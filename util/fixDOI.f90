! Copyright 2018, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=================================
program fixDOI ! tests MLSHDF5 routines
!=================================

  use Dump_0, only: Dump
  use HDF, only: Dfacc_Create, Dfacc_Rdwr, Dfacc_Read
  use HDF5, only: H5fis_HDF5_F, H5gclose_F, H5gopen_F
  use HDFEOS5, only: He5_Swclose, He5_Swopen, MLS_CharType
  use Machine, only: Hp, Getarg, NeverCrash
  use MLSCommon, only: R4, MLSFile_T
  use MLSFILES, only: InitializeMLSFile, MLS_SFStart, MLS_SFEnd, &
    & HDFVersion_5
  use MLSHDF5, only: SaveAsHDF5DS, LoadFromHDF5DS, MakeHDF5Attribute, &
    & GetHDF5Attribute, IsHDF5AttributePresent, MLS_H5open, MLS_H5close
  use MLSHDFEOS, only: He5_Ehrdglatt, MLS_Ehwrglatt, MLS_Swcreate
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMessageConfig
  use MLSFinds, only: FindAll
  use Io_Stuff, only: Get_Lun, Read_TextFile
  use Output_M, only: NewLine, Output
  use HDF5, only: H5gclose_F, H5gopen_F, H5dopen_F, H5dclose_F
  use HDFEOS5, only: He5_Swclose, He5_Swopen, HE5F_ACC_RDWR
   
   implicit none

!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

! Brief description of program
! This program fixes level 2 v4.2x products with erroneous DOI identifiers
! It does so for both l2aux and l2gp files

! What needs correction
! l2aux                                   l2gp
! -----                                   ----
! ATTRIBUTE "identifier_product_DOI" {    ATTRIBUTE "identifier_product_DOI" {
! DATASET "coremetadata.0" {              DATASET "coremetadata.0" {
! DATASET "xmlmetadata" {                 DATASET "xmlmetadata" {

! So they're the same, although their apis are different
! Within the metadata, the fields that need attention are (in odl)
! 
!    OBJECT                 = LOCALVERSIONID
!      NUM_VAL              = 1
!      VALUE                = "c01"
!    END_OBJECT             = LOCALVERSIONID
!
!    OBJECT                 = IDENTIFIER_PRODUCT_DOI
!      NUM_VAL              = 1
!      VALUE                = "10.5067/AURA/MLS/DATA2008"
!    END_OBJECT             = IDENTIFIER_PRODUCT_DOI
!
!   OBJECT                 = LOCALGRANULEID
!     NUM_VAL              = 1
!     VALUE                = "MLS-Aura_L2GP-CO_v04-23-c01_2018d064.he5"
!   END_OBJECT             = LOCALGRANULEID
!
! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"

! Then run it
! cp /data/emls/l2gp/v04.23/2018/064/MLS-Aura_L2GP-CO_v04-23-c01_2018d064.he5 . ; echo "MLS-Aura_L2GP-CO_v04-23-c01_2018d064.he5" | NAG.Linux-6.2/test

! rm coremetadata.0; l2auxdump -d "HDFEOS INFORMATION/coremetadata.0" MLS-Aura_L2GP-BrO_v04-23-c02_2018d064.he5 | sed 's/@/\n/g' > coremetadata.0
! rm xmlmetadata ; l2auxdump -d "HDFEOS INFORMATION/xmlmetadata" MLS-Aura_L2GP-BrO_v04-23-c02_2018d064.he5 | sed 's/@/\n/g' > xmlmetadata
  type options_T
    logical             :: dumpcore           = .false. ! Print (lots) extra
    logical             :: dumpxml            = .false. ! Print (lots) extra
    logical             :: verbose            = .false. ! Print (lots) extra
    character(len=1024) :: corefile           = ''
    character(len=1024) :: xmlfile            = ''
    character(len=255)  :: DOIvalue           = '' 
  end type options_T
  
  type ( options_T ) :: options

! Variables
  integer, parameter ::          MAXFILES = 40
  character(len=255), dimension(MAXFILES) :: filenames
   integer                      :: hdfVersion
   integer                      :: i, nLines
  logical     ::                 is_hdf5
  integer            ::          n_filenames
   integer                      :: nTimes, time, level
   integer, parameter           :: nLevels = 10
   integer                      :: returnStatus, fileID, grp_id
   integer                      :: end_time, start_time
   character(len=1024)          :: Filename
   character(len=*), parameter  :: GlAttrN = 'global attribute'
   character(len=*), parameter  :: GlAttrV = 'GA value'
   character(len=*), parameter  :: DSname = 'xmlmetadata'
   character(len=*), parameter  :: DSname2 = 'coremetadata.0'
   character(len=*), parameter  :: DSUnits = 'test units'
   character(len=*), parameter  :: DimName = 'test dimension'
   character(len=*), parameter  :: DimUnits = 'dimension units'
   character(len=17), dimension(3)    :: myDimensions
   character(len=*), dimension(3), parameter                :: Dimensions = (/&
     &    '001 dimension 001', &
     &    '002 dimension 002', &
     &    '003 dimension 003' /)
   character(len=65535)         :: coremetadata
   character(len=64)            :: DOIValue
   real(r4), dimension(nLevels) :: DimValues
   character(len=*), parameter  :: Charname = 'character data'
   character(len=*), parameter  :: intName = 'integer data'
   real(r4), parameter          :: FILLVALUE=-999.99
   character(len=*), parameter  :: FILLCHAR='*'
   character(len=1000000)       :: textString
   character(len=255), dimension(10000)       :: textArray
   character(len=8), dimension(2) :: int_chars
   integer, dimension(2)        :: the_ints
   integer                      :: itemId
   real(r4)                     :: diff
   logical                      :: is_it
   character(len=1024)          :: textFile
   logical, parameter           :: WRITEATTRIBUTESDURING = .false.
   logical, parameter           :: SKIPWRITINGMIDDLE = .TRUE.
   type(MLSFile_T) ::              l2gp
   ! Executable
   neverCrash = .false.
   call mls_h5open(returnStatus)
   MLSMessageConfig%useToolkit = .false.
   MLSMessageConfig%LogFileUnit = -1
   MLSMessageConfig%crashOnAnyError = .true.
   hdfVersion = HDFVERSION_5
   ! print *, 'Name of file'
   !read(*, '(a80)') Filename
  n_filenames = 0
  do      ! Loop over filenames
     call get_filename( filename, n_filenames, options )
     if ( filename(1:1) == '-' ) cycle
     if ( filename == ' ' ) exit
     call h5fis_hdf5_f( trim(filename), is_hdf5, returnStatus )
     if ( .not. is_hdf5 ) then
       print *, 'Sorry--not recognized as hdf5 file: ', trim( filename )
       cycle
     endif
     n_filenames = n_filenames + 1
     filenames(n_filenames) = filename
  enddo
  if ( n_filenames == 0 ) then
    if ( options%verbose ) print *, 'Sorry no input files supplied'
  elseif ( .not. (options%dumpcore .or. options%dumpxml .or. &
    & len_trim(options%corefile) > 1 .or. len_trim(options%xmlfile) > 1 .or. &
    & len_trim(options%DOIvalue) > 1 ) ) then
    options%dumpcore = .true.
    options%dumpxml  = .true.
  endif
  if ( options%verbose ) call dumpSettings ( options, n_filenames, filenames ) 
   Filename = adjustl(filenames(1))
   textFile = trim(FileName) // '.xml'
   if ( options%verbose ) then
     print *,'hdf version   : ', hdfVersion
     print *,'hdf file name : ', trim(Filename)
     print *,'text file     : ', trim(textFile)
     print *,'len(text file)  : ', len_trim(textFile)
   endif
   ! Is the file plain hdf or hdfeos?
   if ( index( trim(Filename), '.he5' ) < 1 ) then
     fileID = mls_sfstart ( trim(Filename), DFACC_RDWR, &
         &                                          hdfVersion=hdfVersion )
     if ( options%dumpcore ) then
       if ( options%verbose ) print *,'DS name  : ', len_trim(DSname2)
       coremetadata = ' '
       call LoadFromHDF5DS ( fileID, DSname2, coremetadata )
       call output ( trim(coremetadata), advance='yes' )
     endif
     if ( options%dumpxml ) then
       if ( options%verbose ) print *,'DS name  : ', len_trim(DSname)
       coremetadata = ' '
       call LoadFromHDF5DS ( fileID, DSname, coremetadata )
       call output ( trim(coremetadata), advance='yes' )
     endif
     if ( len_trim(options%corefile) > 0 ) then
       call SaveAsHDF5DS ( trim(options%corefile), fileID, DSname2, 4096 )
     endif
     if ( len_trim(options%xmlfile) > 0 ) then
       call SaveAsHDF5DS ( trim(options%xmlfile), fileID, DSname, 4096 )
     endif
     if ( len_trim(options%DOIvalue) > 0 ) then
       call MakeHDF5Attribute ( FileID, "identifier_product_DOI", &
         & trim(options%DOIvalue) )
     endif
     stop

     returnStatus = mls_sfend( fileID, hdfVersion=hdfVersion)
     if ( returnStatus /= 0 ) then
       print *, 'Error in ending hdf access to file'
       stop
     endif
   else
     ! Now let's try to do this with hdfeos files
     returnStatus = InitializeMLSFile( l2gp, access=DFACC_RDWR, &
       & name=Filename )

     fileID = he5_swopen( trim(Filename), HE5F_ACC_RDWR )
     DOIValue = ' '
     returnStatus = HE5_EHRDGLATT( fileID, "identifier_product_DOI", DOIValue )
     if ( options%verbose ) then
       print *,'hdfeos file name : ', trim(Filename)
       print *,'text file     : ', trim(textFile)
       print *,'len(text file)  : ', len_trim(textFile)
       call output ( trim(DOIValue), advance='yes' )
     endif
     if ( len_trim(options%DOIvalue) > 0 ) then
       returnStatus = mls_EHwrglatt( fileID, &
       & 'identifier_product_DOI', MLS_CHARTYPE, 1, &
       &  options%DOIvalue )
     endif
     returnstatus = he5_swclose( fileID )
     
     fileID = mls_sfstart ( trim(Filename), DFACC_RDWR, &
         &                                          hdfVersion=hdfVersion )
     ! print *, 'file ID: ', fileID
     if ( options%dumpcore ) then
       if ( options%verbose ) print *,'DS name  : ', trim(DSname2)
       coremetadata = ' '
       call LoadFromHDF5DS ( fileID, &
         & "HDFEOS INFORMATION/" // DSname2, coremetadata )
       call output ( trim(coremetadata), advance='yes' )
     endif
     if ( options%dumpxml ) then
       if ( options%verbose ) print *,'DS name  : ', trim(DSname)
       coremetadata = ' '
       call LoadFromHDF5DS ( fileID, &
         & "HDFEOS INFORMATION/" // DSname, coremetadata )
       call output ( trim(coremetadata), advance='yes' )
     endif
     if ( len_trim(options%corefile) > 0 ) then
       call SaveAsHDF5DS ( trim(options%corefile), fileID, &
         & "HDFEOS INFORMATION/" // DSname2, 4096, adding_to=.true. )
     endif
     if ( len_trim(options%xmlfile) > 0 ) then
       call SaveAsHDF5DS ( trim(options%xmlfile), fileID, &
         & "HDFEOS INFORMATION/" // DSname, 4096, adding_to=.true. )
     endif
   endif
   call mls_h5close(returnStatus)
contains

!------------------------- get_filename ---------------------
    subroutine get_filename( filename, n_filenames, options )
    ! Added for command-line processing
     character(len=*), intent(out)       :: filename          ! filename
     integer, intent(in)                 :: n_filenames
     type ( options_T ), intent(inout)   :: options
     ! Local variables
     integer ::                         error = 1
     integer, save ::                   i = 1
     character(len=16)       :: number
  ! Get inputfile name, process command-line args
  ! (which always start with -)
    do
      call getarg ( i+hp, filename )
      ! print *, i, ' th Arg: ', trim(filename)
      error = 0
      if ( filename(1:1) /= '-' ) exit
      if ( filename(1:3) == '-h ' ) then
        call print_help
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
        exit
      elseif ( filename(1:4) == '-C ' ) then
        options%dumpcore = .true.
        exit
      elseif ( filename(1:4) == '-X ' ) then
        options%dumpxml = .true.
        exit
      elseif ( filename(1:4) == '-D ' ) then
        options%dumpcore = .true.
        options%dumpxml  = .true.
        exit
      else if ( filename(1:3) == '-a ' ) then
        call getarg ( i+1+hp, options%DOIvalue )
        i = i + 1
        exit
      else if ( filename(1:3) == '-c ' ) then
        call getarg ( i+1+hp, options%corefile )
        i = i + 1
        exit
      else if ( filename(1:3) == '-x ' ) then
        call getarg ( i+1+hp, options%xmlfile )
        i = i + 1
        exit
      else if ( filename(1:3) == '-f ' ) then
        call getarg ( i+1+hp, filename )
        i = i + 1
        exit
      else
        call print_help
      end if
      i = i + 1
    end do
    if ( error /= 0 ) then
      call print_help
    endif
    i = i + 1
    
  end subroutine get_filename
!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'Usage:fixDOI [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options:'
      write (*,*) '  -v              => switch on verbose mode'
      write (*,*) '  -C              => dump just coremetadata.0'
      write (*,*) '  -X              => dump just xmlmetadata'
      write (*,*) '  -D              => dump all metadatasets (default)'
      write (*,*) '  -a value        => set DOI attribute to value'
      write (*,*) '  -c file         => set coremetadata.0 to contents of file'
      write (*,*) '  -x file         => set xmlmetadata to contents of corefile'
      write (*,*) '  -h              => print brief help'
      stop
  end subroutine print_help
!------------------------- dumpSettings ---------------------
    subroutine dumpSettings( options, n_filenames, filenames )
    ! Added for command-line processing
     integer, intent(in)              :: n_filenames
     character(len=255), dimension(:) :: filenames
     type ( options_T ), intent(in)   :: options
     ! Local variables
     integer :: i
     print *, 'dumpcore?            ', options%dumpcore
     print *, 'dumpxml ?            ', options%dumpxml    
     print *, 'verbose ?            ', options%verbose 
     print *, 'corefile             ', trim(options%corefile   )
     print *, 'xmlfile              ', trim(options%xmlfile    )
     print *, 'DOIvalue             ', trim(options%DOIvalue   )
     do i=1, n_filenames
       print *, i, trim(filenames(i))
     enddo
    end subroutine dumpSettings

!==================
end program fixDOI
!==================

! $Log$
! Revision 1.1  2018/05/04 16:36:39  pwagner
! First commit
!
