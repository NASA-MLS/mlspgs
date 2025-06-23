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
program fixAttribute ! tests MLSHDF5 routines
!=================================

  use HDF, only: Dfacc_Rdwr
  use HDF5, only: H5F_ACC_RDWR_F, &
    & H5fis_HDF5_F, H5FOpen_f, H5FClose_f, H5GOpen_f, H5GClose_f
  use HDFEOS5, only: He5_Swclose, He5_Swopen, MLS_CharType
  use Machine, only: Hp, Getarg, NeverCrash
  use MLSCommon, only: MLSFile_T
  use MLSFiles, only: InitializeMLSFile, MLS_SFStart, MLS_SFEnd, &
    & HDFVersion_5
  use MLSHDF5, only: SaveAsHDF5DS, LoadFromHDF5DS, MakeHDF5Attribute, &
    & GetHDF5Attribute, MLS_H5open, MLS_H5close
  use MLSHDFEOS, only: He5_Ehrdglatt, MLS_Ehwrglatt
  use MLSMessageModule, only: MLSMessageConfig
  use Output_M, only: Output
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
! This program does any one of the following:
! It does so for both l2aux and l2gp files
! (1) fixes level 2 umls 1.50 products with erroneous AttributeName identifiers
! (2) overwrites erroneous level 2 metadata
! (3) dumps level 2 metadata

! What needs correction
! l2aux                                   l2gp
! -----                                   ----
! ATTRIBUTE trim(options%AttributeName) {    ATTRIBUTE trim(options%AttributeName) {

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"

! Then run it
! cd ~/mlspgs/tmp
! cp *001.he5* 001 ; /users/pwagner/mlspgs/tests/lib/NAG.Linux-6.2-cool/test -v -A -a "MLS UARS" 001/*.he5
! cp *001.he5* 001 ; /users/pwagner/mlspgs/tests/lib/NAG.Linux-6.2-cool/test -v -C 001/*.he5 > odl.txt
! cp *001*.h5 001 ; /users/pwagner/mlspgs/tests/lib/NAG.Linux-6.2-cool/test -v -A -a "MLS UARS" 001/*.h5
! cp *001.h5* 001 ; /users/pwagner/mlspgs/tests/lib/NAG.Linux-6.2-cool/test -X 001/*.h5 > xml.txt
  type options_T
    logical             :: debug              = .false. ! Print debugging stuff
    logical             :: dumpcore           = .false. ! Print odl-formatted metadata
    logical             :: dumpxml            = .false. ! Print xml-formatted metadata
    logical             :: verbose            = .false. ! Print (lots) extra
    character(len=1024) :: corefile           = ''
    character(len=1024) :: xmlfile            = ''
    character(len=255)  :: AttributeValue     = '' ! 'MLS UARS' 
    character(len=32)   :: AttributeName      = 'InstrumentName' ! or 'PGEVersion'
    character(len=1)    :: AorM               = '' ! 'A' 
  end type options_T
  
  type ( options_T ) :: options

! Variables
  integer, parameter ::          MAXFILES = 40
  character(len=255), dimension(MAXFILES) :: filenames
   integer                      :: hdfVersion
   integer                      :: i
  logical     ::                 is_hdf5
  integer            ::          n_filenames
   integer                      :: returnStatus, fileID
   character(len=1024)          :: Filename
   character(len=*), parameter  :: DSname = '/HDFEOS INFORMATION/xmlmetadata'
   character(len=*), parameter  :: DSname2 = '/HDFEOS INFORMATION/coremetadata.0'
   character(len=65535)         :: coremetadata
   character(len=64)            :: AttributeValue
   character(len=1024)          :: textFile
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
    & len_trim(options%AttributeValue) > 1 ) ) then
    options%dumpcore = .true.
    options%dumpxml  = .true.
  endif
  if ( options%verbose ) call dumpSettings ( options, n_filenames, filenames ) 
  Filename = adjustl(filenames(1))
  textFile = trim(FileName) // '.xml'
  if ( options%debug ) then
    print *,'hdf version   : ', hdfVersion
    print *,'hdf file name : ', trim(Filename)
    print *,'text file     : ', trim(textFile)
    print *,'len(text file)  : ', len_trim(textFile)
  endif
  ! print *, 'Name of file ' // trim(Filename)
  ! Now are we fixing the file level attribute or the metadata?
  ! Or are we dumping metadata?
  if ( options%dumpcore .or. options%dumpxml ) options%AorM = 'D'
  if ( ( len_trim(options%corefile) > 0 .or. &
    &  len_trim(options%xmlfile) > 0 )  .and. &
    &  len_trim(options%AorM) < 1 &
    & ) options%AorM = 'I' ! "I" for "Insert metadat"
  select case ( options%AorM )
  case ( 'A' )
    call TheAttribute
  case ( 'M', 'I' )
    call TheMetadata
  case ( 'D' )
    call Dump
  case default
    print *, 'Sorry; you must run with either -A or -M option set'
    stop
  end select

  call mls_h5close(returnStatus)
contains

!------------------------- TheAttribute ---------------------
  subroutine TheAttribute
   ! Is the file plain hdf or hdfeos?
   if ( index( trim(Filename), '.he5' ) < 1 ) then
     fileID = mls_sfstart ( trim(Filename), DFACC_RDWR, &
         &                                          hdfVersion=hdfVersion )
     if ( len_trim(options%corefile) > 0 ) then
       call SaveAsHDF5DS ( trim(options%corefile), fileID, DSname2, 4096 )
     endif
     if ( len_trim(options%xmlfile) > 0 ) then
       call SaveAsHDF5DS ( trim(options%xmlfile), fileID, DSname, 4096 )
     endif
     call GETHDF5Attribute ( FileID, trim(options%AttributeName), &
       & AttributeValue )
     if ( options%verbose ) then
       print *,'l2aux file name      : ', trim(Filename)
       print *,'OldAttributeValue    : ', trim(AttributeValue)
     endif
     if ( len_trim(options%AttributeValue) > 0 ) then
       call MakeHDF5Attribute ( FileID, trim(options%AttributeName), &
         & trim(options%AttributeValue) )
       call GETHDF5Attribute ( FileID, trim(options%AttributeName), &
         & AttributeValue )
       if ( options%verbose ) &
         & print *,'NewAttributeValue    : ', trim(AttributeValue)
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
     AttributeValue = ' '
     returnStatus = HE5_EHRDGLATT( fileID, trim(options%AttributeName), AttributeValue )
     if ( options%verbose ) then
       print *,'hdfeos file name     : ', trim(Filename)
       print *,'OldAttributeValue    : ', trim(AttributeValue)
     endif
     if ( len_trim(options%AttributeValue) > 0 ) then
       returnStatus = mls_EHwrglatt( fileID, &
       & trim(options%AttributeName), MLS_CHARTYPE, 1, &
       &  options%AttributeValue )
       returnStatus = HE5_EHRDGLATT( fileID, trim(options%AttributeName), AttributeValue )
       if ( options%verbose ) &
       print *,'NewAttributeValue    : ', trim(AttributeValue)
     endif
     returnstatus = he5_swclose( fileID )
     
     fileID = mls_sfstart ( trim(Filename), DFACC_RDWR, &
         &                                          hdfVersion=hdfVersion )
     ! print *, 'file ID: ', fileID
     if ( len_trim(options%corefile) > 0 ) then
       call SaveAsHDF5DS ( trim(options%corefile), fileID, &
         & DSname2, 4096, adding_to=.true., cr2lf=.true. )
     endif
     if ( len_trim(options%xmlfile) > 0 ) then
       call SaveAsHDF5DS ( trim(options%xmlfile), fileID, &
         & DSname, 4096, adding_to=.true., cr2lf=.true. )
     endif
   endif
  end subroutine TheAttribute

!------------------------- TheMetadata ---------------------
  subroutine TheMetadata
    character(len=255) :: DOIValue
    integer            :: grp_id
    integer            :: status
   ! Is the file plain hdf or hdfeos?
   if ( index( trim(Filename), '.he5' ) < 1 ) then
     fileID = mls_sfstart ( trim(Filename), DFACC_RDWR, &
         &                                          hdfVersion=hdfVersion )
     if ( len_trim(options%corefile) > 0 ) then
       call SaveAsHDF5DS ( trim(options%corefile), fileID, DSname2, 4096, &
         & adding_to=.true., cr2lf=.true. )
     endif
     if ( len_trim(options%xmlfile) > 0 ) then
       call SaveAsHDF5DS ( trim(options%xmlfile), fileID, DSname, 4096, &
         & adding_to=.true., cr2lf=.true. )
     endif
     if ( len_trim(options%AttributeValue) > 0 ) then
       call MakeHDF5Attribute ( FileID, "identifier_product_DOI", &
         & trim(options%AttributeValue) )
     endif
     stop

     returnStatus = mls_sfend( fileID, hdfVersion=hdfVersion)
     if ( returnStatus /= 0 ) then
       print *, 'Error in ending hdf access to file'
       stop
     endif
   else
     ! Now let's try to do this with hdfeos files
     if ( .false. ) then
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
     endif
     ! fileID = mls_sfstart ( trim(Filename), DFACC_RDWR, &
     !     &                                          hdfVersion=hdfVersion )
     call h5fopen_f ( trim(FileName), H5F_ACC_RDWR_F, FileID, status )
     if ( status /= 0 ) then
       print *, 'Error in opening hdf access to file'
       stop
     endif
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
     call h5gopen_f( fileID, "HDFEOS INFORMATION", grp_id, status )
     if ( status /= 0 ) then
       print *, 'Error in opening hdf access to group'
       stop
     endif
     if ( len_trim(options%corefile) > 0 ) then
       call SaveAsHDF5DS ( trim(options%corefile), grp_id, &
         & DSname2, 4096 )
     endif
     if ( len_trim(options%xmlfile) > 0 ) then
       call SaveAsHDF5DS ( trim(options%xmlfile), grp_id, &
         & DSname, 4096 )
     endif
     call h5gclose_f( grp_id, status )
!      if ( len_trim(options%AttributeValue) > 0 ) then
!        returnStatus = mls_EHwrglatt( fileID, &
!        & 'identifier_product_DOI', MLS_CHARTYPE, 1, &
!        &  options%AttributeValue )
!      endif
     if ( len_trim(options%AttributeValue) > 0 ) then
       call h5gopen_f( fileID, '/HDFEOS/ADDITIONAL/FILE_ATTRIBUTES', grp_id, status )
       call MakeHDF5Attribute ( grp_ID, "identifier_product_DOI", &
         & trim(options%AttributeValue) )
       call h5gclose_f( grp_id, status )
     endif
     ! returnstatus = he5_swclose( fileID )
     call h5fclose_f(fileID, status)
     
     
   endif
  end subroutine TheMetadata

!------------------------- Dump ---------------------
  subroutine Dump
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
   else
     fileID = mls_sfstart ( trim(Filename), DFACC_RDWR, &
         &                                          hdfVersion=hdfVersion )
     ! print *, 'file ID: ', fileID
     if ( options%dumpcore ) then
       if ( options%verbose ) print *,'DS name  : ', trim(DSname2)
       coremetadata = ' '
       call LoadFromHDF5DS ( fileID, &
         & DSname2, coremetadata )
       call output ( trim(coremetadata), advance='yes' )
     endif
     if ( options%dumpxml ) then
       if ( options%verbose ) print *,'DS name  : ', trim(DSname)
       coremetadata = ' '
       call LoadFromHDF5DS ( fileID, &
         & DSname, coremetadata )
       call output ( trim(coremetadata), advance='yes' )
     endif
     if ( len_trim(options%corefile) > 0 ) then
       call SaveAsHDF5DS ( trim(options%corefile), fileID, &
         & DSname2, 4096, adding_to=.true., cr2lf=.true. )
     endif
     if ( len_trim(options%xmlfile) > 0 ) then
       call SaveAsHDF5DS ( trim(options%xmlfile), fileID, &
         & DSname, 4096, adding_to=.true., cr2lf=.true. )
     endif
   endif
  end subroutine Dump

!------------------------- get_filename ---------------------
    subroutine get_filename( filename, n_filenames, options )
    ! Added for command-line processing
     character(len=*), intent(out)       :: filename          ! filename
     integer, intent(in)                 :: n_filenames
     type ( options_T ), intent(inout)   :: options
     ! Local variables
     integer ::                         error = 1
     integer, save ::                   i = 1
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
      elseif ( filename(1:4) == '-A ' ) then
        options%AorM = 'A'
        exit
      elseif ( filename(1:4) == '-M ' ) then
        options%AorM = 'M'
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
        call getarg ( i+1+hp, options%AttributeValue )
        i = i + 1
        exit
      else if ( filename(1:3) == '-n ' ) then
        call getarg ( i+1+hp, options%AttributeName )
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
      & 'Usage:fixAttribute [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options:'
      write (*,*) '  -A              => operate on attribute AttributeName'
      write (*,*) '  -M              => operate on metadata LocalGranuleID'
      write (*,*) '  -v              => switch on verbose mode'
      write (*,*) '  -C              => dump just coremetadata.0'
      write (*,*) '  -X              => dump just xmlmetadata'
      write (*,*) '  -D              => dump all metadata sets (default)'
      write (*,*) '  -n name         => set attribute "name" to value'
      write (*,*) '  -a value        => set metadata or attribute to value'
      write (*,*) '                     (you must specify either -A or -M)'
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
     print *, 'Attribute or Metadata?  ', options%AorM
     print *, 'dumpcore?               ', options%dumpcore
     print *, 'dumpxml ?               ', options%dumpxml    
     print *, 'verbose ?               ', options%verbose 
     print *, 'corefile                ', trim(options%corefile   )
     print *, 'xmlfile                 ', trim(options%xmlfile    )
     print *, 'Attribute value         ', trim(options%AttributeValue   )
     print *, 'Attribute name          ', trim(options%AttributeName   )
     do i=1, n_filenames
       print *, i, trim(filenames(i))
     enddo
    end subroutine dumpSettings

!==================
end program fixAttribute
!==================

! $Log$
! Revision 1.1  2018/11/14 18:48:03  pwagner
! First commit
!
