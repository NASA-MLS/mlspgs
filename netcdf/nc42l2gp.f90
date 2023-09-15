! Copyright 2023, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=================================
program nc42l2gp ! Rewrite Netcdf4 L2GPData files as hdfeos
!=================================

   use Dump_1, only: Dump
   use Dump_Options, only: SDFormatDefault, DumpDumpOptions
   use HDF, only: Dfacc_Read, Dfacc_Create    
   use HDF5, only: H5F_Acc_RdOnly_F, &
     & H5fclose_F, H5fopen_F, H5gopen_F, H5gclose_F, H5fis_HDF5_F
   use HDFEOS5, only: MLS_Chartype
   use HighOutput, only: OutputNamedValue
   use Intrinsic, only: L_HDF, L_Swath, L_NetCDF4
   use L2GPData, only: L2GPData_T, L2GPnamelen, Maxswathnamesbufsize, Rgp, &
     & Dump, DestroyL2GPcontents, WriteL2GPData
   use Machine, only: Hp, Getarg
   use MLSCommon, only: MLSFile_T, Split_path_name
   use MLSFiles, only: HDFversion_5, InitializeMLSFile, &
     & MLS_CloseFile, MLS_OpenFile, Split_Path_Name
   use MLSHDF5, only: IsHDF5DSPresent, LoadFromHDF5DS, MLS_H5open, MLS_H5close
   use MLSHDFeos, only: MLS_Isglatt, He5_Ehrdglatt, He5_EHWrGlAtt, MLS_EhwrGlAtt
   use MLSMessageModule, only: MLSMSG_Error, MLSMSG_Warning, &
     & MLSMessage
   use MLSNetCDF4, only: MLS_Swrdattr, MLS_SwWrattr, MLS_InqSwath
   use MLSStringLists, only: GetStringElement, NumStringElements, &
     & ReplaceSubstring
   use MLSStrings, only: Lowercase, Reverse
   use NCL2GP, only: ReadNCGlobalAttr, ReadNCL2GPData
   use NetCDF, only: NF90_Char, NF90_Open, NF90_Def_Dim, NF90_Def_Grp, &
     & NF90_Def_Var, NF90_Put_Var, NF90_StrError, NF90_Write
   use Optional_M, only: Default
   use Output_M, only: Blanks, Newline, Output, &
     & ResumeOutput, SuspendOutput, SwitchOutput
   use PCFHdr, only: HE5_WriteGlobalAttr
   use Printit_M, only: Set_Config
   implicit none

!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

! Brief description of program
! This program rewrites L2GPData files as netcdf4
!
! A test of the NCL2GPData module
! To build it, 

! From the root, mlspgs, level, just type
!   make nc42l2gp

! Testing has taken place in
! /users/pwagner/mlspgs/tests/lib/nctest
! 

  integer, parameter :: MetaDataSize = 65535 

  type Options_T
     logical ::             debug         = .false.
     logical ::             verbose       = .false.
     logical ::             attributesToo = .false.
     logical ::             metaDataToo   = .false.
  end type Options_T

  type ( Options_T ) :: options
  character(len=255) :: filename          ! filename
  integer            :: n_filenames
  integer     ::  error ! Counting indices & Error flags
  logical     :: is_hdf5
  logical     :: is_present
  integer, save                   :: numGood = 0
  integer, save                   :: numGoodPrec = 0
  integer, save                   :: numNotUseable = 0
  integer, save                   :: numOddStatus = 0
  integer, save                   :: numPostProcStatus = 0
  real, dimension(3), save        :: numTest = 0.
  integer, parameter              :: PostProcBitIndex = 4
  integer, parameter              :: MAXNUMBITSUSED = 10 !9
  ! The bit number starts at 0: bitNumber[1] = 0
  integer, dimension(MAXNUMBITSUSED), parameter :: bitNumber = &
    & (/ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 /)
  integer, dimension(MAXNUMBITSUSED, 2), save :: bitCounts = 0
  character(len=*), parameter     :: bitNames = &
    & '  dontuse,   bewary,     info,postprocd,' // &
    & '    hicld,    locld,   nogmao,abandoned,   toofew,    crash'
  !   01234567890123456789012345678901234567890123456789012345678901234567890123456789
  real(rgp), dimension(:), pointer :: values => null()
  ! Executable
  call set_config ( useToolkit = .false., logFileUnit = -1 )
  call switchOutput( 'stdout' )
  call mls_h5open(error)
  n_filenames = 0
  do      ! Loop over filenames
     call get_filename(filename, n_filenames, options)
     if ( filename == ' ' ) exit
     n_filenames = n_filenames + 1
     call h5fis_hdf5_f(trim(filename), is_hdf5, error)
     if ( .not. is_hdf5 ) then
       print *, 'Sorry--not recognized as hdf5 file: ', trim(filename)
     endif
     if ( options%verbose ) then
       print *, 'Rewriting swaths in     ', trim(filename)
       print *, 'attributes copied, too? ', options%attributesToo
       print *, 'metadata copied, too?   ', options%metaDataToo
     endif
     call OutputNamedValue( 'HDFEOS5 L2GP File',  trim(filename), &
       & options='--Headline' )
     call rewrite_one_file( trim(filename), options )
     call resumeOutput
  enddo
  call mls_h5close(error)
contains
!------------------------- get_filename ---------------------
    subroutine get_filename( filename, n_filenames, options )
    ! Added for command-line processing
     character(len=255), intent(out) :: filename          ! filename
     integer, intent(in) ::             n_filenames
     type ( Options_T ) :: options
     integer ::                         error = 1
     integer, save ::                   i = 1
     character(len=255) :: argstr
  ! Get inputfile name, process command-line args
  ! (which always start with -)
    do
      call getarg ( i+hp, filename )
      ! print *, i, ' th Arg: ', trim(filename)
      error = 0
      if ( filename(1:1) /= '-' ) exit
      if ( filename(1:3) == '-h ' ) then
        call print_help
      elseif ( filename(1:3) == '-d ' ) then
        options%debug   = .true.
        options%verbose = .true.
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
      elseif ( filename(1:3) == '-a ' ) then
        options%attributesToo = .true.
      elseif ( filename(1:3) == '-m ' ) then
        options%metadataToo = .true.
      else if ( filename(1:3) == '-f ' ) then
        call getarg ( i+1+hp, filename )
        error = 0
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
    if (trim(filename) == ' ' .and. n_filenames == 0) then

    ! Last chance to enter filename
      print *,  "Enter the name of the HDF5 file. " // &
       &  "Datasets in the file will be listed shortly."
      read(*,'(a)') filename
    endif
    
  end subroutine get_filename
!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'l2gpnc4: Converts each netcdf4 file to hdfeos format' &
      & // 'replacing the .nc4 suffix with .he5' &
      & // 'Usage: nc42l2gp [options] [netcdf_filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options:'
      write (*,*) ' -f filename => use filename'
      write (*,*) ' -h          => print brief help'
      write (*,*) ' -a          => copy attributes, too'
      write (*,*) ' -m          => copy metadata, too'
      write (*,*) ' -v          => verbose'
      write (*,*) ' -d          => debug'

      stop
  end subroutine print_help
  
  subroutine Defeat(msg)
    character(len=*), intent(in) :: msg
    call resumeOutput
    call output('Serious error: ' // msg, advance='yes')
    call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'stopping' )
  end subroutine Defeat
  
  subroutine rewrite_one_file ( filename, options )
    use MLSStrings, only:  Asciify
    use PCFHdr, only:  GlobalAttributes_T, HE5_ReadGlobalAttr, &
      & DumpGlobalAttributes, GlobalAttributes
    type (GlobalAttributes_T)            :: gAttributes
    ! Dummy args
    character(len=*), intent(in)         :: filename ! filename
    type ( Options_T )                   :: options
    ! Local variables
    logical, parameter                   :: countEmpty = .true.
    integer :: File1
    integer                              :: i
    integer                              :: listsize
    type (L2GPData_T)                    :: l2gp
    character (len=MAXSWATHNAMESBUFSIZE) :: matches
    type(MLSFile_T)                      :: MLSFile
    type(MLSFile_T)                      :: NC4File
    character(len=255)                   :: ncfilename ! filename
    integer                              :: noSwaths
    integer                              :: status
    character (len=L2GPNameLen)          :: swath
    character (len=MAXSWATHNAMESBUFSIZE) :: SwathList
    character(len=40)                    :: ProcessLevel       
    double precision                     :: TAI93At0zOfGranule 
    integer                              :: DayofYear          
    character (len=L2GPNameLen)          :: HostName
    character (len=MAXSWATHNAMESBUFSIZE) :: MiscNotes
    character (len=L2GPNameLen)          :: DOI
    ! Get swath list
    noSwaths = mls_InqSwath ( filename, SwathList, listSize )
    if ( options%verbose ) then
      call output('swaths in ' // trim(filename), advance='yes')
      call dump(trim(swathList))
      ! return
    endif
    ! Executable code
    noSwaths = NumStringElements(trim(swathList), countEmpty)
    if ( noSwaths < 1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'No swaths to dump--unable to count swaths in ' // trim(swathList) )
    endif
    status = InitializeMLSFile ( NC4File, type=L_NetCDF4, access=DFACC_READ, &
     & name=filename )
    ! Create hdfeos file:
    ! Assume the netcdf file name is some_name.nc4
    ! We want                        some_name.he5
    ! A trick!
    ! Reversing the filenames means simply transform
    ! 4cn.whatever -> 5eh.whatever
    ! which is easy with substrings
    ncfilename = Reverse(trim(filename))
    ncfilename = '5eh' // ncfilename(4:)
    ncfilename = Reverse(trim(ncfilename))
    call outputNamedValue( 'HDFEOS file name', ncfilename )
    status = InitializeMLSFile ( NC4File, type=l_swath, access=DFACC_CREATE, &
     & name=ncfilename, hdfVersion=HDFVERSION_5 )
    if ( options%verbose ) call Dump ( NC4File )
    ! Loop over swaths in file 1
    do i = 1, noSwaths
      call GetStringElement (trim(swathList), swath, i, countEmpty )
      ! HDFEOS INFORMATION isn't a swath
      if ( trim(swath) == 'HDFEOS INFORMATION' ) cycle
    
      ! Allocate and fill l2gp
      ! print *, 'Reading swath from file: ', trim(swath)
      call ReadNCL2GPData ( trim(filename), trim(swath), l2gp )
      ! Rewrite the swath- and file-level attributes
      call output ( 'About to write the l2gp ' // &
        & trim(swath) // '  to the hdfeos file', advance='yes' )
      call WriteL2GPData( l2gp, MLSFile, trim(swath) )
      ! Should we default to writing attributes and metadata?
      ! Maybe not metadata unless we figure out how to change
      ! the suffix ".he5" in LocalGranuleID to ".nc4"
      ! We're trying. We'll see how successful we are.
      if ( options%attributesToo .and. i == 1 ) then
        call MLS_OpenFile( MLSFile )
        call MLS_OpenFile( NC4File )
        file1 = MLSFile%FileID%f_id
        call output ( '(Global Attributes) ', advance='yes' )
        call ReadNCGlobalAttr ( NC4File%FileID%f_id, gAttributes, status )
!        call he5_readglobalattr( file1, gAttributes, &
!         &  ProcessLevel, DayofYear, TAI93At0zOfGranule, &
!         & HostName, MiscNotes, DOI, returnStatus=status )
        GlobalAttributes = gAttributes
        if ( options%verbose ) then
          call DumpGlobalAttributes
          call OutputNamedValue ( 'DOI', trim(asciify( DOI, 'snip' ) ) )
          call OutputNamedValue ( 'HostName', trim(asciify( HostName, 'snip' ) ) )
          call OutputNamedValue ( 'MiscNotes', trim(asciify( MiscNotes, 'snip' ) ) )
          call OutputNamedValue ( 'TAI93', TAI93At0zOfGranule )
        endif
        GlobalAttributes%HostName           = asciify( HostName, 'snip' )
        GlobalAttributes%ProductionLoc      = asciify( HostName, 'snip' )
        GlobalAttributes%MiscNotes          = asciify( MiscNotes, 'snip' )
        GlobalAttributes%DOI                = asciify( DOI, 'snip' )     
        GlobalAttributes%TAI93At0zOfGranule = TAI93At0zOfGranule
        if ( status /= 0 ) then
          call output ('No global attributes found in file', advance='yes')
        endif
        if ( options%verbose ) then
          call dump( file1, l2gp )
          print *, 'Attempting to copy attributes, too'
          print *, 'ProductionLocation: ', trim(GlobalAttributes%ProductionLoc)
        endif
!        call WriteNCGlobalAttr ( NC4File%FileID%f_id, &
!          & DOI=.true., MiscNotes=.true. )
        call he5_writeglobalattr( MLSFile%FileID%f_id, &
         &  DOI=.true., skip_if_already_there=.false. )
        call copyAprioriFileNamesAttr ( MLSFile%FileID%f_id, NC4File%FileID%f_id )
        call MLS_CloseFile ( MLSFile )
        call MLS_CloseFile ( NC4File )
      endif
      if ( options%metaDataToo .and. i == 1 ) then
        ! For this operation we'll treat the hdfeos file like an hdf
        MLSFile%type = l_hdf 
        call output ( '(Meta data) ', advance='yes' )
        call copyMetaData ( MLSFile, nc4File )
      endif
      call DestroyL2GPContents ( l2gp )
    enddo
   end subroutine rewrite_one_file

   subroutine  copyMetaData ( hdf, nc4 )
     ! Copy the metadata datasets
     ! coremetadata.0
     ! xmlmetadata
     ! They will be put into a new group 
     ! "HDFEOS INFORMATION"
     ! Args
     type(MLSFile_T), intent(in)          :: hdf
     type(MLSFile_T), intent(in)          :: nc4
     ! Internal variables
     integer                              :: coredimid
     integer                              :: corevarid
     integer                              :: hdffid
     integer                              :: ncfid
     integer                              :: hdfgrpid
     integer                              :: ncgrpid
     integer                              :: status
     character(len=1024)                  :: path
     character(len=MetaDataSize)          :: strvalue
     character(len=MetaDataSize)          :: tmpstr
     integer, dimension(1)                :: dimids ! For scalar variables
     integer                              :: xmldimid
     integer                              :: xmlvarid
     character(len=128)                   :: hdfname
     character(len=128)                   :: nc4name
     ! Executable
     if ( options%debug ) print *, 'Attempting to copyMetaData'
     call h5fopen_f ( trim(hdf%name), H5F_ACC_RDONLY_F, hdffID, status )
     if ( .not. &
       & IsHDF5DSPresent ( HDF, '/HDFEOS INFORMATION/coremetadata.0' ) &
       & ) then
       call h5fclose_f ( hdffID, status )
       print *, 'cormetadata.0 not found in ', trim(HDF%name)
       return 
     endif
     call split_path_name( trim(hdf%name), path, name=hdfname)
     call split_path_name( trim(nc4%name), path, name=nc4name)
     
     call h5gopen_f ( hdffID, '/HDFEOS INFORMATION', hdfgrpID, status )
     if ( status /= 0 ) call terror ( 'Opening hdf5 group', status )

     status = nf90_open ( nc4%NAME, NF90_Write, ncfID )
     if ( status /= 0 ) call terror ( 'Opening netcdf file', status )
     ! Create nc4 top group
     status = nf90_def_grp( ncfId, 'HDFEOS INFORMATION', ncgrpID )
     if ( status /= 0 ) call terror ( 'Creating hdf5 group HDFEOS INFORMATION', status )
     
     ! Write the 1st metadata
     call LoadFromHDF5DS ( hdfgrpID, 'coremetadata.0', strvalue )
     ! Can we replace hdfname with nc4name? Let's try.
     tmpstr = strvalue
     call ReplaceSubString ( tmpstr, strvalue, trim(hdfname), trim(nc4name) )
     if ( options%debug ) call output( strvalue(1:80), advance='yes' )
     ! Why do we need dims? The metadata are scalars. Blame netcdf
     ! for using c?
     status = nf90_def_dim( ncgrpid, "coresize", len_trim(strvalue), coredimid )
     dimids = coredimid
     status = nf90_def_var( ncgrpID, "coremetadata.0", nf90_char, dimids, corevarid )
     if ( status /= 0 ) call terror ( 'Creating var coremetadata.0', status )

     status = nf90_put_var( ncgrpID, corevarid, trim(strvalue) )
     if ( status /= 0 ) call terror ( 'Writing var coremetadata.0', status )

     ! Write the 2nd metadata
     call LoadFromHDF5DS ( hdfgrpID, 'xmlmetadata', strvalue )
     ! Can we replace hdfname with nc4name? Let's try.
     tmpstr = strvalue
     call ReplaceSubString ( tmpstr, strvalue, trim(hdfname), trim(nc4name) )
     if ( options%debug ) call output( strvalue(1:80), advance='yes' )
     status = nf90_def_dim( ncgrpid, "xmlsize", len_trim(strvalue), xmldimid )
     dimids = xmldimid
     status = nf90_def_var( ncgrpID, "xmlmetadata", nf90_char, dimids, xmlvarid )
     if ( status /= 0 ) call terror ( 'Creating var xmlmetadata', status )

     status = nf90_put_var( ncgrpID, xmlvarid, trim(strvalue) )
     if ( status /= 0 ) call terror ( 'Writing var xmlmetadata.0', status )

   end subroutine  copyMetaData
   subroutine terror ( message, status )
     character(len=*), intent(in) :: message
     integer, intent(in)          :: status
     print *, 'status: ', status
     call output( trim(nf90_strerror(status)), advance='yes' )
     print *, message
     stop
   end subroutine terror

   subroutine  copyAprioriFileNamesAttr ( hdfeosid, nc4id )
     ! Copy the global attributes saying what
     ! kind of meteorology apriori files were used, if any
     ! Args
     integer, intent(in)                  :: hdfeosid
     integer, intent(in)                  :: nc4id
     ! Internal variables
     integer                              :: status
     character(len=Maxswathnamesbufsize)  :: strvalue
     ! Executable
     if ( options%debug ) print *, 'Attempting to copyAprioriNames'
     if ( .not. MLS_Isglatt ( hdfeosID, 'A Priori l2gp' ) ) then
       print *, 'A Priori l2gp attribute not found in ', hdfeosID
       return 
     endif
     status = MLS_Swrdattr ( nc4id, &
     & 'A Priori l2gp', nf90_char, 1, strvalue )
     status = MLS_EhwrGlAtt( hdfeosID, &
      & 'A Priori l2gp', MLS_CHARTYPE, 1, strvalue )

     status = MLS_Swrdattr ( nc4id, &
     & 'A Priori l2aux', nf90_char, 1, strvalue )
     status = MLS_EhwrGlAtt( hdfeosID, &
      & 'A Priori l2aux', MLS_CHARTYPE, 1, strvalue )

     status = MLS_Swrdattr ( nc4id, &
     & 'A Priori ncep', nf90_char, 1, strvalue )
     status = MLS_EhwrGlAtt( hdfeosID, &
      & 'A Priori ncep', MLS_CHARTYPE, 1, strvalue )

     status = MLS_Swrdattr ( nc4id, &
     & 'A Priori gmao', nf90_char, 1, strvalue )
     status = MLS_EhwrGlAtt( hdfeosID, &
      & 'A Priori gmao', MLS_CHARTYPE, 1, strvalue )

     status = MLS_Swrdattr ( nc4id, &
     & 'A Priori geos5', nf90_char, 1, strvalue )
     status = MLS_EhwrGlAtt( hdfeosID, &
      & 'A Priori geos5', MLS_CHARTYPE, 1, strvalue )

     status = MLS_Swrdattr ( nc4id, &
     & 'geos5 type', nf90_char, 1, strvalue )
     status = MLS_EhwrGlAtt( hdfeosID, &
      & 'geos5 type', MLS_CHARTYPE, 1, strvalue )

     status = MLS_Swrdattr ( nc4id, &
      & 'PCF1', nf90_char, 1, strvalue )
     status = MLS_EhwrGlAtt( hdfeosID, &
      & 'PCF1', MLS_CHARTYPE, 1, strvalue )

     ! Can it be that the hdfeos5 library treats strvalue as inout?
     ! strvalue = ' '
     status = MLS_Swrdattr ( nc4id, &
      & 'PCF2', nf90_char, 1, strvalue )
     status = MLS_EhwrGlAtt( hdfeosID, &
      & 'PCF2', MLS_CHARTYPE, 1, strvalue )

   end subroutine  copyAprioriFileNamesAttr

!==================
end program nc42l2gp
!==================

! $Log$
