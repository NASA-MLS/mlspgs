! Copyright 2020, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=================================
program l2gpdiffnc4 ! Diff L2GPData file from NetCDF4
!=================================

   use Dump_1, only: Dump
   use Dump_Options, only: SDFormatDefault, DumpDumpOptions
   use HDF, only: Dfacc_Read, Dfacc_Create    
   use HDF5, only: H5F_ACC_RDONLY_F, &
     & H5fclose_F, H5fopen_F, H5gopen_F, H5gclose_F, H5fis_HDF5_F
   use HighOutput, only: OutputNamedValue
   use Intrinsic, only: L_HDF, L_Swath, L_NetCDF4
   use L2GPData, only: L2GPData_T, L2GPnamelen, Maxswathnamesbufsize, Rgp, &
     & Dump, ReadL2GPData, DestroyL2GPcontents
   use Machine, only: Hp, Getarg
   use MLSCommon, only: MLSFile_T
   use MLSFiles, only: HDFversion_5, InitializeMLSFile, MLS_Inqswath, &
     & MLS_CloseFile, MLS_OpenFile, Split_Path_Name
   use MLSHDF5, only: IsHDF5DSPresent, LoadFromHDF5DS, MLS_H5open, MLS_H5close
   use MLSHDFeos, only: MLS_Isglatt, He5_Ehrdglatt
   use MLSMessageModule, only: MLSMessageConfig, MLSMSG_Error, MLSMSG_Warning, &
     & MLSMessage
   use MLSNetCDF4, only: MLS_SwWrattr
   use MLSStringLists, only: GetStringElement, NumStringElements
   use MLSStrings, only: Lowercase, Reverse
   use NCL2GP, only: WriteNCGlobalAttr, WriteNCL2GPData
   use NetCDF, only: NF90_Char, NF90_Open, NF90_Def_Dim, NF90_Def_Grp, &
     & NF90_Def_Var, NF90_Put_Var, NF90_StrError, NF90_Write
   use Optional_M, only: Default
   use Output_M, only: Blanks, Newline, Output, &
     & ResumeOutput, SuspendOutput, SwitchOutput
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
! This program diffss L2GPData files from netcdf4
!
! A test of the NCL2GPData module
! To build it, 

! From the root, mlspgs, level, just type
!   make l2gpdoffnc4

! To run the resulting executable
! IFC.Linux.ifc17/test -v -a -m nctest/*.he5 nctest/*.nc4
! 

  ! This is just the maximum num of chunks you wish
  ! to dump individually in case you don't want to dump them all
  ! It's not the actual maximum number of chunks.
  integer, parameter :: MetaDataSize = 65535 

  type Options_T
    character(len=80) ::   dumpOptions     = 's'
    integer     ::         Details         = 2
    logical     ::         force           = .false.
    logical     ::         AuBrick         = .true.
    logical     ::         rms             = .false.
    logical     ::         stats           = .false.
    logical     ::         table           = .false.
    logical     ::         silent          = .false.
    logical     ::         timing          = .false.
    logical     ::         verbose         = .false.
    logical     ::         debug           = .false.
    logical     ::         showMissing     = .false.
    logical     ::         ignoreBadChunks = .false.
    integer     ::         numDiffs = 0
  end type Options_T

  type ( Options_T )              :: options
  character(len=255)              :: filename          ! filename
  character(len=255)              :: HEFfilename          ! filename
  character(len=255)              :: NC4filename          ! filename
  integer                         :: n_filenames
  integer                         :: error ! Counting indices & Error flags
  logical                         :: is_hdf5
  logical                         :: is_present
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
  MLSMessageConfig%CrashOnAnyError = .true.
  call mls_h5open(error)
  n_filenames = 0
  do      ! Loop over filenames
     call get_filename( filename, options )
     if ( filename == ' ' ) exit
     n_filenames = n_filenames + 1
     call h5fis_hdf5_f(trim(filename), is_hdf5, error)
     if ( .not. is_hdf5 ) then
       print *, 'Sorry--not recognized as hdf5 file: ', trim(filename)
     endif
     if ( options%verbose ) then
       print *, 'Rewriting swaths in     ', trim(filename)
     endif
     if ( n_filenames == 1 ) HEFFileName = filename
     if ( n_filenames == 2 ) then
       NC4FileName = filename
       exit
     endif
!      call OutputNamedValue( 'HDFEOS5 L2GP File',  trim(filename), &
!        & options='--Headline' )
!      call diff_two_files( trim(filename), options )
!      call resumeOutput
  enddo
  if ( options%ignoreBadChunks ) options%dumpOptions = trim(options%dumpOptions) // 'i'
  if ( options%silent ) options%dumpOptions = trim(options%dumpOptions) // 'm'
  if ( options%rms ) options%dumpOptions = trim(options%dumpOptions) // 'r'
  if ( options%stats ) options%dumpOptions = trim(options%dumpOptions) // 's'
  if ( options%table ) options%dumpOptions = trim(options%dumpOptions) // 'b'
  if ( options%AuBrick ) options%dumpOptions = trim(options%dumpOptions) // '@'
  if ( options%verbose ) options%dumpOptions = trim(options%dumpOptions) // 'v'
  call diff_two_files ( HEFFileName, NC4FileName, options )
  call mls_h5close(error)
contains
!------------------------- get_filename ---------------------
    subroutine get_filename( filename, options )
    ! Added for command-line processing
     character(len=255), intent(out)   :: filename
     type ( Options_T ), intent(inout) :: options
     ! Local variables
     integer ::                         error = 1
     integer, save ::                   i = 1
     character(LEN=160)              :: Chars
  ! Get inputfile name, process command-line args
  ! (which always start with -)
    do
      call getarg ( i+hp, filename )
      ! print *, i, ' th Arg: ', trim(filename)
      error = 0
      if ( filename(1:1) /= '-' ) exit
      if ( filename(1:3) == '-h ' ) then
        call print_help
      elseif ( filename(1:6) == '-force' ) then
        options%force = .true.
        exit
      elseif ( filename(1:8) == '-silent ' ) then
        options%silent = .true.
        exit
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
        exit
      else if ( filename(1:3) == '-f ' ) then
        call getarg ( i+1+hp, filename )
        i = i + 1
        exit
      else if ( filename(1:3) == '-D ' ) then
        call getarg ( i+1+hp, Chars )
        read(Chars, *) options%Details
        i = i + 1
        exit
      else if ( filename(1:3) == '-d ' ) then
        call getarg ( i+1+hp, options%dumpOptions )
        if ( index( options%dumpOptions, '?' ) > 0 ) then
          call DumpDumpOptions( "?" )
          stop
        endif
        i = i + 1
        exit
      else if ( filename(1:4) == '-deb' ) then
        options%debug = .true.
        exit
      else if ( filename(1:4) == '-ign' ) then
        options%ignorebadchunks = .true.
        exit
      else if ( lowercase(filename(1:3)) == '-au' ) then
        options%AuBrick = .true.
        exit
      else if ( filename(1:5) == '-rms ' ) then
        options%rms = .true.
        exit
      else if ( filename(1:3) == '-s ' ) then
        options%stats = .true.
        exit
      else if ( filename(1:2) == '-t' ) then
        options%table = .true.
        exit
      else if ( filename(1:6) == '-miss ' ) then
        options%showMissing = .true.
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
      & 'Usage: l2gpdiffnc4 [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options:'
      write (*,*) ' -f filename => use filename'
      write (*,*) ' -D details  => level of detail'
      write (*,*) ' -d options  => options string passed to diff routines'
      write (*,*) ' -h          => print brief help'
      write (*,*) ' -ls         => dump only swath names'
      write (*,*) ' -deb        => debug'
      write (*,*) ' -v          => verbose'
      write (*,*) ' (details level)'
      write (*,*) ' -0          => diff only scalars, 1-d array'
      write (*,*) ' -1          => diff only scalars'
      write (*,*) ' -2          => diff only swath names (same as -ls)'

      stop
  end subroutine print_help
  
  function IsAttributeInFile( file, attribute ) result(sooDesu)
    use MLSHDF5, only: IsHDF5itempresent
    use HDF5, only: H5fopen_F, H5f_Acc_Rdonly_F
    ! Dummy args
    character(len=*), intent(in) :: file
    character(len=*), intent(in) :: attribute
    logical :: sooDesu
    ! Local variables
    integer :: fileID
    integer :: grpID
    integer :: status
    character(len=len(attribute)) :: path, name
    ! TRUE if attribute in file
    call h5fopen_f ( trim(file), H5F_ACC_RDONLY_F, fileID, status )
    if ( status /= 0 ) call defeat('Opening file')
    call split_path_name ( attribute, path, name )
    call h5gopen_f( fileID, trim(path), grpID, status )
    if ( status /= 0 ) call defeat('Opening group')
    sooDesu = IsHDF5ItemPresent ( grpID, name, '-a' )
    call h5gclose_f(grpID, status)
    if ( status /= 0 ) call defeat('Closing group')
    call h5fclose_f(fileID, status)
    if ( status /= 0 ) call defeat('Closing file')
  end function IsAttributeInFile

  function IsDSInFile( file, DS ) result(sooDesu)
    use MLSHDF5, only: IsHDF5itempresent
    use HDF5, only: H5fopen_F, H5f_Acc_Rdonly_F
    ! Dummy args
    character(len=*), intent(in) :: file
    character(len=*), intent(in) :: DS
    logical :: sooDesu
    ! Local variables
    integer :: fileID
    integer :: status
    integer :: grpID
    character(len=len(DS)) :: path, name
    ! TRUE if DS in file
    call h5fopen_f ( trim(file), H5F_ACC_RDONLY_F, fileID, status )
    call split_path_name ( DS, path, name )
    call h5gopen_f( fileID, trim(path), grpID, status )
    sooDesu = IsHDF5ItemPresent ( grpID, name, '-d' )
    call h5gclose_f(grpID, status)
    call h5fclose_f(fileID, status)
  end function IsDSInFile
  
  subroutine Defeat(msg)
    character(len=*), intent(in) :: msg
    call resumeOutput
    call output('Serious error: ' // msg, advance='yes')
    call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'stopping' )
  end subroutine Defeat
  
   subroutine readHe5Swath ( filename, swath, l2gp )
   use L2GPData, only: L2GPData_T, L2GPnamelen, Maxswathnamesbufsize, Rgp, &
     & Dump, ReadL2GPData, DestroyL2GPcontents
     ! Dummy args
     character(len=*), intent(in)         :: filename ! filename
     character(len=*), intent(in)         :: swath ! swathname
     type (L2GPData_T)                    :: l2gp
     ! Local variables
       ! Allocate and fill l2gp
       ! print *, 'Reading swath from file: ', trim(swath)
       call ReadL2GPData ( trim(filename), trim(swath), l2gp, &
            & hdfVersion=HDFVERSION_5 )
   end subroutine readHe5Swath

   subroutine readNC4Swath ( filename, swath, l2gp )
   use L2GPData, only: L2GPData_T
   use NCL2GP, only: ReadNCL2GPData
     ! Dummy args
     character(len=*), intent(in)         :: filename ! filename
     character(len=*), intent(in)         :: swath ! swathname
     type (L2GPData_T)                    :: l2gp
     ! Local variables
       ! Allocate and fill l2gp
       ! print *, 'Reading swath from file: ', trim(swath)
       call ReadNCL2GPData ( trim(filename), trim(swath), l2gp )
   end subroutine readNC4Swath

   subroutine diff_two_files ( HEFFileName, NC4FileName, options )
     use L2GPData, only: L2GPData_T, Diff
     use PCFHdr, only:  GlobalAttributes_T, HE5_ReadGlobalAttr, &
    & DumpGlobalAttributes, GlobalAttributes
    ! Dummy args
    character(len=*), intent(in)         :: HEFFileName ! filename of hefeos5
    character(len=*), intent(in)         :: NC4FileName ! filename of NetCDF4
    type ( Options_T )                   :: options
    ! Local variables
    logical, parameter                   :: countEmpty = .true.
    integer :: File1
    integer                              :: i
    integer                              :: listsize
    type (L2GPData_T)                    :: he5l2gp
    type (L2GPData_T)                    :: nc4l2gp
    character (len=MAXSWATHNAMESBUFSIZE) :: matches
    type(MLSFile_T)                      :: MLSFile
    type(MLSFile_T)                      :: NC4File
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
    noSwaths = mls_InqSwath ( HEFFileName, SwathList, listSize, &
           & hdfVersion=HDFVERSION_5 )
    if ( options%verbose ) then
      call output('swaths in ' // trim(HEFFileName), advance='yes')
      call dump(trim(swathList))
      ! return
    endif
    ! Executable code
    status = InitializeMLSFile ( MLSFile, type=l_swath, access=DFACC_READ, &
     & name=HEFFileName, HDFVersion=HDFVERSION_5 )
    status = InitializeMLSFile ( NC4File, type=l_netcdf4, access=DFACC_READ, &
     & name=nc4filename )
     call Dump ( MLSFile )
     call Dump ( NC4File )
    noSwaths = NumStringElements(trim(swathList), countEmpty)
    if ( noSwaths < 1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'No swaths to dump--unable to count swaths in ' // trim(swathList) )
    endif
    ! Loop over swaths in file 1
    do i = 1, noSwaths
      call GetStringElement (trim(swathList), swath, i, countEmpty )
      ! Allocate and fill l2gp
      print *, 'Reading He5 swath from file: ', trim(HEFFileName)
      call readHe5Swath ( trim(HEFFileName), trim(swath), He5L2GP )
      print *, 'Min, Max (he5 swath) ', minval(he5L2GP%l2gpvalue), maxval(he5L2GP%l2gpvalue)
      print *, 'shape (he5 swath) ', shape(he5L2GP%l2gpvalue)
      print *, 'Reading NC4 swath from file: ', trim(NC4FileName)
      call readNC4Swath ( trim(NC4filename), trim(swath), NC4L2GP )
      print *, 'Min, Max (nc4 swath) ', minval(NC4L2GP%l2gpvalue), maxval(NC4L2GP%l2gpvalue)
      print *, 'shape (nc4 swath) ', shape(NC4L2GP%l2gpvalue)
      ! Diff the swaths
      call output ( 'About to diff the l2gps ', advance='yes' )
      call Diff( HE5l2gp, NC4l2gp, &
        & details=options%Details, &
        & options=options%dumpOptions )
      call DestroyL2GPContents ( HE5l2gp )
      call DestroyL2GPContents ( NC4l2gp )
    enddo
   end subroutine diff_two_files

!==================
end program l2gpdiffnc4
!==================

! $Log$
