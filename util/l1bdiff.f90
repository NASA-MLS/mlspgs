! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=================================
program l1bdiff ! diffs two l1b or L2AUX files
!=================================

   use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
   use Diff_1, only: Diff, SelfDiff
   use Dump_Options, only: DiffRMSMeansrms, RMSFormat, StatsOnOneLine, &
     & DumpDumpOptions
   use Dump_1, only: Dump
   use HDF, only: Dfacc_Create, Dfacc_Read
   use HDF5, only: H5fis_HDF5_F, &
     & H5gclose_F, H5gopen_F, H5gcreate_F
   use HighOutput, only: OutputnamedValue
   use L1bData, only: L1BData_T, Namelen, &
     & ContractL1BData, DeallocateL1BData, Diff, ReadL1BData
   use Machine, only: Hp, Getarg, NeverCrash
   use MLSFiles, only: FileNotFound, WildcardHDFVersion, &
     & MLS_Exists, MLS_HDF_Version, MLS_Sfstart, MLS_Sfend, &
     & HDFVersion_5
   use MLSHDF5, only: GetAllHDF5DSNames, MLS_H5open, MLS_H5close
   use MLSKinds, only: R8
   use MLSMessageModule, only: MLSMSG_Error, MLSMSG_Warning, &
     & MLSMessageConfig, MLSMessage
   use MLSStringLists, only: GetStringElement, NumStringElements
   use MLSStrings, only: Lowercase, Replace, Streq, WriteIntsToChars
   use Output_M, only: ResumeOutput, SuspendOutput, Output
   use Printit_M, only: Set_Config
   use Time_M, only: Time_Now, Time_Config

   implicit none

!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

! Brief description of program
! diffs two l1b or l2aux files

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"
! Then run it
! LF95.Linux/test [options] [input files]

  type Options_T
    logical     :: halfwaves = .false.
    logical     :: self = .false.
    logical     :: silent = .false.
    logical     :: timing = .false.
    logical     :: unique = .false.
    logical     :: debug = .false.
    logical     :: verbose = .false.
    logical     :: list = .false.
    logical     :: stats = .false.
    logical     :: AuBrick         = .false.
    logical     :: rms = .false.
    logical     :: table = .false.
    logical     :: direct = .true.
    logical     :: oneD = .true.
    logical     :: l2aux = .false.
    logical     :: ascii = .false. ! If true, diff even character fields
    integer     :: hdfVersion = HDFVERSION_5
    integer     :: maf1 = 0
    integer     :: maf2 = 0
    integer     :: moff = 0
    integer     :: numDiffs = 0
    character(len=255) :: group= ''     ! if SDs are within a group, group path
    character(len=255) :: sdList= '*'   ! what SDs to diff
    character(len=255) :: skipList= ''  ! what SDs to skip
    character(len=255) :: referenceFileName= 'default.h5'  ! reference filename
    character(len=80)  :: dumpOptions       = ' '
  end type Options_T

  type ( Options_T ) ::          options
  integer, parameter ::          MAXDS = 300
  integer, parameter ::          MAXSDNAMESBUFSIZE = MAXDS*namelen
  integer, parameter ::          MAXFILES = 100
  ! character(len=8)   ::          options%dumpOptions
  character(len=255) ::          filename          ! input filename
  character(len=255), dimension(MAXFILES) :: filenames
  integer     ::                 i, error ! Counting indices & Error flags
  logical     ::                 is_hdf5
  integer ::                     maf1, maf2
  character (len=MAXSDNAMESBUFSIZE) :: mySdList
  integer            ::          n_filenames
  character(len=16) ::           string
  real        ::                 t1
  real        ::                 t2
  real        ::                 tFile
  !
  call set_config ( useToolkit = .false., logFileUnit = -1 )
  time_config%use_wall_clock = .true.
  DIFFRMSMEANSRMS = .false.
  CALL mls_h5open(error)
  NeverCrash = .false.
  MLSMessageConfig%crashOnAnyError = .true.
  statsOnOneLine = .false.
  n_filenames = 0
  do      ! Loop over filenames
     call get_options(filename, n_filenames, options)
     if ( filename(1:1) == '-' ) cycle
     if ( filename == ' ' ) exit
     call h5fis_hdf5_f(trim(filename), is_hdf5, error)
     if ( options%hdfVersion==HDFVERSION_5 .and. .not. is_hdf5 ) then
       print *, 'Sorry--not recognized as hdf5 file: ', trim(filename)
       cycle
     endif
     n_filenames = n_filenames + 1
     filenames(n_filenames) = filename
  enddo
  if ( n_filenames == 0 ) then
    if ( options%verbose ) print *, 'Sorry no input files supplied'
    stop
  elseif ( options%verbose .and. options%silent ) then
    print *, 'Sorry--either verbose or silent; cant be both'
    stop
  elseif ( options%self ) then
    ! We'll do diffs within each file
  elseif ( options%referenceFileName == 'default.h5' ) then
    options%referenceFileName = filenames(n_filenames)
    n_filenames = n_filenames - 1
  endif
  if ( options%rms .or. options%AuBrick ) rmsFormat = '(1pe9.2)'
  if ( options%silent ) call suspendOutput
  ! options%dumpOptions = '-'
  if ( options%rms ) options%dumpOptions = trim(options%dumpOptions) // 'r'
  if ( options%stats ) options%dumpOptions = trim(options%dumpOptions) // 's'
  if ( options%unique ) options%dumpOptions = trim(options%dumpOptions) // 'u'
  if ( options%silent ) options%dumpOptions = trim(options%dumpOptions) // 'h'
  if ( options%direct ) options%dumpOptions = trim(options%dumpOptions) // 'd'
  if ( options%table ) options%dumpOptions = trim(options%dumpOptions) // 'b'
  if ( options%AuBrick ) options%dumpOptions = trim(options%dumpOptions) // '@'
  call time_now ( t1 )
  if ( options%verbose .and. .not. options%list ) &
    & print *, 'Compare l1b data to: ', trim(options%referenceFileName)
  do i=1, n_filenames
    call time_now ( tFile )
    if ( options%list ) then
      print *, 'DS Names in: ', trim(filenames(i))
      if ( len_trim(options%group) < 1 ) then
        call GetAllHDF5DSNames (trim(filenames(i)), '/', mysdList)
      else
        call GetAllHDF5DSNames (trim(filenames(i)), trim(options%group), mysdList)
      endif
      call dump(mysdList, 'DS names')
    elseif ( options%self )then
      if ( options%verbose ) then
        print *, 'diffing from: ', trim(filenames(i))
      endif
      call mySelfDiff(trim(filenames(i)), options%hdfVersion, options)
      if ( options%timing ) call sayTime('diffing this file', tFile)
    else
      if ( options%verbose ) then
        print *, 'diffing from: ', trim(filenames(i))
      endif
      call myDiff(trim(filenames(i)), &
      & trim(options%referenceFileName), &
      & options%hdfVersion, options)
      if ( options%timing ) call sayTime('diffing this file', tFile)
    endif
  enddo
  if ( options%timing ) call sayTime('diffing all files')
  call resumeOutput
  ! call outputnamedValue( 'num diffs', options%numDiffs )
  if ( options%silent .and. options%numDiffs > 0 ) then
    ! write(string, '(i)') options%numDiffs
    call WriteIntsToChars ( options%numDiffs, string )
    call print_string(string)
  endif
  call mls_h5close(error)
contains
!------------------------- get_options ---------------------
    subroutine get_options( filename, n_filenames, options )
    ! Added for command-line processing of options, filenames
     character(LEN=255), intent(out)       :: filename          ! filename
     integer, intent(in)                   :: n_filenames
     type ( Options_T ), intent(inout)     :: options
     ! Local variables
     integer ::                         error = 1
     integer, save ::                   i = 1
     character(LEN=16)       :: number
  ! Get inputfile name, process command-line args
  ! (which always start with -)
    do
      call getarg ( i+hp, filename )
      ! print *, i, ' th Arg: ', trim(filename)
      error = 0
      if ( filename(1:1) /= '-' ) exit
      if ( filename(1:3) == '-h ' ) then
        call print_help
      elseif ( filename(1:5) == '-maf ' ) then
        call getarg ( i+1+hp, number )
        read(number, *) options%maf1, options%maf2
        i = i + 1
        exit
      elseif ( filename(1:6) == '-moff ' ) then
        call getarg ( i+1+hp, number )
        read(number, *) options%moff
        i = i + 1
        exit
      elseif ( filename(1:3) == '-d ' ) then
        call getarg ( i+1+hp, options%sdList )
        i = i + 1
        exit
      elseif ( filename(1:3) == '-r ' ) then
        call getarg ( i+1+hp, options%referenceFileName )
        i = i + 1
        exit
      elseif ( filename(1:5) == '-half' ) then
        options%halfWaves = .true.
        exit
      else if ( filename(1:4) == '-hdf' ) then
        call getarg ( i+1+hp, number )
        read(number, *) options%hdfVersion
        i = i + 1
        exit
      elseif ( filename(1:4) == '-one' ) then
        statsOnOneLine = .true.
        exit
      elseif ( filename(1:5) == '-opt ' ) then
        call getarg ( i+1+hp, options%dumpOptions )
        if ( index( options%dumpOptions, '?' ) > 0 ) then
          call DumpDumpOptions( "?" )
          stop
        endif
        i = i + 1
        exit
      elseif ( filename(1:6) == '-self ' ) then
        options%self = .true.
        exit
      elseif ( filename(1:8) == '-silent ' ) then
        options%silent = .true.
        exit
      elseif ( filename(1:8) == '-unique ' ) then
        options%unique = .true.
        exit
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
        exit
      elseif ( filename(1:7) == '-ascii ' ) then
        options%ascii = .true.
        exit
      elseif ( filename(1:7) == '-l2aux ' ) then
        options%l2aux = .true.
        exit
      elseif ( filename(1:3) == '-l ' ) then
        options%list = .true.
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
      elseif ( filename(1:6) == '-skip ' ) then
        call getarg ( i+1+hp, options%skipList )
        i = i + 1
        exit
      else if ( filename(1:3) == '-f ' ) then
        call getarg ( i+1+hp, filename )
        i = i + 1
        exit
      else if ( filename(1:3) == '-g ' ) then
        call getarg ( i+1+hp, options%group )
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
      print *,  "Enter the name of the HDF5 l1b or l2aux file. "
      read(*,'(a)') filename
    endif

  end subroutine get_options
!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
    write (*,*) &
    & 'Usage:l1bdiff [options] [filenames]'
    write (*,*) &
    & ' If no filenames supplied, you will be prompted to supply one'
    write (*,*) ' Options: -f filename => add filename to list of filenames'
    write (*,*) '                  (can do the same w/o the -f)'
    write (*,*) '   -d list     => just diff the SDs in list'
    write (*,*) '   -g group    => find SDs under group path'
    write (*,*) '   -ascii      => diff even character-valued fields'
    write (*,*) '   -l2aux      => the files are l2aux, not l1b'
    write (*,*) '   -v          => switch on verbose mode'
    write (*,*) '   -self       => dump successive differences'
    write (*,*) '                  between values in same file'
    write (*,*) '   -half       => (1) (if with -self) '
    write (*,*) '                  show no. of 1/2 waves'
    write (*,*) '                  (2) (otherwise)'
    write (*,*) '                  diff 1/2 of channels (so DACS wont crash)'
    write (*,*) '   -hdf version=> hdf version (default is 5)'
    write (*,*) '   -one        => print name on each line (dont)'
    write (*,*) '   -opt opts   => pass opts to dump routines'
    write (*,*) '                  e.g., "?" to list available ones'
    write (*,*) '   -silent     => switch on silent mode'
    write (*,*) '                 (printing only if diffs found)'
    write (*,*) '   -unique     => dump only unique elements'
    write (*,*) '   -l          => just list sd names in files'
    write (*,*) '   -maf m1,m2  => just diff in the range [m1,m2]'
    write (*,*) '   -moff offset=> 2nd data set starts after 1st'
    write (*,*) '   -au         => format like goldbrick'
    write (*,*) '   -rms        => just print mean, rms'
    write (*,*) '   -s          => just show number of differences'
    write (*,*) '   -skip list  => skip diffing the SDs in list'
    write (*,*) '   -t[able]    => table of % vs. amount of differences (pdf)'
    write (*,*) '   -h          => print brief help'
    stop
  end subroutine print_help
!------------------------- SayTime ---------------------
  subroutine SayTime ( What, startTime )
    character(len=*), intent(in) :: What
    real, intent(in), optional :: startTime
    real :: myt1
    if ( present(startTime) ) then
      myt1 = startTime
    else
      myt1 = t1
    endif
    call time_now ( t2 )
    call output ( "Timing for " // what // " = " )
    call output ( dble(t2 - myt1), advance = 'yes' )
  end subroutine SayTime

  ! ---------------------- myDiff  ---------------------------
  subroutine myDiff( file1, file2, hdfVersion, options )
  !------------------------------------------------------------------------

    ! Given file names file1 and file2,
    ! This routine prints the diffs between them

    ! Arguments

    character (len=*), intent(in) :: file1 ! Name of file 1
    character (len=*), intent(in) :: file2 ! Name of file 2
    integer, intent(in)           :: hdfVersion
    type ( Options_T )            :: options

    ! Local
    logical, parameter            :: countEmpty = .true.
    integer :: data_type
    integer, dimension(7) :: dimsizes
    logical :: file_exists
    integer :: file_access
    integer :: firstChannel
    integer :: grpid
    integer :: i
    integer :: iHalf
    logical :: isl1boa
    type(l1bdata_t) :: L1BDATA  ! Result
    type(l1bdata_t) :: L1BDATA2 ! for diff
    type(l1bdata_t) :: L1BDATAT
    real(r8), dimension(:), pointer :: l1bValues1 => null()
    real(r8), dimension(:), pointer :: l1bValues2 => null()
    integer :: lastChannel
    logical :: mustDiff
    character(len=8) :: myOptions
    character (len=MAXSDNAMESBUFSIZE) :: mySdList
    integer :: nHalves
    integer :: NoMAFs
    integer :: noSds
    integer :: nsize
    integer :: num_Attrs
    integer :: numDiffs
    integer :: rank
    integer :: sdfid1
    integer :: sdfid2
    character (len=80) :: sdName
    integer :: sds_id
    real :: stime
    integer :: status
    integer :: the_hdfVersion

    ! hdf4 externals
    integer :: sfstart, sffinfo, sfselect, sfginfo, sfend
    external :: sfstart, sffinfo, sfselect, sfginfo, sfend

    ! Executable code
    stime = t2
    ! the_hdfVersion = HDFVERSION_5
    the_hdfVersion = hdfVersion
    file_exists = ( mls_exists(trim(File1)) == 0 )
    if ( .not. file_exists ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'File 1 not found; make sure the name and path are correct' &
        & // trim(file1) )
    endif
    if ( the_hdfVersion == WILDCARDHDFVERSION ) then
      the_hdfVersion = mls_hdf_version(File1, hdfVersion)
      if ( the_hdfVersion == FILENOTFOUND ) &
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'File 1 not found; make sure the name and path are correct' &
          & // trim(file1) )
    endif
    if ( options%verbose ) print *, 'the hdf version: ', the_hdfVersion
    if ( the_hdfVersion == HDFVERSION_5 ) then
      if ( len_trim(options%group) < 1 ) then
        call GetAllHDF5DSNames (trim(File1), '/', mysdList)
      else
        call GetAllHDF5DSNames (trim(File1), trim(options%group), mysdList)
      endif
      if ( options%verbose ) then
        call output ( '============ DS names in ', advance='no' )
        call output ( trim(file1) //' ============', advance='yes' )
      endif
      if ( mysdList == ' ' .and. .not. options%silent ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'No way yet to find sdList in ' // trim(File1) )
        return
      else
        if ( options%verbose ) call dump(mysdList, 'DS names')
      endif
      if ( options%sdList /= '*' ) then
        call output( ' Selected SDs to diff', advance='yes' )
        mySDList = options%sdList
        call dump(mysdList, 'DS names')
      endif

      isl1boa = (index(trim(mysdList), '/GHz') > 0)
      file_access = DFACC_READ
      sdfid1 = mls_sfstart(File1, DFACC_READ, hdfVersion=hdfVersion)
      if (sdfid1 == -1 ) then
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        &  'Failed to open l1b file ' // trim(File1) )
      end if
	   call h5gOpen_f (sdfid1,'/', grpID, status)
      if ( status /= 0 .and. .not. options%silent ) then
	     call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	  & 'Unable to open group to read attribute in l2aux file' )
      endif
      sdfId2 = mls_sfstart(trim(file2), file_access, &
                & hdfVersion=hdfVersion)
      if (sdfid2 == -1 ) then
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Failed to open l1b ' // trim(File2) )
      end if
      noSds = NumStringElements(trim(mysdList), countEmpty)
      if ( noSds < 1 .and. .not. options%silent ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'No sdNames cp to file--unable to count sdNames in ' // trim(mysdList) )
      endif
      ! May need to create groups '/sc', '/GHz', and '/THz'
      if ( isL1boa .and. file_access == DFACC_CREATE ) then
        call h5gcreate_f(sdfId2, '/sc', grpID, status)
        if ( status /= 0 ) then
	       call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	  & 'Unable to create group /sc in ' // trim(File2) )
        endif
        call h5gcreate_f(sdfId2, '/GHz', grpID, status)
        if ( status /= 0 ) then
	       call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	  & 'Unable to create group /GHz in ' // trim(File2) )
        endif
        call h5gcreate_f(sdfId2, '/THz', grpID, status)
        if ( status /= 0 ) then
	       call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	  & 'Unable to create group /THz in ' // trim(File2) )
        endif
      endif
    else
      sdfid1 = sfstart( trim(File1), DFACC_READ )
      status = sffinfo( sdfid1, noSds, nsize )
      mysdList = '*'
      sdfid2 = sfstart( trim(File2), DFACC_READ )
      status = sffinfo( sdfid1, noSds, nsize )
    endif
    if ( options%halfwaves ) then
      nHalves = 2
    else
      nHalves = 1
    endif
    ! Loop over sdNames in file 1
    ! (But skip PCF and HDFEOS INFORMATION/coremetadata.0)
    do i = 1, noSds
      numDiffs = 0
      if ( the_hdfversion == HDFVERSION_5 ) then
        call GetStringElement (trim(mysdList), sdName, i, countEmpty )
        if ( len_trim(options%group) > 1 ) sdName = trim(options%group) // sdName
      else
        sds_id = sfselect( sdfid1, i-1 ) ! "c" arrays begin with index 0
        status = sfginfo( sds_id, sdName, rank, dimsizes, data_type, num_attrs )
        ! print *, 'sdName: ', sdName
      endif
      if ( any( &
        & streq( &
        & (/ 'PCF ', 'meta', 'l2cf', 'utcp', 'leap', 'LCF ' /), &
        & sdname, options='-Pw' ) ) .or. &
        &  index(options%skipList, trim(sdName)) > 0 ) cycle
      do ihalf=1, nHalves
        ! Allocate and fill l2aux
        if ( options%debug ) print *, 'About to read ', trim(sdName)
        if ( options%halfwaves ) then
          call ReadL1BData ( sdfid1, trim(sdName), L1bDataT, NoMAFs, status, &
            & hdfVersion=the_hdfVersion, NEVERFAIL=.true., l2aux=options%l2aux )
          if ( noMAFs < 1 ) cycle
          if ( status /= 0 .and. .not. options%silent ) then
	         call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	  & 'Unable to find ' // trim(sdName) // ' in ' // trim(File1) )
            call DeallocateL1BData ( l1bData )
            cycle
          endif
          if ( iHalf == 1 ) then
            firstChannel = 1
            lastChannel = min(L1bDataT%NoAuxInds, L1bDataT%NoAuxInds/2 + 1)
          else
            if ( L1bDataT%NoAuxInds < 3 ) cycle
            firstChannel = L1bDataT%NoAuxInds/2 + 2
            lastChannel = L1bDataT%NoAuxInds
          endif
          call contractL1BData( L1bDataT, L1bData, noMAFs, &
            & firstChannel=firstChannel, lastChannel=lastChannel)
          call deallocateL1BData( L1bDataT )
        else
          call ReadL1BData ( sdfid1, trim(sdName), L1bData, NoMAFs, status, &
            & hdfVersion=the_hdfVersion, NEVERFAIL=.true., l2aux=options%l2aux )
          if ( status /= 0 .and. .not. options%silent ) then
	         call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	  & 'Unable to find ' // trim(sdName) // ' in ' // trim(File1) )
            call DeallocateL1BData ( l1bData )
            cycle
          endif
        endif
        if ( options%timing ) call SayTime( 'Reading l1bdata 1', stime )
        stime = t2
        if ( options%halfwaves ) then
          call ReadL1BData ( sdfid1, trim(sdName), L1bDataT, NoMAFs, status, &
            & hdfVersion=the_hdfVersion, NEVERFAIL=.true., l2aux=options%l2aux )
          call contractL1BData( L1bDataT, L1bData2, noMAFs, &
            & firstChannel=firstChannel, lastChannel=lastChannel )
          call deallocateL1BData( L1bDataT )
        else
          call ReadL1BData ( sdfid2, trim(sdName), L1bData2, NoMAFs, status, &
            & hdfVersion=the_hdfVersion, NEVERFAIL=.true., l2aux=options%l2aux )
          if ( status /= 0 .and. .not. options%silent ) then
	         call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	  & 'Unable to find ' // trim(sdName) // ' in ' // trim(File1) )
            call DeallocateL1BData ( l1bData )
            call DeallocateL1BData ( l1bData2 )
            cycle
          endif
        endif
        if ( options%timing ) call SayTime( 'Reading l1bdata 2', stime )
        stime = t2
        if ( associated(L1bData%charField) .and. .not. options%ascii .and. &
          & .not. options%silent ) then
	       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Skipping diff of char-valued ' // trim(sdName) )
          call DeallocateL1BData ( l1bData )
          call DeallocateL1BData ( l1bData2 )
          cycle
        endif
        maf1 = 1
        if ( options%maf1 > 0 ) maf1 = options%maf1
        maf2 = NoMAFs
        if ( options%maf2 > 0 ) maf2 = options%maf2
        if ( options%verbose ) then
          print *, 'About to diff :: ' // trim(sdName) // ' ::'
        elseif ( .not. options%silent ) then
          print *, ':: ' // trim(sdName) // ' ::'
        endif
        mustDiff = .true.
        if ( associated(L1bData%dpField) .and. associated(L1bData2%dpField) ) then
          if ( all(L1bData%dpField == L1bData2%dpField) ) mustDiff = .false.
        endif
        myOptions = options%dumpOptions
        if ( mustDiff ) myOptions = Replace ( myOptions, 'h', ' ' )
        if ( associated(L1bData%dpField) .and. mustDiff .and. options%debug ) &
          & call outputNamedValue( 'num diffs', &
          & count(L1bData%dpField /= L1bData2%dpField), &
          & advance='yes' )
        if ( .not. mustDiff ) then
          if ( .not. options%silent ) &
            & call output( '(The two arrays are exactly equal)', advance='yes' )
        elseif ( options%oneD .and. associated(L1bData%dpField) ) then
          ! We will store L1BData%dpField in a values array
          nsize = product(shape(L1bData%dpField))
          call dump( L1bData%dpField-L1bData2%dpField, &
            & options=myOptions )
          ! We use dump instead of diff (but why?)
          ! stop
          call allocate_test(l1bValues1, nsize, 'l1bValues1', ModuleName )
          l1bValues1 = reshape( L1bData%dpField, (/nsize/) )
          call deallocate_test( L1bData%dpField, 'l1bData%Values1', ModuleName )
          call allocate_test(l1bValues2, nsize, 'l1bValues2', ModuleName )
          l1bValues2 = reshape( L1bData2%dpField, (/nsize/) )
          call deallocate_test( L1bData2%dpField, 'l1bData%Values2', ModuleName )
          if ( options%timing ) call SayTime( 'Copying values to 1-d arrays', stime )
          stime = t2

          call diff(L1bData, L1bData2, numDiffs=numDiffs, options=myOptions, &
            & l1bValues1=l1bValues1, l1bValues2=l1bValues2 )
          call deallocate_test( L1bValues1, 'l1bValues1', ModuleName )
          call deallocate_test( L1bValues2, 'l1bValues2', ModuleName )
        elseif ( options%direct .or. .not. associated(L1bData%dpField) ) then
          if ( options%verbose ) print *, 'About to do it direct ' // myOptions
          call diff(L1bData, L1bData2, numDiffs=numDiffs, options=myOptions )
        else
          call diff(L1bData, L1bData2, details=0, &
            & numDiffs=numDiffs, options=myOptions )
          if ( .true. .and. associated(L1bData%dpField) .and. &
            & associated(L1bData2%dpField)) then
            if ( options%silent .and. numDiffs < 1 ) then
            elseif ( .false. .and. all(L1bData%dpField(:,:,maf1:maf2) == &
              & L1bData2%dpField(:,:,maf1+options%moff:maf2+options%moff)) ) then
              print *, '(The two fields are exactly equal)'
            elseif( options%moff /= 0 .or. options%maf1 > 0 ) then
              call diff( L1bData%dpField(:,:,maf1:maf2), &
                & '(1)', &
                & L1bData2%dpField(:,:,maf1+options%moff:maf2+options%moff), &
                & '(2)', &
                & options=myOptions )
            else
              call diff( L1bData%dpField, &
                & '(1)', &
                & L1bData2%dpField, &
                & '(2)', &
                & options=myOptions )
            endif
          else
            L1bData%dpField = L1bData%dpField - L1bData2%dpField(:,:,1+options%moff:)
            if ( options%silent ) then
            elseif ( .false. .and. all(L1bData%dpField == 0.d0) ) then
              print *, '(The two fields are exactly equal)'
            elseif ( options%maf1 /= 0 .and. options%maf2 /= 0 ) then
              call dump( L1bData%dpField(:,:,options%maf1:options%maf2), &
                & 'l1bData%dpField', options=myOptions )
            else
              call dump( L1bData%dpField, &
                & 'l1bData%dpField', options=myOptions )
            endif
          endif
        endif
        call DeallocateL1BData ( l1bData )
        call DeallocateL1BData ( l1bData2 )
        options%numDiffs = options%numDiffs + numDiffs
        if ( options%timing ) call SayTime( 'Doing the diff', stime )
        stime = t2
      enddo ! Loop of halves
      if ( options%debug ) print *, 'After ' // trim(sdName) // ' ', options%numDiffs
    enddo ! Loop of datasets
    if ( the_hdfVersion == HDFVERSION_5) call h5gClose_f (grpID, status)
    if ( status /= 0 .and. .not. options%silent ) then
	   call MLSMessage ( MLSMSG_Warning, ModuleName, &
       & 'Unable to close group in l2aux file: ' // trim(File1) // ' after diffing' )
    endif
    if ( the_hdfVersion == HDFVERSION_5 ) then
	   status = mls_sfend(sdfid1, hdfVersion=the_hdfVersion)
      if ( status /= 0 ) &
        call MLSMessage ( MLSMSG_Error, ModuleName, &
         & "Unable to close L2aux file: " // trim(File1) // ' after diffing')
	   status = mls_sfend(sdfid2, hdfVersion=the_hdfVersion)
      if ( status /= 0 ) &
        call MLSMessage ( MLSMSG_Error, ModuleName, &
         & "Unable to close L2aux file: " // trim(File2) // ' after diffing')
    else
      ! print *, 'Ending MyDiff'
      status = sfend( sdfid1 )
      status = sfend( sdfid2 )
    endif
  end subroutine myDiff

  ! ---------------------- mySelfDiff  ---------------------------
  subroutine mySelfDiff( file1, hdfVersion, options )
  !------------------------------------------------------------------------

    ! Given file name file1
    ! This routine prints its selfdiffs

    ! Arguments

    character (len=*), intent(in) :: file1 ! Name of file 1
    integer, intent(in)           :: hdfVersion
    type ( Options_T )            :: options

    ! Local
    logical, parameter            :: countEmpty = .true.
    logical :: file_exists
    integer :: grpid
    integer :: i
    type(l1bdata_t) :: L1BDATA  ! Result
    character (len=MAXSDNAMESBUFSIZE) :: mySdList
    integer :: NoMAFs
    integer :: noSds
    integer :: numDiffs
    integer :: sdfid1
    character (len=80) :: sdName
    integer :: status
    integer :: the_hdfVersion
    logical, parameter :: SKIPIFRANKTOOHIGH = .false.

    ! Executable code
    ! the_hdfVersion = HDFVERSION_5
    the_hdfVersion = hdfVersion
    file_exists = ( mls_exists(trim(File1)) == 0 )
    if ( .not. file_exists ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'File 1 not found; make sure the name and path are correct' &
        & // trim(file1) )
    endif
    if ( the_hdfVersion == WILDCARDHDFVERSION ) then
      the_hdfVersion = mls_hdf_version(File1, hdfVersion)
      if ( the_hdfVersion == FILENOTFOUND ) &
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'File 1 not found; make sure the name and path are correct' &
          & // trim(file1) )
    endif
    call GetAllHDF5DSNames (trim(File1), '/', mysdList)
    if ( options%verbose) then
      call output ( '============ DS names in ', advance='no' )
      call output ( trim(file1) //' ============', advance='yes' )
    endif
    if ( mysdList == ' ' .and. .not. options%silent ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'No way yet to find sdList in ' // trim(File1) )
      return
    else
      if ( options%verbose ) call dump(mysdList, 'DS names')
    endif
    if ( options%sdList /= '*' ) then
      call output( ' Selected SDs to diff', advance='yes' )
      mySDList = options%sdList
      call dump(mysdList, 'DS names')
    endif

    sdfid1 = mls_sfstart(File1, DFACC_READ, hdfVersion=hdfVersion)
    if (sdfid1 == -1 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      &  'Failed to open l1b file ' // trim(File1) )
    end if
	 call h5gOpen_f (sdfid1,'/', grpID, status)
    if ( status /= 0 .and. .not. options%silent ) then
	   call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	& 'Unable to open group to read attribute in l2aux file' )
    endif
    noSds = NumStringElements(trim(mysdList), countEmpty)
    if ( noSds < 1 .and. .not. options%silent ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'No sdNames cp to file--unable to count sdNames in ' // trim(mysdList) )
    endif
    ! Loop over sdNames in file 1
    ! (But skip PCF and HDFEOS INFORMATION/coremetadata.0)
    do i = 1, noSds
      call GetStringElement (trim(mysdList), sdName, i, countEmpty )
      if ( sdName == 'PCF' .or. &
        &  sdName == 'HDFEOS INFORMATION/coremetadata.0' .or. &
        &  sdName == 'l2cf' .or. &
        &  index(options%skipList, trim(sdName)) > 0 ) cycle
      ! Allocate and fill l2aux
      ! if ( options%verbose ) print *, 'About to read ', trim(sdName)
      call ReadL1BData ( sdfid1, trim(sdName), L1bData, NoMAFs, status, &
        & hdfVersion=the_hdfVersion, NEVERFAIL=.true., l2aux=options%l2aux )
      if ( status /= 0 .and. .not. options%silent ) then
	     call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Unable to find ' // trim(sdName) // ' in ' // trim(File1) )
        call DeallocateL1BData ( l1bData )
        cycle
      endif
      if ( associated(L1bData%charField) .and. .not. options%silent ) then
	     call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Skipping diff of char-valued ' // trim(sdName) )
      elseif ( associated(L1bData%intField) .and. .not. options%silent ) then
        if ( ( size(L1bData%intField, 1) > 1 .or. &
          &  size(L1bData%intField, 2) > 1 ) .and. SKIPIFRANKTOOHIGH ) then
	       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Skipping diff high-rank ' // trim(sdName) )
        else
          call selfdiff( L1bData%intField(1, 1, :), trim(sdName), &
            & waves=options%halfWaves, options=options%dumpOptions )
        endif
      elseif ( associated(L1bData%dpField) ) then
        if ( ( size(L1bData%dpField, 1) > 1 .or. &
          &  size(L1bData%dpField, 2) > 1 ) .and. SKIPIFRANKTOOHIGH &
          & .and. .not. options%silent ) then
	       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Skipping diff high-rank ' // trim(sdName) )
        else
          call selfdiff( L1bData%dpField(1, 1, :), trim(sdName), &
            & waves=options%halfWaves, options=options%dumpOptions )
        endif
      endif
      call DeallocateL1BData ( l1bData )
    enddo
    call h5gClose_f (grpID, status)
    if ( status /= 0 .and. .not. options%silent ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & 'Unable to close group in l2aux file: ' // trim(File1) // ' after diffing' )
    endif
	 status = mls_sfend(sdfid1, hdfVersion=the_hdfVersion)
    if ( status /= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Unable to close L2aux file: " // trim(File1) // ' after diffing')
  end subroutine mySelfDiff

!------------------------- print_string ---------------------
  subroutine print_string(string)
    character(len=*), intent(in) :: string
    write(*,'(a)') trim(string)
  end subroutine print_string

!==================
end program l1bdiff
!==================

! $Log$
! Revision 1.37  2017/10/20 20:20:08  pwagner
! Crash with a walkback if an error occurs
!
! Revision 1.36  2017/10/12 18:58:27  pwagner
! CamelCase more use statements
!
! Revision 1.35  2016/10/05 20:14:53  pwagner
! Implemented Au (Gold) option
!
! Revision 1.34  2016/08/09 22:41:40  pwagner
! Consistent with splitting of Dunp_0
!
! Revision 1.33  2016/07/28 01:46:38  vsnyder
! Refactor diff and dump
!
! Revision 1.32  2016/03/23 16:38:28  pwagner
! Added -one commandline option
!
! Revision 1.31  2015/06/30 18:50:33  pwagner
! Corrected error in passing dumpOptions to dump module; -opt '?' now dumps available dump options
!
! Revision 1.30  2014/03/07 21:43:27  pwagner
! Name_Len changed to nameLen; got from MLSCommon
!
! Revision 1.29  2014/01/09 00:31:26  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 1.28  2013/08/23 02:51:47  vsnyder
! Move PrintItOut to PrintIt_m
!
! Revision 1.27  2013/01/09 18:48:05  pwagner
! Dont omit sdname from output
!
! Revision 1.26  2012/09/12 16:40:13  pwagner
! Works more reliably with goldbrick
!
! Revision 1.25  2012/09/06 00:37:02  pwagner
! Works better with goldbrick
!
! Revision 1.24  2012/06/14 00:01:16  pwagner
! Willing to selfdiff higher-rank fields instance-wise
!
! Revision 1.23  2012/04/20 20:55:48  pwagner
! Restored silence to silent option
!
! Revision 1.22  2012/04/20 17:57:51  pwagner
! Fixed syntax error of dumpOptions
!
! Revision 1.21  2012/02/13 23:41:31  pwagner
! -opt opts passes opts to underlying dump routines
!
! Revision 1.20  2011/05/10 18:28:12  pwagner
! Avoid annoying messages unless verbose
!
! Revision 1.19  2010/07/23 17:51:26  pwagner
! Now able to diff hdf4-formatted files
!
! Revision 1.18  2009/11/20 23:00:50  pwagner
! Should not segfault when diffing DACS datasets
!
! Revision 1.17  2009/11/02 19:54:06  pwagner
! Fixed bug causing goldbrick to print excessively
!
! Revision 1.16  2009/10/30 23:06:34  pwagner
! Last change let executable land on bare stop; must fix oneD yet
!
! Revision 1.15  2009/09/11 23:24:01  pwagner
! options now include -t to be consistent with diff apis
!
! Revision 1.14  2009/08/18 20:43:23  pwagner
! Changes needed to diff DACS radiance files
!
! Revision 1.13  2009/08/04 20:45:12  pwagner
! silent option adds 'h' to options arg
!
! Revision 1.12  2009/06/16 22:37:40  pwagner
! Changed api for dump, diff routines; now rely on options for most optional behavior
!
! Revision 1.11  2008/07/10 00:14:43  pwagner
! Added -self, -half, -unique options
!
! Revision 1.10  2008/04/10 20:22:37  pwagner
! Less voluminous output
!
! Revision 1.9  2008/02/28 01:34:07  pwagner
! Normally should skip diffing character-valued fields
!
! Revision 1.8  2007/11/28 19:35:21  pwagner
! May specify that files are l2aux
!
! Revision 1.7  2007/07/18 00:15:25  pwagner
! -rms outputs fractional diffs more clearly
!
! Revision 1.6  2006/11/22 18:33:10  pwagner
! New optional args to diff l1b files with different number MAFs
!
! Revision 1.5  2006/06/14 16:42:38  pwagner
! Should not run out of memory unless direct reset to TRUE
!
! Revision 1.4  2006/01/14 00:58:31  pwagner
! Added -silent option
!
! Revision 1.3  2005/10/29 00:13:56  pwagner
! Removed unused procedures from use statements
!
! Revision 1.2  2005/09/23 21:01:13  pwagner
! use_wall_clock now a component of time_config
!
! Revision 1.1  2005/06/22 19:27:32  pwagner
! Reworded Copyright statement, moved rcs id
!
