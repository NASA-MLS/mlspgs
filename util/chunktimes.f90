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
program chunktimes ! Reads chunk times from l2aux file(s)
!=================================

  use Dump_0, only: Dump
  use Dump_1, only: Dump
  use HDF, only: Dfacc_Rdonly
  use HDF5, only: Hsize_T, H5FIs_HDF5_F, H5GClose_F, H5GOpen_F
  use HighOutput, only: Output_Date_And_Time
  use L1BData, only: Namelen
  use Machine, only: Hp, Getarg
  use MLSKinds, only: R4, R8
  use MLSFiles, only: MLS_Exists, MLS_SFStart, MLS_SFEnd, &
    & HDFVersion_4, HDFVersion_5
  use MLSHDF5, only: &
    & GetAllHDF5DSNames, GetHDF5Attribute, GetHDF5DSDims, &
    & ISHDf5AttributePresent, LoadFromHDF5DS, MLS_H5Open, MLS_H5Close
  use MLSFinds, only: Findall, Findfirst, Findlast, Findnext
  use MLSStats1, only: Stat_T, Statsononeline, &
    & DumpStat=>dump, Statistics
  use MLSStringLists, only: CatLists, GetStringElement, GetUniqueList, &
    & NUMStringElements, StringElement, StringElementNum
  use MLSStrings, only: Lowercase
  use Output_M, only: OutputOptions, &
    & Blanks, Newline, Output
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
! Reads chunk times and failures from list of input files

  type options_T
    logical            :: verbose = .false.
    logical            :: merge = .false.           ! Merge data from files
    logical            :: oneLine = .false.         ! Print on one line
    logical            :: showFailed = .false.      ! show howmany, which failed
    logical            :: showStats = .true.        ! show max, min, mean, etc.
    logical            :: showQManager = .false.    ! show QManager performance
    logical            :: guessFinalPhase = .true.  ! guess how many phases
    logical            :: showWhereFailed = .false. ! show where chunks failed
    character(len=255) :: DSName= 'phase timing'    ! Dataset name
    character(len=255) :: binopts= ' '              ! 'nbins,X1,X2'
    character(len=3)   :: convert= ' '              ! 's2h', 'h2s', ''
    character(len=255) :: details= ' '              ! prnt only detailed params
                                                    ! E.g., 'NumCompletedChunks'
    character(len=512) :: phaseNames= ' '           ! E.g., 'core,core+r3,..'
    character(len=8)   :: tabulate = 'no'           ! tabulate
    integer            :: hdfVersion = HDFVERSION_5
    integer            :: finalPhase = 12           ! phase number ~ total
    integer            :: nHosts = 0                ! number of hosts
    real(r4)           :: longChunks = 0._r4
  end type options_T
  
  type ( options_T )   :: Options
  type(Stat_T)         :: Statistic

  real(r4), parameter  :: Undefinedvalue = -999.99
  logical, parameter   :: Countempty = .true.
  logical, parameter   :: Showdateandtime = .false.
  integer, parameter   :: Maxfiles = 2000
  integer, parameter   :: Maxphases = 50
  integer, parameter   :: Maxchunks = 3600
  character(len=255)   :: filename          ! input filename
  character(len=4096)  :: longChunkList = ''
  character(len=4096)  :: tempChunkList = ''
  integer              :: n_filenames
  integer              :: how_many, i, status, error ! Counting indices & Error flags
  integer              :: fileID
  integer              :: fromgrpID
  integer              :: fileAccess
  logical              :: is_hdf5
  ! logical              :: verbose = .false.
  integer              :: numPhases
  real                 :: t1
  real                 :: t2
  real                 :: tFile
  logical              :: showTimings
  real(r4), dimension(:,:), pointer   :: alltimings => NULL()
  integer(kind=hSize_t), dimension(3) :: dims, maxDims
  character(len=255), dimension(MAXFILES) :: filenames
  real(r4), dimension(:,:,:), pointer :: l2auxValue => NULL()
  real(r4), dimension(:), pointer     :: timings => NULL()
  integer, dimension(MAXCHUNKS)       :: which
  ! 
  call set_config ( useToolkit = .false., logFileUnit = -1 )
  outputoptions%nArrayElmntsPerLine = 40
  time_config%use_wall_clock = .true.
  call mls_h5open(error)
  if ( error /= 0 ) then
    print *, 'Sorry--unable to start hdf5'
  endif
  n_filenames = 0
  do      ! Loop over filenames
     call get_filename(filename, n_filenames, options)
     if ( filename(1:1) == '-' ) cycle
     if ( filename == ' ' ) exit
     if ( mls_exists(trim(filename)) /= 0 ) then
       print *, 'Sorry--file not found: ', trim(filename)
       cycle
     endif
     n_filenames = n_filenames + 1
     filenames(n_filenames) = filename
  enddo
  showTimings = (options%details == ' ')
  statsOnOneLine = options%oneLine
  if ( options%convert == 's2h' ) &
     & outputOptions%sdFormatDefault = '(F11.2)'
  if ( showTimings .and. SHOWDATEANDTIME ) call output_date_and_time(msg='starting chunktimes', &
    & dateFormat='yyyydoy', timeFormat='HH:mm:ss')
  if ( NumStringElements(trim(options%binopts), .true. ) > 0 ) then
    read( options%binopts, *) statistic%nbins, statistic%bounds
    statistic%nbins = max(statistic%nbins, 2)
    allocate(statistic%bincount(statistic%nbins), stat=status)
  endif
  if ( n_filenames == 0 ) then
    if ( options%verbose ) print *, 'Sorry no input files to read'
  else
    ! Check that the hdfversions of the input files accord with default
    status = 0
    do i=1, n_filenames
     call h5fis_hdf5_f(filenames(i), is_hdf5, error)
     select case (options%hdfVersion)
     case (HDFVERSION_4)
       if ( is_hdf5 ) then
         print *, 'Sorry--not recognized as hdf4 file: ', trim(filenames(i))
         status = 1
         cycle
       endif
     case (HDFVERSION_5)
       if ( .not. is_hdf5 ) then
         print *, 'Sorry--not recognized as hdf5 file: ', trim(filenames(i))
         status = 1
         cycle
       endif
     case default
       print *, 'Sorry--unrecognized hdfVersion: ', options%hdfVersion
       status = 1
       cycle
     end select
    enddo
    if ( options%showQManager ) then
      allocate(alltimings(MAXCHUNKS, n_filenames), stat=status)
      alltimings = UNDEFINEDVALUE
    else
      allocate(alltimings(MAXCHUNKS, MAXPHASES), stat=status)
      alltimings = UNDEFINEDVALUE
    endif
    call time_now ( t1 )
    ! if ( options%verbose ) print *, 'Reading chunk l2aux file(s)'
    do i=1, n_filenames
      call time_now ( tFile )
      if ( options%verbose ) then
        if ( options%merge ) then
          print *, 'Merging from: ', trim(filenames(i))
        elseif ( options%oneLine ) then
          call newLine
          call output( trim(filenames(i)) // ': ', advance='no' )
        else
          print *, 'Reading from: ', trim(filenames(i))
        endif
      endif
      fileAccess = DFACC_RDONLY
      fileID = mls_sfstart ( trim(filenames(i)), fileAccess, &
           &                               hdfVersion=options%hdfVersion )
      call GetHDF5DSDims(fileID, trim(options%DSname), dims, maxDims)
      call h5gopen_f(fileID, '/', fromgrpid, status)
      if ( IsHDF5AttributePresent(fromgrpID, 'Phase Names') ) then
        call GetHDF5Attribute (fromgrpID, 'Phase Names', options%phaseNames)
        ! if ( options%verbose ) &
        !   & call output(' Phase Names: ' // trim(options%phaseNames), advance='yes')
      else
        if ( options%verbose ) call output(' attribute not found', advance='yes')
      endif
      call h5gclose_f(fromgrpid, status)
      if ( showTimings ) then
        ! print *, 'dims ', dims
        ! print *, 'maxdims ', maxdims
        ! stop
        allocate(l2auxValue(dims(1), dims(2), dims(3)), stat=status)
        allocate(timings(dims(3)), stat=status)
        l2auxValue = UNDEFINEDVALUE
        call LoadFromHDF5DS (fileID, trim(options%DSname), l2auxValue)
        ! New special feature:
        ! Try to guess how many phases there were by scanning
        ! data for lasst phase with values > 0
        if ( options%guessFinalPhase ) then
          ! print *, 'dims: ', dims
          options%finalPhase = FindLast(l2auxValue(1,:,:) > 0.d0, options='-n' )
          ! print *, 'finalPhase: ', options%finalPhase
        endif
        select case (lowercase(options%convert))
        case ('s2h')
          where (l2auxValue > 0.0)
            l2auxValue = l2auxValue / 3600
          elsewhere
            l2auxValue = UNDEFINEDVALUE
          end where
        case ('h2s')
          where (l2auxValue > 0.0)
            l2auxValue = l2auxValue * 3600
          elsewhere
            l2auxValue = UNDEFINEDVALUE
          end where
        case default
          where (l2auxValue <= 0.0)
            l2auxValue = UNDEFINEDVALUE
          end where
        end select
        timings = l2auxValue(1, options%finalPhase, :)
        if ( options%tabulate /= 'no' ) then
          ! fileID = mls_sfstart ( trim(filenames(1)), fileAccess, &
          !   &                               hdfVersion=options%hdfVersion )
          call tabulate(l2auxValue( 1, 1:options%finalPhase, :), &
            & options%phaseNames, options%tabulate )
        endif
        numPhases = min( MaxPhases, size(l2auxvalue,2) )
        if ( options%showQManager ) then
          alltimings(1:dims(3), i) = timings
        else
          alltimings( 1:size(l2auxvalue,3), 1:numPhases ) = &
            & transpose( l2auxValue(1, 1:numPhases, :) )
        endif
        if ( .not. options%merge ) statistic%count=0
        call statistics(real(timings, r8), statistic, real(UNDEFINEDVALUE, r8))
        if ( options%longChunks > 0._r4 ) then
          call FindAll((timings > options%longChunks), which, how_many=how_many)
          if ( how_many > 0 ) then
            tempChunkList = catLists(longChunkList, which(1:how_many))
            call GetUniqueList( tempChunkList, longChunkList, how_many, &
              & options='-e' )
          endif
        endif
        if ( showTimings .and. .not. options%merge ) then
          if ( options%showStats ) call dumpstat(statistic)
          if ( options%longChunks > 0._r4 ) &
            & call dump(longChunkList, 'list of long chunks')
          longChunkList = ' '
        endif
        deallocate(l2auxValue, stat=status)
      endif
      if ( options%showFailed ) then
        call dumpFailedChunks(fileID)
      endif
      if ( showTimings ) then
        deallocate(timings, stat=status)
      endif
      status = mls_sfend( fileID,hdfVersion=options%hdfVersion )
      if ( options%verbose .and. .not. options%oneLine ) &
        &  call sayTime('reading this file', tFile)
    enddo
    if ( options%merge ) then
      if ( options%showStats) call dumpstat(statistic)
      if ( options%longChunks > 0._r4 ) &
        & call dump(longChunkList, 'list of long chunks')
    endif
    if ( options%verbose ) call sayTime('reading all files')
    if ( options%showQManager ) then
      call QManager(alltimings, options%nHosts, options%verbose)
    endif
  endif
  deallocate(alltimings, stat=status)
  if ( showTimings .and. SHOWDATEANDTIME ) call output_date_and_time(msg='ending chunktimes', &
    & dateFormat='yyyydoy', timeFormat='HH:mm:ss')
  call mls_h5close(error)
contains

!------------------------- get_filename ---------------------
    subroutine get_filename(filename, n_filenames, options)
    ! Added for command-line processing
     character(LEN=255), intent(out) :: filename          ! filename
     integer, intent(in)             :: n_filenames
     type ( options_T ), intent(inout) :: options
     ! character(LEN=*), intent(inout) :: outputFile        ! output filename
     ! logical, intent(inout)          :: verbose
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
      elseif ( filename(1:3) == '-d ' ) then
        call getarg ( i+1+hp, options%DSName )
        i = i + 1
      elseif ( filename(1:2) == '-b' ) then
        call getarg ( i+1+hp, options%binopts )
        i = i + 1
        ! exit
      elseif ( filename(1:6) == '-detai' ) then
        call getarg ( i+1+hp, options%details )
        i = i + 1
        ! exit
      elseif ( filename(1:5) == '-hdf ' ) then
        call getarg ( i+1+hp, filename )
        read(filename, *) options%hdfVersion
        i = i + 1
      elseif ( filename(1:3) == '-l ' ) then
        call getarg ( i+1+hp, filename )
        read(filename, *) options%longChunks
        i = i + 1
      elseif ( filename(1:3) == '-n ' ) then
        call getarg ( i+1+hp, filename )
        read(filename, *) options%finalPhase
        options%guessFinalPhase = .false.
        i = i + 1
      elseif ( filename(1:6) == '-nstat' ) then
        options%showStats = .false.
        exit
      elseif ( filename(1:5) == '-s2h ' ) then
        options%convert = 's2h'
        exit
      elseif ( filename(1:5) == '-h2s ' ) then
        options%convert = 'h2s'
        exit
      elseif ( filename(1:2) == '-m' ) then
        options%merge = .true.
        exit
      elseif ( filename(1:4) == '-one' ) then
        options%oneLine = .true.
        exit
      elseif ( filename(1:5) == '-fail' ) then
        options%showFailed = .true.
        exit
      elseif ( filename(1:3) == '-t' ) then
        options%tabulate = 'yes'
        exit
      elseif ( filename(1:3) == '-tf' ) then
        options%tabulate = 'full'
        exit
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
        exit
      elseif ( filename(1:2) == '-q ' ) then
        call getarg ( i+1+hp, filename )
        options%showQManager = .true.
        read(filename, *) options%nHosts
        i = i + 1
      elseif ( filename(1:5) == '-wher' ) then
        options%showFailed = .true.
        options%showWhereFailed = .true.
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
      print *,  "Enter the name of the hdf5 l2aux-dgm file. " // &
       &  "The default output file name will be used."
      read(*,'(a)') filename
    endif
    
  end subroutine get_filename

!------------------------- deduceFailedChunks ---------------------
  subroutine deduceFailedChunks(number, failedChunks)
    ! Deduce chunks that failed
    ! Based on timings == UNDEFINEDVALUE
    ! Args
    integer, intent(out) :: number
    character(len=*), intent(out) :: failedChunks
    ! Internal variables
    ! Executable
    number = 0
    failedChunks = ' '
    if ( .not. any(timings == UNDEFINEDVALUE) ) return
    call FindAll((timings == UNDEFINEDVALUE), which, how_many=number)
    if ( number > 0 ) then
      tempChunkList = catLists(failedChunks, which(1:number))
      failedChunks = tempChunkList
    endif
  end subroutine deduceFailedChunks

!------------------------- deduceWhereChunksFailed ---------------------
  subroutine deduceWhereChunksFailed( fileID, numberFailed, &
    & whereFailed )
    ! Deduce where chunks failed
    ! Based on dnwt-chiSqRatio-${PhaseName}
    ! Args
    integer, intent(in) :: fileID
    integer, intent(in) :: numberFailed
    integer, dimension(:), intent(inout) :: whereFailed
    ! Internal variables
    integer :: chunk
    integer(kind=hSize_t), dimension(3) :: dims, maxDims
    integer :: i
    integer :: numPhases
    integer :: phase
    integer :: phaseNum
    integer :: total
    integer :: which
    real(r4), dimension(:,:,:), pointer   :: l2auxValue => NULL()
    integer, parameter ::          MAXDS = 1024 ! 500
    integer, parameter ::          MAXSDNAMESBUFSIZE = MAXDS*namelen
    character(len=namelen) :: DSName
    character (len=MAXSDNAMESBUFSIZE) :: mySdList
    character (len=MAXSDNAMESBUFSIZE) :: subList
    real(r4), dimension(MAXCHUNKS, MAXPHASES) :: dcsr
    ! Executable
    call GetAllHDF5DSNames (fileID, '/', mysdList)
    numPhases = NumStringElements( trim(options%phaseNames),  countEmpty )
    subList = ' '
    do i =1,  NumStringElements( trim(mySdList),  countEmpty )
      DSName = StringElement( mySdList, i, countEmpty )
      if ( index( DSName, 'dnwt-chiSqRatio-' ) > 0 ) then
        subList = catLists( subList, DSName )
      endif
    enddo
    ! print *, trim(sublist)
    DSName = StringElement( subList, phase, countEmpty )
    call GetHDF5DSDims(fileID, trim(DSname), dims, maxDims)
    ! print *, 'dims: ', dims
    allocate(l2auxValue(dims(1), dims(2), dims(3)), stat=status)
    dcsr= 0.
    do phase=1, NumStringElements( trim(subList),  countEmpty )
      DSName = StringElement( subList, phase, countEmpty )
      do phaseNum=1, numPhases
        if ( lowercase(DSName) == 'dnwt-chisqratio-' // &
        & trim( &
        & StringElement( lowercase(options%phaseNames), phaseNum, countEmpty ) &
        & ) &
        & ) exit
      enddo
      if ( phaseNum  > numPhases ) phaseNum = 0
      l2auxValue = 0.
      if ( phaseNum > 0 ) &
        & call LoadFromHDF5DS ( fileID, trim(DSname), l2auxValue )
      ! print *, phase, phaseNum, trim(DSName), maxval(l2auxValue), minval(l2auxValue), &
      ! & shape(l2auxValue)
      if ( phaseNum > 0 ) dcsr(:dims(3), PhaseNum) = l2auxValue(1, 1, :)
    enddo
    ! call dump( dcsr(:dims(3),:numPhases), 'dnwt-chi^2 Ratio' )
    total = 0
    do chunk=1, MAXCHUNKS
      which = FindFirst( dcsr(chunk,:), UNDEFINEDVALUE )
      if ( which > 0 .and. which <= options%finalPhase ) then
        total = total + 1
        ! print *, total, chunk, which, trim( StringElement( options%phaseNames, which, countEmpty ) )
        if ( total > numberFailed ) return
        whereFailed(total) = which
      endif
    enddo
  end subroutine deduceWhereChunksFailed

!------------------------- dumpFailedChunks ---------------------
  subroutine dumpFailedChunks(fileID)
  ! Print info on chunks that failed
  ! Args
  integer, intent(in)         :: fileID
  ! Internal variables
  integer                     :: failure
  integer                     :: grp_id
  integer                     :: number
  integer                     :: numberFailed
  integer                     :: returnStatus
  character(len=4096)         :: failedChunks
  character(len=4096)         :: failedMachines
  character(len=32) :: myPhase
  character(len=*), parameter :: NUMCOMPLETEATTRIBUTENAME = 'NumCompletedChunks'
  character(len=*), parameter :: NUMFAILATTRIBUTENAME = 'NumFailedChunks'
  character(len=*), parameter :: FAILATTRIBUTENAME = 'FailedChunks'
  character(len=*), parameter :: MACHATTRIBUTENAME = 'FailedMachines'
  logical :: showAll
  logical :: showThis
  integer, dimension(MAXCHUNKS) :: whereFailed
  ! Executable
  showAll = (options%details == ' ')
  number = -1
  numberFailed = -1
  call h5gopen_f(fileId, '/', grp_id, returnStatus)
  showThis = showAll .or. &
    & StringElementNum(options%details, NUMCOMPLETEATTRIBUTENAME, COUNTEMPTY) > 0
  if ( showThis ) call output(NUMCOMPLETEATTRIBUTENAME, advance='no')
  if ( IsHDF5AttributePresent(grp_id, NUMCOMPLETEATTRIBUTENAME) ) then
    call GetHDF5Attribute (grp_ID, NUMCOMPLETEATTRIBUTENAME, number)
    if ( showThis ) call blanks(4)
    if ( showThis ) call output(number, advance='yes')
  else
    if ( showThis ) call output(' attribute not found', advance='yes')
  endif

  showThis = showAll .or. &
    & StringElementNum(options%details, NUMFAILATTRIBUTENAME, COUNTEMPTY) > 0
  if ( showThis ) call output(NUMFAILATTRIBUTENAME, advance='no')
  if ( IsHDF5AttributePresent(grp_id, NUMFAILATTRIBUTENAME) ) then
    call GetHDF5Attribute (grp_ID, NUMFAILATTRIBUTENAME, number)
    if ( showThis ) call blanks(4)
    if ( showThis ) call output(number, advance='yes')
    numberFailed = number
  else
    if ( showThis ) call output(' attribute not found', advance='yes')
  endif

  showThis = showAll .or. &
    & StringElementNum(options%details, FAILATTRIBUTENAME, COUNTEMPTY) > 0
  if ( IsHDF5AttributePresent(grp_id, FAILATTRIBUTENAME) ) then
    call GetHDF5Attribute (grp_ID, FAILATTRIBUTENAME, failedChunks)
    if ( showThis ) call dump(trim(failedChunks), FAILATTRIBUTENAME)
  else
    if ( showThis ) call output(FAILATTRIBUTENAME, advance='no')
    if ( showThis ) call output(' attribute not found', advance='yes')
  endif

  showThis = showAll .or. &
    & StringElementNum(options%details, MACHATTRIBUTENAME, COUNTEMPTY) > 0
  if ( IsHDF5AttributePresent(grp_id, MACHATTRIBUTENAME) ) then
    call GetHDF5Attribute (grp_ID, MACHATTRIBUTENAME, failedMachines)
    if ( showThis ) call dump(trim(failedMachines), MACHATTRIBUTENAME)
  else
    if ( showThis ) call output(MACHATTRIBUTENAME, advance='no')
    if ( showThis ) call output(' attribute not found', advance='yes')
  endif

  call h5gclose_f(grp_id, returnStatus)
  if ( number < 0 ) then
    call deduceFailedChunks(number, FailedChunks)
    call output(NUMFAILATTRIBUTENAME, advance='no')
    call blanks(4)
    call output(number, advance='yes')
    call dump(trim(failedChunks), FAILATTRIBUTENAME)
    numberFailed = number
  endif
  
  if ( options%showWhereFailed .and. numberFailed > 0 ) then
    whereFailed = -1
    call deduceWhereChunksFailed( fileID, numberFailed, whereFailed )
    call output( 'chunk    phase         phase name', advance='yes' )
    do failure=1, numberFailed
      call GetStringElement( failedChunks, myPhase, failure, .FALSE. )
      call output( trim(myPhase), advance='no' )
      call blanks(6)
      call output( whereFailed(failure), advance='no' )
      call GetStringElement( options%phaseNames, myPhase, whereFailed(failure), .FALSE. )
      call blanks(8)
      call output( trim(myPhase), advance='yes' )
    enddo
  endif
  end subroutine dumpFailedChunks

!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'Usage:chunktimes [options] [filenames]'
      write (*,*) &
      & '  where each filename should be an hdf5 l2aux-dgm file'
      write (*,*) &
      & '(Defaults shown in ())'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options: '
      write (*,*) '-d DSname   => read DSName from list of filenames'
      write (*,*) '                      ("phase timing")'
      write (*,*) '-b "binning options" ("")'
      write (*,*) '               in form "nbins,X1,X2", where'
      write (*,*) '               nbins  => number of bins'
      write (*,*) '               X1,X2  => lower,upper bounds'
      write (*,*) '               the first bin will contain chunks < X1'
      write (*,*) '               the last bin will contain chunks > X2'
      write (*,*) '-details "details" ("")'
      write (*,*) '               in form "param1,param2,..", where'
      write (*,*) '               print only the value[s] of param1[..]'
      write (*,*) '-hdf m      => hdfVersion is m (5)'
      write (*,*) '-l t        => show chunks that took longer than t'
      write (*,*) '-n n        => use phase number n (12)'
      write (*,*) '-nstat      => skip showing statistics'
      write (*,*) '-s2h        => convert from sec to hours; or'
      write (*,*) '-h2s        => convert from hours to sec'
      write (*,*) '-q n        => show queue managers hypothetical'
      write (*,*) '               performance with n hosts (dont)'
      write (*,*) '-v          => switch on verbose mode (off)'
      write (*,*) '-m[erge]    => merge data from all files (dont)'
      write (*,*) '-one        => print statistics on one line (dont)'
      write (*,*) '-fail       => show failed chunks (dont)'
      write (*,*) '-t[abulate] => tabulate data as chunk vs. total time (dont)'
      write (*,*) '-tf         => print full tables showing time for each phase'
      write (*,*) '-where      => show where chunks failed (dont)'
      write (*,*) '-h          => print brief help'
      stop
  end subroutine print_help

!------------------------- QManager ---------------------
  subroutine QManager ( table, nHosts, verbose )
    ! Args
    real(r4), dimension(:,:), intent(in) :: table
    integer, intent(in) :: nHosts
    logical, intent(in) :: verbose
    ! Internal variables
    integer :: chunk
    integer :: nextchunk
    integer :: host
    integer :: nexthost
    integer :: master
    integer :: thismaster
    integer :: nMasters
    integer, dimension(1) :: iarray
    integer, dimension(nHosts) :: chunkNumber
    integer, dimension(nHosts) :: masterNumber
    logical, dimension(nHosts) :: HostBusy
    real(r4), dimension(nHosts) :: HostTimes
    real(r4), dimension(nHosts) :: HostTimesProjected
    real(r4), dimension(size(table, 2)) :: MasterTimes
    real(r4), dimension(size(table, 2)) :: MaxTimes
    logical, dimension(size(table, 1), size(table, 2)) :: ChunkAssigned
    logical, dimension(size(table, 1), size(table, 2)) :: ChunkDone
    real(r4) :: unQTime
    logical, parameter :: DEEBUG = .false.
    ! Executable
    nMasters = size(table, 2)
    if ( verbose ) then
      call output('Starting QManager with ')
      call output(nMasters)
      call output(' Masters ', advance='yes')
    endif
    ChunkDone = (table < 0._r4)
    ChunkAssigned = (table < 0._r4)
    HostBusy = .false.
    HostTimes = 0._r4
    MasterTimes = UNDEFINEDVALUE
    unQTime = 0._r4
    do master=1, nMasters
      MaxTimes(master) = maxval(table(:, master))
      unQTime = unQTime + MaxTimes(master)
    enddo
    ! Loop over assigning chunks to hosts from current master
    ! chunk = 1
    master = 1
    chunk = findFirst(table(:, master) > 0._r4)
    do while(any(.not. ChunkDone))
      do while(any(.not. HostBusy) .and. any(.not. ChunkAssigned))
        ! Find free host
        host = findFirst(.not. HostBusy)
        ! HostTimes(host) = HostTimes(host) + table(chunk, master)
        ChunkAssigned(chunk, master) = .true.
        HostBusy(host) = .true.
        chunkNumber(host) = chunk
        masterNumber(host) = master
        if ( DEEBUG ) print *, 'Assigned ', chunk, master, host
        nextchunk = FindNext(table(:, master) > 0._r4, chunk)
        chunk = nextchunk
        if ( chunk > 0 ) cycle
        if ( master == nMasters ) exit
        master = master + 1
        chunk = findFirst(table(:, master) > 0._r4)
        if ( verbose ) then
          call output('Starting new master ')
          call output(master, advance='yes')
        endif
      enddo
      ! Now find next host, chunk to finish
      HostTimesProjected = UNDEFINEDVALUE
      do host = 1, nHosts
        if ( .not. HostBusy(host) ) cycle
        HostTimesProjected(host) = HostTimes(host) + &
          & table(chunkNumber(host), masterNumber(host))
      enddo
      iarray = minloc(HostTimesProjected, HostTimesProjected > 0._r4)
      nextHost = iarray(1)
      host = nexthost
      HostTimes(host) = HostTimesProjected(host)
      nextchunk = chunkNumber(host)
      thismaster = masterNumber(host)
      HostBusy(host) = .false.
      ChunkDone(nextchunk, thismaster) = .true.
      if ( DEEBUG ) print *, 'finished ', nextchunk, thismaster, host
      ! Now check if this was last chunk for this master
      if ( all(ChunkDone(:, thismaster)) ) then
        MasterTimes(thisMaster) = HostTimes(host)
        if ( verbose ) then
          call output('All chunks done for master ')
          call output(thismaster, advance='yes')
          call output(MasterTimes(thisMaster), advance='yes')
        endif
      endif
    enddo
    ! Print summary of hosttimes
    call output('Summary of QManager performance', advance='yes')
    call output('Time with QManager :', advance='no')
    call output(maxval(hosttimes), advance='yes')
    call output('Time without QManager :', advance='no')
    call output(unQTime, advance='yes')
    !
    call output('Breakdown by master', advance='yes')
    call output('master', advance='no')
    call blanks(6)
    call output('Q time', advance='no')
    call blanks(6)
    call output('noQ time', advance='yes')
    unQTime = 0._r4
    do master=1, nmasters
      unQTime = unQTime + maxTimes(master)
      call output(master)
      call blanks(6)
      call output(masterTimes(master), advance='no')
      call blanks(2)
      call output(unQTime, advance='yes')
    enddo
    !
    call output('Breakdown by host', advance='yes')
    call output('host', advance='no')
    call blanks(6)
    call output('time', advance='yes')
    do host=1, nHosts
      iarray = minloc(hosttimes, hosttimes > 0._r4)
      nexthost = iarray(1)
      call output(nexthost)
      call blanks(6)
      call output(hostTimes(nexthost), advance='yes')
      hostTimes(nexthost) = UNDEFINEDVALUE
    enddo
  end subroutine QManager

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

!------------------------- tabulate ---------------------
  subroutine tabulate ( table, phases, tabulateHow )
    ! Args
    real(r4), dimension(:,:), intent(in) :: table
    character(len=*), intent(in) :: phases
    character(len=*), intent(in) :: tabulateHow
    ! Internal variables
    integer :: chunk
    character(len=32) :: myPhase
    integer :: numPhases
    integer :: phase
    ! Executable
    if ( tabulateHow == 'full' ) then
      numPhases = size(table, 1)
    else
      numPhases = 1
    endif
    call blanks(30, fillchar='*', advance='yes')
    call output('chunk', advance='no')
    call blanks(3)
    do phase=1, numPhases
      if ( tabulateHow /= 'full' ) then
        call output('total', advance='no')
      elseif ( phases /= ' ' ) then
        call GetStringElement( phases, myPhase, phase, .FALSE.)
        call output( myPhase(:14), advance='no' )
      else
        call output('phase', advance='no')
        call output(phase, advance='no')
        call blanks(3)
      endif
    enddo
    call newline
    do chunk=1, size(table, 2)
      if ( any(table(:, chunk) /= UNDEFINEDVALUE) ) then
        call output(chunk, advance='no')
        call blanks(1)
        if ( tabulateHow /= 'full' ) then
          call output( table(size(table, 1), chunk), advance='no' )
        else
          call output(table(:, chunk), advance='no')
          !write(*, '(24F11.2)', advance='no') table(:, chunk)
        endif
        call newline
      endif
    enddo
  end subroutine tabulate

!==================
end program chunktimes
!==================

! $Log$
! Revision 1.31  2018/10/19 01:16:55  pwagner
! CamelCase use statements
!
! Revision 1.30  2016/10/04 22:13:34  pwagner
! Builds properly with some Dumps moved to Dump_1
!
! Revision 1.29  2014/03/07 21:42:59  pwagner
! Name_Len changed to nameLen; got from MLSCommon
!
! Revision 1.28  2014/01/09 00:31:26  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 1.27  2013/08/23 02:51:47  vsnyder
! Move PrintItOut to PrintIt_m
!
! Revision 1.26  2013/08/12 23:50:59  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 1.25  2013/06/05 18:47:39  pwagner
! table should print properly, not wrap
!
! Revision 1.24  2012/07/10 23:14:44  pwagner
! Changed api in accord with GetUniqueList
!
! Revision 1.23  2009/07/08 00:39:34  pwagner
! Should not crash when num chunks is 3500 (as with 1 profile/chunk)
!
! Revision 1.22  2008/08/19 00:39:25  pwagner
! Fixed bug overstepping array bounds
!
! Revision 1.21  2008/06/17 00:06:31  pwagner
! Can force printing stats on single line
!
! Revision 1.20  2007/11/28 18:58:03  pwagner
! Increased limits; should last 5y
!
! Revision 1.19  2007/08/17 00:45:58  pwagner
! Can successfully deduce where chunks failed
!
! Revision 1.18  2007/07/19 00:28:54  pwagner
! Restored ability to print full table of chunks/phases
!
! Revision 1.17  2007/07/18 00:13:35  pwagner
! -where option to show where chunks crashed
!
! Revision 1.16  2007/06/21 22:10:59  pwagner
! -t[abulate] tabulates only total time, not all phases
!
! Revision 1.15  2007/06/14 21:47:01  pwagner
! Should not guessFinalPhase if told so explicitly with -n option
!
! Revision 1.14  2007/03/26 22:55:25  pwagner
! Prints actual Phase Names as column headers when tabulating
!
! Revision 1.13  2006/08/12 00:09:43  pwagner
! Automatically scans timings to guess how many phases
!
! Revision 1.12  2005/09/23 21:01:13  pwagner
! use_wall_clock now a component of time_config
!
! Revision 1.11  2005/07/20 20:38:00  pwagner
! Made defaults consistent with v1.51 (final phase is 12th)
!
! Revision 1.10  2005/06/22 19:27:32  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 1.9  2005/04/18 16:27:18  pwagner
! Mistakenly deallocated timings before possibly needing to use it--fixed
!
! Revision 1.8  2005/04/15 20:08:19  pwagner
! Clarified and corrected type of files required
!
! Revision 1.7  2005/03/24 21:18:52  pwagner
! Avoid printing starting, ending times when useless
!
! Revision 1.6  2005/03/18 01:01:19  pwagner
! -details option narrows output from -fail
!
! Revision 1.5  2005/03/04 18:50:15  pwagner
! Can simulate l2q performance
!
! Revision 1.4  2004/09/28 23:13:18  pwagner
! Added -l, -nstat options; may deduce failed chunks
!
! Revision 1.3  2004/09/16 23:58:46  pwagner
! Reports machine names of failed chunks
!
! Revision 1.2  2004/09/16 00:20:56  pwagner
! May show failed chunks
!
! Revision 1.1  2004/09/14 18:52:37  pwagner
! First commit
!
