! Copyright (c) 2005, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contracts NAS7-1407/NAS7-03001 is acknowledged.

!=================================
program chunktimes ! Reads chunk times from l2aux file(s)
!=================================

   use dump_0, only: dump
   use Hdf, only: DFACC_CREATE, DFACC_RDWR, DFACC_RDONLY, DFACC_READ
   use HDF5, only: HSIZE_T, h5fis_hdf5_f, h5gclose_f, h5gopen_f
   use MACHINE, only: FILSEP, HP, IO_ERROR, GETARG
   use MLSCommon, only: R4, R8
   use MLSFiles, only: mls_exists, MLS_SFSTART, MLS_SFEND, &
     & HDFVERSION_4, HDFVERSION_5
   use MLSHDF5, only: CpHDF5Attribute, GetHDF5Attribute, GetHDF5DSDims, &
     & IsHDF5AttributePresent, LoadFromHDF5DS, mls_h5open, mls_h5close
   use MLSMessageModule, only: MLSMessageConfig
   use MLSSets, only: FindAll, FindFirst, FindNext
   use MLSStats1, only: STAT_T, &
     & ALLSTATS, DUMPSTAT=>DUMP, MLSMIN, MLSMAX, MLSMEAN, MLSSTDDEV, MLSRMS, STATISTICS
   use MLSStringLists, only: catLists, GetStringElement, GetUniqueList, &
     & NumStringElements, StringElementNum
   use MLSStrings, only: lowercase
   use output_m, only: blanks, newline, output, output_date_and_time
   use PCFHdr, only: GlobalAttributes
   use Time_M, only: Time_Now, USE_WALL_CLOCK
   
   implicit none

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! Reads chunk times and failures from list of input files

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it
! LF95.Linux/test [options] [input files]
  type options_T
    logical            :: verbose = .false.
    logical            :: merge = .false.           ! Merge data from files
    logical            :: tabulate = .false.        ! tabulate
    logical            :: showFailed = .false.      ! show howmany, which failed
    logical            :: showStats = .true.        ! show max, min, mean, etc.
    logical            :: showQManager = .false.    ! show QManager performance
    character(len=255) :: DSName= 'phase timing'    ! Dataset name
    character(len=255) :: binopts= ' '              ! 'nbins,X1,X2'
    character(len=255) :: details= ' '              ! prnt only detailed params
                                                    ! E.g., 'NumCompletedChunks'
    character(len=3)   :: convert= ' '              ! 's2h', 'h2s', ''
    integer            :: hdfVersion = HDFVERSION_5
    integer            :: finalPhase = 10           ! phase number ~ total
    integer            :: nHosts = 0                ! number of hosts
    real(r4)           :: longChunks = 0._r4
  end type options_T
  
  type ( options_T ) :: options
  type(STAT_T)       :: statistic

  logical, parameter ::          COUNTEMPTY = .true.
  logical, parameter ::          MODIFYPHASENAMES = .false.
  logical, parameter ::          SHOWDATEANDTIME = .false.
  integer, parameter ::          MAXFILES = 100
  integer, parameter ::          MAXPHASES = 100
  integer, parameter ::          MAXCHUNKS = 360
  character(len=*), parameter :: NEWPHASENAMES = &
    & 'initptan,updateptan,inituth,core,coreplusr2,highcloud,coreplusr3,' // &
    & 'coreplusr4,coreplusr5'
  character(len=255) :: filename          ! input filename
  character(len=4096):: longChunkList = ''
  character(len=4096):: tempChunkList = ''
  character(len=255), dimension(MAXFILES) :: filenames
  integer            :: n_filenames
  integer     ::  how_many, i, count, status, error ! Counting indices & Error flags
  integer(kind=hSize_t), dimension(3) :: dims, maxDims
  integer     ::  fileID, oldfileID
  integer     ::  togrpID, fromgrpID
  integer     ::  fileAccess
  real(r4), dimension(:,:,:), pointer   :: l2auxValue => NULL()
  real(r4), dimension(:), pointer     :: timings => NULL()
  real(r4), dimension(:,:), pointer     :: alltimings => NULL()
  logical     :: is_hdf5
  ! logical     :: verbose = .false.
  real        :: t1
  real        :: t2
  real        :: tFile
  real(r4), parameter :: UNDEFINEDVALUE = -999.99
  integer, dimension(MAXCHUNKS) :: which
  logical :: showTimings
  ! 
  MLSMessageConfig%useToolkit = .false.
  MLSMessageConfig%logFileUnit = -1
  USE_WALL_CLOCK = .true.
  CALL mls_h5open(error)
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
  if ( showTimings .and. SHOWDATEANDTIME ) call output_date_and_time(msg='starting chunktimes', &
    & dateFormat='yyyydoy', timeFormat='HH:mm:ss')
  if ( NumStringElements(trim(options%binopts), .true. ) > 0 ) then
    read(trim(options%binopts), *) statistic%nbins, statistic%bounds
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
    endif
    call time_now ( t1 )
    if ( options%verbose ) print *, 'Reading chunk times'
    do i=1, n_filenames
      call time_now ( tFile )
      if ( options%verbose ) then
        if ( options%merge) then
          print *, 'Merging from: ', trim(filenames(i))
        else
          print *, 'Reading from: ', trim(filenames(i))
        endif
      endif
      if ( MODIFYPHASENAMES .and. i == n_filenames ) then
        fileAccess = DFACC_RDWR
      else
        fileAccess = DFACC_RDONLY
      endif
      fileID = mls_sfstart ( trim(filenames(i)), fileAccess, &
           &                               hdfVersion=options%hdfVersion )
      call GetHDF5DSDims(fileID, trim(options%DSname), dims, maxDims)
      if ( showTimings ) then
        ! print *, 'dims ', dims
        ! print *, 'maxdims ', maxdims
        ! stop
        allocate(l2auxValue(dims(1), dims(2), dims(3)), stat=status)
        allocate(timings(dims(3)), stat=status)
        l2auxValue = UNDEFINEDVALUE
        call LoadFromHDF5DS (fileID, trim(options%DSname), l2auxValue)
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
        if ( options%tabulate ) &
          & call tabulate(l2auxValue(1, 1:options%finalPhase, :))
        if ( options%showQManager ) &
          & alltimings(1:dims(3), i) = timings
        if ( .not. options%merge ) statistic%count=0
        call statistics(real(timings, r8), statistic, real(UNDEFINEDVALUE, r8))
        if ( options%longChunks > 0._r4 ) then
          call FindAll((timings > options%longChunks), which, how_many=how_many)
          if ( how_many > 0 ) then
            tempChunkList = catLists(longChunkList, which(1:how_many))
            call GetUniqueList(tempChunkList, longChunkList, how_many, COUNTEMPTY)
          endif
        endif
        if ( .not. options%merge ) then
          if ( options%showStats) call dumpstat(statistic)
          if ( options%longChunks > 0._r4 ) &
            & call dump(longChunkList, 'list of long chunks')
          longChunkList = ' '
        endif
        deallocate(l2auxValue, timings, stat=status)
      endif
      if ( options%showFailed ) then
        call dumpFailedChunks(fileID)
      endif
      ! The following is just to test whether we can copy phase names
      ! from one file to another
      ! Remove it after testing, please
      if ( MODIFYPHASENAMES  .and. i == n_filenames ) then
        oldfileID = mls_sfstart ( trim(filenames(1)), fileAccess, &
           &                               hdfVersion=options%hdfVersion )
        call h5gopen_f(fileID, '/', togrpid, status)
        call h5gopen_f(oldfileID, '/', fromgrpid, status)
        call CpHDF5Attribute(fromgrpID, togrpID, 'Phase Names')
        call h5gclose_f(togrpid, status)
        call h5gclose_f(fromgrpid, status)
        status = mls_sfend( oldfileID,hdfVersion=options%hdfVersion )
      endif
      status = mls_sfend( fileID,hdfVersion=options%hdfVersion )
      if ( options%verbose )  call sayTime('reading this file', tFile)
    enddo
    if ( options%merge ) then
      if ( options%showStats) call dumpstat(statistic)
      if ( options%longChunks > 0._r4 ) &
        & call dump(longChunkList, 'list of long chunks')
    endif
    if ( options%verbose ) call sayTime('reading all files')
    if ( options%showQManager ) then
      call QManager(alltimings, options%nHosts, options%verbose)
      deallocate(alltimings, stat=status)
    endif
  endif
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
      elseif ( filename(1:5) == '-fail' ) then
        options%showFailed = .true.
        exit
      elseif ( filename(1:2) == '-t' ) then
        options%tabulate = .true.
        exit
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
        exit
      elseif ( filename(1:2) == '-q ' ) then
        call getarg ( i+1+hp, filename )
        options%showQManager = .true.
        read(filename, *) options%nHosts
        i = i + 1
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

!------------------------- dumpFailedChunks ---------------------
  subroutine dumpFailedChunks(fileID)
  ! Print info on chunks that failed
  ! Args
  integer, intent(in)         :: fileID
  ! Internal variables
  integer                     :: grp_id
  integer                     :: number
  integer                     :: returnStatus
  character(len=4096)         :: failedChunks
  character(len=*), parameter :: NUMCOMPLETEATTRIBUTENAME = 'NumCompletedChunks'
  character(len=*), parameter :: NUMFAILATTRIBUTENAME = 'NumFailedChunks'
  character(len=*), parameter :: FAILATTRIBUTENAME = 'FailedChunks'
  character(len=*), parameter :: MACHATTRIBUTENAME = 'FailedMachines'
  logical :: showAll
  logical :: showThis
  ! Executable
  showAll = (options%details == ' ')
  number = -1
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
    call GetHDF5Attribute (grp_ID, MACHATTRIBUTENAME, failedChunks)
    if ( showThis ) call dump(trim(failedChunks), MACHATTRIBUTENAME)
  else
    if ( showThis ) call output(MACHATTRIBUTENAME, advance='no')
    if ( showThis ) call output(' attribute not found', advance='yes')
  endif

  call h5gclose_f(grp_id, returnStatus)
  if ( number > -1 ) return
  call deduceFailedChunks(number, FailedChunks)
  call output(NUMFAILATTRIBUTENAME, advance='no')
  call blanks(4)
  call output(number, advance='yes')
  call dump(trim(failedChunks), FAILATTRIBUTENAME)
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
      write (*,*) ' Options: -d DSname   => read DSName from list of filenames'
      write (*,*) '                         ("phase timing")'
      write (*,*) '          -b "binning options" ("")'
      write (*,*) '                         in form "nbins,X1,X2", where'
      write (*,*) '                         nbins  => number of bins'
      write (*,*) '                         X1,X2  => lower,upper bounds'
      write (*,*) '                         the first bin will contain chunks < X1'
      write (*,*) '                         the last bin will contain chunks > X2'
      write (*,*) '          -details "details" ("")'
      write (*,*) '                         in form "param1,param2,..", where'
      write (*,*) '                         print only the value[s] of param1[..]'
      write (*,*) '          -hdf m      => hdfVersion is m (5)'
      write (*,*) '          -l t        => show chunks that took longer than t'
      write (*,*) '          -n n        => use phase number n (10)'
      write (*,*) '          -nstat      => skip showing statistics'
      write (*,*) '          -s2h        => convert from sec to hours; or'
      write (*,*) '          -h2s        => convert from hours to sec'
      write (*,*) '          -q n        => show queue managers hypothetical'
      write (*,*) '                         performance with n hosts (dont)'
      write (*,*) '          -v          => switch on verbose mode (off)'
      write (*,*) '          -m[erge]    => merge data from all files (dont)'
      write (*,*) '          -fail       => show failed chunks (dont)'
      write (*,*) '          -t[abulate] => print data in tables (dont)'
      write (*,*) '          -h          => print brief help'
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
  subroutine tabulate ( table, phases )
    ! Args
    real(r4), dimension(:,:), intent(in) :: table
    character(len=*), optional, intent(in) :: phases
    ! Internal variables
    integer :: chunk
    character(len=32) :: myPhase
    integer :: numPhases
    integer :: phase
    ! Executable
    numPhases = size(table, 1)
    call blanks(30, fillchar='*', advance='yes')
    call output('chunk', advance='no')
    call blanks(3)
    do phase=1, numPhases
      if ( present(phases) ) then
        call GetStringElement( phases, myPhase, phase, .FALSE.)
        call output(trim(myPhase), advance='no')
        call blanks(3)
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
        call output(table(:, chunk), advance='no')
        call newline
      endif
    enddo
  end subroutine tabulate

!==================
end program chunktimes
!==================

! $Log$
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
