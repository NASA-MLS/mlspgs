! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

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
   use MLSHDF5, only: GetHDF5Attribute, GetHDF5DSDims, LoadFromHDF5DS, &
    &  IsHDF5AttributePresent, mls_h5open, mls_h5close
   use MLSMessageModule, only: MLSMessageConfig
   use MLSStats1, only: STAT_T, &
     & ALLSTATS, DUMPSTAT=>DUMP, MLSMIN, MLSMAX, MLSMEAN, MLSSTDDEV, MLSRMS, STATISTICS
   use MLSStringLists, only: GetStringElement, NumStringElements
   use MLSStrings, only: lowercase
   use output_m, only: blanks, newline, output
   use PCFHdr, only: GlobalAttributes
   use Time_M, only: Time_Now, USE_WALL_CLOCK
   
   implicit none

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! Reads chunk times from list of input files

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
    character(len=255) :: DSName= 'phase timing'    ! Dataset name
    character(len=255) :: binopts= ' '              ! 'nbins,X1,X2'
    character(len=3)   :: convert= 's2h'            ! 's2h', 'h2s', ''
    integer            :: hdfVersion = HDFVERSION_5
    integer            :: finalPhase = 10           ! phase number ~ total
  end type options_T
  
  type ( options_T ) :: options
  type(STAT_T)       :: statistic

  integer, parameter ::          MAXFILES = 100
  integer, parameter ::          MAXPHASES = 100
  integer, parameter ::          MAXCHUNKS = 360
  character(len=255) :: filename          ! input filename
  character(len=255), dimension(MAXFILES) :: filenames
  integer            :: n_filenames
  integer     ::  i, count, status, error ! Counting indices & Error flags
  integer(kind=hSize_t), dimension(3) :: dims, maxDims
  integer     ::  fileID
  real(r4), dimension(:,:,:), pointer   :: l2auxValue => NULL()
  real(r4), dimension(:), pointer     :: timings => NULL()
  logical     :: is_hdf5
  ! logical     :: verbose = .false.
  real        :: t1
  real        :: t2
  real        :: tFile
  real(r4), parameter :: UNDEFINEDVALUE = -999.99
  ! 
  MLSMessageConfig%useToolkit = .false.
  MLSMessageConfig%logFileUnit = -1
  USE_WALL_CLOCK = .true.
  CALL mls_h5open(error)
  if ( status /= 0 ) then
    print *, 'Sorry--unable to allocate timings'
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
      fileID = mls_sfstart ( trim(filenames(i)), DFACC_RDONLY, &
           &                               hdfVersion=options%hdfVersion )
      call GetHDF5DSDims(fileID, trim(options%DSname), dims, maxDims)
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
      if ( .not. options%merge ) statistic%count=0
      call statistics(real(timings, r8), statistic, real(UNDEFINEDVALUE, r8))
      if ( .not. options%merge ) call dumpstat(statistic)
      if ( options%showFailed ) then
        call dumpFailedChunks(fileID)
      endif
      status = mls_sfend( fileID,hdfVersion=options%hdfVersion )
      deallocate(l2auxValue, timings, stat=status)
      if ( options%verbose )  call sayTime('reading this file', tFile)
    enddo
    if ( options%merge ) call dumpstat(statistic)
    if ( options%verbose ) call sayTime('reading all files')
  endif
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
        exit
      elseif ( filename(1:2) == '-b' ) then
        call getarg ( i+1+hp, options%binopts )
        i = i + 1
        exit
      elseif ( filename(1:5) == '-hdf ' ) then
        call getarg ( i+1+hp, filename )
        read(filename, *) options%hdfVersion
        i = i + 1
        exit
      elseif ( filename(1:3) == '-n ' ) then
        call getarg ( i+1+hp, filename )
        read(filename, *) options%finalPhase
        i = i + 1
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
      print *,  "Enter the name of the HDFEOS4 or 5 L2GP file. " // &
       &  "The default output file name will be used."
      read(*,'(a)') filename
    endif
    
  end subroutine get_filename

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
  ! Executable
  call h5gopen_f(fileId, '/', grp_id, returnStatus)
  call output(NUMCOMPLETEATTRIBUTENAME, advance='no')
  if ( IsHDF5AttributePresent(grp_id, NUMCOMPLETEATTRIBUTENAME) ) then
    call GetHDF5Attribute (grp_ID, NUMCOMPLETEATTRIBUTENAME, number)
    call blanks(4)
    call output(number, advance='yes')
  else
    call output(' attribute not found', advance='yes')
  endif

  call output(NUMFAILATTRIBUTENAME, advance='no')
  if ( IsHDF5AttributePresent(grp_id, NUMFAILATTRIBUTENAME) ) then
    call GetHDF5Attribute (grp_ID, NUMFAILATTRIBUTENAME, number)
    call blanks(4)
    call output(number, advance='yes')
  else
    call output(' attribute not found', advance='yes')
  endif

  if ( IsHDF5AttributePresent(grp_id, FAILATTRIBUTENAME) ) then
    call GetHDF5Attribute (grp_ID, FAILATTRIBUTENAME, failedChunks)
    call dump(trim(failedChunks), FAILATTRIBUTENAME)
  else
    call output(FAILATTRIBUTENAME, advance='no')
    call output(' attribute not found', advance='yes')
  endif

  call h5gclose_f(grp_id, returnStatus)
  end subroutine dumpFailedChunks

!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'Usage:chunktimes [options] [filenames]'
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
      write (*,*) '          -hdf m      => hdfVersion is m (5)'
      write (*,*) '          -n n        => use phase number n (10)'
      write (*,*) '          -s2h        => convert from sec to hours; or'
      write (*,*) '          -h2s        => convert from hours to sec (-s2h)'
      write (*,*) '          -v          => switch on verbose mode (off)'
      write (*,*) '          -m[erge]    => merge data from all files (dont)'
      write (*,*) '          -fail       => show failed chunks (dont)'
      write (*,*) '          -t[abulate] => print data in tables (dont)'
      write (*,*) '          -h          => print brief help'
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
! Revision 1.1  2004/09/14 18:52:37  pwagner
! First commit
!
