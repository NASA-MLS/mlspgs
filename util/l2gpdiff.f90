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
program l2gpdiff ! show diffs between swaths in two different files
!=================================

   use Hdf, only: DFACC_CREATE, DFACC_RDWR, DFACC_READ
   use HDF5, only: h5fopen_f, h5fclose_f, h5fis_hdf5_f   
   use HDFEOS5, only: HE5T_NATIVE_CHAR
   use L2GPData, only: L2GPData_T, &
     & Diff, &
     & L2GPNameLen, MAXSWATHNAMESBUFSIZE
   use MACHINE, only: FILSEP, HP, IO_ERROR, GETARG
   use MLSCommon, only: R8
   use MLSFiles, only: mls_exists, &
     & HDFVERSION_4, HDFVERSION_5, MLS_INQSWATH
   use MLSHDF5, only: mls_h5open, mls_h5close
   use MLSMessageModule, only: MLSMessageConfig
   use MLSStringLists, only: ExpandStringRange, &
     & GetStringElement, NumStringElements
   use MLSStrings, only: WriteIntsToChars
   use output_m, only: resumeOutput, suspendOutput, output
   use PCFHdr, only: GlobalAttributes
   use Time_M, only: Time_Now, time_config
   
   implicit none

!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it
! LF95.Linux/test [options] [input files]
  integer, parameter ::          MAXFIELDSLENGTH = 256
  type options_T
    character(len=MAXFIELDSLENGTH) :: chunks = '*' ! wild card means 'all'
    integer     :: Details = 1
    character(len=MAXFIELDSLENGTH) :: fields = '*' ! wild card means 'all'
    logical     :: force = .false.
    logical     :: ignoreBadChunks = .false.
    logical     :: matchTimes = .false.
    logical     :: rms = .false.
    logical     :: showMissing = .false.
    logical     :: stats = .false.
    character(len=MAXSWATHNAMESBUFSIZE) :: swaths1 = '*'
    character(len=MAXSWATHNAMESBUFSIZE) :: swaths2 = '*'
    logical     :: silent = .false.
    logical     :: verbose = .false.
    integer     :: numDiffs = 0
  end type options_T
  
  type ( options_T ) :: options

  integer, parameter ::          MAXFILES = 100
  integer, parameter :: MAXNCHUNKS = 50

  integer, dimension(MAXNCHUNKS) :: chunks
  character(len=10) :: column_name
  character(len=255) :: filename          ! input filename
  character(len=255), dimension(MAXFILES) :: filenames
  integer     ::  i, error ! Counting indices & Error flags
  integer     ::  hdfversion1
  integer     ::  hdfversion2
  logical     :: is_hdf5
  integer :: listSize
  integer            :: n_filenames
  integer :: nChunks
  integer :: numDiffs
  integer :: NUMSWATHSPERFILE
  character(len=16) :: string
  character(len=MAXSWATHNAMESBUFSIZE) :: swathList1
  character(len=MAXSWATHNAMESBUFSIZE) :: swathList2
  character(len=10) :: swath_name
  character(len=*), parameter :: SWATHPREFIX='R3:240.B25D:CO.S1.DACS-1 chisqBinned Core'
  real        :: t1
  real        :: t2
  real        :: tFile
  logical, parameter :: USEALLINPUTSWATHS = .true.
  ! 
  MLSMessageConfig%useToolkit = .false.
  MLSMessageConfig%logFileUnit = -1
  time_config%use_wall_clock = .true.
  CALL mls_h5open(error)
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
  if ( n_filenames < 2 ) then
    if ( options%verbose ) print *, 'Sorry -- need at least 2 input files to diff'
    stop
  elseif ( options%verbose .and. options%silent ) then
    print *, 'Sorry--either verbose or silent; cant be both'
    stop
  endif
  if ( USEALLINPUTSWATHS ) then
    NUMSWATHSPERFILE = mls_InqSwath ( trim(filenames(1)), &
      & swathList1, listSize, hdfVersion=HDFVERSION_5)
    if ( NUMSWATHSPERFILE < 1 ) then
      print *, ' NUMSWATHSPERFILE: ', NUMSWATHSPERFILE
      print *, ' file1: ', trim(filenames(1))
      print *, ' swathList1: ', trim(swathList1)
      print *, ' listSize: ', listSize
      stop
    endif
  else
    NUMSWATHSPERFILE = 2
    swathList1 = ''
  endif
  if ( options%silent ) call suspendOutput
  call time_now ( t1 )
  do i = 2, n_filenames, 2
    call time_now ( tFile )
    if ( options%verbose ) then
      print *, 'diff (1): ', trim(filenames(i-1))
      print *, '     (2): ', trim(filenames(i))
    endif
    if ( options%chunks == '*' ) then
      call diff(trim(filenames(i-1)), trim(filenames(i)), &
      & details=options%Details, stats=options%stats, &
      & rms=options%rms, ignoreBadChunks=options%ignoreBadChunks, &
      & showMissing=options%showMissing, fields=options%fields, &
      & force=options%force, swaths1=options%swaths1, swaths2=options%swaths2, &
      & matchTimes=options%matchTimes, &
      & silent=options%silent, verbose=options%verbose, numDiffs=numDiffs )
    else
      call ExpandStringRange(options%chunks, chunks, nchunks)
      if ( nchunks < 1 ) cycle
      call diff(trim(filenames(i-1)), trim(filenames(i)), &
      & chunks(1:nchunks), &
      & details=options%Details, stats=options%stats, &
      & rms=options%rms, ignoreBadChunks=options%ignoreBadChunks, &
      & showMissing=options%showMissing, fields=options%fields, &
      & force=options%force, swaths1=options%swaths1, swaths2=options%swaths2, &
      & matchTimes=options%matchTimes, &
      & silent=options%silent, verbose=options%verbose, numDiffs=numDiffs )
    endif
    options%numDiffs = options%numDiffs + numDiffs
    call sayTime('diff these files: ' // trim(filenames(i-1)) // &
      & ' and ' // trim(filenames(i)), tFile)
  enddo
  call sayTime('diffing all files')
  call resumeOutput
  if ( options%silent .and. options%numDiffs > 0 ) then
    call WriteIntsToChars ( options%numDiffs, string )
    call print_string(string)
  endif
  call mls_h5close(error)
contains
!------------------------- get_filenames ---------------------
    subroutine get_filename(filename, n_filenames, options)
    ! Added for command-line processing
     character(LEN=255), intent(out) :: filename          ! filename
     integer, intent(in)             :: n_filenames
     type ( options_T ), intent(inout) :: options
     ! Local variables
     integer ::                         error = 1
     integer, save ::                   i = 1
     character(LEN=16)              :: detailChars
  ! Get inputfile name, process command-line args
  ! (which always start with -)
    do
      call getarg ( i+hp, filename )
      ! print *, i, ' th Arg: ', trim(filename)
      error = 0
      if ( filename(1:1) /= '-' ) exit
      if ( filename(1:3) == '-h ' ) then
        call print_help
      else if ( filename(1:6) == '-chunk' ) then
        call getarg ( i+1+hp, options%chunks )
        i = i + 1
        exit
      elseif ( filename(1:6) == '-force' ) then
        options%force = .true.
        exit
      elseif ( filename(1:8) == '-silent ' ) then
        options%silent = .true.
        exit
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
        exit
      else if ( filename(1:6) == '-field' ) then
        call getarg ( i+1+hp, options%fields )
        i = i + 1
        exit
      else if ( filename(1:3) == '-f ' ) then
        call getarg ( i+1+hp, filename )
        i = i + 1
        exit
      else if ( filename(1:3) == '-d ' ) then
        call getarg ( i+1+hp, detailChars )
        ! read(detailChars, '(i)') options%Details
        read(detailChars, *) options%Details
        i = i + 1
        exit
      else if ( filename(1:4) == '-ign' ) then
        options%ignorebadchunks = .true.
        exit
      else if ( filename(1:4) == '-mat' ) then
        options%matchTimes = .true.
        exit
      else if ( filename(1:5) == '-rms ' ) then
        options%rms = .true.
        exit
      else if ( filename(1:3) == '-s ' ) then
        options%stats = .true.
        exit
      else if ( filename(1:3) == '-s1' ) then
        call getarg ( i+1+hp, options%swaths1 )
        i = i + 1
        exit
      else if ( filename(1:3) == '-s2' ) then
        call getarg ( i+1+hp, options%swaths2 )
        i = i + 1
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
      & 'Usage:l2gpdiff [options] [filenames]'
      write (*,*) &
      & ' diffs swaths between paired l2gp files [2k] and [2k - 1]'
      write (*,*) &
      & ' each file contains one or more swaths, each swath many named fields'
      write (*,*) &
      & ' optionally restrict diffs to certain fields, chunks, etc.'
      write (*,*) ' Options: -f filename => add filename to list of filenames'
      write (*,*) '                         (can do the same w/o the -f)'
      write (*,*) '          -chunks "c1,c2,.."'
      write (*,*) '                      =>   diff only chunks c1, c2,..'
      write (*,*) '          -fields "f1,f2,.."'
      write (*,*) '                      =>   diff only fields f1, f2,..'
      write (*,*) '          -d details  => level of details to show'
      write (*,*) '          -force      => force diff even if swath names differ'
      write (*,*) '          -s1 "swath1,swath2,.."'
      write (*,*) '                      =>   diff only swath1,swath2,..'
      write (*,*) '          -s2 "swath1,swath2,.."'
      write (*,*) '                      =>   diff only swath1,swath2,..'
      write (*,*) '                           from even-numbered files in list'
      write (*,*) '          -v          => switch on verbose mode'
      write (*,*) '          -silent         => switch on silent mode'
      write (*,*) '                            (printing only if diffs found)'
      write (*,*) '          -ignore     => ignore bad chunks'
      write (*,*) '          -matchTimes => only matching profile times'
      write (*,*) '          -rms        => just print mean, rms'
      write (*,*) '          -s          => just show statistics'
      write (*,*) '          -miss       => just show which swaths are missing'
      write (*,*) '          -h          => print brief help'
      write (*,*) '    (Notes)'
      write (*,*) ' (1) The -d and -fields options should be mutually exclusive'
      write (*,*) ' (2) -s1 lets you pick out which swaths to diff'
      write (*,*) ' (3) -s2 along with -force diffs swaths despite different names'
      write (*,*) ' (4) the list of chunks may include the range operator "-"'
      stop
  end subroutine print_help

!------------------------- print_string ---------------------
  subroutine print_string(string)
    character(len=*), intent(in) :: string
    write(*,'(a)') trim(string)
  end subroutine print_string

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

!==================
end program l2gpdiff
!==================

! $Log$
! Revision 1.8  2006/03/15 19:18:37  pwagner
! Passes verbose option to L2GPData/diff
!
! Revision 1.7  2006/01/14 00:59:12  pwagner
! Added -silent, -chunks options
!
! Revision 1.6  2005/09/23 21:01:13  pwagner
! use_wall_clock now a component of time_config
!
! Revision 1.5  2005/08/19 23:41:08  pwagner
! May specify fields, swaths
!
! Revision 1.4  2005/06/22 19:27:33  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 1.3  2004/08/07 00:15:55  pwagner
! All stringlist stuff was moved from mlsstrings to mlsstringlists
!
! Revision 1.2  2004/07/22 17:10:46  pwagner
! Added -ignore and -rms options
!
! Revision 1.1  2004/06/16 17:55:06  pwagner
! First commit
!
