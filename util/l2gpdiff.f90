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

   use Dump_Options, only: DumpDumpOptions, DumpTableSide, &
     & PrintNameAtLineEnd, StatsOnOneLine, RmsFormat
   use HighOutput, only: OutputNamedValue
   use Io_Stuff, only: Get_Lun
   use L2GPData, only: Diff, MaxSwathNamesBufSize
   use Machine, only: Hp, Getarg, NeverCrash
   use MLSFiles, only: MLS_Exists, HDFVersion_5, MLS_InqSwath
   use MLSHDF5, only: MLS_H5Open, MLS_H5Close
   use MLSMessageModule, only: MLSMessageConfig, MLSMessage, MLSMSG_Error
   use MLSStringLists, only: CatLists, ExpandStringRange
   use MLSStrings, only: LowerCase, WriteIntsToChars
   use Output_M, only: ResumeOutput, SuspendOutput, Output
   use PrintIt_M, only: Set_Config
   use Time_M, only: Time_Now, Time_Config
   
   implicit none

!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

  integer, parameter ::          MAXFIELDSLENGTH = 256
  type Options_T
    character(len=MAXFIELDSLENGTH) :: chunks = '*' ! wild card means 'all'
    character(len=MAXFIELDSLENGTH) :: pressures = '*' ! wild card means 'all'
    character(len=255) ::  geoBoxNames = '' ! which geolocation names to box
    character(len=MAXFIELDSLENGTH) :: fields = '*' ! wild card means 'all'
    character(len=8) :: HEADSIDE = 'left' ! on which side stats headers printed
    character(len=MAXSWATHNAMESBUFSIZE) :: swaths1 = '*'
    character(len=MAXSWATHNAMESBUFSIZE) :: swaths2 = '*'
    character(len=80) :: dumpOptions       = ' '
    integer            ::  nGeoBoxDims     = 0
    real, dimension(4) ::  geoBoxLowBound
    real, dimension(4) ::  geoBoxHiBound
    integer     ::         Details         = 1
    logical     ::         force           = .false.
    logical     ::         AuBrick         = .false.
    logical     ::         rms             = .false.
    logical     ::         stats           = .false.
    logical     ::         table           = .false.
    logical     ::         silent          = .false.
    logical     ::         timing          = .false.
    logical     ::         verbose         = .false.
    logical     ::         debug           = .false.
    logical     ::         matchTimes      = .false.
    logical     ::         showMissing     = .false.
    logical     ::         ignoreBadChunks = .false.
    integer     ::         numDiffs = 0
  end type Options_T
  
  type ( Options_T ) :: options

  integer, parameter                      :: MAXFILES = 100
  integer, parameter                      :: MAXNCHUNKS = 50

  integer, dimension(MAXNCHUNKS)          :: chunks
  ! character(len=8)                        :: dumpOptions
  character(len=255)                      :: filename          ! input filename
  character(len=255), dimension(MAXFILES) :: filenames
  integer                                 :: i, error ! Counting indices & Error flags
  integer                                 :: listSize
  integer                                 :: nChunks
  integer                                 :: n_filenames
  integer                                 :: nPressures
  integer                                 :: numDiffs
  integer                                 :: NUMSWATHSPERFILE
  real, dimension(MAXNCHUNKS)             :: pressures
  character(len=16)                       :: string
  character(len=MAXSWATHNAMESBUFSIZE)     :: swathList1
  real                                    :: t1
  real                                    :: t2
  real                                    :: tFile
  logical, parameter                      :: USEALLINPUTSWATHS = .true.
  ! 
  call set_config ( useToolkit = .false., logFileUnit = -1 )
  time_config%use_wall_clock = .true.
  CALL mls_h5open(error)
  statsOnOneLine = .false.
  ! PrintNameAtLineEnd = .true.
  n_filenames = 0
  do      ! Loop over filenames, options
     call get_filename(filename, options)
     if ( filename(1:1) == '-' ) cycle
     if ( filename == ' ' ) exit
     if ( mls_exists(trim(filename)) /= 0 ) then
       print *, 'Sorry--file not found: ', trim(filename)
       cycle
     endif
     n_filenames = n_filenames + 1
     filenames(n_filenames) = filename
  enddo
  dumpTableSide = options%headSide ! 'left'
  if ( options%debug ) call dumpProgramOptions
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
  rmsFormat = '(1pe8.1)'
  if ( options%silent ) call suspendOutput
  ! dumpOptions = '-'
  if ( options%ignoreBadChunks ) options%dumpOptions = trim(options%dumpOptions) // 'i'
  if ( options%silent ) options%dumpOptions = trim(options%dumpOptions) // 'm'
  if ( options%rms ) options%dumpOptions = trim(options%dumpOptions) // 'r'
  if ( options%stats ) options%dumpOptions = trim(options%dumpOptions) // 's'
  if ( options%table ) options%dumpOptions = trim(options%dumpOptions) // 'b'
  if ( options%AuBrick ) options%dumpOptions = trim(options%dumpOptions) // '@'
  if ( options%verbose ) options%dumpOptions = trim(options%dumpOptions) // 'v'
  call time_now ( t1 )
  ! print *, 'dumpOptions: ', trim(options%dumpOptions)
  do i = 2, n_filenames, 2
    call time_now ( tFile )
    if ( options%verbose ) then
      print *, 'diff (1): ', trim(filenames(i-1))
      print *, '     (2): ', trim(filenames(i))
    endif
    if ( options%nGeoBoxDims > 0 ) then
      call diff( trim(filenames(i-1)), trim(filenames(i)), &
      & options%geoBoxNames, options%geoBoxLowBound, options%geoBoxHiBound, &
      & details=options%Details, &
      & showMissing=options%showMissing, fields=options%fields, &
      & force=options%force, swaths1=options%swaths1, swaths2=options%swaths2, &
      & matchTimes=options%matchTimes, &
      & numDiffs=numDiffs, &
      & options=options%dumpOptions )
    elseif ( options%chunks == '*' .and. options%pressures == '*' ) then
      if ( options%debug ) call output( 'Standard diff', advance='yes' )
      call diff( trim(filenames(i-1)), trim(filenames(i)), &
      & details=options%Details, &
      & showMissing=options%showMissing, fields=options%fields, &
      & force=options%force, swaths1=options%swaths1, swaths2=options%swaths2, &
      & matchTimes=options%matchTimes, &
      & numDiffs=numDiffs, &
      & options=options%dumpOptions )
    elseif ( options%pressures == '*' ) then
      call ExpandStringRange(options%chunks, chunks, nchunks)
      if ( nchunks < 1 ) cycle
      call diff( trim(filenames(i-1)), trim(filenames(i)), &
      & chunks=chunks(1:nchunks), &
      & details=options%Details, &
      & showMissing=options%showMissing, fields=options%fields, &
      & force=options%force, swaths1=options%swaths1, swaths2=options%swaths2, &
      & matchTimes=options%matchTimes, &
      & numDiffs=numDiffs, &
      & options=options%dumpOptions )
    elseif ( options%chunks == '*' ) then
      call ExpandStringRange(options%pressures, pressures, npressures)
      if ( npressures < 1 ) cycle
      call diff( trim(filenames(i-1)), trim(filenames(i)), &
      & pressures=pressures(1:nPressures), &
      & details=options%Details, &
      & showMissing=options%showMissing, fields=options%fields, &
      & force=options%force, swaths1=options%swaths1, swaths2=options%swaths2, &
      & matchTimes=options%matchTimes, &
      & numDiffs=numDiffs, &
      & options=options%dumpOptions )
    else
      call ExpandStringRange(options%chunks, chunks, nchunks)
      if ( nchunks < 1 ) cycle
      call ExpandStringRange(options%pressures, pressures, npressures)
      if ( npressures < 1 ) cycle
      call diff( trim(filenames(i-1)), trim(filenames(i)), &
      & pressures=pressures(1:nPressures), chunks=chunks(1:nchunks), &
      & details=options%Details, &
      & showMissing=options%showMissing, fields=options%fields, &
      & force=options%force, swaths1=options%swaths1, swaths2=options%swaths2, &
      & matchTimes=options%matchTimes, &
      & numDiffs=numDiffs, &
      & options=options%dumpOptions )
    endif
    options%numDiffs = options%numDiffs + numDiffs
    if ( options%timing ) call sayTime('diff these files: ' // trim(filenames(i-1)) // &
      & ' and ' // trim(filenames(i)), tFile)
  enddo
  if ( options%timing ) call sayTime('diffing all files')
  call resumeOutput
  if ( options%silent .and. options%numDiffs > 0 ) then
    call WriteIntsToChars ( options%numDiffs, string )
    call print_string(string)
  endif
  call mls_h5close(error)
contains
!------------------------- get_filenames ---------------------
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
      else if ( filename(1:6) == '-chunk' ) then
        call getarg ( i+1+hp, options%chunks )
        i = i + 1
        exit
      elseif ( filename(1:6) == '-crash' ) then
        MLSMessageConfig%crashOnAnyError = .true.
        neverCrash = .false.
        exit
      elseif ( filename(1:6) == '-force' ) then
        options%force = .true.
        exit
      else if ( filename(1:4) == '-geo' ) then
        call getarg ( i+1+hp, Chars )
        options%geoBoxNames = catLists( options%geoBoxNames, Chars )
        i = i + 1
        options%nGeoBoxDims = min( options%nGeoBoxDims + 1, 4 )
        call getarg ( i+1+hp, Chars )
        read( Chars, * ) options%geoBoxLowBound(options%nGeoBoxDims), &
          & options%geoBoxHiBound(options%nGeoBoxDims)
        i = i + 1
        exit
      else if ( filename(1:6) == '-press' ) then
        call getarg ( i+1+hp, options%pressures )
        i = i + 1
        exit
      elseif ( filename(1:8) == '-silent ' ) then
        options%silent = .true.
        exit
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
        exit
      elseif ( filename(1:4) == '-one' ) then
        statsOnOneLine = .true.
        exit
      else if ( filename(1:6) == '-field' ) then
        call getarg ( i+1+hp, options%fields )
        i = i + 1
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
      else if ( filename(1:4) == '-mat' ) then
        options%matchTimes = .true.
        exit
      else if ( lowercase(filename(1:3)) == '-au' ) then
        options%AuBrick = .true.
        exit
      else if ( filename(1:5) == '-rms ' ) then
        options%rms = .true.
        exit
      else if ( filename(1:5) == '-side' ) then
        call getarg ( i+1+hp, options%headSide )
        i = i + 1
        exit
      else if ( filename(1:3) == '-s ' ) then
        options%stats = .true.
        exit
      else if ( filename(1:3) == '-S1' ) then
        call getarg ( i+1+hp, Chars )
        call read_swath_names ( Chars, options%swaths1 )
        i = i + 1
        exit
      else if ( filename(1:3) == '-S2' ) then
        call getarg ( i+1+hp, Chars )
        call read_swath_names ( Chars, options%swaths2 )
        i = i + 1
        exit
      else if ( filename(1:3) == '-s1' ) then
        call getarg ( i+1+hp, options%swaths1 )
        i = i + 1
        exit
      else if ( filename(1:3) == '-s2' ) then
        call getarg ( i+1+hp, options%swaths2 )
        i = i + 1
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
!------------------------- dumpProgramOptions ---------------------
  subroutine dumpProgramOptions
    ! dump options
    call outputNamedValue( 'chunks', trim(options%chunks)             )
    call outputNamedValue( 'pressures', trim(options%pressures)       )
    call outputNamedValue( 'geoBoxNames', trim(options%geoBoxNames)   )
    call outputNamedValue( 'fields', trim(options%fields)             )
    call outputNamedValue( 'headSide', trim(options%headSide)         )
    call outputNamedValue( 'swaths1', trim(options%swaths1)           )
    call outputNamedValue( 'swaths2', trim(options%swaths2)           )
    call outputNamedValue( 'nGeoBoxDims    ', options%nGeoBoxDims     )
    call outputNamedValue( 'geoBoxLowBound ', options%geoBoxLowBound  )
    call outputNamedValue( 'geoBoxHiBound  ', options%geoBoxHiBound   )
    call outputNamedValue( 'AuBrick        ', options%AuBrick         )
    call outputNamedValue( 'Details        ', options%Details         )
    call outputNamedValue( 'dumpOptions    ', options%dumpOptions     )
    call outputNamedValue( 'force          ', options%force           )
    call outputNamedValue( 'rms            ', options%rms             )
    call outputNamedValue( 'stats          ', options%stats           )
    call outputNamedValue( 'table          ', options%table           )
    call outputNamedValue( 'silent         ', options%silent          )
    call outputNamedValue( 'timing         ', options%timing          )
    call outputNamedValue( 'verbose        ', options%verbose         )
    call outputNamedValue( 'debug          ', options%debug           )
    call outputNamedValue( 'matchTimes     ', options%matchTimes      )
    call outputNamedValue( 'showMissing    ', options%showMissing     )
    call outputNamedValue( 'ignoreBadChunks', options%ignoreBadChunks )
  end subroutine dumpProgramOptions
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
    write (*,*) 'Options: -f filename => add filename to list of filenames'
    write (*,*) '                 (can do the same w/o the -f)'
    write (*,*) '  -chunks "c1,c2,.."'
    write (*,*) '              =>   diff only chunks c1, c2,..'
    write (*,*) '  -fields "f1,f2,.."'
    write (*,*) '              =>   diff only fields f1, f2,..'
    write (*,*) '  -geo name lo,hi  '
    write (*,*) '              => diff only geobox low <= geo <= hi'
    write (*,*) '                 where geo is one of'
    write (*,*) '                  {latitude, longitude, time, pressure}'
    write (*,*) '                 (may be repeated)'
    write (*,*) '                 if hi < lo then dump is outside geobox'
    write (*,*) '  -pressures "p1,p2,.."'
    write (*,*) '              =>   diff only at pressures p1,p2,..'
    write (*,*) '  -D details  => level of details to show'
    write (*,*) '  -d options  => pass options to dump routines'
    write (*,*) '                 e.g., "?" to list available ones'
    write (*,*) '  -force      => force diff even if swath names differ'
    write (*,*) '  -s1 "swath1,swath2,.."'
    write (*,*) '              =>   diff only swath1,swath2,..'
    write (*,*) '  -s2 "swath1,swath2,.."'
    write (*,*) '              =>   diff only swath1,swath2,..'
    write (*,*) '                   from even-numbered files in list'
    write (*,*) '  -S1 filename'
    write (*,*) '              => read -s1 swath names from filename'
    write (*,*) '  -S2 filename'
    write (*,*) '              => read -s2 swath names from filename'
    write (*,*) '  -v          => switch on verbose mode'
    write (*,*) '  -silent     => switch on silent mode'
    write (*,*) '                    (printing only if diffs found)'
    write (*,*) '  -crash      => crash with walkback on any error'
    write (*,*) '  -debug      => dump options, etc.'
    write (*,*) '  -ignore     => ignore bad chunks'
    write (*,*) '  -matchTimes => only matching profile times'
    write (*,*) '  -au         => format like goldbrick'
    write (*,*) '  -rms        => just print mean, rms'
    write (*,*) '  -s          => just show number of differences'
    write (*,*) '  -one        => print name on each line (dont)'
    write (*,*) '  -side "s"   => print stat headers on one of'
    write (*,*) '                  {"top", "left", "right", "bottom"}'
    write (*,*) '  -t[able]    => table of % vs. amount of differences (pdf)'
    write (*,*) '  -miss       => just show which swaths are missing'
    write (*,*) '  -h          => print brief help'
    write (*,*) '(Notes)'
    write (*,*) '(1) The -D and -fields options should be mutually exclusive'
    write (*,*) '(2) -s1 lets you pick out which swaths to diff'
    write (*,*) '(3) -s2 along with -force diffs swaths despite different names'
    write (*,*) '(4) the list of chunks may include the range operator "-"'
    stop
  end subroutine print_help

!------------------------- print_string ---------------------
  subroutine print_string(string)
    character(len=*), intent(in) :: string
    write(*,'(a)') trim(string)
  end subroutine print_string

!------------------------- read_swath_names ---------------------
  subroutine read_swath_names( filename, swathnames )
    character(len=*), intent(in)    :: filename
    character(len=*), intent(out)   :: swathnames
    !
    character(len=len(swathnames))  :: line
    integer                         :: lun
    integer                         :: stat
    !
    swathnames = ' '
      ! Find a free logical unit number
      call get_lun( lun, msg=.false. )
      if ( lun < 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & "No logical unit numbers available" )
      open ( unit=lun, file=filename,&
        & status='old', form='formatted', &
        & access='sequential', iostat=stat )
      if ( stat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Unable to open swaths file " // filename )
      do
        read ( unit=lun, fmt=*, iostat=stat ) line
        if ( stat < 0 ) exit
        if ( line(1:1) /= '#' ) swathnames = trim(swathnames) // ',' // trim(line)
      enddo
      ! Crude hackery-quackery coming right up
      if ( swathnames(1:1) == ',' ) swathnames = swathnames(2:)

  end subroutine read_swath_names

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
! Revision 1.29  2017/12/14 23:17:06  pwagner
! Fixed syntax error in help screen
!
! Revision 1.28  2017/10/12 20:27:44  pwagner
! Monkeyed with appearance of help page; removed outdated build notes
!
! Revision 1.27  2017/05/13 00:04:09  pwagner
! Added -S1 and -S2 cmdline options
!
! Revision 1.26  2016/09/09 20:38:27  pwagner
! Added Au (Gold) brick option removing some hay from the stack of statistics
!
! Revision 1.25  2016/08/09 22:45:26  pwagner
! Consistent with splitting of Dunp_0
!
! Revision 1.24  2016/03/23 16:38:51  pwagner
! Added -one commandline option
!
! Revision 1.23  2015/07/01 18:10:18  pwagner
! Stray print statement caused headaches for gold brick
!
! Revision 1.22  2015/06/30 18:52:02  pwagner
! -d '?' now dumps available dump otions
!
! Revision 1.21  2014/01/09 00:31:26  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 1.20  2013/08/23 02:51:48  vsnyder
! Move PrintItOut to PrintIt_m
!
! Revision 1.19  2012/02/13 23:42:42  pwagner
! -d opts passes opts to underlying dump routines
!
! Revision 1.18  2010/11/17 01:26:44  pwagner
! debug option dumps commandline options
!
! Revision 1.17  2009/09/11 23:24:47  pwagner
! options now include -t to be consistent with diff apis
!
! Revision 1.16  2009/06/16 22:37:59  pwagner
! Changed api for dump, diff routines; now rely on options for most optional behavior
!
! Revision 1.15  2008/09/25 23:11:54  pwagner
! May confine diffs to a geolocation box
!
! Revision 1.14  2008/04/10 20:23:03  pwagner
! Less voluminous output
!
! Revision 1.13  2007/11/28 18:59:58  pwagner
! May choose where to print stat headers
!
! Revision 1.12  2007/10/12 23:37:31  pwagner
! Removed much unused stuff
!
! Revision 1.11  2007/06/27 19:30:59  pwagner
! -pressures option added
!
! Revision 1.10  2007/06/14 21:48:41  pwagner
! Overrides default rmsFormat
!
! Revision 1.9  2007/02/27 00:05:33  pwagner
! -matchTimes option diffs only profiles with matching times
!
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
