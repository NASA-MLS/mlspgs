! Copyright 2008, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

program wrapLines
   use Io_Stuff, only: Get_Lun
   use Machine, only: Hp, Getarg
   use MLSMacros, only: Dump_Macros, Read_Macros, Expand_Line
   use MLSStringlists, only: Wrap
   use MLSStrings, only: Lentrimtoascii, Nappearances
   use Output_M, only: Output
   implicit none

!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------
  ! wrap lines in stdin (or -i inputFile)
  ! so they do not exceed (by much) 128 chars in length
  
  ! Suggested use:
  ! After m4, pipe l2cf through me before saving as *.l2cf
  ! Then you'll be able to read the resulting l2cf more comfortably
  
  ! Notes and limitations
  ! Input line must not exceed MAXLINELEN (24000) chars
  ! (a bug in NAG makes it even less (1000), so .. don't build with NAG)

  ! Lines will be split at a BREAK (','), unless that line contains quotes
  ! In case you're wondering why we don't wrap lines with quoted strings,
  ! it's because we're not smart enough to tell when a BREAK is embedded
  ! within quotes
  ! which could be trouble if you have directory/file names with BREAKs
  ! in them.
  
  ! E.g., a line with 
  !  ForwardModelGlobal, l2pc='/data/emls/l2cal/l2pc039,not,a,good,idea,to,use,commas,here/MLS-Aura_L2Cal-L2PC-band2-LATSCALARHIRES_v1-7-0-l2pc039_m01.h5'
  ! would be erroneously split thus:
  !  ForwardModelGlobal, l2pc='/data/emls/l2cal/l2pc039,not,a,good,idea,to,use,commas, $
  !  here/MLS-Aura_L2Cal-L2PC-band2-LATSCALARHIRES_v1-7-0-l2pc039_m01.h5'
  
  ! Emerging use:
  ! We hope to replace m4 eventually
  ! Given a file of macros, use their definitions
  ! to replace such occurrences as
  ! (1) "!macro" 
  !      is replaced with value of macro
  ! (2) "!function(args)" 
  !      is replaced with function evaluated for its args
  ! (3) "!forall(function,args)"
  !      is replaced one time for each arg with function evaluated ath that arg
  ! The macros file may contain statements such as
  !   macro_name=value
  !   function(arg)="character string inluding $(arg)"
  !   function()="character string inluding $(n)"

  type options_T
    logical             :: verbose            = .false.
    logical             :: summarize          = .false.
    logical             :: keepTrailingSpaces = .true. ! Less efficient if true
    ! logical             :: wrapAllLines       = .false. ! Even with quoted strs
    logical             :: doWrap             = .true. ! wrap?
    integer             :: printCBL           = 0 ! How many consecutive bls
    integer             :: linefeed           = 10
    integer             :: width              = 128
    character(len=1)    :: BREAK              = ','
    character(len=1)    :: comment            = ';'
    character(len=1)    :: MODE               = 's'
    character(len=3)    :: eol                = ' $'
    character(len=3)    :: quotes             = '''"'
    character(len=255)  :: inputFile          = ' '     ! input filename       
    character(len=255)  :: macrosFile         = 'none'  ! output filename       
  end type options_T
  
  type ( options_T ) :: options

  ! variables
  integer, parameter :: MAXLINELEN = 24000
  integer :: addedLines
  character(len=3) :: advance
  character(len=1) :: c
  character(len=1), dimension(MAXLINELEN) :: cArray
  logical :: containsQuotes
  integer :: i
  integer :: iounit
  character(len=MAXLINELEN) :: lineIn
  character(len=MAXLINELEN) :: lineOut
  integer :: LineLengthRead ! Max encountered
  integer :: nCommentsRead  ! number read
  integer :: nConsectiveBlanks  ! number read
  integer :: nLinesRead     ! number read
  integer :: nQuotesRead    ! number read
  integer :: pos
  character(len=1), dimension(3) :: qArray ! array of quote chars
  integer :: status
  integer :: totalBlankLinesDropped
  integer :: totalLinesAdded
  integer :: totalLinesNeedWrapping
  character(len=12) :: xfmt
  character(len=8) :: xlen
  ! Executable
  ! write (*,*) 'starting'
  call get_options (options)
  ! What format do we use for reading each line?
  xfmt = '(128a1)' ! This is the default
  write( xlen, '(i8)' ) len(lineIn)
  if ( index(xlen, '*') < 1 ) xfmt = '(' // trim(adjustl(xlen)) // 'a1)'
  advance = 'yes'
  if ( .not. options%keepTrailingSpaces ) advance = 'no'
  c = options%comment
  ! write (*,*) 'options%macrosfile ', options%macrosfile
  if ( len_trim(options%macrosfile) > 0 .and. options%macrosfile /= 'none' ) then
    call read_macros ( options%macrosFile, status )
    call dump_macros
  endif
  ! Do we wrap a file or stdin?
  if ( len_trim(options%inputFile) < 1 ) then
    iounit = 5
  else
    call get_lun ( iounit )
    open( iounit, status='OLD', form='FORMATTED', &
      & recl=MAXLINELEN, file=options%inputFile )
  endif
  if ( options%verbose ) call dumpSettings( options )
  LineLengthRead         = 0  
  nCommentsRead          = 0  
  nConsectiveBlanks      = 0
  nLinesRead             = 0  
  nQuotesRead            = 0  
  totalBlankLinesDropped = 0
  totalLinesAdded        = 0  
  totalLinesNeedWrapping = 0
  ! write (*,*) 'options%keepTrailingSpaces ', options%keepTrailingSpaces
  do
    if ( .not. options%keepTrailingSpaces ) then
      read( iounit, '(a)', iostat=status ) lineIn
      ! write (*,*) trim(lineIn)
      if ( options%macrosFile /= 'none' ) then
        call expand_line( LineIn, LineOut )
        LineIn = LineOut
      endif
      LineLengthRead = max( LineLengthRead, len_trim(lineIn) )
    else
      i = 0
      lineIn = ' '
      status = 0
      call null_fill_1d( carray )
      read( iounit, fmt=xfmt, eor=50, end=500, err=50, advance='no' ) cArray
500   status = -1
50    if ( status /= 0 ) exit
      ! write (*,*) 'Read a line starting with: ', cArray(1:1), ichar(cArray(1:1))
      ! Possibly expand given line using macros
      if ( options%macrosFile /= 'none' ) then
        LineIn = transfer( cArray, LineIn )
        pos = index( LineIn, achar(0) )
        ! write(*,*) 'pos: ', pos
        call expand_line( LineIn(:pos-1), LineOut )
        cArray = transfer( LineOut, cArray )
        lineIn = ' '
      endif
      ! write (*,*) 'Expanded, it starts with: ', cArray(1:1), ichar(cArray(1:1))
      oneLine: do pos=1, len(lineIn) - 1
        if ( any(carray(pos:pos+1) == achar(0)) ) exit oneLine
        i = min(i + 1, len(lineIn))
        lineIn(i:i) = carray(pos)
      enddo oneLine
      LineLengthRead = max( LineLengthRead, pos-1 )
      i = min(i + 1, len(lineIn))
      lineIn(i:i) = achar(options%linefeed)
      ! write (*,*) trim(lineIn)
    endif
    ! write (*,*) 'status: ', status
    if ( status /= 0 ) exit
    nLinesRead = nLinesRead + 1
    addedLines = 0
    containsQuotes = .false.
    if ( len_trim(options%quotes) > 0 ) then
      qArray = transfer( options%quotes, qArray )
      containsQuotes = any( &
        & NAppearances( trim(lineIn), qArray(1:len_trim(options%quotes)) ) > 0 &
        & )
    endif
    if ( containsQuotes ) nQuotesRead = nQuotesRead + 1
    if ( lenTrimToAscii(lineIn) < 1 ) then
      nConsectiveBlanks = nConsectiveBlanks + 1
    else
      nConsectiveBlanks = 0
    endif
    if ( 0 < options%printCBL .and. options%printCBL < nConsectiveBlanks ) then
      totalBlankLinesDropped = totalBlankLinesDropped + 1
      cycle
    endif
    if ( index(trim(lineIn), c) > 0 ) then
      lineOut = lineIn
      nCommentsRead = nCommentsRead + 1
    elseif ( .not. options%doWrap .or. containsQuotes ) then
      lineOut = lineIn
    elseif ( len_trim(options%quotes) < 1 ) then
      call wrap( lineIn, lineOut, options%width, &
        & inseparator=trim(options%eol) // achar(options%linefeed), &
        & break=options%break, mode=options%mode, &
        & addedLines=addedLines )
    else
      call wrap( lineIn, lineOut, options%width, &
        & inseparator=trim(options%eol) // achar(options%linefeed), &
        & break=options%break, mode=options%mode, &
        & quotes=trim(options%quotes), addedLines=addedLines )
    endif
    ! write (*,*) trim(lineOut)
    ! Don't double-space unless you mean it
    pos = len_trim( lineOut )
    if ( lineOut(pos:pos) == achar(options%linefeed) ) then
      call output( trim(lineOut), advance='no' )
    else
      call output( trim(lineOut), advance='yes' )
    endif
    totalLinesAdded = totalLinesAdded + addedLines
    if ( addedLines > 0 ) totalLinesNeedWrapping = totalLinesNeedWrapping + 1
  enddo
  if ( options%verbose .or. options%summarize ) then
     print *, c // ' ( I n p u t   f i l e   s u m m a r y )'
     print *, c // 'lines read               ', nLinesRead
     print *, c // 'longest line read        ', LineLengthRead
     print *, c // 'lines with comments      ', nCommentsRead
     print *, c // 'lines with quotes        ', nQuotesRead
     print *, c // 'lines needed wrapping    ', totalLinesNeedWrapping
     print *, c // 'lines added by wrapping  ', totalLinesAdded
     print *, c // 'blank lines dropped      ', totalBlankLinesDropped
  endif
contains
  subroutine null_fill_1d( array, nullChar )
    ! Fill array with null chars
    ! Args
    character(len=*), dimension(:), intent(out) :: array
    character(len=1), optional, intent(in)      :: nullChar
    ! Internal variables
    character(len=1) :: myNull
    integer :: col
    integer :: pos
    ! Executable
    myNull = achar(0)
    if ( present(nullChar) ) myNull = nullChar
    do col=1, size(array)
      do pos=1, len(array(1))
        array(col)(pos:pos) = myNull
      enddo
    enddo
  end subroutine null_fill_1d

!------------------------- dumpSettings ---------------------
    subroutine dumpSettings( options )
    ! Added for command-line processing
     type ( options_T ), intent(in)   :: options
     ! Local variables
     print *, c // ' --- wrapLines settings ---'
     print *, c // 'verbose?                 ', options%verbose
     print *, c // 'summarize?               ', options%summarize
     print *, c // 'wrap?                    ', options%doWrap
     print *, c // 'keep trailing spaces  ?  ', options%keepTrailingSpaces
     print *, c // 'width                    ', options%width
     print *, c // 'print consec blanks      ', options%printCBL
     print *, c // 'break                    ', options%break
     print *, c // 'comment                  ', options%comment
     print *, c // 'mode                     ', options%mode 
     print *, c // 'eol                      ', options%eol  
     print *, c // 'linefeed                 ', options%linefeed
     print *, c // 'quotes                   ', options%quotes
     print *, c // 'xfmt                     ', xfmt
     print *, c // 'advance                  ', advance
     if ( len_trim(options%inputFile) < 1 ) then
     print *, c // 'input  file              ', '<STDIN>'
     else
     print *, c // 'input  file              ', trim(options%inputFile)
     endif
     print *, c // 'macros file              ', trim(options%macrosFile)
     print *, c // ' --- End wrapLines settings ---'
    end subroutine dumpSettings

!------------------------- get_options  ---------------------
    subroutine get_options ( options )
    ! Added for command-line processing
     type ( options_T ), intent(inout) :: options
     ! Local variables
     integer ::                         error = 1
     character(len=255)       :: filename
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
      elseif ( filename(1:3) == '-i ' ) then
        call getarg ( i+1+hp, options%inputFile )
        i = i + 1
      elseif ( filename(1:3) == '-Df' ) then
        call getarg ( i+1+hp, options%macrosFile )
        i = i + 1
      else if ( filename(1:6) == '-blank' ) then
        call getarg ( i+1+hp, number )
        read(number, *) options%printCBL
        i = i + 1
      else if ( filename(1:2) == '-b' ) then
        call getarg ( i+1+hp, options%break )
        options%break = adjustl(options%break)
        i = i + 1
      else if ( filename(1:2) == '-c' ) then
        call getarg ( i+1+hp, options%comment )
        options%comment = adjustl(options%comment)
        i = i + 1
      else if ( filename(1:5) == '-eol ' ) then
        call getarg ( i+1+hp, options%eol )
        options%eol = adjustl(options%eol)
        i = i + 1
      else if ( filename(1:2) == '-q' ) then
        call getarg ( i+1+hp, options%quotes )
        options%quotes = adjustl(options%quotes)
        i = i + 1
      else if ( filename(1:6) == '-mode ' ) then
        call getarg ( i+1+hp, options%mode )
        options%mode = adjustl(options%mode)
        i = i + 1
      else if ( filename(1:7) == '-width ' ) then
        call getarg ( i+1+hp, number )
        read(number, *) options%width
        i = i + 1
      elseif ( filename(1:5) == '-keep' ) then
        options%keepTrailingSpaces = .true.
      elseif ( filename(1:6) == '-nkeep' ) then
        options%keepTrailingSpaces = .false.
      elseif ( filename(1:5) == '-wrap' ) then
        options%doWrap = .true.
      elseif ( filename(1:6) == '-nwrap' ) then
        options%doWrap = .false.
      elseif ( filename(1:2) == '-s' ) then
        options%summarize = .true.
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
      else
        call print_help
      end if
      i = i + 1
    end do
    if ( error /= 0 ) then
      call print_help
    endif
    i = i + 1
    
  end subroutine get_options 
!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'Usage:wrapLines [options]'
      write (*,*) ' Options: -v            => switch on verbose mode'  
      write (*,*) '          -s            => print num of lines, longest'     
      write (*,*) '                            at end of stdout'
      write (*,*) '          -b char       => break lines at char'     
      write (*,*) '                          (default is ",")'         
      write (*,*) '          -blank n      => print no more than n consecutive'
      write (*,*) '                          blank lines (default is infinite)'         
      write (*,*) '          -c char       => treat char as comment'     
      write (*,*) '                          (default is ";")'         
      write (*,*) '          -i file       => wrap file'     
      write (*,*) '                          (default is to wrap stdin)'         
      write (*,*) '          -Df file      => get macros definitions from file'     
      write (*,*) '                          (default is to do no macros expansion)'         
      write (*,*) '          -eol chars    => glue chars to end of broken lines'
      write (*,*) '                            (default is ", $")'
      write (*,*) '          -[n]keep      => do [not] keep trailing spaces'     
      write (*,*) '                            (default is to keep them)'
      write (*,*) '          -[n]wrap      => do [not] wrap long lines'     
      write (*,*) '                            (default is to wrap)'
      write (*,*) '          -q chars      => let chars delimit quoted strings'     
      write (*,*) '                           (a line with quotes will not wrap)'     
      write (*,*) '                          (defaults are ''")'         
      write (*,*) '          -mode mode    => set wrap mode'
      write (*,*) '                            (default is "s" for soft)'
      write (*,*) '          -width width  => wrap lines longer than width'
      write (*,*) '                            (default is 128)'
      write (*,*) '          -h            => print brief help'
      stop
  end subroutine print_help
end program wrapLines
! $Log$
! Revision 1.7  2018/08/13 23:14:16  pwagner
! Use statements made Camel Case
!
! Revision 1.6  2016/12/16 21:58:29  pwagner
! Works with new wrap
!
! Revision 1.5  2012/08/14 21:15:43  pwagner
! l2cf could have uncommented comment lines; fixed
!
! Revision 1.4  2012/08/07 18:04:46  pwagner
! Can wrap a supplied filename instead of stdin
!
! Revision 1.3  2012/08/03 16:54:04  pwagner
! Can now expand some macros using definitions from separate macros file
!
! Revision 1.2  2008/06/17 00:05:15  pwagner
! Can skip consecutive blanks
!
! Revision 1.1  2008/05/22 17:39:37  pwagner
! First commit
!
