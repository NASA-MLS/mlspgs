program F90TEX

!{ Read a Fortran 90 program.  Write a LaTeX file.  Start the output with
! the usual {\tt$\backslash$documentclass} stuff.  In the body of the
! document, wrap code with {\tt$\backslash$begin}\{\emph{verbatim}\} and
! {\tt$\backslash$end}\{\emph{verbatim}\} where \emph{verbatim} is given
! by the variable {\tt CODE\_STYLE}.  When a comment begins with !\{, end
! the {\tt$\backslash$begin}\{\emph{verbatim}\} wrapper, and start it
! again at the end of the comment.

!>> 2000-04-07 F90TEX WV Snyder Original code

! =====     Declarations     ===========================================

  use MACHINE, only: FILSEP, HP, IO_ERROR

  implicit NONE

  logical :: BOX = .true.                         ! Put a box around !{ TeX?
  character(len=8) :: DATE                        ! of program execution
  integer :: I, J                                 ! Subscript / do inductor
  character(len=132) :: IN_FILE, OUT_FILE         ! File names
  integer :: IN_UNIT = 10, OUT_UNIT = 11          ! Unit numbers
  integer :: IOSTAT                               ! I/O status
  character(len=132) :: LINE                      ! From input
  character(len=23) :: NOW                        ! Date and time, formatted
  character :: NUMBER_STEP = '1'                  ! Line number step
  character(len=*), parameter :: START_CODE(2) = &
    (/ 'begin{lstlisting}[texcl]{ }', &
       'begin{verbatim}            ' /)
  integer :: STATE = 0   ! 0 => begin, 1 => doing code, 2 => doing TeX
  character(len=*), parameter :: STOP_CODE(2) = (/ 'end{lstlisting}', &
                                                   'end{verbatim}} ' /)
  integer :: SX = 2      ! index for START_CODE and STOP_CODE, selected
                         ! by the -l or -v options
  character(len=10) :: TIME                       ! of program execution
  character(len=*), parameter :: WIDTH = '6.25in' ! of TeX page
  character(len=*), parameter :: WIDTHP = '6in'   ! of TeX parbox inside fbox

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

! =====     Executable Statements     ==================================

  call date_and_time ( date, time )
  now = date(1:4) // '/' // date(5:6) // '/' // date(7:) // ' ' // &
        time(1:2) // ':' // time(3:4) // ':' // time(5:)
! Get input and output file names

  i = 1
  do
    call getarg ( i+hp, in_file )
    if ( in_file(1:1) /= '-' ) exit
    if ( in_file(1:2) == '- ' ) then
      i = i + 1
      call getarg ( i+hp, in_file )
      exit
    end if
    if ( in_file(1:3) == '-l ' ) then
      sx = 1
    else if ( in_file(1:3) == '-v ' ) then
      sx = 2
    else if ( in_file(1:3) == '-b ' ) then
      box = .true.
    else if ( in_file(1:4) == '-nb ' ) then
      box = .false.
    else if ( in_file(1:3) == '-n ' ) then
      i = i + 1
      call getarg ( i+hp, number_step )
    else
      call getarg ( 0+hp, in_file )
      write (*,*) &
      & 'Usage: ', trim(in_file), ' [options] [ in_file [out_file ]]'
      write (*,*) &
      & ' If in_file and out_file are absent, stdin and stdout are used.'
      write (*,*) &
      & ' If in_file is present but out_file is absent, in_file is used'
      write (*,*) &
      & '  for out_file, after removing the last "." after the last "', &
      & filsep, '"'
      write (*,*) &
      & '  and everything after it, and adding ".tex" on the end.'
      write (*,*) ' Options: -b => put a box around !{ TeX stuff (default)'
      write (*,*) '          -nb => don''t put a box around !{ TeX stuff'
      write (*,*) '          -l => use lstlisting'
      write (*,*) '          -v => use verbatim (default)'
      write (*,*) '          -n x => line number step is x (default ', &
      &                              number_step, ', 0 => no line numbers)'
      write (*,*) '          -<anything else> => this output'
      write (*,*) '          -  => no more options'
      stop
    end if
    i = i + 1
  end do
  call getarg ( i+1+hp, out_file )

  if ( in_file /= ' ' ) then
    open(in_unit, file=in_file, status='old', iostat=iostat)
    if ( iostat /= 0 ) then
      call io_error ( 'Unable to open input file', iostat, in_file )
      stop
    end if
    if ( out_file == ' ' ) then
      out_file = in_file
      j = index( out_file, filsep, back=.true. )
      i = index( out_file, '.', back=.true. )
      if ( i == 0 .or. i < j ) i = len_trim(out_file) + 1
      out_file(i:) = '.tex'
    end if
    open(out_unit, file=out_file, iostat=iostat)
    if ( iostat /= 0 ) then
      call io_error ( 'Unable to open output file', iostat, out_file )
      stop
    end if
  else
    in_unit = -1         ! Standard input
    in_file = 'Standard Input'
    out_unit = -1        ! Standard output
  end if

! Output the beginning stuff

  call output ( 'documentclass[11pt,twoside]{article}', tex=.true. )
  call output ( 'usepackage[fleqn]{amsmath}', tex=.true. )
  if ( sx == 1 ) call output ( 'usepackage{listings}', tex=.true. )
  call output ( 'textwidth ' // width, tex=.true. )
  call output ( 'oddsidemargin -0.25in', tex=.true. )
  call output ( 'evensidemargin -0.25in', tex=.true. )
  call output ( 'topmargin -0.5in', tex=.true. )
  call output ( 'textheight 9.25in', tex=.true. )
  call output ( '' )
  call output ( 'begin{document}', tex=.true. )
  call output ( '' )
  call output ( 'makeatletter', tex=.true. )
  call output ( 'def\@biblabel#1{#1.}', tex=.true. )
  call output ( 'newcommand{\ps@twolines}{%', tex=.true. )
  call output ( '  \renewcommand{\@oddhead}{%' )
  call output ( '    \parbox{' // width // '}{{\bf \Large \hfill%' )
  call output ( trim(in_file) // '%', und=.true. )
  call output ( '    }\newline%' )
  call output ( '    ' // now // &
  &             '\hfill Page \thepage\ of \pageref{lastpage}}}%' )
  call output ( '  \renewcommand{\@evenhead}{%' )
  call output ( '    \parbox{' // width // '}{{\bf \Large %' )
  call output ( trim(in_file) // '%', und=.true. )
  call output ( '    }\newline%' )
  call output ( '    Page \thepage\ of \pageref{lastpage} \hfill ' // &
  &             now // '}}%' )
  call output ( 'renewcommand{\@oddfoot}{}%', tex=.true. )
  call output ( 'renewcommand{\@evenfoot}{}%', tex=.true. )
  call output ( '}%' )
  call output ( 'makeatother', tex=.true. )
  call output ( 'pagestyle{twolines}', tex=.true. )
  call output ( '' )
  call output ( 'parindent 0pt \parskip 5pt', tex=.true. )
  if ( sx == 1 ) then
    call output ( 'lstset{language=[90]Fortran,', tex=.true. )
    call output ( '  basicstyle=\ttfamily\footnotesize,' )
    call output ( '  keywordstyle=\bfseries,' )
    call output ( '  labelstyle=\tiny,labelstep=' // number_step // '}' )
  end if

! Copy the program text
  do
    if ( in_unit >= 0 ) then
      read ( in_unit, '(a)', iostat=iostat ) line
    else
      read ( *, '(a)', iostat=iostat ) line
    end if
    if ( iostat /= 0 ) exit
    i = fnb(line) ! fnb == first nonblank
    select case ( state )
    case ( 0 ) ! begin
      if ( line(i:i+1) /= '!{' ) then
        state = 1
        if ( sx == 2 ) call output ( '{\tt', adv='no' )
        call output ( start_code(sx), tex=.true. )
        call output ( trim(line) )
      else
        if ( box ) &
          call output ( 'framebox[' // width // '][l]{\parbox{' // &
                        widthp // '}{', tex=.true. )
        call output ( trim(line(i+2:)) )
        state = 2
      end if
    case ( 1 ) ! doing code
      if ( line(i:i+1) /= '!{' ) then
        call output ( trim(line) )
      else
        call output ( stop_code(sx), tex=.true. )
        if ( box ) &
          call output ( 'framebox[' // width // '][l]{\parbox{' // &
                        widthp // '}{', tex=.true. )
        call output ( trim(line(i+2:)) )
        state = 2
      end if
    case ( 2 ) ! doing LaTeX
      if ( line(i:i) == '!' ) then
        call output ( trim(line(i+1:)) )
      else
        if ( box ) call output ( '}}' )
        state = 1
        if ( sx == 2 ) call output ( '{\tt', adv='no' )
        call output ( start_code(sx), tex=.true. )
        call output ( trim(line) )
      end if
    end select
  end do

  if ( state == 1 ) call output ( stop_code(sx), tex=.true. )
    
! Output the ending stuff
  call output ( 'label{lastpage}', tex=.true. )
  call output ( 'end{document}', tex=.true. )

  stop

contains
! ----------------------------------------------------     FNB     -----
  integer function FNB ( TEXT )
  ! Find the position of the first nonblank in TEXT
    character(len=*), intent(in) :: TEXT
    integer :: I
    do i = 1, len(text)
      if ( text(i:i) /= ' ' ) then
        fnb = i
        return
      end if
    end do
    fnb = 0
    return
  end function FNB
! -------------------------------------------------     OUTPUT     -----
  subroutine OUTPUT ( TEXT, ADV, UND, TEX )
  ! Output text to OUT_UNIT if it's positive, else to unit *.  Put a "\"
  ! at the beginning if TEX is present (no matter with what value).
    character(len=*), intent(in) :: TEXT
    character(len=*), intent(in), optional :: ADV
    logical, intent(in), optional :: UND, TEX
    character(len=3) ::  ADVANCE   ! for I/O
    integer :: I                   ! Subscript / do inductor
    logical :: UNDER               ! .true. => put "\" before "_"

    advance = 'yes'
    if ( present(adv) ) advance = adv
    under = sx == 1
    i = fnb(text)
    if ( i > 0 ) then
      if ( text(i:i) /= '!' ) under = .false.
    end if
    if ( present(und) ) under = und
    if ( out_unit >= 0 ) then
      if ( present(tex) ) write ( out_unit, '(a)', advance='no' ) '\'
      do i = 1, len(text)
        if ( under .and. text(i:i) == '_' ) &
          write ( out_unit, '(a)', advance='no' ) '\'
        write ( out_unit, '(a)', advance='no' ) text(i:i)
      end do
      write ( out_unit, '(a)' )
    else
      if ( present(tex) ) write ( *, '(a)', advance='no' ) '\'
      do i = 1, len(text)
        if ( under .and. text(i:i) == '_' ) &
          write ( *, '(a)', advance='no' ) '\'
        write ( *, '(a)', advance='no' ) text(i:i)
      end do
      write ( *, '(a)' )
    end if
    return
  end subroutine OUTPUT

end program F90TEX

! $Log$
! Revision 1.2  2001/06/01 21:49:23  vsnyder
! Fix up CVS log variable
!
! Revision 1.1  2001/06/01 21:47:22 vsnyder
! Initial commit
!
