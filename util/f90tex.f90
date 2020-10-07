program F90TEX

!{ Read a Fortran 90 program.  Write a LaTeX file.  Start the output with
! the usual {\tt$\backslash$documentclass} stuff.  In the body of the
! document, wrap code with {\tt$\backslash$begin}\{\emph{verbatim}\} and
! {\tt$\backslash$end}\{\emph{verbatim}\} where \emph{verbatim} is given
! by the variable {\tt CODE\_STYLE}.  When a comment begins with !\{, end
! the {\tt$\backslash$begin}\{\emph{verbatim}\} wrapper, and start it
! again at the end of the comment.  \emph{verbatim} is {\tt verbatim} if
! the {\tt -v} option is specified (default) or {\tt lstlisting} if the
! {\tt -l} option is specified.
!
! A block beginning with !\{ is put in a box if the {\tt -b} option is
! specified (default) and not put in a box if the {\tt -nb} option is
! specified.  If {\tt -nb} is not specified, {\tt$\backslash$newpage} only
! works if it's the first line of the !\{ block, and there's nothing else on
! the line.

!>> 2013-08-08 F90TEX WV Snyder Remove linenumbers stuff using lstlistings
!>> 2013-08-08 F90TEX WV Snyder Remove labelstyle keyword for lstlistings
!>> 2013-08-06 F90TEX WV Snyder Remove dependence on machine module
!>> 2012-02-09 F90TEX WV Snyder Stuff to make \cleardoublepage work
!>> 2010-08-21 F90TEX WV Snyder Added -u option for usepackage
!>> 2003-10-28 F90TEX WV Snyder Stuff to make \newpage work
!>> 2000-04-07 F90TEX WV Snyder Original code

! =====     Declarations     ===========================================

  implicit NONE

  logical :: BOX = .true.                         ! Put a box around !{ TeX?
  character(len=8) :: DATE                        ! of program execution
  character(127) :: Errmsg                        ! from I/O
  character(*), parameter :: FilSep = '/'         ! file separator in path names
  integer :: I, J                                 ! Subscript / do inductor
  character(len=132) :: IN_FILE, OUT_FILE         ! File names
  integer :: IN_UNIT = 10, OUT_UNIT = 11          ! Unit numbers
  integer :: IOSTAT                               ! I/O status
  character(len=132) :: LINE                      ! From input
  integer :: LineNo = 0                           ! Of input
  character(len=23) :: NOW                        ! Date and time, formatted
  character :: NUMBER_STEP = '1'                  ! Line number step
  integer :: Number_Packages = 0                  ! Number of packages to use
  character(len=32) :: Packages(99)               ! LaTeX packages to use
  logical :: Running = .true.                     ! Running line numbers
                                                  !  else pagewise
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

!---------------------------- RCS Module Info --------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!-----------------------------------------------------------------------

! =====     Executable Statements     ==================================

  call date_and_time ( date, time )
  now = date(1:4) // '/' // date(5:6) // '/' // date(7:) // ' ' // &
        time(1:2) // ':' // time(3:4) // ':' // time(5:)
! Get input and output file names

  i = 1
  do
    call get_command_argument ( i, in_file )
    if ( in_file(1:1) /= '-' ) exit
    if ( in_file(1:2) == '- ' ) then
      i = i + 1
      call get_command_argument ( i, in_file )
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
    else if ( in_file(1:2) == '-n' ) then
      if ( in_file(3:) == '' ) then
        i = i + 1
        call get_command_argument ( i, number_step )
      else
        number_step = in_file(3:)
      end if
    else if ( in_file(1:3) == '-p ' ) then
      running = .false.
    else if ( in_file(1:2) == '-u' ) then
      if ( number_packages == ubound(packages,1) ) then
        print *, 'More than ', ubound(packages,1), ' ignored'
      else
        number_packages = number_packages + 1
        if ( in_file(3:) == '' ) then
          i = i + 1
          call get_command_argument ( i, packages(number_packages) )
        else
          packages(number_packages) = in_file(3:)
        end if
      end if
    else
      call get_command_argument ( 0, in_file )
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
      write (*,*) ' Options, [ ] indicates optional blanks:'
      write (*,*) '  -b => Put a box around !{ TeX stuff (default).'
      write (*,*) '  -nb => Don''t put a box around !{ TeX stuff.'
      write (*,*) '  -v => Use verbatim environment for code (default).'
      write (*,*) '  -l => Use lstlisting environment from the listings package.'
      write (*,*) '        This package has numerous customization features;'
      write (*,*) '        You might wish to edit the LaTeX output to exploit them.'
      write (*,*) '        This might not work if your code has "\" in it.'
      write (*,*) '  -n[ ]# => Line number step is # (default ', &
      &                      number_step, ', 0 => no line numbers).'
      write (*,*) '  -p => Pagewise line numbers, default is running numbers.'
      write (*,*) '  -u[ ]package => Use LaTeX package; this option'
      write (*,*) '                  can appear up to 99 times.'
      write (*,*) '  -<anything else> => This output.'
      write (*,*) '  -  => No more options.'
      write (*,*) ' 8 August 2013'
      stop
    end if
    i = i + 1
  end do
  call get_command_argument ( i+1, out_file )

  if ( in_file /= ' ' ) then
    open(in_unit, file=in_file, status='old', iostat=iostat, iomsg=errmsg)
    if ( iostat /= 0 ) then
      print '("Unable to open input file ",a)', trim(in_file)
      print '("Status = ",i0," Message = ")', iostat, trim(errmsg)
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
      print '("Unable to open output file ",a)', trim(out_file)
      print '("Status = ",i0," Message = ")', iostat, trim(errmsg)
      stop
    end if
  else
    in_unit = -1         ! Standard input
    in_file = 'Standard Input'
    out_unit = -1        ! Standard output
  end if

! Replace any \ in the input file name with /, so as not to confuse LaTeX
! when it tries to produce the page header.

  do i = 1, len(in_file)
    if ( in_file(i:i) == '\' ) in_file(i:i) = '/'
  end do

! Output the beginning stuff

  call output ( 'documentclass[11pt,twoside]{article}', tex=.true. )
  call output ( 'usepackage[fleqn]{amsmath}', tex=.true. )
  if ( sx == 1 ) call output ( 'usepackage{listings}', tex=.true. )
  do i = 1, number_packages
    call output ( 'usepackage{' // trim(packages(i)) // '}', tex=.true. )
  end do
  call output ( 'textwidth ' // width, tex=.true. )
  call output ( 'oddsidemargin -0.25in', tex=.true. )
  call output ( 'evensidemargin -0.25in', tex=.true. )
  call output ( 'topmargin -0.5in', tex=.true. )
  call output ( 'textheight 9.25in', tex=.true. )
  if ( sx == 2 .and. number_step /= '0' ) then
    call output ( 'usepackage{lineno}', tex=.true. )
    call output ( 'def\linenumberfont{\normalfont\footnotesize\sffamily}', tex=.true. )
  end if
  call output ( 'usepackage[strings]{underscore}', tex=.true. )
  call output ( '' )
  call output ( 'begin{document}', tex=.true. )
  if ( sx == 2 .and. number_step /= '0' ) then
    if ( running ) then
      call output ( 'runninglinenumbers', tex=.true. )
    else
      call output ( 'pagewiselinenumbers', tex=.true. )
    end if
    call output ( 'modulolinenumbers[' // number_step // ']', tex=.true. )
  end if
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
!   It looks like the labelstyle and labelstep keywords have been removed
!   call output ( '  labelstyle=\tiny,labelstep=' // number_step )
    call output ( '}' )
  end if

! Copy the program text
  do
    lineNo = lineNo + 1
    if ( in_unit >= 0 ) then
      read ( in_unit, '(a)', iostat=iostat ) line
    else
      read ( *, '(a)', iostat=iostat ) line
    end if
    if ( iostat /= 0 ) exit
    i = max(verify(line,' '),1) ! first nonblank, or first character if all blank
    select case ( state )
    case ( 0 ) ! begin
      if ( line(i:i+1) /= '!{' ) then
        state = 1
        if ( sx == 2 ) call output ( '{\tt', adv='no' )
        call output ( start_code(sx), tex=.true. )
        call output ( trim(line) )
      else
        if ( adjustl(line(i+2:)) == '\newpage' .or. &
           & adjustl(line(i+2:)) == '\cleardoublepage' ) then
          state = 3
        else if ( box .and. adjustl(line(i+2:)) /= '\newpage' .or. &
                & box .and. adjustl(line(i+2:)) /= '\cleardoublepage' ) then
          call output ( 'framebox[' // width // '][l]{\parbox{' // &
                        widthp // '}{', tex=.true. )
          state = 2
        end if
        call output ( trim(line(i+2:)) )
      end if
    case ( 1 ) ! doing code
      if ( line(i:i+1) /= '!{' ) then
        call output ( trim(line) )
      else
        call output ( stop_code(sx), tex=.true. )
        if ( sx == 2 .and. number_step /= '0' ) &
          & call output ( 'nolinenumbers', tex=.true. )
        if ( adjustl(line(i+2:)) == '\newpage' .or. &
           & adjustl(line(i+2:)) == '\cleardoublepage' ) then
          state = 3
        else
          if ( box ) &
            & call output ( 'framebox[' // width // '][l]{\parbox{' // &
                          widthp // '}{', tex=.true. )
          state = 2
        end if
        call output ( trim(line(i+2:)) )
      end if
    case ( 2, 3 ) ! doing LaTeX, 3 = first line is \newpage
      if ( line(i:i) == '!' ) then
        if ( state == 3 ) &
          & call output ( 'framebox[' // width // '][l]{\parbox{' // &
                        widthp // '}{', tex=.true. )
        call output ( trim(line(i+1:)) )
        state = 2
      else
        if ( box .and. state == 2 ) call output ( '}}' )
        state = 1
        call output ( '' ) ! without this, the previous paragraph gets line
                           ! numbers if there is more than one paragraph in
                           ! block
        if ( sx == 2 .and. number_step /= '0' ) then
          write ( now, '("linenumbers[",i0,"]")' ) lineNo
          call output ( trim(now), tex=.true. )
        end if
        if ( sx == 2 ) call output ( '{\tt', adv='no' )
        call output ( start_code(sx), tex=.true. )
        call output ( trim(line) )
      end if
    end select
  end do

  if ( state == 1 ) then
    call output ( stop_code(sx), tex=.true. )
  else if ( state > 1 ) then
    if ( box .and. state == 2 ) call output ( '}}' )
    call output ( '' ) ! without this, the previous paragraph gets line
                           ! numbers if there is more than one paragraph in
                           ! block
  end if
    
! Output the ending stuff
! The only effect of the "vspace" line is that "lastpage" is defined
! if a line of text after the last line would be on a new page.
  call output ( 'vspace*{-1in}', tex=.true. )
  call output ( 'label{lastpage}', tex=.true. )
  call output ( 'end{document}', tex=.true. )

  stop

!{\newpage
contains
! -------------------------------------------------     OUTPUT     -----
  subroutine OUTPUT ( TEXT, ADV, UND, TEX )
  ! Output text to OUT_UNIT if it's positive, else to unit *.  Put a
  ! "\" at the beginning if TEX is present (no matter with what value).
    character(len=*), intent(in) :: TEXT
    character(len=*), intent(in), optional :: ADV
    logical, intent(in), optional :: UND, TEX
    character(len=3) ::  ADVANCE   ! for I/O
    integer :: I                   ! Subscript / do inductor
    logical :: UNDER               ! .true. => put "\" before "_"

    advance = 'yes'
    if ( present(adv) ) advance = adv
    under = sx == 1
    i = verify(text,' ') ! first nonblank
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
      write ( out_unit, '(a)', advance=advance )
    else
      if ( present(tex) ) write ( *, '(a)', advance='no' ) '\'
      do i = 1, len(text)
        if ( under .and. text(i:i) == '_' ) &
          write ( *, '(a)', advance='no' ) '\'
        write ( *, '(a)', advance='no' ) text(i:i)
      end do
      write ( *, '(a)', advance=advance )
    end if
    return
  end subroutine OUTPUT

end program F90TEX

! $Log$
! Revision 1.23  2016/03/28 22:18:39  vsnyder
! Keep substring bounds for Line in range
!
! Revision 1.22  2013/08/09 00:40:39  vsnyder
! Some cannonball polishing
!
! Revision 1.21  2013/08/08 20:48:50  vsnyder
! Don't emit linenumbers commands if -l selected
!
! Revision 1.20  2013/08/08 20:11:14  vsnyder
! Remove labelstyle and labelstep from lstlisting style
!
! Revision 1.19  2013/08/07 20:35:42  vsnyder
! Replace back slash in input file name with ordinary slash.  Otherwise,
! LaTeX thinks words in the path after the back slash are commands, and
! gets confused while trying to produce the page heading.
!
! Revision 1.18  2013/08/07 20:19:26  vsnyder
! Repair -nb option
!
! Revision 1.17  2013/08/06 23:40:32  vsnyder
! Explicitly set line number after typeset box
!
! Revision 1.16  2013/08/06 23:14:31  vsnyder
! Remove dependence on machine module
!
! Revision 1.15  2012/02/09 23:50:05  vsnyder
! Go back to using getarg; let machine turn that into get_command_argument
!
! Revision 1.12  2010/08/22 00:45:24  vsnyder
! Add -u option to include LaTeX packages
!
