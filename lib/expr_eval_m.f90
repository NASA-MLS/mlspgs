! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

! Contact: Van Snyder, vsnyder@math.jpl.nasa.gov

!=============================================================================
module Expr_Eval_M

! Evaluate a real scalar expression, given as an array of characters.

! No continuation mark is necessary (or indeed recognized) to proceed from
! one array element to the next.  Otherwise, the syntax is the same as
! free-form Fortran 95.  The numeric intrinsic functions abs, aint, anint,
! ceiling, dim, floor, max, min, mod, modulo and sign are provided, except
! that min and max can only take exactly two arguments.  All of the
! mathematical intrinsic functions are provided; atan can also take two
! arguments. All of the floating-point manipulation functions are
! provided; for scale and set_exponent, nint is applied to the second
! argument.  The numeric inquiry functions are not provided.

! LL(1) Action grammar for expressions parsed by Expr_Eval (actions in {}):

! E ::= '+' Term E_List(T) { E := E_List(Term) } 'end'
!   ::= '-' Term E_List(T) { E := E_List(-Term) } 'end'
!   ::= Term E_List(T) { E := E_List(Term) } 'end'
! E_List(a) ::= '+' Term { b := a+Term } E_List(b) { E_List_1 := E_List_2(b) }
!           ::= '-' Term { b := a-Term } E_List(b) { E_List_1 := E_List_2(b) }
!           ::= \lambda { E_List := a }
! Term ::= Factor T_List(Factor) { Term := T_List(Factor) }
! T_List(a) ::= '*' Factor { b := a*Factor } T_List(b) { T_List_1 := T_List_2(b) }
!           ::= '/' Factor { b := a/Factor } T_List(b) { T_List_1 := T_List_2(b) }
!           ::= \lambda { T_List := a }
! Factor ::= Primary F_List(Primary) { Factor ::= F_List(Primary) }
! F_List(a) ::= '**' Factor { F_List := a**Factor }
!           ::= \lambda { F_List := a }
! Primary ::= 'number' { Primary := value(number) }
!   ::= '(' E ')' { Primary := E }
!   :: = 'name' '(' E [',' E] ')' {Primary := eval function 'name' with E[,E]}

! The module variable "Error" is set if an error occurs.

  implicit NONE

  private

  public :: Expr_Eval

  logical, public :: Error ! True means an error occurred

  logical, public :: Trace = .false.

  integer, parameter :: NameLen = 12
  type :: Token
    double precision :: Value
    character(len=nameLen) :: Name
    integer :: Class
  end type Token

  ! Token classes:
  integer, parameter :: Add = 1
  integer, parameter :: Begin = add + 1
  integer, parameter :: Comma = begin + 1
  integer, parameter :: Divide = Comma + 1
  integer, parameter :: End = divide + 1
  integer, parameter :: Expon = end + 1
  integer, parameter :: LeftPar = expon + 1
  integer, parameter :: Multiply = leftPar + 1
  integer, parameter :: Name = multiply + 1
  integer, parameter :: Number = name + 1
  integer, parameter :: RightPar = number + 1
  integer, parameter :: Subtract = rightPar + 1

  ! Token names, indexed by token classes:
  character(len=2), parameter :: TokenName(add:subtract) = (/ &
    & '+ ', '|-', ', ', '/ ', '-|', '**', '( ', '* ', 'X ', '# ', ') ', '- ' /)

  integer :: Depth  ! to print leading dots for trace

  logical :: NEED   ! "Need a token"

  type(token) :: The_Token    ! Last token examined

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains

  double precision function Expr_Eval ( The_Expr, C, L )
    character(len=*), intent(in) :: The_Expr(:)
    integer, intent(out), optional :: C ! Last character position examined in line L of The_Expr   
    integer, intent(out), optional :: L ! Line of The_Expr currently being examined.               
    integer :: MyC, MyL

    myC = 0
    myL = 1
    depth = 0
    error = .false.
    need = .true.
    the_token = token(0.0d0,'      ',begin)
    expr_eval = expr(the_expr,myC,myL)
    if ( need ) call lexer ( the_token, the_expr, c, l )
    if ( the_token%class /= end ) error = .true.
    if ( present(c) ) c = myC
    if ( present(l) ) l = myL
  end function Expr_Eval

  ! Private procedures in the same order as LHS's in the grammar, except that
  ! lexer, trace_begin and trace_end are at the end.

  recursive double precision function Expr ( The_Expr, C, L ) result ( E )
    character(len=*), intent(in) :: The_Expr(:)
    integer, intent(inout) :: C    ! Last character position examined.
    integer, intent(inout) :: L    ! Line of The_Expr currently being examined.
    double precision :: T          ! Result of Term

    if ( trace ) call trace_begin ( 'Expr', the_expr, c, l )
    if ( need ) call lexer ( the_token, the_expr, c, l )
    if ( .not. error ) then
      select case ( the_token%class )
      case ( add )
        need = .true.
        t = term(the_expr,c,l)
      case ( subtract )
        need = .true.
        t = -term(the_expr,c,l)
      case default
        t = term(the_expr,c,l)
      end select
      e = e_list(the_expr,c,l,t)
    end if
    if ( trace ) call trace_end ( 'Expr', e )
  end function Expr

  recursive double precision function E_List ( The_Expr, C, L, T ) result ( EL )
    character(len=*), intent(in) :: The_Expr(:)
    integer, intent(inout) :: C    ! Last character position examined.
    integer, intent(inout) :: L    ! Line of The_Expr currently being examined.
    double precision, intent(in) :: T   ! Input value
    if ( trace ) call trace_begin ( 'E_List', the_expr, c, l )
    el = t
    do while ( .not. error )
      if ( need ) call lexer ( the_token, the_expr, c, l )
      select case ( the_token%class )
      case ( add )
        need = .true.
        el = el + term(the_expr,c,l)
      case ( subtract )
        need = .true.
        el = el - term(the_expr,c,l)
      case default
        exit
      end select
    end do
    if ( trace ) call trace_end ( 'E_List', el )
  end function E_List

  recursive double precision function Term ( The_Expr, C, L ) result ( T )
    character(len=*), intent(in) :: The_Expr(:)
    integer, intent(inout) :: C    ! Last character position examined.
    integer, intent(inout) :: L    ! Line of The_Expr currently being examined.
    if ( trace ) call trace_begin ( 'Term', the_expr, c, l )
    t = t_list(the_expr,c,l,factor(the_expr,c,l))
    if ( trace ) call trace_end ( 'Term', t )
  end function Term

  recursive double precision function T_List ( The_Expr, C, L, F ) result ( TL )
    character(len=*), intent(in) :: The_Expr(:)
    integer, intent(inout) :: C    ! Last character position examined.
    integer, intent(inout) :: L    ! Line of The_Expr currently being examined.
    double precision, intent(in) :: F   ! Input value
    if ( trace ) call trace_begin ( 'T_List', the_expr, c, l, 'F', f )
    tl = f
    do while ( .not. error )
      if ( need ) call lexer ( the_token, the_expr, c, l )
      select case ( the_token%class )
      case ( divide )
        need = .true.
        tl = tl / factor(the_expr,c,l)
      case ( multiply )
        need = .true.
        tl = tl * factor(the_expr,c,l)
      case default
        exit
      end select
    end do
    if ( trace ) call trace_end ( 'T_List', tl )
  end function T_List

  recursive double precision function Factor ( The_Expr, C, L ) result ( F )
    character(len=*), intent(in) :: The_Expr(:)
    integer, intent(inout) :: C    ! Last character position examined.
    integer, intent(inout) :: L    ! Line of The_Expr currently being examined.
    if ( trace ) call trace_begin ( 'Factor', the_expr, c, l )
    f = f_list(the_expr,c,l,primary(the_expr,c,l))
    if ( trace ) call trace_end ( 'Factor', f )
  end function Factor

  recursive double precision function F_List ( The_Expr, C, L, P ) result ( FL )
    character(len=*), intent(in) :: The_Expr(:)
    integer, intent(inout) :: C    ! Last character position examined.
    integer, intent(inout) :: L    ! Line of The_Expr currently being examined.
    double precision, intent(in) :: P   ! Input value
    if ( trace ) call trace_begin ( 'F_List', the_expr, c, l, 'P', p )
    if ( need ) call lexer ( the_token, the_expr, c, l )
    fl = p
    if ( the_token%class == expon ) then
      need = .true.
      fl = fl ** factor(the_expr,c,l)
    end if
    if ( trace ) call trace_end ( 'F_List', fl )
  end function F_List

  recursive double precision function Primary ( The_Expr, C, L ) result ( P )
    character(len=*), intent(in) :: The_Expr(:)
    integer, intent(inout) :: C    ! Last character position examined.
    integer, intent(inout) :: L    ! Line of The_Expr currently being examined.
    double precision :: E1, E2     ! Function args
    character(len=nameLen) :: Fname     ! Function name
    integer :: I                   ! Loop inductor, subscript
    if ( trace ) call trace_begin ( 'Primary', the_expr, c, l )
    if ( need .and. the_token%class /= end ) &
      & call lexer ( the_token, the_expr, c, l )
    select case ( the_token%class )
    case ( leftPar )
      need = .true.
      p = expr(the_expr,c,l)
      if ( need ) call lexer ( the_token, the_expr, c, l )
      if ( the_token%class /= rightPar ) go to 999
      need = .true.
    case ( number )
      p = the_token%value
      need = .true.
    case ( name )
      fname = the_token%name
      do i = 1, len(fname) ! Convert name to lower case
        if ( fname(i:i) >= 'A' .and. fname(i:i) <= 'Z' ) &
          & fname(i:i) = char(ichar(fname(i:i)) + &
            & ichar('a') - ichar('A'))
      end do
      need = .true.
      call lexer ( the_token, the_expr, c, l )
      if ( the_token%class /= leftPar ) go to 999
      need = .true.
      e1 = expr(the_expr,c,l)
      i = 1 ! one arg
      if ( need ) call lexer ( the_token, the_expr, c, l )
      select case ( the_token%class )
      case ( comma )
        i = 2 ! two args
        need = .true.
        e2 = expr(the_expr,c,l)
      case ( rightPar )
      case default
        go to 999
      end select
      if ( need ) call lexer ( the_token, the_expr, c, l )
      if ( the_token%class /= rightPar ) go to 999
      need = .true.
      select case ( fname )
      case ( 'atan2', 'dim', 'max', 'min', 'mod', 'modulo', 'nearest', 'scale', &
        & 'set_exponent', 'sign' )
        if ( i /= 2 ) go to 999
      case ( 'atan' ) ! don't care, handled below
      case default
        if ( i /= 1 ) go to 999
      end select
      select case ( fname )
      case ( 'abs' )
        p = abs(e1)
      case ( 'acos' )
        p = acos(e1)
      case ( 'aint' )
        p = aint(e1)
      case ( 'anint' )
        p = anint(e1)
      case ( 'asin' )
        p = asin(e1)
      case ( 'atan' )
        if ( i == 1 ) then
          p = atan(e1)
        else
          p = atan2(e1,e2)
        end if
      case ( 'atan2' )
        p = atan2(e1,e2)
      case ( 'ceiling' )
        p = ceiling(e1)
      case ( 'cos' )
        p = cos(e1)
      case ( 'cosh' )
        p = cosh(e1)
      case ( 'dim' )
        p = dim(e1,e2)
      case ( 'exp' )
        p = exp(e1)
      case ( 'exponent' )
        p = exponent(e1)
      case ( 'floor' )
        p = floor(e1)
      case ( 'fraction' )
        p = fraction(e1)
      case ( 'log' )
        p = log(e1)
      case ( 'log10' )
        p = log10(e1)
      case ( 'max' )
        p = max(e1,e2)
      case ( 'min' )
        p = min(e1,e2)
      case ( 'mod' )
        p = mod(e1,e2)
      case ( 'modulo' )
        p = modulo(e1,e2)
      case ( 'nearest' )
        p = nearest(e1,e2)
      case ( 'rrspacing' )
        p = rrspacing(e1)
      case ( 'scale' )
        p = scale(e1,nint(e2))
      case ( 'set_exponent' )
        p = set_exponent(e1,nint(e2))
      case ( 'sign' )
        p = sign(e1,e2)
      case ( 'sin' )
        p = sin(e1)
      case ( 'sinh' )
        p = sinh(e1)
      case ( 'spacing' )
        p = spacing(e1)
      case ( 'sqrt' )
        p = sqrt(e1)
      case ( 'tan' )
        p = tan(e1)
      case ( 'tanh' )
        p = tanh(e1)
      case default
        go to 999
      end select
    case default
      error = .true.
    end select
99  if ( trace ) call trace_end ( 'Primary', p )
    return
999 error = .true.
    p = -huge(p)
    go to 99
 end function Primary

  subroutine Lexer ( The_Token, The_Expr, C, L )
    type(token), intent(out) :: The_Token
    character(len=*), intent(in) ::The_Expr(:)
    integer, intent(inout) :: C    ! Last character position examined.
    integer, intent(inout) :: L    ! Line of The_Expr currently being examined.
    integer :: C1   ! Beginning of token if that's interesting
    integer :: I, IOSTAT
    need = .false.
    the_token = token(0.0d0,'      ',begin) ! just so it always has a value
    do
      c = c + 1
      if ( c > len(the_expr) ) then
        c = 1
        l = l + 1
      end if
      if ( l > size(the_expr) ) then
        the_token%class = end
        go to 99
      end if
      if ( the_expr(l)(c:c) /= ' ' ) exit
    end do
    select case ( the_expr(l)(c:c) )
    case ( '+' )
      the_token%class = add
    case ( '-' )
      the_token%class = subtract
    case ( '*' )
      the_token%class = multiply
      if ( c < len(the_expr) ) then
        if ( the_expr(l)(c+1:c+1) == '*' ) then
          c = c + 1
          the_token%class = expon
        end if
      end if
    case ( '/' )
      the_token%class = divide
    case ( '(' )
      the_token%class = leftPar
    case ( ')' )
      the_token%class = rightPar
    case ( ',' )
      the_token%class = comma
    case ( '.', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' )
      the_token%class = number
      c1 = c
      do
        c = c + 1
        if ( c > len(the_expr) ) go to 88
        select case ( the_expr(l)(c:c) )
        case ( '.', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' )
        case ( 'd', 'D', 'e', 'E' )
          exit
        case default
          go to 88
        end select
      end do
      c = c + 1
      select case ( the_expr(l)(c:c) )
      case ( '+', '-' )
      end select
      do
        c = c + 1
        if ( c > len(the_expr) ) go to 88
        select case ( the_expr(l)(c:c) )
        case ( '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' )
        case default
          exit
        end select
      end do
88    c = c - 1
      read ( the_expr(l)(c1:c), *, iostat=iostat ) the_token%value
      error = iostat /= 0
      if ( error ) the_token%value = -huge(the_token%value)
    case ( 'a':'z', 'A':'Z' )
      the_token%class = name
      i = 1
      the_token%name(i:i) = the_expr(l)(c:c)
      do
        c = c + 1
        if ( c > len(the_expr) ) exit
        select case ( the_expr(l)(c:c) )
        case ( 'a':'z', 'A':'Z', '0':'9', '_' )
          i = i + 1
          if ( i > len(the_token%name) ) exit
          the_token%name(i:i) = the_expr(l)(c:c)
        case default
          exit
        end select
      end do
      c = c - 1
    case default
      the_token%value = -huge(the_token%value)
      error = .true.
    end select
99  if ( trace ) then
      do i = 1, depth
        write ( *, '(a)', advance='no' ) '.'
      end do
      if ( the_token%class == name ) then
        write ( *, * ) 'The_Token%name = ', trim(the_token%name)
      else
        write ( *, * ) 'The_Token%value =', the_token%value, &
          & ', %class = ', trim(tokenName(the_token%class))
      end if
    end if
  end subroutine Lexer

  subroutine Trace_Begin ( Name, The_Expr, C, L, InputName, InputValue )
    character(len=*), intent(in) :: Name
    character(len=*), intent(in) ::The_Expr(:)
    integer, intent(in) :: C       ! Last character position examined.
    integer, intent(in) :: L       ! Line of The_Expr currently being examined.
    character(len=*), intent(in), optional :: InputName
    double precision, intent(in), optional :: InputValue
    integer :: I, MyC, MyL
    myC = c
    myL = l
    if ( l > size(the_expr) ) then
      myL = size(the_expr)
      myC = len(the_expr) + 1
    end if
    do i = 1, depth
      write ( *, '(a)', advance='no' ) '.'
    end do
    depth = depth + 1
    if ( present(inputName) ) then
      write ( *, * ) ' Enter ', trim(name), ' with ', trim(the_expr(myL)), &
        & ' and ', trim(inputName), ' = ', inputValue, ' at ', myC
    else
      write ( *, * ) ' Enter ', trim(name), ' with ', trim(the_expr(myL)), &
        & ' at ', myC
    end if
  end subroutine Trace_Begin

  subroutine Trace_End ( Name, Result )
    character(len=*), intent(in) :: Name
    double precision, intent(in) :: Result
    integer :: I
    depth = depth - 1
    do i = 1, depth
      write ( *, '(a)', advance='no' ) '.'
    end do
    write ( *, * ) ' Exit ', trim(name), ' with result ', result
  end subroutine Trace_End

end module Expr_Eval_M

! $Log$
! Revision 2.1  2001/12/17 23:39:55  vsnyder
! Initial commit
!
