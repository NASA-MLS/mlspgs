! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

! Contact: Van Snyder, vsnyder@math.jpl.nasa.gov

!=============================================================================
program Test_Expr

  use Expr_Eval_m, only: Expr_Eval, Error, Help, Trace

  implicit NONE

  character(len=127) :: Line(1)
  character(len=127) :: Me
  integer :: C, L
  logical :: Many = .false.   ! Set by -m option, overrides "once"
  logical :: Once = .false.   ! Set if an expression is provided on the
                              ! command line

  double precision :: Value

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

! trace = .true.
  ! Process command line arguments.
  ! If you don't have getarg, comment this out.
  call getarg ( 0, me )
  c = 0
  do 
    c = c + 1
    call getarg ( c, line )
    if ( line(1)(1:3) == '-h' .or. line(1)(1:3) == '-? ' ) then
      write ( *, '(3a)' ) 'Usage: ', trim(me), ' [options] [expr]'
      write ( *, '(a)' )  ' Options: -h => This stuff'
      write ( *, '(2a)' ) '          -m => Do Many expressions even if there''s', &
        &                 ' one on the command line'
      write ( *, '(a)' )  '          -t => Trace execution'      
      call help
      stop
    end if
    if ( line(1)(1:3) == '-m ' ) then
      many = .true.
    else if ( line(1)(1:3) == '-t ' ) then
      trace = .true.
    else
      once = .not. many .and. (line(1) /= ' ')
      exit
    end if
  end do
  ! End of command line processing
  do
    if ( line(1) == ' ' ) read ( *, '(a)', end=99 ) line
    value = expr_eval ( line, c, l )
    if ( error ) then
      print *, trim(me), ': Syntax error at line ', l, ', column ', c
    else
      print *, value
    end if
    if ( once ) exit
    line = ' '
  end do
99 stop
end program Test_Expr 

! $Log$
