! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module TRACE_M

  use LEXER_CORE, only: PRINT_SOURCE
  use OUTPUT_M, only: OUTPUT
  use TIME_M, only: TIME_NOW
  use TREE, only: DUMP_TREE_NODE, SOURCE_REF

  private
  public :: Trace_Begin, Trace_End
  integer, public, save :: DEPTH   ! Depth in tree.  Used for trace printing.

  integer, parameter :: ClockStackMax = 100
  real :: ClockStack(0:clockStackMax) = 0.0

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================
! --------------------------------------------------  TRACE_BEGIN  -----
  subroutine TRACE_BEGIN ( NAME, ROOT, INDEX )
  ! Print "ENTER NAME with ROOT = <node_id(root)>" with DEPTH dots in
  ! front.  Increment DEPTH.
    character(len=*), intent(in) :: NAME
    integer, intent(in), optional :: ROOT
    integer, intent(in), optional :: INDEX
    integer :: I              ! Loop inductor
    character(len=10) :: Now  ! For Date_and_time
    if ( present(root) ) then
      call output ( root, 4 ); call output ( ': ' )
    else
      call output ( '      ' )
    end if
    do i = 1, depth
      call output ( '.' )
    end do
    call output ( 'Enter ' ); call output ( name )
    if ( present(index) ) call output( index )
    call date_and_time ( time=now )
    call output ( ' at ' // now(1:2) // ':' // now(3:4) // ':' // now(5:) )
    if ( present(root) ) then
      call output ( ' with ' );
      call dump_tree_node ( root, 0 )
      call output ( ' at ' )
      call print_source ( source_ref(root), advance='yes' )
    else
      call output ( '', advance='yes' )
    end if
    if ( depth >= 0 .and. depth < clockStackMax ) then
      call time_now ( clockStack(depth) )
      clockStack(depth+1) = 0.0
    end if
    depth = depth + 1
  end subroutine TRACE_BEGIN
! --------------------------------------------------    TRACE_END  -----
  subroutine TRACE_END ( NAME, INDEX )
  ! Decrement DEPTH.  Print "EXIT NAME with DEPTH dots in front.
    character(len=*), intent(in) :: NAME
    integer, intent(in), optional :: INDEX
    integer :: I              ! Loop inductor
    character(len=10) :: Now  ! For Date_and_time
    integer :: Values(8)      ! For Date_and_time
    real :: T                 ! For timing
    depth = depth - 1
    call output ( '      ' )
    do i = 1, depth
      call output ( '.' )
    end do
    call output ( 'Exit ' ); call output ( name )
    if ( present(index) ) call output( index )
    call date_and_time ( time=now, values=values )
    call output ( ' at ' // now(1:2) // ':' // now(3:4) // ':' // now(5:) )
    if ( depth >= 0 .and. depth < clockStackMax ) then
      call time_now ( t )
      clockStack(depth) = t - clockStack(depth)
      call output ( ' used ' )
!     call output ( dble(clockStack(depth) - clockStack(depth+1)), &
!       & format='(g10.3)' )
      call output ( dble(clockStack(depth)), format='(g10.3)' )
    end if
    call output ( '', advance='yes' )
  end subroutine TRACE_END
end module TRACE_M

! $Log$
! Revision 2.10  2001/11/09 23:14:08  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.9  2001/09/13 19:36:50  livesey
! Added optional index arguments
!
! Revision 2.8  2001/05/03 02:13:25  vsnyder
! Trying to print time exclusive of calls isn't working
!
! Revision 2.7  2001/05/01 23:53:40  vsnyder
! Print CPU time exclusive of deeper ones at each end_trace
!
! Revision 2.6  2001/04/25 00:08:26  vsnyder
! Use 'fill' argument of 'output' to get leading zeroes on milliseconds
!
! Revision 2.5  2001/04/24 23:33:56  vsnyder
! Add timing
!
! Revision 2.4  2001/04/24 20:12:47  vsnyder
! Emit blanks before 'Enter' if no tree node supplied
!
! Revision 2.3  2001/04/17 20:56:10  vsnyder
! Emit 'depth' dots in trace_begin even if 'root' is absent
!
! Revision 2.2  2001/03/16 21:01:17  vsnyder
! Make ROOT optional in Trace_Begin
!
! Revision 2.1  2000/10/11 18:33:25  vsnyder
! Move from lib/cf_parser to lib; insert copyright notice
!
! Revision 2.0  2000/09/05 17:41:51  dcuddy
! Change revision to 2.0
!
! Revision 1.1  2000/07/06 01:43:12  vsnyder
! Initial check-in
!
