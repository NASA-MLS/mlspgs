! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module TRACE_M

  use LEXER_CORE, only: PRINT_SOURCE
  use OUTPUT_M, only: OUTPUT
  use TREE, only: DUMP_TREE_NODE, SOURCE_REF

  integer, public, save :: DEPTH   ! Depth in tree.  Used for trace printing.

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================
! --------------------------------------------------  TRACE_BEGIN  -----
  subroutine TRACE_BEGIN ( NAME, ROOT )
  ! Print "ENTER NAME with ROOT = <node_id(root)>" with DEPTH dots in
  ! front.  Increment DEPTH.
    character(len=*), intent(in) :: NAME
    integer, intent(in), optional :: ROOT
    integer :: I              ! Loop inductor
    integer :: Values(8)      ! For Date_and_time
    if ( present(root) ) then
      call output ( root, 4 ); call output ( ': ' )
    else
      call output ( '      ' )
    end if
    do i = 1, depth
      call output ( '.' )
    end do
    call output ( 'Enter ' ); call output ( name )
    call date_and_time ( values=values )
    call output ( ' at ' )
    call output ( values(5) ); call output ( ':' ) ! The hour
    call output ( values(6) ); call output ( ':' ) ! The minute
    call output ( values(7) ); call output ( '.' ) ! The second
    call output ( values(8), places=3 )            ! The milliseconds
    if ( present(root) ) then
      call output ( ' with ' );
      call dump_tree_node ( root, 0 )
      call output ( ' at ' )
      call print_source ( source_ref(root), advance='yes' )
    else
      call output ( '', advance='yes' )
    end if
    depth = depth + 1
  end subroutine TRACE_BEGIN
! --------------------------------------------------    TRACE_END  -----
  subroutine TRACE_END ( NAME )
  ! Decrement DEPTH.  Print "EXIT NAME with DEPTH dots in front.
    character(len=*), intent(in) :: NAME
    integer :: I              ! Loop inductor
    integer :: Values(8)      ! For Date_and_time
    depth = depth - 1
    call output ( '      ' )
    do i = 1, depth
      call output ( '.' )
    end do
    call output ( 'Exit ' ); call output ( name )
    call date_and_time ( values=values )
    call output ( ' at ' )
    call output ( values(5) ); call output ( ':' )     ! The hour
    call output ( values(6) ); call output ( ':' )     ! The minute
    call output ( values(7) ); call output ( '.' )     ! The second
    call output ( values(8), advance='yes', places=3 ) ! The milliseconds
  end subroutine TRACE_END
end module TRACE_M

! $Log$
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
