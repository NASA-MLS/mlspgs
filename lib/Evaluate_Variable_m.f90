! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Evaluate_Variable_m

  implicit NONE
  private
  public :: Evaluate_Variable

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================

  subroutine Evaluate_Variable ( Root )

! Evaluate a variable defined by <name> := <expr>.
! Put its value in the Values component of its declaration.

    use Declaration_Table, only: Allocate_Test, Deallocate_Test, Range, &
      & Redeclare, Str_Range, Value_t, Variable
    use Expr_m, only: Expr
    use Toggles, only: Gen, Toggle
    use Trace_m, only: Trace_Begin, Trace_End
    use Tree, only: Nsons, Subtree, Sub_Rosa
    integer, intent(in) :: Root  ! The := tree node of the definition

    integer :: I
    integer :: Me = -1           ! String index for trace cacheing
    integer :: N                 ! For allocation
    integer :: Type
    integer :: Units(2)          ! From expr
    double precision :: Value(2) ! From expr
    type(value_t), allocatable :: Values1(:), Values2(:), Values3(:)

    call trace_begin ( me, 'Evaluate_Variable', root, cond=toggle(gen) )
    do i = 2, nsons(root)
      call expr ( subtree(i,root), units, value, type, values=values2 )
      if ( allocated(values1) ) then
        call allocate_test ( values3, size(values1)+size(values2), &
          & 'Values3', moduleName )
        values3 = [ values1, values2 ]
        call deallocate_test ( values1, 'Values1', moduleName )
        call deallocate_test ( values2, 'Values2', moduleName )
        call move_alloc ( values3, values1 )
      else
        call move_alloc ( values2, values1 )
      end if
    end do
    call redeclare ( sub_rosa(subtree(1,root)), value(1), variable, &
      & values=values1 ) ! Declares it if it wasn't already declared
    call trace_end ( 'Evaluate_Variable', cond=toggle(gen) )

  end subroutine Evaluate_Variable

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Evaluate_Variable_m

! $Log$
! Revision 2.1  2013/10/09 01:06:36  vsnyder
! Initial commit
!
