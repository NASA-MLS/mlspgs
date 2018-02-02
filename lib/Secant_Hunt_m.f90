! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Secant_Hunt_m

  implicit NONE
  private
  public :: Secant_Hunt

  interface Secant_Hunt
    module procedure Secant_Hunt_Int, Secant_Hunt_Real
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Secant_Hunt_Int ( Array, Target, Loc, Probes, First )

    ! Look for Target in Array(1:N) using a secant method, which has log log N
    ! complexity.  If Target < Array(1), Loc = 0.  If Target >= Array(N), Loc =
    ! N.  Else Loc = I such that Array(i) <= Loc <= Array(i+1) with equality in
    ! the upper case only if Array(i) == Array(i+1).  If Array(i) == Array(i+1)
    ! for several values of I, and First is absent or false, any of those values
    ! of I could be reported.  If Array(i) == Array(i+1) for several values of
    ! I, and First is present and true, the first of these is reported.

    ! Assumes Array(:) is ordered.  Will fail silently and badly otherwise.

    integer, intent(in) :: Array(:)          ! The array to search
    integer, intent(in) :: Target            ! The value to find
    integer, intent(out) :: Loc              ! Where Target was found in Array(:)
    integer, intent(out), optional :: Probes ! How many places were examined?
    logical, intent(in), optional :: First   ! Find the first of equals

    include 'Secant_Hunt.f9h'

  end subroutine Secant_Hunt_Int

  subroutine Secant_Hunt_Real ( Array, Target, Loc, Probes, First )

    ! Look for Target in Array(1:N) using a secant method, which has log log N
    ! complexity.  If Target < Array(1), Loc = 0.  If Target >= Array(N), Loc =
    ! N.  Else Loc = I such that Array(i) <= Loc <= Array(i+1) with equality in
    ! the upper case only if Array(i) == Array(i+1).  If Array(i) == Array(i+1)
    ! for several values of I, and First is absent or false, any of those values
    ! of I could be reported.  If Array(i) == Array(i+1) for several values of
    ! I, and First is present and true, the first of these is reported.

    ! Assumes Array(:) is ordered.  Will fail silently and badly otherwise.

    real, intent(in) :: Array(:)             ! The array to search
    real, intent(in) :: Target               ! The value to find
    integer, intent(out) :: Loc              ! Where Target was found in Array(:)
    integer, intent(out), optional :: Probes ! How many places were examined?
    logical, intent(in), optional :: First   ! Find the first of equals

    include 'Secant_Hunt.f9h'

  end subroutine Secant_Hunt_Real

  subroutine Secant_Hunt_Dble ( Array, Target, Loc, Probes, First )

    ! Look for Target in Array(1:N) using a secant method, which has log log N
    ! complexity.  If Target < Array(1), Loc = 0.  If Target >= Array(N), Loc =
    ! N.  Else Loc = I such that Array(i) <= Loc <= Array(i+1) with equality in
    ! the upper case only if Array(i) == Array(i+1).  If Array(i) == Array(i+1)
    ! for several values of I, and First is absent or false, any of those values
    ! of I could be reported.  If Array(i) == Array(i+1) for several values of
    ! I, and First is present and true, the first of these is reported.

    ! Assumes Array(:) is ordered.  Will fail silently and badly otherwise.

    double precision, intent(in) :: Array(:) ! The array to search
    double precision, intent(in) :: Target   ! The value to find
    integer, intent(out) :: Loc              ! Where Target was found in Array(:)
    logical, intent(in), optional :: First   ! Find the first of equals
    integer, intent(out), optional :: Probes ! How many places were examined?

    include 'Secant_Hunt.f9h'

  end subroutine Secant_Hunt_Dble

!-----------------------------------------------------------------------
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Secant_Hunt_m

! $Log$
! Revision 2.1  2018/02/02 03:32:36  vsnyder
! Initial commit
!
