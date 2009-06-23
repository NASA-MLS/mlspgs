
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module HFTI_M ! From Lawson and Hanson, converted to Fortran 90
!=============================================================================

  implicit NONE
  private

  public :: HFTI, SHFTIV, SHFTIM, DHFTIV, DHFTIM

  interface HFTI
    module procedure SHFTIV, SHFTIM, DHFTIV, DHFTIM
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine SHFTIM ( A, B, TAU, KRANK, RNORM, IP )

    !{ Solve $A X \simeq B$.

    use ERMSG_M, only: ERMOR, ERMSG, ERM1, ERV1
    use HT_M, only: HTGEN
    use NRM2_M, only: NRM2

    integer, parameter :: RK = kind(0.0)
    real(rk), intent(inout) :: A(:,:) ! A on input, input for COV2 on output
    real(rk), intent(inout) :: B(:,:) ! B on input, X on output
    real(rk), optional, intent(in) :: TAU       ! Tolerance
    integer, optional, intent(out) :: KRANK     ! Rank of A
    real(rk), optional, intent(out) :: RNORM(:) ! sqrt(sum(residuals**2))
    integer, optional, intent(out) :: IP(:)     ! Column interchanges
    include 'hftim.f9h'

  end subroutine SHFTIM

  subroutine SHFTIV ( A, B, TAU, KRANK, RNORM, IP )

    !{ Solve $A x \simeq b$.

    use ERMSG_M, only: ERMOR, ERMSG, ERM1, ERV1
    use HT_M, only: HTGEN
    use NRM2_M, only: NRM2

    integer, parameter :: RK = kind(0.0)
    real(rk), intent(inout) :: A(:,:) ! A on input, input for COV2 on output
    real(rk), intent(inout) :: B(:)   ! B on input, X on output
    real(rk), optional, intent(in) :: TAU       ! Tolerance
    integer, optional, intent(out) :: KRANK     ! Rank of A
    real(rk), optional, intent(out) :: RNORM    ! sqrt(sum(residuals**2))
    integer, optional, intent(out) :: IP(:)     ! Column interchanges
    include 'hftiv.f9h'

  end subroutine SHFTIV

  subroutine DHFTIM ( A, B, TAU, KRANK, RNORM, IP )

    !{ Solve $A X \simeq B$.

    use ERMSG_M, only: ERMOR, ERMSG, ERM1, ERV1
    use HT_M, only: HTGEN
    use NRM2_M, only: NRM2

    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(inout) :: A(:,:) ! A on input, input for COV2 on output
    real(rk), intent(inout) :: B(:,:) ! B on input, X on output
    real(rk), optional, intent(in) :: TAU       ! Tolerance
    integer, optional, intent(out) :: KRANK     ! Rank of A
    real(rk), optional, intent(out) :: RNORM(:) ! sqrt(sum(residuals**2))
    integer, optional, intent(out) :: IP(:)     ! Column interchanges
    include 'hftim.f9h'

  end subroutine DHFTIM

  subroutine DHFTIV ( A, B, TAU, KRANK, RNORM, IP )

    !{ Solve $A x \simeq b$.

    use ERMSG_M, only: ERMOR, ERMSG, ERM1, ERV1
    use HT_M, only: HTGEN
    use NRM2_M, only: NRM2

    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(inout) :: A(:,:) ! A on input, input for COV2 on output
    real(rk), intent(inout) :: B(:)   ! B on input, X on output
    real(rk), optional, intent(in) :: TAU       ! Tolerance
    integer, optional, intent(out) :: KRANK     ! Rank of A
    real(rk), optional, intent(out) :: RNORM    ! sqrt(sum(residuals**2))
    integer, optional, intent(out) :: IP(:)     ! Column interchanges
    include 'hftiv.f9h'

  end subroutine DHFTIV

!=======================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module HFTI_M

! $Log$
! Revision 2.3  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.2  2006/03/23 03:04:33  vsnyder
! OOPS, need to get ER[VM]1 from ermsg_m
!
! Revision 2.1  2006/03/22 02:07:04  vsnyder
! Initial commit
!
