! Copyright 2015, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module Weight_1D_m
!=============================================================================

  ! Compute 1-dimensional linear interpolation weights.

  use Pure_Hunt_m, only: PureHunt

  implicit NONE
  private

  public :: Weight_1D

  interface Weight_1D
    module procedure Weight_1D_D, Weight_1D_S
  end interface
  
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine Weight_1D_D ( Heights, Height, Jlo, Jhi, Weights )
    ! Interpolate Height in Heights to calculate Weight
    integer, parameter :: RK = kind(1.0d0)
    real(rk), intent(in) :: Heights(:)  ! Array in which to interpolate
    real(rk), intent(in) :: Height      ! Value to which to interpolate
    integer, intent(inout) :: Jlo, Jhi  ! Which array elements used
    real(rk), intent(out) :: Weights(2) ! Interpolation coefficients
    call purehunt ( height, heights, size(heights), jlo, jhi )
    if ( jlo == jhi ) then ! below the bottom or above the top
      weights = 0.5
    else
      weights(1) = ( heights(jhi) - height ) / &
                 & ( heights(jhi) - heights(jlo) )
      weights(2) = 1 - weights(1)
    end if
  end subroutine Weight_1D_D

  subroutine Weight_1D_S ( Heights, Height, Jlo, Jhi, Weights )
    ! Interpolate Height in Heights to calculate Weight
    integer, parameter :: RK = kind(1.0e0)
    real(rk), intent(in) :: Heights(:)  ! Array in which to interpolate
    real(rk), intent(in) :: Height      ! Value to which to interpolate
    integer, intent(inout) :: Jlo, Jhi  ! Which array elements used
    real(rk), intent(out) :: Weights(2) ! Interpolation coefficients
    call purehunt ( height, heights, size(heights), jlo, jhi )
    if ( jlo == jhi ) then ! below the bottom or above the top
      weights = 0.5
    else
      weights(1) = ( heights(jhi) - height ) / &
                 & ( heights(jhi) - heights(jlo) )
      weights(2) = 1 - weights(1)
    end if
  end subroutine Weight_1D_S

!=============================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Weight_1D_m

! $Log$
! Revision 2.2  2016/10/24 22:11:05  vsnyder
! Handle degenerate cases
!
! Revision 2.1  2016/09/09 00:45:20  vsnyder
! Initial commit
!
