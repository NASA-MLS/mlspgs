! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module LinearizedForwardModel_m

  ! Approximate a forward model with a first-order Taylor's series.

  use L2PC_m, only: L2PC_T
  use MatrixModule_1, only: CopyMatrix, DestroyMatrix, Matrix_T, &
    & MultiplyMatrixVectorNoT
  use VectorsModule, only: assignment(=), CopyVector, DestroyVectorInfo, &
    & operator(-), Vector_T

  implicit none
  private
  public :: LinearizedForwardModel

  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), parameter :: IdParm = & 
       "$Id$"
  character(len=len(idParm)) :: Id = idParm
  character(len=*), parameter :: ModuleName = &
    & "$RCSfile$"
  !-----------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  ! -------------------------------------  LinearizedForwardModel  -----
  subroutine LinearizedForwardModel ( L2PC, State, Radiance, Jacobian )
    type(l2PC_T), intent(in) :: L2PC
    type(vector_T), intent(in) :: State
    type(vector_T), intent(inout) :: Radiance
    type(matrix_T), intent(inout) :: Jacobian

    type(matrix_T) :: KP                ! Intermediate Jacobian
    type(vector_T) :: XP, YP            ! Intermediate State and Radiance

    ! XP = Appropriate subset of State
    ! XP = XP - x*
    xp = xp - l2pc%xStar
    ! YP = y* + K* \times ( XP - x*)
    call copyVector ( yp, l2pc%yStar, clone=.true. )
    call MultiplyMatrixVectorNoT ( l2pc%kStar, xp, yp, update=.true. )
    ! KP = K*
    call copyMatrix ( kp, l2pc%kStar )

    ! Interpolate YP to the ptans for X and the correct MAF.
    ! Fill Jacobian for these ptans with free output from the
    !  interpolator.
    ! Fill the rest of Jacobian by interpolating KP to the ptans for
    !  X and the correct MAF.

    call copyVector ( radiance, yp )

    call destroyMatrix ( kp )
    call destroyVectorInfo ( xp )
    call destroyVectorInfo ( yp )
  end subroutine LinearizedForwardModel
end module LinearizedForwardModel_m

! $Log$
