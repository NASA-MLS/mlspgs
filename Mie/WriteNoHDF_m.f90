! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module WriteHDF_m

  implicit NONE

  private

  public :: WriteHDF

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine WriteHDF ( File, R_max, R_min, N_cut, IWC_s, T_s, Theta_s, F_s, &
    & Beta, Eest, nFunc, MaxOrd, P, E_P, nFuncP, MaxOrdP, &
    & WantBeta, WantIWC, WantP, &
    & dBeta_dIWC, E_dBeta_dIWC, dBeta_dT, E_dBeta_dT, &
    & dP_dIWC, E_dP_dIWC, dP_dT, E_dP_dT )

    ! A stub that doesn't write the tables as HDF5 datasets
    use MLSKinds, only: R8

    character(len=*), intent(in) :: File
    real(r8), intent(in) :: R_Max
    real(r8), intent(in) :: R_Min
    integer, intent(in) :: N_Cut
    real(r8), intent(in) :: IWC_s(:)
    real(r8), intent(in) :: T_s(:)
    real(r8), intent(in) :: Theta_s(:)
    real(r8), intent(in) :: F_s(:)
    real(r8), intent(in) :: Beta(:,:,:,:)
    real(r8), intent(in) :: Eest(:,:,:,:)
    integer, intent(in) :: nFunc(:,:,:,:,:)
    integer, intent(in) :: MaxOrd(:,:,:,:)
    real(r8), intent(in) :: P(:,:,:,:)
    real(r8), intent(in) :: E_P(:,:,:,:)
    integer, intent(in) :: nFuncP(:,:,:,:,:,:)
    integer, intent(in) :: MaxOrdP(:,:,:,:,:)
    logical, intent(in) :: WantBeta, WantIWC, WantP
    real(r8), intent(in), optional :: dBeta_dIWC(:,:,:,:)
    real(r8), intent(in), optional :: E_dBeta_dIWC(:,:,:,:)
    real(r8), intent(in), optional :: dBeta_dT(:,:,:,:)
    real(r8), intent(in), optional :: E_dBeta_dT(:,:,:,:)
    real(r8), intent(in), optional :: dP_dIWC(:,:,:,:)
    real(r8), intent(in), optional :: E_dP_dIWC(:,:,:,:)
    real(r8), intent(in), optional :: dP_dT(:,:,:,:)
    real(r8), intent(in), optional :: E_dP_dT(:,:,:,:)

  end subroutine WriteHDF

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, id
  end function not_used_here

end module WriteHDF_m
