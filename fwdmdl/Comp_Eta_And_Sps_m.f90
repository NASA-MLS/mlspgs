! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Comp_Eta_And_Sps_m

  implicit NONE

  private

  public :: Comp_Sps

  interface Comp_Sps
    module procedure Comp_Eta_And_One_Sps, Comp_Eta_And_Several_Sps
    module procedure Comp_Only_Sps
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine Comp_Eta_And_One_Sps ( Grids, Tan_Pt, Z_Path, Phi_Path, &
                                  & Eta, Path_Value )

    ! Interpolate one quantity, which depends only upon Zeta and Phi, and for
    ! which we need to remember any of the Zeta X Phi interpolation coefficients.

    use Comp_Eta_DoCalc_Sparse_m, only: Comp_Eta
    use Comp_Sps_Path_Sparse_m, only: Comp_Sps_Path_Sparse
    use GLNP, only: NGP1
    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use Sparse_Eta_m, only: Sparse_Eta_t

    type(grids_t), intent(in) :: Grids        ! Quantity values
    integer, intent(in) :: Tan_Pt             ! To split Z_Path into two
                                              ! monotone halves
    real(rp), intent(in) :: Z_Path(:)         ! Zetas on the path
    real(rp), intent(in) :: Phi_Path(:)       ! Phis on the path
    type(sparse_eta_t), intent(out) :: Eta    ! Eta Zeta X Phi interpolator
    real(rp), intent(out) :: Path_Value(:,:)  ! Grids%v interpolated to Z X P

    call comp_eta ( grids, tan_pt, z_path, phi_path, eta, NGP1 )
    call comp_sps_path_sparse ( grids, eta, path_value )

  end subroutine Comp_Eta_And_One_Sps

  subroutine Comp_Eta_And_Several_Sps ( Grids, Tan_Pt, Z_Path, Phi_Path, &
                                      & Eta, Path_Value )

    ! Interpolate one quantity, which depends only upon Zeta and Phi, and for
    ! which we need to remember any of the Zeta X Phi interpolation coefficients.

    use Comp_Eta_DoCalc_Sparse_m, only: Comp_Eta
    use Comp_Sps_Path_Sparse_m, only: Comp_Sps_Path_Sparse
    use GLNP, only: NGP1
    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use Sparse_Eta_m, only: Sparse_Eta_t

    type(grids_t), intent(in) :: Grids        ! Quantity values
    integer, intent(in) :: Tan_Pt             ! To split Z_Path into two
                                              ! monotone halves
    real(rp), intent(in) :: Z_Path(:)         ! Zetas on the path
    real(rp), intent(in) :: Phi_Path(:)       ! Phis on the path
    type(sparse_eta_t), intent(out) :: Eta(:) ! Eta Zeta X Phi interpolator
    real(rp), intent(out) :: Path_Value(:,:)  ! Grids%v interpolated to Z X P

    call comp_eta ( grids, tan_pt, z_path, phi_path, eta, NGP1 )
    call comp_sps_path_sparse ( grids, eta, path_value )

  end subroutine Comp_Eta_And_Several_Sps

  subroutine Comp_Only_Sps ( Grids, Tan_Pt, Z_Path, Phi_Path, Path_Value )

    ! Interpolate one quantity, which depends only upon Zeta and Phi, and for
    ! which we do not need to remember any of the interpolation coefficients.

    use Comp_Eta_DoCalc_Sparse_m, only: Comp_Eta
    use Comp_Sps_Path_Sparse_m, only: Comp_Sps_Path_Sparse
    use GLNP, only: NGP1
    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use Sparse_Eta_m, only: Sparse_Eta_t

    type(grids_t), intent(in) :: Grids        ! Quantity values
    integer, intent(in) :: Tan_Pt             ! To split Z_Path into two
                                              ! monotone halves
    real(rp), intent(in) :: Z_Path(:)         ! Zetas on the path
    real(rp), intent(in) :: Phi_Path(:)       ! Phis on the path
    real(rp), intent(out) :: Path_Value(:,:)  ! Grids%v interpolated to Z X P

    type(sparse_eta_t) :: Eta                 ! Eta Zeta X Phi interpolator

    call comp_eta ( grids, tan_pt, z_path, phi_path, eta, NGP1 )
    call comp_sps_path_sparse ( grids, eta, path_value )

  end subroutine Comp_Only_Sps

! ------------------------------------------------  not_used_here  -----
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Comp_Eta_And_Sps_m

! $Log$
! Revision 2.1  2018/09/12 20:43:41  vsnyder
! Initial commit
!
