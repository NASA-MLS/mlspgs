! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Comp_Eta_DoCalc_Sparse_m

  implicit NONE

  private

  public :: Comp_Eta_Docalc_Sparse
  public :: Comp_All_Eta_2D,  Comp_One_Eta_2D
  public :: Comp_All_Eta_QTM, Comp_One_Eta_QTM

  interface Comp_Eta_Docalc_Sparse
    module procedure Comp_All_Eta_2D,  Comp_One_Eta_2D
    module procedure Comp_All_Eta_QTM, Comp_One_Eta_QTM
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine Comp_All_Eta_2D ( Grids_f, Tan_Pt, Z_Path, Eta_Z, &
                             & Phi_Path, Eta_Phi, Eta_zP )

    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use Sparse_m, only: Sparse_t

    type(grids_t), intent(in) :: Grids_f      ! Quantity values
    integer, intent(in) :: Tan_Pt             ! To split Z_Path into two
                                              ! monotone halves
    real(rp), intent(in) :: Z_Path(:)         ! Zetas on the path
    type(sparse_t), intent(inout) :: Eta_Z(:) ! from VMR's zeta to path.
                                              ! InOut so as not to reallocate
                                              ! Eta_Z%Eta
    real(rp), intent(in) :: Phi_Path(:)       ! Phis on the path
    type(sparse_t), intent(out) :: Eta_Phi(:) ! Created here
    class(sparse_t), intent(out) :: Eta_zP(:) ! Created here

    integer :: I

    do i = 1, size(grids_f%mol)
      call comp_one_eta_2d ( grids_f, i, tan_pt, z_path, eta_z(i), phi_path, &
                           & eta_phi(i), eta_zP(i) )
    end do

  end subroutine Comp_All_Eta_2D

  subroutine Comp_All_Eta_QTM ( Grids_f, Z_Path, Eta_Z,  Eta_Path, Eta_zQ )

    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use Sparse_m, only: Sparse_t

    type(grids_t), intent(in) :: Grids_f      ! Quantity values
    real(rp), intent(in) :: Z_Path(:)         ! Zetas on the path
    type(sparse_t), intent(out) :: Eta_Z(:)   ! from VMR's zeta to path.
    type(sparse_t), intent(in) :: Eta_Path(:) ! from Phi to path for
                                              ! all species, because there's
                                              ! only one QTM.
    class(sparse_t), intent(out) :: Eta_zQ(:) ! Created here

    integer :: I

    do i = 1, size(grids_f%mol)
      call comp_one_eta_QTM ( grids_f, i, z_path, eta_z(i), &
                            & eta_path(i), eta_zQ(i) )
    end do

  end subroutine Comp_All_Eta_QTM

  subroutine Comp_One_Eta_2D ( Grids_f, N, Tan_Pt, Z_Path, Eta_Z, &
                             & Phi_Path, Eta_Phi, Eta_zP )

    use intrinsic, only: Lit_Indices
    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use Sparse_m, only: Create_Sparse, Sparse_t
    use Sparse_Eta_m, only: Sparse_Eta_1D, Sparse_Eta_nD

    type(grids_t), intent(in) :: Grids_f      ! Quantity values
    integer, intent(in) :: N                  ! Which quantity
    integer, intent(in) :: Tan_Pt             ! To split Z_Path into two
                                              ! monotone halves
    real(rp), intent(in) :: Z_Path(:)         ! Zetas on the path
    type(sparse_t), intent(out) :: Eta_Z      ! from VMR's Zeta to path.
    real(rp), intent(in) :: Phi_Path(:)       ! Phis on the path, Radians
    type(sparse_t), intent(out) :: Eta_Phi    ! from VMR's Phi to path.
    class(sparse_t), intent(inout) :: Eta_zP  ! Created here

    integer :: P1, P2, Z1, Z2                 ! Boundaries from Grids_f%l_[pz]
    integer :: N_Path

    n_path = size(z_path) ! Assumed == size(p_path)

    p1 = grids_f%l_p(n-1)+1
    p2 = grids_f%l_p(n)
    z1 = grids_f%l_z(n-1)+1
    z2 = grids_f%l_z(n)

    ! Create and compute Eta_Phi
    call sparse_eta_1d ( grids_f%phi_basis(p1:p2), phi_path, eta_phi, &
                       & what=lit_indices(grids_f%mol(n)), resize=.true. )

    ! Create Eta_Z
    call create_sparse ( eta_z, n_path, z2-z1+1, 2*(z2-z1+1), &
                       & what=lit_indices(grids_f%mol(n)) )
    ! Compute the two halves of Eta_Z on either side of the tangent point
    ! separately.  These parts are monotone.  This avoids sorting Z_Path.
    call sparse_eta_1d ( grids_f%zet_basis(z1:z2), z_path, eta_z, &
                       & row1=tan_pt, rowN=1, create=.false. )
    call sparse_eta_1d ( grids_f%zet_basis(z1:z2), z_path, eta_z, &
                       & row1=tan_pt+1, rowN=n_path, create=.false. )

    ! Compute Eta_ZP
    call sparse_eta_nd ( eta_z, eta_phi, eta_zp, &
                       & what=lit_indices(grids_f%mol(n)), resize=.true. )

  end subroutine Comp_One_Eta_2D

  subroutine Comp_One_Eta_QTM ( Grids_f, N, Z_Path, Eta_Z, Eta_Path, Eta_zQ )

    use Sparse_m, only: Sparse_t
    use Sparse_Eta_m, only: Sparse_Eta_1D, Sparse_Eta_nD
    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP

    type(grids_t), intent(in) :: Grids_f      ! Quantity values
    integer, intent(in) :: N                  ! Which quantity
    real(rp), intent(in) :: Z_Path(:)         ! Zetas on the path
    type(sparse_t), intent(out) :: Eta_Z      ! from VMR's zeta to path.
    type(sparse_t), intent(in) :: Eta_Path    ! from Phi to path.
    class(sparse_t), intent(out) :: Eta_zQ    ! Created here

!     integer :: I1, I2                         ! Boundaries from Grids_f%l_*
! 
!     eta_z%n = size(z_path,1)
!     i1 = grids_f%l_z(n-1)+1
!     i2 = grids_f%l_z(n)
!     ! Get the Zeta interpolator list
!     call get_eta_list ( grids_f%zet_basis(i1:i2), z_path, eta_z, sorted=.false. )
!     eta_zQ%n = eta_z%n
!     call get_eta_list ( eta_z%eta(1:eta_z%n), eta_path, eta_zQ )

  end subroutine Comp_One_Eta_QTM

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Comp_Eta_DoCalc_Sparse_m

! $Log$
! Revision 2.4  2017/11/01 19:00:30  vsnyder
! Add tangent point
!
! Revision 2.3  2017/03/11 00:51:18  vsnyder
! Use Grids_F instead of Beta_Group
!
! Revision 2.2  2017/01/14 01:57:09  vsnyder
! Eliminate polymorphic interpolators.  Add arguments to return 1D Etas.
! Assume Z_Path and Phi_Path are not sorted.  Use template%Phi as radians
! because Phi_Path is radians.
!
! Revision 2.1  2016/12/15 02:44:32  vsnyder
! Initial commit
!
