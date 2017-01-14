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

  subroutine Comp_All_Eta_2D ( Beta_Group, Z_Path, Eta_Z, Phi_Path, Eta_Path, &
    & Eta_zP )

    use ForwardModelConfig, only: Beta_Group_t
    use Indexed_Values_m, only: Value_1D_Lists_t, Value_2D_Lists_t
    use MLSKinds, only: RP

    type(beta_group_t), intent(in) :: Beta_Group(:)
    real(rp), intent(in) :: Z_Path(:)         ! Zetas on the path
    type(value_1D_lists_t), intent(inout) :: Eta_Z(:) ! from VMR's zeta to path.
                                              ! InOut so as not to reallocate
                                              ! Eta_Z%Eta
    real(rp), intent(in) :: Phi_Path(:)       ! Phis on the path
    type(value_1D_lists_t), intent(inout) :: Eta_Path(:) ! InOut so as not to
                                              ! reallocate Eta_Path%Eta
    type(value_2D_lists_t), intent(inout) :: Eta_zP(:) ! InOut so as not to
                                              ! reallocate Eta_zP%Eta

    integer :: I

    do i = 1, size(beta_group)
      call comp_one_eta_2d ( beta_group(i), z_path, eta_z(i), phi_path, &
                           & eta_path(i), eta_zP(i) )
    end do

  end subroutine Comp_All_Eta_2D

  subroutine Comp_All_Eta_QTM ( Beta_Group, Z_Path, Eta_Z,  Eta_Path, Eta_zQ )

    use ForwardModelConfig, only: Beta_Group_t
    use Indexed_Values_m, only: Value_1D_Lists_t, Value_QTM_1D_Lists_t, &
      & Value_QTM_2D_Lists_t
    use MLSKinds, only: RP

    type(beta_group_t), intent(in) :: Beta_Group(:)
    real(rp), intent(in) :: Z_Path(:)         ! Zetas on the path
    type(value_1D_lists_t), intent(inout) :: Eta_Z(:) ! from VMR's zeta to path.
                                              ! InOut so as not to reallocate
                                              ! Eta_Z%Eta
    type(value_QTM_1D_lists_t), intent(in) :: Eta_Path(:) ! from VMR's Phi to
                                              ! path.
    type(value_QTM_2D_lists_t), intent(inout) :: Eta_zQ(:) ! InOut so as not to
                                              ! reallocate Eta_zQ%Eta

    integer :: I

    do i = 1, size(beta_group)
      call comp_one_eta_QTM ( beta_group(i), z_path, eta_z(i), &
                            & eta_path(i)%eta, eta_zQ(i)%eta )
    end do

  end subroutine Comp_All_Eta_QTM

  subroutine Comp_One_Eta_2D ( Beta_Group, Z_Path, Eta_Z, Phi_Path, Eta_Path, &
    & Eta_zP )

    use Constants, only: Deg2Rad
    use ForwardModelConfig, only: Beta_Group_t
    use Get_Eta_List_m, only: Get_Eta_List
    use Indexed_Values_m, only: Value_1D_Lists_t, Value_2D_lists_t
    use MLSKinds, only: RP

    type(beta_group_t), intent(in) :: Beta_Group
    real(rp), intent(in) :: Z_Path(:)         ! Zetas on the path
    type(value_1D_lists_t), intent(inout) :: Eta_Z ! from VMR's Zeta to path.
                                              ! InOut so as not to reallocate
                                              ! Eta_Z%Eta
    real(rp), intent(in) :: Phi_Path(:)       ! Phis on the path, Radians
    type(value_1D_lists_t), intent(inout) :: Eta_Path ! from VMR's Phi to path.
                                              ! InOut so as not to reallocate
                                              ! Eta_Z%Eta
    type(value_2D_lists_t), intent(inout) :: Eta_zP ! InOut so as not to
                                              ! reallocate Eta_zP%Eta

    eta_z%n = size(z_path,1)
    call get_eta_list ( beta_group%qty%qty%template%surfs(:,1), z_path, &
      & eta_z, sorted=.false. )
    eta_path%n = size(phi_path,1)
    call get_eta_list ( deg2rad*beta_group%qty%qty%template%phi(1,:), phi_path, &
      & eta_path, sorted=.false. )
    eta_zP%n = eta_z%n
    call get_eta_list ( eta_z%eta(1:eta_z%n), eta_path%eta(1:eta_path%n), &
      & eta_zP%eta(1:eta_zP%n) )

  end subroutine Comp_One_Eta_2D

  subroutine Comp_One_Eta_QTM ( Beta_Group, Z_Path, Eta_Z, Eta_Path, &
    & Eta_zQ )

    use ForwardModelConfig, only: Beta_Group_t
    use Get_Eta_List_m, only: Get_Eta_List
    use Indexed_Values_m, only: Value_1D_Lists_t, Value_1D_Lists_t, &
      & Value_QTM_1D_List_t, Value_QTM_2D_list_t
    use MLSKinds, only: RP

    type(beta_group_t), intent(in) :: Beta_Group
    real(rp), intent(in) :: Z_Path(:)         ! Zetas on the path
    type(value_1D_lists_t), intent(inout) :: Eta_Z ! from VMR's zeta to path.
                                              ! InOut so as not to reallocate
                                              ! Eta_Z%Eta
    type(value_QTM_1D_list_t), intent(in) :: Eta_Path(:) ! from VMR's Phi to path.
    type(value_QTM_2D_list_t), intent(inout) :: Eta_zQ(:)

    eta_z%n = size(z_path,1)
    call get_eta_list ( beta_group%qty%qty%template%surfs(:,1), z_path, &
        & eta_z, sorted=.true. )
    eta_zQ%n = eta_z%n
    call get_eta_list ( eta_z%eta(1:eta_z%n), eta_path, eta_zQ )

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
! Revision 2.2  2017/01/14 01:57:09  vsnyder
! Eliminate polymorphic interpolators.  Add arguments to return 1D Etas.
! Assume Z_Path and Phi_Path are not sorted.  Use template%Phi as radians
! because Phi_Path is radians.
!
! Revision 2.1  2016/12/15 02:44:32  vsnyder
! Initial commit
!
