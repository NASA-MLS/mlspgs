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

  subroutine Comp_All_Eta_2D ( Grids_f, Z_Path, Eta_Z, Phi_Path, Eta_Path, &
    & Eta_zP )

    use Indexed_Values_m, only: Value_1D_Lists_t, Value_2D_Lists_t
    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP

    type(grids_t), intent(in) :: Grids_f      ! Quantity values
    real(rp), intent(in) :: Z_Path(:)         ! Zetas on the path
    type(value_1D_lists_t), intent(inout) :: Eta_Z(:) ! from VMR's zeta to path.
                                              ! InOut so as not to reallocate
                                              ! Eta_Z%Eta
    real(rp), intent(in) :: Phi_Path(:)       ! Phis on the path
    type(value_1D_lists_t), intent(inout) :: Eta_Path(:) ! InOut so as not to
                                              ! reallocate Eta_Path%Eta
    class(value_2D_lists_t), intent(inout) :: Eta_zP(:) ! InOut so as not to
                                              ! reallocate Eta_zP%Eta

    integer :: I

    do i = 1, size(grids_f%mol)
      call comp_one_eta_2d ( grids_f, i, z_path, eta_z(i), phi_path, &
                           & eta_path(i), eta_zP(i) )
    end do

  end subroutine Comp_All_Eta_2D

  subroutine Comp_All_Eta_QTM ( Grids_f, Z_Path, Eta_Z,  Eta_Path, Eta_zQ )

    use Indexed_Values_m, only: Value_1D_Lists_t, Value_QTM_1D_List_t, &
      & Value_QTM_2D_Lists_t
    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP

    type(grids_t), intent(in) :: Grids_f      ! Quantity values
    real(rp), intent(in) :: Z_Path(:)         ! Zetas on the path
    type(value_1D_lists_t), intent(inout) :: Eta_Z(:) ! from VMR's zeta to path.
                                              ! InOut so as not to reallocate
                                              ! Eta_Z%Eta
    type(value_QTM_1D_list_t), intent(in) :: Eta_Path(:) ! from Phi to path for
                                              ! all species, because there's
                                              ! only one QTM.
    class(value_QTM_2D_lists_t), intent(inout) :: Eta_zQ(:) ! InOut so as not to
                                              ! reallocate Eta_zQ%Eta

    integer :: I

    do i = 1, size(grids_f%mol)
      call comp_one_eta_QTM ( grids_f, i, z_path, eta_z(i), &
                            & eta_path, eta_zQ(i)%eta )
    end do

  end subroutine Comp_All_Eta_QTM

  subroutine Comp_One_Eta_2D ( Grids_f, N, Z_Path, Eta_Z, Phi_Path, Eta_Path, &
    & Eta_zP )

    use Get_Eta_List_m, only: Get_Eta_List
    use Indexed_Values_m, only: Value_1D_Lists_t, Value_2D_lists_t
    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP

    type(grids_t), intent(in) :: Grids_f      ! Quantity values
    integer, intent(in) :: N                  ! Which quantity
    real(rp), intent(in) :: Z_Path(:)         ! Zetas on the path
    type(value_1D_lists_t), intent(inout) :: Eta_Z ! from VMR's Zeta to path.
                                              ! InOut so as not to reallocate
                                              ! Eta_Z%Eta
    real(rp), intent(in) :: Phi_Path(:)       ! Phis on the path, Radians
    type(value_1D_lists_t), intent(inout) :: Eta_Path ! from VMR's Phi to path.
                                              ! InOut so as not to reallocate
                                              ! Eta_Z%Eta
    class(value_2D_lists_t), intent(inout) :: Eta_zP ! InOut so as not to
                                              ! reallocate Eta_zP%Eta

    integer :: I1, I2                         ! Boundaries from Grids_f%l_*

    eta_z%n = size(z_path,1)
    i1 = grids_f%l_z(n-1)+1
    i2 = grids_f%l_z(n)
    ! Get the Zeta interpolator list
    call get_eta_list ( grids_f%zet_basis(i1:i2), z_path, eta_z, sorted=.false. )
    eta_path%n = size(phi_path,1)
    i1 = grids_f%l_p(n-1)+1
    i2 = grids_f%l_p(n)
    ! Get the Phi interpolator list
    call get_eta_list ( grids_f%phi_basis(i1:i2), phi_path, eta_path, &
      & sorted=.false. )
    eta_zP%n = eta_z%n
    ! Compute the Zeta X Phi interpolator list from the other two
    call get_eta_list ( eta_z%eta(1:eta_z%n), eta_path%eta(1:eta_path%n), &
      & eta_zP%eta(1:eta_zP%n) )
  end subroutine Comp_One_Eta_2D

  subroutine Comp_One_Eta_QTM ( Grids_f, N, Z_Path, Eta_Z, Eta_Path, Eta_zQ )

    use Get_Eta_List_m, only: Get_Eta_List
    use Indexed_Values_m, only: Value_1D_Lists_t, Value_1D_Lists_t, &
      & Value_QTM_1D_List_t, Value_QTM_2D_list_t
    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP

    type(grids_t), intent(in) :: Grids_f      ! Quantity values
    integer, intent(in) :: N                  ! Which quantity
    real(rp), intent(in) :: Z_Path(:)         ! Zetas on the path
    type(value_1D_lists_t), intent(inout) :: Eta_Z ! from VMR's zeta to path.
                                              ! InOut so as not to reallocate
                                              ! Eta_Z%Eta
    type(value_QTM_1D_list_t), intent(in) :: Eta_Path(:) ! from Phi to path.
    class(value_QTM_2D_list_t), intent(inout) :: Eta_zQ(:)

    integer :: I1, I2                         ! Boundaries from Grids_f%l_*

    eta_z%n = size(z_path,1)
    i1 = grids_f%l_z(n-1)+1
    i2 = grids_f%l_z(n)
    ! Get the Zeta interpolator list
    call get_eta_list ( grids_f%zet_basis(i1:i2), z_path, eta_z, sorted=.false. )
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
