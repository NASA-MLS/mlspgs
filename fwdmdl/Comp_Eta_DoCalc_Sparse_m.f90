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

  subroutine Comp_All_Eta_2D ( Beta_Group, Z_Path, Eta_Path, Eta_Z, Eta_zP )

    use ForwardModelConfig, only: Beta_Group_t
    use Get_Eta_List_m, only: Eta_Lists_t
    use Indexed_Values_m, only: Value_1D_List_t
    use MLSKinds, only: RP

    type(beta_group_t), intent(in) :: Beta_Group(:)
    real(rp), intent(in) :: Z_Path(:)         ! Zetas on the path
    type(value_1D_list_t), intent(in) :: Eta_Path(:)
    type(eta_lists_t), intent(inout) :: Eta_Z(:) ! from VMR's zeta to path.
                                              ! InOut so as not to reallocate
                                              ! Eta_Z%Eta
    type(eta_lists_t), intent(inout) :: Eta_zP(:) ! InOut so as not to
                                              ! reallocate Eta_zP%Eta

    integer :: I

    do i = 1, size(beta_group)
      call comp_one_eta_2d ( beta_group(i), z_path, eta_path, &
                           & eta_z(i), eta_zP(i) )
    end do

  end subroutine Comp_All_Eta_2D

  subroutine Comp_All_Eta_QTM ( Beta_Group, Z_Path, Eta_Path, Eta_Z, Eta_zQ )

    use ForwardModelConfig, only: Beta_Group_t
    use Get_Eta_List_m, only: Eta_Lists_t
    use Indexed_Values_m, only: Value_QTM_1D_List_t
    use MLSKinds, only: RP

    type(beta_group_t), intent(in) :: Beta_Group(:)
    real(rp), intent(in) :: Z_Path(:)         ! Zetas on the path
    type(value_QTM_1D_list_t), intent(in) :: Eta_Path(:)
    type(eta_lists_t), intent(inout) :: Eta_Z(:) ! from VMR's zeta to path.
                                              ! InOut so as not to reallocate
                                              ! Eta_Z%Eta
    type(eta_lists_t), intent(inout) :: Eta_zQ(:) ! InOut so as not to
                                              ! reallocate Eta_zQ%Eta

    integer :: I

    do i = 1, size(beta_group)
      call comp_one_eta_QTM ( beta_group(i), z_path, eta_path, &
                            & eta_z(i), eta_zQ(i) )
    end do

  end subroutine Comp_All_Eta_QTM

  subroutine Comp_One_Eta_2D ( Beta_Group, Z_Path, Eta_Path, Eta_Z, &
    & Eta_zP )

    use ForwardModelConfig, only: Beta_Group_t
    use Get_Eta_List_m, only: Eta_Lists_t, Eta_List_1D
    use Indexed_Values_m, only: Value_1D_List_t, Value_2D_list_t
    use MLSKinds, only: RP

    type(beta_group_t), intent(in) :: Beta_Group
    real(rp), intent(in) :: Z_Path(:)         ! Zetas on the path
    type(value_1D_list_t), intent(in) :: Eta_Path(:)
    type(eta_lists_t), intent(inout) :: Eta_Z ! from VMR's zeta to path.
                                              ! InOut so as not to reallocate
                                              ! Eta_Z%Eta
    type(eta_lists_t), intent(inout) :: Eta_zP ! InOut so as not to
                                              ! reallocate Eta_zP%Eta

    call eta_list_1d ( beta_group%qty%qty%template%surfs(:,1), z_path, &
        & eta_z, sorted=.true. )
    ! Select Type is needed because eta_z%eta and eta_zP%eta are polymorphic,
    ! and can't be type bound because they're not scalars.
    select type ( the_eta_z => eta_z%eta )
    type is ( value_1D_list_t )
      select type ( the_eta_zP => eta_zP%eta )
      type is ( value_2D_list_t )
        call eta_list_1d ( the_eta_z, eta_path, the_eta_zP )
      end select
    end select

  end subroutine Comp_One_Eta_2D

  subroutine Comp_One_Eta_QTM ( Beta_Group, Z_Path, Eta_Path, Eta_Z, &
    & Eta_zQ )

    use ForwardModelConfig, only: Beta_Group_t
    use Get_Eta_List_m, only: Eta_Lists_t, Eta_List_1D
    use Indexed_Values_m, only: Value_1D_List_t, &
      & Value_QTM_1D_List_t, value_QTM_2D_list_t
    use MLSKinds, only: RP

    type(beta_group_t), intent(in) :: Beta_Group
    real(rp), intent(in) :: Z_Path(:)         ! Zetas on the path
    type(value_QTM_1D_list_t), intent(in) :: Eta_Path(:)
    type(eta_lists_t), intent(inout) :: Eta_Z ! from VMR's zeta to path.
                                              ! InOut so as not to reallocate
                                              ! Eta_Z%Eta
    type(eta_lists_t), intent(inout) :: Eta_zQ ! InOut so as not to
                                              ! reallocate Eta_zQ%Eta

    call eta_list_1d ( beta_group%qty%qty%template%surfs(:,1), z_path, &
        & eta_z, sorted=.true. )
    ! Select Type is needed because eta_z%eta and eta_zQ%eta are polymorphic,
    ! and can't be type bound because they're not scalars.
    select type ( the_eta_z => eta_z%eta )
    type is ( value_1D_list_t )
      select type ( the_eta_zQ => eta_zQ%eta )
      type is ( value_QTM_2D_list_t )
        call eta_list_1d ( the_eta_z, eta_path, the_eta_zQ )
      end select
    end select

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
! Revision 2.1  2016/12/15 02:44:32  vsnyder
! Initial commit
!
