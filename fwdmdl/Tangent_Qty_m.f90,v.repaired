! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Tangent_Qty_m          ! Get qty at the tangent point

  ! Interpolate a state-vector quantity, e.g., temperature or water, to MIF
  ! tangent points.

  implicit none

  private

  public :: Tangent_Qty

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Tangent_Qty ( FmConf, UsingQTM, MAF, PTan, PhiTan, Qty, &
                         & WindowStart, WindowFinish, &
                         & TanQty, Eta_ZXP, Eta_P_Out, Eta_Z_Out, LogLin )

    use ForwardModelConfig, only: ForwardModelConfig_T
    use MLSKinds, only: RP
    use Sparse_Eta_m, only: Sparse_Eta_t
    use VectorsModule, only: VectorValue_T

    type (forwardModelConfig_T), intent(inout) :: FmConf ! Configuration options
    logical, intent(in) :: UsingQTM
    integer, intent(in) :: MAF
    type (vectorValue_T), intent(in) :: PTan    ! Tangent pressures
    type (vectorValue_T), intent(in) :: PhiTan  ! Tangent horizontal coordinates
    type (vectorValue_T), intent(in) :: Qty     ! Minor-frame quantity
    integer, intent(in) :: WindowStart, WindowFinish ! Extent of Qty%Values
    real(rp), intent(out) :: TanQty(:)          ! Interpolated result
    type (sparse_eta_t), intent(out) :: Eta_ZXP ! Zeta X Horizontal interpolator
    type(sparse_eta_t), intent(out), optional, target :: Eta_P_Out, Eta_Z_Out
    logical, intent(in), optional :: LogLin     ! Use log-linear interpolation

    type (sparse_eta_t), target :: Eta_P_Local, Eta_Z_Local
    type (sparse_eta_t), pointer :: Eta_P
    type (sparse_eta_t), pointer :: Eta_Z

    ! Is Eta_P to be output?
    eta_p => eta_p_local
    if ( present(eta_p_out) ) eta_p => eta_p_out
    ! Is Eta_Z to be output?
    eta_z => eta_z_local
    if ( present(eta_z_out) ) eta_z => eta_z_out

    ! Get horizontal interpolating coefficients
    if ( usingQTM ) then
      ! Use phitan%template%geodLat(:,fmConf%MAF) and
      ! phitan%template%lon(:,fmConf%MAF) to find QTM
      ! facets and compute interpolation coefficients.
      ! Eta_p%Eta_QTM allocates components -- no need to call eta_p%create
      call eta_p%eta_QTM ( qty%template%the_Hgrid%QTM_Tree, &
                         & phitan%template%lon(:,MAF), &
                         & phitan%template%geodLat(:,MAF), &
                         & what=qty%template%molecule, empty=.true. )
    else
      ! Eta_p%Eta_1D allocates components -- no need to call eta_p%create If
      ! this is changed to pass WindowStart and WindowFinish to the Basis1 and
      ! BasisN arguments of Eta_P%Eta_1D, DO NOT USE WindowStart:WindowFinish
      ! as a section subscript in the second dimension of the "R" argument to
      ! the sparse dot-product routines.
      call eta_p%eta_1d ( qty%template%phi(1,windowStart:windowFinish), &
                        & phitan%values(:,maf), what=qty%template%molecule, &
                        & empty=.true., sorted=.false. )
    end if

    ! Compute Zeta and Zeta X Phi interpolating coefficients
    call eta_z%eta_1d ( qty%template%surfs(:,1), ptan%values(:,maf), &
                      & what=qty%template%molecule, empty=.true., sorted=.false. )
    call eta_zxp%eta_nd ( eta_z, eta_p )

    ! Compute tanQty
    ! We might use logarithmic interpolation here:
    ! tanQty = exp( eta_zxp .dot. log(max(qty%values,1.0e-9_rp)) )
    call eta_zxp%sparse_dot_vec ( &
      & qty%values(:,windowStart:windowFinish), tanQty, explog = logLin )

  end subroutine Tangent_Qty

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Tangent_Qty_m

! $Log$
! Revision 2.2  2020/05/07 21:24:17  vsnyder
! Add comment about section subscript when interpolating
!
! Revision 2.1  2018/12/13 01:13:16  vsnyder
! Initial commit
!
