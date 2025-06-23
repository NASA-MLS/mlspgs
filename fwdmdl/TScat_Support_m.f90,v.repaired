! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module TScat_Support_m

  ! Provide types and procedures for use during calculation and use of TScat
  ! tables by the Full Forward Model.

  implicit NONE
  private
  public :: Find_Scattering_Point
  public :: Get_dB_df, Get_dB_dT, Get_TScat, Get_TScat_Setup
  public :: Get_TScat_Teardown, Interpolate_P_to_theta_e, Mie_Freq_Index
  public :: TScat_Detail_Heading, TScat_Gen_Setup


  integer, public :: Print_TScat_Deriv  ! For debugging, Temp deriv, from -Sdtsct
  integer, public :: Print_TScat_Detail ! For debugging, from -Spsct

  character(len=*), parameter :: TScat_Detail_Heading = &
    & "Scat_Phi   Scat_Ht   Scat_Zeta     Xi       D2     Tan_Phi    Tan_Ht  Begin Phi   End Phi Which Rev Fwd"

!------------------------------ RCS Ident Info -------------------------------
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
!------------------------------------------------------------------------------

contains

  ! --------------------------------------  Find_Scattering_Point  -----
  subroutine Find_Scattering_Point ( Scat_Zeta, Scat_Phi, Scat_Ht, Xi, &
                                   & Scat_Index, Tan_Phi, Npc, Npf, &
                                   & Print_Incopt, Print_IncRad, &
                                   & Print_Path, Req_S, PhiTan, Tan_Ht, &
                                   & FwdModelConf, MAF, Z_Coarse, Vert_Inds, &
                                   & H_Path, Phi_Path, Tan_Pt_C, Tan_Pt_F, &
                                   & Forward, Rev, Which )

    use Constants, only: PI, Rad2Deg
    use Dump_0, only: Dump
    use ForwardModelConfig, only: ForwardModelConfig_T
    use GLNP, only: NG, NGP1
    use MLSKinds, only: RP, R8
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
    use MLSStringLists, only: SwitchDetail
    use Output_m, only: Output
    use Toggles, only: Switches
    use VectorsModule, only: VectorValue_T

    real(rp), intent(in) :: Scat_Zeta
    real(rp), intent(in) :: Scat_Phi         ! of scattering point
    real(rp), intent(in) :: Scat_Ht          ! To check we hit the
                                             ! scattering point
    real(rp), intent(in) :: Xi               ! Scattering angle
    integer, intent(out) :: Scat_Index       ! in coarse grid
    real(rp), intent(in) :: Tan_Phi          ! orbit angle at tangent, radians
    integer, intent(in) :: Npc, Npf          ! Numbers of points on coarse,
                                             ! fine paths

    ! Stuff only for debugging, printing, and error messages
    logical, intent(in) :: Print_Incopt      ! For debugging, from -Sincp
    logical, intent(in) :: Print_IncRad      ! For debugging, from -Sincr
    integer, intent(in) :: Print_Path        ! Nicer format than Print_Incopt,
                                             ! for few molecules
    real(rp), intent(in) :: Req_S            ! Equivalent Earth Radius at
                                             ! height reference surface
    type(vectorValue_t), intent(in) :: PhiTan
    real(rp), intent(in) :: Tan_Ht           ! Geometric tangent height,
                                             ! km from equivalent Earth center
    type(forwardModelConfig_T), intent(in) :: FwdModelConf
    integer, intent(in) :: MAF               ! MAF under consideration

    ! Recomputed path stuff
    real(rp), intent(inout) :: Z_Coarse(:)   ! Zetas on coarse path
    integer, intent(inout) :: Vert_Inds(:)   ! Height indices of fine path in
                                             ! H_Glgrid etc.
    real(rp), intent(inout) :: H_Path(:)     ! Heights on fine path (km)
    real(rp), intent(inout), target :: Phi_Path(:) ! Phi's on fine path, Radians
    integer, intent(inout) :: Tan_Pt_c, Tan_Pt_f ! Index of tangent point in
                                             ! coarse, fine path

    ! Optional stuff
    logical, intent(in), optional :: Forward ! For subsurface rays
                                             ! and Generate_TScat
    logical, intent(in), optional :: Rev     ! Reverse the path
    integer, intent(in), optional :: Which   ! Which TScat call got us here?
                                             ! Negative if we're trying to
                                             ! trace an Earth-intersecting ray.

    character(max(128,len(TScat_Detail_Heading)+10)) :: Line    ! Of output
    logical :: MyRev               ! "Reverse the integration path"
    real(rp), pointer :: Phi_Path_c(:)  ! Phi on the coarse path
    real(rp) :: Scat_D1, Scat_D2   ! To compute scat_index
    integer :: Scat_Temp           ! To compute scat_index
    real(rp), parameter :: Scat_Tol = 1.0 ! max miss of scattering pt, (km/h)**2

    phi_path_c => phi_path(1:npf:ngp1)
    scat_index = 0
    scat_d2 = huge(1.0_rp)
    do scat_temp = 1, npc
      scat_d1 = sin(phi_path_c(scat_temp)-scat_phi)**2 + &
              &  (h_path(scat_temp*ngp1-ng)/scat_ht-1.0_rp)**2
      if ( scat_d1 < scat_d2 ) then
        scat_index = scat_temp
        scat_d2 = scat_d1
      end if
    end do

    if ( present(rev) ) then
      myRev = rev
    else
      myRev = (xi > -pi .and. xi < 0.0) .neqv. (scat_index < tan_pt_c)
    end if
    if ( myRev ) then
      ! Ray from below the scattering point won't see the tangent point, or
      ! Ray from above the scattering point will see the tangent point, so
      ! reverse the path
      z_coarse(:npc) = z_coarse(npc:1:-1)
      vert_inds(:npf) = vert_inds(npf:1:-1)
      h_path(:npf) = h_path(npf:1:-1)
      phi_path(:npf) = phi_path(npf:1:-1)
      scat_index = npc - scat_index + 1
      tan_pt_c = npc - tan_pt_c ! The first one, not the same one
      tan_pt_f = ngp1 * tan_pt_c - ng
    end if
    if ( scat_index == tan_pt_c + 1 ) scat_index = tan_pt_c

    ! Check that we hit it well enough.
    if ( abs(z_coarse(scat_index)-scat_zeta) > 0.05 .or. &
      &  scat_d2 > scat_tol .or. &
      ! or print details if requested
      &  print_TScat_detail > -1 ) then
      phi_path_c => phi_path(1:npf:ngp1)
      if ( print_TScat_detail /= 0 .or. print_TScat_Deriv > -1 .or. &
        & print_incopt .or. print_incrad .or. print_path > -1 ) &
        & call output ( TScat_Detail_Heading, advance="yes" )
      write ( line, "(f7.2,f12.4,f9.3,f11.2,f10.6,f9.2,f11.4,f8.2,i4,f8.2,i4,4x,L1)" ) &
        & rad2deg*scat_phi, scat_ht, scat_zeta, rad2deg*xi, sqrt(scat_d2), & ! km/ht, not (km/ht)**2
        & rad2deg*tan_phi, tan_ht, rad2deg*phi_path(1), scat_index*ngp1-ng, &
        & rad2deg*phi_path_c(scat_index), which, myRev
      if ( present(forward) ) then; line = trim(line) // merge("   T", "   F", forward)
      else ; line = trim(line) // "   T"
      end if
      if ( abs(z_coarse(scat_index)-scat_zeta) > 0.05 .or. &
        &  scat_d2 > scat_tol ) then
        ! Are we trying to trace an Earth-intersecting ray, for TScat
        ! generation, but the detailed Metrics calculation discovered
        ! it's not one?
        line(len(TScat_detail_heading)+1:) = "> abandoned"
        call output ( trim(line), advance="yes" ) ! Too far away from scattering point
        if ( scat_ht <= tan_ht .or. which > 0 .or. print_TScat_detail > -1 ) then
          call output ( req_s, before="Req_s = " )
          call output ( mod(phitan%values(FwdModelConf%TScatMIF,MAF),360.0_r8), &
            & before=", Phi_ref = ", advance="yes" )
          call dump ( z_coarse(:npc), name="Z_Coarse" )
          call dump ( rad2deg*phi_path_c, name="Phi_Path", format="(f14.8)" )
          call dump ( h_path(1:npf:ngp1), name="H_Path", format="(f14.6)" )
        end if
        if ( scat_ht <= tan_ht .or. which > 0 ) &
          call MLSMessage ( merge(MLSMSG_Warning,MLSMSG_Error,switchDetail(switches,'igsc')>0), &
            & moduleName, 'Scattering point appears not to be in path' )
        scat_index = -1
      else
        call output ( trim(line) // " <", advance="yes" ) ! Close enough
        if ( print_TScat_detail > 0 ) then
          call output ( tan_pt_c, before='Tan_Pt_C = ' )
          call dump ( rad2deg*phi_path_c(:scat_index), name=', Phi_Path_C' )
!c              call dump ( h_path(:scat_index*ngp1-ng:ngp1), name='H_Path_C' )
          associate ( h_path_x => h_path(1:npf:ngp1) )
            call dump ( h_path_x(:scat_index), name='H_Path_C' )
          end associate
          call dump ( z_coarse(:scat_index), name='Z_Coarse' )
        end if
      end if
    end if
  end subroutine Find_Scattering_Point

  ! --------------------------------------------------  Get_dB_df  -----
  subroutine Get_dB_df ( Alpha, Beta_c_e, dBeta_c_a_dIWC, &
                       & dBeta_c_s_dIWC, dAlpha_df, W0, Molecule, dB_df, &
                       & dTScat_df )

    use MLSKinds, only: Rp
    use Molecules, only: L_CloudIce

    real(rp), intent(in) :: Alpha(:)     ! on the path
    real(rp), intent(in) :: Beta_c_e(:)  ! on the path
    real(rp), intent(in) :: dBeta_c_a_dIWC(:)  ! on the path, w.r.t. IWC on the path
    real(rp), intent(in) :: dBeta_c_s_dIWC(:)  ! on the path, w.r.t. IWC on the path
    real(rp), intent(in) :: dAlpha_df(:) ! on the path, w.r.t.
                                         ! f == Molecule on the path
    real(rp), intent(in) :: W0(:)        ! on the path
    integer, intent(in) :: Molecule      ! Molecule
    real(rp), intent(out) :: dB_df(:)    ! on the path, w.r.t. f on the path
    real(rp), intent(in), optional :: dTScat_df(:) ! on the path, w.r.t.
                                         ! f == Molecule on the path

!{ \raggedright
!   Compute some or all of {\tt dB_df(i)} =
!   $\frac{\partial \overline{B}_i}{\partial f^k_i}$ on the path,
!           where  $\overline{B}_i = (1-\omega_{0_i}) B_i +
!                                    \omega_{0_i} T_{\text{scat}_i}$
!           and    $\frac{\partial \overline{B_i}}
!                        {\partial f^k_i} =
!    \frac{\partial \omega_{0_i}}{\partial f^k_i}
!     \left( T_{\text{scat}_i} - B_i \right) +
!    \omega_{0_i} \frac{\partial T_{\text{scat}_i}}
!                      {\partial f^k_i}$.
!   The last term is optional, in case it is available.
!   The caller then computes
!   $\frac{\partial \overline{B}_i}{\partial f^k_{lm}}
!   = \frac{\partial \overline{B}_i}{\partial f^k_i}
!     \frac{\partial f^k_i}{\partial f^k_{lm}}
!   = \frac{\partial \overline{B}_i}{\partial f^k_i}
!     \eta^k_{lm}(\zeta_i)$ or
!   $\frac{\partial \overline{B}_i}{\partial f^k_{lm}}
!   = $ {\tt dB_df(i)}$\frac{\partial f^k_i}{\partial f^k_{lm}} +
!   \frac{\partial T_{\text{scat}_i}}{\partial f^k_{lm}} =$
!   {\tt dB_df(i)}$\eta^k_{lm}(\zeta_i) +
!   \frac{\partial T_{\text{scat}_i}}{\partial f^k_{lm}}$.
!   The $k=$ IWC and $k \neq$ IWC cases for
!   $\frac{\partial \omega_{0_i}}{\partial f^k_i}$ simplify
!   differently:
! 
!   \begin{equation}\begin{split}\label{twelve}
!   \frac{\partial \omega_{0_i}}{\partial f^{\text{IWC}}_i}
!   =\,&
!     \omega_{0_i}
!     \left ( \frac1{\beta_{{c\_s}_i}} - \frac1{\beta_{e_i}} \right )
!      \frac{\partial \beta_{{c\_s}_i}}{\partial f^{\text{IWC}}_i}
!     -\frac{\omega_{0_i}}{\beta_{e_i}}
!       \frac{\partial \beta_{{c\_a}_i}}{\partial f^{\text{IWC}}_i}
!   = \frac{1-\omega_{0_i}}{\beta_{e_i}}
!      \frac{\partial \beta_{{c\_s}_i}}{\partial f^{\text{IWC}}_i}
!     -\frac{\omega_{0_i}}{\beta_{e_i}}
!       \frac{\partial \beta_{{c\_a}_i}}{\partial f^{\text{IWC}}_i}
!   \\
!   %
!   \frac{\partial \omega_{0_i}}{\partial f^{k \neq \text{IWC}}_i}
!   =\,&
!    - \frac{\omega_{0_i}}{\beta_{e_i}}
!      \frac{\partial \alpha_{\text{gas}_i}}
!           {\partial f^{k \neq \text{IWC}}_i}\,,
!   \end{split}\end{equation}
!   where $\beta_{e_i} = \alpha_{\text{gas}_i} + \beta_{c\_e_i}$
!   (see {\tt wvs-095}).
!
!   From $\frac{\partial \overline{B}_i}{\partial f^k_{lm}}$ the caller can
!   compute $\frac{\partial\overline{\Delta B}_i}{\partial f^k_{lm}} =
!   \Delta \frac{\partial\overline{B}_i}{\partial f^k_{lm}}$. This is
!   different from {\tt Get_dB_dT}, which computes
!   $\frac{\partial\overline{\Delta B}_i}{\partial T_{lm}}$, because there
!   is only one $T$.

    select case ( molecule )
    case ( l_cloudIce )
      dB_df = ( (1.0-w0) * dBeta_c_s_dIWC - w0 * dBeta_c_a_dIWC ) / &
            & (alpha + beta_c_e)
    case default
      dB_df = -w0 / (alpha + beta_c_e) * dAlpha_df
    end select

    if ( present(dTScat_df) ) dB_df = dB_df + w0 * dTScat_df

  end subroutine Get_dB_df

  ! --------------------------------------------------  Get_dB_dT  -----
  subroutine Get_dB_dT ( Alpha, B, T, Nu, Beta_c_e, dAlpha_dT, &
    &                    dBeta_c_a_dT, dBeta_c_s_dT, TScat, dTScat_dT, W0, &
    &                    Eta_zxp, I_Start, I_Stop, d_delta_B_dT )

!{ \raggedright
!   {\tt d_delta_B_dT(i)} =
!   $\frac{\partial \overline{\Delta B_i}}{\partial T_{lm}}$
!           where  $\overline{B}_i = (1-\omega_{0_i}) B_i +
!                                    \omega_{0_i} T_{\text{scat}_i} $
!           and    $\frac{\partial \overline{B_i}}
!                        {\partial T_{lm}(\zeta_i)} =
!    \frac{\partial \omega_{0_i}}{\partial T_{lm}(\zeta_i)}
!     \left( T_{\text{scat}(\zeta_i)} - B_i \right) +
!     \left( 1 - \omega_{0_i} \right) \frac{\partial B}
!                                          {\partial T_{lm}(\zeta_i)} +
!    \omega_{0_i} \frac{\partial T_{\text{scat}_i}}
!                      {\partial T_{lm}(\zeta_i)}
!   = \left( \frac{\partial \omega_{0_i}}{\partial T}
!     \left( T_{\text{scat}_i} - B_i \right) +
!     \left( 1 - \omega_{0_i} \right) \frac{\partial B}
!                                          {\partial T_i}
!     \right)
!     \frac{\partial T}{\partial T_{lm}(\zeta_i)} +
!     \omega_{0_i} \frac{\partial T_{\text{scat}_i}}{\partial T_{lm}(\zeta_i)}
!   = \left( \frac{\partial \omega_{0_i}}{\partial T}
!      \left( T_{\text{scat}_i} - B_i \right) +
!     \left( 1 - \omega_{0_i} \right) \frac{\partial B}
!                                          {\partial T_i}
!      \right) \eta^T_{lm}(\zeta_i) +
!     \omega_{0_i} \frac{\partial T_{\text{scat}_i}}{\partial T_{lm}(\zeta_i)}$.
!
! $\frac{\partial \omega_{0_i}}{\partial T_i} =
!   \frac{1 - \omega_{0_i}}{\beta_{e_i}}
!    \frac{\partial \beta_{{c\_s}_i}}{\partial T_i} -
!     \frac{\omega_{0_i}}{\beta_{e_i}}
!                   \left( \frac{\partial \beta_{{c\_a}_i}}{\partial T_i} +
!                          \frac{\partial\alpha_{\text{gas}_i}}{\partial T_i}
!                   \right) =
!   \frac1{\beta_{e_i}} \left[
!   \frac{\partial \beta_{{c\_s}_i}}{\partial T_i} -
!     \omega_{0_i} \left( \frac{\partial \beta_{{c\_a}_i}}{\partial T_i} +
!                         \frac{\partial\alpha_{\text{gas}_i}}{\partial T_i} +
!                         \frac{\partial \beta_{{c\_s}_i}}{\partial T_i}
!                   \right) \right]$
! where $\beta_{e_i} = \alpha_{\text{gas}_i} + \beta_{c\_e(\zeta_i)}$
! (see {\tt wvs-095}).
!
! $\frac{\partial \omega_{0_i}}{\partial T_{lm}(\zeta_i)} =
!  \frac{\partial \omega_{0_i}}{\partial T_i}
!  \frac{\partial T_i}{\partial T_{lm}(\zeta_i)} =
!  \frac{\partial \omega_{0_i}}{\partial T_i} \eta^T_{lm}(\zeta_i)$.
!
!  B satisfies the differential equation
!  $\frac{\text{d} B}{\text{d} T} =
!   \frac{B}{T^2} \left( \frac{h\nu}k + B \right)$.

    use d_t_script_dtnp_m, only: dT_Script
    use MLSKinds, only: Rp, R8
    use Physics, only: H_O_K => h_over_k
    use Sparse_m, only: Sparse_t

    real(rp), intent(in) :: Alpha(:)     ! on the path
    real(rp), intent(in) :: B(:)         ! T_Script in some notes
    real(rp), intent(in) :: T(:)         ! Temperature on the path
    real(r8), intent(in) :: Nu           ! Frequency
    real(rp), intent(in) :: Beta_c_e(:)  ! on the path
    real(rp), intent(in) :: dAlpha_dT(:) ! on the path
    real(rp), intent(in) :: dBeta_c_a_dT(:)  ! on the path
    real(rp), intent(in) :: dBeta_c_s_dT(:)  ! on the path
    real(rp), intent(in) :: TScat(:)     ! on the path
    real(rp), intent(in) :: dTScat_dT(:,:) ! dTScat on the path / dT everywhere
    real(rp), intent(in) :: W0(:)        ! on the path
    class(sparse_t), intent(in) :: Eta_zxp ! Interpolating coefficients from
                                           ! state vector to path
    integer, intent(in) :: I_Start, I_Stop ! Range of dB_dT etc to use
    real(rp), intent(out) :: d_delta_B_dT(:,:) ! on the path, w.r.t. T on the grid

    real(rp) :: dB_dT(size(TScat)) ! on the path, w.r.t. T on the path
    integer :: I, Np

    dB_dT = 1.0 / ( alpha + beta_c_e ) * & ! 1 / beta_e
          & ( dBeta_c_s_dT - & 
          &   w0 * ( dAlpha_dT + dBeta_c_a_dT + dBeta_c_s_dT ) ) * & ! d w0 / dT
          & ( TScat - B ) + &
          & ( 1.0 - W0 ) * B * ( h_o_k * nu + B ) / T**2

    call dt_script ( dB_dT, eta_zxp, i_start, i_stop, d_delta_B_dT )

    !{ add $\Delta \omega_{0_i} \frac{\partial T_{\text{scat}_i}}
    !                                {\partial T_{lm}(\zeta_i)}$
    np = size(w0)
    d_delta_B_dT(1,:) = d_delta_B_dT(1,:) + &
                      & 0.5 * ( W0(1) * dTScat_dT(1,:) + W0(2) * dTScat_dT(2,:) )
    do i = 2, np-1  ! Assume size(w0) == size(dTScat_dT,1) ==
                    ! size(d_delta_B_dT,1)
      d_delta_B_dT(i,:) = d_delta_B_dT(i,:) + &
                        & 0.5 * ( W0(i+1) * dTScat_dT(i+1,:) - &
                        &         W0(i-1) * dTScat_dT(i-1,:) )
    end do
    d_delta_B_dT(1,:) = d_delta_B_dT(1,:) - &
                      & 0.5 * ( W0(np-1) * dTScat_dT(np-1,:) + &
                      &         W0(np) * dTScat_dT(np,:) )

  end subroutine Get_dB_dT

  ! --------------------------------------------------  Get_TScat  -----
  subroutine Get_TScat ( FwdModelIn, FwdModelExtra, PhiWindow, MAF, Phitan, &
    &                    Frq, Z_path, Phi_path, Tan_PT, Grids_f, L2PC, dX,  &
    &                    TScat, RadInL2PC, Grids_TScat, TScat_Path,         &
    &                    Atmos_Der, Temp_Der, dTScat_df, dTScat_dT )

  ! Get TScat and its derivatives along the integration path using
  ! the TScat tables and the same strategy as the linear forward model.

    use Comp_Eta_Docalc_Sparse_m, only: Comp_Eta
    use Comp_Sps_Path_Sparse_m, only: Comp_1_Sps_Path_Sparse_Frq
    use Constants, only: Deg2Rad
    use ForwardModelVectorTools, only: GetQuantityForForwardModel
    use Intrinsic, only: L_Temperature
    use L2PC_m, only: L2PC_t
    use Load_Sps_Data_m, only: DestroyGrids_t, Grids_T, Fill_Grids_2, &
      & Load_One_Item_Grid
    use MatrixModule_1, only: Dump, GetVectorFromColumn, &
      & MultiplyMatrixVectorNoT
    use MLSKinds, only: RP, R8
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
    use Output_m, only: Output
    use String_Table, only: Get_String
    use Sparse_Eta_m, only: Sparse_Eta_t
    use Toggles, only: Switches
    use VectorsModule, only: Dump, VectorValue_T, Vector_T

    type(vector_t), intent(in) :: FwdModelIn, FwdModelExtra ! Forward model state vectors.
    real(r8), intent(in) :: PhiWindow   ! Tolerance for TScat phi - state phi
    integer, intent(in) :: MAF          ! Major frame index
    type(vectorValue_T), intent(in) :: Phitan ! Limb tangent geodAngle
    real(r8), intent(in) :: Frq         ! Frequency at which to evaluate TScat
    real(rp), intent(in) :: Z_Path(:)   ! Zeta along the path
    real(rp), intent(in) :: Phi_Path(:) ! Phi along the path
    integer, intent(in) :: Tan_PT       ! Index of tangent point in Z_Path
    type (Grids_T), intent(in) :: Grids_f ! All the coordinates for VMR
    type(L2PC_t), intent(in) :: L2PC    ! The selected L2PC
    type(vector_T), intent(inout) :: dX ! Clone of L2PC%col%vec.  Used to
                                        ! compute value of X - L2PC%col%vec
                                        ! where X is the forward model's
                                        ! state.  Work space to avoid cloning;
                                        ! (inout) to avoid nullifying pointer
                                        ! components.
    type(vector_T), intent(inout) :: TScat ! Same form as L2PC%row%vec. Work
                                        ! space to avoid cloning; (inout) to
                                        ! avoid nullifying pointer components.
    integer, intent(in) :: RadInL2PC    ! Index of qty in TScat, X, dX, L2PC%J
    type(grids_t), intent(inout) :: Grids_TScat ! All the grids for TScat
    logical, intent(in) :: Atmos_Der, Temp_Der  ! Want these derivatives?
    real(rp), intent(out) :: TScat_Path(:) ! TScat on the path.
    real(rp), intent(out) :: dTScat_df(:,:) ! dTScat on path / df everywhere
                                        ! Path X forward model SV
    real(rp), intent(out) :: dTScat_dT(:,:) ! dTScat on path / dT everywhere
                                        ! Path X forward model temperature SV

    real(r8) :: DeltaPhi
    logical :: Do_Calc_FZP(size(z_path), size(grids_TScat%values) )
    logical :: Do_Calc_ZP(size(z_path), grids_TScat%p_len)
    real(rp) :: dTScat_L2PC(size(z_path)) ! dTScat/dX on forward model's path
    logical :: DumpDX, DumpJ, DumpTScat
    type(sparse_eta_t) :: Eta_FZP
    type(sparse_eta_t) :: Eta_ZP
    integer :: F_Inda, F_Indb
    integer :: I               ! Loop index and subscript
    integer :: Instance        ! For which instance in L2PC do we want dTScat?
    real(rp) :: Phi_Offset     ! To shift TScat phi to path phi
    integer :: Qty             ! Loop index and subscript
    type(vectorValue_t), pointer :: SQ ! State vector quantity
    integer :: SQCenter        ! Center instance of SQ
    integer :: SQInstance      ! Instance of SQ
    integer :: SV_f, SV_zp     ! Loop inductor and subscript
    type(vectorValue_t), pointer :: TQ ! L2PC TScat state vector quantity
    integer :: TQCenter        ! center instance of TQ instances
    integer :: TQInstance      ! Instance of TQ
    integer :: V_Ind
    integer :: W_Inda, W_Indb
    character(len=80) :: Word  ! From the string table

    dumpDX = index(switches,"DTSCDX") > 0
    dumpTScat = index(switches,"DTSCAT") > 0
    dumpJ = index(switches,"DTSCJ") > 0

    ! Create empty interpolators
    call eta_zp%create ( size(z_path), grids_TScat%p_len, 2*size(z_path) )
    call eta_fzp%create ( size(z_path), grids_TScat%l_v(1), 2*size(z_path) )

    do qty = 1, size(l2pc%j%col%vec%quantities,1) ! L2PC state quantities
      ! Get quantity in forward model state corresponding to quantity in
      ! L2PC TScat's column vector
      tq => l2pc%j%col%vec%quantities(qty)  ! L2PC state quantity
      sq => getQuantityForForwardModel( fwdModelIn, fwdModelExtra, &
        & quantityType=tq%template%quantityType )
      if ( .not. associated(sq) ) then ! Complain
        call get_string ( tq%template%name, word, strip=.true. )
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          &  "No quantity in state vector to match L2PC TScat quantity "//trim(word) )
      end if
      ! Compute dX using instances of the forward model state vector most nearly
      ! corresponding to instances of the L2PC TScat state vector.
      sqCenter = sq%template%noInstances/2 + 1
      tqCenter = tq%template%noInstances/2 + 1
      do tqInstance = 1, tq%template%noInstances
        sqInstance = max(1,min(sqCenter + tqInstance - tqCenter, &
                        &      sq%template%noInstances ) )
        do ! phiWindowLoop
          deltaPhi = deg2rad*(sq%template%phi(1,sqInstance) - &
                            & tq%template%phi(1,tqInstance) )
          if ( cos(deltaPhi) >= cos(phiWindow) ) exit
          if ( sin(deltaPhi) > 0 ) then
            sqInstance = sqInstance - 1
            if ( sqInstance <= 0 ) exit
          else
            sqInstance = sqInstance + 1
            if ( sqInstance > sq%template%noInstances ) exit
          end if
        end do ! phiWindowLoop
        sqInstance = max(1, min( sqInstance, sq%template%noInstances ) )
        if ( cos(deltaPhi) < cos(phiWindow) ) then
          call output ( sq%template%phi(1,sqInstance), before="State instance phi = " )
          call output ( tq%template%phi(1,tqInstance), &
                      & before=", TScat instance phi = ", advance="yes" )
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
            &  "| State instance phi - TScat instance phi | > phiWindow" )
        end if
        ! Fill the tqInstance part of dX
        dX%quantities(qty)%values(:,tqInstance) = sq%values(:,sqInstance) - &
          &                                       tq%values(:,tqInstance)
      end do ! tqInstance
    end do ! qty
    if ( dumpDX ) call dump ( dX, name='DX' )

    ! Get the TScat quantity for radInL2PC.  We assume TQ is the same quantity
    ! put into Grids_TScat by Get_TScat_Setup.
    tq => TScat%quantities(radInL2PC)
    tq%values = l2pc%j%row%vec%quantities(radInL2PC)%values

    ! Evaluate TScat%quantities(radInL2PC) using first-order approximation.
    ! TScat(out) = TScat(L2PC) + J * dX

    if ( dumpTScat ) call dump ( tq, name='Before' )
    if ( dumpJ ) call dump ( l2pc%J, name="L2PC%J", details=9, row=radInL2PC )
    call MultiplyMatrixVectorNoT ( l2pc%J, dX, TScat, update = .true., row=radInL2PC )
    if ( dumpTScat ) call dump ( tq, name='After' )

    ! Offset so phi in grids_TScat will match up with phi_path.
    phi_offset = sum(phitan%values(:,maf)) / size(phitan%values,1) &
      &        - 0.5 * ( tq%template%phi(1,grids_TScat%windowStart(1)) &
      &                + tq%template%phi(1,grids_TScat%windowFinish(1)) )

    ! Put values and coordinates into Grids_TScat, which was created by
    ! Get_TScat_Setup.  We assume TQ is the same quantity put into
    ! Grids_TScat by Get_TScat_Setup.
    call Fill_Grids_2 ( grids_TScat, 1, tq, setDerivFlags = .true., &
      & phi_offset=phi_offset )

    ! Get eta_zp to put grids_TScat onto TScat_path
    call comp_eta ( grids_TScat, tan_pt, z_path, phi_path, eta_zp )
    ! Get TScat_path.  sideband=0 means absolute frequency, not I.F.
    ! This also returns eta_fzp for putting TScat from the L2PC coordinates
    ! onto the forward model's path and frequency.
    call comp_1_sps_path_sparse_frq ( grids_TScat, frq, eta_zp, eta_fzp, &
      & TScat_path, lo=0.0_r8, sideband=0 )

    ! We have TScat on the path.  Now get the derivatives.

    if ( atmos_der ) dTScat_df = 0.0
    if ( temp_der ) dTScat_dT = 0.0

    if ( atmos_der .or. temp_der ) then
      do qty = 1, size(l2pc%j%col%vec%quantities)
        do instance = 1, size(l2pc%j%col%vec%quantities(qty)%values)
          ! Use TScat vector for a temp, a column of the L2PC Jacobian.
          ! Here, it's dTScat(radInL2PC)(everywhere) / dSomething(qty,phi,zeta)
          call getVectorFromColumn ( l2pc%j, qty, instance, TScat, rowQty=radInL2PC )
          call eta_fzp%sparse_dot_vec ( TScat%quantities(radInL2PC)%value1, &
                                      & dTScat_L2PC )

          ! We have derivatives of TScat w.r.t. one L2PC state vector
          ! quantity and instance on the forward model's path and frequency
          ! but the L2PC state vector.
          ! Now map them into the forward model's state vector
          if ( temp_der .and. &
            &  l2pc%j%col%vec%quantities(qty)%template%quantityType == l_temperature ) then
            dTScat_dT(:,instance) = dTscat_L2PC
          else if ( atmos_der ) then
            do i = 1, size(grids_f%qty)
              if ( l2pc%j%col%vec%quantities(qty)%template%quantityType == &
                 & grids_f%qty(i) ) then
                 dTScat_df(:,grids_f%l_v(i-1)+instance) = dTScat_L2PC
                 exit
               end if
             end do
          end if
        end do ! Instance
      end do ! Qty
    end if ! atmos_der .or. temp_der

  end subroutine Get_TScat

  ! --------------------------------------------  Get_TScat_Setup  -----
  subroutine Get_TScat_Setup ( FmConf, FwdModelIn, FwdModelExtra,       &
    &                          FwdModelOut, Sideband, Frq, MAF, Phitan, &
    &                          L2PC, RadInL2PC, dX, TScat, Grids_TScat )

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test, &
      & Test_Allocate, Test_DeAllocate
    use ForwardModelConfig, only: ForwardModelConfig_t
    use ForwardModelVectorTools, only: GetQuantityForForwardModel
    use Intrinsic, only: L_Radiance, L_TScat
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    use Load_Sps_Data_m, only: Grids_T, Load_One_Item_Grid
    use L2PC_m, only: L2PC_t, L2PCDatabase, PopulateL2PCbin
    use L2PCBins_m, only: SelectL2PCBins
    use MLSKinds, only: R8
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use VectorsModule, only: CloneVector, GetVectorQuantityIndexByType, &
      & VectorValue_T, Vector_T

    type(forwardModelConfig_T), intent(in) :: FMconf
    type(vector_T), intent(in) ::  FWDMODELIN
    type(vector_T), intent(in) ::  FWDMODELEXTRA
    type(vector_T), intent(in) ::  FWDMODELOUT
    integer, intent(in) :: SIDEBAND       ! Which sideband (see below)
    real(r8) :: Frq                       ! The frequency
    integer, intent(in) :: MAF            ! MAF index
    type(vectorValue_T), intent(in) :: Phitan ! Limb tangent geodAngle
    integer :: S                          ! Size in bytes of an object to deallocate
    type(L2PC_t), pointer :: L2PC         ! The selected L2PC
    integer, intent(out) :: RadInL2PC     ! Qty in L2PC for TScat_conf%signals(1)
    type(vector_T), intent(out) :: dX     ! Same form as L2PC%col%vec
    type(vector_T), intent(out) :: TScat  ! Same form as L2PC%row%vec
    type(grids_t), intent(inout) :: Grids_TScat ! All the grids for TScat

    integer, dimension(-1:1) :: L2PCBINS  ! Result
    type (VectorValue_T), pointer :: RADIANCE ! A radiance for the frequency
                                          ! for the TScat we need

    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: Stat

    type(forwardModelConfig_T) :: TScat_Conf

    ! Create an ersatz ForwardModelConfig to select TScat bins

    TScat_conf = fmConf
    TScat_conf%lockBins = .false.
    TScat_conf%xStar = 0
    TScat_conf%yStar = 0

    allocate ( TScat_conf%signals(1), stat=stat )
    addr = 0
    if ( stat == 0 ) addr = transfer(c_loc(TScat_conf%signals(1)), addr)
    call test_allocate ( stat, moduleName, "TScat_Conf%signals", (/1/), (/1/), &
      & storage_size(TScat_conf%signals) / 8, address=addr )

    ! Use the same signals as FmConf%signals(1), but with all channels turned on
    TScat_conf%signals(1) = FmConf%signals(1)
    nullify ( TScat_conf%signals(1)%channels )
    call allocate_test ( TScat_conf%signals(1)%channels, &
      & size(FmConf%signals(1)%channels), "TScat_conf%signals(1)%channels", &
      & moduleName )
    TScat_conf%signals(1)%channels = .true.

    ! Find a radiance for the desired frequency, so we'll have a signal
    ! for selectL2PCBins.  We really want a frequency, but TScat tables
    ! are catalogued by signal, not frequency.
    radiance => GetQuantityForForwardModel ( fwdModelOut, &
      & quantityType=l_radiance, frq=frq )

    ! Determine the best bin to use
    call selectL2PCBins ( TScat_conf, fwdModelIn, fwdModelExtra, &
      & radiance, sideband, maf, l2pcBins )

    ! Load the bin
    call populateL2PCBin ( l2pcBins(sideband), ignoreHessian=.true. )

    ! Get the L2PC from the bin
    l2pc => l2pcDatabase(l2pcBins(sideband))

    radInl2pc = GetVectorQuantityIndexByType ( &
      & l2pc%j%row%vec, quantityType=l_TScat, &
      & signal=TScat_conf%signals(1)%index, sideband=sideband )
    if ( l2pc%j%row%vec%quantities(radInL2PC)%template%noChans /= &
       & radiance%template%noChans ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Channel dimension in l2pc not same as in measurements" )

    call cloneVector ( dX, l2pc%j%col%vec, vectorNameText='__dX' )
    call cloneVector ( TScat, l2pc%j%row%vec, vectorNameText='__TScat' )

    ! Create and fill grids for TScat radiance
    TScat_Conf%phiWindow = size(TScat%quantities(radInl2pc)%template%phi,2)
    call load_one_item_grid ( grids_TScat, TScat%quantities(radInl2pc), &
      & MAF, phitan, TScat_Conf, setDerivFlags = .true. )

    call deallocate_test ( TScat_Conf%signals(1)%channels, &
      & "TScat_Conf%signals%channels", moduleName )
    s = size(TScat_Conf%signals) * storage_size(TScat_Conf%signals) / 8
    addr = 0
    if ( s > 0 ) addr = transfer(c_loc(TScat_Conf%signals(1)), addr)
    deallocate ( TScat_Conf%signals, stat=stat )
    call test_deallocate ( stat, moduleName, "TScat_Conf%signals", s, address=addr )

  end subroutine Get_TScat_Setup

  ! -----------------------------------------  Get_TScat_Teardown  -----
  subroutine Get_TScat_Teardown ( dX, TScat, Grids_TScat )

    ! Clean up the ersatz ForwardModelConfig created by Get_TScat_Setup

    use Load_Sps_Data_m, only: DestroyGrids_t, Grids_T
    use VectorsModule, only: DestroyVectorInfo, Vector_t

    type(vector_t), intent(inout) :: dX     ! Same form as L2PC%col%vec
    type(vector_t), intent(inout) :: TScat  ! Same form as L2PC%row%vec.
    type(grids_t), intent(inout) :: Grids_TScat

    call destroyVectorInfo ( dX )
    call destroyVectorInfo ( TScat ) 

    call destroyGrids_t ( grids_TScat )

  end subroutine Get_TScat_Teardown

  ! ---------------------------------  Interpolate_P_to_theta_e  -----
  subroutine Interpolate_P_to_theta_e ( P, Eta, Theta_e, Beg_Pos_Theta, &
    &                                   Xis, Coeffs, P_on_Xi )

    ! Interpolate P or dP_dT or dP_dIWC from Theta to Xi, for TScat
    ! generation, using periodic spline interpolation.

    use MLSKinds, only: Rp
    use MLSNumerics, only: Coefficients, InterpolateValues
    use Sparse_m, only: Sparse_t

    real(rp), intent(in), contiguous, target :: P(:,:,:) ! P or dP_dT or dP_dIWC
                                       ! on T x IWC x Theta
    class(sparse_t), intent(in) :: Eta ! to interpolate to T, IWC
    real(rp), intent(in) :: Theta_e(:) ! Theta extended to negative values
    integer, intent(in) :: Beg_Pos_Theta ! 1 or 2, depending on whether
                                       ! theta_e(1) is nonzero or zero
    real(rp), intent(in) :: Xis(:) ! Angles on which radiative transfer calculated
    type(coefficients(rp)), intent(in) :: Coeffs ! To interpolate from Theta_e to Xis
    real(rp), intent(out) :: P_on_Xi(:) ! Interpolated result

    real(rp) :: P_T_IWC(size(p,3)) ! P interpolated to T and IWC
    real(rp) :: P_on_Theta_e(2*size(p_t_iwc)+1-beg_pos_theta) ! P * sin(abs(theta))
                             ! interpolated to scattering point IWC and T for
                             ! each theta (including negative theta), intermediary
                             ! to getting P_On_Xi
    real(rp), pointer :: T_X_IWC(:) ! One-D pointer to P(:,:,i)

    integer :: I

    ! Interpolate P to scattering point IWC and temperature, and
    ! multiply that by sin(theta).  This is done one P_T_IWC at a time because
    ! Eta only has one row, because it interpolates to one IWC and
    ! Temperature.
    do i = 1, size(p_t_iwc)
      t_x_iwc(1:size(p,1)*size(p,2)) => p(:,:,i) ! Get 1-D pointer
      ! P_T_IWC(i) = sparse Eta .dot. T_X_IWC
      p_t_iwc(i) = eta%row_dot_vec ( 1, t_x_iwc ) * abs(sin(theta_e(i)))
    end do

    ! The phase function is symmetric in theta, but the radiances
    ! are not.  Unfold p_t_iwc to negative theta.
    p_on_theta_e(1:size(p_t_iwc)) = p_t_iwc(size(p_t_iwc):1:-1)
    p_on_theta_e(size(p_t_iwc)+1:2*size(p_t_iwc)+1-beg_pos_theta) = &
      & p_t_iwc(beg_pos_theta:)

    ! Now interpolate P_On_Theta_e to P_On_Xi using a periodic spline
    call interpolateValues ( coeffs, theta_e, p_on_theta_e, &
      & xis, p_on_xi, method='S', extrapolate='P' )

  end subroutine Interpolate_P_to_theta_e

  ! ---------------------------------------------  Mie_Freq_Index  -----
  integer function Mie_Freq_Index ( Frq, Tol )
    use MLSKinds, only: R8
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use Output_m, only: Output
    use Read_Mie_m, only: F_s
    real(r8), intent(in) :: Frq, Tol
    ! Choose which frequency panel of phase function to use (don't
    ! bother interpolating because the change as a result of the change
    ! of frequency within a single radiometer is small).  Make sure it's
    ! close enough.
    Mie_freq_index = minloc(abs(Frq-f_s),1)
    if ( abs(Frq-f_s(Mie_freq_index)) > tol ) then
      call output ( Mie_freq_index, before="Mie_freq_index = " )
      call output ( f_s(Mie_freq_index), before=", f_s(Mie_freq_index) = " )
      call output ( Frq, before=", Frq = " )
      call output ( tol, before=", tol = ", advance='yes' )
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'In TScat computation, phase function frequency coordinate too far from desired frequency' )
    end if
  end function Mie_Freq_Index

  ! --------------------------------------------  TScat_Gen_Setup  -----
  subroutine TScat_Gen_Setup ( FwdModelConf, FwdModelIn, FwdModelExtra,     &
    &                          FwdModelOut, NoUsedChannels, Temp, Sideband, &
    &                          IWC, ScatteringAngles )

    ! Get ready for computation of TScat table

    use ForwardModelConfig, only: Channels_T, ForwardModelConfig_t
    use ForwardModelVectorTools, only: GetQuantityForForwardModel
    use Intrinsic, only: L_ScatteringAngle, L_TSCAT, L_VMR
    use Molecules, only: L_CLOUDICE
    use Read_Mie_m, only: dP_dT, F_s, P
    use VectorsModule, only: Vector_T, VectorValue_T

    ! Inputs:
    type (forwardModelConfig_t), intent(in) :: FwdModelConf
    type (vector_T), intent(in) ::  FwdModelIn, FwdModelExtra
    type (vector_T), intent(in) :: FwdModelOut  ! Radiances, etc.
    type (vectorValue_T), intent(in) :: TEMP ! Temperature component of state vector
    integer, intent(in) :: NoUsedChannels, Sideband
    ! Outputs:
    type (vectorValue_T), pointer :: IWC  ! IWC at scattering points
    type (vectorValue_T), pointer :: ScatteringAngles ! for TScat computation

    logical :: Error
    integer :: Q ! Loop index
    type (Channels_T), pointer, dimension(:) :: Channels
    type (VectorValue_T), pointer :: TScat, TScat2
    if ( .not. allocated(p) .or. .not. allocated(F_s) .or. &
      &  .not. allocated(dP_dT) ) &
      & call announce_error ( 'TScat table computation requires Mie tables.' )
    if ( FwdModelConf%do_conv ) call announce_error ( &
      & 'Convolution and TScat computation are incompatible.' )
    if ( FwdModelConf%incl_cld ) call announce_error ( &
      & 'Cloud modeling and TScat computation are incompatible.' )
    if ( FwdModelConf%polarized ) call announce_error ( &
      & 'Cannot compute TScat for the polarized model.' )
    if ( FwdModelConf%refract ) call announce_error ( &
      & 'Refractive correction and TScat computation are incompatible.' )
    if ( FwdModelConf%spect_der ) call announce_error ( &
      & 'Spectroscopy derivatives and TScat computation are incompatible.' )
    channels => fwdModelConf%channels
    ! Get IWC field
    iwc => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra,  &
      & quantityType=l_vmr, molecule=l_cloudIce, config=fwdModelConf )
    ! Make sure all the TScat quantities for the selected signals have the
    ! same grids
    TScat => GetQuantityForForwardModel ( &
        & fwdModelOut, quantityType=l_TScat, &
        & signal=fwdModelConf%signals(channels(1)%signal)%index, &
        & sideband=sideband )
    if ( .not. TScat%template%coherent .or. &
      &  .not. TScat%template%stacked ) call announce_error ( &
        & 'TScat coordinates must be stacked and coherent' )
    TScat%values = 0.0
    error = .false.
    do q = 2, noUsedChannels
      TScat2 => GetQuantityForForwardModel ( &
        & fwdModelOut, quantityType=l_TScat, &
        & signal=fwdModelConf%signals(channels(q)%signal)%index, &
        & sideband=sideband )
      error = .not. TScat2%template%coherent .or. &
        &     .not. TScat2%template%stacked .or. &
        &     TScat%template%noChans /= TScat2%template%noChans .or. &
        &     TScat%template%hGridIndex /= TScat2%template%hGridIndex .or. &
        &     TScat%template%vGridIndex /= TScat2%template%vGridIndex
      if ( error ) exit
      TScat2%values = 0.0
    end do
    if ( error ) call announce_error ( "TScat quantities must all have the same grids." )
    ! Make sure IWC and temperature are coherent and stacked
    if ( .not. temp%template%stacked .or. .not. temp%template%coherent .or. &
       & .not.  iwc%template%stacked .or. .not.  iwc%template%coherent ) &
       & call announce_error ( "IWC and temperature must be coherent and stacked." )
    ! Make sure IWC and temperature have the same grids as TScat
    if ( TScat%template%hGridIndex /= temp%template%hGridIndex .or. &
       & TScat%template%vGridIndex /= temp%template%vGridIndex .or. &
       & TScat%template%hGridIndex /=  iwc%template%hGridIndex .or. &
       & TScat%template%vGridIndex /=  iwc%template%vGridIndex ) &
       & call announce_error ( "TScat, IWC and temperature must have the same grids." )
    scatteringAngles => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_scatteringAngle, config=fwdModelConf )

  contains

  ! .............................................  Announce_Error  .....
    subroutine Announce_Error ( Message )
    ! Announce Message using MLSMessage.  Include the configuration name
    ! in the message.  This is the same as in FullForwardModel.
      use MLSMessageModule, only: MLSMessage, MLSMSG_Error
      use String_Table, only: Get_String, String_Length
      character(len=*), intent(in) :: Message
      integer, parameter :: C = len('With config(')
      integer :: L
      character(511) :: Work ! Should be plenty of room
      l = string_length(fwdModelConf%name)
      work(:c) = 'With config('
      call get_string ( fwdModelConf%name, work(c+1:) )
      l = c + l + 1
      work(l:l+3) = '): '
      l = l + 3
      work(l+1:l+len_trim(message)) = message
      call MLSMessage ( MLSMSG_Error, moduleName, work(:l+len_trim(message)) )
    end subroutine Announce_Error

  end subroutine TScat_Gen_Setup

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

end module TScat_Support_m

! $Log$
! Revision 2.18  2018/09/12 23:49:14  vsnyder
! Convert to sparse interpolators
!
! Revision 2.17  2018/05/17 02:15:45  vsnyder
! Use sparse instead of dense interpolation
!
! Revision 2.16  2018/05/15 03:26:25  vsnyder
! Change Mie tables from pointer to allocatable
!
! Revision 2.15  2018/05/14 23:40:58  vsnyder
! Change to sparse eta representation
!
! Revision 2.14  2017/10/31 23:49:36  vsnyder
! Make Coefficients a parameterized type
!
! Revision 2.13  2016/11/11 02:38:41  vsnyder
! Associate the Phi_Path_C pointer before using it
!
! Revision 2.12  2016/09/30 01:26:30  vsnyder
! Create Find_Scattering_Point andd move stuff to it from FullForwardModel
!
! Revision 2.11  2016/05/10 00:03:09  vsnyder
! Cannonball polishing
!
! Revision 2.10  2015/03/28 02:07:09  vsnyder
! Added stuff to trace allocate/deallocate addresses
!
! Revision 2.9  2014/09/05 20:54:48  vsnyder
! More complete and accurate allocate/deallocate size tracking
!
! Revision 2.8  2013/06/12 02:24:01  vsnyder
! Cruft removal
!
! Revision 2.7  2011/07/29 01:53:28  vsnyder
! Make CloudIce a molecule.  Look for CloudIce instead of Cloud_A and
! Cloud_S.  Correct derivative computations w.r.t. CloudIce.
!
! Revision 2.6  2011/03/25 20:46:59  vsnyder
! Delete declarations of unused objects
!
! Revision 2.5  2011/03/23 23:51:38  vsnyder
! Completely revised Get_dB_df, revised some TeXnicalities in Get_dB_dT
!
! Revision 2.4  2011/03/23 23:45:32  vsnyder
! This log entry is bogus.  Check in again to get the right one.
! FOV_Convolve_m.f90
!
! Revision 2.3  2011/01/28 19:16:28  vsnyder
! Lots of stuff for derivatives
!
! Revision 2.2  2010/09/25 01:09:01  vsnyder
! Inching along
!
! Revision 2.1  2010/08/27 23:39:13  vsnyder
! Initial commit -- far from working
!
