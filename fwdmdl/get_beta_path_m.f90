! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module GET_BETA_PATH_M

  implicit NONE
  private
  public :: Create_Beta, Create_Beta_Path, Create_Beta_Path_PFA
  public :: Get_Beta_Path, Get_Beta_Path_Cloud, Get_Beta_Path_PFA
  public :: Get_Beta_Path_Polarized, Get_Beta_Path_Scalar

  interface Get_Beta_Path
    module procedure Get_Beta_Path_Scalar, Get_Beta_Path_Polarized
  end interface

  ! RHi = H2O * exp(RHa - RHb/T) / P
  real, parameter, private :: RHa = 19.6783
  real, parameter, private :: RHb = 6140.4

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: ModuleName= &
    &  "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains

! ---------------------------------------  Get_Beta_Path_Scalar  -----
  subroutine Get_Beta_Path_Scalar ( frq, Flo, p_path, t_path, tanh_path, &
        & beta_group, SX, NoPolarized, gl_slabs, path_inds,         &
        & beta_path, t_der_path_flags, dTanh_dT, VelCor,            &
        & dBeta_dt_path, dBeta_dw_path, dBeta_dn_path, dBeta_dv_path, &
        & dBeta_df_path, Where_dBeta_df, Sps_Path )

    use Dump_0, only: Dump
    use ForwardModelConfig, only: Beta_Group_T, LineCenter, LineWidth, &
      & LineWidth_Tdep
    use Intrinsic, only: Lit_Indices
    use MLSKinds, only: R8, RP, IP
    use MLSStringLists, only: SwitchDetail
    use Molecules, only: L_RHI
    use Output_M, only: Output
    use SLABS_SW_M, only: SLABS_Struct
    use String_Table, only: Display_String
    use Toggles, only: Emit, Levels, Switches, Toggle
    use Trace_M, only: Trace_Begin, Trace_End

! Inputs:

    real(r8), intent(in) :: Frq          ! frequency in MHz
    real(r8), intent(in) :: Flo          ! LO frequency in MHz
    real(rp), intent(in) :: P_path(:)    ! pressures on combined coarse & fine
                                         ! path in hPa.  Use the subset
                                         ! selected by path_inds.
    real(rp), intent(in) :: T_path(:)    ! path temperatures, already selected
                                         ! from combined coarse & fine path by
                                         ! path_inds.
    real(rp), intent(in) :: Tanh_path(:) ! tanh(0.5*h_over_k*frq / t_path)
    type (slabs_struct), dimension(:,:), intent(in) :: Gl_slabs ! Single-line
                                         ! absorption data on the combined
                                         ! coarse & fine path.  Use the subset
                                         ! selected by path_inds.
    integer(ip), intent(in) :: Path_inds(:) ! indices in combined coarse & fine
                                         ! path for reading p_path and gl_slabs.

    type (beta_group_t), intent(in), dimension(:) :: beta_group
    integer, intent(in) :: SX            ! Sideband index, 1 or 2.

    logical, intent(in) :: NoPolarized   ! "Don't work on Zeeman-split lines"

! Optional inputs.

    logical, intent(in) :: t_der_path_flags(:) ! where temperature derivatives
!                               are needed. Only useful for subsetting.
    real(rp), intent(in) :: dTanh_dT(:)    ! dTanh( (-h nu) / (k T) ) / dT on path

! Input for spectroscopy derivatives
    real(rp), intent(in) :: VelCor

! Input for H2O continuum.  Optional so PFA table generation can leave it off.
    real(rp), intent(in), optional :: Sps_Path(:,:) ! Mixing ratios, Path X Species

! Outputs

    real(rp), intent(out) :: beta_path(:,:) ! path beta for each species, km^{-1}

! Optional outputs.  We use size(.)>0 instead of PRESENT so that the
! caller doesn't need multiple branches.

    real(rp), target, intent(out) :: dBeta_dT_path(:,:) ! Temperature
    real(rp), target, intent(out) :: dBeta_dw_path(:,:) ! line width
    real(rp), target, intent(out) :: dBeta_dn_path(:,:) ! line width t dep.
    real(rp), target, intent(out) :: dBeta_dv_path(:,:) ! line position
    real(rp), target, intent(out) :: dBeta_df_path(:,:) ! mixing ratio
    ! Which column of dBeta_df_path to use for a molecule.  Do not compute if
    ! zero.
    integer, intent(in) :: Where_dBeta_df(:)

! Local variables.

    character(len=4), save :: clean
    real(rp), pointer :: dBdf(:), dBdn(:), dBdT(:), dBdv(:), dBdw(:) ! slices of dBeta_d*_path
    integer, save :: DumpLevel ! units digit of lblb switch
                               ! >0 dump beta, >1 dump T, Tanh
    logical, save :: DumpStop  ! hundreds digit of lblb switch > 0
    logical, save :: DumpZeta  ! tens digit of lblb switch > 0, dump Zeta
                               ! instead of P
    real(rp) :: ES(size(t_path)) ! Used for RHi calculation
    logical, save :: First = .true. ! Fist-time flag
    real(rp) :: Sps(size(t_path))
    integer(ip) :: I, N, NP
    integer :: Me = -1         ! String index for trace

! begin the code

    if ( first ) then
      first = .false.
      dumpLevel = switchDetail(switches,'lblb')
      dumpStop = dumpLevel > 99
      dumpZeta = mod(dumpLevel/10,10) > 0
      if ( dumpLevel >= 0 ) dumpLevel = mod(dumpLevel,10)
      ! clean = switchDetail(switches,'clean') > -1
      clean = ' '
      if ( switchDetail(switches,'clean') > -1 ) clean = 'c'
    end if

    call Trace_Begin ( me, 'ForwardModel.Get_Beta_Path_Scalar', &
      & cond=toggle(emit) .and. levels(emit) > 6 )

    nullify ( dBdf, dBdn, dBdT, dBdv, dBdw ) ! Disassociated means "Don't compute it"

    np = size(beta_path,1)

    do i = 1, size(beta_group)
      if ( beta_group(i)%lbl(sx)%spect_der_ix(lineWidth_tDep) /= 0 .and. size(dBeta_dn_path) > 0 ) then
        dBdn => dBeta_dn_path(:np,beta_group(i)%lbl(sx)%spect_der_ix(lineWidth_tDep))
        dBdn = 0.0_rp
      end if
      if ( beta_group(i)%lbl(sx)%spect_der_ix(lineCenter) /= 0 .and. size(dBeta_dv_path) > 0 ) then
        dBdv => dBeta_dv_path(:np,beta_group(i)%lbl(sx)%spect_der_ix (lineCenter))
        dBdv = 0.0_rp
      end if
      if ( beta_group(i)%lbl(sx)%spect_der_ix(lineWidth) /= 0 .and. size(dBeta_dw_path) > 0 ) then
        dBdw => dBeta_dw_path(:np,beta_group(i)%lbl(sx)%spect_der_ix(lineWidth))
        dBdw = 0.0_rp
      end if

      beta_path(:,i) = 0.0_rp
      if ( size(dBeta_dt_path) > 0 ) then
        dBdT => dBeta_dt_path(:np,i)
        dBdT = 0.0_rp
      end if

      if ( size(dBeta_df_path) > 0 ) then
        if ( where_dBeta_df(i) > 0 ) then
          dBdf => dBeta_df_path(:np,where_dBeta_df(i))
          dBdf = 0.0
        end if
      end if

      if ( present(sps_path) ) sps = sps_path(path_inds,i)
      if ( beta_group(i)%molecule == l_rhi ) then
        !{ Get ready for RHi conversion
        ! \begin{equation}\begin{split}
        ! e_s =\,& \frac1P \exp\left(19.6783-\frac{6140.4}T\right)\\
        ! \beta^{\text{RHi}} =\,& \beta^{\text{H}_2\text{O}} e_s\\
        ! \frac{\partial \beta^{\text{RHi}}}{\partial T} =\,&
        !  e_s \frac{\partial \beta^{\text{H}_2\text{O}}}{\partial T} +
        !  \frac{6140.4}{T^2} \beta^{\text{RHi}}
        ! \end{split}\end{equation}
        es = exp(RHa - RHb/t_path) / p_path(path_inds)
        if ( present(sps_path) ) then
          sps = es * sps
          where ( p_path(path_inds) < 50.0 ) sps = 0.0
        end if
      end if

      do n = 1, size(beta_group(i)%lbl(sx)%cat_index)
        if ( present(sps_path) ) then
          call create_beta_path ( path_inds, p_path, t_path, frq, flo,  &
            & beta_group(i)%lbl(sx)%ratio(n),                       &
            & gl_slabs(:,beta_group(i)%lbl(sx)%cat_index(n)),       &
            & tanh_path, noPolarized, velCor,                       &
            & beta_path(:,i), dTanh_dT, t_der_path_flags,           &
            & dBdT, dBdw, dBdn, dBdv, dBdf, sps_path=sps )
        else
          call create_beta_path ( path_inds, p_path, t_path, frq, flo,  &
            & beta_group(i)%lbl(sx)%ratio(n),                       &
            & gl_slabs(:,beta_group(i)%lbl(sx)%cat_index(n)),       &
            & tanh_path, noPolarized, velCor,                       &
            & beta_path(:,i), dTanh_dT, t_der_path_flags,           &
            & dBdT, dBdw, dBdn, dBdv, dBdf )
        end if

        if ( dumpLevel > 0 ) then
          call display_string ( lit_indices(beta_group(i)%lbl(sx)%molecules(1:n)), &
            & before='LBL Betas for' )
          call output ( frq, before=', FRQ = ', advance='yes' )
          call dump ( beta_path(:,i), name='beta', options=clean )
          if ( associated(dBdT) ) call dump ( dBdT, name='dBdT', options=clean )
          if ( associated(dBdf) ) call dump ( dBdf, name='dBdf', options=clean )
        end if
      end do

      if ( beta_group(i)%molecule == l_rhi ) then
        ! Convert specific humidity to relative humidity?
        beta_path(:,i) = beta_path(:,i) * es
        if ( associated(dBdT) ) &
          & dBdT = es * dBdT + RHb / t_path**2 * beta_path(:,i)

        if ( dumpLevel > 0 ) then
          call display_string ( lit_indices(l_rhi), before='LBL Betas for ' )
          call display_string ( lit_indices(beta_group(i)%lbl(sx)%molecules), &
            & before=' using' )
          call dump ( beta_path(:,i), options=clean )
          if ( associated(dBdT) ) call dump ( dBdT, name='dBdT', options=clean )
        end if
      end if

    end do ! i = 1, size(beta_group)

    call Trace_End ( 'ForwardModel.Get_Beta_Path_Scalar', &
      & cond=toggle(emit) .and. levels(emit) > 6 )

    if ( dumpLevel > 1 ) then
      if ( dumpZeta ) then
        call dump ( -log10(p_path(path_inds)), name='Zetas', options=clean )
      else
        call dump ( p_path(path_inds), name='Pressures', options=clean )
      end if
      call dump ( t_path, name='Temperatures', options=clean )
      call dump ( tanh_path, name='tanh(h nu / k T)', options=clean )
    end if
    if ( dumpStop ) stop

  end subroutine Get_Beta_Path_Scalar

  ! ------------------------------------------  Get_Beta_Path_PFA  -----
  subroutine Get_Beta_Path_PFA ( Frq, Frq_i, P_Path, Path_Inds, T_Path, &
    & Beta_Group, SX, Vel_Rel, Sps_Path, Beta_Path, T_Der_Path_Flags, &
    & dBeta_dT_Path, dBeta_dw_Path, dBeta_dn_Path, dBeta_dv_Path, &
    & dBeta_dIWC_Path )

    use Dump_0, only: Dump
    use ForwardModelConfig, only: Beta_Group_T
    use Intrinsic, only: Lit_Indices
    use MLSKinds, only: RP, R8
    use MLSStringLists, only: SwitchDetail
    use Molecules, only: L_Cloudice, L_RHI
    use Output_M, only: Output
    use PFADatabase_M, only: PFAData
    use Read_Mie_M, only: Beta_C_A, DBeta_DIWC_C_A, DBeta_DT_C_A, &
                        & Beta_C_S, DBeta_DIWC_C_S, DBeta_DT_C_S
    use String_Table, only: Display_String
    use Toggles, only: Switches

! Inputs
    real(r8), intent(in) :: Frq         ! Channel center frequency in MHz
    integer, intent(in) :: Frq_i        ! Channel index
    real(rp), intent(in) :: P_path(:)   ! path pressures in hPa!
    integer, intent(in) :: Path_inds(:) ! indicies for reading P_path
    real(rp), intent(in) :: T_path(:)   ! path temperatures
    type(beta_group_t), intent(in) :: Beta_Group(:) ! PFA stuff for the beta group
    integer, intent(in) :: SX           ! Sideband index, 1 or 2.
    real(rp), intent(in) :: Vel_Rel     ! LOS Vel / C
    real(rp), intent(in) :: Sps_Path(:,:) ! VMR's, for nonlinear species

! Output
    real(rp), intent(out) :: beta_path(:,:) ! path beta for each species, km^{-1}

! Optional input.

    logical, intent(in) :: T_Der_Path_Flags(:) ! where temperature derivatives
!                               are needed. Only useful for subsetting.

! Optional outputs.  We use size(...)/=0 instead of PRESENT so that the
! caller doesn't need multiple branches.

    real(rp), target, intent(out) :: dBeta_dT_path(:,:)   ! Temperature
    real(rp), target, intent(out) :: dBeta_dw_path(:,:)   ! line width
    real(rp), target, intent(out) :: dBeta_dn_path(:,:)   ! line width t dep.
    real(rp), target, intent(out) :: dBeta_dv_path(:,:)   ! line position
    real(rp), target, intent(out) :: dBeta_dIWC_path(:,:) ! IWC

    real(rp), pointer :: dBdIWC(:), dBdn(:), dBdT(:), dBdv(:), dBdw(:) ! slices of dBeta_d*_path
    real(rp) :: ES(size(t_path)) ! Used for RHi calculation
    integer :: I, J, NP
    character(len=4), save :: clean
    logical, save :: DumpAll, DumpBeta, DumpStop
    logical, save :: First = .true. ! First-time flag
    integer :: PFAB ! Dump level 0 => Beta, 1 => Temperature + pressure, 2 => stop

    if ( first ) then
      first = .false.
      pfab = switchDetail(switches,'pfab')
      dumpStop = pfab > 1
      dumpAll = pfab > 0
      dumpBeta = pfab >= 0
      ! clean = switchDetail(switches,'clean') > -1
      clean = ' '
      if ( switchDetail(switches,'clean') > -1 ) clean = 'c'
    end if

    nullify ( dBdIWC, dBdT, dBdn, dBdv, dBdw )

    np = size(t_path)

    do i = 1, size(beta_group)
      if ( size(dBeta_dt_path) > 0 ) then
        dBdT => dBeta_dt_path(:np,i)
        dBdT = 0.0_rp
      end if
      if ( size(dBeta_dn_path) > 0 ) then
        dBdn => dBeta_dn_path(:np,i)
        dBdn = 0.0_rp
      end if
      if ( size(dBeta_dv_path) > 0 ) then
        dBdv => dBeta_dv_path(:np,i)
        dBdv = 0.0_rp
      end if
      if ( size(dBeta_dw_path) > 0 ) then
        dBdw => dBeta_dw_path(:np,i)
        dBdw = 0.0_rp
      end if
      if ( size(dBeta_dIWC_path) > 0 ) then
        dBdIWC => dBeta_dIWC_path(:np,i)
        dBdIWC = 0.0_rp
      end if

      beta_path(:,i) = 0.0_rp

      select case ( beta_group(i)%molecule )
      case ( l_cloudIce )
        call create_beta_path_Mie ( frq, t_path, sps_path(:,i), path_inds, &
          & beta_c_a, dBeta_dIWC_c_a, dBeta_dT_c_a, beta_path(:,i), &
          & dBdT, dBdIWC )
        call create_beta_path_Mie ( frq, t_path, sps_path(:,i), path_inds, &
          & beta_c_s, dBeta_dIWC_c_s, dBeta_dT_c_s, beta_path(:,i), &
          & dBdT, dBdIWC )
      case default
        do j = 1, size(beta_group(i)%pfa(sx)%molecules)
          if ( beta_group(i)%pfa(sx)%data(frq_i,j)/= 0 ) then
            call create_beta_path_PFA ( frq, p_path, path_inds, t_path, vel_rel, &
              & PFAData(beta_group(i)%pfa(sx)%data(frq_i,j)),   &
              & beta_group(i)%pfa(sx)%ratio(j), beta_path(:,i), &
              & t_der_path_flags, dBdT, dBdw, dBdn, dBdv )

            if ( dumpBeta ) then
              call display_string ( lit_indices(beta_group(i)%pfa(sx)%molecules(1:j)), &
              & before='PFA Betas for' )
              call output ( frq, before=', FRQ = ', advance='yes' )
              call dump ( beta_path(:,i), name='Beta', options=clean )
              if ( associated(dBdT) ) call dump ( dBdT, name='dBdT', options=clean )
            end if
          end if
        end do ! j = 1, size(beta_group(i)%pfa(sx)%molecules)

        if ( beta_group(i)%molecule == l_rhi ) then
          !{ Convert specific humidity to relative humidity
          ! \begin{equation}\begin{split}
          ! e_s =\,& \frac1P \exp\left(19.6783-\frac{6140.4}T\right)\\
          ! \beta^{\text{RHi}} =\,& \beta^{\text{H}_2\text{O}} e_s\\
          ! \frac{\partial \beta^{\text{RHi}}}{\partial T} =\,&
          !  e_s \frac{\partial \beta^{\text{H}_2\text{O}}}{\partial T} +
          !  \frac{6140.4}{T^2} \beta^{\text{RHi}}
          ! \end{split}\end{equation}
          es = exp(RHa - RHb/t_path) / p_path(path_inds)
          beta_path(:,i) = beta_path(:,i) * es
          if ( associated(dBdT) ) &
            & dBdT = es * dBdT + RHb / t_path**2 * beta_path(:,i)

          if ( dumpBeta ) then
            call display_string ( lit_indices(l_rhi), before='PFA Betas for ' )
            call display_string ( lit_indices(beta_group(i)%pfa(sx)%molecules), &
              & before=' using' )
            call dump ( beta_path(:,i), name='Beta', options=clean )
            if ( associated(dBdT) ) call dump ( dBdT, name='dBdT', options=clean )
          end if
        end if
      end select

    end do ! i = 1, size(beta_group)

    if ( dumpAll ) then
      call dump ( p_path(path_inds), name='Pressures', options=clean )
      call dump ( t_path, name='Temperatures', options=clean )
    end if
    if ( dumpStop ) stop

  end subroutine Get_Beta_Path_PFA

  ! ------------------------------------  Get_Beta_Path_Polarized  -----
  subroutine Get_Beta_Path_Polarized ( Frq, H, Beta_group, GL_slabs, &
                                     & Path_inds, Beta_path, dBeta_path_dT )

    use Dump_0, only: Dump
    use ForwardModelConfig, only: LBL_T
    use Intrinsic, only: Lit_Indices
    use MLSKinds, only: R8, RP, IP
    use MLSStringLists, only: SwitchDetail
    use O2_ABS_CS_M, only: O2_ABS_CS, D_O2_ABS_CS_DT
    use Output_M, only: Output
    use SLABS_SW_M, only: SLABS_Struct
    use String_Table, only: Display_String
    use Toggles, only: Switches

! Inputs:

    real(r8), intent(in) :: Frq       ! Frequency in MHz
    real(rp), intent(in) :: H(:)      ! Magnetic field component in instrument
                                      ! polarization on the path
    type (slabs_struct), dimension(:,:), intent(in) :: GL_slabs
    integer(ip), intent(in) :: Path_inds(:) ! indicies for reading gl_slabs
    type (LBL_T), dimension(:), intent(in) :: Beta_group

! Outputs

!{ The variable {\tt Beta\_Path} isn't really $\beta$.  It lacks a factor of
!  $\tanh\left( \frac{h \nu}{2 k T}\right)$.  This is put in after we
!  compute the weighted average over species, saving as many multiplies as
!  there are species.

    complex(rp), intent(out) :: Beta_path(-1:,:,:) ! path beta for each speciess, km^{-1}
    ! beta_path(-1,:,:) is Sigma_m, beta_path(0,:,:) is Pi,
    ! beta_path(+1,:,:) is Sigma_p

    complex(rp), intent(out) :: dBeta_path_dT(-1:,:,:)

! Local variables..

    integer(ip) :: I, IB, J, K, N, N_PATH
    real(rp) :: RATIO ! Isotope ratio, not mixing ratio
    complex(rp) :: Sigma_m, Pi, Sigma_p
    complex(rp) :: dSigma_m_dT, dPi_dT, dSigma_p_dT
    character(len=4), save :: clean
    logical, save :: DumpBeta, DumpStop
    logical, save :: First = .true. ! First-time flag

    if ( first ) then
      first = .false.
      dumpStop = switchDetail(switches,'POLB') > -1
      dumpBeta = dumpStop .or. ( switchDetail(switches,'polb') > -1 )
      ! clean = switchDetail(switches,'clean') > -1
      clean = ' '
      if ( switchDetail(switches,'clean') > -1 ) clean = 'c'
    end if

! begin the code

    n_path = size(path_inds)

    beta_path = 0.0
    if ( size(dBeta_path_dT) > 0 ) dBeta_path_dT = 0.0

    do i = 1, size(beta_group)
      do n = 1, size(beta_group(i)%cat_index)
        ratio = beta_group(i)%ratio(n)
        ib = beta_group(i)%cat_index(n)

        ! Get beta everywhere on the path
        do j = 1, n_path
          k = path_inds(j)

          if ( size(dBeta_path_dT) == 0 ) then
            call o2_abs_cs ( frq, h(k), gl_slabs(k,ib), &
                           & sigma_p, pi, sigma_m )
          else
            call d_o2_abs_cs_dT ( frq, h(k), gl_slabs(k,ib),    &
                                & sigma_p,     pi,     sigma_m, &
                                & dSigma_p_dT, dPi_dT, dSigma_m_dT )
            dBeta_path_dT(-1,j,i) = dBeta_path_dT(-1,j,i) + ratio * dSigma_m_dT
            dBeta_path_dT( 0,j,i) = dBeta_path_dT( 0,j,i) + ratio * dPi_dT
            dBeta_path_dT(+1,j,i) = dBeta_path_dT(+1,j,i) + ratio * dSigma_p_dT
          end if
          beta_path(-1,j,i) = beta_path(-1,j,i) + ratio * sigma_m
          beta_path( 0,j,i) = beta_path( 0,j,i) + ratio * pi
          beta_path(+1,j,i) = beta_path(+1,j,i) + ratio * sigma_p
        end do ! j

        if ( dumpBeta ) then
          call display_string ( lit_indices(beta_group(i)%molecules(1:n)), &
            & before='Polarized Betas for' )
          call output ( frq, before=', FRQ = ', advance='yes' )
          call dump ( beta_path(:,:,i), name='Beta', options=clean )
          if ( size(dBeta_path_dT) > 0 ) then
            call dump ( real(dBeta_path_dT), name='real(dBeta_path_dT)', options=clean )
            call dump ( aimag(dBeta_path_dT), name='aimag(dBeta_path_dT)', options=clean )
          end if
        end if
      end do ! n
    end do ! i
    if ( dumpStop ) stop

  end subroutine Get_Beta_Path_Polarized

  ! ----------------------------------------  Get_Beta_Path_Cloud  -----
  subroutine Get_Beta_Path_Cloud ( Frq, t_path, tt_path, path_inds, &
        & beta_path_cloud, w0_path, tt_path_c, IPSD, WC, fwdModelConf  )

    use ForwardModelConfig, only: ForwardModelConfig_T
    use Cloud_Extinction, only: Get_Beta_Cloud
    use MLSKinds, only: R8, RP, IP

! Inputs:

    real(r8), intent(in) :: Frq             ! Frequency in MHz
    real(rp), intent(in) :: T_path(:)       ! Path temperatures
    real(rp), intent(in) :: tt_path(:,:)    ! Scating source func on gl grids


    integer(ip), intent(in) :: Path_inds(:) ! indices for reading T_PATH

    type (ForwardModelConfig_T) ,intent(in) :: FWDMODELCONF

    integer, intent(in) :: IPSD(:)
    real(rp), intent(in)  :: WC(:,:)
    real(rp) :: W0       ! SINGLE SCATTERING ALBEDO
    real(rp) :: PHH(fwdModelConf%num_scattering_angles)   ! PHASE FUNCTION

! Outputs

    real(rp), intent(out) :: beta_path_cloud(:) ! cloud extinction, km^{-1}
    real(rp), intent(out) :: w0_path(:)         ! single scattering albedo
    real(rp), intent(out) :: tt_path_c(:)       ! scattering source func coarse grids

! Local variables..

    integer :: NC, NU, NUA, NAB, NR
    integer(ip) :: j, k, n_path
    real(rp) :: cld_ext

! begin the code

    NU  = fwdModelConf%Num_Scattering_Angles
    NUA = fwdModelConf%Num_Azimuth_Angles
    NAB = fwdModelConf%Num_AB_Terms
    NR  = fwdModelConf%Num_Size_Bins

    n_path = size(path_inds)

    beta_path_cloud = 0.0_rp
    w0_path         = 0.0_rp
    tt_path_c       = 0.0_rp

    do j = 1, n_path
      k = path_inds(j)

      call get_beta_cloud ( Frq, t_path(k),                         &
                        &   WC(:,k), IPSD(k), NC, NU, NUA, NAB, NR, &
                        &   cld_ext, W0, PHH                )      

      beta_path_cloud(j) = beta_path_cloud(j) + cld_ext 
      w0_path(j)         = w0_path(j)         + W0 
      tt_path_c(j)       = tt_path_c(j)       + tt_path(k,1)            

    end do

  end subroutine Get_Beta_Path_Cloud

! ----------------------------------------------  Create_beta  ---------

  subroutine Create_beta ( pressure, Temp, Fgr, Flo, slabs_0, tanh1, &
         &                 Beta_Value, NoPolarized, dTanh_dT,   &
         &                 dBeta_dT, dBeta_dw, dBeta_dn, dBeta_dv, Sps, &
         &                 dBeta_df )

!  For a given frequency and height, compute beta_value function.
!  This routine should be called for primary and image separately.

!  If you change Create_Beta, change Create_Beta_Path the same way!
!  Create_Beta and Create_Beta_Path should do the same thing.

!  The reason for the existence of Create_Beta_Path is that by moving the
!  loop for the path down into Create_Beta a substantial improvement in
!  running time was achieved.  Create_Beta is in the inner loop of the
!  forward model.

    use MLSKinds, only: RP, R8, IP
    use Molecules, only: L_N2, L_Extinction, L_ExtinctionV2, L_H2O, &
      & L_MIFExtinction, L_O2
    use SLABS_SW_M, only: DVoigt_Spectral, Voigt_Lorentz, &
      & SLABS_Lines, SLABS_Lines_DT, SLABS_Struct, &
      & SLABSWint_Lines, SLABSWint_Lines_DT
    use SpectroscopyCatalog_M, only: Catalog_T, Lines

! Inputs:
    real(rp), intent(in) :: pressure   ! pressure in hPa
    real(rp), intent(in) :: temp       ! temperature in K
    real(r8), intent(in) :: fgr        ! frequency in MHz
    real(r8), intent(in) :: Flo        ! LO frequency in MHz
    type(slabs_struct), intent(in) :: slabs_0 ! contains, among others:

!    catalog        ! Pointer to spectroscopy catalog
!    v0s(:)         ! pressure shifted line centers
!    x1(:)          ! Doppler width
!    y(:)           ! ratio Pressure to Doppler widths
!    yi(:)          ! Interference coefficients
!    slabs1(:)      ! strengths
!    dslabs1_dv0(:) ! strength derivative wrt line position

    real(rp), intent(in) :: tanh1      ! tanh(h*frq/(2*k*T))

    ! "Don't do line(L) if slabs%catalog%polarized(L)"
    logical, intent(in) :: NoPolarized

! Optional inputs for temperature derivatives
    ! -(h frq) / (2 k T^2) ( tanh( (2 h frq) / (k T) ) - 1/tanh( (h frq) / (2 k T) ) ):
    real(rp), intent(in), optional :: dTanh_dT ! 1/tanh d/dT tanh
! Outputs
    real(rp), intent(out) :: beta_value         ! km^{-1}
! Optional outputs
    real(rp), optional, intent(out) :: DBeta_DT ! Temperature derivative
    real(rp), optional, intent(out) :: DBeta_DW ! line width derivative
    real(rp), optional, intent(out) :: DBeta_DN ! temperature dependence deriv
    real(rp), optional, intent(out) :: DBeta_DV ! line position derivative
! Optional for H2O -- Optional so PFA table generation can ignore it
    real(rp), optional, intent(in) :: SPS       ! Mixing ratio
    real(rp), optional, intent(out) :: DBeta_DF ! mixing ratio derivative

! -----     Local variables     ----------------------------------------

    type(catalog_t), pointer :: Catalog
    real(rp), pointer :: Cont(:)    ! Continuum parameters
    integer(ip) :: LN_I             ! Line index
    integer(ip) :: NL               ! no of lines
    logical :: Spect_Der

    real(rp) :: bv, dNu, dw, dn, ds, dbdt, dbdw, dbdn, dbdv

!----------------------------------------------------------------------------

    spect_der = present(dBeta_dw) .or. present(dBeta_dn) .or. present(dBeta_dv)
    if ( spect_der ) then
      dbdw = 0.0_rp
      dbdn = 0.0_rp
      dbdv = 0.0_rp
    end if

!  Setup absorption coefficients function
!  Now get the beta_value:

    catalog => slabs_0%catalog
    cont => catalog%continuum
    nl = size(catalog%lines)
    select case ( catalog%molecule )
    case ( l_n2 ) ! ...........................................  Dry Air

      if ( present(dBeta_dT) ) then
        call abs_cs_n2_cont_dT ( cont, temp, pressure, fgr, beta_value, dBeta_dT )
      else
        beta_value = abs_cs_n2_cont(cont,Temp,Pressure,Fgr)
      end if

    case ( l_extinction, l_MIFextinction ) ! ...............  Extinction

      beta_value = 1.0_rp
      if ( present(dBeta_dT) ) dBeta_dT = 0.0_rp
      return

    case ( l_extinctionv2 ) ! ............................  ExtinctionV2

      beta_value = ( fgr / flo ) ** 2.0_rp
      if ( present(dBeta_dT) ) dBeta_dT = 0.0_rp
      return

    case ( l_o2 ) ! ................................................  O2

      if ( present(dBeta_dT) ) then
        call abs_cs_o2_cont_dT ( cont, temp, pressure, fgr, beta_value, dBeta_dT )
      else
        beta_value = abs_cs_o2_cont(cont,Temp,Pressure,Fgr)
      end if

    case ( l_h2o ) ! ..............................................  H2O

      if ( present(dBeta_dT) ) then
        call abs_cs_h2o_cont_dT ( cont, temp, pressure, fgr, sps, &
          & beta_value, dBeta_dT, dBeta_df )
      else
        if ( present(dBeta_df) ) then
          call abs_cs_h2o_cont_df ( cont, temp, pressure, fgr, sps, &
          & beta_value, dBeta_df )
        else
          beta_value = abs_cs_h2o_cont(cont,Temp,Pressure,Fgr,sps)
        end if
      end if

    case default ! ..............................................  Other

      if ( present(dBeta_dT) ) then
        call abs_cs_cont_dT ( cont, temp, pressure, fgr, beta_value, dBeta_dT )
      else
        beta_value = abs_cs_cont(cont,Temp,Pressure,Fgr)
      end if

    end select

    if ( nl < 1 ) then
      if ( present(dBeta_dw)) dBeta_dw = 0.0_rp
      if ( present(dBeta_dn)) dBeta_dn = 0.0_rp
      if ( present(dBeta_dv)) dBeta_dv = 0.0_rp
      return
    end if

    if ( spect_der ) then

      do ln_i = 1, nl

        if ( noPolarized ) then
          if ( catalog%polarized(ln_i) ) cycle
        end if

        dNu = Fgr - slabs_0%s(ln_i)%v0s

        if ( abs(slabs_0%s(ln_i)%y)+0.666666_rp*abs(slabs_0%s(ln_i)%x1*dNu) &
        & > 100.0_rp ) then
          call Voigt_Lorentz ( dNu, slabs_0%s(ln_i)%v0s, slabs_0%s(ln_i)%x1,         &
            &  slabs_0%s(ln_i)%yi, slabs_0%s(ln_i)%y, lines(catalog%lines(ln_i))%w,  &
            &  Temp, tanh1, slabs_0%s(ln_i)%slabs1, bv, slabs_0%s(ln_i)%dslabs1_dv0, &
            &  dw, dn, ds )
        else
          call DVoigt_Spectral ( dNu, slabs_0%s(ln_i)%v0s, slabs_0%s(ln_i)%x1,       &
            &  slabs_0%s(ln_i)%yi, slabs_0%s(ln_i)%y, lines(catalog%lines(ln_i))%w,  &
            &  Temp, tanh1, slabs_0%s(ln_i)%slabs1, bv, slabs_0%s(ln_i)%dslabs1_dv0, &
            &  dw, dn, ds )
        end if

        beta_value = beta_value + bv
        dbdw = dbdw + dw
        dbdn = dbdn + dn
        dbdv = dbdv + ds

      end do

      if ( present(dBeta_dw)) dBeta_dw = dbdw
      if ( present(dBeta_dn)) dBeta_dn = dbdn
      if ( present(dBeta_dv)) dBeta_dv = dbdv

    else ! No spectroscopy derivatives required

      if ( .not. present(dBeta_dT) ) then
        if ( slabs_0%useYi ) then
          bv = slabswint_lines ( Fgr, slabs_0, tanh1, noPolarized )
        else
          bv = slabs_lines ( Fgr, slabs_0, tanh1, noPolarized )
        end if
      else ! Temperature derivatives needed
        if ( slabs_0%useYi ) then
          call slabswint_lines_dT ( fgr, slabs_0, tanh1, dTanh_dT, &
            & bv, dbdt, noPolarized )
        else
          call slabs_lines_dT ( fgr, slabs_0, tanh1, dTanh_dT, &
            & bv, dbdt, noPolarized )
        end if
        dBeta_dT = dBeta_dT + dbdt
      end if
      beta_value = beta_value + bv

    end if

  end Subroutine Create_beta

! -----------------------------------------  Create_beta_path  ---------

  pure &
  subroutine Create_beta_path ( Path_inds, Pressure, Temp, Fgr, Flo, Ratio,  &
         &   Slabs_0, Tanh1, NoPolarized, VelCor, Beta_value, dTanh_dT, &
         &   Path_flags, dBeta_dT, dBeta_dw, dBeta_dn, dBeta_dv,        &
         &   dBeta_df, Sps_Path )

!  For a given frequency and height, compute beta_value function. This routine
!  should be called for primary and image separately. Compute dBeta_dT if it's
!  associated.  Compute dBeta_dw, dBeta_dn, dBeta_dv if they're associated. 

    use MLSKinds, only: RP, R8
    use Molecules, only: L_N2, L_Extinction, L_ExtinctionV2, L_H2O, &
      & L_MIFExtinction, L_O2
    use SLABS_SW_M, only: SLABS_Struct, &
      & SLABS_Lines, SLABS_Lines_Dall, SLABS_Lines_DSpectral, SLABS_Lines_DT, &
      & SLABSWint_Lines, &
      & SLABSWint_Lines_Dall, SLABSWint_Lines_DSpectral, SLABSWint_Lines_DT

! Inputs:
    integer, intent(in) :: Path_inds(:)! Which Pressures to use
    real(rp), intent(in) :: Pressure(:)! pressure in hPa on the fine path grid
    real(rp), intent(in) :: Temp(:)    ! temperature in K along the path
    real(r8), intent(in) :: Fgr        ! frequency in MHz
    real(r8), intent(in) :: Flo        ! LO frequency
    real(rp), intent(in) :: Ratio      ! Isotope ratio
    type(slabs_struct), intent(in) :: Slabs_0(:) ! contains, among others:

!    Catalog        ! Pointer to catalog
!    v0s(:)         ! pressure shifted line centers
!    x1(:)          ! Doppler width
!    y(:)           ! ratio Pressure to Doppler widths
!    yi(:)          ! Interference coefficients
!    slabs1(:)      ! strengths
!    dslabs1_dv0(:) ! strength derivative wrt line position

    real(rp), intent(in) :: Tanh1(:)   ! tanh(h*frq/(2*k*T))

    ! "Don't do line L if slabs_0(k)%catalog%Polarized(L)"
    logical, intent(in) :: NoPolarized

    ! For spectroscopy derivatives:
    real(rp), intent(in) :: VelCor

! Outputs
    real(rp), intent(inout) :: Beta_value(:) ! km^{-1}

! Optional inputs for temperature derivatives:
    real(rp), intent(in) :: dTanh_dT(:) ! -h nu / (2 k T^2) 1/tanh(...) dTanh(...)/dT
    logical, intent(in) :: Path_Flags(:) ! to do on fine path -- default true

! Optional outputs.  Disassociated means "do not compute"
    real(rp), pointer :: dBeta_dT(:) ! temperature derivative
    real(rp), pointer :: dBeta_dw(:) ! line width derivative
    real(rp), pointer :: dBeta_dn(:) ! temperature dependence deriv
    real(rp), pointer :: dBeta_dv(:) ! line position derivative
    real(rp), pointer :: dBeta_df(:) ! mixing ratio derivative

! Optional inputs for H2O -- optional so PFA table generation can ignore it
    real(rp), intent(in), optional :: Sps_Path(:) ! Mixing ratios on the
                                     ! current part of the path -- as for
                                     ! temperature, not pressure

! -----     Local variables     ----------------------------------------

    integer :: J, K              ! Subscript, loop inductor
    integer :: NL                ! no of lines
    logical :: Spect_Der         ! Spectroscopy derivatives required
    logical :: Temp_Der          ! Temperature derivatives required

    real(rp) :: bv, dBdf, dBdT, dBdw, dBdn, dBdv

!----------------------------------------------------------------------------

    if ( associated(dBeta_dw) .or. associated(dBeta_dn) .or. associated(dBeta_dv) ) then
      dBdw = 0.0_rp
      dBdn = 0.0_rp
      dBdv = 0.0_rp
    end if

    nl = size(slabs_0(1)%catalog%lines) ! All of the slabs have the same catalog
    spect_der = associated(dBeta_dw) .or. associated(dBeta_dn) .or. &
              & associated(dBeta_dv)

    !ocl temp(k,temp_der,bv,dBdT)
    do j = 1, size(path_inds)
      k = path_inds(j)
      temp_der = associated(dBeta_dT)
      if ( temp_der .and. size(path_flags) > 0 ) temp_der = path_flags(k)

      select case ( slabs_0(k)%catalog%molecule )
      case ( l_n2 ) ! .........................................  Dry Air

        ! This nominally gets multiplied by "ratio**2" but in practice this
        ! function is for all isotopic forms of N2 hence the ratio is one.

        if ( temp_der ) then
          call abs_cs_n2_cont_dT ( slabs_0(k)%catalog%continuum, temp(j), &
            & pressure(k), fgr, bv, dBdT )
          beta_value(j) = beta_value(j) + bv
          dBeta_dT(j) = dBeta_dT(j) + dBdT
        else
          beta_value(j) = beta_value(j) + &
            & abs_cs_n2_cont(slabs_0(k)%catalog%continuum,Temp(j),Pressure(k),Fgr)
        end if

      case ( l_extinction, l_MIFextinction ) ! .............  Extinction

        beta_value(j) = beta_value(j) + ratio
!       if ( temp_der ) dBeta_dT(j) = dBeta_dT(j) + ratio * 0.0
        cycle ! we know there are no spectral lines

      case ( l_extinctionv2 ) ! ............................  Extinction

        beta_value(j) = beta_value(j) + ratio * ( fgr / flo ) ** 2
!       if ( temp_der ) dBeta_dT(j) = dBeta_dT(j) + ratio * 0.0
        cycle ! we know there are no spectral lines

      case ( l_o2 ) ! ..............................................  O2

        if ( temp_der ) then
          call abs_cs_o2_cont_dT ( slabs_0(k)%catalog%continuum, temp(j), &
            & pressure(k), fgr, bv, dBdT )
          beta_value(j) = beta_value(j) + ratio * bv
          dBeta_dT(j) = dBeta_dT(j) + ratio * dBdT
        else
          beta_value(j) = beta_value(j) + ratio * &
            & abs_cs_o2_cont(slabs_0(k)%catalog%continuum,Temp(j),Pressure(k),Fgr)
        end if

      case ( l_h2o ) ! ............................................  H2O

        if ( present(sps_path) ) then
          if ( temp_der ) then
            if ( associated(dBeta_df) ) then
              call abs_cs_h2o_cont_dT ( slabs_0(k)%catalog%continuum, temp(j), &
                & pressure(k), fgr, ratio*sps_path(j), bv, dBdT, dBdf )
              beta_value(j) = beta_value(j) + ratio * bv
              dBeta_dT(j) = dBeta_dT(j) + ratio * dBdT
              dBeta_df(j) = dBeta_df(j) + ratio * dBdf
            else
              call abs_cs_h2o_cont_dT ( slabs_0(k)%catalog%continuum, temp(j), &
                & pressure(k), fgr, ratio*sps_path(j), bv, dBdT )
              beta_value(j) = beta_value(j) + ratio * bv
              dBeta_dT(j) = dBeta_dT(j) + ratio * dBdT
            end if
          else if ( associated(dBeta_df) ) then
            call Abs_CS_H2O_Cont_df ( slabs_0(k)%catalog%continuum, temp(j), &
              & pressure(k), fgr, ratio*sps_path(j), bv, dBdf )
            beta_value(j) = beta_value(j) + ratio * bv
            dBeta_df(j) = dBeta_df(j) + ratio * dBdf
          else
            beta_value(j) = beta_value(j) + ratio * &
              & abs_cs_h2o_cont(slabs_0(k)%catalog%continuum,Temp(j), &
              & Pressure(k),Fgr,ratio*sps_path(j))
          end if
        else ! Same as case default for slabs_0(k)%catalog%molecule
          if ( temp_der ) then
            call abs_cs_cont_dT ( slabs_0(k)%catalog%continuum, temp(j), &
              & pressure(k), fgr, bv, dBdT )
            beta_value(j) = beta_value(j) + ratio * bv
            dBeta_dT(j) = dBeta_dT(j) + ratio * dBdT
          else
            beta_value(j) = beta_value(j) + ratio * &
              & abs_cs_cont(slabs_0(k)%catalog%continuum,Temp(j),Pressure(k),Fgr)
          end if
        end if

      case default ! ............................................  Other

        if ( temp_der ) then
          call abs_cs_cont_dT ( slabs_0(k)%catalog%continuum, temp(j), &
            & pressure(k), fgr, bv, dBdT )
          beta_value(j) = beta_value(j) + ratio * bv
          dBeta_dT(j) = dBeta_dT(j) + ratio * dBdT
        else
          beta_value(j) = beta_value(j) + ratio * &
            & abs_cs_cont(slabs_0(k)%catalog%continuum,Temp(j),Pressure(k),Fgr)
        end if

      end select

      if ( nl < 1 ) then
        if ( associated(dBeta_dw) ) dBeta_dw(j) = 0.0_rp
        if ( associated(dBeta_dn) ) dBeta_dn(j) = 0.0_rp
        if ( associated(dBeta_dv) ) dBeta_dv(j) = 0.0_rp
        cycle
      end if

      if ( .not. temp_der .and. .not. spect_der ) then

        ! Add in sum of betas for all the lines
        if ( slabs_0(k)%useYi ) then
          beta_value(j) = beta_value(j) + &
            & ratio * slabswint_lines ( Fgr, slabs_0(k), tanh1(j), noPolarized )
        else
          beta_value(j) = beta_value(j) + &
            & ratio * slabs_lines ( Fgr, slabs_0(k), tanh1(j), noPolarized )
        end if

      else if ( temp_der .and. .not. spect_der ) then

        if ( slabs_0(k)%useYi ) then
          call slabswint_lines_dT ( fgr, slabs_0(k), tanh1(j), dTanh_dT(j), &
            & bv, dBdT, noPolarized )
        else
          call slabs_lines_dT ( fgr, slabs_0(k), tanh1(j), dTanh_dT(j), &
            & bv, dBdT, noPolarized )
        end if

        beta_value(j) = beta_value(j) + ratio * bv
        dBeta_dT(j) = dBeta_dT(j) + ratio * dBdT

      else if ( .not. temp_der .and. spect_der ) then

        if ( slabs_0(k)%useYi ) then
          call slabswint_lines_dSpectral ( Fgr, slabs_0(k), temp(j), tanh1(j), &
            & velCor, noPolarized, bv, dBdw, dBdn, dBdv )
        else
          call slabs_lines_dSpectral ( Fgr, slabs_0(k), temp(j), tanh1(j), &
            & velCor, noPolarized, bv, dBdw, dBdn, dBdv )
        end if

        beta_value(j) = beta_value(j) + ratio * bv

        if ( associated(dBeta_dw) ) dBeta_dw(j) = dBeta_dw(j) + ratio * dBdw
        if ( associated(dBeta_dn) ) dBeta_dn(j) = dBeta_dn(j) + ratio * dBdn
        if ( associated(dBeta_dv) ) dBeta_dv(j) = dBeta_dv(j) + ratio * dBdv

      else if ( temp_der .and. spect_der ) then

        if ( slabs_0(k)%useYi ) then
          call slabswint_lines_dAll ( Fgr, slabs_0(k), temp(j), tanh1(j), &
            & dTanh_dT(j), velCor, noPolarized, bv, dBdT, dBdw, dBdn, dBdv )
        else
          call slabs_lines_dAll ( Fgr, slabs_0(k), temp(j), tanh1(j), &
            & dTanh_dT(j), velCor, noPolarized, bv, dBdT, dBdw, dBdn, dBdv )
        end if

        beta_value(j) = beta_value(j) + ratio * bv
        dBeta_dT(j) = dBeta_dT(j) + ratio * dBdT

        if ( associated(dBeta_dw) ) dBeta_dw(j) = dBeta_dw(j) + ratio * dBdw
        if ( associated(dBeta_dn) ) dBeta_dn(j) = dBeta_dn(j) + ratio * dBdn
        if ( associated(dBeta_dv) ) dBeta_dv(j) = dBeta_dv(j) + ratio * dBdv

      end if

    end do ! j = 1, size(path_inds)

  end Subroutine Create_beta_path

  ! ---------------------------------------  Create_Beta_Path_Mie  -----
  pure &
  subroutine Create_Beta_Path_Mie ( Frq, T_Path, IWC_Path, Path_Inds, &
    & Beta_Tab, dBeta_dIWC_tab, dBeta_dT_tab, Beta_Path, dBdT, dBdIWC )

  !{Interpolate in Beta\_c\_a and its derivatives to the T and IWC along the
  ! path.
  !
  ! Let $\xi_1 = \frac{x-x_0}{x_1-x_0}$, $\xi_0 = \frac{x_1-x}{x_1-x_0} =
  ! 1-\xi_1$,
  ! $\eta_1 = \frac{y-y_0}{y_1-y_0}$, $\eta_0 = \frac{y_1-y}{y_1-y_0} =
  ! 1-\eta_1$,
  ! $\xi = [\xi_0,\xi_1]^T$, $\eta = [\eta_0,\eta_1]^T$, and
  ! $Z = \left| \begin{array}{cc} z_{00} & z_{01} \\ z_{10} & z_{11}\\
  !             \end{array} \right|$.
  !
  ! The interpolation can be done either by interpolating in $x$ to
  ! $z_c = \xi_0 z_{00} + \xi_1 z_{10}$ and $z_d = \xi_0 z_{01} +
  ! \xi_1 z_{11}$ and then in $y$ to $z = \eta_0 z_c + \eta_1 z_d$,
  ! or interpolating in $y$ to $z_a = \eta_0 z_{00} + \eta_1 z_{01}$
  ! and $z_b = \eta_0 z_{10} + \eta_1 z_{11}$ and then in $x$ to $z =
  ! \xi_0 z_a + \xi_1 z_b$.  The two orders of interpolation are
  ! equivalent, and when expanded give
  !%
  ! $
  ! z(x,y) = \xi_0 \eta_0 z_{00} + \xi_1 \eta_0 z_{10} +
  ! \xi_0 \eta_1 z_{10} + \xi_1 \eta_1 z_{11},
  ! $
  !%
  ! which can be expressed in matrix-vector form as $\xi^T Z \eta$.

    use MLSNumerics, only: PureHunt
    use MLSKinds, only: RP, R8
    use Read_Mie_M, only: F_S, IWC_S, T_S

! Inputs:
    real(r8), intent(in) :: Frq         ! Channel center frequency in MHz
    real(rp), intent(in) :: T_Path(:)   ! Temperature in K along the path
    real(rp), intent(in) :: IWC_Path(:) ! Ice water content along the path
    integer, intent(in) :: Path_Inds(:) ! Which IWC_Path to use
    real(r8), intent(in) :: Beta_Tab(:,:,:)       ! T x IWC x F
    real(r8), intent(in) :: dBeta_dIWC_Tab(:,:,:) ! T x IWC x F
    real(r8), intent(in) :: dBeta_dT_Tab(:,:,:)   ! T x IWC x F

! Outputs:
    real(rp), intent(inout) :: Beta_Path(:)  ! Beta along the path, km^{-1}

! Optional outputs.  These are pointers instead of optional.  Null means
! "don't compute it."  We do this so the caller doesn't need four branches.
    real(rp), pointer :: dBdT(:)       ! dBeta_dT along the path
    real(rp), pointer :: dBdIWC(:)     ! dBeta_dIWC along the path

! -----     Local variables     ----------------------------------------

    real(rp) :: Eta(2,2)  ! Interpolating factors for T x IWC
    integer :: F_I(2)     ! Indices for F_s
    real(rp) :: I_Fac     ! Interpolating factor for IWC
    integer :: I_I(2)     ! Indices in IWC_s
    integer :: J
    real(rp) :: Log10_IWC ! log10(IWC_path(j))
    integer :: N_F        ! Frequency index (3rd) for Beta_c_a etc.
                          ! F_I(1) or F_I(2), whichever is closest to Frq
    real(rp) :: T_Fac     ! Interpolating factor for T
    integer :: T_I(2)     ! Indices in T_s

    f_i = 0 ! Initialize for Hunt
    i_i = 0
    t_i = 0

    call purehunt ( frq, f_s, size(f_s), f_i(1), f_i(2) )
    n_f = merge(f_i(1), f_i(2), frq-f_s(f_i(1)) < frq-f_s(f_i(2)))

    !ocl independent
    !ocl temp
    do j = 1, size(beta_path)
      if ( iwc_path(path_inds(j)) <= tiny(iwc_path(path_inds(j))) ) then
        beta_path(j) = 0.0
        if ( associated(dBdIWC) ) dBdIWC(j) = 0.0
        if ( associated(dBdT) ) dBdT(j) = 0.0
      else
        log10_IWC = log10(iwc_path(path_inds(j)))
        call purehunt ( log10_IWC, iwc_s, size(iwc_s), i_i(1), i_i(2) )
        call purehunt ( t_path(j), t_s, size(t_s), t_i(1), t_i(2) )
        i_fac = (log10_IWC - iwc_s(i_i(1))) / (iwc_s(i_i(2)) - iwc_s(i_i(1)))
        t_fac = (t_path(j) - t_s(t_i(1))) / (t_s(t_i(2)) - t_s(t_i(1)))
        eta = reshape( (/ (1.0 - t_fac) * (1.0 - i_fac), &
                          t_fac         * (1.0 - i_fac), &
                          (1.0 - t_fac) * i_fac        , &
                          t_fac         * i_fac         /), (/ 2, 2 /) )
        beta_path(j) = sum(eta * beta_tab(t_i,i_i,n_f))
        if ( associated(dBdIWC) ) &
          & dBdIWC(j) = sum(eta * dBeta_dIWC_tab(t_i,i_i,n_f))
        if ( associated(dBdT) ) &
          & dBdT(j) = sum(eta * dBeta_dT_tab(t_i,i_i,n_f))
      end if
    end do

  end subroutine Create_Beta_Path_Mie

  ! ---------------------------------------  Create_Beta_Path_PFA  -----
  pure &
  subroutine Create_Beta_Path_PFA ( Frq, P_Path, Path_Inds, T_Path, Vel_Rel, &
    & PFAD, Ratio, Beta_Path, T_Der_Path, dBdT, dBdw, dBdn, dBdv )

    use MLSNumerics, only: PureHunt
    use MLSKinds, only: RP, R8
    use PFADatabase_M, only: PFAData_T

! Inputs:
    real(r8), intent(in) :: Frq         ! Channel center frequency in MHz
    real(rp), intent(in) :: P_Path(:)   ! Log10 ( Pressure in hPa _)
                                        ! on the fine path grid
    integer, intent(in) :: Path_inds(:) ! Which Pressures to use
    real(rp), intent(in) :: T_Path(:)   ! Temperature in K along the path
    real(rp), intent(in) :: Vel_Rel     ! LOS vel/c
    type(PFAData_t), intent(in) :: PFAD ! PFA datum from PFA Database
    real(rp), intent(in) :: Ratio       ! Isotope ratio
    
! Outputs
    real(rp), intent(inout) :: Beta_Path(:)  ! km^{-1}

! Optional inputs for temperature derivatives:
    logical, intent(in) :: T_Der_Path(:)   ! To do on fine path -- default true

! Optional outputs
    real(rp), pointer :: dBdT(:) ! Temperature derivative
    real(rp), pointer :: dBdw(:) ! Line width derivative
    real(rp), pointer :: dBdn(:) ! Temperature dependence deriv
    real(rp), pointer :: dBdv(:) ! Line position derivative

! -----     Local variables     ----------------------------------------

    real(kind(beta_path)) :: BP  ! Temp for one cell of beta_path
    real(rp) :: dBdNu            ! d log Beta / d nu, for Doppler correction
    real(rp) :: Del_T            ! Log Temperature step in tGrid
    real(r8) :: Doppler          ! Doppler corrected frequency offset, MHz
    integer :: J, K
    real(rp) :: LogT             ! Ln ( temperature )
    integer :: P_I1, P_I2        ! Indices in PFAData%vGrid%surfs
    real(rp) :: P_Fac            ! Interpolating factor for Pressure
    logical :: Temp_Der          ! Temperature derivatives required
    integer :: T_I1, T_I2        ! Indices in PFAData%tGrid%surfs
    real(rp) :: T_Fac            ! Interpolating factor for Temperature

    !{ Doppler correction = $\nu_0 \left[ \left( 1 - \frac{v}c \right) -
    !                                     \left( 1 - \frac{v_l}c \right) \right] =
    !                        \nu_0 \left[ \frac{v_l}c - \frac{v}c \right] $

    doppler = frq * ( PFAD%vel_rel - vel_rel )

    p_i1 = 0 ! Initialize for Hunt
    t_i1 = 0

    !ocl independent
    !ocl temp
    do j = 1, size(path_inds)

      k = path_inds(j)
      temp_der = associated(dBdT)
      if ( temp_der .and. size(t_der_path) > 0 ) temp_der = t_der_path(k)

      ! Get interpolating factors
      logT = log(t_path(j))
      call purehunt ( logT, PFAD%tGrid%surfs(:,1), PFAD%tGrid%noSurfs, &
        & t_i1, t_i2 )
      del_t = PFAD%tGrid%surfs(t_i2,1) - PFAD%tGrid%surfs(t_i1,1)
      t_fac = (logT - PFAD%tGrid%surfs(t_i1,1)) / del_t
      call purehunt ( p_path(k), PFAD%vGrid%surfs(:,1), PFAD%vGrid%noSurfs, &
        & p_i1, p_i2 )
      p_fac = (p_path(k) - PFAD%vGrid%surfs(p_i1,1)) / &
        & (PFAD%vGrid%surfs(p_i2,1) - PFAD%vGrid%surfs(p_i1,1))

      ! Interpolate to get d log Beta / d nu.  We need this to Doppler-correct
      ! Beta even if dBdv is not associated.
      dBdNu = &
        & PFAD%dAbsDnu(t_i1  ,p_i1  ) * (1.0-t_fac) * (1.0-p_fac) + &
        & PFAD%dAbsDnu(t_i1+1,p_i1  ) * t_fac       * (1.0-p_fac) + &
        & PFAD%dAbsDnu(t_i1  ,p_i1+1) * (1.0-t_fac) * p_fac       + &
        & PFAD%dAbsDnu(t_i1+1,p_i1+1) * t_fac * p_fac

      ! Interpolate to get log Beta at the linearization velocity, then
      ! exponentiate to get Beta
      bp = ratio * ( 1.0 + doppler * dBdNu ) * exp( &
        & PFAD%absorption(t_i1  ,p_i1  ) * (1.0-t_fac) * (1.0-p_fac) + &
        & PFAD%absorption(t_i1+1,p_i1  ) * t_fac       * (1.0-p_fac) + &
        & PFAD%absorption(t_i1  ,p_i1+1) * (1.0-t_fac) * p_fac       + &
        & PFAD%absorption(t_i1+1,p_i1+1) * t_fac * p_fac )
      beta_path(j) = beta_path(j) + bp

      !{ \raggedright 
      !  $\frac{\partial \beta}{\partial T} \approx
      !   \left. \frac{\partial \beta}{\partial T} \right|_{\nu=\nu_l} +
      !   (\nu-\nu_l)
      !    \left. \frac{\partial^2  \beta}{\partial T \partial \nu}
      !    \right|_{\nu=\nu_l}$;
      !  $\frac{\partial \beta}{\partial T} = \frac\beta{T}
      !   \frac{\partial \log \beta}{\partial \log T}$;
      !  $\frac{\partial^2 \beta}{\partial T \partial \nu} =
      !   \beta \frac{\partial^2 \log \beta}{\partial T \partial \nu} =
      !   \frac\beta{T}
      !    \frac{\partial^2 \log \beta}{\partial \log T \partial \nu}$;
      !  $\frac{\partial \beta}{\partial T} \approx
      !   \beta \left[\frac1T
      !    \left.
      !     \frac{\partial \log \beta}{\partial \log T}
      !    \right|_{\nu=\nu_l} +
      !    (\nu-\nu_l)
      !    \left.
      !     \frac{\partial^2 \log \beta}{\partial T \partial \nu}
      !    \right|_{\nu=\nu_l}
      !   \right] $;
      !  $\frac{\partial \log \beta}{\partial \log T} \approx
      !   \frac{\nabla_T \log \beta}{\Delta \log T}$;
      !  $\frac{\partial^2 \log \beta}{\partial T \partial \nu} \approx
      !   \nabla_T \left(\frac{\partial \log \beta}{\partial \nu} \right)
      !    \frac1{T \Delta \log T}$;
      !  $\frac{\partial \beta}{\partial T} \approx
      !   \frac\beta{T \Delta \log T} \left[
      !    \nabla_T \log \left. \beta \right|_{\nu=\nu_l} +
      !     (\nu-\nu_l) \nabla_T \left(
      !      \left. \frac{\partial \log \beta}{\partial \nu}\right|_{\nu=\nu_l}
      !      \right)
      !   \right]$.\\
      !  $\nabla_T$ means ``Differences in $T$ coordinate, interpolated in
      !  $P$ coordinate.''

      ! Now the derivatives.  Remember that bp includes ratio, so we don't
      ! need to weight by ratio here.

      if ( temp_der ) then
        ! Interpolate to get d^2 log Beta / d log T d Nu, to Doppler-correct
        ! d log Beta d log T.
        ! Interpolate to get d log Beta / d log T, then Doppler correct
        ! and multiply by Beta / T to get dBeta / dT
        dBdT(j) = dBdT(j) + bp / (del_t * t_path(j)) * ( &
          & ( PFAD%absorption(t_i1+1,p_i1  ) - &
          &   PFAD%absorption(t_i1,  p_i1  ) ) * (1.0-p_fac) + &
          & ( PFAD%absorption(t_i1+1,p_i1+1) - &
          &   PFAD%absorption(t_i1,  p_i1+1) ) * p_fac + &
          & doppler * ( ( PFAD%dAbsDnu(t_i1+1,p_i1  ) - &
          &               PFAD%dAbsDnu(t_i1,  p_i1  ) ) * (1.0-p_fac) + &
          &             ( PFAD%dAbsDnu(t_i1+1,p_i1+1) - &
          &               PFAD%dAbsDnu(t_i1,  p_i1+1) ) * p_fac ) )
      end if

      ! Interpolate to get d log Beta / d*, then multiply by Beta to get
      ! dBeta / d*.  We can't Doppler correct these because we don't have
      ! the second partials with respect to d* dNu.
      if ( associated(dBdw) ) dBdw(j) = dBdw(j) + bp * ( &
        & PFAD%dAbsDwc(t_i1  ,p_i1  ) * (1.0-t_fac) * (1.0-p_fac) + &
        & PFAD%dAbsDwc(t_i1+1,p_i1  ) * t_fac       * (1.0-p_fac) + &
        & PFAD%dAbsDwc(t_i1  ,p_i1+1) * (1.0-t_fac) * p_fac       + &
        & PFAD%dAbsDwc(t_i1+1,p_i1+1) * t_fac * p_fac )

      if ( associated(dBdn) ) dBdn(j) = dBdn(j) + bp * ( &
        & PFAD%dAbsDnc(t_i1  ,p_i1  ) * (1.0-t_fac) * (1.0-p_fac) + &
        & PFAD%dAbsDnc(t_i1+1,p_i1  ) * t_fac       * (1.0-p_fac) + &
        & PFAD%dAbsDnc(t_i1  ,p_i1+1) * (1.0-t_fac) * p_fac       + &
        & PFAD%dAbsDnc(t_i1+1,p_i1+1) * t_fac * p_fac )

      if ( associated(dBdv) ) dBdv(j) = dBdv(j) + bp * dBdNu

    end do ! j

  end subroutine Create_Beta_Path_PFA

!     =====     Private Procedures     =================================

  ! ------------------------------------------------  Abs_CS_Cont  -----

  ! Compute the general continuum contribution
  pure &
  function Abs_CS_Cont ( Cont, Temperature, Pressure, Frequency ) &
    & result(Abs_CS_Cont_r)
  ! real(rp) function Abs_CS_Cont ( Cont, Temperature, Pressure, Frequency )

    use MLSKinds, only: RP

    real(rp), intent(in) :: Cont(:)     ! Continuum parameters
    real(rp), intent(in) :: Temperature ! in Kelvin
    real(rp), intent(in) :: Pressure    ! in mbar
    real(rp), intent(in) :: Frequency   ! in MegaHertz
    real(rp) :: Abs_CS_Cont_r

    Abs_CS_Cont_r = cont(1) * pressure * pressure * frequency * frequency * &
      & ( (300.0_rp / temperature)**cont(2) )

  end function Abs_CS_Cont

  ! ---------------------------------------------  Abs_CS_Cont_dT  -----

  ! Compute the general continuum contribution and its temperature derivative
  pure &
  subroutine Abs_CS_Cont_dT ( Cont, Temperature, Pressure, Frequency, &
    & Beta, dBeta_dT )

    use MLSKinds, only: RP

    real(rp), intent(in) :: Cont(:)     ! Continuum parameters
    real(rp), intent(in) :: Temperature ! in Kelvin
    real(rp), intent(in) :: Pressure    ! in mbar
    real(rp), intent(in) :: Frequency   ! in MegaHertz
    real(rp), intent(out) :: Beta, dBeta_dT ! km^{-1}, km^{-1}/K

    real(rp) :: Onedt ! 1/T

!{ Let $\theta = \frac{300}T$.  Then the general continuum contribution is
!  $\beta = c_1 p^2 \nu^2 \theta^{c_2}$.  Noticing that
!  $\frac{\partial \theta}{\partial T} = -\frac{\theta}T$, we have
!  $\frac{\partial \beta}{\partial T} = -\beta \frac{c_2}T$.

    onedt = 1.0 / temperature
    beta = cont(1) * pressure * pressure * frequency * frequency * &
      & ( (300.0_rp * onedt)**cont(2) )

    dBeta_dT = -beta * cont(2) * onedt

  end subroutine Abs_CS_Cont_dT

  ! --------------------------------------------  Abs_CS_H2O_Cont  -----

  ! Compute the general continuum contribution
  pure &
  function Abs_CS_H2O_Cont ( Cont, Temperature, Pressure, Frequency, Sps ) &
    & result(Beta)

    use MLSKinds, only: RP

    real(rp), intent(in) :: Cont(:)     ! Continuum parameters
    real(rp), intent(in) :: Temperature ! in Kelvin
    real(rp), intent(in) :: Pressure    ! in mbar
    real(rp), intent(in) :: Frequency   ! in MegaHertz
    real(rp), intent(in), optional :: SPS ! Mixing ratio
    real(rp) :: Beta                    ! km^{-1}
    real(rp) :: dBdf
    real(rp) :: Psq_Fsq                 ! pressure**2 * frequency**2
    real(rp) :: Temp_Ratio              ! log(300/T)

    !{ Ignoring vmr dependence,
    !  $\beta = c_1 p^2 \nu^2 \left( \frac{300}T \right)^{c_2}$.
    !
    !  With vmr dependence,
    !  $\beta = c_1 p^2 \nu^2 \left( \frac{300}T \right)^{c_2}
    !   + c_3 f p^2 \nu^2 \left( \frac{300}T \right)^{c_4}$.

    psq_fsq = pressure*pressure * frequency*frequency
    temp_ratio = log(300.0_rp / temperature)
    beta = cont(1) * psq_fsq * exp(cont(2) * temp_ratio)
    if ( present(sps) ) then
      dBdf = cont(3) * psq_fsq * exp(cont(4) * temp_ratio)
      beta = beta + dBdf * sps
    end if

  end function Abs_CS_H2O_Cont

  ! -----------------------------------------  Abs_CS_H2O_Cont_df  -----

  ! Compute the general continuum contribution
  pure &
  subroutine Abs_CS_H2O_Cont_df ( Cont, Temperature, Pressure, Frequency, &
    & Sps, Beta, dBeta_df )

    use MLSKinds, only: RP

    real(rp), intent(in) :: Cont(:)     ! Continuum parameters
    real(rp), intent(in) :: Temperature ! in Kelvin
    real(rp), intent(in) :: Pressure    ! in mbar
    real(rp), intent(in) :: Frequency   ! in MegaHertz
    real(rp), intent(in) :: SPS         ! Mixing ratio
    real(rp), intent(out) :: Beta       ! km^{-1}
    real(rp), intent(out) :: dBeta_df
    real(rp) :: Psq_Fsq                 ! pressure**2 * frequency**2
    real(rp) :: Temp_Ratio              ! log(300/T)

    !{ With vmr dependence,
    !  $\beta = c_1 p^2 \nu^2 \left( \frac{300}T \right)^{c_2}
    !   + c_3 f p^2 \nu^2 \left( \frac{300}T \right)^{c_4}$.
    !  $\frac{\partial \beta}{\partial f} =
    !   c_3 p^2 \nu^2 \left( \frac{300}T \right)^{c_4}$.

    psq_fsq = pressure*pressure * frequency*frequency
    temp_ratio = log(300.0_rp / temperature)
    dBeta_df = cont(3) * psq_fsq * exp(cont(4) * temp_ratio )
    beta = cont(1) * psq_fsq * exp(cont(2) * temp_ratio ) + dBeta_df * sps

  end subroutine Abs_CS_H2O_Cont_df

  ! -----------------------------------------  Abs_CS_H2O_Cont_dT  -----

  ! Compute the general continuum contribution and its temperature derivative
  pure &
  subroutine Abs_CS_H2O_Cont_dT ( Cont, Temperature, Pressure, Frequency, &
    & SPS, Beta, dBeta_dT, dBeta_df )

    use MLSKinds, only: RP

    real(rp), intent(in) :: Cont(:)     ! Continuum parameters
    real(rp), intent(in) :: Temperature ! in Kelvin
    real(rp), intent(in) :: Pressure    ! in mbar
    real(rp), intent(in) :: Frequency   ! in MegaHertz
    real(rp), intent(in), optional :: SPS ! Mixing ratio
    real(rp), intent(out) :: Beta, dBeta_dT ! km^{-1}, km^{-1}/K
    real(rp), intent(out), optional :: dBeta_df ! Derivative w.r.t SPS

    real(rp) :: dBdf                    ! Derivative of beta w.r.t SPS
    real(rp) :: Onedt                   ! 1/T
    real(rp) :: Psq_Fsq                 ! pressure**2 * frequency**2
    real(rp) :: Temp_Ratio              ! log(300/T)

!{ Let $\theta = \frac{300}T$.  Then the general continuum contribution
!  (ignoring vmr dependence) is
!  $\beta = \beta_0 = c_1 p^2 \nu^2 \theta^{c_2}$.  Noticing that
!  $\frac{\partial \theta}{\partial T} = -\frac{\theta}T$, we have
!  $\frac{\partial \beta}{\partial T} = -\beta \frac{c_2}T$ and
!  $\frac{\partial \beta}{\partial f} = 0$.
!
!  With vmr dependence,
!  $\beta = c_1 p^2 \nu^2 \theta^{c_2} + f c_3 p^2 \nu^2 \theta^{c_4}$.
!
!  In this case, $\frac{\partial \beta}{\partial T} = -\frac1T \left(
!  c_2 \beta_0 + f c_3 c_4 p^2 \nu^2 \theta^{c_4} \right)$ and
!  $\frac{\partial \beta}{\partial f} = c_3 p^2 \nu^2 \theta^{c_4}$.


    psq_fsq = pressure*pressure * frequency*frequency
    onedt = 1.0 / temperature
    temp_ratio = log(300.0 * onedt)
    beta = cont(1) * psq_fsq * exp(cont(2) * temp_ratio)
    dBeta_dT = cont(2) * beta
    if ( present(sps) ) then
      dBdf = cont(3) * psq_fsq * exp(cont(4) * temp_ratio)
      if ( present(dBeta_df) ) dBeta_df = dBdf
      dBdf = dBdf * sps
      beta = beta + dBdf
      dBeta_dT = dBeta_dT + cont(4) * dBdf
    end if
    dBeta_dT = - dBeta_dT * onedt

  end subroutine Abs_CS_H2O_Cont_dT

  ! ---------------------------------------------  Abs_CS_N2_Cont  -----

  ! Compute the N2 continuum contribution
  pure &
  function Abs_CS_N2_Cont ( Cont, Temperature, Pressure, Frequency ) &
    & result(Abs_CS_N2_Cont_r)
  ! real(rp) Function Abs_CS_N2_cont ( Cont, Temperature, Pressure, Frequency )

    use MLSKinds, only: RP

    real(rp), intent(in) :: Cont(:)     ! Continuum parameters
    real(rp), intent(in) :: Temperature ! in Kelvin
    real(rp), intent(in) :: Pressure    ! in mbar
    real(rp), intent(in) :: Frequency   ! in MegaHertz
    real(rp) :: Abs_CS_N2_Cont_r

    REAL(rp) :: Theta, Fsqr, FSXT

!{ Let $\theta = \frac{300}T$ and $f = p^2 \nu^2 \theta^{c_2}$. Then the N2
!  continuum contribution to $\beta$ is\\
!  $\beta = f ( c_1 e^{-c_3 \nu^2 \theta} +
!   c_4 e^{-c_5  \nu^2 \theta} (c_6^2 + \nu^2) )$.

    theta = 300.0_rp / temperature
    fsqr = frequency * frequency
    fsxt = fsqr * theta
    Abs_CS_N2_Cont_r = pressure * pressure * fsqr * (theta**cont(2)) * &
                   & ( cont(1) * exp(-cont(3) * fsxt) + &
                   &   cont(4) * exp(-cont(5) * fsxt) * &
                   & (cont(6)**2 + fsqr))

  end function Abs_CS_N2_Cont

  ! ------------------------------------------  Abs_CS_N2_Cont_dT  -----

  ! Compute the N2 continuum contribution and its temperature derivative
  pure &
  subroutine Abs_CS_N2_Cont_dT ( Cont, Temperature, Pressure, Frequency, &
    & Beta, dBeta_dT )

    use MLSKinds, only: RP

    real(rp), intent(in) :: Cont(:)     ! Continuum parameters
    real(rp), intent(in) :: Temperature ! in Kelvin
    real(rp), intent(in) :: Pressure    ! in mbar
    real(rp), intent(in) :: Frequency   ! in MegaHertz
    real(rp), intent(out) :: Beta, dBeta_dT ! km^{-1}, km^{-1}/K

    real(rp) :: E1, E2, F, Fsqr, FSXT, OneDT, Theta

!{ Let $\theta = \frac{300}T$ and $f = p^2 \nu^2 \theta^{c_2}$. Then the N2
!  continuum contribution to $\beta$ is\\
!  $\beta = f ( c_1 e^{-c_3 \nu^2 \theta} +
!   c_4 e^{-c_5  \nu^2 \theta} (c_6^2 + \nu^2) )$.  Noticing that
!  $\frac{\partial \theta}{\partial T} = -\frac{\theta}T$, we have\\
!  $\frac{\partial \beta}{\partial T} = \frac1T \left ( -\beta c_2 +
!   f \nu^2 \theta \left( c_3 c_1 e^{-c_3 \nu^2 \theta} +
!    c_5 c_4 (c_6^2 + \nu^2) e^{-c_5  \nu^2 \theta} \right) \right)$.

    onedt = 1.0 / temperature
    theta = 300.0_rp * onedt
    fsqr = frequency * frequency
    fsxt = fsqr * theta
    f = pressure * pressure * fsqr * (theta**cont(2))
    e1 = cont(1) * exp(-cont(3) * fsxt)
    e2 = cont(4) * exp(-cont(5) * fsxt) * (cont(6)**2 + fsqr)
    beta = f * ( e1 + e2 )

    dBeta_dT = onedt * ( f * fsxt * ( cont(3) * e1 + cont(5) * e2 ) &
      &                  - beta * cont(2) )

  end subroutine Abs_CS_N2_Cont_dT

  ! ---------------------------------------------  Abs_CS_O2_Cont  -----

  ! Compute the O2 continuum contribution
  pure &
  function Abs_CS_O2_Cont ( Cont, Temperature, Pressure, Frequency ) &
    & result(Abs_CS_O2_Cont_r)
  ! real(rp) Function ABS_CS_O2_CONT ( Cont, Temperature, Pressure, Frequency )

    use MLSKinds, only: RP

!{ Let $\theta = \frac{300}T$, $f = (c_3 p \exp({c_4} \theta))^2$ and
!  $D = \frac1{\nu^2 + f}$.
!  Then the O2 continuum contribution to beta
!  is $\beta = c_1 p^2 \nu^2 \exp({c_2} \theta) D$.

    real(rp), intent(in) :: Cont(:)     ! Continuum parameters
    real(rp), intent(in) :: Temperature ! in Kelvin
    real(rp), intent(in) :: Pressure    ! in mbar
    real(rp), intent(in) :: Frequency   ! in MegaHertz
    real(rp) :: Abs_CS_O2_Cont_r

    real(rp) :: Theta, Fsqr

    theta = log ( 300.0_rp / temperature )
    fsqr = frequency * frequency
    Abs_CS_O2_Cont_r = cont(1) * pressure * pressure * fsqr * exp(theta*cont(2)) &
                   & / (fsqr + (cont(3) * pressure * exp(theta*cont(4)) )**2 )

  end function Abs_CS_O2_Cont

  ! ------------------------------------------  Abs_CS_O2_Cont_dT  -----

  ! Compute the O2 continuum contribution
  pure &
  subroutine Abs_CS_O2_Cont_dT ( Cont, Temperature, Pressure, Frequency, &
      & Beta, dBeta_dT )

    use MLSKinds, only: RP

    real(rp), intent(in) :: Cont(:)     ! Continuum parameters
    real(rp), intent(in) :: Temperature ! in Kelvin
    real(rp), intent(in) :: Pressure    ! in mbar
    real(rp), intent(in) :: Frequency   ! in MegaHertz
    real(rp), intent(out) :: Beta, dBeta_dT ! km^{-1}, km^{-1}/K

    real(rp) :: D, F, Fsqr, Onedt, Theta

!{ Let $\theta = \frac{300}T$, $f = (c_3 p \theta^{c_4})^2$ and
!  $D = \frac1{\nu^2 + f}$.
!  Then the O2 continuum contribution to beta
!  is $\beta = c_1 p^2 \nu^2 \theta^{c_2} D$. Noticing that
!  $\frac{\partial \theta}{\partial T} = -\frac{\theta}T$, we have
!  $\frac{\partial \beta}{\partial T} = \frac{\beta}T
!   ( 2 c_4 f D - c_2)$.

    onedt = 1.0 / temperature
    theta = log( 300.0_rp * onedt )
    fsqr = frequency * frequency
    f = ( cont(3) * pressure * exp(cont(4)*theta) )**2
    d = 1.0 / (fsqr + f)
    beta = cont(1) * pressure * pressure * fsqr * exp(cont(2)*theta) * d

    dBeta_dT = beta * onedt * ( 2.0 * cont(4) * f * d - cont(2) )

  end subroutine Abs_CS_O2_Cont_dT

!-----------------------------------------------------------------------
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module GET_BETA_PATH_M

! $Log$
! Revision 2.120  2014/04/22 00:09:01  vsnyder
! Remove unused identifiers
!
! Revision 2.119  2014/04/04 19:43:10  vsnyder
! Move computation of QN from get_beta_path to o2_abs_cs.  This avoids the
! need for an array temp.
!
! Revision 2.118  2014/04/02 17:04:53  wgread
! added line dependent qn capability for o2-wgr
!
! Revision 2.117  2013/08/30 03:56:23  vsnyder
! Revise use of trace_begin and trace_end
!
! Revision 2.116  2013/07/26 22:19:05  vsnyder
! Fiddle with dump switches
!
! Revision 2.115  2013/06/12 02:25:18  vsnyder
! Cruft removal
!
! Revision 2.114  2013/05/22 00:08:43  vsnyder
! Repair integer-to-logical assign that NAG caught and Intel missed
!
! Revision 2.113  2013/05/17 22:52:32  vsnyder
! Revise dump switch stuff
!
! Revision 2.112  2011/11/22 00:40:40  vsnyder
! Add units km^{-1} in comments about beta
!
! Revision 2.111  2011/11/11 00:40:41  vsnyder
! Treat Extinction and MIFExtinction the same way
!
! Revision 2.110  2011/11/09 00:15:40  vsnyder
! Consolidate USE statements for slabs_sw_m
!
! Revision 2.109  2011/08/26 17:54:15  pwagner
! purehunt recovers optimized functionality of fwdmdls own hunt
!
! Revision 2.108  2011/08/25 18:17:46  vsnyder
! Remove unused use names
!
! Revision 2.107  2011/07/29 01:57:47  vsnyder
! Only IWC instead of IWC_A and IWC_S
!
! Revision 2.106  2011/05/09 17:56:27  pwagner
! Converted to using switchDetail
!
! Revision 2.105  2011/03/31 19:53:55  vsnyder
! Don't use :np for the first bound in sps_path.  sps_path is fine-path X
! sps.  np is the coarse path length.
!
! Revision 2.104  2011/03/04 03:42:48  vsnyder
! Associate dBdf with dBeta_df_path, not dBeta_dT_path!
!
! Revision 2.103  2011/03/02 02:06:59  vsnyder
! Cannonball polishing
!
! Revision 2.102  2011/02/12 03:57:40  vsnyder
! Add mixing-ratio dependence for H2O derivatives
!
! Revision 2.101  2010/02/09 21:04:18  vsnyder
! Specify extents for derivative section pointers
!
! Revision 2.100  2010/02/05 03:29:09  vsnyder
! Remove unused dummy arguments
!
! Revision 2.99  2010/02/04 23:10:38  vsnyder
! Remove USE for unreferenced names
!
! Revision 2.98  2010/01/23 01:20:15  vsnyder
! Allow betas that depend upon mixing ratio, get beta for Cloud_A and
! Cloud_S from the Mie tables, get kinds from MLSKinds instead of MLSCommon
!
! Revision 2.97  2009/06/23 18:26:11  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.96  2009/06/16 17:37:47  pwagner
! Changed api for dump, diff routines; now rely on options for most optional behavior
!
! Revision 2.95  2009/06/13 01:10:01  vsnyder
! Label Beta array, option to dump zeta
!
! Revision 2.94  2008/10/03 16:29:56  livesey
! Added EXTINCTIONV2
!
! Revision 2.93  2008/04/18 22:51:30  vsnyder
! Correct self-continuum
!
! Revision 2.92  2008/02/29 01:59:39  vsnyder
! Added a separate H2O continuum routine
!
! Revision 2.91  2007/05/23 22:44:24  vsnyder
! New slabs structure
!
! Revision 2.90  2006/12/04 21:17:28  vsnyder
! Reorganize FullForwardModel to use automatic arrays instead of allocating
! pointer arrays.  Requires testing for zero size instead of testing for
! associated in several subsidiary procedures.
!
! Revision 2.89  2006/07/19 22:30:17  vsnyder
! Cannonball polishing
!
! Revision 2.88  2006/06/29 01:44:42  vsnyder
! Add tracing
!
! Revision 2.87  2006/04/11 18:32:41  vsnyder
! Add more dumps
!
! Revision 2.86  2006/04/05 19:16:49  vsnyder
! Use the 'clean' switch in some dumps
!
! Revision 2.85  2006/02/23 01:00:17  vsnyder
! Give spectroscopy vector pointers correct upper bounds
!
! Revision 2.84  2006/02/10 21:50:33  vsnyder
! Move RHi conversion to correct place, spiff up some dumps
!
! Revision 2.83  2006/02/08 21:38:18  vsnyder
! Add relative humidity (RHi) calculation
!
! Revision 2.82  2006/02/08 01:02:22  vsnyder
! More stuff for spectroscopy derivatives
!
! Revision 2.81  2005/10/24 20:19:35  vsnyder
! Corrections for spectroscopy derivatives
!
! Revision 2.80  2005/07/06 02:17:21  vsnyder
! Revisions for spectral parameter derivatives
!
! Revision 2.79  2005/06/09 02:34:16  vsnyder
! Move stuff from l2pc_pfa_structures to slabs_sw_m
!
! Revision 2.78  2005/06/03 01:58:53  vsnyder
! New copyright notice, move Id to not_used_here to avoid cascades,
! Revise PFA data structures.
!
! Revision 2.77  2005/05/24 01:55:18  vsnyder
! Delete unused symbols
!
! Revision 2.76  2005/05/02 23:05:01  vsnyder
! New data structures for PFA
!
! Revision 2.75  2005/03/29 01:58:17  vsnyder
! Make stuff pure
!
! Revision 2.74  2005/03/26 01:26:29  vsnyder
! Make sure dBeta_dw etc get a value even if there are no lines
!
! Revision 2.73  2005/03/25 21:04:57  vsnyder
! Don't clobber continuum in Create_Beat if there are lines
!
! Revision 2.72  2005/03/03 02:07:10  vsnyder
! Move dumps from Create_Beta_Path... to Get_Beta_Path...
!
! Revision 2.71  2005/02/17 02:35:13  vsnyder
! Remove PFA stuff from Channels part of config
!
! Revision 2.70  2005/02/16 23:16:50  vsnyder
! Revise data structures for split-sideband PFA
!
! Revision 2.69  2004/12/13 20:47:52  vsnyder
! Use Slabs_0%UseYi field instead of MaxVal(Abs(...%yi))
!
! Revision 2.68  2004/11/04 03:42:09  vsnyder
! Provide for both LBL_Ratio and PFA_Ratio in beta_group
!
! Revision 2.67  2004/11/01 20:26:36  vsnyder
! Reorganization of representation for molecules and beta groups; PFA may be broken for now
!
! Revision 2.66  2004/10/06 21:21:21  vsnyder
! Change how dumps are done
!
! Revision 2.65  2004/09/04 01:50:31  vsnyder
! get_beta_path_m.f90
!
! Revision 2.64  2004/09/02 18:14:29  vsnyder
! Doppler correct temperature derivative in PFA
!
! Revision 2.63  2004/09/01 01:48:13  vsnyder
! Closing in on PFA
!
! Revision 2.62  2004/08/31 18:32:17  vsnyder
! Move initialization for temp_der into loop in create_beta_path
!
! Revision 2.61  2004/08/05 20:59:02  vsnyder
! Get rid of beta_group%n_elements
!
! Revision 2.60  2004/08/03 22:06:45  vsnyder
! Inching further toward PFA
!
! Revision 2.59  2004/07/08 21:00:23  vsnyder
! Inching toward PFA
!
! Revision 2.58  2004/04/19 21:03:29  vsnyder
! Remove unused stuff; respect tder_path_flags
!
! Revision 2.57  2004/04/17 00:37:00  vsnyder
! Analytic temperature derivatives
!
! Revision 2.56  2004/04/02 00:59:24  vsnyder
! Get catalog from slabs structure
!
! Revision 2.55  2004/03/27 03:35:27  vsnyder
! Add pointer to catalog in slabs_struct.  Use it so as not to need to drag
! line centers and line widths around.  Write slabs_lines and slabswint_lines
! to get sum of beta over all lines; put slabs_struct instead of its components
! in the calling sequence.
!
! Revision 2.54  2004/03/20 04:08:55  vsnyder
! Steps along the way to analytic temperature derivatives
!
! Revision 2.53  2004/03/20 01:15:23  jonathan
!  minor changes
!
! Revision 2.52  2004/03/19 04:07:31  vsnyder
! Fix some blunders re dNu for spectral derivatives
!
! Revision 2.51  2004/03/19 00:47:02  vsnyder
! Use line center instead of pressure-shifted line center in a few places
!
! Revision 2.50  2004/02/27 22:47:32  bill
! fixed bug in n2 continuum calc
!
! Revision 2.49  2004/01/26 21:57:00  vsnyder
! Improve TeXnicalities
!
! Revision 2.48  2003/12/07 19:46:10  jonathan
! update for use in 2D cloud FWM
!
! Revision 2.47  2003/08/20 21:12:39  bill
! fixed tanh1 bug associated with T-ders
!
! Revision 2.46  2003/07/15 17:50:30  vsnyder
! Callers need t_der_path_flags to be a pointer
!
! Revision 2.45  2003/07/14 22:45:09  vsnyder
! Scale BP, BM by isotope ratio in t_power computation
!
! Revision 2.44  2003/07/11 22:43:37  vsnyder
! Multiply the continuum-derived Beta by the isotope ratio
!
! Revision 2.43  2003/07/09 22:47:43  vsnyder
! Make separate branches for the polarized and nonpolarized cases, so
! we don't need to check "if (present(polarized))" inside the loop.  This
! is the inner loop for the full forward model.
!
! Revision 2.42  2003/07/07 19:53:51  vsnyder
! Move newly-public Create_Beta and Create_Beta_Path above the 'Private
! Procedures' comment
!
! Revision 2.41  2003/07/07 19:08:38  vsnyder
! Make Create_Beta and Create_Beta_Path public
!
! Revision 2.40  2003/07/07 16:47:01  pwagner
! Moved declaration of 3 function results inside body to appease NAG
!
! Revision 2.39  2003/07/04 02:47:50  vsnyder
! Move Create_Beta here, add Create_Beta_Path routine
!
! Revision 2.38  2003/06/27 22:09:19  vsnyder
! Move allocation of LineWidths out of loops; simplify exponent calculation
!
! Revision 2.37  2003/06/18 17:23:40  bill
! fixed NAG associated bug
!
! Revision 2.36  2003/06/18 14:44:53  bill
! added subsetting feature for T-ders
!
! Revision 2.35  2003/06/02 22:41:33  vsnyder
! Remove unused symbols
!
! Revision 2.34  2003/05/16 23:51:51  livesey
! Now uses molecules rather than spectags
!
! Revision 2.33  2003/05/15 03:28:52  vsnyder
! Moved some stuff up to FullForwardModel because Get_d_Deltau_pol_dT needs it
!
! Revision 2.32  2003/05/10 00:48:09  vsnyder
! Add TeXnicalities
!
! Revision 2.31  2003/05/09 20:07:07  vsnyder
! Correct kind parameter for H, specify intent for GL_slabs and Beta_group
!
! Revision 2.30  2003/05/09 19:24:38  vsnyder
! Cosmetic change
!
! Revision 2.29  2003/05/05 23:00:25  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.25.2.9  2003/03/24 21:50:44  jonathan
! remove unused 'pressure' in call to cloud_extinction
!
! Revision 2.25.2.8  2003/03/22 04:03:04  vsnyder
! Move Beta_Group_T and Dump_Beta_Group from get_beta_path to Get_Species_Data
!
! Revision 2.25.2.7  2003/03/12 21:35:44  vsnyder
! Add Dump_Beta_Group and generic Dump for it
!
! Revision 2.25.2.6  2003/03/01 03:16:15  vsnyder
! Use 'polarized' to specify size of quantum numbers array
!
! Revision 2.25.2.5  2003/02/27 23:19:23  vsnyder
! Add polarized stuff.  Remove unused z_path_c argument of get_beta_path.
! Remove unused p_path argument to get_beta_cloud.  Cosmetics.
!
! Revision 2.25.2.3  2003/02/14 00:21:42  jonathan
! add singl. scat. albedo W0, ph funct PHH
!
! Revision 2.25.2.2  2003/02/13 22:26:30  jonathan
! changes dimension for beta_path_cloud
!
! Revision 2.25.2.1  2003/02/13 17:34:27  bill
! uses new slabs_sw
!
! Revision 2.25  2003/02/11 00:48:18  jonathan
! changes made after adding get_beta_path_cloud
!
! Revision 2.24  2003/02/07 01:57:19  vsnyder
! Delete USE RHIFromH2O because it's not used
!
! Revision 2.23  2003/02/07 01:08:34  jonathan
! remove ICON option for compute super saturation
!
! Revision 2.22  2003/02/06 22:12:49  jonathan
! fix bug
!
! Revision 2.21  2003/02/06 00:20:16  jonathan
! Add in many stuff to deal with clouds CloudIce, iwc_path, etc
!
! Revision 2.20  2003/02/04 22:03:33  jonathan
! ICON now equal to 0 as default
!
! Revision 2.19  2003/02/04 21:46:27  jonathan
! add ICON options for super saturation and dry cases
!
! Revision 2.18  2003/02/03 22:56:58  vsnyder
! Add Get_bata_path_polarized
!
! Revision 2.17  2003/01/31 18:45:09  jonathan
! use cld_ext only if Incl_Cld is ture
!
! Revision 2.16  2003/01/31 17:53:48  jonathan
! change z_path to z_path_c in passing to get_beta_path
!
! Revision 2.15  2003/01/31 17:16:08  jonathan
! add Inc_Cld, and cld_ext
!
! Revision 2.14  2003/01/30 17:43:04  jonathan
! remove RHtoEV
!
! Revision 2.13  2003/01/30 00:17:42  jonathan
! add z_path to get_beta_path & use Paul's RHIFromH2O to compute VMR from RHi
!
! Revision 2.12  2003/01/14 21:49:33  jonathan
! option for saturation below 100mb
!
! Revision 2.11  2003/01/08 00:17:29  vsnyder
! Use "associated" instead of "present" to control optional computations
!
! Revision 2.10  2002/12/13 02:06:51  vsnyder
! Use a SLABS structure for the slabs quantities
!
! Revision 2.9  2002/10/08 17:08:03  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.8  2002/09/12 23:00:04  vsnyder
! Cosmetic changes, move USEs from module scope to procedure scope
!
! Revision 2.7  2001/12/23 23:30:42  zvi
! Fixing a bug in the dbeta_dt computations
!
! Revision 2.6  2001/12/14 23:43:15  zvi
! Modification for Grouping concept
!
! Revision 2.5  2001/11/15 01:22:01  zvi
! Remove Extiction debug
!
! Revision 2.4  2001/11/10 00:46:40  zvi
! Adding the EXTINCTION capabilitis
!
! Revision 2.3  2001/11/07 22:24:45  zvi
! Further modification for the t-power computations
!
! Revision 2.2  2001/11/07 21:13:48  livesey
! Fixed bug with log(0.0/0.0) for molecues with no lines
! or continua
!
! Revision 2.1  2001/10/16 15:07:18  zvi
! Continuum parameters are now part of Catalog
!
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model
!
! Revision 1.22.2.1  2001/09/10 10:02:32  zvi
! Cleanup..comp_path_entities_m.f90
!
! Revision 1.8  2001/03/09 02:26:11  vsnyder
! More work on deallocation
!
! Revision 1.7  2001/03/09 02:11:28  vsnyder
! Repair deallocating
!
! Revision 1.6  2001/03/05 21:37:20  zvi
! New filter format
!
! Revision 1.1 2001/02/01 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
