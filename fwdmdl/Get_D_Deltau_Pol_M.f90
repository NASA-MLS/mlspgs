! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Get_D_Deltau_Pol_M

  implicit NONE
  private
  public :: Get_D_Deltau_Pol_DT

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    &  "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains

!{Compute {\tt D\_Deltau\_Pol\_DT}.
!
! Assume $\beta$ can be approximated by $\hat\beta = \beta_0
! \left(\frac{T}{T_0}\right)^n$, but we don't know $n$.  So evaluate $\beta$
! for $T$, $T+$ some $\delta T$, and $T-$ some $\delta T$ (not necessarily the
! same $\delta T$).  Taking logarithms and differences, we get three estimates
! $n = \ln ( \beta(T+\delta T)/\beta(T-\delta T) ) /
!      \ln ( (T+\delta T) / (T-\delta T) )$,
! $n = \ln ( \beta(T+\delta T)/\beta(T) ) /
!      \ln ( (T+\delta T) / T)$, and
! $n = \ln ( \beta(T)/\beta(T-\delta T) ) /
!      \ln ( T / (T-\delta T))$.
! (The terms involving $\beta_0$ and $n \ln T_0$ cancel out.)
!
! Now observe that $\frac{\partial \hat\beta}{\partial T} = \frac{n}T \hat\beta$.
! Even better, use $\frac{n}T \beta$.  In any case, use a weighted average of
! the three values of $n$ computed above.  Also remember that
! $T = \sum_{i=1}^N \eta_i T_i$, where $N$ is the number of elements in the
! representation basis and $\eta_i$ are the coefficients, and what we really
! want is
! $\frac{\partial\beta}{\partial T_m} = 
!  \frac{n}T \beta \frac{\partial  T}{\partial T_m} =
!  \frac{n}T \beta \eta_m$, where $m$ is a state vector element index.
!
! Remembering that $\alpha = \sum_{i=1}^s f_i \beta_i$, and
! having $\frac{\partial \beta}{\partial T_m}$, compute
! $\frac{\partial \alpha}{\partial T_m} =
!  \sum_{i=1}^s f_i \frac{\partial\beta_i}{\partial T_m}$ ($s$ is the number
! of species).
!
! Then compute the derivative of incremental optical depth,
! $\frac{\partial \Delta \delta_{i \rightarrow i-1}^k}{\partial T_m}$ from
! $\frac{\partial \alpha}{\partial T_m}$ using the {\tt Opacity} routine.
!
! Finally, compute {\tt D\_Deltau\_Pol\_DT} = $\frac{\partial \bf E}{\partial
! T_m} = \frac\partial{\partial T} \exp ( -\Delta \delta_{i \rightarrow i-1}^k)$
! using $\Delta \delta_{i \rightarrow i-1}^k$,
! $\frac{\partial \Delta \delta_{i \rightarrow i-1}^k}{\partial T_m}$
! and the {\tt dExDt} routine.
! 

  subroutine Get_D_Deltau_Pol_DT ( Frq, H, CT, STCP, STSP, My_Catalog, &
                & Beta_group, GL_slabs_M, GL_slabs_P, &
                & T_Path, T_Path_M, T_Path_P, Tanh_Path, &
                & Beta_Path, SPS_Path, Eta_zxp, Del_S, Path_inds, &
                & Incoptdepth, D_Deltau_Pol_DT )

    use DExdT_m, only: dExdT
    use Get_Beta_Path_m, only: Get_Beta_Path_Polarized
    use Get_Species_Data_m, only: Beta_Group_T
    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
    use MLSCommon, only: R8, RP, IP
    use Opacity_m, only: Opacity
    use Physics, only: H_OVER_K
    use SpectroscopyCatalog_m, only: CATALOG_T

  ! Arguments
    ! SVE == # of state vector elements
    real(r8), intent(in) :: Frq          ! frequency in MHz
    real(rp), intent(in) :: H(:)         ! Magnetic field component in instrument
                                         ! polarization on the path
    real(rp), intent(in) :: CT(:)        ! Cos(Theta), where theta
      ! is the angle between the line of sight and magnetic field vectors.
    real(rp), intent(in) :: STCP(:)      ! Sin(Theta) Cos(Phi) where
      ! theta is as for CT and phi (for this purpose only) is the angle
      ! between the plane defined by the line of sight and the magnetic
      ! field vector, and the "instrument field of view plane polarized"
      ! (IFOVPP) X axis.
    real(rp), intent(in) :: STSP(:)      ! Sin(Theta) Sin(Phi)
    type(catalog_t), intent(in) :: My_Catalog(:)
    type (beta_group_T), intent(in) :: Beta_group(:)
    type (slabs_struct), intent(in) :: GL_slabs_m(:,:), GL_slabs_p(:,:) ! for T -/+ del_T
    real(rp), intent(in) :: T_Path(:)    ! path temperatures
    real(rp), intent(in) :: T_Path_M(:), T_Path_P(:) ! path temperatures -/+ del_temp
    real(rp), intent(in) :: Tanh_Path(:) ! tanh(h \nu / 2 k T_Path)
    complex(rp), intent(in) :: Beta_Path(-1:,:,:) ! -1:1 x path x sps
    real(rp), intent(in) :: SPS_Path(:,:) ! species on path, path x sps
    real(rp), intent(in) :: Eta_zxp(:,:) ! representation basis function
      !                                    main grid.  path x sve
    real(rp), intent(in) :: Del_S(:)     !  ! unrefracted path length.
    integer(ip), intent(in) :: Path_inds(:) ! indicies for reading gl_slabs
    complex(rp), intent(in) :: Incoptdepth(:,:,:) ! 2 x 2 x path

    complex(rp), intent(out) :: D_Deltau_Pol_DT(:,:,:,:) ! 2 x 2 x path x sve

  ! Local variables
    complex(rp) :: Alpha_Path_N(-1:1)    ! alpha_path_n * N
    complex(rp) :: Beta(-1:1), Beta_M(-1:1), Beta_P(-1:1) ! Single elements of
      ! Beta_Path, Beta_Path_M, Beta_Path_P multiplied by Tanh_Path,  Tanh_M, Tanh_P
    complex(rp), dimension(-1:1,size(beta_path,2),size(beta_path,3)) :: Beta_Path_M, Beta_Path_P
    complex(rp):: D_Alpha_DT(-1:1,size(path_inds,1)) ! n/T * Alpha on the path
    complex(rp) :: D_Incoptdepth_dT(2,2,size(path_inds,1))
    real(r8) :: FrqHK                    ! 0.5 * Frq * H_Over_K
    integer :: I
    integer :: J, K
    real(rp) :: L_TTM, L_TPTM, L_TPT     ! Logarithms of temperature ratios
    complex(rp) :: N(-1:1)               ! Exponent of (T/T_0) in
    ! approximation to beta.  One each for Sigma_-, Pi and Sigma_+.
    integer :: N_Path, N_Sps
    real(rp) :: Tanh_M, Tanh_P           ! for T -/+ del_T

    frqhk = 0.5_r8 * Frq * H_Over_K
    n_path = size(path_inds)
    n_sps = size(sps_path,2)

    call get_beta_path_polarized ( frq, h, my_Catalog, beta_group, gl_slabs_m, &
      & path_inds, beta_path_m )
    call get_beta_path_polarized ( frq, h, my_Catalog, beta_group, gl_slabs_p, &
      & path_inds, beta_path_p )

    do i = 1, n_path
      alpha_path_n = (0.0_rp,0.0_rp)
      k = path_inds(i)

      l_ttm = log(t_path(k)/t_path_m(k))
      l_tptm = log(t_path_p(k)/t_path_m(k))
      l_tpt = log(t_path_p(k)/t_path(k))

      tanh_m = tanh( frqhk / t_path_m(k) )
      tanh_p = tanh( frqhk / t_path_p(k) )

      do j = 1, n_sps
        ! Solve for n
        beta = beta_path(:,k,j) * tanh_path(k)
        beta_m = beta_path_m(:,k,j) * tanh_m
        beta_p = beta_path_p(:,k,j) * tanh_p
        n = 0.25 * (      log(beta/beta_m)   / l_ttm +    &
          &         2.0 * log(beta_p/beta_m) / l_tptm +   &
          &               log(beta_p/beta)   / l_tpt )

        ! not quite D alpha, because we haven't multiplied by tanh or
        ! divided by T.
        alpha_path_n = alpha_path_n + n * beta_path(:,k,j) * sps_path(k,j)
      end do ! j
      ! Still not quite D alpha...
      alpha_path_n = alpha_path_n * tanh_path(k)

      ! Now it's more than D alpha, because we've multiplied by del_s, but
      ! this is OK, because OPACITY is linear.
      d_alpha_dT(:,i) = alpha_path_n * del_s(k) / t_path(k)

    end do ! i

    call opacity ( ct, stcp, stsp, d_alpha_dT, d_incoptdepth_dT )

!??? Do we need to have d_incoptdepth_dT for the nonpolarized case, and  ???
!??? subtract half of it from the diagonal of our d_incoptdepth_dT       ???

    do i = 1, n_path            ! along the path
      k = path_inds(i)
      do j = 1, size(eta_zxp,2) ! state vector elements
        if ( eta_zxp(k,j) /= 0.0 ) then
          call dExdT ( - incoptdepth(:,:,k) * eta_zxp(k,j), &
                     & - d_incoptdepth_dT(:,:,i) * eta_zxp(k,j), &
                     & d_deltau_pol_dT(:,:,i,j) ) ! d exp(incoptdepth) / dT
        else
          d_deltau_pol_dT(:,:,i,j) = 0.0_rp
        end if
      end do ! j
    end do ! i

  end subroutine Get_D_Deltau_Pol_DT

!-----------------------------------------------------------------------
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Get_D_Deltau_Pol_M

! $Log$
! Revision 2.1  2003/05/15 03:27:56  vsnyder
! Initial commit
!
