module CONVOLVE_ALL_M
  use MLSCommon, only: I4, R4, R8
  use L2PCDIM, only: Nlvl, Nsps, Nptg, MNP => max_no_phi
  use FOV_CONVOLVE_M, only: FOV_CONVOLVE
  use HYDROSTATIC_INTRP, only: GET_PRESSURES
  use L2PC_FILE_PARAMETERS, only: MXCO => max_no_elmnts_per_sv_component
  use L2PC_PFA_STRUCTURES, only: ATMOS_COMP, LIMB_PRESS, SPECTRO_PARAM, &
                                 K_MATRIX_INFO
  use D_LINTRP_M, only: LINTRP
  use D_CSPLINE_M, only: CSPLINE
  use DCSPLINE_DER_M, only: CSPLINE_DER
  implicit NONE
  private
  public :: CONVOLVE_ALL
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
     "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
     "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------
! This subroutine transfers the derivatives over from the internal
! convolution grid to the users specified points. This module uses
! cubic spline interpolation to do the job.
!
Subroutine convolve_all (ptg_press,atmospheric,n_sps,temp_der,tan_press,   &
           ptg_angles,tan_temp,dx_dt,d2x_dxdt,band,center_angle,fft_pts,   &
           i_raw, k_temp, k_atmos, k_spect_dw, k_spect_dn,k_spect_dnu,     &
           spect_atmos, no_tan_hts, k_info_count, i_star_all,k_star_all,   &
           k_star_info,no_t,no_phi_t,no_phi_f,InDir,Aaap,spectroscopic,    &
           t_z_basis, Ier)
!
    Logical, intent(IN) :: temp_der
!
    integer(i4), intent(IN) :: no_t, n_sps, no_tan_hts, band, &
   &                           fft_pts, no_phi_t
    integer(i4), intent(IN) :: no_phi_f(*), spect_atmos(*)
!
    real(r8), intent(IN) :: CENTER_ANGLE
    real(r8), intent(IN) :: I_RAW(*), T_Z_BASIS(*)
    real(r8), intent(IN) :: TAN_PRESS(*), PTG_ANGLES(*), TAN_TEMP(*)
    real(r8), intent(IN) :: DX_DT(Nptg,*), D2X_DXDT(Nptg,*)

    Real(r4) :: k_temp(Nptg,mxco,mnp)
    Real(r4) :: k_atmos(Nptg,mxco,mnp,Nsps)
    Real(r4) :: k_spect_dw(Nptg,mxco,mnp,Nsps),  &
                k_spect_dn(Nptg,mxco,mnp,Nsps),  &
                k_spect_dnu(Nptg,mxco,mnp,Nsps)
!
    Character(LEN=*), intent(IN) :: InDir, Aaap
!
    type(limb_press), intent(IN) :: PTG_PRESS
    type(atmos_comp), intent(IN) :: ATMOSPHERIC(*)
    type (spectro_param), intent(IN) :: SPECTROSCOPIC(*)
!
! -----     Output Variables   ----------------------------------------
!
    integer(i4), intent(out) :: IER, K_INFO_COUNT
!
    real(r8), intent(OUT) :: I_STAR_ALL(*)
    real(r4), intent(OUT) :: K_STAR_ALL(:,:,:,:)
    type(k_matrix_info), intent(OUT) :: k_star_info(*)
!
! -----     Local Variables     ----------------------------------------
!
    integer(i4) :: FFT_INDEX(2**fft_pts), nz
    integer(i4) :: n,i,is,j,k,nf,Ntr,ptg_i,sv_i,Spectag,ki,kc,jp,si
!
    real(r8) :: FFT_PRESS(2**fft_pts)
    real(r8) :: FFT_ANGLES(2**fft_pts), RAD(2**fft_pts)
    real(r8) :: SC1(Nlvl), TERM(Nlvl), PtP(Nlvl), Q, R
!
    Character(LEN=01) :: CA
!
! -----  Begin the code  -----------------------------------------
!
! Compute the ratio of the strengths
!
! This subroutine is called by channel
!
    Ier = 0
    ntr = 2**fft_pts
    K_INFO_COUNT = 0
    jp = (no_phi_t+1)/2
!
    Rad(1:no_tan_hts) = i_raw(1:no_tan_hts)
!
    kc = 0
    j = ptg_press%no_lin_values
    PtP(1:j) = dble(ptg_press%lin_val(1:j))
!
! Compute the convolution of the mixed radiances
!
    fft_angles(1:no_tan_hts) = ptg_angles(1:no_tan_hts)
    Call fov_convolve(fft_angles,Rad,center_angle,1,no_tan_hts,band, &
   &                  fft_pts,InDir,Aaap,Ier)
    if (Ier /= 0) Return

    si = no_tan_hts - j + 1
    Call Cspline(fft_angles,ptg_angles(si:no_tan_hts),Rad,Sc1,Ntr,j)
    i_star_all(1:j) = Sc1(1:j)
!
! Find out if user wants pointing derivatives
!
    if (ptg_press%der_calc(band)) then
!
!  Get 'Ntr' pressures associated with the fft_angles:
!
      Call get_pressures('a',ptg_angles,tan_temp,tan_press,no_tan_hts, &
                        fft_angles,fft_press,Ntr,Ier)
      if (Ier /= 0) Return
!
! Make sure the fft_press array is MONOTONICALY increasing:
!
      is = 1
      do while (is < Ntr.and.fft_press(is) >= fft_press(is+1))
        is = is + 1
      end do
!
      k = 1
      Rad(k) = Rad(is)
      fft_index(k) = is
      fft_press(k) = fft_press(is)
!
      do ptg_i = is+1, Ntr
        q = fft_press(ptg_i)
        if (q > fft_press(k)) then
          k = k + 1
          fft_press(k) = q
          Rad(k) = Rad(ptg_i)
          fft_index(k) = ptg_i
        endif
      end do
!
      if (k == Ntr) fft_index(1) = -2
!
! Interpolate the output values and store the radiances derivative
! with respect to pointing pressures in: term
!
      j = ptg_press%no_lin_values
      Call Cspline_der(fft_press,PtP,Rad,Sc1,term,k,j)
!
! Derivatives wanted,find index location k_star_all and write the derivative
!
      ki = 1
      kc = kc + 1
      k_star_info(kc)%name = 'PTAN'
      k_star_info(kc)%first_dim_index = ki
      k_star_info(kc)%no_phi_basis = mnp
      do nf = 1, mnp
        k_star_all(ki,1,nf,1:j) = term(1:j)
      end do
!
    endif
!
! Now transfer the other fwd_mdl derivatives to the output pointing
! values
!
! ********************* Temperature derivatives ******************
!
! check to determine if derivative is desired for this parameter
!
    if (temp_der) then
!
! Derivatives needed continue to process
!
      ki = 2
      kc = kc + 1
      k_star_info(kc)%name = 'TEMP'
      k_star_info(kc)%first_dim_index = ki
      k_star_info(kc)%no_zeta_basis = no_t
      k_star_info(kc)%no_phi_basis = no_phi_t
      k_star_info(kc)%zeta_basis(1:no_t) = t_z_basis(1:no_t)
!
      do nf = 1, no_phi_t
!
        do sv_i = 1, no_t
!
! run through representation basis coefficients
!
! Integrand over temperature derivative plus pointing differential
!
          do ptg_i = 1, no_tan_hts
            q = 0.0
            if(nf == jp) q = d2x_dxdt(ptg_i,sv_i)
            Rad(ptg_i) = i_raw(ptg_i) * q + k_temp(ptg_i,sv_i,nf)
          end do
!
! Now, Convolve:
!
          fft_angles(1:no_tan_hts) = ptg_angles(1:no_tan_hts)
          Call fov_convolve(fft_angles,Rad,center_angle,1,no_tan_hts, &
   &           band,fft_pts,InDir,Aaap,Ier)
          if (Ier /= 0) Return
!
          Call Cspline(fft_angles,ptg_angles(si:no_tan_hts),Rad,Sc1,Ntr,j)
          k_star_all(ki,sv_i,nf,1:j) = Sc1(1:j)
!
!  For any index off center Phi, skip the rest of the phi loop ...
!
          if(nf /= jp) CYCLE
!
! Now the convolution of radiance with the derivative antenna field
!
          Rad(1:no_tan_hts) = i_raw(1:no_tan_hts)
!
! Now, Convolve:
!
          fft_angles(1:no_tan_hts) = ptg_angles(1:no_tan_hts)
          Call fov_convolve(fft_angles,Rad,center_angle,2,no_tan_hts, &
   &           band,fft_pts,InDir,Aaap,Ier)
          if (Ier /= 0) Return
!
          Call Cspline(fft_angles,ptg_angles(si:no_tan_hts),Rad,term,Ntr,j)
!
! Transfer dx_dt from convolution grid onto the output grid
!
          Call Lintrp(tan_press,PtP,dx_dt(1:,sv_i),Sc1,no_tan_hts,j)
!
          do ptg_i = 1, j
            r = Sc1(ptg_i) * term(ptg_i)
            q = k_star_all(ki,sv_i,nf,ptg_i)
            k_star_all(ki,sv_i,nf,ptg_i) = q + r
          end do
!
! the convolution of the radiance weighted hydrostatic derivative
! with the antenna derivative
!
          Rad(1:no_tan_hts) = &
              dx_dt(1:no_tan_hts,sv_i) * i_raw(1:no_tan_hts)
!
! Now, convolve:
!
          fft_angles(1:no_tan_hts) = ptg_angles(1:no_tan_hts)
          Call fov_convolve(fft_angles,Rad,center_angle,2,no_tan_hts, &
   &           band,fft_pts,InDir,Aaap,Ier)
          if (Ier /= 0) Return
!
          Call Cspline(fft_angles,ptg_angles(si:no_tan_hts),Rad,term,Ntr,j)
!
          do ptg_i = 1, j
            q = k_star_all(ki,sv_i,nf,ptg_i)
            k_star_all(ki,sv_i,nf,ptg_i) = q - term(ptg_i)
          end do
!
        end do
!
      end do
!
    endif
!
! ****************** atmospheric derivatives ******************
!
    ki = 2
    do is = 1, n_sps
!
      if (atmospheric(is)%der_calc(band)) then
!
        ki = ki + 1
        kc = kc + 1
        nz = atmospheric(is)%no_lin_values
        k_star_info(kc)%name = atmospheric(is)%name
        k_star_info(kc)%first_dim_index = ki
        k_star_info(kc)%no_phi_basis = no_phi_f(is)
        k_star_info(kc)%no_zeta_basis = nz
        k_star_info(kc)%zeta_basis(1:nz) = &
                &  atmospheric(is)%basis_peaks(1:nz)
!
! Derivatives needed continue to process
!
        do nf = 1, no_phi_f(is)
!
          do sv_i = 1, nz
!
! run through representation basis coefficients
!
            Rad(1:no_tan_hts) = k_atmos(1:no_tan_hts,sv_i,nf,is)
!
! Now Convolve the derivative
!
            fft_angles(1:no_tan_hts) = ptg_angles(1:no_tan_hts)
            Call fov_convolve(fft_angles,Rad,center_angle,1,no_tan_hts, &
   &             band,fft_pts,InDir,Aaap,Ier)
            if (Ier /= 0) Return
!
! Interpolate onto the output grid, and store in k_star_all ..
!
            Call Lintrp(fft_angles,ptg_angles(si:no_tan_hts),Rad,Sc1,Ntr,j)
            k_star_all(ki,sv_i,nf,1:j) = Sc1(1:j)
!
          end do
!
        end do
!
      endif
!
    end do
!
! ****************** Spectroscopic derivatives ******************
!
    do is = 1, n_sps
!
      i = spect_atmos(is)
      if(i < 1) CYCLE
      if(.not.spectroscopic(i)%DER_CALC(band)) CYCLE
!
! Derivatives needed continue to process
!
      Spectag = atmospheric(is)%spectag
!
      DO
!
        if(spectroscopic(i)%Spectag /= Spectag) EXIT
        n = spectroscopic(i)%no_phi_values
        nz = spectroscopic(i)%no_zeta_values
        CA = spectroscopic(i)%type
        ki = ki + 1
        kc = kc + 1
        k_star_info(kc)%name = spectroscopic(i)%NAME
        k_star_info(kc)%first_dim_index = ki
        k_star_info(kc)%no_phi_basis = n
        k_star_info(kc)%no_zeta_basis = nz
        k_star_info(kc)%zeta_basis(1:nz) = &
                &  spectroscopic(i)%zeta_basis(1:nz)
!
        do nf = 1, n
!
          do sv_i = 1, nz
!
            select case ( CA )
              case ( 'W' )
                Rad(1:no_tan_hts) = k_spect_dw(1:no_tan_hts,sv_i,nf,i)
              case ( 'N' )
                Rad(1:no_tan_hts) = k_spect_dn(1:no_tan_hts,sv_i,nf,i)
              case ( 'V' )
                Rad(1:no_tan_hts) = k_spect_dnu(1:no_tan_hts,sv_i,nf,i)
              case default
                Ier = -99
                Print *,'** Unknown Spectroscopic element !'
                Return
            end select
!
! Now Convolve the derivative
!
            fft_angles(1:no_tan_hts) = ptg_angles(1:no_tan_hts)
            Call fov_convolve(fft_angles,Rad,center_angle,1,no_tan_hts, &
   &             band,fft_pts,InDir,Aaap,Ier)
            if (Ier /= 0) Return
!
! Interpolate onto the output grid, and store in k_star_all ..
!
            Call Lintrp(fft_angles,ptg_angles(si:no_tan_hts),Rad,Sc1,Ntr,j)
            k_star_all(ki,sv_i,nf,1:j) = Sc1(1:j)
!
          end do        ! sv_i loop
!
        end do          ! nf loop
!
        i = i + 1
        if(i > 3 * n_sps) EXIT
!
      END DO
!
    end do
!
    K_INFO_COUNT = kc
    Return
!
  End Subroutine CONVOLVE_ALL
!
end module CONVOLVE_ALL_M
! $Log$
! Revision 1.2  2001/02/19 22:14:21  zvi
!
! Revision 1.1  2000/06/21 21:56:14  zvi
! First version D.P.
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
