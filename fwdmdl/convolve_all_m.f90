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
Subroutine convolve_all (ich, ptg_press, atmospheric, n_sps, temp_der, &
           tan_press,ptg_angles,tan_temp,dx_dt,d2x_dxdt,band,center_angle, &
           fft_pts, i_raw, k_temp, k_atmos, k_spect_dw, k_spect_dn,    &
           k_spect_dnu, spect_atmos, no_tan_hts, k_info_count, i_star_all, &
           k_star_all, k_star_info, no_t, no_phi_t, no_phi_f, InDir, Aaap, &
           spectroscopic, Ier)
!
    Logical, intent(in) :: temp_der
!
    integer(i4), intent(in) :: no_t, n_sps, no_tan_hts, band, &
   &                           ich, fft_pts, no_phi_t
    integer(i4), intent(in) :: no_phi_f(*), spect_atmos(*)
!
    real(r8), intent(in) :: CENTER_ANGLE
    real(r8), intent(in) :: TAN_PRESS(*), PTG_ANGLES(*), TAN_TEMP(*)
    real(r8), intent(in) :: DX_DT(Nptg,*), D2X_DXDT(Nptg,*)
    real(r8), intent(in) :: I_RAW(*)

    Real(r4) :: k_temp(Nptg,mxco,mnp)
    Real(r4) :: k_atmos(Nptg,mxco,mnp,Nsps)
    Real(r4) :: k_spect_dw(Nptg,mxco,mnp,Nsps),  &
                k_spect_dn(Nptg,mxco,mnp,Nsps),  &
                k_spect_dnu(Nptg,mxco,mnp,Nsps)
!
    Character(LEN=*), intent(in) :: InDir, Aaap
!
    type(limb_press), intent(in) :: PTG_PRESS
    type(atmos_comp), intent(in) :: ATMOSPHERIC(*)
    type (spectro_param), intent(in) :: SPECTROSCOPIC(*)
!
! -----     Output Variables   ----------------------------------------
!
    integer(i4), intent(out) :: IER, K_INFO_COUNT
!
    real(r8), intent(out) :: I_STAR_ALL(*)
    real(r4), intent(out) :: K_STAR_ALL(20,mxco,mnp,Nptg)
    type(k_matrix_info), intent(out) :: k_star_info(20)
!
! -----     Local Variables     ----------------------------------------
!
    integer(i4) :: FFT_INDEX(2**fft_pts)
    integer(i4) :: N, I, IS, J, K, NF, NTR, PTG_I, SV_I, Spectag, ki, kc
!
    real(r8) :: FFT_PRESS(2**fft_pts)
    real(r8) :: FFT_ANGLES(2**fft_pts), RAD(2**fft_pts)
    real(r8) :: SC1(Nlvl), SC2(Nlvl)
    real(r8) :: TERM(Nlvl), PtP(Nlvl), Q, R
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
!
!  Get 'ntr' pressures associated with the fft_angles:
!
    Call get_pressures('a',ptg_angles,tan_temp,tan_press,no_tan_hts, &
                        fft_angles,fft_press,ntr,Ier)
    if (Ier /= 0) Return
!
! Make sure the fft_press array is MONOTONICALY increasing:
!
    is = 1
    do while (is < ntr.and.fft_press(is) >= fft_press(is+1))
      is = is + 1
    end do
!
    k = 1
    Rad(k) = Rad(is)
    fft_index(k) = is
    fft_press(k) = fft_press(is)
!
    do ptg_i = is+1, ntr
      q = fft_press(ptg_i)
      if (q > fft_press(k)) then
        k = k + 1
        fft_press(k) = q
        Rad(k) = Rad(ptg_i)
        fft_index(k) = ptg_i
      endif
    end do
!
    if (k == ntr) fft_index(1) = -2
!
! Interpolate the output values and store the radiances in: i_star_all
!
    j = ptg_press%no_lin_values
    Call Cspline_der(fft_press,PtP,Rad,i_star_all,term,k,j)
!
! Find out if user wants pointing derivatives
!
    if (ptg_press%der_calc(band)) then
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
      k_star_info(kc)%no_phi_basis = no_phi_t
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
            q = d2x_dxdt(ptg_i,sv_i)
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
          if (fft_index(1) > 0) then
            do ptg_i = 1, k
              Rad(ptg_i) = Rad(fft_index(ptg_i))
            end do
          endif
!
          Call Cspline(fft_press,PtP,Rad,Sc1,k,j)
          k_star_all(ki,sv_i,nf,1:j) = Sc1(1:j)
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
          if (fft_index(1) > 0) then
            do ptg_i = 1, k
              Rad(ptg_i) = Rad(fft_index(ptg_i))
            end do
          endif
!
          Call Cspline(fft_press,PtP,Rad,term,k,j)
!
! Transfer dx_dt from convolution grid onto the output grid
!
          Call Lintrp(tan_press,PtP,dx_dt(1:,sv_i),sc2,no_tan_hts,j)
!
          do ptg_i = 1, j
            r = sc2(ptg_i) * term(ptg_i)
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
          if (fft_index(1) > 0) then
            do ptg_i = 1, k
              Rad(ptg_i) = Rad(fft_index(ptg_i))
            end do
          endif
!
          Call Cspline(fft_press,PtP,Rad,term,k,j)
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
! *** DEBUG
!
      if (ich == 61) then
        sv_i = 9
        r = -1.666667
        write(*,910) sv_i,r
        do nf = 1, no_phi_t
          write(*,919) nf,char(92),j
          write(*,909) (k_star_all(ki,sv_i,nf,ptg_i),ptg_i=1,j)
        end do
      endif
!
 909  Format(6(1x,1pe12.5))
 919  Format('di_dtemp_phi_',i1,a1,i4.4,'n')
!
 910  Format(/,21x,'CONVOLVED k_temp dump',/,4x,     &
     &'Derivatrives of Radiance with respect to TEMP',i2.2,  &
     &' (Zeta=',f7.4,')',/)
!
! *** END DEBUG
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
        k_star_info(kc)%name = atmospheric(is)%name
        k_star_info(kc)%first_dim_index = ki
        k_star_info(kc)%no_phi_basis = no_phi_f(is)
!
! Derivatives needed continue to process
!
        do nf = 1, no_phi_f(is)
!
          do sv_i = 1, atmospheric(is)%no_lin_values
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
            if (fft_index(1) > 0) then
              do ptg_i = 1, k
                Rad(ptg_i) = Rad(fft_index(ptg_i))
              end do
            endif
!
! Interpolate onto the output grid, and store in k_star_all ..
!
            Call Lintrp(fft_press,PtP,Rad,Sc1,k,j)
            k_star_all(ki,sv_i,nf,1:j) = Sc1(1:j)
!
          end do
!
        end do
!
! *** DEBUG
!
        if (ich == 61 .and. atmospheric(is)%name == 'H2O') then
          sv_i = 9
          r = -1.666667
          write(*,911) sv_i,r
          do nf = 1, no_phi_f(is)
            write(*,918) nf,char(92),j
            write(*,908) (k_star_all(ki,sv_i,nf,ptg_i),ptg_i=1,j)
          end do
        endif
!
 908  Format(6(1x,1pe12.5))
 918  Format('di_dh2o_phi_',i1,a1,i4.4,'n')
!
 911  Format(/,21x,'CONVOLVED k_atmos_dump',/,4x,       &
     &'Derivatrives of Radiance with respect to H2O',i2.2,   &
     &' (Zeta=',f7.4,')',/)
!
! *** END DEBUG
!
      endif
!
    end do
!
! ****************** Spectroscopic derivatives ******************
!
    do is = 1, n_sps
!
      j = spect_atmos(is)
      if(j < 1) CYCLE
      if(.not.spectroscopic(j)%DER_CALC(band)) CYCLE
!
! Derivatives needed continue to process
!
      Spectag = atmospheric(is)%spectag
!
      DO
!
        if(spectroscopic(j)%Spectag /= Spectag) EXIT
        n = spectroscopic(j)%no_phi_values
        i = spectroscopic(j)%no_zeta_values
        ki = ki + 1
        kc = kc + 1
        k_star_info(kc)%name = spectroscopic(j)%NAME
        k_star_info(kc)%first_dim_index = ki
        k_star_info(kc)%no_phi_basis = n
        CA = spectroscopic(j)%type
!
        do nf = 1, n
!
          do sv_i = 1, i
!
            select case ( CA )
              case ( 'W' )
                Rad(1:no_tan_hts) = k_spect_dw(1:no_tan_hts,sv_i,nf,j)
              case ( 'N' )
                Rad(1:no_tan_hts) = k_spect_dn(1:no_tan_hts,sv_i,nf,j)
              case ( 'V' )
                Rad(1:no_tan_hts) = k_spect_dnu(1:no_tan_hts,sv_i,nf,j)
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
            if (fft_index(1) > 0) then
              do ptg_i = 1, k
                Rad(ptg_i) = Rad(fft_index(ptg_i))
              end do
            endif
!
! Interpolate onto the output grid, and store in k_star_all ..
!
            Call Lintrp(fft_press,PtP,Rad,Sc1,k,j)
            k_star_all(ki,sv_i,nf,1:j) = Sc1(1:j)
!
          end do        ! sv_i loop
!
        end do          ! nf loop
!
        j = j + 1
        if(j > 3 * n_sps) EXIT
!
      END DO
!
! *** DEBUG
!
      if (ich == 61 .and. atmospheric(is)%name == 'H2O') then
        sv_i = 9
        r = -1.666667
        write(*,912) sv_i,r
        do nf = 1, no_phi_f(is)
          write(*,917) nf,char(92),j
          write(*,907) (k_star_all(ki,sv_i,nf,ptg_i),ptg_i=1,j)
        end do
      endif
!
 907  Format(6(1x,1pe12.5))
 917  Format('di_dh2o_w_phi_',i1,a1,i4.4,'n')
!
 912  Format(/,21x,'CONVOLVED k_spectro_dump',/,4x,       &
     &'Derivatrives of Radiance with respect to H2O_W',i2.2,   &
     &' (Zeta=',f7.4,')',/)
!
! *** END DEBUG
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
! Revision 1.1  2000/06/21 21:56:14  zvi
! First version D.P.
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
