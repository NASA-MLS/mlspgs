!*******************  Bill's version ****************
MODULE fov_convolve_v2_m

  USE Allocate_Deallocate, only: allocate_test, deallocate_test
  USE AntennaPatterns_m, only: AntennaPattern_T
  USE fov_convolve_m, only: drft1
  USE MLSNumerics, ONLY: interpolatevalues, hunt
  USE MLSCommon, ONLY: I4, r4, R8, rp
  USE MLSMessageModule, only: MLSMessage, MLSMSG_Error

  IMPLICIT none
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName = &
       "$RCSfile$"
!---------------------------------------------------------------------------
  CONTAINS
! ============================================  convolve_all_v2 =====
! This subprogram adds the effects of antenna smearing to the radiance.
!
  SUBROUTINE fov_convolve_v2(AntennaPattern,chi_in,rad_in,chi_out,rad_out, &
           & req,rsc,earth_frac,surf_angle,di_dt,dx_dt,ddx_dxdt,dx_dt_out, &
           & drad_dt_out,di_df,drad_df_out)
!
! inputs
!
  Type(antennaPattern_T), intent(in) :: AntennaPattern
!
  REAL(rp), INTENT(in) :: chi_in(:)! inputted pointing angles radians
  REAL(rp), INTENT(in) :: rad_in(:)! inputted radiances
  REAL(rp), INTENT(in) :: chi_out(:)! outputted pointing angles radians
!
  REAL(rp), OPTIONAL, INTENT(in) :: req ! equivalent earth radius
  REAL(rp), OPTIONAL, INTENT(in) :: rsc ! spacecraft radius
  REAL(rp), OPTIONAL, INTENT(in) :: earth_frac ! fraction of earth in total
!                                   filled out pattern
! note req, rsc and earth_frac are non critical parameters and don't
! really need to be supplied externally. They are used to partition the
! full fft field between earth and space components.
! stuff for temperature derivatives
  REAL(rp), OPTIONAL, INTENT(in) :: surf_angle ! An angle that defines the
!                                   Earth surface.
  REAL(rp), OPTIONAL, INTENT(in) :: di_dt(:,:) ! derivative of radiance wrt
!                                   temperature on chi_in
  REAL(rp), OPTIONAL, INTENT(in) :: dx_dt(:,:) ! derivative of angle wrt
!                                   temperature on chi_in
  REAL(rp), OPTIONAL, INTENT(in) :: ddx_dxdt(:,:) ! 2nd derivative wrt angle
!                                   temperature on chi_in
  REAL(rp), OPTIONAL, INTENT(in) :: dx_dt_out(:,:) ! derivative of angle wrt
!                                   temperature on chi_out
  REAL(rp), OPTIONAL, INTENT(in) :: di_df(:,:) ! mixing ratio derivatives or
!                                   any parameter where a simple convolution
!                                   will suffice
!
! outputs
!
  REAL(rp), INTENT(out) :: rad_out(:) ! outputted radiances
  REAL(rp), OPTIONAL, INTENT(out) :: drad_dt_out(:,:) ! outputted radiance
!                                    derivatives wrt temperature.
  REAL(rp), OPTIONAL, INTENT(out) :: drad_df_out(:,:) ! outputted radiance
!                                    derivatives for inputted di_df.
! Internal stuff
!
  INTEGER(i4) :: i, j, k, ffth, n_coeffs, zero_out
  INTEGER(i4) :: peak_grad(1)
!
  REAL(r8) :: r_eq, r_sc, e_frac, init_angle, aaap_step
!
  REAL(r8), POINTER :: p(:)
  REAL(r8), POINTER :: dp(:)
  REAL(r8), POINTER :: drad_dt_temp(:)
  REAL(r8), POINTER :: angles(:)
  REAL(r8), POINTER :: rad_fft(:)
  REAL(r8), POINTER :: rad_fft1(:)
!
! some clunky stuff
!
  INTEGER, PARAMETER :: pwr=12, no_fft=2**pwr
  INTEGER(i4), SAVE :: INIT = 0, MS = 0
  REAL(r8), SAVE :: S(no_fft)
!
  r_eq = 6371.0_rp
  r_sc = r_eq + 705.0_rp
  e_frac = 0.185
  IF (PRESENT(req)) r_eq = req
  IF (PRESENT(rsc)) r_sc = rsc
  IF (PRESENT(earth_frac)) e_frac = earth_frac / 2.0
!
! nullify stuff
!
  NULLIFY(p,dp,drad_dt_temp,angles,rad_fft,rad_fft1)
!
! load up the antenna pattern
!
  CALL allocate_test(p,no_fft,'p',modulename)
  CALL allocate_test(dp,no_fft,'dp',modulename)
!
  p = 0.0_r8
  dp = 0.0_r8
  aaap_step = antennaPattern%lambda
  j = MIN(no_fft,SIZE(antennaPattern%aaap))
  p(1:j) = antennaPattern%aaap(1:j)
  dp(1:j) = antennaPattern%d1aap(1:j)
!
! p & dp are really a complex numbers masquerading as a real ones
!
! construct the angles
!
  ffth = no_fft / 2
  CALL allocate_test(angles,no_fft,'angles',modulename)
  angles = (/(i,i=1,no_fft)/) / (no_fft * aaap_step)
  angles = angles - angles(ffth + 1)
  init_angle = ASIN((r_eq - e_frac*SQRT(r_sc**2-r_eq**2)/aaap_step)/r_sc)
!
! set up the radiance array
!
  CALL allocate_test(rad_fft,no_fft,'rad_fft',modulename)
  CALL allocate_test(rad_fft1,no_fft,'rad_fft1',modulename)
  CALL interpolatevalues(chi_in-init_angle,rad_in,angles(ffth:no_fft), &
     & rad_fft(ffth:no_fft),METHOD='S',EXTRAPOLATE='C')
!
! mirror reflect this
!
  rad_fft(1:ffth-1) = (/(rad_fft(no_fft-i),i = 1, ffth-1)/)
!
! I don't know if this step is truly necessary but it rephases the radiances
! identically to the prototype code
!
  rad_fft = CSHIFT(rad_fft,-1)
!
! take fft of radiance array
!
  IF (init > 0 .and. init /= no_fft) ms=0
  CALL drft1(rad_fft,'a',pwr,ms,s)
  if(ms == -2) then
    init=0
    call MLSMessage ( MLSMSG_Error, ModuleName,"Error in drft1")
  endif
!
! apply convolution theorem
!
  rad_fft1(1:2) = rad_fft(1:2) * p(1:2)
  DO i = 3, no_fft - 1, 2
    rad_fft1(i+1) = rad_fft(i) * p(i+1)
    rad_fft1(i)   = rad_fft(i) * p(i)
  ENDDO
!
  CALL drft1(rad_fft1,'s',pwr,ms,s)
  if(ms == -2) then
    init=0
    call MLSMessage ( MLSMSG_Error, ModuleName,"Error in drft1")
  endif
!
! interpolate to output grid
!
  CALL interpolatevalues(angles(ffth:no_fft),rad_fft1(ffth:no_fft), &
     & chi_out-init_angle,rad_out,METHOD='S',EXTRAPOLATE='C')
!
! determine if temperature derivatives are desired
!
  IF (PRESENT(drad_dt_out)) THEN
!
! temperature derivatives calculation
! compute the antenna derivative function
!
    n_coeffs = SIZE(di_dt,dim=2)
!
! find the surface dimension
!
    CALL hunt(angles(ffth:no_fft),surf_angle-init_angle,zero_out)
!
    CALL allocate_test(drad_dt_temp,SIZE(chi_out),'drad_dt_temp',modulename)
!
! third term first (its fft is coefficient independent)
! apply convolution theorem
!
    rad_fft1(1:2) = 0.0_rp
    DO i = 3, no_fft-1, 2
      rad_fft1(i+1) = rad_fft(i) * dp(i+1)
      rad_fft1(i)   = rad_fft(i) * dp(i)
    ENDDO
!
    CALL drft1(rad_fft1,'s',pwr,ms,s)
    if(ms == -2) then
      init=0
      call MLSMessage ( MLSMSG_Error, ModuleName,"Error in drft1")
    endif
!
! interpolate to output grid
!
    CALL interpolatevalues(angles(ffth:no_fft),rad_fft1(ffth:no_fft), &
       & chi_out-init_angle,drad_dt_temp,METHOD='S',EXTRAPOLATE='C')
!
    DO i = 1, n_coeffs
!
! estimate the error compensation
!
      peak_grad = MAXLOC(ddx_dxdt(:,i))
      k = peak_grad(1)
!
! do i*ddx_dxdt piece
!
      CALL interpolatevalues(chi_in-init_angle,(rad_in-rad_in(k)) * &
      &    ddx_dxdt(:,i), angles(ffth+zero_out+1:no_fft),  &
      &    rad_fft(ffth+zero_out+1:no_fft), METHOD='S',EXTRAPOLATE='C')
!
! zero out the subsurface stuff
!
      rad_fft(ffth:ffth+zero_out) = 0.0_rp
!
! add in di_dt part
!
      CALL interpolatevalues(chi_in-init_angle, di_dt(:,i), &
         & angles(ffth:no_fft), rad_fft1(ffth:no_fft), METHOD='S', &
         & EXTRAPOLATE='C')
!
      rad_fft(ffth:no_fft) = rad_fft(ffth:no_fft) + rad_fft1(ffth:no_fft)
!
! resymetrize
!
      rad_fft(1:ffth-1) = (/(rad_fft(no_fft-i),i = 1, ffth-1)/)
!
! I don't know if this step is truly necessary but it rephases the radiances
! identically to the prototype code
!
      rad_fft = CSHIFT(rad_fft,-1)
!
! take fft of radiance array
!
      CALL drft1(rad_fft,'a',pwr,ms,s)
      if(ms == -2) then
        init=0
        call MLSMessage ( MLSMSG_Error, ModuleName,"Error in drft1")
      endif
!
! do the rad_in * dx_dt term
!
      CALL interpolatevalues(chi_in-init_angle,(rad_in-rad_in(k)) * &
       &   dx_dt(:,i), angles(ffth+zero_out+1:no_fft), &
       &   rad_fft1(ffth+zero_out+1:no_fft), METHOD='S', EXTRAPOLATE='C')
!
! zero out array below surf_angle
!
      rad_fft1(ffth:ffth+zero_out) = 0.0_rp
!
! resymetrize
!
      rad_fft1(1:ffth-1) = (/(rad_fft1(no_fft-i),i = 1, ffth-1)/)
!
! I don't know if this step is truly necessary but it rephases the radiances
! identically to the prototype code
!
      rad_fft1 = CSHIFT(rad_fft1,-1)
!
! take fft of radiance array
!
      CALL drft1(rad_fft1,'a',pwr,ms,s)
      if(ms == -2) then
        init=0
        call MLSMessage ( MLSMSG_Error, ModuleName,"Error in drft1")
      endif
!
! apply convolution theorem
!
      rad_fft(1:2) = rad_fft(1:2) * p(1:2)
      DO j = 3, no_fft-1, 2
        rad_fft(j+1) = rad_fft(j) * p(j+1) - rad_fft1(j) * dp(j+1)
        rad_fft(j)   = rad_fft(j) * p(j)   - rad_fft1(j) * dp(j)
      ENDDO
!
! interplolate to chi_out
!
      CALL drft1(rad_fft,'s',pwr,ms,s)
      if(ms == -2) then
        init=0
        call MLSMessage ( MLSMSG_Error, ModuleName,"Error in drft1")
      endif

      CALL interpolatevalues(angles(ffth:no_fft), &
         & rad_fft(ffth:no_fft), chi_out-init_angle, drad_dt_out(:,i), &
         & METHOD='S',EXTRAPOLATE='C')
!
! compute final result
!
      drad_dt_out(:,i) = drad_dt_out(:,i) + dx_dt_out(:,i)*drad_dt_temp
!
    ENDDO               ! On i = 1, n_coeffs
!
    CALL deallocate_test(drad_dt_temp,'drad_dt_temp',modulename)
!
  ENDIF
!
  IF (PRESENT(drad_df_out)) THEN
! nominally the mixing ratio derivatives but can be used for any
! quantity requiring a simple convolution.
    n_coeffs = SIZE(di_df,dim=2)
    DO i = 1, n_coeffs
      CALL interpolatevalues(chi_in-init_angle,di_df(:,i), &
      & angles(no_fft/2:no_fft), rad_fft(no_fft/2:no_fft),METHOD='S', &
      & EXTRAPOLATE='C')
! mirror reflect this
      rad_fft(1:no_fft/2-1) = (/(rad_fft(no_fft-i),i = 1, no_fft/2 - 1)/)
! I don't know if this step is truly necessary but it rephases the radiances
! identically to the prototype code
      rad_fft = CSHIFT(rad_fft,-1)
! take fft of radiance array
      CALL drft1(rad_fft,'a',pwr,ms,s)
! apply convolution theorem
      rad_fft1(1:2) = rad_fft(1:2) * p(1:2)
      DO j = 3, no_fft - 1, 2
        rad_fft1(j+1) = rad_fft(j) * p(j+1)
        rad_fft1(j)   = rad_fft(j) * p(j)
      ENDDO
      CALL drft1(rad_fft1,'s',pwr,ms,s)
! interpolate to output grid
      CALL interpolatevalues(angles(no_fft/2:no_fft), &
      & rad_fft1(no_fft/2:no_fft), chi_out-init_angle,drad_df_out(:,i), &
      & METHOD='S',EXTRAPOLATE='C')
    ENDDO
  ENDIF
!
  init = no_fft

  CALL deallocate_test(p,'p',modulename)
  CALL deallocate_test(dp,'dp',modulename)
  CALL deallocate_test(angles,'angles',modulename)
  CALL deallocate_test(rad_fft,'rad_fft',modulename)
  CALL deallocate_test(rad_fft1,'rad_fft1',modulename)
!
 END SUBROUTINE fov_convolve_v2

END MODULE fov_convolve_v2_m
