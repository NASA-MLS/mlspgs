  MODULE fov_convolve_v2_m
  USE Allocate_Deallocate, only: allocate_test, deallocate_test
  USE AntennaPatterns_m, only: AntennaPattern_T
  USE fov_convolve_m, only: drft1
  USE MLSNumerics, only: interpolatevalues
  USE MLSCommon, ONLY: I4, R8, rp
  IMPLICIT none
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName = &
       "$RCSfile$"
!---------------------------------------------------------------------------
  CONTAINS
! ============================================  fov_convolve_v2 =====
! This subprogram adds the effects of antenna smearing to the radiance.
!
  SUBROUTINE fov_convolve_v2(aaap,aaap_step,chi_in,rad_in,chi_out,rad_out, &
  & req,rsc,earth_frac)
! inputs
  REAL(rp), INTENT(in) :: aaap(:)  ! ft of antenna power pattern
  REAL(rp), INTENT(in) :: aaap_step! aperture step size in wavelengths
  REAL(rp), INTENT(in) :: chi_in(:)! inputted pointing angles radians
  REAL(rp), INTENT(in) :: rad_in(:)! inputted radiances
  REAL(rp), INTENT(in) :: chi_out(:)! outputted pointing angles radians
  REAL(rp), OPTIONAL, INTENT(in) :: req ! equivalent earth radius
  REAL(rp), OPTIONAL, INTENT(in) :: rsc ! spacecraft radius
  REAL(rp), OPTIONAL, INTENT(in) :: earth_frac ! fraction of earth in total
!                                   filled out pattern
! note req, rsc and earth_frac are non critical parameters and don't
! really need to be supplied externally. They are used to partition the
! full fft field between earth and space components.
!  REAL(rp), OPTIONAL, INTENT(in) :: surface  ! An angle that defines the Earth surface.
! outputs
  REAL(rp), INTENT(out) :: rad_out(:) ! outputted radiances
! Internal stuff
  INTEGER(i4) :: i, pwr, no_fft
  REAL(r8) :: r_eq, r_sc, e_frac, init_angle
  REAL(r8), POINTER :: p(:)
  REAL(r8), POINTER :: angles(:)
  REAL(r8), POINTER :: rad_fft(:)
! some clunky stuff
  INTEGER, PARAMETER :: MAXP=12, MAX2P=2**MAXP
  INTEGER(i4), SAVE :: INIT = 0, MS = 0
  REAL(r8), SAVE :: S(MAX2P)
  r_eq = 6371.0_rp
  r_sc = r_eq + 705.0_rp
  e_frac = 0.18
  IF (PRESENT(req)) r_eq = req
  IF (PRESENT(rsc)) r_sc = rsc
  IF (PRESENT(earth_frac)) e_frac = earth_frac / 2.0
! nullify stuff
  NULLIFY(p,angles,rad_fft)
! find size of stuff
  pwr = 12
  no_fft = 2**pwr
! load up the antenna pattern
  CALL allocate_test(p,no_fft,'p',modulename)
  p = 0.0_r8
  p(1:MIN(no_fft,SIZE(aaap))) = aaap(1:MIN(no_fft,SIZE(aaap)))
! p is really a complex number masquerading as a real one
! construct the angles
  CALL allocate_test(angles,no_fft,'angles',modulename)
  angles = (/(i,i=1,no_fft)/) / (no_fft * aaap_step)
  angles = angles - angles(no_fft/2 + 1)
  init_angle = ASIN((r_eq - e_frac*SQRT(r_sc**2-r_eq**2)/aaap_step)/r_sc)
! set up the radiance array
  CALL allocate_test(rad_fft,no_fft,'rad_fft',modulename)
  CALL interpolatevalues(chi_in-init_angle,rad_in,angles(no_fft/2:no_fft), &
  & rad_fft(no_fft/2:no_fft),METHOD='S',EXTRAPOLATE='C')
! mirror reflect this
  rad_fft(1:no_fft/2-1) = (/(rad_fft(no_fft-i),i = 1, no_fft/2 - 1)/)
! I don't know if this step is truly necessary but it rephases the radiances
! identically to the prototype code
  rad_fft = CSHIFT(rad_fft,-1)
! take fft of radiance array
  IF (init > 0 .and. init /= no_fft) ms=0
  CALL drft1(rad_fft,'a',pwr,ms,s)
! apply convolution theorem
  rad_fft(1:2) = rad_fft(1:2) * p(1:2)
  DO i = 3, no_fft - 1, 2
    rad_fft(i+1) = rad_fft(i) * p(i+1)
    rad_fft(i)   = rad_fft(i) * p(i)
  ENDDO
  CALL drft1(rad_fft,'s',pwr,ms,s)
! interpolate to output grid
  CALL interpolatevalues(angles(no_fft/2:no_fft),rad_fft(no_fft/2:no_fft), &
  & chi_out-init_angle,rad_out,METHOD='S',EXTRAPOLATE='C')
  CALL deallocate_test(p,'p',modulename)
  CALL deallocate_test(angles,'angles',modulename)
  CALL deallocate_test(rad_fft,'rad_fft',modulename)
  END SUBROUTINE fov_convolve_v2
  END MODULE fov_convolve_v2_m
