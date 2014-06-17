subroutine Calc_GeodAngle (mif, fmaf, flosf)

  USE oa_file_contents, ONLY: emls_oa
  USE Constants, ONLY: Deg2Rad, Rad2Deg, Pi

  implicit none

  integer, INTENT(IN) :: mif
  real*8, INTENT(IN) :: fmaf(3,3), flosf(3,3)

  REAL*8, parameter :: ellipse(0:1) = (/ 6378.137D0, 6356.7523141D0 /)
  real :: tngt_clat, clat_t, slat_t, clong_t, slong_t, clat, slat
  REAL*8 :: r_e, n_e(3), n_s(3), ABCDEF(3,3), n_fovd(3), n_pp(3), beta_r, &
       del_phi_t, cbeta_r, sbeta_r, mu_r, los_ascdsc, cmu_r, phi, phi_t_inst, &
       gamma, c_e, los_mat(3,3)
  CHARACTER (LEN=*), PARAMETER :: nl = char(10)

  tngt_clat = emls_oa%geoc_lat(mif) * Deg2Rad
  slat_t = SIN (tngt_clat)
  clat_t = COS (tngt_clat)
  slong_t = SIN (emls_oa%lon(mif) * Deg2Rad)
  clong_t = COS (emls_oa%lon(mif) * Deg2Rad)
  r_e = SQRT ((ellipse(0)*clat_t)**2 + (ellipse(1)*slat_t)**2)
  n_e = [ellipse(0)*clat_t*clong_t, ellipse(0)*clat_t*slong_t, &
       ellipse(1)*slat_t] / r_e
  n_s = emls_oa%sc_ECR(:,mif) / SQRT (SUM (emls_oa%sc_ECR(:,mif)**2))
  del_phi_t = ACOS (SUM (n_s * n_e)) * Rad2Deg
  ABCDEF = RESHAPE (emls_oa%ECRtoFOV(:,mif), (/ 3, 3 /))

  n_fovd = ABCDEF(3,:)
  slat = SIN (emls_oa%sc_geoclat(mif) * Deg2Rad)
  clat = COS (emls_oa%sc_geoclat(mif) * Deg2Rad)
  call cross_product (n_s, n_fovd, n_pp)
  n_pp = n_pp / SQRT (SUM (n_pp**2))
  beta_r = ACOS (n_pp(3)) * Rad2Deg
  sbeta_r = SIN (0.5*Pi - beta_r*Deg2Rad)
  cbeta_r = COS (0.5*Pi - beta_r*Deg2Rad)
  mu_r = ASIN (MAX (MIN (-slat*sbeta_r / (clat*cbeta_r), 1.0d0), -1.0d0)) * &
       Rad2Deg

  los_mat = TRANSPOSE (MATMUL ((MATMUL (TRANSPOSE (ABCDEF), fmaf)), &
       TRANSPOSE (flosf)))   ! matches IDL output!!!
  los_ascdsc = los_mat(3,3) !n_fovd(3)
  if (los_ascdsc > 0.0) mu_r = 180.0d0 - mu_r

if (mif == 9916) then
   print *, 'matmuls:'
   print *, 'los_ascdsc: ', los_ascdsc
   !print *, 'los_mat:'
   !print *, los_mat
   !print *, 'mu_r: ', mu_r
endif

  cmu_r = COS (mu_r * Deg2Rad)
  phi = ACOS (cmu_r * clat ) * Rad2Deg
  if (emls_oa%sc_geoclat(mif) >= 0.0) then
     phi = 180.0d0 - phi
  else
     phi = 180.0d0 + phi
  endif
  gamma = phi + del_phi_t
  c_e = SQRT ((ellipse(0)*ellipse(1))**2 / &
       ((ellipse(0) * COS (Deg2Rad * (beta_r - 90.0d0)))**2 &
       + (ellipse(1) * SIN (Deg2Rad * (beta_r - 90.0d0)))**2))
  phi_t_inst = ATAN2 (ellipse(0)**2 * SIN (Deg2Rad * gamma), &
       c_e**2 * COS (Deg2Rad*gamma)) * Rad2Deg
  if (phi_t_inst < 0.0d0) phi_t_inst = 360.0d0 + phi_t_inst

  ! Save the Geodetic Angle:

  emls_oa%geod_ang(mif) = phi_t_inst

  ! Save the OrbIncl:

  emls_oa%sc_orbincl = beta_r

if (mif == 9916) then
   print *, 'c_e: ', c_e
   print *, 'r_e: ', r_e
   print *, 'n_e: ', n_e
   print *, 'ECR: ', emls_oa%ECR(:,mif)
   print *, 'ABCDEF: '
   print *, ABCDEF
   print *, 'phi_t_inst: ', phi_t_inst
   print *, 'beta_r: ', beta_r
   print *, 'del_phi_t: ', del_phi_t
   print *, 'mu_r: ', mu_r
   print *, 'phi: ', phi
   print *, 'gamma: ', gamma
endif

end subroutine Calc_GeodAngle
