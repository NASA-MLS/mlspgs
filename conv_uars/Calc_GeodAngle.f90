! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Calc_GeodAngle_m

  implicit NONE
  private

  public :: Calc_GeodAngle

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains ! ============= Public procedures ===================================


  subroutine Calc_GeodAngle ( mif, fmaf, flosf, emls_oa )

    use OA_File_contents, only: emls_oa_t
    use Constants, only: Deg2Rad, Rad2Deg, Pi
    use Cross_m, only: Cross
    use MLSKinds, only: R8

    implicit none

    ! Args
    integer, intent(in) :: MIF
    real(r8), intent(in) :: FMAF(3,3), FLOSF(3,3)
    type(emls_oa_t), intent(inout) :: EMLS_OA

    real(r8), parameter :: ellipse(0:1) = (/ 6378.137D0, 6356.7523141D0 /)
    real :: tngt_clat, clat_t, slat_t, clong_t, slong_t, clat, slat
    real(r8) :: r_e, n_e(3), n_s(3), ABCDEF(3,3), n_fovd(3), n_pp(3), beta_r, &
         del_phi_t, cbeta_r, sbeta_r, mu_r, los_ascdsc, cmu_r, phi, phi_t_inst, &
         gamma, c_e, los_mat(3,3)
    character (len=*), parameter :: nl = char(10)

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
    n_pp = cross ( n_s, n_fovd, norm=.true. )
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
  end if

    cmu_r = COS (mu_r * Deg2Rad)
    phi = ACOS (cmu_r * clat ) * Rad2Deg
    if (emls_oa%sc_geoclat(mif) >= 0.0) then
       phi = 180.0d0 - phi
    else
       phi = 180.0d0 + phi
    end if
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
  end if

  end subroutine Calc_GeodAngle

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Calc_GeodAngle_m

! $Log$
