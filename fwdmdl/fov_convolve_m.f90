! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Fov_Convolve_m
 
  implicit NONE
  private
  public :: Fov_Convolve, Fov_Convolve_Old

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character(len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = &
    & "$RCSfile$"
!---------------------------------------------------------------------------
 contains
! =================================================  Fov_Convolve  =====
! This subprogram adds the effects of antenna smearing to the radiance.

    subroutine Fov_Convolve ( AntennaPattern, chi_in, rad_in, chi_out, rad_out, &
             & req, rsc, earth_frac, surf_angle, di_dt, dx_dt, ddx_dxdt,        &
             & dx_dt_out, drad_dt_out, di_df, di_df_flag, drad_df_out,          &
             & drad_dx_out )

    use AntennaPatterns_m, only: AntennaPattern_T
    use D_CSPLINE_M, only: CSPLINE
    use DFFT_M, only: DRFT1
    use MLSNumerics, only: interpolateValues, hunt
    use MLSCommon, only: I4, r4, R8, rp
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error

! inputs

    type(antennapattern_t), intent(in) :: AntennaPattern

    real(rp), intent(in) :: chi_in(:)! inputted pointing angles radians
    real(rp), intent(in) :: rad_in(:)! inputted radiances
    real(rp), intent(in) :: chi_out(:)! outputted pointing angles radians

    real(rp), optional, intent(in) :: req ! equivalent earth radius
    real(rp), optional, intent(in) :: rsc ! spacecraft radius
    real(rp), optional, intent(in) :: earth_frac ! fraction of earth in total
!                                   filled out pattern
! note req, rsc and earth_frac are non critical parameters and don't
! really need to be supplied externally. They are used to partition the
! full fft field between earth and space components.
! stuff for temperature derivatives

    real(rp), optional, intent(in) :: surf_angle ! An angle that defines the
!                                     Earth surface.
    real(rp), optional, intent(in) :: di_dt(:,:) ! derivative of radiance wrt
!                                     temperature on chi_in
    real(rp), optional, intent(in) :: dx_dt(:,:) ! derivative of angle wrt
!                                     temperature on chi_in
    real(rp), optional, intent(in) :: ddx_dxdt(:,:) ! 2nd derivative wrt angle
!                                     temperature on chi_in
    real(rp), optional, intent(in) :: dx_dt_out(:,:) ! derivative of angle wrt
!                                     temperature on chi_out
    real(rp), optional, intent(in) :: di_df(:,:) ! mixing ratio derivatives or
!                                     any parameter where a simple convolution
!                                     will suffice
    logical, optional, intent(in) :: di_df_flag(:) ! Flag to indicate which of the 
!                                     above are to be calculated

! outputs

    real(rp), intent(out) :: rad_out(:) ! outputted radiances
    real(rp), optional, intent(out) :: drad_dx_out(:) ! outputted derivative 
!                                      of radiance wrt to Chi_out
    real(rp), optional, intent(out) :: drad_dt_out(:,:) ! outputted radiance
!                                      derivatives wrt temperature.
    real(rp), optional, intent(out) :: drad_df_out(:,:) ! outputted radiance
!                                      derivatives for inputted di_df.
! Internal stuff

    integer(i4) :: i, j, k, ffth, n_coeffs, zero_out

    real(r8) :: r_eq, r_sc, e_frac, init_angle, aaap_step, ang_step

! some clunky stuff

    integer, parameter :: pwr=12, no_fft=2**pwr
    integer(i4), save :: INIT = 0, MS = 0
    real(r8), save :: S(no_fft)

    real(r8), dimension(no_fft) :: p, dp, angles, rad_fft, rad_fft1


    real(r8) :: drad_dt_temp(SIZE(chi_out))

    r_eq = 6371.0_rp
    r_sc = r_eq + 705.0_rp
    e_frac = 0.185
    if ( present(req)) r_eq = req
    if ( present(rsc)) r_sc = rsc
    if ( present(earth_frac)) e_frac = earth_frac / 2.0

! load up the antenna pattern

    p = 0.0_r8
    dp = 0.0_r8
    aaap_step = antennaPattern%lambda
    ang_step = 1.0_r8 / (no_fft * aaap_step)
    j = MIN(no_fft,SIZE(antennaPattern%aaap))
    p(1:j) = antennaPattern%aaap(1:j)
    dp(1:j) = antennaPattern%d1aap(1:j)

! p & dp are really complex numbers masquerading as real ones

! construct the angles

    ffth = no_fft / 2
    angles = (/(i*ang_step,i=1,no_fft)/)
    angles = angles - angles(ffth+1)
    init_angle = ASIN((r_eq - e_frac*SQRT(r_sc**2-r_eq**2)/aaap_step)/r_sc)

! set up the radiance array

    call interpolateValues ( chi_in-init_angle, rad_in, angles(ffth:no_fft),  &
       & rad_fft(ffth:no_fft), METHOD='S', EXTRAPOLATE='C' )

! mirror reflect this

    rad_fft(1:ffth-1) = (/(rad_fft(no_fft-i),i = 1, ffth-1)/)

! I don't know if this step is truly necessary but it rephases the radiances
! identically to the prototype code

    rad_fft = cshift(rad_fft,-1)

! take fft of radiance array

    if ( init > 0 .and. init /= no_fft ) ms=0
    call drft1_t ( rad_fft, 'a' )

! apply convolution theorem

    rad_fft1(1:2) = rad_fft(1:2) * p(1:2)
    do i = 3, no_fft - 1, 2
      rad_fft1(i+1) = rad_fft(i) * p(i+1)
      rad_fft1(i)   = rad_fft(i) * p(i)
    end do

    call drft1_t ( rad_fft1, 's' )

! interpolate to output grid

    if ( present(drad_dx_out) ) then
      call interpolateValues ( angles(ffth:no_fft)-ang_step, &
         & rad_fft1(ffth:no_fft), chi_out-init_angle, rad_out, &
         & METHOD='S', EXTRAPOLATE='C', dyByDx=drad_dx_out )
    else
      call interpolateValues ( angles(ffth:no_fft)-ang_step, &
         & rad_fft1(ffth:no_fft), chi_out-init_angle, rad_out, &
         & METHOD='S', EXTRAPOLATE='C' )
    end if

! determine if temperature derivatives are desired

    if ( present(drad_dt_out) ) then

! temperature derivatives calculation
! compute the antenna derivative function

      n_coeffs = SIZE(di_dt,dim=2)

! find the surface dimension

      call hunt ( angles(ffth:no_fft), surf_angle-init_angle, zero_out )

! third term first (its fft is coefficient independent)
! apply convolution theorem

      rad_fft1(1:2) = 0.0_rp
      do i = 3, no_fft-1, 2
        rad_fft1(i+1) = rad_fft(i) * dp(i+1)
        rad_fft1(i)   = rad_fft(i) * dp(i)
      end do

      call drft1_t ( rad_fft1, 's' )

! interpolate to output grid

      call interpolateValues ( angles(ffth:no_fft)-ang_step, &
         & rad_fft1(ffth:no_fft), chi_out-init_angle, drad_dt_temp, &
         & METHOD='S', EXTRAPOLATE='C' )

      do i = 1, n_coeffs

! estimate the error compensation

        k = maxloc(ddx_dxdt(:,i),1)

! do i*ddx_dxdt piece

        call interpolateValues ( chi_in-init_angle, (rad_in-rad_in(k)) * &
        &    ddx_dxdt(:,i), angles(ffth+zero_out+1:no_fft), &
        &    rad_fft(ffth+zero_out+1:no_fft), METHOD='S', EXTRAPOLATE='C' )

! zero out the subsurface stuff

        rad_fft(ffth:ffth+zero_out) = 0.0_rp

! add in di_dt part

        call interpolateValues ( chi_in-init_angle, di_dt(:, i), &
           & angles(ffth:no_fft), rad_fft1(ffth:no_fft), METHOD='S', &
           & EXTRAPOLATE='C' )

        rad_fft(ffth:no_fft) = rad_fft(ffth:no_fft) + rad_fft1(ffth:no_fft)

! resymetrize

        rad_fft(1:ffth-1) = (/(rad_fft(no_fft-i),i = 1, ffth-1)/)

! I don't know if this step is truly necessary but it rephases the radiances
! identically to the prototype code

        rad_fft = cshift(rad_fft,-1)

! take fft of radiance array

        call drft1_t ( rad_fft, 'a' )

! do the rad_in * dx_dt term

        call interpolateValues ( chi_in-init_angle, (rad_in-rad_in(k)) * &
         &   dx_dt(:,i), angles(ffth+zero_out+1:no_fft), &
         &   rad_fft1(ffth+zero_out+1:no_fft), METHOD='S', EXTRAPOLATE='C' )

! zero out array below surf_angle

        rad_fft1(ffth:ffth+zero_out) = 0.0_rp

! resymetrize

        rad_fft1(1:ffth-1) = (/(rad_fft1(no_fft-i),i=1,ffth-1)/)

! I don't know if this step is truly necessary but it rephases the radiances
! identically to the prototype code

        rad_fft1 = cshift(rad_fft1,-1)

! take fft of radiance array

        call drft1_t ( rad_fft1, 'a' )

! apply convolution theorem

        rad_fft(1:2) = rad_fft(1:2) * p(1:2)
        do j = 3, no_fft-1, 2
          rad_fft(j+1) = rad_fft(j) * p(j+1) - rad_fft1(j) * dp(j+1)
          rad_fft(j)   = rad_fft(j) * p(j)   - rad_fft1(j) * dp(j)
        end do

! interplolate to chi_out

        call drft1_t ( rad_fft, 's' )

        call interpolateValues ( angles(ffth:no_fft)-ang_step, &
           & rad_fft(ffth:no_fft), chi_out-init_angle, drad_dt_out(:, i), &
           & METHOD='S', EXTRAPOLATE='C' )

! compute final result

        drad_dt_out(:,i) = drad_dt_out(:,i) + dx_dt_out(:,i)*drad_dt_temp

      end do               ! On i = 1, n_coeffs

    end if

    if ( present(drad_df_out) ) then

! nominally the mixing ratio derivatives but can be used for any
! quantity requiring a simple convolution.

      n_coeffs = SIZE(di_df,dim=2)

      do i = 1, n_coeffs
        if ( .not. di_df_flag(i) ) cycle
        call interpolateValues ( chi_in-init_angle, di_df(:, i), &
        & angles(ffth:no_fft), rad_fft(ffth:no_fft), METHOD='S', &
        & EXTRAPOLATE='C' )

! mirror reflect this

        rad_fft(1:ffth-1) = (/(rad_fft(no_fft-i),i=1, ffth-1)/)

! I don't know if this step is truly necessary but it rephases the radiances
! identically to the prototype code

        rad_fft = cshift(rad_fft,-1)

! take fft of radiance array

        call drft1_t ( rad_fft, 'a' )

! apply convolution theorem

        rad_fft1(1:2) = rad_fft(1:2) * p(1:2)
        do j = 3, no_fft-1, 2
          rad_fft1(j+1) = rad_fft(j) * p(j+1)
          rad_fft1(j)   = rad_fft(j) * p(j)
        end do

        call drft1_t ( rad_fft1, 's' )

! interpolate to output grid

        call interpolateValues ( angles(ffth:no_fft)-ang_step, &
        & rad_fft1(ffth:no_fft), chi_out-init_angle, drad_df_out(:, i), &
        & METHOD='S', EXTRAPOLATE='C' )

      end do

    end if

    init = no_fft

  contains

    subroutine DRFT1_T ( A, Mode )
    ! Call DRFT1 and test its status flag
      real(r8) :: A(*)
      character :: Mode
      call drft1 ( a, mode, pwr, ms, s )
      if ( ms == -2 ) then
        init = 0
        call MLSMessage ( MLSMSG_Error, ModuleName, "Error in drft1" )
      end if
    end subroutine DRFT1_T

  end subroutine Fov_Convolve

! ===========================================     FOV_CONVOLVE_OLD  =====
! This subprogram adds the effects of antenna smearing to the radiance.

  subroutine FOV_CONVOLVE_OLD ( EIL_ANGLE, RADIANCE, DELTA0, ITYPE, NP, &
    &                           M, AntennaPattern, IER )

    use AntennaPatterns_m, only: AntennaPattern_T
    use MLSCommon, only: I4, R8

    Integer(i4), intent(in) :: ITYPE, NP, M

    Real(r8), intent(inout) :: EIL_ANGLE(:)
    Real(r8), intent(inout) :: RADIANCE(:)
    Real(r8), intent(in) :: DELTA0
    type(antennaPattern_T), intent(in) :: AntennaPattern

    Integer(i4), intent(out) :: IER

    Integer(i4) :: I, J, NTR, IAS

    ias = size(antennaPattern%aaap)/2

    ntr = 2**m
    call ftgrid ( eil_angle, radiance, delta0, antennaPattern%lambda, np, ntr )

    j = ntr/2 + 2
    do i = j, ntr
      radiance(i) = radiance(ntr-i+2)
    end do

    i = itype - 1
    if ( i  ==  0 ) then                               ! Straight data
      call convolve ( Radiance, AntennaPattern%AAAP, M, Ias, ier )
    else if ( i  ==  1 ) then                          ! First derivative
      call convolve ( Radiance, AntennaPattern%D1AAP, M, Ias, ier )
    else if ( i  ==  2 ) then                          ! Second derivative
      call convolve ( Radiance, AntennaPattern%D2AAP, M, Ias, ier )
    end if

    return
  contains
! -----------------------------------------------     CONVOLVE     -----
! This subroutine applies the convolution theorem

    Subroutine CONVOLVE ( RADIANCE, AAAP, M, IAS, IERR )

      use DFFT_M, only: DRFT1

      real(r8), intent(inout) :: RADIANCE(:)
      real(r8), intent(in) :: AAAP(*)
      integer(i4), intent(in) :: M, IAS
      integer(i4), intent(out) :: IERR
      Integer, parameter :: MAXP=12, MAX2P=2**MAXP

      Integer(i4), save :: INIT = 0, MS = 0
      Integer(i4) :: ISTOP, M4, NTR, NTRH, I, J
      Real(r8) :: CR, CI, DBLRAD(MAX2P)
      Real(r8), save :: S(MAX2P)

      ierr = 2
      if ( maxp < m) return


      m4 = m
      ierr = 0
      ntr = 2**m
      ntrh = ntr / 2
      dblrad(:ntr) = radiance(:ntr)
      if ( init > 0 .and. init /= m) ms=0
      call drft1 ( dblrad, 'a', m4, ms, s )
      if ( ms == -2 ) then
        init = 0
        ierr=5
      else
        istop = min(ntrh,ias)
        do i = 2, istop
          j = 2*i-1
          cr = dblrad(j)*aaap(j) -   dblrad(j+1)*aaap(j+1)
          ci = dblrad(j)*aaap(j+1) + dblrad(j+1)*aaap(j)
          dblrad(j) = cr
          dblrad(j+1) = ci
        end do
        if ( ntrh > ias ) then
          do i = ias+1,ntrh
            j = 2*i-1
            dblrad(j) = 0.0d0
            dblrad(j+1) = 0.0d0
          end do
        end if
        dblrad(1) = dblrad(1)*aaap(1)
        dblrad(2) = dblrad(2)*aaap(2)
        call drft1 ( dblrad, 's', m4, ms, s )
        if ( ms == -2 ) then
          init = 0
          ierr = 6
        else
          init=m
          radiance(:ntr) = dblrad(:ntr)
        end if
      end if

      return

    End subroutine CONVOLVE

  ! -------------------------------------------------     FTGRID     -----
  ! This subroutine performs the interpolation onto the computational grid
  ! it uses cubic splines

    subroutine FTGRID ( EIL_ANGLE, RADIANCE, DELTA0, XLAMDA, NP, NTR )

      use D_CSPLINE_M, only: CSPLINE

      Real(r8), intent(inout) :: EIL_ANGLE(:)
      Real(r8), intent(inout) :: RADIANCE(:)
      Real(r8), intent(in) :: DELTA0, XLAMDA
      Integer(i4), intent(in) :: NP, NTR

      Integer(i4) :: I, K1, KN, N, MP

      Real(r8) :: X1, XN, R1,RN, V, OOX, PP, DEL
      Real(r8) :: X(np), R(np)

  !  Make sure the EIL_ANGLE is a Monotonically increasing array:

      mp = 1
      r(1) = radiance(1)
      x(1) = eil_angle(1)
      do i = 2, np
        xn = eil_angle(i)
        if ( xn > x(mp) ) then
          mp = mp + 1
          x(mp) = xn
          r(mp) = radiance(i)
        end if
      end do

      oox = 1.0_r8 / xlamda
      pp = 0.185_r8 * oox
      del = oox / ntr

      do i = 1, ntr
        v = delta0 - pp + del * (i - 1)            ! new code
        eil_angle(i) = v
        radiance(i) = 0.0
      end do

      k1 = 1
      r1 = r(1)
      x1 = x(1)
      do while ( eil_angle(k1) <= x1 .and. k1 < ntr )
        radiance(k1) = r1
        k1 = k1 + 1
      end do

      kn = ntr
      rn = r(mp)
      xn = x(mp)
      do while ( eil_angle(kn) >= xn .and. kn > 1 )
        radiance(kn) = rn
        kn = kn - 1
      end do

      n = kn - k1 + 1
      call cspline ( x, eil_angle(k1:kn), r, radiance(k1:kn), mp, n )

      return
    end subroutine FTGRID
  end subroutine FOV_CONVOLVE_OLD

end module Fov_Convolve_m
! $Log$
! Revision 2.7  2002/09/05 20:48:29  livesey
! Fixed (i.e. added) handling of vmr derivative flags.
!
! Revision 2.6  2002/08/02 00:12:19  bill
! just testing
!
! Revision 2.5  2002/07/08 17:45:40  zvi
! Replace local pointers by automatic arrays
!
! Revision 2.4  2002/06/28 11:06:48  zvi
! Add computing of dI/d(Chi) on output grid
!
! Revision 2.3  2002/06/19 11:00:32  zvi
! Some cosmetic corrections
!
! Revision 2.2 2002/06/17 16:31:21 bill
! Add zvis modification, rename module
!
! Revision 2.0  2001/09/17 20:26:26  livesey
! New forward model
!
! Revision 1.12  2001/05/02 20:49:23  zvi
! Cleaning up code
!
! Revision 1.11  2001/04/10 01:16:34  livesey
! Tidied up convolution
!
! Revision 1.10  2001/04/09 23:32:29  zvi
! Correcting a small error in radiances folding code
!
! Revision 1.9  2001/04/06 01:37:58  zvi
! Put (*) (Assume size) status on CONVOLVE & DFFT arrays..
!
! Revision 1.8  2001/04/05 22:54:39  vsnyder
! Use AntennaPatterns_M
!
! Revision 1.7  2001/03/31 23:40:55  zvi
! Eliminate l2pcdim (dimension parameters) move to allocatable ..
!
! Revision 1.6  2001/03/29 08:51:01  zvi
! Changing the (*) toi (:) everywhere
!
! Revision 1.5  2001/02/26 09:01:16  zvi
! New version - Using "Super-Structures"
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
