MODULE fov_convolve_m
 
  USE Allocate_Deallocate, only: allocate_test, deallocate_test
  USE AntennaPatterns_m, only: AntennaPattern_T
  USE D_CSPLINE_M, only: CSPLINE
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
! ============================================  fov_convolve =====
! This subprogram adds the effects of antenna smearing to the radiance.
!
  SUBROUTINE fov_convolve(AntennaPattern,chi_in,rad_in,chi_out,rad_out, &
           & req,rsc,earth_frac,surf_angle,di_dt,dx_dt,ddx_dxdt,dx_dt_out, &
           & drad_dt_out,di_df,drad_df_out,drad_dx_out)
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
!
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
  REAL(rp), OPTIONAL, INTENT(out) :: drad_dx_out(:) ! outputted derivative 
!                                    of radiance wrt to Chi_out
  REAL(rp), OPTIONAL, INTENT(out) :: drad_dt_out(:,:) ! outputted radiance
!                                    derivatives wrt temperature.
  REAL(rp), OPTIONAL, INTENT(out) :: drad_df_out(:,:) ! outputted radiance
!                                    derivatives for inputted di_df.
! Internal stuff
!
  INTEGER(i4) :: i, j, k, ffth, n_coeffs, zero_out
  INTEGER(i4) :: peak_grad(1)
!
  REAL(r8) :: r_eq, r_sc, e_frac, init_angle, aaap_step, ang_step
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
  ang_step = 1.0_r8 / (no_fft * aaap_step)
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
  angles = (/(i*ang_step,i=1,no_fft)/)
  angles = angles - angles(ffth+1)
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
  end do
!
  CALL drft1(rad_fft1,'s',pwr,ms,s)
  if(ms == -2) then
    init=0
    call MLSMessage ( MLSMSG_Error, ModuleName,"Error in drft1")
  endif
!
! interpolate to output grid
!
  IF ( PRESENT(drad_dx_out) ) THEN
    CALL interpolatevalues(angles(ffth:no_fft)-ang_step, &
       & rad_fft1(ffth:no_fft),chi_out-init_angle,rad_out,&
       & METHOD='S',EXTRAPOLATE='C',dyByDx=drad_dx_out)
  ELSE
    CALL interpolatevalues(angles(ffth:no_fft)-ang_step, &
       & rad_fft1(ffth:no_fft),chi_out-init_angle,rad_out,&
       & METHOD='S',EXTRAPOLATE='C')
  ENDIF
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
    end do
!
    CALL drft1(rad_fft1,'s',pwr,ms,s)
    if(ms == -2) then
      init=0
      call MLSMessage ( MLSMSG_Error, ModuleName,"Error in drft1")
    endif
!
! interpolate to output grid
!
    CALL interpolatevalues(angles(ffth:no_fft)-ang_step, &
       & rad_fft1(ffth:no_fft), chi_out-init_angle,drad_dt_temp, &
       & METHOD='S',EXTRAPOLATE='C')
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
      rad_fft1(1:ffth-1) = (/(rad_fft1(no_fft-i),i=1,ffth-1)/)
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
      end do
!
! interplolate to chi_out
!
      CALL drft1(rad_fft,'s',pwr,ms,s)
      if(ms == -2) then
        init=0
        call MLSMessage ( MLSMSG_Error, ModuleName,"Error in drft1")
      endif
 
      CALL interpolatevalues(angles(ffth:no_fft)-ang_step, &
         & rad_fft(ffth:no_fft), chi_out-init_angle, drad_dt_out(:,i), &
         & METHOD='S',EXTRAPOLATE='C')
!
! compute final result
!
      drad_dt_out(:,i) = drad_dt_out(:,i) + dx_dt_out(:,i)*drad_dt_temp
!
    end do               ! On i = 1, n_coeffs
!
    CALL deallocate_test(drad_dt_temp,'drad_dt_temp',modulename)
!
  ENDIF
!
  IF (PRESENT(drad_df_out)) THEN
!
! nominally the mixing ratio derivatives but can be used for any
! quantity requiring a simple convolution.
!
    n_coeffs = SIZE(di_df,dim=2)
!
    DO i = 1, n_coeffs
      CALL interpolatevalues(chi_in-init_angle,di_df(:,i), &
      & angles(ffth:no_fft), rad_fft(ffth:no_fft),METHOD='S', &
      & EXTRAPOLATE='C')
!
! mirror reflect this
!
      rad_fft(1:ffth-1) = (/(rad_fft(no_fft-i),i=1, ffth-1)/)
!
! I don't know if this step is truly necessary but it rephases the radiances
! identically to the prototype code
!
      rad_fft = CSHIFT(rad_fft,-1)
!
! take fft of radiance array
!
      CALL drft1(rad_fft,'a',pwr,ms,s)
!
! apply convolution theorem
!
      rad_fft1(1:2) = rad_fft(1:2) * p(1:2)
      DO j = 3, no_fft-1, 2
        rad_fft1(j+1) = rad_fft(j) * p(j+1)
        rad_fft1(j)   = rad_fft(j) * p(j)
      end do
!
      CALL drft1(rad_fft1,'s',pwr,ms,s)
!
! interpolate to output grid
!
      CALL interpolatevalues(angles(ffth:no_fft)-ang_step, &
      & rad_fft1(ffth:no_fft), chi_out-init_angle,drad_df_out(:,i), &
      & METHOD='S',EXTRAPOLATE='C')
!
    end do
!
  ENDIF
!
  init = no_fft
 
  CALL deallocate_test(p,'p',modulename)
  CALL deallocate_test(dp,'dp',modulename)
  CALL deallocate_test(angles,'angles',modulename)
  CALL deallocate_test(rad_fft,'rad_fft',modulename)
  CALL deallocate_test(rad_fft1,'rad_fft1',modulename)
!
 END SUBROUTINE fov_convolve

! ===========================================     FOV_CONVOLVE_OLD  =====
! This subprogram adds the effects of antenna smearing to the radiance.
!
  SUBROUTINE FOV_CONVOLVE_OLD(EIL_ANGLE, RADIANCE, DELTA0, ITYPE, NP, &
 &                        M, AntennaPattern, IER )
!
    Integer(i4), intent(in) :: ITYPE, NP, M

    Real(r8), intent(inout) :: EIL_ANGLE(:)
    Real(r8), intent(inout) :: RADIANCE(:)
    Real(r8), intent(in) :: DELTA0
    type(antennaPattern_T), intent(in) :: AntennaPattern

    Integer(i4), intent(out) :: IER

    Integer(i4) :: I, J, NTR, IAS

    ias = size(antennaPattern%aaap)/2
!
    ntr = 2**m
    call ftgrid(eil_angle,radiance,delta0,antennaPattern%lambda,np,ntr)
!
    j = ntr/2 + 2
    do i = j, ntr
      radiance(i) = radiance(ntr-i+2)
    end do
!
    i = itype - 1
    if (i  ==  0) then                               ! Straight data
      call convolve(Radiance,AntennaPattern%AAAP,M,Ias,ier)
    else if (i  ==  1) then                          ! First derivative
      call convolve(Radiance,AntennaPattern%D1AAP,M,Ias,ier)
    else if (i  ==  2) then                          ! Second derivative
      call convolve(Radiance,AntennaPattern%D2AAP,M,Ias,ier)
    end if
!
    Return
  End subroutine FOV_CONVOLVE_OLD
! *****     Private procedures     *************************************
!
! -----------------------------------------------     CONVOLVE     -----
! This subroutine applies the convolution theorem
!
  Subroutine CONVOLVE ( RADIANCE, AAAP, M, IAS, IERR )
    real(r8), intent(inout) :: RADIANCE(:)
    real(r8), intent(in) :: AAAP(*)
    integer(i4), intent(in) :: M, IAS
    integer(i4), intent(out) :: IERR
    Integer, parameter :: MAXP=12, MAX2P=2**MAXP
!
    Integer(i4), save :: INIT = 0, MS = 0
    Integer(i4) :: ISTOP, M4, NTR, NTRH, I, J
    Real(r8) :: CR, CI, DBLRAD(MAX2P)
    Real(r8), save :: S(MAX2P)
!
    ierr = 2
    if(maxp < m) return
!

    m4 = m
    ierr = 0
    ntr = 2**m
    ntrh = ntr / 2
    dblrad(:ntr) = radiance(:ntr)
    if (init > 0 .and. init /= m) ms=0
    call drft1 ( dblrad, 'a', m4, ms, s )
    if(ms == -2) then
      init=0
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
      if(ntrh > ias) then
        do i = ias+1,ntrh
          j = 2*i-1
          dblrad(j) = 0.0d0
          dblrad(j+1) = 0.0d0
        end do
      endif
      dblrad(1) = dblrad(1)*aaap(1)
      dblrad(2) = dblrad(2)*aaap(2)
      call drft1 ( dblrad, 's', m4, ms, s )
      if(ms == -2) then
        init=0
        ierr = 6
      else
        init=m
        radiance(:ntr) = dblrad(:ntr)
      endif
    endif
!
    Return

  End subroutine CONVOLVE
!
! -------------------------------------------------     FTGRID     -----
! This subroutine performs the interpolation onto the computational grid
! it uses cubic splines
!
  Subroutine FTGRID ( EIL_ANGLE, RADIANCE, DELTA0, XLAMDA, NP, NTR )

    Real(r8), intent(inout) :: EIL_ANGLE(:)
    Real(r8), intent(inout) :: RADIANCE(:)
    Real(r8), intent(in) :: DELTA0, XLAMDA
    Integer(i4), intent(in) :: NP, NTR

    Integer(i4) :: I, K1, KN, N, MP
!
    Real(r8) :: X1, XN, R1,RN, V, OOX, PP, DEL
    Real(r8) :: X(np), R(np)
!
!  Make sure the EIL_ANGLE is a Monotonically increasing array:
!
    mp = 1
    r(1) = radiance(1)
    x(1) = eil_angle(1)
    do i = 2, np
      xn = eil_angle(i)
      if(xn > x(mp)) then
        mp = mp + 1
        x(mp) = xn
        r(mp) = radiance(i)
      endif
    end do
!
    oox = 1.0_r8 / xlamda
    pp = 0.185_r8 * oox
    del = oox / ntr
!
    do i = 1, ntr
      v = delta0 - pp + del * (i - 1)            ! new code
      eil_angle(i) = v
      radiance(i) = 0.0
    end do
!
    k1 = 1
    r1 = r(1)
    x1 = x(1)
    do while ( eil_angle(k1) <= x1 .and. k1 < ntr )
      radiance(k1) = r1
      k1 = k1 + 1
    end do
!
    kn = ntr
    rn = r(mp)
    xn = x(mp)
    do while ( eil_angle(kn) >= xn.and.kn > 1 )
      radiance(kn) = rn
      kn = kn - 1
    end do
!
    n = kn - k1 + 1
    call cspline(x, eil_angle(k1:kn), r, radiance(k1:kn), mp, n)
!
    return
  End subroutine FTGRID

!
!=====================================================================
!=================== Coversion of the f77 fft.f to f90 ===============
!  The JPL Double Precision FFT Package.
!---------------------------------------------------------------------
!
      SUBROUTINE DRFT1 (A,MODE,M,MS,S)
!
!>> 1994-11-11 DRFT1  Krogh   Declared all vars.
!>> 1994-10-20 DRFT1 Krogh  Changes to use M77CON
!>> 1989-05-07 DRFT1 FTK & CLL
!>> 1989-04-21 FTK & CLL
!     .  Copyright (C) 1989, California Institute of Technology.
!     .  U. S. Government sponsorship under
!     .  NASA contract NAS7-918 is acknowledged.
!
!     This subroutine computes Fourier transforms of real data using
!     the Cooley-Tukey fast Fourier transform.
!
!     Variables in the calling sequence have the following types
      DOUBLE PRECISION A(*), S(*)
      INTEGER          M, MS
      CHARACTER        MODE
!
!     Programmed by Fred T. Krogh at the Jet Propulsion Laboratory,
!     Pasadena, Calif.   August 1, 1969.
!     Revised for portability by Krogh -- January 22, 1988
!
!     In describing the usage the following notation is used
!     N = 2 ** M
!     W = EXP(2*PI*I/N), where I = SQRT(-1) and PI = 3.14159...
!
!     The usage is as follows
!
! A() is an array of function values if one is doing Fourier analysis,
!  and is an array of Fourier coefficients if one is doing Fourier
!  synthesis.  In our description here we assume that A is a real
!  array with dimension A(N) when A contains the real data, X, and
!  that A is a complex array with dimension A(N/2) when A contains
!  complex Fourier coefficients, C.  (C(k) for k > N/2 need not be
!  saved since for 0 < k < N, C(N-k) = conjugate of C(k).  It is
!  assumed that the imaginary part of a complex number is stored
!  in the cell immediately following its real part, except that
!  A(1) = C(0), and A(2) = C(N/2).  This is possible since these
!  values of C are real and by doing this both X and C require the
!  same storage in A. Of course the storage required for A can be
!  reserved by the user in any way that works.
!
! MODE  Selects Synthesis or Analysis.
!  If MODE = 'A' or 'a', do Fourier analysis, which amounts to setting
!  C(k) = sum for j=0 to N-1 of X(j)*T(M,j,k), for k = 0, N/2
!  with  T(M,j,k) = (1/N) * W ** (-j*k).
!  If MODE = 'S' or 's', do Fourier synthesis, which amounts to setting
!  X(j) = sum for k=0 to N-1 of C(k)*T(M,j,k), for j = 0, N - 1
!  with  T(M,j,k) = W ** (j*k)
!  (Recall that C(N-k) is the conjugate of C(k).)
!
! M is used to indicate N = 2**M, the number of real points.  The
!  number of points must satisfy 1 .le. N .le. 2**21.
!  M = 0 gives an immediate return.
!
! MS gives the state of the sine table.  If MS > 0, there are NT =
!    2 ** (MS-2) good entries in the sine table.  On the initial call,
!    MS must be set to 0, and when the return is made, it will be set
!    to M, which is the value of MS required for computing a
!    transform of size N.  If MS = -1, the sine table will be computed
!    as for the case MS = 0, and THEN a return to the user will be made
!    with MS set as before, but no transform will be computed.  This
!    option is useful if the user would like access to the sine table
!    before computing the FFT.
!    On detected errors the error message subrs are called and
!    execution stops.  If the user overrides the stop to cause
!    continuation, THEN this subr will return with MS = -2.
!
! S() is a vector, S(j) = sin(pi*j/2*NT)), j = 1, 2, ..., NT-1, where
!  NT is defined in the description of MS above.  S is computed by the
!  subroutine if M .gt. MS.  (If S is altered, set MS=0 so that S
!  is recomputed.)
!
!     ------------------------------------------------------------------
!                Notes on COMMON, PARAMETER's, and local variables
!
!     MMAX is the largest value allowed for M
!     The dimension of KE must be at least as large as MMAX-1.
!     The named common CDFFTC is used for communication between this
!     subroutine and the subroutine DFFT which computes a one
!     dimensional complex Fourier transform and computes the sine table.
!     The use of the variables in CDFFTC is contained in the listing
!     of DFFT.
!
!     ANAL = .TRUE. when doing Fourier analysis, and .false. otherwise.
!
!     N1 = 2 ** M
!     ------------------------------------------------------------------
!--D replaces "?": ?RFT1, ?FFT, C?FFTC
!     Both versions use ERMSG, IERM1
!     and need ERFIN, IERV1
!     ------------------------------------------------------------------
      INTEGER MMAX
      INTEGER I, II, II1, II2, IR, IR1, IR2
      INTEGER J, JDIF, JJ
      INTEGER K1, K1N, KN2
      INTEGER L
      INTEGER MA, MSI
      INTEGER N1, N1P, KEDIM
 
      DOUBLE PRECISION FN, HALF
      DOUBLE PRECISION T, TI, TT, TTI, TWO, WI, WR
 
      LOGICAL ANAL
 
      PARAMETER (TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)
      EQUIVALENCE (ILAST, N1)
! Common variables
      PARAMETER (KEDIM=20)
      LOGICAL NEEDST
      INTEGER MT, NT, MM, KS, ILAST, KE(KEDIM), KEE(KEDIM+1)
! Note that KEE(1) is equivalent to ILAST.
      EQUIVALENCE (KE(1), KEE(2))
      COMMON /CDFFTC/ NEEDST, MT, NT, MM, KS, ILAST, KE
      SAVE /CDFFTC/
      PARAMETER (MMAX = KEDIM+1)
!     ------------------------------------------------------------------
!
      if( MODE .eq. 'A' .or. MODE .eq. 'a') THEN
         ANAL = .true.
      elseif( MODE .eq. 'S' .or. MODE .eq. 's') THEN
         ANAL = .false.
      else
         Print *,'** Fatal error in DRFT1 **'
         Print *,'   Bad MODE.  MODE = ',MODE
         MS = -2
         return
      endif
      MA = M
      MSI = MS
      NEEDST = MA .GT. MSI
      IF(NEEDST) THEN
         IF(MA .GT. MMAX .or. MA .lt. 0) THEN
!                               Fatal error, default is to stop in IERM1
            Print *,'** Fatal error in DRFT1 **'
            Print *,'   Require (0 .le. M .le. 21), M=',M
            MS = -2
            RETURN
         ENDIF
         MS = MA
         MT = MA - 2
         CALL DFFT (A, A, S)
         IF(MSI .EQ. -1) RETURN
      ENDIF
      IF(MA .NE. 0) THEN
         MM = MA - 1
         N1 = 2 ** MA
         N1P = N1 + 2
         KN2 = N1 / 2
         JDIF = (4 * NT) / N1
         KS = 2
         IF(ANAL) THEN
!                               Set flags for Fourier analysis
            IR = 2
            II = 1
            FN = HALF / DBLE(N1)
!           Doing Fourier analysis, so multiply by 2 ** M
            DO 10 I = 1, N1
               A(I) = A(I) * FN
   10       CONTINUE
         ELSE
!                              Set flags for Fourier synthesis
            IR = 1
            II = 2
            GOTO 40
         ENDIF
 
!                              Compute complex Fourier transform
   20    DO 30 L = 1, MM
            KEE(L+1) = KEE(L) / 2
   30    CONTINUE
!
         CALL DFFT (A(IR), A(II), S)
!                              End of computing complex transform
!
         IF(.NOT. ANAL) RETURN
!
!        Beginning of calculations relating coefficients of real data
!        with coefficients of associated complex data
!
!        Special case --  K1 = 0
   40    T = A(1) + A(2)
         TI = A(1) - A(2)
         IF(ANAL) THEN
            T = TWO * T
            TI = TWO * TI
         ENDIF
         A(1) = T
         A(2) = TI
         IF(MM .GT. 0) THEN
!                           Special kase -- K1 = N1 / 4
            A(KN2+1) = TWO * A(KN2+1)
            A(KN2+2) = -TWO * A(KN2+2)
            IF(MM .GT. 1) THEN
               J = 0
               DO 50 K1 = 3, KN2, 2
                  K1N = N1P - K1
                  IF(ANAL) THEN
                     IR1 = K1N
                     IR2 = K1
                  ELSE
                     IR1 = K1
                     IR2 = K1N
                  ENDIF
                  II2 = IR2 + 1
                  II1 = IR1 + 1
                  J = J + JDIF
                  WI = S(J)
                  JJ = NT - J
                  WR = S(JJ)
                  T = A(IR1) - A(IR2)
                  TI = A(II1) + A(II2)
                  TT = T * WI + TI * WR
                  TTI = T * WR - TI * WI
                  T = A(IR1) + A(IR2)
                  TI = A(II1) - A(II2)
                  A(IR1) = T - TT
                  A(IR2) = T + TT
                  A(II1) = TTI + TI
                  A(II2) = TTI - TI
   50          CONTINUE
            ENDIF
         ENDIF
         IF(.NOT. ANAL) GOTO 20
      ENDIF
!
      RETURN
!
      END SUBROUTINE DRFT1
!
!==========================================================================
!
      SUBROUTINE DFFT(AR, AI, S)
!
!>> 1994-11-11 DFFT  Krogh   Declared all vars.
!>> 1994-10-20 DFFT  Krogh  Changes to use M77CON
!>> 1989-08-15 DFFT  FTK -- Parameter constants given to more precision.
!>> 1989-04-21 DFFT  FTK & CLL
!     .  Copyright (C) 1989, California Institute of Technology.
!     .  All rights reserved.  U. S. Government sponsorship under
!     .  NASA contract NAS7-918 is acknowledged..
!
!     This subroutine computes one dimensional complex Fourier
!     transforms using the Cooley-Tukey algorithm. It is used by a
!     number of subroutines which compute Fourier transforms.
!
!     Programmed by Fred T. Krogh at the Jet Propulsion Laboratory,
!     Pasadena, Calif.   August 1, 1969.
!     Revised by Krogh at JPL -- January 19, 1988 -- For portability
!
      DOUBLE PRECISION AR(*), AI(*), S(*)
!     Minimum dimensions are AR(ILAST), AI(ILAST), S(NT-1), where ILAST
!     and NT are defined in the common block below.
!
! AR and AI give the arrays used to hold the real and imaginary parts of
!     the Fourier coefficients and data.
! S   contains the sine table required in the calculations.
!
!     -----------------------------------------------------------------
!                Notes on COMMON, PARAMETER's, and local variables
!
!     SPI4 = SIN(PI/4) = 0.5 * SQRT(2)
!     PI4 = PI / 4
!     THE DIMENSION OF KE MUST BE AS LARGE AS THE MAXIMUM VALUE
!     PERMITTED FOR MM.
!
!     WHERE
!     NEEDST= .TRUE. if the sine table must be computed.
!     MT    = base 2 log(NT)
!     NT    = number of entries in the sine table
!     MM    = base 2 log(number of complex fourier coefficients)
!     KS    = distance in memory between successive points.  The i-th
!             coefficient, a(i) is given by AR((I+1)*KS)+AI((I+1)*KS)*
!             sqrt(-1), i=0, 1, ..., (2**MM)-1.
!     ILAST = KS * 2**MM
!     KE(L) = KS * 2**(MM-L), L=1, 2, ..., MM
!     -----------------------------------------------------------------
!--D replaces "?": ?FFT, C?FFTC
!     -----------------------------------------------------------------
      INTEGER I, I1, I2, I3, IJ
      INTEGER J, JDIF, JGO, JI, JI2, JJ, JR
      INTEGER K, KSI
      INTEGER L, L1, L4, LJ, LL
      INTEGER MTC
 
      DOUBLE PRECISION HALF, PI4, SPI4
      DOUBLE PRECISION T, T1, T1I, T2, T2I, T3, T3I, THETA
      DOUBLE PRECISION TI, TP, TP1, TP1I, TPI
      DOUBLE PRECISION WI1, WI2, WI3, WR1, WR2, WR3
 
      PARAMETER (SPI4 = 0.7071067811865475244008443621048490D0)
      PARAMETER (PI4 = 0.7853981633974483096156608458198757D0)
      PARAMETER (HALF = 0.5D0)
 
! Common variables
      INTEGER KEDIM
      PARAMETER (KEDIM=20)
      LOGICAL NEEDST
      INTEGER MT, NT, MM, KS, ILAST, KE(KEDIM), KEE(KEDIM+1)
! Note that KEE(1) is equivalent to ILAST.
      EQUIVALENCE (KE(1), KEE(2))
      COMMON /CDFFTC/ NEEDST, MT, NT, MM, KS, ILAST, KE
      SAVE /CDFFTC/
!     -----------------------------------------------------------------
!
      IF(NEEDST) THEN
!                      Compute the sine table
         NEEDST = .FALSE.
         MTC = MT
         NT = 2**MTC
         IF(MTC .GT. 0) THEN
!                            SET FOR L=1
            J = NT
            JJ = J / 2
            S(JJ) = SPI4
            IF(MTC .GT. 1) THEN
               THETA = PI4
               DO 20 L = 2, MTC
                  THETA = HALF * THETA
                  K = J
                  J = JJ
                  JJ = JJ / 2
!                       At this point J = 2**(MT-L+1) and JJ = 2**(MT-L)
                  S(JJ) = SIN(THETA)
                  L1 = NT - JJ
                  S(L1) = COS(THETA)
                  LL = NT - K
                  IF(LL .GE. J) THEN
                     DO 10 I = J, LL, J
                        I1 = NT - I
                        I2 = I + JJ
                        S(I2) = S(I) * S(L1) + S(JJ) * S(I1)
   10                CONTINUE
                  ENDIF
   20          CONTINUE
            ENDIF
         ENDIF
      ELSE
!                      Compute the transform
!                           Scramble A
!
         IJ = 1
         JI = 1
         L1 = KS
         IF(MM .GT. 1) THEN
   30       IJ = IJ + L1
            DO 40 I = 1, MM
               JI = JI + KE(I)
               KE(I) = -KE(I)
               IF(KE(I) .LT. 0) THEN
                  IF(IJ .LT. JI) THEN
!                       INTERCHANGE THE IJ-TH COEFFICIENT WITH THE JI-TH
                     T = AR(IJ)
                     AR(IJ) = AR(JI)
                     AR(JI) = T
                     T = AI(IJ)
                     AI(IJ) = AI(JI)
                     AI(JI) = T
                  ENDIF
                  GOTO 30
               ENDIF
   40       CONTINUE
         ENDIF
!                          END OF SCRAMBLING A
         JDIF = NT
         IF(MOD(MM, 2) .NE. 0) THEN
!                    SPECIAL CASE -- L = 1,  MM ODD  (RADIX 2 ALGORITHM)
            L1 = L1 + L1
            DO 50 I = 1, ILAST, L1
               KSI = I + KS
               T = AR(I)
               AR(I) = T + AR(KSI)
               AR(KSI) = T - AR(KSI)
               T = AI(I)
               AI(I) = T + AI(KSI)
               AI(KSI) = T - AI(KSI)
   50       CONTINUE
            JDIF = JDIF / 2
         ENDIF
!
         DO 140 L = 2, MM, 2
            L4 = 4 * L1
            LJ = L1 / 2
            J = 0
            JI = 0
!
!           ASSIGN 70 TO JGO
            JGO = 70
!
!           START OF I LOOP  (RADIX 4 ALGORITHM)
   60       IJ = J + 1
            DO 120 I = IJ, ILAST, L4
               I1 = I + L1
               I2 = I1 + L1
               I3 = I2 + L1
!
!              GOTO JGO, (70, 80, 90)
!
               if(JGO.eq.70) then
                 goto 70
               else if(JGO.eq.80) then
                 goto 80
               else if(JGO.eq.90) then
                 goto 90
               endif
!
!              SPECIAL CASE -- J = 0
   70          T = AR(I) + AR(I1)
               T1 = AR(I) - AR(I1)
               TI = AI(I) + AI(I1)
               T1I = AI(I) - AI(I1)
               T2 = AR(I2) + AR(I3)
               T3 = AR(I2) - AR(I3)
               T2I = AI(I2) + AI(I3)
               T3I = AI(I2) - AI(I3)
               GOTO 110
!
!              SPECIAL CASE -- J = L1 / 2
   80          T2 = SPI4 * AR(I2)
               T2I = SPI4 * AI(I2)
               T3 = SPI4 * AR(I3)
               T3I = SPI4 * AI(I3)
               TP = T2 - T2I
               TPI = T2 + T2I
               TP1 = -T3 - T3I
               TP1I = T3 - T3I
               T = AR(I) - AI(I1)
               T1 = AR(I) + AI(I1)
               TI = AI(I) + AR(I1)
               T1I = AI(I) - AR(I1)
               GOTO 100
!
!              USUAL CASE -- J .NE. 0  AND  J .NE. L1 / 2
!
!              WRK AND WIK (K = 1, 2, 3) ARE THE REAL AND IMAGINARY PART
!              RESP. OF EXP(SQRT(-1) * PI * K*(2**(-L-MOD(MM, 2)))*J/KS)
!                             =EXP(SQRT(-1) * PI * K * (J / (4 * L1)))
!
   90          T2 = WR2 * AR(I1) - WI2 * AI(I1)
               T2I = WI2 * AR(I1) + WR2 * AI(I1)
               T = AR(I) + T2
               T1 = AR(I) - T2
               TI = AI(I) + T2I
               T1I = AI(I) - T2I
               TP = WR1 * AR(I2) - WI1 * AI(I2)
               TPI = WI1 * AR(I2) + WR1 * AI(I2)
               TP1 = WR3 * AR(I3) - WI3 * AI(I3)
               TP1I = WI3 * AR(I3) + WR3 * AI(I3)
  100          T2 = TP + TP1
               T3 = TP - TP1
               T2I = TPI + TP1I
               T3I = TPI - TP1I
  110          AR(I) = T + T2
               AI(I) = TI + T2I
               AR(I1) = T1 - T3I
               AI(I1) = T1I + T3
               AR(I2) = T - T2
               AI(I2) = TI - T2I
               AR(I3) = T1 + T3I
               AI(I3) = T1I - T3
  120       CONTINUE
!           END OF I LOOP
!
            IF(J .LT. LJ) THEN
               IF(J .EQ. 0) THEN
!                 ASSIGN 90 TO JGO
                  JGO = 90
                  J = KS
               ELSE
                  J = L1 - J
!                 COMPUTE WR-S AND WI-S FOR J REPLACED BY L1 - J
                  T = WI1
                  WI1 = WR1
                  WR1 = T
                  WR2 = -WR2
                  T = -WI3
                  WI3 = -WR3
                  WR3 = T
                  GOTO 60
               ENDIF
            ELSE IF(J .EQ. LJ) THEN
               GOTO 130
            ELSE
               J = L1 - J + KS
            ENDIF
            IF(J .LT. LJ) THEN
!                              COMPUTE WR-S AND WI-S
               JI = JI + JDIF
               JR = NT - JI
               WR1 = S(JR)
               WI1 = S(JI)
               JI2 = JI + JI
               WI2 = S(JI2)
               JR = NT - JI2
               WR2 = S(JR)
               JI2 = JI + JI2
               IF(JI2 .LE. NT) THEN
                  JR = NT - JI2
                  WR3 = S(JR)
                  WI3 = S(JI2)
                  GOTO 60
               ENDIF
               JR = JI2 - NT
               JI2 = NT - JR
               WI3 = S(JI2)
               WR3 = -S(JR)
               GOTO 60
            ELSE IF(J .EQ. LJ) THEN
!                                    SET FOR J = L1 / 2
!              ASSIGN 80 TO JGO
               JGO = 80
               GOTO 60
            ENDIF
!           END OF COMPUTING WR-S AND WI-S
!
  130       L1 = L4
            JDIF = JDIF / 4
  140    CONTINUE
!        END OF L LOOP
      ENDIF
!
      RETURN
 
      END SUBROUTINE DFFT
!
!=====================================================================
!
END MODULE fov_convolve_m
! $Log$
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
