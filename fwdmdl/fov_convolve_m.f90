module FOV_CONVOLVE_M
  use D_FFT, only: RFT1
  use L2PCdim, only: MAXFFT, NLVL
  use MACHINE, only: IO_ERROR
  use MLSCommon, only: I4, R4, R8
  use S_CSPLINE_M, only: CSPLINE
  implicit NONE
  private
  public :: FOV_CONVOLVE

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------

contains

! ===========================================     FOV_CONVOLVE     =====
! This subprogram adds the effects antenna smearing to the radiance
!
  SUBROUTINE FOV_CONVOLVE( EIL_ANGLE, RADIANCE, DELTA0, ITYPE, NP, &
 &                         MBAND, M, InDir, AntN, IER )
!
    Integer, parameter :: mx3 = 3*maxfft
!
    Real(r4), intent(inout) :: EIL_ANGLE(*)
    Real(r4), intent(inout) :: RADIANCE(*)
    Real(r4), intent(in) :: DELTA0
    Integer(i4), intent(in) :: ITYPE, NP, MBAND, M
    Character(len=*), intent(in) :: InDir, AntN
    Integer(i4), intent(out) :: IER

    Real(r4), save :: AAAP(maxfft,3) = 0.0_r4, D1AAP(maxfft,3)=0.0_r4, &
                      D2AAP(maxfft,3) = 0.0_r4
    Character(len=132) Fn
    integer(i4), save :: IAS, INIT = 0
    Integer(i4) :: I, J, IND, LA, LF, NTR
    real(r4) :: X
    real(r4), save :: XLAMDA

    if ( init == 0 ) then
      fn(1:) = ' '
      la = len_trim(antn)
      lf = len_trim(indir)
      fn = indir(1:lf)//antn(1:la)
      lf = len_trim(fn)
      call antenna ( fn(1:lf), m, xlamda, aaap, d1aap, d2aap, ias, ier )
      if (ier /= 0) return
      init = m
    endif
!
    ntr = 2**m
    call ftgrid(eil_angle,radiance,delta0,xlamda,np,ntr,ier)
    if(ier /= 0) return
!
    j = ntr / 2 + 2
    do i = j,ntr
      radiance(i) = radiance(ntr-i+2)
    end do
!
    x = mband - 1
    ind = 1 + sqrt(x)
!
    i = itype - 1
    if (i  ==  0) then                               ! Straight data
      call convolve(radiance, aaap,m,ias,ind,ier)
    else if (i  ==  1) then                          ! First derivative
      call convolve(radiance,d1aap,m,ias,ind,ier)
    else if (i  ==  2) then                          ! Second derivative
      call convolve(radiance,d2aap,m,ias,ind,ier)
    endif
!
    Return
  End subroutine FOV_CONVOLVE

! *****     Private procedures     *************************************
! ------------------------------------------------     ANTENNA     -----
! This subroutine reads an external antenna aperture autocorrelation
! file
!
  Subroutine ANTENNA ( Fn, M, XLAMDA, AAAP, D1AAP, D2AAP, IAS, IERR )
    use UNITS, only: AAAP_UNIT

    Integer(i4), parameter :: MaxV= 2048
!
    Real(r4), intent(out) :: AAAP(*),D1AAP(*),D2AAP(*),XLAMDA
    integer(i4), intent(out) :: IAS, IERR
!
    Integer(i4) :: M, k, lf, I, J, NTR, L, N
    REAL(r8) :: V(6), VALL(MaxV, 6), dx2p, Q, Q2
!
    Character*(*) Fn
!
    ierr = 0
    lf = len_trim(Fn)
    OPEN(aaap_unit,FILE=Fn(1:lf),action='READ',STATUS='OLD',iostat=k)
    if(k /= 0) then
      ierr = 1
      Print *,'** File: ',Fn(1:lf)
!     Call ErrMsg('** From ANTENNA subroutine',k)
      call io_error ('Opening file in FOV_CONVOLVE%ANTENNA', k, Fn(1:lf) )
      goto 20
    endif
!
    READ(aaap_unit,*,iostat=k) XLAMDA
    if(k /= 0) then
      ierr = 1
      Print *,'** Reading error in file: ',Fn(1:lf)
!     Call ErrMsg('** From ANTENNA subroutine',k)
      call io_error ( 'Reading file in FOV_CONVOLVE%ANTENNA', k, Fn(1:lf) )
      goto 20
    endif
!
    dx2p = 6.28318530717959_r8 * Xlamda         ! 2 * Pi * Lambda
!
    i = 0
    k = 0
!
    do while(k == 0)
!
! This loop MUST exit on an end of file condition
!
      read(aaap_unit,*,iostat=k) v
      if(k == -1) exit
      if(k == 0) then
        if(I == MaxV) then
          ierr = 1
          Print *,'** Error in ANTENNA subroutine !!'
          Print *,'   Too many lines in: ',Fn(1:lf)
          Print *,'   Maximum allowed is:',MaxV
          goto 20
        endif
        i = i + 1
        vall(i,:) = v
      else
        ierr = 1
        Print *,'** Reading error in file: ',Fn(1:lf)
!       Call ErrMsg('** From ANTENNA subroutine',k)
        call io_error ( 'Reading file in FOV_CONVOLVE%ANTENNA', k, Fn(1:lf) )
        goto 20
      endif
!
    end do
!
    ias = i
    ntr = 2**m
!
    do i = 1, ias
!
      j = 2 * i - 1
      l = 2 * ias + j
      n = 4 * ias + j
!
      v = vall(i,:)
!
      aaap(j)   = v(1)
      aaap(j+1) = v(2)
      aaap(l)   = v(3)
      aaap(l+1) = v(4)
      aaap(n)   = v(5)
      aaap(n+1) = v(6)
!
! Below is the  first derivative field:     ! i*Q * F(S), i = Sqrt(-1)
!
      q = (i - 1) * dx2p
      d1aap(j)    = -v(2) * q
      d1aap(j+1)  =  v(1) * q
      d1aap(l)    = -v(4) * q
      d1aap(l+1)  =  v(3) * q
      d1aap(n)    = -v(6) * q
      d1aap(n+1)  =  v(5) * q
!
! Below is the second derivative field:     ! (i*Q)**2 * F(S), i = Sqrt(-1)
!
      q2 = -q * q                         ! (i*q)**2 = -q*q
      d2aap(j)    =  v(1) * q2
      d2aap(j+1)  =  v(2) * q2
      d2aap(l)    =  v(3) * q2
      d2aap(l+1)  =  v(4) * q2
      d2aap(n)    =  v(5) * q2
      d2aap(n+1)  =  v(6) * q2
!
    end do
!
20  close(aaap_unit,iostat=k)
!
    Return
  End subroutine ANTENNA
!
! -----------------------------------------------     CONVOLVE     -----
! This subroutine applies the convolution theorem
!
  Subroutine CONVOLVE ( RADIANCE, AAAP, M, IAS, IND, IERR )

    real(r4), intent(inout) :: RADIANCE(*)
    real(r4), intent(in) :: AAAP(*)
    integer(i4), intent(in) :: M, IAS, IND
    integer(i4), intent(out) :: IERR

    Integer, parameter :: MAXP=12, MAX2P=2**MAXP
!
    Integer(i4), save :: INIT = 0, MS = 0
    Integer(i4) :: ISTOP, M4, NTR, NTRH, I, J, K
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
    dblrad(:ntr) = dble(radiance(:ntr))
    if (init > 0 .and. init /= m) ms=0
    call rft1 ( dblrad, 'a', m4, ms, s )
    if(ms == -2) then
      init=0
      ierr=5
    else
      istop = min(ntrh,ias)
      do i = 2, istop
        j = 2*i-1
        k = 2*ias*(ind-1)+j
        cr = dblrad(j)*aaap(k)-dblrad(j+1)*aaap(k+1)
        ci = dblrad(j)*aaap(k+1)+dblrad(j+1)*aaap(k)
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
      k = 2*ias*(ind-1)+1
      dblrad(1) = dblrad(1)*aaap(k)
      dblrad(2) = dblrad(2)*aaap(k+1)
      call rft1 ( dblrad, 's', m4, ms, s )
      if(ms == -2) then
        init=0
        ierr = 6
      else
        init=m
        radiance(:ntr) = dblrad(:ntr)
      endif
    endif
!
    return
  End subroutine CONVOLVE
!
! -------------------------------------------------     FTGRID     -----
! This subroutine performs the interpolation onto the computational grid
! it uses cubic splines
!
  Subroutine FTGRID ( EIL_ANGLE, RADIANCE, DELTA0, XLAMDA, NP, NTR, IERR )

    Real(r4), intent(inout) :: EIL_ANGLE(*)
    Real(r4), intent(inout) :: RADIANCE(*)
    Real(r4), intent(in) :: DELTA0, XLAMDA
    Integer(i4), intent(in) :: NP, NTR
    Integer(i4), intent(out) :: IERR
    Integer(i4) :: I, K1, KN, N, MP
!
    Real(r4) :: X1, XN, R1,RN
    Real(r4) :: X(Nlvl), R(Nlvl)
    Real(r8) :: V, W, OOX, PP, DEL
!
    ierr = 0
    if (nlvl > np) then
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
      v = dble(xlamda)
      w = dble(delta0)
      oox = 1.0 / v
      pp = 0.185_r8 * oox
      del = oox / ntr
!
      do i = 1, ntr
        v = w - pp + del * (i - 1)             ! new code
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
!     CALL S_LINTRP(X, EIL_ANGLE(K1), R, RADIANCE(K1), MP, N)
      call cspline(x, eil_angle(k1:kn), r, radiance(k1:kn), mp, n)
!
    else
!
      ierr = 4
      Print *,'** Errro in FTGRID subroutine'
      Print *,'   MAXP < = NP (NP too big)'
!
    endif
!
    return
  End subroutine FTGRID
end module FOV_CONVOLVE_M

! $Log$
