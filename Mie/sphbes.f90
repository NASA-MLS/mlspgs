!{Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.
!
! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module SPHBESS

  implicit NONE

  private

  public :: SPHBESS_Z

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine SPHBESS_Z ( X, RF, ACAP, JR, Y, J )
!{ Calculation of spherical Bessel functions
!
! {\tt J(n)} = $j_n( m \chi )$ is calculated for complex argument
!  $m \chi$ = {\tt RF * X} using the logarithmic derivative as first
!  done by Infeld (coding by J. V. Dave)
!
! {\tt JR(n)} = $j_n( \chi )$ is calculated for real argument
!  $\chi$ = {\tt X} using the logarithmic derivative
!
! {\tt Y(n)} = $y_n( \chi )$ is calculated by upward recurrence
!
! {\tt ACAP(n)} = $A_n(m \chi) = -\frac{n}{m \chi} +
!                 \frac{j_{n-1}(m \chi)}{j_n(m \chi)}$
! is calculated by backward recurrence using
! $A_{n-1}(m \chi) = 
! \frac{n}{m \chi} - \left( \frac{n}{m \chi} + A_n(m \chi) \right)^{-1}$.
!
! Adapted (with extensions) in MARCH 1988 from routine DBMIE of
!  J. V. Dave (1968) by C. D. Cantrell, 1988
!
! In particular, Dave used EXP(-IKR)/R for outgoing waves, so he
!  used the Hankel function of the second kind and a complex
!  refractive index with subtracted imaginary part. This routine
!  uses EXP(+IKR)/R for outgoing waves, hence Hankel functions
!  of the first kind and a complex refractive index with
!  added real part.
!
! Adapted from C. D. Cantrell by W. V. Snyder (2007), trivial
! H calculation removed.
!
! In a simplistic comparison with DBESJN and DBESYN from Math77 for
! $\log_{10}(x) = -1(.2)2$, n = 0..100, and RF=1.0, the maximum absolute
! difference for JR was 1.65E-15, the maximum absolute difference
! for REAL(J) was 5.47E-15, and the maximum relative difference for
! Y was 1.33E-12.
    use MLSKinds, only: R8
!   Arguments:
    real(r8), intent(in) :: X             ! Argument \chi of JR and Y, X > 0
    complex(r8), intent(in) :: RF         ! Index of refraction m
    real(r8), intent(out) :: Y(0:)        ! Neumann (Bessel 2nd kind) of X
    real(r8), intent(out) :: JR(0:)       ! Bessel 1st kind of X
    complex(r8), intent(out) :: ACAP (0:) ! A(rf*x)
    complex(r8), intent(out), optional :: J(0:) ! Bessel 1st kind of RF*X
!   Local variables
    real(r8) :: ABSRF                     ! ABS(rf)
    real(r8), save, allocatable :: ACAPR(:)     ! A(x)
    complex(r8), save, allocatable :: ACAPW(:)  ! A(rf*x)
    real(r8) :: Cut                       ! Maximum value for YM1
    real(r8) :: JRM1, RX, YM1, YM2        ! JR_{n-1}, 1/x, Y_{n-1}, Y_{n-2}
    real(r8) :: RXN                       ! rx * ( 2n-1 )
    complex(r8) :: JNM1, RRFX             ! J_{n-1}, 1/(rf*x)
    integer :: MAXN, NMX1, N              ! ubound(y), size(work), loop inductor
!*****Complex index of refraction & its reciprocal
    absrf = abs(rf)
    rx = 1.0_r8 / x
    rrfx = rx * conjg(rf) / absrf**2 ! 1/(rf*x)
!*****Calculate upper limit of order for downward recurrence
    maxn = ubound(y,1)
    nmx1 = max(maxn,150,int(1.1_r8 * absrf * x)) + 1
    if ( allocated(acapr) ) then
      if ( ubound(acapr,1) < nmx1 ) deallocate (acapr, acapw )
    end if
    if ( .not. allocated(acapr) ) allocate ( acapr(0:nmx1), acapw(0:nmx1) )
    cut = 0.25 * huge(ym1) / (maxn * max(rx,1.0_r8))
!{ Use downward recurrence for the logarithmic derivative of J(RF)
!  and JR(X)
!
!  {\tt acapw(n-1)} = $A_{n-1}(m \chi) =
!    \frac{n}{m \chi} - \left( \frac{n}{m \chi} + A_n(m \chi) \right)^{-1}$
!
!  {\tt Acapr(n-1)} = $A_{n-1}(\chi) =
!    \frac{n}{\chi} - \left( \frac{n}{\chi} + A_n(\chi) \right)^{-1}$
!
!  Start with the asymptotic approximation for fixed $z$ as $n$ increases:
!  $A_n(z) \sim \frac{n+1}z$
    acapw(nmx1) = (nmx1+1) * rrfx
    acapr(nmx1) = (nmx1+1) * rx
    do n = nmx1, 1, -1
      acapw(n-1) = n * rrfx - 1.0_r8 / ( n * rrfx + acapw(n) )
      acapr(n-1) = n * rx -   1.0_r8 / ( n * rx +   acapr(n) )
    end do
    acap(:maxn) = acapw(:maxn)
!{ Calculate J(RF) and JR(X) by upward recurrence using the logarithmic
!  derivative
!
! {\tt J(0)} = $j_0(m \chi) = \frac{\sin(m \chi)}{m \chi}$,
! {\tt J(n)} = $j_n(m \chi) = \frac{j_{n-1}(m \chi)}{A_n(m \chi)+\frac{n}{m \chi}}$
!
! {\tt JR(0)} = $j_0(\chi) = \frac{\sin(\chi)}{\chi}$,
! {\tt JR(n)} = $j_n(\chi) = \frac{j_{n-1}(\chi)}{A_n(\chi)+\frac{n}{\chi}}$
    jr(0) = sin (x) * rx    ! sin(x)/x
    jrm1 = jr(0)
    do n = 1, maxn
      jrm1 = jrm1 / (acapr(n) + n * rx)
      jr(n) = jrm1
    end do
    if ( present(j) ) then
      j(0) = sin(rf*x) * rrfx ! sin(rf*x)/(rf*x)
      jnm1 = j(0)
      do n = 1, maxn
        jnm1 = jnm1 / (acapw(n) + n * rrfx)
        j(n) = jnm1
      end do
    end if
!{ Calculate Neumann functions {\tt Y(n)} = $y_n(x)$ for real argument
!  by upward recurrence
!
!  $y_{-1}(x) = \frac{\sin(x)}x$, $y_0(x) = -\frac{\cos(x)}x$,
!  $y_n(x) = \frac{2 n - 1}x y_{n-1}(x) - y_{n-2}(x)$
    ym2 = sin (x) * rx
    ym1 = - cos (x) * rx
    y(0) = ym1
    do n = 1, maxn
      rxn = (2 * n - 1) * rx
      if ( rxn > 1 ) then
        if ( abs(ym1) > cut ) then ! Avoid overflow
          y(n:maxn) = sign(huge(ym1),ym1)
          exit
        end if
      end if
      y(n) = rxn * ym1 - ym2
      ym2 = ym1
      ym1 = y(n)
    end do

  end subroutine SPHBESS_Z

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module SPHBESS

! $Log$
! Revision 1.1  2008/04/19 01:15:27  vsnyder
! Initial commit
!
