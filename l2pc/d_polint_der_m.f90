module D_POLINT_DER_M
  use D_HUNT_M, only: HUNT
  implicit NONE
  private
  public D_POLINT_DER, POLINT_DER
  interface POLINT_DER; module procedure D_POLINT_DER; end interface
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
    "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
    "$RCSfile$"
!---------------------------------------------------------------------------
  integer, private, parameter :: RK = kind(0.0d0)
contains
  subroutine D_POLINT_DER (XX,YY,N,X,Y,DY,D2Y)
!
!   include 'polint_der.f9h'
!
      INTEGER, INTENT(IN) :: N
      REAL(rk), INTENT(IN) :: XX(*),YY(*),X
      REAL(rk), INTENT(OUT) :: Y,DY,D2Y

      INTEGER :: I,J
      REAL(rk) :: DX,YM,YP,D
!
      CALL POLINT(XX,YY,N,X,Y)

      J = 2
      DX = ABS(XX(N)-XX(1))
      DO I = 1, N
        IF(ABS(XX(I)-X) < DX) THEN
          J = I
          DX = ABS(XX(I)-X)
        ENDIF
      END DO

      IF(J == 1) THEN
        D = ABS(XX(2) - XX(1))
      ELSE IF(J == N) THEN
        D = ABS(XX(N) - XX(N-1))
      ELSE
        D = 0.5D0 * ABS(XX(J+1) - XX(J-1))
      ENDIF

      DX = MAX(1.0D-12, 0.05D0 * D)

      CALL POLINT(XX,YY,N,X+DX,YP)
      CALL POLINT(XX,YY,N,X-DX,YM)
      DY = (YP - YM) / (2.0d0 * DX)
      D2Y = (YP - 2.0d0 * Y + YM) / (DX*DX)
!
  end subroutine D_POLINT_DER
!
!---------------------------------------------------------------------
!
  SUBROUTINE POLINT(XA,YA,N,X,Y)
!
      INTEGER, INTENT(IN) :: N
      REAL(rk), INTENT(IN) :: XA(*),YA(*),X
      REAL(rk), INTENT(OUT) :: Y

      INTEGER :: I,M,K,NS,ILO,KNPT,ILOC(1)
      REAL(rk) :: DIF,DIFT,HO,HP,W,DEN,DY

      INTEGER,PARAMETER :: NDEG = 12
      REAL(rk), DIMENSION(NDEG) :: TMPX,TMPY,C,D
!
! Starts processing:
!
      CALL HUNT(X,XA,N,ILO,I)
      I = NDEG / 2
      M = MAX(1, ILO-I)
      K = MIN(M+NDEG-1, N)
      KNPT = K - M + 1

      TMPX(1:KNPT) = XA(M:K)
      TMPY(1:KNPT) = YA(M:K)
!
      C(1:KNPT) = X
      ILOC = MINLOC(ABS(C-TMPX))
      NS = ILOC(1)

      C = TMPY
      D = TMPY
!
      Y = TMPY(NS)
      NS = NS - 1
      DO M = 1,KNPT-1
        DO I = 1,KNPT-M
          HO = TMPX(I)-X
          HP = TMPX(I+M)-X
          W = C(I+1)-D(I)
          DEN = HO-HP
          IF(DEN == 0.0_rk) THEN
            PRINT *,'** BAD XA ARRAY, ?_POLINT BOMBED OUT ..'
            RETURN
          ENDIF
          DEN = W / DEN
          D(I) = HP*DEN
          C(I) = HO*DEN
        END DO
        IF (2*NS < KNPT-M) THEN
          DY = C(NS+1)
        ELSE
          DY = D(NS)
          NS = NS-1
        ENDIF
        Y = Y+DY
      END DO
!
  END SUBROUTINE POLINT
!
end module D_POLINT_DER_M
! $Log$
! Revision 1.1  2000/05/04 18:12:04  Z.Shippony
! Initial conversion to Fortran 90
!
