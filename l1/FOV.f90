! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE FOV   ! Field-Of-View routines/data
!=============================================================================

  USE MLSCommon, ONLY: R8
  USE Units, ONLY: Pi, Deg2Rad

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: InitFOVconsts, CalcMountsToFOV, GroundToFlightMountsGHz, &
       GroundToFlightMountsTHz, ScToGroundMountsGHz, ScToGroundMountsTHz 

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  ! Last dimension equals 2 where 1 = GHz and 2 = THz:

  REAL(r8) :: n0(3,2)
  REAL(r8) :: M(3,3,2)
  REAL(r8) :: S1(3,2)
  REAL(r8), PARAMETER :: H1(3,2) = RESHAPE ( &
       (/ -0.9047748d0,  0.0023903d0, -0.4258837d0, &
          -0.9047769d0, -0.0017294d0, -0.4258824d0 /), (/ 3, 2 /))

  REAL(r8), PARAMETER :: GroundToFlightMountsGHz(3,3) = RESHAPE ((/ &
       1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 /), (/ 3, 3 /))
  REAL(r8), PARAMETER :: GroundToFlightMountsTHz(3,3) = RESHAPE ((/ &
       1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 /), (/ 3, 3 /))
  REAL(r8), PARAMETER :: ScToGroundMountsGHz(3,3) = RESHAPE ((/ &
       1.000000, -0.000401, 0.000428, &
       0.000401, 1.00000, 0.000586, &
       -0.000428, -0.000586, 1.0000000 /), (/ 3, 3 /))
  REAL(r8), PARAMETER :: ScToGroundMountsTHz(3,3) = RESHAPE ((/ &
       1.000000, -0.000376, -0.000090, &
       0.000376, 1.00000, 0.000050, &
       0.000090, -0.000050, 1.0000000 /), (/ 3, 3 /))

CONTAINS

!=============================================================================
  SUBROUTINE InitFOVconsts   ! initialize FOV constants
!=============================================================================

    INTEGER :: i, n

! constants for each Module

    REAL(r8), PARAMETER :: dc(2) = (/ -0.001217d0, 0.000274d0 /)
    REAL(r8), PARAMETER :: thx(2) = (/ 0.000903d0, -0.001574d0 /)
    REAL(r8), PARAMETER :: thz(2) = (/ -0.002263d0, 0.000980d0 /)

! constants for each Rdm

    REAL(r8), PARAMETER :: dx(2) = (/ 0.000025d0, -0.000348d0 /)
    REAL(r8), PARAMETER :: dz(2) = (/ 0.000046d0, 0.001144d0 /)

! other constants

    REAL(r8), PARAMETER :: pi_4 = Pi / 4.0d0

! Initialize other constants

    DO n = 1, 2    ! GHz/THz

       n0(:,n) = (/ COS(pi_4 - dc(n)), SIN(pi_4 - dc(n)), 0.0d0 /)

       M(:,:,n) = RESHAPE ((/ COS(thz(n)), (COS(thx(n)) * SIN(thz(n))), &
            (SIN(thx(n)) * SIN(thz(n))), -SIN(thz(n)), &
            (COS(thx(n)) * COS(thz(n))), (SIN(thx(n)) * COS(thz(n))), 0.0d0, &
            -SIN(thx(n)),  COS(thx(n)) /), (/ 3, 3 /))

       DO i = 1, 3
          S1(i,n)  = DOT_PRODUCT (M(i,:,n), (/ dx(n), &
               -SQRT (1.0 - dx(n)**2 - dz(n)**2), dz(n) /))
       ENDDO
    ENDDO

  END SUBROUTINE InitFOVconsts

!=============================================================================
  SUBROUTINE CalcMountsToFOV (eps_deg, gtindx, MountsToFOV)
!=============================================================================

    REAL, INTENT (IN) :: eps_deg
    INTEGER, INTENT(IN) :: gtindx   ! 1 = GHz, 2 = THz
    REAL(r8), INTENT (OUT) :: MountsToFOV(3,3)

    INTEGER :: i
    REAL(r8) :: E2(3), eps
    REAL(r8) :: S(3,3), n(3), ns(3), h2(3), s2(3)

    eps = eps_deg * Deg2Rad   ! convert degrees to radians

    S = RESHAPE ((/ &
         COS(eps), 0.0d0, SIN(eps), 0.0d0, 1.0d0, 0.0d0, &
         -SIN(eps), 0.0d0, COS(eps) /), (/ 3, 3 /))

    DO i = 1, 3
       ns(i)  = DOT_PRODUCT (s(i,:), n0(:,gtindx))
    ENDDO

    DO i = 1, 3
       n(i)  = DOT_PRODUCT (M(i,:,gtindx), ns)
    ENDDO

    s2 = s1(:,gtindx) - 2 * DOT_PRODUCT (n, s1(:,gtindx)) * n
    h2 = h1(:,gtindx) - 2 * DOT_PRODUCT (n, h1(:,gtindx)) * n

    E2(1) = h2(2) * s2(3) - h2(3) * s2(2)
    E2(2) = h2(3) * s2(1) - h2(1) * s2(3)
    E2(3) = h2(1) * s2(2) - h2(2) * s2(1)

    MountsToFOV(1,:) = E2
    MountsToFOV(2,:) = H2
    MountsToFOV(3,:) = S2

  END SUBROUTINE CalcMountsToFOV

!=============================================================================
END MODULE FOV
!=============================================================================
! $Log$
! Revision 2.2  2004/11/10 15:36:53  perun
! Revise calculations per REC
!
! Revision 2.1  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
