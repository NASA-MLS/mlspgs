! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE FOV   ! Field-Of-View routines/data
!=============================================================================

  USE MLSCommon, ONLY: R8
  USE Units, ONLY: Pi, Deg2Rad

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: InitFOVconsts, CalcMountsToFOV

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  REAL(r8) :: n0(3)
  REAL(r8) :: M(3,3)
  REAL(r8) :: S1(3)
  REAL(r8), PARAMETER :: H1(3) = (/ -0.904784d0, 0.002174d0, -0.425866d0 /)

CONTAINS

!=============================================================================
  SUBROUTINE InitFOVconsts   ! initialize FOV constants
!=============================================================================

    INTEGER :: i

! constants for GHz Module

    REAL(r8), PARAMETER :: dc = -0.001181d0
    REAL(r8), PARAMETER :: thx = 0.000852d0
    REAL(r8), PARAMETER :: thz = -0.002201d0

! constants for each GHz Rdm

    REAL(r8), PARAMETER :: dx = 0.000012d0
    REAL(r8), PARAMETER :: dz = 0.000399d0

! other constants

    REAL(r8), PARAMETER :: pi_4 = Pi / 4.0d0

! Initialize other constants

    n0 = (/ COS(pi_4 - dc), SIN(pi_4 - dc), 0.0d0 /)

    M = RESHAPE ((/ COS(thz), (COS(thx) * SIN(thz)), (SIN(thx) * SIN(thz)),&
         -SIN(thz), (COS(thx) * COS(thz)), (SIN(thx) * COS(thz)), &
         0.0d0, -SIN(thx),  COS(thx) /), (/ 3, 3 /))

    DO i = 1, 3
       S1(i)  = DOT_PRODUCT (M(i,:), (/ dx, -SQRT (1.0 - dx**2 - dz**2), dz /))
    ENDDO

  END SUBROUTINE InitFOVconsts

!=============================================================================
  SUBROUTINE CalcMountsToFOV (eps_deg, MountsToFOV)
!=============================================================================

    REAL, INTENT (IN) :: eps_deg
    REAL(r8), INTENT (OUT) :: MountsToFOV(3,3)

    INTEGER :: i
    REAL(r8) :: E2(3), eps
    REAL(r8) :: S(3,3), n(3), ns(3), h2(3), s2(3)

    eps = eps_deg * Deg2Rad   ! convert degrees to radians

    S = RESHAPE ((/ &
         COS(eps), 0.0d0, SIN(eps), 0.0d0, 1.0d0, 0.0d0, &
         -SIN(eps), 0.0d0, COS(eps) /), (/ 3, 3 /))

    DO i = 1, 3
       ns(i)  = DOT_PRODUCT (s(i,:), n0)
    ENDDO

    DO i = 1, 3
       n(i)  = DOT_PRODUCT (M(i,:), ns)
    ENDDO

    s2 = s1 - 2 * DOT_PRODUCT (n, s1) * n
    h2 = h1 - 2 * DOT_PRODUCT (n, h1) * n

    E2(1) = h2(2) * s2(3) - h2(3) * s2(2)
    E2(2) = h2(3) * s2(1) - h2(1) * s2(3)
    E2(3) = h2(1) * s2(2) - h2(2) * s2(1)

!!$    MountsToFOV(:,1) = E2
!!$    MountsToFOV(:,2) = H2
!!$    MountsToFOV(:,3) = S2

    MountsToFOV(1,:) = E2
    MountsToFOV(2,:) = H2
    MountsToFOV(3,:) = S2

  END SUBROUTINE CalcMountsToFOV

!=============================================================================
END MODULE FOV
!=============================================================================
! $Log$
! Revision 2.1  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
