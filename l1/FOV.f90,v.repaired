! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
MODULE FOV   ! Field-Of-View routines/data
!=============================================================================

  USE Constants, ONLY: Pi, Deg2Rad
  USE MLSKinds, ONLY: R8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: InitFOVconsts, CalcMountsToFOV, GroundToFlightMountsGHz, &
       GroundToFlightMountsTHz, ScToGroundMountsGHz, ScToGroundMountsTHz 

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! Last dimension equals 2 where 1 = GHz and 2 = THz:

  REAL(r8) :: n0(3,2)
  REAL(r8) :: M(3,3,2)
  REAL(r8) :: S1(3,2)

  REAL(r8), PARAMETER :: H1(3,2) = RESHAPE ( &
       (/ -0.9047748d0,  0.0023903d0, -0.4258837d0, &
          -0.9047769d0, -0.0017294d0, -0.4258824d0 /), (/ 3, 2 /))

  REAL(r8), PARAMETER :: GroundToFlightMountsGHz(3,3) = &
       & RESHAPE ((/ 1.0, 0.0, 0.0, &
       &             0.0, 1.0, 0.0, &
       &             0.0, 0.0, 1.0 /), (/ 3, 3 /))

  REAL(r8), PARAMETER :: GroundToFlightMountsTHz(3,3) = &
       & RESHAPE ((/ 1.0, 0.0, 0.0, &
       &             0.0, 1.0, 0.0, &
       &             0.0, 0.0, 1.0 /), (/ 3, 3 /))

  !<whd> 
  ! 
  ! Construct some rotation matrices that handle the small rotations
  ! found while doing integration and test.
  !
  ! </whd>

  REAL(r8), PARAMETER :: ScToGroundMountsGHz(3,3) = RESHAPE ((/ &
       1.000000,  -0.000401, 0.000428, &
       0.000401,   1.00000,  0.000586, &
      -0.000428,  -0.000586, 1.0000000 /), (/ 3, 3 /))

  ! <whd> Same, same for THz antenna </whd>

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
! <whd> 1 for Ghz, 2 for THz </whd>

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

       M(:,:,n) = RESHAPE ((/ &
            & COS(thz(n)), (COS(thx(n)) * SIN(thz(n))), (SIN(thx(n)) * SIN(thz(n))), &
            & -SIN(thz(n)),(COS(thx(n)) * COS(thz(n))), (SIN(thx(n)) * COS(thz(n))), &
            &  0.0d0,                -SIN(thx(n)),              COS(thx(n)) /), &
            &  (/ 3, 3 /))

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

    !<whd> FOV is rotated eps degress about the mounts 'y' axis</whd>
    S = RESHAPE ((/ &
         &  COS(eps),  0.0d0,    SIN(eps), &
         &   0.0d0,    1.0d0,      0.0d0, &
         &  -SIN(eps), 0.0d0,    COS(eps) /), &
             & (/ 3, 3 /))

    DO i = 1, 3
       ns(i)  = DOT_PRODUCT (s(i,:), n0(:,gtindx))
    ENDDO

    DO i = 1, 3
       n(i)  = DOT_PRODUCT (M(i,:,gtindx), ns)
    ENDDO

    !<whd> Snell's law: S = look vector, H is magnetic vector, n the
    !normal to the plane.</whd>
    s2 = s1(:,gtindx) - 2 * DOT_PRODUCT (n, s1(:,gtindx)) * n
    h2 = h1(:,gtindx) - 2 * DOT_PRODUCT (n, h1(:,gtindx)) * n

    ! <whd> E (electric field vector?) = H cross S </whd>
    E2(1) = h2(2) * s2(3) - h2(3) * s2(2)
    E2(2) = h2(3) * s2(1) - h2(1) * s2(3)
    E2(3) = h2(1) * s2(2) - h2(2) * s2(1)

    MountsToFOV(1,:) = E2
    MountsToFOV(2,:) = H2
    MountsToFOV(3,:) = S2

  END SUBROUTINE CalcMountsToFOV

!=============================================================================

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here
END MODULE FOV
!=============================================================================
! $Log$
! Revision 2.5  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.4.6.1  2015/10/09 10:21:38  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.4  2009/05/13 20:33:05  vsnyder
! Get constants from Constants, kinds from MLSKinds
!
! Revision 2.3  2005/06/23 18:41:35  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.2  2004/11/10 15:36:53  perun
! Revise calculations per REC
!
! Revision 2.1  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
