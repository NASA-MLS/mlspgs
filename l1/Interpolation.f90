! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE Interpolation   ! Interpolation routines
!=============================================================================

  USE MLSL1Common, ONLY: r8

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------

  ! This module contains interpolation routines for the MLSL1 program

CONTAINS

  SUBROUTINE QuadInterpW (tVec, apoVec, qualVec, time, nVec, comVec, errmul, &
       status)

    !! Perform weighted quadratic interpolation

    IMPLICIT NONE

    INTEGER :: nVec
    REAL(r8) :: tVec(nVec)
    REAL(r8) :: apoVec(nVec)
    INTEGER :: qualVec(nVec)
    REAL(r8) :: time
    REAL(r8) :: comVec(nVec)
    REAL(r8) :: errmul
    INTEGER :: status

    !! ==================== LOCALS =========================

    INTEGER :: astat
    INTEGER, SAVE :: nsize = 0
    REAL(r8) :: sx1, sx2, sx3, sx4
    REAL(r8), ALLOCATABLE, DIMENSION(:),SAVE :: apoVec_x, tVec_x
    REAL(r8), ALLOCATABLE, DIMENSION(:), SAVE :: x1, x2, x3, x4

    status = 0   ! OK so far

    !! Check for minimum data size:

    IF (nVec <= 2) THEN
       status = 1    ! Not enough points
       RETURN
    ENDIF

    !! Check for minimum valid quality mask

    IF (SUM (qualVec) <= 3) THEN
       status = 2   ! Not enough valid data
       RETURN
    ENDIF

    !! Allocate local arrays:
    
    IF (nsize /= nVec) THEN
       DEALLOCATE (tVec_x, STAT=astat)
       ALLOCATE (tVec_x(nVec))
       DEALLOCATE (apoVec_x, STAT=astat)
       ALLOCATE (apoVec_x(nVec))
       DEALLOCATE (x1, STAT=astat)
       ALLOCATE (x1(nVec))
       DEALLOCATE (x2, STAT=astat)
       ALLOCATE (x2(nVec))
       DEALLOCATE (x3, STAT=astat)
       ALLOCATE (x3(nVec))
       DEALLOCATE (x4, STAT=astat)
       ALLOCATE (x4(nVec))
       nsize = nVec
    ENDIF

    tVec_x = tVec - time
    WHERE (apoVec == 1)
       apoVec_x = 1
    ELSEWHERE
       apoVec_x = 0
    ENDWHERE
    apoVec_x = apoVec * qualVec

    x1 = tVec_x * apoVec_x
    x2 = tVec_x * x1
    x3 = tVec_x * x2
    x4 = tVec_x * x3

    sx1 = SUM (x1)
    sx2 = SUM (x2)
    sx3 = SUM (x3)
    sx4 = SUM (x4)

    comVec = apoVec_x * (sx2 * sx4 - sx3 * sx3) - sx1 * (x1 * sx4 - x2 * sx3) &
         + sx2 * (x1 * sx3 - x2 * sx2)

    comVec = comVec / SUM (comVec)
    errmul = SUM (comVec * comVec)

!!$    print '(A/, 20(i4))', 'tVec: ', INT(tVec)
!!$    print  '(A/, 20(i4))', 'apoVec: ',  INT(apoVec)
!!$    print  '(A/, 20(i4))', 'qualVec: ', INT(qualVec)
!!$stop "interp"
  END SUBROUTINE QuadInterpW

!=============================================================================
END MODULE Interpolation
!=============================================================================

! $Log$
! Revision 2.2  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.1  2001/02/23 20:46:27  perun
! Version 0.5 commit
!
