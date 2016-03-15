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
MODULE Interpolation   ! Interpolation routines
!=============================================================================

  USE MLSL1Common, ONLY: r8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: QuadInterpW

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! This module contains interpolation routines for the MLSL1 program

CONTAINS

  SUBROUTINE QuadInterpW (tVec, &    ! in, time vector
       &                  apoVec, &  ! in, weight vector
       &                  qualVec, & ! in, quality vector
       &                  time, &    ! in, of interpolation
       &                  nVec, &    ! in, number of vectors
       &                  comVec, &  ! out
       &                  errmul, &  ! out
       &                  status)    ! out

    !! Perform weighted quadratic interpolation

    IMPLICIT NONE

    INTEGER,  INTENT(in)  :: nVec          ! dimension of other arrays
    REAL(r8), INTENT(in)  :: tVec(nVec)    ! time, as MIF index
    REAL(r8), INTENT(in)  :: apoVec(nVec)  ! weight
    INTEGER,  INTENT(in)  :: qualVec(nVec) ! quality
    REAL(r8), INTENT(in)  :: time          ! The number of the MIF
    REAL(r8), INTENT(out) :: comVec(nVec) 
    REAL(r8), INTENT(out) :: errmul
    INTEGER,  INTENT(out) :: status

    !! ==================== LOCALS =========================

    REAL(r8) :: sx1, sx2, sx3, sx4
    REAL(r8) :: apoVec_x(nVec), tVec_x(nVec)
    REAL(r8) :: x1(nVec), x2(nVec), x3(nVec), x4(nVec)

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

    tVec_x = tVec - time
    apoVec_x = apoVec * qualVec

    x1 = tVec_x * apoVec_x
    x2 = tVec_x * x1
    x3 = tVec_x * x2
    x4 = tVec_x * x3

    sx1 = SUM (x1)
    sx2 = SUM (x2)
    sx3 = SUM (x3)
    sx4 = SUM (x4)

    ! See equation D.4 in the Level 1 Algorithm Theoretical Basis document
    comVec = apoVec_x * (sx2 * sx4  -  sx3 * sx3) - &
         & sx1 * (x1 * sx4  -  x2 * sx3) + &
         & sx2 * (x1 * sx3  -  x2 * sx2)

    comVec = comVec / SUM (comVec)
    errmul = SUM (comVec * comVec)

  END SUBROUTINE QuadInterpW

!=============================================================================
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here
END MODULE Interpolation
!=============================================================================

! $Log$
! Revision 2.7  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.6.6.1  2015/10/09 10:21:38  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.6  2005/06/23 18:41:35  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.5  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
! Revision 2.4  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.3  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
! Revision 2.2  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.1  2001/02/23 20:46:27  perun
! Version 0.5 commit
!
