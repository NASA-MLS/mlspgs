! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE THzUtils
!=============================================================================

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: LLO_Bias, Bias_err, MaxBias

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  REAL, PARAMETER :: Bias_err = 10.0
  REAL, PARAMETER :: MaxBias = 2.0

CONTAINS

!=============================================================================
  FUNCTION LLO_Bias (LLO_DN, MIF) RESULT (Bias)
!=============================================================================

    CHARACTER (LEN=80), INTENT (IN) :: LLO_DN
    INTEGER, INTENT (IN) :: MIF
    REAL :: Bias

    CHARACTER (LEN=80) :: DN
    CHARACTER (LEN=5), PARAMETER :: ack7 = CHAR(250)//CHAR(250)//CHAR(3)//&
         CHAR(0)//CHAR(7)
    CHARACTER (LEN=5), PARAMETER :: msg01 = CHAR(250)//CHAR(250)//CHAR(11)//&
         CHAR(1)//CHAR(28)
    CHARACTER (LEN=5), PARAMETER :: msg02 = CHAR(250)//CHAR(250)//CHAR(35)//&
         CHAR(2)//CHAR(7)

    DN = LLO_DN
    Bias = Bias_err

    IF (MIF == 141) THEN    ! Could be re-optimizing

       IF (DN(7:11) == msg01 .AND. DN(16:16) == CHAR(1) .AND. DN(27:31) == &
            msg02) DN = DN(21:80)

    ENDIF
    IF (DN(1:5) == ack7) THEN

       Bias = (ICHAR(DN(28:28)) + 256 * ICHAR(DN(29:29))) / 409.6

    ENDIF

  END FUNCTION LLO_Bias

END MODULE THzUtils

! $Log$
! Revision 2.1  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
!
