! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE THzUtils
!=============================================================================

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: LLO_Bias, Bias_err, MaxBias, LLO_Label, ConvertLLO

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  REAL, PARAMETER :: Bias_err = 10.0
  REAL, PARAMETER :: MaxBias = 2.0

  CHARACTER (LEN=*), PARAMETER, DIMENSION(16) :: LLO_Label = (/ &
       "Pump MC PZT      ", "Pump CC PZT      ", "FIR IC PZT       ", &
       "FIR OC PZT       ", "Pump Temperature ", "RF PS Temperature", &
       "FIR Temperature  ", "Pump Thermopile  ", "Bias Signal      ", &
       "RF Forward Power ", "Max FIR Power    ", "RFPA AGC         ", &
       "FIR Thermopile   ", "+5 V Supply      ", "+12 V Supply     ", &
       "-12 V Supply     " /)

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

!=============================================================================
  SUBROUTINE ConvertLLO (DN, EU)
!=============================================================================

    USE MLSL1Utils, ONLY: BigEndianStr, QNan

    CHARACTER (LEN=80), INTENT (IN) :: DN
    REAL, INTENT (OUT) :: EU(16)

! coefficients:

    REAL, DIMENSION(16), PARAMETER :: c0 = (/ &
         0.0, 0.0, 0.0, 0.0, -23.2, -19.655, -23.2, -0.1124, &
         0.0, 0.0, 0.0, -0.1124, -0.02, 0.0, 0.0, -30.0 /)
    REAL, DIMENSION(16), PARAMETER :: c1 = (/ &
         0.1, 0.1, 0.1, 0.1, 16.155, 18.414, 16.155, 3.5636, &
         1.0, 1.0, 1.0, 3.5636, 6.84, 1.0, 2.0, 6.0 /)

    CHARACTER (LEN=2) :: cdat
    CHARACTER (LEN=5), PARAMETER :: ack7 = CHAR(250)//CHAR(250)//CHAR(3)//&
         CHAR(0)//CHAR(7)

    INTEGER :: i, n1
    REAL :: xdn

    EU = QNan()
    IF (DN(1:5) == ack7) THEN

       DO i = 1, 16
          n1 = (i - 1) * 2 + 1
          cdat(1:1) = DN(n1+12:n1+12)
          cdat(2:2) = DN(n1+11:n1+11)
          IF (i == 11) THEN
             EU(i) = BigEndianStr (cdat)
          ELSE IF (i == 12) THEN
             EU(i) = BigEndianStr (DN(n1+11:n1+12))
          ELSE
             xdn = BigEndianStr (cdat)
             IF (i /= 11) xdn = xdn / 409.6
             xdn = c0(i) + c1(i) * xdn
             IF (i /= 6) THEN
                EU(i) = xdn
             ELSE
                EU(i) = xdn + xdn * xdn * (-5.606) + &
                     xdn * xdn * xdn * 1.509
             ENDIF
          ENDIF
       ENDDO

    ENDIF

  END SUBROUTINE ConvertLLO

END MODULE THzUtils

! $Log$
! Revision 2.2  2003/09/15 17:15:54  perun
! Version 1.3 commit
!
! Revision 2.1  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
!
