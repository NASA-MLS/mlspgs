! Copyright 2006, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
MODULE THzUtils
!=============================================================================

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: LLO_Bias, Bias_err, LLO_Label, ConvertLLO

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

  REAL, PARAMETER :: Bias_err = 10.0

  CHARACTER (LEN=*), PARAMETER, DIMENSION(16) :: LLO_Label = (/ &
       "Pump MC PZT      ", "Pump CC PZT      ", "FIR IC PZT       ", &
       "FIR OC PZT       ", "Pump Temperature ", "RF PS Temperature", &
       "FIR Temperature  ", "Pump Thermopile  ", "Bias Signal      ", &
       "RF Forward Power ", "Max FIR Power    ", "RFPA AGC         ", &
       "FIR Thermopile   ", "+5 V Supply      ", "+12 V Supply     ", &
       "-12 V Supply     " /)
  CHARACTER (LEN=5), PARAMETER :: ack7 = CHAR(250)//CHAR(250)//CHAR(3)// &
       CHAR(0)//CHAR(7)

CONTAINS

!=============================================================================
  FUNCTION LLO_Bias (LLO_DN) RESULT (Bias)
!=============================================================================

    CHARACTER (LEN=80), INTENT (IN) :: LLO_DN
    REAL :: Bias

    CHARACTER (LEN=80) :: DN
    CHARACTER (LEN=5), PARAMETER :: msg01 = CHAR(250)//CHAR(250)//CHAR(11)//&
         CHAR(1)//CHAR(28)
    CHARACTER (LEN=5), PARAMETER :: msg02 = CHAR(250)//CHAR(250)//CHAR(35)//&
         CHAR(2)//CHAR(7)
    INTEGER :: dindx

    DN = LLO_DN
    Bias = Bias_err

    dindx = 0
    IF (DN(1:5) == ack7) THEN
       dindx = 28
    ELSE
       IF (DN(7:11) == msg01 .AND. DN(16:16) == CHAR(1) .AND. &
            DN(27:31) == msg02) dindx = 48
    ENDIF

    IF (dindx > 0) Bias = (ICHAR(DN(dindx:dindx)) + 256 * &
         ICHAR(DN(dindx+1:dindx+1))) / 409.6

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

  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here
END MODULE THzUtils

! $Log$
! Revision 2.7  2006/04/05 18:09:12  perun
! Remove unused variables
!
! Revision 2.6  2006/03/24 15:20:53  perun
! Remove MIF 141 test for re-optimizing
!
! Revision 2.5  2005/06/23 18:41:36  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.4  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
! Revision 2.3  2004/01/09 17:46:23  perun
! Version 1.4 commit
!
! Revision 2.2  2003/09/15 17:15:54  perun
! Version 1.3 commit
!
! Revision 2.1  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
!
