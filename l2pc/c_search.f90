!
!------------------------------------------------------------------------
! This module searches an ordered 'c_array' for the 'str' nearest to that
! selected by the user. This is a character routine
!
      SUBROUTINE C_SEARCH(STR,C_ARRAY,NOE,LOW,HIGH)
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*), INTENT(IN) :: STR
      CHARACTER(LEN=*), INTENT(IN) :: C_ARRAY(*)
!
      INTEGER, INTENT(IN) :: NOE
      INTEGER, INTENT(INOUT) :: HIGH,LOW

      INTEGER :: MIDDLE,LS

      LS = LEN_TRIM(STR)
!
      IF(STR(1:LS) <= C_ARRAY(1)(1:LS)) THEN
        LOW = 1
        HIGH = 1
      ELSE IF(STR(1:LS) > C_ARRAY(NOE)(1:LS)) THEN
        LOW = NOE
        HIGH = NOE
      ELSE
        LOW = 1
        HIGH = NOE
        MIDDLE = (HIGH + LOW) / 2
        DO WHILE(HIGH-LOW > 1)
          IF(STR(1:LS) > C_ARRAY(MIDDLE)(1:LS)) THEN
            LOW = MIDDLE
          ELSE IF(STR(1:LS) < C_ARRAY(MIDDLE)(1:LS)) THEN
            HIGH = MIDDLE
          ELSE
            LOW = MIDDLE
            HIGH = MIDDLE
            RETURN
          END IF
          MIDDLE = (LOW + HIGH) / 2
        END DO
      END IF
!
      END SUBROUTINE C_SEARCH
