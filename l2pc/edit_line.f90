!
!------------------------------------------------------------------
!
      SUBROUTINE EDIT_LINE(LINE,L)
!
      IMPLICIT NONE
!
!  Edit the line to eliminate -. or .
!
      CHARACTER(LEN=*), INTENT(INOUT) :: LINE
!
      INTEGER :: L,I
!
      L = LEN_TRIM(LINE)
      IF(L < 1) RETURN
!
      I = INDEX(LINE,' .')
      DO WHILE(I > 0)
        LINE(I:I) = '0'
        I = INDEX(LINE,' .')
      END DO
!
      I = INDEX(LINE,'-.')
      DO WHILE(I > 0)
        LINE(I-1:I+1) = '-0.'
        I = INDEX(LINE,'-.')
      END DO
!
      END SUBROUTINE EDIT_LINE
!
!------------------------------------------------------------------
!
      SUBROUTINE EDIT_LINE_MOD(LINE,L)
!
!  Edit the line to eliminate -. or .
!
      INTEGER, INTENT(INOUT) :: L
      CHARACTER(LEN=*), INTENT(INOUT) :: LINE
!
      INTEGER :: I,J
      CHARACTER(LEN=1) :: CH
      CHARACTER(LEN=254) :: BUFF
!
      L = LEN_TRIM(LINE)
      IF(L < 1) RETURN
!
      IF(LINE(L:L) == '.') THEN
        J = 1
        IF(L > 1) THEN
          CH = LINE(L-1:L-1)
          J = INDEX('-0123456789',CH)
        ENDIF
        IF(J > 0) THEN
          L = L + 1
          LINE(L:L) = '0'
        ENDIF
      ENDIF
!
      I = INDEX(LINE,'. ')
      DO WHILE(I > 0)
        CH = LINE(I-1:I-1)
        J = INDEX('0123456789',CH)
        IF(J > 0) LINE(I+1:I+1) = '0'
        I = INDEX(LINE,'. ')
      END DO
!
      I = INDEX(LINE,'  .')
      DO WHILE(I > 0)
        LINE(I:I+1)=' 0'
        I = INDEX(LINE,'  .')
      END DO
!
      I = INDEX(LINE,' .')
      DO WHILE(I > 1)
        BUFF(1:)=' '
        BUFF(1:I-1)=LINE(1:I-1)
        BUFF(I+1:I+2)='0.'
        BUFF(I+3:L+1) = LINE(I+2:L)
        LINE = BUFF
        L = L + 1
        I = INDEX(LINE,' .')
      END DO
!
      I = INDEX(LINE,' -.')
      DO WHILE(I > 0)
        LINE(I:I+2)='-0.'
        I = INDEX(LINE,' -.')
      END DO
!
      I = INDEX(LINE,'-.')
      DO WHILE(I > 0)
        BUFF(1:)=' '
        BUFF(1:I)=LINE(1:I)
        BUFF(I+1:I+1)='0'
        BUFF(I+2:L+1) = LINE(I+1:L)
        LINE = BUFF
        L = L + 1
        I = INDEX(LINE,'-.')
      END DO
!
      END SUBROUTINE EDIT_LINE_MOD
