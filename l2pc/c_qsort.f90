!
!---------------------------------------------------------------------
!
      SUBROUTINE C_QSORT(N,CHRY)

      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: N

      CHARACTER(LEN=*), INTENT(INOUT) :: CHRY(*)

      INTEGER :: L,IR,I,J
      CHARACTER(LEN=40) :: CA
!
      L = N / 2 + 1
      IR = N
!
      DO WHILE(IR > 1)
        IF(L > 1) THEN
          L = L- 1
          CA = CHRY(L)
        ELSE
          CA = CHRY(IR)
          CHRY(IR) = CHRY(1)
          IR = IR - 1
          IF(IR == 1) THEN
            CHRY(1) = CA
            RETURN
          ENDIF
        ENDIF
        I = L
        J = L + L
        DO WHILE(J <= IR)
          IF(J < IR) THEN
            IF(CHRY(J) < CHRY(J+1)) J = J + 1
          ENDIF
          IF(CA < CHRY(J)) THEN
            CHRY(I) = CHRY(J)
            I = J
            J = J + J
          ELSE
            J = IR + 1
          ENDIF
        END DO
        CHRY(I) = CA
      END DO
!
      END SUBROUTINE C_QSORT
