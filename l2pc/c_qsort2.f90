!
!---------------------------------------------------------------------
!
      SUBROUTINE C_QSORT2(N,CHRY,IREC)

      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(OUT) :: IREC(*)

      CHARACTER(LEN=*), INTENT(INOUT) :: CHRY(*)
!
      CHARACTER(LEN=40) :: CA
      INTEGER :: I,J,L,JA,IR
!
      L = N / 2 + 1
      IR = N
!
      DO WHILE(IR > 1)
        IF(L > 1) THEN
          L = L - 1
          CA = CHRY(L)
          JA = IREC(L)
        ELSE
          CA = CHRY(IR)
          JA = IREC(IR)
          CHRY(IR) = CHRY(1)
          IREC(IR) = IREC(1)
          IR = IR - 1
          IF(IR == 1) THEN
            CHRY(1) = CA
            IREC(1) = JA
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
            IREC(I) = IREC(J)
            I = J
            J = J + J
          ELSE
            J = IR + 1
          ENDIF
        END DO
        CHRY(I) = CA
        IREC(I) = JA
      END DO
!
      END SUBROUTINE C_QSORT2
