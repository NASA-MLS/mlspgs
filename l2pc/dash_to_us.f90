!
!---------------------------------------------------------------------
!
      SUBROUTINE DASH_TO_US(A)
!
      IMPLICIT NONE

      CHARACTER(LEN=*),INTENT(INOUT) :: A
!
      INTEGER :: I,J

      J = LEN_TRIM(A)
      DO I = 1, J
        IF(A(I:I) == '-') A(I:I)='_'
      END DO
!
      END SUBROUTINE DASH_TO_US
