!
!----------------------------------------------------------------------
!  This routine reverse the order in a chracter array
!
      SUBROUTINE C_REVERSE(N,ARY)

      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: N
      CHARACTER(LEN=*), INTENT(INOUT) :: ARY(*)

      INTEGER :: I,K,L
      CHARACTER(LEN=78) :: CA
!
      IF(N < 2) RETURN
!
      DO I = 1, N/2
        CA(1:)=' '
        K = N-I+1
        CA = ARY(K)
        L = LEN_TRIM(CA)
        ARY(K) = ARY(I)
        ARY(I) = CA(1:L)
      END DO
!
      END SUBROUTINE C_REVERSE
