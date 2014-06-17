!***********************************************************************************************************************************
!  CROSS_PRODUCT
!
!  Returns the right-handed vector cross product of two 3-vectors:  C =
!  A x B.
!***********************************************************************************************************************************

SUBROUTINE CROSS_PRODUCT (A, B, C)

  IMPLICIT NONE

  DOUBLE PRECISION, DIMENSION(3), INTENT (IN)    :: A
  DOUBLE PRECISION, DIMENSION(3), INTENT (IN)    :: B
  DOUBLE PRECISION, DIMENSION(3), INTENT (OUT)   :: C

! result: 3-vector cross product

  C(1) = A(2)*B(3) - A(3)*B(2)
  C(2) = A(3)*B(1) - A(1)*B(3)
  C(3) = A(1)*B(2) - A(2)*B(1)

  RETURN

END SUBROUTINE CROSS_PRODUCT
