PROGRAM TestHunt
USE MLSNumerics

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(50) :: list
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: values
INTEGER, DIMENSION(:), ALLOCATABLE :: indices

INTEGER :: I,noValues

DO I=1,SIZE(list)
   list(I)=I
ENDDO

PRINT*,"Enter number of values"
READ*,noValues

ALLOCATE(values(noValues))
ALLOCATE(indices(noValues))

PRINT*,"Enter the values"
READ*,values

PRINT*,"Calling"
CALL Hunt(list,values,indices)
PRINT*,"Done"

PRINT*,indices
END
