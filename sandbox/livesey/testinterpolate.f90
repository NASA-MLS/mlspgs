PROGRAM TestInterpolate
USE MLSNumerics

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(5) :: oldX=(/1,2,3,4,5/)
DOUBLE PRECISION, DIMENSION(5) :: oldY=(/1,3,9,12,10/)
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: newX
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: newY,dyByDX

INTEGER :: I,noValues

PRINT*,"Hello"

PRINT*,"Enter number of values"
READ*,noValues

ALLOCATE(newX(noValues))
ALLOCATE(newY(noValues),dyByDx(noValues))

PRINT*,"Enter the values"
READ*,newX

PRINT*,"Calling"
CALL InterpolateValues(oldX,oldY,newX,newY,'Spline',extrapolate="Constant", &
     & dyByDx=dyByDx)
PRINT*,"Done"

PRINT*,newY
PRINT*,dyByDx
END PROGRAM TestInterpolate

