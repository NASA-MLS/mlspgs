PROGRAM TestUnits
USE MLSCommon
USE Units

CHARACTER(LEN=32) :: line
CHARACTER(LEN=32) :: remainder
INTEGER :: family
REAL(r8) :: value

!PRINT*,"Enter a request"
!READ*,line
line="10.0km 15.0 Hello there 57.3"
PRINT*,"I got:"//line//"|"

CALL ReadNumberWithUnitsFromString(line,value,family,remainder)
PRINT*,"Value: ",value
PRINT*,"Family: ", family
PRINT*,"Remainder: ",remainder
END PROGRAM TestUnits
