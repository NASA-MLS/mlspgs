PROGRAM Test
USE SignalsFile
IMPLICIT NONE

TYPE (MLS_SignalsDatabase_T) :: database
INTEGER :: i

OPEN(UNIT=9,FILE="emls-signals.dat",STATUS="OLD",ACTION="READ")
CALL ReadSignalsDatabase(9,database)
CLOSE(UNIT=9)
DO i=0,database%noRadiometers
   PRINT*,'--------------------'
   PRINT*,database%radiometerInfo(i)%name
   PRINT*,database%radiometerInfo(i)%prefix
   PRINT*,database%radiometerInfo(i)%suffix
   PRINT*,database%radiometerInfo(i)%number
   PRINT*,database%radiometerInfo(i)%modifier
ENDDO

END PROGRAM Test



