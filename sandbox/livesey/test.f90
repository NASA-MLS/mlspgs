PROGRAM Test
USE MLSSignalNomenclature
IMPLICIT NONE

TYPE (MLSSignalsDatabase_T) :: database
TYPE (MLSSignal_T), POINTER, DIMENSION(:) :: signal
INTEGER :: i
CHARACTER (LEN=120) :: request

OPEN(UNIT=9,FILE="../../tables/emls-signals.dat",STATUS="OLD",ACTION="READ")
CALL ReadSignalsDatabase(9,database)
CLOSE(UNIT=9)

loop: DO
   READ*,request
   IF (TRIM(request)=="DONE") EXIT loop
   CALL ParseMLSSignalRequest(request,database,signal)
   DO i=LBOUND(signal,1),UBOUND(signal,1)
      CALL GetFullMLSSignalName(signal(i),request)
      PRINT*,'You mean: '//TRIM(request)
      PRINT*,signal(i)%channelIncluded
   ENDDO
END DO loop
END PROGRAM Test



