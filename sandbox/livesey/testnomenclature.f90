PROGRAM Test
USE MLSCommon
USE MLSSignalNomenclature
USE QuantityTemplates
IMPLICIT NONE

TYPE (MLSSignalsDatabase_T) :: database
TYPE (MLSSignal_T), POINTER, DIMENSION(:) :: signalA,signalB,signalResult
INTEGER :: i
CHARACTER (LEN=120) :: request

PRINT*,MLSInstrumentModule_GHz
PRINT*,MLSInstrumentModule_THz
PRINT*,"Names:",MLSInstrumentModuleNames
PRINT*,"QuantityTypes:",QTYTypeNames

OPEN(UNIT=9,FILE="../../tables/emls-signals.dat",STATUS="OLD",ACTION="READ")
!OPEN(UNIT=9,FILE="umls-signals.dat",STATUS="OLD",ACTION="READ")
CALL ReadSignalsDatabase(9,database)
CLOSE(UNIT=9)

READ*,request
CALL ParseMLSSignalRequest(request,database,signalA)
DO i=1,SIZE(signalA)
   CALL GetFullMLSSignalName(signalA(i),request,showChannel=.TRUE.)
   PRINT*,'SignalA: '//TRIM(request)
ENDDO

READ*,request
CALL ParseMLSSignalRequest(request,database,signalB)
DO i=1,SIZE(signalB)
   CALL GetFullMLSSignalName(signalB(i),request,showChannel=.TRUE.)
   PRINT*,'SignalB: '//TRIM(request)
ENDDO

CALL IntersectionMLSSignalsInfo(signalA,signalB,signalResult)
PRINT*,"Done"

IF (ASSOCIATED(signalResult)) THEN
   DO i=1,SIZE(signalResult)
      CALL GetFullMLSSignalName(signalResult(i),request,showChannel=.TRUE.)
      PRINT*,'Intersection: '//TRIM(request)
   ENDDO
   CALL DestroyMLSSignalsInfo(signalResult)
ENDIF

CALL DestroyMLSSignalsInfo(signalA)
CALL DestroyMLSSignalsInfo(signalB)

CALL DestroySignalsDatabase(database)
END PROGRAM Test



