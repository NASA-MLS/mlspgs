PROGRAM Test
USE MLSSignalNomenclature
IMPLICIT NONE

TYPE (MLSSignalsDatabase_T) :: database
TYPE (MLSSignal_T), POINTER, DIMENSION(:) :: signal
INTEGER :: i
CHARACTER (LEN=120) :: request

OPEN(UNIT=9,FILE="../../tables/emls-signals.dat",STATUS="OLD",ACTION="READ")
!OPEN(UNIT=9,FILE="umls-signals.dat",STATUS="OLD",ACTION="READ")
CALL ReadSignalsDatabase(9,database)
CLOSE(UNIT=9)
! DO i=1,database%noRadiometers
!    PRINT*,'--------------------'
!    PRINT*,database%radiometerInfo(i)%name
!    PRINT*,database%radiometerInfo(i)%lo
!    PRINT*,database%radiometerInfo(i)%prefix
!    PRINT*,database%radiometerInfo(i)%suffix
!    PRINT*,database%radiometerInfo(i)%number
! ENDDO

! DO i=1,database%noBands
!    PRINT*,'--------------------'
!    PRINT*,database%bandInfo(i)%name
!    PRINT*,database%bandInfo(i)%suffix
!    PRINT*,database%bandInfo(i)%number
!    PRINT*,database%bandInfo(i)%centerFreqIF
!    PRINT*,database%bandInfo(i)%spectrometerFamily
!    PRINT*,database%bandInfo(i)%spectrometerFamilyIndex
! ENDDO

! DO i=1,database%noSpectrometers
!    PRINT*,'--------------------'
!    PRINT*,database%spectrometerInfo(i)%name
!    PRINT*,database%spectrometerInfo(i)%fullFamily
!    PRINT*,database%spectrometerInfo(i)%number
!    PRINT*,database%spectrometerInfo(i)%familyIndex
! ENDDO   

loop: DO
   READ*,request
   IF (TRIM(request)=="DONE") EXIT loop
   CALL ParseMLSSignalRequest(request,database,signal)
   DO i=1,SIZE(signal)
      CALL GetFullMLSSignalName(signal(i),request)
      PRINT*,'You mean: '//TRIM(request)
      ! PRINT*,PACK(signal(i)%channelPosition,signal(i)%channelIncluded)
   ENDDO
   CALL DestroyMLSSignalsInfo(signal)
END DO loop
CALL DestroySignalsDatabase(database)
END PROGRAM Test



