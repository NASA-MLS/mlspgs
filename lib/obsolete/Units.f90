! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE Units                    ! Units for physical quantities
!=============================================================================

  USE MLSCommon
  USE MLSStrings
  USE MLSMessageModule

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------


  ! This module defines standard units for physical quantities.

  ! First we define some standard physical quantity types.  Feel free to add
  ! more to these

  INTEGER, PARAMETER :: PHYQ_Invalid=0 ! Invalid unit given by user
  INTEGER, PARAMETER :: PHYQ_Dimensionless=1 ! Dimensionless quantity
  INTEGER, PARAMETER :: PHYQ_Length=2 ! Default meters
  INTEGER, PARAMETER :: PHYQ_Time=3 ! Default seconds
  INTEGER, PARAMETER :: PHYQ_Pressure=4 ! Default mb
  INTEGER, PARAMETER :: PHYQ_Temperature=5 ! Default K
  INTEGER, PARAMETER :: PHYQ_VMR=6 ! Default parts per `one'
  INTEGER, PARAMETER :: PHYQ_Angle=7 ! Default degrees
  INTEGER, PARAMETER :: PHYQ_MAFs=8 ! Default MAFs
  INTEGER, PARAMETER :: PHYQ_MIFs=9 ! Default MIFs
  INTEGER, PARAMETER :: PHYQ_Frequency=10 ! Defaut MHz
  INTEGER, PARAMETER :: PHYQ_Zeta=11 ! -ALOG10(pressure/hPa)

  INTEGER, PRIVATE, PARAMETER :: UnitNameLength=20

  TYPE PhysicalUnit_T
     CHARACTER (LEN=UnitNameLength) :: name ! e.g. km
     INTEGER :: family          ! e.g. PHYQ_Length
     REAL (r8) :: scale         ! e.g. 1e3
  END TYPE PhysicalUnit_T

  ! This parameter is a list of all the possible units, their families and
  ! associated scale factors.  Feel free to add more of these, just remember to:
  !   1 : Enter the names in UPPER CASE only
  !   2 : Increase the UnitsDatabaseSize accordingly.

  INTEGER, PARAMETER, PRIVATE :: UnitsDatabaseSize=36

  TYPE (PhysicalUnit_T), DIMENSION(UnitsDatabaseSize), &
       &    PARAMETER, PRIVATE :: physicalUnitDatabase= (/&

       & PhysicalUnit_T("DIMENSIONLESS", PHYQ_Dimensionless, 1.0), &
       & PhysicalUnit_T("DIMLESS", PHYQ_Dimensionless, 1.0), &
       & PhysicalUnit_T("DL", PHYQ_Dimensionless, 1.0), &

       & PhysicalUnit_T("M", PHYQ_Length, 1.0), &
       & PhysicalUnit_T("METERS", PHYQ_Length, 1.0), &
       & PhysicalUnit_T("KM", PHYQ_Length, 1e3), &

       & PhysicalUnit_T("S", PHYQ_Time, 1.0), &
       & PhysicalUnit_T("SECONDS", PHYQ_Time, 1.0), &
       & PhysicalUnit_T("MINUTES", PHYQ_Time, 60.0), &
       & PhysicalUnit_T("HOURS", PHYQ_Time, 3600.0), &
       & PhysicalUnit_T("DAYS", PHYQ_Time, 3600.0*24.0), &

       & PhysicalUnit_T("MB", PHYQ_Pressure,1.0), &
       & PhysicalUnit_T("HPA", PHYQ_Pressure,1.0), &
       & PhysicalUnit_T("PA", PHYQ_Pressure,1e-2), &

       & PhysicalUnit_T("K", PHYQ_Temperature, 1.0),&

       & PhysicalUnit_T("VMR", PHYQ_vmr, 1.0), &
       & PhysicalUnit_T("PPMV", PHYQ_vmr, 1e-6), &
       & PhysicalUnit_T("PPBV", PHYQ_vmr, 1e-9), &
       & PhysicalUnit_T("PPTV", PHYQ_vmr, 1e-12), &

       & PhysicalUnit_T("DEGREES", PHYQ_Angle, 1.0), &
       & PhysicalUnit_T("DEG", PHYQ_Angle, 1.0), &
       & PhysicalUnit_T("RADIANS", PHYQ_Angle, 57.295779513082323), &
       & PhysicalUnit_T("RAD", PHYQ_Angle, 57.295779513082323), &
       & PhysicalUnit_T("ORBITS", PHYQ_Angle, 360.0), &

       & PhysicalUnit_T("MAFS", PHYQ_MAFs, 1.0), &
       & PhysicalUnit_T("MAF", PHYQ_MAFs, 1.0), &

       & PhysicalUnit_T("MIFS", PHYQ_MIFs, 1.0), &
       & PhysicalUnit_T("MIF", PHYQ_MIFs, 1.0), &

       & PhysicalUnit_T("MHZ", PHYQ_Frequency, 1.0), &
       & PhysicalUnit_T("GHZ", PHYQ_Frequency, 1e3), &
       & PhysicalUnit_T("THZ", PHYQ_Frequency, 1e6), &
       & PhysicalUnit_T("KHZ", PHYQ_Frequency, 1e-3), &
       & PhysicalUnit_T("HZ", PHYQ_Frequency, 1e-6), &

       & PhysicalUnit_T("ZETA", PHYQ_Zeta, 1e-3), &
       & PhysicalUnit_T("LOG10(P)", PHYQ_Zeta, 1e-6), &
       & PhysicalUnit_T("LOGP", PHYQ_Zeta, 1e-6) &
       & /)

CONTAINS

  ! --------------------------------------------------------------------------

  ! This subroutine takes a unit name and returns it's family and scale factor

  SUBROUTINE ParseUnitName(name,family,scale)

    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(OUT) :: family
    REAL(r8), INTENT(OUT) :: scale

    ! Local variables
    CHARACTER (LEN=UnitNameLength) :: tmpName
    TYPE (PhysicalUnit_T), DIMENSION(1) :: thisUnit
    LOGICAL, DIMENSION(UnitsDatabaseSize) :: matches

    ! Executable code

    tmpName=Capitalize(name)
    matches=(tmpName==physicalUnitDatabase%name)

    ! If we don't get a match don't generate an error, just return an invalid
    ! flag.  We'll let the calling code deal with this, as it's going to be
    ! doing it's own tests anyway (e.g. this must be a length).

    IF (COUNT(matches) /= 1) THEN
       family=PHYQ_Invalid
       scale=1.0
    ELSE
       thisUnit=PACK(physicalUnitDatabase,matches)
       family=thisUnit(1)%family
       scale=thisUnit(1)%scale
    ENDIF
  END SUBROUTINE ParseUnitName

  ! --------------------------------------------------------------------------
  
  ! This subroutine takes a string and reads a number with a trailing unit
  ! from it.  It also returns the `remainder' string.  It also copes with cases
  ! where there is no unit.  In these cases, the value is returned, but
  ! the family is set to PHYQ_Invalid.

  SUBROUTINE ReadNumberWithUnitsFromString(line,value,family,remainder,&
       & unscaledValue,scale,unitIsOptional)

    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(IN) :: line ! Input line
    REAL (r8), INTENT(OUT) :: value ! Value of quantity
    INTEGER, INTENT(OUT) :: family ! Family for this quantity eg. PHYQ_Length
    CHARACTER (LEN=*), OPTIONAL, INTENT(OUT) :: remainder ! Rest of string
    REAL (r8), INTENT(OUT), OPTIONAL :: unscaledValue ! Value in file
    REAL (r8), INTENT(OUT), OPTIONAL :: scale ! Scale factor applied
    LOGICAL, INTENT(IN), OPTIONAL :: unitIsOptional ! User may omit unit

    ! Local variables
    CHARACTER (LEN=LEN(line)) :: firstWord
    CHARACTER (LEN=LEN(line)) :: lineRemainder
    CHARACTER (LEN=LEN(line)) :: finalRemainder
    CHARACTER (LEN=LEN(line)) :: numericValue
    CHARACTER (LEN=LEN(line)) :: unitName
    CHARACTER :: thisChar
    INTEGER :: wordLen
    INTEGER :: charNo
    LOGICAL :: unitIsReallyOptional ! unitIsOptional or FALSE if not present

    REAL(r8) :: useScale

    IF (PRESENT(unitIsOptional)) THEN
       unitIsReallyOptional=unitIsOptional
    ELSE
       unitIsReallyOptional=.FALSE.
    END IF

    ! Executable code

    ! First we'll split the first word off this string
    CALL SplitWords(line,firstWord,lineRemainder,delimiter=" ")

    ! Now we look from the end of the first word back and find the first (i.e.
    ! last!) 0-9 or . character

    wordLen=LEN_TRIM(firstWord)
    charNo=wordLen
    unitsReadLoop: DO
       thisChar=firstWord(charNo:charNo)
       IF ( (LGE(thisChar,"0").AND.LLE(thisChar,"9")) .OR. (thisChar==".")) &
            & EXIT unitsReadLoop
       charNo=charNo-1
       IF (charNo==0) EXIT unitsReadLoop
    END DO unitsReadLoop

    IF (charNo==0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & "<Number> <Unit> expected: "//firstWord)

    IF (charNo==wordLen) THEN
       ! The first word is numbers only, the unit may be the next word
       numericValue=firstWord
       CALL SplitWords(lineRemainder,unitName,finalRemainder,delimiter=" ")

       ! Now, check that the unit really is present
       IF ((LEN_TRIM(unitName)==0).OR.&
            & (INDEX("+-0123456789.",unitName(1:1))/=0)) THEN
          IF (unitIsReallyOptional) THEN
             unitName=""           ! No valid unit
             finalRemainder=lineRemainder ! Put this word back in the remainder
          ELSE
             CALL MLSMessage(MLSMSG_Error,ModuleName,"You must give a unit")
          END IF
       END IF
    ELSE
       ! The unit and the word are together split them ourselves
       finalRemainder=lineRemainder
       unitName=firstWord(charNo+1:)
       numericValue=firstWord(1:charNo)
    END IF

    ! Now read the numeric value
    READ (UNIT=numericValue,FMT=*) value
    
    ! Parse the unit name
    CALL ParseUnitName(unitName,family,useScale)

    ! Pass back the unscaled value if requested.
    IF (PRESENT(unscaledValue)) unscaledValue=value

    value=value*useScale

    ! Set the remainder if present
    IF (PRESENT(remainder)) remainder=finalRemainder

    ! Return the scale if requested to
    IF (PRESENT(scale)) scale=useScale

  END SUBROUTINE ReadNumberWithUnitsFromString

!=============================================================================
END MODULE Units
!=============================================================================

!
! $Log$
! Revision 2.0  2000/09/05 17:41:07  dcuddy
! Change revision to 2.0
!
! Revision 1.6  1999/12/16 23:24:43  livesey
! Added zeta family (can use zeta or log10(p) or logp)
!
