
! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE ObtainClimatology !provides subroutines to access climatology files in l3 ascii format

!=============================================================================
USE MLSCommon
USE MLSStrings 
USE MLSMessageModule
USE GriddedData
USE VerticalCoordinate
USE MLSPCF
USE dates_module    

IMPLICIT NONE

PUBLIC

PRIVATE :: Id, ModuleName
!------------------------------- RCS Ident Info ------------------------------
CHARACTER(LEN=130) :: id = & 
   "$Id$"
CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!-----------------------------------------------------------------------------

CONTAINS

!=======================================================
  SUBROUTINE l3ascii_read_field(unit, qty, end_of_file)
!=======================================================
 
    ! ----Arguments ----!
    INTEGER,INTENT(in)::unit
    TYPE(GriddedData_T),INTENT(inout)::qty
    LOGICAL,INTENT(INOUT)::end_of_file

    !-------Local Variables --------!
    CHARACTER(len=*),PARAMETER::dummyyear="1993"
    LOGICAL::opened
    CHARACTER(len=LineLen)::inline
    CHARACTER(len=30)::linetype,axistype,sdstring,edstring
    CHARACTER(len=80)::filename,unitstring
    REAL(kind=r8),POINTER,DIMENSION(:)::tmpaxis,dateStarts,dateEnds
    INTEGER::family, idate , status
    INTEGER,PARAMETER::maxNoDates=30
    REAL(kind=r8),ALLOCATABLE,DIMENSION(:,:,:,:,:,:)::tmpqty

    !---- Executable statements ----! 
    NULLIFY(tmpaxis)

    WRITE(unit=unitstring,fmt="(i4)") unit ! For use in error reporting
    INQUIRE(unit=unit,opened=opened)
    IF (.NOT. opened) THEN
       CALL MLSMessage(MLSMSG_Error,ModuleName,&
       "Unit "//TRIM(unitstring)//" is not connected to a file.")
    ENDIF
    INQUIRE(unit=unit,name=filename) ! find the file name connected to this
    ! unit for use in error messages.

    ! Fix axis arrays: set to default values with length 1 and a sensible 
    ! value. These will be used if the file does not have variation 
    ! along that axis


 
    ALLOCATE(qty%heights(1:1),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
                 MLSMSG_Allocate//" heights")
   

    qty%heights(1)=1000.0
    qty%noHeights=1
    qty%verticalCoordinate=1
    ALLOCATE(qty%lats(1:1),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
                MLSMSG_Allocate//" lats")
   
    qty%lats(1)=0.0
    qty%noLats=1
    qty%equivalentLatitude=.FALSE.
    ALLOCATE(qty%lons(1:1),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
                 MLSMSG_Allocate//" lons")
   
    qty%lons(1)=0.0
    qty%noLons=1
    ALLOCATE(qty%lsts(1:1),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
                 MLSMSG_Allocate//" lsts")
    
    qty%lsts(1)=12.0
    qty%noLsts=1
    ALLOCATE(qty%szas(1:1),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
                 MLSMSG_Allocate//" szas")
    
    qty%szas(1)=30.0
    qty%noSzas=1
    ! Dates are mandatory, so we don't have to give them a default value

    !--- Read qty info and all the axis info ---!

    CALL ReadCompleteLineWithoutComments(unit,inline,eof=end_of_file)
    IF(Capitalize(inline(1:5)) /= "FIELD" ) THEN
      CALL MLSMessage(MLSMSG_Error,ModuleName,&
           "File "//TRIM(filename)// &
           "on unit"//TRIM(unitstring)//" contains no more Fields")
     
    END IF
    
    IF(end_of_file) THEN
       CALL MLSMessage(MLSMSG_Error,ModuleName,&
       "End of File"//TRIM(filename)// &
       " on unit"//TRIM(unitstring))
    END IF

    READ(unit=inline,fmt=*)linetype,qty%quantityName, &
         qty%description, qty%units
    axesloop:DO

        CALL ReadCompleteLineWithoutComments(unit,inline)
        READ(unit=inline,fmt=*)linetype,axistype
        linetype=Capitalize(linetype)
        axistype=Capitalize(axistype)
        IF(linetype(1:4) == "DATE") THEN ! This is always the last "axis"
           EXIT axesloop                 ! and is different from the others
        ENDIF

        CALL ParseVertCoordSpec(axistype, tmpaxis, family)

        IF((linetype(1:8) == "PRESSURE") .OR. (linetype(1:8) == "ALTITUDE") &
           .OR. (linetype(1:3) == "GPH") .OR. (linetype(1:5) == "THETA")) THEN
           qty%noHeights=SIZE(tmpaxis)
           DEALLOCATE(qty%heights)
           ALLOCATE(qty%heights(qty%noHeights),STAT=status)
           IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
                        MLSMSG_Allocate//" heights")
           
           qty%heights=tmpaxis
           qty%verticalCoordinate = ParseVerticalCoordinateName(linetype)
           IF(qty%verticalCoordinate ==  VC_Invalid)THEN
             CALL MLSMessage(MLSMSG_Error,ModuleName,&
                             "Invalid vertical coordinate :"//linetype)
           END IF
        ELSE IF((linetype(1:8) == "LATITUDE") .OR. &
                (linetype(1:8) == "EQUIVLAT")) THEN
           qty%noLats=SIZE(tmpaxis)
           DEALLOCATE(qty%lats)
           ALLOCATE(qty%lats(qty%noLats),STAT=status)
           IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
                        MLSMSG_Allocate//" lats")

           qty%lats=tmpaxis
           IF (linetype(1:8) == "LATITUDE") THEN
              qty%equivalentLatitude=.FALSE.
           ELSE 
              qty%equivalentLatitude=.TRUE.
           ENDIF
        ELSE IF(linetype(1:9) == "LONGITUDE") THEN
           qty%noLons=SIZE(tmpaxis)
           DEALLOCATE(qty%lons)
           ALLOCATE(qty%lons(qty%noLons),STAT=status)
           IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
                        MLSMSG_Allocate//" lons")
   

           qty%lons=tmpaxis
        ELSE IF(linetype(1:9) == "LST") THEN
           qty%noLsts=SIZE(tmpaxis)
           DEALLOCATE(qty%lsts)
           ALLOCATE(qty%lsts(qty%noLsts),STAT=status)
           IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
                        MLSMSG_Allocate//" lsts")
   
           qty%lsts=tmpaxis
        ELSE IF(linetype(1:9) == "SZA") THEN
           qty%noSzas=SIZE(tmpaxis)
           DEALLOCATE(qty%szas)
           ALLOCATE(qty%szas(qty%noSzas),STAT=status)
           IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
                        MLSMSG_Allocate//" szas")
   

           qty%szas=tmpaxis
        ENDIF
    ENDDO axesloop
    DEALLOCATE(tmpaxis)

    ! We already have the first date line read and the axis type extracted
    ! It was the existence of a date line that caused us to exit from 
    ! the loop axesloop. We don't know how many dates there are so we have to 
    ! allocate large arrays and copy their contents to an array of the 
    ! right size 
    qty%noDates=1
    ALLOCATE(tmpqty(1:qty%noHeights,1:qty%noLats,1:qty%noLons, &
         1:qty%noLsts,1:qty%noSzas,1:maxNoDates),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
                 MLSMSG_Allocate//" tmpqty")
    ALLOCATE(dateStarts(1:maxNoDates),dateEnds(1:maxNoDates),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
                 MLSMSG_Allocate//" dateStarts")
    ! Loop to read in the data for the current date and check to see if 
    ! there is another date
    datesloop: DO idate=1,maxNoDates
        !print*,"Datesloop: idate=",idate
        READ(unit=inline,fmt=*)linetype,axistype,sdstring,edstring
        sdstring=ADJUSTL(sdstring)
        edstring=ADJUSTL(edstring)
        ! Date strings can begin with - indicating the year is 
        ! missing and that the file belongs to no year in particular.
        ! To convert dates to SDP toolkit (TAI) times (Seconds since start of 
        ! 1 Jan 1993) we need to stick on a dummy year
        IF (sdstring(1:1) == "-") THEN
           sdstring=dummyyear//sdstring
        ENDIF
        IF (edstring(1:1) == "-") THEN
           edstring=dummyyear//edstring
        ENDIF
        ! ccsds2tai returns days since 1 Jan 1993. 86400==no of secs per day
        dateStarts(idate)=86400*ccsds2tai(sdstring)
        dateEnds(idate)  =86400*ccsds2tai(edstring)
        CALL ReadCompleteLineWithoutComments(unit,inline)
        BACKSPACE(unit=unit)
        !print*,"About to read data: inline=",inline
        !print*,"tmpfield has size:",size(tmpfield),shape(tmpqty)
        
        READ(unit=unit,fmt=*)tmpqty(:,:,:,:,:,idate)
        !print*,"Read data"
        end_of_file=.FALSE.
        CALL ReadCompleteLineWithoutComments(unit,inline,eof=end_of_file)
        !print*,"Next date line:",inline,"EOF=",end_of_file
        IF(end_of_file) THEN 
           ! No more dates and nothing else either
           EXIT datesloop
        ENDIF
        READ(unit=inline,fmt=*)linetype, axistype        
        linetype=Capitalize(linetype)
        axistype=Capitalize(axistype)
        IF(end_of_file .OR. &
             linetype(1:5) == "FIELD" .OR. linetype(1:3)=="END") THEN 
           !Oops! There were no more dates, but there is another qty
           BACKSPACE(unit=unit)
           EXIT datesloop                 
        ENDIF
        IF(linetype(1:4) /= "DATE") THEN ! There should be another date here
           CALL MLSMessage(MLSMSG_Error,ModuleName,&
                "File"//TRIM(filename)//&
                "on unit"//TRIM(unitstring)//" contains a line beginning"//&
                TRIM(linetype)//"where I expected a line beginning Date ")
           RETURN
        ENDIF
        qty%noDates=qty%noDates+1

    ENDDO datesloop

    ALLOCATE(qty%field(1:qty%noHeights,1:qty%noLats,&
             1:qty%noLons,1:qty%noLsts,1:qty%noSzas,1:qty%noDates),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
                 MLSMSG_Allocate//" field")
 
   ALLOCATE(qty%dateStarts(1:qty%noDates),&
           qty%dateEnds(1:qty%noDates),STAT=status)

    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
                 MLSMSG_Allocate//" dateStarts")

    qty%dateStarts=dateStarts(1:qty%noDates)
    qty%dateEnds=dateEnds(1:qty%noDates)
    qty%field=tmpqty(:,:,:,:,:,1:qty%noDates)
    DEALLOCATE(tmpqty,dateStarts,dateEnds)

  END SUBROUTINE l3ascii_read_field

  SUBROUTINE Obtain_Clim(aprioriData)
 
    TYPE (GriddedData_T), DIMENSION(:), POINTER :: aprioriData 
    ! Input a priori database

    TYPE (GriddedData_T):: qty
    CHARACTER (LEN=256) :: msg, mnemonic
    INTEGER:: CliUnit, processCli, returnStatus, version

    LOGICAL :: end_of_file = .FALSE.

    DO CliUnit = mlspcf_l2_cli_start, mlspcf_l2_cli_end

!   Open one Climatology file as a generic file for reading

      returnStatus = Pgs_io_gen_openF (CliUnit, PGSd_IO_Gen_RSeqFrm, 0, &
                                       processCli, version)

      IF (returnStatus /= PGS_S_SUCCESS) THEN

        CALL Pgs_smf_getMsg(returnStatus, mnemonic, msg)
        CALL MLSMessage (MLSMSG_Error, ModuleName, &
                         "Error opening Climatology file:  "//mnemonic//" "//msg)

      ENDIF


      Do while (.not. end_of_file)

        call l3ascii_read_field(processCli, qty, end_of_file)
        CALL AddGridTemplateToDatabase(aprioriData, qty)
        call DestroyGridTemplateContents(qty)

      END DO !(.not. end_of_file)

    END DO !CliUnit = mlspcf_l2_cli_start, mlspcf_l2_cli_end

  RETURN
  END SUBROUTINE Obtain_Clim
============================
END MODULE ObtainClimatology
============================
! $Log:
