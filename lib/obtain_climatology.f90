
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
PUBLIC

PRIVATE :: Id, ModuleName
!------------------------------- RCS Ident Info ------------------------------
CHARACTER(LEN=130) :: id = & 
   "$Id$"
CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!-----------------------------------------------------------------------------


IMPLICIT NONE


CONTAINS


  SUBROUTINE l3ascii_read_field(unit, qty)
    USE dates_module    ! Shoud use SDP Toolkit eventually. 
    ! ----Arguments ----!
    INTEGER,INTENT(in)::unit
    TYPE(GriddedData_T),INTENT(inout)::qty
    !-------Local Variables --------!
    CHARACTER(len=*),PARAMETER::dummyyear="1993"
    LOGICAL::opened
    CHARACTER(len=LineLen)::inline
    CHARACTER(len=30)::linetype,axistype,sdstring,edstring
    CHARACTER(len=80)::filename,unitstring
    REAL(kind=r8),POINTER,DIMENSION(:)::tmpaxis,dateStarts,dateEnds
    INTEGER::tmpaxis_len,idate
    INTEGER,PARAMETER::maxNoDates=30
    REAL(kind=r8),ALLOCATABLE,DIMENSION(:,:,:,:,:,:)::tmpqty
    LOGICAL::end_of_file
    !---- Executable statements ----! 
    NULLIFY(tmpaxis)

    WRITE(unit=unitstring,fmt="(i3)") unit ! For use in error reporting
    INQUIRE(unit=unit,opened=opened)
    IF (.NOT. opened) THEN
       CALL MLSMessage(MLSMSG_Error,ModuleName,&
       " in subroutine l3ascii_read_qty, Unit "//TRIM(unitstring)//&
       "is not connected to a file. Do call l3ascii_open(filename,unit) first")
       RETURN
    ENDIF
    INQUIRE(unit=unit,name=filename) ! find the file name connected to this
    ! unit for use in error messages.

    ! Fix axis arrays: set to default values with length 1 and a sensible 
    ! value. These will be used if the file does not have variation 
    ! along that axis

    ! Unfortunately, we can't tell if this defined type has
    ! been used before by testing a pointer for being associated as for 
    ! a new  struct they are in the undefined state. (This will change with 
    ! Fortran 95, where you can define a pointer to be initialised to the 
    ! nullified state)
    ! We therefore use a special number qty%reusing to test whether 
    ! this defined type has been used before.  
    IF (qty%reusing==313323435) THEN
       !       print*,"This struct has been used before : deallocating"
       DEALLOCATE(qty%heights)
       DEALLOCATE(qty%lats)
       DEALLOCATE(qty%lons)
       DEALLOCATE(qty%lsts)
       DEALLOCATE(qty%szas)
       DEALLOCATE(qty%dateStarts)
       DEALLOCATE(qty%dateEnds)
       DEALLOCATE(qty%field)
    ELSE
       !       print*,"This is a new struct"
       qty%reusing=313323435
    ENDIF
    ALLOCATE(qty%heights(1:1))
    qty%heights(1)=1000.0
    qty%noHeights=1
    qty%verticalCoordinate=1
    ALLOCATE(qty%lats(1:1))
    qty%lats(1)=0.0
    qty%noLats=1
    qty%equivalentLatitude=.FALSE.
    ALLOCATE(qty%lons(1:1))
    qty%lons(1)=0.0
    qty%noLons=1
    ALLOCATE(qty%lsts(1:1))
    qty%lsts(1)=12.0
    qty%noLsts=1
    ALLOCATE(qty%szas(1:1))
    qty%szas(1)=30.0
    qty%noSzas=1
    ! Dates are mandatory, so we don't have to give them a default value

    !--- Read qty info and all the axis info ---!
!    call get_next_noncomment_line(unit,inline)

    CALL ReadCompleteLineWithoutComments(unit,inline,eof=end_of_file)
    IF(Capitalize(inline(1:5)) /= "FIELD" ) THEN
      CALL MLSMessage(MLSMSG_Error,ModuleName,&
           "in subroutine l3ascii_read_field, File "//TRIM(filename)// &
           "on unit"//TRIM(unitstring)//" contains no more Fields")
       RETURN
    ENDIF
    
    IF(end_of_file) THEN
       CALL MLSMessage(MLSMSG_Error,ModuleName,&
       "In subroutine l3ascii_read_field, End of File"//TRIM(filename)// &
       " on unit"//TRIM(unitstring))
       RETURN
    ENDIF

    READ(unit=inline,fmt=*)linetype,qty%quantityName, &
         qty%description,qty%units
    axesloop:DO

        CALL ReadCompleteLineWithoutComments(unit,inline)
        READ(unit=inline,fmt=*)linetype,axistype
        linetype=Capitalize(linetype)
        axistype=Capitalize(axistype)
        IF(linetype(1:4) == "DATE") THEN ! This is always the last "axis"
           EXIT axesloop                 ! and is different from the others
        ENDIF

        CALL ParseVertCoordSpec(inline,tmpaxis, family)

        IF(linetype(1:8) == "PRESSURE" .OR. linetype(1:8) == "ALTITUDE" &
             .OR. linetype(1:3) == "GPH" .OR. linetype(1:5) == "THETA") THEN
           qty%noHeights=len(tmpaxis)
           DEALLOCATE(qty%heights)
           ALLOCATE(qty%heights(qty%noHeights))
           qty%heights=tmpaxis

           IF(linetype(1:8) == "PRESSURE") THEN
              qty%verticalCoordinate=1
           ELSE IF (linetype(1:8) == "ALTITUDE") THEN
              qty%verticalCoordinate=2
           ELSE IF (linetype(1:3) == "GPH") THEN
              qty%verticalCoordinate=3
           ELSE IF (linetype(1:5) == "THETA") THEN
              qty%verticalCoordinate=4
           ENDIF
        ELSE IF(linetype(1:8) == "LATITUDE" .OR. &
             linetype(1:8) == "EQUIVLAT") THEN
           qty%noLats=len(tmpaxis)
           DEALLOCATE(qty%lats)
           ALLOCATE(qty%lats(qty%noLats))
           qty%lats=tmpaxis
           IF (linetype(1:8) == "LATITUDE") THEN
              qty%equivalentLatitude=.FALSE.
           ELSE 
              qty%equivalentLatitude=.TRUE.
           ENDIF
        ELSE IF(linetype(1:9) == "LONGITUDE") THEN
           qty%noLons=len(tmpaxis)
           DEALLOCATE(qty%lons)
           ALLOCATE(qty%lons(qty%noLons))
           qty%lons=tmpaxis
        ELSE IF(linetype(1:9) == "LST") THEN
           qty%noLsts=len(tmpaxis)
           DEALLOCATE(qty%lsts)
           ALLOCATE(qty%lsts(qty%noLsts))
           qty%lsts=tmpaxis
        ELSE IF(linetype(1:9) == "SZA") THEN
           qty%noSzas=len(tmpaxis)
           DEALLOCATE(qty%szas)
           ALLOCATE(qty%szas(qty%noSzas))
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
         1:qty%noLsts,1:qty%noSzas,1:maxNoDates))
    ALLOCATE(dateStarts(1:maxNoDates),dateEnds(1:maxNoDates))

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
        READ(unit=inline,fmt=*)linetype,axistype        
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
                "in subroutine l3ascii_read_field: File"//TRIM(filename)//&
                "on unit"//TRIM(unitstring)//" contains a line beginning"//&
                TRIM(linetype)//"where I expected a line beginning Date ")
           RETURN
        ENDIF
        qty%noDates=qty%noDates+1

    ENDDO datesloop

    ALLOCATE(qty%field(1:qty%noHeights,1:qty%noLats,&
         1:qty%noLons,1:qty%noLsts,1:qty%noSzas,1:qty%noDates))
    ALLOCATE(qty%dateStarts(1:qty%noDates),&
         qty%dateEnds(1:qty%noDates))
    qty%dateStarts=dateStarts(1:qty%noDates)
    qty%dateEnds=dateEnds(1:qty%noDates)
    qty%field=tmpqty(:,:,:,:,:,1:qty%noDates)
    DEALLOCATE(tmpqty,dateStarts,dateEnds)

  END SUBROUTINE l3ascii_read_field

============================
END MODULE ObtainClimatology
============================
! $Log:

