
! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
!===================================================================
MODULE ObtainClimatology 

!provides subroutines to access climatology files in l3 ascii format
!===================================================================
USE MLSCommon
USE MLSStrings 
USE MLSMessageModule
USE GriddedData
USE VerticalCoordinate
USE MLSCF
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
    INTEGER::family, idate , status, tmpaxis_len, word_count
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
!        CALL ParseVertCoordSpec(axistype, tmpaxis, family)
        IF (axistype(1:6) =="LINEAR") THEN
           CALL make_linear_axis(inline,tmpaxis,tmpaxis_len)
        ELSE IF (axistype(1:3) =="LOG") THEN

           CALL make_log_axis(inline,tmpaxis,tmpaxis_len)

        ELSE IF (axistype(1:8) =="EXPLICIT") THEN

           BACKSPACE(unit=unit)
           CALL read_explicit_axis(unit,tmpaxis,tmpaxis_len)

        ELSE
           CALL MLSMessage(MLSMSG_Error,ModuleName,&
                "in subroutine l3ascii_read_field,File"//TRIM(filename)//&
                " on unit"//TRIM(unitstring)//" contains coordinate"//&
                " of invalid type "//TRIM(axistype)//"for axis"//&
                TRIM(linetype))
           RETURN
        ENDIF


        IF(linetype(1:8) == "PRESSURE" .OR. linetype(1:8) == "ALTITUDE" &
             .OR. linetype(1:3) == "GPH" .OR. linetype(1:5) == "THETA") THEN
           qty%noHeights=tmpaxis_len
           DEALLOCATE(qty%heights)
           ALLOCATE(qty%heights(1:tmpaxis_len))
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
           qty%noLats=tmpaxis_len
           DEALLOCATE(qty%lats)
           ALLOCATE(qty%lats(1:tmpaxis_len))
           qty%lats=tmpaxis
           IF (linetype(1:8) == "LATITUDE") THEN
              qty%equivalentLatitude=.FALSE.
           ELSE 
              qty%equivalentLatitude=.TRUE.
           ENDIF
        ELSE IF(linetype(1:9) == "LONGITUDE") THEN
           qty%noLons=tmpaxis_len
           DEALLOCATE(qty%lons)
           ALLOCATE(qty%lons(1:tmpaxis_len))
           qty%lons=tmpaxis
        ELSE IF(linetype(1:3) == "LST") THEN
           qty%noLsts=tmpaxis_len
           DEALLOCATE(qty%lsts)
           ALLOCATE(qty%lsts(1:tmpaxis_len))
           qty%lsts=tmpaxis
        ELSE IF(linetype(1:3) == "SZA") THEN
           qty%noSzas=tmpaxis_len
           DEALLOCATE(qty%szas)
           ALLOCATE(qty%szas(1:tmpaxis_len))
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
       word_count=count_words(inline)
       IF (word_count == 3) THEN
          READ(unit=inline,fmt=*)linetype,axistype,sdstring
          sdstring=ADJUSTL(sdstring)
          edstring=sdstring
       ELSE IF(word_count >= 4 ) THEN
          READ(unit=inline,fmt=*)linetype,axistype,sdstring,edstring
          sdstring=ADJUSTL(sdstring)
          edstring=ADJUSTL(edstring)
       ELSE
          CALL MLSMessage(MLSMSG_Error,ModuleName,&
               "in subroutine l3ascii_read_field: File "//TRIM(filename)//&
               "on unit "//TRIM(unitstring)//" contains a line beginning"//&
               TRIM(linetype)//" Date with too few words ")
       ENDIF
       
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
        
        READ(unit=unit,fmt=*)tmpqty(:,:,:,:,:,idate)
        end_of_file=.FALSE.
        CALL ReadCompleteLineWithoutComments(unit,inline,eof=end_of_file)
        IF(end_of_file) THEN 
           ! No more dates and nothing else either
           EXIT datesloop
        ENDIF

        READ(unit=inline,fmt=*)linetype,axistype        
        linetype=Capitalize(linetype)
        axistype=Capitalize(axistype)
        IF(end_of_file .OR. &
             linetype(1:5) == "FIELD" .OR. linetype(1:3)=="END") THEN 
           !Oops! There were no more dates, but there is another field
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
!==================================
  END SUBROUTINE l3ascii_read_field
!==================================
!===============================================
  SUBROUTINE make_log_axis(inline,axis,axis_len)
!===============================================
    !--------args------------!
    CHARACTER(len=*),INTENT(in)::inline
    REAL(kind=r8),POINTER,DIMENSION(:)::axis
    INTEGER,INTENT(out)::axis_len
    !-------locals--------------!
    CHARACTER(len=30)::linetype,axistype
    REAL(kind=r8)::basepressure
    INTEGER,DIMENSION(:),ALLOCATABLE::n_levs_in_sec,n_levs_per_dec,axints
    INTEGER::nwords,j,nsections,stind,st,i
    REAL(kind=r8)::gridstep
    !-------Executable----------!

    ! Warning: axis must be nullified or associated!
    IF (ASSOCIATED(axis)) THEN 
       DEALLOCATE(axis)
    ENDIF

    !Count words in inline. 
    nwords=1
    DO j=2,LEN(inline)
       IF(inline(j:j) /= " " .AND. inline(j-1:j-1) == " ") THEN
          nwords=nwords+1
       ENDIF
    ENDDO
    nsections=(nwords-3)/2
    ALLOCATE(n_levs_in_sec(1:nsections),n_levs_per_dec(1:nsections),&
         axints(1:nsections*2))
    READ(unit=inline,fmt=*)linetype,axistype,basepressure,axints
    n_levs_in_sec=axints(1:nsections*2-1:2)
    n_levs_per_dec=axints(2:nsections*2:2)

    axis_len=SUM(n_levs_in_sec)

    ALLOCATE(axis(1:axis_len))
    axis(1)=-LOG10(basepressure)
    stind=0
    DO j=1,nsections
       gridstep=1.0_r8/n_levs_per_dec(j)
       IF (j == 1) THEN
          st=2
       ELSE
          st=1
       ENDIF
       DO i=st,n_levs_in_sec(j)
          axis(stind+i)=axis(stind+i-1)+gridstep
       END DO
       stind=stind+n_levs_in_sec(j)
    END DO
    axis=10.0_r8**(-axis)
    DEALLOCATE(n_levs_in_sec, n_levs_per_dec,axints)
!=============================
  END SUBROUTINE make_log_axis
!=============================
!==================================================
  SUBROUTINE make_linear_axis(inline,axis,axis_len)
!==================================================
    !--------args------------!
    CHARACTER(len=*),INTENT(in)::inline
    REAL(kind=r8),POINTER,DIMENSION(:)::axis
    INTEGER,INTENT(out)::axis_len
    !-------locals--------------!
    CHARACTER(len=30)::linetype,axistype
    REAL(kind=r8)::baseval
    INTEGER,DIMENSION(:),ALLOCATABLE::n_levs_in_sec
    REAL(kind=r8),DIMENSION(:),ALLOCATABLE::gridstep,axints
    INTEGER::nwords,j,nsections,stind,st,i
    !-------Executable----------!
    ! Warning: axis must be nullified or associated!
    IF (ASSOCIATED(axis)) THEN 
       DEALLOCATE(axis)
    ENDIF

    !Count words in inline. 

    nwords=1
    DO j=2,LEN(inline)
       IF(inline(j:j) /= " " .AND. inline(j-1:j-1) == " ") THEN
          nwords=nwords+1
       ENDIF
    ENDDO
    nsections=(nwords-3)/2

    ALLOCATE(n_levs_in_sec(1:nsections),gridstep(1:nsections),&
         axints(1:nsections*2))

    READ(unit=inline,fmt=*)linetype,axistype,baseval,axints
    n_levs_in_sec=NINT(axints(1:nsections*2-1:2))
    gridstep=axints(2:nsections*2:2)

    axis_len=SUM(n_levs_in_sec)

    ALLOCATE(axis(1:axis_len))
    axis(1)=baseval
    stind=0
    DO j=1,nsections
       IF (j == 1) THEN
          st=2
       ELSE
          st=1
       ENDIF
       DO i=st,n_levs_in_sec(j)
          axis(stind+i)=axis(stind+i-1)+gridstep(j)
       END DO
       stind=stind+n_levs_in_sec(j)
    END DO
    DEALLOCATE(axints,gridstep,n_levs_in_sec)
!================================
  END SUBROUTINE make_linear_axis
!================================
!==================================================
  SUBROUTINE read_explicit_axis(unit,axis,axis_len)
!==================================================
    !--------args------------!
    INTEGER,INTENT(in)::unit
    REAL(kind=r8),POINTER,DIMENSION(:)::axis
    INTEGER,INTENT(out)::axis_len

    !------- Local vars ---------!
    CHARACTER(len=30)::readitem
    CHARACTER(len=1)::rdchar
    REAL(kind=r8),DIMENSION(1:200)::tmpaxis
    INTEGER::i,iotest
    LOGICAL::foundcb
    !Executables
    ! Warning: axis must be nullified or associated!
    IF (ASSOCIATED(axis)) THEN 
       DEALLOCATE(axis)
    ENDIF

    ! An explicit axis is supplied as a parenthesised list, spread 
    ! over several lines. This is a Royal PIA.
    ! Read chars till we get to the (
    DO 
       READ(unit=unit,fmt="(a)",advance="no")rdchar
  
       IF(rdchar=="(") THEN
          EXIT 
       ENDIF
    ENDDO
    ! now read items and add them to the axis until we get to the )
    axis_len=0
    foundcb=.FALSE.
    itemsloop:DO
       i=1
       charsloop:DO
          READ(unit=unit,fmt="(a)",advance="no",iostat=iotest)rdchar

          IF(rdchar==")") THEN
             foundcb=.TRUE.
             EXIT charsloop
          ENDIF
          IF (rdchar==" ") THEN
             EXIT charsloop
          ENDIF
          readitem(i:i)=rdchar
          i=i+1
       ENDDO charsloop
       IF ((i<=1 .AND. .NOT.foundcb) .OR. (i==2 .AND. readitem(1:1)==" ")) THEN
          CYCLE itemsloop
       ENDIF
       IF ( i <= 1 .AND. foundcb) THEN
          EXIT itemsloop
       ENDIF
       axis_len=axis_len+1
       READ(unit=readitem(1:i-1),fmt=*)tmpaxis(axis_len)
        IF (foundcb) THEN
           EXIT itemsloop
        ENDIF
    END DO itemsloop
    ALLOCATE(axis(1:axis_len))
    axis=tmpaxis(1:axis_len)
!==================================
  END SUBROUTINE read_explicit_axis
!==================================
!===============================================
  SUBROUTINE Obtain_Clim(aprioriData, l2cf_data)
!===============================================
!Arguments 
    TYPE (GriddedData_T), DIMENSION(:), POINTER :: aprioriData 
    ! Input a priori database
    TYPE (MLSCF_T) :: l2cf_data        ! The information from the l2cf file

!Local Variables

    TYPE (GriddedData_T):: qty
    CHARACTER (LEN=256) :: msg, mnemonic
    INTEGER:: CliUnit, processCli, returnStatus, version

    LOGICAL :: end_of_file = .FALSE.

    DO CliUnit = mlspcf_l2clim_start, mlspcf_l2clim_end

!   Open one Climatology file as a generic file for reading
      version = 1
      returnStatus = Pgs_io_gen_openF (CliUnit, PGSd_IO_Gen_RSeqFrm, 0, &
                                       processCli, version)
      IF (returnStatus /= PGS_S_SUCCESS) THEN

        CALL Pgs_smf_getMsg(returnStatus, mnemonic, msg)
        CALL MLSMessage (MLSMSG_Error, ModuleName, &
                         "Error opening Climatology file:  "//mnemonic//" "//msg)

      ENDIF


      DO WHILE (.NOT. end_of_file)

        CALL l3ascii_read_field(processCli, qty, end_of_file)
!        CALL AddGridTemplateToDatabase(aprioriData, qty)
        CALL DestroyGridTemplateContents(qty)

      END DO !(.not. end_of_file)

    END DO !CliUnit = mlspcf_l2clim_start, mlspcf_l2clim_end

  RETURN
!============================
  END SUBROUTINE Obtain_Clim
!============================
!============================
END MODULE ObtainClimatology
!============================
! $Log$










