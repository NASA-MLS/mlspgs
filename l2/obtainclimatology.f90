! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
!===================================================================
module ObtainClimatology 

! Provides subroutines to access climatology files in L3 ascii format
!===================================================================

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use Dates_module, only: CCSDS2TAI
  use GriddedData, only: DestroyGridTemplateContents, GriddedData_T
  use MLSCommon, only: LineLen, R8
  use MLSMessagemodule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error
  use MLSPCF, only: mlspcf_l2clim_end, mlspcf_l2clim_start
  use MLSStrings, only: Capitalize, Count_Words, ReadCompleteLineWithoutComments
  use SDPToolkit, only: Pgs_io_gen_openF, PGS_S_SUCCESS, PGSd_IO_Gen_RSeqFrm
! use VerticalCoordinate

  implicit none
  private
  public :: OBTAIN_CLIM

  private :: Id, moduleName
  !------------------------------- RCS Ident Info ------------------------------
  character(len=130) :: id = & 
     "$Id$"
  character(len=*), parameter :: moduleName="$RCSfile$"
  !-----------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  ! ------------------------------------------------  OBTAIN_CLIM  -----
  !=====================================================================
  subroutine OBTAIN_CLIM ( aprioriData, root )
  !=====================================================================
    !Arguments 
    type (GriddedData_T), dimension(:), pointer :: aprioriData 
    ! Input a priori database
    integer, intent(in) :: root        ! Root of the L2CF abstract syntax tree

    !Local Variables

    type (GriddedData_T):: qty
    character (LEN=256) :: msg, mnemonic
    integer:: CliUnit, processCli, returnStatus, version

    logical :: end_of_file = .FALSE.

    do CliUnit = mlspcf_l2clim_start, mlspcf_l2clim_end

!     Open one Climatology file as a generic file for reading
      version = 1
      returnStatus = Pgs_io_gen_openF ( CliUnit, PGSd_IO_Gen_RSeqFrm, 0, &
                                        processCli, version )
      if ( returnStatus /= PGS_S_SUCCESS ) then

        call Pgs_smf_getMsg ( returnStatus, mnemonic, msg )
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & "Error opening Climatology file:  "//mnemonic//" "//msg)

      end if


      do while (.NOT. end_of_file)

        call l3ascii_read_field ( processCli, qty, end_of_file)
!        call AddGridTemplateToDatabase(aprioriData, qty)
        call DestroyGridTemplateContents ( qty )

      end do !(.not. end_of_file)

    end do ! CliUnit = mlspcf_l2clim_start, mlspcf_l2clim_end

  return
  !============================
  end subroutine OBTAIN_CLIM
  !============================

! =====     Private Procedures     =====================================

  ! -----------------------------------------  L3ASCII_READ_FIELD  -----
  !=====================================================================
  subroutine L3ASCII_READ_FIELD ( unit, qty, end_of_file )
  !=====================================================================
 
    ! ----Arguments ----!
    integer,intent(in) :: unit
    type(GriddedData_T), intent(inout) :: qty
    logical, intent(inout) :: end_of_file

    !-------Local Variables --------!
    character(len=*), parameter :: dummyyear = "1993"
    logical :: opened
    character(len=LineLen) :: inline
    character(len=30) :: linetype, axistype, sdstring, edstring
    character(len=80) :: filename, unitstring
    real(kind=r8), pointer, dimension(:) :: tmpaxis, dateStarts, dateEnds
!   integer :: family
    integer :: idate , status, tmpaxis_len, word_count
    integer, parameter :: maxNoDates = 30
    real(kind=r8), allocatable, dimension(:,:,:,:,:,:) :: tmpqty

    !---- Executable statements ----! 
    nullify ( tmpaxis )

    write ( unit=unitstring, fmt="(i4)" ) unit ! for use in error reporting
    inquire ( unit=unit, opened=opened )
    if ( .NOT. opened ) then
       call MLSMessage ( MLSMSG_Error, moduleName, &
         & "Unit "//TRIM(unitstring)//" is not connected to a file.")
    end if
    inquire ( unit=unit, name=filename ) ! find the file name connected to this
    ! unit for use in error messages.

    ! Fix axis arrays: set to default values with length 1 and a sensible 
    ! value. These will be used if the file does not have variation 
    ! along that axis
 
    call allocate_test ( qty%heights, 1, "qty%heights", ModuleName )
    qty%heights(1) = 1000.0
    qty%noHeights = 1
    qty%verticalCoordinate = 1

    call allocate_test ( qty%lats, 1, "qty%lats", ModuleName )
    qty%lats(1) = 0.0
    qty%noLats = 1
    qty%equivalentLatitude = .FALSE.

    call allocate_test ( qty%lons, 1, "qty%lons", ModuleName )
    qty%lons(1) = 0.0
    qty%noLons = 1

    call allocate_test ( qty%lsts, 1, "qty%lsts", ModuleName )
    qty%lsts(1) = 12.0
    qty%noLsts = 1

    call allocate_test ( qty%szas, 1, "qty%szas", ModuleName )
    qty%szas(1) = 30.0
    qty%noSzas = 1
    ! Dates are mandatory, so we don't have to give them a default value

    !--- read qty info and all the axis info ---!

    call readCompleteLineWithoutComments ( unit, inline, eof=end_of_file )

    if ( Capitalize(inline(1:5)) /= "FIELD" ) then
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & "File "//TRIM(filename)// &
        & "on unit"//TRIM(unitstring)//" contains no more Fields" )
     
    end if
    
    if ( end_of_file ) then
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & "end of File"//TRIM(filename)// &
        & " on unit"//TRIM(unitstring) )
    end if

    read ( unit=inline, fmt=* ) linetype, qty%quantityName, &
      & qty%description, qty%units

    axesloop:do

      call readCompleteLineWithoutComments ( unit, inline )
      read ( unit=inline,fmt=* ) linetype, axistype
      linetype = Capitalize(linetype)
      axistype = Capitalize(axistype)
      if ( linetype(1:4) == "DATE" ) then ! This is always the last "axis"
        exit axesloop                     ! and is different from the others
      end if
!      call ParseVertCoordSpec(axistype, tmpaxis, family)
      if (axistype(1:6) =="LINEAR") then
        call make_linear_axis(inline,tmpaxis,tmpaxis_len)
      else if (axistype(1:3) =="LOG") then

        call make_log_axis(inline,tmpaxis,tmpaxis_len)

      else if (axistype(1:8) =="EXPLICIT") then

        backspace ( unit=unit )
        call read_explicit_axis ( unit, tmpaxis, tmpaxis_len )

      else
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & "in subroutine l3ascii_read_field,File"//TRIM(filename)//&
          & " on unit"//TRIM(unitstring)//" contains coordinate"//&
          & " of invalid type "//TRIM(axistype)//"for axis"//&
          & TRIM(linetype))
         return
      end if


      if ( linetype(1:8) == "PRESSURE" .OR. linetype(1:8) == "ALTITUDE" &
        & .OR. linetype(1:3) == "GPH" .OR. linetype(1:5) == "THETA" ) then
        qty%noHeights = tmpaxis_len
        call deallocate_test ( qty%heights, "qty%heights", ModuleName )
        call allocate_test ( qty%heights, tmpaxis_len, "qty%heights", &
          & ModuleName )
        qty%heights = tmpaxis

        if ( linetype(1:8) == "PRESSURE" ) then
          qty%verticalCoordinate = 1
        else if ( linetype(1:8) == "ALTITUDE" ) then
          qty%verticalCoordinate = 2
        else if ( linetype(1:3) == "GPH" ) then
          qty%verticalCoordinate = 3
        else if ( linetype(1:5) == "THETA" ) then
          qty%verticalCoordinate = 4
        end if
      else if ( linetype(1:8) == "LATITUDE" .OR. &
        & linetype(1:8) == "EQUIVLAT" ) then
        qty%noLats = tmpaxis_len
        call deallocate_test ( qty%lats, "qty%lats", ModuleName )
        call allocate_test ( qty%lats, tmpaxis_len, "qty%lats", ModuleName )
        qty%lats=tmpaxis
        if (linetype(1:8) == "LATITUDE") then
          qty%equivalentLatitude = .FALSE.
        else
          qty%equivalentLatitude = .TRUE.
        end if
      else if ( linetype(1:9) == "LONGITUDE" ) then
        qty%noLons = tmpaxis_len
        call deallocate_test ( qty%lons, "qty%lons", ModuleName )
        call allocate_test ( qty%lons, tmpaxis_len, "qty%lons", ModuleName )
        qty%lons = tmpaxis
      else if ( linetype(1:3) == "LST" ) then
        qty%noLsts = tmpaxis_len
        call deallocate_test ( qty%lsts, "qty%lsts", ModuleName )
        call allocate_test ( qty%lsts, tmpaxis_len, "qty%lsts", ModuleName )
        qty%lsts = tmpaxis
      else if ( linetype(1:3) == "SZA" ) then
        qty%noSzas = tmpaxis_len
        call deallocate_test ( qty%szas, "qty%szas", ModuleName )
        call allocate_test ( qty%szas, tmpaxis_len, "qty%szas", ModuleName )
        qty%szas = tmpaxis
      end if
    end do axesloop

    call deallocate_test ( tmpaxis, "tmpaxis", ModuleName )

    ! We already have the first date line read and the axis type extracted
    ! It was the existence of a date line that caused us to exit from 
    ! the loop axesloop. We don't know how many dates there are so we have to 
    ! allocate large arrays and copy their contents to an array of the 
    ! right size 
    qty%noDates = 1
    allocate ( tmpqty(1:qty%noHeights, 1:qty%noLats, 1:qty%noLons, &
      & 1:qty%noLsts, 1:qty%noSzas, 1:maxNoDates), STAT=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_allocate//" tmpqty" )
    call allocate_test ( dateStarts, maxNoDates, "dateStarts", ModuleName )
    call allocate_test ( dateEnds, maxNoDates, "dateEnds", ModuleName )

    ! Loop to read in the data for the current date and check to see if 
    ! there is another date    

    datesloop: do idate = 1, maxNoDates
      word_count = count_words(inline)
      if ( word_count == 3 ) then
        read ( unit=inline, fmt=* ) linetype, axistype, sdstring
        sdstring = ADJUSTL(sdstring)
        edstring = sdstring
      else if ( word_count >= 4 ) then
        read ( unit=inline, fmt=* ) linetype, axistype, sdstring, edstring
        sdstring = ADJUSTL(sdstring)
        edstring = ADJUSTL(edstring)
      else
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & "in subroutine l3ascii_read_field: File "//TRIM(filename)// &
          & "on unit "//TRIM(unitstring)//" contains a line beginning"// &
          & TRIM(linetype)//" Date with too few words " )
      end if

       ! Date strings can begin with - indicating the year is 
       ! missing and that the file belongs to no year in particular.
       ! To convert dates to SDP toolkit (TAI) times (Seconds since start of 
       ! 1 Jan 1993) we need to stick on a dummy year
       if ( sdstring(1:1) == "-" ) then
         sdstring = dummyyear//sdstring
       end if
       if ( edstring(1:1) == "-" ) then
          edstring = dummyyear//edstring
       end if
       ! ccsds2tai returns days since 1 Jan 1993. 86400==no of secs per day
       dateStarts(idate) = 86400*ccsds2tai(sdstring)
       dateEnds(idate)   = 86400*ccsds2tai(edstring)
       call readCompleteLineWithoutComments ( unit, inline )
       backspace ( unit=unit )

       read ( unit=unit, fmt=* ) tmpqty(:,:,:,:,:,idate)
       end_of_file = .FALSE.
       call readCompleteLineWithoutComments ( unit, inline, eof=end_of_file )
       if ( end_of_file ) then 
         ! No more dates and nothing else either
         exit datesloop
       end if

       read ( unit=inline, fmt=* ) linetype, axistype        
       linetype = Capitalize(linetype)
       axistype = Capitalize(axistype)
       if ( end_of_file .OR. &
         &linetype(1:5) == "FIELD" .OR. linetype(1:3)=="END" ) then 
         !Oops! There were no more dates, but there is another field
         backspace(unit=unit)
         exit datesloop
       end if
       if ( linetype(1:4) /= "DATE" ) then ! There should be another date here
         call MLSMessage ( MLSMSG_Error, moduleName, &
           & "in subroutine l3ascii_read_field: File"//TRIM(filename)// &
           & "on unit"//TRIM(unitstring)//" contains a line beginning"// &
           & TRIM(linetype)//"where I expected a line beginning Date ")
         return
       end if
       qty%noDates = qty%noDates+1

    end do datesloop


    allocate ( qty%field(1:qty%noHeights, 1:qty%noLats, 1:qty%noLons, &
      & 1:qty%noLsts, 1:qty%noSzas, 1:qty%noDates), STAT=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_allocate//" field" )

    call allocate_test ( qty%dateStarts, qty%noDates, "qty%dateStarts", &
      & ModuleName )
    call allocate_test ( qty%dateEnds, qty%noDates, "qty%dateEnds", ModuleName )

    qty%dateStarts = dateStarts(1:qty%noDates)
    qty%dateEnds   = dateEnds(1:qty%noDates)
    qty%field = tmpqty(:,:,:,:,:,1:qty%noDates)
    deallocate ( tmpqty, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_deallocate//" tmpqty" )
    call deallocate_test ( dateStarts, "dateStarts", ModuleName )
    call deallocate_test ( dateEnds, "dateEnds", ModuleName )
  !==================================
  end subroutine l3ascii_read_field
  !==================================

  ! -------------------------------------------  MAKE_LINEAR_AXIS  -----
  !=====================================================================
  subroutine MAKE_LINEAR_AXIS ( inline, axis, axis_len )
  !=====================================================================
    !--------args------------!
    character(len=*), intent(in) :: inline
    real(kind=r8), pointer, dimension(:) :: axis
    integer, intent(out) :: axis_len
    !-------locals--------------!
    character(len=30) :: linetype, axistype
    real(kind=r8) :: baseval
    integer, dimension(:), pointer :: n_levs_in_sec
    real(kind=r8), dimension(:), pointer :: gridstep, axints
    integer :: nwords, j, nsections, stind, st, i

    !-------Executable----------!
    ! Warning: axis must be nullified or associated!
    call deallocate_test ( axis, "axis", ModuleName )

    !Count words in inline. 

    nwords = 1
    do j = 2, LEN(inline)
      if ( inline(j:j) /= " " .AND. inline(j-1:j-1) == " " )  nwords = nwords+1
    end do
    nsections = (nwords-3)/2

    call allocate_test ( n_levs_in_sec, nsections, "n_levs_in_sec", ModuleName )
    call allocate_test ( gridstep, nsections, "gridstep", ModuleName )
    call allocate_test ( axints, nsections*2, "axints", ModuleName )

    read ( unit=inline, fmt=* ) linetype, axistype, baseval, axints
    n_levs_in_sec = NINT(axints(1:nsections*2-1:2))
    gridstep = axints(2:nsections*2:2)

    axis_len = SUM(n_levs_in_sec)

    call allocate_test ( axis, axis_len, "axis", ModuleName )
    axis(1) = baseval
    stind = 0
    st = 2
    do j = 1, nsections
      do i = st, n_levs_in_sec(j)
        axis(stind+i) = axis(stind+i-1)+gridstep(j)
      end do
      stind = stind+n_levs_in_sec(j)
      st = 1
    end do
    call deallocate_test ( axints, "axints", ModuleName )
    call deallocate_test ( gridstep, "gridstep", ModuleName )
    call deallocate_test ( n_levs_in_sec, "n_levs_in_sec", ModuleName )
  !================================
  end subroutine MAKE_LINEAR_AXIS
  !================================

  ! ----------------------------------------------  MAKE_LOG_AXIS  -----
  !=====================================================================
  subroutine MAKE_LOG_AXIS ( inline, axis, axis_len )
  !=====================================================================
    !--------args------------!
    character(len=*),intent(in) :: inline
    REAL(kind=r8), pointer, dimension(:) :: axis
    integer, intent(out) :: axis_len
    !-------locals--------------!
    character(len=30) :: linetype, axistype
    REAL(kind=r8) :: basepressure
    integer, dimension(:), pointer :: n_levs_in_sec, n_levs_per_dec, axints
    integer :: nwords, j, nsections, stind, st, i
    REAL(kind=r8) :: gridstep

    !-------Executable----------!

    ! Warning: axis must be nullified or associated!
    call deallocate_test ( axis, "axis", ModuleName )

    !Count words in inline. 
    nwords = 1
    do j = 2, LEN(inline)
      if (inline(j:j) /= " " .AND. inline(j-1:j-1) == " ") nwords = nwords+1
    end do
    nsections = (nwords-3)/2
    call allocate_test ( n_levs_in_sec, nsections, "n_levs_in_sec", ModuleName )
    call allocate_test ( n_levs_per_dec, nsections, "n_levs_per_dec", &
      & ModuleName )
    call allocate_test ( axints, nsections*2, "axints", ModuleName )
    read ( unit=inline, fmt=* ) linetype, axistype, basepressure, axints
    n_levs_in_sec = axints(1:nsections*2-1:2)
    n_levs_per_dec = axints(2:nsections*2:2)

    axis_len = SUM(n_levs_in_sec)

    call allocate_test ( axis, axis_len, "axis", ModuleName )
    axis(1) = -LOG10(basepressure)
    stind = 0
    st = 2
    do j = 1, nsections
      gridstep = 1.0_r8/n_levs_per_dec(j)
      do i = st, n_levs_in_sec(j)
        axis(stind+i) = axis(stind+i-1)+gridstep
      end do
      stind = stind+n_levs_in_sec(j)
      st = 1
    end do
    axis = 10.0_r8**(-axis)
    call deallocate_test ( n_levs_in_sec, "n_levs_in_sec", ModuleName )
    call deallocate_test ( n_levs_per_dec, "n_levs_per_dec", ModuleName )
    call deallocate_test ( axints, "axints", ModuleName )
  !=============================
  end subroutine MAKE_LOG_AXIS
  !=============================

  ! -----------------------------------------  READ_EXPLICIT_AXIS  -----
  !=====================================================================
  subroutine READ_EXPLICIT_AXIS ( unit,axis,axis_len )
  !=====================================================================
    !--------args------------!
    integer, intent(in) :: unit
    real(kind=r8), pointer, dimension(:) :: axis
    integer, intent(out) :: axis_len

    !------- Local vars ---------!
    character(len=30) :: readitem
    character(len=1) :: rdchar
    real(kind=r8), dimension(1:200) :: tmpaxis
    integer :: i, iotest
    logical :: foundcb

    !Executables
    ! Warning: axis must be nullified or associated!
    call deallocate_test ( axis, "axis", ModuleName )

    ! An explicit axis is supplied as a parenthesised list, spread 
    ! over several lines. This is a Royal PIA.
    ! read chars till we get to the (
    do 
      read ( unit=unit, fmt="(a)", advance="no" ) rdchar

      if ( rdchar=="(" ) exit 
    end do
    ! now read items and add them to the axis until we get to the )
    axis_len = 0
    foundcb = .FALSE.
    itemsloop: do
      i=1
      charsloop: do
        read ( unit=unit, fmt="(a)", advance="no", iostat=iotest ) rdchar

        if ( rdchar==")" ) then
          foundcb = .TRUE.
          exit charsloop
        end if
        if ( rdchar==" " ) exit charsloop
        readitem(i:i) = rdchar
        i = i+1
      end do charsloop
      if ( (i<=1 .AND. .NOT.foundcb) .OR. (i==2 .AND. readitem(1:1)==" ")) then
        cycle itemsloop
      end if
      if ( i <= 1 .AND. foundcb ) exit itemsloop
      axis_len = axis_len+1
      read ( unit=readitem(1:i-1), fmt=* ) tmpaxis(axis_len)
      if (foundcb) exit itemsloop
    end do itemsloop
    call allocate_test ( axis, axis_len, "axis", ModuleName )
    axis = tmpaxis(1:axis_len)
  !==================================
  end subroutine READ_EXPLICIT_AXIS
  !==================================
!============================
end module ObtainClimatology
!============================
! $Log$
! Revision 2.0  2000/09/05 18:57:06  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry
!

