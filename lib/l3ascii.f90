! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module L3ascii ! Collections of Hugh's subroutines to handle TYPE GriddedData_T
!=============================================================================

  use GriddedData, only: GriddedData_T, v_is_pressure, v_is_altitude, &
    & v_is_gph, v_is_theta
  use LEXER_CORE, only: PRINT_SOURCE
  USE MLSCommon, only: R8, LineLen, NameLen
  USE MLSStrings, only: Capitalize, &
    & Count_words, ReadCompleteLineWithoutComments
  USE output_m, only: output
  use TREE, only: DUMP_TREE_NODE, SOURCE_REF

  implicit NONE
  public

  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), private, parameter :: IdParm = & 
    & "$Id$"
  character(len=len(idParm)), private :: Id = idParm
  character(len=*), private, parameter :: ModuleName = &
    & "$RCSfile$"
  !-----------------------------------------------------------------------------

  public::L3ascii_open, L3ascii_read_field, L3ascii_interp_field, Make_log_axis
  public::L3ascii_get_multiplier
  !private::get_next_noncomment_line, 
  private :: Make_linear_axis
  private :: Read_explicit_axis, Ilocate
  integer, private :: ERROR

   ! --------------------------------------------------------------------------

contains

  subroutine L3ascii_open ( filename, unit )
    ! opens a l3ascii file, reads, prints and discards the annoying
    ! header line and returns the unit it chose. No special close routine.
    ! just do close(unit=unit)
    !--------- argument ----------!
    character(len=*),intent(in)::filename
    integer,intent(out)::unit
    !-------locals------!
    logical:: tiedup, found
    integer:: j
    !character(len=LineLen)::headerline
    !----executables----!
    !----- Find first unattached unit -----!
    found = .false.
    error = 0
    do j = 1, 30
      inquire ( unit=j, opened=tiedup )
      if ( .not. tiedup ) then
        found = .true.
        unit = j
        open ( unit=unit, file=filename, status="old", action="read" )
        exit
      end if
    end do
    if ( .not. found ) then
      unit = -1
      call announce_error ( 0,&
           "in subroutine l3ascii_open: No units left" )
!      call MLSMessage ( MLSMSG_Error, ModuleName,&
!           "in subroutine l3ascii_open: No units left" )
    end if

    !First line is not prefaced with ; nor is it of any use. 
    ! we read it to move the file position past it
    !read(unit=unit,fmt="(a)")headerline
    !    print*,headerline
  end subroutine L3ascii_open

  subroutine L3ascii_read_field ( Unit, Field, End_of_file, ErrType )
    use dates_module    ! Shoud use SDP Toolkit eventually. 
    ! ----Arguments ----!
    integer, intent(in) :: unit
    integer, intent(out), optional :: ErrType
    type(GriddedData_T), intent(inout) :: field
    logical , intent(out) :: end_of_file
    !-------Local Variables --------!
    character(len=*),parameter :: dummyyear="1993"
    logical :: opened
    character(len=LineLen) :: inline
    character(len=30) :: linetype, axistype, sdstring, edstring
    character(len=80) :: filename, unitstring
    real(kind=r8), pointer, dimension(:) :: tmpaxis => null()
    real(kind=r8), pointer, dimension(:) :: dateStarts => null()
    real(kind=r8), pointer, dimension(:) :: dateEnds => null()
    integer :: tmpaxis_len, idate, word_count
    integer,parameter :: maxNoDates = 30
    real(kind=r8), allocatable, dimension(:,:,:,:,:,:) :: tmpfield
    !---- Executable statements ----! 
    error = 0
	 
! Abrupt termination--as with an error--means use field at own risk
    if ( present(ErrType) ) then
      ErrType = 1
    end if

!    nullify(tmpaxis)
	 end_of_file = .TRUE.	! Terminate loops based around this on error

    write(unit=unitstring,fmt="(i3)") unit ! For use in error reporting
    inquire(unit=unit,opened=opened)
    if ( .not. opened ) then
        call announce_error(0, &
       " in subroutine l3ascii_read_field, Unit "//trim(unitstring)//&
       "is not connected to a file. Do call l3ascii_open(filename,unit) first")
!        call MLSMessage(MLSMSG_Error,ModuleName,&
!       " in subroutine l3ascii_read_field, Unit "//trim(unitstring)//&
!       "is not connected to a file. Do call l3ascii_open(filename,unit) first")
       return
    end if
    inquire(unit=unit,name=filename) ! find the file name connected to this
    ! unit for use in error messages.

	field%sourceFileName = filename

    ! Fix axis arrays: set to default values with length 1 and a sensible 
    ! value. These will be used if the file does not have variation 
    ! along that axis

    ! Unfortunately, we can't tell if this defined type has
    ! been used before by testing a pointer for being associated as for 
    ! a new  struct they are in the undefined state. (This will change with 
    ! Fortran 95, where you can define a pointer to be initialised to the 
    ! nullified state)
    ! We therefore use a special number field%reusing to test whether 
    ! this defined type has been used before.  
!    if ( field%reusing==313323435 ) then
    if ( associated(field%field) ) then
       !       print*,"This struct has been used before : deallocating"
       deallocate(field%heights)
       deallocate(field%lats)
       deallocate(field%lons)
       deallocate(field%lsts)
       deallocate(field%szas)
       deallocate(field%dateStarts)
       deallocate(field%dateEnds)
       deallocate(field%field)
!    else
       !       print*,"This is a new struct"
!       field%reusing=313323435
    end if
! Automatically create a stub grid template with minimal size
! Each component will be deallocated && reallocated with correct sizes later
    allocate(field%heights(1:1))
    field%heights(1)=1000.0
    field%noHeights=1
    field%verticalCoordinate=1
    allocate(field%lats(1:1))
    field%lats(1)=0.0
    field%noLats=1
    field%equivalentLatitude=.false.
    allocate(field%lons(1:1))
    field%lons(1)=0.0
    field%noLons=1
    allocate(field%lsts(1:1))
    field%lsts(1)=12.0
    field%noLsts=1
    allocate(field%szas(1:1))
    field%szas(1)=30.0
    field%noSzas=1
    allocate(field%dateStarts(1:1))
    field%dateStarts(1)=30.0
    field%noDates=1
    allocate(field%dateEnds(1:1))
    field%dateStarts(1)=30.0
    field%noDates=1
    allocate(field%field(1:1,1:1,1:1,1:1,1:1,1:1))
    field%dateStarts(1)=30.0
    field%noDates=1
    ! Dates are mandatory, so we don't have to give them a default value

    !--- Read field info and all the axis info ---!
!    call get_next_noncomment_line(unit,inline)
    end_of_file=.false.
    call ReadCompleteLineWithoutComments ( unit, inline, eof=end_of_file )
!    print*,"Read line"
!    print*,inline
    if ( Capitalize(inline(1:5)) /= "FIELD" ) then
      call announce_error ( 0, &
        & "in subroutine l3ascii_read_field, File "//trim(filename)// &
        &  "on unit"//trim(unitstring)//" contains no more Fields" )
!      call MLSMessage(MLSMSG_Error,ModuleName,&
!           "in subroutine l3ascii_read_field, File "//trim(filename)// &
!           "on unit"//trim(unitstring)//" contains no more Fields")
	     end_of_file=.true.
      return
    end if
    
    if ( end_of_file ) then
      call announce_error ( 0,&
        & "In subroutine l3ascii_read_field, End of File"//trim(filename)// &
        & " on unit"//trim(unitstring))
!       call MLSMessage(MLSMSG_Error,ModuleName,&
!       "In subroutine l3ascii_read_field, End of File"//trim(filename)// &
!       " on unit"//trim(unitstring))
      return
    end if

    read ( unit=inline, fmt=* ) linetype, field%quantityName, &
         field%description, field%units
axesloop:do

      call ReadCompleteLineWithoutComments ( unit, inline )
      !print*,inline
      read ( unit=inline, fmt=* ) linetype, axistype
      linetype=Capitalize(linetype)
      axistype=Capitalize(axistype)
      if ( linetype(1:4) == "DATE" ) then ! This is always the last "axis"
        exit axesloop                 ! and is different from the others
      end if
      if ( axistype(1:6) =="LINEAR" ) then
!        print*,"Doing linear axis"
        call make_linear_axis ( inline, tmpaxis, tmpaxis_len )
!        print*,"Done linear axis"
      else if ( axistype(1:3) =="LOG" ) then
        !print*,"Doing log axis"
        call make_log_axis ( inline, tmpaxis, tmpaxis_len )
        !print*,"Done log axis"
      else if ( axistype(1:8) =="EXPLICIT" ) then
!       print*,"Doing explicit axis"
        backspace(unit=unit)
        call read_explicit_axis ( unit, tmpaxis, tmpaxis_len )
!         print*,"Done explicit axis"
      else
                   end_of_file=.true.
         call announce_error(0,&
              "in subroutine l3ascii_read_field,File"//trim(filename)//&
              " on unit"//trim(unitstring)//" contains coordinate"//&
              " of invalid type "//trim(axistype)//"for axis"//&
              trim(linetype))
!         call MLSMessage(MLSMSG_Error,ModuleName,&
!              "in subroutine l3ascii_read_field,File"//trim(filename)//&
!              " on unit"//trim(unitstring)//" contains coordinate"//&
!              " of invalid type "//trim(axistype)//"for axis"//&
!              trim(linetype))
         return
      end if

      ! I do not entirely grok what NJL intended verticalCoordinate to be.
      if ( linetype(1:8) == "PRESSURE" .or. linetype(1:8) == "ALTITUDE" &
        &  .or. linetype(1:3) == "GPH" .or. linetype(1:5) == "THETA" ) then
        field%noHeights = tmpaxis_len
        deallocate(field%heights)
        allocate(field%heights(1:tmpaxis_len))
        field%heights = tmpaxis

        if ( linetype(1:8) == "PRESSURE" ) then
           field%verticalCoordinate = v_is_pressure ! 1
        else if ( linetype(1:8) == "ALTITUDE" ) then
           field%verticalCoordinate = v_is_altitude    ! 2
        else if ( linetype(1:3) == "GPH" ) then
           field%verticalCoordinate = v_is_gph ! 3
        else if ( linetype(1:5) == "THETA" ) then
           field%verticalCoordinate = v_is_theta ! 4
        end if
      else if ( linetype(1:8) == "LATITUDE" .or. &
        &  linetype(1:8) == "EQUIVLAT" ) then
        field%noLats = tmpaxis_len
        deallocate(field%lats)
        allocate(field%lats(1:tmpaxis_len))
        field%lats=tmpaxis
        if ( linetype(1:8) == "LATITUDE" ) then
           field%equivalentLatitude = .false.
        else
           field%equivalentLatitude = .true.
        end if
      else if ( linetype(1:9) == "LONGITUDE" ) then
        field%noLons = tmpaxis_len
        deallocate(field%lons)
        allocate(field%lons(1:tmpaxis_len))
        field%lons = tmpaxis
      else if ( linetype(1:9) == "LST" ) then
        field%noLsts = tmpaxis_len
        deallocate(field%lsts)
        allocate(field%lsts(1:tmpaxis_len))
        field%lsts=tmpaxis
      else if ( linetype(1:9) == "SZA" ) then
        field%noSzas = tmpaxis_len
        deallocate(field%szas)
        allocate(field%szas(1:tmpaxis_len))
        field%szas = tmpaxis
      end if
    end do axesloop
    deallocate(tmpaxis)

    ! We already have the first date line read and the axis type extracted
    ! It was the existence of a date line that caused us to exit from 
    ! the loop axesloop. We don't know how many dates there are so we have to 
    ! allocate large arrays and copy their contents to an array of the 
    ! right size 
    field%noDates=1
    allocate(tmpfield(1:field%noHeights,1:field%noLats,1:field%noLons, &
         1:field%noLsts,1:field%noSzas,1:maxNoDates))
    allocate(dateStarts(1:maxNoDates),dateEnds(1:maxNoDates))

    ! Loop to read in the data for the current date and check to see if 
    ! there is another date
datesloop: do idate = 1, maxNoDates
      !print*,"Datesloop: idate=",idate
      word_count = count_words(inline)
      if ( word_count == 3 ) then
        read(unit=inline,fmt=*)linetype,axistype,sdstring
        sdstring = adjustl(sdstring)
        edstring = sdstring
      else if ( word_count >= 4 ) then
        read ( unit=inline, fmt=* ) linetype, axistype, sdstring, edstring
        sdstring = adjustl(sdstring)
        edstring  =adjustl(edstring)
      else
        call announce_error ( 0, &
          & "in subroutine l3ascii_read_field: File"//trim(filename)//&
          & "on unit"//trim(unitstring)//" contains a line beginning"//&
          & trim(linetype)//"Date with too few words " )
!         call MLSMessage(MLSMSG_Error,ModuleName,&
!              "in subroutine l3ascii_read_field: File"//trim(filename)//&
!              "on unit"//trim(unitstring)//" contains a line beginning"//&
!              trim(linetype)//"Date with too few words ")
        end_of_file = .true.
        return
      end if
 
      ! Date strings can begin with - indicating the year is
      ! missing and that the file belongs to no year in particular.
      ! To convert dates to SDP toolkit (TAI) times (Seconds since start of
      ! 1 Jan 1993) we need to stick on a dummy year
      if ( sdstring(1:1) == "-" ) then
        sdstring=dummyyear//sdstring
      end if
      if ( edstring(1:1) == "-" ) then
        edstring=dummyyear//edstring
      end if
      ! ccsds2tai returns days since 1 Jan 1993. 86400==no of secs per day
      dateStarts(idate)=86400*ccsds2tai(sdstring)
      dateEnds(idate)  =86400*ccsds2tai(edstring)
      call ReadCompleteLineWithoutComments ( unit, inline )
      backspace ( unit=unit )
      !print*,"About to read data: inline=",inline
      !print*,"tmpfield has size:",size(tmpfield),shape(tmpfield)
 
      read ( unit=unit, fmt=* ) tmpfield(:,:,:,:,:,idate)
      !print*,"Read data"
      end_of_file = .false.
      call ReadCompleteLineWithoutComments(unit,inline,eof=end_of_file)
      !print*,"Next date line:",inline,"EOF=",end_of_file
      if ( end_of_file ) then
        ! No more dates and nothing else either
        exit datesloop
      end if
      read ( unit=inline, fmt=* ) linetype, axistype
      linetype = Capitalize(linetype)
      axistype = Capitalize(axistype)
      if ( end_of_file .or. &
        & linetype(1:5) == "FIELD" .or. linetype(1:3)=="END" ) then
        !Oops! There were no more dates, but there is another field
        backspace ( unit=unit )
        exit datesloop
      end if
      if ( linetype(1:4) /= "DATE" ) then ! There should be another date here
        call announce_error ( 0, &
          & "in subroutine l3ascii_read_field: File"//trim(filename)//&
          & "on unit"//trim(unitstring)//" contains a line beginning"//&
          & trim(linetype)//"where I expected a line beginning Date " )
!         call MLSMessage(MLSMSG_Error,ModuleName,&
!              "in subroutine l3ascii_read_field: File"//trim(filename)//&
!              "on unit"//trim(unitstring)//" contains a line beginning"//&
!              trim(linetype)//"where I expected a line beginning Date ")
        end_of_file = .true.
        return
      end if
      field%noDates = field%noDates+1

    end do datesloop

    deallocate(field%dateStarts, field%dateEnds, field%field)
!
    allocate(field%field(1:field%noHeights,1:field%noLats,&
      & 1:field%noLons,1:field%noLsts,1:field%noSzas,1:field%noDates))
    allocate(field%dateStarts(1:field%noDates),&
      & field%dateEnds(1:field%noDates))
    field%dateStarts = dateStarts(1:field%noDates)
    field%dateEnds = dateEnds(1:field%noDates)
    field%field = tmpfield(:,:,:,:,:,1:field%noDates)
    deallocate(tmpfield,dateStarts,dateEnds)

! Normal termination--assume field is valid maybe even correct
    if ( present(ErrType) ) then
      ErrType = 0
    end if
  end subroutine L3ascii_read_field

  subroutine L3ascii_interp_field ( field, outval, pressure, lat, lon, lst, &
    & sza, date )
    ! Returns a value in outval containing the value of the 
    ! gridded data set "field" at the pressure, lat, etc specified by 
    ! the other args. Co-ordinates are 
    ! all optional and suitably dull defaults are chosen if no 
    ! arg is supplied. The date is supplied in tai format i.e. seconds since
    ! Midnite, 1 Jan 1993. 
    ! At the moment the height coord has to be pressure.
    !--------------Arguments------------!
    type(GriddedData_T),intent(in)::field
    real(kind=r8),intent(in),optional::pressure
    real(kind=r8),intent(in),optional::lat,lon,lst,sza,date
    real(kind=r8),intent(out)::outval
    !---- local vars: optional arg values ------!
    real(kind=r8):: inlat,inlon,inlst,insza,indate,inpressure,inalt
    !---- local vars: others--------!
    integer:: ilat1,ilat2,ilon1,ilon2,isza1,isza2,ilst1,ilst2,idate1,idate2
    integer:: ialt1,ialt2,i0,i1,i2,i3,i4
    integer,dimension(1:6)::hcshape 
    real(kind=r8),pointer,dimension(:)::tmpalt,tmpdate
    real(kind=r8),allocatable,dimension(:,:,:,:,:,:)::hcube

    !----Executable code ---- !

    ! --- Set default values for all optional parameters
    if ( present(lat) ) then
      inlat = lat
    else
      inlat = 30.0
    end if

    if ( present(pressure) ) then
      inpressure = pressure
    else
      inpressure = 10.0
    end if

    if ( present(lon) ) then
      inlon = lon
    else
      inlon = 0.0
    end if

    if ( present(sza) ) then
      insza = sza
    else
      insza = 30.0
    end if

    if ( present(lst) ) then
      inlst = lst
    else
      inlst = 12.0
    end if

    if ( present(date) ) then
      indate = date
    else
      indate = 15768000.0 ! This should be June (ish)
    end if

    ! OK, that's the optional arguments sorted. 
    ! Now, we have to find which two indices the chosen value lies between
    ! for each of the six values.
    call ilocate ( field%lats, inlat, ilat1, ilat2 )
    call ilocate ( field%lons, inlon, ilon1, ilon2 )
    call ilocate ( field%szas, insza, isza1, isza2 )
    call ilocate ( field%lsts, inlst, ilst1, ilst2 )
    ! We can't interpolate dates and prssures as they are. We construct
    ! a mean date and a log pressure
    allocate(tmpdate(1:field%noDates),tmpalt(1:field%noHeights))
    tmpdate = (field%dateStarts+field%dateEnds)/2.0
    tmpalt = -log10(field%heights)
    inalt = -log10(inpressure)
    call ilocate ( tmpdate, indate, idate1, idate2 )
    call ilocate ( tmpalt, inalt, ialt1, ialt2 )
    ! We now know that the desired point is inside the 6-D hypercuboid 
    ! bin of the array field%field(ialt1:ialt2, ilat1:ilat2, .......)
    ! Time to do the actual interpolation. Yuck. 
    ! There are up to 2^6=64 values
    ! to be taken into account. 

    ! Lets try allocating a cube of the right shape.
    if ( ialt1 == ialt2 ) then !test whether the altitude is between two of the 
      hcshape(1) = 1       ! grid values (so we interpolate) or if there
    else                   ! is only one value (so we just return that value)
      hcshape(1) = 2
    end if
    if ( ilat1 == ilat2 ) then
      hcshape(2) = 1
    else
      hcshape(2) = 2
    end if
    if ( ilon1 == ilon2 ) then
      hcshape(3) = 1
    else
      hcshape(3) = 2
    end if
    if ( ilst1 == ilst2 ) then
      hcshape(4) = 1
    else
      hcshape(4) = 2
    end if
    if ( isza1 == isza2 ) then
      hcshape(5) = 1
    else
      hcshape(5) = 2
    end if
    if ( idate1 == idate2 ) then
      hcshape(6) = 1
    else
      hcshape(6) = 2
    end if


    ! Allocate a 6-D hypercube that our point lies in. Some of the 
    ! dimensions may be "collapsed" 
    allocate(hcube(1:hcshape(1),1:hcshape(2),1:hcshape(3), &
      & 1:hcshape(4),1:hcshape(5),1:hcshape(6)))
    ! Copy data into hypercube
    hcube = field%field(ialt1:ialt2, ilat1:ilat2, ilon1:ilon2, &
      & ilst1:ilst2,isza1:isza2, idate1:idate2)

    ! Now we interpolate along each of the axes where this is needed
    ! Reduce 6-d to 5-d
    if ( hcshape(6) == 2 ) then !Interpolate in date 
      do i0 = 1, hcshape(1)
        do i1 = 1, hcshape(2)
          do i2 = 1, hcshape(3)
            do i3 = 1, hcshape(4)
              do i4 = 1,hcshape(5)
                hcube(i0,i1,i2,i3,i4,1) = hcube(i0,i1,i2,i3,i4,1) + &
                  &  (hcube(i0,i1,i2,i3,i4,2) - &
                  &  hcube(i0,i1,i2,i3,i4,1)) * &
                  &  (indate-tmpdate(idate1)) / &
                  &  (tmpdate(idate2)-tmpdate(idate1))
              end do
            end do
          end do
        end do
      end do
    end if
    ! Reduce  5-d to 4d
    if ( hcshape(5) == 2 ) then !Interpolate in local Solar zenith ang  SZA
      do i0 = 1, hcshape(1)
        do i1 = 1, hcshape(2)
          do i2 = 1, hcshape(3)
            do i3 = 1, hcshape(4)
              hcube(i0,i1,i2,i3,1,1) = hcube(i0,i1,i2,i3,1,1) + &
                &  (hcube(i0,i1,i2,i3,2,1) - &
                &  hcube(i0,i1,i2,i3,1,1)) * &
                &  (insza-field%szas(isza1)) / &
                &  (field%szas(isza2)-field%szas(isza1))
            end do
          end do
        end do
      end do
    end if
    ! Reduce  4-d to 3d
    if ( hcshape(4) == 2 ) then !Interpolate in local solar time
      do i0 = 1, hcshape(1)
        do i1 = 1, hcshape(2)
          do i2 = 1, hcshape(3)
            hcube(i0,i1,i2,1,1,1) = hcube(i0,i1,i2,1,1,1) + &
              &  (hcube(i0,i1,i2,2,1,1) - &
              &  hcube(i0,i1,i2,1,1,1)) * &
              &  (inlst-field%lsts(ilst1)) / &
              &  (field%lsts(ilst2)-field%lsts(ilst1))
          end do
        end do
      end do
    end if
    ! Reduce  3-d to 2d
    if ( hcshape(3) == 2 ) then !Interpolate in longitude
      do i0 = 1, hcshape(1)
        do i1 = 1, hcshape(2)
          hcube(i0,i1,1,1,1,1) = hcube(i0,i1,1,1,1,1) + &
            &  (hcube(i0,i1,2,1,1,1) - &
            &  hcube(i0,i1,1,1,1,1)) * &
            &  (inlon-field%lons(ilon1)) / &
            &  (field%lons(ilon2)-field%lons(ilon1))
        end do
      end do
    end if
    ! Reduce  2-d to 1d
    if ( hcshape(2) == 2 ) then !Interpolate in latitude
      do i0 = 1, hcshape(1)
        hcube(i0,1,1,1,1,1) = hcube(i0,1,1,1,1,1) + &
          &  (hcube(i0,2,1,1,1,1) - &
          &  hcube(i0,1,1,1,1,1)) * &
          &  (inlat-field%lats(ilat1)) / &
          &  (field%lats(ilat2)-field%lats(ilat1))
      end do
    end if

    ! Reduce  1-d to 0d
    if ( hcshape(1) == 2 ) then !interpolate in altitude
      hcube(1,1,1,1,1,1) = hcube(1,1,1,1,1,1) + &
        &  (hcube(2,1,1,1,1,1) - &
        &  hcube(1,1,1,1,1,1)) * &
        &  (inalt-tmpalt(ialt1)) / &
        &  (tmpalt(ialt2)-tmpalt(ialt1))
    end if
    outval = hcube(1,1,1,1,1,1)
    deallocate(tmpdate,tmpalt,hcube)
  end subroutine L3ascii_interp_field

  subroutine Ilocate ( x, xval, ix1, ix2 )
    ! This finds which two elements of x lie on either side of xval
    ! x is assumed to be 1-based. 
    !     *** arguments *** 
    real(kind=r8),pointer,dimension(:):: x
    real(kind=r8),intent(in)::xval
    integer, intent(out)::ix1,ix2
    !     *** other variables ***                                           
    integer ::j,dj,n 
    !     *** executable statements ***                                     
    n=size(x)
    if ( n <= 1 ) then ! x has only one element: no interpolating to be done
      ix1 = 1
      ix2 = 1
    else
      if ( xval > x(n) ) then ! off top . Use end value
        ix1 = n
        ix2 = n
      else if ( xval < x(1) ) then ! of bottom. Use end value
        ix1 = 1
        ix2 = 1
      else ! in range of x. Do binary search.
        j = n/2 
        dj = n/2 
binsearch: do  
          if ( dj > 1 ) then 
            dj = dj/2 
          else 
            dj = 1
          end if
          if ( x(j) > xval .and. x(j+1) > xval ) then 
            j = j-dj 
            cycle
          else  if ( x(j) < xval .and. x(j+1) < xval ) then 
            j = j+dj 
            cycle 
          else 
            exit
          end if
        end do binsearch
        ix1 = j
        ix2 = j+1
      end if
    end if
  end subroutine Ilocate

!  subroutine get_next_noncomment_line(unit,line)
!    !---Arguments----!
!    integer,intent(in)::unit
!    character(len=*),intent(out)::line
!    !---------Local vars------!
!    integer:: ioinfo
!    !---Executable bit -----!
!    line(1:1)=" "
!    ioinfo=0
!    rdloop:do
!       read(unit=unit,fmt="(a)",iostat=ioinfo)line      !read a line
!       if ( ioinfo /= 0 ) then
!          line="End of File Found"
!          exit rdloop
!       else
!          line=adjustl(line)                 ! remove blanks from start
!          if ( line(1:1) /= ";" .and. line(1:1) /= " " ) then !not a comment
!             exit rdloop                     ! so exit loop and return line
!          end if
!       end if
!       !        print*,line(1:80)
!    end do rdloop

    !    print*,"Got non-comment line"
    !    print*,line(1:80)

!  end subroutine get_next_noncomment_line

  subroutine make_log_axis ( inline, axis, axis_len )
    !--------args------------!
    character(len=*),intent(in)::inline
    real(kind=r8),pointer,dimension(:)::axis
    integer,intent(out)::axis_len
    !-------locals--------------!
    character(len=30)::linetype,axistype
    real(kind=r8)::basepressure
    integer,dimension(:),allocatable::n_levs_in_sec,n_levs_per_dec,axints
    integer::nwords,j,nsections,stind,st,i
    real(kind=r8)::gridstep
    !-------Executable----------!

    ! Warning: axis must be nullified or associated!
    if ( associated(axis) ) then 
       deallocate(axis)
    end if

    if ( len(inline) <= 1 ) then
      call announce_error ( 0, &
        & "in make_log_axis: inline, <" // trim(inline)//">, too short" )
    end if
	
    !Count words in inline. 
    nwords = 1
    do j = 2, len(inline)
      if ( inline(j:j) /= " " .and. inline(j-1:j-1) == " " ) then
        nwords = nwords+1
      end if
    end do
    nsections = (nwords-3)/2
!    print*,"Inline=",inline
!    print*,"Nwords=",nwords

    if ( nsections < 1 ) then
      call announce_error ( 0, "in make_log_axis: nsections < 1" )
    end if
	
    allocate(n_levs_in_sec(1:nsections),n_levs_per_dec(1:nsections),&
      &  axints(1:nsections*2))
    read ( unit=inline, fmt=* ) linetype, axistype, basepressure, axints
!    print*,"Linetype=",linetype," axistype=",axistype
!    print*,"basepressure=",basepressure," Axints=",axints
    n_levs_in_sec = axints(1:nsections*2-1:2)
    n_levs_per_dec = axints(2:nsections*2:2)

    axis_len = sum(n_levs_in_sec)

    if ( axis_len < 1 ) then
      call announce_error ( 0, "in make_log_axis: axis_len < 1" )
    end if
	
!    print*,"axis length=", axis_len
    allocate(axis(1:axis_len))
    axis(1) = -log10(basepressure)
    stind = 0
!    print*,"nsections=", nsections
    do j = 1, nsections
      gridstep=1.0_r8/n_levs_per_dec(j)
!    print*,"j=", j
!    print*,"n_levs_per_dec(j)=", n_levs_per_dec(j)
!    print*,"gridstep=", gridstep
      if ( gridstep <= 0.d0 ) then
        call announce_error ( 0, "in make_log_axis: gridstep <= 0" )
        stop
      end if
      if ( j == 1 ) then
        st = 2
      else
        st = 1
      end if
      do i = st, n_levs_in_sec(j)
        axis(stind+i)=axis(stind+i-1)+gridstep
      end do
      stind = stind+n_levs_in_sec(j)
    end do

!    call dump ( axis, &
!          & '  log(log axis) =' )

! not sure why this doesn't always work, but it doesn't for paw
!    axis=10.0_r8**(-axis)
! (possibly a compiler bug for NAG on Linux)
! so instead we'll use the equivalent:
    axis = exp(-log(10.)*axis)
!

    deallocate(n_levs_in_sec, n_levs_per_dec,axints)

!        call dump ( axis, &
!          & '  log axis =' )
  end subroutine Make_log_axis

  subroutine Make_linear_axis ( inline, axis, axis_len )
    !--------args------------!
    character(len=*),intent(in)::inline
    real(kind=r8),pointer,dimension(:)::axis
    integer,intent(out)::axis_len
    !-------locals--------------!
    character(len=30)::linetype,axistype
    real(kind=r8)::baseval
    integer,dimension(:),allocatable::n_levs_in_sec
    real(kind=r8),dimension(:),allocatable::gridstep,axints
    integer::nwords,j,nsections,stind,st,i
    !-------Executable----------!
    ! Warning: axis must be nullified or associated!
    if ( associated(axis) ) then 
      deallocate(axis)
    end if

    !Count words in inline. 

    nwords = 1
    do j = 2,len(inline)
      if ( inline(j:j) /= " " .and. inline(j-1:j-1) == " " ) then
        nwords = nwords+1
      end if
    end do
    nsections = (nwords-3)/2

    allocate(n_levs_in_sec(1:nsections),gridstep(1:nsections),&
      & axints(1:nsections*2))

    read ( unit=inline, fmt=* ) linetype, axistype, baseval, axints
    n_levs_in_sec = nint(axints(1:nsections*2-1:2))
    gridstep = axints(2:nsections*2:2)

    axis_len = sum(n_levs_in_sec)

    allocate(axis(1:axis_len))
    axis(1) = baseval
    stind = 0
    do j = 1, nsections
      if ( j == 1 ) then
        st = 2
      else
        st = 1
      end if
      do i = st, n_levs_in_sec(j)
        axis(stind+i) = axis(stind+i-1)+gridstep(j)
      end do
      stind = stind+n_levs_in_sec(j)
    end do
    deallocate(axints,gridstep,n_levs_in_sec)

  end subroutine Make_linear_axis

  subroutine Read_explicit_axis ( unit, axis, axis_len )
    !--------args------------!
    integer,intent(in)::unit
    real(kind=r8),pointer,dimension(:)::axis
    integer,intent(out)::axis_len
    !------- Local vars ---------!
    integer,parameter::ri_len=30
    character(len=ri_len)::readitem
    character(len=1)::rdchar
    real(kind=r8),dimension(1:200)::tmpaxis
    integer::i,iotest
    logical::foundcb
    !Executables
    ! Warning: axis must be nullified or associated!
    if ( associated(axis) ) then 
      deallocate(axis)
    end if
    ! readitem needs to be initialised or there is a point where it 
    ! can be used before being set.
    readitem=""
    ! An explicit axis is supplied as a parenthesised list, spread 
    ! over several lines. This is a Royal PIA.
    ! Read chars till we get to the (
    do 
      read ( unit=unit, fmt="(a)", advance="no" ) rdchar
       !        write(unit=*,fmt="(a)",advance="no")rdchar 
      if ( rdchar=="(" ) then
        exit 
      end if
    end do
    !    print*,"Got open paren"
    ! now read items and add them to the axis until we get to the )
    axis_len = 0
    foundcb = .false.
itemsloop:do
      i = 1
  charsloop:do
        read ( unit=unit, fmt="(a)", advance="no", iostat=iotest ) rdchar
          !write(unit=*,fmt="(a)",advance="no")rdchar 
        if ( rdchar == ")" ) then
             !print*,"Found ) at end of explicit axis"
          foundcb = .true.
          exit charsloop
        end if
        if ( rdchar = =" " ) then
           exit charsloop
        end if
        readitem(i:i) = rdchar
        i = i+1
      end do charsloop
      if ( (i<=1 .and. .not.foundcb) .or. (i==2 .and. readitem(1:1)==" ") ) then
        cycle itemsloop
      end if
      if ( i <= 1 .and. foundcb ) then
        exit itemsloop
      end if
      !print*,"Axis item is",readitem(1:i-1)
      axis_len = axis_len+1
      read ( unit=readitem(1:i-1), fmt=* ) tmpaxis(axis_len)
      !print*,"Element",axis_len," is ",tmpaxis(axis_len)
      if ( foundcb ) then
        exit itemsloop
      end if
    end do itemsloop
    allocate(axis(1:axis_len))
    axis = tmpaxis(1:axis_len)
  end subroutine Read_explicit_axis

  function L3ascii_get_multiplier ( field ) result ( multiplier )
    ! This function attempts to return the number by which the mixing ratio
    ! in the file was multiplied. If the data are mixing ratio or Kelvin, it 
    ! returns 1, if they are ppmv, it returns 1.e6, if they are ppbv, it
    ! returns 1.e9, if they are pptv it returns 1.e12. This is  just done by
    ! parsing the "units" string in the file. If that is wrong, the function
    ! won't be right either.
    !------Argument)------!
    type(GriddedData_T), intent(in) :: field
    !---Function result-----!
    real(kind=r8) :: multiplier
    !---other vars-------!
    character(len=NameLen) :: ucunits
    integer :: ix
    !----Executable functions---!
    error = 0
    ucunits = Capitalize(field%units)
    ix = index(ucunits,"VMR ")
    if ( ix > 0 ) then
      multiplier = 1.0_r8
      return
    end if
    ix = index(ucunits,"PPM")
    if ( ix > 0 ) then
      multiplier = 1.0e6_r8
      return
    end if
    ix = index(ucunits,"PPB")
    if ( ix > 0 ) then
      multiplier = 1.0e9_r8
      return
    end if
    ix = index(ucunits,"PPT")
    if ( ix > 0 ) then
      multiplier = 1.0e12_r8
      return
    end if
    ix = index(ucunits,"K ")
    if ( ix > 0 ) then
      multiplier = 1.0_r8
      return
    end if
!    call MLSMessage(MLSMSG_Warning,ModuleName,&
!         "in function l3ascii_get_multiplier: Units "//&
!         trim(field%units)//" for field "//trim(field%quantityName)//&
!         "not known. Guessing multiplier=1.0")
    call announce_error ( 0, &
      & "in function l3ascii_get_multiplier: Units "// &
      & trim(field%units)//" for field "//trim(field%quantityName)// &
      & "not known. Guessing multiplier=1.0" )
    !print*,"in function l3ascii_get_multiplier: Units "//&
    !     trim(field%units)//" for field "//trim(field%quantityName)//&
    !     "not known. Guessing multiplier=1.0"
   multiplier = 1.0_r8
  end function L3ascii_get_multiplier

  ! ------------------------------------------------  announce_error  -----
  subroutine Announce_error ( lcf_where, full_message, use_toolkit, &
  & error_number )
  
   ! Arguments
	
    integer, intent(in)    :: lcf_where
    character(LEN=*), intent(in)    :: full_message
    logical, intent(in), optional :: use_toolkit
    integer, intent(in), optional    :: error_number
    ! Local
!    character (len=80) :: msg, mnemonic
!    integer :: status
    logical :: just_print_it
    logical, parameter :: default_output_by_toolkit = .true.
 
    if ( present(use_toolkit) ) then
      just_print_it = use_toolkit
    else if ( default_output_by_toolkit ) then
      just_print_it = .false.
    else
      just_print_it = .true.
    end if
 
    if ( .not. just_print_it ) then
!      CALL Pgs_smf_getMsg(status, mnemonic, msg)
!      CALL MLSMessage (level, ModuleName, &
!                &trim(full_message)//" "//mnemonic//" "//msg)
      error = max(error,1)
      call output ( '***** At ' )

      if ( lcf_where > 0 ) then
          call print_source ( source_ref(lcf_where) )
      else
        call output ( '(no lcf node available)' )
      end if

      call output ( ': ' )
      call output ( "The " );
      if ( lcf_where > 0 ) then
        call dump_tree_node ( lcf_where, 0 )
      else
        call output ( '(no lcf tree available)' )
      end if

      call output ( " Caused the following error:", advance='yes', &
        & from_where=ModuleName )
      call output ( trim(full_message), advance='yes', &
        & from_where=ModuleName )
      if ( present(error_number) ) then
        call output ( 'error number ', advance='no' )
        call output ( error_number, places=9, advance='yes' )
      end if
    else
      print*, '***Error in module ', ModuleName
      print*, trim(full_message)
      if ( present(error_number) ) then
        print*, 'error number ', error_number
      end if
    end if

!===========================
  end subroutine announce_error
!===========================

!=============================================================================
END MODULE L3ascii
!=============================================================================

!
! $Log$
! Revision 2.7  2001/04/12 22:56:11  vsnyder
! Improve an error message, cosmetic changes
!
! Revision 2.6  2001/03/30 00:25:20  pwagner
! Fills sourceFileName
!
! Revision 2.5  2001/03/29 00:51:35  pwagner
! make_log_axis now always works
!
! Revision 2.4  2001/03/28 00:24:38  pwagner
! Some error controls, ErrType added
!
! Revision 2.3  2001/03/27 17:33:30  pwagner
! announce_error replaces MLSMessage
!
! Revision 2.2  2001/03/15 21:40:30  pwagner
! Eliminated unused routines from USE statements
!
! Revision 2.1  2001/03/15 21:27:26  pwagner
! Moved l3ascii methods from GriddedData here
!
!
