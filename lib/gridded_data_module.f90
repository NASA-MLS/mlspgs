! Copyright (c) 1999 the University of Edinburgh. All rights reserved 
! U.K. Government funding under NERC contract F14/6/33 acknowledged. 
!(If Caltech and NASA do this, we're going to as well, so there.)
! U.S. Government Sponsorship under NASA Contract NAS7-1407 also acknowledged.

module gridded_data_module
  use MLSCommon
  use MLSStrings 
  use MLSMessageModule
  implicit none
  ! This module provides subroutines to access climatology files
  ! in the ascii format described by NJL.
  private

  public::l3ascii_open,l3ascii_read_field,l3ascii_interp_field,make_log_axis
  public::l3ascii_get_multiplier
  !private::get_next_noncomment_line, 
  private :: make_linear_axis
  private::read_explicit_axis,ilocate

  character(len=256),parameter,private :: Id = &
       "$Id$"
  character(len=*),parameter,private::ModuleName=&
       "$RCSfile$"

  ! This type contains the entire contents of *one* field from an 
  ! l3ascii file The several fields in a file can be for different 
  ! products and have different numbers of axes, so we deal with 
  ! them quite separately
  type,public:: GriddedData_T
     character(len=NameLen)::quantityName,units
     character(len=LineLen)::description
     integer::verticalCoordinate,noHeights
     real(kind=r8),pointer,dimension(:)::heights
     logical::equivalentLatitude
     integer::noLats
     real(kind=r8),pointer,dimension(:)::lats
     integer::noLons
     real(kind=r8),pointer,dimension(:)::lons
     integer::noLsts
     real(kind=r8),pointer,dimension(:)::lsts
     integer::noSzas
     real(kind=r8),pointer,dimension(:)::szas
     integer::noDates
     real(kind=r8),pointer,dimension(:)::dateStarts
     real(kind=r8),pointer,dimension(:)::dateEnds

     real,pointer,dimension(:,:,:,:,:,:)::field
     integer::reusing
  end type GriddedData_T

contains

  subroutine l3ascii_open(filename,unit)
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
    do j = 1, 30
      inquire ( unit=j, opened=tiedup )
      if (.not. tiedup) then
        found = .true.
        unit = j
        open ( unit=unit, file=filename, status="old", action="read" )
        exit
      end if
    end do
    if ( .not. found ) then
      unit = -1
      call MLSMessage ( MLSMSG_Error, ModuleName,&
           "in subroutine l3ascii_open: No units left" )
    end if

    !First line is not prefaced with ; nor is it of any use. 
    ! we read it to move the file position past it
    !read(unit=unit,fmt="(a)")headerline
    !    print*,headerline
  end subroutine l3ascii_open

  subroutine l3ascii_read_field(unit,field)
    use dates_module    ! Shoud use SDP Toolkit eventually. 
    ! ----Arguments ----!
    integer, intent(in) :: unit
    type(GriddedData_T), intent(inout) :: field
    !-------Local Variables --------!
    character(len=*),parameter :: dummyyear="1993"
    logical :: opened
    character(len=LineLen) :: inline
    character(len=30) :: linetype, axistype, sdstring, edstring
    character(len=80) :: filename, unitstring
    real(kind=r8), pointer, dimension(:) :: tmpaxis, dateStarts, dateEnds
    integer :: tmpaxis_len, idate, word_count
    integer,parameter :: maxNoDates = 30
    real(kind=r8), allocatable, dimension(:,:,:,:,:,:) :: tmpfield
    logical :: end_of_file
    !---- Executable statements ----! 
    nullify(tmpaxis)

    write(unit=unitstring,fmt="(i3)") unit ! For use in error reporting
    inquire(unit=unit,opened=opened)
    if (.not. opened) then
        call MLSMessage(MLSMSG_Error,ModuleName,&
       " in subroutine l3ascii_read_field, Unit "//trim(unitstring)//&
       "is not connected to a file. Do call l3ascii_open(filename,unit) first")
       return
    endif
    inquire(unit=unit,name=filename) ! find the file name connected to this
    ! unit for use in error messages.

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
    if (field%reusing==313323435) then
       !       print*,"This struct has been used before : deallocating"
       deallocate(field%heights)
       deallocate(field%lats)
       deallocate(field%lons)
       deallocate(field%lsts)
       deallocate(field%szas)
       deallocate(field%dateStarts)
       deallocate(field%dateEnds)
       deallocate(field%field)
    else
       !       print*,"This is a new struct"
       field%reusing=313323435
    endif
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
    ! Dates are mandatory, so we don't have to give them a default value

    !--- Read field info and all the axis info ---!
!    call get_next_noncomment_line(unit,inline)
    end_of_file=.false.
    call ReadCompleteLineWithoutComments(unit,inline,eof=end_of_file)
!    print*,"Read line"
!    print*,inline
   if(Capitalize(inline(1:5)) /= "FIELD" ) then
      call MLSMessage(MLSMSG_Error,ModuleName,&
           "in subroutine l3ascii_read_field, File "//trim(filename)// &
           "on unit"//trim(unitstring)//" contains no more Fields")
       return
    endif
    
    if(end_of_file) then
       call MLSMessage(MLSMSG_Error,ModuleName,&
       "In subroutine l3ascii_read_field, End of File"//trim(filename)// &
       " on unit"//trim(unitstring))
       return
    endif

    read(unit=inline,fmt=*)linetype,field%quantityName, &
         field%description,field%units
    axesloop:do

        call ReadCompleteLineWithoutComments(unit,inline)
        !print*,inline
        read(unit=inline,fmt=*)linetype,axistype
        linetype=Capitalize(linetype)
        axistype=Capitalize(axistype)
        if(linetype(1:4) == "DATE") then ! This is always the last "axis"
           exit axesloop                 ! and is different from the others
        endif
        if (axistype(1:6) =="LINEAR") then
!           print*,"Doing linear axis"
           call make_linear_axis(inline,tmpaxis,tmpaxis_len)
!           print*,"Done linear axis"
        else if (axistype(1:3) =="LOG") then
           !print*,"Doing log axis"
           call make_log_axis(inline,tmpaxis,tmpaxis_len)
           !print*,"Done log axis"
        else if (axistype(1:8) =="EXPLICIT") then
!           print*,"Doing explicit axis"
           backspace(unit=unit)
           call read_explicit_axis(unit,tmpaxis,tmpaxis_len)
!           print*,"Done explicit axis"
        else
           call MLSMessage(MLSMSG_Error,ModuleName,&
                "in subroutine l3ascii_read_field,File"//trim(filename)//&
                " on unit"//trim(unitstring)//" contains coordinate"//&
                " of invalid type "//trim(axistype)//"for axis"//&
                trim(linetype))
           return
        endif

        ! I do not entirely grok what NJL intended verticalCoordinate to be. 
        if(linetype(1:8) == "PRESSURE" .or. linetype(1:8) == "ALTITUDE" &
             .or. linetype(1:3) == "GPH" .or. linetype(1:5) == "THETA") then
           field%noHeights=tmpaxis_len
           deallocate(field%heights)
           allocate(field%heights(1:tmpaxis_len))
           field%heights=tmpaxis

           if(linetype(1:8) == "PRESSURE") then
              field%verticalCoordinate=1
           else if (linetype(1:8) == "ALTITUDE") then
              field%verticalCoordinate=2
           else if (linetype(1:3) == "GPH") then
              field%verticalCoordinate=3
           else if (linetype(1:5) == "THETA") then
              field%verticalCoordinate=4
           endif
        else if(linetype(1:8) == "LATITUDE" .or. &
             linetype(1:8) == "EQUIVLAT") then
           field%noLats=tmpaxis_len
           deallocate(field%lats)
           allocate(field%lats(1:tmpaxis_len))
           field%lats=tmpaxis
           if (linetype(1:8) == "LATITUDE") then
              field%equivalentLatitude=.false.
           else 
              field%equivalentLatitude=.true.
           endif
        else if(linetype(1:9) == "LONGITUDE") then
           field%noLons=tmpaxis_len
           deallocate(field%lons)
           allocate(field%lons(1:tmpaxis_len))
           field%lons=tmpaxis
        else if(linetype(1:9) == "LST") then
           field%noLsts=tmpaxis_len
           deallocate(field%lsts)
           allocate(field%lsts(1:tmpaxis_len))
           field%lsts=tmpaxis
        else if(linetype(1:9) == "SZA") then
           field%noSzas=tmpaxis_len
           deallocate(field%szas)
           allocate(field%szas(1:tmpaxis_len))
           field%szas=tmpaxis
        endif
    enddo axesloop
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
    datesloop: do idate=1,maxNoDates
       !print*,"Datesloop: idate=",idate
       word_count=count_words(inline)
       if (word_count == 3) then
          read(unit=inline,fmt=*)linetype,axistype,sdstring
          sdstring=adjustl(sdstring)
          edstring=sdstring
       else if(word_count >= 4 ) then
          read(unit=inline,fmt=*)linetype,axistype,sdstring,edstring
          sdstring=adjustl(sdstring)
          edstring=adjustl(edstring)
       else
          call MLSMessage(MLSMSG_Error,ModuleName,&
               "in subroutine l3ascii_read_field: File"//trim(filename)//&
               "on unit"//trim(unitstring)//" contains a line beginning"//&
               trim(linetype)//"Date with too few words ")
       endif
       
        ! Date strings can begin with - indicating the year is 
        ! missing and that the file belongs to no year in particular.
        ! To convert dates to SDP toolkit (TAI) times (Seconds since start of 
        ! 1 Jan 1993) we need to stick on a dummy year
        if (sdstring(1:1) == "-") then
           sdstring=dummyyear//sdstring
        endif
        if (edstring(1:1) == "-") then
           edstring=dummyyear//edstring
        endif
        ! ccsds2tai returns days since 1 Jan 1993. 86400==no of secs per day
        dateStarts(idate)=86400*ccsds2tai(sdstring)
        dateEnds(idate)  =86400*ccsds2tai(edstring)
        call ReadCompleteLineWithoutComments(unit,inline)
        backspace(unit=unit)
        !print*,"About to read data: inline=",inline
        !print*,"tmpfield has size:",size(tmpfield),shape(tmpfield)
        
        read(unit=unit,fmt=*)tmpfield(:,:,:,:,:,idate)
        !print*,"Read data"
        end_of_file=.false.
        call ReadCompleteLineWithoutComments(unit,inline,eof=end_of_file)
        !print*,"Next date line:",inline,"EOF=",end_of_file
        if(end_of_file) then 
           ! No more dates and nothing else either
           exit datesloop
        endif
        read(unit=inline,fmt=*)linetype,axistype        
        linetype=Capitalize(linetype)
        axistype=Capitalize(axistype)
        if(end_of_file .or. &
             linetype(1:5) == "FIELD" .or. linetype(1:3)=="END") then 
           !Oops! There were no more dates, but there is another field
           backspace(unit=unit)
           exit datesloop                 
        endif
        if(linetype(1:4) /= "DATE") then ! There should be another date here
           call MLSMessage(MLSMSG_Error,ModuleName,&
                "in subroutine l3ascii_read_field: File"//trim(filename)//&
                "on unit"//trim(unitstring)//" contains a line beginning"//&
                trim(linetype)//"where I expected a line beginning Date ")
           return
        endif
        field%noDates=field%noDates+1

    enddo datesloop

    allocate(field%field(1:field%noHeights,1:field%noLats,&
         1:field%noLons,1:field%noLsts,1:field%noSzas,1:field%noDates))
    allocate(field%dateStarts(1:field%noDates),&
         field%dateEnds(1:field%noDates))
    field%dateStarts=dateStarts(1:field%noDates)
    field%dateEnds=dateEnds(1:field%noDates)
    field%field=tmpfield(:,:,:,:,:,1:field%noDates)
    deallocate(tmpfield,dateStarts,dateEnds)

  end subroutine l3ascii_read_field

  subroutine l3ascii_interp_field(field,outval,pressure,lat,lon,lst,sza,date)
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
    if(present(lat)) then
       inlat=lat
    else
       inlat=30.0
    endif

    if(present(pressure)) then
       inpressure=pressure
    else
       inpressure=10.0
    endif

    if(present(lon)) then
       inlon=lon
    else
       inlon=0.0
    endif

    if(present(sza)) then
       insza=sza
    else
       insza=30.0
    endif

    if(present(lst)) then
       inlst=lst
    else
       inlst=12.0
    endif

    if(present(date)) then
       indate=date
    else
       indate=15768000.0 ! This should be June (ish)
    endif

    ! OK, that's the optional arguments sorted. 
    ! Now, we have to find which two indices the chosen value lies between
    ! for each of the six values.
    call ilocate(field%lats,inlat,ilat1,ilat2)
    call ilocate(field%lons,inlon,ilon1,ilon2)
    call ilocate(field%szas,insza,isza1,isza2)
    call ilocate(field%lsts,inlst,ilst1,ilst2)
    ! We can't interpolate dates and prssures as they are. We construct
    ! a mean date and a log pressure
    allocate(tmpdate(1:field%noDates),tmpalt(1:field%noHeights))
    tmpdate=(field%dateStarts+field%dateEnds)/2.0
    tmpalt=-log10(field%heights)
    inalt=-log10(inpressure)
    call ilocate(tmpdate,indate,idate1,idate2)
    call ilocate(tmpalt,inalt,ialt1,ialt2)
    ! We now know that the desired point is inside the 6-D hypercuboid 
    ! bin of the array field%field(ialt1:ialt2, ilat1:ilat2, .......)
    ! Time to do the actual interpolation. Yuck. 
    ! There are up to 2^6=64 values
    ! to be taken into account. 

    ! Lets try allocating a cube of the right shape.
    if (ialt1==ialt2) then !test whether the altitude is between two of the 
       hcshape(1)=1        ! grid values (so we interpolate) or if there
    else                   ! is only one value (so we just return that value)
       hcshape(1)=2
    endif
    if (ilat1==ilat2) then
       hcshape(2)=1
    else
       hcshape(2)=2
    endif
    if (ilon1==ilon2) then
       hcshape(3)=1
    else
       hcshape(3)=2
    endif
    if (ilst1==ilst2) then
       hcshape(4)=1
    else
       hcshape(4)=2
    endif
    if (isza1==isza2) then
       hcshape(5)=1
    else
       hcshape(5)=2
    endif
    if (idate1==idate2) then
       hcshape(6)=1
    else
       hcshape(6)=2
    endif


    ! Allocate a 6-D hypercube that our point lies in. Some of the 
    ! dimensions may be "collapsed" 
    allocate(hcube(1:hcshape(1),1:hcshape(2),1:hcshape(3), &
         1:hcshape(4),1:hcshape(5),1:hcshape(6)))
    ! Copy data into hypercube
    hcube=field%field(ialt1:ialt2, ilat1:ilat2, ilon1:ilon2, &
         ilst1:ilst2,isza1:isza2, idate1:idate2)

    ! Now we interpolate along eahc of the axes where this is needed
    ! Reduce 6-d to 5-d
    if (hcshape(6) == 2 )then !Interpolate in date 
       do i0=1,hcshape(1)
          do i1=1,hcshape(2)
             do i2=1,hcshape(3)
                do i3=1,hcshape(4)
                   do i4=1,hcshape(5)
                      hcube(i0,i1,i2,i3,i4,1)=hcube(i0,i1,i2,i3,i4,1)+ &
                           (hcube(i0,i1,i2,i3,i4,2)- &
                           hcube(i0,i1,i2,i3,i4,1))* &
                           (indate-tmpdate(idate1))/ &
                           (tmpdate(idate2)-tmpdate(idate1))
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif
    ! Reduce  5-d to 4d
    if (hcshape(5) == 2 )then !Interpolate in local Solar zenith ang  SZA
       do i0=1,hcshape(1)
          do i1=1,hcshape(2)
             do i2=1,hcshape(3)
                do i3=1,hcshape(4)
                   hcube(i0,i1,i2,i3,1,1)=hcube(i0,i1,i2,i3,1,1)+ &
                        (hcube(i0,i1,i2,i3,2,1)- &
                        hcube(i0,i1,i2,i3,1,1))* &
                        (insza-field%szas(isza1))/ &
                        (field%szas(isza2)-field%szas(isza1))
                enddo
             enddo
          enddo
       enddo
    endif
    ! Reduce  4-d to 3d
    if (hcshape(4) == 2 )then !Interpolate in local solar time
       do i0=1,hcshape(1)
          do i1=1,hcshape(2)
             do i2=1,hcshape(3)
                hcube(i0,i1,i2,1,1,1)=hcube(i0,i1,i2,1,1,1)+ &
                     (hcube(i0,i1,i2,2,1,1)- &
                     hcube(i0,i1,i2,1,1,1))* &
                     (inlst-field%lsts(ilst1))/ &
                     (field%lsts(ilst2)-field%lsts(ilst1))
             enddo
          enddo
       enddo
    endif
    ! Reduce  3-d to 2d
    if (hcshape(3) == 2 )then !Interpolate in longitude
       do i0=1,hcshape(1)
          do i1=1,hcshape(2)
             hcube(i0,i1,1,1,1,1)=hcube(i0,i1,1,1,1,1)+ &
                  (hcube(i0,i1,2,1,1,1)- &
                  hcube(i0,i1,1,1,1,1))* &
                  (inlon-field%lons(ilon1))/ &
                  (field%lons(ilon2)-field%lons(ilon1))
          enddo
       enddo
    endif
    ! Reduce  2-d to 1d
    if (hcshape(2) == 2 )then !Interpolate in latitude
       do i0=1,hcshape(1)
          hcube(i0,1,1,1,1,1)=hcube(i0,1,1,1,1,1)+ &
               (hcube(i0,2,1,1,1,1)- &
               hcube(i0,1,1,1,1,1))* &
               (inlat-field%lats(ilat1))/ &
               (field%lats(ilat2)-field%lats(ilat1))
       enddo
    endif

    ! Reduce  1-d to 0d
    if (hcshape(1) == 2 )then !interpolate in altitude
       hcube(1,1,1,1,1,1)=hcube(1,1,1,1,1,1)+ &
            (hcube(2,1,1,1,1,1)- &
            hcube(1,1,1,1,1,1))* &
            (inalt-tmpalt(ialt1))/ &
            (tmpalt(ialt2)-tmpalt(ialt1))
    endif
    outval=hcube(1,1,1,1,1,1)
    deallocate(tmpdate,tmpalt,hcube)
  end subroutine l3ascii_interp_field

  subroutine ilocate(x,xval,ix1,ix2)
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
    if (n <=1) then ! x has only one element: no interpolating to be done
       ix1=1
       ix2=1
    else
       if (xval > x(n)) then ! off top . Use end value
          ix1=n
          ix2=n
       else if (xval < x(1)) then ! of bottom. Use end value
          ix1=1
          ix2=1
       else ! in range of x. Do binary search.
          j=n/2 
          dj = n/2 
          binsearch: do  
             if (dj > 1) then 
                dj=dj/2 
             else 
                dj=1 
             end if
             if(x(j) > xval .and. x(j+1) > xval) then 
                j=j-dj 
                cycle
             else  if (x(j) < xval .and. x(j+1) < xval) then 
                j=j+dj 
                cycle 
             else 
                exit
             endif
          enddo binsearch
          ix1=j
          ix2=j+1
       end if
    endif
  end subroutine ilocate

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
!       if (ioinfo /= 0) then
!          line="End of File Found"
!          exit rdloop
!       else
!          line=adjustl(line)                 ! remove blanks from start
!          if(line(1:1) /= ";" .and. line(1:1) /= " ") then !not a comment
!             exit rdloop                     ! so exit loop and return line
!          end if
!       endif
!       !        print*,line(1:80)
!    enddo rdloop

    !    print*,"Got non-comment line"
    !    print*,line(1:80)

!  end subroutine get_next_noncomment_line

  subroutine make_log_axis(inline,axis,axis_len)
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
    if (associated(axis)) then 
       deallocate(axis)
    endif

    !Count words in inline. 
    nwords=1
    do j=2,len(inline)
       if(inline(j:j) /= " " .and. inline(j-1:j-1) == " ") then
          nwords=nwords+1
       endif
    enddo
    nsections=(nwords-3)/2
!    print*,"Inline=",inline
!    print*,"Nwords=",nwords
    allocate(n_levs_in_sec(1:nsections),n_levs_per_dec(1:nsections),&
         axints(1:nsections*2))
    read(unit=inline,fmt=*)linetype,axistype,basepressure,axints
!    print*,"Linetype=",linetype," axistype=",axistype
!    print*,"basepressure=",basepressure," Axints=",axints
    n_levs_in_sec=axints(1:nsections*2-1:2)
    n_levs_per_dec=axints(2:nsections*2:2)

    axis_len=sum(n_levs_in_sec)

    allocate(axis(1:axis_len))
    axis(1)=-log10(basepressure)
    stind=0
    do j=1,nsections
       gridstep=1.0_r8/n_levs_per_dec(j)
       if (j == 1) then
          st=2
       else
          st=1
       endif
       do i=st,n_levs_in_sec(j)
          axis(stind+i)=axis(stind+i-1)+gridstep
       end do
       stind=stind+n_levs_in_sec(j)
    end do
    axis=10.0_r8**(-axis)
    deallocate(n_levs_in_sec, n_levs_per_dec,axints)

  end subroutine make_log_axis

  subroutine make_linear_axis(inline,axis,axis_len)
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
    if (associated(axis)) then 
       deallocate(axis)
    endif

    !Count words in inline. 

    nwords=1
    do j=2,len(inline)
       if(inline(j:j) /= " " .and. inline(j-1:j-1) == " ") then
          nwords=nwords+1
       endif
    enddo
    nsections=(nwords-3)/2

    allocate(n_levs_in_sec(1:nsections),gridstep(1:nsections),&
         axints(1:nsections*2))

    read(unit=inline,fmt=*)linetype,axistype,baseval,axints
    n_levs_in_sec=nint(axints(1:nsections*2-1:2))
    gridstep=axints(2:nsections*2:2)

    axis_len=sum(n_levs_in_sec)

    allocate(axis(1:axis_len))
    axis(1)=baseval
    stind=0
    do j=1,nsections
       if (j == 1) then
          st=2
       else
          st=1
       endif
       do i=st,n_levs_in_sec(j)
          axis(stind+i)=axis(stind+i-1)+gridstep(j)
       end do
       stind=stind+n_levs_in_sec(j)
    end do
    deallocate(axints,gridstep,n_levs_in_sec)

  end subroutine make_linear_axis

  subroutine read_explicit_axis(unit,axis,axis_len)
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
    if (associated(axis)) then 
       deallocate(axis)
    endif
    ! readitem needs to be initialised or there is a point where it 
    ! can be used before being set. The string should be more than ri_len
    ! spaces. 
    readitem="                                              "
    ! An explicit axis is supplied as a parenthesised list, spread 
    ! over several lines. This is a Royal PIA.
    ! Read chars till we get to the (
    do 
       read(unit=unit,fmt="(a)",advance="no")rdchar
       !        write(unit=*,fmt="(a)",advance="no")rdchar 
       if(rdchar=="(") then
          exit 
       endif
    enddo
    !    print*,"Got open paren"
    ! now read items and add them to the axis until we get to the )
    axis_len=0
    foundcb=.false.
    itemsloop:do
       i=1
       charsloop:do
          read(unit=unit,fmt="(a)",advance="no",iostat=iotest)rdchar
          !write(unit=*,fmt="(a)",advance="no")rdchar 
          if(rdchar==")") then
             !print*,"Found ) at end of explicit axis"
             foundcb=.true.
             exit charsloop
          endif
          if (rdchar==" ") then
             exit charsloop
          endif
          readitem(i:i)=rdchar
          i=i+1
       enddo charsloop
       if ((i<=1 .and. .not.foundcb) .or. (i==2 .and. readitem(1:1)==" ")) then
          cycle itemsloop
       endif
       if ( i <= 1 .and. foundcb) then
          exit itemsloop
       endif
       !print*,"Axis item is",readitem(1:i-1)
       axis_len=axis_len+1
       read(unit=readitem(1:i-1),fmt=*)tmpaxis(axis_len)
        !print*,"Element",axis_len," is ",tmpaxis(axis_len)
        if (foundcb) then
           exit itemsloop
        endif
    end do itemsloop
    allocate(axis(1:axis_len))
    axis=tmpaxis(1:axis_len)
  end subroutine read_explicit_axis

  function l3ascii_get_multiplier(field) result(multiplier)
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
    ucunits = Capitalize(field%units)
    ix = index(ucunits,"VMR ")
    if (ix > 0) then
       multiplier=1.0_r8
       return
    endif
    ix=index(ucunits,"PPM")
    if (ix > 0) then
       multiplier=1.0e6_r8
       return
    endif
    ix=index(ucunits,"PPB")
    if (ix > 0) then
       multiplier=1.0e9_r8
       return
    endif
    ix=index(ucunits,"PPT")
    if (ix > 0) then
       multiplier=1.0e12_r8
       return
    endif
    ix=index(ucunits,"K ")
    if (ix > 0) then
       multiplier=1.0_r8
       return
    endif
    call MLSMessage(MLSMSG_Warning,ModuleName,&
         "in function l3ascii_get_multiplier: Units "//&
         trim(field%units)//" for field "//trim(field%quantityName)//&
         "not known. Guessing multiplier=1.0")
    !print*,"in function l3ascii_get_multiplier: Units "//&
    !     trim(field%units)//" for field "//trim(field%quantityName)//&
    !     "not known. Guessing multiplier=1.0"
   multiplier=1.0_r8
  end function l3ascii_get_multiplier

end module gridded_data_module
