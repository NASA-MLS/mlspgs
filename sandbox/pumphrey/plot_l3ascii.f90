program test_l3ascii

  use gridded_data_module
  use plplot_module
  use MLSCommon
  character(len=80)::filename
  integer::unit,counter
  type(GriddedData_T)::field
  logical :: doleak
  real(kind=r8)::outval!,pressure,lat
  ! variables for plotting section
  integer::ic,nperdec,nlevs,nlats,firstlat,lastlat,i,j
  real(kind=PLFLT),allocatable,dimension(:,:)::z
  real(kind=PLFLT),allocatable,dimension(:)::x,y
  real(kind=PLFLT),parameter,dimension(25)::clev=0.2*(/ (ic,ic=0,24) /)+2.0
!  real(kind=PLFLT),parameter,dimension(25)::clev=5*(/ (ic,ic=0,24) /)+180.0

  print*,"Enter l3ascii file name  "
  read "(a)" ,filename

  !open(file="deleteme",unit=1,status="replace",action="write")

  call l3ascii_open(filename,unit)
  print*,"Using unit",unit  
  print*,"Reading file"
  call l3ascii_read_field(unit,field)
  print*,"Quantity: ",field%quantityName
  print*,"Description: ",field%description
  print*,"Units: ",field%units
  print*,"No of Heights: ",field%noHeights
  print*,"Vert coord type: ",field%verticalCoordinate
  print*,"Heights: ",field%heights
  print*,"No of Lats",field%noLats
  print*,"Lats: ",field%lats
  print*,"No of Lons",field%noLons
  print*,"Lons: ",field%lons
  print*,"No of Lsts",field%noLsts
  print*,"Lsts: ",field%lsts
  print*,"No of Szas",field%noSzas
  print*,"Szas: ",field%szas
  print*,"No. of dates",field%noDates
  print*,"start Dates:",field%dateStarts/(24*60*60)
  print*,"end Dates:",field%dateEnds/(24*60*60)
  close(unit=unit)
  
  ! test intrpolate routine
!  lat=30.0
!  pressure=1.0
!  inttest:do 
!      call l3ascii_interp_field(field,outval,lat=lat,pressure=pressure)
!      print*,"Lat=",lat,"pr=",pressure," value=",outval
!      if (lat > 99.0) then 
!         exit inttest
!      endif
!      print*,"Enter another lat and pressure"
!      read*,lat,pressure
!  enddo inttest

  doleak=.false.
  if (doleak) then 
     print*," Applying leakage test. Type Ctrl-C to stop"
     print*," While this is running, run top and check the size"
     print*," of this program's process. If the process keeps getting bigger"
     print*," then there is a memory leak in this program "
     print*," (or maybe your Fortran 90 compiler.......)"
     counter=0
     leak:do counter=1,800
         call l3ascii_open(filename,unit)
         call l3ascii_read_field(unit,field)
         call l3ascii_read_field(unit,field)
!   call l3ascii_interp_field(field,outval,lat=30.0_r8,pressure=45.0_r8)
         close(unit=unit)
         if (modulo(counter,100) == 0) then
            print*,"Done ",counter," open/read/close cycles"
         endif
     end do leak
  end if
 
  ! graphical test section. Uses plplot. Replace with calls to your 
  ! preferred graphics package. 

  ! Test interpolator

  call plsdev("xwin")
  call plinit()
  allocate(z(field%noLats,field%noHeights),x(field%noLats),y(field%noHeights))
  z=transpose(field%field(:,:,1,1,1,1))
  x=field%lats
  y=-log10(field%heights)
  call plfcp(z,x,y,clevels=clev)
  deallocate(z,x,y)

  call plend() ! An obscure X error happens if you don't do this
  call plsdev("xwin")
  call plinit()! but only if you are in true colour mode. These two lines
  !              are unnecessary in 8-bit colour mode

  allocate(z(field%noDates,field%noHeights),x(field%noDates),&
       y(field%noHeights))
  z=transpose(field%field(:,1,1,1,1,:))
  x=field%dateStarts
  y=-log10(field%heights)
  call plfcp(z,x,y,clevels=clev)
  deallocate(z,x,y)
  call plend() ! An obscure X error happens if you don't do this
  call plsdev("xwin")
  call plinit()! but only if you are in true colour mode. These two lines
  !              are unnecessary in 8-bit colour mode(or 16-bit colour mode)
  
  allocate(z(field%noDates,field%noLats),x(field%noDates),&
       y(field%noLats))
  z=transpose(field%field(1,:,1,1,1,:))
  x=field%dateStarts
  y=field%lats
  call plfcp(z,x,y,clevels=clev)
  deallocate(z,x,y)

  
  nlats=85
  firstlat=-88.0
  lastlat=88.0
  nlevs=49
  nperdec=12
  call plend()
  call plsdev("xwin")
  call plinit() ! initiates PLPLOT
  
  allocate(z(1:nlats,1:nlevs),x(1:nlats),y(1:nlevs))
  x=firstlat+(/ (ic,ic=0,nlats-1) /)*(lastlat-firstlat)/real(nlats-1,kind=r8)
  y=-3.0 + (/ (ic/real(nperdec),ic=0,nlevs-1) /)
  y=10.0**(-y)
!  print*,"Lats are",x
  do i=1,nlats
      do j=1,nlevs
          call l3ascii_interp_field(field,outval,lat=x(i),pressure=y(j), &
               date=10.0_r8)
          z(i,j)=outval
      enddo
  enddo
  y=-log10(y)
  call plfcp(z,x,y,clevels=clev) !Makes contour plot
  deallocate(z,x,y)
  call plend()!ends PLPLOT


  print*,"Shape of Field:",shape(field%field)
  print*,"Shape of lats :",shape(field%lats)


  

end program test_l3ascii
