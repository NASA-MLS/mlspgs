program test_l3ascii

  use l3ascii
  use plplot_module
  use MLSCommon
  character(len=80)::filename
  integer::unit
  logical::endoffile
  type(GriddedData_T)::field

  real(kind=r8)::outval,mult!,pressure,lat
  ! variables for plotting section
  integer::ic,nperdec,nlevs,nlats,firstlat,lastlat,i,j
  real(kind=PLFLT),allocatable,dimension(:,:)::z
  real(kind=PLFLT),allocatable,dimension(:)::x,y
  real(kind=PLFLT),parameter,dimension(25)::clev=0.2*(/ (ic,ic=0,24) /)+2.0
  !  real(kind=PLFLT),parameter,dimension(25)::clev=5*(/ (ic,ic=0,24) /)+180.0
  character(len=20)::multstring
  print*,"Enter l3ascii file name  "
  read "(a)" ,filename

  !open(file="deleteme",unit=1,status="replace",action="write")

  call l3ascii_open(filename,unit)
  print*,"Using unit",unit  
  print*,"Reading file"
  fieldsloop:do
     call l3ascii_read_field(unit,field,endoffile)
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

     print*,"Shape of Field:",shape(field%field)
     print*,"Shape of lats :",shape(field%lats)


     ! graphical test section. Uses plplot. Replace with calls to your 
     ! preferred graphics package. 

     ! Test interpolator

     nlats=85
     firstlat=-88.0
     lastlat=88.0
     nlevs=49+24
     nperdec=12
     call plsdev("xwin")
     call plinit() ! initiates PLPLOT

     allocate(z(1:nlats,1:nlevs),x(1:nlats),y(1:nlevs))
     x=firstlat+(/ (ic,ic=0,nlats-1) /)*(lastlat-firstlat)/ &
          real(nlats-1,kind=r8)
     y=-3.0 + (/ (ic/real(nperdec),ic=0,nlevs-1) /)
     y=10.0**(-y)
     !  print*,"Lats are",x
     do i=1,nlats
        do j=1,nlevs
           call l3ascii_interp_field(field,outval,lat=x(i),pressure=y(j), &
                date=0.0_r8)
           z(i,j)=outval
        enddo
     enddo
     z=z/l3ascii_get_multiplier(field)
     y=-log10(y)
     if(maxval(z)  >  1.0e-2) then
        mult=1
        multstring=" "
     else if(maxval(z)  >  1.0e-5) then
        mult=1.0e3
        multstring=" x 1.e3 "
     else if(maxval(z)  >  1.0e-8) then
        mult=1.0e6
        multstring=" x 1.e6 "
     else if(maxval(z)  >  1.0e-11) then
        mult=1.0e9
        multstring=" x 1.e9 "
     else if(maxval(z)  >  1.0e-14) then
        mult=1.0e12
        multstring=" x 1.e12 "
     endif
     z=z*mult
     !Makes contour plot
     call plfcp(z,x,y,xtitle="Latitude",ytitle="Log pressure",&
          title=trim(field%quantityName)//" / "//trim(field%units) & 
          // trim(multstring)) 
     
     deallocate(z,x,y)
     call plend()!ends PLPLOT



  end do fieldsloop
  
  close(unit=unit)

end program test_l3ascii
