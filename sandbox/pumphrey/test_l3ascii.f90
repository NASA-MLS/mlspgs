program test_l3ascii

  use gridded_data_module
  use kinds_module
  character(len=80)::filename
  integer::unit,counter
  type(GriddedData_T)::field
  logical :: doleak
  real(kind=r8)::outval!,pressure,lat
  ! variables for plotting section
  integer::ic,nperdec,nlevs,nlats,firstlat,lastlat,i,j

  print*,"Enter l3ascii file name "
  read "(a)" ,filename

  !open(file="deleteme",unit=1,status="replace",action="write")

  call l3ascii_open(filename,unit)
  field%reusing=0
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


  doleak=.true.
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
!         call l3ascii_read_field(unit,field)
   call l3ascii_interp_field(field,outval,lat=30.0_r8,pressure=45.0_r8)
         close(unit=unit)
         if (modulo(counter,100) == 0) then
            print*,"Done ",counter," open/read/close cycles"
         endif
     end do leak
  end if
 

  print*,"Shape of Field:",shape(field%field)
  print*,"Shape of lats :",shape(field%lats)


  

end program test_l3ascii
