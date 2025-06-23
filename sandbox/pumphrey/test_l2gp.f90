program test_l2gp
  use L2GPData
  use HDFEOS, ONLY: SWOPEN,SWCLOSE,SWINQSWATH

  !What have you done with your mother?
  implicit none !If nun, write "none".

  !This program sucks in the contents of a EOS MLS L2GP file
  ! prints out the contents of the first 3 profiles

  character(len=80):: swathlist
  integer::nswath,strbufsize

  type(L2GPData_T)::swath
  character(len=80)::filename
  integer::swfid,numprofs,j
!  integer,parameter::np_plot=100
  
  print*,"Enter L2GP file name"
  read "(a80)",filename
  print*,"filename=",filename
  nswath=swinqswath(filename,swathlist,strbufsize)
  print*,"This file contains ",nswath,"swaths. A list of their names is:"
  print*,swathlist
  print*,"The string buffer size is ",strbufsize

  ! File has to be opened as a HDF-EOS file
  swfid = swopen(filename,1)
  print*,"File ID=",swfid

  call ReadL2GPData(swfid,trim(swathlist),swath,numprofs)

  do j=1,swath%nLevels
     print*,swath%pressures(j),swath%L2gpValue(1,j,1:2)
  enddo

  print*,"Lat,Long,sza"
  do j=1,10
     print*,  swath%latitude(j),swath%longitude(j),swath%solarZenith(j)
  enddo

  print*,"We Read ",numprofs,"profs"

end program test_l2gp
