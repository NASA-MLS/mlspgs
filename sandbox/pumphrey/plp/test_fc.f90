program test_fc

  ! This program tests  a Fortran 90 convienience layer between the programmer
  ! and the plplot library
  ! for contour plots
  use plplot_module
  !-- Variables--!
  integer,parameter::nx=101,ny=50,nlevel=7
  real,parameter::pi=3.1415926536
  real(kind=PLFLT),dimension(nx)::x
  real(kind=PLFLT),dimension(ny)::y
  real(kind=PLFLT),dimension(nx,ny)::z
  real(kind=PLFLT),dimension(nlevel)::clevel  
  integer::i,ic
  !------Executable code--------------------!
  clevel=-1.1 + (/ (ic*2.2/(nlevel-1),ic=0,nlevel-1) /)

  x=(/ (4*pi*ic/real(nx-1),ic=0,nx-1) /)-2*pi
  y=(/ (2*pi*ic/real(ny-1),ic=0,ny-1) /)-pi
  ! Hack to make x axis not evenly spaced
  x=sin(x/4)*2*pi

  ! make up 2-D array to contour
  do i=1,nx
      z(i,:)=sin(y)*cos(x(i))
  enddo

!  call plscmap0n(5)

  call plinit()
  call plfcp(z,x,y,clevels=clevel)
  print*,clevel
  call plend()

end program test_fc
