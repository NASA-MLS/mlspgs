program test_2d

  ! This program is a first attempt to use the plplot library 
  ! for contour plots
  use plplot_module
  !-- Variables--!
  integer,parameter::nx=101,ny=50,nlevel=21
  real,parameter::pi=3.1415926536
  real(kind=PLFLT),dimension(nx)::x
  real(kind=PLFLT),dimension(ny)::y
  real(kind=PLFLT),dimension(nx,ny)::z
  real(kind=PLFLT),dimension(nlevel)::clevel  
  real(kind=PLFLT)::xmin,xmax,ymin,ymax,shear
  real(kind=PLFLT)::xmult,ymult,sh_color
  real(kind=PLFLT),dimension(6)::tr
  !  character(len=1)::foo
  integer::ic,i
  integer(kind=PLINT)::sh_cmap,min_color,min_width,max_color,max_width,sh_width

  ! Nasty. Very nasty !
  !  common /plplot/ tr

  !------Executable code--------------------!
  clevel=-1 + (/ (ic*2.0/(nlevel-1),ic=0,nlevel-1) /)

  x=(/ (4*pi*ic/real(nx-1),ic=0,nx-1) /)-2*pi
  y=(/ (2*pi*ic/real(ny-1),ic=0,ny-1) /)-pi

  ! Hack to make x axis not evenly spaced
  x=sin(x/4)*2*pi
  do i=1,nx
      z(i,:)=sin(y)*cos(x(i))
  enddo
  call plscmap0n(5)
  call plinit()

  xmin=1
  xmax=nx
  ymin=1
  ymax=ny

  call plfcp(z,x,y)

  call plenv(xmin,xmax,ymin,ymax,0,0)
  call plcol(2)
  call pllab("X Coord","Y Coord","Boring")
  call plcol(4)
  call plcon0(z,nx,ny,1,nx,1,ny,clevel,nlevel)

  xmin=x(1)
  xmax=x(nx)
  ymin=y(1)
  ymax=y(ny)
  shear=0.0
  xmult=(x(nx)-x(1))/(nx-1)
  ymult=(y(ny)-y(1))/(ny-1)
  tr=(/ xmult,shear,xmin,ymult,shear,ymin /)
  call pladv(0)
  call plcol(1)
  call plenv(xmin,xmax,ymin,ymax,0,0)
  call plcol(2)
  call pllab("X Coord","Y Coord","Boring")

  !Do filled contours

  sh_cmap = 1
  min_color = 1
  min_width = 0
  max_color = 0
  max_width = 0
  sh_width=1
  shadeloop:do i=1,nlevel-1
      sh_color=i
      sh_color=sh_color/(nlevel-1)
      ! call plpsty(modulo(i,8)) ! gets different shade patterns 

      !      works if x and y evenly spaced
      !      call plshade0(z,nx,ny," ", xmin,xmax,ymin,ymax,&
      !           clevel(i),clevel(i+1),sh_cmap,sh_color,sh_width,&
      !           min_color,min_width,max_color,max_width)

      call plshade1(z,nx,ny," ", xmin,xmax,ymin,ymax,&
           clevel(i),clevel(i+1),sh_cmap,sh_color,sh_width,&
           min_color,min_width,max_color,max_width,x,y)

  enddo shadeloop
  call plcol(4)
  print*,tr
  call plcon1(z,nx,ny,1,nx,1,ny,clevel,nlevel,x,y)
  call plend()

end program test_2d
