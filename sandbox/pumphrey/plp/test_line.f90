program test_line

  ! This program is a first attempt to use the plplot library
  use plplot_module

  integer,parameter::n=11
  real(kind=PLFLT),dimension(n)::x,y
  integer::ic

  x=(/ (ic/real(n-1),ic=0,n-1) /)
  y=x*x

  call plinit()
  call plenv(x(1),x(n),y(1),y(n),0,1)
  call plline(n,x,y)
  call plend()

end program test_line
