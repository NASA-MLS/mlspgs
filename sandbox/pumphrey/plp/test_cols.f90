program test_cols

use plplot_module
! Experiment to find how color works in plplot -- it isn't documented yet
real(kind=PLFLT)::xmin,xmax,ymin,ymax,disp,pos,just,fpcol
integer,parameter::ncols=10
integer(kind=PLINT),dimension(ncols)::red,green,blue
integer(kind=PLINT)::i
character(len=20)::string
integer::j
xmin=0
xmax=1
ymin=0
ymax=1


call plinit()

! set colours in cmap 0 above 15 to be grey scale
! 

call plscmap0n(ncols)

do j=0,ncols-1
    i=255*real(j)/real(ncols)
    print*,"Setting color",j," to ",i
    red(j+1)=i
    green(j+1)=i
    blue(j+1)=i
!    call plscol0(j,i,i,i)! causes errors
end do

!call plscmap0(red,green,blue,ncols) ! Causes errors

call plenv(xmin,xmax,ymin,ymax,0,0)
just=0.0
pos=0.5
disp=-2

! colours in colour map 0
do i=1,ncols
    pos=real(i)/(ncols+1)
!    call plcol0(i)
    fpcol=real(i)/real(ncols)
    call plcol1(fpcol)
    write(unit=string,fmt=*) "Colour ",i
    print*,string
    call plmtex("lv",disp,pos,just,trim(adjustl(string)))
enddo


call plend()

end program test_cols
