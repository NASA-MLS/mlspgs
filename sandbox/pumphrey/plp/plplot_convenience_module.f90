module plplot_convenience_module

  ! This module provides a few convenience calls to make plplot easier to 
  ! call from Fortran 90 or F. So far, plfcp is the only one. 
  ! I must split the label bar off as a separate subroutine.
  use plplot_params_module
  use plplot_interface_module

  private

  public:: plfcp

contains
  subroutine plfcp(z,x,y,nlevels,clevels,xtitle,ytitle,title)
    ! This subroutine gets you a filled contour plot: with a label bar !
    ! all you need to do is
    ! call plinit(), call plfcp(z,x,y), call plend. I am adding customisation
    ! as OPTIONAL args so the plot can be styled, but without confusing the 
    ! simple use. !
    ! --------------Arguments-------------!
    real(kind=PLFLT),intent(in),dimension(:,:)::z
    real(kind=PLFLT),intent(in),dimension(:)::x,y
    integer(kind=PLINT),optional,intent(in)::nlevels
    real(kind=PLFLT),optional,intent(in),dimension(:)::clevels
    character(len=*),optional,intent(in)::xtitle,ytitle,title
    !---- Local vars ----!
    character(len=40)::my_xtitle,my_ytitle,my_title
    integer(kind=PLINT)::nx,ny
    integer(kind=PLINT)::sh_cmap,min_color,min_width,max_color,max_width
    integer(kind=PLINT)::sh_width,nlevel
    real(kind=PLFLT),parameter::zero=0
    real(kind=PLFLT)::xmin,xmax,ymin,ymax,sh_color!,zmin,zmax
    real(kind=PLFLT)::zmi,zma,disp,pos,just
    real(kind=PLFLT),allocatable,dimension(:)::clevel
    real(kind=PLFLT)::vpx0,vpx1,vpy0,vpy1
    real(kind=PLFLT),dimension(20)::polyx,polyy

    ! Color map
    ! Supposed to be the nice set I use in IDL.
    integer(kind=PLINT),parameter::ncsp=9
    real(kind=PLFLT),dimension(ncsp),parameter::cpos= &
         (/0.0,  1.0,  2.3,  3.0,  3.6,  5.0,  6.0,  6.01,  7.0/) / 7.0
    real(kind=PLFLT),dimension(ncsp),parameter::red = &
         (/0.0,  0.0,  0.0,  0.0,  1.0,  1.0,  1.0,  1.0, 1.0/)
    real(kind=PLFLT),dimension(ncsp),parameter::green=&
         (/0.0,  0.0,  1.0,  1.0,  1.0,  0.0,  0.0,  0.0, 1.0/)
    real(kind=PLFLT),dimension(ncsp),parameter::blue =&
         (/0.0,  1.0,  1.0,  0.0,  0.0,  0.0,  1.0,  1.0, 1.0/)
    integer(kind=PLINT),dimension(ncsp),parameter::rev=&
         (/1,  0,  0,  0,  0,  1,  1,  1,  1/)
!    integer(kind=PLINT),parameter::ncsp=7
!    real(kind=PLFLT),dimension(ncsp),parameter::cpos= &
!         (/0.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0/) / 6.0
!    real(kind=PLFLT),dimension(ncsp),parameter::blue = &
!         (/0.0,  0.0,  0.25,  0.5,  0.75,  1.0,  1.0/)
!    real(kind=PLFLT),dimension(ncsp),parameter::red=&
!         (/0.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0 /)
!    real(kind=PLFLT),dimension(ncsp),parameter::green =& 
!         (/1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  0.0/)
!    integer(kind=PLINT),dimension(ncsp),parameter::rev=&
!         (/0,  0,  0,  0,  0,  0,  0/)

    

    character(len=20)::string
    integer:: ic,i,j,skip
    !--------Executable-------------!
    nx=size(x)
    ny=size(y)! must add check that z is the right size
    if(present(clevels) )then 
       nlevel=size(clevels)
       allocate(clevel(nlevel))
       clevel=clevels
    else
       if(present(nlevels)) then
          nlevel=nlevels
       else
          nlevel=21
       endif
       allocate(clevel(nlevel))
       zmi=minval(z)
       zma=maxval(z)
       clevel=zmi + (/ (ic*(zma-zmi)/(nlevel-1),ic=0,nlevel-1) /)
    endif

    !  text 
    if(present(xtitle)) then
       my_xtitle=xtitle
    else
       my_xtitle="X"
    endif
    if(present(ytitle)) then
       my_ytitle=ytitle
    else
       my_ytitle="X"
    endif
    if(present(title)) then
       my_title=title
    else
       my_title="Contour"
    endif

!    call plscmap0n(5)
    xmin=x(1)
    xmax=x(nx)
    ymin=y(1)
    ymax=y(ny)
    vpx0=0.1
    vpx1=0.9
    vpy0=0.3
    vpy1=0.9
!    print*,"Setting color to 1"
    call plcol(1)
!    call plenv(xmin,xmax,ymin,ymax,0,0)
!    print*,"Advancing"
    call pladv(0)
!    print*,"Setting viewport"
    call plvpor(vpx0,vpx1,vpy0,vpy1)
!    print*,"Setting window"
    call plwind(xmin,xmax,ymin,ymax)
!    print*,"Setting color to 2"
    call plcol(2)
!    print*,"Setting x coord"
    call pllab(trim(my_xtitle),trim(my_ytitle),trim(my_title))
!    print*,"Setting color map 1 for contour shading"
    ! Set colour map of whatever putridity we like
    call plscmap1l(1,ncsp,cpos,red,green,blue,rev)
    !Do filled contours

    sh_cmap = 1   ! I wish I knew what these did. Some of this is 
    min_color = 0 ! cargo-cult programming.
    min_width = 0
    max_color = 0
    max_width = 0
    sh_width=1
!    print*,"Shading contours"
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
!    print*,"Setting color to 4"
    call plcol(4)
!    print*,"Drawing contours"
    call plcon1(z,nx,ny,1,nx,1,ny,clevel,nlevel,x,y)
    ! Draw axes. Wise to do this last for filled plots or the filling
    ! ends up on top of the ticks.
!    print*,"Drawing box"
    call plbox("bncst",zero,0,"bncst",zero,0) 

    ! ---------------- label bar code starts here --------------------!  
    vpx0=0.1
    vpx1=0.9
    vpy0=0.15
    vpy1=0.2
    xmin=0.0
    xmax=1.0
    ymin=0.0
    ymax=1.0
!    print*,"Setting viewport for label bar"
    call plvpor(vpx0,vpx1,vpy0,vpy1)
!    print*,"Setting window for label bar"
    call plwind(xmin,xmax,ymin,ymax)
    polyy(1:5)=(/0.0,1.0,1.0,0.0,0.0/)
    ! loop to draw the coloured boxes in the label bar
!    print*,"Colouring label bar"
    boxloop:do j=1,nlevel-1
       polyx(1:5)=(/j-1,j-1,j,j,j-1/)*1.0 / (nlevel-1)
       sh_color=j
       sh_color=sh_color/(nlevel-1)
       call plcol1(sh_color)
       call plfill(5,polyx(1:5),polyy(1:5))
       call plcol0(4)
    enddo boxloop
    
    ! loop to put lines between the boxes and add text labels
    ! If we ever make the contours different colours, we could do 
    ! that here for the label bar
!    print*,"Setting color to 4"
    call plcol0(4)
    just=0.5
    disp=1.0
    skip=1+nlevel/8
!    print*,"Drawing lines in label bar"
    lineloop:do j=1,nlevel
        polyx(3:4)=(/j-1,j-1/)*1.0 / (nlevel-1)
        pos=real(j-1)/(nlevel-1)
        write(unit=string,fmt="(f6.2)")clevel(j)
!        print*,"Drawing line no.",j
        call plline(2,polyx(3:4),polyy(3:4))
        if (modulo(j,skip)==0) then
!           print*,"Writing label",j," which is",string
           call plmtex("b",disp,pos,just,trim(adjustl(string)))
           
        endif
        
    enddo lineloop
    ! draw box round whole bar 
!    print*,"Should be Setting colour to 4"
    call plcol0(4)
    polyx(1:5)=(/0.0,0.0,1.0,1.0,0.0/)
!    print*,"drawing box round bar"
    call plline(5,polyx(1:5),polyy(1:5))
    ! --------------------- Label Bar finished -------------------!
    deallocate (clevel)
!    print*,"Leaving plfcp"
  end subroutine plfcp

end module plplot_convenience_module
