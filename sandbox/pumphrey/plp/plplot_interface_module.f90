module plplot_interface_module

  ! Module to provide access to plplot from F/fortran 90. 

  ! For PSR VAST/F90 link opts are: -ltcl8.0 -lMatrix -lplplotdX -lplplotdtk
  ! for Fujitsu F90, link opts are: "-L/usr/X11R6/lib \
  ! -L/usr/lib/gcc-lib/i486-linux/egcs-2.91.60/  -ltcl8.0 -lMatrix -lplplotdX \
  ! -lplplotdtk -ltcl8.0  -lg2c "


  public:: plinit,plenv,plline,plend,pllab,plpoin,plsym,pljoin,plptex,plfill
  public:: plcon0,plcon1,pladv,plcol,plscmap0n,plpsty,plshade0,plbox,plshade1
  public:: plvpor,plwind,plscmap0,plscmap1!,plscmap1n,
  public:: plmtex,plsdev


  public:: plscol0,plscmap1l,plcol0,plcol1


  interface 

     subroutine plinit()
     end subroutine plinit

     subroutine plwind(xmin, xmax, ymin, ymax)
       use plplot_params_module
       real(kind=PLFLT),intent(in)::xmin, xmax, ymin, ymax
     end subroutine plwind

     subroutine plmtex(side, disp, pos, just, text) 
       use plplot_params_module
       character(len=*),intent(in)::side,text
       real(kind=PLFLT),intent(in)::disp,pos,just
     end subroutine plmtex

     subroutine plvpor(xmin, xmax, ymin, ymax)
       use plplot_params_module
       real(kind=PLFLT),intent(in)::xmin, xmax, ymin, ymax
     end subroutine plvpor

     subroutine plbox(xopt, xtick, nxsub, yopt, ytick, nysub)
       use plplot_params_module
       character(len=*),intent(in)::xopt,yopt
       real(kind=PLFLT),intent(in)::xtick,ytick
       integer(kind=PLINT),intent(in)::nxsub,nysub
     end subroutine plbox

     subroutine plscol0(icol0,r,g,b) 
       ! Set color icol0 from color map 0 by 8 bit RGB value
       ! Does not result in any additional cells to be allocated.
       use plplot_params_module
       integer(kind=PLINT),intent(in)::icol0,r,g,b
     end subroutine plscol0

     subroutine plscmap0(r,g,b,ncol0)
       ! Set color map 0 colors by 8 bit RGB values.  
       ! This sets the entire color
       ! map -- only as many colors as specified will be allocated.
       use plplot_params_module
       integer(kind=PLINT),intent(in),dimension(*)::r,g,b
       integer(kind=PLINT),intent(in)::ncol0
     end subroutine plscmap0

     subroutine plscmap1(r,g,b,ncol1)
       ! Set color map 1 colors by 8 bit RGB values
       ! This also sets the number of colors.
       use plplot_params_module
       integer(kind=PLINT),intent(in),dimension(*)::r,g,b
       integer(kind=PLINT),intent(in)::ncol1
     end subroutine plscmap1

!     subroutine plscmap1n(ncol1)
       ! Set number of colors in cmap 1, (re-)allocate cmap 1, and set default
       ! values if this is the first allocation.
       !
       ! Note that the driver is allowed to disregard this number.
       ! In particular, most use fewer than we use internally.
!       use plplot_params_module
!       integer(kind=PLINT),intent(in)::ncol1       
!     end subroutine plscmap1n

     subroutine plscmap1l(itype,npts,pos,coord1,coord2,coord3,rev)

       use plplot_params_module
       integer(kind=PLINT),intent(in)::itype,npts
       real(kind=PLFLT),intent(in),dimension(*)::pos,coord1,coord2,coord3
       integer(kind=PLINT),intent(in),dimension(*)::rev


       !-----------------------------------------------------------------
       ! plscmap1l()
       !
       ! Set color map 1 colors using a piece-wise linear relationship between
       ! position in the color map (from 0 to 1) and position in HLS 
       !or RGB colour space.  May be called at any time.
       !
       ! The idea here is to specify a number of control points that
       ! specify the
       ! mapping between HLS (or RGB or CMY) and palette 1 value. Between these
       ! points, linear interpolation is used. By mapping position in the color
       ! map to function value, this gives a smooth variation of color with
       ! intensity.  Any number of control points may be specified, located at
       ! arbitrary positions (intensities), although typically 2-4 are enough.
       ! Another way of stating this is that we are traversing a given no. of
       ! lines through HLS (or RGB) space as we move through cmap1 entries. The
       ! control points at the minimum and maximum intensity (0 and 1) must
       ! always be specified.  By adding more control points you can get more
       ! variation.  One good technique for plotting functions that vary about
       ! some expected average is to use an additional 2 control points in the
       ! center (intensity ~= 0.5) that are the same color as the background
       ! (typically white for paper output, black for crt), and same hue as the
       ! boundary control points.  This allows the highs and lows to be very
       ! easily distinguished.
       !
       ! Each control point must specify the position in cmap 1 as well as 
       ! three
       ! coordinates in HLS or RGB space.  The first point MUST correspond to
       ! position = 0, and the last to position = 1.  
       !
       ! The hue is interpolated around the "front" of the color wheel
       ! (red<->green<->blue<->red) unless the "rev" flag is set, in which case
       ! interpolation proceeds around the back (reverse) side.  Specifying
       ! rev=NULL is equivalent to setting rev[]=0 for every control point.
       !
       ! Bounds on RGB coordinates:
       !	R,G,B		[0, 1]		magnitude
       !
       ! Bounds on HLS coordinates:
       !	hue		[0, 360]	degrees
       !	lightness	[0, 1]		magnitude
       !	saturation	[0, 1]		magnitude
       !
       ! The inputs are:
       !	itype		0: HLS, 1: RGB
       !	npts		number of control points
       !	pos[]		position for each control point
       !	coord1[]	first coordinate for each control point
       !	coord2[]	second coordinate for each control point
       !	coord3[]	third coordinate for each control point 
       !	rev[]		reverse flag for each control point
       !-----------------------------------------------------------
     end subroutine plscmap1l

     subroutine plpsty(n)
       use plplot_params_module
       integer(kind=PLINT),intent(in)::n
     end subroutine plpsty

     subroutine plscmap0n(sub)
       ! Not documented yet. Appears to make "sub" colors (0 to sub-1)
       ! available in color map 0. There is also a color map 1. 
       ! I can't find how to switch to it. The default is sub=16. Colour 15
       ! is white and colour 0 is black.
       use plplot_params_module
       integer(kind=PLINT),intent(in)::sub
     end subroutine plscmap0n

     subroutine pladv(sub)
       use plplot_params_module
       integer(kind=PLINT),intent(in)::sub
     end subroutine pladv

     subroutine plsdev(sub)
       use plplot_params_module
       character(len=*),intent(in)::sub
     end subroutine plsdev

     subroutine plcol(color)
       use plplot_params_module
       integer(kind=PLINT),intent(in)::color
     end subroutine plcol

     subroutine plcol0(color)
       use plplot_params_module
       integer(kind=PLINT),intent(in)::color
     end subroutine plcol0

     subroutine plcol1(color)
       !Switches to a colour from cmap 1, chosen as a value between
       ! 0.0 and 1.0
       use plplot_params_module
       real(kind=PLFLT),intent(in)::color
     end subroutine plcol1

     subroutine plenv(xmin, xmax, ymin, ymax, just, axis)
       use plplot_params_module
       real(kind=PLFLT),intent(in)::xmin,xmax,ymin,ymax
       integer(kind=PLINT),intent(in)::just,axis
     end subroutine plenv

     subroutine plline(n,x,y)
       use plplot_params_module
       real(kind=PLFLT),intent(in),dimension(*)::x,y
       integer(kind=PLINT),intent(in)::n
     end subroutine plline

     subroutine pllab(xlbl,ylbl,toplbl)
       character(len=*),intent(in)::xlbl,ylbl,toplbl
     end subroutine pllab

     subroutine plpoin(n, x, y, code)
       use plplot_params_module
       integer(kind=PLINT),intent(in)::n,code
       real(kind=PLFLT),intent(in),dimension(*)::x,y
     end subroutine plpoin

     subroutine plsym(n, x, y, code)
       use plplot_params_module
       integer(kind=PLINT),intent(in)::n,code
       real(kind=PLFLT),intent(in),dimension(*)::x,y
     end subroutine plsym

     subroutine pljoin(x1, y1, x2, y2)
       use plplot_params_module
       real(kind=PLFLT),intent(in)::x1,y1,x2,y2
     end subroutine pljoin

     subroutine plptex(x, y, dx, dy, just, text)
       use plplot_params_module
       real(kind=PLFLT),intent(in)::x,y,dx,dy,just
       character(len=*),intent(in)::text
     end subroutine plptex

     subroutine plfill(n,x,y)
       use plplot_params_module
       real(kind=PLFLT),intent(in),dimension(*)::x,y
       integer(kind=PLINT),intent(in)::n
     end subroutine plfill

     subroutine plcon0(z, nx, ny, kx, lx, ky, ly, clevel, nlevel)
       use plplot_params_module
       integer(kind=PLINT),intent(in)::nx,ny,kx,lx,ky,ly,nlevel
       real(kind=PLFLT),intent(in),dimension(*)::z
       real(kind=PLFLT),intent(in),dimension(*)::clevel
     end subroutine plcon0

     subroutine plshade0(z, nx, ny, string, xmin,xmax,ymin,ymax, &
          shade_min,shade_max,sh_cmap,sh_color,sh_width,&
          min_color,min_width,max_color,max_width)
       ! Not documented. Appears to fill the space between two contours.
       ! Need to call it once for each pair of contour levels you want filled 
       ! between. requires evenly spaced grid. For uneven grid use plshade1.
       ! z->data, nx,ny->size of z,xmin,xmax,ymin,ymax-> x axis values
       ! shade_min, shade_max: upper and lower boundaries of shaded area
       ! sh_cmap gives a blue-red color map, anything else seems to
       ! give an error. sh_color is the color you want. 0.0 is one end of
       ! the color table, 1.0 is the other
       use plplot_params_module
       integer(kind=PLINT),intent(in)::nx,ny
       character(len=*),intent(in)::string
       real(kind=PLFLT),intent(in):: xmin,xmax,ymin,ymax,shade_min,shade_max
       real(kind=PLFLT),intent(in):: sh_color
       integer(kind=PLINT),intent(in):: sh_cmap,sh_width,min_color
       integer(kind=PLINT),intent(in):: min_width,max_color,max_width
       real(kind=PLFLT),intent(in),dimension(*)::z
     end subroutine plshade0

     subroutine plshade1(z, nx, ny, string, xmin,xmax,ymin,ymax, &
          shade_min,shade_max,sh_cmap,sh_color,sh_width,&
          min_color,min_width,max_color,max_width,xg,yg)
       ! Not documented. Appears to fill the space between two contours.
       ! Need to call it once for each pair of contour levels you want filled 
       ! between. 
       ! z->data, nx,ny->size of z,xmin,xmax,ymin,ymax-> x axis values
       ! shade_min, shade_max: upper and lower boundaries of shaded area
       ! sh_cmap gives a blue-red color map, anything else seems to
       ! give an error. sh_color is the color you want. 0.0 is one end of
       ! the color table, 1.0 is the other
       ! xg -> 1-d array of x values, yg = 1-d array of y values
       use plplot_params_module
       integer(kind=PLINT),intent(in)::nx,ny
       character(len=*),intent(in)::string
       real(kind=PLFLT),intent(in):: xmin,xmax,ymin,ymax,shade_min,shade_max
       real(kind=PLFLT),intent(in):: sh_color
       integer(kind=PLINT),intent(in):: sh_cmap,sh_width,min_color
       integer(kind=PLINT),intent(in):: min_width,max_color,max_width
       real(kind=PLFLT),intent(in),dimension(*)::z,xg,yg
     end subroutine plshade1

     subroutine plcon1(z, nx, ny, kx, lx, ky, ly, clevel, nlevel,x,y)
       use plplot_params_module
       integer(kind=PLINT),intent(in)::nx,ny,kx,lx,ky,ly,nlevel
       real(kind=PLFLT),intent(in),dimension(*)::z,x,y
       real(kind=PLFLT),intent(in),dimension(*)::clevel
     end subroutine plcon1


     subroutine plend()
     end subroutine plend

  end interface

end module plplot_interface_module
