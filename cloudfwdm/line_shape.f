c--------------------------------------------------
	real*8 function myshape(v0,ff,yy,twth0)
	implicit none
        real pi
        parameter (pi=3.1415926)
	real*8 voffm		! frequency difference
	real*8 voffp		! frequency difference
	real*8 ff		! frequency in MHz
	real*8 v0		! line frequency in MHz
	real twth0		! line width in MHz/mb
	real yy 		! interference coeff.
	real*8 twthm, twthp
c
	  voffm = v0 - ff
	  voffp = v0 + ff
	  twthm = twth0 - voffm*yy*1.d0
	  twthp = twth0 - voffp*yy*1.d0
c... Gross function
c	myshape = 4.*ff*v0*twth0/pi/((voffm*voffp)**2+4*(ff*twthm)**2)
c	return

c... Van Vleck-Weisskopf
 	myshape=twthm / (voffm**2 + twth0**2) + 
     +   twthp/(voffp**2 + twth0**2)
 	myshape = myshape * (ff/v0) /pi
c...
	return
	end
