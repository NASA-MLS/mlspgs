	subroutine Miecoeff(ISPI,f,t,nr,r,a,b,nab,nabr,bc)
	implicit none
	include 'constants.h'
        integer ISPI	! cloud type (1=ICE, 2=WATER)
	real f 		! frequency in GHz
	real t 		! Temperature in K 
	real wl 	! wavelength in meters
        integer nr	! no of particle size
	integer nab	! no of a/b terms
	real r(nr)	! particle radius
	integer nabr(nr)	! truncation number for a and b
				! 10um ---> 5; 2000 um ---> 20.
	complex a(nr,nab),b(nr,nab)	! Mie coefficients
	complex m		! refrative index
        complex*16 a0,a1,w9,w0,w1,p1,p2,mx,mx1
        real*8  x,x1
	real ab_err	! relative error in Mie efficiency
	parameter(ab_err = 1.e-3)
	real bc(3,nr)	! single particle Mie efficiencies (abs,scat,ext)
	real dab1, dab2	! a/b contribution to bc at i term
c... working space
	integer i,j

c... initialization
       wl=0.3/f
       do i=1,nr
       do j=1,nab
        a(i,j)=cmplx(0.0)
        b(i,j)=cmplx(0.0)
       end do
       end do
	   
       do 12 j=1,nr

        if(ISPI.eq.1) call ukisub(f,t,m)
        if(ISPI.eq.2) call uksub(f,t,m)

c        x=dreal(2.*pi*r(j)/wl)*1.d-6    !jj

	x=real(2.*pi*r(j)/wl)*1.d-6
	
	bc(2,j) =0.
	bc(3,j) =0.

c... default of nabr is nab, the largest possible
        nabr(j)=nab
       
c... x1 is 1/x
c        x1=dreal(0.5d0*wl/pi/r(j))*1.d6    !jj
	x1=real(0.5d0*wl/pi/r(j))*1.d6

c        mx=dcmplx(m)*x     !jj
	mx=cmplx(m)*x
        mx1=1.0d0/mx

c ....... for CGI     !JJ
c        w9=dcmplx(dcos(x),-dsin(x))
c        w0=dcmplx(dsin(x),dcos(x))
c        a0=cdcos(mx)/cdsin(mx)

c ....   for zvi f95
	w9=cmplx(cos(x),-sin(x))
        w0=cmplx(sin(x),cos(x))
        a0=cos(mx)/sin(mx)


       do 10 i=1,nab
        w1=(2.0d0*i-1.d0)*x1*w0-w9
        a1=-i*mx1+1.0d0/(i*mx1-a0)
        p1=a1/m+i*x1
        p2=m*a1+i*x1

c        a(j,i) = cmplx( (p1*dreal(w1)-dreal(w0))/(p1*w1-w0) )   !jj
c        b(j,i) = cmplx( (p2*dreal(w1)-dreal(w0))/(p2*w1-w0) )   !jj

        a(j,i) = cmplx( (p1*real(w1)-real(w0))/(p1*w1-w0) )
        b(j,i) = cmplx( (p2*real(w1)-real(w0))/(p2*w1-w0) )


c ... the factor of 2./x^2 is left out in dab since nabr is usually
c	important for large x. For x<0.1, usually, nabr=2 is enough.
	dab1=(2.*i+1.)*((cabs(a(j,i)))**2+(cabs(b(j,i)))**2)
	dab2=(2.*i+1.)*real(a(j,i)+b(j,i))
	bc(2,j)=bc(2,j)+dab1
	bc(3,j)=bc(3,j)+dab2

C ... determine cutoff no. for higher order terms
	if(i .ge. 2 ) then
	dab1=dab1/bc(2,j)
	dab2=dab2/bc(3,j)
	if(abs(dab1) .lt. ab_err .and. abs(dab2) .lt. ab_err) then
		nabr(j) = i
		goto 11		! stop looping and use i as nabr
	endif
	endif

        a0=a1
	w9=w0
	w0=w1
 10    continue
 11	bc(2,j)=bc(2,j)*2/x/x
	bc(3,j)=bc(3,j)*2/x/x
 12    continue

       RETURN
       END

! $Log: miecoeff.f,v      



