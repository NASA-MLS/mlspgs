! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MieTheory

! -------------------------------------------------------------------------  
! MODULE TO COMPUTE MIE EFFICIENCIES
! -------------------------------------------------------------------------

        use MLSCommon, only: r8   
	implicit none

        Private
        Public :: MieCoeff

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
 !---------------------------------------------------------------------------
      
contains

	SUBROUTINE MieCoeff(ISPI,f,t,nr,r,a,b,nab,nabr,bc)

	include 'constants.f9h'
        integer :: ISPI	                        ! cloud type (1=ICE, 2=WATER)
	real(r8) :: f 		                ! frequency in GHz
	real(r8) :: t 		                ! Temperature in K 
	real(r8) :: wl 	                        ! wavelength in meters
        integer :: nr	                        ! no of particle size
	integer :: nab	                        ! no of a/b terms
	real (r8) :: r(nr)	                ! particle radius
	integer :: nabr(nr)	                ! truncation number for a and b
				                ! 10um ---> 5; 2000 um ---> 20.
	complex(r8) :: a(nr,nab),b(nr,nab)	! Mie coefficients
	complex(r8) :: m		        ! refrative index
        complex(r8) :: a0,a1,w9,w0,w1,p1,p2,mx,mx1
        real(r8) :: x, x1
	real(r8) :: ab_err	                ! relative error in Mie efficiency
	parameter(ab_err = 1.e-3)
	real(r8) :: bc(3,nr)	! single particle Mie efficiencies (abs,scat,ext)
	real(r8) :: dab1, dab2	! a/b contribution to bc at i term
!... working space
	integer :: i,j

!... initialization
       wl=0.3_r8/f
       do i=1,nr
       do j=1,nab
        a(i,j)=cmplx(0.0)
        b(i,j)=cmplx(0.0)
       end do
       end do
	   
       do 12 j=1,nr

        if(ISPI.eq.1) call ukisub(f,t,m)
        if(ISPI.eq.2) call uksub(f,t,m)

!        x=dreal(2.*pi*r(j)/wl)*1.d-6    !jj

	x=real(2.*pi*r(j)/wl)*1.d-6
	
	bc(2,j) =0.
	bc(3,j) =0.

!... default of nabr is nab, the largest possible
        nabr(j)=nab
       
!... x1 is 1/x
!        x1=dreal(0.5d0*wl/pi/r(j))*1.d6    !jj
	x1=real(0.5*wl/pi/r(j))*1.d6

!        mx=dcmplx(m)*x     !jj
	mx=cmplx(m)*x
        mx1=1.0/mx

! ....... for CGI     !JJ
!        w9=dcmplx(dcos(x),-dsin(x))
!        w0=dcmplx(dsin(x),dcos(x))
!        a0=cdcos(mx)/cdsin(mx)

! ....   for zvi f95
	w9=cmplx(cos(x),-sin(x))
        w0=cmplx(sin(x),cos(x))
        a0=cos(mx)/sin(mx)


       do 10 i=1,nab
        w1=(2.0*i-1.)*x1*w0-w9
        a1=-i*mx1+1.0/(i*mx1-a0)
        p1=a1/m+i*x1
        p2=m*a1+i*x1

        a(j,i) = cmplx( (p1*real(w1)-real(w0))/(p1*w1-w0) )
        b(j,i) = cmplx( (p2*real(w1)-real(w0))/(p2*w1-w0) )

! ... the factor of 2./x^2 is left out in dab since nabr is usually
!	important for large x. For x<0.1, usually, nabr=2 is enough.
	dab1=(2.*i+1.)*((abs(a(j,i)))**2+(abs(b(j,i)))**2)
	dab2=(2.*i+1.)*real(a(j,i)+b(j,i))
	bc(2,j)=bc(2,j)+dab1
	bc(3,j)=bc(3,j)+dab2

! ... determine cutoff no. for higher order terms
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

       END SUBROUTINE MieCoeff

end module MieTheory

! $Log: MieTheory.f90,v      
