! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module RefractiveIndex

!==================================================================
! COMPUTE COMPLEX REFRACTIVE INDEX OF WATER DROPS AND ICE-PARTICLES       
!
!* WI ---- 'W' OR 'I' ( WATER CASE OR ICE CASE )      ++ INPUT    *  
!* T  ---- TEMPERATURE IN  K                          ++ INPUT    *  
!* F  ---- FREQUENCY  IN  GHZ                         ++ INPUT    *
!* E  ---- REFRACTIVE INDEX (E=M**2)  ++ OUTPUT                   *  
!==================================================================

         use MLSCommon, only: r8
 	 implicit none

         Private :: COMX
         Public :: UKSUB, UKISUB

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
 !---------------------------------------------------------------------------
      
contains

         SUBROUTINE UKSUB(F,T,M)
         !-----------------------
         real(r8) :: F, T
         COMPLEX(r8) :: E,M
         CALL COMX('W',T,F,E)
         M=SQRT(E)
	 
	 END SUBROUTINE UKSUB

         SUBROUTINE UKISUB(F,T,M)
         ! ----------------------
         real(r8) :: F, T
         COMPLEX(r8) :: E, M
         CALL COMX('I',T,F,E)
         M=SQRT(E)

	 END SUBROUTINE UKISUB

          subroutine comx(wi,t,f,e) 

          complex(r8) :: e,c
          character*1 wi
          real(r8) :: t, f, th, fp, fs, e0, e1, e2, x, y
          real(r8) :: a, b
		  
          th=300./t
!          wl=30./f
          if ( wi .eq. 'W') then                                          
            goto 100                                                  
          else                                        
            goto 200                                                  
          end if                                                           
                                                         
!... from Liebe (1989)

 100    if (t .lt. 233.15 ) then
           t = 233.15
        endif

        fp=20.09-142.4*(th-1.)+294*(th-1.)**2
	fs=590.-1500.*(th-1.)
	e0=77.66+103.3*(th-1.)
	e1=5.48
	e2=3.51
	x=(e0-e1)/(1+(f/fp)**2)+(e1-e2)/(1+(f/fs)**2)+e2
	y=(e0-e1)*f/fp/(1+(f/fp)**2)+(e1-e2)*f/fs/(1+(f/fs)**2)
	
!..Lu          C=T-273.16 
!	  IF(C.LT.-50.) C=-50.
!  100     E0=5.27137+0.0216474*C-0.00131198*C*C                            
!          A=-16.8129/(C+273.)+0.0609265                                    
!          S=12.5664E+08                                                   
!          R=0.00033836*EXP(2513.98/(C+273.))                               
!          E1=78.54*(1.0-4.579E-3*(C-25.)+
!     1       1.19E-5*(C-25.)**2-2.8E-8*(C-25.)**3)   
!          X1=(R/WL)**(1.0-A)                                                 
!          XX=X1*X1                                                         
!          X2=A*PAI/2                                                       
!          XS=ASIN(X2)                                                       
!          XC=ACOS(X2)                                                       
!          BV=18.8496E+10                                                   
!          O1=(E1-E0)*(1.0+X1*XS)                                             
!          O2=(E1-E0)*X1*XC                                                 
!          O0=1.0+2.0*X1*XS+XX                                                
!          X=E0+O1/O0                                                       
!          Y=S*WL/BV+O2/O0  

        e=x*(1.0,0.0)+y*(0.0,-1.0)

        return              !jj

!... from Hufford (1991) model

  200	a=(50.4+62.*(th-1))*1.e-4*exp(-22.1*(th-1.))
	b=(0.633/th-0.131)*1e-4+(7.36e-4*th/(th-0.9927))**2

!... additional term from Mishima
	c=1.16e-11_r8

	x=3.15
!        y=a/f+b*f
	y=a/f+b*f+c*f**3
!        write(*,*)y,f,t

        e=x*(1.0,0.0)+y*(0.0,-1.0)
        return

        END subroutine comx

end module RefractiveIndex

! $Log$
! Revision 1.4  2002/04/15 22:22:32  jonathan
! check bug
!
! Revision 1.3  2002/04/15 17:13:28  jonathan
! add mishima term
!
! Revision 1.2  2001/09/21 15:51:37  jonathan
! modified F95 version
!
