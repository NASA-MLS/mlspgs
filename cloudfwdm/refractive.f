C******************************************************************
C* COMPLEX REFRACTIVE INDEX OF WATER DROPS AND ICE-PARTICLES      * 
C*  valid for freq=1-1000 GHz
c	taken from Liebe (1989) formulation
C* WI ---- 'W' OR 'I' ( WATER CASE OR ICE CASE )      ++ INPUT    *  
C* T  ---- TEMPERATURE IN  K                          ++ INPUT    *  
C* F  ---- FREQUENCY  IN  GHZ                         ++ INPUT    *
C* E  ---- REFRACTIVE INDEX (E=M**2)  ++ OUTPUT                   *  
C******************************************************************
C                                                                        
          subroutine comx(wi,t,f,e) 
          complex e
          character*1 wi
		  
          th=300./t
c          wl=30./f
          if ( wi .eq. 'W') then                                          
            goto 100                                                  
          else                                        
            goto 200                                                  
          end if                                                           
                                                         
c... from Liebe (1989)

 100    if (t .lt. 233.15 ) then
           t = 233.15
        endif

        fp=20.09-142.4*(th-1)+294*(th-1.)**2
	fs=590.-1500.*(th-1)
	e0=77.66+103.3*(th-1)
	e1=5.48
	e2=3.51
	x=(e0-e1)/(1+(f/fp)**2)+(e1-e2)/(1+(f/fs)**2)+e2
	y=(e0-e1)*f/fp/(1+(f/fp)**2)+(e1-e2)*f/fs/(1+(f/fs)**2)
	
c..Lu          C=T-273.16 
c	  IF(C.LT.-50.) C=-50.
c  100     E0=5.27137+0.0216474*C-0.00131198*C*C                            
c          A=-16.8129/(C+273.)+0.0609265                                    
c          S=12.5664E+08                                                   
c          R=0.00033836*EXP(2513.98/(C+273.))                               
c          E1=78.54*(1.0-4.579E-3*(C-25.)+
c     1       1.19E-5*(C-25.)**2-2.8E-8*(C-25.)**3)   
c          X1=(R/WL)**(1.0-A)                                                 
c          XX=X1*X1                                                         
c          X2=A*PAI/2                                                       
c          XS=ASIN(X2)                                                       
c          XC=ACOS(X2)                                                       
c          BV=18.8496E+10                                                   
c          O1=(E1-E0)*(1.0+X1*XS)                                             
c          O2=(E1-E0)*X1*XC                                                 
c          O0=1.0+2.0*X1*XS+XX                                                
c          X=E0+O1/O0                                                       
c          Y=S*WL/BV+O2/O0  

        e=x*(1.0,0.0)+y*(0.0,-1.0)

        return              !jj

c... from Hufford (1991) model

  200	a=(50.4+62.*(th-1))*1e-4*exp(-22.1*(th-1.))
	b=(0.633/th-0.131)*1e-4+(7.36e-4*th/(th-0.9927))**2
	x=3.15
	y=a/f+b*f

c        write(*,*)y,f,t

        e=x*(1.0,0.0)+y*(0.0,-1.0)
        return

        END
c...
         SUBROUTINE UKSUB(F,T,M)
         COMPLEX E,M
         CALL COMX('W',T,F,E)
         M=CSQRT(E)
c         write(33,*)m,f
	 RETURN
	 END
c...
         SUBROUTINE UKISUB(F,T,M)
         COMPLEX E, M
         CALL COMX('I',T,F,E)
         M=CSQRT(E)
c         write(43,*)m,f
	 RETURN
	 END

! $Log: refractive.f,v      
