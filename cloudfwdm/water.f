c         p ---- mb
c         t,td ---- K
c         RH --- %
c         SD --- g/m3
    
          SUBROUTINE RHtoSD(p,t,RH,SD)
          ES=6.10779*EXP(17.28*(t-273.15)/(t-36.))
          if (ES.lt.0.) ES=0
	  rou=p/2.87/t
          SD=622*ES/(p-ES)*RH*0.01*rou
          RETURN
          END

          SUBROUTINE RHtoQV(p,t,RH,QV)
          ES=6.10779*EXP(17.28*(t-273.15)/(t-36.))
          if (ES.lt.0.) ES=0
          QV=622*ES/(p-ES)*RH*0.01
          RETURN
          END

          SUBROUTINE SDtoRH(p,t,RH,SD)
          ES=6.10779*EXP(17.28*(t-273.15)/(t-36.))                    
          QS=622*ES/(p-ES) 
	  rou=p/2.87/t
	  Q=SD/rou 
          RH=Q/QS*100
          RETURN
          END

          SUBROUTINE TDtoRH(td,t,RH)                                      
          ES=6.10779*EXP(17.28*(t-273.15)/(t-36.))                        
          E=6.10779*EXP(17.28*(td-273.15)/(td-36.))                       
          RH=E/ES*100                                                     
          RETURN                                                          
          END 

          SUBROUTINE TDtoQV(td,t,QV)
          E=6.10779*EXP(17.28*(td-273.15)/(td-36.))
          if (E.lt.0.) E=0
          QV=622*E/(P-E)
          RETURN
          END 

          SUBROUTINE RHtoEV(p,t,RH,EV)
c... relative to water
c         ES=6.10779*EXP(17.28*(t-273.15)/(t-36.))
c... relative to ice

	  es=10**(-2667./t+10.555)
          if (ES.lt.0.) ES=0
          EV=ES*0.01*RH
          RETURN
          END

          SUBROUTINE SDtoEV(p,t,EV,SD)
          ES=6.10779*EXP(17.28*(t-273.15)/(t-36.))                    
          QS=622*ES/(p-ES) 
	  rou=p/2.87/t
	  Q=SD/rou 
          RH=Q/QS
          EV=RH*ES
          RETURN
          END

! $Log: water.f,v      

