! SUBROUTINES FOR CALCULATIONS OF SURFACE EMISSIVITY
       SUBROUTINE SURFACE(F,TS,IS,S,W,X,N,NP,RH,RV)
! IS --- surface type:
!        0 - Ocean
!        1 - Land
!        2 - SeaIce
!        3 - Snow
! TS --- surface temperature in K
!  F --- Frequency in GHz
!  S --- Salinity in per thousand if ocean,  
!        Emissivity value, if land (negative means default,0.9)
!        Ice type, if sea ice, 0=new,1=2nd year,2=multiyear
!        Snow type, if snow, 0=wet,1=dry,2=refrozen
!  W --- Wind Speed in m/s, if ocean
!        Not used, else

	real X(NP)	! scattering angles
	real RH(NP)	! reflectivity for horizontal polarization
	real RV(NP)	! reflectivity for vertical polarization

       if (IS.eq.0) then
!         call ASSEA1(F,TS,S,W,X,N,NP,RH,RV)
!... simple model
	 do i=1,np
	 rh(i) = 0.8
	 rv(i) = 0.8
	 enddo
       else if(IS.lt.4 .and. IS.gt.0) then
         call ASSEAN(F,TS,IS,S,X,N,NP,RH,RV)
       else
         print*,'Check, surface type:',IS
         stop
       end if
       return
       end
    
       SUBROUTINE ASSEA(F,TS,S,X,N,NP,RH,RV)
!      X ----  COS(SITA)
!      S ---- SALINTY PER MILI
       REAL X(NP),RH(NP),RV(NP)
       COMPLEX E,P1,P2

       T=TS-273.16
       D=25.-T
       B=2.033E-2+1.266E-4*D+2.464E-6*D*D      &
     &  -S*(1.849E-5-2.551E-7*D+2.551E-8*D*D)
       S25=S*(0.182521-1.46192E-3*S+2.09324E-5*S*S-   &
     &  1.28205E-7*S*S*S)
       SS=S25*EXP(-D*B)
       E0=8.854E-12
       E9=4.9
       A=0.0
       ES=87.134-1.949E-1*T-1.276E-2*T*T+2.491E-4*T*T*T
       A1=1.+1.613E-5*S*T-3.656E-3*S+3.21E-5*S*S-      &
     &    4.232E-7*S*S*S
       ES=ES*A1
       TA=1.768E-11-6.086E-13*T+1.104E-14*T*T          &
     &          -8.111E-17*T*T*T
       A1=1.+2.282E-5*S*T-7.638E-4*S-                  &
     &    7.76E-6*S*S+1.105E-8*S*S*S
       TA=TA*A1
       PAI=3.1415926
       A1=2.0*PAI*F*1.E+9
       P1=(A1*TA*(0.0,1.0))**(1.-A)+1.0
       P2=(0.0,1.0)*SS/A1/E0
       E=E9+(ES-E9)/P1-P2
       DO I=1,N
        CS=X(I)
        SS2=1.-CS*CS
        RH(I)=(CABS((CS-CSQRT(E-SS2))/      &
     &        (CS+CSQRT(E-SS2))))**2
        RV(I)=(CABS((E*CS-CSQRT(E-SS2))/    &
     &        (E*CS+CSQRT(E-SS2))))**2
       END DO
       RETURN
       END
!CF
       SUBROUTINE ASSEA1(F,TS,S,W,X,N,NP,RH,RV)
       REAL X(NP),RH(NP),RV(NP)
       DATA a,b,c,d,e/0.117,-2.09e-3,7.32e-2,0.115,    &
     &                3.8e-5/

       CALL ASSEA(F,TS,S,X,N,NP,RH,RV)
       IF(W.LE.0) RETURN
       
       DO IN=1,N
        CA=ACOS(X(IN))*180/3.1415926
        DTRV=W*(a+b*exp(c*CA))*F
        DTRH=W*(d+e*CA*CA)*F
        ROUV=RV(IN)-DTRV/TS
        ROUH=RH(IN)-DTRH/TS
        FV=1-9.946e-4*CA+3.218e-5*CA*CA           &
     &      -1.187e-6*CA*CA*CA+7e-20*CA**10
        FH=1-1.748e-3*CA-7.336e-5*CA*CA+1.044e-7*CA**3
        ROUFV=1-(208+F*1.29)/TS*FV
        ROUFH=1-(208+F*1.29)/TS*FH
        FU=2.95e-6*W**3.52
        RH(IN)=ROUH*(1-FU)+ROUFH*FU
        RV(IN)=ROUV*(1-FU)+ROUFV*FU
       END DO
       RETURN
       END
!CF
       SUBROUTINE ASSEA3(F,TS,S,W,X,N,NP,RH,RV)
!     good only for ssm/i channels
       REAL X(NP),RH(NP),RV(NP),C(8,4),DES(8),CG(4)
       DATA CF,W0,CA0,TS0/8.E-3,7.,53.,273./
       DATA C/-0.556,0.406,-0.670,0.479,      &
     & -0.811,0.473,-0.723,0.358,             &
     &  0.357,-0.108,0.455,-0.175,0.551,-0.160,  &
     &  0.404,-0.0351,                           &
     & -0.0312,0.0128,-0.0446,0.0283,            &
     & -0.0365,0.0312,-0.00735,0.0309,           &
     &  0.0106,0.00153,0.0232,-0.0131,           &
     &  0.0149,-0.0150,-0.0126,-0.0121/
!      DATA DES/0.004,0.00354,-0.00451,0.0,
!     1        -0.0124,-0.0125,-0.00189,0.05415/
       DATA DES/8*0.0/
       DATA CG/3.49E-3,3.6E-3,4.85E-3,6.22E-3/
       CALL ASSEA(F,TS,S,X,N,NP,RH,RV)
       IF(W.LE.0) RETURN

       IF(F.LT.20) THEN
         KF=1
       ELSE IF(F.LT.30) THEN
         KF=3
       ELSE IF(F.LT.40) THEN
         KF=5
       ELSE IF(F.LT.100) THEN
         KF=7
       ELSE
         RETURN
       END IF

       TD=TS/TS0
       FU=0
       IF(W.GT.W0) FU=CF*(W-W0)
       CGA=CG(1+INT(KF/2))
       G2=CGA*W
       DO IN=1,N
        CA=ACOS(X(IN))*180/3.1415926
        DC=CA-CA0
        DERV=G2*(C(KF,1)+C(KF,2)*TD+C(KF,3)*DC+C(KF,4)*DC*TD)
        DERH=G2*(C(KF+1,1)+C(KF+1,2)*TD+C(KF+1,3)*DC+        &
     &           C(KF+1,4)*DC*TD)
        E0H=1-RH(IN)
        E0V=1-RV(IN)
        EH=(E0H+DERH)*(1-FU)+FU+DES(1+KF)
        EV=(E0V+DERV)*(1-FU)+FU+DES(KF)
        RH(IN)=1-EH
        RV(IN)=1-EV
       END DO
       RETURN
       END

! emissivity for non-ocean surface
!
       SUBROUTINE ASSEAN(F,TS,IS,S,X,N,NP,RH,RV)

! IS -- 1 = land, w=emissivity of it,0.9 if negative given
!       2 = ice, S---0: New, 1: Second Year, 2: Multiyear
!       3 = snow,S---0: Wet, 1: Dry, 2: Refrozen

       REAL X(NP),RH(NP),RV(NP),CICE(4,3),CSNOW(4,3)
       save CICE,CSNOW
       DATA CICE/0.95,0.95,99,0,   &
     &           0.93,0.83,31,2,   &
     &           0.92,0.64,31,2/
       DATA CSNOW/0.76,0.99, 9,2,  &
     &            0.90,0.75,33,3,  &
     &            0.97,0.53,32,4/

       if (IS.eq.1) then
         e=S
         if(e.lt.0) e=0.9
       else if(IS.eq.2) then
         k=nint(S)+1
         e=(CICE(1,k)+CICE(2,k)*(F/CICE(3,k))**CICE(4,k))  &
     &     /(1+(F/CICE(3,k))**CICE(4,k))
        
       else if(IS.eq.3) then
         k=nint(S)+1
         e=(CSNOW(1,k)+CSNOW(2,k)*(F/CSNOW(3,k))**CSNOW(4,k))  &
     &     /(1+(F/CSNOW(3,k))**CSNOW(4,k))
       end if

       do i=1,n
        RH(i)=1-e
        RV(i)=1-e
       end do

       return
       end

! $Log: surface.f90,v 
