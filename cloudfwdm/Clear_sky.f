C======================================================
C     >>>>>>>>CLEAR-SKY RADIATION SCHEME<<<<<<<<<<
C
C     CALCULATE BACKGROUND ATMOSPHERIC RADIATION
C======================================================

      SUBROUTINE CLEAR_SKY(L,NU,TS,S,LORS,SWIND,XZ,XP,XT,XQ,VMR,
     >                     F,RS,U,T,TAU,Z,TAU100)

C-------------------------------------------------
C     ATMOSPHERIC PROFILE PARAMETERS
C-------------------------------------------------

      INCLUDE 'spectra.h'
      INTEGER L,NU,NP,I
      REAL RS(NU/2),T(L),TAU(L),U(NU),Z(L),TAU100(L)
      REAL XZ(L+1),XP(L+1),XT(L+1),XQ(L+1)
      REAL VMR(5,L+1),VMR1(5)

C-------------------------------------------------
C     SURFACE REFLECTIVITY
C-------------------------------------------------

      PARAMETER(NP=800)                   ! DEMENSION OF WORKING ARRAYS
      REAL RH(NP)                         ! HORIZONTAL
      REAL RV(NP)                         ! VERTICAL
      REAL X(NP)                          ! SCATTERING ANGLES

C------------------------------------------------------------------------
      CALL HEADER(2)

      CALL SETUP_SPECTRA(QLG,V0,GSE,IST,WTH,NTH,DELTA,N1,
     >                   GAMMA,N2,MOL,NMOL,NCNT)

C-------------------------------------------------
C     GET SURFACE REFLECTIVITY
C-------------------------------------------------

      DO I=1,NU/2
         X(I)=U(I)
      ENDDO
      CALL SURFACE(F,TS,LORS,S,WIND,X,NU/2,NP,RH,RV)
      DO I=1,NU/2
         RS(I)=RH(I)
      ENDDO

      DO I=1,L
         Z(I)=XZ(I+1)-XZ(I)
         T(I)=(XT(I+1)+XT(I))*0.5
         P=(XP(I+1)+XP(I))*0.5               !!! NEED TO CHANGE TO P(I)
         DQ=(XQ(I+1)+XQ(I))*0.5

         DO J=1,5
            VMR1(J)=(VMR(J,I+1)+VMR(J,I))*0.5
         ENDDO

         CALL GET_BETA(QLG,V0,GSE,IST,WTH,NTH,DELTA,N1,GAMMA,N2,
     >                 MOL,NMOL,NCNT,T(I),P,F,DQ,VMR1,DR)   
                                             ! HERE DQ IS H2O MIXING RATIO
         TAU(I)=DR*Z(I)

         CALL GET_BETA(QLG,V0,GSE,IST,WTH,NTH,DELTA,N1,GAMMA,N2,
     >                 MOL,NMOL,NCNT,T(I),P,F,100.,VMR1,DR) 
                                             ! HERE DQ IS RELATIVE HUMIDITY!
         TAU100(I)=DR*Z(I)


      ENDDO
      RETURN
      END


