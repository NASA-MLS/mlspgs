
      SUBROUTINE CLEAR_SKY(L,NU,TS,S,LORS,WIND,XZ,XP,XT,XQ,VMR, NS, &
                 &         F,RS,U,T,TAU,Z,TAU100)

!======================================================
!     >>>>>>>>CLEAR-SKY RADIATION SCHEME<<<<<<<<<<
!
!     CALCULATE BACKGROUND ATMOSPHERIC RADIATION
!     LATEST UPDATE, J.JIANG, MAY 18, 2001
!======================================================

      use MLSCommon, only: r8
      INCLUDE 'spectra.f9h' 

!-------------------------------------------------
!     ATMOSPHERIC PROFILE PARAMETERS
!-------------------------------------------------

      INTEGER :: L, NU, NP, I, LORS, NS
      REAL(r8) :: RS(NU/2),T(L),TAU(L),U(NU),Z(L),TAU100(L)
      REAL(r8) :: XZ(L+1),XP(L+1),XT(L+1),XQ(L+1)
      REAL(r8) :: VMR(NS,L+1),VMR1(NS)
      REAL(r8) :: DQ, P, DR, TS, S, WIND, F

!-------------------------------------------------
!     SURFACE REFLECTIVITY
!-------------------------------------------------

      PARAMETER(NP=800)                   ! DEMENSION OF WORKING ARRAYS
      REAL(r8):: RH(NP)                   ! HORIZONTAL
      REAL(r8):: RV(NP)                   ! VERTICAL
      REAL(r8):: X(NP)                    ! SCATTERING ANGLES

!------------------------------------------------------------------------
      CALL HEADER(2)
      
      CALL SETUP_SPECTRA(QLG,V0,GSE,IST,WTH,NTH,DELTA,N1,   &
                  &      GAMMA,N2,MOL,NMOL,NCNT)

!-------------------------------------------------
!     GET SURFACE REFLECTIVITY
!-------------------------------------------------

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

         DO J=1,NS
            VMR1(J)=(VMR(J,I+1)+VMR(J,I))*0.5
         ENDDO

         CALL GET_BETA(QLG,V0,GSE,IST,WTH,NTH,DELTA,N1,GAMMA,N2,  &
              &        MOL,NMOL,NCNT,T(I),P,F,DQ,VMR1,DR,NS)   
                                             ! HERE DQ IS H2O MIXING RATIO
         TAU(I)=DR*Z(I)

         CALL GET_BETA(QLG,V0,GSE,IST,WTH,NTH,DELTA,N1,GAMMA,N2,  &
              &        MOL,NMOL,NCNT,T(I),P,F,100._r8,VMR1,DR,NS) 
                                             ! HERE DQ IS RELATIVE HUMIDITY!
         TAU100(I)=DR*Z(I)

      ENDDO

      RETURN
      END

! $Log: Clear_sky.f90,v      



