
      SUBROUTINE CLEAR_SKY(L,NU,TS,S,LORS,SWIND,XZ,XP,XT,XQ,VMR, &
                 &         F,RS,U,T,TAU,Z,TAU100)

!======================================================
!     >>>>>>>>CLEAR-SKY RADIATION SCHEME<<<<<<<<<<
!
!     CALCULATE BACKGROUND ATMOSPHERIC RADIATION
!     LATEST UPDATE, J.JIANG, MAY 18, 2001
!======================================================

!-------------------------------------------------
!     ATMOSPHERIC PROFILE PARAMETERS
!-------------------------------------------------

!      INCLUDE 'spectra.h' !---------------------------------------------
!... spectral header file 
	integer no_line		! max no. of lines for all species
	integer no_mol		! max no. of molecules
	parameter (no_line= 1000, no_mol = 20)
	real*8 v0(no_line)
	real*4 gse(no_line)
	real*4 ist(no_line)
	real*4 wth(no_line)
	real*4 nth(no_line)
	real*4 qlg(3,no_line)
	real*4 delta(no_line)
	real*4 n1(no_line)
	real*4 gamma(no_line)
	real*4 n2(no_line)
	integer mol		! molecule mass number

	integer nmol		! no. of molecules
	integer ncnt(no_mol)	! no. of lines used for each molecule
! -----------------------------------------------------------------
      INTEGER L,NU,NP,I
      REAL RS(NU/2),T(L),TAU(L),U(NU),Z(L),TAU100(L)
      REAL XZ(L+1),XP(L+1),XT(L+1),XQ(L+1)
      REAL VMR(2,L+1),VMR1(5)

!-------------------------------------------------
!     SURFACE REFLECTIVITY
!-------------------------------------------------

      PARAMETER(NP=800)                   ! DEMENSION OF WORKING ARRAYS
      REAL RH(NP)                         ! HORIZONTAL
      REAL RV(NP)                         ! VERTICAL
      REAL X(NP)                          ! SCATTERING ANGLES

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

         DO J=1,5
            VMR1(J)=(VMR(J,I+1)+VMR(J,I))*0.5
         ENDDO

         CALL GET_BETA(QLG,V0,GSE,IST,WTH,NTH,DELTA,N1,GAMMA,N2,  &
              &        MOL,NMOL,NCNT,T(I),P,F,DQ,VMR1,DR)   
                                             ! HERE DQ IS H2O MIXING RATIO
         TAU(I)=DR*Z(I)

         CALL GET_BETA(QLG,V0,GSE,IST,WTH,NTH,DELTA,N1,GAMMA,N2,  &
              &        MOL,NMOL,NCNT,T(I),P,F,100.,VMR1,DR) 
                                             ! HERE DQ IS RELATIVE HUMIDITY!
         TAU100(I)=DR*Z(I)

      ENDDO

      RETURN
      END

! $Log: Clear_sky.f,v      



