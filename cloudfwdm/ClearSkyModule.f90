! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module ClearSkyModule

! -------------------------------------------------------------------------  
! MLS CLEAR SKY RADIANCE MODEL
! -------------------------------------------------------------------------

      use Bill_GasAbsorption, only: get_beta_bill
      use GasAbsorption, only: GET_BETA
      use MLSCommon, only: r8
      use PrtMsg, only: HEADER
      use SpectraLines, only: SETUP_SPECTRA
      use SurfaceModel, only: SURFACE
      use SpectroscopyCatalog_m, only: CATALOG_T, LINES

      IMPLICIT NONE
      Private
      Public :: CLEAR_SKY

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
 !---------------------------------------------------------------------------
      
contains

      SUBROUTINE CLEAR_SKY(L,NU,TS,S,LORS,WIND,XZ,XP,XT,XQ,VMR, NS, &
                 &         F,RS,U,T,TAU,Z,TAU100, Catalog, Bill_Spectra,LosVel )

!======================================================
!     >>>>>>>>CLEAR-SKY RADIATION SCHEME<<<<<<<<<<
!
!     CALCULATE BACKGROUND ATMOSPHERIC RADIATION
!     LATEST UPDATE, J.JIANG, MAY 18, 2001
!======================================================

      INCLUDE 'spectra.f9h' 

!-------------------------------------------------
!     ATMOSPHERIC PROFILE PARAMETERS
!-------------------------------------------------

      INTEGER :: L, NU, I, J, LORS, NS, n_sps, Spectag

      REAL(r8) :: RS(NU/2),T(L),TAU(L),U(NU),Z(L),TAU100(L)
      REAL(r8) :: XZ(L+1),XP(L+1),XT(L+1),XQ(L+1)
      REAL(r8) :: VMR(NS,L+1),VMR1(NS)
      REAL(r8) :: DQ, P, DR, TS, S, WIND, F, LosVel

!-------------------------------------------------
!     SURFACE REFLECTIVITY
!-------------------------------------------------
! local variables

      REAL(r8):: RH(NU)                   ! HORIZONTAL
      REAL(r8):: RV(NU)                   ! VERTICAL
      REAL(r8):: X(NU)                    ! SCATTERING ANGLES

!-----------------------------------------------------
! Spectra Catalog 
!----------------------------------------------------

      Type(Catalog_T), INTENT(IN) :: Catalog(:)

      Logical :: Bill_Spectra

!-------------------------Begin Execution----------------------------------------

!      CALL HEADER(2)

      CALL SETUP_SPECTRA(QLG,V0,GSE,IST,WTH,NTH,DELTA,N1,   &
                  &      GAMMA,N2,MOL,NMOL,NCNT)

!-------------------------------------------------
!     GET SURFACE REFLECTIVITY
!-------------------------------------------------

      DO I=1,NU/2
         X(I)=U(I)
      ENDDO

      CALL SURFACE(F,TS,LORS,S,WIND,X,NU/2,NU,RH,RV)
      DO I=1,NU/2
         RS(I)=RH(I)
      ENDDO

      DO I=1,L

         Z(I)=XZ(I+1)-XZ(I)
         T(I)=(XT(I+1)+XT(I))*0.5

         P=(LOG10(XP(I+1))+LOG10(XP(I)))*0.5     
         P=10**P

         DQ=(LOG10( max(1.e-39_r8, XQ(I+1)) )+ &
         &         LOG10( max(1.e-39_r8, XQ(I)) ) )*0.5
         DQ= 10**DQ

         DO J=1,NS
            VMR1(J)=(LOG10( max(1.e-29_r8, VMR(J,I+1)) )+ &
                   & LOG10( max(1.e-29_r8, VMR(J,I))  )  )*0.5   !vmr cannot be zero
            VMR1(J)=10**VMR1(J)
         ENDDO

         If ( .not. Bill_Spectra ) then
           !---------------------------------
           ! Using default spectroscopy data
           ! --------------------------------
           CALL GET_BETA(QLG,V0,GSE,IST,WTH,NTH,DELTA,N1,GAMMA,N2,  &
                &        MOL,NMOL,NCNT,T(I),P,F,DQ,VMR1,DR,NS)   
                                             ! HERE DQ IS H2O MIXING RATIO
           TAU(I)=DR*Z(I)

           CALL GET_BETA(QLG,V0,GSE,IST,WTH,NTH,DELTA,N1,GAMMA,N2,  &
                &        MOL,NMOL,NCNT,T(I),P,F,100._r8,VMR1,DR,NS) 
                                             ! HERE DQ IS RELATIVE HUMIDITY!
           TAU100(I)=DR*Z(I)
         
         else
           !--------------------------------
           ! Using bill's spectroscopy data
           !--------------------------------
           call get_beta_bill (T(I), P, F, DQ, VMR1(1), DR, Catalog, LosVel)
           
           TAU(I)=DR*Z(I)

           call get_beta_bill (T(I), P, F, 100._r8, VMR1(1), DR, Catalog, LosVel)

           TAU100(I)=DR*Z(I)

         endif

      ENDDO

      END SUBROUTINE CLEAR_SKY

end module ClearSkyModule

! $Log$
! Revision 1.13  2001/11/15 23:53:02  jonathan
! clean log
!
! Revision 1.12  2001/11/15 23:52:28  jonathan
! add default_spectroscopy
!
! Revision 1.11  2001/11/15 00:54:14  jonathan
! minor changes
!
! Revision 1.10  2001/11/14 00:39:42  jonathan
! switch to use bill'd data
!
! Revision 1.9  2001/11/09 22:06:56  jonathan
! add Bill_GasAbsorption Module
!
! Revision 1.8  2001/11/09 18:07:36  jonathan
! add spectra catalog
!
! Revision 1.7  2001/10/04 23:35:15  dwu
! *** empty log message ***
!
! Revision 1.6  2001/09/21 15:51:37  jonathan
! modified F95 version
!

