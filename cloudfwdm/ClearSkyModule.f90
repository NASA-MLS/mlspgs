! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module ClearSkyModule

! -------------------------------------------------------------------------  
! MLS CLEAR SKY RADIANCE MODEL
! -------------------------------------------------------------------------

      use Bill_GasAbsorption, only: get_beta_bill
      use GasAbsorption, only: GET_BETA
      use MLSCommon, only: r8
      use SpectraLines, only: SETUP_SPECTRA
      use SurfaceModel, only: SURFACE
      use SpectroscopyCatalog_m, only: CATALOG_T
      USE Intrinsic,    only: l_clear

      IMPLICIT NONE
      Private
      Public :: CLEAR_SKY

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
  private :: not_used_here 
 !---------------------------------------------------------------------------
      
contains

      SUBROUTINE CLEAR_SKY(L,NU,TS,S,LORS,WIND,XZ,XP,XT,XQ,VMR, NS, &
                 &         F,RS,U,T,Z,TAU,tau_wet, tau_dry, Catalog, &
                 &         Bill_Spectra,LosVel, i_saturation )

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

      INTEGER, intent(in) :: NS                ! NO. OF CHEMICAL SPECIES
      INTEGER, intent(in) :: NU                ! 2 x NO. OF pts
      INTEGER, intent(in) :: L                 ! NO. OF (levels?)
      INTEGER, intent(in) :: LORS              ! Surface type index (land/sea)
      REAL(r8), intent(in) :: U(NU)            ! Scattering angles
      REAL(r8), intent(in) :: F                ! Frequency
      REAL(r8), intent(in) :: TS               ! Surface temperature
      REAL(r8), intent(in) :: S                ! Salinity
      REAL(r8), intent(in) :: WIND             ! Wind speed
      REAL(r8), intent(in) :: LosVel
      REAL(r8), intent(in) :: XZ(L+1)
      REAL(r8), intent(in) :: XP(L+1)
      REAL(r8), intent(in) :: XT(L+1)
      REAL(r8), intent(in) :: XQ(L+1)
      REAL(r8), intent(in) :: VMR(NS-1,L+1)
      REAL(r8), intent(out) :: RS(NU/2)        ! Surface reflectivity
      REAL(r8), intent(out) :: tau_wet(L)
      REAL(r8), intent(out) :: tau_dry(L)
      REAL(r8), intent(out) :: TAU(L)
      REAL(r8), intent(out) :: T(L)
      REAL(r8), intent(out) :: Z(L)
      Integer, intent(in) :: i_saturation           ! clear-sky flag

!-----------------------------------------------------
! Spectra Catalog 
!----------------------------------------------------

      Type(Catalog_T), INTENT(IN) :: Catalog(:)

      Logical, INTENT(IN) :: Bill_Spectra

!-------------------------------------------------
!     SURFACE REFLECTIVITY
!-------------------------------------------------
! local variables

      ! Reflectivities
      REAL(r8):: RH(NU)                   ! HORIZONTAL
      REAL(r8):: RV(NU)                   ! VERTICAL
      ! Other stuff
      REAL(r8) :: P
      REAL(r8) :: DQ
      REAL(r8) :: DR
      REAL(r8) :: VMR1(NS-1)
      REAL(r8):: X(NU)                    ! SCATTERING ANGLES
      integer :: I
      integer :: J


!-------------------------Begin Execution--------------------------------------

!      CALL HEADER(2)

      QLG   = 0.0_r8
      V0    = 0.0_r8
      GSE   = 0.0_r8
      IST   = 0.0_r8
      WTH   = 0.0_r8
      NTH   = 0.0_r8
      DELTA = 0.0_r8
      N1    = 0.0_r8
      GAMMA = 0.0_r8
      N2    = 0.0_r8
      tau   = 0.0_r8
      tau_wet = 0.0_r8
      tau_dry = 0.0_r8

      If ( .not. Bill_Spectra ) then

        CALL SETUP_SPECTRA(QLG,V0,GSE,IST,WTH,NTH,DELTA,N1,   &
                    &      GAMMA,N2,MOL,NMOL,NCNT)
      Endif
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

         DQ=(LOG10( max(1.e-8_r8, XQ(I+1)) )+ &
         &         LOG10( max(1.e-8_r8, XQ(I)) ) )*0.5
         DQ= 10**DQ

         DO J=1,NS-1
            VMR1(J)=(LOG10( max(1.e-19_r8, VMR(J,I+1)) )+ &
                   & LOG10( max(1.e-19_r8, VMR(J,I))  )  )*0.5   ! vmr cannot
            VMR1(J)=10**VMR1(J)                                  ! be zero
         ENDDO

         If ( .not. Bill_Spectra ) then
           !---------------------------------
           ! Using default spectroscopy data
           ! --------------------------------
           CALL GET_BETA(QLG,V0,GSE,IST,WTH,NTH,DELTA,N1,GAMMA,N2,  &
                &        NMOL,NCNT,T(I),P,F,DQ,VMR1,DR,NS )   
                                             ! HERE DQ IS H2O MIXING RATIO
           TAU(I)=DR*Z(I)

           IF (i_saturation /= l_clear ) then       ! save time for clear sky
           CALL GET_BETA(QLG,V0,GSE,IST,WTH,NTH,DELTA,N1,GAMMA,N2,  &
                &        NMOL,NCNT,T(I),P,F,110._r8,VMR1,DR,NS ) 
                                             ! HERE 110 IS RELATIVE HUMIDITY!
           tau_wet(I)=DR*Z(I)
         
           CALL GET_BETA(QLG,V0,GSE,IST,WTH,NTH,DELTA,N1,GAMMA,N2,  &
                &        NMOL,NCNT,T(I),P,F,0.0_r8,VMR1,DR,NS ) 
                                             ! HERE DQ IS vmr
           tau_dry(I)=DR*Z(I)
           Endif
         
         else
           !--------------------------------
           ! Using bill's spectroscopy data
           !--------------------------------
           call get_beta_bill (T(I), P, F, DQ, VMR1, &
            & NS, DR, Catalog, LosVel)           
           TAU(I)=DR*Z(I)
           
           IF (i_saturation /= l_clear) then    ! save time for clear sky
           call get_beta_bill (T(I), P, F, 110._r8, VMR1, &
            & NS, DR, Catalog, LosVel)
           tau_wet(I)=DR*Z(I)

           call get_beta_bill (T(I), P, F, 0.0_r8, VMR1, &
            & NS, DR, Catalog, LosVel)
           tau_dry(I)=DR*Z(I)
           Endif
         endif

      ENDDO

      END SUBROUTINE CLEAR_SKY

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module ClearSkyModule

! $Log$
! Revision 1.23  2003/04/10 20:25:35  dwu
! make i_saturation as a verbal argument
!
! Revision 1.22  2003/02/01 06:43:06  dwu
! some fixes
!
! Revision 1.21  2003/01/30 23:48:02  dwu
! add icon=-4 for min Tb and some clearups
!
! Revision 1.20  2003/01/23 00:19:09  pwagner
! Some cosmetic only (or so I hope) changes
!
! Revision 1.19  2002/12/18 16:08:58  jonathan
! minor changes
!
! Revision 1.18  2002/10/08 17:08:07  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 1.17  2002/08/22 00:13:29  jonathan
! upgrade to include more molecules
!
! Revision 1.16  2002/08/08 22:45:57  jonathan
! newly improved version
!
! Revision 1.15  2001/11/20 19:36:46  jonathan
! some changes to save CPU
!
! Revision 1.14  2001/11/16 00:40:52  jonathan
! add losVel
!
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

