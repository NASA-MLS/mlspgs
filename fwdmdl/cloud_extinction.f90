! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
 
module Cloud_Extinction

  use MLSCommon,               only: r8, rp

  IMPLICIT NONE
  private
  public ::  get_beta_cloud

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
  private :: not_used_here 
 !---------------------------------------------------------------------------

contains 

  SUBROUTINE get_beta_cloud (frequency, temperature, pressure, &
                          &  WC, IPSD, N, NU, NUA, NAB, NR,    &
                          &  RC_TOT, W0, PHH                   )

!--------------------------------INPUT PARAMETERS----------------------------------

      INTEGER, intent(in) :: IPSD         ! SIZE-DISTRIBUTION FLAG
                                          ! 1000:     MHDISTRIBUTION 
                                          ! 1100:     LIU-CURRY
                                          ! 2000-3900:MODIFIED GAMMA WITH VARIOUS De, alpha
                                          ! 4000-5900:KNOLLENBERG WITH VARIOUS b1
                                          ! 6000:     PSD FOR PSC

      INTEGER, intent(in) :: N            ! NUMBER OF CLOUD SPECIES
      INTEGER, intent(in) :: NU           ! NO. OF SCATTERING ANGLES
      INTEGER, intent(in) :: NUA          ! NO. OF SCAT. AZIMUTH ANGLES
      INTEGER, intent(in) :: NAB          ! MAX. NO. OF A,B, TERMS 
      INTEGER, intent(in) :: NR           ! NUMBER OF SIZE BINS

      REAL(r8),intent(in) :: TEMPERATURE  ! in Kelvin
      REAL(r8),intent(in) :: PRESSURE     ! in mbar
      REAL(r8),intent(in) :: FREQUENCY    ! in MegaHertz
      REAL(r8),intent(in) :: WC(N)        ! CLOUD WATER CONTENT
                                          ! N=1: ICE; N=2: LIQUID

!--------------------------------OUTPUT PARAMETERS--------------------------------

      REAL(r8),intent(out) :: RC_TOT(3)   ! TOTAL ABS.SCAT.EXT COEFFS.
      REAL(r8),intent(out) :: W0(N)       ! SINGLE SCATTERING ALBEDO
      REAL(r8),intent(out) :: PHH(N,NU)   ! PHASE FUNCTION 

!-------------------------------Internal parameters-------------------------------

      REAL(r8) :: CWC                     ! CLOUD WATER CONTENT (g/m3)
      REAL(r8) :: THETA(NU)               ! SCATTERING ANGLES
      REAL(r8) :: U(NU)                   ! COSINE OF SCATTERING ANGLES      
      REAL(r8) :: DU(NU)                  ! DELTA U
      REAL(r8) :: PHI(NUA)                ! SCATTERING AZIMUTH ANGLES
      REAL(r8) :: UA(NUA)                 ! COSINE OF SCAT AZIMUTH ANGLES
      REAL(r8) :: UI(NU,NU,NUA)           ! COSINE OF INCIDENT TB ANGLES
      REAL(r8) :: THETAI(NU,NU,NUA)       ! ANGLES FOR INCIDENT TB
      REAL(r8) :: P11(NU)
      REAL(r8) :: RC11(3)
      REAL(r8) :: DMA
      REAL(r8) :: PH1(NU)                 ! SINGLE PARTICLE PHASE FUNCTION
      REAL(r8) :: P(NAB,NU)               ! LEGENDRE POLYNOMIALS l=1
      REAL(r8) :: DP(NAB,NU)              ! Delt LEGENDRE POLYNOMIALS l=1
      REAL(r8) :: R(NR)                   ! PARTICLE RADIUS
      REAL(r8) :: RN(NR)                  ! NUMBER OF PARTICLES IN EACH BIN
      REAL(r8) :: BC(3,NR)                ! SINGLE PARTICLE ABS/SCAT/EXT COEFFS 
      REAL(r8) :: RC_TMP(N,3)
      REAL(r8) :: DDm(N)                         

      COMPLEX(r8) A(NR,NAB),B(NR,NAB)     ! MIE COEFFICIENCIES

      INTEGER :: NABR(NR)                 ! TRUNCATION NUMBER FOR A, B
      INTEGER :: ISPI, J
!===================================================================================

      CALL ANGLE(THETA,U,DU,NU,PHI,UA,NUA,UI,THETAI)

      DO ISPI=1,N

         CWC = WC(ISPI)
         IF (CWC .ne. 0._r8 ) then           
            CWC = MAX(1.E-9_r8, abs(CWC))

            CALL CLOUDY_SKY (ISPI, CWC, TEMPERATURE, FREQUENCY,       &
                   &         NU, U, DU, P11, RC11, IPSD, DMA,         &
                   &         PH1, NAB, P, DP, NR, R, RN, BC, A, B, NABR)    

                   PHH(ISPI,:)   = P11     ! INTERGRATED PHASE FUNCTION
                   RC_TMP(ISPI,:)= RC11    ! VOLUME EXT/SCAT/ABS COEFFS
                   DDm(ISPI)     = DMA     ! MASS-MEAN-DIAMETER
         ENDIF
      ENDDO
                      
      DO J=1,3                             ! ADD all COEFFICIENTS
         RC_TOT(J) = RC_TMP(1,J)+RC_TMP(2,J)
      ENDDO

      DO ISPI=1,N
         W0(ISPI) = min( 1._r8, RC_TMP(ISPI,2)/RC_TOT(3) )
      ENDDO
     



end Subroutine get_beta_cloud

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Cloud_Extinction

! $Log$
! Revision 2.1  2003/01/31 18:36:25  jonathan
! new module for cloud extinction
!
