! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module SpectraLines

! -------------------------------------------------------------------------  
! SET-UP SPECTRA LINE EMISSIONS
! -------------------------------------------------------------------------

      use MLSCommon, only: r8
      IMPLICIT NONE
      Private
      Public :: SETUP_SPECTRA

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
 !---------------------------------------------------------------------------
      
contains

      SUBROUTINE SETUP_SPECTRA(QLG,V0,GSE,IST,WTH,NTH,DELTA,N1, &
      &                         GAMMA,N2,MOL,NMOL,NCNT)
!-----------------------------------------------------------------------
      INCLUDE 'spectra.f9h'
      INCLUDE 'data.f9h'
      REAL(r8) ::  QQ(3)
      INTEGER :: I, J

! ----------------------------------------------------------------------

      DO I=1, N_O2_LINES
            QLG(1,I) = QQ_O2(1)
            QLG(2,I) = QQ_O2(2)
            QLG(3,I) = QQ_O2(3)
            V0(I)    = V0_O2(I)
            GSE(I)   = GSE_O2(I)
            IST(I)   = IST_O2(I)
            WTH(I)   = WTH_O2(I)
            NTH(I)   = NTH_O2(I)
            DELTA(I) = DELTA_O2(I)
            N1(I)    = N1_O2(I)
            GAMMA(I) = GAMMA_O2(I)
            N2(I)    = N2_O2(I)
      ENDDO

      J=N_O2_LINES
      DO I=J+1, J+N_H2O_LINES
            QLG(1,I)=QQ_H2O(1)
            QLG(2,I)=QQ_H2O(2)
            QLG(3,I)=QQ_H2O(3)
            V0(I)    = V0_H2O(I-J)
            GSE(I)   = GSE_H2O(I-J)
            IST(I)   = IST_H2O(I-J)
            WTH(I)   = WTH_H2O(I-J)
            NTH(I)   = NTH_H2O(I-J)
            DELTA(I) = DELTA_H2O(I-J)
            N1(I)    = N1_H2O(I-J)
            GAMMA(I) = GAMMA_H2O(I-J)
            N2(I)    = N2_H2O(I-J)
      ENDDO

      J=N_O2_LINES+N_H2O_LINES
      DO I=J+1, J+N_O18O_LINES
            QLG(1,I)=QQ_O_18_O(1)
            QLG(2,I)=QQ_O_18_O(2)
            QLG(3,I)=QQ_O_18_O(3)
            V0(I)    = V0_O_18_O(I-J)
            GSE(I)   = GSE_O_18_O(I-J)
            IST(I)   = IST_O_18_O(I-J)
            WTH(I)   = WTH_O_18_O(I-J)
            NTH(I)   = NTH_O_18_O(I-J)
            DELTA(I) = DELTA_O_18_O(I-J)
            N1(I)    = N1_O_18_O(I-J)
            GAMMA(I) = GAMMA_O_18_O(I-J)
            N2(I)    = N2_O_18_O(I-J)
      ENDDO

      J=N_O2_LINES+N_H2O_LINES+N_O18O_LINES
      DO I=J+1, J+N_H2O18_LINES
            QLG(1,I)=QQ_H2O_18(1)
            QLG(2,I)=QQ_H2O_18(2)
            QLG(3,I)=QQ_H2O_18(3)
            V0(I)    = V0_H2O_18(I-J)
            GSE(I)   = GSE_H2O_18(I-J)
            IST(I)   = IST_H2O_18(I-J)
            WTH(I)   = WTH_H2O_18(I-J)
            NTH(I)   = NTH_H2O_18(I-J)
            DELTA(I) = DELTA_H2O_18(I-J)
            N1(I)    = N1_H2O_18(I-J)
            GAMMA(I) = GAMMA_H2O_18(I-J)
            N2(I)    = N2_H2O_18(I-J)
      ENDDO

      J=N_O2_LINES+N_H2O_LINES+N_O18O_LINES+N_H2O18_LINES
      DO I=J+1, J+N_O3_LINES
            QLG(1,I)=QQ_O3(1)
            QLG(2,I)=QQ_O3(2)
            QLG(3,I)=QQ_O3(3)
      ENDDO

      J=N_O2_LINES+N_H2O_LINES+N_O18O_LINES+N_H2O18_LINES
      DO I=J+1, J+NL
            V0(I)    = V0_O3_PT1(I-J)
            GSE(I)   = GSE_O3_PT1(I-J)
            IST(I)   = IST_O3_PT1(I-J)
            WTH(I)   = WTH_O3_PT1(I-J)
            NTH(I)   = NTH_O3_PT1(I-J)
            DELTA(I) = DELTA_O3_PT1(I-J)
            N1(I)    = N1_O3_PT1(I-J)
            GAMMA(I) = GAMMA_O3_PT1(I-J)
            N2(I)    = N2_O3_PT1(I-J)
     
            V0(I+NL)    = V0_O3_PT2(I-J)
            GSE(I+NL)   = GSE_O3_PT2(I-J)
            IST(I+NL)   = IST_O3_PT2(I-J)
            WTH(I+NL)   = WTH_O3_PT2(I-J)
            NTH(I+NL)   = NTH_O3_PT2(I-J)
            DELTA(I+NL) = DELTA_O3_PT2(I-J)
            N1(I+NL)    = N1_O3_PT2(I-J)
            GAMMA(I+NL) = GAMMA_O3_PT2(I-J)
            N2(I+NL)    = N2_O3_PT2(I-J)

            V0(I+2*NL)    = V0_O3_PT3(I-J)
            GSE(I+2*NL)   = GSE_O3_PT3(I-J)
            IST(I+2*NL)   = IST_O3_PT3(I-J)
            WTH(I+2*NL)   = WTH_O3_PT3(I-J)
            NTH(I+2*NL)   = NTH_O3_PT3(I-J)
            DELTA(I+2*NL) = DELTA_O3_PT3(I-J)
            N1(I+2*NL)    = N1_O3_PT3(I-J)
            GAMMA(I+2*NL) = GAMMA_O3_PT3(I-J)
            N2(I+2*NL)    = N2_O3_PT3(I-J)

            V0(I+3*NL)    = V0_O3_PT4(I-J)
            GSE(I+3*NL)   = GSE_O3_PT4(I-J)
            IST(I+3*NL)   = IST_O3_PT4(I-J)
            WTH(I+3*NL)   = WTH_O3_PT4(I-J)
            NTH(I+3*NL)   = NTH_O3_PT4(I-J)
            DELTA(I+3*NL) = DELTA_O3_PT4(I-J)
            N1(I+3*NL)    = N1_O3_PT4(I-J)
            GAMMA(I+3*NL) = GAMMA_O3_PT4(I-J)
            N2(I+3*NL)    = N2_O3_PT4(I-J)

            V0(I+4*NL)    = V0_O3_PT5(I-J)
            GSE(I+4*NL)   = GSE_O3_PT5(I-J)
            IST(I+4*NL)   = IST_O3_PT5(I-J)
            WTH(I+4*NL)   = WTH_O3_PT5(I-J)
            NTH(I+4*NL)   = NTH_O3_PT5(I-J)
            DELTA(I+4*NL) = DELTA_O3_PT5(I-J)
            N1(I+4*NL)    = N1_O3_PT5(I-J)
            GAMMA(I+4*NL) = GAMMA_O3_PT5(I-J)
            N2(I+4*NL)    = N2_O3_PT5(I-J)

            V0(I+5*NL)    = V0_O3_PT6(I-J)
            GSE(I+5*NL)   = GSE_O3_PT6(I-J)
            IST(I+5*NL)   = IST_O3_PT6(I-J)
            WTH(I+5*NL)   = WTH_O3_PT6(I-J)
            NTH(I+5*NL)   = NTH_O3_PT6(I-J)
            DELTA(I+5*NL) = DELTA_O3_PT6(I-J)
            N1(I+5*NL)    = N1_O3_PT6(I-J)
            GAMMA(I+5*NL) = GAMMA_O3_PT6(I-J)
            N2(I+5*NL)    = N2_O3_PT6(I-J)

            V0(I+6*NL)    = V0_O3_PT7(I-J)
            GSE(I+6*NL)   = GSE_O3_PT7(I-J)
            IST(I+6*NL)   = IST_O3_PT7(I-J)
            WTH(I+6*NL)   = WTH_O3_PT7(I-J)
            NTH(I+6*NL)   = NTH_O3_PT7(I-J)
            DELTA(I+6*NL) = DELTA_O3_PT7(I-J)
            N1(I+6*NL)    = N1_O3_PT7(I-J)
            GAMMA(I+6*NL) = GAMMA_O3_PT7(I-J)
            N2(I+6*NL)    = N2_O3_PT7(I-J)
      ENDDO

      J=N_O2_LINES+N_H2O_LINES+N_O18O_LINES+N_H2O18_LINES+7*NL
      DO I=J+1, J+NLL
             V0(I)    = V0_O3_PT8(I-J)
            GSE(I)   = GSE_O3_PT8(I-J)
            IST(I)   = IST_O3_PT8(I-J)
            WTH(I)   = WTH_O3_PT8(I-J)
            NTH(I)   = NTH_O3_PT8(I-J)
            DELTA(I) = DELTA_O3_PT8(I-J)
            N1(I)    = N1_O3_PT8(I-J)
            GAMMA(I) = GAMMA_O3_PT8(I-J)
            N2(I)    = N2_O3_PT8(I-J)        
      ENDDO

      NCNT(1)=N_O2_LINES
      NCNT(2)=N_H2O_LINES
      NCNT(3)=N_O18O_LINES
      NCNT(4)=N_H2O18_LINES
      NCNT(5)=N_O3_LINES
      NMOL=5

      END SUBROUTINE SETUP_SPECTRA

end module SpectraLines

! $Log: SpectraLines.f90,v      

