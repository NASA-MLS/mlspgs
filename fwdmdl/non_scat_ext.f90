! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Non_scat_ext

  use GasAbsorption, only: GET_BETA
  use MLSCommon, only: r8
  use SpectraLines, only: SETUP_SPECTRA
  use SpectroscopyCatalog_m, only: CATALOG_T
  

  IMPLICIT NONE
      Private
      Public :: get_beta_clear

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
  private :: not_used_here 
 !---------------------------------------------------------------------------

contains

      SUBROUTINE GET_BETA_CLEAR ( L, F, T, P, VMR, NU, NS, Ext_coeff )

!========================================================================

      INCLUDE 'spectra.f9h' 

      INTEGER, intent(in) :: NS                ! NO. OF CHEMICAL SPECIES
      INTEGER, intent(in) :: NU                ! 2 x NO. OF pts
      INTEGER, intent(in) :: L                 ! NO. OF (levels?)
      REAL(r8), intent(in) :: F                ! Frequency
      REAL(r8), intent(in) :: T(L)
      REAL(r8), intent(in) :: P(L)
!!     REAL(r8), intent(in) :: Z(L)
      REAL(r8), intent(in) :: VMR(NS,L)
      REAL(r8), intent(out) :: Ext_coeff(L)

! local variables

      REAL(r8) :: VMR_H2O
      REAL(r8) :: Beta
      REAL(r8) :: VMR1(NS-1)
      integer :: I
      integer :: J
!-------------------------------------------------------------------------

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

        CALL SETUP_SPECTRA(QLG,V0,GSE,IST,WTH,NTH,DELTA,N1,      &
                    &      GAMMA,N2,MOL,NMOL,NCNT)

      DO I=1,L

        VMR_H2O = VMR(1,I)
        do J=2,NS
          VMR1(J-1) = VMR(J,I)
        enddo 
!        print*, VMR_H2O,VMR1(1), VMR1(2), T(I), P(I)
        CALL GET_BETA(QLG,V0,GSE,IST,WTH,NTH,DELTA,N1,GAMMA,N2,  &
                    &        NMOL,NCNT,T(I),P(I),F,VMR_H2O,VMR1,Ext_coeff(I),NS )   
!        print*, Ext_Coeff(I), T(I), P(I), Z(I), VMR_H2O, VMR1(1), VMR1(2)

      ENDDO

!      stop

   end subroutine GET_BETA_CLEAR 

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))

  end function not_used_here

end module Non_scat_ext

! $Log$
! Revision 2.1  2003/11/17 17:50:52  jonathan
! initial work module for cloud model
!
