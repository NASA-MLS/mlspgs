! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Non_scat_ext

  use GasAbsorption, only: GET_BETA
  use MLSCommon, only: r8
  use SpectraLines, only: SETUP_SPECTRA
  

  IMPLICIT NONE
      Private
      Public :: get_beta_clear

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

      SUBROUTINE GET_BETA_CLEAR ( L, F, T, P, VMR_H2O, VMR_O3, VmR_N2O, &
        &                         NU, Ext_coeff )

!========================================================================

      INCLUDE 'spectra.f9h' 

      INTEGER, intent(in) :: NU                ! 2 x NO. OF pts
      INTEGER, intent(in) :: L                 ! NO. OF (levels?)
      REAL(r8), intent(in) :: F                ! Frequency
      REAL(r8), intent(in) :: T(L)
      REAL(r8), intent(in) :: P(L)
!!     REAL(r8), intent(in) :: Z(L)
      REAL(r8), intent(in) :: VMR_H2O(:), VMR_O3(:), VMR_N2O(:)
      REAL(r8), intent(out) :: Ext_coeff(L)

! local variables

      integer :: I
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

        CALL GET_BETA(QLG,V0,GSE,IST,WTH,NTH,DELTA,N1,GAMMA,N2,NMOL,NCNT,  &
          &           T(I),P(I),F,VMR_H2O(I),VMR_O3(I),VMR_N2O(I),Ext_coeff(I) )

      ENDDO

!      stop

   end subroutine GET_BETA_CLEAR 

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Non_scat_ext

! $Log$
! Revision 2.7  2013/06/12 02:32:47  vsnyder
! Cruft removal
!
! Revision 2.6  2010/06/07 23:21:30  vsnyder
! Changed some calling sequences
!
! Revision 2.5  2009/06/23 18:26:11  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.4  2008/05/02 23:42:34  vsnyder
! Simplify, remove array copies
!
! Revision 2.3  2007/07/25 20:21:33  vsnyder
! Delete USE for unreferenced entities and declarations for unused variables
!
! Revision 2.2  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.1  2003/11/17 17:50:52  jonathan
! initial work module for cloud model
!
