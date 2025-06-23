! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module PHYSICS
! Physical constants
  use MLSKinds, only: R8, RP
  use Constants, only: LN10
  implicit NONE
  public

  ! NIST as of 2004-04-01: http://physics.nist.gov/cuu/Constants/index.html
  real(r8), parameter :: Avogadro = 6.02214129e23_r8 ! #/mol ± 27e15 NIST 2010
  real(r8), parameter :: Bohr = 1.399624555_r8       ! MHz/Gauss, ± 31e-9 NIST 2010
                                                     ! Bohr magneton
  real(r8), parameter :: G_E = 2.00231930436153_r8   ! Unitless, ± 53e-14 NIST 2010
                                                     ! Electron g factor
  real(rp), parameter :: MMM = 0.028964_rp           ! Kg/mol, Mean molecular mass of atmosphere kg/mol
  real(r8), parameter :: H = 6.62606947e-34_r8       ! J s, ± 29e-42 NIST 2010
  real(r8), parameter :: K = 1.3806488e-23_r8        ! J/K, ± 13e-30 NIST 2010
  real(r8), parameter :: H_OVER_K = H / K * 1.0e6_r8 ! \nu in MHz
! real(r8), parameter :: H_OVER_K = 4.799243276e-5_r8 ! using H and K above
! real(r8), parameter :: H_OVER_K = 4.7992157e-5_r8  ! Zvi's original
  real(rp), parameter :: Boltz = Avogadro * k * ln10 / mmm ! m^2/(K s^2),
! real(rp), parameter :: Boltz = 660.9853899_rp      ! Using above values
! real(rp), parameter :: Boltz = 660.988_rp          ! Zvi's original
  real(r8), parameter :: SpeedOfLight = 299792458.0_r8 ! M/s, exact NIST 2010

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains 
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module PHYSICS
! $Log$
! Revision 2.8  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.7  2007/12/19 04:01:38  vsnyder
! Get R8 and RP from MLSKinds instead of MLSCommon
!
! Revision 2.6  2007/12/19 03:59:39  vsnyder
! Get LN10 from Constants instead of from Units
!
! Revision 2.5  2005/06/22 17:25:50  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.4  2004/04/02 18:23:36  vsnyder
! Changed a few comments
!
! Revision 2.3  2004/04/01 23:59:37  vsnyder
! Get latest values from NIST, add Bohr and G_e
!
! Revision 2.2  2004/03/20 04:05:02  vsnyder
! Add Avogadro and MMM.  Move SpeedOfLight from units.  Define Boltz in terms
! of Avogadro, MMM, k and ln10 (instead of just typing it in).
!
! Revision 2.1  2003/05/05 23:00:05  livesey
! Merged in feb03 newfwm branch
!
! Revision 1.1.2.1  2003/04/16 20:11:30  vsnyder
! Moved from ../fwdmdl
!
! Revision 2.1  2002/10/08 17:08:05  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model
!
! Revision 1.5  2001/06/07 23:39:31  pwagner
! Added Copyright statement
!
! Revision 1.4  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.1  2000/05/04 18:12:06  vsnyder
! Initial conversion to Fortran 90
!
