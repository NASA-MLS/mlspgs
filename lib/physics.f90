! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module PHYSICS
! Physical constants
  use MLSCommon, only: R8, RP
  use Units, only: LN10
  implicit NONE
  public

  ! NIST as of 2004-04-01: http://physics.nist.gov/cuu/Constants/index.html
  real(r8), parameter :: Avogadro = 6.0221415e23_r8  ! #/mol ± 10e16 NIST 2004
  real(r8), parameter :: Bohr = 1.39962458_r8        ! MHz/Gauss, ± 12e-8 NIST 2004
  real(r8), parameter :: G_E = 2.0023193043718_r8    ! Unitless, ± 75e-13 NIST 2004
  real(rp), parameter :: MMM = 0.028964_rp           ! Kg/mol, Mean molecular mass of atmosphere kg/mol
  real(r8), parameter :: H = 6.62606931e-34_r8       ! J s, ± 11e-41 NIST 2004
  real(r8), parameter :: K = 1.3806505e-23_r8        ! J/K, ± 24e-20 NIST 2004
  real(r8), parameter :: H_OVER_K = H / K * 1.0e6_r8 ! \nu in MHz
! real(r8), parameter :: H_OVER_K = 4.79923752e-5_r8 ! using H and K above
! real(r8), parameter :: H_OVER_K = 4.7992157e-5_r8  ! Zvi's original
  real(rp), parameter :: Boltz = Avogadro * k * ln10 / mmm ! m^2/(K s^2),
! real(rp), parameter :: Boltz = 660.9853899_rp      ! Using above values
! real(rp), parameter :: Boltz = 660.988_rp          ! Zvi's original
  real(r8), parameter :: SpeedOfLight = 299792458.0_r8 ! M/s, exact NIST 2004

!---------------------------- RCS Ident Info -------------------------------
  private :: Id, IdParm, ModuleName
  character(len=*), parameter :: IdParm = "$Id$"
  character (len=256) :: Id = idParm
  character (len=*), parameter :: ModuleName= "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains 
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module PHYSICS
! $Log$
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
