module PHYSICS

! Physical constants

  use MLSCommon, only: R8
  implicit NONE
  public

  Real(r8), Parameter :: H = 6.62606891e-34_r8 ! Physics News, 1998
  Real(r8), Parameter :: K = 1.3806503e-23_r8  ! NIST, 1998
  Real(r8), Parameter :: H_OVER_K = H / K * 1.0e6_r8 ! \nu in MHz
! Real(r8), Parameter :: H_OVER_K = 4.7992377e-5_r8 ! using H and K above
! Real(r8), Parameter :: H_OVER_K = 4.7992157e-5_r8 ! Zvi's original

!---------------------------- RCS Ident Info -------------------------------
  private :: Id, ModuleName
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------

end module PHYSICS
! $Log$
