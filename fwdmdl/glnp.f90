! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module GLNP

! **************  Gauss-Legendre 3 point formula ***************
  use MLSKinds, only: R8
  implicit NONE
  private
  public :: GW, GW_All, GX, GX_All, Legendre, Lobatto, NG, NGNEW, NGP1
  integer, parameter :: Ng = 3
  integer, parameter :: Legendre = 0           ! 0 for Legendre, -2 for Lobatto
  logical, parameter :: Lobatto = legendre < 0 ! Lobatto or Legendre?
  integer, parameter :: NGNEW = ng + legendre  ! Number of new points
  integer, parameter :: NGP1 = ngnew + 1       ! Used to allocate grids

! These are the 3-point-Gauss-Legendre abscissae (X-axis) values in (-1,1):

  real(r8), parameter :: Gx_all(Ng) = (/ & ! -sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)
     & -7.74596669241483377036e-1_r8,  0.00000000000000000000e+0_r8, &
     &  7.74596669241483377036e-1_r8 /)

! These are the corresponding 3-point-Gauss-Legendre Weights values:

  real(r8), parameter :: Gw_all(Ng) = (/ & ! 5.0/9.0, 8.0/9.0, 5.0/9.0
     &  5.55555555555555555556e-1_r8,  8.88888888888888888889e-1_r8, &
     &  5.55555555555555555556e-1_r8 /)

! Just the abscissae and weights needed for updating the path integral:
  real(r8), parameter :: Gx(ngnew) = gx_all(1-legendre/2:ng+legendre/2)
  real(r8), parameter :: Gw(ngnew) = gw_all(1-legendre/2:ng+legendre/2)

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

end module GLNP

! $Log$
! Revision 2.7  2009/06/23 18:26:11  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.6  2006/06/29 19:31:21  vsnyder
! Inching toward Lobatto
!
! Revision 2.5  2006/06/16 20:32:31  vsnyder
! Define NGP1 in glnp
!
! Revision 2.4  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.3  2003/09/17 19:59:56  vsnyder
! Correct a comment
!
! Revision 2.2  2002/10/08 17:08:04  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.1  2002/09/06 18:18:43  vsnyder
! Cosmetic changes
!
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model
!
! Revision 1.1  2001/06/21 13:07:09  zvi
! Speed enhancement MAJOR update
!
! Revision 1.1  2001/06/06 00:00:00  Z.Shippony
! Initial conversion to Fortran 90
