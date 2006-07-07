! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Compute_GL_Grid_M

  implicit NONE
  private
  public :: Compute_GL_Grid

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!-----------------------------------------------  Compute_GL_Grid  -----

  subroutine Compute_GL_Grid ( Z_PSIG, NLVL, MaxVert, P_GLgrid, Z_GLgrid )

  ! Compute the pressure and zeta GL grids.

    use Allocate_Deallocate, only: Allocate_Test
    use GLnp, only: GX, NGP1
    use MLSCommon, only: RP, R8

  ! Inputs:
    real(rp), dimension(:), intent(in) :: Z_psig  ! recommended PSIG for
                                       ! radiative transfer calculations
    integer, intent(in) :: NLVL                   ! size(z_psig)

  ! Outputs
    integer, intent(out) :: MaxVert               ! Levels in fine grid

  ! Would be intent(out) if they weren't pointers.
  ! First thing here is to nullify them.
    real(rp), dimension(:), pointer :: P_GLgrid   ! Pressure on glGrid surfs
    real(rp), dimension(:), pointer :: Z_GLgrid   ! Zeta on glGrid surfs

  ! Local variables
    integer :: NLM1                               ! NLVL - 1

    ! New Gauss points (excluding Lobatto end points) with -1 on the left:
    real(kind(gx)), parameter :: G_Grid(ngp1) = (/ -1.0_r8, gx /)

    nullify ( p_glgrid, z_glgrid )

    NLm1 = Nlvl - 1

! Allocate GL grid stuff

    maxVert = NLm1 * Ngp1 + 1
    call allocate_test ( z_glGrid, maxVert, 'z_glGrid', moduleName )
    call allocate_test ( p_glGrid, maxVert, 'p_glGrid', moduleName )

! From the selected integration grid pressures define the GL pressure grid:

    z_glgrid(1:maxVert-1) = reshape ( &
      ! Midpoint of integration grid intervals:
      & spread(0.5_rp * (z_psig(2:Nlvl) + z_psig(1:Nlm1)),1,Ngp1) + &
      ! Half length of integration grid intervals:
      & spread(0.5_rp * (z_psig(2:Nlvl) - z_psig(1:Nlm1)),1,Ngp1) * &
      ! New Gauss points (excludes Lobatto end points) with -1 at front:
      & spread(g_grid,2,NLm1), (/maxVert-1/))
    z_glgrid(maxVert) = z_psig(Nlvl)
    p_glgrid = 10.0_rp**(-z_glgrid)

  end subroutine Compute_GL_Grid

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Compute_GL_Grid_M

! $Log$
! Revision 2.14  2006/07/07 17:55:28  vsnyder
! Move some stuff to Compute_Z_PSIG_m
!
! Revision 2.13  2006/06/29 19:33:44  vsnyder
! Base grid calculations on interior (i.e., new) points in quadrature
! formula in case of Lobatto.
!
! Revision 2.12  2006/06/16 20:32:30  vsnyder
! Define NGP1 in glnp
!
! Revision 2.11  2005/08/25 00:49:01  vsnyder
! Don't look at the size of integration grid if it's not associated
!
! Revision 2.10  2005/08/09 15:15:43  pwagner
! Don't add pointer to z_all if not associated
!
! Revision 2.9  2005/06/22 18:08:18  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.8  2004/11/01 20:23:41  vsnyder
! Moved QtyStuff_t and associated dump routine to ForwardModelConfig
!
! Revision 2.7  2004/07/02 01:35:52  vsnyder
! Correct a comment
!
! Revision 2.6  2004/05/19 18:53:50  vsnyder
! Remove USE for unreferenced symbol
!
! Revision 2.5  2004/02/12 02:21:21  vsnyder
! Cosmetics
!
! Revision 2.4  2003/09/19 18:10:38  vsnyder
! Simplify computation of tangent point indices
!
! Revision 2.3  2003/06/20 19:35:17  pwagner
! Quanities now share grids stored separately in databses
!
! Revision 2.2  2003/05/20 00:06:23  vsnyder
! Remove stuff not used by FullForwardModel
!
! Revision 2.1  2003/05/15 20:25:23  vsnyder
! Initial commit
!
