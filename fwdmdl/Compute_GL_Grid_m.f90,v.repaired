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

  subroutine Compute_GL_Grid ( Z_PSIG, Z_GLgrid, P_GLgrid )

  ! Compute the pressure and zeta GL grids.

    use GLnp, only: GX, NGP1
    use MLSCommon, only: RP, R8

  ! Inputs:
    real(rp), intent(in) :: Z_psig(:)  ! recommended PSIG for

  ! Outputs
    real(rp), dimension(:), intent(out) :: Z_GLgrid   ! Zeta on glGrid surfs
    real(rp), dimension(:), intent(out), optional :: P_GLgrid   ! Pressure on glGrid surfs

  ! Local variables
    integer :: MaxVert     ! Levels in fine grid
    integer :: NLVL        ! size(z_psig)
    integer :: NLM1        ! NLVL - 1

    ! New Gauss points (excluding Lobatto end points) with -1 on the left:
    real(kind(gx)), parameter :: G_Grid(ngp1) = (/ -1.0_r8, gx /)

    nlvl = size(z_psig)
    nlm1 = nlvl - 1
    ! Using "min" is necessary below in case there is no tangent point
    maxVert = min ( nlvl * ngp1, ubound(z_glgrid,1) )
! From the selected integration grid pressures define the GL pressure grid:

    z_glgrid(1:maxVert-1) = reshape ( &
      ! Midpoint of integration grid intervals:
      & spread(0.5_rp * (z_psig(2:Nlvl) + z_psig(1:Nlm1)),1,Ngp1) + &
      ! Half length of integration grid intervals:
      & spread(0.5_rp * (z_psig(2:Nlvl) - z_psig(1:Nlm1)),1,Ngp1) * &
      ! New Gauss points (excludes Lobatto end points) with -1 at front:
      & spread(g_grid,2,NLm1), (/maxVert-1/))
    if ( maxVert > 0 ) z_glgrid(maxVert) = z_psig(Nlvl)
    if ( present(p_glgrid) ) p_glgrid(:maxVert) = 10.0_rp**(-z_glgrid(:maxVert))

  end subroutine Compute_GL_Grid

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Compute_GL_Grid_M

! $Log$
! Revision 2.22  2013/05/17 22:54:30  vsnyder
! Revise nlvl computation
!
! Revision 2.21  2013/04/19 23:56:23  vsnyder
! Don't violate array bounds
!
! Revision 2.19  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.18  2007/02/01 02:45:18  vsnyder
! Exchange order of P_GLGrid, Z_GLGrid, make P_GLGrid optional
!
! Revision 2.17  2006/12/04 21:17:28  vsnyder
! Reorganize FullForwardModel to use automatic arrays instead of allocating
! pointer arrays.  Requires testing for zero size instead of testing for
! associated in several subsidiary procedures.
!
! Revision 2.16  2006/10/10 01:34:59  vsnyder
! Compute p_glgrid correctly in non-allocate case
!
! Revision 2.15  2006/09/20 01:39:09  vsnyder
! Add an optional 'allocate' argument to control whether to allocate the result
!
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
