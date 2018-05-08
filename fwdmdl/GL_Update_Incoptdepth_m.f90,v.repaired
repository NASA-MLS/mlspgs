! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module GL_Update_Incoptdepth_m

  implicit NONE
  private

  public :: GL_Update_Incoptdepth

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

  ! --------------------------------------  GL_Update_Incoptdepth  -----

  subroutine GL_Update_Incoptdepth ( gl_inds, more_inds, i_start, del_zeta, &
                                   & alpha_path, ds_dz_gw, ref_cor, incoptdepth )
    use GLNP, only: NG, NGP1
    use MLSCommon, only: RP, IP

    integer(ip), intent(in) :: gl_inds(:)    ! Gauss-Legendre grid indices
    integer(ip), intent(in) :: more_inds(:)  ! Places in the coarse path
                                             ! where GL is needed
    integer, intent(in) :: I_Start           ! Where in path to start integrating
    real(rp), intent(in) :: del_zeta(:)      ! Path -log(P) differences on the
                                             ! main grid.  This is for the whole
                                             ! coarse path, not just the part up
                                             ! to the black-out
    real(rp), intent(in), target :: alpha_path(:) ! Absorption coefficient on
                                             ! composite coarse & fine grid.
    real(rp), intent(in) :: ds_dz_gw(:)      ! Path length wrt zeta derivative * gw.
    real(rp), intent(in) :: ref_cor(:)       ! Refracted to unrefracted path
                                             ! length ratios.
    real(rp), intent(inout) :: incoptdepth(:) ! Incremental path opacities
                                             ! from one-sided layer calculation
                                             ! on output. It is the full
                                             ! integrated layer opacity.

  ! Internals

    integer :: A, AA, I, II
    real(rp), pointer :: Alpha_Path_c(:)

  ! see if anything needs to be gl-d

    if ( size(gl_inds) > 0 ) then

      !{ Apply Gauss-Legendre quadrature to the panels indicated by
      !  {\tt more\_inds}.  We remove a singularity (which actually only
      !  occurs at the tangent point) by writing
      !  $\int_{\zeta_i}^{\zeta_{i-1}} G(\zeta) \frac{\text{d}s}{\text{d}h}
      !   \frac{\text{d}h}{\text{d}\zeta} \text{d}\zeta =
      !  G(\zeta_i) \int_{\zeta_i}^{\zeta_{i-1}} \frac{\text{d}s}{\text{d}h}
      !   \frac{\text{d}h}{\text{d}\zeta} \text{d}\zeta +
      !  \int_{\zeta_i}^{\zeta_{i-1}} \left[ G(\zeta) - G(\zeta_i) \right]
      !   \frac{\text{d}s}{\text{d}h} \frac{\text{d}h}{\text{d}\zeta}
      !   \text{d}\zeta$.  The first integral is easy -- it's just
      !  $G(\zeta_i) (\zeta_{i-1}-\zeta_i)$.  Here, it is {\tt incoptdepth}.
      !  In the second integral, $G(\zeta)$ is {\tt alpha\_path\_gl} --
      !  which has already been evaluated at the appropriate abscissae -- and
      !  $G(\zeta_i)$ is {\tt alpha\_path\_c}.  The weights are {\tt gw}.

      alpha_path_c => alpha_path(1::ngp1)
      a = 1
      do i = 1, size(more_inds)
        aa = gl_inds(a)           ! First GL point in the panel
        ii = more_inds(i)         ! Index of coarse panel needing GL
        incoptdepth(ii) = incoptdepth(ii) + &
          & del_zeta(ii) * &
          & dot_product( (alpha_path(aa:aa+ng-1) - alpha_path_c(ii)), &
               & ds_dz_gw(aa:aa+ng-1) )
        a = a + ng
      end do ! i
    end if

    incoptdepth(i_start:) = ref_cor(i_start:) * incoptdepth(i_start:)

  end subroutine GL_Update_Incoptdepth

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module GL_Update_Incoptdepth_m

! $Log$
! Revision 2.1  2018/05/08 23:01:41  vsnyder
! Initial commit
!
