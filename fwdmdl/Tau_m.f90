! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Tau_M

  use MLSCommon, only: RP

  implicit NONE
  private

  public :: Destroy_Tau, Dump, Dump_Tau, Get_Tau

  type, public :: Tau_T
    real(rp), pointer :: Tau(:,:) => NULL() ! Path X Frequencies
    integer, pointer :: I_Stop(:) => NULL() ! Amount of path to use, size = #Frequencies
  end type Tau_T

  interface Dump; module procedure Dump_Tau; end interface Dump

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

  ! ----------------------------------------------------  Get_Tau  -----
  subroutine Get_Tau ( Frq_i, gl_inds, more_inds, e_rflty, del_zeta, &
                     &  alpha_path_c, ref_cor, incoptdepth, &
                     &  alpha_path_gl, ds_dz_gw, tau )

  ! Update IncOptDepth with its GL corrections.  Multiply it by the
  ! refraction correction Ref_Cor.  Compute Tau = exp(indefinite sum of
  ! IncOptDepth).  The half-path point is a zero-thickness layer that
  ! doesn't have any IncOptDepth.  Instead, multiply tau(hal_path) by
  ! e_rflty to get tau(half_path+1).

    use GLNP, only: NG
    use MLSCommon, only: RP, IP

  ! inputs

    integer(ip), intent(in) :: Frq_I         ! Which frequency slice in Tau_T?
    integer(ip), intent(in) :: gl_inds(:)    ! Gauss-Legendre grid indices
    integer(ip), intent(in) :: more_inds(:)  ! Places in the coarse path
  !                                            where GL is needed
    real(rp), intent(in) :: e_rflty          ! earth reflectivity value (0--1).
    real(rp), intent(in) :: del_zeta(:)      ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    real(rp), intent(in) :: alpha_path_c(:)  ! absorption coefficient on coarse
  !                                            grid.
    real(rp), intent(in) :: ref_cor(:)       ! refracted to unrefracted path
  !                                            length ratios.
    real(rp), intent(inout) :: incoptdepth(:) ! incremental path opacities
  !                            from one-sided layer calculation on output.
  !                            it is the full integrated layer opacity.
    real(rp), intent(in) :: alpha_path_gl(:) ! absorption coefficient on gl
  !                                            grid.
    real(rp), intent(in) :: ds_dz_gw(:)      ! path length wrt zeta derivative * gw.

  ! outputs

    type(tau_t), intent(inout) :: tau        ! transmission function.  inout so
  !                            as not to clobber the association status of its
  !                            components.  Initial values not used.

  ! Internals

    integer :: A, AA, Half_Path, I, I_Stop, N_Path
    real(rp) :: total_opacity
    real(rp), parameter :: black_out = -15.0_rp

  ! Begin code

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

      a = 1
      do i = 1, size(more_inds)
        aa = gl_inds(a)
        incoptdepth(more_inds(i)) = incoptdepth(more_inds(i)) + &
          & del_zeta(more_inds(i)) * &
          & sum( (alpha_path_gl(a:a+ng-1) - alpha_path_c(more_inds(i))) * &            
               & ds_dz_gw(aa:aa+ng-1) ) 
        a = a + ng
      end do ! i

    end if

    incoptdepth = ref_cor * incoptdepth

  ! Compute Tau = exp(indefinite sum of IncOptDepth)

    n_path = size(incoptdepth)
    half_path = n_path / 2
    tau%tau(1,frq_i) = 1.0_rp
    total_opacity = 0.0_rp

!{ Compute $\tau_i$ for $2 \leq i \leq t$, where $t$ is given by half\_path.
!  $\tau_i = \exp \left \{ - \sum_{j=2}^i \Delta \delta_{j \rightarrow j-1} \right \}$.
!  $\Delta \delta_{j \rightarrow j-1}$ is given by incoptdepth and
!  $- \sum_{j=2}^i \Delta \delta_{j \rightarrow j-1}$ is given by total\_opacity.

    do i_stop = 2, half_path
      total_opacity = total_opacity - incoptdepth(i_stop)
      tau%tau(i_stop,frq_i) = exp(total_opacity)
      if ( total_opacity < black_out ) go to 99
    end do

!{ Account for earth reflectivity at the tangent to the surface:
!  $\tau_{2N - t + 1} = \Upsilon \tau_t$.

    tau%tau(half_path+1,frq_i) = e_rflty * tau%tau(half_path,frq_i)

!{ Compute $\tau_i$ for $i > 2 N - t + 1$, where $t$ is given by half\_path.\\
!  $\tau_i = \tau_{2N - t + 1} \exp \left \{ - \sum_{j=2N - t + 1}^i
!    \Delta \delta_{j-1 \rightarrow j} \right \}$.

! We don't reset total_opacity, so we compute e_rflty * exp(total_opacity)
! instead of tau(half_path) * exp(total_opacity).  i_stop is half_path + 1 here.

    do while ( total_opacity >= black_out .and. i_stop < n_path )
      total_opacity = total_opacity - incoptdepth(i_stop)
      i_stop = i_stop + 1
      tau%tau(i_stop,frq_i) = e_rflty * exp(total_opacity)
    end do

99  continue
    tau%tau(i_stop+1:n_path,frq_i) = 1.0_rp
    tau%i_stop(frq_i) = i_stop

  end subroutine Get_Tau

! --------------------------------------------------  Destroy_Tau  -----
  subroutine Destroy_Tau ( Tau, What, Where )
    use Allocate_Deallocate, only: Deallocate_test
    type(tau_t), intent(inout) :: Tau
    character(len=*), intent(in) :: What, Where
    call deallocate_test ( tau%tau, what//"%Tau", where )
    call deallocate_test ( tau%i_stop, what//"%I_stop", where )
  end subroutine Destroy_Tau

! -----------------------------------------------------  Dump_Tau  -----
  subroutine Dump_Tau ( Tau, What )
    use Dump_0, only: Dump
    use Output_m, only: Output

    type(tau_t), intent(in) :: Tau
    character(len=*), intent(in), optional :: What

    integer :: I

    if ( present(what) ) call output ( what, advance='yes' )

    do i = 1, size(tau%i_stop)
      call dump ( tau%tau(:tau%i_stop(i),i) )
    end do

  end subroutine Dump_Tau

  logical function NOT_USED_HERE()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function NOT_USED_HERE

end module Tau_M

! $Log$
! Revision 2.2  2005/03/03 02:08:32  vsnyder
! Add Destroy, Dump routines
!
! Revision 2.1  2004/10/07 17:34:15  vsnyder
! Initial commit
!
