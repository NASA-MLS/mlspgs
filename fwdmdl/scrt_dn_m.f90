! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module SCRT_DN_M

! "Scalar Condensed Radiative Transfer Down."  "Down" is an anachronism.

  implicit NONE
  private
  public :: SCRT, SCRT_PFA, DSCRT_DX, DSCRT_DT

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains

  ! -------------------------------------------------------  SCRT  -----
  subroutine SCRT ( Tan_pt, T_SCRIPT, E_RFLTY, INCOPTDEPTH, TAU, RADIANCE, &
                  & INC_RAD_PATH, I_STOP )
    use MLSKinds, only: IP, RP

!{ Accumulate the incremental opacities multiplied by the differential
!  temperatures to get radiative transfer:
!  $I(\mathbf{x}) = \sum_{i=1}^{2N} \Delta B_i \tau_i$
!  In the code below, $\Delta B$ is called t_script.

! inputs:

    integer, intent(in) :: Tan_pt          ! Tangent point index in t_script
    real(rp), intent(in) :: t_script(:)    ! differential temperatures (K).
    real(rp), intent(in) :: e_rflty        ! Earth surface reflectivity (0--1).
    real(rp), intent(in) :: incoptdepth(:) ! layer incremental optical depth,
!                           this comes in as a positive number.
! outputs:

    real(rp), intent(out) :: tau(:)        ! transmission function
    real(rp), intent(out) :: inc_rad_path(:) ! incremental radiance along the
                                           ! path.  t_script * tau.
    real(rp), intent(out) :: radiance      ! radiance(K).
    integer(ip), intent(out) :: i_stop     ! integration stop index

! internals

    integer(ip) :: n_path
    real(rp) :: total_opacity
    real(rp), parameter :: black_out = -15.0_rp

! begin code

    n_path = size(t_script)
    tau(1) = 1.0_rp
    inc_rad_path(1) = t_script(1)
    total_opacity = 0.0_rp
    radiance = t_script(1)

!{ Compute $\tau_i$ for $2 \leq i \leq t$, where $t$ is given by half_path.
!  $\tau_i = \exp \left \{ - \sum_{j=2}^i \Delta \delta_{j \rightarrow j-1} \right \}$.
!  $\Delta \delta_{j \rightarrow j-1}$ is given by incoptdepth and
!  $- \sum_{j=2}^i \Delta \delta_{j \rightarrow j-1}$ is given by total_opacity.

    do i_stop = 2, tan_pt
      total_opacity = total_opacity - incoptdepth(i_stop)
      tau(i_stop) = exp(total_opacity)
      inc_rad_path(i_stop) = t_script(i_stop) * tau(i_stop)
      radiance = radiance + inc_rad_path(i_stop)
      if ( total_opacity < black_out ) then
        tau(i_stop+1:n_path) = 0.0_rp
        inc_rad_path(i_stop+1:n_path) = 0.0_rp
        return
      end if
    end do

!{ Account for earth reflectivity at the tangent to the surface:
!  $\tau_{2N - t + 1} = \Upsilon \tau_t$.

    tau(tan_pt+1) = e_rflty * tau(tan_pt)
    inc_rad_path(tan_pt+1) = t_script(tan_pt+1) * tau(tan_pt+1)
    radiance = radiance + inc_rad_path(tan_pt+1)

!{ Compute $\tau_i$ for $i > 2 N - t + 1$, where $t$ is given by half_path.\\
!  $\tau_i = \tau_{2N - t + 1} \exp \left \{ - \sum_{j=2N - t + 1}^i
!    \Delta \delta_{j-1 \rightarrow j} \right \}$.

! We don't reset total_opacity, so we compute e_rflty * exp(total_opacity)
! instead of tau(tan_pt) * exp(total_opacity).  i_stop is tan_pt + 1 here.

    do while ( total_opacity >= black_out .and. i_stop < n_path )
      total_opacity = total_opacity - incoptdepth(i_stop)
      i_stop = i_stop + 1
      tau(i_stop) = e_rflty * exp(total_opacity)
      inc_rad_path(i_stop) = t_script(i_stop) * tau(i_stop)
      radiance = radiance + inc_rad_path(i_stop)
    end do

    tau(i_stop+1:n_path) = 0.0_rp
    inc_rad_path(i_stop+1:n_path) = 0.0_rp

  end subroutine SCRT

!------------------------------------------------------  SCRT_PFA  -----

  ! Combine the Tau's for PFA and non-PFA models.
  ! We don't have to worry about the tangent point, because SCRT
  ! already did when computing the Tau's.

  subroutine SCRT_PFA ( Channel, Tau_LBL, Tau_PFA, T_Script_PFA, RadV )

    use MLSKinds, only: RP
    use Tau_M, only: Tau_T

    integer, intent(in) :: Channel              ! Which channel in Tau_PFA?
                                                ! Index in channels stru, not chan#
    type(tau_t), intent(in) :: Tau_LBL, Tau_PFA ! Tau structures
    real(rp), intent(in) :: T_Script_PFA(:,:)   ! Delta B, Path X Channels
    real(rp), intent(out) :: RadV(:)            ! Radiances at frequencies
                                                ! within the channel

    integer :: FRQ_I                            ! Frequency index
    integer :: N_Tau_PFA                        ! Tau_PFA%i_stop(channel)
    integer :: PATH_I                           ! Path index

    n_tau_PFA = tau_PFA%i_stop(channel)

    do frq_i = 1, size(radV)
      radV(frq_i) = 0.0
      do path_i = 1, min(tau_LBL%i_stop(frq_i), n_tau_PFA)
        radV(frq_i) = radV(frq_i) + &
          & t_script_pfa(path_i,channel) * tau_lbl%tau(path_i,frq_i) * &
          &                                tau_pfa%tau(path_i,channel)
      end do ! path_i
      ! Tau's after i_stop are 0.0.
    end do ! frq_i

  end subroutine SCRT_PFA

!------------------------------------------------------  DSCRT_DT  -----
! Compute the scalarized condensed radiative transfer derivatives,
! with derivatives of the differential Temperatures.

  subroutine DSCRT_DT ( TAN_PT, D_DELTA_DT, TAU, INC_RAD_PATH, DT_SCRIPT_DT, &
                      & I_START, I_END, DRAD_DT )
    use MLSKinds, only: IP, RP

!{ $\frac{\text{d}I(\mathbf{x})}{\text{d}x} = \sum_{i=1}^{2N} Q_i \tau_i$,
!  where $Q_i = \frac{\partial \Delta B_i}{\partial x} - \Delta B_i W_i$
!  and $W_i$ is given below. The only difference from {\tt DSCRT_DX} is
!  that $\frac{\partial \Delta B_i}{\partial x}$ is not assumed to be zero
!  here. This works for any $x$, not just temperature.  $\Delta B$ is
!  called T-script in some notes, so {\tt DT_SCRIPT_DT} is
!  $\frac{\partial \Delta B_i}{\partial x}$.

    integer, intent(in) :: Tan_pt           ! Tangent point index in inc_rad_path
    integer(ip), intent(in) :: i_start      ! where non-zeros on the path begin
    integer(ip), intent(in) :: i_end        ! where non-zeros on the path end
    real(rp), intent(in) :: d_delta_dt(:)   ! path opacity derivatives wrt sve
    real(rp), intent(in) :: tau(:)          ! path transmission
    real(rp), intent(in) :: inc_rad_path(:) ! incremental radiance along the
                                            ! path.  t_script * tau.
    real(rp), intent(in) :: dt_script_dt(:) ! derivatives of differential temps
    real(rp), intent(out) :: drad_dt        ! radiance derivative wrt temperature

    integer(ip) :: i, n_path
    real(rp) :: w

    w = 0.0_rp
    drad_dt = dt_script_dt(1)

    n_path = size(inc_rad_path)

!{ $-W_i = -\sum_{j=2}^i \frac{\partial \Delta \delta_{j \rightarrow j - 1}}
!                           {\partial x}$,
!  where the derivative is given by {\tt D_DELTA_DT}.
!
! $\frac{\text{d}I(\mathbf{x})}{\text{d}x} :=
!  \frac{\text{d}I(\mathbf{x})}{\text{d}x}
!  + \frac{\partial \Delta B_i}{\partial x} \mathcal{T}_i +
!  \Delta B_i \mathcal{T}_i W_i$, where $\Delta B_i \mathcal{T}_i$ is given
!  by {\tt INC_RAD_PATH(i)}. See {\tt wvs-024}.

    do i = max(2,i_start), min(i_end,tan_pt)
      w = w - d_delta_dt(i)
      drad_dt = drad_dt + dt_script_dt(i) * tau(i) + inc_rad_path(i) * w
    end do

    if ( i_end <= tan_pt ) return

! Tangent point or bounce off the earth's surface

    drad_dt = drad_dt + dt_script_dt(tan_pt+1) * tau(tan_pt+1)  &
            + inc_rad_path(tan_pt+1) * w

! Same as above, but don't add in the zero-emission layer at the 
! tangent point or earth's surface

    do i = tan_pt+2, min(i_end,n_path)
      w = w - d_delta_dt(i-1)
      drad_dt = drad_dt + dt_script_dt(i) * tau(i) + inc_rad_path(i) * w
    end do

  end subroutine DSCRT_DT

!------------------------------------------------------  DSCRT_DX  -----
! Compute the scalarized condensed radiative transfer derivatives,
! without derivatives of the differential Temperatures.

  subroutine DSCRT_DX ( TAN_PT, D_DELTA_DX, INC_RAD_PATH, I_START, I_END, DRAD_DX )
    use MLSKinds, only: IP, RP

!{ $\frac{\partial I(\mathbf{x})}{\partial x_k} = \sum_{i=1}^{2N} Q_i \tau_i$,
!  where $Q_i = \frac{\partial \Delta B_i}{\partial x_k} - \Delta B_i W_i$
!  and $W_i$ is given below. The only difference from {\tt DSCRT_DT} is
!  that $\frac{\partial \Delta B_i}{\partial x}$ is assumed to be zero here.

! Inputs

    integer, intent(in) :: Tan_pt         ! Tangent point index in inc_rad_path
    real(rp), intent(in) :: d_delta_dx(:) ! path opacity derivatives
    real(rp), intent(in) :: inc_rad_path(:) ! incremental radiance along the
                                          ! path.  t_script * tau.

    integer(ip), intent(in) :: i_start    ! where non-zeros on the path begin
    integer(ip), intent(in) :: i_end      ! where non-zeros on path end

! Output

    real(rp), intent(out) :: drad_dx      ! radiance derivative wrt x

! internals

    integer(ip) :: i, n_path
    real(rp) :: w

    drad_dx = 0.0_rp
    w = 0.0_rp
    n_path = size(inc_rad_path)

!{ $-W_i = -\sum_{j=2}^i \frac{\partial \Delta \delta_{j \rightarrow j - 1}}
!                           {\partial x_k}$,
!  where the derivative is given by d_delta_dx.

    do i = max(2,i_start), min(i_end,tan_pt)
      w = w - d_delta_dx(i)
      drad_dx = drad_dx + inc_rad_path(i) * w
    end do

    if ( i_end <= tan_pt ) return

! Tangent point or bounce off the earth's surface

    if ( i_start <= tan_pt ) drad_dx = drad_dx + inc_rad_path(tan_pt+1) * w

! Same as above, but don't add in the zero-emission layer at the 
! tangent point or earth's surface

    do i = max(i_start,tan_pt+2) , min(n_path,i_end)
      w = w - d_delta_dx(i-1)
      drad_dx = drad_dx + inc_rad_path(i) * w
    end do

  end subroutine DSCRT_DX

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module SCRT_DN_M
! $Log$
! Revision 2.21  2011/08/31 20:03:15  vsnyder
! TeXnicalities
!
! Revision 2.20  2011/07/29 01:48:13  vsnyder
! Remove tabs (Intel compiler complains)
!
! Revision 2.19  2011/07/08 20:51:22  yanovsky
! Minor comments changes
!
! Revision 2.18  2011/06/02 22:38:01  yanovsky
! Add opacity second derivatives D2_DELTA_DX2_qr to subroutine D2SCRT_DX2
! to allow for computations of analytical Hessians in logarithmic basis
!
! Revision 2.17  2010/08/19 02:00:02  vsnyder
! Polish up TeXnicalities
!
! Revision 2.16  2010/06/09 22:09:04  yanovsky
! d2scrt_dx2 subroutine
!
! Revision 2.15  2010/03/09 00:17:39  vsnyder
! Repair TeXnicalities
!
! Revision 2.14  2009/06/23 18:26:11  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.13  2009/01/16 23:40:03  vsnyder
! Correct case of starting point on path past the tangent point.
! add PRINT statement to not_used_here.
!
! Revision 2.12  2006/12/13 02:32:03  vsnyder
! Drag the tangent point around instead of assuming it's the middle one
!
! Revision 2.11  2005/11/21 22:57:27  vsnyder
! PFA derivatives stuff
!
! Revision 2.10  2005/11/01 23:02:21  vsnyder
! PFA Derivatives
!
! Revision 2.9  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.8  2004/10/06 21:17:36  vsnyder
! Make blacked-out Taus 1.0 instead of 0.0
!
! Revision 2.7  2004/09/04 01:49:46  vsnyder
! Handle earth reflectivity at tangent point correctly in derivatives.
!
! Revision 2.6  2003/06/03 23:58:19  vsnyder
! Correct a comment
!
! Revision 2.5  2002/10/10 19:50:50  vsnyder
! Mostly cosmetic changes, including adding LaTeX equations
!
! Revision 2.4  2002/10/09 22:47:58  vsnyder
! Don't set tau(i_stop) to zero at the end
!
! Revision 2.3  2002/10/08 17:08:06  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.2  2002/10/02 20:08:16  vsnyder
! Insert copyright notice, other cosmetic changes
!
! Revision 2.1  2002/05/13 05:15:40  zvi
! Adding missing intent on dummy variables
!
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model
!
! Revision 1.6.2.1  2001/09/10 10:02:32  zvi
! Cleanup..comp_path_entities_m.f90
!
! Revision 1.1  2000/05/04 18:12:06  vsnyder
! Initial conversion to Fortran 90
