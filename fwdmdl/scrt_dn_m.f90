! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module SCRT_DN_M

  implicit NONE
  private
  public :: SCRT_DN, GET_DSCRT_NO_T_DN, GET_DSCRT_DN

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains

  ! ----------------------------------------------------  SCRT_DN  -----
  subroutine SCRT_DN ( T_SCRIPT, E_RFLTY, INCOPTDEPTH, TAU, RADIANCE, I_STOP )
    use MLSCommon, only: IP, RP

!{ Accumulate the incremental opacities multiplied by the differential
!  temperatures to get radiative transfer:
!  $I(\mathbf{x}) = \sum_{i=1}^{2N} \Delta B_i \tau_i$

! inputs:

    real(rp), intent(in) :: t_script(:)    ! differential temperatures (K).
    real(rp), intent(in) :: e_rflty        ! Earth surface reflectivity (0--1).
    real(rp), intent(in) :: incoptdepth(:) ! layer incremental optical depth,
!                           this comes in as a positive number.
! outputs:

    real(rp), intent(out) :: tau(:)        ! transmission function
    real(rp), intent(out) :: radiance      ! radiance(K).
    integer(ip), intent(out) :: i_stop     ! integration stop index

! internals

    integer(ip) :: n_path, half_path
    real(rp) :: total_opacity
    real(rp), parameter :: black_out = -15.0_rp

! begin code

    n_path = size(t_script)
    half_path = n_path / 2
    tau(1) = 1.0_rp
    total_opacity = 0.0_rp
    radiance = t_script(1)
    i_stop = 2

    do while ( total_opacity >= black_out .and. i_stop <= half_path )
      total_opacity = total_opacity - incoptdepth(i_stop)
      tau(i_stop) = exp(total_opacity)
      radiance = radiance + t_script(i_stop) * tau(i_stop)
      i_stop = i_stop + 1
    end do

    if ( i_stop <= half_path ) then
      tau(i_stop:n_path) = 0.0_rp
      i_stop = i_stop - 1
      return
    end if

! the tangent

    tau(i_stop) = e_rflty * tau(i_stop-1)
    radiance = radiance + t_script(i_stop) * tau(i_stop)

    do while ( total_opacity >= black_out .and. i_stop < n_path )
      total_opacity = total_opacity - incoptdepth(i_stop)
      i_stop = i_stop + 1
      tau(i_stop) = exp(total_opacity)
      radiance = radiance + t_script(i_stop) * tau(i_stop)
    end do

    tau(i_stop+1:n_path) = 0.0_rp

  end subroutine SCRT_DN
!-------------------------------------------  GET_DSCRT_NO_T_DN  -------
! Compute the scalarized condensed radiative transfer derivatives,
! without Temperature derivatives

  subroutine GET_DSCRT_NO_T_DN ( D_DELTA_DX, T_SCRIPT, TAU, I_START, &
    &                            I_END, DRAD_DX )
    use MLSCommon, only: IP, RP

! Inputs

    real(rp), intent(in) :: d_delta_dx(:) ! path opacity derivatives
    real(rp), intent(in) :: t_script(:) ! differential temparatures
    real(rp), intent(in) :: tau(:) ! path transmission path

    integer(ip), intent(in) :: i_start ! where non-zeros on the path begin
    integer(ip), intent(in) :: i_end   ! where non-zeros on path end

! Output

    real(rp), intent(out) :: drad_dx ! radiance derivative wrt x

! internals

    integer(ip) :: i, n_path, half_path
    real(rp) :: w

    w = 0.0_rp
    drad_dx = 0.0_rp
    n_path = size(t_script)
    half_path = n_path/2

    do i = max(2,i_start),min(i_end,half_path)
      w = w - d_delta_dx(i)
      drad_dx = drad_dx + t_script(i) * w * tau(i)
    end do

    if ( i_end <= half_path ) return

    drad_dx = drad_dx + t_script(half_path+1) * w * tau(half_path)

    do i = half_path+2 , min(n_path,i_end)
      w = w - d_delta_dx(i-1)
      drad_dx = drad_dx + t_script(i) * w * tau(i)
    end do

  end subroutine GET_DSCRT_NO_T_DN
!--------------------------------------------------  GET_DSCRT_DN  -----
! Compute the scalarized condensed radiative transfer derivatives

  subroutine GET_DSCRT_DN ( DDER, T_SCRIPT, TAU, DT_SCRIPT_DX, I_START, &
                  &         I_END, DRAD_DX)
    use MLSCommon, only: IP, RP

    integer(ip), intent(in) :: i_start, i_end
    real(rp), intent(in) :: DDER(:), T_SCRIPT(:), TAU(:), DT_SCRIPT_DX(:)
    real(rp), intent(out) :: drad_dx

    integer(ip) :: i, n_path, half_path
    real(rp) :: w

    w = 0.0_rp
    drad_dx = dt_script_dx(1)

    n_path = size(t_script)
    half_path = n_path/2

    do i = max(2,i_start),min(i_end,half_path)
      w = w - dder(i)
      drad_dx = drad_dx + (dt_script_dx(i) + t_script(i) * w) * tau(i)
    end do

    if ( i_end <= half_path ) return

    drad_dx = drad_dx + (dt_script_dx(half_path+1)  &
            + t_script(half_path+1) * w) * tau(half_path)

    do i = half_path+2 , min(n_path,i_end)
      w = w - dder(i-1)
      drad_dx = drad_dx + (dt_script_dx(i) + t_script(i) * w) * tau(i)
    end do

  end subroutine GET_DSCRT_DN

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module SCRT_DN_M
! $Log$
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
