! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module CCRT_M

  implicit NONE
  private
  public :: CCRT

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains

  ! ----------------------------------------------------  SCRT_DN  -----
  subroutine CCRT ( T_SCRIPT, E_RFLTY, INCOPTDEPTH, TAU, RADIANCE, I, TSCAT, W0 )

    use MLSCommon, only: IP, RP

! "Scalar Radiative Transfer Down."  "Down" is an anachronism.

!{ Accumulate the incremental opacities multiplied by the differential
!  temperatures to get radiative transfer:
!  $I(\mathbf{x}) = \sum_{i=1}^{2N} \Delta B_i \tau_i$
!  In the code below, $\Delta B$ is called t\_script.

! inputs:

    real(rp), intent(in) :: t_script(:)    ! differential temperatures (K).
    real(rp), intent(in) :: e_rflty        ! Earth surface reflectivity (0--1).
    real(rp), intent(in) :: incoptdepth(:) ! layer incremental optical depth,
                                           ! this comes in as a positive number.
    real(rp), intent(in) :: tscat(:)       ! scattering source function
    real(rp), intent(in) :: w0(:)          ! single scattering albedo    

! outputs:

    real(rp), intent(out) :: tau(:)        ! transmission function
    real(rp), intent(out) :: radiance      ! radiance(K).
    integer(ip), intent(out) :: i          ! integration stop index

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

    do i = 2, half_path
      total_opacity = total_opacity - incoptdepth(i)
      tau(i) = exp(total_opacity)
      radiance = radiance + t_script(i) * tau(i) 
      if ( total_opacity < black_out ) then
        tau(i+1:n_path) = 0.0_rp
        return
      end if
      call update_rad ( radiance, t_script(i), tscat(i), w0(i) )
    end do

    tau(i) = e_rflty * tau(i-1)
    radiance = radiance + t_script(i) * tau(i)
    call update_rad ( radiance, t_script(i), tscat(i), w0(i) )

    do while ( total_opacity >= black_out .and. i < n_path )
      total_opacity = total_opacity - incoptdepth(i)
      i = i + 1
      tau(i) = e_rflty * exp(total_opacity)
      radiance = radiance + t_script(i) * tau(i)
      call update_rad ( radiance, t_script(i), tscat(i), w0(i) )
    end do

    tau(i+1:n_path) = 0.0_rp

  end subroutine CCRT

  subroutine Update_rad ( Radiance, Scalar, ScatSource, ScatAlb )
 
  use MLSCommon, only: Rk => Rp 

  real(rk), intent(in) :: Scalar
  real(rk), intent(in) :: ScatSource
  real(rk), intent(in) :: ScatAlb
  real(rk), intent(inout) :: Radiance

!----------------------------------------------------------------------

      Radiance = (1-ScatAlb)*Radiance + ScatAlb*ScatSource 
     
      !! this is still under construction ! 

  end subroutine Update_rad

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module CCRT_M
! $Log$
! Revision 2.1  2003/12/08 17:51:03  jonathan
! compute scattering in radiative transfer
!

