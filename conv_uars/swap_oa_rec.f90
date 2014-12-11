! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Swap_OA_Rec_m

  implicit NONE
  private

  public :: Swap_OA_Rec

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains ! ============= Public procedures ===================================


  subroutine Swap_OA_Rec ( Limb_OA )

    USE Rad_File_Contents, only : limb_oa_t
    USE SwapEndian, only: SwapBig

    type(limb_oa_t), intent(inout) :: Limb_OA

    integer :: i, j

    limb_oa%oa_ephem_status = SwapBig (limb_oa%oa_ephem_status)
    do i = 1, 32
       limb_oa%oa_limb_calc_status(i) = SwapBig (limb_oa%oa_limb_calc_status(i))
       limb_oa%oa_sat_att_status(i) = SwapBig (limb_oa%oa_sat_att_status(i))
       limb_oa%oa_sat_orb_status(i) = SwapBig (limb_oa%oa_sat_orb_status(i))
       limb_oa%earth_geod_rad(i) = SwapBig (limb_oa%earth_geod_rad(i))

       limb_oa%ptg_fov_azim_offset(i) = SwapBig (limb_oa%ptg_fov_azim_offset(i))
       limb_oa%ptg_fov_elev_offset(i) = SwapBig (limb_oa%ptg_fov_elev_offset(i))
       limb_oa%rollrate_uars(i) = SwapBig (limb_oa%rollrate_uars(i))
       limb_oa%roll_uars(i) = SwapBig (limb_oa%roll_uars(i))
       limb_oa%sat_gcrad(i) = SwapBig (limb_oa%sat_gcrad(i))
       limb_oa%tngt_geod_alt(i) = SwapBig (limb_oa%tngt_geod_alt(i))
       limb_oa%tngt_geod_lat(i) = SwapBig (limb_oa%tngt_geod_lat(i))
       limb_oa%tngt_long(i) = SwapBig (limb_oa%tngt_long(i))
    enddo

  !  limb_oa%X = SwapBig (limb_oa%X)

    limb_oa%ptg_fov_bo_diag_mmif_fst = SwapBig (limb_oa%ptg_fov_bo_diag_mmif_fst)
    limb_oa%ptg_fov_bo_diag_mmif_lst = SwapBig (limb_oa%ptg_fov_bo_diag_mmif_lst)
    limb_oa%ptg_fov_bo_diag_mmif_num = SwapBig (limb_oa%ptg_fov_bo_diag_mmif_num)
    limb_oa%ref_mmif = SwapBig (limb_oa%ref_mmif)
    limb_oa%ref_solar_illum = SwapBig (limb_oa%ref_solar_illum)
    do i = 1, 2
       limb_oa%ref_time(i) = SwapBig (limb_oa%ref_time(i))
       limb_oa%ptg_fov_azim_thm(i) = SwapBig (limb_oa%ptg_fov_azim_thm(i))
       limb_oa%ptg_fov_elev_thm(i) = SwapBig (limb_oa%ptg_fov_elev_thm(i))
       limb_oa%ptg_inst2mac_elev(i) = SwapBig (limb_oa%ptg_inst2mac_elev(i))
    enddo
    limb_oa%sat_geod_status = SwapBig (limb_oa%sat_geod_status)
    limb_oa%grnw_sid_time = SwapBig (limb_oa%grnw_sid_time)
    limb_oa%ptg_fov_bo_diag_azimdif = SwapBig (limb_oa%ptg_fov_bo_diag_azimdif)
    limb_oa%ptg_fov_bo_diag_elevdif = SwapBig (limb_oa%ptg_fov_bo_diag_elevdif)
    do i = 1, 7
       limb_oa%ptg_fov_bo_diag_mmaf(i) = SwapBig (limb_oa%ptg_fov_bo_diag_mmaf(i))
    enddo
    do i = 1, 3
       limb_oa%ptg_limb_pt(i) = SwapBig (limb_oa%ptg_limb_pt(i))
       limb_oa%sat_vel(i) = SwapBig (limb_oa%sat_vel(i))
       limb_oa%ypr(i) = SwapBig (limb_oa%ypr(i))
       limb_oa%ypr_rate(i) = SwapBig (limb_oa%ypr_rate(i))
    enddo
    limb_oa%ref_earth_radius = SwapBig (limb_oa%ref_earth_radius)
    limb_oa%ref_lat = SwapBig (limb_oa%ref_lat)
    limb_oa%ref_long = SwapBig (limb_oa%ref_long)
    limb_oa%ref_solar_time = SwapBig (limb_oa%ref_solar_time)
    limb_oa%ref_solar_zen = SwapBig (limb_oa%ref_solar_zen)
    limb_oa%sat_geod_alt = SwapBig (limb_oa%sat_geod_alt)
    limb_oa%sat_geod_lat = SwapBig (limb_oa%sat_geod_lat)
    limb_oa%sat_long = SwapBig (limb_oa%sat_long)
    do i = 1, 3
       do j = 1, 3
          limb_oa%trans_inst2eci(i,j) = SwapBig(limb_oa%trans_inst2eci(i,j))
       enddo
    enddo

  end subroutine Swap_OA_Rec

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Swap_OA_Rec_m

! $Log$
