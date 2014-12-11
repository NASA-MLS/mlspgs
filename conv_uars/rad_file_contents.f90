! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

MODULE rad_file_contents

IMPLICIT NONE

TYPE lvl1_hdr_t
   sequence
   character (len=1) :: rec_type
   character (len=1) :: qualifier
   character (len=1) :: qualifier_pad(2)
   integer :: recordno
   integer :: write_time(2)
   character (len=4) :: lvl1_versionno
   integer :: nummmaf
   integer :: uars_day
   integer :: start_time(2)
   integer :: stop_time(2)
   integer :: mls_status_day
END TYPE

TYPE limb_hdr_t
   sequence
   character (len=1) :: rec_type
   character (len=1) :: band_bank
   character (len=2) :: band_bank_pad
   integer :: recordno
   integer :: mmaf_time(2)
   integer :: mmafno
END TYPE

TYPE limb_stat_t
   sequence
   real :: prd_temps(16)
   integer :: dgap_mmaf
   integer :: maneuver_stat
   integer :: mls_status
   integer(kind=2) :: wall_mmaf(90)
   integer(kind=2) :: window_red_refs(90)
   integer(kind=2) :: window_red_sz(90)
   character (len=1) :: mmif_stat(32)
   character (len=1) :: mmaf_stat
   logical(kind=1) :: hga_interfer
   character (len=1) :: hga_pad(2)
END TYPE

TYPE limb_rad_t
   sequence
   integer(kind=2) :: rad_l1(90,2,32)
END TYPE

TYPE limb_oa_t
   sequence
   logical(kind=1) :: ptg_fov_bo_diag_map(7)
   character (len=1) :: ptg_fov_pad
   character (len=12) :: oa_att_retrn
   character (len=12) :: oa_orb_retrn
   integer :: oa_ephem_status
   integer :: oa_limb_calc_status(32)
   integer :: oa_sat_att_status(32)
   integer :: oa_sat_orb_status(32)
   integer :: ptg_fov_bo_diag_mmif_fst
   integer :: ptg_fov_bo_diag_mmif_lst
   integer :: ptg_fov_bo_diag_mmif_num
   integer :: ref_mmif
   integer :: ref_solar_illum
   integer :: ref_time(2)
   integer :: sat_geod_status
   real :: earth_geod_rad(32)
   real :: grnw_sid_time
   real :: ptg_fov_azim_offset(32)
   real :: ptg_fov_azim_thm(2)
   real :: ptg_fov_bo_diag_azimdif
   real :: ptg_fov_bo_diag_elevdif
   real :: ptg_fov_bo_diag_mmaf(7)
   real :: ptg_fov_elev_offset(32)
   real :: ptg_fov_elev_thm(2)
   real :: ptg_inst2mac_elev(2)
   real :: ptg_limb_pt(3)
   real :: ref_earth_radius
   real :: ref_lat
   real :: ref_long
   real :: ref_solar_time
   real :: ref_solar_zen
   real :: rollrate_uars(32)
   real :: roll_uars(32)
   real :: sat_gcrad(32)
   real :: sat_geod_alt
   real :: sat_geod_lat
   real :: sat_long
   real :: sat_vel(3)
   real :: tngt_geod_alt(32)
   real :: tngt_geod_lat(32)
   real :: tngt_long(32)
   real :: trans_inst2eci(3,3)
   real :: ypr(3)
   real :: ypr_rate(3)

END TYPE

END MODULE rad_file_contents

! $Log$
