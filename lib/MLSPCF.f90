! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
! NOTE: This module is automatically created by the makemlspcfdatamodule
!       perl script.  Do *NOT* attempt to modify this by hand.
!
MODULE MLSPCF
   INTEGER, PARAMETER :: mlspcf_l1b_oa_start = 20000                 
   INTEGER, PARAMETER :: mlspcf_l1b_oa_end = 20039
   !
   INTEGER, PARAMETER :: mlspcf_l1b_rad_start = 20040
   INTEGER, PARAMETER :: mlspcf_l1b_rad_end = 20079
   !
   INTEGER, PARAMETER :: mlspcf_l1b_log_start = 20080
   INTEGER, PARAMETER :: mlspcf_l1b_log_end = 20080
   !
   INTEGER, PARAMETER :: mlspcf_l1b_eng_start = 20081
   INTEGER, PARAMETER :: mlspcf_l1b_eng_end = 20081
   !
   INTEGER, PARAMETER :: mlspcf_l2gp_start = 21000                 
   INTEGER, PARAMETER :: mlspcf_l2gp_end = 21050
   !
   INTEGER, PARAMETER :: mlspcf_MCF_start = 40260 
   INTEGER, PARAMETER :: mlspcf_MCF_end = 40276
   !
   INTEGER, PARAMETER :: mlspcf_l2aux_start = 21440
   INTEGER, PARAMETER :: mlspcf_l2aux_end = 21440
   !
   INTEGER, PARAMETER :: mlspcf_l2fwm_start = 21480
   INTEGER, PARAMETER :: mlspcf_l2fwm_end = 21480
   !
   INTEGER, PARAMETER :: mlspcf_l2ncep_start = 21500
   INTEGER, PARAMETER :: mlspcf_l2ncep_end = 21500

   INTEGER, PARAMETER :: mlspcf_l2dao_start = 21600
   INTEGER, PARAMETER :: mlspcf_l2dao_end = 21600

   INTEGER, PARAMETER :: mlspcf_l2clim_start = 21700
   INTEGER, PARAMETER :: mlspcf_l2clim_end = 21700

   INTEGER, PARAMETER :: mlspcf_l2cf_start = 21482
   INTEGER, PARAMETER :: mlspcf_l2cf_end = 21482
   !
   INTEGER, PARAMETER :: mlspcf_nomen_start = 23000
   INTEGER, PARAMETER :: mlspcf_nomen_end = 23000
   !
   INTEGER, PARAMETER :: mlspcf_l3_param_InputVersion = 3000
   INTEGER, PARAMETER :: mlspcf_l3_param_OutputVersion = 3001
   INTEGER, PARAMETER :: mlspcf_l3_param_Cycle = 3002
   INTEGER, PARAMETER :: mlspcf_l3_param_L2DayRange = 3003
   INTEGER, PARAMETER :: mlspcf_l3_param_MinDays = 3004
   INTEGER, PARAMETER :: mlspcf_l3_param_RangDays = 3005
   !
   INTEGER, PARAMETER :: mlspcf_l3cf_start = 30000
   INTEGER, PARAMETER :: mlspcf_l3cf_end = 30002
   !
   INTEGER, PARAMETER :: mlspcf_l3dm_start = 30100
   INTEGER, PARAMETER :: mlspcf_l3dm_end = 30199 
   !
   INTEGER, PARAMETER :: mlspcf_mcf_l3log_start = 50001
   INTEGER, PARAMETER :: mlspcf_mcf_l3log_end = 50001
   !
   INTEGER, PARAMETER :: mlspcf_mcf_l3dm_start = 50100
   INTEGER, PARAMETER :: mlspcf_mcf_l3dm_end = 50199
END MODULE MLSPCF

!$Log $
