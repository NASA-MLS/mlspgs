! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
! NOTE: This module is automatically created by the makemlspcfmodule
!       perl script.  Do *NOT* attempt to modify this by hand.
!
MODULE MLSPCF2
   INTEGER, PARAMETER :: mlspcf_l2_param_InputVersion = 2000
   INTEGER, PARAMETER :: mlspcf_l2_param_PGEVersion = 2001
   INTEGER, PARAMETER :: mlspcf_l2_param_Cycle = 2002
   INTEGER, PARAMETER :: mlspcf_l2_param_CCSDSStartId = 2003
   INTEGER, PARAMETER :: mlspcf_l2_param_CCSDSEndId = 2004
   INTEGER, PARAMETER :: mlspcf_l2_param_spec_keys = 2005
   INTEGER, PARAMETER :: mlspcf_l2_param_spec_hash = 2006
   !
   INTEGER, PARAMETER :: mlspcf_pcf_start = 900
   INTEGER, PARAMETER :: mlspcf_pcf_end = 900
   !
   INTEGER, PARAMETER :: mlspcf_l2cf_start = 901
   INTEGER, PARAMETER :: mlspcf_l2cf_end = 905
   !
   INTEGER, PARAMETER :: mlspcf_nomen_start = 20000
   INTEGER, PARAMETER :: mlspcf_nomen_end = 20000
   !
   INTEGER, PARAMETER :: mlspcf_l2ncep_start = 21000
   INTEGER, PARAMETER :: mlspcf_l2ncep_end = 21019
   !
   INTEGER, PARAMETER :: mlspcf_l2dao_start = 21020
   INTEGER, PARAMETER :: mlspcf_l2dao_end = 21039
   !
   INTEGER, PARAMETER :: mlspcf_l2clim_start = 21040
   INTEGER, PARAMETER :: mlspcf_l2clim_end = 21049
   !
   INTEGER, PARAMETER :: mlspcf_l1b_rad_start = 21050
   INTEGER, PARAMETER :: mlspcf_l1b_rad_end = 21109
   !
   INTEGER, PARAMETER :: mlspcf_l1b_oa_start = 21110
   INTEGER, PARAMETER :: mlspcf_l1b_oa_end = 21139
   !
   INTEGER, PARAMETER :: mlspcf_l2gp_start = 30000
   INTEGER, PARAMETER :: mlspcf_l2gp_end = 30569
   !
   INTEGER, PARAMETER :: mlspcf_l2dgg_start = 30570
   INTEGER, PARAMETER :: mlspcf_l2dgg_end = 30599
   !
   INTEGER, PARAMETER :: mlspcf_l2dgm_start = 30600
   INTEGER, PARAMETER :: mlspcf_l2dgm_end = 30629
   !
   INTEGER, PARAMETER :: mlspcf_l2fwm_full_start = 30630
   INTEGER, PARAMETER :: mlspcf_l2fwm_full_end = 30659
   !
   INTEGER, PARAMETER :: mlspcf_mcf_l2log_start = 4000
   INTEGER, PARAMETER :: mlspcf_mcf_l2log_end = 4000
   !
   INTEGER, PARAMETER :: mlspcf_mcf_l2gp_start = 4001
   INTEGER, PARAMETER :: mlspcf_mcf_l2gp_end = 4021
   !
   INTEGER, PARAMETER :: mlspcf_mcf_l2dgg_start = 4022
   INTEGER, PARAMETER :: mlspcf_mcf_l2dgg_end = 4022
   !
   INTEGER, PARAMETER :: mlspcf_mcf_l2dgm_start = 4023
   INTEGER, PARAMETER :: mlspcf_mcf_l2dgm_end = 4023
   !

   ! Change the following to 1 before delivering to sids;
   ! when set to 0, it allows program to run w/o creating metadata
   integer, parameter :: PENALTY_FOR_NO_METADATA = 0
END MODULE MLSPCF2
