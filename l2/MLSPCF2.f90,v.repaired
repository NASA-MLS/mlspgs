! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.
! NOTE: This module is automatically created by the makemlspcfmodule
!       perl script.  Do *NOT* attempt to modify this by hand.
!

MODULE MLSPCF2
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

   INTEGER, PARAMETER :: mlspcf_l2_param_InputVersion = 2000
   INTEGER, PARAMETER :: mlspcf_l2_param_PGEVersion = 2001
   INTEGER, PARAMETER :: mlspcf_l2_param_Cycle = 2002
   INTEGER, PARAMETER :: mlspcf_l2_param_CCSDSStartId = 2003
   INTEGER, PARAMETER :: mlspcf_l2_param_CCSDSEndId = 2004
   INTEGER, PARAMETER :: mlspcf_l2_param_spec_keys = 2005
   INTEGER, PARAMETER :: mlspcf_l2_param_spec_mcfnames = 2006
   INTEGER, PARAMETER :: mlspcf_l2_param_switches = 2007
   INTEGER, PARAMETER :: mlspcf_l2_param_col_spec_keys = 2008
   INTEGER, PARAMETER :: mlspcf_l2_param_col_spec_mcfnames = 2009
   INTEGER, PARAMETER :: mlspcf_l2_param_col_spec_doinames = 2010
   INTEGER, PARAMETER :: mlspcf_l2_param_CrashMsg = 2011
   !
   INTEGER, PARAMETER :: mlspcf_pcf_start = 900
   INTEGER, PARAMETER :: mlspcf_pcf_end = 900
   !
   INTEGER, PARAMETER :: mlspcf_mobile_l2cf_start = 901
   INTEGER, PARAMETER :: mlspcf_mobile_l2cf_end = 905
   !
   INTEGER, PARAMETER :: mlspcf_nomen_start = 20000
   INTEGER, PARAMETER :: mlspcf_nomen_end = 20000
   !
   INTEGER, PARAMETER :: mlspcf_antpats_start = 20001
   INTEGER, PARAMETER :: mlspcf_antpats_end = 20001
   !
   INTEGER, PARAMETER :: mlspcf_filtshps_start = 20002
   INTEGER, PARAMETER :: mlspcf_filtshps_end = 20002
   !
   INTEGER, PARAMETER :: mlspcf_dacsfltsh_start = 20003
   INTEGER, PARAMETER :: mlspcf_dacsfltsh_end = 20003
   !
   INTEGER, PARAMETER :: mlspcf_ptggrids_start = 20004
   INTEGER, PARAMETER :: mlspcf_ptggrids_end = 20004
   !
   INTEGER, PARAMETER :: mlspcf_pfa_start = 20005
   INTEGER, PARAMETER :: mlspcf_pfa_end = 20024
   !
   INTEGER, PARAMETER :: mlspcf_IGRF_start = 20025
   INTEGER, PARAMETER :: mlspcf_IGRF_end = 20025
   !
   INTEGER, PARAMETER :: mlspcf_ElevOffset_start = 20026
   INTEGER, PARAMETER :: mlspcf_ElevOffset_end = 20026
   !
   INTEGER, PARAMETER :: mlspcf_IsoTopeRatio_start = 20027
   INTEGER, PARAMETER :: mlspcf_IsoTopeRatio_end = 20027
   !
   INTEGER, PARAMETER :: mlspcf_SBFraction_start = 20028
   INTEGER, PARAMETER :: mlspcf_SBFraction_end = 20028
   !
   INTEGER, PARAMETER :: mlspcf_SysTemperature_start = 20029
   INTEGER, PARAMETER :: mlspcf_SysTemperature_end = 20029
   !
   INTEGER, PARAMETER :: mlspcf_NoiseBW_start = 20030
   INTEGER, PARAMETER :: mlspcf_NoiseBW_end = 20030
   !
   INTEGER, PARAMETER :: mlspcf_Misc_start = 20031
   INTEGER, PARAMETER :: mlspcf_Misc_end = 20031
   !
   INTEGER, PARAMETER :: mlspcf_MieTables_start = 20032
   INTEGER, PARAMETER :: mlspcf_MieTables_end = 20032
   !
   INTEGER, PARAMETER :: mlspcf_spectroscopy_start = 20033
   INTEGER, PARAMETER :: mlspcf_spectroscopy_end = 20033
   !
   INTEGER, PARAMETER :: mlspcf_surfaceHeight_start = 20034
   INTEGER, PARAMETER :: mlspcf_surfaceHeight_end = 20034
   !
   INTEGER, PARAMETER :: mlspcf_l2pc_start = 20035
   INTEGER, PARAMETER :: mlspcf_l2pc_end = 20099
   !
   INTEGER, PARAMETER :: mlspcf_l2apriori_start = 21000
   INTEGER, PARAMETER :: mlspcf_l2apriori_end = 21049
   !
   INTEGER, PARAMETER :: mlspcf_mobile_l1b_rad_start = 21050
   INTEGER, PARAMETER :: mlspcf_mobile_l1b_rad_end = 21109
   !
   INTEGER, PARAMETER :: mlspcf_mobile_l1b_oa_start = 21110
   INTEGER, PARAMETER :: mlspcf_mobile_l1b_oa_end = 21139
   !
   INTEGER, PARAMETER :: mlspcf_polygon_start = 21150
   INTEGER, PARAMETER :: mlspcf_polygon_end = 21150
   !
   INTEGER, PARAMETER :: mlspcf_l2ncep_start = 21900
   INTEGER, PARAMETER :: mlspcf_l2ncep_end = 21999
   !
   INTEGER, PARAMETER :: mlspcf_l2geos5_start = 22000
   INTEGER, PARAMETER :: mlspcf_l2geos5_end = 22099
   !
   INTEGER, PARAMETER :: mlspcf_l2dao_start = 22100
   INTEGER, PARAMETER :: mlspcf_l2dao_end = 22199
   !
   INTEGER, PARAMETER :: mlspcf_l2clim_start = 22200
   INTEGER, PARAMETER :: mlspcf_l2clim_end = 22249
   !
   INTEGER, PARAMETER :: mlspcf_l2neurnet_start = 22250
   INTEGER, PARAMETER :: mlspcf_l2neurnet_end = 22299
   !
   INTEGER, PARAMETER :: mlspcf_l2ascii_start = 22300
   INTEGER, PARAMETER :: mlspcf_l2ascii_end = 22399
   !
   INTEGER, PARAMETER :: mlspcf_mobile_l2gp_start = 30000
   INTEGER, PARAMETER :: mlspcf_mobile_l2gp_end = 30569
   !
   INTEGER, PARAMETER :: mlspcf_mobile_l2dgg_start = 30570
   INTEGER, PARAMETER :: mlspcf_mobile_l2dgg_end = 30599
   !
   INTEGER, PARAMETER :: mlspcf_mobile_l2dgm_start = 30600
   INTEGER, PARAMETER :: mlspcf_mobile_l2dgm_end = 30629
   !
   INTEGER, PARAMETER :: mlspcf_mobile_l2fwm_full_start = 30630
   INTEGER, PARAMETER :: mlspcf_mobile_l2fwm_full_end = 30659
   !
   INTEGER, PARAMETER :: mlspcf_mobile_mcf_l2log_start = 4000
   INTEGER, PARAMETER :: mlspcf_mobile_mcf_l2log_end = 4000
   !
   INTEGER, PARAMETER :: mlspcf_mobile_mcf_l2gp_start = 4001
   INTEGER, PARAMETER :: mlspcf_mobile_mcf_l2gp_end = 4021
   !
   INTEGER, PARAMETER :: mlspcf_mobile_mcf_l2dgg_start = 4022
   INTEGER, PARAMETER :: mlspcf_mobile_mcf_l2dgg_end = 4022
   !
   INTEGER, PARAMETER :: mlspcf_mobile_mcf_l2dgm_start = 4023
   INTEGER, PARAMETER :: mlspcf_mobile_mcf_l2dgm_end = 4023
   !
   include 'pcf2mob.f9h'
contains 
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

END MODULE MLSPCF2

! $Log$
! Revision 2.30  2022/03/16 22:20:55  pwagner
! Remove conflict between mlspcf_l2_param_CrashMsg and mlspcf_l2_param_col_spec_keys
!
! Revision 2.29  2020/12/22 22:21:01  pwagner
! Reserved id nums for Neural Net Coefficient Files
!
! Revision 2.28  2017/07/10 18:53:23  pwagner
! Added mlspcf_l2ascii_start, mlspcf_l2ascii_end
!
! Revision 2.27  2016/01/29 00:48:54  pwagner
! Added Polygon file
!
! Revision 2.26  2014/03/26 17:45:49  pwagner
! Added ProductionLocation, identifier_product_DOI to attributes
!
! Revision 2.25  2012/04/05 20:13:01  pwagner
! Added ids for igrf and other parameter files
!
! Revision 2.24  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.23  2008/12/18 22:00:01  pwagner
! Added print to not_used_here
!
! Revision 2.22  2008/05/08 19:18:22  pwagner
! added Mie Tables entries
!
! Revision 2.21  2008/01/08 00:18:48  pwagner
! Levels 1 and 2 can use same shared PCF now
!
! Revision 2.20  2007/01/10 00:02:54  pwagner
! Made PCF entries for spectroscopy and surfaceHeight files
!
! Revision 2.19  2006/06/15 17:35:32  pwagner
! Increased max num of l2pc from 30 to 65
!
! Revision 2.18  2006/06/13 18:19:08  pwagner
! Added pcfids for geos5; moved ncep pcfids backwards to make room
!
! Revision 2.17  2006/01/19 00:31:02  pwagner
! Set aside range for pfa files, col_spec_keys and _hash
!
! Revision 2.16  2005/08/25 20:30:20  pwagner
! Restored larger mlspcf_l2gp_end--better way to limit warnings found
!
! Revision 2.15  2005/06/22 18:57:02  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.14  2003/08/01 20:05:09  pwagner
! Lowered mlspcf_l2gp_end to limit toolkit warnings
!
! Revision 2.13  2003/07/16 21:50:53  pwagner
! Added mlspcf_dacsfltsh_start, _end
!
! Revision 2.12  2003/05/29 17:54:28  pwagner
! Able to read dao, ncep files w/o knowing name fragment
!
! Revision 2.11  2002/10/08 17:36:22  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.10  2002/10/01 22:25:16  pwagner
! Fixed RCS Ident Block
!
! Revision 2.9  2002/10/01 20:19:06  bwknosp
! Added Id, RCS, and Log info
!
