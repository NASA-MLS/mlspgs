! Copyright 2006, by the California Institute of Technology. ALL
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

MODULE MLSPCF1

  IMPLICIT NONE

  PUBLIC

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

  INTEGER, PARAMETER :: mlspcf_l1_param_StartUTC = 1001
  INTEGER, PARAMETER :: mlspcf_l1_param_EndUTC = 1002
  INTEGER, PARAMETER :: mlspcf_l1_param_OutputVersion = 1003
  INTEGER, PARAMETER :: mlspcf_l1_param_Cycle = 1004
  INTEGER, PARAMETER :: mlspcf_l1_param_doinames = 1010
!
  INTEGER, PARAMETER :: mlspcf_pcf_start = 900
  INTEGER, PARAMETER :: mlspcf_pcf_end = 900
!
  INTEGER, PARAMETER :: mlspcf_APID1732_start = 21000
  INTEGER, PARAMETER :: mlspcf_APID1732_end = 21019
!
  INTEGER, PARAMETER :: mlspcf_APID1734_start = 21020
  INTEGER, PARAMETER :: mlspcf_APID1734_end = 21039
!
  INTEGER, PARAMETER :: mlspcf_APID1736_start = 21040
  INTEGER, PARAMETER :: mlspcf_APID1736_end = 21059
!
  INTEGER, PARAMETER :: mlspcf_APID1738_start = 21060
  INTEGER, PARAMETER :: mlspcf_APID1738_end = 21079
!
  INTEGER, PARAMETER :: mlspcf_APID1740_start = 21080
  INTEGER, PARAMETER :: mlspcf_APID1740_end = 21099
!
  INTEGER, PARAMETER :: mlspcf_APID1742_start = 21100
  INTEGER, PARAMETER :: mlspcf_APID1742_end = 21119
!
  INTEGER, PARAMETER :: mlspcf_APID1744_start = 21120
  INTEGER, PARAMETER :: mlspcf_APID1744_end = 21139
!
  INTEGER, PARAMETER :: mlspcf_APID1746_start = 21140
  INTEGER, PARAMETER :: mlspcf_APID1746_end = 21159
!
  INTEGER, PARAMETER :: mlspcf_APID1748_start = 21160
  INTEGER, PARAMETER :: mlspcf_APID1748_end = 21179
!
  INTEGER, PARAMETER :: mlspcf_l1cf_start = 901
  INTEGER, PARAMETER :: mlspcf_l1cf_end = 901
!
  INTEGER, PARAMETER :: mlspcf_nomen_start = 902
  INTEGER, PARAMETER :: mlspcf_nomen_end = 902
!
  INTEGER, PARAMETER :: mlspcf_engtbl_start = 903
  INTEGER, PARAMETER :: mlspcf_engtbl_end = 903
!
  INTEGER, PARAMETER :: mlspcf_defltgains_start = 904
  INTEGER, PARAMETER :: mlspcf_defltgains_end = 904
!
  INTEGER, PARAMETER :: mlspcf_defltzeros_start = 905
  INTEGER, PARAMETER :: mlspcf_defltzeros_end = 906
!
  INTEGER, PARAMETER :: mlspcf_dacsconst_start = 907
  INTEGER, PARAMETER :: mlspcf_dacsconst_end = 907
!
  INTEGER, PARAMETER :: mlspcf_sidebandfrac_start = 908
  INTEGER, PARAMETER :: mlspcf_sidebandfrac_end = 908
!
  INTEGER, PARAMETER :: mlspcf_spilloverloss_start = 909
  INTEGER, PARAMETER :: mlspcf_spilloverloss_end = 909
!
  INTEGER, PARAMETER :: mlspcf_defltchi2_start = 910
  INTEGER, PARAMETER :: mlspcf_defltchi2_end = 910
!
  INTEGER, PARAMETER :: mlspcf_defltbaselineAC_start = 911
  INTEGER, PARAMETER :: mlspcf_defltbaselineAC_end = 911
!
  INTEGER, PARAMETER :: mlspcf_bandalts_start = 912
  INTEGER, PARAMETER :: mlspcf_bandalts_end = 912
!
  INTEGER, PARAMETER :: mlspcf_bandsw_start = 913
  INTEGER, PARAMETER :: mlspcf_bandsw_end = 913
!
  INTEGER, PARAMETER :: mlspcf_strayrad_start = 914
  INTEGER, PARAMETER :: mlspcf_strayrad_end = 914
!
  INTEGER, PARAMETER :: mlspcf_sciMAF_start = 920
  INTEGER, PARAMETER :: mlspcf_sciMAF_end = 920
!
  INTEGER, PARAMETER :: mlspcf_engMAF_start = 921
  INTEGER, PARAMETER :: mlspcf_engMAF_end = 921
!
  INTEGER, PARAMETER :: mlspcf_MAF_data_start = 922
  INTEGER, PARAMETER :: mlspcf_MAF_data_end = 922
!
  INTEGER, PARAMETER :: mlspcf_l1b_radf_start = 30000
  INTEGER, PARAMETER :: mlspcf_l1b_radf_end = 30000
!
  INTEGER, PARAMETER :: mlspcf_l1b_radd_start = 30001
  INTEGER, PARAMETER :: mlspcf_l1b_radd_end = 30001
!
  INTEGER, PARAMETER :: mlspcf_l1b_oa_start = 30002
  INTEGER, PARAMETER :: mlspcf_l1b_oa_end = 30002
!
  INTEGER, PARAMETER :: mlspcf_l1b_eng_start = 30003
  INTEGER, PARAMETER :: mlspcf_l1b_eng_end = 30003
!
  INTEGER, PARAMETER :: mlspcf_l1b_diag_start = 30004
  INTEGER, PARAMETER :: mlspcf_l1b_diag_end = 30004
!
  INTEGER, PARAMETER :: mlspcf_l1b_radt_start = 30005
  INTEGER, PARAMETER :: mlspcf_l1b_radt_end = 30005
!
  INTEGER, PARAMETER :: mlspcf_l1b_log_start = 30006
  INTEGER, PARAMETER :: mlspcf_l1b_log_end = 30006
!
  INTEGER, PARAMETER :: mlspcf_l1b_diagT_start = 30007
  INTEGER, PARAMETER :: mlspcf_l1b_diagT_end = 30007
!
  INTEGER, PARAMETER :: mlspcf_mcf_l1log_start = 4000
  INTEGER, PARAMETER :: mlspcf_mcf_l1log_end = 4000
!
  INTEGER, PARAMETER :: mlspcf_mcf_l1boa_start = 4001
  INTEGER, PARAMETER :: mlspcf_mcf_l1boa_end = 4001
!
  INTEGER, PARAMETER :: mlspcf_mcf_l1bradf_start = 4002
  INTEGER, PARAMETER :: mlspcf_mcf_l1bradf_end = 4002
!
  INTEGER, PARAMETER :: mlspcf_mcf_l1bradd_start = 4003
  INTEGER, PARAMETER :: mlspcf_mcf_l1bradd_end = 4003
!
  INTEGER, PARAMETER :: mlspcf_mcf_l1bradt_start = 4004
  INTEGER, PARAMETER :: mlspcf_mcf_l1bradt_end = 4004

CONTAINS

  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here
END MODULE MLSPCF1

! $Log$
! Revision 2.14  2014/04/11 16:52:54  pwagner
! reserved PCFId for DOI
!
! Revision 2.13  2006/06/14 13:47:23  perun
! Add pcf numbers for stray radiance table
!
! Revision 2.12  2006/03/24 15:14:13  perun
! Add pcf numbers for BandAlts and BandSwitches tables
!
! Revision 2.11  2005/06/23 18:41:36  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.10  2004/11/10 15:36:34  perun
! Add pcf number for default baselineAC table file
!
! Revision 2.9  2004/08/12 13:51:50  perun
! Version 1.44 commit
!
! Revision 2.8  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
! Revision 2.7  2004/01/09 17:46:22  perun
! Version 1.4 commit
!
! Revision 2.6  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
