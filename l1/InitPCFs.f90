! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE InitPCFs ! Init PCF data used by MLSL1 program
!=============================================================================

  USE SDPToolkit, ONLY: PGS_PC_GetReference, PGS_S_SUCCESS
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------
  
! This data type is used to store User-defined Runtime Parameters and other
! information taken from the PCF.

   TYPE PCFData_T
      CHARACTER(LEN=27) :: StartUTC
      CHARACTER(LEN=27) :: EndUTC
      CHARACTER(LEN=5) :: OutputVersion
      CHARACTER(LEN=5) :: Cycle
   END TYPE PCFData_T

   TYPE (PCFData_T) :: L1PCF

CONTAINS

!------------------------------------
   SUBROUTINE GetPCFParameters
!------------------------------------

     USE MLSPCF1, ONLY: mlspcf_l1_param_StartUTC, mlspcf_l1_param_EndUTC, &
          mlspcf_l1_param_OutputVersion, mlspcf_l1_param_Cycle

! This subroutine retrieves the User-Defined Runtime Parameters from the PCF
! and stores them in variables used in the MLSL1 code.

! Parameters
     
     CHARACTER (LEN=*), PARAMETER :: UDRP_ERR = 'Failed to get runtime &
          &parameter from PCF:  '

! Functions
     
     INTEGER, EXTERNAL :: PGS_PC_GetConfigData

! Variables

     CHARACTER (LEN=480) :: msr
     INTEGER :: returnStatus
     
     returnStatus = PGS_PC_GetConfigData (mlspcf_l1_param_startUTC, &
          L1PCF%startUTC)
     IF (returnStatus /= PGS_S_SUCCESS) THEN
        msr = UDRP_ERR // 'startUTC'
        CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
     ENDIF
     
     returnStatus = PGS_PC_GetConfigData (mlspcf_l1_param_endUTC, L1PCF%endUTC)
     IF (returnStatus /= PGS_S_SUCCESS) THEN
        msr = UDRP_ERR // 'endUTC'
        CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
     ENDIF
     
     returnStatus = PGS_PC_GetConfigData (mlspcf_l1_param_OutputVersion, &
          L1PCF%OutputVersion)
     IF (returnStatus /= PGS_S_SUCCESS) THEN
        msr = UDRP_ERR // 'OutputVersion'
        CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
     ENDIF
     
     returnStatus = PGS_PC_GetConfigData (mlspcf_l1_param_Cycle, L1PCF%Cycle)
     IF (returnStatus /= PGS_S_SUCCESS) THEN
        msr = UDRP_ERR // 'Cycle'
        CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
     ENDIF

   END SUBROUTINE GetPCFParameters

END MODULE InitPCFs

! $Log$
! Revision 2.2  2001/03/12 16:15:12  perun
! Add OutputVersion and Cycle parameters
!
! Revision 2.1  2001/02/23 20:44:56  perun
! Version 0.5 commit
!
