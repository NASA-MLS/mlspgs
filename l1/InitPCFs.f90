! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
MODULE InitPCFs ! Init PCF data used by MLSL1 program
!=============================================================================

  USE SDPToolkit, ONLY: PGS_PC_GetReference, PGS_S_SUCCESS
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: L1PCF, PCFData_T, GetPCFParameters
  
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
  
! This data type is used to store User-defined Runtime Parameters and other
! information taken from the PCF.

   TYPE PCFData_T
      CHARACTER(LEN=27) :: StartUTC
      CHARACTER(LEN=27) :: EndUTC
      CHARACTER(LEN=80) :: OutputVersion
      CHARACTER(LEN=80) :: Cycle
      CHARACTER(LEN=132) :: PCF_filename
      CHARACTER(LEN=132) :: L1CF_filename
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

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here
END MODULE InitPCFs

! $Log$
! Revision 2.5  2005/06/23 18:41:35  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.4  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
! Revision 2.3  2002/07/17 15:25:42  perun
! Increase length of outputversion and cycle
!
! Revision 2.2  2001/03/12 16:15:12  perun
! Add OutputVersion and Cycle parameters
!
! Revision 2.1  2001/02/23 20:44:56  perun
! Version 0.5 commit
!
