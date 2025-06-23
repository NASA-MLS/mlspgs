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

  use SDPToolkit, ONLY: PGS_S_SUCCESS
  use MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Info, &
    & MLSMessageConfig
  use MLSL1Common, ONLY : FileNameLen
  use Machine, ONLY : NeverCrash
  implicit none

  private

  public :: L1PCF, PCFData_T, GetPCFParameters

  ! -------------------------------------------------------------------
  ! The next parameter optionally sets level 1 to crash with a walkback
  ! if it logs a message containing a fatal string
  ! E.g., put the next line in your PCF
  ! 1005|CrashMsg|Drop. Dead.
  ! and your run will automatically crash at the point where it logs
  ! any message containing the string "Drop. Dead."
  
  integer, parameter :: mlspcf_l1_param_CrashMsg = 1005
  ! -------------------------------------------------------------------
  
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
      CHARACTER(LEN=FileNameLen) :: PCF_filename
      CHARACTER(LEN=FileNameLen) :: L1CF_filename
   END TYPE PCFData_T

   TYPE (PCFData_T) :: L1PCF

CONTAINS

!------------------------------------
   SUBROUTINE GetPCFParameters
!------------------------------------

   use MLSPCF1, only: MLSPcf_L1_Param_StartUTC, MLSPcf_L1_Param_EndUTC, &
     & MLSPcf_L1_Param_OutputVersion, MLSPcf_L1_Param_Cycle

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
     INTEGER :: msgStatus

     returnStatus = PGS_PC_GetConfigData (mlspcf_l1_param_startUTC, &
          L1PCF%startUTC)
     IF (returnStatus /= PGS_S_SUCCESS) THEN
        msgStatus=MLSMSG_Error
        msr = UDRP_ERR // 'startUTC'
     ELSE
        msgStatus=MLSMSG_Info
        msr = 'Start UTC : '//L1PCF%startUTC
     ENDIF 
     CALL MLSMessage (msgStatus, ModuleName, msr)
     
     returnStatus = PGS_PC_GetConfigData (mlspcf_l1_param_endUTC, L1PCF%endUTC)
     IF (returnStatus /= PGS_S_SUCCESS) THEN
        msr = UDRP_ERR // 'endUTC'
        msgStatus=MLSMSG_Error
     ELSE
        msgStatus=MLSMSG_Info
        msr = 'End UTC   : '//L1PCF%endUTC
     ENDIF
     CALL MLSMessage (msgStatus, ModuleName, msr)
     
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

     ! -----------------------------------------------------------------
     ! These next relate to message logging and to whether level 1
     ! will be peermitted to crash with walkback on any error
     ! -----------------------------------------------------------------
     ! Default settings (not yet configurable by the user)
     NeverCrash                       = .false.
     MLSMessageConfig%crashOnAnyError = .true.
     returnStatus = PGS_PC_GetConfigData ( mlspcf_l1_param_CrashMsg, msr )
     ! We'll allow you to either omit the PCFid or leave it blank
     ! to side-step this hair-trigger option
     IF ( returnStatus /= PGS_S_SUCCESS ) return
     if ( len_trim(msr) < 1 ) return
     NeverCrash = .false.
     MLSMessageConfig%CrashIfMsgSays = msr
     ! -----------------------------------------------------------------

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
! Revision 2.9  2023/06/06 22:36:04  pwagner
! May crash with walkback if appropriate
!
! Revision 2.8  2017/03/24 22:57:06  pwagner
! Made new PCFid 1005 that tells level 1 to crash with walkback if special msg logged
!
! Revision 2.7  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.6.6.1  2015/10/09 10:21:38  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.6  2006/04/05 18:11:02  perun
! Remove unused variables
!
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
