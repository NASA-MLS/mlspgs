! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!==============================================================================
MODULE MLSMessageModule         ! Basic messaging for the MLSPGS suite
!==============================================================================

  USE SDPToolkit

  IMPLICIT NONE
  
  PRIVATE :: Id,ModuleName
! ------------------------------- RCS Ident Info ------------------------------
CHARACTER (LEN=130) :: Id = &
     & "$Id$"
CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
! -----------------------------------------------------------------------------

! This module provides low level messaging for the MLSPGS suite.  The main
! routine is MLSMessage, which generates log messages as directed by the user.
! The MLSMessage routine logs a message using the SDPToolkit routine
! PGS_SMF_GenerateStatusReport.  This writes a string to the `LogReport'
! file (PCF# 10101) in the toolkit.  In the Toolkit `substitute' it just does
! a simple print.

! The user can also choose to log the messages to a seperate file when running
! under the toolkit.  This is setup by MLSMessageSetup and closed by
! MLSMessageClose.  The cataloging of such a file is left up to the calling
! code.

! ---------------------------------------------------------------------------

! Define some low level parameters.  These are used by the calling code to
! indicate the severity or otherwise of the messages.

INTEGER, PARAMETER :: MLSMSG_Debug=1
INTEGER, PARAMETER :: MLSMSG_Info=2
INTEGER, PARAMETER :: MLSMSG_Warning=3
INTEGER, PARAMETER :: MLSMSG_Error=4

PRIVATE :: SeverityNames
CHARACTER (LEN=*), DIMENSION(4), PARAMETER :: SeverityNames = &
     & (/"Debug  ","Info   ","Warning","Error  "/)

! This set of parameters are simple prefixes for common messages

CHARACTER (LEN=*), PARAMETER :: MLSMSG_Allocate = &
     & "Allocation failed: "
CHARACTER (LEN=*), PARAMETER :: MLSMSG_Fileopen = &
     & "Failed to open file: "
CHARACTER (LEN=*), PARAMETER :: MLSMSG_Keyword = &
     & "Unrecognized configuration file keyword: "
CHARACTER (LEN=*), PARAMETER :: MLSMSG_L1BRead = &
     & "Unable to read L1B data item: "
CHARACTER (LEN=*), PARAMETER :: MLSMSG_Duplicate = &
     & "There is already an entry with the name "

! This datatype describes the configuration of the messaging suite

PRIVATE MLSMSG_PrefixLen
INTEGER, PARAMETER :: MLSMSG_PrefixLen=32

TYPE MLSMessageConfig_T
   LOGICAL :: suppressDebugs
   INTEGER :: logFileUnit
   CHARACTER (LEN=MLSMSG_PrefixLen) :: prefix
END TYPE MLSMessageConfig_T

! This private variable describes this configuration

PRIVATE MLSMessageConfig
TYPE (MLSMessageConfig_T) :: MLSMessageConfig=MLSMessageConfig_T(.TRUE.,-1,"")

CONTAINS

  ! -------------------------------------------------------------------------

  ! This first routine is the main `messaging' code.

  SUBROUTINE MLSMessage(severity,moduleNameIn,message)

    ! Dummy arguments
    INTEGER, INTENT(IN) :: severity ! e.g. MLSMSG_Error
    CHARACTER (LEN=*), INTENT(IN) :: moduleNameIn ! Name of module (see below)
    CHARACTER (LEN=*), INTENT(IN) :: message ! Line of text

    ! Local variables
    INTEGER :: dummy
    CHARACTER (LEN=512) :: line ! Line to output, should be long enough

    ! Executable code

    ! The moduleNameIn is <dollar>RCSFile: <filename>,v <dollar>

    IF ((.NOT. MLSMessageConfig%suppressDebugs).OR. &
         & (severity/=MLSMSG_Debug)) THEN
       
       ! Assemble a full message line

       line=TRIM(SeverityNames(severity))// &
            & " ("//moduleNameIn(11:(LEN_TRIM(moduleNameIn)-8)) &
            &  //"): "//message

       ! Log the message using the toolkit routine

       dummy=PGS_SMF_GenerateStatusReport(TRIM(MLSMessageConfig%prefix)// &
            & TRIM(line))

       ! Now, if we're also logging to a file then write to that too.

       IF (MLSMessageConfig%logFileUnit /= -1) &
            & WRITE (UNIT=MLSMessageConfig%logFileUnit,FMT=*) TRIM(line)

    ENDIF

    ! Now if it's an error, then try to close log file if any and quit

    IF (severity==MLSMSG_Error) THEN
       IF (MLSMessageConfig%logFileUnit /= -1) &
            & CLOSE(MLSMessageConfig%logFileUnit)
       STOP
    ENDIF
  END SUBROUTINE MLSMessage

  ! ----------------------------------------------------------------------

  ! This routine sets up the MLSMessage suite.  The defaults are of course
  ! sensible, but the user may wish to chage things.

  SUBROUTINE MLSMessageSetup(suppressDebugs,logFileUnit,prefix)

    ! Dummy arguments
    LOGICAL, OPTIONAL, INTENT(IN) :: suppressDebugs
    INTEGER, OPTIONAL, INTENT(IN) :: logFileUnit
    CHARACTER (LEN=*), OPTIONAL, INTENT(IN) :: prefix

    ! Local variables
    INTEGER :: fileUnit, status

    ! Executable code

    IF (PRESENT(suppressDebugs)) &
         & MLSMessageConfig%suppressDebugs=suppressDebugs

    IF (PRESENT(prefix)) &
         & MLSMessageConfig%prefix=prefix

    IF (PRESENT(logFileUnit)) THEN
       IF (MLSMessageConfig%logFileUnit/=-1) CALL MLSMessage(MLSMSG_Error, &
            & ModuleName,"Already writing to a log file")
       MLSMessageConfig%logFileUnit=logFileUnit
    END IF
  END SUBROUTINE MLSMessageSetup

  ! ----------------------------------------------------------------------

  ! This routine simply closes the MLSMessage log file if there is one.

  SUBROUTINE MLSMessageClose
    ! Executable code
    IF (MLSMessageConfig%logFileUnit/=-1) &
         & MLSMessageConfig%logFileUnit=-1
  END SUBROUTINE MLSMessageClose

!===========================================================================
END MODULE MLSMessageModule
!===========================================================================

!
! $Log$
! Revision 2.1  2002/03/29 20:23:48  perun
! Version 1.0 commit
!
! Revision 1.8  1999/12/17 21:38:17  livesey
! Added MLSMSG_Duplicate
!
! Revision 1.7  1999/12/16 17:52:38  livesey
! Added MLSMSG_L1BRead
!
! Revision 1.6  1999/12/16 00:15:06  livesey
! Added MLSMSG_Keyword string constant.
!
! Revision 1.5  1999/12/10 18:24:41  nakamura
! Removed declaration for PGS_SMF_GenerateStatusReport (redundant with INTERFACE in SDPToolkit.f90).
!
! Revision 1.4  1999/12/08 17:59:53  nakamura
! Added function declaration for PGS_SMF_GenerateStatusReport.
!
! Revision 1.3  1999/12/02 23:47:03  livesey
! Added the file logging capability, but not properly tested yet.
!
! Revision 1.2  1999/12/02 01:56:11  livesey
! Changed filenames to mixed case, redid makefiles.
!
! Revision 1.1  1999/12/01 23:01:40  livesey
! Before renaming things to upper/lower case
!
!
