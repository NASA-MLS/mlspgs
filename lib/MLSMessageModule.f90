! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!==============================================================================
module MLSMessageModule         ! Basic messaging for the MLSPGS suite
!==============================================================================

  use HDF, only: DFACC_CREATE, DFACC_RDONLY, DFACC_RDWR
  use MACHINE, only: CRASH_BURN, EXIT_WITH_STATUS, NEVERCRASH
  use MLSCommon, only: MLSFile_T
  use MLSStrings, only: Capitalize
  use PVM, only: InfoTag, &
    & PVMDATADEFAULT, PVMFInitSend, PVMF90Pack, SIG_AboutToDie
  use SDPToolkit, only: PGS_S_SUCCESS, &
    & PGS_SMF_MASK_LEV_N, PGS_SMF_MASK_LEV_E, PGS_SMF_MASK_LEV_F, &
    & PGS_SMF_MASK_LEV_W, PGS_SMF_MASK_LEV_M, PGS_SMF_MASK_LEV_S, &
    & UseSDPToolkit, &
    & PGS_SMF_GenerateStatusReport, PGS_SMF_TestStatusLevel

  implicit none
  private

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! This module provides low level messaging for the MLSPGS suite.  The main
  ! routine is MLSMessage, which generates log messages as directed by the
  ! user. The MLSMessage routine logs a message using the SDPToolkit routine
  ! PGS_SMF_GenerateStatusReport.  This writes a string to the `LogReport'
  ! file (PCF# 10101) in the toolkit.  In the Toolkit `substitute' it just
  ! does a simple print.
  
  ! Alternate entries for special circumstances are PVMErrorMessage
  ! (to log a PVM Error)
  ! and ReportTKStatus
  
  ! Another choice to report an error is StopWithErrorMsg
  ! which lets you dump a calling stack you create using MLSMessageCalls
  ! (to report on the severity of a PGS Toolkit return status)
  
  ! Yet another mode is MLSMessage_, useful if you overload MLSMessage with
  ! another module's MLSMessage subroutine accepting extra args
  ! which can then turn around and call MLSMessage_

  ! The user can also choose to log the messages to a seperate file when
  ! running under the toolkit.  This is setup by MLSMessageSetup and closed
  ! by MLSMessageClose.  The cataloging of such a file is left up to the
  ! calling code.

  ! A lighter-weight substitute is MLSMessageSubstitute.f90 which dispenses
  ! with most of the toolkit panoply, needing only the modules
  ! Machine
  ! MLSKinds
  ! MLSCommon

  include 'MLSMessage.f9h'

  ! --------------------------------------------  ReportTKStatus  -----

  ! Report on severity represented by value returned by toolkit function
  ! Some functions may return non-zero values that are mere warnings or notices
  ! This routine converts these values to a severity level
  ! By default, it then prints only if severity is MLSMSG_Error or higher
  ! but by setting threshold this can be changed
  ! E.g., setting threshold to 0 will print every status, even success

  subroutine ReportTKStatus( status, ModuleNameIn, Message, Threshold )
    ! Dummy arguments
    integer, intent(in) :: status ! e.g. PGS_TD_NOLEAPSECFILE
    character (len=*), intent(in) :: ModuleNameIn ! Name of module (see below)
    character (len=*), intent(in) :: Message ! Line of text
    integer, intent(in), optional :: Threshold ! Min severity to log message

    ! Internal variables
    integer :: levelmask
    integer :: myThreshold
    integer :: severity

    ! Executable code
    myThreshold = MLSMSG_Error
    if ( present(threshold) ) myThreshold = threshold
    levelMask = pgs_smf_teststatuslevel(status)
    severity = level2severity(levelMask)
    if ( severity < myThreshold ) return
    call MLSMessage( severity, ModuleNameIn, Message )
  end subroutine ReportTKStatus

  ! Private procedures
  !-----------------------------------------  accessDFACCToStr  -----
  function accessDFACCToStr ( dfacc ) result(str)

    ! This routine converts an hdf access type
    ! like DFACC_RDONLY into a string like 'rdonly'
    ! If access type is unrecognized, returns 'unknown'
    ! Args
    integer, intent(in)           :: dfacc
    character(len=8)              :: str
    ! Executable
    select case (dfacc)
    case (DFACC_CREATE)
      str = 'create'
    case (DFACC_RDONLY)
      str = 'rdonly'
    case (DFACC_RDWR)
      str = 'rdwrite'
    case default
      str = 'unknown' ! Why not ' '? Or '?'
    end select
  end function accessDFACCToStr

  ! -------------------- LastGasp -------------------
  ! We're a slave and we're about to expire
  ! Before we do, however, try to tell the master why
  subroutine LastGasp ( ModuleNameIn, Message )
    character (len=*), intent(in) :: ModuleNameIn ! Name of module (see below)
    character (len=*), intent(in) :: Message ! Line of text
    ! Local variables
    integer :: BUFFERID                 ! ID for buffer to send
    integer :: INFO                     ! Flag from PVM
    ! Executable code
    call PVMFInitSend ( PvmDataDefault, bufferID  )
    call PVMF90Pack ( SIG_AboutToDie, info  )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing about-to-die signal'  )
    call PVMF90Pack ( ModuleNameIn // trim(message), info  )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing last gasp message'  )
    call PVMFSend ( MLSMessageConfig%masterTid, InfoTag, info  )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'sending last gasp'  )
  end subroutine LastGasp

  !-----------------------------------------  level2severity  -----
  function level2severity ( level ) result(severity)

    ! This routine converts a toolkit levelmask to an mls severity
    ! Args
    integer, intent(in)           :: level
    integer :: severity
    select case (level)
    case (PGS_SMF_MASK_LEV_S)
      severity = MLSMSG_Success
    case (PGS_SMF_MASK_LEV_N)
      severity = MLSMSG_Debug
    case (PGS_SMF_MASK_LEV_E)
      severity = MLSMSG_Error
    case (PGS_SMF_MASK_LEV_F)
      severity = MLSMSG_Error
    case (PGS_SMF_MASK_LEV_W)
      severity = MLSMSG_Warning
    case (PGS_SMF_MASK_LEV_M)
      severity = MLSMSG_Info
    case default
      severity = MLSMSG_Info
    end select
  end function level2severity

  ! --------------------------------------------  PRINTITOUT  -----
  subroutine PRINTITOUT ( INLINE, SEVERITY, LINE_LEN, NOPREFIX  )
    ! In any way we're asked
    ! Args
    character(len=*), intent(in) :: INLINE
    integer, intent(in) :: SEVERITY
    integer, optional, intent(in) :: LINE_LEN
    logical, optional, intent(in) :: NOPREFIX
    ! Local variables
    character(len=len(inline)) :: Line
    logical :: log_it
    integer :: loggedLength
    character(len=len(inline)+len(MLSMessageConfig%prefix)) :: loggedLine
    integer :: ioerror
    integer :: maxLineLength
    logical :: myNoPrefix
    ! Executable
    if ( MLSMessageConfig%AsciifyMessages ) then
      line = asciify(inLine)
    else
      line = inLine
    endif
    loggedLength = len_trim(line)
    if ( present(line_len) ) loggedLength = line_len
    myNoPrefix = .false.
    if ( present(noPrefix) ) myNoPrefix = noPrefix
    loggedLine = line
    if ( TRIM(MLSMessageConfig%prefix) /= ' ' .and. .not. myNoPrefix ) then
      loggedLength = loggedLength + len_trim(MLSMessageConfig%prefix)
      loggedLine = TRIM(MLSMessageConfig%prefix) // &
           & TRIM(line)
    endif
    maxLineLength = min( loggedLength, len(loggedLine) )
    log_it = &
    & (MLSMessageConfig%useToolkit .and. UseSDPToolkit) &
    & .or. &
    & severity >= MLSMSG_Severity_to_quit
    if( log_it .and. loggedLength > 0 .and. MLSMessageConfig%useToolkit) then
      ioerror = PGS_SMF_GenerateStatusReport ( loggedLine(1:maxLineLength) )
    end if

    ! Now, if we're also logging to a file then write to that too.
    select case ( MLSMessageConfig%logFileUnit  )
    case ( 0 :  )
      write ( UNIT=max(MLSMessageConfig%logFileUnit,1), FMT=* ) TRIM(line)
    case ( STDOUTLOGUNIT  )
      if ( USEDEFAULTFORMATSTDOUT ) then
        write ( UNIT=*, FMT=* ) TRIM(line)
      else
        write ( UNIT=*, FMT='(a)' ) TRIM(line)
      endif
    case default ! DEFAULTLOGUNIT
    end select

  end subroutine PRINTITOUT

!=======================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MLSMessageModule
!=======================================================================

!
! $Log$
! Revision 2.42  2012/08/16 17:38:07  pwagner
! Refers to module variable constant instead of '-1'
!
! Revision 2.41  2011/10/10 23:56:02  pwagner
! Prevent writing non-ascii chars to stdout
!
! Revision 2.40  2010/05/23 03:11:06  honghanh
! Fix a bug on line 200, to check if useToolkit is true before we log
!
! Revision 2.39  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.38  2009/06/16 17:10:12  pwagner
! Can Capitalize messages in warnings summaries
!
! Revision 2.37  2007/11/08 00:02:08  pwagner
! asciify not used any more; removed
!
! Revision 2.36  2007/08/29 19:51:23  pwagner
! Worked around Intel quirk that wraps stdout when 'FMT=*'
!
! Revision 2.35  2007/08/27 23:53:44  pwagner
! Fixed many small bugs; now used MLSMessage.f9h
!
! Revision 2.34  2007/08/23 22:13:55  pwagner
! Fixed small bugs; add MLSMSG_Severity_to_walkback
!
! Revision 2.33  2007/08/17 00:29:32  pwagner
! MLSMessageCalls commands include 'depth', 'length', 'remain'
!
! Revision 2.32  2007/08/13 17:10:45  pwagner
! Implement MLSCallStack for printing walkback
!
! Revision 2.31  2007/03/23 00:15:34  pwagner
! Last-ditch attempt to prevent error when outputting extra-long lines
!
! Revision 2.30  2007/01/23 17:12:30  pwagner
! Restore ability to suppress Debugs
!
! Revision 2.29  2007/01/13 01:47:07  pwagner
! Added MLSMessageInternalFile function to return what would be logged
!
! Revision 2.28  2005/12/12 19:54:16  pwagner
! Correctly returns error when given PGS_SMF_MASK_LEV_E
!
! Revision 2.27  2005/12/10 00:23:21  pwagner
! Added ReportTKStatus
!
! Revision 2.26  2005/07/18 17:44:19  pwagner
! Bug allowed substring index past actual length; fixed
!
! Revision 2.25  2005/07/15 20:37:40  pwagner
! Handles warning messages longer than 80 chars better
!
! Revision 2.24  2005/07/15 20:03:21  pwagner
! A work-around for Lahey memory leak doing any(strarray == str)
!
! Revision 2.23  2005/06/22 17:25:50  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.22  2005/06/14 20:32:51  pwagner
! Many changes to accommodate the new fields in MLSFile_T
!
! Revision 2.21  2005/05/31 17:48:26  pwagner
! Added MLSFile as optional arg to be dumped on error
!
! Revision 2.20  2005/05/03 15:56:37  pwagner
! NAG doesnt like space between / and ) in array constructor
!
! Revision 2.19  2005/05/02 22:56:32  vsnyder
! Make MLSMessage generic
!
! Revision 2.18  2005/03/15 23:47:52  pwagner
! Slaves given last chance to send error message to master
!
! Revision 2.17  2005/03/12 00:46:41  pwagner
! limits to warnings now work correctly
!
! Revision 2.16  2005/03/10 00:29:20  pwagner
! Limits Warnings to 1st 1000 each message
!
! Revision 2.15  2004/08/19 00:18:16  pwagner
! New way to respond to severe errors-kaBOOM
!
! Revision 2.14  2002/10/08 00:09:12  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.13  2001/08/08 23:50:35  pwagner
! Added farewell optional arg to MLSMessageExit
!
! Revision 2.12  2001/07/16 23:44:10  pwagner
! Added MLSMessageQuit
!
! Revision 2.11  2001/05/15 20:32:35  pwagner
! Fix wrong unit nums in write and close; reset default
!
! Revision 2.10  2001/05/14 23:43:59  pwagner
! Added severity for reason to write to log
!
! Revision 2.9  2001/05/11 23:42:27  pwagner
! Prevented unwanted double printing
!
! Revision 2.8  2001/05/09 23:30:13  pwagner
! Detachable from toolkit
!
! Revision 2.7  2001/05/04 23:26:01  vsnyder
! Call Exit_with_status with nonzero status to terminate
!
! Revision 2.6  2001/04/20 20:43:15  pwagner
! Check severity against MLSMSG_Severity_to_quit
!
! Revision 2.5  2001/03/16 19:44:18  vsnyder
! Don't stop until advance='yes' -- i.e. not before the message is complete
!
! Revision 2.4  2001/02/23 00:14:54  vsnyder
! Maybe the coordination of output_m and MLSMessageModule is OK now...
!
! Revision 2.3  2001/02/22 23:27:16  vsnyder
! Correct routing of output through MLSMessage
!
! Revision 2.2  2000/10/04 18:06:39  vsnyder
! Added an optional "advance" argument to MLSMessage
!
! Revision 2.1  2000/10/03 01:34:10  vsnyder
! Corrected a spelling error, simplified MLSMessageClose, standardized
! some spelling and spacing.
!
! Revision 2.0  2000/09/05 17:41:06  dcuddy
! Change revision to 2.0
!
! Revision 1.10  2000/06/23 01:08:48  vsnyder
! Delete unused variables (except ID) to keep NAG f95 happy
!
