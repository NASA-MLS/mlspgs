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

  use Call_stack_m, only: dump_stack
  use Highoutput, only: banner
  use Machine, only: crash_burn_rude=>crash_burn, exit_with_status, nevercrash
  use MLSCommon, only: MLSFile_t, MLSDebug, MLSVerbose, &
    & MLSDebugsticky, MLSVerboseSticky, dontCrashHere
  use MLSStrings, only: capitalize
  use Printit_m, only: assembleFullLine, get_config, logUnitName, prefixLen, &
    & MLSMSG_allocate, MLSMSG_deallocate, &
    & MLSMSG_crash, MLSMSG_debug, MLSMSG_error, MLSMSG_info, MLSMSG_success, &
    & MLSMSG_testwarning, MLSMSG_warning, MLSMessageconfig_t, &
    & Defaultlogunit, invalidlogunit, prefixlen, &
    & Printitout, sniprcsfrom, &
    & Stdoutlogunit, MLSMessageconfig, &
    & MLSMSG_severity_so_far, MLSMSG_severity_to_quit, MLSMSG_severity_to_walkback
  implicit none

  private
  
  public :: MLSMSG_ALLOCATE, MLSMSG_DEALLOCATE, &
    & MLSMSG_CRASH, MLSMSG_DEBUG, MLSMSG_ERROR, MLSMSG_INFO, MLSMSG_SUCCESS, &
    & MLSMSG_TESTWARNING, MLSMSG_WARNING, MLSMESSAGECONFIG_T, &
    & DEFAULTLOGUNIT, INVALIDLOGUNIT, PREFIXLEN, &
    & STDOUTLOGUNIT, MLSMESSAGECONFIG, &
    & MLSMSG_SEVERITY_SO_FAR, MLSMSG_SEVERITY_TO_QUIT, MLSMSG_SEVERITY_TO_WALKBACK

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! A low-weight substitute for the full module MLSMessageModule.f90.
  ! which provides low level messaging for the MLSPGS suite.  The main
  ! routine is MLSMessage, which generates log messages as directed by the
  ! user. In the high-fat module, which needs the toolkit,
  ! MLSMessage routine logs a message using the SDPToolkit routine
  ! PGS_SMF_GenerateStatusReport.  This writes a string to the `LogReport'
  ! file (PCF# 10101) in the toolkit.  
  
  ! Here, however, the Toolkit `substitute' just does a simple print.
  
  ! For alternate procedures see the toc in MLSMessage.f9h
  
  ! We dispense with most of the toolkit panoply, needing only the modules
  ! Machine
  ! MLSKinds
  ! MLSCommon
  ! CALL_STACK_M
  ! PRINTIT_M
  ! MLSSTRINGS

  integer, public, parameter ::  PGS_S_SUCCESS = 0

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

    ! Executable code
    call MLSMessage( MLSMSG_Info, ModuleNameIn, Message )
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
    case default
      str = 'unknown' ! Why not ' '? Or '?'
    end select
  end function accessDFACCToStr

  ! -------------------- LastGasp -------------------
  ! In the full module this is called when:
  ! We're a slave and we're about to expire
  ! Before we do, however, try to tell the master why
  subroutine LastGasp ( ModuleNameIn, Message )
    character (len=*), intent(in) :: ModuleNameIn ! Name of module (see below)
    character (len=*), intent(in) :: Message ! Line of text
    call MLSMessage( MLSMSG_Info, ModuleNameIn, 'No pvm to send ' // Message )
  end subroutine LastGasp

  !-----------------------------------------  level2severity  -----
  function level2severity ( level ) result(severity)

    ! In the full module:
    ! This routine converts a toolkit levelmask to an mls severity
    ! Args
    integer, intent(in)           :: level
    integer :: severity
    severity = MLSMSG_Info
  end function level2severity

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
! Revision 2.17  2015/06/30 18:41:29  pwagner
! Added a conditional crash_burn
!
! Revision 2.16  2014/04/22 16:29:33  pwagner
! Added bummer--writes error message as an eye-catching banner
!
! Revision 2.15  2013/11/18 22:23:34  pwagner
! Sticky versions of verbose, debug available
!
! Revision 2.14  2013/11/15 00:03:51  pwagner
! Comments confess need for 3 more modules; LastGasp gasps its inability
!
! Revision 2.13  2013/11/13 21:41:35  pwagner
! Compatible with PRINTIT_M, MLSMessage.f9h
!
! Revision 2.12  2013/08/28 00:35:39  pwagner
! Moved more stuff from MLSMessage down to PrintIt module
!
! Revision 2.11  2013/08/23 02:51:04  vsnyder
! Move PrintItOut to PrintIt_m
!
! Revision 2.10  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.9  2009/06/16 17:10:12  pwagner
! Can Capitalize messages in warnings summaries
!
! Revision 2.8  2007/11/08 00:02:08  pwagner
! asciify not used any more; removed
!
! Revision 2.7  2007/08/29 19:51:30  pwagner
! Worked around Intel quirk that wraps stdout when 'FMT=*'
!
! Revision 2.6  2007/08/27 23:53:37  pwagner
! Fixed many small bugs; now used MLSMessage.f9h
!
! Revision 2.5  2007/08/23 22:14:07  pwagner
! Fixed small bugs; add MLSMSG_Severity_to_walkback
!
! Revision 2.4  2007/08/17 00:29:32  pwagner
! MLSMessageCalls commands include 'depth', 'length', 'remain'
!
! Revision 2.3  2007/08/13 17:11:07  pwagner
! Implement MLSCallStack for printing walkback
!
! Revision 2.2  2007/01/23 17:18:04  pwagner
! Fixed an obvious bug; now compiles successfully
!
! Revision 2.1  2007/01/12 00:25:26  pwagner
! First commit
!
