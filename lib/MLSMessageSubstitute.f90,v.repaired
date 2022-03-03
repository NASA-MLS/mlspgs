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

  ! use HighOutput, only: Banner
  ! use Intrinsic, only: L_HDFeos, L_HDF, L_Swath, L_Zonalavg, Lit_Indices
  use Machine, only: Crash_Burn_Rude=>crash_Burn, Exit_With_Status, Nevercrash
  use MLSCommon, only: MLSFile_T, MLSDebug, MLSVerbose, &
    & MLSDebugsticky, MLSVerboseSticky, DontCrashHere
  use MLSStrings_0, only: Capitalize, Lowercase
  use Printit_M, only: AssembleFullLine, Get_Config, LogUnitName, PrefixLen, &
    & MLSMSG_Allocate, MLSMSG_Deallocate, &
    & MLSMSG_Crash, MLSMSG_Debug, MLSMSG_Error, MLSMSG_Info, MLSMSG_Success, &
    & MLSMSG_Testwarning, MLSMSG_Warning, MLSMessageconfig_T, &
    & Bothlogunit, Defaultlogunit, Invalidlogunit, Prefixlen, &
    & Printitout, Sniprcsfrom, &
    & Stdoutlogunit, MLSMessageconfig, &
    & MLSMSG_Severity_So_Far, MLSMSG_Severity_To_Quit, MLSMSG_Severity_To_Walkback
  use SDPToolkit, only: PGSd_Pc_File_Path_Max, UseSDPToolkit
  implicit none

  private
  
  public :: MLSMSG_Allocate, MLSMSG_Deallocate, &
   & MLSMSG_Crash, MLSMSG_Debug, MLSMSG_Error, MLSMSG_Info, MLSMSG_Success, &
   & MLSMSG_Testwarning, MLSMSG_Warning, MLSMessageconfig_T, &
   & Defaultlogunit, Invalidlogunit, Prefixlen, &
   & Stdoutlogunit, MLSMessageconfig, &
   & MLSMSG_Severity_So_Far, MLSMSG_Severity_To_Quit, MLSMSG_Severity_To_Walkback

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
  integer, parameter :: L_HDFeos      = 0
  integer, parameter :: L_HDF         = L_HDFeos + 1
  integer, parameter :: L_Swath       = L_HDF + 1
  integer, parameter :: L_Zonalavg    = L_Swath + 1
  integer, dimension(1) :: Lit_Indices = 1

  interface AddRow
    module procedure AddRow_Chars
    module procedure AddRow_int
    module procedure AddRow_intarray
    module procedure AddRow_log
  end interface

  interface Banner
    module procedure Banner_Chars
    module procedure Banner_Chararray
  end interface

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
  
  subroutine AddRow_chars ( chars, how )
    character(len=*), intent(in) :: chars
    character(len=*), intent(in) :: how
  end subroutine AddRow_chars

  subroutine AddRow_intarray ( chars, how )
    character(len=*), intent(in) :: chars
    integer, dimension(:), intent(in) :: how
  end subroutine AddRow_intarray

  subroutine AddRow_int ( chars, how )
    character(len=*), intent(in) :: chars
    integer, intent(in) :: how
  end subroutine AddRow_int

  subroutine AddRow_log ( chars, how )
    character(len=*), intent(in) :: chars
    logical, intent(in) :: how
  end subroutine AddRow_log

  subroutine AddRow_header ( sep, border )
    character(len=*), intent(in) :: sep
    character(len=*), intent(in) :: border
  end subroutine AddRow_header

  subroutine AddRow_divider ( chars )
    character(len=*), intent(in) :: chars
  end subroutine AddRow_divider

  subroutine startTable
  end subroutine startTable

  subroutine outputTable ( sep, border )
    character(len=*), intent(in) :: sep
    character(len=*), intent(in) :: border
  end subroutine outputTable

  ! -----------------------------------------------------  BANNER  -----
  ! We put a simpler version of thsi family of routines here to avoid
  ! dragging in highOutput and everything else that entails
  !
  ! Surround your message with stars and stripes; e.g.,
  ! *-----------------------------------------------*
  ! *            Your message here                  *
  ! *-----------------------------------------------*
  ! proclaiming its great importance to an uncaring world.
  ! For multiline messages, you may divide them into elements of
  ! a character array, or else a longer character scalar and
  ! supply LineLength asking the routine to wrap at word boundaries
  subroutine BANNER_chars ( chars, &
    & columnRange, alignment, skips, lineLength, mode, pattern )
    character(len=*), intent(in)                :: CHARS
    ! If columnRange(1) < 1, just use starting columns; otherwise move to
    integer, dimension(2), optional, intent(in) :: COLUMNRANGE
    character(len=1), intent(in), optional      :: ALIGNMENT ! L, R, C, or J
    integer, optional, intent(in)               :: SKIPS ! How many spaces between chars
    integer, optional, intent(in)               :: LINELENGTH
    character (len=*), optional, intent(in)     :: mode ! if not 'hard'
    character (len=1), optional, intent(in)     :: pattern ! if not stripes
    !
    call PrintItOut ( '*-----------------------------------------------*', MLSMSG_Info )
    call PrintItOut ( chars, MLSMSG_Info )
    call PrintItOut ( '*-----------------------------------------------*', MLSMSG_Info )
  end subroutine BANNER_chars

  subroutine Banner_chararray ( charArray, &
    & columnRange, alignment, skips, pattern )
    character(len=*), dimension(:), intent(in)  :: CHARARRAY
    ! If columnRange(1) < 1, just use starting columns; otherwise move to
    integer, dimension(2), optional, intent(in) :: COLUMNRANGE
    character(len=1), intent(in), optional      :: ALIGNMENT ! L, R, C, or J
    integer, optional, intent(in)               :: SKIPS ! How many spaces between chars
    character (len=1), optional, intent(in)     :: pattern ! if not stripes
    ! Internal variables
    integer :: i
    ! Executable
    call PrintItOut ( '*-----------------------------------------------*', MLSMSG_Info )
    do i = 1, size(chararray)
      call PrintItOut ( charArray(i), MLSMSG_Info )
    enddo
    call PrintItOut ( '*-----------------------------------------------*', MLSMSG_Info )
  end subroutine Banner_chararray

  ! -------------------- Dump_Stack -------------------
  ! In the full module this is called when:
  ! We're a slave and we're about to expire
  ! Before we do, however, try to tell the master why
  subroutine Dump_Stack ( Top, Before, Where, Size, SysSize, CPU, DoDepth, &
    & Rev, Index, String, StringIndex, ShowTime, Used, Advance, &
    & PrintMemoryReport, StackIsEmpty )
    logical, intent(in), optional :: Top   ! Dump only the top frame
    character(len=*), intent(in), optional :: Before ! first thing output
    logical, intent(in), optional :: Where ! Dump tree location
    logical, intent(in), optional :: Size  ! Dump memory size (default true)
    logical, intent(in), optional :: SysSize ! Dump memory size, as the system
                                           ! accounts for it, in kB (default
                                           ! Show_Sys_Memory)
    logical, intent(in), optional :: CPU   ! Print CPU (default false)
    logical, intent(in), optional :: DoDepth ! Print "depth" dots (default true)
    logical, intent(in), optional :: Rev   ! Print in reverse order (default false)
    integer, intent(in), optional :: Index ! Print this instead of from stack
    character(len=*), optional, intent(in) :: String
    integer, optional, intent(in) :: StringIndex
    logical, intent(in), optional :: ShowTime   ! Show time when we dumped
    character(len=*), optional :: Used
    character(len=*), intent(in), optional :: Advance ! Default 'yes'
    logical, intent(in), optional :: PrintMemoryReport! Report on mem usage so far
    logical, intent(out), optional :: StackIsEmpty ! Return true if stack is empty
    !
    call MLSMessageCalls ( 'dump' )
  end subroutine Dump_Stack

  ! -------------------- Get_String -------------------
  ! In the full module this is called when:
  ! We're a slave and we're about to expire
  ! Before we do, however, try to tell the master why
  subroutine Get_String ( strIndex, string, strip )
    integer, intent(in)              :: strIndex
    character(len=*), intent(out)    :: string
    logical, intent(in), optional    :: strip
    !
    string = ' '
  end subroutine Get_String

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
! Revision 2.24  2022/03/03 22:25:42  pwagner
! Added args mising from Dump_Stack
!
! Revision 2.23  2019/04/09 20:34:50  pwagner
! Moved some procedures from MLSStrings to new MLSStrings_0
!
! Revision 2.22  2019/03/18 22:04:40  pwagner
! Try harder to Print Errors to stdout
!
! Revision 2.21  2019/01/24 18:30:21  pwagner
! Reorganized modules that print to simplify toolkit-free builds
!
! Revision 2.20  2018/03/27 22:57:39  pwagner
! updated api for Dump_Stack; Freed from use-ing modules string_table and intrinsic
!
! Revision 2.19  2018/03/15 16:40:07  pwagner
! Moved 'Use' statement to .f90 where make can see it
!
! Revision 2.18  2016/10/10 22:42:06  pwagner
! Avoid the Call_Stack module if no toolkit
!
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
