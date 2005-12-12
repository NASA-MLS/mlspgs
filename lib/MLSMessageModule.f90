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
  use Machine, only: CRASH_BURN, Exit_with_status
  use MLSCommon, only: MLSFile_T
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
  ! (to report on the severity of a PGS Toolkit return status)

  ! The user can also choose to log the messages to a seperate file when
  ! running under the toolkit.  This is setup by MLSMessageSetup and closed
  ! by MLSMessageClose.  The cataloging of such a file is left up to the
  ! calling code.

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (parameters)
! MLSMSG_SUCCESS           status returned when all went well
! MLSMSG_DEBUG             should print only if debugging turned on
! MLSMSG_INFO              fyi only
! MLSMSG_WARNING           not fatal, but deserving of attention
! MLSMSG_ERROR             quits after printing
! MLSMSG_CRASH             should give traceback before quitting
! MLSMSG_Severity_to_quit  severity level needed to quit
! MLSMSG_Allocate          mesg prefix for this type of error
! MLSMSG_Fileopen          mesg prefix for this type of error
! MLSMSG_Keyword           mesg prefix for this type of error
! MLSMSG_L1BRead           mesg prefix for this type of error
! MLSMSG_Duplicate         mesg prefix for this type of error
! MLSMSG_DeAllocate        mesg prefix for this type of error
! MLSMSG_PVM               mesg prefix for this type of error
! MLSMessageConfig         configuration controlling where to print, etc.

!     (subroutines and functions)
! MLSMessage               main messaging routine
! MLSMessageSetup          routine interface to change some parts of MLSMessageConfig
! MLSMessageClose          close MLSMessage log file; but see MLSMessageExit
! MLSMessageExit           recommended way to finish main program
! MLSMessageReset          reset flags, counters, etc. during runtime
! PVMErrorMessage          log a PVM error
! ReportTKStatus           converts SDP status to severity, prints if needed

! === (end of toc) ===

! === (start of api) ===
! MLSMessage ( int Severity, char* ModuleNameIn, char* Message, 
!      [char* Advance], [MLSFile_T MLSFile] ) 
! MLSMessageSetup ( [log SuppressDebugs], [int LogFileUnit], [char* Prefix],
!      [log useToolkit], [log CrashOnAnyError] )
! MLSMessageExit ( [int status], [char* farewell] )
! MLSMessageReset ( [int logFileUnit], [log CrashOnAnyError], [log Warnings] )
! PVMErrorMessage ( int INFO, char* PLACE  )
! ReportTKStatus( int status, char* ModuleNameIn, char* Message, 
!      [int Threshold] )
! === (end of api) ===
  ! ---------------------------------------------------------------------------

  ! Define some low level parameters.  These are used by the calling code to
  ! indicate the severity or otherwise of the messages.

  integer, public, parameter :: MLSMSG_Success = PGS_S_SUCCESS ! == 0
  integer, public, parameter :: MLSMSG_Debug   = MLSMSG_Success + 1
  integer, public, parameter :: MLSMSG_Info    = MLSMSG_Debug + 1
  integer, public, parameter :: MLSMSG_Warning = MLSMSG_Info + 1
  integer, public, parameter :: MLSMSG_Error   = MLSMSG_Warning + 1
  ! Warning--a Crash may not properly close files opened by your run
  ! Use it only for specific debugging where you need a walkback
  ! See also MLSMessageConfig%crashOnAnyError
  integer, public, parameter :: MLSMSG_Crash   = MLSMSG_Error + 1

  ! MLSMSG_Severity_to_quit can be reset in a main program to cause us
  ! to become more lenient (set it higher) or strict (set it lower )
  integer, public            :: MLSMSG_Severity_to_quit = MLSMSG_Error

  private :: SeverityNames
  character (len=*), dimension(MLSMSG_Success:MLSMSG_Crash), parameter :: &
     & SeverityNames = (/&
     & "Success", &
     & "Debug  ", &
     & "Info   ", &
     & "Warning", &
     & "Error  ", &
     & "Crash  " &
     /)

  ! So that we may limit the number of times warnings printed, messagewise
  character(len=*), parameter :: WARNINGSSUPPRESSED = '(No more warnings of this)'
  integer, parameter :: MAXNUMWARNINGS = 40 ! was 80, but most tests < 10
  integer, parameter :: WARNINGMESSLENGTH = 80
  character(len=WARNINGMESSLENGTH), dimension(MAXNUMWARNINGS), save :: &
    &                   warningmessages = ' '
  integer, dimension(MAXNUMWARNINGS), save :: timeswarned = 0
  integer, save :: numwarnings = 0
  ! This set of parameters are simple prefixes for common messages

  character (len=*), public, parameter :: MLSMSG_Allocate = &
     & "Allocation failed: "
  character (len=*), public, parameter :: MLSMSG_Fileopen = &
     & "Failed to open file: "
  character (len=*), public, parameter :: MLSMSG_Keyword = &
     & "Unrecognized configuration file keyword: "
  character (len=*), public, parameter :: MLSMSG_L1BRead = &
     & "Unable to read L1B data item: "
  character (len=*), public, parameter :: MLSMSG_Duplicate = &
     & "There is already an entry with the name "
  character (len=*), public, parameter :: MLSMSG_DeAllocate = &
     & "Deallocation failed: "
  character (len=*), public, parameter :: MLSMSG_PVM = &
     & "PVM Error: "
  ! This datatype describes the configuration of the messaging suite

  integer, private, parameter :: MLSMSG_PrefixLen = 32

   ! May get some of these from MLSLibOptions? 
  type, public :: MLSMessageConfig_T
    ! We log messages by toolkit (if useToolkit and UseSDPToolkit are TRUE )
    ! In the following, values would have the effect of adding logged messages:
    ! -2: none added
    ! -1: to stdout
    !  n: to ftn unit n
    integer :: logFileUnit                     = -2
    ! In the following, values would have the effect on identical warnings of:
    ! -1: Print every one without suppression
    !  0: Suppress every one
    !  1: Print every one only once
    integer :: limitWarnings                   = 1000 ! Max number each warning
    integer :: masterTID                       = -1 ! Where to send error msg
    character (len=MLSMSG_PrefixLen) :: prefix = ''   ! Prefix to every msg
    logical :: suppressDebugs                  = .false.
    logical :: useToolkit                      = .true.
    logical :: CrashOnAnyError                 = .false. ! See crash warning
    logical :: SendErrMsgToMaster              = .false. ! Whether to send last
  end type MLSMessageConfig_T

  ! This variable describes the configuration

  type (MLSMessageConfig_T), public, save :: MLSMessageConfig
  
  ! Public procedures
  public :: MLSMessage, MLSMessage_, MLSMessageSetup, MLSMessageClose
  public :: MLSMessageExit, MLSMessageReset, PVMErrorMessage
  public :: ReportTKStatus

  interface MLSMessage
    module procedure MLSMessage_
  end interface

contains

  ! ------------------------------------------------  MLSMessage_  -----

  ! This first routine is the main `messaging' code.

  subroutine MLSMessage_ ( Severity, ModuleNameIn, Message, Advance, MLSFile )

    ! Dummy arguments
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module (see below)
    character (len=*), intent(in) :: Message ! Line of text
    character (len=*), intent(in), optional :: Advance ! Do not advance
    !                                 if present and the first character is 'N'
    !                                 or 'n'
    type(MLSFile_T), intent(in), optional :: MLSFile

    ! Local variables
    character (len=512), save :: Line   ! Line to output, should be long enough
    integer, save :: Line_len=0         ! Number of saved characters in line.
    !                                     If nonzero, do not insert prefix.
    integer :: msgLength                  
    logical :: My_adv
    logical :: nosubsequentwarnings
    logical :: newwarning
    integer :: warning_index

    ! Executable code

    my_adv = .true.
    if ( present(advance) ) &
      & my_adv = advance(1:1) /= 'n' .and. advance(1:1) /= 'N'
    ! This is the smaller of the actual length and what we can check for repeats
    msgLength = min(len(message), WARNINGMESSLENGTH)

    ! Here's where we suppress warning messages beyond a limit for each
    nosubsequentwarnings = .false.
    if ( severity == MLSMSG_Warning .and. MLSMessageConfig%limitWarnings > -1 &
      & .and. numwarnings <= MAXNUMWARNINGS .and. message /= ' ' ) then
      ! See if we have seen this message before
      ! newwarning = .not. any(warningmessages == trim(message))
      newwarning = .true.
      do warning_index = 1, numwarnings
        newwarning = newwarning .and. &
          & ( warningmessages(warning_index) /= trim(message(1:msgLength)) )
      enddo
      if ( newwarning .and. numwarnings >= MAXNUMWARNINGS ) then
      else if ( newwarning .or. &
        & numwarnings < 1 ) then
        numwarnings = numwarnings + 1
        warningmessages(numwarnings) = message
        timeswarned(numwarnings) = timeswarned(numwarnings) + 1
        if ( timeswarned(numwarnings) > MLSMessageConfig%limitWarnings ) return
        timeswarned(numwarnings) = min(timeswarned(numwarnings) + 1, &
          & MLSMessageConfig%limitWarnings + 1 )
        nosubsequentwarnings = &
          & (timeswarned(numwarnings) >= MLSMessageConfig%limitWarnings)
      else
        do warning_index = 1, numwarnings
          if ( warningmessages(warning_index) == message(1:msgLength) ) exit
        end do
        if ( warning_index > numwarnings ) return
        if ( timeswarned(warning_index) > MLSMessageConfig%limitWarnings ) return
        timeswarned(warning_index) = min(timeswarned(warning_index) + 1, &
          & MLSMessageConfig%limitWarnings + 1 )
        nosubsequentwarnings = &
          & (timeswarned(warning_index) >= MLSMessageConfig%limitWarnings)
      end if
    end if
      
    if ( (.not. MLSMessageConfig%suppressDebugs).OR. &
         & (severity /= MLSMSG_Debug) ) then
       
       ! Assemble a full message line

       if ( line_len == 0 ) then
         if ( severity > MLSMSG_Success-1 .and. severity < MLSMSG_Crash+1 ) then
           line = SeverityNames(severity)
         else
           line = 'Unknown'
         end if
         line_len = len_trim(line)
         line(line_len+1:line_len+2) = ' ('
         if ( moduleNameIn(1:1) == '$' ) then
         ! The moduleNameIn is <dollar>RCSFile: <filename>,v <dollar>
           line(line_len+3:) = moduleNameIn(11:(LEN_TRIM(moduleNameIn)-8))
         else
           line(line_len+3:) = moduleNameIn
         end if
         line_len = len_trim(line) + 3
         line(line_len-2:line_len-1) = '):'
       end if
       if ( nosubsequentwarnings ) then
         line(line_len+1:) = WARNINGSSUPPRESSED // message
         line_len = line_len + len(WARNINGSSUPPRESSED) + len_trim(message)
       else
         line(line_len+1:) = message
         line_len = line_len + len(message) ! Not len-trim, so we can get
         ! trailing blanks into a part of a message.  If there are trailing
         ! blanks remaining when my_adv is true, they'll be trimmed off.
       end if

       ! Log the message using the toolkit routine
       ! (or its substitute )
       ! if either using toolkit or severity is sufficient to
       ! quit (which means we might have been called directly
       ! rather than from output module )

       if ( my_adv ) then
         call printitout(line, severity)
         line_len = 0
         line = ' '
       end if

    end if

    ! Now if it's an error, and the message is complete, then try to close
    ! log file if any and quit (or crash)

    if ( my_adv .and. severity >= MLSMSG_Severity_to_quit ) then
      if ( present(MLSFile) ) call dumpFile(MLSFile)
      if ( MLSMessageConfig%SendErrMsgToMaster .and. &
        & MLSMessageConfig%masterTID > 0 ) call LastGasp(ModulenameIn, Message )
      if ( MLSMessageConfig%logFileUnit > 0 ) &
        & close ( MLSMessageConfig%logFileUnit  )
      if ( severity >= MLSMSG_Crash .or. MLSMessageConfig%CrashOnAnyError ) &
        & call crash_burn
      call exit_with_status ( 1  )
    end if
  end subroutine MLSMessage_

  ! --------------------------------------------  MLSMessageSetup  -----

  ! This routine sets up the MLSMessage suite.  The defaults are of course
  ! sensible, but the user may wish to change things.

  subroutine MLSMessageSetup ( SuppressDebugs, LogFileUnit, Prefix, useToolkit, &
    & CrashOnAnyError  )

    ! Dummy arguments
    logical, optional, intent(in) :: SuppressDebugs
    integer, optional, intent(in) :: LogFileUnit
    character (len=*), optional, intent(in) :: Prefix
    logical, optional, intent(in) :: useToolkit
    logical, optional, intent(in) :: CrashOnAnyError

    ! Local variables

    ! Executable code

    if ( present(suppressDebugs) ) &
      & MLSMessageConfig%suppressDebugs=suppressDebugs

    if ( present(prefix) ) &
      & MLSMessageConfig%prefix=prefix

    if ( present(logFileUnit) ) then
      if ( MLSMessageConfig%logFileUnit /= -1 ) call MLSMessage ( &
        & MLSMSG_Error, ModuleName,"Already writing to a log file" )
      MLSMessageConfig%logFileUnit = logFileUnit
    end if

    if ( present(useToolkit) ) &
      & MLSMessageConfig%useToolkit=useToolkit
    if ( present(CrashOnAnyError) ) &
      & MLSMessageConfig%CrashOnAnyError=CrashOnAnyError

  end subroutine MLSMessageSetup

  ! --------------------------------------------  MLSMessageClose  -----

  ! This routine simply closes the MLSMessage log file if there is one.

  subroutine MLSMessageClose
    ! Executable code
    MLSMessageConfig%logFileUnit=-1
  end subroutine MLSMessageClose

  ! --------------------------------------------  MLSMessageExit  -----

  ! This routine (optionally) logs farewell, advances
  ! (hopefully) gracefully ends logging, and exits 
  ! (optionally with status )
  ! if farewell present, and non-blank, logs it
  ! if farewell present,  but blank, logs default message
  ! if farewell absent, does not log
  subroutine MLSMessageExit ( status, farewell )
  integer, optional, intent(in) :: STATUS
  character(LEN=*), optional, intent(in) :: FAREWELL
  CHARACTER(LEN=36) :: mesg

    ! Executable code
    if(present(status)) then
      if(present(farewell)) then
        if(farewell == ' ') then
          write(mesg, '(A29, I2, A1)') 'Exiting with status (', &
          status, ')'
          call MLSMessage ( MLSMSG_Info, ModuleName, mesg, advance='y' )
        else
          call MLSMessage ( MLSMSG_Info, ModuleName, farewell, advance='y' )
        end if
      end if
      call MLSMessageClose
      call exit_with_status ( status  )
    else
      if(present(farewell)) then
        if(farewell == ' ') then
          mesg='Exiting normally with "stop"'
          call MLSMessage ( MLSMSG_Info, ModuleName, mesg, advance='y' )
        else
          call MLSMessage ( MLSMSG_Info, ModuleName, farewell, advance='y' )
        end if
      end if
      call MLSMessageClose
      stop
    end if
  end subroutine MLSMessageExit

  ! --------------------------------------------  MLSMessageReset  -----

  ! This routine allows you to reset flags, counters, etc. during runtime

  subroutine MLSMessageReset ( logFileUnit, CrashOnAnyError, Warnings )
    ! Args
    integer, intent(in), optional :: logFileUnit
    logical, intent(in), optional :: CrashOnAnyError
    logical, intent(in), optional :: Warnings
    character(len=6) :: logname
    ! Executable code
    if ( present(logFileUnit) ) then
      if ( logFileUnit /= MLSMessageConfig%logFileUnit ) then
        write(logname, '(i6)') MLSMessageConfig%logFileUnit
        call MLSMessage ( MLSMSG_Info, ModuleName, &
          & 'Closing output on' // logname )
        call MLSMessageClose
        MLSMessageConfig%logFileUnit = logFileUnit
        write(logname, '(i6)') MLSMessageConfig%logFileUnit
        call MLSMessage ( MLSMSG_Info, ModuleName, &
          & 'Opening output on' // logname )
      end if
    end if
    if ( present(CrashOnAnyError) ) MLSMessageConfig%CrashOnAnyError = CrashOnAnyError
    if ( present(Warnings) ) then
      numwarnings = 0
      timeswarned = 0
      warningmessages = ' '
    end if
  end subroutine MLSMessageReset

  ! --------------------------------------------  PVMERRORMESSAGE  -----
  subroutine PVMErrorMessage ( INFO, PLACE  )
    ! This routine is called to log a PVM error
    integer, intent(in) :: INFO
    character (LEN=*) :: PLACE

    character (LEN=132) :: LINE

    write (line, * ) info
    call MLSMessage ( MLSMSG_Error, Place, MLSMSG_PVM // &
      & ' Info='//trim(adjustl(line)))
  end subroutine PVMErrorMessage

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

  ! --------------------------------------------  dumpFile  -----
  subroutine dumpFile ( MLSFile  )
    ! Show everything about it
    type(MLSFile_T) :: MLSFile
    integer, parameter :: SEVERE = MLSMSG_Error
    ! Executable code
    call printitout ( 'MLS File Info: ', MLSMSG_Error )                                  
    call dump ( '(name) ', charValue=trim(MLSFile%Name))                                  
    call dump ( 'short name ', charValue=trim(MLSFile%shortName))                                  
    call dump ( '    Type (int)   : ', MLSFile%Type)
    call dump ( '    Type         : ', charValue=trim(MLSFile%TypeStr))
    call dump ( '    Access       : ', charValue=trim(accessDFACCToStr(MLSFile%access)))
    call dump ( '    content      : ', charValue=trim(MLSFile%content))
    call dump ( '    last Operatn : ', charValue=trim(MLSFile%lastOperation))
    call dump ( '    File ID      : ', MLSFile%FileId%f_id)
    call dump ( '    Group ID     : ', MLSFile%FileId%grp_id)
    call dump ( '    DataSet ID   : ', MLSFile%FileId%sd_id)
    call dump ( '    PCF ID       : ', MLSFile%PCFId)
    call dump ( '    PCF Range    : ', MLSFile%PCFidRange%Bottom)
    call dump ( '                 : ', MLSFile%PCFidRange%Top)
    call dump ( '    hdf version  : ', MLSFile%HDFVersion)
    call dump ( '    record length: ', MLSFile%recordLength)
    call dump ( '    Open?        : ', logValue= MLSFile%StillOpen)
    call dump ( '    error code   : ', MLSFile%errorCode)
  end subroutine dumpFile

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

  ! --------------------------------------------  dump  -----
  subroutine dump ( name, intValue, charValue, logValue  )
    ! In any way we're asked
    character(len=*), intent(in) :: name
    integer, intent(in), optional :: intValue
    character(len=*), intent(in), optional :: charValue
    logical, intent(in), optional :: logValue
    !
    character(len=132) :: line
    character(len=32) :: value
    value = ''
    if ( present(intValue) ) then
      write(value, '(i10)') intValue
    elseif ( present(logValue) ) then
      write(value, '(l10)') logValue
    endif
    if ( present(charValue) ) then
      line = trim(name) // ' : ' // trim(charvalue)
    else
      line = trim(name) // ' : ' // trim(value)
    endif
    call printitout(trim(line), MLSMSG_Error)
  end subroutine dump

  ! --------------------------------------------  PRINTITOUT  -----
  subroutine PRINTITOUT ( LINE, SEVERITY  )
    ! In any way we're asked
    character(len=*), intent(in) :: LINE
    integer, intent(in) :: SEVERITY
    logical :: log_it
    integer :: ioerror
      log_it = &
      & (MLSMessageConfig%useToolkit .and. UseSDPToolkit) &
      & .or. &
      & severity >= MLSMSG_Severity_to_quit
      if( log_it ) &
      & ioerror = PGS_SMF_GenerateStatusReport&
      & (TRIM(MLSMessageConfig%prefix)// &
           & TRIM(line) )

      ! Now, if we're also logging to a file then write to that too.

      select case ( MLSMessageConfig%logFileUnit  )
      case ( 0 :  )
        write ( UNIT=max(MLSMessageConfig%logFileUnit,1), FMT=* ) TRIM(line)
      case ( -1  )
        write ( UNIT=*, FMT=* ) TRIM(line)
      case default
      end select

  end subroutine PRINTITOUT

!=======================================================================
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MLSMessageModule
!=======================================================================

!
! $Log$
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
