! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!==============================================================================
module MLSMessageModule         ! Basic messaging for the MLSPGS suite
!==============================================================================

  use Machine, only: Exit_with_status
  use SDPToolkit, only: PGS_SMF_GenerateStatusReport, UseSDPToolkit

  implicit none

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  ! This module provides low level messaging for the MLSPGS suite.  The main
  ! routine is MLSMessage, which generates log messages as directed by the
  ! user. The MLSMessage routine logs a message using the SDPToolkit routine
  ! PGS_SMF_GenerateStatusReport.  This writes a string to the `LogReport'
  ! file (PCF# 10101) in the toolkit.  In the Toolkit `substitute' it just
  ! does a simple print.

  ! The user can also choose to log the messages to a seperate file when
  ! running under the toolkit.  This is setup by MLSMessageSetup and closed
  ! by MLSMessageClose.  The cataloging of such a file is left up to the
  ! calling code.

  ! ---------------------------------------------------------------------------

  ! Define some low level parameters.  These are used by the calling code to
  ! indicate the severity or otherwise of the messages.

  integer, parameter :: MLSMSG_Debug=1
  integer, parameter :: MLSMSG_Info=2
  integer, parameter :: MLSMSG_Warning=3
  integer, parameter :: MLSMSG_Error=4

  ! MLSMSG_Severity_to_quit can be reset in a main program to cause us
  ! to become more lenient (set it higher) or strict (set it lower)
  integer            :: MLSMSG_Severity_to_quit = MLSMSG_Error

  private :: SeverityNames
  character (len=*), dimension(4), parameter :: SeverityNames = &
     & (/"Debug  ","Info   ","Warning","Error  "/)

  ! This set of parameters are simple prefixes for common messages

  character (len=*), parameter :: MLSMSG_Allocate = &
     & "Allocation failed: "
  character (len=*), parameter :: MLSMSG_Fileopen = &
     & "Failed to open file: "
  character (len=*), parameter :: MLSMSG_Keyword = &
     & "Unrecognized configuration file keyword: "
  character (len=*), parameter :: MLSMSG_L1BRead = &
     & "Unable to read L1B data item: "
  character (len=*), parameter :: MLSMSG_Duplicate = &
     & "There is already an entry with the name "
  character (len=*), parameter :: MLSMSG_DeAllocate = &
     & "Deallocation failed: "
  ! This datatype describes the configuration of the messaging suite

  integer, private, parameter :: MLSMSG_PrefixLen = 32

   ! May get some of these from MLSLibOptions? 
  type MLSMessageConfig_T
    integer :: logFileUnit                     = -1
    character (len=MLSMSG_PrefixLen) :: prefix = ''
    logical :: suppressDebugs                  = .false.
    logical :: useToolkit                      = .true.
  end type MLSMessageConfig_T

  ! This variable describes the configuration

  type (MLSMessageConfig_T), public, save :: MLSMessageConfig

contains

  ! -------------------------------------------------  MLSMessage  -----

  ! This first routine is the main `messaging' code.

  subroutine MLSMessage ( Severity, ModuleNameIn, Message, Advance )

    ! Dummy arguments
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module (see below)
    character (len=*), intent(in) :: Message ! Line of text
    character (len=*), intent(in), optional :: Advance ! Do not advance
    !                                 if present and the first character is 'N'
    !                                 or 'n'

    ! Local variables
    integer :: Dummy
    character (len=512), save :: Line   ! Line to output, should be long enough
    integer, save :: Line_len=0         ! Number of saved characters in line.
    !                                     If nonzero, do not insert prefix.
    logical :: log_it
    logical :: My_adv

    ! Executable code

    my_adv = .true.
    if ( present(advance) ) &
      & my_adv = advance(1:1) /= 'n' .and. advance(1:1) /= 'N'

    if ( (.not. MLSMessageConfig%suppressDebugs).OR. &
         & (severity /= MLSMSG_Debug) ) then
       
       ! Assemble a full message line

       if ( line_len == 0 ) then
         line_len = len_trim(SeverityNames(severity))
         line = SeverityNames(severity)
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
       line(line_len+1:) = message
       line_len = line_len + len(message) ! Not len-trim, so we can get
       ! trailing blanks into a part of a message.  If there are trailing
       ! blanks remaining when my_adv is true, they'll be trimmed off.

       ! Log the message using the toolkit routine
       ! (or its substitute)
       ! if either using toolkit or severity is sufficient to
       ! quit (which means we might have been called directly
       ! rather than from output module)

       if ( my_adv ) then
         log_it = &
         & (MLSMessageConfig%useToolkit .and. UseSDPToolkit) &
         & .or. &
         & severity >= MLSMSG_Severity_to_quit
         if( log_it ) &
         & dummy = PGS_SMF_GenerateStatusReport&
         & (TRIM(MLSMessageConfig%prefix)// &
              & TRIM(line))

         ! Now, if we're also logging to a file then write to that too.

         select case ( MLSMessageConfig%logFileUnit )
         case ( 0 : )
           write ( UNIT=MLSMessageConfig%logFileUnit, FMT=* ) TRIM(line)
         case ( -1 )
           write ( UNIT=*, FMT=* ) TRIM(line)
         case default
         end select

         line_len = 0
         line = ' '
       end if

    end if

    ! Now if it's an error, and the message is complete, then try to close
    ! log file if any and quit

    if ( my_adv .and. severity >= MLSMSG_Severity_to_quit ) then
      if ( MLSMessageConfig%logFileUnit /= -1 ) &
        & close ( MLSMessageConfig%logFileUnit )
      call exit_with_status ( 1 )
    end if
  end subroutine MLSMessage

  ! --------------------------------------------  MLSMessageSetup  -----

  ! This routine sets up the MLSMessage suite.  The defaults are of course
  ! sensible, but the user may wish to change things.

  subroutine MLSMessageSetup ( SuppressDebugs, LogFileUnit, Prefix, useToolkit )

    ! Dummy arguments
    logical, optional, intent(in) :: SuppressDebugs
    integer, optional, intent(in) :: LogFileUnit
    character (len=*), optional, intent(in) :: Prefix
    logical, optional, intent(in) :: useToolkit

    ! Local variables

    ! Executable code

    if ( present(suppressDebugs) ) &
      & MLSMessageConfig%suppressDebugs=suppressDebugs

    if ( present(prefix) ) &
      & MLSMessageConfig%prefix=prefix

    if ( present(logFileUnit) ) then
      if ( MLSMessageConfig%logFileUnit /= -1 ) call MLSMessage ( &
        & MLSMSG_Error, ModuleName,"Already writing to a log file")
      MLSMessageConfig%logFileUnit = logFileUnit
    end if

    if ( present(useToolkit) ) &
      & MLSMessageConfig%useToolkit=useToolkit

  end subroutine MLSMessageSetup

  ! --------------------------------------------  MLSMessageClose  -----

  ! This routine simply closes the MLSMessage log file if there is one.

  subroutine MLSMessageClose
    ! Executable code
    MLSMessageConfig%logFileUnit=-1
  end subroutine MLSMessageClose

!=======================================================================
end module MLSMessageModule
!=======================================================================

!
! $Log$
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
