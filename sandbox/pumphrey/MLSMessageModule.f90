! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
! Hugh's private version
!==============================================================================
module MLSMessageModule         ! Basic messaging for the MLSPGS suite
!==============================================================================

  !use SDPToolkit
  !private
  !implicit none
  public::MLSMessage,MLSMessageSetup,MLSMessageClose

! ------------------------------- RCS Ident Info ------------------------------
character (len=130),private,parameter :: Id = &
      "$Id$"
character (len=*), parameter,private::&
     ModuleName= "$RCSfile$"
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

integer, parameter,public :: MLSMSG_Debug=1
integer, parameter,public  :: MLSMSG_Info=2
integer, parameter ,public :: MLSMSG_Warning=3
integer, parameter,public  :: MLSMSG_Error=4

character (len=*), dimension(4), parameter,private :: SeverityNames = &
     (/"Debug  ","Info   ","Warning","Error  "/)

! This set of parameters are simple prefixes for common messages

character (len=*), parameter ,public :: MLSMSG_Allocate = &
      "Allocation failed: "
character (len=*), parameter,public  :: MLSMSG_Fileopen = &
      "Failed to open file: "
character (len=*), parameter,public  :: MLSMSG_Keyword = &
      "Unrecognized configuration file keyword: "
character (len=*), parameter,public  :: MLSMSG_L1BRead = &
      "Unable to read L1B data item: "
character (len=*), parameter,public  :: MLSMSG_Duplicate = &
      "There is already an entry with the name "

! This datatype describes the configuration of the messaging suite

integer, parameter,private :: MLSMSG_PrefixLen=32

type,public:: MLSMessageConfig_T
   logical :: suppressDebugs
   integer :: logFileUnit
   character (len=MLSMSG_PrefixLen) :: prefix
end type MLSMessageConfig_T

! This private variable describes this configuration
type(MLSMessageConfig_T),public :: &
     MLSMessageConfig! =MLSMessageConfig_T(.true.,-1,"")

contains

  ! -------------------------------------------------------------------------

  ! This first routine is the main `messaging' code.

  function MLSMessage(severity,moduleNameIn,message) result(dummy)

    ! Dummy arguments
    integer, intent(in) :: severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: moduleNameIn ! Name of module (see below)
    character (len=*), intent(in) :: message ! Line of text

    ! Local variables
    integer :: dummy
    character (len=512) :: line ! Line to output, should be long enough

    ! Executable code

    ! The moduleNameIn is <dollar>RCSFile: <filename>,v <dollar>

    if ((.not. MLSMessageConfig%suppressDebugs).or. &
         (severity/=MLSMSG_Debug)) then

       ! Assemble a full message line

       line=trim(SeverityNames(severity))// &
            " ("//moduleNameIn(11:(len_trim(moduleNameIn)-8)) &
            //"): "//message

       ! Log the message using the toolkit routine

       print * ,trim(MLSMessageConfig%prefix)// &
            trim(line)

    endif

    ! Now if it's an error, then try to close log file if any and quit

    if (severity==MLSMSG_Error) then
       stop
    endif
    dummy=0
  end function MLSMessage

  ! ----------------------------------------------------------------------

  ! This routine sets up the MLSMessage suite.  The defaults are of course
  ! sensible, but the user may wish to chage things.

  subroutine MLSMessageSetup(suppressDebugs,logFileUnit,prefix)

    ! Dummy arguments
    logical, optional, intent(in) :: suppressDebugs
    integer, optional, intent(in) :: logFileUnit
    character (len=*), optional, intent(in) :: prefix

    ! Local variables
    !integer :: fileUnit, status
    integer::dummy
    ! Executable code

    if (present(suppressDebugs)) then
       MLSMessageConfig%suppressDebugs=suppressDebugs
    endif
    if (present(prefix)) then
       MLSMessageConfig%prefix=prefix
    endif
    if (present(logFileUnit)) then
       if (MLSMessageConfig%logFileUnit/=-1) then
          dummy=MLSMessage(MLSMSG_Error, &
               ModuleName,"Already writing to a log file")
       endif
       MLSMessageConfig%logFileUnit=logFileUnit
    end if
  end subroutine MLSMessageSetup

  ! ----------------------------------------------------------------------

  ! This routine simply closes the MLSMessage log file if there is one.

  subroutine MLSMessageClose()
    ! Executable code
    if (MLSMessageConfig%logFileUnit/=-1) then
       MLSMessageConfig%logFileUnit=-1
    endif
  end subroutine MLSMessageClose

!===========================================================================
end module MLSMessageModule
!===========================================================================

!
! $Log$
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
