! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

program tellMasterToQuit
  use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
  use dates_module, only: DATEFORM, REFORMATDATE
  use L2PARINFO, only: PARALLEL, INITPARALLEL
  use L2ParInfo, only: MACHINE_T, PARALLEL, &
    & PETITIONTAG, GIVEUPTAG, GRANTEDTAG, NOTIFYTAG, &
    & SIG_FINISHED, SIG_REGISTER, SIG_SWEARALLEGIANCE, SIG_SWITCHALLEGIANCE, &
    & SIG_HOSTDIED, SIG_RELEASEHOST, SIG_REQUESTHOST, SIG_THANKSHOST, &
    & MACHINENAMELEN, GETMACHINENAMES, &
    & DUMP, ADDMACHINETODATABASE
  use MACHINE ! At least HP for command lines, and maybe GETARG, too
  use MLSCOMMON, only: FILENAMELEN
  use MLSL2Options, only: CURRENT_VERSION_ID
  use MLSMessageModule, only: MLSMessage, MLSMessageConfig, MLSMessageExit, &
    & MLSMSG_Allocate, MLSMSG_DeAllocate, MLSMSG_Debug, MLSMSG_Error, &
    & MLSMSG_Info, MLSMSG_Success, MLSMSG_Warning, PVMERRORMESSAGE
  use MLSFINDS, only: FINDFIRST
  use MLSSTRINGLISTS, only: CATLISTS, GETSTRINGELEMENT, NUMSTRINGELEMENTS, &
    & STRINGELEMENTNUM
  use MLSSTRINGS, only: LOWERCASE, READINTSFROMCHARS, STREQ
  use OUTPUT_M, only: BLANKS, NEWLINE, &
    & OUTPUT, OUTPUT_DATE_AND_TIME, outputNamedValue, OutputOptions, &
    & TIMESTAMP
  use PVM, only: PVMOK, &
    & ClearPVMArgs, GETMACHINENAMEFROMTID, &
    & PVMDATADEFAULT, PVMFINITSEND, PVMF90PACK, PVMFKILL, PVMFMYTID, &
    & PVMF90UNPACK, PVMFPSTAT, &
    & PVMFSEND, PVMFNOTIFY, PVMTASKEXIT, &
    & PVMFFREEBUF
  use Sort_M, only: SORT
  use Time_M, only: Time_Now, time_config
  use TOGGLES, only: GEN, LEVELS, &
    & TOGGLE

  ! === (start of toc) ===
  !     c o n t e n t s
  !     - - - - - - - -

  ! Main program to tell master to quit
  ! === (end of toc) ===
  ! (It is assumed that pvm is already up and running)
  ! Usage:
  ! tellMasterToQuit [options] [--] [<] [list]
  ! where list is a list of master tids
  
  ! For a list of options 
  ! tellMasterToQuit --help
  

  implicit none

  integer :: BUFFERID
  character(len=2048) :: command_line ! All the opts
  logical, parameter :: COUNTEMPTY = .true.
  logical, parameter :: DEEBUG = .false.
  logical, parameter :: DUMPDBSONDEBUG = .false.
  integer :: ERROR
  integer :: I
  integer :: INFO
  integer :: INUNIT = -1       ! Input unit, * if < 0
  character(len=2048) :: LINE      ! Into which is read the command args
  integer, parameter :: LIST_UNIT = 20  ! Unit # for hosts file if not stdin
  character(len=MachineNameLen), dimension(:), pointer :: MACHINENAMES => null()
  integer, parameter :: MAXNUMMASTERS = 100 ! Mas num running simultaneously
  integer, parameter :: MAXNUMMULTIPROCS = 8 ! For some architectures > 1000 
  integer :: RECL = 10000          ! Record length for list
  integer :: STATUS                ! From OPEN
  logical :: SWITCH                ! "First letter after -- was not n"
  real :: T0, T1, T2, T_CONVERSION ! For timing
  integer :: TAG
  integer :: TID                   ! TID of master
  character(len=32) :: TIDSTR

  character(len=*), parameter :: GROUPNAME = "mlsl2"
  character(len=*), parameter :: LISTNAMEEXTENSION = ".txt"
  integer, parameter          :: AVOIDSELECTEDHOSTSTAG = GIVEUPTAG - 1
  integer, parameter          :: CHECKREVIVEDHOSTSTAG = AVOIDSELECTEDHOSTSTAG - 1
  integer, parameter          :: CHECKSELECTEDHOSTSTAG = CHECKREVIVEDHOSTSTAG - 1
  integer, parameter          :: CLEANMASTERDBTAG = CHECKSELECTEDHOSTSTAG - 1
  integer, parameter          :: DUMPDBTAG = CLEANMASTERDBTAG - 1
  integer, parameter          :: DUMPMASTERSDBTAG = DUMPDBTAG - 1
  integer, parameter          :: DUMPHOSTSDBTAG = DUMPMASTERSDBTAG - 1
  integer, parameter          :: FREEANYHOSTSTAG = DUMPHOSTSDBTAG - 1
  integer, parameter          :: FREEHOSTSTAG = FREEANYHOSTSTAG - 1
  integer, parameter          :: KILLMASTERSTAG = FREEHOSTSTAG - 1
  integer, parameter          :: SUICIDETAG = KILLMASTERSTAG - 1
  integer, parameter          :: SWITCHDUMPFILETAG = SUICIDETAG - 1
  integer, parameter          :: TURNREVIVALSONTAG = SWITCHDUMPFILETAG - 1
  integer, parameter          :: TURNREVIVALSOFFTAG = TURNREVIVALSONTAG - 1

  ! These are special tid values
  integer, parameter          :: UNASSIGNED = -1
  integer, parameter          :: AWAITINGREVIVAL = UNASSIGNED - 1
  integer, parameter          :: AWAITINGTHANKS = AWAITINGREVIVAL - 1

  integer, parameter          :: DUMPUNIT = LIST_UNIT + 1
  integer, parameter          :: TEMPUNIT = DUMPUNIT + 1
  integer, parameter          :: HDBUNIT = TEMPUNIT + 1
  integer, parameter          :: MDBUNIT = HDBUNIT + 1
!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it
! LF95.Linux/test [options] [input files]

  ! Our data type for the master tasks we'll be communicating with via pvm
  type master_T
    character(len=MACHINENAMELEN) :: name = ' '
    character(len=16            ) :: date = ' '
    integer                       :: tid = 0
    integer                       :: numChunks = 0
    integer                       :: numHosts = 0
    integer                       :: numFreed = 0
    integer, dimension(:), pointer:: hosts => null() ! hosts assigned to this master
    logical                       :: needs_host = .false.
    logical                       :: owes_thanks = .false.
    logical                       :: finished = .false.
  end type master_T

  !
  MLSMessageConfig%useToolkit = .false.
  MLSMessageConfig%logFileUnit = -1
  MLSMessageConfig%CrashOnAnyError = .true.
  time_config%use_wall_clock = .true.
  parallel%master = .true.  ! not merely master, but master of masters
  parallel%slaveFilename = 'pvm' ! for later cures only
  ! call get_options
  print *, 'tid of master you wish to quit'
  read *, tid
  call PVMFInitSend ( PvmDataDefault, bufferID )
  call PVMF90Pack ( GROUPNAME, info )
  ! By default we won't to quit if an error occurs
  ! (unless you want us to)
  call PVMFSend ( tid, GiveUpTag, info )
  call output('Already-running tellMasterToQuit tid: ', advance='no')
  call timestamp(tid, advance='yes')
  call timestamp(' commanded to quit', advance='yes')
end program tellMasterToQuit

! $Log$
! Revision 1.1  2010/04/13 20:30:46  pwagner
! First commit
!
