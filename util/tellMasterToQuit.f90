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
  use HighOutput, only: Timestamp
  use L2parinfo, only: Parallel, &
    & Giveuptag, Machinenamelen
  use Machine ! At Least Hp For Command Lines, And Maybe Getarg, Too
  use MLSMessagemodule, only: MLSMessageconfig
  use Output_M, only: Output
  use Printit_M, only: Set_Config
  use Pvm, only: PvmDatadefault, Pvmfinitsend, Pvmf90pack, &
    & Pvmfsend
  use Time_M, only: Time_Config

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
  integer :: INFO
  integer :: TID                   ! TID of master

  character(len=*), parameter :: GROUPNAME = "mlsl2"
!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

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
  call set_config ( useToolkit = .false., logFileUnit = -1 )
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
! Revision 1.5  2014/04/29 23:19:15  pwagner
! Fixed bug in last commit
!
! Revision 1.4  2014/04/29 23:15:16  pwagner
! Builds successfully with highOutput module
!
! Revision 1.3  2013/08/23 02:51:48  vsnyder
! Move PrintItOut to PrintIt_m
!
! Revision 1.2  2013/08/12 23:50:59  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 1.1  2010/04/13 20:30:46  pwagner
! First commit
!
