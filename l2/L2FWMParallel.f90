! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module L2FWMParallel
  ! This module is an alternative approach to parallel processing in
  ! mlsl2.  The idea is that the user can run a single chunk with
  ! different slave tasks computing the forward model for different mafs

  ! In general as this is for online studies, it does not need to be as robust
  ! as the main parallel method, so I've put in less stuff to check for tasks
  ! dying etc.

  implicit none
  private

  public :: LaunchFWMSlaves, L2FWMSlaveTask

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  !---------------------------------------------------------------------------

  integer, pointer, dimension(:), private, save :: slaveTIDs => NULL()
  ! We keep a record of the slaves

contains
  
  ! ------------------------------------------- LaunchFwmSlaves ----
  subroutine LaunchFWMSlaves ( chunk )
    use Output_m, only: Output
    use Allocate_Deallocate, only: Allocate_Test, Deallocate_test
    use Machine, only: SHELL_COMMAND
    use MLSCommon, only: MLSChunk_T
    use L2ParInfo, only: PARALLEL, GETMACHINENAMES, MACHINENAMELEN, &
      & SLAVEARGUMENTS, SIG_REGISTER, INFOTAG, MAFTAG
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use Toggles, only: SWITCHES
    use PVM, only: PVMDATADEFAULT, MYPVMSPAWN, PVMFCATCHOUT, PVMERRORMESSAGE, &
      & PVMFRECV, PVMFBUFINFO, PVMF90UNPACK, PVMFINITSEND, PVMFSEND, PVMF90PACK
    ! Dummy arguments
    type (MLSChunk_T), intent(in) :: CHUNK ! The chunk we're processing

    ! Local variables
    integer :: BUFFERID                 ! From PVM
    integer :: BYTES                    ! Used for call to bufinfo
    integer :: INFO                     ! Flag from PVM
    integer :: MAF                      ! Loop counter
    integer :: MSGTAG                   ! Message tag
    integer :: NOMACHINES               ! Number of machines
    integer :: NOMAFS                   ! Number of MAFs to process
    integer :: PVMTASKHOST              ! For bufInfo
    integer :: SIGNAL                   ! Flag from slave
    integer :: SLAVETID                 ! ID of one slave

    logical :: USINGSUBMIT              ! Set if using a submit mechanism

    integer, dimension(1) :: TID1       ! For MyPVMSpawn
    logical, pointer, dimension(:) :: HEARDFROMSLAVE ! Which slaves have we heard from

    character(len=16) :: MAFNOSTR       ! MAF number as a string
    character(len=2048) :: COMMANDLINE
    character(len=MachineNameLen), dimension(:), pointer :: MACHINENAMES

    ! Executable code
    usingSubmit = trim ( parallel%submit ) /= ''
    noMAFs = chunk%lastMAFIndex - chunk%firstMAFIndex + 1

    ! Work out the information on our virtual machine
    if ( .not. usingSubmit ) then
      nullify ( machineNames )
      call GetMachineNames ( machineNames )
      noMachines = size(machineNames)
      if ( noMachines < noMAFs ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Too few machines for fwmParallel mode' )
    else
      noMachines = 0
    end if
    call Allocate_test ( slaveTids, noMAFS, 'slaveTids', ModuleName )
    nullify ( heardFromSlave )
    call Allocate_test ( heardFromSlave, noMAFS, 'heardFromSlave', ModuleName )
    
    if ( index ( switches, 'mas' ) /= 0 ) &
      & call output ( 'Launching forward model slaves', advance='yes' )
    ! Now we're going to launch the slaves
    if ( usingSubmit ) then ! --------------------- Using a batch system
      do maf = 1, noMAFs
        write ( mafNoStr, '(i0)' ) maf
        commandLine = &
          & trim(parallel%submit) // ' ' // &
          & trim(parallel%executable) // ' ' // &
          & ' --slaveMAF ' // trim(mafNoStr) // ' ' // &
          & trim(slaveArguments)
        call shell_command ( trim(commandLine) )
      end do
    else ! ----------------------------------------- Start job using pvmspawn
      commandLine = 'mlsl2'
      if ( index(switches,'slv') /= 0 ) then
        call PVMFCatchOut ( 1, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, "calling catchout" )
      end if
      do maf = 1, noMAFs
        info = myPVMSpawn ( trim(commandLine), PvmTaskHost, &
          trim(machineNames(maf)), 1, tid1 )
        slaveTids(maf) = tid1(1)
        ! Did this launch work
        if ( info /= 1 ) then
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Unable to launch fwmSlave on '//trim(machineNames(maf)) )
        end if
      end do
    end if

    ! Now wait for them to get in touch, we behave differntly in the different
    ! modes (using submit or not)
    if ( usingSubmit ) maf = 1
    heardFromSlave = .false.
    if ( index ( switches, 'mas' ) /= 0 ) &
      & call output ( 'Waiting to hear from slaves', advance='yes' )
    contactLoop: do
      call PVMFRecv ( -1, InfoTag, bufferID )
      ! Got a message who sent this
      call PVMFBufInfo ( bufferID, bytes, msgTag, slaveTid, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, "calling PVMFBufInfo" )
      call PVMF90Unpack ( signal, info )
      if ( info /= 0 ) then
        call PVMErrorMessage ( info, "unpacking signal" )
      endif
      if ( signal /= sig_register ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, 'Expected registration message from fwmSlave' )
      if ( usingSubmit ) then
        slaveTids ( maf ) = slaveTid
        maf = maf + 1
      end if
      ! Otherwise, we know who this is.
      heardFromSlave ( maf ) = .true.
      ! Send back the maf number as a confirmation.
      call PVMFinitSend ( PVMDataDefault, bufferID )
      call PVMF90Pack ( maf, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'packing maf' )
      call PVMFSend ( slaveTid, MAFTag, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'sending maf' )
      if ( all ( heardFromSlave ) ) exit contactLoop
    end do contactLoop

    if ( index ( switches, 'mas' ) /= 0 ) &
      & call output ( 'All slaves ready', advance='yes' )

    call Deallocate_test ( machineNames, 'MachineNames', ModuleName )
    call Deallocate_test ( heardFromSlave, 'heardFromSlave', ModuleName )

  end subroutine LaunchFWMSlaves

  ! ------------------------------------------------ L2FWMSlaveTask -----
  subroutine L2FWMSlaveTask ( mifGeolocation )
    use QuantityTemplates, only: QUANTITYTEMPLATE_T
    type (QuantityTemplate_T), dimension(:), pointer :: mifGeolocation
  end subroutine L2FWMSlaveTask

end module L2FWMParallel

! $Log$
! Revision 2.1  2002/10/05 02:06:28  livesey
! First version, should have checked in earlier
!
