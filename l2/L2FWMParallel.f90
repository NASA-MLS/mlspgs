! Copyright (c) 2005, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contracts NAS7-1407/NAS7-03001 is acknowledged.

module L2FWMParallel
  ! This module is an alternative approach to parallel processing in
  ! mlsl2.  The idea is that the user can run a single chunk with
  ! different slave tasks computing the forward model for different mafs

  ! In general as this is for online studies, it does not need to be as robust
  ! as the main parallel method, so I've put in less stuff to check for tasks
  ! dying etc.

  implicit none
  private

  public :: LaunchFWMSlaves, L2FWMSlaveTask, SetupFWMSlaves, TriggerSlaveRun
  public :: RequestSlavesOutput, ReceiveSlavesOutput

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

  ! Local parameters
  ! These are the three vectors / templates we're after
  integer, parameter :: FWMIN = 1
  integer, parameter :: FWMEXTRA = FWMIn + 1
  integer, parameter :: FWMOUT = FWMExtra + 1
  integer, parameter :: NOVECTORS = 3

  ! We keep a record of the slaves
  integer, pointer, dimension(:), save :: slaveTIDs => NULL()

  ! Need a global flag here as the slave is called many times
  logical, save :: finished = .false.

contains
  
  ! --------------------------------------------  LaunchFwmSlaves  -----
  subroutine LaunchFWMSlaves ( Chunk )
    use Allocate_Deallocate, only: Allocate_Test, Deallocate_test
    use Chunks_m, only: MLSChunk_T
    use L2ParInfo, only: PARALLEL, GETMACHINENAMES, MACHINENAMELEN, &
      & SLAVEARGUMENTS, SIG_REGISTER, NOTIFYTAG, GETNICETIDSTRING
    use Machine, only: SHELL_COMMAND
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, PVMERRORMESSAGE
    use MLSSets, only: FINDFIRST
    use Output_m, only: Output
    use PVM, only: INFOTAG, MYPVMSPAWN, PVMFCATCHOUT, &
      & PVMFBUFINFO, PVMF90UNPACK, PVMFINITSEND, PVMFSEND, PVMF90PACK, &
      & PVMTASKHOST, PVMTASKEXIT, PVMRAW
    use Toggles, only: SWITCHES
    ! Dummy arguments
    type (MLSChunk_T), intent(in) :: CHUNK ! The chunk we're processing

    ! Local variables
    integer :: BUFFERID                 ! From PVM
    integer :: BYTES                    ! Used for call to bufinfo
    integer :: INFO                     ! Flag from PVM
    integer :: MAF                      ! Loop counter
    integer :: MACHINEIND                  ! Loop counter
    integer :: MSGTAG                   ! Message tag
    integer :: NOMACHINES               ! Number of machines
    integer :: NOMAFS                   ! Number of MAFs to process
    integer :: NOTASKS                  ! min ( noMachines, noMAFs )
    integer :: SIGNAL                   ! Flag from slave
    integer :: SLAVETID                 ! ID of one slave

    logical :: USINGSUBMIT              ! Set if using a submit mechanism

    integer, dimension(1) :: TID1       ! For MyPVMSpawn
    logical, pointer, dimension(:) :: HEARDFROMSLAVE ! Which slaves have we heard from

    character(len=2048) :: COMMANDLINE
    character(len=MachineNameLen), dimension(:), pointer :: MACHINENAMES

    ! Executable code
    usingSubmit = trim ( parallel%submit ) /= ''
    if ( usingSubmit ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Cannot use submit in fwmParallel mode' )
    noMAFs = chunk%lastMAFIndex - chunk%firstMAFIndex + 1

    ! Work out the information on our virtual machine
    nullify ( machineNames )
    call GetMachineNames ( machineNames )
    noMachines = size(machineNames)
    nullify ( slaveTids, heardFromSlave )
    noTasks = min ( noMAFs, noMachines )
    parallel%noFWMSlaves = noTasks
    call Allocate_test ( slaveTids, noMAFs, 'slaveTids', ModuleName )
    call Allocate_test ( heardFromSlave, noTasks, 'heardFromSlave', ModuleName )
    
    if ( index ( switches, 'mas' ) /= 0 ) &
      & call output ( 'Launching forward model slaves', advance='yes' )
    ! Now we're going to launch the slaves
    commandLine = trim ( parallel%executable )
    if ( index(switches,'slv') /= 0 ) then
      call PVMFCatchOut ( 1, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "calling catchout" )
    end if
    do machineInd = 1, noTasks
      info = myPVMSpawn ( trim(commandLine), PvmTaskHost, &
        trim(machineNames(machineInd)), 1, tid1 )
      ! Did this launch work
      if ( info /= 1 ) then
        call PVMErrorMessage ( tid1(1), &
          & 'Unable to launch fwmSlave on '//trim(machineNames(machineInd)) )
      end if
      slaveTids ( machineInd ) = tid1(1)
    end do
    do maf = noTasks+1, noMAFs
      slaveTids ( maf ) = slaveTids ( mod ( maf-1, noTasks ) + 1 )
    end do
    
    ! Now wait for them to get in touch, we behave differntly in the different
    ! modes (using submit or not)
    heardFromSlave = .false.
    if ( index ( switches, 'mas' ) /= 0 ) &
      & call output ( 'Waiting to hear from slaves', advance='yes' )
    contactLoop: do
      call IntelligentPVMFRecv ( -1, InfoTag, bufferID )
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
      machineInd = FindFirst ( slaveTids, slaveTid )
      if ( machineInd == 0 ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, 'Heard from an unknown forward model slave' )
      heardFromSlave ( machineInd ) = .true.
      if ( index ( switches, 'mas' ) /= 0 ) then
        call output ( 'Heard from ' )
        call output ( trim ( GetNiceTidString ( slaveTid ) ) )
        call output ( ', now heard from ' )
        call output ( count ( heardFromSlave ) )
        call output ( ' / ' )
        call output ( noTasks, advance='yes' )
      end if
      call PVMFNotify ( PVMTaskExit, NotifyTag, 1, (/ slaveTids(machineInd) /), info )
      if ( all ( heardFromSlave ) ) exit contactLoop
    end do contactLoop

    if ( index ( switches, 'mas' ) /= 0 ) &
      & call output ( 'All slaves started', advance='yes' )

    call Deallocate_test ( machineNames, 'MachineNames', ModuleName )
    call Deallocate_test ( heardFromSlave, 'heardFromSlave', ModuleName )

  end subroutine LaunchFWMSlaves

  ! ------------------------------------------------ L2FWMSlaveTask -----
  subroutine L2FWMSlaveTask ( mifGeolocation )
    ! This is the core routine for the 'slave mode' of the L2Fwm parallel stuff
    use Dump_0, only: Dump
    use Allocate_Deallocate, only: ALLOCATE_TEST
    use QuantityTemplates, only: QUANTITYTEMPLATE_T, DESTROYQUANTITYTEMPLATEDATABASE, &
      & INFLATEQUANTITYTEMPLATEDATABASE
    use VectorsModule, only: VECTORTEMPLATE_T, VECTOR_T, DESTROYVECTORDATABASE, &
      & DESTROYVECTORTEMPLATEDATABASE, CREATEVECTOR, CREATEMASK, &
      & CONSTRUCTVECTORTEMPLATE
    use ForwardModelConfig, only: FORWARDMODELCONFIG_T, DESTROYFWMCONFIGDATABASE, &
      & PVMUNPACKFWMCONFIG
    use ForwardModelIntermediate, only: FORWARDMODELINTERMEDIATE_T, FORWARDMODELSTATUS_T
    use L2ParInfo, only: PARALLEL, SIG_FINISHED, SIG_NEWSETUP, SIG_RUNMAF, &
      & SIG_SENDRESULTS, NOTIFYTAG
    use MorePVM, only: PVMUNPACKSTRINGINDEX
    use PVM, only: INFOTAG, PVMFINITSEND, &
      & PVMF90UNPACK, PVMRAW, PVMFFREEBUF
    use PVMIDL, only: PVMIDLUNPACK, PVMIDLPACK
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Allocate, &
      & PVMERRORMESSAGE
    use MatrixModule_0, only: M_Absent, M_Banded, M_Column_Sparse, M_Full, MatrixElement_T
    use MatrixModule_1, only: Matrix_T, CREATEEMPTYMATRIX, CLEARMATRIX, &
      & DESTROYMATRIX
    use QuantityPVM, only: PVMRECEIVEQUANTITY
    use ForwardModelWrappers, only: FORWARDMODEL
    use String_table, only: DISPLAY_STRING
    use Output_M, only: OUTPUT

    ! Dummy argument
    type (QuantityTemplate_T), dimension(:), pointer :: mifGeolocation

    ! Local variables
    integer :: I,J                      ! Loop counters
    integer :: INFO                     ! From pvm
    integer :: NAME                     ! An enumerated name
    integer :: NOFWMCONFIGS             ! Number of forward model confis
    integer :: NOQUANTITIES             ! Number of quantity templates we'll need
    integer :: NOQUANTITIESINVECTOR     ! Number of quantities in the vector
    integer :: RECVBUFFERID             ! For pvm
    integer :: SENDBUFFERID             ! For pvm
    integer :: SIGNAL                   ! Signal code from master
    logical :: FLAG                     ! A flag sent via pvm
    logical, dimension(2) :: L2         ! Two flags sent by pvm
    integer, dimension(:), pointer :: QTINDS ! Index of relevant quantities

    type (ForwardModelIntermediate_T) :: FMW
    type (ForwardModelStatus_T) :: FMSTAT
    type (QuantityTemplate_T), dimension(:), pointer :: QUANTITIES
    type (VectorTemplate_T), dimension(:), pointer :: VECTORTEMPLATES
    type (Vector_T), dimension(:), pointer :: VECTORS
    type (ForwardModelConfig_T), dimension(:), pointer :: FWMCONFIGS
    type (Matrix_T) :: JACOBIAN
    type (MatrixElement_T), pointer :: B ! A block from the jacobian matrix

    type (Vector_T), pointer :: V

    ! Executable code

    ! If the finished flag is set (by an earlier call to this routine) exit
    if ( finished ) return

    ! Some setup
    nullify ( quantities, vectorTemplates, vectors, fwmConfigs, qtInds )

    mainLoop: do
      ! We'll actually do a non-blocking check here to listen out for
      ! dead masters
      call IntelligentPVMFrecv ( parallel%masterTid, InfoTag, recvBufferId )

      ! We get from this just an integer (at least at first), based on which we do
      ! various tasks
      call PVMF90Unpack ( signal, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, 'Unpacking signal from master' )

      ! There are three possible signals
      ! SIG_NewSetup - a new set of state vectors and measurement vectors are to be used
      ! SIG_RunMAF - run a forward model for a given MAF and a given set of values of the
      !              state vector.
      ! SIG_Finished - That's all we require
      select case ( signal )
      case ( SIG_Finished )
        finished = .true.
      case ( SIG_NewSetup )
        ! For this we have to get the templates for the main state and
        ! measurement vectors, and the complete values for the auxilliary
        ! state vector.  We also get all the forward model configurations.
        ! First destroy all our old information
        call ClearSetup

        ! Get the forward model configs we'll need
        call PVMIDLUnpack ( noFWMConfigs, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, 'Unpacking noFwmConfigs' )
        ! Setup the structure for them
        allocate ( fwmConfigs ( noFwmConfigs ), STAT=info )
        if ( info /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_Allocate//'fwmConfigs' )
        do i = 1, noFwmConfigs
          call PVMUnpackFWMConfig ( fwmConfigs(i) )
        end do

        ! Now get all the quantity templates we'll need
        call PVMIDLUnpack ( noQuantities, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, 'Unpacking noQuantities' )
        info = InflateQuantityTemplateDatabase ( quantities, noQuantities )
        do i = 1, noQuantities
          ! Unpack a flag, if true this quantity is worth getting
          call PVMIDLUnpack ( flag, info )
          if ( info /= 0 ) call PVMErrorMessage ( info, 'Unpacking quantity relevant flag' )
          if ( flag ) then
            call PVMReceiveQuantity ( quantities(i), justUnpack=.true., &
              & mifGeolocation=mifGeolocation )
          end if
        end do

        ! Allocate the vector templates and vectors
        allocate ( vectorTemplates ( noVectors ), STAT=info )
        if ( info /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_Allocate//'vectorTemplates' )
        allocate ( vectors ( noVectors ), STAT=info )
        if ( info /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_Allocate//'vectors' )

        ! Now get the three vector templates and vectors
        do i = 1, noVectors
          v => vectors ( i )
          call PVMIDLUnpack ( noQuantitiesInVector, info )
          if ( info /= 0 ) call PVMErrorMessage ( info, 'noQuantitiesInVector' )
          call Allocate_Test ( qtInds, noQuantitiesInVector, 'qtInds', ModuleName )
          call PVMIDLUnpack ( qtInds, info )
          if ( info /= 0 ) call PVMErrorMessage ( info, 'qtInds' )
          call PVMUnpackStringIndex ( name, info )
          if ( info /= 0 ) call PVMErrorMessage ( info, 'vector template name' )
          call ConstructVectorTemplate ( name, quantities, qtInds, &
            & vectorTemplates ( i ) )
          ! Create the vector
          call PVMUnpackStringIndex ( name, info )
          if ( info /= 0 ) call PVMErrorMessage ( info, 'vector name' )
          vectors ( i ) = CreateVector ( name, vectorTemplates(i), &
            & quantities )
          ! Get the vector values in some circumstances
          if ( i == FWMExtra .or. i == FWMOut ) then
            do j = 1, noQuantitiesInVector
              if ( i == FWMExtra ) then
                call PVMIDLUnpack ( vectors(i)%quantities(j)%values, info )
                if ( info /= 0 ) call PVMErrorMessage ( info, 'extra/out vector values' )
              end if
              call PVMIDLUnpack ( flag, info )
              if ( info /= 0 ) call PVMErrorMessage ( info, 'extra/out vector mask flag' )
              if ( flag ) then
                call CreateMask ( vectors(i)%quantities(j) )
                call PVMIDLUnpack ( vectors(i)%quantities(j)%mask, info )
                if ( info /= 0 ) call PVMErrorMessage ( info, 'extra/out vector mask' )
              end if
            end do
          end if
        end do

        ! OK, we now have vectors and vector templates for fwmIn, fwmExtra and fwmOut
        ! We have also filled fwmExtra and got the forward model configs.
        ! Create a jacobian, get its information from the master
        call PVMIDLUnpack ( l2, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, '2 logical flags' )
        call CreateEmptyMatrix ( jacobian, 0, vectors(fwmOut), vectors(fwmIn), &
          & .not. l2(1), .not. l2(2) )
        call Allocate_test ( fmStat%rows, jacobian%row%nb, 'fmStat%rows', ModuleName )

        ! Free up this rather bulky receive buffer
        call PVMFFreeBuf ( recvBufferID, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, 'Freeing receive buffer' )

      case ( SIG_RunMAF )
        ! Get the state vector for this iteration
        do j = 1, size ( vectors(fwmIn)%quantities )
          call PVMIDLUnpack ( vectors(fwmIn)%quantities(j)%values, info )
          if ( info /= 0 ) call PVMErrorMessage ( info, 'state vector values' )
          call PVMIDLUnpack ( flag, info )
          if ( info /= 0 ) call PVMErrorMessage ( info, 'state vector mask flag' )
          if ( flag ) then
            call CreateMask ( vectors(fwmIn)%quantities(j) )
            call PVMIDLUnpack ( vectors(fwmIn)%quantities(j)%mask, info )
            if ( info /= 0 ) call PVMErrorMessage ( info, 'state vector mask' )
          end if
        end do
        ! Setup the fmStat stuff
        call PVMIDLUnpack ( fmStat%maf, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, 'fmStat%maf' )

        ! Loop over the configs and call the forward model
        fmStat%rows = .false.
        do i = 1, size ( fwmConfigs )
          call ForwardModel ( fwmConfigs(i), &
            & vectors(fwmIn), vectors(fwmExtra), vectors(fwmOut), fmw, fmStat, jacobian )
        end do

        ! Pack up our results in anticipation of sending them off
        call PVMFInitSend ( PVMRAW, sendBufferID )
        if ( sendBufferID <= 0 ) &
          & call PVMErrorMessage ( sendBufferID, 'Setting up results buffer' )
        call PVMIDLPack ( fmStat%rows, info ) 
        if ( info /= 0 ) call PVMErrorMessage ( info, 'fmStat%rows' )
        do i = 1, jacobian%row%nb
          ! Send corresponding values of fwmOut
          if ( fmStat%rows ( i ) ) then
            call PVMIDLPack ( vectors(fwmOut)%quantities ( &
              & jacobian%row%quant(i) ) % values ( :, &
              & jacobian%row%inst(i) ), info )
            if ( info /= 0 ) call PVMErrorMessage ( info, 'row of fwmOut' )
            ! Send non empty blocks of jacobian
            do j = 1, jacobian%col%nb
              b => jacobian%block ( i, j )
              call PVMIDLPack ( b%kind, info )
              if ( info /= 0 ) call PVMErrorMessage ( info, 'b%kind' )
              if ( b%kind == M_Banded .or. b%kind == M_Column_sparse ) then
                call PVMIDLPack ( size ( b%values, 1 ), info )
                if ( info /= 0 ) call PVMErrorMessage ( info, 'noValues for block' )
                call PVMIDLPack ( b%r1, info )
                if ( info /= 0 ) call PVMErrorMessage ( info, 'b%r1' )
                call PVMIDLPack ( b%r2, info )
                if ( info /= 0 ) call PVMErrorMessage ( info, 'b%r2' )
              end if
              if ( b%kind /= M_Absent ) then
                call PVMIDLPack ( b%values, info )
                if ( info /= 0 ) call PVMErrorMessage ( info, 'b%values' )
              end if
            end do
          end if                      ! Any rows here?
        end do

        ! OK, now wait patiently for the request to send the results off
        call IntelligentPVMFrecv ( parallel%masterTid, infoTag, recvBufferID )
        if ( recvBufferID <= 0 ) &
          & call PVMErrorMessage ( recvBufferID, 'Receiveing go-ahead from master' )
        call PVMF90Unpack ( signal, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, 'Unpacking signal from master' )
        if ( signal == SIG_SendResults ) then
          ! OK, now send the results off
          call PVMFSend ( parallel%masterTid, InfoTag, info )
          if ( info /= 0 ) call PVMErrorMessage ( info, 'Sending results' )
          ! OK, now free up the send buffer to save space
          call PVMFFreeBuf ( sendBufferID, info )
          if ( info /= 0 ) call PVMErrorMessage ( info, 'Freeing up send buffer' )
        else if ( signal == SIG_Finished ) then
          finished = .true.
        else
          call output ( 'Signal: ' )
          call output ( signal, advance='yes' )
          call MLSMessage ( MLSMSG_Error, ModuleName, 'Got unexpected message from master' )
        end if
        ! OK, we've run our forward models and sent our results
        call ClearMatrix ( jacobian )
      end select
      if ( finished ) exit mainLoop
    end do mainLoop

    ! Finished, tidy up
    call PVMFFreeBuf ( sendBufferID, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Freeing up send buffer' )
    call ClearSetup
  contains

    subroutine ClearSetup
      call DestroyMatrix ( jacobian )
      call DestroyFWMConfigDatabase ( fwmConfigs, deep=.true. )
      call DestroyVectorDatabase ( vectors )
      call DestroyVectorTemplateDatabase ( vectorTemplates )
      call DestroyQuantityTemplateDatabase ( quantities )
    end subroutine ClearSetup

  end subroutine L2FWMSlaveTask

  ! ------------------------------------------------ ReceiveSlavesOutput ---
  subroutine ReceiveSlavesOutput ( outVector, fmStat, jacobian )
    ! The master uses this routine to get the output from each forward model
    ! slave.
    use VectorsModule, only: VECTOR_T
    use ForwardModelIntermediate, only: FORWARDMODELSTATUS_T
    use MatrixModule_1, only: MATRIX_T
    use MLSMessageModule, only: PVMERRORMESSAGE
    use PVM, only: INFOTAG, PVMFFREEBUF
    use PVMIDL, only: PVMIDLUNPACK
    use L2ParInfo, only: GETNICETIDSTRING
    use MatrixModule_1, only: CREATEBLOCK
    use MatrixModule_0, only: M_ABSENT, M_BANDED, M_COLUMN_SPARSE, MATRIXELEMENT_T
    use Toggles, only: SWITCHES
    use Output_m, only: OUTPUT
    type (Vector_T), intent(inout) :: OUTVECTOR
    type (ForwardModelStatus_T), intent(inout) :: FMSTAT
    type (Matrix_T), intent(inout) :: JACOBIAN
    ! Local variables
    integer :: KIND                     ! Kind for block
    integer :: NOVALUES                 ! Number of values for banded/sparse
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! Flag from PVM
    integer :: I, J                     ! Loop counters
    type ( MatrixElement_T), pointer :: B ! A block from the jacobian

    ! Executable code
    if ( index ( switches, 'mas' ) /= 0 ) then
      call output ( 'Recieving results packet from ' )
      call output ( trim ( GetNiceTidString ( slaveTids(fmStat%maf) ) ) )
      call output ( ' MAF ' )
      call output ( fmStat%maf )
      call output ( ' ...' )
    end if
    call IntelligentPVMFrecv ( slaveTids ( fmStat%maf ), infoTag, bufferID )
    if ( bufferID <= 0 ) &
      & call PVMErrorMessage ( bufferID, 'Receiveing results from slave' )
    if ( index ( switches, 'mas' ) /= 0 ) &
      & call output ( ' Done', advance='yes' )
    call PVMIDLUnpack ( fmStat%rows, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Unpacking fmStat%rows' )
    do i = 1, jacobian%row%nb
      if ( fmStat%rows(i) ) then
        call PVMIDLUnpack ( outVector%quantities ( &
          & jacobian%row%quant(i) ) % values ( :, &
          & jacobian%row%inst(i) ), info )
        if ( info /= 0 ) call PVMErrorMessage ( info, 'unpacking values for fwmOut' )
        ! Get blocks of jacobian
        do j = 1, jacobian%col%nb
          call PVMIDLUnpack ( kind, info )
          if ( info /= 0 ) call PVMErrorMessage ( info, 'unpacking block kind' )
          if ( kind == M_Banded .or. kind == M_Column_sparse ) then
            call PVMIDLUnpack ( noValues, info )
            if ( info /= 0 ) call PVMErrorMessage ( info, 'unpacking noValues for block' )
          else
            noValues = 1                ! Not really used
          end if
          call CreateBlock ( jacobian, i, j, kind, noValues )
          b => jacobian%block ( i, j )
          if ( kind == M_Banded .or. kind == M_Column_sparse ) then
            call PVMIDLUnpack ( b%r1, info )
            if ( info /= 0 ) call PVMErrorMessage ( info, 'unpacking b%r1' )
            call PVMIDLUnpack ( b%r2, info )
            if ( info /= 0 ) call PVMErrorMessage ( info, 'unpacking b%r2' )
          end if
          if ( b%kind /= M_Absent ) then
            call PVMIDLUnpack ( b%values, info )
            if ( info /= 0 ) call PVMErrorMessage ( info, 'unpacking b%values' )
          end if
        end do
      end if
    end do
    ! Now free up our buffer to save space
    call PVMFFreeBuf ( bufferID, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'freeing receive buffer' )

  end subroutine ReceiveSlavesOutput

  ! ----------------------------------------------- RequestSlavesOutput ---
  subroutine RequestSlavesOutput ( maf )
    ! The master uses this routine to ask a slave to pack its output up
    use MLSMessageModule, only: PVMERRORMESSAGE
    use PVM, only: INFOTAG, PVMFINITSEND, PVMF90PACK, PVMFSEND, &
      & PVMRAW
    use L2ParInfo, only: PARALLEL, SIG_SENDRESULTS, GETNICETIDSTRING
    use Toggles, only: SWITCHES
    use Output_m, only: OUTPUT
    integer, intent(in) :: MAF
    ! Local variables
    integer :: INFO
    integer :: BUFFERID
    ! Executable code
    if ( index ( switches, 'mas' ) /= 0 ) then 
      call output ( 'Requesting output from ' )
      call output ( trim(GetNiceTidString ( slaveTids(maf) ) ) )
      call output ( ' MAF ' )
      call output ( maf, advance='yes' )
    end if
    call PVMFInitSend ( PVMRAW, bufferID )
    if ( bufferID <= 0 ) call PVMErrorMessage ( info, 'Setting up output request' )
    call PVMF90Pack ( SIG_SendResults, info )
    if ( bufferID <= 0 ) call PVMErrorMessage ( info, 'Packing output request signal' )
    call PVMFSend ( slaveTIDs ( maf ), InfoTag, info )
    if ( bufferID <= 0 ) call PVMErrorMessage ( info, 'Sending output request' )
  end subroutine RequestSlavesOutput

  ! ------------------------------------------------ SetupFWMSlaves ----
  subroutine SetupFWMSlaves ( configs, inVector, extraVector, outVector, jacobian )
    ! The master uses this routine to send the core information on the
    ! state vector layout etc. to each forward model slave.
    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use ForwardModelConfig, only: FORWARDMODELCONFIG_T, PVMPACKFWMCONFIG
    use VectorsModule, only: VECTOR_T
    use MLSMessageModule, only: PVMERRORMESSAGE
    use PVM, only: INFOTAG, PVMFINITSEND, PVMF90PACK, &
      & PVMFBCAST, PVMRAW, PVMFFREEBUF
    use PVMIDL, only: PVMIDLPACK, PVMIDLUNPACK
    use L2ParInfo, only: SIG_NEWSETUP, FWMSLAVEGROUP
    use QuantityPVM, only: PVMSENDQUANTITY
    use MorePVM, only: PVMPACKSTRINGINDEX
    use Toggles, only: SWITCHES
    use Output_m, only: OUTPUT
    use MatrixModule_1, only: MATRIX_T
    type (ForwardModelConfig_T), dimension(:), intent(in) :: CONFIGS
    type (Vector_T), target, intent(in) :: INVECTOR
    type (Vector_T), target, intent(in) :: EXTRAVECTOR
    type (Vector_T), target, intent(in) :: OUTVECTOR
    type (Matrix_T), intent(in) :: JACOBIAN

    ! Local variables
    integer :: BUFFERID                 ! For PVM
    integer :: INFO                     ! Flag from PVM
    integer :: NOQUANTITIES             ! Number of quantities that are relevant
    integer :: I,J                      ! Loop counters
    integer, dimension(:), pointer :: VECINDS
    integer, dimension(:), pointer :: QTYINDS
    type (Vector_T), pointer :: V

    ! Executable code
    if ( index ( switches, 'mas' ) /= 0 ) &
      & call output ( 'Setting up forward model slaves', advance='yes' )
    call PVMFInitSend ( PVMRaw, bufferID )
    if ( bufferID <= 0 ) &
      & call PVMErrorMessage ( bufferID, 'Setting up buffer for FWMSlaveSetup' )
    call PVMF90Pack ( SIG_NewSetup, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Packing signal for FWMSlaveSetup' )

    ! Now we send the forward model configs
    call PVMIDLPack ( size ( configs ), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Packing no fwm configs' )
    do i = 1, size ( configs )
      call PVMPackFWMConfig ( configs ( i ) )
    end do

    ! Now we work out how many quantities are relevant
    noQuantities = max ( &
      & maxval ( inVector%template%quantities ), &
      & maxval ( extraVector%template%quantities ), &
      & maxval ( outVector%template%quantities ) )
    nullify ( vecInds, qtyInds )
    call Allocate_test ( vecInds, noQuantities, 'vecInds', ModuleName )
    call Allocate_test ( qtyInds, noQuantities, 'qtyInds', ModuleName )
    vecInds = 0
    qtyInds = 0
    ! Note that some quantities may be in two vectors, it doesn't
    ! matter as long as the slave gets the template from somewhere
    do i = 1, noVectors
      select case ( i )
      case ( FWMIn )
        v => inVector
      case ( FWMExtra )
        v => extraVector
      case  ( FWMOut )
        v => outVector
      end select
      do j = 1, size ( v%quantities )
        vecInds ( v%template%quantities(j) ) = i
        qtyInds ( v%template%quantities(j) ) = j
      end do
    end do
    ! Now pack up the relevant ones.
    call PVMIDLpack ( noQuantities, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'packing noQuantities' )
    do j = 1, noQuantities
      call PVMIDLPack ( vecInds ( j ) /= 0 , info )
      if ( info /= 0 ) call PVMErrorMessage ( info, 'Packing quantity relevant flag' )
      if ( vecInds ( j ) /= 0 ) then
        select case ( vecInds(j) )
        case ( FWMIn )
          v => inVector
        case ( FWMExtra )
          v => extraVector
        case  ( FWMOut )
          v => outVector
        end select
        call PVMSendQuantity ( v%quantities(qtyInds(j)), justPack=.true., &
          & noValues=.true., noMask=.true., skipMIFGeolocation=.true. )
      end if
    end do

    ! Now send the three templates, and the vector values in some circumstances
    do i = 1, noVectors
      select case ( i )
      case ( FWMIn )
        v => inVector
      case ( FWMExtra )
        v => extraVector
      case  ( FWMOut )
        v => outVector
      end select
      ! Send the vector template information
      call PVMIDLPack ( size ( v%quantities ), info )
      if ( info /= 0 ) call PVMErrorMessage ( info, 'packing size(quantities)' )
      call PVMIDLPack ( v%template%quantities, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, 'packing template%quantities' )
      call PVMPackStringIndex ( v%template%name, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, 'packing vector template name' )
      ! Send the vector information
      call PVMPackStringIndex ( v%name, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, 'packing vector name' )
      ! Pack the vector values in some circumstances
      if ( i == FWMExtra .or. i == FWMOut ) then
        do j = 1, size ( v%quantities )
          if ( i == FWMExtra ) then
            call PVMIDLPack ( v%quantities(j)%values, info )
            if ( info /= 0 ) call PVMErrorMessage ( info, 'packing vector values' )
          end if
          call PVMIDLPack ( associated ( v%quantities(j)%mask ), &
            & info )
          if ( info /= 0 ) call PVMErrorMessage ( info, 'packing quantity mask flag' )
          if ( associated ( v%quantities(j)%mask ) ) then
            call PVMIDLPack ( v%quantities(j)%mask, info )
            if ( info /= 0 ) call PVMErrorMessage ( info, 'packing quantity mask' )
          end if                        ! Send mask?
        end do                          ! Loop over quantities
      end if                            ! Extra or output
    end do                              ! Loop over vectors

    ! Finally, just send the 'instance first' flags
    call PVMIDLPack ( (/ jacobian%row%instFirst, jacobian%col%instFirst /), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, &
      & 'packing jacobian instance first flags' )

    ! That's it, let's get this information sent off
    call PVMFBcast ( FWMSlaveGroup, InfoTag, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Broadcasting setup information' )

    ! Free up the send buffer
    call PVMFFreeBuf ( bufferID, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'freeing send buffer' )
      
  end subroutine SetupFWMSlaves

  ! ------------------------------------------------ TriggerSlaveRun ---
  subroutine TriggerSlaveRun ( state, maf )
    ! This routine is used by the master to launch one run
    use VectorsModule, only: VECTOR_T
    use MLSMessageModule, only: PVMERRORMESSAGE
    use PVM, only: INFOTAG, PVMFINITSEND, PVMFSEND, PVMRAW, PVMF90PACK
    use PVMIDL, only: PVMIDLPACK
    use L2ParInfo, only: SIG_RUNMAF, GETNICETIDSTRING
    use Toggles, only: SWITCHES
    use Output_m, only: OUTPUT
    type (Vector_T), intent(in) :: STATE
    integer, intent(in) :: MAF
    ! Local variables
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! Flag from PVM
    integer :: J                        ! Loop counter
    integer :: TASK                     ! Task index

    ! Executable code
    call PVMFInitSend ( PVMRaw, bufferID )
    if ( bufferID <= 0 ) call PVMErrorMessage ( bufferID, &
      & 'Setting up buffer for TriggerSlaveRun' )
    call PVMF90Pack ( SIG_RunMAF, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Packing sig_RunMAF' )
    ! Send the state vector values
    do j = 1, size ( state%quantities )
      call PVMIDLPack ( state%quantities(j)%values, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, &
        & 'Packing state quantity values' )
      call PVMIDLPack ( associated ( state%quantities(j)%mask ), info )
      if ( info /= 0 ) call PVMErrorMessage ( info, &
        & 'Packing state quantity mask flag' )
      if ( associated ( state%quantities(j)%mask ) ) then
        call PVMIDLPack ( state%quantities(j)%mask, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, &
          & 'Packing state quantity mask' )
      end if
    end do
    ! Send the fmStat stuff
    call PVMIDLPack ( maf, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'packing maf' )
    ! OK, send this off
    task = mod ( maf-1, size(slaveTids) ) + 1
    if ( index ( switches, 'mas' ) /= 0 ) then
      call output ( 'Triggering ' // trim ( GetNiceTidString(slaveTids(task)) ) )
      call output ( ' to do MAF ' )
      call output ( task, advance='yes' )
    end if
      
    call PVMFSend ( slaveTids ( task ), InfoTag, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Sending trigger packet' )
  end subroutine TriggerSlaveRun

  ! ================================================ Private procedures

  subroutine IntelligentPVMFRecv ( tid, tag, bufferID )
    use L2PARINFO, only: NOTIFYTAG
    use MLSMessageModule, only: MLSMSG_ERROR, MLSMESSAGE, PVMERRORMESSAGE
    use PVM, only: PVMFRECV
    ! Dummy arguments
    integer, intent(in) :: TID
    integer, intent(in) :: TAG
    integer, intent(out) :: BUFFERID
    ! Parameters etc.
    integer, parameter :: DELAY = 200000  ! For Usleep, no. microsecs
    external :: USLEEP

    ! Local variables
    listenLoop: do
      call PVMFNRecv ( -1, NotifyTag, bufferID )
      if ( bufferID > 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'A task in the pvm system died' )
      call PVMFNRecv ( tid, tag, bufferID )
      if ( bufferID > 0 ) exit listenLoop
      if ( bufferID < 0 ) call PVMErrorMessage ( bufferID, 'Listening for message' )
      call usleep ( delay )
    end do listenLoop
  end subroutine IntelligentPVMFRecv

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module L2FWMParallel

! $Log$
! Revision 2.19  2005/03/15 23:50:15  pwagner
! PVMERRORMESSAGE now part of MLSMessageModule
!
! Revision 2.18  2004/06/10 00:58:45  vsnyder
! Move FindFirst, FindNext from MLSCommon to MLSSets
!
! Revision 2.17  2004/05/19 19:16:11  vsnyder
! Move MLSChunk_t to Chunks_m
!
! Revision 2.16  2003/06/20 19:38:25  pwagner
! Allows direct writing of output products
!
! Revision 2.15  2003/01/13 20:59:02  livesey
! Added calls to PVMFFreeBuf
!
! Revision 2.14  2003/01/13 20:15:49  livesey
! Slight changes to slave launching, uses parllel%executable
!
! Revision 2.13  2002/12/11 02:14:08  livesey
! Removed extra dumps/outputs
!
! Revision 2.12  2002/12/11 02:08:17  livesey
! Bug fix with handling of fmStat%rows?
!
! Revision 2.11  2002/12/11 01:59:27  livesey
! Slightly new approach, also some extra diagnostics that will have to go soon
!
! Revision 2.10  2002/12/06 18:43:26  livesey
! New approach to sharing out the work load.  Doesn't require there to be
! a full complement of 'machines' one for each MAF
!
! Revision 2.9  2002/12/05 02:21:29  livesey
! Changes to improve performance (hopefully)
!
! Revision 2.8  2002/10/08 20:34:02  livesey
! Added notify stuff for fwm parallel
!
! Revision 2.7  2002/10/08 17:40:35  livesey
! Lots of debugging
!
! Revision 2.6  2002/10/08 17:36:21  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.5  2002/10/07 01:23:51  livesey
! OK, all the code written and compiles, but not tested.
!
! Revision 2.4  2002/10/06 22:22:20  livesey
! Nearly all routines filled out now, only one more to go
!
! Revision 2.3  2002/10/06 02:04:31  livesey
! More progress all routines at least stubbed out.
!
! Revision 2.2  2002/10/06 01:10:31  livesey
! More code added, still not complete.
!
! Revision 2.1  2002/10/05 02:06:28  livesey
! First version, should have checked in earlier
!
