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

  public :: LaunchFWMSlaves, L2FWMSlaveTask, SetupFWMSlaves, TriggerSlaveRun
  public :: RequestSlavesOutput, ReceiveSlavesOutput

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
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
  
  ! ------------------------------------------- LaunchFwmSlaves ----
  subroutine LaunchFWMSlaves ( chunk )
    use Output_m, only: Output
    use Allocate_Deallocate, only: Allocate_Test, Deallocate_test
    use Machine, only: SHELL_COMMAND
    use MLSCommon, only: MLSChunk_T
    use L2ParInfo, only: PARALLEL, GETMACHINENAMES, MACHINENAMELEN, &
      & SLAVEARGUMENTS, SIG_REGISTER, INFOTAG
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
      commandLine = &
        & trim(parallel%submit) // ' ' // &
        & trim(parallel%executable) // ' ' // &
        & trim(slaveArguments)
      do maf = 1, noMAFs
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
    end do contactLoop

    if ( index ( switches, 'mas' ) /= 0 ) &
      & call output ( 'All slaves ready', advance='yes' )

    call Deallocate_test ( machineNames, 'MachineNames', ModuleName )
    call Deallocate_test ( heardFromSlave, 'heardFromSlave', ModuleName )

  end subroutine LaunchFWMSlaves

  ! ------------------------------------------------ L2FWMSlaveTask -----
  subroutine L2FWMSlaveTask ( mifGeolocation )
    ! This is the core routine for the 'slave mode' of the L2Fwm parallel stuff
    use Allocate_Deallocate, only: ALLOCATE_TEST
    use QuantityTemplates, only: QUANTITYTEMPLATE_T, DESTROYQUANTITYTEMPLATEDATABASE, &
      & INFLATEQUANTITYTEMPLATEDATABASE
    use VectorsModule, only: VECTORTEMPLATE_T, VECTOR_T, DESTROYVECTORDATABASE, &
      & DESTROYVECTORTEMPLATEDATABASE, CREATEVECTOR, CREATEMASK, &
      & CONSTRUCTVECTORTEMPLATE
    use ForwardModelConfig, only: FORWARDMODELCONFIG_T, DESTROYFWMCONFIGDATABASE, &
      & PVMUNPACKFWMCONFIG
    use ForwardModelIntermediate, only: FORWARDMODELINTERMEDIATE_T, FORWARDMODELSTATUS_T
    use L2ParInfo, only: PARALLEL, SIG_FINISHED, SIG_NEWSETUP, SIG_RUNMAF, INFOTAG, &
      & SIG_SENDRESULTS
    use MorePVM, only: PVMUNPACKSTRINGINDEX
    use PVM, only: PVMFRECV, PVMERRORMESSAGE, PVMDATADEFAULT, PVMFINITSEND, &
      & PVMF90UNPACK
    use PVMIDL, only: PVMIDLUNPACK, PVMIDLPACK
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Allocate
    use MatrixModule_0, only: M_Absent, M_Banded, M_Column_Sparse, M_Full, MatrixElement_T
    use MatrixModule_1, only: Matrix_T, CREATEEMPTYMATRIX, CLEARMATRIX
    use QuantityPVM, only: PVMRECEIVEQUANTITY
    use ForwardModelWrappers, only: FORWARDMODEL
    ! Dummy argument
    type (QuantityTemplate_T), dimension(:), pointer :: mifGeolocation

    ! Local variables
    integer :: INFO                     ! From pvm
    integer :: BUFFERID                 ! For pvm
    integer :: SIGNAL                   ! Signal code from master
    integer :: NOFWMCONFIGS             ! Number of forward model confis
    integer :: NOQUANTITIES             ! Number of quantity templates we'll need
    integer :: NOQUANTITIESINVECTOR     ! Number of quantities in the vector
    integer :: I,J                      ! Loop counters
    integer :: NAME                     ! An enumerated name
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

      ! We basically wait for some communication from the master task.
      call PVMFrecv ( parallel%masterTid, InfoTag, bufferID )
      if ( bufferID <= 0 ) &
        & call PVMErrorMessage ( bufferID, 'Receiveing information from master' )
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
        call DestroyFWMConfigDatabase ( fwmConfigs, deep=.true. )
        call DestroyVectorDatabase ( vectors )
        call DestroyVectorTemplateDatabase ( vectorTemplates )
        call DestroyQuantityTemplateDatabase ( quantities, &
          & ignoreMinorFrame=.true. )

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
                if ( info /= 0 ) call PVMErrorMessage ( info, 'vector values' )
              end if
              call PVMIDLUnpack ( flag, info )
              if ( info /= 0 ) call PVMErrorMessage ( info, 'vector mask flag' )
              if ( flag ) then
                call CreateMask ( vectors(i)%quantities(j) )
                call PVMIDLUnpack ( vectors(i)%quantities(j)%mask, info )
                if ( info /= 0 ) call PVMErrorMessage ( info, 'vector mask' )
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
          & l2(1), l2(2) )
        call Allocate_test ( fmStat%rows, jacobian%row%nb, 'fmStat%rows', ModuleName )

        ! OK, I think we're ready to go

      case ( SIG_RunMAF )
        ! Get the state vector for this iteration
        do j = 1, size ( vectors(fwmIn)%quantities )
          call PVMIDLUnpack ( vectors(fwmIn)%quantities(j)%values, info )
          if ( info /= 0 ) call PVMErrorMessage ( info, 'vector values' )
          call PVMIDLUnpack ( flag, info )
          if ( info /= 0 ) call PVMErrorMessage ( info, 'vector mask flag' )
          if ( flag ) then
            call CreateMask ( vectors(fwmIn)%quantities(j) )
            call PVMIDLUnpack ( vectors(fwmIn)%quantities(j)%mask, info )
            if ( info /= 0 ) call PVMErrorMessage ( info, 'vector mask' )
          end if
        end do
        ! Setup the fmStat stuff
        call PVMIDLUnpack ( fmStat%maf, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, 'fmStat%maf' )

        ! Loop over the configs and call the forward model
        do i = 1, size ( fwmConfigs )
          call ForwardModel ( fwmConfigs(i), &
            & vectors(fwmIn), vectors(fwmExtra), vectors(fwmOut), fmw, fmStat, jacobian )
        end do

        ! Now we sit patiently and wait for an instruction to 'dump' our results
        call PVMFrecv ( parallel%masterTid, infoTag, bufferID )
        if ( bufferID <= 0 ) &
          & call PVMErrorMessage ( bufferID, 'Receiveing go-ahead from master' )
        call PVMF90Unpack ( signal, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, 'Unpacking signal from master' )
        if ( signal == SIG_SendResults ) then
          call PVMFInitSend ( PVMDataDefault, bufferID )
          if ( bufferID <= 0 ) &
            & call PVMErrorMessage ( bufferID, 'Setting up results buffer' )
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
          ! OK, now send the results off
          call PVMFSend ( parallel%masterTid, InfoTag, info )
          if ( info /= 0 ) call PVMErrorMessage ( info, 'Sending results' )
        else if ( signal == SIG_Finished ) then
          finished = .true.
        else
          call MLSMessage ( MLSMSG_Error, ModuleName, 'Got unexpected message from master' )
        end if

        ! OK, we've run our forward models and sent our results
        call ClearMatrix ( jacobian )
      end select
      if ( finished ) exit mainLoop
    end do mainLoop
  end subroutine L2FWMSlaveTask

  ! ------------------------------------------------ ReceiveSlavesOutput ---
  subroutine ReceiveSlavesOutput ( outVector, fmStat, jacobian )
    ! The master uses this routine to get the output from each forward model
    ! slave.
    use VectorsModule, only: VECTOR_T
    use ForwardModelIntermediate, only: FORWARDMODELSTATUS_T
    use MatrixModule_1, only: MATRIX_T
    type (Vector_T), intent(inout) :: OUTVECTOR
    type (ForwardModelStatus_T), intent(inout) :: FMSTAT
    type (Matrix_T), intent(inout) :: JACOBIAN
  end subroutine ReceiveSlavesOutput

  ! ----------------------------------------------- RequestSlavesOutput ---
  subroutine RequestSlavesOutput ( maf )
    ! The master uses this routine to ask a slave to pack its output up
    use PVM, only: PVMFINITSEND, PVMF90PACK, PVMFSEND, PVMDATADEFAULT, &
      & PVMERRORMESSAGE
    use L2ParInfo, only: PARALLEL, SIG_SENDRESULTS, INFOTAG
    integer, intent(in) :: MAF
    ! Local variables
    integer :: INFO
    integer :: BUFFERID
    ! Executable code
    call PVMFInitSend ( PVMDataDefault, bufferID )
    if ( bufferID <= 0 ) call PVMErrorMessage ( info, 'Setting up output request' )
    call PVMF90Pack ( SIG_SendResults, info )
    if ( bufferID <= 0 ) call PVMErrorMessage ( info, 'Packing output request signal' )
    call PVMFSend ( slaveTIDs ( maf ), InfoTag, info )
    if ( bufferID <= 0 ) call PVMErrorMessage ( info, 'Sending output request' )
  end subroutine RequestSlavesOutput

  ! ------------------------------------------------ SetupFWMSlaves ----
  subroutine SetupFWMSlaves ( configs, inVector, extraVector, outVector )
    ! The master uses this routine to send the core information on the
    ! state vector layout etc. to each forward model slave.
    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use ForwardModelConfig, only: FORWARDMODELCONFIG_T, PVMPACKFWMCONFIG
    use VectorsModule, only: VECTOR_T
    use PVM, only: PVMFINITSEND, PVMDATADEFAULT, PVMERRORMESSAGE, PVMF90PACK, &
      & PVMFBCAST
    use PVMIDL, only: PVMIDLPACK, PVMIDLUNPACK
    use L2ParInfo, only: SIG_NEWSETUP, FWMSLAVEGROUP, INFOTAG
    use QuantityPVM, only: PVMSENDQUANTITY
    use MorePVM, only: PVMPACKSTRINGINDEX
    type (ForwardModelConfig_T), dimension(:), intent(in) :: CONFIGS
    type (Vector_T), target, intent(in) :: INVECTOR
    type (Vector_T), target, intent(in) :: EXTRAVECTOR
    type (Vector_T), target, intent(in) :: OUTVECTOR

    ! Local variables
    integer :: BUFFERID                 ! For PVM
    integer :: INFO                     ! Flag from PVM
    integer :: NOQUANTITIES             ! Number of quantities that are relevant
    integer :: I,J                      ! Loop counters
    integer, dimension(:), pointer :: VECINDS
    integer, dimension(:), pointer :: QTYINDS
    type (Vector_T), pointer :: V

    ! Executable code
    call PVMFInitSend ( PVMDataDefault, bufferID )
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
      call PVMIDLUnpack ( v%template%quantities, info )
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
            call PVMPack ( v%quantities(j)%values, info )
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

    ! That's it, let's get this information sent off
    call PVMFBcast ( FWMSlaveGroup, InfoTag, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Broadcasting setup information' )
      
  end subroutine SetupFWMSlaves

  ! ------------------------------------------------ TriggerSlaveRun ---
  subroutine TriggerSlaveRun ( state, maf )
    ! This routine is used by the master to launch one run
    use VectorsModule, only: VECTOR_T
    use PVM, only: PVMFINITSEND, PVMFSEND, PVMDATADEFAULT, PVMF90PACK
    use L2ParInfo, only: SIG_RUNMAF, INFOTAG
    type (Vector_T), intent(in) :: STATE
    integer, intent(in) :: MAF
    ! Local variables
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! Flag from PVM
    integer :: J                        ! Loop counter

    ! Executable code
    call PVMFInitSend ( PVMDataDefault, bufferID )
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
    call PVMFSend ( slaveTids ( maf ), InfoTag, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Sending trigger packet' )
  end subroutine TriggerSlaveRun

end module L2FWMParallel

! $Log$
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
