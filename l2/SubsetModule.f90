! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module SubsetModule

  ! This module deals with all the user requests to 'subset' quantities
  ! through the mask field.  In the first instance, it was abstracted
  ! out of the retrieval module.

  implicit none
  private

  public :: RestrictRange, SetupSubset, SetupFlagCloud, UpdateMask

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------
  
  ! Local parameters
  character(len=32), private, parameter :: WRONGUNITS = 'The wrong units were supplied'

contains ! ========= Public Procedures ============================

  ! -------------------------------------------------- RestrictRange ---
  subroutine RestrictRange ( key, vectors )
    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use Intrinsic, only: PHYQ_DIMENSIONLESS
    use MLSCommon, only: R8
    use Expr_M, only: EXPR
    use VectorsModule, only: VECTOR_T, VECTORVALUE_T, GETVECTORQTYBYTEMPLATEINDEX, &
      & M_LINALG, GETVECTORQUANTITYBYTYPE, SETMASK, CREATEMASK
    use Init_Tables_Module, only: F_QUANTITY, F_PTANQUANTITY, F_BASISFRACTION, &
      & F_MINCHANNELS, F_SIGNALS, F_MEASUREMENTS, F_MASK
    use Init_Tables_Module, only: FIELD_FIRST, FIELD_LAST
    use Init_Tables_Module, only: L_ZETA, L_RADIANCE
    use Tree, only: NSONS, SUBTREE, DECORATION, SUB_ROSA
    use MoreTree, only: GET_FIELD_ID
    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR, MLSMSG_ALLOCATE, &
      & MLSMSG_DEALLOCATE
    use MLSSignals_M, only: Signal_T
    use Parse_Signal_M, only: EXPAND_SIGNAL_LIST
    use String_Table, only: GET_STRING
    use ManipulateVectorQuantities, only: FINDONECLOSESTINSTANCE
    use MLSNumerics, only: HUNT

    ! Dummy arguments
    integer, intent(in) :: KEY          ! Tree node
    type (Vector_T), dimension(:), intent(inout), target :: VECTORS

    ! Local variables

    integer :: CHANNEL                  ! Loop counter / channel index
    integer :: FIELD                    ! ID for field in l2cf line
    integer :: GSON                     ! Tree node
    integer :: INDEX                    ! Loop counter
    integer :: INDEX0                   ! Array index
    integer :: INDEX1                   ! Array index
    integer :: INSTANCE                 ! Instance index
    integer :: J                        ! Loop counter
    integer :: MAF                      ! Major frame index
    integer :: MASK                     ! Bits to mask
    integer :: MEASQTY                  ! Loop counter
    integer :: MINCHANNELS              ! Value of argument of same name
    integer :: NOGOODMIFS               ! How many MIFs meet our criteria for one chan
    integer :: NOINSTANCES              ! Number of instances
    integer :: NOMAFS                   ! Number of major frames
    integer :: NOMEASUREMENTS           ! Number of relevant measurements
    integer :: NOMIFS                   ! Number of minor frames
    integer :: NOVALIDCHANNELS          ! How many channels meet our criteria
    integer :: QUANTITYINDEX            ! Index in database, not vector
    integer :: SIDEBAND                 ! Returned by parse_signal
    integer :: SON                      ! Tree node
    integer :: STATUS                   ! Flag from allocate etc.
    integer :: SURF                     ! Loop counter
    integer :: VECTORINDEX              ! Vector index for quantity

    integer, dimension(2) :: EXPRUNIT   ! Units for expression
    integer, dimension(:), pointer :: MAFFORINSTANCE ! Maps instances to mafs
    integer, dimension(:), pointer :: MAFSFORINSTANCE ! Which maf relevant for each inst.
    integer, dimension(:), pointer :: MEASQTYINDS ! Indices
    integer, dimension(:), pointer :: MIFPOINTINGS ! Which basis for each mif?
    integer, dimension(:), pointer :: SIGNALINDS ! Returned by parse_signal

    logical :: ERRORFLAG                ! Set if problem
    logical :: FOUNDONE                 ! Flag for instance identification
    logical :: Got(field_first:field_last)   ! "Got this field already"

    character (len=132) :: SIGNALSTRING ! One signal
    character(len=1), dimension(:), pointer :: MIFMASKS ! Part of rad%mask

    real(r8) :: BASISFRACTION           ! Value of argument of same name
    real(r8), dimension(2) :: EXPRVALUE ! Vaue for expression
    real(r8), dimension(:), pointer :: SURFS ! Surfaces for quantity

    type (Signal_T), dimension(:), pointer :: SIGNALS
    type (VectorValue_T), pointer :: PTAN ! The ptan quantity
    type (VectorValue_T), pointer :: QTY ! The quantity to subset
    type (VectorValue_T), pointer :: RAD ! A radiance quantity
    type (Vector_T), pointer :: MEASUREMENTS ! Measurement vector

    ! Executable code
    nullify ( signals )
    ! Setup defaults
    basisFraction = 0.5
    minChannels = 1
    mask = m_linAlg
    got = .false.
    ! Loop over arguments to this l2cf command
    do j = 2, nsons ( key )
      son = subtree ( j, key )
      field = get_field_id ( son )
      if ( nsons ( son )  > 1 ) gson = subtree ( 2, son )
      select case ( field )
      case ( f_quantity )
        vectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        qty => GetVectorQtyByTemplateIndex ( vectors(vectorIndex), quantityIndex )
      case ( f_ptanQuantity )
        vectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        ptan => GetVectorQtyByTemplateIndex ( vectors(vectorIndex), quantityIndex )
      case ( f_minChannels )
        call expr ( gson, exprUnit, exprValue )
        if ( exprUnit(1) /= PHYQ_Dimensionless ) call AnnounceError ( key, &
          WrongUnits, f_minChannels )
        minChannels = nint ( exprValue(1) )
      case ( f_basisFraction )
        call expr ( gson, exprUnit, exprValue )
        if ( exprUnit(1) /= PHYQ_Dimensionless ) call AnnounceError ( key, &
          WrongUnits, f_basisFraction )
        basisFraction = exprValue(1)
      case ( f_measurements )
        vectorIndex = decoration(decoration(gson))
        measurements => vectors ( vectorIndex )
      case ( f_signals )
        call Expand_Signal_List ( son, signals, errorFlag )
        if ( errorFlag ) call AnnounceError ( key, 'Bad signal spec.' )
      case ( f_mask )
        mask = GetMaskBit ( son )
      end select
    end do
    ! The parser has ensured that we 'got' all the fields we required,
    ! and we supplied defaults for the rest, so we're ready to go now.
    ! Some quick error checking
    if ( qty%template%verticalCoordinate /= l_zeta ) &
      & call AnnounceError ( key, 'Quantity is not on log pressure surfaces' )
    if ( .not. associated ( qty%mask ) ) call CreateMask ( qty )

    ! Identify all the measurement quantities we'll be using
    noMeasurements = size(signals)
    nullify ( measQtyInds )
    call Allocate_test ( measQtyInds, noMeasurements, 'measQtyInds', ModuleName )
    do measQty = 1, noMeasurements
      rad => GetVectorQuantityByType ( measurements, quantityType=l_radiance, &
        & signal=signals(measQty)%index, sideband=signals(measQty)%sideband )
      if ( rad%template%instrumentModule /= ptan%template%instrumentModule ) &
        & call AnnounceError ( key, 'Ptan/radiance not from the same module' )
      measQtyInds(measQty) = rad%index
    end do

    ! Now identify the relevant MAFs for each profile of our quantity
    ! This kind of involves running FindOneClosestInstance backwards
    noInstances = qty%template%noInstances
    noMAFS = ptan%template%noInstances
    nullify ( mafsForInstance )
    call Allocate_test ( mafsForInstance, noInstances, 'mafsForInstance', ModuleName )
    mafsForInstance = 0
    do maf = 1, noMAFs
      instance = FindOneClosestInstance ( qty, ptan, maf )
      mafsForInstance ( instance ) = maf
    enddo
    ! Deal with any we missed (probably at the ends)
    foundOne = .false.
    do instance = 1, noInstances
      if ( mafsForInstance ( instance ) /= 0 ) then
        foundOne = .true.
      else
        if ( foundOne ) then
          ! OK, make this gap contain the previous value,
          ! we know instance>1 here.
          mafsForInstance ( instance ) = mafsForInstance ( instance - 1 )
        else
          ! OK, still not seen on, make this use the first maf
          mafsForInstance ( instance ) = 1
        end if
      end if
    end do

    ! Setup some more arrrays
    nullify ( mifPointings )
    call Allocate_test ( mifPointings, ptan%template%noSurfs, &
      & 'mifPointings', ModuleName )

    ! Now we get down to the main work, we'll go through this instance by instance
    do instance = 1, noInstances
      maf = MAFsForInstance ( instance )
      if ( qty%template%coherent ) then
        surfs => qty%template%surfs ( :, 1 )
      else
        surfs => qty%template%surfs ( :, instance )
      end if
      ! Work out where each MIF points
      call Hunt ( surfs, ptan%values(:,maf), mifPointings, allowTopValue=.true. )
      ! Go through each surface/channel in the quantity and work out whether
      ! We're going to use it.
      surfLoop: do surf = 1, qty%template%noSurfs
        index0 = ( surf-1 ) * qty%template%noChans + 1
        index1 = index0 + qty%template%noChans - 1
        ! If were retrieving data for any 'channel' at this surface
        ! think about changing our minds, otherwise don't bother.
        if ( any ( iand ( ichar(qty%mask(index0:index1,instance)), &
          & m_linAlg ) == 0 ) ) then
          ! How many MIFs total above this surface?
          noMIFs = count ( mifPointings == surf )
          ! Now, how many channels have enough good channels (i.e. enough channels
          ! with more valid MIFs in this basis range than goodMIFs*basisFraction)
          noValidChannels = 0
          do measQty = 1, noMeasurements
            rad => measurements%quantities ( measQtyInds ( measQty ) )
            if ( associated ( rad%mask ) ) then
              channelLoop: do channel = 1, rad%template%noChans
                if ( .not. signals(measQty)%channels(channel) ) cycle channelLoop
                ! Point to the relevant masks for this radiance channel only
                mifMasks => rad%mask ( channel : rad%template%instanceLen : &
                  & rad%template%noChans, instance )
                noGoodMIFs = count ( mifPointings == surf .and. &
                  & iand( ichar(mifMasks), m_linAlg ) == 0 )
                if ( (1.0*noGoodMIFs)/noMIFs > basisFraction ) &
                  & noValidChannels = noValidChannels + 1
              end do channelLoop
            else
              noValidChannels = noValidChannels + 1
            end if
          end do ! Signal loop
          if ( noValidChannels < minChannels ) then
            ! We don't have enough information to retrieve this value,
            ! So mark it as not to be retrieved (or whatever bits have been supplied)
            do index = index0, index1
              call SetMask ( qty%mask(:,instance), (/ index /), what=mask )
            end do
          else
            ! We can retrieve this surface.  Not only that, I want this code
            ! to only consider the 'bottom' limit case, above that I don't
            ! want it to bother, so I'm going to escape from the surface loop here
            exit surfLoop
          end if                        ! Is this surface ok?
        end if                          ! Was planning to retrieve here
      end do surfLoop                   ! Surface loop
    end do                              ! Instance loop

    call Deallocate_test ( mifPointings, 'mifPointings', ModuleName )
    call Deallocate_test ( MAFsForInstance, 'MAFsForInstance', ModuleName )
    call Deallocate_test ( measQtyInds, 'measQtyInds', ModuleName )
  end subroutine RestrictRange

  ! -------------------------------------------------- SetupSubset ---
  subroutine SetupSubset ( key, vectors )

    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use Expr_m, only: EXPR, GETINDEXFLAGSFROMLIST
    use MLSCommon, only: R8
    use Declaration_table, only: NUM_VALUE
    use Intrinsic, only: PHYQ_LENGTH, PHYQ_PRESSURE, PHYQ_INVALID, PHYQ_DIMENSIONLESS
    use Init_Tables_Module, only: FIELD_FIRST, FIELD_LAST
    use Init_Tables_Module, only: F_QUANTITY, F_PTANQUANTITY, F_CHANNELS, F_HEIGHT, &
      & F_MASK, F_OPTICALDEPTH, F_OPTICALDEPTHCUTOFF, F_MAXVALUE, F_MINVALUE, &
      & F_IGNORE, F_RESET, F_ADDITIONAL
    use Init_Tables_Module, only: L_RADIANCE, L_OPTICALDEPTH, L_NONE, L_ZETA, L_PRESSURE
    use Tree_Types, only: N_COLON_LESS, N_LESS_COLON, &
      & N_LESS_COLON_LESS, N_NAMED
    use VectorsModule, only: GETVECTORQTYBYTEMPLATEINDEX, SETMASK, VECTORVALUE_T, &
      & VECTOR_T, CLEARMASK, CREATEMASK, M_LINALG, DUMPMASK
    use Tree, only: NSONS, SUBTREE, DECORATION, NODE_ID
    use MoreTree, only: GET_FIELD_ID, GET_BOOLEAN
    use String_Table, only: DISPLAY_STRING
    use Toggles, only: SWITCHES
    use Output_M, only: OUTPUT
    integer, intent(in) :: KEY        ! Tree node
    type (Vector_T), dimension(:) :: VECTORS

    ! Local variables
    integer :: CHANNEL                ! Loop index
    integer :: CHANNELSNODE           ! Tree node for channels values
    integer :: COORDINATE             ! Vertical coordinate type
    integer :: FIELD                  ! Field type from tree
    integer :: GSON                   ! Tree node
    integer :: HEIGHT                 ! Loop counter
    integer :: HEIGHTNODE             ! Tree node for height values
    integer :: HEIGHTUNIT             ! Unit for heights command
    integer :: I, J                   ! Subscripts, loop inductors
    integer :: IND                    ! An array index
    integer :: INSTANCE               ! Loop counter
    integer :: INSTANCEOR1            ! For coherent quantities
    integer :: MaskBit                ! Bits corresponding to Mask
    integer :: MAINVECTORINDEX        ! Vector index of quantity to subset
    integer :: NROWS                  ! Loop limit dumping mask
    integer :: ODCUTOFFHEIGHT         ! `First' index optically thick
    integer :: QUANTITYINDEX          ! Index
    integer :: RANGEID                ! nodeID of a range
    integer :: RANGE_LOW, RANGE_HI    ! Bounds of a range
    integer :: ROW                    ! Row index dumping mask
    integer :: S1(1), S2(1)           ! Results of minloc intrinsic
    integer :: SCANDIRECTION          ! +/-1 for up or down
    integer :: SON                    ! Tree node
    integer :: STATUS                 ! Flag
    integer :: TESTUNIT               ! Either vector%globalUnit or qty tmplt unit
    integer :: TYPE                   ! Type of value returned by expr
    integer :: UNITS(2)               ! Units returned by expr
    integer :: VECTORINDEX            ! Index
    integer :: MINUNIT                ! Units for minValue
    integer :: MAXUNIT                ! Units for maxValue

    real(r8) :: HeightMin, HeightMax
    real(r8), dimension(:), pointer :: THESEHEIGHTS ! Subset of heights
    real(r8) :: VALUE(2)              ! Value returned by expr
    real(r8) :: OPTICALDEPTHCUTOFF    ! Maximum value of optical depth to allow
    real(r8) :: MAXVALUE, MINVALUE    ! Cutoff ranges
    type (VectorValue_T), pointer :: QTY ! The quantity to mask
    type (VectorValue_T), pointer :: PTAN ! The ptan quantity if needed
    type (VectorValue_T), pointer :: OPTICALDEPTH ! The opticalDepth quantity if needed
    logical :: Got(field_first:field_last)   ! "Got this field already"
    logical, dimension(:), pointer :: CHANNELS ! Are we dealing with these channels
    logical :: IGNORE                 ! Flag
    logical :: RESET                  ! Flag
    logical :: ADDITIONAL             ! Flag
    logical :: DOTHISCHANNEL          ! Flag
    logical :: DOTHISHEIGHT           ! Flag
    character(len=1), dimension(:), pointer :: ORIGINALMASK

    ! Executable code
    nullify ( channels, qty, ptan, opticalDepth )
    got = .false.
    ignore = .false.
    reset = .false.
    additional = .false.
    maskBit = m_linalg
    minUnit = 0
    maxUnit = 0
    do j = 2, nsons(key) ! fields of the "subset" specification
      son = subtree(j, key)
      field = get_field_id(son)   ! tree_checker prevents duplicates
      if (nsons(son) > 1 ) gson = subtree(2,son) ! Gson is value
      select case ( field )
      case ( f_quantity )
        mainVectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        qty => GetVectorQtyByTemplateIndex(vectors(mainVectorIndex), quantityIndex)
      case ( f_ptanquantity )
        vectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        ptan => GetVectorQtyByTemplateIndex(vectors(vectorIndeX), quantityIndex)
      case ( f_channels )
        channelsNode = son
      case ( f_height )
        heightNode = son
      case ( f_mask )
        if ( .not. got(f_mask) ) maskBit = 0 ! clear default first time
        maskBit = GetMaskBit ( son )
      case ( f_opticalDepth )
        vectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        opticalDepth => GetVectorQtyByTemplateIndex(vectors(vectorIndeX), quantityIndex)
      case ( f_opticalDepthCutoff )
        call expr ( subtree (2, son), units, value, type )
        if ( units(1) /= phyq_dimensionless ) &
          & call announceError ( key, WrongUnits, &
          & f_opticalDepthCutoff )
        opticalDepthCutoff = value(1)
      case ( f_maxValue )
        call expr ( subtree (2, son), units, value, type )
        maxUnit = units(1)
        maxValue = value(1)
      case ( f_minValue )
        call expr ( subtree (2, son), units, value, type )
        minUnit = units(1)
        minValue = value(1)
      case ( f_ignore )
        ignore = Get_Boolean ( son )
      case ( f_reset )
        reset = Get_Boolean ( son )
      case ( f_additional )
        additional = Get_Boolean ( son )
      case default
        ! Shouldn't get here if the type checker worked
      end select
      got(field) = .true.
    end do ! j = 2, nsons(key)

    ! Do some error checking for the optical depth issues
    if ( any(got((/ f_opticalDepth, f_opticalDepthCutoff /))) ) then
      if ( .not. all(got((/ f_opticalDepth, f_opticalDepthCutoff /))) ) &
        & call AnnounceError ( key, &
        & 'Must supply both opticalDepth and opicalDepthCutoff' )
      if ( qty%template%quantityType /= l_radiance .or. &
        &  opticalDepth%template%quantityType /= l_opticalDepth ) &
        & call AnnounceError ( key, 'Supplied quantity is not optical depth' )
      if ( qty%template%signal /= opticalDepth%template%signal .or. &
        &  qty%template%sideband /= opticalDepth%template%sideband ) &
        & call AnnounceError ( key, 'Optical depth does not match subsetted quantity' )
    endif

    if ( vectors(mainVectorIndex)%globalUnit /= phyq_invalid ) then
      testUnit = vectors(mainVectorIndex)%globalUnit
    else
      testUnit = qty%template%unit
    end if
    if ( got ( f_minValue ) .and. ( minUnit /= testUnit ) ) &
      & call AnnounceError ( key, WrongUnits, f_minValue )
    if ( got ( f_maxValue ) .and. ( maxUnit /= testUnit ) ) &
      & call AnnounceError ( key, WrongUnits, f_maxValue )

    ! Process the channels field.
    if ( qty%template%frequencyCoordinate /= l_none ) then
      call Allocate_test ( channels, qty%template%noChans, &
        & 'channels', ModuleName )
      if ( got(f_channels) ) then     ! This subset is only for some channels
        call GetIndexFlagsFromList ( channelsNode, channels, status, &
          & lower=lbound(channels,1) )
        if ( status /= 0 ) call announceError ( key, &
          & 'There was a problem with the channels field' )
      else
        channels = .true.             ! Apply this to all channels
      end if
    end if

    ! Check that got one of ignore, height, reset
    if ( count ( got ( (/ f_ignore, f_height, f_reset /) ) ) /= 1 ) &
      & call announceError ( key, 'Subset must be one o ignore, height or reset' )

    ! Preprocess the height stuff.  
    heightUnit = phyq_dimensionless
    if ( got(f_height) ) then
      do j = 2, nsons(heightNode)
        call expr ( subtree(j,heightNode), units, value, type )
        ! Make sure the range has non-dimensionless units -- the type
        ! checker only verifies that they're consistent.  We need to
        ! check each range separately, because the units determine the
        ! scaling of the values.
        if ( all(units == phyq_dimensionless) ) call announceError ( &
          & key, WrongUnits, f_height )
        ! Check consistency of units -- all the same, or dimensionless. The
        ! type checker verifies the consistency of units of ranges, but not
        ! of array elements.
        do i = 1, 2
          if ( heightUnit == phyq_dimensionless ) then
            heightUnit = units(i)
          else if ( units(i) /= phyq_dimensionless .and. &
            &       units(i) /= heightUnit ) then
            call announceError ( key, WrongUnits, f_height )
          end if
        end do
      end do
      ! Check for correct units
      if ( heightUnit == phyq_pressure ) then
        if ( qty%template%minorFrame .and. .not. got(f_ptanQuantity) ) &
          & call announceError ( key, WrongUnits, f_height, f_ptanQuantity )
      else if ( heightUnit /= phyq_length ) then
        call announceError ( key, WrongUnits, f_height )
      end if
    end if

    ! Create the mask if it doesn't exist
    if ( .not. associated( qty%mask ) ) call CreateMask ( qty )

    ! Make a space to save the original values if doing an additional mask
    if ( additional ) then
      nullify ( originalMask )
      call Allocate_test ( originalMask, qty%template%instanceLen, &
        & 'originalMask', ModuleName )
    end if

    ! Now we loop over the instances
    do instance = 1, qty%template%noInstances

      ! Possibly save original mask
      if ( additional ) originalMask = qty%mask ( :, instance )

      instanceOr1 = instance
      if ( qty%template%coherent ) then
        theseHeights => qty%template%surfs(:,1)
        coordinate = qty%template%verticalCoordinate
        instanceOr1 = 1
      else if ( qty%template%minorFrame .and. heightUnit == phyq_pressure ) then
        if ( ptan%template%instrumentModule /= qty%template%instrumentModule ) &
          & call AnnounceError ( key, WrongUnits, f_ptanQuantity, &
          & f_quantity )
        theseHeights => ptan%values(:,instance)
        coordinate = l_zeta
      else
        theseHeights => qty%template%surfs(:,instance)
        coordinate = qty%template%verticalCoordinate
      end if

      ! Now, make sure for the channels we're considering that the
      ! default is to ignore all, unless we've not got a heights or ignore
      ! in which case we default to use all
      if ( got(f_height) .or. ignore ) then
        do channel = 1, qty%template%noChans
          doThisChannel = .true.
          if ( associated(channels) ) doThisChannel = channels(channel)
          if ( doThisChannel ) then
            do height = 1, qty%template%noSurfs
              !??? Make sure mask bit numbers begin at 1, even when
              !??? channel numbers don't.
              call SetMask ( qty%mask(:,instance), &
                & (/ channel+qty%template%noChans*(height-1) /), &
                & what=maskBit )
            end do                    ! Height loop
          end if                      ! Do this channel
        end do                        ! Channel loop
      else
        ! If not got heights and/or ignore, default must be 
        if ( reset ) then
          do channel = 1, qty%template%noChans
            doThisChannel = .true.
            if ( associated(channels) ) doThisChannel = channels(channel)
            if ( doThisChannel ) then
              do height = 1, qty%template%noSurfs
                !??? Make sure mask bit numbers begin at 1, even when
                !??? channel numbers don't.
                call ClearMask ( qty%mask(:,instance), &
                  & (/ channel+qty%template%noChans*(height-1) /), &
                  & what=maskBit )
              end do                  ! Height loop
            end if                    ! Do this channel
          end do                      ! Channel loop
        end if                        ! Reset specified
      end if                          ! Got heights or ignore

      ! Now go and `unmask' the ones we want to consider.  For coherent
      ! quantities we can simply loop over the indices of each height range.
      ! For incoherent ones we have to go through and mark each point
      ! appropriately.
      if ( got(f_height) ) then
        do j = 2, nsons(heightNode)
          ! Get values for this ragne
          son = subtree ( j, heightNode )
          rangeId = node_id ( son )
          call expr ( son, units, value, type )
          ! Now maybe do something nasty to value to get in right units.
          if ( coordinate == l_zeta .and. heightUnit == phyq_pressure ) then
            value = -log10(value)
          else
            if ( coordinate /= qty%template%verticalCoordinate ) &
              & call AnnounceError ( key, WrongUnits, f_height )
          end if

          ! Do special things for coherent quantities
          if ( qty%template%coherent ) then
            if ( coordinate == l_pressure ) then
              s1 = minloc ( abs ( -log10(theseHeights) + log10(value(1)) ) )
              s2 = minloc ( abs ( -log10(theseHeights) + log10(value(2)) ) )
            else
              s1 = minloc ( abs ( theseHeights - value(1) ) )
              s2 = minloc ( abs ( theseHeights - value(2) ) )
            end if

            ! Now consider the open range issue
            select case ( rangeId )
            case ( n_colon_less )
              s1 = min ( s1 + 1, qty%template%noSurfs )
            case ( n_less_colon )
              s2 = max ( s2 - 1, 1 )
            case ( n_less_colon_less )
              s1 = min ( s1 + 1, qty%template%noSurfs )
              s2 = max ( s2 - 1, 1 )
            end select
          end if

          scanDirection = 0
          do channel = 1, qty%template%noChans
            doThisChannel = .true.
            if ( associated(channels) ) doThisChannel = channels(channel)
            if ( doThisChannel ) then

              ! Think about optical depth.  We want the cutoff to apply once
              ! and definitively, not have channels come in and out of use as
              ! a function of tangent height.  To do this we need to work out
              ! where we 'first' go optically thick. To do this we need to
              ! work out whether we're scanning up or down.
              if ( associated ( opticalDepth ) ) then
                ! Don't bother deducing direction if we alredy know
                if ( scanDirection == 0 ) then
                  ind = qty%template%noChans + channel
                  do height = 2, qty%template%noSurfs
                    scanDirection = scanDirection + merge ( 1, -1, &
                      & opticalDepth%template%surfs(height,instanceOr1) > &
                      & opticalDepth%template%surfs(height-1,instanceOr1) )
                  end do
                  ! Now convert it to +/-1.
                  if ( scanDirection == 0 ) scanDirection = 1 ! Default upscan
                  scanDirection = scanDirection / abs(scanDirection) 
                end if
                ! Now find `first' optically thick radiance look down from top
                if ( scanDirection == 1 ) then ! Scanning up
                  ind = qty%template%noChans*(qty%template%noSurfs-1) + channel
                  do odCutoffHeight = qty%template%noSurfs, 1, - 1 
                    if ( opticalDepth%values ( ind, instance ) > &
                      & opticalDepthCutoff ) exit
                    ind = ind - qty%template%noChans
                  end do
                else ! else scanning down
                  ind = qty%template%noChans + channel
                  do odCutoffHeight = 1, qty%template%noSurfs
                    if ( opticalDepth%values ( ind, instance ) > &
                      & opticalDepthCutoff ) exit
                    ind = ind + qty%template%noChans
                  end do
                end if
              end if

              if ( qty%template%coherent ) then
                ! For coherent quantities simply do a loop over a range.
                do height = s1(1), s2(1)
                  !??? Make sure mask bit numbers begin at 1, even when
                  !??? channel numbers don't.
                  ind = channel + qty%template%noChans*(height-1)
                  doThisHeight = .true.
                  if ( got ( f_minValue ) ) doThisHeight = doThisHeight .and. &
                    & qty%values( ind, instance ) > minValue
                  if ( got ( f_maxValue ) ) doThisHeight = doThisHeight .and. &
                    & qty%values( ind, instance ) < maxValue
                  if ( doThisHeight ) then
                    if ( associated( opticalDepth ) ) then
                      if ( scanDirection * ( height - odCutoffHeight ) > 0 ) &
                        & call ClearMask ( qty%mask(:,instance), (/ind/), what=maskBit )
                    else
                      call ClearMask ( qty%mask(:,instance), (/ind/), what=maskBit )
                    end if
                  end if
                end do                ! Height loop
              else
                ! For Incoherent quantities check each height individually
                ! We might as well consider open and non open ranges, but
                ! it's really rather unnecessary, as the chance of a ptan
                ! being exactly equal to the range given is remote, and the
                ! difference probably doesn't matter anyway.
                do height = 1, qty%template%noSurfs
                  ind = channel + qty%template%noChans*(height-1)
                  doThisHeight = .true.
                  if (any(rangeID==(/ n_less_colon,n_less_colon_less /))) then
                    doThisHeight = doThisHeight .and. theseHeights(height) > value(1)
                  else
                    doThisHeight = doThisHeight .and. theseHeights(height) >= value(1)
                  end if
                  if (any(rangeID==(/ n_colon_less,n_less_colon_less /))) then
                    doThisHeight = doThisHeight .and. theseHeights(height) < value(2)
                  else
                    doThisHeight = doThisHeight .and. theseHeights(height) <= value(2)
                  end if
                  if ( associated ( opticalDepth ) ) then
                    doThisHeight = doThisHeight .and. &
                      & scanDirection * ( height - odCutoffHeight ) > 0
                  end if
                  if ( got ( f_minValue ) ) doThisHeight = doThisHeight .and. &
                    & qty%values( ind, instance ) > minValue
                  if ( got ( f_maxValue ) ) doThisHeight = doThisHeight .and. &
                    & qty%values( ind, instance ) < maxValue
                  if ( doThisHeight ) call ClearMask ( qty%mask(:,instance), &
                    & (/ ind /), what=maskBit )
                end do                ! Height loop
              end if                  ! Coherent
            end if                    ! Do this channel
          end do                      ! Channel loop
        end do                        ! Height entries in l2cf
      end if                          ! Got a height entry

      ! If this is supposed to be an 'additional' mask, merge in the
      ! original value
      if ( additional ) &
        & qty%mask(:,instance) = char ( ior ( &
        & ichar ( qty%mask(:,instance) ), ichar ( originalMask ) ) )
    end do                            ! Instance loop

    ! Tidy up
    call Deallocate_test ( channels, 'channels', ModuleName )
    if ( additional ) call Deallocate_test ( originalMask, &
      & 'originalMask', ModuleName )

    if ( index(switches,'msk') /= 0 ) then
      call output ( 'Dumping mask' )
      if ( qty%template%name /= 0 ) then
        call output ( ' for quantity ' )
        call display_string ( qty%template%name )
      end if
      call output ( '', advance='yes' )
      call DumpMask ( qty )
    end if
  end subroutine SetupSubset

  ! -------------------------------------------------- FlagCloud ---
  ! m_cloud is the default if no maskbit is given

  subroutine SetupFlagCloud ( key, vectors )
    
    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use EXPR_M, only: EXPR, GETINDEXFLAGSFROMLIST
    use MLSCommon, only: R8
    use Init_Tables_Module, only: FIELD_FIRST, FIELD_LAST
    use Init_Tables_Module, only: F_QUANTITY, F_PTANQUANTITY, &
      & F_HEIGHT, F_CLOUDHEIGHT, F_CHANNELS, F_CLOUDCHANNELS, F_CLOUDRADIANCE, &
      & F_CLOUDRADIANCECUTOFF, F_MASK
    use Init_Tables_Module, only: L_RADIANCE, L_CLOUDINDUCEDRADIANCE, L_ZETA
    use Declaration_table, only: NUM_VALUE
    use Intrinsic, only: PHYQ_PRESSURE, PHYQ_TEMPERATURE, &
      & PHYQ_DIMENSIONLESS
    use MoreTree, only: GET_FIELD_ID
    use VectorsModule, only: ClearMask, CreateMask, &
      & GetVectorQtyByTemplateIndex, SetMask, VectorValue_T, Vector_T, &
      & M_LINALG, m_cloud
    use Tree, only: NSONS, SUBTREE, DECORATION

    integer, intent(in) :: KEY        ! Tree node
    type (Vector_T), dimension(:) :: VECTORS

    ! Local variables
    integer :: CHANNEL                ! Loop index
    integer :: CHANNELSNODE           ! Tree node for channels values
    integer :: COORDINATE             ! Vertical coordinate type
    integer :: FIELD                  ! Field type from tree
    integer :: GSON                   ! Tree node
    integer :: HEIGHT                 ! Loop counter
    integer :: HEIGHT1                 ! Loop counter
    integer :: CLOUDHEIGHTNODE             ! Tree node for cloudheight values
    integer :: HEIGHTNODE             ! Tree node for height values
    integer :: CLOUDHEIGHTUNIT             ! Unit for cloudheights command
    integer :: HEIGHTUNIT             ! Unit for heights command
    integer :: IND,IND1,J, I          ! Aarray indices
    integer :: INSTANCE               ! Loop counter
    integer :: MaskBit                ! Bits corresponding to Mask
    integer :: QUANTITYINDEX          ! Index
    integer :: RANGEID                ! nodeID of a range
    integer :: SON                    ! Tree node
    integer :: STATUS                 ! Flag
    integer :: TYPE                   ! Type of value returned by expr
    integer :: UNITS(2)               ! Units returned by expr
    integer :: VECTORINDEX            ! Index
    integer :: USETHISCHANNEL         ! cloud radiance channel

    real(r8), dimension(:), pointer :: THESEHEIGHTS ! Subset of heights
    real(r8) :: VALUE(2)              ! Value returned by expr
    real(r8) :: HeightVALUE(2)              ! Value returned by height expr
    real(r8) :: CHeightVALUE(2)              ! Value returned by cloudHeight expr
    real(r8) :: cloudRadianceCutOff              ! threshold for flagging cloud
    type (VectorValue_T), pointer :: QTY ! The quantity to mask
    type (VectorValue_T), pointer :: PTAN ! The ptan quantity if needed
    type (VectorValue_T), pointer :: cloudRadiance ! cloud radiances

    logical :: Got(field_first:field_last)   ! "Got this field already"
    logical, dimension(:), pointer :: CHANNELS ! Are we dealing with these channels
    logical :: DOTHISCLOUDHEIGHT           ! Flag
    logical :: DOTHISHEIGHT           ! Flag
    logical :: DOTHISCHANNEL          ! Flag
    logical :: ISCLOUD                ! Flag

    ! Executable code
    nullify ( channels, qty, ptan, cloudRadiance )
    got = .false.
    maskBit = m_cloud

    do j = 2, nsons(key) ! fields of the "Flagcloud" specification
      son = subtree(j, key)
      field = get_field_id(son)   ! tree_checker prevents duplicates
      if (nsons(son) > 1 ) gson = subtree(2,son) ! Gson is value
      select case ( field )
      case ( f_quantity )
        vectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        qty => GetVectorQtyByTemplateIndex(vectors(vectorIndeX), quantityIndex)
      case ( f_ptanquantity )
        vectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        ptan => GetVectorQtyByTemplateIndex(vectors(vectorIndeX), quantityIndex)
      case ( f_height )
        heightNode = son
      case ( f_cloudHeight )
        cloudHeightNode = son
      case ( f_channels )
        channelsNode = son
      case ( f_cloudChannels )
        if ( nsons(son) /= 2 ) call AnnounceError ( key, &
          & 'Too many cloudChannels' )
        call expr ( subtree(2,son), units, value )
        useThisChannel = nint ( value(1) )
      case ( f_cloudRadiance )
        vectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        cloudRadiance => GetVectorQtyByTemplateIndex(vectors(vectorIndeX), quantityIndex)
      case ( f_cloudRadianceCutoff )
        call expr ( subtree (2, son), units, Value )
        if ( units(1) /= PHYQ_temperature ) &
          & call announceError ( key, WrongUnits, f_cloudRadianceCutoff )
        cloudRadianceCutoff = value(1)
      case ( f_mask )
        maskBit = GetMaskBit ( son )
      case default
        ! Shouldn't get here if the type checker worked
      end select
      got(field) = .true.
    end do ! j = 2, nsons(key)

    ! Quantity must be radiance
    if ( (.not. any(qty%template%quantityType==(/ l_radiance,l_cloudinducedradiance /))) &
      & .or. (cloudRadiance%template%quantityType /= l_cloudInducedRadiance &
      & .and. cloudRadiance%template%quantityType /= l_radiance)) &
      & call AnnounceError ( key, 'You supplied the wrong quantities to flagCloud' )

    if( got(f_channels)) then
      call Allocate_test ( Channels, qty%template%noChans, &
        & 'channels', ModuleName )
      call GetIndexFlagsFromList ( channelsNode, channels, status, &
        & lower=lbound(channels,1) )
    end if

    ! Preprocess the height stuff.
    heightUnit = phyq_dimensionless
    if ( got(f_height) ) then
      if (nsons(heightNode) .gt. 2) call AnnounceError ( key, &
        & 'Only one height range allowed', f_height )

      call expr ( subtree(2,heightNode), units, value, type )
      ! Make sure the range has non-dimensionless units -- the type
      ! checker only verifies that they're consistent.  We need to
      ! check each range separately, because the units determine the
      ! scaling of the values.
      if ( all(units == phyq_dimensionless) ) call announceError ( &
        & key, WrongUnits, f_height )
      ! Check consistency of units -- all the same, or dimensionless. The
      ! type checker verifies the consistency of units of ranges, but not
      ! of array elements.
      do i = 1, 2
        if ( heightUnit == phyq_dimensionless ) then
          heightUnit = units(i)
        else if ( units(i) /= phyq_dimensionless .and. &
          &       units(i) /= heightUnit ) then
          call announceError ( key, WrongUnits, f_height )
        end if
      end do

      ! convert pressure unit
      Heightvalue = value
      if ( heightUnit == phyq_pressure ) Heightvalue = -log10(value)

    end if
    ! Preprocess the height stuff.
    cloudheightUnit = phyq_dimensionless
    if ( got(f_cloudheight) ) then
      if (nsons(cloudheightNode) .gt. 2) call AnnounceError ( key, &
        & 'Can only supply one height range' ) 

      call expr ( subtree(2,cloudheightNode), units, value, type )
      if ( all(units == phyq_dimensionless) ) call announceError ( &
        & key, WrongUnits, f_cloudheight )
      do i = 1, 2
        if ( cloudheightUnit == phyq_dimensionless ) then
          cloudheightUnit = units(i)
        else if ( units(i) /= phyq_dimensionless .and. &
          &       units(i) /= cloudheightUnit ) then
          call announceError ( key, 'Problem with units', f_cloudheight )
        end if
      end do

      ! convert pressure unit
      CHeightvalue = value
      if ( cloudheightUnit == phyq_pressure ) CHeightvalue = -log10(value)          
    end if

    ! ----- finish checking ------

    ! Create the mask if it doesn't exist
    if ( .not. associated( qty%mask ) ) call CreateMask ( qty )

    ! Now loop over the instances
    do instance = 1, qty%template%noInstances
      theseHeights => ptan%values(:,instance)
      coordinate = l_zeta

      isCloud = .false.

      ! use the given height range to find cloud flag
      if( got(f_cloudheight) ) then 
        do height1 = 1,cloudRadiance%template%noSurfs
          doThisCloudHeight = .true.
          doThisCloudHeight = doThisCloudHeight .and. &
            & theseHeights(height1) > Cheightvalue(1)
          doThisCloudHeight = doThisCloudHeight .and. &
            & theseHeights(height1) < Cheightvalue(2)
          ind1 = useThisChannel + cloudRadiance%template%noChans*(height1-1)
          if (cloudRadianceCutoff > 0._r8) then
             if ( doThisCloudHeight .and. &
               & cloudRadiance%values ( ind1, instance ) > cloudRadianceCutoff) &
               & isCloud = .true.
          else
             if ( doThisCloudHeight .and. &
               & cloudRadiance%values ( ind1, instance ) < cloudRadianceCutoff) &
               & isCloud = .true.
          endif
        enddo
      endif
      do height = 1, qty%template%noSurfs
        doThisHeight = .true.
        doThisHeight = doThisHeight .and. &
          & theseHeights(height) > heightvalue(1)
        doThisHeight = doThisHeight .and. &
          & theseHeights(height) < heightvalue(2)

        ! determine cloud flag
        ! if cloud flag is not from a height range then do comparison at each minor frame
        if ( .not. got(f_cloudheight) ) then 
          ! use the same mif radiance to find cloud flag
          isCloud = .false.
          ind1 = useThisChannel + cloudRadiance%template%noChans*(height-1)
          if (cloudRadianceCutoff > 0._r8) then
            if ( cloudRadiance%values ( ind1, instance ) > cloudRadianceCutoff) &
            & isCloud = .true.
          else
            if ( cloudRadiance%values ( ind1, instance ) < cloudRadianceCutoff) &
            & isCloud = .true.
          endif
        endif

        ! assign cloud flag
        do channel = 1, qty%template%noChans
          doThisChannel = .true.
          if ( associated(channels) ) doThisChannel = channels(channel)
          if ( doThisChannel ) then
            ind = channel + qty%template%noChans*(height-1)
            if ( doThisHeight .and. isCloud )  &
              &     call SetMask ( qty%mask(:,instance), (/ ind /), &
              &     what=maskBit )
          end if                        ! do this channel
        end do                          ! Channel loop
      end do                            ! height loop
    end do                              ! Instance loop

    ! Tidy up
    if ( associated(channels) ) &
      & call Deallocate_test ( channels, 'channels', ModuleName )

  end subroutine SetupFlagCloud

  ! ---------------------------------------------- UpdateMask -----
  subroutine UpdateMask ( key, vectors )

    use Init_Tables_Module, only: FIELD_FIRST, FIELD_LAST
    use Init_Tables_Module, only: F_QUANTITY, F_SOURCEQUANTITY, &
      & F_OPERATION, F_SOURCEMASK, F_MASK
    use Init_Tables_Module, only: L_INVERT, L_COPY, L_ANDMASKS, L_ORMASKS
    use VectorsModule, only: GETVECTORQTYBYTEMPLATEINDEX, SETMASK, VECTORVALUE_T, &
      & VECTOR_T, CLEARMASK, CREATEMASK, M_LINALG, DUMPMASK
    use Tree, only: NSONS, SUBTREE, DECORATION, NODE_ID
    use MoreTree, only: GET_FIELD_ID
    use Toggles, only: SWITCHES
    use Output_M, only: OUTPUT

    integer, intent(in) :: KEY          ! Tree node
    type (Vector_T), dimension(:) :: VECTORS

    ! Local variables
    logical :: GOT(field_first:field_last)   ! "Got this field already"
    integer :: SON                      ! Tree node
    integer :: GSON                     ! Another tree node
    integer :: MASKBIT                  ! Which bits
    integer :: SOURCEMASKBIT            ! Which bits from source
    integer :: J                        ! Loop counter
    integer :: FIELD                    ! A field id
    integer :: VECTORINDEX              ! In database
    integer :: QUANTITYINDEX            ! In database not vector
    integer :: OPERATION                ! What operation to perform
    type (VectorValue_T), pointer :: QUANTITY
    type (VectorValue_T), pointer :: SOURCEQUANTITY

    ! Executable code
    ! First work out what we've been asked to do

    ! Executable code
    nullify ( sourceQuantity )
    got = .false.
    do j = 2, nsons(key) ! fields of the "subset" specification
      son = subtree(j, key)
      field = get_field_id(son)   ! tree_checker prevents duplicates
      if (nsons(son) > 1 ) gson = subtree(2,son) ! Gson is value
      select case ( field )
      case ( f_quantity )
        vectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        quantity => GetVectorQtyByTemplateIndex(vectors(vectorIndex), quantityIndex)
      case ( f_sourceQuantity )
        vectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        sourceQuantity => GetVectorQtyByTemplateIndex(vectors(vectorIndeX), quantityIndex)
      case ( f_mask )
        if ( .not. got(f_mask) ) maskBit = 0 ! clear default first time
        maskBit = GetMaskBit ( son, single=.true. )
      case ( f_sourceMask )
        if ( .not. got(f_sourceMask) ) sourceMaskBit = 0 ! clear default first time
        sourceMaskBit = GetMaskBit ( son, single=.true. )
      case ( f_operation )
        operation = decoration(gson)
      case default
        ! Shouldn't get here if the type checker worked
      end select
      got(field) = .true.
    end do ! j = 2, nsons(key)

    ! Special checking for invert case
    if ( operation == l_invert ) then
      if ( associated(sourceQuantity) ) call AnnounceError ( key, &
        & 'sourceQuantity not appropriate for invert updateMask' )
      if ( got(f_sourceMask) ) call AnnounceError ( key, &
        & 'sourceMask not appropriate for invert updateMask' )
    end if

    ! If we've got a source quantity then check it out.
    if ( associated ( sourceQuantity ) ) then
      if ( sourceQuantity%template%name /= quantity%template%name ) &
        & call AnnounceError ( key, 'sourceQuantity and quantity not compatible' )
      if ( .not. got(f_sourceMask) ) sourceMaskBit = maskBit
    else
      if ( operation /= l_invert ) then
        if ( .not. got(f_sourceMask) ) call AnnounceError ( key, &
          & 'sourceMask required for this operation' )
        if ( sourceMaskBit == maskBit ) call AnnounceError ( key, &
          & 'mask and sourceMask are the same information' )
      end if
      sourceQuantity => quantity
    end if

    ! If the result quantity has no mask then create one.
    if ( .not. associated ( quantity%mask ) ) call CreateMask ( quantity )
    ! Create one in the source too if it doesn't have one, potentially time wasting
    ! later on, but unlikely
    if ( .not. associated ( sourceQuantity%mask ) ) call CreateMask ( sourceQuantity )

    ! Now work out what to do.
    select case ( operation )
    case ( l_invert )
      quantity%mask = char ( &
        & ior ( iand ( ichar(quantity%mask), not(maskBit) ), & ! Other bits
        &       iand ( not ( iand ( ichar(quantity%mask), maskBit ) ), maskBit ) ) )
    case ( l_copy )
      quantity%mask = char ( &
        & ior ( iand ( ichar(quantity%mask), not(maskBit) ), & ! Other bits
        &       merge ( maskBit, 0, &
        &              iand ( ichar(sourceQuantity%mask), sourceMaskBit ) /= 0 ) ) )
    case ( l_andMasks )
      ! Note, because we 'invert' the logic in storing the mask
      ! this operation may look wrong at first glance
      quantity%mask = char ( &
        & ior ( iand ( ichar(quantity%mask), not(maskBit) ), & ! Other bits
        &       merge ( 0, maskBit, &
        &               iand ( ichar(sourceQuantity%mask), sourceMaskBit ) == 0 .and. &
        &               iand ( ichar(quantity%mask), maskBit ) == 0  ) ) )
    case ( l_orMasks )
      ! Note, because we 'invert' the logic in storing the mask
      ! this operation may look wrong at first glance
      quantity%mask = char ( &
        & ior ( iand ( ichar(quantity%mask), not(maskBit) ), & ! Other bits
        &       merge ( 0, maskBit, &
        &               iand ( ichar(sourceQuantity%mask), sourceMaskBit ) == 0 .or. &
        &               iand ( ichar(quantity%mask), maskBit ) == 0  ) ) )
    end select

  end subroutine UpdateMask
 
  ! ===================================== Private procedures ============

  ! --------------------------------------- AnnounceError -----------
  
  subroutine AnnounceError ( NODE, STRING, FIELDINDEX, ANOTHERFIELDINDEX )
    
    use Intrinsic, only: FIELD_INDICES, SPEC_INDICES
    use Lexer_Core, only: PRINT_SOURCE
    use String_Table, only: DISPLAY_STRING
    use Tree, only: SOURCE_REF
    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR
    use Output_M, only: OUTPUT
    
    integer, intent(in) :: NODE
    character(len=*), intent(in) :: STRING
    integer, intent(in), optional :: FIELDINDEX, ANOTHERFIELDINDEX ! f_...
    
    call output ( '***** At ' )
    call print_source ( source_ref(node) )
    call output ( ', SubsetModule complained: ', advance='yes' )
    call output ( trim(string), advance='yes' )
    if ( present ( fieldIndex ) ) then
      call output ( 'Offending field(s): ' )
      call output ( field_indices(fieldIndex) )
      if ( present ( anotherFieldIndex ) ) then
        call output ( ', ' )
        call output ( field_indices(anotherFieldIndex) )
      end if
      call output ( '', advance='yes' )
    end if
    call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Problem in SubsetModule, see above' )
  end subroutine AnnounceError

  ! -------------------------------------- GetMaskBit -------------

  function GetMaskBit ( node, single ) result ( maskBit )
    use Tree, only: NSONS, DECORATION, SUBTREE
    use Init_Tables_Module, only: L_FILL, L_FULL_DERIVATIVES, &
      & L_TIKHONOV, L_LINALG, L_SPARE, L_CLOUD
    use VectorsModule, only: M_FILL, M_FULLDERIVATIVES, M_LINALG, &
      & M_TIKHONOV, M_SPARE, M_CLOUD

    ! This routine parses the mask field in an l2cf command
    ! Dummy arguments
    integer, intent(in) :: NODE         ! Tree node
    logical, optional, intent(in) :: SINGLE ! Raise error if more than one
    integer :: MASKBIT                  ! Result
    integer :: MASK                     ! Tree decoration
    integer :: I
    if ( present ( single ) ) then
      if ( single .and. nsons(node) > 2 ) &
        & call AnnounceError ( node, 'Only one mask expected' )
    end if
    maskBit = 0
    do i = 2, nsons(node)
      mask = decoration(subtree(i,node))
      select case ( mask )
      case ( l_Cloud )
        maskBit = ior(maskBit, m_Cloud)
      case ( l_fill )
        maskBit = ior(maskBit, m_fill)
      case ( l_full_derivatives )
        maskBit = ior(maskBit, m_fullDerivatives)
      case ( l_linAlg )
        maskBit = ior(maskBit, m_linAlg)
      case ( l_Spare )
        maskBit = ior(maskBit, m_Spare)
      case ( l_Tikhonov )
        maskBit = ior(maskBit, m_Tikhonov)
      end select
    end do
  end function GetMaskBit

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module SubsetModule
 
! $Log$
! Revision 2.12  2004/05/28 00:57:50  vsnyder
! Move GetIndexFlagsFromList from MoreTree to Expr_m
!
! Revision 2.11  2003/08/06 17:24:16  livesey
! Essentially cosmetic change
!
! Revision 2.10  2003/05/14 22:02:08  dwu
! fix a bug in flagcloud
!
! Revision 2.9  2003/05/12 18:55:01  dwu
! allow cloud flag for cloudinducedradiance
!
! Revision 2.8  2003/04/14 22:21:26  dwu
! make cloudradiancecutoff sign dependent: it is 'greater than' if positive; 'less than' if negative
!
! Revision 2.7  2003/04/10 18:30:57  dwu
! make m_cloud as default in FlagCloud
!
! Revision 2.6  2003/04/08 23:12:18  dwu
! add m_cloud in setupFlagCloud
!
! Revision 2.5  2003/04/04 22:01:46  livesey
! Added UpdateMask
!
! Revision 2.4  2003/03/19 01:58:36  livesey
! Added nullify on signals
!
! Revision 2.3  2003/03/07 03:17:12  livesey
! Added RestrictRange
!
! Revision 2.2  2003/03/06 00:47:04  livesey
! Added logging information
!
