! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE ConstructQuantityTemplates ! Construct templates from user supplied info
!=============================================================================

  ! This module has various functionality for constructing quantity templates.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use EXPR_M, only: EXPR
  use FGrid, only: fGrid_T
  use HGrid, only: hGrid_T
  use INIT_TABLES_MODULE, only: F_GEODANGLE, F_FGRID, F_HGRID, &
    & F_LOGBASIS, F_MINVALUE, F_MODULE, F_MOLECULE, F_NOMIFS, F_RADIOMETER, &
    & F_SIGNAL, F_SGRID, F_TYPE, F_UNIT, F_VGRID, F_IRREGULAR, &
    & F_SOLARTIME, F_SOLARZENITH
  use INIT_TABLES_MODULE, only: &
    FIRST_LIT, LAST_LIT, L_BASELINE, L_BOUNDARYPRESSURE, &
    L_CHANNEL, L_CHISQCHAN, L_CHISQMMAF, L_CHISQMMIF, L_CHUNK, L_CLOUDICE, &
    L_CLOUDEXTINCTION, L_CLOUDWATER, &
    L_TOTALEXTINCTION, L_MASSMEANDIAMETERICE, &
    L_CLOUDINDUCEDRADIANCE, L_CLOUDRADSENSITIVITY, &
    L_COLUMNABUNDANCE, L_DNWT_AJN, L_DNWT_AXMAX, &
    L_DNWT_CAIT, L_DNWT_CHISQMINNORM, L_DNWT_CHISQNORM, L_DNWT_DIAG, &
    L_DNWT_DXDX, L_DNWT_DXDXL, L_DNWT_DXN, L_DNWT_DXNL, L_DNWT_FLAG, &
    L_DNWT_FNMIN, L_DNWT_FNORM, L_DNWT_GDX, L_DNWT_GFAC, L_DNWT_GRADN, &
    L_DNWT_SQ, L_DNWT_SQT, &
    L_EARTHREFL, L_EARTHRADIUS, L_EFFECTIVEOPTICALDEPTH, &
    L_ELEVOFFSET, L_EXTINCTION, L_FREQUENCY, L_GEODALTITUDE, L_GPH, &
    L_HEIGHTOFFSET, L_ITERATION, L_JACOBIAN_COLS, L_JACOBIAN_ROWS, &
    L_LOSTRANSFUNC, L_LOSVEL, L_MAGNETICFIELD, &
    L_MAF, L_MASSMEANDIAMETERICE, L_MASSMEANDIAMETERWATER, L_MIF, &
    L_NOISEBANDWIDTH, L_NONE, L_NUMJ, L_ORBITINCLINATION, L_OPTICALDEPTH, &
    L_PHITAN, L_PTAN, L_RADIANCE, L_RHI, &
    L_REFGPH, L_SCANRESIDUAL, L_SCECI, L_SCGEOCALT, L_SCVEL, &
    L_SCVELECI, L_SCVELECR, L_SIDEBANDRATIO, L_SIZEDISTRIBUTION, &
    L_SPACERADIANCE, L_SURFACETYPE, L_SYSTEMTEMPERATURE, &
    L_TEMPERATURE, L_TNGTECI, L_TNGTGEOCALT, L_TNGTGEODALT, &
    L_TRUE,&
    L_VMR, L_XYZ, PHYQ_ANGLE, PHYQ_DIMENSIONLESS, PHYQ_EXTINCTION, &
    PHYQ_FREQUENCY, PHYQ_DOBSONUNITS, PHYQ_IceDensity, PHYQ_LENGTH, PHYQ_PCTRHI, &
    PHYQ_PRESSURE, PHYQ_TEMPERATURE, PHYQ_VELOCITY, PHYQ_VMR, &
    PHYQ_ZETA, PHYQ_GAUSS, PHYQ_TIME
  use L1BData, only: L1BData_T, READL1BDATA, DEALLOCATEL1BDATA, &
    & FindL1BData, AssembleL1BQtyName, PRECISIONSUFFIX
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: L1BInfo_T, MLSChunk_T, NameLen, R8
  use MLSFiles, only: mls_hdf_version
  use MLSL2Options, only: LEVEL1_HDFVERSION
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_L1BRead
  use MLSSignals_m, only:  IsModuleSpacecraft, &
    & GetModuleFromRadiometer, GetModuleFromSignal, GetModuleName, &
    & GetRadiometerFromSignal, GetSignal,&
    & GetSignalName, Signal_T, SIGNALS, MODULES
  use OUTPUT_M, only: OUTPUT
  use Parse_Signal_m, only: PARSE_SIGNAL
  use QuantityTemplates, only: QuantityTemplate_T,SetupNewQuantityTemplate, &
    & QuantityTemplateCounter, NullifyQuantityTemplate
  use STRING_TABLE, only: GET_STRING, DISPLAY_STRING
  use TOGGLES, only: GEN, TOGGLE, SWITCHES, LEVELS
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATION, NODE_ID, NSONS, SOURCE_REF, SUB_ROSA, SUBTREE
  use TREE_TYPES, only: N_SET_ONE
  use VGridsDatabase, only: VGrid_T

  implicit none
  private
  public :: ConstructMinorFrameQuantity, CreateQtyTemplateFromMLSCFInfo
  public :: ForgeMinorFrames

! -----     Private declarations     -----------------------------------

  integer :: ERROR
  logical, parameter :: DEEBUG = .FALSE.                 ! Normally FALSE
  logical, parameter :: ALWAYSFIRSTSIGNAL = .FALSE.      ! Normally FALSE
  logical, parameter :: BE_WHINY_ABOUT_IT = .FALSE.      ! Normally FALSE
  integer, parameter :: HOWSEVEREISNOGOODSIGNAL = 0      ! Normally 0

! Error codes for "announce_error"
  integer, parameter :: No_Error_Code = 0
  integer, parameter :: BadUnitMessage =         No_Error_Code+1
  integer, parameter :: InappropriateQuantity =  BadUnitMessage+1
  integer, parameter :: NeedGrid =               InappropriateQuantity+1
  integer, parameter :: NoQuantityType =         NeedGrid+1
  integer, parameter :: UnnecessaryGrid =        NoQuantityType+1
  integer, parameter :: noModule=                UnnecessaryGrid+1

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================
  ! -----------------------------  CreateQtyTemplateFromMLSCFInfo  -----
  type (QuantityTemplate_T) function CreateQtyTemplateFromMLSCFInfo ( &
    & Name, Root, FGrids, VGrids, HGrids, L1bInfo, Chunk, MifGeolocation, &
    & returnStatus ) &
    result ( QTY )

  ! This routine constructs a vector quantity template based on instructions
  ! passed in an mlscf line.

    ! Dummy arguments
    integer, intent(in) :: NAME              ! Sub-rosa index of name
    integer, intent(in) :: ROOT              ! Root of QuantityTemplate subtree
    type (FGrid_T), dimension(:), pointer :: FGrids
    type (VGrid_T), dimension(:), pointer :: VGrids
    type (HGrid_T), dimension(:), pointer :: HGrids
    type (l1bInfo_T), intent(in) :: L1bInfo
    type (MLSChunk_T), intent(in) :: Chunk
    type (QuantityTemplate_T), dimension(:), intent(in), optional :: &
      & MifGeolocation
    integer, intent(out) :: returnStatus      ! 0 unless trouble

    ! Local variables
    integer, dimension(2) :: EXPR_UNITS
    real(r8), dimension(2) :: EXPR_VALUE
    integer :: Family
    integer :: FGridIndex
    logical, save :: First = .true.
    integer :: FrequencyCoordinate
    integer :: HGridIndex
    integer :: I                        ! Loop counter
    integer :: IERR                     ! So display_string doesn't crash and burn
    integer :: InstrumentModule         ! Database index
    integer :: Key            ! Field name, F_... from Init_Tables_Module
    character(len=127) :: SIGNALSTRING
    logical :: LOGBASIS                 ! To place in quantity
    integer :: Molecule
    character(len=127) :: moleculeString
    integer, save :: Natural_Units(first_lit:last_lit)
    integer :: NoInstances
    integer :: NoSurfs
    real(r8) :: MinValue                ! Minimumvalue
    logical :: MinorFrame               ! Is a minor frame quantity
    logical :: MajorFrame               ! Is a major frame quantity
    integer :: NoChans
    integer :: QuantityType
    integer :: Radiometer               ! Database index
    logical :: Regular                  ! Flag
    real(r8) :: ScaleFactor
    integer :: Sideband
    integer :: Signal                   ! Database index
    integer :: s_index
    integer, dimension(:), pointer :: SignalInds ! From parse signal
    type (signal_T) :: SignalInfo       ! Details of the appropriate signal
    integer :: Son                      ! A Son of Root -- an n_assign node
    integer :: Type_Field               ! Index in subtree of "type"
    integer :: Value                    ! Node index of value of field of spec
    integer :: VGridIndex
    integer :: SGridIndex

    ! Executable code

    if ( toggle(gen) .and. levels(gen) > 0 ) &
      & call trace_begin ( "CreateQtyTemplateFromMLSCFInfo", root )

! ??? Do we need a GOT_FIELD check like in VGrid, e.g. ???

    call nullifyQuantityTemplate ( qty ) ! for Sun's rubbish compiler
    nullify ( signalInds )
    error = 0
    family = 0
    hGridIndex = 0
    instrumentModule = 0
    logBasis = .false.
    regular = .true.
    minValue = -huge ( 0.0_r8 )
    molecule = 0

    if ( first ) then
      ! Fill NATURAL_UNITS on first call only (it's a SAVE variable)
      first = .false.
      natural_units = 0
      natural_units(l_baseline) =                PHYQ_Temperature
      natural_units(l_boundaryPressure) =        PHYQ_Pressure
      natural_units(l_chisqchan) =               PHYQ_Dimensionless
      natural_units(l_chisqmmaf) =               PHYQ_Dimensionless
      natural_units(l_chisqmmif) =               PHYQ_Dimensionless
      natural_units(l_columnAbundance) =         PHYQ_DobsonUnits
      natural_units(l_cloudice) =                PHYQ_IceDensity
      natural_units(l_cloudextinction) =         PHYQ_Dimensionless
      natural_units(l_totalextinction) =         PHYQ_Dimensionless
      natural_units(l_massmeandiameterice) =     PHYQ_Dimensionless
      natural_units(l_earthRefl) =               PHYQ_Dimensionless
      natural_units(l_elevOffset) =              PHYQ_Angle
      natural_units(l_extinction) =              PHYQ_Extinction
      natural_units(l_gph) =                     PHYQ_Length
      natural_units(l_heightOffset ) =           PHYQ_Length
      natural_units(l_losTransFunc) =            PHYQ_Dimensionless
      natural_units(l_losVel) =                  PHYQ_Velocity
      natural_units(l_orbitInclination) =        PHYQ_Angle
      natural_units(l_noiseBandwidth) =          PHYQ_Frequency
      natural_units(l_phitan) =                  PHYQ_Angle
      natural_units(l_ptan) =                    PHYQ_Zeta
      natural_units(l_radiance) =                PHYQ_Temperature
      natural_units(l_cloudinducedradiance) =    PHYQ_Temperature
      natural_units(l_cloudradsensitivity) =     PHYQ_Temperature
      natural_units(l_effectiveopticaldepth) =   PHYQ_Dimensionless
      natural_units(l_earthradius) =             PHYQ_Length
      natural_units(l_refGPH) =                  PHYQ_Length
      natural_units(l_rhi) =                     PHYQ_PctRHI
      natural_units(l_scGeocAlt ) =              PHYQ_Length
      natural_units(l_scVel) =                   PHYQ_Velocity
      natural_units(l_scVelECI) =                PHYQ_Velocity   
      natural_units(l_scVelECR) =                PHYQ_Velocity   
      natural_units(l_scanResidual ) =           PHYQ_Length
      natural_units(l_spaceRadiance) =           PHYQ_Temperature
      natural_units(l_systemTemperature) =       PHYQ_Temperature
      natural_units(l_temperature) =             PHYQ_Temperature
      natural_units(l_tngtGeocAlt) =             PHYQ_Length
      natural_units(l_tngtGeodAlt) =             PHYQ_Length
      natural_units(l_vmr) =                     PHYQ_Vmr
      natural_units(l_magneticField) =           PHYQ_Gauss
    end if

    noChans = 1
    quantitytype = 0
    radiometer = 0
    scaleFactor = 1.0
    sideband = 0
    signal = 0
    fGridIndex = 0
    vGridIndex = 0
    sGridIndex = 0
    signalString = ' '
    moleculeString = ' '

    ! First we'll loop over the MLSCF keys.

    do i = 2, nsons(root)
      son = subtree(i,root)
      key = subtree(1,son)
      if ( node_id(son) == n_set_one ) then
        value = l_true
      else
        value = decoration(subtree(2,son))
      end if

      select case ( decoration(key) )
      case ( f_fgrid )
        fGridIndex = decoration(value)
      case ( f_hgrid )
        hGridIndex = decoration(value) ! node_id(value) == n_spec_args
      case ( f_logBasis )
        logBasis = (value == l_true)
      case ( f_irregular )
        regular = (value /= l_true)
      case ( f_minValue )
        call expr ( subtree(2,son), expr_units, expr_value )
        minValue = expr_value(1)
      case ( f_module)
        instrumentModule = decoration(decoration(subtree(2,son)))
!      case ( f_molecule );          molecule = value  ! I don't understand this
      case ( f_molecule )
        molecule = value
        call get_string( sub_rosa(subtree(2,son)), moleculeString, strip=.true. )
      case ( f_radiometer )
        radiometer = decoration(decoration(subtree(2,son)))
        instrumentModule = GetModuleFromRadiometer(radiometer)
      case ( f_signal )
        !??? For the moment it is simple, later we'll be more intelligent here
        !??? for example, letting the user choose either R1A or R1B.
        call get_string( sub_rosa(subtree(2,son)), signalString, strip=.true. )
        !??? Here we would do intelligent stuff to work out which bands
        !??? are present, for the moment choose the first
        call parse_Signal ( signalString, signalInds, &
          & tree_index=son, sideband=sideband )
        if ( .not. associated(signalInds) ) then ! A parse error occurred
          call MLSMessage ( MLSMSG_Error, ModuleName,&
            & 'Unable to parse signal string' )
        end if
        ! print *, 'Parsed ', trim(signalString), ' for ', size(signalInds), &
        !  & ' bands'
        if ( size(signalInds) == 1 .or. ALWAYSFIRSTSIGNAL &
          & .or. .not. associated(L1bInfo%L1BRADIds) ) then
          signal = signalInds(1)
        else
          ! Seek a signal with any precision values !< 0
          do s_index=1, size(signalInds)
            if ( any_good_signaldata ( signalInds(s_index), sideband, &
              & l1bInfo, Chunk) ) exit
          enddo
          if ( s_index > size(signalInds) ) then
            if (BE_WHINY_ABOUT_IT ) call announce_error (root, No_Error_Code, &
              & 'Warning: no good signal data found for ' &
              & // trim(signalString), HOWSEVEREISNOGOODSIGNAL)
            signal = signalInds(1)
          else
            signal = signalInds(s_index)
          endif
        endif
        call deallocate_test ( signalInds, 'signalInds', ModuleName )
        instrumentModule = GetModuleFromSignal(signal)
        radiometer = GetRadiometerFromSignal(signal)
      case ( f_sgrid )
        sGridIndex = decoration(value) ! node_id(value) == n_spec_args
      case ( f_type )
        quantityType = value
        type_field = son
      case ( f_unit );              scaleFactor = value
      case ( f_vgrid )
        vGridIndex = decoration(value) ! node_id(value) == n_spec_args
      end select
    end do

    ! Now we know what the user asked for, try to make sense of it.
    ! First see that the `type' has been defined

    if ( quantityType == 0 ) call announce_error ( root, noQuantityType )
 
    ! Set defaults for other parameters
    if ( fGridIndex /= 0 ) then
      frequencyCoordinate = fGrids(fGridIndex)%frequencyCoordinate
      noChans = fGrids(fGridIndex)%noChans
    else
      frequencyCoordinate = L_None
      noChans = 1
    endif

    ! Now, depending on the type, check out stuff and see if it's ok to
    ! first order.

    if ( family == 0 ) family = natural_units(quantityType)
    minorFrame = any(quantityType == (/ l_phiTan, l_Ptan, l_Radiance, &
      & l_cloudInducedRadiance, l_cloudRADSensitivity, l_effectiveOpticalDepth, &
      & l_tngtECI, l_tngtGeodAlt, l_tngtGeocAlt, l_scECI, l_scGeocAlt,&
      & l_scVel, l_scVelECI, l_scVelECR, l_losVel, l_heightOffset, &
      & l_scanResidual, l_chisqmmif, l_opticalDepth, l_orbitInclination /) )

    majorFrame = any(quantityType == (/ l_chisqchan, l_chisqmmaf /) )
 
    ! Here the code splits, for minor frame quantities, we take the information
    ! from the previously constructed MIFGeolocation information.  Otherwise,
    ! we'll probably need to use any supplied vGrid/hGrid information.

    if ( minorFrame ) then

      ! This is a minor frame type quantity.
      if ( (hGridIndex /= 0) .OR. (vGridIndex /= 0 )) then
        call announce_error ( root, unnecessaryGrid )
      end if
      if ( instrumentModule == 0 ) then 
        call announce_error ( root, noModule )
      end if

      ! Work out the channel information
      if ( signal /= 0 ) then
        signalInfo = GetSignal(signal)
        noChans = size(signalInfo%frequencies)
        frequencyCoordinate = l_channel
      end if
    
      ! For some cases we know the quantity is an xyz vector
      if ( any(quantityType == &
       & (/ l_tngtECI, l_scECI, l_scVel, l_scVelECI, l_scVelECR /)) &
       &  ) then
        noChans = 3
        frequencyCoordinate = l_xyz
      end if

      ! Construct an empty quantity
      call ConstructMinorFrameQuantity ( l1bInfo, chunk, instrumentModule, &
        & qty, noChans=noChans, mifGeolocation=mifGeolocation, regular=regular )

      ! Make absolutely certain template's dimensions are what we want
      if ( quantityType == l_chiSqMMIF ) then
        qty%noChans = 1
        qty%instanceLen = qty%NoSurfs
        qty%frequencyCoordinate = l_none
      endif
        
    elseif ( majorFrame ) then

      ! This is a major frame type quantity.
      if ( vGridIndex/=0 ) then
        call announce_error ( root, No_Error_Code, &
        &  'No vGrid should be specified for a major frame quantity' )
      endif

      ! Work out the channel information
      if ( signal == 0 ) then
        call announce_error ( root, No_Error_Code, &
        &  'A signal is required for every major frame quantity' )
      elseif ( instrumentModule == 0 ) then 
        call announce_error ( root, noModule )
      else
        signalInfo = GetSignal(signal)
        noChans = size(signalInfo%frequencies)
        frequencyCoordinate = l_channel
      end if

      ! Make absolutely certain template's dimensions are what we want
      if ( quantityType == l_chisqMMAF .or. quantityType == l_chisqMMIF ) then
        noChans = 1
        frequencyCoordinate = l_none
      end if

      ! Construct an empty quantity
      call ConstructMajorFrameQuantity ( chunk, instrumentModule, &
        & qty, noChans, mifGeolocation )
      qty%frequencyCoordinate = frequencyCoordinate
        
   ! for losTransFunc type of quantity 
   elseif (quantityType == l_losTransFunc) then
      ! replace noChans with no of grid along path which is specified from sGrid
        noChans = VGrids(sGridIndex)%noSurfs
        frequencyCoordinate = l_losTransFunc

       ! Construct an empty quantity
       call ConstructMinorFrameQuantity ( l1bInfo, chunk, instrumentModule, &
        & qty, noChans=noChans )
       call Allocate_test ( qty%frequencies, noChans, 'qty%frequencies', ModuleName )
       
       qty%frequencies = VGrids(sGridIndex)%surfs
        
   else

      ! This is not a minor frame quantity, set it up from FGrids, VGrids and HGrids

      if ( hGridIndex/=0 ) then
        noInstances=hGrids(hGridIndex)%noProfs
      else
        noInstances=1
      end if

      if ( vGridIndex/=0 ) then
        noSurfs=vGrids(vGridIndex)%noSurfs
      else
        noSurfs=1
      end if
      
      ! Some special cases for certain quantities
      select case (quantityType)
      case ( l_SidebandRatio, l_NoiseBandwidth, l_SystemTemperature )
        frequencyCoordinate = l_channel
        signalInfo = GetSignal(signal)
        noChans = size ( signalInfo%frequencies ) 
      case ( l_magneticField )
        frequencyCoordinate = l_xyz
        noChans = 3
      case default
      end select
      if ( .not. regular ) call announce_error ( root, no_error_code, &
        & 'Inappropriate irregular quantity request' )
      call SetupNewQuantityTemplate ( qty, noInstances=noInstances, &
        & noSurfs=noSurfs, coherent=.TRUE., stacked=.TRUE., regular=.TRUE.,&
        & noChans=noChans )
      ! ??? Note in later versions we'll need to think about channels here

      if ( hGridIndex /=0 ) then
        call CopyHGridInfoIntoQuantity ( hGrids(hGridIndex), qty )
      else                      ! Set `empty' values
        qty%phi = 0.0
        qty%geodLat = 0.0
        qty%lon = 0.0
        qty%time = 0.0
        qty%solarTime = 0.0
        qty%solarZenith = 0.0
        qty%losAngle = 0.0
      end if

      if ( quantityType == l_columnAbundance &
      & .or. &
      &   quantityType == l_boundaryPressure) then
        if ( vGridIndex/=0 ) then
          call announce_error( root, no_error_code, &
          & 'No vGrid should be specified for column abundance ' &
          & // 'boundary pressure')
        endif
        frequencyCoordinate = L_None
        qty%verticalCoordinate = L_None
      elseif ( vGridIndex /= 0 ) then
        call CopyVGridInfoIntoQuantity ( vGrids(vGridIndex), qty )
      else
        qty%surfs = 0.0
        qty%verticalCoordinate = L_None
      end if

      if ( fGridIndex /= 0 ) then
        call Allocate_test ( qty%frequencies, qty%noChans, 'qty%frequencies', &
          & ModuleName )
        qty%frequencies = fGrids(fGridIndex)%values
      end if

    end if

    ! Now fill up the remaining items, e.g. name etc.

    qty%badValue = -999.99              ! Think more about this later NJL !????
    qty%frequencyCoordinate = frequencyCoordinate
    qty%instrumentmodule = instrumentmodule
    qty%logBasis = logBasis
    qty%minValue = minValue
    qty%molecule = molecule
    qty%name = name
    qty%quantityType = quantityType
    qty%radiometer = radiometer
    qty%scaleFactor = scaleFactor
    qty%sideband = sideband
    qty%signal = signal
    qty%unit = family

    if ( majorFrame ) then
      qty%minorFrame = .false.
      qty%majorFrame = .true.
    endif
    if ( DEEBUG .or. index(switches,'qtmp') /= 0) then
      call output( 'Template name: ', advance='no' )
      call display_string ( qty%name, advance='no', ierr=ierr )
      call output( '   quantityType: ', advance='no' )
      call output( qty%quantityType, advance='no' )
      call output( '   major frame? ', advance='no' )
      call output( qty%majorFrame, advance='no' )
      call output( '   minor frame? ', advance='no' )
      call output( qty%minorFrame, advance='yes' )
      call output( '   signal string: ', advance='no' )
      call output( trim(signalString), advance='no' )
      call output( '    signal: ', advance='no' )
      if ( qty%signal < 1 ) then
        call output( ' (none available) ', advance='yes' )
      else
        call display_string ( signals(qty%signal)%name, advance='yes', ierr=ierr )
      endif
      call output( '   Instrument module: ', advance='no' )
      if ( qty%instrumentModule < 1 ) then
        call output( ' (none available) ', advance='yes' )
      else
        call display_string ( modules(qty%instrumentModule)%name, advance='yes', &
          & ierr=ierr )
      endif
!      call output( '    Molecule: ', advance='no' )
!      call display_string ( qty%Molecule, advance='yes', ierr=ierr )
      call output( '   molecule string: ', advance='no' )
      call output( trim(moleculeString), advance='yes' )
      call output ( '   noChans = ' )
      call output ( qty%noChans, advance='no' )
      call output ( ' noSurfs = ' )
      call output ( qty%noSurfs, advance='no' )
      call output ( ' noInstances = ' )
      call output ( qty%noInstances, advance='no' )
      call output ( ' instanceLen = ' )
      call output ( qty%instanceLen, advance='yes' )
    endif
    if ( toggle(gen) .and. levels(gen) > 0 ) &
      & call trace_end ( "CreateQtyTemplateFromMLSCFInfo" )
    returnStatus = error

  end function CreateQtyTemplateFromMLSCFInfo

! =====     Private Procedures     =====================================

! -----------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( WHERE, CODE, ExtraMessage, severity )
    integer, intent(in) :: WHERE   ! Tree node where error was noticed
    integer, intent(in) :: CODE    ! Code for error message
    character (LEN=*), intent(in), optional :: ExtraMessage
    integer, intent(in), optional :: severity

    if (present(severity)) then
      error = max(error,severity)
    else
      error = max(error,1)
    endif
    call output ( '***** At ' )
    if ( where > 0 ) then
      call print_source ( source_ref(where) )
    else
      call output ( '(no lcf tree available)' )
    end if
    call output ( ': ' )
    select case ( code )
    case ( badUnitMessage )
      call output ( "Incorrect or absent unit.", advance='yes' )
    case ( InappropriateQuantity )
      call output ( "A quantity inappropriate for this version is specified.", &
                    advance='yes' )
    case ( needGrid )
      call output ( "Quantity needs a vGrid or hGrid.", advance='yes' )
    case ( noQuantityType )
      call output ( "No quantity type specified.", advance='yes' )
    case ( unnecessaryGrid )
      call output ( "Quantity has an unnecessary vGrid or hGrid.", advance='yes' )
    case ( noModule )
      call output ( "No module given or deducible.", advance='yes' )
    case default
      call output ( " command caused an unrecognized programming error", advance='yes' )
    end select
    if ( present(ExtraMessage) ) then
      call output(ExtraMessage, advance='yes')
    end if
  end subroutine ANNOUNCE_ERROR

! -----------------------------------------------  ANY_GOOD_SIGNALDATA  -----
  function ANY_GOOD_SIGNALDATA ( signal, sideband, l1bInfo, Chunk )  result (answer)
  ! Read precision of signal
  ! if all values < 0.0, return FALSE
  ! if no precision data in file, return FALSE
  ! otherwise return true
  ! Arguments
    integer, intent(in) :: signal
    integer, intent(in) :: sideband
    logical             :: answer
    type (MLSChunk_T), intent(in) :: Chunk
    type (l1bInfo_T), intent(in) :: L1bInfo
  ! Private
    integer :: FileID, flag, noMAFs
    character(len=127)  :: namestring
    type (l1bData_T) :: L1BDATA
    integer :: hdfVersion

  ! Executable
    hdfVersion = mls_hdf_version(trim(l1bInfo%L1BOAFileName), LEVEL1_HDFVERSION)
    if ( hdfversion <= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Illegal hdf version for l1boa file (file missing or non-hdf?)' )
    call GetSignalName ( signal, nameString, &                   
    & sideband=sideband, noChannels=.TRUE. )                     
    nameString = AssembleL1BQtyName ( nameString, hdfVersion, .false. )
    nameString = trim(nameString) // PRECISIONSUFFIX
    fileID = FindL1BData (l1bInfo%l1bRadIDs, nameString, hdfVersion )
    if ( fileID <= 0 ) then
      answer = .false.
      return
    endif
    ! print *, 'About to read ', trim(nameString)
    ! print *, 'From Fileid ', fileID
      call ReadL1BData ( fileID , nameString, l1bData, noMAFs, flag, &
        & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
        & NeverFail= .true., hdfVersion=hdfVersion )
      if (flag == 0) then
        answer = .not. all (l1bData%DpField < 0._r8)
        call deallocate_test(l1bData%DpField, trim(nameString), ModuleName)
      else
        answer = .false.
      endif
  end function ANY_GOOD_SIGNALDATA

  ! --------------------------------  ConstructMajorFrameQuantity  -----
  subroutine ConstructMajorFrameQuantity( chunk, instrumentModule, qty, noChans, &
    & mifGeolocation )
    ! Dummy arguments
    type (MLSChunk_T), intent(in) :: CHUNK
    integer, intent(in) :: INSTRUMENTMODULE
    type (QuantityTemplate_T), intent(out) :: QTY
    type (QuantityTemplate_T), dimension(:), target :: MIFGEOLOCATION
    integer, intent(in) :: NOCHANS
    ! Local variables
    type (QuantityTemplate_T), pointer :: source
    
    ! Executable code
    source => mifGeolocation(instrumentModule)
    call SetupNewQuantityTemplate ( qty, noInstances=source%noInstances, &
      & noSurfs=1, coherent=.true., stacked=.true., regular=.true., &
      & noChans=noChans )
    qty%phi => source%phi(1:1,:)
    qty%geodLat => source%geodLat(1:1,:)
    qty%lon => source%lon(1:1,:)
    qty%time => source%time(1:1,:)
    qty%solarTime => source%solarTime(1:1,:)
    qty%solarZenith => source%solarZenith(1:1,:)
    qty%losAngle => source%losAngle(1:1,:)
    qty%mafIndex => source%mafIndex
    qty%mafCounter => source%mafCounter

    qty%verticalCoordinate = l_none
    qty%majorFrame = .true.
    qty%minorFrame = .false.

    qty%noInstancesLowerOverlap = chunk%noMAFsLowerOverlap
    qty%noInstancesUpperOverlap = chunk%noMAFsUpperOverlap
  end subroutine ConstructMajorFrameQuantity

  ! --------------------------------  ConstructMinorFrameQuantity  -----
  subroutine ConstructMinorFrameQuantity ( l1bInfo, chunk, instrumentModule, &
    & qty, noChans, regular, instanceLen, mifGeolocation )

  ! This routine constructs a minor frame based quantity.

    ! Dummy arguments
    type (L1BInfo_T), intent(in) :: l1bInfo  ! File handles for l1bdata
    type (MLSChunk_T), intent(in) :: chunk   ! The chunk under consideration
    integer, intent(in) :: instrumentModule  ! Database index
    type (QuantityTemplate_T), intent(out) :: qty ! Resulting quantity
    integer, intent(in), optional :: noChans
    logical, intent(in), optional :: regular
    integer, intent(in), optional :: instanceLen
    type (QuantityTemplate_T), intent(in), dimension(:), optional :: &
         & mifGeolocation

    ! Local parameters
    real(r8), parameter :: SIXTH = 1.0_r8 / 6.0_r8
    ! Note the similarity to CreateHGridFromMLSCFInfo:
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! IF MODIFYING THIS SECTION PLEASE TAKE CARE, SEE BELOW!
    integer, parameter :: NoL1BItemsToRead=7
    character (len=15), dimension(NoL1BItemsToRead), parameter :: &
         & L1bItemsToRead = &
         & (/"MAFStartTimeTAI","tpGeodLat      ","tpLon          ",&
         &   "tpGeodAngle    ","tpSolarZenith  ","tpSolarTime    ",&
         &   "tpLosAngle     "/)
    integer, parameter :: TransitionToModularItems = 2
    ! Entries in the above array below TransitionToModularItems are prefixed
    ! with either GHz or THz.  The layout of the above array is critically
    ! bound to the select case(l1bItem) code below.  So TAKE CARE! when
    ! modifing it.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! Local variables

    type (L1BData_T) :: l1bField
    character (len=NameLen) :: l1bItemName
    integer :: noMAFs, l1bFlag, l1bItem, mafIndex, mifIndex, hdfVersion

    ! Executable code. There are basically two cases here. If we have a
    ! MIFGeolocation argument this conveys all the geolocation for this
    ! quantity.  Otherwise, we have to read it all from the l1boa file
    ! ourselves.

    hdfVersion = mls_hdf_version(trim(l1bInfo%L1BOAFileName), LEVEL1_HDFVERSION)
    if ( hdfversion <= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Illegal hdf version for l1boa file (file missing or non-hdf?)' )
    if ( present(mifGeolocation) ) then
      ! -------------------------------------- Got mifGeolocation ------------
      if ( .not. (present(noChans)) ) &
         call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'You must supply NoChans to reuse geolocation information' )
      ! We have geolocation information, setup the quantity as a clone of that.
      qty = mifGeolocation(instrumentModule)
      if ( present ( regular ) ) qty%regular = regular
      qty%noChans = noChans
      if ( qty%regular ) then
        qty%instanceLen = qty%noChans*qty%noSurfs
      else
        qty%instanceLen = 0
      end if
      ! Increment the id counter and set the id field
      quantityTemplateCounter = quantityTemplateCounter + 1
      qty%id = quantityTemplateCounter

    else
      ! -------------------------------------- Not Got mifGeolocation ------------
      ! We have no geolocation information, we have to read it ourselves
      ! from the L1BOA file.

      ! First we read xxGeodalt to get the size of the quantity.
      if ( IsModuleSpacecraft(instrumentModule) ) then
        l1bItemName = 'scGeocAlt' 
      else
        call GetModuleName ( instrumentModule, l1bItemName )
        l1bItemName = TRIM(l1bItemName) // "." // "tpGeodAlt"
      endif
      l1bItemName = AssembleL1BQtyName ( l1bItemName, hdfVersion, .false. )

      call ReadL1BData ( l1bInfo%l1boaid, l1bItemName, l1bField, noMAFs, &
        & l1bFlag, firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
        & hdfVersion=hdfVersion )
      if ( l1bFlag==-1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_L1BRead//l1bItemName )
      
      ! Now noMAFs qty%noInstances, l1bField%maxMIFs is no surfs.
      !print *, 'About to SetupNewQuantityTemplate'
      !print *, 'noInstances ', noMAFs
      !print *, 'noSurfs ', l1bField%maxMIFs
      !if ( present(noChans) ) then
      !  print *, 'noChans ', noChans
      !else
      !  print *, 'noChans not present'
      !endif
      call SetupNewQuantityTemplate ( qty, noInstances=noMAFs, &
        & noSurfs=l1bField%maxMIFs, noChans=noChans, coherent=.false., &
        & stacked=.false., regular=regular, instanceLen=instanceLen, &
        & minorFrame=.true. )
      qty%noInstancesLowerOverlap = chunk%noMAFsLowerOverlap
      qty%noInstancesUpperOverlap = chunk%noMAFsUpperOverlap

      if ( .not. IsModuleSpacecraft(instrumentModule) ) then
        call GetModuleName ( instrumentModule, l1bItemName )
        l1bItemName = TRIM(l1bItemName) // "." // "tpGeodAlt"
        
        l1bItemName = AssembleL1BQtyName ( l1bItemName, hdfVersion, .false. )
        call ReadL1BData ( l1bInfo%l1boaid, l1bItemName, l1bField, noMAFs, &
          & l1bFlag, firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
          & hdfVersion=hdfVersion )
        if ( l1bFlag==-1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_L1BRead//l1bItemName )
        
        ! Now noMAFs qty%noInstances, l1bField%maxMIFs is no surfs.
        
        call SetupNewQuantityTemplate ( qty, noInstances=noMAFs, &
          & noSurfs=l1bField%maxMIFs, noChans=noChans, coherent=.FALSE., &
          & stacked=.FALSE., regular=regular, instanceLen=instanceLen, &
          & minorFrame=.TRUE. )
        
        qty%noInstancesLowerOverlap = chunk%noMAFsLowerOverlap
        qty%noInstancesUpperOverlap = chunk%noMAFsUpperOverlap
        
        ! Now we're going to deal with a VGrid for this quantity
        
        qty%verticalCoordinate = l_geodAltitude
        qty%surfs = l1bField%dpField(1,:,:)  ! Vert coord is tpGeodAlt read above.
        qty%mafCounter = l1bField%counterMAF
        do mafIndex = chunk%firstMAFIndex, chunk%lastMAFIndex
          qty%mafIndex(mafIndex-chunk%firstMAFIndex+1) = mafIndex
        end do
        call DeallocateL1BData(l1bfield)

        ! Now we're going to fill in the hGrid information
        
        do l1bItem = 1, NoL1BItemsToRead
          ! Get the name of the item to read
          l1bItemName = l1bItemsToRead(l1bItem)
          if ( l1bItem >= TransitionToModularItems ) then
            call GetModuleName ( instrumentModule, l1bItemName )
            l1bItemName = trim(l1bItemName)//'.'//l1bItemsToRead(l1bItem)
          else
            l1bItemName = l1bItemsToRead(l1bItem)
          end if
          
          ! Read it from the l1boa file
          l1bItemName = AssembleL1BQtyName ( l1bItemName, hdfVersion, .false. )
          call ReadL1BData ( l1bInfo%l1boaid, l1bItemName, l1bField, noMAFs, &
            & l1bFlag, firstMAF=chunk%firstMafIndex, &
            & lastMAF=chunk%lastMafIndex, hdfVersion=hdfVersion )
          if ( l1bFlag == -1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & MLSMSG_L1BRead//l1bItemName )
          
          ! Now we have to save this field in the quantity data.
          ! This is rather a kludgy way of doing it but this worked out the
          ! least boring way to write the code.  See the definition of
          ! L1BItemsToRead above for reference.
          
          select case(l1bItem)
          case ( 1 )
            !??? For time we have to do something a little more complicated.
            !??? This is a real kludge, and we have to find a way
            !??? to do it better in 0.5. Probably simply have time as a minor
            !??? frame quantity in L1, or MIF duration. !???????
            !??? Also note that it fills in times even for non existant MIFs
            do mafIndex = 1, noMAFs
              do mifIndex = 1, qty%noSurfs 
                qty%time(mifIndex,mafIndex) = &
                  & l1bField%dpField(1,1,mafIndex) + &
                  & (mifIndex-1) * sixth
              end do
            end do
          case ( 2 )
            qty%geodLat = l1bField%dpField(1,:,:)
          case ( 3 )
            qty%lon = l1bField%dpField(1,:,:)
          case ( 4 )
            qty%phi = l1bField%dpField(1,:,:)
          case ( 5 )
            qty%solarZenith = l1bField%dpField(1,:,:)
          case ( 6 )
            qty%solarTime = l1bField%dpField(1,:,:)
          case ( 7 )
            qty%losAngle = l1bField%dpField(1,:,:)
          end select

          call DeallocateL1BData ( l1bField )
        end do                          ! Loop over l1b quantities
      else                              ! Spacecraft module
        ! just zero stuff out.
        qty%surfs = 0.0
        qty%phi = 0.0
        qty%geodLat = 0.0
        qty%lon = 0.0
        qty%time = 0.0
        qty%solarTime = 0.0
        qty%solarZenith = 0.0
        qty%losAngle = 0.0
        
        qty%mafCounter = l1bField%counterMAF
        do mafIndex = chunk%firstMAFIndex, chunk%lastMAFIndex
          qty%mafIndex(mafIndex-chunk%firstMAFIndex+1) = mafIndex
        end do
        
        call DeallocateL1BData(l1bfield)
        
      end if
      
    end if
    qty%frequencyCoordinate = L_None
    qty%instrumentModule = instrumentModule

    ! In later versions we'll probably need to think about FILL_VALUEs and
    ! setting things to the badData flag.

  end subroutine ConstructMinorFrameQuantity

  ! ------------------------------------ ForgeMinorFrames --------------
  subroutine ForgeMinorFrames ( root, chunk, mifGeolocation )
    ! This routine is used when we're in the mode of creating l2pc files
    ! and we want to invent a set of minor frame quantities with no
    ! reference to the l1 files
    
    ! Dummy arguments
    integer, intent(in) :: ROOT         ! Tree vertex
    type (MLSChunk_T), intent(in) :: CHUNK ! This chunk
    type (QuantityTemplate_T), dimension(:), intent(inout) :: MIFGEOLOCATION

    ! Local variables
    integer :: GEODANGLENODE            ! Tree vertex
    integer :: I                        ! Loop counter
    integer :: INSTRUMENTMODULE         ! Database index
    integer :: KEY                      ! Tree vertex
    integer :: MAF                      ! Loop counter
    integer :: NODE                     ! A tree node
    integer :: NOMAFS                   ! Dimension
    integer :: NOMIFS                   ! Dimension
    integer :: PARAM                    ! Loop counter
    integer :: SOLARTIMENODE            ! Tree vertex
    integer :: SOLARZENITHNODE          ! Tree vertex
    integer :: SON                      ! Tree vertex
    integer :: UNITS                    ! Units for node

    real(r8), dimension(:,:), pointer :: VALUES ! An array to fill
    integer, dimension(2) :: EXPR_UNITS ! From tree
    real(r8), dimension(2) :: EXPR_VALUE ! From tree

    ! Executable code
    solarTimeNode = 0
    solarZenithNode = 0
    geodAngleNode = 0
    do i = 2, nsons( root )
      son = subtree(i,root)
      key = subtree(1,son)

      select case ( decoration(key) )
      case ( f_module )
        instrumentModule = decoration(decoration(subtree(2,son)))
      case ( f_geodAngle )
        geodAngleNode = son
      case ( f_solarTime )
        solarTimeNode = son
      case ( f_solarZenith )
        solarZenithNode = son
      case ( f_noMIFs )
        call expr ( subtree(2,son), expr_units, expr_value )
        noMIFs = expr_value(1)
      case default ! Can't get here if tree_checker works correctly
      end select
    end do

    ! Work out how many MAFs we're after
    if ( geodAngleNode /= 0 ) noMAFs = nsons ( geodAngleNode ) - 1
    if ( solarTimeNode /= 0 ) then
      if ( noMAFs /= 0 ) then
        if ( nsons ( solarTimeNode ) - 1 /= noMAFs ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Inconsistent explicit Forge definition' )
      else
        noMAFs = nsons ( solarTimeNode ) - 1
      end if
    end if
    if ( solarZenithNode /= 0 ) then
      if ( noMAFs /= 0 ) then
        if ( nsons ( solarZenithNode ) - 1 /= noMAFs ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Inconsistent explicit Forge definition' )
      else
        noMAFs = nsons ( solarZenithNode ) - 1
      end if
    end if

    ! Setup the minor frame quantity template
    call SetupNewQuantityTemplate ( mifGeolocation(instrumentModule), &
      & noInstances=noMAFs, noSurfs=noMIFs, noChans=1,&
      & coherent=.false., stacked=.false., regular=.true.,&
      & minorFrame=.true. )

    ! Put zeros etc. in the appropriate places, may overwrite these later
    mifGeolocation(instrumentModule)%instrumentModule = instrumentModule
    mifGeolocation(instrumentModule)%noInstancesLowerOverlap = 0
    mifGeolocation(instrumentModule)%noInstancesUpperOverlap = 0
    mifGeolocation(instrumentModule)%verticalCoordinate = l_none
    mifGeolocation(instrumentModule)%time = 0.0
    mifGeolocation(instrumentModule)%phi = 0.0
    mifGeolocation(instrumentModule)%geodLat = 0.0
    mifGeolocation(instrumentModule)%solarTime = 0.0
    mifGeolocation(instrumentModule)%solarZenith = 0.0
    mifGeolocation(instrumentModule)%lon = 0.0
    mifGeolocation(instrumentModule)%losAngle = 0.0
    mifGeolocation(instrumentModule)%mafIndex = 0.0
    mifGeolocation(instrumentModule)%mafCounter = 0.0
    mifGeolocation(instrumentModule)%surfs = 0.0

    ! Loop over the parameters we might have
    do param = 1, 3
      select case ( param )
      case ( 1 )
        values => mifGeolocation(instrumentModule)%phi
        node = geodAngleNode
        units = phyq_angle
      case ( 2 )
        values => mifGeolocation(instrumentModule)%solarTime
        node = solarTimeNode
        units = phyq_time
      case ( 3 )
        values => mifGeolocation(instrumentModule)%solarZenith
        node = solarZenithNode
        units = phyq_angle
      end select

      if ( node /= 0 ) then
        do maf = 1, noMAFs
          call expr ( subtree ( maf+1, node), expr_units, expr_value )
          if ( all ( expr_units(1) /= (/ phyq_dimensionless, units /) ) ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Invalid units for explicit hGrid' )
          values(:,maf) = expr_value(1)
        end do
      end if
    end do
    ! Make the latitude the same as the geod angle
    ! This is a bit of a hack, but it should be OK in all cases.
    mifGeolocation(instrumentModule)%geodLat = &
      & mifGeolocation(instrumentModule)%phi
    
  end subroutine ForgeMinorFrames

  ! ----------------------------------  CopyHGridInfoIntoQuantity  -----
  subroutine CopyHGridInfoIntoQuantity ( hGrid, qty )

  ! This routine copies HGrid information into an already defined quantity

    ! Dummy arguments
    type (hGrid_T), intent(in) :: hGrid
    type (QuantityTemplate_T), intent(inout) :: QTY

    ! Executable code

    if ( qty%noInstances/=hGrid%noProfs ) call MLSMessage ( MLSMSG_Error,&
      & ModuleName, "Size of HGrid not compatible with size of quantity" )

    if ( .NOT. qty%stacked ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Cannot copy hGrids into unstacked quantities")

    qty%phi(1,:) = hGrid%phi
    qty%geodLat(1,:) = hGrid%geodLat
    qty%lon(1,:) = hGrid%lon
    qty%time(1,:) = hGrid%time
    qty%solarTime(1,:) = hGrid%solarTime
    qty%solarZenith(1,:) = hGrid%solarZenith
    qty%losAngle(1,:) = hGrid%losAngle
    qty%noInstancesLowerOverlap = hGrid%noProfsLowerOverlap
    qty%noInstancesUpperOverlap = hGrid%noProfsUpperOverlap

  end subroutine CopyHGridInfoIntoQuantity

  ! ----------------------------------  CopyVGridInfoIntoQuantity  -----
  subroutine CopyVGridInfoIntoQuantity ( vGrid, qty )

  ! This similar routine copies VGrid information into an already 
  ! defined quantity

    ! Dummy arguments
    type (VGrid_T), intent(in) :: vGrid
    type (QuantityTemplate_T), intent(inout) :: qty

    ! Executable code

    if ( vGrid%noSurfs /= qty%noSurfs ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, "Size of vGrid not compatible with size of quantity" )

    if ( .NOT. qty%coherent ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Cannot copy vGrid information into incoherent quantities")

    qty%verticalCoordinate = vGrid%verticalCoordinate
    qty%surfs(:,1) = vGrid%surfs

  end subroutine CopyVGridInfoIntoQuantity


!=============================================================================
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module ConstructQuantityTemplates
!=============================================================================

!
! $Log$
! Revision 2.87  2003/02/12 21:53:48  pwagner
! Fix to errant ANY_GOOD_SIGNALDATA in case of hdf5 files
!
! Revision 2.86  2003/02/07 03:42:50  vsnyder
! Fill NATURAL_UNITS on first call only -- it's a SAVE variable now
!
! Revision 2.85  2003/02/07 03:38:05  vsnyder
! Cosmetic change
!
! Revision 2.84  2003/02/07 00:41:41  livesey
! Bug fix/workaround
!
! Revision 2.83  2003/02/06 23:31:00  livesey
! New approach to Forge
!
! Revision 2.82  2003/01/14 00:40:29  pwagner
! Moved GetQuantityAttributes to L2AUXData
!
! Revision 2.81  2003/01/09 00:09:15  pwagner
! routine GetQuantityAttributes added
!
! Revision 2.80  2003/01/08 23:50:44  livesey
! Added irregular argument
!
! Revision 2.79  2003/01/07 23:46:53  livesey
! Added magnetic field stuff
!
! Revision 2.78  2002/12/11 22:17:05  pwagner
! Added error checks on hdf version
!
! Revision 2.77  2002/11/26 23:37:50  livesey
! Better handling of major frame quantities
!
! Revision 2.76  2002/11/22 12:16:08  mjf
! Added nullify routine(s) to get round Sun's WS6 compiler not
! initialising derived type function results.
!
! Revision 2.75  2002/11/13 01:05:03  pwagner
! Actually reads hdf5 radiances
!
! Revision 2.74  2002/10/08 17:36:19  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.73  2002/09/26 18:03:06  livesey
! Changed extinction to a vmr
!
! Revision 2.72  2002/09/25 20:07:55  livesey
! Made -g less verbose
!
! Revision 2.71  2002/09/24 21:37:44  livesey
! Added minValue stuff
!
! Revision 2.70  2002/09/24 00:27:16  pwagner
! Wont bomb if no l1brads; nor whine if no good signals
!
! Revision 2.69  2002/09/18 22:48:32  pwagner
! Chooses signals first band with any good data
!
! Revision 2.68  2002/08/04 16:01:19  mjf
! Added some nullify statements for Sun's rubbish compiler.
!
! Revision 2.67  2002/06/14 16:39:49  livesey
! Orbital inclination now minor frame
!
! Revision 2.66  2002/06/04 22:07:35  livesey
! Added phiTan as a state vector element
!
! Revision 2.65  2002/05/22 19:06:32  jonathan
! added units for cloudextinction(m-1), totalextinction(m-1), and
! massmeandiameterice(micron m) as dimentionless for now, may define
! more clear later
!
! Revision 2.64  2002/05/14 00:27:42  livesey
! New code for system temperatures and noise bandwidths
!
! Revision 2.63  2002/05/07 20:02:54  livesey
! Added noise bandwidth
!
! Revision 2.62  2002/04/10 17:44:22  pwagner
! Added rhi quantity (but is this enough?)
!
! Revision 2.61  2002/03/19 00:51:32  pwagner
! Added new scVel quantity types
!
! Revision 2.60  2002/02/09 21:35:52  livesey
! Added optical depth stuff
!
! Revision 2.59  2001/11/08 00:13:38  livesey
! Sorted out extinction stuff
!
! Revision 2.58  2001/10/31 19:07:25  livesey
! Added fGrid stuff
!
! Revision 2.57  2001/10/12 23:15:05  pwagner
! Fixed biggest erros in diagnostic quantity templates
!
! Revision 2.56  2001/10/02 23:12:50  pwagner
! More chi^2 fixes
!
! Revision 2.55  2001/10/02 20:50:54  livesey
! No longer asserts baseline to be minor frame
!
! Revision 2.54  2001/09/17 23:11:50  pwagner
! Tiny changes for chi^2
!
! Revision 2.53  2001/09/17 21:58:50  livesey
! Added allocate of frequencies if needed
!
! Revision 2.52  2001/09/14 23:30:33  pwagner
! Now constructs major frame quantity templates
!
! Revision 2.51  2001/08/01 00:04:29  dwu
! add qty%frequencies = VGrids(sGridIndex)%surfs for quantity l_losTransFunc
!
! Revision 2.50  2001/07/30 23:28:38  pwagner
! Added columnAbundances scaffolding--needs fleshing out
!
! Revision 2.49  2001/07/26 17:34:25  jonathan
! add DTcir, etc, jonathan
!
! Revision 2.48  2001/07/19 22:17:44  jonathan
! add cloud stuff , jonathan/dwu
!
! Revision 2.47  2001/07/19 17:42:31  dwu
! add sGrid field
!
! Revision 2.46  2001/07/18 23:17:30  dwu
! rename l_radiusofearth as l_earthradius
!
! Revision 2.45  2001/07/18 18:42:19  dwu
! add radiusofearth quantity type
!
! Revision 2.44  2001/07/17 23:23:05  dwu
! make l_losTransFunc as non-minorframe but minorframe-like quantity
!
! Revision 2.43  2001/07/16 18:24:45  dwu
! add feature for losTransFunc type of quantities
!
! Revision 2.42  2001/07/13 18:41:59  dwu
! fix problem after adding losTransFunc
!
! Revision 2.41  2001/07/13 18:10:03  dwu
! add quantity losTransFunc
!
! Revision 2.40  2001/07/11 21:40:00  livesey
! More bug fixes
!
! Revision 2.39  2001/07/10 23:45:16  jonathan
! added cloudicedensity and template for cloudsfwm, paul/jonathan
!
! Revision 2.38  2001/07/09 22:37:23  livesey
! Embarassing memory leak caught!  It's our old friend
! mifGeolocation again.  I'm going to regret trying
! to be this clever!
!
! Revision 2.37  2001/05/31 19:53:56  livesey
! Whoops, debug stuff left in.
!
! Revision 2.36  2001/05/30 23:59:51  livesey
! Thought I'd made this change already.  I'm confused
!
! Revision 2.35  2001/05/30 23:55:28  livesey
! Previous one was debug version, this is correct one.
!
! Revision 2.34  2001/05/30 23:53:01  livesey
! For new version of L1Bdata
!
! Revision 2.33  2001/05/29 23:21:35  livesey
! Changed l_orbitincline to l_orbitinclination
!
! Revision 2.32  2001/05/26 00:20:20  livesey
! Cosmetic changes
!
! Revision 2.31  2001/05/10 01:08:53  livesey
! Changed hGrids and vGrids to pointers, rather than intent(in)
! to allow them to be empty.
!
! Revision 2.30  2001/05/03 23:08:36  livesey
! Added stuff to support scan model items.
!
! Revision 2.29  2001/05/03 20:30:09  vsnyder
! Add a 'nullify' and some cosmetic changes
!
! Revision 2.28  2001/04/26 02:45:25  vsnyder
! Fix up CVS stuff
!
! Revision 2.27  2001/04/26 02:44:17  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.26  2001/04/25 20:33:07  livesey
! Minor bug fix, Forge now also zeros surfs.
!
! Revision 2.25  2001/04/25 19:29:49  livesey
! Fixed bug in forge, now sets mafCounter and mafIndex correctly.
!
! Revision 2.24  2001/04/25 00:01:23  livesey
! Bug fix, no default units for scGeocAlt
!
! Revision 2.23  2001/04/24 22:21:17  livesey
! Gave up on latitude for forge
!
! Revision 2.22  2001/04/23 23:25:10  livesey
! Fixed bug in forge
!
! Revision 2.21  2001/04/20 23:11:39  livesey
! Added forge stuff for minor frames
!
! Revision 2.20  2001/04/19 20:30:06  livesey
! Added specific stuff for sideband ratio
!
! Revision 2.19  2001/04/12 23:25:29  vsnyder
! Give "Sideband" an initial value
!
! Revision 2.18  2001/04/12 21:41:42  livesey
! Signal now a string.
!
! Revision 2.17  2001/04/10 22:27:47  vsnyder
! Nullify explicitly instead of with <initialization> so as not to give
! pointers the SAVE attribute.  <initialization> is NOT executed on each
! entry to a procedure.
!
! Revision 2.16  2001/04/07 01:50:48  vsnyder
! Move some of VGrid to lib/VGridsDatabase.  Move ForwardModelConfig_T and
! some related stuff to fwdmdl/ForwardModelConfig.
!
! Revision 2.15  2001/03/28 23:48:13  livesey
! Bug fixes zero out some stuff.
!
! Revision 2.14  2001/03/28 18:33:19  livesey
! Fixed bug with logBasis (wasn't initialised!)
!
! Revision 2.13  2001/03/21 02:13:30  livesey
! Bug with logBasis, put in a work around. Will need to fix later
!
! Revision 2.12  2001/03/17 02:23:55  livesey
! Added logBasis (and set value for badData)
!
! Revision 2.11  2001/03/15 21:07:47  vsnyder
! Cross-references between databases are by database index, not tree index
!
! Revision 2.10  2001/03/15 18:41:17  livesey
! Tidied up the frequency coordinate stuff.
!
! Revision 2.9  2001/03/08 21:49:26  livesey
! Added elev_offset
!
! Revision 2.8  2001/03/03 00:08:09  livesey
! Lots of changes mostly with minor frame quantities
!
! Revision 2.7  2001/03/02 01:28:23  livesey
! New quantity types etc.
!
! Revision 2.6  2001/02/28 01:17:04  livesey
! Interim version, on the way to using proper signals stuff
!
! Revision 2.5  2001/02/22 23:37:24  livesey
! Really removed all references to firstIndexChannel
!
! Revision 2.4  2001/02/21 01:09:00  livesey
! Allowed for quantities with no h/v grid
!
! Revision 2.3  2001/02/20 18:43:50  livesey
! Removed all references to firstIndexChannel
!
! Revision 2.2  2001/02/14 00:12:45  livesey
! Removed firstIndexChannel
!
! Revision 2.1  2001/02/09 00:38:22  livesey
! Various updates
!
! Revision 2.0  2000/09/05 18:57:02  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:03  vsnyder
! Initial entry
!
! Revision 1.13  2000/05/17 23:19:07  lungu
! Added "." between MLSInstrumentModuleNames and l1bItemName.
! Made hGridIndex=0 and vGridIndex=0 upon entry, so it does not "inherit" attributes from previous call.
! Made caseInsensitive=.TRUE. for all searches.
! Added type for QTY_Gph.
! Made stacked=.TRUE. so that CopyHGridInfoIntoQuantity works.
!
