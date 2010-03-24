module CFM_QuantityTemplate_m
   use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test, Test_Allocate
   use FGrid, only: fGrid_T
   use HGridsDatabase, only: hGrid_T
   use VGridsDatabase, only: VGrids, VGrid_T
   use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
   use QuantityTemplates, only: QuantityTemplate_T, NullifyQuantityTemplate, &
                                DestroyQuantityTemplateDatabase, &
                                AddQuantityTemplateToDatabase, Dump
   use Init_Tables_Module, only: FIRST_LIT, LAST_LIT, L_ADOPTED, &
      L_BASELINE, L_BOUNDARYPRESSURE, L_CALSIDEBANDFRACTION, &
      L_CHISQBINNED, L_CHISQCHAN, L_CHISQMMAF, L_CHISQMMIF, L_CLOUDICE, &
      L_CLOUDINDUCEDRADIANCE, L_CLOUDEXTINCTION, L_CLOUDMINMAX, L_CLOUDRADSENSITIVITY, &
      L_CLOUDTEMPERATURE, L_CLOUDWATER, L_COLUMNABUNDANCE, &
      L_DNWT_ABANDONED, L_DNWT_AJN, L_DNWT_AXMAX, L_DNWT_CAIT, &
      L_DNWT_CHISQMINNORM, L_DNWT_CHISQNORM, L_DNWT_CHISQRATIO, &
      L_DNWT_COUNT, L_DNWT_DIAG, L_DNWT_DXDX, L_DNWT_DXDXL, &
      L_DNWT_DXN, L_DNWT_DXNL, L_DNWT_FLAG, L_DNWT_FNMIN, &
      L_DNWT_FNORM, L_DNWT_GDX, L_DNWT_GFAC, &
      L_DNWT_GRADN, L_DNWT_SQ, L_DNWT_SQT,&
      L_EARTHRADIUS, L_EARTHREFL, L_ECRTOFOV, L_EFFECTIVEOPTICALDEPTH, &
      L_ELEVOFFSET, L_EXTINCTION, L_EXTINCTIONV2, &
      L_FIELDAZIMUTH, L_FIELDELEVATION, L_FIELDSTRENGTH, &
      L_GEODALTITUDE, L_GEOLOCATION, L_GPH, L_HEIGHTOFFSET, L_ISOTOPERATIO, L_IWC, &
      L_JACOBIAN_COLS, L_JACOBIAN_ROWS, &
      L_L1BMAFBASELINE, L_L1BMIF_TAI, L_LIMBSIDEBANDFRACTION, &
      L_LineCenter, L_LineWidth, L_LineWidth_TDep, &
      L_LOSTRANSFUNC, L_LOSVEL, &
      L_MASSMEANDIAMETERICE, L_MASSMEANDIAMETERWATER, L_MAGNETICFIELD, &
      L_MIFDEADTIME, L_ZETA, L_XYZ, L_MATRIX3X3, L_Channel, &
      L_NOISEBANDWIDTH, L_NORADSPERMIF, L_NORADSBINNED, &
      L_NUMGRAD, L_NUMJ, L_NUMNEWT, L_OPTICALDEPTH, L_ORBITINCLINATION, &
      L_PHASETIMING, L_PHITAN, L_PTAN, L_QUALITY, L_RADIANCE, &
      L_REFGPH, L_REFLTEMP, L_REFLTRANS, L_REFLREFL, L_REFLSPILL, &
      L_RHI, L_SINGLECHANNELRADIANCE, L_SIZEDISTRIBUTION, &
      L_SCANRESIDUAL, L_SCATTERINGANGLE, L_SCECI, L_SCVEL, L_SCVELECI, &
      L_SCVELECR, L_SCGEOCALT, L_SPACERADIANCE, L_STATUS, &
      L_STRAYRADIANCE, L_SurfaceHeight, L_SURFACETYPE, L_SYSTEMTEMPERATURE, &
      L_TEMPERATURE, L_TNGTECI, L_TNGTGEODALT, L_TNGTGEOCALT, &
      L_TOTALPOWERWEIGHT, L_TSCAT, L_VMR, L_NONE, &
      PHYQ_ANGLE, PHYQ_COLMABUNDANCE, &
      PHYQ_DIMENSIONLESS, PHYQ_EXTINCTION, PHYQ_FREQUENCY,&
      PHYQ_GAUSS, PHYQ_IceDensity, PHYQ_LENGTH, &
      PHYQ_PRESSURE, PHYQ_TEMPERATURE, PHYQ_TIME, PHYQ_VELOCITY, &
      PHYQ_VMR, PHYQ_ZETA
    use Intrinsic, only: LIT_INDICES
    use Output_M, only: OUTPUT
    use String_Table, only: DISPLAY_STRING
    use MLSCommon, only: r8
    use MLSSignals_m, only: GetModuleIndex
    use Molecules, only: GetMoleculeIndex
    use Chunks_m, only: MLSChunk_T
    use MLSCommon, only: MLSFile_T
    use ConstructQuantityTemplates, only : AnyGoodSignalData, GetQtyTypeIndex, &
            ConstructMajorFrameQuantity, ConstructMinorFrameQuantity, &
            InitQuantityTemplates, &
            unitstable, propertyTable, noProperties, &
            p_majorFrame, p_minorFrame, p_mustBeZeta, p_fGrid, p_hGrid, &
            p_vgrid, p_radiometer, p_radiometerOptional, p_scModule, &
            p_signal, p_signalOptional, p_xyz, p_matrix3x3, p_flexibleVHGrid, &
            p_suppressChannels, p_module, p_molecule, p_fGridOptional
    use Construct, only: ConstructMIFGeolocation
    use Parse_Signal_m, only: PARSE_SIGNAL
    use MLSSignals_m, only: GetModuleFromSignal, GetRadiometerFromSignal, &
                            GetSignal, Signal_T, IsModuleSpacecraft, &
                            GetRadiometerIndex, GetModuleFromRadiometer

   implicit none

   public :: CreateQtyTemplate, DestroyQuantityTemplateDatabase
   public :: AddQuantityTemplateToDatabase
   public :: QuantityTemplate_T, Dump
   public :: InitQuantityTemplates, ConstructMIFGeolocation

   private

   character(len=20) :: ModuleName = "CFM_QuantityTemplate"

   contains

   ! This subroutine design is prone to bugs, I'll fix it when I got some time. -haley

   ! Creating a quantity based on the optional input this subroutine is provided with.
   type(QuantityTemplate_T) function CreateQtyTemplate (type, filedatabase, chunk, &
        avgrid, ahgrid, afgrid, qInstModule, qMolecule, qLogBasis, qMinValue, qSignal, &
        qRadiometer, qBadValue, mifGeolocation) result(qty)
      character(len=*), intent(in) :: type
      type (MLSFile_T), dimension(:), pointer, optional ::     FILEDATABASE
      type (MLSChunk_T), intent(in), optional :: Chunk
      type(VGrid_T), intent(in), optional :: avgrid
      type(HGrid_T), intent(in), optional :: ahgrid
      type(FGrid_T), intent(in), optional :: afgrid
      character(len=*), optional :: qInstModule
      character(len=*), optional :: qMolecule
      character(len=*), optional :: qSignal
      logical, optional :: qLogBasis
      real(r8), optional :: qMinValue
      character(len=*), optional :: qRadiometer
      real(r8), optional :: qBadValue
      type (QuantityTemplate_T), dimension(:), intent(in), optional, target :: &
      & MifGeolocation ! TARGET attribute needed so pointers to MifGeolocation
                       ! created at lower levels in the call tree don't become
                       ! undefined when those procedures return.

      logical, dimension(noProperties) :: PROPERTIES ! Properties for this quantity type
      character(len=127) :: signalString
      integer, dimension(:), pointer :: SignalInds ! From parse signal
      logical, pointer :: Channels(:)     ! From Parse_Signal
      integer :: s_index, temp
      type(Signal_T) :: signalInfo
      integer :: quantityType, noChans, sideband, signal, radiometer

      ! Executables
      call NullifyQuantityTemplate(qty)
      qty%instrumentModule = 0
      call GetQtyTypeIndex(type, quantityType)
      qty%name = 0
      qty%fGridIndex = 0 ! we're not getting fgrid from fgridDatabase
      qty%hGridIndex = 0 ! we're not getting hgrid from hgridDatabase
      qty%vGridIndex = 0 ! we're not getting vgrid from vgridDatabase
      qty%logBasis = .false.
      qty%minValue = -huge(0.0_r8)
      qty%molecule = 0
      noChans = 1
      qty%frequencyCoordinate = l_none
      radiometer = 0
      qty%reflector = 0
      qty%regular = .true.
      sideband = 0
      signal = 0
      qty%vGridIndex = 0
      qty%noInstances = 1
      qty%noSurfs = 1
      qty%regular = .true. ! irregular quantity not supported
      qty%coherent = .true.
      qty%stacked = .true.
      qty%minorFrame = .false.
      qty%majorFrame = .false.
      qty%sharedVGrid = .true.
      qty%sharedHGrid = .true.
      qty%sharedFGrid = .true.
      qty%badValue = huge(0.0_r8)   !Default bad value
      ! These field are not used in cfm
      qty%noInstancesLowerOverlap = 0
      qty%noInstancesUpperOverlap = 0

      properties = propertyTable(:, quantityType)

      !Sanity check
      if (.not. properties (p_flexibleVHGrid)) then
         if ( present ( ahGrid ) .neqv. properties ( p_hGrid ) ) &
           & call Announce_error ( 0, trim ( merge ( 'unexpected', 'need      ', &
           & present(ahGrid) ) ) // ' hGrid for quantity type ', quantityType, &
           & severity='nonfatal' )
         if ( present ( avGrid ) .neqv. properties ( p_vGrid ) ) &
           & call Announce_error ( 0, trim ( merge ( 'unexpected', 'need      ', &
           & present(avGrid) ) ) // ' vGrid for quantity type ', quantityType )
      end if
      if ( present ( qMolecule ) .neqv. properties ( p_Molecule ) ) &
         & call Announce_error ( 0, trim ( merge ( 'unexpected', 'need      ', &
         & present(qMolecule) ) ) // ' molecule for quantity type ', quantityType )
      if ( present ( qInstModule ) .neqv. ( properties ( p_Module ) .or. properties ( p_scModule) ) ) then
         call Announce_error ( 0, trim ( merge ( 'unexpected', 'need      ', &
         & present(qInstModule) ) ) // ' module for quantity type ', quantityType )
      end if
      if ( .not. properties ( p_signalOptional ) ) then
         if ( present ( qSignal ) .neqv. properties ( p_Signal ) ) &
            & call Announce_error ( 0, trim ( merge ( 'unexpected', 'need      ', &
            & present(qSignal) ) ) // ' signal for quantity type ', quantityType )
      end if
      if ( properties ( p_mustBeZeta ) .and. present(avGrid) ) then
         if ( avGrid%verticalCoordinate /= l_zeta ) &
           & call Announce_error ( 0, 'Expecting log pressure coordinates for', &
           & qty%quantityType )
      end if
      if (.not. properties(p_fGridOptional)) then
         if ( present ( aFGrid ) .neqv. properties ( p_fGrid ) ) &
           & call Announce_error ( 0, trim ( merge ( 'unexpected', 'need      ', &
           & present(aFGrid) ) ) // ' fGrid for quantity type ', quantityType )
      end if
      if ( .not. properties ( p_radiometerOptional ) ) then
         if ( present ( qRadiometer ) .neqv. properties ( p_Radiometer ) ) &
           & call Announce_error ( 0, trim ( merge ( 'unexpected', 'need      ', &
           & present(qRadiometer) ) ) // ' radiometer for quantity type ', quantityType )
      end if

     ! Set up a non major/minor frame quantity
      qty%minorFrame = .not. present(ahgrid)

      if (present(qInstModule)) then
         call GetModuleIndex (qInstModule, qty%instrumentModule)
      end if
      if (present (qRadiometer)) then
         call GetRadiometerIndex(qRadiometer, radiometer)
         ! Every radiometer pairs with only one instrument module
         qty%instrumentModule = GetModuleFromRadiometer(radiometer)
         ! If the user happen to pass in the wrong instrumentModule,
         ! silently correct it, for simplicity
      end if

      if (present(qLogBasis)) qty%logBasis = qLogBasis

      if (present(qMinValue)) qty%minValue = qMinValue

     if (present(qSignal)) then
          if (.not. (present(filedatabase) .and. present(chunk))) &
            call MLSMessage(MLSMSG_Error, ModuleName, &
            'Need filedatabase and chunk to check for good signal data')
          nullify ( channels, signalInds )
          signalString = qSignal
          !??? Here we would do intelligent stuff to work out which bands
          !??? are present, for the moment choose the first
          call parse_Signal ( signalString, signalInds, &
            & sideband=sideband, channels=channels )
          if ( .not. associated(signalInds) ) then ! A parse error occurred
            call MLSMessage ( MLSMSG_Error, ModuleName,&
              & 'Unable to parse signal string' )
          end if
          if ( size(signalInds) == 1 .or. .not. associated(filedatabase) ) then
            signal = signalInds(1)
          else
            ! Seek a signal with any precision values !< 0
            do s_index=1, size(signalInds)
              if ( AnyGoodSignalData ( signalInds(s_index), sideband, &
                & filedatabase, chunk) ) exit
            end do
            if ( s_index > size(signalInds) ) then
              signal = signalInds(1)
            else
              signal = signalInds(s_index)
            end if
          end if
          call deallocate_test ( signalInds, 'signalInds', ModuleName )
          ! if the user has pass in the wrong instrumentModule or radiometer
          ! just silently correct it, for simplicity.
          qty%instrumentModule = GetModuleFromSignal(signal)
          radiometer = GetRadiometerFromSignal(signal)
     end if

      ! After all the possible places where instrumentModule can be set
      if ( properties ( p_scModule ) ) then
         if ( .not. IsModuleSpacecraft ( qty%instrumentModule ) ) &
           & call Announce_error ( 0, 'Module must be spacecraft' )
      end if

      ! Now establish the frequency coordinate system
      if (present(afgrid)) then
         qty%frequencyCoordinate = afgrid%frequencyCoordinate
         qty%frequencies => afgrid%values
         noChans = afgrid%noChans
      else if (properties(p_xyz)) then
         ! XYZ quantity (e.g. ECI/ECR stuff)
         qty%frequencyCoordinate = l_xyz
         noChans = 3
      else if (properties(p_matrix3x3)) then
         ! XYZ^2 quantity (e.g. rotation matrix)
         qty%frequencyCoordinate = l_matrix3x3
         noChans = 9
      else if (present(qSignal) .and. .not. properties(p_suppressChannels)) then
         ! This is a channel based quantity
         signalInfo = GetSignal(signal)
         qty%frequencyCoordinate = l_channel
         noChans = size(signalInfo%frequencies)
      end if

      ! Now do the setup for the different families of quantities
      if (qty%minorFrame) then
      if (.not. present(mifGeoLocation) .and. &
          .not. (present(filedatabase) .and. present(chunk))) &
            call MLSMessage(MLSMSG_Error, ModuleName, &
            'Need either mifGeolocation or both filedatabase and ' &
            // 'chunk to create minor frame quantity')
         call ConstructMinorFrameQuantity (filedatabase, chunk, qty%instrumentModule, &
              qty, noChans=noChans, regular=qty%regular, mifGeoLocation=mifGeoLocation)
      else if (properties(p_majorFrame)) then
         if (.not. present(mifGeoLocation)) &
            call MLSMessage(MLSMSG_Error, moduleName, "Missing needed mifGeoLocation")
         if (.not. present(chunk)) &
            call MLSMessage(MLSMSG_Error, moduleName, "Need chunk to create major frame quantity")
         ! Setup a major frame quantity
         call ConstructMajorFrameQuantity (chunk, qty%instrumentModule, qty, &
              noChans, mifGeoLocation)
      else
         ! Setup the quantity template
         if (present(ahgrid)) then
            qty%noInstances = ahgrid%noprofs
            call PointQuantityToHGrid (ahgrid, qty)
         else
            call SetupEmptyHGridForQuantity(qty)
         end if

         if (present(avGrid)) then
            qty%verticalCoordinate = avGrid%verticalCoordinate
            qty%surfs => avGrid%surfs
            qty%noSurfs = avgrid%noSurfs
         else
            call SetupEmptyVGridForQuantity(qty)
         end if

         qty%instanceLen = qty%noSurfs * noChans
     end if

     if (present(qBadValue)) qty%badValue = qBadValue

     if (present(qMolecule)) then
         call GetMoleculeIndex(qMolecule, qty%molecule)
     end if

     qty%quantityType = quantityType
     qty%unit = unitsTable(quantityType)
     qty%noChans = noChans
     qty%sideband = sideband
     qty%signal = signal
     qty%radiometer = radiometer

   end function

  ! ----------------------------------  PointQuantityToHGrid  -----
  subroutine PointQuantityToHGrid ( hGrid, qty )

  ! This routine copies HGrid information into an already defined quantity

    ! Dummy arguments
    type (hGrid_T), intent(in) :: HGRID
    type (QuantityTemplate_T), intent(inout) :: QTY

    ! Executable code
    if ( qty%noInstances/=hGrid%noProfs ) call MLSMessage ( MLSMSG_Error,&
      & ModuleName, "Size of HGrid not compatible with size of quantity" )
    if ( .not. qty%stacked ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Cannot copy hGrids into unstacked quantities")

    qty%phi => hGrid%phi
    qty%geodLat => hGrid%geodLat
    qty%lon => hGrid%lon
    qty%time => hGrid%time
    qty%solarTime => hGrid%solarTime
    qty%solarZenith => hGrid%solarZenith
    qty%losAngle => hGrid%losAngle
    qty%noInstancesLowerOverlap = hGrid%noProfsLowerOverlap
    qty%noInstancesUpperOverlap = hGrid%noProfsUpperOverlap

  end subroutine PointQuantityToHGrid

  ! ---------------------------------- SetupEmptyHGridForQuantity
  subroutine SetupEmptyHGridForQuantity ( qty )
    ! Dummy arguments
    type ( QuantityTemplate_T ), intent(inout) :: QTY
    ! Executable code
    qty%sharedHGrid = .false.
    qty%hGridIndex = 0
    nullify( qty%frequencies, qty%geodLat, qty%lon, qty%time, qty%solarTime, &
      & qty%solarZenith, qty%losAngle ) ! Lest we deallocate a database entry
    call Allocate_test ( qty%phi, 1, 1, 'qty%phi(1,1)', ModuleName )
    call Allocate_test ( qty%geodLat, 1, 1, 'qty%geodLat(1,1)', ModuleName )
    call Allocate_test ( qty%lon, 1, 1, 'qty%lon(1,1)', ModuleName )
    call Allocate_test ( qty%time, 1, 1, 'qty%time(1,1)', ModuleName )
    call Allocate_test ( qty%solarTime, 1, 1, 'qty%solarTime(1,1)', ModuleName )
    call Allocate_test ( qty%solarZenith, 1, 1, 'qty%solarZenith(1,1)', ModuleName )
    call Allocate_test ( qty%losAngle, 1, 1, 'qty%losAngle(1,1)', ModuleName )
    qty%phi = 0.0
    qty%geodLat = 0.0
    qty%lon = 0.0
    qty%time = 0.0
    qty%solarTime = 0.0
    qty%solarZenith = 0.0
    qty%losAngle = 0.0
  end subroutine SetupEmptyHGridForQuantity

  ! ---------------------------------- SetupEmptyVGridForQuantity
  subroutine SetupEmptyVGridForQuantity ( qty )
    ! Dummy arguments
    type ( QuantityTemplate_T ), intent(inout) :: QTY
    ! Executable code
    qty%sharedVGrid = .false.
    qty%vGridIndex = 0
    qty%verticalCoordinate = l_none
    nullify(qty%surfs) ! Lest we deallocate a database entry
    call Allocate_test ( qty%surfs, 1, 1, 'qty%surfs(1,1)', ModuleName )
    qty%surfs = 0. ! We used to have impossible values for bnd. prs.
  end subroutine SetupEmptyVGridForQuantity

  ! -----------------------------------------------  Announce_Error  -----
  subroutine Announce_Error ( where, message, extra, severity )

    use LEXER_CORE, only: PRINT_SOURCE
    use OUTPUT_M, only: BLANKS, OUTPUT
    use TREE, only: SOURCE_REF
    use Intrinsic, only: LIT_INDICES
    use String_Table, only: DISPLAY_STRING
    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR

    integer, intent(in) :: WHERE   ! Tree node where error was noticed
    character (LEN=*), intent(in) :: MESSAGE
    integer, intent(in), optional :: EXTRA
    character(len=*), intent(in), optional :: SEVERITY ! 'nonfatal' or 'fatal'
    character(len=8) :: mySeverity

    mySeverity = 'fatal'
    if ( present(severity) ) then
      if (index('WwNn', severity(1:1)) > 0 ) mySeverity='warning'
    end if
    if ( mySeverity /= 'fatal' ) then
      call blanks(5)
      call output ( ' (warning) ' )
    else
      call blanks(5, fillChar='*')
      call output ( ' (fatal) ' )
    end if
    call output ( 'At ' )
    if ( where > 0 ) then
      call print_source ( source_ref(where) )
    else
      call output ( '(no lcf tree available)' )
    end if
    call output ( ': ' )
    call output ( message )
    if ( present ( extra ) ) call display_string ( lit_indices ( extra ), strip=.true. )
    call output ( '', advance='yes' )
    if ( mySeverity == 'fatal') &
      & call MLSMessage ( MLSMSG_Error, ModuleName, 'Problem in creating quantity' )
  end subroutine Announce_Error

end module
