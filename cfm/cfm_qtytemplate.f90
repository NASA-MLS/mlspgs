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
      PHYQ_VMR, PHYQ_ZETA, l_ghz
    use Output_M, only: OUTPUT
    use String_Table, only: DISPLAY_STRING
    use MLSCommon, only: r8
    use MLSSignals_m, only: GetModuleIndex
    use Molecules, only: GetMoleculeIndex
    use Chunks_m, only: MLSChunk_T
    use MLSCommon, only: MLSFile_T
    use ConstructQuantityTemplates, only : AnyGoodSignalData, &
            ConstructMajorFrameQuantity, ConstructMinorFrameQuantity, &
            InitQuantityTemplates, unitstable, propertyTable, noProperties, &
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

!---------------------------- RCS Ident Info -------------------------------
   character(len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
   private :: not_used_here
!---------------------------------------------------------------------------

   contains

   ! Creating a quantity based on the optional inputs this subroutine is provided with.
   type(QuantityTemplate_T) function CreateQtyTemplate (quantityType, filedatabase, chunk, &
        avgrid, ahgrid, afgrid, qInstModule, qMolecule, qLogBasis, qMinValue, qSignal, &
        qRadiometer, qBadValue) result(qty)
      ! an integer representing a quantity type, see CFM documentation appendix
      integer, intent(in) :: quantityType
      ! is an array of open files (see CFM_MLSSetup)
      type (MLSFile_T), dimension(:), pointer, optional ::     FILEDATABASE
      ! an input holder, storing the time range of the data to be read
      ! (see CFM_MLSSetup)
      type (MLSChunk_T), intent(in), optional :: Chunk
      ! the z-coordinate samples of the spacecraft's path
      type(VGrid_T), intent(in), optional :: avgrid
      ! the (x,y) coordinate samples of the spacecraft's path
      type(HGrid_T), intent(in), optional :: ahgrid
      ! the frequency in which the data is gathered
      type(FGrid_T), intent(in), optional :: afgrid
      ! instrument module, a case-insensitive string, either "THz" or "GHz"
      character(len=*), optional :: qInstModule
      ! an integer representing one of the molecules listed in the CFM document.
      integer, optional :: qMolecule
      ! if type is "radiance", then a qSignal must be provided
      character(len=*), optional :: qSignal
      ! tells whether the value of avgrid and ahgrid
      ! is on logarithmic scale. The default value is .false..
      logical, optional :: qLogBasis
      ! is used only qLogBasis is .true., in which case
      ! it is the threshold to which greater or equal coordinate values are recorded
      ! as their logarithmic equivalent. The default value is 1.0.
      real(r8), optional :: qMinValue
      character(len=*), optional :: qRadiometer
      ! a number to indicate that a data point shouldn't be used. Default
      ! is -huge(0.0_r8)
      real(r8), optional :: qBadValue

      type (QuantityTemplate_T), dimension(:), pointer :: MifGeolocation
      logical, dimension(noProperties) :: PROPERTIES ! Properties for this quantity type
      character(len=127) :: signalString
      integer, dimension(:), pointer :: SignalInds ! From parse signal
      logical, pointer :: Channels(:)     ! From Parse_Signal
      integer :: s_index, temp
      type(Signal_T) :: signalInfo
      integer :: noChans, sideband, signal, radiometer, instrumentModule
      logical :: isMinorFrame

      ! Executables
      call NullifyQuantityTemplate(qty)
      instrumentModule = 0
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
      ! Because some quantities are provided by MLS, and some are supplied
      ! by GMAO, right now for simplicity, let's not share the grids
      qty%sharedVGrid = .false.
      qty%sharedHGrid = .false.
      qty%sharedFGrid = .false.
      qty%badValue = huge(0.0_r8)   !Default bad value
      ! These field are not used in cfm
      qty%noInstancesLowerOverlap = 0
      qty%noInstancesUpperOverlap = 0

      properties = propertyTable(:, quantityType)

      !Sanity check
      if (.not. properties (p_flexibleVHGrid)) then
         if ( present ( ahGrid ) .neqv. properties ( p_hGrid ) ) &
           & call Announce_error ( trim ( merge ( 'unexpected', 'need      ', &
           & present(ahGrid) ) ) // ' hGrid for quantity type ', quantityType, &
           & severity='nonfatal' )
         if ( present ( avGrid ) .neqv. properties ( p_vGrid ) ) &
           & call Announce_error ( trim ( merge ( 'unexpected', 'need      ', &
           & present(avGrid) ) ) // ' vGrid for quantity type ', quantityType )
         isMinorFrame = properties(p_minorFrame)
      else
         if (present(ahGrid) .neqv. present(avgrid)) &
            call Announce_Error ( 'Must supply both or neither vGrid and hGrid ' &
            // 'for quantity type ', quantityType, severity='fatal' )
         isMinorFrame = .not. present(ahgrid)
      end if
      if ( present ( qMolecule ) .neqv. properties ( p_Molecule ) ) &
         & call Announce_error ( trim ( merge ( 'unexpected', 'need      ', &
         & present(qMolecule) ) ) // ' molecule for quantity type ', quantityType )
      if ( present ( qInstModule ) .neqv. ( properties ( p_Module ) .or. properties ( p_scModule) ) ) then
         call Announce_error ( trim ( merge ( 'unexpected', 'need      ', &
         & present(qInstModule) ) ) // ' module for quantity type ', quantityType )
      end if
      if ( .not. properties ( p_signalOptional ) ) then
         if ( present ( qSignal ) .neqv. properties ( p_Signal ) ) &
            & call Announce_error ( trim ( merge ( 'unexpected', 'need      ', &
            & present(qSignal) ) ) // ' signal for quantity type ', quantityType )
      end if
      if ( properties ( p_mustBeZeta ) .and. present(avGrid) ) then
         if ( avGrid%verticalCoordinate /= l_zeta ) &
           & call Announce_error ( 'Expecting log pressure coordinates for', &
           & quantityType )
      end if
      if (.not. properties(p_fGridOptional)) then
         if ( present ( aFGrid ) .neqv. properties ( p_fGrid ) ) &
           & call Announce_error ( trim ( merge ( 'unexpected', 'need      ', &
           & present(aFGrid) ) ) // ' fGrid for quantity type ', quantityType )
      end if
      if ( .not. properties ( p_radiometerOptional ) ) then
         if ( present ( qRadiometer ) .neqv. properties ( p_Radiometer ) ) &
           & call Announce_error ( trim ( merge ( 'unexpected', 'need      ', &
           & present(qRadiometer) ) ) // ' radiometer for quantity type ', quantityType )
      end if

      if (present(qInstModule)) then
         call GetModuleIndex (qInstModule, instrumentModule)
      end if
      if (present (qRadiometer)) then
         call GetRadiometerIndex(qRadiometer, radiometer)
         ! Every radiometer pairs with only one instrument module
         instrumentModule = GetModuleFromRadiometer(radiometer)
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
          instrumentModule = GetModuleFromSignal(signal)
          radiometer = GetRadiometerFromSignal(signal)
     end if

      ! After all the possible places where instrumentModule can be set
      if ( properties ( p_scModule ) ) then
         if ( .not. IsModuleSpacecraft ( instrumentModule ) ) &
           & call Announce_error ( 'Module must be spacecraft' )
      end if

      ! Now establish the frequency coordinate system
      if (present(afgrid)) then
         qty%frequencyCoordinate = afgrid%frequencyCoordinate
         call allocate_test(qty%frequencies, size(afgrid%values), &
         "qty%frequencies", moduleName)
         qty%frequencies = afgrid%values
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
      if (isMinorFrame) then
         if (.not. (present(filedatabase) .and. present(chunk))) &
            call MLSMessage(MLSMSG_Error, ModuleName, &
            'Need either mifGeolocation or both filedatabase and ' &
            // 'chunk to create minor frame quantity')
         call ConstructMinorFrameQuantity (instrumentModule, &
              qty, noChans=noChans, regular=qty%regular, &
              filedatabase=filedatabase, chunk=chunk)
      else if (properties(p_majorFrame)) then
         if (.not. present(chunk) .or. .not. present(filedatabase)) &
            call MLSMessage(MLSMSG_Error, moduleName, "Need both chunk and filedatabase &
            to create major frame quantity")
         ! We have to construct this because ConstructMajorFrameQuantity
         ! doesn't work with mifGeolocation. The proper fix should be
         ! ConstructMajorFrameQuantity, but I've got no time to test that.
         call ConstructMIFGeolocation (mifGeoLocation, filedatabase, chunk)
         ! Setup a major frame quantity
         call ConstructMajorFrameQuantity (chunk, instrumentModule, qty, &
              noChans, mifGeoLocation)
         call DestroyQuantityTemplateDatabase (mifGeoLocation)
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
            call allocate_test(qty%surfs, size(avGrid%surfs,1), &
            size(avGrid%surfs,2), "qty%surfs", moduleName)
            qty%surfs = avGrid%surfs
            qty%noSurfs = avgrid%noSurfs
         else
            call SetupEmptyVGridForQuantity(qty)
         end if

         qty%instanceLen = qty%noSurfs * noChans
     end if

     if (present(qBadValue)) qty%badValue = qBadValue

     if (present(qMolecule)) then
         qty%molecule = qMolecule
     end if

     qty%quantityType = quantityType
     qty%unit = unitsTable(quantityType)
     qty%noChans = noChans
     qty%sideband = sideband
     qty%signal = signal
     qty%radiometer = radiometer
     qty%instrumentModule = instrumentModule
     qty%minorFrame = isMinorFrame

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

    call allocate_test(qty%phi, size(hGrid%phi,1), size(hGrid%phi,2), &
    "qty%phi", moduleName)
    call allocate_test(qty%geodLat, size(hGrid%geodLat,1), size(hGrid%geodLat,2), &
    "qty%geodLat", moduleName)
    call allocate_test(qty%lon, size(hGrid%lon,1), size(hGrid%lon,2), &
    "qty%lon", moduleName)
    call allocate_test(qty%time, size(hGrid%time,1), size(hGrid%time,2), &
    "qty%time", moduleName)
    call allocate_test(qty%solarTime, size(hGrid%solarTime,1), &
    size(hGrid%solarTime,2), "qty%solarTime", moduleName)
    call allocate_test(qty%solarZenith, size(hGrid%solarZenith,1), &
    size(hGrid%solarZenith,2), "qty%solarZenith", moduleName)
    call allocate_test(qty%losAngle, size(hGrid%losAngle,1), &
    size(hGrid%losAngle,2), "qty%losAngle", moduleName)
    qty%phi = hGrid%phi
    qty%geodLat = hGrid%geodLat
    qty%lon = hGrid%lon
    qty%time = hGrid%time
    qty%solarTime = hGrid%solarTime
    qty%solarZenith = hGrid%solarZenith
    qty%losAngle = hGrid%losAngle
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
  subroutine Announce_Error ( message, extra, severity )

    use LEXER_CORE, only: PRINT_SOURCE
    use OUTPUT_M, only: BLANKS, OUTPUT
    use TREE, only: SOURCE_REF
    use Intrinsic, only: LIT_INDICES
    use String_Table, only: DISPLAY_STRING
    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR

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
    call output ( '(no lcf tree available)' )
    call output ( ': ' )
    call output ( message )
    if ( present ( extra ) ) &
       call display_string ( lit_indices ( extra ), strip=.true. )
    call output ( '', advance='yes' )
    if ( mySeverity == 'fatal') &
      & call MLSMessage ( MLSMSG_Error, ModuleName, 'Problem in creating quantity' )
  end subroutine Announce_Error

  !--------------------------- end bloc --------------------------------------
   logical function not_used_here()
   character (len=*), parameter :: IdParm = &
       "$Id$"
   character (len=len(idParm)) :: Id = idParm
      not_used_here = (id(1:1) == ModuleName(1:1))
      print *, Id ! .mod files sometimes change if PRINT is added
   end function not_used_here
!---------------------------------------------------------------------------

end module
