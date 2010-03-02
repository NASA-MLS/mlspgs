module CFM_QuantityTemplate
   use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test, Test_Allocate
   use FGrid, only: fGrid_T
   use HGridsDatabase, only: hGrid_T
   use VGridsDatabase, only: VGrids, VGrid_T
   use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
   use QuantityTemplates, only: QuantityTemplate_T, NullifyQuantityTemplate
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
      L_MIFDEADTIME, &
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
    use ConstructQuantityTemplates, only : AnyGoodSignalData
    use Parse_Signal_m, only: PARSE_SIGNAL
    use MLSSignals_m, only: GetModuleFromSignal, GetRadiometerFromSignal
 
   implicit none

   public :: CreateQtyTemplate, InitQuantityTemplates
   private

  ! The various properties has/can have
  integer, parameter :: NEXT = -1
  integer, parameter :: P_CHUNKED            = 1
  integer, parameter :: P_MAJORFRAME         = P_CHUNKED + 1
  integer, parameter :: P_MINORFRAME         = P_MAJORFRAME + 1
  integer, parameter :: P_MUSTBEZETA         = P_MINORFRAME + 1
  integer, parameter :: P_FGRID              = P_MUSTBEZETA + 1
  integer, parameter :: P_FGRIDOPTIONAL      = P_FGRID + 1
  integer, parameter :: P_FLEXIBLEVHGRID     = P_FGRIDOPTIONAL + 1
  integer, parameter :: P_HGRID              = P_FLEXIBLEVHGRID + 1
  integer, parameter :: P_MODULE             = P_HGRID + 1
  integer, parameter :: P_MOLECULE           = P_MODULE + 1
  integer, parameter :: P_SGRID              = P_MOLECULE + 1
  integer, parameter :: P_VGRID              = P_SGRID + 1
  integer, parameter :: P_RADIOMETER         = P_VGRID + 1
  integer, parameter :: P_RADIOMETEROPTIONAL = P_RADIOMETER + 1
  integer, parameter :: P_REFLECTOR          = P_RADIOMETEROPTIONAL + 1
  integer, parameter :: P_SCMODULE           = P_REFLECTOR + 1
  integer, parameter :: P_SIGNAL             = P_SCMODULE + 1
  integer, parameter :: P_SIGNALOPTIONAL     = P_SIGNAL + 1
  integer, parameter :: P_SUPPRESSCHANNELS   = P_SIGNALOPTIONAL + 1
  integer, parameter :: P_XYZ                = P_SUPPRESSCHANNELS + 1
  integer, parameter :: P_MATRIX3X3          = P_XYZ + 1

  integer, parameter :: NOPROPERTIES = P_MATRIX3X3


   character(len=20) :: ModuleName = "CFM_QuantityTemplate"
   integer, save, dimension (first_lit : last_lit) :: unitsTable
   logical, dimension(noProperties, first_lit : last_lit), save :: propertyTable

   contains

   type(QuantityTemplate_T) function CreateQtyTemplate (type, filedatabase, chunk, &
        avgrid, ahgrid, qInstModule, qMolecule, qLogBasis, qMinValue, qSignal) &
        & result(qty)
      character(len=*), intent(in) :: type
      type (MLSFile_T), dimension(:), pointer, optional ::     FILEDATABASE
      type (MLSChunk_T), intent(in), optional :: Chunk
      type(VGrid_T), intent(in), optional :: avgrid
      type(HGrid_T), intent(in), optional :: ahgrid
      character(len=*), optional :: qInstModule
      character(len=*), optional :: qMolecule
      character(len=*), optional :: qSignal
      logical, optional :: qLogBasis
      real(r8), optional :: qMinValue

      logical, dimension(noProperties) :: PROPERTIES ! Properties for this quantity type
      character(len=127) :: signalString
      integer, dimension(:), pointer :: SignalInds ! From parse signal
      logical, pointer :: Channels(:)     ! From Parse_Signal
      integer :: s_index

      call NullifyQuantityTemplate(qty)
      qty%instrumentModule = 0
      call GetQtyTypeIndex(type, qty%quantityType)
      qty%name = 0
      qty%fGridIndex = 0
      qty%hGridIndex = 0 ! we're not getting hgrid from hgridDatabase
      qty%vGridIndex = 0 ! we're not getting vgrid from vgridDatabase
      qty%logBasis = .false.
      qty%minValue = -huge(0.0_r8)
      qty%molecule = 0
      qty%noChans = 1
      qty%frequencyCoordinate = l_none
      qty%radiometer = 0
      qty%reflector = 0
      qty%regular = .true.
      qty%sideband = 0
      qty%signal = 0
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
      ! These field are not used in cfm
      qty%noInstancesLowerOverlap = 0
      qty%noInstancesUpperOverlap = 0

      properties = propertyTable(:, qty%quantityType)

      !Sanity check
      if (.not. properties (p_flexibleVHGrid)) then
         if ( present ( ahGrid ) .neqv. properties ( p_hGrid ) ) &
           & call Announce_error ( 0, trim ( merge ( 'unexpected', 'need      ', &
           & present(ahGrid) ) ) // ' hGrid for quantity type ', qty%quantityType, &
           & severity='nonfatal' )
         if ( present ( avGrid ) .neqv. properties ( p_vGrid ) ) &
           & call Announce_error ( 0, trim ( merge ( 'unexpected', 'need      ', &
           & present(avGrid) ) ) // ' vGrid for quantity type ', qty%quantityType )
      end if
      if ( present ( qMolecule ) .neqv. properties ( p_Molecule ) ) &
         & call Announce_error ( 0, trim ( merge ( 'unexpected', 'need      ', &
         & present(qMolecule) ) ) // ' molecule for quantity type ', qty%quantityType )
      if ( present ( qInstModule ) .neqv. ( properties ( p_Module ) .or. properties ( p_scModule) ) ) then
         call Announce_error ( 0, trim ( merge ( 'unexpected', 'need      ', &
         & present(qInstModule) ) ) // ' module for quantity type ', qty%quantityType )
      end if
      if ( .not. properties ( p_signalOptional ) ) then
         if ( present ( qSignal ) .neqv. properties ( p_Signal ) ) &
            & call Announce_error ( 0, trim ( merge ( 'unexpected', 'need      ', &
            & present(qSignal) ) ) // ' signal for quantity type ', qty%quantityType )
      end if

    ! Set up a non major/minor frame quantity
      if (present(avGrid)) then
         qty%verticalCoordinate = avGrid%verticalCoordinate
         qty%surfs => avGrid%surfs
         qty%noSurfs = avgrid%noSurfs
      else
         call SetupEmptyVGridForQuantity(qty)
      end if 

      if (present(ahgrid)) then
         qty%noInstances = ahgrid%noprofs
         call PointQuantityToHGrid (ahgrid, qty)
      else
         call SetupEmptyHGridForQuantity(qty)
      end if

      qty%minorFrame = .not. present(ahgrid)

      qty%unit = unitsTable(qty%quantityType)

      if (present(qInstModule)) then
         call GetModuleIndex (qInstModule, qty%instrumentModule)
      end if

      if (present(qMolecule)) then
         call GetMoleculeIndex(qMolecule, qty%molecule)
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
            & sideband=qty%sideband, channels=channels )
          if ( .not. associated(signalInds) ) then ! A parse error occurred
            call MLSMessage ( MLSMSG_Error, ModuleName,&
              & 'Unable to parse signal string' )
          end if
          if ( size(signalInds) == 1 .or. .not. associated(filedatabase) ) then
            qty%signal = signalInds(1)
          else
            ! Seek a signal with any precision values !< 0
            do s_index=1, size(signalInds)
              if ( AnyGoodSignalData ( signalInds(s_index), qty%sideband, &
                & filedatabase, chunk) ) exit
            end do
            if ( s_index > size(signalInds) ) then
              qty%signal = signalInds(1)
            else
              qty%signal = signalInds(s_index)
            end if
          end if
          call deallocate_test ( signalInds, 'signalInds', ModuleName )
          qty%instrumentModule = GetModuleFromSignal(qty%signal)
          qty%radiometer = GetRadiometerFromSignal(qty%signal)
      end if
   end function

  ! ----------------------------------------------  GetQtyTypeIndex  -----
  subroutine GetQtyTypeIndex(string_text, QtyType)
    use Intrinsic, only: LIT_INDICES
    use MLSStrings, only: LowerCase
    use STRING_TABLE, only: GET_STRING, DISPLAY_STRING
    ! Returns the lit index,  given QtyType name in mixed case
    ! Returns 0 if QtyType name not found
    ! (inverse function: GetQtyTypeName)
    integer, intent(out) :: QtyType
    character (len=*), intent(in) :: string_text
    ! Local variables
    character (len=len(string_text))             :: string_test
    do QtyType=first_lit, last_lit
      ! Skip if not quantity type
      if ( .not. any(PROPERTYTABLE(:, QtyType) ) ) cycle
      call get_string ( lit_indices(QtyType), string_test, strip=.true. )
      if ( LowerCase(trim(string_text)) == LowerCase(trim(string_test))) &
        & return
    end do
    QtyType = 0
  end subroutine GetQtyTypeIndex

  ! ----------------------------------------------- InitQuantityTemplates ----
  subroutine InitQuantityTemplates
    ! This routine initializes the quantity template properties
    ! This is the routine one needs to update when one introduces a new quantity type.

    ! Local variables
    integer :: I                        ! Loop counter
    integer, dimension(0), parameter :: NONE = (/ ( 0, i=1,0 ) /)
    logical :: VALID                    ! Flag
    character(len=132) :: MESSAGE       ! An error message

    ! Executable code
    ! Basically here, we're going to go through and populate the various tables

    propertyTable = .false.
    unitsTable = 0

    call DefineQtyTypes ( (/ &
      l_adopted, phyq_dimensionless, none, next, &
      l_baseline, phyq_temperature, p_flexibleVHGrid, p_fGrid, p_radiometerOptional, &
                  p_signalOptional, p_suppressChannels, p_mustBeZeta, next, &
      l_boundaryPressure, phyq_pressure, p_hGrid, next, &
      l_calSidebandFraction, phyq_dimensionless, p_signal, next, &
      l_chisqBinned, phyq_dimensionless, p_hGrid, p_vGrid, &
                     p_signal, p_suppressChannels, p_mustBeZeta, next, &
      l_chisqChan, phyq_dimensionless, p_majorFrame, p_signal, next, &
      l_chisqMMAF, phyq_dimensionless, p_majorFrame, p_signal, &
                   p_suppressChannels, next, &
      l_chisqMMIF, phyq_dimensionless, p_minorFrame, p_signal, &
                   p_suppressChannels, next, &
      l_cloudExtinction, phyq_dimensionless, p_hGrid, p_vGrid, p_mustBeZeta, next, &
      l_cloudIce, phyq_IceDensity, p_hGrid, p_vGrid, p_mustBeZeta, next, &
      l_cloudInducedRadiance, phyq_temperature, p_minorFrame, p_signal, next, &
      l_cloudMinMax, phyq_temperature, p_hGrid, p_vGrid, p_mustbezeta, &
                     p_signal, p_suppressChannels, next, &
      l_cloudRadSensitivity, phyq_temperature, p_minorFrame, p_signal, next, &
      l_cloudTemperature, phyq_temperature, p_hGrid, p_vGrid, p_mustBeZeta, next, &
      l_cloudWater, phyq_dimensionless, p_hGrid, p_vGrid, p_mustBeZeta, next, &
      l_columnAbundance, phyq_colmabundance, p_hGrid, p_molecule, next, &
      l_dnwt_abandoned, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_ajn, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_axmax, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_cait, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_chisqminnorm, phyq_dimensionless, p_vGrid /) )

    call DefineQtyTypes ( (/ &
      l_dnwt_chisqnorm, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_chisqratio, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_count, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_diag, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_dxdx, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_dxdxl, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_dxn, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_dxnl, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_flag, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_fnmin, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_fnorm, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_gdx, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_gfac, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_gradn, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_sq, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_sqt, phyq_dimensionless, p_vGrid, next, &
      l_geolocation, phyq_dimensionless, p_majorframe, p_module, next, &
      l_earthRadius, phyq_length, p_hGrid, next, &
      l_earthRefl, phyq_dimensionless, none /) )

    call DefineQtyTypes ( (/ &
      l_ecrToFOV, phyq_dimensionless, p_minorFrame, p_module, p_matrix3x3, next, &
      l_effectiveOpticalDepth, phyq_dimensionless, p_minorFrame, p_signal, next, &
      l_elevOffset, phyq_angle, p_signal, next, &
      l_extinction, phyq_extinction, p_hGrid, p_vGrid, p_fGrid, p_radiometer, &
                    p_mustBeZeta, next, &
      l_extinctionv2, phyq_extinction, p_hGrid, p_vGrid, p_fGrid, p_radiometer, &
                      p_mustBeZeta, next, &
      l_fieldAzimuth, phyq_angle, p_hGrid, p_vGrid, next, &
      l_fieldElevation, phyq_angle, p_hGrid, p_vGrid, next, &
      l_fieldStrength, phyq_gauss, p_hGrid, p_vGrid, next, &
      l_geodAltitude, phyq_length, p_minorFrame, p_module, next, &
      l_gph, phyq_length, p_hGrid, p_vGrid, p_mustBeZeta, next, &
      l_heightOffset, phyq_length, p_hGrid, p_vGrid, next, &
      l_isotopeRatio, phyq_dimensionless, p_molecule, next, &
      l_iwc, phyq_icedensity, p_hGrid, p_vGrid, next, &
      l_jacobian_cols, phyq_dimensionless, p_vGrid, next, &
      l_jacobian_rows, phyq_dimensionless, p_vGrid, next /) )

    call DefineQtyTypes ( (/ &
      l_l1bMAFBaseline, phyq_temperature, p_majorFrame, p_signal, next, &
      l_l1bMIF_TAI, phyq_time, p_minorFrame, p_scmodule, next, &
      l_limbSidebandFraction, phyq_dimensionless, p_signal, next, &
      l_lineCenter, phyq_frequency, p_hGrid, p_vGrid, p_molecule, next, &
      l_lineWidth,  phyq_frequency, p_hGrid, p_vGrid, p_molecule, next, &
      l_lineWidth_TDep, phyq_dimensionless, p_hGrid, p_vGrid, p_molecule, next, &
      p_hGrid, p_vGrid, next, &
      l_losTransFunc, phyq_dimensionless, p_minorFrame, p_sGrid, p_module, next, &
      l_losVel, phyq_dimensionless, p_minorFrame, p_module, next, &
      l_magneticField, phyq_gauss, p_vGrid, p_hGrid, p_xyz, p_mustBeZeta, next, &
      l_massMeanDiameterIce, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_massMeanDiameterWater, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_mifDeadTime, phyq_time, next, &
      l_noRadsBinned, phyq_dimensionless, p_vGrid, p_hGrid, &
                      p_signal, p_suppressChannels, p_mustBeZeta, next, &
      l_noRadsPerMIF, phyq_dimensionless, p_minorFrame, p_signal, &
                      p_suppressChannels, next, &
      l_noiseBandwidth, phyq_frequency, p_signal, next, &
      l_numGrad, phyq_dimensionless, p_vGrid, next, &
      l_numJ, phyq_dimensionless, p_vGrid, next, &
      l_numNewt, phyq_dimensionless, p_vGrid, next /) )

    call DefineQtyTypes ( (/ &
      l_opticalDepth, phyq_dimensionless, p_minorFrame, p_signal, next, &
      l_orbitInclination, phyq_angle, p_minorFrame, p_scModule, next, &
      l_phaseTiming, phyq_dimensionless, p_vGrid, next, &
      l_phiTan, phyq_angle, p_minorFrame, p_module, next, &
      l_ptan, phyq_zeta, p_minorFrame, p_module, next, &
      l_quality, phyq_dimensionless, p_hGrid, next, &
      l_radiance, phyq_temperature, p_minorFrame, p_signal, next, &
      l_refltemp, phyq_temperature, p_majorFrame, p_reflector, p_module, next, &
      l_refltrans, phyq_dimensionless, p_signal, p_reflector, next, &
      l_reflrefl, phyq_dimensionless, p_signal, p_reflector, next, &
      l_reflspill, phyq_temperature, p_signal, p_majorframe, p_reflector, next, &
      l_refGPH, phyq_length, p_hGrid, p_vGrid, p_mustBeZeta, next, &
      l_rhi, phyq_dimensionless, p_hGrid, p_vGrid, p_molecule, p_mustBeZeta, next, &
      l_scECI, phyq_length, p_minorFrame, p_scModule, p_xyz, next, &
      l_scGeocAlt, phyq_length, p_minorFrame, p_scModule, next, &
      l_scVel, phyq_velocity, p_minorFrame, p_scModule, p_xyz, next, &
      l_scVelECI, phyq_velocity, p_minorFrame, p_scModule, p_xyz, next, &
      l_scVelECR, phyq_velocity, p_minorFrame, p_scModule, p_xyz, next, &
      l_scanResidual, phyq_length, p_minorFrame, p_module, next, &
      l_scatteringAngle, phyq_angle, p_vGrid, next, &
      l_singleChannelRadiance, phyq_temperature, p_minorFrame, p_signal, &
                               p_suppressChannels, next, &
      l_sizeDistribution, phyq_dimensionless, p_hGrid, p_vGrid, p_mustBeZeta, next, &
      l_spaceRadiance, phyq_temperature, none, next, &
      l_status, phyq_dimensionless, p_hGrid, next, &
      l_strayRadiance, phyq_temperature, p_signal, p_majorFrame, next, &
      l_surfaceHeight, phyq_length, p_hGrid, next, &
      l_surfaceType, phyq_dimensionless, p_hGrid, next, &
      l_systemTemperature, phyq_temperature, p_signal, next, &
      l_temperature, phyq_temperature, p_hGrid, p_vGrid, p_mustbezeta, next, &
      l_totalPowerWeight, phyq_dimensionless, p_signal, next, &
      l_tngtECI, phyq_length, p_minorFrame, p_module, p_xyz, next, &
      l_tngtGeocAlt, phyq_length, p_minorFrame, p_module, next, &
      l_tngtGeodAlt, phyq_length, p_minorFrame, p_module, next, &
      l_TScat, phyq_temperature, p_hGrid, p_signal, p_vGrid, next, &
      l_vmr, phyq_vmr, p_hGrid, p_vGrid, p_fGridOptional, &
             p_molecule, p_radiometerOptional, p_mustbezeta, next /) )

    ! Do a bit of checking
    do i = first_lit, last_lit
      valid = .true.
      message =  ''
      ! Check it's not both major and minor frame
      if ( count ( propertyTable ( (/ p_minorFrame, p_majorFrame /), i ) ) > 1 ) then
        valid = .false.
        message = "Quantity cannot be both major and minor frame"
      end if
      ! Check that we can identify the module for major/minor frame
      if ( ( propertyTable ( p_minorFrame, i ) .or. &
        &    propertyTable ( p_majorFrame, i ) ) .and. &
        & .not. any ( propertyTable ( &
        & (/ p_radiometer, p_module, p_scModule, p_signal /), i ) ) ) then
        valid = .false.
        message = "Badly defined major/minor frame quantity"
      end if
      ! Check that mustBeZeta quantities have a vGrid
      if ( propertyTable ( p_mustBeZeta, i ) .and. .not. &
        & ( propertyTable ( p_vGrid, i ) .or. propertyTable ( p_flexibleVHGrid, i ) ) ) then
        valid = .false.
        message = "Quantity must have vGrid if it must be on log-pressure"
      end if
      ! Print out any message
      if ( .not. valid ) then
        call output ( "Offending quantity: " )
        call display_string ( lit_indices ( i ), strip=.true., advance='yes' )
        call MLSMessage ( MLSMSG_Error, ModuleName, message )
      end if
    end do
  contains
    ! --------------------------- Internal subroutine
    subroutine DefineQtyTypes ( info )
      integer, dimension(:), intent(in) :: INFO
      ! Local variables
      integer :: I                      ! Location
      integer :: QTYTYPE                ! Index
      ! Executable code
      qtyType = 0
      i = 1
      do while ( i <= size(info) )
        if ( qtyType == 0 ) then
          qtyType = info ( i )
          propertyTable ( :, qtyType ) = .false.
          i = i + 1
          if ( i > size(info) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Malformed call to DefineQtyTypes' )
          unitsTable ( qtyType ) = info ( i )
        else
          if ( info(i) /= next ) then
            propertyTable ( info(i), qtyType ) = .true.
          else
            qtyType = 0
          end if
        end if
        i = i + 1
      end do
    end subroutine DefineQtyTypes

  end subroutine InitQuantityTemplates

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
      & call MLSMessage ( MLSMSG_Error, ModuleName, 'Problem in Construct' )
  end subroutine Announce_Error

end module
