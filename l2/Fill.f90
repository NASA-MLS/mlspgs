! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module Fill                     ! Create vectors and fill them.
  !=============================================================================

  ! This module performs the Fill operation in the Level 2 software.  
  ! This takes a vector template, and creates and fills an appropriate vector

  implicit none
  private
  public :: MLSL2Fill

! === (start of toc) ===
! MLSL2Fill          given a vector template, and creates and fills a vector
! === (end of toc) ===

! === (start of api) ===
! MLSL2Fill (int root, l1binfo_t l1binfo, *griddedData_T GriddedDataBase(:),
!        *vectorTemplate_T VectorTemplates(:), &
!        *vector_t Vectors(:), *quantityTemplate_T QtyTemplates(:),
!        *matrix_database_T Matrices(:), VGrid_T vGrids(:),
!        *l2GPData_T L2GPDatabase(:), *l2AUXData_T L2AUXDatabase(:),
!        *mlSChunk_T Chunks(:), int ChunkNo )      
! === (end of api) ===
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  !---------------------------------------------------  MLSL2Fill  -----

  subroutine MLSL2Fill ( Root, L1bInfo, GriddedDataBase, VectorTemplates, &
    & Vectors, QtyTemplates , Matrices, vGrids, L2GPDatabase, L2AUXDatabase, &
    & Chunks, ChunkNo )

    ! This is the main routine for the module.  It parses the relevant lines
    ! of the l2cf and works out what to do.

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use Expr_M, only: EXPR
    use GriddedData, only: GriddedData_T
    ! We need many things from Init_Tables_Module.  First the fields:
    use INIT_TABLES_MODULE, only: F_A, F_APRIORIPRECISION, F_B, F_BOUNDARYPRESSURE, &
      & F_COLUMNS, F_DESTINATION, F_DIAGONAL, F_dontMask, F_EARTHRADIUS, &
      & F_EXPLICITVALUES, F_EXTINCTION, &
      & F_FRACTION, F_GEOCALTITUDEQUANTITY, F_GPHQUANTITY, F_HIGHBOUND, F_H2OQUANTITY, &
      & F_H2OPRECISIONQUANTITY, &
      & F_IGNORENEGATIVE, F_IGNOREZERO, F_INSTANCES, F_INTEGRATIONTIME, &
      & F_INTERPOLATE, F_INVERT, F_INTRINSIC, F_ISPRECISION, &
      & F_LENGTHSCALE, F_LOSQTY, F_LOWBOUND, F_LSB, F_LSBFRACTION, &
      & F_MANIPULATION, F_MATRIX, F_MAXITERATIONS, F_MEASUREMENTS, F_METHOD, &
      & F_MODEL, F_MULTIPLIER, F_NOFINEGRID, F_NOISE, F_NOISEBANDWIDTH, &
      & F_OFFSETAMOUNT, F_ORBITINCLINATION, F_PHITAN, &
      & F_PHIWINDOW, F_PRECISION, F_PRECISIONFACTOR, &
      & F_PROFILE, F_PROFILEVALUES, F_PTANQUANTITY, &
      & F_QUANTITY, F_RADIANCEQUANTITY, F_RATIOQUANTITY, &
      & F_REFRACT, F_REFGPHQUANTITY, F_REFGPHPRECISIONQUANTITY, F_RESETSEED, &
      & F_RHIQUANTITY, F_Rows, &
      & F_SCECI, F_SCVEL, F_SCVELECI, F_SCVELECR, F_SEED, F_SKIPMASK, &
      & F_SOURCE, F_SOURCEGRID, F_SOURCEL2AUX, F_SOURCEL2GP, &
      & F_SOURCEQUANTITY, F_SOURCEVGRID, F_SPREAD, F_SUPERDIAGONAL, &
      & F_SYSTEMTEMPERATURE, F_TEMPERATUREQUANTITY, F_TEMPPRECISIONQUANTITY, &
      & F_TEMPLATE, F_TNGTECI, &
      & F_TYPE, F_USB, F_USBFRACTION, F_VECTOR, F_VMRQUANTITY, &
      & FIELD_FIRST, FIELD_LAST
    ! Now the literals:
    use INIT_TABLES_MODULE, only: L_ADDNOISE, L_BOUNDARYPRESSURE, L_CHISQCHAN, &
      & L_CHISQMMAF, L_CHISQMMIF, L_CHOLESKY, L_COLUMNABUNDANCE, &
      & L_ESTIMATEDNOISE, L_EXPLICIT, L_FOLD, L_GEODALTITUDE, &
      & L_GPH, L_GPHPRECISION, L_GRIDDED, L_H2OFROMRHI, &
      & L_HEIGHT, &
      & L_HYDROSTATIC, L_ISOTOPE, L_ISOTOPERATIO, L_KRONECKER, L_L1B, L_L2GP, &
      & L_L2AUX, L_LOSVEL, L_MAGNETICFIELD, L_MAGNETICMODEL, &
      & L_MANIPULATE, L_NEGATIVEPRECISION, L_NOISEBANDWIDTH, L_NONE, &
      & L_OFFSETRADIANCE, L_ORBITINCLINATION, L_PHITAN, &
      & L_PLAIN, L_PRESSURE, L_PROFILE, L_PTAN, &
      & L_RADIANCE, L_RECTANGLEFROMLOS, L_REFGPH, L_REFRACT, L_RHI, &
      & L_RHIFROMH2O, L_RHIPRECISIONFROMH2O, &
      & L_SCECI, L_SCGEOCALT, L_SCVEL, L_SCVELECI, L_SCVELECR, &
      & L_SIDEBANDRATIO, L_SPD, L_SPECIAL, L_SYSTEMTEMPERATURE, &
      & L_TEMPERATURE, L_TNGTECI, L_TNGTGEODALT, &
      & L_TNGTGEOCALT, L_TRUE, L_VECTOR, L_VGRID, L_VMR, L_XYZ, L_ZETA
    ! Now the specifications:
    use INIT_TABLES_MODULE, only: S_DESTROY, S_DUMP, S_FILL, S_FILLCOVARIANCE, &
      & S_FILLDIAGONAL, S_MATRIX,  S_NEGATIVEPRECISION, S_SNOOP, S_TIME, &
      & S_TRANSFER, S_VECTOR
    ! Now some arrays
    use Intrinsic, only: FIELD_INDICES
    use Intrinsic, only: &
      & PHYQ_Dimensionless, PHYQ_Invalid, PHYQ_Temperature, &
      & PHYQ_Time, PHYQ_Length, PHYQ_Pressure, PHYQ_Zeta, PHYQ_Angle, PHYQ_Profiles
    use L1BData, only: DeallocateL1BData, Dump, FindL1BData, L1BData_T, &
      & PRECISIONSUFFIX, ReadL1BData, AssembleL1BQtyName
    use L2GPData, only: L2GPData_T
    use L2AUXData, only: L2AUXData_T
    use L3ASCII, only: L3ASCII_INTERP_FIELD
    use LEXER_CORE, only: PRINT_SOURCE
    use ManipulateVectorQuantities, only: DOFGRIDSMATCH, DOHGRIDSMATCH, &
      & DOVGRIDSMATCH, DOQTYSDESCRIBESAMETHING
    use MatrixModule_0, only: Sparsify, MatrixInversion
    use MatrixModule_1, only: AddToMatrixDatabase, CreateEmptyMatrix, &
      & DestroyMatrix, Dump, GetActualMatrixFromDatabase, GetDiagonal, &
      & FindBlock, GetKindFromMatrixDatabase, GetFromMatrixDatabase, K_SPD, &
      & Matrix_Cholesky_T, Matrix_Database_T, Matrix_Kronecker_T, Matrix_SPD_T, &
      & Matrix_T, NullifyMatrix, UpdateDiagonal
    ! NOTE: If you ever want to include defined assignment for matrices, please
    ! carefully check out the code around the call to snoop.
    use MLSCommon, only: FileNameLen, L1BInfo_T, MLSChunk_T, R8, RM, RV
    use MLSFiles, only: mls_hdf_version, HDFVERSION_5, &
      & ERRORINH5FFUNCTION, WRONGHDFVERSION
    use MLSL2Options, only: LEVEL1_HDFVERSION
    use MLSL2Timings, only: SECTION_TIMES, TOTAL_TIMES
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Allocate, MLSMSG_Deallocate
    use MLSNumerics, only: InterpolateValues
    use MLSRandomNumber, only: drang, mls_random_seed, MATH77_RAN_PACK
    use MLSSignals_m, only: GetSignalName, GetModuleName, IsModuleSpacecraft
    use Molecules, only: L_H2O
    use MoreTree, only: Get_Boolean, Get_Field_ID, Get_Spec_ID, GetIndexFlagsFromList
    use OUTPUT_M, only: BLANKS, OUTPUT
    use QuantityTemplates, only: Epoch, QuantityTemplate_T
    use RHIFromH2O, only: RHIFromH2O_Factor, RHIPrecFromH2O
    use ScanModelModule, only: GetBasisGPH, Get2DHydrostaticTangentPressure, GetGPHPrecision
    use SnoopMLSL2, only: SNOOP
    use String_Table, only: Display_String
    use Time_M, only: Time_Now
    use TOGGLES, only: GEN, LEVELS, SWITCHES, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, &
      & SOURCE_REF, SUB_ROSA, SUBTREE
    use TREE_TYPES, only: N_NAMED, N_SET_ONE
    use UNITS, only: PI
    use VectorsModule, only: AddVectorToDatabase, &
      & ClearUnderMask, CopyVector, CreateMask, CreateVector, &
      & DestroyVectorInfo, Dump, &
      & GetVectorQtyByTemplateIndex, isVectorQtyMasked, MaskVectorQty, &
      & rmVectorFromDatabase, ValidateVectorQuantity, Vector_T, &
      & VectorTemplate_T, VectorValue_T, M_Fill, M_LinAlg
    use VGridsDatabase, only: VGRID_T, GETUNITFORVERTICALCOORDINATE

    ! Dummy arguments
    integer, intent(in) :: ROOT    ! Of the FILL section in the AST
    type (l1BInfo_T), intent(in) :: L1bInfo
    type (griddedData_T), dimension(:), pointer :: GriddedDataBase
    type (vectorTemplate_T), dimension(:), pointer :: VectorTemplates
    type (vector_T), dimension(:), pointer :: Vectors
    type (quantityTemplate_T), dimension(:), pointer :: QtyTemplates
    type (matrix_database_T), dimension(:), pointer :: Matrices
    type (VGrid_T), dimension(:), intent(in) :: vGrids
    type (l2GPData_T), dimension(:), pointer :: L2GPDatabase
    type (l2AUXData_T), dimension(:), pointer :: L2AUXDatabase
    type (mlSChunk_T), dimension(:), pointer :: Chunks
    integer, intent(in) :: ChunkNo

    ! -----     Declarations for Fill and internal subroutines     -------

    integer :: ERROR
    logical, parameter :: DEEBUG = .FALSE.                 ! Usually FALSE

    ! Error codes for "announce_error"  
    integer, parameter :: No_Error_code = 0
    integer, parameter :: Wrong_Number = No_Error_code+1     ! of fields of a VECTOR command
    integer, parameter :: UnknownQuantityName = wrong_number + 1
    integer, parameter :: Source_not_in_db = unknownQuantityName + 1
    integer, parameter :: ZeroProfilesFound = source_not_in_db + 1
    integer, parameter :: ZeroGeodSpan = zeroProfilesFound + 1
    integer, parameter :: CantFillFromL2AUX = ZeroGeodSpan + 1
    integer, parameter :: VectorWontMatchPDef = cantFillFromL2AUX + 1
    integer, parameter :: CantFillFromL1B = vectorWontMatchPDef + 1

    ! Error codes for "Matrix" specification
    integer, parameter :: MissingField = cantFillFromL1B + 1

    ! More Error codes relating to FillVector
    integer, parameter :: NumInstancesisZero = missingField + 1
    integer, parameter :: NumSurfsisZero = numInstancesisZero + 1
    integer, parameter :: NumChansisZero = numSurfsisZero + 1
    integer, parameter :: ObjIsFullRank3 = numChansisZero + 1
    integer, parameter :: OtherErrorInFillVector = objIsFullRank3 + 1
    integer, parameter :: NoSourceGridGiven= otherErrorInFillVector + 1
    integer, parameter :: NoSourceL2GPGiven= noSourceGridGiven + 1
    integer, parameter :: NoSourceL2AUXGiven= noSourceL2GPGiven + 1
    integer, parameter :: NoExplicitValuesGiven= noSourceL2AUXGiven + 1
    integer, parameter :: NoSourceQuantityGiven= noExplicitValuesGiven + 1
    integer, parameter :: InvalidExplicitFill= noSourceQuantityGiven + 1
    integer, parameter :: BadUnitsForExplicit= invalidExplicitFill + 1
    integer, parameter :: BadUnitsForIntegrationTime = badUnitsForExplicit + 1
    integer, parameter :: BadUnitsForSystemTemperature = badUnitsForIntegrationTime + 1
    integer, parameter :: BadIsotopeFill = badUnitsForSystemTemperature + 1
    integer, parameter :: BadlosGridFill = badIsotopeFill + 1
    integer, parameter :: CantInterpolate3d = badlosGridFill + 1

    ! Error codes resulting from FillCovariance
    integer, parameter :: NotSPD = CantInterpolate3D + 1
    integer, parameter :: NotImplemented = notSPD + 1
    integer, parameter :: BothFractionAndLength = NotImplemented + 1

    ! Error codes resulting from squeeze
    integer, parameter :: N1_is_zero = BothFractionAndLength + 1
    integer, parameter :: N2_is_zero = n1_is_zero + 1
    integer, parameter :: N3_is_zero = n2_is_zero + 1
    integer, parameter :: M1_too_small = n3_is_zero + 1
    integer, parameter :: M2_too_small = m1_too_small + 1
    integer, parameter :: Not_permutation = m2_too_small + 1
    integer, parameter :: Allocation_err = not_permutation + 1
    integer, parameter :: Deallocation_err = allocation_err + 1

    ! Miscellaneous
    integer, parameter :: Miscellaneous_err = deallocation_err + 1
    integer, parameter :: ErrorReadingL1B = miscellaneous_err + 1
    integer, parameter :: NeedTempREFGPH = errorReadingL1B + 1
    integer, parameter :: NeedH2O = needTempRefGPH + 1
    integer, parameter :: NeedOrbitInclination = needH2O + 1
    integer, parameter :: NeedGeocAltitude = needOrbitInclination + 1
    integer, parameter :: NeedGeodAltitude = needGeocAltitude + 1
    integer, parameter :: BadGeocAltitudeQuantity = needGeodAltitude + 1
    integer, parameter :: BadTemperatureQuantity = badGeocAltitudeQuantity + 1
    integer, parameter :: BadREFGPHQuantity = badTemperatureQuantity + 1
    integer, parameter :: NonConformingHydrostatic = badREFGPHQuantity + 1
    integer, parameter :: BadUnitsForMaxIterations = nonConformingHydrostatic + 1
    integer, parameter :: NoSpecialFill = badUnitsForMaxIterations + 1
    integer, parameter :: BadlosVelFill = noSpecialFill + 1
    integer, parameter :: NotZetaForGrid = BadLosVelFill + 1
    integer, parameter :: BadEstNoiseFill = NotZetaForGrid + 1
    integer, parameter :: BadRefractFill = BadEstNoiseFill + 1

    real, parameter ::    UNDEFINED_VALUE = -999.99 ! Same as %template%badvalue

    ! Local variables

    type (vectorValue_T), pointer :: APRIORIPRECISION
    type (vectorValue_T), pointer :: AQUANTITY
    type (vectorValue_T), pointer :: BNDPRESSQTY
    type (vectorValue_T), pointer :: BQUANTITY
    type (vectorValue_T), pointer :: EARTHRADIUSQTY
    type (vectorValue_T), pointer :: GEOCALTITUDEQUANTITY
    type (vectorValue_T), pointer :: GPHQUANTITY
    type (vectorValue_T), pointer :: H2OPRECISIONQUANTITY
    type (vectorValue_T), pointer :: H2OQUANTITY
    type (vectorValue_T), pointer :: LOSQTY
    type (vectorValue_T), pointer :: LSB
    type (vectorValue_T), pointer :: LSBFRACTION
    type (vectorValue_T), pointer :: MEASQTY
    type (vectorValue_T), pointer :: MODELQTY
    type (vectorValue_T), pointer :: NBWQUANTITY
    type (vectorValue_T), pointer :: NOISEQTY
    type (vectorValue_T), pointer :: ORBITINCLINATIONQUANTITY
    type (vectorValue_T), pointer :: PRECISIONQUANTITY
    type (vectorValue_T), pointer :: QUANTITY ! Quantity to be filled
    type (vectorValue_T), pointer :: RADIANCEQUANTITY
    type (vectorValue_T), pointer :: RATIOQUANTITY
    type (vectorValue_T), pointer :: REFGPHQUANTITY
    type (vectorValue_T), pointer :: REFGPHPRECISIONQUANTITY
    type (vectorValue_T), pointer :: SCECIQUANTITY
    type (vectorValue_T), pointer :: SCVELQUANTITY
    type (vectorValue_T), pointer :: SOURCEQUANTITY
    type (vectorValue_T), pointer :: SYSTEMPQUANTITY
    type (vectorValue_T), pointer :: TEMPERATUREQUANTITY
    type (vectorValue_T), pointer :: TEMPPRECISIONQUANTITY
    type (vectorValue_T), pointer :: TNGTECIQUANTITY
    type (vectorValue_T), pointer :: PHITANQUANTITY
    type (vectorValue_T), pointer :: PTANQUANTITY
    type (vectorValue_T), pointer :: USB
    type (vectorValue_T), pointer :: USBFRACTION
    type (vectorValue_T), pointer :: VMRQTY

    type (Matrix_T), dimension(:), pointer :: SNOOPMATRICES
    type (Matrix_T), pointer :: ONEMATRIX

    integer :: APRPRECQTYINDEX          ! Index of apriori precision quantity    
    integer :: APRPRECVCTRINDEX         ! Index of apriori precision vector
    integer :: AQTYINDEX                ! Index of a quantity in vector
    integer :: AVECINDEX                ! Index of a vector
    integer :: BNDPRESSQTYINDEX
    integer :: BNDPRESSVCTRINDEX
    integer :: BQTYINDEX                ! Index of a quantity in vector
    integer :: BVECINDEX                ! Index of a vector
    integer :: COLVECTOR                ! Vector defining columns of Matrix
    type(matrix_SPD_T), pointer :: Covariance
    integer :: DESTINATIONVECTORINDEX   ! For transfer commands
    !                                     -- for FillCovariance
    integer :: EARTHRADIUSQTYINDEX
    integer :: EARTHRADIUSVECTORINDEX
    integer :: Diagonal                 ! Index of diagonal vector in database
    !                                     -- for FillCovariance
    logical :: DONTMASK                 ! Use even masked values if TRUE
    integer :: ERRORCODE                ! 0 unless error; returned by called routines
    logical :: Extinction               ! Flag for cloud extinction calculation
    integer :: FIELDINDEX               ! Entry in tree
    integer :: FieldValue               ! Value of a field in the L2CF
    integer :: FILLMETHOD               ! How will we fill this quantity
    integer :: FRACTION                 ! Index of fraction vector in database
    integer :: GEOCALTITUDEQUANTITYINDEX    ! In the source vector
    integer :: GEOCALTITUDEVECTORINDEX      ! In the vector database
    integer :: GPHQUANTITYINDEX         ! In the source vector
    integer :: GPHVECTORINDEX           ! In the vector database
    logical, dimension(field_first:field_last) :: GOT
    integer :: GRIDINDEX                ! Index of requested grid
    integer :: GSON                     ! Descendant of Son
    logical :: HIGHBOUND                ! Flag
    integer :: H2OQUANTITYINDEX         ! in the quantities database
    integer :: H2OVECTORINDEX           ! In the vector database
    integer :: H2OPRECISIONQUANTITYINDEX         ! in the quantities database
    integer :: H2OPRECISIONVECTORINDEX           ! In the vector database
    integer :: I, J                     ! Loop indices for section, spec, expr
    integer :: GLOBALUNIT               ! To go into the vector
    logical :: IGNOREZERO               ! Don't sum chi^2 at values of noise = 0
    logical :: IGNORENEGATIVE           ! Don't sum chi^2 at values of noise < 0
    real(r8) :: INTEGRATIONTIME         ! For estimated noise
    logical :: INTERPOLATE              ! Flag for l2gp etc. fill
    integer :: INSTANCE                 ! Loop counter
    integer :: INSTANCESNODE            ! Tree node
    logical :: INVERT                   ! "Invert the specified covariance matrix"
    logical :: ISPRECISION              ! l1b precision, not radiances if TRUE
    integer :: KEY                      ! Definitely n_named
    integer :: L2AUXINDEX               ! Index into L2AUXDatabase
    integer :: L2GPINDEX                ! Index into L2GPDatabase
    integer :: LENGTHSCALE              ! Index of lengthscale vector in database
    integer :: LOSVECTORINDEX           ! index in vector database
    integer :: LOSQTYINDEX              ! index in QUANTITY database
    logical :: LOWBOUND                 ! Flag
    integer :: LSBVECTORINDEX           ! Index in vector database
    integer :: LSBQUANTITYINDEX         ! Index in vector database
    integer :: LSBFRACTIONVECTORINDEX   ! Index in vector database
    integer :: LSBFRACTIONQUANTITYINDEX ! Index in vector database
    type(matrix_Cholesky_T) :: MatrixCholesky
    type(matrix_Kronecker_T) :: MatrixKronecker
    type(matrix_SPD_T) :: MatrixSPD
    type(matrix_T) :: MatrixPlain
    integer :: MANIPULATION             ! String index
    integer :: MATRIXTOFILL             ! Index in database
    integer :: MATRIXTOKILL             ! Index in database
    integer :: MATRIXTYPE               ! Type of matrix, L_... from init_tables
    integer :: MAXITERATIONS            ! For hydrostatic fill
    real, dimension(2) :: MULTIPLIER   ! To scale source,noise part of addNoise
    integer :: MULTIPLIERNODE           ! For the parser
    integer :: NBWVECTORINDEX           ! In vector database
    integer :: NBWQUANTITYINDEX         ! In vector database
    integer :: NOFINEGRID               ! no of fine grids for cloud extinction calculation
    integer :: NOSNOOPEDMATRICES        ! No matrices to snoop
    real(rv) :: OFFSETAMOUNT            ! For offsetRadiance
    integer :: ORBITINCLINATIONVECTORINDEX ! In the vector database
    integer :: ORBITINCLINATIONQUANTITYINDEX ! In the quantities database
    integer :: PHITANVECTORINDEX        ! In the vector database
    integer :: PHITANQUANTITYINDEX      ! In the quantities database
    real(r8) :: PHIWINDOW               ! For hydrostatic ptan guesser
    integer :: PHIWINDOWUNITS           ! For hydrostatic ptan guesser
    real(r8) :: PRECISIONFACTOR         ! For setting -ve error bars
    integer :: PRECISIONQUANTITYINDEX   ! For precision quantity
    integer :: PRECISIONVECTORINDEX     ! In the vector database
    integer :: PROFILE                  ! A single profile
    integer :: PTANVECTORINDEX          ! In the vector database
    integer :: PTANQUANTITYINDEX        ! In the quantities database
    integer :: QUANTITYINDEX            ! Within the vector
    integer :: RADIANCEQUANTITYINDEX    ! For radiance quantity
    integer :: RADIANCEVECTORINDEX      ! For radiance quantity
    integer :: RATIOQUANTITYINDEX       ! in the quantities database
    integer :: RATIOVECTORINDEX         ! In the vector database
    logical :: REFRACT                  ! Do refraction in phiTan fill
    integer :: REFGPHQUANTITYINDEX      ! in the quantities database
    integer :: REFGPHVECTORINDEX        ! In the vector database
    integer :: REFGPHPRECISIONQUANTITYINDEX      ! in the quantities database
    integer :: REFGPHPRECISIONVECTORINDEX        ! In the vector database
    logical :: RESETSEED                ! Let mls_random_seed choose new seed
    integer :: ROWVECTOR                ! Vector defining rows of Matrix
    integer :: SCECIQUANTITYINDEX       ! In the quantities database
    integer :: SCECIVECTORINDEX         ! In the vector database
    integer :: SCVELQUANTITYINDEX       ! In the quantities database
    integer :: SCVELVECTORINDEX         ! In the vector database
    integer, dimension(2) :: SEED       ! integers used by random_numbers
    integer :: SEEDNODE                 ! For the parser
    logical :: SKIPMASK                 ! Flag for transfer
    integer :: SON                      ! Of root, an n_spec_args or a n_named
    integer :: SOURCEQUANTITYINDEX      ! in the quantities database
    integer :: SOURCEVECTORINDEX        ! In the vector database
    logical :: SPREADFLAG               ! Do we spread values accross instances in explict
    integer :: STATUS                   ! Flag from allocate etc.
    integer :: SUPERDIAGONAL            ! Index of superdiagonal matrix in database
    logical :: Switch2intrinsic         ! Have mls_random_seed call intrinsic
    !                                     -- for FillCovariance
    real :: T1, T2                      ! for timing 
    integer :: SYSTEMPQUANTITYINDEX     ! in the quantities database
    integer :: SYSTEMPVECTORINDEX       ! in the vector database
    integer :: TEMPERATUREQUANTITYINDEX ! in the quantities database
    integer :: TEMPERATUREVECTORINDEX   ! In the vector database
    integer :: TEMPPRECISIONQUANTITYINDEX ! in the quantities database
    integer :: TEMPPRECISIONVECTORINDEX   ! In the vector database
    integer :: TEMPLATEINDEX            ! In the template database
    logical :: TIMING
    integer :: TNGTECIQUANTITYINDEX     ! In the quantities database
    integer :: TNGTECIVECTORINDEX       ! In the vector database
    integer, dimension(2) :: UNITASARRAY ! From expr
    integer :: USBVECTORINDEX           ! Inddex in vector database
    integer :: USBQUANTITYINDEX         ! Inddex in vector database
    integer :: USBFRACTIONVECTORINDEX   ! Index in vector database
    integer :: USBFRACTIONQUANTITYINDEX ! Index in vector database
    real(r8), dimension(2) :: VALUEASARRAY ! From expr
    integer :: VALUESNODE               ! For the parser
    integer :: VECTORINDEX              ! In the vector database
    integer :: VECTORNAME               ! Name of vector to create
    integer :: VGRIDINDEX               ! Index of sourceVGrid
    integer :: vmrQtyIndex
    integer :: vmrQtyVctrIndex
    integer :: measQtyIndex
    integer :: measVectorIndex
    integer :: modelQtyIndex
    integer :: modelVectorIndex
    integer :: noiseQtyIndex
    integer :: noiseVectorIndex
    logical :: old_math77_ran_pack      ! To restore math77_ran_pack

    ! Executable code
    timing = section_times
    if ( timing ) call time_now ( t1 )
    old_math77_ran_pack = math77_ran_pack

    if ( toggle(gen) ) call trace_begin ( "MLSL2Fill", root )

    error = 0
    templateIndex = -1
    vectorIndex = -1
    maxIterations = 4

    ! Loop over the lines in the configuration file

    do i = 2, nsons(root)-1 ! Skip the section name at begin and end
      son = subtree(i,root)
      if ( node_id(son) == n_named ) then ! Is spec labeled?
        key = subtree(2,son)
        vectorName = sub_rosa(subtree(1,son))
      else
        key = son
        vectorName = 0
      end if
      dontMask = .false.
      extinction = .false.
      got= .false.
      ignoreZero = .false.
      ignoreNegative = .false.
      interpolate = .false.
      invert = .false.
      isPrecision = .false.
      resetSeed = .false.
      refract = .false.
      spreadFlag = .false.
      switch2intrinsic = .false.
      seed = 0
      noFineGrid = 1
      precisionFactor = 0.5
      offsetAmount = 1000.0             ! Default to 1000K
      profile = -1
      phiWindow = 4
      phiWindowUnits = phyq_angle

      ! Node_id(key) is now n_spec_args.

      select case( get_spec_id(key) )
      case ( s_vector ) ! ===============================  Vector  =====
        got = .false.
        globalUnit = PHYQ_Invalid
        lowBound = .false.
        highBound = .false.
        do j = 2, nsons(key)
          son = subtree(j,key)              ! The field
          fieldIndex = get_field_id(son)
          got(fieldIndex) = .true.
          if ( nsons(son) > 1 ) then
            fieldValue = decoration(subtree(2,son)) ! The field's value
          else
            fieldValue = son
          end if
          select case ( fieldIndex )
          case ( f_template )
            templateIndex = decoration(decoration(subtree(2,son)))
          case ( f_lengthScale )
            if ( get_boolean(fieldValue) ) globalUnit = phyq_length
          case ( f_fraction )
            if ( get_boolean(fieldValue) ) globalUnit = phyq_dimensionless
          case ( f_lowBound )
            lowBound = get_boolean ( fieldValue )
          case ( f_highBound )
            highBound = get_boolean ( fieldValue )
          end select
        end do

        if ( all(got((/f_lengthScale,f_fraction/))) ) &
          & call Announce_Error ( key, bothFractionAndLength )
        ! Create the vector, and add it to the database.
        call decorate ( key, AddVectorToDatabase ( vectors, &
          & CreateVector ( vectorName, vectorTemplates(templateIndex), &
          & qtyTemplates, globalUnit=globalUnit, highBound=highBound, lowBound=lowBound ) ) )

        ! That's the end of the create operation

      case ( s_dump ) ! ============================== Dump ==========
        ! Currently the only field here is the quantity, but let's
        ! do a loop anyway, as I may add more later
        do j = 2, nsons(key)
          gson = subtree(j,key) ! The argument
          fieldIndex = get_field_id(gson)
          if (nsons(gson) > 1) gson = subtree(2,gson) ! Now value of said argument
          got(fieldIndex) = .true.
          select case ( fieldIndex )
          case ( f_vector )
            vectorIndex = decoration(decoration(gson))
          end select
        end do
        call dump ( vectors(vectorIndex) )

      case ( s_matrix ) ! ===============================  Matrix  =====
        got = .false.
        matrixType = l_plain
        do j = 2, nsons(key)
          gson = subtree(j,key)              ! The field
          fieldIndex = get_field_id(gson)
          got(fieldIndex) = .true.
          if ( nsons(gson) > 1 ) &
            & fieldValue = decoration(subtree(2,gson)) ! The field's value
          select case ( fieldIndex )
          case ( f_columns )
            colVector = decoration(fieldValue)
          case ( f_rows )
            rowVector = decoration(fieldValue)
          case ( f_type )
            matrixType = fieldValue
          end select
        end do
        if ( got(f_columns) .and. (got(f_rows) .or. matrixType == l_spd) ) then
          select case ( matrixType )
          case ( l_cholesky )
            call NullifyMatrix ( matrixCholesky%m )
            call createEmptyMatrix ( matrixCholesky%m, vectorName, &
              & vectors(rowVector), vectors(colVector) )
            call decorate ( key, addToMatrixDatabase(matrices, matrixCholesky) )
          case ( l_kronecker )
            call NullifyMatrix ( matrixKronecker%m )
            call createEmptyMatrix ( matrixKronecker%m, vectorName, &
              & vectors(rowVector), vectors(colVector) )
            call decorate ( key, addToMatrixDatabase(matrices, matrixKronecker) )
          case ( l_plain )
            call NullifyMatrix ( matrixPlain )
            call createEmptyMatrix ( matrixPlain, vectorName, vectors(rowVector), &
              vectors(colVector) )
            call decorate ( key, addToMatrixDatabase(matrices, matrixPlain) )
          case ( l_spd )
            call NullifyMatrix ( matrixSPD%m )
            call createEmptyMatrix ( matrixSPD%m, vectorName, &
              & vectors(colVector), vectors(colVector) )
            call decorate ( key, addToMatrixDatabase(matrices, matrixSPD) )
          end select
        else
          call announce_error ( key, missingField, &
            & extraInfo = (/ f_columns, f_rows /) )
        end if

      case ( s_fill ) ! ===================================  Fill  =====
        ! Now we're on actual Fill instructions.
        ! Loop over the instructions to the Fill command

        do j = 2, nsons(key)
          gson = subtree(j,key) ! The argument
          fieldIndex = get_field_id(gson)
          if (nsons(gson) > 1) gson = subtree(2,gson) ! Now value of said argument
          got(fieldIndex)=.TRUE.
          select case ( fieldIndex )
          case ( f_a )
            aVecIndex = decoration(decoration(subtree(1,gson)))
            aQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_aprioriPrecision )
            aprPrecVctrIndex = decoration(decoration(subtree(1,gson)))
            aprPrecQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_b )
            bVecIndex = decoration(decoration(subtree(1,gson)))
            bQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_boundaryPressure )     ! For special fill of columnAbundance
            bndPressVctrIndex = decoration(decoration(subtree(1,gson)))
            bndPressQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_earthRadius ) ! For losGrid fill
            earthRadiusVectorIndex = decoration(decoration(subtree(1,gson)))
            earthRadiusQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_noise )   ! Only used for chi^2 special fills or addnoise
            noiseVectorIndex = decoration(decoration(subtree(1,gson)))
            noiseQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_explicitValues ) ! For explicit fill
            valuesNode=subtree(j,key)
          case ( f_extinction ) ! For cloud extinction fill
            if ( node_id(gson) == n_set_one ) then
              extinction=.TRUE.
            else
              extinction = decoration(subtree(2,gson)) == l_true
            end if
          case ( f_geocAltitudeQuantity ) ! For hydrostatic
            geocAltitudeVectorIndex = decoration(decoration(subtree(1,gson)))
            geocAltitudeQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_gphQuantity ) ! For magnetic field fill
            gphVectorIndex = decoration(decoration(subtree(1,gson)))
            gphQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_h2oQuantity ) ! For hydrostatic or rhi
            h2oVectorIndex = decoration(decoration(subtree(1,gson)))
            h2oQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_h2oPrecisionQuantity ) ! For rhi precision
            h2oPrecisionVectorIndex = decoration(decoration(subtree(1,gson)))
            h2oPrecisionQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_dontMask )
            dontMask = get_boolean ( gson )
          case ( f_ignoreZero )
            ignoreZero = get_boolean ( gson )
          case ( f_ignoreNegative )
            ignoreNegative = get_boolean ( gson )
          case ( f_instances )
            instancesNode = subtree(j,key)
          case ( f_integrationTime )
            call expr ( gson , unitAsArray, valueAsArray )
            if ( all (unitAsArray /= (/PHYQ_Time, PHYQ_Invalid/) ) ) &
              call Announce_error ( key, badUnitsForIntegrationtime )
            integrationTime = valueAsArray(1)
          case ( f_interpolate ) ! For l2gp etc. fill
            if ( node_id(gson) == n_set_one ) then
              interpolate=.TRUE.
            else
              interpolate = decoration(subtree(2,gson)) == l_true
            end if
          case ( f_intrinsic )
            switch2intrinsic = get_boolean ( gson )
!         case ( f_invert )
!           invert = get_boolean ( gson )
          case ( f_isPrecision )
            isPrecision = get_boolean ( gson )
          case ( f_losQty ) ! For losGrid fill
            losVectorIndex = decoration(decoration(subtree(1,gson)))
            losQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_lsb ) ! For folding
            lsbVectorIndex = decoration(decoration(subtree(1,gson)))
            lsbQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_lsbFraction ) ! For folding
            lsbFractionVectorIndex = decoration(decoration(subtree(1,gson)))
            lsbFractionQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_manipulation )
            manipulation = sub_rosa ( gson )
          case ( f_maxIterations )      ! For hydrostatic fill
            call expr ( subtree(2,subtree(j,key)), unitAsArray,valueAsArray )
            if ( all(unitAsArray(1) /= (/PHYQ_Dimensionless,PHYQ_Invalid/)) ) &
              & call Announce_error ( key, badUnitsForMaxIterations )
            maxIterations = valueAsArray(1)
          case ( f_measurements )   ! Only used for diagnostic special fills
            measVectorIndex = decoration(decoration(subtree(1,gson)))
            measQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_method )   ! How are we going to fill it?
            fillMethod = decoration(gson)
          case ( f_model )   ! Only used for diagnostic special fills
            modelVectorIndex = decoration(decoration(subtree(1,gson)))
            modelQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_multiplier ) ! For scaling noise part of addnoise
            multiplierNode=subtree(j,key)
          case ( f_noFineGrid )      ! For cloud extinction fill
            call expr ( subtree(2,subtree(j,key)), unitAsArray,valueAsArray )
            if ( all(unitAsArray(1) /= (/PHYQ_Dimensionless,PHYQ_Invalid/)) ) &
              & call Announce_error ( key, badUnitsForMaxIterations )
            noFineGrid = valueAsArray(1)
          case ( f_noiseBandwidth )
            nbwVectorIndex = decoration(decoration(subtree(1,gson)))
            nbwQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_offsetAmount )    ! For marking unused radiances
            call expr ( subtree(2,subtree(j,key)), unitAsArray, valueAsArray )
            if ( unitAsArray(1) /= PHYQ_Temperature ) &
              & call Announce_error ( key, No_Error_code, &
              & 'Bad units for offsetAmount' )
            offsetAmount = valueAsArray(1)
          case ( f_orbitInclination ) ! For hydrostatic fill
            orbitinclInationVectorIndex = &
              & decoration(decoration(subtree(1,gson)))
            orbitinclInationQuantityIndex = &
              & decoration(decoration(decoration(subtree(2,gson))))
          case ( f_precision )      ! For masking l1b radiances
            precisionVectorIndex = decoration(decoration(subtree(1,gson)))
            precisionQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_precisionFactor )    ! For setting negative errors
            call expr ( subtree(2,subtree(j,key)), unitAsArray, valueAsArray )
            if ( unitAsArray(1) /= PHYQ_Dimensionless ) &
              & call Announce_error ( key, No_Error_code, &
              & 'Bad units for precisionFactor' )
            precisionFactor = valueAsArray(1)
          case ( f_profile )
            call expr ( gson , unitAsArray, valueAsArray )
            if ( all (unitAsArray /= (/PHYQ_Dimensionless, PHYQ_Invalid/) ) ) &
              call Announce_error ( key, 0, 'Bad units for profile' )
            profile = valueAsArray(1)
          case ( f_profileValues )
            valuesNode = subtree(j,key)
          case ( f_PtanQuantity ) ! For losGrid fill
            PtanVectorIndex = decoration(decoration(subtree(1,gson)))
            PtanQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_Phitan ) ! For losGrid fill
            PhitanVectorIndex = decoration(decoration(subtree(1,gson)))
            PhitanQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_phiWindow )
            call expr ( gson, unitAsArray, valueAsArray )
            phiWindow = valueAsArray(1)
            if ( all ( unitAsArray(1) /= (/ PHYQ_Profiles, PHYQ_Angle /) ) ) &
              call Announce_Error ( key, 0, 'Bad units for phiWindow' )
            phiWindowUnits = unitAsArray(1)
          case ( f_quantity )   ! What quantity are we filling quantity=vector.quantity
            vectorIndex = decoration(decoration(subtree(1,gson)))
            quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_radianceQuantity )      ! For estimated noise
            radianceVectorIndex = decoration(decoration(subtree(1,gson)))
            radianceQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_ratioQuantity )      ! For isotope ratio
            ratioVectorIndex = decoration(decoration(subtree(1,gson)))
            ratioQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_refract )
            refract = get_boolean ( gson )
          case ( f_refGPHQuantity ) ! For hydrostatic or rhi
            refGPHVectorIndex = decoration(decoration(subtree(1,gson)))
            refGPHQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))    
          case ( f_refGPHPrecisionQuantity ) ! For GPH precision
            refGPHPrecisionVectorIndex = decoration(decoration(subtree(1,gson)))
            refGPHPrecisionQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))    
          case ( f_resetSeed )
            resetSeed = get_boolean ( gson )
          case ( f_rhiQuantity ) ! For h2o from rhi
            sourceVectorIndex = decoration(decoration(subtree(1,gson)))
            sourceQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))    
          case ( f_scECI )                ! For special fill of losVel
            scECIVectorIndex = decoration(decoration(subtree(1,gson)))
            scECIQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_scVel )                ! For special fill of losVel
            scVelVectorIndex = decoration(decoration(subtree(1,gson)))
            scVelQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_scVelECI )                ! For special fill of losVel
            scVelVectorIndex = decoration(decoration(subtree(1,gson)))
            scVelQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_scVelECR )                ! For special fill of losVel
            scVelVectorIndex = decoration(decoration(subtree(1,gson)))
            scVelQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_seed ) ! For explicitly setting mls_random_seed
            seedNode=subtree(j,key)
          case ( f_sourceL2AUX )          ! Which L2AUXDatabase entry to use
            l2auxIndex = decoration(decoration(gson))
          case ( f_sourceL2GP )           ! Which L2GPDatabase entry to use
            l2gpIndex=decoration(decoration(gson))
          case ( f_sourceGrid )
            gridIndex=decoration(decoration(gson))
          case ( f_sourceQuantity )       ! When filling from a vector, what vector/quantity
            sourceVectorIndex = decoration(decoration(subtree(1,gson)))
            sourceQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_sourceVGrid )
            vGridIndex=decoration(decoration(gson))
          case ( f_spread ) ! For explicit fill, note that gson here is not same as others
            if ( node_id(gson) == n_set_one ) then
              spreadFlag = .true.
            else
              spreadFlag = decoration(subtree(2,gson)) == l_true
            end if
          case ( f_systemTemperature )
            sysTempVectorIndex = decoration(decoration(subtree(1,gson)))
            sysTempQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_tngtECI )              ! For special fill of losVel
            tngtECIVectorIndex = decoration(decoration(subtree(1,gson)))
            tngtECIQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_temperatureQuantity ) ! For hydrostatic or rhi
            temperatureVectorIndex = decoration(decoration(subtree(1,gson)))
            temperatureQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_tempPrecisionQuantity ) ! For rhi precision
            tempPrecisionVectorIndex = decoration(decoration(subtree(1,gson)))
            tempPrecisionQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_usb ) ! For folding
            usbVectorIndex = decoration(decoration(subtree(1,gson)))
            usbQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_usbFraction ) ! For folding
            usbFractionVectorIndex = decoration(decoration(subtree(1,gson)))
            usbFractionQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_vmrQuantity )     ! For special fill of columnAbundance
            vmrQtyVctrIndex = decoration(decoration(subtree(1,gson)))
            vmrQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          end select
        end do                  ! Loop over arguments to fill instruction

        ! Now call various routines to do the filling
        quantity => GetVectorQtyByTemplateIndex( &
          & vectors(vectorIndex), quantityIndex )

        select case ( fillMethod )
        case ( l_addNoise ) ! ----- Add random noise to source Quantity -------
          if (DEEBUG) call output('add noise method', advance='yes')
          if (.not. all(got( (/f_Quantity, f_sourceQuantity, f_noise/) ) ) ) &
            call Announce_error ( key, No_Error_code, &
             'Missing a required field to add noise'  )
          Quantity => GetVectorQtyByTemplateIndex( &
            & vectors(VectorIndex), QuantityIndex )
          sourceQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(sourceVectorIndex), sourceQuantityIndex )
          noiseQty => GetVectorQtyByTemplateIndex( &
            & vectors(noiseVectorIndex), noiseQtyIndex)
          math77_ran_pack = .not. switch2intrinsic
          if (DEEBUG) then
            call output('Switch to intrinsic? ', advance='no')
            call output(switch2intrinsic, advance='yes')
            call output('resetSeed? ', advance='no')
            call output(resetSeed, advance='yes')
            call output('got(f_seed)? ', advance='no')
            call output(got(f_seed), advance='yes')
          end if
          if ( resetSeed ) then
            call mls_random_seed(new_seed=seed(1:))
            if (DEEBUG) then
              call output('Letting mls choose new seed ', advance='no')
              call output(seed, advance='yes')
            end if
          else if ( got(f_seed) ) then
            do j=1, min(nsons(seedNode)-1, 2)
              call expr(subtree(j+1,seedNode),unitAsArray,valueAsArray)
              seed(j) = int(valueAsArray(1))+chunkNo
            end do
            if ( seed(1) /= 0 .and. seed(2) /= 0 ) then
              call mls_random_seed(pput=seed(1:))
              if (DEEBUG) then
                call output('Setting new seed ', advance='no')
                call output(seed, advance='yes')
              end if
            else
              call mls_random_seed(new_seed=seed(1:))
              if (DEEBUG) then
                call output('Letting mls choose new seed ', advance='no')
                call output(seed, advance='yes')
              end if
            end if
          else
            if (DEEBUG) then
              call mls_random_seed(gget=seed(1:))
              call output('Reusing current seed ', advance='no')
              call output(seed, advance='yes')
            end if
          end if

          ! Either multiplier = [a, b] or multiplier = b are possible
          if ( got(f_multiplier) ) then
            multiplier = UNDEFINED_VALUE
            do j=1, min(nsons(multiplierNode)-1, 2)
              call expr(subtree(j+1,multiplierNode),unitAsArray,valueAsArray)
              multiplier(j) = valueAsArray(1)
            end do
          else
            multiplier = 1.
          end if

          if (DEEBUG) then
            call output('Using multipliers: ', advance='no')
            call output(multiplier, advance='yes')
          end if
          call addGaussianNoise ( key, quantity, sourceQuantity, &
            & noiseQty, multiplier )

        case ( l_gphPrecision) ! -------------  GPH precision  -----
          ! Need a tempPrecision and a refgphPrecision quantity
          if ( .not.all(got( (/ f_refGPHPrecisionQuantity, f_tempPrecisionQuantity /))) ) &
            call Announce_Error ( key,needTempREFGPH )

          tempPrecisionQuantity => GetVectorQtyByTemplateIndex( &
            &  vectors(tempPrecisionVectorIndex), tempPrecisionQuantityIndex)
          if ( tempPrecisionQuantity%template%quantityType /= l_Temperature ) &
            & call Announce_Error ( key, badTemperatureQuantity )

          refGPHPrecisionQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(refGPHPrecisionVectorIndex), refGPHPrecisionQuantityIndex)
          if ( refGPHPrecisionQuantity%template%quantityType /= l_refGPH ) &
            & call Announce_Error ( key, badrefGPHQuantity )

          call FillGPHPrecision ( key, quantity, tempPrecisionQuantity, &
            & refGPHPrecisionQuantity )          

        case ( l_hydrostatic ) ! -------------  Hydrostatic fills  -----
          ! Need a temperature and a refgph quantity
          if ( .not.all(got( (/ f_refGPHQuantity, f_temperatureQuantity /))) ) &
            call Announce_Error ( key,needTempREFGPH )

          temperatureQuantity => GetVectorQtyByTemplateIndex( &
            &  vectors(temperatureVectorIndex), temperatureQuantityIndex)
          if ( temperatureQuantity%template%quantityType /= l_Temperature ) &
            & call Announce_Error ( key, badTemperatureQuantity )

          refGPHQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(refGPHVectorIndex), refGPHQuantityIndex)
          if ( refGPHQuantity%template%quantityType /= l_refGPH ) &
            & call Announce_Error ( key, badrefGPHQuantity )

          if ( quantity%template%quantityType==l_ptan ) then
            if ( .not. got(f_geocAltitudeQuantity) ) &
              & call Announce_Error ( key, needGeocAltitude )
            geocAltitudeQuantity => GetVectorQtyByTemplateIndex( &
              & vectors(geocAltitudeVectorIndex), geocAltitudeQuantityIndex)
            if ( geocAltitudeQuantity%template%quantityType /= l_tngtgeocAlt ) &
              & call Announce_Error ( key, badGeocAltitudeQuantity )

            if ( .not. got(f_phiTan) ) &
              & call Announce_Error ( key, 0, 'Needs phiTan quantity' )
            phiTanQuantity => GetVectorQtyByTemplateIndex( &
              & vectors(phiTanVectorIndex), phiTanQuantityIndex)
            if ( phiTanQuantity%template%quantityType /= l_phiTan ) &
              & call Announce_Error ( key, 0, 'Has a bad phiTan quantity' )

            if ( .not. got(f_h2oQuantity) ) &
              & call Announce_Error ( key, needH2O )
            h2oQuantity => GetVectorQtyByTemplateIndex( &
              & vectors(h2oVectorIndex), h2oQuantityIndex)
            if ( .not. ValidateVectorQuantity(h2oQuantity, &
              & quantityType=(/l_vmr/), molecule=(/l_h2o/)) )&
              & call Announce_Error ( key, badGeocAltitudeQuantity )

            if ( .not. got(f_orbitInclination) ) &
              & call Announce_Error ( key, needOrbitInclination )
            orbitInclinationQuantity => GetVectorQtyByTemplateIndex( &
              & vectors(orbitInclinationVectorIndex), orbitInclinationQuantityIndex)

          else
            nullify ( geocAltitudeQuantity, h2oQuantity )
          end if
          call FillVectorQtyHydrostatically ( key, quantity, temperatureQuantity, &
            & refGPHQuantity, h2oQuantity, orbitInclinationQuantity, &
            & phiTanQuantity, geocAltitudeQuantity, maxIterations, &
            & phiWindow, phiWindowUnits )          

        case ( l_isotope ) ! --------------- Isotope based fills -------
          if (.not. all(got( (/f_ratioQuantity, f_sourceQuantity/) ) ) ) &
            & call Announce_Error ( key, badIsotopeFill )
          ratioQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(ratioVectorIndex), ratioQuantityIndex )
          sourceQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(sourceVectorIndex), sourceQuantityIndex )
          call FillVectorQtyFromIsotope ( key, quantity, sourceQuantity, &
            & ratioQuantity )

        case ( l_manipulate ) ! ---------------------------- Manipulate --
          if ( .not. got ( f_a ) ) &
            & call Announce_error ( key, 0 ,'aQuantity not supplied' )
          if ( .not. got ( f_manipulation ) ) &
            & call Announce_error ( key, 0 ,'manipulation not supplied' )
          aQuantity => GetVectorQtyByTemplateIndex ( &
            & vectors(aVecIndex), aQtyIndex )
          if ( got ( f_b ) ) then
            bQuantity => GetVectorQtyByTemplateIndex ( &
              & vectors(bVecIndex), bQtyIndex )
          else
            nullify ( bQuantity )
          end if
          call FillQuantityByManipulation ( quantity, aQuantity, bQuantity, &
            & manipulation, key )

        case ( l_magneticModel ) ! --------------------- Magnetic Model --
          if ( .not. got ( f_gphQuantity ) ) then
            call Announce_Error ( key, 0, 'Need gphQuantity for magnetic model' )
          else
            GPHQuantity => GetVectorQtyByTemplateIndex( &
              & vectors(GPHVectorIndex), GPHQuantityIndex)
            call FillQuantityUsingMagneticModel ( quantity, gphQuantity, key )          
          end if
          
        case ( l_offsetRadiance ) ! ------------------- Offset radiance --
          if ( .not. got ( f_radianceQuantity ) ) &
            & call Announce_error ( key, 0, 'radianceQuantity not supplied' )
          radianceQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(radianceVectorIndex), radianceQuantityIndex )
          call OffsetRadianceQuantity ( quantity, radianceQuantity, offsetAmount )

        case ( l_profile ) ! ------------------------ Profile fill -------
          if ( .not. got ( f_profileValues ) ) &
            call Announce_error ( key, 0, 'profileValues not supplied' )
          if ( .not. got ( f_instances ) ) instancesNode = 0
          call FillVectorQtyFromProfile ( key, quantity, valuesNode, &
            & instancesNode, vectors(vectorIndex)%globalUnit, dontMask )

        case ( l_refract )              ! --------- refraction for phiTan -----
          if ( refract ) then 
            if ( .not. all ( got ( &
            & (/ f_temperatureQuantity, f_h2oQuantity, f_ptanQuantity /) ) ) ) &
            & call Announce_error ( key, badRefractFill )
            temperatureQuantity => GetVectorQtyByTemplateIndex( &
              & vectors(temperatureVectorIndex), temperatureQuantityIndex)
            h2oQuantity => GetVectorQtyByTemplateIndex( &
              & vectors(h2oVectorIndex), h2oQuantityIndex)
            ptanQuantity => GetVectorQtyByTemplateIndex( &
              & vectors(ptanVectorIndex), ptanQuantityIndex)
          else
            temperatureQuantity => NULL()
            h2oQuantity => NULL()
            ptanQuantity => NULL()
          end if
          call FillPhiTanWithRefraction ( key, quantity, ptanQuantity, temperatureQuantity, &
            & h2oQuantity, refract )

        case ( l_rectanglefromlos ) ! -------fill from losGrid quantity -------
          if (.not. all(got((/f_losQty,f_earthRadius,f_PtanQuantity/))))&
            & call Announce_Error ( key, badlosGridFill )
          earthRadiusQty => GetVectorQtyByTemplateIndex( &
            & vectors(earthRadiusVectorIndex), earthRadiusQtyIndex )
          PtanQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(PtanVectorIndex), PtanQuantityIndex )
          losQty => GetVectorQtyByTemplateIndex( &
            & vectors(losVectorIndex), losQtyIndex )
          call FillQuantityFromLosGrid ( key, Quantity, losQty, &
            & ptanQuantity, earthRadiusQty, &
            & noFineGrid, extinction, errorCode )

        case ( l_RHIfromH2O ) ! -------fill RHI from H2O quantity -------
            if ( .not. any(got( &
             & (/f_h2oquantity, f_temperatureQuantity/) &
             & )) ) then
              call Announce_error ( key, No_Error_code, &
              & 'Missing a required field to fill rhi'  )
            else ! value
                h2oQuantity => GetVectorQtyByTemplateIndex( &
                  & vectors(h2oVectorIndex), h2oQuantityIndex)
                temperatureQuantity => GetVectorQtyByTemplateIndex( &
                  & vectors(temperatureVectorIndex), temperatureQuantityIndex)
                if ( .not. ValidateVectorQuantity(h2oQuantity, &
                  & quantityType=(/l_vmr/), molecule=(/l_h2o/)) ) then
                  call Announce_Error ( key, No_Error_code, &
                    & 'The h2oQuantity is not a vmr for the H2O molecule'  )
                else if ( .not. ValidateVectorQuantity(temperatureQuantity, &
                  & quantityType=(/l_temperature/)) ) then
                  call Announce_Error ( key, No_Error_code, &
                    & 'The temperatureQuantity is not a temperature'  )
                else
                  call FillRHIFromH2O ( key, quantity, &
                    & h2oQuantity, temperatureQuantity, &
                    & dontMask, ignoreZero, ignoreNegative, interpolate, &
                    & .true., &   ! Mark Undefined values?
                    & invert )    ! invert rather than convert?
                end if
            end if
        case ( l_H2OfromRHI ) ! -------fill H2O from RHI quantity -------
            if ( .not. any(got( &
             & (/f_rhiQuantity, f_temperatureQuantity/) &
             & )) ) then
              call Announce_error ( key, No_Error_code, &
              & 'Missing a required field to fill h2o from rhi'  )
            else
              sourceQuantity => GetVectorQtyByTemplateIndex( &
                & vectors(sourceVectorIndex), sourceQuantityIndex)
              temperatureQuantity => GetVectorQtyByTemplateIndex( &
                & vectors(temperatureVectorIndex), temperatureQuantityIndex)
              if ( .not. ValidateVectorQuantity(sourceQuantity, &
                & quantityType=(/l_rhi/)) ) then
                call Announce_Error ( key, No_Error_code, &
                & 'The rhiQuantity is not an rhi'  )
              else if ( .not. ValidateVectorQuantity(Quantity, &
                & quantityType=(/l_vmr/), molecule=(/l_h2o/)) ) then
                call Announce_Error ( key, No_Error_code, &
                & 'The Quantity is not a vmr for the H2O molecule'  )
              else if ( .not. ValidateVectorQuantity(temperatureQuantity, &
                & quantityType=(/l_temperature/)) ) then
                call Announce_Error ( key, No_Error_code, &
                & 'The temperatureQuantity is not a temperature'  )
              else
                invert = .true.
                call FillRHIFromH2O ( key, quantity, &
                & sourceQuantity, temperatureQuantity, &
                & dontMask, ignoreZero, ignoreNegative, interpolate, &
                & .true., &   ! Mark Undefined values?
                & invert )    ! invert rather than convert?
              end if
            end if
        case ( l_RHIPrecisionfromH2O ) ! --fill RHI prec. from H2O quantity --
                if ( .not. any(got( &
                  & (/f_h2oquantity, f_temperatureQuantity, &
                  & f_h2oPrecisionquantity, f_tempPrecisionQuantity/) &
                  & )) ) then
                  call Announce_error ( key, No_Error_code, &
                    & 'Missing a required field to fill rhi precision'  )
                else
                  h2oQuantity => GetVectorQtyByTemplateIndex( &
                    & vectors(h2oVectorIndex), h2oQuantityIndex)
                  temperatureQuantity => GetVectorQtyByTemplateIndex( &
                    & vectors(temperatureVectorIndex), temperatureQuantityIndex)
                  h2oPrecisionQuantity => GetVectorQtyByTemplateIndex( &
                    & vectors(h2oPrecisionVectorIndex), h2oPrecisionQuantityIndex)
                  tempPrecisionQuantity => GetVectorQtyByTemplateIndex( &
                    & vectors(tempPrecisionVectorIndex), tempPrecisionQuantityIndex)
                  if ( .not. ValidateVectorQuantity(h2oQuantity, &
                    & quantityType=(/l_vmr/), molecule=(/l_h2o/)) ) then
                    call Announce_Error ( key, No_Error_code, &
                      & 'The h2oQuantity is not a vmr for the H2O molecule'  )
                  else if ( .not. ValidateVectorQuantity(temperatureQuantity, &
                    & quantityType=(/l_temperature/)) ) then
                    call Announce_Error ( key, No_Error_code, &
                      & 'The temperatureQuantity is not a temperature'  )
                  else if ( .not. ValidateVectorQuantity(h2oPrecisionQuantity, &
                    & quantityType=(/l_vmr/), molecule=(/l_h2o/)) ) then
                    call Announce_Error ( key, No_Error_code, &
                      & 'The h2oPrecisionQuantity is not a vmr for the H2O molecule'  )
                  else if ( .not. ValidateVectorQuantity(tempPrecisionQuantity, &
                    & quantityType=(/l_temperature/)) ) then
                    call Announce_Error ( key, No_Error_code, &
                      & 'The tempPrecisionQuantity is not a temperature'  )
                    ! This is not the right way to do an invert fill
                    ! The first quantity named on the fill line must
                    ! _always_ be the one getting filled according to
                    ! the pattern '[verb] [direct-object] [modifier(s)]'
                    ! see case l_vmr below for how to do this
                    !else if ( invert ) then
                    !  call FillRHIFromH2O ( key, h2oquantity, &
                    !  & Quantity, temperatureQuantity, &
                    !  & dontMask, ignoreZero, ignoreNegative, interpolate, &
                    !  & .true., &   ! Mark Undefined values?
                    !  & invert )    ! invert rather than convert?
                  else
                    call FillRHIPrecisionFromH2O ( key, quantity, &
                      & h2oPrecisionQuantity, tempPrecisionQuantity, h2oQuantity, temperatureQuantity, &
                      & dontMask, ignoreZero, ignoreNegative, interpolate, &
                      & .true., &   ! Mark Undefined values?
                      & invert )    ! invert rather than convert?
                  end if
                end if
        case ( l_special ) ! -  Special fills for some quantities  -----
          ! Either multiplier = [a, b] or multiplier = b are possible
          if ( got(f_multiplier) ) then
            multiplier = UNDEFINED_VALUE
            do j=1, min(nsons(multiplierNode)-1, 2)
              call expr(subtree(j+1,multiplierNode),unitAsArray,valueAsArray)
              multiplier(j) = valueAsArray(1)
            end do
          else
            multiplier = 1.
          end if

          if (DEEBUG) then
            call output('Using multipliers: ', advance='no')
            call output(multiplier, advance='yes')
          end if

          select case ( quantity%template%quantityType )
          case ( l_losVel )
            if ( .not. any(got( &
            & (/f_tngtECI, f_scECI, f_scVel, f_scVelECI, f_scVelECR/) )) ) then
              call Announce_error ( key, badlosVelFill )
            else
              tngtECIQuantity => GetVectorQtyByTemplateIndex( &
                & vectors(tngtECIVectorIndex), tngtECIQuantityIndex)
              scECIQuantity => GetVectorQtyByTemplateIndex( &
                & vectors(scECIVectorIndex), scECIQuantityIndex)
              scVelQuantity => GetVectorQtyByTemplateIndex( &
                & vectors(scVelVectorIndex), scVelQuantityIndex)
              call FillLOSVelocity ( key, quantity, tngtECIQuantity, &
                & scECIquantity, scVelQuantity )
            end if
          case ( l_columnAbundance )
            if ( .not. any(got( (/f_vmrQuantity, f_boundaryPressure/) )) ) then
              call Announce_error ( key, No_Error_code, &
              & 'Missing a required field to fill column abundance'  )
            else
              bndPressQty => GetVectorQtyByTemplateIndex( &
                & vectors(bndPressVctrIndex), bndPressQtyIndex)
              vmrQty => GetVectorQtyByTemplateIndex( &
                & vectors(vmrQtyVctrIndex), vmrQtyIndex)
              call FillColAbundance ( key, quantity, &
                & bndPressQty, vmrQty )
            end if
          case ( l_chiSQChan )
            if ( .not. any(got( (/f_measurements, f_model, f_noise/) )) ) then
              call Announce_error ( key, No_Error_code, &
              & 'Missing a required field to fill chi^2 on channels'  )
            else
              measQty => GetVectorQtyByTemplateIndex( &
                & vectors(measVectorIndex), measQtyIndex)
              modelQty => GetVectorQtyByTemplateIndex( &
                & vectors(modelVectorIndex), modelQtyIndex)
              noiseQty => GetVectorQtyByTemplateIndex( &
                & vectors(noiseVectorIndex), noiseQtyIndex)
              call FillChiSqChan ( key, quantity, &
                & measQty, modelQty, noiseQty, &
                & dontMask, ignoreZero, ignoreNegative, multiplier )
            end if
          case ( l_chiSQMMaf )
            if ( .not. any(got( (/f_measurements, f_model, f_noise/) )) ) then
              call Announce_error ( key, No_Error_code, &
              & 'Missing a required field to fill chi^2 on MAFs'  )
            else
              measQty => GetVectorQtyByTemplateIndex( &
                & vectors(measVectorIndex), measQtyIndex)
              modelQty => GetVectorQtyByTemplateIndex( &
                & vectors(modelVectorIndex), modelQtyIndex)
              noiseQty => GetVectorQtyByTemplateIndex( &
                & vectors(noiseVectorIndex), noiseQtyIndex)
              call FillChiSqMMaf ( key, quantity, &
                & measQty, modelQty, noiseQty, &
                & dontMask, ignoreZero, ignoreNegative, multiplier )
            end if
          case ( l_chiSQMMif )
            if ( .not. any(got( (/f_measurements, f_model, f_noise/) )) ) then
              call Announce_error ( key, No_Error_code, &
              & 'Missing a required field to fill chi^2 on MIFs'  )
            else
              measQty => GetVectorQtyByTemplateIndex( &
                & vectors(measVectorIndex), measQtyIndex)
              modelQty => GetVectorQtyByTemplateIndex( &
                & vectors(modelVectorIndex), modelQtyIndex)
              noiseQty => GetVectorQtyByTemplateIndex( &
                & vectors(noiseVectorIndex), noiseQtyIndex)
              call FillChiSqMMif ( key, quantity, &
                & measQty, modelQty, noiseQty, &
                & dontMask, ignoreZero, ignoreNegative, multiplier )
            end if
          case default
            call Announce_error ( key, noSpecialFill )
          end select

        case ( l_negativePrecision ) ! ------------ Set output SD -ve wrt apriori
          if ( .not. got ( f_aprioriPrecision ) ) &
            & call Announce_Error ( key, No_Error_code, &
            & 'Missing aprioriPrecision field for negativePrecision fill' )
          aprioriPrecision => GetVectorQtyByTemplateIndex ( &
            & vectors(aprPrecVctrIndex), aprPrecQtyIndex )
          where ( quantity%values >= aprioriPrecision%values*precisionFactor )
            quantity%values = - quantity%values
          end where

        case ( l_vGrid ) ! ---------------------  Fill from vGrid  -----
          if (.not. ValidateVectorQuantity(quantity, &
            & quantityType=(/l_ptan/), &
            & frequencyCoordinate=(/l_none/) ) ) &
            & call Announce_Error ( key, No_Error_code, &
            &   'vGrids can only be used to fill ptan quantities' )
          if ( vGrids(vGridIndex)%verticalCoordinate /= l_zeta ) &
            & call Announce_Error ( key, No_Error_code, &
            &  'Vertical coordinate in vGrid is not zeta' )
          if ( vGrids(vGridIndex)%noSurfs /= quantity%template%noSurfs )&
            & call Announce_Error ( key, No_Error_code, &
            &  'VGrid is not of the same size as the quantity' )
          quantity%values = spread ( vGrids(vGridIndex)%surfs, 2, &
            & quantity%template%noInstances )

        case ( l_fold ) ! --------------- Fill by sideband folding -----
          lsb => GetVectorQtyByTemplateIndex ( &
            & vectors(lsbVectorIndex), lsbQuantityIndex )
          lsbFraction => GetVectorQtyByTemplateIndex ( &
            & vectors(lsbFractionVectorIndex), lsbFractionQuantityIndex )
          usb => GetVectorQtyByTemplateIndex ( &
            & vectors(usbVectorIndex), usbQuantityIndex )
          usbFraction => GetVectorQtyByTemplateIndex ( &
            & vectors(usbFractionVectorIndex), usbFractionQuantityIndex )
          call FillFoldedRadiance ( quantity, lsb, usb, lsbFraction, usbFraction, key )

        case ( l_gridded ) ! ------------  Fill from gridded data  -----
          if ( .not. got(f_sourceGrid) ) &
            & call Announce_Error ( key,noSourceGridGiven )
          call FillVectorQuantityFromGrid &
            & ( quantity, griddedDataBase(gridIndex), errorCode )
          if ( errorCode /= 0 ) call Announce_error ( key, errorCode )

        case ( l_l2gp ) ! --------------  Fill from L2GP quantity  -----
          if ( .NOT. got(f_sourceL2GP) ) &
            & call Announce_Error ( key, noSourceL2GPGiven )
          call FillVectorQuantityFromL2GP &
            & ( quantity, l2gpDatabase(l2gpIndex), interpolate, profile, errorCode )
          if ( errorCode /= 0 ) call Announce_error ( key, errorCode )

        case ( l_l2aux ) ! ------------  Fill from L2AUX quantity  -----
          if ( .NOT. got(f_sourceL2AUX) ) &
            & call Announce_Error ( key, noSourceL2AUXGiven )
          call FillVectorQuantityFromL2AUX(quantity,l2auxDatabase(l2auxIndex),errorCode)
          if ( errorCode /= 0 ) call Announce_error ( key, errorCode )

        case ( l_vector ) ! ---------------- Fill from another qty.
          ! This is VERY PRELIMINARY, A more fancy one needs to be written
          ! before too long.
          ! Note that this does *NOT* copy the mask (at least for the moment)
          ! It is assumed that the original one (e.g. inherited from transfer)
          ! is still relevant.
          if ( .not. got(f_sourceQuantity) ) &
            & call Announce_Error ( key, No_Error_Code, &
            & 'Missing a source field for vector fill' )
          Quantity => GetVectorQtyByTemplateIndex( &
            & vectors(VectorIndex), QuantityIndex )
          sourceQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(sourceVectorIndex), sourceQuantityIndex )
          if ( quantity%template%name /= sourceQuantity%template%name ) then
            if ( .not. interpolate) then
              call Announce_Error ( key, No_Error_Code, &
                & 'Quantity and sourceQuantity do not have the same template' )
            else
              call FillQtyFromInterpolatedQty ( quantity, sourceQuantity, key )
            end if
          else
            ! Just a straight copy
            ! If we have a mask and we're going to obey it then do so
            if ( associated(quantity%mask) .and. .not. dontMask ) then
              where ( iand ( ichar(quantity%mask(:,:)), m_Fill ) == 0 )
                quantity%values(:,:) = sourceQuantity%values(:,:)
              end where
            else ! Otherwise, just blindly copy
              quantity%values = sourceQuantity%values
            end if
          end if

        case ( l_estimatedNoise ) ! ----------- Fill with estimated noise ---
          if (.not. all(got( (/ f_radianceQuantity, &
            & f_systemTemperature, f_integrationTime /)))) &
            & call Announce_Error ( key, badEstNoiseFill )
          radianceQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(radianceVectorIndex), radianceQuantityIndex )
          sysTempQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(sysTempVectorIndex), sysTempQuantityIndex )
          if ( got ( f_noiseBandwidth ) ) then
            nbwQuantity => GetVectorQtyByTemplateIndex( &
              & vectors(nbwVectorIndex), nbwQuantityIndex )
          else
            nbwQuantity => NULL()
          end if
          call FillVectorQtyWithEstNoise ( &
            & quantity, radianceQuantity, sysTempQuantity, nbwQuantity, &
            & integrationTime )

        case ( l_explicit ) ! ---------  Explicity fill from l2cf  -----
          if ( .not. got(f_explicitValues) ) &
            & call Announce_Error ( key, noExplicitValuesGiven )
          call ExplicitFillVectorQuantity ( quantity, valuesNode, spreadFlag, &
            & vectors(vectorIndex)%globalUnit, dontmask )

        case ( l_l1b ) ! --------------------  Fill from L1B data  -----
          if ( got(f_precision) ) then
            precisionQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(precisionVectorIndex), precisionQuantityIndex )
            call FillVectorQuantityFromL1B ( key, quantity, chunks(chunkNo), &
              & l1bInfo, isPrecision, precisionQuantity )
          else
            call FillVectorQuantityFromL1B ( key, quantity, chunks(chunkNo), &
              & l1bInfo, isPrecision )
          end if

        case default
          call Announce_error ( key, 0, 'This fill method not yet implemented' )
        end select

      case ( s_FillCovariance ) ! ===============  FillCovariance  =====
        invert = .false. ! Default if the field isn't present
        lengthScale = 0
        fraction = 0
        do j = 2, nsons(key)
          gson = subtree(j,key) ! The argument
          fieldIndex = get_field_id(gson)
          if (nsons(gson) > 1) &
            & gson = decoration(decoration(subtree(2,gson))) ! Now value of said argument
          got(fieldIndex)=.true.
          select case ( fieldIndex )
          case ( f_matrix )
            matrixToFill = gson
            if ( getKindFromMatrixDatabase(matrices(matrixToFill)) /= k_spd ) &
              call announce_error ( key, notSPD )
          case ( f_diagonal )
            diagonal = gson
          case ( f_lengthScale )
            lengthScale = gson
          case ( f_fraction )
            fraction = gson
            call announce_error ( key, notImplemented, "Decay" ) !???
          case ( f_invert )
            invert = get_boolean ( subtree(j,key) )
          case ( f_superDiagonal )
            superDiagonal = gson
            call announce_error ( key, notImplemented, "SuperDiagonal" ) !???
          end select
        end do

        call getFromMatrixDatabase ( matrices(matrixToFill), covariance )
        call FillCovariance ( covariance, vectors, diagonal, lengthScale, fraction, &
          & invert )

      case ( s_FillDiagonal ) ! ===============  FillDiagonal  =====
        do j = 2, nsons(key)
          gson = subtree(j,key) ! The argument
          fieldIndex = get_field_id(gson)
          if (nsons(gson) > 1) &
            & gson = decoration(decoration(subtree(2,gson))) ! Now value of said argument
          got(fieldIndex)=.true.
          select case ( fieldIndex )
          case ( f_matrix )
            matrixToFill = gson
            if ( getKindFromMatrixDatabase(matrices(matrixToFill)) /= k_spd ) &
              call announce_error ( key, notSPD )
          case ( f_diagonal )
            diagonal = gson
          end select
        end do

        call getFromMatrixDatabase ( matrices(matrixToFill), covariance )
        call getDiagonal( covariance%m, vectors(diagonal))

      ! End of fill operations

      case ( s_destroy ) ! ===============================  Destroy ==
        if (DEEBUG) call output('Destroy vector/matrix instruction', &
        &  advance='no')
        ! Here we're to try to shrink the vector database by destroying a vector
        ! or the matrix database database by destroying a matrix
        ! Loop over the instructions
        ! (Shall we allow multiple vectors on a single line? Maybe later)
        do j = 2, nsons(key)
          son = subtree(j,key)  ! The argument
          fieldIndex = get_field_id(son)
          got(fieldIndex)=.true.
          if ( nsons(son) > 1 ) then
            fieldValue = decoration(subtree(2,son)) ! The field's value
          else
            fieldValue = son
          end if
          select case ( fieldIndex )
          case ( f_matrix )
            matrixToKill = decoration(decoration(subtree(2,son)))
          case ( f_vector )
            sourceVectorIndex = decoration(decoration(subtree(2,son)))
          case default ! Can't get here if type checker worked
          end select
        end do
        if ( got(f_vector) ) then
          if (DEEBUG) then
            if ( vectors(sourceVectorIndex)%name /= 0 ) then
              call output ( '   Vector Name = ' )
              call display_string ( vectors(sourceVectorIndex)%name )
            end if
            if ( vectors(sourceVectorIndex)%template%name /= 0 ) then
              call output ( ' Template_Name = ' )
              call display_string ( vectors(sourceVectorIndex)%template%name )
              call output ( ' ', advance='yes' )
            end if
            call output ( ' -- vector database before removal --', advance='yes' )
            call dump(vectors, details=-2)
          end if

          call DestroyVectorInfo ( vectors(sourceVectorIndex) )
     !    vectorindex = rmVectorFromDatabase ( vectors, vectors(sourceVectorIndex) )
          if (DEEBUG) then
            call output ( ' -- vector database after removal --', advance='yes' )
            call dump(vectors, details=-2)
          end if
        else if ( got(f_matrix) ) then
          if (DEEBUG) then
            ! if ( matrices(matrixToKill)%matrix%name /= 0 ) then
             ! call output ( '   Matrix Name = ' )
             ! call display_string ( matrices(matrixToKill)%matrix%name )
             call dump(matrices(matrixToKill), -1)
            ! end if
          end if

          call DestroyMatrix ( matrices(matrixToKill) )
        end if

      case ( s_negativePrecision ) ! ===============================  Transfer ==
        ! Here we're on a setNegative instruction
        ! Loop over the instructions
        skipMask = .false.
        do j = 2, nsons(key)
          gson = subtree(j,key)  ! The argument
          fieldIndex = get_field_id(gson)
          if ( nsons(gson) > 1 ) then
            gson = subtree(2,gson) ! Now the value of said argument
            fieldValue = decoration(gson) ! The field's value
          else
            fieldValue = gson
          end if
          select case ( fieldIndex )
          case ( f_precision )
            precisionVectorIndex = decoration(fieldValue)
          case ( f_aprioriPrecision )
            aprPrecVctrIndex =  decoration(fieldValue)
          case ( f_precisionFactor )
            call expr ( gson, unitAsArray, valueAsArray )
            if ( unitAsArray(1) /= PHYQ_Dimensionless ) &
              & call Announce_error ( key, No_Error_code, &
              & 'Bad units for precisionFactor' )
            precisionFactor = valueAsArray(1)
          case default ! Can't get here if type checker worked
          end select
        end do
        if ( vectors(precisionVectorIndex)%template%name /= &
          & vectors(aprPrecVctrIndex)%template%name ) then
          call Announce_error ( key, no_error_code, &
            & 'Incompatible vectors in negativePrecision statement' )
        else
          do j = 1, vectors(precisionVectorIndex)%template%noQuantities
            where ( vectors(precisionVectorIndex)%quantities(j)%values > &
              & vectors(aprPrecVctrIndex)%quantities(j)%values * precisionFactor )
              vectors(precisionVectorIndex)%quantities(j)%values = &
                & - vectors(precisionVectorIndex)%quantities(j)%values
            end where
          end do
        end if
      
      case ( s_transfer ) ! ===============================  Transfer ==
        ! Here we're on a transfer instruction
        ! Loop over the instructions
        skipMask = .false.
        do j = 2, nsons(key)
          gson = subtree(j,key)  ! The argument
          fieldIndex = get_field_id(gson)
          if ( nsons(gson) > 1 ) then
            gson = subtree(2,gson) ! Now the value of said argument
            fieldValue = decoration(gson) ! The field's value
          else
            fieldValue = gson
          end if
          select case ( fieldIndex )
          case ( f_source )
            sourceVectorIndex = decoration(fieldValue)
          case ( f_destination )
            destinationVectorIndex = decoration(fieldValue)
          case ( f_skipMask )
            skipMask = get_boolean ( fieldValue )
          case default ! Can't get here if type checker worked
          end select
        end do
        call TransferVectors ( vectors(sourceVectorIndex), &
          & vectors(destinationVectorIndex), skipMask )

      case ( s_time ) ! ===================================  Time  =====
        if ( timing ) then
          call sayTime
        else
          call time_now ( t1 )
          timing = .true.
        end if

      case ( s_snoop )
        ! Generate a matrix database
        noSnoopedMatrices = 0
        if ( associated ( matrices ) ) then
          do j = 1, size ( matrices )
            call GetActualMatrixFromDatabase ( matrices(j), oneMatrix )
            if ( associated ( oneMatrix ) ) &
              & noSnoopedMatrices = noSnoopedMatrices + 1
          end do
        end if
        allocate ( snoopMatrices ( noSnoopedMatrices ), STAT=status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_Allocate//'snoopMatrices' )
        if ( associated ( matrices ) ) then
          noSnoopedMatrices = 1
          do j = 1, size ( matrices )
            call GetActualMatrixFromDatabase ( matrices(j), oneMatrix )
            if ( associated ( oneMatrix ) ) then 
              snoopMatrices(noSnoopedMatrices) = oneMatrix
              noSnoopedMatrices = noSnoopedMatrices + 1
            end if
          end do
        end if
        call Snoop ( key=key, vectorDatabase=vectors, matrixDatabase=snoopMatrices )
        deallocate ( snoopMatrices, STAT=status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_Deallocate//'snoopMatrices' )

      case default ! Can't get here if tree_checker worked correctly
      end select
    end do

    if ( ERROR /= 0 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, 'Problem with Fill section' )
    end if

    if ( toggle(gen) ) then
      if ( levels(gen) > 0 ) then
        call dump ( vectors, details=levels(gen)-1 )
        call dump ( matrices, details=levels(gen)-1 )
      end if
      call trace_end ( "MLSL2Fill" )
    end if
    math77_ran_pack = old_math77_ran_pack
    if ( timing ) call sayTime

  ! =====     Internal Procedures     ==================================

  contains

    ! --------------------------------------------------  SayTime  -----
    subroutine SayTime
      call time_now ( t2 )
      if ( total_times ) then
        call output ( "Total time = " )
        call output ( dble(t2), advance = 'no' )
        call blanks ( 4, advance = 'no' )
      end if
      call output ( "Timing for MLSL2Fill = " )
      call output ( dble(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime

    ! ------------------------------------------- addGaussianNoise ---
    subroutine addGaussianNoise ( key, quantity, sourceQuantity, &
              & noiseQty, multiplier )
      ! A special fill: quantity = sourceQuantity + g() noiseQty
      ! where g() is a random number generator with mean 0 and std. dev. 1
      ! Generalized into ( a sourceQuantity + b g() noiseQty )
      ! where a and b are multipliers)
      ! Formal arguments
      integer, intent(in) :: KEY
      type (VectorValue_T), intent(inout) ::      quantity
      type (VectorValue_T), intent(in) ::         sourceQuantity
      type (VectorValue_T), intent(in) ::         noiseQty
      real, dimension(:), intent(in), optional :: multiplier

      ! Local variables
      integer                          ::    ROW, COLUMN
      real                             ::    a, b

      ! Executable code
      ! First check that things are OK.
      if (.not. FillableChiSq ( quantity, &
        & sourceQuantity, noiseQty ) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Incompatibility among vector quantities adding noise'  )
        return
      end if

     ! Either multiplier = [a, b] or multiplier = b are possible
      if ( .not. present(multiplier) ) then
        a = 1.
        b = 1.
      else if (&
      & multiplier(1) == UNDEFINED_VALUE .and. multiplier(2) == UNDEFINED_VALUE &
      & ) then
        a = 1.
        b = 1.
      else if ( multiplier(2) == UNDEFINED_VALUE ) then
        a = 1.
        b = multiplier(1)
      else
        a = multiplier(1)
        b = multiplier(2)
      end if

      do column=1, size(quantity%values(1, :))
        do row=1, size(quantity%values(:, 1))
          quantity%values(row, column) = &
            & sourceQuantity%values(row, column) * a &
            & + &
            & drang() * noiseQty%values(row, column) * b
        end do
      end do

    end subroutine addGaussianNoise

    !------------------------------------- FillCovariance ------------

    subroutine FillCovariance ( covariance, vectors, diagonal, &
      & lengthScale, fraction, invert )
      ! This routine fills a covariance matrix from a given set of vectors
      type (Matrix_SPD_T), intent(inout) :: COVARIANCE ! The matrix to fill
      type (Vector_T), dimension(:), intent(in), target :: VECTORS ! The vector database
      integer, intent(in) :: DIAGONAL     ! Index of vector describing diagonal
      integer, intent(in) :: LENGTHSCALE  ! Index of vector describing length scale
      integer, intent(in) :: FRACTION     ! Index of vector describing fraction
      logical, intent(in) :: INVERT       ! We actually want the inverse

      ! Local parameters
      real(r8), parameter :: DECADE = 16000.0 ! Number of meters per decade.

      ! Local variables
      type (VectorValue_T), pointer :: d  ! Diagonal
      type (VectorValue_T), pointer :: l  ! Length
      type (VectorValue_T), pointer :: f  ! Fraction
      type (QuantityTemplate_T), pointer :: qt ! One quantity template
      type (Vector_T) :: DMASKED ! Masked diagonal
      type (Vector_T) :: LMASKED ! Masked length scale
      integer :: B                        ! Block index
      integer :: I                        ! Instance index
      integer :: J                        ! Loop index
      integer :: K                        ! Loop index
      integer :: N                        ! Size of matrix block
      integer :: Q                        ! Quantity index
      real (rm), dimension(:,:), pointer :: M ! The matrix being filled
      real (r8), dimension(:), pointer :: SURFS ! The vertical coordinate
      real (r8) :: distance               ! Distance between two points
      real (r8) :: meanLength             ! Geometric mean length scale
      real (r8) :: meanDiag               ! Geometric mean diagonal value
      real (r8) :: thisFraction           ! Geometric mean diagonal value
      logical, dimension(:), pointer :: condition ! Condition

      ! Executable code

      ! Apply mask to diagonal
      nullify ( m, condition )
      call CopyVector ( Dmasked, vectors(diagonal), clone=.true., &
        & vectorNameText='_Dmasked' )
      call ClearUnderMask ( Dmasked )

      if ( lengthScale == 0 ) then
        call updateDiagonal ( covariance, vectors(diagonal), square=.true., &
          & invert=invert )

      else ! Do a more complex fill.

        ! Setup some stuff
        call CopyVector ( Lmasked, vectors(lengthScale), clone=.true., &
          & vectorNameText='_Lmasked' ) 
        call ClearUnderMask ( Lmasked )

        ! Check the validity of the supplied vectors
        if ( covariance%m%row%vec%template%id /= &
          & dMasked%template%id ) call MLSMessage ( MLSMSG_Error, &
          & ModuleName, "diagonal and covariance not compatible in fillCovariance" )
        if ( covariance%m%row%vec%template%id /= &
          & lMasked%template%id ) call MLSMessage ( MLSMSG_Error, &
          & ModuleName, "lengthScale and covariance not compatible in fillCovariance" )
        if ( lMasked%globalUnit /= phyq_length ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "length vector does not have dimensions of length" )

        if ( fraction /= 0 ) then
          if ( covariance%m%row%vec%template%id /= &
            & vectors(fraction)%template%id ) call MLSMessage ( MLSMSG_Error, &
            & ModuleName, "fraction and covariance not compatible in fillCovariance" )
          if ( vectors(fraction)%globalUnit /= phyq_dimensionless ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "fraction vector is not dimensionless" )
        else
          thisFraction = 1.0
        end if

        ! Now loop over the quantities
        do q = 1, covariance%m%col%vec%template%noQuantities

          ! Setup pointers etc.
          d => dMasked%quantities(q)
          qt => d%template
          l => lMasked%quantities(q)
          if (fraction /=0) f => vectors(fraction)%quantities(q)
          n = qt%instanceLen
          if ( qt%coherent ) surfs => qt%surfs(:,1)
          if ( .not. qt%regular ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "Unable to handle irregular quantity in FillCovariance" )
          call Allocate_test ( m, n, n, 'M', ModuleName )

          ! Loop over the instances
          do i = 1, qt%noInstances
            if ( .not. qt%coherent ) surfs => qt%surfs(:,i)

            ! Clear the working matrix and load the diagonal
            m = 0.0_rm
            do j = 1, n
              m(j,j) = d%values(j,i) ** 2.0
            end do

            ! Now if appropriate add off diagonal terms.
            if ( any( qt%verticalCoordinate == (/ l_height, l_pressure, l_zeta /) ) ) then
              ! Loop over off diagonal terms
              do j = 1, n
                do k = j+1, n
                  meanLength = sqrt ( l%values(j,i) * l%values(k,i) )
                  meanDiag = sqrt ( m(j,j) * m(k,k) ) 
                  if ( fraction /= 0) thisFraction = f%values(j,i)
                  select case (qt%verticalCoordinate)
                  case ( l_height )
                    distance = abs ( surfs ( (j-1)/qt%noChans + 1 ) - &
                      & surfs ( (k-1)/qt%noChans + 1) )
                  case ( l_zeta )
                    distance = abs ( surfs ( (j-1)/qt%noChans + 1 ) - &
                      & surfs ( (k-1)/qt%noChans + 1 ) ) * decade
                  case ( l_pressure )
                    distance = abs ( -log10 ( surfs( (j-1)/qt%noChans + 1) ) + &
                      &               log10 ( surfs( (k-1)/qt%noChans + 1) ) ) / decade
                  end select
                  if ( meanLength > 0.0 ) &
                    & m(j,k) = meanDiag*thisFraction*exp(-distance/meanLength)
                end do                    ! Loop over k (in M)
              end do                      ! Loop over j (in M)
            end if                        ! An appropriate vertical coordinate

            ! Now we may need to invert this, if so we need to be clever.
            if ( invert ) then
              call Allocate_test ( condition, n, 'condition', ModuleName )
              condition = d%values(:,i) <= 0.0_rv
              do j = 1, n
                if ( condition(j) ) M(j,j) = 1.0_rm
              end do
              call MatrixInversion(M, upper=.true.)
              do j = 1, n
                if ( condition(j) ) M(j,j) = 0.0_rm
              end do
              call Deallocate_test ( condition, 'condition', ModuleName )
            end if

            b = FindBlock ( covariance%m%col, q, i )
            call Sparsify ( M, covariance%m%block(b,b) )
          end do                          ! Loop over instances
          call Deallocate_test ( m, 'M', ModuleName )
        end do                            ! Loop over quantities
      end if                              ! A non diagonal fill

      call DestroyVectorInfo ( DMasked )
      call DestroyVectorInfo ( LMasked )

    end subroutine FillCovariance

    !=============================== FillVectorQuantityFromGrid ============
    subroutine FillVectorQuantityFromGrid(quantity,grid, errorCode)
      ! Dummy arguments
      type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      type (GriddedData_T), intent(in) :: GRID ! Grid to fill it from
      integer, intent(out) :: ERRORCODE   ! Error code (one of constants defined above)

      ! Local variables
      integer :: instance,surf            ! Loop counter
      integer :: instIndex,surfIndex      ! Indices

      ! Executable code
      errorCode = 0

      if (quantity%template%verticalCoordinate /= l_zeta) then
        errorCode=NotZetaForGrid
        return
      end if

      instIndex=1
      surfIndex=1

      do instance = 1, quantity%template%noInstances
        if (.not. quantity%template%stacked) instIndex=instance

        do surf = 1, quantity%template%noSurfs
          if (.not. quantity%template%coherent) surfIndex=surf
          call l3ascii_interp_field(grid, quantity%values(surf,instance), &
            & pressure=10.0**(-quantity%template%surfs(surf,instIndex)), &
            & lat=quantity%template%geodLat(surfIndex,instance), &
            & lon=quantity%template%lon(surfIndex,instance), &
            & lst=quantity%template%solarTime(surfIndex,instance), &
            & sza=quantity%template%solarZenith(surfIndex,instance), &
            & date=quantity%template%time(surfIndex,instance))
        end do                            ! End surface loop
      end do                              ! End instance loop
    end subroutine FillVectorQuantityFromGrid

    !=============================== FillVectorQuantityFromL2GP ==========
    subroutine FillVectorQuantityFromL2GP ( quantity,l2gp, interpolate, profile, &
      & errorCode )
      use MLSNumerics, only: COEFFICIENTS_R8, INTERPOLATEARRAYSETUP, &
        & INTERPOLATEARRAYTEARDOWN

      ! If the times, pressures, and geolocations match, fill the quantity with
      ! the appropriate subset of profiles from the l2gp

      ! Dummy arguments
      type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      type (L2GPData_T), intent(in) :: L2GP ! L2GP to fill from
      logical, intent(in) :: interpolate  ! Flag
      integer, intent(in) :: profile    ! Single profile to use or -1 for default
      integer, intent(out) :: errorCode ! Error code

      ! Local parameters
      real(r8), parameter :: FTOL = 1.0e-3 ! 1 kHz
      real(r8), parameter :: TOLERANCE=0.05 ! Tolerence for angles
      real(r8), parameter :: TIMETOL=5.0  ! Tolerence for time (not sure why this
      ! needs to be so big !????????? NJL)

      ! Local variables
      integer ::    FIRSTPROFILE, LASTPROFILE
      integer, dimension(1) :: FIRSTPROFILEASARRAY
      integer :: INSTANCE               ! Loop counter 
      integer :: THISPROFILE            ! Index
      type (Coefficients_R8) :: COEFFS  ! For interpolation
      real (r8), dimension(quantity%template%noSurfs) :: outZeta

      errorCode=0
      ! Make sure this quantity is appropriate
      if (.not. ValidateVectorQuantity(quantity, coherent=.TRUE., stacked=.TRUE., &
        & verticalCoordinate= (/ l_pressure, l_zeta /) ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Quantity to fill is not on pressure or zeta coordinates' )

      if ( (quantity%template%noChans/=l2gp%nFreqs) .and. &
        &  ((quantity%template%noChans/=1) .or. (l2gp%nFreqs/=0)) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Quantity and l2gp have different number of channels' )
      if ( associated ( quantity%template%frequencies ) ) then
        if ( any ( abs ( l2gp%frequency - &
          & quantity%template%frequencies ) > fTol ) ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Quantity and l2pg have different frequency grids' )
      end if

      if ( quantity%template%noSurfs /= l2gp%nLevels .and. (.not. interpolate) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Quantity and l2gp have different number of surfaces (set interpolate?)' )

      if (.not. interpolate) then 
        if ( quantity%template%verticalCoordinate == l_pressure ) then
          if ( any(ABS(-LOG10(quantity%template%surfs(:,1))+ &
            & LOG10(l2gp%pressures)) > TOLERANCE) ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Quantity and l2gp are on different surfaces (set interpolate?)' )
        else                                ! Must be l_zeta
          if ( any(ABS(quantity%template%surfs(:,1)+ &
            & LOG10(l2gp%pressures)) > TOLERANCE) ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Quantity and l2gp are on different surfaces (set interpolate?)' )
        end if
      end if

      ! Skip the position checks if we're forcing in a particular profile.
      if ( profile == -1 ) then
        ! Attempt to match up the first location
        firstProfileAsArray=minloc(abs(quantity%template%phi(1,1)-l2gp%geodAngle))
        firstProfile=firstProfileAsArray(1)
        
        ! Well, the last profile has to be noInstances later, check this would be OK
        lastProfile=firstProfile+quantity%template%noInstances-1
        if (lastProfile > l2gp%nTimes ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Quantity has profiles beyond the end of the l2gp' )

        ! Now check that geodAngle's are a sufficient match
        if (any(abs(l2gp%geodAngle(firstProfile:lastProfile)-&
          &         quantity%template%phi(1,:)) > tolerance) ) then
          call dump ( l2gp%geodAngle(firstProfile:lastProfile), 'L2GP geodetic angle' )
          call dump ( quantity%template%phi(1,:), 'Quantity Geodetic angle' )
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Quantity has profiles that mismatch l2gp in geodetic angle' )
        end if

        ! Now check that the times match
        if (any(abs(l2gp%time(firstProfile:lastProfile)- &
          &         quantity%template%time(1,:)) > timeTol) ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Quantity has profiles that mismatch l2gp in time' )
        
        ! Currently the code cannot interpolate in 3 dimensions, wouldn't
        ! be hard to code up, but no need as yet.
        if (interpolate .and. quantity%template%noChans /= 1) then
          errorCode=cantInterpolate3D
          return
        end if
      else
        ! Given a specific profile, check it's legal
        if ( profile == 0 .or. profile > l2gp%nTimes ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Illegal profile request in l2gp fill' )
      end if

      ! OK, now do the filling, it's easier if we don't have to interpolate
      if (.not. interpolate) then
        if ( profile == -1 ) then
          ! Not forcing a particular profile to all instances
          quantity%values=reshape(l2gp%l2gpValue(:,:,firstProfile:lastProfile),&
            & (/quantity%template%noChans*quantity%template%noSurfs,&
            &   quantity%template%noInstances/))
        else
          ! Spread one profile onto all instances
          quantity%values = spread ( reshape ( l2gp%l2gpValue(:,:,profile), &
            &   (/ quantity%template%noChans*quantity%template%noSurfs/) ), &
            & 2, quantity%template%noInstances )
        end if
      else
        if ( quantity%template%verticalCoordinate == l_pressure ) then
          outZeta = -log10 ( quantity%template%surfs(:,1) )
        else
          outZeta = quantity%template%surfs(:,1)
        end if
        ! Setup the interpolation we'll be doing.
        call InterpolateArraySetup ( -log10(real(l2gp%pressures, r8)), &
          & outZeta, 'Linear', coeffs, extrapolate='Clamp' )
        do instance = 1, quantity%template%noInstances
          if ( profile == -1 ) then
            thisProfile = firstProfile + instance - 1
          else
            thisProfile = profile
          end if
          ! Guess I don't really need the loop here, but it
          ! does make the spread stuff much easier, and the setup/teardown
          ! at least makes it more efficient.
          call InterpolateValues ( coeffs, &
            & -log10(real(l2gp%pressures, r8)), &  ! Old X
            & real(l2gp%l2gpValue(1,:,thisProfile), r8), & ! OldY
            & outZeta, & ! New X
            & quantity%values(:,instance), & ! New Y
            & method='Linear', extrapolate='Clamp' )
        end do
        call InterpolateArrayTeardown ( coeffs )
      end if

    end subroutine FillVectorQuantityFromL2GP

    ! -------------------------------------- FillVectorQuantityFromProfile --
    subroutine FillVectorQtyFromProfile ( key, quantity, valuesNode, &
      & instancesNode, globalUnit, dontMask )
      use MLSNumerics, only: HUNT
      ! This fill is slightly complicated.  Given a values array along
      ! the lines of [ 1000mb : 1.0K, 100mb : 1.0K,  10mb : 2.0K] etc. it
      ! does the linear interpolation appropriate to perform the fill.
      integer, intent(in) :: KEY          ! Tree node
      type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      integer, intent(in) :: VALUESNODE   ! Tree node for values
      integer, intent(in) :: INSTANCESNODE ! Tree node for instances
      integer, intent(in) :: GLOBALUNIT   ! Possible global unit
      logical, intent(in) :: DONTMASK     ! If set don't follow the fill mask

      ! Local variables
      integer :: C                      ! Channel loop counter
      integer :: HEIGHTUNIT             ! Unit for height
      integer :: NOPOINTS               ! Number of points supplied
      integer :: I,J                    ! Loop counters / indices
      integer :: S                      ! Surface loop counter
      integer :: STATUS                 ! Flag
      integer :: TESTUNIT               ! Unit for value
      logical :: LOCALOUTHEIGHTS ! Set if out heights is our own variable
      real (r8), dimension(:), pointer :: HEIGHTS ! Heights for the points
      real (r8), dimension(:), pointer :: VALUES ! Values for the points
      real (r8), dimension(:), pointer :: OUTHEIGHTS ! Heights for output
      real (r8), dimension(:), pointer :: OUTVALUES ! Single profile for output
      logical, dimension(:), pointer :: INSTANCES ! Flags
      real (r8), dimension(2) :: EXPRVALUE ! Value of expression
      integer, dimension(2) :: EXPRUNIT   ! Unit for expression
      integer, dimension(:), pointer :: ININDS ! Indices

      ! Executable code

      ! Check the quantity is amenable to this type of fill
      if ( .not. ValidateVectorQuantity ( quantity, &
        & coherent=.true. ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'The quantity is not amenable to a profile fill' )

      ! Check the units
      testUnit = quantity%template%unit
      if ( globalUnit /= phyq_Invalid ) testUnit = globalUnit

      ! Set some stuff up
      noPoints = nsons ( valuesNode ) - 1
      nullify ( heights, values, outHeights, outValues, instances )
      call Allocate_test ( heights, noPoints, 'heights', ModuleName )
      call Allocate_test ( values, noPoints, 'values', ModuleName )
      call Allocate_test ( outValues, quantity%template%noSurfs, &
        & 'outValues', ModuleName )
      call Allocate_test ( instances, quantity%template%noInstances, &
        & 'instances', ModuleName )

      ! Loop over the values
      do i = 1, noPoints
        ! Get value from tree
        call expr ( subtree ( i+1, valuesNode ), exprUnit, exprValue )
        ! Check height unit OK
        heightUnit = GetUnitForVerticalCoordinate ( quantity%template%verticalCoordinate )
        if ( exprUnit(1) /= heightUnit .and. exprUnit(1) /= PHYQ_Dimensionless &
          & .and. .not. ( exprUnit(1) == PHYQ_Pressure .and. heightUnit == PHYQ_Zeta ) ) &
          & call Announce_error ( valuesNode, 0, 'Bad height units for profile fill' )
        ! Store height
        if ( heightUnit == PHYQ_Zeta ) then
          ! Assume zeta coordinates are expressed in mb
          heights(i) = -log10 ( exprValue(1) )
        else
          heights(i) = exprValue(1)
        end if
        ! Check value unit OK
        if ( all ( exprUnit(2) /= (/ testUnit, PHYQ_Dimensionless /) ) ) &
          & call Announce_error ( valuesNode, 0, 'Bad units for profile fill' )
        ! Store value
        values ( i ) = exprValue(2)
      end do

      ! Get the appropriate height coordinate for output, for pressure take log.
      if ( quantity%template%verticalCoordinate == l_pressure ) then
        localOutHeights = .true.
        call Allocate_test ( outHeights, quantity%template%noSurfs, &
          & 'outHeights', ModuleName )
        outHeights = -log10 ( quantity%template%surfs(:,1) )
      else
        localOutHeights = .false.
        outHeights => quantity%template%surfs(:,1)
      end if

      ! Now, if the quantity is coherent, let's assume the user wanted the
      ! 'nearest' values
      if ( quantity%template%coherent ) then
        nullify ( inInds )
        call allocate_test ( inInds, noPoints, 'inInds', ModuleName )
        call hunt ( outHeights, heights, inInds, &
          & nearest=.true., allowTopValue=.true. )
        heights = outHeights ( inInds )
        call deallocate_test ( inInds, 'inInds', ModuleName )
      end if

      ! Now do the interpolation for the first instance
      call InterpolateValues ( heights, values, outHeights, &
        & outValues, 'Linear', extrapolate='Constant' )

      ! Work out what instances we're going to spread this to
      instances = .true.
      if ( instancesNode /= 0 ) then
        call GetIndexFlagsFromList ( instancesNode, instances, &
          & status, noError = .true. )
        if ( status /= 0 ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to parse instances (wrong units?)' )
      end if

      ! Spread it into the other instances, worry about the mask and instances
      do i = 1, quantity%template%noInstances
        if ( instances(i) ) then
          j = 1
          do s = 1, quantity%template%noSurfs
            do c = 1, quantity%template%noChans
              if ( associated(quantity%mask) ) then
                if ( iand ( ichar(quantity%mask(j,i)), m_fill ) == 0 ) &
                  & quantity%values(j,i) = outValues(s)
              else
                quantity%values(j,i) = outValues(s)
              end if
              j = j + 1
            end do
          end do
        end if
      end do

      ! Finish off
      if ( localOutHeights ) call Deallocate_test ( outHeights, &
        & 'outHeights', ModuleName )
      call Deallocate_test ( heights, 'heights', ModuleName )
      call Deallocate_test ( values, 'values', ModuleName )
      call Deallocate_test ( outValues, 'outValues', ModuleName )
      call Deallocate_test ( instances, 'instances', ModuleName )
    end subroutine FillVectorQtyFromProfile

    ! ------------------------------------------- FillLOSVelocity ---
    subroutine FillLOSVelocity ( key, qty, tngtECI, scECI, scVel)
      ! A special fill from geometry arguments
      use Geometry, only: OMEGA => W
      integer, intent(in) :: KEY
      type (VectorValue_T), intent(inout) :: QTY
      type (VectorValue_T), intent(in) :: TNGTECI
      type (VectorValue_T), intent(in) :: SCECI
      type (VectorValue_T), intent(in) :: SCVEL


      ! Local variables
      integer :: MAF                      ! Loop counter
      integer :: MIF                      ! Loop counter
      integer :: noMAFs                   ! Number of major frames
      integer :: noMIFs                   ! Number of minor frames for this module
      integer :: x,y,z                    ! Indicies into the vectors

      real (r8), dimension(3) :: tngtVel   ! Due to rotation of earth
      real (r8), dimension(3) :: los       ! Normalised line of sight vector

      ! Executable code
      ! First check that things are OK.
      if ( .not. ValidateVectorQuantity ( qty, &
        & quantityType=(/l_losVel/), &
        & minorFrame=.true., &
        & frequencyCoordinate=(/l_none/) ) ) call Announce_Error ( key, No_Error_Code, &
        & 'Quantity to fill is not a valid LOS Velocity' )
      if ( .not. ValidateVectorQuantity ( tngtECI, &
        & quantityType=(/l_tngtECI/), &
        & minorFrame=.true., &
        & frequencyCoordinate=(/l_xyz/) ) ) call Announce_Error ( key, No_Error_Code, &
        & 'Tangent ECI quantity is not of an appropriate form' )
      if ( .not. ValidateVectorQuantity ( scECI, &
        & quantityType=(/l_scECI/), &
        & minorFrame=.true., &
        & frequencyCoordinate=(/l_xyz/) ) ) call Announce_Error ( key, No_Error_Code, &
        & 'Spacecraft ECI quantity is not of an approriate form' )
      if ( qty%template%instrumentModule /= tngtECI%template%instrumentModule ) &
        & call Announce_Error ( key, No_Error_Code, &
        & 'LOS Velocity and Tangent ECI quantities are not for the same module' )
      if ( .not. IsModuleSpacecraft ( scECI%template%instrumentModule ) ) &
        & call Announce_Error ( key, No_Error_Code, &
        & 'Spacecraft ECI quantity is not for the spacecraft' )

      noMAFs = qty%template%noInstances
      noMIFs = qty%template%noSurfs

      do maf = 1, noMAFs
        do mif = 1, noMIFs

          ! First compute the tangent point velocity in ECI coordinates due 
          ! to the rotation of the earth.  This no doubt makes approximations
          ! due to the slight non alignment between the earth's rotation axis and
          ! the ECI z axis, but I'm going to ignore this.

          ! Work out the indices in 3*mif,maf space
          x = 1 + 3*(mif-1)
          y = x+1
          z = x+2

          tngtVel= omega* (/ -tngtECI%values(y,maf), &
            &                 tngtECI%values(x,maf), 0.0_r8 /)

          ! Now compute the line of sight direction normal
          los = tngtECI%values(x:z,maf) - scECI%values(x:z,maf)
          los = los / sqrt(sum(los**2))

          ! Now compute the net velocity in this direction.  For the moment I'll
          ! assume +ve means the sc and tp are moving apart, and -ve that they're
          ! getting closer.

          qty%values(mif,maf) = dot_product(tngtVel, los) - &
            &                   dot_product(scVel%values(x:z,maf), los)

          ! Note that even though x,y,z have been used up to now for a GHz/THz
          ! minor frame quantity, they're OK with this sc one too.
        end do
      end do
    end subroutine FillLOSVelocity

    ! ------------------------------------------- FillableChiSq ---
    function FillableChiSq ( qty, measQty, modelQty, noiseQty ) result ( aok )
      ! Purpose (A)
      !   Check whether we may proceed with special fill of chi squared
      !   case where all VectorValue_T args present
      ! Purpose (B)
      !   Check whether we may proceed with special fill of addNoise
      !   case where missing noiseQty arg
      type (VectorValue_T), intent(in) ::              QTY
      type (VectorValue_T), intent(in) ::              modelQty
      type (VectorValue_T), intent(in) ::              measQty
      type (VectorValue_T), optional, intent(in) ::    noiseQty
      LOGICAL ::                                       AOK

      ! What we will check is that (for the args we have been given):
      ! (0) all quantities have associated values
      ! (1) all quantities have same molecule
      ! (2) all quantities have same signal
      ! (3) all quantities have same HGrid
      ! (4A) all but qty (chiSq) have same VGrid
      ! (4B) all have same VGrid
      ! Ah, but radiances have no VGrid or HGrid,
      ! and molecule should be part of signal, so
      ! if radiances need only check on signal
      ! if vmr, check on others

     ! Local variables
      LOGICAL ::       minorFrame   ! TRUE if radiances, FALSE if vmr

      aok = .true.

      ! (0)
      if ( present(noiseQty) ) then
        aok = associated(noiseQty%values)
        if ( DEEBUG .and. .not. associated(noiseQty%values) ) &
          & call announce_error( 0, No_Error_code, &
          & 'Noise values unassociated in FillableChiSq')
      end if
      aok = aok .and. &
        & associated(qty%values) .and. &
        & associated(measQty%values) .and. &
        & associated(modelQty%values)

      if ( DEEBUG ) then
        if ( .not. associated(qty%values) ) &
          & call announce_error( 0, No_Error_code, &
          & 'Quantity values unassociated in FillableChiSq')
        if ( .not. associated(measQty%values) ) &
          & call announce_error( 0, No_Error_code, &
          & 'Measurements values unassociated in FillableChiSq')
        if ( .not. associated(modelQty%values) ) &
          & call announce_error( 0, No_Error_code, &
          & 'Model values unassociated in FillableChiSq')
      end if

      if ( .not. aok ) return

      minorFrame = qty%template%minorFrame .or. qty%template%majorFrame
      ! (1)
      if (.not. minorFrame ) then
        aok = aok .and. &
          & (qty%template%molecule == measQty%template%molecule) &
          & .and. &
          & (qty%template%molecule == modelQty%template%molecule)
        if ( present(noiseQty) ) aok = aok &
          & .and. &
          & (qty%template%molecule == noiseQty%template%molecule)
      end if

      ! (2)
      if (minorFrame ) then
        aok = aok .and. &
          & (qty%template%signal == measQty%template%signal) &
          & .and. &
          & (qty%template%signal == modelQty%template%signal)
        if ( present(noiseQty) ) aok = aok &
          & .and. &
          & (qty%template%signal == noiseQty%template%signal)
      end if

      ! (3)
      if (.not. minorFrame ) then
        aok = aok .and. &
          & DoHgridsMatch( qty, measQty ) &
          & .and. &
          & DoHgridsMatch( qty, modelQty )
        if ( present(noiseQty) ) aok = aok &
          & .and. &
          & DoHgridsMatch( qty, noiseQty )

        ! (4)
        aok = aok .and. &
          & DoVgridsMatch( measqty, modelQty )
        if ( present(noiseQty) ) then
          aok = aok &
          & .and. &
          & DoVgridsMatch( measqty, noiseQty )
        else
          aok = aok &
          & .and. &
          & DoVgridsMatch( measqty, Qty )
        end if
      end if

      return
    end function FillableChiSq

    ! ------------------------------------------- FillChiSqChan ---
    subroutine FillChiSqChan ( key, qty, measQty, modelQty, noiseQty, &
    & dontMask, ignoreZero, ignoreNegative, multiplier, &
    & firstInstance, lastInstance )
      ! A special fill of chi squared 
      ! broken out according to channels
      ! Formal arguments
      integer, intent(in) :: KEY
      type (VectorValue_T), intent(inout) :: QTY
      type (VectorValue_T), intent(in) ::    modelQty
      type (VectorValue_T), intent(in) ::    measQty
      type (VectorValue_T), intent(in) ::    noiseQty
      logical, intent(in)           ::       dontMask    ! Use even masked values
      logical, intent(in)           ::       ignoreZero  ! Ignore 0 values of noiseQty
      logical, intent(in)           ::       ignoreNegative  ! Ignore <0 values of noiseQty
      real, dimension(:), intent(in), optional :: multiplier

      integer, intent(in), optional ::       firstInstance, lastInstance
      ! The last two are set if only part (e.g. overlap regions) of the quantity
      ! is to be stored in qty

      ! Local variables
      real(r8), dimension(:), pointer  ::    VALUES => NULL()
      integer ::                             UseFirstInstance, UseLastInstance, &
      &                                      NoOutputInstances
      integer ::                             C           ! Channel loop counter
      integer ::                             S           ! Surface loop counter
      integer ::                             I           ! Instances
      integer ::                             QINDEX
      integer ::                             NOCHANS
      integer ::                             N           ! Num. of summed values
      logical ::                             skipMe
      real                             ::    a, b

      ! Executable code

     ! Either multiplier = [a, b] or multiplier = 1/a if a=b are possible
      if ( .not. present(multiplier) ) then
        a = 1.
        b = 1.
      else if (&
      & multiplier(1) == UNDEFINED_VALUE .and. multiplier(2) == UNDEFINED_VALUE &
      & ) then
        a = 1.
        b = 1.
      else if ( multiplier(2) == UNDEFINED_VALUE ) then
        a = 0.
        if ( multiplier(1) /= 0.0 ) a = 1.0/multiplier(1)
        b = a
      else
        a = multiplier(1)
        b = multiplier(2)
      end if

      ! First check that things are OK.
      if (.not. ValidateVectorQuantity ( qty, &
        & quantityType=(/l_chiSqChan/), majorFrame=.true.) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Attempting to fill wrong quantity with chi^2 channelwise'  )
        if ( DEEBUG ) then
          call output('major frame? ', advance = 'no')
          call output(qty%template%majorFrame, advance = 'no')
          call output('   quantity type ', advance = 'no')
          call output(qty%template%quantityType, advance = 'no')
          call output('   compared with ', advance = 'no')
          call output(l_chiSqChan, advance = 'yes')
        end if
        return
      else if (.not. FillableChiSq ( qty, measQty, modelQty, noiseQty ) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Incompatibility among vector quantities filling chi^2 channelwise'  )
        return
      else if (any ( noiseQty%values == 0.0) .and. &
        & .not. (ignoreZero .or. .not. dontMask) ) then
        call Announce_error ( key, No_Error_code, &
        & 'A vanishing error filling chi^2 channelwise'  )
        return
      end if

      ! Work out what to do with the first and last Instance information

      if ( PRESENT(firstInstance) ) then
        useFirstInstance = firstInstance
      else
        useFirstInstance = 1
      end if

      if ( PRESENT(lastInstance) ) then
        useLastInstance = lastInstance
      else
        useLastInstance = qty%template%noInstances
      end if
      noOutputInstances = useLastInstance-useFirstInstance+1
      ! If we've not been asked to output anything then don't carry on
      if ( noOutputInstances < 1 ) return

      call allocate_test(values, measQty%template%noSurfs, &
        & 'chi^2 unsummed', ModuleName)
      noChans = qty%template%noChans
      do i=useFirstInstance, useLastInstance
        do c=1, noChans
          N = 0
          values = 0.0
          do s=1, measQty%template%noSurfs
            qIndex = c + (s-1)*nochans
            skipMe = &
            & .not. dontMask .and. ( &
            &   isVectorQtyMasked(measQty, qIndex, i) .or. &
            &   isVectorQtyMasked(modelQty, qIndex, i) .or. &
            &   isVectorQtyMasked(noiseQty, qIndex, i) ) &
            & .or. (ignoreNegative .and. noiseQty%values(qIndex, i) < 0.0 ) &
            & .or. (ignoreZero .and. noiseQty%values(qIndex, i) == 0.0 )
            if ( .not. skipMe ) then
              values(s) = ( &
              & (a*measQty%values(qIndex, i) - b*modelQty%values(qIndex, i)) &
              & / &
              & noiseQty%values(qIndex, i) &
              &  ) ** 2
              N = N + 1
            end if
          end do
          if ( N > 0 ) then
            qty%values(c, i) = sum(values) / N
          else
            qty%values(c, i) = 0.
          end if
        end do
      end do
      call deallocate_test(values, &
        & 'chi^2 unsummed', ModuleName)
    end subroutine FillChiSqChan

    ! ------------------------------------------- FillChiSqMMaf ---
    subroutine FillChiSqMMaf ( key, qty, measQty, modelQty, noiseQty, &
    & dontMask, ignoreZero, ignoreNegative, multiplier, &
    & firstInstance, lastInstance )
      ! A special fill of chi squared 
      ! broken out according to major frames
      ! Formal arguments
      integer, intent(in) :: KEY
      type (VectorValue_T), intent(inout) :: QTY
      type (VectorValue_T), intent(in) ::    modelQty
      type (VectorValue_T), intent(in) ::    measQty
      type (VectorValue_T), intent(in) ::    noiseQty
      logical, intent(in)           ::       dontMask    ! Use even masked values
      logical, intent(in)           ::       ignoreZero  ! Ignore 0 values of noiseQty
      logical, intent(in)           ::       ignoreNegative  ! Ignore <0 values of noiseQty
      real, dimension(:), intent(in), optional :: multiplier

      integer, intent(in), optional ::       firstInstance, lastInstance
      ! The last two are set if only part (e.g. overlap regions) of the quantity
      ! is to be stored in the qty

      ! Local variables
      real(r8), dimension(:), pointer  ::    VALUES => NULL()
      integer ::                             UseFirstInstance, UseLastInstance, &
      &                                      NoOutputInstances
      integer ::                             I           ! Instances
      integer ::                             ROW         ! Running 1st coord
      integer ::                             INSTANCELEN ! Num of rows
      integer ::                             N           ! Num. of summed values
      logical ::                             skipMe

      ! Executable code
      real                             ::    a, b

      ! Executable code

     ! Either multiplier = [a, b] or multiplier = 1/a if a=b are possible
      if ( .not. present(multiplier) ) then
        a = 1.
        b = 1.
      else if (&
      & multiplier(1) == UNDEFINED_VALUE .and. multiplier(2) == UNDEFINED_VALUE &
      & ) then
        a = 1.
        b = 1.
      else if ( multiplier(2) == UNDEFINED_VALUE ) then
        a = 0.
        if ( multiplier(1) /= 0.0 ) a = 1.0/multiplier(1)
        b = a
      else
        a = multiplier(1)
        b = multiplier(2)
      end if

      ! First check that things are OK.
      if (.not. ValidateVectorQuantity ( qty, &
        & quantityType=(/l_chiSqMMaf/), majorFrame=.true.) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Attempting to fill wrong quantity with chi^2 MMAFwise'  )
        if ( DEEBUG ) then
          call output('major frame? ', advance = 'no')
          call output(qty%template%majorFrame, advance = 'no')
          call output('   quantity type ', advance = 'no')
          call output(qty%template%quantityType, advance = 'no')
          call output('   compared with ', advance = 'no')
          call output(l_chiSqMMaf, advance = 'yes')
        end if
        return
      else if (.not. FillableChiSq ( qty, measQty, modelQty, noiseQty ) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Incompatibility among vector quantities filling chi^2 MMAFwise'  )
        return
      else if (any ( noiseQty%values == 0.0) .and. &
        & .not. (ignoreZero .or. .not. dontMask) ) then
        call Announce_error ( key, No_Error_code, &
        & 'A vanishing noise filling chi^2 MMAFwise'  )
        return
      end if

      ! Work out what to do with the first and last Instance information

      if ( PRESENT(firstInstance) ) then
        useFirstInstance = firstInstance
      else
        useFirstInstance = 1
      end if

      if ( PRESENT(lastInstance) ) then
        useLastInstance = lastInstance
      else
        useLastInstance = qty%template%noInstances
      end if
      noOutputInstances = useLastInstance-useFirstInstance+1
      ! If we've not been asked to output anything then don't carry on
      if ( noOutputInstances < 1 ) return

      instanceLen = measQty%template%noChans * measQty%template%noSurfs
      call allocate_test(values, instanceLen, &
        & 'chi^2 unsummed', ModuleName)
      do i=useFirstInstance, useLastInstance
        if( .not. (.not. dontMask .or. ignoreNegative .or. ignoreZero )) then
            values = ( &
            & (measQty%values(:, i) - modelQty%values(:, i)) &
            & / &
            & noiseQty%values(:, i) &
            &  ) ** 2
            qty%values(1, i) = sum(values) / instanceLen
        else
          N = 0
          values = 0.0
          do row = 1, instanceLen
            skipMe = &
            & .not. dontMask .and. ( &
            &   isVectorQtyMasked(measQty, row, i) .or. &
            &   isVectorQtyMasked(modelQty, row, i) .or. &
            &   isVectorQtyMasked(noiseQty, row, i) ) &
            & .or. (ignoreNegative .and. noiseQty%values(row, i) < 0.0 ) &
            & .or. (ignoreZero .and. noiseQty%values(row, i) == 0.0 )
            if ( .not. skipMe ) then
              values(row) = ( &
              & (a*measQty%values(row, i) - b*modelQty%values(row, i)) &
              & / &
              & noiseQty%values(row, i) &
              &  ) ** 2
              N = N + 1
            end if
          end do
          if ( N > 0 ) then
            qty%values(1, i) = sum(values) / N
          else
            qty%values(1, i) = 0.
          end if
        end if
      end do
      call deallocate_test(values, &
        & 'chi^2 unsummed', ModuleName)
    end subroutine FillChiSqMMaf

    ! ------------------------------------------- FillChiSqMMif ---
    subroutine FillChiSqMMif ( key, qty, measQty, modelQty, noiseQty, &
    & dontMask, ignoreZero, ignoreNegative, multiplier, &
    & firstInstance, lastInstance )
      ! A special fill of chi squared 
      ! broken out according to Mifs
      ! Formal arguments
      integer, intent(in) :: KEY
      type (VectorValue_T), intent(inout) :: QTY
      type (VectorValue_T), intent(in) ::    modelQty
      type (VectorValue_T), intent(in) ::    measQty
      type (VectorValue_T), intent(in) ::    noiseQty
      logical, intent(in)           ::       dontMask    ! Use even masked values
      logical, intent(in)           ::       ignoreZero  ! Ignore 0 values of noiseQty
      logical, intent(in)           ::       ignoreNegative  ! Ignore <0 values of noiseQty
      real, dimension(:), intent(in), optional :: multiplier

      integer, intent(in), optional ::       firstInstance, lastInstance
      ! The last two are set if only part (e.g. overlap regions) of the quantity
      ! is to be stored in the qty

      ! Local variables
      real(r8), dimension(:), pointer  ::    VALUES => NULL()
      integer ::                             UseFirstInstance, UseLastInstance, &
      &                                      NoOutputInstances
      integer ::                             C           ! Channel loop counter
      integer ::                             S           ! Surface loop counter
      integer ::                             I           ! Instances
      integer ::                             QINDEX
      integer ::                             NOMIFS
      integer ::                             N           ! Num. of summed values
      logical ::                             skipMe

      ! Executable code
      real                             ::    a, b

      ! Executable code

     ! Either multiplier = [a, b] or multiplier = 1/a if a=b are possible
      if ( .not. present(multiplier) ) then
        a = 1.
        b = 1.
      else if (&
      & multiplier(1) == UNDEFINED_VALUE .and. multiplier(2) == UNDEFINED_VALUE &
      & ) then
        a = 1.
        b = 1.
      else if ( multiplier(2) == UNDEFINED_VALUE ) then
        a = 0.
        if ( multiplier(1) /= 0.0 ) a = 1.0/multiplier(1)
        b = a
      else
        a = multiplier(1)
        b = multiplier(2)
      end if

      ! First check that things are OK.
      if (.not. ValidateVectorQuantity ( qty, &
        & quantityType=(/l_chiSqMMif/), minorFrame=.true.) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Attempting to fill wrong quantity with chi^2 MMIFwise'  )
        if ( DEEBUG ) then
          call output('minor frame? ', advance = 'no')
          call output(qty%template%minorFrame, advance = 'no')
          call output('   quantity type ', advance = 'no')
          call output(qty%template%quantityType, advance = 'no')
          call output('   compared with ', advance = 'no')
          call output(l_chiSqMMif, advance = 'yes')
        end if
        return
      else if (.not. FillableChiSq ( qty, measQty, modelQty, noiseQty ) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Incompatibility among vector quantities filling chi^2 MMIFwise'  )
        return
      else if (any ( noiseQty%values == 0.0) .and. &
        & .not. (ignoreZero .or. .not. dontMask) ) then
        call Announce_error ( key, No_Error_code, &
        & 'A vanishing noise filling chi^2 MMIFwise'  )
        return
      end if

      ! Work out what to do with the first and last Instance information

      if ( PRESENT(firstInstance) ) then
        useFirstInstance = firstInstance
      else
        useFirstInstance = 1
      end if

      if ( PRESENT(lastInstance) ) then
        useLastInstance = lastInstance
      else
        useLastInstance = qty%template%noInstances
      end if
      noOutputInstances = useLastInstance-useFirstInstance+1
      ! If we've not been asked to output anything then don't carry on
      if ( noOutputInstances < 1 ) return

      call allocate_test(values, measQty%template%noChans, &
        & 'chi^2 unsummed', ModuleName)
      noMIFs = measQty%template%noSurfs
      do i=useFirstInstance, useLastInstance
        do s=1, noMIFs
          N = 0
          values = 0.0
          do c=1, measQty%template%noChans
            qIndex = c + (s-1)*measQty%template%noChans
            skipMe = &
            & .not. dontMask .and. ( &
            &   isVectorQtyMasked(measQty, qIndex, i) .or. &
            &   isVectorQtyMasked(modelQty, qIndex, i) .or. &
            &   isVectorQtyMasked(noiseQty, qIndex, i) ) &
            & .or. (ignoreNegative .and. noiseQty%values(qIndex, i) < 0.0 ) &
            & .or. (ignoreZero .and. noiseQty%values(qIndex, i) == 0.0 )
            if ( .not. skipMe ) then
              values(c) = ( &
              & (a*measQty%values(qIndex, i) - b*modelQty%values(qIndex, i)) &
              & / &
              & noiseQty%values(qIndex, i) &
              &  ) ** 2
              N = N + 1
            end if
          end do
          if ( N > 0 ) then
            qty%values(s, i) = sum(values) / N
          else
            qty%values(s, i) = 0.
          end if
        end do
      end do
      call deallocate_test(values, &
        & 'chi^2 unsummed', ModuleName)
    end subroutine FillChiSqMMif

    ! ------------------------------------------- FillColAbundance ---
    subroutine FillColAbundance ( key, qty, bndPressQty, vmrQty, &
    & firstInstance, lastInstance )
      ! A special fill according to Appendix A of EOS MLS ATBD
      ! JPL D-16159
      ! EOS MLS DRL 601 (part 3)
      ! ATBD-MLS-03
      ! (Livesey and Wu)

      ! Assumptions:
      ! (1)This fill operation is triggered by a command
      !    such as the following in the lcf
      !      Fill, state.columnO3, method=special, vmrQuantity=state.o3, $
      !      boundaryPressure=state.tpPressure
      ! (2)the vmr is in units of PHYQ_VMR and not, say, ppmv;
      !    it is in fact identical to the coefficients of the mls basis functions
      ! (3)The pressure surfaces are in hPa, but not all necessarily at the
      !    same logarithmic distance from one another
      ! (4)The tropospheric boundary pressure is somewhere in between the surfs
      ! (5)Unless first,last instances are args, fill all instances
      !    (unlike join which has to worry about chunks and overlaps)
      integer, intent(in) :: KEY
      type (VectorValue_T), intent(inout) :: QTY
      type (VectorValue_T), intent(in) ::    bndPressQty
      type (VectorValue_T), intent(in) ::    vmrQty
      integer, intent(in), optional ::       firstInstance, lastInstance
      ! The last two are set if only part (e.g. overlap regions) of the quantity
      ! is to be stored in the column data.

      ! Local variables
      real (r8) :: AoverMg         ! A/(M g) from Appendix A
      logical ::   zeta_surfs      ! If true, surfs are zeta-type; else pressures
      integer ::   status
      integer ::   surface
      integer ::   instance
      integer ::   surfaceInstance
      integer ::   firstSurface
      integer ::   UseFirstInstance, UseLastInstance, &
      &            NoOutputInstances
      real (r8) :: columnSum
      real (r8) :: Delta_p_plus    ! p[j+1] - p[j]
      real (r8) :: Delta_p_minus   ! p[j-1] - p[j]
      real (r8) :: Delta_log_plus  ! ln p[j+1] - ln p[j]
      real (r8) :: Delta_log_minus ! ln p[j-1] - ln p[j]

      real (r8), allocatable, dimension(:) :: p         ! p[i] in hPa

      ! Executable code
      ! First check that things are OK.
      if ( (qty%template%quantityType /= l_columnAbundance) .or. &
        &  (bndPressQty%template%quantityType /= l_boundaryPressure) .or. &
        &  (vmrQty%template%quantityType /= l_vmr) ) then
            call Announce_error ( key, No_Error_code, &
                & 'Wrong quantity type found while filling column abundance'  )
        return
      else if ( qty%template%molecule /= vmrQty%template%molecule) then
            call Announce_error ( key, No_Error_code, &
                & 'Attempt to fill column abundance with different molecule'  )
        return
      else if ( &
      & .not. ( &
      & DoHgridsMatch( qty, vmrQty ) &
      & .and. &
      & DoHgridsMatch( qty, bndPressQty ) &
      & ) &
      & ) then
            call Announce_error ( key, No_Error_code, &
                & 'Attempt to fill column abundance with different HGrids'  )
        return
      else if ( .not. &
      & any(vmrQty%template%verticalCoordinate == (/l_pressure, l_zeta/)) &
      & ) then
            call Announce_error ( key, No_Error_code, &
                & 'Fill column abundance, but vmr not on [log]pressure surfs.'  )
        return
      end if

      ! Work out what to do with the first and last Instance information

      if ( PRESENT(firstInstance) ) then
        useFirstInstance = firstInstance
      else
        useFirstInstance = 1
      end if

      if ( PRESENT(lastInstance) ) then
        useLastInstance = lastInstance
      else
        useLastInstance = qty%template%noInstances
      end if
      noOutputInstances = useLastInstance-useFirstInstance+1
      ! If we've not been asked to output anything then don't carry on
      if ( noOutputInstances < 1 ) return
     !    AoverMg = 4.12e25 / (2.687e20 * 0.192)
     !    This assumes that
     ! (1) p is in hPa
     ! (2) f is in PHYQ_vmr (*not* ppmv)
     AoverMg = 4.12e5 / (2.687 * 0.192)

     zeta_surfs = vmrQty%template%verticalCoordinate == l_zeta
     allocate(p(vmrQty%template%noSurfs), stat=status)
     if(status /= 0) then
            call Announce_error ( key, No_Error_code, &
                & 'Error in allocating p'  )
        return
     end if   
  !   do instance=1, vmrQty%template%noInstances
     do instance=useFirstInstance, useLastInstance
        if(vmrQty%template%coherent) then
           surfaceInstance=1
        else
           surfaceInstance=instance
        end if

        if(zeta_surfs) then
           ! Invert zeta = -log10(p)
           do surface=1, vmrQty%template%noSurfs
              p(surface) = exp(-log(10.)* &
              & vmrQty%template%surfs(surface, surfaceInstance))
           end do
        else
           do surface=1, vmrQty%template%noSurfs
              p(surface) = vmrQty%template%surfs(surface, surfaceInstance)
           end do
        end if

        if(p(1) &
           &  < bndPressQty%values(1, instance)) then
            call Announce_error ( key, No_Error_code, &
                & 'Fill column abundance, but tropopause below VGrid'  )
        end if

     ! Find 1st surface at or above tropopause
     ! (i.e., at a pressure equal to or less than boundaryPressure)
        do surface=1, vmrQty%template%noSurfs
           if(p(surface) &
           &  <= bndPressQty%values(1, instance)) exit
        end do
        firstSurface = surface
        if(firstSurface > vmrQty%template%noSurfs-2) then
            call Announce_error ( key, No_Error_code, &
                & 'Fill column abundance, but tropopause above VGrid'  )
        end if
     ! Do summation
     ! Initialize sum, Deltas
        columnSum = 0.
        Delta_p_plus = p(firstSurface+1) - p(firstSurface)
        Delta_log_plus = log(p(firstSurface+1)) - log(p(firstSurface))
     ! Loop over surfaces from tropoause+1 to uppermost-1
        do surface = firstSurface+1, vmrQty%template%noSurfs-1
           Delta_p_minus = - Delta_p_plus
           Delta_log_minus = - Delta_log_plus
           Delta_p_plus = p(surface+1) - p(surface)
           Delta_log_plus = log(p(surface+1)) - log(p(surface))
           columnSum = columnSum + &
           & vmrQty%values(surface, instance)* &
           & ( &
           & Delta_p_minus/Delta_log_minus &
           & - &
           & Delta_p_plus/Delta_log_plus &
           & )          
        end do
        qty%values(1, instance) = AoverMg * columnSum
     end do

     deallocate(p, stat=status)
     if(status /= 0) then
            call Announce_error ( key, No_Error_code, &
                & 'Error in deallocating p'  )
     end if   

    end subroutine FillColAbundance

    ! ------------------------------------- FillFoldedRadiance ---
    subroutine FillFoldedRadiance ( radiance, lsb, usb, &
      & lsbFraction, usbFraction, key )
      type (VectorValue_T), intent(inout) :: RADIANCE
      type (VectorValue_T), intent(in) :: USB
      type (VectorValue_T), intent(in) :: LSB
      type (VectorValue_T), intent(in) :: USBFRACTION
      type (VectorValue_T), intent(in) :: LSBFRACTION
      integer, intent(in) :: KEY

      ! Local variables
      integer :: C                        ! Channel loop inductor
      integer :: I                        ! Array index
      integer :: MIF                      ! Minor frame loop inductor

      ! Executable code
      ! First some sanity checks
      if (.not. ValidateVectorQuantity ( radiance, quantityType=(/l_radiance/), &
        & sideband=(/0/), minorFrame=.true. )) &
        & call Announce_Error ( key, 0, 'Inappropriate radiance quantity to fill' )
      if (.not. ValidateVectorQuantity ( lsb, quantityType=(/l_radiance/), &
        & sideband=(/-1/), signal=(/radiance%template%signal/), minorFrame=.true. )) &
        & call Announce_Error ( key, 0, 'Inappropriate lsb radiance quantity for fill' )
      if (.not. ValidateVectorQuantity ( usb, quantityType=(/l_radiance/), &
        & sideband=(/1/), signal=(/radiance%template%signal/), minorFrame=.true. )) &
        & call Announce_Error ( key, 0, 'Inappropriate usb radiance quantity for fill' )
      if (.not. ValidateVectorQuantity ( lsbFraction, quantityType=(/l_sidebandRatio/), &
        & signal=(/radiance%template%signal/), sideband=(/-1/) ) ) &
        & call Announce_Error ( key, 0, 'Inappropriate lsbFraction quantity for fill' )
      if (.not. ValidateVectorQuantity ( usbFraction, quantityType=(/l_sidebandRatio/), &
        & signal=(/radiance%template%signal/), sideband=(/-1/) ) ) &
        & call Announce_Error ( key, 0, 'Inappropriate usbFraction quantity for fill' )

      ! Now do the work
      i = 1                               ! Use i as a composit mif,channel index
      do mif = 1, radiance%template%noSurfs
        do c = 1, radiance%template%noChans
          radiance%values(i,:) = &
            & lsbFraction%values(c,1) * lsb%values(i,:) + &
            & usbFraction%values(c,1) * usb%values(i,:)
          i = i + 1
        end do
      end do

    end subroutine FillFoldedRadiance

    ! ------------------------------------ FillPhiTanWithRefraction --
    subroutine FillPhiTanWithRefraction ( key, quantity, ptanQuantity, &
      & temperatureQuantity, h2oQuantity, refract )
      integer, intent(in) :: KEY          ! Tree node
      type (VectorValue_T), intent(inout) :: QUANTITY ! PhiTan quantity to fill
      type (VectorValue_T), pointer :: PTANQUANTITY ! Ptan for same module
      type (VectorValue_T), pointer :: TEMPERATUREQUANTITY ! Temperature
      type (VectorValue_T), pointer :: H2OQUANTITY ! Water vapor
      logical, intent(in) :: REFRACT      ! Do refraction or not

      ! Executable code
      ! First check sanity
      if ( .not. ValidateVectorQuantity ( quantity, &
        & quantityType=(/l_phiTan/), minorFrame=.true. ) ) &
        & call Announce_error ( key, 0, 'Quantity to fill is not phiTan' )
      if ( refract ) then
        ! More sanity checks
        if ( .not. ValidateVectorQuantity ( ptanQuantity, &
          & quantityType=(/l_ptan/), minorFrame=.true. ) ) &
          & call Announce_error ( key, 0, 'Problem with ptan quantity for phiTan fill' )
        if ( quantity%template%instrumentModule /= &
          &  ptanQuantity%template%instrumentModule ) &
          & call Announce_error ( key, 0, 'phiTan and ptan quantities not for same module' )
        if ( .not. ValidateVectorQuantity ( temperatureQuantity, &
          & quantityType=(/l_temperature/), coherent=.true., stacked=.true., &
          & frequencyCoordinate=(/l_none/), verticalCoordinate=(/l_zeta/) ) ) &
          & call Announce_error ( key, 0, 'Problem with temperature quantity for phiTan fill' )
        if ( .not. ValidateVectorQuantity ( h2oQuantity, &
          & quantityType=(/l_vmr/), molecule=(/l_h2o/), coherent=.true., stacked=.true., &
          & frequencyCoordinate=(/l_none/), verticalCoordinate=(/l_zeta/) ) ) &
          & call Announce_error ( key, 0, 'Problem with temperature quantity for phiTan fill' )

        ! OK, do the refraction calculation
        call Announce_error ( key, 0, 'Refract=true, not yet supported' )
      else
        ! Just copy it from the template
        quantity%values = quantity%template%phi
      end if

    end subroutine FillPhiTanWithRefraction

      ! ------------------------------------- FillRHIFromH2O ----
    subroutine FillRHIFromH2O ( key, quantity, &
     & sourceQuantity, temperatureQuantity, &
     & dontMask, ignoreZero, ignoreNegative, interpolate, &
     & markUndefinedValues, invert )
      ! Convert h2o vmr to %RHI for all instances, channels, surfaces
      ! (See Eq. 9 from "UARS Microwave Limb Sounder upper tropospheric
      !  humidity measurement: Method and validation" Read et. al. 
      !  J. Geoph. Res. Dec. 2001 (106) D23)

      !  Method:
      ! (1) straight convert--all quantities must have the same shape
      !     (strictly we assume they have _all_ the same geolocations)
      ! (2) interpolate--all quantities may have different shapes
      !     (the interpolation will be along the vertical coordinate only)
      !     I.e., for xQuantity (where x can be h2o or temperature)
      !     if NoChans(xQuantity) /= NoChans(Quantity)
      !        => use only xQuantity(channel==1)
      !     if NoInstances(xQuantity) /= NoInstances(Quantity)
      !        => use only xQuantity(instance==1)
      !
      ! (3) if invert is TRUE, like (1) but its inverse: %RHI to h2o vmr
      integer, intent(in) :: key          ! For messages
      ! Actually, the meaning of the next two is reversed if invert is TRUE)
      type (VectorValue_T), intent(inout) :: QUANTITY ! (rhi) Quantity to fill
      type (VectorValue_T), intent(in) :: sourceQuantity ! vmr (unless invert)
      type (VectorValue_T), intent(in) :: temperatureQuantity ! T(zeta)
  !   type (VectorValue_T), intent(in) :: refGPHQuantity ! zeta
      logical, intent(in)           ::    dontMask    ! Use even masked values
      logical, intent(in)           ::    ignoreZero  ! Ignore 0 values of h2o
      logical, intent(in)           ::    ignoreNegative  ! Ignore <0 values
      logical, intent(in)           ::    interpolate ! If VGrids or HGrids differ
      logical, intent(in)           ::    markUndefinedValues ! as UNDEFINED_VALUE
      logical, intent(in)           ::    invert      ! %RHI -> vmr if TRUE

      ! Local variables
      integer ::                          Channel     ! Channel loop counter    
      integer ::                          Chan_h2o    ! Channel loop counter    
      integer ::                          Chan_T      ! Channel loop counter    
      logical, parameter ::               DEEBUG_RHI = .false.
      integer                          :: dim
      integer ::                          I           ! Instances
      integer ::                          I_H2O       ! Instance num for values
      integer ::                          I_T         ! Instance num for values
      integer ::                          invs        ! 1 if invert, else -1
      integer ::                          QINDEX                                
      integer ::                          N           ! Num. of summed values   
      logical                          :: matched_h2o_channels
      logical                          :: matched_h2o_instances
      logical                          :: matched_sizes
      logical                          :: matched_surfs
      logical                          :: matched_T_channels
      logical                          :: matched_T_instances
      integer ::                          S           ! Surface loop counter    
      integer ::                          S_H2O       ! Instance num for surfs
      integer ::                          S_RHI       ! Instance num for surfs
      integer ::                          S_T         ! Instance num for surfs
      logical ::                          skipMe                                
      real (r8) ::                        T
      character(len=*), parameter ::      VMR_UNITS = 'vmr'
      integer ::                          VMR_UNIT_CNV
      logical ::                          wereAnySkipped
      ! These automatic arrays could cause trouble later
      ! You may consider declaring them as pointers and
      ! calling allocate_test and deallocate_test
      real (r8), dimension(quantity%template%noSurfs) :: &
       &                                  zeta, TofZeta, H2OofZeta
      real (r8), dimension(Temperaturequantity%template%noSurfs) :: &
       &                                  zetaTemperature, oldTemperature
      real (r8), dimension(sourceQuantity%template%noSurfs) :: &
       &                                  zetaH2o, oldH2o
      ! Executable statements
      ! Let any undefined values be so marked (but not necessarily masked)
      if ( markUndefinedValues ) Quantity%values = UNDEFINED_VALUE
      ! Will we convert %RHI to vmr?
      if ( invert ) then
        invs = 1
      else
        invs = -1
      end if
      ! Do we need to internally convert the vmr units?
      if ( VMR_UNITS == 'ppmv' ) then
        vmr_unit_cnv = 6
      else if ( VMR_UNITS == 'ppbv' ) then
        vmr_unit_cnv = 9
      else
        vmr_unit_cnv = 0
      end if
      ! Check that all is well
      if ( invert .and. interpolate ) then
       call Announce_Error ( key, No_Error_code, &
        & ' FillRHIFromH2O unable to invert and interpolate simultaneously' )
       return
      end if
      matched_sizes = .true.
      do dim=1, 2
        matched_sizes = matched_sizes .and. &
        & .not. any( size(Quantity%values,dim) /= &
        &(/ size(sourceQuantity%values,dim), &
        & size(temperatureQuantity%values,dim) /)&
        & )
      end do
      if ( .not. (matched_sizes .or. interpolate) ) then
       call Announce_Error ( key, No_Error_code, &
        & 'Incompatible quantities in FillRHIFromH2O--' //&
        & '(unless interpolating, all must have same shape)' )
       return
      end if
      matched_surfs = .true.
      matched_surfs = matched_surfs .and. &
       & .not. any( Quantity%template%noSurfs /= &
       &(/ sourceQuantity%template%noSurfs, &
       & temperatureQuantity%template%noSurfs /)&
       & )
      if ( .not. (matched_surfs .or. interpolate) ) then
       call Announce_Error ( key, No_Error_code, &
        & 'Different vertical coords in FillRHIFromH2O--' //&
        & '(unless interpolating, all must be on the same VGrid)' )
       return
      end if
      matched_h2o_channels = &
       &   (sourceQuantity%template%noChans == Quantity%template%noChans)
      matched_h2o_instances = &
       &   (sourceQuantity%template%noInstances == Quantity%template%noInstances)
      matched_T_channels = &
       &   (TemperatureQuantity%template%noChans == Quantity%template%noChans)
      matched_T_instances = &
       &   (TemperatureQuantity%template%noInstances == Quantity%template%noInstances)
      wereAnySkipped = .false.
      ! Now let's do the actual conversion
      do i=1, quantity%template%noInstances
        if ( quantity%template%coherent ) then
          s_rhi = 1
        else
          s_rhi = i
        end if
        if ( sourceQuantity%template%coherent ) then
          s_h2o = 1
        else
          s_h2o = i
        end if
        if ( temperaturequantity%template%coherent ) then
          s_t = 1
        else
          s_t = i
        end if
        ! zeta must be in log(hPa) units
        if ( quantity%template%verticalCoordinate == l_pressure ) then 
          zeta = -log10 ( quantity%template%surfs(:,s_rhi) )             
        else
          zeta = quantity%template%surfs(:,s_rhi)                        
        end if
        if ( interpolate .and. .not. matched_h2o_instances ) then
          i_h2o = 1
        else
          i_h2o = i
        end if
        if ( interpolate .and. .not. matched_T_instances ) then
          i_T = 1
        else
          i_T = i
        end if
        if ( sourceQuantity%template%verticalCoordinate == l_pressure ) then 
          zetah2o = -log10 ( sourceQuantity%template%surfs(:,s_h2o) )             
        else                                                           
          zetah2o = sourceQuantity%template%surfs(:,s_h2o)            
        end if
        if ( Temperaturequantity%template%verticalCoordinate == l_pressure ) then 
          zetaTemperature = -log10 ( Temperaturequantity%template%surfs(:,s_T) )             
        else                                                           
          zetaTemperature = Temperaturequantity%template%surfs(:,s_T)            
        end if
        N = 0
        do Channel=1, quantity%template%noChans
          if ( interpolate .and. .not. matched_h2o_channels ) then
            Chan_h2o = 1
          else
            Chan_h2o = Channel
          end if
          if ( interpolate .and. .not. matched_T_channels ) then
            Chan_T = 1
          else
            Chan_T = Channel
          end if
          if ( interpolate ) then
            do s=1, sourceQuantity%template%noSurfs
              qIndex = Chan_h2o + (s-1)*sourceQuantity%template%noChans
              oldH2o(s) = sourceQuantity%values(qIndex, i_h2o)
            end do
            ! Know the following about the procedure we will call:
            ! First pair of args are old(X,Y), next pair are new(X,Y)
            ! We want newY(newX) via linear interp. w/o extrapolating
            ! and mark undefined values among oldY with UNDEFINED_VALUE
            call InterpolateValues( zetah2o, oldH2o, &
             & zeta, H2OofZeta, &
             & 'Linear', extrapolate='Constant', &
             & badValue=real(UNDEFINED_VALUE, r8), &
             & missingRegions=.TRUE. )
            do s=1, Temperaturequantity%template%noSurfs
              qIndex = Chan_T + (s-1)*Temperaturequantity%template%noChans
              oldTemperature(s) = TemperatureQuantity%values(qIndex, i_T)
            end do
            call InterpolateValues( zetaTemperature, oldTemperature, &
             & zeta, TofZeta, &
             & 'Linear', extrapolate='Constant', &
             & badValue=real(UNDEFINED_VALUE, r8), &
             & missingRegions=.TRUE. )
          else
            do s=1, quantity%template%noSurfs
              N = N + 1
              qIndex = Channel + (s-1)*quantity%template%noChans
              H2OofZeta(s) = sourceQuantity%values(qIndex, i)
              TofZeta(s) = TemperatureQuantity%values(qIndex, i)
            end do
          end if
          do s=1, quantity%template%noSurfs
            N = N + 1
            qIndex = Channel + (s-1)*quantity%template%noChans
            skipMe = .false.
            if ( .not. interpolate) then
             skipMe = skipMe .or. &
             & .not. dontMask .and. ( &
             &   isVectorQtyMasked(sourceQuantity, qIndex, i) .or. &
             &   isVectorQtyMasked(temperatureQuantity, qIndex, i) &
             & )
            end if
            skipMe = skipMe .or. &
            & .not. dontMask .and. ( &
            & (ignoreNegative .and. H2OofZeta(s) < 0.0 ) &
            & .or. (ignoreZero .and. H2OofZeta(s) == 0.0 ) &
            & )
            ! But skip no matter what else if temperature illegal
            skipMe = skipMe .or. TofZeta(s) <= 0.0
            if ( .not. skipMe ) then
              T = TofZeta(s)
              Quantity%values(qIndex, i) = &
               & H2OofZeta(s) &
               & * &
               & RHIFromH2O_Factor(T, zeta(qIndex), vmr_unit_cnv, invert)
! >                & exp(invs*( &
! >                & (C(T)+zeta(qIndex)+vmr_unit_cnv) * log(10.) &
! >                & + &
! >                & 3.56654*log(T/273.16) &
! >                & ))
            end if
            wereAnySkipped = wereAnySkipped .or. skipMe
          end do
        end do
      end do
      if ( DEEBUG_RHI ) then
        call output('rhi Num. instances: ', advance='no')
        if ( invert ) then
          call output(sourceQuantity%template%noInstances, advance='yes')
        else
          call output(quantity%template%noInstances, advance='yes')
        end if
        call output('  size(surfs,2) ', advance='no')
        call output(size(quantity%template%surfs,2), advance='yes')
        call output('Were any rhi left undefined? ', advance='no')
        call output(wereAnySkipped, advance='yes')
        call dump(zeta, 'zeta(-log hPa)')
        do s=1, quantity%template%noSurfs
          if ( invert ) then
            zeta(s) = 1000000*Quantity%values(s,1)
          else
            zeta(s) = 1000000*sourceQuantity%values(s,1)
          end if
        end do
        call dump(zeta, 'H2O(ppmv)')
        call dump(TemperatureQuantity%values(:,1), 'Temperature(K)')
        if ( invert ) then
          call dump(sourceQuantity%values(:,1), 'RHI(%)')
        else
          call dump(Quantity%values(:,1), 'RHI(%)')
        end if
      end if
    end subroutine FillRHIFromH2O

      ! Properly belongs in FillRHIFromH2O, but that would make two levels
! >       function C ( T )
! >         ! As found in ref.
! >         real(r8), intent(in)   :: T
! >         real(r8)               :: C
! >         ! Local
! >         real(r8), parameter    :: a0 = -1.2141649d0
! >         real(r8), parameter    :: a1 = 9.09718d0
! >         real(r8), parameter    :: a2 = 0.876793d0
! >         real, parameter        :: ILLEGALTEMP = UNDEFINED_VALUE
! >         !
! >         if ( T > 0.d0 ) then
! >           C = a0 - a1*(273.16/T -1.0d0) + a2*(1.0d0 - T/273.16)
! >         else
! >           C = ILLEGALTEMP
! >         end if
! >       end function C

!MJF
      ! ------------------------------------- FillRHIPrecisionFromH2O ----
    subroutine FillRHIPrecisionFromH2O ( key, quantity, &
     & sourcePrecisionQuantity, tempPrecisionQuantity, sourceQuantity, temperatureQuantity, &
     & dontMask, ignoreZero, ignoreNegative, interpolate, &
     & markUndefinedValues, invert )
      ! For precisions:
      ! Convert h2o vmr to %RHI for all instances, channels, surfaces
      ! (See Eq. 9 from "UARS Microwave Limb Sounder upper tropospheric
      !  humidity measurement: Method and validation" Read et. al. 
      !  J. Geoph. Res. Dec. 2001 (106) D23)

      !  Method:
      ! (1) straight convert--all quantities must have the same shape
      !     (strictly we assume they have _all_ the same geolocations)
      ! (2) interpolate--all quantities may have different shapes
      !     (the interpolation will be along the vertical coordinate only)
      !     I.e., for xQuantity (where x can be h2o or temperature)
      !     if NoChans(xQuantity) /= NoChans(Quantity)
      !        => use only xQuantity(channel==1)
      !     if NoInstances(xQuantity) /= NoInstances(Quantity)
      !        => use only xQuantity(instance==1)
      !
      ! (3) Don't know how to invert precisions
      integer, intent(in) :: key          ! For messages
      ! Actually, the meaning of the next two is reversed if invert is TRUE)
      type (VectorValue_T), intent(inout) :: QUANTITY ! (rhi) Quantity to fill
      type (VectorValue_T), intent(in) :: sourcePrecisionQuantity ! vmr (unless invert)
      type (VectorValue_T), intent(in) :: tempPrecisionQuantity ! T(zeta)
      type (VectorValue_T), intent(in) :: sourceQuantity ! vmr (unless invert)
      type (VectorValue_T), intent(in) :: temperatureQuantity ! T(zeta)
  !   type (VectorValue_T), intent(in) :: refGPHQuantity ! zeta
      logical, intent(in)           ::    dontMask    ! Use even masked values
      logical, intent(in)           ::    ignoreZero  ! Ignore 0 values of h2o
      logical, intent(in)           ::    ignoreNegative  ! Ignore <0 values
      logical, intent(in)           ::    interpolate ! If VGrids or HGrids differ
      logical, intent(in)           ::    markUndefinedValues ! as UNDEFINED_VALUE
      logical, intent(in)           ::    invert      ! %RHI -> vmr if TRUE

      ! Local variables
      integer ::                          Channel     ! Channel loop counter    
      integer ::                          Chan_h2oPrecision    ! Channel loop counter    
      integer ::                          Chan_TPrecision      ! Channel loop counter    
      integer ::                          Chan_h2o    ! Channel loop counter    
      integer ::                          Chan_T      ! Channel loop counter    
      logical, parameter ::               DEEBUG_RHI = .false.
      integer                          :: dim
      integer ::                          I           ! Instances
      integer ::                          I_H2OPrecision       ! Instance num for values
      integer ::                          I_TPrecision         ! Instance num for values
      integer ::                          I_H2O       ! Instance num for values
      integer ::                          I_T         ! Instance num for values
      integer ::                          invs        ! 1 if invert, else -1
      integer ::                          QINDEX                                
      integer ::                          N           ! Num. of summed values   
      logical                          :: matched_h2oPrecision_channels
      logical                          :: matched_h2oPrecision_instances
      logical                          :: matched_h2o_channels
      logical                          :: matched_h2o_instances
      logical                          :: matched_sizes
      logical                          :: matched_surfs
      logical                          :: matched_TPrecision_channels
      logical                          :: matched_TPrecision_instances
      logical                          :: matched_T_channels
      logical                          :: matched_T_instances
      integer ::                          S           ! Surface loop counter    
      integer ::                          S_H2OPrecision       ! Instance num for surfs
      integer ::                          S_H2O       ! Instance num for surfs
      integer ::                          S_RHI       ! Instance num for surfs
      integer ::                          S_TPrecision         ! Instance num for surfs
      integer ::                          S_T         ! Instance num for surfs
      logical ::                          skipMe                                
!     real (r8) ::                        TPrecision
!     real (r8) ::                        T
      character(len=*), parameter ::      VMR_UNITS = 'vmr'
      integer ::                          VMR_UNIT_CNV
      logical ::                          wereAnySkipped
!     real (r8) ::                        df_db       ! RHi deriv wrt H2O
!     real (r8) ::                        df_dT       ! RHi deriv wrt T
      real (r8) ::                        rhi_precision
      ! These automatic arrays could cause trouble later
      ! You may consider declaring them as pointers and
      ! calling allocate_test and deallocate_test
      real (r8), dimension(quantity%template%noSurfs) :: &
       &                                  zeta, TPrecisionofZeta, H2OPrecisionofZeta, TofZeta, H2OofZeta
      real (r8), dimension(TempPrecisionquantity%template%noSurfs) :: &
       &                                  zetaTempPrecision, oldTempPrecision
      real (r8), dimension(sourcePrecisionQuantity%template%noSurfs) :: &
       &                                  zetaH2oPrecision, oldH2oPrecision
      real (r8), dimension(Temperaturequantity%template%noSurfs) :: &
       &                                  zetaTemperature, oldTemperature
      real (r8), dimension(sourceQuantity%template%noSurfs) :: &
       &                                  zetaH2o, oldH2o
      ! Executable statements
      ! Let any undefined values be so marked (but not necessarily masked)
      if ( markUndefinedValues ) Quantity%values = UNDEFINED_VALUE
      ! Will we convert %RHI to vmr?
      if ( invert ) then
       call Announce_Error ( key, No_Error_code, &
        & ' FillRHIPrecisionFromH2O unable to invert' )
       return
      end if
!!$      if ( invert ) then
!!$        invs = 1
!!$      else
        invs = -1
!!$      end if
      ! Do we need to internally convert the vmr units?
      if ( VMR_UNITS == 'ppmv' ) then
        vmr_unit_cnv = 6
      else if ( VMR_UNITS == 'ppbv' ) then
        vmr_unit_cnv = 9
      else
        vmr_unit_cnv = 0
      end if
      ! Check that all is well
      if ( invert .and. interpolate ) then
       call Announce_Error ( key, No_Error_code, &
        & ' FillRHIPrecisionFromH2O unable to invert and interpolate simultaneously' )
       return
      end if
      matched_sizes = .true.
      do dim=1, 2
        matched_sizes = matched_sizes .and. &
        & .not. any( size(Quantity%values,dim) /= &
        &(/ size(sourcePrecisionQuantity%values,dim), &
        & size(tempPrecisionQuantity%values,dim), &
        & size(sourceQuantity%values,dim), &
        & size(temperatureQuantity%values,dim) /)&
        & )
      end do
      if ( .not. (matched_sizes .or. interpolate) ) then
       call Announce_Error ( key, No_Error_code, &
        & 'Incompatible quantities in FillRHIPrecisionFromH2O--' //&
        & '(unless interpolating, all must have same shape)' )
       return
      end if
      matched_surfs = .true.
      matched_surfs = matched_surfs .and. &
       & .not. any( Quantity%template%noSurfs /= &
       &(/ sourcePrecisionQuantity%template%noSurfs, &
       & tempPrecisionQuantity%template%noSurfs, &
       & sourceQuantity%template%noSurfs, &
       & temperatureQuantity%template%noSurfs /)&
       & )
      if ( .not. (matched_surfs .or. interpolate) ) then
       call Announce_Error ( key, No_Error_code, &
        & 'Different vertical coords in FillRHIPrecisionFromH2O--' //&
        & '(unless interpolating, all must be on the same VGrid)' )
       return
      end if
      matched_h2oPrecision_channels = &
       &   (sourcePrecisionQuantity%template%noChans == Quantity%template%noChans)
      matched_h2oPrecision_instances = &
       &   (sourcePrecisionQuantity%template%noInstances == Quantity%template%noInstances)
      matched_TPrecision_channels = &
       &   (TempPrecisionQuantity%template%noChans == Quantity%template%noChans)
      matched_TPrecision_instances = &
       &   (TempPrecisionQuantity%template%noInstances == Quantity%template%noInstances)
      matched_h2o_channels = &
       &   (sourceQuantity%template%noChans == Quantity%template%noChans)
      matched_h2o_instances = &
       &   (sourceQuantity%template%noInstances == Quantity%template%noInstances)
      matched_T_channels = &
       &   (TemperatureQuantity%template%noChans == Quantity%template%noChans)
      matched_T_instances = &
       &   (TemperatureQuantity%template%noInstances == Quantity%template%noInstances)
      wereAnySkipped = .false.
      ! Now let's do the actual conversion
      do i=1, quantity%template%noInstances
        if ( quantity%template%coherent ) then
          s_rhi = 1
        else
          s_rhi = i
        end if
        if ( sourcePrecisionQuantity%template%coherent ) then
          s_h2oPrecision = 1
        else
          s_h2oPrecision = i
        end if
        if ( tempPrecisionquantity%template%coherent ) then
          s_tPrecision = 1
        else
          s_tPrecision = i
        end if
        if ( sourceQuantity%template%coherent ) then
          s_h2o = 1
        else
          s_h2o = i
        end if
        if ( temperaturequantity%template%coherent ) then
          s_t = 1
        else
          s_t = i
        end if
        ! zeta must be in log(hPa) units
        if ( quantity%template%verticalCoordinate == l_pressure ) then 
          zeta = -log10 ( quantity%template%surfs(:,s_rhi) )             
        else
          zeta = quantity%template%surfs(:,s_rhi)                        
        end if
        if ( interpolate .and. .not. matched_h2oPrecision_instances ) then
          i_h2oPrecision = 1
        else
          i_h2oPrecision = i
        end if
        if ( interpolate .and. .not. matched_TPrecision_instances ) then
          i_TPrecision = 1
        else
          i_TPrecision = i
        end if
        if ( interpolate .and. .not. matched_h2o_instances ) then
          i_h2o = 1
        else
          i_h2o = i
        end if
        if ( interpolate .and. .not. matched_T_instances ) then
          i_T = 1
        else
          i_T = i
        end if
        if ( sourcePrecisionQuantity%template%verticalCoordinate == l_pressure ) then 
          zetah2oPrecision = -log10 ( sourcePrecisionQuantity%template%surfs(:,s_h2o) )             
        else                                                           
          zetah2oPrecision = sourcePrecisionQuantity%template%surfs(:,s_h2o)            
        end if
        if ( TempPrecisionquantity%template%verticalCoordinate == l_pressure ) then 
          zetaTempPrecision = -log10 ( TempPrecisionquantity%template%surfs(:,s_T) )             
        else                                                           
          zetaTempPrecision = TempPrecisionquantity%template%surfs(:,s_T)            
        end if
        if ( sourceQuantity%template%verticalCoordinate == l_pressure ) then 
          zetah2o = -log10 ( sourceQuantity%template%surfs(:,s_h2o) )             
        else                                                           
          zetah2o = sourceQuantity%template%surfs(:,s_h2o)            
        end if
        if ( Temperaturequantity%template%verticalCoordinate == l_pressure ) then 
          zetaTemperature = -log10 ( Temperaturequantity%template%surfs(:,s_T) )             
        else                                                           
          zetaTemperature = Temperaturequantity%template%surfs(:,s_T)            
        end if
        N = 0
        do Channel=1, quantity%template%noChans
          if ( interpolate .and. .not. matched_h2oPrecision_channels ) then
            Chan_h2oPrecision = 1
          else
            Chan_h2oPrecision = Channel
          end if
          if ( interpolate .and. .not. matched_TPrecision_channels ) then
            Chan_TPrecision = 1
          else
            Chan_TPrecision = Channel
          end if
          if ( interpolate .and. .not. matched_h2o_channels ) then
            Chan_h2o = 1
          else
            Chan_h2o = Channel
          end if
          if ( interpolate .and. .not. matched_T_channels ) then
            Chan_T = 1
          else
            Chan_T = Channel
          end if
          if ( interpolate ) then
            do s=1, sourcePrecisionQuantity%template%noSurfs
              qIndex = Chan_h2oPrecision + (s-1)*sourcePrecisionQuantity%template%noChans
              oldH2oPrecision(s) = sourcePrecisionQuantity%values(qIndex, i_h2oPrecision)
            end do
            ! Know the following about the procedure we will call:
            ! First pair of args are old(X,Y), next pair are new(X,Y)
            ! We want newY(newX) via linear interp. w/o extrapolating
            ! and mark undefined values among oldY with UNDEFINED_VALUE
            call InterpolateValues( zetah2oPrecision, oldH2oPrecision, &
             & zeta, H2OPrecisionofZeta, &
             & 'Linear', extrapolate='Constant', &
             & badValue=real(UNDEFINED_VALUE, r8), &
             & missingRegions=.TRUE. )
            do s=1, TempPrecisionquantity%template%noSurfs
              qIndex = Chan_TPrecision + (s-1)*TempPrecisionquantity%template%noChans
              oldTempPrecision(s) = TempPrecisionQuantity%values(qIndex, i_TPrecision)
            end do
            call InterpolateValues( zetaTempPrecision, oldTempPrecision, &
             & zeta, TPrecisionofZeta, &
             & 'Linear', extrapolate='Constant', &
             & badValue=real(UNDEFINED_VALUE, r8), &
             & missingRegions=.TRUE. )
            do s=1, sourceQuantity%template%noSurfs
              qIndex = Chan_h2o + (s-1)*sourceQuantity%template%noChans
              oldH2o(s) = sourceQuantity%values(qIndex, i_h2o)
            end do
            ! Know the following about the procedure we will call:
            ! First pair of args are old(X,Y), next pair are new(X,Y)
            ! We want newY(newX) via linear interp. w/o extrapolating
            ! and mark undefined values among oldY with UNDEFINED_VALUE
            call InterpolateValues( zetah2o, oldH2o, &
             & zeta, H2OofZeta, &
             & 'Linear', extrapolate='Constant', &
             & badValue=real(UNDEFINED_VALUE, r8), &
             & missingRegions=.TRUE. )
            do s=1, Temperaturequantity%template%noSurfs
              qIndex = Chan_T + (s-1)*Temperaturequantity%template%noChans
              oldTemperature(s) = TemperatureQuantity%values(qIndex, i_T)
            end do
            call InterpolateValues( zetaTemperature, oldTemperature, &
             & zeta, TofZeta, &
             & 'Linear', extrapolate='Constant', &
             & badValue=real(UNDEFINED_VALUE, r8), &
             & missingRegions=.TRUE. )
          else
            do s=1, quantity%template%noSurfs
              N = N + 1
              qIndex = Channel + (s-1)*quantity%template%noChans
              H2OPrecisionofZeta(s) = sourcePrecisionQuantity%values(qIndex, i)
              TPrecisionofZeta(s) = TempPrecisionQuantity%values(qIndex, i)
              H2OofZeta(s) = sourceQuantity%values(qIndex, i)
              TofZeta(s) = TemperatureQuantity%values(qIndex, i)
            end do
          end if
          do s=1, quantity%template%noSurfs
            N = N + 1
            qIndex = Channel + (s-1)*quantity%template%noChans
            skipMe = .false.
            if ( .not. interpolate) then
             skipMe = skipMe .or. &
             & .not. dontMask .and. ( &
             &   isVectorQtyMasked(sourcePrecisionQuantity, qIndex, i) .or. &
             &   isVectorQtyMasked(tempPrecisionQuantity, qIndex, i) .or. &
             &   isVectorQtyMasked(sourceQuantity, qIndex, i) .or. &
             &   isVectorQtyMasked(temperatureQuantity, qIndex, i) &
             & )
            end if
            skipMe = skipMe .or. &
            & .not. dontMask .and. ( &
            & (ignoreNegative .and. H2OofZeta(s) < 0.0 ) &
            & .or. (ignoreZero .and. H2OofZeta(s) == 0.0 ) &
            & )
!!$            skipMe = skipMe .or. &
!!$            & .not. dontMask .and. ( &
!!$            & (ignoreNegative .and. H2OPrecisionofZeta(s) < 0.0 ) &
!!$            & .or. (ignoreZero .and. H2OPrecisionofZeta(s) == 0.0 ) &
!!$            & .or. (ignoreNegative .and. H2OofZeta(s) < 0.0 ) &
!!$            & .or. (ignoreZero .and. H2OofZeta(s) == 0.0 ) &
!!$            & )
            ! But skip no matter what else if temperature illegal
            skipMe = skipMe .or. TofZeta(s) <= 0.0
            if ( .not. skipMe ) then
!               T = TofZeta(s)
!               df_db = exp(invs*( &
!                & (C(T)+zeta(qIndex)+vmr_unit_cnv) * log(10.) &
!                & + &
!                & 3.56654*log(T/273.16) &
!                & ))
!               df_dT = H2OofZeta(s) * exp(invs*( &
!                & (C(T)+zeta(qIndex)+vmr_unit_cnv) * log(10.) &
!                & + &
!                & 3.56654*log(T/273.16) &
!                & )) &
!                & * invs * ( dC_dT(T) * log(10.) + 3.56654 / T )
!               Quantity%values(qIndex, i) = sqrt (&
!                & ( H2OPrecisionofZeta(s) * df_db )**2 &
!                & + ( TPrecisionofZeta(s) * df_dT )**2 &
!                & )
              call RHIPrecFromH2O(H2OofZeta(s), &
               & TofZeta(s), zeta(qIndex), vmr_unit_cnv, &
               & H2OPrecisionofZeta(s), TPrecisionofZeta(s), &
               & rhi_precision)
              Quantity%values(qIndex, i) = rhi_precision
            end if
            wereAnySkipped = wereAnySkipped .or. skipMe
          end do
        end do
      end do
      if ( DEEBUG_RHI ) then
        call output('rhi Num. instances: ', advance='no')
        if ( invert ) then
          call output(sourceQuantity%template%noInstances, advance='yes')
        else
          call output(quantity%template%noInstances, advance='yes')
        end if
        call output('  size(surfs,2) ', advance='no')
        call output(size(quantity%template%surfs,2), advance='yes')
        call output('Were any rhi left undefined? ', advance='no')
        call output(wereAnySkipped, advance='yes')
        call dump(zeta, 'zeta(-log hPa)')
        do s=1, quantity%template%noSurfs
          if ( invert ) then
            zeta(s) = 1000000*Quantity%values(s,1)
          else
            zeta(s) = 1000000*sourceQuantity%values(s,1)
          end if
        end do
        call dump(zeta, 'H2O(ppmv)')
        call dump(TemperatureQuantity%values(:,1), 'Temperature(K)')
        if ( invert ) then
          call dump(sourceQuantity%values(:,1), 'RHI(%)')
        else
          call dump(Quantity%values(:,1), 'RHI(%)')
        end if
      end if
    end subroutine FillRHIPrecisionFromH2O

      ! Properly belongs in FillRHIPrecisionFromH2O, but that would make two levels
!       function dC_dT ( T )
!         ! As found in ref.
!         real(r8), intent(in)   :: T
!         real(r8)               :: dC_dT
!         ! Local
!         real(r8), parameter    :: a0 = -1.2141649d0
!         real(r8), parameter    :: a1 = 9.09718d0
!         real(r8), parameter    :: a2 = 0.876793d0
!         real, parameter        :: ILLEGALTEMP = UNDEFINED_VALUE
!         !
!         if ( T > 0.d0 ) then
!           dC_dT = a1*(273.16/T**2) - a2/273.16
! !!$          C = a0 - a1*(273.16/T -1.0d0) + a2*(1.0d0 - T/273.16)
!         else
!           dC_dT = ILLEGALTEMP
! !!$          C = ILLEGALTEMP
!         end if
!       end function dC_dT

!MJF
    ! ---------------------------------- FillVectorQuantityWithEsimatedNoise ---
    subroutine FillVectorQtyWithEstNoise ( quantity, radiance, &
      & sysTemp, nbw, integrationTime )

      use MLSSignals_m, only: signals

      ! Dummy arguments
      type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      type (VectorValue_T), intent(in) :: RADIANCE ! Radiances to use in calculation
      type (VectorValue_T), intent(in) :: SYSTEMP ! System temperature
      type (VectorValue_T), pointer :: NBW ! Noise bandwidth
      real(r8), intent(in) :: INTEGRATIONTIME ! Integration time in seconds

      ! Local variables
      integer :: C                        ! Channel loop counter
      integer :: S                        ! Surface loop counter
      integer :: I                        ! Index into first dimension of values

      real (r8), dimension(:), pointer :: WIDTH ! Channel widths in MHz

      ! Executable code

      if (.not. ValidateVectorQuantity ( quantity, &
        & quantityType=(/l_radiance/), minorFrame=.true.) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Invalid quantity for estimated noise fill.")

      if (.not. ValidateVectorQuantity ( radiance, &
        & quantityType=(/l_radiance/), minorFrame=.true.) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Invalid quantity for estimated noise fill.")

      if ( ( radiance%template%signal /= quantity%template%signal ) .or. &
        &  ( radiance%template%sideband /= quantity%template%sideband ) ) &
        & call MLSMEssage ( MLSMSG_Error, ModuleName, &
        & "Quantity and radiances not same signal for estimated noise fill.")

      if ( ( systemp%template%signal /= quantity%template%signal ) .or. &
        &  ( systemp%template%sideband /= quantity%template%sideband ) ) &
        & call MLSMEssage ( MLSMSG_Error, ModuleName, &
        & "Quantity and system temperature not same signal for estimated noise fill." )

      if (.not. ValidateVectorQuantity ( &
        & sysTemp, &
        & quantityType=(/l_systemTemperature/), &
        & verticalCoordinate=(/l_none/), &
        & noInstances=(/1/) ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Invalid system temperature quantity for estimated noise fill.")

      if ( associated ( nbw ) ) then 
        if (.not. ValidateVectorQuantity ( &
          & nbw, &
          & quantityType=(/l_noiseBandwidth/), &
          & verticalCoordinate=(/l_none/), &
          & noInstances=(/1/) ) ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Invalid noise bandwidth quantity for estimated noise fill.")
        if ( ( nbw%template%signal /= quantity%template%signal ) .or. &
          &    ( nbw%template%sideband /= quantity%template%sideband ) ) &
          & call MLSMEssage ( MLSMSG_Error, ModuleName, &
          & "Quantity and noise bandwidth not same signal for estimated noise fill." )
        width => nbw%values(:,1)
      else
        width => signals(radiance%template%signal)%widths(:)
      end if

      i = 1
      do s = 1, quantity%template%noSurfs
        do c = 1, quantity%template%noChans
          quantity%values(i,:) = ( radiance%values(i,:) + sysTemp%values(c,1) ) / &
            & sqrt ( integrationTime * 1e6 * width(c) )
          i = i + 1
        end do
      end do

    end subroutine FillVectorQtyWithEstNoise

    ! ------------------------------------- FillVectorHydrostatically ----
    subroutine FillVectorQtyHydrostatically ( key, quantity, &
      & temperatureQuantity, refGPHQuantity, h2oQuantity, &
      & orbitInclinationQuantity, phiTanQuantity, geocAltitudeQuantity, maxIterations, &
      & phiWindow, phiWindowUnits )
      ! Various hydrostatic fill operations
      integer, intent(in) :: key          ! For messages
      type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      type (VectorValue_T), intent(in) :: TEMPERATUREQUANTITY
      type (VectorValue_T), intent(in) :: REFGPHQUANTITY
      type (VectorValue_T), pointer :: H2OQUANTITY
      type (VectorValue_T), pointer :: ORBITINCLINATIONQUANTITY
      type (VectorValue_T), pointer :: PHITANQUANTITY
      type (VectorValue_T), pointer :: GEOCALTITUDEQUANTITY
      integer, intent(in) :: MAXITERATIONS
      real(r8), intent(in) :: PHIWINDOW
      integer, intent(in) :: PHIWINDOWUNITS
      ! H2OQuantity and GeocAltitudeQuantity have to be pointers
      ! as they may be absent.

      ! Local variables

      ! Executable code

      if ( toggle(gen) .and. levels(gen) > 0 ) &
        & call trace_begin ( "FillVectorQtyHydrostatically", key )

      select case ( quantity%template%quantityType )
      case ( l_gph )
        if ( (temperatureQuantity%template%noSurfs /= &
          &   quantity%template%noSurfs) .or. &
          &  (refGPHQuantity%template%noInstances /= &
          &   quantity%template%noInstances) .or. &
          &  (temperatureQuantity%template%noInstances /= &
          &   quantity%template%noInstances) ) then
          call Announce_Error ( key, nonConformingHydrostatic, &
            & "case l_gph failed first test" )
	  if ( toggle(gen) .and. levels(gen) > 0 ) &
            & call trace_end ( "FillVectorQtyHydrostatically")
          return
        end if
        if ( (any(quantity%template%surfs /= temperatureQuantity%template%surfs)) .or. &
          & (any(quantity%template%phi /= temperatureQuantity%template%phi)) .or. &
          & (any(quantity%template%phi /= refGPHQuantity%template%phi)) ) then
          call Announce_Error ( key, nonConformingHydrostatic, &
            &  "case l_gph failed second test" )
	  if ( toggle(gen) .and. levels(gen) > 0 ) &
            & call trace_end ( "FillVectorQtyHydrostatically")
          return
        end if
        call GetBasisGPH ( temperatureQuantity, refGPHQuantity, quantity%values )
      case ( l_ptan )
        if ( (temperatureQuantity%template%noInstances /= &
          &   refGPHquantity%template%noInstances) .or. &
          &  (temperatureQuantity%template%noInstances /= &
          &   h2oQuantity%template%noInstances) ) then
          call Announce_Error ( key, nonConformingHydrostatic, &
            & "case l_ptan failed first test" )
	  if ( toggle(gen) .and. levels(gen) > 0 ) &
            & call trace_end ( "FillVectorQtyHydrostatically")
          return
        end if
        if ( (any(refGPHquantity%template%phi /= temperatureQuantity%template%phi)) .or. &
          & (any(h2oQuantity%template%phi /= temperatureQuantity%template%phi)) ) then
          call Announce_Error ( key, nonConformingHydrostatic, &
            & "case l_ptan failed second test" )
	  if ( toggle(gen) .and. levels(gen) > 0 )&
            &  call trace_end ( "FillVectorQtyHydrostatically")
          return
        end if
        if ( (.not. ValidateVectorQuantity(quantity, minorFrame=.true.) ) .or. &
          &  (.not. ValidateVectorQuantity(geocAltitudeQuantity, minorFrame=.true.) ) .or. &
          &  (quantity%template%instrumentModule /= &
          &   geocAltitudeQuantity%template%instrumentModule) )  then
          call Announce_Error ( key, nonConformingHydrostatic, &
            & "case l_ptan failed third test" )
  !        print *, 'ValidateVectorQuantity(quantity, minorFrame=.true.) ', &
  !        &  ValidateVectorQuantity(quantity, minorFrame=.true., sayWhyNot=.true.)
  !        print *, 'ValidateVectorQuantity(geocAltitudeQuantity, minorFrame=.true.) ', &
  !        & ValidateVectorQuantity(geocAltitudeQuantity, minorFrame=.true., sayWhyNot=.true.)
  !        print *, 'quantity%template%instrumentModule ', &
  !        & quantity%template%instrumentModule
  !        print *, 'geocAltitudeQuantity%template%instrumentModule ', &
  !        & geocAltitudeQuantity%template%instrumentModule
           call output('ValidateVectorQuantity(quantity, minorFrame=.true.) ')
           call blanks(3)
           call output( &
           & ValidateVectorQuantity(quantity, minorFrame=.true., sayWhyNot=.true.), &
           & advance='yes')

           call output('ValidateVectorQuantity(geocAltitudeQuantity, minorFrame=.true.) ')
           call blanks(3)
           call output( &
           & ValidateVectorQuantity(geocAltitudeQuantity, minorFrame=.true., sayWhyNot=.true.), &
           & advance='yes')

           call output('quantity%template%instrumentModule ')
           call blanks(3)
           call output( &
           & quantity%template%instrumentModule, &
           & advance='yes')

           call output('geocAltitudeQuantity%template%instrumentModule ')
           call blanks(3)
           call output( &
           & geocAltitudeQuantity%template%instrumentModule, &
           & advance='yes')
	  if ( toggle(gen) .and. levels(gen) > 0 ) &
            & call trace_end ( "FillVectorQtyHydrostatically")
          return
        end if
        call Get2DHydrostaticTangentPressure ( quantity, temperatureQuantity,&
          & refGPHQuantity, h2oQuantity, orbitInclinationQuantity, &
          & phiTanQuantity, geocAltitudeQuantity, maxIterations, &
          & phiWindow, phiWindowUnits )
      case default
        call Announce_error ( 0, 0, 'No such fill yet' )
      end select

      if ( toggle(gen) .and. levels(gen) > 0 ) &
        & call trace_end ( "FillVectorQtyHydrostatically" )

    end subroutine FillVectorQtyHydrostatically

    ! ------------------------------------- FillVectorHydrostatically ----
    subroutine FillGPHPrecision ( key, quantity, &
      & tempPrecisionQuantity, refGPHPrecisionQuantity )
      ! Fill the GPH precision from the temperature and refGPH precision,
      ! ignoring the of diagonal elements (not available outside
      ! RetrievalModule anyway).
      integer, intent(in) :: key          ! For messages
      type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      type (VectorValue_T), intent(in) :: TEMPPRECISIONQUANTITY
      type (VectorValue_T), intent(in) :: REFGPHPRECISIONQUANTITY

      ! Local variables

      ! Executable code

      if ( toggle(gen) .and. levels(gen) > 0 ) &
        & call trace_begin ( "FillGPHPrecision", key )

      select case ( quantity%template%quantityType )
      case ( l_gph )
        if ( (tempPrecisionQuantity%template%noSurfs /= &
          &   quantity%template%noSurfs) .or. &
          &  (refGPHPrecisionQuantity%template%noInstances /= &
          &   quantity%template%noInstances) .or. &
          &  (tempPrecisionQuantity%template%noInstances /= &
          &   quantity%template%noInstances) ) then
          call Announce_Error ( key, nonConformingHydrostatic, &
            & "case l_gph failed first test" )
	  if ( toggle(gen) .and. levels(gen) > 0 ) &
            & call trace_end ( "FillGPHPrecision")
          return
        end if
        if ( (any(quantity%template%surfs /= tempPrecisionQuantity%template%surfs)) .or. &
          & (any(quantity%template%phi /= tempPrecisionQuantity%template%phi)) .or. &
          & (any(quantity%template%phi /= refGPHPrecisionQuantity%template%phi)) ) then
          call Announce_Error ( key, nonConformingHydrostatic, &
            &  "case l_gph failed second test" )
	  if ( toggle(gen) .and. levels(gen) > 0 ) &
            & call trace_end ( "FillGPHPrecision")
          return
        end if
        call GetGPHPrecision ( tempPrecisionQuantity, refGPHPrecisionQuantity, quantity%values )
      case default
        call Announce_error ( 0, 0, 'GPH precision needed for result of FillGPHPrecision' )
      end select

      if ( toggle(gen) .and. levels(gen) > 0 ) &
        & call trace_end ( "FillGPHPrecision" )

    end subroutine FillGPHPrecision

    ! -------------------------------------- FillVectorQtyFromIsotope -----------

    subroutine FillVectorQtyFromIsotope ( key, quantity, sourceQuantity, &
              & ratioQuantity )
      ! This routine fills one vector from another, given an appropriate
      ! isotope ratio.

      integer, intent(in) :: KEY          ! Tree node
      type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      type (VectorValue_T), intent(in) :: SOURCEQUANTITY ! Quantity to take vmr from
      type (VectorValue_T), intent(in) :: RATIOQUANTITY ! Isotope ratio information

      ! Local variables
      real (r8) :: FACTOR                 ! Multiplier to apply to sourceQuantity

      ! Executable code

      if (.not. ValidateVectorQuantity ( quantity, &
        & quantityType=(/ l_vmr /), frequencyCoordinate=(/ l_none /) ) ) &
        &   call MLSMessage ( MLSMSG_Error, ModuleName, &
        &      "Inappropriate quantity for isotope fill")

      if (.not. ValidateVectorQuantity ( sourceQuantity, &
        & quantityType=(/ l_vmr /), frequencyCoordinate=(/ l_none /) ) ) &
        &   call MLSMessage ( MLSMSG_Error, ModuleName, &
        &      "Inappropriate source quantity for isotope fill")

      if (.not. ValidateVectorQuantity ( ratioQuantity, &
        & quantityType=(/ l_isotopeRatio /), frequencyCoordinate=(/ l_none /), &
        & noInstances=(/1/), noSurfs=(/1/) ) ) &
        &   call MLSMessage ( MLSMSG_Error, ModuleName, &
        &     "Inappropriate form/quantity for isotope ratio")


      if ( .not. DoHGridsMatch ( quantity, sourceQuantity ) .or. &
        &  .not. DoVGridsMatch ( quantity, sourceQuantity ) ) &
        &    call MLSMessage ( MLSMSG_Error, ModuleName, &
        &      "Quantity and source quantity don't match for isotope fill" )

      if ( quantity%template%molecule == sourceQuantity%template%molecule ) &
        & call MLSMessage( MLSMSG_Error, ModuleName, &
        &   "Source and quantity both describe same molecule in isotope fill")

      if ( ratioQuantity%template%molecule == quantity%template%molecule ) then
        ! Going from parent to isotope
        factor = ratioQuantity%values(1,1)
      else if ( ratioQuantity%template%molecule == sourceQuantity%template%molecule ) then
        ! Going from isotope to parent
        factor = 1.0/ratioQuantity%values(1,1)
      else 
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Unable to understand isotope fill" )
      end if

      quantity%values = sourceQuantity%values * factor

    end subroutine FillVectorQtyFromIsotope

    !=============================================== ExplicitFillVectorQuantity ==
    subroutine ExplicitFillVectorQuantity(quantity, valuesNode, spreadFlag, globalUnit, &
      & dontmask)

      ! This routine is called from MLSL2Fill to fill values from an explicit
      ! fill command line

      ! Dummy arguments
      type (VectorValue_T), intent(inout) :: QUANTITY ! The quantity to fill
      integer, intent(in) :: VALUESNODE   ! Tree node
      logical, intent(in) :: SPREADFLAG   ! One instance given, spread to all
      integer, intent(in) :: GLOBALUNIT   ! From parent vector
      logical, intent(in) :: DONTMASK     ! Don't bother with the mask

      ! Local variables
      integer :: K                        ! Loop counter
      integer :: I,J                      ! Other indices
      integer, DIMENSION(2) :: unitAsArray ! Unit for value given
      real (r8), DIMENSION(2) :: valueAsArray ! Value given
      integer :: TestUnit                 ! Unit to use

      ! Executable code
      testUnit = quantity%template%unit
      if ( globalUnit /= phyq_Invalid ) testUnit = globalUnit

      if (spreadFlag) then ! 1 instance/value given, spread to all instances

        ! Check we have the right number of values
        if ( ( (nsons(valuesNode)-1 /= quantity%template%instanceLen) .and. &
          &    (nsons(valuesNode)-1 /= 1) ) .or. &
          &  (.not. quantity%template%regular)) &
          & call Announce_error ( valuesNode, invalidExplicitFill )

        ! Loop over the values
        do k=1,nsons(valuesNode)-1
          ! Get value from tree
          call expr(subtree(k+1,valuesNode),unitAsArray,valueAsArray)
          ! Check unit OK
          if ( (unitAsArray(1) /= testUnit) .and. &
            &  (unitAsArray(1) /= PHYQ_Dimensionless) ) &
            & call Announce_error ( valuesNode, badUnitsForExplicit )
          ! Store value
          if ( .not. dontMask .and. associated ( quantity%mask ) ) then
            if ( nsons(valuesNode)-1 == 1 ) then
              where ( iand ( ichar(quantity%mask), m_Fill ) == 0 )
                quantity%values=valueAsArray(1)
              end where
            else
              where ( iand ( ichar(quantity%mask(k,:)), m_Fill ) == 0 )
                quantity%values(k,:)=valueAsArray(1)
              end where
            end if
          else
            ! No mask to worry about
            if (nsons(valuesNode)-1 == 1) then
              quantity%values=valueAsArray(1)
            else
              quantity%values(k,:)=valueAsArray(1)
            end if
          end if
        end do

      else                              ! Not spread, fill all values

        ! Check we have the right number of values
        if (nsons(valuesNode)-1 /= &
          & quantity%template%noInstances*quantity%template%instanceLen) &
          & call Announce_error ( valuesNode, invalidExplicitFill )

        ! Loop over values
        do k=1,nsons(valuesNode)-1
          ! Get value from tree
          call expr(subtree(k+1,valuesNode),unitAsArray,valueAsArray)
          ! Check unit OK
          if ( (unitAsArray(1) /= testUnit) .and. &
            &  (unitAsArray(1) /= PHYQ_Dimensionless) ) &
            & call Announce_error ( valuesNode, badUnitsForExplicit )
          ! Store value
          i = mod(k-1,quantity%template%instanceLen) + 1
          j = (k-1) / quantity%template%instanceLen + 1
          if ( .not. dontMask .and. associated ( quantity%mask ) ) then
            if ( iand ( ichar(quantity%mask(i,j)), m_Fill ) == 0 ) &
              & quantity%values(i,j) = valueAsArray(1)
          else
            quantity%values(i,j) = valueAsArray(1)
          end if
        end do
      end if
    end subroutine ExplicitFillVectorQuantity

    ! ----------------------------------------- FillVectorQuantityFromL1B ----
    subroutine FillVectorQuantityFromL1B ( root, quantity, chunk, l1bInfo, &
      & isPrecision, PrecisionQuantity )
      integer, intent(in) :: root
      type (VectorValue_T), INTENT(INOUT) ::        QUANTITY
      type (MLSChunk_T), INTENT(IN) ::              CHUNK
      type (l1bInfo_T), INTENT(IN) ::               L1BINFO
      logical, intent(in)               ::          ISPRECISION
      type (VectorValue_T), INTENT(IN), optional :: PRECISIONQUANTITY

      ! Local variables
      character (len=80) :: NAMESTRING
      character (len=FileNameLen) :: FILENAMESTRING
      integer :: fileID, FLAG, NOMAFS
      type (l1bData_T) :: L1BDATA
      integer :: ROW, COLUMN
      integer :: this_hdfVersion

      ! Executable code

      if ( toggle(gen) .and. levels(gen) > 0 ) &
        & call trace_begin ("FillVectorQuantityFromL1B",root)
      fileID=l1bInfo%l1bOAID
      if ( quantity%template%quantityType /= l_radiance ) then
        filenamestring = l1bInfo%L1BOAFileName
      else
        filenamestring = l1bInfo%L1BRADFileNames(1)
      endif
      this_hdfVersion = mls_hdf_version(trim(filenamestring), LEVEL1_HDFVERSION)
      if ( this_hdfVersion == ERRORINH5FFUNCTION ) then
        call Announce_Error ( root, No_Error_code, &
          & 'Error in finding hdfversion of l1b file' )
      elseif ( this_hdfVersion == WRONGHDFVERSION ) then
        call Announce_Error ( root, No_Error_code, &
          & 'Wrong hdfversion declared/coded for l1b file' )
      endif
      
      select case ( quantity%template%quantityType )
      case ( l_ptan )
        call GetModuleName( quantity%template%instrumentModule,nameString )
!        nameString=TRIM(nameString)//'.ptan'
        nameString = AssembleL1BQtyName('ptan', this_hdfVersion, .FALSE., &
          & trim(nameString))
      case ( l_radiance )
        call GetSignalName ( quantity%template%signal, nameString, &
          & sideband=quantity%template%sideband, noChannels=.TRUE. )
        nameString = AssembleL1BQtyName(nameString, this_hdfVersion, .FALSE.)
        fileID = FindL1BData (l1bInfo%l1bRadIDs, nameString, this_hdfVersion )
      case ( l_tngtECI )
        call GetModuleName( quantity%template%instrumentModule,nameString )
!        nameString=TRIM(nameString)//'.tpECI'
        nameString = AssembleL1BQtyName('ECI', this_hdfVersion, .TRUE., &
          & trim(nameString))
      case ( l_tngtGeodAlt )
        call GetModuleName( quantity%template%instrumentModule,nameString )
!        nameString=TRIM(nameString)//'.tpGeodAlt'
        nameString = AssembleL1BQtyName('GeodAlt', this_hdfVersion, .TRUE., &
          & trim(nameString))
      case ( l_tngtGeocAlt )
        call GetModuleName( quantity%template%instrumentModule,nameString )
!        nameString=TRIM(nameString)//'.tpGeocAlt'
        nameString = AssembleL1BQtyName('GeocAlt', this_hdfVersion, .TRUE., &
          & trim(nameString))
      case ( l_scECI )
!        nameString='scECI'
        nameString = AssembleL1BQtyName('ECI', this_hdfVersion, .FALSE., 'sc')
      case ( l_scVel )
!        nameString='scVel'
        nameString = AssembleL1BQtyName('Vel', this_hdfVersion, .FALSE., 'sc')
      case ( l_scVelECI )
!        nameString='scVelECI'
        nameString = AssembleL1BQtyName('VelECI', this_hdfVersion, .FALSE., &
          & 'sc')
      case ( l_scVelECR )
!        nameString='scVelECR'
        nameString = AssembleL1BQtyName('VelECR', this_hdfVersion, .FALSE., &
          & 'sc')
      case ( l_scGeocAlt )
!        nameString='scGeocAlt'
        nameString = AssembleL1BQtyName('GeocAlt', this_hdfVersion, .FALSE., &
          & 'sc')
      case ( l_orbitInclination )
!        nameString='scOrbIncl'
        nameString = AssembleL1BQtyName('OrbIncl', this_hdfVersion, .FALSE., &
          & 'sc')
      case default
        call Announce_Error ( root, cantFillFromL1B )
      end select

      if ( isPrecision ) nameString = trim(nameString) // PRECISIONSUFFIX

      call ReadL1BData ( fileID , nameString, l1bData, noMAFs, flag, &
        & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
        & NeverFail= .false., hdfVersion=this_hdfVersion )
      ! We'll have to think about `bad' values here .....
      if ( flag /= 0 ) then
        call Announce_Error ( root, errorReadingL1B )
        if ( toggle(gen) .and. levels(gen) > 0 ) &
          & call trace_end ( "FillVectorQuantityFromL1B")
        return
      end if
      quantity%values = RESHAPE(l1bData%dpField, &
        & (/ quantity%template%instanceLen, quantity%template%noInstances /) )
      if ( isPrecision ) then
        do column=1, size(quantity%values(1, :))
          do row=1, size(quantity%values(:, 1))
            if ( quantity%values(row, column) < 0.d0 ) &
              & call MaskVectorQty(quantity, row, column)
          end do
        end do
      else if ( present(precisionQuantity) ) then
        do column=1, size(quantity%values(1, :))
          do row=1, size(quantity%values(:, 1))
            if ( isVectorQtyMasked(precisionQuantity, row, column) ) &
              & call MaskVectorQty(quantity, row, column)
          end do
        end do
      end if

      if ( index(switches, 'l1b') /= 0 ) &
        & call Dump( l1bData )
      call DeallocateL1BData(l1bData)

      if (toggle(gen) .and. levels(gen) > 0 ) call trace_end( "FillVectorQuantityFromL1B" )
    end subroutine FillVectorQuantityFromL1B

    ! ------------------------------------------- FillVectorQuantityFromL2AUX --
    subroutine FillVectorQuantityFromL2AUX ( qty, l2aux, errorCode )
      type ( VectorValue_T), intent(inout) :: QTY
      type ( L2AUXData_T), intent(in) :: L2AUX
      integer, intent(inout) :: ERRORCODE

      ! Executable code
      qty%values = reshape ( l2aux%values ( :, :,  &
        & qty%template%mafIndex(1)+1 : &
        & qty%template%mafIndex(qty%template%noInstances)+1 ), &
        & (/ qty%template%instanceLen, qty%template%noInstances /) )
    end subroutine FillVectorQuantityFromL2AUX

    ! --------------------------------------- FillQuantityUsingMagneticModel --
    subroutine FillQuantityUsingMagneticModel ( qty, gphQty, key )
      use Geometry, only: SecPerYear
      use IGRF_INT, only: FELDC, FELDCOF, To_Cart
      type (VectorValue_T), intent(inout) :: QTY
      type (VectorValue_T), intent(inout) :: GPHQTY
      integer, intent(in) :: KEY
      ! Local variables
      real :: B(3)                      ! Magnetic field
      integer :: INSTANCE               ! Loop counter
      integer :: SURF                   ! Loop counter
      integer :: SURFOR1                ! Index
      real    :: XYZ(3)                 ! lat, lon, height for to_cart

      ! Executable code

      if ( .not. ValidateVectorQuantity ( qty, quantityType=(/l_magneticField/), &
        & frequencyCoordinate=(/ l_xyz /) ) ) then
        call Announce_Error ( key, 0, &
          & 'Quantity does not describe magnetic field' )
        return
      end if
      if ( .not. ValidateVectorQuantity ( gphQty, quantityType=(/l_gph/), &
        & frequencyCoordinate=(/ l_none /), verticalCoordinate=(/l_zeta/) ) ) then
        call Announce_Error ( key, 0, &
          & 'GPH quantity does not describe gph field' )
        return
      end if
      if ( .not. DoHGridsMatch ( qty, gphQty ) ) then
        call Announce_Error ( key, 0, &
          & 'Quantity and GPHQuanity do not share the same horizontal basis' )
        return
      end if
      if ( .not. DoVGridsMatch ( qty, gphQty ) ) then
        call Announce_Error ( key, 0, &
          & 'Quantity and GPHQuantity do not share the same vertical basis' )
        return
      end if

      ! Assume the time is close enough to constant that one call to
      ! FELDCOF is accurate enough.

      call feldcof ( real(qty%template%time(1,1)/secPerYear + epoch) )

      do instance = 1, qty%template%noInstances
        do surf = 1, qty%template%noSurfs
          if ( qty%template%stacked ) then
            surfOr1 = 1
          else
            surfOr1 = surf
          end if
          ! Convert (/ lat(deg), lon(deg), height(km) /) to cartesian (e-radius)
          call to_cart ( real( (/ qty%template%geodLat(surfOr1,instance), &
            &                     qty%template%lon(surfOr1,instance), &
            &                     gphQty%values(surf,instance)*1.0e-3 /) ), xyz )
          ! Compute the field at and w.r.t. cartesian coordinates
          call feldc ( xyz, b )
          qty%values ( surf*3-2 : surf*3, instance) = b
        end do
      end do

    end subroutine FillQuantityUsingMagneticModel

    ! ------------------------------------------- FillQtyFromInterpolatedQty
    subroutine FillQtyFromInterpolatedQty ( qty, source, key )
      type (VectorValue_T), intent(inout) :: QTY
      type (VectorValue_T), intent(in) :: SOURCE
      integer, intent(in) :: KEY

      ! Local variables
      real (r8), dimension(:), pointer :: oldSurfs, newSurfs
      real (r8), dimension(:,:), pointer :: newValues
      logical :: mySurfs, myNewValues

      ! Executable code
      if ( .not. DoQtysDescribeSameThing ( qty, source ) ) then
        call Announce_error ( key, no_error_code, &
          & 'Mismatch in quantities' )
        return
      end if
      if ( .not. doHGridsMatch ( qty, source ) ) then
        call Announce_error ( key, no_error_code, &
          & 'Mismatch in horizontal grid' )
        return
      end if
      if ( .not. doFGridsMatch ( qty, source ) ) then
        call Announce_error ( key, no_error_code, &
          & 'Mismatch in frequency grid' )
        return
      end if

      ! These checks are for cases the code can't (yet) handle, 
      ! may add this functionality later. 
      if ( qty%template%noChans /= 1 ) then
        call Announce_error ( key, no_error_code, &
          & 'Code cannot (yet?) interpolate multi channel quantities' )
        return
      end if
      if ( .not. all ( (/ qty%template%coherent, source%template%coherent /) ) ) then
        call Announce_error ( key, no_error_code, &
          & 'Code cannot (yet?) interpolate incoherent quantities' )
        return
      end if

      ! Work out vertical coordinate issues
      if ( qty%template%verticalCoordinate == l_pressure ) then
        nullify ( oldSurfs, newSurfs )
        call Allocate_test ( oldSurfs, source%template%noSurfs, 'oldSurfs', ModuleName )
        call Allocate_test ( newSurfs, qty%template%noSurfs, 'newSurfs', ModuleName )
        oldSurfs = -log10 ( source%template%surfs(:,1) )
        newSurfs = -log10 ( qty%template%surfs(:,1) )
        mySurfs = .true.
      else
        oldSurfs => source%template%surfs ( :, 1 )
        newSurfs => qty%template%surfs ( :, 1 )
        mySurfs = .false.
      end if

      ! Work out if we have to obey the mask
      myNewValues = .false.
      if ( associated ( qty%mask ) ) &
        & myNewValues = any ( iand ( ichar(qty%mask(:,:)), m_fill ) /= 0 )
      if ( myNewValues ) then
        nullify ( newValues )
        call Allocate_test ( newValues, qty%template%instanceLen, &
          & qty%template%noInstances, 'myNewValues', ModuleName )
      else
        newValues => qty%values
      end if
      
      ! OK, do the work
      if ( qty%template%logBasis ) then
        call InterpolateValues ( &
          & oldSurfs, log ( max ( source%values, sqrt(tiny(0.0_r8)) ) ), &
          & newSurfs, newValues, &
          & method='Linear', extrapolate='Constant' )
        newValues = exp ( newValues )
      else
        call InterpolateValues ( &
          & oldSurfs, source%values, &
          & newSurfs, newValues, &
          & method='Linear', extrapolate='Constant' )
      end if

      if ( myNewValues ) then
        where ( iand ( ichar(qty%mask(:,:)),m_fill) == 0 )
          qty%values = newValues
        end where
        call Deallocate_test ( newValues, 'myNewValues', ModuleName )
      endif

      ! Tidy up
      if ( mySurfs ) then
        call Deallocate_test ( oldSurfs, 'oldSurfs', ModuleName )
        call Deallocate_test ( newSurfs, 'newSurfs', ModuleName )
      end if

    end subroutine FillQtyFromInterpolatedQty

    !=============================== FillQuantityFromLosGrid ====
    subroutine FillQuantityFromLosGrid ( key, Qty, LOS, &
      & Ptan, Re, noFineGrid, extinction, errorCode )

      ! This is to fill a l2gp type of quantity with a los grid type of quantity.
      ! The los quantity is a vector quantity that has dimension of (s, mif, maf),
      ! where s is the path along los.

      ! Linear interpolation is used to fill l2gp grids and unfilled grids are
      ! marked with the baddata flag (-999.)

      ! Dummy arguments
      integer, intent(in) :: key          ! For messages    
      type (VectorValue_T), intent(in) :: LOS ! Vector quantity to fill from
      type (VectorValue_T), intent(in) :: Ptan ! tangent pressure
      type (VectorValue_T), intent(in) :: Re ! Earth's radius
      type (VectorValue_T), INTENT(INOUT) :: QTY ! Quantity to fill
      integer, intent(in) :: noFineGrid   ! make finer sGrid with this number
      logical, intent(in) :: extinction  ! Flag for extinction fill and calculation
      integer, intent(out) :: errorCode ! Error code

      ! Local variables
      integer :: i, j, maf, mif                ! Loop counter
      integer :: maxZ, minZ                    ! pressure range indices of sGrid
      integer :: noMAFs,noMIFs,noDepths
      integer, dimension(qty%template%noSurfs,qty%template%noInstances) :: cnt
      real (r8), dimension(qty%template%noSurfs,qty%template%noInstances) :: out
      real (r8), dimension(qty%template%noSurfs) :: outZeta, phi_out, beta_out
      real (r8), dimension(los%template%noChans) :: x_in, y_in, sLevel
      real (r8), dimension(los%template%noSurfs) :: zt
      real (r8), dimension(los%template%noChans*noFineGrid) :: betaFine, TransFine, SFine
      real (r8), dimension(los%template%noChans, &
        & los%template%noSurfs,los%template%noInstances) :: beta
      real (r8) :: ds, ColTrans

      if ( toggle(gen) ) call trace_begin ( "FillQuantityFromLosGrid", key )

      errorCode=0

      ! Make sure this quantity is appropriate
      !    if (.not. ValidateVectorQuantity(qty, coherent=.TRUE., stacked=.TRUE., &
      !      & verticalCoordinate= (/ l_pressure, l_zeta /) ) ) then
      !      call output ( " quantity vertical grid in FillQuantityFromLOSgrid is not valid")
      !    end if

      if ( qty%template%verticalCoordinate == l_pressure ) then
        outZeta = -log10 ( qty%template%surfs(:,1) )
      else
        outZeta = qty%template%surfs(:,1)
      end if

      noMAFs=los%template%noInstances
      noMIFs=los%template%noSurfs
      ! Now, we use frequency coordinate as sGrid along the path
      noDepths=los%template%noChans
      sLevel = los%template%frequencies

      ! the input losQty is the increment of cloud transmission function by default.
      ! it is converted to cloud extinction if extinction flag is on.
      if(extinction) then
        ! both sGrid and sFineGrid are expected to be evenly spaced at present
        ds = sLevel(2)-sLevel(1)  
        do i=1,noDepths
          do j=1,noFineGrid
            Sfine(j+(i-1)*noFineGrid) = sLevel(i)+(j-1._r8)*ds/noFineGrid
          end do
        end do

        do maf=1,noMafs
          do mif=1,noMifs
            do i=1,noDepths
              ! convert increments to derivatives
              y_in(i) = los%values(i+(mif-1)*noDepths,maf)/ds
            end do
            call InterpolateValues(sLevel,y_in,sFine,TransFine,method='Linear')  
            ! calculate column transmission function by integrating
            ! the derivatives on fine grid
            do i=1,noFineGrid*noDepths
              betaFine(i) = 0._r8
              colTrans=0._r8
              do j=1,i
                colTrans=colTrans + transFine(j)*ds/noFineGrid
              end do
              colTrans = 1._r8 - colTrans
              if(colTrans > 0.02_r8) betaFine(i)= transFine(i)/colTrans
            end do
            ! interpolate betaFine back to the coarser sGrid
            call InterpolateValues(sFine,betaFine,sLevel,beta(:,mif,maf),method='Linear')

          end do
        end do

      else

        do maf=1,noMafs
          do mif=1,noMifs
            do i=1,noDepths
              beta(i,mif,maf)=los%values(i+(mif-1)*noDepths,maf)
            end do
          end do
        end do

      end if

      ! initialize quantity
      do j = 1, qty%template%noInstances
        do i = 1, qty%template%noSurfs
          qty%values(i,j)=qty%template%badValue
          cnt(i,j)=0
          out(i,j)=0._r8
        end do
      end do

      do maf=1,noMAFs  
        zt = ptan%values(:,maf)   ! noChans=1 for ptan
        zt = (zt+3.)*16.                      ! converted to height in km
        do mif=1,noMIFs
          if (ptan%values(mif,maf) .gt. -2.5) cycle ! for testing
          ! find altitude of each s grid
          x_in = sLevel**2/2./(re%values(1,maf)*0.001_r8 + zt(mif))
          ! converted to zeta
          x_in = x_in/16. + ptan%values(mif,maf)
          ! find minimum and maximum pressures indices in sGrid
          do i = 2,qty%template%noSurfs-1
            if (ptan%values(mif,maf) < (outZeta(i)+outZeta(i+1))/2. .and. &
              & ptan%values(mif,maf) > (outZeta(i)+outZeta(i-1))/2.) &
              & minZ = i
          end do
          if (ptan%values(mif,maf) < (outZeta(1)+outZeta(2))/2.) minZ=1
          if (ptan%values(mif,maf) > outZeta(qty%template%noSurfs)) cycle ! goto next mif

          do i = 2,qty%template%noSurfs-1
            if (x_in(noDepths) < (outZeta(i)+outZeta(i+1))/2. .and. &
              & x_in(noDepths) > (outZeta(i)+outZeta(i-1))/2.) &
              & maxZ = i
          end do
          if (x_in(noDepths) < (outZeta(1)+outZeta(2))/2.) cycle    ! goto next mif
          if (x_in(noDepths) > outZeta(qty%template%noSurfs)) maxZ=qty%template%noSurfs

          ! get phi along path for each mif (phi is in degree)
          y_in = los%template%phi(mif,maf) &
            & - atan(sLevel/(re%values(1,maf)*0.001_r8 + zt(mif)))*180._r8/Pi
          ! interpolate phi onto standard vertical grids     
          call InterpolateValues(x_in,y_in,outZeta(minZ:maxZ),phi_out(minZ:maxZ), &
            & method='Linear')
          ! interpolate quantity to standard vertical grids      
          y_in = beta(:,mif,maf)
          call InterpolateValues(x_in,y_in,outZeta(minZ:maxZ),beta_out(minZ:maxZ), &
            & method='Linear')
          ! interpolate quantity to standard phi grids
          do i=minZ,maxZ  
            do j = 2, qty%template%noInstances-1
              if(phi_out(i) .lt. &     
                & (qty%template%phi(1,j)+qty%template%phi(1,j+1))/2. &
                & .and. phi_out(i) .ge. &  
                & (qty%template%phi(1,j-1)+qty%template%phi(1,j))/2. ) then
                out(i,j)=out(i,j) + beta_out(i)
                cnt(i,j)=cnt(i,j)+1       !  counter
              end if
            end do
          end do
        end do                            ! End surface loop
      end do                              ! End instance loop
      ! average all non-zero bins
      where (cnt > 0) qty%values = out/cnt

    end subroutine FillQuantityFromLosGrid

    ! --------------------------------------------- FillQuantityByManipulation ---
    subroutine FillQuantityByManipulation ( quantity, a, b, manipulation, key )
      use String_table, only: GET_STRING
      type (VectorValue_T), intent(inout) :: QUANTITY
      type (VectorValue_T), intent(inout) :: A
      type (VectorValue_T), pointer :: B
      integer, intent(in) :: MANIPULATION
      integer, intent(in) :: KEY        ! Tree node
      ! Local variables
      character (len=128) :: MSTR
      ! Executable code
      call get_string ( manipulation, mstr, strip=.true. )
      if ( mstr /= 'a+b' .and. mstr /= 'a-b' ) then
        call Announce_Error ( key, 0, &
          & 'Only a+b or a-b allowed for manipulation at the moment' )
        return
      end if
      ! Check that this operation makes sense. Most of the time this means that
      ! the quantities have the same template. In some cases however, we can be
      ! a little more lenient
      if ( a%template%name /= quantity%template%name .and. .not. ( &
          &  ( a%template%minorFrame .and. quantity%template%minorFrame ) .and. &
          &  ( a%template%signal == quantity%template%signal ) .and. &
          &  ( a%template%sideband == quantity%template%sideband ) .and. &
          &  ( a%template%frequencyCoordinate == &
          &      quantity%template%frequencyCoordinate ) ) ) then
        call Announce_Error ( key, 0, &
          & 'a is not of the same (or close enough) type as quantity' )
        return
      end if
      if ( associated ( b ) ) then
        if ( b%template%name /= quantity%template%name .and. .not. ( &
          &  ( b%template%minorFrame .and. quantity%template%minorFrame ) .and. &
          &  ( b%template%signal == quantity%template%signal ) .and. &
          &  ( b%template%sideband == quantity%template%sideband ) .and. &
          &  ( b%template%frequencyCoordinate == &
          &      quantity%template%frequencyCoordinate ) ) ) then
          call Announce_Error ( key, 0, &
            & 'b is not of the same (or close enough) type as quantity' )
          return
        end if
      else
        ! Later we'll be more tolerant of this
        call Announce_Error ( key, 0, &
          & 'You did not supply a b quantity' )
        return
      end if
      ! OK do the simple work for now
      ! Later we'll do fancy stuff to parse the manipulation.
      if(mstr .eq. 'a+b') then
        if ( .not. associated ( quantity%mask ) ) then
          quantity%values = a%values + b%values
        else
          where ( iand ( ichar(quantity%mask(:,:)), m_fill ) == 0 )
            quantity%values = a%values + b%values
          end where
        end if
      end if
      if(mstr .eq. 'a-b') then
        if ( .not. associated ( quantity%mask ) ) then
          quantity%values = a%values - b%values
        else
          where ( iand ( ichar(quantity%mask(:,:)), m_fill ) == 0 )
            quantity%values = a%values - b%values
          end where
        end if
      end if
    end subroutine FillQuantityByManipulation

    ! ----------------------------------------------- OffsetRadianceQuantity ---
    subroutine OffsetRadianceQuantity ( quantity, radianceQuantity, amount )
      type (VectorValue_T), intent(inout) :: QUANTITY
      type (VectorValue_T), intent(in) :: RADIANCEQUANTITY
      real (rv), intent(in) :: AMOUNT

      ! Executable code
      if ( .not. ValidateVectorQuantity ( quantity, &
        & quantityType=(/l_radiance/) ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Quantity for offsetRadiance fill is not radiance' )
      if ( .not. ValidateVectorQuantity ( quantity, &
        & quantityType=(/l_radiance/) ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Radiance quantity for offsetRadiance fill is not radiance' )
      if ( quantity%template%signal /= radianceQuantity%template%signal .or. &
        & quantity%template%sideband /= radianceQuantity%template%sideband ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Quantity and rad. qty. in offsetRadiance fill different signal/sideband' )
      if ( .not. associated ( radianceQuantity%mask ) ) return
      where ( iand ( ichar(radianceQuantity%mask), m_linAlg ) == 1 )
        quantity%values = quantity%values + amount
      end where
    end subroutine OffsetRadianceQuantity

    ! ---------------------------------------------- TRANSFERVECTORS -----
    subroutine TransferVectors ( source, dest, skipMask )
      ! Copy common items in source to those in dest
      type (Vector_T), intent(in) :: SOURCE
      type (Vector_T), intent(inout) :: DEST
      logical, intent(in) :: SKIPMASK

      ! Local variables
      type (VectorValue_T), pointer :: DQ ! Destination quantity
      type (VectorValue_T), pointer :: SQ ! Source quantity
      integer :: SQI                      ! Quantity index in source

      ! Executable code

      ! First copy those things in source, loop over them
      dest%globalUnit = source%globalUnit
      do sqi = 1, size ( source%quantities )
        ! Try to find this in dest
        sq => source%quantities(sqi)
        dq => GetVectorQtyByTemplateIndex ( dest, source%template%quantities(sqi) )
        if ( associated ( dq ) ) then
          dq%values = sq%values
          if ( .not. skipMask ) then
            if (associated(sq%mask)) then
              if (.not. associated(dq%mask)) call CreateMask ( dq )
              dq%mask = sq%mask
            else
              if ( associated(dq%mask) ) &
                & call Deallocate_test ( dq%mask, 'dq%mask', ModuleName )
            end if
          end if
        end if
      end do
    end subroutine TransferVectors

    ! ---------------------------------------------  ANNOUNCE_ERROR  -----
    subroutine ANNOUNCE_ERROR ( where, CODE , ExtraMessage, ExtraInfo )
      integer, intent(in) :: where   ! Tree node where error was noticed
      integer, intent(in) :: CODE    ! Code for error message
      character (LEN=*), intent(in), optional :: ExtraMessage
      integer, intent(in), dimension(:), optional :: ExtraInfo

      integer :: I

      error = max(error,1)
      call output ( '***** At ' )

      if ( where > 0 ) then
        call print_source ( source_ref(where) )
      else
        call output ( '(no lcf tree available)' )
      end if

      call output ( ': The' );

      select case ( code )
      case ( allocation_err )
        call output ( " command caused an allocation error in squeeze.", advance='yes' )
      case ( badEstNoiseFill )
        call output ( " missing information for estimated noise fill", advance='yes' )
      case ( badGeocAltitudeQuantity )
        call output ( " geocAltitudeQuantity is not geocAltitude", advance='yes' )
      case ( badlosGridfill )
        call output ( " incomplete/incorrect information for los Grid fill", advance='yes' )
      case ( badlosvelfill )
        call output ( " incomplete/incorrect information for los velocity", advance='yes' )
      case ( badIsotopeFill )
        call output ( " incomplete/incorrect information for isotope fill", advance='yes' )
      case ( badREFGPHQuantity )
        call output ( " refGPHQuantity is not refGPH", advance='yes' )
      case ( badTemperatureQuantity )
        call output ( " temperatureQuantity is not temperature", advance='yes' )
      case ( badUnitsForExplicit )
        call output ( " explitictValues field has inappropriate " // &
          & "units for Fill instruction.", advance='yes' )
      case ( badUnitsForIntegrationTime )
        call output ( " has inappropriate units for integration time.", advance='yes' )
      case ( badUnitsForSystemTemperature )
        call output ( " has inappropriate units for system temperature.", advance='yes' )
      case ( badUnitsForMaxIterations )
        call output ( " maxIterations should be dimensionless", advance='yes' )
      case ( bothFractionAndLength )
        call output ( " specified both fraction and lengthScale", advance='yes' )
      case ( cantFillFromL1B )
        call output ( " command could not be filled from L1B.", advance='yes' )
      case ( cantFillFromL2AUX )
        call output ( " command could not be filled from L2AUX.", advance='yes' )
      case ( cantInterpolate3D )
        call output ( " cannot interpolate 3d quantities (yet).", advance='yes' )
      case ( deallocation_err )
        call output ( " command caused an deallocation error in squeeze.", advance='yes' )
      case ( errorReadingL1B )
        call output ( " L1B file could not be read.", advance='yes' )
      case ( invalidExplicitFill )
        call output ( " has inappropriate dimensionality for explicit fill.", advance='yes' )
      case ( m1_too_small )
        call output ( " command caused a m1 too small error in squeeze.", advance='yes' )
      case ( m2_too_small )
        call output ( " command caused a m2 too small error in squeeze.", advance='yes' )
      case ( missingField )
        call output ( " fields " )
        do i = 1, size(extraInfo)
          call display_string ( field_indices(i) )
          if ( i == size(extraInfo) ) then
            call output ( " and " )
          else
            call output ( ", " )
          end if
        end do
        call output ( " are required.", advance='yes' )
      case ( needGeocAltitude )
        call output ( " needs geocAltitudeQuantity.", advance='yes' )
      case ( needGeodAltitude )
        call output ( " vertical coordinate should be geoditic altitude.", &
          & advance='yes' )
      case ( needH2O )
        call output ( " needs H2OQuantity.", advance='yes' )
      case ( needOrbitInclination )
        call output ( " needs OrbitalInclination.", advance='yes' )
      case ( needTempREFGPH )
        call output ( " needs temperatureQuantity and refGPHquantity.", advance='yes' )
      case ( noExplicitValuesGiven )
        call output ( " no explicit values given for explicit fill.", advance='yes' )
      case ( nonConformingHydrostatic )
        call output ( " quantities needed for hydrostatic fill do not conform", advance='yes' )
      case ( noSourceGridGiven )
        call output ( " no sourceGrid field given for gridded fill.", advance='yes' )
      case ( noSourceL2AUXGiven )
        call output ( " no sourceL2AUX field given for L2AUX fill.", advance='yes' )
      case ( noSourceL2GPGiven )
        call output ( " no sourceL2GP field given for L2GP fill.", advance='yes' )
      case ( noSourceQuantityGiven )
        call output ( " no sourceQuantity field given for vector fill.", advance='yes' )
      case ( noSpecialFill )
        call output ( " invalid special fill", advance='yes' )
      case ( notImplemented )
        call output ( extraMessage )
        call output ( " is not implemented yet.", advance='yes' )
      case ( notSPD )
        call output ( " is not a SPD matrix.", advance='yes' )
      case ( not_permutation )
        call output ( " command caused an illegal permutation in squeeze.", advance='yes' )
      case ( numInstancesisZero )
        call output ( " command has zero instances.", advance='yes' )
      case ( numSurfsisZero )
        call output ( " command has zero surfaces.", advance='yes' )
      case ( n1_is_zero )
        call output ( " command caused an n1=0 error in squeeze.", advance='yes' )
      case ( n2_is_zero )
        call output ( " command caused an n2=0 error in squeeze.", advance='yes' )
      case ( n3_is_zero )
        call output ( " command caused an n3=0 error in squeeze.", advance='yes' )
      case ( objIsFullRank3 )
        call output ( " command array is full rank 3.", advance='yes' )
      case ( otherErrorInFillVector )
        call output ( " command caused an error in FillVector.", advance='yes' )
      case ( source_not_in_db )
        call output ( " source was not found in the db.", advance='yes' )
      case ( unknownQuantityName )
        call output ( " quantity was not found in the vector", advance='yes' )
      case ( vectorWontMatchPDef )
        call output ( " command found new and prev. vectors unmatched.", advance='yes' )
      case ( wrong_number )
        call output ( " command does not have exactly one field.", advance='yes' )
      case ( zeroGeodSpan )
        call output ( " command found zero geod. ang. span.", advance='yes' )
      case ( zeroProfilesFound )
        call output ( " command found zero profiles.", advance='yes' )
      case ( badRefractFill )
        call output ( " missing information for phiTan refract fill", advance='yes' )
      case default
        call output ( " command caused an unrecognized programming error", advance='yes' )
      end select
      if ( present(ExtraMessage) ) then
        call output(ExtraMessage, advance='yes')
      end if
    end subroutine ANNOUNCE_ERROR
  end subroutine MLSL2Fill

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Fill
!=============================================================================

!
! $Log$
! Revision 2.184  2003/02/18 23:59:06  livesey
! Added phiWindow for hydrostatic fill.
!
! Revision 2.183  2003/02/15 00:35:03  livesey
! Added error checking for range of profile on L2GP fill.
!
! Revision 2.182  2003/02/13 21:42:12  livesey
! Added specific profile stuff to fill from l2gp.
!
! Revision 2.181  2003/01/29 01:59:19  livesey
! Changed some MLSMessages to Announce_Errors to get a line number.
!
! Revision 2.180  2003/01/28 21:53:07  pwagner
! RHI H2O conversions moved to fwdmdl from l2/Fill
!
! Revision 2.179  2003/01/16 21:48:58  vsnyder
! Fix a comment
!
! Revision 2.178  2003/01/15 23:29:41  pwagner
! Compatible with adjustable data types L2GPData
!
! Revision 2.177  2003/01/15 02:49:06  vsnyder
! Get SecPerYear from Geometry module
!
! Revision 2.176  2003/01/14 23:53:10  livesey
! Bug fix in magnetic model
!
! Revision 2.175  2003/01/14 22:39:25  livesey
! Bug fixes in magnetic stuff
!
! Revision 2.174  2003/01/14 21:58:59  vsnyder
! OOPS, Left out 'template' in reference to verticalCoordinate
!
! Revision 2.173  2003/01/14 21:34:09  vsnyder
! More work on magnetic field vector quantity
!
! Revision 2.172  2003/01/12 07:34:05  dwu
! with some fix for a-b manipulation
!
! Revision 2.171  2003/01/12 05:13:50  dwu
! add a-b manipulation (under the same conditions of the a+b case
!
! Revision 2.170  2003/01/08 23:52:16  livesey
! Bug fix in offset radiance
!
! Revision 2.169  2003/01/07 23:46:38  livesey
! Added magentic model
!
! Revision 2.168  2002/11/29 22:46:15  livesey
! Tidyup on l2aux fill
!
! Revision 2.167  2002/11/27 22:59:21  livesey
! Made the checking in FillQuantityByManipulation a little more lenient
!
! Revision 2.166  2002/11/27 22:18:10  dwu
! Change the error handling in the new manipulation feature. Instead of quitting, just send off a warning
!
! Revision 2.165  2002/11/27 19:25:45  livesey
! Added manipulate method and bug fix in snooper when no matrices
!
! Revision 2.164  2002/11/21 01:18:11  livesey
! Added negativePrecision command (as distinct from fill method of same
! name).
!
! Revision 2.163  2002/11/14 17:28:01  livesey
! Made the profile fill do a 'nearest' on the input heights in the case of
! coherent quantities.
!
! Revision 2.162  2002/11/13 01:06:42  pwagner
! Fixed small bug
!
! Revision 2.161  2002/11/06 02:01:05  livesey
! Changes to fill from l2aux
!
! Revision 2.160  2002/10/26 00:02:51  livesey
! Another typo! Going too fast!
!
! Revision 2.159  2002/10/26 00:00:07  livesey
! Typo
!
! Revision 2.158  2002/10/25 23:59:45  livesey
! Forgot to include parameter f_offsetamount
!
! Revision 2.157  2002/10/25 23:56:14  livesey
! Added the offsetAmount argument default 1000K
!
! Revision 2.156  2002/10/17 18:18:50  livesey
! Added low/high bound to vector creation
!
! Revision 2.155  2002/10/16 20:15:27  mjf
! Added GPH precision.
!
! Revision 2.154  2002/10/10 23:52:57  pwagner
! Added code to fill from L1bdata with hdf5; untested
!
! Revision 2.153  2002/10/08 17:36:20  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.152  2002/10/03 13:44:04  mjf
! Renamed temperaturePrecisionQuantitiy to tempPrecisionQuantity to get
! names < 31 long.
!
! Revision 2.151  2002/10/02 23:04:47  pwagner
! RHI now separate fill methods
!
! Revision 2.150  2002/10/01 18:28:04  mjf
! New Fill for RHi precision including T error.
!
! Revision 2.149  2002/09/26 20:41:04  vsnyder
! Get Omega from Geometry instead of Units -- that's where it really belongs
!
! Revision 2.148  2002/09/25 20:08:14  livesey
! Made -g less verbose
!
! Revision 2.147  2002/09/13 18:10:10  pwagner
! May change matrix precision rm from r8
!
! Revision 2.146  2002/09/12 22:07:05  livesey
! Added masking to interpolated vector fill
!
! Revision 2.145  2002/09/10 20:50:33  livesey
! Added interpolated vector fill
!
! Revision 2.144  2002/09/10 01:00:32  livesey
! Added calls to NullifyMatrix
!
! Revision 2.143  2002/09/06 00:54:13  livesey
! Tiny change
!
! Revision 2.142  2002/09/05 21:48:54  livesey
! More fixes in transfer skipmask
!
! Revision 2.141  2002/09/05 20:49:38  livesey
! Added skipMask to transfer
!
! Revision 2.140  2002/08/28 01:13:52  livesey
! Added OffsetRadianceQuantity
!
! Revision 2.139  2002/08/26 20:01:09  livesey
! Added instances argument to profile fill
!
! Revision 2.138  2002/08/21 23:06:12  livesey
! Bug fix in profile fill
!
! Revision 2.137  2002/08/20 22:10:49  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.136  2002/08/20 19:30:24  livesey
! Typo!
!
! Revision 2.135  2002/08/20 19:19:41  livesey
! Tidied up FillVectorQtyFromProfile
!
! Revision 2.134  2002/08/16 16:08:58  livesey
! Bug fix in the Matrix fill (covariance with frequency variation).
!
! Revision 2.133  2002/08/15 03:52:52  livesey
! Added profile fill, still not complete though.
!
! Revision 2.132  2002/08/03 20:41:40  livesey
! Added snooping of matrices
!
! Revision 2.131  2002/06/26 01:26:21  livesey
! Added 2D pressure guesser
!
! Revision 2.130  2002/06/18 20:00:36  livesey
! Removed unwanted print
!
! Revision 2.129  2002/06/14 16:40:00  livesey
! Orbital inclination can now be filled from l1b
!
! Revision 2.128  2002/06/04 23:22:36  livesey
! Bug fixes on phiTan fill, other cleanups
!
! Revision 2.127  2002/06/04 22:40:44  livesey
! Added framework for phiTan fill
!
! Revision 2.126  2002/05/28 17:08:42  livesey
! More explicit error message still in l2gp fill.
!
! Revision 2.125  2002/05/23 20:51:53  livesey
! Bug fix, checking wrong in special fill for los velocity.
!
! Revision 2.124  2002/05/17 17:55:48  livesey
! Added sideband folding fill
!
! Revision 2.123  2002/05/14 00:26:25  livesey
! New code for system temperatures etc.
!
! Revision 2.122  2002/05/06 22:30:51  livesey
! Tidied up l2gp fill.
!
! Revision 2.121  2002/04/25 20:47:02  livesey
! Added channel dependent system temperature, and removed
! embarassing bug whereby radiance noise was sqrt(10) too small!
!
! Revision 2.120  2002/04/18 20:14:52  pwagner
! Fills either rhi from h2o or inverse; passes non-interpolating test
!
! Revision 2.119  2002/04/16 23:27:43  pwagner
! FillRHI testing begun; incomplete
!
! Revision 2.118  2002/04/13 00:31:46  pwagner
! More flesh on FillrhiFromH2o; still untested
!
! Revision 2.117  2002/04/11 23:51:28  pwagner
! Fleshed out FillRHIFromH2O; untested yet
!
! Revision 2.116  2002/04/10 17:45:44  pwagner
! Added FillRHI from h2oquantity (just a placeholder)
!
! Revision 2.115  2002/04/04 16:32:42  livesey
! Added negative error bar stuff
!
! Revision 2.114  2002/03/27 17:37:57  livesey
! Minor changes to random number stuff.
! Now seed incremented with chunk number
!
! Revision 2.113  2002/03/19 00:52:40  pwagner
! SOme new checks added to FillLOSVelocity
!
! Revision 2.112  2002/03/14 17:29:59  pwagner
! Fixed check in FillLOSVelocity
!
! Revision 2.111  2002/03/14 01:01:17  pwagner
! Can fill scVelECI and scVelECR from l1b
!
! Revision 2.110  2002/03/13 22:01:41  livesey
! Added masking of fill from vector quantity, also changed
! from m_explicitfill to m_fill
!
! Revision 2.109  2002/03/08 08:07:00  livesey
! Added explicit fill mask
!
! Revision 2.105  2002/02/06 01:35:29  livesey
! Added ability to load ptan from l1b
!
! Revision 2.104  2002/02/05 01:45:21  livesey
! Added preliminary vector/vector fill, not yet tested.
!
! Revision 2.103  2002/01/18 00:24:21  livesey
! Added sideband argument to call to GetSignalName when reading
! L1B radiances (so can now read split signals).
!
! Revision 2.102  2002/01/16 18:06:31  livesey
! Changed neverFail flag in ReadL1B to false.  I want to see
! the error message.
!
! Revision 2.101  2002/01/09 00:00:04  pwagner
! Replaced write or print statements with calls to output
!
! Revision 2.100  2001/12/06 23:45:07  livesey
! Dealt with move of Omega to Units
!
! Revision 2.99  2001/11/09 23:17:22  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.98  2001/10/26 23:23:05  pwagner
! Complies with l1b data dump
!
! Revision 2.97  2001/10/26 18:15:18  livesey
! Added upper argument to MatrixInversion in FillCovariance
!
! Revision 2.96  2001/10/25 23:32:11  pwagner
! Responds to l1b switch by dumping l1b quantity during Fill
!
! Revision 2.95  2001/10/24 22:35:33  dwu
! add FillDiagonal
!
! Revision 2.94  2001/10/23 16:38:09  pwagner
! Fill from l1b can fill precision, set mask
!
! Revision 2.93  2001/10/19 23:42:50  pwagner
! Applies multipliers to chi^2 fills, too
!
! Revision 2.92  2001/10/19 22:32:44  pwagner
! Can destroy a vector or a matrix; def. multiplier(s) in addNoise
!
! Revision 2.91  2001/10/19 00:00:36  pwagner
! Replaced remove with destroy
!
! Revision 2.90  2001/10/18 23:29:57  pwagner
! Fixes in addNoise, remove
!
! Revision 2.89  2001/10/18 23:07:40  livesey
! Removed debug dump statement
!
! Revision 2.88  2001/10/18 23:06:05  livesey
! Tidied up some bugs in transfer, must have been asleep or something!
!
! Revision 2.87  2001/10/18 22:30:30  livesey
! Added s_dump and more functionality to fillCovariance
!
! Revision 2.86  2001/10/18 03:51:46  livesey
! Just some tidying up.
!
! Revision 2.85  2001/10/18 00:46:32  livesey
! Bug fixes in Transfer
!
! Revision 2.84  2001/10/17 23:40:08  pwagner
! Exploits visibility of MATH77_ran_pack
!
! Revision 2.83  2001/10/16 23:34:05  pwagner
! intrinsic, resetseed, seed fields added to addnoise method
!
! Revision 2.82  2001/10/16 00:07:48  livesey
! Got smoothing working.
!
! Revision 2.81  2001/10/15 22:10:42  livesey
! Interim version with smoothing stubbed
!
! Revision 2.80  2001/10/02 23:12:50  pwagner
! More chi^2 fixes
!
! Revision 2.79  2001/09/28 23:59:20  pwagner
! Fixed various timing problems
!
! Revision 2.78  2001/09/28 17:50:30  pwagner
! MLSL2Timings module keeps timing info
!
! Revision 2.77  2001/09/24 17:28:15  pwagner
! Gets drang from MLSRandomNumber
!
! Revision 2.76  2001/09/21 23:23:35  pwagner
! Stiff fails to add noise properly
!
! Revision 2.75  2001/09/20 20:57:25  pwagner
! Fleshed out adding noise thing
!
! Revision 2.74  2001/09/19 23:42:29  pwagner
! New Remove command, ignore() fields
!
! Revision 2.73  2001/09/18 23:53:08  pwagner
! Replaced error field name with noise; began addNoise Fill method
!
! Revision 2.72  2001/09/17 23:12:21  pwagner
! Fleshed out FillChi.. some more
!
! Revision 2.71  2001/09/14 23:34:13  pwagner
! Now should allow special fill of chi^2..
!
! Revision 2.70  2001/08/03 23:13:52  pwagner
! Began testing; at least now exits normally again
!
! Revision 2.69  2001/08/02 00:17:06  pwagner
! Mostly done with column fill; untested
!
! Revision 2.68  2001/08/01 00:05:25  dwu
! remove f_sourceSGrid
!
! Revision 2.67  2001/07/31 23:53:37  dwu
! remove sGrid source
!
! Revision 2.66  2001/07/31 23:24:17  pwagner
! column abundance calculation more fleshed out--not tested
!
! Revision 2.65  2001/07/30 23:28:38  pwagner
! Added columnAbundances scaffolding--needs fleshing out
!
! Revision 2.64  2001/07/26 20:33:40  vsnyder
! Eliminate the 'extra' field of the 'matrix' spec
!
! Revision 2.63  2001/07/20 20:03:30  dwu
! fix problems in cloud extinction calculation
!
! Revision 2.62  2001/07/20 19:25:03  dwu
! add cloud extinction calculation
!
! Revision 2.61  2001/07/19 21:45:33  dwu
! some fixes for FillQuantityFromLOS
!
! Revision 2.60  2001/07/19 18:05:42  dwu
! add sourceSGRID
!
! Revision 2.59  2001/07/19 00:56:27  dwu
! fix bugs in FillQuantityFromLos
!
! Revision 2.58  2001/07/19 00:19:42  dwu
! add new method=rectanglefromlos
!
! Revision 2.57  2001/06/22 05:37:27  livesey
! First version of transfer command
!
! Revision 2.56  2001/06/13 20:41:35  vsnyder
! Make sure 'invert' has a value
!
! Revision 2.55  2001/05/30 23:56:39  livesey
! Changed for new L1BData
!
! Revision 2.54  2001/05/30 20:16:26  vsnyder
! Add 'invert' field to 'fillCovariance' spec
!
! Revision 2.53  2001/05/23 04:38:16  livesey
! Changes chunks to pointer rather than intent(in), so it gets the right indices
!
! Revision 2.52  2001/05/19 00:15:39  livesey
! OK, that should have been square!!! (fool!)
!
! Revision 2.51  2001/05/19 00:13:41  livesey
! Made fillcovariance apply square root option to update diagonal
!
! Revision 2.50  2001/05/18 22:39:45  livesey
! Bug fix.
!
! Revision 2.49  2001/05/18 19:50:50  livesey
! Added interpolate option for l2gp fills
!
! Revision 2.48  2001/05/16 19:44:16  livesey
! Added estimated noise stuff
!
! Revision 2.47  2001/05/10 23:25:12  livesey
! Added isotope scaling stuff, tidied up some old code.
!
! Revision 2.46  2001/05/08 20:34:26  vsnyder
! Cosmetic changes
!
! Revision 2.45  2001/05/03 20:30:29  vsnyder
! Add a 'nullify' and some cosmetic changes
!
! Revision 2.44  2001/04/28 01:43:21  vsnyder
! Improved the timing message
!
! Revision 2.43  2001/04/26 02:44:17  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.42  2001/04/24 23:12:02  livesey
! Made spread more flexible
!
! Revision 2.41  2001/04/23 23:26:05  livesey
! Removed some unnecessary logic
!
! Revision 2.40  2001/04/20 17:12:24  livesey
! Add fill from l2gp
!
! Revision 2.39  2001/04/19 00:09:34  pwagner
! Longer error messages; halts if problem in fill
!
! Revision 2.38  2001/04/10 23:45:17  vsnyder
! Construct matrix properly
!
! Revision 2.37  2001/04/10 20:04:17  livesey
! Now does fill from Grid.
!
! Revision 2.36  2001/04/10 00:02:19  vsnyder
! Implement 'matrix' spec in Fill section
!
! Revision 2.35  2001/04/07 00:12:05  pwagner
! Calls trace_end if needed before every return
!
! Revision 2.34  2001/04/05 23:45:39  pwagner
! Deleted all MLSMessages
!
! Revision 2.33  2001/03/29 19:12:40  livesey
! Added gridded data fill
!
! Revision 2.32  2001/03/15 23:28:23  livesey
! Bug fix
!
! Revision 2.31  2001/03/15 21:18:57  vsnyder
! Use Get_Spec_ID instead of decoration(subtree...
!
! Revision 2.30  2001/03/15 21:12:11  livesey
! Added special fill for losVel, and dealt with new MLSSignals_m
!
! Revision 2.29  2001/03/15 18:40:38  livesey
! Added some more l1b fill options.
!
! Revision 2.28  2001/03/14 05:33:39  livesey
! Added snoop option
!
! Revision 2.27  2001/03/07 22:42:23  livesey
! Got pressure guesser to work
!
! Revision 2.26  2001/03/06 22:41:07  livesey
! New L2AUX stuff
!
! Revision 2.25  2001/03/06 00:34:46  livesey
! Regular commit.
!
! Revision 2.24  2001/03/05 01:20:14  livesey
! Regular commit, hydrostatic stuff in place.
!
! Revision 2.23  2001/03/03 05:54:29  livesey
! Started hydrostic stuff
!
! Revision 2.22  2001/03/03 00:10:14  livesey
! Removed debuging dump.
!
! Revision 2.21  2001/03/03 00:07:40  livesey
! Added fill from l1b
!
! Revision 2.20  2001/02/27 17:39:03  livesey
! Tidied stuff up a bit.
!
! Revision 2.19  2001/02/27 01:25:15  livesey
! Used new ValidateVectorQuantity routine
!
! Revision 2.18  2001/02/27 00:50:53  livesey
! Made sure verticalCoordinate=l_zeta worked for filling from L2GP
!
! Revision 2.17  2001/02/23 18:16:26  livesey
! Regular commit
!
! Revision 2.16  2001/02/21 01:07:34  livesey
! Got the explicit fill to work.
!
! Revision 2.15  2001/02/08 01:17:41  vsnyder
! Simplify access to abstract syntax tree.
!
! Revision 2.14  2001/01/26 00:11:12  pwagner
! Can fill from prev. defd. vector
!
! Revision 2.13  2001/01/24 23:31:00  pwagner
! Using announce_error, simplified
!
! Revision 2.12  2001/01/10 21:47:45  pwagner
! Chunk bounds determined by geodet.ang.
!
! Revision 2.11  2001/01/03 18:15:13  pwagner
! Changed types of t1, t2 to real
!
! Revision 2.10  2001/01/03 17:49:49  pwagner
! Accounts for chunking when filling from old L2GP
!
! Revision 2.9  2000/12/07 00:41:46  pwagner
! added whatquantitynumber
!
! Revision 2.8  2000/12/06 00:01:20  pwagner
! Completed FillOL2AUXData; changed squeeze, nearby
!
! Revision 2.7  2000/12/05 00:40:50  pwagner
! Added FillOL2AUXVector
!
! Revision 2.6  2000/11/30 00:22:52  pwagner
! functions properly moved to read a priori
!
! Revision 2.5  2000/11/16 02:15:25  vsnyder
! Implement timing.
!
! Revision 2.4  2000/11/13 23:02:21  pwagner
! Adapted for rank2 vectorsModule
!
! Revision 2.3  2000/10/06 22:18:47  pwagner
! Fills from old l2gp data
!
! Revision 2.2  2000/09/11 19:52:51  ahanzel
! Removed old log entries in file.
!
! Revision 2.1  2000/09/08 22:55:56  vsnyder
! Revised to use the tree output by the parser
!
!
