! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module Fill                     ! Create vectors and fill them.
  !=============================================================================

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use Expr_M, only: EXPR
  use GriddedData, only: GriddedData_T
  ! We need many things from Init_Tables_Module.  First the fields:
  use INIT_TABLES_MODULE, only: F_BOUNDARYPRESSURE, &
    & F_COLUMNS, F_DESTINATION, F_DIAGONAL, F_NOISE, &
    & F_GEOCALTITUDEQUANTITY, F_EARTHRADIUS, F_EXPLICITVALUES, &
    & F_EXTINCTION, F_FRACTION, F_H2OQUANTITY, F_LOSQTY,&
    & F_dontMask, F_IGNORENEGATIVE, F_IGNOREZERO, &
    & F_INTEGRATIONTIME, F_INTERPOLATE, F_INVERT, F_INTRINSIC, F_LENGTHSCALE, &
    & F_MATRIX, F_MAXITERATIONS, F_METHOD, &
    & F_MEASUREMENTS, F_MODEL, F_NOFINEGRID, F_PTANQUANTITY, &
    & F_QUANTITY, F_RADIANCEQUANTITY, F_RATIOQUANTITY, F_REFGPHQUANTITY, &
    & F_RESETSEED, F_Rows, F_SCECI, F_SCVEL, F_SEED, F_SOURCE, F_SOURCEGRID, &
    & F_SOURCEL2AUX, F_SOURCEL2GP, F_SOURCEQUANTITY, F_SOURCEVGRID, &
    & F_SPREAD, F_SUPERDIAGONAL, &
    & F_SYSTEMTEMPERATURE, F_TEMPERATUREQUANTITY, F_TEMPLATE, F_TNGTECI, &
    & F_TYPE, F_VECTOR, F_VMRQUANTITY, FIELD_FIRST, FIELD_LAST
  ! Now the literals:
  use INIT_TABLES_MODULE, only: L_ADDNOISE, L_BOUNDARYPRESSURE, L_CHISQCHAN, &
    & L_CHISQMMAF, L_CHISQMMIF, L_CHOLESKY, &
    & L_COLUMNABUNDANCE, L_ESTIMATEDNOISE, L_EXPLICIT, L_GPH, L_GRIDDED, L_HEIGHT, &
    & L_HYDROSTATIC, L_ISOTOPE, L_ISOTOPERATIO, L_KRONECKER, L_L1B, L_L2GP, L_L2AUX, &
    & L_RECTANGLEFROMLOS, L_LOSVEL, L_NONE, L_PLAIN, &
    & L_PRESSURE, L_PTAN, L_RADIANCE, &
    & L_REFGPH, L_SCECI, L_SCGEOCALT, L_SCVEL, &
    & L_SPD, L_SPECIAL, L_TEMPERATURE, L_TNGTECI, L_TNGTGEODALT, &
    & L_TNGTGEOCALT, L_TRUE, L_VGRID, L_VMR, L_ZETA
  ! Now the specifications:
  use INIT_TABLES_MODULE, only: S_DUMP, S_FILL, S_FILLCOVARIANCE, S_MATRIX, S_REMOVE, &
    & S_SNOOP, S_TIME, S_TRANSFER, S_VECTOR
  ! Now some arrays
  use Intrinsic, only: Field_Indices
  use Intrinsic, only: &
    & PHYQ_Dimensionless, PHYQ_Invalid, PHYQ_Temperature, &
    & PHYQ_Time, PHYQ_Length
  use L1BData, only: DeallocateL1BData, FindL1BData, L1BData_T, ReadL1BData
  use L2GPData, only: L2GPData_T
  use L2AUXData, only: L2AUXData_T
  use L3ASCII, only: L3ASCII_INTERP_FIELD
  use LEXER_CORE, only: PRINT_SOURCE
  use ManipulateVectorQuantities, only: DOHGRIDSMATCH, DOVGRIDSMATCH
  use MatrixModule_0, only: Sparsify, MatrixInversion
  use MatrixModule_1, only: AddToMatrixDatabase, CreateEmptyMatrix, &
    & Dump, FindBlock, GetKindFromMatrixDatabase, GetFromMatrixDatabase, K_SPD, &
    & Matrix_Cholesky_T, Matrix_Database_T, Matrix_Kronecker_T, Matrix_SPD_T, &
    & Matrix_T, UpdateDiagonal
  use MLSCommon, only: L1BInfo_T, MLSChunk_T, R8
  use MLSL2Timings, only: SECTION_TIMES, TOTAL_TIMES
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use MLSNumerics, only: InterpolateValues
  use MLSRandomNumber, only: drang, mls_random_seed, MATH77_RAN_PACK
  use MLSSignals_m, only: GetSignalName, GetModuleName
  use Molecules, only: L_H2O
  use MoreTree, only: Get_Boolean, Get_Field_ID, Get_Spec_ID
  use OUTPUT_M, only: BLANKS, OUTPUT
  use QuantityTemplates, only: QuantityTemplate_T
  use ScanModelModule, only: GetBasisGPH, GetHydrostaticTangentPressure, OMEGA
  use SnoopMLSL2, only: SNOOP
  use String_Table, only: Display_String
  use TOGGLES, only: GEN, LEVELS, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, &
    & SOURCE_REF, SUB_ROSA, SUBTREE
  use TREE_TYPES, only: N_NAMED, N_SET_ONE
  use UNITS
  use VectorsModule, only: AddVectorToDatabase, &
    & ClearUnderMask, CopyVector, CreateMask, CreateVector, &
    & DestroyVectorInfo, Dump, &
    & GetVectorQtyByTemplateIndex, isVectorQtyMasked, &
    & rmVectorFromDatabase, ValidateVectorQuantity, Vector_T, &
    & VectorTemplate_T, VectorValue_T
  use VGridsDatabase, only: VGRID_T

  implicit none
  private
  public :: MLSL2Fill

  ! -----     Private declarations     ---------------------------------

  integer, private :: ERROR
  logical, parameter :: DEEBUG = .FALSE.                 ! Usually FALSE

  ! Error codes for "announce_error"  
  integer, parameter :: No_Error_code = 0
  integer, parameter :: Wrong_Number = No_Error_code+1     ! of fields of a VECTOR command
  integer, parameter :: UnknownQuantityName = wrong_number + 1
  integer, parameter :: Source_not_in_db = unknownQuantityName + 1
  integer, parameter :: ZeroProfilesFound = source_not_in_db + 1
  integer, parameter :: ZeroGeodSpan = zeroProfilesFound + 1
  integer, parameter :: VectorWontMatchL2GP = zeroGeodSpan + 1
  integer, parameter :: CantFillFromL2AUX = vectorWontMatchL2GP + 1
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
  integer, parameter :: NeedGeocAltitude = needH2O + 1
  integer, parameter :: BadGeocAltitudeQuantity = needGeocAltitude + 1
  integer, parameter :: BadTemperatureQuantity = badGeocAltitudeQuantity + 1
  integer, parameter :: BadREFGPHQuantity = badTemperatureQuantity + 1
  integer, parameter :: NonConformingHydrostatic = badREFGPHQuantity + 1
  integer, parameter :: BadUnitsForMaxIterations = nonConformingHydrostatic + 1
  integer, parameter :: NoSpecialFill = badUnitsForMaxIterations + 1
  integer, parameter :: BadlosVelFill = noSpecialFill + 1
  integer, parameter :: NotZetaForGrid = BadLosVelFill + 1

  !  integer, parameter :: s_Fill = 0   ! to be replaced by entry in init_tables_module
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  ! This module performs the Fill operation in the Level 2 software.  
  ! This takes a vector template, and creates and fills an appropriate vector

contains ! =====     Public Procedures     =============================

  !---------------------------------------------------  MLSL2Fill  -----

  subroutine MLSL2Fill ( Root, L1bInfo, GriddedData, VectorTemplates, Vectors, &
    & QtyTemplates , Matrices, vGrids, L2GPDatabase, L2AUXDatabase, Chunks, ChunkNo )

    ! This is the main routine for the module.  It parses the relevant lines
    ! of the l2cf and works out what to do.

    ! Dummy arguments
    integer, intent(in) :: ROOT    ! Of the FILL section in the AST
    type (l1BInfo_T), intent(in) :: L1bInfo
    type (griddedData_T), dimension(:), pointer :: GriddedData
    type (vectorTemplate_T), dimension(:), pointer :: VectorTemplates
    type (vector_T), dimension(:), pointer :: Vectors
    type (quantityTemplate_T), dimension(:), pointer :: QtyTemplates
    type (matrix_database_T), dimension(:), pointer :: Matrices
    type (VGrid_T), dimension(:), intent(in) :: vGrids
    type (l2GPData_T), dimension(:), pointer :: L2GPDatabase
    type (l2AUXData_T), dimension(:), pointer :: L2AUXDatabase
    type (mlSChunk_T), dimension(:), pointer :: Chunks
    integer, intent(in) :: ChunkNo

    ! Local variables

    type (vectorValue_T), pointer :: QUANTITY ! Quantity to be filled
    type (vectorValue_T), pointer :: GEOCALTITUDEQUANTITY
    type (vectorValue_T), pointer :: H2OQUANTITY
    type (vectorValue_T), pointer :: RADIANCEQUANTITY
    type (vectorValue_T), pointer :: RATIOQUANTITY
    type (vectorValue_T), pointer :: REFGPHQUANTITY
    type (vectorValue_T), pointer :: SCECIQUANTITY
    type (vectorValue_T), pointer :: SCVELQUANTITY
    type (vectorValue_T), pointer :: SOURCEQUANTITY
    type (vectorValue_T), pointer :: TEMPERATUREQUANTITY
    type (vectorValue_T), pointer :: TNGTECIQUANTITY
    type (vectorValue_T), pointer :: TNGTPRESQUANTITY
    type (vectorValue_T), pointer :: earthRadiusQty
    type (vectorValue_T), pointer :: losQty
    type (vectorValue_T), pointer :: bndPressQty
    type (vectorValue_T), pointer :: vmrQty
    type (vectorValue_T), pointer :: measQty
    type (vectorValue_T), pointer :: modelQty
    type (vectorValue_T), pointer :: noiseQty

    integer :: bndPressQtyIndex
    integer :: bndPressVctrIndex
    integer :: ColVector                ! Vector defining columns of Matrix
    type(matrix_SPD_T), pointer :: Covariance
    integer :: DESTINATIONVECTORINDEX   ! For transfer commands
    !                                     -- for FillCovariance
    integer :: EARTHRADIUSQTYINDEX
    integer :: EARTHRADIUSVECTORINDEX
    integer :: Diagonal                 ! Index of diagonal vector in database
    !                                     -- for FillCovariance
    integer :: ERRORCODE                ! 0 unless error; returned by called routines
    logical :: Extinction               ! Flag for cloud extinction calculation
    integer :: FIELDINDEX               ! Entry in tree
    integer :: FieldValue               ! Value of a field in the L2CF
    integer :: FILLMETHOD               ! How will we fill this quantity
    integer :: FRACTION                 ! Index of fraction vector in database
    integer :: GEOCALTITUDEQUANTITYINDEX    ! In the source vector
    integer :: GEOCALTITUDEVECTORINDEX      ! In the vector database
    logical, dimension(field_first:field_last) :: GOT
    integer :: GRIDINDEX                ! Index of requested grid
    integer :: GSON                     ! Descendant of Son
    integer :: H2OQUANTITYINDEX         ! in the quantities database
    integer :: H2OVECTORINDEX           ! In the vector database
    integer :: I, J                     ! Loop indices for section, spec, expr
  !  The next three are FALSE by default
    logical :: DONTMASK                 ! Use even masked values if TRUE
    integer :: GLOBALUNIT               ! To go into the vector
    logical :: IGNOREZERO               ! Don't sum chi^2 at values of noise = 0
    logical :: IGNORENEGATIVE           ! Don't sum chi^2 at values of noise < 0
    real(r8) :: INTEGRATIONTIME         ! For estimated noise
    logical :: INTERPOLATE              ! Flag for l2gp etc. fill
    integer :: INSTANCE                 ! Loop counter
    logical :: INVERT                   ! "Invert the specified covariance matrix"
    integer :: KEY                      ! Definitely n_named
    integer :: L2AUXINDEX               ! Index into L2AUXDatabase
    integer :: L2GPINDEX                ! Index into L2GPDatabase
    integer :: LENGTHSCALE              ! Index of lengthscale vector in database
    integer :: LOSVECTORINDEX           ! index in vector database
    integer :: LOSQTYINDEX              ! index in QUANTITY database
    type(matrix_Cholesky_T) :: MatrixCholesky
    type(matrix_Kronecker_T) :: MatrixKronecker
    type(matrix_SPD_T) :: MatrixSPD
    type(matrix_T) :: MatrixPlain
    integer :: MatrixToFill             ! Index in database
    integer :: MatrixType               ! Type of matrix, L_... from init_tables
    integer :: MAXITERATIONS            ! For hydrostatic fill
    integer :: NoFineGrid               ! no of fine grids for cloud extinction calculation
    integer :: PTANVECTORINDEX          !
    integer :: PTANQTYINDEX             !
    integer :: QUANTITYINDEX            ! Within the vector
    integer :: RADIANCEQUANTITYINDEX    ! For radiance quantity
    integer :: RADIANCEVECTORINDEX    ! For radiance quantity
    integer :: RATIOQUANTITYINDEX       ! in the quantities database
    integer :: RATIOVECTORINDEX         ! In the vector database
    integer :: REFGPHQUANTITYINDEX      ! in the quantities database
    integer :: REFGPHVECTORINDEX        ! In the vector database
    logical :: ResetSeed                ! Let mls_random_seed choose new seed
    integer :: RowVector                ! Vector defining rows of Matrix
    integer :: SCECIQUANTITYINDEX       ! In the quantities database
    integer :: SCECIVECTORINDEX         ! In the vector database
    integer :: SCVELQUANTITYINDEX       ! In the quantities database
    integer :: SCVELVECTORINDEX         ! In the vector database
    integer, dimension(2) :: SEED       ! integers used by random_numbers
    integer :: SON                      ! Of root, an n_spec_args or a n_named
    integer :: SOURCEQUANTITYINDEX      ! in the quantities database
    integer :: SOURCEVECTORINDEX        ! In the vector database
    logical :: SPREAD                   ! Do we spread values accross instances in explict
    integer :: SUPERDIAGONAL            ! Index of superdiagonal matrix in database
    logical :: Switch2intrinsic         ! Have mls_random_seed call intrinsic
    !                                     -- for FillCovariance
    real(r8) :: SYSTEMTEMPERATURE       ! For estimated noise
    real :: T1, T2                      ! for timing
    integer :: TEMPERATUREQUANTITYINDEX ! in the quantities database
    integer :: TEMPERATUREVECTORINDEX   ! In the vector database
    integer :: TEMPLATEINDEX            ! In the template database
    logical :: TIMING
    integer :: TNGTECIQUANTITYINDEX     ! In the quantities database
    integer :: TNGTECIVECTORINDEX       ! In the vector database
    integer, dimension(2) :: UNITASARRAY ! From expr
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

    ! Executable code
    timing = section_times
    if ( timing ) call cpu_time ( t1 )
    dontMask = .false.
    ignoreZero = .false.
    ignoreNegative = .false.

    if ( toggle(gen) ) call trace_begin ( "MLSL2Fill", root )

    ! Logical id of file(s) holding old L2GP data
    !    OL2FileHandle = mlspcf_ol2gp_start

    ! starting quantities number for *this* vector; what if we have more?
    !    qtiesStart = 1
    !   Calculate qtiesStart for the specific quantity below

    error = 0
    templateIndex = -1
    vectorIndex = -1
    spread = .false.
    interpolate = .false.
    extinction = .false.
    maxIterations = 4
    noFineGrid = 1

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
      got= .false.
      resetSeed = .false.
      switch2intrinsic = .false.
      seed = 0

      ! Node_id(key) is now n_spec_args.

      select case( get_spec_id(key) )
      case ( s_vector ) ! ===============================  Vector  =====
        got = .false.
        globalUnit = PHYQ_Invalid
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
          end select
        end do

        if ( all(got((/f_lengthScale,f_fraction/))) ) &
          & call Announce_Error ( key, bothFractionAndLength )
        ! Create the vector, and add it to the database.
        call decorate ( key, AddVectorToDatabase ( vectors, &
          & CreateVector ( vectorName, vectorTemplates(templateIndex), &
          & qtyTemplates, globalUnit=globalUnit ) ) )

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
            call createEmptyMatrix ( matrixCholesky%m, vectorName, &
              & vectors(rowVector), vectors(colVector) )
            call decorate ( key, addToMatrixDatabase(matrices, matrixCholesky) )
          case ( l_kronecker )
            call createEmptyMatrix ( matrixKronecker%m, vectorName, &
              & vectors(rowVector), vectors(colVector) )
            call decorate ( key, addToMatrixDatabase(matrices, matrixKronecker) )
          case ( l_plain )
            call createEmptyMatrix ( matrixPlain, vectorName, vectors(rowVector), &
              vectors(colVector) )
            call decorate ( key, addToMatrixDatabase(matrices, matrixPlain) )
          case ( l_spd )
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
          case ( f_h2oQuantity ) ! For hydrostatic
            h2oVectorIndex = decoration(decoration(subtree(1,gson)))
            h2oQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_dontMask )
            dontMask = get_boolean ( gson )
          case ( f_ignoreZero )
            ignoreZero = get_boolean ( gson )
          case ( f_ignoreNegative )
            ignoreNegative = get_boolean ( gson )
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
          case ( f_losQty ) ! For losGrid fill
            losVectorIndex = decoration(decoration(subtree(1,gson)))
            losQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
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
          case ( f_noFineGrid )      ! For cloud extinction fill
            call expr ( subtree(2,subtree(j,key)), unitAsArray,valueAsArray )
            if ( all(unitAsArray(1) /= (/PHYQ_Dimensionless,PHYQ_Invalid/)) ) &
              & call Announce_error ( key, badUnitsForMaxIterations )
            noFineGrid = valueAsArray(1)
          case ( f_PtanQuantity ) ! For losGrid fill
            PtanVectorIndex = decoration(decoration(subtree(1,gson)))
            PtanQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_quantity )   ! What quantity are we filling quantity=vector.quantity
            vectorIndex = decoration(decoration(subtree(1,gson)))
            quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_radianceQuantity )      ! For estimated noise
            radianceVectorIndex = decoration(decoration(subtree(1,gson)))
            radianceQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_ratioQuantity )      ! For isotope ratio
            ratioVectorIndex = decoration(decoration(subtree(1,gson)))
            ratioQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_refGPHQuantity ) ! For hydrostatic
            refGPHVectorIndex = decoration(decoration(subtree(1,gson)))
            refGPHQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))    
          case ( f_resetSeed )
            resetSeed = get_boolean ( gson )
          case ( f_sourceQuantity )       ! When filling from a vector, what vector/quantity
            sourceVectorIndex = decoration(decoration(subtree(1,gson)))
            sourceQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_scECI )                ! For special fill of losVel
            scECIVectorIndex = decoration(decoration(subtree(1,gson)))
            scECIQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_scVel )                ! For special fill of losVel
            scVelVectorIndex = decoration(decoration(subtree(1,gson)))
            scVelQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_seed ) ! For explicit fill
            valuesNode=subtree(j,key)
          case ( f_sourceL2AUX )          ! Which L2AUXDatabase entry to use
            l2auxIndex = decoration(decoration(gson))
          case ( f_sourceL2GP )           ! Which L2GPDatabase entry to use
            l2gpIndex=decoration(decoration(gson))
          case ( f_sourceGrid )
            gridIndex=decoration(decoration(gson))
          case ( f_sourceVGrid )
            vGridIndex=decoration(decoration(gson))
          case ( f_spread ) ! For explicit fill, note that gson here is not same as others
            if ( node_id(gson) == n_set_one ) then
              spread=.TRUE.
            else
              spread = decoration(subtree(2,gson)) == l_true
            end if
          case ( f_systemTemperature )
            call expr ( gson , unitAsArray, valueAsArray )
            if ( all (unitAsArray /= (/PHYQ_Temperature, PHYQ_Invalid/) ) ) &
              call Announce_error ( key, badUnitsForSystemTemperature )
            systemTemperature = valueAsArray(1)
          case ( f_tngtECI )              ! For special fill of losVel
            tngtECIVectorIndex = decoration(decoration(subtree(1,gson)))
            tngtECIQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_temperatureQuantity ) ! For hydrostatic
            temperatureVectorIndex = decoration(decoration(subtree(1,gson)))
            temperatureQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
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
          endif
          if ( resetSeed ) then
            call mls_random_seed(new_seed=seed(1:))
            if (DEEBUG) then
              call output('Letting mls choose new seed ', advance='no')
              call output(seed, advance='yes')
            endif
          elseif ( got(f_seed) ) then
            do j=1, nsons(valuesNode)-1
              call expr(subtree(j+1,valuesNode),unitAsArray,valueAsArray)
              seed(j) = int(valueAsArray(1))
            enddo
            if ( seed(1) /= 0 .and. seed(2) /= 0 ) then
              call mls_random_seed(pput=seed(1:))
              if (DEEBUG) then
                call output('Setting new seed ', advance='no')
                call output(seed, advance='yes')
              endif
            else
              call mls_random_seed(new_seed=seed(1:))
              if (DEEBUG) then
                call output('Letting mls choose new seed ', advance='no')
                call output(seed, advance='yes')
              endif
            endif
          else
            call mls_random_seed(gget=seed(1:))
            if (DEEBUG) then
              call output('Reusing current seed ', advance='no')
              call output(seed, advance='yes')
            endif
          endif
          call addGaussianNoise ( key, quantity, sourceQuantity, &
            & noiseQty )

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
            if ( .not. got(f_h2oQuantity) ) &
              & call Announce_Error ( key, needH2O )
            h2oQuantity => GetVectorQtyByTemplateIndex( &
              & vectors(h2oVectorIndex), h2oQuantityIndex)
            if ( .not. ValidateVectorQuantity(h2oQuantity, &
              & quantityType=(/l_vmr/), molecule=(/l_h2o/)) )&
              & call Announce_Error ( key, badGeocAltitudeQuantity )
          else
            nullify ( geocAltitudeQuantity, h2oQuantity )
          end if
          call FillVectorQtyHydrostatically ( key, quantity, temperatureQuantity, &
            & refGPHQuantity, h2oQuantity, geocAltitudeQuantity, maxIterations )          

        case ( l_isotope ) ! --------------- Isotope based fills -------
          if (.not. all(got( (/f_ratioQuantity, f_sourceQuantity/) ) ) ) &
            & call Announce_Error ( key, badIsotopeFill )
          ratioQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(ratioVectorIndex), ratioQuantityIndex )
          sourceQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(sourceVectorIndex), sourceQuantityIndex )
          call FillVectorQtyFromIsotope ( key, quantity, sourceQuantity, &
            & ratioQuantity )

        case ( l_rectanglefromlos ) ! -------fill from losGrid quantity -------
          if (.not. all(got((/f_losQty,f_earthRadius,f_PtanQuantity/))))&
            & call Announce_Error ( key, badlosGridFill )
          earthRadiusQty => GetVectorQtyByTemplateIndex( &
            & vectors(earthRadiusVectorIndex), earthRadiusQtyIndex )
          TngtPresQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(PtanVectorIndex), PtanQtyIndex )
          losQty => GetVectorQtyByTemplateIndex( &
            & vectors(losVectorIndex), losQtyIndex )
          call FillQuantityFromLosGrid ( key, Quantity, losQty, &
            & tngtPresQuantity, earthRadiusQty, &
            & noFineGrid, extinction, errorCode )

        case ( l_special ) ! -  Special fills for some quantities  -----
          select case ( quantity%template%quantityType )
          case ( l_losVel )
            if ( .not. any(got( (/f_tngtECI, f_scECI, f_scVel/) )) ) then
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
            endif
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
                & dontMask, ignoreZero, ignoreNegative )
            endif
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
                & dontMask, ignoreZero, ignoreNegative )
            endif
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
                & dontMask, ignoreZero, ignoreNegative )
            endif
          case default
            call Announce_error ( key, noSpecialFill )
          end select

        case ( l_vGrid ) ! ---------------------  Fill from vGrid  -----
          if (.not. ValidateVectorQuantity(quantity, &
            & quantityType=(/l_ptan/), &
            & frequencyCoordinate=(/l_none/) ) ) &
            & call MLSMessage ( MLSMSG_Error,ModuleName,&
            &   'vGrids can only be used to fill ptan quantities' )
          if ( vGrids(vGridIndex)%verticalCoordinate /= l_zeta ) &
            & call MLSMessage ( MLSMSG_Error,ModuleName, &
            &  'Vertical coordinate in vGrid is not zeta' )
          if ( vGrids(vGridIndex)%noSurfs /= quantity%template%noSurfs )&
            & call MLSMessage ( MLSMSG_Error,ModuleName, &
            &  'VGrid is not of the same size as the quantity' )
          do instance = 1, quantity%template%noInstances
            quantity%values(:,instance) = vGrids(vGridIndex)%surfs
          end do
          !quantity%values = spread ( vGrids(vGridIndex)%surfs, 2, &
          !  & quantity%template%noInstances )

        case ( l_gridded ) ! ------------  Fill from gridded data  -----
          if ( .not. got(f_sourceGrid) ) &
            & call Announce_Error ( key,noSourceGridGiven )
          call FillVectorQuantityFromGrid &
            & ( quantity, griddedData(gridIndex), errorCode )
          if ( errorCode /= 0 ) call Announce_error ( key, errorCode )

        case ( l_l2gp ) ! --------------  Fill from L2GP quantity  -----
          if ( .NOT. got(f_sourceL2GP) ) &
            & call Announce_Error ( key, noSourceL2GPGiven )
          call FillVectorQuantityFromL2GP &
            & ( quantity, l2gpDatabase(l2gpIndex), interpolate, errorCode )
          if ( errorCode /= 0 ) call Announce_error ( key, errorCode )

        case ( l_l2aux ) ! ------------  Fill from L2AUX quantity  -----
          if ( .NOT. got(f_sourceL2AUX) ) &
            & call Announce_Error ( key, noSourceL2AUXGiven )
!          call FillVectorQuantityFromL2AUX(quantity,l2auxDatabase(l2auxIndex),errorCode)
          if ( errorCode /= 0 ) call Announce_error ( key, errorCode )

        case ( l_estimatedNoise ) ! ----------- Fill with estimated noise ---
          radianceQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(radianceVectorIndex), radianceQuantityIndex )
          call FillVectorQtyWithEstNoise ( &
            & quantity, radianceQuantity, systemTemperature, integrationTime )

        case ( l_explicit ) ! ---------  Explicity fill from l2cf  -----
          if ( .not. got(f_explicitValues) ) &
            & call Announce_Error ( key, noExplicitValuesGiven )
          call ExplicitFillVectorQuantity ( quantity, valuesNode, spread, &
            & vectors(vectorIndex)%globalUnit )

        case ( l_l1b ) ! --------------------  Fill from L1B data  -----
          call FillVectorQuantityFromL1B ( key, quantity, chunks(chunkNo), l1bInfo )

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

      ! End of fill operations

      case ( s_remove ) ! ===============================  Remove ==
        if (DEEBUG) call output('Remove vector instruction', advance='yes')
        ! Here we're to try to shrink the vector database by removing a vector
        ! Loop over the instructions
        ! (Shall we allow multiple rms on a single line? Maybe later)
        do j = 2, nsons(key)
          gson = subtree(j,key)  ! The argument
          fieldIndex = get_field_id(gson)
          gson = subtree(2,gson) ! Now the value of said argument
          select case ( fieldIndex )
          case ( f_source )
            sourceVectorIndex = decoration(gson)
          case default ! Can't get here if type checker worked
          end select
        end do

        vectorindex = rmVectorFromDatabase ( vectors, vectors(sourceVectorIndex) )

      case ( s_transfer ) ! ===============================  Transfer ==
        ! Here we're on a transfer instruction
        ! Loop over the instructions
        do j = 2, nsons(key)
          gson = subtree(j,key)  ! The argument
          fieldIndex = get_field_id(gson)
          gson = subtree(2,gson) ! Now the value of said argument
          select case ( fieldIndex )
          case ( f_source )
            sourceVectorIndex = decoration(decoration(gson))
          case ( f_destination )
            destinationVectorIndex = decoration(decoration(gson))
          case default ! Can't get here if type checker worked
          end select
        end do

        call TransferVectors ( vectors(sourceVectorIndex), &
          & vectors(destinationVectorIndex) )

      case ( s_time ) ! ===================================  Time  =====
        if ( timing ) then
          call sayTime
        else
          call cpu_time ( t1 )
          timing = .true.
        end if

      case ( s_snoop )
        call Snoop ( key=key, vectorDatabase=vectors )

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
    if ( timing ) call sayTime

  contains
    subroutine SayTime
      call cpu_time ( t2 )
      if ( total_times ) then
        call output ( "Total time = " )
        call output ( dble(t2), advance = 'no' )
        call blanks ( 4, advance = 'no' )
      endif
      call output ( "Timing for MLSL2Fill = " )
      call output ( dble(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime
  end subroutine MLSL2Fill

  ! =====     Private Procedures     =====================================


  ! ------------------------------------------- addGaussianNoise ---
  subroutine addGaussianNoise ( key, quantity, sourceQuantity, &
            & noiseQty )
    ! A special fill: quantity = sourceQuantity + g() noiseQty
    ! where g() is a random number generator with mean 0 and std. dev. 1
    ! Formal arguments
    integer, intent(in) :: KEY
    type (VectorValue_T), intent(inout) ::   quantity
    type (VectorValue_T), intent(in) ::      sourceQuantity
    type (VectorValue_T), intent(in) ::      noiseQty

    ! Local variables
    integer                          ::    ROW, COLUMN

    ! Executable code
    ! First check that things are OK.
    if (.not. FillableChiSq ( quantity, &
      & sourceQuantity, noiseQty ) ) then
      call Announce_error ( key, No_Error_code, &
      & 'Incompatibility among vector quantities adding noise'  )
      return
    endif
    
    do column=1, size(quantity%values(1, :))
      do row=1, size(quantity%values(:, 1))
        quantity%values(row, column) = sourceQuantity%values(row, column) &
          & + &
          & drang() * noiseQty%values(row, column)
      enddo
    enddo

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
    real (r8), dimension(:,:), pointer :: M ! The matrix being filled
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
      
      call Dump ( LMasked )

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
          m = 0.0
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
                  distance = abs ( surfs ( j/qt%noChans ) - surfs ( k/qt%noChans ) )
                case ( l_zeta )
                  distance = abs ( surfs ( j/qt%noChans ) - surfs ( k/qt%noChans ) ) * decade
                case ( l_pressure )
                  distance = abs ( -log10 ( surfs(j/qt%noChans) ) + &
                    &               log10 ( surfs(k/qt%noChans) ) ) / decade
                end select
                if ( meanLength > 0.0 ) &
                  & m(j,k) = meanDiag*thisFraction*exp(-distance/meanLength)
              end do                    ! Loop over k (in M)
            end do                      ! Loop over j (in M)
          end if                        ! An appropriate vertical coordinate

          ! Now we may need to invert this, if so we need to be clever.
          if ( invert ) then
            call Allocate_test ( condition, n, 'condition', ModuleName )
            condition = d%values(:,i) <= 0.0
            do j = 1, n
              if ( condition(j) ) M(j,j) = 1.0
            end do
            call MatrixInversion(M)
            do j = 1, n
              if ( condition(j) ) M(j,j) = 0.0
            end do
            call Deallocate_test ( condition, 'condition', ModuleName )
          endif

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
    endif

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
  subroutine FillVectorQuantityFromL2GP ( quantity,l2gp, interpolate, errorCode )

    ! If the times, pressures, and geolocations match, fill the quantity with
    ! the appropriate subset of profiles from the l2gp

    ! Dummy arguments
    type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
    type (L2GPData_T), intent(in) :: L2GP ! L2GP to fill from
    logical, intent(in) :: interpolate  ! Flag
    integer, intent(out) :: errorCode ! Error code

    ! Local parameters
    real(r8), parameter:: TOLERANCE=0.05 ! Tolerence for angles

    ! Local variables
    integer ::    FIRSTPROFILE, LASTPROFILE
    integer, dimension(1) :: FIRSTPROFILEASARRAY
    integer :: INSTANCE                 ! Loop counter 

    real (r8), dimension(quantity%template%noSurfs) :: outZeta

    errorCode=0
    ! Make sure this quantity is appropriate
    if (.not. ValidateVectorQuantity(quantity, coherent=.TRUE., stacked=.TRUE., &
      & verticalCoordinate= (/ l_pressure, l_zeta /) ) ) then
      errorCode=vectorWontMatchL2GP
      return
    end if

    if ( (quantity%template%noChans/=l2gp%nFreqs) .and. &
      &  ((quantity%template%noChans/=1) .or. (l2gp%nFreqs/=0)) ) then
      errorCode=vectorWontMatchL2GP
      return
    end if

    if ( quantity%template%noSurfs /= l2gp%nLevels .and. (.not. interpolate) ) then
      errorCode=vectorWontMatchL2GP
      return
    end if

    if (.not. interpolate) then 
      if ( quantity%template%verticalCoordinate == l_pressure ) then
        if ( any(ABS(-LOG10(quantity%template%surfs(:,1))+ &
          & LOG10(l2gp%pressures)) > TOLERANCE) ) then
          errorCode=vectorWontMatchL2GP
          return
        end if
      else                                ! Must be l_zeta
        if ( any(ABS(quantity%template%surfs(:,1)+ &
          & LOG10(l2gp%pressures)) > TOLERANCE) ) then
          errorCode=vectorWontMatchL2GP
          return
        end if
      end if
    end if

    ! Attempt to match up the first location
    firstProfileAsArray=MINLOC(ABS(quantity%template%phi(1,1)-l2gp%geodAngle))
    firstProfile=firstProfileAsArray(1)

    ! Well, the last profile has to be noInstances later, check this would be OK
    lastProfile=firstProfile+quantity%template%noInstances-1
    if (lastProfile > l2gp%nTimes ) then
      errorCode=vectorWontMatchL2GP
      return
    end if

    ! Now check that geodAngle's are a sufficient match
    if (any(abs(l2gp%geodAngle(firstProfile:lastProfile)-&
      &         quantity%template%phi(1,:)) > tolerance) ) then
      errorCode=vectorWontMatchL2GP
      return
    end if

    if (any(abs(l2gp%time(firstProfile:lastProfile)- &
      &         quantity%template%time(1,:)) > tolerance) ) then
      errorCode=vectorWontMatchL2GP
      return
    end if

    if (interpolate .and. quantity%template%noChans /= 1) then
      errorCode=cantInterpolate3D
      return
    endif

    if (.not. interpolate) then
      quantity%values=RESHAPE(l2gp%l2gpValue(:,:,firstProfile:lastProfile),&
        & (/quantity%template%noChans*quantity%template%noSurfs,&
        &   quantity%template%noInstances/))
    else
      if ( quantity%template%verticalCoordinate == l_pressure ) then
        outZeta = -log10 ( quantity%template%surfs(:,1) )
      else
        outZeta = quantity%template%surfs(:,1)
      endif
      do instance = 1, quantity%template%noInstances
        call InterpolateValues ( &
          & -log10(l2gp%pressures), &  ! Old X
          & l2gp%l2gpValue(1,:,firstProfile+instance-1), & ! OldY
          & outZeta, & ! New X
          & quantity%values(:,instance), & ! New Y
          & method='Linear', extrapolate='Clamp' )
      enddo
    endif

  end subroutine FillVectorQuantityFromL2GP

  ! ------------------------------------------- FillLOSVelocity ---
  subroutine FillLOSVelocity ( key, qty, tngtECI, scECI, scVel)
    ! A special fill from geometry arguments
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
    if ( (qty%template%quantityType /= l_losVel) .or. &
      &  (tngtECI%template%quantityType /= l_tngtECI) .or. &
      &  (scECI%template%quantityType /= l_scECI) .or. &
      &  (scVel%template%quantityType /= l_scVel) ) then
      call Announce_Error ( key, badLOSVelFill )
      return
    end if

    if ( qty%template%instrumentModule /= tngtECI%template%instrumentModule ) then
      call Announce_Error ( key, badLOSVelFill )
      return
    end if

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
    endif
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
    endif

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
    endif

    ! (2)
    if (minorFrame ) then
      aok = aok .and. &
        & (qty%template%signal == measQty%template%signal) &
        & .and. &
        & (qty%template%signal == modelQty%template%signal)
      if ( present(noiseQty) ) aok = aok &
        & .and. &
        & (qty%template%signal == noiseQty%template%signal)
    endif

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
      endif
    endif
    
    return
  end function FillableChiSq

  ! ------------------------------------------- FillChiSqChan ---
  subroutine FillChiSqChan ( key, qty, measQty, modelQty, noiseQty, &
  & dontMask, ignoreZero, ignoreNegative, firstInstance, lastInstance )
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

    ! Executable code
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
      endif
      return
    elseif (.not. FillableChiSq ( qty, measQty, modelQty, noiseQty ) ) then
      call Announce_error ( key, No_Error_code, &
      & 'Incompatibility among vector quantities filling chi^2 channelwise'  )
      return
    elseif (any ( noiseQty%values == 0.0) .and. &
      & .not. (ignoreZero .or. .not. dontMask) ) then
      call Announce_error ( key, No_Error_code, &
      & 'A vanishing error filling chi^2 channelwise'  )
      return
    endif

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
            & (measQty%values(qIndex, i) - modelQty%values(qIndex, i)) &
            & / &
            & noiseQty%values(qIndex, i) &
            &  ) ** 2
            N = N + 1
          endif
        enddo
        if ( N > 0 ) then
          qty%values(c, i) = sum(values) / N
        else
          qty%values(c, i) = 0.
        endif
      enddo
    enddo
    call deallocate_test(values, &
      & 'chi^2 unsummed', ModuleName)
  end subroutine FillChiSqChan

  ! ------------------------------------------- FillChiSqMMaf ---
  subroutine FillChiSqMMaf ( key, qty, measQty, modelQty, noiseQty, &
  & dontMask, ignoreZero, ignoreNegative, firstInstance, lastInstance )
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
      endif
      return
    elseif (.not. FillableChiSq ( qty, measQty, modelQty, noiseQty ) ) then
      call Announce_error ( key, No_Error_code, &
      & 'Incompatibility among vector quantities filling chi^2 MMAFwise'  )
      return
    elseif (any ( noiseQty%values == 0.0) .and. &
      & .not. (ignoreZero .or. .not. dontMask) ) then
      call Announce_error ( key, No_Error_code, &
      & 'A vanishing noise filling chi^2 MMAFwise'  )
      return
    endif

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
            & (measQty%values(row, i) - modelQty%values(row, i)) &
            & / &
            & noiseQty%values(row, i) &
            &  ) ** 2
            N = N + 1
          endif
        enddo
        if ( N > 0 ) then
          qty%values(1, i) = sum(values) / N
        else
          qty%values(1, i) = 0.
        endif
      endif
    enddo
    call deallocate_test(values, &
      & 'chi^2 unsummed', ModuleName)
  end subroutine FillChiSqMMaf

  ! ------------------------------------------- FillChiSqMMif ---
  subroutine FillChiSqMMif ( key, qty, measQty, modelQty, noiseQty, &
  & dontMask, ignoreZero, ignoreNegative, firstInstance, lastInstance )
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
      endif
      return
    elseif (.not. FillableChiSq ( qty, measQty, modelQty, noiseQty ) ) then
      call Announce_error ( key, No_Error_code, &
      & 'Incompatibility among vector quantities filling chi^2 MMIFwise'  )
      return
    elseif (any ( noiseQty%values == 0.0) .and. &
      & .not. (ignoreZero .or. .not. dontMask) ) then
      call Announce_error ( key, No_Error_code, &
      & 'A vanishing noise filling chi^2 MMIFwise'  )
      return
    endif

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
            & (measQty%values(qIndex, i) - modelQty%values(qIndex, i)) &
            & / &
            & noiseQty%values(qIndex, i) &
            &  ) ** 2
            N = N + 1
          endif
        enddo
        if ( N > 0 ) then
          qty%values(s, i) = sum(values) / N
        else
          qty%values(s, i) = 0.
        endif
      enddo
    enddo
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
    elseif ( qty%template%molecule /= vmrQty%template%molecule) then
          call Announce_error ( key, No_Error_code, &
              & 'Attempt to fill column abundance with different molecule'  )
      return
    elseif ( &
    & .not. ( &
    & DoHgridsMatch( qty, vmrQty ) &
    & .and. &
    & DoHgridsMatch( qty, bndPressQty ) &
    & ) &
    & ) then
          call Announce_error ( key, No_Error_code, &
              & 'Attempt to fill column abundance with different HGrids'  )
      return
    elseif ( .not. &
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
   endif   
!   do instance=1, vmrQty%template%noInstances
   do instance=useFirstInstance, useLastInstance
      if(vmrQty%template%coherent) then
         surfaceInstance=1
      else
         surfaceInstance=instance
      endif

      if(zeta_surfs) then
         ! Invert zeta = -log10(p)
         do surface=1, vmrQty%template%noSurfs
            p(surface) = exp(-log(10.)* &
            & vmrQty%template%surfs(surface, surfaceInstance))
         enddo
      else
         do surface=1, vmrQty%template%noSurfs
            p(surface) = vmrQty%template%surfs(surface, surfaceInstance)
         enddo
      endif
 
      if(p(1) &
         &  < bndPressQty%values(1, instance)) then
          call Announce_error ( key, No_Error_code, &
              & 'Fill column abundance, but tropopause below VGrid'  )
      endif

   ! Find 1st surface at or above tropopause
   ! (i.e., at a pressure equal to or less than boundaryPressure)
      do surface=1, vmrQty%template%noSurfs
         if(p(surface) &
         &  <= bndPressQty%values(1, instance)) exit
      enddo
      firstSurface = surface
      if(firstSurface > vmrQty%template%noSurfs-2) then
          call Announce_error ( key, No_Error_code, &
              & 'Fill column abundance, but tropopause above VGrid'  )
      endif
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
      enddo
      qty%values(1, instance) = AoverMg * columnSum
   enddo

   deallocate(p, stat=status)
   if(status /= 0) then
          call Announce_error ( key, No_Error_code, &
              & 'Error in deallocating p'  )
   endif   

  end subroutine FillColAbundance

  ! ---------------------------------- FillVectorQuantityWithEsimatedNoise ---
  subroutine FillVectorQtyWithEstNoise ( quantity, radiance, &
    & systemTemperature, integrationTime )

    use MLSSignals_m, only: signals

    ! Dummy arguments
    type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
    type (VectorValue_T), intent(in) :: RADIANCE ! Radiances to use in calculation
    real(r8), intent(in) :: SYSTEMTEMPERATURE ! System temperature in K
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

    if ( radiance%template%signal /= quantity%template%signal ) &
      & call MLSMEssage ( MLSMSG_Error, ModuleName, &
      & "Quantity and radiances not same signal for estimated noise fill.")

    if ( DEEBUG .and. .not. associated(quantity%values) ) &
      & call announce_error( 0, No_Error_code, &
      & 'quantity values unassociated in FillVectorQtyWithEstNoise')

    width => signals(radiance%template%signal)%widths

    i = 1
    do s = 1, quantity%template%noSurfs
      do c = 1, quantity%template%noChans
        quantity%values(i,:) = ( radiance%values(i,:) + systemTemperature ) / &
          & sqrt ( integrationTime * 10e6 * width(c) )
        i = i + 1
      end do
    end do

  end subroutine FillVectorQtyWithEstNoise

  ! ------------------------------------- FillVectorHydrostatically ----
  subroutine FillVectorQtyHydrostatically ( key, quantity, &
    & temperatureQuantity, refGPHQuantity, h2oQuantity, &
    & geocAltitudeQuantity, maxIterations )
    ! Various hydrostatic fill operations
    integer, intent(in) :: key          ! For messages
    type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
    type (VectorValue_T), intent(in) :: TEMPERATUREQUANTITY
    type (VectorValue_T), intent(in) :: REFGPHQUANTITY
    type (VectorValue_T), pointer :: H2OQUANTITY
    type (VectorValue_T), pointer :: GEOCALTITUDEQUANTITY
    integer, intent(in) :: MAXITERATIONS
    ! H2OQuantity and GeocAltitudeQuantity have to be pointers
    ! as they may be absent.

    ! Local variables

    ! Executable code

    if ( toggle(gen) ) call trace_begin ( "FillVectorQtyHydrostatically", key )

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
	if ( toggle(gen) ) call trace_end ( "FillVectorQtyHydrostatically")
        return
      end if
      if ( (any(quantity%template%surfs /= temperatureQuantity%template%surfs)) .or. &
        & (any(quantity%template%phi /= temperatureQuantity%template%phi)) .or. &
        & (any(quantity%template%phi /= refGPHQuantity%template%phi)) ) then
        call Announce_Error ( key, nonConformingHydrostatic, &
          &  "case l_gph failed second test" )
	if ( toggle(gen) ) call trace_end ( "FillVectorQtyHydrostatically")
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
	if ( toggle(gen) ) call trace_end ( "FillVectorQtyHydrostatically")
        return
      end if
      if ( (any(refGPHquantity%template%phi /= temperatureQuantity%template%phi)) .or. &
        & (any(h2oQuantity%template%phi /= temperatureQuantity%template%phi)) ) then
        call Announce_Error ( key, nonConformingHydrostatic, &
          & "case l_ptan failed second test" )
	if ( toggle(gen) ) call trace_end ( "FillVectorQtyHydrostatically")
        return
      end if
      if ( (.not. ValidateVectorQuantity(quantity, minorFrame=.true.) ) .or. &
        &  (.not. ValidateVectorQuantity(geocAltitudeQuantity, minorFrame=.true.) ) .or. &
        &  (quantity%template%instrumentModule /= &
        &   geocAltitudeQuantity%template%instrumentModule) )  then
        call Announce_Error ( key, nonConformingHydrostatic, &
          & "case l_ptan failed third test" )
        print *, 'ValidateVectorQuantity(quantity, minorFrame=.true.) ', &
        &  ValidateVectorQuantity(quantity, minorFrame=.true., sayWhyNot=.true.)
        print *, 'ValidateVectorQuantity(geocAltitudeQuantity, minorFrame=.true.) ', &
        & ValidateVectorQuantity(geocAltitudeQuantity, minorFrame=.true., sayWhyNot=.true.)
        print *, 'quantity%template%instrumentModule ', &
        & quantity%template%instrumentModule
        print *, 'geocAltitudeQuantity%template%instrumentModule ', &
        & geocAltitudeQuantity%template%instrumentModule
	if ( toggle(gen) ) call trace_end ( "FillVectorQtyHydrostatically")
        return
      end if
      call GetHydrostaticTangentPressure ( quantity, temperatureQuantity,&
        & refGPHQuantity, h2oQuantity, geocAltitudeQuantity, maxIterations )
    case default
      call Announce_error ( 0, 0, 'No such fill yet' )
    end select

    if ( toggle(gen) ) call trace_end ( "FillVectorQtyHydrostatically" )

  end subroutine FillVectorQtyHydrostatically

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
  subroutine ExplicitFillVectorQuantity(quantity, valuesNode, spread, globalUnit)

    ! This routine is called from MLSL2Fill to fill values from an explicit
    ! fill command line

    ! Dummy arguments
    type (VectorValue_T), intent(inout) :: QUANTITY ! The quantity to fill
    integer, intent(in) :: VALUESNODE ! Tree node
    logical, intent(in) :: SPREAD ! One instance given, spread to all
    integer, intent(in) :: GLOBALUNIT   ! From parent vector

    ! Local variables
    integer :: K                ! Loop counter
    integer, DIMENSION(2) :: unitAsArray ! Unit for value given
    real (r8), DIMENSION(2) :: valueAsArray ! Value given
    integer :: TestUnit                 ! Unit to use
    
    ! Executable code
    testUnit = quantity%template%unit
    if ( globalUnit /= phyq_Invalid ) testUnit = globalUnit

    if (spread) then      ! 1 instance/value given, spread to all instances

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
        if (nsons(valuesNode)-1 == 1) then
          quantity%values=valueAsArray(1)
        else
          quantity%values(k,:)=valueAsArray(1)
        endif
      end do

    else                  ! Not spread, fill all values

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
        quantity%values(mod(k-1,quantity%template%instanceLen)+1,&
          &             (k-1)/quantity%template%instanceLen+1)=&
          & valueAsArray(1)
      end do
    end if
  end subroutine ExplicitFillVectorQuantity

  ! ----------------------------------------- FillVectorQuantityFromL1B ----
  subroutine FillVectorQuantityFromL1B ( root, quantity, chunk, l1bInfo )
    integer, intent(in) :: root
    type (VectorValue_T), INTENT(INOUT) :: QUANTITY
    type (MLSChunk_T), INTENT(IN) :: CHUNK
    type (l1bInfo_T), INTENT(IN) :: L1BINFO

    ! Local variables
    character (len=80) :: NAMESTRING
    integer :: fileID, FLAG, NOMAFS
    type (l1bData_T) :: L1BDATA

    ! Executable code

    if ( toggle(gen) ) call trace_begin ("FillVectorQuantityFromL1B",root)

    fileID=l1bInfo%l1bOAID
    select case ( quantity%template%quantityType )
    case ( l_radiance )
      call GetSignalName ( quantity%template%signal, nameString, noChannels=.TRUE. )
      fileID = FindL1BData (l1bInfo%l1bRadIDs, nameString )
    case ( l_tngtECI )
      call GetModuleName( quantity%template%instrumentModule,nameString )
      nameString=TRIM(nameString)//'.tpECI'
    case ( l_tngtGeodAlt )
      call GetModuleName( quantity%template%instrumentModule,nameString )
      nameString=TRIM(nameString)//'.tpGeodAlt'
    case ( l_tngtGeocAlt )
      call GetModuleName( quantity%template%instrumentModule,nameString )
      nameString=TRIM(nameString)//'.tpGeocAlt'
    case ( l_scECI )
      nameString='scECI'
    case ( l_scVel )
      nameString='scVel'
    case ( l_scGeocAlt )
      nameString='scGeocAlt'
    case default
      call Announce_Error ( root, cantFillFromL1B )
    end select

    call ReadL1BData ( fileID , nameString, l1bData, noMAFs, flag, &
      & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex )
    ! We'll have to think about `bad' values here .....
    if ( flag /= 0 ) then
      call Announce_Error ( root, errorReadingL1B )
      if ( toggle(gen) ) call trace_end ( "FillVectorQuantityFromL1B")
      return
    end if
    quantity%values = RESHAPE(l1bData%dpField, &
      & (/ quantity%template%instanceLen, quantity%template%noInstances /) )
    call DeallocateL1BData(l1bData)

    if (toggle(gen) ) call trace_end( "FillVectorQuantityFromL1B" )
  end subroutine FillVectorQuantityFromL1B

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

  ! ---------------------------------------------- TRANSFERVECTORS -----
  subroutine TransferVectors ( source, dest )
    ! Copy common items in source to those in dest
    type (Vector_T), intent(in) :: SOURCE
    type (Vector_T), intent(inout) :: DEST

    ! Local variables
    integer :: Dummy                   ! Dummy integer
    type (VectorValue_T), pointer :: DQ ! Destination quantity
    type (VectorValue_T), pointer :: SQ ! Source quantity
    integer :: SQI                      ! Quantity index in source

    ! Executable code

    ! First copy those things in source, loop over them
    dest%globalUnit = source%globalUnit
    do sqi = 1, size ( source%quantities )
      ! Try to find this in dest
      sq => source%quantities(sqi)
      dq => GetVectorQtyByTemplateIndex ( dest, source%template%quantities(sqi), dummy )
      if ( associated ( dq ) ) then
        dq%values = sq%values
        if (associated(sq%mask)) then
          if (.not. associated(dq%mask)) call CreateMask ( dq )
          dq%mask = sq%mask
        else
          if ( associated(dq%mask) ) &
            & call Deallocate_test ( dq%mask, 'dq%mask', ModuleName )
        endif
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

    call output ( ': ' )
    call output ( "The " );

    select case ( code )
    case ( allocation_err )
      call output ( " command caused an allocation error in squeeze.", advance='yes' )
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
    case ( needH2O )
      call output ( " needs H2OQuantity.", advance='yes' )
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
    case ( vectorWontMatchL2GP )
      call output ( " command found no match of vetor and L2GP (set interpolate?).",&
        & advance='yes' )
    case ( wrong_number )
      call output ( " command does not have exactly one field.", advance='yes' )
    case ( zeroGeodSpan )
      call output ( " command found zero geod. ang. span.", advance='yes' )
    case ( zeroProfilesFound )
      call output ( " command found zero profiles.", advance='yes' )
    case default
      call output ( " command caused an unrecognized programming error", advance='yes' )
    end select
    if ( present(ExtraMessage) ) then
      call output(ExtraMessage, advance='yes')
    end if
  end subroutine ANNOUNCE_ERROR


end module Fill
!=============================================================================

!
! $Log$
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
