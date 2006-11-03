! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module Fill                     ! Create vectors and fill them.
  !=============================================================================

  use MLSCommon, only: MLSFile_T, R4, R8, RM, RV, &
    & DEFAULTUNDEFINEDVALUE

  ! This module performs the Fill operation in the Level 2 software.
  ! This takes a vector template, and creates and fills an appropriate vector

  implicit none
  private
  public :: MLSL2Fill

! === (start of toc) ===
! MLSL2Fill          given a vector template, and creates and fills a vector
! === (end of toc) ===

! === (start of api) ===
! MLSL2Fill (int root, *MLSFile_T fileDataBase(:), 
!        *griddedData_T GriddedDataBase(:),
!        *vectorTemplate_T VectorTemplates(:),
!        *vector_t Vectors(:), *quantityTemplate_T QtyTemplates(:),
!        *matrix_database_T Matrices(:),
!        *l2GPData_T L2GPDatabase(:), *l2AUXData_T L2AUXDatabase(:),
!        *ForwardModelConfig_T FWModelConfig(:),
!        *mlSChunk_T Chunks(:), int ChunkNo )
! === (end of api) ===
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

  logical, parameter :: DONTPAD = .false.
  logical, parameter :: USEREICHLER = .true.

  type arrayTemp_T
     real(rv), dimension(:,:), pointer :: VALUES => NULL() ! shaped like a
  end type arrayTemp_T
  
  type(arrayTemp_T), dimension(:), save, pointer :: primitives => null()

contains ! =====     Public Procedures     =============================

  !---------------------------------------------------  MLSL2Fill  -----

  subroutine MLSL2Fill ( Root, filedatabase, GriddedDataBase, VectorTemplates, &
    & Vectors, QtyTemplates , Matrices, L2GPDatabase, L2AUXDatabase, &
    & FWModelConfig, Chunks, ChunkNo )

    ! This is the main routine for the module.  It parses the relevant lines
    ! of the l2cf and works out what to do.

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test, Test_Allocate
    use Chunks_m, only: MLSChunk_T
    use DestroyCommand_m, only: DestroyCommand
    use DumpCommand_m, only: DumpCommand
    use Expr_M, only: EXPR, EXPR_CHECK, GetIndexFlagsFromList
    use ForwardModelConfig, only: ForwardModelConfig_T
    use ForwardModelSupport, only: FillFwdModelTimings
    use GLOBAL_SETTINGS, only: BrightObjects
    use GriddedData, only: GriddedData_T, WrapGriddedData
    ! We need many things from Init_Tables_Module.  First the fields:
    use INIT_TABLES_MODULE, only: F_A, F_ADDITIONAL, F_ALLOWMISSING, &
      & F_APRIORIPRECISION, F_AVOIDBRIGHTOBJECTS, &
      & F_B, F_BADRANGE, F_BASELINEQUANTITY, F_BIN, F_BOUNDARYPRESSURE, &
      & F_BOXCARMETHOD, &
      & F_C, F_CENTERVERTICALLY, F_CHANNEL, F_COLUMNS, &
      & F_DESTINATION, F_DIAGONAL, F_DONTMASK,&
      & F_ECRTOFOV, F_EARTHRADIUS, F_EXCLUDEBELOWBOTTOM, F_EXPLICITVALUES, &
      & F_EXTINCTION, F_FIELDECR, F_FILE, F_FLAGS, F_FORCE, &
      & F_FRACTION, F_FROMPRECISION, F_GEOCALTITUDEQUANTITY, F_GPHQUANTITY, &
      & F_HEIGHT, F_HIGHBOUND, F_H2OQUANTITY, F_H2OPRECISIONQUANTITY, &
      & F_IFMISSINGGMAO, F_IGNORENEGATIVE, F_IGNOREGEOLOCATION, F_IGNOREZERO, &
      & F_INSTANCES, F_INTEGRATIONTIME, F_INTERNALVGRID, &
      & F_INTERPOLATE, F_INVERT, F_INTRINSIC, F_ISPRECISION, &
      & F_LENGTHSCALE, F_LOGSPACE, F_LOSQTY, F_LOWBOUND, F_LSB, F_LSBFRACTION, &
      & F_MANIPULATION, F_MATRIX, F_MAXITERATIONS, F_MAXVALUE, F_MEASUREMENTS, &
      & F_METHOD, F_MINNORMQTY, F_MINVALUE, F_MODEL, F_MULTIPLIER, &
      & F_NOFINEGRID, F_NOISE, F_NOISEBANDWIDTH, F_NORMQTY, &
      & F_OFFSETAMOUNT, F_ORBITINCLINATION, F_PHITAN, &
      & F_PHIWINDOW, F_PHIZERO, F_PRECISION, F_PRECISIONFACTOR, &
      & F_PROFILE, F_PROFILEVALUES, F_PTANQUANTITY, &
      & F_QUADRATURE, F_QUANTITY, F_RADIANCEQUANTITY, F_RATIOQUANTITY, &
      & F_REFRACT, F_REFGPHQUANTITY, F_REFGPHPRECISIONQUANTITY, F_RESETSEED, &
      & F_RHIQUANTITY, F_ROWS, F_SCALE, &
      & F_SCALEINSTS, F_SCALERATIO, F_SCALESURFS, &
      & F_SCECI, F_SCVEL, F_SCVELECI, F_SCVELECR, F_SEED, F_SKIPMASK, &
      & F_SOURCE, F_SOURCEGRID, F_SOURCEL2AUX, F_SOURCEL2GP, &
      & F_SOURCEQUANTITY, F_SOURCEVGRID, F_SPREAD, F_STATUS, F_SUFFIX, F_SUPERDIAGONAL, &
      & F_SYSTEMTEMPERATURE, F_TEMPERATUREQUANTITY, F_TEMPPRECISIONQUANTITY, &
      & F_TEMPLATE, F_TNGTECI, F_TERMS, &
      & F_TYPE, F_UNIT, F_USB, F_USBFRACTION, F_VECTOR, F_VMRQUANTITY, &
      & F_WHEREFILL, F_WHERENOTFILL, F_WIDTH, &
      & FIELD_FIRST, FIELD_LAST
    ! Now the literals:
    use INIT_TABLES_MODULE, only: L_ADDNOISE, L_APPLYBASELINE, L_ASCIIFILE, &
      & L_BINMAX, L_BINMEAN, L_BINMIN, L_BINTOTAL, &
      & L_BOUNDARYPRESSURE, L_BOXCAR, L_CHISQBINNED, L_CHISQCHAN, &
      & L_CHISQMMAF, L_CHISQMMIF, L_CHISQRATIO, L_CHOLESKY, &
      & L_cloudice, L_cloudextinction, L_cloudInducedRADIANCE, &
      & L_COMBINECHANNELS, L_COLUMNABUNDANCE, L_CONVERGENCERATIO, &
      & l_dnwt_flag, l_dnwt_chiSqMinNorm, l_dnwt_chiSqNorm, l_dnwt_chiSqRatio, &
      & L_DOBSONUNITS, L_DU, &
      & L_ECRTOFOV, L_ESTIMATEDNOISE, L_EXPLICIT, L_EXTRACTCHANNEL, L_FOLD, &
      & L_FWDMODELTIMING, L_FWDMODELMEAN, L_FWDMODELSTDDEV, &
      & L_GEOCALTITUDE, L_GEODALTITUDE, L_GPH, L_GPHPRECISION, L_GRIDDED, L_H2OFROMRHI, &
      & L_HEIGHT, L_HYDROSTATIC, L_ISOTOPE, L_ISOTOPERATIO, &
      & L_IWCFROMEXTINCTION, L_KRONECKER, &
      & L_L1B, L_L1BMAFBASELINE, L_L1BMIF_TAI, L_L2GP, L_L2AUX, &
      & L_LIMBSIDEBANDFRACTION, L_LOSVEL, &
      & L_LSGLOBAL, L_LSLOCAL, L_LSWEIGHTED, &
      & L_MAGAZEL, L_MAGNETICFIELD, L_MAGNETICMODEL, &
      & L_MANIPULATE, L_MAX, L_MEAN, L_MIN, L_MOLCM2, L_NEGATIVEPRECISION, &
      & L_NOISEBANDWIDTH, L_NONE, &
      & L_NORADSPERMIF, L_OFFSETRADIANCE, L_ORBITINCLINATION, &
      & L_PHASETIMING, L_PHITAN, &
      & L_PLAIN, L_PRESSURE, L_PROFILE, L_PTAN,  L_QUALITY, &
      & L_RADIANCE, L_RECTANGLEFROMLOS, L_REFGPH, L_REFRACT, &
      & L_REFLECTORTEMPMODEL, L_REFLTEMP, L_RESETUNUSEDRADIANCES, L_RHI, &
      & L_RHIFROMH2O, L_RHIPRECISIONFROMH2O, L_ROTATEFIELD, &
      & L_SCALEOVERLAPS, L_SCECI, L_SCGEOCALT, L_SCVEL, L_SCVELECI, L_SCVELECR, &
      & L_SECTIONTIMING, &
      & L_SINGLECHANNELRADIANCE, L_SPD, L_SPECIAL, L_SPREADCHANNEL, &
      & L_SPLITSIDEBAND, L_STATUS, L_SYSTEMTEMPERATURE, &
      & L_TEMPERATURE, L_TNGTECI, L_TNGTGEODALT, &
      & L_TNGTGEOCALT, L_VECTOR, L_VGRID, L_VMR, L_WMOTROPOPAUSE, &
      & L_XYZ, L_ZETA
    ! Now the specifications:
    use INIT_TABLES_MODULE, only: S_LOAD, S_DESTROY, S_DUMP, S_FILL, S_FILLCOVARIANCE, &
      & S_FILLDIAGONAL, S_FLUSHL2PCBINS, S_FLUSHPFA, S_MATRIX,  S_NEGATIVEPRECISION, &
      & S_PHASE, S_POPULATEL2PCBIN, S_SNOOP, S_TIME, &
      & S_TRANSFER, S_VECTOR, S_SUBSET, S_FLAGCLOUD, S_RESTRICTRANGE, S_UPDATEMASK
    ! Now some arrays
    use Intrinsic, only: lit_indices, &
      & PHYQ_Dimensionless, PHYQ_Invalid, PHYQ_Temperature, &
      & PHYQ_Time, PHYQ_Length, PHYQ_Pressure, PHYQ_Zeta, PHYQ_Angle, PHYQ_Profiles
    use L1BData, only: DeallocateL1BData, Dump, GetL1BFile, L1BData_T, &
      & PRECISIONSUFFIX, ReadL1BData, AssembleL1BQtyName
    use L2GPData, only: L2GPData_T, COL_SPECIES_HASH, COL_SPECIES_KEYS
    use L2AUXData, only: L2AUXData_T
    use L2PC_m, only: POPULATEL2PCBINBYNAME, LOADMATRIX, LOADVECTOR
    use LinearizedForwardModel_m, only: FLUSHLOCKEDBINS
    use L3ASCII, only: L3ASCII_INTERP_FIELD
    use ManipulateVectorQuantities, only: DOFGRIDSMATCH, DOHGRIDSMATCH, &
      & DOVGRIDSMATCH, DOQTYSDESCRIBESAMETHING, FILLWITHCOMBINEDCHANNELS
    use MatrixModule_0, only: Sparsify, MatrixInversion
    use MatrixModule_1, only: AddToMatrixDatabase, CreateEmptyMatrix, &
      & Dump, GetActualMatrixFromDatabase, GetDiagonal, &
      & FindBlock, GetKindFromMatrixDatabase, GetFromMatrixDatabase, K_Plain, K_SPD, &
      & Matrix_Cholesky_T, Matrix_Database_T, Matrix_Kronecker_T, Matrix_SPD_T, &
      & Matrix_T, NullifyMatrix, UpdateDiagonal
    ! NOTE: If you ever want to include defined assignment for matrices, please
    ! carefully check out the code around the call to snoop.
    use MLSFiles, only: GetMLSFileByType
    use MLSL2Options, only: SPECIALDUMPFILE
    use MLSL2Timings, only: SECTION_TIMES, TOTAL_TIMES, &
      & addPhaseToPhaseNames, fillTimings, finishTimings
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning, &
      & MLSMSG_Allocate, MLSMSG_Deallocate
    use MLSNumerics, only: InterpolateValues, Hunt
    use MLSRandomNumber, only: drang, mls_random_seed, MATH77_RAN_PACK
    use MLSSets, only: FindFirst, FindLast
    use MLSSignals_m, only: GetFirstChannel, GetSignalName, GetModuleName, IsModuleSpacecraft, &
      & GetSignal, Signal_T
    use MLSStringLists, only: catLists, GetHashElement, GetStringElement, &
      & NumStringElements, PutHashElement, &
      & ReplaceSubString, StringElement, StringElementNum, switchDetail
    use MLSStrings, only: lowerCase, SplitNest
    use Molecules, only: L_H2O
    use MoreTree, only: Get_Boolean, Get_Field_ID, Get_Spec_ID
    use OUTPUT_M, only: BLANKS, NEWLINE, OUTPUT, output_name_v_pair, &
      & revertoutput, switchOutput
    use PFAData_m, only: Flush_PFAData
    use QuantityTemplates, only: Epoch, QuantityTemplate_T
    use readAPriori, only: APrioriFiles
    use RHIFromH2O, only: RHIFromH2O_Factor, RHIPrecFromH2O
    use ScanModelModule, only: GetBasisGPH, Get2DHydrostaticTangentPressure, GetGPHPrecision
    use SnoopMLSL2, only: SNOOP
    use String_Table, only: Display_String, get_string
    use SubsetModule, only: SETUPSUBSET, SETUPFLAGCLOUD, RESTRICTRANGE, UPDATEMASK
    use Time_M, only: Time_Now
    use TOGGLES, only: GEN, LEVELS, SWITCHES, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, &
      & SOURCE_REF, SUB_ROSA, SUBTREE
    use TREE_TYPES, only: N_NAMED
    use UNITS, only: Deg2Rad, Rad2Deg
    use VectorsModule, only: AddVectorToDatabase, &
      & ClearUnderMask, CopyVector, CreateMask, CreateVector, &
      & DestroyVectorInfo, Dump, &
      & GetVectorQtyByTemplateIndex, isVectorQtyMasked, MaskVectorQty, &
      & ValidateVectorQuantity, Vector_T, &
      & VectorTemplate_T, VectorValue_T, M_Fill, M_LinAlg, M_Cloud
    use VGridsDatabase, only: VGRIDS, GETUNITFORVERTICALCOORDINATE

    ! Dummy arguments
    integer, intent(in) :: ROOT    ! Of the FILL section in the AST
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    type (griddedData_T), dimension(:), pointer :: GriddedDataBase
    type (vectorTemplate_T), dimension(:), pointer :: VectorTemplates
    type (vector_T), dimension(:), pointer :: Vectors
    type (quantityTemplate_T), dimension(:), pointer :: QtyTemplates
    type (matrix_database_T), dimension(:), pointer :: Matrices
    type (l2GPData_T), dimension(:), pointer :: L2GPDatabase
    type (l2AUXData_T), dimension(:), pointer :: L2AUXDatabase
    type(ForwardModelConfig_T), dimension(:), pointer :: FWModelConfig
    type (mlSChunk_T), dimension(:), pointer :: Chunks
    integer, intent(in) :: ChunkNo

    ! -----     Declarations for Fill and internal subroutines     -------

    integer :: ERROR
    logical, parameter :: DEEBUG = .FALSE.                 ! Usually FALSE

    ! Error codes for "announce_error"
    integer, parameter :: No_Error_code = 0
    integer, parameter :: CantFillFromL2AUX = No_Error_code + 1
    integer, parameter :: CantFillFromL1B = cantFillFromL2AUX + 1

    ! Error codes for "Matrix" specification
    integer, parameter :: MissingField = cantFillFromL1B + 1

    ! More Error codes relating to FillVector
    integer, parameter :: NumChansisZero = missingField + 1
    integer, parameter :: NoSourceGridGiven= numChansisZero + 1
    integer, parameter :: NoSourceL2GPGiven= noSourceGridGiven + 1
    integer, parameter :: NoSourceL2AUXGiven= noSourceL2GPGiven + 1
    integer, parameter :: NoExplicitValuesGiven= noSourceL2AUXGiven + 1
    integer, parameter :: InvalidExplicitFill = noExplicitValuesGiven + 1
    integer, parameter :: BadIsotopeFill = invalidExplicitFill + 1
    integer, parameter :: BadlosGridFill = badIsotopeFill + 1
    integer, parameter :: CantInterpolate3d = badlosGridFill + 1
    integer, parameter :: WrongUnits = CantInterpolate3d + 1

    ! Error codes resulting from FillCovariance
    integer, parameter :: NotSPD = WrongUnits + 1
    integer, parameter :: NotPlain = NotSPD + 1
    integer, parameter :: NotImplemented = notPlain + 1
    integer, parameter :: BothFractionAndLength = NotImplemented + 1

    ! Miscellaneous
    integer, parameter :: Miscellaneous_err = BothFractionAndLength + 1
    integer, parameter :: ErrorReadingL1B = miscellaneous_err + 1
    integer, parameter :: NeedTempREFGPH = errorReadingL1B + 1
    integer, parameter :: NeedH2O = needTempRefGPH + 1
    integer, parameter :: NeedOrbitInclination = needH2O + 1
    integer, parameter :: NeedGeocAltitude = needOrbitInclination + 1
    integer, parameter :: BadGeocAltitudeQuantity = needGeocAltitude + 1
    integer, parameter :: BadTemperatureQuantity = badGeocAltitudeQuantity + 1
    integer, parameter :: BadREFGPHQuantity = badTemperatureQuantity + 1
    integer, parameter :: NonConformingHydrostatic = badREFGPHQuantity + 1
    integer, parameter :: NoSpecialFill = nonConformingHydrostatic + 1
    integer, parameter :: BadlosVelFill = noSpecialFill + 1
    integer, parameter :: NotZetaForGrid = BadLosVelFill + 1
    integer, parameter :: BadEstNoiseFill = NotZetaForGrid + 1
    integer, parameter :: BadRefractFill = BadEstNoiseFill + 1
    integer, parameter :: MissingDataInGrid = BadRefractFill + 1
    integer, parameter :: EmptyGridForFill = MissingDataInGrid + 1

    ! -999.99 ! Same as %template%badvalue
    real, parameter ::    UNDEFINED_VALUE = DEFAULTUNDEFINEDVALUE

    ! Local variables

    type (vectorValue_T), pointer :: APRIORIPRECISION
    type (vectorValue_T), pointer :: AQUANTITY
    type (vectorValue_T), pointer :: BASELINEQUANTITY
    type (vectorValue_T), pointer :: BNDPRESSQTY
    type (vectorValue_T), pointer :: BQUANTITY
    type (vectorValue_T), pointer :: EARTHRADIUSQTY
    type (vectorValue_T), pointer :: ECRTOFOV
    type (vectorValue_T), pointer :: FIELDECR
    type (vectorValue_T), pointer :: FLAGQTY
    type (vectorValue_T), pointer :: GEOCALTITUDEQUANTITY
    type (vectorValue_T), pointer :: GPHQUANTITY
    type (vectorValue_T), pointer :: H2OPRECISIONQUANTITY
    type (vectorValue_T), pointer :: H2OQUANTITY
    type (vectorValue_T), pointer :: LOSQTY
    type (vectorValue_T), pointer :: LSB
    type (vectorValue_T), pointer :: LSBFRACTION
    type (vectorValue_T), pointer :: MEASQTY
    type (vectorValue_T), pointer :: MINNORMQTY
    type (vectorValue_T), pointer :: MODELQTY
    type (vectorValue_T), pointer :: NBWQUANTITY
    type (vectorValue_T), pointer :: NOISEQTY
    type (vectorValue_T), pointer :: NORMQTY
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

    logical :: ADDITIONAL               ! Flag for binMax/binMin
    logical :: ALLOWMISSING             ! Flag from l2cf
    integer :: APRPRECQTYINDEX          ! Index of apriori precision quantity
    integer :: APRPRECVCTRINDEX         ! Index of apriori precision vector
    integer :: AQTYINDEX                ! Index of a quantity in vector
    integer :: AVECINDEX                ! Index of a vector
    character(len=256) :: AVOIDOBJECTS  ! Which bright objects to avoid
    real(r8) :: BADRANGE(2)             ! Range for 'missing' data value
    integer :: BASELINEQTYINDEX
    integer :: BASELINEVCTRINDEX
    integer :: BINNAME                  ! Name of an l2pc bin
    integer :: BNDPRESSQTYINDEX
    integer :: BNDPRESSVCTRINDEX
    integer :: BOMASK                   ! Set Prec. neg. if which BO in FOV
    integer :: BONUM
    integer :: BQTYINDEX                ! Index of a quantity in vector
    integer :: BVECINDEX                ! Index of a vector
    integer :: BOXCARMETHOD             ! l_min, l_max, l_mean
    real(rv) :: C                       ! constant "c" in manipulation
    logical :: CENTERVERTICALLY         ! For bin based fills
    integer :: CHANNEL                  ! For spreadChannels fill
    integer :: COLMABUNITS              ! l_DOBSONUNITS, or l_MOLCM2
    integer :: COLVECTOR                ! Vector defining columns of Matrix
    type(matrix_SPD_T), pointer :: Covariance
    integer :: DESTINATIONVECTORINDEX   ! For transfer commands
    !                                     -- for FillCovariance
    integer :: EARTHRADIUSQTYINDEX
    integer :: EARTHRADIUSVECTORINDEX
    character(len=256) :: EXTRAOBJECTS  ! Which bright objects to avoid
    integer :: Diagonal                 ! Index of diagonal vector in database
    !                                     -- for FillCovariance
    logical :: DONTMASK                 ! Use even masked values if TRUE
    integer :: ECRTOFOVQUANTITYINDEX    ! Rotation matrix
    integer :: ECRTOFOVVECTORINDEX      ! Rotation matirx
    integer :: ERRORCODE                ! 0 unless error; returned by called routines
    logical :: EXCLUDEBELOWBOTTOM       ! If set in binmax/binmin does not consider stuff below bottom
    character(len=16)  :: explicitUnit  ! E.g., DU
    logical :: Extinction               ! Flag for cloud extinction calculation
    integer :: FIELDINDEX               ! Entry in tree
    integer :: FieldValue               ! Value of a field in the L2CF
    integer :: FIELDECRQUANTITYINDEX    ! Magnetic field
    integer :: FIELDECRVECTORINDEX      ! Magnetic field
    integer :: FILLMETHOD               ! How will we fill this quantity
    integer :: FILENAME                 ! String index for ascii fill
    integer :: FLAGQTYINDEX
    integer :: FLAGVECTORINDEX
    logical :: FORCE                    ! Bypass checks on some operations
    integer :: FRACTION                 ! Index of fraction vector in database
    logical :: FROMPRECISION            ! Fill from l2gpPrecision not l2gpValue
    integer :: GEOCALTITUDEQUANTITYINDEX    ! In the source vector
    integer :: GEOCALTITUDEVECTORINDEX      ! In the vector database
    integer :: GPHQUANTITYINDEX         ! In the source vector
    integer :: GPHVECTORINDEX           ! In the vector database
    logical, dimension(field_first:field_last) :: GOT
    integer :: GRIDINDEX                ! Index of requested grid
    integer :: GSON                     ! Descendant of Son
    integer :: HEIGHTNODE               ! Descendant of son
    logical :: HIGHBOUND                ! Flag
    integer :: H2OQUANTITYINDEX         ! in the quantities database
    integer :: H2OVECTORINDEX           ! In the vector database
    integer :: H2OPRECISIONQUANTITYINDEX         ! in the quantities database
    integer :: H2OPRECISIONVECTORINDEX           ! In the vector database
    integer :: I, J                     ! Loop indices for section, spec, expr
    integer :: GLOBALUNIT               ! To go into the vector
    integer :: IBO
    logical :: IGNOREZERO               ! Don't sum chi^2 at values of noise = 0
    logical :: IGNORENEGATIVE           ! Don't sum chi^2 at values of noise < 0
    logical :: IGNOREGEOLOCATION        ! Don't copy geolocation to vector qua
    integer :: INTERNALVGRIDINDEX       ! Index for internal vgrid (wmoTrop)
    real(r8) :: INTEGRATIONTIME         ! For estimated noise
    logical :: INTERPOLATE              ! Flag for l2gp etc. fill
    integer :: INSTANCESNODE            ! Tree node
    logical :: INVERT                   ! "Invert the specified covariance matrix"
    logical :: ISPRECISION              ! l1b precision, not radiances if TRUE
    integer :: KEY                      ! Definitely n_named
    integer :: L2AUXINDEX               ! Index into L2AUXDatabase
    integer :: L2GPINDEX                ! Index into L2GPDatabase
    integer :: LENGTHSCALE              ! Index of lengthscale vector in database
    logical :: LOGSPACE                 ! Interpolate in log space?
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
    real(r8) :: MAXVALUE                 ! Value of f_maxValue field
    integer :: MAXVALUEUNIT              ! Unit for f_maxValue field
    integer :: MANIPULATION             ! String index
    type(matrix_T), pointer :: MATRIX
    integer :: MATRIXTOFILL             ! Index in database
    integer :: MATRIXTYPE               ! Type of matrix, L_... from init_tables
    integer :: MAXITERATIONS            ! For hydrostatic fill
    integer :: MEASQTYINDEX
    integer :: MEASVECTORINDEX
    character(len=80) :: MESSAGE        ! Possible error message
    integer :: MINNORMQTYINDEX
    integer :: MINNORMVECTORINDEX
    integer :: MODELQTYINDEX
    integer :: MODELVECTORINDEX
    character(len=16)  :: MOL  ! E.g., H2O
    real, dimension(2) :: MULTIPLIER   ! To scale source,noise part of addNoise
    real(r8) :: MINVALUE                 ! Value of f_minValue field
    integer :: MINVALUEUNIT              ! Unit for f_minValue field
    logical :: MISSINGGMAO              ! Only if missing GMAO
    integer :: MULTIPLIERNODE           ! For the parser
    integer :: NBO
    integer :: NBWVECTORINDEX           ! In vector database
    integer :: NBWQUANTITYINDEX         ! In vector database
    integer :: NEEDEDCOORDINATE         ! For vGrid fills
    integer :: NOFINEGRID               ! no of fine grids for cloud extinction calculation
    integer :: NOISEQTYINDEX
    integer :: NOISEVECTORINDEX
    integer :: NORMQTYINDEX
    integer :: NORMVECTORINDEX
    integer :: NOSNOOPEDMATRICES        ! No matrices to snoop
    real(rv) :: OFFSETAMOUNT            ! For offsetRadiance
    logical :: OLD_MATH77_RAN_PACK      ! To restore math77_ran_pack
    integer :: ORBITINCLINATIONVECTORINDEX ! In the vector database
    integer :: ORBITINCLINATIONQUANTITYINDEX ! In the quantities database
    integer :: PHITANVECTORINDEX        ! In the vector database
    integer :: PHITANQUANTITYINDEX      ! In the quantities database
    real(r8) :: PHIWINDOW               ! For hydrostatic ptan guesser
    real(r8) :: PHIZERO                 ! For hydrostatic ptan guesser
    integer :: PHIWINDOWUNITS           ! For hydrostatic ptan guesser
    real(r8) :: PRECISIONFACTOR         ! For setting -ve error bars
    integer :: PRECISIONQUANTITYINDEX   ! For precision quantity
    integer :: PRECISIONVECTORINDEX     ! In the vector database
    integer :: PROFILE                  ! A single profile
    integer :: PTANVECTORINDEX          ! In the vector database
    integer :: PTANQUANTITYINDEX        ! In the quantities database
    logical :: QUADRATURE               ! Apply baseline in quadarture
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
    real(r8) :: SCALE                   ! Scale factor
    real(r8) :: SCALEINSTANCES          ! Scale factor for weighted LS fill
    real(r8) :: SCALERATIO              ! Scale factor for weighted LS fill
    real(r8) :: SCALESURFS              ! Scale factor for weighted LS fill
    integer :: SCECIQUANTITYINDEX       ! In the quantities database
    integer :: SCECIVECTORINDEX         ! In the vector database
    integer :: SCVELQUANTITYINDEX       ! In the quantities database
    integer :: SCVELVECTORINDEX         ! In the vector database
    integer, dimension(2) :: SEED       ! integers used by random_numbers
    integer :: SEEDNODE                 ! For the parser
    logical :: SKIPMASK                 ! Flag for transfer
    integer :: SON                      ! Of root, an n_spec_args or a n_named
    integer :: SOURCE                   ! l_rows or l_colums for adoption
    integer :: SOURCEQUANTITYINDEX      ! in the quantities database
    integer :: SOURCEVECTORINDEX        ! In the vector database
    logical :: SKIPFILL                 ! Don't execute Fill command
    logical :: SPREADFLAG               ! Do we spread values accross instances in explict
    integer :: STATUS                   ! Flag from allocate etc.
    integer :: STATUSVALUE              ! Vaue of f_status
    integer :: SUPERDIAGONAL            ! Index of superdiagonal matrix in database
    logical :: Switch2intrinsic         ! Have mls_random_seed call intrinsic
    !                                     -- for FillCovariance
    real :: T1, T2                      ! for timing
    integer :: SYSTEMPQUANTITYINDEX     ! in the quantities database
    integer :: SYSTEMPVECTORINDEX       ! in the vector database
    integer :: SUFFIX                   ! Possible suffix in reading L1BData
    integer :: TEMPERATUREQUANTITYINDEX ! in the quantities database
    integer :: TEMPERATUREVECTORINDEX   ! In the vector database
    integer :: TEMPPRECISIONQUANTITYINDEX ! in the quantities database
    integer :: TEMPPRECISIONVECTORINDEX   ! In the vector database
    integer :: TEMPLATEINDEX            ! In the template database
    integer :: TERMSNODE                ! Tree index
    logical :: TIMING
    integer :: TNGTECIQUANTITYINDEX     ! In the quantities database
    integer :: TNGTECIVECTORINDEX       ! In the vector database
    integer, dimension(2) :: UNITASARRAY ! From expr
    logical :: UNITSERROR               ! From expr
    integer :: USBVECTORINDEX           ! Inddex in vector database
    integer :: USBQUANTITYINDEX         ! Inddex in vector database
    integer :: USBFRACTIONVECTORINDEX   ! Index in vector database
    integer :: USBFRACTIONQUANTITYINDEX ! Index in vector database
    real(r8), dimension(2) :: VALUEASARRAY ! From expr
    integer :: VALUESNODE               ! For the parser
    integer :: VECTORINDEX              ! In the vector database
    integer :: VECTORNAME               ! Name of vector to create
    integer :: VGRIDINDEX               ! Index of sourceVGrid
    integer :: VMRQTYINDEX
    integer :: VMRQTYVCTRINDEX
    logical :: WHEREFILL                ! Replace only fill values
    logical :: WHERENOTFILL             ! Don't replace fill values
    integer :: WIDTH                    ! Width of boxcar

    ! Executable code
    timing = section_times
    if ( timing ) call time_now ( t1 )
    old_math77_ran_pack = math77_ran_pack

    if ( toggle(gen) ) call trace_begin ( "MLSL2Fill", root )
    if ( specialDumpFile /= ' ' ) &
      & call switchOutput( specialDumpFile, keepOldUnitOpen=.true. )

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
      additional = .false.
      allowMissing = .false.
      boxCarMethod = l_mean
      c = 0.
      channel = 0
      colmabunits = l_molcm2 ! default units for column abundances
      dontMask = .false.
      excludeBelowBottom = .false.
      extinction = .false.
      centerVertically = .false.
      force = .false.
      fromPrecision = .false.
      got= .false.
      MissingGMAO = .false.
      ignoreZero = .false.
      ignoreNegative = .false.
      ignoreGeolocation = .false.
      interpolate = .false.
      invert = .false.
      isPrecision = .false.
      maxValue = huge(0.0_r8)
      maxValueUnit = 0
      minValue = -huge(0.0_r8)
      minValueUnit = 0
      logSpace = .false.
      resetSeed = .false.
      refract = .false.
      scale = 0.0
      scaleInstances = -1.0
      scaleRatio = 1.0
      scaleSurfs = -1.0
      spreadFlag = .false.
      switch2intrinsic = .false.
      suffix = 0
      seed = 0
      noFineGrid = 1
      precisionFactor = 0.5
      offsetAmount = 1000.0             ! Default to 1000K
      profile = -1
      phiWindow = 4
      phiWindowUnits = phyq_angle
      phiZero = 0.0
      quadrature = .false.
      skipFill = .false.
      statusValue = 0
      heightNode = 0
      whereFill = .false.
      whereNotFill = .false.

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
          & qtyTemplates, globalUnit=globalUnit, highBound=highBound, lowBound=lowBound, &
          & where=source_ref(key) ) ) )

        ! That's the end of the create operation

      case ( s_dump ) ! ============================== Dump ==========
        ! Handle disassociated pointers by allocating them with zero size
        status = 0
        if ( .not. associated(qtyTemplates) ) allocate ( qtyTemplates(0), stat=status )
        call test_allocate ( status, moduleName, 'QtyTemplates', (/0/), (/0/) )
        if ( .not. associated(vectorTemplates) ) allocate ( vectorTemplates(0), stat=status )
        call test_allocate ( status, moduleName, 'VectorTemplates', (/0/), (/0/) )
        if ( .not. associated(vectors) ) allocate ( vectors(0), stat=status )
        call test_allocate ( status, moduleName, 'Vectors', (/0/), (/0/) )
        call dumpCommand ( key, qtyTemplates, vectorTemplates, vectors )

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
              & vectors(rowVector), vectors(colVector), where=source_ref(key) )
            call decorate ( key, addToMatrixDatabase(matrices, matrixCholesky) )
          case ( l_kronecker )
            call NullifyMatrix ( matrixKronecker%m )
            call createEmptyMatrix ( matrixKronecker%m, vectorName, &
              & vectors(rowVector), vectors(colVector), where=source_ref(key) )
            call decorate ( key, addToMatrixDatabase(matrices, matrixKronecker) )
          case ( l_plain )
            call NullifyMatrix ( matrixPlain )
            call createEmptyMatrix ( matrixPlain, vectorName, vectors(rowVector), &
              vectors(colVector), where=source_ref(key) )
            call decorate ( key, addToMatrixDatabase(matrices, matrixPlain) )
          case ( l_spd )
            call NullifyMatrix ( matrixSPD%m )
            call createEmptyMatrix ( matrixSPD%m, vectorName, &
              & vectors(colVector), vectors(colVector), where=source_ref(key) )
            call decorate ( key, addToMatrixDatabase(matrices, matrixSPD) )
          end select
        else
          call announce_error ( key, missingField, &
            & extraInfo = (/ f_columns, f_rows /) )
        end if

      case ( s_subset )
        if ( toggle(gen) .and. levels(gen) > 0 ) &
          & call trace_begin ( "Fill.subset", root )
        call SetupSubset ( key, vectors )
        if ( toggle(gen) .and. levels(gen) > 0 ) &
          & call trace_end ( "Fill.subset" )

      case ( s_restrictRange )
        if ( toggle(gen) .and. levels(gen) > 0 ) &
          & call trace_begin ( "Fill.RestrictRange", root )
        call RestrictRange ( key, vectors )
        if ( toggle(gen) .and. levels(gen) > 0 ) &
          & call trace_end ( "Fill.RestrictRange" )

      case ( s_flushL2PCBins )
        call FlushLockedBins

      case ( s_flushPFA )
        call flush_PFAData ( key, status )
        error = max(error,status)

      case ( s_load )
        got = .false.
        do j = 2, nsons(key)
          gson = subtree(j,key)              ! The field
          fieldIndex = get_field_id(gson)
          got(fieldIndex) = .true.
          select case ( fieldIndex )
          case ( f_matrix )
            matrixToFill = decoration(decoration( subtree(2, gson) ))
            if ( getKindFromMatrixDatabase(matrices(matrixToFill)) /= k_plain ) &
              call announce_error ( key, notPlain )
          case ( f_vector )
            vectorIndex = decoration(decoration( subtree(2, gson) ))
          case ( f_bin )
            binName = sub_rosa ( subtree(2,gson) )
          case ( f_source )
            source = decoration ( subtree(2,gson) )
          end select
        end do
        if ( got ( f_matrix ) ) then
          if ( got ( f_vector ) ) call Announce_Error ( key, no_Error_Code, &
            & 'Cannot load both vector and matrix' )
          call GetFromMatrixDatabase ( matrices(matrixToFill), matrix )
          call LoadMatrix ( matrix, binName, message )
        else if ( got ( f_vector ) ) then
          if ( .not. got ( f_source ) ) call Announce_Error ( key, no_Error_Code, &
            & 'Must supply source=rows/columns for vector adoption' )
          call LoadVector ( vectors(vectorIndex), binName, source, message )
        else
          call Announce_Error ( key, no_Error_Code, 'Must supply either matrix or vector to adopt' )
        end if
        if ( len_trim(message) /= 0 ) &
          & call Announce_Error ( key, no_Error_Code, trim(message) )

      case ( s_populateL2PCBin )
        got = .false.
        do j = 2, nsons(key)
          gson = subtree(j,key)              ! The field
          fieldIndex = get_field_id(gson)
          got(fieldIndex) = .true.
          select case ( fieldIndex )
          case ( f_bin )
            call PopulateL2PCBinByName ( sub_rosa ( subtree(2,gson) ) )
          end select
        end do

      case ( s_updateMask )
        if ( toggle(gen) .and. levels(gen) > 0 ) &
          & call trace_begin ( "Fill.UpdateMask", root )
        call UpdateMask ( key, vectors )
        if ( toggle(gen) .and. levels(gen) > 0 ) &
          & call trace_end ( "Fill.UpdateMask" )

      case ( s_flagCloud )
        if ( toggle(gen) .and. levels(gen) > 0 ) &
          & call trace_begin ( "Fill.flagCloud", root )
        call SetupflagCloud ( key, vectors )
        if ( toggle(gen) .and. levels(gen) > 0 ) &
          & call trace_end ( "Fill.flagCloud" )

      case ( s_fill ) ! ===================================  Fill  =====
        ! Now we're on actual Fill instructions.
        ! Loop over the instructions to the Fill command
        BOMask = 0
        AvoidObjects = ' '
        do j = 2, nsons(key)
          gson = subtree(j,key) ! The argument
          fieldIndex = get_field_id(gson)
          if ( nsons(gson) > 1) gson = subtree(2,gson) ! Now value of said argument
          got(fieldIndex)=.TRUE.
          select case ( fieldIndex )
          case ( f_a )
            aVecIndex = decoration(decoration(subtree(1,gson)))
            aQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_additional )
            additional = get_boolean ( gson )
          case ( f_allowMissing )
            allowMissing = get_boolean ( gson )
          case ( f_aprioriPrecision )
            aprPrecVctrIndex = decoration(decoration(subtree(1,gson)))
            aprPrecQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_avoidBrightObjects )
            call get_string( sub_rosa(gson), extraObjects, strip=.true. )
            avoidObjects = catLists( avoidObjects, extraObjects )
          case ( f_b )
            bVecIndex = decoration(decoration(subtree(1,gson)))
            bQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_badRange )
            call expr ( gson , unitAsArray, valueAsArray )
            badRange = valueAsArray
          case ( f_baselineQuantity )
            baselineVctrIndex = decoration(decoration(subtree(1,gson)))
            baselineQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_boxCarMethod )
            boxCarMethod = decoration(gson)
          case ( f_boundaryPressure )     ! For special fill of columnAbundance
            bndPressVctrIndex = decoration(decoration(subtree(1,gson)))
            bndPressQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          case(f_c)
            call expr ( gson , unitAsArray, valueAsArray )
            c = valueAsArray(1)
          case ( f_channel )
            call expr_check ( gson , unitAsArray, valueAsArray, &
              & (/PHYQ_Dimensionless/), unitsError )
            if ( unitsError ) call Announce_error ( subtree(j,key), wrongUnits, &
              & extraInfo=(/unitAsArray(1), PHYQ_Dimensionless/) )
            channel = valueAsArray(1)
          case ( f_centerVertically )
            centerVertically = get_boolean ( gson )
          case ( f_earthRadius ) ! For losGrid fill
            earthRadiusVectorIndex = decoration(decoration(subtree(1,gson)))
            earthRadiusQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_excludeBelowBottom )
            excludeBelowBottom = get_boolean ( gson )
          case ( f_ECRToFOV ) ! For hydrostatic
            ecrToFOVVectorIndex = decoration(decoration(subtree(1,gson)))
            ecrToFOVQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_explicitValues ) ! For explicit fill
            valuesNode = subtree(j,key)
          case ( f_extinction ) ! For cloud extinction fill
            extinction = get_boolean ( gson )
          case ( f_fieldECR ) ! For hydrostatic
            fieldECRVectorIndex = decoration(decoration(subtree(1,gson)))
            fieldECRQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_file ) ! For asciifile
            filename = sub_rosa ( gson )
          case ( f_flags ) ! For chi^2 ratio
            flagVectorIndex = decoration(decoration(subtree(1,gson)))
            flagQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_force )
            force = get_boolean ( gson )
          case ( f_fromPrecision )
            fromPrecision = get_boolean ( gson )
          case ( f_geocAltitudeQuantity ) ! For hydrostatic
            geocAltitudeVectorIndex = decoration(decoration(subtree(1,gson)))
            geocAltitudeQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_gphQuantity ) ! For magnetic field fill
            gphVectorIndex = decoration(decoration(subtree(1,gson)))
            gphQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_height )
            heightNode = subtree(j,key)
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
          case ( f_ignoreGeolocation ) ! For l2gp etc. fill
            ignoreGeolocation =get_boolean ( gson )
          case ( f_ignoreNegative )
            ignoreNegative = get_boolean ( gson )
          case ( f_ifMissingGMAO )
            MissingGMAO = get_boolean ( gson ) .and. &
              & ( APrioriFiles%dao // AprioriFiles%ncep  // AprioriFiles%geos5 &
              &   == ' ' )
          case ( f_instances )
            instancesNode = subtree(j,key)
          case ( f_integrationTime )
            call expr_check ( gson , unitAsArray, valueAsArray, &
              & (/PHYQ_Time/), unitsError )
            if ( unitsError ) call Announce_error ( subtree(j,key), wrongUnits, &
              & extraInfo=(/unitAsArray(1), PHYQ_Time/) )
            integrationTime = valueAsArray(1)
          case ( f_internalVGrid )
            internalVGridIndex=decoration(decoration(gson))
          case ( f_interpolate ) ! For l2gp etc. fill
            interpolate =get_boolean ( gson )
          case ( f_intrinsic )
            switch2intrinsic = get_boolean ( gson )
!         case ( f_invert )
!           invert = get_boolean ( gson )
          case ( f_isPrecision )
            isPrecision = get_boolean ( gson )
          case ( f_logSpace )
            logSpace = get_boolean ( gson )
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
            call expr_check ( gson , unitAsArray, valueAsArray, &
              & (/PHYQ_Dimensionless/), unitsError )
            if ( unitsError ) call Announce_error ( subtree(j,key), wrongUnits, &
              & extraInfo=(/unitAsArray(1), PHYQ_Dimensionless/) )
            maxIterations = valueAsArray(1)
          case ( f_maxValue )      ! For status fill
            call expr ( gson, unitAsArray, valueAsArray )
            maxValueUnit = unitAsArray(1)
            maxValue = valueAsArray(1)
          case ( f_measurements )   ! Only used for diagnostic special fills
            measVectorIndex = decoration(decoration(subtree(1,gson)))
            measQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_method )   ! How are we going to fill it?
            fillMethod = decoration(gson)
          case ( f_minNormQty )   ! Only used for chi^2 ratio fills
            minNormVectorIndex = decoration(decoration(subtree(1,gson)))
            minNormQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_minValue )      ! For status fill
            call expr ( gson, unitAsArray, valueAsArray )
            minValueUnit = unitAsArray(1)
            minValue = valueAsArray(1)
          case ( f_model )   ! Only used for diagnostic special fills
            modelVectorIndex = decoration(decoration(subtree(1,gson)))
            modelQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_multiplier ) ! For scaling noise part of addnoise
            multiplierNode=subtree(j,key)
          case ( f_noFineGrid )      ! For cloud extinction fill
            call expr_check ( gson , unitAsArray, valueAsArray, &
              & (/PHYQ_Dimensionless/), unitsError )
            if ( unitsError ) call Announce_error ( subtree(j,key), wrongUnits, &
              & extraInfo=(/unitAsArray(1), PHYQ_Dimensionless/) )
            noFineGrid = valueAsArray(1)
          case ( f_noise )   ! Only used for chi^2 special fills or addnoise
            noiseVectorIndex = decoration(decoration(subtree(1,gson)))
            noiseQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_noiseBandwidth )
            nbwVectorIndex = decoration(decoration(subtree(1,gson)))
            nbwQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_normQty )   ! Only used for chi^2 ratio fills
            normVectorIndex = decoration(decoration(subtree(1,gson)))
            normQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_offsetAmount )    ! For marking unused radiances
            call expr_check ( gson , unitAsArray, valueAsArray, &
              & (/PHYQ_Temperature/), unitsError )
            if ( unitsError ) call Announce_error ( subtree(j,key), wrongUnits, &
              & extraInfo=(/unitAsArray(1), PHYQ_Temperature/) )
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
            call expr_check ( gson , unitAsArray, valueAsArray, &
              & (/PHYQ_Dimensionless/), unitsError )
            if ( unitsError ) call Announce_error ( subtree(j,key), wrongUnits, &
              & extraInfo=(/unitAsArray(1), PHYQ_Dimensionless/) )
            precisionFactor = valueAsArray(1)
          case ( f_profile )
            call expr_check ( gson , unitAsArray, valueAsArray, &
              & (/PHYQ_Dimensionless/), unitsError )
            if ( unitsError ) call Announce_error ( subtree(j,key), wrongUnits, &
              & extraInfo=(/unitAsArray(1), PHYQ_Dimensionless/) )
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
            call expr_check ( gson , unitAsArray, valueAsArray, &
              & (/ PHYQ_Profiles, PHYQ_Angle /), unitsError )
            if ( unitsError ) call Announce_error ( subtree(j,key), wrongUnits, &
              & extraInfo=(/unitAsArray(1), PHYQ_Profiles, PHYQ_Angle/) )
            phiWindow = valueAsArray(1)
            phiWindowUnits = unitAsArray(1)
          case ( f_phiZero )
            call expr_check ( gson , unitAsArray, valueAsArray, &
              & (/PHYQ_Angle/), unitsError )
            if ( unitsError ) call Announce_error ( subtree(j,key), wrongUnits, &
              & extraInfo=(/unitAsArray(1), PHYQ_Angle/) )
            phiZero = valueAsArray(1)
          case ( f_quadrature )
            quadrature = get_boolean ( gson )
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
          case ( f_scale, f_scaleInsts, f_scaleRatio, f_scaleSurfs )
            call expr_check ( gson , unitAsArray, valueAsArray, &
              & (/PHYQ_Dimensionless/), unitsError )
            if ( unitsError ) call Announce_error ( subtree(j,key), wrongUnits, &
              & extraInfo=(/unitAsArray(1), PHYQ_Dimensionless/) )
            select case ( fieldIndex )
            case ( f_scale )
              scale = valueAsArray(1)
            case ( f_scaleInsts )
              scaleInstances = valueAsArray(1)
            case ( f_scaleRatio )
              scaleRatio = valueAsArray(1)
            case ( f_scaleSurfs )
              scaleSurfs = valueAsArray(1)
            end select
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
            spreadFlag = get_boolean ( gson )
          case ( f_status )
            valuesNode = subtree(j,key)
            call expr_check ( gson , unitAsArray, valueAsArray, &
              & (/PHYQ_Dimensionless/), unitsError )
            if ( unitsError ) call Announce_error ( valuesNode, wrongUnits, &
              & extraInfo=(/unitAsArray(1), PHYQ_Dimensionless/) )
            statusValue = nint ( valueAsArray(1) )
          case ( f_suffix )
            suffix = sub_rosa ( gson )
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
          case ( f_terms )
            termsNode = subtree(j,key)
          case ( f_unit ) ! For folding
            colmabunits = decoration(gson) ! decoration(subtree(2,gson))
          case ( f_usb ) ! For folding
            usbVectorIndex = decoration(decoration(subtree(1,gson)))
            usbQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_usbFraction ) ! For folding
            usbFractionVectorIndex = decoration(decoration(subtree(1,gson)))
            usbFractionQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_vmrQuantity )     ! For special fill of columnAbundance
            vmrQtyVctrIndex = decoration(decoration(subtree(1,gson)))
            vmrQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
          case ( f_whereFill )
            whereFill = get_boolean ( gson )
          case ( f_whereNotFill )
            whereNotFill = get_boolean ( gson )
          case ( f_width )
            call expr_check ( gson , unitAsArray, valueAsArray, &
              & (/PHYQ_Dimensionless/), unitsError )
            if ( unitsError ) call Announce_error ( subtree(j,key), wrongUnits, &
              & extraInfo=(/unitAsArray(1), PHYQ_Dimensionless/) )
            width = valueAsArray(1)
          end select
        end do                  ! Loop over arguments to fill instruction

        ! Put various conditions under which you would want to skip this fill
        if ( got(f_ifMissingGMAO) .and. .not. MissingGMAO ) skipFill = .true.
        if ( skipFill ) fillMethod = -1 ! We'll assume no l_value can be this

        ! Now call various routines to do the filling
        quantity => GetVectorQtyByTemplateIndex( &
          & vectors(vectorIndex), quantityIndex )

        select case ( fillMethod )
        case ( l_addNoise ) ! ----- Add random noise to source Quantity -------
          if ( DEEBUG) call output('add noise method', advance='yes')
          if ( .not. all(got( (/f_Quantity, f_sourceQuantity, f_noise/) ) ) ) &
            call Announce_error ( key, No_Error_code, &
             'Missing a required field to add noise'  )
          Quantity => GetVectorQtyByTemplateIndex( &
            & vectors(VectorIndex), QuantityIndex )
          sourceQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(sourceVectorIndex), sourceQuantityIndex )
          noiseQty => GetVectorQtyByTemplateIndex( &
            & vectors(noiseVectorIndex), noiseQtyIndex)
          math77_ran_pack = .not. switch2intrinsic
          if ( DEEBUG ) then
            call output('Switch to intrinsic? ', advance='no')
            call output(switch2intrinsic, advance='yes')
            call output('resetSeed? ', advance='no')
            call output(resetSeed, advance='yes')
            call output('got(f_seed)? ', advance='no')
            call output(got(f_seed), advance='yes')
          end if
          if ( resetSeed ) then
            call mls_random_seed(new_seed=seed(1:))
            if ( DEEBUG ) then
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
              if ( DEEBUG ) then
                call output('Setting new seed ', advance='no')
                call output(seed, advance='yes')
              end if
            else
              call mls_random_seed(new_seed=seed(1:))
              if ( DEEBUG ) then
                call output('Letting mls choose new seed ', advance='no')
                call output(seed, advance='yes')
              end if
            end if
          else
            if ( DEEBUG ) then
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

          if ( DEEBUG ) then
            call output('Using multipliers: ', advance='no')
            call output(multiplier, advance='yes')
          end if
          call addGaussianNoise ( key, quantity, sourceQuantity, &
            & noiseQty, multiplier )

        case ( l_applyBaseline )
          if ( .not. got ( f_baselineQuantity ) ) &
            & call Announce_Error ( key, no_Error_Code, &
            & 'Need baselineQuantity for applyBaseline fill' )
          baselineQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(baselineVctrIndex), baselineQtyIndex )
          call ApplyBaseline ( key, quantity, baselineQuantity, &
            & quadrature, dontmask )

        case ( l_asciiFile )
          if ( .not. got ( f_file ) ) &
            & call Announce_Error ( key, no_Error_Code, &
            & 'Need filename for asciiFile fill' )
          if ( got ( f_badRange ) ) then
            call FillQuantityFromASCIIFile ( key, quantity, filename, badRange )
          else
            call FillQuantityFromASCIIFile ( key, quantity, filename )
          end if

        case ( l_binMax, l_binMean, l_binMin, l_binTotal, &
             & l_lsLocal, l_lsGlobal, l_lsWeighted )
          if ( .not. got ( f_sourceQuantity ) ) &
            & call Announce_Error ( key, no_Error_Code, &
            & 'Need source quantity for bin fill or least-squares fill' )
          sourceQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(sourceVectorIndex), sourceQuantityIndex )
          if ( got ( f_ptanQuantity ) ) then
            ptanQuantity => GetVectorQtyByTemplateIndex( &
              & vectors(ptanVectorIndex), ptanQuantityIndex)
          else
            nullify ( ptanQuantity )
          end if

          if ( sourceQuantity%template%verticalCoordinate /= &
            & quantity%template%verticalCoordinate ) then
            if ( .not. sourceQuantity%template%minorFrame .or. &
              &  .not. got ( f_ptanQuantity ) .or. &
              &  quantity%template%verticalCoordinate /= l_zeta ) then
                call Announce_Error ( key, no_Error_Code, &
                & ' vertical coordinates mismatch, perhaps supply tangent pressure?' )
            else
              if ( ptanQuantity%template%instrumentModule .ne. &
                & sourceQuantity%template%instrumentModule ) &
                & call Announce_Error ( key, no_Error_Code, &
                & 'Instrument module mismatch between ptan and source quantity' )
            end if
          end if

          if ( error == 0 ) then
            select case ( fillMethod )
            case ( l_binMax, l_binMean, l_binMin, l_binTotal )
              call FillWithBinResults ( key, quantity, sourceQuantity, ptanQuantity, &
                & channel, fillMethod, additional, excludeBelowBottom, centerVertically )
            case default
              call FillUsingLeastSquares ( key, quantity, sourceQuantity, ptanQuantity, &
                & channel, fillMethod, scaleInstances, scaleRatio, scaleSurfs )
            end select
          end if

        case ( l_boxcar )
          if ( .not. got ( f_sourceQuantity ) ) &
            & call Announce_Error ( key, no_Error_Code, &
            & 'Need source quantity for boxcar fill' )
          sourceQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(sourceVectorIndex), sourceQuantityIndex )
          if ( .not. got ( f_width ) ) call Announce_Error ( key, no_Error_Code, &
            & 'Must supply width for boxcar fill' )
          call FillWithBoxcarFunction ( key, quantity, sourceQuantity, width, &
            & boxCarMethod )

        case ( l_chiSqRatio ) ! ----------- Fill with convergence ratio ---
          if ( .not. all(got( (/ f_normQty, &
            & f_minNormQty, f_flags /)))) &
            & call Announce_Error ( key, no_Error_Code, &
              & 'Missing required fields for chi^2 ratio' )
          normQty => GetVectorQtyByTemplateIndex( &
            & vectors(normVectorIndex), normQtyIndex )
          minNormQty => GetVectorQtyByTemplateIndex( &
            & vectors(minNormVectorIndex), minNormQtyIndex )
          flagQty => GetVectorQtyByTemplateIndex( &
            & vectors(flagVectorIndex), flagQtyIndex )
          call FillChiSqRatio ( key, &
            & quantity, normQty, minNormQty, flagQty, dontMask )

        case ( l_combineChannels )
          if ( .not. got ( f_sourceQuantity ) ) &
            & call Announce_Error ( key, no_Error_Code, &
            & 'Need source quantity for combine channels fill' )
          sourceQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(sourceVectorIndex), sourceQuantityIndex )
          call FillWithCombinedChannels ( quantity, sourceQuantity, message )
          if ( message /= '' ) call Announce_Error ( key, no_Error_Code, trim(message) )

        case ( l_convergenceRatio )
          if ( .not. all ( got ( (/ f_sourceQuantity, f_scale /) ) ) ) &
            call Announce_error ( key, no_error_code, &
            & 'Need sourceQuanitity and scale for quality fill' )
          sourceQuantity => GetVectorQtyByTemplateIndex ( vectors(sourceVectorIndex), &
            & sourceQuantityIndex )
          call FillConvergenceFromChisq ( key, quantity, sourceQuantity, scale )

        case ( l_estimatedNoise ) ! ----------- Fill with estimated noise ---
          if ( .not. all(got( (/ f_radianceQuantity, &
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

        case ( l_explicit ) ! ---------  Explicitly fill from l2cf  -----
          if ( .not. got(f_explicitValues) ) &
            & call Announce_Error ( key, noExplicitValuesGiven )
          call ExplicitFillVectorQuantity ( quantity, valuesNode, spreadFlag, &
            & vectors(vectorIndex)%globalUnit, dontmask )

        case ( l_extractChannel )
          if ( .not. all(got ( (/f_sourceQuantity,f_channel/)))) &
            & call Announce_Error ( key, no_Error_Code, &
              & 'Need sourceQuantity and channel for this fill' )
          sourceQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(sourceVectorIndex), sourceQuantityIndex )
          call ExtractSingleChannel ( key, quantity, sourceQuantity, channel )

        case ( l_fold ) ! --------------- Fill by sideband folding -----
          nullify ( lsb, lsbFraction, usb, usbFraction )
          if ( got ( f_lsb ) ) lsb => GetVectorQtyByTemplateIndex ( &
            & vectors(lsbVectorIndex), lsbQuantityIndex )
          if ( got ( f_lsbFraction ) ) lsbFraction => GetVectorQtyByTemplateIndex ( &
            & vectors(lsbFractionVectorIndex), lsbFractionQuantityIndex )
          if ( got ( f_usb ) ) usb => GetVectorQtyByTemplateIndex ( &
            & vectors(usbVectorIndex), usbQuantityIndex )
          if ( got ( f_usbFraction ) ) usbFraction => GetVectorQtyByTemplateIndex ( &
            & vectors(usbFractionVectorIndex), usbFractionQuantityIndex )
          call FillFoldedRadiance ( quantity, lsb, usb, lsbFraction, usbFraction, key )

        case ( l_fwdModelTiming ) ! --- Fill timings for forward model  -----
          call FillFwdModelTimings (quantity%values(:,1), FWModelConfig, 'fwdTiming')
        case ( l_fwdModelMean ) ! --- Fill mean for forward model  -----
          call FillFwdModelTimings (quantity%values(:,1), FWModelConfig, 'mean')
        case ( l_fwdModelStdDev ) ! --- Fill std_dev for forward model  -----
          call FillFwdModelTimings (quantity%values(:,1), FWModelConfig, 'stdDev')
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

        case ( l_gridded ) ! ------------  Fill from gridded data  -----
          if ( .not. got(f_sourceGrid) ) &
            & call Announce_Error ( key, noSourceGridGiven )
          ! call output( 'Filling quantity from grid', advance='yes' )
          call FillVectorQuantityFromGrid &
            & ( quantity, griddedDataBase(gridIndex), allowMissing, errorCode )
          ! call output_name_v_pair( 'error code', errorCode )
          if ( errorCode /= 0 ) call Announce_error ( key, errorCode )

        case ( l_l1b ) ! --------------------  Fill from L1B data  -----
          if ( got(f_precision) ) then
            precisionQuantity => GetVectorQtyByTemplateIndex( &
              & vectors(precisionVectorIndex), precisionQuantityIndex )
            call FillVectorQuantityFromL1B ( key, quantity, chunks(chunkNo), &
              & filedatabase, isPrecision, suffix=suffix, &
              & precisionQuantity=precisionQuantity )
          elseif ( got(f_avoidbrightobjects) ) then
            avoidObjects = lowerCase(avoidObjects)
            nBO = NumStringElements( avoidObjects, .true. )
            do iBO = 1, nBO
              BOnum = StringElementNum(lowercase(BrightObjects), &
                & trim(StringElement(avoidObjects, iBO, .true.)), &
                & .true.)
              if ( BOnum > 0 ) BOMask = ibset( BOMask, BOnum )
            end do
            ! Special case: moon in space port
            if ( index(avoidObjects, 'mooninsp') > 0 ) &
              & BOMask = ibset( BOMask, 0 )
            call FillVectorQuantityFromL1B ( key, quantity, chunks(chunkNo), &
              & filedatabase, isPrecision, suffix=suffix, BOMask=BOMask )
          else
            call FillVectorQuantityFromL1B ( key, quantity, chunks(chunkNo), &
              & filedatabase, isPrecision, suffix=suffix )
          end if

        case ( l_l2gp ) ! --------------  Fill from L2GP quantity  -----
          if ( .NOT. got(f_sourceL2GP) ) &
            & call Announce_Error ( key, noSourceL2GPGiven )
          call FillVectorQuantityFromL2GP &
            & ( quantity, l2gpDatabase(l2gpIndex), interpolate, profile, errorCode, &
            & ignoreGeolocation, fromPrecision  )
          if ( errorCode /= 0 ) call Announce_error ( key, errorCode )

        case ( l_l2aux ) ! ------------  Fill from L2AUX quantity  -----
          if ( .NOT. got(f_sourceL2AUX) ) &
            & call Announce_Error ( key, noSourceL2AUXGiven )
          call FillVectorQuantityFromL2AUX(quantity,l2auxDatabase(l2auxIndex),errorCode)
          if ( errorCode /= 0 ) call Announce_error ( key, errorCode )

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
                & quantityType=(/l_rhi, l_vmr/)) ) then
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

        case ( l_hydrostatic ) ! -------------  Hydrostatic fills  -----
          ! Need a temperature and a refgph quantity
          if ( .not.all(got( (/ f_refGPHQuantity, f_temperatureQuantity /))) ) &
            call Announce_Error ( key, needTempREFGPH )

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
              & call Announce_Error ( key, no_Error_Code, 'Needs phiTan quantity' )
            phiTanQuantity => GetVectorQtyByTemplateIndex( &
              & vectors(phiTanVectorIndex), phiTanQuantityIndex)
            if ( phiTanQuantity%template%quantityType /= l_phiTan ) &
              & call Announce_Error ( key, no_Error_Code, 'Has a bad phiTan quantity' )

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
            & phiWindow, phiWindowUnits, chunkNo )

        case ( l_isotope ) ! --------------- Isotope based fills -------
          if ( .not. all(got( (/f_ratioQuantity, f_sourceQuantity/) ) ) ) &
            & call Announce_Error ( key, badIsotopeFill )
          ratioQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(ratioVectorIndex), ratioQuantityIndex )
          sourceQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(sourceVectorIndex), sourceQuantityIndex )
          call FillVectorQtyFromIsotope ( quantity, sourceQuantity, &
            & ratioQuantity )

        case ( l_IWCfromExtinction ) ! -------fill IWC from cloud extinction -------
            if ( .not. any(got( &
             & (/f_extinction, f_temperatureQuantity/) &
             & )) ) then
              call Announce_error ( key, No_Error_code, &
              & 'Missing a required field to fill iwc from cloudextinction'  )
            else
              sourceQuantity => GetVectorQtyByTemplateIndex( &
                & vectors(sourceVectorIndex), sourceQuantityIndex)
              temperatureQuantity => GetVectorQtyByTemplateIndex( &
                & vectors(temperatureVectorIndex), temperatureQuantityIndex)
              if ( .not. ValidateVectorQuantity(sourceQuantity, &
                & quantityType=(/l_cloudextinction/)) ) then
                call Announce_Error ( key, No_Error_code, &
                & 'The extinctionQuantity is not an cloudextinction'  )
              else if ( .not. ValidateVectorQuantity(Quantity, &
                & quantityType=(/l_cloudice/)) ) then
                call Announce_Error ( key, No_Error_code, &
                & 'The filled Quantity is not a type of cloudice '  )
              else if ( .not. ValidateVectorQuantity(temperatureQuantity, &
                & quantityType=(/l_temperature/)) ) then
                call Announce_Error ( key, No_Error_code, &
                & 'The temperatureQuantity is not a temperature'  )
              else
                call FillIWCFromExtinction ( quantity, &
                  & sourceQuantity, temperatureQuantity)
              end if
            end if

        case ( l_manipulate ) ! ---------------------------- Manipulate --
          if ( .not. got ( f_a ) ) &
            & call Announce_error ( key, no_Error_Code,'aQuantity not supplied' )
          if ( .not. got ( f_manipulation ) ) &
            & call Announce_error ( key, no_Error_Code,'manipulation not supplied' )
          aQuantity => GetVectorQtyByTemplateIndex ( &
            & vectors(aVecIndex), aQtyIndex )
          if ( got ( f_b ) ) then
            bQuantity => GetVectorQtyByTemplateIndex ( &
              & vectors(bVecIndex), bQtyIndex )
          else
            nullify ( bQuantity )
          end if
          call FillQuantityByManipulation ( quantity, aQuantity, bQuantity, &
            & manipulation, key, force, c )

        case ( l_magAzEl ) ! -- Magnetic Explicit from stren, azim, elev --
          if ( .not. got(f_explicitValues) ) &
            & call Announce_Error ( key, noExplicitValuesGiven )
          call ExplicitFillVectorQuantity ( quantity, valuesNode, spreadFlag, &
            & vectors(vectorIndex)%globalUnit, dontmask, azEl=.true. )

        case ( l_magneticModel ) ! --------------------- Magnetic Model --
          if ( .not. got ( f_gphQuantity ) ) then
            call Announce_Error ( key, no_Error_Code, 'Need gphQuantity for magnetic model' )
          else
            GPHQuantity => GetVectorQtyByTemplateIndex( &
              & vectors(GPHVectorIndex), GPHQuantityIndex)
            call FillQuantityUsingMagneticModel ( quantity, gphQuantity, key )
          end if

        case ( l_negativePrecision ) ! ------------ Set output SD -ve wrt apriori
          if ( .not. got ( f_aprioriPrecision ) ) &
            & call Announce_Error ( key, No_Error_code, &
            & 'Missing aprioriPrecision field for negativePrecision fill' )
          aprioriPrecision => GetVectorQtyByTemplateIndex ( &
            & vectors(aprPrecVctrIndex), aprPrecQtyIndex )
          where ( quantity%values >= aprioriPrecision%values*precisionFactor .and. &
            & aprioriPrecision%values > 0.0_r8 )
            quantity%values = - quantity%values
          end where

        case ( l_offsetRadiance ) ! ------------------- Offset radiance --
          if ( .not. got ( f_radianceQuantity ) ) &
            & call Announce_error ( key, no_Error_Code, 'radianceQuantity not supplied' )
          radianceQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(radianceVectorIndex), radianceQuantityIndex )
          call OffsetRadianceQuantity ( quantity, radianceQuantity, offsetAmount )

        case ( l_phaseTiming ) ! ---------  Fill timings for phases  -----
          call finishTimings('phases', returnStatus=status)
          if ( status /= 0 ) then
            call MLSMessage ( MLSMSG_Warning, ModuleName, &
              & 'Unable to finish phases timings' )
          else
            call FillTimings ( quantity%values(:,1), 'phases', 'all', .true. )
            ! call dump( quantity%values(:,1), 'phases' )
          end if

        case ( l_sectionTiming ) ! ---------  Fill timings for sections  -----
          call finishTimings('sections', returnStatus=status)
          if ( status /= 0 ) then
            call MLSMessage ( MLSMSG_Warning, ModuleName, 'Unable to finish sections timings' )
          else
            call FillTimings ( quantity%values(:,1), 'sections', 'all', .true. )
            ! call dump( quantity%values(:,1), 'sections' )
          end if

        case ( l_profile ) ! ------------------------ Profile fill -------
          if ( .not. got ( f_profileValues ) ) &
            call Announce_error ( key, no_Error_Code, 'profileValues not supplied' )
          if ( .not. got ( f_instances ) ) instancesNode = 0
          if ( got ( f_logSpace ) ) then
            call FillVectorQtyFromProfile ( quantity, valuesNode, &
              & instancesNode, vectors(vectorIndex)%globalUnit, dontMask, logSpace=logSpace )
          else
            call FillVectorQtyFromProfile ( quantity, valuesNode, &
              & instancesNode, vectors(vectorIndex)%globalUnit, dontMask )
          end if

        case ( l_refract )              ! --------- refraction for phiTan -----
          ! More sanity checks
          if ( .not. ValidateVectorQuantity ( quantity, &
            & quantityType=(/l_phiTan/), minorFrame=.true. ) ) &
            & call Announce_error ( key, no_Error_Code, 'Quantity to fill is not phiTan' )
          ! Start off with a copy from the template
          quantity%values = quantity%template%phi
          if ( refract ) then
            if ( .not. all ( got ( (/ f_h2oQuantity, f_orbitinclination, &
            & f_ptanQuantity, f_refGPHquantity, f_temperatureQuantity /) ) ) ) then
              call Announce_error ( key, badRefractFill )
              call announce_error ( key, missingField, extraInfo= &
                & (/ f_h2oQuantity, f_orbitinclination, &
                & f_ptanQuantity, f_refGPHquantity, f_temperatureQuantity /) )
            end if
            h2oQuantity => GetVectorQtyByTemplateIndex( &
              & vectors(h2oVectorIndex), h2oQuantityIndex)
            orbitInclinationQuantity => GetVectorQtyByTemplateIndex( &
              & vectors(orbitInclinationVectorIndex), orbitInclinationQuantityIndex)
            ptanQuantity => GetVectorQtyByTemplateIndex( &
              & vectors(ptanVectorIndex), ptanQuantityIndex)
            refGPHquantity => GetVectorQtyByTemplateIndex( &
              & vectors(refGPHVectorIndex), refGPHQuantityIndex)
            temperatureQuantity => GetVectorQtyByTemplateIndex( &
              & vectors(temperatureVectorIndex), temperatureQuantityIndex)
            call FillPhiTanWithRefraction ( key, quantity, h2oQuantity, &
              & orbitInclinationQuantity, ptanQuantity, refGPHquantity, temperatureQuantity )
          end if

        case ( l_reflectorTempModel ) ! --------------- Reflector temperature model
          call FillWithReflectorTemperature ( key, quantity, phiZero, termsNode )

        case ( l_resetUnusedRadiances )
          call ResetUnusedRadiances ( quantity, offsetAmount )

        case ( l_rectanglefromlos ) ! -------fill from losGrid quantity -------
          if ( .not. all(got((/f_losQty,f_earthRadius,f_PtanQuantity/))))&
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

        case ( l_quality )
          if ( .not. all ( got ( (/ f_sourceQuantity, f_scale /) ) ) ) &
            call Announce_error ( key, no_error_code, &
            & 'Need sourceQuanitity and scale for quality fill' )
          sourceQuantity => GetVectorQtyByTemplateIndex ( vectors(sourceVectorIndex), &
            & sourceQuantityIndex )
          call FillQualityFromChisq ( key, quantity, sourceQuantity, scale, heightNode )

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
            else
              call FillRHIPrecisionFromH2O ( key, quantity, &
                & h2oPrecisionQuantity, tempPrecisionQuantity, h2oQuantity, &
                & temperatureQuantity, dontMask, ignoreZero, &
                & ignoreNegative, interpolate, &
                & .true., &   ! Mark Undefined values?
                & invert )    ! invert rather than convert?
            end if
          end if

        case ( l_rotateField )
          if ( .not. all ( got ( (/ f_fieldECR, f_ecrtofov /) ) ) ) then
            call Announce_Error ( key, no_error_code, &
              & 'Must supply field and ecrToFov for rotateField fill' )
          else
            fieldECR => GetVectorQtyByTemplateIndex ( &
              & vectors(fieldECRVectorIndex), fieldECRQuantityIndex )
            ecrToFOV => GetVectorQtyByTemplateIndex ( &
              & vectors(ecrToFOVVectorIndex), ecrToFovQuantityIndex )
            call RotateMagneticField ( key, quantity, fieldECR, ecrToFov )
          end if

        case ( l_scaleOverlaps )
          if ( .not. got ( f_multiplier ) ) then
            call Announce_Error ( key, no_error_code, &
              & 'Must supply multipler for scaleOverlaps fill' )
          else
            call ScaleOverlaps ( quantity, multiplierNode, dontMask )
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

          if ( DEEBUG ) then
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
            elseif ( .not. &
              & any( colmabunits == (/l_dobsonUnits, l_DU, l_molcm2/) ) ) then
              call Announce_error ( key, No_Error_code, &
              & 'Wrong units to fill column abundance'  )
            else
              bndPressQty => GetVectorQtyByTemplateIndex( &
                & vectors(bndPressVctrIndex), bndPressQtyIndex)
              vmrQty => GetVectorQtyByTemplateIndex( &
                & vectors(vmrQtyVctrIndex), vmrQtyIndex)
              if ( got(f_unit) ) then
                ! Switch column species hash to explicit unit
                call get_string( lit_indices(colmabunits), explicitUnit, &
                  & strip=.true. )
                call get_string( lit_indices(quantity%template%molecule), mol, &
                  & strip=.true. )
                call PutHashElement( col_species_keys, col_species_hash, &
                  & lowerCase(mol), ExplicitUnit, countEmpty=.true. )
                ! print *,'col_species_keys: '
                ! print *,col_species_keys
                ! print *,'col_species_hash: '
                ! print *,col_species_hash
                ! call dump( .true., col_species_keys, col_species_hash, &
                !   & 'column species units' )
              else
                call get_string( lit_indices(quantity%template%molecule), mol, &
                  & strip=.true. )
                call GetHashElement (col_species_keys, &
                  & col_species_hash, trim(lowercase(mol)), &
                  & ExplicitUnit, .true.)
                if ( index(lowerCase(ExplicitUnit), 'd') > 0 ) colmabunits = l_DU
              end if
              ! print *, 'species, column unit: ', mol, ExplicitUnit
              call FillColAbundance ( key, quantity, &
                & bndPressQty, vmrQty, colmAbUnits )
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
          case ( l_noRadsPerMIF )
            if ( .not. got ( f_measurements ) ) then
              call Announce_error ( key, No_Error_code, &
              & 'Missing a required field to fill noRadsPerMIF on MIFs'  )
            else
              measQty => GetVectorQtyByTemplateIndex( &
                & vectors(measVectorIndex), measQtyIndex)
              call FillNoRadsPerMif ( key, quantity, measQty )
            end if
          case default
            call Announce_error ( key, noSpecialFill )
          end select ! quantity types in special fill cases

        case ( l_splitSideband ) ! --------------- Split the sidebands
          if ( .not. got(f_sourceQuantity) ) &
            & call Announce_Error ( key, No_Error_Code, &
            & 'Missing a source field for vector fill' )
          if ( .not. all(got( (/f_lsbFraction,f_usbFraction/) ))) &
            & call Announce_Error ( key, No_Error_Code, &
            & 'Missing a usb/lsb fraction field for vector fill' )
          if ( .not. got ( f_channel ) ) call Announce_Error ( key, &
            & no_error_code, 'Must supply channel for spreadChannel fill' )
          sourceQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(sourceVectorIndex), sourceQuantityIndex )
          lsbFraction => GetVectorQtyByTemplateIndex ( &
            & vectors(lsbFractionVectorIndex), lsbFractionQuantityIndex )
          usbFraction => GetVectorQtyByTemplateIndex ( &
            & vectors(usbFractionVectorIndex), usbFractionQuantityIndex )
          if ( got ( f_usb ) ) then
               usb => GetVectorQtyByTemplateIndex ( &
               & vectors(usbVectorIndex), usbQuantityIndex )
          else
            nullify ( usb )
          end if

          call FillFromSplitSideband ( quantity, sourceQuantity, &
            & lsbFraction, usbFraction, spreadFlag, usb, channel, key )

        case ( l_spreadChannel )
          if ( .not. got ( f_channel ) .and. .not. got( f_sourceQuantity ) ) &
            & call Announce_Error ( key, &
            & no_error_code, 'Must supply channel or sourcequantity for spreadChannel fill' )
          if ( got(f_sourceQuantity) ) then
          sourceQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(sourceVectorIndex), sourceQuantityIndex )
            if ( .not. got(f_channel) ) channel = 1
            call SpreadChannelFill ( quantity, channel, dontMask, key, &
              & sourceQuantity )
          else
            call SpreadChannelFill ( quantity, channel, dontMask, key )
          endif

        case ( l_status )
          if ( got(f_ifMissingGMAO) ) then
            if ( MissingGMAO ) call ExplicitFillVectorQuantity ( quantity, &
              & valuesNode, .true., phyq_Invalid, .true., options='-v' )
          elseif ( .not. all ( got ( (/ f_sourceQuantity, f_status /) ) ) ) then
            call Announce_Error ( key, no_error_code, &
            & 'Need sourceQuantity and status fields for status fill' )
          elseif ( .not. any ( got ( (/ f_minValue, f_maxValue /) ) ) ) then
            call Announce_Error ( key, no_error_code, &
            & 'Need one or both of maxValue, minValue for status fill' )
          else
            sourceQuantity => &
              &GetVectorQtyByTemplateIndex ( vectors(sourceVectorIndex), &
              & sourceQuantityIndex )
            if ( got ( f_maxValue) .and. &
              &  all ( maxValueUnit /= &
              & (/ sourceQuantity%template%unit, PHYQ_Dimensionless/) ) ) &
              & call Announce_Error ( key, no_error_code, &
              & 'Bad unit for maxValue' )
            if ( got ( f_minValue) .and. &
              &  all ( minValueUnit /= &
              & (/ sourceQuantity%template%unit, PHYQ_Dimensionless/) ) ) &
              & call Announce_Error ( key, no_error_code, &
              & 'Bad unit for minValue' )
            if ( all ( got ( (/ f_minValue, f_maxValue /) ) ) .and. &
              &  maxValue <= minValue ) call Announce_Error ( key, no_error_code, &
              & 'Bad combination of max/min values' )
            call FillStatusQuantity ( key, quantity, &
              & sourceQuantity, statusValue, &
              & minValue, maxValue, heightNode, additional )
          end if

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
            if ( .not. interpolate .and. .not. force ) then
              call Announce_Error ( key, No_Error_Code, &
                & 'Quantity and sourceQuantity do not have the same template' )
            else
              call FillQtyFromInterpolatedQty ( quantity, sourceQuantity, force, key )
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

        case ( l_vGrid ) ! ---------------------  Fill from vGrid  -----
          select case ( quantity%template%quantityType )
          case ( l_ptan )
            neededCoordinate = l_zeta
          case ( l_tngtGeodAlt )
            neededCoordinate = l_geodAltitude
          case ( l_tngtGeocAlt )
            neededCoordinate = l_geocAltitude
          end select
          if ( vGrids(vGridIndex)%verticalCoordinate /= neededCoordinate ) &
            & call Announce_Error ( key, No_Error_code, &
            &  'Vertical coordinate vGrid does not match fill' )
          if ( vGrids(vGridIndex)%noSurfs /= quantity%template%noSurfs )&
            & call Announce_Error ( key, No_Error_code, &
            &  'VGrid is not of the same size as the quantity' )
          quantity%values = spread ( vGrids(vGridIndex)%surfs(:,1), 2, &
            & quantity%template%noInstances )

        case ( l_wmoTropopause ) ! ---------------- Fill with wmo tropopause --
!           if ( .not. all(got( (/ f_temperatureQuantity, f_refGPHquantity, &
!             & f_internalVGrid /) ))) &
!             & call Announce_Error ( key, no_error_code, &
!             & 'wmoTropopause fill needs temperatureQuantity, refGPHQuantity ' // &
!             & 'and internalVGrid' )
!           if ( vGrids(internalVGridIndex)%verticalCoordinate /= l_zeta ) &
!             & call Announce_Error ( key, No_Error_code, &
!             &  'Vertical coordinate in internal vGrid is not log pressure' )
          if ( .not. all(got( (/ f_temperatureQuantity, f_refGPHquantity /) ))) &
            & call Announce_Error ( key, no_error_code, &
            & 'wmoTropopause fill needs temperatureQuantity, refGPHQuantity' )
          if ( .not. ValidateVectorQuantity ( quantity, &
            & quantityType = (/l_boundaryPressure/), verticalCoordinate=(/l_none/), &
            & coherent=.true., stacked=.true. ) ) &
            & call Announce_Error ( key, no_error_code, &
            & 'Quantity is not a boundary pressure' )

          temperatureQuantity => GetVectorQtyByTemplateIndex( &
            &  vectors(temperatureVectorIndex), temperatureQuantityIndex)
          if ( .not. ValidateVectorQuantity ( temperatureQuantity, &
            & quantityType = (/l_temperature/), verticalCoordinate=(/l_zeta/), &
            & coherent=.true., stacked=.true. ) ) &
            & call Announce_Error ( key, badTemperatureQuantity )

          refGPHQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(refGPHVectorIndex), refGPHQuantityIndex)
          if ( .not. ValidateVectorQuantity ( refGPHQuantity, &
            & quantityType = (/l_refGPH/), verticalCoordinate=(/l_zeta/), &
            & coherent=.true., stacked=.true., noSurfs=(/1/) ) ) &
            & call Announce_Error ( key, badrefGPHQuantity )

          if ( .not. DoHGridsMatch ( quantity, temperatureQuantity ) .or. &
            &  .not. DoHGridsMatch ( quantity, refGPHQuantity ) ) &
            & call Announce_Error ( key, no_error_code, &
            & 'Horizontal coordinates for temperature/refGPH and/or quantity disagree' )

          ! OK, we must be ready to go
          if ( .not. USEREICHLER ) then
            call FillQtyWithWMOTropopause ( quantity, &
            & temperatureQuantity, refGPHQuantity, vGrids(internalVGridIndex) )
          else
            call FillQtyWithReichlerWMOTP ( quantity, &
            & temperatureQuantity, refGPHQuantity )
          end if
        case (-1)
          ! We must have decided to skip this fill
        case default
          call Announce_error ( key, no_Error_Code, 'This fill method not yet implemented' )
        end select      ! s_method

      case ( s_FillCovariance ) ! ===============  FillCovariance  =====
        invert = .false. ! Default if the field isn't present
        lengthScale = 0
        fraction = 0
        do j = 2, nsons(key)
          gson = subtree(j,key) ! The argument
          fieldIndex = get_field_id(gson)
          if ( nsons(gson) > 1) &
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
          if ( nsons(gson) > 1) &
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
        call destroyCommand ( key, matrices, vectors )

      case ( s_negativePrecision ) ! =======================  negativePrecision ==
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
            call expr_check ( gson , unitAsArray, valueAsArray, &
              & (/PHYQ_Dimensionless/), unitsError )
            if ( unitsError ) call Announce_error ( subtree(j,key), wrongUnits, &
              & extraInfo=(/unitAsArray(1), PHYQ_Dimensionless/) )
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

      case ( s_phase ) ! ===============================  Phase ==
        ! Set the name for this phase
        call addPhaseToPhaseNames ( vectorname, key )

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

!   User can put dumpBlocks, details=whatever, /allMatrices
!   and dump, details=whatever, /allVectors at the end of the Fill section.
!     if ( toggle(gen) ) then
!       if ( levels(gen) > 0 ) then
!         call dump ( vectors, details=levels(gen)-1 )
!         call dump ( matrices, details=levels(gen)-1 )
!       end if
!     end if
    if ( specialDumpFile /= ' ' ) call revertOutput
    if ( toggle(gen) ) call trace_end ( "MLSL2Fill" )
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
      if ( .not. FillableChiSq ( quantity, &
        & sourceQuantity, noiseQty ) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Incompatibility among vector quantities adding noise'  )
        return
      end if

     ! Either multiplier = [a, b] or multiplier = b are possible
      if ( .not. present(multiplier) ) then
        a = 1.
        b = 1.
      else if ( &
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

    ! ------------------------------------------- ApplyBaseline ----------
    subroutine ApplyBaseline ( key, quantity, baselineQuantity, &
      & quadrature, dontmask )
      integer, intent(in) :: KEY        ! Tree node
      type (VectorValue_T), intent(inout) :: QUANTITY ! Radiance quantity to modify
      type (VectorValue_T), intent(in) :: BASELINEQUANTITY ! L1B MAF baseline to use
      logical, intent(in) :: QUADRATURE ! If set add in quadrature (for noise)
      logical, intent(in) :: DONTMASK ! If set ignore baselinequantity mask
      ! Local variables
      integer :: MIF
      integer :: CHAN
      integer :: IND                    ! Combined MIF/CHAN
      integer :: i
      integer :: numProfs
      logical :: skipMe
      ! Executable code
      ! call output( 'Now in applyBaseline', advance='yes' )
      if ( quantity%template%quantityType /= l_radiance ) &
        & call Announce_Error ( key, no_Error_Code, &
        &   'Quantity to fill must be a radiance' )
      if ( baselineQuantity%template%quantityType /= l_l1bMAFBaseline ) &
        & call Announce_Error ( key, no_Error_Code, &
        &   'Quantity to fill must be a L1BMAFBaseline' )
      if ( baselineQuantity%template%signal /= quantity%template%signal .or. &
        &  baselineQuantity%template%sideband /= quantity%template%sideband ) &
        & call Announce_Error ( key, no_Error_Code, &
        &   'Quantity and baselineQuantity must have matching signal/sideband' )
      ind = 1
      numProfs = size( quantity%values ( ind, : ) )
      if ( quadrature ) then
        do mif = 1, quantity%template%noSurfs
          do chan = 1, quantity%template%noChans
            if ( .not. dontMask .and. associated(baselineQuantity%mask) ) then
              do i=1, numProfs
                skipMe = .not. dontMask .and. &
                  &  isVectorQtyMasked(baselineQuantity, chan, i)
                if ( .not. skipMe )  &
                & quantity%values ( ind, i ) = sqrt ( &
                  & quantity%values ( ind, i )**2 + &
                  & baselineQuantity%values ( chan, i )**2 )
              enddo
            else
              quantity%values ( ind, : ) = sqrt ( quantity%values ( ind, : )**2 + &
                & baselineQuantity%values ( chan, : )**2 )
            endif
            ind = ind + 1
          end do
        end do
      else
        do mif = 1, quantity%template%noSurfs
          do chan = 1, quantity%template%noChans
            if ( .not. dontMask .and. associated(baselineQuantity%mask) ) then
              do i=1, numProfs
                skipMe = .not. dontMask .and. &
                  &  isVectorQtyMasked(baselineQuantity, chan, i)
                if ( .not. skipMe )  &
                & quantity%values ( ind, i ) = &
                  & quantity%values ( ind, i ) + &
                  & baselineQuantity%values ( chan, i )
              enddo
            else
              quantity%values ( ind, : ) = quantity%values ( ind, : ) + &
                & baselineQuantity%values ( chan, : )
            endif
            ind = ind + 1
          end do
        end do
      end if
    end subroutine ApplyBaseline

    ! ------------------------------------- deallocateStuff ---
    subroutine deallocateStuff(Zetab, Zetac, Zetai, Pb, Pc, Pi)
      real (r8), pointer, dimension(:) :: Zetab
      real (r8), pointer, dimension(:) :: Zetac
      real (r8), pointer, dimension(:) :: Zetai
      real (r8), pointer, dimension(:) :: Pb         ! p[i] in hPa
      real (r8), pointer, dimension(:) :: Pc         ! p[i] in hPa
      real (r8), pointer, dimension(:) :: Pi         ! p[i] in hPa
      ! call Deallocate_test ( pa, 'pa', ModuleName )
      call Deallocate_test ( pb, 'pb', ModuleName )
      call Deallocate_test ( pc, 'pc', ModuleName )
      ! call Deallocate_test ( pd, 'pd', ModuleName )
      call Deallocate_test ( pi, 'pi', ModuleName )
      ! call Deallocate_test ( Zetaa, 'Zetaa', ModuleName )
      call Deallocate_test ( Zetab, 'Zetab', ModuleName )
      call Deallocate_test ( Zetac, 'Zetac', ModuleName )
      ! call Deallocate_test ( Zetad, 'Zetad', ModuleName )
      call Deallocate_test ( Zetai, 'Zetai', ModuleName )
    end subroutine DeallocateStuff

    ! ------------------------------------------- ExtractSingleChannel ---
    subroutine ExtractSingleChannel ( key, quantity, sourceQuantity, channel )
      integer, intent(in) :: KEY        ! Tree node
      type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      type (VectorValue_T), intent(in) :: SOURCEQUANTITY ! Source quantity for radiances
      integer, intent(in) :: CHANNEL    ! Channel number
      ! Local variables
      integer :: CHANIND                ! Channel index
      integer :: MIF                    ! Minor frame index
      ! Executable code
      if ( quantity%template%quantityType /= l_singleChannelRadiance ) &
        & call announce_error ( key, no_Error_Code, 'Quantity to fill must be of type singleChannelRadiance' )
      if ( all ( sourceQuantity%template%quantityType /= (/ l_cloudInducedRadiance, l_radiance /) ) ) &
        & call announce_error ( key, no_Error_Code, 'source quantity for fill must be of type [cloudInduced]radiance' )
      if ( quantity%template%signal /= sourceQuantity%template%signal .or. &
        & quantity%template%sideband /= sourceQuantity%template%sideband ) &
        & call announce_error ( key, no_Error_Code, 'quantity/sourceQuantity must be same signal/sideband' )
      if ( .not. sourceQuantity%template%regular ) &
        & call Announce_Error ( key, no_Error_Code, 'source quantity must be regular' )
      chanInd = channel - GetFirstChannel ( quantity%template%signal ) + 1
      do mif = 1, quantity%template%noSurfs
        quantity%values ( mif, : ) = &
          & sourceQuantity%values ( chanInd + ( mif - 1 ) * sourceQuantity%template%noChans, : )
      end do
    end subroutine ExtractSingleChannel

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
      real (r8) :: thisLength             ! Geometric mean length scale
      real (r8) :: meanDiag               ! Geometric mean diagonal value
      real (r8) :: thisFraction           ! Geometric mean diagonal value
      logical, dimension(:), pointer :: condition ! Condition
      logical :: ANYOFFDIAG             ! Flag to indicate presence of off diagonal elements

      ! Executable code

      ! Apply mask to diagonal
      nullify ( m, condition )
      call CopyVector ( Dmasked, vectors(diagonal), clone=.true., &
        & vectorNameText='_Dmasked' )
      call ClearUnderMask ( Dmasked )

      if ( lengthScale == 0 ) then
        call updateDiagonal ( covariance, vectors(diagonal), square=.true.,&
          & invert=invert, forgiveZeros=.true. )
      else
        ! Do a more complex fill, either we're doing non-diagonal, or there might
        ! be zeros to 'invert'

        ! Setup some stuff
        if ( lengthScale /= 0 ) then
          call CopyVector ( Lmasked, vectors(lengthScale), clone=.true., &
            & vectorNameText='_Lmasked' )
          call ClearUnderMask ( Lmasked )
        end if

        ! Check the validity of the supplied vectors
        if ( covariance%m%row%vec%template%name /= &
          & dMasked%template%name ) call MLSMessage ( MLSMSG_Error, &
          & ModuleName, "diagonal and covariance not compatible in fillCovariance" )
        if ( lengthScale /= 0 ) then    ! Check length if supplied
          if ( covariance%m%row%vec%template%name /= &
            & lMasked%template%name ) call MLSMessage ( MLSMSG_Error, &
            & ModuleName, "lengthScale and covariance not compatible in fillCovariance" )
          if ( lMasked%globalUnit /= phyq_length ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "length vector does not have dimensions of length" )
        else
          thisLength = 0.0
        end if
        if ( fraction /= 0 ) then       ! Check fraction if supplied
          if ( covariance%m%row%vec%template%name /= &
            & vectors(fraction)%template%name ) call MLSMessage ( MLSMSG_Error, &
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
          if ( lengthScale /= 0 ) l => lMasked%quantities(q)
          if ( fraction /=0 ) f => vectors(fraction)%quantities(q)
          n = qt%instanceLen
          if ( qt%coherent ) surfs => qt%surfs(:,1)
          if ( .not. qt%regular ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "Unable to handle irregular quantity in FillCovariance" )
          call Allocate_test ( m, n, n, 'M', ModuleName )

          ! Loop over the instances
          do i = 1, qt%noInstances
            anyOffDiag = .false.
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
                do k = 1, j-1
                  meanDiag = sqrt ( m(j,j) * m(k,k) )
                  if ( lengthScale /= 0 ) &
                    & thisLength = sqrt ( l%values(j,i) * l%values(k,i) )
                  if ( fraction /= 0 ) thisFraction = f%values(j,i)
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
                  if ( thisLength > 0.0 .and. thisFraction > 0.0 ) then
                    m(j,k) = meanDiag*thisFraction*exp(-distance/thisLength)
                    anyOffDiag = .true.
                  end if
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
              if ( anyOffDiag ) then
                call MatrixInversion(M, upper=.true.)
              else
                do j = 1, n
                  m(j,j) = 1.0 / m(j,j)
                end do
              end if
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
    subroutine FillVectorQuantityFromGrid(quantity, grid, allowMissing, errorCode)
      ! Dummy arguments
      type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      type (GriddedData_T), intent(inout) :: GRID ! Grid to fill it from
      ! Needs to be inout because we wrap it.
      logical, intent(in) :: ALLOWMISSING ! If set missing data in grid ok
      integer, intent(out) :: ERRORCODE   ! Error code (one of constants defined above)

      ! Local variables
      integer :: instance,surf            ! Loop counter
      integer :: instIndex,surfIndex      ! Indices
      real(rv) :: newValue

      ! Executable code
      errorCode = 0
      if ( grid%empty ) then
        if ( index(lowercase(grid%description), 'tropopause') > 0 ) then
          ! Must allow this as missing gmao files are a possibility
          ! to be handled with grace and aplomb
          call MLSMessage ( MLSMSG_Warning, moduleName, &
            & 'No tropopause values in grid--filling with missing values' )
          quantity%values = grid%missingValue
        else
          errorCode=EmptyGridForFill
        end if
        return
      end if

      if ( quantity%template%verticalCoordinate /= l_zeta .and. &
        & quantity%template%noInstances > 1 .and. grid%noHeights > 1 ) then
        errorCode=NotZetaForGrid
        return
      end if

      instIndex=1
      surfIndex=1

      ! Wrap the grid to be sure that we can interpolate it in longitude
      ! This will skip out if it's already been done.
      call WrapGriddedData ( grid )

      do instance = 1, quantity%template%noInstances
        if ( .not. quantity%template%stacked) instIndex=instance

        do surf = 1, quantity%template%noSurfs
          if ( .not. quantity%template%coherent) surfIndex=surf
          call l3ascii_interp_field(grid, newValue, &
            & pressure=10.0**(-quantity%template%surfs(surf,instIndex)), &
            & lat=quantity%template%geodLat(surfIndex,instance), &
            & lon=quantity%template%lon(surfIndex,instance), &
            & lst=quantity%template%solarTime(surfIndex,instance), &
            & sza=quantity%template%solarZenith(surfIndex,instance), &
            & date=quantity%template%time(surfIndex,instance))
          if ( newValue >= nearest ( grid%missingValue, -1.0 ) .and. &
            &  newValue <= nearest ( grid%missingValue,  1.0 ) .and. &
            & .not. allowMissing ) errorCode = MissingDataInGrid
          quantity%values(surf,instance) = newValue
        end do                            ! End surface loop
      end do                              ! End instance loop
    end subroutine FillVectorQuantityFromGrid

    !=============================== FillVectorQuantityFromL2GP ==========
    subroutine FillVectorQuantityFromL2GP ( quantity,l2gp, interpolate, profile, &
      & errorCode, ignoreGeolocation, fromPrecision )
      use MLSNumerics, only: COEFFICIENTS_R8, INTERPOLATEARRAYSETUP, &
        & INTERPOLATEARRAYTEARDOWN

      ! If the times, pressures, and geolocations match, fill the quantity with
      ! the appropriate subset of profiles from the l2gp

      ! Dummy arguments
      type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      type (L2GPData_T), intent(in), target :: L2GP ! L2GP to fill from
      logical, intent(in) :: interpolate  ! Flag
      integer, intent(in) :: profile    ! Single profile to use or -1 for default
      integer, intent(out) :: errorCode ! Error code
      logical, intent(in) :: ignoreGeolocation  ! Flag
      logical, intent(in) :: fromPrecision ! Flag

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
      real (r4), dimension(:,:,:), pointer :: SOURCE

      errorCode=0
      ! Make sure this quantity is appropriate
      if ( .not. ValidateVectorQuantity(quantity, coherent=.TRUE., stacked=.TRUE., &
        & verticalCoordinate= (/ l_pressure, l_zeta /) ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Quantity to fill is not on pressure or zeta coordinates' )

      if ( (quantity%template%noChans/=l2gp%nFreqs) .and. &
        &  ((quantity%template%noChans/=1) .or. (l2gp%nFreqs/=0)) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Quantity and l2gp have different number of channels' )
!       call dump ( l2gp%frequency, 'l2gp' )
!       call dump ( quantity%template%frequencies, 'quantity' )
!       if ( associated ( quantity%template%frequencies ) ) then
!         if ( any ( abs ( l2gp%frequency - &
!           & quantity%template%frequencies ) > fTol ) ) &
!           & call MLSMessage ( MLSMSG_Error, ModuleName, &
!           & 'Quantity and l2gp have different frequency grids' )
!       end if

      if ( quantity%template%noSurfs /= l2gp%nLevels .and. (.not. interpolate) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Quantity and l2gp have different number of surfaces (set interpolate?)' )

      if ( .not. interpolate ) then
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
      firstProfile = 1
      lastProfile=firstProfile+quantity%template%noInstances-1
      if ( (profile == -1) .and. (.not. ignoreGeolocation) ) then
        ! Attempt to match up the first location
        firstProfileAsArray=minloc(abs(quantity%template%phi(1,1)-l2gp%geodAngle))
        firstProfile=firstProfileAsArray(1)

        ! Well, the last profile has to be noInstances later, check this would be OK
        lastProfile=firstProfile+quantity%template%noInstances-1
        if ( lastProfile > l2gp%nTimes ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Quantity has profiles beyond the end of the l2gp' )

        ! Now check that geodAngle's are a sufficient match
        if ( any(abs(l2gp%geodAngle(firstProfile:lastProfile)-&
          &         quantity%template%phi(1,:)) > tolerance) ) then
          call dump ( l2gp%geodAngle(firstProfile:lastProfile), 'L2GP geodetic angle' )
          call dump ( quantity%template%phi(1,:), 'Quantity Geodetic angle' )
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Quantity has profiles that mismatch l2gp in geodetic angle' )
        end if

        ! Now check that the times match
        if ( any(abs(l2gp%time(firstProfile:lastProfile)- &
          &         quantity%template%time(1,:)) > timeTol) ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Quantity has profiles that mismatch l2gp in time' )

        ! Currently the code cannot interpolate in 3 dimensions, wouldn't
        ! be hard to code up, but no need as yet.
        if ( interpolate .and. quantity%template%noChans /= 1 ) then
          errorCode=cantInterpolate3D
          return
        end if
      else
        ! Given a specific profile, check it's legal
        if ( profile == 0 .or. profile > l2gp%nTimes ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Illegal profile request in l2gp fill' )
      end if

      ! OK, things seem to be OK, so start getting the data
      if ( fromPrecision ) then
        source => l2gp%l2gpPrecision
      else
        source => l2gp%l2gpValue
      end if

      ! OK, now do the filling, it's easier if we don't have to interpolate
      if ( .not. interpolate ) then
        if ( profile == -1 ) then
          ! Not forcing a particular profile to all instances
          quantity%values=reshape(source(:,:,firstProfile:lastProfile),&
            & (/quantity%template%noChans*quantity%template%noSurfs,&
            &   quantity%template%noInstances/))
        else
          ! Spread one profile onto all instances
          quantity%values = spread ( reshape ( source(:,:,profile), &
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
            & real(source(1,:,thisProfile), r8), & ! OldY
            & outZeta, & ! New X
            & quantity%values(:,instance), & ! New Y
            & method='Linear', extrapolate='Clamp' )
        end do
        call InterpolateArrayTeardown ( coeffs )
      end if

    end subroutine FillVectorQuantityFromL2GP

    ! -------------------------------------- FillVectorQuantityFromProfile --
    subroutine FillVectorQtyFromProfile ( quantity, valuesNode, &
      & instancesNode, globalUnit, dontMask, logSpace )
      use MLSNumerics, only: HUNT
      ! This fill is slightly complicated.  Given a values array along
      ! the lines of [ 1000mb : 1.0K, 100mb : 1.0K,  10mb : 2.0K] etc. it
      ! does the linear interpolation appropriate to perform the fill.
      type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      integer, intent(in) :: VALUESNODE   ! Tree node for values
      integer, intent(in) :: INSTANCESNODE ! Tree node for instances
      integer, intent(in) :: GLOBALUNIT   ! Possible global unit
      logical, intent(in) :: DONTMASK     ! If set don't follow the fill mask
      logical, intent(in), optional :: LOGSPACE ! Interpolate in logspace

      ! Local variables
      integer :: C                      ! Channel loop counter
      integer :: HEIGHTUNIT             ! Unit for height
      integer :: NOPOINTS               ! Number of points supplied
      integer :: NOUNIQUE               ! Number of unique heights supplied
      integer :: I,J                    ! Loop counters / indices
      integer :: S                      ! Surface loop counter
      integer :: STATUS                 ! Flag
      integer :: TESTUNIT               ! Unit for value
      logical :: Fail                   ! Status from Hunt
      logical :: LOCALOUTHEIGHTS ! Set if out heights is our own variable
      logical :: MYLOGSPACE             ! Interpolate in log space?
      real (r8), dimension(:), pointer :: HEIGHTS ! Heights for the points
      real (r8), dimension(:), pointer :: VALUES ! Values for the points
      real (r8), dimension(:), pointer :: OUTHEIGHTS ! Heights for output
      real (r8), dimension(:), pointer :: OUTVALUES ! Single profile for output
      logical, dimension(:), pointer :: INSTANCES ! Flags
      logical, dimension(:), pointer :: DUPLICATED ! Flags
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
      myLogSpace = quantity%template%logBasis
      if ( present ( logSpace ) ) myLogSpace = logSpace
      noPoints = nsons ( valuesNode ) - 1
      nullify ( heights, values, duplicated, outHeights, outValues, instances )
      call Allocate_test ( heights, noPoints, 'heights', ModuleName )
      call Allocate_test ( values, noPoints, 'values', ModuleName )
      call Allocate_test ( duplicated, noPoints, 'duplicated', ModuleName )
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
          & call Announce_error ( valuesNode, no_Error_Code, 'Bad height units for profile fill' )
        ! Store height
        if ( heightUnit == PHYQ_Zeta ) then
          ! Assume zeta coordinates are expressed in mb
          heights(i) = -log10 ( exprValue(1) )
        else
          heights(i) = exprValue(1)
        end if
        ! Check value unit OK
        if ( all ( exprUnit(2) /= (/ testUnit, PHYQ_Dimensionless /) ) ) &
          & call Announce_error ( valuesNode, no_error_code, 'Bad units for profile fill' )
        ! Store value
        values ( i ) = exprValue(2)
      end do

      if ( myLogSpace .and. any ( values <= 0.0 ) ) then
        call Announce_Error ( valuesNode, no_error_code, &
          & 'Non-positive input data in log profile fill' )
        return
      end if
      if ( myLogSpace ) values = log ( values )

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
          & nearest=.true., allowTopValue=.true., fail=fail )
        if ( fail ) then
          call Announce_Error ( valuesNode, no_error_code, &
          & 'Problem in Hunt' )
        return
      end if
        duplicated = .false.
        do i = 1, noPoints - 1
          do j = i + 1, noPoints
            if ( inInds(i) == inInds(j) ) then
              duplicated ( j ) = .true.
            end if
          end do
        end do
        noUnique = count ( .not. duplicated )
        inInds(1:noUnique) = pack ( inInds, .not. duplicated )
        heights(1:noUnique) = outHeights ( inInds(1:noUnique) )
        values(1:noUnique) = pack ( values, .not. duplicated )
        call deallocate_test ( inInds, 'inInds', ModuleName )
      end if

      ! Now do the interpolation for the first instance
      call InterpolateValues ( heights(1:noUnique), values(1:noUnique), outHeights, &
        & outValues, 'Linear', extrapolate='Constant' )

      if ( myLogSpace ) outValues = exp ( outValues )

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
              if ( associated(quantity%mask) .and. .not. dontMask ) then
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
      call Deallocate_test ( duplicated, 'duplicated', ModuleName )
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
          & call announce_error( no_Error_Code, No_Error_code, &
          & 'Noise values unassociated in FillableChiSq')
      end if
      aok = aok .and. &
        & associated(qty%values) .and. &
        & associated(measQty%values) .and. &
        & associated(modelQty%values)

      if ( DEEBUG ) then
        if ( .not. associated(qty%values) ) &
          & call announce_error( no_Error_Code, No_Error_code, &
          & 'Quantity values unassociated in FillableChiSq')
        if ( .not. associated(measQty%values) ) &
          & call announce_error( no_Error_Code, No_Error_code, &
          & 'Measurements values unassociated in FillableChiSq')
        if ( .not. associated(modelQty%values) ) &
          & call announce_error( no_Error_Code, No_Error_code, &
          & 'Model values unassociated in FillableChiSq')
      end if

      if ( .not. aok ) return

      minorFrame = qty%template%minorFrame .or. qty%template%majorFrame
      ! (1)
      if ( .not. minorFrame ) then
        aok = aok .and. &
          & (qty%template%molecule == measQty%template%molecule) &
          & .and. &
          & (qty%template%molecule == modelQty%template%molecule)
        if ( present(noiseQty) ) aok = aok &
          & .and. &
          & (qty%template%molecule == noiseQty%template%molecule)
      end if

      ! (2)
      if ( minorFrame ) then
        aok = aok .and. &
          & (qty%template%signal == measQty%template%signal) &
          & .and. &
          & (qty%template%signal == modelQty%template%signal)
        if ( present(noiseQty) ) aok = aok &
          & .and. &
          & (qty%template%signal == noiseQty%template%signal)
      end if

      ! (3)
      if ( .not. minorFrame ) then
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
      else if ( &
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
      if ( .not. ValidateVectorQuantity ( qty, &
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
      else if ( .not. FillableChiSq ( qty, measQty, modelQty, noiseQty ) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Incompatibility among vector quantities filling chi^2 channelwise'  )
        return
      else if ( any ( noiseQty%values == 0.0) .and. &
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
      else if ( &
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
      if ( .not. ValidateVectorQuantity ( qty, &
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
      else if ( .not. FillableChiSq ( qty, measQty, modelQty, noiseQty ) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Incompatibility among vector quantities filling chi^2 MMAFwise'  )
        return
      else if ( any ( noiseQty%values == 0.0) .and. &
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
        if ( .not. (.not. dontMask .or. ignoreNegative .or. ignoreZero ) ) then
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
      else if ( &
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
      if ( .not. ValidateVectorQuantity ( qty, &
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
      else if ( .not. FillableChiSq ( qty, measQty, modelQty, noiseQty ) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Incompatibility among vector quantities filling chi^2 MMIFwise'  )
        return
      else if ( any ( noiseQty%values == 0.0) .and. &
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

    ! ------------------------------------------- FillChiSqRatio ---
    subroutine FillChiSqRatio ( key, qty, normQty, minNormQty, flagQty, &
    & dontMask, firstInstance, lastInstance )
      ! A special fill of the ratio
      !  chi squared Norm
      ! ----------------     [iter_n, *]
      ! chi squared Min Norm
      ! where iter_n is the final iteration number
      
      ! Note the following tricks:
      ! The number of surfaces is the maximum allowed number of iterations
      ! The actual number of iterations will be less than this
      ! After the last iteration, all "surfaces" above this are zero-filled
      
      ! The number of instances will be the number of chunks
      ! (yes, an unfortunate fact)
      integer, intent(in) :: KEY
      type (VectorValue_T), intent(inout) :: QTY
      type (VectorValue_T), intent(inout) ::    normQty
      type (VectorValue_T), intent(inout) ::    minNormQty
      type (VectorValue_T), intent(inout) ::    flagQty
      logical, intent(in)           ::       dontMask    ! Use even masked values

      integer, intent(in), optional ::       firstInstance, lastInstance
      ! The last two are set if only part (e.g. overlap regions) of the quantity
      ! is to be stored in qty

      ! Local variables
      integer ::                             UseFirstInstance, UseLastInstance, &
      &                                      NoOutputInstances
      integer ::                             I           ! Instances
      integer ::                             ITER        ! Instances
      integer ::                             QINDEX
      integer ::                             NOCHANS
      logical ::                             skipMe
      logical, parameter ::                  FakeData = .false.

      ! Executable code

      ! First check that things are OK.
      if ( .not. ValidateVectorQuantity ( qty, &
        & quantityType=(/l_dnwt_chiSqRatio/) ) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Attempting to fill wrong quantity with chi^2 ratio'  )
        return
      elseif ( .not. ValidateVectorQuantity ( normqty, &
        & quantityType=(/l_dnwt_chiSqNorm/) ) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Attempting to fill using wrong norm quantity with chi^2 ratio'  )
        return
      elseif ( .not. ValidateVectorQuantity ( minnormqty, &
        & quantityType=(/l_dnwt_chiSqMinNorm/) ) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Attempting to fill using wrong min norm quantity with chi^2 ratio'  )
        return
      elseif ( .not. ValidateVectorQuantity ( flagqty, &
        & quantityType=(/l_dnwt_flag/) ) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Attempting to fill using wrong flag quantity with chi^2 ratio'  )
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

      noChans = qty%template%noChans
      do i=useFirstInstance, useLastInstance
        if ( FakeData ) then
          ! Let's just fake up some data here
          flagQty%values(:,i) = 0._rv
          normQty%values(:,i) = 0._rv
          minNormQty%values(:,i) = 0._rv
          ! Say we converged at iteration number 5
          qIndex = 5
          flagQty%values(:qIndex,i) = -2._rv
          minNormQty%values(:qIndex,i) = -2._rv
          do iter = 1, qIndex
            minNormQty%values(iter, i) = 1.1
            normQty%values(iter, i) = 1.1 + 0.05*(qIndex+1-iter)
          enddo
        endif
        ! Now find the iteration number
        qIndex = findLast( flagQty%values(:,i) /= 0._rv )
        if ( qIndex == 0 .or. qIndex >= qty%template%noSurfs ) cycle
        skipMe = &
          & .not. dontMask .and. ( &
          &   isVectorQtyMasked(normQty, qIndex, i) .or. &
          &   isVectorQtyMasked(minNormQty, qIndex, i) .or. &
          &   minNormQty%values(qIndex, i) == 0. &
          & )
          qty%values(:,i) = &
            & normQty%values(qIndex, i) / minNormQty%values(qIndex, i)
      end do
    end subroutine FillChiSqRatio

    ! ------------------------------------------- FillColAbundance ---
    subroutine FillColAbundance ( key, qty, bndPressQty, vmrQty, colmAbUnits, &
      & firstInstance, lastInstance )
      ! A special fill according to W.R.Read's idl code
      ! Similar to his hand-written notes, but with a small correction

      ! Assumptions:
      ! (See above)
      integer, intent(in) :: KEY
      type (VectorValue_T), intent(inout) :: QTY
      type (VectorValue_T), intent(in) :: BNDPRESSQTY
      type (VectorValue_T), intent(in) :: VMRQTY
      integer, intent(in) :: colmAbUnits
      integer, intent(in), optional :: FIRSTINSTANCE
      integer, intent(in), optional :: LASTINSTANCE
      ! The last two are set if only part (e.g. overlap regions) of the quantity
      ! is to be stored in the column data.

      ! Local parameters
      real(r8) :: LN10         ! = LOG(10.d0)
      real(r8), parameter :: INVERMGPPMV = 0.789 ! in DU / (ppmv hPa)
      real(r8), parameter :: INVERMGDU = INVERMGPPMV*1.d6 ! in DU / (vmr hPa)
      real(r8), parameter :: INVERMGMOLCM2 = INVERMGDU*2.687d16 ! in mol / cm^2

      ! Local variables
      logical :: printMe
      integer :: SURFACE
      integer :: INSTANCE
      integer :: FIRSTSURFACEABOVE
      integer :: FIRSTSURFACEBELOW
      integer :: USEFIRSTINSTANCE
      integer :: USELASTINSTANCE
      integer :: N
      integer :: NOOUTPUTINSTANCES
      real (r8) :: THISBNDPRESS
      real (r8) :: zeta            ! -log(THISBNDPRESS)
      real (r8) :: zetaTmp
      real (r8) :: COLUMNSUM
      real (r8) :: TRAPEZOIDSUM
      real (r8) :: DELTAZETA      ! Zetai[s+1] - Zetai[s]
      real (r8) :: INVERMG
      real (r8)                        :: Zetaa
      real (r8), pointer, dimension(:) :: Zetab
      real (r8), pointer, dimension(:) :: Zetac
      real (r8)                        :: Zetad
      real (r8), pointer, dimension(:) :: Zetai
      real (r8)                        :: Pa         ! p[i] in hPa
      real (r8), pointer, dimension(:) :: Pb         ! p[i] in hPa
      real (r8), pointer, dimension(:) :: Pc         ! p[i] in hPa
      real (r8)                        :: Pd         ! p[i] in hPa
      real (r8), pointer, dimension(:) :: Pi         ! p[i] in hPa

      ! Executable code
      LN10 = LOG(10.d0)  ! Compiler won't let this be a parameter
      ! First check that things are OK.
      if ( (qty%template%quantityType /= l_columnAbundance) .or. &
        &  (bndPressQty%template%quantityType /= l_boundaryPressure) .or. &
        &  (vmrQty%template%quantityType /= l_vmr) ) then
        call Announce_error ( key, No_Error_code, &
          & 'Wrong quantity type found while filling column abundance'  )
        return
      else if ( qty%template%molecule /= vmrQty%template%molecule ) then
        call Announce_error ( key, No_Error_code, &
          & 'Attempt to fill column abundance with different molecule'  )
        return
      else if ( .not. ( DoHgridsMatch( qty, vmrQty ) .and. &
        & DoHgridsMatch( qty, bndPressQty ) ) ) then
        call Announce_error ( key, No_Error_code, &
          & 'Attempt to fill column abundance with different HGrids'  )
        return
      else if ( .not. any(vmrQty%template%verticalCoordinate == &
        & (/l_zeta/)) ) then
        call Announce_error ( key, No_Error_code, &
          & 'Fill column abundance, but vmr not on zeta surfs.'  )
        return
      else if ( vmrQty%template%noSurfs < 2 ) then
        call Announce_error ( key, No_Error_code, &
          & 'Fill column abundance, but too few vmr surfaces'  )
        return
      end if
      select case (colmAbUnits)
      case (l_DobsonUnits)
        INVERMG = INVERMGDU
      case (l_DU)
        INVERMG = INVERMGDU
      case (l_molcm2)
        INVERMG = INVERMGMOLCM2
      case default
        call Announce_error ( key, No_Error_code, &
          & 'Fill column abundance, but wrong units.'  )
      end select
      ! Work out what to do with the first and last Instance information
      useFirstInstance = 1
      useLastInstance = qty%template%noInstances
      if ( present ( firstInstance ) ) useFirstInstance = firstInstance
      if ( present ( lastInstance ) ) useLastInstance = lastInstance
      noOutputInstances = useLastInstance-useFirstInstance+1

      ! If we've not been asked to output anything then don't carry on
      if ( noOutputInstances < 1 ) return

      nullify ( Zetab, Zetac, Zetai, pb, pc, pi )
      call allocate_test ( Zetab, vmrQty%template%noSurfs, 'Zetab', ModuleName )
      call allocate_test ( Zetac, vmrQty%template%noSurfs, 'Zetac', ModuleName )
      call allocate_test ( Zetai, vmrQty%template%noSurfs, 'Zetai', ModuleName )
      call allocate_test ( pb, vmrQty%template%noSurfs, 'pb', ModuleName )
      call allocate_test ( pc, vmrQty%template%noSurfs, 'pc', ModuleName )
      call allocate_test ( pi, vmrQty%template%noSurfs, 'pi', ModuleName )
      Zetai = vmrQty%template%surfs(:,1)
      pi = 10.0 ** ( -Zetai  )
      N = vmrQty%template%noSurfs
      do instance = useFirstInstance, useLastInstance
        printMe = ( instance == useFirstInstance) .and. &
         & (switchDetail(switches, 'column') > -1 )
        if ( printMe ) print *, 'switches: ', trim(switches)
        if ( printMe ) &
          & print *, 'switchDetail(switches, column) ', switchDetail(switches, 'column')
        ! Find 1st surface at or above tropopause
        ! (i.e., at a pressure equal to or less than boundaryPressure)
        ! This next check should be unnecessary--the HGrids were already matched
        if ( instance > size(bndPressQty%values, 2) ) then
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Cant fill column--instance outside b.pres. range' )
          call deallocateStuff(Zetab, Zetac, Zetai, Pb, Pc, Pi)
          return
        end if
        thisBndPress = bndPressQty%values(1,instance)
        ! In case where WMO algorithm failed, use bottom of basis
        if ( thisBndPress <= 0.0 ) &
          & thisBndPress = 10.0 ** ( - vmrQty%template%surfs(1,1) )
        if ( thisBndPress <= 0._r8 ) then
          call Announce_error ( key, No_Error_code, &
          & 'Fill column abundance, illegal bound. pr. at this instance' )
        end if
        zeta = -log10 ( thisBndPress )
        Zetaa = max(Zetai(N), zeta)
        do surface=1, N
          Zetab(surface) = Zetai(min(surface+1,N))
        end do
        do surface=1, N
          Zetab(surface) = min(max(zeta, Zetai(surface)), Zetab(surface))
          zetaTmp = zeta
          if ( surface > 1 ) zetaTmp = max(zeta, Zetai(surface-1))
          Zetac(surface) = min(zetaTmp, Zetai(surface))
        end do
        Zetad = min(Zetai(1), zeta)
        Pa = 10. ** (-Zetaa)
        Pb = 10. ** (-Zetab)
        Pc = 10. ** (-Zetac)
        Pd = 10. ** (-Zetad)
        if ( printMe ) then
          call output ( 'zeta: ', advance='no')
          call output ( zeta, advance='yes')
          call output ( 'zetaA: ', advance='no')
          call output ( zetaA, advance='yes')
          call output ( 'zetaD: ', advance='no')
          call output ( zetaD, advance='yes')
          call output ( 'PA: ', advance='no')
          call output ( PA, advance='yes')
          call output ( 'PD: ', advance='no')
          call output ( PD, advance='yes')
          call output ( 'zetaI   zetaB     zetaC      Pi   Pb   Pc', advance='yes')
          do surface=1, N
            call output( (/ zetai(surface), zetab(surface), zetac(surface), &
             & Pi(surface), Pb(surface), Pc(surface)/) , advance='yes')
          end do
        end if
        ! Find 1st surface immediately above tropopause
        firstSurfaceAbove = FindFirst (Pi < thisBndPress)
        if ( firstSurfaceAbove < 1 ) then
          call output ( 'tropopause: ', advance='no')
          call output ( thisBndPress, advance='yes')
          call output ( 'p0 ', advance='no')
          call output ( Pi(1), advance='yes')
          call output ( 'pTop ', advance='no')
          call output ( Pi(N), advance='yes')
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Filling column, but tropopause outside pres. surfaces' )
          firstSurfaceBelow = 1
        else if ( firstSurfaceAbove == 1 ) then
          ! Nothing special
          firstSurfaceBelow = 1
        else if ( firstSurfaceAbove > N ) then
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Cant column, tropopause above top surface' )
          call deallocateStuff(Zetab, Zetac, Zetai, Pb, Pc, Pi)
          return
        else  ! Nothing special
          firstSurfaceBelow = firstSurfaceAbove - 1
        end if
        if ( printMe ) then
          print *, 'thisBndPress: ', thisBndPress
          print *, 'firstSurfaceAbove: ', firstSurfaceAbove
          print *, 'firstSurfaceBelow: ', firstSurfaceBelow
        end if
        ! Do summation
        columnSum = vmrQty%values(N, instance) * Pa  ! Initialize sum
        if ( printMe ) then
          print *, 'columnSum: ', columnSum
          trapezoidSum = 0.
          do surface=firstSurfaceAbove, N-1
            trapezoidsum = trapezoidsum + 0.5 * ( &
              & vmrQty%values(surface, instance) + vmrQty%values(surface+1, instance) &
              & ) * ( Pi(surface+1)-Pi(surface))
          end do
        end if
        ! Loop over surfaces from 1 to uppermost-1
        do surface = 1, N-1
          deltaZeta = Zetai(surface+1) - Zetai(surface)
          if ( deltaZeta == 0._r8 ) cycle
          columnSum = columnSum + &
            & vmrQty%values(surface, instance) / deltaZeta * &
            & ( &
            & Pb(surface)*(Zetai(surface+1)-Zetab(surface)) &
            & + &
            & (Pi(surface+1) - Pb(surface))/ln10 &
            & )
        if ( printMe ) print *, 'columnSum: ', columnSum, surface
        end do
        ! Loop over surfaces from 2 to uppermost
        do surface = 2, N
          deltaZeta = Zetai(surface) - Zetai(surface-1)
          if ( deltaZeta == 0._r8 ) cycle
          columnSum = columnSum + &
            & vmrQty%values(surface, instance) / deltaZeta * &
            & ( &
            & Pc(surface)*(Zetac(surface)-Zetai(surface)) &
            & + &
            & (Pc(surface) - Pi(surface))*(1./ln10 + deltaZeta) &
            & )
          ! columnSum = columnSum + &
          ! & vmrQty%values(surface, instance) / deltaZeta * &
          ! & ( &
          ! & Pc(surface)*(Zetac(surface)-Zetai(surface)) &
          ! & + &
          ! & (Pc(surface) - Pi(surface))/ln10 &
          ! & )
        if ( printMe ) print *, 'columnSum: ', columnSum, surface
        end do
        columnSum = columnSum + vmrQty%values(1, instance) * (Pd - Pi(1))
        if ( printMe ) print *, 'columnSum: ', columnSum
        if ( printMe ) print *, 'trapezoid: ', -trapezoidSum
        qty%values ( 1, instance ) = InverMg * columnSum
      end do
      call deallocateStuff(Zetab, Zetac, Zetai, Pb, Pc, Pi)
    end subroutine FillColAbundance

    ! ------------------------------------- FillFoldedRadiance ---
    subroutine FillFoldedRadiance ( radiance, lsb, usb, &
      & lsbFraction, usbFraction, key )
      type (VectorValue_T), intent(inout) :: RADIANCE
      type (VectorValue_T), pointer :: USB
      type (VectorValue_T), pointer :: LSB
      type (VectorValue_T), pointer :: USBFRACTION
      type (VectorValue_T), pointer :: LSBFRACTION
      integer, intent(in) :: KEY

      ! Local variables
      integer :: C                        ! Channel loop inductor
      integer :: I                        ! Array index
      integer :: MIF                      ! Minor frame loop inductor

      ! Executable code
      ! First some sanity checks
      if ( .not. ValidateVectorQuantity ( radiance, quantityType=(/l_radiance/), &
        & sideband=(/0/), minorFrame=.true. )) &
        & call Announce_Error ( key, no_error_code, 'Inappropriate radiance quantity to fill' )
      if ( ( associated ( lsb ) .neqv. associated ( lsbFraction ) ) .or. &
        &  ( associated ( usb ) .neqv. associated ( usbFraction ) ) ) then
        call Announce_Error ( key, no_error_code, 'Must supply sidebands and fractions together' )
        return
      end if
      if ( associated ( lsb ) ) then
        if ( .not. ValidateVectorQuantity ( lsb, quantityType=(/l_radiance/), &
          & sideband=(/-1/), signal=(/radiance%template%signal/), minorFrame=.true. )) &
          & call Announce_Error ( key, no_error_code, 'Inappropriate lsb radiance quantity for fill' )
        if ( .not. ValidateVectorQuantity ( lsbFraction, &
          & quantityType=(/l_limbSidebandFraction/), &
          & signal=(/radiance%template%signal/), sideband=(/-1/) ) ) &
          & call Announce_Error ( key, no_error_code, 'Inappropriate lsbFraction quantity for fill' )
      end if
      if ( associated ( usb ) ) then
        if ( .not. ValidateVectorQuantity ( usb, quantityType=(/l_radiance/), &
          & sideband=(/1/), signal=(/radiance%template%signal/), minorFrame=.true. )) &
          & call Announce_Error ( key, no_error_code, 'Inappropriate usb radiance quantity for fill' )
        if ( .not. ValidateVectorQuantity ( usbFraction, &
          & quantityType=(/l_limbSidebandFraction/), &
          & signal=(/radiance%template%signal/), sideband=(/1/) ) ) &
          & call Announce_Error ( key, no_error_code, 'Inappropriate usbFraction quantity for fill' )
      end if

      if ( .not. associated ( lsb ) .and. .not. associated ( usb ) ) &
        & call Announce_Error ( key, no_error_code, 'Must supply one or both of lsb/usb' )

      ! Now do the work
      ! Note that this is a very inefficient loop, but I don't really care as it's
      ! Never used in routine processing.
      i = 1                               ! Use i as a composit mif,channel index
      do mif = 1, radiance%template%noSurfs
        do c = 1, radiance%template%noChans
          radiance%values(i,:) = 0.0
          if ( associated ( lsb ) ) radiance%values(i,:) = radiance%values(i,:) + &
            & lsbFraction%values(c,1) * lsb%values(i,:)
          if ( associated ( usb ) ) radiance%values(i,:) = radiance%values(i,:) + &
            & usbFraction%values(c,1) * usb%values(i,:)
          i = i + 1
        end do
      end do

    end subroutine FillFoldedRadiance

    ! ------------------------------------ FillPhiTanWithRefraction --
    subroutine FillPhiTanWithRefraction ( key, quantity, &
      & h2o, orbIncline, ptan, refGPH, temperature )

      use Geometry, only: EarthRadA, EarthRadB, GEODTOGEOCLAT
      use Hydrostatic_M, only: HYDROSTATIC
      use MLSKinds, only: RP
      use MLSNumerics, only: InterpolateValues
      use Phi_Refractive_Correction_m, only: Phi_Refractive_Correction_Up
      use Refraction_m, only: REFRACTIVE_INDEX
      use Units, only: DEG2RAD, RAD2DEG

      integer, intent(in) :: KEY          ! Tree node, for error messages
      type (VectorValue_T), intent(inout) :: QUANTITY ! PhiTan quantity to update
      type (VectorValue_T), intent(in) :: H2O         ! Water vapor
      type (VectorValue_T), intent(in) :: OrbIncline  ! Orbital inclination
      type (VectorValue_T), intent(in) :: PTan        ! Tangent pressure
      type (VectorValue_T), intent(in) :: RefGPH      ! Reference GPH
      type (VectorValue_T), intent(in) :: TEMPERATURE ! Temperature

      real(rp), dimension(quantity%template%noInstances) :: CP2, CSQ, REQ, SP2
      real(rp) :: PhiCorrs(temperature%template%noInstances,temperature%template%noSurfs)
      real(rp), dimension(temperature%template%noInstances) :: REQS
      real(rv), dimension(temperature%template%noSurfs) :: Heights, N, PhiCorr, PS
      integer :: I, J ! Subscripts, loop inductors

      ! Executable code
      ! More sanity checks
      if ( quantity%template%instrumentModule /= ptan%template%instrumentModule ) &
        & call Announce_Error ( key, No_Error_Code, &
        & 'PHITan and PTan quantities are not for the same module' )
      if ( any(shape(quantity%values)/=shape(ptan%values)) ) &
        & call Announce_Error ( key, No_Error_Code, &
        & 'PHITan and PTan quantities are not the same size' )
      if ( .not. ValidateVectorQuantity ( temperature, &
        & quantityType=(/l_temperature/), coherent=.true., stacked=.true., &
        & frequencyCoordinate=(/l_none/), verticalCoordinate=(/l_zeta/) ) ) &
        & call Announce_error ( key, no_error_code, 'Problem with temperature quantity for phiTan fill' )
      if ( .not. ValidateVectorQuantity ( h2o, &
        & quantityType=(/l_vmr/), molecule=(/l_h2o/), coherent=.true., stacked=.true., &
        & frequencyCoordinate=(/l_none/), verticalCoordinate=(/l_zeta/) ) ) &
        & call Announce_error ( key, no_error_code, 'Problem with temperature quantity for phiTan fill' )
      if ( .not. ValidateVectorQuantity ( refGPH, &
        & quantityType = (/l_refGPH/), coherent=.true., stacked=.true., &
        & verticalCoordinate=(/l_zeta/), frequencyCoordinate=(/l_none/), noSurfs=(/1/) ) ) &
        & call Announce_Error ( key, badrefGPHQuantity )

      i = 0
      if ( .not. DoHGridsMatch ( temperature, refGPH ) ) i = i + 1
      if ( .not. DoHGridsMatch ( temperature, h2o ) ) i = i + 2
      if ( quantity%template%noInstances /= orbincline%template%noInstances ) i = i + 4
      if ( .not. DoVGridsMatch ( temperature, h2o ) ) i = i + 8
      if ( i /= 0 ) call Announce_Error ( key, no_error_code, &
        & ' coordinates for temperature/refGPH/H2O disagree', &
        & extraInfo = (/i/) )

      ! compute equivalent earth radius at uncorrected tangent phi nearest surface
      csq = (earthrada * earthradb)**2 / &
            & ((earthrada**2-earthradb**2)*sin(orbincline%values(1,:)*deg2rad)**2 + &
            &   earthradb**2)
      cp2 = cos(quantity%values(1,:)*deg2rad)**2
      sp2 = 1.0_rp - cp2
      req = 0.001_rp*sqrt((earthrada**4*sp2 + csq**2*cp2) / &
                        & (earthrada**2*cp2 + csq*sp2))

      ! Interpolate REQ to temperature/h2o/refGPH hGrid
      call InterpolateValues ( quantity%values(1,:), req, &
        & temperature%template%phi(1,:), reqs, 'Linear', extrapolate='Constant' )

      ! Get pressures corresponding to Temperature etc
      ps = 10.0**(-temperature%template%surfs(:,1))

      ! OK, do the refraction calculation on Temperature's grids
      do i = 1, temperature%template%noInstances

        ! Get heights.  Temperature and RefGPH are on same hGrids.
        ! RefGPH is in meters, but Hydrostatic wants it in km.
        call Hydrostatic ( GeodToGeocLat ( temperature%template%geodLat(1,i) ), &
          & temperature%template%surfs(:,1), temperature%values(:,i), &
          & temperature%template%surfs(:,1), &
          & refGPH%template%surfs(1,1), 0.001*refGPH%values(1,i), &
          & heights )
        heights = heights + reqs(i)
        ! Get refractive indices.  Temperature and H2O are on same grids.
        call refractive_index ( ps, temperature%values(:,i), n, h2o%values(:,i) )
        do j = 1, temperature%template%noSurfs
          call phi_refractive_correction_up ( n(j:), heights(j:), phiCorr(j:) )
          phiCorrs(i,j) = phiCorr(temperature%template%noSurfs) * rad2deg
        end do
      end do

      ! Now interpolate Phi corrections onto the PhiTan grid and apply them.
      ! Its zeta grid is the value of PTan; its phi grid is its own value.

      do j = 1, size(quantity%values,2)
        call interpolateValues ( h2o%template%phi(1,:), quantity%values(:,j), &
          &                      h2o%template%surfs(:,1), ptan%values(:,j), &
          &                      phiCorrs, quantity%values(:,j), update=.true. )
      end do ! j

    end subroutine FillPhiTanWithRefraction

      ! ------------------------------------- FillIWCFromExtinction ----
    subroutine FillIWCFromExtinction ( quantity, &
     & sourceQuantity, temperatureQuantity)
      ! Actually, the meaning of the next two is reversed if invert is TRUE)
      type (VectorValue_T), intent(inout) :: QUANTITY ! (IWC) Quantity to fill
      type (VectorValue_T), intent(in) :: sourceQuantity ! cloud extinction
      type (VectorValue_T), intent(in) :: temperatureQuantity ! T(zeta)

      ! local variables
      real (r8), dimension(Temperaturequantity%template%noSurfs,quantity%template%noInstances) :: temp
      real (r8), dimension(Temperaturequantity%template%noSurfs) :: zt, yt
      real (r8), dimension(sourceQuantity%template%noSurfs) :: ze, ye
      real (r8), dimension(quantity%template%noSurfs) :: z, Tz, Ez, iwc0
      real (r8), dimension(quantity%template%noInstances) :: x2, y2
      real (r8), dimension(Temperaturequantity%template%noInstances) :: x1, y1
      integer :: i

      if ( .not. (quantity%template%coherent .and. sourceQuantity%template%coherent &
         .and. Temperaturequantity%template%coherent)) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "one of the quantities is not coherent")
      if ( sourceQuantity%template%noInstances /= Quantity%template%noInstances) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "IWC and Extinction profile numbers are not matched")

      if ( TemperatureQuantity%template%noInstances /= Quantity%template%noInstances ) then
         x1 = TemperatureQuantity%template%phi(1,:)
         x2 = quantity%template%phi(1,:)
         do i=1,Temperaturequantity%template%noSurfs
           y1 = TemperatureQuantity%values(i, :)
           call InterpolateValues( x1, y1, x2, y2, method='Linear' )
           temp(i,:) = y2
         end do
      else
         temp = TemperatureQuantity%values
      end if

      if ( quantity%template%verticalCoordinate == l_pressure ) then
        z = -log10 ( quantity%template%surfs(:,1) )
      else
        z = quantity%template%surfs(:,1)
      end if

      if ( sourceQuantity%template%verticalCoordinate == l_pressure ) then
        ze = -log10 ( sourceQuantity%template%surfs(:,1) )
      else
        ze = sourceQuantity%template%surfs(:,1)
      end if

      if ( Temperaturequantity%template%verticalCoordinate == l_pressure ) then
        zt = -log10 ( Temperaturequantity%template%surfs(:,1) )
      else
        zt = Temperaturequantity%template%surfs(:,1)
      end if

      do i=1, quantity%template%noInstances
        yt = Temp(:, i)
        call InterpolateValues( zt, yt, z,  Tz, method='Linear' )
        ye = sourceQuantity%values(:, i)*1000._r8  ! converted to 1/km
        call InterpolateValues( ze, ye, z, Ez, method='Linear' )
        ! see ATBD for the conversion based on the size distribution from
        ! McFarquhar and Heymsfield [1996]
        !
        iwc0 = 10**(-2.77+0.01*(Tz-273.15))
        where(Ez < 0._r8) Ez = 0._r8
        Ez = Ez**(1./1.4)

        quantity%values(:,i) = iwc0*Ez
      end do

    end subroutine FillIWCFromExtinction

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
            if ( .not. interpolate ) then
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
!MJF
    ! ------------------------------------- FillNoRadsPerMIF -----
    subroutine FillNoRadsPerMif ( key, quantity, measQty )
      integer, intent(in) :: KEY
      type(VectorValue_T), intent(inout) :: QUANTITY
      type(VectorValue_T), intent(in) :: MEASQTY
      ! Local variables
      integer :: MIF, MAF               ! Loop counters
      integer :: I0, I1                 ! Indices

      ! Executable code
      ! Do some fairly limited checking.
      if ( .not. ValidateVectorQuantity ( measQty, quantityType=(/l_radiance/), &
        & signal=(/quantity%template%signal/), sideband=(/quantity%template%sideband/) ) ) &
        & call Announce_Error ( key, no_error_code, &
        & 'Quantity and measurement quantity disagree' )
      if ( associated ( measQty%mask ) ) then
        do maf = 1, quantity%template%noInstances
          do mif = 1, quantity%template%noSurfs
            i0 = 1 +  ( mif-1 ) * measQty%template%noChans
            i1 = i0 + measQty%template%noChans - 1
            quantity%values ( mif, maf ) = count ( &
              & iand ( ichar ( measQty%mask ( i0:i1, maf ) ), m_linAlg ) == 0 )
          end do
        end do
      else
        quantity%values = measQty%template%noChans
      end if
    end subroutine FillNoRadsPerMIF

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
      !
      ! Note: If either of TempPrecision or H2OPrecision is negative, then
      ! set Precision to negative, too (if NEGATIVETOO is TRUE)
      integer, intent(in) :: key          ! For messages
      ! Actually, the meaning of the next two is reversed if invert is TRUE)
      type (VectorValue_T), intent(inout) :: QUANTITY ! (rhi) Quantity to fill
      type (VectorValue_T), intent(in) :: sourcePrecisionQuantity ! vmr (unless invert)
      type (VectorValue_T), intent(in) :: tempPrecisionQuantity ! T(zeta)
      type (VectorValue_T), intent(in) :: sourceQuantity ! vmr (unless invert)
      type (VectorValue_T), intent(in) :: temperatureQuantity ! T(zeta)
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
      integer ::                          N           ! Num. of summed values
      logical, parameter ::               NEGATIVETOO = .true.
      integer ::                          QINDEX
      real (r8) ::                        rhi_precision
      integer ::                          S           ! Surface loop counter
      integer ::                          S_H2O       ! Instance num for surfs
      integer ::                          S_RHI       ! Instance num for surfs
      integer ::                          S_T         ! Instance num for surfs
      logical ::                          skipMe
      character(len=*), parameter ::      VMR_UNITS = 'vmr'
      integer ::                          VMR_UNIT_CNV
      logical ::                          wereAnySkipped
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
            if ( .not. interpolate ) then
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
            ! But skip no matter what else if temperature illegal
            skipMe = skipMe .or. TofZeta(s) <= 0.0
            if ( .not. skipMe ) then
              call RHIPrecFromH2O( H2OofZeta(s), &
               & TofZeta(s), zeta(qIndex), vmr_unit_cnv, &
               & H2OPrecisionofZeta(s), TPrecisionofZeta(s), &
               & rhi_precision, negativeToo )
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

      if ( .not. ValidateVectorQuantity ( quantity, &
        & quantityType=(/l_radiance/), minorFrame=.true.) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Invalid quantity for estimated noise fill.")

      if ( .not. ValidateVectorQuantity ( radiance, &
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

      if ( .not. ValidateVectorQuantity ( &
        & sysTemp, &
        & quantityType=(/l_systemTemperature/), &
        & verticalCoordinate=(/l_none/), &
        & noInstances=(/1/) ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Invalid system temperature quantity for estimated noise fill.")

      if ( associated ( nbw ) ) then
        if ( .not. ValidateVectorQuantity ( &
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

    ! ------------------------------------- FillVectorQtyHydrostatically ----
    subroutine FillVectorQtyHydrostatically ( key, quantity, &
      & temperatureQuantity, refGPHQuantity, h2oQuantity, &
      & orbitInclinationQuantity, phiTanQuantity, geocAltitudeQuantity, maxIterations, &
      & phiWindow, phiWindowUnits, chunkNo )
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
      integer, intent(in), optional :: chunkNo
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
          &   geocAltitudeQuantity%template%instrumentModule) ) then
          call Announce_Error ( key, nonConformingHydrostatic, &
            & "case l_ptan failed third test" )
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
          & phiWindow, phiWindowUnits, chunkNo )
      case default
        call Announce_error ( 0, 0, 'No such fill yet' )
      end select

      if ( toggle(gen) .and. levels(gen) > 0 ) &
        & call trace_end ( "FillVectorQtyHydrostatically" )

    end subroutine FillVectorQtyHydrostatically

    ! ------------------------------------- FillFromSplitSideband ----
    subroutine FillFromSplitSideband ( quantity, sourceQuantity, &
      & lsbFraction, usbFraction, spreadFlag, usb, channel, key )

      type (VectorValue_T), intent(inout) :: QUANTITY
      type (VectorValue_T), intent(in) :: SOURCEQUANTITY
      type (VectorValue_T), intent(in) :: LSBFRACTION
      type (VectorValue_T), intent(in) :: USBFRACTION
      type (VectorValue_T), pointer :: USB
      logical, intent(in) :: SPREADFLAG   ! One instance given, spread to all
      integer, intent(in) :: CHANNEL
      integer, intent(in) :: KEY
      ! Local variables
      integer :: MYCHANNEL              ! Possibly offset channel
      integer :: i, mif, maf
      type (Signal_T) :: signalIn, signalOut, signalRef
      real(r8), dimension(:), pointer :: freq1, freqL1, freqU1
      real(r8), dimension(:), pointer :: freq2, freqL2, freqU2
      real(r8), dimension(:), pointer :: freq, freqL, freqU
      real(r8) :: ratio1, ratio2    ! signal sideband fractions
      real(r8) :: scaledRad   ! scaled radiance according to the f^4 law

      ! Executable code
      nullify(freq, freqL, freqU, freq1, freqL1, freqU1, freq2, freqL2, freqU2)

      ! check for qualified quantity
      if ( .not. ValidateVectorQuantity ( quantity, quantityType=(/l_cloudInducedRadiance/), &
        & sideband=(/-1,1/), minorFrame=.true. )) &
        & call Announce_Error ( key, no_error_code, 'Quantity must be cloud-induced-radiances to fill' )
      if ( .not. associated(quantity%mask)) &
        & call Announce_Error ( key, no_error_code, 'Quantity must be a subset to fill' )

      signalIn = GetSignal ( sourceQuantity%template%signal )     ! sideband info gets lost
      signalOut = GetSignal ( Quantity%template%signal )          ! sideband info gets lost

      myChannel = channel

      call allocate_test ( freq1, size(signalIn%frequencies), 'frequencies', ModuleName )
      call allocate_test ( freqL1, size(signalIn%frequencies), 'LSBfrequencies', ModuleName )
      call allocate_test ( freqU1, size(signalIn%frequencies), 'USBfrequencies', ModuleName )
      call allocate_test ( freq2, size(signalOut%frequencies), 'frequencies', ModuleName )
      call allocate_test ( freqL2, size(signalOut%frequencies), 'LSBfrequencies', ModuleName )
      call allocate_test ( freqU2, size(signalOut%frequencies), 'USBfrequencies', ModuleName )

      ! find input signal frequencies
      freqL1 = signalIn%lo - signalIn%centerFrequency - &
        & signalIn%direction*signalIn%frequencies    ! lower sideband freq
      freqU1 = signalIn%lo + signalIn%centerFrequency + &
        & signalIn%direction*signalIn%frequencies     ! upper sideband freq
      if ( sourceQuantity%template%sideband == -1) freq1 = freqL1
      if ( sourceQuantity%template%sideband == 1) freq1 = freqU1
      ! find output signal frequencies
      freqL2 = signalOut%lo - signalOut%centerFrequency - &
        & signalOut%direction*signalOut%frequencies      ! lower sideband freq
      freqU2 = signalOut%lo + signalOut%centerFrequency + &
        & signalOut%direction*signalOut%frequencies    ! upper sideband freq
      if ( quantity%template%sideband == -1) freq2 = freqL2
      if ( quantity%template%sideband == 1) freq2 = freqU2

      quantity%values=0._r8
      if ( spreadFlag ) then
        ! spread a cloudy radiance to other bands according to the f^4 law
        ! The source cloudy radiance must be from a single sideband signal
        if ( .not. ValidateVectorQuantity ( sourceQuantity, &
          & quantityType=(/l_cloudInducedRadiance/), &
          & sideband=(/-1,1/), minorFrame=.true. )) &
          & call Announce_Error ( key, no_error_code, 'Inappropriate sourceQuantity radiance for fill' )
        do i=1,size(signalOut%frequencies)
          do maf=1, quantity%template%noInstances
            do mif=1, quantity%template%noSurfs
              scaledRad = sourceQuantity%values ( MyChannel + &
                & (mif-1) * size(signalIn%frequencies), maf) * &
                & freq2(i)**4/freq1(MyChannel)**4
              if ( iand ( ichar ( Quantity%mask(i+(mif-1) * &
                & size(signalOut%frequencies), maf)), m_cloud) /= 0 ) &
                & quantity%values(i+(mif-1)*size(signalOut%frequencies), maf) = scaledRad
            end do
          end do
        end do

      else
        ! split two sideband radiances according to the f^4 law
        ! The source cloudy radiance must be from the same signal of double sideband

        if ( .not. ValidateVectorQuantity ( sourceQuantity, &
          & quantityType=(/l_cloudInducedRadiance/), &
          & sideband=(/0/), signal=(/quantity%template%signal/), minorFrame=.true. )) &
          & call Announce_Error ( key, no_error_code, 'Inappropriate sourceQuantity radiance for fill' )
        if ( .not. ValidateVectorQuantity ( lsbFraction, &
          & quantityType=(/l_limbSidebandFraction/), &
          & signal=(/quantity%template%signal/), sideband=(/-1/) ) ) &
          & call Announce_Error ( key, no_error_code, 'Inappropriate lsbFraction quantity for fill' )
        if ( .not. ValidateVectorQuantity ( usbFraction, &
          & quantityType=(/l_limbSidebandFraction/), &
          & signal=(/quantity%template%signal/), sideband=(/1/) ) ) &
          & call Announce_Error ( key, no_error_code, 'Inappropriate usbFraction quantity for fill' )

        ! This method is only appliable to the cloud-induced radiances. One may assume
        ! that the scattering radiances are proportional to f^4. And this operation
        ! is only applied to maskbit = m_cloud.

        if ( .not. associated(usb) ) then
          ! If both sidebands are within 20GHz and have similar penetration depths
          ! in the attenuative atmosphere. We can neglect the attenuation and split cloudy
          ! radiances assuming that they are fully due to scattering and obey the f^4 law.
          do i=1,size(signalOut%frequencies)
            do maf=1, quantity%template%noInstances
              do mif=1, quantity%template%noSurfs
                if ( iand ( ichar ( Quantity%mask( i + &
                  & (mif-1) * size(signalOut%frequencies), maf)), m_cloud) /= 0 ) then
                  quantity%values(i+(mif-1)*size(signalOut%frequencies), maf ) = &
                    & sourceQuantity%values ( myChannel + &
                    & (mif-1)*size(signalIn%frequencies), maf) *freq2(i)**4/ &
                    & (lsbFraction%values ( myChannel, 1 ) * freqL1(MyChannel)**4 + &
                    & usbFraction%values ( myChannel, 1 ) * freqU1(MyChannel)**4)
                end if
              end do
            end do
          end do
        else
          ! If both sidebands are within 20GHz but have very different penetration depths,
          ! where one is optically thick and one is optically thin. We may use a reference
          ! cloud radiance near the optically-thin sideband to split cloud radiance in the
          ! optically-thick sideband. The reference radiance must be a single sideband
          ! radiance in this case, namely usb (upper sideband in most cases). The usb is
          ! scaled to the sourceQty upperside (thin) via the f^4 law, and the rest cloud
          ! radiance is for the lower sideband.
          !
          ! In this case, channel will for the reference signal and sideband frac remains for
          ! source signal.
          !
          if ( .not. ValidateVectorQuantity ( usb, &
            & quantityType=(/l_cloudInducedRadiance/), sideband=(/-1,1/), &
            & minorFrame=.true. )) call Announce_Error ( key, no_error_code, &
            & 'Inappropriate reference radiance for splitting' )

          signalRef = GetSignal ( usb%template%signal )

          call allocate_test ( freq, size(signalRef%frequencies), 'frequencies', ModuleName )
          call allocate_test ( freqL, size(signalRef%frequencies), &
            & 'LSBfrequencies', ModuleName )
          call allocate_test ( freqU, size(signalRef%frequencies), &
            & 'USBfrequencies', ModuleName )
          freq = signalRef%centerFrequency + signalRef%direction*signalRef%frequencies
          freqL = signalOut%lo - freq    ! lower sideband freq
          freqU = signalOut%lo + freq    ! upper sideband freq
          if ( usb%template%sideband == -1) freq = freqL
          if ( usb%template%sideband == 1) freq = freqU
          do i=1,size(signalOut%frequencies)

            ! need to scale the opposite sideband for the output signal
            if ( quantity%template%sideband == 1 ) then
              ratio1 = lsbFraction%values(i,1)
              ratio2 = usbFraction%values(i,1)
              freq2 = freqL2
            end if
            if ( quantity%template%sideband == -1 ) then
              ratio2 = lsbFraction%values(i,1)
              ratio1 = usbFraction%values(i,1)
              freq2 = freqU2
            end if

            do maf=1, quantity%template%noInstances
              do mif=1, quantity%template%noSurfs
                scaledRad = usb%values(MyChannel+(mif-1)*size(signalIn%frequencies), maf) * &
                  & freq2(i)**4/freq(MyChannel)**4

                if ( iand ( ichar ( Quantity%mask( i + &
                  & (mif-1)*size(signalOut%frequencies), maf)), m_cloud) /= 0 ) then
                  quantity%values(i+(mif-1)*size(signalOut%frequencies), maf) = &
                    &  (sourceQuantity%values(i+(mif-1)*size(signalIn%frequencies), maf) - &
                    &  ratio1*scaledRad)/ratio2
                end if
              end do
            end do
          end do

          call deallocate_test ( freq, 'frequencies', ModuleName )
          call deallocate_test ( freqL,'LSBfrequencies', ModuleName )
          call deallocate_test ( freqU,'USBfrequencies', ModuleName )
        end if
      end if

      call deallocate_test ( freq1, 'frequencies', ModuleName )
      call deallocate_test ( freqL1,'LSBfrequencies', ModuleName )
      call deallocate_test ( freqU1,'USBfrequencies', ModuleName )
      call deallocate_test ( freq2, 'frequencies', ModuleName )
      call deallocate_test ( freqL2,'LSBfrequencies', ModuleName )
      call deallocate_test ( freqU2,'USBfrequencies', ModuleName )
    end subroutine FillFromSplitSideband

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
        call Announce_error ( 0, no_error_code, 'GPH precision needed for result of FillGPHPrecision' )
      end select

      if ( toggle(gen) .and. levels(gen) > 0 ) &
        & call trace_end ( "FillGPHPrecision" )

    end subroutine FillGPHPrecision

    ! -------------------------------------- FillVectorQtyFromIsotope -----------

    subroutine FillVectorQtyFromIsotope ( quantity, sourceQuantity, &
              & ratioQuantity )
      ! This routine fills one vector from another, given an appropriate
      ! isotope ratio.

      type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      type (VectorValue_T), intent(in) :: SOURCEQUANTITY ! Quantity to take vmr from
      type (VectorValue_T), intent(in) :: RATIOQUANTITY ! Isotope ratio information

      ! Local variables
      real (r8) :: FACTOR                 ! Multiplier to apply to sourceQuantity

      ! Executable code

      if ( .not. ValidateVectorQuantity ( quantity, &
        & quantityType=(/ l_vmr /), frequencyCoordinate=(/ l_none /) ) ) &
        &   call MLSMessage ( MLSMSG_Error, ModuleName, &
        &      "Inappropriate quantity for isotope fill")

      if ( .not. ValidateVectorQuantity ( sourceQuantity, &
        & quantityType=(/ l_vmr /), frequencyCoordinate=(/ l_none /) ) ) &
        &   call MLSMessage ( MLSMSG_Error, ModuleName, &
        &      "Inappropriate source quantity for isotope fill")

      if ( .not. ValidateVectorQuantity ( ratioQuantity, &
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

    ! ----------------------------------- FillQuantityFromASCIIFile --------
    subroutine FillQuantityFromAsciiFile ( key, quantity, filename, badRange )
      use IO_stuff, only: GET_LUN
      use Machine, only: IO_Error
      use MoreMessage, only: MLSMessage
      integer, intent(in) :: KEY        ! Tree node
      type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      integer, intent(in) :: FILENAME   ! ASCII filename to read from
      real(r8), dimension(2), optional, intent(in) :: BADRANGE ! Range for missing data
      ! Local variables
      integer :: LUN                    ! Unit number
      integer :: STATUS                 ! Flag from open/close/read etc.
      character(len=1024) :: FILENAMESTR
      ! Executable code
      call get_lun ( lun, msg=.false. )
      if ( lun < 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "No logical unit numbers available" )
      call Get_String ( filename, filenameStr, strip=.true. )
      open ( unit=lun, file=filenameStr, status='old', form='formatted', &
        & access='sequential', iostat=status )
      if ( status /= 0 ) then
        call io_Error ( "Unable to open ASCII input file ", status, filenameStr )
        call MLSMessage( MLSMSG_Error, ModuleName, 'Error opening ASCII file at %l', &
          & (/ source_ref(key) /) )
      end if
      read ( unit=lun, fmt=*, iostat=status ) quantity%values
      if ( status /= 0 ) then
        call io_Error ( "Unable to read ASCII input file ", status, filenameStr )
        call MLSMessage( MLSMSG_Error, ModuleName, 'Error reading ASCII file %l', &
          & (/ source_ref(key) /) )
      end if
      close ( unit=lun, iostat=status )
      if ( status /= 0 ) then
        call io_Error ( "Unable to close ASCII input file ", status, filenameStr )
        call MLSMessage( MLSMSG_Error, ModuleName, 'Error closing ASCII file at %l', &
          & (/ source_ref(key) /) )
      end if
      if ( present ( badRange ) ) then
        if ( .not. associated ( quantity%mask ) ) call CreateMask ( quantity )
        where ( quantity%values >= badRange(1) .and. &
          & quantity%values <= badRange(2) )
          quantity%mask = char(ior(ichar(quantity%mask),m_linAlg))
        end where
      end if
    end subroutine FillQuantityFromAsciiFile

    ! --------------------------------------------- RotateMagneticField ----
    subroutine RotateMagneticField ( key, qty, fieldECR, ecrToFOV )
      use Intrinsic, only: L_FIELDAZIMUTH, L_FIELDELEVATION, L_FIELDSTRENGTH
      use Units, only: RAD2DEG
      integer, intent(in) :: KEY        ! Where are we in the l2cf?
      type (VectorValue_T), intent(inout) :: QTY ! The quantity to fill
      type (VectorValue_T), intent(in) :: FIELDECR ! The input field
      type (VectorValue_T), intent(in) :: ECRTOFOV ! The rotation matrix
      ! Local variables
      integer :: INSTANCE               ! Loop counter
      integer :: SURFACE                ! Loop counter
      integer :: MAF(1)                 ! Which MAF is the best match to this instance
      real(r8) :: THISFIELD(3)          ! A magnetic field vector
      real(r8) :: ROTATION(3,3)         ! One rotation matrix
      real(r8) :: STRENGTH              ! The field strength

      ! Executable code
      ! Do some sanity checks
      if ( .not. any ( qty%template%quantityType == &
        & (/ l_fieldStrength, l_fieldAzimuth, l_fieldElevation /) ) ) then
        call Announce_Error ( key, no_error_code, 'Inappropriate quantity for this fill' )
        return
      end if
      if ( .not. DoHGridsMatch ( qty, fieldECR ) ) then
        call Announce_Error ( key, no_error_code, &
          & 'Field and result quantity must have matching hGrids' )
        return
      end if
      if ( .not. DoVGridsMatch ( qty, fieldECR ) ) then
        call Announce_Error ( key, no_error_code, &
          & 'Field and result quantity must have matching hGrids' )
        return
      end if

      do instance = 1, qty%template%noInstances
        ! Identify the relevant MAF and pull out the first MIF's rotation matrix
        ! (first MIF is probably close enough to all the others to be useful).
        maf = minloc ( abs ( qty%template%phi(1,instance) - &
          & ecrToFOV%template%phi(1,:) ) )
        rotation = reshape ( ecrToFOV%values( 1:9, maf), (/ 3, 3 /) )
        ! Now loop over the pressure levels in the input and output field information
        do surface = 1, qty%template%noSurfs
          ! Get the field in ECR coordinates
          thisField = fieldECR%values ( (surface-1) * 3 + 1 : surface*3, &
            & instance )
          ! Now rotate it into IFOVPP coordinates
          thisField = matmul ( rotation, thisField )
          ! Now work out the strength / angles as appropriate.
          strength = sqrt ( sum ( thisField ** 2 ) )
          select case ( qty%template%quantityType )
          case ( l_fieldStrength )
            qty%values(surface,instance) = strength
          case ( l_fieldElevation )
            if ( strength /= 0.0_r8 ) then
              ! Elevation is constrained to 0--90 degrees instead of 0--180 degrees because
              ! radiative transfer Physics is symmetric.  We save half of l2pc bins.
              qty%values(surface,instance) = acos ( abs ( thisField(3) / strength ) ) * rad2deg
            else
              qty%values(surface,instance) = 0.0
            end if
          case ( l_fieldAzimuth )
            if ( thisField(1) /= 0.0_r8 ) then
              qty%values(surface,instance) = atan ( thisField(2) / thisField(1) ) * rad2deg
            else
              qty%values(surface,instance) = merge ( 90.0_r8, -90.0_r8, &
                & thisField(1) > 0.0_r8 )
            end if
          end select
        end do
      end do

    end subroutine RotateMagneticField

    !=============================================== ExplicitFillVectorQuantity ==
    subroutine ExplicitFillVectorQuantity ( quantity, valuesNode, spreadFlag, &
      & globalUnit, dontmask, &
      & AzEl, options, FillValue )

      ! This routine is called from MLSL2Fill to fill values from an explicit
      ! fill command line

      ! Dummy arguments
      type (VectorValue_T), intent(inout) :: QUANTITY ! The quantity to fill
      integer, intent(in) :: VALUESNODE   ! Tree node
      logical, intent(in) :: SPREADFLAG   ! One instance given, spread to all
      integer, intent(in) :: GLOBALUNIT   ! From parent vector
      logical, intent(in) :: DONTMASK     ! Don't bother with the mask
      logical, intent(in), optional :: AzEl ! Values are in [Mag, Az, El]; the
        ! desired quantity is components of Mag in the coordinate system to
        ! which Az and El are referenced.  So the number of values has to be
        ! a multiple of 3.
                                          ! (defaults to replacing all)
      ! The options are peculiar to this procedure, apart from verbose
      ! option           meaning
      ! ------           -------
      !   v              verbose
      !   e              replace only values in quantity == FillValue
      !   n              replace only values in quantity != FillValue
      !                   (defaults to replacing all)
      character (len=*), optional, intent(in) :: options ! E.g., '-v'
      real(r8), intent(in), optional :: FillValue

      ! Local variables
      integer :: K                        ! Loop counter
      integer :: I,J                      ! Other indices
      logical :: MyAzEl
      real(kind(quantity%values)) :: myFillValue
      character (len=8) :: myOptions
      integer :: NoValues
      integer :: TestUnit                 ! Unit to use
      integer, dimension(2) :: unitAsArray ! Unit for value given
      real (r8), pointer, dimension(:) :: VALUES
      real (r8), dimension(2) :: valueAsArray ! Value given
      logical :: Verbose
      character(len=2) :: whichToReplace ! '/=' (.ne. fillValue), '==', or ' ' (always)

      ! Executable code
      myAzEl = .false.
      if ( present(azEl) ) myAzEl = azEl
      myOptions = ' '
      if ( present(options) ) myOptions = options

      testUnit = quantity%template%unit
      if ( globalUnit /= phyq_Invalid ) testUnit = globalUnit
      noValues = nsons(valuesNode) - 1

      myFillValue = 0.
      if ( present(FillValue) ) myFillValue = FillValue

      whichToReplace = ' '
      if ( index(myOptions, 'e') > 0 ) then
        whichToReplace = '=='
      elseif ( index(myOptions, 'n') > 0 ) then
        whichToReplace = '/='
      end if
      verbose = ( index(myOptions, 'v') > 0 )
      ! Check the dimensions work out OK
      if ( myAzEl .and. mod(noValues,3) /= 0 ) &
          & call Announce_Error ( valuesNode, invalidExplicitFill )
      if ( spreadFlag ) then
        if ( noValues /= quantity%template%instanceLen .and. &
          & noValues /= quantity%template%noChans .and. &
          & noValues /= 1 ) &
          & call Announce_Error ( valuesNode, invalidExplicitFill )
      else
        if ( noValues /= &
          & quantity%template%instanceLen * quantity%template%noInstances ) &
          & call Announce_Error ( valuesNode, invalidExplicitFill, &
            & extraInfo = (/ &
              & quantity%template%instanceLen * quantity%template%noInstances /) )
      end if

      ! Get the values the user asked for, checking their units
      nullify ( values )
      call Allocate_test ( values, noValues, 'values', ModuleName )
      if ( .not. myAzEl ) then
        do k = 1, noValues
          call expr_check ( subtree(k+1,valuesNode) , unitAsArray, valueAsArray, &
            & (/testUnit, PHYQ_Dimensionless/), unitsError )
          if ( unitsError ) call Announce_error ( valuesNode, wrongUnits, &
            & extraInfo=(/unitAsArray(1), testUnit, PHYQ_Dimensionless/) )
          values ( k ) = valueAsArray(1)
        end do
      else
        ! Convert from Mag, Az, El to 3-D projections
        do k = 1, noValues, 3
          call expr_check ( subtree(k+1,valuesNode) , unitAsArray, valueAsArray, &
            & (/testUnit, PHYQ_Dimensionless/), unitsError )
          if ( unitsError ) call Announce_error ( valuesNode, wrongUnits, &
            & extraInfo=(/unitAsArray(1), testUnit, PHYQ_Dimensionless/) )
          values ( k ) = valueAsArray(1)
          ! Next two quantities have to be angles
          call expr_check ( subtree(k+2,valuesNode) , unitAsArray, valueAsArray, &
            & (/PHYQ_Angle/), unitsError )
          if ( unitsError ) call Announce_error ( valuesNode, wrongUnits, &
            & extraInfo=(/unitAsArray(1), PHYQ_Angle/) )
          values (k+1) = deg2rad * valueAsArray(1)
          call expr_check ( subtree(k+3,valuesNode) , unitAsArray, valueAsArray, &
            & (/PHYQ_Angle/), unitsError )
          if ( unitsError ) call Announce_error ( valuesNode, wrongUnits, &
            & extraInfo=(/unitAsArray(1), PHYQ_Angle/) )
          values (k+2) = deg2rad * valueAsArray(1)
          values(k:k+2) = values(k) * (/ cos(values(k+1))*cos(values(k+2)), &
                                         sin(values(k+1))*cos(values(k+2)), &
                                         sin(values(k+2)) /)
        end do
      end if

      if ( verbose ) then
        call output('Explicit fill of values for ', advance='no')
        call display_string ( quantity%template%name )
        call newline
        call output(values)
        call newline
      end if
      ! Now loop through the quantity
      k = 0
      do i = 1, quantity%template%noInstances
        do j = 1, quantity%template%instanceLen
          k = k + 1
          if ( .not. dontMask .and. associated ( quantity%mask ) ) then
            if ( iand ( ichar(quantity%mask(j,i)), m_Fill ) /= 0 ) cycle
          end if
          select case (whichToReplace)
          case ('/=')
            if ( quantity%values(j,i) == myFillValue ) cycle
          case ('==')
            if ( quantity%values(j,i) /= myFillValue ) cycle
          end select
          quantity%values(j,i) = values ( mod ( k-1, noValues ) + 1 )
        end do
      end do

      if ( verbose ) then
        call output(quantity%values(1,:))
        call newline
      end if
      ! Tidy up
      call Deallocate_test ( values, 'values', ModuleName )

    end subroutine ExplicitFillVectorQuantity

    ! ----------------------------------------- FillVectorQuantityFromL1B ----
    subroutine FillVectorQuantityFromL1B ( root, quantity, chunk, filedatabase, &
      & isPrecision, suffix, PrecisionQuantity, BOMask )
      use BitStuff, only: NegativeIfBitPatternSet
      use MLSFiles, only: HDFVERSION_5
      use MLSStrings, only: lowercase
      use output_m, only: blanks, output
      integer, intent(in) :: root
      type (VectorValue_T), INTENT(INOUT) ::        QUANTITY
      type (MLSChunk_T), INTENT(IN) ::              CHUNK
      type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
      logical, intent(in)               ::          ISPRECISION
      integer, intent(in), optional :: SUFFIX
      type (VectorValue_T), INTENT(IN), optional :: PRECISIONQUANTITY
      integer, intent(in), optional :: BOMask ! A pattern of bits--
                                              ! set prec. neg. if matched
      ! Local variables
      integer                               :: BO_error
      type (l1bData_T)                      :: BO_stat
      integer                               :: channel
      character (len=132)                   :: MODULENAMESTRING
      character (len=132)                   :: NAMESTRING
      integer                               :: FLAG, NOMAFS, maxMIFs
      type (l1bData_T)                      :: L1BDATA
      type (MLSFile_T), pointer             :: L1BFile
      type (MLSFile_T), pointer             :: L1BOAFile
      integer                               :: ROW, COLUMN
      integer                               :: myBOMask
      integer                               :: this_hdfVersion

      ! Executable code
      myBOMask = 0
      if ( present(BOMask) ) myBOMask = BOMask
      if ( toggle(gen) .and. levels(gen) > 0 ) &
        & call trace_begin ("FillVectorQuantityFromL1B",root)
      ! print *, 'Filling vector quantity from l1b'
      L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
      L1BOAFile => GetMLSFileByType(filedatabase, content='l1boa')
      this_hdfVersion = L1BFile%HDFVersion
      ! fileID = L1BFile%FileID%f_id

      select case ( quantity%template%quantityType )
      case ( l_ECRtoFOV )
        call GetModuleName( quantity%template%instrumentModule, nameString )
        nameString = AssembleL1BQtyName('ECRtoFOV', this_hdfVersion, .TRUE., &
          & trim(nameString))
      case ( l_L1BMAFBaseline )
        call GetSignalName ( quantity%template%signal, nameString, &
          & sideband=quantity%template%sideband, noChannels=.TRUE. )
        nameString = AssembleL1BQtyName(nameString, this_hdfVersion, .FALSE.)
      case ( l_l1bMIF_TAI )
        if ( this_hdfVersion == HDFVERSION_5 ) then
          call GetModuleName ( quantity%template%instrumentModule, nameString )
          nameString = AssembleL1BQtyName ('MIF_TAI', this_hdfVersion, .FALSE., &
            & trim(nameString) )
        else ! ??? MIF_TAI is goofy in HDF4 files -- no sc, no tp, no GHz....
          nameString = 'MIF_TAI'
        end if
      case ( l_LosVel )
        call GetModuleName ( quantity%template%instrumentModule, nameString )
        nameString = AssembleL1BQtyName ('LosVel', this_hdfVersion, .TRUE., &
          & trim(nameString) )
      case ( l_orbitInclination )
        nameString = AssembleL1BQtyName('OrbIncl', this_hdfVersion, .FALSE., &
          & 'sc')
      case ( l_ptan )
        call GetModuleName( quantity%template%instrumentModule,nameString )
        nameString = AssembleL1BQtyName('ptan', this_hdfVersion, .FALSE., &
          & trim(nameString))
      case ( l_radiance )
        call GetSignalName ( quantity%template%signal, nameString, &
          & sideband=quantity%template%sideband, noChannels=.TRUE. )
        nameString = AssembleL1BQtyName(nameString, this_hdfVersion, .FALSE.)
        L1BFile => GetL1bFile(filedatabase, namestring)
        ! fileID = FindL1BData (filedatabase, nameString )
      case ( l_scECI )
        nameString = AssembleL1BQtyName('ECI', this_hdfVersion, .FALSE., 'sc')
      case ( l_scGeocAlt )
        nameString = AssembleL1BQtyName('GeocAlt', this_hdfVersion, .FALSE., &
          & 'sc')
      case ( l_scVel )
        nameString = AssembleL1BQtyName('Vel', this_hdfVersion, .FALSE., 'sc')
      case ( l_scVelECI )
        nameString = AssembleL1BQtyName('VelECI', this_hdfVersion, .FALSE., &
          & 'sc')
      case ( l_scVelECR )
        nameString = AssembleL1BQtyName('VelECR', this_hdfVersion, .FALSE., &
          & 'sc')
      case ( l_tngtECI )
        call GetModuleName( quantity%template%instrumentModule,nameString )
        nameString = AssembleL1BQtyName('ECI', this_hdfVersion, .TRUE., &
          & trim(nameString))
      case ( l_tngtGeocAlt )
        call GetModuleName( quantity%template%instrumentModule,nameString )
        nameString = AssembleL1BQtyName('GeocAlt', this_hdfVersion, .TRUE., &
          & trim(nameString))
      case ( l_tngtGeodAlt )
        call GetModuleName( quantity%template%instrumentModule,nameString )
        nameString = AssembleL1BQtyName('GeodAlt', this_hdfVersion, .TRUE., &
          & trim(nameString))
      case default
        call Announce_Error ( root, cantFillFromL1B )
      end select

      ! Perhaps will need to read bright object status from l1bOA file
      if ( isPrecision .and. myBOMask /= 0 ) then
        call GetModuleName ( quantity%template%instrumentModule, moduleNameString )
        moduleNameString = AssembleL1BQtyName('BO_stat', this_hdfVersion, .TRUE., &
          & trim(moduleNameString))
        call ReadL1BData ( L1BOAFile, moduleNameString, BO_stat, noMAFs, &
          & flag=BO_error, firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
          & NeverFail= .true., &
          & dontPad=DONTPAD )
      end if

      if ( present ( suffix ) ) then
        if ( suffix /= 0 ) then
          call Get_String ( suffix, &
            & nameString(len_trim(nameString)+1:), strip=.true. )
          ! Look for the field again
          L1BFile => GetL1bFile(filedatabase, namestring)
          ! fileID = FindL1BData (filedatabase, nameString )
          ! Note we won't find it if it's in the OA file, I'm going to ignore that
          ! possibility for the moment.

        end if
      end if

      ! If the quantity exists (or it doesn't exist but it's not a radiance)
      if ( index(lowercase(namestring), 'baseline') > 0 .and. .false. ) then
        call output('namestring: ', advance='no')
        call output(trim(namestring), advance='yes')
        call output('associated(L1BFile): ', advance='no')
        call output(associated(L1BFile), advance='yes')
        call output('qty type: ', advance='no')
        call output(quantity%template%quantityType, advance='no')
        call blanks(3)
        call output(l_radiance, advance='no')
        call blanks(3)
        call output(l_L1BMAFBaseline, advance='yes')
      endif
      if ( associated(L1BFile) .or. ( quantity%template%quantityType /= l_radiance .and. &
        & quantity%template%quantityType /= l_L1BMAFBaseline ) ) then
        if ( isPrecision ) nameString = trim(nameString) // PRECISIONSUFFIX
        L1BFile => GetL1bFile(filedatabase, namestring)

        call ReadL1BData ( L1BFile, nameString, l1bData, noMAFs, flag, &
          & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
          & NeverFail= .false., &
          & dontPad=DONTPAD )
        ! If it didn't exist in the not-a-radiance case, then we'll fail here.
        if ( flag /= 0 ) then
          call Announce_Error ( root, errorReadingL1B )
          if ( toggle(gen) .and. levels(gen) > 0 ) &
            & call trace_end ( "FillVectorQuantityFromL1B")
          return
        end if
        if ( quantity%template%noInstances /= size ( l1bData%dpField, 3 ) .or. &
          &  quantity%template%instanceLen /= &
          &   size ( l1bData%dpField, 1 ) * size ( l1bData%dpField, 2 ) ) then
          call output ( 'Quantity shape:' )
          call output ( quantity%template%instanceLen )
          call output ( ' ( ' )
          call output ( quantity%template%noChans )
          call output ( ', ' )
          call output ( quantity%template%noSurfs )
          call output ( ' ), ' )
          call output ( quantity%template%noInstances, advance='yes' )
          call output ( 'L1B shape:' )
          call output ( size ( l1bData%dpField, 1 ) )
          call output ( ', ' )
          call output ( size ( l1bData%dpField, 2 ) )
          call output ( ', ' )
          call output ( size ( l1bData%dpField, 3 ), advance='yes' )
          call Announce_Error ( root, no_error_code, 'L1B data is wrong shape' )
          if ( toggle(gen) .and. levels(gen) > 0 ) &
            & call trace_end ( "FillVectorQuantityFromL1B")
          return
        end if

        if ( isPrecision .and. myBOMask /= 0 .and. BO_error == 0 ) then
          noMAFs = size(l1bData%dpField, 3)
          maxMIFs = l1bData%maxMIFs
          if ( switchDetail(switches, 'glob') > 0 ) then ! e.g., 'glob1'
            call output ( 'Quantity shape:' )
            call output ( quantity%template%instanceLen )
            call output ( ' ( ' )
            call output ( quantity%template%noChans )
            call output ( ', ' )
            call output ( quantity%template%noSurfs )
            call output ( ' ), ' )
            call output ( quantity%template%noInstances, advance='yes' )
            call output ( 'L1B shape:' )
            call output ( size ( l1bData%dpField, 1 ) )
            call output ( ', ' )
            call output ( size ( l1bData%dpField, 2 ) )
            call output ( ', ' )
            call output ( size ( l1bData%dpField, 3 ), advance='yes' )
            call output_name_v_pair( 'shape' // trim(namestring), shape(l1bData%dpField) )
            call output_name_v_pair( 'shape(BO_stat)', shape(BO_stat%intField) )
            call output_name_v_pair( 'noMAFs', noMAFs )
            call output_name_v_pair( 'maxMIFs', maxMIFs )
          endif
          do channel = 1, quantity%template%noChans
          l1bData%dpField(channel,:,:) = &
            & NegativeIfBitPatternSet( l1bData%dpField(channel,:,:), &
            & BO_stat%intField(1, 1:maxMIFs, 1:noMAFs), myBOMask )
          enddo
          call DeallocateL1BData(BO_stat)
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

        if ( switchDetail(switches, 'l1b') > -1 ) &
          & call Dump( l1bData )
        call DeallocateL1BData(l1bData)
      else
        ! This is the case where it's a radiance we're after and it's missing
        quantity%values = DEFAULTUNDEFINEDVALUE ! -1.0
        do column=1, size(quantity%values(1,:))
          do row=1, size(quantity%values(:,1))
            call MaskVectorQty ( quantity, row, column )
          end do
        end do
      end if
      if ( toggle(gen) .and. levels(gen) > 0 ) call trace_end( "FillVectorQuantityFromL1B" )
    end subroutine FillVectorQuantityFromL1B

    ! ------------------------------------------- FillVectorQuantityFromL2AUX --
    subroutine FillVectorQuantityFromL2AUX ( qty, l2aux, errorCode )
      type ( VectorValue_T), intent(inout) :: QTY
      type ( L2AUXData_T), intent(in) :: L2AUX
      integer, intent(inout) :: ERRORCODE
      ! Local variables
      integer :: FIRSTPROFILE
      integer :: LASTPROFILE
      ! Executable code
      errorCode = 0
      ! Work out which profile in the l2aux this belongs to
      firstProfile = qty%template%instanceOffset - qty%template%noInstancesLowerOverlap
      lastProfile = firstProfile + qty%template%noInstances - 1
      ! In the case of minor/major frame quanties, while instanceOffset is
      ! Numbered from zero, as in L1B, our array starts from 1.
      if ( qty%template%minorFrame .or. qty%template%majorFrame ) then
        firstProfile = firstProfile + 1
        lastProfile = lastProfile + 1
      end if
      ! Check that the dimensions are appropriate
      if ( firstProfile < lbound ( l2aux%values, 3 ) ) then
        errorCode = CantFillFromL2AUX
        return
      end if
      if ( lastProfile > ubound ( l2aux%values, 3 ) ) then
        errorCode = CantFillFromL2AUX
        return
      end if
      if ( size ( l2aux%values, 1 ) /= qty%template%noChans .or. &
        &  size ( l2aux%values, 2 ) /= qty%template%noSurfs ) then
        errorCode = CantFillFromL2AUX
        return
      end if
      ! Do the fill
      qty%values = reshape ( l2aux%values ( :, :,  &
        & firstProfile : lastProfile ), &
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
        call Announce_Error ( key, no_error_code, &
          & 'Quantity does not describe magnetic field' )
        return
      end if
      if ( .not. ValidateVectorQuantity ( gphQty, quantityType=(/l_gph/), &
        & frequencyCoordinate=(/ l_none /), verticalCoordinate=(/l_zeta/) ) ) then
        call Announce_Error ( key, no_error_code, &
          & 'GPH quantity does not describe gph field' )
        return
      end if
      if ( .not. DoHGridsMatch ( qty, gphQty ) ) then
        call Announce_Error ( key, no_error_code, &
          & 'Quantity and GPHQuanity do not share the same horizontal basis' )
        return
      end if
      if ( .not. DoVGridsMatch ( qty, gphQty ) ) then
        call Announce_Error ( key, no_error_code, &
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

      if ( switchDetail(switches,'mag') > -1 ) &
        & call dump ( qty, clean=(switchDetail(switches,'clean') > -1) )

    end subroutine FillQuantityUsingMagneticModel

    ! ------------------------------------------- FillQtyFromInterpolatedQty
    subroutine FillQtyFromInterpolatedQty ( qty, source, force, key )
      type (VectorValue_T), intent(inout) :: QTY
      type (VectorValue_T), intent(in) :: SOURCE
      logical, intent(in) :: FORCE
      integer, intent(in) :: KEY

      ! Local variables
      real (r8), dimension(:), pointer :: oldSurfs, newSurfs
      real (r8), dimension(:,:), pointer :: newValues
      logical :: mySurfs, myNewValues

      ! Executable code
      if ( .not. DoQtysDescribeSameThing ( qty, source ) .and. .not. force ) then
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

      ! Two cases here, one where we have to interpolate vertically (has to be
      ! single channel quantity).  The other where we don't.  Most of the latter
      ! cases can be handled by the code that calls this routine.  The exception
      ! is when we've used the force option to e.g. copy one band into another.
      if ( .not. doVGridsMatch ( qty, source ) ) then
        ! This quantity needs vertical interpolation
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
        end if

        ! Tidy up
        if ( mySurfs ) then
          call Deallocate_test ( oldSurfs, 'oldSurfs', ModuleName )
          call Deallocate_test ( newSurfs, 'newSurfs', ModuleName )
        end if
      else ! ------------------------
        ! No interpolation needed, more like the
        ! case handled in the calling code, except we're more lenient
        ! If we have a mask and we're going to obey it then do so
        if ( associated(quantity%mask) .and. .not. dontMask ) then
          where ( iand ( ichar(quantity%mask(:,:)), m_Fill ) == 0 )
            quantity%values(:,:) = sourceQuantity%values(:,:)
          end where
        else ! Otherwise, just blindly copy
          quantity%values = sourceQuantity%values
        end if
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
      !    if ( .not. ValidateVectorQuantity(qty, coherent=.TRUE., stacked=.TRUE., &
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
      if ( extinction ) then
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
              if ( colTrans > 0.02_r8) betaFine(i)= transFine(i)/colTrans
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
          if ( ptan%values(mif,maf) .gt. -2.5) cycle ! for testing
          ! find altitude of each s grid
          x_in = sLevel**2/2./(re%values(1,maf)*0.001_r8 + zt(mif))
          ! converted to zeta
          x_in = x_in/16. + ptan%values(mif,maf)
          ! find minimum and maximum pressures indices in sGrid
          do i = 2,qty%template%noSurfs-1
            if ( ptan%values(mif,maf) < (outZeta(i)+outZeta(i+1))/2. .and. &
              & ptan%values(mif,maf) > (outZeta(i)+outZeta(i-1))/2.) &
              & minZ = i
          end do
          if ( ptan%values(mif,maf) < (outZeta(1)+outZeta(2))/2.) minZ=1
          if ( ptan%values(mif,maf) > outZeta(qty%template%noSurfs)) cycle ! goto next mif

          do i = 2,qty%template%noSurfs-1
            if ( x_in(noDepths) < (outZeta(i)+outZeta(i+1))/2. .and. &
              & x_in(noDepths) > (outZeta(i)+outZeta(i-1))/2.) &
              & maxZ = i
          end do
          if ( x_in(noDepths) < (outZeta(1)+outZeta(2))/2.) cycle    ! goto next mif
          if ( x_in(noDepths) > outZeta(qty%template%noSurfs)) maxZ=qty%template%noSurfs

          ! get phi along path for each mif ( phi is in degree)
          y_in = los%template%phi(mif,maf) &
            & - atan(sLevel/(re%values(1,maf)*0.001_r8 + zt(mif)))*rad2deg
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
              if ( phi_out(i) .lt. &
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
    subroutine FillQuantityByManipulation ( quantity, a, b, manipulation, &
      & key, force, c )
      use String_table, only: GET_STRING
      type (VectorValue_T), intent(inout) :: QUANTITY
      type (VectorValue_T), pointer :: A
      type (VectorValue_T), pointer :: B
      integer, intent(in) :: MANIPULATION
      integer, intent(in) :: KEY        ! Tree node
      logical, intent(in) :: FORCE      ! If set throw caution to the wind
      real(rv) :: C                     ! constant "c" in manipulation
      ! Local parameters
      integer, parameter :: NO2WAYMANIPULATIONS = 8
      character(len=*), parameter :: VALID2WAYMANIPULATIONS ( NO2WAYMANIPULATIONS ) = (/ &
        & 'a+b     ', &
        & '(a+b)/2 ', &
        & 'a-b     ', &
        & 'a*b     ', &
        & 'a>b     ', &
        & 'a<b     ', &
        & 'a|b     ', &
        & 'a/b     ' /)
      integer, parameter :: NO1WAYMANIPULATIONS = 7
      character(len=*), parameter :: VALID1WAYMANIPULATIONS ( NO1WAYMANIPULATIONS ) = (/ &
        & '-a      ', &
        & '1/a     ', &
        & 'abs(a)  ', &
        & 'sign(a) ', &
        & 'exp(a)  ', &
        & 'log(a)  ', &
        & 'log10(a)' /)
      ! Local variables
      character (len=128) :: MSTR
      character (len=1) :: ABNAME
      logical :: OKSOFAR
      logical :: OneWay
      logical :: TwoWay
      logical :: USESC
      integer :: I
      logical :: NEEDSB
      integer :: NUMWAYS ! 1 or 2
      type (VectorValue_T), pointer :: AORB
      ! Executable code

      ! Currently we have a rather brain dead approach to this, so
      ! check that what the user has asked for, we can supply.
      call get_string ( manipulation, mstr, strip=.true. )
      mstr = lowercase(mstr)
      
      OneWay = any ( mstr == valid1WayManipulations )
      TwoWay = any ( mstr == valid2WayManipulations )
      usesc  = .not. ( OneWay .or. TwoWay ) .and. &
        & index(mstr, 'c') > 0
      if ( .not. ( OneWay .or. TwoWay .or. usesc ) ) then
        call Announce_Error ( key, no_error_code, 'Invalid manipulation' )
        return
      end if
      
      needsB = TwoWay .or. (usesC .and. (index(mstr, 'b') > 0) )

      ! Now check the sanity of the request.
      if ( OneWay .or. usesc ) then
        numWays = 1
        if ( .not. associated ( a ) ) then
          call Announce_Error ( key, no_error_code, &
            & 'You must supply the a quantity' )
          return
        end if
      elseif ( TwoWay ) then
        numWays = 2
        if ( .not. associated ( a ) .or. .not. associated ( b ) ) then
          call Announce_Error ( key, no_error_code, &
            & 'You must supply both a and b quantities' )
          return
        end if
      end if

      okSoFar = .true.
      do i = 1, numWays ! 2
        if ( i == 1 ) then
          aorb => a
          abName = 'a'
        else
          aorb => b
          abName = 'b'
        end if

        ! For minor frame quantities, check that we're on the same page
        if ( .not. force ) then
          if ( quantity%template%minorFrame ) then
            okSoFar = okSoFar .and. aorb%template%minorFrame .and. &
              & quantity%template%signal == aorb%template%signal .and. &
              & quantity%template%sideband == aorb%template%sideband .and. &
              & quantity%template%frequencyCoordinate == aorb%template%frequencyCoordinate
          else if ( mstr == 'a*b' ) then
            ! In this case, just check that the coordinate systems for these quantities match
            okSoFar = okSoFar .and. &
              & DoHGridsMatch ( quantity, aorb ) .and. &
              & DoVGridsMatch ( quantity, aorb ) .and. &
              & DoFGridsMatch ( quantity, aorb, sizeOnly=.true. )
          else
            ! For a+/-b these quantities must share a template
            okSoFar = okSoFar .and. quantity%template%name == aorb%template%name
          end if
        else
          okSoFar = okSoFar &
            & .and. quantity%template%noInstances == aorb%template%noInstances &
            & .and. quantity%template%instanceLen == aorb%template%instanceLen
        end if

        if ( .not. okSoFar ) then
          call Announce_Error ( key, no_error_code, &
            & abName // ' is not of the same (or close enough) type as quantity' )
          return
        end if
      end do

      ! OK do the simple work for now
      ! Later we'll do fancy stuff to parse the manipulation.
      select case ( mstr )
      ! The binary manipulations
      case ( 'a+b' )
        if ( .not. associated ( quantity%mask ) ) then
          quantity%values = a%values + b%values
        else
          where ( iand ( ichar(quantity%mask(:,:)), m_fill ) == 0 )
            quantity%values = a%values + b%values
          end where
        end if
      case ( '(a+b)/2' )
        if ( .not. associated ( quantity%mask ) ) then
          quantity%values = 0.5 * ( a%values + b%values )
        else
          where ( iand ( ichar(quantity%mask(:,:)), m_fill ) == 0 )
            quantity%values = 0.5 * ( a%values + b%values )
          end where
        end if
      case ( 'a-b' )
        if ( .not. associated ( quantity%mask ) ) then
          quantity%values = a%values - b%values
        else
          where ( iand ( ichar(quantity%mask(:,:)), m_fill ) == 0 )
            quantity%values = a%values - b%values
          end where
        end if
      case ( 'a*b' )
        if ( .not. associated ( quantity%mask ) ) then
          quantity%values = a%values * b%values
        else
          where ( iand ( ichar(quantity%mask(:,:)), m_fill ) == 0 )
            quantity%values = a%values * b%values
          end where
        end if
      case ( 'a>b' )
        if ( .not. associated ( quantity%mask ) ) then
          quantity%values = max ( a%values, b%values )
        else
          where ( iand ( ichar(quantity%mask(:,:)), m_fill ) == 0 )
            quantity%values = max ( a%values, b%values )
          end where
        end if
      case ( 'a<b' )
        if ( .not. associated ( quantity%mask ) ) then
          quantity%values = min ( a%values, b%values )
        else
          where ( iand ( ichar(quantity%mask(:,:)), m_fill ) == 0 )
            quantity%values = min ( a%values, b%values )
          end where
        end if
      case ( 'a|b' )
        if ( .not. associated ( quantity%mask ) ) then
          quantity%values = ior ( nint(a%values), nint(b%values) )
        else
          where ( iand ( ichar(quantity%mask(:,:)), m_fill ) == 0 )
            quantity%values = ior ( nint(a%values), nint(b%values) )
          end where
        end if
      case ( 'a/b' )
        if ( .not. associated ( quantity%mask ) ) then
          where ( b%values /= 0._rv )
            quantity%values = a%values / b%values
          end where
        else
          where ( iand ( ichar(quantity%mask(:,:)), m_fill ) == 0 .and. &
            & ( b%values /= 0._rv ) )
            quantity%values = a%values / b%values
          end where
        end if
      ! The unary manipulations
      case ( '-a'  )
        if ( .not. associated ( quantity%mask ) ) then
            quantity%values = -a%values
        else
          where ( iand ( ichar(quantity%mask(:,:)), m_fill ) == 0  )
            quantity%values = -a%values
          end where
        end if
      case ( '1/a'  )
        if ( .not. associated ( quantity%mask ) ) then
          where ( a%values /= 0._rv )
            quantity%values = 1./a%values
          end where
        else
          where ( iand ( ichar(quantity%mask(:,:)), m_fill ) == 0 .and. &
            & ( a%values /= 0._rv ) )
            quantity%values = 1./a%values
          end where
        end if
      case ( 'abs(a)'  )
        if ( .not. associated ( quantity%mask ) ) then
            quantity%values = abs(a%values)
        else
          where ( iand ( ichar(quantity%mask(:,:)), m_fill ) == 0  )
            quantity%values = abs(a%values)
          end where
        end if
      case ( 'sign(a)' )
        if ( .not. associated ( quantity%mask ) ) then
          where ( a%values /= 0._rv )
            quantity%values = sign(1._rv, a%values)
          end where
        else
          where ( iand ( ichar(quantity%mask(:,:)), m_fill ) == 0 .and. &
            & ( a%values /= 0._rv ) )
            quantity%values = sign(1._rv, a%values)
          end where
        end if
      case ( 'exp(a)'  )
        if ( .not. associated ( quantity%mask ) ) then
            quantity%values = exp(a%values)
        else
          where ( iand ( ichar(quantity%mask(:,:)), m_fill ) == 0  )
            quantity%values = exp(a%values)
          end where
        end if
      case ( 'log(a)' )
        if ( .not. associated ( quantity%mask ) ) then
          where ( a%values > 0._rv )
            quantity%values = log(a%values)
          end where
        else
          where ( iand ( ichar(quantity%mask(:,:)), m_fill ) == 0 .and. &
            & ( a%values > 0._rv ) )
            quantity%values = log(a%values)
          end where
        end if
      case ( 'log10(a)' )
        if ( .not. associated ( quantity%mask ) ) then
          where ( a%values > 0._rv )
            quantity%values = log10(a%values)
          end where
        else
          where ( iand ( ichar(quantity%mask(:,:)), m_fill ) == 0 .and. &
            & ( a%values > 0._rv ) )
            quantity%values = log10(a%values)
          end where
        end if
      case default
        ! This should be one of the cases which use the constant "c"
        call SimpleExprWithC( quantity, a, b, c, mstr )
      end select

    end subroutine FillQuantityByManipulation

      subroutine SimpleExprWithC( quantity, a, b, c, mstr )
      type (VectorValue_T), intent(inout) :: QUANTITY
      type (VectorValue_T), pointer :: A
      type (VectorValue_T), pointer :: B
      real(rv) :: C                     ! constant "c" in manipulation
      character (len=128) :: MSTR
        ! Evaluate mstr assuming it's of the form
        ! expr1 [op1 expr2]
        ! where each expr is either a primitive 'x' (one of {a, b, or c})
        ! or else ['('] 'x op y' [')']
        ! where 'op' is one of {+, -, *, /,<,>}

        ! Method:
        ! Progressively collapse all the '(..)' pairs into their values
        ! (stored in primitives database)
        ! until only primitives remain
        ! Then evaluate the primitives
        !
        ! Limitations:
        ! Does not check for unmatched parens or other illegal syntax
        integer, parameter :: MAXSTRLISTLENGTH = 128
        integer, parameter :: MAXNESTINGS=64 ! Max number of '(..)' pairs
        character(len=MAXSTRLISTLENGTH) :: collapsedstr
        integer :: level
        integer :: np ! number of primitives
        character(len=MAXSTRLISTLENGTH) :: part1
        character(len=MAXSTRLISTLENGTH) :: part2
        character(len=MAXSTRLISTLENGTH) :: part3
        character(len=4) :: vchar
        logical, parameter :: DEEBUG = .false.
        ! Executable
        if ( DeeBUG ) print *, 'mstr: ', trim(mstr)
        nullify(primitives)
        np = 0
        
        ! We're unable to ensure operator precedence
        ! so we'll attempt to identify multiplications and divisions
        ! and surround such subexpressions with extra parentheses
        
        call reorderPrecedence(mstr, collapsedstr)
        if ( DeeBUG ) then
          print *, 'incoming ', mstr
          print *, 'after reordering precedence ', collapsedstr
        endif
        mstr = collapsedstr
        
        ! 1st--make sure spaces surround each operator
        ! (It takes two steps for each to avoid threat of infinite loop)
        call ReplaceSubString( mstr, collapsedstr, '+', ' & ', &
          & which='all', no_trim=.true. )
        call ReplaceSubString( collapsedstr, mstr, '&', '+', &
          & which='all', no_trim=.false. )

        call ReplaceSubString( mstr, collapsedstr, '*', ' & ', &
          & which='all', no_trim=.true. )
        call ReplaceSubString( collapsedstr, mstr, '&', '*', &
          & which='all', no_trim=.false. )

        call ReplaceSubString( mstr, collapsedstr, '-', ' & ', &
          & which='all', no_trim=.true. )
        call ReplaceSubString( collapsedstr, mstr, '&', '-', &
          & which='all', no_trim=.false. )

        call ReplaceSubString( mstr, collapsedstr, '/', ' & ', &
          & which='all', no_trim=.true. )
        call ReplaceSubString( collapsedstr, mstr, '&', '/', &
          & which='all', no_trim=.false. )

        call ReplaceSubString( mstr, collapsedstr, '<', ' & ', &
          & which='all', no_trim=.true. )
        call ReplaceSubString( collapsedstr, mstr, '&', '<', &
          & which='all', no_trim=.false. )

        call ReplaceSubString( mstr, collapsedstr, '>', ' & ', &
          & which='all', no_trim=.true. )
        call ReplaceSubString( collapsedstr, mstr, '&', '>', &
          & which='all', no_trim=.false. )

        collapsedstr = lowerCase(mstr)
        ! Collapse every sub-formula nested within parentheses
        do level =1, MAXNESTINGS ! To prevent endlessly looping if ill-formed
          if ( index( collapsedstr, '(' ) < 1 ) exit
          call SplitNest ( collapsedstr, part1, part2, part3 )
          ! Now evaluate the part2
          if ( DeeBUG ) then
            print *, 'part1 ', part1
            print *, 'part2 ', part2
            print *, 'part3 ', part3
          endif
          if ( part2 == ' ' ) then
            ! This should never happen with well-formed formulas
            collapsedstr = part1
            cycle
          else
            np = evaluatePrimitive( trim(part2), &
              & a, b, c )
            write(vChar, '(i4)') np
          endif
          ! And substitute its value for the spaces it occupied
          if (  part1 == ' ' ) then
            collapsedstr = trim(vChar) // ' ' // part3
          elseif ( part3 == ' ' ) then
            collapsedstr = trim(part1) // ' ' // vChar
          else
            collapsedstr = trim(part1) // ' ' // trim(vChar) // &
              & ' ' // part3
          endif
          if ( DeeBUG ) then
            print *, 'collapsedstr ', collapsedstr
          endif
        enddo
        ! Presumably we have collapsed all the nested '(..)' pairs by now
        np = evaluatePrimitive( trim(collapsedstr), &
              & a, b, c )
        if ( DeeBUG ) then
          print *, 'np ', np
          print *, 'size(database) ', size(primitives)
        endif
        quantity%values = 0.
        if ( np < 1 .or. np > size(primitives) ) then
          print *, 'np ', np
          print *, 'size(database) ', size(primitives)
          call Announce_Error ( key, no_error_code, &
            & 'Illegal index for primitives array' )
          return
        endif
        if ( .not. associated ( quantity%mask ) ) then
          quantity%values = primitives(np)%values
        else
          where ( iand ( ichar(quantity%mask(:,:)), m_fill ) == 0 )
            quantity%values = primitives(np)%values
          end where
        end if
        if ( DeeBUG ) call dumpPrimitives(primitives)
        call destroyPrimitives(primitives)
      end subroutine SimpleExprWithC

      subroutine reorderPrecedence(mstr, collapsedstr)
        ! Identify all the terms where each term are separated by
        ! the lower-precedence operators {+, -,<,>}
        ! If any terms contain higher-precedence operators {*, /}
        ! then surround them by parentheses
        character(len=*), intent(in)  :: mstr
        character(len=*), intent(out) :: collapsedstr
        ! Internal variables
        logical, parameter :: COUNTEMPTY = .true.
        integer :: i
        integer :: n
        character(len=(len(mstr)+3)) :: element
        character(len=(len(mstr)+3)) :: temp
        ! Executable
        ! 1st -- replace each '-' with '+-'
        ! (Don't worry--we'll undo this before returning)
        ! (It takes two steps for each to avoid threat of infinite loop)
        call ReplaceSubString( mstr, collapsedstr, '-', '&', &
          & which='all', no_trim=.true. )
        call ReplaceSubString( collapsedstr, temp, '&', '+-', &
          & which='all', no_trim=.false. )

        call ReplaceSubString( temp, collapsedstr, '<', '&', &
          & which='all', no_trim=.true. )
        call ReplaceSubString( collapsedstr, temp, '&', '+<', &
          & which='all', no_trim=.false. )

        call ReplaceSubString( temp, collapsedstr, '>', '&', &
          & which='all', no_trim=.true. )
        call ReplaceSubString( collapsedstr, temp, '&', '+>', &
          & which='all', no_trim=.false. )
        ! Now loop over terms
        n = NumStringElements( temp, COUNTEMPTY, inseparator='+' )
        if ( n < 1 ) then
          call ReplaceSubString( temp, collapsedstr, '+-', '-', &
            & which='all', no_trim=.false. )
          call ReplaceSubString( collapsedstr, temp, '+<', '<', &
            & which='all', no_trim=.false. )
          call ReplaceSubString( temp, collapsedstr, '+>', '>', &
            & which='all', no_trim=.false. )
          return
        endif
        collapsedstr = ' '
        do i=1, n
          call GetStringElement ( temp, element, i, countEmpty, inseparator='+' )
          ! Surround term with parentheses if it's a product or quotient
          ! but not if it's not simple
          if ( ( index(element, '*') > 0 .or. index(element, '/') > 0 ) .and. &
            & index(element, ')')  < 1 ) then
            element = '(' // trim(element) // ')'
          endif
          collapsedstr = catLists( collapsedstr, element, inseparator='+' )
        enddo
        
        ! Now undo change by reverting all '+-'
        ! (including any that may have been split by parentheses
        call ReplaceSubString( collapsedstr, temp, '+-', '-', &
          & which='all', no_trim=.false. )
        call ReplaceSubString( temp, collapsedstr, '+(-', '-(', &
          & which='all', no_trim=.false. )

        call ReplaceSubString( collapsedstr, temp, '+<', '<', &
          & which='all', no_trim=.false. )
        call ReplaceSubString( temp, collapsedstr, '+(<', '<(', &
          & which='all', no_trim=.false. )

        call ReplaceSubString( collapsedstr, temp, '+>', '>', &
          & which='all', no_trim=.false. )
        call ReplaceSubString( temp, collapsedstr, '+(>', '>(', &
          & which='all', no_trim=.false. )
      end subroutine reorderPrecedence

      subroutine destroyPrimitives(primitives)
        ! deallocate all the arrays we created
        type(arrayTemp_T), dimension(:), pointer :: primitives
        integer :: i
        if ( .not. associated(primitives) ) return
        if ( size(primitives) < 1 ) return
        do i=1, size(primitives)
          call deallocate_test( primitives(i)%values, &
            & 'values', ModuleName // '/destroyPrimitives' )
        enddo
      end subroutine destroyPrimitives

      subroutine dumpAPrimitive(primitive)
        ! dump all the values in the array
        type(arrayTemp_T), intent(in) :: primitive
        if ( .not. associated(primitive%values) ) then
          call output( 'values not associated ', advance='yes' )
          return
        endif
        if ( size(primitive%values) < 1 ) then
          call output( 'values array is of 0 size ', advance='yes' )
          return
        endif
        call dump( primitive%values, 'values' )
      end subroutine dumpAPrimitive

      subroutine dumpPrimitives(primitives)
        ! dump all the arrays we created
        type(arrayTemp_T), dimension(:), pointer :: primitives
        integer :: i
        if ( .not. associated(primitives) ) then
          call output( 'database not associated ', advance='yes' )
          return
        endif
        if ( size(primitives) < 1 ) then
          call output( 'empty database ', advance='yes' )
          return
        endif
        call output( 'size of primitives database: ', advance='no' )
        call output( size(primitives), advance='yes' )
        do i=1, size(primitives)
          call output ( ' index of primitive: ', advance='no' )
          call output ( i, advance='yes' )
          call dumpAPrimitive( primitives(i) )
        enddo
      end subroutine dumpPrimitives

      integer function AddPrimitiveToDatabase( DATABASE, ITEM )

        ! This function adds a primitive data type to a database of said types,
        ! creating a new database if it doesn't exist.  The result value is
        ! the size -- where it was put.

        ! Dummy arguments
        type (arrayTemp_T), dimension(:), pointer :: DATABASE
        type (arrayTemp_T), intent(in) :: ITEM

        ! Local variables
        type (arrayTemp_T), dimension(:), pointer :: tempDatabase
        !This include causes real trouble if you are compiling in a different
        !directory.
        include "addItemToDatabase.f9h" 

        AddPrimitiveToDatabase = newSize
      end function AddPrimitiveToDatabase

      function evaluatePrimitive( str, a, b, c ) result(value)
        ! Evaluate an expression composed entirely of
        ! (0) constants ('c')
        ! (1) primitives (e.g., '2')
        ! (2) unary operators ('-')
        ! (3) binary operators {'+', '-', '*', '/','<','>'}
        ! Dummy args
        character(len=*)                :: str
        integer                         :: value
        type (VectorValue_T), pointer   :: A
        type (VectorValue_T), pointer   :: B
        real(rv) :: C                     ! constant "c" in manipulation
        ! Internal variables
        logical                         :: done
        integer                         :: elem
        logical                         :: hit
        character(len=3)                :: lastOp ! {'+', '-', '*', '/'}
        integer                         :: n
        logical                         :: negating
        integer                         :: partID
        integer, dimension(2)           :: shp
        character(len=32)               :: variable
        type (arrayTemp_T)              :: newone
        type (arrayTemp_T)              :: part
        logical, parameter              :: DEEBUG = .false.
        ! Executable
        shp = shape(a%values)
        call allocate_test( newone%values, shp(1), shp(2), &
          & 'newone', ModuleName // '/evaluatePrimitive' )
        call allocate_test( part%values, shp(1), shp(2), &
          & 'part', ModuleName // '/evaluatePrimitive' )

        if ( deeBug ) then
          print *, 'Complete dump of database'
          call dumpPrimitives(primitives)
        endif

        done = .false.
        negating = .false.
        elem = 0
        lastOp = 'nul' ! 'or'
        newone%values = 0.
        n = NumStringElements( trim(str), countEmpty=.false., &
          & inseparator=' ' )
          if ( DeeBUG ) then
            print *, n, ' str: ', trim(str)
          endif
        do
          ! go through the elements, re-evaluating every time we "hit" a primitive
          ! Otherwise revising our lastOp or negating status
          elem = elem + 1
          call GetStringElement ( trim(str), variable, elem, &
            & countEmpty=.false., inseparator=' ' )
          if ( DeeBUG ) then
            print *, elem, ' variable: ', trim(variable)
          endif
          select case( trim(variable) )
          case ('a')
            partID = -1
            part%values = a%values
            hit = .true.
          case ('b')
            partID = -2
            part%values = b%values
            hit = .true.
          case ('c')
            partID = -3
            part%values = c
            hit = .true.
          case ('+')
            lastOp = '+'
            hit = .false.
          case ('*')
            lastOp = '*'
            hit = .false.
          case ('/')
            lastOp = '/'
            hit = .false.
          case ('-') ! could be unary or binary; how do we tell?
            if ( hit ) then ! already have a primitive; looking for an op
              lastOp = '-'
              hit = .false.
            else
              ! case ('not', '~')
              negating = .true.
              hit = .false.
            endif
          case ('<')
            lastOp = '<'
            hit = .false.
          case ('>')
            lastOp = '>'
            hit = .false.
          case default
            read( variable, * ) partID
            if ( partID < 1 ) then
              print *, 'partID: ', partID
              call Announce_Error ( key, no_error_code, 'partID too small' )
              return
            elseif( partID > size(primitives) ) then
              print *, 'partID: ', partID
              call Announce_Error ( key, no_error_code, 'partID too big' )
              return
            endif
            part%values = primitives(partID)%values
            hit = .true.
            if ( deeBug ) then
              print *, 'part"s values after ' // trim(lastOp) // trim(variable)
              call dumpAPrimitive(part)
              print *, 'based on'
              call dumpAPrimitive(primitives(partID))
            endif
          end select
          if ( hit ) then
            if ( negating ) part%values = -part%values
            select case(lastOp)
            case ('nul')
                newone%values = part%values
            case ('+')
                newone%values = newone%values + part%values
            case ('-')
                newone%values = newone%values - part%values
            case ('*')
                newone%values = newone%values * part%values
            case ('/')
              where ( part%values /= 0._rv )
                newone%values = newone%values / part%values
              end where
            case ('<')
                newone%values = min( newone%values, part%values )
            case ('>')
                newone%values = max( newone%values, part%values )
            case default
              ! How could this happen?
                call MLSMessage( MLSMSG_Error, ModuleName, &
                  & lastOp // ' not a legal binary op in evaluatePrimitive' )
            end select
            negating = .false.
            if ( deeBug ) then
              print *, 'newone"s values after ' // trim(lastOp) // trim(variable)
              call dumpAPrimitive(newone)
            endif
          endif
          if ( DeeBUG ) then
            print *, 'variable ', variable
            print *, 'partID ', partID
            print *, 'hit ', hit
            print *, 'negating ', negating
            print *, 'lastOp ', lastOp
          endif
          done = ( elem >= n )
          if ( done ) exit
        enddo
        value = AddPrimitiveToDatabase( primitives, newone )
!         call deallocate_test( newone%values, &
!           & 'newone', ModuleName // '/evaluatePrimitive' )
        call deallocate_test( part%values, &
          & 'part', ModuleName // '/evaluatePrimitive' )
        if ( .not. DEEBUG ) return
        print *, 'value ', value
        print *, 'newone"s values ' // trim(str)
        call dumpAPrimitive(newone)
        
        print *, 'values stored in db '
        call dumpAPrimitive(primitives(value))
      end function evaluatePrimitive

    ! ----------------------------------------- FillWithReflectorTemperature ---
    subroutine FillWithReflectorTemperature ( key, quantity, phiZero, termsNode )
      use Units, only: DEG2RAD
      integer, intent(in) :: KEY         ! Tree node for messages
      type (VectorValue_T), intent(inout) :: QUANTITY ! The quantity to fill
      real(r8), intent(in) :: PHIZERO   ! Offset term
      integer, intent(in) :: TERMSNODE

      ! Local variables
      integer :: I                      ! Loop counter
      integer, DIMENSION(2) :: UNITASARRAY ! Unit for value given
      real (r8), DIMENSION(2) :: VALUEASARRAY ! Value give

      ! Executable code
      if ( quantity%template%quantityType /= l_reflTemp ) &
        & call Announce_Error ( key, no_error_code, &
        & 'Inappropriate quantity for reflector temperature fill' )

      ! Loop over fourier terms, coefficients are son i+2
      do i = 0, nsons ( termsNode ) - 2
        call expr_check ( subtree(i+2,termsNode), unitAsArray, valueAsArray, &
          & (/PHYQ_Temperature/), unitsError )
        if ( unitsError ) call Announce_error ( termsNode, wrongUnits, &
          & extraInfo=(/unitAsArray(1), PHYQ_Temperature/) )
        ! Add in this coefficient
        if ( i == 0 ) then
          quantity%values = valueAsArray(1)
        else if ( mod ( i, 2 ) == 1 ) then
          ! A sine term
          quantity%values(1,:) = quantity%values(1,:) + 2 * valueAsArray(1) * &
            & cos ( Deg2Rad * ((i+1)/2) * ( quantity%template%phi(1,:) - phiZero ) )
        else
          ! A cosine term
          quantity%values(1,:) = quantity%values(1,:) - 2 * valueAsArray(1) * &
            & sin ( Deg2Rad * ((i+1)/2) * ( quantity%template%phi(1,:) - phiZero ) )
        end if
      end do

    end subroutine FillWithReflectorTemperature

    ! ----------------------------------------- FillQtyWithReichlerWMOTP -------------
    subroutine FillQtyWithReichlerWMOTP ( tpPres, temperature, refGPH )
      use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
      use dump_0, only: dump
      use MLSFillValues, only: IsFillValue, RemoveFillValues
      use output_m, only: output
      
      use wmoTropopause, only: ExtraTropics, twmo
      ! Implements the algorithm published in GRL
      ! Loosely called the "Reichler" algorithm
      ! Ideas the same as in FillQtyWithWMOTropopause
      ! But implemented differently
      ! 
      type (VectorValue_T), intent(inout) :: TPPRES ! Result
      type (VectorValue_T), intent(in) :: TEMPERATURE
      type (VectorValue_T), intent(in) :: REFGPH
      ! Local variables
      integer :: instance
      integer :: invert
      real :: MISSINGVALUE
      integer :: nlev
      integer :: nvalid
      real, dimension( size(temperature%values, 1) ) :: p ! Pa
      real, parameter :: pliml = 65.*100 ! in Pa
      real, parameter :: plimlex = 65.*100 ! in Pa
      real, parameter :: plimu = 550.*100 ! in Pa
      real, dimension( size(temperature%values, 1) ) :: t
      real :: trp
      real, dimension(:), pointer :: xyTemp, xyPress
      logical :: alreadyDumped
      logical, parameter :: DEEBUG = .false.
      ! Executable
      nlev = size(temperature%values, 1)
      MISSINGVALUE = REAL( DEFAULTUNDEFINEDVALUE )
      ! Loop over the instances
      tpPres%values = MISSINGVALUE
      if ( DEEBUG ) then
      call output( 'num levels', advance='no' )
      call output( nlev, advance='yes' )
      call output( 'num instances', advance='no' )
      call output( temperature%template%noInstances, advance='yes' )
      call output( 'shape of tpPress values', advance='no' )
      call output( shape(tpPres%values), advance='yes' )
      endif
      if ( -temperature%template%surfs(1,1) .gt.&
        & -temperature%template%surfs(2,1) ) then
        invert=1
        ! p = refGPH%values(nlev:1:-1,instance)*100.  ! hPa > Pa
        p = 10.0**(-temperature%template%surfs(nlev:1:-1,1))
        p = p * 100.  ! hPa > Pa
      else
        invert=0
        ! p=refGPH%values(:,instance)*100.         ! hPa > Pa
        p = 10.0**(-temperature%template%surfs(:,1)) 
        p = p * 100.  ! hPa > Pa
      end if
      instanceLoop: do instance = 1, temperature%template%noInstances
      alreadyDumped = .false.
        ! Check the temperatures are in a sensible range
        if ( any ( &
          & temperature%values(:,instance) < 10.0 .or. &
          & temperature%values(:,instance) > 1000.0 ) ) then
          call output('invalid temperatures in this instance', advance='yes' )
          cycle instanceLoop
        end if
        ! check vertical orientation of data
        ! (twmo expects ordered from top downward)
        ! but surfs are -1.0 * log10(press)
        if ( invert == 1 ) then
          t = temperature%values(nlev:1:-1,instance)
        else
          t = temperature%values(:,instance)
        end if
        where ( t < 0. .or. t > 100000. )
          t = MissingValue
        end where
        nvalid = count( .not. isFillValue(t) )
        if ( nvalid < 2 ) then
          if ( DEEBUG ) then
          call output('not enough valid temperatures in this instance', advance='yes' )
          call output(count( .not. isFillValue(t) ), advance='yes' )
          call output(count( isFillValue(t) ), advance='yes' )
          call dump(temperature%values(:,instance), 'temperature%values')
          endif
          cycle
        endif
        call Allocate_test (xyTemp, nvalid, 'xyTemp', ModuleName )
        call Allocate_test (xyPress, nvalid, 'xyPress', ModuleName )
        call RemoveFillValues( t, MISSINGVALUE, xyTemp, &
          & p, xyPress )
        call twmo(nvalid, xyTemp, xyPress, plimu, pliml, trp)
        if ( .not. alreadyDumped .and. DEEBUG ) then
           alreadyDumped = .true.
           call dump(p, 'p ')
           call dump(t, 't ')
           call dump(xyTemp, 'xyTemp')
           call dump(xyPress, 'xyPress')
           call output( 'plimu, pliml, trp', advance='no' )
           call output( (/ plimu, pliml, trp /), advance='yes' )
        endif
        ! Don't let tropopause sink too low in "extra tropics"
        if ( trp < plimlex .and. &
          & extraTropics(temperature%template%geodLat(1, instance) ) ) then
          if ( DEEBUG ) then
          call output( 'trp too low in extra tropics', advance='no' )
          call output( (/ plimlex, trp /), advance='yes' )
          endif
          trp = MISSINGVALUE
        endif
        if ( trp > 0. .and. trp < 100000000. ) tpPres%values(1, instance) = trp/100
        call Deallocate_test ( xyTemp, 'xyTemp', ModuleName )
        call Deallocate_test ( xyPress, 'xyPress', ModuleName )
      end do instanceLoop
      if ( DEEBUG ) call dump( tpPres%values, 'tpPres%values' )
    end subroutine FillQtyWithReichlerWMOTP

    ! ----------------------------------------- FillQtyWithWMOTropopause ------
    subroutine FillQtyWithWMOTropopause ( tpPres, temperature, refGPH, grid )
      use Geometry, only: GEODTOGEOCLAT
      use Hydrostatic_M, only: HYDROSTATIC
      use MLSNumerics, only: HUNT
      use VGridsDatabase, only: VGRID_T

      type (VectorValue_T), intent(inout) :: TPPRES ! Result
      type (VectorValue_T), intent(in) :: TEMPERATURE
      type (VectorValue_T), intent(in) :: REFGPH
      type (VGrid_T), intent(in) :: GRID

      ! Local variables
      real(rv), dimension(grid%noSurfs) :: TFINE ! Temperature on fine grid
      real(rv), dimension(grid%noSurfs) :: HFINE ! Temperature on fine grid
      real(rv), dimension(grid%noSurfs) :: DTDH ! d(tFine)/d(hFine)
      real(rv), dimension(grid%noSurfs) :: DUMMYT ! Extra arg. for interpolateValues
      integer :: I                      ! Instance counter
      integer :: S500                   ! 500mb surface?
      integer :: LOWCANDIDATE           ! Possbile tb but below 500hPa
      integer :: S                      ! Surface counter
      integer :: S0                     ! Surface to start looking from
      integer :: DS                     ! No. surfs less than 2km above s
      logical :: VALIDTP                ! Flag

      ! Now the rules for the WMO tropopause are:

      ! 1. The "first tropopause" is defined as the lowest height at which the
      ! lapse rate decreases to 2K per kilometer or less, provided also that
      ! the average lapse rate between this height and all higher altitudes
      ! within 2 kilometers does not exceed 2K per kilometer.

      ! 2. If, above the first tropopause, the average lapse rate between any
      ! height and all higher altitudes within a 1-kilometer intend exceeds 3K
      ! per kilometer, then another tropopause is defined by the same criteria
      ! as under 1 above. This second tropopause may be within or above the
      ! 1-kilometer layer.

      ! There are also two qualifying remarks attached to the selection
      ! criteria. They are as follows:

      ! 1. A height below the 500-mb level is not designated as a tropopause
      ! unless the sounding reaches the 200-mb level and the height is the only
      ! height satisfying the above definitions.

      ! 2. When the second or higher tropopauses are being determined, the
      ! 1-kilometer interval with an average lapse rate of 3K per kilometer can
      ! occur at any height above the conventional tropopause and not only at a
      ! height more than 2 kilometers above the first tropopause.

      ! Stricly melding the WMO tropopause definition with our own
      ! 'linear interpolation between tie points' mantra is problematic.
      ! If we are strict we'd only ever give the one of our basis points as
      ! the solution.  So we're going to interpolate to a higher resolution
      ! using a spline which pop's out the derivatives too.

      ! Note here we assume that the quantities supplied are all valid
      ! (ie describe the right quantities, and all on the same horizontal
      ! grid).

      ! Loop over the instances
      tpPres%values = 0.0
      instanceLoop: do i = 1, temperature%template%noInstances
        ! Check the temperatures are in a sensible range
        if ( any ( &
          & temperature%values(:,i) < 10.0 .or. &
          & temperature%values(:,i) > 1000.0 ) ) then
          cycle instanceLoop
        end if
        ! Perhaps think later about keeping track of the 'top'.
        ! Currently we're going to assert that we have data extending to
        ! pressures lower than 200mb.

        ! Interpolate temperature onto 'finer' pressure grid
        call InterpolateValues ( &
          & temperature%template%surfs(:,1), temperature%values(:,i), &
          & grid%surfs(:,1), tFine, method='Spline' )
        ! Now get the height for this.
        call Hydrostatic ( GeodToGeocLat ( temperature%template%geodLat(1,i) ), &
          & grid%surfs(:,1), tFine, grid%surfs(:,1), &
          & refGPH%template%surfs(1,1), refGPH%values(1,i), hFine )
        ! Note while much of the software thinks in meters, the hydrostatic
        ! routine works in km.
        ! Now do another spline 'interpolation' to get dTdH
        call InterpolateValues ( hFine, tFine, hFine, dummyT, &
          & dYbyDx=dtdh, method='Spline' )
        ! Negate dTdH to turn it into lapse rate
        dTdH = - dTdH

        ! Locate the 500mb surface
        call Hunt ( grid%surfs(:,1), -log10(500.0_r8), s500 )

        ! Find the 'first' tropopause.  We're not going to bother
        ! with the 'second'
        lowCandidate = 0
        validTP = .false.
        s0 = 1
        tpHunt: do
          ! Find first place where lapse rate less than 2K/km
          s = FindFirst ( dtdh(s0:) < 2.0 ) + s0 - 1
          ! If not found such a place give up looking
          if ( s == s0 - 1 ) exit tpHunt
          ! Find last surface within 2km of this one
          ds = FindFirst ( hFine(s:) > hFine(s) + 2.0 ) - 1
          ! If not got data 2km above s0 (FindFirst gave 0, so ds=-1), give up
          ! Also, if data on such coarse resolution that the next point above
          ! s0 is more than 2km away (FindFirst gave 1 so ds=0), also give up
          if ( ds <= 0 ) exit tpHunt
          ! Must have mean lapse rate less than 2K/km within 2km of s.
          validTp = sum ( dTdH ( s : s + ds - 1 ) / ds ) < 2.0
          ! If this is below 500mb, keep an eye on it and keep looking
          if ( s < s500 .and. validTP ) then
            validTP = .false.
            lowCandidate = s
          end if
          if ( validTP ) exit tpHunt
          ! Otherwise keep looking higher up
          s0 = s + 1
          if ( s0 > grid%noSurfs ) exit tpHunt
        end do tpHunt

        ! Now pick up the pieces of that complex logic
        if ( .not. validTP .and. lowCandidate > 0 ) then
          s = lowCandidate
          validTP = .true.
        end if

        ! If we never found one, leave 0 in the result and move onto
        ! the next instance
        if ( .not. validTP ) cycle instanceLoop

        ! OK, our tropopause is below surface s
        ! Now we do some interpolation to get
        ! the value we really want.
        if ( s > 1 ) then
          tpPres%values(1,i) = 10.0 ** ( -( &
            & grid%surfs(s-1,1) + &
            & ( grid%surfs(s,1) - grid%surfs(s-1,1) ) * ( dTdH(s-1) - 2.0 ) / &
            & ( dTdH(s-1) - dTdH(s) ) ) )
        else
          tpPres%values(1,i) = 10.0 ** ( - grid%surfs(s,1) )
        end if
      end do instanceLoop
    end subroutine FillQtyWithWMOTropopause

    ! -------------------------------------------- FillWithBinResults -----
    subroutine FillWithBinResults ( key, quantity, sourceQuantity, ptanQuantity, &
      & channel, method, additional, excludeBelowBottom, centerVertically )
      ! This fills a coherent quantity with the max/min binned value of
      ! a typically incoherent one.  The bins are centered horizontally
      ! on the profiles in quantity, but vertically the bins run between one
      ! surface and the next one up.
      integer, intent(in) :: KEY        ! Tree node
      type (VectorValue_T), intent(inout) :: QUANTITY
      type (VectorValue_T), intent(in) :: SOURCEQUANTITY
      type (VectorValue_T), pointer :: PTANQUANTITY
      integer, intent(in) :: CHANNEL
      integer, intent(in) :: METHOD
      logical, intent(in) :: ADDITIONAL
      logical, intent(in) :: EXCLUDEBELOWBOTTOM
      logical, intent(in) :: CENTERVERTICALLY

      ! Local variables
      real(r8), dimension(:,:), pointer :: SOURCEHEIGHTS ! might be ptan.
      real(r8), dimension(:), pointer :: PHIBOUNDARIES
      real(r8) :: V                     ! A value
      integer, dimension(:,:), pointer :: SURFS ! Surface mapping source->quantity
      integer, dimension(:,:), pointer :: INSTS ! Instance mapping source->quantity
      integer :: QS, QI, SS, SI                   ! Loop counters
      integer :: MYCHANNEL              ! Channel or 1

      ! Executable code

      ! Check the output quantity
      if ( .not. ValidateVectorQuantity ( quantity, &
        & coherent=.true., stacked=.true. ) ) &
        & call Announce_Error ( key, no_error_code, &
        & 'Illegal quantity for bin min/max/total fill' )

      ! Also should have the condition:
      !  quantityType = (/ sourceQuantity%template%quantityType /)
      ! However, the code does not yet really support the ability to do that.

      ! Work out source vertical coordinate
      if ( associated ( ptanQuantity ) .and. sourceQuantity%template%minorFrame ) then
        nullify ( sourceHeights )
        call Allocate_test ( sourceHeights, sourceQuantity%template%nosurfs, &
          & sourceQuantity%template%noinstances, 'sourceHeights', ModuleName )
        sourceHeights = ptanQuantity%values
        if ( quantity%template%verticalCoordinate /= l_zeta ) &
          & call Announce_Error ( key, no_error_code, &
          & 'Vertical coordinate in quantity to fill is not zeta' )
      else
        sourceHeights => sourceQuantity%template%surfs
        if ( sourceQuantity%template%verticalCoordinate /= quantity%template%verticalCoordinate ) &
          & call Announce_Error ( key, no_error_code, &
          & 'Vertical coordinates in binned fill do not match' )
      end if

      ! Setup index arrays
      nullify ( surfs, insts )
      call Allocate_test ( surfs, sourceQuantity%template%nosurfs, &
        & sourceQuantity%template%noinstances, 'surfs', ModuleName )
      call Allocate_test ( insts, sourceQuantity%template%nosurfs, &
        & sourceQuantity%template%noinstances, 'insts', ModuleName )

      ! Work out the vertical mapping, a function of instance for
      ! incoherent quantities
      if ( sourceQuantity%template%coherent ) then
        call Hunt ( quantity%template%surfs(:,1), sourceHeights(:,1), surfs(:,1), &
          & allowTopValue=.true., allowBelowValue=excludeBelowBottom, nearest=centerVertically )
        surfs = spread ( surfs(:,1), 2, sourceQuantity%template%noInstances )
      else
        do si = 1, sourceQuantity%template%noInstances
          call Hunt ( quantity%template%surfs(:,1), sourceHeights(:,si), surfs(:,si), &
            & allowTopValue=.true., allowBelowValue=excludeBelowBottom, nearest=centerVertically )
        end do
      end if

!      call dump ( surfs, 'surfs' )

      ! Work out the horizontal mapping, a function of height for unstacked quantities.
      ! Bin to the center of the profiles rather than the edges.
      if ( quantity%template%noInstances > 1 ) then
        ! Setup and work out the phiBoundaries
        nullify ( phiBoundaries )
        call Allocate_test ( phiBoundaries, quantity%template%noInstances-1, &
          & 'phiBoundaries', ModuleName )
        phiBoundaries = 0.5 * ( &
          & quantity%template%phi(1,1:quantity%template%noInstances-1) + &
          & quantity%template%phi(1,2:quantity%template%noInstances) )
        if ( sourceQuantity%template%stacked ) then
          call Hunt ( phiBoundaries, sourceQuantity%template%phi(1,:), &
            & insts(1,:), allowTopValue=.true., allowBelowValue=.true. )
          insts = spread ( surfs(1,:), 1, sourceQuantity%template%noSurfs ) + 1
        else
          do ss = 1, sourceQuantity%template%noSurfs
            call Hunt ( phiBoundaries, sourceQuantity%template%phi(ss,:), &
              & insts(ss,:), allowTopValue=.true., allowBelowValue=.true. )
          end do
          insts = insts + 1
        end if
        call Deallocate_test ( phiBoundaries, 'phiBoundaries', ModuleName )
      else
        insts = 1
      end if

!      call dump ( insts, 'insts' )

      ! Modify the surfs index to account for the channel, so now it's more
      ! of a 'values' index
      if ( sourceQuantity%template%frequencyCoordinate /= l_none .and. &
        & channel == 0 ) then
        call Announce_Error ( key, no_error_code, &
          & 'Must supply channel for this bin max/min fill' )
        return
      end if
      myChannel = channel
      if ( channel == 0 ) myChannel = 1

      ! Now loop over the output quantity points and work out the information
      do qi = 1, quantity%template%noInstances
        do qs = 1, quantity%template%noSurfs
          if ( count ( surfs == qs .and. insts == qi ) > 0 ) then
            select case ( method )
            case  ( l_binMax )
              ! Compute the maximum value in the bin
              v = maxval ( pack ( sourceQuantity%values ( &
                & myChannel : sourceQuantity%template%instanceLen : &
                &   sourceQuantity%template%noChans, : ), &
                & surfs == qs .and. insts == qi ) )
              if ( additional ) then
                quantity%values(qs,qi) = max ( quantity%values(qs,qi), v )
              else
                quantity%values(qs,qi) = v
              end if
            case ( l_binMean )
              ! Compute the average in the bin, be careful about dividing by zero
              v = sum ( pack ( sourceQuantity%values ( &
                & myChannel : sourceQuantity%template%instanceLen : &
                &   sourceQuantity%template%noChans, : ), &
                & surfs == qs .and. insts == qi ) ) / &
                & max ( count ( surfs == qs .and. insts == qi ), 1 )
              if ( additional ) then
                quantity%values(qs,qi) = 0.5 * ( quantity%values(qs,qi) + v )
              else
                quantity%values(qs,qi) = v
              end if
            case ( l_binMin )
              ! Compute the minimum value in the bin
              v = minval ( pack ( sourceQuantity%values ( &
                & myChannel : sourceQuantity%template%instanceLen : &
                &   sourceQuantity%template%noChans, : ), &
                & surfs == qs .and. insts == qi ) )
              if ( additional ) then
                quantity%values(qs,qi) = min ( quantity%values(qs,qi), v )
              else
                quantity%values(qs,qi) = v
              end if
            case ( l_binTotal )
              ! Compute the total in the bin
              v = sum ( pack ( sourceQuantity%values ( &
                & myChannel : sourceQuantity%template%instanceLen : &
                &   sourceQuantity%template%noChans, : ), &
                & surfs == qs .and. insts == qi ) )
              if ( additional ) then
                quantity%values(qs,qi) = quantity%values(qs,qi) + v
              else
                quantity%values(qs,qi) = v
              end if
            end select
          else
            quantity%values(qs,qi) = 0.0
          end if
        end do
      end do

      ! Now tidy up
      call Deallocate_test ( surfs, 'surfs', ModuleName )
      call Deallocate_test ( insts, 'insts', ModuleName )
      if ( associated ( ptanQuantity ) .and. sourceQuantity%template%minorFrame ) &
        & call Deallocate_test ( sourceHeights, 'sourceHeights', ModuleName )
    end subroutine FillWithBinResults

    ! --------------------------------------------- FillWithBoxcarAvergage  ----
    subroutine FillWithBoxcarFunction ( key, quantity, sourceQuantity, width, method )
      integer, intent(in) :: KEY        ! Key for tree node
      type (VectorValue_T), intent(inout) :: QUANTITY
      type (VectorValue_T), intent(in), target :: SOURCEQUANTITY
      integer, intent(in) :: WIDTH
      integer, intent(in) :: METHOD     ! L_MEAN, L_MAX, L_MIN
      ! Local variables
      integer :: I, I1, I2              ! Instance indices
      integer :: HALFWIDTH
      real(r8), dimension(:,:), pointer :: OLDVALUES

      ! Executable code
      if ( quantity%template%name /= sourceQuantity%template%name ) then
        call Announce_Error ( key, no_error_code, 'Quantity and source quantity do not match' )
        return
      end if
      if ( width <= 1 .or. mod ( width, 2 ) == 0 ) then
        call Announce_Error ( key, no_error_code, 'width must be greater than 1 and odd' )
        return
      end if

      if ( associated ( quantity%values, sourceQuantity%values ) ) then
        nullify ( oldValues )
        call Allocate_test ( oldValues, quantity%template%instanceLen, quantity%template%noInstances, &
          & 'oldValues', ModuleName )
        oldValues = sourceQuantity%values
      else
        oldValues => sourceQuantity%values
      end if

      halfWidth = width/2
      select case ( method )
      case ( l_mean )
        do i = 1, quantity%template%noInstances
          i1 = max ( i - halfWidth, 1 )
          i2 = min ( i + halfWidth, quantity%template%noInstances )
          quantity%values(:,i) = sum ( oldValues(:,i1:i2), dim=2 ) / (i2-i1+1)
        end do
      case ( l_min )
        do i = 1, quantity%template%noInstances
          i1 = max ( i - halfWidth, 1 )
          i2 = min ( i + halfWidth, quantity%template%noInstances )
          quantity%values(:,i) = minval ( oldValues(:,i1:i2), dim=2 )
        end do
      case ( l_max )
        do i = 1, quantity%template%noInstances
          i1 = max ( i - halfWidth, 1 )
          i2 = min ( i + halfWidth, quantity%template%noInstances )
          quantity%values(:,i) = maxval ( oldValues(:,i1:i2), dim=2 )
        end do
      end select

      if ( associated ( quantity%values, sourceQuantity%values ) ) then
        call Deallocate_test ( oldValues, 'oldValues', ModuleName )
      end if

    end subroutine FillWithBoxcarFunction

    ! -------------------------------------------- FillStatusQuantity --------
    subroutine FillStatusQuantity ( key, quantity, sourceQuantity, statusValue, &
      & minValue, maxValue, heightNode, additional )
      integer, intent(in) :: KEY        ! Tree node
      type ( VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      type ( VectorValue_T), intent(in) :: SOURCEQUANTITY ! Chisq like quantity on which it's based
      integer, intent(in) :: STATUSVALUE
      real(r8), intent(in) :: MINVALUE     ! A scale factor
      real(r8), intent(in) :: MAXVALUE     ! A scale factor
      integer, intent(in) :: HEIGHTNODE ! What heights
      logical, intent(in) :: ADDITIONAL ! Is this an additional flag or a fresh start?
      ! Local variables
      integer, dimension(2) :: UNITASARRAY ! From expr
      real(r8), dimension(2) :: VALUEASARRAY ! From expr
      real(r8) :: HEIGHT                ! The height to consider
      integer :: SURFACE                ! Surface index
      ! Executable code
      ! Do some sanity checking
      if ( quantity%template%quantityType /= l_status ) call Announce_error ( key, no_error_code, &
        & 'Quality quantity must be quality' )
      if ( .not. DoHGridsMatch ( quantity, sourceQuantity ) ) call Announce_error ( &
        & key, no_error_code, 'quantity and sourceQuantity do not have matching hGrids' )

      ! Work out the height
      if ( heightNode /= 0 ) then
        if ( nsons ( heightNode ) /= 2 ) call Announce_Error ( key, no_error_code, &
          & 'Only one height can be supplied for status fill' )
        if ( sourceQuantity%template%verticalCoordinate /= l_zeta ) &
          & call Announce_Error ( key, no_error_code, 'Bad vertical coordinate for sourceQuantity' )
        call expr_check ( subtree(2,heightNode) , unitAsArray, valueAsArray, &
          & (/PHYQ_Pressure/), unitsError )
        if ( unitsError ) call Announce_error ( heightNode, wrongUnits, &
          & extraInfo=(/unitAsArray(1), PHYQ_Pressure/) )
        height = - log10 ( valueAsArray(1) )
        call Hunt ( sourceQuantity%template%surfs(:,1), height, surface, nearest=.true. )
      else
        surface = 1
      end if
      if ( .not. additional ) quantity%values = 0.0_r8
        ! quantity%values = iand ( nint ( quantity%values ), not ( statusValue ) )
      where ( sourceQuantity%values(surface,:) > maxValue .or. sourceQuantity%values(surface,:) < minValue )
        quantity%values(1,:) = ior ( nint ( quantity%values(1,:) ), statusValue )
      end where
    end subroutine FillStatusQuantity

    ! -------------------------------------------- FillQualityFromChisq --------
    subroutine FillQualityFromChisq ( key, quantity, sourceQuantity, scale, heightNode )
      integer, intent(in) :: KEY        ! Tree node
      type ( VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      type ( VectorValue_T), intent(in) :: SOURCEQUANTITY ! Chisq like quantity on which it's based
      real(r8), intent(in) :: SCALE     ! A scale factor
      integer, intent(in) :: HEIGHTNODE ! What heights
      ! Local variables
      integer, dimension(2) :: UNITASARRAY ! From expr
      real(r8), dimension(2) :: VALUEASARRAY ! From expr
      real(r8) :: HEIGHT                ! The height to consider
      integer :: SURFACE                ! Surface index
      ! Executable code
      ! Do some sanity checking
      if ( quantity%template%quantityType /= l_quality ) call Announce_error ( key, no_error_code, &
        & 'Quality quantity must be quality' )
      if ( sourceQuantity%template%quantityType /= l_chisqBinned ) call Announce_error ( &
        & key, no_error_code, 'sourceQuantity must be of type chisqBinned' )
      if ( .not. DoHGridsMatch ( quantity, sourceQuantity ) ) call Announce_error ( &
        & key, no_error_code, 'quantity and sourceQuantity do not have matching hGrids' )

      ! Work out the height
      if ( heightNode /= 0 ) then
        if ( nsons ( heightNode ) /= 2 ) call Announce_Error ( key, no_error_code, &
          & 'Only one height can be supplied for quality fill' )
        call expr_check ( subtree(2,heightNode) , unitAsArray, valueAsArray, &
          & (/PHYQ_Pressure/), unitsError )
        if ( unitsError ) call Announce_error ( heightNode, wrongUnits, &
          & extraInfo=(/unitAsArray(1), PHYQ_Pressure/) )
        height = - log10 ( valueAsArray(1) )
        call Hunt ( sourceQuantity%template%surfs(:,1), height, surface, nearest=.true. )
      else
        surface = 1
      end if
      quantity%values(1,:) = 0.0_r8
      where ( sourceQuantity%values(surface,:) /= 0.0_r8 )
        quantity%values(1,:) = scale / sourceQuantity%values(surface,:)
      end where
    end subroutine FillQualityFromChisq

    ! -------------------------------------------- FillConvergenceFromChisq --------
    subroutine FillConvergenceFromChisq ( key, quantity, sourceQuantity, scale )
      integer, intent(in) :: KEY        ! Tree node
      type ( VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      type ( VectorValue_T), intent(in) :: SOURCEQUANTITY ! dnwt_ChisqRatio quantity on which it's based
      real(r8), intent(in) :: SCALE     ! A scale factor
      ! Local variables
      ! Executable code
      ! Do some sanity checking
      if ( quantity%template%quantityType /= l_quality ) call Announce_error ( key, no_error_code, &
        & 'Convergence quantity must be quality' )
      if ( sourceQuantity%template%quantityType /= l_dnwt_chisqRatio ) call Announce_error ( &
        & key, no_error_code, 'sourceQuantity must be of type chisqRatio' )

      quantity%values(1,:) = scale * sourceQuantity%values(1,1)
    end subroutine FillConvergenceFromChisq

    ! ------------------------------------------ FillUsingLeastSquares -----
    subroutine FillUsingLeastSquares  ( key, Quantity, SourceQuantity, ptanQuantity, &
      & channel, method, scaleInstances, scaleRatio, scaleSurfs )
      ! This fills a coherent Quantity from a a typically incoherent
      ! SourceQuantity using a least-squares approximation to a first-order
      ! Taylor series.

      use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
      use HFTI_M, only: HFTI
      use MLSKinds, only: R8

      integer, intent(in) :: KEY        ! Tree node
      type (VectorValue_T), intent(inout) :: QUANTITY
      type (VectorValue_T), intent(in) :: SOURCEQUANTITY
      type (VectorValue_T), pointer :: PTANQUANTITY
      integer, intent(in) :: CHANNEL
      integer, intent(in) :: METHOD
      real(r8), intent(inout) :: ScaleInstances, ScaleRatio, ScaleSurfs

      ! Local variables
      real(r8), dimension(:,:), pointer :: INSTS ! Instance coordinates for source
      real(r8), dimension(:,:), pointer :: LSMatrix ! Least Squares problem matrix
      real(r8), dimension(:), pointer :: RHS
      real(r8), dimension(:,:), pointer :: SOURCEHEIGHTS ! might be ptan.
      real(r8), dimension(:,:), pointer :: SURFS ! Surface coordinates for source
      real(r8) :: W, WMAX               ! Weight if method == l_lsWeighted
      integer :: QS, QI, SS, SI         ! Loop counters
      integer :: KRANK                  ! Rank of least-squares solution
      integer :: MYCHANNEL              ! Channel or 1
      integer :: NRows
      integer :: NSourceQuant  ! Number of values in SourceQuantity

      ! Executable code

      ! Check the output quantity
      if ( .not. ValidateVectorQuantity ( quantity, &
        & coherent=.true., stacked=.true. ) ) &
        & call Announce_Error ( key, no_error_code, &
        & 'Quantity for least-squares fill must be stacked and coherent' )

      ! Also should have the condition:
      !  quantityType = (/ sourceQuantity%template%quantityType /)
      ! However, the code does not yet really support the ability to do that.

      ! Work out source vertical coordinate
      if ( associated ( ptanQuantity ) .and. sourceQuantity%template%minorFrame ) then
        nullify ( sourceHeights )
        call Allocate_test ( sourceHeights, sourceQuantity%template%nosurfs, &
          & sourceQuantity%template%noinstances, 'sourceHeights', ModuleName )
        sourceHeights = ptanQuantity%values
        if ( quantity%template%verticalCoordinate /= l_zeta ) &
          & call Announce_Error ( key, no_error_code, &
          & 'Vertical coordinate in quantity to fill is not zeta' )
      else
        sourceHeights => sourceQuantity%template%surfs
        if ( sourceQuantity%template%verticalCoordinate /= quantity%template%verticalCoordinate ) &
          & call Announce_Error ( key, no_error_code, &
          & 'Vertical coordinates in least-squares fill do not match' )
      end if

      ! Get the channel number
      if ( sourceQuantity%template%frequencyCoordinate /= l_none .and. &
        & channel == 0 ) then
        call Announce_Error ( key, no_error_code, &
          & 'Must supply channel for this least-squares fill' )
        return
      end if
      myChannel = max(channel,1)

      ! Setup coordinate arrays
      nullify ( surfs, insts )
      call Allocate_test ( surfs, sourceQuantity%template%nosurfs, &
        & sourceQuantity%template%noinstances, 'surfs', ModuleName )
      call Allocate_test ( insts, sourceQuantity%template%nosurfs, &
        & sourceQuantity%template%noinstances, 'insts', ModuleName )

      ! Get the vertical coordinates, a function of instance for
      ! incoherent quantities
      do ss = 1, sourceQuantity%template%noSurfs
        if ( sourceQuantity%template%coherent ) then
          surfs(ss,:) = sourceHeights(ss,1)
        else
          surfs(ss,:) = sourceHeights(ss,:)
        end if
      end do

!      call dump ( surfs, 'surfs' )

      ! Get the horizontal coordinates, a function of height for unstacked quantities.
      do si = 1, sourceQuantity%template%noInstances
        if ( sourceQuantity%template%coherent ) then
          insts(:,si) = sourceQuantity%template%phi(1,si)
        else
          insts(:,si) = sourceQuantity%template%phi(:,si)
        end if
      end do

!      call dump ( insts, 'insts' )

      ! Allocate the arrays
      nullify ( lsMatrix, RHS )
      if ( sourceQuantity%template%regular ) then
        nSourceQuant = sourceQuantity%template%noInstances * &
          &            sourceQuantity%template%noSurfs
      else
        nSourceQuant = sourceQuantity%template%noInstances * &
          &                sourceQuantity%template%instanceLen
      end if
      call allocate_test ( lsMatrix, nSourceQuant, 4, 'LSMatrix', moduleName )
      call allocate_test ( RHS, nSourceQuant, 'RHS', moduleName )

      if ( method == l_lsWeighted ) then
        if ( scaleInstances < 0.0 ) scaleInstances = &
          ( quantity%template%phi(1,quantity%template%noInstances) - &
            quantity%template%phi(1,1) ) / &
          ( quantity%template%noInstances - 1 )
        if ( scaleSurfs < 0.0 ) scaleSurfs = &
          ( quantity%template%surfs(quantity%template%noSurfs,1) - &
            quantity%template%surfs(1,1) ) / &
          ( quantity%template%noSurfs - 1 )
        scaleSurfs = scaleSurfs * scaleRatio
      end if
      ! Now loop over the output quantity points and set up the least-squares
      ! problem.
      do qi = 1, quantity%template%noInstances
        do qs = 1, quantity%template%noSurfs
          nRows = 0
          wMax = -1.0
          do si = 1, sourceQuantity%template%noInstances
            do ss = 1, sourceQuantity%template%noSurfs
              if ( method == l_lslocal ) then
                if ( surfs(ss,si) < quantity%template%surfs(max(qs-1,1),1) ) cycle
                if ( surfs(ss,si) > quantity%template%surfs(min(qs+1,quantity%template%noSurfs),1) ) cycle
                if ( insts(ss,si) < quantity%template%phi(1,max(qi-1,1)) ) cycle
                if ( insts(ss,si) > quantity%template%phi(1,min(qi+1,quantity%template%noInstances)) ) cycle
              end if
              nRows = nRows + 1
              lsMatrix(nRows,:3) = (/ 1.0_r8, &
                & insts(ss,si) - quantity%template%phi(1,qi), &
                & surfs(ss,si) - quantity%template%surfs(qs,1) /)
              lsMatrix(nRows,4) = 0.5 * lsMatrix(nRows,2) * lsMatrix(nRows,3)
              rhs(nRows) = sourceQuantity%values(myChannel+(ss-1)*sourceQuantity%template%noChans,si)
              if ( method == l_lsWeighted ) then
                w = exp( - (lsMatrix(nRows,2)/scaleInstances)**2 &
                       & - (lsMatrix(nRows,3)/scaleSurfs)**2 )
                if ( w > wMax ) then
                  wMax = w
                else if ( w < sqrt(epsilon(w))*wMax ) then
                  nRows = nRows - 1
                  cycle
                end if
                lsMatrix(nRows,:) = lsMatrix(nRows,:) * w
                rhs(nRows) = rhs(nRows) * w
              end if
            end do
          end do
          if ( nRows < 4 ) then
            quantity%values(qs,qi) = quantity%template%badValue
          else
            ! Solve the least squares problem
            call hfti ( lsMatrix(:nRows,:), rhs(:nRows), krank=krank )
            if ( krank < 4 ) then
              quantity%values(qs,qi) = quantity%template%badValue
            else
              quantity%values(qs,qi) = rhs(1)
            end if
          end if
        end do
      end do

      ! Now tidy up
      call deallocate_test ( surfs,    'surfs',    moduleName )
      call deallocate_test ( insts,    'insts',    moduleName )
      call deallocate_test ( lsMatrix, 'LSMatrix', moduleName )
      call deallocate_test ( RHS,      'RHS',      moduleName )
      if ( associated ( ptanQuantity ) .and. sourceQuantity%template%minorFrame ) &
        & call Deallocate_test ( sourceHeights, 'sourceHeights', ModuleName )

    end subroutine FillUsingLeastSquares

    ! ----------------------------------------- offsetradiancequantity -----
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
      where ( iand ( ichar(radianceQuantity%mask), m_linAlg ) /= 0 )
        quantity%values = quantity%values + amount
      end where
    end subroutine OffsetRadianceQuantity

    ! ---------------------------------------------- ResetUnusedRadiances --
    subroutine ResetUnusedRadiances ( quantity, amount )
      type (VectorValue_T), intent(inout) :: QUANTITY
      real (rv), intent(in) :: AMOUNT
      ! Executable code
      if ( .not. ValidateVectorQuantity ( quantity, &
        & quantityType=(/l_radiance/) ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Quantity for resetUnusedRadiances fill is not radiance' )
      where ( quantity%values > amount*0.9 )
        quantity%values = quantity%values - amount
      end where
    end subroutine ResetUnusedRadiances

    ! ---------------------------------------------- ScaleOverlaps -------
    subroutine ScaleOverlaps ( quantity, multiplierNode, dontMask )
      type (VectorValue_T), intent(inout) :: QUANTITY
      integer, intent(in) :: MULTIPLIERNODE ! Tree node for factors
      logical :: DONTMASK               ! Flag
      ! Local variables
      real (r8) :: exprValue(2)         ! Tree expression
      integer :: exprUnit(2)            ! Tree units
      integer :: i,j                    ! Loop counters / indices

      ! Executable code
      scaleLowerLoop: do i = 1, quantity%template%noInstancesLowerOverlap
        if ( i+1 > nsons ( multiplierNode ) ) exit scaleLowerLoop
        call expr_check ( subtree( i+1, multiplierNode ), exprUnit, exprValue, &
          & (/PHYQ_Dimensionless/), unitsError )
        if ( unitsError ) then
          call Announce_error ( multiplierNode, wrongUnits, &
            & extraMessage="for scaleOverlaps fill", &
            & extraInfo=(/exprUnit(1), PHYQ_Dimensionless/) )
          return
        end if
        if ( associated ( quantity%mask ) .and. .not. dontMask ) then
          where ( iand ( ichar(quantity%mask(:,i)), m_Fill ) == 0 )
            quantity%values ( :, i ) = quantity%values ( :, i ) * exprValue(1)
          end where
        else
          quantity%values ( :, i ) = quantity%values ( :, i ) * exprValue(1)
        end if
      end do scaleLowerLoop

      scaleUpperLoop: do i = quantity%template%noInstances, &
        & quantity%template%noInstances - &
        & quantity%template%noInstancesUpperOverlap + 1, - 1
        j = quantity%template%noInstances - i + 1
        if ( j+1 > nsons ( multiplierNode ) ) exit scaleUpperLoop
        call expr_check ( subtree( j+1, multiplierNode ), exprUnit, exprValue, &
          & (/PHYQ_Dimensionless/), unitsError )
        if ( unitsError ) then
          call Announce_error ( multiplierNode, wrongUnits, &
            & extraMessage="for scaleOverlaps fill", &
            & extraInfo=(/exprUnit(1), PHYQ_Dimensionless/) )
          return
        end if
        if ( associated ( quantity%mask ) .and. .not. dontMask ) then
          where ( iand ( ichar(quantity%mask(:,i)), m_Fill ) == 0 )
            quantity%values ( :, i ) = quantity%values ( :, i ) * exprValue(1)
          end where
        else
          quantity%values ( :, i ) = quantity%values ( :, i ) * exprValue(1)
        end if
      end do scaleUpperLoop
    end subroutine ScaleOverlaps

    ! ----------------------------------------- SpreadChannelFill --------
    subroutine SpreadChannelFill ( quantity, channel, dontMask, key, &
      & sourceQuantity )
      type(VectorValue_T), intent(inout) :: QUANTITY
      integer, intent(in) :: CHANNEL
      logical, intent(in) :: DONTMASK
      integer, intent(in) :: KEY
      type(VectorValue_T), intent(in), optional :: SOURCEQUANTITY
      ! Local variables
      integer :: I                      ! Instance loop counter
      integer :: C                      ! Channel loop counter
      integer :: S                      ! Surface loop counter
      integer :: J                      ! Destination index
      integer :: MYCHANNEL              ! Possibly offset channel
      type (Signal_T) ::  signal        ! Signal for this quantity

      ! Exectuable code
      if ( present(sourceQuantity) ) then
        do i = 1, quantity%template%noInstances
          do s = 1, quantity%template%noSurfs
            do c = 1, quantity%template%noChans
              j = (s-1)*quantity%template%noChans + c
              if ( associated ( quantity%mask ) .and. .not. dontMask ) then
                if ( iand ( ichar(quantity%mask(j,i)), m_Fill ) == 1 ) cycle
              end if
              quantity%values ( j, i ) = &
                & sourceQuantity%values ( s, i )
            end do
          end do
        end do
        return
      endif
      ! Deal with any channel numbering issues.
      signal = GetSignal ( quantity%template%signal )
      myChannel = channel - lbound ( signal%frequencies, 1 ) + 1
      if ( myChannel < 1 .or. myChannel > quantity%template%noChans ) &
        & call Announce_Error ( key, no_error_code, &
        & 'Invalid channel for spread channel' )
      do i = 1, quantity%template%noInstances
        do s = 1, quantity%template%noSurfs
          do c = 1, quantity%template%noChans
            if ( c == myChannel ) cycle
            j = (s-1)*quantity%template%noChans + c
            if ( associated ( quantity%mask ) .and. .not. dontMask ) then
              if ( iand ( ichar(quantity%mask(j,i)), m_Fill ) == 1 ) cycle
            end if
            quantity%values ( j, i ) = &
              & quantity%values ( (s-1)*quantity%template%noChans + myChannel, i )
          end do
        end do
      end do
    end subroutine SpreadChannelFill

    ! ---------------------------------------------- TransferVectors -----
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
            if ( associated(sq%mask) ) then
              if ( .not. associated(dq%mask)) call CreateMask ( dq )
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

      use Dump_0, only: Dump
      use Intrinsic, only: Field_indices, PHYQ_Indices
      use MoreTree, only: Get_Field_Id, StartErrorMessage
      use String_Table, only: Display_String

      integer, intent(in) :: where   ! Tree node where error was noticed
      integer, intent(in) :: CODE    ! Code for error message
      character (LEN=*), intent(in), optional :: ExtraMessage
      integer, intent(in), dimension(:), optional :: ExtraInfo

      integer :: I

      error = max(error,1)
      call StartErrorMessage ( where )
      if ( code  > no_Error_Code ) call output ( 'The' );

      select case ( code )
      case ( badEstNoiseFill )
        call output ( " estimated noise fill is missing information", advance='yes' )
      case ( badGeocAltitudeQuantity )
        call output ( " geocAltitudeQuantity is not geocAltitude", advance='yes' )
      case ( badlosGridfill )
        call output ( " information for los Grid fill is incomplete/incorrect", advance='yes' )
      case ( badlosvelfill )
        call output ( " information for los velocity is incomplete/incorrect", advance='yes' )
      case ( badIsotopeFill )
        call output ( " information for isotope fill is incomplete/incorrect", advance='yes' )
      case ( badREFGPHQuantity )
        call output ( " refGPHQuantity is not refGPH", advance='yes' )
      case ( badRefractFill )
        call output ( " phiTan refract fill is missing information", advance='yes' )
      case ( badTemperatureQuantity )
        call output ( " temperatureQuantity is not temperature", advance='yes' )
      case ( bothFractionAndLength )
        call output ( " config specifies both fraction and lengthScale", advance='yes' )
      case ( cantFillFromL1B )
        call output ( " command could not be filled from L1B.", advance='yes' )
      case ( cantFillFromL2AUX )
        call output ( " command could not be filled from L2AUX.", advance='yes' )
      case ( cantInterpolate3D )
        call output ( " program cannot interpolate 3d quantities (yet).", advance='yes' )
      case ( emptyGridForFill )
        call output ( " config specifies an empty grid for the fill", advance='yes' )
      case ( errorReadingL1B )
        call output ( " L1B file could not be read.", advance='yes' )
      case ( invalidExplicitFill )
        call output ( " value has inappropriate dimensionality for explicit fill.", advance='yes' )
        call output ( " Should have " )
        call output ( extraInfo )
        call output ( " elements.", advance='yes' )
      case ( missingDataInGrid )
        call output ( " grid for fill has missing/bad data points", advance='yes' )
      case ( missingField )
        call output ( " fields " )
        do i = 1, size(extraInfo)
          call display_string ( field_indices(extraInfo(i)) )
          if ( i == size(extraInfo) ) then
            call output ( " and " )
          else
            call output ( ", " )
          end if
        end do
        call output ( " are required.", advance='yes' )
      case ( needGeocAltitude )
        call output ( " fill needs geocAltitudeQuantity.", advance='yes' )
      case ( needH2O )
        call output ( " fill needs H2OQuantity.", advance='yes' )
      case ( needOrbitInclination )
        call output ( " fill needs OrbitalInclination.", advance='yes' )
      case ( needTempREFGPH )
        call output ( " needs temperatureQuantity and refGPHquantity.", advance='yes' )
      case ( no_Error_Code ) ! Handled at the bottom
      case ( noExplicitValuesGiven )
        call output ( " explicit fill requires explicit values.", advance='yes' )
      case ( nonConformingHydrostatic )
        call output ( " quantities needed for hydrostatic fill do not conform", advance='yes' )
      case ( noSourceGridGiven )
        call output ( " gridded fill requires a sourceGrid field.", advance='yes' )
      case ( noSourceL2AUXGiven )
        call output ( " L2AUX fill requires a sourceL2AUX field.", advance='yes' )
      case ( noSourceL2GPGiven )
        call output ( " L2GP fill requires a sourceL2GP field.", advance='yes' )
      case ( noSpecialFill )
        call output ( " special fill is invalid", advance='yes' )
      case ( notImplemented )
        call output ( extraMessage )
        call output ( " method is not implemented yet.", advance='yes' )
      case ( notSPD )
        call output ( " matrix is not a SPD matrix.", advance='yes' )
      case ( notPlain )
        call output ( " matrix is not a plain matrix.", advance='yes' )
      case ( NotZetaForGrid )
        call output ( " Quantity not on zeta surfaces.", advance='yes' )
      case ( wrongUnits )
        call display_string ( field_indices(Get_Field_Id(where)), &
          & before=" values of the " )
        call output ( " field have the wrong units", advance="yes" )
        if ( present(extraInfo) ) then
          i = 1
          if ( size(extraInfo) > 1 ) then
            call display_String ( phyq_indices(extraInfo(1)), &
              & before="Units are " )
            call output ( ", " )
            i = 2
          end if
          call display_String ( phyq_indices(extraInfo(i:)), advance='yes', &
            & before="Units should be" )
        end if
      case default
        call output_name_v_pair( ' error code', Code, advance='no' )
        call output ( " command caused an unrecognized programming error", advance='yes' )
      end select
      if ( present(ExtraMessage) )  call output(ExtraMessage, advance='yes')
      if ( code == no_Error_Code .and. present(extraInfo) ) &
        & call dump ( extraInfo, name='Extra info' )
    end subroutine ANNOUNCE_ERROR

  end subroutine MLSL2Fill

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Fill
!=============================================================================

!
! $Log$
! Revision 2.347  2006/11/03 00:26:07  pwagner
! Fixed bug in tropopause calculation
!
! Revision 2.346  2006/10/11 22:56:11  pwagner
! May fill convergence from dnwt_chisqRatio
!
! Revision 2.345  2006/10/03 20:24:11  pwagner
! Optional test, tweaks to FillChiSqRatio
!
! Revision 2.344  2006/10/02 23:05:03  pwagner
! May Fill chi^2 ratio to measure convergence
!
! Revision 2.343  2006/08/03 01:58:03  vsnyder
! Better error message if profile is out of order
!
! Revision 2.342  2006/08/02 19:52:29  vsnyder
! Move destroy processing to DestroyCommand_m
!
! Revision 2.341  2006/07/28 01:55:20  vsnyder
! Improve error detection and reporting for some allocations.
! Comment out matrix and vector dumps conditioned on the -g level since
! these dumps can be gotten by dump and dumpBlocks commands.  Inserted
! a comment to that effect.
!
! Revision 2.340  2006/07/27 23:07:38  pwagner
! Attempt to enforce conformity between column unit attributes and how we compute them
!
! Revision 2.339  2006/07/27 03:54:11  vsnyder
! Include source_ref in created vectors and matrices
!
! Revision 2.338  2006/07/12 20:41:16  pwagner
! Fixed BO size mismatch that only NAG caught
!
! Revision 2.337  2006/07/07 23:08:53  pwagner
! Fixed bug in filling from GEOS5-derived grid
!
! Revision 2.336  2006/06/13 22:14:18  pwagner
! Recover gracefully if l1boa file lacks BO_stat dataset
!
! Revision 2.335  2006/06/12 19:28:52  pwagner
! Fallback to climatology noted only if all of ncep, goes4/5 missing
!
! Revision 2.334  2006/06/08 17:29:27  dwu
! add option to allow sourceQuantity in spreadChannel
!
! Revision 2.333  2006/06/03 01:43:36  vsnyder
! Allow multiple fields and multiple vectors/matrices per field on destroy
!
! Revision 2.332  2006/05/22 21:56:00  pwagner
! Fixed bug besetting manipulation fills
!
! Revision 2.331  2006/05/19 00:00:13  pwagner
! Added min, max operators ('<', '>') to manipulation fills with c
!
! Revision 2.330  2006/05/03 22:18:26  pwagner
! Sets mask reading quantities missing from l1b file
!
! Revision 2.329  2006/03/23 03:06:35  vsnyder
! Use HFTI instead of Cholesky for FillUsingLeastSquares, for stability
!
! Revision 2.328  2006/03/23 00:38:57  vsnyder
! Make SolveLS a little more general, in case somebody else wants to use it
!
! Revision 2.327  2006/03/22 23:47:37  vsnyder
! More work on least-squares fill
!
! Revision 2.326  2006/03/15 23:55:11  pwagner
! Use reorderPrecedence to ensure higher precedes lower in manipulations
!
! Revision 2.325  2006/03/13 23:42:37  pwagner
! Added c=numeric type field to Fill via manipulation
!
! Revision 2.324  2006/03/09 16:25:01  pwagner
! Fixed a bug only NAG caught
!
! Revision 2.323  2006/03/08 21:30:32  pwagner
! Added new manipulations: a/b, 1/a, -a, abs(a), sign(a), log(a)
!
! Revision 2.322  2006/02/21 19:16:10  pwagner
! GetHashElement is now a generic
!
! Revision 2.321  2006/02/15 00:01:19  pwagner
! Shoule let you fill H2O from new RHIvmr quantity
!
! Revision 2.320  2006/02/10 21:18:48  pwagner
! Added code for wmoTropopause (unused); dumps may go to special dumpfile
!
! Revision 2.319  2006/01/11 17:04:32  pwagner
! May specify unit when filling column abundances
!
! Revision 2.318  2006/01/06 01:16:34  pwagner
! silent boolean field can silence selected phases
!
! Revision 2.317  2006/01/05 03:48:17  vsnyder
! Use Interp_Bilinear_2d_1d (as InterpolateValues)
!
! Revision 2.316  2006/01/05 00:04:45  vsnyder
! Implement refractive correction for PhiTan fill, correct some error
! messages and rephrase others to make sense.
!
! Revision 2.315  2005/12/22 21:05:20  vsnyder
! Use Hydrostatic_no_der (by way of generic)
!
! Revision 2.314  2005/11/17 20:12:51  pwagner
! Can now read BO_stat successfully
!
! Revision 2.313  2005/11/15 00:21:12  pwagner
! Removed space between Is and Precision
!
! Revision 2.312  2005/11/11 21:49:41  pwagner
! May set l1b precisions negative if avoidBrightObjects set
!
! Revision 2.311  2005/10/18 16:56:43  pwagner
! Negative RHIPrecision when either T or H2O Precisions are
!
! Revision 2.310  2005/09/21 23:21:58  pwagner
! Use of single arg options in ExplicitFillVectorQuantity replaces three
!
! Revision 2.309  2005/08/04 03:28:50  vsnyder
! Correct fill for goofy L1BMIF_TAI, lots of cannonball polishing
!
! Revision 2.308  2005/08/03 18:09:38  vsnyder
! Cannon ball polishing, scan averaging
!
! Revision 2.307  2005/07/21 23:42:31  pwagner
! Repaired bugs in Fill status for missingGMAO; extras for explicit fill
!
! Revision 2.306  2005/07/12 17:40:52  pwagner
! May fill status with condition that no gmaos found
!
! Revision 2.305  2005/06/21 23:56:24  livesey
! Added forgiveZeros handling to FillCovariance for efficiency.
!
! Revision 2.304  2005/06/03 02:05:29  vsnyder
! New copyright notice, move Id to not_used_here to avoid cascades,
! get VGrids from VGridsDatabase instead of passing as an argument.
!
! Revision 2.303  2005/05/31 18:11:45  pwagner
! Restored 2.300 revisions mistakenly omitted from 2.302
!
! Revision 2.302  2005/05/31 17:51:17  pwagner
! Began switch from passing file handles to passing MLSFiles
!
! Revision 2.301  2005/05/28 03:25:40  vsnyder
! Cannonball polishing
!
! Revision 2.300  2005/05/27 20:03:06  vsnyder
! Dissassociated -> zero size before dump
!
! Revision 2.299  2005/03/24 21:23:46  pwagner
! Removed buggy, unused FillColAbundance
!
! Revision 2.298  2005/03/12 00:50:27  pwagner
! May restart warnings counter at each phase
!
! Revision 2.297  2004/12/01 01:24:44  livesey
! Handles missing L1BMAFBaselines in the same manner as missing radiances.
!
! Revision 2.296  2004/11/30 01:41:58  livesey
! Make folded sideband fill cope with absence of one sideband (for R1A
! case).
!
! Revision 2.295  2004/11/29 21:53:33  livesey
! Bug fix for reading L1B data when no radiances files specified.
! Also, added (a+b)/2 manipulation and changed definition of quality to be
! 1.0/chisq.
!
! Revision 2.294  2004/11/24 22:51:37  livesey
! Bug fix in off line sideband folding
!
! Revision 2.293  2004/11/08 21:57:00  livesey
! Added handling of 'badRange' in ASCII fill
!
! Revision 2.292  2004/10/21 19:32:55  livesey
! Got the Fill from ASCII file working.
!
! Revision 2.291  2004/10/16 17:26:08  livesey
! Added stub of fill from ascii file
!
! Revision 2.290  2004/10/13 02:25:11  livesey
! Changes to fill from vGrid
!
! Revision 2.289  2004/09/28 22:26:46  livesey
! Bug fix in applyBaseline, had quadrature handled the wrong way round.
!
! Revision 2.288  2004/09/27 20:11:29  livesey
! Added stuff for reading and applying L1BMAFBaseline.  This includes new
! suffix argument to L1B reading.
!
! Revision 2.287  2004/09/25 00:16:31  livesey
! Removed 'key' argument in CombineChannels call
!
! Revision 2.286  2004/09/24 17:55:57  livesey
! Moved FillWithCombinedChannels into ManipulateVectorQuantitites
!
! Revision 2.285  2004/09/24 03:38:26  livesey
! Added optional mapping matrix output to combine channels fill
!
! Revision 2.284  2004/09/21 22:59:18  livesey
! Another change to the chi-squared to quality conversion.
!
! Revision 2.283  2004/09/21 19:17:23  livesey
! Changed the 1-tanh to an exp in converting from chi squared to quality.
!
! Revision 2.282  2004/09/16 23:55:04  livesey
! Stopped binned fill from being so fussy.
!
! Revision 2.281  2004/09/10 23:53:10  livesey
! Added centerVertically option for binmean/max/min fill
!
! Revision 2.280  2004/08/24 23:24:46  pwagner
! Asks ReadL1BData to pad, contract--partly tested
!
! Revision 2.279  2004/08/23 21:59:42  pwagner
! Disabled debugging dumps of section, phase timings
!
! Revision 2.278  2004/08/03 18:01:14  pwagner
! Gets DEFAULTUNDEFINEDVALUE from MLSCommon
!
! Revision 2.277  2004/07/30 00:17:22  livesey
! Changed some errors to warnings
!
! Revision 2.276  2004/07/22 20:39:14  cvuu
! Now can fill ForwardModel time, mean and std_dev
!
! Revision 2.275  2004/06/29 18:06:28  pwagner
! May fill phase, section timings
!
! Revision 2.274  2004/06/10 00:58:45  vsnyder
! Move FindFirst, FindNext from MLSCommon to MLSSets
!
! Revision 2.273  2004/05/28 00:57:49  vsnyder
! Move GetIndexFlagsFromList from MoreTree to Expr_m
!
! Revision 2.272  2004/05/19 20:38:04  vsnyder
! Remove unreferenced symbols, polish some cannonballs
!
! Revision 2.271  2004/05/19 19:16:09  vsnyder
! Move MLSChunk_t to Chunks_m
!
! Revision 2.270  2004/05/04 01:03:56  livesey
! Added excludeBelowBottom flag for binmax/binmin fill
!
! Revision 2.269  2004/05/01 04:04:36  vsnyder
! Use DumpCommand
!
! Revision 2.268  2004/04/28 00:30:58  livesey
! Added a|b option in manipulate fill.
!
! Revision 2.267  2004/04/19 21:04:02  livesey
! Bug fix in ExtractSingleChannel
!
! Revision 2.266  2004/04/16 00:49:03  livesey
! Added extractChannel fill
!
! Revision 2.265  2004/04/13 21:19:10  livesey
! Bug fix in negative precision flagging.
!
! Revision 2.264  2004/04/02 01:06:46  livesey
! Got the status filling working.
!
! Revision 2.263  2004/03/22 18:25:25  livesey
! Added CombineChannels fill (may actually replace this before too long).
!
! Revision 2.262  2004/03/18 17:41:31  livesey
! Bug fix in quality fill.
!
! Revision 2.261  2004/03/17 17:16:11  livesey
! New status fill and new manipulations.
!
! Revision 2.260  2004/03/10 22:19:51  livesey
! Added quality fill method
!
! Revision 2.259  2004/03/03 22:40:59  livesey
! Added a<b and a>b to manipulations
!
! Revision 2.258  2004/03/03 19:26:38  pwagner
! More printing if tropopause outside pressure grid
!
! Revision 2.257  2004/02/20 00:43:27  pwagner
! Clarified warning message when tropopause too big/little
!
! Revision 2.256  2004/02/19 23:59:00  pwagner
! Integrates column abundances using WReads method
!
! Revision 2.255  2004/02/17 14:08:27  livesey
! Added functionality to the box car and binning fills
!
! Revision 2.254  2004/02/06 01:01:40  livesey
! Added boxcar fill method
!
! Revision 2.253  2004/01/30 23:28:33  livesey
! Insist on loading a plain matrix
!
! Revision 2.252  2004/01/29 03:32:42  livesey
! Made FillCovariance (temporarily?) fill both sides of the digaonal (in
! any case was wrongly doing upper).
!
! Revision 2.251  2004/01/23 19:07:35  livesey
! Finished off the adoption stuff
!
! Revision 2.250  2004/01/23 05:47:38  livesey
! Added the adoption stuff
!
! Revision 2.249  2004/01/20 20:26:03  livesey
! Added the binMean fill
!
! Revision 2.248  2003/12/04 22:19:32  livesey
! Added ability to fill from l2gpPrecision field
!
! Revision 2.247  2003/11/25 21:54:47  livesey
! Made the column filling algorithm much less fussy.
!
! Revision 2.246  2003/11/05 18:37:25  pwagner
! Now can dump either entire vector or a single quantity
!
! Revision 2.245  2003/10/22 21:17:06  pwagner
! aPhaseName: Phase added to Fill, Construct sections to time phases
!
! Revision 2.244  2003/10/15 23:12:08  livesey
! Added ResetUnusedRadiances
!
! Revision 2.243  2003/10/07 15:44:27  cvuu
! add new flag ignoreGeolocation in subroutine FillVectorQuantityFromL2GP
!
! Revision 2.242  2003/09/25 16:41:12  michael
! magnetic field Elevation angle is constrained to 0-90 degrees.
!
! Revision 2.241  2003/09/09 22:06:38  livesey
! Added resilency to missing radiances.
!
! Revision 2.240  2003/08/28 00:44:54  livesey
! Made the a*b manipulation even more lenient
!
! Revision 2.239  2003/08/21 16:07:17  livesey
! Now calls FlushLockedBins (LinearizedForwardModel) rather than
! FlushL2PCBins (which is called by the former), as we want it also to
! forget which bin it chose.
!
! Revision 2.238  2003/08/20 20:05:42  livesey
! Added the a*b possibility to the manipulation fill.
!
! Revision 2.237  2003/08/16 00:29:37  vsnyder
! Correct a blunder: deg2rad should have been rad2deg
!
! Revision 2.236  2003/08/15 23:58:48  vsnyder
! Add MagAzEl fill method
!
! Revision 2.235  2003/08/13 19:23:45  vsnyder
! Make an error message more informative
!
! Revision 2.234  2003/08/08 23:07:02  livesey
! Added the rotate field fill.
!
! Revision 2.233  2003/07/16 22:39:47  livesey
! Bug fix in fill from l2aux
!
! Revision 2.232  2003/07/08 00:17:46  livesey
! Bug fix in column filling
!
! Revision 2.231  2003/06/24 19:59:42  livesey
! Fixed something in OffsetRadianceQuantity that might have become a bug
! one day (assumed m_linAlg=1).
!
! Revision 2.230  2003/06/20 19:37:06  pwagner
! Quanities now share grids stored separately in databses
!
! Revision 2.229  2003/06/05 22:08:55  livesey
! Cosmetic and superficial changes to FillFromSplitSideband
!
! Revision 2.228  2003/06/03 19:23:51  livesey
! Added flushL2PCBins
!
! Revision 2.227  2003/05/29 20:01:55  livesey
! Added reflector temperature model.
!
! Revision 2.226  2003/05/29 16:41:56  livesey
! Renamed sideband fraction
!
! Revision 2.225  2003/05/28 06:00:06  livesey
! Bug fix in profile fill where the 'latching' to the output surfaces was
! leading to redundancy in the interpolation
!
! Revision 2.224  2003/05/26 06:32:50  livesey
! Various mainly cosmetic changes to the column stuff
!
! Revision 2.223  2003/05/22 02:23:15  livesey
! Rewrite of explicit fill to make spread option more flexible.
!
! Revision 2.222  2003/05/22 00:26:33  dwu
! fix a problem in iwcfromextinction
!
! Revision 2.221  2003/05/21 18:58:57  dwu
! allow temperature and iwc on different hGrids in iwcFromExtinction
!
! Revision 2.220  2003/05/21 18:04:30  livesey
! Added a bit more intelligence to FillCovariance
!
! Revision 2.219  2003/05/20 23:10:24  dwu
! complete the addition of fill IWC from extinction
!
! Revision 2.218  2003/05/20 20:20:01  dwu
! add IWCfromExtinction
!
! Revision 2.217  2003/05/15 19:09:17  dwu
! changes in splitsideband
!
! Revision 2.216  2003/05/14 23:14:00  dwu
! nullify pointers in splitsideband
!
! Revision 2.215  2003/05/12 23:53:55  dwu
! fix a bug in splitsideband
!
! Revision 2.214  2003/05/12 22:11:07  dwu
! add more checkpoints in splitsideband
!
! Revision 2.213  2003/05/11 00:05:06  livesey
! Informative error message when L1B data wrong size
!
! Revision 2.212  2003/05/10 23:40:27  livesey
! Bug fixes in bining, other general tidyups.
!
! Revision 2.211  2003/05/10 22:20:12  livesey
! Made wmo tropopause resilient to being given stupid (i.e. 0) profiles.
!
! Revision 2.210  2003/05/10 01:07:58  livesey
! Added binTotal and noRads fills
!
! Revision 2.209  2003/05/08 19:34:58  dwu
! add more options to splitsideband, and tidy up
!
! Revision 2.208  2003/05/07 00:16:52  livesey
! Added the dump for magnetic field results.
!
! Revision 2.207  2003/05/06 21:00:23  livesey
! Bug fix on L1B stuff
!
! Revision 2.206  2003/04/30 22:07:14  pwagner
! Always sets errorCode to 0 in return from FillVectorQuantityFromL2AUX
!
! Revision 2.205  2003/04/24 22:17:02  dwu
! remove dump statement in fill binMinMax
!
! Revision 2.204  2003/04/24 00:35:14  dwu
! modify splitSideband to allow the splitted sideband cloud radiances being spread to other bands assuming the f**4 law
!
! Revision 2.203  2003/04/24 00:29:45  dwu
! modify splitSideband to allow the splitted sideband cloud radiances being spread to other bands assuming the f**4 law
!
! Revision 2.202  2003/04/23 17:06:36  livesey
! Added binmax binmin fills
!
! Revision 2.201  2003/04/11 23:15:09  livesey
! Added force option to vector fill, and spreadChannel fill method.
!
! Revision 2.200  2003/04/11 21:56:40  livesey
! Added wmo tropopause stuff
!
! Revision 2.199  2003/04/08 23:13:01  dwu
! an update on splitsideband
!
! Revision 2.198  2003/04/07 06:37:42  dwu
! implement splitsideband for cloud radiance
!
! Revision 2.197  2003/04/05 00:26:47  livesey
! Bug fix in sideband splitting stub
!
! Revision 2.196  2003/04/05 00:05:37  livesey
! Added call to getSignal in split sideband
!
! Revision 2.195  2003/04/04 23:53:57  livesey
! Added skeleton for split sideband fill
!
! Revision 2.194  2003/04/04 22:01:59  livesey
! Added call to updateMask
!
! Revision 2.193  2003/04/04 00:08:06  livesey
! Added the wrapping of the gridded data before it's used in fill.
!
! Revision 2.192  2003/03/27 20:45:02  livesey
! Added logSpace argument to profile fill, and made it obey the dontMask
! flag
!
! Revision 2.191  2003/03/26 21:23:47  livesey
! Added ScaleOverlaps stuff
!
! Revision 2.190  2003/03/19 19:22:24  pwagner
! Passes chunkNo around more widely
!
! Revision 2.189  2003/03/07 03:16:12  livesey
! Added RestrictRange
!
! Revision 2.188  2003/03/06 00:46:30  livesey
! Added ability to do subset and flagCloud
!
! Revision 2.187  2003/03/05 19:11:11  livesey
! Added allowMissing capability to gridded fill.
!
! Revision 2.186  2003/02/28 02:26:23  livesey
! Added checking for bad/missing data in fill from gridded data.
!
! Revision 2.185  2003/02/27 00:38:52  livesey
! Better handling of missing length scale in FillCovariance
!
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
