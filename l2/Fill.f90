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

  use MLSCOMMON, only: MLSFILE_T, DEFAULTUNDEFINEDVALUE
  use MLSKINDS, only: R8, RV
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

  logical, parameter :: USEREICHLER = .true.
contains ! =====     Public Procedures     =============================

  !---------------------------------------------------  MLSL2Fill  -----

  subroutine MLSL2FILL ( ROOT, FILEDATABASE, GRIDDEDDATABASE, VECTORTEMPLATES, &
    & VECTORS, QTYTEMPLATES, MATRICES, HESSIANS, L2GPDATABASE, L2AUXDATABASE, &
    & FWMODELCONFIG, CHUNKS, CHUNKNO )

    ! This is the main routine for the module.  It parses the relevant lines
    ! of the l2cf and works out what to do.

    use ALLOCATE_DEALLOCATE, only: TEST_ALLOCATE
    use CHUNKS_M, only: MLSCHUNK_T
    use DESTROYCOMMAND_M, only: DESTROYCOMMAND
    use DUMPCOMMAND_M, only: BOOLEANFROMANYGOODRADIANCES, &
      & BOOLEANFROMANYGOODVALUES, &
      & BOOLEANFROMCATCHWARNING, BOOLEANFROMCOMPARINGQTYS, BOOLEANFROMFORMULA, &
      & DUMPCOMMAND, INITIALIZEREPEAT, NEXTREPEAT, &
      & MLSCASE, MLSENDSELECT, MLSSELECT, MLSSELECTING, &
      & REPEAT=>SKIP, SKIP
    use EXPR_M, only: EXPR, EXPR_CHECK
    use FILLUTILS_1, only: ADDGAUSSIANNOISE, APPLYBASELINE, AUTOFILLVECTOR, &
      & COMPUTETOTALPOWER, DERIVATIVEOFSOURCE, FILLCOVARIANCE, &
      & EXTRACTSINGLECHANNEL, FILLERROR, FROMANOTHER, FROMGRID, &
      & FROML2GP, FROMPROFILE, GATHER, LOSVELOCITY, &
      & CHISQCHAN, CHISQMMAF, CHISQMMIF, CHISQRATIO, &
      & COLABUNDANCE, FOLDEDRADIANCE, PHITANWITHREFRACTION, &
      & IWCFROMEXTINCTION, RHIFROMORTOH2O, NORADSPERMIF, &
      & RHIPRECISIONFROMORTOH2O, WITHESTNOISE, &
      & HYDROSTATICALLY, FROMSPLITSIDEBAND, GPHPRECISION, &
      & FROMISOTOPE, FROMASCIIFILE, ROTATEMAGNETICFIELD, &
      & EXPLICIT, FROML1B, &
      & FROML2AUX, USINGMAGNETICMODEL, &
      & FROMINTERPOLATEDQTY, FROMLOSGRID, &
      & BYMANIPULATION, MANIPULATEVECTORS, WITHREFLECTORTEMPERATURE, &
      & WITHREICHLERWMOTP, &
      & WITHWMOTROPOPAUSE, WITHBINRESULTS, WITHBOXCARFUNCTION, &
      & STATUSQUANTITY, QUALITYFROMCHISQ, CONVERGENCEFROMCHISQ, &
      & USINGLEASTSQUARES, OFFSETRADIANCEQUANTITY, RESETUNUSEDRADIANCES, &
      & SCALEOVERLAPS, SCATTER, SPREADCHANNELFILL, TRANSFERVECTORS, &
      & TRANSFERVECTORSBYMETHOD, UNCOMPRESSRADIANCE, &
      & QTYFROMFILE, VECTORFROMFILE, ANNOUNCE_ERROR, &
      ! codes for ANNOUNCE_ERROR:
      & BADESTNOISEFILL, BADGEOCALTITUDEQUANTITY, BADISOTOPEFILL, &
      & BADLOSGRIDFILL, BADLOSVELFILL, BADREFGPHQUANTITY, &
      & BADREFRACTFILL, BADTEMPERATUREQUANTITY, BOTHFRACTIONANDLENGTH, &
      & MISSINGFIELD, &
      & NEEDGEOCALTITUDE, NEEDH2O, NEEDORBITINCLINATION, &
      & NEEDTEMPREFGPH, NOCODEFOR, NO_ERROR_CODE, NOEXPLICITVALUESGIVEN, &
      & NOSOURCEGRIDGIVEN, NOSOURCEL2AUXGIVEN, NOSOURCEL2GPGIVEN, &
      & NOTIMPLEMENTED, NOTPLAIN, NOTSPD, &
      & WRONGUNITS
    use FORWARDMODELCONFIG, only: FORWARDMODELCONFIG_T
    use FORWARDMODELSUPPORT, only: FILLFWDMODELTIMINGS
    use GLOBAL_SETTINGS, only: BRIGHTOBJECTS
    use GRIDDEDDATA, only: GRIDDEDDATA_T
    use HESSIANMODULE_1, only: ADDHESSIANTODATABASE, CREATEEMPTYHESSIAN, &
      & STREAMLINEHESSIAN, HESSIAN_T
    ! We need many things from init_tables_module.  first the fields:
    use INIT_TABLES_MODULE, only: F_A, F_ADDITIONAL, F_ALLOWMISSING, &
      & F_APRIORIPRECISION, F_ASPERCENTAGE, F_AUTOFILL, F_AVOIDBRIGHTOBJECTS, &
      & F_B, F_BADRANGE, F_BASELINEQUANTITY, F_BIN, F_BLOCK, &
      & F_BOUNDARYPRESSURE, F_BOXCARMETHOD, &
      & F_C, F_CENTERVERTICALLY, F_CHANNEL, F_CHANNELS, F_COLUMNS, F_COUNT, &
      & F_DESTINATION, F_DIAGONAL, F_DIMLIST, &
      & F_DONTMASK, &
      & F_ECRTOFOV, F_EARTHRADIUS, F_EXACT, F_EXCLUDEBELOWBOTTOM, &
      & F_EXPLICITVALUES, F_EXPR, F_EXTINCTION, &
      & F_FIELDECR, F_FILE, F_FLAGS, F_FORCE, F_SHAPE, &
      & F_FRACTION, F_FROMPRECISION, &
      & F_GEOCALTITUDEQUANTITY, F_Geolocation, F_GPHQUANTITY, &
      & F_HEIGHT, F_HDFVERSION, F_HEIGHTRANGE, F_HESSIAN, &
      & F_HIGHBOUND, F_H2OQUANTITY, F_H2OPRECISIONQUANTITY, &
      & F_IFMISSINGGMAO, &
      & F_IGNORENEGATIVE, F_IGNOREGEOLOCATION, F_IGNORETEMPLATE, F_IGNOREZERO, &
      & F_INSTANCES, F_INTEGRATIONTIME, F_INTERNALVGRID, &
      & F_INTERPOLATE, F_INVERT, F_INTRINSIC, F_ISPRECISION, &
      & F_LENGTHSCALE, F_LOGSPACE, F_LOSQTY, F_LOWBOUND, F_LSB, F_LSBFRACTION, &
      & F_MANIPULATION, F_MATRIX, F_MAXITERATIONS, F_MAXVALUE, F_MEASUREMENTS, &
      & F_METHOD, F_MINNORMQTY, F_MINVALUE, F_MODEL, F_MULTIPLIER, &
      & F_NOFINEGRID, F_NOISE, F_NOISEBANDWIDTH, F_NORMQTY, &
      & F_OFFSETAMOUNT, F_OPTIONS, F_ORBITINCLINATION, F_PHITAN, &
      & F_PHIWINDOW, F_PHIZERO, F_PRECISION, F_PRECISIONFACTOR, &
      & F_PROFILE, F_PROFILEVALUES, F_PTANQUANTITY, &
      & F_QUADRATURE, F_QUANTITY, F_QUANTITYNAMES, &
      & F_RADIANCEQUANTITY, F_RATIOQUANTITY, &
      & F_REFRACT, F_REFGPHQUANTITY, F_REFGPHPRECISIONQUANTITY, F_RESETSEED, &
      & F_RHIPRECISIONQUANTITY, F_RHIQUANTITY, F_ROWS, F_SCALE, &
      & F_SCALEINSTS, F_SCALERATIO, F_SCALESURFS, F_SCECI, &
      & F_SCVEL, F_SCVELECI, F_SCVELECR, F_SDNAME, F_SEED, F_SKIPMASK, &
      & F_SOURCE, F_SOURCEGRID, F_SOURCEL2AUX, F_SOURCEL2GP, &
      & F_SOURCEQUANTITY, F_SOURCEVGRID, F_SPREAD, F_START, F_STATUS, F_STRIDE, &
      & F_SUFFIX, F_SUPERDIAGONAL, F_SURFACE, &
      & F_SYSTEMTEMPERATURE, F_TEMPERATUREQUANTITY, F_TEMPPRECISIONQUANTITY, &
      & F_TEMPLATE, F_TNGTECI, F_TERMS, F_TOTALPOWERQUANTITY, &
      & F_TYPE, F_UNIT, F_USB, F_USBFRACTION, F_VECTOR, F_VMRQUANTITY, &
      & F_WHEREFILL, F_WHERENOTFILL, F_WIDTH, &
      & FIELD_FIRST, FIELD_LAST
    ! Now the literals:
    use INIT_TABLES_MODULE, only: L_ADDNOISE, L_APPLYBASELINE, L_ASCIIFILE, &
      & L_BINMAX, L_BINMEAN, L_BINMIN, L_BINTOTAL, &
      & L_BOUNDARYPRESSURE, L_BOXCAR, L_CHISQCHAN, &
      & L_CHISQMMAF, L_CHISQMMIF, L_CHISQRATIO, L_CHOLESKY, &
      & L_CLOUDICE, L_CLOUDEXTINCTION, &
      & L_COMBINECHANNELS, L_COLUMNABUNDANCE, L_CONVERGENCERATIO, &
      & L_DERIVATIVE, L_DOBSONUNITS, L_DU, &
      & L_ESTIMATEDNOISE, L_EXPLICIT, L_EXTRACTCHANNEL, L_FOLD, &
      & L_FWDMODELTIMING, L_FWDMODELMEAN, L_FWDMODELSTDDEV, &
      & L_GATHER, L_GEOLOCATION, &
      & L_GEOCALTITUDE, L_GEODALTITUDE, L_GPHPRECISION, L_GRIDDED, &
      & L_H2OFROMRHI, L_H2OPRECISIONFROMRHI, L_HYDROSTATIC, L_ISOTOPE, &
      & L_IWCFROMEXTINCTION, L_KRONECKER, &
      & L_L1B, L_L2GP, L_L2AUX, &
      & L_LOSVEL, L_LSGLOBAL, L_LSLOCAL, L_LSWEIGHTED, &
      & L_MAGAZEL, L_MAGNETICMODEL, &
      & L_MANIPULATE, L_MEAN, L_MODIFYTEMPLATE, L_MOLCM2, &
      & L_NEGATIVEPRECISION, L_NONE, &
      & L_NORADSPERMIF, L_OFFSETRADIANCE, &
      & L_PHASETIMING, L_PHITAN, &
      & L_PLAIN, L_PROFILE, L_PTAN,  L_QUALITY, &
      & L_RECTANGLEFROMLOS, L_REFGPH, L_REFRACT, &
      & L_REFLECTORTEMPMODEL, L_RESETUNUSEDRADIANCES, L_RHI, &
      & L_RHIFROMH2O, L_RHIPRECISIONFROMH2O, L_ROTATEFIELD, L_SCALEOVERLAPS, &
      & L_SECTIONTIMING, L_SCATTER, L_SPD, L_SPREADCHANNEL, &
      & L_SPLITSIDEBAND, L_STATUS, L_SWAPVALUES, &
      & L_TEMPERATURE, L_TNGTGEODALT, &
      & L_TNGTGEOCALT, L_UNCOMPRESSRADIANCE, L_VECTOR, L_VGRID, L_VMR, L_WMOTROPOPAuse, &
      & L_ZETA
    ! Now the specifications:
    use INIT_TABLES_MODULE, only: S_ANYGOODVALUES, S_ANYGOODRADIANCES, &
      & S_CASE, S_CATCHWARNING, S_COMPARE, S_COMPUTETOTALPOWER, S_DESTROY, &
      & S_DIFF, S_DIRECTREAD, S_DUMP, S_ENDSELECT, S_FILL, S_FILLCOVARIANCE, &
      & S_FILLDIAGONAL, S_FLAGCLOUD, S_FLUSHL2PCBINS, S_FLUSHPFA, S_HESSIAN, &
      & S_LOAD, S_MATRIX,  S_NEGATIVEPRECISION, S_PHASE, S_POPULATEL2PCBIN, &
      & S_REEVALUATE, S_REPEAT, S_RESTRICTRANGE, &
      & S_SELECT, S_SKIP, S_SNOOP, S_STREAMLINEHESSIAN, S_SUBSET, &
      & S_TIME, S_TRANSFER, S_UPDATEMASK, S_VECTOR
    ! Now some arrays
    use INTRINSIC, only: LIT_INDICES, &
      & PHYQ_DIMENSIONLESS, PHYQ_INVALID, PHYQ_TEMPERATURE, &
      & PHYQ_TIME, PHYQ_LENGTH, PHYQ_ANGLE, PHYQ_PROFILES
    use L2GPDATA, only: L2GPDATA_T, COL_SPECIES_HASH, COL_SPECIES_KEYS
    use L2AUXDATA, only: L2AUXDATA_T
    use L2PC_M, only: POPULATEL2PCBINBYNAME, LOADHESSIAN, LOADMATRIX, LOADVECTOR
    use L2PCBINS_M, only: FLUSHLOCKEDBINS
    use MANIPULATEVECTORQUANTITIES, only: DOHGRIDSMATCH, &
      & FILLWITHCOMBINEDCHANNELS
    use MATRIXMODULE_1, only: ADDTOMATRIXDATABASE, CREATEEMPTYMATRIX, &
      & GETACTUALMATRIXFROMDATABASE, GETDIAGONAL, &
      & GETKINDFROMMATRIXDATABASE, GETFROMMATRIXDATABASE, K_PLAIN, K_SPD, &
      & MATRIX_CHOLESKY_T, MATRIX_DATABASE_T, MATRIX_KRONECKER_T, MATRIX_SPD_T, &
      & MATRIX_T, NULLIFYMATRIX
    ! NOTE: if you ever want to include defined assignment for matrices, please
    ! carefully check out the code around the call to snoop.
    use MLSL2OPTIONS, only: DEFAULT_HDFVERSION_READ, L2CFNODE, &
      & SKIPRETRIEVAL, SPECIALDUMPFILE, MLSMESSAGE
    use MLSL2TIMINGS, only: SECTION_TIMES, TOTAL_TIMES, &
      & ADDPHASETOPHASENAMES, FILLTIMINGS, FINISHTIMINGS
    use MLSMESSAGEMODULE, only: MLSMSG_ERROR, MLSMSG_WARNING, &
      & MLSMSG_ALLOCATE, MLSMSG_DEALLOCATE, MLSMESSAGERESET
    use MLSPCF2, only: MLSPCF_L2APRIORI_START, MLSPCF_L2APRIORI_END
    use MLSRANDOMNUMBER, only: MLS_RANDOM_SEED, MATH77_RAN_PACK
    use MLSSTRINGLISTS, only: CATLISTS, GETHASHELEMENT, &
      & NUMSTRINGELEMENTS, PUTHASHELEMENT, &
      & STRINGELEMENT, STRINGELEMENTNUM
    use MLSSTRINGS, only: LOWERCASE
    use MOLECULES, only: L_H2O
    use MORETREE, only: GET_BOOLEAN, GET_FIELD_ID, GET_SPEC_ID
    use OUTPUT_M, only: BLANKS, OUTPUT, OUTPUTNAMEDVALUE, &
      & REVERTOUTPUT, SWITCHOUTPUT
    use PFADATA_M, only: FLUSH_PFADATA
    use QUANTITYTEMPLATES, only: QUANTITYTEMPLATE_T, &
      & MODIFYQUANTITYTEMPLATE
    use READAPRIORI, only: APRIORIFILES
    use SNOOPMLSL2, only: SNOOP
    use STRING_TABLE, only: GET_STRING
    use SUBSETMODULE, only: APPLYMASKTOQUANTITY, RESTRICTRANGE, &
      & SETUPFLAGCLOUD, SETUPSUBSET, UPDATEMASK
    use TIME_M, only: TIME_NOW
    use TOGGLES, only: GEN, LEVELS, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, &
      & SOURCE_REF, SUB_ROSA, SUBTREE
    use TREE_TYPES, only: N_NAMED
    use VECTORSMODULE, only: ADDVECTORTODATABASE, &
      & CLEARMASK, CLONEVECTORQUANTITY, CREATEVECTOR, &
      & DESTROYVECTORQUANTITYVALUE, DUMPQUANTITYMASK, &
      & GETVECTORQTYBYTEMPLATEINDEX, &
      & VALIDATEVECTORQUANTITY, VECTOR_T, &
      & VECTORTEMPLATE_T, VECTORVALUE_T, M_FILL
    use VGRIDSDATABASE, only: VGRIDS

    ! Dummy arguments
    integer, intent(in)                               :: ROOT ! Of the FILL section in the AST
    type (MLSFile_T), dimension(:), pointer           :: FILEDATABASE
    type (griddedData_T), dimension(:), pointer       :: GRIDDEDDATABASE
    type (vectorTemplate_T), dimension(:), pointer    :: VECTORTEMPLATES
    type (vector_T), dimension(:), pointer            :: VECTORS
    type (quantityTemplate_T), dimension(:), pointer  :: QTYTEMPLATES
    type (matrix_database_T), dimension(:), pointer   :: MATRICES
    type (Hessian_T), dimension(:), pointer           :: HESSIANS
    type (l2GPData_T), dimension(:), pointer          :: L2GPDATABASE
    type (l2AUXData_T), dimension(:), pointer         :: L2AUXDATABASE
    type(ForwardModelConfig_T), dimension(:), pointer :: FWMODELCONFIG
    type (mlSChunk_T), dimension(:), pointer          :: CHUNKS
    integer, intent(in)                               :: CHUNKNO

    ! -----     Declarations for Fill and internal subroutines     -------

    logical, parameter :: DEEBUG = .FALSE.                 ! Usually FALSE

    ! -999.99 ! Same as %template%badvalue
    real, parameter ::    UNDEFINED_VALUE = DEFAULTUNDEFINEDVALUE

    ! Local variables

    type (vectorValue_T), pointer :: APRIORIPRECISION
    type (vectorValue_T), pointer :: AQUANTITY
    type (vectorValue_T), pointer :: BASELINEQUANTITY
    type (vectorValue_T), pointer :: BNDPRESSQTY
    type (vector_T),      pointer :: DESTVECTOR
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
    type (vector_T),      pointer :: MEASVECTOR
    type (vectorValue_T), pointer :: MINNORMQTY
    type (vectorValue_T), pointer :: MODELQTY
    type (vector_T),      pointer :: MODELVECTOR
    type (vectorValue_T), pointer :: NBWQUANTITY
    type (vectorValue_T), pointer :: NOISEQTY
    type (vector_T),      pointer :: NOISEVECTOR
    type (vectorValue_T), pointer :: NORMQTY
    type (vectorValue_T), pointer :: NULLQUANTITY
    type (vectorValue_T), pointer :: ODQUANTITY              ! optical depth
    type (vectorValue_T), pointer :: ORBITINCLINATIONQUANTITY
    type (vectorValue_T), pointer :: PHITANQUANTITY
    type (vectorValue_T), pointer :: PRECISIONQUANTITY
    type (vectorValue_T), pointer :: PTANQUANTITY
    type (vectorValue_T), pointer :: RADQUANTITY
    type (vectorValue_T), pointer :: QUANTITY ! Quantity to be filled
    type (vectorValue_T), pointer :: RADIANCEQUANTITY
    type (vectorValue_T), pointer :: RATIOQUANTITY
    type (vectorValue_T), pointer :: REFGPHQUANTITY
    type (vectorValue_T), pointer :: REFGPHPRECISIONQUANTITY
    type (vectorValue_T), pointer :: RHIPRECISIONQUANTITY
    type (vectorValue_T), pointer :: SCECIQUANTITY
    type (vectorValue_T), pointer :: SCVELQUANTITY
    type (vectorValue_T), pointer :: SOURCEQUANTITY
    type (vector_T),      pointer :: SOURCEVECTOR
    type (vectorValue_T), pointer :: SYSTEMPQUANTITY
    type (vectorValue_T), pointer :: TEMPERATUREQUANTITY
    type (vectorValue_T), pointer :: TEMPPRECISIONQUANTITY
    type (vectorValue_T)          :: TEMPSWAPQUANTITY
    type (vectorValue_T), pointer :: TOTALPOWERQUANTITY
    type (vectorValue_T), pointer :: TNGTECIQUANTITY
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
    logical :: ASPERCENTAGE             ! Flag for noRadsPerMIF
    logical :: AUTOFILL                 ! Flag for automatically filling qtys
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
    character(len=80) :: BOOLEANNAME    ! E.g., 'BQTYS'
    integer :: BOXCARMETHOD             ! l_min, l_max, l_mean
    integer :: BQTYINDEX                ! Index of a quantity in vector
    integer :: BVECINDEX                ! Index of a vector
    real(rv) :: C                       ! constant "c" in manipulation
    logical :: CENTERVERTICALLY         ! For bin based fills
    integer :: CHANNEL                  ! For spreadChannels fill
    integer :: CHANNELSNODE
    integer :: COLMABUNITS              ! l_DOBSONUNITS, or l_MOLCM2
    integer :: COLVECTOR                ! Vector defining columns of Matrix
    type(matrix_SPD_T), pointer :: Covariance
    integer :: DESTINATIONVECTORINDEX   ! For transfer commands
    !                                     -- for Covariance
    integer :: EARTHRADIUSQTYINDEX
    integer :: EARTHRADIUSVECTORINDEX
    character(len=256) :: EXTRAOBJECTS  ! Which bright objects to avoid
    integer :: Diagonal                 ! Index of diagonal vector in database
    !                                     -- for Covariance
    character(len=16) :: DIMLIST        ! 's', 'c', or 'i' in manipulation's shift
    logical :: DONTMASK                 ! Use even masked values if TRUE
    integer :: ECRTOFOVQUANTITYINDEX    ! Rotation matrix
    integer :: ECRTOFOVVECTORINDEX      ! Rotation matirx
    integer :: ERRORCODE                ! 0 unless error; returned by called routines
    logical :: EXACT                    ! Set STATUS to value exactly, don't OR bits
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
    logical :: FORCE                    ! Fill in as many (instances) as will fit
    integer :: SHAPENODE              ! For the parser
    integer :: FRACTION                 ! Index of fraction vector in database
    logical :: FROMPRECISION            ! Fill from l2gpPrecision not l2gpValue
    integer :: GEOCALTITUDEQUANTITYINDEX    ! In the source vector
    integer :: GEOCALTITUDEVECTORINDEX      ! In the vector database
    integer :: GLOBALUNIT               ! To go into the vector
    character(len=16) :: GLSTR          ! geo. loc. in manipulation='..'
    integer :: GPHQUANTITYINDEX         ! In the source vector
    integer :: GPHVECTORINDEX           ! In the vector database
    logical, dimension(field_first:field_last) :: GOT
    integer :: GRIDINDEX                ! Index of requested grid
    integer :: GSON                     ! Descendant of Son
    integer :: HEIGHTNODE               ! Descendant of son
    character(len=8) :: HEIGHTRANGE     ! 'above', 'below', or ' '
    type(hessian_T), pointer :: hessian
    integer :: HESSIANINDEX             ! An index into the hessians
    integer :: HESSIANTOFILL             ! Index in database
    logical :: HIGHBOUND                ! Flag
    integer :: H2OQUANTITYINDEX         ! in the quantities database
    integer :: H2OVECTORINDEX           ! In the vector database
    integer :: H2OPRECISIONQUANTITYINDEX         ! in the quantities database
    integer :: H2OPRECISIONVECTORINDEX           ! In the vector database
    integer :: I, J                     ! Loop indices for section, spec, expr
    integer :: IBO
    logical :: IGNOREGEOLOCATION        ! Don't copy geolocation to vector qua
    logical :: IGNORENEGATIVE           ! Don't sum chi^2 at values of noise < 0
    logical :: IGNORETEMPLATE           ! Don't check compatible--just fill values
    logical :: IGNOREZERO               ! Don't sum chi^2 at values of noise = 0
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
    type(matrix_Cholesky_T) :: MATRIXCHOLESKY
    type(matrix_Kronecker_T) :: MATRIXKRONECKER
    type(matrix_SPD_T) :: MATRIXSPD
    type(matrix_T) ::     MATRIXPLAIN
    real(r8) :: MAXVALUE                ! Value of f_maxValue field
    integer :: MAXVALUEUNIT             ! Unit for f_maxValue field
    integer :: MANIPULATION             ! String index
    type(matrix_T), pointer :: MATRIX
    integer :: MATRIXTOFILL             ! Index in database
    integer :: MATRIXTYPE               ! Type of matrix, L_... from init_tables
    integer :: MAXITERATIONS            ! For hydrostatic fill
    integer :: Me = -1                  ! String index for trace
    integer :: Me_FillCommand = -1      ! String index for trace
    integer :: Me_FlagCloud = -1        ! String index for trace
    integer :: Me_RestrictRange = -1    ! String index for trace
    integer :: Me_Spec = -1             ! String index for trace
    integer :: Me_Subset = -1           ! String index for trace
    integer :: Me_UpdateMask = -1       ! String index for trace
    integer :: Me_Transfer = -1         ! String index for trace
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
    type (MLSFile_T), pointer   :: MLSFILE
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
    character(len=16) :: OPTIONS
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
    logical :: QUADRATURE               ! Apply baseline in quadrature
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
    logical :: REPEATLOOP               ! Do we repeat this section?
    logical :: RESETSEED                ! Let mls_random_seed choose new seed
    integer :: RHIPRECISIONQUANTITYINDEX         ! in the quantities database
    integer :: RHIPRECISIONVECTORINDEX           ! In the vector database
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
    integer, dimension(2) :: SHP
    logical :: SKIPMASK                 ! Flag for transfer
    integer :: SON                      ! Of root, an n_spec_args or a n_named
    integer :: SOURCE                   ! l_rows or l_colums for adoption
    integer :: SOURCEQUANTITYINDEX      ! in the quantities database
    integer :: SOURCEVECTORINDEX        ! In the vector database
    logical :: SKIPFILL                 ! Don't execute Fill command
    logical :: SPREADFLAG               ! Do we spread values accross instances in explict
    integer :: STATUS                   ! Flag from allocate etc.
    integer :: STATUSVALUE              ! Vaue of f_status
!   logical :: STRICT                   ! Maximize checking
    integer :: SUPERDIAGONAL            ! Index of superdiagonal matrix in database
    integer :: SURFNODE                 ! Descendant of son
    logical :: SWITCH2INTRINSIC         ! Have mls_random_seed call intrinsic
    !                                     -- for Covariance
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
    integer :: TOTALPOWERQUANTITYINDEX    ! In the quantities database
    integer :: TOTALPOWERVECTORINDEX      ! In the vector database
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
    integer, parameter :: whereRange  = 0
    integer :: WIDTH                    ! Width of boxcar

    ! Executable code
    timing = section_times
    if ( timing ) call time_now ( t1 )
    old_math77_ran_pack = math77_ran_pack

    call trace_begin ( me, "MLSL2Fill", root, cond=toggle(gen) )
    if ( specialDumpFile /= ' ' ) &
      & call switchOutput( specialDumpFile, keepOldUnitOpen=.true. )

    fillerror = 0
    templateIndex = -1
    vectorIndex = -1
    maxIterations = 4
    repeatLoop = .false. ! By default, we will not repeat
    call InitializeRepeat

    repeat_loop: do ! RepeatLoop
    ! Loop over the lines in the configuration file

    do i = 2, nsons(root)-1 ! Skip the section name at begin and end
      son = subtree(i,root)
      call trace_begin ( me_spec, "Fill.spec", son, &
        & cond=toggle(gen) .and. levels(gen) > 0 )
      if ( node_id(son) == n_named ) then ! Is spec labeled?
        key = subtree(2,son)
        vectorName = sub_rosa(subtree(1,son))
      else
        key = son
        vectorName = 0
      end if
      L2CFNODE = key
      if ( MLSSelecting .and. &
        & .not. any( get_spec_id(key) == (/ s_endselect, s_select, s_case /) ) ) cycle
      additional = .false.
      allowMissing = .false.
      asPercentage = .false.
      autoFill = .false.
      boxCarMethod = l_mean
      c = 0.
      centerVertically = .false.
      channel = 0
      channelsNode = 0
      colmabunits = l_molcm2 ! default units for column abundances
      dontMask = .false.
      exact = .false.
      excludeBelowBottom = .false.
      extinction = .false.
      fillMethod = l_none
      force = .false.
      fromPrecision = .false.
      GLStr = ' '
      got= .false.
      heightNode = 0
      heightRange = ' '
      ignoreZero = .false.
      ignoreNegative = .false.
      ignoreGeolocation = .false.
      ignoreTEMPLATE = .false.
      instancesNode = 0
      interpolate = .false.
      invert = .false.
      isPrecision = .false.
      maxValue = huge(0.0_r8)
      maxValueUnit = 0
      minValue = -huge(0.0_r8)
      minValueUnit = 0
      MissingGMAO = .false.
      logSpace = .false.
      options = ' '
      resetSeed = .false.
      refract = .false.
      scale = 0.0
      scaleInstances = -1.0
      scaleRatio = 1.0
      scaleSurfs = -1.0
      spreadFlag = .false.
!     strict = .false.
      surfNode = 0
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
      whereFill = .false.
      whereNotFill = .false.
      nullify ( nullQuantity, ptanQuantity, radQuantity, odQuantity )

      ! Node_id(key) is now n_spec_args.

      if ( get_spec_id(key) /= s_catchWarning ) &
        & call MLSMessageReset( clearLastWarning=.true. )
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
          case ( f_autoFill )
            autoFill = get_boolean ( fieldValue )
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
        ! Do we want to autofill any of its qtys?
        if ( autoFill .and. associated(vectors) ) &
          & call AutoFillVector ( vectors(size(vectors)) )

      case ( s_anygoodvalues )
        call decorate ( key, &
          & BooleanFromAnyGoodValues ( key, vectors ) )
      case ( s_anygoodradiances )
        call decorate ( key, &
          & BooleanFromAnyGoodRadiances ( key, chunks(chunkNo), filedatabase ) )
      case ( s_select ) ! ============ Start of select .. case ==========
        ! We'll start seeking a matching case
        call MLSSelect (key)
      case ( s_case ) ! ============ seeking matching case ==========
        ! We'll continue seeking a match unless the case is TRUE
        call MLSCase (key)
      case ( s_endSelect ) ! ============ End of select .. case ==========
        ! We'done with seeking a match
        call MLSEndSelect (key)
      case ( s_catchWarning )
        call decorate ( key,  BooleanFromCatchWarning ( key ) )
      case ( s_compare )
        call decorate ( key,  BooleanFromComparingQtys ( key, vectors ) )
      case ( s_computeTotalPower )
        call ComputeTotalPower ( key, vectors )
      case ( s_Reevaluate )
        call decorate ( key,  BooleanFromFormula ( 0, key, vectors ) )
      case ( s_diff, s_dump ) ! ======================== Diff, Dump ==========
        ! Handle disassociated pointers by allocating them with zero size
        status = 0
        if ( .not. associated(qtyTemplates) ) allocate ( qtyTemplates(0), stat=status )
        call test_allocate ( status, moduleName, 'QtyTemplates', (/0/), (/0/) )
        if ( .not. associated(vectorTemplates) ) allocate ( vectorTemplates(0), stat=status )
        call test_allocate ( status, moduleName, 'VectorTemplates', (/0/), (/0/) )
        if ( .not. associated(vectors) ) allocate ( vectors(0), stat=status )
        call test_allocate ( status, moduleName, 'Vectors', (/0/), (/0/) )
        call dumpCommand ( key, qtyTemplates, vectorTemplates, vectors, &
          & GriddedDataBase=GriddedDataBase, FileDataBase=FileDataBase, &
          & MatrixDatabase=Matrices, HessianDatabase=Hessians )
      case ( s_hessian ) ! ===============================  Hessian  =====
        got = .false.
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
          end select
        end do
        if ( got(f_columns) .and. got(f_rows) ) then
          ! We decorate with a negative number to differentiate hessians from
          ! Jacobians when it comes time to write them out
          ! *******************************************
          ! The above is an enormously bad idea
          ! The pronounced justification fails because
          ! Hessians and Jacobians are stored in separate databases
          ! databases.
          ! Long term: Don't decorate with a negative number
          ! which was an unacceptably crude hack even if it had been
          ! necessary, which it was not
          ! *******************************************
          call decorate ( key, -addHessianToDatabase(hessians, &
            & createEmptyHessian ( vectorName, &
            & vectors(rowVector), vectors(colVector), where=source_ref(key) )) )
        else
          call announce_error ( key, missingField, &
            & extraInfo = (/ f_columns, f_rows /) )
        end if

      case ( s_Repeat ) ! ============================== Repeat ==========
        ! We'll Repeat the section as long as the Boolean cond'n is TRUE
        RepeatLoop = Repeat(key)
        ! call outputNamedValue ( 'repeat loop?', repeatLoop )
        if ( .not. RepeatLoop ) exit repeat_loop
        call nextRepeat
      case ( s_skip ) ! ============================== Skip ==========
        ! We'll skip the rest of the section if the Boolean cond'n is TRUE
        if ( Skip(key) ) exit repeat_loop
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
        call trace_begin ( me_subset, "Fill.subset", root, &
          & cond=toggle(gen) .and. levels(gen) > 1 )
        call SetupSubset ( key, vectors )
        call trace_end ( cond=toggle(gen) .and. levels(gen) > 1 )

      case ( s_restrictRange )
        call trace_begin ( me_restrictRange, "Fill.RestrictRange", root, &
          & cond=toggle(gen) .and. levels(gen) > 1 )
        call RestrictRange ( key, vectors )
        call trace_end ( cond=toggle(gen) .and. levels(gen) > 1 )

      case ( s_flushL2PCBins )
        call FlushLockedBins

      case ( s_flushPFA )
        call flush_PFAData ( key, status )
        fillerror = max(fillerror,status)

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
          case ( f_hessian )
            hessianToFill = decoration(decoration( subtree(2, gson) ))
          case ( f_vector )
            vectorIndex = decoration(decoration( subtree(2, gson) ))
          case ( f_bin )
            binName = sub_rosa ( subtree(2,gson) )
          case ( f_source )
            source = decoration ( subtree(2,gson) )
          end select
        end do
        if ( got ( f_hessian ) ) then
          if ( got ( f_vector ) ) call Announce_Error ( key, no_Error_Code, &
            & 'Cannot load both vector and hessian (Why not?)' )
          Hessian => hessians(hessianToFill)
          call LoadHessian ( Hessian, binName, message )
          if ( len_trim(message) > 0 ) call output( '* * *' // message, advance='yes' )
          if ( got ( f_matrix ) ) then
            call GetFromMatrixDatabase ( matrices(matrixToFill), matrix )
            call LoadMatrix ( matrix, binName, message )
          endif
        elseif ( got ( f_matrix ) ) then
          if ( got ( f_vector ) ) call Announce_Error ( key, no_Error_Code, &
            & 'Cannot load both vector and matrix (Why not?)' )
          call GetFromMatrixDatabase ( matrices(matrixToFill), matrix )
          call LoadMatrix ( matrix, binName, message )
          if ( len_trim(message) > 0 ) call output( '* * *' // message, advance='yes' )
        else if ( got ( f_vector ) ) then
          if ( .not. got ( f_source ) ) call Announce_Error ( key, no_Error_Code, &
            & 'Must supply source=rows/columns for vector adoption' )
          call LoadVector ( vectors(vectorIndex), binName, source, message )
        else
          call Announce_Error ( key, no_Error_Code, &
            & 'Must supply hessian, matrix or vector to adopt' )
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
        call trace_begin ( me_updateMask, "Fill.UpdateMask", root, &
          & cond=toggle(gen) .and. levels(gen) > 1 )
        call UpdateMask ( key, vectors )
        call trace_end ( cond=toggle(gen) .and. levels(gen) > 1 )

      case ( s_flagCloud )
        call trace_begin ( me_flagCloud, "Fill.flagCloud", root, &
          & cond=toggle(gen) .and. levels(gen) > 1 )
        call SetupflagCloud ( key, vectors )
        call trace_end ( cond=toggle(gen) .and. levels(gen) > 1 )

      case ( s_directRead ) ! =======================  directRead  =====
        call directReadCommand
      case ( s_fill ) ! ===================================  Fill  =====
        call trace_begin ( me_fillCommand, "Fill.fillCommand", key, &
          & cond=toggle(gen) .and. levels(gen) > 1 )
        call fillCommand
        call trace_end ( cond=toggle(gen) .and. levels(gen) > 1 )
      case ( s_fillcovariance ) ! ===================  Covariance  =====
        invert = .false. ! Default if the field isn't present
        lengthScale = 0
        fraction = 0
        do j = 2, nsons(key)
          gson = subtree(j,key) ! The argument
          fieldIndex = get_field_id(gson)
          if ( nsons(gson) > 1) gson =  &
            & decoration(decoration(subtree(2,gson))) ! Now value of said argument
          select case ( fieldIndex )
          case ( f_matrix )
            matrixToFill = gson
            if ( getKindFromMatrixDatabase(matrices(matrixToFill)) /= k_spd ) &
              call announce_error ( key, notSPD )
          case ( f_diagonal )
            diagonal = gson
          case ( f_ignoreTemplate )
            ignoreTemplate = get_boolean ( gson )
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
        if ( SKIPRETRIEVAL ) then
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Unable to fill covariance when skipping retrievals' )
        else
          call FillCovariance ( covariance, &
            & vectors, diagonal, lengthScale, fraction, invert, ignoreTemplate )
        endif

      case ( S_FILLDIAGONAL ) ! ===============  Diagonal  =====
        do j = 2, nsons(key)
          gson = subtree(j,key) ! The argument
          fieldIndex = get_field_id(gson)
          if ( nsons(gson) > 1) gson = &
            & decoration(decoration(subtree(2,gson))) ! Now value of said argument
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
        call destroyCommand ( key, matrices, vectors, griddedDataBase )

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
        call trace_begin ( me_transfer, "Fill.Transfer", key, &
          & cond=toggle(gen) .and. levels(gen) > 1 )
        nullify( destVector, measvector, modelvector, noisevector, sourceVector )
        interpolate = .false.
        skipMask = .false.
        booleanName = ' '
        do j = 2, nsons(key)
          gson = subtree(j,key)  ! The argument
          fieldIndex = get_field_id(gson)
          got(fieldIndex) = .true.
          if ( nsons(gson) > 1 ) then
            gson = subtree(2,gson) ! Now the value of said argument
            fieldValue = decoration(gson) ! The field's value
          else
            fieldValue = gson
          end if
          select case ( fieldIndex )
          case ( f_a )
            aVecIndex = decoration(fieldValue)
          case ( f_b )
            bVecIndex = decoration(fieldValue)
          case(f_c)
            call expr ( gson , unitAsArray, valueAsArray )
            c = valueAsArray(1)
          case ( f_source )
            sourceVectorIndex = decoration(fieldValue)
            sourceVector => vectors( sourceVectorIndex )
          case ( f_destination )
            destinationVectorIndex = decoration(fieldValue)
            destVector => vectors( destinationVectorIndex )
          case ( f_dontMask )
            dontMask = get_boolean ( gson )
          case ( f_ignoreNegative )
            ignoreNegative = get_boolean ( gson )
          case ( f_ignoreZero )
            ignoreZero = get_boolean ( gson )
          case ( f_interpolate )
            interpolate = get_boolean ( fieldValue )
          case ( f_manipulation )
            manipulation = sub_rosa ( gson )
          case ( f_method )   ! How ? (if not default copy)
            fillMethod = decoration(gson)
          case ( f_measurements )   ! Only used for diagnostic special fills
            measVectorIndex = decoration(fieldValue)
            measVector => vectors( measVectorIndex )
          case ( f_model )   ! Only used for diagnostic special fills
            modelVectorIndex = decoration(fieldValue)
            modelVector => vectors( modelVectorIndex )
          case ( f_noise )   ! Only used for chi^2 special fills or addnoise
            noiseVectorIndex = decoration(fieldValue)
            noiseVector => vectors( noiseVectorIndex )
          case ( f_PtanQuantity ) ! For minorframe qty
            PtanVectorIndex = decoration(decoration(subtree(1,gson)))
            PtanQuantityIndex = &
              & decoration(decoration(decoration(subtree(2,gson))))
          case (f_quantityNames)
            call get_string ( sub_rosa(gson), booleanName, strip=.true. )
            booleanName = lowerCase(booleanName)
          case ( f_skipMask )
            skipMask = get_boolean ( fieldValue )
          case default ! Can't get here if type checker worked
          end select
        end do
        if ( got(f_a) ) then
          if ( .not. got(f_b) ) bVecIndex = aVecIndex
          if ( got(f_c) ) then
            call ManipulateVectors ( manipulation, &
              & vectors(destinationVectorIndex), &
              & vectors(aVecIndex), vectors(bVecIndex), c, &
              & booleanName=booleanName )
          else
            call ManipulateVectors ( manipulation, &
              & vectors(destinationVectorIndex), &
              & vectors(aVecIndex), vectors(bVecIndex), &
              & booleanName=booleanName )
          endif
        else if ( fillMethod /= l_none ) then
          if ( got ( f_ptanQuantity )  ) then
            ptanQuantity => GetVectorQtyByTemplateIndex( &
              & vectors(ptanVectorIndex), ptanQuantityIndex)
          endif
          call TransferVectorsByMethod ( key, destvector, &
            & sourceVector, fillMethod, dontMask, interpolate, &
            & ignorenegative, ignoreZero, measVector, modelVector, &
            & noiseVector, ptanQuantity, booleanName )
        else if ( got(f_source) ) then
          call TransferVectors ( vectors(sourceVectorIndex), &
            & vectors(destinationVectorIndex), skipMask, interpolate, &
            & booleanName )
        else
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Transfer command requires either source or a to be present' )
        end if
        call trace_end ( "Fill.Transfer", cond=toggle(gen) .and. levels(gen) > 1 )

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
        call Snoop ( key=key, vectorDatabase=vectors, &
          & matrixDatabase=snoopMatrices )
        deallocate ( snoopMatrices, STAT=status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_Deallocate//'snoopMatrices' )

      case ( s_StreamlineHessian ) ! =============== StreamlineHessian =====
        call doStreamline

      ! End of fill operations

      case default ! Can't get here if tree_checker worked correctly
      end select
      call trace_end ( "Fill.spec", cond=toggle(gen) .and. levels(gen) > 0 )
    end do ! End of loop of specs
    if ( .not. repeatLoop ) exit ! Otherwise, we will repeat the section
    end do  repeat_loop!  RepeatLoop

    if ( fillError /= 0 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, 'Problem with Fill section' )
    end if

    if ( specialDumpFile /= ' ' ) call revertOutput
    call trace_end ( "MLSL2Fill", cond=toggle(gen) )
    math77_ran_pack = old_math77_ran_pack
    if ( timing ) call sayTime

  ! =====     Internal Procedures     ==================================

  contains

    ! ---------------------------------------------- DoStreamline  -----
    subroutine DoStreamline
      use INIT_TABLES_MODULE, only: F_GEODANGLE, F_HESSIAN, &
        & F_SCALEHEIGHT, F_SURFACE, F_THRESHOLD
      real(r8) :: GEODANGLE        ! For StreamlineHessian
      real(r8) :: SCALEHEIGHT      ! Scale height for StreamlineHessian
      integer :: SURFACE           ! Number of surfaces for StreamlineHessian
      real(r8) :: THRESHOLD
      geodAngle = -1.0   ! means "not specified"
      scaleHeight = -1.0 ! means "not specified"
      surface = -1       ! means "not specified"
      threshold = -1.0   ! means "not specified"
      do j = 2, nsons(key)
        gson = subtree(j,key) ! The argument
        fieldIndex = get_field_id(gson)
        select case ( fieldIndex )
        case ( f_geodAngle )
          call expr_check ( subtree(2,gson), unitAsArray, valueAsArray, &
            & (/PHYQ_Angle/), unitsError )
          if ( unitsError ) call Announce_error ( subtree(j,key), wrongUnits, &
            & extraInfo=(/unitAsArray(1), PHYQ_Angle/) )
          geodAngle = valueAsArray(1)
        case ( f_hessian )
          hessianIndex = -decoration(decoration(subtree(2,gson)))
        case ( f_scaleHeight )
          call expr_check ( subtree(2,gson), unitAsArray, valueAsArray, &
            & (/PHYQ_Dimensionless/), unitsError )
          if ( unitsError ) call Announce_error ( subtree(j,key), wrongUnits, &
            & extraInfo=(/unitAsArray(1), PHYQ_Dimensionless/) )
          scaleHeight = valueAsArray(1)
        case ( f_surface )
          call expr_check ( subtree(2,gson), unitAsArray, valueAsArray, &
            & (/PHYQ_Dimensionless/), unitsError )
          if ( unitsError ) call Announce_error ( subtree(j,key), wrongUnits, &
            & extraInfo=(/unitAsArray(1), PHYQ_Dimensionless/) )
          surface = valueAsArray(1)
        case ( f_threshold )
          call expr_check ( subtree(2,gson), unitAsArray, valueAsArray, &
            & (/PHYQ_Dimensionless/), unitsError )
          if ( unitsError ) call Announce_error ( subtree(j,key), wrongUnits, &
            & extraInfo=(/unitAsArray(1), PHYQ_Dimensionless/) )
          threshold = valueAsArray(1)
        end select
      end do
      call StreamlineHessian ( hessians ( hessianIndex ), &
        & surface, scaleHeight, geodAngle, threshold )
    end subroutine DoStreamline

    ! ------------------------------------------------ directReadCommand -----
    subroutine directReadCommand
    use L2PC_M, only: READCOMPLETEHDF5L2PCFILE
    ! Now we're on actual directRead instructions.
    logical, parameter :: DEEBUG = .true.
      integer :: EXPRUNITS(2)             ! From expr
      real (r8) :: EXPRVALUE(2)           ! From expr
      integer :: fieldIndex
      integer :: file                     ! Index into string table
      integer :: fileType
      character(len=8) :: fileTypeStr
      integer :: gson
      integer :: hdfversion
      logical :: interpolate
      integer :: j
      character(len=32) :: sdname
      logical :: spread
      type (vector_T), pointer :: Vector
      ! Loop over the instructions to the directRead command
      got = .false.
      hdfVersion = DEFAULT_HDFVERSION_READ
      file = 0
      options = ' '
      binname = 0
      sdname = ' '
      spread = .false.
      interpolate = .false.
      if ( DEEBUG ) call output ( 'In DirectReadCommand', advance='yes' )
      do j = 2, nsons(key)
        gson = subtree(j,key) ! The argument
        fieldIndex = get_field_id(gson)
        if ( nsons(gson) > 1) gson = subtree(2,gson) ! Now value of said argument
        got(fieldIndex)=.TRUE.
        select case ( fieldIndex )
        case ( f_quantity )
          vectorIndex = decoration(decoration(subtree(1,gson)))
          quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        case ( f_vector )
          vectorIndex = decoration(decoration(gson))
          if ( DEEBUG ) call output ( 'Processing vector field', advance='yes' )
        case ( f_matrix )
          matrixToFill = decoration(decoration(gson))
          if ( DEEBUG ) call output ( 'Processing matrix field', advance='yes' )
        case ( f_hessian )
          hessianToFill = decoration(decoration(gson))
          if ( DEEBUG ) call output ( 'Processing hessian field', advance='yes' )
        case ( f_hdfVersion )
          if ( DEEBUG ) call output ( 'Begin Processing hdfversion field', advance='yes' )
          call expr ( gson, exprUnits, exprValue )
          if ( DEEBUG ) call output ( 'Processing hdfversion field', advance='yes' )
          if ( exprUnits(1) /= phyq_dimensionless ) &
            & call Announce_error ( son, NO_ERROR_CODE, &
            & 'No units allowed for hdfVersion: just integer 4 or 5')
          hdfVersion = exprValue(1)
        case ( f_bin )
          binName = sub_rosa ( gson )
        case ( f_file )
          file = sub_rosa(gson)
          if ( DEEBUG ) call outputNamedValue ( 'Processing file field', file )
        case ( f_options )
          if ( DEEBUG ) call output ( 'Begin Processing options field', advance='yes' )
          call get_string ( sub_rosa(gson), options, strip=.true. )
          if ( DEEBUG ) call output ( 'Processing options field', advance='yes' )
        case ( f_sdName )
          if ( DEEBUG ) call output ( 'Begin Processing sdName field', advance='yes' )
          call get_string ( sub_rosa(gson), sdName, strip=.true. )
          if ( DEEBUG ) call output ( 'Processing sdName field', advance='yes' )
        case ( f_interpolate )
            interpolate = get_boolean ( gson )
        case ( f_spread )
            spread = get_boolean ( gson )
        case ( f_type )
          if ( DEEBUG ) call output ( 'Begin Processing type field', advance='yes' )
          fileType = decoration(gson)
          call get_string ( lit_indices(fileType), fileTypeStr, strip=.true. )
          if ( DEEBUG ) call output ( 'Processing type field', advance='yes' )
        end select
      end do
      if ( DEEBUG ) call output ( 'Done processing fields', advance='yes' )
      if ( .not. got(f_file) ) call Announce_error ( key, 0, 'No file supplied' )
      if ( .not. got(f_type) ) call Announce_error ( key, 0, 'No type supplied' )
      if ( .not. any( (/ &
        & got(f_quantity), got(f_vector), got(f_matrix), got(f_hessian) /) ) )&
        & call Announce_error ( key, 0, &
        & 'Must supply one of matrix, quantity or vector' )
      if ( DEEBUG ) call output ( 'About to  get_file_name', advance='yes' )
      call get_file_name ( file, mlspcf_l2apriori_start, &
            & fileType, filedatabase, MLSFile, &
            & 'DirectRead File not found in PCF', mlspcf_l2apriori_end )
      if ( DEEBUG ) call output ( 'Done with get_file_name', advance='yes' )
      if ( got(f_quantity) ) then
        quantity => GetVectorQtyByTemplateIndex( &
          & vectors(vectorIndex), quantityIndex )
        call QtyFromFile ( key, quantity, MLSFile, &
          & filetypestr, options, sdName, spread, interpolate )
      elseif ( got ( f_hessian ) ) then
        Hessian => hessians(hessianToFill)
        call ReadCompleteHDF5L2PCFile ( MLSFile, key, shallow=.false. )
        call LoadHessian ( Hessian, binName, message )
        if ( len_trim(message) > 0 ) call output( '* * *' // message, advance='yes' )
        if ( got ( f_matrix ) ) then
          call GetFromMatrixDatabase ( matrices(matrixToFill), matrix )
          call LoadMatrix ( matrix, binName, message )
        endif
      elseif ( got(f_matrix) ) then
        if ( .not. got(f_bin) ) call Announce_error ( key, 0, 'No bin supplied' )
        call GetFromMatrixDatabase ( matrices(matrixToFill), matrix )
        call ReadCompleteHDF5L2PCFile ( MLSFile, key, shallow=.false. )
        call LoadMatrix ( matrix, binName, message )
        if ( len_trim(message) > 0 ) call output( '* * *' // message, advance='yes' )
      else
        vector => vectors(vectorIndex)
        call VectorFromFile ( key, vector, MLSFile, &
          & filetypestr, options, spread, interpolate )
      endif
    end subroutine directReadCommand

    ! ------------------------------------------------ fillCommand -----
    subroutine fillCommand
    ! Now we're on actual Fill instructions.
      use Init_Tables_Module, only: L_None
      integer :: GEOLOCATION
      integer :: JJ
      integer :: MUL
      integer, dimension(3) :: START, COUNT, STRIDE, BLOCK
      logical :: QTYWASMASKED
      type(vectorValue_T) :: TEMPQUANTITY  ! For storing original qty's mask
      ! Executable
      ! Loop over the instructions to the Fill command
      BOMask = 0
      AvoidObjects = ' '
      dimList = 'c' ! defaults to shift or slip by channel, or surface if noFreqs < 2
      geolocation = l_none
      got = .false.
      multiplier = 1.
      start = 0
      count = 0
      stride = 0
      block = 0
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
        case ( f_asPercentage )
          asPercentage = get_boolean ( gson )
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
        case ( f_channels )
          channelsNode = son
        case ( f_centerVertically )
          centerVertically = get_boolean ( gson )
        case ( f_dimList )
          ! dimList = sub_rosa ( gson )
          call get_string ( sub_rosa ( gson ), dimList, strip=.true. )
          dimList = lowercase( dimList )
        case ( f_earthRadius ) ! For losGrid fill
          earthRadiusVectorIndex = decoration(decoration(subtree(1,gson)))
          earthRadiusQtyIndex = decoration(decoration(decoration(subtree(2,gson))))
        case ( f_excludeBelowBottom )
          excludeBelowBottom = get_boolean ( gson )
        case ( f_ECRToFOV ) ! For hydrostatic
          ecrToFOVVectorIndex = decoration(decoration(subtree(1,gson)))
          ecrToFOVQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        case ( f_exact )
          exact = get_boolean ( gson )
        case ( f_explicitValues ) ! For explicit fill
          valuesNode = subtree(j,key)
        case ( f_expr )
          call Announce_Error ( gson, noCodeFor, extraInfo=(/ f_expr /) )
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
        case ( f_shape )
          shapeNode = subtree(j,key)
        case ( f_fromPrecision )
          fromPrecision = get_boolean ( gson )
        case ( f_geocAltitudeQuantity ) ! For hydrostatic
          geocAltitudeVectorIndex = decoration(decoration(subtree(1,gson)))
          geocAltitudeQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        case ( f_geolocation )
          geolocation = decoration(gson)
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
        case ( f_heightRange )
          manipulation = sub_rosa ( gson )
          heightRange = ' '
          ! If heightRange field was present, it should have been one of
          ! 'a[bove]' meaning fill heights above supplied value
          ! 'b[elow]' meaning fill heights below supplied value
          call get_string ( manipulation, heightRange, strip=.true. )
          select case ( heightRange(1:1) )
          case ( 'a' )
            heightRange = 'above'
          case ( 'b' )
            heightRange = 'below'
          case ( ' ' )
            heightRange = ' '
          case default
            call Announce_Error ( key, no_Error_Code, &
            & 'invalid heightRange: ' // trim(heightRange) )
          end select
        case ( f_ignoreZero )
          ignoreZero = get_boolean ( gson )
        case ( f_ignoreGeolocation ) ! For l2gp etc. fill
          ignoreGeolocation =get_boolean ( gson )
        case ( f_ignoreNegative )
          ignoreNegative = get_boolean ( gson )
        case ( f_ignoreTemplate )
          ignoreTemplate = get_boolean ( gson )
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
          ! Either multiplier = [a, b] or multiplier = b are possible
          multiplier = UNDEFINED_VALUE
          do jj=1, min(nsons(multiplierNode)-1, 2)
            call expr(subtree(jj+1,multiplierNode),unitAsArray,valueAsArray)
            multiplier(jj) = valueAsArray(1)
          end do
          if ( DEEBUG ) then
            call output('Using multipliers: ', advance='no')
            call output(multiplier, advance='yes')
          end if
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
        case ( f_PtanQuantity ) ! For minorframe qty
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
        case ( f_rhiPrecisionQuantity ) ! For converting to h2o precision
          rhiPrecisionVectorIndex = decoration(decoration(subtree(1,gson)))
          rhiPrecisionQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
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
        case ( f_start, f_count, f_stride, f_block ) ! For selecting hyperslab
          multiplierNode = subtree(j,key)
          ! Either start = [a, b] or start = b are possible
          do jj=1, min(nsons(multiplierNode)-1, 2)
            call expr(subtree(jj+1,multiplierNode),unitAsArray,valueAsArray)
            mul = valueAsArray(1)
            select case ( fieldIndex )
            case ( f_start )
              start(jj) = mul
            case ( f_count )
              count(jj) = mul
            case ( f_stride )
              stride(jj) = mul
            case ( f_block )
              block(jj) = mul
            end select
          end do
          if ( DEEBUG ) then
            call output('index start: ', advance='no')
            call output(valueAsArray, advance='yes')
          end if
        case ( f_status )
          valuesNode = subtree(j,key)
          call expr_check ( gson , unitAsArray, valueAsArray, &
            & (/PHYQ_Dimensionless/), unitsError )
          if ( unitsError ) call Announce_error ( valuesNode, wrongUnits, &
            & extraInfo=(/unitAsArray(1), PHYQ_Dimensionless/) )
          statusValue = nint ( valueAsArray(1) )
        ! case ( f_strict )
        !   strict = get_boolean ( gson )
        case ( f_suffix )
          suffix = sub_rosa ( gson )
        case ( f_surface )
          surfNode = subtree(j,key)
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
        case ( f_totalPowerQuantity )
          totalpowerVectorIndex = decoration(decoration(subtree(1,gson)))
          totalpowerQuantityIndex = decoration(decoration(decoration(subtree(2,gson))))
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
      ! This is the actual quantity we shall Fill
      quantity => GetVectorQtyByTemplateIndex( &
        & vectors(vectorIndex), quantityIndex )
      if ( got ( f_ptanQuantity )  ) then
        ptanQuantity => GetVectorQtyByTemplateIndex( &
          & vectors(ptanVectorIndex), ptanQuantityIndex)
      endif

      ! However, we will first mask the quantity to account for
      ! any height, instances, etc. constraints you may have set
      qtyWasMasked = associated( quantity%mask )
      call CloneVectorQuantity ( tempQuantity, quantity )
      ! Ignore original mask if /dontMask 
      if ( qtyWasMasked .and. dontMask .and. fillMethod /= l_l1b ) &
        & quantity%mask = char(0)
      if ( any( fillMethod == (/ &
        & l_status, l_quality /) ) ) then
        call ApplyMaskToQuantity( quantity, &
          & radQuantity, ptanQuantity, odQuantity, nullQuantity, 0._r8, &
          & maxValue, minValue, heightRange, whereRange, &
          & .false., .false., additional=.true., reset=.false., &
          & maskBit=M_Fill, heightNode=0, surfNode=0, &
          & instancesNode=instancesNode, channelsNode=0 )
      else
        call ApplyMaskToQuantity( quantity, &
          & radQuantity, ptanQuantity, odQuantity, nullQuantity, 0._r8, &
          & maxValue, minValue, heightRange, whereRange, &
          & .false., .false., additional=.true., reset=.false., &
          & maskBit=M_Fill, heightNode=heightNode, surfNode=surfNode, &
          & instancesNode=instancesNode, channelsNode=channelsNode )
      endif
      if ( heightNode /= 0 .and. DEEBUG ) call dumpQuantityMask ( quantity )
      ! Then we Fill the newly masked quantity
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
          & noiseQty, multiplier, spreadFlag, ignoreTemplate )

      case ( l_applyBaseline )
        if ( .not. got ( f_baselineQuantity ) ) &
          & call Announce_Error ( key, no_Error_Code, &
          & 'Need baselineQuantity for applyBaseline fill' )
        baselineQuantity => GetVectorQtyByTemplateIndex( &
          & vectors(baselineVctrIndex), baselineQtyIndex )
        call ApplyBaseline ( key, quantity, baselineQuantity, &
          & quadrature, dontMask, ignoreTemplate )

      case ( l_asciiFile )
        if ( .not. got ( f_file ) ) &
          & call Announce_Error ( key, no_Error_Code, &
          & 'Need filename for asciiFile fill' )
        if ( got ( f_badRange ) ) then
          call FromASCIIFile ( key, quantity, filename, badRange )
        else
          call FromASCIIFile ( key, quantity, filename )
        end if

      case ( l_binMax, l_binMean, l_binMin, l_binTotal, &
           & l_lsLocal, l_lsGlobal, l_lsWeighted )
        if ( .not. got ( f_sourceQuantity ) ) &
          & call Announce_Error ( key, no_Error_Code, &
          & 'Need source quantity for bin fill or least-squares fill' )
        sourceQuantity => GetVectorQtyByTemplateIndex( &
          & vectors(sourceVectorIndex), sourceQuantityIndex )

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

        if ( fillerror == 0 ) then
          select case ( fillMethod )
          case ( l_binMax, l_binMean, l_binMin, l_binTotal )
            call WithBinResults ( key, quantity, sourceQuantity, ptanQuantity, &
              & channel, fillMethod, additional, excludeBelowBottom, centerVertically )
          case default
            call UsingLeastSquares ( key, quantity, sourceQuantity, ptanQuantity, &
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
        if ( width == 1 ) then
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Boxcar Fill with width=1 results in straightforward copy' )
          call FromAnother ( quantity, sourceQuantity, ptanQuantity, &
            & key, ignoreTemplate=.true., spreadflag=.false., &
            & interpolate=.false., force=.false. )
        elseif (  mod ( width, 2 ) == 0 ) then
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Boxcar Fill called with even width: adding one to make it odd' )
          call WithBoxcarFunction ( key, quantity, sourceQuantity, width+1, &
            & boxCarMethod, ignoreTemplate )
        else
          call WithBoxcarFunction ( key, quantity, sourceQuantity, width, &
            & boxCarMethod, ignoreTemplate )
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
          call ChiSqChan ( key, quantity, &
            & measQty, modelQty, noiseQty, &
            & dontMask, ignoreZero, ignoreNegative, ignoreTemplate, multiplier )
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
          call ChiSqMMaf ( key, quantity, &
            & measQty, modelQty, noiseQty, &
            & dontMask, ignoreZero, ignoreNegative, ignoreTemplate, multiplier )
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
          call ChiSqMMif ( key, quantity, &
            & measQty, modelQty, noiseQty, &
            & dontMask, ignoreZero, ignoreNegative, ignoreTemplate, multiplier )
        end if
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
        call ChiSqRatio ( key, &
          & quantity, normQty, minNormQty, flagQty, dontMask, ignoreTemplate )

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
          else
            call get_string( lit_indices(quantity%template%molecule), mol, &
              & strip=.true. )
            call GetHashElement (col_species_keys, &
              & col_species_hash, trim(lowercase(mol)), &
              & ExplicitUnit, .true.)
            if ( index(lowerCase(ExplicitUnit), 'd') > 0 ) colmabunits = l_DU
          end if
          call ColAbundance ( key, quantity, &
            & bndPressQty, vmrQty, colmAbUnits, ignoreTemplate )
        end if

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
        call ConvergenceFromChisq ( key, quantity, sourceQuantity, scale, &
          & ignoreTemplate )

      case ( l_derivative )
        if ( .not. any ( got( &
          & (/ f_sourceQuantity, f_geocAltitudeQuantity, f_dimList /) &
          & ) ) ) &
          & call Announce_Error ( key, no_Error_Code, &
          & 'Need source quantity, geocAltitudeQuantity, dimList for ' // &
          & 'derivative fill' )
        sourceQuantity => GetVectorQtyByTemplateIndex( &
          & vectors(sourceVectorIndex), sourceQuantityIndex )
        geocAltitudeQuantity => GetVectorQtyByTemplateIndex( &
          & vectors(geocAltitudeVectorIndex), geocAltitudeQuantityIndex)
        call DerivativeOfSource ( quantity, sourceQuantity, geocAltitudeQuantity, &
          & dimList, ignoreTemplate )

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
        call WithEstNoise ( &
          & quantity, radianceQuantity, sysTempQuantity, nbwQuantity, &
          & integrationTime )

      case ( l_explicit ) ! ---------  Explicitly fill from l2cf  -----
        if ( .not. got(f_explicitValues) ) &
          & call Announce_Error ( key, noExplicitValuesGiven )
        call Explicit ( quantity, valuesNode, spreadFlag, force, &
          & vectors(vectorIndex)%globalUnit, channel, &
          & .false., options=heightRange(1:1) )
      case ( l_extractChannel )
        if ( .not. all(got ( (/f_sourceQuantity,f_channel/)))) &
          & call Announce_Error ( key, no_Error_Code, &
            & 'Need sourceQuantity and channel for this fill' )
        sourceQuantity => GetVectorQtyByTemplateIndex( &
          & vectors(sourceVectorIndex), sourceQuantityIndex )
        call ExtractSingleChannel ( key, quantity, sourceQuantity, channel, &
          & ignoreTemplate )

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
        call FoldedRadiance ( quantity, lsb, usb, lsbFraction, usbFraction, key )

      case ( l_fwdModelTiming ) ! --- Fill timings for forward model  -----
        call FillFwdModelTimings (quantity%values(:,1), FWModelConfig, 'fwdTiming')
      case ( l_fwdModelMean ) ! --- Fill mean for forward model  -----
        call FillFwdModelTimings (quantity%values(:,1), FWModelConfig, 'mean')
      case ( l_fwdModelStdDev ) ! --- Fill std_dev for forward model  -----
        call FillFwdModelTimings (quantity%values(:,1), FWModelConfig, 'stdDev')
      case ( l_geoLocation )
        ! otherwise hard-to-obtain geo location data, settings, switches, etc.
        if ( .not. got ( f_manipulation ) ) &
          & call Announce_error ( key, no_Error_Code,'manipulation not supplied' )
        call get_string ( manipulation, GLStr, strip=.true. )
        select case (lowercase( trim(GLStr) ))
        ! Things about our chunk
        case ( 'abandoned' )
          quantity%values = 0
          if ( Chunks(ChunkNo)%abandoned ) quantity%values = 1
        case ( 'chunk' )
          quantity%values = ChunkNo
        case ( '1stmaf' )
          quantity%values = Chunks(ChunkNo)%firstMAFIndex
        ! Things about our Quantity template
        case ( 'frequencies' )
          quantity%values(:,1) = quantity%template%frequencies
        case ( 'geodlat' )
          quantity%values = quantity%template%geodLat
        case ( 'lon' )
          quantity%values = quantity%template%lon
        case ( 'losangle' )
          quantity%values = quantity%template%losangle
        case ( 'phi' )
          quantity%values = quantity%template%phi
        case ( 'solartime' )
          quantity%values = quantity%template%solartime
        case ( 'solarzenith' )
          quantity%values = quantity%template%solarzenith
        case ( 'surfs' )
          quantity%values = quantity%template%surfs
        case ( 'time' )
          quantity%values = quantity%template%time
        case default
          call Announce_error ( key, no_Error_Code, trim(GLStr) // &
            & 'geolocation not recognized' )
        end select

      case ( l_gather, l_scatter ) ! ------------  Gather  -----
        if ( .not. all(got( &
          &(/ f_sourceQuantity, f_start, f_count, f_stride, f_block /) &
          & )) ) &
          & call Announce_Error ( key, no_Error_Code, &
            & 'Need sourceQuantity, start, count, stride, and block for this fill' )
        sourceQuantity => GetVectorQtyByTemplateIndex( &
          & vectors(sourceVectorIndex), sourceQuantityIndex )
        if ( fillMethod == l_gather ) then
          call Gather( quantity, sourceQuantity, start, count, stride, block )
        else
          call Scatter( quantity, sourceQuantity, start, count, stride, block )
        endif

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

        call GPHPrecision ( key, quantity, tempPrecisionQuantity, &
          & refGPHPrecisionQuantity )

      case ( l_gridded ) ! ------------  Fill from gridded data  -----
        if ( .not. got(f_sourceGrid) ) &
          & call Announce_Error ( key, noSourceGridGiven )
        ! call output( 'Filling quantity from grid', advance='yes' )
        call FromGrid &
          & ( quantity, griddedDataBase(gridIndex), allowMissing, errorCode )
        ! call outputNamedValue( 'error code', errorCode )
        if ( errorCode /= 0 ) call Announce_error ( key, errorCode )

      case ( l_l1b ) ! --------------------  Fill from L1B data  -----
        if ( any( got( &
          & (/ f_height, f_heightRange, f_surface, f_instances, &
          & f_maxValue, f_minValue /) &
          & ) ) ) then
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Mask bits set during l1b Fill are sticky--use Subset to clear them' )
        endif
        if ( got(f_precision) ) then
          precisionQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(precisionVectorIndex), precisionQuantityIndex )
          call FromL1B ( key, quantity, chunks(chunkNo), &
            & filedatabase, isPrecision, suffix=suffix, geolocation=geolocation, &
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
          call FromL1B ( key, quantity, chunks(chunkNo), &
            & filedatabase, isPrecision, suffix=suffix, geolocation=geolocation, &
            & BOMask=BOMask )
        else
          call FromL1B ( key, quantity, chunks(chunkNo), &
            & filedatabase, isPrecision, suffix=suffix, geolocation=geolocation )
        end if

      case ( l_l2gp ) ! --------------  Fill from L2GP quantity  -----
        if ( .NOT. got(f_sourceL2GP) ) &
          & call Announce_Error ( key, noSourceL2GPGiven )
        call FromL2GP &
          & ( quantity, l2gpDatabase(l2gpIndex), interpolate, profile, errorCode, &
          & ignoreGeolocation, fromPrecision  )
        if ( errorCode /= 0 ) call Announce_error ( key, errorCode )

      case ( l_l2aux ) ! ------------  Fill from L2AUX quantity  -----
        if ( .NOT. got(f_sourceL2AUX) ) &
          & call Announce_Error ( key, noSourceL2AUXGiven )
        call FromL2AUX(quantity,l2auxDatabase(l2auxIndex),errorCode)
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
              call RHIFromOrToH2O ( key, quantity, &
                & sourceQuantity, temperatureQuantity, &
                & dontMask, ignoreZero, ignoreNegative, interpolate, &
                & .false., &   ! Mark Undefined values?
                & invert )    ! invert rather than convert?
            end if
          end if

      case ( l_H2OPrecisionfromRHi ) ! --fill H2O prec. from RHi quantity --
        if ( .not. any(got( &
          & (/f_h2oquantity, f_temperatureQuantity, &
          & f_rhiPrecisionquantity, f_tempPrecisionQuantity/) &
          & )) ) then
          call Announce_error ( key, No_Error_code, &
            & 'Missing a required field to fill h2o precision'  )
        else
          h2oQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(h2oVectorIndex), h2oQuantityIndex)
          temperatureQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(temperatureVectorIndex), temperatureQuantityIndex)
          rhiPrecisionQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(rhiPrecisionVectorIndex), rhiPrecisionQuantityIndex)
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
          else if ( .not. ValidateVectorQuantity(tempPrecisionQuantity, &
            & quantityType=(/l_temperature/)) ) then
            call Announce_Error ( key, No_Error_code, &
              & 'The tempPrecisionQuantity is not a temperature'  )
          else
            call RHIPrecisionFromOrToH2O ( key, quantity, &
              & rhiPrecisionQuantity, tempPrecisionQuantity, h2oQuantity, &
              & temperatureQuantity, dontMask, ignoreZero, &
              & ignoreNegative, interpolate, &
              & .true., &   ! Mark Undefined values?
              & invert=.true. )    ! convert RHiPrecisionToH2O
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
        call Hydrostatically ( key, quantity, temperatureQuantity, &
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
        call FromIsotope ( quantity, sourceQuantity, &
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
              call IWCFromExtinction ( quantity, &
                & sourceQuantity, temperatureQuantity)
            end if
          end if

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
          call LOSVelocity ( key, quantity, tngtECIQuantity, &
            & scECIquantity, scVelQuantity )
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
        ! We can make this simpler by imitating Van's use of pointers in
        ! place of optional args
        ! (Later, for now just see if you can make it work at all)
        if ( got(f_c) ) then
          if ( any( got( &
            & (/ f_height, f_channel /) &
            & ) ) ) then
            call CloneVectorQuantity( tempswapquantity, quantity )
            call ByManipulation ( tempswapquantity, aQuantity, bQuantity, &
              & manipulation, key, ignoreTemplate, &
              & spreadflag, dimList, &
              & c )
            call Explicit ( quantity, valuesNode, spreadFlag, force, &
              & vectors(vectorIndex)%globalUnit, channel, &
              & .false., options=heightRange(1:1), extraQuantity=tempswapquantity )
            call destroyVectorQuantityValue ( tempswapquantity, destroyMask=.true., &
              forWhom = 'tempswapquantity' )
          else
            call ByManipulation ( quantity, aQuantity, bQuantity, &
              & manipulation, key, ignoreTemplate, &
              & spreadflag, dimList, &
              & c )
          endif
        else
          if ( any( got( &
            & (/ f_height, f_channel /) &
            & ) ) ) then
            call CloneVectorQuantity( tempswapquantity, quantity )
            call ByManipulation ( tempswapquantity, aQuantity, bQuantity, &
              & manipulation, key, ignoreTemplate, &
              & spreadflag, dimList )
            call Explicit ( quantity, valuesNode, spreadFlag, force, &
              & vectors(vectorIndex)%globalUnit, channel, &
              & .false., options=heightRange(1:1), &
              & extraQuantity=tempswapquantity )
            call destroyVectorQuantityValue ( tempswapquantity, destroyMask=.true., &
              forWhom = 'tempswapquantity' )
          else
            call ByManipulation ( quantity, aQuantity, bQuantity, &
              & manipulation, key, ignoreTemplate, &
              & spreadflag, dimList )
          endif
        endif

      case ( l_magAzEl ) ! -- Magnetic Explicit from stren, azim, elev --
        if ( .not. got(f_explicitValues) ) &
          & call Announce_Error ( key, noExplicitValuesGiven )
        call Explicit ( quantity, valuesNode, spreadFlag, force, &
          & vectors(vectorIndex)%globalUnit, channel, &
          & azEl=.true. )

      case ( l_magneticModel ) ! --------------------- Magnetic Model --
        if ( .not. got ( f_gphQuantity ) ) then
          call Announce_Error ( key, no_Error_Code, 'Need gphQuantity for magnetic model' )
        else
          GPHQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(GPHVectorIndex), GPHQuantityIndex)
          call UsingMagneticModel ( quantity, gphQuantity, key )
        end if

      case ( l_modifyTemplate )
        shp = 0
        ! override qty template fileds: time, phi, longitude, etc.
        if ( .not. got ( f_manipulation ) ) &
          & call Announce_error ( key, no_Error_Code,'manipulation not supplied' )
        call get_string ( manipulation, GLStr, strip=.true. )
        if ( got( f_shape ) ) then
          do j = 2, nsons(shapeNode)
            call expr ( subtree ( j, shapeNode ), unitAsArray, valueAsArray )
            shp(j-1) = valueAsArray(1)
          enddo
        endif
        if ( got( f_sourceQuantity ) ) then
          sourceQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(sourceVectorIndex), sourceQuantityIndex )
          call ModifyQuantityTemplate  ( quantity%template, GLStr, &
            & sourceQuantity%values, spreadFlag )
        elseif ( got( f_explicitValues ) ) then
          if ( .not. got ( f_shape ) ) &
            & call Announce_error ( key, no_Error_Code,'shape not supplied' )
          call ModifyQuantityTemplate  ( quantity%template, GLStr, &
            & shp, valuesNode, spreadFlag )
        elseif ( .not. got( f_c ) ) then
          call Announce_error ( key, no_Error_Code, &
          & 'Must supply c for new value(s)' )
        else
          call ModifyQuantityTemplate  ( quantity%template, GLStr, c )
        endif

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

      case ( l_noRadsPerMIF )
        if ( .not. got ( f_measurements ) ) then
          call Announce_error ( key, No_Error_code, &
          & 'Missing a required field to fill noRadsPerMIF on MIFs'  )
        else
          measQty => GetVectorQtyByTemplateIndex( &
            & vectors(measVectorIndex), measQtyIndex)
          call NoRadsPerMif ( key, quantity, measQty, asPercentage )
        end if

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
          call fillTimings ( quantity%values(:,1), 'phases', 'all', .true. )
          ! call dump( quantity%values(:,1), 'phases' )
        end if

      case ( l_sectionTiming ) ! ---------  Fill timings for sections  -----
        call finishTimings('sections', returnStatus=status)
        if ( status /= 0 ) then
          call MLSMessage ( MLSMSG_Warning, ModuleName, 'Unable to finish sections timings' )
        else
          call fillTimings ( quantity%values(:,1), 'sections', 'all', .true. )
          ! call dump( quantity%values(:,1), 'sections' )
        end if

      case ( l_profile ) ! ------------------------ Profile fill -------
        if ( .not. got ( f_profileValues ) ) &
          call Announce_error ( key, no_Error_Code, 'profileValues not supplied' )
        if ( .not. got ( f_instances ) ) instancesNode = 0
        ! The next bit may seem odd - why not just pass logSpace on and let it default to false?
        ! The problem is, it should default to quantity%template%logSpace so absent means
        ! "don't care", not "logSpace=.false."
        if ( got ( f_logSpace ) ) then
          call FromProfile ( quantity, valuesNode, &
            & instancesNode, vectors(vectorIndex)%globalUnit, &
            & ptanQuantity, logSpace=logSpace )
        else
          call FromProfile ( quantity, valuesNode, &
            & instancesNode, vectors(vectorIndex)%globalUnit, &
            & ptanQuantity )
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
          call PhiTanWithRefraction ( key, quantity, h2oQuantity, &
            & orbitInclinationQuantity, &
            & ptanQuantity, refGPHquantity, temperatureQuantity, &
            & ignoreTemplate )
        end if

      case ( l_reflectorTempModel ) ! --------------- Reflector temperature model
        call WithReflectorTemperature ( key, quantity, phiZero, termsNode )

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
        call FromLosGrid ( key, Quantity, losQty, &
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
              call RHIFromOrToH2O ( key, quantity, &
                & h2oQuantity, temperatureQuantity, &
                & dontMask, ignoreZero, ignoreNegative, interpolate, &
                & .false., &   ! Mark Undefined values?
                & invert )    ! invert rather than convert?
            end if
        end if

      case ( l_quality )
        if ( .not. all ( got ( (/ f_sourceQuantity, f_scale /) ) ) ) &
          call Announce_error ( key, no_error_code, &
          & 'Need sourceQuanitity and scale for quality fill' )
        sourceQuantity => GetVectorQtyByTemplateIndex ( vectors(sourceVectorIndex), &
          & sourceQuantityIndex )
        call QualityFromChisq ( key, quantity, sourceQuantity, &
          & scale, heightNode, ignoreTemplate )

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
            call RHIPrecisionFromOrToH2O ( key, quantity, &
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
          call ScaleOverlaps ( quantity, multiplierNode )
        end if

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

        call FromSplitSideband ( quantity, sourceQuantity, &
          & lsbFraction, usbFraction, spreadFlag, usb, channel, key )

      case ( l_spreadChannel )
        if ( .not. got ( f_channel ) .and. .not. got( f_sourceQuantity ) ) &
          & call Announce_Error ( key, no_error_code, &
          & 'Must supply channel or sourcequantity for spreadChannel fill' )
        if ( got(f_sourceQuantity) ) then
        sourceQuantity => GetVectorQtyByTemplateIndex( &
          & vectors(sourceVectorIndex), sourceQuantityIndex )
          if ( .not. got(f_channel) ) channel = 1
          call SpreadChannelFill ( quantity, channel, key, &
            & sourceQuantity )
        else
          call SpreadChannelFill ( quantity, channel, key )
        endif

      case ( l_status )
        if ( got(f_ifMissingGMAO) ) then
          if ( MissingGMAO ) call Explicit ( quantity, &
            & valuesNode, .true., .false., phyq_Invalid, &
            & channel, .false., options=heightRange(1:1) )
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
          call StatusQuantity ( key, quantity, &
            & sourceQuantity, statusValue, &
            & minValue, maxValue, heightNode, &
            & additional, exact, ignoreTemplate )
        end if

      case ( l_swapvalues )
        if ( .not. got( f_sourceQuantity ) ) &
          & call Announce_Error ( key, no_error_code, &
          & 'Must supply sourcequantity for swapValues fill' )
        sourceQuantity => GetVectorQtyByTemplateIndex( &
          & vectors(sourceVectorIndex), sourceQuantityIndex )
        ! nullify( tempswapquantity )
        call CloneVectorQuantity( tempswapquantity, quantity )
        if ( quantity%template%name /= sourceQuantity%template%name ) then
          if ( .not. interpolate .and. .not. ignoreTemplate ) then
            call Announce_Error ( key, No_Error_Code, &
              & 'Quantity and sourceQuantity do not have the same template' )
          else
            call FromInterpolatedQty ( tempswapquantity, sourceQuantity, &
              & key, ignoreTemplate )
            call FromInterpolatedQty ( sourceQuantity, quantity, &
              & key, ignoreTemplate )
            call FromInterpolatedQty ( quantity, tempswapquantity, &
              & key, ignoreTemplate )
          end if
        else
          ! Just a straight copy
          ! If we have a mask and we're going to obey it then do so
          if ( associated(quantity%mask) ) then
            where ( iand ( ichar(quantity%mask(:,:)), m_Fill ) == 0 )
              tempswapquantity%values(:,:) = sourceQuantity%values(:,:)
              sourceQuantity%values(:,:) = Quantity%values(:,:)
              quantity%values(:,:) = tempswapquantity%values(:,:)
            end where
          else ! Otherwise, just blindly copy
            tempswapquantity%values = sourceQuantity%values
            sourceQuantity%values = Quantity%values
            quantity%values = tempswapquantity%values
          end if
        end if
        call destroyVectorQuantityValue ( tempswapquantity, destroyMask=.true., &
          forWhom = 'tempswapquantity' )
      case ( l_vector ) ! ---------------- Fill from another qty.
        if ( .not. got(f_sourceQuantity) ) &
          & call Announce_Error ( key, No_Error_Code, &
          & 'Missing a source field for vector fill' )
        Quantity => GetVectorQtyByTemplateIndex( &
          & vectors(VectorIndex), QuantityIndex )
        sourceQuantity => GetVectorQtyByTemplateIndex( &
          & vectors(sourceVectorIndex), sourceQuantityIndex )
        call FromAnother ( quantity, sourceQuantity, ptanQuantity, &
          & key, ignoreTemplate, spreadflag, interpolate, force )
      case ( l_uncompressRadiance )
        if ( .not. got(f_systemTemperature) ) &
          & call Announce_Error ( key, No_Error_Code, &
          & 'Missing a systemTemperature field for uncompress radiance fill' )
        if ( .not. got(f_totalPowerQuantity) ) &
          & call Announce_Error ( key, No_Error_Code, &
          & 'Missing a totalPowerQuantity field for uncompress radiance fill' )
        if ( .not. got(f_terms) ) &
          & call Announce_Error ( key, No_Error_Code, &
          & 'Missing a terms field for uncompress radiance fill' )

        sysTempQuantity => GetVectorQtyByTemplateIndex( &
          & vectors(sysTempVectorIndex), sysTempQuantityIndex )
        totalpowerQuantity => GetVectorQtyByTemplateIndex ( &
          & vectors(totalpowerVectorIndex), totalpowerQuantityIndex )
        call UncompressRadiance ( key, quantity, totalPowerQuantity, sysTempQuantity, termsNode )

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
          call WithWMOTropopause ( quantity, &
          & temperatureQuantity, refGPHQuantity, vGrids(internalVGridIndex) )
        else
          call WithReichlerWMOTP ( quantity, &
          & temperatureQuantity )
        end if
      case (-1)
        ! We must have decided to skip this fill
      case default
        call Announce_error ( key, no_Error_Code, 'This fill method not yet implemented' )
      end select      ! s_method
      ! Before leaving, we must restore the original mask (or lack of one)
      ! (unless the Fill method is l1b because we may be using the
      ! radiance precision to set the radiance mask)
      if ( fillMethod == l_l1b ) then
        ! For this Fill method the new mask is sticky
        qtyWasMasked = associated(quantity%mask)
      elseif ( qtyWasMasked ) then
        quantity%mask = tempQuantity%mask
      else
        if ( associated(quantity%mask) ) call ClearMask( quantity%mask )
      endif
      ! Housekeeping
      call destroyVectorQuantityValue ( tempQuantity, &
        & destroyMask=.true., destroyTemplate=.false. )
    end subroutine fillCommand

    ! ............................................  Get_File_Name  .....
    subroutine Get_File_Name ( fileIndex, pcfCode, &
      & fileType, fileDataBase, MLSFile, MSG, pcfEndCode )
      use FILLUTILS_1, only: ANNOUNCE_ERROR
      use HDF, only: DFACC_RDONLY
      use INTRINSIC, only: L_ASCII, L_HDF
      use MLSCOMMON, only: MLSFILE_T
      use MLSFILES, only: HDFVERSION_5, &
        & ADDINITIALIZEMLSFILE, GETPCFROMREF, SPLIT_PATH_NAME
      use MLSL2OPTIONS, only: TOOLKIT
      use SDPTOOLKIT, only: PGS_PC_GETREFERENCE
      use STRING_TABLE, only: GET_STRING
      ! Dummy args
      integer, intent(in) :: fileIndex ! index into string table
      integer, intent(in) :: pcfCode
      integer, intent(in) :: fileType ! l_hdf, etc.
      type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
      type (MLSFile_T), pointer   :: MLSFile
      character(len=*), intent(in) :: MSG ! in case of error
      integer, intent(in), optional :: pcfEndCode
      ! Internal variables
      character(len=255) :: fileName, fileTypeStr, PCFFileName, path, shortName
      integer :: lun
      integer :: mypcfEndCode
      integer :: returnStatus             ! non-zero means trouble
      integer :: version
      ! Executable
      mypcfEndCode = 0
      lun = 0
      version = 1
      call outputNamedValue( 'About to get_string for file name', fileIndex )
      call get_string ( fileIndex, shortName, strip=.true. )
      fileName = shortName
      call outputNamedValue( 'Result', fileName )
      call outputNamedValue( 'About to get_string for file type', lit_indices(fileType) )
      call get_string ( lit_indices(fileType), fileTypeStr, strip=.true. )
      call outputNamedValue( 'Result', fileTypeStr )
      if ( TOOLKIT .and. index (fileName, '/') < 1 ) then
        mypcfEndCode = pcfCode
        if ( present(pcfEndCode) ) mypcfEndCode = pcfEndCode
        if ( fileName == ' ' ) then
          returnStatus = Pgs_pc_getReference(pcfCode, version, &
            & fileName)
          lun = pcfCode
        else
          PCFFileName = fileName
          call split_path_name ( PCFFileName, path, fileName )
          lun = GetPCFromRef(fileName, pcfCode, &
            & mypcfEndCode, &
            & TOOLKIT, returnStatus, Version, DEEBUG, &
            & exactName=PCFFileName)
          if ( returnStatus /= 0 ) then
            call Announce_Error ( 0, son, extraMessage=MSG )
          else
            fileName = PCFFileName
          end if
        end if
      end if
      if ( fileType == l_hdf ) then
        MLSFile => AddInitializeMLSFile(filedatabase, &
          & content=fileTypeStr, &
          & name=Filename, shortName=shortName, &
          & type=l_hdf, access=dfacc_rdonly, HDFVersion=HDFVERSION_5)
      else
        MLSFile => AddInitializeMLSFile(filedatabase, &
          & content=fileTypeStr, &
          & name=Filename, shortName=shortName, &
          & type=l_ascii, access=dfacc_rdonly)
      endif      
      MLSFile%PCFId = lun
    end subroutine Get_File_Name

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

  end subroutine MLSL2Fill

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------
end module Fill
!=============================================================================

!
! $Log$
! Revision 2.427  2013/09/17 22:43:45  pwagner
! Added Scatter, Gather methods
!
! Revision 2.426  2013/09/17 00:52:38  vsnyder
! Correct 'no_code_for' error message
!
! Revision 2.425  2013/08/30 02:45:37  vsnyder
! Revise calls to trace_begin and trace_end
!
! Revision 2.424  2013/07/12 23:26:51  vsnyder
! Mistakenly commited a bunch of debugging stuff
!
! Revision 2.423  2013/07/12 23:25:28  vsnyder
! Remove unreferenced error messages
!
! Revision 2.422  2013/06/14 18:49:22  vsnyder
! Decruftification
!
! Revision 2.421  2013/05/31 00:42:12  vsnyder
! Add geolocation field to fill, used only if method=l1b
!
! Revision 2.420  2013/05/17 00:49:29  pwagner
! May constrain Transfer command command to quantitynames by r/t Boolean
!
! Revision 2.419  2013/04/24 00:37:42  pwagner
! Added InitRepeat and NextRepeat calls to set/increment r/t Boolean count
!
! Revision 2.418  2013/04/22 17:51:32  pwagner
! Reevaluate may store a literal instead of a Boolean value
!
! Revision 2.417  2013/04/17 00:05:07  pwagner
! Added new Repeat control structure to Fill, Retrieve sections
!
! Revision 2.416  2013/01/17 20:02:02  pwagner
! New where field in Subset command
!
! Revision 2.415  2013/01/02 21:40:59  pwagner
! Added derivative method to Fill command; Transfer can do Fill methods, too
!
! Revision 2.414  2012/11/14 20:04:11  pwagner
! Fix boxcar bug added with last commit
!
! Revision 2.413  2012/11/14 00:58:42  pwagner
! Use dimList for choosing which of {csi} to average over; finished treating width=1
!
! Revision 2.412  2012/11/08 23:21:30  pwagner
! dimList field lets us specifiy whether to shift by [c,s,i] during manipulate
!
! Revision 2.411  2012/11/05 19:02:48  pwagner
! Fixed various bugs related to last changes
!
! Revision 2.410  2012/10/31 00:08:40  pwagner
! Must treat Fill method L1B specially when restoring quantity mask
!
! Revision 2.409  2012/10/29 17:18:24  pwagner
! Made consistent with FillUtils_1 api
!
! Revision 2.408  2012/10/27 00:27:29  pwagner
! Temporary fixes; cause of error not yet understood
!
! Revision 2.407  2012/10/22 18:14:50  pwagner
! Many Subset operations now available in Fill
!
! Revision 2.406  2012/10/09 00:48:30  pwagner
! New ignoreTemplate, changed force meaning in Fill
!
! Revision 2.405  2012/08/16 17:57:16  pwagner
! Exploit level 2-savvy MLSMessage
!
! Revision 2.404  2012/07/31 00:48:16  vsnyder
! Use DestroyVectorQuantityValue abstraction
!
! Revision 2.403  2012/05/08 17:49:04  pwagner
! Added Select .. Case .. EndSelect control structure
!
! Revision 2.402  2012/02/24 21:20:44  pwagner
! DirectRead may /interpolate vertically
!
! Revision 2.401  2012/02/13 23:28:45  pwagner
! Corrected error message
!
! Revision 2.400  2012/02/10 23:45:43  vsnyder
! Add more tracing, change some tracing levels
!
! Revision 2.399  2012/02/02 01:19:05  pwagner
! Can DirectRead matrix or hessian
!
! Revision 2.398  2012/01/25 01:19:01  pwagner
! Removed unused args to ..FromFile calls
!
! Revision 2.397  2012/01/18 20:38:59  vsnyder
! Check consistency of covariance matrix and diagonal vector
!
! Revision 2.396  2011/12/15 01:49:43  pwagner
! Added sdName and /spread fields to DirectRead
!
! Revision 2.395  2011/11/04 00:28:18  pwagner
! Added autoFill flag to Vector Spec to start it off with non-zero values in appropriated quantities
!
! Revision 2.394  2011/10/07 00:06:02  pwagner
! May dump Matrices, Hessians from Fill, Join
!
! Revision 2.393  2011/06/16 20:52:22  vsnyder
! Get codes for Announce_Error from FillUtils.  Add f_expr, but with a
! noCodeFor error message -- to be implemented later.
!
! Revision 2.392  2011/04/20 16:46:37  pwagner
! Removed unused declaration
!
! Revision 2.391  2011/03/22 23:47:54  pwagner
! May now reshape qty template field while filling explicitly
!
! Revision 2.390  2011/03/15 22:51:58  pwagner
! May now modify quantity template fields with fill method
!
! Revision 2.389  2010/11/20 00:01:14  pwagner
! May specifiy surfaces gap beyond which to zero out in Streamline
!
! Revision 2.388  2010/08/06 23:08:48  pwagner
! Pass Hessians, matrices to DumpCommand
!
! Revision 2.387  2010/07/22 17:42:24  pwagner
! Replaced method=special fills with unique names
!
! Revision 2.386  2010/07/06 16:06:06  pwagner
! Better error checking in Transfer
!
! Revision 2.385  2010/07/01 00:49:19  pwagner
! Transfer between vectors may now also manipulate
!
! Revision 2.384  2010/05/19 23:06:45  pwagner
! Shorten most Fill routine names
!
! Revision 2.383  2010/05/19 17:53:19  pwagner
! Removed unused stuff
!
! Revision 2.382  2010/04/28 16:24:11  pwagner
! May specify instances range in explicit Fill
!
! Revision 2.381  2010/04/22 23:36:46  pwagner
! May fill num rads/MIF as a percentage
!
! Revision 2.380  2010/04/13 01:43:09  vsnyder
! Move FlushLockedBins from LinearizedForwardModel_m to L2PCBins_m
!
! Revision 2.379  2010/03/26 23:16:56  vsnyder
! Add Threshold to StreamlineHessian
!
! Revision 2.378  2010/03/25 01:50:25  vsnyder
! Make sure GeodAngle and ScaleHeight get values in StreamlineHessian
!
! Revision 2.377  2010/02/25 18:37:51  pwagner
! Adds support for new Hessian data type
!
! Revision 2.376  2009/10/26 17:11:28  pwagner
! Added Diff command to be used like Dump in l2cf
!
! Revision 2.375  2009/08/24 20:13:47  pwagner
! May Fill H2O precision from RHI precision
!
! Revision 2.374  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.373  2009/04/29 23:12:54  pwagner
! Manipulation Fills can be restricted by height and heightRange
!
! Revision 2.372  2009/04/28 20:02:50  pwagner
! No longer sets undefined values in RHI to -999.99
!
! Revision 2.371  2009/04/16 21:55:05  pwagner
! /exact keyword in status Fill to fix radiance bug
!
! Revision 2.370  2009/04/13 20:45:57  pwagner
! heightRange in explicit Fill can fill above or below specified height
!
! Revision 2.369  2009/03/05 18:37:59  pwagner
! May specifiy height, channel with explicit Fill
!
! Revision 2.368  2008/12/18 21:14:05  pwagner
! May now dump an l2pc or allL2PCs (use with caution)
!
! Revision 2.367  2008/11/06 21:50:46  pwagner
! Fill method swapValues swaps values between two quantities
!
! Revision 2.366  2008/09/24 16:46:37  livesey
! Tidy up handling of ptan in profile fill
!
! Revision 2.365  2008/09/19 23:55:05  pwagner
! May now Destroy GriddedData
!
! Revision 2.364  2008/09/16 22:29:05  pwagner
! pass optional field ptanQuantity to profile, vector methods
!
! Revision 2.363  2008/08/14 20:59:00  pwagner
! /interpolate now possible field in Transfer command
!
! Revision 2.362  2008/05/28 21:52:48  pwagner
! geo location Fill method to fill chunk numbers[maf]
!
! Revision 2.361  2008/04/26 00:40:07  livesey
! Added total power stuff
!
! Revision 2.360  2008/04/11 01:17:22  livesey
! Added uncompressRadiance fill
!
! Revision 2.359  2007/12/07 01:13:39  pwagner
! Lets us catch warnings and assign to runtime Booleans
!
! Revision 2.358  2007/11/15 22:53:16  pwagner
! May set runtimeBooleans by anyGood.., Compare, Reevaluate commands
!
! Revision 2.357  2007/11/05 18:38:09  pwagner
! May Skip remaining lines in Fill, Join, Retrieve sections depending on Boolean
!
! Revision 2.356  2007/09/27 22:00:45  pwagner
! Much moved into new FillUtils_1 module; Intel compiler can now optimize
!
! Revision 2.355  2007/08/27 23:56:34  pwagner
! /spread flag now affects statistical manipulation Fills
!
! Revision 2.354  2007/08/23 22:17:05  pwagner
! manipulation Fills can now use statistical functions
!
! Revision 2.353  2007/08/20 22:04:48  pwagner
! Many procedures now push their names onto MLSCallStack
!
! Revision 2.352  2007/07/06 17:09:16  pwagner
! Reverted to former ApplyBaseline
!
! Revision 2.351  2007/06/21 22:34:32  pwagner
! Fixed inconsequential bug in adding baseline when quantity was masked
!
! Revision 2.350  2007/03/23 00:26:09  pwagner
! More unused debugging; skips filling covariance matrix if skipping retrievals
!
! Revision 2.349  2007/01/12 00:34:04  pwagner
! Renamed routine outputNamedValue
!
! Revision 2.348  2006/11/03 19:40:17  pwagner
! Fixed unassociated pointers NAG caught
!
! Revision 2.347  2006/11/03 00:26:07  pwagner
! Fixed bug in tropopause calculation
!
! Revision 2.346  2006/10/11 22:56:11  pwagner
! May fill convergence from dnwt_chisqRatio
!
! Revision 2.345  2006/10/03 20:24:11  pwagner
! Optional test, tweaks to ChiSqRatio
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
! Use HFTI instead of Cholesky for UsingLeastSquares, for stability
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
! Use of single arg options in Explicit replaces three
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
! Added forgiveZeros handling to Covariance for efficiency.
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
! Removed buggy, unused ColAbundance
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
! Moved WithCombinedChannels into ManipulateVectorQuantitites
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
! Made Covariance (temporarily?) fill both sides of the digaonal (in
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
! add new flag ignoreGeolocation in subroutine FromL2GP
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
! Cosmetic and superficial changes to FromSplitSideband
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
! Added a bit more intelligence to Covariance
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
! Always sets errorCode to 0 in return from FromL2AUX
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
! Better handling of missing length scale in Covariance
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
! Made the checking in ByManipulation a little more lenient
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
! Tidied up FromProfile
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
! RHI testing begun; incomplete
!
! Revision 2.118  2002/04/13 00:31:46  pwagner
! More flesh on FillrhiFromH2o; still untested
!
! Revision 2.117  2002/04/11 23:51:28  pwagner
! Fleshed out RHIFromH2O; untested yet
!
! Revision 2.116  2002/04/10 17:45:44  pwagner
! Added RHI from h2oquantity (just a placeholder)
!
! Revision 2.115  2002/04/04 16:32:42  livesey
! Added negative error bar stuff
!
! Revision 2.114  2002/03/27 17:37:57  livesey
! Minor changes to random number stuff.
! Now seed incremented with chunk number
!
! Revision 2.113  2002/03/19 00:52:40  pwagner
! Some new checks added to FillLOSVelocity
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
! add Diagonal
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
! some fixes for FillFromLOS
!
! Revision 2.60  2001/07/19 18:05:42  dwu
! add sourceSGRID
!
! Revision 2.59  2001/07/19 00:56:27  dwu
! fix bugs in FillFromLos
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
! Completed FillL2AUXData; changed squeeze, nearby
!
! Revision 2.7  2000/12/05 00:40:50  pwagner
! Added FillL2AUXVector
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
