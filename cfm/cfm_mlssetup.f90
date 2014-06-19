! Copyright 2010, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.
module CFM_MLSSetup_m

    use ForwardModelConfig, only: ForwardModelConfig_T
    use QuantityTemplates, only: QuantityTemplate_T, Dump
    use LEXER_CORE, only: INIT_LEXER
    use DECLARATION_TABLE, only: ALLOCATE_DECL, DEALLOCATE_DECL
    use TREE, only: ALLOCATE_TREE, DEALLOCATE_TREE
    use H5LIB, ONLY: h5open_f, h5close_f
    use MLSMessageModule, only: MLSMessageConfig, &
                                MLSMessage, MLSMSG_Error
    use SYMBOL_TABLE, only: DESTROY_SYMBOL_TABLE
    use MLSCommon, only: MLSFile_T, TAI93_Range_T, r8
    use Intrinsic, only: l_hdf
    use INIT_TABLES_MODULE, only: phyq_pressure, l_logarithmic, l_zeta, &
                                  phyq_vmr, l_vmr, l_earthRefl, l_losVel, &
                                  l_scgeocalt, l_spaceradiance, l_o2, &
                                  l_elevOffset, l_limbsidebandFraction
    use ConstructQuantityTemplates, only: InitQuantityTemplates
    use ChunkDivide_m, only: ChunkDivideConfig_T, &
                             GetChunkFromTimeRange => CFM_ChunkDivide
    use Chunks_m, only: MLSChunk_T ! To also be referenced outside
    use VGridsDatabase, only: VGrid_T
    use HGridsDatabase, only: HGrid_T
    use VectorsModule, only: GetVectorQtyByTemplateIndex, & ! being used in many functions
                             VectorTemplate_T, Dump, Vector_T, VectorValue_T
    use CFM_Vector_m, only: CreateValue4AgileVector, AddValue2Vector
    use CFM_Constants_m
    use CFM_QuantityTemplate_m, only: CreateQtyTemplate
    use CFM_Fill_M, only: FillVectorQtyFromProfile, ExplicitFillVectorQuantity, &
        FillVectorQuantityFromL1B

    implicit none

    private
    public :: CFM_MLSSetup, CFM_MLSCleanup, CreateMLSValue_O2
    public :: CreateMLSValue_EarthReflectivity, CreateMLSValue_LSF
    public :: CreateMLSValue_FromL1BOA, CreateMLSValue_SpaceRadiance
    public :: GetConstantQuantities, CreateMLSValue_ElevationOffset
    public :: timeRange2MafRange

    interface CFM_MLSSetup
        module procedure CFM_MLSSetup_Obsolete, CFM_MLSSetup_Compact
    end interface

    interface CFM_MLSCleanup
        module procedure CFM_MLSCleanup_Obsolete, CFM_MLSCleanup_Compact
    end interface

!---------------------------- RCS Ident Info -------------------------------
    character(len=*), private, parameter :: ModuleName= &
        "$RCSfile$"
    private :: not_used_here
!---------------------------------------------------------------------------

    type(ChunkDivideConfig_T), save :: chunkDivideConfig  ! Using default options

    contains

    ! Initializes CFM library and HDF5 library. Also parses the forward model
    ! configuration file.
    subroutine CFM_MLSSetup_Compact (signalFileName, configFileName, ForwardModelConfigDatabase)
        use CFM_Tree_Walker_m, only : Walk_Tree
        use EmpiricalGeometry, only: CFM_InitEmpiricalGeometry
        use string_table, only: AddInUnit
        use io_stuff, only: get_lun
        use INIT_TABLES_MODULE, only: INIT_TABLES
        use Parser, only: Clean_Up_Parser, Configuration
        use Parser_Table_m, only:  Destroy_Parser_Table, Parser_Table_t
        use Parser_Tables_L2CF, only: Init_Parser_Table
        use TREE_CHECKER, only: CHECK_TREE

        ! The name of the signal database file
        character(len=*), intent(in) :: signalFileName
        ! The name of the forward model configuration file
        character(len=*), intent(in) :: configFileName
        ! Output: to store all ForwardModelConfig_T objects declared in
        ! the config file.
        ! This argument will be nullified inside this subroutine.
        type (ForwardModelConfig_T), pointer :: ForwardModelConfigDatabase(:)
        type(Parser_Table_t) :: Parser_Table

        integer :: Root
        integer :: First_Section
        integer :: error
        integer :: signalIn, configIn

        ! Set this, so we don't have the missing log error
        MLSMessageConfig%useToolkit = .false.
        MLSMessageConfig%logFileUnit = -1

        ! Set up empirical geometry before reading the L2CFs, which include
        ! global_setting, where vGrids can be defined.
        call CFM_InitEmpiricalGeometry (empiricalGeometry_noIterations, &
                                        empiricalGeometry_terms)

        call get_lun(signalIn)
        open (unit=signalIn, file=trim(signalFileName), status='OLD', &
        iostat=error)
        if (error /= 0) &
            call MLSMessage (MLSMSG_Error, moduleName, &
            'Error opening ' // trim(signalFileName))

        call get_lun(configIn)
        open(unit=configIn, file=trim(configFileName), status='OLD', &
        iostat=error)
        if (error /= 0) &
            call MLSMessage (MLSMSG_Error, moduleName, &
            'Error opening ' // trim(configFileName))

        call init_lexer ( n_chars=80000, n_symbols=4000, hash_table_size=611957 )
        call allocate_decl ( ndecls=8000 )
        call allocate_tree ( n_tree=2000000 )
        call init_tables

        call AddInUnit(signalIn)
        call AddInUnit(configIn)

        ! Parse the L2CF, producing an abstract syntax tree
        call init_parser_table ( parser_table )
        call configuration ( root, parser_table )
        call destroy_parser_table ( parser_table )
        call clean_up_parser
        if (Root <= 0) then
            call MLSMessage (MLSMSG_Error, moduleName, &
            'A syntax error occurred -- there is no abstract syntax tree')
        end if

        call check_tree ( root, error, first_section )
        if (error /= 0) &
            call MLSMessage (MLSMSG_Error, moduleName, "Error in tree")

        nullify(ForwardModelConfigDatabase)
        call Walk_Tree ( Root, First_Section, &
            ForwardModelConfigDatabase=ForwardModelConfigDatabase )

        close (signalIn, iostat=error)
        if (error /= 0) &
            call MLSMessage (MLSMSG_Error, moduleName, &
            "Error closing " // trim(signalFileName))

        close (configIn, iostat=error)
        if (error /= 0) &
            call MLSMessage (MLSMSG_Error, moduleName, &
           'Error closing ' // trim(configFileName))

        ! We have to call this before opening any HDF5 file
        call h5open_f(error)
        if (error /= 0) &
            call MLSMessage (MLSMSG_Error, moduleName, "Error in initialize hdf5 library")

        ! Have to initialize before we start creating quantity templates
        call InitQuantityTemplates
    end subroutine

    ! Read the signal file and populate the signal database
    ! Read the config file for forward model configuration(s),
    ! then put them in ForwardModelConfigDatabase.
    ! Reads L1BOA file and put it into the filedatabase.
    ! Uses the startTime, endTime and leapsec file to create a
    ! MLSChunk_T object that uses by other subroutines that read
    ! data from L1BOA.
    ! Construct the quantity template databases, which will include
    ! the quantities inside stateVectorExtra
    ! Create stateVectorExtra that is the second input to the forward
    ! model.
    ! Note: the use of a filedatabase, instead of just one MLSFile_T
    ! object to store L1BOA is for convenience when wrapping
    ! existing code in MLSPGS.
    ! CFM_MLSSetup is to be called only once before CFM_MLSCleanup is called.
    subroutine CFM_MLSSetup_Obsolete (startTime, endTime, l1boa, leapsecFile, signalFileName, &
        configFileName, filedatabase, qtyTemplates, chunk, ForwardModelConfigDatabase, &
        stateVectorExtra)
        use CFM_Tree_Walker_m, only : Walk_Tree
        use VectorsModule, only: Vector_T
        use string_table, only: AddInUnit
        use MLSFiles, only: AddFileToDatabase, InitializeMLSFile, mls_openFile
        use io_stuff, only: get_lun
        use SDPToolkit, only: mls_utctotai
        use EmpiricalGeometry, only: CFM_InitEmpiricalGeometry
        use Hdf, only: DFACC_RDONLY
        use INIT_TABLES_MODULE, only: INIT_TABLES, l_ghz, phyq_mafs, l_orbital
        use Parser, only: Clean_Up_Parser, Configuration
        use Parser_Table_m, only:  Destroy_Parser_Table, Parser_Table_t
        use Parser_Tables_L2CF, only: Init_Parser_Table
        use TREE_CHECKER, only: CHECK_TREE

        ! The start time of the data to be read in the format
        ! yyyy-doyThh:mm:ss.zzzz
        character(len=CCSDSlen), intent(in) :: startTime
        ! The stop time of the data to be read in the same format as startTime.
        character(len=CCSDSlen), intent(in) :: endTime
        ! The file name of L1BOA file.
        character(len=*), intent(in) :: l1boa
        ! The name of leapsec file
        character(len=*), intent(in) :: leapsecFile
        ! The name of the signal database file
        character(len=*), intent(in) :: signalFileName
        ! The name of the forward model configuration file
        character(len=*), intent(in) :: configFileName
        ! Output: An array of MLSFile_T object, representing open file(s),
        ! including L1BOA.
        type (MLSFile_T), dimension(:), pointer :: filedatabase
        ! The quantity template database, which is nullified by this
        ! subroutine.
        type(QuantityTemplate_T), dimension(:), pointer :: qtyTemplates
        ! Output: A data holder that holds the startTime and endTime
        ! in a way that can be understood by other subroutines.
        type (MLSChunk_T), intent(out) :: chunk
        ! Output: to store all ForwardModelConfig_T objects declared in
        ! the config file.
        ! This argument will be nullified inside this subroutine.
        type (ForwardModelConfig_T), pointer :: ForwardModelConfigDatabase(:)
        type(Parser_Table_t) :: Parser_Table
        ! A vector filled with quantities that can be supplied by MLS
        ! and needed for the forward model.
        type (Vector_T), intent(out) :: stateVectorExtra

        integer :: Root
        integer :: First_Section
        integer :: error
        type (MLSFile_T), target :: l1bfile
        type (TAI93_Range_T) :: processingRange
        integer :: signalIn, configIn

        !Executables
        nullify(filedatabase)

        ! Set this, so we don't have the missing log error
        MLSMessageConfig%useToolkit = .false.
        MLSMessageConfig%logFileUnit = -1

        ! Set up empirical geometry before reading the L2CFs, which include
        ! global_setting, where vGrids can be defined.
        call CFM_InitEmpiricalGeometry (empiricalGeometry_noIterations, &
                                        empiricalGeometry_terms)

        call get_lun(signalIn)
        open (unit=signalIn, file=trim(signalFileName), status='OLD', &
        iostat=error)
        if (error /= 0) &
            call MLSMessage (MLSMSG_Error, moduleName, &
            'Error opening ' // trim(signalFileName))

        call get_lun(configIn)
        open(unit=configIn, file=trim(configFileName), status='OLD', &
        iostat=error)
        if (error /= 0) &
            call MLSMessage (MLSMSG_Error, moduleName, &
            'Error opening ' // trim(configFileName))

        call init_lexer ( n_chars=80000, n_symbols=4000, hash_table_size=611957 )
        call allocate_decl ( ndecls=8000 )
        call allocate_tree ( n_tree=2000000 )
        call init_tables

        call AddInUnit(signalIn)
        call AddInUnit(configIn)

        ! Parse the L2CF, producing an abstract syntax tree
        call init_parser_table ( parser_table )
        call configuration ( root, parser_table )
        call destroy_parser_table ( parser_table )
        call clean_up_parser
        if (Root <= 0) then
            call MLSMessage (MLSMSG_Error, moduleName, &
            'A syntax error occurred -- there is no abstract syntax tree')
        end if

        call check_tree ( root, error, first_section )
        if (error /= 0) &
            call MLSMessage (MLSMSG_Error, moduleName, "Error in tree")

        nullify(ForwardModelConfigDatabase)
        call Walk_Tree ( Root, First_Section, &
            ForwardModelConfigDatabase=ForwardModelConfigDatabase )

        close (signalIn, iostat=error)
        if (error /= 0) &
            call MLSMessage (MLSMSG_Error, moduleName, &
            "Error closing " // trim(signalFileName))

        close (configIn, iostat=error)
        if (error /= 0) &
            call MLSMessage (MLSMSG_Error, moduleName, &
            'Error closing ' // trim(configFileName))

        ! We have to call this before opening any HDF5 file
        call h5open_f(error)
        if (error /= 0) &
            call MLSMessage (MLSMSG_Error, moduleName, "Error in initialize hdf5 library")

        error = InitializeMLSFile(l1bfile, content='l1boa', &
        name=trim(l1boa), shortName='L1BOA', type=l_hdf, access=DFACC_RDONLY)
        if (error /= 0) &
            call MLSMessage (MLSMSG_Error, moduleName, &
            "Error initializing " // trim(l1boa))

        call mls_openFile(L1BFile, error)
        if (error /= 0 ) &
            call MLSMessage (MLSMSG_Error, moduleName, &
            "Error opening " // trim(l1boa))

        ! Just something to capture the output of the function
        ! in order to not get compiling error
        error = AddFileToDatabase(filedatabase, l1bfile)

        ! Get the start time and end time encoded
        error = mls_utctotai (trim(leapsecFile), startTime, processingrange%starttime)
        if (error /= 0) &
            call MLSMessage (MLSMSG_Error, moduleName, &
            "Could not convert UTC Start time to TAI")

        error = mls_utctotai (trim(leapsecFile), endtime, processingrange%endtime)
        if (error /= 0) &
            call MLSMessage (MLSMSG_Error, moduleName, &
            "Could not convert UTC End time to TAI")

        ! Initialize ChunkDivideConfig
        chunkDivideConfig%method = l_orbital
        chunkDivideConfig%maxLength = 3500
        chunkDivideConfig%maxLengthFamily = phyq_mafs
        chunkDivideConfig%skipL1BCheck = .true.
        chunkDivideConfig%homeModule = l_ghz
        chunkDivideConfig%criticalModules = l_ghz
        chunkDivideConfig%homeGeodAngle = 0.0_r8
        chunkDivideConfig%maxGap = 2.0_r8
        chunkDivideConfig%maxGapFamily = phyq_mafs
        chunkDivideConfig%maxOrbY = 50000.0_r8 ! Unit is meter
        chunkDivideConfig%scanLowerLimit = (/-20000.0_r8, 10000.0_r8 /) ! unit is meter
        chunkDivideConfig%scanUpperLimit = (/ 40000.0_r8, 200000.0_r8 /) ! unit is meter
        chunkDivideConfig%lowerOverlap = 0.0
        chunkDivideConfig%upperOverlap = 0.0
        chunkDivideConfig%lowerOverlapFamily = phyq_mafs
        chunkDivideConfig%upperOverlapFamily = phyq_mafs
        chunkDivideConfig%noChunks = 1
        ! Create a fake chunk out of the start time, end time, and L1BOA
        chunk = GetChunkFromTimeRange(processingRange, filedatabase, chunkDivideConfig)

        ! Have to initialize before we start creating quantity templates
        call InitQuantityTemplates

        call CreateStateVectorExtra(filedatabase, chunk, qtyTemplates, stateVectorExtra)

    end subroutine

    ! Clean up allocated memory, close opened files.
    ! This subroutine will deallocate the forward model config database.
    subroutine CFM_MLSCleanup_Compact (forwardModelConfigDatabase)
        use ForwardModelConfig, only: DestroyFWMConfigDatabase
        use STRING_TABLE, only: DESTROY_CHAR_TABLE, DESTROY_HASH_TABLE, &
                                DESTROY_STRING_TABLE
        use MLSSignals_m, only: modules, bands, radiometers, spectrometerTypes, &
                                signals, DestroyBandDatabase, DestroySignalDatabase, &
                                DestroySpectrometerTypeDatabase, DestroyModuleDatabase, &
                                DestroyRadiometerDatabase
        use EmpiricalGeometry, only: CFM_ResetEmpiricalGeometry

        ! The forward model configurations created by CFM_MLSSetup
        type (ForwardModelConfig_T), pointer :: ForwardModelConfigDatabase(:)

        integer :: error

        call DestroyFWMConfigDatabase(forwardModelConfigDatabase)

        ! Clean up for the tree
        call destroy_char_table
        call destroy_hash_table
        call destroy_string_table
        call destroy_symbol_table
        call deallocate_decl
        call deallocate_tree

        ! Destroy signal and related database
        call DestroyBandDatabase(bands)
        call DestroySpectrometerTypeDatabase(spectrometerTypes)
        call DestroySignalDatabase(signals)
        call DestroyRadiometerDatabase(radiometers)
        call DestroyModuleDatabase(modules)

        ! Reset EmpiricalGeometry because something is allocated
        call CFM_ResetEmpiricalGeometry

        ! Have to call this, after we stops using HDF5 library
        call h5close_f (error)
        if (error /= 0) then
            print *, "Error in finishing up hdf5 library"
            return
        end if
    end subroutine

    ! Clean up allocated memory, close opened files
    subroutine CFM_MLSCleanup_Obsolete (filedatabase, qtyTemplates, forwardModelConfigDatabase, &
                                        stateVectorExtra)
        use ForwardModelConfig, only: DestroyFWMConfigDatabase
        use MLSFiles, only: mls_closeFile
        use STRING_TABLE, only: DESTROY_CHAR_TABLE, DESTROY_HASH_TABLE, &
                                DESTROY_STRING_TABLE
        use QuantityTemplates, only: DestroyQuantityTemplateDatabase
        use VectorsModule, only: Vector_T, DestroyVectorInfo
        use VectorsModule, only: DestroyVectorTemplateInfo
        use MLSSignals_m, only: modules, bands, radiometers, spectrometerTypes, &
                                signals, DestroyBandDatabase, DestroySignalDatabase, &
                                DestroySpectrometerTypeDatabase, DestroyModuleDatabase, &
                                DestroyRadiometerDatabase
        use EmpiricalGeometry, only: CFM_ResetEmpiricalGeometry

        ! The filedatabase created by CFM_MLSSetup
        type (MLSFile_T), dimension(:), pointer :: filedatabase
        ! The quantity template database created by CFM_MLSSetup
        type(QuantityTemplate_T), dimension(:), pointer :: qtyTemplates
        ! The forward model configurations created by CFM_MLSSetup
        type (ForwardModelConfig_T), pointer :: ForwardModelConfigDatabase(:)
        ! The vector created by CFM_MLSSetup
        type(Vector_T) :: stateVectorExtra

        integer :: error
        integer :: i

        call DestroyFWMConfigDatabase(forwardModelConfigDatabase)

        ! Clean up for the tree
        call destroy_char_table
        call destroy_hash_table
        call destroy_string_table
        call destroy_symbol_table
        call deallocate_decl
        call deallocate_tree

        do i = 1, size(filedatabase)
            call mls_closefile(filedatabase(i))
        end do

        deallocate(filedatabase)
        nullify(filedatabase)

        call DestroyVectorTemplateInfo(stateVectorExtra%template)
        call DestroyVectorInfo(stateVectorExtra)
        call DestroyQuantityTemplateDatabase (qtyTemplates)

        ! Destroy signal and related database
        call DestroyBandDatabase(bands)
        call DestroySpectrometerTypeDatabase(spectrometerTypes)
        call DestroySignalDatabase(signals)
        call DestroyRadiometerDatabase(radiometers)
        call DestroyModuleDatabase(modules)

        ! Reset EmpiricalGeometry because something is allocated
        call CFM_ResetEmpiricalGeometry

        ! Have to call this, after we stops using HDF5 library
        call h5close_f (error)
        if (error /= 0) then
            print *, "Error in finishing up hdf5 library"
            return
        end if

    end subroutine CFM_MLSCleanup_Obsolete

    ! Creates a VectorValue_T object and fills it with constants values for O2.
    ! The constants are declared in cfm_constants.f90.
    type(VectorValue_T) function CreateMLSValue_O2 (avgrid, ahgrid, qname) result(o2)
        type(VGrid_T), intent(in) :: avgrid
        type(HGrid_T), intent(in) :: ahgrid
        character(len=*), intent(in), optional :: qname

        type(QuantityTemplate_T) :: o2template

        o2template = CreateQtyTemplate(l_vmr, avgrid=avgrid, ahgrid=ahgrid, qMolecule=l_o2, qname=qname)

        o2 = CreateValue4AgileVector(o2template)
        call FillVectorQtyFromProfile (o2, .false., o2_heights, o2_values, phyq_vmr)
    end function

    ! Creates a VectorValue_T object and fills it with values for earth 
    ! reflectivity. The values are declared in cfm_constants.f90.
    type(VectorValue_T) function CreateMLSValue_EarthReflectivity () result(earthRefl)
        type(QuantityTemplate_T) :: ertemplate

        ertemplate = CreateQtyTemplate(l_earthRefl)
        earthRefl = CreateValue4AgileVector(ertemplate, value=earthRefl_values)
    end function

    ! Creates a VectorValue_T object and fills it with values for space
    ! radiance. The values are declared in cfm_constants.f90.
    type(VectorValue_T) function CreateMLSValue_SpaceRadiance () result(spaceRad)
        type(QuantityTemplate_T) :: srtemplate

        srtemplate = CreateQtyTemplate(l_spaceradiance)
        spaceRad = CreateValue4AgileVector(srtemplate, value=spaceRad_values)
    end function

    ! Creates a VectorValue_T object and fills it with values from an L1BOA file.
    ! This is possible because every quantity in L1BOA is of a different type.
    type(VectorValue_T) function CreateMLSValue_FromL1BOA (qType, instrumentModule, filedatabase, &
    firstL1Maf, lastL1Maf, qname) result (vv)
        ! The type of the quantity to be read. This type is used to 
        ! identify the quantity in the L1BOA file.
        integer, intent(in) :: qType
        ! One of "GHz", "THz" and "sc"
        character(len=*), intent(in) :: instrumentModule
        ! The filedatabase with L1BOA already opened.
        type (MLSFile_T), dimension(:), pointer :: filedatabase
        ! The start L1 Maf of the value range to be read.
        integer, intent(in) :: firstL1Maf
        ! The end L1 Maf of the value range to be read.
        integer, intent(in) :: lastL1Maf
        ! Optional name to give to the quantity. 
        character(len=*), intent(in), optional :: qname

        type(MLSChunk_T) :: chunk
        type(QuantityTemplate_T) :: template

        chunk%firstMafIndex = firstL1Maf
        chunk%lastMafIndex = lastL1Maf

        template = CreateQtyTemplate(qType, filedatabase=filedatabase, chunk=chunk, &
                                     qInstModule=instrumentModule, qname=qname)
        vv = CreateValue4AgileVector(template)
        call FillVectorQuantityFromL1B (vv, chunk, filedatabase, .false.)
    end function

    ! Creates a VectorValue_T object for elevation offset of a given upper/lower band
    ! and gives it an optional name. This quantity's values are MLS-specified constants
    ! which are declared in cfm_constants.f90.
    type(VectorValue_T) function CreateMLSValue_ElevationOffset (band, upper, qname) result (vv)
        integer, intent(in) :: band
        logical, intent(in) :: upper ! if upper is false, then it's lower
        character(len=*), intent(in), optional :: qname

        type(QuantityTemplate_T) :: template
        character(len=32) :: signal
        real(r8) :: spreadvalue

        select case(band)
        case (1)
            if (upper) then
                call MLSMessage(MLSMSG_Error, moduleName, "Band 1U does not exist.")
            else
                signal = band1L
                spreadvalue = velev1L
            endif
        case (2)
            if (upper) then
                signal = band2U
                spreadvalue = velev2U
            else
                signal = band2L
                spreadvalue = velev2L
            endif
        case (3)
            if (upper) then
                signal = band3U
                spreadvalue = velev3U
            else
                signal = band3L
                spreadvalue = velev3L
            endif
        case (4)
            if (upper) then
                signal = band4U
                spreadvalue = velev4U
            else
                signal = band4L
                spreadvalue = velev4L
            endif
        case (5)
            if (upper) then
                signal = band5U
                spreadvalue = velev5U
            else
                signal = band5L
                spreadvalue = velev5L
            endif
        case (6)
            if (upper) then
                signal = band6U
                spreadvalue = velev6U
            else
                signal = band6L
                spreadvalue = velev6L
            endif
        case (7)
            if (upper) then
                signal = band7U
                spreadvalue = velev7U
            else
                signal = band7L
                spreadvalue = velev7L
            endif
        case (8)
            if (upper) then
                signal = band8U
                spreadvalue = velev8U
            else
                signal = band8L
                spreadvalue = velev8L
            endif
        case (9)
            if (upper) then
                signal = band9U
                spreadvalue = velev9U
            else
                signal = band9L
                spreadvalue = velev9L
            endif
        case (10)
            if (upper) then
                signal = band10U
                spreadvalue = velev10U
            else
                signal = band10L
                spreadvalue = velev10L
            endif
        case (11)
            if (upper) then
                signal = band11U
                spreadvalue = velev11U
            else
                signal = band11L
                spreadvalue = velev11L
            endif
        case (12)
            if (upper) then
                signal = band12U
                spreadvalue = velev12U
            else
                signal = band12L
                spreadvalue = velev12L
            endif
        case (13)
            if (upper) then
                signal = band13U
                spreadvalue = velev13U
            else
                signal = band13L
                spreadvalue = velev13L
            endif
        case (14)
            if (upper) then
                signal = band14U
                spreadvalue = velev14U
            else
                signal = band14L
                spreadvalue = velev14L
            endif
        case (15)
            if (upper) then
                signal = band15U
                spreadvalue = velev15U
            else
                signal = band15L
                spreadvalue = velev15L
            endif
        case (16)
            if (upper) then
                signal = band16U
                spreadvalue = velev16U
            else
                signal = band16L
                spreadvalue = velev16L
            endif
        case (17)
            if (upper) then
                signal = band17U
                spreadvalue = velev17U
            else
                signal = band17L
                spreadvalue = velev17L
            endif
        case (18)
            if (upper) then
                signal = band18U
                spreadvalue = velev18U
            else
                signal = band18L
                spreadvalue = velev18L
            endif
        case (19)
            if (upper) then
                signal = band19U
                spreadvalue = velev19U
            else
                signal = band19L
                spreadvalue = velev19L
            endif
        case (20)
            if (upper) then
                signal = band20U
                spreadvalue = velev20U
            else
                signal = band20L
                spreadvalue = velev20L
            endif
        case (21)
            if (upper) then
                call MLSMessage(MLSMSG_Error, moduleName, "Band 21U does not exist.")
            else
                signal = band21L
                spreadvalue = velev21L
            endif
        case (22)
            if (upper) then
                call MLSMessage(MLSMSG_Error, moduleName, "Band 22U does not exist.")
            else
                signal = band22L
                spreadvalue = velev22L
            endif
        case (23)
            if (upper) then
                signal = band23U
                spreadvalue = velev23U
            else
                signal = band23L
                spreadvalue = velev23L
            endif
        case (24)
            if (upper) then
                signal = band24U
                spreadvalue = velev24U
            else
                signal = band24L
                spreadvalue = velev24L
            endif
        case (25)
            if (upper) then
                signal = band25U
                spreadvalue = velev25U
            else
                signal = band25L
                spreadvalue = velev25L
            endif
        case (26)
            if (upper) then
                call MLSMessage(MLSMSG_Error, moduleName, "Band 26U does not exist.")
            else
                signal = band26L
                spreadvalue = velev26L
            endif
        case (27)
            if (upper) then
                signal = band27U
                spreadvalue = velev27U
            else
                signal = band27L
                spreadvalue = velev27L
            endif
        case (28)
            if (upper) then
                signal = band28U
                spreadvalue = velev28U
            else
                signal = band28L
                spreadvalue = velev28L
            endif
        case (29)
            if (upper) then
                signal = band29U
                spreadvalue = velev29U
            else
                signal = band29L
                spreadvalue = velev29L
            endif
        case (30)
            if (upper) then
                signal = band30U
                spreadvalue = velev30U
            else
                signal = band30L
                spreadvalue = velev30L
            endif
        case (31)
            if (upper) then
                signal = band31U
                spreadvalue = velev31U
            else
                signal = band31L
                spreadvalue = velev31L
            endif
        case (32)
            if (upper) then
                call MLSMessage(MLSMSG_Error, moduleName, "Band 32U does not exist.")
            else
                signal = band32L
            endif
        case (33)
            if (upper) then
                signal = band33U
            else
                signal = band33L
            endif
        case (34)
            if (upper) then
                call MLSMessage(MLSMSG_Error, moduleName, "Band 34U does not exist.")
            else
                signal = band34L
            endif
        case default
            call MLSMessage (MLSMSG_Error, moduleName, "Band does not exist.")
        end select

        template = CreateQtyTemplate (l_elevOffset, qSignal=signal, qName=qname)

        if (upper) then
            if (band == 33) then
                vv = CreateValue4AgileVector(template, value=velev33U)
            else
                vv = CreateValue4AgileVector(template, spreadvalue=spreadvalue)
            endif
        else
            select case(band)
            case (32)
                vv = CreateValue4AgileVector(template, value=velev32L)
            case (33)
                vv = CreateValue4AgileVector(template, value=velev33L)
            case (34)
                vv = CreateValue4AgileVector(template, value=velev34L)
            case default
                vv = CreateValue4AgileVector(template, spreadvalue=spreadvalue)
            end select
        endif
    end function

    ! Create a vector value object that contains the value for limb sideband 
    ! fraction of the specified upper or lower band, and give it an optional name.
    ! The value for this quantity are MLS-specified constants, which are declared
    ! in cfm_constants.f90.
    type(VectorValue_T) function CreateMLSValue_LSF (band, upper, qname) result (vv)
        integer, intent(in) :: band
        logical, intent(in) :: upper ! if upper is false, then it's lower
        character(len=*), intent(in), optional :: qname

        type(QuantityTemplate_T) :: template
        character(len=32) :: signal
        real(r8) :: spreadvalue

        select case(band)
        case (1)
            if (upper) then
                call MLSMessage(MLSMSG_Error, moduleName, "Band 1U does not exist.")
            else
                signal = band1L
                spreadvalue = vlimbSidebandFraction1L
            endif
        case (2)
            if (upper) then
                signal = band2U
            else
                signal = band2L
            endif
        case (3)
            if (upper) then
                signal = band3U
            else
                signal = band3L
            endif
        case (4)
            if (upper) then
                signal = band4U
            else
                signal = band4L
            endif
        case (5)
            if (upper) then
                signal = band5U
            else
                signal = band5L
            endif
        case (6)
            if (upper) then
                signal = band6U
            else
                signal = band6L
            endif
        case (7)
            if (upper) then
                signal = band7U
            else
                signal = band7L
            endif
        case (8)
            if (upper) then
                signal = band8U
            else
                signal = band8L
            endif
        case (9)
            if (upper) then
                signal = band9U
            else
                signal = band9L
            endif
        case (10)
            if (upper) then
                signal = band10U
                spreadvalue = vlimbSidebandFraction10U
            else
                signal = band10L
                spreadvalue = vlimbSidebandFraction10L
            endif
        case (11)
            if (upper) then
                signal = band11U
                spreadvalue = vlimbSidebandFraction11U
            else
                signal = band11L
                spreadvalue = vlimbSidebandFraction11L
            endif
        case (12)
            if (upper) then
                signal = band12U
                spreadvalue = vlimbSidebandFraction12U
            else
                signal = band12L
                spreadvalue = vlimbSidebandFraction12L
            endif
        case (13)
            if (upper) then
                signal = band13U
                spreadvalue = vlimbSidebandFraction13U
            else
                signal = band13L
                spreadvalue = vlimbSidebandFraction13L
            endif
        case (14)
            if (upper) then
                signal = band14U
                spreadvalue = vlimbSidebandFraction14U
            else
                signal = band14L
                spreadvalue = vlimbSidebandFraction14L
            endif
        case (15)
            if (upper) then
                signal = band15U
            else
                signal = band15L
            endif
        case (16)
            if (upper) then
                signal = band16U
            else
                signal = band16L
            endif
        case (17)
            if (upper) then
                signal = band17U
            else
                signal = band17L
            endif
        case (18)
            if (upper) then
                signal = band18U
            else
                signal = band18L
            endif
        case (19)
            if (upper) then
                signal = band19U
            else
                signal = band19L
            endif
        case (20)
            if (upper) then
                signal = band20U
            else
                signal = band20L
            endif
        case (21)
            if (upper) then
                call MLSMessage(MLSMSG_Error, moduleName, "Band 21U does not exist.")
            else
                signal = band21L
                spreadvalue = vlimbSidebandFraction21L
            endif
        case (22)
            if (upper) then
                call MLSMessage(MLSMSG_Error, moduleName, "Band 22U does not exist.")
            else
                signal = band22L
                spreadvalue = vlimbSidebandFraction22L
            endif
        case (23)
            if (upper) then
                signal = band23U
            else
                signal = band23L
            endif
        case (24)
            if (upper) then
                signal = band24U
            else
                signal = band24L
            endif
        case (25)
            if (upper) then
                signal = band25U
            else
                signal = band25L
            endif
        case (26)
            if (upper) then
                call MLSMessage(MLSMSG_Error, moduleName, "Band 26U does not exist.")
            else
                signal = band26L
                spreadvalue = vlimbSidebandFraction26L
            endif
        case (27)
            if (upper) then
                signal = band27U
            else
                signal = band27L
            endif
        case (28)
            if (upper) then
                signal = band28U
            else
                signal = band28L
                spreadvalue = vlimbSidebandFraction28L
            endif
        case (29)
            if (upper) then
                signal = band29U
            else
                signal = band29L
                spreadvalue = vlimbSidebandFraction29L
            endif
        case (30)
            if (upper) then
                signal = band30U
            else
                signal = band30L
                spreadvalue = vlimbSidebandFraction30L
            endif
        case (31)
            if (upper) then
                signal = band31U
            else
                signal = band31L
                spreadvalue = vlimbSidebandFraction31L
            endif
        case (32)
            if (upper) then
                call MLSMessage(MLSMSG_Error, moduleName, "Band 32U does not exist.")
            else
                signal = band32L
            endif
        case (33)
            if (upper) then
                signal = band33U
            else
                signal = band33L
            endif
        case (34)
            if (upper) then
                call MLSMessage(MLSMSG_Error, moduleName, "Band 34U does not exist.")
            else
                signal = band34L
            endif
        case default
            call MLSMessage (MLSMSG_Error, moduleName, "Band does not exist.")
        end select

        template = CreateQtyTemplate (l_limbsidebandFraction, qSignal=signal, qName=qname)

        if (upper) then
            select case(band)
            case (2)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction2U)
            case (3)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction3U)
            case (4)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction4U)
            case (5)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction5U)
            case (6)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction6U)
            case (7)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction7U)
            case (8)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction8U)
            case (9)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction9U)
            case (15)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction15U)
            case (16)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction16U)
            case (17)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction17U)
            case (18)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction18U)
            case (19)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction19U)
            case (20)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction20U)
            case (23)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction23U)
            case (24)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction24U)
            case (25)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction25U)
            case (27)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction27U)
            case (33)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction33U)
            case default
                vv = CreateValue4AgileVector(template, spreadvalue=spreadvalue)
            end select
        else
            select case(band)
            case (2)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction2L)
            case (3)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction3L)
            case (4)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction4L)
            case (5)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction5L)
            case (6)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction6L)
            case (7)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction7L)
            case (8)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction8L)
            case (9)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction9L)
            case (15)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction15L)
            case (16)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction16L)
            case (17)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction17L)
            case (18)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction18L)
            case (19)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction19L)
            case (20)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction20L)
            case (23)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction23L)
            case (24)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction24L)
            case (25)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction25L)
            case (27)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction27L)
            case (32)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction32L)
            case (33)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction33L)
            case (34)
                vv = CreateValue4AgileVector(template, value=vlimbSidebandFraction34L)
            case default
                vv = CreateValue4AgileVector(template, spreadvalue=spreadvalue)
            end select
        endif
    end function

    ! Add constant quantities to stateVectorExtra.
    ! The constant quantities are: earth reflectivity,
    ! space radiance, elevation offset for all bands, 
    ! and limb sideband fraction for all bands.
    subroutine GetConstantQuantities (stateVectorExtra)
        type (Vector_T), intent(inout) :: stateVectorExtra

        type(VectorValue_T) :: qty
        integer :: band

        qty = CreateMLSValue_EarthReflectivity()
        call AddValue2Vector(stateVectorExtra, qty)

        qty = CreateMLSValue_SpaceRadiance()
        call AddValue2Vector(stateVectorExtra, qty)

        do band=1,34
            if (band .ne. 1 .and. band .ne. 21 .and. band .ne. 22 &
            .and. band .ne. 26 .and. band .ne. 32) then
                qty = CreateMLSValue_ElevationOffset (band, .true.)
                call AddValue2Vector(stateVectorExtra, qty)

                qty = CreateMLSValue_lsf (band, .true.)
                call AddValue2Vector(stateVectorExtra, qty)
            end if

            qty = CreateMLSValue_ElevationOffset (band, .false.)
            call AddValue2Vector(stateVectorExtra, qty)

            qty = CreateMLSValue_lsf (band, .false.)
            call AddValue2Vector(stateVectorExtra, qty)
        enddo
    end subroutine

    ! Create stateVectorExtra with the given filedatabase, chunk, and quantity database.
    ! StateVectorExtra will include the following quantities: O2, earth reflectivity,
    ! LOS velocity of GHz module, spacecraft geocentric altitude, space radiance,
    ! elevation offset, limb sideband fraction.
    ! See CFM_Constants_m for a list of bands.
    subroutine CreateStateVectorExtra (filedatabase, chunk, qtyTemplates, stateVectorExtra)
        use CFM_VGrid_m, only: CreateVGrid
        use VGridsDatabase, only: DestroyVGridContents
        use CFM_HGrid_m, only: CreateRegularHGrid
        use HGridsDatabase, only: DestroyHGridContents
        use QuantityTemplates, only: AddQuantityTemplateToDatabase
        use CFM_VectorTemplate_m, only: CreateVectorTemplate
        use CFM_Vector_m, only: CreateVector
        use Allocate_Deallocate, only: allocate_test, deallocate_test
        use CFM_LSF_M, only: CreateLimbSidebandFractions, FillLimbSidebandFractions
        use CFM_EO_M, only: CreateElevationOffsets, FillElevationOffsets

        type (MLSFile_T), dimension(:), pointer :: filedatabase
        type (MLSChunk_T), intent(in) :: chunk
        type (QuantityTemplate_T), dimension(:), pointer :: qtyTemplates
        type (Vector_T), intent(out) :: stateVectorExtra

        type(QuantityTemplate_T) :: O2, earthRefl, losVelGHz, scGeocAlt, &
                                    spaceRadiance
        type(VectorTemplate_T) :: stateTemplateExtra
        type(VGrid_T) :: vGridStandard
        type(HGrid_T) :: hGridStandard
        integer :: O2_index, earthRefl_index, losVelGHz_index, scGeocAlt_index, &
                   spaceRadiance_index
        type(VectorValue_T), pointer :: quantity
        integer :: i
        integer, dimension(:), pointer :: selected

        ! Executables

        nullify(qtyTemplates)
        ! Create O2 for use in the forward model
        vGridStandard = CreateVGrid (l_zeta, phyq_pressure, l_logarithmic, &
                                     start=1000.0d0, formula="37:6")

        ! Have insetoverlaps, and not single
        hGridStandard = CreateRegularHGrid("GHz", 0.0_r8, 1.5_r8, .true., &
                                           filedatabase, chunk)

        O2 = CreateQtyTemplate(l_vmr, filedatabase=filedatabase, chunk=chunk, &
                               avgrid=vGridStandard, ahgrid=hGridStandard, &
                               qMolecule=l_o2)
        earthRefl = CreateQtyTemplate(l_earthRefl)
        losVelGHz = CreateQtyTemplate(l_losVel, filedatabase=filedatabase, chunk=chunk, &
                                      qInstModule="GHz")
        scGeocAlt = CreateQtyTemplate(l_scgeocalt, filedatabase=filedatabase, chunk=chunk, &
                                      qInstModule="sc")
        spaceRadiance = CreateQtyTemplate(l_spaceradiance)

        o2_index = AddQuantityTemplateToDatabase(qtyTemplates, o2)
        earthRefl_index = AddQuantityTemplateToDatabase(qtyTemplates, earthRefl)
        losVelGHz_index = AddQuantityTemplateToDatabase(qtyTemplates, losVelGHz)
        scGeocAlt_index = AddQuantityTemplateToDatabase(qtyTemplates, scGeocAlt)
        spaceRadiance_index = AddQuantityTemplateToDatabase(qtyTemplates, spaceRadiance)
        ! Add sideband fraction
        call CreateLimbSidebandFractions (chunk, filedatabase, qtyTemplates)
        call CreateElevationOffsets (chunk, filedatabase, qtyTemplates)

        ! Don't need the grids anymore
        call DestroyVGridContents(vGridStandard)
        call DestroyHGridContents(hGridStandard)

        ! Everything in the the quantity template database is selected
        nullify(selected)
        call allocate_test (selected, size(qtyTemplates), "selected", moduleName)
        do i = 1, size(selected)
            selected(i) = i
        end do

        stateTemplateExtra = CreateVectorTemplate(qtyTemplates, selected)

        call deallocate_Test(selected, "selected", moduleName)

        stateVectorExtra = CreateVector(stateTemplateExtra, qtyTemplates, name='stateExtra')

        quantity => GetVectorQtyByTemplateIndex (stateVectorExtra, o2_index)
        call FillVectorQtyFromProfile (quantity, .false., o2_heights, &
                                       o2_values, phyq_vmr)
        !print *, "O2 value"
        !call dump(quantity, details=3)
        quantity => GetVectorQtyByTemplateIndex (stateVectorExtra, earthRefl_index)
        call ExplicitFillVectorQuantity (quantity, earthRefl_values)
        !print *, "Earth Reflectivity value"
        !call dump(quantity, details=3)
        quantity => GetVectorQtyByTemplateIndex (stateVectorExtra, spaceRadiance_index)
        call ExplicitFillVectorQuantity (quantity, spaceRad_values)
        !print *, "Space radiance"
        !call dump(quantity, details=3)
        quantity => GetVectorQtyByTemplateIndex (stateVectorExtra, losVelGHz_index)
        call FillVectorQuantityFromL1B (quantity, chunk, filedatabase, .false.)
        !print *, "LOS Velocity GHz"
        !call dump(quantity, details=3)
        quantity => GetVectorQtyByTemplateIndex (stateVectorExtra, scGeocAlt_index)
        call FillVectorQuantityFromL1B (quantity, chunk, filedatabase, .false.)
        !print *, "scGeocAlt value"
        !call dump(quantity, details=3)
        call FillLimbSidebandFractions(stateVectorExtra)
        call FillElevationOffsets(stateVectorExtra)
    end subroutine

    ! This subroutine convert a time range to a pair of L1 mafs from 
    ! a L1BOA file.
    subroutine timeRange2MafRange (startTime, endTime, leapsecfile, filedatabase, &
    startL1Maf, endL1Maf)
        use INIT_TABLES_MODULE, only: l_ghz, phyq_mafs, l_orbital
        use SDPToolkit, only: mls_utctotai

        ! The start time of the data to be read in the format
        ! yyyy-doyThh:mm:ss.zzzz
        character(len=CCSDSlen), intent(in) :: startTime
        ! The stop time of the data to be read in the same format as startTime.
        character(len=CCSDSlen), intent(in) :: endTime
        ! The leap second file to correct time calculation.
        character(len=*), intent(in) :: leapsecfile
        ! The file database that contains an already opened L1BOA.
        type (MLSFile_T), dimension(:), pointer :: filedatabase
        ! The output L1 maf that correspond to startTime
        integer, intent(out) :: startL1Maf
        ! The output L1 maf that correspond to stopTime
        integer, intent(out) :: endL1Maf

        type (TAI93_Range_T) :: processingRange
        integer :: error
        type (MLSChunk_T) :: chunk

        ! Get the start time and end time encoded
        error = mls_utctotai (trim(leapsecFile), startTime, processingrange%starttime)
        if (error /= 0) &
            call MLSMessage (MLSMSG_Error, moduleName, &
            "Could not convert UTC Start time to TAI")

        error = mls_utctotai (trim(leapsecFile), endtime, processingrange%endtime)
        if (error /= 0) &
            call MLSMessage (MLSMSG_Error, moduleName, &
            "Could not convert UTC End time to TAI")

        ! Initialize ChunkDivideConfig
        chunkDivideConfig%method = l_orbital
        chunkDivideConfig%maxLength = 3500
        chunkDivideConfig%maxLengthFamily = phyq_mafs
        chunkDivideConfig%skipL1BCheck = .true.
        chunkDivideConfig%homeModule = l_ghz
        chunkDivideConfig%criticalModules = l_ghz
        chunkDivideConfig%homeGeodAngle = 0.0_r8
        chunkDivideConfig%maxGap = 2.0_r8
        chunkDivideConfig%maxGapFamily = phyq_mafs
        chunkDivideConfig%maxOrbY = 50000.0_r8 ! Unit is meter
        chunkDivideConfig%scanLowerLimit = (/-20000.0_r8, 10000.0_r8 /) ! unit is meter
        chunkDivideConfig%scanUpperLimit = (/ 40000.0_r8, 200000.0_r8 /) ! unit is meter
        chunkDivideConfig%lowerOverlap = 0.0
        chunkDivideConfig%upperOverlap = 0.0
        chunkDivideConfig%lowerOverlapFamily = phyq_mafs
        chunkDivideConfig%upperOverlapFamily = phyq_mafs
        chunkDivideConfig%noChunks = 1
        ! Create a fake chunk out of the start time, end time, and L1BOA
        chunk = GetChunkFromTimeRange(processingRange, filedatabase, chunkDivideConfig)

        startL1Maf = chunk%firstMafIndex
        endL1Maf = chunk%lastMafIndex
    end subroutine

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

! $Log$
! Revision 1.34  2011/12/15 18:27:44  honghanh
! Documentation and code clean up, including removing unused and broken
! subroutines.
!
! Revision 1.33  2011/12/15 16:53:24  honghanh
! Correct the name of CreateMLSValue_EarthReflectivity
!
! Revision 1.32  2011/12/14 22:54:18  honghanh
! Add timeRange2MafRange method in CFM.
!
! Revision 1.30  2011/11/08 16:13:57  honghanh
! Add 'qname' parameter to CreateMLSValue_O2 subroutine.
!
! Revision 1.29  2011/10/31 20:04:59  honghanh
! Change CreateMLSValue_FromL1B to CreateMLSValue_FromL1BOA.
!
! Revision 1.28  2011/10/20 00:23:07  honghanh
! Add elevation offset creation subroutine to public API
!
! Revision 1.27  2011/10/19 11:35:28  honghanh
! Add extra APIs for creating MLS vector values.
!
! Revision 1.26  2011/10/18 17:13:35  honghanh
! Use earth reflectivity and space radiance values from CFM_Constants_m.
!
! Revision 1.25  2011/10/18 17:04:01  honghanh
! Move InitQuantityTemplates to setup subroutines.
! (Previously in create statevectorExtra).
!
! Revision 1.24  2011/10/17 20:41:02  honghanh
! Add a more concise CFM_MLSSetup/CFM_MLSCleanup subroutines
! and extract literal constants from CFM_MLSSetup_m to put in
! CFM_Constants_m
!
! Revision 1.23  2011/03/23 20:09:46  honghanh
! nullify selected array before calling allocate_test
!
! Revision 1.22  2010/11/03 20:17:01  honghanh
! Add name as an optional argument to CreateVector.
!
! Revision 1.21  2010/09/17 16:47:45  honghanh
! Fix memory leak bug by calling CFM_ResetEmpiricalGeometry
!
! Revision 1.18  2010/07/08 21:39:16  honghanh
! Add ApplyBaseline to cfm_fill_m
!
! Revision 1.17  2010/06/29 17:02:47  honghanh
! Change the identifier 'fakeChunk' to 'chunk' because
! since it is created with ChunkDivide, it's as real as a chunk
! can get.
!
! Revision 1.16  2010/06/29 16:40:23  honghanh
! Remove all function/subroutine and user type forwarding from
! all CFM modules except for from cfm.f90
!
! Revision 1.15  2010/06/29 15:53:45  honghanh
! Add copyright comments and support for CVS log in the file
!
