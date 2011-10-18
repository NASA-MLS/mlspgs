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
                                  l_scgeocalt, l_spaceradiance, l_o2
    use ConstructQuantityTemplates, only: InitQuantityTemplates
    use ChunkDivide_m, only: ChunkDivideConfig_T, &
                             GetChunkFromTimeRange => CFM_ChunkDivide
    use Chunks_m, only: MLSChunk_T ! To also be referenced outside
    use VGridsDatabase, only: VGrid_T
    use HGridsDatabase, only: HGrid_T
    use VectorsModule, only: GetVectorQtyByTemplateIndex, & ! being used in many functions
                             VectorTemplate_T, Dump, Vector_T, VectorValue_T
    use CFM_Vector_m, only: CreateValue4AgileVector
    use CFM_Constants_m
    use CFM_QuantityTemplate_m, only: CreateQtyTemplate
    use CFM_Fill_M, only: FillVectorQtyFromProfile, ExplicitFillVectorQuantity, &
        FillVectorQuantityFromL1B

    implicit none

    private
    public :: CFM_MLSSetup, CFM_MLSCleanup

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

    subroutine CFM_MLSSetup_Compact (signalFileName, configFileName, ForwardModelConfigDatabase)
        use CFM_Tree_Walker_m, only : Walk_Tree
        use EmpiricalGeometry, only: CFM_InitEmpiricalGeometry
        use string_table, only: AddInUnit
        use io_stuff, only: get_lun
        use INIT_TABLES_MODULE, only: INIT_TABLES
        use PARSER, only: CONFIGURATION
        use TREE_CHECKER, only: CHECK_TREE

        ! The name of the signal database file
        character(len=*), intent(in) :: signalFileName
        ! The name of the forward model configuration file
        character(len=*), intent(in) :: configFileName
        ! Output: to store all ForwardModelConfig_T objects declared in
        ! the config file.
        ! This argument will be nullified inside this subroutine.
        type (ForwardModelConfig_T), pointer :: ForwardModelConfigDatabase(:)

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

        call configuration(Root)
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
        use PARSER, only: CONFIGURATION
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

        call configuration(Root)
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

    ! Clean up allocated memory, close opened files
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

    type(VectorValue_T) function CreateMLSValue_O2 (avgrid, ahgrid) result(o2)
        type(VGrid_T), intent(in) :: avgrid
        type(HGrid_T), intent(in) :: ahgrid

        type(QuantityTemplate_T) :: o2template

        o2template = CreateQtyTemplate(l_vmr, avgrid=avgrid, ahgrid=ahgrid, qMolecule=l_o2)

        o2 = CreateValue4AgileVector(o2template)
        call FillVectorQtyFromProfile (o2, .false., o2_heights, o2_values, phyq_vmr)
    end function

    type(VectorValue_T) function CreateMLSValue_EarthReflectiviy () result(earthRefl)
        type(QuantityTemplate_T) :: ertemplate

        ertemplate = CreateQtyTemplate(l_earthRefl)
        earthRefl = CreateValue4AgileVector(ertemplate, value=earthRefl_values)
    end function

    type(VectorValue_T) function CreateMLSValue_SpaceRadiance () result(spaceRad)
        type(QuantityTemplate_T) :: srtemplate

        srtemplate = CreateQtyTemplate(l_spaceradiance)
        spaceRad = CreateValue4AgileVector(srtemplate, value=spaceRad_values)
    end function

    type(VectorValue_T) function CreateMLSValue_FromL1B (qType, instrumentModule, filedatabase, &
    firstL1Maf, lastL1Maf) result (vv)
        integer, intent(in) :: qType
        character(len=*), intent(in) :: instrumentModule
        type (MLSFile_T), dimension(:), pointer :: filedatabase
        integer, intent(in) :: firstL1Maf
        integer, intent(in) :: lastL1Maf

        type(MLSChunk_T) :: chunk
        type(QuantityTemplate_T) :: template

        chunk%firstMafIndex = firstL1Maf
        chunk%lastMafIndex = lastL1Maf

        template = CreateQtyTemplate(qType, filedatabase=filedatabase, chunk=chunk, &
                                     qInstModule=instrumentModule)
        vv = CreateValue4AgileVector(template)
        call FillVectorQuantityFromL1B (vv, chunk, filedatabase, .false.)
    end function

    ! Create stateVectorExtra with the given filedatabase, chunk, and quantity database.
    ! StateVectorExtra will include the following quantities: O2, earth reflectivity,
    ! LOS velocity of GHz module, spacecraft geocentric altitude, space radiance,
    ! elevation offset, limb sideband fraction, reference GPH, phitan of GHz module
    ! for every band. See CFM documentation for a list of bands.
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
        call ExplicitFillVectorQuantity (quantity, (/0.05_r8/))
        !print *, "Earth Reflectivity value"
        !call dump(quantity, details=3)
        quantity => GetVectorQtyByTemplateIndex (stateVectorExtra, spaceRadiance_index)
        call ExplicitFillVectorQuantity (quantity, (/2.735_r8/))
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
