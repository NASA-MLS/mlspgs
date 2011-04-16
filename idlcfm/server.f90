program server
    use QuantityPVM
    use CFM
    use IDLCFM2_m
    use MLSCommon, only: r8
    use DECLARATION_TABLE, only: ALLOCATE_DECL, DEALLOCATE_DECL
    use TREE, only: ALLOCATE_TREE, DEALLOCATE_TREE
    use SYMBOL_TABLE, only: DESTROY_SYMBOL_TABLE
    use MLSMessageModule, only: MLSMessageConfig, PVMERRORMESSAGE
    use PVMIDL, only: PVMIDLPACK, PVMIDLUNPACK
    use PVM, only: PVMFRECV, PVMFBUFINFO
    use H5LIB, ONLY: h5open_f, h5close_f
    use EmpiricalGeometry, only: CFM_InitEmpiricalGeometry

    implicit none

!---------------------------- RCS Ident Info ------------------------------
    character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
    character (len=*), parameter :: IdParm = &
       "$Id$"
    character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------
    integer, parameter :: SIG_SETUP = 0
    integer, parameter :: SIG_CLEANUP = 1
    integer, parameter :: SIG_FWDMDL = 2

    integer :: tid, info, bufid, bytes, msgtag, cid
    integer :: reqcode
    type (ForwardModelConfig_T), pointer :: ForwardModelConfigDatabase(:)
    logical :: locked = .false.

    MLSMessageConfig%logFileUnit = -1
    MLSMessageConfig%useToolkit = .false.

    call PVMFMyTid(tid)

    if (tid <= 0) &
        call MLSMessage (MLSMSG_Error, moduleName, "Can't contact PVM daemon")

    print *, "Server start with tid ", tid

    do while (.true.)
        reqcode = -1 ! init to some number that is not one of the signals
        call PVMFrecv ( -1, QtyMsgTag, bufid )
        call PVMIDLUnpack (reqcode, info)
        if ( info /= 0 ) then
            call PVMErrorMessage ( info, "unpacking reqcode." )
            cycle
        endif
        call pvmfbufinfo(bufid, bytes, msgtag, cid, info)
        if (info /= 0) then
            call PVMErrorMessage (info, "get bufinfo")
            cycle
        endif
        print *, "cid is ", cid

        select case (reqcode)
        case (SIG_SETUP)
            call ICFM_Setup (info)
        case (SIG_CLEANUP)
            call ICFM_Cleanup
        case (SIG_FWDMDL)
            call ICFM_ForwardModel (cid, info)
        end select
    enddo

    call PVMFExit(tid, info)

    contains

    subroutine ICFM_Setup (info)
        use LEXER_CORE, only: INIT_LEXER
        use ConstructQuantityTemplates, only: InitQuantityTemplates
        use string_table, only: AddInUnit
        use io_stuff, only: get_lun
        use PARSER, only: CONFIGURATION
        use TREE_CHECKER, only: CHECK_TREE
        use CFM_Tree_Walker_m, only : Walk_Tree
        use INIT_TABLES_MODULE, only: INIT_TABLES

        integer, intent(out) :: info

        character(len=128) :: signalfile, configFile, spectroscopy, antennaPatterns, filterShapes
        character(len=128) :: dacsFilterShapes, pointingGrids, pfa, l2pc

        integer :: signalIn, configIn, root, error, first_section
        integer :: numFiles, i

        integer, parameter :: empiricalGeometry_noIterations = 10
        real(r8), dimension(21), parameter :: empircalGeometry_terms = (/ &
        -1.06863, 43.0943, -16.2062, 8.12730, -4.58416, 2.75786, -1.72880, &
        1.11523, -0.733464, 0.489792, -0.331852, 0.227522, -0.156428, &
        0.108031, -0.0757825, 0.0536980, -0.0375161, 0.0260555, &
        -0.0188811, 0.0138453, -0.00959350 /)

        print *, "call ICFM_Setup"

        if (locked) then
            call MLSMessage (MLSMSG_Warning, moduleName, "a session is already started")
            return
        endif

        call CFM_InitEmpiricalGeometry (empiricalGeometry_noIterations, &
                                        empircalGeometry_terms)

        call PVMIDLUnpack (signalfile, info)
        if (info /= 0) then
            call PVMErrorMessage ( info, "unpacking signalfile." )
            return
        endif
        print *, "signalfile ", signalfile

        call PVMIDLUnpack (configfile, info)
        if (info /= 0) then
            call PVMErrorMessage ( info, "unpacking configfilename." )
            return
        endif
        print *, "configfile ", configfile

        ! Read signal database
        call get_lun(signalIn)
        open (unit=signalIn, file=trim(signalfile), status='OLD', iostat=error)
        if (error /= 0) call MLSMessage (MLSMSG_Error, moduleName, &
            'Error opening ' // trim(signalfile))

        call get_lun(configIn)
        open(unit=configIn, file=trim(configFile), status='OLD', iostat=error)
        if (error /= 0) &
            call MLSMessage (MLSMSG_Error, moduleName, &
            'Error opening ' // trim(configFile))

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
        if (error /= 0) call MLSMessage (MLSMSG_Error, moduleName, "Error in tree")

        nullify(ForwardModelConfigDatabase)
        call Walk_Tree ( Root, First_Section, &
            ForwardModelConfigDatabase=ForwardModelConfigDatabase )

        close (signalIn, iostat=error)
        if (error /= 0) call MLSMessage (MLSMSG_Error, moduleName, &
            "Error closing " // trim(signalfile))
        close (configIn, iostat=error)
        if (error /= 0) &
            call MLSMessage (MLSMSG_Error, moduleName, &
            'Error closing ' // trim(configFile))

        ! init property table
        call InitQuantityTemplates

        ! We have to call this before opening any HDF5 file
        call h5open_f(error)
        if (error /= 0) &
            call MLSMessage (MLSMSG_Error, moduleName, "Error in initialize hdf5 library")

        call PVMIDLUnpack (spectroscopy, info)
        if (info /= 0) then
            call PVMErrorMessage ( info, "unpacking spectroscopy." )
            return
        endif
        print *, "spectroscopy ", spectroscopy
        call Read_Spectroscopy (spectroscopy, 'HDF5')

        call PVMIDLUnpack (antennaPatterns, info)
        if (info /= 0) then
            call PVMErrorMessage ( info, "unpacking antennaPatterns." )
            return
        endif
        print *, "antennaPatterns ", antennaPatterns
        call ReadAntennaPatterns (antennaPatterns)

        call PVMIDLUnpack (filterShapes, info)
        if (info /= 0) then
            call PVMErrorMessage ( info, "unpacking filterShapes." )
            return
        endif
        print *, "filterShapes ", filterShapes
        call ReadFilterShapes(filterShapes)

        call PVMIDLUnpack (dacsFilterShapes, info)
        if (info /= 0) then
            call PVMErrorMessage ( info, "unpacking dacsFilterShapes." )
            return
        endif
        print *, "dacsFilterShapes ", dacsFilterShapes
        call ReadDACSFilterShapes (DACSFilterShapes)

        call PVMIDLUnpack (pointinggrids, info)
        if (info /= 0) then
            call PVMErrorMessage ( info, "unpacking pointinggrids." )
            return
        endif
        print *, "pointinggrids ", pointinggrids
        call ReadPointingGrids (pointingGrids)

        call PVMIDLUnpack (numFiles, info)
        if (info /= 0) then
            call PVMErrorMessage ( info, "unpacking numFiles." )
            return
        endif
        print *, "numFiles ", numFiles

        do i = 1, numFiles
            call PVMIDLUnpack(pfa, info)
            if (info /= 0) then
                call PVMErrorMessage ( info, "unpacking pfa." )
                return
            endif
            print *, "pfa ", pfa
            call ReadPFAFile (pfa)
        enddo

        call PVMIDLUnpack (numFiles, info)
        if (info /= 0) then
            call PVMErrorMessage ( info, "unpacking numFiles." )
            return
        endif
        print *, "numFiles ", numFiles

        do i = 1, numFiles
            call PVMIDLUnpack(l2pc, info)
            if (info /= 0) then
                call PVMErrorMessage ( info, "unpacking l2pc." )
                return
            endif
            print *, "l2pc ", l2pc
            call ReadHDF5L2PC (l2pc)
        enddo

        locked = .true.

    end subroutine

    subroutine ICFM_Cleanup
        use STRING_TABLE, only: DESTROY_CHAR_TABLE, DESTROY_HASH_TABLE, DESTROY_STRING_TABLE
        use ForwardModelConfig, only: DestroyFWMConfigDatabase
        use MLSSignals_m, only: modules, bands, radiometers, spectrometerTypes, &
                              signals, DestroyBandDatabase, DestroySignalDatabase, &
                              DestroySpectrometerTypeDatabase, DestroyModuleDatabase, &
                              DestroyRadiometerDatabase
        use EmpiricalGeometry, only: CFM_ResetEmpiricalGeometry

        integer :: error

        print *, "call ICFM_Cleanup"

        if (.not. locked) then
            call MLSMessage (MLSMSG_Warning, moduleName, "there is no session to clean up")
            return
        endif

        call Destroy_DACS_Filter_Database
        call Destroy_Filter_Shapes_Database
        call Destroy_Ant_Patterns_Database
        call Destroy_SpectCat_Database
        call Destroy_Line_Database
        call Destroy_Pointing_Grid_Database
        call DestroyL2PCDatabase
        call Destroy_PFADataBase

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

        locked = .false.

    end subroutine

    subroutine ICFM_ForwardModel (tid, info)
        use ForwardModelIntermediate, only: FORWARDMODELSTATUS_T
        use ForwardModelWrappers, only: ForwardModel

        integer, intent(in) :: tid
        integer, intent(out) :: info

        type(QuantityTemplate_T), dimension(:), pointer :: qtydb
        type(Vector_T) :: state, stateExtra, measurement
        integer :: mafNo, i
        type(forwardModelStatus_t) :: FMSTAT ! Reverse comm. stuff

        if (.not. locked) then
            call MLSMessage (MLSMSG_Warning, moduleName, "must run set up first")
            return
        endif

        nullify(qtydb)

        call ICFMReceiveVector (state, qtydb)
        call ICFMReceiveVector (stateExtra, qtydb)
        call ICFMReceiveVector (measurement, qtydb)
        !call dump (state, details=3)
        !call dump(stateextra, details=3)
        call PVMIDLUnpack(mafNo, info)
        if (info /= 0) then
            call PVMErrorMessage ( info, "unpacking mafNo." )
        endif

        if (info == 0) then
!           print *, "mafNo ", mafNo ! mafNo is 0-based

            fmStat%newScanHydros = .true.
            fmStat%maf = mafNo + 1 ! in fortran it 1-based

            do i=1, size(ForwardModelConfigDatabase)
                call ForwardModel (ForwardModelConfigDatabase(i), state, stateExtra, &
                measurement, fmStat)
            enddo

            call dump(measurement)

            call ICFMSendVector(measurement, tid, info)
            if (info /= 0) then
                call PVMErrorMessage ( info, "send measurement" )
            endif
        end if

        call DestroyVectorTemplateInfo(stateExtra%template)
        call DestroyVectorInfo(stateExtra)
        call DestroyVectorTemplateInfo(state%template)
        call DestroyVectorInfo(state)
        call DestroyVectorTemplateInfo(measurement%template)
        call DestroyVectorInfo(measurement)
        call DestroyQuantityTemplateDatabase (qtydb)

    end subroutine

    subroutine senddummyvector (tid)
        use input
        use SDPToolkit, only: mls_utctotai
        use MLSFiles, only: AddFileToDatabase, InitializeMLSFile, mls_openFile
        use Hdf, only: DFACC_RDONLY
        use MLSCommon, only: MLSFile_T, TAI93_Range_T, r8
        use ChunkDivide_m, only: ChunkDivideConfig_T, &
                            GetChunkFromTimeRange => CFM_ChunkDivide
        use Chunks_m, only: MLSChunk_T
        use INIT_TABLES_MODULE, only: l_ghz, &
                                    l_none, phyq_mafs, phyq_time, phyq_angle, &
                                    l_orbital, l_fixed, l_both, l_either

        integer, intent(in) :: tid

        type(QuantityTemplate_T) :: temperature, GPH, H2O, O3, ptanGHz, &
                                    geodAltitude, orbincl, geocAlt, refGPH, phitanGHz
        type(QuantityTemplate_T), dimension(:), pointer :: qtyTemplates
        type (MLSChunk_T) :: chunk
        type (MLSFile_T), dimension(:), pointer :: filedatabase
        type (MLSFile_T), target :: l1bfile
        type (TAI93_Range_T) :: processingRange
        integer :: error
        type(ChunkDivideConfig_T) :: chunkDivideConfig
        integer :: CCSDSLen = 27
        character(len=3) :: GHz = "GHz"
        character(len=2) :: sc = "sc"
        integer :: temperature_index, h2o_index
        integer :: o3_index, ptanGHz_index, phitanGHz_index
        integer :: geodAlt_index, orbincl_index, gph_index
        integer :: geocAlt_index, refGPH_index
        type(VGrid_T) :: vGridStandard, vGridRefGPH
        type(HGrid_T) :: hGridStandard
        integer :: stateSelected(4)
        type(Vector_T) :: state
        type(VectorTemplate_T) :: stateTemplate
        type(VectorValue_T), pointer :: quantity, h2o_vv, orbincl_vv, geocAlt_vv, &
                                        ptanG_vv, temperature_vv, refGPH_vv, &
                                        phitan_vv

        nullify(filedatabase)
        error = InitializeMLSFile(l1bfile, content='l1boa', &
        name=trim(l1boa), shortName='L1BOA', type=l_hdf, access=DFACC_RDONLY)
        call mls_openFile(L1BFile, error)
        error = AddFileToDatabase(filedatabase, l1bfile)
        error = mls_utctotai (trim(leapsecFile), startTime, processingrange%starttime)
        error = mls_utctotai (trim(leapsecFile), endtime, processingrange%endtime)
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

        vGridStandard = CreateVGrid (l_zeta, phyq_pressure, l_logarithmic, &
                                     start=1000.0d0, formula="37:6")

        vGridRefGPH = CreateVGrid (l_zeta, phyq_pressure, l_explicit, &
                                   values=(/100.0_r8/))

        ! Have insetoverlaps, and not single
        hGridStandard = CreateRegularHGrid(GHz, 0.0_r8, 1.5_r8, .true., &
                                           filedatabase, chunk)

        temperature = CreateQtyTemplate(l_temperature, filedatabase=filedatabase, &
                                        chunk=chunk, qName='temperature', &
                                        avgrid=vGridStandard, ahgrid=hGridStandard)
        GPH = CreateQtyTemplate(l_gph, filedatabase=filedatabase, chunk=chunk, &
                                avgrid=vGridStandard, ahgrid=hGridStandard, qName='GPH')
        O3 = CreateQtyTemplate(l_vmr, filedatabase=filedatabase, chunk=chunk, &
                               avgrid=vGridStandard, ahgrid=hGridStandard, qMolecule=l_o3, &
                               qName='O3')
        H2O = CreateQtyTemplate(l_vmr, filedatabase=filedatabase, chunk=chunk, &
                                avgrid=vGridStandard, ahgrid=hGridStandard, qMolecule=l_h2o, &
                                qLogBasis=.true., qMinValue=0.1E-6_r8, qName='H2O')
        ptanGHz = CreateQtyTemplate(l_ptan, filedatabase=filedatabase, &
                                    chunk=chunk, qInstModule=GHz, qName='ptanGHz')

        nullify(qtyTemplates)
        temperature_index = AddQuantityTemplateToDatabase(qtyTemplates, temperature)
        gph_index = AddQuantityTemplateToDatabase(qtyTemplates, GPH)
        o3_index = AddQuantityTemplateToDatabase(qtyTemplates, O3)
        h2o_index = AddQuantityTemplateToDatabase(qtyTemplates, H2O)
        ptanGHz_index = AddQuantityTemplateToDatabase(qtyTemplates, ptanGHz)

        stateSelected = (/temperature_index,o3_index,h2o_index, ptanGHz_index/)
        stateTemplate = CreateVectorTemplate(qtyTemplates, stateSelected)

        state = CreateVector(stateTemplate, qtyTemplates, name='state')
        temperature_vv => GetVectorQtyByTemplateIndex(state, temperature_index)
        call ExplicitFillVectorQuantity(temperature_vv, TemperatureInput)

        h2o_vv => GetVectorQtyByTemplateIndex(state, h2o_index)
        call ExplicitFillVectorQuantity(h2o_vv, H2OInput)

        quantity => GetVectorQtyByTemplateIndex(state, o3_index)
        call ExplicitFillVectorQuantity(quantity, O3Input)

        call ICFMSendVector(state, tid, error)
        print *, "error ", error

    end subroutine

end program

! $Log$
! Revision 1.1  2011/03/15 15:23:51  honghanh
! Initial imports
!
