! Copyright 2010, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or proVersionIding access to foreign persons.

! This program is meant to serve as an example, and proof that the CFM library
! is working. Consequently, the design of this program does not follow good
! software design principles. This program should not be used as a part in
! any programs or software suite meant for long-term use.
program MLS_CFM_Main

       use CFM
       use input
       use machine, only: getarg


       ! To convert MAF to Profile, Pranjit Saha 
       use global_settings, only: L1MAFToL2Profile, L2ProfileToL1MAF
       use MLSHDF5, only: MLS_H5OPEN, MLS_H5CLOSE
       use HDF5, only:  h5open_f, h5close_f
       use Allocate_Deallocate, only: Allocate_Test, DeAllocate_Test
       
       

       implicit none


!---------------------------- RCS Ident Info ------------------------------
       character (len=*), parameter :: ModuleName= &
           "$RCSfile$"
       character (len=*), parameter :: IdParm = &
           "$Id$"
       character (len=len(idParm)) :: Id = idParm
       character (len=*) , parameter :: VersionId = 'MLSCFM1.57_001' ! version  ID for this executable
!---------------------------------------------------------------------------

       integer :: i, j ,ierr
       ! To write values in file, read MLS L1BOA and L2GP files, Pranjit Saha
       CHARACTER(LEN=20), PARAMETER :: FMT1 = "(F16.5)"
       CHARACTER(LEN=20), PARAMETER :: FMT2 = "(E18.8)"
       CHARACTER(LEN=256) :: string_buffer
       
       
       integer :: FileIndex2
       type (MLSFile_T), dimension(:), pointer :: FileDatabase2 => null()
       type (MLSFile_T) :: L1BFile2, L2GPFile2
       integer :: startL1Maf_new
       integer :: endL1Maf_new
       ! For reading ptan Values from L2AUX-DGM file, Pranjit Saha
       real(r8), dimension(125) :: ptanValuesRead
 
       
       type(ForwardModelConfig_T), dimension(:), pointer :: forwardModelConfigDatabase
       type(MLSFile_T), dimension(:), pointer :: filedatabase
       type(MLSFile_T) :: l1bfile
       type(VGrid_T) :: vGridStandard55, vGridStandard37, vGridTESCO, vGridRefGPH
       type(VGrid_T) :: vGridExtinction
       type(HGrid_T) :: hGridStandard
       type(FGrid_T) :: fGridExtinctionConstant
       type(QuantityTemplate_T) :: qtemp, qCO, qso2, qhno3, qo3, qextinctionv2r3
       type(QuantityTemplate_T) :: qphitanGhz, qrefGPH, qgph, qh2o, newer_qptanGHz
       type(QuantityTemplate_T) :: qband9
       type(QuantityTemplate_T) :: qptanGHz
       type(QuantityTemplate_T) :: qbaseline9
       type(Vector_T) :: state, stateExtra
       type(Vector_T) :: radiance, diffVector
       type(Vector_T) :: observed, obsPrecision
       character(len=3) :: GHz = "GHz"
       character(len=2) :: sc = "sc"
       type(VectorValue_T) :: temperature, co, o2, so2, hno3, o3, extinctionv2r3, ptanGHz
       type(VectorValue_T) :: phitanGhz, refGPH, gph, h2o, newer_ptanGHz 
       type(VectorValue_T) :: limbSidebandFraction9L
       type(VectorValue_T) :: limbSidebandFraction9U, elev9L, elev9U, earthReflectivity
       type(VectorValue_T) :: orbitInclination, spaceRadiance, scGeocAlt, tngtGeocAltGHz
       type(VectorValue_T) :: losVelGHz, band9, precision9
       type(VectorValue_T) :: correction9, noise9
       ! configuration file name
       integer, parameter :: strlen = 256
       character(len=strlen) :: signalFileName, configFileName
       
       ! Input file stored in an HDF5 file datafile name in inputlst, Zheng Qu
       character(len=strlen) ::  inputlst='mlscfm_osp.lst', inputH5fileName
       character(len=strlen) ::  outputH5fileName, chbuff, &
        logProcess, logError
       character(len=strlen) :: MLS_CFM_OSP_PATH, MLS_CFM_OUTPUT_PATH, &
        MLS_INPUT_FILE_PATH, TES_MLS_OSP_PATH

       type(Matrix_T) :: jacobian ! stores intermediate Jacobians matrices for CO, O3, Extinction etc. , Zheng Qu 
       integer :: error
       character(len=1024) cbuff

       LOGICAL :: exists ! for testing file/dir existance

!=========== Customized initialization, added by Zheng Qu =====>      
        
       !Get command line arguement for stateVector file path , Zheng Qu
       call getarg(1, inputH5fileName)

        ! now construct log file name from stateVector file name and open them
        ! if stateVector file (inputH5fileName) is like '/[path]/TES_0000010310_0767_004_MLS_1679_2009d041_Version_1_StateVector.hdf'
        ! then log file should be:
        ! TES_0000010310_0767_004_MLS_1679_2009d041_[VersionId]_run.log
        ! TES_0000010310_0767_004_MLS_1679_2009d041_[VersionId]_error.log
        ![VersionId] is CFM version number
        
       call  openLogFiles(trim(adjustl(inputH5fileName)), VersionId, strlen, logProcess, logError)
       
       call runlog(VersionId, 'MLS CFM starts.')
       call runlog(VersionId, 'Preparing for input, output and log files.')
       
       !CALL GET_ENVIRONMENT_VARIABLE('MLS_CFM_OSP_PATH', MLS_CFM_OSP_PATH)
       CALL GET_ENVIRONMENT_VARIABLE('TES_MLS_OSP_PATH', TES_MLS_OSP_PATH)
       MLS_CFM_OSP_PATH = path_join(TES_MLS_OSP_PATH,'MLS_CFM')
       CALL GET_ENVIRONMENT_VARIABLE('MLS_CFM_OUTPUT_PATH', MLS_CFM_OUTPUT_PATH)
       CALL GET_ENVIRONMENT_VARIABLE('MLS_INPUT_FILE_PATH', MLS_INPUT_FILE_PATH )

       ! Testing if directory MLS_CFM_OSP_PATH exists -- please note this is IFORT implementation only
       inquire( DIRECTORY=trim(MLS_CFM_OSP_PATH), exist=exists )
       if (.NOT. exists) then
            call errlog ( moduleName, &
         &(/'Environment variable TES_MLS_OSP_PATH was not properly set?',&
            '$TES_MLS_OSP_PATH/MLS_CFM not found!'/)) 
       endif
       
       
       inquire( DIRECTORY=trim(MLS_CFM_OUTPUT_PATH), exist=exists )
       if (.NOT. exists) then
            call errlog ( moduleName, &
           & (/'Environment variable MLS_CFM_OUTPUT_PATH was not properly set.'/) ) 
       endif

       
       inquire( DIRECTORY=trim(MLS_INPUT_FILE_PATH ), exist=exists )
       if (.NOT. exists) then
            call errlog ( moduleName, &
           & (/'Environment variable MLS_INPUT_FILE_PATH  was not properly set.'/) ) 
       endif
       
       nullify(filedatabase)

       ! initialize hdf5 fortran interface
       CALL h5open_f(error)
       if (error /=0) call errlog ( ModuleName, &
            & (/'Unable to initialize hdf5 fortran interface.'/) )

       call runlog(VersionId, 'Reading OSP files.')
       call read_txtinputdata(MLS_CFM_OSP_PATH, inputlst, spectroscopy, &
        leapsecFile,  antennaPatterns, filterShapes, &
        DACSFilterShapes, pointingGrids, pfaFiles, l2pc , &
            signalFileName, configFileName, &
            vGridStandard37Start ,&
            vGridExtinctionStart ,& 
            hGridStandardVal1 , &
            hGridStandardVal2 , &
            qH2OMinValue , &
            vGridStandard37formula , &
            vGridExtinctionFormula, &
            outputTxtFile, &
            H2O_Sample_Vals)
               
       call runlog(VersionId, 'CFM_MLSSetup -- CFM initialization.')
       call CFM_MLSSetup(signalFileName, configFileName, forwardModelConfigDatabase)

       
       call MLSMessageSetup ( SuppressDebugs=.TRUE., &
        LogFileUnit=errorLogFileUnit, &
        & CrashOnAnyError=.FALSE.  )
       

         
       ! Reading state vector file, it also constructs H2OInput using 
       call runlog(VersionId, 'Reading stateVector file: ' &
          //trim(inputH5fileName))
       call read_H5inputdata (MLS_INPUT_FILE_PATH , trim(inputH5fileName), &
        TemperatureInput, H2OInput , &
        O3Input, SO2Input, HNO3Input, COInput, extinctionV2R3Input, &
        PressureCOInput, PressureStandardInput, &
        REFGPHINPUT, VGRIDREFGPHVALS, MLS_Year, MLS_Day, &
        MLS_ProfileNumber, Tes_run, Tes_scan, TES_sequnce, &
        l1boa, l1brad, l2GP,l2dgm, H2O_Sample_Vals)
   
       !write(*,*)  'TemperatureInput = &   ' 
       !write(*,'(5(E14.7, ",")," &")')  TemperatureInput
       
       !H2OInput = H2OInput*2.
       !write(*,*)  ' H2OInput = & ' 
       !write(*,'(5(E14.7, ",")," &")')  H2OInput  
       
       
       !write(*,*)  'SO2Input= & '
       !write(*,'(5(E14.7, ",")," &")')  SO2Input
       
       !write(*,*)  'HNO3Input= & '
       !write(*,'(5(E14.7, ",")," &")')  HNO3Input
       
       
       !write(*,*)  ' COInput = &  '
       !write(*,'(5(E14.7, ",")," &")')  COInput 
       
       
       !write(*,*)  ' extinctionV2R3Input = &'
       !write(*,'(5(E14.7, ",")," &")')  extinctionV2R3Input 
       

       !write(*,*) 'O3Input= &'
       !write(*,'(5(E14.7, ",")," &")') O3Input 
       
       !write(*,*)  ' PressureCOInput = &'
       !write(*,'(5(E14.7, ",")," &")')  PressureCOInput 
       
       !write(*,*)  'PressureStandardInput = &'
       !write(*,'(5(E14.7, ",")," &")')  PressureStandardInput 
       
       !TemperatureInput = &
       !    (/   299.456_r8,       289.667_r8,       280.912_r8,       272.268_r8,       263.273_r8, &
       !         253.692_r8,       245.674_r8,       231.690_r8,       219.471_r8,       210.275_r8, &
       !         205.069_r8,       201.019_r8,       197.023_r8,       192.212_r8,       192.879_r8, &
       !         196.030_r8,       201.309_r8,       206.541_r8,       208.083_r8,       208.639_r8, &
       !         210.903_r8,       215.163_r8,       219.846_r8,       223.787_r8,       227.240_r8, &
       !         231.812_r8,       237.576_r8,       242.384_r8,       244.396_r8,       243.992_r8, &
       !         242.964_r8,       244.642_r8,       248.765_r8,       255.249_r8,       263.369_r8, &
       !         271.624_r8,       274.746_r8,       271.979_r8,       260.006_r8,       247.311_r8, &
       !         239.643_r8,       239.779_r8,       239.759_r8,       213.600_r8,       201.630_r8, &
       !         190.568_r8,       182.280_r8,       177.571_r8,       175.984_r8,       177.844_r8, &
       !         183.204_r8,       201.225_r8,       264.396_r8,       368.598_r8,       368.598_r8 /)

                
        ! MLS Profile number read from statevector is 0 based
        ! Now we change this number to 1 based.
        ! Reasons:
        ! (From Ming's email 10/8/2012) The reason for me to use '1' based is because in 'mockup.f90', 
        !  all codes created by Paul Wagner and Pranjit assume that 'profile number' starts from '1', 
        !  and all Haley's codes assume the 'MAF number' starts from '0'.  
        
        startProfile = MLS_ProfileNumber + 1 
        endProfile = MLS_ProfileNumber + 1
       

        !File name pattern (for log and output files): 'TES_0000010310_0768_003_MLS_1681_2009d041_'
        
       write(chbuff, "('TES_',I10.10,'_',I4.4,'_',I3.3," // &
          "'_MLS_',I4.4,'_',I4.4,'d',I3.3, '_')") &
          Tes_run, TES_sequnce, &
          Tes_scan, MLS_ProfileNumber, MLS_Year, MLS_Day
          
          
       if (index(trim(inputH5fileName), trim(chbuff)) .le. 0 ) then

           call errlog ( ModuleName, &
               & (/ 'StateVector file', &
               trim(inputH5fileName), &
               'does NOT have the pattern', &
               trim(chbuff) , &
               'which was derived from TES run-seq-scan and MLS profile#, year and day', & 
               'read from the stateVector file itself!' /))
       endif
       
       
       outputH5fileName = trim(chbuff)//VersionId//'_RadianceJacobians.h5'

       
       !log input variables
       call log_inputdata( TES_MLS_OSP_PATH, MLS_CFM_OUTPUT_PATH, &
            MLS_INPUT_FILE_PATH , inputlst, &
            trim(inputH5fileName),spectroscopy, &
            leapsecFile,  antennaPatterns, filterShapes, &
            DACSFilterShapes, pointingGrids, pfaFiles, l2pc , &
            signalFileName, configFileName, &
            vGridStandard37Start ,&
            vGridExtinctionStart ,& 
            hGridStandardVal1 , &
            hGridStandardVal2 , &
            qH2OMinValue , &
            vGridStandard37formula , &
            vGridExtinctionFormula, &
            outputTxtFile, &
            H2O_Sample_Vals, &
            REFGPHINPUT, VGRIDREFGPHVALS, MLS_Year, MLS_Day, &
            MLS_ProfileNumber, Tes_run, Tes_scan, TES_sequnce, &
            l1boa, l1brad, l2GP,l2dgm)

!==========================Start doing the real job!=====================================================
        ! Convert profile to MAF (Major Frame) number, Pranjit Saha
       call runlog(VersionId, 'Converting profile to MAF (Major Frame) number.')
       FileIndex2 = InitializeMLSFile( L1BFile2, content = 'l1boa', &
         & name=trim(l1boa), shortName='L1BOA', type=l_hdf, access=DFACC_RDONLY )
       FileIndex2 = AddFileToDataBase( FileDatabase2, L1BFile2 ) 

       FileIndex2 = InitializeMLSFile( L2GPFile2, content = 'l2gp', &
         & name=trim(L2GP), shortName='L2GP', &
         & type=l_hdfeos, access=DFACC_RDONLY )
       FileIndex2 = AddFileToDataBase( FileDatabase2, L2GPFile2 )

       startL1Maf_new = L2ProfileToL1MAF( startProfile, FileDatabase2 )
       endL1Maf_new = L2ProfileToL1MAF( endProfile, FileDatabase2 )
       
       
       ! Read ptanGHz values from  MLS *DGM* file, Pranjit Saha
       call runlog(VersionId, 'Reading ptanGHz values from  MLS *DGM* file')
       call Read_ptan(trim(l2dgm), 'HDF5', startL1Maf_new, ptanValuesRead)    

       if (minval(ptanValuesRead) .LT. -3.3) then
            write(cbuff, '(F10.3)') minval(ptanValuesRead)
            call MLSMessage ( MLSMSG_Error, moduleName, &
                & 'min(ptanGHZ)='//trim(cbuff)//', <-3.3, MLSCFM failed.' )
       endif
    ! Check ptanGHZ values for -999 for bad cases 
    !  call check_patanGHZ(newer_ptanGHz)
       !========================= Run the forward model ==========================
       call runlog(VersionId, 'Starting running forward model.')
       call runlog(VersionId, 'CFM Initialization--prepare database, vector & templates')
       ! Read L1BOA file
       error = InitializeMLSFile (l1bfile, content='l1boa', name=trim(l1boa), &
       shortName='L1BOA', type=l_hdf, access=DFACC_RDONLY)
       if (error /= 0) call MLSMessage (MLSMSG_Error, moduleName, &
         "Error initializing " // trim(l1boa))
       call mls_openfile (l1bfile, error)
       if (error /= 0) call MLSMessage (MLSMSG_Error, moduleName, &
         "Error opening " // trim(l1boa))
       ! don't care about return value of the following function
       i = AddFileToDatabase(filedatabase, l1bfile)

   ! Create necessary template and vector for ptanGHz already read from file, Pranjit Saha
       newer_qptanGHz = CreateQtyTemplate(l_ptan, startL1Maf_new, endL1Maf_new, &
            & filedatabase=filedatabase, qInstModule='Ghz')
       newer_ptanGHz  = CreateValue4AgileVector(newer_qptanGHz, value=ptanValuesRead)


    ! From Ming: Haley assumes 0-based indexing for L1Maf
       startL1Maf_new = startL1Maf_new - 1
       endL1Maf_new = endL1Maf_new - 1

       ! read MLS input data files
       call runlog(VersionId, 'Reading MLS input data files.')
       
       call Read_Spectroscopy (spectroscopy, 'HDF5')
       call ReadAntennaPatterns (antennaPatterns)
       call ReadFilterShapes(filterShapes)
       call ReadDACSFilterShapes (DACSFilterShapes)
       call ReadPointingGrids (pointingGrids)

       do i = 1, size(pfaFiles)
           call ReadPFAFile (pfaFiles(i))
       end do

       do i = 1, size(l2pc)
           call ReadHDF5L2PC (l2pc(i))
       end do

       vGridStandard37 = CreateVGrid (l_zeta, phyq_pressure, l_logarithmic, &
                                      start=vGridStandard37Start, &
                                      formula=trim(vGridStandard37formula))
       !!vGridStandard37Start = 1000.0d0
       !!vGridStandard37formula = "37:6"
       
       vGridRefGPH = CreateVGrid (l_zeta, phyq_pressure, l_explicit, &
                                  values=vGridRefGPHVals)
       !! vGridRefGPHVals = (/100.0_r8/)
       vGridStandard55 = CreateVGrid (l_zeta, phyq_pressure, l_explicit, &
       values= PressureStandardInput)

       vGridTESCO = CreateVGrid (l_zeta, phyq_pressure, l_explicit, &
       values = PressureCOInput)

       vGridExtinction = CreateVGrid (l_zeta, phyq_pressure, l_logarithmic, &
       start=vGridExtinctionStart, formula=trim(vGridExtinctionFormula))
       !!vGridExtinctionStart=1000.0d0
       !!vGridExtinctionFormula="21:12,14:6,12:3"
       ! Have insetoverlaps, and not single
       hGridStandard = CreateRegularHGrid(GHz, hGridStandardVal1, hGridStandardVal2, &
                                          .true., &
                                          filedatabase, startL1Maf_new, endL1Maf_new)
       !!hGridStandardVal1 = 0.0_r8
       !!hGridStandardVal2 = 1.5_r8
       fGridExtinctionConstant = CreateFGrid(L_IntermediateFrequency, (/0.0_r8/))

       ! Construct state vector
       call runlog(VersionId, 'Constructing the state vector.')
       state = CreateAgileVector(name='state')

       qtemp = CreateQtyTemplate(l_temperature, qName='temperature', &
                                 avgrid=vGridStandard55, ahgrid=hGridStandard)
       temperature = CreateValue4AgileVector(qtemp, value=TemperatureInput)
       call AddValue2Vector(state, temperature)

       qCO = CreateQtyTemplate(l_vmr, avgrid=vGridTESCO, ahgrid=hGridStandard, &
       qMolecule=l_co, qName='CO')
       co = CreateValue4AgileVector(qco, value=COInput)
       call AddValue2Vector(state, co)

       o2 = CreateMLSValue_O2 (vGridStandard37, hGridStandard, qname='O2')
       call AddValue2Vector(state, o2)

       qSO2 = CreateQtyTemplate(l_vmr, avgrid=vGridStandard37, ahgrid=hGridStandard, &
       qMolecule=l_so2, qName='SO2')
       so2 = CreateValue4AgileVector(qSO2, value=SO2Input)
       call AddValue2Vector(state, so2)

       qHNO3 = CreateQtyTemplate(l_vmr, avgrid=vGridStandard37, ahgrid=hGridStandard, &
       qMolecule=l_hno3, qName='HNO3')
       hno3 = CreateValue4AgileVector(qhno3, value=HNO3Input)
       call AddValue2Vector(state, hno3)

       qO3 = CreateQtyTemplate(l_vmr, avgrid=vGridStandard55, ahgrid=hGridStandard, &
       qMolecule=l_o3, qName='O3')
       o3 = CreateValue4AgileVector(qo3, value=O3Input)
       call AddValue2Vector(state, o3)

       qExtinctionv2r3 = CreateQtyTemplate(l_vmr, avgrid=vGridExtinction, &
       ahgrid=hGridStandard, afgrid=fgridExtinctionConstant, qRadiometer="R3", &
       qMolecule=l_extinctionv2)
       extinctionv2r3 = CreateValue4AgileVector(qExtinctionv2r3, value=extinctionV2R3Input)
       call AddValue2Vector(state, extinctionv2r3)

       !Construct stateExtra vector.
       call runlog(VersionId, 'Constructing the stateExtra vector.')
       stateExtra = CreateAgileVector(name='stateExtra')

       ! qPtanGHz = CreateQtyTemplate(l_ptan, startL1Maf_new, endL1Maf_new, filedatabase, &
       ! qInstModule=GHz, qName='ptanGHz')
       ! ptanGHz = CreateValue4AgileVector(qPtanGhz)

       qPhitanGHz = CreateQtyTemplate(l_phitan, startL1Maf_new, endL1Maf_new, qInstModule=GHz, &
       filedatabase=filedatabase, qName='phitanGHz')
       phitanGhz = CreateValue4AgileVector(qPhitanGhz)
       call FillPhitanQuantity(phitanGhz)
       call AddValue2Vector(stateExtra, phitanGhz)

       qRefGPH = CreateQtyTemplate(l_refGPH, avgrid=vGridRefGPH, &
       ahgrid=hGridStandard, qName='refGPH')
       refGPH = CreateValue4AgileVector(qRefGPH, spreadvalue=refGPHInput) ! unit is meter
       call AddValue2Vector(stateExtra, refGPH)

       qGPH = CreateQtyTemplate(l_gph, avgrid=vGridStandard55, ahgrid=hGridStandard, &
                                qName='GPH')
       gph = CreateValue4AgileVector(qGPH)
       call AddValue2Vector(stateExtra, gph)

       qH2O = CreateQtyTemplate(l_vmr, avgrid=vGridStandard55, ahgrid=hGridStandard, &
       qMolecule=l_h2o, qLogBasis=.true., qMinValue=qH2OMinValue, qName='H2O')
       !!qH2OMinValue = 0.1E-6_r8
       h2o = CreateValue4AgileVector(qH2O, value=H2OInput)
       call AddValue2Vector(stateExtra, H2O)

       ! Fill orbit inclination, tangent geocentric altitude with
       ! data from MLS L1B file, and use them, along with other
       ! quantities to calculate ptan
       
       call runlog(VersionId, 'Computing ptan.')
       limbSidebandFraction9L = CreateMLSValue_LSF(9, .false., qname='lsf9L')
       call AddValue2Vector(stateExtra, limbSidebandFraction9L)

       limbSidebandFraction9U = CreateMLSValue_LSF(9, .true., qname='lsf9U')
       call AddValue2Vector(stateExtra, limbSidebandFraction9U)

       elev9L = CreateMLSValue_ElevationOffset(9, .false.)
       call AddValue2Vector(stateExtra, elev9L)

       elev9U = CreateMLSValue_ElevationOffset(9, .true.)
       call AddValue2Vector(stateExtra, elev9U)

       earthReflectivity = CreateMLSValue_EarthReflectivity()
       call AddValue2Vector(stateExtra, earthReflectivity)

       orbitInclination = CreateMLSValue_FromL1BOA (l_orbitInclination, sc, &
       filedatabase, startL1Maf_new, endL1Maf_new)
       call AddValue2Vector(stateExtra, orbitInclination)

       spaceRadiance = CreateMLSValue_SpaceRadiance()
       call AddValue2Vector(stateExtra, spaceRadiance)

       scGeocAlt = CreateMLSValue_FromL1BOA(l_scgeocalt, sc, filedatabase, &
       startL1Maf_new, endL1Maf_new)
       call AddValue2Vector(stateExtra, scGeocAlt)

       tngtGeocAltGHz = CreateMLSValue_FromL1BOA (l_tngtgeocalt, GHz, &
       filedatabase, startL1Maf_new, endL1Maf_new)
       call AddValue2Vector(stateExtra, tngtGeocAltGHz)

       losVelGHz = CreateMLSValue_FromL1BOA(l_losVel, GHz, filedatabase, &
       startL1Maf_new, endL1Maf_new)
       call AddValue2Vector(stateExtra, losVelGHz)

       
       call runlog(VersionId, 'Constructing Radiance vector.')
       radiance = CreateAgileVector(name='simulatedRadiance')

       qband9 = CreateQtyTemplate(l_radiance, startL1Maf_new, endL1Maf_new, filedatabase, &
                                  qSignal="R3:240.B9F:CO", qName='band9')
       band9 = CreateValue4AgileVector(qband9)
       call AddValue2Vector(radiance, band9)

       ! calculate ptan
       !call FillPtanQuantity (ptanGHz, temperature, refGPH, h2o, &
       !orbitInclination, phitanGhz, tngtGeocAltGHz)
       !call AddValue2Vector(stateExtra, ptanGhz)

       ! Add ptan_GHz values read from file in the stateExtra, Pranjit Saha
       
       call runlog(VersionId, 'Adding ptan_GHz values read from file to stateExtra.')
       call AddValue2Vector(stateExtra, newer_ptanGHz)

       ! We no longer need vGrid because the quantity templates have copied it
       call DestroyVGridContents(vGridStandard55)
       call DestroyVGridContents(vGridStandard37)
       call DestroyVGridContents(vGridTESCO)
       ! No long need hGrid, fGrid either
       call DestroyHGridContents(hGridStandard)
       call DestroyFGridContents(fGridExtinctionConstant)

       ! GPH is filled by the forward model

       ! Create jacobian
       jacobian = CreatePlainMatrix(radiance, state)
       
       ! Call the forward model
       call runlog(VersionId, 'Calling the forward model')
       call ForwardModel2 (0, forwardModelConfigDatabase, state, &
                           stateExtra, radiance, jacobian)
						
       ! Write 'state', 'radiance', 'jacobian', 'ptanGHz' in separate files, Pranjit Saha
       call runlog(VersionId, 'Writing state, radiance, ptanGHz etc. to file and concatenate Jacobians.')
       if (outputTxtFile) call Write_To_File1 (state, radiance, jacobian, newer_ptanGHz)

       ! from radiance (template) and Jacobian (template),copy Radiance_Calculated and concatenate Jacobians. Zheng Qu

       !call Write_To_HDF5('mlscfm.old_out.h5', state, radiance, jacobian, newer_ptanGHz)
       call Copy_RadianceJacobian ( radiance, jacobian, newer_ptanGHz, &
        Radiance_Calculated, Jacobian_Concatenated, ptanGHz_save)
        

       !=================== Finish running the forward model =====================
        call runlog(VersionId, 'Finishing running the forward model and doing cleaning up.')
       !== Clean up anything that is not related to reading observed radiance ====
       call DestroyMatrix(jacobian)
       call DestroyAgileVectorContent (state)
       call DestroyAgileVectorContent (stateExtra)
       call Destroy_DACS_Filter_Database
       call Destroy_Filter_Shapes_Database
       call Destroy_Ant_Patterns_Database
       call Destroy_SpectCat_Database
       call Destroy_Line_Database
       call Destroy_Pointing_Grid_Database
       call DestroyL2PCDatabase
       call Destroy_PFADataBase

       ! Destroy all quantity templates that goes in state and stateExtra
       call DestroyQuantityTemplateContents(qtemp)
       call DestroyQuantityTemplateContents(qCO)
       call DestroyQuantityTemplateContents(qSO2)
       call DestroyQuantityTemplateContents(qHNO3)
       call DestroyQuantityTemplateContents(qO3)
       call DestroyQuantityTemplateContents(qExtinctionv2r3)
       call DestroyQuantityTemplateContents(qPtanGHz)
       call DestroyQuantityTemplateContents(qPhitanGHz)
       call DestroyQuantityTemplateContents(qRefGPH)
       call DestroyQuantityTemplateContents(qGPH)
       call DestroyQuantityTemplateContents(qH2O)
       ! even quantities created by CFM subroutines has templates
       call DestroyQuantityTemplateContents(o2%template)
       call DestroyQuantityTemplateContents(limbSidebandFraction9L%template)
       call DestroyQuantityTemplateContents(limbSidebandFraction9U%template)
       call DestroyQuantityTemplateContents(elev9L%template)
       call DestroyQuantityTemplateContents(elev9U%template)
       call DestroyQuantityTemplateContents(earthReflectivity%template)
       call DestroyQuantityTemplateContents(orbitInclination%template)
       call DestroyQuantityTemplateContents(spaceRadiance%template)
       call DestroyQuantityTemplateContents(scGeocAlt%template)
       call DestroyQuantityTemplateContents(tngtGeocAltGHz%template)
       call DestroyQuantityTemplateContents(losVelGHz%template)
       !========================== finish cleaning ===============================

       !====================== Read observed radiance ============================
       ! Open l1brad
       call runlog(VersionId, 'Reading the observed radiance/noise.')
       error = InitializeMLSFile(l1bfile, content='l1brad', &
       name=trim(l1brad), shortName='L1BRAD', type=l_hdf, access=DFACC_RDONLY)
       if (error /= 0) call MLSMessage (MLSMSG_Error, moduleName, &
       "Error initializing " // trim(l1brad))

       call mls_openFile(l1bfile, error)
       if (error /= 0 ) call MLSMessage (MLSMSG_Error, moduleName, &
       "Error opening " // trim(l1brad))

       ! Add it to the filedatabase
       ! AddFileToDatabase doesn't return an error.
       ! I don't care about the return value of AddFileToDatabase,
       ! but Fortran dictate that the return value has to be captured,
       ! so error is being used as a dummy variable.
       error = AddFileToDatabase(filedatabase, l1bfile)

       observed = CreateAgileVector(name='observedRadiance')
       obsPrecision = CreateAgileVector(name='observedRadiancePrecision')

       ! need to read precision before reading quantity
       band9 = CreateValue4AgileVector(qband9)
       call AddValue2Vector(observed, band9)
       ! precision of a quantity has the same template as the quantity
       precision9 = CreateValue4AgileVector(qband9)
       call AddValue2Vector(obsPrecision, precision9)

       ! need to fill precision first
       call FillVectorQuantityFromL1B(precision9, startL1Maf_new, endL1Maf_new, &
       filedatabase, .true.)
       ! then fill the quantity and setting precision at the same time
       call FillVectorQuantityFromL1B(band9, startL1Maf_new, endL1Maf_new, &
       filedatabase, .false., precisionQuantity=precision9)

       ! then we need to get the baseline correction and noise
       qbaseline9 = CreateQtyTemplate(l_l1bmafbaseline, startL1Maf_new, endL1Maf_new, &
       filedatabase, qSignal="R3:240.B9F:CO", qname='baseline 9')

       correction9 = CreateValue4AgileVector (qbaseline9)
       ! the space in ' Baseline' is very important
       call FillVectorQuantityFromL1B(correction9, startL1Maf_new, endL1Maf_new, &
       filedatabase, .false., suffix=' Baseline')

       ! apply correction to the quantity
       call ApplyBaseline (band9, correction9, .false., .false.)

       ! get the noise
       noise9 = CreateValue4AgileVector(qbaseline9) ! same template as baseline
       ! again, the string must match exactly
       call FillVectorQuantityFromL1B(noise9, startL1Maf_new, endL1Maf_new, &
       filedatabase, .false., suffix=' Baseline precision')

       ! apply the noise to the precision
       call ApplyBaseline (precision9, noise9, .true., .false.)
       !call dump(observed, details=1)
       
       !  Filter Radiance and Jacobians using ptanGHz and noise mask values, Zheng Qu
       ! 
       call runlog(VersionId, 'Filtering Radiance and Jacobians using ptanGHz and noise mask values.') 
       call Filter_RadianceJacobian(band9, precision9, ptanGHz_save,&
            Radiance_Calculated, &
            Jacobian_Concatenated , &
            Radiance_Noise, &
            Radiance_Observed, &
            Radiance_Calculated_out, &
            Jacobian_Concatenated_out , &
            Radiance_Noise_out, &
            Radiance_Observed_out)
       ! Write Jacobians and Radiances in HDF5 output file, Zheng Qu
       !call Write2HDF5_old('mlscfm.old_out.h5',band9, precision9)

       ! band9 -> Radiance_Observed , Radiance -> Radiance_Calculated, Presision9 -> Radiance_Noise
       call runlog(VersionId, 'Writing radiance+Jacobians to output file.')
       call Write2HDF5(path_join(MLS_CFM_OUTPUT_PATH,outputH5fileName, novalidation=.TRUE.),&
            Radiance_Calculated_out, &
            Jacobian_Concatenated_out , &
            Radiance_Noise_out, &
            Radiance_Observed_out)
       ! Write 'band9' and 'precision9' in separate text files, Pranjit Saha
       if (outputTxtFile) call Write_To_File2 (band9, precision9)
       ! read output H5 file for validation

       !if (outputTxtFile) call read_H5outputdata(trim(outputH5fileName))
       !==================== Finish reading observed radiance ====================
       call runlog(VersionId, 'Final cleaning up.')
       !============= At this point we don't need L1B file anymore ===============
       do i = 1, size(filedatabase)
           call mls_closefile(filedatabase(i))
       end do


       deallocate(filedatabase)
       
       deallocate( filedatabase2)
       

       !================== Done closing and clean up file objects ================

       diffVector = observed - radiance
       !call dump(diffVector, details=1)

       !===== Clean up calculated and observed radiance and related vectors ======
       call DestroyAgileVectorContent (radiance)

       call DestroyAgileVectorContent (observed) !band9 also destroyed
       call DestroyAgileVectorContent (diffVector)
       call DestroyAgileVectorContent (obsPrecision) !precision9 also destroyed


       call DestroyVectorValueContent (noise9)
       call DestroyVectorValueContent (correction9)

       call DestroyQuantityTemplateContents(qband9)
       call DestroyQuantityTemplateContents(qbaseline9)


       !========================================================

       call CFM_MLSCleanup(forwardModelConfigDatabase)

       
       ! clean up other variables/ vectors/ templates
       
       call DestroyQuantityTemplateContents(newer_qPtanGHz)
       
       call Deallocate_test ( TemperatureInput, &
                & 'Unable to deallocate array for TemperatureInput', ModuleName )

       
       call Deallocate_test ( H2OInput, &
                & 'Unable to deallocate array for H2OInput', ModuleName )

         
         
       call Deallocate_test ( O3Input, &
                & 'Unable to deallocate array for O3Input', ModuleName )
  
         
       call Deallocate_test ( SO2Input, &
                & 'Unable to deallocate array for SO2Input', ModuleName )

                
         
       call Deallocate_test ( HNO3Input, &
                & 'Unable to deallocate array for HNO3Input', ModuleName )

       call Deallocate_test ( COInput, &
                & 'Unable to deallocate array for COInput', ModuleName )

       call Deallocate_test ( extinctionV2R3Input, &
                & 'Unable to deallocate array for extinctionV2R3Input', ModuleName )

       call Deallocate_test ( PressureCOInput, &
                & 'Unable to deallocate array for PressureCOInput', ModuleName )

       call Deallocate_test ( PressureStandardInput, &
                & 'Unable to deallocate array for PressureStandardInput', ModuleName )
         
       call Deallocate_test ( Radiance_Calculated, &
                & 'Unable to deallocate array for Radiance_Calculated', ModuleName )

       call Deallocate_test ( Radiance_Observed, &
                & 'Unable to deallocate array for Radiance_Observed', ModuleName )
  
       call Deallocate_test ( Radiance_Noise, &
                & 'Unable to deallocate array for Radiance_Noise', ModuleName )

  
       call Deallocate_test ( Jacobian_Concatenated, &
                & 'Unable to deallocate array for Jacobian_Concatenated', ModuleName )
      
       call Deallocate_test ( Radiance_Calculated_out, &
                & 'Unable to deallocate array for Radiance_Calculated_out', ModuleName )
  
       call Deallocate_test ( Radiance_Observed_out, &
                & 'Unable to deallocate array for Radiance_Observed_out', ModuleName )
  
       call Deallocate_test ( Radiance_Noise_out, &
                & 'Unable to deallocate array for Radiance_Noise_out', ModuleName )
  
       call Deallocate_test ( Jacobian_Concatenated_out, &
                & 'Unable to deallocate array for Jacobian_Concatenated_out', ModuleName )

       call Deallocate_test ( ptanGHz_save, &
                & 'Unable to deallocate array for ptanGHz_save', ModuleName )

       write(*, '("outputfile=",A)') trim(path_join(MLS_CFM_OUTPUT_PATH, &
            outputH5fileName, novalidation=.TRUE.))
       call runlog(VersionId, 'MLSCFM complated successfully.')
       call runlog(VersionId, "outputfile="// &
            trim(path_join(MLS_CFM_OUTPUT_PATH,outputH5fileName,&
            novalidation=.TRUE.)))
       close(runLogFileUnit) 
       
       close(errorLogFileUnit)
       ! delete errorLog file if empty
       cbuff = '/bin/rm `/usr/bin/find '//trim(logError)// &
        ' -size 0` > /dev/null 2>&1'

       call system(trim(cbuff))
       
end program

        ! $Log$
        ! Revision 1.57  2011/12/24 18:39:13  honghanh
        ! Clean up unused imports and variables
        !
        ! Revision 1.56  2011/12/23 22:56:12  honghanh
        ! Add AutoFillVector call to ForwardModel2 subroutine,
        ! to automatically add and fill isotope ratio in beta group
        ! (according to L2CF configuration file)
        !
        ! Revision 1.55  2011/12/15 16:53:24  honghanh
        ! Correct the name of CreateMLSValue_EarthReflectivity
        !
        ! Revision 1.54  2011/12/14 22:54:18  honghanh
        ! Add timeRange2MafRange method in CFM.
        !
        ! Revision 1.53  2011/11/10 17:07:07  honghanh
        ! Add 'qname' optional argument to a few CreateMLSValue_*
        ! subroutines.
        !
        ! Revision 1.52  2011/11/08 16:13:57  honghanh
        ! Add 'qname' parameter to CreateMLSValue_O2 subroutine.
        !
        ! Revision 1.51  2011/11/03 17:44:32  honghanh
        ! Bug fix in mockup
        !
        ! Revision 1.50  2011/11/02 04:20:16  honghanh
        ! Mockup for new CFM API.
        !
        ! Revision 1.48  2011/03/24 15:16:46  honghanh
        ! Add new interfaces for creating vector and vector values without going through quantity template databases
        !
        ! Revision 1.44  2010/11/03 20:17:01  honghanh
        ! Add name as an optional argument to CreateVector.
        !
        ! Revision 1.42  2010/09/28 14:42:42  honghanh
        ! Add call to forwardModel with jacobian
        !
        ! Revision 1.39  2010/08/06 14:15:06  honghanh
        ! Call dump on diff vector instead of measurement vector.
        !
        ! Revision 1.38  2010/08/05 16:23:03  honghanh
        ! Added Jacobian to forwardModel subroutine
        !
        ! Revision 1.36  2010/07/08 21:39:16  honghanh
        ! Add ApplyBaseline to cfm_fill_m
        !
        ! Revision 1.34  2010/06/29 17:02:47  honghanh
        ! Change the identifier 'fakeChunk' to 'chunk' because
        ! since it is created with ChunkDiVersionIde, it's as real as a chunk
        ! can get.
        !
        ! Revision 1.33  2010/06/29 15:29:33  honghanh
        ! Develop FillPtanQuantity to compute ptan, instead of using
        ! Get2DHydrostaticTangentPressure
        !
        ! Revision 1.32  2010/06/29 02:28:17  honghanh
        ! Change mockup to import functions and literals from CFM module
