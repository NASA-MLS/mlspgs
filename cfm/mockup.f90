! Copyright 2010, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

! This program is meant to serve as an example, and proof that the CFM library
! is working. Consequently, the design of this program does not follow good
! software design principles. This program should not be used as a part in
! any programs or software suite meant for long-term use.

program mockup

    use CFM
    use highOutput, only: outputNamedValue
    use input
    use machine, only: getarg

    ! To convert MAF to Profile, Pranjit Saha 
    use global_settings, only: L1MAFToL2Profile, L2ProfileToL1MAF
    use MLSL2Options, only: Toolkit
    use MLSHDF5, only: MLS_H5OPEN, MLS_H5CLOSE

    implicit none

!---------------------------- RCS Ident Info ------------------------------
    character (len=*), parameter :: ModuleName= &
        "$RCSfile$"
    character (len=*), parameter :: IdParm = &
        "$Id$"
    character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

    integer :: i

    ! To write values in file, read MLS L1BOA and L2GP files, Pranjit Saha
    CHARACTER(LEN=20), PARAMETER :: FMT1 = "(F16.5)"
    CHARACTER(LEN=20), PARAMETER :: FMT2 = "(E18.8)"
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
    type(QuantityTemplate_T) :: qphitanGhz, qrefGPH, qgph, qisotoperatioO_18_O
    type(QuantityTemplate_T) :: qisotoperatioO3_ASYM_O_18, qisotoperatioO3_V2
    type(QuantityTemplate_T) :: qisotoperatioHNO3, qisotoperatioCO, qband7
    type(QuantityTemplate_T) :: qisotoperatioO3, qisotoperatioS_32_O2, newer_qptanGHz
    type(QuantityTemplate_T) :: qisotoperatioO3_SYM_O_18, qbaseline7
    type(Vector_T) :: state, stateExtra
    type(Vector_T) :: radiance, diffVector
    type(Vector_T) :: observed, obsPrecision
    type(Vector_T) :: corrections, correctionNoise
    character(len=3) :: GHz = "GHz"
    character(len=2) :: sc = "sc"
    type(VectorValue_T) :: temperature, co, o2, so2, hno3, o3, extinctionv2r3, newer_ptanGHz 
    type(VectorValue_T) :: phitanGhz, refGPH, gph, isotoperatioO_18_O, isotoperatioO3
    type(VectorValue_T) :: isotoperatioO3_ASYM_O_18, isotoperatioO3_V2, isotoperatioS_32_O2
    type(VectorValue_T) :: isotoperatioHNO3, isotoperatioCO, limbSidebandFraction7L
    type(VectorValue_T) :: limbSidebandFraction7U, elev7L, elev7U, earthReflectivity
    type(VectorValue_T) :: orbitInclination, spaceRadiance, scGeocAlt, tngtGeocAltGHz
    type(VectorValue_T) :: losVelGHz, band7, isotoperatioO3_SYM_O_18, precision7
    type(VectorValue_T) :: correction7, noise7
    character(len=256) :: signalFileName, configFileName
    type(Matrix_T) :: jacobian
    integer :: error

    call getarg(1, signalFileName)
    call getarg(2, configFileName)
    Toolkit = .false.

    nullify(filedatabase)

    call CFM_MLSSetup(signalFileName, configFileName, forwardModelConfigDatabase)

    ! Convert profile to MAF (Major Frame) number, Pranjit Saha
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
    call Read_ptan(l2dgm, 'HDF5', startL1Maf_new, ptanValuesRead)    

    !========================= Run the forward model ==========================

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

    ! Ming: Haley assumes 0-based indexing for L1Maf
    startL1Maf_new = startL1Maf_new - 1
    endL1Maf_new = endL1Maf_new - 1

    ! read MLS input data file
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
                                   start=1000.0d0, formula="37:6")

    vGridRefGPH = CreateVGrid (l_zeta, phyq_pressure, l_explicit, &
                               values=(/100.0_r8/))
    vGridStandard55 = CreateVGrid (l_zeta, phyq_pressure, l_explicit, &
    values=(/1000.00_r8,  825.404_r8,  681.292_r8,  562.341_r8,  464.159_r8,  &
             383.119_r8,  316.228_r8,  261.016_r8,  215.443_r8,  177.828_r8,  &
             146.780_r8,  121.153_r8,  100.000_r8,  82.5404_r8,  68.1292_r8,  &
             56.2341_r8,  46.4159_r8,  38.3119_r8,  31.6228_r8, 26.1016_r8,  &
             21.5443_r8,  17.7828_r8,  14.6780_r8,  12.1153_r8,  10.0000_r8, &
             8.25404_r8,  6.81292_r8,  5.62341_r8,  4.64159_r8,  3.83119_r8,  &
             3.16228_r8,  2.61016_r8,  2.15443_r8,  1.77828_r8,  1.46780_r8,  &
             1.21153_r8,  1.00000_r8,  0.681292_r8, 0.464159_r8,  0.316228_r8, &
             0.215443_r8, 0.146780_r8,  0.100000_r8,  0.0464159_r8, 0.0215443_r8, &
             0.01000_r8,  0.00464159_r8, 0.00215443_r8, 0.00100_r8, 0.000464159_r8,&
             0.000215443_r8, 0.000100_r8, 4.64159e-05_r8, 2.15443e-05_r8, 1.00000e-05_r8 /))

    vGridTESCO = CreateVGrid (l_zeta, phyq_pressure, l_explicit, &
   values= (/   1.00829E+03_r8,   1.00000E+03_r8,   9.08514E+02_r8,   8.25402E+02_r8,   7.49893E+02_r8, &
            6.81291E+02_r8,   6.18966E+02_r8,   5.62342E+02_r8,   5.10898E+02_r8,   4.64160E+02_r8, &
            4.21698E+02_r8,   3.83117E+02_r8,   3.48069E+02_r8,   3.16227E+02_r8,   2.87298E+02_r8, &
            2.61016E+02_r8,   2.37137E+02_r8,   2.15444E+02_r8,   1.95735E+02_r8,   1.77829E+02_r8, &
            1.61561E+02_r8,   1.46779E+02_r8,   1.33352E+02_r8,   1.21152E+02_r8,   1.10069E+02_r8, &
            1.00000E+02_r8,   9.08518E+01_r8,   8.25406E+01_r8,   7.49896E+01_r8,   6.81295E+01_r8, &
            6.18963E+01_r8,   5.62339E+01_r8,   5.10896E+01_r8,   4.64158E+01_r8,   4.21696E+01_r8, &
            3.83119E+01_r8,   3.48071E+01_r8,   3.16229E+01_r8,   2.87299E+01_r8,   2.61017E+01_r8, &
            2.37136E+01_r8,   2.15443E+01_r8,   1.95734E+01_r8,   1.77828E+01_r8,   1.61560E+01_r8, &
            1.46780E+01_r8,   1.33352E+01_r8,   1.21153E+01_r8,   1.10070E+01_r8,   1.00000E+01_r8, &
            9.08514E+00_r8,   8.25402E+00_r8,   6.81291E+00_r8,   5.10898E+00_r8,   4.64160E+00_r8, &
            3.16227E+00_r8,   2.61016E+00_r8,   2.15443E+00_r8,   1.61560E+00_r8,   1.33352E+00_r8, &
            1.00000E+00_r8,   6.81292E-01_r8,   3.83118E-01_r8,   2.15443E-01_r8,   1.00000E-01_r8 /))

    vGridExtinction = CreateVGrid (l_zeta, phyq_pressure, l_logarithmic, &
    start=1000.0d0, formula="21:12,14:6,12:3")

    ! Have insetoverlaps, and not single
    hGridStandard = CreateRegularHGrid(GHz, 0.0_r8, 1.5_r8, .true., &
                                       filedatabase, startL1Maf_new, endL1Maf_new)

    fGridExtinctionConstant = CreateFGrid(L_IntermediateFrequency, (/0.0_r8/))

    ! Construct state vector
    state = CreateAgileVector(name='state')

    qtemp = CreateQtyTemplate(l_temperature, qName='temperature', &
                              avgrid=vGridStandard55, ahgrid=hGridStandard)
    temperature = CreateValue4AgileVector(qtemp, value=TemperatureInput)
    ! call outputNamedValue ( 'is T qtyTemplate%crossAngles associated?', &
    !   & associated(qtemp%crossAngles) )
    ! call outputNamedValue ( 'is T%qtyTemplate%crossAngles associated?', &
    !   & associated(temperature%template%crossAngles) )
    call AddValue2Vector(state, temperature)

    qCO = CreateQtyTemplate(l_vmr, avgrid=vGridTESCO, ahgrid=hGridStandard, &
    qMolecule=l_co, qName='CO')
    !co = CreateValue4AgileVector(qco, value=COInput)
    !call AddValue2Vector(state, co)

    o2 = CreateMLSValue_O2 (vGridStandard37, hGridStandard)
    call AddValue2Vector(state, o2)

    qSO2 = CreateQtyTemplate(l_vmr, avgrid=vGridStandard37, ahgrid=hGridStandard, &
    qMolecule=l_so2, qName='SO2')
    so2 = CreateValue4AgileVector(qSO2, value=SO2Input)
    call AddValue2Vector(state, so2)

    qHNO3 = CreateQtyTemplate(l_vmr, avgrid=vGridStandard37, ahgrid=hGridStandard, &
    qMolecule=l_hno3, qName='HNO3')
    hno3 = CreateValue4AgileVector(qhno3, value=HNO3Input)
    call AddValue2Vector(state, hno3)

    qO3 = CreateQtyTemplate(l_vmr, avgrid=vGridTESCO, ahgrid=hGridStandard, &
    ! qO3 = CreateQtyTemplate(l_vmr, avgrid=vGridStandard55, ahgrid=hGridStandard, &
    qMolecule=l_o3, qName='O3')
    o3 = CreateValue4AgileVector(qo3, value=O3Input)
    call AddValue2Vector(state, o3)

    qExtinctionv2r3 = CreateQtyTemplate(l_vmr, avgrid=vGridStandard55, &
    ! qExtinctionv2r3 = CreateQtyTemplate(l_vmr, avgrid=vGridExtinction, &
    ahgrid=hGridStandard, afgrid=fgridExtinctionConstant, qRadiometer="R3", &
    qMolecule=l_extinctionv2)
    extinctionv2r3 = CreateValue4AgileVector(qExtinctionv2r3, value=extinctionV2R3Input)
    call AddValue2Vector(state, extinctionv2r3)
    call outputnamedValue( 'size(state)', size(state%quantities) )
    call Dump( state, details=-1 )

    stateExtra = CreateAgileVector(name='stateExtra')

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

    qIsotoperatioO_18_O = CreateQtyTemplate(l_isotoperatio, &
    qName='isotoperatioO_18_O', qmolecule=l_o_18_o)
    isotoperatioO_18_O = CreateValue4AgileVector(qIsotoperatioO_18_O, &
    spreadvalue=0.00409000_r8)
    call AddValue2Vector(stateExtra, isotoperatioO_18_O)

    qIsotoperatioO3 = CreateQtyTemplate(l_isotoperatio, qName='isotoperatioO3', &
    qmolecule=l_o3)
    isotoperatioO3 = CreateValue4AgileVector(qIsotoperatioO3, &
    spreadvalue=0.99290103_r8)
    call AddValue2Vector(stateExtra, isotoperatioO3)

    qIsotoperatioO3_ASYM_O_18 = CreateQtyTemplate(l_isotoperatio, &
    qmolecule=l_o3_asym_o_18, qName='isotoperatioO3_ASYM_O_18')
    isotoperatioO3_ASYM_O_18 = CreateValue4AgileVector(qIsotoperatioO3_ASYM_O_18, &
    spreadvalue=0.00398194_r8)
    call AddValue2Vector(stateExtra, isotoperatioO3_ASYM_O_18)

    qIsotoperatioO3_SYM_O_18 = CreateQtyTemplate(l_isotoperatio, &
    qmolecule=l_o3_sym_o_18, qName='isotoperatioO3_SYM_O_18')
    isotoperatioO3_SYM_O_18 = CreateValue4AgileVector(qIsotoperatioO3_SYM_O_18, &
    spreadvalue=0.00199097_r8)
    call AddValue2Vector(stateExtra, isotoperatioO3_SYM_O_18)

    qIsotoperatioO3_V2 = CreateQtyTemplate(l_isotoperatio, qmolecule=l_o3_v2, &
    qName='isotoperatioO3_V2')
    isotoperatioO3_V2 = CreateValue4AgileVector(qIsotoperatioO3_V2, &
    spreadvalue=0.99290103_r8)
    call AddValue2Vector(stateExtra, isotoperatioO3_V2)

    qIsotoperatioS_32_O2 = CreateQtyTemplate(l_isotoperatio, qmolecule=l_s_32_o2, &
    qName='isotoperatioS_32_O2')
    isotoperatioS_32_O2 = CreateValue4AgileVector(qIsotoperatioS_32_O2, &
    spreadvalue=0.94568002_r8)
    call AddValue2Vector(stateExtra, isotoperatioS_32_O2)

    qIsotoperatioHNO3 = CreateQtyTemplate(l_isotoperatio, qmolecule=l_hno3, &
    qName='isotoperatioHNO3')
    isotoperatioHNO3 = CreateValue4AgileVector(qIsotoperatioHNO3, &
    spreadvalue=0.98910999_r8)
    call AddValue2Vector(stateExtra, isotoperatioHNO3)

    qIsotoperatioCO = CreateQtyTemplate(l_isotoperatio, qmolecule=l_co, &
    qName='isotoperatioCO')
    isotoperatioCO = CreateValue4AgileVector(qIsotoperatioCO, &
    spreadvalue=0.98654002_r8)
    call AddValue2Vector(stateExtra, isotoperatioCO)

    ! Fill orbit inclination, tangent geocentric altitude with
    ! data from MLS L1B file, and use them, along with other
    ! quantities to calculate ptan
    limbSidebandFraction7L = CreateMLSValue_LSF(7, .false.)
    call AddValue2Vector(stateExtra, limbSidebandFraction7L)

    limbSidebandFraction7U = CreateMLSValue_LSF(7, .true.)
    call AddValue2Vector(stateExtra, limbSidebandFraction7U)

    elev7L = CreateMLSValue_ElevationOffset(7, .false.)
    call AddValue2Vector(stateExtra, elev7L)

    elev7U = CreateMLSValue_ElevationOffset(7, .true.)
    call AddValue2Vector(stateExtra, elev7U)

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

    radiance = CreateAgileVector(name='simulatedRadiance')

    qband7 = CreateQtyTemplate(l_radiance, startL1Maf_new, endL1Maf_new, filedatabase, &
                               qSignal="R3:240.B7F:O3", qName='band7')
    band7 = CreateValue4AgileVector(qband7)
    call AddValue2Vector(radiance, band7)

    ! Add ptan_GHz values read from file in the stateExtra, Pranjit Saha
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
    call ForwardModel2 (0, forwardModelConfigDatabase, state, &
                        stateExtra, radiance, jacobian)

    ! Write 'state', 'radiance', 'jacobian', 'ptanGHz' in separate files, Pranjit Saha
    call Write_To_File1 (state, radiance, jacobian, newer_ptanGHz)

    !=================== Finish running the forward model =====================

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
    if ( .false. ) then
    call DestroyQuantityTemplateContents(qtemp)
    call DestroyQuantityTemplateContents(qCO)
    call DestroyQuantityTemplateContents(qSO2)
    call DestroyQuantityTemplateContents(qHNO3)
    call DestroyQuantityTemplateContents(qO3)
    call DestroyQuantityTemplateContents(qExtinctionv2r3)
    call DestroyQuantityTemplateContents(newer_qPtanGHz)
    call DestroyQuantityTemplateContents(qPhitanGHz)
    call DestroyQuantityTemplateContents(qRefGPH)
    call DestroyQuantityTemplateContents(qGPH)
    call DestroyQuantityTemplateContents(qIsotoperatioO_18_O)
    call DestroyQuantityTemplateContents(qIsotoperatioO3)
    call DestroyQuantityTemplateContents(qIsotoperatioO3_ASYM_O_18)
    call DestroyQuantityTemplateContents(qIsotoperatioO3_SYM_O_18)
    call DestroyQuantityTemplateContents(qIsotoperatioO3_V2)
    call DestroyQuantityTemplateContents(qIsotoperatioS_32_O2)
    call DestroyQuantityTemplateContents(qIsotoperatioHNO3)
    call DestroyQuantityTemplateContents(qIsotoperatioCO)
    ! even quantities created by CFM subroutines has templates
    call DestroyQuantityTemplateContents(o2%template)
    call DestroyQuantityTemplateContents(limbSidebandFraction7L%template)
    call DestroyQuantityTemplateContents(limbSidebandFraction7U%template)
    call DestroyQuantityTemplateContents(elev7L%template)
    call DestroyQuantityTemplateContents(elev7U%template)
    call DestroyQuantityTemplateContents(earthReflectivity%template)
    call DestroyQuantityTemplateContents(orbitInclination%template)
    call DestroyQuantityTemplateContents(spaceRadiance%template)
    call DestroyQuantityTemplateContents(scGeocAlt%template)
    call DestroyQuantityTemplateContents(tngtGeocAltGHz%template)
    call DestroyQuantityTemplateContents(losVelGHz%template)
    endif
    ! stop
    !========================== finish cleaning ===============================

    !====================== Read observed radiance ============================
    ! Open l1brad
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
    band7 = CreateValue4AgileVector(qband7)
    call AddValue2Vector(observed, band7)
    ! precision of a quantity has the same template as the quantity
    precision7 = CreateValue4AgileVector(qband7)
    call AddValue2Vector(obsPrecision, precision7)

    ! need to fill precision first
    call FillVectorQuantityFromL1B(precision7, startL1Maf_new, endL1Maf_new, &
    filedatabase, .true.)
    ! then fill the quantity and setting precision at the same time
    call FillVectorQuantityFromL1B(band7, startL1Maf_new, endL1Maf_new, &
    filedatabase, .false., precisionQuantity=precision7)

    ! then we need to get the baseline correction and noise
    qbaseline7 = CreateQtyTemplate(l_l1bmafbaseline, startL1Maf_new, endL1Maf_new, &
    filedatabase, qSignal="R3:240.B7F:O3", qname='baseline 7')

    correction7 = CreateValue4AgileVector (qbaseline7)
    ! the space in ' Baseline' is very important
    call FillVectorQuantityFromL1B(correction7, startL1Maf_new, endL1Maf_new, &
    filedatabase, .false., suffix=' Baseline')

    ! apply correction to the quantity
    call ApplyBaseline (band7, correction7, .false., .false.)

    ! get the noise
    noise7 = CreateValue4AgileVector(qbaseline7) ! same template as baseline
    ! again, the string must match exactly
    call FillVectorQuantityFromL1B(noise7, startL1Maf_new, endL1Maf_new, &
    filedatabase, .false., suffix=' Baseline precision')

    ! apply the noise to the precision
    call ApplyBaseline (precision7, noise7, .true., .false.)

    ! Write 'band7' and 'precision7' in separate text files, Pranjit Saha
    call Write_To_File2 (band7, precision7)

    !==================== Finish reading observed radiance ====================

    !============= At this point we don't need L1B file anymore ===============
    do i = 1, size(filedatabase)
        call mls_closefile(filedatabase(i))
    end do

    deallocate(filedatabase)
    nullify(filedatabase)
    !================== Done closing and clean up file objects ================

    diffVector = observed - radiance

    !===== Clean up calculated and observed radiance and related vectors ======
    call DestroyAgileVectorContent (radiance)
    call DestroyAgileVectorContent (observed)
    call DestroyAgileVectorContent (diffVector)
    call DestroyAgileVectorContent (obsPrecision)

    call DestroyVectorValueContent (noise7)
    call DestroyVectorValueContent (correction7)

    call DestroyQuantityTemplateContents(qband7)
    call DestroyQuantityTemplateContents(qbaseline7)
    !========================== Done cleaning up ==============================

    call CFM_MLSCleanup(forwardModelConfigDatabase)

end program

! $Log$
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
! since it is created with ChunkDivide, it's as real as a chunk
! can get.
!
! Revision 1.33  2010/06/29 15:29:33  honghanh
! Develop FillPtanQuantity to compute ptan, instead of using
! Get2DHydrostaticTangentPressure
!
! Revision 1.32  2010/06/29 02:28:17  honghanh
! Change mockup to import functions and literals from CFM module
