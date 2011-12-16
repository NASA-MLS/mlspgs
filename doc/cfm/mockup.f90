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
    use CFM, only: ForwardModelConfig_T 
    use input
    use machine, only: getarg

    implicit none

!---------------------------- RCS Ident Info ------------------------------
    character (len=*), parameter :: ModuleName= &
        "$RCSfile$"
    character (len=*), parameter :: IdParm = &
        "$Id$"
    character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

    integer :: i
    type(ForwardModelConfig_T), dimension(:), pointer :: forwardModelConfigDatabase
    type(MLSFile_T), dimension(:), pointer :: filedatabase
    type(MLSFile_T) :: l1bfile
    type(VGrid_T) :: vGridStandard55, vGridStandard37, vGridTESCO, vGridRefGPH
    type(VGrid_T) :: vGridExtinction
    type(HGrid_T) :: hGridStandard
    type(FGrid_T) :: fGridExtinctionConstant
    type(QuantityTemplate_T) :: qtemp, qCO, qso2, qhno3, qo3, qextinctionv2r3
    type(QuantityTemplate_T) :: qphitanGhz, qrefGPH, qgph, qh2o, qisotoperatioO_18_O
    type(QuantityTemplate_T) :: qisotoperatioO3_ASYM_O_18, qisotoperatioO3_V2
    type(QuantityTemplate_T) :: qisotoperatioHNO3, qisotoperatioCO, qband9
    type(QuantityTemplate_T) :: qptanGHz, qisotoperatioO3, qisotoperatioS_32_O2
    type(QuantityTemplate_T) :: qisotoperatioO3_SYM_O_18, qbaseline9
    type(Vector_T) :: state, stateExtra
    type(Vector_T) :: radiance, diffVector
    type(Vector_T) :: observed, obsPrecision
    type(Vector_T) :: corrections, correctionNoise
    character(len=3) :: GHz = "GHz"
    character(len=2) :: sc = "sc"
    type(VectorValue_T) :: temperature, co, o2, so2, hno3, o3, extinctionv2r3, ptanGHz
    type(VectorValue_T) :: phitanGhz, refGPH, gph, h2o, isotoperatioO_18_O, 
    type(VectorValue_T) :: isotoperatioO3, isotoperatioS_32_O2
    type(VectorValue_T) :: isotoperatioO3_ASYM_O_18, isotoperatioO3_V2
    type(VectorValue_T) :: isotoperatioHNO3, isotoperatioCO, limbSidebandFraction9L
    type(VectorValue_T) :: limbSidebandFraction9U, elev9L, elev9U, earthReflectivity
    type(VectorValue_T) :: orbitInclination, spaceRadiance, scGeocAlt, tngtGeocAltGHz
    type(VectorValue_T) :: losVelGHz, band9, isotoperatioO3_SYM_O_18, precision9
    type(VectorValue_T) :: correction9, noise9
    character(len=256) :: signalFileName, configFileName
    type(Matrix_T) :: jacobian
    integer :: error

    call getarg(1, signalFileName)
    call getarg(2, configFileName)

    nullify(filedatabase)

    call CFM_MLSSetup(signalFileName, configFileName, forwardModelConfigDatabase)

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
             0.000215443_r8, 0.000100_r8, 4.64159e-05_r8, 2.15443e-05_r8, &
             1.00000e-05_r8 /))

    vGridTESCO = CreateVGrid (l_zeta, phyq_pressure, l_explicit, &
    values=(/1012.93_r8,  1000.00_r8,  908.514_r8,  825.402_r8,  749.893_r8,  &
             681.291_r8,  618.966_r8,  562.342_r8, 510.898_r8,  464.160_r8,  421.698_r8, &
             383.117_r8,  348.069_r8,  316.227_r8,  287.298_r8,  261.016_r8, 237.137_r8, &
             215.444_r8,  195.735_r8,  177.829_r8,  161.561_r8,  146.779_r8,  133.352_r8, &
             121.152_r8, 110.069_r8,  100.000_r8,  90.8518_r8,  82.5406_r8,  74.9896_r8, &
             68.1295_r8,  61.8963_r8,  56.2339_r8, 51.0896_r8,  46.4158_r8,  42.1696_r8, &
             38.3119_r8,  34.8071_r8,  31.6229_r8,  28.7299_r8,  26.1017_r8, 23.7136_r8, &
             21.5443_r8,  19.5734_r8,  17.7828_r8,  16.1560_r8,  14.6780_r8,  13.3352_r8, &
             12.1153_r8, 11.0070_r8,  10.0000_r8,  9.08514_r8,  8.25402_r8,  6.81291_r8, &
             5.10898_r8,  4.64160_r8,  3.16227_r8, 2.61016_r8,  2.15443_r8,  1.61560_r8, &
             1.33352_r8,  1.00000_r8,  0.681292_r8,  0.383118_r8,  0.215443_r8, &
             0.100000_r8 /))

    vGridExtinction = CreateVGrid (l_zeta, phyq_pressure, l_logarithmic, &
    start=1000.0d0, formula="21:12,14:6,12:3")

    ! Have insetoverlaps, and not single
    hGridStandard = CreateRegularHGrid(GHz, 0.0_r8, 1.5_r8, .true., &
                                       filedatabase, startL1Maf, endL1Maf)

    fGridExtinctionConstant = CreateFGrid(L_IntermediateFrequency, (/0.0_r8/))

    ! Construct state vector
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
    extinctionv2r3 = CreateValue4AgileVector(qExtinctionv2r3, &
                     value=extinctionV2R3Input)
    call AddValue2Vector(state, extinctionv2r3)

    stateExtra = CreateAgileVector(name='stateExtra')

    qPtanGHz = CreateQtyTemplate(l_ptan, startL1Maf, endL1Maf, filedatabase, &
    qInstModule=GHz, qName='ptanGHz')
    ptanGHz = CreateValue4AgileVector(qPtanGhz)

    qPhitanGHz = CreateQtyTemplate(l_phitan, startL1Maf, endL1Maf, qInstModule=GHz, &
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
    qMolecule=l_h2o, qLogBasis=.true., qMinValue=0.1E-6_r8, qName='H2O')
    h2o = CreateValue4AgileVector(qH2O, value=H2OInput)
    call AddValue2Vector(stateExtra, H2O)

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
    filedatabase, startL1Maf, endL1Maf)
    call AddValue2Vector(stateExtra, orbitInclination)

    spaceRadiance = CreateMLSValue_SpaceRadiance()
    call AddValue2Vector(stateExtra, spaceRadiance)

    scGeocAlt = CreateMLSValue_FromL1BOA(l_scgeocalt, sc, filedatabase, &
    startL1Maf, endL1Maf)
    call AddValue2Vector(stateExtra, scGeocAlt)

    tngtGeocAltGHz = CreateMLSValue_FromL1BOA (l_tngtgeocalt, GHz, &
    filedatabase, startL1Maf, endL1Maf)
    call AddValue2Vector(stateExtra, tngtGeocAltGHz)

    losVelGHz = CreateMLSValue_FromL1BOA(l_losVel, GHz, filedatabase, &
    startL1Maf, endL1Maf)
    call AddValue2Vector(stateExtra, losVelGHz)

    radiance = CreateAgileVector(name='simulatedRadiance')

    qband9 = CreateQtyTemplate(l_radiance, startL1Maf, endL1Maf, filedatabase, &
                               qSignal="R3:240.B9F:CO", qName='band9')
    band9 = CreateValue4AgileVector(qband9)
    call AddValue2Vector(radiance, band9)

    ! calculate ptan
    call FillPtanQuantity (ptanGHz, temperature, refGPH, h2o, &
    orbitInclination, phitanGhz, tngtGeocAltGHz)
    call AddValue2Vector(stateExtra, ptanGhz)

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
!    call ForwardModel2 (0, forwardModelConfigDatabase, state, &
!                        stateExtra, radiance, jacobian)
!    call dump(radiance, details=1)

    !call dump(jacobian, details=3)

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
    call FillVectorQuantityFromL1B(precision9, startL1Maf, endL1Maf, &
    filedatabase, .true.)
    ! then fill the quantity and setting precision at the same time
    call FillVectorQuantityFromL1B(band9, startL1Maf, endL1Maf, &
    filedatabase, .false., precisionQuantity=precision9)

    ! then we need to get the baseline correction and noise
    qbaseline9 = CreateQtyTemplate(l_l1bmafbaseline, startL1Maf, endL1Maf, &
    filedatabase, qSignal="R3:240.B9F:CO", qname='baseline 9')

    correction9 = CreateValue4AgileVector (qbaseline9)
    ! the space in ' Baseline' is very important
    call FillVectorQuantityFromL1B(correction9, startL1Maf, endL1Maf, &
    filedatabase, .false., suffix=' Baseline')

    ! apply correction to the quantity
    call ApplyBaseline (band9, correction9, .false., .false.)

    ! get the noise
    noise9 = CreateValue4AgileVector(qbaseline9) ! same template as baseline
    ! again, the string must match exactly
    call FillVectorQuantityFromL1B(noise9, startL1Maf, endL1Maf, &
    filedatabase, .false., suffix=' Baseline precision')

    ! apply the noise to the precision
    call ApplyBaseline (precision9, noise9, .true., .false.)
    call dump(observed, details=1)

    !==================== Finish reading observed radiance ====================

    !============= At this point we don't need L1B file anymore ===============
    do i = 1, size(filedatabase)
        call mls_closefile(filedatabase(i))
    end do

    deallocate(filedatabase)
    nullify(filedatabase)
    !================== Done closing and clean up file objects ================

    diffVector = observed - radiance
    !call dump(diffVector, details=1)

    !===== Clean up calculated and observed radiance and related vectors ======
    call DestroyAgileVectorContent (radiance)
    call DestroyAgileVectorContent (observed)
    call DestroyAgileVectorContent (diffVector)
    call DestroyAgileVectorContent (obsPrecision)

    call DestroyVectorValueContent (noise9)
    call DestroyVectorValueContent (correction9)

    call DestroyQuantityTemplateContents(qband9)
    call DestroyQuantityTemplateContents(qbaseline9)
    !========================== Done cleaning up ==============================

    call CFM_MLSCleanup(forwardModelConfigDatabase)

end program

