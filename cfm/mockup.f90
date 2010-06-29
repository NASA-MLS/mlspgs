program mockup
   use input   ! Provides hard-coded input for testing purposes only

   use CFM_MLSSetup_m, only: CFM_MLSCleanup, MLSChunk_T, CFM_MLSSetup
   use CFM_IO_M, only: Read_Spectroscopy, ReadDACSFilterShapes, &
                     ReadAntennaPatterns, ReadFilterShapes, &
                     ReadPointingGrids, ReadPFAFile, ReadHDF5L2PC, &
                     Destroy_DACS_Filter_Database, &
                     Destroy_Filter_Shapes_Database, &
                     Destroy_Ant_Patterns_Database, &
                     Destroy_SpectCat_Database, &
                     Destroy_Line_Database, ReadObservedRadiances, &
                     Destroy_Pointing_Grid_Database, &
                     DestroyL2PCDatabase, Destroy_PFADataBase
   use CFM_VGrid_m, only: CreateVGrid, DestroyVGridContents, &
                        VGrid_T, Dump
   use CFM_HGrid_m, only: CreateRegularHGrid, HGrid_T, &
                        DestroyHGridContents, Dump
   use CFM_FGrid_m, only: CreateFGrid, FGrid_T, DestroyFGridContents, &
                        Dump
   use CFM_QuantityTemplate_m, only: CreateQtyTemplate, Dump, &
                        AddQuantityTemplateToDatabase, &
                        DestroyQuantityTemplateDatabase, &
                        QuantityTemplate_T
   use CFM_VectorTemplate_m, only: CreateVectorTemplate, Dump, &
                        VectorTemplate_T, DestroyVectorTemplateInfo
   use CFM_Vector_m, only: CreateVector, Dump, &
                         Vector_T, VectorValue_T, &
                         DestroyVectorInfo, GetVectorQtyByTemplateIndex
   use CFM_Fill_m, only: ExplicitFillVectorQuantity, FillPhitanQuantity, &
                         FillVectorQuantityFromL1B, SpreadFillVectorQuantity
   use CFM_FWDMDL_M, only: ForwardModel, FORWARDMODELSTATUS_T, &
                         ForwardModelConfig_T
   use MLSCommon, only: MLSFile_T, r8
   use Init_tables_module, only: l_logarithmic, l_zeta, l_temperature, &
                                 L_IntermediateFrequency, l_vmr, l_gph, &
                                 l_ptan, l_radiance, l_orbitInclination, &
                                 l_tngtgeodalt, l_tngtgeocalt, l_o3, l_refGPH, &
                                 phyq_pressure, phyq_angle, l_h2o, l_explicit, &
                                 l_phitan
   use MLSFiles, only: GetMLSFileByType, InitializeMLSFile, mls_openFile, &
                       AddFileToDatabase
   use ScanModelModule, only: Get2DHydrostaticTangentPressure
   use machine, only: getarg
   use Hdf, only: DFACC_RDONLY
   use Intrinsic, only: l_hdf
   use MLSMessageModule, only: MLSMessage, MLSMSG_Error

   implicit none

!---------------------------- RCS Ident Info ------------------------------
   character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
   character (len=*), parameter :: IdParm = &
       "$Id$"
   character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

   integer :: error, i
   type(ForwardModelConfig_T), dimension(:), pointer :: forwardModelConfigDatabase
   type(MLSFile_T), dimension(:), pointer :: filedatabase
   type(MLSChunk_T) :: fakeChunk
   type(VGrid_T) :: vGridStandard, vGridRefGPH
   type(HGrid_T) :: hGridStandard
   type(FGrid_T) :: fGridStandard
   type(QuantityTemplate_T) :: temperature, GPH, H2O, O3, ptanGHz, band7, &
                               geodAltitude, orbincl, geocAlt, refGPH, band2, &
                               band8, phitanGHz
   type(QuantityTemplate_T), dimension(:), pointer :: qtyTemplates
   type(VectorTemplate_T) :: stateTemplate, measurementTemplate
   type(Vector_T) :: state, measurement, stateExtra, observed, obsPrecision
   character(len=3) :: GHz = "GHz"
   character(len=2) :: sc = "sc"
   integer :: stateSelected(10), measurementSelected(3)
   type(VectorValue_T), pointer :: quantity, h2o_vv, orbincl_vv, geocAlt_vv, &
                          ptanG_vv, temperature_vv, refGPH_vv, precQty
   integer :: temperature_index, h2o_index, band2_index
   integer :: o3_index, ptanGHz_index, band7_index, phitanGHz_index
   integer :: geodAlt_index, orbincl_index, gph_index
   integer :: geocAlt_index, band8_index, refGPH_index
   character(len=256) :: signalFileName, configFileName
   type (MLSFile_T) :: l1bfile

   call getarg(1, signalFileName)
   call getarg(2, configFileName)

   call CFM_MLSSetup(startTime, endTime, l1boa, leapsecFile, signalFileName, &
   configFileName, filedatabase, qtyTemplates, fakeChunk, &
   forwardModelConfigDatabase, stateExtra)

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

   vGridStandard = CreateVGrid (l_zeta, phyq_pressure, l_logarithmic, &
                                start=1000.0d0, formula="37:6")
   !call dump(vGridStandard, details=2)

   vGridRefGPH = CreateVGrid (l_zeta, phyq_pressure, l_explicit, &
                              values=(/100.0_r8/))

   ! Have insetoverlaps, and not single
   hGridStandard = CreateRegularHGrid(GHz, 0.0_r8, 1.5_r8, .true., &
        filedatabase, fakeChunk)
   !call dump(hGridStandard)

   fGridStandard = CreateFGrid(L_IntermediateFrequency, (/0.0_r8/))
   !call dump(fGridStandard)

   temperature = CreateQtyTemplate(l_temperature, filedatabase=filedatabase, &
                                   chunk=fakeChunk, &
                                   avgrid=vGridStandard, ahgrid=hGridStandard)
   GPH = CreateQtyTemplate(l_gph, filedatabase=filedatabase, chunk=fakeChunk, &
                           avgrid=vGridStandard, ahgrid=hGridStandard)
   O3 = CreateQtyTemplate(l_vmr, filedatabase=filedatabase, chunk=fakeChunk, &
                          avgrid=vGridStandard, ahgrid=hGridStandard, qMolecule=l_o3)
   H2O = CreateQtyTemplate(l_vmr, filedatabase=filedatabase, chunk=fakeChunk, &
                           avgrid=vGridStandard, ahgrid=hGridStandard, qMolecule=l_h2o, &
                           qLogBasis=.true., qMinValue=0.1E-6_r8)
   ptanGHz = CreateQtyTemplate(l_ptan, filedatabase=filedatabase, &
                               chunk=fakeChunk, qInstModule=GHz)
   ! band 2,7,8 is the band whose radiances are to be computed
   ! see CFM document for a list of signals corresponding to bands
   band7 = CreateQtyTemplate(l_radiance, filedatabase=filedatabase, chunk=fakeChunk, &
                             qSignal="R3:240.B7F:O3")
   band2 = CreateQtyTemplate(l_radiance, filedatabase=filedatabase, chunk=fakeChunk, &
                             qSignal="R2:190.B2F:H2O")
   band8 = CreateQtyTemplate(l_radiance, filedatabase=filedatabase, chunk=fakeChunk, &
                             qSignal="R3:240.B8F:PT")
   geodAltitude = CreateQtyTemplate(l_tngtgeodalt, filedatabase=filedatabase, &
                                    chunk=fakeChunk, qInstModule=GHz)
   geocAlt = CreateQtyTemplate(l_tngtgeocalt, filedatabase=filedatabase, &
                               chunk=fakeChunk, qInstModule=GHz)
   orbincl = CreateQtyTemplate(l_orbitInclination, filedatabase=filedatabase, &
                               chunk=fakeChunk, qInstModule=sc)
   refGPH = CreateQtyTemplate(l_refGPH, avgrid=vGridRefGPH, ahgrid=hGridStandard)
   phitanGHz = CreateQtyTemplate(l_phitan, qInstModule="GHz", filedatabase=filedatabase, &
                                 chunk=fakeChunk)

   temperature_index = AddQuantityTemplateToDatabase(qtyTemplates, temperature)
   gph_index = AddQuantityTemplateToDatabase(qtyTemplates, GPH)
   o3_index = AddQuantityTemplateToDatabase(qtyTemplates, O3)
   h2o_index = AddQuantityTemplateToDatabase(qtyTemplates, H2O)
   ptanGHz_index = AddQuantityTemplateToDatabase(qtyTemplates, ptanGHz)
   geodAlt_index = AddQuantityTemplateToDatabase(qtyTemplates, geodAltitude)
   geocAlt_index = AddQuantityTemplatetoDatabase(qtyTemplates, geocAlt)
   orbincl_index = AddQuantityTemplateToDatabase(qtyTemplates, orbIncl)
   refGPH_index = AddQuantityTemplateToDatabase(qtyTemplates, refGPH)
   phitanGHz_index = AddQuantityTemplateToDatabase(qtyTemplates, phitanGHz)
   band7_index = AddQuantityTemplateToDatabase(qtyTemplates, band7)
   band2_index = AddQuantityTemplateToDatabase(qtyTemplates, band2)
   band8_index = AddQuantityTemplateToDatabase(qtyTemplates, band8)

   ! We no longer need vGrid because the quantity templates have copied it
   call DestroyVGridContents(vGridStandard)
   ! No long need hGrid, fGrid either
   call DestroyHGridContents(hGridStandard)
   call DestroyFGridContents(fGridStandard)

!   call dump(ptanGHz, details=2)
!   call dump(ptanTHz, details=2)
!   call dump(temperature, details=2)
!   call dump(GPH, details=2)
!   call dump(O3, details=2)
!   call dump(H2O, details=2)
!   call dump(band7, details=2)
!   call dump(geodAltitude, details=2)

   ! The numbers are the order that quantities template were added
   stateSelected = (/temperature_index,o3_index,h2o_index, phitanGHz_index, &
                     ptanGHz_index,geodAlt_index, geocAlt_index, &
                     orbincl_index, gph_index, refGPH_index/)
   stateTemplate = CreateVectorTemplate(qtyTemplates, stateSelected)
!   call dump(stateTemplate, details=2, quantities=qtyTemplates)

   measurementSelected = (/band7_index, band2_index, band8_index/)
   measurementTemplate = CreateVectorTemplate(qtyTemplates, measurementSelected)
!   call dump(measurementTemplate, details=2, quantities=qtyTemplates)

   state = CreateVector(stateTemplate, qtyTemplates)
   measurement = CreateVector(measurementTemplate, qtyTemplates)

   refGPH_vv => GetVectorQtyByTemplateIndex (state, refGPH_index) ! refGPH
   call SpreadFillVectorQuantity (refGPH_vv, refGPHInput) ! unit is meter

   ! supply temperature, GPH, H2O, and O3 data
   temperature_vv => GetVectorQtyByTemplateIndex(state, temperature_index)
   call ExplicitFillVectorQuantity(temperature_vv, TemperatureInput)

   h2o_vv => GetVectorQtyByTemplateIndex(state, h2o_index)
   call ExplicitFillVectorQuantity(h2o_vv, H2OInput)

   quantity => GetVectorQtyByTemplateIndex(state, o3_index)
   call ExplicitFillVectorQuantity(quantity, O3Input)

   quantity => GetVectorQtyByTemplateIndex(state, geodAlt_index)
   call FillVectorQuantityFromL1B(quantity, fakeChunk, filedatabase, &
      .false.)

   ! Fill orbit inclination, tangent geocentric altitude with
   ! data from MLS L1B file, and use them, along with other
   ! quantities to calculate ptan
   orbincl_vv => GetVectorQtyByTemplateIndex (state, orbincl_index)
   call FillVectorQuantityFromL1B(orbincl_vv, fakeChunk, filedatabase, &
      .false.)

   geocAlt_vv => GetVectorQtyByTemplateIndex (state, geocAlt_index)
   call FillVectorQuantityFromL1B(geocAlt_vv, fakeChunk, filedatabase, &
      .false.)

   ptanG_vv => GetVectorQtyByTemplateIndex (state, ptanGHz_index)

   quantity => GetVectorQtyByTemplateIndex (state, phitanGHz_index)
   call FillPhitanQuantity(quantity)

   ! calculate ptan
   !call dump(temperature_vv, details=3)
   call Get2DHydrostaticTangentPressure(ptanG_vv, temperature_vv, refGPH_vv, &
             h2o_vv, orbincl_vv, quantity, geocAlt_vv, 4, 0.0_r8, phyq_angle)

   ! GPH is filled by the forward model

   ! Call the forward model
   call ForwardModel (fakeChunk, forwardModelConfigDatabase, state, &
                      stateExtra, measurement)

   call dump(measurement, details=3)

!   ! Create an identical vector as simulated radiance vector for observed radiances
!   observed = CreateVector(measurementTemplate, qtyTemplates)
!   obsPrecision = CreateVector(measurementTemplate, qtyTemplates)
!
!   ! Open l1brad
!   error = InitializeMLSFile(l1bfile, content='l1brad', &
!   name=trim(l1brad), shortName='L1BRAD', type=l_hdf, access=DFACC_RDONLY)
!   if (error /= 0) &
!      call MLSMessage (MLSMSG_Error, moduleName, &
!      "Error initializing " // trim(l1brad))
!
!   call mls_openFile(l1bfile, error)
!   if (error /= 0 ) &
!      call MLSMessage (MLSMSG_Error, moduleName, &
!      "Error opening " // trim(l1brad))
!
!   ! Add it to the filedatabase
!   error = AddFileToDatabase(filedatabase, l1bfile)
!
!   ! Fill band 2,7,8
!   quantity => GetVectorQtyByTemplateIndex(observed, band2_index)
!   precQty => GetVectorQtyByTemplateIndex(obsPrecision, band2_index)
!
!   ! However, only band 9 and 25 have BOMask=1
!   ! Because these bands have bright object status read from L1BOA file,
!   ! we always have to have L1BOA file in the filedatabase.
!   ! You have to fill the precision quantity first
!   call FillVectorQuantityFromL1B(precQty, fakeChunk, filedatabase, &
!   .true.)
!   call FillVectorQuantityFromL1B(quantity, fakeChunk, filedatabase, &
!   .false., precisionQuantity=precQty)
!
!   quantity => GetVectorQtyByTemplateIndex(observed, band7_index)
!   precQty => GetVectorQtyByTemplateIndex(obsPrecision, band7_index)
!   call FillVectorQuantityFromL1B(precQty, fakeChunk, filedatabase, &
!   .true.)
!   call FillVectorQuantityFromL1B(quantity, fakeChunk, filedatabase, &
!   .false., precisionQuantity=precQty)
!
!   quantity => GetVectorQtyByTemplateIndex(observed, band8_index)
!   precQty => GetVectorQtyByTemplateIndex(obsPrecision, band8_index)
!   call FillVectorQuantityFromL1B(precQty, fakeChunk, filedatabase, &
!   .true.)
!   call FillVectorQuantityFromL1B(quantity, fakeChunk, filedatabase, &
!   .false., precisionQuantity=precQty)
!
!   ! CFM_MLSCleanup will close all files in the filedatabase
!
!   call dump(observed, details=3)
!
!   ! Clean up allocated memory for creating observed radiance vector
!   call DestroyVectorInfo(observed)
!   call DestroyVectorInfo(obsPrecision)

   ! Clean up memory
   call DestroyVectorInfo (state)
   call DestroyVectorInfo (measurement)
   call DestroyVectorTemplateInfo(stateTemplate)
   call DestroyVectorTemplateInfo(measurementTemplate)
   call Destroy_DACS_Filter_Database
   call Destroy_Filter_Shapes_Database
   call Destroy_Ant_Patterns_Database
   call Destroy_SpectCat_Database
   call Destroy_Line_Database
   call Destroy_Pointing_Grid_Database
   call DestroyL2PCDatabase
   call Destroy_PFADataBase

   call CFM_MLSCleanup(filedatabase, qtyTemplates, &
   forwardModelConfigDatabase, stateExtra)

end program
