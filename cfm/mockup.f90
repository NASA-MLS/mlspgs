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
   use input   ! Provides hard-coded input for testing purposes only

   use CFM, only: & ! Alphabetized, case-insensitive list (underscores are ignored though)
      AddFileToDatabase, AddQuantityTemplateToDatabase, ApplyBaseline, &
      CFM_MLSCleanup, CFM_MLSSetup, CreateFGrid, CreatePlainMatrix, &
      CreateQtyTemplate, CreateRegularHGrid, CreateVector, &
      CreateVectorTemplate, CreateVGrid, &
      Destroy_Ant_Patterns_Database, DestroyFGridContents, &
      Destroy_DACS_Filter_Database, Destroy_Filter_Shapes_Database, &
      DestroyHGridContents, DestroyL2PCDatabase, Destroy_Line_Database, &
      Destroy_PFADataBase, Destroy_Pointing_Grid_Database, &
      DestroyQuantityTemplateDatabase, Destroy_SpectCat_Database, &
      DestroyVectorInfo, DestroyVectorTemplateInfo, DestroyVGridContents, &
      DFACC_RDONLY, Dump, &
      ExplicitFillVectorQuantity, &
      FGrid_T, FillPhitanQuantity, FillPtanQuantity, &
      FillVectorQuantityFromL1B, ForwardModel, ForwardModelConfig_T, &
      FORWARDMODELSTATUS_T, &
      GetMLSFileByType, GetVectorQtyByTemplateIndex, &
      HGrid_T, &
      InitializeMLSFile, &
      l_explicit, l_gph, l_h2o, l_hdf, L_IntermediateFrequency, &
      l_logarithmic, l_o3, l_orbitInclination, l_phitan, l_ptan, &
      l_radiance, l_refGPH, l_temperature, l_tngtgeocalt, &
      l_tngtgeodalt, l_vmr, l_zeta, l_l1bMAFBaseline, &
      Matrix_T, MLSChunk_T, MLSFile_T, MLSMessage, MLSMSG_Error, &
      mls_openFile, &
      operator(+), operator(-), &
      phyq_angle, phyq_pressure, &
      QuantityTemplate_T, &
      r8, ReadAntennaPatterns, ReadDACSFilterShapes, ReadFilterShapes, &
      ReadHDF5L2PC, ReadPFAFile, ReadPointingGrids, Read_Spectroscopy, &
      SpreadFillVectorQuantity, &
      Vector_T, VectorTemplate_T, VectorValue_T, VGrid_T
   use machine, only: getarg

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
   type(MLSChunk_T) :: chunk
   type(VGrid_T) :: vGridStandard, vGridRefGPH
   type(HGrid_T) :: hGridStandard
   type(FGrid_T) :: fGridStandard
   type(QuantityTemplate_T) :: temperature, GPH, H2O, O3, ptanGHz, band7, &
                               geodAltitude, orbincl, geocAlt, refGPH, band2, &
                               band8, phitanGHz, baseline7, baseline2, baseline8
   type(QuantityTemplate_T), dimension(:), pointer :: qtyTemplates
   type(VectorTemplate_T) :: stateTemplate, measurementTemplate, correctionTemplate
   type(Vector_T) :: state, measurement, stateExtra, observed, obsPrecision
   type(Vector_T) :: diffVector, corrections, correctionNoise
   character(len=3) :: GHz = "GHz"
   character(len=2) :: sc = "sc"
   integer :: stateSelected(10), measurementSelected(3), baselineSelected(3)
   type(VectorValue_T), pointer :: quantity, h2o_vv, orbincl_vv, geocAlt_vv, &
                          ptanG_vv, temperature_vv, refGPH_vv, precQty
   type(VectorValue_T), pointer :: band2L1BMAFBaseline, band7L1BMAFBaseline, band8L1BMAFBaseline
   integer :: temperature_index, h2o_index, band2_index
   integer :: o3_index, ptanGHz_index, band7_index, phitanGHz_index
   integer :: geodAlt_index, orbincl_index, gph_index
   integer :: geocAlt_index, band8_index, refGPH_index
   integer :: baseline7_index, baseline2_index, baseline8_index
   character(len=256) :: signalFileName, configFileName
   type (MLSFile_T) :: l1bfile
   type(Matrix_T) :: jacobian

   call getarg(1, signalFileName)
   call getarg(2, configFileName)

   call CFM_MLSSetup(startTime, endTime, l1boa, leapsecFile, signalFileName, &
   configFileName, filedatabase, qtyTemplates, chunk, &
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
        filedatabase, chunk)
   !call dump(hGridStandard)

   fGridStandard = CreateFGrid(L_IntermediateFrequency, (/0.0_r8/))
   !call dump(fGridStandard)

   temperature = CreateQtyTemplate(l_temperature, filedatabase=filedatabase, &
                                   chunk=chunk, &
                                   avgrid=vGridStandard, ahgrid=hGridStandard)
   GPH = CreateQtyTemplate(l_gph, filedatabase=filedatabase, chunk=chunk, &
                           avgrid=vGridStandard, ahgrid=hGridStandard)
   O3 = CreateQtyTemplate(l_vmr, filedatabase=filedatabase, chunk=chunk, &
                          avgrid=vGridStandard, ahgrid=hGridStandard, qMolecule=l_o3)
   H2O = CreateQtyTemplate(l_vmr, filedatabase=filedatabase, chunk=chunk, &
                           avgrid=vGridStandard, ahgrid=hGridStandard, qMolecule=l_h2o, &
                           qLogBasis=.true., qMinValue=0.1E-6_r8)
   ptanGHz = CreateQtyTemplate(l_ptan, filedatabase=filedatabase, &
                               chunk=chunk, qInstModule=GHz)
   ! band 2,7,8 is the band whose radiances are to be computed
   ! see CFM document for a list of signals corresponding to bands
   band7 = CreateQtyTemplate(l_radiance, filedatabase=filedatabase, chunk=chunk, &
                             qSignal="R3:240.B7F:O3")
   band2 = CreateQtyTemplate(l_radiance, filedatabase=filedatabase, chunk=chunk, &
                             qSignal="R2:190.B2F:H2O")
   band8 = CreateQtyTemplate(l_radiance, filedatabase=filedatabase, chunk=chunk, &
                             qSignal="R3:240.B8F:PT")
   geodAltitude = CreateQtyTemplate(l_tngtgeodalt, filedatabase=filedatabase, &
                                    chunk=chunk, qInstModule=GHz)
   geocAlt = CreateQtyTemplate(l_tngtgeocalt, filedatabase=filedatabase, &
                               chunk=chunk, qInstModule=GHz)
   orbincl = CreateQtyTemplate(l_orbitInclination, filedatabase=filedatabase, &
                               chunk=chunk, qInstModule=sc)
   refGPH = CreateQtyTemplate(l_refGPH, avgrid=vGridRefGPH, ahgrid=hGridStandard)
   phitanGHz = CreateQtyTemplate(l_phitan, qInstModule="GHz", filedatabase=filedatabase, &
                                 chunk=chunk)

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
   call FillVectorQuantityFromL1B(quantity, chunk, filedatabase, &
      .false.)

   ! Fill orbit inclination, tangent geocentric altitude with
   ! data from MLS L1B file, and use them, along with other
   ! quantities to calculate ptan
   orbincl_vv => GetVectorQtyByTemplateIndex (state, orbincl_index)
   call FillVectorQuantityFromL1B(orbincl_vv, chunk, filedatabase, &
      .false.)

   geocAlt_vv => GetVectorQtyByTemplateIndex (state, geocAlt_index)
   call FillVectorQuantityFromL1B(geocAlt_vv, chunk, filedatabase, &
      .false.)

   ptanG_vv => GetVectorQtyByTemplateIndex (state, ptanGHz_index)

   quantity => GetVectorQtyByTemplateIndex (state, phitanGHz_index)
   call FillPhitanQuantity(quantity)

   ! calculate ptan
   !call dump(temperature_vv, details=3)
   call FillPtanQuantity (ptanG_vv, temperature_vv, refGPH_vv, &
                          h2o_vv, orbincl_vv, quantity, geocAlt_vv)

   ! GPH is filled by the forward model

   ! Create jacobian
   jacobian = CreatePlainMatrix(measurement, state)

   ! Call the forward model
   call ForwardModel (chunk, forwardModelConfigDatabase, state, &
                      stateExtra, measurement, jacobian)

   !call dump(measurement, details=3)

   ! Create an identical vector as simulated radiance vector for observed radiances
   observed = CreateVector(measurementTemplate, qtyTemplates)
   obsPrecision = CreateVector(measurementTemplate, qtyTemplates)

   ! Open l1brad
   error = InitializeMLSFile(l1bfile, content='l1brad', &
   name=trim(l1brad), shortName='L1BRAD', type=l_hdf, access=DFACC_RDONLY)
   if (error /= 0) &
      call MLSMessage (MLSMSG_Error, moduleName, &
      "Error initializing " // trim(l1brad))

   call mls_openFile(l1bfile, error)
   if (error /= 0 ) &
      call MLSMessage (MLSMSG_Error, moduleName, &
      "Error opening " // trim(l1brad))

   ! Add it to the filedatabase
   error = AddFileToDatabase(filedatabase, l1bfile)

   ! Fill band 2,7,8
   quantity => GetVectorQtyByTemplateIndex(observed, band2_index)
   precQty => GetVectorQtyByTemplateIndex(obsPrecision, band2_index)

   ! However, only band 9 and 25 have BOMask=1
   ! Because these bands have bright object status read from L1BOA file,
   ! we always have to have L1BOA file in the filedatabase.
   ! You have to fill the precision quantity first
   call FillVectorQuantityFromL1B(precQty, chunk, filedatabase, &
   .true.)
   call FillVectorQuantityFromL1B(quantity, chunk, filedatabase, &
   .false., precisionQuantity=precQty)

   quantity => GetVectorQtyByTemplateIndex(observed, band7_index)
   precQty => GetVectorQtyByTemplateIndex(obsPrecision, band7_index)
   call FillVectorQuantityFromL1B(precQty, chunk, filedatabase, &
   .true.)
   call FillVectorQuantityFromL1B(quantity, chunk, filedatabase, &
   .false., precisionQuantity=precQty)

   quantity => GetVectorQtyByTemplateIndex(observed, band8_index)
   precQty => GetVectorQtyByTemplateIndex(obsPrecision, band8_index)
   call FillVectorQuantityFromL1B(precQty, chunk, filedatabase, &
   .true.)
   call FillVectorQuantityFromL1B(quantity, chunk, filedatabase, &
   .false., precisionQuantity=precQty)

   ! For applying baseline corrections
   baseline2 = CreateQtyTemplate(l_L1BMAFBaseline, filedatabase=filedatabase, &
                                 chunk=chunk, qSignal="R2:190.B2F:H2O")
   baseline7 = CreateQtyTemplate(l_L1BMAFBaseline, filedatabase=filedatabase, &
                                 chunk=chunk, qSignal="R3:240.B7F:O3")
   baseline8 = CreateQtyTemplate(l_L1BMAFBaseline, filedatabase=filedatabase, &
                                 chunk=chunk, qSignal="R3:240.B8F:PT")
   baseline2_index = AddQuantityTemplateToDatabase(qtyTemplates, baseline2)
   baseline7_index = AddQuantityTemplateToDatabase(qtyTemplates, baseline7)
   baseline8_index = AddQuantityTemplateToDatabase(qtyTemplates, baseline8)
   baselineSelected = (/baseline2_index, baseline7_index, baseline8_index/)
   correctionTemplate = CreateVectorTemplate(qtyTemplates, baselineSelected)
   corrections = CreateVector(correctionTemplate, qtyTemplates)
   correctionNoise = CreateVector(correctionTemplate, qtyTemplates)

   band2L1BMAFBaseline => GetVectorQtyByTemplateIndex(corrections, baseline2_index)
   band7L1BMAFBaseline => GetVectorQtyByTemplateIndex(corrections, baseline7_index)
   band8L1BMAFBaseline => GetVectorQtyByTemplateIndex(corrections, baseline8_index)
   call FillVectorQuantityFromL1B(band2L1BMAFBaseline, chunk, filedatabase, &
   .false., suffix=' Baseline')
   call FillVectorQuantityFromL1B(band7L1BMAFBaseline, chunk, filedatabase, &
   .false., suffix=' Baseline')
   call FillVectorQuantityFromL1B(band8L1BMAFBaseline, chunk, filedatabase, &
   .false., suffix=' Baseline')

   quantity => GetVectorQtyByTemplateIndex(observed, band2_index)
   call ApplyBaseline(quantity, band2L1BMAFBaseline, .false., .false.)
   quantity => GetVectorQtyByTemplateIndex(observed, band7_index)
   call ApplyBaseline(quantity, band7L1BMAFBaseline, .false., .false.)
   quantity => GetVectorQtyByTemplateIndex(observed, band8_index)
   call ApplyBaseline(quantity, band8L1BMAFBaseline, .false., .false.)

   band2L1BMAFBaseline => GetVectorQtyByTemplateIndex(correctionNoise, baseline2_index)
   band7L1BMAFBaseline => GetVectorQtyByTemplateIndex(correctionNoise, baseline7_index)
   band8L1BMAFBaseline => GetVectorQtyByTemplateIndex(correctionNoise, baseline8_index)
   call FillVectorQuantityFromL1B(band2L1BMAFBaseline, chunk, filedatabase, &
   .false., suffix=' Baseline precision')
   call FillVectorQuantityFromL1B(band7L1BMAFBaseline, chunk, filedatabase, &
   .false., suffix=' Baseline precision')
   call FillVectorQuantityFromL1B(band8L1BMAFBaseline, chunk, filedatabase, &
   .false., suffix=' Baseline precision')

   ! quadrature is true because this is precision
   quantity => GetVectorQtyByTemplateIndex(obsPrecision, band2_index)
   call ApplyBaseline(quantity, band2L1BMAFBaseline, .true., .false.)
   quantity => GetVectorQtyByTemplateIndex(obsPrecision, band7_index)
   call ApplyBaseline(quantity, band7L1BMAFBaseline, .true., .false.)
   quantity => GetVectorQtyByTemplateIndex(obsPrecision, band8_index)
   call ApplyBaseline(quantity, band8L1BMAFBaseline, .true., .false.)

   ! CFM_MLSCleanup will close all files in the filedatabase

   diffVector = observed - measurement

   call dump(diffVector, details=3)

   ! Clean up allocated memory for creating observed radiance vector
   call DestroyVectorInfo(observed)
   call DestroyVectorInfo(obsPrecision)
   ! the subtraction create internal field for diffVector
   ! and we need to clean that up
   call DestroyVectorInfo(diffVector)

   ! Clean up For baseline
   call DestroyVectorInfo(corrections)
   call DestroyVectorInfo(correctionNoise)
   call DestroyVectorTemplateInfo(correctionTemplate)

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

! $Log$
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
!
