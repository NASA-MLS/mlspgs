program mockup
   use input   ! Provides hard-coded input for testing purposes only

   use CFM_MLSSetup_m, only: CFM_MLSSetup, CFM_MLSCleanup, MLSChunk_T
   use CFM_VGrid_m, only: CreateVGrid, DestroyVGridContents, &
                        VGrid_T, Dump
   use CFM_HGrid_m, only: CreateRegularHGrid, HGrid_T, &
                        DestroyHGridContents, Dump
   use CFM_FGrid_m, only: CreateFGrid, FGrid_T, DestroyFGridContents, &
                        Dump
   use CFM_QuantityTemplate_m, only: CreateQtyTemplate, Dump, &
                        AddQuantityTemplateToDatabase, &
                        DestroyQuantityTemplateDatabase, &
                        QuantityTemplate_T, InitQuantityTemplates, &
                        ConstructMIFGeolocation
   use CFM_VectorTemplate_m, only: CreateVectorTemplate, Dump, &
                        VectorTemplate_T, DestroyVectorTemplateInfo
   use CFM_Vector_m, only: CreateVector, Dump, &
                         Vector_T, VectorValue_T, &
                         DestroyVectorInfo, GetVectorQtyByTemplateIndex
   use CFM_Fill_m, only: ExplicitFillVectorQuantity
   use MLSCommon, only: MLSFile_T, r8

   use ForwardModelConfig, only: ForwardModelConfig_T
   use Init_tables_module, only: l_logarithmic, l_zeta, &
                                 phyq_pressure, L_IntermediateFrequency

   implicit none

!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

   integer :: error, i, numQty
   type(ForwardModelConfig_T), dimension(:), pointer :: forwardModelConfigDatabase
   type(MLSFile_T), dimension(:), pointer :: filedatabase
   type(MLSChunk_T) :: fakeChunk
   type(VGrid_T) :: vGridStandard
   type(HGrid_T) :: hGridStandard
   type(FGrid_T) :: fGridStandard
   type(QuantityTemplate_T) :: temperature, GPH, H2O, O3, ptanGHz, ptanTHz, band7, &
                               geodAltitude
   type(QuantityTemplate_T), dimension(:), pointer :: mifGeoLocation => NULL()
   type(QuantityTemplate_T), dimension(:), pointer :: qtyTemplates
   type(VectorTemplate_T) :: stateTemplate, measurementTemplate
   type(Vector_T) :: state, measurement
   character(len=3) :: GHz = "GHz"
   character(len=3) :: THz = "THz"
   integer :: stateSelected(7), measurementSelected(1)
   type(VectorValue_T) :: quantity

   ! Reads L2CF from standard input, populates ForwardModelConfigDatabase, filedatabase,
   ! and fakeChunk
   call CFM_MLSSetup(startTime, endTime, l1boa, error, filedatabase, fakeChunk, &
                     ForwardModelConfigDatabase)
   if (error /= 0 ) stop

   vGridStandard = CreateVGrid (l_zeta, l_logarithmic, 1000.0d0, "37:6", phyq_pressure)
   !call dump(vGridStandard, details=2)

   ! Have insetoverlaps, and not single
   hGridStandard = CreateRegularHGrid(GHz, 0.0_r8, 1.5_r8, .true., &
        filedatabase, fakeChunk)
   !call dump(hGridStandard)

   fGridStandard = CreateFGrid(L_IntermediateFrequency, (/0.0_r8/))
   !call dump(fGridStandard)

   ! Have to initialize before we start creating quantity templates
   call InitQuantityTemplates
   ! Reading mifGeolocation from L1BOA file
   call ConstructMIFGeolocation(mifGeoLocation, filedatabase, fakeChunk)
   temperature = CreateQtyTemplate("temperature", filedatabase=filedatabase, &
                                   avgrid=vGridStandard, ahgrid=hGridStandard, &
                                   mifGeolocation=mifGeolocation)
   GPH = CreateQtyTemplate("GPH", filedatabase=filedatabase, avgrid=vGridStandard, &
                           ahgrid=hGridStandard, mifGeolocation=mifGeolocation)
   O3 = CreateQtyTemplate("vmr", filedatabase=filedatabase, avgrid=vGridStandard, &
                          ahgrid=hGridStandard, qMolecule="O3", mifGeolocation=mifGeolocation)
   H2O = CreateQtyTemplate("vmr", filedatabase=filedatabase, avgrid=vGridStandard, &
                           ahgrid=hGridStandard, qMolecule="H2O", &
                           qLogBasis=.true., qMinValue=0.1_r8, mifGeolocation=mifGeolocation)
   ptanGHz = CreateQtyTemplate("ptan", filedatabase=filedatabase, qInstModule=GHz, &
                               mifGeoLocation=mifGeolocation)
   ptanTHz = CreateQtyTemplate("ptan", filedatabase=filedatabase, qInstModule=THz, &
                               mifGeolocation=mifGeolocation)
   band7 = CreateQtyTemplate("radiance", filedatabase=filedatabase, chunk=fakeChunk, &
                             qSignal="R3:240.B7F:O3", mifGeolocation=mifGeolocation)
   geodAltitude = CreateQtyTemplate("geodAltitude", filedatabase=filedatabase, &
                                    qInstModule=GHz, mifGeolocation=mifGeolocation)

   numQty = AddQuantityTemplateToDatabase(qtyTemplates, temperature)
   numQty = AddQuantityTemplateToDatabase(qtyTemplates, GPH)
   numQty = AddQuantityTemplateToDatabase(qtyTemplates, O3)
   numQty = AddQuantityTemplateToDatabase(qtyTemplates, H2O)
   numQty = AddQuantityTemplateToDatabase(qtyTemplates, ptanGHz)
   numQty = AddQuantityTemplateToDatabase(qtyTemplates, ptanTHz)
   numQty = AddQuantityTemplateToDatabase(qtyTemplates, band7)
   numQty = AddQuantityTemplateToDatabase(qtyTemplates, geodAltitude)

!   call dump(ptanGHz, details=2)
!   call dump(ptanTHz, details=2)
!   call dump(temperature, details=2)
!   call dump(GPH, details=2)
!   call dump(O3, details=2)
!   call dump(H2O, details=2)
!   call dump(band7, details=2)
!   call dump(geodAltitude, details=2)

   stateSelected = (/5,6,1,2,4,3,8/)   ! The numbers are the order that
                                       ! quantities template were added
   stateTemplate = CreateVectorTemplate(qtyTemplates, stateSelected)
!   call dump(stateTemplate, details=2, quantities=qtyTemplates)

   measurementSelected = (/7/)
   measurementTemplate = CreateVectorTemplate(qtyTemplates, measurementSelected)
!   call dump(measurementTemplate, details=2, quantities=qtyTemplates)

   state = CreateVector(stateTemplate, qtyTemplates)
   measurement = CreateVector(measurementTemplate, qtyTemplates)

!   call dump(state, details=2)
!   call dump(measurement, details=2)

   quantity = GetVectorQtyByTemplateIndex(state, 1)
   call ExplicitFillVectorQuantity(quantity, TemperatureInput)

   quantity = GetVectorQtyByTemplateIndex(state, 2)
   call ExplicitFillVectorQuantity(quantity, GPHInput)

   quantity = GetVectorQtyByTemplateIndex(state, 4)
   call ExplicitFillVectorQuantity(quantity, H2OInput)

   quantity = GetVectorQtyByTemplateIndex(state, 3)
   call ExplicitFillVectorQuantity(quantity, O3Input)

   call dump(state, details=2)

   ! Temporarily take out the call to ForwardModel because this step is excluded
   ! for v0.1 delivery
   ! call ForwardModel (fmConfig, stateVector, stateVectorExtra, radiances, &
   ! & fmStatus, jacobian=jacobianMatrix)

   ! Clean up memory
   call DestroyVectorInfo (state)
   call DestroyVectorInfo (measurement)
   call DestroyVectorTemplateInfo(stateTemplate)
   call DestroyVectorTemplateInfo(measurementTemplate)
   call DestroyQuantityTemplateDatabase (qtyTemplates)
   call DestroyHGridContents(hGridStandard)
   call DestroyVGridContents(vGridStandard)
   call DestroyFGridContents(fGridStandard)
   call CFM_MLSCleanup
end program
