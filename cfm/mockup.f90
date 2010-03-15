program mockup
   use input
   use CFM_MLSSetup_m, only: CFM_MLSSetup, CFM_MLSCleanup
   use CFM_VGrid, only: CreateVGrid
   use CFM_HGrid, only: CreateRegularHGrid
   use CFM_FGrid, only: CreateFGrid
   use CFM_QuantityTemplate, only: CreateQtyTemplate
   use CFM_VectorTemplate, only: CreateVectorTemplate
   use CFM_Vector, only: CreateVector
   use CFM_Fill, only: ExplicitFillVectorQuantity

   use Chunks_m, only: MLSChunk_T
   use ForwardModelConfig, only: ForwardModelConfig_T
   use MLSCommon, only: MLSFile_T, r8
   use VGridsDatabase, only: VGrid_T, DestroyVGridContents
   use HGridsDatabase, only: HGrid_T, DestroyHGridContents, Dump
   use Intrinsic, only: L_IntermediateFrequency
   use FGrid, only: FGrid_T, DestroyFGridContents
   use QuantityTemplates, only: QuantityTemplate_T, &
         AddQuantityTemplateToDatabase, DestroyQuantityTemplateDatabase
   use Init_tables_module, only: l_logarithmic, l_zeta, &
                                 phyq_pressure
   use VectorsModule, only: VectorTemplate_T, Vector_T, VectorValue_T, &
                            DestroyVectorTemplateInfo, DestroyVectorInfo, &
                            GetVectorQtyByTemplateIndex
   use Construct, only: ConstructMIFGeolocation
   use ConstructQuantityTemplates, only: InitQuantityTemplates

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

   call CFM_MLSSetup(error, ForwardModelConfigDatabase, filedatabase, fakeChunk)
   if (error /=0) stop

   vGridStandard = CreateVGrid(l_zeta, l_logarithmic, &
                               1000.0d0, "37:6", phyq_pressure)
   hGridStandard = CreateRegularHGrid(GHz, 0.0_r8, 1.5_r8, filedatabase, fakeChunk)
   call dump(hGridStandard)
   fGridStandard = CreateFGrid(L_IntermediateFrequency, (/0.0_r8/))

   ! Have to initialize before we start creating quantity templates
   call ConstructMIFGeolocation(mifGeoLocation, filedatabase, fakeChunk)
   call InitQuantityTemplates
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

   stateSelected = (/1,2,3,4,5,6,8/)
   stateTemplate = CreateVectorTemplate(qtyTemplates, stateSelected)

   measurementSelected = (/7/)
   measurementTemplate = CreateVectorTemplate(qtyTemplates, measurementSelected)

   state = CreateVector(stateTemplate, qtyTemplates)
   measurement = CreateVector(measurementTemplate, qtyTemplates)

   quantity = GetVectorQtyByTemplateIndex(state, 1)
   call ExplicitFillVectorQuantity(quantity, TemperatureInput)   

   quantity = GetVectorQtyByTemplateIndex(state, 2)
   call ExplicitFillVectorQuantity(quantity, GPHInput)

   quantity = GetVectorQtyByTemplateIndex(state, 4)
   call ExplicitFillVectorQuantity(quantity, H2OInput)

   quantity = GetVectorQtyByTemplateIndex(state, 3)
   call ExplicitFillVectorQuantity(quantity, O3Input)

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
