program mockup
   use CFM_MLSSetup_m, only: CFM_MLSSetup, CFM_MLSCleanup
   use CFM_VGrid, only: CreateVGrid
   use CFM_HGrid, only: CreateRegularHGrid
   use CFM_FGrid, only: CreateFGrid
   use CFM_QuantityTemplate, only: CreateQtyTemplate, InitQuantityTemplates
   use CFM_VectorTemplate, only: CreateVectorTemplate

   use Chunks_m, only: MLSChunk_T
   use ForwardModelConfig, only: ForwardModelConfig_T
   use MLSCommon, only: MLSFile_T, r8
   use VGridsDatabase, only: VGrid_T, DestroyVGridContents
   use HGridsDatabase, only: HGrid_T, DestroyHGridContents
   use Intrinsic, only: phyq_pressure, l_zeta, L_IntermediateFrequency
   use FGrid, only: FGrid_T, DestroyFGridContents
   use QuantityTemplates, only: QuantityTemplate_T, &
         AddQuantityTemplateToDatabase, DestroyQuantityTemplateDatabase
   use Init_tables_module, only: l_logarithmic
   use VectorsModule, only: VectorTemplate_T, Vector_T, &
                            DestroyVectorTemplateInfo, Dump

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
   type(QuantityTemplate_T), dimension(:), pointer :: qtyTemplates
   type(VectorTemplate_T) :: stateTemplate, measurementTemplate
   type(Vector_T) :: state
   character(len=3) :: GHz = "GHz"
   character(len=3) :: THz = "THz"
   integer :: stateSelected(7), measurementSelected(1)

   call CFM_MLSSetup(error, ForwardModelConfigDatabase, filedatabase, fakeChunk)
   if (error /=0) stop

   vGridStandard = CreateVGrid(l_zeta, l_logarithmic, &
                               1000.0d0, "37:6", phyq_pressure)
   hGridStandard = CreateRegularHGrid(GHz, 0.0_r8, 1.5_r8, filedatabase, fakeChunk)
   fGridStandard = CreateFGrid(L_IntermediateFrequency, (/0.0_r8/))

   ! Have to initialize before we start creating quantity templates
   call InitQuantityTemplates
   temperature = CreateQtyTemplate("temperature", avgrid=vGridStandard, ahgrid=hGridStandard)
   GPH = CreateQtyTemplate("GPH", avgrid=vGridStandard, ahgrid=hGridStandard)
   O3 = CreateQtyTemplate("vmr", avgrid=vGridStandard, &
                          ahgrid=hGridStandard, qMolecule="O3")
   H2O = CreateQtyTemplate("vmr", avgrid=vGridStandard, &
                           ahgrid=hGridStandard, qMolecule="H2O", &
                           qLogBasis=.true., qMinValue=0.1_r8)
   ptanGHz = CreateQtyTemplate("ptan", qInstModule=GHz)
   ptanTHz = CreateQtyTemplate("ptan", qInstModule=THz)
   band7 = CreateQtyTemplate("radiance", filedatabase, fakeChunk, qSignal="R3:240.B7F:O3")
   geodAltitude = CreateQtyTemplate("geodAltitude", qInstModule=GHz)

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

   call dump(stateTemplate, quantities=qtyTemplates)
   call dump(measurementTemplate, quantities=qtyTemplates)

   call DestroyVectorTemplateInfo(stateTemplate)
   call DestroyVectorTemplateInfo(measurementTemplate)
   call DestroyQuantityTemplateDatabase (qtyTemplates)
   call DestroyHGridContents(hGridStandard)
   call DestroyVGridContents(vGridStandard)
   call DestroyFGridContents(fGridStandard)
   call CFM_MLSCleanup
end program
