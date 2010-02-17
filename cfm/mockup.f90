program mockup
   use CFM_MLSSetup_m, only: CFM_MLSSetup, CFM_MLSCleanup
   use CFM_VGrid, only: CreateVGrid
   use CFM_HGrid, only: CreateRegularHGrid
   use CFM_FGrid, only: CreateFGrid
   use CFM_QuantityTemplate, only: CreateQtyTemplate, InitQuantityTemplates

   use Chunks_m, only: MLSChunk_T
   use ForwardModelConfig, only: ForwardModelConfig_T
   use MLSCommon, only: MLSFile_T, r8
   use VGridsDatabase, only: VGrid_T
   use HGridsDatabase, only: HGrid_T
   use Intrinsic, only: phyq_pressure, l_zeta, L_IntermediateFrequency
   use Init_Tables_Module, only: l_logarithmic, l_temperature, l_gph, l_vmr
   use FGrid, only: FGrid_T
   use QuantityTemplates, only: QuantityTemplate_T, Dump

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
   type(VGrid_T) :: vGridStandard
   type(HGrid_T) :: hGridStandard
   type(FGrid_T) :: fGridStandard
   type(QuantityTemplate_T) :: temperature, GPH, H2O, O3
   type(QuantityTemplate_T) :: qtyTemplates(4)
   character(len=3) :: GHz = "GHz"

   call CFM_MLSSetup(error, ForwardModelConfigDatabase, filedatabase, fakeChunk)
   if (error /=0) stop

   vGridStandard = CreateVGrid(l_zeta, l_logarithmic, &
                               1000.0d0, "37:6", phyq_pressure)
   hGridStandard = CreateRegularHGrid(GHz, 0.0_r8, 1.5_r8, filedatabase, fakeChunk)
   fGridStandard = CreateFGrid(L_IntermediateFrequency, (/0.0_r8/))

   ! Have to initialize before we start creating quantity templates
   call InitQuantityTemplates

   temperature = CreateQtyTemplate(l_temperature, vGridStandard, hGridStandard)
   GPH = CreateQtyTemplate(l_gph, vGridStandard, hGridStandard)
   O3 = CreateQtyTemplate(l_vmr, vGridStandard, hGridStandard, qMolecule="O3")
   H2O = CreateQtyTemplate(l_vmr, vGridStandard, hGridStandard, qMolecule="H2O", &
                           qLogBasis=.true., qMinValue=0.1_r8)
   qtyTemplates(1) = temperature
   qtyTemplates(2) = GPH
   qtyTemplates(3) = O3
   qtyTemplates(4) = H2O

   do i = 1, size(qtyTemplates)
      call dump(qtyTemplates(i), details=2)
   end do

   call CFM_MLSCleanup
end program
