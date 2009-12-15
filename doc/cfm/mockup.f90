program mockup

   use CFM_MLSSetup_m, only: CFM_MLSSetup, CFM_MLSCleanup
   use ForwardModelConfig, only: ForwardModelConfig_T
   use ForwardModelWrappers, only: ForwardModel
   use VectorsModule, only: Vector_T, VectorTemplate_T, CreateVector
   use ForwardModelIntermediate, only: FORWARDMODELSTATUS_T
   use MatrixModule_1, only: Matrix_T
   use QuantityTemplates, only: QuantityTemplate_T, &
          AddQuantityTemplateToDatabase
   use FGrid, only: FGrid_T
   use HGridsDatabase, only: HGrid_T
   use MLSCommon, only: MLSFile_T

   use Dummy, only: dummyHGrid => dummyHGrid1

   implicit none

   integer :: error 
   type(ForwardModelConfig_T) :: aForwardModelConfig
   type (ForwardModelConfig_T), pointer :: ForwardModelConfigDatabase(:) => NULL()
   type(Vector_T) :: fwdModelIn, fwdModelExtra, fwdModelOut
   type(VectorTemplate_T) :: inputVectorTemplate
   type(QuantityTemplate_T), dimension(:), pointer :: qtyTemplates => NULL()
   type(ForwardModelStatus_t) :: fwdStatus
   type(Vector_T), dimension(:), pointer :: vectors => NULL()
   type(Matrix_T) :: aJacobian
   type(QuantityTemplate_T) :: qtyTemplate
   type(FGrid_T), dimension(:), pointer :: fgrids => NULL()
   type(HGrid_T), dimension(:), pointer :: hgrids => NULL()
   type(MLSFile_T), dimension(:), pointer :: filedatabase => NULL()
   real :: a

   !Executables
   a = .0
   print *, "Haley: a", a
   call CFM_MLSSetup(error, ForwardModelConfigDatabase)

   if (error /=0) stop

   aForwardModelConfig = ForwardModelConfigDatabase(1)

   call CreateQuantityTemplate(fgrids, hgrids, filedatabase, qtyTemplate)
   allocate (qtyTemplates(1), stat=error)
   if (error /= 0) stop
   error = AddQuantityTemplateToDatabase (qtyTemplates, qtyTemplate)

   call VectorFill (inputVectorTemplate, qtyTemplates, fwdModelIn)

   call ForwardModel (aForwardModelConfig, fwdModelIn, fwdModelExtra, fwdModelOut, &
                     & fwdStatus, jacobian=aJacobian, vectors=vectors)

   call CFM_MLSCleanup

   contains 

      subroutine VectorFill (inputVectorTemplate, qtyTemplates, vector)
         type(VectorTemplate_T), intent(in) :: inputVectorTemplate
         type(QuantityTemplate_T), dimension(:), pointer :: qtyTemplates
         type(Vector_T), intent(out) :: vector

         vector = CreateVector(-1, inputVectorTemplate, qtyTemplates)

      end subroutine VectorFill

      subroutine CreateQuantityTemplate (fgrids, hgrids, filedatabase, &
                     qtyTemplate) 
         type(FGrid_T), dimension(:), pointer :: fgrids
         type(HGrid_T), dimension(:), pointer :: hgrids
         type(MLSFile_T), dimension(:), pointer :: filedatabase
         type(QuantityTemplate_T), intent(out) :: qtyTemplate

         qtyTemplate%name = -1

      end subroutine CreateQuantityTemplate

end program
