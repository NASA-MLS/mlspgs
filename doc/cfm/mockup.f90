program mockup

   use CFM 

   implicit none

   integer :: error
   type(ForwardModelConfig_T) :: fmConfig
   type(MLSFile_T), dimension(:), pointer :: filedatabase => NULL()
   type(MLSChunk_T) :: fakeChunk
   ! type(HGrid_T) ... more attributes for HGrid_T object(s)
   ! type(VGrid_T) ... more attributes for VGrid_T object(s)
   type(QuantityTemplate_T), dimension(:), pointer :: qtyTemplateDB => NULL()
   logical, dimension(:), pointer :: stateSelected => NULL()
   logical, dimension(:), pointer :: stateExtraSelected => NULL()
   type(VectorTemplate_T) :: stateVectorTemplate
   type(VectorTemplate_T) :: stateVectorExtraTemplate
   type(Vector_T) :: stateVector
   type(Vector_T) :: stateVectorExtra
   type(Vector_T) :: radiances
   type(ForwardModelStatus_T) :: fmStatus
   type(Matrix_T) :: jacobianMatrix

   !Executables

   ! Create fmConfig
   call CFM_MLSSetup(error, fmConfig, filedatabase, fakeChunk)
   if (error /= 0) then
      ! Print your choice of error message
      stop
   end if

   ! You need to create vgrid and hgrid for a qtyTemplate

   ! To create qtyTemplateDB, see the example for creating a quantity
   ! template snippet
   ! Initialize stateSelected and stateExtraSelected as appropriate

   ! Create vector template for stateVectorIn
   stateVectorInTemplate = CreateVectorTemplate (qtyTemplateDB, stateSelected)
   ! Create stateVectorIn
   stateVectorIn = CreateVector(stateVectorTemplate, qtyTemplateDB)
   ! Fill stateVectorIn%quantities, see VectorValue_T

   ! Create vector template for stateVectorExtra
   stateVectorInTemplate = CreateVectorTemplate (qtyTemplateDB, &
           stateExtraSelected)
   ! Create stateVectorExtra
   stateVectorExtra = CreateVector(stateVectorExtraTemplate, &
                                   qtyTemplateDB)
   ! Fill stateVectorExtra%quantities, see VectorValue_T
 
   ! Invoke the forward model  
   call ForwardModel (fmConfig, stateVectorIn, stateVectorExtra, radiances, &
                     & fmStatus, jacobian=jacobianMatrix)

   call CFM_MLSCleanup

end program
