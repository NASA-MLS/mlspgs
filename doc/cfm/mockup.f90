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
   logical, dimension(:), pointer :: radianceSelected => NULL()
   type(VectorTemplate_T) :: stateVectorTemplate
   type(VectorTemplate_T) :: stateVectorExtraTemplate
   type(VectorTemplate_T) :: radianceTemplate
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
   ! Initialize stateSelected, stateExtraSelected, radianceSelected 
   ! as appropriate

   ! Create vector template for stateVectorIn
   stateVectorTemplate = CreateVectorTemplate (qtyTemplateDB, stateSelected)
   ! Create stateVector
   stateVector = CreateVector(stateVectorTemplate, qtyTemplateDB)
   ! Fill stateVector%quantities, see VectorValue_T

   ! Create vector template for stateVectorExtra
   stateVectorExtraTemplate = CreateVectorTemplate (qtyTemplateDB, &
           stateExtraSelected)
   ! Create stateVectorExtra
   stateVectorExtra = CreateVector(stateVectorExtraTemplate, &
                                   qtyTemplateDB)
   ! Fill stateVectorExtra%quantities, see VectorValue_T

   ! You have to create an empty radiances vector
   radianceTemplate = CreateVectorTemplate(qtyTemplateDB, radianceSelected)
   radiances = CreateVector(radianceTemplate, qtyTemplateDB)
 
   ! Invoke the forward model  
   call ForwardModel (fmConfig, stateVector, stateVectorExtra, radiances, &
                     & fmStatus, jacobian=jacobianMatrix)

   call CFM_MLSCleanup

end program
