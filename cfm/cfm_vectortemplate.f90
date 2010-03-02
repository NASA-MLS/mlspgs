module CFM_VectorTemplate
   use VectorsModule, only: VectorTemplate_T, NullifyVectorTemplate
   use QuantityTemplates, only: QuantityTemplate_T
   use Allocate_Deallocate, only: Allocate_test

   implicit none

   private

   public :: CreateVectorTemplate

   character(len=20) :: moduleName = "CFM_VectorTemplate"

   contains

   type(VectorTemplate_T) function CreateVectorTemplate (qtyTemplateDB, selectedQty) &
         result (vecTemplate)

      integer, intent(in) :: selectedQty(:)
      type(QuantityTemplate_T), intent(in) :: qtyTemplateDB(:)

      call nullifyVectorTemplate(vecTemplate)
      vecTemplate%name = 0
      vecTemplate%noQuantities = size(selectedQty)
      vecTemplate%totalInstances = sum(qtyTemplateDB(selectedQty)%noInstances)
      vecTemplate%totalElements = &
          sum(qtyTemplateDB(selectedQty)%noInstances & 
          *qtyTemplateDB(selectedQty)%instanceLen)
      call allocate_test(vecTemplate%quantities, vecTemplate%noQuantities, &
            'vecTemplate%quantities', moduleName)
      vecTemplate%quantities = selectedQty
   end function
end module
