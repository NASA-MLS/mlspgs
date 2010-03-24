module CFM_VectorTemplate_m
   use VectorsModule, only: VectorTemplate_T, NullifyVectorTemplate, &
                            DestroyVectorTemplateInfo, Dump
   use QuantityTemplates, only: QuantityTemplate_T
   use Allocate_Deallocate, only: Allocate_test

   implicit none

   private

   public :: CreateVectorTemplate, DestroyVectorTemplateInfo
   public :: Dump
   public :: VectorTemplate_T

   character(len=20) :: moduleName = "CFM_VectorTemplate"

   contains

   ! qtyTemplateDB is a 1-based, 1-dimensional array of QuantityTemplate_T
   ! selectedQty is a list of the array indices of the chosen quantity templates.
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
