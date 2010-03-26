module CFM_VectorTemplate_m
   use VectorsModule, only: VectorTemplate_T, NullifyVectorTemplate, &
                            DestroyVectorTemplateInfo, Dump
   use QuantityTemplates, only: QuantityTemplate_T
   use Allocate_Deallocate, only: Allocate_test

   implicit none

   private

   public :: CreateVectorTemplate, DestroyVectorTemplateInfo
   ! Please see mockup.f90 for examples of using Dump subroutines.
   ! Details level typically range from -3 to 3, for most subroutines
   ! the domain for the details argument can be smaller. However,
   ! passing a too big or too small number won't cause the subroutine
   ! to crash. This applies to all Dump subroutines from all modules.
   public :: Dump
   public :: VectorTemplate_T

   character(len=20) :: moduleName = "CFM_VectorTemplate"

   contains

   ! Create a VectorTemplate_T, which references one or many quantities.
   ! To ensure correctness, the same qtyTemplateDB, must be passed to
   ! CreateVector function in conjunction with the returned VectorTemplate_T.
   type(VectorTemplate_T) function CreateVectorTemplate (qtyTemplateDB, selectedQty) &
         result (vecTemplate)
      ! a list of the array indices of the chosen quantity templates.
      integer, intent(in) :: selectedQty(:)
      ! a 1-based, 1-dimensional array of all QuantityTemplate_T objects
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
