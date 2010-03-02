module CFM_Vector
   use VectorsModule, only: Vector_T, VectorTemplate_T
   use QuantityTemplates, only: QuantityTemplate_T
   use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, &
         MLSMSG_Error, MLSMSG_Warning

   implicit none
   public :: CreateVector
   private
   character(len=20) :: moduleName="CFM_Vector"

   contains

   type(Vector_T) function CreateVector (vectorTemplate, qtyDatabase) &
             result (vector)
      type (VectorTemplate_T), intent (in), target :: vectorTemplate
      type (QuantityTemplate_T), dimension(:), intent(in), target :: qtyDatabase

      integer :: quantity
      integer :: status

      vector%template = vectorTemplate
      allocate (vector%quantities(vectorTemplate%noQuantities), stat=status)
      if (status /= 0) call MLSMessage (MLSMSG_Error, moduleName, &
          MLSMSG_Allocate // "Vector quantities")

      do quantity = 1, vectorTemplate%noQuantities
         vector%quantities(quantity)%index = quantity
         vector%quantities(quantity)%template = &
            qtyDatabase(vectorTemplate%quantities(quantity))
      end do
   end function

end module
