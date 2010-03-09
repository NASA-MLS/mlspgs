module CFM_Vector
   use VectorsModule, only: Vector_T, VectorTemplate_T
   use QuantityTemplates, only: QuantityTemplate_T
   use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, &
         MLSMSG_Error, MLSMSG_Warning
   use Allocate_Deallocate, only: Allocate_Test
   use MLSCommon, only: r8

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
      call CreateValues(vector)
   end function

   ! =====     Private Procedures     =====================================
   subroutine CreateValues ( Vector)
      ! Allocate space for the values of a vector.
      type(Vector_T), intent(inout) :: Vector
      integer :: QTY
      real(r8), parameter :: MYHUGE = 1.0e15

      do qty = 1, size(vector%quantities)
         call allocate_test ( vector%quantities(qty)%values, &
           & vector%quantities(qty)%template%noChans * &
           & vector%quantities(qty)%template%noSurfs, &
           & vector%quantities(qty)%template%noInstances, &
           & "vector%quantities(qty)%%values", ModuleName )
         vector%quantities(qty)%values = 0.0_r8
    end do
  end subroutine


end module
