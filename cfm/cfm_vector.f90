! Copyright 2010, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.
module CFM_Vector_m
   use VectorsModule, only: Vector_T, VectorTemplate_T, VectorValue_T
   use QuantityTemplates, only: QuantityTemplate_T
   use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, &
         MLSMSG_Error, MLSMSG_Warning
   use Allocate_Deallocate, only: Allocate_Test
   use MLSCommon, only: r8
   use String_table, only: create_string

   implicit none

   public :: CreateVector

!---------------------------- RCS Ident Info -------------------------------
   character(len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
   private :: not_used_here
!---------------------------------------------------------------------------

   private

   contains

   type(Vector_T) function CreateVector (vectorTemplate, qtyDatabase, name) &
             result (vector)
      ! template listing the quantities to be stored in this vector
      type (VectorTemplate_T), intent (in), target :: vectorTemplate
      ! quantity template database to retrieve the template for the quantities
      ! to be stored in this vector
      type (QuantityTemplate_T), dimension(:), intent(in), target :: qtyDatabase
      ! name of the vector as string
      character(len=*), optional :: name

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

      if (present(name)) then
         vector%name = create_string(name)
      else
         vector%name = 0
      end if
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

!--------------------------- end bloc --------------------------------------
   logical function not_used_here()
   character (len=*), parameter :: IdParm = &
       "$Id$"
   character (len=len(idParm)) :: Id = idParm
      not_used_here = (id(1:1) == ModuleName(1:1))
      print *, Id ! .mod files sometimes change if PRINT is added
   end function not_used_here
!---------------------------------------------------------------------------

end module

! $Log$
! Revision 1.6  2010/06/29 16:40:23  honghanh
! Remove all function/subroutine and user type forwarding from
! all CFM modules except for from cfm.f90
!
! Revision 1.5  2010/06/29 15:53:45  honghanh
! Add copyright comments and support for CVS log in the file
!
