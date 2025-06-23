! Copyright 2010, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.
module CFM_VectorTemplate_m
   use VectorsModule, only: VectorTemplate_T, NullifyVectorTemplate
   use QuantityTemplates, only: QuantityTemplate_T
   use Allocate_Deallocate, only: Allocate_test

   implicit none

   private

   public :: CreateVectorTemplate

!---------------------------- RCS Ident Info -------------------------------
   character(len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
   private :: not_used_here
!---------------------------------------------------------------------------

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
! Revision 1.5  2010/06/29 15:53:45  honghanh
! Add copyright comments and support for CVS log in the file
!
