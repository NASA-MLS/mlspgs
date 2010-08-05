! Copyright 2010, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.
module cfm_matrix_m
   use MatrixModule_1, only: Matrix_T
   use VectorsModule, only: Vector_T

   implicit none

!---------------------------- RCS Ident Info -------------------------------
   character(len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
   private :: not_used_here
!---------------------------------------------------------------------------

   private
   public :: CreatePlainMatrix

   contains

   type(Matrix_T) function CreatePlainMatrix (rows, columns) result(matrix)
      use MatrixModule_1, only: CreateEmptyMatrix, NullifyMatrix

      type(Vector_T), intent(in) :: rows, columns

      ! Executables
      call NullifyMatrix(matrix)
      call CreateEmptyMatrix(matrix, 0, rows, columns)
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
