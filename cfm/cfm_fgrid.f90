! Copyright 2010, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.
module CFM_FGrid_m
   use FGrid, only: FGrid_T, nullifyFGrid
   use Allocate_Deallocate, only : allocate_test, Deallocate_test
   use Intrinsic, only: l_none
   use MLSCommon, only: r8

   implicit none

   public :: CreateFGrid

!---------------------------- RCS Ident Info -------------------------------
   character(len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
   private :: not_used_here
!---------------------------------------------------------------------------

   private
   contains

   ! Creates and fills an FGrid_T given a coordinate and a list of values
   type(FGrid_T) function CreateFGrid (frequencyCoordinate, values) &
   result (fGrid)
      integer, intent(in) :: frequencyCoordinate
      real(r8), dimension(:), intent(in) :: values

      integer :: i

      call nullifyFGrid(fGrid)
      fGrid%name = -1
      fGrid%frequencyCoordinate = frequencyCoordinate
      fGrid%noChans = size(values)
      call allocate_test(fGrid%values, fGrid%noChans, 'fgrid%values', &
           moduleName)
      do i = 1, fgrid%nochans
         fgrid%values(i) = values(i)
      end do

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
