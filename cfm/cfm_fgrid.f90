module CFM_FGrid
   use FGrid, only: FGrid_T, nullifyFGrid
   use Allocate_Deallocate, only : allocate_test, Deallocate_test
   use Intrinsic, only: l_none
   use MLSCommon, only: r8

   implicit none
   public :: CreateFGrid, DestroyFGridContents
   private
   character(len=20), parameter :: moduleName="CFM_FGrid"
   contains

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

   subroutine DestroyFGridContents (fGrid)
      type(FGrid_T) :: fGrid

      call deallocate_test (fGrid%values, 'fgrid%values', modulename)
      fgrid%name = 0
      fgrid%noChans = 0
      fGrid%frequencyCoordinate = l_none
   end subroutine
end module
