module CFM_VGrid
   use VGridsDatabase, only: VGrid_T, NullifyVGrid
   use MLSCommon, only: r8
   use Intrinsic, only: l_zeta, phyq_pressure
   use MLSStrings, only: ReadIntsFromChars
   use MLSStringLists, only: List2Array, NumStringElements
   use Init_tables_module, only: l_logarithmic
   use Allocate_Deallocate, only: allocate_test, deallocate_test

   implicit none
   public :: CreateVGrid

   private
   character(len=20) :: moduleId = "CFM_VGrid.f90"

   contains

   type(VGrid_T) function CreateVGrid (coordinate, type, start, &
                 formula, unit) result (vGrid)

      real(r8), intent(in) :: start
      character(len=*), intent(in) :: formula
      integer, intent(in) :: coordinate
      integer, intent(in) :: type
      integer, intent(in) :: unit

      character(len=10) :: pair(2)
      character(len=20), dimension(:), pointer :: tokens => NULL()
      integer :: temp, currentsurf, numTokens
      integer :: i,j
      real(r8) :: step

      call nullifyVGrid(vGrid) ! for sun's still useless compiler
      vgrid%name = 0
      vgrid%noSurfs = 0
      vGrid%verticalCoordinate = coordinate

      numTokens = NumStringElements(formula, .false., ',')
      call allocate_test(tokens, numTokens, "tokens", moduleId)
      call List2Array(formula, tokens, .false., ',', .true.)
      if (type == l_logarithmic) then
         ! get the number of surfaces
         do i = 1, size(tokens)
            call List2Array (tokens(i), pair, .false., ':', .true.)
            call ReadIntsFromChars(pair(1), temp)
            vgrid%nosurfs = vgrid%nosurfs + temp
         end do
         call allocate_test(vgrid%surfs, vgrid%nosurfs, 1, "vgrid%surfs", &
               moduleId)
         currentSurf = 1
         vgrid%surfs(1,1) = start
         do i = 1, size(tokens)
            call List2Array (tokens(i), pair, .false., ':', .true.)
            call ReadIntsFromChars(pair(2), temp)
            step = 10.0 ** 1.0d0 / temp
            call ReadIntsFromChars(pair(1), temp)
            do j = currentSurf, (currentSurf + temp - 1)
               vgrid%surfs(j,1) = vgrid%surfs(j-1,1) * step
            end do
            currentSurf = currentSurf + temp
         end do
      end if
      if (coordinate == l_zeta .and. unit == phyq_pressure) then
         vgrid%surfs = -Log10(vgrid%surfs)
      end if
      call deallocate_test(tokens, "tokens", moduleId)

   end function

end module
