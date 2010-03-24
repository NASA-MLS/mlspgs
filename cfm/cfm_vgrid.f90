module CFM_VGrid
   use VGridsDatabase, only: VGrid_T, NullifyVGrid, &
                             DestroyVGridContents, Dump
   use MLSCommon, only: r8
   use MLSStrings, only: ReadIntsFromChars
   use MLSStringLists, only: List2Array, NumStringElements
   use Intrinsic, only: l_zeta, &
                        phyq_dimensionless, phyq_pressure
   use Init_tables_module, only: l_logarithmic
   use Allocate_Deallocate, only: allocate_test, deallocate_test
   use MLSMessageModule, only: MLSMessage, MLSMSG_Error

   implicit none
   public :: CreateVGrid, DestroyVGridContents, Dump
   public :: VGrid_T

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

      call nullifyVGrid(vGrid) ! for sun's still useless compiler
      vgrid%name = 0
      vGrid%verticalCoordinate = coordinate

      if (type == l_logarithmic) then
         call LogarithmicFormula(vgrid, -1.0d0, formula, start)
      else
         call MLSMessage (MLSMSG_Error, moduleId, "VGrid's type not supported")
      end if

      if (coordinate == l_zeta .and. &
         (unit == phyq_pressure .or. unit == phyq_dimensionless)) then
         vgrid%surfs = -Log10(vgrid%surfs)
      end if
   end function

   ! not working right now
   subroutine LinearFormula (vgrid, stepSign, formula, start)
      real(r8), intent(in) :: stepSign
      type(VGrid_T) :: vgrid
      character(len=*), intent(in) :: formula
      real(r8), intent(in) :: start

      real(r8) :: step
      character(len=10) :: pair(2)
      character(len=20), dimension(:), pointer :: tokens => NULL()
      integer :: temp, numTokens
      integer :: i,j,k,n

      vgrid%noSurfs = 0
      numTokens = NumStringElements(formula, .false., ',')
      call allocate_test(tokens, numTokens, "tokens", moduleId)
      call List2Array(formula, tokens, .false., ',', .true.)

      ! get the number of surfaces
      do i = 1, size(tokens)
         call List2Array (tokens(i), pair, .false., ':', .true.)
         call ReadIntsFromChars(pair(1), temp)
         vgrid%nosurfs = vgrid%nosurfs + temp
      end do
      call allocate_test(vgrid%surfs, vgrid%nosurfs, 1, "vgrid%surfs", &
            moduleId)

      n = 0
      do i = 1, numTokens
         call List2Array (tokens(i), pair, .false., ':', .true.)
         call ReadIntsFromChars(pair(2), temp)
         n = n + 1
         vgrid%surfs(n,1) = temp

      end do
   end subroutine

   subroutine LogarithmicFormula (vgrid, stepSign, formula, start)
      real(r8), intent(in) :: stepSign
      type(VGrid_T) :: vgrid
      character(len=*), intent(in) :: formula
      real(r8), intent(in) :: start

      real(r8) :: step
      character(len=10) :: pair(2)
      character(len=20), dimension(:), pointer :: tokens => NULL()
      integer :: temp, numTokens
      integer :: i,j,k,n

      vgrid%noSurfs = 0
      numTokens = NumStringElements(formula, .false., ',')
      call allocate_test(tokens, numTokens, "tokens", moduleId)
      call List2Array(formula, tokens, .false., ',', .true.)

      ! get the number of surfaces
      do i = 1, size(tokens)
         call List2Array (tokens(i), pair, .false., ':', .true.)
         call ReadIntsFromChars(pair(1), temp)
         vgrid%nosurfs = vgrid%nosurfs + temp
      end do
      call allocate_test(vgrid%surfs, vgrid%nosurfs, 1, "vgrid%surfs", &
            moduleId)

      k = 1
      n = 1 ! One less surface the first time, since we have one at the start.
      vgrid%surfs(1,1) = start
      do i = 1, size(tokens)
         call List2Array (tokens(i), pair, .false., ':', .true.)
         call ReadIntsFromChars(pair(2), temp)
         step = 10.0 ** (stepSign / temp)
         call ReadIntsFromChars(pair(1), temp)
         do j = 1, temp - n
            k = k+1
            vgrid%surfs(k,1) = vgrid%surfs(k-1,1) * step
         end do
         n = 0
      end do

      call deallocate_test(tokens, "tokens", moduleId)

   end subroutine

end module
