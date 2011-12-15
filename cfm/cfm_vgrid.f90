! Copyright 2010, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.
module CFM_VGrid_m
   use VGridsDatabase, only: VGrid_T, NullifyVGrid
   use MLSCommon, only: r8
   use MLSStrings, only: ReadIntsFromChars
   use MLSStringLists, only: List2Array, NumStringElements
   use Intrinsic, only: l_zeta, &
                        phyq_dimensionless, phyq_pressure
   use Init_tables_module, only: l_logarithmic, l_explicit
   use Allocate_Deallocate, only: allocate_test, deallocate_test
   use MLSMessageModule, only: MLSMessage, MLSMSG_Error

   implicit none
   public :: CreateVGrid

   private

!---------------------------- RCS Ident Info -------------------------------
   character(len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
   private :: not_used_here
!---------------------------------------------------------------------------

   contains

   ! Create a VGrid, from a given formula, a starting value, a type
   ! (currently only supported logarithmic type)
   type(VGrid_T) function CreateVGrid (coordinate, unit, type, start, &
                 formula, values) result (vGrid)
      ! Currently only support l_zeta
      integer, intent(in) :: coordinate
      ! Must be either phyq_dimensionless or phyq_pressure
      integer, intent(in) :: unit
      ! Currently only support l_logarithmic (log-based vgrid)
      integer, intent(in) :: type

      ! If type is l_logarithmic, then start and formula
      ! should be provided
      ! The first value of surfaces
      real(r8), intent(in), optional :: start
      ! The format of formula is
      ! "number_of_surfaces:num_decades_between_surfaces".
      ! For example, "37:6" or "25:8,12:6".
      character(len=*), intent(in), optional :: formula
      ! If type is l_explicit, this must be present
      real(r8), dimension(:), intent(in), optional :: values

      call nullifyVGrid(vGrid) ! for sun's still useless compiler
      vgrid%name = 0
      vGrid%verticalCoordinate = coordinate

      if (type == l_logarithmic) then
         call LogarithmicFormula(vgrid, -1.0d0, formula, start)
      else if (type == l_explicit) then
         if (.not. present(values)) &
            call MLSMessage(MLSMSG_Error, ModuleName, &
            "Need values to create explicit VGrid")
         call CreateExplicitVGrid(vgrid, values)
      else
         call MLSMessage (MLSMSG_Error, ModuleName, "VGrid's type not supported")
      end if

      if (coordinate == l_zeta .and. &
         (unit == phyq_pressure .or. unit == phyq_dimensionless)) then
         vgrid%surfs = -Log10(vgrid%surfs)
      end if
   end function

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
      call allocate_test(tokens, numTokens, "tokens", ModuleName)
      call List2Array(formula, tokens, .false., ',', .true.)

      ! get the number of surfaces
      do i = 1, size(tokens)
         call List2Array (tokens(i), pair, .false., ':', .true.)
         call ReadIntsFromChars(pair(1), temp)
         vgrid%nosurfs = vgrid%nosurfs + temp
      end do
      call allocate_test(vgrid%surfs, vgrid%nosurfs, 1, "vgrid%surfs", &
            ModuleName)

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

      call deallocate_test(tokens, "tokens", ModuleName)

   end subroutine

   ! Fill the given vgrid with the given values
   subroutine CreateExplicitVGrid (vgrid, values)
      ! The vgrid to be filled
      type(VGrid_T) :: vgrid
      ! The values to fill
      real(r8), dimension(:), intent(in) :: values

      integer :: i

      vgrid%noSurfs = size(values)
      call allocate_test(vgrid%surfs, vgrid%nosurfs, 1, "vgrid%surfs", &
            ModuleName)

      do i = 1, vgrid%nosurfs
         vgrid%surfs(i,1) = values(i)
      end do
      vGrid%surfs = 10.0**( nint ( log10(vGrid%surfs) * values(1)) / values(1))

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
! Revision 1.10  2010/06/29 16:40:23  honghanh
! Remove all function/subroutine and user type forwarding from
! all CFM modules except for from cfm.f90
!
! Revision 1.9  2010/06/29 15:53:45  honghanh
! Add copyright comments and support for CVS log in the file
!
