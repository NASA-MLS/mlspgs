! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module FGrid                    ! Frequency grid information

  use Allocate_Deallocate, only: Allocate_test, Deallocate_test
  use EXPR_M, only: EXPR
  use Intrinsic, only: L_Frequency, L_IntermediateFrequency, &
    & L_LSBFrequency, L_None, L_USBFrequency
  use Init_tables_module, only: F_Coordinate, F_Values
  use MLSCommon, only: r8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Error, &
    & MLSMSG_Deallocate
  use Tree, only: DECORATION, NSONS, SUBTREE
  use Units, only: PHYQ_DIMENSIONLESS, PHYQ_FREQUENCY

  implicit none
  private

  public :: FGrid_T, AddFGridToDatabase, CreateFGridFromMLSCFInfo, &
    & DestroyFGridContents, DestroyFGridDatabase

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  !---------------------------------------------------------------------------

  ! This module contains code for defining and maintaining 'FGrids'. These
  ! define frequency coordinates for use in quantity templates and thence
  ! vector quantities.  Note however that in perhaps the most obvious example
  ! of a quantity having a frequency coordinate - the radiances - the fGrid
  ! information is not used, in these cases the information given in the signal
  ! field describes the frequency.  FGrids are mainly used for quantities such
  ! as baseline and extinction

  ! This is the main data type
  type FGrid_T
    integer :: name                     ! String index of name
    integer :: noChans                  ! Number of `channels'
    integer :: frequencyCoordinate      ! Literal, see below
    real(r8), dimension(:), pointer :: values => NULL() ! The frequencies
  end type FGrid_T
  ! Note - frequency coordinate is a literal, examples are:
  !        l_frequency, l_usbFrequency, l_lsbFrequency, l_intermediateFrequency
  !      - values *must* remain r8, as we're after kHz precision in THz values.


contains ! ===================================== Public procedures =====

  ! ----------------------------------------------- AddFGridToDatabase
  integer function AddFGridToDatabase ( database, item )
    ! Dummy arguments
    type (FGrid_T), dimension(:), pointer :: database
    type (FGrid_T), intent(in) :: item

    ! Local variables
    type (FGrid_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddFGridToDatabase = newSize
  end function AddFGridToDatabase

  ! -------------------------------------------- CreateFGridFromMLSCFInfo
  type (FGrid_T) function CreateFGridFromMLSCFInfo ( name, root ) &
    & result ( fGrid )
    ! Dummy arguments
    integer, intent(in) :: NAME         ! String index of name
    integer, intent(in) :: ROOT         ! Tree root
    
    ! Local variables
    integer :: I,J                      ! Loop counter
    integer :: SON                      ! Tree node
    integer :: UNITS(2)                 ! From expr
    real(r8) :: VALUES(2)               ! From expr

    ! Executable code
    fGrid%name = name
    do i = 2, nsons(root)
      son = subtree(i,root)
      select case ( decoration( subtree(1,son) ) )
      case ( f_coordinate )
        fGrid%frequencyCoordinate = decoration( subtree(2,son) )
      case ( f_values )
        fGrid%noChans = nsons ( son ) - 1
        call Allocate_test ( fGrid%values, fGrid%noChans, 'fGrid%values', &
          & ModuleName )
        do j = 1, fGrid%noChans
          call expr ( subtree ( j+1, son), units, values )
          if ( all ( units(1) /= &
            & (/ phyq_frequency, phyq_dimensionless /) ) ) then
            call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Bad units for fGrid' )
          end if
          fGrid%values(j) = values(1)
        end do
      end select
    end do
    ! Because the parser does such a great job, that's all we need to do here!
  end function CreateFGridFromMLSCFInfo

  ! -------------------------------------------- DestroyFGridContents
  subroutine DestroyFGridContents ( fGrid )
    ! Dummy arguments
    type (FGrid_T), intent(inout) :: fGrid

    ! Executable code
    call Deallocate_test ( fGrid%values, 'fGrid%values', ModuleName )
    fGrid%name = 0
    fGrid%noChans = 0
    fGrid%frequencyCoordinate = l_none
  end subroutine DestroyFGridContents

  ! ----------------------------------------- Destroy FGridDatabase
  subroutine DestroyFGridDatabase ( database )
    ! Dummy arguments
    type (FGrid_T), dimension(:), pointer :: database

    ! Local variables
    integer :: I                        ! Loop counter
    integer :: STATUS                   ! From deallocate

    ! Executable code

    do i = 1, size ( database )
      call DestroyFGridContents ( database(i) )
    end do

    deallocate ( database, STAT=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//'database' )
  end subroutine DestroyFGridDatabase

end module FGrid

! $Log$
! Revision 2.1  2001/10/31 18:36:19  livesey
! First version
!
