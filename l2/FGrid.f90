! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module FGrid                    ! Frequency grid information

  use Allocate_Deallocate, only: Allocate_test, Deallocate_test
  use dump_0, only: dump
  use EXPR_M, only: EXPR
  use Intrinsic, only: L_Frequency, L_IntermediateFrequency, &
    & L_LSBFrequency, L_None, L_USBFrequency, PHYQ_DIMENSIONLESS, PHYQ_FREQUENCY, &
    & L_CHANNEL, LIT_INDICES
  use Init_tables_module, only: F_Coordinate, F_Values
  use MLSCommon, only: r8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Error, &
    & MLSMSG_Deallocate
  use Output_M, only: Output
  use STRING_TABLE, only: DISPLAY_STRING
  use Tree, only: DECORATION, NSONS, SUBTREE

  implicit none
  private

  public :: FGrid_T, AddFGridToDatabase, CreateFGridFromMLSCFInfo, &
    & DestroyFGridContents, DestroyFGridDatabase, NullifyFgrid, dump

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
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
  interface DUMP
    module procedure DUMPFGRID
    module procedure DUMPFGRIDS
  end interface


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
    call nullifyFGrid ( fGrid ) ! for Sun's still useless compiler
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
    if ( fGrid%frequencyCoordinate == l_channel ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Cannot actually use channel as fGrid coordinate' )

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
    if (.not. associated(database) ) return
    do i = 1, size ( database )
      call DestroyFGridContents ( database(i) )
    end do

    deallocate ( database, STAT=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//'database' )
  end subroutine DestroyFGridDatabase

  ! -------------------------------------------- DumpFGrid
  subroutine DumpFGrid ( fGrid )
    ! Dummy arguments
    type (FGrid_T), intent(in) :: fGrid

    ! Executable code
    call output('FGrid name: ', advance='no')
    call display_string ( fgrid%name, advance='yes' )
    call output('Number channels: ', advance='no')
    call output(fGrid%noChans, advance='yes')
    call output('Frequency coord: ', advance='no')
    call display_string ( lit_indices(fGrid%frequencyCoordinate), &
        &             strip=.true., advance='yes' )
    if ( associated(fgrid%values) ) then
      call dump(fgrid%values, 'fgrid values')
    else
      call output('(values not associated)', advance='yes')
    endif
  end subroutine DumpFGrid

  ! -------------------------------------------- DumpFGrids
  subroutine DumpFGrids ( fGrids, destroy )
    ! Dummy arguments
    type (FGrid_T), dimension(:), intent(inout) :: fGrids
    logical, optional, intent(in) :: destroy
    integer :: I
    logical :: myDestroy
    myDestroy = .false.
    if ( present(destroy) ) myDestroy = destroy
    call output ( size(fgrids), before='FGRIDS: SIZE = ', advance='yes' )
    do i = 1, size(fgrids)
      call output ( i, 4, after=': ' )
      call dump ( fgrids(i) )
    end do
    if ( .not. myDestroy ) return
    do i = 1, size(fgrids)
      call output ( i, 4, after=': ' )
      call destroyFGridContents ( fgrids(i) )
    end do
  end subroutine DumpFGrids

  ! ----------------------------------------NullifyFGrid -----
  subroutine NullifyFGrid ( F )
    ! Given a fGrid, nullify all the pointers associated with it
    type ( FGrid_T ), intent(out) :: F

    ! Executable code
    nullify ( f%values )
  end subroutine NullifyFGrid

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module FGrid

! $Log$
! Revision 2.7  2004/08/26 18:51:48  pwagner
! Added dump methods
!
! Revision 2.6  2004/01/23 05:38:06  livesey
! Tidied up handling of l_channel
!
! Revision 2.5  2003/08/15 23:58:20  vsnyder
! Get PHYQ_... directly from Intrinsic instead of indirectly via Units
!
! Revision 2.4  2002/11/22 12:17:48  mjf
! Added nullify routine(s) to get round Sun's WS6 compiler not
! initialising derived type function results.
!
! Revision 2.3  2002/10/08 17:36:20  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.2  2001/11/01 00:23:14  livesey
! Add check in DestroyFGridDatabase to not do anything if it's empty
!
! Revision 2.1  2001/10/31 18:36:19  livesey
! First version
!
