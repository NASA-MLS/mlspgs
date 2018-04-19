! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module FGrid                    ! Frequency grid information

  use Allocate_Deallocate, only: Allocate_test, Deallocate_test
  use dump_0, only: dump
  use EXPR_M, only: EXPR
  use Intrinsic, only:  L_None, PHYQ_DIMENSIONLESS, PHYQ_FREQUENCY, &
    & L_CHANNEL, LIT_INDICES
  use Init_tables_module, only: F_Coordinate, F_Values
  use MLSKinds, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
  use Output_M, only: Output
  use STRING_TABLE, only: DISPLAY_STRING
  use Tree, only: DECORATION, NSONS, SUBTREE

  implicit none
  private

  public :: FGrid_T, AddFGridToDatabase, CreateFGridFromMLSCFInfo, &
    & DestroyFGridContents, DestroyFGridDatabase, NullifyFgrid, dump

!---------------------------- RCS Module Info ------------------------------
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

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

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
    use TOGGLES, only: GEN, LEVELS, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    ! Dummy arguments
    integer, intent(in) :: NAME         ! String index of name
    integer, intent(in) :: ROOT         ! Tree root

    ! Local variables
    integer :: I,J                      ! Loop counter
    integer :: Me = -1                  ! String index for trace
    integer :: SON                      ! Tree node
    integer :: UNITS(2)                 ! From expr
    real(r8) :: VALUES(2)               ! From expr

    ! Executable code
    call trace_begin ( me, "CreateFGridFromMLSCFInfo", &
      & cond=toggle(gen) .and. levels(gen) > 0 )
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
    if ( fGrid%frequencyCoordinate == l_channel ) call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & 'Did you mean to use channel as fGrid coordinate?' )

    ! Because the parser does such a great job, that's all we need to do here!
    call trace_end ( "CreateFGridFromMLSCFInfo", &
      & cond=toggle(gen) .and. levels(gen) > 0 )
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

    use Allocate_Deallocate, only: Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type (FGrid_T), dimension(:), pointer :: database

    ! Local variables
    integer(c_intptr_t) :: Addr         ! for tracing
    integer :: I                        ! Loop counter
    integer :: S                        ! Size in bytes of object to deallocate
    integer :: STATUS                   ! From deallocate

    ! Executable code
    if (.not. associated(database) ) return
    do i = 1, size ( database )
      call DestroyFGridContents ( database(i) )
    end do

    s = size(database) * storage_size(database) / 8
    addr = 0
    if ( s > 0 ) addr = transfer(c_loc(database(1)), addr)
    deallocate ( database, STAT=status )
    call test_deallocate ( status, ModuleName, 'database', s, address=addr )

  end subroutine DestroyFGridDatabase

  ! -------------------------------------------- DumpFGrid
  subroutine DumpFGrid ( fGrid )
    ! Dummy arguments
    type (FGrid_T), intent(in) :: fGrid
    integer :: ierr

    ! Executable code
    call output('FGrid name: ', advance='no')
    if ( fgrid%name > 0 ) then
        call display_string ( fgrid%name, advance='yes', ierr=ierr )
        if ( ierr /= 0 ) call output ( '(not found in string table)')
    else
        call output('(unknown)', advance='yes')
    endif
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

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module FGrid

! $Log$
! Revision 2.17  2018/04/19 01:14:16  vsnyder
! Remove USE statements for unused names
!
! Revision 2.16  2015/03/28 02:31:30  vsnyder
! Added stuff to trace allocate/deallocate addresses
!
! Revision 2.15  2014/09/05 00:53:20  vsnyder
! More complete and accurate allocate/deallocate size tracking.  Add some
! tracing.
!
! Revision 2.14  2014/09/05 00:49:06  vsnyder
! EmpiricalGeometry.f90 -- Wrong comment
!
! Revision 2.13  2011/08/20 01:04:33  vsnyder
! IERR needs a value before it can be referenced
!
! Revision 2.12  2010/02/17 20:27:26  honghanh
! Fix Dump subroutine to print FGrid when it has no name
!
! Revision 2.11  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.10  2009/03/05 16:20:17  pwagner
! May now use channel as FGrid coordinate
!
! Revision 2.9  2005/06/22 18:57:01  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.8  2004/12/27 23:05:14  vsnyder
! Remove unreferenced use names
!
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
