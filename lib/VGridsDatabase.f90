! Copyright (c) 1999, California Institute of Technology. ALL RIGHTS RESERVED.
! U.S. Government sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module VGridsDatabase
!=============================================================================

  use Allocate_Deallocate, only: Deallocate_Test
  use Dump_0, only: DUMP
  use MLSCommon, only: R8
  use MLSMessageModule, only: & ! Message logging
    & MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, MLSMSG_Error
  use OUTPUT_M, only: OUTPUT
  use STRING_TABLE, only: DISPLAY_STRING

  implicit NONE
  private

  ! Define the vGrid data type.  This is used to store all the vGrid
  ! information. Note that this is only relevant for coherent quantities. 
  ! Incoherent ones deal with vGrids seperately.

  type, public :: VGrid_T
    integer:: Name                 ! String index of name
    integer :: VerticalCoordinate  ! One of t_vGridCoordinate's literals, or
                                   ! 0 if empty
    integer :: NoSurfs             ! Number of surfaces
    real(r8), dimension(:), pointer :: Surfs => NULL()  ! Array of surfaces
                                   ! (actually dimensioned 1:noSurfs)
  end type VGrid_T

  ! Public procedures:
  interface Dump
    module procedure Dump_VGrids
  end interface Dump

  public :: AddVGridToDatabase, DestroyVGridContents, DestroyVGridDatabase, Dump
  public :: Dump_VGrids

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains

  !------------------------------------------  AddVGridToDatabase  -----
  integer function AddVGridToDatabase ( DATABASE, ITEM )

  ! This routine adds a vGrid to a database of vGrids, creating the database
  ! if necessary.

    ! Dummy arguments
    type (VGrid_T), dimension(:), pointer :: DATABASE
    type (VGrid_T), intent(in) :: ITEM

    ! Local variables
    type (VGrid_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddVGridToDatabase = newSize
  end function AddVGridToDatabase

  !----------------------------------------  DestroyVGridContents  -----
  subroutine DestroyVGridContents ( vGrid )

  ! This routine destroys the array information created with the vGrid

    ! Dummy arguments

    type (vGrid_T), intent(inout) :: vGrid

    ! Executable code

    vGrid%noSurfs = 0
    vGrid%verticalCoordinate = 0

    call deallocate_test ( vGrid%surfs, "vGrid%surfs", ModuleName )

  end subroutine DestroyVGridContents

  ! ---------------------------------------  DestroyVGridDatabase  -----
  subroutine DestroyVGridDatabase ( DATABASE )

  ! This subroutine destroys a vGrid database

    ! Dummy argument
    type (VGrid_T), dimension(:), pointer :: DATABASE

    ! Local variables
    integer :: vgridIndex, Status

    if ( associated(database) ) then
      do vgridIndex = 1, SIZE(database)
        call DestroyVGridContents ( database(vgridIndex) )
      end do
      deallocate ( database, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Deallocate // "database" )
    end if
  end subroutine DestroyVGridDatabase

  ! ------------------------------------------------  Dump_VGrids  -----
  subroutine Dump_VGrids ( VGrids, Lit_Indices, Details )
    type(vGrid_T), intent(in) :: VGrids(:)             ! The database
    integer, intent(in), dimension(:) :: Lit_Indices   ! From init_tables
    integer, intent(in), optional :: Details ! <= 0 => Don't dump arrays
    !                                        ! >0   => Do dump arrays
    !                                        ! Default 1
    integer :: I, MyDetails
    myDetails = 1
    if ( present(details) ) myDetails = details
    call output ( 'VGRIDS: SIZE = ' )
    call output ( size(vgrids), advance='yes' )
    do i = 1, size(vgrids)
      call output ( i, 4 )
      call output ( ': Name = ' )
      call display_string ( vgrids(i)%name )
      call output ( ' noSurfs = ' )
      call output ( vgrids(i)%noSurfs )
      call output ( ' verticalCoordinate = ' )
      call display_string ( lit_indices(vgrids(i)%verticalCoordinate) )
      if ( details > 0 ) call dump ( vgrids(i)%surfs, ' Surfs = ' )
    end do
  end subroutine Dump_VGrids

end module VGridsDatabase

! $Log$
! Revision 2.2  2001/04/11 00:03:52  vsnyder
! Improve 'dump'
!
! Revision 2.1  2001/04/07 01:54:08  vsnyder
! Initial Commit
!
