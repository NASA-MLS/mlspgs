! Copyright (c) 1999, California Institute of Technology. ALL RIGHTS RESERVED.
! U.S. Government sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module VGridsDatabase
!=============================================================================

  use Allocate_Deallocate, only: Deallocate_Test
  use MLSCommon, only: R8
  use MLSMessageModule, only: & ! Message logging
    & MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, MLSMSG_Error

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
    real(r8), dimension(:,:), pointer :: Surfs => NULL()  ! Array of surfaces
                                   ! (actually dimensioned 1:noSurfs)
  end type VGrid_T

  ! Public procedures:
  interface Dump
    module procedure Dump_VGrids, Dump_a_VGrid
  end interface Dump

  public :: AddVGridToDatabase, DestroyVGridContents, DestroyVGridDatabase
  public :: Dump, Dump_a_VGrid, Dump_VGrids, GetUnitForVerticalCoordinate
  public :: NullifyVGrid
  public :: PVMPackVGrid, PVMUnpackVGrid

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  private :: not_used_here 
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

  ! -----------------------------------------------  Dump_a_VGrid  -----
  subroutine Dump_a_VGrid ( VGrid, Details )
    use Dump_0, only: DUMP
    use Intrinsic, only: Lit_Indices
    use OUTPUT_M, only: OUTPUT
    use STRING_TABLE, only: DISPLAY_STRING
    type(vGrid_T), intent(in) :: VGrid
    integer, intent(in), optional :: Details ! <= 0 => Don't dump arrays
    !                                        ! >0   => Do dump arrays
    !                                        ! Default 1
    integer :: MyDetails
    myDetails = 1
    if ( present(details) ) myDetails = details
    call output ( ' Name = ' )
    call display_string ( vgrid%name )
    call output ( ' noSurfs = ' )
    call output ( vgrid%noSurfs )
    call output ( ' verticalCoordinate = ' )
    call display_string ( lit_indices(vgrid%verticalCoordinate) )
    if ( myDetails > 0 ) call dump ( vgrid%surfs(:,1), ' Surfs = ' )
  end subroutine Dump_a_VGrid


  ! ------------------------------------------------  Dump_VGrids  -----
  subroutine Dump_VGrids ( VGrids, Details, Where )

    use MoreTree, only: StartErrorMessage
    use OUTPUT_M, only: OUTPUT
    type(vGrid_T), pointer :: VGrids(:)             ! The database
    integer, intent(in), optional :: Details ! <= 0 => Don't dump arrays
    !                                        ! >0   => Do dump arrays
    !                                        ! Default 1
    integer, intent(in), optional :: Where   ! Tree node index

    integer :: I
    if ( associated(vGrids) ) then
      call output ( 'VGRIDS: SIZE = ' )
      call output ( size(vgrids), advance='yes' )
      do i = 1, size(vgrids)
        call output ( i, 4 )
        call dump_a_vGrid ( vgrids(i), details )
      end do
    else
      if ( present(where) ) call startErrorMessage ( where )
      call output ( 'No VGrids to dump.', advance='yes' )
    end if
  end subroutine Dump_VGrids

  ! ------------------------------------------------  GetUnitForVGrid  -----
  integer function GetUnitForVerticalCoordinate ( coordinate )
    use intrinsic, only: L_ANGLE, L_GEODALTITUDE, L_GPH, L_NONE, &
      & L_PRESSURE, L_THETA, L_ZETA, PHYQ_Pressure, PHYQ_Zeta, PHYQ_Temperature, &
      & PHYQ_Length, PHYQ_Angle, PHYQ_Invalid
    integer, intent(in) :: coordinate
    ! Excutable code
    select case ( coordinate )
    case ( l_angle )
      GetUnitForVerticalCoordinate = PHYQ_Angle
    case ( l_geodAltitude, l_gph )
      GetUnitForVerticalCoordinate = PHYQ_Length
    case ( l_none )
      GetUnitForVerticalCoordinate = PHYQ_Invalid
    case ( l_pressure ) 
      GetUnitForVerticalCoordinate = PHYQ_Pressure
    case ( l_theta ) 
      GetUnitForVerticalCoordinate = PHYQ_Temperature
    case ( l_zeta )
      GetUnitForVerticalCoordinate = PHYQ_Zeta
    case default
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Invalid vGrid type' )
    end select
  end function GetUnitForVerticalCoordinate

  ! ------------------------------------------------  PVMPackVGrid ----
  subroutine PVMPackVgrid ( VGRID )
    use PVMIDL, only: PVMIDLPack
    use MorePVM, only: PVMPackStringIndex, PVMPackLitIndex
    use PVM, only: PVMErrorMessage

    ! Dummy argument
    type ( VGrid_T), intent(in) :: VGRID

    ! Local variables
    integer :: INFO                     ! Flag from PVM

    ! Executable code

    call PVMPackStringIndex ( vGrid%name, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Packing Vgrid name' )
    call PVMPackLitIndex ( vGrid%verticalCoordinate, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Packing Vgrid coordinate' )

    call PVMIDLPack ( vGrid%noSurfs, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Packing Vgrid size' )
    if ( associated ( vGrid%surfs ) ) then
      call PVMIDLPack ( vGrid%surfs, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, 'Packing Vgrid surfaces' )
    else
      if ( vGrid%noSurfs > 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Asked to pack a strange vGrid' )
    end if
    
  end subroutine PVMPackVgrid

  ! ------------------------------------------------  PVMUnpackVGrid ----
  subroutine PVMUnpackVgrid ( VGRID )
    use PVMIDL, only: PVMIDLUnpack
    use MorePVM, only: PVMUnpackStringIndex, PVMUnpackLitIndex
    use PVM, only: PVMErrorMessage
    use Allocate_Deallocate, only: Allocate_test

    ! Dummy argument
    type ( VGrid_T), intent(out) :: VGRID

    ! Local variables
    integer :: INFO                     ! Flag from PVM

    ! Executable code

    call PVMUnpackStringIndex ( vGrid%name, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Unpacking Vgrid name' )
    call PVMUnpackLitIndex ( vGrid%verticalCoordinate, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Unpacking Vgrid coordinate' )

    call PVMIDLUnpack ( vGrid%noSurfs, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Unpacking Vgrid size' )
    call Allocate_test ( vGrid%surfs, vGrid%noSurfs, 1, 'vGrid%surfs', ModuleName )
    if ( vGrid%noSurfs > 0 ) then
      call PVMIDLUnpack ( vGrid%surfs, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, 'Unpacking Vgrid surfaces' )
    end if
    
  end subroutine PVMUnpackVgrid

  ! ----------------------------------------NullifyVGrid -----
  subroutine NullifyVGrid ( V )
    ! Given a vGrid, nullify all the pointers associated with it
    type ( VGrid_T ), intent(out) :: V

    ! Executable code
    nullify ( v%surfs )
  end subroutine NullifyVGrid

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module VGridsDatabase

! $Log$
! Revision 2.12  2004/05/29 02:46:39  vsnyder
! Fix a bug in Dump_a_VGrid, some cannonball-polishing
!
! Revision 2.11  2004/05/22 02:27:07  vsnyder
! Add Dump_a_VGrid
!
! Revision 2.10  2003/09/15 23:19:04  vsnyder
! Remove unused local variables and USEs
!
! Revision 2.9  2003/08/16 00:32:50  vsnyder
! Get PHYQ_... directly from Intrinsic instead of indirectly via Units
!
! Revision 2.8  2003/06/20 19:33:53  pwagner
! Quanities now share grids stored separately in databses
!
! Revision 2.7  2002/11/22 12:55:12  mjf
! Added nullify routine(s) to get round Sun's WS6 compiler not
! initialising derived type function results.
!
! Revision 2.6  2002/10/08 00:09:15  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.5  2002/10/05 00:41:12  livesey
! Added pvm pack and unpack stuff
!
! Revision 2.4  2002/08/20 19:19:32  livesey
! Added GetUnitForVerticalCoordinate
!
! Revision 2.3  2001/04/26 02:33:03  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.2  2001/04/11 00:03:52  vsnyder
! Improve 'dump'
!
! Revision 2.1  2001/04/07 01:54:08  vsnyder
! Initial Commit
!
