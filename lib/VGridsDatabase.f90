! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module VGridsDatabase
!=============================================================================

  use Allocate_Deallocate, only: Deallocate_Test
  use HyperSlabs, only: EssentiallyEqual
  use Intrinsic, only: L_Angle, L_Eta, L_Geodaltitude, L_Gph, L_Integer, &
    & L_None, L_Pressure, L_Theta, L_Zeta, &
    & PHYQ_Dimensionless, PHYQ_Pressure, PHYQ_Zeta, PHYQ_Temperature, &
    & PHYQ_Length, PHYQ_Angle, PHYQ_Invalid
  use MLSKinds, only: RS => R8 ! Real Kind For Surfs
  use MLSMessageModule, only: & ! Message Logging
    & MLSMessage, MLSMSG_Error, &
    & PVMErrorMessage

  implicit none
  private

  ! Define the vGrid data type.  This is used to store all the vGrid
  ! information. Note that this is only relevant for coherent quantities.
  ! Incoherent ones deal with vGrids seperately.

  type, public :: VGrid_T
    integer:: Name                 ! String index of name
    integer :: VerticalCoordinate  ! One of t_vGridCoordinate's literals, or
                                   ! 0 if empty
    integer :: NoSurfs             ! Number of surfaces
    real(rs), dimension(:,:), pointer :: Surfs => NULL()  ! Array of surfaces
                                   ! (actually dimensioned (1:noSurfs,1))
  end type VGrid_T

  ! The VGrids database:
  type(VGrid_t), pointer, save, public :: VGrids(:) => NULL()

  ! Public procedures:
  interface ConvertVGrid
    module procedure ConvertVGrid_inout, ConvertVGrid_sngl, ConvertVGrid_dbl
  end interface ConvertVGrid

  interface DoVGridsMatch
    module procedure DoVGridsMatch_VG
  end interface

  interface Dump
    module procedure Dump_VGrids, Dump_a_VGrid
  end interface Dump

  public :: AddVgridIfNecessary, AddVGridToDatabase
  public :: ConvertVGrid
  public :: DestroyVGridContents, DestroyVGridDatabase
  public :: DoVGridsMatch, DoVGridsMatch_VG
  public :: Dump, Dump_a_VGrid, Dump_VGrids, GetUnitForVerticalCoordinate
  public :: NullifyVGrid
  public :: PVMPackVGrid, PVMUnpackVGrid
  public :: RS

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  private :: not_used_here
  !---------------------------------------------------------------------------

contains


  !-----------------------------------------  AddVGridIfNecessary  -----
  integer function AddVgridIfNecessary ( VGrid, RelErr )
    type(vGrid_t), intent(inout) :: VGrid
    real(rs), intent(in), optional :: RelErr
    ! If there is a vGrid in vGrids that matches vGrid according to
    ! doVGridsMatch, destroy vGrid and return the index of the matching one.
    ! Otherwise, add VGrid to the database of VGrids and return the index
    ! of the added one.

    integer :: I

    if ( associated ( vGrids ) ) then
      do i = 1, size(vGrids)
        if ( doVGridsMatch(vGrid,vGrids(i),relerr) ) then
          call destroyVGridContents ( vGrid )
          addVgridIfNecessary = i
          return
        end if
      end do
    end if

    addVgridIfNecessary = addVGridToDatabase ( vGrids, vGrid )

  end function AddVgridIfNecessary

  !------------------------------------------  AddVGridToDatabase  -----
  integer function AddVGridToDatabase ( DATABASE, ITEM )

  ! This routine adds a vGrid to a database of vGrids, creating the database
  ! if necessary.

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

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

  !----------------------------------------  ConvertVGrid_inout  -----
  subroutine ConvertVGrid_inout ( vGrid, newVerticalCoordinate )

  ! This routine converts the surfaces of a VGrid between vertical coordinates

    ! Dummy arguments

    type (vGrid_T), intent(inout) :: vGrid
    integer, intent(in)           :: newVerticalCoordinate

    ! Executable code
    if ( VGrid%verticalCoordinate == newVerticalCoordinate ) return
    select case (newVerticalCoordinate)
    case (l_GeodAltitude)
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--not able to convert VGrids between vertical coordinates' )
    case (l_GPH)
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--not able to convert VGrids between vertical coordinates' )
    case (l_eta)
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--not able to convert VGrids between vertical coordinates' )
    case (l_theta)
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--not able to convert VGrids between vertical coordinates' )
    case (l_pressure)
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--not able to convert VGrids between vertical coordinates' )
    case default
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--not able to convert VGrids between vertical coordinates' )
    end select

  end subroutine ConvertVGrid_inout

  !----------------------------------------  ConvertVGrid_sngl  -----
  subroutine ConvertVGrid_sngl ( vGrid, newVerticalCoordinate, surfs )

  ! This routine converts the surfaces of a VGrid between vertical coordinates

    ! Dummy arguments

    type (vGrid_T), intent(in) :: vGrid
    integer, intent(in)        :: newVerticalCoordinate
    real, dimension(:)         :: surfs


    ! Executable code
    surfs = vGrid%surfs(:,1)
    if ( VGrid%verticalCoordinate == newVerticalCoordinate ) return
    select case (newVerticalCoordinate)
    case (l_GeodAltitude)
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--not able to convert VGrids between vertical coordinates' )
    case (l_GPH)
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--not able to convert VGrids between vertical coordinates' )
    case (l_eta)
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--not able to convert VGrids between vertical coordinates' )
    case (l_theta)
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--not able to convert VGrids between vertical coordinates' )
    case (l_pressure)
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--not able to convert VGrids between vertical coordinates' )
    case default
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--not able to convert VGrids between vertical coordinates' )
    end select

  end subroutine ConvertVGrid_sngl

  !----------------------------------------  ConvertVGrid_dbl  -----
  subroutine ConvertVGrid_dbl ( vGrid, newVerticalCoordinate, surfs )

  ! This routine converts the surfaces of a VGrid between vertical coordinates

    ! Dummy arguments

    type (vGrid_T), intent(in)     :: vGrid
    integer, intent(in)            :: newVerticalCoordinate
    double precision, dimension(:) :: surfs


    ! Executable code
    surfs = vGrid%surfs(:,1)
    if ( VGrid%verticalCoordinate == newVerticalCoordinate ) return
    select case (newVerticalCoordinate)
    case (l_GeodAltitude)
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--not able to convert VGrids between vertical coordinates' )
    case (l_GPH)
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--not able to convert VGrids between vertical coordinates' )
    case (l_eta)
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--not able to convert VGrids between vertical coordinates' )
    case (l_theta)
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--not able to convert VGrids between vertical coordinates' )
    case (l_pressure)
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--not able to convert VGrids between vertical coordinates' )
    case default
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--not able to convert VGrids between vertical coordinates' )
    end select

  end subroutine ConvertVGrid_dbl

  ! ---------------------------------------  DestroyVGridDatabase  -----
  subroutine DestroyVGridDatabase ( DATABASE )

  ! This subroutine destroys a vGrid database

    use Allocate_Deallocate, only: Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc

    ! Dummy argument
    type (VGrid_T), dimension(:), pointer :: DATABASE

    ! Local variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: S, Status, vgridIndex

    if ( associated(database) ) then
      do vgridIndex = 1, SIZE(database)
        call DestroyVGridContents ( database(vgridIndex) )
      end do
      s = size(database) * storage_size(database) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(database(1)), addr)
      deallocate ( database, stat=status )
      call test_deallocate ( status, ModuleName, "database", s, address=addr )
    end if

  end subroutine DestroyVGridDatabase

  ! -------------------------------------------  DoVGridsMatch_VG  -----
  logical function DoVGridsMatch_VG ( A, B, RelErr )
    ! Returns true if A and B are essentially the same VGrid.
    type (vGrid_T), intent(in) :: A
    type (vGrid_T), intent(in) :: B
    real(rs), intent(in), optional :: RelErr ! "essentially equal" means
      ! "all(abs(a%surfs-b%surfs) <=
      ! relerr*max(maxval(abs(a%surfs)),maxval(abs(b%surfs))))"
    real(rs) :: Test

    ! Executable code
    doVGridsMatch_VG = .false.
    if ( a%verticalCoordinate /= b%verticalCoordinate ) return
    if ( a%noSurfs /= b%noSurfs ) return
    if ( .not. present(relErr) ) then
      if ( any ( .not. essentiallyEqual ( a%surfs, b%surfs ) ) ) return
    else
      test = relerr*max(maxval(abs(a%surfs)),maxval(abs(b%surfs)))
      if ( any(abs(a%surfs-b%surfs) > test) ) return
    end if
    doVGridsMatch_VG = .true.

  end function DoVGridsMatch_VG

  ! -----------------------------------------------  Dump_a_VGrid  -----
  subroutine Dump_a_VGrid ( VGrid, Details, What )
    use Dump_0, only: Dump
    use Intrinsic, only: Lit_Indices
    use Output_M, only: Newline, Output
    use String_Table, only: Display_String
    type(vGrid_T), intent(in) :: VGrid
    integer, intent(in), optional :: Details ! <= 0 => Don't dump arrays
    !                                        ! >0   => Do dump arrays
    !                                        ! Default 1
    character(len=*), intent(in), optional :: What ! Prints ' '//what//':'
    integer :: MyDetails
    myDetails = 1
    if ( present(details) ) myDetails = details
    if ( present(what) ) call output ( ' ' // trim(what) // ':' )
    call output ( ' Name = ' )
    if ( vGrid%name /= 0 ) then
      call display_string ( vgrid%name )
    else
      call output ( '<none>' )
    end if
    call output ( ' noSurfs = ' )
    call output ( vgrid%noSurfs )
    call output ( ' Coordinate = ' )
    call display_string ( lit_indices(vgrid%verticalCoordinate) )
    call newline
    if ( myDetails > 0 ) then
      call dump ( vgrid%surfs(:,1), ' Surfs = ' )
      call newline
    end if
  end subroutine Dump_a_VGrid


  ! ------------------------------------------------  Dump_VGrids  -----
  subroutine Dump_VGrids ( VGrids, Details, Where )

    use MoreTree, only: StartErrorMessage
    use Output_M, only: Output
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

  ! ----------------------------------  GetUnitForVerticalCoordinate  -----
  integer function GetUnitForVerticalCoordinate ( coordinate )
    integer, intent(in) :: coordinate
    ! Excutable code
    select case ( coordinate )
    case ( l_angle )
      GetUnitForVerticalCoordinate = PHYQ_Angle
    case ( l_geodAltitude, l_gph )
      GetUnitForVerticalCoordinate = PHYQ_Length
    case ( l_integer )
      GetUnitForVerticalCoordinate = PHYQ_Dimensionless
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

  ! ----------------------------------------NullifyVGrid -----
  subroutine NullifyVGrid ( IntentionallyNotUsed )
    ! Given a vGrid, nullify all the pointers associated with it
    type ( VGrid_T ), intent(out) :: IntentionallyNotUsed

    ! Executable code isn't necessary because IntentionallyNotUsed is
    ! intent(out) and IntentionallyNotUsed%surfs has default initialization
    ! to NULL()
  end subroutine NullifyVGrid

  ! ------------------------------------------------  PVMPackVGrid ----
  subroutine PVMPackVgrid ( VGRID )
    use PVMIDL, only: PVMIDLPack
    use MorePVM, only: PVMPackStringIndex, PVMPackLitIndex

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

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module VGridsDatabase

! $Log$
! Revision 2.30  2018/04/19 02:00:36  vsnyder
! Compute address for allocate/deallocate tracking.  Remove USE statements for
! unused names.
!
! Revision 2.29  2017/11/03 20:02:31  pwagner
! Most array gymnastics moved from MLSFillValues to HyperSlabs module
!
! Revision 2.28  2015/03/28 01:43:46  vsnyder
! Added stuff to trace allocate/deallocate addresses
!
! Revision 2.27  2014/09/05 00:19:47  vsnyder
! More complete and accurate allocate/deallocate size tracking.  Some
! cannonball polishing.
!
! Revision 2.26  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.25  2008/08/27 19:58:30  vsnyder
! Add PRINT to not_used_here
!
! Revision 2.24  2008/06/06 22:52:21  pwagner
! EssentiallyEqual moved to MLSFillValues
!
! Revision 2.23  2008/06/06 01:54:43  vsnyder
! Get kinds from MLSKinds, not MLSCommon
!
! Revision 2.22  2006/05/02 18:59:50  pwagner
! Added ConvertVGrid, though mostly non-functional for now
!
! Revision 2.21  2006/02/08 21:34:46  vsnyder
! Add a 'what' argument, that just prints, to Dump_a_VGrid
!
! Revision 2.20  2005/06/03 01:54:36  vsnyder
! Make VGrids a public module variable, some cannonball polishing,
! new copyright notice, move Id to not_used_here to avoid cascades.
!
! Revision 2.19  2005/04/19 20:23:57  livesey
! Bug fix in AddVGridIfNecessary for case when this is the first vGrid
!
! Revision 2.18  2005/03/15 23:48:55  pwagner
! PVMERRORMESSAGE now part of MLSMessageModule
!
! Revision 2.17  2005/01/12 03:06:08  vsnyder
! Added AddVGridIfNecessary and relative error test option to DoVGridsMatch.
! Don't try to dump VGrid's name if its string index is zero.
!
! Revision 2.16  2005/01/07 00:38:53  vsnyder
! Call the kind for the Surfs field RS
!
! Revision 2.15  2004/12/13 20:29:36  vsnyder
! Added DoVGridsMatch generic with DoVGridsMatch_VG specific
!
! Revision 2.14  2004/06/17 22:35:10  pwagner
! Added new integer type for vertical coordinate
!
! Revision 2.13  2004/06/03 22:57:48  vsnyder
! Cosmetic changes to account for using VGrid struct for TGrids too
!
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
