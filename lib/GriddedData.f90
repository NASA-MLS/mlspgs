! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module GriddedData ! Contains the derived TYPE GriddedData_T

  use MLSCommon, only: R8, LINELEN, NAMELEN, RP
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ALLOCATE, MLSMSG_ERROR, &
    & MLSMSG_DEALLOCATE
  use Toggles, only: TOGGLE, GEN
  use Output_m, only: OUTPUT
  use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
  use Trace_M, only: TRACE_BEGIN, TRACE_END

  implicit none
  public

  private :: Id,ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  character(LEN=130) :: id = & 
    "$Id$"
  character(LEN=*), parameter :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  public::GriddedData_T, SetupNewGriddedData, DestroyGriddedData, &
    & AddGriddedDataToDatabase, DestroyGriddedDataDatabase, Dump

    logical, private, parameter :: MAYDUMPFIELDVALUES = .false.

  interface DUMP
    module procedure DumpGriddedData
    module procedure DumpGriddedDatabase
  end interface

  ! These are 'enumerated types' consistent with hph's
  ! work in l3ascii_read_field
  public::v_is_pressure,v_is_altitude,v_is_gph,v_is_theta

  integer, parameter :: v_is_pressure = 1
  integer, parameter :: v_is_altitude = v_is_pressure+1
  integer, parameter :: v_is_gph = v_is_altitude+1
  integer, parameter :: v_is_theta = v_is_gph+1

  ! This type reflects the format of the Level 3 ASCII files, though note that
  ! these files can store multiple quantities such as these.

  type GriddedData_T

    ! First the comment line(s) from the relevant input file
    character (LEN=LineLen), pointer, dimension(:) :: fileComments => NULL()

    ! Now the name, description and units information
    character (LEN=LineLen) :: sourceFileName ! Input file name
    character (LEN=NameLen) :: quantityName ! From input file
    character (LEN=LineLen) :: description ! Quantity description
    character (LEN=NameLen) :: units ! Units for quantity

    ! Now define the various coordinate systems, first vertical
    integer :: verticalCoordinate ! An 'enumerated' type
    integer :: noHeights         ! Number of surfaces
    real (rp), pointer, dimension(:) :: heights  => NULL()
    ! Surfaces (e.g. pressures etc.) [noHeights]

    ! Now the latitudinal coordinate
    logical :: equivalentLatitude       ! If set, coordinate is equivalent latitude
    integer :: noLats                   ! Number of latitudes
    real (rp), pointer, dimension(:) :: lats => NULL() ! Latitudes [noLats]
    integer :: noLons                   ! Number of longitudes
    real (rp), pointer, dimension(:) :: lons => NULL() ! Longitudes [noLons]
    integer noLsts                      ! Number of local times
    real (rp), pointer, dimension(:) :: lsts => NULL() ! Local times [noLsts]
    integer noSzas                      ! Number of solar zenith angles
    real (rp), pointer, dimension(:) :: szas => NULL() ! Zenith angles [noSzas]
    integer noDates                     ! Number of dates in data
    real (rp), pointer, dimension(:) :: dateStarts => NULL()
    ! Starting dates in SDP toolkit format
    real (rp), pointer, dimension(:) :: dateEnds => NULL()
    ! Ending dates in SDP toolkit format

    ! The data itself.  This is stored as
    !  [noHeights, noLats, noLons, noLsts, noSzas, noDates]
    real (rp), pointer, dimension(:,:,:,:,:,:) :: field => NULL()

  end type GriddedData_T

! ============================================================================
contains

  ! ---------------------------------------------------AddGridTemplateToDatabase --
  integer function AddGriddedDataToDatabase(database,item)
  ! This subroutine adds a quantity template to a database, or creates the
  ! database if it doesn't yet exist

    ! Dummy arguments
    type (GriddedData_T), dimension(:), pointer :: database
    type (GriddedData_T), intent(in) :: item

    ! Local variables
    type (GriddedData_T), dimension(:), pointer :: tempDatabase

    ! Executable code

    include "addItemToDatabase.f9h"
    AddGriddedDataToDatabase = newSize

  end function AddGriddedDataToDatabase

  ! ------------------------------------------------ DestroyGriddedData ------

  subroutine DestroyGriddedData(qty)
    ! This subroutine destroys a quantity template

    ! Dummy argument
    type (GriddedData_T), intent(INOUT) :: qty

    ! Local variables
    integer :: STATUS

    ! Executable code
    call Deallocate_test ( qty%heights, "qty%heights", ModuleName )
    call Deallocate_test ( qty%lats, "qty%lats", ModuleName )
    call Deallocate_test ( qty%lons, "qty%lons", ModuleName )
    call Deallocate_test ( qty%lsts, "qty%lsts", ModuleName )
    call Deallocate_test ( qty%szas, "qty%szas", ModuleName )

    ! Now the temporal coordinates
    call Deallocate_test ( qty%dateStarts, "qty%dateStarts", ModuleName )
    call Deallocate_test ( qty%dateEnds, "qty%dateEnds", ModuleName )

    ! Now the data itself
    deallocate(qty%field, STAT=status)

    if (status /= 0) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//"field")

  end subroutine DestroyGriddedData

  ! ------------------------------------------- DestroyGridTemplateDatabase --
  subroutine DestroyGriddedDataDatabase(database)
  ! This subroutine destroys a quantity template database

    ! Dummy argument
    type (GriddedData_T), dimension(:), pointer :: database

    ! Local variables
    integer :: qtyIndex, status

    if ( toggle(gen) ) call trace_begin ( "DestroyGridTemplateDatabase" )

    if (associated(database)) then
      do qtyIndex=1,size(database)
        call DestroyGriddedData(database(qtyIndex))
      enddo
      deallocate(database, stat=status)
      if (status /= 0) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate//"database")
    endif
    if ( toggle(gen) ) then
      call trace_end ( "DestroyGridTemplateDatabase" )
    end if
  end subroutine DestroyGriddedDataDatabase

  ! --------------------------------  DumpGriddedDatabase  -----
  subroutine DumpGriddedDatabase(GriddedData, Details)
    use Dump_0, only: Dump
    ! Imitating what dump_pointing_grid_database does, but for gridded data
    ! which may come from climatology, ncep, dao

    type (GriddedData_T), dimension(:), pointer :: GriddedData 

    integer, intent(in), optional :: DETAILS
    
    ! Local Variables
    integer            :: i

    if ( .not. associated(GriddedData)) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, 'Gridded database still null')

    call output ( '============ Gridded Data Base ============', advance='yes' )
    call output ( ' ', advance='yes' )
    call output ( 'database: a priori grids: SIZE = ' )
    call output ( size(GriddedData), advance='yes' )
    if ( size(GriddedData) < 1 ) return
    do i = 1, size(GriddedData)

      call output ( 'item number ' )
      call output ( i, advance='yes' )

      call DumpGriddedData(GriddedData(i), Details)
    end do ! i
  end subroutine DumpGriddedDatabase

  ! --------------------------------  DumpGriddedData  -----
  subroutine DumpGriddedData(GriddedData, Details)
    use Dump_0, only: Dump

    ! Imitating what dump_pointing_grid_database does, but for gridded data
    ! which may come from climatology, ncep, dao
    type (GriddedData_T) :: GriddedData 
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump multidim arrays
    !                                        ! -1 Skip even 1-d arrays
    !                                        ! -2 Skip all but name
    !                                        ! >0 Dump even multi-dim arrays
    !                                        ! Default 1

    ! Local Variables
    integer :: MYDETAILS

    ! Executable code
    myDetails = 1
    if ( present(details) ) myDetails = details

    call output('Gridded quantity name ' // GriddedData%quantityName, advance='yes')
      if ( myDetails < -1 ) return
    call output('description ' // GriddedData%description, advance='yes')
    call output('units ' // GriddedData%units, advance='yes')

    call output ( ' ************ Geometry ********** ' ,advance='yes')

    call output ( ' Vertical coordinate = ' )
    call output ( GriddedData%verticalCoordinate, advance='yes' )
    call output ( ' No. of heights = ' )
    call output ( GriddedData%noHeights, advance='yes' )
    if ( myDetails >= 0 ) call dump ( GriddedData%heights, &
      & '    Heights =' )

    call output ( ' Equivalent latitude = ' )
    call output ( GriddedData%equivalentLatitude, advance='yes' )
    call output ( ' No. of latitudes = ' )
    call output ( GriddedData%noLats, advance='yes' )
    if ( myDetails >= 0 ) call dump ( GriddedData%lats, &
      & '    latitudes =' )

    call output ( ' No. of longitudes = ' )
    call output ( GriddedData%noLons, advance='yes' )
    if ( myDetails >= 0 ) call dump ( GriddedData%lons, &
      & '    longitudes =' )

    call output ( ' No. of local times = ' )
    call output ( GriddedData%noLsts, advance='yes' )
    if ( myDetails >= 0 ) call dump ( GriddedData%lsts, &
      & '    local times =' )

    call output ( ' No. of solar zenith angles = ' )
    call output ( GriddedData%noSzas, advance='yes' )
    if ( myDetails >= 0 ) call dump ( GriddedData%szas, &
      & '    solar zenith angles =' )

    call output ( ' No. of dates = ' )
    call output ( GriddedData%noDates, advance='yes' )
    if ( myDetails >= 0 ) call dump ( GriddedData%dateStarts, &
      & '    starting dates =' )
    if ( myDetails >= 0 ) call dump ( GriddedData%dateEnds, &
      & '    ending dates =' )

    if ( MAYDUMPFIELDVALUES .and. myDetails > 0 ) then
      call output ( ' ************ tabulated field values ********** ' ,advance='yes')

      ! No dump for 6-dimensional double arrays yet, anyway
      !     call dump ( GriddedData%field, &
      !      & '    gridded field values =' )
      call output ( ' *(Sorry, dump_6d_double not yet coded)* ' ,advance='yes')
    endif

  end subroutine DumpGriddedData

  ! ------------------------------------------- SetupNewGriddedData ---------
  subroutine SetupNewGriddedData(qty, source, noHeights, noLats, &
    & noLons, noLsts, noSzas, noDates)
  ! This first routine sets up a new quantity template according to the user
  ! input.  This may be based on a previously supplied template (with possible
  ! modifications), or created from scratch.
    ! Dummy arguments
    type (GriddedData_T), intent(OUT) :: QTY ! Result
    type (GriddedData_T), optional, intent(in) :: SOURCE ! Template
    integer, optional, intent(in) :: NOHEIGHTS, NOLATS, NOLONS, NOLSTS, NOSZAS, NODATES

    ! Local variables
    integer :: status           ! Status from allocates etc.

    ! Executable code

    ! First, if we have a template setup according to that
    if (present(source)) then
      qty%noHeights=source%noHeights
      qty%noLats=source%noLats
      qty%noLons=source%noLons
      qty%noLsts=source%noLsts
      qty%noSzas=source%noSzas
      qty%noDates=source%noDates
    else ! We have no template, setup a very bare quantity
      qty%noHeights=1
      qty%noLats=1
      qty%noLons=1
      qty%noLsts=1
      qty%noSzas=1
      qty%noDates=1
    endif

    ! Now, see if the user asked for modifications to this
    if (present(noHeights)) qty%noHeights=noHeights
    if (present(noLats)) qty%noLats=noLats
    if (present(noLons)) qty%noLons=noLons
    if (present(noLsts)) qty%noLsts=noLsts
    if (present(noSzas)) qty%noSzas=noSzas
    if (present(noDates)) qty%noDates=noDates

    ! First the vertical/horizontal coordinates
    call Allocate_test ( qty%heights, qty%noHeights, "qty%heights", ModuleName )
    call Allocate_test ( qty%lats, qty%noLats, "qty%lats", ModuleName )
    call Allocate_test ( qty%lons, qty%noLons, "qty%lons", ModuleName )
    call Allocate_test ( qty%lsts, qty%noLsts, "qty%lsts", ModuleName )
    call Allocate_test ( qty%szas, qty%noSzas, "qty%szas", ModuleName )

    ! Now the temporal coordinates
    call Allocate_test ( qty%dateStarts, qty%noDates, "qty%dateStarts", ModuleName )
    call Allocate_test ( qty%dateEnds, qty%noDates, "qty%dateEnds", ModuleName )

    ! Now the data itself
    allocate(qty%field(qty%noHeights, qty%noLats, qty%noLons,  &
      qty%noLsts, qty%noSzas, qty%noDates), STAT=status)
    if ( status /= 0 ) call MLSMessage(MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'qty%field')
  
  end subroutine SetupNewGriddedData

end module GriddedData

!
! $Log$
! Revision 2.13  2001/10/26 23:17:14  pwagner
! Provides a single dump module interface and details
!
! Revision 2.12  2001/09/10 23:36:56  livesey
! Expanded and tidied up.  Stuff that was in ncep_dao is now here
!
! Revision 2.11  2001/04/10 20:04:54  livesey
! Changed some NameLens to LineLens
!
! Revision 2.10  2001/03/30 00:24:40  pwagner
! Added sourceFileName
!
! Revision 2.9  2001/03/15 21:28:08  pwagner
! Now only the data type
!
! Revision 2.8  2001/03/15 00:37:13  pwagner
! Still not complete; missing gdrdfld
!
! Revision 2.7  2001/03/14 00:32:47  pwagner
! More changes--still wrong, though
!
! Revision 2.6  2001/03/10 00:33:16  pwagner
! Some corrections in ReadGriddedData
!
! Revision 2.5  2001/03/09 01:02:55  pwagner
! Fixed announce_error
!
! Revision 2.4  2001/03/08 01:08:35  pwagner
! Added announce_error
!
! Revision 2.3  2001/03/07 01:03:19  pwagner
! ReadGriddedData added
!
! Revision 2.2  2001/02/21 00:36:43  pwagner
! l3ascii_read_field now has eof as intent(out) arg
!
! Revision 2.1  2001/02/20 21:51:39  pwagner
! Functions absorbed from gridded_data_module
!
! Revision 2.0  2000/09/05 17:41:05  dcuddy
! Change revision to 2.0
!
! Revision 1.5  2000/06/20 22:19:22  lungu
! Changed DOUBBLE PRECISION to REAL (r8).
!
