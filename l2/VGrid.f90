! Copyright (c) 1999, California Institute of Technology. ALL RIGHTS RESERVED.
! U.S. Government sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module vGrid                    ! Definitions for vGrids in vector quantities
!=============================================================================

  use EXPR_M, only: EXPR
  use INIT_TABLES_MODULE, only: F_COORDINATE, F_NUMBER, F_PER_DECADE, F_START, &
    & F_STOP, F_TYPE, F_VALUES, FIELD_FIRST, FIELD_INDICES, FIELD_LAST, &
    & L_ANGLE, L_EXPLICIT, L_GEODALTITUDE, L_GPH, L_LINEAR, L_LOGARITHMIC, &
    & L_NONE, L_PRESSURE, L_THETA, L_ZETA, LIT_INDICES, PHYQ_Angle, &
    & PHYQ_Dimensionless, PHYQ_Length, PHYQ_Pressure, PHYQ_Temperature
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: R8       ! General constants etc.
  use MLSMessageModule, only: & ! Message logging
    & MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, MLSMSG_Error
  use OUTPUT_M, only: OUTPUT
  use STRING_TABLE, only: DISPLAY_STRING
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TOGGLES, only: GEN, TOGGLE
  use TREE, only: DECORATION, DUMP_TREE_NODE, NSONS, SOURCE_REF, SUBTREE

  implicit none
  public

  private :: ID, ModuleName
  !------------------------------- RCS Ident Info ---------------------------
  character (len=130) :: Id= "$Id$"
  character (len=*), parameter :: ModuleName="$RCSfile$"
  !--------------------------------------------------------------------------

  ! Define the vGrid data type.  This is used to store all the vGrid
  ! information. Note that this is only relevant for coherent quantities. 
  ! Incoherent ones deal with vGrids seperately.

  integer, public, parameter :: VC_Invalid = 0
  type vGrid_T
    integer:: NAME                 ! String index of name
    integer :: verticalCoordinate  ! One of t_vGridCoordinate's literals, or
                                   ! VC_Invalid if empty
    integer :: noSurfs             ! Number of surfaces
    real(r8), dimension(:), pointer :: surfs  ! Array of surfaces
                                   ! (actually dimensioned noSurfs)
  end type vGrid_T

! -----     Private declarations     ---------------------------------

  integer, private :: ERROR

! Error codes for "announce_error"
  integer, private, parameter :: DuplicateField = 1
  integer, private, parameter :: ExtraIf = 2
  integer, private, parameter :: InconsistentSizes = 3
  integer, private, parameter :: InconsistentUnits = 4
  integer, private, parameter :: NoCoordField = 5
  integer, private, parameter :: NotPositive = 6
  integer, private, parameter :: NoTypeField = 7
  integer, private, parameter :: RequiredIf = 8
  integer, private, parameter :: RequireExplicit = 9
  integer, private, parameter :: StartStopUnits = 10
  integer, private, parameter :: TooFew = 11
  integer, private, parameter :: Unitless = 12
  integer, private, parameter :: UnitsPressure = 13
  integer, private, parameter :: WrongUnits = 14

contains ! =====     Public Procedures     =============================

  !------------------------------------  CreateVGridFromMLSCFInfo  -----
  type(vGrid_T) function CreateVGridFromMLSCFInfo ( name, root ) &
    & result ( vGrid )

    ! This routine creates a vGrid according to user supplied information
    ! in the l2cf.

    ! Dummy arguments
    integer, intent(in) :: NAME    ! String index of name
    integer, intent(in) :: ROOT    ! Root of vGrid subtree in abstract syntax

    ! Local parameters
    character (len=*), parameter :: UnitsMessage= &
         & "Inappropriate units for vertical coordinates"

    ! Local variables
    integer :: CoordType           ! One of t_vGridCoord's literals
    integer :: FIELD               ! First son of Son
    integer :: FIELD_INDEX         ! F_..., see Init_Tables_Module
    logical :: GOT_FIELD(field_first:field_last)
    integer :: I, J, K, L, N       ! Loop counter or limit
    integer :: NUMBER              ! Index in tree of Number field
    integer :: PER_DECADE          ! Index in tree of Per_decade field
    integer :: PREV_UNITS          ! Units of previous element of Value field
    integer :: SON                 ! Son of Root
    integer :: START               ! Index in tree of start field
    integer :: STATUS              ! From Allocate
    double precision :: STEP       ! Step for linear grid
    double precision :: STOP(2)    ! Value of Stop field
    integer :: STOP_UNITS(2)       ! Units of Stop field
    integer :: UNITS(2)            ! Output from Expr
    integer :: VALUE               ! Index in tree of value tree for a field
    integer :: VALUE_FIELD         ! Index in tree of value field
    double precision :: VALUES(2)  ! Output from Expr

    ! Executable code

    if ( toggle(gen) ) call trace_begin ( "CreateVGridFromMLSCFInfo", root )

    coordType = 0
    error = 0
    got_field = .false.
    number = 0
    vGrid%name = name
    vGrid%noSurfs = 0
    vGrid%verticalCoordinate = VC_Invalid

    do i = 2, nsons(root)
      son = subtree(i,root)
      field = subtree(1,son)
      value = subtree(2,son)
      field_index = decoration(field)
      if ( got_field(field_index) ) &
        & call announce_error ( field, duplicateField )
      got_field(field_index) = .true.
      select case ( field_index )
      case ( f_coordinate )
        vGrid%verticalCoordinate = decoration(value)
      case ( f_number )
        number = son
      case ( f_per_decade )
        per_decade = son
      case ( f_start )
        start = son
      case ( f_stop )
        call expr ( value, stop_units, stop )
      case ( f_type )
        coordType = decoration(value)
      case ( f_values )
        value_field = son
      case default ! Can't get here if tree_checker works properly
      end select
    end do

    ! Now check that this is a sensible vGrid; first the obvious stuff.

    if ( .not. got_field(f_coordinate) ) &
      & call announce_error ( root, noCoordField )
    if ( .not. got_field(f_type) ) call announce_error ( root, noTypeField )

    select case ( coordType )
    case ( l_explicit )
      call check_fields ( root, l_explicit, got_field, &
        & required=(/ f_values /), &
        & extra=(/ f_number, f_per_decade, f_start, f_stop /) )
      vgrid%noSurfs = nsons(value_field)-1
      allocate ( vgrid%surfs(vgrid%noSurfs), stat=status )
      if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//"vGrid%surfs" )
      if ( got_field(f_values) ) &
        & prev_units = check_units ( value_field, f_values, vgrid%surfs )
    case ( l_linear )
      call check_fields ( root, l_linear, got_field, &
        & required=(/ f_number, f_start, f_stop /), &
        & extra=(/ f_per_decade, f_values /) )
      if ( got_field(f_number) ) then
        call expr ( subtree(2,number), units, values )
        if ( units(1) /= phyq_dimensionless ) &
          & call announce_error ( subtree(1,number), unitless )
        vgrid%noSurfs = nint(values(1))
        if ( vgrid%noSurfs < 2 ) call announce_error ( root, tooFew )
        allocate ( vgrid%surfs(vgrid%noSurfs), stat=status )
        if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_Allocate//"vGrid%surfs" )
      end if
      if ( got_field(f_start) ) then
        call expr ( subtree(2,start), units, values )
        prev_units = units(1)
        vgrid%surfs(1) = values(1)
        if ( prev_units /= stop_units(1) ) &
          & call announce_error ( root, startStopUnits )
      end if
      if ( got_field(f_stop) ) vgrid%surfs(vgrid%noSurfs) = stop(1)
      if ( error == 0 ) then
        step = ( stop(1) - values(1) ) / ( number-1 )
        do i = 2, number-1
          vgrid%surfs(i) = vgrid%surfs(1) + (i-1) * step
        end do
      end if
    case ( l_logarithmic )
      call check_fields ( root, l_logarithmic , got_field, &
        & required=(/ f_number, f_per_decade, f_start /), &
        & extra=(/ f_stop, f_values /) )
      vgrid%noSurfs = 0
      if ( got_field(f_number) ) then
        j = nsons(number)
        do i = 2, j ! Compute total number of surfaces
          call expr ( subtree(i,number), units, values )
          if ( units(1) /= phyq_dimensionless ) &
            call announce_error ( subtree(1,number), wrongUnits )
          vgrid%noSurfs = vgrid%noSurfs + nint(values(1))
        end do
        if ( got_field(f_start) ) prev_units = check_units ( start, f_start )
        if ( error == 0 ) then
          if ( nsons(per_decade) /= j ) &
            & call announce_error ( root, inconsistentSizes )
        end if
      end if
      if ( error == 0 ) then
        allocate ( vgrid%surfs(vgrid%noSurfs), stat=status )
        if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
         & MLSMSG_Allocate//"vGrid%surfs" )
        k = 1
        n = 1 ! One less surface the first time, since we have one at the start.
        call expr ( subtree(2,start), units, values )
        vgrid%surfs(1) = values(1)
        do i = 2, j
          call expr ( subtree(i,per_decade), units, values )
          if ( values(1) <= 0.0d0 ) then
            call announce_error ( subtree(1,per_decade), notPositive )
            exit
          end if
          step = 10.0 ** (-1.0d0/values(1))
          if ( units(1) /= phyq_dimensionless ) then
            call announce_error ( subtree(1,per_decade), wrongUnits )
            exit
          end if
          call expr ( subtree(i,number), units, values )
          if ( units(1) /= phyq_dimensionless ) then
            call announce_error ( subtree(1,number), wrongUnits )
            exit
          end if
          do l = 1, nint(values(1)) - n
            k = k + 1
            vgrid%surfs(k) = vgrid%surfs(k-1) * step
          end do
          n = 0 ! Do all of the surfaces after the first time.
        end do
      end if
    end select

    ! Check that the given surfaces are in an appropriate unit

    if ( got_field(f_start) .and. &
      &( prev_units == PHYQ_Pressure .neqv. coordType == l_logarithmic) ) &
      & call announce_error ( subtree(1,start), unitsPressure )

    select case ( vGrid%verticalCoordinate )
    case (l_none)
      if ( coordType /= l_explicit ) &
        & call announce_error ( root, requireExplicit )
      if ( prev_units /= PHYQ_Dimensionless ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, UnitsMessage )
    case (l_pressure)
      if ( prev_units /= PHYQ_Pressure) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, UnitsMessage )
    case (l_zeta)
      if ( prev_units /= PHYQ_Dimensionless ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, UnitsMessage )
    case (l_geodAltitude)
      if ( prev_units /= PHYQ_Length) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, UnitsMessage )
    case (l_gph)
      if ( prev_units /= PHYQ_Length) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, UnitsMessage )
    case (l_theta)
      if ( prev_units /= PHYQ_Temperature) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, UnitsMessage )
    case (l_angle)
      if ( prev_units /= PHYQ_Angle) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, UnitsMessage )
    end select

    if ( toggle(gen) ) call trace_end ( "CreateVGridFromMLSCFInfo" )

  end function CreateVGridFromMLSCFInfo

  !----------------------------------------  DestroyVGridContents  -----
  subroutine DestroyVGridContents ( vGrid )

  ! This routine destroys the array information created with the vGrid

    ! Dummy arguments

    type (vGrid_T), intent(inout) :: vGrid

    ! Local Variables
    integer :: Status    ! From deallocate

    ! Executable code

    vGrid%noSurfs = 0
    vGrid%verticalCoordinate = VC_Invalid

    deallocate ( vGrid%surfs, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate // "vGrid%surfs" )

  end subroutine DestroyVGridContents

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

  ! ---------------------------------------  DestroyVGridDatabase  -----
  subroutine DestroyVGridDatabase ( DATABASE )

  ! This subroutine destroys a vGrid database

    ! Dummy argument
    type (VGrid_T), dimension(:), pointer :: DATABASE

    ! Local variables
    integer :: vgridIndex, Status

    if (associated(database)) then
      do vgridIndex = 1, SIZE(database)
        call DestroyVGridContents ( database(vgridIndex) )
      end do
      deallocate ( database, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Deallocate // "database" )
    end if
  end subroutine DestroyVGridDatabase

! =====     Private Procedures     =====================================

! -----------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( WHERE, CODE, FIELD_INDEX, LIT_INDEX )
    integer, intent(in) :: WHERE   ! Tree node where error was noticed
    integer, intent(in) :: CODE    ! Code for error message
    integer, intent(in), optional :: FIELD_INDEX, LIT_INDEX ! Extra stuff

    error = max(error,1)
    call output ( '***** At ' )
    call print_source ( source_ref(where) )
    call output ( ': ' )
    select case ( code )
    case ( duplicateField )
      call output ( "The " )
      call dump_tree_node ( where, 0 )
      call output ( " field appears more than once.", advance='yes' )
    case ( extraIf )
      call output ( "If the type is '" )
      call display_string ( lit_indices(lit_index) )
      call output ( "' the '" )
      call display_string ( field_indices(field_index) )
      call output ( "' field is not permitted.", advance='yes' )
    case ( inconsistentSizes )
      call output ( &
        & "The 'number' and 'per_decade' fields are", &
        & advance='yes' )
      call output ( "required to have the same sizes", advance='yes' )
    case ( inconsistentUnits )
      call output ( "Elements of the '" )
      call display_string ( field_indices(field_index) )
      call output ( "' field have inconsistent units", advance='yes' )
    case ( noCoordField )
      call output ( "The required field 'coordinate' is not present.", &
        & advance='yes' )
    case ( notPositive )
      call output ( "The value of the '" )
      call dump_tree_node ( where, 0 )
      call output ( "' field is required to be positive.", advance='yes' )
    case ( noTypeField )
      call output ( "The required field 'type' is not present.", advance='yes' )
    case ( requiredIf )
      call output ( "The '" )
      call display_string ( field_indices(field_index) )
      call output ( "' field is required if the type is '" )
      call display_string ( lit_indices(lit_index) )
      call output ( ".", advance='yes' )
    case ( requireExplicit )
      call output ( "If the 'coordinate' field is 'none', the 'type' field", &
        & advance='yes' )
      call output ( "shall be 'explicit'.", advance='yes' )
    case ( startStopUnits )
      call output ( &
        & "The 'start' and 'stop' fields do not have the same units.", &
        & advance='yes' )
    case ( tooFew )
      call output ( "If 'type' is linear 'number' >= 2 is required.", &
        & advance='yes' )
    case ( unitless )
      call output ( "The value for the " )
      call dump_tree_node ( where, 0 )
      call output ( " field is required to be unitless.", advance='yes' )
    case ( unitsPressure )
      call output ( "Pressure units are are required for logarithmic grids,", &
        & advance='yes' )
      call output ( "and prohibited otherwise.", advance='yes' )
    case ( wrongUnits )
      call output ( "The " )
      call dump_tree_node ( where, 0 )
      call output ( " field has incorrect units.", advance='yes' )
    end select
    end subroutine ANNOUNCE_ERROR

  ! -----------------------------------------------  CHECK_FIELDS  -----
  subroutine CHECK_FIELDS ( ROOT, GRID_TYPE, GOT_FIELD, REQUIRED, EXTRA )

  ! If any fields of Root that are listed in Extra are marked .true.
  ! in Got_field, or any fields of Root that are listed in Required are marked
  ! .false. in Got_field, announce an error.

    integer, intent(in) :: ROOT         ! Root of vGrid, used for message
    integer, intent(in) :: GRID_TYPE    ! Type field of vGrid, used for message
    logical, intent(in) :: Got_field(:) ! Which fields were present?
    integer, intent(in) :: REQUIRED(:)  ! List of required fields
    integer, intent(in) :: EXTRA(:)     ! List of prohibited fields

    integer :: I

    do i = 1, size(required)
      if ( .not. got_field(required(i)) ) &
        & call announce_error ( root, requiredIf, required(i), grid_type )
    end do
    do i = 1, size(extra)
      if ( got_field(extra(i)) ) &
        & call announce_error ( root, extraIf, extra(i), grid_type )
    end do
  end subroutine CHECK_FIELDS

  ! ------------------------------------------------  CHECK_UNITS  -----
  integer function CHECK_UNITS ( ROOT, FIELD_INDEX, FIELD_VALUES )

  ! Check that subtrees 2-n of Root have the same units, or that all
  ! but one are PHYQ_Dimensionless

    integer, intent(in) :: ROOT         ! Root of the subtree
    integer, intent(in) :: FIELD_INDEX  ! F_... From Init_Tables_Module
    real(r8), intent(out), optional :: FIELD_VALUES(:)  ! Values of the fields

    integer :: I                        ! Loop counter
    integer :: UNITS(2)                 ! Units of an expression
    double precision :: VALUES(2)       ! Values of an expression

    check_units = phyq_dimensionless
    do i = 2, nsons(root)
      call expr ( subtree(i,root), units, values )
      if ( present(field_values) ) field_values(i-1) = values(1)
      if ( check_units /= phyq_dimensionless .and. &
        &  units(1) /= phyq_dimensionless .and. &
        &  check_units /= units(1) ) &
        &  call announce_error ( subtree(1,field_index), inconsistentUnits, &
        &  field_index )
      if ( units(1) /= phyq_dimensionless ) check_units = units(1)
    end do
  end function CHECK_UNITS

end module vGrid

!
! $Log$
! Revision 2.1  2000/12/04 23:34:38  vsnyder
! Move more of addItemToDatabase into the include.
!
! Revision 2.0  2000/09/05 18:57:05  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry
!

