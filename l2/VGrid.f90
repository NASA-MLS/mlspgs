! Copyright (c) 1999, California Institute of Technology. ALL RIGHTS RESERVED.
! U.S. Government sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module vGrid                    ! Definitions for vGrids in vector quantities
!=============================================================================

  use Allocate_Deallocate, only: Allocate_Test
  use EXPR_M, only: EXPR
  use INIT_TABLES_MODULE, only: F_COORDINATE, F_FORMULA, F_NUMBER, F_START, &
    & F_STOP, F_TYPE, F_VALUES, FIELD_FIRST, FIELD_INDICES, FIELD_LAST, &
    & L_ANGLE, L_EXPLICIT, L_GEODALTITUDE, L_GPH, L_LINEAR, L_LOGARITHMIC, &
    & L_NONE, L_PRESSURE, L_THETA, L_ZETA, LIT_INDICES, PHYQ_Angle, &
    & PHYQ_Dimensionless, PHYQ_Length, PHYQ_Pressure, PHYQ_Temperature
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: R8       ! General constants etc.
  use MLSMessageModule, only: & ! Message logging
    & MLSMessage, MLSMSG_Allocate, MLSMSG_Error
  use OUTPUT_M, only: OUTPUT
  use STRING_TABLE, only: DISPLAY_STRING
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TOGGLES, only: GEN, TOGGLE
  use TREE, only: DECORATION, DUMP_TREE_NODE, NSONS, SOURCE_REF, SUBTREE
  use VGridsDatabase, only: AddVGridToDatabase, Dump, VGrid_T

  implicit none
  private

  public :: CreateVGridFromMLSCFInfo, Dump

  interface Dump
    module procedure MyDump_VGrids
  end interface Dump

  !------------------------------- RCS Ident Info ---------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), parameter, private :: ModuleName="$RCSfile$"
  !--------------------------------------------------------------------------

! -----     Private declarations     ---------------------------------

  integer, private :: ERROR

! Error codes for "announce_error"
  integer, private, parameter :: ExtraIf = 1
  integer, private, parameter :: InconsistentUnits = ExtraIf + 1
  integer, private, parameter :: NotPositive = InconsistentUnits + 1
  integer, private, parameter :: RequiredIf = NotPositive + 1
  integer, private, parameter :: RequireExplicit = RequiredIf + 1
  integer, private, parameter :: StartStopUnits = RequireExplicit + 1
  integer, private, parameter :: TooFew = StartStopUnits + 1
  integer, private, parameter :: Unitless = TooFew + 1
  integer, private, parameter :: UnitsPressure = Unitless + 1
  integer, private, parameter :: WrongUnits = UnitsPressure + 1

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
    integer :: FORMULA             ! Index in tree of formula field
    logical :: GOT_FIELD(field_first:field_last)
    integer :: I, J, K, L, N       ! Loop counter or limit
    integer :: NUMBER              ! Index in tree of Number field
    integer :: PREV_UNITS          ! Units of previous element of Value field
    integer :: SON                 ! Son of Root
    integer :: START               ! Index in tree of start field
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
    vGrid%verticalCoordinate = L_None

    do i = 2, nsons(root)
      son = subtree(i,root)
      field = subtree(1,son)
      value = subtree(2,son)
      field_index = decoration(field)
      got_field(field_index) = .true.
      select case ( field_index )
      case ( f_coordinate )
        vGrid%verticalCoordinate = decoration(value)
      case ( f_formula )
        formula = son
      case ( f_number )
        number = son
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

    select case ( coordType )
    case ( l_explicit )
      call check_fields ( root, l_explicit, got_field, &
        & required=(/ f_values /), &
        & extra=(/ f_formula, f_start, f_stop /) )
      vgrid%noSurfs = nsons(value_field)-1
      call allocate_test ( vgrid%surfs, vgrid%noSurfs, "vGrid%surfs", &
        & ModuleName )
      if ( got_field(f_values) ) &
        & prev_units = check_units ( value_field, f_values, vgrid%surfs )
    case ( l_linear )
      call check_fields ( root, l_linear, got_field, &
        & required=(/ f_number, f_start, f_stop /), &
        & extra=(/ f_formula, f_values /) )
      if ( got_field(f_number) ) then
        call expr ( subtree(2,number), units, values )
        if ( units(1) /= phyq_dimensionless ) &
          & call announce_error ( subtree(1,number), unitless )
        vgrid%noSurfs = nint(values(1))
        if ( vgrid%noSurfs < 2 ) call announce_error ( root, tooFew )
        call allocate_test ( vgrid%surfs, vgrid%noSurfs, "vGrid%surfs", &
          & ModuleName )
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
        & required=(/ f_formula, f_start /), &
        & extra=(/ f_stop, f_values /) )
      vgrid%noSurfs = 0
      if ( got_field(f_formula) ) then
        j = nsons(formula)
        do i = 2, j ! Compute total number of surfaces
          call expr ( subtree(i,formula), units, values )
          if ( units(1) /= phyq_dimensionless .or. &
            &  units(2) /= phyq_dimensionless ) &
            call announce_error ( subtree(1,formula), wrongUnits )
          vgrid%noSurfs = vgrid%noSurfs + nint(values(1))
        end do
        if ( got_field(f_start) ) prev_units = check_units ( start, f_start )
      end if
      if ( error == 0 ) then
        call allocate_test ( vgrid%surfs, vgrid%noSurfs, "vGrid%surfs", &
          & ModuleName )
        k = 1
        n = 1 ! One less surface the first time, since we have one at the start.
        call expr ( subtree(2,start), units, values )
        vgrid%surfs(1) = values(1)
        do i = 2, j
          call expr ( subtree(i,formula), units, values )
          if ( values(2) <= 0.0d0 ) then
            call announce_error ( subtree(1,formula), notPositive )
            exit
          end if
          step = 10.0 ** (-1.0d0/values(2))
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
      select case ( prev_units )
      case ( PHYQ_Dimensionless )     ! OK, do nothing
      case ( PHYQ_Pressure )          ! Need to take log
        vgrid%surfs= -LOG10(vgrid%surfs)
      case default
        call MLSMessage ( MLSMSG_Error, ModuleName, UnitsMessage )
      end select
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
    case ( extraIf )
      call output ( "If the type is '" )
      call display_string ( lit_indices(lit_index) )
      call output ( "' the '" )
      call display_string ( field_indices(field_index) )
      call output ( "' field is not permitted.", advance='yes' )
    case ( inconsistentUnits )
      call output ( "Elements of the '" )
      call display_string ( field_indices(field_index) )
      call output ( "' field have inconsistent units", advance='yes' )
    case ( notPositive )
      call output ( "The value of the '" )
      call dump_tree_node ( where, 0 )
      call output ( "' field is required to be positive.", advance='yes' )
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

  ! ----------------------------------------------  MyDump_VGrids  -----
  subroutine MyDump_VGrids ( VGrids )
    type(vGrid_T), intent(in), dimension(:) :: VGrids
    call dump ( vGrids, lit_indices )
  end subroutine MyDump_VGrids

end module vGrid

!
! $Log$
! Revision 2.8  2001/04/07 01:50:49  vsnyder
! Move some of VGrid to lib/VGridsDatabase.  Move ForwardModelConfig_T and
! some related stuff to fwdmdl/ForwardModelConfig.
!
! Revision 2.7  2001/03/28 03:03:38  vsnyder
! Remove use, only's that aren't used
!
! Revision 2.6  2001/03/28 01:25:38  vsnyder
! Move DUMP_VGRIDS from dumper.f90 to VGrid.f90
!
! Revision 2.5  2001/02/28 17:21:05  livesey
! Allowed user to specify zeta grids in pressure as well as log pressure space.
!
! Revision 2.4  2001/02/22 23:44:50  livesey
! Got rid of VC_Invalid
!
! Revision 2.3  2001/02/22 23:43:43  livesey
! Nullified vGrid_T%surfs by default
!
! Revision 2.2  2001/02/09 19:30:16  vsnyder
! Move checking for required and duplicate fields to init_tables_module
!
! Revision 2.1  2000/12/04 23:34:38  vsnyder
! Move more of addItemToDatabase into the include.
!
! Revision 2.0  2000/09/05 18:57:05  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry
!

