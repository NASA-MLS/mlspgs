! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module SpectroscopyCatalog_m

! Process the Spectroscopy section.  Read the "old format" spectroscopy catalog

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use String_Table, only: Display_String
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error
  use MLSSignals_m, only: MaxSigLen, Signals
  use Output_m, only: Blanks, Output
  use Toggles, only: Gen, Levels, Switches, Toggle
  use Trace_M, only: Trace_begin, Trace_end

  ! More USEs below in each procedure, if they're only used therein.

  private
  ! Public procedures:
  public :: Spectroscopy
  public :: Destroy_Line_Database, Destroy_SpectCat_Database
  public :: Dump_Lines_Database, Dump_SpectCat_Database

  ! Public types:
  type, public :: Line_T           ! One line in the spectrum for a species
    integer :: Line_Name           ! Sub_rosa index
    Real(r8) :: DELTA              ! Delta interference coefficient at 300K 1/mb
    Real(r8) :: EL                 ! Lower state energy cm-1
    Real(r8) :: GAMMA              ! Gamma interference coefficient at 300K 1/mb
    Real(r8) :: N                  ! Temperature power dependence of w
    Real(r8) :: N1                 ! Temperature dependency of delta
    Real(r8) :: N2                 ! Temperature dependency of gamma
    Real(r8) :: PS                 ! Pressure shift parameter in MHz/mbar
    Real(r8) :: STR                ! Integrated spectral intensity
                                   ! Log(nm**2 MHz) at 300 K
    Real(r8) :: V0                 ! Line center frequency MHz
    Real(r8) :: W                  ! Collision broadening parameter
                                   ! MHz/mbar at 300 K
  end type Line_T

  type, public :: Catalog_T        ! Catalog entry for a species
    integer :: Species_Name        ! Sub_rosa index
    integer :: Spec_Tag            ! Spectroscopy tag
    integer, pointer :: Lines(:)   ! Indices in Lines database
    integer :: Molecule            ! L_...
    real(r8) :: Qlog(3)            ! Logarithm of the partition function
                                   ! At 300 , 225 , and 150 K
  end type Catalog_T

  ! Public Variables:
  ! The lines database:
  type(line_T), public, pointer, dimension(:), save :: Lines => NULL()

  ! The spectroscopy database:
  type(catalog_T), public, pointer, dimension(:), save :: Catalog => NULL()

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter, private :: ModuleName = &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains ! =====  Public Procedures  ===================================

  ! -----------------------------------------------  Spectroscopy  -----
  subroutine Spectroscopy ( Root )
  ! Process the spectroscopy section.
    ! We need a lot of names from Init_Spectroscopy_Module.  First, the spec
    ! ID's:
    use Init_Spectroscopy_M, only: S_Line, S_Spectra
    ! Now the Fields:
    use Init_Spectroscopy_M, only: F_Delta, F_El, F_Gamma, F_Lines, &
      & F_Molecule, F_N, F_N1, F_N2, F_Ps, F_Qlog, F_Str, F_V0, F_W
    use Intrinsic, only: Phyq_Dimless => Phyq_Dimensionless, Phyq_Frequency, &
      & S_Time
    use Molecules, only: Spec_Tags
    use MoreTree, only: Get_Field_Id, Get_Spec_Id
    use Tree, only: Decorate, Decoration, Node_ID, NSons, Sub_Rosa, Subtree
    use Tree_Types, only: N_Named

    ! Dummy argument
    integer, intent(in) :: Root         ! Of the AST for the section

    ! Local Variables
    integer :: Error                    ! /= 0 => An error occured
    integer :: I, J, K                  ! Loop inductors, Subscripts
    type(line_t) :: OneLine             ! To be added to the database
    type(catalog_t) :: OneSpecies       ! To be added to the database
    integer :: Son                      ! Of root or key
    integer :: Key                      ! Index of spec_arg
    integer :: Name                     ! Index of name in string table
    logical :: TIMING                   ! For S_Time
    real :: T1, T2                      ! For S_Time

    ! Error message codes
    integer, parameter :: WrongSize = 1                ! Wrong number of elements
    integer, parameter :: WrongUnits = wrongSize + 1   ! Wrong physical units

    if ( toggle(gen) ) call trace_begin ( "Spectroscopy", root )

    error = 0
    timing = .false.

    do i = 2, nsons(root)-1             ! Skip names of section
      son = subtree(i,root)
      if ( node_id(son) == n_named ) then
        name = sub_rosa(subtree(1,son))
        key = subtree(2,son)
      else
        name = 0
        key = son
      end if
      select case ( get_spec_id(key) )
      case ( s_line ) ! ...................................  LINE  .....
        oneLine%line_Name = name
        do j = 2, nsons(key)
          son = subtree(j,key)
          select case ( get_field_id(son) )
          case ( f_delta )
            call expr_check ( subtree(2,son), oneLine%delta, phyq_dimless )
          case ( f_el )
            call expr_check ( subtree(2,son), oneLine%el, phyq_dimless )
          case ( f_gamma )
            call expr_check ( subtree(2,son), oneLine%gamma, phyq_dimless )
          case ( f_n )
            call expr_check ( subtree(2,son), oneLine%n, phyq_dimless )
          case ( f_n1 )
            call expr_check ( subtree(2,son), oneLine%n1, phyq_dimless )
          case ( f_n2 )
            call expr_check ( subtree(2,son), oneLine%n2, phyq_dimless )
          case ( f_ps )
            call expr_check ( subtree(2,son), oneLine%ps, phyq_dimless )
          case ( f_str )
            call expr_check ( subtree(2,son), oneLine%str, phyq_dimless )
          case ( f_v0 )
            call expr_check ( subtree(2,son), oneLine%v0, phyq_frequency )
          case ( f_w )
            call expr_check ( subtree(2,son), oneLine%w, phyq_dimless )
          case default
            ! Can't get here if the type checker worked
          end select
        end do
        call decorate ( key, addLineToDatabase ( lines, oneLine ) )
      case ( s_spectra ) ! .............................  SPECTRA  .....
        nullify(oneSpecies%lines) ! So that Allocate_Test doesn't deallocate
        ! the most recently filled one
        oneSpecies%species_Name = name
        do j = 2, nsons(key)
          son = subtree(j,key)
          select case ( get_field_id(son) )
          case ( f_lines )
            k = nsons(son)
            call allocate_test ( oneSpecies%lines, k-1, "OneSpecies%Lines", &
              & moduleName )
            do k = 2, k
              oneSpecies%lines(k-1) = decoration(decoration(subtree(k,son)))
            end do
          case ( f_molecule )
            oneSpecies%molecule = decoration(subtree(2,son))
            oneSpecies%spec_Tag = spec_Tags(oneSpecies%molecule)
          case ( f_qlog )
            if ( nsons(son) /= 4 ) call announce_error ( son, wrongSize, 3 )
            do k = 2, 4
              call expr_check ( subtree(k,son), oneSpecies%qlog(k-1), &
                & phyq_dimless )
            end do
          end select
        end do
        call decorate ( key, addSpeciesToCatalog ( Catalog, oneSpecies ) )
      case ( s_time ) ! ...................................  TIME  .....
        if ( timing ) then
          call sayTime
        else
          call cpu_time ( t1 )
          timing = .true.
        end if
      end select
    end do ! i

    if ( toggle(gen) ) then
      if ( levels(gen) > 0 .or. index(switches,'C') /= 0 ) &
        & call dump_SpectCat_database
      call trace_end ( "Spectroscopy" )
    end if
    if ( timing ) call sayTime

    if ( error > 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & 'Error(s) in input for spectroscopy database' )
  contains
    ! ...........................................  Announce_Error  .....
    subroutine Announce_Error ( Where, Code, More )
      use Intrinsic, only: Phyq_Indices
      use Lexer_Core, only: Print_Source
      use Tree, only: Source_Ref
      integer, intent(in) :: Where      ! In the tree
      integer, intent(in) :: Code       ! The error code
      integer, intent(in), optional :: MOre  ! In case some error messages need

      error = max(error, 1)
      call output ( '***** At ' )
      call print_source ( source_ref ( where ) )
      call output ( ' Spectroscopy complained: ' )
      select case ( code )
      case ( wrongSize )
        call output ( 'The field is required to have ' )
        call output ( more )
        call output ( ' elements.', advance='yes' )
      case ( wrongUnits )
        call output ( "The field's units ought to be " )
        call display_string ( phyq_indices(more), advance='yes' )
      end select
    end subroutine Announce_Error

    ! ...............................................  Expr_Check  .....
    subroutine Expr_Check ( Root, Value, NeededUnits )
    ! Evaluate the expression at Root giving Value.  Make sure its units
    ! are NeededUnits
      use Expr_M, only: Expr
      integer, intent(in) :: Root
      real(r8), intent(out) :: Value
      integer, intent(in) :: NeededUnits
      integer :: Units(2)               ! From expr
      double precision :: Values(2)     ! From expr
      call expr ( root, units, values )
      if ( units(1) /= neededUnits ) &
        & call announce_error ( root, wrongUnits, neededUnits )
      value = values(1)
    end subroutine Expr_Check

    ! ..................................................  SayTime  .....
    subroutine SayTime
      call cpu_time ( t2 )
      call output ( "Timing for Spectroscopy = " )
      call output ( dble(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime
  end subroutine Spectroscopy

  ! ------------------------------------------  AddLineToDatabase  -----
  integer function AddLineToDatabase ( Database, Item )
  ! Add a line to the Lines database, creating the database
  ! if necessary.

    ! Dummy arguments
    type(line_T), pointer, dimension(:) :: Database
    type(line_T), intent(in) :: Item


    ! Local variables
    type(line_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddLineToDatabase = newSize
  end function AddLineToDatabase

  ! ----------------------------------------  AddSpeciesToCatalog  -----
  integer function AddSpeciesToCatalog ( Database, Item )
  ! Add a line to the Lines database, creating the database
  ! if necessary.

    ! Dummy arguments
    type(catalog_t), pointer, dimension(:) :: Database
    type(catalog_t), intent(in) :: Item


    ! Local variables
    type(catalog_t), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddSpeciesToCatalog = newSize
  end function AddSpeciesToCatalog

  ! --------------------------------------  Destroy_Line_Database  -----
  subroutine Destroy_Line_Database
    integer :: Status                   ! From deallocate
    deallocate ( lines, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_DeAllocate // "Lines" )
  end subroutine Destroy_Line_Database

  ! ----------------------------------  Destroy_SpectCat_Database  -----
  subroutine Destroy_SpectCat_Database
    integer :: I, Status
    do i = 1, size(catalog)
      call deallocate_test ( catalog(i)%lines, moduleName, "Catalog(i)%lines" )
    end do ! i
    deallocate ( catalog, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      MLSMSG_DeAllocate // "Catalog" )
  end subroutine Destroy_SpectCat_Database

  ! ----------------------------------------  Dump_Lines_Database  -----
  subroutine Dump_Lines_Database ( Start, End, Number )
    use Dump_0, only: Dump
    integer, intent(in), optional :: Start, End
    logical, intent(in), optional :: Number
    integer :: I                   ! Subscript, loop inductor
    integer :: MyStart, MyEnd
    logical :: MyNumber

    myStart = 1
    if ( present(start) ) myStart = start
    myEnd = size(lines)
    if ( present(end) ) myEnd = end
    myNumber = .true.
    if ( present(number) ) myNumber = number
    if ( .not. present(start) .and. .not. present(end) ) then
      call output ('Spectroscopy lines database: SIZE = ' )
      call output ( size(lines), advance='yes' )
    end if
    do i = myStart, myEnd
      if ( myNumber ) then
        call output ( i, 4 )
        call output ( ': ' )
      else
        call blanks ( 6 )
      end if
      if ( lines(i)%line_name /= 0 ) then
        call output ( 'Name = ' )
        call display_string ( lines(i)%line_name )
        call output ( ', ' )
      end if
      call output ( 'V0 = ' )
      call output ( lines(i)%v0 )
      call output ( ', El = ' )
      call output ( lines(i)%el )
      call output ( ', Str =' )
      call output ( lines(i)%str )
      call output ( ', W =' )
      call output ( lines(i)%str, advance='yes' )
      call blanks ( 6 )
      call output ( 'Ps = ' )
      call output ( lines(i)%ps )
      call output ( ', N = ' )
      call output ( lines(i)%n )
      call output ( ': Delta = ' )
      call output ( lines(i)%delta )
      call output ( ', N1 = ' )
      call output ( lines(i)%n1, advance='yes' )
      call blanks ( 6 )
      call output ( 'Gamma = ' )
      call output ( lines(i)%gamma )
      call output ( ', N2 = ' )
      call output ( lines(i)%n2, advance='yes' )
    end do
  end subroutine Dump_Lines_Database
  ! -------------------------------------  Dump_SpectCat_Database  -----
  subroutine Dump_SpectCat_Database
    use Dump_0, only: Dump
    use Intrinsic, only: Lit_indices

    integer :: I, J                ! Subscript, loop inductor

    call output ( 'Spectroscopy catalog: SIZE = ' )
    call output ( size(catalog), advance='yes' )
    do i = 1, size(catalog)
      call output ( i, 4 )
      call output ( ': ' )
      if ( catalog(i)%species_name /= 0 ) then
        call display_string ( catalog(i)%species_name )
        call output ( ', ' )
      end if
      call output ( 'Species = ' )
      call display_string ( lit_indices(catalog(i)%molecule) )
      call output ( ', SpecTag = ' )
      call output ( catalog(i)%spec_tag )
      call output ( ', Qlog = [ ' )
      do j = 1, 3
        call output ( catalog(i)%qlog(j) )
        if ( j < 3 ) call output ( ', ' )
      end do
      call output ( ' ]', advance='yes' )
      call blanks ( 6 + int(log10(i+0.0)) )
      call output ( 'Lines:', advance='yes' )
      do j = 1, size(catalog(i)%lines)
        call dump_lines_database ( catalog(i)%lines(j), catalog(i)%lines(j), &
          & .false. )
      end do
    end do ! i
  end subroutine Dump_SpectCat_Database

end module SpectroscopyCatalog_m

! $Log$
! Revision 1.6  2001/04/26 02:36:52  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 1.5  2001/04/23 23:16:16  vsnyder
! Add 'time' command
!
! Revision 1.4  2001/04/20 17:26:31  vsnyder
! Remove arguments from Destroy..., publish ...Lines...
!
! Revision 1.3  2001/04/04 23:56:46  zvi
! Correfting Typo error in SpectCat  names
!
! Revision 1.2  2001/04/04 23:21:46  vsnyder
! Add comments for fields of Lines_T and Catalog_T
!
! Revision 1.1  2001/04/04 02:09:16  vsnyder
! Initial Commit
!
