! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module SpectroscopyCatalog_m

! Process the Spectroscopy section.  Read the "old format" spectroscopy catalog

  use Intrinsic, only: L_none
  use MLSCommon, only: R8

  ! More USEs below in each procedure.

  implicit none

  private
  ! Public procedures:
  public :: Spectroscopy, SearchByQn
  public :: Destroy_Line_Database, Destroy_SpectCat_Database
  public :: Dump_Lines_Database, Dump_SpectCat_Database, Dump

  interface DUMP
    module procedure Dump_SpectCat_Database, Dump_SpectCat_Database_2d
  end interface

  ! Private parameters
  integer, parameter :: MaxContinuum = 6

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
    Real(r8) :: NS                 ! Pressure shift on temperature dependency
    Real(r8) :: STR                ! Integrated spectral intensity
                                   ! Log(nm**2 MHz) at 300 K
    Real(r8) :: V0                 ! Line center frequency MHz
    Real(r8) :: W                  ! Collision broadening parameter
                                   ! MHz/mbar at 300 K
    integer, dimension(:), pointer :: QN=>NULL()      ! Optional quantum numbers
    integer, dimension(:), pointer :: Signals=>NULL() ! List of signal indices for line
    integer, dimension(:), pointer :: Sidebands=>NULL() ! Sidebands for above bands (-1,0,1)
    logical, dimension(:), pointer :: Polarized=>NULL() ! Process this signal and
                                   ! sideband using the polarized model
  end type Line_T

  type, public :: Catalog_T        ! Catalog entry for a species
    real(r8) :: continuum(MaxContinuum)      ! Continuum coefficients
    integer, pointer :: Lines(:)=>NULL() ! Indices in Lines database
    real(r8) :: Mass               ! Molecular mass in AMU
    integer :: Molecule            ! L_...
    logical, pointer :: Polarized(:)=>NULL() ! Used only in catalog extract to
                                   ! indicate that the lines(:) are to be
                                   ! processed with the polarized model
    real(r8) :: Qlog(3)            ! Logarithm of the partition function
                                   ! At 300 , 225 , and 150 K
    integer :: Species_Name        ! Sub_rosa index
  end type Catalog_T

  type(catalog_T), public, parameter :: Empty_Cat = catalog_t ( &
    & 0.0_r8, NULL(), 0.0_r8, l_none, NULL(), 0.0_r8, 0 )

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
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains ! =====  Public Procedures  ===================================

  ! -----------------------------------------------  Spectroscopy  -----
  subroutine Spectroscopy ( Root )
  ! Process the spectroscopy section.
    ! We need a lot of names from Init_Spectroscopy_Module.  First, the spec
    ! ID's:

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use Init_Spectroscopy_M, only: S_Line, S_Spectra
    ! Now the Fields:
    use Init_Spectroscopy_M, only: F_Continuum, F_Delta, F_El, F_Gamma, F_Lines, &
      & F_Mass, F_Molecule, F_N, F_N1, F_N2, F_Ns, F_Ps, F_Qlog, F_QN, F_Str, F_V0, F_W, &
      & F_EMLSSIGNALS, F_EMLSSIGNALSPOL, F_MLS1SIGNALS, F_UMLSSIGNALS
    use Intrinsic, only: Phyq_Dimless => Phyq_Dimensionless, Phyq_Frequency, &
      & S_Time, L_EMLS, L_UMLS, L_MLS1
    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, &
      & MLSMSG_DeAllocate, MLSMSG_Error
    use MoreTree, only: Get_Field_Id, Get_Spec_Id
    use Parse_Signal_m, only: PARSE_SIGNAL
    use String_Table, only: Get_string
    use Time_M, only: Time_Now
    use Toggles, only: Gen, Switches, Toggle
    use Trace_M, only: Trace_begin, Trace_end
    use Tree, only: Decorate, Decoration, Node_ID, NSons, Sub_Rosa, Subtree
    use Tree_Types, only: N_Named
    use MLSSignals_m, only: Instrument

    ! Dummy argument
    integer, intent(in) :: Root         ! Of the AST for the section

    ! Local Variables
    integer :: Error                    ! /= 0 => An error occured
    logical :: GotLines, GotMass        ! Got a "lines" or "mass" field
    integer :: I, J, K, L               ! Loop inductors, Subscripts
    integer :: Key                      ! Index of spec_arg
    integer :: Name                     ! Index of name in string table
    integer :: NOSIGNALS                ! For the bands part
    integer :: NumCatalog               ! Number of species in catalog
    integer :: NumLines                 ! Number of lines in catalog
    integer :: OffsetCatalog            ! Number of species previously in catalog
    integer :: OffsetLines              ! Number of lines previously in catalog
    real(r8) :: QN                      ! for call to expr_check for QN field
    integer :: SIDEBAND                 ! A single sideband
    integer :: SignalsNode              ! Tree node for emls/umls bands
    integer :: SignalsNodePol           ! Tree node for emls/umls bands for
                                        ! Zeeman-split lines
    integer, dimension(:), pointer :: SIGINDS ! From Parse_signal
    character ( len=80 ) :: SIGNAME     ! The signal
    integer :: Son                      ! Of root or key
    integer :: Status                   ! From Allocate or Deallocate
    type(line_t), pointer, dimension(:) :: TempLines
    type(catalog_t), pointer, dimension(:) :: TempCatalog
    integer :: TheSignal                ! SubRosa for a signal
    integer :: THISMANY                 ! Conted up to noSignals
    logical :: TIMING                   ! For S_Time
    real :: T1, T2                      ! For S_Time

    ! Error message codes
    integer, parameter :: No_Mass = 1   ! Lines field but no mass field
    integer, parameter :: NotInt = no_mass + 1         ! QN not an integer
    integer, parameter :: NotListedSignal = NotInt + 1 ! Polarized signal is not
                                        ! listed as emlsSignal
    integer, parameter :: TooBig = NotListedSignal + 1 ! Too many elements
    integer, parameter :: WrongSize = TooBig + 1       ! Wrong number of elements
    integer, parameter :: WrongUnits = wrongSize + 1   ! Wrong physical units

    if ( toggle(gen) ) call trace_begin ( "Spectroscopy", root )

    error = 0
    timing = .false.

    ! Determine sizes of created or expanded lines and species databases.
    offsetLines = 0
    if ( associated(lines) ) offsetLines = size(lines)
    numLines = offsetLines
    offsetCatalog = 0
    if ( associated(catalog) ) offsetCatalog = size(catalog)
    numCatalog = offsetCatalog
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
        numLines = numLines + 1
      case ( s_spectra ) ! .............................  SPECTRA  .....
        numCatalog = numCatalog + 1
      case default ! Don't care about the others
      end select
    end do

    ! Create or expand the Lines database
    tempLines => Lines
    allocate ( lines(numLines), stat=status )
    if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "Lines")
    if ( associated(tempLines) ) then
      lines(:offsetLines) = tempLines
      deallocate ( tempLines, stat=status )
      if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate // "TempLines")
    end if

    ! Create or expand the Species Catalog
    tempCatalog => Catalog
    allocate ( catalog(numCatalog), stat=status )
    if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "Catalog")
    if ( associated(tempCatalog) ) then
      catalog(:offsetCatalog) = tempCatalog
      deallocate ( tempCatalog, stat=status )
      if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate // "TempCatalog")
    end if

    numLines = offsetLines
    numCatalog = offsetCatalog
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
        numLines = numLines + 1
        lines(numLines)%line_Name = name
        signalsNode = 0
        signalsNodePol = 0
        do j = 2, nsons(key)
          son = subtree(j,key)
          select case ( get_field_id(son) )
          case ( f_delta )
            call expr_check ( subtree(2,son), lines(numLines)%delta, phyq_dimless )
          case ( f_el )
            call expr_check ( subtree(2,son), lines(numLines)%el, phyq_dimless )
          case ( f_emlsSignals )
            if ( instrument == l_emls ) signalsNode = son
          case ( f_emlsSignalsPol )
            if ( instrument == l_emls ) signalsNodePol = son
          case ( f_gamma )
            call expr_check ( subtree(2,son), lines(numLines)%gamma, phyq_dimless )
          case ( f_mls1Signals )
            if ( instrument == l_mls1 ) signalsNode = son
          case ( f_n )
            call expr_check ( subtree(2,son), lines(numLines)%n, phyq_dimless )
          case ( f_n1 )
            call expr_check ( subtree(2,son), lines(numLines)%n1, phyq_dimless )
          case ( f_n2 )
            call expr_check ( subtree(2,son), lines(numLines)%n2, phyq_dimless )
          case ( f_ns )
            call expr_check ( subtree(2,son), lines(numLines)%ns, phyq_dimless )
          case ( f_ps )
            call expr_check ( subtree(2,son), lines(numLines)%ps, phyq_dimless )
          case ( f_qn )
            call allocate_test ( lines(numLines)%qn, nsons(son)-1, 'qn', ModuleName )
            do k = 2, nsons(son)
              call expr_check ( subtree(k,son), qn, phyq_dimless )
              lines(numLines)%qn(k-1) = nint(qn)
              if ( abs(qn - lines(numLines)%qn(k-1)) > 0.1_r8 ) &
                & call announce_error ( subtree(k,son), notInt )
            end do
          case ( f_str )
            call expr_check ( subtree(2,son), lines(numLines)%str, phyq_dimless )
          case ( f_umlsSignals )
            if ( instrument == l_umls ) signalsNode = son
          case ( f_v0 )
            call expr_check ( subtree(2,son), lines(numLines)%v0, phyq_frequency )
          case ( f_w )
            call expr_check ( subtree(2,son), lines(numLines)%w, phyq_dimless )
          case default
            ! Can't get here if the type checker worked
          end select
        end do
        if ( signalsNode /= 0 ) then
          ! First work out how many signals we're dealing with
          noSignals = 0
          nullify ( sigInds )
          do j = 2, nsons(signalsNode)    ! Skip name
            call get_string ( sub_rosa ( subtree (j, signalsNode) ), &
              & sigName, strip=.true. )
            call Parse_Signal ( sigName, sigInds, signalsNode, onlyCountEm=thisMany )
            noSignals = noSignals + thisMany
          end do

          ! Compile list of all the bands/sidebands named
          call Allocate_test ( lines(numLines)%signals, noSignals, 'signals', ModuleName )
          call Allocate_test ( lines(numLines)%sidebands, noSignals, 'sidebands', ModuleName )
          if ( signalsNodePol /= 0 ) then
            call allocate_test ( lines(numLines)%polarized, noSignals, 'Polarized', moduleName )
            lines(numLines)%polarized = .false.
          end if
          nullify ( sigInds )
          k = 1
          do j = 2, nsons(signalsNode)    ! Skip name
            theSignal = sub_rosa ( subtree (j, signalsNode) )
            call get_string ( theSignal, sigName, strip=.true. )
            call Parse_Signal ( sigName, sigInds, signalsNode, sideband=sideband )
            if ( .not. associated(sigInds) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Invalid signal in spectroscopy' )
            
            lines(numLines)%signals ( k:k+size(sigInds)-1 ) = sigInds
            lines(numLines)%sidebands ( k:k+size(sigInds)-1 ) = sideband
            if ( signalsNodePol /= 0 ) then
              do l = 2, nsons(signalsNodePol)
                if ( sub_rosa ( subtree (l, signalsNodePol) ) == theSignal ) then
                  lines(numLines)%polarized ( k:k+size(sigInds)-1 ) = .true.
              exit
                end if
              end do
            end if
            k = k + size(sigInds)
          end do
          call Deallocate_test ( sigInds, 'sigInds', ModuleName )
        end if
        if ( signalsNodePol /= 0 ) then
          if ( signalsNode == 0 ) then
            call announce_Error ( signalsNodePol, notListedSignal )
          else
            k = nsons(signalsNode)
            do j = 2, nsons(signalsNodePol)
              theSignal = sub_rosa ( subtree (j, signalsNodePol) )
              do l = 2, k
                if ( theSignal == sub_rosa ( subtree (l, signalsNode) ) ) exit
              end do
              if ( l > k ) then
                call announce_Error ( subtree (j, signalsNodePol), notListedSignal )
            cycle
              end if
            end do
          end if
        end if
        call decorate ( key, numLines )
      case ( s_spectra ) ! .............................  SPECTRA  .....
        numCatalog = numCatalog + 1
        catalog(numCatalog)%species_Name = name
        catalog(numCatalog)%continuum = 0.0
        gotLines = .false.; gotMass = .false.
        do j = 2, nsons(key)
          son = subtree(j,key)
          select case ( get_field_id(son) )
          case ( f_lines )
            gotLines = .true.
            call allocate_test ( catalog(numCatalog)%lines, nsons(son)-1, &
              & "catalog(numCatalog)%Lines", moduleName )
            do k = 2, nsons(son)
              catalog(numCatalog)%lines(k-1) = decoration(decoration(subtree(k,son)))
            end do
          case ( f_continuum )
            if ( nsons(son) > MaxContinuum + 1 ) &
              & call announce_error ( son, tooBig, MaxContinuum )
            do k = 2, nsons(son)
              call expr_check ( subtree(k,son), catalog(numCatalog)%continuum(k-1), &
                & phyq_dimless )
            end do
          case ( f_mass )
            gotMass = .true.
            if ( nsons(son) /= 2 ) call announce_error ( son, wrongSize, 1 )
            call expr_check ( subtree(2,son), catalog(numCatalog)%mass, &
              & phyq_dimless )
          case ( f_molecule )
            catalog(numCatalog)%molecule = decoration(subtree(2,son))
          case ( f_qlog )
            if ( nsons(son) /= 4 ) call announce_error ( son, wrongSize, 3 )
            do k = 2, 4
              call expr_check ( subtree(k,son), catalog(numCatalog)%qlog(k-1), &
                & phyq_dimless )
            end do
          end select
        end do
        if ( gotLines .and. .not. gotMass ) call announce_error ( key, no_mass )
        call decorate ( key, numCatalog )
      case ( s_time ) ! ...................................  TIME  .....
        if ( timing ) then
          call sayTime
        else
          call time_now ( t1 )
          timing = .true.
        end if
      end select
    end do ! i

    if ( index(switches,'spec') /= 0 ) call dump_SpectCat_database ( catalog )
    if ( toggle(gen) ) then
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
      use Output_m, only: Output
      use String_Table, only: Display_String
      use Tree, only: Source_Ref
      integer, intent(in) :: Where      ! In the tree
      integer, intent(in) :: Code       ! The error code
      integer, intent(in), optional :: MOre  ! In case some error messages need

      error = max(error, 1)
      call output ( '***** At ' )
      call print_source ( source_ref ( where ) )
      call output ( ' Spectroscopy complained: ' )
      select case ( code )
      case ( no_mass )
        call output ( &
          & 'A "mass" field is required if the "lines" field is present.', &
          & advance='yes' )
      case ( notInt )
        call output ( 'The field is too far from being an integer.',& 
          & advance='yes' )
      case ( notListedSignal )
        call output ( 'The Polarized signal is not listed as an EMLS signal.', &
          & advance='yes' )
      case ( tooBig )
        call output ( 'The field cannot have more than ' )
        call output ( more )
        call output ( ' elements.', advance='yes' )
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
      use Output_m, only: Output
      call time_now ( t2 )
      call output ( "Timing for Spectroscopy = " )
      call output ( dble(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime
  end subroutine Spectroscopy

  ! ------------------------------------------  AddLineToDatabase  -----
  integer function AddLineToDatabase ( Database, Item )
  ! Add a line to the Lines database, creating the database
  ! if necessary.

    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
      & MLSMSG_Error

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

    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
      & MLSMSG_Error

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

    use MLSMessageModule, only: MLSMessage, MLSMSG_DeAllocate, MLSMSG_Error

    integer :: Status                   ! From deallocate
    if ( .not. associated(lines) ) return
    deallocate ( lines, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_DeAllocate // "Lines" )
  end subroutine Destroy_Line_Database

  ! ----------------------------------  Destroy_SpectCat_Database  -----
  subroutine Destroy_SpectCat_Database
    use Allocate_Deallocate, only: Deallocate_Test
    use MLSMessageModule, only: MLSMessage, MLSMSG_DeAllocate, MLSMSG_Error
    integer :: I, Status
    if ( .not. associated(catalog) ) return
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
    use Output_m, only: Blanks, Output
    use String_Table, only: Display_String
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
      call output ( 'Ns = ' )
      call output ( lines(i)%ns )
      call output ( ', Ps = ' )
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
      if ( associated(lines(i)%qn) ) call dump ( lines(i)%qn, '      QN' )
      if ( associated(lines(i)%signals) ) &
        & call dump ( lines(i)%signals, '      Signals' )
      if ( associated(lines(i)%sidebands) ) &
        & call dump ( lines(i)%sidebands, '      Sidebands' )
      if ( associated(lines(i)%polarized) ) &
        & call dump ( lines(i)%polarized, '      Polarized' )
    end do
  end subroutine Dump_Lines_Database
  ! ----------------------------------  Dump_SpectCat_Database_2D  -----
  subroutine Dump_SpectCat_Database_2d ( Catalog, Name, Details )
    type(catalog_T), pointer :: Catalog(:,:)
    character(len=*), intent(in), optional :: Name
    integer, optional, intent(in) :: Details ! <= 0 => Don't dump lines, default 0
    integer :: Sideband
    ! Executable code
    do sideband = lbound(catalog,1), ubound(catalog,1), 2
      call Dump ( catalog(sideband,:), name, sideband, details )
    end do
  end subroutine Dump_SpectCat_Database_2d

  ! -------------------------------------  Dump_SpectCat_Database  -----
  subroutine Dump_SpectCat_Database ( Catalog, Name, Sideband, Details )
    use Dump_0, only: Dump
    use Intrinsic, only: Lit_indices
    use Output_m, only: Blanks, NewLine, Output
    use String_Table, only: Display_String

    type(catalog_T), intent(in) :: Catalog(:)
    character(len=*), intent(in), optional :: Name
    integer, intent(in), optional :: Sideband
    integer, optional, intent(in) :: Details ! <= 0 => Don't dump lines, default 0

    integer :: I, J                ! Subscript, loop inductor
    integer :: MyDetails
    character(len=15) :: Print

    myDetails = 0
    if ( present(details) ) myDetails = details
    call output ( 'Spectroscopy catalog' )
    if ( present(name) ) call output ( ' '//trim(name) )
    if ( present(sideband) ) call output ( sideband, before=' for sideband ' )
    call output ( ': SIZE = ' )
    call output ( size(catalog), advance='yes' )
    do i = 1, size(catalog)
      if ( catalog(i)%molecule == l_none ) then
        if ( .not. associated(catalog(i)%lines) ) cycle
        if ( size(catalog(i)%lines) == 0 ) cycle
      end if
      call output ( i, 4 )
      call output ( ': ' )
      if ( catalog(i)%species_name /= 0 ) then
        call display_string ( catalog(i)%species_name )
        call output ( ', ' )
      end if
      call output ( 'Species = ' )
      call display_string ( lit_indices(catalog(i)%molecule) )
      write ( print, '(f10.3)' ) catalog(i)%mass
      call output ( ' Mass = ' // trim(adjustl(print)) // ', Qlog = [ ' )
      do j = 1, 3
        write ( print, '(f10.4)' ) catalog(i)%qlog(j)
        call output ( trim(adjustl(print)) )
        if ( j < 3 ) call output ( ', ' )
      end do
      call output ( ' ]', advance='yes' )
      call blanks ( 6 )
      call output ( 'continuum = [ ' )
      do j = 1, MaxContinuum
        write ( print, '(g15.3)' ) catalog(i)%continuum(j)
        call output ( trim(adjustl(print)) )
        if ( j < MaxContinuum ) call output ( ', ' )
      end do
      call output ( ' ]', advance='yes' )
      if ( associated(catalog(i)%polarized) ) &
        & call dump ( catalog(i)%polarized, '      Polarized:' )
      if ( myDetails > 0 ) then
        call blanks ( 6 )
        call output ( 'Lines:' )
        if ( associated(catalog(i)%lines) ) then
          if ( size(catalog(i)%lines) > 0 ) then
            call newLine
            do j = 1, size(catalog(i)%lines)
              call dump_lines_database ( catalog(i)%lines(j), catalog(i)%lines(j), &
                & .false. )
            end do
          else
            call output ( ' none', advance='yes' )
          end if
        else
          call output ( ' none', advance='yes' )
        end if
      end if
    end do ! i
  end subroutine Dump_SpectCat_Database

! ---------------------------------------------------  SearchByQn  -----
  integer function SearchByQn ( molecule, QNs )
  ! Find an element in the Catalog database that has a molecule.
  ! Then in its Lines database find an object that has its QN field
  ! equal to QNs.  If found, return the index of that line in the lines
  ! database.  Return zero if no such object can be found.

    integer, intent(in) :: MOLECULE
    integer, intent(in) :: QNs(:)

    integer :: I, J

    do i = 1, size(catalog)
      if ( catalog(i)%molecule == molecule ) then
        if ( associated(catalog(i)%lines) ) then
          do j = 1, size(catalog(i)%lines)
            if ( associated(lines(catalog(i)%lines(j))%qn) ) then
              if ( size(lines(catalog(i)%lines(j))%qn) == size(qns) ) then
                if ( all(lines(catalog(i)%lines(j))%qn == qns) ) then
                  searchByQn = catalog(i)%lines(j)
                  return
                end if
              end if
            end if
          end do
        end if
      end if
    end do
    searchByQn = 0 ! not found

  end function SearchByQn

! ------------------------------------------------  not_used_here  -----
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module SpectroscopyCatalog_m

! $Log$
! Revision 2.25  2004/08/07 00:34:57  vsnyder
! Correct subscript error in Dump_SpectCat_Database_2d
!
! Revision 2.24  2004/08/03 02:27:05  vsnyder
! Add 'Details' argument to dump
!
! Revision 2.23  2004/07/08 02:47:44  vsnyder
! Fix up the dump routines
!
! Revision 2.22  2004/04/02 23:58:33  vsnyder
! Require a MASS field in SPECTRA if a LINE field appears
!
! Revision 2.21  2004/04/02 01:00:01  vsnyder
! Cosmetic change
!
! Revision 2.20  2004/01/09 07:25:20  livesey
! Added the fictitious instrument mls1
!
! Revision 2.19  2003/07/15 18:17:31  livesey
! Added a 2D dump
!
! Revision 2.18  2003/05/21 22:14:53  vsnyder
! Add 'new' fields to the dump routines
!
! Revision 2.17  2003/05/19 19:58:07  vsnyder
! Remove USEs for unreferenced symbols, remove unused local variables
!
! Revision 2.16  2003/05/16 23:50:42  livesey
! Removed spec_tag, added mass.
!
! Revision 2.15  2003/05/10 22:20:57  livesey
! Tried to calm down -g1..
!
! Revision 2.14  2003/05/05 23:00:25  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.13  2003/02/27 18:09:04  bill
! Reverted from wrongly-committed newfwm branch version
!
! Revision 2.11.2.8  2003/04/24 21:56:27  vsnyder
! Determine size and allocate the database, then fill it
!
! Revision 2.11.2.7  2003/03/13 00:27:30  vsnyder
! Add a NAME argument to Dump_SpectCat_Database
!
! Revision 2.11.2.6  2003/03/12 21:24:54  vsnyder
! Add DUMP generic for Dump_SpectCat_Database.
! Add Catalog argument to Dump_SpectCat_Database.
!
! Revision 2.11.2.5  2003/03/01 03:16:51  vsnyder
! Nullify polarized field, so it doesn't end up in all lines after one gets one
!
! Revision 2.11.2.4  2003/02/27 23:20:00  vsnyder
! Add 'polarized' component to catalog_t type
!
! Revision 2.12  2003/02/27 03:25:06  vsnyder
! Add Polarized field for catalog extracts in Get_Species_Data
!
! Revision 2.11.2.3  2003/02/26 02:32:38  vsnyder
! Remove PUBLIC declaration for recently removed generic DUMP
!
! Revision 2.11.2.2  2003/02/26 00:01:55  vsnyder
! Remove ambiguous 'dump' interface that nobody used anyway
!
! Revision 2.11.2.1  2003/02/22 00:49:58  vsnyder
! Add EMLSSignalsPol field to Line
!
! Revision 2.11  2002/12/03 01:26:41  vsnyder
! Add QN field to lines spec, add SearchByQn function
!
! Revision 2.10  2002/10/08 17:08:06  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.9  2002/01/08 01:02:03  livesey
! Nullify lines by default
!
! Revision 2.8  2001/11/09 23:20:17  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.7  2001/10/18 23:53:03  livesey
! Tidied up dump, added new stuff
!
! Revision 2.6  2001/10/15 18:10:37  livesey
! Added continuum
!
! Revision 2.5  2001/10/09 22:38:41  livesey
! Added stuff for ns
!
! Revision 2.4  2001/09/19 04:38:48  livesey
! Lines per band stuff works now
!
! Revision 2.3  2001/09/18 01:25:48  livesey
! Changed emls/umls bands to emls/umls signals
!
! Revision 2.2  2001/09/18 01:23:34  livesey
! Changed bands to signals
!
! Revision 2.1  2001/09/18 00:08:25  livesey
! Added the bands information stuff
!
! Revision 2.0  2001/09/17 20:26:26  livesey
! New forward model
!
! Revision 1.8  2001/05/26 00:22:19  livesey
! Made destroy stuff return if nothing to do.
!
! Revision 1.7  2001/05/03 22:10:15  vsnyder
! Add copyright notice, make databases SAVE
!
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
