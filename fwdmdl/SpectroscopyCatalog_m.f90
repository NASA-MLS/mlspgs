! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module SpectroscopyCatalog_m

! Process the Spectroscopy section.  Read the "old format" spectroscopy catalog

  use Intrinsic, only: L_none
  use MLSCommon, only: R8
  use Molecules, only: First_Molecule, Last_Molecule

  ! More USEs below in each procedure.

  implicit none

  private
  ! Public procedures:
  public :: Spectroscopy, SearchByQn
  public :: Destroy_Line_Database, Destroy_SpectCat_Database
  public :: Dump_Lines_Database, Dump_SpectCat_Database, Dump_SpectCat_Item, Dump
  public :: Read_Spectroscopy, Write_Spectroscopy

  interface DUMP
    module procedure Dump_Line, Dump_SpectCat_Database, Dump_SpectCat_Database_2d
    module procedure Dump_SpectCat_Item
  end interface

  ! Private parameters
  integer, parameter :: MaxContinuum = 6

  ! Public types:
  type, public :: Line_T           ! One line in the spectrum for a species
    integer :: Line_Name           ! Sub_rosa index
    real(r8) :: DELTA              ! Delta interference coefficient at 300K 1/mb
    real(r8) :: EL                 ! Lower state energy cm-1
    real(r8) :: GAMMA              ! Gamma interference coefficient at 300K 1/mb
    real(r8) :: N                  ! Temperature power dependence of w
    real(r8) :: N1                 ! Temperature dependency of delta
    real(r8) :: N2                 ! Temperature dependency of gamma
    real(r8) :: NS                 ! Pressure shift dependency on temperature
    real(r8) :: PS                 ! Pressure shift parameter in MHz/mbar
    real(r8) :: STR                ! Integrated spectral intensity
                                   ! Log(nm**2 MHz) at 300 K
    real(r8) :: V0                 ! Line center frequency MHz
    real(r8) :: W                  ! Collision broadening parameter
                                   ! MHz/mbar at 300 K
    logical :: UseYi               ! Delta /= 0 .or. Gamma /= 0
    integer, dimension(:), pointer :: QN=>NULL()      ! Optional quantum numbers
    integer, dimension(:), pointer :: Signals=>NULL() ! List of signal indices for line
    integer, dimension(:), pointer :: Sidebands=>NULL() ! Sidebands for above bands (-1,0,1)
    logical, dimension(:), pointer :: Polarized=>NULL() ! Process this signal and
                                   ! sideband using the polarized model
  end type Line_T

  type, public :: Catalog_T        ! Catalog entry for a species
    real(r8) :: Continuum(MaxContinuum)  ! Continuum coefficients
    real :: DefaultIsotopeRatio = 1.0
    integer, pointer :: Lines(:)=>NULL() ! Indices in Lines database
    real(r8) :: Mass               ! Molecular mass in AMU
    integer :: Molecule = l_none   ! L_...; l_none => no catalog entry
    logical, pointer :: Polarized(:)=>NULL() ! Used only in catalog extract to
                                   ! indicate that the lines(:) are to be
                                   ! processed with the polarized model
    real(r8) :: Qlog(3)            ! Logarithm of the partition function
                                   ! at 300, 225, and 150 K
    integer :: Species_Name        ! Sub_rosa index of s_spectra label
  end type Catalog_T

  type(catalog_T), public, parameter :: Empty_Cat = catalog_t ( &
    & 0.0_r8, 1.0, NULL(), 0.0_r8, l_none, NULL(), 0.0_r8, -1 )

  ! Public Variables:
  ! The spectroscopy database:
  type(catalog_T), public, save :: Catalog(first_molecule:last_molecule)

  ! Greatest number of lines in any catalog entry:
  integer, public, save :: MostLines = 0

  ! The lines database:
  type(line_T), public, pointer, dimension(:), save :: Lines => NULL()

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
    use Init_Spectroscopy_M, only: S_Line, S_Spectra, S_ReadSpectroscopy, &
      & S_WriteSpectroscopy
    ! Now the Fields:
    use Init_Spectroscopy_M, only: F_Continuum, F_Delta, F_DefaultIsotopeRatio, &
      & F_El, F_EMLSSIGNALS, F_EMLSSIGNALSPOL, F_Gamma, F_Lines, F_Mass, &
      & F_Molecule, F_MLS1SIGNALS, F_N, F_N1, F_N2, F_Ns, F_Ps, F_Qlog, F_QN, &
      & F_Str, F_UMLSSIGNALS, F_V0, F_W
    use Intrinsic, only: L_EMLS, L_MLS1, L_UMLS, &
      & Phyq_Dimless => Phyq_Dimensionless, Phyq_Frequency, S_Time
    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, &
      & MLSMSG_DeAllocate, MLSMSG_Error
    use MoreTree, only: Get_Field_Id, Get_Spec_Id
    use Parse_Signal_m, only: PARSE_SIGNAL
    use String_Table, only: Get_string
    use Time_M, only: Time_Now
    use Toggles, only: Gen, Switches, Toggle
    use Trace_M, only: Trace_begin, Trace_end
    use Tree, only: Decorate, Decoration, Node_ID, NSons, Sub_Rosa, Subtree
    use Tree_Types, only: N_Named, N_String
    use MLSSignals_m, only: Instrument

    ! Dummy argument
    integer, intent(in) :: Root         ! Of the AST for the section

    ! Local Variables
    integer :: Error                    ! /= 0 => An error occured
    character(len=1023) :: FileName     ! For WriteSpectroscopy
    character(len=31) :: FileType
    logical :: GotLines, GotMass        ! Got a "lines" or "mass" field
    integer :: I, J, K, L               ! Loop inductors, Subscripts
    integer :: Key                      ! Index of spec_arg
    integer :: Molecule                 ! Molecule for which the catalog applies
    integer :: Name                     ! Index of name in string table
    integer :: NOSIGNALS                ! For the bands part
    integer :: NumLines                 ! Number of lines in catalog
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
    integer :: TheSignal                ! SubRosa for a signal
    integer :: THISMANY                 ! Conted up to noSignals
    logical :: TIMING                   ! For S_Time
    real :: T1, T2                      ! For S_Time
    real(r8) :: VALUE                   ! From Expr_Check

    ! Error message codes
    integer, parameter :: DupSpectra = 1               ! Duplicate s_spectra
    integer, parameter :: Negative = dupSpectra + 1    ! Parameter is negative
    integer, parameter :: No_Mass = negative + 1       ! Lines field but no mass field
    integer, parameter :: NotInt = no_mass + 1         ! QN not an integer
    integer, parameter :: NotListedSignal = NotInt + 1 ! Polarized signal is not
                                        ! listed as emlsSignal
    integer, parameter :: TooBig = NotListedSignal + 1 ! Too many elements
    integer, parameter :: WrongSize = TooBig + 1       ! Wrong number of elements
    integer, parameter :: WrongUnits = wrongSize + 1   ! Wrong physical units

    if ( toggle(gen) ) call trace_begin ( "Spectroscopy", root )

    error = 0
    timing = .false.

    ! Determine size of created or expanded lines database.
    offsetLines = 0
    if ( associated(lines) ) offsetLines = size(lines)
    numLines = offsetLines
    do i = 2, nsons(root)-1             ! Skip names of section
      son = subtree(i,root)
      if ( node_id(son) == n_named ) then
        key = subtree(2,son)
      else
        key = son
      end if
      if ( get_spec_id(key) == s_line ) numLines = numLines + 1
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

    numLines = offsetLines
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
            if ( lines(numLines)%el < 0.0_r8 ) call announce_error ( son, negative )
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
        lines(numLines)%useYi = abs(lines(numLines)%delta) > 0.0 .or. &
          &                     abs(lines(numLines)%gamma) > 0.0
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
      case ( s_readSpectroscopy ) ! ..........  ReadSpectroscopy  .....
        ! Only one field -- f_file -- allowed
        j = subtree(2,subtree(2,key))
        fileType = 'HDF5'
        if ( node_id(j) /= n_string ) then ! must be n_*colon
          call get_string ( sub_rosa(subtree(2,j)), fileType, strip=.true. )
          j = subtree(1,j)
        end if
        call get_string ( sub_rosa(j), fileName, strip=.true. )
        call Read_Spectroscopy ( j, fileName, fileType )
      case ( s_spectra ) ! .............................  SPECTRA  .....
        ! Get the molecule
        do j = 2, nsons(key)
          son = subtree(j,key)
          if ( get_field_id(son) == f_molecule ) then
            molecule = decoration(subtree(2,son))
            if ( catalog(molecule)%molecule /= l_none ) &
              & call announce_error ( son, dupSpectra, molecule )
            catalog(molecule)%molecule = molecule
        exit
          end if
        end do
        catalog(molecule)%species_Name = name
        catalog(molecule)%continuum = 0.0
        gotLines = .false.; gotMass = .false.
        do j = 2, nsons(key)
          son = subtree(j,key)
          select case ( get_field_id(son) )
          case ( f_lines )
            gotLines = .true.
            call allocate_test ( catalog(molecule)%lines, nsons(son)-1, &
              & "catalog(molecule)%Lines", moduleName )
            mostLines = max(mostLines, size(catalog(molecule)%lines) )
            do k = 2, nsons(son)
              catalog(molecule)%lines(k-1) = decoration(decoration(subtree(k,son)))
            end do
          case ( f_continuum )
            if ( nsons(son) > MaxContinuum + 1 ) &
              & call announce_error ( son, tooBig, MaxContinuum )
            do k = 2, nsons(son)
              call expr_check ( subtree(k,son), catalog(molecule)%continuum(k-1), &
                & phyq_dimless )
            end do
          case ( f_defaultIsotopeRatio )
            call expr_check ( subtree(2,son), value, phyq_dimless )
            catalog(molecule)%defaultIsotopeRatio = value
          case ( f_mass )
            gotMass = .true.
            if ( nsons(son) /= 2 ) call announce_error ( son, wrongSize, 1 )
            call expr_check ( subtree(2,son), catalog(molecule)%mass, &
              & phyq_dimless )
        ! case ( f_molecule ) ! Already done above
          case ( f_qlog )
            if ( nsons(son) /= 4 ) call announce_error ( son, wrongSize, 3 )
            do k = 2, 4
              call expr_check ( subtree(k,son), catalog(molecule)%qlog(k-1), &
                & phyq_dimless )
            end do
          end select
        end do
        if ( gotLines .and. .not. gotMass ) call announce_error ( key, no_mass )
        call decorate ( key, molecule )
      case ( s_time ) ! ...................................  TIME  .....
        if ( timing ) then
          call sayTime
        else
          call time_now ( t1 )
          timing = .true.
        end if
      case ( s_writeSpectroscopy ) ! .........  WriteSpectroscopy  .....
        ! Only one field -- f_file -- allowed
        j = subtree(2,subtree(2,key))
        fileType = 'HDF5'
        if ( node_id(j) /= n_string ) then ! must be n_*colon
          call get_string ( sub_rosa(subtree(2,j)), fileType, strip=.true. )
          j = subtree(1,j)
        end if
        call get_string ( sub_rosa(j), fileName, strip=.true. )
        call Write_Spectroscopy ( j, fileName, fileType )
      end select
    end do ! i

    if ( index(switches,'speC') /= 0 ) then
      call dump_SpectCat_database ( catalog )
      stop
    end if
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
      use Intrinsic, only: Field_Indices, Lit_Indices, Phyq_Indices
      use MoreTree, only: StartErrorMessage
      use Output_m, only: Output
      use String_Table, only: Display_String
      integer, intent(in) :: Where      ! In the tree
      integer, intent(in) :: Code       ! The error code
      integer, intent(in), optional :: More  ! In case some error messages need

      error = max(error, 1)
      call startErrorMessage ( where )
      call output ( ' Spectroscopy complained: ' )
      select case ( code )
      case ( dupSpectra )
        call display_string ( lit_indices(more), &
          & before='Duplicate SPECTRA specification for ', advance='yes' )
      case ( negative )
        call display_string ( field_indices(get_field_id(where)), before='The ' )
        call output ( ' parameter shall not be negative.', advance='yes' )
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

  ! --------------------------------------  Destroy_Line_Database  -----
  subroutine Destroy_Line_Database

    use Allocate_Deallocate, only: Deallocate_Test
    use MLSMessageModule, only: MLSMessage, MLSMSG_DeAllocate, MLSMSG_Error

    integer :: I, Status                   ! From deallocate
    if ( .not. associated(lines) ) return
    do i = 1, size(lines)
      call deallocate_test ( lines(i)%qn, "Lines(i)%QN", moduleName )
      call deallocate_test ( lines(i)%signals, "Lines(i)%Signals", moduleName )
      call deallocate_test ( lines(i)%sidebands, "Lines(i)%Sidebands", moduleName )
      call deallocate_test ( lines(i)%polarized, "Lines(i)%Polarized", moduleName )
    end do
    deallocate ( lines, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_DeAllocate // "Lines" )
  end subroutine Destroy_Line_Database

  ! ----------------------------------  Destroy_SpectCat_Database  -----
  subroutine Destroy_SpectCat_Database
    use Allocate_Deallocate, only: Deallocate_Test
    integer :: I
    do i = first_molecule, last_molecule
      call deallocate_test ( catalog(i)%lines, "Catalog(i)%lines", moduleName )
      call deallocate_test ( catalog(i)%polarized, "Catalog(i)%polarized", moduleName )
    end do ! i
    mostLines = 0
    catalog%molecule = l_none ! Clobber them all
  end subroutine Destroy_SpectCat_Database

  ! --------------------------------------------------  Dump_Line  -----
  subroutine Dump_Line ( Line )
    use Dump_0, only: Dump
    use Output_m, only: Blanks, Output
    use String_Table, only: Display_String
    type(line_t), intent(in) :: Line

    if ( line%line_name /= 0 ) then
      call display_string ( line%line_name, before='Name = ' )
      call output ( ', ' )
    end if
    call output ( line%v0,    before='V0 = ' )
    call output ( line%el,    before=', El = ' )
    call output ( line%str,   before=', Str = ', advance='yes' )
    call blanks ( 6 )
    call output ( line%w,     before='W = ' )
    call output ( line%ns,    before=', Ns = ' )
    call output ( line%ps,    before=', Ps = ' )
    call output ( line%n,     before=', N = ', advance='yes' )
    call blanks ( 6 )
    call output ( line%delta, before='Delta = ' )
    call output ( line%n1,    before=', N1 = ' )
    call output ( line%gamma, before=', Gamma = ' )
    call output ( line%n2,    before=', N2 = ', advance='yes' )
    if ( associated(line%qn) ) call dump ( line%qn, '      QN' )
    if ( associated(line%signals) ) &
      & call dump ( line%signals, '      Signals' )
    if ( associated(line%sidebands) ) &
      & call dump ( line%sidebands, '      Sidebands' )
    if ( associated(line%polarized) ) &
      & call dump ( line%polarized, '      Polarized' )

  end subroutine Dump_Line

  ! ----------------------------------------  Dump_Lines_Database  -----
  subroutine Dump_Lines_Database ( Start, End, Number )
    use Output_m, only: Blanks, Output
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
      call dump_line ( lines(i) )
    end do
  end subroutine Dump_Lines_Database

  ! -------------------------------------  Dump_SpectCat_Database  -----
  subroutine Dump_SpectCat_Database ( Catalog, Name, Sideband, Details )
    use Output_m, only: NewLine, Output

    type(catalog_T), intent(in) :: Catalog(:)
    character(len=*), intent(in), optional :: Name
    integer, intent(in), optional :: Sideband
    integer, optional, intent(in) :: Details ! <= 0 => Don't dump lines, default 0

    integer :: I                   ! Subscript, loop inductor
    character(len=3), parameter :: SB(-1:1) = (/ 'low', '   ', 'upp' /)

    call output ( 'Spectroscopy catalog' )
    if ( present(name) ) call output ( ' '//trim(name) )
    if ( present(sideband) ) call output ( ' for ' // sb(sideband) // 'er sideband' )
    call newLine
    do i = lbound(catalog,1), ubound(catalog,1)
      if ( catalog(i)%molecule == l_none ) cycle
      call output ( i, 4 )
      call output ( ': ' )
      call dump_spectCat_item ( catalog(i), details )
    end do ! i
  end subroutine Dump_SpectCat_Database

  ! -----------------------------------------  Dump_SpectCat_Item  -----
  subroutine Dump_SpectCat_Item ( Catalog, Details )
    use Dump_0, only: Dump
    use Intrinsic, only: Lit_indices
    use Output_m, only: Blanks, NewLine, Output
    use String_Table, only: Display_String, String_Length

    type(catalog_T), intent(in) :: Catalog
    integer, optional, intent(in) :: Details ! <= 0 => Don't dump lines,
      !                                      ! == 1 => Dump line names,
      !                                      !  > 1 => Dump lines, default 0

    integer :: J                   ! Subscript, loop inductor
    integer :: MyDetails
    character(len=15) :: Print
    integer :: W

    myDetails = 0
    if ( present(details) ) myDetails = details
    if ( catalog%species_name > 0 ) then
      call display_string ( catalog%species_name )
      call output ( ', ' )
    end if
    call display_string ( lit_indices(catalog%molecule), before='Species = ' )
    write ( print, '(f10.3)' ) catalog%mass
    call output ( ', Mass = ' // trim(adjustl(print)) // ', Qlog = [ ' )
    do j = 1, 3
      write ( print, '(f10.4)' ) catalog%qlog(j)
      call output ( trim(adjustl(print)) )
      if ( j < 3 ) call output ( ', ' )
    end do
    call output ( ' ]', advance='yes' )
    call blanks ( 6 )
    call output ( 'Continuum = [ ' )
    do j = 1, MaxContinuum
      write ( print, '(g15.3)' ) catalog%continuum(j)
      call output ( trim(adjustl(print)) )
      if ( j < MaxContinuum ) call output ( ', ' )
    end do
    call output ( ' ]', advance='yes' )
    call blanks ( 6 )
    call output ( catalog%defaultIsotopeRatio, before='Default Isotope Ratio = ', &
      & advance='yes' )
    if ( myDetails > 0 ) then
      call blanks ( 6 )
      call output ( 'Lines:' )
      if ( associated(catalog%lines) ) then
        if ( size(catalog%lines) > 0 ) then
          if ( details == 1 ) then
            w = 12
            do j = 1, size(catalog%lines)
              if ( w > 72 ) then
                call newLine
                call blanks ( 12 )
                w = 12
              end if
              call blanks ( 1 )
              call display_string ( lines(catalog%lines(j))%line_name )
              w = w + string_length ( lines(catalog%lines(j))%line_name )
            end do
            call newLine
          else
            call newLine
            do j = 1, size(catalog%lines)
              call dump_lines_database ( catalog%lines(j), catalog%lines(j), &
                & .false. )
            end do
          end if
        else
          call output ( ' none', advance='yes' )
        end if
      else
        call output ( ' none', advance='yes' )
      end if
      if ( associated(catalog%polarized) ) &
        & call dump ( catalog%polarized, '      Polarized:' )
    end if ! myDetails > 0
  end subroutine Dump_SpectCat_Item

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

! --------------------------------------------  Read_Spectroscopy  -----
  subroutine Read_Spectroscopy ( Where, FileName, FileType )
    use Allocate_Deallocate, only: Allocate_Test, DeAllocate_Test
!   use Declaration_Table, only: Declare, Decls, Get_Decl, Label
!   use Intrinsic, only: Phyq_Invalid
    use IO_Stuff, only: Get_Lun
    use Machine, only: IO_Error
    use MLSHDF5, only: GetHDF5DSDims, LoadFromHDF5DS, LoadPtrFromHDF5DS
    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
      & MLSMSG_Error
    use MLSSignals_m, only: MaxSigLen
    use MLSStrings, only: Capitalize
    use MoreTree, only: GetLitIndexFromString, GetStringIndexFromString
    use Parse_Signal_m, only: Parse_Signal
    use Tree, only: Null_Tree
    use HDF5, only: H5F_ACC_RDONLY_F, H5FOpen_F, H5FClose_F, HSize_T

    integer, intent(in) :: Where ! in the parse tree
    character(len=*), intent(in) :: FileName, FileType

    type(catalog_t) :: CatalogItem
    character(len=maxSigLen), pointer :: CatNames(:)
    real(r8), pointer :: Continuum(:,:)
!   type(decls) :: Decl
    logical :: Error
    integer :: FileID            ! HDF5
    integer :: I, IOSTAT, J, L, Line1, LineN, LUN, L2
    type(line_t) :: Line
    integer, pointer :: LineIndices(:)
    integer, pointer :: LineList(:)
    character(len=maxSigLen) :: LineName
    character(len=maxSigLen), pointer :: LineNames(:)
    type(catalog_t), pointer :: MyCatalog(:)
    type(line_t), pointer :: MyLines(:)
    character(len=63) :: MoleculeName
    character(len=63), pointer :: MoleculeNames(:)
    integer :: NCat              ! Number of catalog items
    integer :: NLines            ! Number of lines for a species
    integer :: NPol              ! Number of polarized flags for a line
                                 ! (zero or nSig )
    integer :: NQN, NSig         ! Numbers of quantum numbers, signals for a line
    integer, pointer :: PolarizedIndices(:) ! PolarizedIndices(i) is index in
                                 ! PolarizedList of last Polarized for line I.
    logical, pointer :: PolarizedList(:) ! Concatenation from all lines
    real(r8), pointer :: Qlog(:,:)
    integer, pointer :: QNIndices(:) ! QNIndices(i) is index in
                                 ! QNList of last QN for line I.
    integer, pointer :: QNList(:) ! Concatenation from all lines
    integer(hsize_t) :: Shp(1), Shp2(2) ! To get the shapes of datasets HDF
    integer, pointer :: SidebandList(:) ! Concatenation from all lines
    integer, dimension(:), pointer :: SigInds ! From Parse_signal
    integer , pointer:: SignalIndices(:) ! signalIndices(i) is index in
                                 ! SidebandList and SignalList of last signal
                                 ! for line I.
    integer, pointer :: SignalList(:) ! Concatenation from all lines
    character(len=MaxSigLen) :: SignalName
    character(len=MaxSigLen), pointer :: SignalNames(:)
    character(len=63) :: SpeciesName
    character(len=5) :: What

    error = .false.
    if ( capitalize(fileType) == 'HDF5' ) then
      call h5fopen_f ( trim(fileName), H5F_ACC_RDONLY_F, fileID, iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName,&
        & 'Unable to open HDF5 Spectroscopy file ' // trim(fileName) // '.' )
      ! Create or expand the Lines database
      call getHDF5DSDims ( fileID, 'Delta', Shp )
      nLines = shp(1)
      line1 = 0
      if ( associated(lines) ) line1 = size(lines)
      lineN = line1 + nLines
      allocate ( myLines(lineN), stat=iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Allocate // 'MyLines' )
      if ( associated(lines) ) then
        myLines(:line1) = lines
        deallocate ( lines, stat=iostat )
        if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_DeAllocate // 'Lines' )
        lines => myLines
      end if
      ! Fill in the expanded part
      nullify ( lineNames, polarizedIndices, polarizedList, &
        & qnIndices, qnList, signalIndices, signalList, &
        & sidebandList, sigInds, signalNames )
      call loadPtrFromHDF5DS ( fileID, 'LineNames', lineNames )
      call loadPtrFromHDF5DS ( fileID, 'SignalNames', signalNames )
      call loadPtrFromHDF5DS ( fileID, 'PolarizedList', polarizedList )
      call loadPtrFromHDF5DS ( fileID, 'PolarizedIndices', polarizedIndices, &
        & lowBound=line1 )
      call loadPtrFromHDF5DS ( fileID, 'QNList', qnList )
      call loadPtrFromHDF5DS ( fileID, 'QNIndices', qnIndices, lowBound=line1 )
      call loadPtrFromHDF5DS ( fileID, 'SignalList', signalList )
      call loadPtrFromHDF5DS ( fileID, 'SidebandList', sidebandList )
      call loadPtrFromHDF5DS ( fileID, 'SignalIndices', signalIndices, &
        & lowBound=line1 )
      call loadFromHDF5DS ( fileID, 'Delta', lines(line1+1:lineN)%delta )
      call loadFromHDF5DS ( fileID, 'EL', lines(line1+1:lineN)%el )
      call loadFromHDF5DS ( fileID, 'Gamma', lines(line1+1:lineN)%gamma )
      call loadFromHDF5DS ( fileID, 'N', lines(line1+1:lineN)%n )
      call loadFromHDF5DS ( fileID, 'N1', lines(line1+1:lineN)%n1 )
      call loadFromHDF5DS ( fileID, 'N2', lines(line1+1:lineN)%n2 )
      call loadFromHDF5DS ( fileID, 'NS', lines(line1+1:lineN)%ns )
      call loadFromHDF5DS ( fileID, 'PS', lines(line1+1:lineN)%ps )
      call loadFromHDF5DS ( fileID, 'Str', lines(line1+1:lineN)%str )
      call loadFromHDF5DS ( fileID, 'V0', lines(line1+1:lineN)%v0 )
      call loadFromHDF5DS ( fileID, 'W', lines(line1+1:lineN)%w )
      call loadFromHDF5DS ( fileID, 'UseYi', lines(line1+1:lineN)%useYi )
      do i = line1+1, lineN
        lines(i)%line_name = 0
        if ( lineNames(i) /= '' ) &
          & lines(i)%line_name = getStringIndexFromString(trim(lineNames(i)))
        ! Don't need to nullify QN, Polarized, Sidebands or Signals fields:
        ! They spring into existence nullified.
        call allocate_test ( lines(i)%qn, qnIndices(i)-qnIndices(i-1), &
          & 'Lines(i)%QN', moduleName )
        lines(i)%qn = qnList(qnIndices(i-1)+1:qnIndices(i))
        if ( signalIndices(i) /= signalIndices(i-1) ) then
          call allocate_test ( lines(i)%signals, signalIndices(i)-signalIndices(i-1), &
            & 'Lines(i)%Signals', moduleName )
          do j = 1, size(lines(i)%signals)
            call parse_signal ( trim(signalNames(signalList(signalIndices(i-1)+j))), &
              & sigInds, null_tree )
            if ( .not. associated(sigInds) ) then
              call announceError ( &
              & 'The string ' // trim(signalNames(lines(i)%signals(j))) // &
              & ' is not a signal name.' )
            else
              lines(i)%signals(j) = sigInds(1)
            end if
            nullify ( sigInds )
          end do ! j = 1, size(lines(i)%signals)
          call allocate_test ( lines(i)%sidebands, signalIndices(i)-signalIndices(i-1), &
            & 'Lines(i)%Sidebands', moduleName )
          lines(i)%sidebands = sidebandList(signalIndices(i-1)+1:signalIndices(i))
          if ( polarizedIndices(i) /= polarizedIndices(i-1) ) then
            call allocate_test ( lines(i)%polarized, polarizedIndices(i)-polarizedIndices(i-1), &
              & 'Lines(i)%Polarized', moduleName )
            lines(i)%polarized = polarizedList(polarizedIndices(i-1)+1:polarizedIndices(i))
          end if
        end if
      end do
      call deallocate_test ( lineNames, 'LineNames', moduleName )
      call deallocate_test ( signalNames, 'SignalNames', moduleName )
      call deallocate_test ( polarizedIndices, 'PolarizedIndices', moduleName )
      call deallocate_test ( polarizedList, 'PolarizedList', moduleName )
      call deallocate_test ( qnIndices, 'QNIndices', moduleName )
      call deallocate_test ( qnList, 'QNList', moduleName )
      call deallocate_test ( signalIndices, 'SignalIndices', moduleName )
      call deallocate_test ( signalList, 'SignalList', moduleName )
      call deallocate_test ( sidebandList, 'SidebandList', moduleName )
      ! Fill the catalog
      call getHDF5DSDims ( fileID, 'Continuum', shp2 )
      if ( shp2(2) /= maxContinuum ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Second dimension of continuum field of catalog has changed.' )
      allocate ( myCatalog(shp2(1)), stat=iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Allocate // 'MyCatalog' )
      nullify ( catNames, continuum, lineIndices, lineList, moleculeNames, qlog )
      call loadPtrFromHDF5DS ( fileID, 'CatNames', catNames )
      call loadPtrFromHDF5DS ( fileID, 'Continuum', continuum )
      call loadPtrFromHDF5DS ( fileID, 'LineList', lineList )
      call loadPtrFromHDF5DS ( fileID, 'LineIndices', lineIndices, lowBound=0 )
      call loadPtrFromHDF5DS ( fileID, 'MoleculeNames', moleculeNames  )
      call loadPtrFromHDF5DS ( fileID, 'Qlog', qlog )
      call loadFromHDF5DS ( fileID, 'Mass', myCatalog%mass )
      call loadFromHDF5DS ( fileID, 'IsotopeRatio', myCatalog%defaultIsotopeRatio )
      call loadFromHDF5DS ( fileID, 'Molecule', myCatalog%molecule )
      do i = 1, size(myCatalog)
        myCatalog(i)%species_name = 0
        if ( catNames(i) /= '' ) myCatalog(i)%species_name = &
          & GetStringIndexFromString(trim(catNames(i)))
        myCatalog(i)%continuum = continuum(i,:)
        myCatalog(i)%qlog = qlog(i,:)
        if ( myCatalog(i)%molecule <= 0 ) then
          myCatalog(i)%molecule = l_none
        else
          j = getLitIndexFromString(trim(moleculeNames(myCatalog(i)%molecule)))
          if ( j < first_molecule .or. j > last_molecule ) call MLSMessage ( &
            & MLSMSG_Error, moduleName, 'The string ' // &
            & trim(moleculeNames(myCatalog(i)%molecule)) // ' is not a molecule name.' )
          myCatalog(i)%molecule = j
          call allocate_test ( myCatalog(i)%lines, lineIndices(i)-lineIndices(i-1), &
            & 'MyCatalog(i)%lines', moduleName )
          myCatalog(i)%lines = lineList(lineIndices(i-1)+1:lineIndices(i)) + line1
          mostLines = max(mostLines, size(myCatalog(i)%lines))
          catalog(myCatalog(i)%molecule) = myCatalog(i)
        end if
      end do
      deallocate ( myCatalog, stat=iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Deallocate // 'MyCatalog' )
      call deallocate_test ( catNames,      'CatNames', moduleName )
      call deallocate_test ( continuum,     'Continuum', moduleName )
      call deallocate_test ( lineList,      'LineList', moduleName )
      call deallocate_test ( lineIndices,   'LineIndices', moduleName )
      call deallocate_test ( moleculeNames, 'MoleculeNames', moduleName )
      call deallocate_test ( qLog,          'QLog', moduleName )
      call H5FClose_F ( fileID, iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName,&
        & 'Unable to close HDF5 Spectroscopy file ' // trim(fileName) // '.' )
    else if ( capitalize(fileType) == 'UNFORMATTED' ) then
      ! Open the file
      call get_lun ( lun )
      what = 'open'
      open ( lun, file=trim(fileName), status='old', form='unformatted', &
        & access='sequential', iostat=iostat, err=9 )
      ! Read the catalog
      what = 'read'
      read ( lun, iostat=iostat, err=9 ) nCat, nLines
      do i = 1, nCat
        catalogItem = empty_cat
        catalogItem%species_name = 0
        read ( lun, iostat=iostat, err=9 ) l, speciesName(:l), &
          & l2, moleculeName(:l2), catalogItem%continuum, catalogItem%mass, &
          & catalogItem%defaultIsotopeRatio, catalogItem%qlog, nLines
        catalogItem%molecule = getLitIndexFromString ( moleculeName(:l2) )
        if ( catalogItem%molecule < first_molecule .or. &
          &  catalogItem%molecule > last_molecule ) &
          & call announceError ( moleculeName(:l2)// ' is not a molecule.' )
        if ( l /= 0 ) then
          catalogItem%species_name = getStringIndexFromString ( speciesName(:l) )
! This isn't of any use unless we provide a way to turn off type checking,
! and do it ourselves where we want to use these labels in, say, "dump."
!         decl = get_decl(catalogItem%species_name,label)
!         if ( decl%type == label ) then
!           call announceError ( 'The catalog item label ' // &
!             & speciesName(:l) // ' is already used as a label' )
!         else
!           call declare ( catalogItem%species_name, 0.0d0, label, phyq_invalid, &
!             & null_tree )
!         end if
        end if
        if ( nLines > 0 ) then
          call allocate_test ( catalogItem%lines, nLines, &
          & 'catalogItem%lines', moduleName )
          mostLines = max(mostLines, nLines)
        end if
        ! Read the lines for the catalog entry
        do j = 1, nLines
          nullify ( line%qn, line%signals, line%sidebands, line%polarized )
          read ( lun, iostat=iostat, err=9 ) l2, lineName(:l2), line%delta, &
            & line%el, line%gamma, line%n, line%n1, line%n2, &
            & line%ps, line%ns, line%str, line%v0, line%w, &
            & line%useYi, npol, nqn, nsig
          line%line_name = getStringIndexFromString ( lineName(:l2) )
! This isn't of any use unless we provide a way to turn off type checking,
! and do it ourselves where we want to use these labels in, say, "dump."
!         decl = get_decl(line%line_name,label)
!         if ( decl%type == label ) then
!           call announceError ( 'The line label ' // &
!             & speciesName(:l) // ' is already used as a label' )
!         else
!           call declare ( line%line_name, 0.0d0, label, phyq_invalid, null_tree )
!         end if
          if ( nQN /= 0 ) then
            call allocate_test ( line%qn, nqn, 'line%qn', moduleName )
            read ( lun, iostat=iostat, err=9 ) line%qn
          end if
          if ( nSig /= 0 ) then
            call allocate_test ( line%signals, nSig, 'line%signals', moduleName )
            call allocate_test ( line%sidebands, nSig, 'line%sidebands', moduleName )
          end if
          nullify ( sigInds )
          do l2 = 1, nSig
            read ( lun, iostat=iostat, err=9 ) l, signalName(:l), line%sidebands(l2)
            call parse_signal ( signalName(:l), sigInds, null_tree )
            if ( .not. associated(sigInds) ) then
              call announceError ( &
              & 'The string ' // signalName(:l) // ' is not a signal name.' )
            else
              line%signals(l2) = sigInds(1)
            end if
          end do
          call deallocate_test ( sigInds, 'SigInds', moduleName )
          if ( nPol /= 0 ) then ! nPol required to be == 0 or nSig
            call allocate_test ( line%polarized, nSig, 'line%polarized', moduleName )
            read ( lun, iostat=iostat, err=9 ) line%polarized
          end if
          if ( .not. error ) catalogItem%lines(j) = &
            & addLineToDatabase ( lines, line )
        end do ! j = 1, nLines
        catalog(catalogItem%molecule) = catalogItem
      end do ! i = 1, size(catalog)
      what = 'close'
      close ( lun, iostat=iostat, err=9 )
    else
      call announceError ( ' File type ' // trim(fileType) // ' is not supported.' )
      go to 99
    end if
    if ( error ) go to 99
    return

9   call io_error ( 'Unable to ' // trim(what) // ' output file ', iostat, trim(fileName) )
99  call MLSMessage ( MLSMSG_Error, moduleName, &
      & 'Error reading spectroscopy database.' )

  contains
    subroutine AnnounceError ( What )
      use MoreTree, only: StartErrorMessage
      use Output_m, only: Output
      character(len=*), intent(in) :: What
      call startErrorMessage ( where )
      call output ( what, advance='yes' )
      error = .true.
    end subroutine AnnounceError
  end subroutine Read_Spectroscopy 

  ! -------------------------------------------------  SearchByQn  -----
  integer function SearchByQn ( molecule, QNs )
  ! If the Catalog database for a molecule has an object in its Lines
  ! database that has its QN field equal to QNs, return the index of that
  ! line in the lines database.  Otherwise return zero.

    integer, intent(in) :: MOLECULE
    integer, intent(in) :: QNs(:)

    integer :: I

    if ( catalog(molecule)%molecule == molecule ) then
      if ( associated(catalog(molecule)%lines) ) then
        do i = 1, size(catalog(molecule)%lines)
          if ( associated(lines(catalog(molecule)%lines(i))%qn) ) then
            if ( size(lines(catalog(molecule)%lines(i))%qn) == size(qns) ) then
              if ( all(lines(catalog(molecule)%lines(i))%qn == qns) ) then
                searchByQn = catalog(molecule)%lines(i)
                return
              end if
            end if
          end if
        end do
      end if
    end if
    searchByQn = 0 ! not found

  end function SearchByQn

! -------------------------------------------  Write_Spectroscopy  -----
  subroutine Write_Spectroscopy ( Where, FileName, FileType )
    use Allocate_Deallocate, only: Allocate_Test, DeAllocate_Test
    use Intrinsic, only: Lit_Indices
    use IO_Stuff, only: Get_Lun
    use Machine, only: IO_Error
    use MLSHDF5, only: SaveAsHDF5DS
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSSignals_m, only: GetSignalName, MaxSigLen, Signals
    use MLSStrings, only: Capitalize
    use MoreTree, only: StartErrorMessage
    use Output_m, only: Output
    use String_Table, only: Get_String, String_Length
    use HDF5, only: H5FCREATE_F, H5FClose_F, H5F_ACC_TRUNC_F

    integer, intent(in) :: Where ! in the parse tree
    character(len=*), intent(in) :: FileName, FileType

    character(len=maxSigLen) :: CatNames(first_molecule:last_molecule)
    real(r8) :: Continuum(first_molecule:last_molecule,maxContinuum)
    integer :: FileID   ! For HDF
    integer :: I, IOSTAT, J, L, LUN, L2
    integer :: LineIndices(first_molecule-1:last_molecule)
    integer, pointer :: LineList(:)
    character(len=63) :: LineName
    character(len=maxSigLen) :: LineNames(size(lines))
    character(len=63) :: MoleculeName
    character(len=63) :: MoleculeNames(first_molecule:last_molecule)
    integer :: Molecules(first_molecule:last_molecule)
    integer :: NLines            ! Number of lines for a species
    integer :: NPol              ! Number of polarized flags for a line
                                 ! (zero or nSig )
    integer :: NQN, NSig         ! Numbers of quantum numbers, signals for a line
    integer :: PolarizedIndices(0:size(lines)) ! PolarizedIndices(i) is index in
                                 ! PolarizedList of last Polarized for line I.
    logical, pointer :: PolarizedList(:) ! Concatenation from all lines
    real(r8) :: Qlog(first_molecule:last_molecule,3)
    integer :: QNIndices(0:size(lines)) ! QNIndices(i) is index in
                                 ! QNList of last QN for line I.
    integer, pointer :: QNList(:) ! Concatenation from all lines
    integer, pointer :: SidebandList(:) ! Concatenation from all lines
    integer :: SignalIndices(0:size(lines)) ! signalIndices(i) is index in
                                 ! SidebandList and SignalList of last signal
                                 ! for line I.
    integer, pointer :: SignalList(:) ! Concatenation from all lines
    character(len=127) :: SignalName
    character(len=MaxSigLen) :: SignalNames(size(signals))
    character(len=63) :: SpeciesName
    character(len=5) :: What

    if ( .not. associated(lines) ) then
      call startErrorMessage ( where )
      call output ( ' The Lines database doesn''t exist yet.', advance='yes' )
      go to 99
    end if
    if ( .not. associated(signals) ) then
      call startErrorMessage ( where )
      call output ( ' The Signals database doesn''t exist yet.', advance='yes' )
      go to 99
    end if
    if ( capitalize(fileType) == 'HDF5' ) then

      ! Data sets for the signals:
      ! SignalNames: 1-D array of MaxSigLen characters

      ! Data sets for the lines:
      ! Delta: 1-D array real, from all lines
      ! EL: 1-D array real, from all lines
      ! Gamma: 1-D array real, from all lines
      ! LineNames: 1-D array of MaxSigLen characters
      ! N: 1-D array real, from all lines
      ! N1: 1-D array real, from all lines
      ! N2: 1-D array real, from all lines
      ! NS: 1-D array real, from all lines
      ! PolarizedIndices: 1-D array of integer, index in PolarizedList for last
      !   Polarized in lines(i)
      ! PolarizedList: 1-D array of logical, concatenated from all lines
      ! PS: 1-D array real, from all lines
      ! QNIndices: 1-D array of integer, index in QNList for last QN in lines(i)
      ! QNList: 1-D array of integer, concatenated from all lines
      ! SidebandList: 1-D array of integer, concatenated from all lines, indices
      !   in SignalNames DS
      ! SignalIndices: 1-D array of integer, index in SidebandList and SignalList
      !   for last Signal in lines(i)
      ! SignalList: 1-D array of integer, concatenated from all lines, indices
      !   in SignalNames dataset
      ! Str: 1-D array real, from all lines
      ! V0: 1-D array real, from all lines
      ! W: 1-D array real, from all lines
      ! UseYi: 1-D array real, from all lines

      ! Data sets for the molecules:
      ! CatNames: 1-d array of len=maxSigLen for catalog%species_name
      ! MoleculeNames: 1-d array of len=63 characters
      ! Continuum: 2-d array (molecules X MaxContinuum) of real(r8)
      ! LineIndices: 1-d array of integer, index in LineList for last Line in
      !   catalog(i)
      ! LineList: 1-d array of integer, concatenated from all Catalog
      ! Mass: 1-d array of real(r8)
      ! Molecule: 1-d array of integer
      ! Qlog: 2-d array (molecules x 3) of real(r8)

      ! Create the file
      call H5FCreate_F ( trim(fileName), H5F_ACC_TRUNC_F, fileID, iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create hdf5 Spectroscopy file ' // trim(fileName) )
      ! Create the SignalNames data set
      do i = 1, size(signals)
        call getSignalName ( i, signalNames(i) )
      end do
      call saveAsHDF5DS ( fileID, 'SignalNames', signalNames )
      ! Concatenate the arrays from the lines structures
      polarizedIndices(0) = 0
      qnIndices(0) = 0
      signalIndices(0) = 0
      do i = 1, size(lines)
        if ( associated(lines(i)%polarized) ) then
          polarizedIndices(i) = polarizedIndices(i-1) + size(lines(i)%polarized)
        else
          polarizedIndices(i) = polarizedIndices(i-1)
        end if
        if ( associated(lines(i)%qn) ) then
          qnIndices(i) = qnIndices(i-1) + size(lines(i)%qn)
        else
          qnIndices(i) = qnIndices(i-1)
        end if
        if ( associated(lines(i)%signals) ) then
          signalIndices(i) = signalIndices(i-1) + size(lines(i)%signals)
        else
          signalIndices(i) = signalIndices(i-1)
        end if
      end do
      nullify ( lineList, polarizedList, qnList, sidebandList, signalList )
      call allocate_test ( polarizedList, polarizedIndices(size(lines)), &
        & 'PolarizedList', moduleName )
      call allocate_test ( qnList, qnIndices(size(lines)), &
        & 'QNList', moduleName )
      call allocate_test ( sidebandList, signalIndices(size(lines)), &
        & 'SidebandList', moduleName )
      call allocate_test ( signalList, signalIndices(size(lines)), &
        & 'SignalList', moduleName )
      do i = 1, size(lines)
        if ( associated(lines(i)%polarized) ) &
          & polarizedList(polarizedIndices(i-1)+1:polarizedIndices(i)) = lines(i)%polarized
        if ( associated(lines(i)%qn) ) &
          & qnList(qnIndices(i-1)+1:qnIndices(i)) = lines(i)%qn
        if ( associated(lines(i)%signals) ) then
          sidebandList(signalIndices(i-1)+1:signalIndices(i)) = lines(i)%sidebands
          signalList(signalIndices(i-1)+1:signalIndices(i)) = lines(i)%signals
        end if
      end do
      call saveAsHDF5DS ( fileID, 'PolarizedList', polarizedList )
      call saveAsHDF5DS ( fileID, 'PolarizedIndices', polarizedIndices )
      call saveAsHDF5DS ( fileID, 'QNList', qnList )
      call saveAsHDF5DS ( fileID, 'QNIndices', qnIndices )
      call saveAsHDF5DS ( fileID, 'SignalList', signalList )
      call saveAsHDF5DS ( fileID, 'SidebandList', sidebandList )
      call saveAsHDF5DS ( fileID, 'SignalIndices', signalIndices )
      call saveAsHDF5DS ( fileID, 'Delta', lines%delta )
      call saveAsHDF5DS ( fileID, 'EL', lines%el )
      call saveAsHDF5DS ( fileID, 'Gamma', lines%gamma )
      call saveAsHDF5DS ( fileID, 'N', lines%n )
      call saveAsHDF5DS ( fileID, 'N1', lines%n1 )
      call saveAsHDF5DS ( fileID, 'N2', lines%n2 )
      call saveAsHDF5DS ( fileID, 'NS', lines%ns )
      call saveAsHDF5DS ( fileID, 'PS', lines%ps )
      call saveAsHDF5DS ( fileID, 'Str', lines%str )
      call saveAsHDF5DS ( fileID, 'V0', lines%v0 )
      call saveAsHDF5DS ( fileID, 'W', lines%w )
      call saveAsHDF5DS ( fileID, 'UseYi', lines%useYi )
      call deallocate_test ( polarizedList, 'PolarizedList', moduleName )
      call deallocate_test ( qnList, 'QNList', moduleName )
      call deallocate_test ( sidebandList, 'SidebandList', moduleName )
      call deallocate_test ( signalList, 'SignalList', moduleName )
      ! Create the LineNames data set
      do i = 1, size(lines)
        call get_string ( lines(i)%line_name, lineNames(i) )
      end do
      call saveAsHDF5DS ( fileID, 'LineNames', lineNames )
      ! Create datasets for the catalog
      lineIndices(first_molecule-1) = 0
      do i = first_molecule, last_molecule
        catNames(i) = ''
        if ( catalog(i)%species_name /= 0 ) &
          & call get_string ( catalog(i)%species_name, catNames(i) )
        if ( catalog(i)%molecule /= l_none .and. associated(catalog(i)%lines) ) then
          lineIndices(i) = lineIndices(i-1) + size(catalog(i)%lines)
        else
          lineIndices(i) = lineIndices(i-1)
        end if
        call get_string ( lit_indices(i), moleculeNames(i) )
      end do
      call saveAsHDF5DS ( fileID, 'CatNames', catNames )
      call saveAsHDF5DS ( fileID, 'MoleculeNames', moleculeNames )
      call allocate_test ( lineList, lineIndices(last_molecule), &
        & 'LineList', moduleName )
      do i = first_molecule, last_molecule
        if ( catalog(i)%molecule /= l_none .and. associated(catalog(i)%lines) ) &
          & lineList(lineIndices(i-1)+1:lineIndices(i)) = catalog(i)%lines
        continuum(i,:) = catalog(i)%continuum
        qlog(i,:) = catalog(i)%qlog
      end do
      call saveAsHDF5DS ( fileID, 'Continuum', continuum )
      call saveAsHDF5DS ( fileID, 'LineList', lineList )
      call saveAsHDF5DS ( fileID, 'LineIndices', lineIndices )
      call saveAsHDF5DS ( fileID, 'Mass', catalog%mass )
      call saveAsHDF5DS ( fileID, 'IsotopeRatio', catalog%defaultIsotopeRatio )
      ! Molecule indexes MoleculeNames, not the lits -- the lits may change.
      molecules = catalog%molecule
      where ( molecules == l_none )
        molecules = -1
      elsewhere
        molecules = molecules - first_molecule+1
      end where
      call saveAsHDF5DS ( fileID, 'Molecule', molecules )
      call saveAsHDF5DS ( fileID, 'Qlog', qlog )
      call deallocate_test ( lineList, 'LineList', moduleName )
      ! Close the file
      call H5FClose_F ( fileID, iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName,&
        & 'Unable to close hdf5 Spectroscopy file ' // trim(fileName) // '.' )
    else if ( capitalize(fileType) == 'UNFORMATTED' ) then
      ! Open the file
      call get_lun ( lun )
      what = 'open'
      open ( lun, file=trim(fileName), status='unknown', form='unformatted', &
        & access='sequential', iostat=iostat, err=9 )
      ! Write the catalog
      what = 'write'
      write ( lun, iostat=iostat, err=9 ) count(catalog%molecule/=l_none), size(lines)
      do i = first_molecule, last_molecule
        if ( catalog(i)%molecule == l_none ) cycle
        l = 0
        if ( catalog(i)%species_name /= 0 ) then
          l = string_length(catalog(i)%species_name)
          call get_string ( catalog(i)%species_name, speciesName )
        end if
        l2 = string_length(lit_indices(catalog(i)%molecule))
        call get_string ( lit_indices(catalog(i)%molecule), moleculeName )
        nLines = 0
        if ( associated(catalog(i)%lines) ) nLines = size(catalog(i)%lines)
        write ( lun, iostat=iostat, err=9 ) l, speciesName(:l), &
          & l2, moleculeName(:l2), catalog(i)%continuum, catalog(i)%mass, &
          & catalog(i)%defaultIsotopeRatio, catalog(i)%qlog, nLines
        ! Write the lines for the catalog entry
        do j = 1, nLines
          l = catalog(i)%lines(j)
          l2 = string_length(lines(l)%line_name)
          call get_string ( lines(l)%line_name, lineName )
          nQN = 0
          if ( associated(lines(l)%qn) ) nqn = size(lines(l)%qn)
          nSig = 0
          if ( associated(lines(l)%signals) ) nSig = size(lines(l)%signals)
          nPol = 0
          if (  associated(lines(l)%polarized) ) nPol = nSig
          write ( lun, iostat=iostat, err=9 ) l2, lineName(:l2), lines(l)%delta, &
            & lines(l)%el, lines(l)%gamma, lines(l)%n, lines(l)%n1, lines(l)%n2, &
            & lines(l)%ps, lines(l)%ns, lines(l)%str, lines(l)%v0, lines(l)%w, &
            & lines(l)%useYi, nPol, nQN, nSig
          if ( nQN /= 0 ) write ( lun, iostat=iostat, err=9 ) lines(l)%qn
          do l2 = 1, nSig
            call getSignalName ( lines(l)%signals(l2), signalName )
            write ( lun, iostat=iostat, err=9 ) len_trim(signalName), &
              & trim(signalName), lines(l)%sidebands(l2)
          end do
          if ( nPol /= 0 ) write ( lun, iostat=iostat, err=9 ) lines(l)%polarized
        end do ! j = 1, nLines
      end do ! i = 1, size(catalog)
      what = 'close'
      close ( lun, iostat=iostat, err=9 )
    else
      call startErrorMessage ( where )
      call output ( ' File type  ' // trim(fileType) // ' is not supported.', &
        & advance='yes' )
      go to 99
    end if
    return

9   call io_error ( 'Unable to ' // trim(what) // ' output file ', iostat, trim(fileName) )
99  call MLSMessage ( MLSMSG_Error, moduleName, &
      & 'Error writing spectroscopy database.' )
  end subroutine Write_Spectroscopy 

! ------------------------------------------------  not_used_here  -----
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module SpectroscopyCatalog_m

! $Log$
! Revision 2.36  2005/01/20 02:28:45  vsnyder
! Cannonball polishing
!
! Revision 2.35  2005/01/13 01:32:24  vsnyder
! Add DefaultIsotopeRatio to Catalog_t
!
! Revision 2.34  2005/01/12 03:10:00  vsnyder
! Add error message in Read_Spectroscopy if the file cannot be opened
!
! Revision 2.33  2005/01/03 18:56:31  pwagner
! Moved use hdf5 statements to avoid Lahey internal compiler error
!
! Revision 2.32  2004/12/31 02:40:44  vsnyder
! Read/Write HDF Spectroscopy catalog
!
! Revision 2.31  2004/12/13 20:44:22  vsnyder
! Added Dump_Line and Dump_SpectCat_Item, added UseYi field to Line_T.
! Improved error checking.
!
! Revision 2.30  2004/11/16 03:05:18  vsnyder
! Not using correct bounds for dump
!
! Revision 2.29  2004/11/04 03:56:35  vsnyder
! Correct a blunder
!
! Revision 2.28  2004/11/04 03:40:42  vsnyder
! Index spectroscopy catalog by molecule instead of searching
!
! Revision 2.27  2004/11/01 20:26:35  vsnyder
! Reorganization of representation for molecules and beta groups; PFA may be broken for now
!
! Revision 2.26  2004/09/16 18:35:02  vsnyder
! Pass 'details' through Dump_SpectCat_Database_2d
!
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
