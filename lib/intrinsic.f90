! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module INTRINSIC

! Intrinsic constants needed by Init_Tables_Module, DeclarationTable, etc.

! Declaring the definitions is handled by the tree walker.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test

  implicit NONE
  public
  private :: Allocate_Test, Deallocate_Test ! may make .mod files smaller

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

! A "spec_def" vertex may be decorated with (sums of) the following flags:
  integer, parameter :: NO_DUP = 1        ! Duplicate fields prohibited
  integer, parameter :: ALL_FIELDS = 2    ! All fields required
  integer, parameter :: NO_POSITIONAL = 4 ! Positional fields prohibited
! A "field_type", "field_spec" or "dot" vertex may be decorated with the
! following flag:
  integer, parameter :: REQ_FLD = 1       ! Required field

! Data types that don't have enumerated literals:
  integer, parameter :: T_FIRST             = 1
  integer, parameter :: T_NUMERIC           = t_first
  integer, parameter :: T_NUMERIC_RANGE     = t_numeric + 1
  integer, parameter :: T_STRING            = t_numeric_range + 1
! Enumeration types:
  integer, parameter :: T_BOOLEAN           = t_string + 1
  integer, parameter :: T_INSTRUMENT        = t_boolean + 1
  integer, parameter :: LAST_INTRINSIC_TYPE = t_instrument

! We don't define any fields here, but here's the first index:
  integer, parameter :: Field_First = 1

! Abstract physical quantities:
  integer, parameter :: PHYQ_INVALID =         0 ! Invalid unit given by user
  integer, parameter :: PHYQ_DIMENSIONLESS =   phyq_invalid+1     ! Dimensionless quantity
  integer, parameter :: PHYQ_LENGTH =          phyq_dimensionless+1  ! Default meters
  integer, parameter :: PHYQ_TIME =            phyq_length+1         ! Default seconds
  integer, parameter :: PHYQ_PRESSURE =        phyq_time+1        !  Default millibars
  integer, parameter :: PHYQ_TEMPERATURE =     phyq_pressure+1    ! Default Kelvins
  integer, parameter :: PHYQ_VMR =             phyq_temperature+1 ! Default parts-per-one
  integer, parameter :: PHYQ_ANGLE =           phyq_vmr+1         ! Default degrees
  integer, parameter :: PHYQ_MAFS =            phyq_angle+1       ! Default MAFs
  integer, parameter :: PHYQ_MIFS =            phyq_mafs+1        ! Default MIFs
  integer, parameter :: PHYQ_FREQUENCY =       phyq_mifs+1        ! Default MHz
  integer, parameter :: PHYQ_ZETA =            phyq_frequency+1   ! log10(pressure/hPa)
  integer, parameter :: PHYQ_VELOCITY =        phyq_zeta+1        ! Default meters/second
  integer, parameter :: PHYQ_EXTINCTION =      phyq_velocity+1    ! Default 1/meters
  integer, parameter :: PHYQ_ICEDENSITY =      phyq_extinction+1  ! Default g/meters^3
  integer, parameter :: PHYQ_DOBSONUNITS =     phyq_icedensity+1  ! 1 DU = 2.687e20 molecules/m^2
  integer, parameter :: FIRST_PHYQ = phyq_invalid, LAST_PHYQ = PHYQ_DobsonUnits
  integer :: PHYQ_INDICES(first_phyq:last_phyq)

! Enumeration literals:
  integer, parameter :: FIRST_LIT       = 1
! Don't edit the following file directly--it is generated automatically
! based on the file lit_names.txt
  include 'lit_parm.f9h'

  ! Specifications
  integer, parameter :: Spec_First = 1
  integer, parameter :: S_TIME          = Spec_First
  integer, parameter :: Last_Intrinsic_Spec = S_Time

  ! The following parameters are for building trees:
  integer, parameter :: BEGIN = -1       ! Start of a tree
  integer, parameter :: D = 1000000      ! Decoration
  integer, parameter :: F = 1000         ! Field index
  integer, parameter :: L = 2000         ! Lit index
  integer, parameter :: N = 0            ! Tree index
  integer, parameter :: NADP = n+d*(all_fields+no_dup+no_positional)
  integer, parameter :: ND = n+d*no_dup
  integer, parameter :: NDP = n+d*(no_dup+no_positional)
  integer, parameter :: NP = n+d*no_positional
  integer, parameter :: NR = n+d*req_fld
  integer, parameter :: P = 3000         ! Parameter index
  integer, parameter :: S = 4000         ! Spec index
  integer, parameter :: T = 5000         ! Type index
  integer, parameter :: Z = 6000         ! Section index

  ! Tables used for type checking:
  integer, save, pointer, dimension(:) :: DATA_TYPE_INDICES=>NULL()
  integer, save, pointer, dimension(:) :: FIELD_INDICES=>NULL()
  integer, save, pointer, dimension(:) :: LIT_INDICES=>NULL()
  integer, save, pointer, dimension(:) :: PARM_INDICES=>NULL()
  integer, save, pointer, dimension(:) :: SECTION_INDICES=>NULL()
  integer, save, pointer, dimension(:) :: SPEC_INDICES=>NULL()

contains ! =====     Public procedures     =============================
! -----------------------------------------------  INIT_INTRINSIC  -----
  subroutine INIT_INTRINSIC ( N_DATA_TYPE_INDICES, N_FIELD_INDICES, &
    & N_LIT_INDICES, FIRST_PARM_INDEX, LAST_PARM_INDEX, N_SECTION_INDICES, &
    & N_SPEC_INDICES )

    ! This really belongs in make_tree, but "make depends" can't see it there
    ! (because of the "include"):
    use TREE, only: BUILD_TREE, PUSH_PSEUDO_TERMINAL
    use TREE_TYPES, only: N_DT_DEF

    integer, intent(in) :: N_DATA_TYPE_INDICES
    integer, intent(in) :: N_FIELD_INDICES
    integer, intent(in) :: N_LIT_INDICES
    integer, intent(in) :: FIRST_PARM_INDEX, LAST_PARM_INDEX
    integer, intent(in) :: N_SECTION_INDICES
    integer, intent(in) :: N_SPEC_INDICES

    ! Allocate the string index tables for the various categories of
    ! names
    call allocate_test ( data_type_indices, n_data_type_indices, &
      & 'DATA_TYPE_INDICES', moduleName )
    call allocate_test ( field_indices, n_field_indices, &
      & 'FIELD_INDICES', moduleName )
    call allocate_test ( lit_indices, n_lit_indices, &
      & 'LIT_INDICES', moduleName )
    call allocate_test ( parm_indices, last_parm_index, &
      & 'PARM_INDICES', moduleName, lowBound=first_parm_index )
    call allocate_test ( section_indices, n_section_indices, &
      & 'SECTION_INDICES', moduleName )
    call allocate_test ( spec_indices, n_spec_indices, &
      & 'SPEC_INDICES', moduleName )

    ! Put intrinsic predefined identifiers into the symbol table.

    ! Put intrinsic non-enumeration type names into symbol table
    data_type_indices(t_numeric) =         add_ident ( 'numeric' )
    data_type_indices(t_numeric_range) =   add_ident ( 'numeric_range' )
    data_type_indices(t_string) =          add_ident ( 'string' )
    ! Put intrinsic enumeration type names into the symbol table
    data_type_indices(t_boolean) =         add_ident ( 'boolean' )
    data_type_indices(t_instrument) =      add_ident ( 'instrument' )
    ! Put intrinsic enumeration literals into the symbol table:
! Don't edit the following file directly--it is generated automatically
! based on the file lit_names.txt
    include 'lit_add.f9h'

    ! Put spec names into the symbol table
    spec_indices(s_time) =                 add_ident ( 'time' )

    ! Put abstract physical quantities into the symbol table
    phyq_indices(phyq_invalid) =           add_ident ( 'invalid' )
    phyq_indices(phyq_dimensionless) =     add_ident ( 'dimensionless' )
    phyq_indices(phyq_length) =            add_ident ( 'length' )
    phyq_indices(phyq_time) =              add_ident ( 'time' )
    phyq_indices(phyq_pressure) =          add_ident ( 'pressure' )
    phyq_indices(phyq_temperature) =       add_ident ( 'temperature' )
    phyq_indices(phyq_vmr) =               add_ident ( 'vmr' )
    phyq_indices(phyq_angle) =             add_ident ( 'angle' )
    phyq_indices(phyq_mafs) =              add_ident ( 'mafs' )
    phyq_indices(phyq_mifs) =              add_ident ( 'mifs' )
    phyq_indices(phyq_frequency) =         add_ident ( 'frequency' )
    phyq_indices(phyq_velocity) =          add_ident ( 'velocity' )
    phyq_indices(phyq_zeta) =              add_ident ( 'zeta' )
    phyq_indices(phyq_extinction) =        add_ident ( 'extinction' )
    phyq_indices(phyq_icedensity) =        add_ident ( 'icedensity' )
    phyq_indices(phyq_dobsonunits) =       add_ident ( 'dobsonunits' )

  ! Definitions are represented by trees.  The notation in the comments
  ! for the trees is < root first_son ... last_son >.  This is sometimes
  ! called "Cambridge Polish Notation."  It was developed to represent
  ! LISP by McCarthy et. al. at MIT (in Cambridge, MA).

  ! Notice that in the argument for make_tree, the tree node id is at
  ! the END of the subtree, while in Cambridge Polish Notation it is at
  ! the BEGINNING of the subtree!

  ! Put the definition trees into the tree space before the parser runs.
  ! After the parsing is done, they're automatically "glued in" to the
  ! "left" of the trees that represent the input.  The tree-walker
  ! stumbles upon them in its normal course of operation, never really
  ! realizing they're special (because by then they're not).

  ! Start with the definitions of types. These are represented by trees of
  ! the form  < n_dt_def t_type_name l_lit ... l_lit >

  ! Start with the definitions of types. These are represented by trees of
  ! the form  < n_dt_def t_type_name l_lit ... l_lit >

    ! Define the intrinsic data types
    call make_tree ( (/ &
      begin, t+t_numeric, n+n_dt_def, &
      begin, t+t_numeric_range, n+n_dt_def, &
      begin, t+t_string, n+n_dt_def /) )
    ! Define the enumerated types
    call make_tree ( (/ &
      begin, t+t_boolean, l+l_true, l+l_false, n+n_dt_def /) )
    call make_tree ( (/ &
      begin, t+t_instrument, l+l_emls, l+l_umls, n+n_dt_def /) )

  contains
    ! ................................................  MAKE_TREE  .....
    include "make_tree.f9h"

  end subroutine INIT_INTRINSIC

  ! --------------------------------------------------  Add_Ident  -----
  integer function ADD_IDENT ( TEXT )
    use SYMBOL_TABLE, only: ENTER_TERMINAL
    use SYMBOL_TYPES, only: T_IDENTIFIER
    character(len=*), intent(in) :: TEXT
    add_ident = enter_terminal ( text, t_identifier )
  end function ADD_IDENT

  ! -----------------------------------  DestroyTypeCheckerTables  -----
  subroutine DestroyTypeCheckerTables
    call deallocate_test ( data_type_indices, 'DATA_TYPE_INDICES', moduleName )
    call deallocate_test ( field_indices,     'FIELD_INDICES',     moduleName )
    call deallocate_test ( lit_indices,       'LIT_INDICES',       moduleName )
    call deallocate_test ( parm_indices,      'PARM_INDICES',      moduleName )
    call deallocate_test ( section_indices,   'SECTION_INDICES',   moduleName )
    call deallocate_test ( spec_indices,      'SPEC_INDICES',      moduleName )
  end subroutine

end module INTRINSIC

! $Log$
! Revision 2.42  2001/10/04 22:12:57  pwagner
! Now includes files lit_add.f9h and lit_parm.f9h
!
! Revision 2.41  2001/10/03 18:32:44  vsnyder
! OOPS, forgot some lit_indices(...) = add_ident(...)
!
! Revision 2.40  2001/10/03 17:38:11  vsnyder
! OOPS, defined L_DNWT_SQ twice
!
! Revision 2.39  2001/10/03 17:36:47  vsnyder
! Add lits for DNWT quantities
!
! Revision 2.38  2001/10/02 23:39:39  vsnyder
! Add L_DegreesOfFreedom
!
! Revision 2.37  2001/09/17 23:14:14  livesey
! Bug fix, added name for t_instrument
!
! Revision 2.36  2001/09/17 22:53:23  livesey
! Added t_instrument, l_emls and l_umls
!
! Revision 2.35  2001/07/30 23:28:38  pwagner
! Added columnAbundances scaffolding--needs fleshing out
!
! Revision 2.34  2001/07/18 23:15:31  dwu
! rename l_radiusofearth as l_earthradius
!
! Revision 2.33  2001/07/17 18:53:31  jonathan
! remove earthradius (redundant quantity),jonathan/wu
!
! Revision 2.32  2001/07/13 19:04:56  dwu
! fix problem after adding lostransfunc
!
! Revision 2.31  2001/07/13 18:20:18  dwu
! add quantity losTransFunc
!
! Revision 2.30  2001/07/10 23:46:35  jonathan
! added l_icedensity, paul/jonathan
!
! Revision 2.29  2001/07/06 18:55:40  jonathan
! Modified for cloud model, Paul/Jonathan
!
! Revision 2.28  2001/05/31 22:07:33  livesey
! Updated cloud quantity types
!
! Revision 2.27  2001/05/31 20:27:37  livesey
! New vector type associated with cloud quantities.
!
! Revision 2.26  2001/05/29 22:45:45  livesey
! Moved some state vector component literals in from l2/init_tables_module.f90
!
! Revision 2.25  2001/05/10 23:26:45  livesey
! Added isotope ratio vector quantity type.
!
! Revision 2.24  2001/05/03 22:27:04  livesey
! Added l_heightoffset
!
! Revision 2.23  2001/04/26 16:18:39  livesey
! Fixed bug. All those nice integer arrays weren't initially nullified, whoops!
!
! Revision 2.22  2001/04/26 02:33:03  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.21  2001/04/23 20:57:42  vsnyder
! Move the first spec (time) to 'intrinsic'
!
! Revision 2.20  2001/04/09 20:59:35  vsnyder
! Add C (for Celsius) unit and l_c name for it
!
! Revision 2.19  2001/04/09 14:56:04  perun
! Corrected enumeration literal definition for l_temperature
!
! Revision 2.18  2001/04/04 17:56:42  vsnyder
! Insert "USE TREE" because "make depends" can't see the one in "make_tree"
! (because of the "include").
!
! Revision 2.17  2001/04/04 17:21:12  pwagner
! Added extra use tree line to tweak dependencies
!
! Revision 2.16  2001/04/03 19:09:12  vsnyder
! Change the order of initialization to intrinsic, Molecules, MLSSignals.
! Use the revised make_tree.f9h, which requires revision of init...
! calling sequences.
!
! Revision 2.15  2001/03/17 02:23:40  livesey
! Bug fix, defined phyq_indices(phyq_velocity)
!
! Revision 2.14  2001/03/15 18:41:04  livesey
! Added some more, losvel etc.
!
! Revision 2.13  2001/03/15 07:37:35  livesey
! Added l_frequency
!
! Revision 2.12  2001/03/14 02:05:52  vsnyder
! Moved MLSSignals_m to mlspgs/lib.
!
! Revision 2.11  2001/03/06 22:41:59  livesey
! Minor changes
!
! Revision 2.10  2001/03/02 01:32:21  livesey
! Added some new PHYQs
!
! Revision 2.9  2001/02/22 23:57:38  vsnyder
! Remove ", public" from parameters, because default accessibility is public
!
! Revision 2.8  2001/02/09 18:37:37  vsnyder
! Add REQ_FLD flag for specification definitions
!
! Revision 2.7  2001/02/09 01:05:18  livesey
! Thought I'd done this
!
! Revision 2.6  2001/02/08 21:11:13  vsnyder
! Move "theta" from init_tables_module to intrinsic.
!
! Revision 2.5  2001/02/05 21:18:57  vsnyder
! Add parameters for type checking rules.
!
! Revision 2.4  2001/02/01 20:18:50  vsnyder
! Correct index and spelling for gph_precision and temperature_precision
!
! Revision 2.3  2001/02/01 01:23:18  vsnyder
! Account for the Molecules module
!
! Revision 2.2  2001/01/31 23:32:31  vsnyder
! Moved l_temperature l_temperature_prec l_ptan l_tangentheight l_sidebandratio
! l_scvel l_orbitinclination l_geodaltitude l_radiance l_scanresidual l_gph
! l_gph_precision l_refgph l_baseline l_extinction l_linewidth from L2's
! init_tables_module
!
! Revision 2.1  2000/10/11 18:24:39  vsnyder
! Initial entry
!
