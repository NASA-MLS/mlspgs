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
  integer, parameter :: LAST_INTRINSIC_TYPE = t_boolean

! We don't define any fields here, but here's the first index:
  integer, parameter :: Field_First = 1

! Abstract physical quantities:
  integer, parameter :: PHYQ_INVALID = 0 ! Invalid unit given by user
  integer, parameter :: PHYQ_DIMENSIONLESS = 1 ! Dimensionless quantity
  integer, parameter :: PHYQ_LENGTH = 2        ! Default meters
  integer, parameter :: PHYQ_TIME = 3          ! Default seconds
  integer, parameter :: PHYQ_PRESSURE = 4      ! Default millibars
  integer, parameter :: PHYQ_TEMPERATURE = 5   ! Default Kelvins
  integer, parameter :: PHYQ_VMR = 6           ! Default parts-per-one
  integer, parameter :: PHYQ_ANGLE = 7         ! Default degrees
  integer, parameter :: PHYQ_MAFS = 8          ! Default MAFs
  integer, parameter :: PHYQ_MIFS = 9          ! Default MIFs
  integer, parameter :: PHYQ_FREQUENCY = 10    ! Default MHz
  integer, parameter :: PHYQ_ZETA = 11         ! log10(pressure/hPa)
  integer, parameter :: PHYQ_VELOCITY = 12     ! Default meters/second
  integer, parameter :: PHYQ_EXTINCTION =13    ! Default 1/meters
  integer, parameter :: FIRST_PHYQ = phyq_invalid, LAST_PHYQ = phyq_extinction
  integer :: PHYQ_INDICES(first_phyq:last_phyq)

! Enumeration literals:
  integer, parameter :: FIRST_LIT       = 1
  integer, parameter :: L_BASELINE      = first_lit
  integer, parameter :: L_C             = l_baseline + 1
  integer, parameter :: L_CHANNEL       = l_c + 1
  integer, parameter :: L_CLOUDINDUCEDRADIANCE = l_channel + 1
  integer, parameter :: L_CLOUDOPTICALDEPTH    = l_cloudInducedRadiance + 1
  integer, parameter :: L_CLOUDSENSITIVITY     = l_cloudOpticalDepth + 1
  integer, parameter :: L_DAYS          = l_cloudSensitivity + 1
  integer, parameter :: L_DEG           = l_days + 1
  integer, parameter :: L_DEGREES       = l_deg + 1
  integer, parameter :: L_DIMENSIONLESS = l_degrees + 1
  integer, parameter :: L_DIMLESS       = l_dimensionless + 1
  integer, parameter :: L_DL            = l_dimless + 1
  integer, parameter :: L_EARTHREFL     = l_dl + 1
  integer, parameter :: L_EFFECTIVEOPTICALDEPTH = l_earthRefl + 1
  integer, parameter :: L_ELEVOFFSET    = l_effectiveOpticalDepth + 1
  integer, parameter :: L_EXTINCTION    = l_elevOffset + 1
  integer, parameter :: L_FALSE         = l_extinction + 1
  integer, parameter :: L_FREQUENCY     = l_false + 1
  integer, parameter :: L_GEODALTITUDE  = l_frequency + 1
  integer, parameter :: L_GEODANGLE     = l_geodaltitude + 1
  integer, parameter :: L_GHZ           = l_geodangle + 1
  integer, parameter :: L_GPH           = l_ghz + 1
  integer, parameter :: L_GPH_PRECISION = l_gph + 1
  integer, parameter :: L_HEIGHTOFFSET  = l_gph_precision + 1
  integer, parameter :: L_HOURS         = l_heightOffset + 1
  integer, parameter :: L_HPA           = l_hours + 1
  integer, parameter :: L_HZ            = l_hpa + 1
  integer, parameter :: L_INTERMEDIATEFREQUENCY= l_hz + 1
  integer, parameter :: L_ISOTOPERATIO  = l_intermediatefrequency + 1
  integer, parameter :: L_K             = l_isotopeRatio + 1
  integer, parameter :: L_KHZ           = l_k  + 1
  integer, parameter :: L_KM            = l_khz + 1
  integer, parameter :: L_LINEWIDTH     = l_km + 1
  integer, parameter :: L_LOGP          = l_linewidth + 1
  integer, parameter :: L_LOSVEL        = l_logp + 1
  integer, parameter :: L_LSBFREQUENCY  = l_losvel + 1
  integer, parameter :: L_M             = l_lsbfrequency + 1
  integer, parameter :: L_MAF           = l_m + 1
  integer, parameter :: L_MAFS          = l_maf + 1
  integer, parameter :: L_MASSMEANDIAMETERICE   = l_mafs + 1
  integer, parameter :: L_MASSMEANDIAMETERWATER = l_massMeanDiameterIce + 1
  integer, parameter :: L_MB            = l_massMeanDiameterWater + 1
  integer, parameter :: L_METERS        = l_mb + 1
  integer, parameter :: L_MHZ           = l_meters  + 1
  integer, parameter :: L_MIF           = l_mhz + 1
  integer, parameter :: L_MIFS          = l_mif + 1
  integer, parameter :: L_MINUTES       = l_mifs + 1
  integer, parameter :: L_NONE          = l_minutes + 1
  integer, parameter :: L_ORBITINCLINATION = l_none + 1
  integer, parameter :: L_ORBITS        = l_orbitinclination + 1
  integer, parameter :: L_PA            = l_orbits + 1
  integer, parameter :: L_PPBV          = l_pa + 1
  integer, parameter :: L_PPMV          = l_ppbv + 1
  integer, parameter :: L_PPTV          = l_ppmv + 1
  integer, parameter :: L_PTAN          = l_pptv + 1
  integer, parameter :: L_RAD           = l_ptan + 1
  integer, parameter :: L_RADIANCE      = l_rad + 1
  integer, parameter :: L_RADIANS       = l_radiance + 1
  integer, parameter :: L_REFGPH        = l_radians + 1
  integer, parameter :: L_S             = l_refgph + 1
  integer, parameter :: L_SCANRESIDUAL  = l_s + 1
  integer, parameter :: L_SCECI         = l_scanresidual + 1
  integer, parameter :: L_SCGEOCALT     = l_scECI + 1
  integer, parameter :: L_SCVEL         = l_scGeocAlt + 1
  integer, parameter :: L_SECONDS       = l_scvel + 1
  integer, parameter :: L_SIDEBANDRATIO = l_seconds + 1
  integer, parameter :: L_SPACERADIANCE = l_sidebandratio + 1
  integer, parameter :: L_TEMPERATURE   = l_spaceRadiance + 1
  integer, parameter :: L_TEMPERATURE_PREC = l_temperature + 1
  integer, parameter :: L_THETA         = l_temperature_prec + 1
  integer, parameter :: L_THZ           = l_theta + 1
  integer, parameter :: L_TIME          = l_thz + 1
  integer, parameter :: L_TNGTECI       = l_time + 1
  integer, parameter :: L_TNGTGEOCALT   = l_tngteci + 1
  integer, parameter :: L_TNGTGEODALT   = l_tngtgeocalt + 1
  integer, parameter :: L_TOTALOPTICALDEPTH =  l_tngtgeodalt + 1
  integer, parameter :: L_TRUE         =  l_totalOpticalDepth + 1
  integer, parameter :: L_USBFREQUENCY  = l_true + 1
  integer, parameter :: L_VMR           = l_usbfrequency + 1
  integer, parameter :: L_XYZ           = l_vmr + 1
  integer, parameter :: L_ZETA          = l_xyz + 1
  integer, parameter :: LAST_INTRINSIC_LIT = l_zeta

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
    ! Put intrinsic enumeration literals into the symbol table:
    lit_indices(l_baseline) =              add_ident ( 'baseline' )
    lit_indices(l_c) =                     add_ident ( 'C' )
    lit_indices(l_channel) =               add_ident ( 'channel' )
    lit_indices(l_cloudInducedRadiance) =  add_ident ( 'cloudInducedRadiance' )
    lit_indices(l_cloudOpticalDepth) =     add_ident ( 'cloudOpticalDepth' )
    lit_indices(l_cloudSensitivity) =      add_ident ( 'cloudSensitivity' )
    lit_indices(l_days) =                  add_ident ( 'days' )
    lit_indices(l_deg) =                   add_ident ( 'deg' )
    lit_indices(l_degrees) =               add_ident ( 'degrees' )
    lit_indices(l_dimensionless) =         add_ident ( 'dimensionless' )
    lit_indices(l_dimless) =               add_ident ( 'dimless' )
    lit_indices(l_dl) =                    add_ident ( 'dl' )
    lit_indices(l_earthRefl) =             add_ident ( 'earthRefl' )
    lit_indices(l_elevOffset) =            add_ident ( 'elevOffset' )
    lit_indices(l_effectiveOpticalDepth) = add_ident ( 'effectiveOpticalDepth' )
    lit_indices(l_extinction) =            add_ident ( 'extinction' )
    lit_indices(l_false) =                 add_ident ( 'false' )
    lit_indices(l_frequency) =             add_ident ( 'frequency' )
    lit_indices(l_geodaltitude) =          add_ident ( 'geodAltitude' )
    lit_indices(l_geodangle) =             add_ident ( 'geodAngle' )
    lit_indices(l_ghz) =                   add_ident ( 'GHz' )
    lit_indices(l_gph) =                   add_ident ( 'gph' )
    lit_indices(l_gph_precision) =         add_ident ( 'gph_precision' )
    lit_indices(l_heightOffset) =          add_ident ( 'heightOffset' )
    lit_indices(l_hours) =                 add_ident ( 'hours' )
    lit_indices(l_hpa) =                   add_ident ( 'hPa' )
    lit_indices(l_hz) =                    add_ident ( 'Hz' )
    lit_indices(l_intermediatefrequency) = add_ident ( 'intermediatefrequency' )
    lit_indices(l_isotopeRatio) =          add_ident ( 'isotopeRatio' )
    lit_indices(l_k) =                     add_ident ( 'K' )
    lit_indices(l_khz) =                   add_ident ( 'KHz' )
    lit_indices(l_km) =                    add_ident ( 'km' )
    lit_indices(l_linewidth) =             add_ident ( 'linewidth' )
    lit_indices(l_logp) =                  add_ident ( 'logp' )
    lit_indices(l_losVel) =                add_ident ( 'LOSVel' )
    lit_indices(l_lsbfrequency) =          add_ident ( 'LSBFrequency' )
    lit_indices(l_m) =                     add_ident ( 'm' )
    lit_indices(l_maf) =                   add_ident ( 'maf' )
    lit_indices(l_mafs) =                  add_ident ( 'mafs' )
    lit_indices(l_massMeanDiameterIce) =   add_ident ( 'massMeanDiameterIce' )
    lit_indices(l_massMeanDiameterWater) = add_ident ( 'massMeanDiameterWater' )
    lit_indices(l_mb) =                    add_ident ( 'mb' )
    lit_indices(l_meters) =                add_ident ( 'meters' )
    lit_indices(l_mhz) =                   add_ident ( 'MHz' )
    lit_indices(l_mif) =                   add_ident ( 'mif' )
    lit_indices(l_mifs) =                  add_ident ( 'mifs' )
    lit_indices(l_minutes) =               add_ident ( 'minutes' )
    lit_indices(l_none) =                  add_ident ( 'none' )
    lit_indices(l_orbitinclination) =      add_ident ( 'orbitInclination' )
    lit_indices(l_orbits) =                add_ident ( 'orbits' )
    lit_indices(l_pa) =                    add_ident ( 'pa' )
    lit_indices(l_ptan) =                  add_ident ( 'ptan' )
    lit_indices(l_ppbv) =                  add_ident ( 'ppbv' )
    lit_indices(l_ppmv) =                  add_ident ( 'ppmv' )
    lit_indices(l_pptv) =                  add_ident ( 'pptv' )
    lit_indices(l_rad) =                   add_ident ( 'rad' )
    lit_indices(l_radiance) =              add_ident ( 'radiance' )
    lit_indices(l_radians) =               add_ident ( 'radians' )
    lit_indices(l_refgph) =                add_ident ( 'refGPH' )
    lit_indices(l_s) =                     add_ident ( 's' )
    lit_indices(l_scanresidual) =          add_ident ( 'scanResidual' )
    lit_indices(l_scECI) =                 add_ident ( 'scECI' )
    lit_indices(l_scGeocAlt) =             add_ident ( 'scGeocAlt' )
    lit_indices(l_scvel) =                 add_ident ( 'scVel' )
    lit_indices(l_seconds) =               add_ident ( 'seconds' )
    lit_indices(l_sidebandratio) =         add_ident ( 'sidebandRatio' )
    lit_indices(l_spaceRadiance) =         add_ident ( 'spaceRadiance' )
    lit_indices(l_temperature) =           add_ident ( 'temperature' )
    lit_indices(l_temperature_prec) =      add_ident ( 'temperature_precision' )
    lit_indices(l_theta) =                 add_ident ( 'theta' )
    lit_indices(l_thz) =                   add_ident ( 'THz' )
    lit_indices(l_time) =                  add_ident ( 'time' )
    lit_indices(l_tngteci) =               add_ident ( 'tngteci' )
    lit_indices(l_tngtgeocalt) =           add_ident ( 'tngtgeocalt' )
    lit_indices(l_tngtgeodalt) =           add_ident ( 'tngtgeodalt' )
    lit_indices(l_totalOpticalDepth) =     add_ident ( 'totalOpticalDepth' )
    lit_indices(l_true) =                  add_ident ( 'true' )
    lit_indices(l_usbfrequency) =          add_ident ( 'USBFrequency')
    lit_indices(l_vmr) =                   add_ident ( 'vmr' )
    lit_indices(l_xyz) =                   add_ident ( 'xyz' )
    lit_indices(l_zeta) =                  add_ident ( 'zeta' )

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
