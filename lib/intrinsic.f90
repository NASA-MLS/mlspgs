! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module INTRINSIC

! Intrinsic constants needed by Init_Tables_Module, DeclarationTable, etc.

! Declaring the definitions is handled by the tree walker.

  use MOLECULES ! everything, in particular FIRST_MOLECULE, INIT_MOLECULES,
  !               and LAST_MOLECULE.  There is no "only" clause so as to
  !               make all of the literals, e.g. l_h2o, available here, too.
  use SYMBOL_TABLE, only: ENTER_TERMINAL
  use SYMBOL_TYPES, only: T_IDENTIFIER

  implicit NONE
  public

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

! Data types that don't have enumerated literals:
  integer, public, parameter :: T_FIRST          = 1
  integer, public, parameter :: T_NUMERIC        = t_first
  integer, public, parameter :: T_NUMERIC_RANGE  = t_numeric + 1
  integer, public, parameter :: T_STRING         = t_numeric_range + 1
! Enumeration types:
  integer, public, parameter :: T_BOOLEAN        = t_string + 1
  integer, public, parameter :: T_LAST_INTRINSIC = t_boolean

! Abstract physical quantities:
  integer, public, parameter :: PHYQ_INVALID = 0 ! Invalid unit given by user
  integer, public, parameter :: PHYQ_DIMENSIONLESS = 1 ! Dimensionless quantity
  integer, public, parameter :: PHYQ_LENGTH = 2        ! Default meters
  integer, public, parameter :: PHYQ_TIME = 3          ! Default seconds
  integer, public, parameter :: PHYQ_PRESSURE = 4      ! Default millibars
  integer, public, parameter :: PHYQ_TEMPERATURE = 5   ! Default Kelvins
  integer, public, parameter :: PHYQ_VMR = 6           ! Default parts-per-one
  integer, public, parameter :: PHYQ_ANGLE = 7         ! Default degrees
  integer, public, parameter :: PHYQ_MAFS = 8          ! Default MAFs
  integer, public, parameter :: PHYQ_MIFS = 9          ! Default MIFs
  integer, public, parameter :: PHYQ_FREQUENCY = 10    ! Default MHz
  integer, public, parameter :: PHYQ_ZETA = 11         ! log10(pressure/hPa)
  integer, public, parameter :: FIRST_PHYQ = phyq_invalid, LAST_PHYQ = phyq_zeta
  integer, public :: PHYQ_INDICES(first_phyq:last_phyq)

! Enumeration literals:
  integer, public, parameter :: FIRST_LIT       = first_molecule
  integer, public, parameter :: L_BASELINE      = last_molecule+1
  integer, public, parameter :: L_DAYS          = l_baseline + 1
  integer, public, parameter :: L_DEG           = l_days + 1
  integer, public, parameter :: L_DEGREES       = l_deg + 1
  integer, public, parameter :: L_DIMENSIONLESS = l_degrees + 1
  integer, public, parameter :: L_DIMLESS       = l_dimensionless + 1
  integer, public, parameter :: L_DL            = l_dimless + 1
  integer, public, parameter :: L_EXTINCTION    = l_dl + 1
  integer, public, parameter :: L_FALSE         = l_extinction + 1
  integer, public, parameter :: L_GEODALTITUDE  = l_false + 1
  integer, public, parameter :: L_GHZ           = l_geodaltitude + 1
  integer, public, parameter :: L_GPH           = l_ghz + 1
  integer, public, parameter :: L_GPH_PRECISION = l_gph + 1
  integer, public, parameter :: L_HOURS         = l_gph_precision + 1
  integer, public, parameter :: L_HPA           = l_hours + 1
  integer, public, parameter :: L_HZ            = l_hpa + 1
  integer, public, parameter :: L_K             = l_hz + 1
  integer, public, parameter :: L_KHZ           = l_k  + 1
  integer, public, parameter :: L_KM            = l_khz + 1
  integer, public, parameter :: L_LINEWIDTH     = l_km + 1
  integer, public, parameter :: L_LOGP          = l_linewidth + 1
  integer, public, parameter :: L_M             = l_logp + 1
  integer, public, parameter :: L_MAF           = l_m + 1
  integer, public, parameter :: L_MAFS          = l_maf + 1
  integer, public, parameter :: L_MB            = l_mafs + 1
  integer, public, parameter :: L_METERS        = l_mb + 1
  integer, public, parameter :: L_MHZ           = l_meters  + 1
  integer, public, parameter :: L_MIF           = l_mhz + 1
  integer, public, parameter :: L_MIFS          = l_mif + 1
  integer, public, parameter :: L_MINUTES       = l_mifs + 1
  integer, public, parameter :: L_ORBITINCLINATION = l_minutes + 1
  integer, public, parameter :: L_ORBITS        = l_orbitinclination + 1
  integer, public, parameter :: L_PA            = l_orbits + 1
  integer, public, parameter :: L_PPBV          = l_pa + 1
  integer, public, parameter :: L_PPMV          = l_ppbv + 1
  integer, public, parameter :: L_PPTV          = l_ppmv + 1
  integer, public, parameter :: L_PTAN          = l_pptv + 1
  integer, public, parameter :: L_RAD           = l_ptan + 1
  integer, public, parameter :: L_RADIANCE      = l_rad + 1
  integer, public, parameter :: L_RADIANS       = l_radiance + 1
  integer, public, parameter :: L_REFGPH        = l_radians + 1
  integer, public, parameter :: L_S             = l_refgph + 1
  integer, public, parameter :: L_SCANRESIDUAL  = l_s + 1
  integer, public, parameter :: L_SCVEL         = l_scanresidual + 1
  integer, public, parameter :: L_SECONDS       = l_scvel + 1
  integer, public, parameter :: L_SIDEBANDRATIO = l_seconds + 1
  integer, public, parameter :: L_TANGENTHEIGHT = l_sidebandratio + 1
  integer, public, parameter :: L_TEMPERATURE   = l_tangentheight + 1
  integer, public, parameter :: L_TEMPERATURE_PREC = l_temperature + 1
  integer, public, parameter :: L_THZ           = l_temperature_prec + 1
  integer, public, parameter :: L_TRUE          = l_thz + 1
  integer, public, parameter :: L_VMR           = l_true + 1
  integer, public, parameter :: L_ZETA          = l_vmr + 1
  integer, public, parameter :: LAST_INTRINSIC_LIT = l_zeta

contains ! =====     Public procedures     =============================
! -----------------------------------------------  INIT_INTRINSIC  -----
  subroutine INIT_INTRINSIC ( DATA_TYPE_INDICES, LIT_INDICES )
    integer, intent(inout) :: DATA_TYPE_INDICES(:)
    integer, intent(inout) :: LIT_INDICES(:)

  ! Put molecules into the symbol table.

    call init_molecules ( lit_indices )

  ! Put intrinsic predefined identifiers into the symbol table.

    ! Put intrinsic non-enumeration type names into symbol table
    data_type_indices(t_numeric) =         add_ident ( 'numeric' )
    data_type_indices(t_numeric_range) =   add_ident ( 'numeric_range' )
    data_type_indices(t_string) =          add_ident ( 'string' )
    ! Put intrinsic enumeration type names into the symbol table
    data_type_indices(t_boolean) =         add_ident ( 'boolean' )
    ! Put intrinsic enumeration literals into the symbol table:
    lit_indices(l_baseline) =              add_ident ( 'baseline' )
    lit_indices(l_days) =                  add_ident ( 'days' )
    lit_indices(l_deg) =                   add_ident ( 'deg' )
    lit_indices(l_degrees) =               add_ident ( 'degrees' )
    lit_indices(l_dimensionless) =         add_ident ( 'dimensionless' )
    lit_indices(l_dimless) =               add_ident ( 'dimless' )
    lit_indices(l_dl) =                    add_ident ( 'dl' )
    lit_indices(l_extinction) =            add_ident ( 'extinction' )
    lit_indices(l_false) =                 add_ident ( 'false' )
    lit_indices(l_geodaltitude) =          add_ident ( 'geodAltitude' )
    lit_indices(l_ghz) =                   add_ident ( 'GHz' )
    lit_indices(l_gph) =                   add_ident ( 'gph' )
    lit_indices(l_gph) =                   add_ident ( 'gph_precision' )
    lit_indices(l_hours) =                 add_ident ( 'hours' )
    lit_indices(l_hpa) =                   add_ident ( 'hPa' )
    lit_indices(l_hz) =                    add_ident ( 'Hz' )
    lit_indices(l_k) =                     add_ident ( 'k' )
    lit_indices(l_khz) =                   add_ident ( 'KHz' )
    lit_indices(l_km) =                    add_ident ( 'km' )
    lit_indices(l_linewidth) =             add_ident ( 'linewidth' )
    lit_indices(l_logp) =                  add_ident ( 'logp' )
    lit_indices(l_m) =                     add_ident ( 'm' )
    lit_indices(l_maf) =                   add_ident ( 'maf' )
    lit_indices(l_mafs) =                  add_ident ( 'mafs' )
    lit_indices(l_mb) =                    add_ident ( 'mb' )
    lit_indices(l_meters) =                add_ident ( 'meters' )
    lit_indices(l_mhz) =                   add_ident ( 'MHz' )
    lit_indices(l_mif) =                   add_ident ( 'mif' )
    lit_indices(l_mifs) =                  add_ident ( 'mifs' )
    lit_indices(l_minutes) =               add_ident ( 'minutes' )
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
    lit_indices(l_scvel) =                 add_ident ( 'scVel' )
    lit_indices(l_seconds) =               add_ident ( 'seconds' )
    lit_indices(l_sidebandratio) =         add_ident ( 'sidebandRatio' )
    lit_indices(l_tangentheight) =         add_ident ( 'tangentHeight' )
    lit_indices(l_temperature) =           add_ident ( 'temperature' )
    lit_indices(l_temperature) =           add_ident ( 'temperature_prec' )
    lit_indices(l_thz) =                   add_ident ( 'THz' )
    lit_indices(l_true) =                  add_ident ( 'true' )
    lit_indices(l_vmr) =                   add_ident ( 'vmr' )
    lit_indices(l_zeta) =                  add_ident ( 'zeta' )

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
    phyq_indices(phyq_zeta) =              add_ident ( 'zeta' )

  contains

    integer function ADD_IDENT ( TEXT )
      character(len=*), intent(in) :: TEXT
      add_ident = enter_terminal ( text, t_identifier )
    end function ADD_IDENT

  end subroutine INIT_INTRINSIC

end module INTRINSIC

! $Log$
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
