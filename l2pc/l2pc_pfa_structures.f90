module L2PC_PFA_STRUCTURES
  use L2PCDim, only: MAX_NO_PHI, NLVL
  use L2PC_File_Parameters, only: MAX_NO_BANDS, &
                                  MAX_NO_ELMNTS_PER_SV_COMPONENT, &
                                  MAX_NO_POINTINGS
  use MLSCommon, only: I4, R4, R8
  implicit NONE
!  This has the standard structure constructions used for
!    computing the data for the l2pc_xx.dat file.
!  Originally written by W. G. Read on 8/13/1990.
!  Modified, Mar/10/2000, Z. Shippony, to fit for EOS concept
!---------------------------- RCS Ident Info -------------------------------
  PRIVATE :: Id, ModuleName
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName = &
       "$RCSfile$"
!---------------------------------------------------------------------------
!_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
! Parameter declarations
  Integer(i4), parameter :: MAXSPS = 20
  Integer(i4), parameter :: MAXSUBSPS = 15
  Integer(i4), parameter :: MAXGEOM = 13
  Integer(i4), parameter :: MAXGEOPHYS = 2
! maxsps + maxgeom + maxgeophys + 4 (for pointing and 3 X field strength
! computations) should equal max_no_sv_components.
! More structures :
! These are for catagories of state vector types. This catagorization is
! because the method of differentiation differs according to catagory.
! This structure contains the a priori limb point pressure
  type :: LIMB_PRESS
    character(len=8) :: NAME
    Integer(i4) :: NO_LIN_VALUES
    Logical*1 DER_CALC(max_no_bands)
    Real(r4) :: LIN_VAL(max_no_pointings)
  end type LIMB_PRESS
! This structure contains the linearization values for magnetic field strength.
  type :: MAG_FIELD
    character(len=8) :: NAMES(3)
    Logical*1 DER_CALC(max_no_bands)
    Real(r4) :: TOTAL_MAG_FIELD
    Real(r4) :: PROP_DIR_ANGLE
    Real(r4) :: POLAR_ANGLE
  end type MAG_FIELD
! This structure contains geometric and any other parameters characterized
! with a single value
  type :: GEOM_PARAM
    character(len=8) :: NAME
    Logical*1 DER_CALC(max_no_bands)
    Real(r4) :: LIN_VAL
  end type GEOM_PARAM
! This structure contains spectroscopic information about the species
  type :: SPECTRO_PARAM
    character(len=1) :: TYPE
    character(len=8) :: NAME
    Integer(i4) :: SPECTAG
    Integer(i4) :: NO_PHI_VALUES
    Integer(i4) :: NO_ZETA_VALUES
    Logical*1 DER_CALC(max_no_bands)
    Real(r4) :: PHI_BASIS(max_no_phi+2)
    Real(r4) :: ZETA_BASIS(max_no_elmnts_per_sv_component+2)
  end type SPECTRO_PARAM
! This structure contains geophysical information about the atmosphere
  type :: GEOPHYS_PARAM
    character(len=8) :: NAME
    Integer(i4) :: NO_LIN_VALUES
    Logical*1 DER_CALC(max_no_bands)
    Real(r4) :: LIN_VAL(max_no_elmnts_per_sv_component)
    Real(r4) :: BASIS_PEAKS(max_no_elmnts_per_sv_component)
  end type GEOPHYS_PARAM
! This structure contains atmospheric composition of the atmosphere
  type :: ATMOS_COMP
    character(len=8) :: NAME
    Integer(i4) :: SPECTAG
    Integer(i4) :: NO_LIN_VALUES
    Logical*1 FWD_CALC(max_no_bands)
    Logical*1 DER_CALC(max_no_bands)
    Real(r4) :: LIN_VAL(max_no_elmnts_per_sv_component)
    Real(r4) :: BASIS_PEAKS(max_no_elmnts_per_sv_component+2)
  end type ATMOS_COMP
! The cross section structure
  type :: CS_DATA
    Integer(i4) :: SPECTAG
    Integer(i4) :: NO_HTS
    Real(r8) FREQ
    Real(r4) :: ABS_CS(nlvl)
    Real(r4) :: DABS_CS_DF(nlvl)
    Real(r4) :: TEMP_DEP(nlvl)
  end type CS_DATA
!-----------------------------------------------------------------------
! Physical Constants
  Real(r4), parameter :: EARTH_MINOR = 6356.755
  Real(r4), parameter :: EARTH_MAJOR = 6378.140
!-----------------------------------------------------------------------
!  NEW PFA parameters file, Z. Shippony  Jun/5/92
!  Modified, Aug/18/92, Z. Shippony (Add new type of structure for new
!  "PQM" approach)
!  Modified, Mar/12/00, Z. Shippony (Set maxaitkenpts back to 30)
!  Modified, Jul/03/97, Z. Shippony (Add pressure shift parameter)
!  Modified, Jun/21/97, Z. Shippony (Up maxlines from: 30  to: 35)
!  Modified, Jul/28/95, Z. Shippony (Up maxpfach from: 25  to: 45)
!  Modified, Apr/08/94, Z. Shippony (Up maxaitkenpts from: 25  to: 201)
  Integer(i4), parameter :: MAXPFALINES = 6
  Integer(i4), parameter :: MAXPFACH = 45
  Integer(i4), parameter :: MAXFILTPTS = 161
  Integer(i4), parameter :: MAXAITKENPTS = 30
  Integer(i4), parameter :: MAXLINES = 35
  Integer(i4), parameter :: MAXRAT = 11
  type :: PFA_SLAB
    Integer(i4) :: NO_SPS
    Integer(i4) :: NO_LINES
    Integer(i4) :: NRAT(nlvl)
    Integer(i4) :: SPS_SPECTAG
    Character(len=8) :: SPS_NAME
    Real(r8) :: VARM(nlvl)
    Real(r8) :: SPS_V0(maxlines)
    Real(r8) :: SPS_EL(maxlines)
    Real(r8) :: SPS_STR(maxlines)
    Real(r8) :: SPS_W(maxlines)
    Real(r8) :: SPS_PS(maxlines)
    Real(r8) :: SPS_N(maxlines)
    Real(r8) :: SPS_N1(maxlines)
    Real(r8) :: SPS_N2(maxlines)
    Real(r8) :: SPS_GAMMA(maxlines)
    Real(r8) :: SPS_DELTA(maxlines)
    Real(r8) :: SPS_PART(3,maxlines)
    Real(r8) :: XX(maxrat,nlvl)
    Real(r8) :: YY(maxrat,nlvl)
    Real(r8) :: DY(maxrat,nlvl)
  end type PFA_SLAB
end module L2PC_PFA_STRUCTURES
! $Log$
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
!
