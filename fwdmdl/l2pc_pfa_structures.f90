module L2PC_PFA_STRUCTURES
  use MLSCommon, only: I4, R4, R8
  use L2PC_File_Parameters, only: MAX_NO_ELMNTS_PER_SV_COMPONENT, &
                                  MAX_NO_POINTINGS
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
  type LIMB_PRESS
    character(len=8) :: NAME
    Integer(i4) :: NO_LIN_VALUES
    Logical DER_CALC(6)
    Real(r4) :: LIN_VAL(max_no_pointings)
  end type LIMB_PRESS

! This structure contains geometric and any other parameters characterized
! with a single value
  type GEOM_PARAM
    character(len=8) :: NAME
    Logical DER_CALC(6)
    Real(r4) :: LIN_VAL
  end type GEOM_PARAM

! This structure contains spectroscopic information about the species
  type SPECTRO_PARAM
    character(len=1) :: TYPE
    character(len=8) :: NAME
    Integer(i4) :: SPECTAG
    Integer(i4) :: NO_PHI_VALUES
    Integer(i4) :: NO_ZETA_VALUES
    Logical DER_CALC(6)
    Real(r4) :: PHI_BASIS(12)
    Real(r4) :: ZETA_BASIS(52)
  end type SPECTRO_PARAM

! This structure contains geophysical information about the atmosphere
  type GEOPHYS_PARAM
    character(len=8) :: NAME
    Integer(i4) :: NO_LIN_VALUES
    Logical DER_CALC(6)
    Real(r4) :: LIN_VAL(50)
    Real(r4) :: BASIS_PEAKS(52)
  end type GEOPHYS_PARAM

! This structure contains atmospheric composition of the atmosphere
  type ATMOS_COMP
    character(len=8) :: NAME
    Integer(i4) :: SPECTAG
    Integer(i4) :: NO_LIN_VALUES
    Logical FWD_CALC(6)
    Logical DER_CALC(6)
    Real(r4) :: LIN_VAL(50)
    Real(r4) :: BASIS_PEAKS(52)
  end type ATMOS_COMP

! This NEW structure contains info about the K matrix
  type K_MATRIX_INFO
    character(len=8) :: NAME
    Integer(i4) :: FIRST_DIM_INDEX
    Integer(i4) :: NO_ZETA_BASIS
    Integer(i4) :: NO_PHI_BASIS
    Real(r4)    :: PHI_BASIS(12)
    Real(r4)    :: ZETA_BASIS(52)
  end type K_MATRIX_INFO

!-----------------------------------------------------------------------
! Physical Constants
  Real(r4), parameter :: EARTH_MINOR = 6356.755
  Real(r4), parameter :: EARTH_MAJOR = 6378.140
!-----------------------------------------------------------------------
!  NEW PFA parameters file, Z. Shippony  Jun/5/92
!  Modified, Aug/18/92, Z. Shippony (Add new type of structure for new
!  "PQM" approach)
!  Modified, Jul/03/97, Z. Shippony (Add pressure shift parameter)
!  Modified, Jun/21/97, Z. Shippony (Up maxlines from: 30  to: 35)
!  Modified, Jul/28/95, Z. Shippony (Up maxpfach from: 25  to: 45)
  Integer(i4), parameter :: MAXPFALINES = 6
  Integer(i4), parameter :: MAXPFACH = 45
  Integer(i4), parameter :: MAXFILTPTS = 161
  Integer(i4), parameter :: MAXLINES = 35

  type :: PFA_SLAB
    Integer(i4) :: NO_SPS
    Integer(i4) :: NO_LINES
    Integer(i4) :: SPS_SPECTAG
    Character(len=8) :: SPS_NAME
    Real(r8) :: SPS_QLOG(3)
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
  end type PFA_SLAB

end module L2PC_PFA_STRUCTURES
! $Log$
! Revision 1.7  2001/02/19 22:20:40  zvi
! Latest modification: Conv/NoConv
!
! Revision 1.6  2001/02/19 22:14:21  zvi
!
! Revision 1.1  2000/06/21 21:56:15  zvi
! First version D.P.
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
