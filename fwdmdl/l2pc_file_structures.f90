module L2PC_FILE_STRUCTURES
  use L2PC_File_Parameters ! everything
  use MLSCommon, only: I4, R4
  implicit NONE

!  This has the standard structure constructions used for
!    reading and writing the l2pc_xx.dat file.
!  Written by W. G. Read on 8/13/90.

!  Modified by W. G. Read 3/09/90:
!    Includes new structure for use in reading pseudo key access file

!  Modified by Z. Shippony, Jul/10/96:
!    Update to ver 5.0 of l2pc code. Add l2pc_header_tri  structure

!  Last Modified by Z. Shippony, Oct/10/96:
!    Update to ver 5.0 of l2pc code.

!  Last Modified by Z. Shippony, May/5/97:
!    Update to ver 5.1 of l2pc code.

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

!---------------------------- RCS Ident Info -------------------------------
  PRIVATE :: Id, ModuleName
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------

! Header structures

  type :: L2PC_HEADER_ONE
    Character(len=80) :: LINE1
    Character(len=80) :: LINE2
    Character(len=80) :: LINE3
    Character(len=8) :: AVAIL_KEYS(max_no_avail_keys)
    Character(len=8) :: SV_COMPONENTS(max_no_sv_components)
    Logical*1 SV_RTRVL_BY_BAND(max_no_sv_components,max_no_bands)
    Integer(i4) :: NO_ELMNTS_PER_SV_COMPONENT(max_no_sv_components)
    Integer(i4) :: SV_COMPONENT_FIRST_ELMNT_INDEX(max_no_sv_components)
    Integer(i4) :: NO_BANDS
    Integer(i4) :: NO_CHANNELS_PER_BAND
    Integer(i4) :: NO_POINTINGS
    Integer(i4) :: NO_AVAIL_KEYS
    Integer(i4) :: NO_SV_COMPONENTS
    Integer(i4) :: NO_COEFF_PER_COMPONENT
    Integer(i4) :: NO_K_RECORDS_PER_BIN
    Integer(i4) :: NO_MAG_FIELDS
    Integer(i4) :: NO_B_THETA_LIN_VAL
    Real(r4) :: POINTINGS(max_no_pointings)
    Real(r4) :: B_FIELDS(max_no_mag_fields)
    Real(r4) :: B_PHI_LIN_VAL
    Real(r4) :: B_THETA_LIN_VAL(max_no_theta_val)
  end type L2PC_HEADER_ONE

!**********  Size of(l2pc_header_one) = 1676 bytes (419 words)

  type :: L2PC_HEADER_TWO
    Integer(i4) :: NO_SV_ELMNTS
    Real(r4) :: TRI_BASIS_VERT_GRID(max_no_sv_elmnts)
  end type L2PC_HEADER_TWO

!**********  Size of(l2pc_header_two) = 2316 bytes (579 words)

  type :: L2PC_HEADER_TRI
    Integer(i4) :: MATDIM
    Character(len=1) :: SECOND_DER_MATRIX_BANDS(max_table_2d)
    Character(len=2) :: SECOND_DER_MATRIX_NAMID(max_table_2d)
  end type L2PC_HEADER_TRI

!**********  Size of(l2pc_header_tri) = 1588 bytes (397 words)

! Spectral power (Radiances) structure

  type :: I_STAR
    Integer(i4) :: FIRST_CHANNEL_NUMBER
    Real(r4) :: RADIANCES(max_no_pointings,max_no_channels_per_band)
  end type I_STAR

!**********  Size of(i_star) = 2588 bytes ( 647 words)

! State vector structure

  type :: X_STAR
    Integer(i4) :: NO_POINTINGS
    Integer(i4) :: NO_SV_ELMNTS
    Real(r4) :: POINTINGS(max_no_pointings)
    Real(r4) :: SV_ELMNTS(max_no_sv_elmnts)
  end type X_STAR

!**********  Size of(x_star) = 4024 bytes (1006 words)

! The derivative of the forward model with respect to a state vector element

  type :: K_STAR
    Character(len=8) :: SV_NAME1
    Character(len=8) :: SV_NAME2
    Real(r4) :: LOG_PSV1
    Real(r4) :: LOG_PSV2
    Integer(i4) :: FIRST_CHANNEL_NUMBER
    Real(r4) :: DERIVATIVES(max_no_pointings,max_no_channels_per_band)
  end type K_STAR

!**********  Size of(k_star) = 4048 bytes (1012 words)

! This Structure is used in collating the k_star records

  type :: DUMP_K_STAR
    Character(len=40) :: KEY
    Character(len=8) :: SV_NAME1
    Character(len=8) :: SV_NAME2
    Real(r4) :: LOG_PSV1
    Real(r4) :: LOG_PSV2
    Integer(i4) :: FIRST_CHANNEL_NUMBER
    Real(r4) :: DERIVATIVES(max_no_pointings)
  end type DUMP_K_STAR

!**********  Size of(dump_k_star) = 336 bytes (84 words)

! This Structure goes into a companion file which is used to find keys

  type :: L2PC_KEYS
    Character(len=32) :: L2PC_KEY
    Integer(i4) :: REC_NO
  end type L2PC_KEYS

!**********  Size of(l2pc_keys) = 36 bytes (9 words)

end module L2PC_FILE_STRUCTURES

! $Log$
