module L2PC_FILE_PARAMETERS
  use MLSCommon, only: I4, R4
  use L2PCDIM, only: NSPS
  implicit NONE
  public

!---------------------------- RCS Ident Info -------------------------------
  private ID, ModuleName
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------

! * * * * * * * * * * *  L2PC_FILE_PARAMETERS V51 * * * * * * * * * * * * *
!
!  This include file has the parameters values for
!    reading and writing the l2pc_xx.dat file.
!
!  Written by W. G. Read  on 8/13/90.
!
!  Modified to Version 5.0 of l2pc code, by Z. Shippony, Jul/10/96
!  Adding the second derivative table dimension Parameters and the
!  third header's key
!
!  Updated by Z. Shippony on May/5/97.  (For version 5.1)
!      setting    max_no_sv_derivatives to: 350
!    (Increasing  max_no_elmnts_per_sv_component from 22 to: 43)
!     Increasing  max_no_elmnts_per_sv_component from 22 to: 34
!    (Increasing  max_no_pointings from: 43  to: 85)
!     Increasing  max_no_pointings from: 43  to: 67
!
!  Updated by Z. Shippony on Oct/16/97
!     resetting   max_no_sv_derivatives to: 400  ON SGI ONLY !!
!
! Declarations for Parameters.
!
! If these parameters are changed be sure to modify the record length
! in the l2pc file OPEN statement to match the new length. Keep
! max_char_len_l2pc_key divisable by 4

  integer(i4), parameter :: MAX_CHAR_LEN_L2PC_KEY = 4*10
  integer(i4), parameter :: MAX_CHAR_LEN_AVAIL_KEYS = 8
  integer(i4), parameter :: MAX_CHAR_LEN_SV_COMPONENTS = 8
  integer(i4), parameter :: MAX_NO_AVAIL_KEYS = 48
  integer(i4), parameter :: MAX_NO_BANDS = 6
  integer(i4), parameter :: MAX_NO_CHANNELS_PER_BAND = 15
! integer(i4), parameter :: MAX_NO_POINTINGS = 85
  integer(i4), parameter :: MAX_NO_POINTINGS = 67
! integer(i4), parameter :: MAX_NO_ELMNTS_PER_SV_COMPONENT = 43
  integer(i4), parameter :: MAX_NO_ELMNTS_PER_SV_COMPONENT = 34
! integer(i4), parameter :: MAX_NO_KEY_ADDR = 100000
  integer(i4), parameter :: MAX_NO_KEY_ADDR = 5000
  integer(i4), parameter :: MAX_NO_THETA_VAL = 6
  integer(i4), parameter :: MAX_NO_MAG_FIELDS = 4
  integer(i4), parameter :: MAX_NO_SV_COMPONENTS = 32

  integer(i4), parameter :: MAX_NO_SV_ELMNTS =                          &
     &                       MAX_NO_ELMNTS_PER_SV_COMPONENT * NSPS

! integer(i4), parameter :: MAX_NO_SV_DERIVATIVES = 350
! integer(i4), parameter :: MAX_NO_SV_DERIVATIVES = 385
  integer(i4), parameter :: MAX_NO_SV_DERIVATIVES = 400

  integer(i4), parameter :: MAX_CHAR_LEN_L2PC_KEY_DES =                 &
 &                            max_char_len_l2pc_key - max_char_len_avail_keys
  integer(i4), parameter :: MAX_NO_CHANNELS =                           &
 &                            max_no_bands * max_no_channels_per_band
  integer(i4), parameter :: MAX_REC_LEN_KEY_FILE =                      &
 &                            1 + max_char_len_l2pc_key_des/4
  character(len=40), parameter :: L2PC_HEADER_KEY1 =                    &
 &                                 'First header record.    '
  character(len=40), parameter :: L2PC_HEADER_KEY2 =                    &
 &                                 'Second header record.   '
  character(len=40), parameter :: L2PC_HEADER_KEY3 =                    &
 &                                 'Third header record.    '

  integer(i4), parameter :: J2D=max_no_sv_components,                   &
 &                          MAX_TABLE_2D=(j2d*(j2d+1))/2

  ! Conversion factor degrees to radians.
  real(r4), parameter :: DEG2RAD = 1.745329252E-2 ! [Radians / Degree]

  ! Conversion factor radians to degrees.
  real(r4), parameter :: RAD2DEG = 57.29577951308 ! [Degrees / Radian]
!
end module L2PC_FILE_PARAMETERS

! $Log$
