module MDBETA
  use MLSCommon, only: I4, R4, R8
  implicit NONE
  private
  public :: MAX_NO_ZETA, MAX_NO_FREQ, NO_T_PHI, MDB_BETA_REC
!  This contains all the parameters and structure definitions
!  needed by the EOS Master DataBase routines
!  Written by: Z. Shippony, June/24/1998
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
  Integer(i4), parameter :: MAX_NO_ZETA = 337
  Integer(i4), parameter :: MAX_NO_FREQ= 15
  Integer(i4), parameter :: NO_T_PHI = 5
  type :: MDB_BETA_REC
    Integer(i4) :: NO_FREQ
    Integer(i4) :: SPECTAG
    Real(r8) :: FREQ_P(max_no_freq)
    Real(r8) :: FREQ_I(max_no_freq)
    Real(r8) :: BETA_P(max_no_zeta,no_t_phi,max_no_freq)
    Real(r8) :: BETA_I(max_no_zeta,no_t_phi,max_no_freq)
  end type MDB_BETA_REC
end module MDBETA
! $Log$
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
!
