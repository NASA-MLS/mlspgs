module EOS_MDB
  use MLSCommon, only: I4, R4, R8
  implicit NONE
  public
! This file contains all the parameters and structure definitions
! needed by the EOS Master DataBase routines
! Written  by: Z. Shippony, Aug/18/1999
! Modified by: Z. Shippony, Oct/12/1999
!---------------------------- RCS Ident Info -------------------------------
  private ID, ModuleName
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
!
  Real(r4), parameter :: TEMP_LO = 170.0  ! Lower limit of temperature
  Real(r4), parameter :: TTEMP_HI = 320.0 ! Upper limit of temperature
  Real(r4), parameter :: ZETA_LO = -3.0   ! Lower limit of Zeta
  Real(r4), parameter :: ZETA_HI = 4.0    ! Upper limit of Zeta
!
  Integer(i4), parameter :: NDEC=7        ! Number of Zeta decades
  Integer(i4), parameter :: NZPD=12       ! Number of Zeta points (subdivision)
                                          ! per decade
! Integer(i4), parameter :: MAX_ZETA = ndec * nzpd * (Ng + 1) + 1)
  Integer(i4), parameter :: MAX_ZETA = ndec * nzpd + 1
  Integer(i4), parameter :: MAX_TEMP=15, MAX_NO_LINES=20, MAX_FREQ=30
!
  Type EOS_MDB_HDR
    Integer(i4) :: Spectag
    Integer(i4) :: no_lines
    Integer(i4) :: no_f_grid(max_no_lines)
    Real(r8) :: el(max_no_lines),log_i(max_no_lines),n(max_no_lines),   &
 &         w(max_no_lines),delta(max_no_lines),n1(max_no_lines),        &
 &         n2(max_no_lines),gamma(max_no_lines),v0(max_no_lines),       &
 &         ps(max_no_lines),q_log(3)
    Real(r4) :: Zeta(max_zeta)
    Real(r4) :: Log_Temp(max_temp)
    Real(r8) :: x_grid(max_freq,max_no_lines)
  End Type EOS_MDB_HDR
!
  Type EOS_MDB_REC
    Real(r4) :: Log_beta(max_zeta,max_temp,max_freq)
    Real(r4) :: dLog_beta_dw(max_zeta,max_temp,max_freq)
    Real(r4) :: dLog_beta_dn(max_zeta,max_temp,max_freq)
    Real(r4) :: dLog_beta_dNu0(max_zeta,max_temp,max_freq)
  End type EOS_MDB_REC
!
end module EOS_MDB
! $Log$
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
!
