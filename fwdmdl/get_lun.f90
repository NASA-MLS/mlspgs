! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module GET_LUN
  implicit NONE
  public
!  This is used with sid_fwd_mdl routines and uth routines.
!  It contains parameters that are the unit numbers for various files opened
!  by both programs, and it centralized the allocation so the above programs
!  will not contradict each other or write/read over each other's files.
!  Written by: Z. Shippony,    Aug/12/1997
  integer, parameter :: mdb_unit_dat=37     ! master database unit number
  integer, parameter :: mdb_unit_ndx=38     ! master database index unit number
  integer, parameter :: mdb_pqm_unit_dat=33 ! master pqm database unit number
  integer, parameter :: mdb_pqm_unit_ndx=34 ! master pqm database index unit num
  integer, parameter :: p_vs_h_unit=27      ! p_vs_h file unit number
  integer, parameter :: filter_unit=47      ! filter shapes file(s) unit number
  integer, parameter :: aaap_unit=72        ! Antenna file unit number
  integer, parameter :: cs_tbl_unit=31      ! ascii cs table file(s) unit number
  integer, parameter :: ascii_pqm_unit=17   ! ascii pqm file(s) unit number
  integer, parameter :: disp_pqm_unit=43    ! ascii disposition pqm file(s) unit
  integer, parameter :: velcor_unit=19      ! Velocity correction file unit
  integer, parameter :: model_unit=62       ! climatology file(s) unit
!---------------------------- RCS Ident Info -------------------------------
  private ID, ModuleName
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
end module GET_LUN
! $Log$
! Revision 1.4  2001/06/07 23:30:34  pwagner
! Added Copyright statement
!
! Revision 1.3  2001/03/05 23:35:30  zvi
! Recovering get_lun.f90 from Paul
!
! Revision 1.1  2001/01/31 23:57:43  zvi
! Various changes..
!
! Revision 1.4  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.1  2000/05/04 18:12:06  vsnyder
! Initial conversion to Fortran 90
!
