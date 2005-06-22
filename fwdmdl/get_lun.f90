! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

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
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains 
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module GET_LUN
! $Log$
! Revision 2.1  2002/10/08 17:08:04  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model
!
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
