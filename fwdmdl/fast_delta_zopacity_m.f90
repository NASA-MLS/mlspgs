! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module FAST_DELTA_ZOPACITY_M
  use MLSCommon, only: I4, R8
  use ELLIPSE_M, only: ELLIPSE
  use RAD_DELTA_INTEGRAL_M, only: RAD_DELTA_INTEGRAL
  use PATH_ENTITIES_M, only: PATH_VECTOR, PATH_BETA
  implicit NONE
  private
  public :: FAST_DELTA_ZOPACITY
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
  "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
  "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------
!  This routine computes the d_delta function. Integration is done
!  using the Gauss-Legendre method.
!
!  ** NOTE: This routine integrate in ZETA Space !
!
  Subroutine FAST_DELTA_ZOPACITY(mid,brkpt,no_ele,z_path,h_path,     &
 &           beta_path,dHdz_path,spsfunc_path,n_sps,N_lvls,ref_corr, &
 &           elvar,midval_ndx,no_midval_ndx,gl_ndx,no_gl_ndx,        &
 &           midval_delta,del_opacity,Ier)
!
    Integer(i4), intent(in) :: N_SPS, N_LVLS, MID, BRKPT, NO_ELE
    Integer(i4), intent(in) :: midval_ndx(:,:), no_midval_ndx
    Integer(i4), intent(in) :: gl_ndx(:,:), no_gl_ndx

    Real(r8), intent(in) :: REF_CORR(:)
    Real(r8), intent(in) :: midval_delta(:,:)    ! (N2lvl,Nsps)

    Real(r8), intent(out) :: DEL_OPACITY(:)

    Integer(i4), intent(out) :: IER

    Type(ELLIPSE), intent(in out) :: elvar

    Type(path_beta), intent(in) :: BETA_PATH(:)      ! (Nsps)

    Type(path_vector), intent(in) :: Z_PATH, H_PATH, DHDZ_PATH
    Type(path_vector), intent(in) :: SPSFUNC_PATH(:)
!
! -----     Local variables     ----------------------------------------
!
    Integer(i4) :: I, J
!
    Real(r8), ALLOCATABLE, DIMENSION(:) :: Integrand
!
! -----     Executable statements     ----------------------------------
!
    Ier = 0
    DEALLOCATE(Integrand, STAT=i)
!
    ALLOCATE(Integrand(no_ele), STAT=ier)
    IF(ier /= 0) THEN
      Print *,'** ALLOCATION error in FAST_DELTA_ZOPACITY routine ..'
      Print *,'   Allocation STAT =',ier
      goto 99
    endif

    integrand(1:no_ele) = 0.0_r8
!
    do j = 1, n_sps
      integrand(1:no_ele) = integrand(1:no_ele)              +  &
     &                      spsfunc_path(j)%values(1:no_ele) *  &
     &                      beta_path(j)%values(1:no_ele)
    end do
!
!  Integrate:
!
    Call rad_delta_integral(mid,brkpt,no_ele,z_path,h_path,dhdz_path,   &
      &  N_lvls,ref_corr,integrand,elvar,gl_ndx,no_gl_ndx,del_opacity,Ier)
    IF(ier /= 0) goto 99
!
! Now complete the del_opacity for the case of midval delta
!
    if(no_midval_ndx > 0) then
      do j = 1, no_midval_ndx
        i = midval_ndx(j,1)
        del_opacity(i) = SUM(midval_delta(i,1:n_sps))
      end do
    endif

 99  DEALLOCATE(Integrand, STAT=i)

    Return
!
  End Subroutine FAST_DELTA_ZOPACITY
!
End module FAST_DELTA_ZOPACITY_M
! $Log$
! Revision 1.1  2001/06/21 13:07:09  zvi
! Speed enhancement MAJOR update
!
! Revision 1.6  2001/06/07 23:30:34  pwagner
! Added Copyright statement
!
! Revision 1.5  2001/03/31 23:40:55  zvi
! Eliminate l2pcdim (dimension parameters) move to allocatable ..
!
! Revision 1.4  2001/03/30 20:28:21  zvi
! General fix-up to get rid of COMMON BLOCK (ELLIPSE)
!
! Revision 1.3  2001/03/29 08:51:01  zvi
! Changing the (*) toi (:) everywhere
!
! Revision 1.2  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.1  2000/06/21 21:56:13  zvi
! First version D.P.
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
