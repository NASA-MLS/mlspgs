module GET_FILTERS_M
  use MLSCommon, only: I4, R4, R8
  use FILTER_SW_M, only: FILTER
  use DSIMPSON_MODULE, only: SIMPS
  Implicit NONE
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
  "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
  "$RCSfile$"
!---------------------------------------------------------------------------
contains
!----------------------------------------------------------------------

SUBROUTINE get_filters(no_pfa_ch,no_filt_pts,pfa_ch,f_grid_filter, &
               &       freqs,filter_func,InDir,ld,primag,ier)

!  ===============================================================
!  Declaration of variables for sub-program: get_filters
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: pfa_ch(*)
Integer(i4), INTENT(IN) :: no_pfa_ch, no_filt_pts

Integer(i4), INTENT(OUT) :: ier, ld

Real(r8), INTENT(IN) :: freqs(*)

Real(r8), INTENT(OUT) :: filter_func(:,:), f_grid_filter(:,:)

Character (LEN=*), INTENT(IN) :: InDir, primag

!  ----------------
!  Local variables:
!  ----------------

Integer(i4) :: j, ch_i, mch

Real(r8) :: xlhs, xrhs, df, q, area, frq

! Begin code:

  ier = 0

! Find the species index in the l2pc mixing ratio database:

  DO ch_i = 1, no_pfa_ch

    mch = pfa_ch(ch_i)
    frq = freqs(mch)
    IF(frq < 1.0D0) THEN
      ier = 1
      WRITE(6,900) mch
      Return
    END IF

! Set up filter's response function

    q = 0.0
    df = Filter(q,mch,xlhs,xrhs,area,ier,InDir,primag,ld)
    IF(ier /= 0) GO TO 99

    df = (xrhs-xlhs)/(no_filt_pts-1)
    DO j = 1, no_filt_pts
      q = xlhs + (j - 1) * df
      f_grid_filter(j,ch_i) = frq + q
      filter_func(j,ch_i) = Filter(q)
    END DO

!  Normalize the filter's response array:

    CALL Simps(filter_func(1:,ch_i),df,no_filt_pts,q)
    DO j = 1, no_filt_pts
      filter_func(j,ch_i) = filter_func(j,ch_i) / q
    END DO

  END DO                         ! On ch_i

 99 IF(ier /= 0) THEN
      PRINT *,' ** Error in get_filters subroutine **'
      PRINT *,'    After calling subroutine: Get_Filter_Param'
      CALL errmsg(' ',ier)
    END IF

  900  FORMAT(' ** Error in get_filters subroutine **',/, &
      '    Inconsistant User Input.',/, &
      '    PFA Channel:',i3,' not among the non-PFA channels !')

  Return

END SUBROUTINE get_filters
end module GET_FILTERS_M
! $Log$
! Revision 1.1 2000/06/09 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
