!
module GET_RAD_M
  use MLSCommon, only: I4, R4, R8
  use D_QSORT_M, only: QSORT
  use D_LINTRP_M, only: LINTRP
  use TRAPEZ_INT_M, only: TRAPEZ_INT
  Implicit NONE
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

Real(r8) FUNCTION get_rad(f_grid,f_grid_fltr,fltr_func,rad,nvr,nfp,ier)

!  ===============================================================
!  Declaration of variables for sub-program: get_rad
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: nvr, nfp
Integer(i4), INTENT(OUT) :: ier

Real(r8), INTENT(IN) :: rad(:), f_grid(:), f_grid_fltr(:), fltr_func(:)
!  ----------------------
!  PARAMETER Declaration:
!  ----------------------
Integer(i4), PARAMETER :: m = 370
!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: j, k, n

Real(r8) :: rxf(m), tmpary(m), tmpf(m), tmpflt(m), q, fr, pf
!
! Begin code
!
  ier = 2000
  IF(nfp > m) THEN
    PRINT *,'** Error in subroutine Get_Rad !'
    PRINT *,'   Number of fltr points too large:',nfp
    PRINT *,'   Maximum allowed:',m
    RETURN
  END IF

  n = 0
  ier = 0
  DO k = 1, nfp
    n = n + 1
    tmpf(n) = f_grid_fltr(k)
  END DO

  DO k = 1, nvr
    n = n + 1
    tmpf(n) = f_grid(k)
    rxf(k) = rad(k)
  END DO

  CALL Qsort(n,tmpf)

  j = n
  n = 0
  pf = -1.0D0
  DO k = 1, j
    fr = tmpf(k)
    IF(fr > pf) THEN
      pf = fr
      n = n + 1
      tmpf(n) = fr
    END IF
  END DO

  CALL lintrp(f_grid_fltr,tmpf,fltr_func,tmpflt,nfp,n)
  CALL lintrp(f_grid,tmpf,rxf,tmpary,nvr,n)

  DO k = 1, n
    rxf(k) = tmpary(k) * tmpflt(k)
  END DO

  CALL trapez_int(n,tmpf,rxf,q)
  get_rad = q

  RETURN
END FUNCTION get_rad
end module GET_RAD_M
! $Log$
! Revision 1.1 2000/06/09 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
