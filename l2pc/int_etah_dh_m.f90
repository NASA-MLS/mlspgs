!
module INT_ETAH_DH_M
  use MLSCommon, only: I4, R4, R8
  Implicit NONE
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------
! This subroutine compute the integral of{Eta(iq,h)*h*dh/Sqrt(h*h-Ht*Ht)}
! Where Eta(iq,h) is the Eta function for the iq'th coefficient at h

SUBROUTINE int_etah_dh(h_lo,h_hi,iq,no_tc,ht2,hqm1,hq,hqp1,eh)

!  ===============================================================
!  Declaration of variables for sub-program: int_etah_dh
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: iq, no_tc

Real(r8), INTENT(IN) :: ht2, h_lo, h_hi, hqm1, hq, hqp1

Real(r8), INTENT(OUT) :: eh
!  ----------------
!  Local variables:
!  ----------------
Real(r8) :: a, b, r, p, t, dh, sa, sb, deh
!
! Begin code
!
  deh = 0.0D0

  IF(iq > 1) THEN            ! do standard uprise integral

    a = MAX(h_lo,hqm1)
    b = MIN(h_hi,hq)

    IF (a < b) THEN          ! non zero contribution here
      dh = hq - hqm1
      sa = SQRT(a*a-ht2)
      sb = SQRT(b*b-ht2)
      r = sa - sb
      p = b*sb - a*sa
      t = LOG((a+sa)/(b+sb))
      deh = (hqm1 * r + 0.5D0 * (p - ht2 * t)) / dh
    END IF

  END IF

  IF(iq < no_tc) THEN        ! do standard down rise integral

    a = MAX(h_lo,hq)
    b = MIN(h_hi,hqp1)

    IF(a < b) THEN           ! non zero contribution here
      dh = hqp1 - hq
      sa = SQRT(a*a-ht2)
      sb = SQRT(b*b-ht2)
      r = sb - sa
      p = b*sb - a*sa
      t = LOG((a+sa)/(b+sb))
      deh = deh + (hqp1 * r + 0.5D0 * (ht2 * t - p)) / dh
    END IF

  END IF

! Now look for special cases

  r = 0.0D0
  IF(iq == 1) THEN           ! up rising on lowest coefficient

    a = h_lo
    b = MIN(h_hi,hq)
    if(a < b) r = SQRT(b*b-ht2) - SQRT(a*a-ht2)   ! Non-zero contribution

  ELSE IF(iq == no_tc) THEN      ! down rising on topmost coefficient

    b = h_hi
    a = MAX(h_lo,hq)
    IF(a < b) r = SQRT(b*b-ht2) - SQRT(a*a-ht2)  ! Non-zero contribution

  END IF

  eh = deh + r

  RETURN
END SUBROUTINE int_etah_dh
end module INT_ETAH_DH_M
! $Log$
! Revision 1.1 2000/06/09 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
