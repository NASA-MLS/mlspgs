!
module GET_ANGLES_M
  use MLSCommon, only: I4, R4, R8
  use D_HUNT_M, only: HUNT
  Implicit NONE
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains

!----------------------------------------------------------------------
! Set up array of pointing angles

SUBROUTINE get_angles(h_grid,conv_hts,n_grid,a_grid,ptg_angle, &
                      roc,h_obs,n_lvls,n_obs,no_conv_hts,ier)

!  ===============================================================
!  Declaration of variables for sub-program: get_angles
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: n_lvls, n_obs, no_conv_hts
Integer(i4), INTENT(OUT) :: ier

Real(r8), INTENT(IN) :: roc, h_obs, h_grid(:), n_grid(:), conv_hts(:)

Real(r8), INTENT(OUT) :: ptg_angle(:), a_grid(:)
!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: i, j, k

Real(r8) :: q, ht, ang

! Get 'N_lvls' pointing angles

  ier = 0
  DO i = 1, n_lvls
    ht = h_grid(i)
    q = 1.0 + n_grid(i)
    ang = q * (ht + roc) / h_obs
    IF(ABS(ang) > 1.0) THEN
      ier = 1
      PRINT *,'*** Error in get_angles subroutine !'
      PRINT *,'    arg > 1.0 in ArcSin(arg) ..'
      RETURN
    END IF
    a_grid(i) = DASIN(ang)
  END DO

! Get 'No_conv_hts' convolution angles

  j = -1
  DO i = 1, no_conv_hts
    ht = conv_hts(i)
    CALL hunt(ht,h_grid,n_lvls+1,j,k)
    q = 1.0 + n_grid(j) - n_grid(n_obs)
    ang = q * (ht + roc) / h_obs
    IF(ABS(ang) > 1.0) THEN
      ier = 1
      PRINT *,'*** Error in get_angles subroutine !'
      PRINT *,'    arg > 1.0 in ArcSin(arg) ..'
      RETURN
    END IF
    ptg_angle(i) = DASIN(ang)
  END DO

  RETURN
END SUBROUTINE get_angles
end module GET_ANGLES_M
! $Log$
! Revision 1.1 2000/06/09 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
