!
module CH_RUN_SETUP_M
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
!-----------------------------------------------------------------------

SUBROUTINE ch_run_setup(conv_hts,no_conv_hts,h_grid,t_grid,n_lvls, &
                        n_tan,t_tan)

!  ===============================================================
!  Declaration of variables for sub-program: ch_run_setup
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: n_lvls, no_conv_hts

Integer(i4), INTENT(OUT) :: n_tan(:)

Real(r8), INTENT(IN) :: conv_hts(:), t_grid(:), h_grid(:)

Real(r8), INTENT(OUT) :: t_tan(:)
!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: h_i, ptr, ptg_i

Real(r8) :: h, h1, h_tan, dh, ht

!  Find the tangent height for pointing ptg_i

  ptr = -1
  h1 = h_grid(1)
  DO ptg_i = 1, no_conv_hts

    h = conv_hts(ptg_i)
    CALL hunt(h,h_grid,n_lvls,ptr,h_i)
    IF(ABS(h-h_grid(h_i)) < ABS(h-h_grid(ptr))) ptr = h_i

    IF(ptr < n_lvls) THEN

      h_tan = h
      n_tan(ptg_i) = ptr
      ht = MAX(h_tan,h1)
      dh = (ht - h_grid(ptr)) / (h_grid(ptr+1) - h_grid(ptr))
      t_tan(ptg_i) = t_grid(ptr)+(t_grid(ptr+1)-t_grid(ptr)) * dh

    ELSE

      n_tan(ptg_i) = n_lvls
      t_tan(ptg_i) = t_grid(n_lvls)

    END IF

  END DO

  RETURN
END SUBROUTINE ch_run_setup
end module CH_RUN_SETUP_M
! $Log$
! Revision 1.1 2000/06/09 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
