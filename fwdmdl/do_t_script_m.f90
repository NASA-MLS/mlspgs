module DO_T_SCRIPT_M
  use MLSCommon, only: IP, R8, RP
  use D_STAT_TEMP_M, only: STAT_TEMP
  implicit NONE
  private
  PUBLIC :: DO_T_SCRIPT, TWO_D_T_SCRIPT
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
     "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
     "$RCSfile$"
!---------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------
! This routine builds the t_script array
!
  Subroutine DO_T_SCRIPT (Ngp1,Frq,sp_tmp,brkpt,no_ele,t_path,mid,t_script)
!
! 'sp_tmp' is space temperature.
!
    Integer(ip), intent(in) :: Ngp1, brkpt, no_ele
!
    Real(rp), intent(in) :: Frq, sp_tmp

    Real(rp), intent(in) :: T_PATH(:)

    Integer(ip), intent(out) :: mid
    Real(rp), intent(out) :: t_script(:)
!
    Integer(ip) :: i, j, m
    real(rp)    :: Tb(Size(t_script))

    tb(1:) = 0.0_rp
    t_script(1:) = 0.0_rp

! 'brkpt' is the index of the path break-point (when it change from
!         incoming ray to outgoing ray)
! 'no_ele' is the total number of entries in t_path(1...no_ele)

    m = 0
    j = 1 - Ngp1

    do
      j = j + Ngp1
      if (j > brkpt) EXIT
      m = m + 1
      Tb(m) = stat_temp(t_path(j),Frq)
    end do
!
    mid = m
    j = brkpt + 1 - Ngp1

    do
      j = j + Ngp1
      if (j > no_ele) EXIT
      m = m + 1
      Tb(m) = stat_temp(t_path(j),Frq)
    end do
!
    t_script(1) = 0.5_rp * (Tb(1) + Tb(2))
    do i = 2, mid-1
      t_script(i) = 0.5_rp * (Tb(i+1) - Tb(i-1))
    end do
!
    t_script(mid  ) = 0.5_rp * (Tb(mid+1) - Tb(mid-1))
    t_script(mid+1) = 0.5_rp * (Tb(mid+2) - Tb(mid  ))
!
! Note that Tb(mid) = Tb(mid+1) because this index is redundant
! at the tangent vertical (or center phi)
!
    do i = mid+2, m-1
      t_script(i) = 0.5_rp * (Tb(i+1) - Tb(i-1))
    end do
!
    t_script(m) = stat_temp(sp_tmp,Frq) - 0.5_rp * (Tb(m-1) + Tb(m))
!
    Return
  End Subroutine DO_T_SCRIPT
!
!---------------------------------------------------------------
!
  SUBROUTINE two_d_t_script(t_grid,t_space,nu,t_scr)
!
! INPUTS:
!
  REAL(rp), INTENT(in) :: t_grid(:) ! path temperatures
  REAL(r8), INTENT(in) :: nu ! Frequency
  REAL(rp), INTENT(in) :: t_space    ! space temperature
!
! OUTPUTS:
!
  REAL(rp), INTENT(out) :: t_scr(:) ! path differential temperatures
!
! Internal stuff
!
  INTEGER :: n_path
  Real(rp) :: a, b
  Real(rp) :: V1(SIZE(t_grid))
!
! get path information
!
  n_path = SIZE(t_grid)
  V1 = stat_temp(t_grid,nu)
  b = -stat_temp(t_grid(1),nu)
  a = 2.0_rp*stat_temp(t_space,nu)-stat_temp(t_grid(n_path),nu)
!
  t_scr = 0.5_rp * (EOSHIFT(V1,1,a) - EOSHIFT(V1,-1,b))
!
  END SUBROUTINE two_d_t_script
!
END module DO_T_SCRIPT_M
! $Log$
! Revision 1.7.2.1  2001/09/10 10:02:32  zvi
! Cleanup..comp_path_entities_m.f90
!
! Revision 1.1  2000/08/30 18:12:04  vsnyder
! Initial conversion to Fortran 90
