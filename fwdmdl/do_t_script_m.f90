module DO_T_SCRIPT_M
  use L2PCdim, only: N2LVL
  use MLSCommon, only: I4, R4, R8
  use S_STAT_TEMP_M, only: STAT_TEMP
  implicit NONE
  private
  public :: DO_T_SCRIPT

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
! This routine builds the t_script array
!
  Subroutine DO_T_SCRIPT ( N_lvls, Ng, Frq, sp_tmp, t_path, IndxR,   &
 &                         IndxL, path_brkpt, t_script )
!
! sp_tmp is space temperature. The begining of the ray has index 1, and
! the end of the ray has index: 2*n. The uppermost index of the ray on the
! right is: IndxR, and the lowermost index ot the ray on the left is: IndxL
!
    Integer(i4), intent(in) :: N_lvls, Ng
    Real(r8), intent(in) :: Frq
    Real(r4), intent(in) :: sp_tmp, t_path(*)
    Integer(i4), intent(in) :: IndxR, IndxL, path_brkpt(*)
    Real(r4), intent(out) :: t_script(*)

    Integer(i4) :: Ngp1,i,j,k
    real(r4) :: Tb(N2lvl)

    tb = 0.0
    t_script(1:n2lvl) = 0.0
!
    j = -Ng
    Ngp1 = Ng + 1
    do i = 1, IndxR
      j = j + Ngp1
      if ( j <= path_brkpt(1) ) Tb(i) = stat_temp(t_path(j),Frq)
    end do
!
    k = 2 * N_lvls
    j = path_brkpt(2) - Ngp1
    do i = IndxL, k
      j = j + Ngp1
      Tb(i) = stat_temp(t_path(j),Frq)
    end do
!
    t_script(1) = 0.5 * (Tb(1) + Tb(2))
    do i = 2, IndxR-1
      t_script(i) = 0.5 * (Tb(i+1) - Tb(i-1))
    end do
!
    t_script(IndxR) = 0.5 * (Tb(IndxL) - Tb(IndxR-1))
    t_script(IndxL) = 0.5 * (Tb(IndxL+1) - Tb(IndxR))
!
! Note that Tb(IndxR) = Tb(IndxL) because this index is redundant
! at the tangent vertical (or center phi)
!
    do i = IndxL+1, k-1
      t_script(i) = 0.5 * (Tb(i+1) - Tb(i-1))
    end do
!
    t_script(k) = stat_temp(sp_tmp,Frq) - 0.5 * (Tb(k-1) + Tb(k))
!
    Return
  End Subroutine DO_T_SCRIPT
end module DO_T_SCRIPT_M

! $Log$
