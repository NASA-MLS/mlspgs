module DO_T_SCRIPT_M
  use MLSCommon, only: I4, R8
  use PATH_ENTITIES_M, only: PATH_VECTOR
  use D_STAT_TEMP_M, only: STAT_TEMP
  implicit NONE
  private
  public :: DO_T_SCRIPT
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
    Integer(i4), intent(in) :: Ngp1, brkpt, no_ele
!
    Real(r8), intent(in) :: Frq, sp_tmp

    Type(path_vector), intent(in) :: T_PATH

    Integer(i4), intent(out) :: mid
    Real(r8), intent(out) :: t_script(:)
!
    Integer(i4) :: i, j, m
    Real(r8) :: Tb(size(t_script))
!
!  Begin code:
!
    tb(:) = 0.0
    t_script(:) = 0.0

! 'brkpt' is the index of the path break-point (when it change from
!         incoming ray to outgoing ray)
! 'no_ele' is the total number of entries in t_path%values(1...no_ele)

    m = 0
    j = 1 - Ngp1

    do
      j = j + Ngp1
      if (j > brkpt) EXIT
      m = m + 1
      Tb(m) = stat_temp(t_path%values(j),Frq)
    end do
!
    mid = m
    j = brkpt + 1 - Ngp1

    do
      j = j + Ngp1
      if (j > no_ele) EXIT
      m = m + 1
      Tb(m) = stat_temp(t_path%values(j),Frq)
    end do
!
    t_script(1) = 0.5 * (Tb(1) + Tb(2))
    do i = 2, mid-1
      t_script(i) = 0.5 * (Tb(i+1) - Tb(i-1))
    end do
!
    t_script(mid  ) = 0.5 * (Tb(mid+1) - Tb(mid-1))
    t_script(mid+1) = 0.5 * (Tb(mid+2) - Tb(mid  ))
!
! Note that Tb(mid) = Tb(mid+1) because this index is redundant
! at the tangent vertical (or center phi)
!
    do i = mid+2, m-1
      t_script(i) = 0.5 * (Tb(i+1) - Tb(i-1))
    end do
!
    t_script(m) = stat_temp(sp_tmp,Frq) - 0.5 * (Tb(m-1) + Tb(m))
!
    Return
  End Subroutine DO_T_SCRIPT
end module DO_T_SCRIPT_M
! $Log$
! Revision 1.5  2001/03/29 08:51:01  zvi
! Changing the (*) toi (:) everywhere
!
! Revision 1.4  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.1  2000/08/30 18:12:04  vsnyder
! Initial conversion to Fortran 90
!
