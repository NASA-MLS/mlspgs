! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module DO_T_SCRIPT_M

  implicit NONE
  private
  public :: DO_T_SCRIPT, TWO_D_T_SCRIPT

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!----------------------------------------------------  DO_T_SCRIPT  -----

  subroutine DO_T_SCRIPT ( Ngp1, Frq, sp_tmp, brkpt, no_ele, t_path, mid, &
    &                      t_script)

! This routine builds the t_script array

    use D_STAT_TEMP_M, only: STAT_TEMP
    use MLSCommon, only: IP, RP

! 'sp_tmp' is space temperature.

    integer(ip), intent(in) :: Ngp1, brkpt, no_ele

    real(rp), intent(in) :: Frq, sp_tmp

    real(rp), intent(in) :: T_PATH(:)

    integer(ip), intent(out) :: mid
    real(rp), intent(out) :: t_script(:)

    integer(ip) :: i, j, m
    real(rp)    :: Tb(Size(t_script))

    tb(1:) = 0.0_rp
    t_script(1:) = 0.0_rp

! 'brkpt' is the index of the path break-point (when it changes from
!         incoming ray to outgoing ray)
! 'no_ele' is the total number of entries in t_path(1...no_ele)

    m = 0
    j = 1 - Ngp1

    do
      j = j + Ngp1
      if ( j > brkpt ) exit
      m = m + 1
      Tb(m) = stat_temp(t_path(j),Frq)
    end do

    mid = m
    j = brkpt + 1 - Ngp1

    do
      j = j + Ngp1
      if ( j > no_ele ) exit
      m = m + 1
      Tb(m) = stat_temp(t_path(j),Frq)
    end do

    t_script(1) = 0.5_rp * (Tb(1) + Tb(2))
    do i = 2, mid-1
      t_script(i) = 0.5_rp * (Tb(i+1) - Tb(i-1))
    end do

    t_script(mid  ) = 0.5_rp * (Tb(mid+1) - Tb(mid-1))
    t_script(mid+1) = 0.5_rp * (Tb(mid+2) - Tb(mid  ))

! Note that Tb(mid) = Tb(mid+1) because this index is redundant
! at the tangent vertical (or center phi)

    do i = mid+2, m-1
      t_script(i) = 0.5_rp * (Tb(i+1) - Tb(i-1))
    end do

    t_script(m) = stat_temp(sp_tmp,Frq) - 0.5_rp * (Tb(m-1) + Tb(m))

    return
  end subroutine DO_T_SCRIPT

!------------------------------------------------  Two_D_T_Script  -----

  subroutine Two_D_T_Script ( t_grid, t_space, nu, t_scr )

    use D_STAT_TEMP_M, only: STAT_TEMP
    use MLSCommon, only: R8, RP

! Inputs:

    real(rp), intent(in) :: t_grid(:) ! path temperatures
    real(r8), intent(in) :: nu        ! Frequency
    real(rp), intent(in) :: t_space   ! farside boundary temperature
    !                                   usually cosmic space (2.7K).

! Outputs:

    real(rp), intent(out) :: t_scr(:) ! path differential temperatures

! Internal stuff

    integer :: n_path
    real(rp) :: a, b
    real(rp) :: V1(size(t_grid))

! get path information

    n_path = size(t_grid)
    V1 = stat_temp(t_grid,nu)
    b = -stat_temp(t_grid(1),nu)
    a = 2.0_rp * stat_temp(t_space,nu)-stat_temp(t_grid(n_path),nu)

    t_scr = 0.5_rp * (eoshift(V1,1,a) - eoshift(V1,-1,b))

  end subroutine Two_D_T_Script

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module DO_T_SCRIPT_M

! $Log$
! Revision 2.2  2002/10/08 17:08:02  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.1  2002/09/27 00:04:15  vsnyder
! Insert copyright notice, cosmetic changes
!
! Revision 2.0  2001/09/17 20:26:26  livesey
! New forward model
!
! Revision 1.7.2.1  2001/09/10 10:02:32  zvi
! Cleanup..comp_path_entities_m.f90
!
! Revision 1.1  2000/08/30 18:12:04  vsnyder
! Initial conversion to Fortran 90
