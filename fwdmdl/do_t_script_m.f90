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

  subroutine DO_T_SCRIPT ( t_path, sp_tmp, Frq, t_script )

! This routine builds the t_script array.  In some notes it's called Delta B.

    use D_STAT_TEMP_M, only: STAT_TEMP
    use MLSCommon, only: IP, RP

    real(rp), intent(in) :: T_PATH(:)    ! path temperatures            
    real(rp), intent(in) :: Frq          ! Frequency                    
    real(rp), intent(in) :: sp_tmp       ! farside boundary temperature 
    !                                      usually cosmic space (2.7K). 

    real(rp), intent(out) :: t_script(:) ! path differential temperatures

    integer(ip) :: I, N_path
    real(rp)    :: Tb(Size(t_script))

    n_path = size(t_path)

    tb = stat_temp(t_path,Frq)

    t_script(1) = 0.5_rp * (Tb(1) + Tb(2))
    do i = 2, n_path-1
      t_script(i) = 0.5_rp * (Tb(i+1) - Tb(i-1))
    end do
    t_script(n_path) = stat_temp(sp_tmp,Frq) - 0.5_rp * (Tb(n_path-1) + Tb(n_path))

    return
  end subroutine DO_T_SCRIPT

!------------------------------------------------  Two_D_T_Script  -----

  subroutine Two_D_T_Script ( t_grid, t_scat, w0, t_space, nu, t_scr )

! This routine builds the t_script array.  In some notes it's called Delta B.

    use D_STAT_TEMP_M, only: STAT_TEMP
    use MLSCommon, only: R8, RP

! Inputs:

    real(rp), intent(in) :: t_grid(:) ! path temperatures
    real(rp), intent(in) :: t_scat(:) ! path scattering
    real(r8), intent(in) :: w0(:)     ! single scattering albedo
    real(r8), intent(in) :: nu        ! Frequency
    real(rp), intent(in) :: t_space   ! farside boundary temperature
    !                                   usually cosmic space (2.7K).

    real(rp) :: t_tmp(size(t_grid)) ! path temperatures modified by scattering

! Outputs:

    real(rp), intent(out) :: t_scr(:) ! path differential temperatures

! Internal stuff

    integer :: n_path, i
    real(rp) :: a, b
    real(rp) :: V1(size(t_grid))

! get path information

    n_path = size(t_grid)

    do i=1, n_path-1
       t_tmp(i) = (1.0_rp-w0(i))*t_grid(i) + w0(i)*t_scat(i)
        print*, t_tmp(i), t_grid(i), grid(i), w0(i), t_scat(i)
!       t_grid(i) = t_tmp(i)  !will add in later soon
    enddo
    stop

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
! Revision 2.4  2003/06/06 22:17:26  vsnyder
! Simplify do_t_script (which isn't used anywhere)
!
! Revision 2.3  2003/06/06 20:38:25  vsnyder
! Add some comments about dummy arguments
!
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
