! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module DO_T_SCRIPT_M

  implicit NONE
  private
  public :: DO_T_SCRIPT, TWO_D_T_SCRIPT, TWO_D_T_SCRIPT_CLOUD

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains

!----------------------------------------------------  DO_T_SCRIPT  -----
  subroutine DO_T_SCRIPT ( t_path, sp_tmp, Frq, t_script )

! This routine builds the t_script array.  In some notes it's called Delta B.

    use D_STAT_TEMP_M, only: STAT_TEMP
    use MLSKinds, only: RP, R8

    real(rp), intent(in) :: T_PATH(:)    ! path temperatures            
    real(r8), intent(in) :: Frq          ! Frequency                    
    real(rp), intent(in) :: sp_tmp       ! farside boundary temperature 
    !                                      usually cosmic space (2.7K). 

    real(rp), intent(out) :: t_script(:) ! path differential temperatures

    integer  :: I, N_path
    real(rp) :: Tb(Size(t_script))

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
  subroutine Two_D_T_Script ( t_grid, t_space, nu, t_scr, B )


!{ This routine builds the {\tt t_script} array.
!  In some notes and reports it's called $\Delta B_i =
!  \frac12 ( B_{i+1} - B_{i-1} )$
!  for $i = 2 \dots n-1$, $\Delta B_1 = \frac12 ( B_2 + B_1 )$, and
!  $\Delta B_{n} = B_{\text{space}} - \frac12 ( B_{n} + B_{n-1} )$,
!  where $n$ = {\tt n_path}.

    use D_STAT_TEMP_M, only: STAT_TEMP
    use MLSKinds, only: R8, RP

! Inputs:

    real(rp), intent(in) :: t_grid(:) ! path temperatures
    real(r8), intent(in) :: nu        ! Frequency
    real(rp), intent(in) :: t_space   ! farside boundary temperature
    !                                   usually cosmic space (2.7K).

!    real(rp) :: t_tmp(size(t_grid))  ! path temperatures modified by scattering

! Outputs:

    real(rp), intent(out) :: t_scr(:) ! path differential temperatures
    real(rp), intent(out) :: B(:)     ! Planck radiation function on path

! Internal stuff

    integer :: n_path

! get path information

    n_path = size(t_grid)
  
    B = stat_temp(t_grid,nu)
    t_scr = 0.5_rp * (eoshift(B,1,2.0_rp * stat_temp(t_space,nu) - B(n_path)) - &
      &               eoshift(B,-1,-B(1)))

  end subroutine Two_D_T_Script

!------------------------------------------  Two_D_T_Script_Cloud  -----
  subroutine Two_D_T_Script_Cloud ( t_grid, t_scat, w0, t_space, nu, t_scr, B )

!{ This routine builds the {\tt t_script} array.  In some notes it's called
!  Delta B. This one combines B with TScat.  In the cloud scattering case,
!  $\Delta B_i = \frac12 \left[ (1-\omega_{0_{i+1}}) B_{i+1} +
!                             \omega_{0_{i+1}} T_{\text{scat}_{i+1}} -
!                             (1-\omega_{0_{i-1}}) B_{i-1} -
!                             \omega_{0_{i-1}} T_{\text{scat}_{i-1}} \right]$,
!  for $i = 2 \dots\ n-1$,
!  $\Delta B_1 = \frac12 \left[ (1-\omega_{0_1}) B_1 +
!                             \omega_{0_1} T_{\text{scat}_1} +
!                             (1-\omega_{0_2}) B_2 +
!                             \omega_{0_2} T_{\text{scat}_2} \right]$,
!  and
!  $\Delta B_n = B(\text{\tt t_space}) -
!                \frac12 \left[ (1-\omega_{0_n}) B_n +
!                             \omega_{0_n} T_{\text{scat}_n} +
!                             (1-\omega_{0_{n-1}}) B_{n-1} +
!                             \omega_{0_{n-1}} T_{\text{scat}_{n-1}} \right]$,


    use D_STAT_TEMP_M, only: STAT_TEMP
    use MLSKinds, only: R8, RP

! Inputs:

    real(rp), intent(in) :: t_grid(:) ! path temperatures
    real(rp), intent(in) :: t_scat(:) ! path scattering
    real(rp), intent(in) :: w0(:)     ! single scattering albedo
    real(r8), intent(in) :: nu        ! Frequency
    real(rp), intent(in) :: t_space   ! farside boundary temperature
    !                                   usually cosmic space (2.7K).

! Outputs:

    real(rp), intent(out) :: t_scr(:) ! path differential temperatures
    real(rp), intent(out) :: B(:)     ! Planck radiation function on path

! Internal stuff

    integer :: n_path
    real(rp) :: V1(size(t_grid))

! get path information

    n_path = size(t_grid)

    B = stat_temp(t_grid,nu)                 ! clear-sky Planck  function on path
    V1 = (1.0_rp-w0) * B + w0 * t_scat       !(1-w0)*B + w0*Tscat

    t_scr = 0.5_rp * (eoshift(V1,1,2.0_rp * stat_temp(t_space,nu) - B(n_path)) - &
      &               eoshift(V1,-1,-B(1)))

  end subroutine Two_D_T_Script_Cloud

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module DO_T_SCRIPT_M

! $Log$
! Revision 2.14  2010/12/15 21:36:37  vsnyder
! TeXnicalities
!
! Revision 2.13  2010/08/19 02:04:50  vsnyder
! Cannonball polishing
!
! Revision 2.12  2009/06/23 18:26:11  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.11  2005/11/01 23:02:21  vsnyder
! PFA Derivatives
!
! Revision 2.10  2005/06/22 18:08:18  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.9  2004/03/31 20:35:17  jonathan
! correction on the cloud effect term
!
! Revision 2.8  2004/03/30 00:45:09  vsnyder
! Remove unused variable declaration
!
! Revision 2.7  2004/03/27 01:36:56  jonathan
! now cloud effects are included
!
! Revision 2.6  2004/03/20 01:25:06  jonathan
! minor changes
!
! Revision 2.5  2004/03/20 01:15:10  jonathan
! add in scattering correction term in two t_script
!
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
