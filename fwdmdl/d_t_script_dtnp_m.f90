! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module D_T_SCRIPT_DTNP_M

  implicit NONE
  private
  public :: DT_SCRIPT_DT, DT_SCRIPT_DT_W0, DT_SCRIPT_DF_W0, DT_SCRIPT

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

  ! -----------------------------------------------  DT_SCRIPT_DT  -----
  subroutine DT_SCRIPT_DT ( T_PATH, B, ETA_ZXP, NZ_ZXP, NNZ_ZXP, &
    &                       NU, DT_SCR_DT )
    ! Build the derivative of the B array w.r.t to T_np
    ! (The 'n' zeta coeff. and the 'p' phi coefficient)

!{ \parskip 5pt
!  Given $T$, $\nu$ and the Planck function
!  $B = \frac{\frac{h \nu}k}{exp\left(\frac{h \nu}{k T}\right)-1}$,
!  compute $\frac{\text{d} B}{\text{d} T}$.
!
!  B satisfies the differential equation
!  $\frac{\text{d} B}{\text{d} T} =
!   \frac{B}{T^2} \left( \frac{h\nu}k + B \right)$.
!
!  The old way was to compute $\frac{\text{d} B}{\text{d} T} =
!   \exp \left ( \frac{h \nu}{k T} \right )
!   \left ( \frac{h \nu}
!     {k T \left [ \exp \left ( \frac{h \nu}{k T} \right ) - 1 \right ]}
!     \right ) ^2$.  Dividing numerator and denominator by
!     $\exp \left ( \frac{h \nu}{k T} \right )$ gives
!     $\frac{\left ( \frac{h \nu}{k T} \right )^2}
!           { 2 \left [ \cosh \left ( \frac{h \nu}{k T} \right ) - 1 \right ]}$,
!   which is substantially more expensive than using the differential equation.
!   $2 ( \cosh x \, - \, 1) / x^2$ has substantial cancellation for small $x$,
!   so a specially-developed procedure should be used to evaluate it.
!
!  From $\frac{\text{d} B}{\text{d} T}$ compute
!  $\frac{\text{d} B}{\text{d} T_{np}} =
!   \frac{\text{d} B}{\text{d} T} \eta_n \eta_p$ where
!  $\eta_n$ = Eta(n,zeta) and $\eta_p$ = Eta(p,phi).
!  Then difference $\frac{\text{d} B}{\text{d} T_{np}}$ along
!  the path to get {\tt dT\_scr\_dT} = $\frac{\text{d} \Delta B}{\text{d} T}$,
!  where $\Delta B = \frac12 \left( B_{i+1} - B_{i-1} \right)$.
!  $\Delta B$  is called ``T script'' (not Tau) in many notes and reports.

    use MLSCommon, only: R8, RP
    use PHYSICS, only: H_O_K => h_over_k

! inputs

    real(rp), intent(in) :: t_path(:)     ! path temperatures K
    real(rp), intent(in) :: B(:)          ! Planck function
    real(rp), intent(in) :: eta_zxp(:,:)  ! path eta functions
    integer, intent(in) :: nz_zxp(:,:)    ! Where eta_zxp is not zero
    integer, intent(in) :: nnz_zxp(:)     ! Numbers of rows in nz_zxp
    real(r8), intent(in) :: nu            ! calculation frequency (MHz)

! output

    real(rp), intent(out) :: dT_scr_dT(:,:) ! path dt_script temperature
!                                           derivative

! internals

    real(rp) :: dBdT(size(B))

    dBdT = B * ( h_o_k * nu + B ) / t_path**2
    call dT_script ( dBdT, eta_zxp, nz_zxp, nnz_zxp, dT_scr_dT )

  end subroutine DT_SCRIPT_DT

  ! --------------------------------------------  DT_SCRIPT_DT_W0  -----
  subroutine DT_SCRIPT_DT_W0 ( T_PATH, B, W0, dW0_dT, &
    &                          ETA_ZXP, NZ_ZXP, NNZ_ZXP, &
    &                          NU, DT_SCR_DT )
    ! Build the derivative of the B array w.r.t to T_np
    ! (The 'n' zeta coeff. and the 'p' phi coefficient)

!{ \parskip 5pt
!  Given $T$, $\nu$, $\omega_0$, $\frac{\partial \omega_0}{\partial T}$,
!  and the Planck function
!  $B = \frac{\frac{h \nu}k}{exp\left(\frac{h \nu}{k T}\right)-1}$,
!  compute \\
!  $\frac{\partial (1-\omega_0) B}{\partial T} =
!  -B \frac{\partial \omega_0}{\partial T} +
!     (1-\omega_0) \frac{\partial B}{\partial T}$.
!
!  B satisfies the differential equation
!  $\frac{\text{d} B}{\text{d} T} =
!   \frac{B}{T^2} \left( \frac{h\nu}k + B \right)$, so\\
!  $\frac{\partial (1-\omega_0) B}{\partial T} =
!   B \left( \frac{1-\omega_0}{T^2} \left(\frac{h \nu}k + B \right) -
!   \frac{\partial \omega_0}{\partial T} \right)$.
!
!  From $\frac{\partial (1-\omega_0) B}{\partial T}$ compute
!  $\frac{\partial (1-\omega_0) B}{\partial T_{np}} =
!   \frac{\partial (1-\omega_0) B}{\partial T} \eta_n \eta_p$ where
!  $\eta_n$ = Eta(n,zeta) and $\eta_p$ = Eta(p,phi).
!  Then difference $\frac{\partial (1-\omega_0) B}{\partial T_{np}}$ along
!  the path to get {\tt dT\_scr\_dT} =
!  $\frac{\partial \Delta (1-\omega_0) B}{\partial T}$, where
!  $\Delta (1-\omega_0) B = \frac12 \left( (1-\omega_{0_{i+1}}) B_{i+1} -
!                                          (1-\omega_{0_{i-1}}) B_{i-1} \right)$.
!  $\Delta B$  is called ``T script'' (not Tau) in many notes and reports.

    use MLSCommon, only: R8, RP
    use PHYSICS, only: H_O_K => h_over_k

! inputs

    real(rp), intent(in) :: T_path(:)     ! path temperatures K
    real(rp), intent(in) :: B(:)          ! Planck function
    real(rp), intent(in) :: W0(:)         ! \omega_0
    real(rp), intent(in) :: dW0_dT(:)     ! \frac{\partial \omega_0}{\partial T}
    real(rp), intent(in) :: eta_zxp(:,:)  ! path eta functions
    integer, intent(in) :: nz_zxp(:,:)    ! Where eta_zxp is not zero
    integer, intent(in) :: nnz_zxp(:)     ! Numbers of rows in nz_zxp
    real(r8), intent(in) :: nu            ! calculation frequency (MHz)

! output

    real(rp), intent(out) :: dT_scr_dT(:,:) ! path dt_script temperature
!                                           derivative

! internals

    real(rp) :: dBdT(size(B))

    dBdT = B * ( ( 1.0_rp - w0 ) * &
         &       ( h_o_k * nu + B ) / t_path**2 - dW0_dT )
    call dT_script ( dBdT, eta_zxp, nz_zxp, nnz_zxp, dT_scr_dT )

  end subroutine DT_SCRIPT_DT_W0

  ! --------------------------------------------  DT_SCRIPT_DF_W0  -----
  subroutine DT_SCRIPT_DF_W0 ( B, dW0_dX, ETA_ZXP, NZ_ZXP, NNZ_ZXP, &
    &                          DT_SCR_DF )
    ! Build the derivative of the B array w.r.t to mixing ratio F_np
    ! (The 'n' zeta coeff. and the 'p' phi coefficient)

!{ \parskip 5pt
!  Given $\frac{\partial \omega_0}{\partial f}$ and the Planck function
!  $B = \frac{\frac{h \nu}k}{exp\left(\frac{h \nu}{k T}\right)-1}$,
!  compute \\
!  $\frac{\partial (1-\omega_0) B}{\partial f} =
!  -B \frac{\partial \omega_0}{\partial f}$
!  ($B$ does not depend upon $f$).
!
!  From $\frac{\partial (1-\omega_0) B}{\partial f}$ compute
!  $\frac{\partial (1-\omega_0) B}{\partial f_{np}} =
!   \frac{\partial (1-\omega_0) B}{\partial f} \eta_n \eta_p$ where
!  $\eta_n$ = Eta(n,zeta) and $\eta_p$ = Eta(p,phi).
!  Then difference $\frac{\partial (1-\omega_0) B}{\partial f_{np}}$ along
!  the path to get {\tt dT\_scr\_df} =
!  $\frac{\partial \Delta (1-\omega_0) B}{\partial f}$, where
!  $\Delta (1-\omega_0) B = \frac12 \left( (1-\omega_{0_{i+1}}) B_{i+1} -
!                                          (1-\omega_{0_{i-1}}) B_{i-1} \right)$.
!  $\Delta B$  is called ``T script'' (not Tau) in many notes and reports.

    use MLSCommon, only: RP

! inputs

    real(rp), intent(in) :: B(:)          ! Planck function
    real(rp), intent(in) :: dW0_dX(:)     ! \frac{\partial \omega_0}{\partial X}
    real(rp), intent(in) :: eta_zxp(:,:)  ! path eta functions
    integer, intent(in) :: nz_zxp(:,:)    ! Where eta_zxp is not zero
    integer, intent(in) :: nnz_zxp(:)     ! Numbers of rows in nz_zxp

! output

    real(rp), intent(out) :: DT_SCR_DF(:,:) ! path dt_script temperature
!                                           derivative

    call dT_script ( -B * dW0_dX, eta_zxp, nz_zxp, nnz_zxp, dT_scr_df )

  end subroutine DT_SCRIPT_DF_W0

  ! --------------------------------------------------  DT_SCRIPT  -----
  subroutine DT_SCRIPT ( dBdx, ETA_ZXP, NZ_ZXP, NNZ_ZXP, DT_SCR, W0, dTScat )

!{ From {\tt dBdx} = $\frac{\partial B}{\partial x}$ compute
!  $\frac{\partial B}{\partial x_{np}} =
!   \frac{\partial B}{\partial x} \eta_n \eta_p$ where
!  $\eta_n$ = Eta(n,zeta), $\eta_p$ = Eta(p,phi), and $\eta_n \eta_p$ =
!  {\tt ETA\_ZXP}.
!  Then difference $\frac{\partial B}{\partial x_{np}}$ along
!  the path to get {\tt dT\_scr} =
!  $\frac{\partial \Delta B}{\partial x_{np}}$, where
!  $\Delta B = \frac12 \left( B_{i+1} - B_{i-1} \right)$.
!  $\Delta B$ is called ``T script'' (not Tau) in many notes and reports.

    use MLSCommon, only: RP

! inputs

    real(rp), intent(in) :: dBdx(:)       ! Derivative of Planck function
    real(rp), intent(in) :: eta_zxp(:,:)  ! path eta functions
    integer, intent(in) :: nz_zxp(:,:)    ! Where eta_zxp is not zero
    integer, intent(in) :: nnz_zxp(:)     ! Numbers of rows in nz_zxp
    real(rp), intent(in), optional :: W0(:)       ! On the path
    real(rp), intent(in), optional :: dTScat(:,:) ! Path X SV.  On the path,
                                          ! w.r.t. molecules in the grid.

! output

    real(rp), intent(out) :: dt_scr(:,:)  ! path dt_script dx derivative

! internals

    integer :: j, n_path, n_sv, p_i, sv_i

    real(rp) :: dT_x_eta(1:size(dBdx),1:size(eta_zxp,2))

    n_path = size(dBdx)
    n_sv = size(eta_zxp,dim=2)

    ! Make sure the elements of dT_x_eta that we use have defined values
    dT_x_eta(1:2,:) = 0.0
    dT_x_eta(n_path-1:n_path,:) = 0.0
    ! Fill other nonzero elements of dT_x_eta
    do sv_i = 1, n_sv
      dT_x_eta(nz_zxp(:nnz_zxp(sv_i),sv_i),sv_i) = &
        & 0.5 * dbdx(nz_zxp(:nnz_zxp(sv_i),sv_i)) * &
        & eta_zxp(nz_zxp(:nnz_zxp(sv_i),sv_i),sv_i)
    end do

    if ( present(w0) .and. present(dTScat) ) then
      do sv_i = 1, n_sv
        dT_x_eta(nz_zxp(:nnz_zxp(sv_i),sv_i),sv_i) = &
          & dT_x_eta(nz_zxp(:nnz_zxp(sv_i),sv_i),sv_i) + &
          & 0.5 * w0 * dTScat(nz_zxp(:nnz_zxp(sv_i),sv_i),sv_i)
      end do
    end if

    do sv_i = 1 , n_sv
      dt_scr(1,sv_i) = dT_x_eta(1,sv_i) + dT_x_eta(2,sv_i)
      dt_scr(n_path,sv_i) = -dT_x_eta(n_path-1,sv_i) - dT_x_eta(n_path,sv_i)
! This is what would be happening if we didn't pay attention to the nonzeros:
!     dt_scr(2:n_path-1,sv_i) = dT_x_eta(3:n_path,sv_i) - &
!         &                        dT_x_eta(1:n_path-2,sv_i)
! But we pay attention to the nonzeros to improve performance.
! If we kept track of the nonzeros in dt_scr, we wouldn't neet the next line.
      dt_scr(2:n_path-1,sv_i) = 0.0
      ! Now fill the nonzeros in column sv_i
      do j = 1, nnz_zxp(sv_i)
        p_i = nz_zxp(j,sv_i)
        if ( p_i > 2 ) dt_scr(p_i-1,sv_i) = dT_x_eta(p_i,sv_i)
      end do
      do j = 1, nnz_zxp(sv_i)
        p_i = nz_zxp(j,sv_i)
        if ( p_i >= n_path-1 ) exit
        dt_scr(p_i+1,sv_i) = dt_scr(p_i+1,sv_i) - dT_x_eta(p_i,sv_i)
      end do
    end do

  end subroutine DT_SCRIPT

  ! ----------------------------------------------  NOT_USED_HERE  -----
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module D_T_SCRIPT_DTNP_M
! $Log$
! Revision 2.11  2010/09/25 01:12:37  vsnyder
! Add dT_script_dT and dT_script_df, some cannonball polishing
!
! Revision 2.10  2010/06/23 02:38:19  vsnyder
! Improve TeXnicalities
!
! Revision 2.9  2010/06/23 02:28:39  vsnyder
! Correct nomenclature: T_Script should be B
!
! Revision 2.8  2009/06/23 18:26:11  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.7  2008/08/27 19:56:51  vsnyder
! Add PRINT to not_used_here
!
! Revision 2.6  2007/06/26 00:38:55  vsnyder
! Use column-sparse eta
!
! Revision 2.5  2005/11/01 23:02:21  vsnyder
! PFA Derivatives
!
! Revision 2.4  2005/06/22 18:08:18  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.3  2002/10/11 21:26:56  vsnyder
! use COSH1 to avoid errors for small nu
!
! Revision 2.2  2002/10/10 19:49:04  vsnyder
! Move USE statements from module scope to procedure scope.  Remove an
! unnecessary one.  Get rid of some array temps.  Cosmetic changes.
!
! Revision 2.1  2002/10/08 17:08:02  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.0  2001/09/17 20:26:26  livesey
! New forward model
!
! Revision 1.8.2.2  2001/09/10 23:57:40  zvi
! Import from nasty..
!
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
