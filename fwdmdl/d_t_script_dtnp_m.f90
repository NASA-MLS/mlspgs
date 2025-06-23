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
  public :: dT_Script_dT, dT_Script_dT_Sparse, dT_Script
!   public :: dT_Script_dT_W0, dT_Script_dF_W0

  interface dT_Script
    module procedure dT_Script_1D, dT_Script_2D
    module procedure dT_Script_1D_Sparse, dT_Script_2D_Sparse
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

  ! -----------------------------------------------  dT_Script_dT  -----
  subroutine dT_Script_dT ( T_Path, B, Eta_ZXP, NZ_ZXP, NNZ_ZXP, NU, &
    &                       DT_SCR_DT )
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

  end subroutine dT_Script_dT

  ! ----------------------------------------  dT_Script_dT_Sparse  -----
  subroutine dT_Script_dT_Sparse ( T_Path, B, Eta_ZXP, I_Start, I_Stop, &
                                 & NU, dT_Scr_dT, Skip )
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
    use Physics, only: H_O_K => h_over_k
    use Sparse_m, only: Sparse_t

! inputs

    real(rp), intent(in) :: t_path(:)      ! path temperatures K
    real(rp), intent(in) :: B(:)           ! Planck function
    class(sparse_t), intent(in) :: Eta_ZXP ! Interpolating coefficients from
                                           ! state vector to path
    integer, intent(in) :: I_Start, I_Stop ! Range of dBdx etc to use
    real(r8), intent(in) :: nu             ! calculation frequency (MHz)

! output

    real(rp), intent(out) :: dT_scr_dT(:,:) ! path dt_script temperature
                                          ! derivative

! optional

    integer, intent(in), optional :: Skip ! Process every skip'th row, default 1
                                          ! This allows, e.g., to process only
                                          ! elements on the coarse path, in
                                          ! which case Skip should == NGP1
! internals

    real(rp) :: dBdT(size(B))

    dBdT = B * ( h_o_k * nu + B ) / t_path**2
    call dT_script ( dBdT, eta_zxp, i_start, i_stop, dT_scr_dT, skip )

  end subroutine dT_Script_dT_Sparse

!   ! --------------------------------------------  dT_Script_dT_W0  -----
!   subroutine dT_Script_dT_W0 ( T_PATH, B, W0, dW0_dT, &
!     &                          ETA_ZXP, NZ_ZXP, NNZ_ZXP, &
!     &                          NU, DT_SCR_DT )
!     ! Build the derivative of the B array w.r.t to T_np
!     ! (The 'n' zeta coeff. and the 'p' phi coefficient)
! 
! !{ \parskip 5pt
! !  Given $T$, $\nu$, $\omega_0$, $\frac{\partial \omega_0}{\partial T}$,
! !  and the Planck function
! !  $B = \frac{\frac{h \nu}k}{exp\left(\frac{h \nu}{k T}\right)-1}$,
! !  compute \\
! !  $\frac{\partial (1-\omega_0) B}{\partial T} =
! !  -B \frac{\partial \omega_0}{\partial T} +
! !     (1-\omega_0) \frac{\partial B}{\partial T}$.
! !
! !  B satisfies the differential equation
! !  $\frac{\text{d} B}{\text{d} T} =
! !   \frac{B}{T^2} \left( \frac{h\nu}k + B \right)$, so\\
! !  $\frac{\partial (1-\omega_0) B}{\partial T} =
! !   B \left( \frac{1-\omega_0}{T^2} \left(\frac{h \nu}k + B \right) -
! !   \frac{\partial \omega_0}{\partial T} \right)$.
! !
! !  From $\frac{\partial (1-\omega_0) B}{\partial T}$ compute
! !  $\frac{\partial (1-\omega_0) B}{\partial T_{np}} =
! !   \frac{\partial (1-\omega_0) B}{\partial T} \eta_n \eta_p$ where
! !  $\eta_n$ = Eta(n,zeta) and $\eta_p$ = Eta(p,phi).
! !  Then difference $\frac{\partial (1-\omega_0) B}{\partial T_{np}}$ along
! !  the path to get {\tt dT\_scr\_dT} =
! !  $\frac{\partial \Delta (1-\omega_0) B}{\partial T}$, where
! !  $\Delta (1-\omega_0) B = \frac12 \left( (1-\omega_{0_{i+1}}) B_{i+1} -
! !                                          (1-\omega_{0_{i-1}}) B_{i-1} \right)$.
! !  $\Delta B$  is called ``T script'' (not Tau) in many notes and reports.
! 
!     use MLSCommon, only: R8, RP
!     use PHYSICS, only: H_O_K => h_over_k
! 
! ! inputs
! 
!     real(rp), intent(in) :: T_path(:)     ! path temperatures K
!     real(rp), intent(in) :: B(:)          ! Planck function
!     real(rp), intent(in) :: W0(:)         ! \omega_0
!     real(rp), intent(in) :: dW0_dT(:)     ! \frac{\partial \omega_0}{\partial T}
!     real(rp), intent(in) :: eta_zxp(:,:)  ! path eta functions
!     integer, intent(in) :: nz_zxp(:,:)    ! Where eta_zxp is not zero
!     integer, intent(in) :: nnz_zxp(:)     ! Numbers of rows in nz_zxp
!     real(r8), intent(in) :: nu            ! calculation frequency (MHz)
! 
! ! output
! 
!     real(rp), intent(out) :: dT_scr_dT(:,:) ! path dt_script temperature
! !                                           derivative
! 
! ! internals
! 
!     real(rp) :: dBdT(size(B))
! 
!     dBdT = B * ( ( 1.0_rp - w0 ) * &
!          &       ( h_o_k * nu + B ) / t_path**2 - dW0_dT )
!     call dT_script ( dBdT, eta_zxp, nz_zxp, nnz_zxp, dT_scr_dT )
! 
!   end subroutine dT_Script_dT_W0
! 
!   ! --------------------------------------------  dT_Script_dF_W0  -----
!   subroutine dT_Script_dF_W0 ( B, dW0_dX, ETA_ZXP, NZ_ZXP, NNZ_ZXP, &
!     &                          DT_SCR_DF )
!     ! Build the derivative of the B array w.r.t to mixing ratio F_np
!     ! (The 'n' zeta coeff. and the 'p' phi coefficient)
! 
! !{ \parskip 5pt
! !  Given $\frac{\partial \omega_0}{\partial f}$ and the Planck function
! !  $B = \frac{\frac{h \nu}k}{exp\left(\frac{h \nu}{k T}\right)-1}$,
! !  compute \\
! !  $\frac{\partial (1-\omega_0) B}{\partial f} =
! !  -B \frac{\partial \omega_0}{\partial f}$
! !  ($B$ does not depend upon $f$).
! !
! !  From $\frac{\partial (1-\omega_0) B}{\partial f}$ compute
! !  $\frac{\partial (1-\omega_0) B}{\partial f_{np}} =
! !   \frac{\partial (1-\omega_0) B}{\partial f} \eta_n \eta_p$ where
! !  $\eta_n$ = Eta(n,zeta) and $\eta_p$ = Eta(p,phi).
! !  Then difference $\frac{\partial (1-\omega_0) B}{\partial f_{np}}$ along
! !  the path to get {\tt dT\_scr\_df} =
! !  $\frac{\partial \Delta (1-\omega_0) B}{\partial f}$, where
! !  $\Delta (1-\omega_0) B = \frac12 \left( (1-\omega_{0_{i+1}}) B_{i+1} -
! !                                          (1-\omega_{0_{i-1}}) B_{i-1} \right)$.
! !  $\Delta B$  is called ``T script'' (not Tau) in many notes and reports.
! 
!     use MLSCommon, only: RP
! 
! ! inputs
! 
!     real(rp), intent(in) :: B(:)          ! Planck function
!     real(rp), intent(in) :: dW0_dX(:)     ! \frac{\partial \omega_0}{\partial X}
!     real(rp), intent(in) :: eta_zxp(:,:)  ! path eta functions
!     integer, intent(in) :: nz_zxp(:,:)    ! Where eta_zxp is not zero
!     integer, intent(in) :: nnz_zxp(:)     ! Numbers of rows in nz_zxp
! 
! ! output
! 
!     real(rp), intent(out) :: DT_SCR_DF(:,:) ! path dt_script / dT
! 
!     call dT_script ( -B * dW0_dX, eta_zxp, nz_zxp, nnz_zxp, dT_scr_df )
! 
!   end subroutine dT_Script_dF_W0

  ! -----------------------------------------------  dT_Script_1D  -----
  subroutine dT_Script_1D ( dBdx, ETA_ZXP, NZ_ZXP, NNZ_ZXP, DT_SCR, W0, dTScat )

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
    real(rp), intent(in) :: Eta_ZXP(:)    ! Interpolating coefficients from
                                          ! state vector to path
    integer, intent(in) :: Nz_ZXP(:)      ! Where eta_zxp is not zero
    integer, intent(in) :: Nnz_ZXP        ! Numbers of rows in nz_zxp
    real(rp), intent(in), optional :: W0(:)       ! On the path
    real(rp), intent(in), optional :: dTScat(:)   ! On the path, w.r.t.
                                          ! molecules in the grid.

! output

    real(rp), intent(out) :: dt_scr(:)    ! path dt_script / dx

! internals

    integer :: j, n_path, p_i

    real(rp) :: dT_x_eta(1:size(dBdx,1))

    n_path = size(dBdx)

    ! Make sure the elements of dT_x_eta that we use have defined values
    dT_x_eta(1:2) = 0.0
    dT_x_eta(n_path-1:n_path) = 0.0
    ! Fill other nonzero elements of dT_x_eta
    dT_x_eta(nz_zxp(:nnz_zxp)) = &
      & 0.5 * dbdx(nz_zxp(:nnz_zxp)) * eta_zxp(nz_zxp(:nnz_zxp))

    if ( present(w0) .and. present(dTScat) ) then
      dT_x_eta(nz_zxp(:nnz_zxp)) = &
        & dT_x_eta(nz_zxp(:nnz_zxp)) + &
        & 0.5 * w0(nz_zxp(:nnz_zxp)) * dTScat(nz_zxp(:nnz_zxp))
    end if

   dt_scr(1) = dT_x_eta(1) + dT_x_eta(2)
   dt_scr(n_path) = -dT_x_eta(n_path-1) - dT_x_eta(n_path)
! This is what would be happening if we didn't pay attention to the nonzeros:
!     dt_scr(2:n_path-1) = dT_x_eta(3:n_path) - dT_x_eta(1:n_path-2)
! But we pay attention to the nonzeros to improve performance.
! If we kept track of the nonzeros in dt_scr, we wouldn't neet the next line.
    dt_scr(2:n_path-1) = 0.0
    ! Now fill the nonzeros
    do j = 1, nnz_zxp
      p_i = nz_zxp(j)
      if ( p_i > 2 ) dt_scr(p_i-1) = dT_x_eta(p_i)
    end do
    do j = 1, nnz_zxp
      p_i = nz_zxp(j)
      if ( p_i >= n_path-1 ) exit
      dt_scr(p_i+1) = dt_scr(p_i+1) - dT_x_eta(p_i)
    end do

  end subroutine dT_Script_1D

  ! ----------------------------------------  dT_Script_1D_Sparse  -----
  subroutine dT_Script_1D_Sparse ( dBdx, Eta_ZXP, Column, I_Start, I_Stop, &
                                 & dT_SCR, Skip, &
                                 & W0, dTScat )

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
    use Sparse_m, only: Sparse_t

! inputs

    real(rp), intent(in) :: dBdx(:)       ! Derivative of Planck function
    class(sparse_t), intent(in) :: Eta_ZXP ! Interpolating coefficients from
                                          ! state vector to path
    integer, intent(in) :: Column         ! Column of Eta_ZXP to use
    integer, intent(in) :: I_Start, I_Stop  ! Range of dBdx etc to use

! output

    real(rp), intent(out) :: dt_scr(:)    ! path dt_script / dx

! optional

    integer, intent(in), optional :: Skip ! Process every skip'th row, default 1
                                          ! This allows, e.g., to process only
                                          ! elements on the coarse path, in
                                          ! which case Skip should == NGP1
    real(rp), intent(in), optional :: W0(:)       ! On the path
    real(rp), intent(in), optional :: dTScat(:)   ! On the path, w.r.t.
                                          ! molecules in the grid.


! internals

    integer :: j, mySkip, n_path, r, rs

    real(rp) :: dT_x_eta(1:size(dBdx,1))

    mySkip = 1
    if ( present(skip) ) mySkip = skip

    j = eta_zxp%cols(column) ! Last element in column c
    if ( j == 0 ) then
      dt_scr = 0
      return
    end if

    n_path = size(dBdx)

    ! Make sure the elements of dT_x_eta that we use have defined values
    dT_x_eta(1:2) = 0.0
    dT_x_eta(n_path-1:n_path) = 0.0
    ! Fill other nonzero elements of dT_x_eta
    do
      j = eta_zxp%e(j)%nc ! Element in next row of column c
      r = eta_zxp%e(j)%r  ! Row subscript of that element
      if ( mod(r-1,mySkip) == 0 ) then
        rs = (r-1) / skip + 1
        if ( rs >= i_start .and. rs <= i_stop ) &
          & dT_x_eta(rs) = 0.5 * dbdx(rs) * eta_zxp%e(j)%v
      end if
      if ( j == eta_zxp%cols(column) ) exit ! Processed last element in column
    end do

    if ( present(w0) .and. present(dTScat) ) then
      j = eta_zxp%cols(column) ! Last element in column c
      do
        j = eta_zxp%e(j)%nc ! Element in next row of column c
        r = eta_zxp%e(j)%r  ! Row subscript of that element
        if ( mod(r-1,mySkip) == 0 ) then
          rs = (r-1) / skip + 1
          if ( rs >= i_start .and. rs <= i_stop ) &
            & dT_x_eta(rs) = dT_x_eta(rs) + 0.5 * w0(rs) * dTScat(rs)
        end if
        if ( j == eta_zxp%cols(column) ) exit ! Processed last element in column
      end do
    end if

   dt_scr(1) = dT_x_eta(1) + dT_x_eta(2)
   dt_scr(n_path) = -dT_x_eta(n_path-1) - dT_x_eta(n_path)
! This is what would be happening if we didn't pay attention to the nonzeros:
!     dt_scr(2:n_path-1) = dT_x_eta(3:n_path) - dT_x_eta(1:n_path-2)
! But we pay attention to the nonzeros to improve performance.
! If we kept track of the nonzeros in dt_scr, we wouldn't neet the next line.
    dt_scr(2:n_path-1) = 0.0
    ! Now fill the nonzeros
    j = eta_zxp%cols(column) ! Last element in column c
    do
      j = eta_zxp%e(j)%nc ! Element in next row of column c
      r = eta_zxp%e(j)%r  ! Row subscript of that element
      if ( mod(r-1,mySkip) == 0 ) then
        rs = (r-1) / skip + 1
        if ( rs > 2 .and. rs >= i_start .and. rs <= i_stop ) &
          & dt_scr(rs-1) = dT_x_eta(rs)
      end if
      if ( j == eta_zxp%cols(column) ) exit ! Processed last element in column
    end do
    j = eta_zxp%cols(column) ! Last element in column c
    do
      j = eta_zxp%e(j)%nc ! Element in next row of column c
      r = eta_zxp%e(j)%r  ! Row subscript of that element
      if ( mod(r-1,mySkip) == 0 ) then
        rs = (r-1) / skip + 1
        if ( rs < n_path-1 .and. rs >= i_start .and. rs <= i_stop ) &
          & dt_scr(rs+1) = dt_scr(rs+1) - dT_x_eta(rs)
      end if
      if ( j == eta_zxp%cols(column) ) exit ! Processed last element in column
    end do

  end subroutine dT_Script_1D_Sparse

  ! -----------------------------------------------  dT_Script_2D  -----
  subroutine dT_Script_2D ( dBdx, ETA_ZXP, NZ_ZXP, NNZ_ZXP, DT_SCR, W0, dTScat )

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

    real(rp), intent(out) :: dt_scr(:,:)  ! path dt_script / dx

! internals

    integer :: sv_i

    if ( present(w0) .and. present(dTScat) ) then
      do sv_i = 1, size(eta_zxp,dim=2)
        call dt_script ( dBdx, eta_zxp(:,sv_i), nz_zxp(:,sv_i), nnz_zxp(sv_i), &
                       & dt_scr(:,sv_i), w0, dTscat(:,sv_i) )
      end do
    else
      do sv_i = 1, size(eta_zxp,dim=2)
        call dt_script ( dBdx, eta_zxp(:,sv_i), nz_zxp(:,sv_i), nnz_zxp(sv_i), &
                       & dt_scr(:,sv_i) )
      end do
    end if

  end subroutine dT_Script_2D

  ! ----------------------------------------  dT_Script_2D_Sparse  -----
  subroutine dT_Script_2D_Sparse ( dBdx, Eta_ZXP, I_Start, I_Stop, &
                                 & dT_Scr, Skip, W0, dTScat )

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
    use Sparse_m, only: Sparse_t

! inputs

    real(rp), intent(in) :: dBdx(:)        ! Derivative of Planck function
    class(sparse_t), intent(in) :: Eta_ZXP ! Interpolating coefficients from
                                           ! state vector to path
    integer, intent(in) :: I_Start, I_Stop ! Range of dBdx etc to use

! output

    real(rp), intent(out) :: dt_scr(:,:)  ! path dt_script / dx

! optional

    integer, intent(in), optional :: Skip
    real(rp), intent(in), optional :: W0(:)       ! On the path
    real(rp), intent(in), optional :: dTScat(:,:) ! Path X SV.  On the path,
                                          ! w.r.t. molecules in the grid.

! internals

    integer :: sv_i

    if ( present(w0) .and. present(dTScat) ) then
      do sv_i = 1, size(eta_zxp%cols)
        call dt_script ( dBdx, eta_zxp, sv_i, i_start, i_stop, dt_scr(:,sv_i), &
                       & skip, w0, dTscat(:,sv_i) )
      end do
    else
      do sv_i = 1, size(eta_zxp%cols)
        call dt_script ( dBdx, eta_zxp, sv_i, i_start, i_stop, dt_scr(:,sv_i), &
                       & skip )
      end do
    end if

  end subroutine dT_Script_2D_Sparse

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
! Revision 2.13  2011/03/23 23:52:30  vsnyder
! Add ability to do dTScat update in dT_Script
!
! Revision 2.12  2011/03/23 23:45:32  vsnyder
! This log entry is bogus.  Check in again to get the right one.
! FOV_Convolve_m.f90
!
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
