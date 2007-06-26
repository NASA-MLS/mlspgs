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
  public :: DT_SCRIPT_DT

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------
! Build the derivative of the t_script array w.r.t to Tnp
! (The 'n' zeta coeff. and the 'p' phi coefficient)

  ! -----------------------------------------------  DT_SCRIPT_DT  -----
  subroutine DT_SCRIPT_DT ( T_PATH, T_SCRIPT, ETA_ZXP, NZ_ZXP, NNZ_ZXP, &
    &                       NU, DT_SCR_DT )

!{ \parskip 5pt
!  Given $T$, $\nu$ and
!  $B = \frac{\frac{h \nu}k}{exp\left(\frac{h \nu}{k T}-1\right)}$,
!  compute $\frac{\text{d} B}{\text{d} T}$.  $B$ is called ``T script''
!  in many notes and reports.
!
! B satisfies the differential equation
! $\frac{\text{d} B}{\text{d} T} = B \frac{\frac{h\nu}k + B}{T^2}$.
!
! The old way was to compute $\frac{\text{d} B}{\text{d} T} =
!   \exp \left ( \frac{h \nu}{k T} \right )
!   \left ( \frac{h \nu}
!     {k T \left [ \exp \left ( \frac{h \nu}{k T} \right ) - 1 \right ]}
!     \right ) ^2$.  Dividing numerator and denominator by
!     $\exp \left ( \frac{h \nu}{k T} \right )$ gives
!     $\frac{\left ( \frac{h \nu}{k T} \right )^2}
!           { 2 \left [ \cosh \left ( \frac{h \nu}{k T} \right ) - 1 \right ]}$.
!   $2 ( \cosh x \, - \, 1) / x^2$ has substantial cancellation near $x = 0$,
!   so a specially-developed procedure should be used to evaluate it.
!
! From $\frac{\text{d} B}{\text{d} T}$ compute
! $\frac{\text{d} B}{\text{d} T_{np}} =
!  \frac{\text{d} B}{\text{d} T} \eta_n \eta_p$ where
! $\eta_n$ = Eta(n,zeta) and $\eta_p$ = Eta(p,phi).

!   use DCOSH1_M, only: COSH1             ! In case RP is double
    use MLSCommon, only: R8, RP, IP
    use PHYSICS, only: H_O_K => h_over_k
!   use SCOSH1_M, only: COSH1             ! In case RP is single

! inputs

    real(rp), intent(in) :: t_path(:)     ! path temperatures K
    real(rp), intent(in) :: t_script(:)   ! B, called "T script" in many notes
    real(rp), intent(in) :: eta_zxp(:,:)  ! path eta functions
    integer, intent(in) :: nz_zxp(:,:)    ! Where eta_zxp is not zero
    integer, intent(in) :: nnz_zxp(:)     ! Numbers of rows in nz_zxp
    real(r8) :: nu                        ! calculation frequency (MHz)

! output

    real(rp), intent(out) :: dt_scr_dt(:,:) ! path dt_script temperature
!                                           derivative
! internals

    integer(ip) :: j, n_path, n_sv, p_i, sv_i

    real(rp) :: dstdt(1:size(t_path))     ! d "script T" / dT
    real(rp) :: dT_x_eta(1:size(t_path),1:size(eta_zxp,2))

! Do this inefficiently for now because it is quick and easy

    n_path = size(t_path)
    n_sv = size(eta_zxp,dim=2)

!   dstdt = 0.5_rp / cosh1(h_o_k * nu / t_path) !{ \frac{a^2}{2 ( \cosh a - 1 )}
    dstdt = 0.5 * t_script * ( h_o_k * nu + t_script ) / t_path**2

    ! Make sure the elements of dT_x_eta that we use have defined values
    dT_x_eta(1:2,:) = 0.0
    dT_x_eta(n_path-1:n_path,:) = 0.0
    ! Fill other nonzero elements of dT_x_eta
    do sv_i = 1, n_sv
      dT_x_eta(nz_zxp(:nnz_zxp(sv_i),sv_i),sv_i) = &
        & dstdt(nz_zxp(:nnz_zxp(sv_i),sv_i)) * eta_zxp(nz_zxp(:nnz_zxp(sv_i),sv_i),sv_i)
    end do

    do sv_i = 1 , n_sv
      dt_scr_dt(1,sv_i) = dT_x_eta(1,sv_i) + dT_x_eta(2,sv_i)
      dt_scr_dt(n_path,sv_i) = -dT_x_eta(n_path-1,sv_i) - dT_x_eta(n_path,sv_i)
! This is what would be happening if we didn't pay attention to the nonzeros:
!     dt_scr_dt(2:n_path-1,sv_i) = dT_x_eta(3:n_path,sv_i) - &
!         &                        dT_x_eta(1:n_path-2,sv_i)
! But we pay attention to the nonzeros to improve performance.
! If we kept track of the nonzeros in dt_scr_dt, we wouldn't neet the next line.
      dt_scr_dt(2:n_path-1,sv_i) = 0.0
      ! Now fill the nonzeros in column sv_i
      do j = 1, nnz_zxp(sv_i)
        p_i = nz_zxp(j,sv_i)
        if ( p_i > 2 ) dt_scr_dt(p_i-1,sv_i) = dT_x_eta(p_i,sv_i)
      end do
      do j = 1, nnz_zxp(sv_i)
        p_i = nz_zxp(j,sv_i)
        if ( p_i >= n_path-1 ) exit
        dt_scr_dt(p_i+1,sv_i) = dt_scr_dt(p_i+1,sv_i) - dT_x_eta(p_i,sv_i)
      end do
    end do

  end subroutine DT_SCRIPT_DT

  ! ----------------------------------------------  NOT_USED_HERE  -----
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module D_T_SCRIPT_DTNP_M
! $Log$
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
