module D_T_SCRIPT_DTNP_M

  implicit NONE
  private
  public :: DT_SCRIPT_DT

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = &
    & "$RCSfile$"
!---------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------
! Build the derivative of the t_script array w.r.t to Tnp
! (The 'n' zeta coeff. and the 'p' phi coefficient)

  ! -----------------------------------------------  DT_SCRIPT_DT  -----
  subroutine DT_SCRIPT_DT ( T_PATH, ETA_ZXP, NU, DT_SCR_DT )

! Let h_o_k = h / k = (Planck's constant) / (Boltzmann's constant)
! Let hxf = h_o_k * f  (f = Frq.)
! Let B = stat_temp(T,f) = hxf / (exp(hxf/T) - 1.0)
! Then: dB/dT = ((B/T)**2)*exp(hxf/T)
! And : dB/dTnp = dB/dT * Etan * Etap
! Where: Etan = Eta(n,zeta), Etap = Eta(p,phi)

    use MLSCommon, only: R8, RP, IP
    use PHYSICS, only: H_O_K => h_over_k

! inputs

    real(rp), intent(in) :: t_path(:)     ! path temperatures K
    real(rp), intent(in) :: eta_zxp(:,:)  ! path eta functions
!    logical, intent(in) :: not_zero(:,:) ! where eta is not zero
    real(r8) :: nu                        ! calculation frequency (MHz)

! output

    real(rp), intent(out) :: dt_scr_dt(:,:) ! path dt_script temparature
!                                           derivative
! internals

    integer(ip) :: n_path, n_sv, sv_i

    real(rp) :: dstdt(1:size(t_path))
    real(rp) :: a, b

! Do this inefficiently for now because it is quick and easy

    n_path = size(t_path)
    n_sv = size(eta_zxp,dim=2)

! Avoid some array temps by not using whole-array operations

!{ Compute $\frac{\text{d} B}{\text{d} T} =
!   \exp \left ( \frac{h \nu}{k T} \right )
!   \left ( \frac{h \nu}
!     {k T \left [ \exp \left ( \frac{h \nu}{k T} \right ) - 1 \right ]} \right ) ^2$

    do sv_i = 1, n_path
      a = h_o_k * nu / t_path(sv_i)        !{ \frac{h \nu}{k T}
      b = exp(a)                           !{ \exp \left ( \frac{h \nu}{k T} \right )
      dstdt(sv_i) = b * (a/(b-1.0_rp))**2  !{ \frac{\text{d} B}{\text{d} T}
    end do

! This part Zvi will hate

    do sv_i = 1 , n_sv
      a = -dstdt(1) * eta_zxp(1,sv_i)
      b = -dstdt(n_path) * eta_zxp(n_path,sv_i)
      dt_scr_dt(:,sv_i) = 0.5_rp * (  &
          &     eoshift(dstdt * eta_zxp(:,sv_i), 1, b) -   &
          &     eoshift(dstdt * eta_zxp(:,sv_i),-1, a) )
    end do

  end subroutine DT_SCRIPT_DT

  ! ----------------------------------------------  NOT_USED_HERE  -----
  logical function NOT_USED_HERE()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function NOT_USED_HERE

end module D_T_SCRIPT_DTNP_M
! $Log$
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
