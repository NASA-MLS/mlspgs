module D_T_SCRIPT_DTNP_M
  use D_STAT_TEMP_M, only: STAT_TEMP
  use MLSCommon, only: R8, RP, IP
  use PHYSICS, only: H_O_K => h_over_k
  implicit NONE
  private
  public :: DT_SCRIPT_DT
!
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
   "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
   "$RCSfile$"
!---------------------------------------------------------------------------
contains
!----------------------------------------------------------------------
! This routine builds the derivative of the t_script array w.r.t to Tnp
! (The 'n' zeta coeff. and the 'p' phi coefficient)

  SUBROUTINE dt_script_dt(t_path,eta_zxp,nu,dt_scr_dt)
!
! Let h_o_k = h / k = (Planck's constant) / (Boltzmann's constant)
! Let hxf = h_o_k * f  (f = Frq.)
! Let B = stat_temp(T,f) = hxf / (exp(hxf/T) - 1.0)
! Then: dB/dT = ((B/T)**2)*exp(hxf/T)
! And : dB/dTnp = dB/dT * Etan * Etap
! Where: Etan = Eta(n,zeta), Etap = Eta(p,phi)
!
! inputs
!
  REAL(rp), INTENT(in) :: t_path(:) ! path temperatures K
  REAL(rp), INTENT(in) :: eta_zxp(:,:) ! path eta functions
!  LOGICAL, INTENT(in) :: not_zero(:,:) ! where eta is not zero
  REAL(r8) :: nu ! calculation frequency (MHz)
!
! output
!
  REAL(rp), INTENT(out) :: dt_scr_dt(:,:) ! path dt_script temparature
!                                           derivative
! internals
!
  INTEGER(ip) :: n_path, n_sv, sv_i

  REAL(rp), ALLOCATABLE :: dstdt(:)
  REAL(rp) :: a, b
!
! I will do this inefficiently for now because it is quick and easy
!
  n_path = SIZE(t_path)
  n_sv = SIZE(eta_zxp,DIM=2)
!
  ALLOCATE(dstdt(1:n_path))
!
  dstdt = EXP(h_o_k * nu / t_path)
  dstdt = dstdt * (h_o_k * nu / (t_path * (dstdt - 1.0_rp)))**2
!
! This part Zvi will hate
!
  DO sv_i = 1 , n_sv
    a = -dstdt(1) * eta_zxp(1,sv_i)
    b = -dstdt(n_path) * eta_zxp(n_path,sv_i)
    dt_scr_dt(:,sv_i) = 0.5_rp * (  &
        &     EOSHIFT(dstdt * eta_zxp(:,sv_i), 1, b) -   &
        &     EOSHIFT(dstdt * eta_zxp(:,sv_i),-1, a) )
  ENDDO
!
  DEALLOCATE(dstdt)
!
  END SUBROUTINE dt_script_dt
!
END module D_T_SCRIPT_DTNP_M
! $Log$
! Revision 1.8.2.2  2001/09/10 23:57:40  zvi
! Import from nasty..
!
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
