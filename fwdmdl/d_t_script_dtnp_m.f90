module D_T_SCRIPT_DTNP_M
  use D_STAT_TEMP_M, only: STAT_TEMP
  use L2PCdim, only: N2LVL
  use MLSCommon, only: I4, R4, R8
  use PHYSICS, only: H_O_K => h_over_k
  use S_GET_ONE_ETA_M, only: GET_ONE_ETA
  implicit NONE
  private
  public :: D_T_SCRIPT_DTNP

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
!
! Let h_o_k = h / k = (Planck's constant) / (Boltzmann's constant)
! Let hxf = h_o_k * f  (f = Frq.)
! Let B = stat_temp(T,f) = hxf / (exp(hxf/T) - 1.0)
! Then: dB/dT = ((B/T)**2)*exp(hxf/T)
! And : dB/dTnp = dB/dT * Etan * Etap
! Where: Etan = Eta(n,zeta), Etap = Eta(p,phi)
!
! The begining of the ray has index 1, and the end of the ray has index: 2*n.
! The uppermost index of the ray on the right is: IndxR, and the
! lowermost index ot the ray on the left is: IndxL
!
  Subroutine D_T_SCRIPT_DTNP ( Frq, t_basis, phi_basis, t_path, z_path,    &
 &           phi_path, N_lvls, Ng, path_brkpt, IndxR, IndxL, in, ip, no_t, &
 &           no_phi, dt_script_dtnp )

    real(r8), intent(in) :: FRQ
    real(r4), intent(in) :: T_BASIS(*)
    real(r4), intent(in) :: PHI_BASIS(*)
    real(r4), intent(in) :: T_PATH(*)
    real(r4), intent(in) :: Z_PATH(*)
    real(r4), intent(in) :: PHI_PATH(*)
    integer(i4), intent(in) :: N_lvls
    integer(i4), intent(in) :: Ng
    integer(i4), intent(in) :: PATH_BRKPT(*)
    integer(i4), intent(in) :: INDXR, INDXL
    integer(i4), intent(in) :: IN
    integer(i4), intent(in) :: IP
    integer(i4), intent(in) :: NO_T
    integer(i4), intent(in) :: NO_PHI
    real(r4), intent(out) :: DT_SCRIPT_DTNP(*)

    Integer(i4) :: I, J, K, NGP1
    real(r8) :: DPST, HXF, S, T1, T2, DT_DNP(N2lvl)
    real(r4) :: ETAN, ETAP
!
    hxf = h_o_k * Frq
!
!  Initialize the array:
!
    dt_dnp(1:N2lvl) = 0.0
    dt_script_dtnp(1:N2lvl) = 0.0
!
    i = 0
    j = -Ng
    Ngp1 = Ng + 1
    do while (i < IndxR .and. j < path_brkpt(1) )
      i = i + 1
      j = j + Ngp1
      Call get_one_eta(z_path(j),t_basis,no_t,in,etan)
      if(etan > 0.0) then
        Call get_one_eta(phi_path(j),phi_basis,no_phi,ip,etap)
        if(etap > 0.0) then
          t1 = dble(t_path(j))
          t2 = stat_temp(t1,frq)
          s = t2 / t1
          dt_dnp(i) = s * s * exp(hxf/t1) * dble(etan) * dble(etap)
        endif
      endif
    end do
!
    k = 2 * N_lvls
    j = path_brkpt(2) - Ngp1
    do i = IndxL, k
      j = j + Ngp1
      Call get_one_eta(z_path(j),t_basis,no_t,in,etan)
      if(etan > 0.0) then
        Call get_one_eta(phi_path(j),phi_basis,no_phi,ip,etap)
        if(etap > 0.0) then
          t1 = dble(t_path(j))
          t2 = stat_temp(t1,frq)
          s = t2 / t1
          dt_dnp(i) = s * s * exp(hxf/t1) * dble(etan) * dble(etap)
        endif
      endif
    end do
!
    dt_script_dtnp(1) = 0.5_r8 * (dt_dnp(1) + dt_dnp(2))
    do i = 2, IndxR-1
      dt_script_dtnp(i) = 0.5_r8 * (dt_dnp(i+1) - dt_dnp(i-1))
    end do
!
    dt_script_dtnp(IndxR) = 0.5_r8 * (dt_dnp(IndxL) - dt_dnp(IndxR-1))
    dt_script_dtnp(IndxL) = 0.5_r8 * (dt_dnp(IndxL+1) - dt_dnp(IndxR))
!
! Note that dt_dnp(IndxR) = dt_dnp(IndxL) because this index
! is redundant at the tangent vertical (or center phi)
!
    do i = IndxL+1, k-1
      dt_script_dtnp(i) = 0.5_r8 * (dt_dnp(i+1) - dt_dnp(i-1))
    end do
!
    dt_script_dtnp(k) = 0.5_r8 * (dt_dnp(k-1) + dt_dnp(k))
!
    Return
  End subroutine D_T_SCRIPT_DTNP
end module D_T_SCRIPT_DTNP_M

! $Log$
