module D_T_SCRIPT_DTNP_M
  use D_STAT_TEMP_M, only: STAT_TEMP
  use MLSCommon, only: I4, R8
  use PHYSICS, only: H_O_K => h_over_k
  use PATH_ENTITIES_M, only: PATH_VECTOR
  use D_GET_ONE_ETA_M, only: GET_ONE_ETA
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
  Subroutine D_T_SCRIPT_DTNP ( Frq, t_basis, phi_basis, brkpt, no_ele,  &
 &           z_path, t_path, phi_path, Ng, in, ip, no_t, no_phi,   &
 &           dt_script_dtnp )

    real(r8), intent(in) :: FRQ
    real(r8), intent(in) :: T_BASIS(:), PHI_BASIS(:)

    integer(i4), intent(in) :: IN,IP,NO_T,NO_PHI,Ng,brkpt,no_ele

    Type(path_vector), intent(in) :: T_PATH, Z_PATH, PHI_PATH

    real(r8), intent(out) :: DT_SCRIPT_DTNP(:)

    integer(i4) :: i, j, m, mid, Ngp1
    real(r8) :: hxf, s, t1, t2, etan, etap
    real(r8) :: dt_dnp(size(dt_script_dtnp))
!
    hxf = h_o_k * Frq
!
!  Initialize the array:
!
    dt_dnp(:) = 0.0
    dt_script_dtnp(:) = 0.0

! 'brkpt' is the index of the path break-point (when it change from
!         incoming ray to outgoing ray)
! 'no_ele' is the total number of entries in t_path%values(1...no_ele)

    m = 0
    Ngp1 = Ng + 1
    j = 1 - Ngp1

    do
      j = j + Ngp1
      if (j > brkpt) EXIT
      m = m + 1
      Call get_one_eta(z_path%values(j),t_basis,no_t,in,etan)
      if(etan > 0.0) then
        Call get_one_eta(phi_path%values(j),phi_basis,no_phi,ip,etap)
        if(etap > 0.0) then
          t1 = t_path%values(j)
          t2 = stat_temp(t1,frq)
          s = t2 / t1
          dt_dnp(m) = s * s * exp(hxf/t1) * etan * etap
        endif
      endif
    end do
!
    mid = m
    j = brkpt + 1 - Ngp1

    do
      j = j + Ngp1
      if (j > no_ele) EXIT
      m = m + 1
      Call get_one_eta(z_path%values(j),t_basis,no_t,in,etan)
      if(etan > 0.0) then
        Call get_one_eta(phi_path%values(j),phi_basis,no_phi,ip,etap)
        if(etap > 0.0) then
          t1 = t_path%values(j)
          t2 = stat_temp(t1,frq)
          s = t2 / t1
          dt_dnp(m) = s * s * exp(hxf/t1) * etan * etap
        endif
      endif
    end do
!
    dt_script_dtnp(1) = 0.5_r8 * (dt_dnp(1) + dt_dnp(2))
    do i = 2, mid-1
      dt_script_dtnp(i) = 0.5_r8 * (dt_dnp(i+1) - dt_dnp(i-1))
    end do
!
    dt_script_dtnp(mid  ) = 0.5_r8 * (dt_dnp(mid+1) - dt_dnp(mid-1))
    dt_script_dtnp(mid+1) = 0.5_r8 * (dt_dnp(mid+2) - dt_dnp(mid))
!
! Note that dt_dnp(mid) = dt_dnp(mid+1) because this index
! is redundant at the tangent vertical (or center phi)
!
    do i = mid+2, m-1
      dt_script_dtnp(i) = 0.5_r8 * (dt_dnp(i+1) - dt_dnp(i-1))
    end do
!
    dt_script_dtnp(m) = -0.5_r8 * (dt_dnp(m-1) + dt_dnp(m))
!
    Return
  End subroutine D_T_SCRIPT_DTNP
end module D_T_SCRIPT_DTNP_M
! $Log$
! Revision 1.6  2001/03/31 23:40:55  zvi
! Eliminate l2pcdim (dimension parameters) move to allocatable ..
!
! Revision 1.5  2001/03/29 08:51:01  zvi
! Changing the (*) toi (:) everywhere
!
! Revision 1.4  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
