module HYDROSTATIC_INTRP
  use D_SOLVE_QUAD_M, only: SOLVE_QUAD
  use MLSCommon, only: I4, R8
  use D_HUNT_M, only: HUNT
  implicit NONE
  private
  public :: GET_HEIGHTS, GET_PRESSURES
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------     GET_HEIGHTS     -----
! Given an array of pressures, get set of pointing angles or heights:
!
  Subroutine GET_HEIGHTS ( ha, v_grid, t_grid, z_grid, n_i, p_out, v_out, &
 &                         n_o, ier )
!
! Given that: T(z)*dz = m*G*dh/K           (z=-Log(P))
! Given that: h =  D*Sin(a)   (Where D=Rs/(1 + ngrid), constant)
!             dh = D*Cos(a)*da
! Thus      : T(z)*dz = m*G*D*Cos(a)*da/K
! Assume that in the in the interval (z1,z2) we have:
!             T(z) = T1+w*(z-z1)     (w=(T2-T1)/(z2-z1))
! Then integrating the above over (z1,z2) gives (for delz=z2-z1):
!            delz*(T1+0.5*w*delz) = C*(Sin(a2)-Sin(a1))
! Where C=m*G*D/K, assume constant over the interval. Solving for C gives:
!       C=delz*(T1+0.5*w*delz)/(Sin(a2)-Sin(a1))
! Or:   C=delz*(T1+0.5*w*delz)/(h2-h1)
! Now, given 'z' inside (z1,z2), what is the 'a' ('a' inside (a1,a2)) ?
! let delz=z-z1, given C for that interval (z1,z2),
! Compute the integral: vl = delz*(T1+0.5*w*delz)/C
! And now h = h1 + vl   Or: a = a1 + vl
    character, intent(in) :: HA         ! ha='h' if v_grid=h_grid
                                        ! ha='a' if v_grid=ptg_angles
    real(r8), intent(in) :: V_GRID(*)   ! Input grid Angles or Heights
    real(r8), intent(in) :: T_GRID(*)   ! Input grid temperatures
    real(r8), intent(in) :: Z_GRID(*)   ! Input grid Log pressures
    integer(i4), intent(in) :: N_I      ! Number of input points
    real(r8), intent(in) :: P_OUT(*)    ! Output grid Log pressures
    real(r8), intent(out) :: V_OUT(*)   ! Output grid Angles or Heights
    integer(i4), intent(in) :: N_O      ! Number of output points
    integer(i4), intent(out) :: IER     ! Error flag
!
    character :: CH
    real(r8) :: COEFF(N_I), DZ
    real(r8), parameter :: EPS = 1.0e-5_r8
    integer :: I, KHI, KLO
    real(r8) :: P, RR, SA, T, V, VL, VR, W
!
    ier = 0
    ch = '@'
    if (ha == 'h' .or. ha == 'H') ch = 'h'
    if (ha == 'a' .or. ha == 'A') ch = 'a'
    if (ch == '@') then
      ier = 1
      Print *,'** Error in get_heights subroutine ..'
      Print *,'   variable: ha must be: "h" or "a"'
      Return
    end if
!
    call compute_coeff ( ch, n_i, v_grid, z_grid, t_grid, coeff )
!
    klo = -1
    do i = 1, n_o
      rr = p_out(i)
      Call hunt ( rr, z_grid, n_i, klo, khi )
      p = z_grid(klo)
      v = v_grid(klo)
      if (abs(rr-p) < eps) then
        v_out(i) = v
      else if (abs(rr-z_grid(khi)) < eps) then
        v_out(i) = v_grid(khi)
      else
        t = t_grid(klo)
        sa = v
        if (ch == 'a') sa = Sin(v)
        dz = z_grid(khi) - p
        w = (t_grid(khi) - t) / dz
        dz = rr - p
        vl = dz * (t + 0.5 * w * dz) / coeff(klo)
        vr = sa + vl
        if (ch == 'a') then
          if (abs(vr) <= 1.0) then
            vr = Asin(vr)
          else
            ier = 1
            Print *,'** Error in get_heights subroutine ..'
            Print *,'   Asking for ArcSin(arg) when abs(arg) > 1.0'
            Return
          end if
        end if
        v_out(i) = vr
      end if
    end do
!
    Return
  End Subroutine GET_HEIGHTS
!
!-------------------------------------------     GET_PRESSURES     -----
! Get an array of pressures from a given set of pointing angles or heights:
!
  Subroutine GET_PRESSURES ( ha, v_grid, t_grid, z_grid, n_i, v_out, p_out, &
 &                           n_o, ier )
!
! Given that: T(z)*dz = m*G*dh/K           (z=-Log(P))
! Given that: h =  H*Sin(a)   (Where H=Rs/(1 + ngrid), constant)
!             dh = H*Cos(a)*da
! Thus      : T(z)*dz = m*G*H*Cos(a)*da/K
! Assume that in the in the interval (z1,z2) we have:
!             T(z) = T1+w*(z-z1)     (w=(T2-T1)/(z2-z1))
! Then integrating the above over (z1,z2) gives (for delz=z2-z1):
!            delz*(T1+0.5*w*delz) = C*(Sin(a2)-Sin(a1))
! Where C=m*G*H/K, assume constant over the interval. Solving for C gives:
!       C=delz*(T1+0.5*w*delz)/(Sin(a2)-Sin(a1))
! Now, given a inside (a1,a2), what is the z (z inside (z1,z2)) ?
! let delz=z-z1, given C for that interval (z1,z2),
! solve the quadratic in delz: 0.5*w*delz*delz+T1*delz-C*(Sin(a)-Sin(a1))
! and get z=z1+delz which goes with the given a.
    character, intent(in) :: HA         ! ha='h' if v_grid=h_grid
                                        ! ha='a' if v_grid=ptg_angles
    real(r8), intent(in) :: V_GRID(*)   ! Input grid Angles or Heights
    real(r8), intent(in) :: T_GRID(*)   ! Input grid temperatures
    real(r8), intent(in) :: Z_GRID(*)   ! Input grid Log pressures
    integer(i4), intent(in) :: N_I      ! Number of input points
    real(r8), intent(in) :: V_OUT(*)    ! Output grid Angles or Heights
    real(r8), intent(out) :: P_OUT(*)   ! Output grid Log pressures
    integer(i4), intent(in) :: N_O      ! Number of output points
    integer(i4), intent(out) :: IER     ! Error flag
!
    real(r8) :: A, C
    character :: CH
    real(r8) :: COEFF(N_I)
    real(r8), parameter :: EPS = 1.0e-5_r8
    integer :: I, KHI, KLO
    real(r8) :: P, RR, S4, SA
    real(r8) :: T
    real(r8) :: V, W
    real(r8) :: X1, X2
!
    ier = 0
    ch = '@'
    if (ha == 'h' .or. ha == 'H') ch = 'h'
    if (ha == 'a' .or. ha == 'A') ch = 'a'
    if (ch == '@') then
      ier = 1
      Print *,'** Error in get_pressures subroutine ..'
      Print *,'   variable: ha must be: "h" or "a"'
      Return
    end if
!
    call compute_coeff ( ch, n_i, v_grid, z_grid, t_grid, coeff )
!
! For each output angle, solve for the zeta.
! In solving the quadratic equation, choose the smaller root always.
!
    klo = -1
    do i = 1, n_o
      rr = v_out(i)
      Call hunt ( rr, v_grid, n_i, klo, khi )
      v = v_grid(klo)
      p = z_grid(klo)
      if (abs(rr-v) < eps) then
        p_out(i) = p
      else if (abs(rr-v_grid(khi)) < eps) then
        p_out(i) = z_grid(khi)
      else
        t = t_grid(klo)
        sa = v
        if (ch == 'a') sa = Sin(sa)
        w = (t_grid(khi) - t) / (z_grid(khi) - p)
        a = 0.5 * w
        s4 = rr
        if (ch == 'a') s4 = Sin(rr)
        c = -coeff(klo) * (s4 - sa)
        Call Solve_Quad ( a, t, c, x1, x2 )
        p_out(i) = p + x2            ! Take the smaller root always
      end if
    end do
!
    Return
  End Subroutine GET_PRESSURES
! *****     Private procedure     **************************************
! ------------------------------------------     COMPUTE_COEFF     -----
  subroutine COMPUTE_COEFF ( ch, n_i, v_grid, z_grid, t_grid, coeff )
    character, intent(in) :: CH         ! ch='h' if v_grid=h_grid
                                        ! ch='a' if v_grid=ptg_angles
    integer(i4), intent(in) :: N_I      ! Number of input points
    real(r8), intent(in) :: V_GRID(*)   ! Input grid Angles or Heights
    real(r8), intent(in) :: T_GRID(*)   ! Input grid temperatures
    real(r8), intent(in) :: Z_GRID(*)   ! Input grid Log pressures
    real(r8), intent(out) :: COEFF(N_I)
    real(r8) :: DZ
    integer :: I
    real(r8) :: SAI, SAIP1, VL, W
!
! Computes the coefficients (C) for each sub-interval
!
    saip1 = v_grid(1)
    if (ch == 'a') saip1 = Sin(saip1)
    do i = 1, n_i - 1
      sai = saip1
      saip1 = v_grid(i+1)
      if (ch == 'a') saip1 = Sin(saip1)
      dz = z_grid(i+1) - z_grid(i)
      w = (t_grid(i+1) - t_grid(i)) / dz   ! ??? Do we really want to divide
      vl = dz * (t_grid(i) + 0.5 * w * dz) ! ??? by DZ, and then multiply?
      coeff(i) = vl / (saip1 - sai)
    end do
  end subroutine COMPUTE_COEFF
end module HYDROSTATIC_INTRP
! $Log$
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
!
