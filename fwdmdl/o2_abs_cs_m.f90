! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module O2_Abs_CS_M

  implicit NONE
  private
  public :: O2_Abs_CS, Get_QN_By_Frequency

  integer :: O2_in_catalog = -1   ! Index in Spectroscopy catalog of O2 line.

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    &  "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

! ----------------------------------------------------  O2_Abs_CS  -----
  subroutine O2_Abs_CS ( Freq, Qn, H, Slabs_0, Polarized, &
    &                    Slabs_nonres, Y_nonres, &
    &                    Sigma_p, Pi, Sigma_m )

! Compute the complex absorption cross section.
! Modified to use Voigt with interfered lineshape

    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
    use MLSCommon, only: IP, R8, Rk => RP
    use Units, only: SqrtPi

    real(r8), intent(in) :: Freq              ! Observation frequency
    integer(ip), intent(in) :: Qn(:)          ! Quantum numbers
    real(rk), intent(in) :: H                 ! Magnetic field
    type(slabs_struct), intent(in) :: Slabs_0 ! contains, among others:
!    v0s(:)         ! pressure shifted line centers
!    x1(:)          ! Doppler width
!    y(:)           ! ratio Pressure to Doppler widths
!    yi(:)          ! Interference coefficients
!    slabs1(:)      ! strengths
    logical, intent(in) :: Polarized(:)       ! Which lines to use.  Same
                                              ! size as Slabs_0%v0s.
    real(rk), intent(in) :: Slabs_nonres      ! Nonresonant absorption
    real(rk), intent(in) :: Y_nonres          ! Nonresonant ratio
    complex(rk), intent(out) :: Sigma_p       ! Output Beta values
    complex(rk), intent(out) :: Pi
    complex(rk), intent(out) :: Sigma_m

    integer(ip) :: J, No_lines

    real(rk) :: F_o_v0, X, Z, U, V, Denomm

    complex(rk) :: Wing

    no_lines = size(polarized)

    sigma_p = 0.0_rk
    pi = 0.0_rk
    sigma_m = 0.0_rk
    wing = 0.0_rk

! Do magnetic calculation even if h = 0 for all of the Zeeman-split lines

    do j = 1, no_lines

      if ( .not. polarized(j) ) cycle

      call mag_o2_abs_cs ( qn(j), freq, slabs_0%x1(j), slabs_0%slabs1(j),  &
        &                  slabs_0%y(j), slabs_0%yi(j), h, slabs_0%v0s(j), &
        &                  sigma_p, pi, sigma_m )

! Fill in negative frequency part of VVW lineshape for jth line

      f_o_v0 = freq / slabs_0%v0s(j)
      z = slabs_0%x1(j) * (slabs_0%v0s(j) + freq)
      denomm = sqrtPi * (z*z + slabs_0%y(j)*slabs_0%y(j))
      wing = wing + ( 0.5_rk * slabs_0%slabs1(j) * f_o_v0 ) * cmplx( &
        & f_o_v0 * (slabs_0%y(j) - slabs_0%yi(j)*z) / denomm, & ! Real part
        & (z + slabs_0%y(j) * (slabs_0%y(j)/slabs_0%x1(j) -   & ! Imaginary part...
        &   freq*slabs_0%yi(j)) / slabs_0%v0s(j)) / denomm -  &
        &   2.0_rk/(slabs_0%x1(j) * slabs_0%v0s(j)) &
        & )

    end do

    sigma_p = sigma_p + 0.5_rk * wing
    pi = pi + wing
    sigma_m = sigma_m + 0.5_rk * wing

! Contribution from non-Zeeman-split lines is done in get_beta_path_scalar.

! Non resonant contribution is done in get_beta_path_scalar too.

  end subroutine O2_Abs_CS

! ------------------------------------------------  Mag_O2_Abs_CS  -----
  subroutine Mag_O2_Abs_CS ( n, nu, x1, s, w, y, h, v0, sigma_p, pi, sigma_m )

! Compute the frequency dependent absorption cross section for magnetic o2.

! Other notes:
! Document refers to "MLS Spectroscopic Data Base" W. G. Read, Version 1.0
! October 19, 1990.

    use MLSCommon, only: IP, R8, Rk => RP
    use SLabs_SW_M, only: Simple_Voigt

    integer, intent(in) :: N    ! rotational quantum number, sign indicates delta J
    real(r8), intent(in) :: Nu  ! transmission frequency in MHz
    real(rk), intent(in) :: X1  ! Doppler width factor sqrt(ln 2) / D_width
    real(rk), intent(in) :: S   ! strength factor slabs1 from slabs routine
    real(rk), intent(in) :: W   ! collision to doppler width ratio in the
                                ! document corrected for temperature and pressure
    real(rk), intent(in) :: Y   ! interference coefficient in the document
                                ! corrected for pressure and temperature
    real(rk), intent(in) :: H   ! magnetic field in Gauss
    real(rk), intent(in) :: V0  ! zero magnetic field line position

    complex(rk), intent(inout) :: Sigma_P ! absorption coefficient at unity mixing
                                ! ratio for Delta M = +1
    complex(rk), intent(inout) :: Pi ! absorption coefficient at unity mixing
                                ! ratio for Delta M = 0
    complex(rk), intent(inout) :: Sigma_M ! absorption coefficient at unity mixing
                                ! ratio for Delta M = -1
    integer(ip) :: M

    real(rk) :: Nu_offst, Del_nu, Xi, Denom1, Denom2, Z, F_o_v0, U, V, Zr, Zi

! Compute the absorption coefficient at unity mixing ratio for N transition

    if ( n == 0 ) return

    f_o_v0 = nu / v0
    del_nu = v0 - nu

    if ( n > 0 ) then

! Delta J = +1

      denom1 = n * (n + 1)
      denom2 = (n + 1) * (2 * n +1) * (2 * n + 3)

      do m = -n , n

! pi transition

        nu_offst = x1*(del_nu + 2.803_rk*m*(1-n)*h / denom1)
        call simple_voigt ( nu_offst, w, u, v )
        xi = 3.0_rk * ((n + 1) * (n + 1) - m * m) / denom2
        z  = 0.5_rk * s * xi * f_o_v0
        zr = z * f_o_v0 * (u - y*v)
        zi = z * (v + u * f_o_v0 * (w / (x1*nu) + y))
        pi = pi + cmplx(zr, zi)

! sigma_p transition

        nu_offst = x1*(del_nu + 2.803_rk*(m*(1-n)-n)*h / denom1)
        call simple_voigt ( nu_offst, w, u, v )
        xi = 3.0_rk * (n + m + 1) * (n + m + 2) / (4.0_rk * denom2)
        z  = 0.5_rk * s * xi * f_o_v0
        zr = z * f_o_v0 * (u - y * v)
        zi = z * (v + u * f_o_v0 * (w/(x1 * nu) + y))
        sigma_p = sigma_p + cmplx(zr, zi)

! sigma_m transition

        nu_offst = x1*(del_nu + 2.803_rk*(m*(1-n)+n)*h / denom1)
        call simple_voigt ( nu_offst, w, u, v )
        xi = 3.0_rk * (n - m + 1) * (n - m + 2) / (4.0_rk * denom2)
        z  = 0.5_rk * s * xi * f_o_v0
        zr = z * f_o_v0 * (u - y * v)
        zi = z * (v + u * f_o_v0 * (w/(x1*nu) + y))
        sigma_m = sigma_m + cmplx(zr, zi)

      end do

    else

! Delta J = -1 n is a negative number

      denom1 =  n * (n - 1)
      denom2 = -n * (4 * n * n - 1)

      do m = n+2 , -(n+2)

! pi transition

        nu_offst = x1*(del_nu + 2.803_rk*m*(2-n)*h / denom1)
        call simple_voigt ( nu_offst, w, u, v )
        xi = 3.0_rk * (n * n - m * m) / denom2
        z  = 0.5_rk * s * xi * f_o_v0
        zr = z * f_o_v0 * (u - y*v)
        zi = z * (v + u * f_o_v0 * (w/(x1*nu) + y))
        pi = pi + cmplx(zr, zi)

! sigma_p transition

        nu_offst = x1*(del_nu + 2.803_rk*((m+1)*(2-n)-1)*h / denom1)
        call simple_voigt ( nu_offst, w, u, v )
        xi = 3.0_rk * (n + m) * (n + m + 1) / (4.0_rk * denom2)
        z  = 0.5_rk * s * xi * f_o_v0
        zr = z * f_o_v0 * (u - y*v)
        zi = z * (v + u * f_o_v0 * (w/(x1*nu) + y))
        sigma_p = sigma_p + cmplx(zr, zi)

! sigma_m transition

        nu_offst = x1*(del_nu + 2.803_rk*((m-1)*(2-n)+1)*h / denom1)
        call simple_voigt ( nu_offst, w, u, v )
        xi = 3.0_rk * (m - n) * (m - n - 1) / (4.0_rk * denom2)
        z  = 0.5_rk * s * xi * f_o_v0
        zr = z * f_o_v0 * (u - y*v)
        zi = z * (v + u * f_o_v0 * (w/(x1*nu) + y))
        sigma_m = sigma_m + cmplx(zr, zi)

      end do

! Finish off special transition cases for pi and sigmas
! m = -n

      m = n

! sigma_p transition

      nu_offst = x1*(del_nu + 2.803_rk*((m+1)*(2-n)-1)*h / denom1)
      call simple_voigt ( nu_offst, w, u, v )
      xi = 3.0_rk * (n + m) * (n + m + 1) / (4.0_rk * denom2)
      z  = 0.5_rk * s * xi * f_o_v0
      zr = z * f_o_v0 * (u - y*v)
      zi = z * (v + u * f_o_v0 * (w/(x1*nu) + y))
      sigma_p = sigma_p + cmplx(zr, zi)

! m = n

      m = -n

! sigma_m transition

      nu_offst = x1*(del_nu + 2.803_rk*((m-1)*(2-n)+1)*h / denom1)
      call simple_voigt ( nu_offst, w, u, v )
      xi = 3.0_rk * (m - n) * (m - n - 1) / (4.0_rk * denom2)
      z  =  0.5_rk * s * xi * f_o_v0
      zr = z * f_o_v0 * (u - y*v)
      zi = z * (v + u * f_o_v0 * (w/(x1*nu) + y))
      sigma_m = sigma_m + cmplx(zr, zi)

      if ( n /= -1 ) then

! m = -(n-1)

        m = n + 1

! pi transition

        nu_offst = x1*(del_nu + 2.803_rk*m*(2-n)*h / denom1)
        call simple_voigt ( nu_offst, w, u, v )
        xi = 3.0_rk * (n * n - m * m) / denom2
        z  = 0.5_rk * s * xi * f_o_v0
        zr = z * f_o_v0 * (u - y*v)
        zi = z * (v + u * f_o_v0 * (w/(x1*nu) + y))
        pi = pi + cmplx(zr, zi)

! sigma_p transition

        nu_offst = x1*(del_nu + 2.803_rk*((m+1)*(2-n)-1)*h / denom1)
        call simple_voigt ( nu_offst, w, u, v )
        xi = 3.0_rk * (n + m) * (n + m + 1) / (4.0_rk * denom2)
        z  = 0.5_rk * s * xi * f_o_v0
        zr = z * f_o_v0 * (u - y*v)
        zi = z * (v + u * f_o_v0 * (w/(x1*nu) + y))
        sigma_p = sigma_p + cmplx(zr, zi)

! m = n - 1

        m = -(n + 1)

! pi transition

        nu_offst = x1*(del_nu + 2.803_rk*m*(2-n)*h / denom1)
        call simple_voigt ( nu_offst, w, u, v )
        xi = 3.0_rk * (n * n - m * m) / denom2
        z  = 0.5_rk * s * xi * f_o_v0
        zr = z * f_o_v0 * (u - y*v)
        zi = z * (v + u * f_o_v0 * (w/(x1*nu) + y))
        pi = pi + cmplx(zr, zi)

! sigma_m transition

        nu_offst = x1*(del_nu + 2.803_rk*((m-1)*(2-n)+1)*h / denom1)
        call simple_voigt ( nu_offst, w, u, v )
        xi = 3.0_rk * (m - n) * (m - n - 1) / (4.0_rk * denom2)
        z  = 0.5_rk * s * xi * f_o_v0
        zr = z * f_o_v0 * (u - y*v)
        zi = z * (v + u * f_o_v0 * (w/(x1*nu) + y))
        sigma_m = sigma_m + cmplx(zr, zi)

      else

        m = n + 1

! pi transition

        nu_offst = x1*(del_nu + 2.803_rk*m*(2-n)*h / denom1)
        call simple_voigt ( nu_offst, w, u, v )
        xi = 3.0_rk * (n * n - m * m) / denom2
        z  = 0.5_rk * s * xi * f_o_v0
        zr = z * f_o_v0 * (u - y*v)
        zi = z * (v + u * f_o_v0 * (w/(x1*nu) + y))
        pi = pi + cmplx(zr, zi)

      end if

    end if

  end subroutine Mag_O2_Abs_CS

! ------------------------------------------------------  Find_O2  -----
  subroutine Find_O2
  ! Find the O2 in the spectroscopy catalog

    use Molecules, only: l_o2
    use SpectroscopyCatalog_m, only: Catalog

    do o2_in_catalog = 1, size(catalog)
      if ( catalog(o2_in_catalog)%molecule == l_o2 ) return
    end do

    o2_in_catalog = -1 ! Not found
  end subroutine Find_O2

! ------------------------------------------  Get_QN_By_Frequency  -----
  subroutine Get_QN_By_Frequency ( V0, N )

    ! Get the quantum number for the O2 line that has a center frequency
    ! nearest to V0.  The quantum number of interest is given by the
    ! difference between the fourth and second numbers in the QN field
    ! of the line specification of the spectroscopy database.

    use MLSCommon, only: IP, R8, Rk => RP
    use SpectroscopyCatalog_M, only: Catalog, Lines

    real(r8), intent(in) :: V0
    integer(ip), intent(out) :: N

    integer(ip) :: I, J, K
    real(rk) :: Q, R

    ! Find O2 line the spectroscopy catalog
    if ( o2_in_catalog < 0 ) call find_o2
    n = -1
    ! Find the O2 line that has a center frequency nearest to V0
    j = catalog(o2_in_catalog)%lines(1)
    r = abs(v0-lines(j)%v0)
    do i = 2, size(catalog(o2_in_catalog)%lines)
      k = catalog(o2_in_catalog)%lines(i)
      q = abs(v0-lines(k)%v0)
      if ( q < r ) then
        r = q
        j = k
      end if
    end do
    if ( associated(lines(j)%qn) ) then
      if ( size(lines(j)%qn) >= 4 ) n = lines(j)%qn(4) - lines(j)%qn(2)
    end if
  end subroutine Get_QN_By_Frequency


! ------------------------------------------------  not_used_here  -----
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module O2_Abs_CS_M

! $Log$
! Revision 2.5  2003/05/17 00:30:39  pwagner
! Ousted last bogus ref to sp_o2
!
! Revision 2.4  2003/05/16 23:52:53  livesey
! Now uses molecule indices rather than spectags
!
! Revision 2.3  2003/05/16 02:45:08  vsnyder
! Removed USE's for unreferenced symbols
!
! Revision 2.2  2003/05/05 23:00:26  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.1.2.5  2003/03/05 03:29:59  vsnyder
! Add in the wing calculation
!
! Revision 2.1.2.4  2003/03/01 03:15:19  vsnyder
! Finish deleting Get_QN -- form the PUBLIC statement
!
! Revision 2.1.2.3  2003/03/01 03:13:47  vsnyder
! Delete unused procedure Get_QN so we won't need o2_DBase
!
! Revision 2.1.2.2  2003/03/01 03:11:02  vsnyder
! Use 'polarized' array for size, delete nonresonant computation
!
! Revision 2.1.2.1  2003/02/27 23:35:24  vsnyder
! Process all Zeeman-split lines, and no others
!
! Revision 2.1  2003/02/03 22:55:26  vsnyder
! Initial commit
!
