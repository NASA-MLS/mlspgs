! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module O2_DBase

! Module creates the frequency independent oxygen data base for
! cross section calculations.

!**********************  Revision History  *****************************
!                                                                      *
!     early  1991  W. G. Read     Program created with initial data    *
!                                 based on Liebe 1989 MPM              *
!     17 Jun 1991  Z.    Shippony Increased broadening parameter by 5% *
!                                 based on phone conversation with     *
!                                 H. Liebe                             *
!     21 Feb 1992  W. G. Read     Widths and interferences are modified*
!                                 to agree with Liebe 1991 MPM         *
!                                 Removed 5% broadening parameter      *
!                                 increase                             *
!      5 Mar 1992  W. G. Read     Increase data base to include more   *
!                                 lines and more flexible temperature  *
!                                 dependences per line                 *
!      2 Dec 1994  Z.    Shippony Added 10 more lines (34 to 44)       *
!                                                                      *
!***********************************************************************

  use MLSCommon, only: Rk => R8

  implicit NONE
  public

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    &  "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  Integer i, j

  integer, parameter:: maxlines = 44

! Load in the data base

  real(rk) :: o2_freq(maxlines) = &
    & (/ 118750.341e0_rk, 56264.777e0_rk, 62486.253e0_rk, 58446.589e0_rk, &
    &     60306.057e0_rk, 59590.982e0_rk, 59164.204e0_rk, 60434.775e0_rk, &
    &     58323.874e0_rk, 61150.558e0_rk, 57612.481e0_rk, 61800.152e0_rk, &
    &     56968.180e0_rk, 62411.212e0_rk, 56363.387e0_rk, 62997.971e0_rk, &
    &     55783.800e0_rk, 63568.520e0_rk, 55221.365e0_rk, 64127.764e0_rk, &
    &     54671.157e0_rk, 64678.900e0_rk, 54129.999e0_rk, 65224.067e0_rk, &
    &     53595.748e0_rk, 65764.769e0_rk, 53066.906e0_rk, 66302.088e0_rk, &
    &     52542.303e0_rk, 66836.827e0_rk, 52021.409e0_rk, 67369.595e0_rk, &
    &     51503.350e0_rk, 67900.862e0_rk, 50987.749e0_rk, 68431.005e0_rk, &
    &     50474.238e0_rk, 68960.311e0_rk,368498.350e0_rk,424763.124e0_rk, &
    &    487249.370e0_rk,715393.150e0_rk,773839.675e0_rk,834145.330e0_rk /)

  real(rk), private, parameter :: O2_Data_1_17(7,17) = reshape( (/ &
! N Quantum  G.S.Energy Int Intensity Width Interference Coefficients TD
!      Number    cm^{-1}   nm^{2} MHz MHz/mBar   1/mBar     1/mBar t_dep_w
    & -1.0e0_rk,   0.2458e0_rk,-6.5323e0_rk,1.630e0_rk, -0.310e-4_rk, 0.080e-4_rk, 0.8e0_rk, &
    &  1.0e0_rk,   2.3301e0_rk,-7.0948e0_rk,1.646e0_rk,  3.390e-4_rk,-0.980e-4_rk, 0.8e0_rk, &
    & -3.0e0_rk,  16.4987e0_rk,-6.6075e0_rk,1.468e0_rk, -4.330e-4_rk, 0.840e-4_rk, 0.8e0_rk, &
    &  3.0e0_rk,  16.6334e0_rk,-6.6544e0_rk,1.449e0_rk,  6.500e-4_rk,-1.270e-4_rk, 0.8e0_rk, &
    & -5.0e0_rk,  42.4460e0_rk,-6.4769e0_rk,1.382e0_rk, -6.130e-4_rk, 0.700e-4_rk, 0.8e0_rk, &
    &  5.0e0_rk,  42.4698e0_rk,-6.4846e0_rk,1.360e0_rk,  6.650e-4_rk,-0.780e-4_rk, 0.8e0_rk, &
    & -7.0e0_rk,  79.8529e0_rk,-6.4315e0_rk,1.319e0_rk, -6.280e-4_rk, 2.310e-4_rk, 0.8e0_rk, &
    &  7.0e0_rk,  79.8105e0_rk,-6.4120e0_rk,1.297e0_rk,  6.060e-4_rk,-2.820e-4_rk, 0.8e0_rk, &
    & -9.0e0_rk, 128.7379e0_rk,-6.4410e0_rk,1.266e0_rk, -1.780e-4_rk, 0.440e-4_rk, 0.8e0_rk, &
    &  9.0e0_rk, 128.6437e0_rk,-6.3993e0_rk,1.248e0_rk,  0.900e-4_rk,-0.580e-4_rk, 0.8e0_rk, &
    &-11.0e0_rk, 189.0990e0_rk,-6.4933e0_rk,1.221e0_rk, -5.330e-4_rk, 6.060e-4_rk, 0.8e0_rk, &
    & 11.0e0_rk, 188.9594e0_rk,-6.4319e0_rk,1.207e0_rk,  4.960e-4_rk,-6.620e-4_rk, 0.8e0_rk, &
    &-13.0e0_rk, 260.9285e0_rk,-6.5824e0_rk,1.181e0_rk, -3.620e-4_rk, 6.450e-4_rk, 0.8e0_rk, &
    & 13.0e0_rk, 260.7469e0_rk,-6.5028e0_rk,1.171e0_rk,  3.130e-4_rk,-6.760e-4_rk, 0.8e0_rk, &
    &-15.0e0_rk, 344.2156e0_rk,-6.7047e0_rk,1.144e0_rk, -2.580e-4_rk, 6.550e-4_rk, 0.8e0_rk, &
    & 15.0e0_rk, 343.7484e0_rk,-6.6076e0_rk,1.139e0_rk,  2.080e-4_rk,-6.680e-4_rk, 0.8e0_rk, &
    &-17.0e0_rk, 438.9474e0_rk,-6.8580e0_rk,1.110e0_rk, -1.440e-4_rk, 6.130e-4_rk, 0.8e0_rk /), &
    & (/ 7, 17 /) )

  real(rk), private, parameter :: O2_Data_18_34(7,17) = reshape( (/ &
    & 17.0e0_rk, 438.4418e0_rk,-6.7441e0_rk,1.110e0_rk,  0.940e-4_rk,-6.140e-4_rk, 0.8e0_rk, &
    &-19.0e0_rk, 545.1087e0_rk,-7.0408e0_rk,1.080e0_rk,  2.240e-4_rk, 2.950e-4_rk, 0.8e0_rk, &
    & 19.0e0_rk, 544.8117e0_rk,-6.9105e0_rk,1.080e0_rk, -2.700e-4_rk,-2.890e-4_rk, 0.8e0_rk, &
    &-21.0e0_rk, 662.6827e0_rk,-7.2518e0_rk,1.050e0_rk,  3.250e-4_rk, 2.650e-4_rk, 0.8e0_rk, &
    & 21.0e0_rk, 662.3489e0_rk,-7.1054e0_rk,1.050e0_rk, -3.660e-4_rk,-2.590e-4_rk, 0.8e0_rk, &
    &-23.0e0_rk, 791.6504e0_rk,-7.4903e0_rk,1.020e0_rk,  2.910e-4_rk, 3.750e-4_rk, 0.8e0_rk, &
    & 23.0e0_rk, 791.2803e0_rk,-7.3279e0_rk,1.020e0_rk, -3.260e-4_rk,-3.680e-4_rk, 0.8e0_rk, &
    &-25.0e0_rk, 931.9908e0_rk,-7.7557e0_rk,1.000e0_rk,  2.000e-4_rk, 5.080e-4_rk, 0.8e0_rk, &
    & 25.0e0_rk, 931.5849e0_rk,-7.5775e0_rk,1.000e0_rk, -2.320e-4_rk,-5.000e-4_rk, 0.8e0_rk, &
    &-27.0e0_rk,1083.6814e0_rk,-8.0474e0_rk,0.970e0_rk,  1.140e-4_rk, 6.210e-4_rk, 0.8e0_rk, &
    & 27.0e0_rk,1083.2400e0_rk,-7.8535e0_rk,0.970e0_rk, -1.460e-4_rk,-6.090e-4_rk, 0.8e0_rk, &
    &-29.0e0_rk,1246.6976e0_rk,-8.3650e0_rk,0.940e0_rk,  1.180e-4_rk, 6.530e-4_rk, 0.8e0_rk, &
    & 29.0e0_rk,1246.2208e0_rk,-8.1554e0_rk,0.940e0_rk, -1.470e-4_rk,-6.390e-4_rk, 0.8e0_rk, &
    &-31.0e0_rk,1421.0129e0_rk,-8.7081e0_rk,0.920e0_rk,  1.440e-4_rk, 6.640e-4_rk, 0.8e0_rk, &
    & 31.0e0_rk,1420.5009e0_rk,-8.4830e0_rk,0.920e0_rk, -1.740e-4_rk,-6.470e-4_rk, 0.8e0_rk, &
    &-33.0e0_rk,1606.5990e0_rk,-9.0766e0_rk,0.890e0_rk,  1.710e-4_rk, 6.730e-4_rk, 0.8e0_rk, &
    & 33.0e0_rk,1606.0520e0_rk,-8.8359e0_rk,0.890e0_rk, -1.980e-4_rk,-6.550e-4_rk, 0.8e0_rk /), &
    & (/ 7, 17 /) )

  real(rk), private, parameter :: O2_Data_35_44(7,10) = reshape( (/ &
    &-35.0e0_rk,1803.4400e0_rk,-9.4701e0_rk,0.870e0_rk,  1.900e-4_rk, 6.800e-4_rk, 0.8e0_rk, &
    & 35.0e0_rk,1802.8600e0_rk,-9.2139e0_rk,0.870e0_rk, -2.100e-4_rk,-6.600e-4_rk, 0.8e0_rk, &
    &-37.0e0_rk,2011.4800e0_rk,-9.8883e0_rk,0.850e0_rk,  2.100e-4_rk, 6.850e-4_rk, 0.8e0_rk, &
    & 37.0e0_rk,2010.8100e0_rk,-9.6166e0_rk,0.850e0_rk, -2.200e-4_rk,-6.650e-4_rk, 0.8e0_rk, &
    &  0.0e0_rk,   4.2070e0_rk,-7.1852e0_rk,1.920e0_rk,  0.000e+0_rk, 0.000e+0_rk, 0.2e0_rk, &
    &  0.0e0_rk,   2.3302e0_rk,-6.1476e0_rk,1.926e0_rk,  0.000e+0_rk, 0.000e+0_rk, 0.2e0_rk, &
    &  0.0e0_rk,   2.3302e0_rk,-6.5168e0_rk,1.920e0_rk,  0.000e+0_rk, 0.000e+0_rk, 0.2e0_rk, &
    &  0.0e0_rk,  18.5832e0_rk,-6.7340e0_rk,1.810e0_rk,  0.000e+0_rk, 0.000e+0_rk, 0.2e0_rk, &
    &  0.0e0_rk,  16.6336e0_rk,-5.9339e0_rk,1.810e0_rk,  0.000e+0_rk, 0.000e+0_rk, 0.2e0_rk, &
    &  0.0e0_rk,  16.6336e0_rk,-6.3965e0_rk,1.810e0_rk,  0.000e+0_rk, 0.000e+0_rk, 0.2e0_rk /), &
    & (/ 7, 10 /) )

  real(rk), parameter :: O2_Data(7,maxlines) = reshape( &
    & (/ o2_data_1_17, o2_data_18_34, o2_data_35_44 /), (/ 7, maxlines /) )

! o2 data for partition function

  real(rk), parameter :: o2_logq(3) = (/ 2.3398e0_rk, 2.215e0_rk, 2.0398e0_rk /)

! o2 data for temperature dependences

  real(rk), parameter :: temp_dep_d = 0.8e0_rk, temp_dep_g = 1.8e0_rk
  real(rk), parameter :: temp_dep_wnr = 0.8e0_rk

! o2 mass

  real(rk), parameter :: o2_mass = 32.0e0_rk

! Non resonant amplitude and broadening constants

  real(rk), parameter :: debye = 1.10816e-3_rk, w_nonres = 0.48e0_rk

contains

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module O2_DBase

! $Log$
! Revision 2.1  2003/02/03 23:08:09  vsnyder
! Initial commit
!
