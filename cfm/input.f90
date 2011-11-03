! Copyright 2010, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.
module Input
   use MLSCommon, only: r8

   implicit none

!Ming: Haley is working on startTime, endTime 
   character(len=27) :: startTime = "2009-040T00:00:00.0000"
   character(len=27) :: endTime = "2009-040T00:00:00.9999"

   ! Replace startTime and endTime with L1 maf number in L1BOA
   integer, parameter :: startL1Maf = 7
   integer, parameter :: endL1Maf = 7

!Ming: change the example day to 2009 day040
   character(len=100) :: l1boa = "/data/emls/l1b/v02.23/2009/040/MLS-Aura_L1BOA_v02-23-c01_2009d040.h5"
   character(len=100) :: l1brad = "/data/emls/l1b/v02.23/2009/040/MLS-Aura_L1BRADG_v02-23-c01_2009d040.h5"
!Ming: I assume we keep these
   character(len=100) :: spectroscopy = '/data/emls/l2cal/MLS-Aura_L2Cal-Spectroscopy-PFA_v3-0-5i_0000d000.h5'
   character(len=100) :: leapsecFile = '/software/toolkit/LF6.1/toolkit/database/common/TD/leapsec.dat'
   character(len=100) :: antennaPatterns = '/data/emls/l2cal/MLS-Aura_L2Cal-AAAP_v2-0-0_0000d000.txt'
   character(len=100) :: filterShapes = '/data/emls/l2cal/MLS-Aura_L2Cal-Filters_v3-0-2_0000d000.txt'
   character(len=100) :: DACSFilterShapes = '/data/emls/l2cal/MLS-Aura_L2Cal-DACSFilters_v1-5-1_0000d000.txt'
   character(len=100) :: pointingGrids = '/data/emls/l2cal/MLS-Aura_L2Cal-PFG_v3-0-0_0000d000.txt'
   character(len=100), dimension(8) :: &
   pfaFiles= (/ &
   '/data/emls/l2cal/PFA_R5V_ME-03.h5                                                                   ', &
   '/data/emls/l2cal/PFA_R5H_ME-03.h5                                                                   ', &
   '/data/emls/l2cal/PFA_R4_ME-03.h5                                                                    ', &
   '/data/emls/l2cal/PFA_R3_ME-03.h5                                                                    ', &
   '/data/emls/l2cal/PFA_R2_ME-03.h5                                                                    ', &
   '/data/emls/l2cal/PFA_R1B_ME-03.h5                                                                   ', &
   '/data/emls/l2cal/PFA_R1A_ME-03.h5                                                                   ', &
   '/data/emls/l2cal/PFA_DACS_ME-03.h5                                                                  ' /)
   character(len=100), dimension(1) :: &
   l2pc = (/ &
   '/data/emls/l2cal/l2pc_30H6/MLS-Aura_L2Cal-L2PC-band9-LATSCALARHIRESO3HR_v3-00-HO-06_m02.h5          '/)
!Ming: I changed the following profiles to consist Temperature, ExtinctionV2, O3, 3_V1_3, O3_V2 at 55 levels;
!      HNO3, S_32_O2 at 37 levels; CO at <67 levels.  Note: is O_18_O set somewhere else? Extingction? O3_V1_3 etc?
   real(r8), dimension(55) :: &
   TemperatureInput= (/       297.158_r8,       289.724_r8,       284.356_r8,       276.756_r8,       265.709_r8, &
                254.132_r8,       249.935_r8,       236.775_r8,       225.171_r8,       214.553_r8, &
                206.161_r8,       199.582_r8,       194.390_r8,       191.274_r8,       192.612_r8, &
                197.358_r8,       203.675_r8,       209.069_r8,       212.745_r8,       215.404_r8, &
                218.021_r8,       221.393_r8,       225.156_r8,       228.391_r8,       230.946_r8, &
                233.755_r8,       237.216_r8,       241.305_r8,       245.444_r8,       249.706_r8, &
                254.621_r8,       259.924_r8,       264.561_r8,       267.950_r8,       269.553_r8, &
                268.795_r8,       264.352_r8,       257.105_r8,       248.433_r8,       242.248_r8, &
                235.518_r8,       226.208_r8,       221.316_r8,       212.292_r8,       201.969_r8, &
                192.598_r8,       178.324_r8,       172.662_r8,       169.580_r8,       172.959_r8, &
                185.157_r8,       212.902_r8,       284.879_r8,       382.289_r8,       382.289_r8 /)
   real(r8), dimension(55) :: &
   H2OInput        = (/    0.00241224_r8,    0.00172104_r8,    0.00139789_r8,   0.000935927_r8,   0.000450439_r8, &
                0.000190159_r8,   8.86722e-05_r8,   3.51943e-06_r8,   0.000137363_r8,   5.49711e-05_r8, &
                2.05114e-05_r8,   2.01636e-07_r8,   3.24751e-06_r8,   4.51751e-06_r8,   3.87043e-06_r8, &
                3.14470e-06_r8,   3.48467e-06_r8,   3.93302e-06_r8,   3.98300e-06_r8,   4.19780e-06_r8, &
                4.75373e-06_r8,   5.03954e-06_r8,   4.90524e-06_r8,   4.78210e-06_r8,   4.91733e-06_r8, &
                5.10569e-06_r8,   5.20030e-06_r8,   5.26312e-06_r8,   5.53712e-06_r8,   5.87348e-06_r8, &
                5.91152e-06_r8,   5.69914e-06_r8,   5.67516e-06_r8,   6.01006e-06_r8,   6.38157e-06_r8, &
                6.61089e-06_r8,   6.13885e-06_r8,   7.04676e-06_r8,   8.08649e-06_r8,   6.43817e-06_r8, &
                6.35526e-06_r8,   6.83005e-06_r8,   7.75657e-06_r8,   7.30504e-06_r8,   5.95156e-06_r8, &
                5.06159e-06_r8,   2.75447e-06_r8,   7.92934e-07_r8,   1.28573e-07_r8,   5.27341e-07_r8, &
                5.27341e-07_r8,   5.27341e-07_r8,   5.27341e-07_r8,   5.27341e-07_r8,   5.27341e-07_r8 /)
   real(r8), dimension(55) :: &
   O3Input         = (/  1.25014e-08_r8,   2.18536e-08_r8,   3.12059e-08_r8,   3.77980e-08_r8,   4.43901e-08_r8, &
                4.81388e-08_r8,   7.53838e-08_r8,  -4.84879e-08_r8,   3.47201e-08_r8,   2.25237e-07_r8, &
                2.32295e-07_r8,  -5.50474e-09_r8,   3.22532e-08_r8,   3.27599e-07_r8,   5.12197e-07_r8, &
                1.07709e-06_r8,   1.64005e-06_r8,   2.69266e-06_r8,   3.79579e-06_r8,   5.10495e-06_r8, &
                6.27606e-06_r8,   7.29860e-06_r8,   8.48605e-06_r8,   8.85052e-06_r8,   9.08717e-06_r8, &
                8.92521e-06_r8,   8.94099e-06_r8,   8.51976e-06_r8,   8.06099e-06_r8,   7.45707e-06_r8, &
                6.75821e-06_r8,   5.83119e-06_r8,   4.78850e-06_r8,   4.03797e-06_r8,   3.37683e-06_r8, &
                3.00758e-06_r8,   2.80204e-06_r8,   2.37477e-06_r8,   1.83857e-06_r8,   1.47623e-06_r8, &
                1.02984e-06_r8,   4.47357e-07_r8,   5.67443e-07_r8,  -3.06701e-08_r8,  -1.79143e-06_r8, &
               -1.57880e-06_r8,   2.02320e-06_r8,   2.33698e-06_r8,   2.86198e-06_r8,   1.94295e-06_r8, &
                1.94295e-06_r8,   1.94295e-06_r8,   1.94295e-06_r8,   1.94295e-06_r8,   1.94295e-06_r8 /)
   real(r8), dimension(37) :: &
   SO2Input         = (/0.00000_r8,       0.00000_r8,       0.00000_r8,   3.32686e-08_r8,  -2.60494e-08_r8, &
               -5.61417e-09_r8,   2.56858e-09_r8,   1.93806e-09_r8,  -1.18218e-09_r8,  -4.93486e-09_r8, &
                1.64919e-08_r8,  -9.11367e-10_r8,  -4.31667e-10_r8,  -4.67996e-09_r8,   2.79592e-09_r8, &
                1.19575e-08_r8,   5.17797e-09_r8,   2.99216e-09_r8,   6.40992e-09_r8,   2.37867e-09_r8, &
               -9.46628e-10_r8,  -1.31347e-09_r8,   6.71434e-11_r8,   6.60937e-10_r8,  -4.07461e-11_r8, &
                5.35449e-11_r8,  -3.21519e-10_r8,  -3.27190e-10_r8,  -4.46835e-11_r8,   3.85534e-11_r8, &
               -2.90539e-12_r8,       0.00000_r8,       0.00000_r8,       0.00000_r8,       0.00000_r8, &
                0.00000_r8,       0.00000_r8 /)

   real(r8), dimension(37) :: &
   HNO3Input         = (/ 1.09330e-11_r8,   4.07050e-11_r8,   6.55338e-11_r8,  -1.02235e-09_r8,  -4.70137e-10_r8, &
                4.37410e-10_r8,  -1.24584e-10_r8,   8.45327e-10_r8,   4.55798e-09_r8,   4.16154e-09_r8, &
                4.20141e-09_r8,   3.71814e-09_r8,   2.73836e-09_r8,   2.15569e-09_r8,   1.99514e-10_r8, &
                1.70946e-11_r8,  -1.84183e-09_r8,  -8.75599e-10_r8,   1.65315e-09_r8,   2.08126e-09_r8, &
               -1.37495e-09_r8,   8.02555e-10_r8,  -4.83877e-10_r8,  -7.17038e-09_r8,  -1.03526e-08_r8, &
               -9.59712e-09_r8,  -4.17081e-09_r8,  -1.30650e-09_r8,   2.10400e-15_r8,   1.43523e-17_r8, &
                1.93479e-19_r8,   1.93479e-19_r8,   1.93479e-19_r8,   1.93479e-19_r8,   1.93479e-19_r8, &
                1.93479e-19_r8,   1.93479e-19_r8 /)
   real(r8), dimension(65) :: &
   COInput         = (/ 6.80256e-08_r8,    6.79563e-08_r8,    6.74404e-08_r8,    6.69285e-08_r8,    6.64205e-08_r8, &
                6.59163e-08_r8,    6.48891e-08_r8,    6.38780e-08_r8,    6.28826e-08_r8,    6.30514e-08_r8, &
                6.32205e-08_r8,    6.33901e-08_r8,    6.34645e-08_r8,    6.35390e-08_r8,    6.36135e-08_r8, &
                6.32930e-08_r8,    6.29741e-08_r8,    6.26568e-08_r8,    6.16546e-08_r8,    6.06685e-08_r8, &
                5.96981e-08_r8,    5.77026e-08_r8,    5.57741e-08_r8,    5.39099e-08_r8,    5.26075e-08_r8, &
                5.13365e-08_r8,    5.00962e-08_r8,    4.41271e-08_r8,    3.88691e-08_r8,    3.42377e-08_r8, &
                3.01577e-08_r8,    2.65643e-08_r8,    2.33990e-08_r8,    2.17289e-08_r8,    2.01780e-08_r8, &
                1.87377e-08_r8,    1.74003e-08_r8,    1.61583e-08_r8,    1.50050e-08_r8,    1.56685e-08_r8, &
                1.63615e-08_r8,    1.70851e-08_r8,    1.78407e-08_r8,    1.86296e-08_r8,    1.94535e-08_r8, &
                2.03139e-08_r8,    2.12122e-08_r8,    2.21503e-08_r8,    2.31298e-08_r8,    2.41528e-08_r8, &
                2.52210e-08_r8,    2.63364e-08_r8,    2.87173e-08_r8,    3.26982e-08_r8,    3.41443e-08_r8, &
                3.42490e-08_r8,    3.43014e-08_r8,    3.43540e-08_r8,    3.44330e-08_r8,    3.44857e-08_r8, &
                3.45650e-08_r8,    3.46710e-08_r8,    3.48306e-08_r8,    3.49909e-08_r8,    3.52058e-08_r8 /)

    real(r8), dimension(47) :: &
    extinctionV2R3Input = (/ 1.10044_r8,      0.481280_r8,     0.158084_r8,     5.825650E-02_r8, 2.094034E-02_r8, &
                        1.018303E-02_r8, 5.922833E-03_r8, 2.936215E-03_r8, 1.437336E-03_r8, 7.931669E-04_r8, &
                        4.991122E-04_r8, 3.534503E-04_r8, 2.518090E-04_r8, 1.649636E-04_r8, 1.088406E-04_r8, &
                        6.841273E-05_r8, 4.305684E-05_r8, 2.769214E-05_r8, 1.783691E-05_r8, 1.167463E-05_r8, &
                        7.642561E-06_r8, 3.310648E-06_r8, 1.435593E-06_r8, 6.105654E-07_r8, 2.500542E-07_r8, &
                        1.025592E-07_r8, 4.297000E-08_r8, 1.877126E-08_r8, 8.611204E-09_r8, 4.155832E-09_r8, &
                        2.095296E-09_r8, 1.083296E-09_r8, 5.644514E-10_r8, 2.880354E-10_r8, 1.470881E-10_r8, &
                        3.675219E-11_r8, 8.514243E-12_r8, 1.923638E-12_r8, 4.633676E-13_r8, 1.155358E-13_r8, &
                        2.669294E-14_r8, 5.691287E-15_r8, 1.124997E-15_r8, 1.468009E-16_r8, 1.615710E-17_r8, &
                        2.040401E-18_r8, 1.493425E-19_r8 /)


!Ming: ask Bill/Paul what is this and if we need this:
   real(r8) :: refGPHInput = 16000.0_r8
end module
