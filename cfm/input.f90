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

   character(len=27) :: startTime = "2005-037T00:00:00.0000"
   character(len=27) :: endTime = "2005-037T00:00:00.9999"
   character(len=100) :: l1boa = "/data/emls/l1b/v02.23/2005/037/MLS-Aura_L1BOA_v02-23-c01_2005d037.h5"
   character(len=100) :: l1brad = "/data/emls/l1b/v02.23/2005/037/MLS-Aura_L1BRADG_v02-23-c01_2005d037.h5"
   character(len=100) :: spectroscopy = '/data/emls/l2cal/MLS-Aura_L2Cal-Spectroscopy-PFA_v3-0-4_0000d000.h5'
   character(len=100) :: leapsecFile = '/software/toolkit/LF6.1/toolkit/database/common/TD/leapsec.dat'
   character(len=100) :: antennaPatterns = '/data/emls/l2cal/MLS-Aura_L2Cal-AAAP_v2-0-0_0000d000.txt'
   character(len=100) :: filterShapes = '/data/emls/l2cal/MLS-Aura_L2Cal-Filters_v3-0-2_0000d000.txt'
   character(len=100) :: DACSFilterShapes = '/data/emls/l2cal/MLS-Aura_L2Cal-DACSFilters_v1-5-1_0000d000.txt'
   character(len=100) :: pointingGrids = '/data/emls/l2cal/MLS-Aura_L2Cal-PFG_v3-0-0_0000d000.txt'
   character(len=100), dimension(8) :: &
   pfaFiles= (/ &
   '/data/emls/l2cal/PFA_R5V_FS-04.h5                                                                   ', &
   '/data/emls/l2cal/PFA_R5H_FS-04.h5                                                                   ', &
   '/data/emls/l2cal/PFA_R4_FS-04.h5                                                                    ', &
   '/data/emls/l2cal/PFA_R3_FS-04.h5                                                                    ', &
   '/data/emls/l2cal/PFA_R2_FS-04.h5                                                                    ', &
   '/data/emls/l2cal/PFA_R1B_FS-04.h5                                                                   ', &
   '/data/emls/l2cal/PFA_R1A_FS-04.h5                                                                   ', &
   '/data/emls/l2cal/PFA_DACS_FS-04.h5                                                                  ' /)
   character(len=100), dimension(3) :: &
   l2pc = (/ &
   '/data/emls/l2cal/l2pc_30H6/MLS-Aura_L2Cal-L2PC-band7-LATSCALARHIRESO3HR_v3-00-HO-06_m02.h5          ', &
   '/data/emls/l2cal/l2pc_30H6/MLS-Aura_L2Cal-L2PC-band2-LATSCALARHIRES_v3-00-HO-06_m02.h5              ', &
   '/data/emls/l2cal/l2pc_30H6/MLS-Aura_L2Cal-L2PC-band8-LATSCALARHIRESO3HR_v3-00-HO-06_m02.h5          '/)
   real(r8), dimension(37) :: &
   TemperatureInput= (/ 2.9128e+02_r8,  2.9149e+02_r8,  2.8439e+02_r8,  2.7649e+02_r8,  2.6942e+02_r8, &
                2.5727e+02_r8,  2.4453e+02_r8,  2.2941e+02_r8,  2.1723e+02_r8,  2.1070e+02_r8, &
                2.0637e+02_r8,  2.0195e+02_r8,  1.9559e+02_r8,  1.8838e+02_r8,  1.9129e+02_r8, &
                2.0175e+02_r8,  2.1001e+02_r8,  2.1052e+02_r8,  2.1138e+02_r8,  2.1788e+02_r8, &
                2.2112e+02_r8,  2.2371e+02_r8,  2.2438e+02_r8,  2.3447e+02_r8,  2.4688e+02_r8, &
                2.4973e+02_r8,  2.5171e+02_r8,  2.6056e+02_r8,  2.7136e+02_r8,  2.6748e+02_r8, &
                1.7861e+02_r8,  1.6997e+02_r8,  1.8009e+02_r8,  2.0964e+02_r8,  2.8234e+02_r8, &
                3.8117e+02_r8,  3.8117e+02_r8 /)
   real(r8), dimension(37) :: &
   H2OInput        = (/ 2.3632e-02_r8,  2.9074e-02_r8,  2.0823e-02_r8,  1.3613e-02_r8,  9.2022e-03_r8, &
                3.7962e-03_r8,  3.1691e-03_r8,  3.1161e-04_r8,  3.4557e-05_r8,  1.8799e-05_r8, &
                9.2256e-06_r8,  5.4310e-06_r8,  4.3650e-06_r8,  3.5293e-06_r8,  3.4100e-06_r8, &
                3.5061e-06_r8,  4.0785e-06_r8,  4.3385e-06_r8,  2.7219e-06_r8,  5.1609e-06_r8, &
                4.3181e-06_r8,  4.9063e-06_r8,  4.7031e-06_r8,  5.3227e-06_r8,  5.4470e-06_r8, &
                5.8672e-06_r8,  5.3741e-06_r8,  5.9369e-06_r8,  6.7043e-06_r8,  6.9746e-06_r8, &
                6.9707e-06_r8,  6.4823e-06_r8,  6.3464e-06_r8,  6.7609e-06_r8,  6.9218e-06_r8, &
                5.3200e-07_r8,  5.3200e-07_r8 /)
   real(r8), dimension(37) :: &
   O3Input         = (/ 1.3001e-08_r8,  3.2795e-08_r8,  4.4579e-08_r8,  4.0661e-07_r8,  8.0144e-08_r8, &
                9.1738e-08_r8,  1.2884e-07_r8,  3.1488e-07_r8,  1.9350e-06_r8,  3.7353e-06_r8, &
                6.8007e-06_r8,  8.5696e-06_r8,  9.0027e-06_r8,  9.2723e-06_r8,  8.1097e-06_r8, &
                6.7720e-06_r8,  5.1827e-06_r8,  3.8637e-06_r8,  3.2787e-06_r8,  2.5215e-06_r8, &
                1.8182e-06_r8,  1.4148e-06_r8,  1.2598e-06_r8,  1.3303e-06_r8,  1.1685e-06_r8, &
                1.2587e-06_r8,  4.0577e-07_r8,  1.1117e-06_r8,  2.3770e-06_r8,  1.4602e-06_r8, &
                1.6944e-06_r8,  1.9545e-06_r8,  1.9545e-06_r8,  1.9545e-06_r8,  1.9545e-06_r8, &
                1.9545e-06_r8,  1.9545e-06_r8 /)
   real(r8), dimension(37) :: &
   TemperatureInput2= (/ 2.7927e+02_r8,  2.7695e+02_r8,  2.7325e+02_r8,  2.6708e+02_r8,  2.5718e+02_r8, &
                2.4276e+02_r8,  2.2883e+02_r8,  2.1503e+02_r8,  2.0760e+02_r8,  2.1141e+02_r8, &
                2.2078e+02_r8,  2.2327e+02_r8,  2.1580e+02_r8,  2.0642e+02_r8,  2.0583e+02_r8, &
                2.1095e+02_r8,  2.1502e+02_r8,  2.1524e+02_r8,  2.1338e+02_r8,  2.1086e+02_r8, &
                2.1983e+02_r8,  2.1632e+02_r8,  2.1611e+02_r8,  2.2542e+02_r8,  2.2380e+02_r8, &
                2.2027e+02_r8,  2.3052e+02_r8,  2.3974e+02_r8,  2.4960e+02_r8,  2.4963e+02_r8, &
                2.4025e+02_r8,  2.3605e+02_r8,  2.4090e+02_r8,  2.4880e+02_r8,  2.3582e+02_r8, &
                3.5728e+02_r8,  3.5728e+02_r8 /)
   real(r8), dimension(37) :: &
   H2OInput2        = (/ 9.1453e-03_r8,  9.2172e-03_r8,  8.2678e-03_r8,  5.9595e-03_r8,  2.9752e-03_r8, &
                8.7081e-04_r8,  4.3487e-05_r8,  1.7336e-05_r8,  6.6668e-06_r8,  4.3221e-06_r8, &
                3.6797e-06_r8,  4.7577e-06_r8,  4.8634e-06_r8,  4.0305e-06_r8,  3.8945e-06_r8, &
                4.1516e-06_r8,  4.7683e-06_r8,  5.2437e-06_r8,  3.4030e-06_r8,  6.0271e-06_r8, &
                5.4048e-06_r8,  5.9748e-06_r8,  7.4699e-06_r8,  7.9878e-06_r8,  6.5561e-06_r8, &
                4.8183e-06_r8,  3.3121e-06_r8,  2.8018e-06_r8,  2.7419e-06_r8,  3.2127e-06_r8, &
                1.0000e-07_r8,  3.5347e-07_r8,  3.5347e-07_r8,  3.5347e-07_r8,  3.5347e-07_r8, &
                3.5347e-07_r8,  3.5347e-07_r8 /)
   real(r8), dimension(37) :: &
   O3Input2         = (/ 4.9074e-08_r8,  5.6807e-08_r8,  7.0802e-08_r8,  1.0482e-07_r8,  1.5493e-07_r8, &
                4.2987e-07_r8,  4.3618e-07_r8,  1.2536e-06_r8,  3.2734e-06_r8,  5.3652e-06_r8, &
                5.0849e-06_r8,  5.7286e-06_r8,  5.6696e-06_r8,  6.6899e-06_r8,  6.4951e-06_r8, &
                7.4254e-06_r8,  7.7712e-06_r8,  5.4279e-06_r8,  4.3685e-06_r8,  2.2732e-06_r8, &
                1.3391e-06_r8,  1.0155e-06_r8,  9.3310e-07_r8,  8.8999e-07_r8,  3.8050e-07_r8, &
                2.3795e-07_r8, -9.1263e-07_r8,  9.1669e-07_r8,  1.3863e-06_r8, -1.8748e-06_r8, &
                1.7927e-06_r8,  2.9872e-06_r8,  2.9872e-06_r8,  2.9872e-06_r8,  2.9872e-06_r8, &
                2.9872e-06_r8,  2.9872e-06_r8 /)

   real(r8) :: refGPHInput = 16000.0_r8
end module
