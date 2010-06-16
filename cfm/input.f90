module Input
   use MLSCommon, only: r8

   implicit none

   character(len=27) :: startTime = "2005-037T00:00:00.0000"
   character(len=27) :: endTime = "2005-037T00:00:00.9999"
   character(len=100) :: l1boa = "/data/emls/l1b/v02.23/2005/037/MLS-Aura_L1BOA_v02-23-c01_2005d037.h5"
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
end module
