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
   pfaFiles= (/'/data/emls/l2cal/PFA_R5V_FS-04.h5', &
    '/data/emls/l2cal/PFA_R5H_FS-04.h5', &
    '/data/emls/l2cal/PFA_R4_FS-04.h5 ', &
    '/data/emls/l2cal/PFA_R3_FS-04.h5 ', &
    '/data/emls/l2cal/PFA_R2_FS-04.h5 ', &
    '/data/emls/l2cal/PFA_R1B_FS-04.h5', &
    '/data/emls/l2cal/PFA_R1A_FS-04.h5', &
    '/data/emls/l2cal/PFA_DACS_FS-04.h5' /)
   character(len=100), dimension(1) :: &
   l2pc = (/'/data/emls/l2cal/l2pc_30H6/MLS-Aura_L2Cal-L2PC-band7-LATSCALARHIRESO3HR_v3-00-HO-06_m02.h5'/)
   real(r8), dimension(37) :: &
   TemperatureInput= (/ 2.9128e+02,  2.9149e+02,  2.8439e+02,  2.7649e+02,  2.6942e+02, &
                2.5727e+02,  2.4453e+02,  2.2941e+02,  2.1723e+02,  2.1070e+02, &
                2.0637e+02,  2.0195e+02,  1.9559e+02,  1.8838e+02,  1.9129e+02, &
                2.0175e+02,  2.1001e+02,  2.1052e+02,  2.1138e+02,  2.1788e+02, &
                2.2112e+02,  2.2371e+02,  2.2438e+02,  2.3447e+02,  2.4688e+02, &
                2.4973e+02,  2.5171e+02,  2.6056e+02,  2.7136e+02,  2.6748e+02, &
                1.7861e+02,  1.6997e+02,  1.8009e+02,  2.0964e+02,  2.8234e+02, &
                3.8117e+02,  3.8117e+02 /)
   real(r8), dimension(37) :: &
   GPHInput        = (/ 3.0347e+02,  1.9401e+03,  3.5574e+03,  5.1326e+03,  6.6657e+03, &
                8.1449e+03,  9.5541e+03,  1.0885e+04,  1.2139e+04,  1.3341e+04, &
                1.4513e+04,  1.5659e+04,  1.6776e+04,  1.7854e+04,  1.8920e+04, &
                2.0024e+04,  2.1181e+04,  2.2362e+04,  2.3546e+04,  2.4752e+04, &
                2.5985e+04,  2.8483e+04,  3.1000e+04,  3.3577e+04,  3.6281e+04, &
                3.9070e+04,  4.1887e+04,  4.4764e+04,  4.7752e+04,  5.0778e+04, &
                5.3718e+04,  5.6560e+04,  5.9348e+04,  6.2052e+04,  6.4649e+04, &
                1.1777e+05,  1.2727e+05 /)
   real(r8), dimension(37) :: &
   H2OInput        = (/ 2.3632e-02,  2.9074e-02,  2.0823e-02,  1.3613e-02,  9.2022e-03, &
                3.7962e-03,  3.1691e-03,  3.1161e-04,  3.4557e-05,  1.8799e-05, &
                9.2256e-06,  5.4310e-06,  4.3650e-06,  3.5293e-06,  3.4100e-06, &
                3.5061e-06,  4.0785e-06,  4.3385e-06,  2.7219e-06,  5.1609e-06, &
                4.3181e-06,  4.9063e-06,  4.7031e-06,  5.3227e-06,  5.4470e-06, &
                5.8672e-06,  5.3741e-06,  5.9369e-06,  6.7043e-06,  6.9746e-06, &
                6.9707e-06,  6.4823e-06,  6.3464e-06,  6.7609e-06,  6.9218e-06, &
                5.3200e-07,  5.3200e-07 /)
   real(r8), dimension(37) :: &
   O3Input         = (/ 1.3001e-08,  3.2795e-08,  4.4579e-08,  4.0661e-07,  8.0144e-08, &
                9.1738e-08,  1.2884e-07,  3.1488e-07,  1.9350e-06,  3.7353e-06, &
                6.8007e-06,  8.5696e-06,  9.0027e-06,  9.2723e-06,  8.1097e-06, &
                6.7720e-06,  5.1827e-06,  3.8637e-06,  3.2787e-06,  2.5215e-06, &
                1.8182e-06,  1.4148e-06,  1.2598e-06,  1.3303e-06,  1.1685e-06, &
                1.2587e-06,  4.0577e-07,  1.1117e-06,  2.3770e-06,  1.4602e-06, &
                1.6944e-06,  1.9545e-06,  1.9545e-06,  1.9545e-06,  1.9545e-06, &
                1.9545e-06,  1.9545e-06 /)
end module
