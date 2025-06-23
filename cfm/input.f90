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

   ! StartProfile and endProfile will be converted 
   ! to MAF number, Pranjit Saha
   integer, parameter :: startProfile =  1259
   integer, parameter :: endProfile =  1259
 
   character(len=100) :: l1boa = "/data/emls/l1b/v04.20/2007/207/MLS-Aura_L1BOA_v04-20-c01_2007d207.h5"
   character(len=100) :: l1brad = "/data/emls/l1b/v04.20/2007/207/MLS-Aura_L1BRADG_v04-20-c01_2007d207.h5"
 
   ! To convert profile to MAF number, Pranjit Saha
   character(len=100) :: l2GP = "/data/emls/l2gp/v04.20/2007/207/MLS-Aura_L2GP-O3_v04-20-c01_2007d207.he5"
   ! To read ptan values, Pranjit Saha
   character(len=100) :: l2dgm = "/data/emls/l2aux/v04.20/2007/207/MLS-Aura_L2AUX-DGM_v04-20-c01_2007d207.h5"
 
   character(len=100) :: spectroscopy = '/data/emls/l2cal/MLS-Aura_L2Cal-Spectroscopy-PFA_v4-0-3_0000d000.h5'
   character(len=100) :: leapsecFile = '/software/toolkit/LF6.1/toolkit/database/common/TD/leapsec.dat'
   character(len=100) :: antennaPatterns = '/data/emls/l2cal/MLS-Aura_L2Cal-AAAP_v2-0-0_0000d000.txt'
   character(len=100) :: filterShapes = '/data/emls/l2cal/MLS-Aura_L2Cal-Filters_v3-0-2_0000d000.txt'
   character(len=100) :: DACSFilterShapes = '/data/emls/l2cal/MLS-Aura_L2Cal-DACSFilters_v1-5-1_0000d000.txt'
   character(len=100) :: pointingGrids = '/data/emls/l2cal/MLS-Aura_L2Cal-PFG_v3-0-0_0000d000.txt'
   character(len=100), dimension(8) :: &
   pfaFiles= (/ &
   '/data/emls/l2cal/PFA/PFA_R5V_CA-08.h5                                                                   ', &
   '/data/emls/l2cal/PFA/PFA_R5H_CA-08.h5                                                                   ', &
   '/data/emls/l2cal/PFA/PFA_R4_CA-08.h5                                                                    ', &
   '/data/emls/l2cal/PFA/PFA_R3_CA-08.h5                                                                    ', &
   '/data/emls/l2cal/PFA/PFA_R2_CA-08.h5                                                                    ', &
   '/data/emls/l2cal/PFA/PFA_R1B_CA-08.h5                                                                   ', &
   '/data/emls/l2cal/PFA/PFA_R1A_CA-08.h5                                                                   ', &
   '/data/emls/l2cal/PFA/PFA_DACS_CA-08.h5                                                                  ' /)
   character(len=100), dimension(1) :: &
   l2pc = (/ &
   '/data/emls/l2cal/l2pc_040CA08/MLS-Aura_L2Cal-L2PC-band7-LATSCALARHIRESO3HR_v4-00-CA-08_m02.h5       '/)
   real(r8), dimension(55) :: &
   TemperatureInput= (/       300.150_r8,       289.918_r8,       281.060_r8,       272.966_r8,       265.034_r8, &
                255.851_r8,       237.844_r8,       226.169_r8,       215.611_r8,       207.645_r8, &
                204.001_r8,       201.055_r8,       198.594_r8,       198.365_r8,       200.322_r8, &
                203.076_r8,       206.759_r8,       209.636_r8,       211.445_r8,       213.604_r8, &
                217.958_r8,       223.254_r8,       226.517_r8,       229.135_r8,       231.166_r8, &
                232.312_r8,       232.046_r8,       233.874_r8,       235.821_r8,       238.425_r8, &
                242.082_r8,       248.006_r8,       254.357_r8,       259.690_r8,       263.520_r8, &
                265.752_r8,       264.551_r8,       259.498_r8,       252.170_r8,       246.055_r8, &
                242.513_r8,       239.993_r8,       230.524_r8,       212.998_r8,       210.112_r8, &
                191.580_r8,       173.331_r8,       161.116_r8,       151.806_r8,       151.632_r8, &
                169.061_r8,       202.305_r8,       259.684_r8,       356.545_r8,       356.545_r8 /)
   real(r8), dimension(65) :: &
   O3Input= (/   3.26427E-08_r8,   3.28156E-08_r8,   3.48940E-08_r8,   3.71038E-08_r8,   3.68163E-08_r8, &
            3.65310E-08_r8,   3.59035E-08_r8,   3.52867E-08_r8,   3.53855E-08_r8,   3.54845E-08_r8, &
            3.54941E-08_r8,   3.55037E-08_r8,   3.56235E-08_r8,   3.57436E-08_r8,   3.64826E-08_r8, &
            3.72368E-08_r8,   3.89248E-08_r8,   4.06891E-08_r8,   4.50273E-08_r8,   4.98279E-08_r8, &
            5.51404E-08_r8,   7.05450E-08_r8,   9.02488E-08_r8,   1.15459E-07_r8,   1.57574E-07_r8, &
            2.15050E-07_r8,   2.93494E-07_r8,   3.92605E-07_r8,   5.25188E-07_r8,   7.02539E-07_r8, &
            9.09776E-07_r8,   1.17811E-06_r8,   1.52559E-06_r8,   1.92305E-06_r8,   2.42406E-06_r8, &
            3.05557E-06_r8,   3.57256E-06_r8,   4.17701E-06_r8,   4.88376E-06_r8,   5.47605E-06_r8, &
            6.14028E-06_r8,   6.88496E-06_r8,   7.36978E-06_r8,   7.88875E-06_r8,   8.44427E-06_r8, &
            8.68432E-06_r8,   8.93120E-06_r8,   9.18508E-06_r8,   9.16914E-06_r8,   9.15323E-06_r8, &
            9.13734E-06_r8,   9.02868E-06_r8,   8.81521E-06_r8,   8.67289E-06_r8,   8.62597E-06_r8, &
            6.76798E-06_r8,   5.99496E-06_r8,   5.31021E-06_r8,   4.42692E-06_r8,   3.76995E-06_r8, &
            2.96270E-06_r8,   2.14861E-06_r8,   1.55006E-06_r8,   1.11826E-06_r8,   7.23547E-07_r8 /)
   real(r8), dimension(37) :: &
   SO2Input= (/   0.00000E+00_r8,   0.00000E+00_r8,   0.00000E+00_r8,  -1.64908E-08_r8,  -6.87417E-09_r8, &
           -1.37411E-08_r8,   2.64650E-09_r8,  -4.60870E-09_r8,   2.23168E-09_r8,   8.17165E-09_r8, &
           -5.64156E-09_r8,   9.52680E-09_r8,  -1.34418E-08_r8,   8.53389E-09_r8,   8.53559E-09_r8, &
            4.26519E-09_r8,  -6.27082E-09_r8,   2.57639E-09_r8,  -2.73861E-10_r8,  -7.10464E-09_r8, &
           -1.48396E-09_r8,  -2.60092E-09_r8,  -5.06087E-09_r8,  -2.80760E-09_r8,  -2.50622E-09_r8, &
           -2.79323E-09_r8,  -1.26484E-09_r8,  -5.31186E-10_r8,  -2.73478E-10_r8,  -1.27632E-10_r8, &
           -6.47193E-11_r8,   0.00000E+00_r8,   0.00000E+00_r8,   0.00000E+00_r8,   0.00000E+00_r8, &
            0.00000E+00_r8,   0.00000E+00_r8 /)
   real(r8), dimension(37) :: &
   HNO3Input= (/   2.02937E-11_r8,   3.94648E-11_r8,   5.72974E-11_r8,  -1.93186E-09_r8,  -6.88858E-10_r8, &
            1.36816E-09_r8,   4.68370E-10_r8,   5.37981E-10_r8,   3.25219E-09_r8,   3.10499E-09_r8, &
            3.83946E-09_r8,   3.17680E-09_r8,   1.85473E-09_r8,   1.80085E-09_r8,   8.60910E-10_r8, &
            8.23724E-10_r8,   1.03454E-09_r8,  -6.18602E-10_r8,   1.91962E-09_r8,  -1.21803E-09_r8, &
           -3.40837E-09_r8,  -2.13270E-09_r8,  -7.14335E-09_r8,  -8.39571E-09_r8,  -3.17500E-09_r8, &
           -2.79427E-09_r8,  -2.93694E-09_r8,  -3.21508E-10_r8,   1.23244E-14_r8,   3.11581E-16_r8, &
            2.06287E-18_r8,   2.06286E-18_r8,   2.06286E-18_r8,   2.06286E-18_r8,   2.06286E-18_r8, &
            2.06286E-18_r8,   2.06286E-18_r8 /)
   real(r8), dimension(55) :: &
   extinctionV2R3Input= (/   1.00000E-02_r8,   1.00000E-02_r8,   9.99725E-03_r8,   9.00421E-03_r8,   7.88025E-03_r8, &
            9.37069E-03_r8,   5.20747E-03_r8,   2.90742E-03_r8,   1.55356E-03_r8,   8.82107E-04_r8, &
            5.40313E-04_r8,   3.55975E-04_r8,   2.34633E-04_r8,   1.51644E-04_r8,   9.79357E-05_r8, &
            6.20004E-05_r8,   4.12940E-05_r8,   2.62466E-05_r8,   1.55667E-05_r8,   9.49930E-06_r8, &
            7.08288E-06_r8,   5.16998E-06_r8,   2.21251E-06_r8,   1.70082E-06_r8,  -7.08658E-07_r8, &
            1.29696E-06_r8,  -1.68620E-06_r8,  -8.84385E-07_r8,  -1.17012E-06_r8,  -4.95608E-07_r8, &
            1.85951E-07_r8,   7.51485E-07_r8,  -3.71730E-07_r8,  -6.62393E-07_r8,  -1.10324E-06_r8, &
           -5.16387E-07_r8,  -6.32101E-07_r8,  -8.84112E-08_r8,  -2.76626E-07_r8,  -5.19296E-07_r8, &
           -2.95020E-07_r8,  -3.76369E-08_r8,  -6.60811E-08_r8,  -8.74369E-09_r8,   1.78002E-10_r8, &
           -2.84114E-10_r8,   2.23855E-10_r8,   1.03013E-10_r8,   1.03013E-10_r8,   1.03013E-10_r8, &
            1.03013E-10_r8,   1.03013E-10_r8,   1.03013E-10_r8,   1.03013E-10_r8,   1.03013E-10_r8 /)
   real(r8) :: refGPHInput = 16514.4_r8
   real(r8), dimension(:), pointer ::  H2OInput
end module
