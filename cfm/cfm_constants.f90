! Copyright 2011, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

! This module is to contain literal constants, which are data MLS uses in
! its retrieval and forward model. Other modules of CFM will use this
! module, and users of CFM are also welcome to use values out of this module
! if they don't wish to supply their own values.
module CFM_Constants_m
    use MLSCommon, only: r8

    implicit none

!---------------------------- RCS Ident Info -------------------------------
    character(len=*), private, parameter :: ModuleName= &
        "$RCSfile$"
    private :: not_used_here
!---------------------------------------------------------------------------

    public

!------------------------- Date-Time constants -----------------------------
    integer, parameter :: CCSDSLen = 27
!---------------------------------------------------------------------------

!------------------- Emperical Geometry constants --------------------------
    integer, parameter :: empiricalGeometry_noIterations = 10
    real(r8), dimension(21), parameter :: empiricalGeometry_terms = (/ &
    -1.06863, 43.0943, -16.2062, 8.12730, -4.58416, 2.75786, -1.72880, &
     1.11523, -0.733464, 0.489792, -0.331852, 0.227522, -0.156428, &
     0.108031, -0.0757825, 0.0536980, -0.0375161, 0.0260555, &
    -0.0188811, 0.0138453, -0.00959350 /)
!---------------------------------------------------------------------------

!------------------------------ O2 constants -------------------------------
    real(r8), dimension(17), parameter :: o2_heights = &
        (/1.0e+03_r8, 8.0131e-03_r8, 5.8925e-03_r8, 4.3241e-03_r8, 3.1594e-03_r8, &
        2.2961e-03_r8, 1.6581e-03_r8, 1.1874e-03_r8, 8.4392e-04_r8, 5.9869e-04_r8, &
        4.2472e-04_r8, 3.0332e-04_r8, 2.1863e-04_r8, 1.5948e-04_r8, 1.1809e-04_r8, &
        8.8552e-05_r8, 6.6696e-05_r8/)
    real(r8), dimension(17), parameter :: o2_values = &
        (/0.2095_r8, 0.2095_r8, 0.2092_r8, 0.2089_r8, 0.2086_r8, 0.2083_r8, &
        0.2080_r8, 0.2070_r8, 0.2061_r8, 0.2051_r8, 0.2042_r8, 0.2032_r8, &
        0.1915_r8, 0.1798_r8, 0.1681_r8, 0.1564_r8, 0.1447_r8/)
!---------------------------------------------------------------------------

!--------------------- Earth reflectivity constants ------------------------
    real(r8), dimension(1), parameter :: earthRefl_values = (/0.05_r8/)
!---------------------------------------------------------------------------

!----------------------- Space radiance constants --------------------------
    real(r8), dimension(1), parameter :: spaceRad_values = (/2.735_r8/)
!---------------------------------------------------------------------------

!------------------------ Band signal constants ----------------------------
    character(len=20), parameter :: band1L = "R1A:118.B1LF:PT     "
    character(len=20), parameter :: band2L = "R2:190.B2LF:H2O     "
    character(len=20), parameter :: band3L = "R2:190.B3LF:N2O     "
    character(len=20), parameter :: band4L = "R2:190.B4LF:HNO3    "
    character(len=20), parameter :: band5L = "R2:190.B5LF:ClO     "
    character(len=20), parameter :: band6L = "R2:190.B6LF:O3      "
    character(len=20), parameter :: band7L = "R3:240.B7LF:O3      "
    character(len=20), parameter :: band8L = "R3:240.B8LF:PT      "
    character(len=20), parameter :: band9L = "R3:240.B9LF:CO      "
    character(len=20), parameter :: band10L = "R4:640.B10LF:ClO    "
    character(len=20), parameter :: band11L = "R4:640.B11LF:BrO    "
    character(len=20), parameter :: band12L = "R4:640.B12LF:N2O    "
    character(len=20), parameter :: band13L = "R4:640.B13LF:HCl    "
    character(len=20), parameter :: band14L = "R4:640.B14LF:O3     "
    character(len=20), parameter :: band15L = "R5H:2T5.B15LF:OH    "
    character(len=20), parameter :: band16L = "R5H:2T5.B16LF:OH    "
    character(len=20), parameter :: band17L = "R5H:2T5.B17LF:PT    "
    character(len=20), parameter :: band18L = "R5V:2T5.B18LF:OH    "
    character(len=20), parameter :: band19L = "R5V:2T5.B19LF:OH    "
    character(len=20), parameter :: band20L = "R5V:2T5.B20LF:PT    "
    character(len=20), parameter :: band21L = "R1B:118.B21LF:PT    "
    character(len=20), parameter :: band22L = "R1A:118.B22LD:PT    "
    character(len=20), parameter :: band23L = "R2:190.B23LD:H2O    "
    character(len=20), parameter :: band24L = "R3:240.B24LD:O3     "
    character(len=20), parameter :: band25L = "R3:240.B25LD:CO     "
    character(len=20), parameter :: band26L = "R1B:118.B26LD:PT    "
    character(len=20), parameter :: band27L = "R2:190.B27LM:HCN    "
    character(len=20), parameter :: band28L = "R4:640.B28LM:HO2    "
    character(len=20), parameter :: band29L = "R4:640.B29LM:HOCl   "
    character(len=20), parameter :: band30L = "R4:640.B30LM:HO2    "
    character(len=20), parameter :: band31L = "R4:640.B31LM:BrO    "
    character(len=20), parameter :: band32L = "R1A:118.B32LW:PT    "
    character(len=20), parameter :: band33L = "R3:240.B33LW:O3     "
    character(len=20), parameter :: band34L = "R1B:118.B34LW:PT    "

    character(len=20), parameter :: band2U = "R2:190.B2UF:H2O     "
    character(len=20), parameter :: band3U = "R2:190.B3UF:N2O     "
    character(len=20), parameter :: band4U = "R2:190.B4UF:HNO3    "
    character(len=20), parameter :: band5U = "R2:190.B5UF:ClO     "
    character(len=20), parameter :: band6U = "R2:190.B6UF:O3      "
    character(len=20), parameter :: band7U = "R3:240.B7UF:O3      "
    character(len=20), parameter :: band8U = "R3:240.B8UF:PT      "
    character(len=20), parameter :: band9U = "R3:240.B9UF:CO      "
    character(len=20), parameter :: band10U = "R4:640.B10UF:ClO    "
    character(len=20), parameter :: band11U = "R4:640.B11UF:BrO    "
    character(len=20), parameter :: band12U = "R4:640.B12UF:N2O    "
    character(len=20), parameter :: band13U = "R4:640.B13UF:HCl    "
    character(len=20), parameter :: band14U = "R4:640.B14UF:O3     "
    character(len=20), parameter :: band15U = "R5H:2T5.B15UF:OH    "
    character(len=20), parameter :: band16U = "R5H:2T5.B16UF:OH    "
    character(len=20), parameter :: band17U = "R5H:2T5.B17UF:PT    "
    character(len=20), parameter :: band18U = "R5V:2T5.B18UF:OH    "
    character(len=20), parameter :: band19U = "R5V:2T5.B19UF:OH    "
    character(len=20), parameter :: band20U = "R5V:2T5.B20UF:PT    "
    character(len=20), parameter :: band23U = "R2:190.B23UD:H2O    "
    character(len=20), parameter :: band24U = "R3:240.B24UD:O3     "
    character(len=20), parameter :: band25U = "R3:240.B25UD:CO     "
    character(len=20), parameter :: band27U = "R2:190.B27UM:HCN    "
    character(len=20), parameter :: band28U = "R4:640.B28UM:HO2    "
    character(len=20), parameter :: band29U = "R4:640.B29UM:HOCl   "
    character(len=20), parameter :: band30U = "R4:640.B30UM:HO2    "
    character(len=20), parameter :: band31U = "R4:640.B31UM:BrO    "
    character(len=20), parameter :: band33U = "R3:240.B33UW:O3     "
!---------------------------------------------------------------------------

!---------------------- Elevation offset constants -------------------------
    ! spread fill value
    real(r8), parameter :: velev1L = 0.0187_r8
    real(r8), parameter :: velev2L = 0.0151_r8
    real(r8), parameter :: velev3L = 0.0150_r8
    real(r8), parameter :: velev4L = 0.0148_r8
    real(r8), parameter :: velev5L = 0.0146_r8
    real(r8), parameter :: velev6L = 0.0144_r8
    real(r8), parameter :: velev7L = -0.0002_r8
    real(r8), parameter :: velev8L = -0.0000_r8
    real(r8), parameter :: velev9L = -0.0004_r8
    real(r8), parameter :: velev10L = 0.0077_r8
    real(r8), parameter :: velev11L = 0.0077_r8
    real(r8), parameter :: velev12L = 0.0077_r8
    real(r8), parameter :: velev13L = 0.0076_r8
    real(r8), parameter :: velev14L = 0.0076_r8
    real(r8), parameter :: velev15L = 0.0000_r8
    real(r8), parameter :: velev16L = 0.0000_r8
    real(r8), parameter :: velev17L = 0.0000_r8
    real(r8), parameter :: velev18L = 0.0000_r8
    real(r8), parameter :: velev19L = 0.0000_r8
    real(r8), parameter :: velev20L = 0.0000_r8
    real(r8), parameter :: velev21L = 0.0056_r8
    real(r8), parameter :: velev22L = 0.0187_r8
    real(r8), parameter :: velev23L = 0.0151_r8
    real(r8), parameter :: velev24L = -0.0002_r8
    real(r8), parameter :: velev25L = -0.0004_r8
    real(r8), parameter :: velev26L = 0.0056_r8
    real(r8), parameter :: velev27L = 0.0144_r8
    real(r8), parameter :: velev28L = 0.0077_r8
    real(r8), parameter :: velev29L = 0.0077_r8
    real(r8), parameter :: velev30L = 0.0076_r8
    real(r8), parameter :: velev31L = 0.0076_r8
    ! explitcit fill values
    real(r8), dimension(4), parameter :: velev32L = (/0.0200_r8, 0.0215_r8, 0.0186_r8, 0.0172_r8/)
    real(r8), dimension(4), parameter :: velev33L = (/-0.0002_r8, 0.0003_r8, 0.0003_r8, 0.0001_r8/)
    real(r8), dimension(4), parameter :: velev34L = (/0.0057_r8, 0.0057_r8, 0.0060_r8, 0.0063_r8/)

    ! spread fill value
    real(r8), parameter :: velev2U = 0.0141_r8
    real(r8), parameter :: velev3U = 0.0140_r8
    real(r8), parameter :: velev4U = 0.0138_r8
    real(r8), parameter :: velev5U = 0.0133_r8
    real(r8), parameter :: velev6U = 0.0130_r8
    real(r8), parameter :: velev7U = 0.0001_r8
    real(r8), parameter :: velev8U = 0.0003_r8
    real(r8), parameter :: velev9U = -0.0003_r8
    real(r8), parameter :: velev10U = 0.0078_r8
    real(r8), parameter :: velev11U = 0.0081_r8
    real(r8), parameter :: velev12U = 0.0089_r8
    real(r8), parameter :: velev13U = 0.0079_r8
    real(r8), parameter :: velev14U = 0.0078_r8
    real(r8), parameter :: velev15U = 0.0000_r8
    real(r8), parameter :: velev16U = 0.0000_r8
    real(r8), parameter :: velev17U = 0.0000_r8
    real(r8), parameter :: velev18U = 0.0000_r8
    real(r8), parameter :: velev19U = 0.0000_r8
    real(r8), parameter :: velev20U = 0.0000_r8
    real(r8), parameter :: velev23U = 0.0141_r8
    real(r8), parameter :: velev24U = 0.0001_r8
    real(r8), parameter :: velev25U = -0.0003_r8
    real(r8), parameter :: velev27U = 0.0130_r8
    real(r8), parameter :: velev28U = 0.0079_r8
    real(r8), parameter :: velev29U = 0.0080_r8
    real(r8), parameter :: velev30U = 0.0078_r8
    real(r8), parameter :: velev31U = 0.0078_r8
    ! explicit fill values
    real(r8), dimension(4), parameter :: velev33U = (/0.0001_r8, 0.0002_r8, -0.0000_r8, -0.0001_r8/)
!---------------------------------------------------------------------------

!------------------- Limb sideband fraction constants ----------------------
    real(r8), parameter :: vlimbSidebandFraction1L = 0.97333_r8
    real(r8), dimension(25), parameter :: vlimbSidebandFraction2L = &
    (/0.54213_r8, 0.53927_r8, 0.53572_r8, 0.53244_r8, 0.52971_r8, 0.52693_r8, 0.52451_r8, &
      0.52279_r8, 0.52160_r8, 0.52076_r8, 0.52018_r8, 0.51977_r8, 0.51948_r8, 0.51919_r8, &
      0.51878_r8, 0.51822_r8, 0.51742_r8, 0.51632_r8, 0.51479_r8, 0.51273_r8, 0.51052_r8, &
      0.50846_r8, 0.50613_r8, 0.50368_r8, 0.50160_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction3L = &
    (/0.52287_r8, 0.51885_r8, 0.51507_r8, 0.51216_r8, 0.50998_r8, 0.50797_r8, 0.50635_r8, &
      0.50526_r8, 0.50455_r8, 0.50405_r8, 0.50372_r8, 0.50348_r8, 0.50332_r8, 0.50316_r8, &
      0.50293_r8, 0.50263_r8, 0.50220_r8, 0.50163_r8, 0.50088_r8, 0.49992_r8, 0.49898_r8, &
      0.49819_r8, 0.49741_r8, 0.49676_r8, 0.49637_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction4L = &
    (/0.50272_r8, 0.50216_r8, 0.50153_r8, 0.50096_r8, 0.50049_r8, 0.50002_r8, 0.49958_r8, &
      0.49929_r8, 0.49908_r8, 0.49893_r8, 0.49882_r8, 0.49875_r8, 0.49870_r8, 0.49864_r8, &
      0.49858_r8, 0.49848_r8, 0.49833_r8, 0.49814_r8, 0.49787_r8, 0.49752_r8, 0.49715_r8, &
      0.49683_r8, 0.49652_r8, 0.49629_r8, 0.49626_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction5L = &
    (/0.49469_r8, 0.49534_r8, 0.49601_r8, 0.49659_r8, 0.49706_r8, 0.49753_r8, 0.49794_r8, &
      0.49824_r8, 0.49845_r8, 0.49860_r8, 0.49870_r8, 0.49877_r8, 0.49883_r8, 0.49888_r8, &
      0.49895_r8, 0.49906_r8, 0.49920_r8, 0.49941_r8, 0.49970_r8, 0.50010_r8, 0.50054_r8, &
      0.50099_r8, 0.50149_r8, 0.50207_r8, 0.50259_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction6L = &
    (/0.49089_r8, 0.49044_r8, 0.49002_r8, 0.48970_r8, 0.48946_r8, 0.48925_r8, 0.48909_r8, &
      0.48898_r8, 0.48892_r8, 0.48889_r8, 0.48887_r8, 0.48885_r8, 0.48884_r8, 0.48883_r8, &
      0.48881_r8, 0.48879_r8, 0.48877_r8, 0.48875_r8, 0.48874_r8, 0.48877_r8, 0.48888_r8, &
      0.48907_r8, 0.48946_r8, 0.49023_r8, 0.49136_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction7L = &
    (/0.48358_r8, 0.48470_r8, 0.48580_r8, 0.48669_r8, 0.48739_r8, 0.48807_r8, 0.48865_r8, &
      0.48905_r8, 0.48932_r8, 0.48950_r8, 0.48964_r8, 0.48973_r8, 0.48979_r8, 0.48985_r8, &
      0.48995_r8, 0.49007_r8, 0.49024_r8, 0.49047_r8, 0.49078_r8, 0.49116_r8, 0.49153_r8, &
      0.49178_r8, 0.49191_r8, 0.49171_r8, 0.49100_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction8L = &
    (/0.47597_r8, 0.47481_r8, 0.47369_r8, 0.47277_r8, &
      0.47205_r8, 0.47135_r8, 0.47076_r8, 0.47034_r8, 0.47006_r8, &
      0.46987_r8, 0.46972_r8, 0.46962_r8, 0.46956_r8, 0.46949_r8, &
      0.46940_r8, 0.46927_r8, 0.46907_r8, 0.46882_r8, 0.46846_r8, &
      0.46799_r8, 0.46747_r8, 0.46700_r8, 0.46645_r8, 0.46589_r8, &
      0.46543_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction9L = &
    (/0.47767_r8, 0.47701_r8, 0.47627_r8, 0.47565_r8, &
      0.47517_r8, 0.47471_r8, 0.47433_r8, 0.47408_r8, 0.47390_r8, &
      0.47379_r8, 0.47371_r8, 0.47365_r8, 0.47361_r8, 0.47357_r8, &
      0.47352_r8, 0.47344_r8, 0.47333_r8, 0.47319_r8, 0.47299_r8, &
      0.47274_r8, 0.47248_r8, 0.47224_r8, 0.47195_r8, 0.47165_r8, &
      0.47136_r8/)
    real(r8), parameter :: vlimbSidebandFraction10L = 0.50882_r8
    real(r8), parameter :: vlimbSidebandFraction11L = 0.50882_r8
    real(r8), parameter :: vlimbSidebandFraction12L = 0.50690_r8
    real(r8), parameter :: vlimbSidebandFraction13L = 0.51563_r8
    real(r8), parameter :: vlimbSidebandFraction14L = 0.51563_r8
    real(r8), dimension(25), parameter :: vlimbSidebandFraction15L = &
    (/0.95750_r8, 0.95800_r8, 0.95850_r8, 0.95890_r8, &
      0.95920_r8, 0.95960_r8, 0.95990_r8, 0.96010_r8, 0.96020_r8, &
      0.96030_r8, 0.96040_r8, 0.96040_r8, 0.96050_r8, 0.96050_r8, &
      0.96060_r8, 0.96060_r8, 0.96070_r8, 0.96090_r8, 0.96110_r8, &
      0.96140_r8, 0.96170_r8, 0.96210_r8, 0.96250_r8, 0.96300_r8, &
      0.96350_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction16L = &
    (/0.80400_r8, 0.80870_r8, 0.81290_r8, 0.81670_r8, &
      0.81960_r8, 0.82270_r8, 0.82520_r8, 0.82710_r8, 0.82840_r8, &
      0.82940_r8, 0.83000_r8, 0.83050_r8, 0.83080_r8, 0.83110_r8, &
      0.83160_r8, 0.83220_r8, 0.83310_r8, 0.83440_r8, 0.83630_r8, &
      0.83890_r8, 0.84170_r8, 0.84490_r8, 0.84840_r8, 0.85320_r8, &
      0.85720_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction17L = &
    (/0.42500_r8, 0.41920_r8, 0.41350_r8, 0.40850_r8, &
      0.40450_r8, 0.40030_r8, 0.39680_r8, 0.39440_r8, 0.39270_r8, &
      0.39150_r8, 0.39060_r8, 0.39000_r8, 0.38950_r8, 0.38910_r8, &
      0.38850_r8, 0.38760_r8, 0.38640_r8, 0.38460_r8, 0.38220_r8, &
      0.37860_r8, 0.37490_r8, 0.37090_r8, 0.36590_r8, 0.36010_r8, &
      0.35400_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction18L = &
    (/0.94650_r8, 0.94680_r8, 0.94700_r8, 0.94730_r8, &
      0.94740_r8, 0.94760_r8, 0.94780_r8, 0.94790_r8, 0.94800_r8, &
      0.94800_r8, 0.94800_r8, 0.94810_r8, 0.94810_r8, 0.94810_r8, &
      0.94810_r8, 0.94820_r8, 0.94820_r8, 0.94830_r8, 0.94840_r8, &
      0.94860_r8, 0.94870_r8, 0.94890_r8, 0.94910_r8, 0.94940_r8, &
      0.94970_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction19L = &
    (/0.78500_r8, 0.78980_r8, 0.79410_r8, 0.79800_r8, &
      0.80110_r8, 0.80400_r8, 0.80670_r8, 0.80860_r8, 0.80990_r8, &
      0.81090_r8, 0.81160_r8, 0.81200_r8, 0.81240_r8, 0.81270_r8, &
      0.81320_r8, 0.81380_r8, 0.81470_r8, 0.81600_r8, 0.81810_r8, &
      0.82060_r8, 0.82370_r8, 0.82670_r8, 0.83040_r8, 0.83510_r8, &
      0.83960_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction20L = &
    (/0.41870_r8, 0.41310_r8, 0.40650_r8, 0.40200_r8, &
      0.39800_r8, 0.39370_r8, 0.39040_r8, 0.38790_r8, 0.38620_r8, &
      0.38500_r8, 0.38410_r8, 0.38350_r8, 0.38300_r8, 0.38260_r8, &
      0.38200_r8, 0.38120_r8, 0.37990_r8, 0.37810_r8, 0.37570_r8, &
      0.37200_r8, 0.36850_r8, 0.36450_r8, 0.35930_r8, 0.35380_r8, &
      0.34760_r8/)
    real(r8), parameter :: vlimbSidebandFraction21L = 0.94708_r8
    real(r8), parameter :: vlimbSidebandFraction22L = 0.97333_r8
    real(r8), dimension(129), parameter :: vlimbSidebandFraction23L = &
    (/0.51922_r8, 0.51923_r8, 0.51923_r8, 0.51923_r8, &
      0.51924_r8, 0.51924_r8, 0.51924_r8, 0.51925_r8, 0.51925_r8, &
      0.51925_r8, 0.51926_r8, 0.51926_r8, 0.51926_r8, 0.51927_r8, &
      0.51927_r8, 0.51927_r8, 0.51928_r8, 0.51928_r8, 0.51928_r8, &
      0.51929_r8, 0.51929_r8, 0.51929_r8, 0.51930_r8, 0.51930_r8, &
      0.51930_r8, 0.51931_r8, 0.51931_r8, 0.51931_r8, 0.51932_r8, &
      0.51932_r8, 0.51932_r8, 0.51933_r8, 0.51933_r8, 0.51933_r8, &
      0.51934_r8, 0.51934_r8, 0.51934_r8, 0.51935_r8, 0.51935_r8, &
      0.51935_r8, 0.51936_r8, 0.51936_r8, 0.51936_r8, 0.51937_r8, &
      0.51937_r8, 0.51937_r8, 0.51938_r8, 0.51938_r8, 0.51938_r8, &
      0.51939_r8, 0.51939_r8, 0.51939_r8, 0.51940_r8, 0.51940_r8, &
      0.51940_r8, 0.51941_r8, 0.51941_r8, 0.51941_r8, 0.51942_r8, &
      0.51942_r8, 0.51942_r8, 0.51943_r8, 0.51943_r8, 0.51943_r8, &
      0.51944_r8, 0.51944_r8, 0.51944_r8, 0.51945_r8, 0.51945_r8, &
      0.51945_r8, 0.51946_r8, 0.51946_r8, 0.51946_r8, 0.51947_r8, &
      0.51947_r8, 0.51947_r8, 0.51948_r8, 0.51948_r8, 0.51948_r8, &
      0.51949_r8, 0.51949_r8, 0.51949_r8, 0.51950_r8, 0.51950_r8, &
      0.51950_r8, 0.51951_r8, 0.51951_r8, 0.51951_r8, 0.51952_r8, &
      0.51952_r8, 0.51952_r8, 0.51953_r8, 0.51953_r8, 0.51953_r8, &
      0.51954_r8, 0.51954_r8, 0.51954_r8, 0.51955_r8, 0.51955_r8, &
      0.51955_r8, 0.51956_r8, 0.51956_r8, 0.51957_r8, 0.51957_r8, &
      0.51957_r8, 0.51958_r8, 0.51958_r8, 0.51958_r8, 0.51959_r8, &
      0.51959_r8, 0.51959_r8, 0.51960_r8, 0.51960_r8, 0.51960_r8, &
      0.51961_r8, 0.51961_r8, 0.51961_r8, 0.51962_r8, 0.51962_r8, &
      0.51962_r8, 0.51963_r8, 0.51963_r8, 0.51963_r8, 0.51964_r8, &
      0.51964_r8, 0.51964_r8, 0.51965_r8, 0.51965_r8, 0.51965_r8/)
    real(r8), dimension(129), parameter :: vlimbSidebandFraction24L = &
    (/0.48982_r8, 0.48982_r8, 0.48982_r8, 0.48982_r8, &
      0.48982_r8, 0.48982_r8, 0.48982_r8, 0.48982_r8, 0.48982_r8, &
      0.48982_r8, 0.48982_r8, 0.48982_r8, 0.48982_r8, 0.48982_r8, &
      0.48982_r8, 0.48982_r8, 0.48982_r8, 0.48981_r8, 0.48981_r8, &
      0.48981_r8, 0.48981_r8, 0.48981_r8, 0.48981_r8, 0.48981_r8, &
      0.48981_r8, 0.48981_r8, 0.48981_r8, 0.48981_r8, 0.48981_r8, &
      0.48981_r8, 0.48981_r8, 0.48981_r8, 0.48981_r8, 0.48981_r8, &
      0.48981_r8, 0.48981_r8, 0.48981_r8, 0.48981_r8, 0.48981_r8, &
      0.48981_r8, 0.48980_r8, 0.48980_r8, 0.48980_r8, 0.48980_r8, &
      0.48980_r8, 0.48980_r8, 0.48980_r8, 0.48980_r8, 0.48980_r8, &
      0.48980_r8, 0.48980_r8, 0.48980_r8, 0.48980_r8, 0.48980_r8, &
      0.48980_r8, 0.48980_r8, 0.48980_r8, 0.48980_r8, 0.48980_r8, &
      0.48980_r8, 0.48980_r8, 0.48980_r8, 0.48980_r8, 0.48979_r8, &
      0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, &
      0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, &
      0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, &
      0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, &
      0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, &
      0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, &
      0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, &
      0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, &
      0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, &
      0.48978_r8, 0.48978_r8, 0.48978_r8, 0.48978_r8, 0.48978_r8, &
      0.48978_r8, 0.48978_r8, 0.48978_r8, 0.48978_r8, 0.48978_r8, &
      0.48978_r8, 0.48978_r8, 0.48978_r8, 0.48978_r8, 0.48978_r8, &
      0.48978_r8, 0.48978_r8, 0.48978_r8, 0.48978_r8, 0.48978_r8/)
    real(r8), dimension(129), parameter :: vlimbSidebandFraction25L = &
    (/0.47355_r8, 0.47355_r8, 0.47355_r8, 0.47355_r8, 0.47355_r8, 0.47355_r8, 0.47355_r8, &
      0.47355_r8, 0.47355_r8, 0.47355_r8, 0.47356_r8, 0.47356_r8, 0.47356_r8, 0.47356_r8, &
      0.47356_r8, 0.47356_r8, 0.47356_r8, 0.47356_r8, 0.47356_r8, 0.47356_r8, 0.47356_r8, &
      0.47356_r8, 0.47356_r8, 0.47356_r8, 0.47357_r8, 0.47357_r8, 0.47357_r8, 0.47357_r8, &
      0.47357_r8, 0.47357_r8, 0.47357_r8, 0.47357_r8, 0.47357_r8, 0.47357_r8, 0.47357_r8, &
      0.47357_r8, 0.47357_r8, 0.47358_r8, 0.47358_r8, 0.47358_r8, 0.47358_r8, 0.47358_r8, &
      0.47358_r8, 0.47358_r8, 0.47358_r8, 0.47358_r8, 0.47358_r8, 0.47358_r8, 0.47358_r8, &
      0.47358_r8, 0.47359_r8, 0.47359_r8, 0.47359_r8, 0.47359_r8, 0.47359_r8, 0.47359_r8, &
      0.47359_r8, 0.47359_r8, 0.47359_r8, 0.47359_r8, 0.47359_r8, 0.47359_r8, 0.47359_r8, &
      0.47360_r8, 0.47360_r8, 0.47360_r8, 0.47360_r8, 0.47360_r8, 0.47360_r8, 0.47360_r8, &
      0.47360_r8, 0.47360_r8, 0.47360_r8, 0.47360_r8, 0.47360_r8, 0.47360_r8, 0.47361_r8, &
      0.47361_r8, 0.47361_r8, 0.47361_r8, 0.47361_r8, 0.47361_r8, 0.47361_r8, 0.47361_r8, &
      0.47361_r8, 0.47361_r8, 0.47361_r8, 0.47361_r8, 0.47361_r8, 0.47362_r8, 0.47362_r8, &
      0.47362_r8, 0.47362_r8, 0.47362_r8, 0.47362_r8, 0.47362_r8, 0.47362_r8, 0.47362_r8, &
      0.47362_r8, 0.47362_r8, 0.47362_r8, 0.47362_r8, 0.47362_r8, 0.47363_r8, 0.47363_r8, &
      0.47363_r8, 0.47363_r8, 0.47363_r8, 0.47363_r8, 0.47363_r8, 0.47363_r8, 0.47363_r8, &
      0.47363_r8, 0.47363_r8, 0.47363_r8, 0.47363_r8, 0.47364_r8, 0.47364_r8, 0.47364_r8, &
      0.47364_r8, 0.47364_r8, 0.47364_r8, 0.47364_r8, 0.47364_r8, 0.47364_r8, 0.47364_r8, &
      0.47364_r8, 0.47364_r8, 0.47364_r8/)
    real(r8), parameter :: vlimbSidebandFraction26L = 0.94708_r8
    real(r8), dimension(11), parameter :: vlimbSidebandFraction27L = &
    (/0.48916_r8, 0.48929_r8, 0.48940_r8, 0.48948_r8, 0.48955_r8, 0.48959_r8, 0.48963_r8, &
      0.48971_r8, 0.48981_r8, 0.48998_r8, 0.49023_r8/)
    real(r8), parameter :: vlimbSidebandFraction28L = 0.50882_r8
    real(r8), parameter :: vlimbSidebandFraction29L = 0.50882_r8
    real(r8), parameter :: vlimbSidebandFraction30L = 0.51563_r8
    real(r8), parameter :: vlimbSidebandFraction31L = 0.51563_r8
    real(r8), dimension(4), parameter :: vlimbSidebandFraction32L = &
    (/0.97293_r8, 0.97144_r8, 0.97015_r8, 0.96945_r8/)
    real(r8), dimension(4), parameter :: vlimbSidebandFraction33L = &
    (/0.48936_r8, 0.48045_r8, 0.46610_r8, 0.46746_r8/)
    real(r8), dimension(4), parameter :: vlimbSidebandFraction34L = &
    (/0.95080_r8, 0.94541_r8, 0.94532_r8, 0.94825_r8/)

    real(r8), dimension(25), parameter :: vlimbSidebandFraction2U = &
    (/0.43731_r8, 0.44019_r8, 0.44373_r8, 0.44702_r8, 0.44975_r8, 0.45254_r8, 0.45496_r8, &
      0.45668_r8, 0.45787_r8, 0.45871_r8, 0.45929_r8, 0.45970_r8, 0.45999_r8, 0.46028_r8, &
      0.46069_r8, 0.46126_r8, 0.46205_r8, 0.46316_r8, 0.46469_r8, 0.46675_r8, 0.46896_r8, &
      0.47102_r8, 0.47335_r8, 0.47581_r8, 0.47789_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction3U = &
    (/0.45660_r8, 0.46062_r8, 0.46440_r8, 0.46732_r8, 0.46950_r8, 0.47151_r8, 0.47314_r8, &
      0.47422_r8, 0.47494_r8, 0.47544_r8, 0.47577_r8, 0.47601_r8, 0.47616_r8, 0.47633_r8, &
      0.47656_r8, 0.47686_r8, 0.47729_r8, 0.47786_r8, 0.47861_r8, 0.47958_r8, 0.48052_r8, &
      0.48130_r8, 0.48208_r8, 0.48274_r8, 0.48312_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction4U = &
    (/0.47677_r8, 0.47733_r8, 0.47796_r8, 0.47853_r8, 0.47900_r8, 0.47948_r8, 0.47991_r8, &
      0.48020_r8, 0.48042_r8, 0.48056_r8, 0.48067_r8, 0.48074_r8, 0.48079_r8, 0.48085_r8, &
      0.48092_r8, 0.48102_r8, 0.48116_r8, 0.48136_r8, 0.48162_r8, 0.48198_r8, 0.48235_r8, &
      0.48266_r8, 0.48298_r8, 0.48320_r8, 0.48324_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction5U = &
    (/0.48466_r8, 0.48401_r8, 0.48333_r8, 0.48276_r8, 0.48229_r8, 0.48181_r8, 0.48140_r8, &
      0.48110_r8, 0.48089_r8, 0.48075_r8, 0.48064_r8, 0.48057_r8, 0.48051_r8, 0.48046_r8, &
      0.48039_r8, 0.48029_r8, 0.48014_r8, 0.47993_r8, 0.47964_r8, 0.47924_r8, 0.47880_r8, &
      0.47836_r8, 0.47785_r8, 0.47727_r8, 0.47675_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction6U = &
    (/0.48836_r8, 0.48881_r8, 0.48923_r8, 0.48956_r8, 0.48979_r8, 0.49001_r8, 0.49016_r8, &
      0.49027_r8, 0.49033_r8, 0.49037_r8, 0.49039_r8, 0.49041_r8, 0.49042_r8, 0.49043_r8, &
      0.49045_r8, 0.49047_r8, 0.49049_r8, 0.49051_r8, 0.49052_r8, 0.49049_r8, 0.49038_r8, &
      0.49018_r8, 0.48979_r8, 0.48903_r8, 0.48789_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction7U = &
    (/0.48639_r8, 0.48527_r8, 0.48417_r8, 0.48328_r8, 0.48258_r8, 0.48190_r8, 0.48132_r8, &
      0.48092_r8, 0.48065_r8, 0.48046_r8, 0.48033_r8, 0.48024_r8, 0.48017_r8, 0.48012_r8, &
      0.48002_r8, 0.47990_r8, 0.47973_r8, 0.47949_r8, 0.47918_r8, 0.47881_r8, 0.47844_r8, &
      0.47818_r8, 0.47806_r8, 0.47826_r8, 0.47897_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction8U = &
    (/0.49391_r8, 0.49506_r8, 0.49619_r8, 0.49711_r8, &
      0.49783_r8, 0.49853_r8, 0.49912_r8, 0.49954_r8, 0.49982_r8, &
      0.50001_r8, 0.50016_r8, 0.50025_r8, 0.50032_r8, 0.50039_r8, &
      0.50048_r8, 0.50061_r8, 0.50081_r8, 0.50106_r8, 0.50142_r8, &
      0.50189_r8, 0.50241_r8, 0.50288_r8, 0.50343_r8, 0.50399_r8, &
      0.50446_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction9U = &
    (/0.49235_r8, 0.49301_r8, 0.49375_r8, 0.49437_r8, &
      0.49485_r8, 0.49531_r8, 0.49569_r8, 0.49594_r8, 0.49612_r8, &
      0.49623_r8, 0.49631_r8, 0.49637_r8, 0.49641_r8, 0.49645_r8, &
      0.49651_r8, 0.49658_r8, 0.49669_r8, 0.49684_r8, 0.49703_r8, &
      0.49728_r8, 0.49754_r8, 0.49779_r8, 0.49807_r8, 0.49837_r8, &
      0.49866_r8/)
    real(r8), parameter :: vlimbSidebandFraction10U = 0.46045_r8
    real(r8), parameter :: vlimbSidebandFraction11U = 0.46045_r8
    real(r8), parameter :: vlimbSidebandFraction12U = 0.46232_r8
    real(r8), parameter :: vlimbSidebandFraction13U = 0.45360_r8
    real(r8), parameter :: vlimbSidebandFraction14U = 0.45360_r8
    real(r8), dimension(25), parameter :: vlimbSidebandFraction15U = &
    (/0.04250_r8, 0.04200_r8, 0.04150_r8, 0.04110_r8, &
      0.04080_r8, 0.04040_r8, 0.04010_r8, 0.03990_r8, 0.03980_r8, &
      0.03970_r8, 0.03960_r8, 0.03960_r8, 0.03950_r8, 0.03950_r8, &
      0.03940_r8, 0.03940_r8, 0.03930_r8, 0.03910_r8, 0.03890_r8, &
      0.03860_r8, 0.03830_r8, 0.03790_r8, 0.03750_r8, 0.03700_r8, &
      0.03650_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction16U = &
    (/0.19600_r8, 0.19130_r8, 0.18710_r8, 0.18330_r8, &
      0.18040_r8, 0.17730_r8, 0.17480_r8, 0.17290_r8, 0.17160_r8, &
      0.17060_r8, 0.17000_r8, 0.16950_r8, 0.16920_r8, 0.16890_r8, &
      0.16840_r8, 0.16780_r8, 0.16690_r8, 0.16560_r8, 0.16370_r8, &
      0.16110_r8, 0.15830_r8, 0.15510_r8, 0.15160_r8, 0.14680_r8, &
      0.14280_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction17U = &
    (/0.57500_r8, 0.58080_r8, 0.58650_r8, 0.59150_r8, &
      0.59550_r8, 0.59970_r8, 0.60320_r8, 0.60560_r8, 0.60730_r8, &
      0.60850_r8, 0.60940_r8, 0.61000_r8, 0.61050_r8, 0.61090_r8, &
      0.61150_r8, 0.61240_r8, 0.61360_r8, 0.61540_r8, 0.61780_r8, &
      0.62140_r8, 0.62510_r8, 0.62910_r8, 0.63410_r8, 0.63990_r8, &
      0.64600_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction18U = &
    (/0.05350_r8, 0.05320_r8, 0.05300_r8, 0.05270_r8, &
      0.05260_r8, 0.05240_r8, 0.05220_r8, 0.05210_r8, 0.05200_r8, &
      0.05200_r8, 0.05200_r8, 0.05190_r8, 0.05190_r8, 0.05190_r8, &
      0.05190_r8, 0.05180_r8, 0.05180_r8, 0.05170_r8, 0.05160_r8, &
      0.05140_r8, 0.05130_r8, 0.05110_r8, 0.05090_r8, 0.05060_r8, &
      0.05030_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction19U = &
    (/0.21500_r8, 0.21020_r8, 0.20590_r8, 0.20200_r8, &
      0.19890_r8, 0.19600_r8, 0.19330_r8, 0.19140_r8, 0.19010_r8, &
      0.18910_r8, 0.18840_r8, 0.18800_r8, 0.18760_r8, 0.18730_r8, &
      0.18680_r8, 0.18620_r8, 0.18530_r8, 0.18400_r8, 0.18190_r8, &
      0.17940_r8, 0.17630_r8, 0.17330_r8, 0.16960_r8, 0.16490_r8, &
      0.16040_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction20U = &
    (/0.58130_r8, 0.58690_r8, 0.59350_r8, 0.59800_r8, &
      0.60200_r8, 0.60630_r8, 0.60960_r8, 0.61210_r8, 0.61380_r8, &
      0.61500_r8, 0.61590_r8, 0.61650_r8, 0.61700_r8, 0.61740_r8, &
      0.61800_r8, 0.61880_r8, 0.62010_r8, 0.62190_r8, 0.62430_r8, &
      0.62800_r8, 0.63150_r8, 0.63550_r8, 0.64070_r8, 0.64620_r8, &
      0.65240_r8/)
    real(r8), dimension(129), parameter :: vlimbSidebandFraction23U = &
    (/0.46025_r8, 0.46024_r8, 0.46024_r8, 0.46024_r8, &
      0.46023_r8, 0.46023_r8, 0.46023_r8, 0.46022_r8, 0.46022_r8, &
      0.46022_r8, 0.46021_r8, 0.46021_r8, 0.46021_r8, 0.46020_r8, &
      0.46020_r8, 0.46020_r8, 0.46019_r8, 0.46019_r8, 0.46019_r8, &
      0.46018_r8, 0.46018_r8, 0.46018_r8, 0.46017_r8, 0.46017_r8, &
      0.46017_r8, 0.46016_r8, 0.46016_r8, 0.46016_r8, 0.46015_r8, &
      0.46015_r8, 0.46015_r8, 0.46014_r8, 0.46014_r8, 0.46014_r8, &
      0.46013_r8, 0.46013_r8, 0.46013_r8, 0.46012_r8, 0.46012_r8, &
      0.46012_r8, 0.46011_r8, 0.46011_r8, 0.46011_r8, 0.46010_r8, &
      0.46010_r8, 0.46010_r8, 0.46009_r8, 0.46009_r8, 0.46009_r8, &
      0.46008_r8, 0.46008_r8, 0.46008_r8, 0.46007_r8, 0.46007_r8, &
      0.46007_r8, 0.46006_r8, 0.46006_r8, 0.46006_r8, 0.46005_r8, &
      0.46005_r8, 0.46005_r8, 0.46004_r8, 0.46004_r8, 0.46004_r8, &
      0.46003_r8, 0.46003_r8, 0.46003_r8, 0.46002_r8, 0.46002_r8, &
      0.46002_r8, 0.46001_r8, 0.46001_r8, 0.46001_r8, 0.46000_r8, &
      0.46000_r8, 0.46000_r8, 0.45999_r8, 0.45999_r8, 0.45999_r8, &
      0.45998_r8, 0.45998_r8, 0.45998_r8, 0.45997_r8, 0.45997_r8, &
      0.45997_r8, 0.45996_r8, 0.45996_r8, 0.45996_r8, 0.45995_r8, &
      0.45995_r8, 0.45995_r8, 0.45994_r8, 0.45994_r8, 0.45994_r8, &
      0.45994_r8, 0.45994_r8, 0.45994_r8, 0.45993_r8, 0.45993_r8, &
      0.45993_r8, 0.45992_r8, 0.45992_r8, 0.45991_r8, 0.45991_r8, &
      0.45991_r8, 0.45990_r8, 0.45990_r8, 0.45990_r8, 0.45989_r8, &
      0.45989_r8, 0.45989_r8, 0.45988_r8, 0.45988_r8, 0.45988_r8, &
      0.45987_r8, 0.45987_r8, 0.45987_r8, 0.45986_r8, 0.45986_r8, &
      0.45986_r8, 0.45985_r8, 0.45985_r8, 0.45985_r8, 0.45984_r8, &
      0.45984_r8, 0.45984_r8, 0.45983_r8, 0.45983_r8, 0.45983_r8/)
    real(r8), dimension(129), parameter :: vlimbSidebandFraction24U = &
    (/0.48014_r8, 0.48014_r8, 0.48014_r8, 0.48014_r8, &
      0.48014_r8, 0.48014_r8, 0.48014_r8, 0.48014_r8, 0.48014_r8, &
      0.48014_r8, 0.48014_r8, 0.48014_r8, 0.48014_r8, 0.48014_r8, &
      0.48014_r8, 0.48014_r8, 0.48014_r8, 0.48015_r8, 0.48015_r8, &
      0.48015_r8, 0.48015_r8, 0.48015_r8, 0.48015_r8, 0.48015_r8, &
      0.48015_r8, 0.48015_r8, 0.48015_r8, 0.48015_r8, 0.48015_r8, &
      0.48015_r8, 0.48015_r8, 0.48015_r8, 0.48015_r8, 0.48015_r8, &
      0.48015_r8, 0.48015_r8, 0.48015_r8, 0.48015_r8, 0.48015_r8, &
      0.48015_r8, 0.48016_r8, 0.48016_r8, 0.48016_r8, 0.48016_r8, &
      0.48016_r8, 0.48016_r8, 0.48016_r8, 0.48016_r8, 0.48016_r8, &
      0.48016_r8, 0.48016_r8, 0.48016_r8, 0.48016_r8, 0.48016_r8, &
      0.48016_r8, 0.48016_r8, 0.48016_r8, 0.48016_r8, 0.48016_r8, &
      0.48016_r8, 0.48016_r8, 0.48016_r8, 0.48016_r8, 0.48017_r8, &
      0.48017_r8, 0.48017_r8, 0.48017_r8, 0.48017_r8, 0.48017_r8, &
      0.48017_r8, 0.48017_r8, 0.48017_r8, 0.48017_r8, 0.48017_r8, &
      0.48017_r8, 0.48017_r8, 0.48017_r8, 0.48017_r8, 0.48017_r8, &
      0.48017_r8, 0.48017_r8, 0.48017_r8, 0.48017_r8, 0.48017_r8, &
      0.48017_r8, 0.48017_r8, 0.48018_r8, 0.48018_r8, 0.48018_r8, &
      0.48018_r8, 0.48018_r8, 0.48018_r8, 0.48018_r8, 0.48018_r8, &
      0.48018_r8, 0.48018_r8, 0.48018_r8, 0.48018_r8, 0.48018_r8, &
      0.48018_r8, 0.48018_r8, 0.48018_r8, 0.48018_r8, 0.48018_r8, &
      0.48018_r8, 0.48018_r8, 0.48018_r8, 0.48018_r8, 0.48018_r8, &
      0.48019_r8, 0.48019_r8, 0.48019_r8, 0.48019_r8, 0.48019_r8, &
      0.48019_r8, 0.48019_r8, 0.48019_r8, 0.48019_r8, 0.48019_r8, &
      0.48019_r8, 0.48019_r8, 0.48019_r8, 0.48019_r8, 0.48019_r8, &
      0.48019_r8, 0.48019_r8, 0.48019_r8, 0.48019_r8, 0.48019_r8/)
    real(r8), dimension(129), parameter :: vlimbSidebandFraction25U = &
    (/0.49647_r8, 0.49647_r8, 0.49647_r8, 0.49647_r8, &
      0.49647_r8, 0.49647_r8, 0.49647_r8, 0.49647_r8, 0.49647_r8, &
      0.49647_r8, 0.49646_r8, 0.49646_r8, 0.49646_r8, 0.49646_r8, &
      0.49646_r8, 0.49646_r8, 0.49646_r8, 0.49646_r8, 0.49646_r8, &
      0.49646_r8, 0.49646_r8, 0.49646_r8, 0.49646_r8, 0.49646_r8, &
      0.49645_r8, 0.49645_r8, 0.49645_r8, 0.49645_r8, 0.49645_r8, &
      0.49645_r8, 0.49645_r8, 0.49645_r8, 0.49645_r8, 0.49645_r8, &
      0.49645_r8, 0.49645_r8, 0.49645_r8, 0.49644_r8, 0.49644_r8, &
      0.49644_r8, 0.49644_r8, 0.49644_r8, 0.49644_r8, 0.49644_r8, &
      0.49644_r8, 0.49644_r8, 0.49644_r8, 0.49644_r8, 0.49644_r8, &
      0.49644_r8, 0.49643_r8, 0.49643_r8, 0.49643_r8, 0.49643_r8, &
      0.49643_r8, 0.49643_r8, 0.49643_r8, 0.49643_r8, 0.49643_r8, &
      0.49643_r8, 0.49643_r8, 0.49643_r8, 0.49643_r8, 0.49642_r8, &
      0.49642_r8, 0.49642_r8, 0.49642_r8, 0.49642_r8, 0.49642_r8, &
      0.49642_r8, 0.49642_r8, 0.49642_r8, 0.49642_r8, 0.49642_r8, &
      0.49642_r8, 0.49642_r8, 0.49641_r8, 0.49641_r8, 0.49641_r8, &
      0.49641_r8, 0.49641_r8, 0.49641_r8, 0.49641_r8, 0.49641_r8, &
      0.49641_r8, 0.49641_r8, 0.49641_r8, 0.49641_r8, 0.49641_r8, &
      0.49640_r8, 0.49640_r8, 0.49640_r8, 0.49640_r8, 0.49640_r8, &
      0.49640_r8, 0.49640_r8, 0.49640_r8, 0.49640_r8, 0.49640_r8, &
      0.49640_r8, 0.49640_r8, 0.49640_r8, 0.49640_r8, 0.49639_r8, &
      0.49639_r8, 0.49639_r8, 0.49639_r8, 0.49639_r8, 0.49639_r8, &
      0.49639_r8, 0.49639_r8, 0.49639_r8, 0.49639_r8, 0.49639_r8, &
      0.49639_r8, 0.49639_r8, 0.49638_r8, 0.49638_r8, 0.49638_r8, &
      0.49638_r8, 0.49638_r8, 0.49638_r8, 0.49638_r8, 0.49638_r8, &
      0.49638_r8, 0.49638_r8, 0.49638_r8, 0.49638_r8, 0.49638_r8/)
    real(r8), dimension(11), parameter :: vlimbSidebandFraction27U = &
    (/0.49004_r8, 0.48991_r8, 0.48981_r8, 0.48973_r8, &
      0.48966_r8, 0.48962_r8, 0.48957_r8, 0.48949_r8, 0.48939_r8, &
      0.48923_r8, 0.48897_r8/)
    real(r8), parameter :: vlimbSidebandFraction28U = 0.46045_r8
    real(r8), parameter :: vlimbSidebandFraction29U = 0.46045_r8
    real(r8), parameter :: vlimbSidebandFraction30U = 0.45360_r8
    real(r8), parameter :: vlimbSidebandFraction31U = 0.45360_r8
    real(r8), dimension(4), parameter :: vlimbSidebandFraction33U = &
    (/0.48061_r8, 0.48947_r8, 0.50383_r8, 0.50247_r8/)
!---------------------------------------------------------------------------

    contains
    logical function not_used_here()
        character (len=*), parameter :: IdParm = &
        "$Id$"
        character (len=len(idParm)) :: Id = idParm
        not_used_here = (id(1:1) == ModuleName(1:1))
        print *, Id ! .mod files sometimes change if PRINT is added
    end function not_used_here

end module CFM_Constants_m

! $Log$
! Revision 1.3  2011/10/18 16:56:47  honghanh
! Refractoring limb sideband constants to CFM_Constants_m.
!
