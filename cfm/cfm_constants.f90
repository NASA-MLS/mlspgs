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

    contains
    logical function not_used_here()
        character (len=*), parameter :: IdParm = &
        "$Id$"
        character (len=len(idParm)) :: Id = idParm
        not_used_here = (id(1:1) == ModuleName(1:1))
        print *, Id ! .mod files sometimes change if PRINT is added
    end function not_used_here

end module CFM_Constants_m
