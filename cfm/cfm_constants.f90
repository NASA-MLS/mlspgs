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
    real(r8), dimension(21), parameter :: empircalGeometry_terms = (/ &
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

    contains
    logical function not_used_here()
        character (len=*), parameter :: IdParm = &
        "$Id$"
        character (len=len(idParm)) :: Id = idParm
        not_used_here = (id(1:1) == ModuleName(1:1))
        print *, Id ! .mod files sometimes change if PRINT is added
    end function not_used_here

end module CFM_Constants_m
