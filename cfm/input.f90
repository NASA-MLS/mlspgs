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
   use MLSCommon, only: r8, i4, rv

   implicit none

   ! StartProfile and endProfile will be converted 
   ! to MAF number, Pranjit Saha
       integer :: startProfile , endProfile 
     
       character(len=512) :: l1boa, l1brad, l2GP,l2dgm
       character(len=256) :: spectroscopy, &
        leapsecFile,  antennaPatterns, filterShapes, &
        DACSFilterShapes, pointingGrids
       character(len=256), dimension(8) :: pfaFiles
       !character(len=256), dimension(:), pointer :: pfaFiles
       character(len=256), dimension(1) :: l2pc 
       !character(len=256), dimension(:), pointer :: l2pc 
       
       real(r8) :: refGPHInput
       !! Some variables used by mockup with default values
       real(r8) :: vGridStandard37Start = 1000.0_r8
       real(r8) :: vGridExtinctionStart=1000.0_r8
       real(r8) :: hGridStandardVal1 = 0.0_r8
       real(r8) :: hGridStandardVal2 = 1.5_r8
       real(r8) :: qH2OMinValue = 0.1E-6_r8
       character(len=64) :: vGridStandard37formula = "37:6"
       real(r8), dimension(1) :: vGridRefGPHVals = (/100.0_r8/)
       real(r8), dimension(55) :: H2O_Sample_Vals
       
       character(len=64) :: vGridExtinctionFormula="21:12,14:6,12:3"
       logical :: outputTxtFile = .TRUE.
       integer(i4)  MLS_Year, MLS_Day, MLS_ProfileNumber, &
           Tes_run, Tes_scan, TES_sequnce
           
       real(rv), dimension(:,:), pointer :: &
            Radiance_Calculated => NULL(), &
            Radiance_Observed => NULL(),  &
            Radiance_Noise => NULL(),  &
            Jacobian_Concatenated => NULL(), &
            Jacobian_Concatenated_out => NULL(), &
            ptanGHz_save => NULL()
       
      real(rv), dimension(:), pointer :: &
            Radiance_Calculated_out=> NULL(), &
            Radiance_Noise_out=> NULL(), &
            Radiance_Observed_out=> NULL()
            

!Ming: I changed the following profiles to consist Temperature, ExtinctionV2, O3, 3_V1_3, O3_V2 at 55 levels;
!      HNO3, S_32_O2 at 37 levels; CO at <67 levels.  Note: is O_18_O set somewhere else? Extingction? O3_V1_3 etc?
    real(r8), dimension(:), pointer ::  TemperatureInput, H2OInput , &
    O3Input, SO2Input, HNO3Input, COInput, extinctionV2R3Input , &
    PressureCOInput, PressureStandardInput 

end module
