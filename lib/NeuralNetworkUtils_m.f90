! Copyright 2021, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use  must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.
!
! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

  ! ----------------------------------------------------------------------------
  ! This module contains the utility code to be used when calculating
  ! species using Frank Werner's Neural Network solutions. At the time
  ! of this writing (2021/01/20) only Temperature is being retrieved,
  ! but the code should be able to do any species, provided the correct
  ! training data is provided. 

  ! To adapt to other species, we must revisit the hardcoded use  of some
  ! signal names, arrays sizes, etc. both here and in l2/NeuralNet_m.f90
  !
  ! Note that it currently supports only the retrieval of Temperature
  ! Also it assumes Bands 1, 8, and 22 are supplied
  ! The current weights file format is hard-coded here. 
  ! It would require changes and testing if that format should ever change.
  ! We also assume the weights are stored in a manner consistent
  ! with the following order for the "collapsed" measurement vector:
  !   [b1c1m1,b1c2m1,..,b1c1m2,b1c2m2,..,b2c1m1,..,..,b22c1m1,..]
  ! where
  !   bk is Band k
  !   cj is channel j
  !   mi is MIF i
  ! ----------------------------------------------------------------------------

MODULE NeuralNetUtils_M

  use  Dump_0, only: Dump
  use  HighOutput, only: OutputNamedValue
  use  MLSCommon, only: MLSFile_T, UndefinedValue
  use  MLSHDF5, only: SaveAsHDF5DS
  use  MLSFinds, only: FindFirst
  use  MLSKinds, only: R4, R8, Rv
  use  MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use  Output_M, only: Output 

  implicit none

  type Radiance_T
     character(len=64) :: signal
     !integer :: NumMAFS ! <----- it's only ever going to be one MAF, right?
     integer :: NumChannels
     integer :: NumMIFs
     real(r8), dimension(:,:), allocatable :: values
  end type Radiance_T

  type NeuralNetCoeffs_T
     character(len=32)                             :: Activation_Function
     ! Will be allocated in the caller. 
     ! [3]  Number of Bands
     character(len=32),  dimension(:), allocatable :: Bands
     
     ! The first index is M, the max num of channels used (probably 51)
     ! The second index is the number of Bands (currently 3)
     ! We need this info because  the Neural Net doesn't currently use
     ! every channel in the DACS 
     ! (and someday may even skip some in the non-DACS)
     ! So (M,3)
     integer,   dimension(:,:), allocatable :: Channels_In_Each_Band
     ! Variabls of this type will be allocated in the caller. 

     ! these are the dimensions for each part (assuming that we keep
     ! with Frank set up for Temperature retrieval)

     ! [42], number of pressure levels
     real(r8),  dimension(:), allocatable :: Intercepts_Hidden_Labels_Layer
     real(r8),  dimension(:), allocatable :: Output_Pressure_Levels
     integer,   dimension(:), allocatable :: Output_Pressure_Levels_Indices

     ! [5078]. Number of Neurons in hidden layers
     real(r8),  dimension(:),allocatable :: Intercepts_Hidden_Layer_1
     real(r8),  dimension(:),allocatable :: Intercepts_Hidden_Layer_2

     ! [5078, 42] = [nNeurons, nLevs]
     real(r8),  dimension(:,:),allocatable :: Weights_Hidden_Labels_Layer
     ! [7575,5078] = [nVars, nNeurons]
     real(r8),  dimension(:,:),allocatable :: Weights_Hidden_Layer_1
     ! [nNeurons, nNeurons]
     real(r8),  dimension(:,:),allocatable :: Weights_Hidden_Layer_2


     ! [42]. Number of levels
     real(r8),  dimension(:),allocatable :: Normalization_Labels_Max
     real(r8),  dimension(:),allocatable :: Normalization_Labels_Min

     
     ! [75]. Number of MIFs
     integer,  dimension(:), allocatable :: MIFs

     ! [7575] number of 'variables', i.e. 2*25*75 + 51*75
     ! This is band 1[channel X MIF] + band 8 [channel X MIF ]  and 
     ! band 22 [channel X MIF]
     real(r8),  dimension(:),allocatable :: Standardization_Brightness_Temperature_Mean
     real(r8),  dimension(:),allocatable :: Standardization_Brightness_Temperature_Std

     ! Just for debugging
     real(r8),  dimension(:,:), allocatable :: Stddevs
     real(r8),  dimension(:,:), allocatable :: Means
  end type NeuralNetCoeffs_T


  type NeuralNetInputData_T
     ! I'm assuming that will have been corrected for the baseline and
     ! reduced to the appropriate MIFs, 22-96 for all bands, all 25
     ! channels for Bands 1 and 8 and channels 40-90 for Band
     ! 22. 
     !
     ! All of these will end up having 
     ! These will be allocated in the caller.

     type (Radiance_T) :: Band_1_Radiances 
     type (Radiance_T) :: Band_8_Radiances
     type (Radiance_T) :: Band_22_Radiances 

  end type NeuralNetInputData_T

  ! Paul; I used this type to pass info out of the routine. In the
  ! final configutation, all you'll need is `prediction', so you could
  ! modify the type to use  just that variable, or comment it out
  ! entirely and just put `prediction' directly into the interface for
  ! NeuralNetFit. 
  type NeuralNetOutputData_T
     real(r8),  dimension(:),allocatable::prediction, &
          &  working_space, &
          & neuronValues, &
          & neuronValues2, &
          & y_pred
  end type NeuralNetOutputData_T

  interface Dump
    module procedure Dump_Coeffs !, Dump_Radiance, Dump_InputData
  end interface


  ! In case we're checking against Frank's results
  ! These are purely for debugging against matched results
  ! Comment them out, and references to them, when you're
  ! satisfied our results match.
  real(r8),  dimension(7575), public, save :: matchedStdRadiances = 0
  integer, public, save                    :: matchedbinNum       = 0
  integer, public, save                    :: matchedMAF          = 0

  public :: CheckMAFs, Dump, &
       & NeuralNetInputData_T, & 
       & NeuralNetCoeffs_T, & 
       & NeuralNetFit, & 
       & Radiance_T, &
       & StandardizeRadiances

  
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------


contains 
  subroutine CheckTemperatures ( TempFileID, ProfileRange, TemperatureValues )
    ! Compare with the Temperature values written to the TempFile by
    ! Frank's python script
    use  Intrinsic, only: L_HDF
    use  MLSHDF5, only: LoadFromHDF5DS
    use  MLSFiles, only: HDFVersion_5
    ! Dummy args
    integer, intent(in)                  :: TempFileID ! For optional comparisons
    integer,  dimension(:), intent(in)    :: ProfileRange ! For optional comparisons
    real(rv),   dimension(:)              :: TemperatureValues
    ! Internal variables
    integer                                :: profile
    type (MLSFile_T)                       :: TempFile
    real(r4), allocatable,  dimension(:,:)  :: values2
    ! Executable
    call outputNamedValue( 'Range of profiles to read Temperatures', ProfileRange )
    TempFile%FileID%f_id = TempFileID
    TempFile%StillOpen = .true.
    TempFile%type = l_hdf
    TempFile%hdfVersion = HDFVersion_5
    allocate ( values2(55, 3495) )
    call LoadFromHDF5DS ( TempFile, &
      & "ANN_Prediction", &
      & values2 )
    do profile = ProfileRange(1), ProfileRange(2)
      TemperatureValues = values2(:, profile)
    enddo
  end subroutine CheckTemperatures 

  function CheckMafs ( TempFileID, StdRadiances ) result ( MAF )
    ! Compare with the standardized rads written to the TempFile by
    ! Frank's python script
    ! This comparison seeks to find which MAF he most likely used
    use  Intrinsic, only: L_HDF
    use  MLSHDF5, only: LoadFromHDF5DS
    use  MLSFiles, only: HDFVersion_5
    ! Dummy args
    integer, intent(in)                     :: TempFileID ! For optional comparisons
    real(r8),  dimension(:), intent(in)     :: StdRadiances
    integer                                 :: MAF ! Remember--starts at 0
    ! Internal variables
    integer                                 :: j
    type (MLSFile_T)                        :: TempFile
    real(r4), allocatable,  dimension(:,:)  :: ratios2
    real(r4), allocatable,  dimension(:)    :: recalc
    ! Executable
    TempFile%FileID%f_id = TempFileID
    TempFile%StillOpen = .true.
    TempFile%type = l_hdf
    TempFile%hdfVersion = HDFVersion_5
    ! Recalculate ratios using Frank's own radiances, means, and stddevs
    ! Then check against the ratios he wrote out to TempFile
    allocate ( ratios2(7575, 3495) )
    call LoadFromHDF5DS ( TempFile, &
      & "Standardized_Brightness_Temps_Matrix", &
      & ratios2 )
    MAF = FindFirst ( &
      & (StdRadiances(1)-ratios2(1,:))**2 &
      &  + &
      & (StdRadiances(2)-ratios2(2,:))**2 &
      & < 1.e-6 ) - 1
    call OutputNamedValue ( 'min std rads', minval(StdRadiances) )
    call OutputNamedValue ( 'max std rads', maxval(StdRadiances) )
    call OutputNamedValue ( 'min matched rads', minval(matchedStdRadiances) )
    call OutputNamedValue ( 'max matched rads', maxval(matchedStdRadiances) )
    call OutputNamedValue ( 'min rads(matchedMAF)', minval(ratios2(:,matchedMAF+1)) )
    call OutputNamedValue ( 'max rads(matchedMAF)', maxval(ratios2(:,matchedMAF+1)) )
  end function CheckMafs

  subroutine CheckBinNums ( TempFileID, MAFRange, nnCoeffs )
    ! Compare with the standardized rads written to the TempFile by
    ! Frank's python script
    ! This comparison seeks to find which bin number he most likely used
    use  Intrinsic, only: L_HDF
    use  MLSHDF5, only: LoadFromHDF5DS
    use  MLSFiles, only: HDFVersion_5
    ! Dummy args
    integer, intent(in)                    :: TempFileID ! For optional comparisons
    integer,  dimension(:), intent(in)     :: MAFRange ! For optional comparisons
    type (NeuralNetCoeffs_T), intent(in)   :: nnCoeffs
    ! Internal variables
    integer                                :: binNum
    integer                                :: binNum2
    integer                                :: j
    integer                                :: MAF ! Remember--MAFs start at 0
    type (MLSFile_T)                       :: TempFile
    real(r4), allocatable,  dimension(:,:) :: ratios2 
    real(r4), allocatable,  dimension(:,:) :: recalc  
    real(r4), allocatable,  dimension(:,:) :: values2 
    ! Executable
    call outputNamedValue( 'Range of MAFs to check bin numbers', MAFRange )
    TempFile%FileID%f_id = TempFileID
    TempFile%StillOpen = .true.
    TempFile%type = l_hdf
    TempFile%hdfVersion = HDFVersion_5
    ! Recalculate ratios using Frank's own radiances, means, and stddevs
    ! Then check against the ratios he wrote out to TempFile
    allocate ( ratios2(7575, 3495) )
    allocate ( values2(7575, 3495) )
    call LoadFromHDF5DS ( TempFile, &
      & "Brightness_Temps_Matrix", &
      & values2 )
    call LoadFromHDF5DS ( TempFile, &
      & "Standardized_Brightness_Temps_Matrix", &
      & ratios2 )
    do MAF = MAFRange(1), min(MAFRange(2), size(values2, 2) - 1)
      allocate ( recalc(18,2) )
      do binNum=1, 18
        recalc(binNum,1:2) = &
          & (values2(1:2, MAF+1) - nnCoeffs%means(binNum,1:2)) &
          & / &
          & nnCoeffs%stddevs(binNum,1:2)
      enddo
      ! Now for the musical question:
      ! which binNum most closely approximates Frank's ratio?
      binNum = FindFirst ( (ratios2(1,MAF+1)-recalc(:,1))**2 + &
        &  (ratios2(2,MAF+1)-recalc(:,2))**2 < 1.e-8 )
      if ( binNum < 1 ) then
        call output ( 'No matching binNum found in Franks file', advance='yes' )
      else
        call outputNamedValue ( 'MAF, Bin num matched', (/MAF, BinNum /) )
        call outputNamedValue ( 'Frank,me', (/ratios2(1,MAF+1),recalc(binNum,1) /) )
        j = 1 + mod(binNum,18)
        call outputNamedValue ( 'if wrong Bin', recalc(j,2)  )
        ! Check if a second bin number also matches
        if ( binNum < 18 ) then
          binNum2 = FindFirst ( (ratios2(1,MAF+1)-recalc(binNum+1:,1))**2 + &
            & (ratios2(2,MAF+1)-recalc(binNum+1:,2))**2 < 1.e-8 )
          if ( binNum2 > 0 ) then
            call outputNamedValue ( 'A 2nd match found after 1st', BinNum+BinNum2 )
            ! stop
          endif
        endif
        matchedBinNum = binNum
        ! We match at the first MIF, channel, and Band. How about all the rest?
        deallocate ( recalc )
        allocate ( recalc(7575,1) )
        do j=1, 7575
          recalc(j,1) = &
            & (values2(j, MAF+1) - nnCoeffs%means(binNum,j)) &
            & / &
            & nnCoeffs%stddevs(binNum,j)
        enddo
        call outputNamedValue ( 'max diff over all MIFs, etc.', &
          & maxval( abs(ratios2(:,MAF+1)-recalc(:,1)) ) )
      endif
      deallocate ( recalc )
    enddo

  end subroutine CheckBinNums

  subroutine CompareStandardizedRads ( TempFileID, MAFRange, Rads, &
    & standardized, MatchingMAF, Recalculate, mean, stddev )
    ! Compare with the standardized rads written to the TempFile by
    ! Frank's python script
    use  Intrinsic, only: L_HDF
    use  MLSHDF5, only: LoadFromHDF5DS
    use  MLSFiles, only: HDFVersion_5
    ! Dummy args
    integer, intent(in)                 :: TempFileID ! For optional comparisons
    integer,  dimension(:), intent(in)  :: MAFRange ! For optional comparisons
    real(r8), intent(in),  dimension(:) :: Rads
    logical, intent(in)                 :: standardized
    integer, intent(out), optional      :: MatchingMAF ! starts at 0
    logical, intent(in), optional       :: Recalculate
    real(r8),  dimension(:), intent(in), optional :: Mean
    real(r8),  dimension(:), intent(in), optional :: StdDev
    !
    type (MLSFile_T)                        :: TempFile
    real(r4), allocatable,  dimension(:,:)  :: ratios2
    real(r4), allocatable,  dimension(:)    :: recalc
    real(r4), allocatable,  dimension(:,:)  :: values2
    real(r4),  dimension(size(Rads))        :: diffs
    character(len=32)                       :: radianceType
    integer                                 :: j
    integer                                 :: MAF ! Remember--starts at 0
    integer                                 :: MIF
    logical                                 :: DEEBug
    logical                                 :: myRecalculate
    ! Executable
    if ( present(MatchingMAF) ) MatchingMAF = 0
    myRecalculate = .false.
    if ( present(Recalculate) ) myRecalculate = Recalculate
    call outputNamedValue( 'Range of MAFs to compare', MAFRange )
    TempFile%FileID%f_id = TempFileID
    TempFile%StillOpen = .true.
    TempFile%type = l_hdf
    TempFile%hdfVersion = HDFVersion_5
    allocate ( ratios2(7575, 3495) )
    call LoadFromHDF5DS ( TempFile, &
      & "Standardized_Brightness_Temps_Matrix", &
      & ratios2 )
    if ( myRecalculate ) then
      call output ( 'Recalculating the standardized brightness Temps in Franks file', &
        & advance='yes' )
      ! Recalculate ratios using Frank's own radiances, means, and stddevs
      ! Then check against the ratios he wrote out to TempFile
      allocate ( recalc(7575) )
      allocate ( values2(7575, 3495) )
      call LoadFromHDF5DS ( TempFile, &
        & "Brightness_Temps_Matrix", &
        & values2 )
      do MAF = MAFRange(1), min(MAFRange(2), size(values2, 2) - 1)
        do j=1, size(mean)
          recalc(j) = (values2(j, MAF+1) - mean(j)) / stddev(j)
        enddo
        diffs = recalc - ratios2(:,MAF+1)
        call ShowDiffs ( MAF+1 )
      enddo
      return
    endif
    allocate ( values2(7575, 3495) )
    if ( standardized ) then
      call LoadFromHDF5DS ( TempFile, &
        & "Standardized_Brightness_Temps_Matrix", &
        & values2 )
      radiancetype = 'standardized'
      DEEBug = .true. ! .false.
    else
      call LoadFromHDF5DS ( TempFile, &
        & "Brightness_Temps_Matrix", &
        & values2 )
      radiancetype = 'original'
      DEEBug = .true.
    endif
!     call output( '     ------------------', advance='yes' )
!     call output( '   (Full matrix for last MAF)', advance='yes' )
!     call Dump ( diffs, 'diffs' )
    if ( .not. DEEBug ) return
    do MAF = MAFRange(1), min(MAFRange(2), size(values2, 2) - 1)
      call ShowDiffs ( MAF )
    enddo
    call output( '     ------------------', advance='yes' )
    ! MAF = FindFirst ( abs(Rads(1)-values2(1,:)) < 1.e-4     )
    MAF = FindFirst ( (Rads(1)-values2(1,:))**2 + &
      & (Rads(2)-values2(2,:))**2 < 1.e-8     ) - 1 ! Remember--MAFs start at 0
    if ( present(MatchingMAF) ) MatchingMAF = MAF
    if ( MAF < 0 ) then
      call output( 'No matching MAF found', advance='yes' )
      call outputNamedValue( 'Closest we come is', &
        & minval(abs(Rads(1)-values2(1,:))) )
      call outputNamedValue( 'at MAF', &
        & minloc(abs(Rads(1)-values2(1,:))) )
      return
    else
      call outputNamedValue( 'Matching MAF in Franks file', MAF )
      call output ( (/ real(rads(1), r4), values2(1,MAF+1) /), advance='yes' )
      call output( '     ------------------', advance='yes' )
      call ShowDiffs ( MAF )
      ! Could the MIF indexes be jumbled?
      MIF = FindFirst ( abs(Rads(2)-values2(:,MAF+1)) < 1.e-4     )
      if ( MIF < 1 ) then
        call output( 'No matching MIF found for our 2', advance='yes' )
      else
        call outputNamedValue( 'Matching MIF in Franks file for MIF=2', MIF )
        call output ( (/ real(rads(2), r4), values2(MIF,MAF+1) /), advance='yes' )
        call output( '     ------------------', advance='yes' )
      endif
    endif
    ! Are we saving the matched valus for comparisons and debugging?
    matchedMAF = MAF
    matchedStdRadiances = ratios2(:,MAF+1)
    call OutputNamedValue ( 'min matched rads', minval(matchedStdRadiances) )
    call OutputNamedValue ( 'max matched rads', maxval(matchedStdRadiances) )
    if ( MAF+1 > size(values2,2) ) return
    ! Check for a 2nd matching MAF
    MAF = FindFirst ( abs(Rads(1)-values2(1,MAF+2:)) < 1.e-4     )
    if ( MAF > 0 ) then
      call outputNamedValue( '2nd Matching MAF in Franks file', MAF )
    else
      call output( 'No 2nd matching MAF found', advance='yes' )
    endif    
    contains
      subroutine ShowDiffs ( MAF )
        integer, intent(in)          :: MAF
        ! Executable
        if ( .not. present(Recalculate) ) then
          call output( '     ------------------', advance='yes' )
          call outputNamedValue ( 'MAF used for comparison', MAF )
          call output( '     (As we compute them)', advance='yes' )
          call outputNamedValue ( trim(radianceType) // 'Rad min', minval(Rads) )
          call outputNamedValue ( trim(radianceType) // 'Rad max', maxval(Rads) )
          if ( MAF+1 > size(values2,2) ) then
            call outputNamedValue ( 'MAF too big', MAF )
            call output( '     ------------------', advance='yes' )
            return
          else
          call output( '     (As read from Franks file)', advance='yes' )
          call outputNamedValue ( trim(radianceType) // 'Rad min', minval(values2(:, MAF+1)) )
          call outputNamedValue ( trim(radianceType) // 'Rad max', maxval(values2(:, MAF+1)) )
          diffs = Rads - values2(:, MAF+1)
          endif
        endif
        call outputNamedValue ( 'min diff', minval(diffs) )
        call outputNamedValue ( 'max diff', maxval(diffs) )
        call Dump ( diffs, 'diffs', options='@' )
        call output( '     ------------------', advance='yes' )
      end subroutine ShowDiffs
      
  end subroutine CompareStandardizedRads

  subroutine Dump_Coeffs ( Coeffs, Details )

    use  HighOutput, only: OutputNamedValue
    type ( NeuralNetCoeffs_T )                 :: Coeffs
    integer, intent(in), optional              :: Details
    !                                          
    integer                                    :: myDetails
    ! Executable
    myDetails = 0
    if ( present(Details) ) myDetails = Details
    
    call outputNamedValue ( 'Activation Function', trim(Coeffs%Activation_Function) )
    ! call outputNamedValue ( 'Number of Bands', size(Coeffs%Bands) )
    ! if ( myDetails > 0 ) then
    ! call dump ( Coeffs%Bands, 'Bands' )
    ! call dump ( Coeffs%Channels_In_Each_Band, 'Channels' )
    ! endif

    ! call outputNamedValue ( 'Number of MIFs', size(Coeffs%MIFs) )
    ! if ( myDetails > 0 ) &
    ! call dump ( Coeffs%MIFs, 'MIFs' )
    
    call outputNamedValue ( 'Number of pressure levels', &
         & SIZE(Coeffs%Intercepts_Hidden_Labels_Layer) )
    IF ( myDetails > 0 ) THEN
      !call dump ( Coeffs%Output_Pressure_Levels, 'Output_Pressure_Levels' )
      !call dump ( Coeffs%Output_Pressure_Levels_Indices, 'Output_Pressure_Indices' )
      call dump ( Coeffs%Intercepts_Hidden_Labels_Layer, 'Intercepts layer' )
      call dump ( Coeffs%Normalization_Labels_Max, 'Normalization Labels Max' )
      call dump ( Coeffs%Normalization_Labels_Min, 'Normalization Labels Max' )
    ENDIF

    call outputNamedValue ( 'Number of neurons', size(Coeffs%Intercepts_Hidden_Layer_1) )
    if ( myDetails < 0 ) return    
    if ( myDetails > 0 ) then
    call dump ( Coeffs%Intercepts_Hidden_Layer_1, 'Intercepts layer 1' )
    
    call dump ( Coeffs%Intercepts_Hidden_Layer_2, 'Intercepts layer 2' )
    endif
    
    call outputNamedValue ( 'shape(WHLL)', shape(Coeffs%Weights_Hidden_Labels_Layer ) )
    call outputNamedValue ( 'shape(WHL1)', shape(Coeffs%Weights_Hidden_Layer_1      ) )
    call outputNamedValue ( 'shape(WHL2)', shape(Coeffs%Weights_Hidden_Layer_2      ) )
    if ( myDetails < 1 ) return
    call dump ( Coeffs%Intercepts_Hidden_Layer_1, 'Intercepts layer 1' )
    
    call dump ( Coeffs%Weights_Hidden_Labels_Layer, 'Weights_Hidden_Labels_Layer' )
    call dump ( Coeffs%Weights_Hidden_Layer_1     , 'Weights_Hidden_Layer_1     ' )
    call dump ( Coeffs%Weights_Hidden_Layer_2     , 'Weights_Hidden_Layer_2     ' )
    

    call dump ( Coeffs%Standardization_Brightness_Temperature_Mean, 'Std_Bright_Temp_Mean' )
    call dump ( Coeffs%Standardization_Brightness_Temperature_Std , 'Std_Bright_Temp_Std ' )
  end subroutine Dump_Coeffs

  subroutine NeuralNetFit( nnInputData, nnCoeffs, nHL, prediction, precision, &
    & RadFileID, TempFileID, MAF, profile, Debugging, StdRadiances, &
    & NNOutputData )

    ! Fortran version of the IDL code in `example_idl.pro' from Frank

    type (NeuralNetInputData_T),intent(in):: nnInputData
    type (NeuralNetCoeffs_T), intent(in)  :: nnCoeffs
    integer,intent(in) :: nHL ! number of hidden layers. Will probably always be 2
    ! These will already have been allocated by the caller
    real(r8),  dimension(:) :: prediction
    real(r8),  dimension(:) :: precision
    integer, optional, intent(in) :: RadFileID ! For optionally saving Datasets
    integer, optional, intent(in) :: TempFileID ! For optional comparisons
    integer, optional, intent(in) :: MAF ! For optional comparisons
    integer, optional, intent(in) :: profile ! For optional comparisons
    logical, optional, intent(in) :: Debugging
    real(r8), optional,  dimension(:), intent(in) :: StdRadiances
    type(NeuralNetOutputData_T), optional :: NNOutputData



    ! Local variables
    logical,save :: first = .TRUE.
    integer, save:: nMIFs
    integer, save,  dimension(2) :: nChans
    integer, save :: nSurfs, nVars, nNeurons
    character(len=8) :: type ! One of 'tanh', 'relu' or 'sigmoid'

    ! Various working space arrays
    real(r8), dimension(:), allocatable :: working_space
    real(r8), dimension(:), allocatable :: NeuronValues, NeuronValues2

    real(r8) :: y_pred2

    integer :: c, n, m, jj, s, v ! various counters
    integer :: istat
    integer :: status
    integer, dimension(2) :: dims2
    logical :: debug
    real(r8) :: spread

    ! ----------- executable statements -----------------
     type = nnCoeffs%Activation_Function
     status = 0
     debug = .false.
     IF (PRESENT(debugging)) debug=debugging
     IF (debug) THEN 
       ! open an output file (for debugging purposes)
       OPEN(unit=7,file='intermediate_results.dat',status='replace',&
            & form='unformatted', iostat=istat, &
            & access='stream')
       IF (istat /= 0) THEN 
         print *,'bad open! istat = ',istat
         return
       ENDIF
     ENDIF

    ! make sure type is correct.
    IF ( index( 'relu,sigmoid,tanh', trim(type) ) < 1 ) then 
      print *,"input var type must be one of 'relu', 'tanh' or 'sigmoid'"
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Illegal activation type ' // trim(type) )
    ENDIF
    IF (NHL /=  1 .AND. NHL /= 2) THEN 
      print *,"NHL must equal either 1 or 2"
     ENDIF


    ! These are fixed for the run
    IF (FIRST) THEN 
      nMIFs = nnInputData%Band_1_Radiances%numMIFs
      nChans(1)=nnInputData%Band_1_Radiances%numChannels
      nChans(2)=nnInputData%Band_22_Radiances%numChannels
      nSurfs = SIZE(nnCoeffs%Intercepts_HIdden_Labels_Layer)
      dims2= size( nnCoeffs%Intercepts_HIdden_Layer_1 )
      nNeurons = dims2(1)
      dims2=SIZE( nnCoeffs%Standardization_Brightness_Temperature_Mean ) 
      nVars = dims2(1)
      first=.FALSE.
    ENDIF


    IF (debug) THEN 
      print *,'nMIFs=',nMIFs
      print *,'nChans=',nChans
      print *,'nSurfs=',nSurfs
      print *,'nNeurons=',nNeurons
      print *,'nVars=',nVars
    ENDIF

    ! allocate space                        
    ! nVars = 2*Chanels(bands 1&8) *  nMIFs + nChans(22)*nMIFs
    ! = 2*25*75 + 51*75  = 7575
    ! nNeurons=5078
    ! nLevs=42
    allocate( working_space(nVars )) 
    allocate(neuronValues(nNeurons))
    ! These will have already have been allocated in the caller
!     allocate(prediction(nSurfs))
!     allocate(precision(nSurfs))

    allocate(neuronValues2(nNeurons))

    ! <Paul> See note above about this output user type.
    if ( present(nnOutputData) ) then
      allocate(nnOutputData%working_space( nVars )) 
      allocate(nnOutputData%neuronValues(nNeurons))
      allocate(nnOutputData%neuronValues2(nNeurons))
      allocate(nnOutputData%prediction(nSurfs))
      allocate(nnOutputData%y_pred(nSurfs)) 
    endif


    neuronValues = 0.0
    prediction   = 0.0
    precision    = 0.0

    ! load the data into the working_space
    jj = 1
    ! Load band1
    DO m=1,nMIFs
      DO c=1,nChans(1)
        working_space( jj ) = nnInputData%Band_1_Radiances%Values(c,m)
        jj = jj + 1
      end DO
    end DO

    ! Load band8
    DO m=1,nMIFs
      DO c=1,nChans(1)
        working_space( jj ) = nnInputData%Band_8_Radiances%Values(c,m)
        jj = jj + 1
      end DO
    end DO

    ! load band 22
    DO m=1,nMIFs
      DO c=1,nChans(2)
        working_space( jj ) = nnInputData%Band_22_Radiances%Values(c,m)
        jj = jj + 1
      end DO
    end DO

    IF (debug) THEN 
      print *,'before normalization: size(working_space) = ',SIZE(working_space)
      WRITE (7,iostat=istat) working_space

    ! ============ normalize the brightness temperatures.  ===========

    ! <Paul> I don't recall whether this is done when running in the L2
    ! environment. If not, uncomment these lines
    ! working_space = (working_space - &
    !      & nnCoeffs%Standardization_Brightness_Temperature_Mean ) / &
    !      & nnCoeffs%Standardization_Brightness_Temperature_Std

    call OutputNamedValue ( 'nVars', nVars )
    call OutputNamedValue ( 'shape(WHL1)', &
         SHAPE(nnCoeffs%Weights_Hidden_Layer_1) )

      print *,'sbtm, sbts'
      write (7,iostat=istat) nnCoeffs%Standardization_Brightness_Temperature_Mean 
      write (7,iostat=istat) nnCoeffs%Standardization_Brightness_Temperature_Std
      print *,'after normalization: size(working_space) e= ',SIZE(working_space)
      write (7,iostat=istat) working_space

      print *,'whl1 [nNeurons, nVars], whl2 [nNeurons, nNeurons]'
      print *,'size(Weights_Hidden_Layer_1)=',SIZE(nnCoeffs%Weights_Hidden_Layer_1)
      print *,'size(Weights_Hidden_Layer_2) = ',SIZE(nnCoeffs%Weights_Hidden_Layer_2)
      write (7,iostat=istat) nnCoeffs%Weights_Hidden_Layer_1, & 
           & nnCoeffs%Weights_Hidden_Layer_2
    ENDIF

    neuronValues=0.0
    ! working_space is [nVars]
    ! weights_hl_1 is [nNeurons,nVars]
    ! weights_hl_2 is [nNeurons,nNeurons]
    ! neuronValues is [nNeurons]
    ! Weights_Hidden_layer_1 is [numNeurons, numVars]
    ! Weights_Hidden_layer_2 is [numNeurons,numNeurons ]
    DO v=1,nVars
      neuronValues = neuronValues + &
           & working_space(v) * nnCoeffs%Weights_Hidden_Layer_1(v,:) 
!           & working_space(v) * nnCoeffs%Weights_Hidden_Layer_1(:,v) 
    end DO
    IF (debug) THEN 
      print *,'after calculation with whl1: size(neuronValues) = ',SIZE(neuronValues)
      WRITE(7,iostat=istat) neuronValues
    ENDIF
    !print *,'k = ',k ! 8
    ! k=k+1
    neuronValues = Activation ( type, neuronValues + & 
           & nnCoeffs%Intercepts_Hidden_Layer_1 )
    IF (debug) THEN 
      print *,'After ihl1'
      print*,'intercepts_hidden_layer_[12]',size(nnCoeffs%Intercepts_Hidden_Layer_1)
      WRITE (7) nnCoeffs%Intercepts_Hidden_Layer_1, nnCoeffs%Intercepts_Hidden_Layer_2
      print *,'ihll '
      WRITE (7) nnCoeffs%Intercepts_Hidden_Labels_Layer
      print *,'whll'
      WRITE (7) nnCoeffs%Weights_Hidden_Labels_Layer
      print *,'after ihl1: size(neuronValues) = ',SIZE(neuronValues)
      WRITE(7,iostat=istat) neuronValues
    ENDIF
    if ( debug ) then
      call OutputNamedValue ( 'num neurons', nNeurons )
      call OutputNamedValue ( 'shape(WHLL)', shape(nnCoeffs%Weights_Hidden_Labels_Layer) )
    endif


    !print *,'k = ',k ! 9
    ! k=k+1


    !print *,'nHL = ',nHL
    IF (nHL .EQ. 1) THEN 
      ! 1 hidden layer
      ! Loop over surfaces
      
      DO s=1,nSurfs 

!        y_pred2 = DOT( neuronValues, &
        y_pred2 = dot_product( neuronValues, &
             & nnCoeffs%Weights_Hidden_Labels_Layer(:,s)  ) + &
             &   nnCoeffs%Intercepts_Hidden_Labels_Layer(s)

        if ( present(nnOutputData) ) then
        nnOutputData%y_pred(s)=y_pred2

        ! Denormalize the 'labels'
        nnOutputData%prediction(s) = y_pred2 * &
             & ( nnCoeffs%Normalization_Labels_Max(s) - &
             &   nnCoeffs%Normalization_Labels_Min(s)) + &
             &   nnCoeffs%Normalization_Labels_Min(s)
        endif
        prediction(s) = y_pred2 * &
             & ( nnCoeffs%Normalization_Labels_Max(s) - &
             &   nnCoeffs%Normalization_Labels_Min(s)) + &
             &   nnCoeffs%Normalization_Labels_Min(s)

        !print '(a14,",",i2,",",f0.2,",",f0.2)','s, yp2,pred = ',&
        !     &              s,y_pred2,nnOutputData%prediction(s)
      ENDDO ! loop over pressure surfaces

    ELSEIF (nHL .EQ. 2) THEN 


      ! 2 hidden layers
      neuronValues2=0.0
      DO n=1,nNeurons
        neuronValues2 = neuronValues2 + &
             & neuronValues(n) * nnCoeffs%Weights_Hidden_Layer_2(n,:)
      end DO

      IF (debug) THEN 
        print *,'Writing size(neuronValues2) = ',SIZE(neuronValues2)
        WRITE(7,iostat=istat) neuronValues2
      ENDIF
      neuronValues2 = Activation ( type, neuronValues2 + & 
           & nnCoeffs%Intercepts_Hidden_Layer_1 )
      IF (debug) THEN 
        print *,'after applying ihl2: Writing size(neuronValues2) = ',SIZE(neuronValues2)
        WRITE(7,iostat=istat) neuronValues2
      ENDIF

      ! Calculate the final product
      DO s=1,nSurfs 

!        y_pred2 = DOT( neuronValues2, &
        y_pred2 = dot_product( neuronValues2, &
             & nnCoeffs%Weights_Hidden_Labels_Layer(:,s)  ) + &
             &   nnCoeffs%Intercepts_Hidden_Labels_Layer(s)

        if ( present(nnOutputData) ) then
        nnOutputData%y_pred(s)=y_pred2

        ! Denormalize the 'labels'
        nnOutputData%prediction(s) = y_pred2 * &
             & ( nnCoeffs%Normalization_Labels_Max(s) - &
             &   nnCoeffs%Normalization_Labels_Min(s)) + &
             &   nnCoeffs%Normalization_Labels_Min(s)
        endif
        prediction(s) = y_pred2 * &
             & ( nnCoeffs%Normalization_Labels_Max(s) - &
             &   nnCoeffs%Normalization_Labels_Min(s)) + &
             &   nnCoeffs%Normalization_Labels_Min(s)

        !print '(a14,",",i2,",",f0.2,",",f0.2)','s, yp2,pred = ',&
        !     &              s,y_pred2,nnOutputData%prediction(s)
      ENDDO ! loop over pressure surfaces
      IF (debug) THEN 
        print *,'Normalization labels min/max'
        write (7,iostat=istat) nnCoeffs%Normalization_Labels_Min, &
             & nnCoeffs%Normalization_Labels_Max
        if ( present(nnOutputData) ) then
          print *,'size(nnOutputData%y_pred) = ',SIZE(nnOutputData%y_pred)
          write(7,iostat=istat) nnOutputData%y_pred
          print *,'size(nnOutputData%prediction) = ',SIZE(nnOutputData%prediction)
          write(7,iostat=istat) nnOutputData%prediction
        endif
      endif

    endif ! NHl .eq. 2

    ! Check that the prediction is within the range
    !    [min-spread[s], max+spread[s]]
    ! where spread[s] = max - min
    ! If any level strays outisde the range, set all predictions to -999.99
    do s=1,nSurfs 
      spread = &
        & nnCoeffs%Normalization_Labels_Max(s) &
        & - &
        & nnCoeffs%Normalization_Labels_Min(s)
        if ( &
          & (prediction(s) < nnCoeffs%Normalization_Labels_Min(s) - spread) &
          & .or. & 
          & (prediction(s) > nnCoeffs%Normalization_Labels_Max(s) + spread) &
          & ) &
          & precision = UndefinedValue
    enddo




    if ( present(nnOutputData) ) then
      nnOutputData%neuronValues  = neuronValues
      nnOutputData%neuronValues2 = neuronValues2
      nnOutputData%working_space = working_space
    endif

    status=1

    !print *,'k = ',k ! 10
    ! k=k+1

    if (debug) CLOSE(7)
    contains
      function Activation ( type, args ) result ( values )
        ! Apply the activation function type to the args
        ! returning the result in values
        character(len=*), intent(in)           :: type
        real(r8), dimension(:), intent(in)     :: args
        real(r8), dimension(size(args))        :: values
        ! Executable
        select case ( trim(type) )
        case ( 'tanh' )  
          ! In Frank's code he calls this `neuronValuesActivation', but I'm
          ! just going to reuse  the variable.
          print *, '*** tanh *** '
          Values = TANH( args )
        case ( 'sigmoid' )
          print *, '*** sigmoid *** '
          call sigmoid( args, Values )

        case ( 'relu' )
          print *, '*** relu *** '
          Values = max( 0._r8, args )

        end select
      end function Activation

  end subroutine NeuralNetFit

  subroutine StandardizeRadiances ( nnInputData, &
              & Mean, &
              & StdDev, &
              & TempFileID, MAF, nnCoeffs, Debugging )
    ! Compute the ratio
    !       values - mean
    !       --------------
    !           StdDev
    ! and use  it to replace corresponding values in NNMeasurements.
    type (NeuralNetInputData_T),intent(inout)              :: nnInputData
    real(r8),  dimension(:), intent(in)                    :: Mean
    real(r8),  dimension(:), intent(in)                    :: StdDev
    integer, optional, intent(in) :: TempFileID ! For optional comparisons
    integer, optional, intent(in) :: MAF ! For optional comparisons; starts at 0
    type (NeuralNetCoeffs_T), optional, intent(in)        :: nnCoeffs
    logical, optional, intent(in)       :: Debugging
    ! Internal variables
    real(r8),  dimension(:), allocatable                   :: ratio
    real(r8),  dimension(:), allocatable                   :: values
    integer                                               :: channel
    logical                                               :: DeeBug
    integer                                               :: j ! index into mean
    integer                                               :: MatchingMAF
    integer                                               :: MIF
    integer                                               :: n ! how many overall
    ! Executable
    DeeBug = .false.
    if ( present(Debugging) ) DeeBug = Debugging
    if ( size(mean) /= size(stddev) ) then
         CALL MLSMessage( MLSMSG_error, ModuleName, &
             & "input arrays mean and stddev must have the same size" )
    endif
    allocate( values(size(mean) ) )
    allocate( ratio(size(mean) ) )
    ! Gather values
    j = 0
    ! Band _1_
    do MIF=1, nnInputData%Band_1_Radiances%NumMIFs
      do channel=1, nnInputData%Band_1_Radiances%NumChannels
        j = j + 1
        values(j) = nnInputData%Band_1_Radiances%values(channel, MIF)
      enddo
    enddo
    ! Band _8_
    do MIF=1, nnInputData%Band_8_Radiances%NumMIFs
      do channel=1, nnInputData%Band_8_Radiances%NumChannels
        j = j + 1
        values(j) = nnInputData%Band_8_Radiances%values(channel, MIF)
      enddo
    enddo
    ! Band _22_
    do MIF=1, nnInputData%Band_22_Radiances%NumMIFs
      do channel=1, nnInputData%Band_22_Radiances%NumChannels
        j = j + 1
        values(j) = nnInputData%Band_22_Radiances%values(channel, MIF)
      enddo
    enddo
    n = j
    if ( n /= size(mean) ) then
         CALL MLSMessage(MLSMSG_error, ModuleName, &
             & "n must equal size of mean array")
    endif
    if ( DeeBug ) print *, 'n: ', n
    do j=1, n
      ratio(j) = ( values(j) - mean(j) ) / stddev(j)
    enddo
    if ( .false. ) then
      call Dump( values, 'values' )
      call Dump( mean ,  'mean  ' )
      call Dump( stddev, 'stddev' )
      call Dump( ratio , 'ratio ' )
    endif
    if ( DeeBug ) then
      call OutputNamedValue ( '1st(values)', values(1) )
      call OutputNamedValue ( 'max(values)', maxval(values) )
      call OutputNamedValue ( '1st(mean  )', mean(1) )
      call OutputNamedValue ( 'max(mean  )', maxval(mean) )
      call OutputNamedValue ( '1st(stddev)', stddev(1) )
      call OutputNamedValue ( 'max(stddev)', maxval(stddev) )
      call OutputNamedValue ( '1st(Ratio )', ratio(1) )
      call OutputNamedValue ( 'min(Ratio)', minval(ratio) )
      call OutputNamedValue ( 'max(Ratio)', maxval(ratio) )
    endif
    if ( present(TempFileID) ) then
      call CompareStandardizedRads ( TempFileID, &
      & (/ MAF, MAF /), values, &
      & standardized=.false., MatchingMAF=MatchingMAF )
      call OutputNamedValue ( 'Matching MAF', MatchingMAF )
      call OutputNamedValue ( 'Compared to', MAF )
    endif
    ! Now use  these standardized radiances to replace the measurement's values
    j = 0
    ! Band _1_
    do MIF=1, nnInputData%Band_1_Radiances%NumMIFs
      do channel=1, nnInputData%Band_1_Radiances%NumChannels
        j = j + 1
        nnInputData%Band_1_Radiances%values(channel, MIF) = ratio(j)
      enddo
    enddo
!     print *, 'j: ', j
!     print *, ratio(j)
    ! Band _8_
    do MIF=1, nnInputData%Band_8_Radiances%NumMIFs
      do channel=1, nnInputData%Band_8_Radiances%NumChannels
        j = j + 1
        nnInputData%Band_8_Radiances%values(channel, MIF) = ratio(j)
      enddo
    enddo
!     print *, 'j: ', j
!     print *, ratio(j)
    ! Band _22_
    do MIF=1, nnInputData%Band_22_Radiances%NumMIFs
      do channel=1, nnInputData%Band_22_Radiances%NumChannels
        j = j + 1
        nnInputData%Band_22_Radiances%values(channel, MIF) = ratio(j)
      enddo
    enddo
!     print *, 'j: ', j
!     print *, ratio(j)
    !
    if ( .not. present(TempFileID) ) return
    call CompareStandardizedRads ( TempFileID, &
    & (/ MAF, MAF /), values, &
    & standardized=.true., recalculate=.true., mean=mean, stddev=stddev )
    
    ! What if the binNum we used was different from Frank's?
    if ( present(nnCoeffs) ) &
      & call CheckBinNums ( TempFileID, (/ MAF, MAF /), nnCoeffs )
    
  end subroutine StandardizeRadiances

  ! Dot product.
  ! (Not used)
  real(r8) function dot(X,Y)
    real(r8), dimension(:), intent(in):: x,y
    integer :: i,nn
    dot=0
    nn=size(x)
    IF (nn .NE. SIZE(y)) THEN 
      CALL MLSMessage(MLSMSG_error, ModuleName, &
           & 'The two vectors are not the same size!')
    ENDIF

    DO i=1,SIZE(x)
      dot = dot + x(i)*y(i)
    end DO
  end function dot

  ! Sigmoid subroutine
  subroutine sigmoid(x,s) 
    real(r8), dimension(:),intent(in) :: x
    real(r8), dimension(:) :: s
    s = 1.0/(1.0+EXP(-x))
  end subroutine sigmoid

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if print is added
  end function not_used_here
!---------------------------------------------------------------------------

end module NeuralNetUtils_M
! $Log$
! Revision 2.10  2021/07/28 23:41:37  pwagner
! Fix bug preventing full use of relu activation function
!
! Revision 2.9  2021/07/08 23:29:11  pwagner
! New api for NeuralNetFit; new 'relu' type; housekeeping
!
! Revision 2.8  2021/06/18 15:20:33  pwagner
! Refine search for matching bins; avoid array bounds error
!
! Revision 2.7  2021/05/27 23:41:00  pwagner
! Corrected errors traced to MAF index starting at 0; removed undefined variable 'k'
!
! Revision 2.6  2021/05/18 15:52:46  pwagner
! Many bugs fixed
!
! Revision 2.5  2021/04/01 23:46:21  pwagner
! many more (too many?) debugging options
!
! Revision 2.4  2021/03/18 23:47:41  pwagner
! Fixed some more errors; added more debugging aids
!
! Revision 2.3  2021/03/05 00:53:34  pwagner
! Some progress but still wrong
!
! Revision 2.2  2021/02/19 00:29:46  pwagner
! repaired many bugs; still unsatisfactory imo
!
! Revision 2.1  2021/02/05 05:14:40  pwagner
! First commit
!
