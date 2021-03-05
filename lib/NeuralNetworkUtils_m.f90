! Copyright 2021, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.
!
! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

! This module contains the utility code to be used when calculating
! species using Frank Werner's Neural Network solutions. At the time
! of this writing (2021/01/20) only Temperature is being retrieved,
! but the code should be able to do any species, provided the correct
! training data is provided. 

! To adapt to other species, we must revisit the hardcoded use of some
! signal names, arrays sizes, etc. both here and in l2/NeuralNet_m.f90

MODULE NeuralNetUtils_M

  use Dump_0, only: Dump
  use HighOutput, only: OutputNamedValue
  use MLSCommon, only: MLSFile_T
  use MLSHDF5, only: SaveAsHDF5DS
  use MLSKinds, only: R4, R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use Output_M, only: Output
  implicit none

  TYPE radiance_T
     character(len=64) :: signal
     !INTEGER :: NumMAFS ! <----- it's only ever going to be one MAF, right?
     !
     ! Instead of the original quantity's larger number of channels and MIFs,
     ! we store only the useful ones, by downsampling using integer
     ! arrays in the NeuralNetCoeffs (see below)
     integer :: NumChannels
     integer :: NumMIFs
     real(R8), dimension(:,:), allocatable :: values
  END TYPE radiance_T

  TYPE NeuralNetCoeffs_T
     ! Will be allocated in the caller. 
     ! [3]  Number of Bands
     character(len=32), dimension(:), allocatable :: Bands
     
     ! The first index is M, the max num of channels used (probably 51)
     ! The second index is the number of Bands (currently 3)
     ! We need this info because the Neural Net doesn't currently use
     ! every channel in the DACS 
     ! (and someday may even skip some in the non-DACS)
     ! So (M,3)
     integer,  dimension(:,:), allocatable :: Channels_In_Each_Band

     ! these are the dimensions for each part (assuming that we keep
     ! with Frank set up for Temperature retrieval)

     ! [42], number of pressure levels
     real(R8), dimension(:), allocatable :: Intercepts_Hidden_Labels_Layer
     real(R8), dimension(:), allocatable :: Output_Pressure_Levels
     integer,  dimension(:), allocatable :: Output_Pressure_Levels_Indices

     ! [5078]. Number of Neurons in hidden layers
     real(R8), dimension(:), allocatable :: Intercepts_Hidden_Layer_1
     real(R8), dimension(:), allocatable :: Intercepts_Hidden_Layer_2

     ! [5078, 42] = [nNeurons, nLevels]
     real(R8), dimension(:,:), allocatable :: Weights_Hidden_Labels_Layer
     real(R8), dimension(:,:), allocatable :: Weights_Hidden_Layer_1
     real(R8), dimension(:,:), allocatable :: Weights_Hidden_Layer_2

     ! [42]. Number of levels
     real(R8), dimension(:), allocatable :: Normalization_Labels_Max
     real(R8), dimension(:), allocatable :: Normalization_Labels_Min
     
     ! [75]. Number of MIFs
     integer, dimension(:), allocatable :: MIFs

     ! [7575] number of 'variables', i.e. 2*25*75 + 51*75
     ! This is band 1[channel X MIF] + band 8 [channel X MIF ]  and 
     ! band 22 [channel X MIF]
     real(R8), dimension(:), allocatable :: Standardization_Brightness_Temperature_Mean
     real(r8), dimension(:), allocatable :: Standardization_Brightness_Temperature_Std

  END TYPE NeuralNetCoeffs_T

  TYPE NeuralNetInputData_T
     ! I'm assuming that will have been corrected for the baseline and
     ! reduced to the appropriate MIFs, 22-96 for all bands, all 25
     ! channels for Bands 1 and 8 and channels 40-90 for Band
     ! 22. 
     !
     ! All of these will end up having 
     ! These will be allocated in the caller.

     TYPE (Radiance_T) :: Band_1_Radiances 
     TYPE (Radiance_T) :: Band_8_Radiances
     TYPE (Radiance_T) :: Band_22_Radiances 

  END TYPE NeuralNetInputData_T


  interface Dump
    module procedure Dump_Coeffs !, Dump_Radiance, Dump_InputData
  end interface

  PUBLIC :: Dump, &
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


CONTAINS 
  subroutine CompareStandardizedRads ( TempFileID, MAF, Rads, standardized )
    ! Compare with the standardized rads written to the TempFile by
    ! Frank's python script
    use Intrinsic, only: L_HDF
    use MLSHDF5, only: LoadFromHDF5DS
    use MLSFiles, only: HDFVersion_5
    ! Dummy args
    integer, intent(in) :: TempFileID ! For optional comparisons
    integer, intent(in) :: MAF ! For optional comparisons
    real(r8), intent(in), dimension(:) :: Rads
    logical, intent(in) :: standardized
    !
    type (MLSFile_T)                       :: TempFile
    real(r4), allocatable, dimension(:,:)  :: values2
    real(r4), dimension(size(Rads))        :: diffs
    character(len=32)                      :: radianceType
    ! Executable
    TempFile%FileID%f_id = TempFileID
    TempFile%StillOpen = .true.
    TempFile%type = l_hdf
    TempFile%hdfVersion = HDFVersion_5
    allocate ( values2(7575, 3495) )
    if ( standardized ) then
      call LoadFromHDF5DS ( TempFile, &
        & "Standardized_Brightness_Temps_Matrix", &
        & values2 )
      radianceType = 'standardized'
    else
      call LoadFromHDF5DS ( TempFile, &
        & "Brightness_Temps_Matrix", &
        & values2 )
      radianceType = 'original'
    endif
    diffs = Rads - values2(:, MAF)
    call output( '     ------------------', advance='yes' )
    call output( '     (As we compute them)', advance='yes' )
    call outputNamedValue ( trim(radianceType) // 'Rad min', minval(Rads) )
    call outputNamedValue ( trim(radianceType) // 'Rad max', maxval(Rads) )
    call output( '     (As read from Franks file)', advance='yes' )
    call outputNamedValue ( trim(radianceType) // 'Rad min', minval(values2(:, MAF)) )
    call outputNamedValue ( trim(radianceType) // 'Rad max', maxval(values2(:, MAF)) )
    call outputNamedValue ( 'min diff', minval(diffs) )
    call outputNamedValue ( 'max diff', maxval(diffs) )
    call Dump ( diffs, 'diffs', options='@' )
    call output( '     ------------------', advance='yes' )
  end subroutine CompareStandardizedRads

  subroutine Dump_Coeffs ( Coeffs, Details )
    type ( NeuralNetCoeffs_T )                 :: Coeffs
    integer, intent(in), optional              :: Details
    !                                          
    integer                                    :: myDetails
    ! Executable
    myDetails = 0
    if ( present(Details) ) myDetails = Details
    
    call outputNamedValue ( 'Number of Bands', size(Coeffs%Bands) )
    if ( myDetails > 0 ) then
    call dump ( Coeffs%Bands, 'Bands' )
    call dump ( Coeffs%Channels_In_Each_Band, 'Channels' )
    endif

    call outputNamedValue ( 'Number of MIFs', size(Coeffs%MIFs) )
    if ( myDetails > 0 ) &
    call dump ( Coeffs%MIFs, 'MIFs' )
    
    call outputNamedValue ( 'Number of pressure levels', &
         & SIZE(Coeffs%Intercepts_Hidden_Labels_Layer) )
    if ( myDetails > 0 ) then
    call dump ( Coeffs%Output_Pressure_Levels, 'Output_Pressure_Levels' )
    call dump ( Coeffs%Output_Pressure_Levels_Indices, 'Output_Pressure_Indices' )
    call dump ( Coeffs%Intercepts_Hidden_Labels_Layer, 'Intercepts layer' )
    call dump ( Coeffs%Normalization_Labels_Max, 'Normalization Labels Max' )
    call dump ( Coeffs%Normalization_Labels_Min, 'Normalization Labels Max' )
    endif
    
    call outputNamedValue ( 'Number of neurons', size(Coeffs%Intercepts_Hidden_Layer_1) )
    
    if ( myDetails > 0 ) then
    call dump ( Coeffs%Intercepts_Hidden_Layer_1, 'Intercepts layer 1' )
    
    call dump ( Coeffs%Intercepts_Hidden_Layer_2, 'Intercepts layer 2' )
    endif
    
    call outputNamedValue ( 'Number of neurons', size(Coeffs%Intercepts_Hidden_Layer_1) )
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

  FUNCTION NeuralNetFit( nnInputData, nnCoeffs, nHL, TYPE, &
    & RadFileID, TempFileID, MAF, profile ) RESULT(prediction)

    ! Fortran version of the IDL code in `example_idl.pro' from Frank

    TYPE (NeuralNetInputData_T),intent(in):: nnInputData
    TYPE (NeuralNetCoeffs_T), intent(in)  :: nnCoeffs
    integer,intent(in) :: nHL ! number of hidden layers. Will probably always be 2
    CHARACTER(len=*),intent(in) :: TYPE ! type  may be either 'tanh' or 'sigmoid'
    real(r8), dimension(:), allocatable :: prediction
    integer, optional, intent(in) :: RadFileID ! For optionally saving Datasets
    integer, optional, intent(in) :: TempFileID ! For optional comparisons
    integer, optional, intent(in) :: MAF ! For optional comparisons
    integer, optional, intent(in) :: profile ! For optional comparisons

    ! Local variables
    logical,save :: first = .TRUE.
    integer, save:: nMIFs
    integer, save, dimension(2) :: nChans
    integer, save :: nSurfs, nVars, nNeurons

    ! Various working space arrays
    real(r8),dimension(:), allocatable :: working_space
    real(r8), dimension(:), allocatable  :: NeuronValues, NeuronValues2

    real(r8) :: y_pred

    integer :: c,n,m,jj,s,v! various counters
    integer, dimension(2) :: dims2


    ! ----------- executable statements -----------------

    ! make sure TYPE is correct.
    IF ( TRIM(TYPE) .NE. 'tanh' .AND. &
         & TRIM(TYPE) .NE. 'sigmoid' ) THEN 
       CALL MLSMessage(MLSMSG_error, ModuleName, &
           & "input var TYPE must equal either 'tanh' or 'sigmoid'")
    ENDIF
    IF (NHL /=  1 .AND. NHL /= 2) THEN 
       CALL MLSMessage(MLSMSG_error, ModuleName, &
           & "NHL must equal either 1 or 2")
     ENDIF


    call Dump ( nnCoeffs )

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

    ! allocate space
    ! nMIFs=75 for all bands (so far)
    ! nVars = 2*Chanels(bands 1&8) *  nMIFs + nChans(22)*nMIFs
    ALLOCATE( working_space(nVars )) 
    ALLOCATE(neuronValues(nNeurons))
    ALLOCATE(neuronValues2(nNeurons))
    ALLOCATE(prediction(nSurfs))

    neuronValues = 0.0
    prediction = 0.0

    ! load the data into the working_space. The data is arranged as
    ! the radiances in MIF, then channel order for band 1, band 8 and
    ! then band 22
    jj=1

    ! Load band1
    DO c=1,nChans(1)
      DO m=1,nMIFs
        working_space( jj ) = nnInputData%Band_1_Radiances%Values(c,m)
        jj = jj + 1
      END DO
    END DO

    ! Load band8
    DO c=1,nChans(1)
      DO m=1,nMIFs
        working_space( jj ) = nnInputData%Band_8_Radiances%Values(c,m)
        jj = jj + 1
      END DO
    END DO

    ! load band 22
    DO c=1,nChans(2)
      DO m=1,nMIFs
        working_space( jj ) = nnInputData%Band_22_Radiances%Values(c,m)
        jj = jj + 1
      END DO
    END DO


    call OutputNamedValue ( 'nVars', nVars )
    call OutputNamedValue ( 'shape(WHL1)', shape(nnCoeffs%Weights_Hidden_Layer_1) )
    call OutputNamedValue ( 'min(Rads)', minval(working_space(1:jj-1)) )
    call OutputNamedValue ( 'max(Rads)', maxval(working_space(1:jj-1)) )
    
    ! Optionally save this dataset for comparison with
    ! PyThon or iDl or whaTeVer you guys use these daYs
    if ( present(RadFileID) ) &
      & call SaveAsHDF5DS ( RadFileID, 'Radiances', working_space(1:jj-1) )

    ! Optionally save this dataset for comparison with
    ! PyThon or iDl or whaTeVer you guys use these daYs
    if ( present(TempFileID) ) &
      & call CompareStandardizedRads ( TempFileID, MAF, working_space(1:jj-1), &
      & standardized=.true. )

    ! unlike the IDL code, this routine will only ever see 1 MAF
    ! (sample) at a time, so we don't need the IDL loop over `samples'
    ! working_space is [nVars]
    ! weights is [nVars,nNeurons]
    ! neuronValues is [nNeurons]
    DO v=1,nVars
      neuronValues = neuronValues + &
           & SUM( working_space(v) * nnCoeffs%Weights_Hidden_Layer_1(v,:))
    END DO


    IF (TRIM(TYPE) .EQ. 'tanh') THEN 
      neuronValues=TANH(neuronValues + & 
           & nnCoeffs%Intercepts_Hidden_Layer_1)
    ELSE IF (TRIM(TYPE) .EQ. 'sigmoid') THEN 
      CALL SIGMOID(neuronValues + & 
           & nnCoeffs%Intercepts_Hidden_Layer_1, neuronValues)

    ENDIF
    
    call OutputNamedValue ( 'num neurons', nNeurons )
    call OutputNamedValue ( 'shape(WHLL)', shape(nnCoeffs%Weights_Hidden_Labels_Layer) )


    neuronValues2=0.0
    ! Loop over surfaces
    DO s=1,nSurfs 
      ! neuron values 'activation' 
      IF (nHL .EQ. 1) THEN 

        IF (TRIM(TYPE) .EQ. 'tanh') THEN
          neuronValues = TANH(neuronValues + &
               & nnCoeffs%Weights_Hidden_Labels_Layer(:,s)) + &
               & nnCoeffs%Intercepts_Hidden_labels_layer(s)
        ELSE IF(trim(TYPE) .EQ. 'sigmoid') THEN 
          CALL SIGMOID(neuronValues + &
               & nnCoeffs%Weights_Hidden_Labels_Layer(:,s), neuronValues)
               neuronValues = neuronValues + nnCoeffs%Intercepts_Hidden_Labels_Layer(s)
        ENDIF

      ELSEIF (nHL .EQ. 2) THEN 

        DO n=1,nNeurons
          neuronValues2 = neuronValues2 + &
               & SUM( neuronValues(n) * nnCoeffs%Weights_Hidden_Layer_2(n,:))
        END DO


        IF (TRIM(TYPE) .EQ. 'tanh') THEN
          neuronValues2 = TANH(neuronValues2 + nnCoeffs%Intercepts_Hidden_Layer_2)
        ELSE IF(trim(TYPE) .EQ. 'sigmoid') THEN 
          call SIGMOID(neuronValues2 + nnCoeffs%Intercepts_Hidden_Layer_2,neuronValues)
        ENDIF
      ENDIF

      y_pred = DOT( neuronValues2, &
           & nnCoeffs%Weights_Hidden_Labels_Layer(:,s) ) + &
           &   nnCoeffs%Intercepts_Hidden_Labels_Layer(s)

      ! Calculate the final product
      prediction(s) = y_pred * &
           & ( nnCoeffs%Normalization_Labels_Max(s) - &
           &   nnCoeffs%Normalization_Labels_Min(s)) + &
           &   nnCoeffs%Normalization_Labels_Min(s)

    ENDDO ! loop over pressure surfaces


  END function NeuralNetFit
  
  subroutine StandardizeRadiances ( nnInputData, &
              & Mean, &
              & StdDev, &
              & TempFileID, MAF )
    ! Compute the ratio
    !       values - mean
    !       --------------
    !           StdDev
    ! and use it to replace corresponding values in NNMeasurements.
    type (NeuralNetInputData_T),intent(inout)             :: nnInputData
    real(r8), dimension(:), intent(in)                    :: Mean
    real(r8), dimension(:), intent(in)                    :: StdDev
    integer, optional, intent(in) :: TempFileID ! For optional comparisons
    integer, optional, intent(in) :: MAF ! For optional comparisons
    ! Internal variables
    real(r8), dimension(:), allocatable                   :: ratio
    real(r8), dimension(:), allocatable                   :: values
    integer                                               :: channel
    integer                                               :: j ! index into mean
    integer                                               :: MIF
    integer                                               :: n ! how many overall
    ! Executable
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
    print *, 'n: ', n
    do j=1, n
      ratio(j) = ( values(j) - mean(j) ) / stddev(j)
    enddo
    call Dump( values, 'values' )
    call Dump( mean ,  'mean  ' )
    call Dump( stddev, 'stddev' )
    call Dump( ratio , 'ratio ' )
    call OutputNamedValue ( 'min(Ratio)', minval(ratio) )
    call OutputNamedValue ( 'max(Ratio)', maxval(ratio) )
    if ( present(TempFileID) ) &
      & call CompareStandardizedRads ( TempFileID, MAF, values, &
      & standardized=.false. )
    ! Now use these standardized radiances to replace the measurement's values
    j = 0
    ! Band _1_
    do MIF=1, nnInputData%Band_1_Radiances%NumMIFs
      do channel=1, nnInputData%Band_1_Radiances%NumChannels
        j = j + 1
        nnInputData%Band_1_Radiances%values(channel, MIF) = ratio(j)
      enddo
    enddo
    print *, 'j: ', j
    print *, ratio(j)
    ! Band _8_
    do MIF=1, nnInputData%Band_8_Radiances%NumMIFs
      do channel=1, nnInputData%Band_8_Radiances%NumChannels
        j = j + 1
        nnInputData%Band_8_Radiances%values(channel, MIF) = ratio(j)
      enddo
    enddo
    print *, 'j: ', j
    print *, ratio(j)
    ! Band _22_
    do MIF=1, nnInputData%Band_22_Radiances%NumMIFs
      do channel=1, nnInputData%Band_22_Radiances%NumChannels
        j = j + 1
        nnInputData%Band_22_Radiances%values(channel, MIF) = ratio(j)
      enddo
    enddo
    print *, 'j: ', j
    print *, ratio(j)
  end subroutine StandardizeRadiances

  ! Dot product.
  real(r8) FUNCTION dot(X,Y)
    real(r8), dimension(:),intent(in):: x,y
    integer :: i,nn
    dot=0
    nn=size(x)
    IF (nn .NE. SIZE(y)) THEN 
      CALL MLSMessage(MLSMSG_error, ModuleName, &
           & 'The two vectors are not the same size!')
    ENDIF

    DO i=1,SIZE(x)
      dot = dot + x(i)*y(i)
    END DO
  END FUNCTION dot

  ! Sigmoid subroutine
  SUBROUTINE sigmoid(x,s) 
    real(r8),dimension(:),intent(in) :: x
    real(r8),dimension(:) :: s
    s = 1.0/(1.0+EXP(-x))
  END SUBROUTINE sigmoid

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

END MODULE NeuralNetUtils_M
! $Log$
! Revision 2.2  2021/02/19 00:29:46  pwagner
! repaired many bugs; still unsatisfactory imo
!
! Revision 2.1  2021/02/05 05:14:40  pwagner
! First commit
!
