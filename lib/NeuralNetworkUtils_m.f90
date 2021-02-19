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

MODULE NeuralNetUtils_M

  USE Dump_0, only: Dump
  USE HighOutput, only: OutputNamedValue
  USE MLSHDF5, only: SaveAsHDF5DS
  USE MLSKinds, only: R8
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error
  IMPLICIT NONE

  TYPE radiance_T
     character(len=64) :: signal
     !INTEGER :: NumMAFS ! <----- it's only ever going to be one MAF, right?
     INTEGER :: NumChannels
     INTEGER :: NumMIFs
     REAL(R8), dimension(:,:), allocatable :: values
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
     Integer,  dimension(:,:), allocatable :: Channels_In_Each_Band

     ! these are the dimensions for each part (assuming that we keep
     ! with Frank set up for Temperature retrieval)

     ! [42], number of pressure levels
     REAL(R8), dimension(:), allocatable :: Intercepts_Hidden_Labels_Layer
     REAL(R8), dimension(:), allocatable :: Output_Pressure_Levels
     Integer,  dimension(:), allocatable :: Output_Pressure_Levels_Indices

     ! [5078]. Number of Neurons in hidden layers
     REAL(R8), dimension(:), allocatable :: Intercepts_Hidden_Layer_1
     REAL(R8), dimension(:), allocatable :: Intercepts_Hidden_Layer_2

     ! [5078, 42] = [nNeurons, nLevels]
     REAL(R8), dimension(:,:), allocatable :: Weights_Hidden_Labels_Layer
     REAL(R8), dimension(:,:), allocatable :: Weights_Hidden_Layer_1
     REAL(R8), dimension(:,:), allocatable :: Weights_Hidden_Layer_2

     ! [42]. Number of levels
     REAL(R8), dimension(:), allocatable :: Normalization_Labels_Max
     REAL(R8), dimension(:), allocatable :: Normalization_Labels_Min
     
     ! [75]. Number of MIFs
     integer, dimension(:), allocatable :: MIFs

     ! [7575] number of 'variables', i.e. 2*25*75 + 51*75
     ! This is band 1[channel X MIF] + band 8 [channel X MIF ]  and 
     ! band 22 [channel X MIF]
     REAL(R8), dimension(:), allocatable :: Standardization_Brightness_Temperature_Mean
     REAL(r8), dimension(:), allocatable :: Standardization_Brightness_Temperature_Std

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
       & Radiance_T

  
  
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------


CONTAINS 
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

  FUNCTION NeuralNetFit(nnInputData, nnCoeffs, nHL, TYPE, &
    & FileID ) RESULT(prediction)

    ! Fortran version of the IDL code in `example_idl.pro' from Frank

    TYPE (NeuralNetInputData_T),intent(IN):: nnInputData
    TYPE (NeuralNetCoeffs_T), INTENT(IN)  :: nnCoeffs
    INTEGER,intent(in) :: nHL ! number of hidden layers. Will probably always be 2
    CHARACTER(len=*),INTENT(IN) :: TYPE ! type  may be either 'tanh' or 'sigmoid'
    REAL(r8), dimension(:), ALLOCATABLE :: prediction
    integer, optional, intent(in) :: FileID ! For optionally saving Datasets



    ! Local variables
    LOGICAL,SAVE :: first = .TRUE.
    INTEGER, SAVE:: nMIFs
    INTEGER, SAVE, dimension(2) :: nChans
    INTEGER, SAVE :: nSurfs, nVars, nNeurons

    ! Various working space arrays
    REAL(r8),dimension(:), allocatable :: working_space
    REAL(r8), dimension(:), ALLOCATABLE  :: NeuronValues, NeuronValues2

    real(r8) :: y_pred

    INTEGER :: c,n,m,jj,s,v! various counters
    INTEGER, dimension(2) :: dims2


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
    
    ! Optionally save this dataset for comparison with
    ! PyThon or iDl or whaTeVer you guys use these daYs
    if ( present(FileID) ) &
      & call SaveAsHDF5DS ( FileID, 'Radiances', working_space(1:jj-1) )

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


  ! Dot product.
  REAL(r8) FUNCTION dot(X,Y)
    REAL(r8), dimension(:),INTENT(IN):: x,y
    INTEGER :: i,nn
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
    REAL(r8),dimension(:),INTENT(in) :: x
    REAL(r8),dimension(:) :: s
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
! Revision 2.1  2021/02/05 05:14:40  pwagner
! First commit
!
