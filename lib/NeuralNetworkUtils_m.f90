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
  USE MLSKinds, only: R8
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error
  IMPLICIT NONE

  TYPE radiance_T
     character(len=64) :: signal
     !INTEGER :: NumMAFS ! <----- it's only ever going to be one MAF, right?
     INTEGER :: NumChannels
     INTEGER :: NumMIFs
     REAL(R8), DIMENSION(:,:),ALLOCATABLE :: values
  END TYPE radiance_T

  ! According to Frank, the 
  TYPE NeuralNetCoeffs_T
     ! Will be allocated in the caller. 
     ! [3]  Number of Bands
     character(len=32), dimension(:), allocatable :: Bands

     ! these are the dimensions for each part (assuming that we keep
     ! with Frank set up for Temperature retrieval)

     ! [42], number of pressure levels
     REAL(R8), DIMENSION(:),allocatable :: Intercepts_Hidden_Labels_Layer

     ! [5078]. Number of Neurons in hidden layers
     REAL(R8), DIMENSION(:),allocatable :: Intercepts_Hidden_Layer_1
     REAL(R8), DIMENSION(:),allocatable :: Intercepts_Hidden_Layer_2

     ! [5078, 41] = [nNeurons, nLevels]
     REAL(R8), DIMENSION(:,:),allocatable :: Weights_Hidden_Labels_Layer
     REAL(R8), DIMENSION(:,:),allocatable :: Weights_Hidden_Layer_1
     REAL(R8), DIMENSION(:,:),allocatable :: Weights_Hidden_Layer_2

     ! [41]. Number of levels
     REAL(R8), DIMENSION(:),allocatable :: Normalization_Labels_Max
     REAL(R8), DIMENSION(:),allocatable :: Normalization_Labels_Min
     
     ! [75]. Number of MIFs
     integer, dimension(:), allocatable :: MIFs

     ! [7575] number of 'variables', i.e. 2*25*75 + 51*75
     ! This is band 1[channel X MIF] + band 8 [channel X MIF ]  and 
     ! band 22 [channel X MIF]
     REAL(R8), DIMENSION(:),allocatable :: Standardization_Brightness_Temperature_Mean
     REAL(r8), DIMENSION(:),allocatable :: Standardization_Brightness_Temperature_Std

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
    if ( myDetails > 0 ) &
    call dump ( Coeffs%Bands, 'Bands' )

    call outputNamedValue ( 'Number of MIFs', size(Coeffs%MIFs) )
    if ( myDetails > 0 ) &
    call dump ( Coeffs%MIFs, 'MIFs' )
    
    call outputNamedValue ( 'Number of pressure levels', size(Coeffs%Intercepts_Hidden_Labels_Layer) )
    if ( myDetails > 0 ) then
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

  FUNCTION NeuralNetFit(nnInputData, nnCoeffs, nHL,TYPE) RESULT(prediction)

    ! Fortran version of the IDL code in `example_idl.pro' from Frank

    TYPE (NeuralNetInputData_T),intent(IN):: nnInputData
    TYPE (NeuralNetCoeffs_T), INTENT(IN)  :: nnCoeffs
    INTEGER,intent(in) :: nHL ! number of hidden layers. Will probably always be 2
    CHARACTER(len=*),INTENT(IN) :: TYPE ! type  may be either 'tanh' or 'sigmoid'
    REAL(r8), DIMENSION(:), ALLOCATABLE :: prediction



    ! Local variables
    LOGICAL,SAVE :: first = .TRUE.
    INTEGER, SAVE:: nMIFs
    INTEGER, SAVE, DIMENSION(2) :: nChans
    INTEGER, SAVE :: nSurfs, nVars, nNeurons

    ! Various working space arrays
    REAL(r8),DIMENSION(:),ALLOCATABLE :: working_space
    REAL(r8), DIMENSION(:), ALLOCATABLE  :: NeuronValues

    real(r8) :: y_pred

    INTEGER :: c,n,m,jj,s ! various counters
    INTEGER, dimension(2) :: dims2


    ! ----------- executable statements -----------------
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
    ALLOCATE( working_space(nVars )) ! nVars = 2*Chanels(bands 1&8) *  nMIFs + nChans(22)*nMIFs
    ALLOCATE(neuronValues(nNeurons))
    ALLOCATE(prediction(nSurfs))

    neuronValues = 0.0
    prediction = 0.0

    ! load the data into the working_space
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
        working_space( jj ) = nnInputData%Band_8_Radiances%Values(c,m)
        jj = jj + 1
      END DO
    END DO


    call OutputNamedValue ( 'nVars', nVars )
    call OutputNamedValue ( 'shape(WHL1)', shape(nnCoeffs%Weights_Hidden_Layer_1) )
      
    DO n=1,nNeurons
      neuronValues(n) = DOT( working_space, nnCoeffs%Weights_Hidden_Layer_1(n,:))
    END DO
    
    call OutputNamedValue ( 'num neurons', nNeurons )
    call OutputNamedValue ( 'shape(WHLL)', shape(nnCoeffs%Weights_Hidden_Labels_Layer) )
      
    DO s=1,nSurfs 

      ! neuron values 'activation' 
      IF (nHL .EQ. 1) THEN 

        IF (TRIM(TYPE) .EQ. 'tanh') THEN
          neuronValues = TANH(neuronValues + nnCoeffs%Intercepts_Hidden_Layer_1)
        ELSE IF(trim(TYPE) .EQ. 'sigmoid') THEN 
          call SIGMOID(neuronValues + nnCoeffs%Intercepts_Hidden_Layer_2,neuronValues)
        ENDIF

      ELSEIF (nHL .EQ. 2) THEN 

        IF (TRIM(TYPE) .EQ. 'tanh') THEN
          neuronValues = TANH(neuronValues + nnCoeffs%Intercepts_Hidden_Layer_2)
        ELSE IF(trim(TYPE) .EQ. 'sigmoid') THEN 
          call SIGMOID(neuronValues + nnCoeffs%Intercepts_Hidden_Layer_2,neuronValues)
        ENDIF

      ELSE 

        CALL MLSMessage(MLSMSG_error, ModuleName, &
             & "input var TYPE must equal either 'tanh' or 'sigmoid'")

      ENDIF

      y_pred = DOT( neuronValues,nnCoeffs%Weights_Hidden_Labels_Layer(:,s)) + &
           nnCoeffs%Intercepts_Hidden_Labels_Layer(s)
      

      ! Calculate the final product
       prediction(s) = y_pred * &
            & ( nnCoeffs%Normalization_Labels_Max(s) - &
            &   nnCoeffs%Normalization_Labels_Min(s)) + &
            &   nnCoeffs%Normalization_Labels_Min(s)

    ENDDO ! loop over pressure surfaces

  END function NeuralNetFit


  ! Dot product.
  REAL(r8) FUNCTION dot(X,Y)
    REAL(r8), DIMENSION(:),INTENT(IN):: x,y
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

  ! Sigmoid function 
  SUBROUTINE sigmoid(x,s) 
    REAL(r8),DIMENSION(:),INTENT(in) :: x
    REAL(r8),DIMENSION(:) :: s
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
