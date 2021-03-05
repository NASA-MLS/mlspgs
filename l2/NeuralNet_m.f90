! Copyright 2021, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module NeuralNet_m              ! Use Neural Net Model to Retrieve State
!=============================================================================

  use Chunks_m, only: Dump
  use Global_Settings, only: L1MAFToL2Profile, L2ProfileToL1MAF
  use HDF, only: DFAcc_Create, DFAcc_RDOnly
  use HighOutput, only: BeVerbose, Dump, LetsDebug, OutputNamedValue
  use MLSCommon, only: MLSChunk_T, MLSFile_T
  use MLSFiles, only: MLS_CloseFile, MLS_OpenFile, MLS_SFEnd, MLS_SFStart, &
    & AddInitializeMLSFile, Dump, HDFVersion_5
  use MLSFinds, only: FindFirst, FindLast
  use MLSKinds, only: R4, R8, Rv
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
  use MLSSignals_M, only: DisplayRadiometer, GetSignalName
  use MLSStats1, only: MLSMax, MLSMin
  use MLSStrings, only: Capitalize, ReplaceNonAscii
  use MoreTree, only: Get_Field_Id
  use NeuralNetUtils_m, only: NeuralNetCoeffs_T, NeuralNetInputData_T, &
    & Dump, NeuralNetFit, StandardizeRadiances
  use Output_M, only: Output
  use QuantityTemplates, only: Dump
  use String_Table, only: Get_String
  use Toggles, only: Gen, Levels, Switches, Toggle
  use Trace_M, only: Trace_Begin, Trace_End
  use Tree, only: Decoration, Sub_Rosa, Subtree, Nsons, Subtree
  use VectorsModule, only: &
    & Dump, &
    & Vector_T, &
    & VectorValue_T
  ! This module performs the NeuralNet operation in the Level 2 software.
  ! This takes a measurement vector, 
  ! then returns a state vector with values calculated
  ! using a file of weight coefficients.

  implicit none
  private
  public :: NeuralNet

! === (start of toc) ===
! NeuralNet          Given a measurement vector and a file of coefficients
!                      calculate the resulting state vector
! === (end of toc) ===

! === (start of api) ===
! NeuralNet ( int key, 
!        *Vector_t VectorsDatabase(:),
!        *MLSChunk_T Chunk, 
!        *MLSFile_T FileDatabase(:) )
! === (end of api) ===
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

  ! The following must someday also be read from the file of n-n coefficients
  integer, parameter              :: NumHiddenLayers  = 2
  character(len=*), parameter     :: switchOver       = 'tanh' ! or 'sigmoid'
  integer, parameter              :: NumChnnelsBand1  = 25
  integer, parameter              :: NumChnnelsBand8  = 25
  integer, parameter              :: NumChnnelsBand22 = 51
  integer, parameter              :: NumMIFs          = 75
  integer, parameter              :: NumLevels        = 42
  
  Logical, parameter              :: StandardizeRadiancesHere = .true.

contains ! =====     Public Procedures     =============================

  !---------------------------------------------------  NeuralNet  -----


  subroutine NeuralNet ( Key, VectorsDatabase, Chunk, FileDatabase )
    use Init_Tables_Module, only: L_HDF, L_Temperature, L_Radiance
    use Init_Tables_Module, only: F_State, F_Measurements, F_File
    integer, intent(in)                        :: Key
    type(vector_T), dimension(:), target       :: VectorsDatabase
    type (MLSChunk_T), intent(in)              :: Chunk
    type (MLSFile_T), dimension(:), pointer    :: FileDatabase
    integer                   :: Me = -1 ! String index for trace

    ! Local variables
    logical, dimension(3) :: NeededBands
    integer :: BinNum                 ! Loop counter
    logical :: DeeBug
    integer :: FIELD                  ! Entry in tree
    integer :: File
    integer :: FirstBin
    integer :: FirstProfile
    integer :: I                      ! Loop counter
    integer :: J                      ! Tree index
    integer :: LastBin
    integer :: LastProfile
    integer :: MAF
    real(r4):: MaxLat
    real(r4):: MinLat
    integer :: MyFile                 ! FileDatabase index
    integer :: Profile
    integer :: RadfileID
    integer :: SON                    ! Tree node
    integer :: SIGNAL                 ! the signal we're looking for
    character(len=32) :: SignalStr    ! its string value
    integer :: Status
    real(Rv), allocatable, dimension(:) :: TemperatureValues
    integer :: TempfileID
    integer :: thisMAF
    integer :: thisProfile

    type (Vector_T), pointer        :: Measurements ! The measurement vector
    type (Vector_T), pointer        :: State ! The state vector to fill

    type (MLSFile_T), pointer       :: CoeffsFile
    type (MLSFile_T), pointer       :: RadiancesFile
    type (VectorValue_T), pointer   :: Radiances ! The radiances for this band
    type (VectorValue_T), pointer   :: Temperature ! The quantity to fill
    
    type(NeuralNetCoeffs_T)         :: Coeffs
    type(NeuralNetInputData_T)      :: NNMeasurements
    character (len=1024)            :: FileName
    character (len=*), parameter    :: DefaultFileName = &
      & '/users/fwerner/Documents/database/trained_neural_nets/temperature/v1.14/1/' &
      & // &
      & 'Temperature_trained_neural_net.h5' ! For coefficients
    character (len=*), parameter    :: DefaultRadiances = &
      & '/users/pwagner/' &
      & // &
      & 'Radiances_used_in_nn.h5' ! For radiances
    character (len=*), parameter    :: DefaultTemperatures = &
      & '/users/fwerner/Documents/database/trained_neural_nets/test_data/temperature/v1.13/' &
      & // &
      & 'Temperature_20190101.h5' ! For results to compare with
    ! Beware if the following ever change
    ! For we will then need to read them from the file of coefficients
    ! For now we have 18 bins, each of size 10, spanning latitudes from the
    ! South Pole to the North
    integer, parameter              :: NumBins = 18
    real(r4), dimension(NumBins), parameter :: BinArray = &
      & (/ -90., -80., -70., -60., -50., -40., -30., -20., -10., &
      &      0.,  10.,  20.,  30.,  40.,  50.,  60.,  70.,  80. /)
    real(r4), parameter             :: BinSz = 10.
        
    ! Executable code
    DEEBUG = LetsDebug ( 'neu', 0 ) .or. .true.
    call trace_begin ( me, 'NeuralNet_m.NeuralNet', key, &
      & cond=toggle(gen) .and. levels(gen) > 1 )
    NeededBands = .false.

    RadfileID = mls_sfstart ( trim(DefaultRadiances), DFAcc_Create, &
        &                                          hdfVersion=HDFVersion_5 )
    TempfileID = mls_sfstart ( trim(DefaultTemperatures), DFAcc_RDOnly, &
        &                                          hdfVersion=HDFVersion_5 )
    print *, 'Opened RadFile id ', RadfileID
    print *, 'Opened TempFile id ', TempfileID
    if ( RadfileID < 0 ) then
      call announce_error ( key, &
        & 'Failed to Open/Create Radiances File: ' // trim(DefaultRadiances) )
    endif
!       RadiancesFile => AddInitializeMLSFile ( FileDatabase, name=DefaultRadiances, &
!         & shortName='Radiances', access=Dfacc_Create,  &
!         & type=l_hdf, content='NNRadiances', HDFVersion=HDFVersion_5 )

    do i = 2, nsons(key)
      son = subtree(i,key)
      field = get_field_id(son)
      select case ( field )
      case ( f_measurements )
        measurements => VectorsDatabase(decoration(decoration(subtree(2,son))))
      case ( f_state )
        state => VectorsDatabase(decoration(decoration(subtree(2,son))))
      case ( f_file )
        file = sub_rosa(subtree(2,son))
      end select
    end do ! i = 2, nsons(key)

    if ( file > 0 ) then
      call get_string ( file, filename, strip=.true. )
      myFile =  FileNameToID( trim(filename), FileDataBase ) 
    else
      myFile = 0
    end if
    if ( myFile < 1 ) then
      filename = defaultFileName
    endif
    if ( DeeBug ) then
      call outputNamedValue ( 'filename', trim(filename) )
    endif
    ! Loop over the quantities in the vectors

    do i = 1, state%template%noQuantities
      if ( state%quantities(i)%template%quantityType /= l_temperature ) then
        ! call announce_error ( key, &
        !   & "state vector must contain only Temperature" )
        print *, 'Skipping non-Temperature state quantity'
        cycle
      end if
      Temperature => state%quantities(i)
      Temperature%values = 0.0_rv
      ! Now go through all the bands in the measurement vector that are in this radiometer
    end do                          ! End loop over bands

    do j = 1, measurements%template%noQuantities
      if ( measurements%quantities(j)%template%quantityType /= l_radiance ) then
        print *, 'measurement quantity type is not a radiance'
        cycle
      endif
!       We discover (to our spine-tingling horror)
!       that not all quantity templates have a radiometer; or bother to store
!       their radiometer index in the appropriate component.
!       if ( measurements%quantities(j)%template%radiometer /= Temperature%template%radiometer ) then
!         print *, 'measurement and Tmperature radiometers do not match '
!         call displayRadiometer ( measurements%quantities(j)%template%radiometer )
!         call displayRadiometer ( Temperature%template%radiometer )
!         cycle
!       endif
      radiances => measurements%quantities(j)
      signal = measurements%quantities(j)%template%signal
      call GetSignalName ( signal, signalStr )
      signalStr = ReplaceNonAscii( signalStr, ' ' )
      select case (Capitalize(signalStr) )
      case ('R1A:118.B1F:PT.S0.FB25-1')
        print *, 'Got Band 1'
        NeededBands(1) = .true.
        NNMeasurements%Band_1_Radiances%NumChannels = NumChnnelsBand1
        NNMeasurements%Band_1_Radiances%NumMIFs     = NumMIFs
        allocate( &
          & NNMeasurements%Band_1_Radiances%&
          &   values(NumChnnelsBand1, NumMIFs) )
      case ('R3:240.B8F:PT.S3.FB25-8')
        print *, 'Got Band 8'
        NeededBands(2) = .true.
        NNMeasurements%Band_8_Radiances%NumChannels = NumChnnelsBand8
        NNMeasurements%Band_8_Radiances%NumMIFs     = NumMIFs
        allocate( &
          & NNMeasurements%Band_8_Radiances%&
          &   values(NumChnnelsBand8, NumMIFs) )
      case ('R1A:118.B22D:PT.S0.DACS-4')
        print *, 'Got Band 22'
        NeededBands(3) = .true.
        NNMeasurements%Band_22_Radiances%NumChannels = NumChnnelsBand22
        NNMeasurements%Band_22_Radiances%NumMIFs     = NumMIFs
        allocate( &
          & NNMeasurements%Band_22_Radiances%&
          &   values(NumChnnelsBand22, NumMIFs) )
      case default
        call output ( 'Unrecognized or unwanted signal string:', advance='no' )
        call output ( Capitalize(trim(signalStr)), advance='yes' )
        call announce_error ( key, &
          & 'Signal error' )
      end select
    end do                          ! End loop over bands

    ! Sanity checks
    ! Must have all needed Bands
    if ( .not. any(NeededBands) ) then
      call announce_error ( key, &
        & 'Must have all 3 needed bands; check your l2cf' )
    endif
    if ( myFile > 0 ) then
      CoeffsFile => FileDatabase(myFile)
    else
      CoeffsFile => AddInitializeMLSFile ( FileDatabase, name=fileName, &
        & shortName='Temperature', &
        & type=l_hdf, content='NNCoeffs', HDFVersion=HDFVersion_5 )
    endif
    if ( DeeBug .and. .false. ) call Dump ( FileDatabase )
    call MLS_OpenFile( CoeffsFile, Status )
    ! call Dump( RadiancesFile )
!     call MLS_OpenFile( RadiancesFile, Status )
!     if ( Status /= 0 ) then
!       call announce_error ( key, &
!         & 'Failed to Open/Create Radiances File' )
!     endif
    ! These are absolute profile numbers, i.e. they start at '1' only
    ! for the first chunk
    call Dump ( Chunk )
    firstProfile = L1MAFToL2Profile ( Chunk%firstMAFIndex, FileDatabase )
    lastProfile  = L1MAFToL2Profile ( Chunk%lastMAFIndex , FileDatabase )
    print *, 'firstProfile ', firstProfile
    print *, 'lastProfile  ',  lastProfile
    allocate(TemperatureValues(NumLevels))
    
    ! Here's something new: must find first and last latitude bins
    minLat = mlsmin( temperature%template%geodLat(1,:) )
    maxLat = mlsmax( temperature%template%geodLat(1,:) )
    firstBin = FindLast ( BinArray - minLat < 0._r4 )
    lastBin  = FindFirst ( BinArray - maxLat > 0._r4 )
    print *, 'firstBin ', firstBin
    print *, 'lastBin  ',  lastBin
    do BinNum = firstBin, lastBin
      call ReadCoeffsFile ( CoeffsFile, BinNum, &
        & Coeffs, &
        & BinNum==firstBin ) ! MustAllocate Coeffs arrays on the 1st time through
      do profile = 1, temperature%template%NoInstances ! firstProfile, lastProfile
!         call InitializeNNMeasurements ( NNMeasurements )
        ! Does this profile fall within this bin num?
        print *, 'lat of profile ', profile, ' ', temperature%template%geodLat(1,profile)
        print *, 'bin range ', BinArray(BinNum), BinArray(BinNum)+BinSz
        if ( &
          & temperature%template%geodLat(1,profile) < BinArray(BinNum) &
          & .or. &
          & temperature%template%geodLat(1,profile) > BinArray(BinNum) + BinSz &
          & ) &
          & cycle
        thisProfile = profile + firstProfile - 1
        thisMAF = L2ProfileToL1MAF ( thisProfile, fileDatabase )
        MAF = thisMAF - Chunk%firstMAFIndex + 1
        print *, 'binNum, profile, MAF ', binNum, profile, MAF
        print *, 'thisprofile, thisMAF ', thisprofile, thisMAF
        print *, 'Num measurement quantities ', measurements%template%noQuantities
        print *, 'Temperature radiometer ', Temperature%template%radiometer
        ! MAF and profile are indices inside the chunk, not absolute indices
        do j = 1, measurements%template%noQuantities
          call Dump ( measurements%quantities(j)%template )
          if ( measurements%quantities(j)%template%quantityType /= l_radiance ) cycle
          !
          print *, 'Qty radiometer ', measurements%quantities(j)%template%radiometer
          ! Because there is no radiometer defined for Temperature,
          ! we can hardly compare it to whatever radiometer the measurement
          ! quantity uses
          ! if ( measurements%quantities(j)%template%radiometer /= Temperature%template%radiometer ) cycle
          !
          radiances => measurements%quantities(j)
          signal = measurements%quantities(j)%template%signal
          print *, 'Calling AssembleNNMeasurement'
          call AssembleNNMeasurement ( NNMeasurements, &
            & radiances, signal, &
            & Coeffs%MIFs, Coeffs%Channels_In_Each_Band, MAF )
        end do     
        ! End loop over bands
        if ( all(NNMeasurements%Band_1_Radiances%values == 0._r8) ) &
          print *, 'All band 1 radiances vanish'
        if ( all(NNMeasurements%Band_8_Radiances%values == 0._r8) ) &
          print *, 'All band 8 radiances vanish'
        if ( all(NNMeasurements%Band_22_Radiances%values == 0._r8) ) &
          print *, 'All band 22 radiances vanish'
        if ( StandardizeRadiancesHere ) then
          call StandardizeRadiances ( NNMeasurements, &
            & Coeffs%Standardization_Brightness_Temperature_Mean, &
            & Coeffs%Standardization_Brightness_Temperature_Std, &
            & TempFileID, thisMAF)
        endif
!         call RunNeuralNet ( NNMeasurements, c, MyNeuralNet )
        print *, 'Calling NeuralNetFit'
        print *, 'Rad FileID: ', RadfileID
        ! Temperature%values(8:49, profile) = NeuralNetFit ( NNMeasurements, &
        !   & Coeffs, NumHiddenLayers, switchOver )
        TemperatureValues = NeuralNetFit ( NNMeasurements, &
          & Coeffs, NumHiddenLayers, switchOver, &
          & RadfileID, TempFileID, thisMAF, thisProfile )
        do j=1, NumLevels
          Temperature%values(Coeffs%Output_Pressure_Levels_Indices(j), profile) = &
            & TemperatureValues(j)
        enddo
        call Dump ( Temperature%values(:, profile), 'Temps (nn)' )
      enddo
    enddo
    call MLS_CloseFile( CoeffsFile, Status )

    Status = MLS_SFEnd( RadfileID, hdfVersion=HDFVersion_5 )
    print *, 'Closed RadFile id ', RadfileID
    print *, Status
    if ( Status /= 0 ) then
      print *, 'Error in ending hdf access to file'
      stop
    endif

    Status = MLS_SFEnd( TempfileID, hdfVersion=HDFVersion_5 )
    print *, 'Closed TempFile id ', TempfileID
    print *, Status
    if ( Status /= 0 ) then
      print *, 'Error in ending hdf access to file'
      stop
    endif

!debug
!call dump(Temperature%values, 'beforedivide')



    if ( BeVerbose( 'neu', 0 ) .or. DeeBug ) then
      call Dump( Temperature%values )
      call Dump( Coeffs )
    endif
    call trace_end ( 'NeuralNet_m.NeuralNet', &
      & cond=toggle(gen) .and. levels(gen) > 1 )
    stop
  end subroutine NeuralNet

  ! ------------------ Private ---------------------------------
    ! ---------------------------------------------  ANNOUNCE_ERROR  -----
    subroutine announce_error ( wherewasit, &
      & extramessage, qty, extrainfo )

      use Moretree, only: StartErrorMessage

      integer, intent(in) :: WhereWasIt   ! Tree node WhereWasIt error was noticed
      character (len=*), intent(in), optional     :: EXTRAMESSAGE
      type (VectorValue_T), optional, intent(in)  :: QTY
      integer, intent(in), dimension(:), optional :: EXTRAINFO

      if ( present(extraMessage) ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & trim(extraMessage) )
      else
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Calling ANNOUNCE_ERROR' )
      end if
      call StartErrorMessage ( WhereWasIt )
      if ( present(ExtraMessage) )  call output(ExtraMessage, advance='yes')
      if ( present(extraInfo) ) &
        & call dump ( extraInfo, name='Extra info' )
      if ( present(qty) ) call dump( qty )
      call MLSMessage ( MLSMSG_Error, ModuleName, &
          & ' ' )
    end subroutine announce_error
  
  !------------------------------------------  AssembleNNMeasurements  -----
    subroutine AssembleNNMeasurement ( NNMeasurements, &
             & radiances, signal, MIFs, ChannelNums, MAF )
      ! Assemble masurements array to be used by RunNeuralNet
      ! real(r4), allocatable, dimension(:,:) :: NNMeasurements
      type (VectorValue_T), pointer   :: Radiances ! The radiances for this band
      type(NeuralNetInputData_T)      :: NNMeasurements
      integer                         :: SIGNAL    ! the signal we're looking for
      integer, dimension(:)           :: MIFs      ! MIF numbers in radiances
      integer, dimension(:,:)         :: ChannelNums! channel numbers in radiances
      integer                         :: MAF
      ! Local variables
      integer                         :: channel
      integer                         :: chanNum
      integer                         :: MIF
      character(len=32)               :: SignalStr  ! its string value
      ! Executable
      call GetSignalName ( signal, signalStr )
      signalStr = ReplaceNonAscii( signalStr, ' ' )
      select case (Capitalize(signalStr) )
      case ('R1A:118.B1F:PT.S0.FB25-1')
        do MIF = 1, NumMIFs
          do channel = 1, NumChnnelsBand1
            chanNum = ChannelNums(channel, 1)
            NNMeasurements%Band_1_Radiances%values(channel,MIF) = &
              & radiances%value3(chanNum, MIFs(MIF), MAF)
          enddo
        enddo
        
      case ('R3:240.B8F:PT.S3.FB25-8')
        do MIF = 1, NumMIFs
          do channel = 1, NumChnnelsBand8
            chanNum = ChannelNums(channel, 2)
            NNMeasurements%Band_8_Radiances%values(channel,MIF) = &
              & radiances%value3(chanNum, MIFs(MIF), MAF)
          enddo
        enddo
      case ('R1A:118.B22D:PT.S0.DACS-4')
        do MIF = 1, NumMIFs
          do channel = 1, NumChnnelsBand22
            chanNum = ChannelNums(channel, 3)
            NNMeasurements%Band_22_Radiances%values(channel,MIF) = &
              & radiances%value3(chanNum, MIFs(MIF), MAF)
          enddo
        enddo
      case  default
        print *, 'Failed to recognize: ', Capitalize(signalStr)
      end select
    end subroutine AssembleNNMeasurement

  !------------------------------------------  InitializeNNMeasurements  -----
    subroutine InitializeNNMeasurements ! ( NNMeasurements )
      ! allocate and initialize masurements array to be used by RunNeuralNet
      ! real(r4), allocatable, dimension(:,:) :: NNMeasurements
    end subroutine InitializeNNMeasurements

  !------------------------------------------  ReadCoeffsFile  -----
    subroutine ReadCoeffsFile ( CoeffsFile, &
      & binNum, &
      & Coeffs, &
      & MustAllocate )
      use MLSHDF5, only: LoadFromHDF5DS
      
      type (MLSFile_T)                        :: CoeffsFile
      integer, intent(in)                     :: BinNum
      ! real, dimension(:,:,:), allocatable     :: Coeffs
      type(NeuralNetCoeffs_T)                 :: Coeffs
      logical, intent(in)                     :: MustAllocate
      ! Local variables
      integer, dimension(2)                   :: shp    ! So we can reshape
      integer, dimension(3)                   :: start
      integer, dimension(3)                   :: count
      integer, dimension(3)                   :: stride
      integer, dimension(3)                   :: block
      character(len=35), dimension(:), allocatable     :: bands
      integer, dimension(:), allocatable      :: Channels
      integer, dimension(:), allocatable      :: Levels
      real(r4), dimension(:), allocatable     :: values1
      real(r4), dimension(:,:), allocatable   :: values2
      real(r4), dimension(:,:,:), allocatable :: values3
      ! Executable
      print *, 'Now in ReadCoeffsFile'
      ! stop
      ! Allocate the max space we'll need
      allocate ( values1(5078) )
      allocate ( values2(1, 7575) )
      allocate ( values3(1, 5078, 7575) )
      if ( MustAllocate ) then
        ! Allocate the fields of our datatype to be read
        print *, 'Allocating ..'
        allocate ( Coeffs%Channels_In_Each_Band(51,3) )
        allocate ( Coeffs%Intercepts_Hidden_Layer_1(5078) )
        allocate ( Coeffs%Intercepts_Hidden_Layer_2(5078) )
        allocate ( Coeffs%Weights_Hidden_Labels_Layer(5078, 42) )
        allocate ( Coeffs%Weights_Hidden_Layer_1     (7575, 5078) )
        allocate ( Coeffs%Weights_Hidden_Layer_2     (5078, 5078) )
        allocate ( Coeffs%Normalization_Labels_Max(42) )
        allocate ( Coeffs%Normalization_Labels_Min(42) )
        allocate ( Coeffs%Output_Pressure_Levels(42) )
        allocate ( Coeffs%Output_Pressure_Levels_Indices(42) )
        allocate ( Coeffs%Standardization_Brightness_Temperature_Mean(7575) )
        allocate ( Coeffs%Standardization_Brightness_Temperature_Std (7575) )
        allocate ( Coeffs%MIFs(75) )
        allocate ( Coeffs%Bands(3) )
      endif
      if ( .not. allocated(Coeffs%Channels_In_Each_Band) ) &
        & call announce_error(0, &
        & 'Allocation failed somehow' )
      ! Now read each dataset, 
      ! using a rank 2 or rank 3 temporary as appropriate
      ! Sadly, our implementation of hyperslabs seems to be buggy.
      !
      !  For an example of how to do it right, see
      ! /software/toolkit/ifc17/hdf5/fortran/examples//hyperslab.f90
      !
      ! Until we repair the hyperslab in LoadFromHDF5DS, we must resort to
      ! reading all 18 latitude bins every time (oof!)
      
      ! Bands
      allocate( Bands(3) )
      call LoadFromHDF5DS ( CoeffsFile, &
        & "Bands", &
        & Bands )
        Coeffs%Bands = &
          Bands

      ! MIFs
      allocate( Levels(75) )
      call LoadFromHDF5DS ( CoeffsFile, &
        & "MIFs", &
        & Levels )
        Coeffs%MIFs = &
          Levels

      ! Output pressure levels (indices)
      deallocate( Levels )
      allocate( Levels(42) )
      call LoadFromHDF5DS ( CoeffsFile, &
        & "Output_Pressure_Levels_Indices", &
        & Levels )
        Coeffs%Output_Pressure_Levels_Indices = &
          Levels

      ! Channels in each Band (indices)
      allocate( Channels(25) )
      call LoadFromHDF5DS ( CoeffsFile, &
        & "Channels_Band_#1", &
        & Channels )
        Coeffs%Channels_In_Each_Band(1:25,1) = &
          Channels(1:25)

      call LoadFromHDF5DS ( CoeffsFile, &
        & "Channels_Band_#2", &
        & Channels )
        Coeffs%Channels_In_Each_Band(1:25,2) = &
          Channels(1:25)

      deallocate( Channels )
      allocate( Channels(51) )
      call LoadFromHDF5DS ( CoeffsFile, &
        & "Channels_Band_#3", &
        & Channels )
        Coeffs%Channels_In_Each_Band(1:51,3) = &
          Channels(1:51)
      print *, 'Channels in each band'

      !
      start = (/ binNum, 1, 0 /)
      stride = (/ 1, 1, 1 /)
      count = (/ 1, 5078, 0 /)
      block = (/ 1, 1, 0 /)
!       print *, start(1:1), count(1:1), stride(1:1), block(1:1)
!       call LoadFromHDF5DS ( CoeffsFile, &
!         & "Intercepts_Hidden_Layer_1", &
!         & values2(:,1:count(2)), &
!         & start(1:2), count(1:2) )
! !        & start(1:2), count(1:2), stride(1:2), block(1:2) )
!         Coeffs%Intercepts_Hidden_Layer_1 = &
!           values2(1,1:count(2))
      deallocate ( values2 )

      allocate ( values2(18, 42) )
      call LoadFromHDF5DS ( CoeffsFile, &
        & "Intercepts_Hidden_Labels_Layer", &
        & values2 )
        Coeffs%Intercepts_Hidden_Labels_Layer = &
          values2(binNum,:)
      print *, 'Intercepts_Hidden_Labels_Layer'

      deallocate ( values1 )
      allocate ( values1(42) )
      call LoadFromHDF5DS ( CoeffsFile, &
        & "Output_Pressure_Levels", &
        & values1 )
        Coeffs%Output_Pressure_Levels = &
          values1
      print *, 'Intercepts_Hidden_Labels_Layer'

      deallocate ( values2 )

      allocate ( values2(18, 5078) )
      call LoadFromHDF5DS ( CoeffsFile, &
        & "Intercepts_Hidden_Layer_1", &
        & values2 )
        Coeffs%Intercepts_Hidden_Layer_1 = &
          values2(binNum,:)

      print *, 'Loaded Intercepts_Hidden_Layer_1'
!       call LoadFromHDF5DS ( CoeffsFile, &
!         & "Intercepts_Hidden_Layer_2", &
!         & values1, &
!         & start, count, stride, block )
!         Coeffs%Intercepts_Hidden_Layer_2 = &
!           (:)
! 
      call LoadFromHDF5DS ( CoeffsFile, &
        & "Intercepts_Hidden_Layer_2", &
        & values2 )
        Coeffs%Intercepts_Hidden_Layer_2 = &
          values2(binNum,:)

      print *, 'Loaded Intercepts_Hidden_Layer_2'

      deallocate ( values2 )

      allocate ( values2(18, 7575) )

      call LoadFromHDF5DS ( CoeffsFile, &
        & "Standardization_Brightness_Temperatures_Mean", &
        & values2 )
        Coeffs%Standardization_Brightness_Temperature_Mean = &
          values2(binNum,:)

      print *, 'Standardization_Brightness_Temperatures_Mean'

      call LoadFromHDF5DS ( CoeffsFile, &
        & "Standardization_Brightness_Temperatures_Std", &
        & values2 )
        Coeffs%Standardization_Brightness_Temperature_Std = &
          values2(binNum,:)
      print *, 'Standardization_Brightness_Temperatures_Std'

      deallocate ( values3 )

      allocate ( values3(18, 42,5078) )
      shp = (/ 5078, 42 /)
      call LoadFromHDF5DS ( CoeffsFile, &
        & "Weights_Hidden_Labels_Layer", &
        & values3 )
        Coeffs%Weights_Hidden_Labels_Layer = &
          Reshape( values3(binNum,:,:), shp, order=(/2,1/) )
      print *, 'Weights_Hidden_Labels_Layer'
      print *, shape(values3)
      print *, shape(Coeffs%Weights_Hidden_Labels_Layer)
      ! stop

      deallocate ( values3 )

      allocate ( values3(18, 5078, 7575) )
      shp = (/ 7575, 5078 /)

      call LoadFromHDF5DS ( CoeffsFile, &
        & "Weights_Hidden_Layer_1", &
        & values3 )
        Coeffs%Weights_Hidden_Layer_1 = &
          Reshape( values3(binNum,:,:), shp, order=(/2,1/) )
      print *, 'Weights_Hidden_Layer_1'

      deallocate ( values3 )

      allocate ( values3(18, 5078, 5078) )
      shp = (/ 5078, 5078 /)
      call LoadFromHDF5DS ( CoeffsFile, &
        & "Weights_Hidden_Layer_2", &
        & values3 )
        Coeffs%Weights_Hidden_Layer_2 = &
          Reshape( values3(binNum,:,:), shp, order=(/2,1/) )
      print *, 'Weights_Hidden_Layer_2'

      deallocate ( values2 )

      allocate ( values2(18, 42) )
      call LoadFromHDF5DS ( CoeffsFile, &
        & "Normalization_Labels_Max", &
        & values2 )
        Coeffs%Normalization_Labels_Max = &
          values2(binNum, :)
      print *, 'Normalization_Labels_Max'

      call LoadFromHDF5DS ( CoeffsFile, &
        & "Normalization_Labels_Min", &
        & values2 )
        Coeffs%Normalization_Labels_Min = &
          values2(binNum, :)
      print *, 'Normalization_Labels_Min'
    end subroutine ReadCoeffsFile

  !------------------------------------------  FileNameToID  -----
  function FileNameToID ( fileName, DataBase )  result(ID)

    ! Given filename, returns index; if name not found in db, returns 0

    ! Dummy arguments
    type (MLSFile_T), dimension(:), pointer :: DATABASE
    character(len=*), intent(in)            :: FileName
    integer                                 :: ID

    ! Local variables
    ! Executable
    id = 0
    if ( .not. associated(dataBase) .or. len_trim(filename) < 1 ) return
    do id =1, size(database)
      if ( fileName == dataBase(id)%Name ) exit
      if ( fileName == dataBase(id)%ShortName ) exit
    enddo
    if ( id > size(database) ) id = 0
  end function FileNameToID

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module NeuralNet_m
!=============================================================================

!
! $Log$
! Revision 2.4  2021/03/05 00:55:30  pwagner
! A tiny bit of progress
!
! Revision 2.3  2021/02/19 00:28:08  pwagner
! repaired many bugs; still unsatisfactory imo
!
! Revision 2.2  2021/02/05 05:16:37  pwagner
! Repaired many errors; others doubtless remain
!
! Revision 2.1  2021/01/22 00:22:18  pwagner
! First commit
!
