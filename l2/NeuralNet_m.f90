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
  use HGridsDatabase, only: L1BGeolocation
  use HighOutput, only: BeVerbose, Dump, LetsDebug, OutputNamedValue
  use Hunt_M, only: Hunt
  use MLSCommon, only: MLSChunk_T, MLSFile_T, UndefinedValue
  use MLSFiles, only: MLS_CloseFile, MLS_OpenFile, MLS_SFEnd, MLS_SFStart, &
    & AddInitializeMLSFile, Dump, HDFVersion_5
  use MLSFinds, only: FindFirst, FindLast
  use MLSKinds, only: R4, R8, Rv, Rk => R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
  use MLSSignals_M, only: GetSignalName
  use MLSStats1, only: MLSMax, MLSMin
  use MLSStrings, only: Capitalize, ReplaceNonAscii
  use MoreTree, only: Get_Field_Id
  use NeuralNetUtils_m, only: NeuralNetCoeffs_T, &
    & NeuralNetInputData_T, NeuralNetInputData_2_T, &
    & MatchedBinNum, MatchedMAF, MatchedStdRadiances, &
    & CheckMAFs, CheckTemperatures, Dump, NeuralNetFit, StandardizeRadiances
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
  ! ----------------------------------------------------------------------------
  ! This module performs the NeuralNet operation in the Level 2 software.
  ! This takes a measurement vector, 
  ! then returns a state vector with values calculated
  ! using a file of weight coefficients.
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
  !
  ! Another wrinkle is the infamous decision that MAFs start at 0 rather than 1
  ! ----------------------------------------------------------------------------

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

  ! The following must someday also be read from the file of weight coefficients
  !  integer, parameter              :: NumHiddenLayers  = 2
  ! The type of activation function is now read from the weights file
  ! (and anyway, it could be and in fact ultimately is 'relu'
  ! character(len=*), parameter     :: switchOver       = 'tanh' ! or 'sigmoid'
  ! Instead of plopping these constants into parameter declarations
  ! we'll either
  ! (1) Read them from the coefficients file; or
  ! (2) Make them fields of the l2cf's NeuralNet specification
  integer, parameter              :: NumChnnelsBand1  = 25
  integer, parameter              :: NumChnnelsBand8  = 25
  integer, parameter              :: NumChnnelsBand22 = 51
  integer, parameter              :: NumMIFs          = 75
  integer, parameter              :: NumLevels        = 42
  
  Logical, parameter              :: StandardizeRadiancesHere = .true.
  Logical, parameter              :: CheckTempsAndRads        = .false.
  Logical, parameter              :: UseMatchedRadiances      = .false.
  Logical, parameter              :: UseMatchedTemps          = .false.

contains ! =====     Public Procedures     =============================

  !---------------------------------------------------  NeuralNet  -----


  subroutine NeuralNet ( Key, VectorsDatabase, Chunk, FileDatabase )
    use Init_Tables_Module, only: L_H2O, L_HDF, L_Temperature, L_Radiance
    use Init_Tables_Module, only: F_State, F_Measurements, F_File, F_OutputSD
    integer, intent(in)                        :: Key
    type(vector_T), dimension(:), target       :: VectorsDatabase
    type (MLSChunk_T), intent(in)              :: Chunk
    type (MLSFile_T), dimension(:), pointer    :: FileDatabase
    integer                   :: Me = -1 ! String index for trace

    ! Local variables
    logical, dimension(23) :: NeededBands
    integer :: BinNum                 ! Loop counter
    logical :: DeeBug
    integer :: FIELD                  ! Entry in tree
    integer :: File
    integer :: FirstBin
    integer :: FirstProfile
    integer :: GSON                    ! Tree node
    ! These set how we match MAF numbers to the profile
    ! 3  methods are possible: Hunt, minloc, and FindFirst
    integer, parameter :: MMinLoc    = 1
    integer, parameter :: MHunt      = 2
    integer, parameter :: MFindFirst = 3
    integer, parameter :: howMatched = MMinLoc
    integer :: I                      ! Loop counter
    integer :: J                      ! Tree index
    integer :: LastBin
    integer :: LastProfile
    real(r4):: Lat
    integer :: MAF
    real(r4):: MaxLat
    real(r4):: MinLat
    integer :: MyFile                 ! FileDatabase index
    integer :: NumHiddenLayers        ! = 2
    integer :: NumMAFs                ! In this chunk
    real(r4):: Phi
    integer :: Profile
    integer :: RadfileID
    integer :: SON                    ! Tree node
    integer :: SIGNAL                 ! the signal we're looking for
    character(len=32) :: SignalStr    ! its string value
    integer :: Status
    real(Rv), allocatable, dimension(:) :: FranksValues
    real(Rv), allocatable, dimension(:) :: PrecisionValues
    real(Rv), allocatable, dimension(:) :: Values
    integer :: TempfileID
    integer, dimension(3) :: theseMAFs
    integer :: thisMAF
    integer :: thisProfile
    logical :: verbose

    type (Vector_T), pointer        :: Measurements ! The measurement vector
    type (Vector_T), pointer        :: State ! The state vector to fill
    type (Vector_T), pointer        :: OutputSD ! The Precision vector to fill

    type (MLSFile_T), pointer       :: CoeffsFile
    type (MLSFile_T), pointer       :: RadiancesFile
    type (VectorValue_T), pointer   :: Radiances ! The radiances for this band
    character(len=16)               :: RetrievedQty ! 'Temperature' or 'H2O'
    type (VectorValue_T), pointer   :: H2O ! The quantity to fill
    type (VectorValue_T), pointer   :: H2OPrecision ! Another quantity to fill
    type (VectorValue_T), pointer   :: Temperature ! The quantity to fill
    type (VectorValue_T), pointer   :: TemperaturePrecision ! Another quantity to fill
    type (VectorValue_T), pointer   :: Qty ! The quantity to fill
    type (VectorValue_T), pointer   :: Precision ! Another quantity to fill
    
    type(NeuralNetCoeffs_T)         :: Coeffs
    type(NeuralNetInputData_T)      :: NNMeasurements
    type(NeuralNetInputData_2_T)    :: NNMeasurements_2
    character (len=1024)            :: FileName
    real(rk), dimension(:,:), pointer :: GeodAngle
    real(rk), dimension(:,:), pointer :: GeodLat
    ! The next file holds the coefficients
    ! although we should try to supply the filename in the l2cf/PCF
    character (len=*), parameter    :: DefaultFileName = &
      & '/users/fwerner/Documents/database/trained_neural_nets/temperature/v1.14/1/' &
      & // &
      & 'Temperature_trained_neural_net.h5' ! For coefficients
    character (len=*), parameter    :: DefaultRadiances = &
      & 'Radiances_used_in_nn.h5' ! For radiances
    character (len=*), parameter    :: DefaultTemperatures = &
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
    DEEBUG = LetsDebug ( 'neu', 0 )
    verbose = BeVerbose ( 'neu', 0 )
    call trace_begin ( me, 'NeuralNet_m.NeuralNet', key, &
      & cond=toggle(gen) .and. levels(gen) > 1 )
    NeededBands = .false.
    nullify ( GeodAngle, GeodLat )

    do i = 2, nsons(key)
      son = subtree(i,key)
      gson = subtree(i,key) ! The argument
      if ( nsons(gson) > 1) gson = subtree(2,gson) ! Now value of said argument
      field = get_field_id(son)
      select case ( field )
      case ( f_measurements )
        measurements => VectorsDatabase(decoration(decoration(subtree(2,son))))
      case ( f_state )
        state => VectorsDatabase(decoration(decoration(subtree(2,son))))
      case ( f_outputSD )
        outputSD => VectorsDatabase(decoration(decoration(subtree(2,son))))
      case ( f_file )
        file = sub_rosa(gson)
        call outputNamedValue ( 'Processing file field', file )
      end select
    end do ! i = 2, nsons(key)

    if ( file > 0 ) then
      call get_string ( file, filename, strip=.true. )
      if ( verbose ) print *, 'filename read: ', trim(filename)
      myFile =  FileNameToID( trim(filename), FileDataBase ) 
    else
      myFile = 0
      if ( verbose ) print *, 'No filename read'
    end if
    if ( myFile < 1 ) then
      filename = defaultFileName
      if ( verbose ) print *, 'Must use default filename: ', trim(defaultFileName)
    endif
    if ( DeeBug ) then
      call outputNamedValue ( 'filename', trim(filename) )
    endif
    ! Loop over the quantities in the vectors

    if ( DeeBug ) then
      call outputNamedValue ( 'num quantities in state', state%template%noQuantities )
    endif
    do i = 1, state%template%noQuantities
      if ( state%quantities(i)%template%quantityType == l_temperature ) then
        TemperaturePrecision => OutputSD%quantities(i)
        Temperature => state%quantities(i)
        Temperature%values = 0.0_rv
        RetrievedQty = 'Temperature'
        NumHiddenLayers = 2
        if ( .not. associated(Temperature%BinNumber) ) &
          & allocate(Temperature%BinNumber(Temperature%template%NoInstances) )
        if ( .not. associated(Temperature%MAF) ) &
          & allocate(Temperature%MAF(Temperature%template%NoInstances) )
        Qty => Temperature
        Precision => TemperaturePrecision
      elseif ( state%quantities(i)%template%quantityType == l_h2o ) then
        H2OPrecision => OutputSD%quantities(i)
        H2O => state%quantities(i)
        H2O%values = 0.0_rv
        RetrievedQty = 'H2O'
        NumHiddenLayers = 4
        if ( .not. associated(H2O%BinNumber) ) &
          & allocate(H2O%BinNumber(H2O%template%NoInstances) )
        if ( .not. associated(H2O%MAF) ) &
          & allocate(H2O%MAF(H2O%template%NoInstances) )
        Qty => H2O
        Precision => H2OPrecision
      else
        if ( verbose ) &
          & print *, 'Skipping non-Temperature, non-H2O state quantity'
        cycle
      end if
    end do                          ! End loop over state quantities

      ! Now go through all the bands in the measurement vector that are in this radiometer
    do j = 1, measurements%template%noQuantities
      if ( measurements%quantities(j)%template%quantityType /= l_radiance ) then
        if ( verbose ) print *, 'measurement quantity type is not a radiance'
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
      NumMAFs = measurements%quantities(j)%template%NoInstances
      call GetSignalName ( signal, signalStr )
      signalStr = ReplaceNonAscii( signalStr, ' ' )
      select case (Capitalize(signalStr) )
      case ('R1A:118.B1F:PT.S0.FB25-1')
        if ( verbose ) print *, 'Got Band 1'
        NeededBands(1) = .true.
        ! Band 1 radiances are used in retrieving both T and H2O
        if ( RetrievedQty == 'Temperature' ) then
          NNMeasurements%Band_1_Radiances%NumChannels = NumChnnelsBand1
          NNMeasurements%Band_1_Radiances%NumMIFs     = NumMIFs
          allocate( &
            & NNMeasurements%Band_1_Radiances%&
            &   values(NumChnnelsBand1, NumMIFs) )
        else
          NNMeasurements_2%Band_1_Radiances%NumChannels = NumChnnelsBand1
          NNMeasurements_2%Band_1_Radiances%NumMIFs     = NumMIFs
          allocate( &
            & NNMeasurements_2%Band_1_Radiances%&
            &   values(NumChnnelsBand1, NumMIFs) )
        endif
      case ('R2:190.B2F:H2O.S0.FB25-1')
        if ( verbose ) print *, 'Got Band 2'
        NeededBands(2) = .true.
        NNMeasurements_2%Band_2_Radiances%NumChannels = NumChnnelsBand1
        NNMeasurements_2%Band_2_Radiances%NumMIFs     = NumMIFs
        allocate( &
          & NNMeasurements_2%Band_2_Radiances%&
          &   values(NumChnnelsBand1, NumMIFs) )
      case ('R2:190.B3F:N2O.S0.FB25-1')
        if ( verbose ) print *, 'Got Band 3'
        NeededBands(3) = .true.
        NNMeasurements_2%Band_3_Radiances%NumChannels = NumChnnelsBand1
        NNMeasurements_2%Band_3_Radiances%NumMIFs     = NumMIFs
        allocate( &
          & NNMeasurements_2%Band_3_Radiances%&
          &   values(NumChnnelsBand1, NumMIFs) )
      case ('R3:240.B8F:PT.S3.FB25-8')
        if ( verbose ) print *, 'Got Band 8'
        NeededBands(8) = .true.
        NNMeasurements%Band_8_Radiances%NumChannels = NumChnnelsBand8
        NNMeasurements%Band_8_Radiances%NumMIFs     = NumMIFs
        allocate( &
          & NNMeasurements%Band_8_Radiances%&
          &   values(NumChnnelsBand8, NumMIFs) )
      case ('R1A:118.B22D:PT.S0.DACS-4')
        if ( verbose ) print *, 'Got Band 22'
        NeededBands(22) = .true.
        NNMeasurements%Band_22_Radiances%NumChannels = NumChnnelsBand22
        NNMeasurements%Band_22_Radiances%NumMIFs     = NumMIFs
        allocate( &
          & NNMeasurements%Band_22_Radiances%&
          &   values(NumChnnelsBand22, NumMIFs) )
      case ('R2:190.B23D:H2O.S0.DACS-4')
        if ( verbose ) print *, 'Got Band 23'
        NeededBands(23) = .true.
        NNMeasurements_2%Band_23_Radiances%NumChannels = NumChnnelsBand22
        NNMeasurements_2%Band_23_Radiances%NumMIFs     = NumMIFs
        allocate( &
          & NNMeasurements_2%Band_23_Radiances%&
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
    if ( RetrievedQty == 'Temperature') then
      if ( .not. &
        & (NeededBands(1) .and. NeededBands(8) .and. NeededBands(22) ) &
        &  ) then
        call announce_error ( key, &
        & 'Must have all 3 needed bands to retrieve T; check your l2cf' )
      endif
    else ! H2O
      if ( .not. &
        & (NeededBands(1) .and. NeededBands(2) .and. &
        &  NeededBands(3) .and. NeededBands(23) ) &
        &  ) then
        call announce_error ( key, &
        & 'Must have all 4 needed bands to retrieve H2O; check your l2cf' )
      endif
    endif
    if ( myFile > 0 ) then
      CoeffsFile => FileDatabase(myFile)
    else
      CoeffsFile => AddInitializeMLSFile ( FileDatabase, name=fileName, &
        & shortName=RetrievedQty, &
        & type=l_hdf, content='NNCoeffs', HDFVersion=HDFVersion_5 )
    endif
    call MLS_OpenFile( CoeffsFile, Status )
    ! These are absolute profile numbers, i.e. they start at '1' only
    ! for the first chunk
    call Dump ( Chunk )
    call L1BGeoLocation ( filedatabase, 'GHz/GeodAngle    ', &
      & 'GHz', values2d=GeodAngle )
    call L1BGeoLocation ( filedatabase, 'GHz/GeodLat    ', &
      & 'GHz', values2d=GeodLat )
    ! We'll use the Geod angles directly
    call Hunt ( Qty%template%phi(1,:),  &
      & GeodAngle(36,Chunk%firstMAFIndex+1), firstProfile, &
      & allowTopValue=.true. )
    call Hunt ( Qty%template%phi(1,:),  &
      & GeodAngle(36,Chunk%lastMAFIndex+1), lastProfile, &
      & allowTopValue=.true. )
    ! Now we actually want the profile numbers in the entire
    ! precessing range, not just this chunk
    firstProfile = firstprofile + chunk%hGridOffsets(1)
    lastProfile  = lastprofile  + chunk%hGridOffsets(1)
    if ( verbose ) print *, 'chunkNumber  ', chunk%ChunkNumber
    if ( verbose ) print *, 'HGrid offset  ', chunk%hGridOffsets(1)
    if ( verbose ) print *, 'firstProfile ', firstProfile
    if ( verbose ) print *, 'lastProfile  ',  lastProfile
    allocate(Values(NumLevels))
    allocate(PrecisionValues(NumLevels))
    allocate(FranksValues(Qty%template%NoSurfs))
    
    ! Here's something new: must find first and last latitude bins
    minLat = mlsmin( Qty%template%geodLat(1,:) )
    maxLat = mlsmax( Qty%template%geodLat(1,:) )
    firstBin = FindLast ( BinArray - minLat < 0._r4 )
    lastBin  = FindFirst ( BinArray - maxLat > 0._r4 )
    ! if ( lastBin < firstBin ) lastBin = firstBin ! In case firstBin == NumBins
    ! *** We ought to be careful to handle all anomalous cases ***
    if ( firstBin < 1 ) then
      firstBin = NumBins
      lastBin = NumBins
    elseif ( lastBin < 1 ) then
      lastBin = NumBins
    elseif ( firstBin > lastBin ) then
      lastBin = firstBin
    endif
    if ( verbose ) print *, 'minLat ', minLat
    if ( verbose ) print *, 'maxLat ', maxLat
    if ( verbose ) print *, 'firstBin ', firstBin
    if ( verbose ) print *, 'lastBin  ',  lastBin
    do BinNum = firstBin, lastBin
      if ( verbose ) call outputNamedValue ( 'BinNum', BinNum )
      call ReadCoeffsFile ( CoeffsFile, BinNum, &
        & Coeffs, &
        & BinNum==firstBin ) ! MustAllocate Coeffs arrays on the 1st time through
      if ( verbose ) print *, 'Done reading Coeffs file'
      ! call OutputNamedValue ( 'Num profiles', Qty%template%NoInstances )
      do profile = 1, Qty%template%NoInstances ! firstProfile, lastProfile
        if ( verbose ) call outputNamedValue ( 'instance number in chunk ', profile )
        ! Does this profile fall within this bin num?
        if ( &
          & Qty%template%geodLat(1,profile) < BinArray(BinNum) &
          & .or. &
          & Qty%template%geodLat(1,profile) > BinArray(BinNum) + BinSz &
          & ) &
          & cycle
        if ( verbose ) print *, 'lat of profile ', profile, ' ', Qty%template%geodLat(1,profile)
        if ( verbose ) print *, 'bin range ', BinArray(BinNum), BinArray(BinNum)+BinSz
        thisProfile = profile + firstProfile - 1
        ! Find thisMAF matching profile
        ! Choose which method by which to match MAF to profile
        call Hunt ( GeodAngle(36,:), &
          & Qty%template%phi(1,profile), thisMaf, &
          & allowTopValue=.true. )
        if ( verbose ) call OutputNamedValue ( 'Hunt returned MAF', thisMAF )
        theseMAFs(2) = thisMAF
        ! Do we need to scale?
        theseMAFs(1:1) = minloc( &
          & abs(GeodAngle(36,:)-Qty%template%phi(1,profile)) &
          & + &
          & abs(GeodLat(36,:)-Qty%template%GeodLat(1,profile)) &
          & ) - 1
        thisMAF = theseMAFs(1) ! minloc insists on returning an array 
        if ( verbose ) call OutputNamedValue ( 'minloc returned MAF', thisMAF )
        thisMAF = FindFirst( &
          & (abs(GeodAngle(36,:)-Qty%template%phi(1,profile)) < 2.) &
          & .and. &
          & (abs(GeodLat(36,:)-Qty%template%GeodLat(1,profile)) < 1.) &
          & )
        if ( verbose ) call OutputNamedValue ( 'FindFirst returned MAF', thisMAF )
        theseMAFs(3) = thisMAF
        thisMAF = theseMAFs(howMatched)
        MAF = thisMAF - Chunk%firstMAFIndex ! + 1
        ! Must constrain MAF to be within range
        !   [0, NumMAFs-1]
        MAF = max( MAF, 0 )
        MAF = min( MAF, NumMAFs - 1 )
        thisMAF = MAF + Chunk%firstMAFIndex
        if ( verbose ) print *, 'thisprofile, thisMAF ', thisprofile, thisMAF
        if ( verbose ) print *, 'firstProfile, firstMAFIndex ', firstProfile, Chunk%firstMAFIndex
        if ( verbose ) print *, 'binNum, profile, MAF, chunk ', binNum, profile, MAF, Chunk%ChunkNumber
        if ( verbose ) print *, 'binNum, thisprofile, thisMAF, chunk ', binNum, thisprofile, thisMAF, Chunk%ChunkNumber
        if ( verbose ) print *, 'Num measurement quantities ', measurements%template%noQuantities
        if ( verbose ) print *, 'Qty radiometer ', Qty%template%radiometer
        Qty%BinNumber(profile) = binNum
        Qty%MAF(profile)       = thisMAF
        ! MAF and profile are indices inside the chunk, not absolute indices
        do j = 1, measurements%template%noQuantities
          ! call Dump ( measurements%quantities(j)%template )
          if ( measurements%quantities(j)%template%quantityType /= l_radiance ) cycle
          !
          ! Because there is no radiometer defined for Temperature,
          ! we can hardly compare it to whatever radiometer the measurement
          ! quantity uses
          ! if ( measurements%quantities(j)%template%radiometer /= Temperature%template%radiometer ) cycle
          !
          radiances => measurements%quantities(j)
          signal = measurements%quantities(j)%template%signal
          ! print *, 'Calling AssembleNNMeasurement for MAF ', thisMAF
          if ( RetrievedQty == 'Temperature' ) then
            call AssembleNNMeasurement ( NNMeasurements, &
              & radiances, signal, &
              & Coeffs%MIFs, Coeffs%Channels_In_Each_Band, MAF, DeeBug )
          else ! H2O
            call AssembleNNMeasurement_2 ( NNMeasurements_2, &
              & radiances, signal, &
              & Coeffs%MIFs, Coeffs%Channels_In_Each_Band, MAF, DeeBug )
          endif
        end do     
        ! End loop over bands
        if ( RetrievedQty == 'Temperature' ) then
          if ( all(NNMeasurements%Band_1_Radiances%values == 0._r8) ) &
            print *, 'All band 1 radiances vanish'
          if ( all(NNMeasurements%Band_8_Radiances%values == 0._r8) ) &
            print *, 'All band 8 radiances vanish'
          if ( all(NNMeasurements%Band_22_Radiances%values == 0._r8) ) &
            print *, 'All band 22 radiances vanish'
          if ( StandardizeRadiancesHere ) then
            call StandardizeRadiances ( NNMeasurements, &
              & Coeffs%Standardization_Brightness_Temperature_Mean, &
              & Coeffs%Standardization_Brightness_Temperature_Std )
          endif
          ! stop
          ! print *, 'Calling NeuralNetFit'
          if ( verbose ) print *, 'Using NeuralNetFit with our own raw Radiances'
          call NeuralNetFit ( NNMeasurements, &
            & Coeffs, NumHiddenLayers, Values, PrecisionValues, &
            & Debugging=.false. )
          do j=1, NumLevels
            Temperature%values(Coeffs%Output_Pressure_Levels_Indices(j), profile) = &
              & Values(j)
          enddo
          if ( any(PrecisionValues == UndefinedValue) ) &
            & TemperaturePrecision%values(:,profile) = UndefinedValue
          ! if ( StandardizeRadiancesHere ) &
          !  & call Dump ( FranksValues, 'Temps (Franks ANN)' )
        else ! H2O
          if ( all(NNMeasurements_2%Band_1_Radiances%values == 0._r8) ) &
            print *, 'All band 1 radiances vanish'
          if ( all(NNMeasurements_2%Band_2_Radiances%values == 0._r8) ) &
            print *, 'All band 2 radiances vanish'
          if ( all(NNMeasurements_2%Band_3_Radiances%values == 0._r8) ) &
            print *, 'All band 3 radiances vanish'
          if ( all(NNMeasurements_2%Band_23_Radiances%values == 0._r8) ) &
            print *, 'All band 23 radiances vanish'
          if ( StandardizeRadiancesHere ) then
            call StandardizeRadiances ( NNMeasurements_2, &
              & Coeffs%Standardization_Brightness_Temperature_Mean, &
              & Coeffs%Standardization_Brightness_Temperature_Std )
          endif
          ! stop
          ! print *, 'Calling NeuralNetFit'
          if ( verbose ) print *, 'Using NeuralNetFit with our own raw Radiances'
          call NeuralNetFit ( NNMeasurements_2, &
            & Coeffs, NumHiddenLayers, Values, PrecisionValues, &
            & Debugging=.false. )
          do j=1, NumLevels
            Qty%values(Coeffs%Output_Pressure_Levels_Indices(j), profile) = &
              & Values(j)
          enddo
          if ( any(PrecisionValues == UndefinedValue) ) &
            & Precision%values(:,profile) = UndefinedValue
          ! if ( StandardizeRadiancesHere ) &
          !  & call Dump ( FranksValues, 'Temps (Franks ANN)' )
        endif
      enddo ! Loop of profiles
    enddo ! Loop of BinNums
    call MLS_CloseFile( CoeffsFile, Status )


    ! stop
!debug
!call dump(Temperature%values, 'beforedivide')

    call Dump( Temperature%BinNumber, 'BinNumbers' )
    call Dump( Temperature%MAF, 'MAFs' )

    if ( BeVerbose( 'neu', 0 ) .or. DeeBug ) then
      call Dump( Temperature%values )
      call Dump( Coeffs )
      call output ( '*** Neural Net complete ***', advance='yes' )
    endif
    call trace_end ( 'NeuralNet_m.NeuralNet', &
      & cond=toggle(gen) .and. levels(gen) > 1 )
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
    ! Assemble masurements array to be used by RunNeuralNet
    subroutine AssembleNNMeasurement ( NNMeasurements, &
             & radiances, signal, MIFs, ChannelNums, MAF, DeeBug )
      type (VectorValue_T), pointer   :: Radiances ! The radiances for this band
      type(NeuralNetInputData_T)      :: NNMeasurements
      integer                         :: SIGNAL    ! the signal we're looking for
      integer, dimension(:)           :: MIFs      ! MIF numbers in radiances
      integer, dimension(:,:)         :: ChannelNums! channel numbers in radiances
      integer                         :: MAF ! Remember--1st MAF is 0
      logical, intent(in)             :: DeeBug
      ! Local variables
      integer                         :: channel
      integer                         :: chanNum
      integer                         :: MIF
      character(len=32)               :: SignalStr  ! its string value
      ! logical, parameter              :: DeeBug = .true.
      ! Executable
      call GetSignalName ( signal, signalStr )
      signalStr = ReplaceNonAscii( signalStr, ' ' )
      select case (Capitalize(signalStr) )
      case ('R1A:118.B1F:PT.S0.FB25-1')
        do MIF = 1, NumMIFs
          do channel = 1, NumChnnelsBand1
            chanNum = ChannelNums(channel, 1)
            NNMeasurements%Band_1_Radiances%values(channel,MIF) = &
              & radiances%value3(chanNum, MIFs(MIF), MAF+1)
          enddo
        enddo
        if ( DeeBug ) call OutputNamedValue ( 'Band 1 rad', &
          & NNMeasurements%Band_1_Radiances%values(1,1) )
      case ('R3:240.B8F:PT.S3.FB25-8')
        do MIF = 1, NumMIFs
          do channel = 1, NumChnnelsBand8
            chanNum = ChannelNums(channel, 2)
            NNMeasurements%Band_8_Radiances%values(channel,MIF) = &
              & radiances%value3(chanNum, MIFs(MIF), MAF+1)
          enddo
        enddo
        if ( DeeBug ) call OutputNamedValue ( 'Band 8 rad', &
          & NNMeasurements%Band_8_Radiances%values(1,1) )
      case ('R1A:118.B22D:PT.S0.DACS-4')
        do MIF = 1, NumMIFs
          do channel = 1, NumChnnelsBand22
            chanNum = ChannelNums(channel, 3)
            NNMeasurements%Band_22_Radiances%values(channel,MIF) = &
              & radiances%value3(chanNum, MIFs(MIF), MAF+1)
          enddo
        enddo
        if ( DeeBug ) call OutputNamedValue ( 'Band 22 rad', &
          & NNMeasurements%Band_22_Radiances%values(1,1) )
      case  default
        print *, 'Failed to recognize: ', Capitalize(signalStr)
      end select
    end subroutine AssembleNNMeasurement

    subroutine AssembleNNMeasurement_2 ( NNMeasurements, &
             & radiances, signal, MIFs, ChannelNums, MAF, DeeBug )
      ! real(r4), allocatable, dimension(:,:) :: NNMeasurements
      type (VectorValue_T), pointer   :: Radiances ! The radiances for this band
      type(NeuralNetInputData_2_T)    :: NNMeasurements
      integer                         :: SIGNAL    ! the signal we're looking for
      integer, dimension(:)           :: MIFs      ! MIF numbers in radiances
      integer, dimension(:,:)         :: ChannelNums! channel numbers in radiances
      integer                         :: MAF ! Remember--1st MAF is 0
      logical, intent(in)             :: DeeBug
      ! Local variables
      integer                         :: channel
      integer                         :: chanNum
      integer                         :: MIF
      character(len=32)               :: SignalStr  ! its string value
      ! logical, parameter              :: DeeBug = .true.
      ! Executable
      call GetSignalName ( signal, signalStr )
      signalStr = ReplaceNonAscii( signalStr, ' ' )
      select case (Capitalize(signalStr) )
      case ('R1A:118.B1F:PT.S0.FB25-1')
        do MIF = 1, NumMIFs
          do channel = 1, NumChnnelsBand1
            chanNum = ChannelNums(channel, 1)
            NNMeasurements%Band_1_Radiances%values(channel,MIF) = &
              & radiances%value3(chanNum, MIFs(MIF), MAF+1)
          enddo
        enddo
        if ( DeeBug ) call OutputNamedValue ( 'Band 1 rad', &
          & NNMeasurements%Band_1_Radiances%values(1,1) )
      case ('R2:190.B2F:H2O.S3.FB25-8')
        do MIF = 1, NumMIFs
          do channel = 1, NumChnnelsBand8
            chanNum = ChannelNums(channel, 2)
            NNMeasurements%Band_2_Radiances%values(channel,MIF) = &
              & radiances%value3(chanNum, MIFs(MIF), MAF+1)
          enddo
        enddo
        if ( DeeBug ) call OutputNamedValue ( 'Band 2 rad', &
          & NNMeasurements%Band_2_Radiances%values(1,1) )
      case ('R2:190.B3F:N2O.S3.FB25-8')
        do MIF = 1, NumMIFs
          do channel = 1, NumChnnelsBand8
            chanNum = ChannelNums(channel, 2)
            NNMeasurements%Band_3_Radiances%values(channel,MIF) = &
              & radiances%value3(chanNum, MIFs(MIF), MAF+1)
          enddo
        enddo
        if ( DeeBug ) call OutputNamedValue ( 'Band 3 rad', &
          & NNMeasurements%Band_3_Radiances%values(1,1) )
      case ('R1A:118.B23D:H2O.S0.DACS-4')
        do MIF = 1, NumMIFs
          do channel = 1, NumChnnelsBand22
            chanNum = ChannelNums(channel, 3)
            NNMeasurements%Band_23_Radiances%values(channel,MIF) = &
              & radiances%value3(chanNum, MIFs(MIF), MAF+1)
          enddo
        enddo
        if ( DeeBug ) call OutputNamedValue ( 'Band 23 rad', &
          & NNMeasurements%Band_23_Radiances%values(1,1) )
      case  default
        print *, 'Failed to recognize: ', Capitalize(signalStr)
      end select
    end subroutine AssembleNNMeasurement_2

  !------------------------------------------  InitializeNNMeasurements  -----
    subroutine InitializeNNMeasurements ! ( NNMeasurements )
      ! allocate and initialize masurements array to be used by RunNeuralNet
      ! real(r4), allocatable, dimension(:,:) :: NNMeasurements
    end subroutine InitializeNNMeasurements

  !------------------------------------------  ReadCoeffsFile  -----
    subroutine ReadCoeffsFile ( CoeffsFile, &
      & binNum, &
      & Coeffs, &
      & MustAllocate, &
      & Debugging )
      use MLSHDF5, only: IsHDF5DSPresent, LoadFromHDF5DS
      use MLSStrings, only: Asciify
      
      type (MLSFile_T)                        :: CoeffsFile
      integer, intent(in)                     :: BinNum
      ! real, dimension(:,:,:), allocatable     :: Coeffs
      type(NeuralNetCoeffs_T)                 :: Coeffs
      logical, intent(in)                     :: MustAllocate
      logical, optional, intent(in)           :: Debugging
      ! Local variables
      integer, dimension(2)                   :: shp    ! So we can reshape
      integer, dimension(3)                   :: start
      integer, dimension(3)                   :: count
      integer, dimension(3)                   :: stride
      integer, dimension(3)                   :: block
      character(len=35), dimension(:), allocatable     :: bands
      character(len=35), dimension(:), allocatable     :: charvalues
      integer, dimension(:), allocatable      :: Channels
      logical                                 :: DeeBug
      integer, dimension(:), allocatable      :: Levels
      real(r4), dimension(:), allocatable     :: values1
      real(r4), dimension(:,:), allocatable   :: values2
      real(r4), dimension(:,:,:), allocatable :: values3
      ! Executable
      DeeBug = .false.
      if ( present(Debugging) ) DeeBug = Debugging
      if ( DeeBug ) print *, 'Now in ReadCoeffsFile with BinNum ', BinNum
      ! stop
      ! Allocate the max space we'll need
      allocate ( values1(5078) )
      allocate ( values2(1, 7575) )
      allocate ( values3(1, 5078, 7575) )
      if ( MustAllocate ) then
        ! Allocate the fields of our datatype to be read
        if ( DeeBug ) print *, 'Allocating ..'
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
        ! Just for debugging
        allocate ( Coeffs%Means(18,7575) )
        allocate ( Coeffs%Stddevs (18,7575) )
      endif
      if ( .not. allocated(Coeffs%Channels_In_Each_Band) ) &
        & call announce_error(0, &
        & 'Allocation failed somehow' )
      ! Now read each dataset, 
      ! using a rank 2 or rank 3 temporary as appropriate
      ! Remember: 
      ! hyperslabs in hdf5 treat the "start" array as if it means "offset"
      !
      
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
      if ( DeeBug ) print *, 'Channels in each band'

      !
!       start = (/ binNum-1, 0, 0 /)
!       stride = (/ 1, 1, 1 /)
!       count = (/ 1, 5078, 0 /)
!       block = (/ 1, 1, 0 /)
!       print *, start(1:1), count(1:1), stride(1:1), block(1:1)
!       call LoadFromHDF5DS ( CoeffsFile, &
!         & "Intercepts_Hidden_Layer_1", &
!         & values2(:,1:count(2)), &
!         & start(1:2), count(1:2) )
! !        & start(1:2), count(1:2), stride(1:2), block(1:2) )
!         Coeffs%Intercepts_Hidden_Layer_1 = &
!           values2(1,1:count(2))
      deallocate ( values2 )

      start = (/ binNum-1, 0, 0 /)
      stride = (/ 1, 1, 1 /)
      count = (/ 1, 42, 0 /)
      block = (/ 1, 1, 0 /)
      allocate ( values2(1, 42) )
      call LoadFromHDF5DS ( CoeffsFile, &
        & "Intercepts_Hidden_Labels_Layer", &
        & values2, &
        & start(1:2), count(1:2), stride(1:2), block(1:2) )
!       & values2 )
        Coeffs%Intercepts_Hidden_Labels_Layer = &
          values2(1,:)
      if ( DeeBug ) print *, 'Intercepts_Hidden_Labels_Layer'

      deallocate ( values1 )
      allocate ( values1(42) )
      call LoadFromHDF5DS ( CoeffsFile, &
        & "Output_Pressure_Levels", &
        & values1 )
        Coeffs%Output_Pressure_Levels = &
          values1
      if ( DeeBug ) print *, 'Intercepts_Hidden_Labels_Layer'

      deallocate ( values2 )

      start = (/ binNum-1, 0, 0 /)
      stride = (/ 1, 1, 1 /)
      count = (/ 1, 5078, 0 /)
      block = (/ 1, 1, 0 /)
      allocate ( values2(1, 5078) )
      call LoadFromHDF5DS ( CoeffsFile, &
        & "Intercepts_Hidden_Layer_1", &
        & values2, &
        & start(1:2), count(1:2), stride(1:2), block(1:2) )
        Coeffs%Intercepts_Hidden_Layer_1 = &
          values2(1,:)

      if ( DeeBug ) print *, 'Loaded Intercepts_Hidden_Layer_1'
!       call LoadFromHDF5DS ( CoeffsFile, &
!         & "Intercepts_Hidden_Layer_2", &
!         & values1, &
!         & start, count, stride, block )
!         Coeffs%Intercepts_Hidden_Layer_2 = &
!           (:)
! 
      call LoadFromHDF5DS ( CoeffsFile, &
        & "Intercepts_Hidden_Layer_2", &
        & values2, &
        & start(1:2), count(1:2), stride(1:2), block(1:2) )
        Coeffs%Intercepts_Hidden_Layer_2 = &
          values2(1,:)

      if ( DeeBug ) print *, 'Loaded Intercepts_Hidden_Layer_2'

      deallocate ( values2 )

      ! Just for debugging
      start = (/ 0, 0, 0 /)
      stride = (/ 1, 1, 1 /)
      count = (/ 18, 7575, 0 /)
      block = (/ 1, 1, 0 /)
      allocate ( values2(18, 7575) )

      call LoadFromHDF5DS ( CoeffsFile, &
        & "Standardization_Brightness_Temperatures_Mean", &
        & values2, &
        & start(1:2), count(1:2), stride(1:2), block(1:2) )
        Coeffs%Means = &
          values2

      call LoadFromHDF5DS ( CoeffsFile, &
        & "Standardization_Brightness_Temperatures_Std", &
        & values2, &
        & start(1:2), count(1:2), stride(1:2), block(1:2) )
        Coeffs%Stddevs = &
          values2

      if ( DeeBug ) print *, 'Loaded means, Stddevs for debugging'

      ! Back to stuff just for this bin
      deallocate ( values2 )
      start = (/ binNum-1, 0, 0 /)
      stride = (/ 1, 1, 1 /)
      count = (/ 1, 7575, 0 /)
      block = (/ 1, 1, 0 /)
      allocate ( values2(1, 7575) )

      call LoadFromHDF5DS ( CoeffsFile, &
        & "Standardization_Brightness_Temperatures_Mean", &
        & values2, &
        & start(1:2), count(1:2), stride(1:2), block(1:2) )
        Coeffs%Standardization_Brightness_Temperature_Mean = &
          values2(1,:)

      if ( DeeBug ) print *, 'Standardization_Brightness_Temperatures_Mean'

      call LoadFromHDF5DS ( CoeffsFile, &
        & "Standardization_Brightness_Temperatures_Std", &
        & values2, &
        & start(1:2), count(1:2), stride(1:2), block(1:2) )
        Coeffs%Standardization_Brightness_Temperature_Std = &
          values2(1,:)
      if ( DeeBug ) print *, 'Standardization_Brightness_Temperatures_Std'

      deallocate ( values3 )

      start = (/ binNum-1, 0, 0 /)
      stride = (/ 1, 1, 1 /)
      count = (/ 1, 42, 5078 /)
      block = (/ 1, 1, 1 /)
      allocate ( values3(1, 42, 5078) )
      shp = (/ 5078, 42 /)
      call LoadFromHDF5DS ( CoeffsFile, &
        & "Weights_Hidden_Labels_Layer", &
        & values3, start, count, stride, block )
        Coeffs%Weights_Hidden_Labels_Layer = &
          Reshape( values3(1,:,:), shp, order=(/2,1/) )
      if ( DeeBug ) then
        print *, 'Weights_Hidden_Labels_Layer'
        print *, shape(values3)
        print *, shape(Coeffs%Weights_Hidden_Labels_Layer)
      endif
      ! stop

      deallocate ( values3 )

      start = (/ binNum-1, 0, 0 /)
      stride = (/ 1, 1, 1 /)
      count = (/ 1, 5078, 7575 /)
      block = (/ 1, 1, 1 /)
      allocate ( values3(1, 5078, 7575) )
      shp = (/ 7575, 5078 /)

      call LoadFromHDF5DS ( CoeffsFile, &
        & "Weights_Hidden_Layer_1", &
        & values3, start, count, stride, block )
        Coeffs%Weights_Hidden_Layer_1 = &
          Reshape( values3(1,:,:), shp, order=(/2,1/) )
      if ( DeeBug ) print *, 'Weights_Hidden_Layer_1'

      deallocate ( values3 )

      start = (/ binNum-1, 0, 0 /)
      stride = (/ 1, 1, 1 /)
      count = (/ 1, 5078, 5078 /)
      block = (/ 1, 1, 1 /)
      allocate ( values3(1, 5078, 5078) )
      shp = (/ 5078, 5078 /)
      call LoadFromHDF5DS ( CoeffsFile, &
        & "Weights_Hidden_Layer_2", &
        & values3, start, count, stride, block )
        Coeffs%Weights_Hidden_Layer_2 = &
          Reshape( values3(1,:,:), shp, order=(/2,1/) )
      if ( DeeBug ) print *, 'Weights_Hidden_Layer_2'

      deallocate ( values2 )

      start = (/ binNum-1, 0, 0 /)
      stride = (/ 1, 1, 1 /)
      count = (/ 1, 42, 0 /)
      block = (/ 1, 1, 0 /)
      allocate ( values2(1, 42) )
      call LoadFromHDF5DS ( CoeffsFile, &
        & "Normalization_Labels_Max", &
        & values2, start(1:2), count(1:2), &
        & stride(1:2), block(1:2) )
        Coeffs%Normalization_Labels_Max = &
          values2(1, :)
      if ( DeeBug ) print *, 'Normalization_Labels_Max'

      call LoadFromHDF5DS ( CoeffsFile, &
        & "Normalization_Labels_Min", &
        & values2, start(1:2), count(1:2), &
        & stride(1:2), block(1:2) )
        Coeffs%Normalization_Labels_Min = &
          values2(1, :)
      if ( DeeBug ) print *, 'Normalization_Labels_Min'
      
      ! The newer weights files include the dataset Activation_Function
      ! Older ones did not
      if ( .not. &
        & IsHDF5DSPresent ( CoeffsFile%fileID%f_id, "Activation_Function" ) &
        & ) then
        Coeffs%Activation_Function = 'tanh'
        if ( DeeBug ) print *, 'Activation_Function (default): ', &
          & trim(Coeffs%Activation_Function)
        return
      endif

      start = (/ 0, 0, 0 /)
      stride = (/ 1, 1, 1 /)
      count = (/ 1, 42, 0 /)
      block = (/ 1, 1, 0 /)
      allocate ( charvalues(1) )
      call LoadFromHDF5DS ( CoeffsFile, &
        & "Activation_Function", &
        & charvalues, start(1:1), count(1:1), &
        & stride(1:1), block(1:1) )
      if ( DeeBug ) then
        print *, 'Activation_Function (read): ', &
          & trim(charvalues(1))
          Coeffs%Activation_Function = &
            Asciify( charvalues(1), how='snip' )
        print *, 'Activation_Function (asciified): ', &
          & trim(Coeffs%Activation_Function)
      endif
    end subroutine ReadCoeffsFile

  !------------------------------------------  FileNameToID  -----
  ! Given a file name or a fragment of a file name found in the PCF
  ! return its index number into the file database.
  ! If the file is not found in the database, add it to the database.
  function FileNameToID ( fileName, DataBase )  result( ID )
    use Hdf, only: Dfacc_Rdonly
    use Intrinsic, only: L_HDF
    use MLSCommon, only: MLSFile_T
    use MLSFiles, only: HDFVersion_5, &
      & AddInitializeMLSFile, GetPCFromRef
    use MLSL2Options, only: Toolkit
    use MLSPCF2, only: MLSPCF_L2NeurNet_Start, MLSPCF_L2NeurNet_End

    ! Dummy arguments
    type (MLSFile_T), dimension(:), pointer :: DATABASE
    character(len=*), intent(in)            :: FileName ! full name or fragment
    integer                                 :: ID ! Index of file in database

    ! Local variables
    integer :: lun
    type (MLSFile_T), pointer   :: MLSFile
    integer :: mypcfEndCode
    integer :: PCBottom
    integer :: PCTop
    character(len=255) :: PCFFileName
    integer :: returnStatus
    integer :: Version
    ! Executable
    mypcfEndCode = 0
    lun = 0
    version = 1
    id = 0
    PCBottom = MLSPCF_L2NeurNet_Start
    PCTop    = MLSPCF_L2NeurNet_End
    print *, 'associated database: ', associated(dataBase)
    print *, 'len_trim(filename): ', len_trim(filename)
    ! if ( .not. associated(dataBase) .or. len_trim(filename) < 1 ) return
    if ( associated(database) ) then
      do id =1, size(database)
        if ( fileName == dataBase(id)%Name ) exit
        if ( fileName == dataBase(id)%ShortName ) exit
      enddo
      if ( id > size(database) ) id = 0
      if ( id > 0 ) return
    endif
    ! Must add file to database
    lun = GetPCFromRef( fileName, PCBottom, &
      & PCTop, &
      & TOOLKIT, returnStatus, Version, .false., &
      & exactName=PCFFileName )
    MLSFile => AddInitializeMLSFile( database, &
      & content='NeuralNetCoeffs', &
      & name=PCFFilename, shortName=fileName, &
      & type=l_hdf, access=dfacc_rdonly, HDFVersion=HDFVERSION_5 )
    call Dump ( MLSFile )
    ! Because we just now added a new item to the database, id is its new size
    id = size(database)
    print *, 'id: ', id
    print *, 'len_trim(exact filename): ', len_trim(PCFfilename)
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
! Revision 2.17  2022/04/28 20:46:32  pwagner
! First steps toward a 2nd n-n retrieval
!
! Revision 2.16  2022/04/13 21:33:03  pwagner
! Removed a lot of unneeded debugging stuff
!
! Revision 2.15  2021/10/14 22:25:25  pwagner
! Changes to acommodate setting precisions
!
! Revision 2.14  2021/07/28 23:43:12  pwagner
! Take care not to read Activation_Function unless it is present
!
! Revision 2.13  2021/07/22 23:18:26  pwagner
! Strip out extraordinary diagnostics relying on separate rad,temp files
!
! Revision 2.12  2021/07/08 23:33:03  pwagner
! Read Activation_Function from weights file; obey new api for NeuralNetFit; housekeeping
!
! Revision 2.11  2021/06/24 23:31:17  pwagner
! Coded 3 different approaches to matching profile to MAF
!
! Revision 2.10  2021/06/18 15:18:25  pwagner
! Distinguish StandardizeRadiancesHere (should always) from CheckTempsAndRads (optional)
!
! Revision 2.9  2021/06/10 23:49:11  pwagner
! Store BinNumber and MAF for retrieved qty
!
! Revision 2.8  2021/05/27 23:52:30  pwagner
! Now gets coeffs file from PCF; MAF index starts at 0; avoid use of L1MAF To and From functions in favor of L1BGeoLocation
!
! Revision 2.7  2021/05/18 15:53:56  pwagner
! Many bugs fixed; still in debug mode
!
! Revision 2.6  2021/04/01 23:52:57  pwagner
! Debugging til the cows come home (moo)
!
! Revision 2.5  2021/03/18 23:48:18  pwagner
! Fixed some more errors; added more debugging aids
!
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
