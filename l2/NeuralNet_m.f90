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

  use Global_Settings, only: L1MAFToL2Profile, L2ProfileToL1MAF
  use HighOutput, only: BeVerbose, Dump, LetsDebug, OutputNamedValue
  use MLSCommon, only: MLSChunk_T, MLSFile_T
  use MLSFiles, only: MLS_CloseFile, MLS_OpenFile
  use MLSFinds, only: FindFirst
  use MLSKinds, only: Rv
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
  use MLSStats1, only: MLSMax, MLSMin
  use MoreTree, only: Get_Field_Id
  ! use NeuralNetUtils_m, only: NeuralNet_T, RunNeuralNet
  use Output_M, only: Output
  use QuantityTemplates, only: Rt
  use String_Table, only: Display_String, Get_String
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
! NeuralNet          Given a measurement vector and file of coefficients
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

  logical, parameter :: countEmpty = .true. ! Except where overriden locally

contains ! =====     Public Procedures     =============================

  !---------------------------------------------------  NeuralNet  -----


  subroutine NeuralNet ( Key, VectorsDatabase, Chunk, FileDatabase )
    use Init_Tables_Module, only: L_Temperature, L_Radiance
    use Init_Tables_Module, only: F_State, F_Measurements, F_File
    integer, intent(in)                        :: Key
    type(vector_T), dimension(:), target       :: VectorsDatabase
    type (MLSChunk_T), intent(in)              :: Chunk
    type (MLSFile_T), dimension(:), pointer    :: FileDatabase
    integer                   :: Me = -1 ! String index for trace

    ! Local variables
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
    real(rt):: MaxLat
    real(rt):: MinLat
    integer :: MyFile                 ! FileDatabase index
    integer :: Profile
    integer :: SON                    ! Tree node
    integer :: SIGNAL                 ! the signal we're looking for
    integer :: Status
    integer :: thisMAF
    integer :: thisProfile

    type (Vector_T), pointer        :: Measurements ! The measurement vector
    type (Vector_T), pointer        :: State ! The state vector to fill

    type (MLSFile_T), pointer       :: CoeffsFile
    type (VectorValue_T), pointer   :: Radiances ! The radiances for this band
    type (VectorValue_T), pointer   :: Temperature ! The quantity to fill
    
    real, dimension(:,:,:), allocatable     :: Coeffs
    ! type (NeuralNet_T)              :: MyNeuralNet
    character (len=1024)            :: FileName
    character (len=*), parameter    :: DefaultFileName = &
      & '/users/fwerner/Documents/database/trained_neural_nets/temperature/v1.13/1/' &
      & // &
      & 'Temperature_trained_neural_net.h5'
    ! Beware if the following ever change
    ! For we will then need to read them from the file of coefficients
    ! For now we have 18 bins, each of size 10, spanning latitudes from the
    ! South Pole to the North
    integer, parameter              :: NumBins = 18
    real(rt), dimension(NumBins), parameter :: BinArray = &
      & (/ -90., -80., -70., -60., -50., -40., -30., -20., -10., &
      &      0.,  10.,  20.,  30.,  40.,  50.,  60.,  70.,  80. /)
    real(rt), parameter             :: BinSz = 10.
        
    ! Executable code
    DEEBUG = LetsDebug ( 'neu', 0 )
    call trace_begin ( me, 'NeuralNet_m.NeuralNet', key, &
      & cond=toggle(gen) .and. levels(gen) > 1 )

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
    else
      filename = defaultFileName
    end if
    myFile =  FileNameToID( trim(filename), FileDataBase ) 
    if ( DeeBug ) then
      call outputNamedValue ( 'filename', trim(filename) )
    endif
    ! Loop over the quantities in the vectors

    do i = 1, state%template%noQuantities
      if ( state%quantities(i)%template%quantityType /= l_temperature ) then
        call announce_error ( key, &
          & "state vector must contain only Temperature" )
      end if
      Temperature => state%quantities(i)
      Temperature%values = 0.0_rv
      ! Now go through all the bands in the measurement vector that are in this radiometer
    end do                          ! End loop over bands

    do j = 1, measurements%template%noQuantities
      if ( measurements%quantities(j)%template%quantityType /= l_radiance ) cycle
      if ( measurements%quantities(j)%template%radiometer /= Temperature%template%radiometer ) cycle
      radiances => measurements%quantities(j)
      signal = measurements%quantities(j)%template%signal
    end do                          ! End loop over bands

    CoeffsFile => FileDatabase(myFile)
    call MLS_OpenFile( CoeffsFile, Status )
    ! These are absolute profile numbers, i.e. they start at '1' only
    ! for the first chunk
    firstProfile = L1MAFToL2Profile ( Chunk%firstMAFIndex, FileDatabase )
    lastProfile  = L1MAFToL2Profile ( Chunk%lastMAFIndex , FileDatabase )
    
    ! Here's something new: must find first and last latitude bins
    minLat = mlsmin( temperature%template%geodLat(1,:) )
    maxLat = mlsmax( temperature%template%geodLat(1,:) )
    firstBin = FindFirst ( BinArray - minLat > 0._rt )
    lastBin  = FindFirst ( BinArray - maxLat > 0._rt )
    do BinNum = firstBin, lastBin
      call ReadCoeffsFile ( CoeffsFile, BinNum, &
        & Coeffs ) ! Should be & MyNeuralNet
      do profile = 1, temperature%template%NoInstances ! firstProfile, lastProfile
!         call InitializeNNMeasurements ( NNMeasurements )
        ! Does this profile fall within this bin num?
        if ( &
          & temperature%template%geodLat(1,profile) < BinArray(BinNum) &
          & .or. &
          & temperature%template%geodLat(1,profile) > BinArray(BinNum) + BinSz &
          & ) &
          & cycle
        thisProfile = profile + firstProfile - 1
        thisMAF = L2ProfileToL1MAF ( thisProfile, fileDatabase )
        MAF = thisMAF - Chunk%firstMAFIndex + 1
        ! MAF and profile are indices inside the chunk, not absolute indices
        do j = 1, measurements%template%noQuantities
          if ( measurements%quantities(j)%template%quantityType /= l_radiance ) cycle
          if ( measurements%quantities(j)%template%radiometer /= Temperature%template%radiometer ) cycle
          radiances => measurements%quantities(j)
          signal = measurements%quantities(j)%template%signal
!           call AssembleNNMeasurement ( NNMeasurements, &
!             & radiances%value3(:,:,MAF), &
!             & MyNeuralNet )
        end do                          ! End loop over bands
!         call RunNeuralNet ( NNMeasurements, Temperature%values, MyNeuralNet )
      enddo
    enddo
    call MLS_CloseFile( CoeffsFile, Status )
!debug
!call dump(Temperature%values, 'beforedivide')



    if ( BeVerbose( 'neu', 0 ) ) then
      call Dump( Temperature%values )
      call Dump( Coeffs )
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

      integer :: I

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
    subroutine AssembleNNMeasurements ! ( NNMeasurements, &
!             & radiances%value3(:,:,MAF), &
!             & MyNeuralNet )
      ! Assemble masurements array to be used by RunNeuralNet
      ! real(rt), allocatable, dimension(:,:) :: NNMeasurements
    end subroutine AssembleNNMeasurements

  !------------------------------------------  InitializeNNMeasurements  -----
    subroutine InitializeNNMeasurements ! ( NNMeasurements )
      ! allocate and initialize masurements array to be used by RunNeuralNet
      ! real(rt), allocatable, dimension(:,:) :: NNMeasurements
    end subroutine InitializeNNMeasurements

  !------------------------------------------  ReadCoeffsFile  -----
    subroutine ReadCoeffsFile ( CoeffsFile, &
      & binNum, &
      & Coeffs )
      use MLSHDF5, only: LoadFromHDF5DS
      
      type (MLSFile_T)                        :: CoeffsFile
      integer, intent(in)                     :: BinNum
      real, dimension(:,:,:), allocatable     :: Coeffs
      ! Local variables
      integer, dimension(3)                   :: start
      integer, dimension(3)                   :: count
      integer, dimension(3)                   :: stride
      integer, dimension(3)                   :: block
      real(rt), dimension(:,:), allocatable   :: values2
      real(rt), dimension(:,:,:), allocatable :: values3
      ! Executable
      ! Allocate the max space we'll need
      allocate ( values2(7575, 1) )
      allocate ( values3(7575, 5078, 1) )
      ! Now read each dataset, 
      ! using a rank 2 or rank 3 temporary as appropriate
      start = (/ 1, binNum, 0 /)
      stride = (/ 1, 1, 1 /)
      count = (/ 5078, 1, 0 /)
      block = (/ 1, 1, 0 /)
      call LoadFromHDF5DS ( CoeffsFile, &
        & "Standardization_Brightness_Temperatures_Mean", &
        & values2, &
        & start(1:2), count(1:2), stride(1:2), block(1:2) )
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
! Revision 2.1  2021/01/22 00:22:18  pwagner
! First commit
!
