! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module QuantityTemplates         ! Quantities within vectors
!=============================================================================

  ! This module defines the `quantities' that make up vectors and their
  ! template information.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use HyperSlabs, only: Rerank
  use Dump_0, only: Dump
  use Expr_M, only: Expr_Check
  use HGridsDatabase, only: HGrid_T, Dump
  use HighOutput, only: OutputNamedValue
  use Intrinsic, only: L_Geodetic, L_GeodAltitude, L_None, L_Phitan, &
    & L_VMR, L_QTM, Lit_Indices, Phyq_Angle, Phyq_Dimensionless, &
    & Phyq_Frequency, Phyq_Indices, Phyq_Time, Phyq_Vmr
  use MLSKinds, only: Rt => R8 ! Rt Is "kind Of Real Components Of Template"
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
  use MLSFillValues, only: NaNFunction
  use MLSFinds, only: FindFirst
  use MLSStringLists, only: SwitchDetail
  use MLSStrings, only: Lowercase, WriteIntsToChars
  use Output_M, only: Output
  use String_Table, only: IsStringInTable
  use Toggles, only: Switches
  use Tree, only: Nsons, Subtree

  implicit NONE
  private

  public :: Epoch, QuantityTemplate_T, RT
  public :: AddQuantityTemplateToDatabase, InflateQuantityTemplateDatabase
  public :: CheckIntegrity, CopyQuantityTemplate, CreateGeolocationFields, &
    & CreateLatitudeFields, CreateLongitudeFields, &
    & DestroyGeolocationFields, DestroyQuantityTemplateContents, &
    & DestroyQuantityTemplateDatabase, Dump, GetHGridFromQuantity, &
    & ModifyQuantityTemplate, &
    & NullifyQuantityTemplate, PointQuantityToHGrid, QuantitiesAreCompatible, &
    & SetupNewQuantityTemplate, WriteAttributes

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

  ! Define some global parameters and data types.

  real(rt), parameter :: Epoch = 1993.0 ! Starting point for time references

  type QuantityTemplate_T

    ! Some administrative stuff

    integer :: Name = 0        ! Sub-rosa index of quantity name

    ! This integer is of an enumerated type describing what kind of
    ! quantity this is -- one of the l_lits of type t_quantityType
    ! in Init_Tables_Module, e.g. l_Temperature.

    integer :: QuantityType

    ! The dimensions of this quantity

    integer :: NoInstances     ! Number of along-track horizontal instances in
                               ! this quantity
    integer :: NoSurfs         ! Number of surfaces per instance
    integer :: NoChans = 1     ! Number of channels
    integer :: NoCrossTrack = 1! Number of cross-track horizontal instances in
                               ! this quantity (from xGrid)

    ! Flags describing the quantity

    logical :: Coherent        ! Do instances have same vertical coordinates?
    logical :: Stacked         ! Are instances true vertical profiles?
    logical :: Regular         ! Are all channels/heights represented

    ! This next one allows software using the vector quantities to be somewhat
    ! lazy and, for example, avoid interpolation.  Minor frame quantities are
    ! incoherent and unstacked, but may be regular or irregular.  However, not
    ! all incoherent unstacked quantities are minor frame quantities.

    logical :: MinorFrame      ! Is this a minor frame quantity.
    logical :: MajorFrame      ! Is this a major frame quantity.

    ! This one indicates whether log or linear interpolation should be used
    logical :: LogBasis                 ! If set use log
    real(rt) :: MinValue                ! Minimum value to consider if using log

    ! This information describes how much of the data is in the overlap
    ! regions if any.

    integer :: NoInstancesLowerOverlap
    integer :: NoInstancesUpperOverlap

    ! Misc. information
    real(rt) :: BadValue       ! Value used to flag bad/missing data
    integer :: Unit = PHYQ_VMR ! Unit quantity is in when scaled as below,
                               ! an l_lit of the type t_units.  Units are
                               ! defined in units.f90, but their names are
                               ! declared in intrinsic.f90, and their membership
                               ! in the type t_units is defined in init_tables_module.

    ! For regular quantities the number of elements of each instance
    ! is simply noSurfs*noChans.  For irregular ones it is less, but it is
    ! constant from instance to instance; this is that number.
    integer :: InstanceLen

    ! Vertical coordinate
    integer :: VerticalCoordinate=l_geodAltitude ! The vertical coordinate
                                  ! used.  These are l_lits of the type
                                  ! t_VGridCoord defined in Init_Tables_Module.
    logical :: SharedVGrid        ! Set if surfs is a pointer not a copy
    integer :: VGridIndex         ! Index of any vGrid used

    ! Surfs is dimensioned (noSurfs,noCrossTrack) for coherent quantities and
    ! (noSurfs, noInstances*noCrossTrack) for incoherent ones.  Pretending the
    ! values are dimensioned (noChans, noSurfs, noInstances), or
    ! (noChans, noSurfs, noInstances, noCrossTrack), the SURFS coordinate
    ! for the (:,i,j,:) values is surfs(i,1) for a coherent quantity or
    ! surfs(i,j) for an incoherent one.
    ! Surfs is allocatable instead of a pointer because the quantity template
    ! is copied into the vector value, not targeted therein.  If the VGrid
    ! is shared, deallocating a pointer would result in a dangling pointer.
    ! If noCrossTrack > 1, use the Surfs3 function to access it.
    real(rt), allocatable :: Surfs(:,:) ! zeta or meters, depending on verticalCoordinate

    ! Horizontal coordinates in the orbit plane.
    integer :: HorizontalCoordinate = l_phiTan ! The horizontal coordinate used.
                                        ! Either l_phiTan or l_time
    integer :: GrandTotalInstances      ! Total number of instances in destination output file
    ! for example MAF index, or profile index.
    integer :: hGridIndex               ! Index of any hGrid used
    type(hGrid_t), pointer :: The_HGrid => NULL()
    integer :: InstanceOffset           ! Ind of 1st non overlapped instance in output
    integer :: XGridIndex = 0           ! Index of any xGrid used

    ! First subscript values for GeoLocation component
    integer :: OrbitCoordinateIndex = 1 ! For spacecraft position
    integer :: LOSCoordinateIndex = 2   ! For line of sight
    real(rt), dimension(:,:,:), pointer :: Geolocation => NULL()

    ! Geolocation is dimensioned (*, 1, noInstances) for stacked quantities and
    ! (*, noSurfs, noInstances) for unstacked ones.  The Geolocation coordinate
    ! for the (*,i,j) value is geolocation(*,1,j) for a stacked quantity and
    ! geolocation(*,i,j) for an unstacked one.  The "*" is taken from either
    ! the OrbitCoordinateIndex or LOSCoordinateIndex component.

    real(rt), allocatable :: ECR(:,:,:) ! Meters

    ! ECR is dimensioned (3,1,noInstances*noCrossTrack) for stacked quantities
    ! and (3, noSurfs, noInstances*noCrossTrack) for unstacked ones.

    real(rt), allocatable :: Phi(:,:)   ! Degrees

    ! Phi is dimensioned (1, noInstances*noCrossTrack) for stacked quantities
    ! and (noSurfs, noInstances*noCrossTrack) for unstacked ones.  The PHI
    ! coordinate for the (i,j) value is phi(1,j) for a stacked quantity and
    ! phi(i,j) for an unstacked one.  Phi is either taken from or derived from
    ! Geolocation.  If noCrossTrack > 1, use the Phi3 function to access it.

    ! These other coordinates are dimensioned in the same manner as Phi.
    ! GeodLat and Lon are allocatable instead of pointers because the
    ! quantity template is copied into the vector value, not targeted
    ! therein.  If geolocations are computed, as for the magnetic field,
    ! if they were pointers the result would be dangling pointers.
    real(rt), allocatable :: GeodLat(:,:)            ! Degrees
    real(rt), allocatable :: Lon(:,:)                ! Degrees
    real(rt), pointer :: Time(:,:) => NULL()         ! Seconds since EPOCH
    real(rt), pointer :: SolarTime(:,:) => NULL()
    real(rt), pointer :: SolarZenith(:,:) => NULL()  ! Degrees
    real(rt), pointer :: LosAngle(:,:) => NULL()     ! Degrees

    ! Notwithstanding that the name of the latitude coordinate is GeodLat,
    ! it is sometimes geocentric latitude.  LatitudeCoordinate is either
    ! L_Geocentric or L_Geodetic.

    integer :: LatitudeCoordinate = L_Geodetic

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! For quantities having cross-track extent other than one, CrossAngles
    ! gives the angles from Phi, in the direction of the instrument.  Positive
    ! values are away from the instrument, negative values are toward the
    ! instrument.  The values come from the quantity's xGrid, which is
    ! reqired to be an explicit hGrid.
    real(rt), dimension(:), pointer :: CrossAngles => NULL()   ! Degrees

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! For quantities containing `channels' the following information may or
    ! may not be useful.

    ! Some quantities are on arbitrary freqency grids; these quantities refer
    ! to those.
    integer, pointer :: ChanInds(:) => NULL() ! Indices of values of Channels
                                        ! that are true
    logical, pointer :: Channels(:) => NULL() ! If /keepChannels is set
    integer :: fGridIndex               ! Index of any fGrid Index used
    real(rt), dimension(:), pointer :: frequencies => NULL() ! List of frequencies
                                        ! for Channels(ChanInds)
    integer :: FrequencyCoordinate      ! An enumerated type, e.g. FG_USBFreq
    real(rt) :: Lo                      ! Local oscillator frequency, MHz
    logical :: SharedFGrid              ! Set of frequencies are a pointer not a copy
    integer :: Sideband                 ! Associated sideband -1, 0, +1
    integer :: Signal                   ! Index into signals database

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Some families of quantities require special additional information.
    ! This is given here if needed.

    integer :: InstrumentModule ! Index in the Modules database in MLSSignals_m
    integer :: Radiometer       ! For ptan etc., index into radiometers database
    integer :: Reflector        ! For reflector efficiency etc. terms
    integer :: Molecule ! What molecule does this refer to? (One of the l_...
                        ! lits of type t_molecule in Molecules.)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! For irregular quantities, we have these arrays to
    ! help us navigate around the quantity.

    integer, dimension(:,:), pointer :: SurfIndex => NULL()
    integer, dimension(:,:), pointer :: ChanIndex => NULL()
    ! These are actually dimensioned (instanceLen, noInstances)
  contains
    procedure :: GeocLat => GetGeocLat
    procedure :: GeocLat3 => GetGeocLat3
    procedure :: GeodLat3 => GetLat3
    procedure :: IsQTM
    procedure :: Lon3 => GetLon3
    procedure :: Phi3 => GetPhi3
    procedure :: PutGeocLat
    procedure :: PutGeocLat3
    procedure :: PutLat
    procedure :: PutLat3
    procedure :: PutLon3
    procedure :: PutPhi3
    procedure :: PutSurfs3
    procedure :: Surfs3 => GetSurfs3
  end type QuantityTemplate_T

  interface Dump
    module procedure Dump_Quantity_Template, Dump_Quantity_Templates
  end interface

  interface ModifyQuantityTemplate
    module procedure ModifyQuantityTemplate_allocate
    module procedure ModifyQuantityTemplate_array, ModifyQuantityTemplate_sca
  end interface

  interface ReadAttributes
    module procedure ReadAttributes_QuantityTemplate
  end interface

  interface WriteAttributes
    module procedure WriteAttributes_QuantityTemplate
  end interface

  interface CHECKINTEGRITY
    module procedure CheckIntegrity_QuantityTemplate
  end interface

  ! Local procedures
  interface myValuesToField
    module procedure myValuesToField_2d_real
    module procedure myValuesToField_1d_dble
    module procedure myValuesToField_2d_dble
  end interface

  integer, parameter :: NUMMODS = 9
  character(*), dimension(NUMMODS), parameter :: MODIFIABLEFIELDS = (/&
    & 'surfs      ', &
    & 'phi        ', &
    & 'geodlat    ', &
    & 'lon        ', &
    & 'time       ', &
    & 'solartime  ', &
    & 'solarzenith', &
    & 'losangle   ', &
    & 'frequencies' &
    & /)

contains

 ! =====     Public Procedures     =============================

  ! Subroutines to deal with these quantitites

  ! ------------------------------  AddQuantityTemplateToDatabase  -----
  integer function AddQuantityTemplateToDatabase ( database, item )

  ! Add a quantity template to a database, or create the database if it
  ! doesn't yet exist

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

    ! Dummy arguments
    type (QuantityTemplate_T), dimension(:), pointer :: database
    type (QuantityTemplate_T), intent(in) :: item

    ! Local variables
    type (QuantityTemplate_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddQuantityTemplateToDatabase = newSize
  end function AddQuantityTemplateToDatabase

  ! ----------------------------  CopyQuantityTemplate  -----
  subroutine CopyQuantityTemplate ( Z, A, dontDestroy )
    ! This routine does a 'deep' copy of a quantity template.
    ! We don't need to do if often as typically only a shallow copy
    ! is required.  Note that this also follows any 'links' to h/v/xGrids
    ! and expands them too.
    use DeepCopy_m, only: DeepCopy
    type (QuantityTemplate_T), intent(inout) :: Z
    type (QuantityTemplate_T), intent(in)    :: A
    logical, optional, intent(in)            :: dontDestroy
    ! Internal variables
    logical                                  :: mustDestroy
    ! Executable code
    mustDestroy = .true.
    if ( present(dontDestroy) ) mustDestroy = .not. dontDestroy
    ! Destroy result
    if ( mustDestroy ) then
      call DestroyQuantityTemplateContents ( z )
    else
      nullify( z%the_HGrid )
    endif
    ! Setup result
    z%name = a%name
    call SetupNewQuantityTemplate ( z, noInstances=a%noInstances, &
      & noSurfs=a%noSurfs, noChans=a%noChans, noCrossTrack=a%noCrossTrack, &
      & coherent=a%coherent, stacked=a%stacked, regular=a%regular, &
      & instanceLen=a%instanceLen, minorFrame=a%minorFrame, &
      & majorFrame=a%majorFrame )
    ! Copy each other component -- tedious, but a shallow copy
    ! would lose newly allocated arrays
    ! 1st, scalars
    z%name                         = a%name
    z%quantityType                 = a%quantityType
    z%noInstances                  = a%noInstances            
    z%noSurfs                      = a%noSurfs                
    z%noChans                      = a%noChans                
    z%coherent                     = a%coherent               
    z%stacked                      = a%stacked                
    z%regular                      = a%regular                
    z%minorFrame                   = a%minorFrame             
    z%majorFrame                   = a%majorFrame             
    z%logBasis                     = a%logBasis               
    z%minValue                     = a%minValue               
    z%noInstancesLowerOverlap      = a%noInstancesLowerOverlap
    z%noInstancesUpperOverlap      = a%noInstancesUpperOverlap
    z%badValue                     = a%badValue               
    z%unit                         = a%unit                   
    z%instanceLen                  = a%instanceLen
    z%verticalCoordinate           = a%verticalCoordinate     
    z%sharedVGrid                  = a%sharedVGrid            
    z%vGridIndex                   = a%vGridIndex             
    z%xGridIndex                   = a%xGridIndex             
    z%horizontalCoordinate         = a%horizontalCoordinate             
    z%hGridIndex                   = a%hGridIndex              
    z%instanceOffset               = a%instanceOffset                    
    z%grandTotalInstances          = a%grandTotalInstances               
    z%fGridIndex                   = a%fGridIndex                
    z%frequencyCoordinate          = a%frequencyCoordinate             
    z%lo                           = a%lo                              
    z%sharedFGrid                  = a%sharedFGrid                     
    z%sideband                     = a%sideband                        
    z%signal                       = a%signal                          
    z%instrumentModule             = a%instrumentModule
    z%radiometer                   = a%radiometer      
    z%reflector                    = a%reflector       
    z%molecule                     = a%molecule        
    if ( associated(a%the_HGrid) ) &
  & z%the_HGrid                    => a%the_HGrid        

    ! Next, arrays
    if ( allocated(z%surfs)        .and. allocated(a%surfs) )        z%surfs =       a%surfs
    if ( allocated(z%phi)          .and. allocated(a%phi) )          z%phi =         a%phi
    if ( allocated(z%geodLat)      .and. allocated(a%geodLat) )      z%geodLat =     a%geodLat
    if ( allocated(z%lon)          .and. allocated(a%lon) )          z%lon =         a%lon
    if ( associated(z%time)        .and. associated(a%time) )        z%time =        a%time
    if ( associated(z%solarTime)   .and. associated(a%solarTime) )   z%solarTime =   a%solarTime
    if ( associated(z%solarZenith) .and. associated(a%solarZenith) ) z%solarZenith = a%solarZenith
    if ( associated(z%losAngle)    .and. associated(a%losAngle) )    z%losAngle =    a%losAngle
    if ( associated(z%chanInds)    .and. associated(a%chanInds) )    z%chanInds =    a%chanInds
    if ( associated(z%channels)    .and. associated(a%channels) )    z%channels =    a%channels
    call deepCopy ( z%frequencies, a%frequencies )
    if ( .not. z%regular ) then
      z%surfIndex = a%surfIndex
      z%chanIndex = a%chanIndex
    end if

  end subroutine CopyQuantityTemplate

  ! ------------------------------------  CreateGeolocationFields  -----
  subroutine CreateGeolocationFields ( Qty, NoSurfsToAllocate, What, Phi )
    ! Allocate Qty%GeodLat and Qty%Lon.

    type (QuantityTemplate_T), intent(inout) :: QTY
    integer, intent(in) :: NoSurfsToAllocate
    character(len=*), intent(in) :: What
    logical, intent(in), optional :: Phi

    call createLatitudeFields ( Qty, NoSurfsToAllocate, What )
    Qty%geodLat = NaNFunction(0._rt)
    call createLongitudeFields ( Qty, NoSurfsToAllocate, What )
    Qty%Lon = NaNFunction(0._rt)
    call createSurfsFields ( Qty, What )
    if ( present(phi) ) then
      if ( phi ) call createPhiFields ( Qty, What )
      Qty%phi = NaNFunction(0._rt)
    end if

  end subroutine CreateGeolocationFields

  ! ---------------------------------------  CreateLatitudeFields  -----
  subroutine CreateLatitudeFields ( Qty, NoSurfsToAllocate, What )
    ! Allocate Qty%GeodLat.

    type (QuantityTemplate_T), intent(inout) :: QTY
    integer, intent(in) :: NoSurfsToAllocate
    character(len=*), intent(in) :: What

    integer :: NumAlloc

    numAlloc = qty%noInstances * qty%noCrossTrack
    call allocate_test ( qty%geodLat, noSurfsToAllocate, numAlloc, trim(what) // "%geodLat", &
      & moduleName )

  end subroutine CreateLatitudeFields

  ! --------------------------------------  CreateLongitudeFields  -----
  subroutine CreateLongitudeFields ( Qty, NoSurfsToAllocate, What )
    ! Allocate Qty%Lon.

    type (QuantityTemplate_T), intent(inout) :: QTY
    integer, intent(in) :: NoSurfsToAllocate
    character(len=*), intent(in) :: What

    integer :: NumAlloc

    numAlloc = qty%noInstances * qty%noCrossTrack
    call allocate_test ( qty%lon, noSurfsToAllocate, numAlloc, trim(what) // "%lon", &
      & moduleName )

  end subroutine CreateLongitudeFields

  ! ------------------------------------------  CreatePhiFields  -----
  subroutine CreatePhiFields ( Qty, What )
    ! Allocate Qty%Phi

    type (QuantityTemplate_T), intent(inout) :: QTY
    character(len=*), intent(in) :: What

    call allocate_test ( qty%phi, merge(1,qty%noSurfs,qty%stacked), &
      & qty%noInstances*qty%noCrossTrack, moduleName, trim(what) // "%Phi" )

  end subroutine CreatePhiFields

  ! ------------------------------------------  CreateSurfsFields  -----
  subroutine CreateSurfsFields ( Qty, What )
    ! Allocate Qty%Surfs

    type (QuantityTemplate_T), intent(inout) :: QTY
    character(len=*), intent(in) :: What

    call allocate_test ( qty%surfs, qty%noSurfs, &
      & merge(1,qty%noInstances,qty%coherent), moduleName, &
      & trim(what) // "%surfs" )

  end subroutine CreateSurfsFields

  ! -----------------------------------  DestroyGeolocationFields  -----
  subroutine DestroyGeolocationFields ( Qty )
    ! Deallocate the latitude and longitude fields.  They are no longer
    ! shared with corresponding fields of HGrids, so we don't need to worry
    ! about shared HGrids here.
    use Allocate_Deallocate, only: Deallocate_Test
    ! Args
    type (QuantityTemplate_T), intent(inout) :: QTY
    ! Executable
    call deallocate_test ( qty%geodLat, 'Qty%GeodLat', moduleName )
    call deallocate_test ( qty%lon, 'Qty%Lon', moduleName )

  end subroutine DestroyGeolocationFields

  ! ----------------------------  DestroyQuantityTemplateContents  -----
  subroutine DestroyQuantityTemplateContents ( qty )
    use String_Table, only: Get_String
    ! Dummy argument
    type (QuantityTemplate_T), intent(inout) :: QTY
    ! Local variables
    character (len=80) :: typeStr
    logical :: Verbose
    character(63) :: What

    ! Executable code
    verbose = ( switchDetail(switches, 'qtmp' ) > -1 .or. switchDetail(switches, 'destroy' ) > -1 )
    ! May not destroy GeoLocations (until we discover why not)
    if ( qty%quantityType < 1 .or. qty%quantityType > size(lit_indices) ) return
    call get_string( lit_indices(qty%quantityType), typeStr )
    if ( lowercase( typeStr ) == 'geolocation' ) then
      call output( 'Unable to destroy this geoLocation quantity', advance='yes' )
      return
    elseif ( lowercase( typeStr ) == 'lostransfunc' ) then
      call output( 'Unable to destroy this LOSTransFunc quantity', advance='yes' )
      return
    endif
    if ( verbose ) then
      call dump( qty )
      call output( 'About to destroy this quantity', advance='yes' )
    end if
    if ( qty%name == 0 ) then
      what = "qty"
    else if ( verbose ) then
      call myGetString ( qty%name, what )
    else
      what = ''
    end if
    if ( index( lowercase( what ), 'phitan' ) > 0 .or. &
      & index( lowercase( what ), 'mif' ) > 0 ) then
      call output( 'Unable to destroy this particular quantity', advance='yes' )
      return
    endif

    if ( .not. qty%sharedVGrid ) then
      if ( verbose ) call output( 'About to deallocate surfs', advance='yes' )
      call deallocate_test ( qty%surfs, trim(what) // "%surfs", ModuleName )
    end if

    nullify ( qty%the_hGrid )

    ! If not stacked, there is no Phi.
    if ( qty%stacked ) then
      if ( verbose ) call output( 'About to deallocate stacked phi', advance='yes' )
      call deallocate_test ( qty%phi, trim(what) // '%phi', ModuleName )
    end if
    if ( verbose ) call output( 'About to deallocate geosdlat', advance='yes' )
    call deallocate_test ( qty%geodLat, trim(what) // '%geodLat', ModuleName )
    if ( verbose ) call output( 'About to deallocate lons', advance='yes' )
    call deallocate_test ( qty%lon, trim(what) // '%lon', ModuleName )
    if ( verbose ) call output( 'About to deallocate times', advance='yes' )
    call deallocate_test ( qty%time, trim(what) // '%time', ModuleName )
    if ( verbose ) call output( 'About to deallocate solartime', advance='yes' )
    call deallocate_test ( qty%solarTime, trim(what) // '%solarTime', ModuleName )
    if ( verbose ) call output( 'About to deallocate solarzenits', advance='yes' )
    call deallocate_test ( qty%solarZenith, trim(what) // '%solarZenith', ModuleName )
    if ( verbose ) call output( 'About to deallocate losangle', advance='yes' )
    call deallocate_test ( qty%losAngle, trim(what) // '%losAngle', ModuleName )
    if ( verbose ) call output( 'About to deallocate crossAngles', advance='yes' )
    call deallocate_test ( qty%crossAngles, trim(what) // '%crossAngles', ModuleName )

    if ( .not. qty%sharedFGrid ) then
      if ( verbose ) call output( 'About to deallocate chanInds', advance='yes' )
      call deallocate_test ( qty%chanInds, trim(what) // "%chanInds", ModuleName )
      if ( verbose ) call output( 'About to deallocate channels', advance='yes' )
      call deallocate_test ( qty%channels, trim(what) // "%channels", ModuleName )
      if ( verbose ) call output( 'About to deallocate freqs', advance='yes' )
      call deallocate_test ( qty%frequencies, trim(what) // "%frequencies", ModuleName )
    else
      nullify ( qty%chanInds, qty%channels, qty%frequencies )
    end if

    if ( .not. qty%regular ) then
      if ( verbose ) call output( 'About to deallocate surfindex', advance='yes' )
      call deallocate_test ( qty%surfIndex, trim(what) // "%surfIndex", ModuleName )
      if ( verbose ) call output( 'About to deallocate chanindex', advance='yes' )
      call deallocate_test ( qty%chanIndex, trim(what) // '%chanIndex', ModuleName )
    else
      nullify ( qty%surfIndex, qty%chanIndex )
    end if

  end subroutine DestroyQuantityTemplateContents

  ! ----------------------------  DestroyQuantityTemplateDatabase  -----
  subroutine DestroyQuantityTemplateDatabase ( database )
    use Output_m, only: Blanks
    ! Dummy argument

    use Allocate_Deallocate, only: Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc

    type (QuantityTemplate_T), dimension(:), pointer :: DATABASE
    ! Local variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: qtyIndex, s, status
    logical :: verbose

    ! Executable code
    verbose = ( switchDetail(switches, 'qtmp' ) > -1 .or. switchDetail(switches, 'destroy' ) > -1 )
    if ( associated(database) ) then
      if ( verbose ) call outputNamedValue( 'size(qty db)', size ( database ) )
      do qtyIndex = 1, size ( database )
        if ( verbose ) then
          call blanks ( 9, FillChar = '-' )
          call output ( ' Quantity index ' )
          call output ( qtyIndex )
          call blanks ( 1 )
          call blanks ( 9, FillChar = '-', advance='yes' )
        endif
        call DestroyQuantityTemplateContents ( database(qtyIndex) )
      end do
      s = size(database) * storage_size(database) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(database(1)), addr)
      deallocate ( database, stat=status )
      call test_deallocate ( status, ModuleName, "database", s, address=addr )
    end if
  end subroutine DestroyQuantityTemplateDatabase

  ! -------------------------------------  Dump_Quantity_Template  -----
  subroutine Dump_Quantity_Template ( Qty, Details, NoL2CF, What )

    use MLSSignals_M, only: Signals, Dump, GetRadiometerName, GetModuleName
    use Output_M, only: Blanks, NewLine
    use VGridsDatabase, only: Dump

    type(QuantityTemplate_T), intent(in) :: Qty
    integer, intent(in), optional :: Details ! <= 0 => Don't dump arrays
                                             ! >0   => Do signal, phi, surfs
                                             !         and frequency
                                             ! >1   => Dump all arrays
                                             ! Default 1
    logical, intent(in), optional :: NoL2CF  ! if TRUE => Don't dump L2-specific
                                             !  stuff
    character(*), intent(in), optional :: What ! In case you want to label it
    ! Local variables
    integer :: MyDetails
    logical :: myNoL2CF
    character (len=80) :: Str

    myDetails = 1
    if ( present(details) ) myDetails = details
    myNoL2CF = switchDetail(switches, 'nl2cf') > -1 ! .false.
    if ( present(NoL2CF) ) myNoL2CF = NoL2CF
    call output ( ' Name = ' )
    if ( isStringInTable( qty%name ) ) &
          & call myDisplayString ( qty%name )
    if ( .not. myNoL2CF .and. qty%quantityType > 0 ) then
      call output ( ' quantityType = ' )
      if ( isStringInTable(qty%quantityType, lit_indices) ) &
        &  call myDisplayString ( lit_indices(qty%quantityType) )
    else
      call output ( ' unknown quantityType' )
    end if
    if ( present(what) ) call output ( ' ' // trim(what) )
    call newline
    call Blanks( 5 )
    call output ( qty%noChans,              before=' NoChans = ' )
    call output ( qty%noSurfs,              before=' NoSurfs = ' )
    call output ( qty%noInstances,          before=' NoInstances = ' )
    call output ( qty%grandtotalinstances,  before=' grandtotalinstances = ' )
    call output ( qty%noCrossTrack,         before=' NoCrossTrack = ' )
    call newLine
    call blanks( 6 )
    call output ( trim(merge('  ','in',qty%coherent)) // 'coherent ' )
    call output ( trim(merge('   ','non',qty%stacked)) // 'stacked ' )
    call output ( trim(merge('  ','ir',qty%regular)) // 'regular ' )
    call output ( trim(merge('log-   ','linear-',qty%logBasis)) // 'basis ' )
    call output ( trim(merge('   ','non',qty%minorFrame)) // 'minorFrame ' )
    call output ( trim(merge('   ','non',qty%majorFrame)) // 'majorFrame', &
      & advance='yes' )
    call output ( '      NoInstancesLowerOverlap = ' )
    call output ( qty%noInstancesLowerOverlap )
    call output ( ' NoInstancesUpperOverlap = ' )
    call output ( qty%noInstancesUpperOverlap, advance='yes' )
    if ( .not. myNoL2CF .and. isStringInTable( qty%unit, phyq_indices) ) then
      call myDisplayString ( phyq_indices(qty%unit), before='      Unit = ' )
    else
      call blanks ( 5 )
    end if
    call output ( qty%badValue, before=' BadValue = ' )
    call output ( qty%InstanceLen, before=' InstanceLen = ', advance='yes' )
    if ( .not. myNoL2CF .and. &
      & isStringInTable( qty%horizontalCoordinate, lit_indices) ) &
      & call myDisplayString ( lit_indices(qty%horizontalCoordinate), &
      & before=   '      horizontal coordinate = ' )
    if ( isStringInTable( qty%horizontalCoordinate, lit_indices ) ) &
      & call myDisplayString ( lit_indices(qty%horizontalCoordinate), &
      & before=' latitude coordinate = ', advance='yes' )
    call output ( qty%hGridIndex, before='      hGridIndex = ' )
    call output ( qty%xGridIndex, before=' xGridIndex = ' )
    call newline
    call output ( qty%sharedVGrid, before='      sharedVGrid = ' )
    if ( qty%sharedVGrid ) then
      call output ( ' vGridIndex = ' )
      call output ( qty%vGridIndex )
    end if
    if ( .not. myNoL2CF .and. &
      & isStringInTable( qty%verticalCoordinate, lit_indices ) ) &
      & call myDisplayString ( lit_indices(qty%verticalCoordinate), &
      & before=' vertical coordinate = ' )
    call newLine
    call output ( qty%sharedFGrid, before='      sharedFGrid = ' )
    if ( qty%sharedFGrid ) then
      call output ( ' fGridIndex = ' )
      call output ( qty%fGridIndex )
    end if
    if ( isStringInTable( qty%frequencyCoordinate, lit_indices ) ) &
      & call myDisplayString ( lit_indices(qty%frequencyCoordinate), &
      & before=   ' frequency coordinate = ', advance='yes' )
    if ( qty%radiometer /= 0 .and. .not. myNoL2CF ) then
      call output ( '      Radiometer = ' )
      call GetRadiometerName ( qty%radiometer, str )
      call output ( trim(str), advance='yes' )
    end if
    if ( .not. myNoL2CF ) then
      if ( qty%quantityType == l_vmr .or. qty%instrumentModule /= 0 ) &
        & call output ( '     ' )
      if ( qty%quantityType == l_vmr ) then
        call output ( ' Molecule = ' )
        if ( isStringInTable( qty%molecule, lit_indices ) ) &
          & call myDisplayString ( lit_indices(qty%molecule) )
      end if
      if ( qty%instrumentModule /= 0 ) then
        call output ( ' Instrument Module = ' )
        call GetModuleName ( qty%instrumentModule, str )
        call output ( trim(str) )
      end if
      if ( qty%quantityType == l_vmr .or. qty%instrumentModule /= 0 ) &
        & call newLine
    end if

    if ( associated (qty%the_HGrid) ) then
      call Dump( qty%the_HGrid )
    else
      call output ( '      No The_HGrid', advance='yes' )
    end if
    if ( myDetails > 0 ) then
      if ( qty%signal /= 0 ) then
        call output ( '      Signal ' )
        call output ( qty%signal )
        call output ( ':', advance='yes' )
        if ( associated(qty%channels) ) then
          call dump ( signals(qty%signal), &
            & otherChannels=qty%channels )
        else
          call dump ( signals(qty%signal) )
        end if
      end if
      if ( allocated(qty%phi) ) then
        if ( qty%noCrossTrack == 1 ) then
          call dump ( reshape ( qty%phi, &
            & [ size(qty%geodLat,1),qty%noInstances] ), &
            & '      Phi' )
        else
          call dump ( reshape ( qty%phi, &
            & [ size(qty%geodLat,1),qty%noInstances,qty%noCrossTrack] ), &
            & '      Phi' )
        end if
      else
        call output ( '      No Phi' )
      end if

      if ( allocated(qty%surfs) ) then
        call dump ( qty%surfs, '      Surfs' )
      else
        call output ( '      No Surfs' )
      end if
      if ( myDetails > 1 ) then
        call maybe_dump_2_I ( qty%surfIndex, 'SurfIndex' )
        call maybe_dump_2_I ( qty%chanIndex, 'ChanIndex' )
        if ( allocated(qty%geodlat) ) then
          call dump ( reshape ( qty%geodlat, &
            & [ size(qty%geodLat,1),qty%noInstances,qty%noCrossTrack] ), &
            & '      GeodLat' )
        else
          call output ( '      No GeodLat' )
        end if
        if ( allocated(qty%lon) ) then
          call dump ( reshape ( qty%lon, &
            & [ size(qty%geodlat,1),qty%noInstances,qty%noCrossTrack] ), &
            & '      Lon' )
        else
          call output ( '      No Lon' )
        end if
        call maybe_dump_2_rt ( qty%time, 'Time' )
        call maybe_dump_2_rt ( qty%solarTime, 'SolarTime' )
        call maybe_dump_2_rt ( qty%solarZenith, 'SolarZenith' )
        call maybe_dump_2_rt ( qty%losAngle, 'LosAngle' )
        call maybe_dump_1_rt ( qty%crossAngles, 'CrossAngles' )
      end if
      if ( associated(qty%frequencies)  .and. .not. myNoL2CF ) then
        if ( isStringInTable( qty%frequencyCoordinate, lit_indices ) ) &
          & call myDisplayString ( lit_indices(qty%frequencyCoordinate), &
          & before='      FrequencyCoordinate = ', advance='yes' )
        call dump ( qty%frequencies, ' Frequencies ' )
      end if
    else
      if ( associated(qty%frequencies)  .and. .not. myNoL2CF .and. &
        & isStringInTable( qty%frequencyCoordinate, lit_indices ) ) &
        & call myDisplayString ( lit_indices(qty%frequencyCoordinate), &
        & before='      FrequencyCoordinate = ', advance='yes' )
    end if
  
  contains
  
    subroutine Maybe_Dump_2_I ( Array, Name )
      integer, intent(in), pointer :: Array(:,:)
      character(*), intent(in) :: Name
      if ( associated(array) ) then
        call dump ( array, '      ' // trim(name) // ' = ' )
      else
        call output ( '      No ' // trim(name), advance='yes' )
      end if
    end subroutine Maybe_Dump_2_I

    subroutine Maybe_Dump_1_RT ( Array, Name )
      real(rt), intent(in), pointer :: Array(:)
      character(*), intent(in) :: Name
      if ( associated(array) ) then
        call dump ( array, '      ' // trim(name) // ' = ' )
      else
        call output ( '      No ' // trim(name), advance='yes' )
      end if
    end subroutine Maybe_Dump_1_RT

    subroutine Maybe_Dump_2_RT ( Array, Name )
      real(rt), intent(in), pointer :: Array(:,:)
      character(*), intent(in) :: Name
      if ( associated(array) ) then
        call dump ( array, '      ' // trim(name) )
      else
        call output ( '      No ' // trim(name), advance='yes' )
      end if
    end subroutine Maybe_Dump_2_RT

  end subroutine Dump_Quantity_Template

  ! ------------------------------------  Dump_Quantity_Templates  -----
  subroutine Dump_Quantity_Templates ( Quantity_Templates, Details, NoL2CF, What )

    type(QuantityTemplate_T), intent(in) :: Quantity_Templates(:)
    integer, intent(in), optional :: Details ! <= 0 => Don't dump arrays
    !                                        ! >0   => Do dump arrays
    !                                        ! Default 1
    logical, intent(in), optional :: NoL2CF  ! if TRUE => Don't dump l2-specific
    character(*), intent(in), optional :: What ! In case you want to label it

    integer :: I

    call output ( size(quantity_templates), &
                & before='Quantity_Templates: SIZE = ', advance='yes' )
    do i = 1, size(quantity_templates)
      call output ( i, 4 )
      call output ( ':' )
      call dump_quantity_template ( quantity_templates(i), details, nol2cf, What )
    end do

  end subroutine Dump_Quantity_Templates

  ! ------------------------------ InflateQuantityTemplateDatabase -----
  integer function InflateQuantityTemplateDatabase ( database, extra )
    ! Make a quantity template database bigger by extra
    ! Return index of first new element

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

    ! Dummy arguments
    type (QuantityTemplate_T), dimension(:), pointer :: DATABASE
    integer, intent(in) :: EXTRA

    ! Local variables
    type (QuantityTemplate_T), dimension(:), pointer :: TEMPDATABASE

    include "inflateDatabase.f9h"
    InflateQuantityTemplateDatabase = firstNewItem
  end function InflateQuantityTemplateDatabase

! --------------------------------------------------------  IsQTM  -----
  pure logical function IsQTM ( Qty )
    class(QuantityTemplate_T), intent(in) :: Qty
    IsQTM = associated(qty%the_Hgrid)
    if ( IsQTM ) IsQTM = qty%the_Hgrid%type == L_QTM
    if ( isQTM ) IsQTM = allocated ( qty%the_Hgrid%QTM_Tree )
  end function IsQTM

  ! ------------------------------------  ModifyQuantityTemplate   -----
  ! This family modifies a quantity template's fields according to
  ! specified input
  ! This is something of a hack: normally we create a quantity template
  ! and its fields are filled according to its type and any vgrids,
  ! hgrids, and fgrids specified at that time
  ! However, in order to reverse-engineer an l2cf we may wish to
  ! go back later and override these fields by a Fill command

  subroutine ModifyQuantityTemplate_allocate  ( Z, FIELD, SHP, VALUESNODE, &
    & spread )
    ! This routine modifies any field whose name matches the field
    ! so that it takes the new values supplied by the source array
    
    ! Note that we assume the destination field is large enough
    ! to accomodate the source array
    
    ! How would you go about changing integer or l_ -valued fields?
    type (QuantityTemplate_T), intent(inout) :: Z
    character(len=*), intent(in)             :: field
    integer, dimension(:), intent(in)        :: SHP
    integer, intent(in)                      :: VALUESNODE   ! Tree node for values
    logical, intent(in)                      :: spread
    ! Executable
    if ( findFirst(MODIFIABLEFIELDS, lowercase(field)) < 1 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName // &
        & 'ModifyQuantityTemplate_allocate', &
        & trim(field) // " not a modifiable field" )
    end if
    select case(lowercase(field))
    case ( 'surfs' )
      if ( any(shape(z%surfs) /= shp) ) then
        call deallocate_test( z%surfs, 'template surfs', &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
        call allocate_test ( z%surfs, shp(1), shp(2), &
          & "template surfs", &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
      end if
      call myValuesToField( z%surfs, SHP, VALUESNODE, spread, PHYQ_Dimensionless )
    case ( 'phi' )
      if ( any(shape(z%phi) /= shp) ) then
        call deallocate_test( z%phi, 'template phi', &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
        call allocate_test ( z%phi, shp(1), shp(2), &
          & "template phi", &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
      end if
      call myValuesToField( z%phi, SHP, VALUESNODE, spread, phyq_angle )
    case ( 'geodlat' )
      if ( any(shape(z%geodlat) /= shp) ) then
        call deallocate_test( z%geodLat, 'template geodLat', &
          & ModuleName // '%ModifyQuantityTemplate_allocate' )
        call createLatitudeFields ( z, size(z%the_hGrid%geodLat,1), &
          & 'ModifyQuantityTemplate' )
      end if
      call myValuesToField( z%geodLat, SHP, VALUESNODE, spread, phyq_angle )
    case ( 'lon' )
      if ( any(shape(z%lon) /= shp) ) then
        call deallocate_test( z%lon, 'template lon', &
          & ModuleName // '%ModifyQuantityTemplate_allocate' )
        call createLongitudeFields ( z, size(z%the_hGrid%Lon,1), &
          & 'ModifyQuantityTemplate' )
      end if
      call myValuesToField( z%lon, SHP, VALUESNODE, spread, phyq_angle )
    case ( 'time' )
      if ( any(shape(z%time) /= shp) ) then
        call deallocate_test( z%time, 'template time', &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
        call allocate_test ( z%time, shp(1), shp(2), &
          & "template time", &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
      end if
      call myValuesToField( z%time, SHP, VALUESNODE, spread, phyq_time )
    case ( 'solartime' )
      if ( any(shape(z%solartime) /= shp) ) then
        call deallocate_test( z%solartime, 'template solartime', &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
        call allocate_test ( z%solartime, shp(1), shp(2), &
          & "template solartime", &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
      end if
      call myValuesToField( z%solartime, SHP, VALUESNODE, spread, phyq_time )
    case ( 'solarzenith' )
      if ( any(shape(z%solarzenith) /= shp) ) then
        call deallocate_test( z%solarzenith, 'template solarzenith', &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
        call allocate_test ( z%solarzenith, shp(1), shp(2), &
          & "template solarzenith", &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
      end if
      call myValuesToField( z%solarzenith, SHP, VALUESNODE, spread, phyq_angle )
    case ( 'losangle' )
      if ( any(shape(z%losangle) /= shp) ) then
        call deallocate_test( z%losangle, 'template losangle', &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
        call allocate_test ( z%losangle, shp(1), shp(2), &
          & "template losangle", &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
      end if
      call myValuesToField( z%losangle, SHP, VALUESNODE, spread, phyq_angle )
    case ( 'frequencies' )
      if ( .not. associated(z%frequencies) ) then
        call allocate_test ( z%frequencies, shp(1), &
          & "template frequencies", &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
      else if ( size(z%frequencies) /= shp(1) ) then
        call deallocate_test( z%frequencies, 'template frequencies', &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
        call allocate_test ( z%frequencies, shp(1), &
          & "template frequencies", &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
      end if
      call myValuesToField( z%frequencies, VALUESNODE, spread, phyq_frequency )
    case default
    end select
    if ( lowercase(field) == 'surfs' ) then
      z%NoSurfs = shp(1)
      if ( .not. z%coherent ) z%noInstances = shp(2)
    else if ( lowercase(field) == 'frequencies' ) then
      z%NoChans = shp(1)
    else
      z%noInstances = shp(2)
      if ( .not. z%stacked ) z%NoSurfs = shp(1)
    end if
  end subroutine ModifyQuantityTemplate_allocate

  subroutine ModifyQuantityTemplate_array  ( Z, FIELD, array, spread )
    ! This routine modifies any field whose name matches the field
    ! so that it takes the new values supplied by the source array
    
    ! Note that we assume the destination field is large enough
    ! to accomodate the source array
    
    ! How would you go about changing integer or l_ -valued fields?
    type (QuantityTemplate_T), intent(inout) :: Z
    character(len=*), intent(in)             :: field
    real(rt), dimension(:,:), intent(in)     :: array
    logical, intent(in)                      :: spread
    ! Local variables
    integer :: shp(2)
    ! Executable
    if ( findFirst(MODIFIABLEFIELDS, lowercase(field)) < 1 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName // &
        & 'ModifyQuantityTemplate_allocate', &
        & trim(field) // " not a modifiable field" )
    end if
    shp = shape(array)
    if ( spread .and. shp(1) == 1 ) then
      select case(lowercase(field))
      case ( 'surfs' )
        z%surfs(:,1:shp(2)) = array
      case ( 'phi' )
        z%phi(:,1:shp(2)) = array
      case ( 'geodlat' )
        z%geodLat(:,1:shp(2)) = array
      case ( 'lon' )
        z%lon(:,1:shp(2)) = array
      case ( 'time' )
        z%time(:,1:shp(2)) = array
      case ( 'solartime' )
        z%solartime(:,1:shp(2)) = array
      case ( 'solarzenith' )
        z%solarzenith(:,1:shp(2)) = array
      case ( 'losangle' )
        z%losangle(:,1:shp(2)) = array
      case ( 'frequencies' )
        z%frequencies(:) = array(1,1)
      case default
      end select
    else if ( spread .and. shp(2) == 1 ) then
      select case(lowercase(field))
      case ( 'surfs' )
        z%surfs(1:shp(1),:) = array
      case ( 'phi' )
        z%phi(1:shp(1),:) = array
      case ( 'geodlat' )
        z%geodLat(1:shp(1),:) = array
      case ( 'lon' )
        z%lon(1:shp(1),:) = array
      case ( 'time' )
        z%time(1:shp(1),:) = array
      case ( 'solartime' )
        z%solartime(1:shp(1),:) = array
      case ( 'solarzenith' )
        z%solarzenith(1:shp(1),:) = array
      case ( 'losangle' )
        z%losangle(1:shp(1),:) = array
      case ( 'frequencies' )
        z%frequencies(:) = array(1,1)
      case default
      end select
    else
      select case(lowercase(field))
      case ( 'surfs' )
        z%surfs(1:shp(1),1:shp(2)) = array
      case ( 'phi' )
        z%phi(1:shp(1),1:shp(2)) = array
      case ( 'geodlat' )
        z%geodLat(1:shp(1),1:shp(2)) = array
      case ( 'lon' )
        z%lon(1:shp(1),1:shp(2)) = array
      case ( 'time' )
        z%time(1:shp(1),1:shp(2)) = array
      case ( 'solartime' )
        z%solartime(1:shp(1),1:shp(2)) = array
      case ( 'solarzenith' )
        z%solarzenith(1:shp(1),1:shp(2)) = array
      case ( 'losangle' )
        z%losangle(1:shp(1),1:shp(2)) = array
      case ( 'frequencies' )
        z%frequencies(1:shp(1)) = array(:,1)
      case default
      end select
    end if
  end subroutine ModifyQuantityTemplate_array

  subroutine ModifyQuantityTemplate_sca  ( Z, FIELD, NEWVALUE )
    ! This routine modifies any field whose name matches the field
    ! so that it takes the new value
    
    ! How would you go about changing integer or l_ -valued fields?
    type (QuantityTemplate_T), intent(inout) :: Z
    character(len=*), intent(in)             :: field
    real(rt), intent(in)                     :: newvalue
    if ( findFirst(MODIFIABLEFIELDS, lowercase(field)) < 1 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName // &
        & 'ModifyQuantityTemplate_allocate', &
        & trim(field) // " not a modifiable field" )
    end if
    select case(lowercase(field))
    case ( 'surfs' )
      z%surfs = newvalue
    case ( 'phi' )
      z%phi = newvalue
    case ( 'geodlat' )
      z%geodLat = newvalue
    case ( 'lon' )
      z%lon = newvalue
    case ( 'time' )
      z%time = newvalue
    case ( 'solartime' )
      z%solartime = newvalue
    case ( 'solarzenith' )
      z%solarzenith = newvalue
    case ( 'losangle' )
      z%losangle = newvalue
    case ( 'frequencies' )
      z%frequencies = newvalue
    case default
    end select
  end subroutine ModifyQuantityTemplate_sca

  ! ------------------------------------  NullifyQuantityTemplate  -----
  subroutine NullifyQuantityTemplate ( IntentionallyNotUsed )
    ! Given a quantity template, nullify all the pointers within it
    type ( QuantityTemplate_T ), intent(out) :: IntentionallyNotUsed

    ! Executable code not needed because IntentionallyNotUsed is intent(out),
    ! and therefore undergoes default initialization as a consequence of
    ! argument association. Since all pointers within it have default
    ! initialization, they therefore become nullified.

    ! All non-pointer components of course become undefined and so must
    ! be explicitly defined after the call.

  end subroutine NullifyQuantityTemplate

  ! ---------------------------------------  GetHGridFromQuantity  -----
  subroutine GetHGridFromQuantity ( Qty )

  ! This routine associates HGrid information from an already defined quantity,
  ! except that if Qty%NoCrossTrack > 1, it creates a copy.

    use HGridsDatabase, only: HGrid_T, CreateEmptyHGrid
    use Toggles, only: Gen, Levels, Toggle
    use Trace_m, only: Trace_Begin, Trace_End

    ! Dummy argument
    type (QuantityTemplate_T), intent(inout) :: Qty

    ! Local variables
    integer :: Me = -1                  ! String index for trace
    type (hGrid_T), pointer :: HGrid

    ! Executable code

    call trace_begin ( me, "GetHGridFromQuantity", &
      & cond=toggle(gen) .and. levels(gen) > 2 )
    nullify ( HGrid )
    qty%the_hGrid             => HGrid              
    HGrid%noProfs             = qty%noInstances     
    HGrid%noProfsLowerOverlap = qty%noInstancesLowerOverlap
    HGrid%noProfsUpperOverlap = qty%noInstancesUpperOverlap
    call CreateEmptyHGrid ( hGrid )
    HGrid%time                => qty%time              
    HGrid%solarTime           => qty%solarTime         
    HGrid%solarZenith         => qty%solarZenith       
    HGrid%losAngle            => qty%losAngle          
    call trace_end ( "GetHGridFromQuantity", &
      & cond=toggle(gen) .and. levels(gen) > 2 )

  end subroutine GetHGridFromQuantity

  ! ---------------------------------------  PointQuantityToHGrid  -----
  subroutine PointQuantityToHGrid ( Qty )

  ! This routine associates HGrid information into an already defined quantity,
  ! except that if Qty%NoCrossTrack > 1, it creates a copy.

    use HGridsDatabase, only: HGrid_T
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
    use MoreMessage, only: MLSMessage
    use Toggles, only: Gen, Levels, Toggle
    use Trace_m, only: Trace_Begin, Trace_End

    ! Dummy argument
    type (QuantityTemplate_T), intent(inout) :: Qty

    ! Local variables
    integer :: I, J, K
    integer :: Me = -1                  ! String index for trace
    type (hGrid_T), pointer :: HGrid

    ! Executable code

    call trace_begin ( me, "PointQuantityToHGrid", &
      & cond=toggle(gen) .and. levels(gen) > 2 )

    hGrid => qty%the_hGrid

    if ( qty%noInstances/=hGrid%noProfs ) call MLSMessage ( MLSMSG_Error,&
      & ModuleName, "Size of HGrid not compatible with size of quantity" )

    call allocate_test ( qty%phi, merge(1,qty%noSurfs, qty%stacked), &
      & qty%noInstances*qty%noCrossTrack, 'qty%phi', ModuleName )
    if ( qty%stacked ) then
      if ( qty%noCrossTrack == 1 ) then
        qty%phi(:,:) = hGrid%phi
      else
        do k = 1, qty%noCrossTrack
          do j = 1, size(hGrid%phi,2)
            do i = 1, size(hGrid%phi,1)
              call qty%putPhi3 ( i, j, k, hGrid%phi(i,j) )
            end do
          end do
        end do
      end if
      call CreateGeolocationFields ( qty, size(hGrid%geodLat,1), 'Qty' )
      do k = 1, qty%noCrossTrack
        do j = 1, qty%noInstances ! size(qty%geodLat,2) would include cross-track
          do i = 1, size(qty%geodLat,1)
            call qty%putLat3 ( i, j, k, real(hGrid%geodLat(i,j),rt) )
            call qty%putLon3 ( i, j, k, real(hGrid%lon(i,j),rt) )
          end do
        end do
      end do
    else
      call CreateGeolocationFields ( qty, qty%noSurfs, 'Qty' )
      if ( qty%name /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & "Cannot copy hGrids into unstacked quantity %S" // &
        & "; assume Lat, Lon, Phi computed somehow.", qty%name)
      else
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & "Cannot copy hGrids into unstacked quantities " // &
          & "; assume Lat, Lon, Phi computed somehow.")
      end if
    end if
    qty%time                    => hGrid%time
    qty%solarTime               => hGrid%solarTime
    qty%solarZenith             => hGrid%solarZenith
    qty%losAngle                => hGrid%losAngle
    qty%noInstancesLowerOverlap = hGrid%noProfsLowerOverlap
    qty%noInstancesUpperOverlap = hGrid%noProfsUpperOverlap

    call trace_end ( "PointQuantityToHGrid", &
      & cond=toggle(gen) .and. levels(gen) > 2 )

  end subroutine PointQuantityToHGrid

  ! ------------------------------------  QuantitiesAreCompatible  -----
  logical function QuantitiesAreCompatible ( Qty_1, Qty_2, DifferentTypeOK, &
                                           & DifferentChansOK )
    use Intrinsic, only: L_QTM
    type(quantityTemplate_t), intent(in) :: Qty_1, Qty_2
    logical, intent(in), optional :: DifferentTypeOK ! Quantity, not hGrid, type
    logical, intent(in), optional :: DifferentChansOK

    QuantitiesAreCompatible = qty_1%quantityType == qty_2%quantityType
    if ( present(differentTypeOK) ) then
      if ( differentTypeOK ) QuantitiesAreCompatible = .true.
    end if
    if ( .not. QuantitiesAreCompatible ) return
    QuantitiesAreCompatible = associated(qty_1%the_hGrid) .eqv. &
                            & associated(qty_2%the_hGrid)
    if ( .not. QuantitiesAreCompatible ) return
    if ( associated(qty_1%the_hGrid) ) then
      QuantitiesAreCompatible = qty_1%the_hGrid%type == qty_2%the_hGrid%type
      if ( .not. QuantitiesAreCompatible ) return
      if ( qty_1%the_hGrid%type .eq. l_QTM ) then
        QuantitiesAreCompatible = allocated(qty_1%the_hGrid%QTM_tree%ZOT_in) .and. &
                                & allocated(qty_2%the_hGrid%QTM_tree%ZOT_in)
        if ( QuantitiesAreCompatible ) &
          & QuantitiesAreCompatible = size(qty_1%the_hGrid%QTM_tree%ZOT_in) == &
                                    & size(qty_2%the_hGrid%QTM_tree%ZOT_in)
      else
        QuantitiesAreCompatible = .not. allocated(qty_1%the_hGrid%QTM_tree%ZOT_in) .and. &
                                & .not. allocated(qty_2%the_hGrid%QTM_tree%ZOT_in)
      end if
    end if
    if ( .not. QuantitiesAreCompatible ) return
    QuantitiesAreCompatible = qty_1%noChans == qty_2%noChans
    if ( present(differentChansOK) ) then
      if ( differentChansOK ) QuantitiesAreCompatible = .true.
    end if
    if ( .not. QuantitiesAreCompatible ) return

    QuantitiesAreCompatible = &
      & qty_1%noInstances == qty_2%noInstances .and. &
      & qty_1%noSurfs == qty_2%noSurfs .and. &
      & qty_1%noCrossTrack == qty_2%noCrossTrack .and. &
      & qty_1%coherent .eqv. qty_2%coherent .and. &
      & qty_1%stacked .eqv. qty_2%stacked .and. &
      & qty_1%regular .eqv. qty_2%regular .and. &
      & qty_1%minorFrame .eqv. qty_2%minorFrame .and. &
      & qty_1%majorFrame .eqv. qty_2%majorFrame .and. &
      & qty_1%logBasis .eqv. qty_2%logBasis .and. &
      & qty_1%instanceLen == qty_2%instanceLen .and. &
      & qty_1%verticalCoordinate == qty_2%verticalCoordinate .and. &
      & qty_1%grandTotalInstances == qty_2%grandTotalInstances .and. &
      & qty_1%frequencyCoordinate == qty_2%frequencyCoordinate

  end function QuantitiesAreCompatible

  ! ----------------------------  ReadAttributes_QuantityTemplate  -----
  subroutine ReadAttributes_QuantityTemplate ( dsID, QT )
    ! Note:
    ! Most or all of the character-valued attributes are to be stored in the 
    ! quantity template as string table indexes or other indexes, e.g. signals
    ! Therefore we must do a bit of table lookups
    use Declaration_Table, only: Decls, Get_Decl, Phys_Unit_Name
    use Intrinsic, only: L_Dimensionless
    use MLSHDF5, only: GetHDF5Attribute
    use MLSSignals_m, only: GetRadiometerIndex, GetModuleIndex, &
      & GetSignalIndex

    ! Arguments
    integer, intent(in) :: dsID
    type(QuantityTemplate_T), intent(inout) :: qt

    type(decls) :: Decl
    character (len=80) :: Str
    ! Executable
    call GetHDF5AttrAsStrID ( dsID, 'TemplateName', qt%name )
    call GetHDF5AttrAsLitID ( dsID, 'tempQtyType', qt%quantityType )
    call GetHDF5Attribute ( dsID, 'noInstances', qt%noInstances )
    call GetHDF5Attribute ( dsID, 'noChans    ', qt%noChans     )
    call GetHDF5Attribute ( dsID, 'noSurfs    ', qt%noSurfs     )
    call GetHDF5Attribute ( dsID, 'coherent   ', qt%coherent     )
    call GetHDF5Attribute ( dsID, 'stacked    ', qt%stacked     )
    call GetHDF5Attribute ( dsID, 'regular    ', qt%regular     )
    call GetHDF5Attribute ( dsID, 'logBasis    ', qt%logBasis     )
    call GetHDF5Attribute ( dsID, 'minorFrame    ', qt%minorFrame     )
    call GetHDF5Attribute ( dsID, 'badValue    ', qt%badValue     )
    call GetHDF5Attribute ( dsID, 'tempQtyUnit', str )
    decl = get_decl ( str, phys_unit_name )
    if ( decl%type == phys_unit_name ) then
      qt%unit = decl%tree
    else ! ??? Should we emit an error message here ???
      qt%unit = l_dimensionless
    end if
    call GetHDF5Attribute ( dsID, 'instanceLen    ', qt%instanceLen     )
    call GetHDF5AttrAsLitID ( dsID, 'verticalCoordinate', qt%verticalCoordinate )
    call GetHDF5Attribute ( dsID, 'radiometer', str )
    call GetRadiometerIndex ( str, qt%radiometer )
    call GetHDF5AttrAsLitID ( dsID, 'molecule', qt%molecule )
    call GetHDF5Attribute ( dsID, 'instrumentModule', str )
    call GetModuleIndex ( str, qt%instrumentModule )
    call GetHDF5Attribute ( dsID, 'signal', str )
    call GetSignalIndex( str, qt%signal )
  end subroutine ReadAttributes_QuantityTemplate

  ! -----------------------------------  SetupNewQuantityTemplate  -----
  subroutine SetupNewQuantityTemplate ( qty, noInstances, noSurfs, &
    & noChans, coherent, stacked, regular, instanceLen, noCrossTrack, &
    & minorFrame, majorFrame, &
    & sharedVGrid, sharedFGrid, badValue, verticalCoordinate )

    use Intrinsic, only: L_Dimensionless

  ! Set up a new quantity template according to the user input.  This may
  ! be based on a previously supplied template (with possible
  ! modifications), or created from scratch.  The name isn't set here; the
  ! caller is expected to do it.

    ! Dummy arguments
    type (QuantityTemplate_T), intent(inout) :: qty ! Result

    integer, intent(in), optional :: noInstances
    integer, intent(in), optional :: noSurfs
    integer, intent(in), optional :: noChans
    logical, intent(in), optional :: coherent
    logical, intent(in), optional :: stacked
    logical, intent(in), optional :: regular
    integer, intent(in), optional :: instanceLen
    integer, intent(in), optional :: noCrossTrack
    logical, intent(in), optional :: minorFrame
    logical, intent(in), optional :: majorFrame
    logical, intent(in), optional :: sharedVGrid
    logical, intent(in), optional :: sharedFGrid
    real(rt), intent(in), optional :: badValue
    integer, intent(in), optional, value :: verticalCoordinate ! VALUE in case
                                        ! actual arg is qty%verticalCoordinate

    ! Local variables
    integer :: noSurfsToAllocate        ! For allocations
    integer :: noInstancesToAllocate    ! For allocations
    logical :: Verbose
    logical :: Verboser
    character(63) :: What

    ! Executable code
    if ( qty%name > 0 ) then
      call mygetString ( qty%name, what )
    else
      what = "qty"
    end if
    verbose = ( switchDetail(switches, 'qtmp' ) > -1 )
    verboser = ( switchDetail(switches, 'qtmp' ) > 0 )
    qty%quantityType = 0
    qty%noChans = 1
    qty%noInstances = 1
    qty%noSurfs = 1
    qty%noCrossTrack = 1
    qty%coherent = .true.
    qty%stacked = .true.
    qty%regular = .true.
    qty%minorFrame = .false.
    qty%majorFrame = .false.
    qty%logBasis = .false.
    qty%minValue = - huge ( 0.0_rt )
    qty%noInstancesLowerOverlap = 0
    qty%noInstancesUpperOverlap = 0
    qty%badValue = huge ( 0.0_rt )
    qty%unit = l_dimensionless
    qty%instanceLen = 1
    qty%verticalCoordinate = l_none
    qty%sharedVGrid = .false.
    qty%vGridIndex = 0
    qty%xGridIndex = 0
    qty%horizontalCoordinate = l_phiTan
    qty%hGridIndex = 0
    qty%instanceOffset = 0
    qty%frequencyCoordinate = l_none
    qty%sharedFGrid = .false.
    qty%fGridIndex = 0
    qty%lo = 0.0_rt
    qty%signal = 0
    qty%sideband = 0
    qty%instrumentModule = 0
    qty%radiometer = 0
    qty%reflector = 0
    qty%molecule = 0

    ! Now, see if the user asked for modifications to this
    if ( present (noChans) )            qty%noChans = noChans
    if ( present (noInstances) )        qty%noInstances = noInstances
    if ( present (noSurfs) )            qty%noSurfs = noSurfs
    if ( present (noCrossTrack) )       qty%noCrossTrack = noCrossTrack
    if ( present (regular) )            qty%regular = regular
    if ( present (minorFrame) )         qty%minorFrame = minorFrame
    if ( present (majorFrame) )         qty%majorFrame = majorFrame
    if ( present (sharedVGrid) )        qty%sharedVGrid = sharedVGrid
    if ( present (sharedFGrid) )        qty%sharedFGrid = sharedFGrid
    if ( present (badValue) )           qty%badValue = badValue
    if ( present (verticalCoordinate) ) qty%verticalCoordinate = verticalCoordinate

    if ( qty%minorFrame ) then
      if ( present(coherent) ) then
        if ( coherent ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Minor frame quantities must be incoherent" )
      end if
      qty%coherent = .FALSE.
      if ( present(stacked) ) then
        if ( stacked ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Minor frame quantities must be unstacked" )
      end if
      qty%stacked = .FALSE.
    else
      if ( present(coherent) ) qty%coherent = coherent
      if ( present(stacked) ) qty%stacked = stacked
    end if

    ! Now think about instanceLen
    if ( .not. qty%regular ) then
      if ( present(instanceLen) ) then
        qty%instanceLen = instanceLen
      else
        qty%instanceLen = 0
      end if
    else
      qty%instanceLen = qty%noSurfs * qty%noChans * qty%noCrossTrack
    end if

    ! Now we allocate all the arrays we're going to need if necessary
    if ( qty%coherent ) then
      noInstancesToAllocate = 1
    else
      noInstancesToAllocate = qty%noInstances
    end if

    if ( qty%stacked ) then
      noSurfsToAllocate = 1
    else
      noSurfsToAllocate = qty%noSurfs
    end if

    ! First the vertical coordinates
    if ( .not. qty%sharedVGrid ) then
      what = trim(what) // "%surfs"
      call allocate_test ( qty%surfs, qty%noSurfs, noInstancesToAllocate, &
        & what(1:len_trim(what)), ModuleName )
    end if

    ! Now the horizontal coordinates

    call allocate_test ( qty%crossAngles, qty%noCrossTrack, &
      & trim(what) // "%crossAngles", ModuleName )
    qty%crossAngles = 0.0 ! In case there actually is no xGrid
    call allocate_test ( qty%phi, noSurfsToAllocate, qty%noInstances, &
      & trim(what) // "%phi", ModuleName )
    ! Create GeodLat and Lon fields.
    call createGeolocationFields ( qty, noSurfsToAllocate, what )
    call allocate_test ( qty%time, noSurfsToAllocate, qty%noInstances, &
      & trim(what) // "%time", ModuleName )
    call allocate_test ( qty%solarTime, noSurfsToAllocate, qty%noInstances, &
      & trim(what) // "%solarTime", ModuleName )
    call allocate_test ( qty%solarZenith, noSurfsToAllocate, qty%noInstances, &
      & trim(what) // "%solarZenith", ModuleName )
    call allocate_test ( qty%losAngle, noSurfsToAllocate, qty%noInstances, &
      & trim(what) // "%losAngle", ModuleName )

    if ( .not. qty%regular ) then        !
      call allocate_test ( qty%surfIndex, qty%instanceLen, qty%noInstances, &
        & trim(what) // "%surfIndex", ModuleName )
      call allocate_test ( qty%chanIndex, qty%instanceLen, qty%noInstances, &
        & trim(what) // "%chanIndex", ModuleName )
    else
      nullify ( qty%surfIndex, qty%chanIndex )
    end if
    if ( verbose ) &
      & call outputnamedvalue ( 'grandTotalInstances (from setup)', qty%grandTotalInstances )

    if ( verboser )  call dump(qty, details=0, noL2CF=.true.)
  end subroutine SetupNewQuantityTemplate

  ! ---------------------------  WriteAttributes_QuantityTemplate  -----
  subroutine WriteAttributes_QuantityTemplate ( dsID, NAME, &
    & QT, NOL2CF )

    use MLSHDF5, only: MakeHDF5Attribute
    use MLSSignals_m, only: GetRadiometerName, GetModuleName, &
      & GetSignalName

    ! Arguments
    integer, intent(in) :: dsID
    character(len=*), intent(in) :: name
    type(QuantityTemplate_T), intent(in) :: qt
    logical, intent(in), optional :: NOL2CF  ! if TRUE => Skip l2cf-dependent
                                             ! attributes
    character (len=80) :: Str
    logical :: myNoL2CF
    ! Executable
    myNoL2CF = switchDetail(switches, 'nl2cf') > -1 ! .false.
    if ( present(NoL2CF) ) myNoL2CF = NoL2CF
    call myGetString ( qt%name, str, strip=.true. )
    call MakeHDF5Attribute ( dsID, name, 'TemplateName', str )
    if ( .not. myNoL2CF ) then
      call myGetString ( lit_indices(qt%quantityType), str, strip=.true. )
      call MakeHDF5Attribute ( dsID, name, 'tempQtyType', str )
    end if
    call MakeHDF5Attribute ( dsID, name, 'noInstances', qt%noInstances )
    call MakeHDF5Attribute ( dsID, name, 'noChans    ', qt%noChans     )
    call MakeHDF5Attribute ( dsID, name, 'noSurfs    ', qt%noSurfs     )
    call MakeHDF5Attribute ( dsID, name, 'coherent   ', qt%coherent     )
    call MakeHDF5Attribute ( dsID, name, 'stacked    ', qt%stacked     )
    call MakeHDF5Attribute ( dsID, name, 'regular    ', qt%regular     )
    call MakeHDF5Attribute ( dsID, name, 'logBasis    ', qt%logBasis     )
    call MakeHDF5Attribute ( dsID, name, 'minorFrame    ', qt%minorFrame     )
    call MakeHDF5Attribute ( dsID, name, 'badValue    ', qt%badValue     )
    if ( .not. myNoL2CF ) then
      call myGetString ( lit_indices(qt%unit), str, strip=.true. )
      call MakeHDF5Attribute ( dsID, name, 'tempQtyUnit', str )
    end if
    call MakeHDF5Attribute ( dsID, name, 'instanceLen    ', qt%instanceLen     )
    if ( .not. myNoL2CF ) then
      call myGetString ( lit_indices(qt%verticalCoordinate), str, strip=.true. )
      call MakeHDF5Attribute ( dsID, name, 'verticalCoordinate', str )
    end if
    if ( qt%radiometer /= 0 .and. .not. myNoL2CF ) then
      call GetRadiometerName ( qt%radiometer, str )
      call MakeHDF5Attribute ( dsID, name, 'radiometer', str )
    end if
    if ( qt%molecule + &
      &  qt%instrumentModule /= 0 .and. .not. myNoL2CF ) then
      if ( qt%molecule /= 0 ) then
        call myGetString ( lit_indices(qt%molecule), str, strip=.true. )
        call MakeHDF5Attribute ( dsID, name, 'molecule', str )
      end if
      if ( qt%instrumentModule /= 0 ) then
        call GetModuleName ( qt%instrumentModule, str )
        call MakeHDF5Attribute ( dsID, name, 'tempQtyInstrumentModule', str )
      end if
    end if
    if ( qt%signal /= 0 ) then
      call GetSignalName ( qt%signal, str )
      call MakeHDF5Attribute ( dsID, name, 'signal', str )
    end if
    ! Should we always write these, or only when specifically requested?
    if ( allocated(qt%surfs) ) &
      & call MakeHDF5Attribute ( dsID, name, 'surfs', qt%surfs(:,1) )
  end subroutine WriteAttributes_QuantityTemplate

  ! =====      Private procedures     ==================================

  ! ----------------------------  CheckIntegrity_QuantityTemplate  -----
  logical function CheckIntegrity_QuantityTemplate ( qty, noError )
    type (QuantityTemplate_T), intent(in) :: QTY
    logical, intent(in), optional :: NOERROR

    ! Local variables
    integer :: NOINSTANCESOR1           ! Test value
    integer :: NOSURFSOR1               ! Test value

    integer :: MESSAGETYPE
    character ( len=132 ) :: NAME

    ! Executable code
    messageType = MLSMSG_Error
    if ( present ( noError ) ) then
      if ( noError ) messageType = MLSMSG_Warning
    end if

    ! Now check the integrity of the template
    if ( qty%coherent ) then
      noInstancesOr1 = 1
    else
      noInstancesOr1 = qty%noInstances
    end if

    if ( qty%stacked ) then
      noSurfsOr1 = 1
    else
      noSurfsOr1 = qty%noSurfs
    end if

    if ( qty%name > 0 ) then
      call myGetString ( qty%name, name, strip=.true. )
    else
      name = '<no name>'
    end if

    CheckIntegrity_QuantityTemplate = .true.

    ! Check the instances / overlap stuff
    if ( qty%noInstances < 0 ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad noInstances for quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    if ( ( qty%noInstancesLowerOverlap < 0 ) .or. &
      &  ( qty%noInstancesLowerOverlap > qty%noInstances ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Inappropriate noInstancesLowerOverlap for quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    if ( ( qty%noInstancesUpperOverlap < 0 ) .or. &
      &  ( qty%noInstancesUpperOverlap > qty%noInstances ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Inappropriate noInstancesUpperOverlap for quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    if ( qty%noInstancesLowerOverlap + qty%noInstancesUpperOverlap > &
      & qty%noInstances ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Too much overlap for quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    ! Check the surfaces stuff
    if ( qty%noSurfs < 0 ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad noSurfs for quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    ! Check the channels stuff
    if ( qty%noChans < 0 ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad noChans for quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    ! Check the instanceLen
    if ( qty%regular .and. &
      & qty%InstanceLen /= qty%noSurfs * qty%noChans ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'The quantity template '//trim(name)//' does not have the right instanceLen' )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    ! Check the arrays are associated.  Note these have to be errors, as later
    ! tests will fail otherwise.
    if ( .not. allocated ( qty%surfs ) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'The quantity template '//trim(name)// ' does not have surfs associated' )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    if ( .not. allocated ( qty%phi ) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'The quantity template '//trim(name)// ' does not have phi allocated' )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    if ( .not. allocated ( qty%geodLat ) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'The quantity template '//trim(name)// ' does not have geodLat associated' )
      CheckIntegrity_QuantityTemplate = .false.
    end if
    if ( .not. allocated ( qty%lon ) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'The quantity template '//trim(name)// ' does not have lon associated' )
      CheckIntegrity_QuantityTemplate = .false.
    end if
    if ( .not. associated ( qty%time ) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'The quantity template '//trim(name)// ' does not have time associated' )
      CheckIntegrity_QuantityTemplate = .false.
    end if
    if ( .not. associated ( qty%solarTime ) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'The quantity template '//trim(name)// ' does not have solarTime associated' )
      CheckIntegrity_QuantityTemplate = .false.
    end if
    if ( .not. associated ( qty%solarZenith ) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'The quantity template '//trim(name)// ' does not have solarZenith associated' )
      CheckIntegrity_QuantityTemplate = .false.
    end if
    if ( .not. associated ( qty%losAngle ) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'The quantity template '//trim(name)// ' does not have losAngle associated' )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    ! Check the array lower bounds
    if ( any ( lbound ( qty%surfs ) /= 1 ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad lbound for surfs array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    if ( any ( lbound ( qty%phi ) /= 1 ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad lbound for phi array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    if ( any ( lbound ( qty%geodLat ) /= 1 ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad lbound for geodLat array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if
    if ( any ( lbound ( qty%lon ) /= (/1,1/) ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad lbound for lon array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if
    if ( any ( lbound ( qty%time ) /= 1 ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad lbound for time array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if
    if ( any ( lbound ( qty%solarTime ) /= 1 ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad lbound for solarTime array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if
    if ( any ( lbound ( qty%solarZenith ) /= 1 ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad lbound for solarZenith array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if
    if ( any ( lbound ( qty%losAngle ) /= 1 ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad lbound for losAngle array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    ! Check the array upper bounds
    if ( any ( ubound ( qty%surfs ) /= (/qty%noSurfs, noInstancesOr1/) ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad ubound for surfs array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    if ( any ( ubound ( qty%phi ) /= (/noSurfsOr1, qty%noInstances/) ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad ubound for phi array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    if ( any ( ubound ( qty%geodLat ) /= (/noSurfsOr1, qty%noInstances/) ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad ubound for geodLat array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if
    if ( any ( ubound ( qty%lon ) /= (/noSurfsOr1, qty%noInstances/) ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad ubound for lon array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if
    if ( any ( ubound ( qty%time ) /= (/noSurfsOr1, qty%noInstances/) ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad ubound for time array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if
    if ( any ( ubound ( qty%solarTime ) /= (/noSurfsOr1, qty%noInstances/) ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad ubound for solarTime array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if
    if ( any ( ubound ( qty%solarZenith ) /= (/noSurfsOr1, qty%noInstances/) ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad ubound for solarZenith array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if
    if ( any ( ubound ( qty%losAngle ) /= (/noSurfsOr1, qty%noInstances/) ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad ubound for losAngle array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    ! Check irregular stuff
    if ( .not. qty%regular ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'The quantity '//trim(name)//' appears to be irregular' )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    ! Could check channels stuff here, but not sure what to do.

  end function CheckIntegrity_QuantityTemplate

  ! -----------------------------------------  GetHDF5AttrAsLitID  -----
  subroutine GetHDF5AttrAsLitID ( dsID, attrName, LitID )
    ! Given a DS, File or GroupID, find the character-valued attribute
    ! for the attribute named attrName of the dataset name
    ! Look up its id in the lit indices and return that id as LitID
    use Intrinsic, only: First_Lit, Last_Auto_Lit
    use MLSHDF5, only: GetHDF5Attribute
    use String_Table, only: Add_Char, Lookup
    ! Args
    integer, intent(in)           :: dsID      ! dataset, file or group ID
    character(len=*), intent(in)  :: attrName  ! attribute name
    integer, intent(out)          :: LitID     ! where to find attr's value
    ! Internal variables
    logical :: found
    character(len=64) :: str
    integer :: strID
    ! litID = -1 ! meaning not found
    call GetHDF5Attribute ( dsID, attrname, str )
    call add_char( trim(str) )
    call lookup ( strID, found, caseless=.true., debug=0 )
    do litID=first_lit, Last_auto_lit
      if ( lit_indices(litID) == strID ) return
    end do
    litID = -1 ! Still not found
  end subroutine GetHDF5AttrAsLitID

  ! -----------------------------------------  GetHDF5AttrAsStrID  -----
  subroutine GetHDF5AttrAsStrID ( dsID, attrName, strID )
    ! Given a DS, File or GroupID, find the character-valued attribute
    ! for the attribute named attrName of the dataset name
    ! Look up its id in the string table and return that id as strID
    use MLSHDF5, only: GetHDF5Attribute
    use String_Table, only: Add_Char, Lookup
    use Tree_Types, only: Add_Char
    ! Args
    integer, intent(in)           :: dsID      ! dataset, file or group ID
    character(len=*), intent(in)  :: attrName  ! attribute name
    integer, intent(out)          :: strID     ! where to find attr's value
    ! Internal variables
    logical :: found
    character(len=64) :: str
    strID = -1 ! meaning not found
    call GetHDF5Attribute ( dsID, attrname, str )
    call add_char( trim(str) )
    call lookup ( strID, found, caseless=.true., debug=0 )
    
  end subroutine GetHDF5AttrAsStrID

  ! -------------------------------------------------  GetGeocLat  -----
  pure real(rt) function GetGeocLat ( Qty, Surf, Inst )
    ! Return the geocentric latitude in degrees.
    ! This is intended to be used only as a type-bound function with a
    ! binding named GeocLat.  It's OK to use a surface index without
    ! checking whether the quantity is stacked -- it's checked here.
    use Constants, only: Rad2Deg
    use Geometry, only: GeodToGeocLat
    class(quantityTemplate_t), intent(in) :: Qty
    integer, intent(in) :: Surf, Inst
    integer :: SurfOr1
    surfOr1 = merge(1,surf,qty%stacked)
    if ( qty%latitudeCoordinate == l_geodetic ) then
      getGeocLat = geodToGeocLat ( qty%geodlat(surfOr1, inst ) ) * rad2deg
    else
      getGeocLat = qty%geodlat(surfOr1, inst )
    end if
  end function GetGeocLat

  ! ------------------------------------------------  GetGeocLat3  -----
  pure real(rt) function GetGeocLat3 ( Qty, Surf, Inst, CrossIndex )
    ! Return the geocentric latitude in degrees.
    ! This is intended to be used only as a type-bound function with a
    ! binding named GeocLat.  It's OK to use a surface index without
    ! checking whether the quantity is stacked -- it's checked here.
    use Constants, only: Rad2Deg
    use Geometry, only: GeodToGeocLat
    class(quantityTemplate_t), intent(in) :: Qty
    integer, intent(in) :: Surf, Inst, CrossIndex
    integer :: SurfOr1
    surfOr1 = merge(1,surf,qty%stacked)
    if ( qty%latitudeCoordinate == l_geodetic ) then
      getGeocLat3 = geodToGeocLat ( qty%geodlat(surfOr1, &
        & inst + qty%noInstances * ( crossIndex - 1 ) ) ) * rad2deg
    else
      getGeocLat3 = qty%geodlat(surfOr1, inst + qty%noInstances * ( crossIndex - 1 ) )
    end if
  end function GetGeocLat3

  ! ----------------------------------------------------  GetLat3  -----
  pure real(rt) function GetLat3 ( Qty, Surf, Inst, CrossIndex )
    ! Return the geodetic latitude in degrees.
    ! This is intended to be used only as a type-bound function with a
    ! binding named GeodLat3.  It's OK to use a surface index without
    ! checking whether the quantity is stacked -- it's checked here.
    use Geometry, only: GeocToGeodLat
    class(quantityTemplate_t), intent(in) :: Qty
    integer, intent(in) :: Surf, Inst, CrossIndex
    integer :: SurfOr1
    surfOr1 = merge(1,surf,qty%stacked)
    if ( qty%latitudeCoordinate == l_geodetic ) then
      getLat3 = qty%geodlat(surfOr1, inst + qty%noInstances * ( crossIndex - 1 ) )
    else
      getLat3 = geocToGeodLat ( qty%geodlat(surfOr1, inst + qty%noInstances * ( crossIndex - 1 ) ) )
    end if
  end function GetLat3

  ! ----------------------------------------------------  GetLon3  -----
  pure real(rt) function GetLon3 ( Qty, Surf, Inst, CrossIndex )
    ! This is intended to be used only as a type-bound function with a
    ! binding named Lon3.  It's OK to use a surface index without
    ! checking whether the quantity is stacked -- it's checked here.
    class(quantityTemplate_t), intent(in) :: Qty
    integer, intent(in) :: Surf, Inst, CrossIndex
    integer :: SurfOr1
    surfOr1 = merge(1,surf,qty%stacked)
    getLon3 = qty%lon(surfOr1, inst + qty%noInstances * ( crossIndex - 1 ) )
  end function GetLon3

  ! ----------------------------------------------------  GetPhi3  -----
  pure real(rt) function GetPhi3 ( Qty, Surf, Inst, CrossIndex )
    ! This is intended to be used only as a type-bound function with a
    ! binding named Phi3.  It's OK to use a surface index without
    ! checking whether the quantity is stacked -- it's checked here.
    class(quantityTemplate_t), intent(in) :: Qty
    integer, intent(in) :: Surf, Inst, CrossIndex
    integer :: SurfOr1
    surfOr1 = merge(1,surf,qty%stacked)
    getPhi3 = qty%phi(surfOr1, inst + qty%noInstances * ( crossIndex - 1 ) )
  end function GetPhi3

  ! --------------------------------------------------  GetSurfs3  -----
  pure real(rt) function GetSurfs3 ( Qty, Surf, Inst, CrossIndex )
    ! This is intended to be used only as a type-bound function with a
    ! binding named Surfs3.  It's OK to use a surface index without
    ! checking whether the quantity is stacked -- it's checked here.
    class(quantityTemplate_t), intent(in) :: Qty
    integer, intent(in) :: Surf, Inst, CrossIndex
    integer :: SurfOr1
    surfOr1 = merge(1,surf,qty%stacked)
    getSurfs3 = qty%Surfs(surfOr1, inst + qty%noInstances * ( crossIndex - 1 ) )
  end function GetSurfs3

  ! --------------------------------------------  myDisplayString  -----
  subroutine myDisplayString ( index, advance, before )
    ! Given a string index, display the string or an error message
    use String_Table, only: Display_String, How_Many_Strings
    integer, intent(in) :: index
    character(len=*), intent(in), optional :: advance
    character(len=*), intent(in), optional :: before

    ! Executable code
    if ( present(before) ) call output ( before )
    if ( index < 1 ) then
      call output ( '(string index < 1)', advance=advance )
    else if ( index > how_many_strings() ) then ! How can an integer be a NaN ?????
      call output ( how_many_strings(), before='(string index > ', after=')', &
        & advance=advance )
    else
      call display_string ( index, advance )
    end if
  end subroutine myDisplayString

  ! ------------------------------------------------  myGetString  -----
  subroutine myGetString ( index, what, strip )
    ! Given a string index, Get the string or an error message
    use String_Table, only: Get_String, How_Many_Strings
    integer, intent(in)           :: index
    character(len=*), intent(out) :: what
    logical, intent(in), optional :: strip

    ! Executable code
    if ( index < 1 ) then
      what = '(string index < 1)'
    else if ( index > how_many_strings() ) then
      call writeIntsToChars( how_many_strings(), what )
      what = '(string index >' // trim(what) // ')'
    else
      call get_string ( index, what, strip )
    end if
  end subroutine myGetString

  ! --------------------------------------------  myValuesToField  -----
  ! This family of subroutines assigns from the values field
  ! explicitly to the template's own field
  ! Unless spread is TRUE, we assume there are exactly enough values
  subroutine myValuesToField_1d_dble ( TFIELD, VALUESNODE, SPREAD, TESTUNIT )
    double precision, dimension(:), intent(out)        :: tField ! Template's own field
    integer, intent(in)                      :: VALUESNODE   ! Tree node for values
    logical, intent(in)                      :: spread
    integer, intent(in) :: TestUnit                 ! Unit to use
    ! Internal variables
    integer :: k
    integer :: noValues
    integer, dimension(2) :: unitAsArray ! Unit for value given
    logical :: UNITSERROR               ! From expr
    real (rt), dimension(2) :: valueAsArray ! Value given
    ! Executable code
    if ( valuesNode < 1 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'illegal valuesNode in template field modify' )
    end if
    noValues = nsons(valuesNode) - 1
    if ( noValues < 1 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Too few values in template field modify' )
    end if
    do k = 1, noValues
      call expr_check ( subtree(k+1,valuesNode) , unitAsArray, valueAsArray, &
        & (/testUnit, PHYQ_Dimensionless/), unitsError )
      if ( unitsError ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'No units allowed for values in template field modify' )
      if ( spread ) then
        tField = valueAsArray(1)
        return
      else
        if ( k < size(tField) ) &
          & tField( k ) = valueAsArray(1)
      end if
    end do
  end subroutine myValuesToField_1d_dble

  subroutine myValuesToField_2d_real ( TFIELD, SHP, VALUESNODE, SPREAD, TESTUNIT )
    real, dimension(:,:), intent(out)        :: tField ! Template's own field
    integer, dimension(:), intent(in)        :: SHP
    integer, intent(in)                      :: VALUESNODE   ! Tree node for values
    logical, intent(in)                      :: spread
    integer, intent(in) :: TestUnit                 ! Unit to use
    ! Internal variables
    integer, dimension(2) :: indices
    integer :: k
    integer :: noValues
    integer, dimension(2) :: unitAsArray ! Unit for value given
    logical :: UNITSERROR               ! From expr
    real (rt), dimension(2) :: valueAsArray ! Value given

    ! Executable code
    if ( valuesNode < 1 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'illegal valuesNode in template field modify' )
    end if
    noValues = nsons(valuesNode) - 1
    if ( noValues < 1 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Too few values in template field modify' )
    end if
    ! call outputNamedValue( 'Number of values', NoValues )
    do k = 1, noValues
      call expr_check ( subtree(k+1,valuesNode) , unitAsArray, valueAsArray, &
        & (/testUnit, PHYQ_Dimensionless/), unitsError )
      if ( unitsError ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'No units allowed for values in template field modify' )
      ! call outputNamedValue( 'value', valueAsArray(1) )
      if ( spread ) then
        tField = valueAsArray(1)
        return
      else
        call rerank( k, shp, indices )
        ! call outputNamedValue( 'indices', indices )
        if ( any( indices < 1 ) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to rerank values in template field modify' )
        if ( all(indices < shp) ) &
          & tField( indices(1), indices(2) ) = valueAsArray(1)
      end if
    end do
  end subroutine myValuesToField_2d_real

  subroutine myValuesToField_2d_dble ( TFIELD, SHP, VALUESNODE, SPREAD, TESTUNIT )
    double precision, dimension(:,:), intent(out)        :: tField ! Template's own field
    integer, dimension(:), intent(in)        :: SHP
    integer, intent(in)                      :: VALUESNODE   ! Tree node for values
    logical, intent(in)                      :: spread
    integer, intent(in) :: TestUnit                 ! Unit to use
    ! Internal variables
    integer, dimension(2) :: indices
    integer :: k
    integer :: noValues
    integer, dimension(2) :: unitAsArray ! Unit for value given
    logical :: UNITSERROR               ! From expr
    real (rt), dimension(2) :: valueAsArray ! Value given

    ! Executable code
    if ( valuesNode < 1 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'illegal valuesNode in template field modify' )
    end if
    noValues = nsons(valuesNode) - 1
    if ( noValues < 1 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Too few values in template field modify' )
    end if
    ! call outputNamedValue( 'Number of values', NoValues )
    do k = 1, noValues
      call expr_check ( subtree(k+1,valuesNode) , unitAsArray, valueAsArray, &
        & (/testUnit, PHYQ_Dimensionless/), unitsError )
      if ( unitsError ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'No units allowed for values in template field modify' )
      ! call outputNamedValue( 'value', valueAsArray(1) )
      if ( spread ) then
        tField = valueAsArray(1)
        return
      else
        call rerank( k, shp, indices )
        ! call outputNamedValue( 'indices', indices )
        if ( any( indices < 1 ) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to rerank values in template field modify' )
        if ( all(indices <= shp) ) &
          & tField( indices(1), indices(2) ) = valueAsArray(1)
      end if
    end do
  end subroutine myValuesToField_2d_dble

  ! -------------------------------------------------  PutGeocLat  -----
  pure subroutine PutGeocLat ( Qty, Surf, Inst, Lat )
    ! Store a geocentric latitude as if the GeodLat component were a rank-3
    ! array with extents (surfs,insts,cross).  It's OK to use a surface index
    ! without checking whether the quantity is stacked -- it's checked here.
    use Geometry, only: GeocToGeodLat
    class(quantityTemplate_t), intent(inout) :: Qty
    integer, intent(in) :: Surf, Inst
    real(rt), intent(in) :: Lat ! geocentric, degrees
    integer :: SurfOr1
    surfOr1 = merge(1,surf,qty%stacked)
    if ( qty%latitudeCoordinate == l_geodetic ) then
      qty%geodlat(surfOr1, inst ) = geocToGeodLat ( lat )
    else
      qty%geodlat(surfOr1, inst ) = lat
    end if
  end subroutine PutGeocLat

  ! ------------------------------------------------  PutGeocLat3  -----
  pure subroutine PutGeocLat3 ( Qty, Surf, Inst, CrossIndex, Lat )
    ! Store a geocentric latitude as if the GeodLat component were a rank-3
    ! array with extents (surfs,insts,cross).  It's OK to use a surface index
    ! without checking whether the quantity is stacked -- it's checked here.
    use Geometry, only: GeocToGeodLat
    class(quantityTemplate_t), intent(inout) :: Qty
    integer, intent(in) :: Surf, Inst, CrossIndex
    real(rt), intent(in) :: Lat ! geocentric, degrees
    integer :: SurfOr1
    surfOr1 = merge(1,surf,qty%stacked)
    if ( qty%latitudeCoordinate == l_geodetic ) then
      qty%geodlat(surfOr1, inst + qty%noInstances * ( crossIndex - 1 ) ) = &
        & geocToGeodLat ( lat )
    else
      qty%geodlat(surfOr1, inst + qty%noInstances * ( crossIndex - 1 ) ) = lat
    end if
  end subroutine PutGeocLat3

  ! -----------------------------------------------------  PutLat  -----
  pure subroutine PutLat ( Qty, Surf, Inst, Lat )
    ! Store a geodetic latitude as if the GeodLat component were a rank-3 array
    ! with extents (surfs,insts,cross).  It's OK to use a surface index
    ! without checking whether the quantity is stacked -- it's checked here.
    use Constants, only: Rad2Deg
    use Geometry, only: GeocToGeodLat
    class(quantityTemplate_t), intent(inout) :: Qty
    integer, intent(in) :: Surf, Inst
    real(rt), intent(in) :: Lat ! degrees
    integer :: SurfOr1
    surfOr1 = merge(1,surf,qty%stacked)
    if ( qty%latitudeCoordinate == l_geodetic ) then
      qty%geodlat(surfOr1, inst ) = lat
    else
      qty%geodlat(surfOr1, inst ) = geocToGeodLat ( lat ) * rad2deg
    end if
  end subroutine PutLat

  ! ----------------------------------------------------  PutLat3  -----
  pure subroutine PutLat3 ( Qty, Surf, Inst, CrossIndex, Lat )
    ! Store a geodetic latitude as if the GeodLat component were a rank-3 array
    ! with extents (surfs,insts,cross).  It's OK to use a surface index
    ! without checking whether the quantity is stacked -- it's checked here.
    use Constants, only: Rad2Deg
    use Geometry, only: GeocToGeodLat
    class(quantityTemplate_t), intent(inout) :: Qty
    integer, intent(in) :: Surf, Inst, CrossIndex
    real(rt), intent(in) :: Lat ! degrees
    integer :: SurfOr1
    surfOr1 = merge(1,surf,qty%stacked)
    if ( qty%latitudeCoordinate == l_geodetic ) then
      qty%geodlat(surfOr1, inst + qty%noInstances * ( crossIndex - 1 ) ) = lat
    else
      qty%geodlat(surfOr1, inst + qty%noInstances * ( crossIndex - 1 ) ) = &
        & geocToGeodLat ( lat ) * rad2deg
    end if
  end subroutine PutLat3

  ! ----------------------------------------------------  PutLon3  -----
  pure subroutine PutLon3 ( Qty, Surf, Inst, CrossIndex, Lon )
    ! Store a Lonitude as if the Lon component were a rank-3 array
    ! with extents (surfs,insts,cross).  It's OK to use a surface index
    ! without checking whether the quantity is stacked -- it's checked here.
    class(quantityTemplate_t), intent(inout) :: Qty
    integer, intent(in) :: Surf, Inst, CrossIndex
    real(rt), intent(in) :: Lon
    integer :: SurfOr1
    surfOr1 = merge(1,surf,qty%stacked)
    qty%lon(surfOr1, inst + qty%noInstances * ( crossIndex - 1 ) ) = lon
  end subroutine PutLon3

  ! ----------------------------------------------------  PutPhi3  -----
  pure subroutine PutPhi3 ( Qty, Surf, Inst, CrossIndex, Phi )
    ! Store a Phi angle as if the Phi component were a rank-3 array
    ! with extents (surfs,insts,cross).  It's OK to use a surface index
    ! without checking whether the quantity is stacked -- it's checked here.
    class(quantityTemplate_t), intent(inout) :: Qty
    integer, intent(in) :: Surf, Inst, CrossIndex
    real(rt), intent(in) :: Phi
    integer :: SurfOr1
    surfOr1 = merge(1,surf,qty%stacked)
    qty%phi(surfOr1, inst + qty%noInstances * ( crossIndex - 1 ) ) = phi
  end subroutine PutPhi3

  ! ----------------------------------------------------  PutSurfs3  -----
  pure subroutine PutSurfs3 ( Qty, Surf, Inst, CrossIndex, SurfValue )
    ! Store a Surfs angle as if the Surfs component were a rank-3 array
    ! with extents (surfs,insts,cross).  It's OK to use a surface index
    ! without checking whether the quantity is stacked -- it's checked here.
    class(quantityTemplate_t), intent(inout) :: Qty
    integer, intent(in) :: Surf, Inst, CrossIndex
    real(rt), intent(in) :: SurfValue
    integer :: SurfOr1
    surfOr1 = merge(1,surf,qty%stacked)
    qty%surfs(surfOr1, inst + qty%noInstances * ( crossIndex - 1 ) ) = surfValue
  end subroutine PutSurfs3

!=============================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module QuantityTemplates
!=============================================================================

!
! $Log$
! Revision 2.121  2018/04/19 02:00:36  vsnyder
! Compute address for allocate/deallocate tracking.  Remove USE statements for
! unused names.
!
! Revision 2.120  2017/11/03 19:57:01  pwagner
! Most array gymnastics moved from MLSFillValues to HyperSlabs module
!
! Revision 2.119  2017/09/18 19:30:41  vsnyder
! Spiff the dump
!
! Revision 2.118  2017/01/25 17:19:26  pwagner
! Now require setting verboser for some Dumps
!
! Revision 2.117  2016/10/01 01:37:28  vsnyder
! Make QTM_Tree component of HGrid_t allocatable
!
! Revision 2.116  2016/08/30 20:27:24  vsnyder
! Add IsQTM type-bound function
!
! Revision 2.115  2016/08/23 00:42:43  vsnyder
! Components within or adjacent to the polygon are now within the QTM_Tree_t
! structure instead of the HGrid_t structure.
!
! Revision 2.114  2016/07/28 01:36:34  vsnyder
! Remove unreferenced USE and local variables
!
! Revision 2.113  2016/05/25 00:21:01  vsnyder
! Optionally allow different numbers of channels in QuantitiesAreCompatible.
! Check that HGrids are the same type in QuantitiesAreCompatible.
!
! Revision 2.112  2016/05/24 01:24:53  vsnyder
! Add checking for hGrid associated and of same type to QuantitiesAreCompatibie
!
! Revision 2.111  2016/05/19 23:28:07  pwagner
! Try harder not to crash on invalid strings; CopyQuantityTemplate takes optional arg dontDestroy
!
! Revision 2.110  2016/05/18 01:34:37  vsnyder
! HGridsDatabase.f90
!
! Revision 2.109  2016/05/12 15:22:46  pwagner
! Added GetHGridFromQuantity
!
! Revision 2.108  2016/05/04 18:34:16  pwagner
! Changed default unit to a valid phyq; some misunderstandings still remain however
!
! Revision 2.107  2015/09/25 02:12:25  vsnyder
! Add an optional verticalCoordinate argument to SetupNewQuantityTemplate.
! It has the VALUE attribute so qty%verticalCoordinate can be the actual
! argument without violating the standard, and without getting clobbered
! before it's needed.
!
! Revision 2.106  2015/09/23 22:39:09  vsnyder
! Add type-bound procedures to access and store latitude as either
! geocentric or geodetic latitude, except for the rank-2 GeodLat
! component.  Whether the latter is geocentric or geodetic still depends
! upon the value of the LatitudeCoordinate component, which users are
! expected to examine.
!
! Revision 2.105  2015/09/22 23:15:01  vsnyder
! Add 3D Phi and Surfs, spiff the dump
!
! Revision 2.104  2015/08/31 17:26:04  pwagner
! Fixed error in displaying qty%unit
!
! Revision 2.103  2015/08/26 01:08:17  vsnyder
! Yet more dump spiffing
!
! Revision 2.102  2015/08/25 18:36:03  vsnyder
! More dump spiffing
!
! Revision 2.101  2015/08/21 01:00:47  vsnyder
! Spiff a dump
!
! Revision 2.100  2015/07/31 20:42:22  pwagner
! Improved Dump
!
! Revision 2.99  2015/07/29 00:27:28  vsnyder
! Convert Phi from pointer to allocated
!
! Revision 2.98  2015/07/23 23:45:57  vsnyder
! qty%unit should be index in phyq_indices, not in lit_indices
!
! Revision 2.97  2015/06/04 03:13:16  vsnyder
! Make Surfs component of quantity template allocatable
!
! Revision 2.96  2015/06/03 23:09:00  pwagner
! Tried to prevent end-of-run crashes
!
! Revision 2.95  2015/06/02 23:53:00  vsnyder
! Add type-bound procedures to do rank-3 reference and update for latitude
! and longitude, instead of using rank-remapped pointers.
!
! Revision 2.94  2015/05/28 20:32:40  vsnyder
! Deallocate the correct component in DestroyGeolocationFields
!
! Revision 2.93  2015/05/27 22:41:46  vsnyder
! Move PointQuantityToHGrid here from HGridDatabase.  Eliminate shared
! HGrids.
!
! Revision 2.92  2015/05/01 02:09:28  vsnyder
! Spiff a dump
!
! Revision 2.91  2015/04/29 00:53:28  vsnyder
! Spiff the dump
!
! Revision 2.90  2015/04/07 02:51:50  vsnyder
! Remove CONTIGUOUS attribute from 2-D geolocation components.  Non-
! contiguous sections of them are targets in ConstructMajorFrameQuantity in
! ConstructQuantityTemplates.
!
! Revision 2.89  2015/03/28 01:40:45  vsnyder
! Added NoCrossTrack.  Added Unit (for values).  Added xGridIndex.  Added
! 1-d and 3-d geolocation quantities; made them contiguous.  Added
! LatitudeCoordinate (geocentric or geodetic).  Added CrossAngles.  Checked
! pointer components before copying them in CopyQuantityTemplate.  Added
! CreateGeolocationFields.  Fiddled with shared HGrids, but probably didn't
! improve anything.  Spiffed a dump, dump units, dump CrossAngles, dump
! LatitudeCorrdinate.  Adjusted QuantitiesAreCompatible to allow different
! types if requested.  Added NoCrossTrack to SetupNewQuantityTemplate.
! Added stuff to trace allocate/deallocate addresses -- some commented out
! because NAG build 1017 doesn't yet allow arrays as arguments to C_LOC.
!
! Revision 2.88  2014/10/29 23:04:29  vsnyder
! Specified units of several components
!
! Revision 2.87  2014/09/05 00:17:16  vsnyder
! More complete and accurate allocate/deallocate size tracking
!
! Revision 2.86  2014/08/19 00:28:33  vsnyder
! Make sure 'what' is always defined in DestroyQuantityTemplate
!
! Revision 2.85  2014/08/07 22:45:15  vsnyder
! Default HorizontalCoordinate to L_Phi_Tan instead of undefined
!
! Revision 2.84  2014/08/06 23:23:02  vsnyder
! Forgot to save from editor before commiting last time
!
! Revision 2.83  2014/08/06 23:22:28  vsnyder
! Combine several USE statements for the same module.  Remove declaration
! of unused parameter DEEBUG.
!
! Revision 2.82  2014/04/24 23:49:25  pwagner
! Added horizontalCoordinate component
!
! Revision 2.81  2014/03/20 01:39:47  vsnyder
! Unified types in Intrinsic
!
! Revision 2.80  2014/01/09 00:24:29  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.79  2013/12/12 01:57:17  vsnyder
! Change type of debug from logical to integer
!
! Revision 2.78  2013/09/19 23:31:08  vsnyder
! Use MyDisplayString more, add QuantitiesAreCompatible
!
! Revision 2.77  2013/08/16 02:27:04  vsnyder
! Declare RT named constant = R8, and use it for REAL components
!
! Revision 2.76  2013/08/12 23:47:25  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.75  2013/07/12 23:57:42  vsnyder
! Added Geolocation component
!
! Revision 2.74  2013/06/12 02:13:40  vsnyder
! Cruft removal
!
! Revision 2.73  2012/10/30 22:06:14  pwagner
! Fixed some obscure bugs when modifying templates
!
! Revision 2.72  2012/10/29 17:41:16  pwagner
! Attempted a more complete CopyQuantityTemplate
!
! Revision 2.71  2012/08/08 20:00:21  vsnyder
! Honest! I only changed some comments!
!
! Revision 2.70  2012/07/10 03:53:49  vsnyder
! Use DeepCopy
!
! Revision 2.69  2012/02/24 21:11:50  pwagner
! Include surfs when writing quantity attributes
!
! Revision 2.68  2012/02/23 00:08:35  vsnyder
! Don't dump molecule names if quantity type is not vmr
!
! Revision 2.67  2012/02/13 23:22:31  pwagner
! Print moleccule when dumping template
!
! Revision 2.66  2012/01/05 01:17:50  pwagner
! Added ReadAttributes; improved WriteAttributes
!
! Revision 2.65  2011/10/25 18:07:02  pwagner
! Added WriteAttributes to attach qty template attributes when writing datasets
!
! Revision 2.64  2011/05/09 17:26:03  pwagner
! Converted to using switchDetail
!
! Revision 2.63  2011/03/31 18:30:02  pwagner
! Corrected spelling in MODIFIABLEFIELDS
!
! Revision 2.62  2011/03/23 00:42:08  pwagner
! Tried to fix some of the more obvious bugs in ModifyQuantityTemplate_allocate
!
! Revision 2.61  2011/03/22 23:39:50  pwagner
! May now change both shape and values of qtytemplate field
!
! Revision 2.60  2011/03/15 22:43:52  pwagner
! Added ModifyQuantityTemplate; defaults to private
!
! Revision 2.59  2011/02/18 17:54:49  pwagner
! Prevented crashes when run w/o l2cf
!
! Revision 2.58  2010/09/25 01:16:35  vsnyder
! Add ChanInds component, some cannonball polishing
!
! Revision 2.57  2010/09/17 00:04:54  pwagner
! Workaround for obscure crashes when called from outside mlsl2
!
! Revision 2.56  2010/08/31 02:05:10  vsnyder
! Deallocate channels component in DestroyQuantityTemplateContents
!
! Revision 2.55  2010/02/04 23:08:00  vsnyder
! Remove USE or declaration for unused names
!
! Revision 2.54  2009/09/25 02:42:07  vsnyder
! Added badValue to SetupNewQuantityTemplate, dump channels field of
! quantity if it's associated.
!
! Revision 2.53  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.52  2008/09/30 22:28:03  vsnyder
! Remove AuxGrids -- didn't need them after all
!
! Revision 2.51  2008/06/06 01:54:08  vsnyder
! Aux grids have to be vGrids, not indices in vGridsDatabase, else clients
! will have to have the database.
! Make sure to deallocate the auxGrids.  Dump auxGrids.
!
! Revision 2.50  2008/06/05 02:05:53  vsnyder
! Added Aux grids
!
! Revision 2.49  2007/09/12 00:16:12  vsnyder
! Default initialize name component of QuantityTemplate_T to zero
!
! Revision 2.48  2007/03/23 00:11:52  pwagner
! qtmp switch now warns while destroying quantitytemplates
!
! Revision 2.47  2006/08/04 20:54:09  pwagner
! get_string for quantity name only if positive
!
! Revision 2.46  2006/08/04 01:54:16  vsnyder
! Use >0 instead of ==0 to test the name string
!
! Revision 2.45  2006/08/03 01:10:06  vsnyder
! Put l2cf names in leak track database
!
! Revision 2.44  2006/03/22 23:49:20  vsnyder
! Change the name of a dummy argument, add some comments
!
! Revision 2.43  2006/03/22 02:15:18  vsnyder
! Spiff up a dump
!
! Revision 2.42  2006/01/05 03:47:28  vsnyder
! Add some stuff to the dump
!
! Revision 2.41  2005/08/04 02:57:27  vsnyder
! Cannonball polishing
!
! Revision 2.40  2005/06/22 17:25:50  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.39  2004/08/26 18:54:49  pwagner
! qtmp switch now dumps in CreateQtyTemplateFromMLSCFInfo, not setup..
!
! Revision 2.38  2004/08/16 17:07:11  pwagner
! qtmp switch dumps quantity template after setup
!
! Revision 2.37  2004/05/01 04:07:44  vsnyder
! Rearranged some dumping stuff
!
! Revision 2.36  2004/04/15 20:51:51  pwagner
! Added DUMP_QUANTITY_TEMPLATES (found in l2/dumper)
!
! Revision 2.35  2004/01/24 01:02:43  livesey
! Added CopyQuantityTemplate
!
! Revision 2.34  2003/07/01 19:29:00  livesey
! Added grandTotalInstances
!
! Revision 2.33  2003/06/20 19:33:53  pwagner
! Quanities now share grids stored separately in databses
!
! Revision 2.32  2003/05/29 16:36:41  livesey
! Added the reflector item
!
! Revision 2.31  2003/01/14 21:35:53  vsnyder
! Add EPOCH and a comment about it in 'time' component
!
! Revision 2.30  2003/01/08 21:39:55  livesey
! Minor change in irregular quantity handling
!
! Revision 2.29  2002/11/27 01:06:26  livesey
! Better handling of major frame quantities
!
! Revision 2.28  2002/11/22 12:54:34  mjf
! Added nullify routine(s) to get round Sun's WS6 compiler not
! initialising derived type function results.
!
! Revision 2.27  2002/10/08 00:09:13  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.26  2002/09/24 21:36:42  livesey
! Added minValue
!
! Revision 2.25  2002/08/28 20:42:11  livesey
! Added InflateQuantityTemplateDatabase
!
! Revision 2.24  2002/07/22 03:26:05  livesey
! Added checkIntegrity
!
! Revision 2.23  2002/07/01 23:51:07  vsnyder
! Plug a memory leak
!
! Revision 2.22  2001/10/12 23:09:25  pwagner
! More debugging statements
!
! Revision 2.21  2001/10/03 17:42:27  pwagner
! reset DEEBUG to FALSE
!
! Revision 2.20  2001/10/02 23:12:50  pwagner
! More chi^2 fixes
!
! Revision 2.19  2001/09/17 21:59:26  livesey
! Removed allocate of frequencies, it's deferred to later in the code
!
! Revision 2.18  2001/09/13 19:59:43  pwagner
! Added majorframe as possible quantity type
!
! Revision 2.17  2001/07/31 23:39:12  dwu
! allocate and deallocate qty%frequencies
!
! Revision 2.16  2001/07/11 21:41:16  livesey
! Made quantityTemplateCounter public
!
! Revision 2.15  2001/07/02 17:25:30  livesey
! Some changes to comments, following walk through
!
! Revision 2.14  2001/05/23 20:38:35  livesey
! Updated a comment
!
! Revision 2.13  2001/04/23 23:52:16  livesey
! Sorry, should have put comment in one below.  Now has optional ignoreMinorFrame
! argument to DestroyQuantityTemplateDatabase
!
! Revision 2.12  2001/04/23 23:50:41  livesey
! *** empty log message ***
!
! Revision 2.11  2001/04/12 21:43:06  livesey
! Added sideband field
!
! Revision 2.10  2001/04/10 22:37:49  vsnyder
! Fix a type
!
! Revision 2.9  2001/03/24 00:31:12  pwagner
! USEs output in case we replace MLSMessage with output in additem..
!
! Revision 2.8  2001/03/17 02:23:18  livesey
! Added log basis field
!
! Revision 2.7  2001/03/15 20:20:59  vsnyder
! Correct the description of 'InstrumentModule'
!
! Revision 2.6  2001/03/02 01:34:03  livesey
! New signals stuff
!
! Revision 2.5  2001/02/23 17:47:01  livesey
! Nullified pointers.
!
! Revision 2.4  2001/02/14 00:12:34  livesey
! Removed firstIndexChannel
!
! Revision 2.3  2001/02/09 00:38:56  livesey
! Various changes
!
! Revision 2.2  2000/12/04 23:43:59  vsnyder
! Move more of addItemToDatabase into the include
!
! Revision 2.1  2000/10/13 00:00:37  vsnyder
! Moved from mlspgs/l2 to mlspgs/lib
!
! Revision 2.0  2000/09/05 18:57:04  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry
!
