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

  use ALLOCATE_DEALLOCATE, only: ALLOCATE_TEST, DEALLOCATE_TEST
  use DUMP_0, only: DUMP
  use EXPR_M, only: EXPR_CHECK
  use INTRINSIC, only: PHYQ_ANGLE, PHYQ_DIMENSIONLESS, PHYQ_FREQUENCY, &
    & PHYQ_TIME
  use MLSFILLVALUES, only: RERANK
  use MLSKINDS, only: RT => R8 ! RT is "kind of Real components of template"
  use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ALLOCATE, MLSMSG_DEALLOCATE, &
    & MLSMSG_ERROR, MLSMSG_WARNING
  use INTRINSIC, only: L_NONE, L_VMR, LIT_INDICES, PHYQ_INDICES
  use MLSFINDS, only: FINDFIRST
  use MLSSTRINGLISTS, only: SWITCHDETAIL
  use MLSSTRINGS, only: LOWERCASE, WRITEINTSTOCHARS
  use OUTPUT_M, only: OUTPUT, OUTPUTNAMEDVALUE
  use TOGGLES, only: SWITCHES
  use TREE, only: NSONS, SUBTREE

  implicit none
  private

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

  logical, parameter, private :: DEEBUG = .FALSE.           ! Usually FALSE

  ! Define some global parameters and data types.

  real(rt), parameter :: EPOCH = 1993.0 ! Starting point for time references

  type QuantityTemplate_T

    ! Some administrative stuff

    integer :: name = 0        ! Sub-rosa index of quantity name

    ! This integer is of an enumerated type describing what kind of
    ! quantity this is -- one of the l_lits of type t_quantityType
    ! in Init_Tables_Module, e.g. l_Temperature.

    integer :: quantityType

    ! The dimensions of this quantity

    integer :: noInstances     ! Number of horizontal instances in this quantity
    integer :: noSurfs         ! Number of surfaces per instance
    integer :: noChans         ! Number of channels

    ! Flags describing the quantity

    logical :: coherent        ! Do instances have same vertical coordinates?
    logical :: stacked         ! Are instances true vertical profiles?
    logical :: regular         ! Are all channels/heights represented

    ! This next one allows software using the vector quantities to be somewhat
    ! lazy and, for example, avoid interpolation.  Minor frame quantities are
    ! incoherent and unstacked, but may be regular or irregular.  However, not
    ! all incoherent unstacked quantities are minor frame quantities.

    logical :: minorFrame      ! Is this a minor frame quantity.
    logical :: majorFrame      ! Is this a major frame quantity.

    ! This one indicates whether log or linear interpolation should be used
    logical :: logBasis                 ! If set use log
    real(rt) :: minValue                ! Minimum value to consider if using log

    ! This information describes how much of the data is in the overlap
    ! regions if any.

    integer :: noInstancesLowerOverlap
    integer :: noInstancesUpperOverlap

    ! Misc. information
    real(rt) :: badValue      ! Value used to flag bad/missing data
    integer :: unit           ! Unit quantity is in when scaled as below,
                              ! an l_lit of the type t_units in Units.f90.

    ! For regular quantities the number of elements of each instance
    ! is simply noSurfs*noChans.  For irregular ones it is less, but it is
    ! constant from instance to instance; this is that number.
    integer :: instanceLen

    ! Vertical coordinate
    integer :: verticalCoordinate ! The vertical coordinate used.  These
                                  ! are l_lits of the type t_VGridCoord
                                  ! defined in Init_Tables_Module.
    logical :: sharedVGrid        ! Set if surfs is a pointer not a copy
    integer :: vGridIndex         ! Index of any vGrid used
    real(rt), dimension(:,:), pointer :: surfs => NULL()

    ! This is dimensioned (noSurfs,1) for coherent quantities and
    ! (noSurfs, noInstances) for incoherent ones.  Pretending the values are
    ! dimensioned (noChans, noSurfs, noInstances), the SURFS coordinate
    ! for the (:,i,j) values is surfs(i,1) for a coherent quantity or
    ! surfs(i,j) for an incoherent one.

    ! Horizontal coordinates
    logical :: sharedHGrid              ! Set if horiz coord is a pointer not a copy
    integer :: hGridIndex               ! Index of any hGrid used
    integer :: instanceOffset           ! Ind of 1st non overlapped instance in output
    integer :: grandTotalInstances      ! Total number of instances in destination output file
    ! for example MAF index, or profile index.

    ! First subscript values for GeoLocation component
    integer :: OrbitCoordinateIndex = 1 ! For spacecraft position
    integer :: LOSCoordinateIndex = 2   ! For line of sight
    real(rt), dimension(:,:,:), pointer :: Geolocation => NULL()

    ! Geolocation is dimensioned (*,1, noInstances) for stacked quantities and
    ! (*,noSurfs, noInstances) for unstacked ones.  The Geolocation coordinate
    ! for the (*,i,j) value is geolocation(*,1,j) for a stacked quantity and
    ! geolocation(*,i,j) for an unstacked one.  The "*" is taken from either
    ! the OrbitCoordinateIndex or LOSCoordinateIndex component.

    real(rt), dimension(:,:), pointer :: phi => NULL()

    ! Phi is dimensioned (1, noInstances) for stacked quantities and
    ! (noSurfs, noInstances) for unstacked ones.  The PHI coordinate for the
    ! (i,j) value is phi(1,j) for a stacked quantity and phi(i,j) for an
    ! unstacked one.  Phi is either taken from or derived from Geolocation.

    ! These other coordinates are dimensioned in the same manner as Phi:
    real(rt), dimension(:,:), pointer :: geodLat => NULL()
    real(rt), dimension(:,:), pointer :: lon => NULL()
    real(rt), dimension(:,:), pointer :: time => NULL() ! Seconds since EPOCH
    real(rt), dimension(:,:), pointer :: solarTime => NULL()
    real(rt), dimension(:,:), pointer :: solarZenith => NULL()
    real(rt), dimension(:,:), pointer :: losAngle => NULL()

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
    integer :: frequencyCoordinate      ! An enumerated type, e.g. FG_USBFreq
    real(rt) :: lo                      ! Local oscillator frequency, MHz
    logical :: sharedFGrid              ! Set of frequencies are a pointer not a copy
    integer :: sideband                 ! Associated sideband -1, 0, +1
    integer :: signal                   ! Index into signals database

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Some families of quantities require special additional information.
    ! This is given here if needed.

    integer :: instrumentModule ! Index in the Modules database in MLSSignals_m
    integer :: radiometer       ! For ptan etc., index into radiometers database
    integer :: reflector        ! For reflector efficiency etc. terms
    integer :: molecule ! What molecule does this refer to? (One of the l_...
                        ! lits of type t_molecule in Molecules.)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! For irregular quantities, we have these arrays to
    ! help us navigate around the quantity.

    integer, dimension(:,:), pointer :: surfIndex => NULL()
    integer, dimension(:,:), pointer :: chanIndex => NULL()
    ! These are actually dimensioned (instanceLen, noInstances)
  end type QuantityTemplate_T

  interface DUMP
    module procedure Dump_Quantity_Template, Dump_Quantity_Templates
  end interface

  interface MODIFYQUANTITYTEMPLATE
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

  public :: EPOCH, QuantityTemplate_T, RT
  public :: AddQuantityTemplateToDatabase, InflateQuantityTemplateDatabase
  public :: CheckIntegrity, CopyQuantityTemplate, &
    & DestroyQuantityTemplateContents, DestroyQuantityTemplateDatabase, &
    & Dump, ModifyQuantityTemplate, NullifyQuantityTemplate, &
    & QuantitiesAreCompatible, SetupNewQuantityTemplate, WriteAttributes
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

    ! Dummy arguments
    type (QuantityTemplate_T), dimension(:), pointer :: database
    type (QuantityTemplate_T), intent(in) :: item

    ! Local variables
    type (QuantityTemplate_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddQuantityTemplateToDatabase = newSize
  end function AddQuantityTemplateToDatabase

  ! ----------------------------  CopyQuantityTemplate  -----
  subroutine CopyQuantityTemplate ( Z, A )
    ! This routine does a 'deep' copy of a quantity template.
    ! We don't need to do if often as typically only a shallow copy
    ! is required.  Note that this also follows any 'links' to h/vGrids and expands
    ! them too.
    use DeepCopy_m, only: DeepCopy
    type (QuantityTemplate_T), intent(inout) :: Z
    type (QuantityTemplate_T), intent(in) :: A

    ! Executable code
    ! Destroy result
    call DestroyQuantityTemplateContents ( z )
    ! Setup result
    z%name = a%name
    call SetupNewQuantityTemplate ( z, a%noInstances, a%noSurfs, a%noChans, &
      & a%coherent, a%stacked, a%regular, a%instanceLen, a%minorFrame, &
      & a%majorFrame )
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
    z%sharedHGrid                  = a%sharedHGrid             
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

    ! Next, arrays
    z%surfs       = a%surfs
    z%phi         = a%phi
    z%geodLat     = a%geodLat
    z%lon         = a%lon
    z%time        = a%time
    z%solarTime   = a%solarTime
    z%solarZenith = a%solarZenith
    z%losAngle    = a%losAngle
    z%ChanInds    = a%ChanInds
    z%Channels    = a%Channels
    call deepCopy ( z%frequencies, a%frequencies )
    if ( .not. z%regular ) then
      z%surfIndex = a%surfIndex
      z%chanIndex = a%chanIndex
    end if

  end subroutine CopyQuantityTemplate

  ! ----------------------------  DestroyQuantityTemplateContents  -----
  subroutine DestroyQuantityTemplateContents ( qty )
    ! Dummy argument
    type (QuantityTemplate_T), intent(inout) :: QTY

    logical :: Verbose
    character(63) :: What

    ! Executable code
    verbose = ( switchDetail(switches, 'qtmp' ) > -1 )
    if ( verbose ) then
      call dump( qty )
      call output( 'About to destroy this quantity', advance='yes' )
    end if
    if ( qty%name == 0 ) then
      what = "qty"
    elseif ( verbose ) then
      call myGetString ( qty%name, what )
    end if

    if ( .not. qty%sharedVGrid ) then
      if ( verbose ) call output( 'About to deallocate surfs', advance='yes' )
      call deallocate_test ( qty%surfs, trim(what) // "%surfs", ModuleName )
    end if

    if ( .not. qty%sharedHGrid ) then
      if ( verbose ) call output( 'About to deallocate phi', advance='yes' )
      call deallocate_test ( qty%phi, trim(what) // '%phi', ModuleName )
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
    end if

    if ( .not. qty%sharedFGrid ) then
      if ( verbose ) call output( 'About to deallocate chanInds', advance='yes' )
      call deallocate_test ( qty%chanInds, trim(what) // "%chanInds", ModuleName )
      if ( verbose ) call output( 'About to deallocate channels', advance='yes' )
      call deallocate_test ( qty%channels, trim(what) // "%channels", ModuleName )
      if ( verbose ) call output( 'About to deallocate freqs', advance='yes' )
      call deallocate_test ( qty%frequencies, trim(what) // "%frequencies", ModuleName )
    end if

    if ( .not. qty%regular ) then
      if ( verbose ) call output( 'About to deallocate surfindex', advance='yes' )
      call deallocate_test ( qty%surfIndex, trim(what) // "%surfIndex", ModuleName )
      if ( verbose ) call output( 'About to deallocate chanindex', advance='yes' )
      call deallocate_test ( qty%chanIndex, trim(what) // '%chanIndex', ModuleName )
    end if

  end subroutine DestroyQuantityTemplateContents

  ! ----------------------------  DestroyQuantityTemplateDatabase  -----
  subroutine DestroyQuantityTemplateDatabase ( database )
    ! Dummy argument
    type (QuantityTemplate_T), dimension(:), pointer :: DATABASE

    ! Local variables
    integer :: qtyIndex, status
    logical :: verbose

    ! Executable code
    verbose = ( switchDetail(switches, 'qtmp' ) > -1 )
    if ( associated(database) ) then
      if ( verbose ) call outputNamedValue( 'size(qty db)', size ( database ) )
      do qtyIndex = 1, size ( database )
        call DestroyQuantityTemplateContents ( database(qtyIndex) )
      end do
      deallocate ( database, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate // "database" )
    end if
  end subroutine DestroyQuantityTemplateDatabase

  ! -------------------------------------  DUMP_QUANTITY_TEMPLATE  -----
  subroutine DUMP_QUANTITY_TEMPLATE ( QUANTITY_TEMPLATE, DETAILS, NOL2CF )

    use MLSSIGNALS_M, only: SIGNALS, DUMP, GETRADIOMETERNAME, GETMODULENAME
    use OUTPUT_M, only: NEWLINE
    use VGRIDSDATABASE, only: DUMP

    type(QuantityTemplate_T), intent(in) :: QUANTITY_TEMPLATE
    integer, intent(in), optional :: DETAILS ! <= 0 => Don't dump arrays
                                             ! >0   => Do signal, phi, surfs
                                             !         and frequency
                                             ! >1   => Dump all arrays
                                             ! Default 1
    logical, intent(in), optional :: NOL2CF  ! if TRUE => Don't dump L2-specific
                                             !  stuff
    integer :: MyDetails
    character (len=80) :: Str
    logical :: myNoL2CF
    myDetails = 1
    if ( present(details) ) myDetails = details
    myNoL2CF = switchDetail(switches, 'nl2cf') > -1 ! .false.
    if ( present(NoL2CF) ) myNoL2CF = NoL2CF
    call output ( ' Name = ' )
    call myDisplayString ( quantity_template%name )
    if ( .not. myNoL2CF ) then
      call output ( ' quantityType = ' )
      call myDisplayString ( lit_indices(quantity_template%quantityType), &
        & advance='yes' )
    end if
    call output ( quantity_template%noInstances, before='      NoInstances = ' )
    call output ( quantity_template%noChans,     before=' noChans = ' )
    call output ( quantity_template%noSurfs,     before=' NoSurfs = ', advance='yes' )
    call output ( '      ' )
    if ( .not. quantity_template%coherent ) call output ( 'in' )
    call output ( 'coherent ' )
    if ( .not. quantity_template%stacked ) call output ( 'non' )
    call output ( 'stacked ' )
    if ( .not. quantity_template%regular ) call output ( 'ir' )
    call output ( 'regular ' )
    if ( quantity_template%logBasis ) then
      call output ('log-')
    else
      call output ('linear-')
    end if
    call output ('basis ' )
    if ( .not. quantity_template%minorFrame ) call output ( 'non' )
    call output ( 'minorFrame', advance='yes' )
    call output ( '      NoInstancesLowerOverlap = ' )
    call output ( quantity_template%noInstancesLowerOverlap )
    call output ( ' NoInstancesUpperOverlap = ' )
    call output ( quantity_template%noInstancesUpperOverlap, advance='yes' )
    call output ( '      BadValue = ' )
    call output ( quantity_template%badValue )
    if ( .not. myNoL2CF ) then
      call output ( ' Unit = ' )
      call myDisplayString ( phyq_indices(quantity_template%unit) )
    end if
    call output ( ' InstanceLen = ' )
    call output ( quantity_template%InstanceLen, advance='yes' )
    call output ( '      sharedHGrid = ' )
    call output ( quantity_template%sharedHGrid, advance='no' )
    if ( quantity_template%sharedHGrid ) then
      call output ( ' hGridIndex = ' )
      call output ( quantity_template%hGridIndex, advance='yes' )
    else
      call newline
    end if
    call output ( '      sharedVGrid = ' )
    call output ( quantity_template%sharedVGrid, advance='no' )
    if ( quantity_template%sharedVGrid ) then
      call output ( ' vGridIndex = ' )
      call output ( quantity_template%vGridIndex )
    end if
    if ( .not. myNoL2CF ) call myDisplayString ( lit_indices(quantity_template%verticalCoordinate), &
      & before=' vertical coordinate = ', advance='yes' )
    call output ( '      sharedFGrid = ' )
    call output ( quantity_template%sharedFGrid, advance='no' )
    if ( quantity_template%sharedFGrid ) then
      call output ( ' fGridIndex = ' )
      call output ( quantity_template%fGridIndex, advance='yes' )
    else
      call newline
    end if
    if ( quantity_template%radiometer /= 0 .and. .not. myNoL2CF ) then
      call output ( '      Radiometer = ' )
      call GetRadiometerName ( quantity_template%radiometer, str )
      call output ( trim(str), advance='yes' )
    end if
    if ( .not. myNoL2CF ) then
      call output ( '     ' )
      if ( quantity_template%quantityType == l_vmr ) then
        call output ( ' Molecule = ' )
        call myDisplayString ( lit_indices(quantity_template%molecule) )
      end if
      if ( quantity_template%instrumentModule /= 0 ) then
        call output ( ' Instrument Module = ' )
        call GetModuleName ( quantity_template%instrumentModule, str )
        call output ( trim(str) )
      end if
      call newLine
    end if
    if ( myDetails > 0 ) then
      if ( quantity_template%signal /= 0 ) then
        call output ( '   Signal ' )
        call output ( quantity_template%signal )
        call output ( ':', advance='yes' )
        if ( associated(quantity_template%channels) ) then
          call dump ( signals(quantity_template%signal), &
            & otherChannels=quantity_template%channels )
        else
          call dump ( signals(quantity_template%signal) )
        end if
      end if
      if ( associated(quantity_template%phi) ) &
        & call dump ( quantity_template%phi,           '      Phi = ' )
      if ( associated(quantity_template%surfs) ) &
        & call dump ( quantity_template%surfs,         '      Surfs = ' )
      if ( myDetails > 1 ) then
        if ( associated(quantity_template%surfIndex) ) &
          & call dump ( quantity_template%surfIndex,   '      SurfIndex = ' )
        if ( associated(quantity_template%chanIndex) ) &
          & call dump ( quantity_template%chanIndex,   '      ChanIndex = ' )
        if ( associated(quantity_template%geodLat) ) &
          & call dump ( quantity_template%geodLat,     '      GeodLat = ' )
        if ( associated(quantity_template%lon) ) &
          & call dump ( quantity_template%lon,         '      Lon = ' )
        if ( associated(quantity_template%time) ) &
          & call dump ( quantity_template%time,        '      Time = ' )
        if ( associated(quantity_template%solarTime) ) &
          & call dump ( quantity_template%solarTime,   '      SolarTime = ' )
        if ( associated(quantity_template%solarZenith) ) &
          & call dump ( quantity_template%solarZenith, '      SolarZenith = ' )
        if ( associated(quantity_template%losAngle) ) &
          & call dump ( quantity_template%losAngle,    '      LosAngle = ' )
      end if
      if ( associated(quantity_template%frequencies)  .and. .not. myNoL2CF ) then
        call myDisplayString ( lit_indices(quantity_template%frequencyCoordinate), &
          & before='      FrequencyCoordinate = ', advance='yes' )
        call dump ( quantity_template%frequencies, ' Frequencies = ' )
      end if
    else
      if ( associated(quantity_template%frequencies)  .and. .not. myNoL2CF ) &
        & call myDisplayString ( lit_indices(quantity_template%frequencyCoordinate), &
          & before='      FrequencyCoordinate = ', advance='yes' )
    end if
  end subroutine DUMP_QUANTITY_TEMPLATE

  ! ------------------------------------  DUMP_QUANTITY_TEMPLATES  -----
  subroutine DUMP_QUANTITY_TEMPLATES ( QUANTITY_TEMPLATES, DETAILS, NOL2CF )

    type(QuantityTemplate_T), intent(in) :: QUANTITY_TEMPLATES(:)
    integer, intent(in), optional :: DETAILS ! <= 0 => Don't dump arrays
    !                                        ! >0   => Do dump arrays
    !                                        ! Default 1
    logical, intent(in), optional :: NOL2CF  ! if TRUE => Don't dump l2-specific

    integer :: I

    call output ( 'QUANTITY_TEMPLATES: SIZE = ' )
    call output ( size(quantity_templates), advance='yes' )
    do i = 1, size(quantity_templates)
      call output ( i, 4 )
      call output ( ':' )
      call dump_quantity_template ( quantity_templates(i), details, nol2cf )
    end do

  end subroutine DUMP_QUANTITY_TEMPLATES

  ! ------------------------------ InflateQuantityTemplateDatabase -----
  integer function InflateQuantityTemplateDatabase ( database, extra )
    ! Make a quantity template database bigger by extra
    ! Return index of first new element

    ! Dummy arguments
    type (QuantityTemplate_T), dimension(:), pointer :: DATABASE
    integer, intent(in) :: EXTRA

    ! Local variables
    type (QuantityTemplate_T), dimension(:), pointer :: TEMPDATABASE

    include "inflateDatabase.f9h"
    InflateQuantityTemplateDatabase = firstNewItem
  end function InflateQuantityTemplateDatabase

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
    endif
    select case(lowercase(field))
    case ( 'surfs' )
      if ( any(shape(z%surfs) /= shp) ) then
        call deallocate_test( z%surfs, 'template surfs', &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
        call allocate_test ( z%surfs, shp(1), shp(2), &
          & "template surfs", &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
      endif
      call myValuesToField( z%surfs, SHP, VALUESNODE, spread, PHYQ_Dimensionless )
    case ( 'phi' )
      if ( any(shape(z%phi) /= shp) ) then
        call deallocate_test( z%phi, 'template phi', &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
        call allocate_test ( z%phi, shp(1), shp(2), &
          & "template phi", &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
      endif
      call myValuesToField( z%phi, SHP, VALUESNODE, spread, phyq_angle )
    case ( 'geodlat' )
      if ( any(shape(z%geodlat) /= shp) ) then
        call deallocate_test( z%geodLat, 'template geodLat', &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
        call allocate_test ( z%geodLat, shp(1), shp(2), &
          & "template geodLat", &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
      endif
      call myValuesToField( z%geodLat, SHP, VALUESNODE, spread, phyq_angle )
    case ( 'lon' )
      if ( any(shape(z%lon) /= shp) ) then
        call deallocate_test( z%lon, 'template lon', &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
        call allocate_test ( z%lon, shp(1), shp(2), &
          & "template lon", &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
      endif
      call myValuesToField( z%lon, SHP, VALUESNODE, spread, phyq_angle )
    case ( 'time' )
      if ( any(shape(z%time) /= shp) ) then
        call deallocate_test( z%time, 'template time', &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
        call allocate_test ( z%time, shp(1), shp(2), &
          & "template time", &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
      endif
      call myValuesToField( z%time, SHP, VALUESNODE, spread, phyq_time )
    case ( 'solartime' )
      if ( any(shape(z%solartime) /= shp) ) then
        call deallocate_test( z%solartime, 'template solartime', &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
        call allocate_test ( z%solartime, shp(1), shp(2), &
          & "template solartime", &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
      endif
      call myValuesToField( z%solartime, SHP, VALUESNODE, spread, phyq_time )
    case ( 'solarzenith' )
      if ( any(shape(z%solarzenith) /= shp) ) then
        call deallocate_test( z%solarzenith, 'template solarzenith', &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
        call allocate_test ( z%solarzenith, shp(1), shp(2), &
          & "template solarzenith", &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
      endif
      call myValuesToField( z%solarzenith, SHP, VALUESNODE, spread, phyq_angle )
    case ( 'losangle' )
      if ( any(shape(z%losangle) /= shp) ) then
        call deallocate_test( z%losangle, 'template losangle', &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
        call allocate_test ( z%losangle, shp(1), shp(2), &
          & "template losangle", &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
      endif
      call myValuesToField( z%losangle, SHP, VALUESNODE, spread, phyq_angle )
    case ( 'frequencies' )
      if ( .not. associated(z%frequencies) ) then
        call allocate_test ( z%frequencies, shp(1), &
          & "template frequencies", &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
      elseif ( size(z%frequencies) /= shp(1) ) then
        call deallocate_test( z%frequencies, 'template frequencies', &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
        call allocate_test ( z%frequencies, shp(1), &
          & "template frequencies", &
          & ModuleName // 'ModifyQuantityTemplate_allocate' )
      endif
      call myValuesToField( z%frequencies, VALUESNODE, spread, phyq_frequency )
    case default
    end select
    if ( lowercase(field) == 'surfs' ) then
      z%NoSurfs = shp(1)
      if ( .not. z%coherent ) z%noInstances = shp(2)
    elseif ( lowercase(field) == 'frequencies' ) then
      z%NoChans = shp(1)
    else
      z%noInstances = shp(2)
      if ( .not. z%stacked ) z%NoSurfs = shp(1)
    endif
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
    endif
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
    elseif ( spread .and. shp(2) == 1 ) then
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
    endif
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
    endif
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

  ! ------------------------------------  QuantitiesAreCompatible  -----
  logical function QuantitiesAreCompatible ( Qty_1, Qty_2 )
    type(quantityTemplate_t), intent(in) :: Qty_1, Qty_2
    QuantitiesAreCompatible = &
      & qty_1%quantityType == qty_2%quantityType .and. &
      & qty_1%noInstances == qty_2%noInstances .and. &
      & qty_1%noSurfs == qty_2%noSurfs .and. &
      & qty_1%noChans == qty_2%noChans .and. &
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
    use INTRINSIC, only: GET_PHYQ
    use MLSHDF5, only: GETHDF5ATTRIBUTE
    use MLSSIGNALS_M, only: GETRADIOMETERINDEX, GETMODULEINDEX, &
      & GETSIGNALINDEX

    ! Arguments
    integer, intent(in) :: dsID
    type(QuantityTemplate_T), intent(inout) :: qt
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
    qt%unit = get_phyq( str )
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
    & noChans, coherent, stacked, regular, instanceLen, minorFrame, majorFrame, &
    & sharedVGrid, sharedHGrid, sharedFGrid, badValue )

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
    logical, intent(in), optional :: minorFrame
    logical, intent(in), optional :: majorFrame
    logical, intent(in), optional :: sharedVGrid
    logical, intent(in), optional :: sharedHGrid
    logical, intent(in), optional :: sharedFGrid
    real(rt), intent(in), optional :: badValue

    ! Local variables
    integer :: noSurfsToAllocate        ! For allocations
    integer :: noInstancesToAllocate    ! For allocations

    character(63) :: What

    ! Executable code
    if ( qty%name > 0 ) then
      call mygetString ( qty%name, what )
    else
      what = "qty"
    end if

    qty%quantityType = 0
    qty%noChans = 1
    qty%noInstances = 1
    qty%noSurfs = 1
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
    qty%unit = 0
    qty%instanceLen = 1
    qty%verticalCoordinate = l_none
    qty%sharedVGrid = .false.
    qty%vGridIndex = 0
    qty%sharedHGrid = .false.
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
    if ( present (noChans) )      qty%noChans = noChans
    if ( present (noInstances) )  qty%noInstances = noInstances
    if ( present (noSurfs) )      qty%noSurfs = noSurfs
    if ( present (regular) )      qty%regular = regular
    if ( present (minorFrame) )   qty%minorFrame = minorFrame
    if ( present (majorFrame) )   qty%majorFrame = majorFrame
    if ( present (sharedVGrid) )  qty%sharedVGrid = sharedVGrid
    if ( present (sharedHGrid) )  qty%sharedHGrid = sharedHGrid
    if ( present (sharedFGrid) )  qty%sharedFGrid = sharedFGrid
    if ( present (badValue) )     qty%badValue = badValue

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
      qty%instanceLen = qty%noSurfs * qty%noChans
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
    if ( qty%sharedVGrid ) then
      nullify ( qty%surfs )
    else
      call allocate_test ( qty%surfs, qty%noSurfs, noInstancesToAllocate, &
        & trim(what) // "%surfs", ModuleName )
    end if

    ! Now the horizontal coordinates

    if ( qty%sharedHGrid ) then
      nullify ( qty%phi, qty%geodLat, qty%lon, qty%time, qty%solarTime, &
        & qty%solarZenith, qty%losAngle )
    else
      call allocate_test ( qty%phi, noSurfsToAllocate, qty%noInstances, &
        & trim(what) // "%phi", ModuleName )
      call allocate_test ( qty%geodLat, noSurfsToAllocate, qty%noInstances, &
        & trim(what) // "%geodLat", ModuleName )
      call allocate_test ( qty%lon, noSurfsToAllocate, qty%noInstances, &
        & trim(what) // "%lon", ModuleName )
      call allocate_test ( qty%time, noSurfsToAllocate, qty%noInstances, &
        & trim(what) // "%time", ModuleName )
      call allocate_test ( qty%solarTime, noSurfsToAllocate, qty%noInstances, &
        & trim(what) // "%solarTime", ModuleName )
      call allocate_test ( qty%solarZenith, noSurfsToAllocate, qty%noInstances, &
        & trim(what) // "%solarZenith", ModuleName )
      call allocate_test ( qty%losAngle, noSurfsToAllocate, qty%noInstances, &
        & trim(what) // "%losAngle", ModuleName )
    end if

    if ( .not. qty%regular ) then        !
      call allocate_test ( qty%surfIndex, qty%instanceLen, qty%noInstances, &
        & trim(what) // "%surfIndex", ModuleName )
      call allocate_test ( qty%chanIndex, qty%instanceLen, qty%noInstances, &
        & trim(what) // "%chanIndex", ModuleName )
    else
      nullify ( qty%surfIndex, qty%chanIndex )
    end if

    ! if ( switchDetail(switches, 'qtmp') > -1 ) call dump(qty, details=0, noL2CF=.true.)
  end subroutine SetupNewQuantityTemplate

  ! ---------------------------  WriteAttributes_QuantityTemplate  -----
  subroutine WriteAttributes_QuantityTemplate ( dsID, NAME, &
    & QT, NOL2CF )

    use MLSHDF5, only: MAKEHDF5ATTRIBUTE
    use MLSSIGNALS_M, only: GETRADIOMETERNAME, GETMODULENAME, &
      & GETSIGNALNAME

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
      call myGetString ( phyq_indices(qt%unit), str, strip=.true. )
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
    if ( associated(qt%surfs) ) &
      & call MakeHDF5Attribute ( dsID, name, 'surfs', qt%surfs(:,1) )
  end subroutine WriteAttributes_QuantityTemplate

  ! --------- Private procedures ---
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
    if ( .not. associated ( qty%surfs ) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'The quantity template '//trim(name)// ' does not have surfs associated' )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    if ( .not. associated ( qty%phi ) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'The quantity template '//trim(name)// ' does not have phi associated' )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    if ( .not. associated ( qty%geodLat ) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'The quantity template '//trim(name)// ' does not have geodLat associated' )
      CheckIntegrity_QuantityTemplate = .false.
    end if
    if ( .not. associated ( qty%lon ) ) then
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
    use INTRINSIC, only: FIRST_LIT, LAST_AUTO_LIT
    use MLSHDF5, only: GETHDF5ATTRIBUTE
    use STRING_TABLE, only: ADD_CHAR, LOOKUP
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
    call lookup ( strID, found, caseless=.true., debug=.false. )
    do litID=first_lit, Last_auto_lit
      if ( lit_indices(litID) == strID ) return
    enddo
    litID = -1 ! Still not found
  end subroutine GetHDF5AttrAsLitID

  ! -----------------------------------------  GetHDF5AttrAsStrID  -----
  subroutine GetHDF5AttrAsStrID ( dsID, attrName, strID )
    ! Given a DS, File or GroupID, find the character-valued attribute
    ! for the attribute named attrName of the dataset name
    ! Look up its id in the string table and return that id as strID
    use MLSHDF5, only: GETHDF5ATTRIBUTE
    use STRING_TABLE, only: ADD_CHAR, LOOKUP
    use TREE_TYPES, only: ADD_CHAR
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
    call lookup ( strID, found, caseless=.true., debug=.false. )
    
  end subroutine GetHDF5AttrAsStrID

  ! --------------------------------------------  myDisplayString  -----
  subroutine myDisplayString ( index, advance, before )
    ! Given a string index, display the string or an error message
    use String_Table, only: DISPLAY_STRING, HOW_MANY_STRINGS
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
    use STRING_TABLE, only: GET_STRING, HOW_MANY_STRINGS
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
    endif
    noValues = nsons(valuesNode) - 1
    if ( noValues < 1 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Too few values in template field modify' )
    endif
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
      endif
    enddo
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
    endif
    noValues = nsons(valuesNode) - 1
    if ( noValues < 1 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Too few values in template field modify' )
    endif
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
      endif
    enddo
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
    endif
    noValues = nsons(valuesNode) - 1
    if ( noValues < 1 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Too few values in template field modify' )
    endif
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
      endif
    enddo
  end subroutine myValuesToField_2d_dble

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
