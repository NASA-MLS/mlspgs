! Copyright (c) 1999, California Institute of Technology. ALL RIGHTS RESERVED.
! U.S. Government sponsorship under NASA Contract NAS7407 is acknowledged.

!=============================================================================
module QuantityTemplates         ! Quantities within vectors
!=============================================================================

  ! This module defines the `quantities' that make up vectors and their
  ! template information.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use DUMP_0, only: DUMP
  use MLSCommon, only: NameLen, R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error, MLSMSG_Warning
  use Intrinsic, only: L_None, LIT_INDICES, PHYQ_INDICES
  use Output_m, only: Output
  use String_Table, only: DISPLAY_STRING, Get_String

  implicit none
  public

  !------------------------------- RCS Ident Info ------------------------------
  character(len=130), private :: id = & 
       "$Id$"
  character(len=*), parameter, private :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here 
  !-----------------------------------------------------------------------------

  logical, parameter, private :: DEEBUG = .FALSE.           ! Usually FALSE

  ! Define some global parameters and data types.

  real(r8), parameter :: EPOCH = 1993.0 ! Starting point for time references

  type QuantityTemplate_T

    ! Some administrative stuff

    integer :: name            ! Sub-rosa index of quantity name

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
    real(r8) :: minValue                ! Minimum value to consider if using log

    ! This information describes how much of the data is in the overlap
    ! regions if any.

    integer :: noInstancesLowerOverlap
    integer :: noInstancesUpperOverlap

    ! Misc. information
    real(r8) :: badValue      ! Value used to flag bad/missing data
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
    real(r8), dimension(:,:), pointer :: surfs => NULL()

    ! This is dimensioned (noSurfs,1) for coherent quantities and
    ! (noSurfs, noInstances) for incoherent ones.

    ! Horizontal coordinates
    logical :: sharedHGrid              ! Set if horiz coords a pointer not a copy
    integer :: hGridIndex               ! Index of any hGrid used
    integer :: instanceOffset           ! Ind of 1st non overlapped instance in output
    integer :: grandTotalInstances      ! Total number of instances in destination output file
    ! for example MAF index, or profile index.
    real(r8), dimension(:,:), pointer :: phi => NULL()

    ! This is dimensioned (1, noInstances) for stacked quantities and
    ! (noSurfs, noInstances) for unstacked ones.
    
    ! These other coordinates are dimensioned in the same manner:
    real(r8), dimension(:,:), pointer :: geodLat => NULL()
    real(r8), dimension(:,:), pointer :: lon => NULL()
    real(r8), dimension(:,:), pointer :: time => NULL() ! Seconds since EPOCH
    real(r8), dimension(:,:), pointer :: solarTime => NULL()
    real(r8), dimension(:,:), pointer :: solarZenith => NULL()
    real(r8), dimension(:,:), pointer :: losAngle => NULL()

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! For quantities containing `channels' the following information may or
    ! may not be useful.

    ! Some quantities are on abritrary freqency grids; these quantities refer
    ! to those.
    integer :: frequencyCoordinate ! An enumerated type, e.g. FG_USBFreq
    logical :: sharedFGrid              ! Set of frequencies are a pointer not a copy
    integer :: fGridIndex               ! Index of any fGrid Index used
    real(r8), dimension(:), pointer :: frequencies => NULL() ! List of frequencies
                                                   ! (noChans)
    real(r8) :: lo                      ! Local oscillator
    integer :: signal                   ! Index into signals database
    integer :: sideband                 ! Associated sideband -1, 0, +1

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Some families of quantities require special additional information.
    ! This is given here if needed.

    integer :: instrumentModule ! Index in the Modules database in MLSSignals_m
    integer :: radiometer       ! For ptan etc., index into radiometers database
    integer :: reflector        ! For reflector efficiency etc. terms
    integer :: molecule ! What molecule does this refer to? (One of the l_...
                        ! lits of type t_molecule in Molecules.)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! For irregular quantities, instead of using the we have these arrays to
    ! help us navigate around the quantity.

    integer, dimension(:,:), pointer :: surfIndex => NULL()
    integer, dimension(:,:), pointer :: chanIndex => NULL()
    ! These are actually dimensioned (instanceLen, noInstances)
  end type QuantityTemplate_T

  ! Local procedures
  interface CheckIntegrity
    module procedure CheckIntegrity_QuantityTemplate
  end interface

  private :: CheckIntegrity_QuantityTemplate
  
contains ! =====     Public Procedures     =============================

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

  ! ---------------------------- CheckIntegrity_QuantityTemplate -------
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

    if ( qty%name /= 0 ) then
      call get_string ( qty%name, name, strip=.true. )
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

  ! ----------------------------  CopyQuantityTemplate  -----
  subroutine CopyQuantityTemplate ( Z, A )
    ! This routine does a 'deep' copy of a quantity template.
    ! We don't need to do if often as typically only a shallow copy
    ! is required.  Note that this also follows any 'links' to h/vGrids and expands
    ! them too.
    type (QuantityTemplate_T), intent(inout) :: Z
    type (QuantityTemplate_T), intent(in) :: A

    ! Executable code
    ! Destroy result
    call DestroyQuantityTemplateContents ( z )
    ! Setup result
    call SetupNewquantityTemplate ( z, a%noInstances, a%noSurfs, a%noChans, &
      & a%coherent, a%stacked, a%regular, a%instanceLen, a%minorFrame, a%majorFrame )
    ! Copy each other component -- tedious, but can't do a shallow copy
    ! as would lose newly allocated arrays
    z%name = a%name
    z%quantityType = a%quantityType
    z%logBasis = a%logBasis
    z%minValue = a%minValue
    z%noInstancesLowerOverlap = a%noInstancesLowerOverlap
    z%noInstancesUpperOverlap = a%noInstancesUpperOverlap
    z%badValue = a%badValue
    z%unit = a%unit
    z%verticalCoordinate = a%verticalCoordinate
    z%surfs = a%surfs
    z%instanceOffset = a%instanceOffset
    z%grandTotalInstances = a%grandTotalInstances
    z%phi = a%phi
    z%geodLat = a%geodLat
    z%lon = a%lon
    z%time = a%time
    z%solarTime = a%solarTime
    z%solarZenith = a%solarZenith
    z%losAngle = a%losAngle
    z%frequencyCoordinate = a%frequencyCoordinate
    if ( associated ( a%frequencies ) ) then
      call Allocate_test ( z%frequencies, size(a%frequencies), 'z%frequencies', ModuleName )
      z%frequencies = a%frequencies
    end if
    z%lo = a%lo
    z%signal = a%signal
    z%sideband = a%sideband
    z%instrumentModule = a%instrumentModule
    z%radiometer = a%radiometer
    z%reflector = a%reflector
    z%molecule = a%molecule
    if ( .not. z%regular ) then
      z%surfIndex = a%surfIndex
      z%chanIndex = a%chanIndex
    end if
   
  end subroutine CopyQuantityTemplate

  ! ----------------------------  DestroyQuantityTemplateContents  -----
  subroutine DestroyQuantityTemplateContents ( qty )
    ! Dummy argument
    type (QuantityTemplate_T), intent(inout) :: QTY

    ! Executable code
    if ( .not. qty%sharedVGrid ) then
      call deallocate_test ( qty%surfs, "qty%surfs", ModuleName )
    end if

    if ( .not. qty%sharedHGrid ) then
      call deallocate_test ( qty%phi, "qty%phi", ModuleName )
      call deallocate_test ( qty%geodLat, "qty%geodLat", ModuleName )
      call deallocate_test ( qty%lon, "qty%lon", ModuleName )
      call deallocate_test ( qty%time, "qty%time", ModuleName )
      call deallocate_test ( qty%solarTime, "qty%solarTime", ModuleName )
      call deallocate_test ( qty%solarZenith, "qty%solarZenith", ModuleName )
      call deallocate_test ( qty%losAngle, "qty%losAngle", ModuleName )
    end if
    
    if ( .not. qty%sharedFGrid ) then
      call deallocate_test ( qty%frequencies, "qty%frequencies", ModuleName )
    end if

    if ( .not. qty%regular) then
      call deallocate_test ( qty%surfIndex, "qty%surfIndex", ModuleName )
      call deallocate_test ( qty%chanIndex, "qty%chanIndex", ModuleName )
    end if
  end subroutine DestroyQuantityTemplateContents

  ! ----------------------------  DestroyQuantityTemplateDatabase  -----
  subroutine DestroyQuantityTemplateDatabase ( database )
    ! Dummy argument
    type (QuantityTemplate_T), dimension(:), pointer :: DATABASE

    ! Local variables
    integer :: qtyIndex, status

    if ( associated(database) ) then
      do qtyIndex = 1, size ( database )
        call DestroyQuantityTemplateContents ( database(qtyIndex) )
      end do
      deallocate ( database, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate // "database" )
    end if
  end subroutine DestroyQuantityTemplateDatabase

  ! ------------------------------------  DUMP_QUANTITY_TEMPLATES  -----
  subroutine DUMP_QUANTITY_TEMPLATES ( QUANTITY_TEMPLATES, DETAILS, NOL2CF )

    use MLSSignals_m, only: signals, DUMP, GetRadiometerName, GetModuleName

    type(QuantityTemplate_T), intent(in) :: QUANTITY_TEMPLATES(:)
    integer, intent(in), optional :: DETAILS ! <= 0 => Don't dump arrays
    !                                        ! >0   => Do dump arrays
    !                                        ! Default 1
    logical, intent(in), optional :: NOL2CF  ! if TRUE => Don't dump l2-specific
    integer :: I, MyDetails                   !  stuff
    character (len=80) :: Str
    logical :: myNoL2CF
    myDetails = 1
    if ( present(details) ) myDetails = details
    myNoL2CF = .false.
    if ( present(NoL2CF) ) myNoL2CF = NoL2CF
    call output ( 'QUANTITY_TEMPLATES: SIZE = ' )
    call output ( size(quantity_templates), advance='yes' )
    do i = 1, size(quantity_templates)
      call output ( i, 4 )
      call output ( ': Name = ' )
      call myDisplayString ( quantity_templates(i)%name )
      if ( .not. myNoL2CF ) then
      call output ( ' quantityType = ' )
      call myDisplayString ( lit_indices(quantity_templates(i)%quantityType), &
        & advance='yes' )
      endif
      call output ( '      NoInstances = ' )
      call output ( quantity_templates(i)%noInstances )
      call output ( ' NoSurfs = ' )
      call output ( quantity_templates(i)%noSurfs )
      call output ( ' noChans = ' )
      call output ( quantity_templates(i)%noChans, advance='yes' )
      call output ( '      ' )
      if ( .not. quantity_templates(i)%coherent ) call output ( 'in' )
      call output ( 'coherent ' )
      if ( .not. quantity_templates(i)%stacked ) call output ( 'non' )
      call output ( 'stacked ' )
      if ( .not. quantity_templates(i)%regular ) call output ( 'ir' )
      call output ( 'regular ' )
      if ( quantity_templates(i)%logBasis ) then
        call output ('log-')
      else
        call output ('linear-')
      endif
      call output ('basis ' )  
      if ( .not. quantity_templates(i)%minorFrame ) call output ( 'non' )
      call output ( 'minorFrame', advance='yes' )
      call output ( '      NoInstancesLowerOverlap = ' )
      call output ( quantity_templates(i)%noInstancesLowerOverlap )
      call output ( ' NoInstancesUpperOverlap = ' )
      call output ( quantity_templates(i)%noInstancesUpperOverlap, advance='yes' )
      call output ( '      BadValue = ' )
      call output ( quantity_templates(i)%badValue )
      if ( .not. myNoL2CF ) then
      call output ( ' Unit = ' )
      call myDisplayString ( phyq_indices(quantity_templates(i)%unit) )
      endif
      call output ( ' InstanceLen = ' )
      call output ( quantity_templates(i)%InstanceLen, advance='yes' )
      if ( myDetails < 0 ) then
        call dump ( quantity_templates(i)%surfs, '  Surfs = ' )
        call dump ( quantity_templates(i)%phi, '      Phi = ' )
        call dump ( quantity_templates(i)%geodLat, '      GeodLat = ' )
        call dump ( quantity_templates(i)%lon, '      Lon = ' )
        call dump ( quantity_templates(i)%time, '      Time = ' )
        call dump ( quantity_templates(i)%solarTime, '      SolarTime = ' )
        call dump ( quantity_templates(i)%solarZenith, '      SolarZenith = ' )
        call dump ( quantity_templates(i)%losAngle, '      LosAngle = ' )
        if ( associated(quantity_templates(i)%frequencies) ) then
          call output ( '      FrequencyCoordinate = ' )
          call output ( quantity_templates(i)%frequencyCoordinate )
          call dump ( quantity_templates(i)%frequencies, ' Frequencies = ' )
        end if
      end if
      if ( quantity_templates(i)%radiometer /= 0 .and. .not. myNoL2CF ) then
        call output ( '      Radiometer = ' )
        call GetRadiometerName ( quantity_templates(i)%radiometer, str )
        call output ( trim(str), advance='yes' )
      end if
      if ( quantity_templates(i)%molecule + &
        &  quantity_templates(i)%instrumentModule /= 0 .and. .not. myNoL2CF ) then
        call output ( '     ' )
        if ( quantity_templates(i)%molecule /= 0 ) then
          call output ( ' Molecule = ' )
          call myDisplayString ( lit_indices(quantity_templates(i)%molecule) )
        end if
        if ( quantity_templates(i)%instrumentModule /= 0 ) then
          call output ( ' Instrument Module = ' )
          call GetModuleName ( quantity_templates(i)%instrumentModule, str )
          call output ( trim(str) )
        end if
        call output ( '', advance = 'yes')
      end if
      if ( myDetails > 0 ) then
        if ( quantity_templates(i)%signal /= 0 ) then
          call dump ( signals( (/ quantity_templates(i)%signal /) ) )
        end if
        if ( quantity_templates(i)%radiometer + &
          &  quantity_templates(i)%molecule /= 0 ) &
          &  call output ( '', advance='yes' )
        if ( associated(quantity_templates(i)%surfIndex) ) then
          call dump ( quantity_templates(i)%surfIndex, '      SurfIndex = ' )
        end if
        if ( associated(quantity_templates(i)%chanIndex) ) then
          call dump ( quantity_templates(i)%chanIndex, '      ChanIndex = ' )
        end if
      end if
    end do
  end subroutine DUMP_QUANTITY_TEMPLATES

  ! ---------------------------------- InflateQuantityTemplateDatabase --
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

  ! ---------------------------------------- myDisplayString -----
  subroutine myDisplayString ( index, advance )
    ! Given a quantity template, nullify all the pointers associated with it
    integer, intent(in) :: index
    character(len=*), intent(in), optional :: advance

    ! Executable code
    if ( index < 1 ) then
      call output ( '(string index < 1)' )
    elseif ( index /= index ) then
      call output ( '(string index is NaN)' )
    else
      call display_string ( index, advance )
    endif
  end subroutine myDisplayString

  ! ----------------------------------------NullifyQuantityTemplate -----
  subroutine NullifyQuantityTemplate ( Q )
    ! Given a quantity template, nullify all the pointers associated with it
    type ( QuantityTemplate_T ), intent(out) :: Q

    ! Executable code
    nullify ( q%surfs )
    nullify ( q%phi )
    nullify ( q%geodLat )
    nullify ( q%lon )
    nullify ( q%time )
    nullify ( q%solarTime )
    nullify ( q%solarZenith )
    nullify ( q%losAngle )
    nullify ( q%frequencies )
    nullify ( q%surfIndex )
    nullify ( q%chanIndex )
  end subroutine NullifyQuantityTemplate

  ! -----------------------------------  SetupNewQuantityTemplate  -----
  subroutine SetupNewQuantityTemplate ( qty, noInstances, noSurfs, &
    & noChans, coherent, stacked, regular, instanceLen, minorFrame, majorFrame, &
    & sharedVGrid, sharedHGrid, sharedFGrid )

  ! Set up a new quantity template according to the user input.  This may
  ! be based on a previously supplied template (with possible
  ! modifications), or created from scratch.

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

    ! Local variables
    integer :: noSurfsToAllocate        ! For allocations
    integer :: noInstancesToAllocate    ! For allocations

    ! Executable code
    qty%name = 0
    qty%quantityType = 0
    qty%noInstances = 1
    qty%noSurfs = 1
    qty%noChans = 1
    qty%coherent = .true.
    qty%stacked = .true.
    qty%regular = .true.
    qty%minorFrame = .false.
    qty%majorFrame = .false.
    qty%logBasis = .false.
    qty%minValue = - huge ( 0.0_r8 )
    qty%noInstancesLowerOverlap = 0
    qty%noInstancesUpperOverlap = 0
    qty%badValue = huge ( 0.0_r8 )
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
    qty%lo = 0.0_r8
    qty%signal = 0
    qty%sideband = 0
    qty%instrumentModule = 0
    qty%radiometer = 0
    qty%reflector = 0
    qty%molecule = 0

    ! Now, see if the user asked for modifications to this
    if ( present ( noInstances) )  qty%noInstances = noInstances
    if ( present ( noSurfs) )      qty%noSurfs = noSurfs
    if ( present ( noChans) )      qty%noChans = noChans
    if ( present ( regular) )      qty%regular = regular
    if ( present ( minorFrame) )   qty%minorFrame = minorFrame
    if ( present ( majorFrame) )   qty%majorFrame = majorFrame
    if ( present ( sharedVGrid ) ) qty%sharedVGrid = sharedVGrid
    if ( present ( sharedHGrid ) ) qty%sharedHGrid = sharedHGrid
    if ( present ( sharedFGrid ) ) qty%sharedFGrid = sharedFGrid

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
      qty%instanceLen = qty%noSurfs*qty%noChans
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
        & "qty%surfs", ModuleName )
    end if

    ! Now the horizontal coordinates

    if ( qty%sharedHGrid ) then
      nullify ( qty%phi, qty%geodLat, qty%lon, qty%time, qty%solarTime, &
        & qty%solarZenith, qty%losAngle )
    else
      call allocate_test ( qty%phi, noSurfsToAllocate, qty%noInstances, &
        & "qty%phi", ModuleName )
      call allocate_test ( qty%geodLat, noSurfsToAllocate, qty%noInstances, &
        & "qty%geodLat", ModuleName )
      call allocate_test ( qty%lon, noSurfsToAllocate, qty%noInstances, &
        & "qty%lon", ModuleName )
      call allocate_test ( qty%time, noSurfsToAllocate, qty%noInstances, &
        & "qty%time", ModuleName )
      call allocate_test ( qty%solarTime, noSurfsToAllocate, qty%noInstances, &
        & "qty%solarTime", ModuleName )
      call allocate_test ( qty%solarZenith, noSurfsToAllocate, qty%noInstances, &
        & "qty%solarZenith", ModuleName )
      call allocate_test ( qty%losAngle, noSurfsToAllocate, qty%noInstances, &
        & "qty%losAngle", ModuleName )
    end if

    if ( .not. qty%regular) then        !
      call allocate_test ( qty%surfIndex, qty%instanceLen, qty%noInstances, &
        & "qty%surfIndex", ModuleName )
      call allocate_test ( qty%chanIndex, qty%instanceLen, qty%noInstances, &
        & "qty%chanIndex", ModuleName )
    else
      nullify ( qty%surfIndex, qty%chanIndex )
    end if

  end subroutine SetupNewQuantityTemplate

!=============================================================================
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module QuantityTemplates
!=============================================================================

!
! $Log$
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
