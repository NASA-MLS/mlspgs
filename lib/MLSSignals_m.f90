! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MLSSignals_M

  ! Process the MLSSignals section of the L2 configuration file and deal with
  ! parsing signal request strings.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use DUMP_0, only: DUMP
  use Expr_M, only: Expr
  use Init_MLSSignals_m ! Everything
  use Intrinsic, only: Field_First, Field_indices, &
    & PHYQ_Dimensionless, PHYQ_Frequency, PHYQ_Indices, S_Time, l_emls
  use Lexer_Core, only: Print_Source
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error
  use MoreTree, only: Get_Boolean
  use Output_M, only: Output
  use String_Table, only: Display_String, Get_String
  use Toggles, only: Gen, Levels, Switches, Toggle
  use Trace_M, only: Trace_begin, Trace_end
  use Tree, only: Decorate, Decoration, Node_ID, Nsons, Source_Ref, Sub_Rosa, &
    & Subtree
  use Tree_Types, only: N_named, N_Set_One

  implicit none

  private ! So as not to re-export everything accessed by USE association.

  ! Public procedures and interfaces:
  public :: AddBandToDatabase, AddModuleToDatabase, AddRadiometerToDatabase
  public :: AddSignalToDatabase, AddSpectrometerTypeToDatabase, AreSignalsSuperset
  public :: DestroyBandDatabase, DestroyModuleDatabase
  public :: DestroyRadiometerDatabase, DestroySignal, DestroySignalDatabase
  public :: DestroySpectrometerType, DestroySpectrometerTypeDatabase, Dump
  public :: Dump_Bands, Dump_Radiometers, Dump_Signals, Dump_Spectrometertypes
  public :: GetAllModules, GetBandName, GetModuleFromRadiometer
  public :: GetModuleFromSignal, GetModuleName, GetNameOfSignal
  public :: GetRadiometerFromSignal, GetRadiometerName, GetSignal, GetSignalName
  public :: GetSpectrometerTypeName, IsModuleSpacecraft, MatchSignal, MLSSignals

  integer, public, parameter :: MaxSigLen = 80 ! Maximum length of a signal name

  ! =====     Defined Operators and Generic Identifiers     ==============
  
  interface Dump
    module procedure Dump_Bands, Dump_Radiometers, Dump_Signals, &
      & Dump_SpectrometerTypes
  end interface

  ! This boring type defines a module
  type, public :: Module_T
    integer :: Name                     ! Sub_rosa index of declaration's label
    integer :: Node                     ! Node of tree where module declared
    logical :: spaceCraft               ! Set if `module' is in fact s/c
  end type Module_T

  ! This type defines a radiometer.

  type, public :: Radiometer_T
    real(r8) :: LO                      ! Local oscillator in MHz
    integer :: InstrumentModule         ! Index in Modules database
    integer :: Prefix                   ! Sub_rosa index of declaration's label
    integer :: Suffix                   ! Sub_rosa index
  end type Radiometer_T

  ! The second type describes a band within that radiometer

  type, public :: Band_T
    real(r8) :: CenterFrequency         ! Negative if not present (wide filter)
    integer :: Prefix                   ! Sub_rosa index of declaration's label
    integer :: Radiometer               ! Index in Radiometers database
    integer :: SpectrometerType         ! Index in SpectrometerTypes database
    integer :: Suffix                   ! Sub_rosa index
  end type Band_T

  ! This type gives the information for specific spectrometer families.  For
  ! all apart from the WF4 spectrometers, we list frequencies and widths. 
  ! Otherwise, the arrays are empty.

  type, public :: SpectrometerType_T
    real(r8), pointer, dimension(:) :: Frequencies => NULL(), Widths => NULL()
    integer :: Name                     ! Sub_rosa index of declaration's label
    logical :: Deferred=.false.         ! "Frequencies/widths are deferred"
  end type SpectrometerType_T

  ! This is the key type; it describes a complete signal (one band, or a
  ! subset of the channels in one band).

  type, public :: Signal_T
    real(r8) :: CenterFrequency         ! Band local oscillator
    real(r8) :: LO                      ! Radiometer local oscillator
    real(r8), pointer, dimension(:) :: Frequencies=>NULL() ! Mainly a shallow copy
    real(r8), pointer, dimension(:) :: Widths=>NULL()  ! Mainly a shallow copy
    logical, pointer, dimension(:) :: Channels=>NULL() ! The ones actually used
    !                                   (This is for derived signals.  Every
    !                                   element is true for "real" signals).

    integer :: Band                     ! Index in Bands database
    logical :: Deferred = .false.       ! "Frequencies/widths are deferred"
    integer :: Index                    ! Index into master signals database
    integer :: InstrumentModule         ! Index in Modules database
    integer :: Name                     ! Sub_rosa index of declaration's label
    integer :: Radiometer               ! Index in Radiometers database
    integer :: SideBand                 ! -1=lower, +1=upper, 0=folded
    integer :: Spectrometer             ! Just a spectrometer number
    integer :: SpectrometerType         ! Index in SpectrometerTypes database
    integer :: Switch                   ! Just a switch number
  end type Signal_T

  ! Now some databases, the first are fairly obvious.
  !??? Should these be public ???

  type(module_T), public, save, pointer, dimension(:) :: Modules => NULL()
  type(band_T), public, save, pointer, dimension(:) :: Bands => NULL()
  type(radiometer_T), public, save, pointer, dimension(:) :: Radiometers => NULL()
  type(spectrometerType_T), public, save, pointer, dimension(:) ::&
    & SpectrometerTypes => NULL()

  ! This array is the signals database.  The first entries are the official
  ! `valid' signals in the instrument.  Later one can derive things from that.
  ! for subsets of channels etc.
  type(signal_T), public, save, pointer, dimension(:) :: Signals => NULL()
  integer, public, save :: Instrument = l_emls

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter, private :: ModuleName = &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains

  ! -------------------------------------------------  MLSSignals  -----
  subroutine MLSSignals ( ROOT )
    ! Process the MLSSignals section of the L2 configuration file.
    integer, intent(in) :: ROOT         ! The "cf" vertex for the section

    type(band_T) :: Band                ! To be added to the database
    integer :: Channels                 ! subtree index of field
    logical :: Deferred                 ! "Frequencies/widths are deferred"
    integer :: Error                    ! Error level seen so far
    integer :: Field                    ! Field index -- f_something
    integer :: First                    ! "first" field of "spectrometer"
    logical :: Got(field_first:last_Signal_Field)   ! "Got this field already"
    integer :: Gson                     ! Son of Son.
    integer :: I, J, K                  ! Subscript and loop inductor.
    integer :: Key                      ! Indexes the spec_args vertex.
    integer :: Last                     ! "last" field of "spectrometer"
    integer :: Name                     ! sub_rosa of label of specification,
                                        ! if any, else zero.
    type(radiometer_T) :: Radiometer    ! To be added to the database
    integer :: Son                      ! Some subtree of root.
    type(signal_T) :: Signal            ! To be added to the database
    type(spectrometerType_T) :: SpectrometerType ! To be added to the database
    real(r8) :: Start                   ! "start" field of "spectrometer"
    real(r8) :: Step                    ! "step" field of "spectrometer"
    type(module_T) :: thisModule        ! To be added to database
    logical :: TIMING                   ! For S_Time
    real :: T1, T2                      ! For S_Time
    integer :: Units(2)                 ! of an expression
    double precision :: Value(2)        ! of an expression
    real(r8) :: Width                   ! "width" field of "spectrometer"

    ! Error message codes
    integer, parameter :: AllOrNone = 1 ! All of a set of fields, or none
    integer, parameter :: AtLeastOne = allOrNone + 1   ! At least one
    integer, parameter :: DeferredChannels = atLeastone + 1 ! Deferred .neqv channels
    integer, parameter :: BadMix = deferredChannels + 1 ! Disallowed mixture of
    !                                     fields.
    integer, parameter :: WrongUnits = badMix + 1 ! Field has wrong units

    error = 0
    timing = .false.
    if ( toggle(gen) ) call trace_begin ( "MLSSignals", root )
    do i = 2, nsons(root)-1 ! skip "MLSSignals" at "begin" and "end"
      son = subtree(i,root) ! A spec_args vertex now
      if ( node_id(son) == n_named ) then
        name = sub_rosa(subtree(1, son))
        key = subtree(2, son)
      else
        name = 0
        key = son
      end if
      ! node_id(key) is now n_spec_args

      got = .false.
      select case ( decoration(subtree(1,decoration(subtree(1,key)))) )


      case ( s_band ) ! ...................................  BAND  .....
        band%prefix = name
        band%centerFrequency = -1.0_r8 ! "The 'frequency' field is absent"
        do j = 2, nsons(key)
          son = subtree(j,key)
          field = decoration(subtree(1,son))
          gson = subtree(2,son)
          got(field) = .true.
          select case ( field )
          case ( f_centerFrequency )
            call expr_check ( gson, units, value, field, phyq_frequency )
            band%centerFrequency = value(1)
          case ( f_suffix )
            band%suffix = sub_rosa(gson)
          case ( f_radiometer )
            band%radiometer = decoration(decoration(gson))
          case ( f_spectrometerType )
            band%spectrometerType = decoration(decoration(gson))
          case default
            ! Shouldn't get here if the type checker worked
          end select
        end do ! i = 2, nsons(key)
        call decorate ( key, addBandToDatabase ( bands, band ) )


      case ( s_module ) ! ............................  MODULE  ........
        thisModule%name = name
        thisModule%spaceCraft = .false.
        thisModule%node = decoration(name)
        do j = 2,nsons(key)
          son = subtree(j,key)
          field = decoration(subtree(1,son))
          if ( nsons(son) > 1 ) then
            gson = subtree(2,son)
          else
            gson = son
          end if
          got(field) = .true.
          select case ( field )
          case (f_spaceCraft)
            thisModule%spacecraft = get_boolean(son)
          case default
            ! Shouldn't get here if parser worked
          end select
        end do
        call decorate ( key, AddModuleToDatabase (modules, thisModule ) )


      case ( s_radiometer ) ! .......................  RADIOMETER  .....
        radiometer%prefix = name
        do j = 2, nsons(key)
          son = subtree(j,key)
          field = decoration(subtree(1,son))
          gson = subtree(2,son)
          select case ( field )
          case ( f_lo )
            call expr_check ( gson, units, value, field, phyq_frequency )
            radiometer%lo = value(1)
          case ( f_suffix )
            radiometer%suffix = sub_rosa(gson)
          case ( f_module )
            radiometer%instrumentModule = decoration(decoration(gson))
          case default
            ! Shouldn't get here if the type checker worked
          end select
        end do ! i = 2, nsons(key)
        call decorate ( key, addRadiometerToDatabase ( radiometers, radiometer ) )


      case ( s_signal ) ! ..........................  VALIDSIGNAL  .....
        signal%sideband = 0
        do j = 2, nsons(key)
          son = subtree(j,key)
          field = decoration(subtree(1,son))
          gson = subtree(2,son)
          got(field) = .true.
          select case ( field )
          case ( f_band )
            signal%band = decoration(decoration(gson))
          case ( f_channels )
            channels = son
          case ( f_spectrometer )
            call expr_check ( gson, units, value, field, phyq_dimensionless )
            signal%spectrometer = value(1)
          case ( f_switch )
            call expr_check ( gson, units, value, field, phyq_dimensionless )
            signal%switch = value(1)
          case default
            ! Shouldn't get here if the type checker worked
          end select
        end do ! i = 2, nsons(key)
        signal%name = name
        ! Set default values for remaining parameters
        signal%radiometer = bands(signal%band)%radiometer
        signal%lo = radiometers(signal%radiometer)%lo
        signal%instrumentModule = radiometers(signal%radiometer)%instrumentModule
        signal%spectrometerType = bands(signal%band)%spectrometerType
        signal%centerFrequency = bands(signal%band)%centerFrequency
        signal%deferred = spectrometerTypes(signal%spectrometerType)%deferred
        if ( signal%deferred .neqv. got(f_channels) ) &
          & call announceError ( deferredChannels )
        ! For the wide filters, we specify frequency etc. here.
        if ( got(f_channels) ) then
          if ( error == 0 ) then
            call allocate_Test ( signal%frequencies, nsons(channels)-1, &
              & 'signal%frequencies', moduleName )
            call allocate_Test ( signal%widths, nsons(channels)-1, &
              & 'signal%widths', moduleName)
            call allocate_Test ( signal%channels, nsons(channels)-1, &
              & 'signal%channels', moduleName)
            do k = 2, nsons(channels)
              call expr ( subtree(k,channels), units, value )
              signal%frequencies(k-1) = value(1)
              signal%widths(k-1) = value(2)
            end do
            if ( any(units /= phyq_frequency) ) &
              & call announceError ( wrongUnits, f_channels, &
                & (/ phyq_frequency /) )
            signal%channels = .true.
          end if
        else
          signal%frequencies => spectrometerTypes(signal%spectrometerType)% &
            & frequencies
          signal%widths => spectrometerTypes(signal%spectrometerType)%widths
        end if
        call decorate ( key, addSignalToDatabase ( signals, signal ) )
        signals(size(signals))%index = size(signals)
        ! Now nullify pointers so they don't get hosed later by allocate_test
        nullify ( signal%channels )
        nullify ( signal%frequencies )
        nullify ( signal%widths )


      case ( s_spectrometerType ) ! ...........  SPECTROMETERTYPE  .....
        spectrometerType%name = name
        deferred = .false.
        first = 0
        do j = 2, nsons(key)
          son = subtree(j,key)
          field = decoration(subtree(1,son))
          if ( nsons(son) > 1 ) then
            gson = subtree(2,son)
          else
            gson = son                  ! For case of /deferred
          end if
          got(field) = .true.
          select case ( field )
          case ( f_channels )
            channels = son
          case ( f_deferred )
            deferred = get_boolean(son)
          case ( f_first, f_last )
            call expr_check ( gson, units, value, field, phyq_dimensionless )
            select case ( field )
            case ( f_first )
              first = value(1)
            case ( f_last )
              last = value(1)
            end select
          case ( f_start, f_step, f_width )
            call expr_check ( gson, units, value, field, phyq_frequency )
            select case ( field )
            case ( f_start )
              start = value(1)
            case ( f_step )
              step = value(1)
            case ( f_width )
              width = value(1)
            end select
          case default
            ! Shouldn't get here if the type checker worked
          end select
        end do ! i = 2, nsons(key)
        if ( got(f_channels) ) then
          if ( any(got( (/ f_last, f_start, f_step, f_width /) )) ) &
            & call announceError ( badMix, f_channels, &
            & (/ f_last, f_start, f_step, f_width /) )
          if ( error == 0 ) then
            call allocate_Test ( spectrometerType%frequencies, &
              & nsons(channels)-first, 'spectrometerType%frequencies', &
              & moduleName, lowBound = first )
            call allocate_Test ( spectrometerType%widths, &
              & nsons(channels)-first, 'spectrometerType%widths', &
              & moduleName, lowBound = first )
            do k = 2, nsons(channels)
              call expr ( subtree(k,channels), units, value )
              spectrometerType%frequencies(k-2+first) = value(1)
              spectrometerType%widths(k-2+first) = value(2)
            end do
            if ( any(units /= phyq_frequency) ) &
              & call announceError ( wrongUnits, f_channels, &
                & (/ phyq_frequency /) )
          end if
        end if
        if ( got(f_start) ) then
          if ( .not. all( got((/ f_last, f_step, &
          & f_width /)) ) ) call announceError ( allOrNone, f_start, &
          & (/ f_last, f_step, f_width /) )
          if ( error == 0 ) then
            k = last - first 
            call allocate_Test ( spectrometerType%frequencies, k-first, &
              & 'spectrometerType%frequencies', moduleName, lowBound = first )
            call allocate_Test ( spectrometerType%widths, k-first, &
              & 'spectrometerType%widths', moduleName, lowBound = first )
            spectrometerType%widths = width
            spectrometerType%frequencies(first) = start
            do k = first+1, last
              spectrometerType%frequencies(k) = start + (k-first) * step
            end do ! k
          end if
        end if
        if ( deferred ) then            ! For deferred types, wait till later
          if ( any(got( (/ f_channels, f_last, f_start, f_step, f_width /) )) ) &
            & call announceError ( badMix, f_channels, &
            & (/ f_last, f_start, f_step, f_width /) )
          if ( error == 0 ) then
            spectrometerType%deferred = deferred
            nullify(spectrometerType%frequencies)
            nullify(spectrometerType%widths)
          end if
        end if
        if ( .not. any(got( (/ f_channels, f_deferred, f_start /) )) ) &
          & call announceError ( atLeastOne, f_channels, (/ f_deferred, f_start /) )
        if ( error == 0 ) call decorate ( key, addSpectrometerTypeToDatabase ( &
          & spectrometerTypes, spectrometerType ) )

        ! Nullify pointers to temporary stuff so it doesn't get hosed later
        nullify ( spectrometerType%frequencies )
        nullify ( spectrometerType%widths )

      case ( s_time ) ! ...................................  TIME  .....
        if ( timing ) then
          call sayTime
        else
          call cpu_time ( t1 )
          timing = .true.
        end if
      case default
        ! Shouldn't get here if the type checker worked
      end select

    end do

    if ( error > 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "Unable to create MLSSignals database" )

    if ( toggle(gen) ) then
      if ( levels(gen) > 0 .or. index(switches, 'S') /= 0 ) then
        call dump ( radiometers )
        call dump ( spectrometerTypes )
        call dump ( bands )
        call dump ( signals )
      end if
      call trace_end ( "MLSSignals" )
    end if
    if ( timing ) call sayTime

    contains
    ! --------------------------------------------  AnnounceError  -----
    subroutine AnnounceError ( Code, FieldIndex, MoreFields )
      integer, intent(in) :: Code       ! Index of error message
      integer, intent(in), optional :: FieldIndex ! f_...
      integer, intent(in), optional :: MoreFields(:)

      integer :: I
      integer :: Source

      error = max(error,1)
      source = source_ref ( son )
      call output ( 'At ' )
      call print_source ( source )
      call output ( ' MLSSignals complained: ' )
      select case ( code )
      case ( allOrNone )
        call output ( 'Either all of the fields ' )
        call display_string ( field_indices(fieldIndex) )
        do i = 1, size(moreFields)
          if ( i == size(moreFields) ) then
            call output ( ' and ' )
          else
            call output ( ', ' )
          end if
          call display_string ( field_indices(moreFields(i)) )
        end do ! i
        call output ( ' shall appear, or none of them shall.', advance = 'yes' )
      case ( atLeastOne )
        call output ( 'At least one of the fields ' )
        call display_string ( field_indices(fieldIndex) )
        do i = 1, size(moreFields)
          if ( i == size(moreFields) ) then
            call output ( ' or ' )
          else
            call output ( ', ' )
          end if
          call display_string ( field_indices(moreFields(i)) )
        end do ! i
        call output ( ' shall appear.', advance='yes' )
      case ( badMix )
        call output ( 'If the ' )
        call display_string ( field_indices(fieldIndex) )
        call output ( ' field appears, no other field except ' )
        do i = 1, size(moreFields)
          if ( i == size(moreFields) ) then
            if ( i > 1 ) call output ( ' or ' )
          else if ( i > 1 ) then
            call output ( ', ' )
          end if
          call display_string ( field_indices(moreFields(i)) )
        end do ! i
        call output ( ' shall appear.', advance='yes' )
      case ( deferredChannels )
        call output ( "Channels shall be specified if and only if the band's" )
        call output ( ' radiometer has deferred channels', advance='yes' )
      case ( wrongUnits )
        call output ( 'The values of the ' )
        call display_string ( field_indices(fieldIndex) )
        call output ( ' field have the wrong units -- ' )
        call display_string ( phyq_indices(moreFields(1)) )
        call output ( ' required.', advance='yes' )
      end select
    end subroutine AnnounceError

    ! -----------------------------------------------  Expr_Check  -----
    subroutine Expr_Check ( Root, Units, Value, Field, NeededUnits )
    ! Evaluate an expression and check its units.
      integer, intent(in) :: Root, Field, NeededUnits
      integer, intent(out) :: Units(2)
      real(r8), intent(out) :: Value(2)
      call expr ( root, units, value )
      if ( units(1) /= neededUnits ) &
        & call announceError ( wrongUnits, field, (/ neededUnits /) )
    end subroutine Expr_Check

    ! --------------------------------------------------  SayTime  -----
    subroutine SayTime
      call cpu_time ( t2 )
      call output ( "Timing for MLSSignals = " )
      call output ( dble(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime

  end subroutine MLSSignals

  ! ------------------------------------------  AddBandToDatabase  -----
  integer function AddBandToDatabase ( Database, Item )
    type(band_T), dimension(:), pointer :: Database
    type(band_T), intent(in) :: Item

    ! Local variables
    type (Band_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddBandToDatabase = newSize
  end function AddBandToDatabase

  ! ----------------------------------------  AddModuleToDatabase  -----
  integer function AddModuleToDatabase ( Database, Item )
    type(module_T), dimension(:), pointer :: Database
    type(module_T), intent(in) :: Item

    ! Local variables
    type (Module_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddModuleToDatabase = newSize
  end function AddModuleToDatabase

  ! ------------------------------------  AddRadiometerToDatabase  -----
  integer function AddRadiometerToDatabase ( Database, Item )
    type(radiometer_T), dimension(:), pointer :: Database
    type(radiometer_T), intent(in) :: Item

    ! Local variables
    type (Radiometer_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddRadiometerToDatabase = newSize
  end function AddRadiometerToDatabase

  ! ----------------------------------------  AddSignalToDatabase  -----
  integer function AddSignalToDatabase ( Database, Item )
    type(signal_T), dimension(:), pointer :: Database
    type(signal_T), intent(in) :: Item

    ! Local variables
    type (Signal_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddSignalToDatabase = newSize
  end function AddSignalToDatabase

  ! ------------------------------  AddSpectrometerTypeToDatabase  -----
  integer function AddSpectrometerTypeToDatabase ( Database, Item )
    type(spectrometerType_T), dimension(:), pointer :: Database
    type(spectrometerType_T), intent(in) :: Item

    ! Local variables
    type (spectrometerType_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddSpectrometerTypeToDatabase = newSize
  end function AddSpectrometerTypeToDatabase

  ! ----------------------------------------  DestroyBandDatabase  -----
  subroutine DestroyBandDatabase ( Bands )
    type(band_T), dimension(:), pointer :: Bands
    integer :: Status
    if ( associated(bands) ) then
      deallocate ( bands, stat = status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_DeAllocate // 'Band database' )
    end if
  end subroutine DestroyBandDatabase

  ! --------------------------------------  DestroyModuleDatabase  -----
  subroutine DestroyModuleDatabase ( Modules )
    type(module_T), dimension(:), pointer :: Modules
    integer :: Status
    if ( associated(modules) ) then
      deallocate ( modules, stat = status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_DeAllocate // 'Module database' )
    end if
  end subroutine DestroyModuleDatabase

  ! ----------------------------------  DestroyRadiometerDatabase  -----
  subroutine DestroyRadiometerDatabase ( Radiometers )
    type(radiometer_T), dimension(:), pointer :: Radiometers
    integer :: Status
    if ( associated(radiometers) ) then
      deallocate ( radiometers, stat = status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_DeAllocate // 'Radiometer database' )
    end if
  end subroutine DestroyRadiometerDatabase

  ! ----------------------------------------------  DestroySignal  -----
  subroutine DestroySignal ( Signal )
    ! Destroy one signal in the signals database (or elsewhere)
    type(signal_T), intent(inout) :: Signal

    call deallocate_test ( signal%channels, 'Signal%channels', moduleName )
    ! Don't destroy Frequencies or Widths unless signal%Deferred.  Those
    ! fields are shallow copies here.  They're destroyed in
    ! DestroySpectrometerType.
    if ( signal%deferred ) then
      call deallocate_test ( signal%frequencies, "signal%frequencies", &
        & moduleName )
      call deallocate_test ( signal%widths, "signal%widths", moduleName )
    end if
  end subroutine DestroySignal

  ! --------------------------------------  DestroySignalDatabase  -----
  subroutine DestroySignalDatabase ( Signals )
    type(signal_T), dimension(:), pointer :: Signals
    integer :: I, Status

    ! Executable code
    if ( associated(signals) ) then
      do i = 1, size(signals)
        call destroySignal ( signals(i) )
      end do
      deallocate ( signals, stat = status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_DeAllocate // 'Signal database' )
    end if
  end subroutine DestroySignalDatabase

  ! ------------------------------------  DestroySpectrometerType  -----
  subroutine DestroySpectrometerType ( SpectrometerType )
    ! Destroy one SpectrometerType object.
    type(spectrometerType_T), intent(inout) :: SpectrometerType
    call deallocate_Test ( spectrometerType%frequencies, &
      & 'Spectrometer%frequencies', moduleName )
    call deallocate_Test ( spectrometerType%widths, 'Spectrometer%widths', &
      &  moduleName )
  end subroutine DestroySpectrometerType

  ! ----------------------------  DestroySpectrometerTypeDatabase  -----
  subroutine DestroySpectrometerTypeDatabase ( SpectrometerTypes )
    type(spectrometerType_T), dimension(:), pointer :: spectrometerTypes
    integer :: I, Status
    if ( associated(spectrometerTypes) ) then
      do i = 1, size(spectrometerTypes)
        call destroySpectrometerType ( spectrometerTypes(i) )
      end do
      deallocate ( spectrometerTypes, stat = status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_DeAllocate // 'Spectrometer database' )
    end if
  end subroutine DestroySpectrometerTypeDatabase

  ! --------------------------------------------------  DumpBands  -----
  subroutine DUMP_BANDS ( BANDS )
    type (Band_T), intent(in) :: BANDS(:)
    integer :: i
    call output ( 'BANDS: SIZE = ')
    call output ( size(bands), advance='yes' )
    do i = 1, size(bands)
      call output ( i, 2 )
      call output ( ': ' )
      call display_string (bands(i)%prefix)
      call output ( ':' )
      call display_string (bands(i)%suffix, advance='yes', strip=.true. )
      call output ( '   Radiometer: ')
      call output ( bands(i)%radiometer )
      call output ( ' - ' )
      call display_string ( radiometers(bands(i)%radiometer)%prefix )
      call output ( '   SpectrometerType: ' )
      call output ( bands(i)%spectrometerType )
      call output ( ' - ' )
      call display_string ( spectrometerTypes(bands(i)%spectrometerType)%name, &
        & advance='yes' )
      call output ( '   Frequency: ')
      call output ( bands(i)%centerFrequency )
    end do
  end subroutine DUMP_BANDS

  ! -------------------------------------------  Dump_Radiometers  -----
  subroutine DUMP_RADIOMETERS ( RADIOMETERS )
    type (Radiometer_T), intent(in) :: RADIOMETERS(:)
    integer :: i
    call output ( 'RADIOMETERS: SIZE = ')
    call output ( size(radiometers), advance='yes' )
    do i = 1, size(radiometers)
      call output ( i,1 )
      call output ( ': ')
      call display_string (radiometers(i)%prefix)
      call output ( ':' )
      call display_string (radiometers(i)%suffix, advance='yes', strip=.true. )
      call output ( '   Module: ')
      call output ( radiometers(i)%instrumentModule )
      call output ( ' - ' )
      call display_string ( modules(radiometers(i)%instrumentModule)%name, advance='yes' ) 
      call output ( '   LO: ')
      call output ( radiometers(i)%lo,advance='yes' )
    end do
  end subroutine DUMP_RADIOMETERS

  ! ------------------------------------------------  DumpSignals  -----
  subroutine DUMP_SIGNALS ( SIGNALS, DETAILS )
    type (signal_T), intent(in) :: SIGNALS(:)
    logical, intent(in), optional :: Details ! false => don't dump frequencies
    integer :: I
    logical :: My_Details
    character (len=80) :: Str
    my_details = .true.
    if ( present(details) ) my_details = details
    call output ( 'SIGNALS: SIZE = ')
    call output ( size(signals), advance='yes' )
    do i = 1, size(signals)
      call output ( 'Signal ' )
      call output ( i )
      if ( signals(i)%name > 0 ) then
        call output ( ': ' )
        call display_string ( signals(i)%name, advance='yes' )
      else
        call output ( '', advance='yes' )
      end if
      call output ( '   Module: ')
      call output ( signals(i)%instrumentModule )
      call output ( ' - ' )
      if ( associated(modules) ) then
        call display_string ( modules(signals(i)%instrumentModule)%name, &
          & advance='yes' )
      else
        call output ( 'Cannot get module name', advance='yes' )
      end if
      call output ( '   Radiometer: ')
      call output ( signals(i)%radiometer )
      call output ( ' - ' )
      if ( associated(radiometers) ) then
        call getRadiometerName ( signals(i)%radiometer, str )
        call output ( TRIM(str) )
      else
        call output ( 'Cannot get radiometer name', advance='yes' )
      end if
      call output ( '   First LO: ')
      call output ( signals(i)%lo, advance='yes')
      call output ( '   Band: ')
      call output ( signals(i)%band )
      call output (' - ')
      if ( associated(bands) ) then
        call getBandName ( signals(i)%band, str )
        call output ( TRIM(str) )
      else
        call output ( 'Cannot get band name', advance='yes' )
      end if
      call output ( '   Band center frequency: ')
      call output ( signals(i)%centerFrequency, advance='yes')
      call output ( '   SpectrometerType: ')
      call output ( signals(i)%spectrometerType )
      call output ( ' - ' )
      if ( associated(spectrometerTypes) ) then
        call display_string ( spectrometerTypes(signals(i)%spectrometerType)%name )
      else
        call output ( 'Cannot get spectrometer type name', advance='yes' )
      end if
      call output ( '   Channels: ' )
      call output ( lbound(signals(i)%frequencies,1), 3 )
      call output ( ':' )
      call output ( ubound(signals(i)%frequencies,1), 3, advance='yes' )
      call output ( '   Sideband: ' )
      call output ( signals(i)%sideband, advance='yes')
      if ( my_details ) then
        call output ( '   Frequencies' )
        if ( signals(i)%deferred ) call output ( '(deferred)' )
        call output ( ':', advance='yes' )
        call dump ( signals(i)%frequencies )
        call output ( '   Widths:' )
        if ( signals(i)%deferred ) call output ( '(deferred)' )
        call output ( ':', advance='yes' )
        call dump ( signals(i)%widths )
      else
        call output ( '   Frequencies and widths are' )
        if ( .not. signals(i)%deferred ) call output ( ' not' )
        call output ( ' deferred', advance='yes' )
      end if ! my_details
      if (associated(signals(i)%channels)) then
        call output ( '   Channel Flags:', advance='yes' )
        call dump( signals(i)%channels )
      else
        call output ( '   All channels selected', advance='yes' )
      end if
    end do
  end subroutine DUMP_SIGNALS

  ! --------------------------------------  DumpSpectrometerTypes  -----
  subroutine DUMP_SPECTROMETERTYPES ( SPECTROMETERTYPES )
    type (SpectrometerType_T), intent(in) :: SPECTROMETERTYPES(:)
    integer :: i
    call output ( 'SPECTROMETERTYPES: SIZE = ')
    call output ( size(spectrometerTypes), advance='yes' )
    do i = 1, size(spectrometerTypes)
      call output ( i,1 )
      call output ( ': ')
      call display_string ( spectrometerTypes(i)%name, advance='yes' )
      if ( associated(spectrometerTypes(i)%frequencies) ) then
        call output ( '  Channels: ' )
        call output ( lbound(spectrometerTypes(i)%frequencies,1), 3 )
        call output ( ':' )
        call output ( ubound(spectrometerTypes(i)%frequencies,1), 3, &
          & advance='yes' )
        call output ( '  Frequencies:', advance='yes' )
        call dump ( spectrometerTypes(i)%frequencies )
        call output ( '  Widths:', advance='yes' )
        call dump ( spectrometerTypes(i)%widths )
      else
        call output ('   Frequencies and widths deferred.', advance='yes' )
      end if
    end do
  end subroutine DUMP_SPECTROMETERTYPES

  ! ----------------------------------------------  GetAllModules  -----
  subroutine GetAllModules(moduleNodes)
    ! Return tree nodes for all modules
    integer, dimension(:), pointer :: moduleNodes
    
    call Allocate_Test ( moduleNodes, size(modules), 'ModuleNodes', ModuleName )
    moduleNodes = modules%node
  end subroutine GetAllModules

  ! ------------------------------------------------  GetBandName  -----
  subroutine GetBandName(band, string_text, sideband, noSuffix)
    ! Given an index in the Bands database, place band name in string
    integer, intent(in) :: BAND                   ! Database index
    character(len=*), intent(out) :: STRING_TEXT  ! Result
    integer, intent(in), optional :: SIDEBAND     ! L_Folded, L_Lower, L_Upper
    logical, intent(in), optional :: NOSUFFIX     ! Omit suffix if present and true

    ! Local variables
    logical :: MY_NOSUFFIX
    integer :: MY_SIDEBAND
    character (len=1) :: SB_CHAR        ! U/L for sideband
    character (len=1) :: SB_CHARS(-1:1) = (/ 'L', ' ', 'U' /)

    ! Executable code
    my_noSuffix = .false.
    my_sideband = 0
    if ( present(noSuffix) ) my_noSuffix = noSuffix
    if ( present(sideband) ) my_sideband = sideband

    if ( my_sideband < -1 .or. my_sideband > 1 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, 'Illegal sideband' )
    sb_char = sb_chars(my_sideband)
    
    call get_string ( bands(band)%prefix, string_text, cap=.true. )
    string_text = string_text(1:LEN_TRIM(string_text)-1) // TRIM(sb_char) // &
      & string_text(LEN_TRIM(string_text):LEN_TRIM(string_text))
    if ( (.not. my_noSuffix) .and. &
      &  (len_trim(string_text) < len(string_text)) ) then
      string_text = TRIM(string_text) // ':'
      call get_string ( bands(band)%suffix,&
        & string_text(LEN_TRIM(string_text)+1:), cap=.true., strip=.true. )
    endif

  end subroutine GetBandName

  ! ------------------------------------  GetModuleFromRadiometer  -----
  integer function GetModuleFromRadiometer(radiometer)
    ! Returns module field from given radiometer given as database index
    integer, intent(in) :: radiometer
    GetModuleFromRadiometer = radiometers(radiometer)%instrumentModule
  end function GetModuleFromRadiometer

  ! ----------------------------------------  GetModuleFromSignal  -----
  integer function GetModuleFromSignal(signal)
    ! Returns module field from given signal given as database index
    integer, intent(in) :: signal
    GetModuleFromSignal = signals(signal)%instrumentModule
  end function GetModuleFromSignal

  ! ----------------------------------------------  GetModuleName  -----
  subroutine GetModuleName(instrumentModule, string_text)
    ! Given the index in the module database, returns module name in mixed case
    integer, intent(in) :: instrumentModule
    character (len=*), intent(out) :: string_text
    call Get_String ( modules(instrumentModule)%name, string_text )
  end subroutine GetModuleName

  ! --------------------------------------------  GetNameOfSignal  -----
  subroutine GetNameOfSignal ( Signal, String_text, NoRadiometer, NoBand, &
    & NoSwitch, NoSpectrometer, NoChannels, NoSuffix, sideband )
    ! Given a signal object, this routine constructs a full signal name.
    type(signal_T), intent(in) :: SIGNAL
    character(len=*), intent(inout) :: STRING_TEXT
    logical, intent(in), optional :: NORADIOMETER
    logical, intent(in), optional :: NOBAND
    logical, intent(in), optional :: NOSWITCH
    logical, intent(in), optional :: NOSPECTROMETER
    logical, intent(in), optional :: NOCHANNELS
    logical, intent(in), optional :: NOSUFFIX
    integer, intent(in), optional :: SIDEBAND

    ! Local variables
    logical :: First     ! First channel in signal text
    integer :: I, J, L
    logical :: MY_NORADIOMETER, MY_NOBAND, MY_NOSWITCH
    logical :: MY_NOSPECTROMETER, MY_NOCHANNELS
    integer :: MY_SIDEBAND
    character (len=8) :: word

    ! Executable code
    string_text       = ''
    my_noRadiometer   = .false.
    my_noBand         = .false.
    my_noSwitch       = .false.
    my_noSpectrometer = .false.
    my_noChannels     = .false.
    my_sideband       = signal%sideband

    if ( present(noRadiometer) )   my_noRadiometer =   noRadiometer
    if ( present(noBand) )         my_noBand =         noBand
    if ( present(noSwitch) )       my_noSwitch =       noSwitch
    if ( present(noSpectrometer) ) my_noSpectrometer = noSpectrometer
    if ( present(noChannels) )     my_noChannels =     noChannels
    if ( present(sideband) )       my_sideband =       sideband

    if ( .not. my_noRadiometer ) call GetRadiometerName ( signal%radiometer, &
      & string_text, noSuffix=noSuffix )

    if ( .not. my_noBand ) then
      if ( (len_trim(string_text) /= 0) .and. &
        &  (len_trim(string_text)<len(string_text)) ) &
        &  string_text = TRIM(string_text) // '.'
      call GetBandName ( signal%band, &
        & string_text(LEN_TRIM(string_text)+1:), sideband=my_sideband, &
        & noSuffix=noSuffix )
    end if

    if ( .not. my_noSwitch ) then
      if ( (len_trim(string_text) /= 0) .and. &
        &  (len_trim(string_text)+1<len(string_text)) ) &
        &  string_text = TRIM(string_text) // '.S'
      write (word,'(I8)') signal%switch
      word = adjustl(word)
      if ( len_trim(string_text)+len_trim(word) < len(string_text) )&
        & string_text = TRIM(string_text) // TRIM(word)
    end if

    if ( .not. my_noSpectrometer ) then
      if ( (len_trim(string_text) /= 0) .and. &
        &  (len_trim(string_text)<len(string_text)) ) &
        &  string_text = TRIM(string_text) // '.'
      call GetSpectrometerTypeName ( signal%spectrometerType, &
        & signal%spectrometer, string_text(LEN_TRIM(string_text)+1:) )
    end if

    if ( .not. my_noChannels .and. associated(signal%channels) ) then
      l = len_trim(string_text)
      call addToSignalString ( '.C' )
      i = lbound(signal%channels, 1)
      first = .true.
oc:   do
        do
          if ( i > ubound(signal%channels, 1) ) exit oc
          if ( signal%channels(i) ) exit
          i = i + 1
        end do
        if ( .not. first ) call addToSignalString ( '+' )
        first = .false.
        j = i
        do while ( j < ubound(signal%channels, 1) )
          if ( .not. signal%channels(j+1) ) exit
        end do
        if ( j > i ) then
          write ( word, * ) i
          call addToSignalString ( word )
          call addToSignalString ( ':' )
          write ( word, * ) j
          call addToSignalString ( word )
        else
          write ( word, * ) i
          call addToSignalString ( word )
        end if
        i = j + 1
      end do oc
    end if

  contains
    subroutine AddToSignalString ( Text )
      ! Assumes L = len_trim(string_text) on entry, preserves it
      character(len=*), intent(in) :: Text
      integer :: Lt
      lt = len_trim(adjustl(text))
      if ( l+lt <= len(string_text) ) &
        & string_text(l+1:l+lt) = trim(adjustl(text))
      l = l + lt
    end subroutine AddToSignalString
  end subroutine GetNameOfSignal

  ! ------------------------------------  GetRadiometerFromSignal  -----
  integer function GetRadiometerFromSignal(signal)
    ! Returns radiometer field from given signal given as database index
    integer, intent(in) :: signal
    GetRadiometerFromSignal = signals(signal)%radiometer
  end function GetRadiometerFromSignal

  ! ------------------------------------------  GetRadiometerName  -----
  subroutine GetRadiometerName(radiometer, string_text, noSuffix)
    ! Given an index in the Radiometers database, place radiometer name
    ! in string
    integer, intent(in) :: RADIOMETER             ! Database index
    character(len=*), intent(out) :: STRING_TEXT  ! Result
    logical, intent(in), optional :: NOSUFFIX     ! Omit suffix if present and true

    ! Local variables
    logical :: MY_NOSUFFIX

    ! Executable code
    my_noSuffix = .false.
    if ( present(noSuffix) ) my_noSuffix = noSuffix

    call get_string ( radiometers(radiometer)%prefix, string_text, cap=.true., &
      & strip=.true. )
    if ( (.not. my_noSuffix) .and. &
      &  (len_trim(string_text) < len(string_text)) ) then
      string_text = TRIM(string_text) // ':'
      call get_string ( radiometers(radiometer)%suffix, &
        & string_text(LEN_TRIM(string_text)+1:), cap=.true., strip=.true. )
    endif

  end subroutine GetRadiometerName

  ! --------------------------------------------------- GetSignal ------
  type (Signal_T) function GetSignal(signal)
    ! Given the database index, this routine returns the signal data structure
    integer, intent(in) :: SIGNAL       ! Requested signal
    
    GetSignal = signals(signal)
  end function GetSignal    

  ! ----------------------------------------------  GetSignalName  -----
  subroutine GetSignalName ( Signal, String_text, NoRadiometer, NoBand, &
    & NoSwitch, NoSpectrometer, NoChannels, NoSuffix, SIDEBAND )
    ! Given an index in the signals database, this routine constructs a
    ! full signal name.
    integer, intent(in) :: SIGNAL                 ! Database index
    character(len=*), intent(inout) :: STRING_TEXT
    logical, intent(in), optional :: NORADIOMETER
    logical, intent(in), optional :: NOBAND
    logical, intent(in), optional :: NOSWITCH
    logical, intent(in), optional :: NOSPECTROMETER
    logical, intent(in), optional :: NOCHANNELS
    logical, intent(in), optional :: NOSUFFIX
    integer, intent(in), optional :: SIDEBAND

    call getNameOfSignal ( signals(signal), string_text, noRadiometer, noBand, &
      & noSwitch, noSpectrometer, noChannels, noSuffix, sideband )
  end subroutine GetSignalName 

  ! ----------------------------------------  GetSpectrometerName  -----
  subroutine GetSpectrometerTypeName(spectrometerType, number, string_text)
    ! Place spectrometer name and number in string
    integer, intent(in) :: SPECTROMETERTYPE
    integer, intent(in) :: NUMBER
    character (len=*), intent(out) :: STRING_TEXT

    ! Local variables
    character (len=8) :: word

    ! Executable code
    call get_string ( spectrometerTypes(spectrometerType)%name, string_text )
    if ( len_trim(string_text) < len(string_text) ) &
      & string_text = TRIM(string_text) // '-'
    write ( word,'(I8)' ) number
    word = adjustl(word)
    if ( len_trim(string_text)+len_trim(word) < len(string_text) ) &
      & string_text = TRIM(string_text) // TRIM(word)
  end subroutine GetSpectrometerTypeName

  ! -----------------------------------------  IsModuleSpacecraft  -----
  logical function IsModuleSpacecraft(thisModule)
    ! Returns true if the module is really the spacecraft
    integer, intent(in) :: thisModule
    IsModuleSpacecraft = modules(thisModule)%spacecraft
  end function IsModuleSpacecraft

  ! ------------------------------------------------  MatchSignal  -----
  integer function MatchSignal ( Signals, Probe, sideband, channel, matchFlags )
    ! Given an array Signals, find the one in the array that provides
    ! the smallest superset of features of the signal Probe.  The result
    ! is zero if no signals match.
    ! If sideband or channel are present they are used for the probe, instead
    ! of the values in probe.

    type(signal_T), dimension(:), intent(in) :: Signals
    type(signal_T), intent(in) :: Probe
    integer, intent(in), optional :: sideband     ! Use this instead of probe%sideband
    logical, dimension(size(signals)), intent(out), optional :: matchFlags
    integer, intent(in), optional :: CHANNEL ! Just this channel

    integer :: BestMatch                ! The smallest number of 
    integer :: I                        ! Loop inductors, subscripts
    logical :: Match                    ! Channels in probe are in signal
    integer :: NumChannelsMatch
    integer :: mySideband               ! Either sideband or probe%sideband

    if ( present(matchFlags) ) matchFlags = .false.

    if ( present(sideband) ) then
      mySideband = sideband
    else
      mySideband = probe%sideband
    end if

    bestMatch = huge(bestMatch)
    matchSignal = 0
    do i = 1, size(signals)
      ! First, the signal must have the same band, instrument module,
      ! radiometer, spectrometer, spectrometer type and switch number as
      ! the probe signal
      if ( signals(i)%band /= probe%band .or. &
        &  signals(i)%instrumentModule /= probe%instrumentModule.or. &
        &  signals(i)%radiometer /= probe%radiometer .or. &
        &  signals(i)%sideband /= mySideband .or. &
        &  signals(i)%spectrometer /= probe%spectrometer .or. &
        &  signals(i)%spectrometerType /= probe%spectrometerType .or. &
        &  signals(i)%switch /= probe%switch ) cycle
      ! Now the channels in Probe all have to be present in Signals(i)
      if (present(channel)) then        ! User asked for a specific channel
        if ( associated(signals(i)%channels) ) then
          match = signals(i)%channels(channel)
        else
          match = .true.
        endif
      else
        match = (.not. associated(probe%channels)) .or. &
          & (.not. associated(signals(i)%channels) )
        if ( .not. match ) match = all( (probe%channels .and. &
          & signals(i)%channels(lbound(probe%channels,1):ubound(probe%channels,1)) ) &
          & .eqv. probe%channels )
      endif
      if ( match ) then
        if ( present( matchFlags ) ) matchFlags(i) = .true.
        if ( associated(signals(i)%channels) ) then
          numChannelsMatch = count(signals(i)%channels)
        else
          numChannelsMatch = size( signals(i)%frequencies )
        endif
        if ( numChannelsMatch < bestMatch ) then
          matchSignal = i
          bestMatch = numChannelsMatch
        end if
      end if
    end do
  end function MatchSignal

  ! ----------------------------------------- AreSignalsSubset ----------
  integer function AreSignalsSuperset ( signals, probe, sideband, channel )
    ! This is related to MatchSignal.  Given two arrays of signals: signals and
    ! probe, it returns -1 if probe is not a subset of signals, otherwise, it
    ! returns the number of channels by which signals is a superset of probe.
    ! If sideband and channel are present, they are taken to override the
    ! values in probe, if that makes sense.

    ! Dummy arguments
    type (Signal_T), dimension(:), intent(in) :: SIGNALS ! Signals to search
    type (Signal_T), dimension(:), intent(in) :: PROBE ! Potential subset
    integer, intent(in), optional :: SIDEBAND ! If present overrides one in probe
    integer, intent(in), optional :: CHANNEL ! If present overrides one in probe

    ! Local variables
    integer :: I                        ! Loop counter
    integer :: MATCH                    ! Result of call to match signal
    integer :: NOPBCHANNELS             ! Total channels in probe
    integer :: NOSIGCHANNELS            ! Total channels in signals

    ! Executable code

    AreSignalsSuperset = -1               ! Default to no
    noPbChannels = 0
    noSigChannels = 0
    
    ! First we identify whether they are a superset or not
    do i = 1, size(probe)
      ! Is this probe in signals?
      match = MatchSignal ( signals, probe(i), sideband=sideband, channel=channel )
      if ( match == 0 ) return

      ! OK, now how big is this probe
      if ( associated(probe(i)%channels ) ) then
        noPbChannels = noPbChannels + count( probe(i)%channels )
      else
        noPbChannels = noPbChannels + size ( probe(i)%frequencies )
      end if
    end do
  
    ! OK, if we got here, then we know we're a superset
    ! How many channels in signals
    do i = 1, size(signals)
      if ( associated(signals(i)%channels ) ) then
        noSigChannels = noSigChannels + count( signals(i)%channels )
      else
        noSigChannels = noSigChannels + size ( signals(i)%frequencies )
      end if
    end do

    ! Now return the difference
    AreSignalsSuperset = noSigChannels - noPbChannels

  end function AreSignalsSuperset

end module MLSSignals_M

! $Log$
! Revision 2.40  2001/09/17 22:53:46  livesey
! Added Instrument variable
!
! Revision 2.39  2001/05/16 23:05:06  livesey
! Added channel argument to AreSignalsSuperset
!
! Revision 2.38  2001/05/16 01:24:06  livesey
! Added AreSignalsSuperset routine
!
! Revision 2.37  2001/05/03 02:06:07  vsnyder
! Insert copyright notice
!
! Revision 2.36  2001/05/02 21:50:26  livesey
! Added initialization to deferred
!
! Revision 2.35  2001/05/02 19:12:11  vsnyder
! Nullify signal%channels so it won't get hosed by subsequent allocate_test
! Spiffify dump_signals and make it work if modules etc. are already gone
! Get label into name field of signals.
!
! Revision 2.34  2001/04/26 19:33:25  livesey
! Add sideband field to getSignalName etc.
!
! Revision 2.33  2001/04/26 02:33:03  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.32  2001/04/24 22:36:27  vsnyder
! Correct the 'what' argument to deallocate_test
!
! Revision 2.31  2001/04/23 23:12:45  vsnyder
! Finish adding 'time' command
!
! Revision 2.30  2001/04/23 23:10:37  vsnyder
! Add 'time' command
!
! Revision 2.29  2001/04/21 01:06:37  vsnyder
! Make Signal%Deferred initially false
!
! Revision 2.28  2001/04/21 01:05:26  vsnyder
! Deallocate Frequencies and Widths in DestroySignal if it's deferred
!
! Revision 2.27  2001/04/20 22:41:25  vsnyder
! Publish DestroySignal
!
! Revision 2.26  2001/04/19 20:29:24  livesey
! Added sideband argument to MatchSignal and made MatchSignal look at it or
! probe%sideband
!
! Revision 2.25  2001/04/13 23:54:06  livesey
! Change intent inout to intent out for matchFlags in MatchSignal
!
! Revision 2.24  2001/04/13 23:28:59  livesey
! Tidied up some gothcas in MatchSignal.
!
! Revision 2.23  2001/04/13 20:40:23  vsnyder
! Create GetNameOfSignal that returns the name of a Signal_T object (the
! already-present subroutine GetSignalName takes a signals-database index).
! Add DestroySignal to destroy one Signal_T object.
!
! Revision 2.22  2001/04/12 18:14:19  vsnyder
! Get the right vertex for 'deferred'
!
! Revision 2.21  2001/04/11 20:19:27  vsnyder
! Undo changes to 'deferred'
!
! Revision 2.20  2001/04/11 19:57:55  vsnyder
! OOPS! More work on 'deferred' spectrometers
!
! Revision 2.19  2001/04/11 19:54:00  vsnyder
! More work on 'deferred' spectrometers
!
! Revision 2.18  2001/04/11 18:31:04  vsnyder
! Change 'deferred' from boolean to numeric
!
! Revision 2.17  2001/04/10 23:20:05  livesey
! Added index field
!
! Revision 2.16  2001/04/10 17:59:53  vsnyder
! Remove sideband field from signal
!
! Revision 2.15  2001/04/09 20:30:46  vsnyder
! More work on MatchSignal
!
! Revision 2.14  2001/04/09 20:16:40  vsnyder
! Correct numChannelsMatch calculation
!
! Revision 2.13  2001/04/07 01:52:58  vsnyder
! Initial cut at MatchSignal
!
! Revision 2.12  2001/04/03 01:48:25  vsnyder
! No need to declare Make_Tree private -- It's internal!
!
! Revision 2.11  2001/03/29 23:52:31  vsnyder
! Added MaxSigLen parameter
!
! Revision 2.10  2001/03/28 19:51:58  vsnyder
! Remove "frequencies" and "widths" fields.  Use range instead in a "channels"
! field.  Put in units checking for every numeric field.
!
! Revision 2.9  2001/03/16 21:32:23  vsnyder
! Add 'Switches' test for dumping
!
! Revision 2.8  2001/03/16 02:10:32  vsnyder
! Put the sideband character in the correct place
!
! Revision 2.7  2001/03/16 01:54:47  vsnyder
! Use enumerated type instead of numbers for sideband; add a field for it
! in the "signal" spec.
!
! Revision 2.6  2001/03/16 00:30:49  vsnyder
! Make way for the pointing grid
!
! Revision 2.5  2001/03/15 21:02:07  vsnyder
! Cross-references between databases are by database index, not tree index
!
! Revision 2.4  2001/03/15 18:42:58  livesey
! Added GetSignal
!
! Revision 2.3  2001/03/15 18:39:42  vsnyder
! Periodic commit
!
! Revision 2.2  2001/03/14 23:44:47  vsnyder
! Correct a comment, other cosmetic changes in comments
!
! Revision 2.1  2001/03/14 02:05:52  vsnyder
! Moved MLSSignals_m to mlspgs/lib.
!
! Revision 2.11  2001/03/03 00:08:23  livesey
! Minor changes to module_t
!
! Revision 2.10  2001/03/02 01:29:31  livesey
! Added option of spacecraft module
!
! Revision 2.9  2001/02/28 21:44:21  livesey
! Moved back into L2, whoops!
!

! Tried to move this into lib, but failed, needs init_tables_module
! ----------------------------------------------------------
! Revision 2.7  2001/02/28 21:36:02  livesey
! Another plateau
!
! Revision 2.6  2001/02/28 01:16:11  livesey
! Interim version
!
! Revision 2.5  2001/02/27 01:28:34  vsnyder
! Rearranged a comment
!
! Revision 2.4  2001/02/27 00:23:48  livesey
! Very interim version
!
! Revision 2.3  2001/02/15 22:03:00  vsnyder
! Remove checking for ranges -- the type checker does it
!
