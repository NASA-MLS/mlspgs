module MLSSignals_M

  ! Process the MLSSignals section of the L2 configuration file and deal with
  ! parsing signal request strings.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use Declaration_Table, only: Range
  use DUMP_0, only: DUMP
  use Expr_M, only: Expr
  use Init_Tables_Module, only: f_band, f_centerFrequency, f_channel, &
    & f_channels, f_deferred, f_first, f_frequencies, f_frequency, &
    & f_last, f_lo, f_module, f_radiometer, f_spacecraft, f_spectrometer, &
    & f_spectrometerType, f_start, f_step, f_suffix, f_switch, &
    & f_width, f_widths, field_first, field_indices, field_last, &
    & s_band,  s_module, s_radiometer, s_spectrometerType, s_validSignal
  use Intrinsic, only: L_TRUE
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
             & MLSMSG_Error
  use Output_M, only: Output
  use String_Table, only: Display_String
  use Toggles, only: Gen, Toggle, Levels
  use Trace_M, only: Trace_begin, Trace_end
  use Tree, only: Decorate, Decoration, Node_ID, Nsons, Source_Ref, Sub_Rosa, &
    & Subtree
  use Tree_Types, only: N_named, N_Set_One

  implicit none

  ! =====     Defined Operators and Generic Identifiers     ==============
  
  interface DUMP
    module procedure DUMP_BANDS, DUMP_RADIOMETERS, DUMP_SIGNALS, &
      & DUMP_SPECTROMETERTYPES
  end interface

  ! This boring type defines a module
  type Module_T
    integer :: Name                     ! Sub_rosa index
    logical :: spaceCraft               ! Set if `module' is in fact s/c
  end type Module_T

  ! This type defines a radiometer.

  type Radiometer_T
    integer :: InstrumentModule         ! Tree index
    integer :: Prefix                   ! Sub_rosa index
    integer :: Suffix                   ! Sub_rosa index
    real(r8) :: LO                      ! Local oscillator in MHz
  end type Radiometer_T

  ! The second type describes a band within that radiometer

  type Band_T
    integer :: Prefix                   ! Sub_rosa index
    integer :: Radiometer               ! Tree index
    integer :: SpectrometerType         ! Tree index
    integer :: Suffix                   ! Sub_rosa index
    real(r8) :: CenterFrequency         ! Negative if not present (wide filter)
  end type Band_T

  ! This type gives the information for specific spectrometer families.  For
  ! all apart from the WF4 spectrometers, we list frequencies and widths. 
  ! Otherwise, the arrays are empty.

  type SpectrometerType_T
    integer :: Name                     ! Name for spectrometer type
    real(r8), pointer, dimension(:) :: Frequencies => NULL(), Widths => NULL()
  end type SpectrometerType_T

  ! This is the key type, it describes a complete signal (one band, or a
  ! subset of the channels in one band).  We have to main variables of this
  ! type flying around, see later.

  type Signal_T
    integer :: Band                     ! Tree index
    integer :: InstrumentModule         ! Tree index
    integer :: ParentSignal             ! Tree index of parent if derived, else 0
    integer :: Radiometer               ! Tree index
    integer :: SideBand                 ! -1:lower, +1:upper, 0:folded
    integer :: Spectrometer             ! Just a spectrometer number
    integer :: SpectrometerType         ! Tree index
    integer :: Switch                   ! Just a switch number

    real(r8) :: LO                      ! Radiometer local oscillator
    real(r8) :: CenterFrequency         ! Band local oscillator
    real(r8), POINTER, DIMENSION(:) :: Frequencies=>NULL() ! Mainly a shallow copy
    real(r8), POINTER, DIMENSION(:) :: Widths=>NULL() ! Mainly a shallow copy
    logical, POINTER, DIMENSION(:) :: Selected ! A bit field, or absent==all.
  end type Signal_T

  ! Now some databases, the first are fairly obvious.

  type(module_T), save, pointer, dimension(:) :: Modules => NULL()
  type(band_T), save, pointer, dimension(:) :: Bands => NULL()
  type(radiometer_T), save, pointer, dimension(:) :: Radiometers => NULL()
  type(spectrometerType_T), save, pointer, dimension(:) ::&
    & SpectrometerTypes => NULL()

  ! This array is the signals database.  The first entries are the official
  ! `valid' signals in the instrument.  Later one can derive things from that.
  ! for subsets of channels etc.
  type(signal_T), save, pointer, dimension(:) :: signals => NULL()

  !---------------------------- RCS Ident Info -------------------------------
  character (len=130), private :: Id = &
    & "$Id$"
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains
  ! -------------------------------------------------  MLSSignals  -----
  subroutine MLSSignals ( ROOT )
    ! Process the MLSSignals section of the L2 configuration file.
    integer, intent(in) :: ROOT         ! The "cf" vertex for the section

    type(band_T) :: Band                ! To be added to the database
    integer :: MyChannels               ! subtree index of field
    logical :: Deferred                 ! Set if frequencies/widths deferred
    integer :: Error                    ! Error level seen so far
    integer :: Field                    ! Field index -- f_something
    integer :: First                    ! "first" field of "spectrometer"
    integer :: Frequencies              ! subtree index of field
    logical :: Got(field_first:field_last)   ! "Got this field already"
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
    integer :: Units(2)                 ! of an expression
    double precision :: Value(2)        ! of an expression
    real(r8) :: Width                   ! "width" field of "spectrometer"
    integer :: Widths                   ! subtree index of field

    ! Error message codes
    integer, parameter :: AllOrNone = 1 ! All of a set of fields, or none
    integer, parameter :: AtLeastOne = allOrNone + 1   ! At least one
    integer, parameter :: BadMix = atLeastone + 1 ! Disallowed mixture of
    !                                     fields.
    integer, parameter :: Sizes = badMix + 1 ! Fields don't have same sizes

    error = 0
    if ( toggle(gen) ) call trace_begin ( "MLSSignals", root )
    do i = 2, nsons(root)-1 ! skip "MLSSignals" at "begin" and "end"
      son = subtree(i,root) ! A spec_args vertex now
      if ( node_id(son) == n_named ) then
        name = subtree(1, son)
        key = subtree(2, son)
      else
        name = 0
        key = son
      end if
      ! node_id(key) is now n_spec_args

      got = .false.
      select case ( decoration(subtree(1,decoration(subtree(1,key)))) )


      case ( s_band ) ! ...................................  BAND  .....
        band%prefix= sub_rosa(name)
        band%centerFrequency = -1.0_r8 ! "The 'frequency' field is absent"
        do j = 2, nsons(key)
          son = subtree(j,key)
          field = decoration(subtree(1,son))
          gson = subtree(2,son)
          got(field) = .true.
          select case ( field )
          case ( f_centerFrequency )
            call expr ( gson, units, value )
            band%centerFrequency = value(1)
          case ( f_suffix )
            band%suffix = sub_rosa(gson)
          case ( f_radiometer )
            band%radiometer = gson
          case ( f_spectrometerType )
            band%spectrometerType= gson
          case default
            ! Shouldn't get here if the type checker worked
          end select
        end do ! i = 2, nsons(key)
        call decorate ( key, addBandToDatabase ( bands, band ) )


      case ( s_module ) ! .............................. MODULE ........
        thisModule%name = sub_rosa(name)
        thisModule%spaceCraft = .false.
        do j = 2,nsons(key)
          son = subtree(j,key)
          field = decoration(subtree(1,son))
          if (nsons(son) > 1) then
            gson = subtree(2,son)
          else
            gson = son
          end if
          got(field) = .true.
          select case ( field )
          case (f_spaceCraft)
            if (node_id(son) == n_set_one) then
              deferred=.true.
            else
              deferred=decoration(gson) == l_true
            endif
          case default
            ! Shouldn't get here if parser worked
          end select
        end do
        call decorate ( key, AddModuleToDatabase (modules, thisModule ) )


      case ( s_radiometer ) ! .......................  RADIOMETER  .....
        radiometer%prefix= sub_rosa(name)
        do j = 2, nsons(key)
          son = subtree(j,key)
          field = decoration(subtree(1,son))
          gson = subtree(2,son)
          select case ( field )
          case ( f_lo )
            call expr ( gson, units, value )
            radiometer%lo = value(1)
          case ( f_suffix )
            radiometer%suffix = sub_rosa(gson)
          case ( f_module )
            radiometer%instrumentModule=gson
          case default
            ! Shouldn't get here if the type checker worked
          end select
        end do ! i = 2, nsons(key)
        call decorate ( key, addRadiometerToDatabase ( radiometers, radiometer ) )


      case ( s_validSignal ) ! ..........................  VALIDSIGNAL  .....
        do j = 2, nsons(key)
          son = subtree(j,key)
          field = decoration(subtree(1,son))
          gson = subtree(2,son)
          got(field)=.true.
          select case ( field )
          case ( f_band )
            signal%band = gson
          case ( f_spectrometer )
            call expr ( gson, units, value )
            signal%spectrometer = value(1)
          case ( f_switch )
            call expr ( gson, units, value )
            signal%switch = value(1)
          case ( f_frequencies )
            frequencies=son
          case ( f_widths ) 
            widths=son
          case default
            ! Shouldn't get here if the type checker worked
          end select
        end do ! i = 2, nsons(key)
        ! Set default values for remaining parameters
        signal%sideband=0
        signal%parentSignal=0
        signal%radiometer=bands( &
          & decoration(decoration(signal%band)))%radiometer
        signal%lo=radiometers( &
          & decoration(decoration(signal%radiometer)))%lo
        signal%instrumentModule=radiometers( &
          & decoration(decoration(signal%radiometer)))%instrumentModule
        signal%spectrometerType=bands( &
          & decoration(decoration(signal%band)))%spectrometerType
        signal%centerFrequency=bands( &
          & decoration(decoration(signal%band)))%centerFrequency
        nullify(signal%selected)
        ! For the wide filters, we specify frequency etc. here.
        if ( got(f_frequencies) .neqv. got(f_widths) ) call announceError ( &
          & allOrNone, f_frequencies, (/ f_widths /) )
        if ( got(f_frequencies) ) then
          if ( nsons(frequencies) /= nsons(widths) ) &
            call announceError ( sizes, f_frequency, (/ f_widths /) )
          if ( error == 0 ) then
            call allocate_Test ( signal%frequencies, nsons(son)-1, &
              & 'signal%frequencies', moduleName )
            call allocate_Test ( signal%widths, nsons(son)-1, &
              & 'signal%widths', moduleName)
            do k = 2, nsons(frequencies)
              call expr ( subtree(k,frequencies), units, value )
              signal%frequencies(k-1) = value(1)
              call expr ( subtree(k,widths), units, value )
              signal%widths(k-1) = value(1)
            end do
          end if
        else
          signal%frequencies=> spectrometerTypes( &
            & decoration(decoration(signal%spectrometerType)))%frequencies
          signal%widths=> spectrometerTypes( &
            & decoration(decoration(signal%spectrometerType)))%widths
        end if
        call decorate ( key, addSignalToDatabase ( signals, signal ) )
        ! Now nullify pointers so they don't get hosed later
        nullify ( signal%frequencies )
        nullify ( signal%widths )
        nullify ( signal%selected )


      case ( s_spectrometerType ) ! .............  SPECTROMETERTYPE .....
        spectrometerType%name=sub_rosa(name)
        deferred=.false.
        first = 0
        do j = 2, nsons(key)
          son = subtree(j,key)
          field = decoration(subtree(1,son))
          if (nsons(son) > 1) then
            gson = subtree(2,son)
          else
            gson = son                  ! For case of /deferred
          end if
          got(field) = .true.
          select case ( field )
          case ( f_first, f_last, f_start, f_step, f_width )
            call expr ( gson, units, value )
            select case ( field )
            case ( f_first )
              first = value(1)
            case ( f_last )
              last = value(1)
            case ( f_start )
              start = value(1)
            case ( f_step )
              step = value(1)
            case ( f_width )
              width = value(1)
            end select
          case ( f_frequencies )   ! Postpone processing until later, so that
            frequencies = son      ! we can verify that Frequencies and
          case ( f_widths )        ! Widths have the same number of values
            widths = son
          case ( f_deferred )
            if (node_id(son) == n_set_one) then
              deferred=.true.
            else
              deferred=decoration(gson) == l_true
            endif
          case default
            ! Shouldn't get here if the type checker worked
          end select
        end do ! i = 2, nsons(key)
        if ( got(f_frequencies) .neqv. got(f_widths) ) call announceError ( &
          & allOrNone, f_frequencies, (/ f_widths /) )
        if ( got(f_frequencies) ) then
          if ( any(got( (/ f_last, f_start, f_step, f_width /) )) ) &
            & call announceError ( badMix, f_frequencies, &
            & (/ f_last, f_start, f_step, f_width /) )
          if ( nsons(frequencies) /= nsons(widths) ) &
            call announceError ( sizes, f_frequency, (/ f_widths /) )
          if ( error == 0 ) then
            call allocate_Test ( spectrometerType%frequencies, nsons(son)-first, &
              & 'spectrometerType%frequencies', moduleName, lowBound = first )
            call allocate_Test ( spectrometerType%widths, nsons(son)-first, &
              & 'spectrometerType%widths', moduleName, lowBound = first )
            do k = 2, nsons(frequencies)
              call expr ( subtree(k,frequencies), units, value )
              spectrometerType%frequencies(k-2+first) = value(1)
              call expr ( subtree(k,widths), units, value )
              spectrometerType%widths(k-2+first) = value(1)
            end do
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
          nullify(spectrometerType%frequencies)
          nullify(spectrometerType%widths)
        end if
        if ( .not. any(got( (/ f_frequencies, f_start /) )) .and. &
          & (.not. deferred) ) &
          & call announceError ( atLeastOne, f_frequencies, (/ f_start /) )
        if ( error == 0 ) call decorate ( key, addSpectrometerTypeToDatabase ( &
          & spectrometerTypes, spectrometerType ) )

        ! Nullify pointers to temporary stuff so it doesn't get hosed later
        nullify ( spectrometerType%frequencies )
        nullify ( spectrometerType%widths )
      case default
        ! Shouldn't get here if the type checker worked
      end select

    end do

    if ( toggle(gen) ) call trace_end ( "MLSSignals" )

    if ( levels(gen) > 0 ) then
      call dump(radiometers)
      call dump(spectrometerTypes)
      call dump(bands)
      call dump(signals)
    endif

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
      call output ( 'At line '  )
      call output ( mod(source,256) )
      call output ( ', column ' )
      call output ( source/256 )
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
        call output ( ' shall appear, or none of them shall.', advance='yes' )
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
          CALL display_string ( field_indices(moreFields(i)) )
        end do ! i
        call output ( ' shall appear.', advance='yes' )
      case ( sizes )
        call output ( 'The fields ' )
        do i = 1, size(moreFields)
          if ( i == size(moreFields) ) then
            call output ( ' and ' )
          else
            call output ( ', ' )
          end if
          call display_string ( field_indices(moreFields(i)) )
        end do ! i
        call output ( ' shall all have the same size.', advance='yes' )
      end select
    end subroutine AnnounceError

  end subroutine MLSSignals

  ! ----------------------------------  AddBandToDatabase  -----
  integer function AddBandToDatabase ( Database, Item )
    type(band_T), dimension(:), pointer :: Database
    type(band_T), intent(in) :: Item

    ! Local variables
    type (Band_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddBandToDatabase = newSize
  end function AddBandToDatabase

  ! ----------------------------------  AddModuleToDatabase  -----
  integer function AddModuleToDatabase ( Database, Item )
    type(module_T), dimension(:), pointer :: Database
    type(module_T), intent(in) :: Item

    ! Local variables
    type (Module_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddModuleToDatabase = newSize
  end function AddModuleToDatabase

  ! ----------------------------------  AddRadiometerToDatabase  -----
  integer function AddRadiometerToDatabase ( Database, Item )
    type(radiometer_T), dimension(:), pointer :: Database
    type(radiometer_T), intent(in) :: Item

    ! Local variables
    type (Radiometer_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddRadiometerToDatabase = newSize
  end function AddRadiometerToDatabase

  ! ----------------------------------  AddSignalToDatabase  -----
  integer function AddSignalToDatabase ( Database, Item )
    type(signal_T), dimension(:), pointer :: Database
    type(signal_T), intent(in) :: Item

    ! Local variables
    type (Signal_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddSignalToDatabase = newSize
  end function AddSignalToDatabase

  ! -------------------------------  AddSpectrometerTypeToDatabase  -----
  integer function AddSpectrometerTypeToDatabase ( Database, Item )
    type(spectrometerType_T), dimension(:), pointer :: Database
    type(spectrometerType_T), intent(in) :: Item

    ! Local variables
    type (spectrometerType_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddSpectrometerTypeToDatabase = newSize
  end function AddSpectrometerTypeToDatabase

  ! --------------------------------  DestroyBandDatabase  -----
  subroutine DestroyBandDatabase ( Bands )
    type(band_T), dimension(:), pointer :: Bands
    integer :: Status
    if ( associated(bands) ) then
      deallocate ( bands, stat = status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_DeAllocate // 'Band database' )
    end if
  end subroutine DestroyBandDatabase

  ! --------------------------------  DestroyModuleDatabase  -----
  subroutine DestroyModuleDatabase ( Modules )
    type(module_T), dimension(:), pointer :: Modules
    integer :: Status
    if ( associated(modules) ) then
      deallocate ( modules, stat = status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_DeAllocate // 'Module database' )
    end if
  end subroutine DestroyModuleDatabase

  ! --------------------------------  DestroyRadiometerDatabase  -----
  subroutine DestroyRadiometerDatabase ( Radiometers )
    type(radiometer_T), dimension(:), pointer :: Radiometers
    integer :: Status
    if ( associated(radiometers) ) then
      deallocate ( radiometers, stat = status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_DeAllocate // 'Radiometer database' )
    end if
  end subroutine DestroyRadiometerDatabase

  ! --------------------------------  DestroySignalDatabase  -----
  subroutine DestroySignalDatabase ( Signals )
    type(signal_T), dimension(:), pointer :: Signals
    integer :: Status
    if ( associated(signals) ) then
      deallocate ( signals, stat = status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_DeAllocate // 'Signal database' )
    end if
  end subroutine DestroySignalDatabase

  ! ------------------------------------  DestroySpectrometerType  -----
  subroutine DestroySpectrometerType ( SpectrometerType )
    type(spectrometerType_T) :: SpectrometerType
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

  ! ------------------------------------------------ DumpBands -----
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
      call output ( decoration(bands(i)%radiometer))
      call output ( bands(i)%radiometer )
      call output ( ' - ' )
      call display_string ( sub_rosa(bands(i)%radiometer), advance='yes' )
      call output ( '   SpectrometerType: ')
      call output ( decoration(bands(i)%spectrometerType) )
      call output ( ' - ' )
      call display_string ( sub_rosa(bands(i)%spectrometerType),advance='yes' )
      call output ( '   Frequency: ')
      call output ( bands(i)%centerFrequency,advance='yes' )
    end do
  end subroutine DUMP_BANDS

  ! ----------------------------------------------- Dump_Radiometers --
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
      call output ( decoration(radiometers(i)%instrumentModule) )
      call output ( ' - ' )
      call display_string ( sub_rosa(radiometers(i)%instrumentModule), advance='yes' ) 
      call output ( '   LO: ')
      call output ( radiometers(i)%lo,advance='yes' )
    end do
  end subroutine DUMP_RADIOMETERS

  ! ------------------------------------------- DumpSpectrometerTypes --
  subroutine DUMP_SIGNALS ( SIGNALS )
    type (signal_T), intent(in) :: SIGNALS(:)
    integer :: i
    call output ( 'SIGNALS: SIZE = ')
    call output ( size(signals), advance='yes' )
    do i = 1, size(signals)
      call output ( i, advance='yes' )
      call output ( '   Module: ')
      call output ( decoration(signals(i)%instrumentModule ) )
      call output ( ' - ' )
      call display_string ( sub_rosa(signals(i)%instrumentModule), advance='yes' )
      call output ( '   Radiometer: ')
      call output ( decoration( signals(i)%radiometer ) )
      call output ( ' - ' )
      call display_string ( sub_rosa(signals(i)%radiometer), advance='yes' )
      call output ( '   First LO: ')
      call output ( signals(i)%lo, advance='yes')
      call output ( '   Band: ')
      call output ( decoration(signals(i)%band ) )
      call output (' - ')
      call display_string ( sub_rosa(signals(i)%band),advance='yes' )
      call output ( '   Band center frequency: ')
      call output ( signals(i)%centerFrequency, advance='yes')
      call output ( '   SpectrometerType: ')
      call output ( decoration( signals(i)%spectrometerType ) )
      call output ( ' - ' )
      call display_string ( sub_rosa(signals(i)%spectrometerType),advance='yes' )
      call output ( '   Channels: ' )
      call output ( lbound(signals(i)%frequencies,1), 3 )
      call output ( ':' )
      call output ( ubound(signals(i)%frequencies,1), 3, advance='yes' )
      call output ( '   Frequencies:', advance='yes' )
      call dump ( signals(i)%frequencies )
      call output ( '   Widths:', advance='yes' )
      call dump ( signals(i)%widths )
      if (associated ( signals(i)%selected) ) then
        call output ( '   Channel selected flags:', advance='yes' )
        call dump ( signals(i)%selected )
      else
        call output ( '   All channels selected.', advance='yes' )
      end if
    end do
  end subroutine DUMP_SIGNALS

  ! ------------------------------------------- DumpSpectrometerTypes --
  subroutine DUMP_SPECTROMETERTYPES ( SPECTROMETERTYPES )
    type (SpectrometerType_T), intent(in) :: SPECTROMETERTYPES(:)
    integer :: i
    call output ( 'SPECTROMETERTYPES: SIZE = ')
    call output ( size(spectrometerTypes), advance='yes' )
    do i = 1, size(spectrometerTypes)
      call output ( i,1 )
      call output ( ': ')
      call display_string (spectrometerTypes(i)%name,advance='yes')
      if (associated(spectrometerTypes(i)%frequencies)) then
        call output ( '  Channels: ' )
        call output ( lbound(spectrometerTypes(i)%frequencies,1), 3 )
        call output ( ':' )
        call output ( ubound(spectrometerTypes(i)%frequencies,1), 3, advance='yes' )
        call output ( '  Frequencies:', advance='yes' )
        call dump ( spectrometerTypes(i)%frequencies )
        call output ( '  Widths:', advance='yes' )
        call dump ( spectrometerTypes(i)%widths )
      else
        call output ('   Frequencies and widths deferred.', advance='yes' )
      end if
    end do
  end subroutine DUMP_SPECTROMETERTYPES
end module MLSSignals_M

! $Log$
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
