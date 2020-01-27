! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module MLSSignals_M

  ! Process the MLSSignals section of the L2 configuration file and deal with
  ! parsing signal request strings.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use Dump_0, only: Dump
  use Expr_M, only: Expr
  use HighOutput, only: HeadLine, OutputNamedValue, OutputTable, Tab
  use Init_MLSSignals_M ! Everything
  use Intrinsic, only: Field_First, Field_Indices, Lit_Indices, &
    & Phyq_Dimensionless, Phyq_Frequency, Phyq_Indices, S_Time, L_A, L_EMLS
  use MLSKinds, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, PVMErrorMessage
  use MLSStrings, only: Lowercase, Capitalize
  use Output_M, only: Blanks, Newline, Output
  use String_Table, only: Display_String, Get_String

  implicit none

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

! AddBandToDatabase               As the name implies
! AddModuleToDatabase             ...
! AddRadiometerToDatabase         ...
! AddSignalToDatabase             ...
! AddSpectrometerTypeToDatabase   ...
! AreSignalsSuperset              ...
! DestroyBandDatabase             ...
! DestroyModuleDatabase           ...
! DestroyRadiometerDatabase       ...
! DestroySignal                   ...
! DestroySignalDatabase           ...
! DestroySpectrometerType         ...
! DestroySpectrometerTypeDatabase ...
! DisplayRadiometer               ... given its index in radiometer database
! DisplaySignalName               ...
! DisplaySignalName_index         ... given a signal index
! DisplaySignalName_signal        ... given a signal structure
! Dump                            ...
! Dump_All                        ...
! Dump_Bands                      ...
! Dump_Modules                    ...
! Dump_Radiometers                ...
! Dump_Signal                     ...
! Dump_Signals                    ...
! Dump_Spectrometertypes          ...
! GetAllModules                   Return tree nodes for all modules
! GetBandName                     Given an index in the Bands database, place band name in string
! GetModuleFromRadiometer         Returns module field from given radiometer given as database index
! GetModuleIndex                  Returns the index in the module database, given module name in mixed case
! GetModuleFromSignal             Returns module field from given signal given as database index
! GetModuleName                   Given the index in the module database, returns module name in mixed case
! GetNameOfSignal                 Given a signal object, this routine constructs a full signal name
! GetRadiometerFromSignal         Returns radiometer field from given signal given as database index
! GetRadiometerName               Given an index in the Radiometers database, place radiometer name
! GetRadiometerIndex              Returns index in the Radiometers database, given radiometer name
! GetSidebandStartStop            Given a signal, compute SidebandStart and SidebandStop
! GetSignal                       Given the database index, this routine returns the signal data structure
! GetSignalIndex                  Returns the index in the signals database, given signal name in mixed case
! GetSignalName                   Given an index in the signals database, constructs full signal name.
! GetSpectrometerTypeName         Place spectrometer name and number in string
! IsAnyModuleSpacecraft           Returns true if any module is really the spacecraft
! IsModuleSpacecraft              Returns true if the module is really the spacecraft
! IsSpacecraftAura                Returns true if the s/c is really Aura
! MatchSignal                     Given an array Signals, find the matching one
! MatchSignalPair                 Determine whether two signals match
! MLSSignals                      Process the MLSSignals section of the L2 configuration file
! === (end of toc) ===

! === (start of api) ===
! === (end of api) ===

  private ! So as not to re-export everything accessed by USE association.

  ! Public procedures and interfaces:
  public :: AddbandtoDatabase, AddModuletoDatabase, AddradiometertoDatabase
  public :: AddsignaltoDatabase, AddspectrometertypetoDatabase, Aresignalssuperset
  public :: DestroybandDatabase, DestroyModuleDatabase
  public :: DestroyradiometerDatabase, Destroysignal, DestroysignalDatabase
  public :: Destroyspectrometertype, DestroyspectrometertypeDatabase
  public :: Displayradiometer, DisplaysignalName
  public :: DisplaysignalName_Index, DisplaysignalName_Signal
  public :: Dump, Dump_Bands, Dump_Radiometers, Dump_Signal, Dump_Signals
  public :: Dump_All, Dump_Modules, Dump_Spectrometertypes
  public :: GetallModules, GetbandName, GetfirstChannel, GetModulefromradiometer
  public :: GetModuleindex, Getsidebandloop, Getsidebandstartstop, Getsignalindex
  public :: GetModulefromsignal, GetModuleName, GetNameofsignal
  public :: Getradiometerfromsignal, GetradiometerName, Getradiometerindex
  public :: Getsignal, GetsignalName, GetspectrometertypeName
  public :: Isspacecraftaura, IsanyModulespacecraft, IsModulespacecraft
  public :: Matchsignal, Matchsignalpair, MLSSignals
  public :: Pvmpacksignal, Pvmunpacksignal

  integer, public, parameter :: MaxSigLen = 80 ! Maximum length of a signal name

  interface DisplaySignalName
    module procedure DisplaySignalName_index, DisplaySignalName_signal
  end interface

  interface MatchSignal
    module procedure MatchSignals, MatchSignalPair
  end interface

  ! =====     Defined Operators and Generic Identifiers     ==============
  
  interface Dump
    module procedure Dump_Bands, Dump_Modules, Dump_Oneradiometer, &
      & Dump_Radiometers, Dump_Signal, Dump_Signals, &
      & Dump_Spectrometertype, Dump_Spectrometertypes
  end interface

  ! This boring type defines a module
  type, public :: Module_T
    integer :: Name                     ! Sub_rosa index of declaration's label
    integer :: Node                     ! Node of tree where module declared
    logical :: spaceCraft               ! Set if module is in fact s/c
    integer :: supportedModule          ! See discssion below
    logical :: Aura  = .true.           ! Set if s/c is in fact Aura
    character(len=32) :: NameString = ' '  ! If not get_string (e.g., camelCase)
  end type Module_T
  ! A-SMLS does things in a more complicated way, as the timing of the
  ! different spectrometer modules varies.  This means that there need
  ! to be multiple Spacecraft modules, with timings that are on the
  ! same cadence as the corresponding spectrometer.  In these cases,
  ! supportedModule is an index into the corresponding module.  That
  ! is to say for each Spectrometer there will be two Module_Ts: one
  ! containing the tangent pointing information as usual, the other
  ! containing the spacecraft information on the same timing, for
  ! which the "supportedModule" entry will point to the tangent
  ! information

  ! This type defines a radiometer.

  type, public :: Radiometer_T
    real(r8) :: LO                      ! Local oscillator in MHz
    integer :: InstrumentModule         ! Index in Modules database
    integer :: Polarization             ! L_A or L_B, default is L_A
    integer :: Prefix                   ! Sub_rosa index of declaration's label
    integer :: Suffix                   ! Sub_rosa index
    integer :: SingleSideband           ! +/-1 indicates indicates single sideband 0 folded
  end type Radiometer_T

  ! The second type describes a band within that radiometer

  type, public :: Band_T
    real(r8) :: CenterFrequency         ! Zero if not present (wide filter)
    integer :: Prefix                   ! Sub_rosa index of declaration's label
    integer :: Radiometer               ! Index in Radiometers database, or none if deferred
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
    logical :: DACS=.false.             ! Set if this spectrometer is a DACS
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
    logical :: DACS = .false.           ! This signal is a DACS
    logical :: Deferred = .false.       ! "Frequencies/widths are deferred"
    integer :: Direction                ! +1 Channel 1 closest to lo, -1 reverse.
    integer :: Index                    ! Index into master signals database
    integer :: InstrumentModule         ! Index in Modules database
    integer :: Name                     ! Sub_rosa index of declaration's label
    integer :: Radiometer               ! Index in Radiometers database
    integer :: SideBand                 ! -1=lower, +1=upper, 0=folded
    integer :: SingleSideband           ! Indicates only sideband for this radiometer
    integer :: Spectrometer             ! Just a spectrometer number
    integer :: SpectrometerType         ! Index in SpectrometerTypes database
    integer :: Switch                   ! Just a switch number
  end type Signal_T

  ! Now some databases, the first are fairly obvious.
  !??? Should these be public ???

  type(module_T), public, save, pointer, dimension(:)  :: Modules => null()   
  type(band_T), public, save, pointer, dimension(:)    :: Bands => null()     
  type(radiometer_T), public, save, pointer, dimension(:) &
    &                                                  :: Radiometers => null()
  type(spectrometerType_T), public, save, pointer, dimension(:) &
    &                                                  :: Spectrometertypes => null()

  ! This array is the signals database.  The first entries are the official
  ! `valid' signals in the instrument.  Later one can derive things from that.
  ! for subsets of channels etc.
  type(signal_T), public, save, pointer, dimension(:)  :: Signals => null()           
  integer, public, save                                :: Instrument = l_emls         
  integer, parameter                                   :: Maxradiometernamelen = 16   

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! -------------------------------------------------  MLSSignals  -----
  subroutine MLSSignals ( ROOT )
    ! Process the MLSSignals section of the L2 configuration file.

    use MLSStringLists, only: Switchdetail
    use MoreTree, only: Get_Boolean, Get_Label_And_Spec, StarterrorMessage
    use Next_Tree_Node_M, only: Next_Tree_Node, Next_Tree_Node_State
    use Time_M, only: Time_Now
    use Toggles, only: Gen, Levels, Switches, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    use Tree, only: Decorate, Decoration, Nsons, Sub_Rosa, Subtree

    integer, intent(in) :: ROOT         ! The "cf" vertex for the section

    type(band_T) :: Band                ! To be added to the database
    integer :: Channels                 ! subtree index of field
    logical :: DACS                     ! This spectrometer type is a DACS
    logical :: Deferred                 ! "Frequencies/widths are deferred"
    integer :: Error                    ! Error level seen so far
    integer :: Field                    ! Field index -- f_something
    integer :: First                    ! "first" field of "spectrometer"
    logical :: Got(field_first:last_Signal_Field)   ! "Got this field already"
    integer :: Gson                     ! Son of Son.
    integer :: J, K                     ! Subscript and loop inductor.
    integer :: Key                      ! Indexes the spec_args vertex.
    integer :: Last                     ! "last" field of "spectrometer"
    integer :: Me = -1                  ! String index for trace
    integer :: Me_Band = -1             ! String index for trace
    integer :: Me_Module = -1           ! String index for trace
    integer :: Me_Radiometer = -1       ! String index for trace
    integer :: Me_Signal = -1           ! String index for trace
    integer :: Me_SpectrometerType = -1 ! String index for trace
    integer :: Name                     ! sub_rosa of label of specification,
                                        ! if any, else zero.
    type(radiometer_T) :: Radiometer    ! To be added to the database
    integer :: Son                      ! Some subtree of root.
    type(signal_T) :: Signal            ! To be added to the database
    type(spectrometerType_T) :: SpectrometerType ! To be added to the database
    real(r8) :: Start                   ! "start" field of "spectrometer"
    type(next_tree_node_state) :: State ! of tree traverser
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
    integer, parameter :: UnneededRadiometer = WrongUnits + 1 ! Radiomter spec in wrong place
    integer, parameter :: NoDeferredRadiometer = UnneededRadiometer + 1 ! Radiometer spec needed
    integer, parameter :: PriorTrouble = NoDeferredRadiometer + 1
    integer, parameter :: noModuleSpecified = PriorTrouble + 1 ! No module in either signal or radimoeter
    integer, parameter :: conflictingModules = noModuleSpecified + 1 ! Module in signal and radiometer differ
    integer, parameter :: notSpacecraft = conflictingModules + 1 ! Had supportedModule set innappropriately

    error = 0
    timing = .false.
    call trace_begin ( me, "MLSSignals", root, cond=toggle(gen) )
    do
      son = next_tree_node(root,state)
      if ( son == 0 ) exit
      call get_label_and_spec ( son, name, key )
      ! node_id(key) is now n_spec_args

      got = .false.
      select case ( decoration(subtree(1,decoration(subtree(1,key)))) )

      case ( s_band ) ! ...................................  BAND  .....
        call trace_begin ( me_band, "MLSSignals.band", son, &
          & cond=toggle(gen) .and. levels(gen) > 0 )
        band%prefix = name
        band%centerFrequency = 0.0_r8 ! "The 'frequency' field is absent"
        band%radiometer = 0             ! No radiometer defined at first.
        do j = 2, nsons(key)
          son = subtree(j,key)
          field = decoration(subtree(1,son))
          gson = subtree(2,son)
          got(field) = .true.
          select case ( field )
          case ( f_centerFrequency )
            ! Front end should be checking units now, but just in case....
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
        end do ! j = 2, nsons(key)
        call decorate ( key, addBandToDatabase ( bands, band ) )
        call trace_end ( "MLSSignals.band", &
          & cond=toggle(gen) .and. levels(gen) > 0 )

      case ( s_module ) ! ............................  MODULE  ........
        call trace_begin ( me_module, "MLSSignals.module", son, &
          & cond=toggle(gen) .and. levels(gen) > 0 )
        thisModule%name = name
        thisModule%spaceCraft = .false.
        thisModule%node = decoration(name)
        thisModule%supportedModule = 0
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
          case (f_Aura)
            thisModule%Aura = get_boolean(son)
          case (f_instrument)
            instrument = decoration(gson)
          case (f_spaceCraft)
            thisModule%spacecraft = get_boolean(son)
          case (f_supportedModule)
            thisModule%supportedModule = decoration(decoration(gson))
          case default
            ! Shouldn't get here if parser worked
          end select
        end do
        if ( .not. thisModule%spacecraft .and. thisModule%supportedModule /= 0 ) &
          call announceError ( notSpacecraft )
        call decorate ( key, AddModuleToDatabase (modules, thisModule ) )
        call trace_end ( "MLSSignals.module", &
          & cond=toggle(gen) .and. levels(gen) > 0 )

      case ( s_radiometer ) ! .......................  RADIOMETER  .....
        call trace_begin ( me_radiometer, "MLSSignals.radiometer", son, &
          & cond=toggle(gen) .and. levels(gen) > 0 )
        radiometer%polarization = l_a
        radiometer%prefix = name
        radiometer%singleSideband = 0
        radiometer%instrumentModule = 0
        do j = 2, nsons(key)
          son = subtree(j,key)
          field = decoration(subtree(1,son))
          gson = subtree(2,son)
          select case ( field )
          case ( f_lo )
            ! Front end should be checking units now, but just in case....
            call expr_check ( gson, units, value, field, phyq_frequency )
            radiometer%lo = value(1)
          case ( f_module )
            radiometer%instrumentModule = decoration(decoration(gson))
          case ( f_polarization )
            radiometer%polarization = decoration(gson)
          case ( f_singlesideband )
            ! Front end should be checking units now, but just in case....
            call expr_check ( gson, units, value, field, phyq_dimensionless )
            radiometer%singleSideband = nint ( value(1) )
          case ( f_suffix )
            radiometer%suffix = sub_rosa(gson)
          case default
            ! Shouldn't get here if the type checker worked
          end select
        end do ! j = 2, nsons(key)
        call decorate ( key, addRadiometerToDatabase ( radiometers, radiometer ) )
        call trace_end ( "MLSSignals.radiometer", &
          & cond=toggle(gen) .and. levels(gen) > 0 )

      case ( s_signal ) ! ..........................  VALIDSIGNAL  .....
        call trace_begin ( me_signal, "MLSSignals.signal", son, &
          & cond=toggle(gen) .and. levels(gen) > 0 )
        signal%sideband = 0
        signal%radiometer = 0
        signal%instrumentModule = 0
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
          case ( f_direction )
            ! Front end should be checking units now, but just in case....
            call expr_check ( gson, units, value, field, phyq_dimensionless )
            signal%direction = nint(value(1))
          case ( f_radiometer )
            signal%radiometer = decoration(decoration(gson))
          case ( f_module )
            signal%instrumentModule = decoration(decoration(gson))
          case ( f_spectrometer )
            ! Front end should be checking units now, but just in case....
            call expr_check ( gson, units, value, field, phyq_dimensionless )
            signal%spectrometer = value(1)
          case ( f_switch )
            ! Front end should be checking units now, but just in case....
            call expr_check ( gson, units, value, field, phyq_dimensionless )
            signal%switch = value(1)
          case default
            ! Shouldn't get here if the type checker worked
          end select
        end do ! j = 2, nsons(key)
        signal%name = name
        ! Did this band have a radiometer spec deferred to this point?
        if ( bands(signal%band)%radiometer == 0 ) then
          ! The radiometer was deferred in the band definition
          if ( signal%radiometer == 0 ) call announceError ( noDeferredRadiometer )
        else
          ! Didn't need to supply radiometer then.
          if ( got(f_radiometer) ) call announceError ( unneededRadiometer )
          signal%radiometer = bands(signal%band)%radiometer
        end if
        if ( associated(radiometers) .and. associated(bands) .and. &
           & associated(spectrometerTypes) ) then
          signal%lo = radiometers(signal%radiometer)%lo
          ! The instrument module may have been supplied here, or as part of the radiometer specification
          ! either way, make sure they're consistent
          if ( signal%instrumentModule == 0 ) then
            if ( radiometers(signal%radiometer)%instrumentModule == 0 ) &
              & call announceError ( noModuleSpecified )
            signal%instrumentModule = radiometers(signal%radiometer)%instrumentModule
          else
            if ( ( radiometers(signal%radiometer)%instrumentModule /= 0 ) .and. &
                 ( radiometers(signal%radiometer)%instrumentModule /= signal%instrumentModule ) ) &
              & call announceError ( conflictingModules )
          end if
          signal%spectrometerType = bands(signal%band)%spectrometerType
          signal%singleSideband = radiometers(signal%radiometer)%singleSideband
          signal%centerFrequency = bands(signal%band)%centerFrequency
          signal%deferred = spectrometerTypes(signal%spectrometerType)%deferred
          signal%dacs = spectrometerTypes(signal%spectrometerType)%dacs
          if ( signal%deferred .neqv. got(f_channels) ) &
            & call announceError ( deferredChannels )
          ! For the wide filters, we specify frequency etc. here.
          if ( got(f_channels) ) then
            if ( error == 0 ) then
              call allocate_Test ( signal%frequencies, nsons(channels)-1, &
                & 'signal%frequencies', moduleName )
              call allocate_Test ( signal%widths, nsons(channels)-1, &
                & 'signal%widths', moduleName)
              do k = 2, nsons(channels)
                call expr ( subtree(k,channels), units, value )
                signal%frequencies(k-1) = value(1)
                signal%widths(k-1) = value(2)
              end do
              if ( any(units /= phyq_frequency) ) &
                ! Front end should be checking units now, but just in case....
                & call announceError ( wrongUnits, f_channels, &
                  & (/ phyq_frequency /) )
            end if
          else
            signal%frequencies => spectrometerTypes(signal%spectrometerType)% &
              & frequencies
            signal%widths => spectrometerTypes(signal%spectrometerType)%widths
          end if
          call decorate ( key, addSignalToDatabase ( signals, signal ) )
          signals(size(signals))%index = size(signals)
        else
          call announceError ( priorTrouble )
        end if
        ! Now nullify pointers so they don't get hosed later by allocate_test
        nullify ( signal%frequencies )
        nullify ( signal%widths )
        call trace_end ( "MLSSignals.signal", &
          & cond=toggle(gen) .and. levels(gen) > 0 )

      case ( s_spectrometerType ) ! ...........  SPECTROMETERTYPE  .....
        call trace_begin ( me_spectrometertype, "MLSSignals.spectrometerType", &
          & son, cond=toggle(gen) .and. levels(gen) > 0 )
        spectrometerType%name = name
        deferred = .false.
        dacs = .false.
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
          case ( f_dacs )
            dacs = get_boolean(son)
          case ( f_first, f_last )
            ! Front end should be checking units now, but just in case....
            call expr_check ( gson, units, value, field, phyq_dimensionless )
            select case ( field )
            case ( f_first )
              first = value(1)
            case ( f_last )
              last = value(1)
            end select
          case ( f_start, f_step, f_width )
            ! Front end should be checking units now, but just in case....
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
        end do ! j = 2, nsons(key)
        if ( got(f_channels) ) then
          if ( any(got( (/ f_last, f_start, f_step, f_width /) )) ) &
            & call announceError ( badMix, f_channels, &
            & (/ f_last, f_start, f_step, f_width /) )
          if ( error == 0 ) then
            call allocate_Test ( spectrometerType%frequencies, &
              & nsons(channels)-2+first, 'spectrometerType%frequencies', &
              & moduleName, lowBound = first )
            call allocate_Test ( spectrometerType%widths, &
              & nsons(channels)-2+first, 'spectrometerType%widths', &
              & moduleName, lowBound = first )
            do k = 2, nsons(channels)
              call expr ( subtree(k,channels), units, value )
              spectrometerType%frequencies(k-2 +first) = value(1)
              spectrometerType%widths(k-2+first) = value(2)
            end do
            if ( any(units /= phyq_frequency) ) &
              ! Front end should be checking units now, but just in case....
              & call announceError ( wrongUnits, f_channels, &
                & (/ phyq_frequency /) )
          end if
        end if
        if ( got(f_start) ) then
          if ( .not. all( got((/ f_last, f_step, &
          & f_width /)) ) ) call announceError ( allOrNone, f_start, &
          & (/ f_last, f_step, f_width /) )
          if ( error == 0 ) then
            call allocate_Test ( spectrometerType%frequencies, last, &
              & 'spectrometerType%frequencies', moduleName, lowBound = first )
            call allocate_Test ( spectrometerType%widths, last, &
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
        spectrometerType%dacs = dacs
        if ( error == 0 ) call decorate ( key, addSpectrometerTypeToDatabase ( &
          & spectrometerTypes, spectrometerType ) )

        ! Nullify pointers to temporary stuff so it doesn't get hosed later
        nullify ( spectrometerType%frequencies )
        nullify ( spectrometerType%widths )
        call trace_end ( "MLSSignals.spectrometerType", &
          & cond=toggle(gen) .and. levels(gen) > 0 )

      case ( s_time ) ! ...................................  TIME  .....
        if ( timing ) then
          call sayTime
        else
          call time_now ( t1 )
          timing = .true.
        end if
      case default
        ! Shouldn't get here if the type checker worked
      end select

    end do

    if ( error > 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "Unable to create MLSSignals database" )

    if ( switchDetail(switches, 'sig') > -1 ) then
      call dump ( radiometers )
      call dump ( spectrometerTypes )
      call dump ( bands )
      call dump ( signals )
    end if
    call trace_end ( "MLSSignals", cond=toggle(gen) )
    if ( timing ) call sayTime

  contains
    ! --------------------------------------------  AnnounceError  -----
    subroutine AnnounceError ( Code, FieldIndex, MoreFields )
      integer, intent(in) :: Code       ! Index of error message
      integer, intent(in), optional :: FieldIndex ! f_...
      integer, intent(in), optional :: MoreFields(:)

      integer :: I

      error = max(error,1)
      call startErrorMessage ( son ) ! print '***** At <where>: '
      call output ( 'MLSSignals complained: ' )
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
      case ( priorTrouble )
        call output ( 'Unable to finish Signals because of prior trouble', &
          & advance='yes' )
      case ( wrongUnits )
        call output ( 'The values of the ' )
        call display_string ( field_indices(fieldIndex) )
        call output ( ' field have the wrong units -- ' )
        call display_string ( phyq_indices(moreFields(1)) )
        call output ( ' required.', advance='yes' )
      case ( unneededRadiometer )
        call output ( 'Radiometer field supplied when not necessary', advance='yes' )
      case ( noDeferredRadiometer )
        call output ( 'Radiometer field neeed for this signal as not specified in band', advance='yes' )
      case ( noModuleSpecified )
        call output ( 'No module found in either signal or radiometer definition', advance='yes' )
      case ( conflictingModules )
        call output ( 'Module given in signal and in radiometer conflict', advance='yes' )
      case ( notSpacecraft )
        call output ( 'Should only supply supportedModule for spacecraft modules', advance='yes' )
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
      call time_now ( t2 )
      call output ( "Timing for MLSSignals = " )
      call output ( dble(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime

  end subroutine MLSSignals

  ! ------------------------------------------  AddBandToDatabase  -----
  integer function AddBandToDatabase ( Database, Item )

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

    type(band_T), dimension(:), pointer :: Database
    type(band_T), intent(in) :: Item

    ! Local variables
    type (Band_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddBandToDatabase = newSize
  end function AddBandToDatabase

  ! ----------------------------------------  AddModuleToDatabase  -----
  integer function AddModuleToDatabase ( Database, Item )

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

    type(module_T), dimension(:), pointer :: Database
    type(module_T), intent(in) :: Item

    ! Local variables
    type (Module_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddModuleToDatabase = newSize
  end function AddModuleToDatabase

  ! ------------------------------------  AddRadiometerToDatabase  -----
  integer function AddRadiometerToDatabase ( Database, Item )

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

    type(radiometer_T), dimension(:), pointer :: Database
    type(radiometer_T), intent(in) :: Item

    ! Local variables
    type (Radiometer_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddRadiometerToDatabase = newSize
  end function AddRadiometerToDatabase

  ! ----------------------------------------  AddSignalToDatabase  -----
  integer function AddSignalToDatabase ( Database, Item )

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

    type(signal_T), dimension(:), pointer :: Database
    type(signal_T), intent(in) :: Item

    ! Local variables
    type (Signal_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddSignalToDatabase = newSize
  end function AddSignalToDatabase

  ! ------------------------------  AddSpectrometerTypeToDatabase  -----
  integer function AddSpectrometerTypeToDatabase ( Database, Item )

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

    type(spectrometerType_T), dimension(:), pointer :: Database
    type(spectrometerType_T), intent(in) :: Item

    ! Local variables
    type (spectrometerType_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddSpectrometerTypeToDatabase = newSize
  end function AddSpectrometerTypeToDatabase

  ! ----------------------------------------  DestroyBandDatabase  -----
  subroutine DestroyBandDatabase ( Bands )

    use Allocate_Deallocate, only: Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type(band_T), dimension(:), pointer :: Bands
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: S, Status
    if ( associated(bands) ) then
      s = size(bands) * storage_size(bands) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(bands(1)), addr)
      deallocate ( bands, stat = status )
      call test_deallocate ( status, moduleName, 'Band database', s, address=addr )
    end if

  end subroutine DestroyBandDatabase

  ! --------------------------------------  DestroyModuleDatabase  -----
  subroutine DestroyModuleDatabase ( Modules )

    use Allocate_Deallocate, only: Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type(module_T), dimension(:), pointer :: Modules
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: S, Status
    if ( associated(modules) ) then
      s = size(modules) * storage_size(modules) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(modules(1)), addr)
      deallocate ( modules, stat = status )
      call test_deallocate ( status, moduleName, 'Modules database', s, address=addr )
    end if
  end subroutine DestroyModuleDatabase

  ! ----------------------------------  DestroyRadiometerDatabase  -----
  subroutine DestroyRadiometerDatabase ( Radiometers )

    use Allocate_Deallocate, only: Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type(radiometer_T), dimension(:), pointer :: Radiometers
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: S, Status
    if ( associated(radiometers) ) then
      s = size(radiometers) * storage_size(radiometers) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(radiometers(1)), addr)
      deallocate ( radiometers, stat = status )
      call test_deallocate ( status, moduleName, 'Radiometer database', s, address=addr )
    end if
  end subroutine DestroyRadiometerDatabase

  ! ----------------------------------------------  DestroySignal  -----
  subroutine DestroySignal ( Signal, justChannels )
    ! Destroy one signal in the signals database (or elsewhere)
    ! If justChannels is set, don't destroy frequencies/widths, even
    ! for deferred, as it is a copy of the master.
    type(signal_T), intent(inout) :: SIGNAL
    logical, intent(in), optional :: JUSTCHANNELS

    logical :: MYJUSTCHANNELS

    myJustChannels = .false.
    if ( present(justChannels) ) myJustChannels = justChannels
    call deallocate_test ( signal%channels, 'Signal%channels', moduleName )
    ! Don't destroy Frequencies or Widths unless signal%Deferred.  Those
    ! fields are shallow copies here.  They're destroyed in
    ! DestroySpectrometerType.
    if ( signal%deferred .and. .not. myJustChannels ) then
      call deallocate_test ( signal%frequencies, "signal%frequencies", &
        & moduleName )
      call deallocate_test ( signal%widths, "signal%widths", moduleName )
    end if
  end subroutine DestroySignal

  ! --------------------------------------  DestroySignalDatabase  -----
  subroutine DestroySignalDatabase ( Signals, justChannels )

    use Allocate_Deallocate, only: Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type(signal_T), dimension(:), pointer :: Signals
    logical, intent(in), optional :: JUSTCHANNELS
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: I, S, Status

    ! Executable code
    if ( associated(signals) ) then
      do i = 1, size(signals)
        call destroySignal ( signals(i), justChannels )
      end do
      s = size(signals) * storage_size(signals) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(signals(1)), addr)
      deallocate ( signals, stat = status )
      call test_deallocate ( status, moduleName, 'Signal database', s, address=addr )
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

    use Allocate_Deallocate, only: Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type(spectrometerType_T), dimension(:), pointer :: spectrometerTypes
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: I, S, Status
    if ( associated(spectrometerTypes) ) then
      do i = 1, size(spectrometerTypes)
        call destroySpectrometerType ( spectrometerTypes(i) )
      end do
      s = size(spectrometerTypes) * storage_size(spectrometerTypes) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(spectrometerTypes(1)), addr)
      deallocate ( spectrometerTypes, stat = status )
      call test_deallocate ( status, moduleName, 'Spectrometer database', s, &
        & address=addr )
    end if
  end subroutine DestroySpectrometerTypeDatabase

  ! ------------------------------------------  DisplayRadiometer  -----
  subroutine DisplayRadiometer ( Index )
    integer, intent(in) :: Index ! in radiometer database
    call display_string ( radiometers(index)%prefix )
    call output ( ':' )
    call display_string ( radiometers(index)%suffix, strip=.true. )
  end subroutine

  ! ------------------------------------  DisplaySignalName_index  -----
  subroutine DisplaySignalName_index ( Signal, Advance, Before, Sideband, &
                                     & Channel, OtherChannels )
    ! Given a signal object, this routine displays a full signal name.
    integer, intent(in) :: SIGNAL
    character(len=*), intent(in), optional :: Advance, Before
    integer, intent(in), optional :: Sideband ! Use this instead of Signal's
    integer, intent(in), optional :: Channel  ! Use this instead of Signal's
    logical, pointer, optional :: OtherChannels(:) ! Use these instead of Signal's
    call displaySignalName ( signals(signal), advance, before, sideband, channel, &
                           & otherChannels )
   end subroutine DisplaySignalName_index

  ! -----------------------------------  DisplaySignalName_signal  -----
  subroutine DisplaySignalName_signal ( Signal, Advance, Before, Sideband, &
                                      & Channel, OtherChannels )
    ! Given a signal object, this routine displays a full signal name.
    use String_Table, only: Display_String
    type(signal_T), intent(in) :: SIGNAL
    character(len=*), intent(in), optional :: Advance, Before
    integer, intent(in), optional :: Sideband ! Use this instead of Signal's
    integer, intent(in), optional :: Channel  ! Use this instead of Signal's
    logical, pointer, optional :: OtherChannels(:) ! Use these instead of Signal's
    character(len=15) :: BandName
    logical, pointer :: Channels(:)
    logical :: First
    integer :: I, J, SB

    if ( present(before) ) call output ( before )
    if ( present(sideband) ) then
      sb = sideband
    else
      sb = signal%sideband
    end if
    call displayRadiometer ( signal%radiometer )
    call output ( '.' )
    call GetBandName ( signal%band, bandName, sideband=sb )
    call output ( trim(bandName) )
    call output ( signal%switch, before='.S', after='.' )
    call display_string ( spectrometerTypes(signal%spectrometerType)%name )
    call output ( signal%spectrometer, before='-' )
    if ( present(channel) ) then
      call output ( channel, before='.C' )
    else
      channels => signal%channels
      if ( present(otherChannels) ) channels => otherChannels
      if ( associated(channels) ) then
        if ( .not. all(channels) .or. &
          & lbound(channels,1) /= lbound(signal%frequencies,1) .or. &
          & ubound(channels,1) /= ubound(signal%frequencies,1) ) then
          first = .true.
          call output ( '.C' )
          i = lbound(channels, 1)
oc:       do
            do
              if ( i > ubound(channels, 1) ) exit oc
              if ( channels(i) ) exit
              i = i + 1
            end do
            if ( .not. first ) call output ( '+' )
            first = .false.
            j = i
            do while ( j < ubound(channels, 1) )
              if ( .not. channels(j+1) ) exit
              j = j + 1
            end do
            if ( j > i ) then
              call output ( i, after = ':' )
              call output ( j )
            else
              call output ( i )
            end if
            i = j + 1
          end do oc
        end if
      end if
    end if
    call output ( '', advance=advance )
  end subroutine DisplaySignalName_signal

  ! -------------------------------------------------  Dump_Bands  -----
  subroutine DUMP_BANDS ( BANDS )
    type (Band_T), intent(in) :: BANDS(:)
    integer :: i
    call output ( 'BANDS: SIZE = ')
    call output ( size(bands), advance='yes' )
    do i = 1, size(bands)
      call output ( i )
      call output ( ': ' )
      call display_string (bands(i)%prefix)
      call output ( ':' )
      call display_string (bands(i)%suffix, strip=.true. )
      call output ( '   Radiometer: ')
      if ( bands(i)%radiometer /= 0 ) then
        call output ( bands(i)%radiometer )
        call output ( ' - ' )
        call display_string ( radiometers(bands(i)%radiometer)%prefix )
      else
        call output ( ' deferred' )
      end if
      call output ( '   SpectrometerType: ' )
      call output ( bands(i)%spectrometerType )
      call output ( ' - ' )
      call display_string ( spectrometerTypes(bands(i)%spectrometerType)%name, &
        & advance='yes' )
      call output ( '   Frequency: ')
      call output ( bands(i)%centerFrequency, advance='yes' )
    end do
  end subroutine DUMP_BANDS

  ! -------------------------------------------  Dump_All  -----
  ! Dump all databases related to Signals, Modules, etc.
  subroutine Dump_All
    !
    call headline( 'All databases related to Signals', &
      & fillChar='-', before='*', after='*' )
    call output ( '   Instrument: ', advance='no' )
    call display_string ( lit_indices(instrument) )
    call newLine
    call dump_modules
    call dump_radiometers ( radiometers )
    call dump_bands( bands )
    call dump_spectrometerTypes( spectrometerTypes )
    call Dump_Signals( signals )
  end subroutine Dump_All

  ! -------------------------------------------  Dump_Modules  -----
  subroutine Dump_Modules
    !
    character(len=16), dimension(:,:), pointer :: array
    integer :: i
    integer :: n
    ! Executable
    if ( .not. associated(modules) ) then
      call output( '(modules not associated)', advance='yes' )
      return
    endif
    call headline( 'modules' )
    nullify( array )
    n = size(modules)
    allocate( array(n+1, 4 ) )
    array(1,1) = 'name'
    array(1,2) = 's/c'
    array(1,3) = 'Aura'
    array(1,4) = 'Supports'
    array(2:n+1,1) = modules(:)%nameString
    do i=1, n
      if ( len_trim(modules(i)%nameString) < 1 ) &
        & call get_string( modules(i)%name, array(i+1,1) )
      if ( modules(i)%supportedModule == 0 ) then
        array(i+1,4) = '<none>'
      else
        call get_string ( modules(modules(i)%supportedModule)%name,&
          & array(i+1,4) )
      end if
    enddo
    array(2:n+1,2) = merge( 'T', 'F', modules(:)%spaceCraft )
    array(2:n+1,3) = merge( 'T', 'F', modules(:)%Aura )
    call outputTable( array, border='-', headliner='-' )
    deallocate( array )
  end subroutine Dump_Modules

  ! -------------------------------------------  Dump_OneRadiometer  -----
  subroutine Dump_OneRadiometer ( RADIOMETER )
    type (Radiometer_T), intent(in) :: RADIOMETER
    call display_string (radiometer%prefix)
    call output ( ':' )
    call display_string (radiometer%suffix, advance='yes', strip=.true. )
    call output ( '   Polarization: ' )
    call display_string (lit_indices(radiometer%polarization))
    call output ( '   Module: ')
    if ( radiometer%instrumentModule /= 0 ) then
      call output ( radiometer%instrumentModule )
      call output ( ' - ' )
      call display_string ( modules(radiometer%instrumentModule)%name, advance='yes' ) 
    else
      call output ( '<undefined>', advance='yes' )
    end if
    call output ( '   LO: ')
    call output ( radiometer%lo )
    call output ( '   Single sideband: ' )
    call output ( radiometer%singleSideband, advance='yes' )
  end subroutine Dump_OneRadiometer

  ! -------------------------------------------  Dump_Radiometers  -----
  subroutine DUMP_RADIOMETERS ( RADIOMETERS )
    type (Radiometer_T), intent(in) :: RADIOMETERS(:)
    integer :: i
    call output ( 'RADIOMETERS: SIZE = ')
    call output ( size(radiometers), advance='yes' )
    do i = 1, size(radiometers)
      call output ( i,1 )
      call output ( ': ')
      call dump_OneRadiometer( radiometers(i) )
    end do
  end subroutine DUMP_RADIOMETERS

  ! ------------------------------------------------  Dump_Signal  -----
  subroutine Dump_Signal ( Signal, Details, OtherChannels )
    type (signal_T), intent(in) :: SIGNAL
    integer, intent(in), optional :: Details ! 0  => don't dump frequencies
                                             ! -1 => dump only name
    logical, pointer, optional :: OtherChannels(:)
    logical, pointer :: Channels(:)
    integer :: My_Details
    character (len=80) :: Str
    my_details = 0
    if ( present(details) ) my_details = details
    channels => signal%channels
    if ( present(otherChannels) ) channels => otherChannels
    if ( signal%name > 0 ) call display_string ( signal%name, advance='yes' )
    call output ( '   Signal: ' )
    call displaySignalName ( signal, otherChannels=otherChannels )
    if ( My_Details < 0 ) then
      call NewLine 
      return
    endif
    call output ( '   Module: ')
    call output ( signal%instrumentModule )
    call output ( ' - ' )
    if ( associated(modules) ) then
      call display_string ( modules(signal%instrumentModule)%name, &
        & advance='yes' )
    else
      call output ( 'Cannot get module name', advance='yes' )
    end if
    call output ( '   Radiometer: ')
    call output ( signal%radiometer )
    call output ( ' - ' )
    if ( associated(radiometers) ) then
      call getRadiometerName ( signal%radiometer, str )
      call output ( TRIM(str) )
    else
      call output ( 'Cannot get radiometer name', advance='yes' )
    end if
    call output ( '   First LO: ')
    call output ( signal%lo, after=' MHz', advance='yes' )
    call output ( '   Band: ')
    call output ( signal%band )
    call output (' - ')
    if ( associated(bands) ) then
      call getBandName ( signal%band, str )
      call output ( TRIM(str) )
    else
      call output ( 'Cannot get band name', advance='yes' )
    end if
    call output ( '   Band center frequency: ')
    call output ( signal%centerFrequency, after=' MHz', advance='yes' )
    call output ( '   SpectrometerType: ')
    call output ( signal%spectrometerType )
    call output ( ' - ' )
    if ( associated(spectrometerTypes) ) then
      call display_string ( spectrometerTypes(signal%spectrometerType)%name )
    else
      call output ( 'Cannot get spectrometer type name', advance='yes' )
    end if
    call output ( '   DACS?: ' )
    call output ( signal%dacs )
    call output ( '   Channels: ' )
    call output ( lbound(signal%frequencies,1), 3 )
    call output ( ':' )
    call output ( ubound(signal%frequencies,1), 3, advance='yes' )
    call output ( '   Sideband: ' )
    call output ( signal%sideband )
    call output ( '   Single Sideband: ' )
    call output ( signal%singleSideband, advance='yes')
    if ( my_details > 0) then
      call output ( '   Frequencies (MHz)' )
      if ( signal%deferred ) call output ( ' (deferred)' )
      call output ( ':', advance='yes' )
      call dump ( signal%frequencies )
      call output ( '   Widths' )
      if ( signal%deferred ) call output ( '(deferred)' )
      call output ( ':', advance='yes' )
      call dump ( signal%widths )
    else
      call output ( '   Frequencies and widths are' )
      if ( .not. signal%deferred ) call output ( ' not' )
      call output ( ' deferred', advance='yes' )
    end if ! my_details
    if (associated(channels)) then
      call output ( '   Channel Flags:', advance='yes' )
      call dump ( channels, lbound=lbound(channels,1) )
    else
      call output ( '   All channels selected', advance='yes' )
    end if
  end subroutine Dump_Signal

  ! -----------------------------------------------  Dump_Signals  -----
  subroutine Dump_Signals ( Signals, isSignalCritical, Details, OtherChannels )
    type (signal_T), intent(in)                 :: Signals(:)
    logical, dimension(:), intent(in), optional :: isSignalCritical
    integer, intent(in), optional               :: Details ! How verbose?
    logical, pointer, optional                  :: OtherChannels(:)
    integer :: I
!   character(len=*), parameter                 :: The80 = &
!     & '12345678901234567890123456789012345678901234567890123456789012345678901234567890'
!   character(len=*), parameter                 :: TheDecades = &
!     & '         1         2         3         4         5         6         7         8'
    character(len=2) :: str
    integer :: myDetails
    ! Executable
    myDetails = 0
    if ( present(Details) ) myDetails = Details
    call outputNamedValue ( 'Size(signals)', size(signals), options='--Headline' )
    if ( myDetails < 0 ) then
      do i=1, size(signals)
        call output( i, advance='no' )
        call blanks ( 2 )
        call Dump_Signal ( signals(i), Details )
      enddo
      return
    endif
    ! call output( The80, advance='yes' )
    ! call output( TheDecades, advance='yes' )
    ! Column headers: 4 columns
    ! Index   signal  band  critical module
    call output( 'Index', advance='no' )
    call tab ( 3 )
    call output( 'signal', advance='no' )
    call tab ( 7 )
    call output( 'band', advance='no' )
    call tab ( 9 )
    call output( 'critical', advance='no' )
    call tab ( 11 )
    call output( 'module', advance='no' )
    if ( myDetails > 0 ) then
      call tab ( 13 )
      call output( 'DACS', advance='no' )
      call tab ( 15 )
      call output( 'Sideband', advance='no' )
      call tab ( 17 )
      call output( 'Spectrometer', advance='no' )
      call tab ( 19 )
      call output( 'Channels', advance='no' )
    endif
    call NewLine
    do i = 1, size(signals)
      call output ( i )
      call tab
      call displaySignalName ( signals(i), otherChannels=otherChannels )
      call tab ( 7 )
      call output ( signals(i)%Band )
      if ( present( isSignalCritical ) ) then
        str = merge( 'y', 'n', isSignalCritical(i) )
      else
        str = '?'
      endif
      call tab ( 9 )
      call output ( str )
      call tab ( 11 )
      call display_string ( modules(signals(i)%instrumentModule)%name, &
        & advance='no' )
      if ( myDetails > 0 ) then
        call tab ( 13 )
        call output ( signals(i)%dacs )
        call tab ( 15 )
        call output ( signals(i)%sideband )
        call tab ( 17 )
        call display_string ( spectrometerTypes(signals(i)%spectrometerType)%name )
        call tab ( 19 )
        call output( lbound(signals(i)%frequencies,1) )
        call output( ':' )
        call output( ubound(signals(i)%frequencies,1) )
      endif
      call NewLine
    end do
  end subroutine Dump_Signals

  ! --------------------------------------  Dump_SpectrometerType  -----
  subroutine DUMP_SPECTROMETERTYPE ( SPECTROMETERTYPE, N )
    type (SpectrometerType_T), intent(in) :: SPECTROMETERTYPE
    integer, intent(in), optional :: N

    if ( present(n) ) then
      call output ( n, 1 )
      call output ( ': ')
    end if
    call display_string ( spectrometerType%name, advance='yes' )
    if ( associated(spectrometerType%frequencies) ) then
      call output ( '  Channels: ' )
      call output ( lbound(spectrometerType%frequencies,1), 3 )
      call output ( ':' )
      call output ( ubound(spectrometerType%frequencies,1), 3 )
      call output ( '  DACS?: ' )
      call output ( spectrometerType%dacs, advance='yes' )
      call output ( '  Frequencies:', advance='yes' )
      call dump ( spectrometerType%frequencies )
      call output ( '  Widths:', advance='yes' )
      call dump ( spectrometerType%widths )
    else
      call output ('   Frequencies and widths deferred.' )
      call output ( '  DACS?: ' )
      call output ( spectrometerType%dacs, advance='yes' )
    end if

  end subroutine DUMP_SPECTROMETERTYPE

  ! -------------------------------------  Dump_SpectrometerTypes  -----
  subroutine DUMP_SPECTROMETERTYPES ( SPECTROMETERTYPES )
    type (SpectrometerType_T), intent(in) :: SPECTROMETERTYPES(:)
    integer :: i
    call output ( 'SPECTROMETERTYPES: SIZE = ')
    call output ( size(spectrometerTypes), advance='yes' )
    do i = 1, size(spectrometerTypes)
      call dump ( spectrometerTypes(i), i )
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
    integer :: ENDOFFIRSTNUMBER
    logical :: FOUNDEND
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
    foundEnd = .false.
    do endOfFirstNumber = 2, len_trim ( string_text )
      if ( index ( '0123456789', string_text(endOfFirstNumber:endOfFirstNumber) ) == 0 ) &
        & foundEnd = .true.
      if ( foundEnd ) exit
    end do
    if ( .not. foundEnd ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, 'Cannot understand band name' )
    string_text = string_text(1:endOfFirstNumber-1) // TRIM(sb_char) // &
      & string_text(endOfFirstNumber:LEN_TRIM(string_text))
    if ( (.not. my_noSuffix) .and. &
      &  (len_trim(string_text) < len(string_text)) ) then
      string_text = TRIM(string_text) // ':'
      call get_string ( bands(band)%suffix,&
        & string_text(LEN_TRIM(string_text)+1:), cap=.true., strip=.true. )
    end if

  end subroutine GetBandName

  ! ------------------------------------  GetFirstChannel --------------
  integer function GetFirstChannel ( signal )
    integer, intent(in) :: SIGNAL
    GetFirstchannel = lbound ( signals(signal)%frequencies, 1 )
  end function GetFirstChannel

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

  ! ----------------------------------------------  GetModuleIndex  -----
  subroutine GetModuleIndex(string_text, instrumentModule)
    ! Returns the index in the module database, given module name in mixed case
    ! Returns 0 if module name not found
    ! (inverse function: GetModuleName)
    integer, intent(out) :: instrumentModule
    character (len=*), intent(in) :: string_text
    ! Local variables
    character (len=len(string_text))             :: string_test
    if ( size(modules) < 1 ) then
      instrumentModule = 0
      return
    end if
    do instrumentModule=1, size(modules)
      if ( modules(instrumentModule)%name > 0 ) then
        call Get_String ( modules(instrumentModule)%name, string_test )
        if ( LowerCase(trim(string_text)) == LowerCase(trim(string_test))) &
          & return
      elseif ( len_trim(modules(instrumentModule)%nameString) > 0 ) then
        string_test = modules(instrumentModule)%nameString
        if ( LowerCase(trim(string_text)) == LowerCase(trim(string_test))) &
          & return
      end if
    end do
    instrumentModule = 0
  end subroutine GetModuleIndex

  ! ----------------------------------------------  GetModuleName  -----
  subroutine GetModuleName(instrumentModule, string_text)
    ! Given the index in the module database, returns module name in mixed case
    ! (inverse function: GetModuleIndex)
    integer, intent(in) :: instrumentModule
    character (len=*), intent(out) :: string_text
    if ( modules(instrumentModule)%name < 1 .or. &
      & len_trim(modules(instrumentModule)%nameString) > 0 ) then
      string_text =  modules(instrumentModule)%nameString
    else
      call Get_String ( modules(instrumentModule)%name, string_text )
    endif
    ! call outputNamedValue ( 'got module name', trim(string_text) )
  end subroutine GetModuleName

  ! ----------------------------------------------  GetSignalIndex  -----
  subroutine GetSignalIndex(string_text, signal_index)
    ! Returns the index in the signals database, given signal name in mixed case
    ! Returns 0 if signal name not found
    ! (inverse function: GetSignalName)
    integer, intent(out) :: signal_index
    character (len=*), intent(in) :: string_text
    ! Local variables
    character (len=len(string_text))             :: string_test
    if ( size(signals) < 1 ) then
      signal_index = 0
      return
    end if
    do signal_index=1, size(signals)
      if ( signals(signal_index)%name > 0 .and. .false. ) then
        call Get_String ( signals(signal_index)%name, string_test )
        if ( LowerCase(trim(string_text)) == LowerCase(trim(string_test))) &
          & return
      else
        call GetNameOfSignal ( signals(signal_index), string_test )
        if ( LowerCase(trim(string_text)) == LowerCase(trim(string_test))) &
          & return
      end if
    end do
    signal_index = 0
  end subroutine GetSignalIndex

  ! --------------------------------------------  GetNameOfSignal  -----
  subroutine GetNameOfSignal ( Signal, String_text, NoRadiometer, NoBand, &
    & NoSwitch, NoSpectrometer, NoChannels, NoSuffix, sideband, channel, &
    & OtherChannels )
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
    integer, intent(in), optional :: CHANNEL ! Only this channel, noChannels overrides it
    logical, pointer, optional :: OtherChannels(:) ! instead of from Signal,
                                                   ! Channel overrides this

    ! Local variables
    logical :: First     ! First channel in signal text
    logical, pointer :: Channels(:) ! From signal or OtherChannels
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

    if ( .not. my_noChannels ) then
      channels => signal%channels
      if ( present(otherChannels) ) channels => otherChannels
      if ( present(channel) ) then
        l = len_trim(string_text)
        call addToSignalString ( '.C' )
        write ( word,'(I8)') channel
        call addToSignalString ( word )
      else if ( associated(channels) ) then
        if ( .not. all(channels) .or. &
          & lbound(channels,1) /= lbound(signal%frequencies,1) .or. &
          & ubound(channels,1) /= ubound(signal%frequencies,1) ) then
          l = len_trim(string_text)
          call addToSignalString ( '.C' )
          i = lbound(channels, 1)
          first = .true.
          oc: do
            do
              if ( i > ubound(channels, 1) ) exit oc
              if ( channels(i) ) exit
              i = i + 1
            end do
            if ( .not. first ) call addToSignalString ( '+' )
            first = .false.
            j = i
            do while ( j < ubound(channels, 1) )
              if ( .not. channels(j+1) ) exit
              j = j + 1
            end do
            if ( j > i ) then
              write ( word,'(I8)') i
              call addToSignalString ( word )
              call addToSignalString ( ':' )
              write ( word,'(I8)') j
              call addToSignalString ( word )
            else
              write ( word,'(I8)') i
              call addToSignalString ( word )
            end if
            i = j + 1
          end do oc
        end if
      end if
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

  ! ----------------------------------------------  GetRadiometerIndex  -----
  subroutine GetRadiometerIndex(string_text, radiometer)
    ! Returns the index in the radiometer database, 
    ! given radiometer name in mixed case
    ! Returns 0 if radiometer name not found
    ! Examples of string_text:
    ! 'R1A:118', 'R3:240', 'R3',  (yes, we may omit suffix or prefix)
    ! Note: beware of an ambiguous string_tetxt;
    ! e.g., if you just use a suffix, e.g., '118', you will end up with
    ! the first match, 'R1A:118' and not 'R1B:118'
    ! (inverse function: GetradiometerName)
    integer, intent(out) :: radiometer
    character (len=*), intent(in)        :: string_text
    ! Local variables
    character (len=MAXRADIOMETERNAMELEN) :: prefix
    character (len=MAXRADIOMETERNAMELEN) :: suffix
    ! Executable
    radiometer = 0 ! If no matching radiometer found
    if ( size(radiometers) < 1 ) then
      return
    end if
    do radiometer=1, size(radiometers)
      ! Return first match
      if ( radiometers(radiometer)%prefix > 0 ) then
        call Get_String ( radiometers(radiometer)%prefix, prefix )
        ! Did we just try the prefix?
        if ( LowerCase(trim(string_text)) == LowerCase(trim(prefix))) &
          & return
        call get_string ( radiometers(radiometer)%suffix, &
          & suffix, cap=.true., strip=.true. )
        ! Did we just try the suffix?
        if ( LowerCase(trim(string_text)) == LowerCase(trim(suffix))) &
          & return
        ! Did we go for prefix:suffix?
        if ( LowerCase( trim(string_text) ) == &
          & LowerCase( trim(prefix) // ':' // trim(suffix) ) ) &
          & return
      end if
    end do
    ! Oops, sorry no match found
  end subroutine GetRadiometerIndex

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
    end if

  end subroutine GetRadiometerName

  ! --------------------------------------------  GetSidebandLoop  -----
  subroutine GetSidebandLoop ( signal, sideband, split, &
    & sidebandStart, sidebandStop, sidebandStep )
    ! This routine gets the loop limits for a loop over sidebands from
    ! a signal in the database.
    integer, intent(in) :: SIGNAL       ! Index into signals
    integer, intent(in) :: SIDEBAND     ! -1,0,1
    logical, intent(in) :: SPLIT        ! If set do a split sideband loop for folded
    integer, intent(out) :: SIDEBANDSTART ! Loop lower limit
    integer, intent(out) :: SIDEBANDSTOP ! Loop upper limit
    integer, intent(out) :: SIDEBANDSTEP ! Loop step

    ! Executable code
    if ( sideband == 0 ) then
      if ( split ) then
        if ( signals(signal)%singleSideband /= 0 ) then
          sidebandStart = signals(signal)%singleSideband
          sidebandStop = sidebandStart
          sidebandStep = 1
        else
          sidebandStart = -1
          sidebandStop = 1
          sidebandStep = 2
        end if
      else
        sidebandStart = 0
        sidebandStop = sidebandStart
        sidebandStep = 1
      end if
    else
      sidebandStart = sideband
      sidebandStop = sidebandStart
      sidebandStep = 1
    end if
  end subroutine GetSidebandLoop

  ! ---------------------------------------  GetSidebandStartStop  -----
  subroutine GetSidebandStartStop ( Signal, SidebandStart, SidebandStop )
    ! This routine also gets the sideband start and stop, but from any
    ! signal, including one parsed from a config.  We don't bother to
    ! compute the step, because if SidebandStart == SidebandStop it doesn't
    ! matter what the step is, and otherwise the step should be 2. So users
    ! should always just use 2.
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    type(signal_t), intent(in) :: Signal
    integer, intent(out) :: SidebandStart, SidebandStop
    if ( ( signal%sideband == 0 ) .and. &
      &  ( signal%singleSideband == 0 ) ) then
      ! Do a folded measurement
      sidebandStart = -1
      sidebandStop = 1
    else
      ! It's either a single sideband radiometer, or the user requested a
      ! specific sideband.
      ! Check sanity, if they are both non zero they should be the same.
      if ( ( signal%singleSideband /= 0 ) .and. &
        &  ( signal%sideband /= 0 ) .and. &
        &  ( signal%singleSideband /= &
        &    signal%sideband ) ) call MLSMessage ( &
        &      MLSMSG_Error, ModuleName, &
        &      "User requested a sideband that doesn't exist" )
      ! OK, use whichever one is given
      if ( signal%singleSideband /= 0 ) then
        sidebandStart = signal%singleSideband
      else
        sidebandStart = signal%sideband
      end if
      sidebandStop = sidebandStart
    end if
  end subroutine GetSidebandStartStop

  ! --------------------------------------------------  GetSignal  -----
  type (Signal_T) function GetSignal(signal)
    ! Given the database index, this routine returns the signal data structure
    integer, intent(in) :: SIGNAL       ! Requested signal
    
    GetSignal = signals(signal)
  end function GetSignal    

  ! ----------------------------------------------  GetSignalName  -----
  subroutine GetSignalName ( Signal, String_text, NoRadiometer, NoBand, &
    & NoSwitch, NoSpectrometer, NoChannels, NoSuffix, Sideband, Channel, &
    & OtherChannels )
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
    integer, intent(in), optional :: CHANNEL ! Only this channel, noChannels overrides this
    logical, pointer, optional :: OtherChannels(:) ! instead of from Signal,
                                                   ! Channel overrides this

    call getNameOfSignal ( signals(signal), string_text, noRadiometer, noBand, &
      & noSwitch, noSpectrometer, noChannels, noSuffix, sideband, channel, &
      & otherChannels )
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
    word = adjustl(Capitalize(word))
    if ( len_trim(string_text)+len_trim(word) < len(string_text) ) &
      & string_text = TRIM(string_text) // TRIM(word)
    string_text = Capitalize(string_text)
  end subroutine GetSpectrometerTypeName

  ! -----------------------------------------  IsAnyModuleSpacecraft  -----
  logical function IsAnyModuleSpacecraft()
    ! Returns true if any module is really the spacecraft
    integer :: i
    IsAnyModuleSpacecraft = modules(1)%spacecraft
    do i = 2, size(modules)
      IsAnyModuleSpacecraft = IsAnyModuleSpacecraft &
        & .or. modules(i)%spacecraft
    enddo
  end function IsAnyModuleSpacecraft

  ! -----------------------------------------  IsModuleSpacecraft  -----
  logical function IsModuleSpacecraft(thisModule)
    ! Returns true if the module is really the spacecraft
    integer, intent(in) :: thisModule
    IsModuleSpacecraft = modules(thisModule)%spacecraft
  end function IsModuleSpacecraft

  ! -----------------------------------------  IsSpaceCraftAura  -----
  logical function IsSpaceCraftAura()
    ! Returns true if the s/c is really the Aura spacecraft
    IsSpaceCraftAura = associated( modules )
    if ( IsSpaceCraftAura ) IsSpaceCraftAura = all( modules%Aura )
  end function IsSpaceCraftAura

  ! -----------------------------------------------  MatchSignals  -----
  integer function MatchSignals ( Signals, Probe, sideband, channel, matchFlags, &
    & NoMatchFails, FromWhere, DSBSSB )
    ! Given an array Signals, find the one in the array that provides
    ! the smallest superset of features of the signal Probe.  The result
    ! is zero if no signals match.
    ! If sideband or channel are present they are used for the probe, instead
    ! of the values in probe.

    type(signal_T), dimension(:), intent(in) :: Signals
    type(signal_T), intent(in) :: Probe
    integer, intent(in), optional :: sideband     ! Use this instead of probe%sideband
    logical, dimension(size(signals)), intent(out), optional :: matchFlags
    integer, intent(in), optional :: CHANNEL      ! Just this channel
    logical, intent(in), optional :: NoMatchFails ! Fail if no match
    character(len=*), intent(in), optional :: FromWhere ! For an error message
    logical, intent(in), optional :: DSBSSB       ! OK if one sig is DSB, other is SSB

    integer :: BestMatch                ! The smallest number of 
    integer :: I                        ! Loop inductors, subscripts
    integer :: NumChannelsMatch

    if ( present(matchFlags) ) matchFlags = .false.

    bestMatch = huge(bestMatch)
    matchSignals = 0
    do i = 1, size(signals)
      numChannelsMatch = matchSignalPair ( signals(i), probe, sideband, &
        & channel, DSBSSB=DSBSSB )
      if ( numChannelsMatch > 0 ) then
        if ( present( matchFlags ) ) matchFlags(i) = .true.
        if ( numChannelsMatch < bestMatch ) then
          matchSignals = i
          bestMatch = numChannelsMatch
        end if
      end if
    end do
    if ( present(noMatchFails) ) then
      if ( noMatchFails .and. matchSignals == 0 ) then
        if ( present(fromWhere) ) then
          call MLSMessage ( MLSMSG_Error, moduleName, &
            & 'No match for requested signal from ' // trim(fromWhere) )
        else
          call MLSMessage ( MLSMSG_Error, moduleName, &
            & 'No match for requested signal' )
        end if
      end if
    end if
  end function MatchSignals

  ! --------------------------------------------  MatchSignalPair  -----
  integer function MatchSignalPair ( Signal, Probe, sideband, channel, &
    & DSBSSB, NoMatchFails, FromWhere )
    ! Given a Signal and a Probe, determine whether they match.
    ! If Sideband or Channel are present they are used for the probe, instead
    ! of the values in Probe.
    ! The result is -1 if not, else the number of matching channels.

    type(signal_T), intent(in) :: Signal, Probe
    integer, intent(in), optional :: sideband     ! Use this instead of probe%sideband
    integer, intent(in), optional :: CHANNEL      ! Just this channel
    logical, intent(in), optional :: DSBSSB       ! OK if one DSB, one SSB
    logical, intent(in), optional :: NoMatchFails ! Fail if no match
    character(len=*), intent(in), optional :: FromWhere ! For an error message

    logical :: Match
    logical :: MyDSBSSB
    integer :: MySideband               ! Either sideband or probe%sideband
    logical :: SidebandOK

    matchSignalPair = -1

    mySideband = probe%sideband
    if ( present(sideband) ) mySideband = sideband

    myDSBSSB = .false.
    if ( present(DSBSSB) ) myDSBSSB = DSBSSB

    ! First, the signal must have the same band, instrument module,
    ! radiometer, spectrometer, spectrometer type and switch number as
    ! the probe signal
    sidebandOK = (signal%sideband == mySideband) .or. &
      &          myDSBSSB .and. (signal%sideband * mySideband == 0)
    match = signal%band == probe%band .and. &
      &  signal%instrumentModule == probe%instrumentModule.and. &
      &  signal%radiometer == probe%radiometer .and. &
      &  sidebandOK .and. &
      &  signal%spectrometer == probe%spectrometer .and. &
      &  signal%spectrometerType == probe%spectrometerType .and. &
      &  signal%switch == probe%switch
    if ( match ) then
      ! Now the channels in Probe all have to be present in signal
      if (present(channel)) then        ! User asked for a specific channel
        if ( associated(signal%channels) ) then
          if ( channel < lbound(signal%channels,1) ) then
            match = .false.
          else if ( channel > ubound(signal%channels,1) ) then
            match = .false.
          else
            match = signal%channels(channel)
          end if
        else
          match = .true.
        end if
      else
        match = (.not. associated(probe%channels)) .or. &
          & (.not. associated(signal%channels) )
        if ( .not. match ) match = all( (probe%channels .and. &
          & signal%channels(lbound(probe%channels,1):ubound(probe%channels,1)) ) &
          & .eqv. probe%channels )
      end if
      if ( match ) then
        if ( associated(signal%channels) ) then
          matchSignalPair = count(signal%channels)
        else
          matchSignalPair = size( signal%frequencies )
        end if
      end if
    end if
    if ( .not. match .and. present(noMatchFails) ) then
      if ( noMatchFails ) then
        if ( present(fromWhere) ) then
          call MLSMessage ( MLSMSG_Error, moduleName, &
            & 'No match for requested signal from ' // trim(fromWhere) )
        else
          call MLSMessage ( MLSMSG_Error, moduleName, &
            & 'No match for requested signal' )
        end if
      end if
    end if
  end function MatchSignalPair

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
      if ( present(channel) ) then
        noPbChannels = 1
      else
        if ( associated(probe(i)%channels ) ) then
          noPbChannels = noPbChannels + count( probe(i)%channels )
        else
          noPbChannels = noPbChannels + size ( probe(i)%frequencies )
        end if
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

  ! ------------------------------------------------ PVMPackSignal ---
  subroutine PVMPackSignal ( signal )
    use PVMIDL, only: PVMIDLpack
    ! Dummy arguments
    type ( Signal_T ), intent(in) :: SIGNAL

    ! Local variables
    integer :: INFO                     ! Flag from PVM

    ! Executable code
    ! Let's pack things by type do scalars first
    call PVMIDLPack ( (/ signal%centerFrequency, signal%LO /), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Packing signal lo/cf' )
    call PVMIDLPack ( (/ signal%band, signal%direction, signal%index, &
      & signal%instrumentModule, signal%name, signal%radiometer, signal%sideband, &
      & signal%singleSideband, signal%spectrometer, signal%spectrometerType, &
      & signal%switch /), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Packing integers from signal' )
    call PVMIDLPack ( signal%deferred, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Packing signal deferred flag' )

    ! Now do each vector
    if ( associated ( signal%frequencies ) ) then
      call PVMIDLPack ( size(signal%frequencies), info )
      if ( info == 0 ) call PVMIDLPack ( signal%frequencies, info )
    else
      call PVMIDLPack ( 0, info )
    end if
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Packing signal frequencies/size' )

    if ( associated ( signal%widths ) ) then
      call PVMIDLPack ( size(signal%widths), info )
      if ( info == 0 ) call PVMIDLPack ( signal%widths, info )
    else
      call PVMIDLPack ( 0, info )
    end if
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Packing signal widths' )

    if ( associated ( signal%channels ) ) then
      call PVMIDLPack ( size(signal%channels), info )
      if ( info == 0 ) call PVMIDLPack ( signal%channels, info )
    else
      call PVMIDLPack ( 0, info )
    end if
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Packing signal channels' )

    ! That's it
  end subroutine PVMPackSignal
    
  ! ------------------------------------------------ PVMUnpackSignal ---
  subroutine PVMUnpackSignal ( signal )
    use PVMIDL, only: PVMIDLunpack
    ! Dummy arguments
    type ( Signal_T ), intent(out) :: SIGNAL

    ! Local variables
    integer :: INFO                     ! Flag from PVM
    integer, dimension(11) :: i11       ! Temporary array
    real(r8), dimension(2) :: f2        ! Temporary array
    integer :: SZ

    ! Executable code
    ! Unpack in mirror image to above pack
    call PVMIDLUnpack ( f2, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Unpacking signal lo/cf' )
    signal%centerFrequency = f2(1)
    signal%lo = f2(2)

    call PVMIDLUnpack ( i11, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Unpacking integers from signal' )
    signal%band = i11(1)
    signal%direction = i11(2)
    signal%index = i11(3)
    signal%instrumentModule = i11(4)
    signal%name = i11(5)
    signal%radiometer = i11(6)
    signal%sideband = i11(7)
    signal%singleSideband = i11(8)
    signal%spectrometer = i11(9)
    signal%spectrometerType = i11(10)
    signal%switch = i11(11)
    call PVMIDLUnpack ( signal%deferred, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Unpacking signal deferred flag' )

    ! Now do each vector
    call PVMIDLUnpack ( sz, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Unpacking size of frequencies' )
    call Allocate_test ( signal%frequencies, sz, 'signal%frequencies', ModuleName )
    if ( sz > 0 ) then
      call PVMIDLUnpack ( signal%frequencies, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, 'Unpacking frequencies' )
    end if
    
    call PVMIDLUnpack ( sz, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Unpacking size of widths' )
    call Allocate_test ( signal%widths, sz, 'signal%widths', ModuleName )
    if ( sz > 0 ) then
      call PVMIDLUnpack ( signal%widths, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, 'Unpacking widths' )
    end if
    
    call PVMIDLUnpack ( sz, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Unpacking size of channels' )
    call Allocate_test ( signal%channels, sz, 'signal%channels', ModuleName )
    if ( sz > 0 ) then
      call PVMIDLUnpack ( signal%channels, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, 'Unpacking channels' )
    end if
    
    ! That's it
  end subroutine PVMUnpackSignal

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MLSSignals_M

! $Log$
! Revision 2.116  2020/01/27 21:32:47  pwagner
! Dump_Signals may now dump only signal names if Details negative
!
! Revision 2.115  2018/04/19 00:53:09  vsnyder
! Remove USE statements and declarations for unused names.  Add OtherChannels
! to call to displaySignalName in DisplaySignalName_index.  Remove Channels
! pointer in Dump_Signals because it's not used.
!
! Revision 2.114  2018/03/05 19:24:39  pwagner
! Improved DUMP_SIGNALS
!
! Revision 2.113  2018/02/27 00:50:00  livesey
! Added the supportedModule functionality to support A-SMLS
!
! Revision 2.112  2017/09/15 15:44:18  livesey
! Updated to allow modules to be defferred until signal definition
!
! Revision 2.111  2016/09/21 00:38:09  pwagner
! Change appearance of Dump_Modules
!
! Revision 2.110  2016/07/27 22:14:59  pwagner
! Added NameString component to module type
!
! Revision 2.109  2016/07/25 23:15:41  pwagner
! Added IsAnyModuleSpacecraft function
!
! Revision 2.108  2015/03/28 01:18:43  vsnyder
! Added stuff to trace allocate/deallocate addresses
!
! Revision 2.107  2014/09/05 00:11:11  vsnyder
! More complete and accurate allocate/deallocate size tracking.  Get
! kinds from MLSKinds instead of MLSCommon.
!
! Revision 2.106  2014/06/02 23:15:35  livesey
! Fix bug with radiometer%singleSideband not being initialized
!
! Revision 2.105  2014/05/24 01:30:17  vsnyder
! Calculate upper bound of frequencies and widths correctly
!
! Revision 2.104  2014/02/07 02:28:06  vsnyder
! Fail gracefully in case a database didn't get allocated due to prior
! trouble.  Add more tracing at -g1 level.
!
! Revision 2.103  2014/01/11 01:41:02  vsnyder
! Decruftification
!
! Revision 2.102  2014/01/09 00:24:29  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.101  2013/12/12 02:08:36  vsnyder
! Use iterator to handle variables, and IF and SELECT constructs
!
! Revision 2.100  2013/11/06 01:47:26  pwagner
! May read instrument field of module; e.g. emls
!
! Revision 2.99  2013/10/09 01:08:03  vsnyder
! Add call to Evaluate_Variable
!
! Revision 2.98  2013/08/30 03:56:02  vsnyder
! Revise use of trace_begin and trace_end
!
! Revision 2.97  2013/08/23 23:27:48  pwagner
! Added function to return whether or not s/c is Aura
!
! Revision 2.96  2013/06/28 18:10:45  pwagner
! Correct erroneous GetSignalIndex; do we ever need signal_T%name?
!
! Revision 2.95  2013/06/12 02:12:23  vsnyder
! Cruft removal
!
! Revision 2.94  2012/05/01 20:52:09  vsnyder
! Remove unused assignment
!
! Revision 2.93  2011/06/02 19:22:35  pwagner
! Fixed bug in getRadiometerIndex
!
! Revision 2.92  2011/05/09 17:24:55  pwagner
! Converted to using switchDetail
!
! Revision 2.91  2011/01/29 00:51:38  vsnyder
! Add comments about units checking in front end
!
! Revision 2.90  2010/02/04 23:08:00  vsnyder
! Remove USE or declaration for unused names
!
! Revision 2.89  2009/10/26 17:08:02  pwagner
! Added GetRadiometerIndex
!
! Revision 2.88  2009/10/01 00:49:13  vsnyder
! Add OtherSignals to DisplaySignalName_*
!
! Revision 2.87  2009/09/25 02:44:14  vsnyder
! Added OtherChannels to dump and getSignalName routines to use specific
! channels instead of the ones in the data structure (because sometimes
! we have a signals database index plus channels from some other source)
!
! Revision 2.86  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.85  2008/09/29 22:56:40  vsnyder
! Add PRINT statement in Not_Used_Here to reduce compilation cascades
!
! Revision 2.84  2008/09/29 22:55:17  vsnyder
! Add DisplayRadiometer subroutine
!
! Revision 2.83  2007/05/22 02:28:28  vsnyder
! Don't use list-directed formatting for internal writes.  The standard
! allows a processor to insert any number of blanks it wishes, which might
! overflow the length of the internal file.
!
! Revision 2.82  2007/04/26 20:30:53  pwagner
! Bugfix for way ifc writes ints to strings
!
! Revision 2.81  2006/04/20 01:07:57  vsnyder
! Display a signal name given either a signal structure or an index
!
! Revision 2.80  2005/06/22 17:25:50  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.79  2005/04/06 23:15:31  vsnyder
! Add Before, Sideband and Channel arguments for Display_Signal_Name
!
! Revision 2.78  2005/03/15 23:48:55  pwagner
! PVMERRORMESSAGE now part of MLSMessageModule
!
! Revision 2.77  2005/01/12 03:07:37  vsnyder
! Added CHANNEL argument to GetNameOfSignal and GetSignalName
!
! Revision 2.76  2004/10/30 00:24:34  vsnyder
! Add GetSidebandStartStop
!
! Revision 2.75  2004/10/06 21:16:26  vsnyder
! Add MatchSignalPair, use it for MatchSignals
!
! Revision 2.74  2004/08/03 21:49:05  vsnyder
! Inching toward PFA
!
! Revision 2.73  2004/07/23 19:48:12  vsnyder
! Some cannonball polishing
!
! Revision 2.72  2004/07/23 18:35:17  vsnyder
! Dump low bound of channels array
!
! Revision 2.71  2004/05/29 02:45:28  vsnyder
! Add DisplaySignalName, fix a bug in GetNameOfSignal
!
! Revision 2.70  2004/05/01 04:07:44  vsnyder
! Rearranged some dumping stuff
!
! Revision 2.69  2004/04/16 00:44:24  livesey
! Added GetFirstChannel
!
! Revision 2.68  2004/03/24 23:08:58  livesey
! Bug fix in getBandName
!
! Revision 2.67  2004/03/24 01:02:17  livesey
! Slight change in GetBandName to make smls stuff easier.
!
! Revision 2.66  2004/03/22 18:22:59  livesey
! Bug fixes, only relevant for smls
!
! Revision 2.65  2004/02/11 02:24:18  livesey
! Added (probably unnecessary) initialization for DACS in
! SpectromterType_T
!
! Revision 2.64  2004/01/28 02:10:07  vsnyder
! Polish up some dump routines, other cosmetics
!
! Revision 2.63  2004/01/28 01:17:36  vsnyder
! Do spectrometerType%dacs = dacs BEFORE putting spectrometerType in the database
!
! Revision 2.62  2004/01/16 21:36:28  livesey
! Added the ability to defer the connection between bands and radiometers
! until you define the signal.  This is to support some SMLS related
! research work.
!
! Revision 2.61  2003/08/16 01:14:03  vsnyder
! Add optional 'polarization' field to 'radiometer' spec
!
! Revision 2.60  2003/07/23 21:51:36  pwagner
! Tried to fix problem with lower case dacs
!
! Revision 2.59  2003/07/23 18:04:32  livesey
! Ensure that the spectrometer names are capitalized when outputing them.
!
! Revision 2.58  2003/07/18 20:23:34  livesey
! Added DACS flag
!
! Revision 2.57  2003/05/16 02:44:18  vsnyder
! Removed USE's for unreferenced symbols
!
! Revision 2.56  2003/05/10 22:21:12  livesey
! Tried to calm down -g1..
!
! Revision 2.55  2003/03/07 03:17:50  livesey
! Add optional justchannels argument to DestroySignalDatabase
!
! Revision 2.54  2002/10/08 17:42:10  livesey
! Bug fixes in pack/unpack
!
! Revision 2.53  2002/10/08 00:09:12  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.52  2002/10/05 00:40:28  livesey
! Added pvm packing and unpacking of signals, and deep option on destroy
!
! Revision 2.51  2002/10/03 05:38:11  livesey
! Infinite loop fix
!
! Revision 2.50  2002/09/05 20:27:39  livesey
! Got rid of print statement left over from long ago.
!
! Revision 2.49  2002/07/17 06:00:40  livesey
! Added GetSidebandLoop routine
!
! Revision 2.48  2002/05/15 17:29:48  livesey
! Fixed channels related bug in MatchSignal
!
! Revision 2.47  2002/05/14 22:31:59  livesey
! Added singleSideband
!
! Revision 2.46  2002/05/03 22:38:41  livesey
! Added direction field to bands.
!
! Revision 2.45  2002/02/14 23:00:45  livesey
! Added justChannels optional argument to destroySignal
!
! Revision 2.44  2002/02/14 18:37:40  livesey
! Big fix for AreSignalsSuperset
!
! Revision 2.43  2002/02/13 23:57:34  livesey
! Tidied up a bit.  Channels doesn't need to be set for
! defered signals.  Made getSignalName skip channels field
! if appropriate.  However, possible bug lurking in channel
! printing loop.
!
! Revision 2.42  2001/11/09 23:14:08  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.41  2001/10/12 23:07:23  pwagner
! Added two inverse functions to signal (module) indexes
!
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
