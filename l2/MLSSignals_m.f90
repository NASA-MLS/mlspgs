module MLSSignals_M

  ! Process the MLSSignals section of the L2 configuration file and deal with
  ! parsing signal request strings.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use Declaration_Table, only: Range
  use Expr_M, only: Expr
  use Init_Tables_Module, only: f_band, f_channel, f_channels, f_first, &
    & f_frequencies, f_frequency, f_last, f_lo, f_radiometer, &
    & f_spectrometer, f_start, f_step, f_suffix, &
    & f_switch, f_width, f_widths, field_first, field_indices, field_last, &
    & s_band, s_channel,  s_radiometer,  s_signal, s_spectrometerType
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
             & MLSMSG_Error
  use Output_M, only: Output
  use String_Table, only: Display_String
  use Toggles, only: Gen, Toggle
  use Trace_M, only: Trace_begin, Trace_end
  use Tree, only: Decorate, Decoration, Node_ID, Nsons, Source_Ref, Sub_Rosa, &
    & Subtree
  use Tree_Types, only: N_named

  ! This first type defines a radiometer.

  type Radiometer_T
    real(r8) :: LO              ! Local oscillator in MHz
    integer :: Suffix           ! Sub_rosa index
    integer :: InstrumentModule ! Literal value L_THz, L_GHz
  end type Radiometer_T

  ! The second type describes a band within that radiometer

  type Band_T
    real(r8) :: Frequency               ! Negative if not present (wide filter)
    integer :: Suffix                   ! Sub_rosa index
  end type Band_T

  ! This type gives the information for specific spectrometer families.  For
  ! all apart from the WF4 spectrometers, we list frequencies and widths. 
  ! Otherwise, the arrays are empty.

  type SpectrometerType_T
    real(r8), pointer, dimension(:) :: Frequencies => NULL(), Widths => NULL()
  end type SpectrometerType_T

  ! This is the key type, it describes a complete signal (one band, or a
  ! subset of the channels in one band).  We have to main variables of this
  ! type flying around, see later.

  type Signal_T
    integer :: Band                     ! Index in Bands data base
    integer :: Radiometer               ! Index in Radiometers data base
    integer :: SideBand                 ! -1:lower, +1:upper, 0:folded
    integer :: Spectrometer             ! Just a spectrometer number
    integer :: SpectrometerType       ! Index in SpectrometerType data base
    integer :: Switch                   ! Just a switch number
    real(r8), POINTER, DIMENSION(:) :: Frequencies=>NULL() ! Mainly a shallow copy
    real(r8), POINTER, DIMENSION(:) :: Widths=>NULL() ! Mainly a shallow copy
    logical, POINTER, DIMENSION(:) :: Selected ! A bit field, or absent==all.
  end type Signal_T

  ! Now some databases, the first are fairly obvious.

  type(band_T), save, pointer, dimension(:) :: Bands => NULL()
  type(radiometer_T), save, pointer, dimension(:) :: Radiometers => NULL()
  type(spectrometerType_T), save, pointer, dimension(:) ::&
    & SpectrometerTypes => NULL()

  ! The next two are important, the first describes all the `allowed' signals
  ! in the instrument.  The second is signals defined in the l2cf for later use
  ! in identifying bands or subsets of a band.

  type(signal_T), save, pointer, dimension(:) :: ValidSignals => NULL()

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
    integer :: Error                    ! Error level seen so far
    integer :: Field                    ! Field index -- f_something
    integer :: First                    ! "first" field of "spectrometer"
    integer :: Frequencies              ! subtree index of field
    logical :: Got(field_first:field_last)   ! "Got this field already"
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
        band%frequency = -1.0_r8 ! "The 'frequency' field is absent"
        do j = 2, nsons(key)
          son = subtree(j,key)
          field = decoration(subtree(1,son))
          got(field) = .true.
          select case ( field )
          case ( f_frequency )
            call expr ( subtree(2,son), units, value )
            band%frequency = value(1)
          case ( f_suffix )
            band%suffix = sub_rosa(son)
          case default
            ! Shouldn't get here if the type checker worked
          end select
        end do ! i = 2, nsons(key)
        call decorate ( key, addBandToDatabase ( bands, band ) )
      case ( s_spectrometerType ) ! ......................... SPECTROMETER TYPE  .....
        do j = 2, nsons(key)
          son = subtree(j,key)
          field = decoration(subtree(1,son))
          got(field) = .true.
          select case ( field )
          case ( f_frequencies )
            frequencies = son
          case ( f_widths )
            widths = son
          case default
            ! Shouldn't get here if the type checker worked
          end select
        end do ! i = 2, nsons(key)
        if ( nsons(frequencies) /= nsons(widths) ) &
          call announceError ( sizes, f_frequency, (/ f_widths /) )
        if ( error == 0 ) then
          call allocate_Test ( spectrometerType%frequencies, nsons(son)-1, &
            & 'spectrometerType%frequencies', moduleName, lowBound = first )
          call allocate_Test ( spectrometerType%widths, nsons(son)-1, &
            & 'spectrometerType%widths', moduleName, lowBound = first )
          do k = 2, nsons(frequencies)
            call expr ( subtree(k,frequencies), units, value )
            spectrometerType%frequencies(k-2+first) = value(1)
            call expr ( subtree(k,widths), units, value )
            spectrometerType%widths(k-2+first) = value(1)
          end do
          call decorate ( key, addSpectrometerTypeToDatabase (&
            & spectrometerTypes, spectrometerType ) )
          call destroySpectrometerType ( SpectrometerType )
        end if
      case ( s_radiometer ) ! .......................  RADIOMETER  .....
        do j = 2, nsons(key)
          son = subtree(j,key)
          field = decoration(subtree(1,son))
          select case ( field )
          case ( f_lo )
            call expr ( subtree(2,son), units, value )
            radiometer%lo = value(1)
          case ( f_suffix )
            radiometer%suffix = sub_rosa(son)
          case default
            ! Shouldn't get here if the type checker worked
          end select
        end do ! i = 2, nsons(key)
        call decorate ( key, addRadiometerToDatabase ( radiometers, radiometer ) )
!       case ( s_signal ) ! ...............................  SIGNAL  .....
!         do j = 2, nsons(key)
!           son = subtree(j,key)
!           field = decoration(subtree(1,son))
!           select case ( field )
!           case ( f_band )
!             signal%band = decoration(subtree(2,son))
!           case ( f_spectrometerType )
!             call expr ( subtree(2,son), units, value )
!             signal%spectrometerType = value(1)
!           case ( f_radiometer )
!             signal%radiometer = decoration(subtree(2,son))
!           case ( f_spectrometer )
!             signal%spectrometer = decoration(subtree(2,son))
!           case ( f_switch )
!             call expr ( subtree(2,son), units, value )
!             signal%switch = value(1)
!           case default
!             ! Shouldn't get here if the type checker worked
!           end select
!         end do ! i = 2, nsons(key)
!         call decorate ( key, addSignalToDatabase ( signals, signal ) )
!       case ( s_spectrometerType ) ! .............  SPECTROMETERTYPE .....
!         first = 0
!         do j = 2, nsons(key)
!           son = subtree(j,key)
!           field = decoration(subtree(1,son))
!           got(field) = .true.
!           select case ( field )
!           case ( f_spectrometerTypes )
!             mySpectrometerTypes = son
!           case ( f_first, f_last, f_start, f_step, f_width )
!             call expr ( subtree(2,son), units, value )
!             select case ( field )
!             case ( f_first )
!               first = value(1)
!             case ( f_last )
!               last = value(1)
!             case ( f_start )
!               start = value(1)
!             case ( f_step )
!               step = value(1)
!             case ( f_width )
!               width = value(1)
!             end select
!           case ( f_frequencies )   ! Postpone processing until later, so that
!             frequencies = son      ! we can verify that Frequencies and
!           case ( f_widths )        ! Widths have the same number of values
!             widths = son
!           case default
!             ! Shouldn't get here if the type checker worked
!           end select
!         end do ! i = 2, nsons(key)
!         if ( got(f_frequencies) .neqv. got(f_widths) ) call announceError ( &
!           & allOrNone, f_frequencies, (/ f_widths /) )
!         if ( got(f_spectrometerTypes) ) then
!           if ( any(got( (/ f_frequencies, f_last, f_start, &
!             & f_step, f_width, f_widths /) )) ) call announceError ( badMix, &
!             & f_spectrometerTypes, (/ f_frequencies, f_last, f_start, f_step, f_width, &
!             &                f_widths /) )
!           if ( error == 0 ) then
!             call allocate_Test ( spectrometer%spectrometerTypes, nsons(mySpectrometerTypes)-1, &
!               & 'spectrometer%spectrometerTypes', moduleName, lowBound = first )
!             do k = first, nsons(mySpectrometerTypes)-2+first
!               spectrometer%spectrometerTypes(k) = decoration(subtree(k,mySpectrometerTypes))
!             end do
!           end if
!         end if
!         if ( got(f_frequencies) ) then
!           if ( any(got( (/ f_last, f_start, f_step, f_width /) )) ) &
!             & call announceError ( badMix, f_frequencies, &
!             & (/ f_last, f_start, f_step, f_width /) )
!           if ( nsons(frequencies) /= nsons(widths) ) &
!             call announceError ( sizes, f_frequency, (/ f_widths /) )
!           if ( error == 0 ) then
!             call allocate_Test ( spectrometer%frequencies, nsons(son)-1, &
!               & 'spectrometer%frequencies', moduleName, lowBound = first )
!             call allocate_Test ( spectrometer%widths, nsons(son)-1, &
!               & 'spectrometer%widths', moduleName, lowBound = first )
!             do k = 2, nsons(frequencies)
!               call expr ( subtree(k,frequencies), units, value )
!               spectrometer%frequencies(k-2+first) = value(1)
!               call expr ( subtree(k,widths), units, value )
!               spectrometer%frequencies(k-2+first) = value(1)
!             end do
!           end if
!         end if
!         if ( got(f_start) ) then
!           if ( .not. all( got((/ f_last, f_step, &
!           & f_width /)) ) ) call announceError ( allOrNone, f_start, &
!           & (/ f_last, f_step, f_width /) )
!           if ( error == 0 ) then
!             k = last - first + 1
!             call allocate_Test ( spectrometer%frequencies, k, &
!               & 'spectrometer%frequencies', moduleName, lowBound = first )
!             call allocate_Test ( spectrometer%widths, k, &
!               & 'spectrometer%widths', moduleName, lowBound = first )
!             spectrometer%widths = width
!             spectrometer%frequencies(first) = start
!             do k = first+1, last
!               spectrometer%frequencies(k) = start + (k-first) * step
!             end do ! k
!           end if
!         end if
!         if ( .not. any(got( (/ f_spectrometerTypes, f_frequencies, f_start /) )) ) &
!           & call announceError ( atLeastOne, f_frequencies, (/ f_start /) )
!         if ( error == 0 ) call decorate ( key, addSpectrometerToDatabase ( &
!           & spectrometers, spectrometer ) )
!         call destroySpectrometer ( spectrometer )
      case default
        ! Shouldn't get here if the type checker worked
      end select

    end do

    if ( toggle(gen) ) call trace_end ( "MLSSignals" )

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

  ! ----------------------------------  AddSpectrometerToDatabase  -----
  integer function AddSpectrometerTypeToDatabase ( Database, Item )
    type(spectrometerType_T), dimension(:), pointer :: Database
    type(spectrometerType_T), intent(in) :: Item

    ! Local variables
    type (spectrometerType_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddSpectrometerToDatabase = newSize
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


end module MLSSignals_M

! $Log$
! Revision 2.4  2001/02/27 00:23:48  livesey
! Very interim version
!
! Revision 2.3  2001/02/15 22:03:00  vsnyder
! Remove checking for ranges -- the type checker does it
!
