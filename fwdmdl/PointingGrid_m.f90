module PointingGrid_m

  ! Read the pointing grid file.  Make a database of pointing grids.
  ! Link them to and from the Signals database in MLSSignals_m

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use Dump_0, only: Dump
  use Machine, only: IO_Error
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error, MLSMSG_Info
  use MLSSignals_m, only: Bands, Dump_Bands, Dump_Signals, GetSignalName, &
    & Signals
  use Output_m, only: Blanks, MLSMSG_Level, Output, PrUnit
  use Parse_Signal_m, only: Parse_Signal
  use String_Table, only: Display_String
  use Toggles, only: Gen, Levels, Switches, Toggle
  use Trace_M, only: Trace_begin, Trace_end

  private
  ! Public procedures:
  public :: Open_Pointing_Grid_File, Read_Pointing_Grid_File
  public :: Close_Pointing_Grid_File, Destroy_Pointing_Grid_Database
  public :: Dump_Pointing_Grid_Database

  type, public :: OneGrid_T
    real(r8) :: Height                  ! Zeta, actually
    real(r8), pointer, dimension(:) :: Frequencies => NULL()
  end type OneGrid_T

  type, public :: PointingGrid_T
    integer, pointer, dimension(:) :: Signals => NULL() ! Database indices
    real(r8) :: CenterFrequency         ! Should be gotten from Bands database
    !??? Maybe not.  The one in the pointing grid file appears to have been
    !??? Doppler shifted.
    type(oneGrid_t), pointer, dimension(:) :: OneGrid => NULL()
  end type PointingGrid_T

  integer, public :: ExtraHeights = 0   ! How many extra spaces in each
    !                                     frequency grid?
  type(pointingGrid_t), public, pointer, dimension(:) :: &
    & PointingGrids => NULL()

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter, private :: ModuleName = &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains

  ! ------------------------------------  Open_Pointing_Grid_File  -----
  subroutine Open_Pointing_Grid_File ( Filename, Lun )
    character(len=*), intent(in) :: Filename ! Name of the pointing grid file
    integer, intent(out) :: Lun              ! Logical unit number to read it
    logical :: Exist, Opened
    integer :: Status
    do lun = 20, 99
      inquire ( unit=lun, exist=exist, opened=opened )
      if ( exist .and. .not. opened ) exit
    end do
    if ( opened .or. .not. exist ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "No logical unit numbers available" )
    open ( unit=lun, file=filename, status='old', form='formatted', &
      & access='sequential', iostat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "Unable to open pointing grid file " // Filename )
  end subroutine Open_Pointing_Grid_File

  ! ------------------------------------  Read_Pointing_Grid_File  -----
  subroutine Read_Pointing_Grid_File ( Lun, Spec_Indices )
    integer, intent(in) :: Lun               ! Logical unit number to read it
    integer, intent(in) :: Spec_Indices(:)   ! Needed by Parse_Signal, q.v.

    real(r8) :: Frequency                    ! Center frequency for the grid.
    !                                          Read from the input file.
    !                                          Should ultimately come from the
    !                                          signals database.
    real(r8) :: Height                       ! Zeta, actually, from the file
    integer, dimension(size(signals)) :: HowManyGrids  ! per radiometer
    integer :: HowManyRadiometers            ! gotten by counting the file
    integer :: I, N                          ! Loop inductor, subscript, temp
    character(len=80) :: Line                ! From the input file
    integer :: NumHeights                    ! Read from the input
    integer :: Status                        ! From read or allocate
    integer, pointer, dimension(:) :: Signal_Indices   ! From Parse_Signal, q.v.

    if ( toggle(gen) ) call trace_begin ( "Read_Pointing_Grid_File" )

    ! First, read through the file and count how much stuff is there.
    read ( lun, '(a)', iostat=status ) line  ! Skip the first radiometer spec
    howManyRadiometers = 0
outer1: do
      howManyRadiometers = howManyRadiometers + 1
      if ( howManyRadiometers > size(howManyGrids) ) &
        & call MLSMessage ( MLSMSG_Error, moduleName, &
          & "More radiometers in the file than signals in the database" )
      read ( lun, *, err=99 ) Frequency
      howManyGrids(howManyRadiometers) = 0
      do
        read ( lun, '(a)', iostat=status ) line
        if ( status < 0 ) exit outer1
        if ( status > 0 ) go to 98
        if ( verify(line(1:1), ' 0123456789.+-') /= 0 ) exit ! not a number
        howManyGrids(howManyRadiometers) = howManyGrids(howManyRadiometers) + 1
        read ( line, *, err=99 ) height, numHeights
        read ( lun, *, err=99, end=98 ) ( height, i = 1, numHeights )
      end do
    end do outer1

    rewind ( lun )
    ! Now that we know how much is there, allocate the data structures.
    allocate ( pointingGrids(howManyRadiometers), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_Allocate // "PointingGrids" )
    read ( lun, '(a)', iostat=status ) line  ! Read the first radiometer spec
    howManyRadiometers = 0
outer2: do
      howManyRadiometers = howManyRadiometers + 1
      nullify ( signal_indices )
      call parse_signal ( line, signal_indices, spec_indices )
      if ( .not. associated(signal_indices) ) call MLSMessage ( MLSMSG_Error, &
        & moduleName, "Improper signal specification" )
      ! Check that pointing grids have not already been specified for any
      ! bands implied by the present signal string.
      do i = 1, size(signal_indices)
        if ( signals(signal_indices(i))%pointingGrid /= 0 ) then
          prunit = -2 ! To do output via MLSMessage
          MLSMSG_Level = MLSMSG_Info
          call dump_signals ( (/ signals(signal_indices(i)) /) )
          MLSMSG_Level = MLSMSG_Error
          call output ( "More than one pointing grid specified for signal " )
          call display_string ( signals(signal_indices(i))%name, advance='yes' )
        end if
        if ( bands(signals(signal_indices(i))%band)%pointingGrid /= 0 ) then
          prunit = -2 ! To do output via MLSMessage
          MLSMSG_Level = MLSMSG_Info
          call dump_bands ( (/ bands(signals(signal_indices(i))%band) /) )
          MLSMSG_Level = MLSMSG_Error
          call output ( "More than one pointing grid specified for band " )
          call display_string ( bands(signals(signal_indices(i))%band)%prefix, &
            & advance='yes' )
        end if
      end do
      ! But allow the present signal string to specify the grid for a band
      ! several times.
      signals(signal_indices)%pointingGrid = howManyRadiometers
      bands(signals(signal_indices)%band)%pointingGrid = howManyRadiometers
      pointingGrids(howManyRadiometers)%signals => signal_indices
      ! The next thing should be gotten from the signals database.
      !??? Maybe not.  The one in the pointing grid file appears to have been
      !??? Doppler shifted.
      read ( lun, *, err=99 ) pointingGrids(howManyRadiometers)%centerFrequency
      allocate ( pointingGrids(howManyRadiometers)% &
        & oneGrid(howManyGrids(howManyRadiometers)), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Allocate // "PointingGrids(?)%OneGrid" )
      n = 0
      do
        read ( lun, '(a)', iostat=status ) line
        if ( status < 0 ) exit outer2
        if ( status > 0 ) go to 98
        if ( verify(line(1:1), ' 0123456789.+-') /= 0 ) exit ! not a number
        n = n + 1
        read ( line, *, iostat=status, err=99 ) & 
          & pointingGrids(howManyRadiometers)%oneGrid(n)%height, numHeights
        call allocate_test ( &
          & pointingGrids(howManyRadiometers)%oneGrid(n)%frequencies, &
          & numHeights + extraHeights, &
          & "PointingGrids(?)%oneGrid(?)%frequencies", moduleName )
          pointingGrids(howManyRadiometers)%oneGrid(n)% &
            & frequencies(:extraHeights) = -huge(1.0_r8)
        read ( lun, *, iostat=status, err=99, end=98 ) &
          & pointingGrids(howManyRadiometers)%oneGrid(n)% &
            & frequencies(extraHeights+1:)
      end do
    end do outer2

    if ( toggle(gen) ) then
      if ( levels(gen) > 0 .or. index(switches,'P') /= 0 ) &
        & call dump_pointing_grid_database
      call trace_end ( "Read_Pointing_Grid_File" )
    end if

    return
  98 call MLSMessage ( MLSMSG_Error, moduleName, "Unexpected end-of-file" )
  99 call io_error ( "While reading the pointing grid file", status )
     call MLSMessage ( MLSMSG_Error, moduleName, "Input error" )
  end subroutine Read_Pointing_Grid_File

  ! -----------------------------------  Close_Pointing_Grid_File  -----
  subroutine Close_Pointing_Grid_File ( Lun )
    close ( lun )
  end subroutine Close_Pointing_Grid_File

  ! -----------------------------  Destroy_Pointing_Grid_Database  -----
  subroutine Destroy_Pointing_Grid_Database
    integer :: I, J, Status
    do i = 1, size(pointingGrids)
      call deallocate_test ( pointingGrids(i)%signals, &
        & "PointingGrids(?)%signals", moduleName )
      do j = 1, size(pointingGrids(i)%oneGrid)
        call deallocate_test ( pointingGrids(i)%oneGrid(j)%frequencies, &
          & "PointingGrids(?)%oneGrid(?)%frequencies", moduleName )
      end do ! j
      deallocate ( pointingGrids(i)%oneGrid, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        MLSMSG_DeAllocate // "PointingGrids(?)%oneGrid(?)" )
    end do ! i
    deallocate ( pointingGrids, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      MLSMSG_DeAllocate // "PointingGrids" )
  end subroutine Destroy_Pointing_Grid_Database

  ! --------------------------------  Dump_Pointing_Grid_Database  -----
  subroutine Dump_Pointing_Grid_Database
    integer :: I, J                ! Subscripts, loop inductors
    character(len=80) :: SigName   ! From GetSignalName
    call output ( 'Pointing Grids: SIZE = ' )
    call output ( size(pointingGrids), advance='yes' )
    do i = 1, size(pointingGrids)
      call output ( i, 4 )
      call output ( ':    Signals =', advance='yes' )
      do j = 1, size(pointingGrids(i)%signals)
        call blanks ( 6 )
        call getSignalName ( pointingGrids(i)%signals(j), sigName )
        call output ( trim(sigName), advance='yes' )
      end do ! j = 1, size(pointingGrids(i)%signals)
      call output ( ' Center Frequency = ' )
      call output ( pointingGrids(i)%centerFrequency, advance='yes' )
      do j = 1, size(pointingGrids(i)%oneGrid)
        call output ( j, 4 )
        call output ( ':: Height = ' )
        call output ( pointingGrids(i)%oneGrid(j)%height )
        call dump ( pointingGrids(i)%oneGrid(j)%frequencies, &
          & '    Frequencies =' )
      end do ! j = 1, size(pointingGrids(i)%oneGrid)
    end do ! i
  end subroutine Dump_Pointing_Grid_Database

end module PointingGrid_m

! $Log$
