module PointingGrid_m

  ! Read the pointing grid file.  Make a database of pointing grids.
  ! Link them to and from the Signals database in MLSSignals_m

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use D_HUNT_M, only: HUNT
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error, MLSMSG_Info
  use MLSSignals_m, only: Bands, Dump_Bands, Dump_Signals, GetSignalName, &
    & Signals
  use Output_m, only: Blanks, MLSMSG_Level, Output, PrUnit

  ! More USEs below in each procedure, if they're only used therein.

  private
  ! Public procedures:
  public :: Open_Pointing_Grid_File, Read_Pointing_Grid_File
  public :: Close_Pointing_Grid_File, Destroy_Pointing_Grid_Database
  public :: Dump_Pointing_Grid_Database, Get_Grids_Near_Tan_Pressures
  public :: Get_Nearest_Tan_Pressure

  type, public :: OneGrid_T
    real(r8) :: Height                            ! Zeta, actually
    real(r8) :: NearestTanPress = huge(0.0_r8)/4  ! Nearest tangent pressure
    integer :: WhichTanPress = -1                 ! Which one is nearest
    real(r8), pointer, dimension(:) :: Frequencies => NULL()
  end type OneGrid_T

  type, public :: PointingGrid_T
    integer, pointer, dimension(:) :: Signals => NULL() ! Database indices
    real(r8) :: CenterFrequency         ! Should be gotten from Bands database
    !??? Maybe not.  The one in the pointing grid file appears to have been
    !??? Doppler shifted.
    type(oneGrid_t), pointer, dimension(:) :: OneGrid => NULL()
  end type PointingGrid_T

  type(pointingGrid_t), public, pointer, dimension(:) :: &
    & PointingGrids => NULL()

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter, private :: ModuleName = &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

  ! Has Get_Nearest_Tan_Pressure been called?  (Get_Grids_Near_Tan_Pressures
  ! needs it to have been done.)
  logical, private, save :: Got_Nearest_Tan_Press = .false.

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
    use Machine, only: IO_Error
    use Parse_Signal_m, only: Parse_Signal
    use String_Table, only: Display_String
    use Toggles, only: Gen, Levels, Switches, Toggle
    use Trace_M, only: Trace_begin, Trace_end

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

    if ( associated(pointingGrids) ) call destroy_pointing_grid_database

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
          & numHeights, "PointingGrids(?)%oneGrid(?)%frequencies", moduleName )
        read ( lun, *, iostat=status, err=99, end=98 ) &
          & pointingGrids(howManyRadiometers)%oneGrid(n)%frequencies
        ! The frequencies are relative to the band center.  Make them
        ! absolute
        pointingGrids(howManyRadiometers)%oneGrid(n)%frequencies = &
            & pointingGrids(howManyRadiometers)%oneGrid(n)%frequencies + &
                & pointingGrids(howManyRadiometers)%centerFrequency
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
    Got_Nearest_Tan_Press = .false.
  end subroutine Destroy_Pointing_Grid_Database

  ! --------------------------------  Dump_Pointing_Grid_Database  -----
  subroutine Dump_Pointing_Grid_Database
    use Dump_0, only: Dump

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

  ! -------------------------------  Get_Grids_Near_Tan_Pressures  -----
  subroutine Get_Grids_Near_Tan_Pressures ( Which_Pointing_Grid, Tan_Press, &
    & Tol, Grids )
  ! Get the index of the grid nearest to each tangent pressure, and
  ! within TOL of it.  Assumes that Get_Nearest_Tan_Pressure has been
  ! called already.

    integer, intent(in) :: Which_Pointing_Grid    ! Database index
    real(r8), intent(in), dimension(:) :: Tan_Press
    real(r8), intent(in) :: Tol                   ! How close?
    integer, intent(out) :: Grids(:)              ! Database indices

    integer :: I

    call Get_Nearest_Tan_Pressure ( tan_press ) ! Just in case

    grids = 0
    do i = 1, size(pointingGrids(which_pointing_grid)%oneGrid)
      if ( abs(pointingGrids(which_pointing_grid)%oneGrid(i)%height - &
        &      tan_press(pointingGrids(which_pointing_grid)%oneGrid(i)% &
               & whichTanPress) ) <= tol ) then
        if ( grids(pointingGrids(which_pointing_grid)%oneGrid(i)% &
          & whichTanPress) /= 0 ) &
            & call MLSMessage ( MLSMSG_Error, moduleName, &
              & "More than one frequency grid for a tangent height" )
        grids(pointingGrids(which_pointing_grid)%oneGrid(i)%whichTanPress) = i
      end if
    end do ! i = 1, size(pointingGrids(which_pointing_grid))
  end subroutine Get_Grids_Near_Tan_Pressures
  
  ! -----------------------------------  Get_Nearest_Tan_Pressure  -----
  subroutine Get_Nearest_Tan_Pressure ( Tan_Press )
  ! Find the nearest tangent pressure to each PointingGrids%oneGrid

    real(r8), intent(in), dimension(:) :: Tan_Press

    integer :: I, J, K, L

    if ( got_Nearest_Tan_Press ) return ! No point in doing it twice

    Got_Nearest_Tan_Press = .true.

    do i = 1, size(pointingGrids)
      do j = 1, size(pointingGrids(i)%oneGrid)
        k = -1
        call Hunt ( pointingGrids(i)%oneGrid(j)%height, tan_press, &
          & size(tan_press), k, l )
        if ( abs( pointingGrids(i)%oneGrid(j)%height - tan_press(l) ) < &
             abs( pointingGrids(i)%oneGrid(j)%height - tan_press(k) ) ) k = l
        pointingGrids(i)%oneGrid(j)%nearestTanPress = tan_press(k)
        pointingGrids(i)%oneGrid(j)%whichTanPress = k
      end do ! j = 1, size(pointingGrids%oneGrid)
    end do ! i = 1, size(pointingGrids)

  end subroutine Get_Nearest_Tan_Pressure
end module PointingGrid_m

! $Log$
! Revision 1.4  2001/03/17 02:34:03  vsnyder
! Get rid of "ExtraHeights" -- That's related to tan_press, not frq_grid
!
! Revision 1.3  2001/03/17 01:22:23  vsnyder
! Add Get_Grids_Near_Tan_Pressures
!
! Revision 1.2  2001/03/17 00:24:56  vsnyder
! Add Get_Nearest_Tan_Pressure
!
! Revision 1.1  2001/03/16 23:12:14  vsnyder
! Initial commit
!
