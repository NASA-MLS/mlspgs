! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module PointingGrid_m

  ! Read the pointing grid file.  Make a database of pointing grids.
  ! Link them to and from the Signals database in MLSSignals_m

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use D_HUNT_M, only: HUNT
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error, MLSMSG_Info
  use MLSSignals_m, only: Bands, DestroySignalDatabase, GetNameOfSignal, &
    & MaxSigLen, Signals, Signal_T
  use Output_m, only: Blanks, MLSMSG_Level, Output, PrUnit

  ! More USEs below in each procedure, if they're only used therein.

  private
  ! Public procedures:
  public :: Open_Pointing_Grid_File, Read_Pointing_Grid_File
  public :: Close_Pointing_Grid_File, Destroy_Pointing_Grid_Database
  public :: Dump_Pointing_Grid_Database

  type, public :: OneGrid_T
    real(r8) :: Height                            ! Zeta, actually
    real(r8), pointer, dimension(:) :: Frequencies => NULL()
  end type OneGrid_T

  type, public :: PointingGrid_T
    type(signal_T), pointer, dimension(:) :: Signals => NULL()
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

contains

  ! ------------------------------------  Open_Pointing_Grid_File  -----
  subroutine Open_Pointing_Grid_File ( Filename, Lun )

    character(len=*), intent(in) :: Filename ! Name of the pointing grid file
    integer, intent(out) :: Lun              ! Logical unit number to read it

    logical :: Exist, Opened
    integer :: Status

    do lun = 20, 99
      inquire ( unit=lun, exist=exist, opened=opened )
      if ( exist .and. .not. opened ) EXIT
    end do
    if ( opened .or. .not. exist ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "No logical unit numbers available" )
    open ( unit=lun, file=filename, status='old', form='formatted', &
      & access='sequential', iostat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "Unable to open pointing grid file " // Filename )
  end subroutine Open_Pointing_Grid_File

  ! ------------------------------------  Read_Pointing_Grid_File  -----
  subroutine Read_Pointing_Grid_File ( Lun )
    use Machine, only: IO_Error
    use Parse_Signal_m, only: Parse_Signal
    use String_Table, only: Display_String
    use Toggles, only: Gen, Levels, Switches, Toggle
    use Trace_M, only: Trace_begin, Trace_end

    integer, intent(in) :: Lun               ! Logical unit number to read it

    logical, pointer, dimension(:) :: Channels ! Specified in a signal
    real(r8) :: Frequency                    ! Center frequency for the grid.
    !                                          Read from the input file.
    !                                          Should ultimately come from the
    !                                          signals database.
    real(r8) :: Height                       ! Zeta, actually, from the file
    integer, pointer, dimension(:) :: HowManyGrids  ! per radiometer batch
    integer, pointer, dimension(:) :: HowManySignals ! per radiometer batch
    integer :: HowManyRadiometers            ! gotten by counting the file
    integer :: I, N                          ! Loop inductor, subscript, temp
    character(len=MaxSigLen) :: Line         ! From the input file
    integer :: NumHeights                    ! Read from the input
    integer :: Sideband                      ! Specified in a signal
    integer :: SignalCount                   ! From Parse_Signal
    integer :: Status                        ! From read or allocate
    integer, pointer, dimension(:) :: Signal_Indices   ! From Parse_Signal, q.v.

    if ( toggle(gen) ) call trace_begin ( "Read_Pointing_Grid_File" )

    if ( associated(pointingGrids) ) call destroy_pointing_grid_database

    nullify ( howManyGrids, howManySignals )
    call allocate_test ( howManyGrids, size(signals), 'HowManyGrids', &
      & moduleName )
    call allocate_test ( howManySignals, size(signals), 'HowManySignals', &
      & moduleName )

    ! First, read through the file and count how much stuff is there.
    read ( lun, '(A)', end=98, err=99, iostat=status ) line
    howManyRadiometers = 0
outer1: do
      howManyRadiometers = howManyRadiometers + 1
      if ( howManyRadiometers > size(howManyGrids) ) then ! Double table sizes
        signal_indices => howManyGrids
        nullify ( howManyGrids )
        call allocate_test ( howManyGrids, 2*size(signal_indices), &
          & 'HowManyGrids', moduleName )
        howManyGrids(:size(signal_indices)) = signal_indices
        call deallocate_test ( signal_indices, 'Old HowManyGrids', moduleName )
        signal_indices => howManySignals
        nullify ( howManySignals )
        call allocate_test ( howManySignals, 2*size(signal_indices), &
          & 'HowManySignals', moduleName )
        howManySignals(:size(signal_indices)) = signal_indices
        call deallocate_test ( signal_indices, 'Old HowManySignals', moduleName )
      end if
      howManySignals(howManyRadiometers) = 0
      nullify(signal_indices)
      do
        line = adjustl(line)
        call parse_signal ( line, signal_indices, onlyCountEm = signalCount )
        if ( signalCount == 0 ) call MLSMessage ( MLSMSG_Error, &
          & moduleName, "Improper signal specification: " // trim(line) )
        howManySignals(howManyRadiometers) = &
          & howManySignals(howManyRadiometers) + signalCount
        read ( lun, '(A)', end=98, err=99, iostat=status ) line
        if ( verify(line(1:1), ' 0123456789.+-') == 0 ) EXIT ! a number
      end do
      read ( line, *, err=99 ) Frequency
      howManyGrids(howManyRadiometers) = 0
      do
        read ( lun, '(A)', err=99, iostat=status ) line
        if ( status < 0 ) EXIT outer1
        line = adjustl(line)
        if ( verify(line(1:1), ' 0123456789.+-') /= 0 ) EXIT ! not a number
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
    read ( lun, '(A)', iostat=status ) line  ! Read the first radiometer spec
    if ( status > 0 ) go to 99
    howManyRadiometers = 0
outer2: do
      howManyRadiometers = howManyRadiometers + 1
      allocate ( pointingGrids(howManyRadiometers)%signals( &
        & howManySignals(howManyRadiometers) ), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Allocate // "PointingGrids(?)%signals" )
      n = 0 ! Counter in pointingGrids(howManyRadiometers)%signals
      nullify ( signal_indices )
      nullify ( channels )
      do
        call parse_signal ( line, signal_indices, &
          & sideband=sideband, channels=channels )
        ! Errors were checked during the "counting" phase above
        do i = 1, size(signal_indices)
          n = n + 1
          pointingGrids(howManyRadiometers)%signals(n) = &
            & signals(signal_indices(i))
          pointingGrids(howManyRadiometers)%signals(n)%sideband = sideband
          pointingGrids(howManyRadiometers)%signals(n)%channels => channels
        end do
        call deallocate_test ( signal_indices, 'Signal_Indices', moduleName )
        read ( lun, '(A)', err=99, iostat=status ) line
        if ( verify(line(1:1), ' 0123456789.+-') == 0 ) EXIT ! a number
      end do
      !??? Should the centerFrequency be gotten from the signals database?
      !??? Maybe not.  The one in the pointing grid file appears to have been
      !??? Doppler shifted.
      read ( line, *, err=99, iostat=status ) &
        & pointingGrids(howManyRadiometers)%centerFrequency
      allocate ( pointingGrids(howManyRadiometers)% &
        & oneGrid(howManyGrids(howManyRadiometers)), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Allocate // "PointingGrids(?)%OneGrid" )
      n = 0
      do
        read ( lun, '(A)', err=99, iostat=status ) line
        if ( status < 0 ) EXIT outer2
        if ( verify(line(1:1), ' 0123456789.+-') /= 0 ) EXIT ! not a number
        n = n + 1
        read ( line, *, iostat=status, err=99 ) height, numHeights
        if ( status < 0 ) EXIT outer2
        if ( status /= 0 ) goto 99
        pointingGrids(howManyRadiometers)%oneGrid(n)%height = height
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

    call deallocate_test ( howManySignals, 'HowManySignals', moduleName )
    call deallocate_test ( howManyGrids, 'HowManyGrids', moduleName )

    if ( toggle(gen) ) then
      if ( levels(gen) > 0 .or. index(switches,'P') /= 0 ) &
        & call dump_pointing_grid_database
      call trace_end ( "Read_Pointing_Grid_File" )
    end if

    Return

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
    if (.not. associated(pointingGrids) ) return
    do i = 1, size(pointingGrids)
      ! It used to do this.
      !  call destroySignalDatabase ( pointingGrids(i)%signals )
      ! but that seemed to be a little overzelous, so I've replaced
      ! it with this, until I understand the issues better.
      do j = 1, size(pointingGrids(i)%signals)
        call Deallocate_Test ( pointingGrids(i)%signals(j)%channels, &
          & 'pointingGrids(?)%signals(?)%channels', ModuleName )
      end do
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
    use Dump_0, only: Dump

    integer :: I, J                     ! Subscripts, loop inductors
    character(len=MaxSigLen) :: SigName ! From GetSignalName
    call output ( 'Pointing Grids: SIZE = ' )
    call output ( size(pointingGrids), advance='yes' )
    do i = 1, size(pointingGrids)
      call output ( i, 4 )
      call output ( ':    Signals =', advance='yes' )
      do j = 1, size(pointingGrids(i)%signals)
        call blanks ( 6 )
        call getNameOfSignal ( pointingGrids(i)%signals(j), sigName )
        call output ( trim(sigName), advance='yes' )
      end do ! j = 1, size(pointingGrids(i)%signals)
      call output ( ' Center Frequency = ' )
      call output ( pointingGrids(i)%centerFrequency, advance='yes' )
      do j = 1, size(pointingGrids(i)%oneGrid)
        call output ( j, 4 )
        call output ( ':: Zeta = ' )
        call output ( pointingGrids(i)%oneGrid(j)%height )
        call dump ( pointingGrids(i)%oneGrid(j)%frequencies, &
          & '    Frequencies =' )
      end do ! j = 1, size(pointingGrids(i)%oneGrid)
    end do ! i
  end subroutine Dump_Pointing_Grid_Database

end module PointingGrid_m

! $Log$
! Revision 2.0  2001/09/17 20:26:26  livesey
! New forward model
!
! Revision 1.17.2.2  2001/09/14 20:14:49  livesey
! Rewrote Destroy_Pointing_Grid_Database to handle the signals
! differently.  Seemed to fix some memory stomping I didn't
! completely understand.
!
! Revision 1.17.2.1  2001/09/09 03:06:08  livesey
! Minor change
!
! Revision 1.17  2001/06/07 23:30:33  pwagner
! Added Copyright statement
!
! Revision 1.16  2001/05/04 00:49:43  livesey
! Let destroy quit if nothing to destroy
!
! Revision 1.15  2001/04/26 02:36:52  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 1.14  2001/04/20 02:55:56  zvi
! Get rid of the 1/48 griiding assumption in the file reader
!
! Revision 1.13  2001/04/19 06:48:14  zvi
! Fixing memory leaks..
!
! Revision 1.12  2001/04/13 22:50:27  livesey
! Tidied up, added some nullifies etc.
!
! Revision 1.11  2001/04/13 21:39:15  vsnyder
! Replace 'indices of signals' by 'array of signals'
!
! Revision 1.10  2001/03/30 01:34:33  vsnyder
! More work on I/O status handling
!
! Revision 1.9  2001/03/30 01:28:13  vsnyder
! Repair handling of I/O status values
!
! Revision 1.8  2001/03/29 23:53:21  vsnyder
! Use MaxSigLen parameter from MLSSignals
!
! Revision 1.7  2001/03/28 00:42:37  livesey
! Got rid of obsolete routines.
!
! Revision 1.6  2001/03/20 02:31:32  livesey
! Modified so always returns at least lowest grid rather than 0
!
! Revision 1.5  2001/03/17 19:17:57  vsnyder
! Destroy database if it exists before reading a new one
!
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
