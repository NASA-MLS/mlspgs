! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module PointingGrid_m

  ! Read the pointing grid file.  Make a database of pointing grids.
  ! Link them to and from the Signals database in MLSSignals_m
  
  ! ---------------------- TBD ----------------------------
  
  ! Should point to any relevant documentation on Pointing Grids
  ! How are the files formatted?
  ! How are they created?
  ! ---------------------- TBD ----------------------------

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use MLSKinds, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use MLSSignals_M, only: DisplaysignalName, Maxsiglen, Signals, Signal_T

  ! More USEs below in each procedure, if they're only used therein.

  private
  ! Public procedures:
  public :: Open_Pointing_Grid_File, Read_Pointing_Grid_File
  public :: Close_Pointing_Grid_File, Destroy_Pointing_Grid_Database
  public :: Dump_Pointing_Grid_Database, Dump_Pointing_Grid

  type, public :: OneGrid_T
    real(r8) :: Height                            ! Zeta, actually
    real(r8), allocatable, dimension(:) :: Frequencies
  end type OneGrid_T

  type, public :: PointingGrid_T
    type(signal_T), allocatable, dimension(:) :: Signals
    real(r8) :: CenterFrequency         ! Should be gotten from Bands database
    !??? Maybe not.  The one in the pointing grid file appears to have been
    !??? Doppler shifted.
    type(oneGrid_t), allocatable, dimension(:) :: OneGrid
  end type PointingGrid_T

  type(pointingGrid_t), public, pointer, dimension(:) :: &
    & PointingGrids => NULL()

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ------------------------------------  Open_Pointing_Grid_File  -----
  subroutine Open_Pointing_Grid_File ( Filename, Lun )

    character(len=*), intent(in) :: Filename ! Name of the pointing grid file
    integer, intent(out) :: Lun              ! Logical unit number to read it

    logical :: Exist, Opened
    character(127) :: IOMSG
    integer :: Status

    do lun = 20, 99
      inquire ( unit=lun, exist=exist, opened=opened )
      if ( exist .and. .not. opened ) EXIT
    end do
    if ( opened .or. .not. exist ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "No logical unit numbers available" )
    open ( unit=lun, file=filename, status='old', form='formatted', &
      & access='sequential', iostat=status, iomsg=iomsg )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "Unable to open pointing grid file " // Filename // ": " // trim(iomsg) )
  end subroutine Open_Pointing_Grid_File

  ! ------------------------------------  Read_Pointing_Grid_File  -----
  ! Should document format of Pointing_Grid_File or else
  ! point to location of file format.
  subroutine Read_Pointing_Grid_File ( Lun, Where )
    use Allocate_Deallocate, only: Test_Allocate
    use, Intrinsic :: Iso_C_Binding, only: C_Intptr_T, C_Loc
    use Machine, only: Io_Error
    use MLSStringlists, only: Switchdetail
    use Parse_Signal_M, only: Parse_Signal
    use Toggles, only: Gen, Levels, Switches, Toggle
    use Trace_M, only: Trace_Begin, Trace_End

    integer, intent(in) :: Lun               ! Logical unit number to read it
    integer, intent(in) :: Where             ! In the L2CF tree, for tracing

    integer(c_intptr_t) :: Addr         ! For tracing
    logical, pointer, dimension(:) :: Channels ! Specified in a signal
    real(r8) :: Frequency                    ! Center frequency for the grid.
    !                                          Read from the input file.
    !                                          Should ultimately come from the
    !                                          signals database.
    real(r8) :: Height                       ! Zeta, actually, from the file
    integer, allocatable, dimension(:) :: HowManyGrids  ! per radiometer batch
    integer, allocatable, dimension(:) :: HowManySignals ! per radiometer batch
    integer :: HowManyRadiometers            ! gotten by counting the file
    integer :: I, N                          ! Loop inductor, subscript, temp
    character(127) :: IOMSG
    character(len=MaxSigLen) :: Line         ! From the input file
    integer :: LineNo                        ! In the input file
    integer :: Me = -1                       ! String index for trace
    integer :: NumHeights                    ! Read from the input
    integer :: Sideband                      ! Specified in a signal
    integer :: SignalCount                   ! From Parse_Signal
    integer :: Status                        ! From read or allocate
    integer, pointer, dimension(:) :: Signal_Indices   ! From Parse_Signal, q.v.
    integer, allocatable :: Temp(:)          ! For increasing sizes
    character(3) :: Which                    ! Which I/O statement caused error

    call trace_begin ( me, "Read_Pointing_Grid_File", where, &
      & cond=toggle(gen) .and. levels(gen) > 0 )

    if ( associated(pointingGrids) ) call destroy_pointing_grid_database

    call allocate_test ( howManyGrids, size(signals), 'HowManyGrids', &
      & moduleName )
    call allocate_test ( howManySignals, size(signals), 'HowManySignals', &
      & moduleName )

    ! First, read through the file and count how much stuff is there.
    which = ' A ' ! in case of error
    lineNo = 1
    read ( lun, '(A)', end=98, err=99, iostat=status, iomsg=iomsg ) line
    howManyRadiometers = 0
outer1: do
      howManyRadiometers = howManyRadiometers + 1
      if ( howManyRadiometers > size(howManyGrids) ) then ! Double table sizes
        call allocate_test ( temp, 2*size(howManyGrids), &
          & 'HowManyGrids', moduleName )
        temp(:size(howManyGrids)) = howManyGrids
        call move_alloc ( temp, howManyGrids )
        call allocate_test ( temp, 2*size(howmanySignals), &
          & 'HowmanySignals', moduleName )
        temp(:size(howmanySignals)) = howmanySignals
        call move_alloc ( temp, howmanySignals )
      end if
      howManySignals(howManyRadiometers) = 0
      nullify(signal_indices)
      line = adjustl(line)
      do
        call parse_signal ( line, signal_indices, onlyCountEm = signalCount )
        if ( signalCount == 0 ) call MLSMessage ( MLSMSG_Error, &
          & moduleName, "Improper signal specification: " // trim(line) )
        howManySignals(howManyRadiometers) = &
          & howManySignals(howManyRadiometers) + signalCount
        which = ' B ' ! in case of error
        lineNo = lineNo + 1
        read ( lun, '(A)', end=98, err=99, iostat=status, iomsg=iomsg ) line
        line = adjustl(line)
        if ( verify(line(1:1), ' 0123456789.+-') == 0 ) EXIT ! a number
      end do
      which = ' C ' ! in case of error
      read ( line, *, err=99, iostat=status, iomsg=iomsg ) Frequency
      howManyGrids(howManyRadiometers) = 0
      do
        which = ' D ' ! in case of error
        lineNo = lineNo + 1
        read ( lun, '(A)', err=99, iostat=status, iomsg=iomsg ) line
        if ( status < 0 ) EXIT outer1
        line = adjustl(line)
        if ( verify(line(1:1), ' 0123456789.+-') /= 0 ) EXIT ! not a number
        howManyGrids(howManyRadiometers) = howManyGrids(howManyRadiometers) + 1
        which = ' E ' ! in case of error
        read ( line, *, err=99, iostat=status, iomsg=iomsg ) height, numHeights
        which = ' F ' ! in case of error
        lineNo = lineNo + 1
        read ( lun, *, err=99, end=98, iostat=status, iomsg=iomsg ) ( height, i = 1, numHeights )
      end do
    end do outer1

    rewind ( lun )
    ! Now that we know how much is there, allocate the data structures.
    allocate ( pointingGrids(howManyRadiometers), stat=status )
    addr = 0
    if ( status == 0 .and. howManyRadiometers > 0 ) &
      & addr = transfer(c_loc(pointingGrids(1)), addr)
    call test_allocate ( status, moduleName, "PointingGrids", &
      & uBounds = howManyRadiometers, elementSize = storage_size(pointingGrids) / 8, &
      & address=addr )
    which = ' G ' ! in case of error
    lineNo = 1
    read ( lun, '(A)', iostat=status, iomsg=iomsg ) line  ! Read the first radiometer spec
    if ( status > 0 ) go to 99
    howManyRadiometers = 0
outer2: do
      howManyRadiometers = howManyRadiometers + 1
      allocate ( pointingGrids(howManyRadiometers)%signals( &
        & howManySignals(howManyRadiometers) ), stat=status )
      addr = 0
      if ( status == 0 .and. howManySignals(howManyRadiometers) > 0 ) &
        & addr = transfer(c_loc(pointingGrids(howManyRadiometers)%signals(1)), addr)
      call test_allocate ( status, moduleName, "PointingGrids(?)%signals", &
        & uBounds = howManySignals(howManyRadiometers), &
        & elementSize = storage_size(pointingGrids(howManyRadiometers)%signals) / 8, &
        & address=addr )
      n = 0 ! Counter in pointingGrids(howManyRadiometers)%signals
      nullify ( signal_indices, channels )
      line = adjustl(line)
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
        which = ' H ' ! in case of error
        lineNo = lineNo + 1
        read ( lun, '(A)', err=99, iostat=status, iomsg=iomsg ) line
        line = adjustl(line)
        if ( verify(line(1:1), ' 0123456789.+-') == 0 ) EXIT ! a number
      end do
      !??? Should the centerFrequency be gotten from the signals database?
      !??? Maybe not.  The one in the pointing grid file appears to have been
      !??? Doppler shifted.
      which = ' I ' ! in case of error
      read ( line, *, err=99, iostat=status, iomsg=iomsg ) &
        & pointingGrids(howManyRadiometers)%centerFrequency
      allocate ( pointingGrids(howManyRadiometers)% &
        & oneGrid(howManyGrids(howManyRadiometers)), stat=status )
      addr = 0
      if ( status == 0 .and. howManyGrids(howManyRadiometers) > 0 ) &
        & addr = transfer(c_loc(pointingGrids(howManyRadiometers)%oneGrid(1)), addr)
      call test_allocate ( status, moduleName, "PointingGrids(?)%OneGrid", &
        & uBounds = howManyGrids(howManyRadiometers), &
        & elementSize = storage_size(pointingGrids(howManyRadiometers)%oneGrid) / 8, &
        & address=addr )
      n = 0
      do
        which = ' J ' ! in case of error
        lineNo = lineNo + 1
        read ( lun, '(A)', err=99, iostat=status, iomsg=iomsg ) line
        line = adjustl(line)
        if ( status < 0 ) EXIT outer2
        if ( verify(line(1:1), ' 0123456789.+-') /= 0 ) EXIT ! not a number
        n = n + 1
        which = ' K ' ! in case of error
        read ( line, *, iostat=status, err=99, iomsg=iomsg ) height, numHeights
        if ( status < 0 ) EXIT outer2
        if ( status /= 0 ) goto 99
        pointingGrids(howManyRadiometers)%oneGrid(n)%height = height
        call allocate_test ( &
          & pointingGrids(howManyRadiometers)%oneGrid(n)%frequencies, &
          & numHeights, "PointingGrids(?)%oneGrid(?)%frequencies", moduleName )
        which = ' L ' ! in case of error
        lineNo = lineNo + 1
        read ( lun, *, iostat=status, err=99, end=98, iomsg=iomsg ) &
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

    if ( switchDetail(switches,'point') > -1 ) call dump_pointing_grid_database
    call trace_end ( "Read_Pointing_Grid_File", &
      & cond=toggle(gen) .and. levels(gen) > 0 )

    return

  98 call MLSMessage ( MLSMSG_Error, moduleName, "Unexpected end-of-file" )
  99 write ( line, '(a,i0)' ) ' at line ', lineNo
     call io_error ( "While reading the pointing grid file" // trim(line) // &
                   & ": " // which // trim(iomsg), status )
     call MLSMessage ( MLSMSG_Error, moduleName, "Input error" // trim(line) //&
                     & which // trim(iomsg))
  end subroutine Read_Pointing_Grid_File

  ! -----------------------------------  Close_Pointing_Grid_File  -----
  subroutine Close_Pointing_Grid_File ( Lun )
    close ( lun )
  end subroutine Close_Pointing_Grid_File

  ! -----------------------------  Destroy_Pointing_Grid_Database  -----
  subroutine Destroy_Pointing_Grid_Database
    use Allocate_Deallocate, only: Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: I, J, S, Status
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
      s = size(pointingGrids(i)%oneGrid) * &
        & storage_size(pointingGrids(i)%oneGrid) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(pointingGrids(i)%oneGrid(1)), addr)
      deallocate ( pointingGrids(i)%oneGrid, stat=status )
      call test_deallocate ( status, moduleName, "PointingGrids(?)%oneGrid(?)", s, &
        & address=addr )
    end do ! i
    s = size(pointingGrids) * storage_size(pointingGrids) / 8
    addr = 0
    if ( s > 0 ) addr = transfer(c_loc(pointingGrids(1)), addr)
    deallocate ( pointingGrids, stat=status )
    call test_deallocate ( status, moduleName, "PointingGrids", s, address=addr )
  end subroutine Destroy_Pointing_Grid_Database

  ! --------------------------------  Dump_Pointing_Grid_Database  -----
  !
  ! This family of routines Dumps all or just one pointing frequency grid
  subroutine Dump_Pointing_Grid_Database ( where, details )
    use Dump_0, only: Dump
    use Moretree, only: StarterrorMessage
    use Output_M, only: Blanks, Output

    integer, intent(in), optional :: Where   ! Tree node index
    integer, intent(in), optional :: Details ! Show heights, freqs if > 0

    integer :: I, J                     ! Subscripts, loop inductors
    integer :: myDetails
    ! Executable
    myDetails = 0
    if ( present(Details) ) myDetails = Details
    if ( associated(pointingGrids) ) then
      call output ( 'Pointing Grids: SIZE = ' )
      call output ( size(pointingGrids), advance='yes' )
      do i = 1, size(pointingGrids)
        call output ( i, 4 )
        call output ( ':    Signals =', advance='yes' )
        call Dump_Pointing_Grid ( pointingGrids(i), Details )
      end do ! i
    else
      if ( present(where) ) call startErrorMessage ( where )
      call output ( 'No pointing grids database to dump.', advance='yes' )
    end if
  end subroutine Dump_Pointing_Grid_Database

  ! --------------------------------  Dump_Pointing_Grid  -----
  subroutine Dump_Pointing_Grid ( Grid, details )
    use Dump_0, only: Dump
    use Moretree, only: StarterrorMessage
    use Output_M, only: Blanks, Output

    type(PointingGrid_T), intent(in) :: Grid
    integer, intent(in), optional    :: Details ! Show heights, freqs if > 0

    integer :: J                     ! Subscripts, loop inductors
    integer :: myDetails
    ! Executable
    myDetails = 0
    if ( present(Details) ) myDetails = Details
    do j = 1, size(Grid%signals)
      call blanks ( 6 )
      call DisplaySignalName ( Grid%signals(j), advance='yes' )
    end do ! j = 1, size(Grid%signals)
    if ( myDetails < 1 ) return
    call output ( ' Center Frequency = ' )
    call output ( Grid%centerFrequency, advance='yes' )
    do j = 1, size(Grid%oneGrid)
      call output ( j, 4 )
      call output ( ':: Zeta = ' )
      call output ( Grid%oneGrid(j)%height )
      call dump ( Grid%oneGrid(j)%frequencies, &
        & '    Frequencies =' )
    end do ! j = 1, size(Grid%oneGrid)
  end subroutine Dump_Pointing_Grid

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module PointingGrid_m

! $Log$
! Revision 2.18  2020/01/27 18:17:26  pwagner
! Noted that docs are sorely needed; Dump_Pointing_Grid_Database now takes optional arg details
!
! Revision 2.17  2016/09/08 20:52:59  vsnyder
! Make components allocatable instead of pointers
!
! Revision 2.16  2015/03/28 02:00:28  vsnyder
! Added stuff to trace allocate/deallocate addresses
!
! Revision 2.15  2014/11/06 00:00:09  vsnyder
! Robustify: left adjust LINE before checking for numbers
!
! Revision 2.14  2014/09/05 20:51:33  vsnyder
! More complete and accurate allocate/deallocate size tracking
!
! Revision 2.13  2013/08/30 03:56:23  vsnyder
! Revise use of trace_begin and trace_end
!
! Revision 2.12  2011/05/09 17:52:15  pwagner
! Converted to using switchDetail
!
! Revision 2.11  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.10  2007/10/03 23:58:26  vsnyder
! Add 'where' for tracing
!
! Revision 2.9  2007/05/23 22:40:20  vsnyder
! Change tracing level of Read_Pointing_Grid_File from zero to one
!
! Revision 2.8  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.7  2004/05/29 02:49:51  vsnyder
! Simplifications from using DisplaySignalName
!
! Revision 2.6  2004/05/26 23:54:14  vsnyder
! Don't dump the database if it's not allocated
!
! Revision 2.5  2003/05/19 19:58:07  vsnyder
! Remove USEs for unreferenced symbols, remove unused local variables
!
! Revision 2.4  2003/05/10 22:20:57  livesey
! Tried to calm down -g1..
!
! Revision 2.3  2002/10/08 17:08:05  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.2  2002/09/13 22:58:22  vsnyder
! Cosmetic changes
!
! Revision 2.1  2002/05/08 08:53:43  zvi
! All radiometers grid concept implementation
!
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
