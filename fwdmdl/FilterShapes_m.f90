module FilterShapes_m

  ! Read the filter shapes file.  Attach each array of filter shapes to
  ! the database for the specified signal from the Signals database in
  ! MLSSignals_m

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error
  use MLSSignals_m, only: MaxSigLen, Signals
  use Output_m, only: Output

  ! More USEs below in each procedure, if they're only used therein.

  private
  ! Public procedures:
  public :: Open_Filter_Shapes_File, Read_Filter_Shapes_File
  public :: Close_Filter_Shapes_File, Destroy_Filter_Shapes_Database
  public :: Dump_Filter_Shapes_Database

  type, public :: FilterShape_T
    real(r8) :: LHS, RHS
    real(r8), dimension(:), pointer :: FilterShape => NULL()
    character(len=MaxSigLen) :: Signal
  end type FilterShape_T

  ! The filter shape database:
  type(filterShape_T), dimension(:), pointer, public :: FilterShapes => NULL()

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter, private :: ModuleName = &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains

  ! ------------------------------------  Open_Filter_Shapes_File  -----
  subroutine Open_Filter_Shapes_File ( Filename, Lun )

    character(len=*), intent(in) :: Filename ! Name of the filter shape file
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
      & "Unable to open filter shapes file " // Filename )
  end subroutine Open_Filter_Shapes_File

  ! ------------------------------------  Read_Filter_Shapes_File  -----
  subroutine Read_Filter_Shapes_File ( Lun, Spec_Indices )
    use Machine, only: IO_Error
    use Parse_Signal_m, only: Parse_Signal
    use Toggles, only: Gen, Levels, Switches, Toggle
    use Trace_M, only: Trace_begin, Trace_end

    integer, intent(in) :: Lun               ! Logical unit number to read it
    integer, intent(in) :: Spec_Indices(:)   ! Needed by Parse_Signal, q.v.

    integer :: DataBaseSize                  ! How many filter shapes?
    integer :: I                             ! Loop inductors, subscripts
    integer :: NumChannels                   ! For the signal
    integer :: NumFilterPts                  ! How many points in each filter
    !                                          shape array -- all the same
    !                                          for each signal.
    character(len=MaxSigLen) :: SigName      ! Signal Name
    integer :: Status                        ! From read or allocate
    integer, pointer, dimension(:) :: Signal_Indices   ! From Parse_Signal, q.v.
    type(filterShape_T), dimension(:), pointer :: TempFilterShapes => NULL()

    if ( toggle(gen) ) call trace_begin ( "Read_Filter_Shapes_File" )

    if ( associated(filterShapes) ) call destroy_filter_shapes_database

    do
      read ( lun, *, iostat=status ) numFilterPts, sigName
      if ( status > 0 ) go to 99
      if ( status < 0 ) exit
      call parse_signal ( sigName, signal_indices, spec_indices )
      if ( .not. associated(signal_indices) ) &
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & trim(sigName) // " is not a valid signal." )
      numChannels = size(signals(signal_indices(1))%frequencies)
      do i = 1, size(signal_indices)
        if ( size(signals(signal_indices(i))%frequencies) /= numChannels ) &
          call MLSMessage ( MLSMSG_Error, moduleName, &
            & "The signals implied by " // trim(sigName) // &
            & " do not all have the same number of frequencies" )
      end do
      dataBaseSize = 0
      if ( associated(filterShapes) ) dataBaseSize = size(filterShapes)
      tempFilterShapes => filterShapes
      allocate ( filterShapes(dataBaseSize + numChannels), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Allocate // 'FilterShapes' )
      if ( dataBaseSize > 0 ) then
        filterShapes(:dataBaseSize) = tempFilterShapes
        deallocate ( tempFilterShapes, stat=status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
          & MLSMSG_DeAllocate // 'TempFilterShapes' )
      end if
      do i = dataBaseSize + 1, dataBaseSize + numChannels
        filterShapes(i)%signal = sigName
        call allocate_test ( filterShapes(i)%filterShape, numFilterPts, &
          & "filterShapes(i)%filterShape", moduleName )
        call read_one_filter ( filterShapes(i)%lhs, filterShapes(i)%rhs, &
          filterShapes(i)%filterShape )
        if ( status < 0 ) go to 99
        if ( status > 0 ) go to 98
      end do
    end do

    if ( toggle(gen) ) then
      if ( levels(gen) > 0 .or. index(switches,'F') /= 0 ) &
        & call dump_filter_shapes_database
      call trace_end ( "Read_Filter_Shapes_File" )
    end if

    return
  98 call MLSMessage ( MLSMSG_Error, moduleName, "Unexpected end-of-file" )
  99 call io_error ( "While reading the filter shape file", status )
     call MLSMessage ( MLSMSG_Error, moduleName, "Input error" )

  contains
    ! ..........................................  Read_One_Filter  .....
    subroutine Read_One_Filter ( LHS, RHS, ArgFilterShape )
      real(r8), intent(out) :: LHS, RHS, ArgFilterShape(:)
      integer :: Channel           ! Only for its documentary value in the file
      real(r8) :: FilterShape(255)
      integer :: N
      namelist / Filter / Channel, FilterShape, LHS, RHS
      n = size(argFilterShape)
      read ( lun, filter, iostat=status )
      argFilterShape(:n) = filterShape(:n)
    end subroutine Read_One_Filter
  end subroutine Read_Filter_Shapes_File

  ! -----------------------------------  Close_Filter_Shapes_File  -----
  subroutine Close_Filter_Shapes_File ( Lun )
    close ( lun )
  end subroutine Close_Filter_Shapes_File

  ! -----------------------------  Destroy_Filter_Shapes_Database  -----
  subroutine Destroy_Filter_Shapes_Database
    integer :: I, Status
    do i = 1, size(filterShapes)
      call deallocate_test ( filterShapes(i)%filterShape, &
        & "FilterShapes(?)%filterShape", moduleName )
    end do ! i
    deallocate ( filterShapes, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      MLSMSG_DeAllocate // "FilterShapes" )
  end subroutine Destroy_Filter_Shapes_Database

  ! --------------------------------  Dump_Filter_Shapes_Database  -----
  subroutine Dump_Filter_Shapes_Database
    use Dump_0, only: Dump

    integer :: I                   ! Subscripts, loop inductors
    call output ( 'Filter Shapes: SIZE = ' )
    call output ( size(filterShapes), advance='yes' )
    do i = 1, size(filterShapes)
      call output ( i, 4 )
      call output ( ':    Signal =' )
      call output ( trim(filterShapes(i)%signal), advance='yes' )
      call output ( ' LHS = ' )
      call output ( filterShapes(i)%lhs )
      call output ( '  RHS = ' )
      call output ( filterShapes(i)%rhs, advance='yes' )
      call dump ( filterShapes(i)%filterShape, name='FilterShape' )
    end do ! i
  end subroutine Dump_Filter_Shapes_Database

end module FilterShapes_m

! $Log$
! Revision 1.2  2001/03/29 21:57:31  vsnyder
! NAG actually compiles this one
!
! Revision 1.1  2001/03/29 21:42:41  vsnyder
! Changed name to FilterShapes_m
!
! Revision 1.3  2001/03/29 21:27:05  vsnyder
! This one may actually work...
!
! Revision 1.2  2001/03/29 19:49:18  vsnyder
! At least it compiles...
!
