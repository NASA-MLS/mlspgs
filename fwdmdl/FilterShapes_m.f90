! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module FilterShapes_m

  ! Read the filter shapes file.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error
  use MLSSignals_m, only: GetSignalName, MaxSigLen, Signals, Signal_T
  
  implicit none

  ! More USEs below in each procedure, if they're only used therein.

  private
  ! Public procedures:
  public :: Open_Filter_Shapes_File, Read_Filter_Shapes_File
  public :: Close_Filter_Shapes_File, Destroy_Filter_Shapes_Database
  public :: Dump_Filter_Shapes_Database

  type, public :: FilterShape_T
    real(r8), dimension(:,:), pointer :: FilterGrid => NULL()      ! Abscissae
    real(r8), dimension(:,:), pointer :: FilterShape => NULL()     ! Ordinates
    type (Signal_T) :: Signal
  end type FilterShape_T

  ! The filter shape database:
  type(filterShape_T), dimension(:), pointer, public, save :: &
    & FilterShapes => NULL()

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
  subroutine Read_Filter_Shapes_File ( Lun )
    use Machine, only: IO_Error
    use Parse_Signal_m, only: Parse_Signal
    use Toggles, only: Gen, Levels, Switches, Toggle
    use Trace_M, only: Trace_begin, Trace_end

    integer, intent(in) :: Lun          ! Logical unit number to read it

    integer :: DataBaseSize             ! How many filter shapes?
    real(r8) :: DX                      ! To compute FilterGrid
    real(r8) :: LHS, RHS                ! For computing grid
    integer :: I, J                     ! Loop inductors, subscripts
    integer :: NumChannels              ! For the signal
    integer :: NumFilterPts             ! How many points in each filter
    integer :: Sideband                 ! From parse signal
    !                                          shape array -- all the same
    !                                          for each signal.
    character(len=MaxSigLen) :: SigName      ! Signal Name
    integer :: Status                        ! From read or allocate
    integer, pointer, dimension(:) :: Signal_Indices   ! From Parse_Signal, q.v.
    type(filterShape_T) :: thisShape

    if ( toggle(gen) ) call trace_begin ( "Read_Filter_Shapes_File" )

    nullify ( signal_indices )
    ! if ( associated(filterShapes) ) call destroy_filter_shapes_database

    do ! Loop over filter shapes
      read ( lun, *, iostat=status ) numFilterPts, sigName
      if ( status > 0 ) go to 99
      if ( status < 0 ) exit
      call parse_signal ( sigName, signal_indices, sideband=sideband )
      if ( .not. associated(signal_indices) ) &
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & trim(sigName) // " is not a valid signal." )
      if ( size ( signal_indices ) /= 1 ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & trim(sigName) // " is ambiguous." )

      thisShape%signal = signals(signal_indices(1))
      thisShape%signal%sideband = sideband
      numChannels = size(thisShape%signal%frequencies)
      call deallocate_test ( signal_indices, "Signal_Indices", moduleName )
      
      call allocate_test ( thisShape%filterGrid,&
        & numChannels, numFilterPts, &
        & 'thisShape%filterGrid', ModuleName )
      call allocate_test ( thisShape%filterShape,&
        & numChannels, numFilterPts, &
        & 'thisShape%filterShape', ModuleName )

      do i = 1, numChannels
        call read_one_filter ( lhs, rhs, thisShape%filterShape(i,:) )
        if ( status < 0 ) go to 99
        if ( status > 0 ) go to 98
        dx = ( rhs - lhs ) / (numFilterPts - 1)
        do j = 1, numFilterPts
          thisShape%filterGrid(i,j) = lhs + (j-1) * dx
        end do ! j
      end do ! i

      call AddFilterShapeToDatabase ( filterShapes, thisShape )
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
    integer, intent(in) :: lun
    close ( lun )
  end subroutine Close_Filter_Shapes_File

  ! ----------------------------------- AddFilterShapeToDatabase -------
  integer function AddFilterShapeToDatabase ( database, item )

    ! Add a quantity template to a database, or create the database if it
    ! doesn't yet exist
    
    ! Dummy arguments
    type (FilterShape_T), dimension(:), pointer :: database
    type (FilterShape_T), intent(in) :: item

    ! Local variables
    type (FilterShape_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddFilterShapeToDatabase = newSize
  end function AddFilterShapeToDatabase

  ! -----------------------------  Destroy_Filter_Shapes_Database  -----
  subroutine Destroy_Filter_Shapes_Database
    integer :: I, Status
    if ( .not. associated(filterShapes) ) return
    do i = 1, size(filterShapes)
      call deallocate_test ( filterShapes(i)%filterGrid, &
        & "FilterShapes(?)%filterGrid", moduleName )
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
    use Output_m, only: Output

    integer :: I                   ! Subscripts, loop inductors
    character(len=MaxSigLen) :: sigName
    call output ( 'Filter Shapes: SIZE = ' )
    call output ( size(filterShapes), advance='yes' )
    do i = 1, size(filterShapes)
      call output ( i, 4 )
      call output ( ':    Signal = ' )
      call GetSignalName ( filterShapes(i)%signal%index, sigName )
      call output ( sigName, advance='yes' )
      call dump ( filterShapes(i)%filterShape, name='FilterShape' )
      call dump ( filterShapes(i)%filterGrid, name='FilterGrid' )
    end do ! i
  end subroutine Dump_Filter_Shapes_Database

end module FilterShapes_m

! $Log$
! Revision 1.12  2001/05/04 00:49:06  livesey
! Let destroy quit if nothing to destroy
!
! Revision 1.11  2001/05/03 22:05:22  vsnyder
! Add a nullify, make database SAVE,
!
! Revision 1.10  2001/04/26 02:36:52  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 1.9  2001/04/21 01:21:11  vsnyder
! Fix a memory leak
!
! Revision 1.8  2001/04/20 17:19:05  vsnyder
! Deallocate FilterGrid component in Destroy...
!
! Revision 1.7  2001/04/02 20:56:56  vsnyder
! Add FilterGrid field and compute it
!
! Revision 1.6  2001/03/30 02:10:12  vsnyder
! Improve 'dump' routine
!
! Revision 1.5  2001/03/30 01:12:29  vsnyder
! Correct some comments, move "use Output" to "dump_filter_shapes_database"
!
! Revision 1.4  2001/03/30 00:02:06  livesey
! Nullified another pointer
!
! Revision 1.3  2001/03/29 23:53:06  vsnyder
! This one seems to work (not just compile)
!
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
