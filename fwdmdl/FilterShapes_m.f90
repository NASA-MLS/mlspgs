! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module FilterShapes_m

  ! Read the filter shapes file.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error
  use MLSSignals_m, only: GetNameOfSignal, MaxSigLen, Signals, Signal_T
  use Dump_0, only: Dump
  
  implicit none

  ! More USEs below in each procedure, if they're only used therein.

  private
  ! Public procedures:
  public :: Open_Filter_Shapes_File, Read_Filter_Shapes_File
  public :: Close_Filter_Shapes_File, Destroy_Filter_Shapes_Database
  public :: Dump_Filter_Shapes_Database

  ! There is a separate FilterShape_T object for each COMPLETE signal
  ! specification, including the channel number.  It isn't necessary
  ! for all of the filter shapes to have the same size.
  type, public :: FilterShape_T
    real(r8), dimension(:), pointer :: FilterGrid => NULL()      ! Abscissae
    real(r8), dimension(:), pointer :: FilterShape => NULL()     ! Ordinates
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

    real(r8) :: DX                      ! To compute FilterGrid
    real(r8) :: LHS, RHS                ! For computing grid
    integer :: I                        ! Loop inductor, subscript
    character(80) :: Line               ! From the file
    integer :: Number_in_shape          ! How many points in each filter
    integer :: Sideband                 ! From parse signal
    !                                          shape array -- all the same
    !                                          for each signal.
    character(len=MaxSigLen) :: SigName      ! Signal Name
    integer :: Status                        ! From read or allocate
    integer :: dummy                    ! Result of add to database
    integer, pointer, dimension(:) :: Signal_Indices   ! From Parse_Signal, q.v.
    type(filterShape_T) :: thisShape

    namelist /Filter/ lhs, rhs, number_in_shape

    if ( toggle(gen) ) call trace_begin ( "Read_Filter_Shapes_File" )

    nullify ( signal_indices )
    ! if ( associated(filterShapes) ) call destroy_filter_shapes_database

    do ! Loop over filter shapes
      read ( lun, '(a)', iostat=status ) line
      if ( status > 0 ) go to 99
      if ( status < 0 ) exit
      if ( line == '' ) cycle           ! Skip blank lines
      line = adjustl(line)
      if ( line(1:1) == '!' ) cycle     ! Skip comments
      sigName = line
      nullify ( thisShape%signal%channels )
      call parse_signal ( sigName, signal_indices, sideband=sideband, &
        channels=thisShape%signal%channels )
      if ( .not. associated(signal_indices) ) &
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & trim(sigName) // " is not a valid signal." )
      ! Just take the first one.
      thisShape%signal = signals(signal_indices(1))
      thisShape%signal%sideband = sideband
      call deallocate_test ( signal_indices, "Signal_Indices", moduleName )

      ! Read the lhs, rhs and num_in_shape
      read ( lun, filter, iostat=status )
      if ( status /= 0 ) go to 99

      ! Need to nullify so this add doesn't hose any previous work
      nullify ( thisShape%filterShape, thisShape%filterGrid )
      call allocate_test ( thisShape%filterGrid,&
        & number_in_shape, 'thisShape%filterGrid', ModuleName )
      call allocate_test ( thisShape%filterShape,&
        & number_in_shape, 'thisShape%filterShape', ModuleName )

      ! Read the shape array and calculate the associate abscissae
      read ( lun, *, iostat=status ) thisShape%filterShape
      if ( status /= 0 ) go to 99
      dx = ( rhs - lhs ) / (number_in_shape - 1)
      do i = 1, number_in_shape
        thisShape%filterGrid(i) = lhs + (i-1) * dx
      end do ! i

      dummy = AddFilterShapeToDatabase ( filterShapes, thisShape )

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
      call deallocate_test ( filterShapes(i)%signal%channels, &
        & "FilterShapes(?)%signal%channels", moduleName )
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
      call GetNameOfSignal ( filterShapes(i)%signal, sigName )
      call output ( sigName, advance='yes' )
      call dump ( filterShapes(i)%filterShape, name='FilterShape' )
      call dump ( filterShapes(i)%filterGrid, name='FilterGrid' )
    end do ! i
  end subroutine Dump_Filter_Shapes_Database

end module FilterShapes_m

! $Log$
! Revision 2.1  2002/05/10 00:21:39  vsnyder
! Revise to cope with new filter shape file.  filterShapes%filterGrid
! and filterShapes%filterShape are now one-dimensional.
! filterShapes%signal%channel has exactly one .true. element that
! indicates the channel to which the shape applies.
!
! Revision 2.0  2001/09/17 20:26:25  livesey
! New forward model
!
! Revision 1.15  2001/05/17 01:00:22  livesey
! Odd bug, was allowing me to call AddFilterShapeToDatabase as a
! subroutine when it was in fact a function.
!
! Revision 1.14  2001/05/16 23:04:29  livesey
! Bug fix.
!
! Revision 1.13  2001/05/16 01:25:02  livesey
! New version.  Stores thing differently.
!
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
