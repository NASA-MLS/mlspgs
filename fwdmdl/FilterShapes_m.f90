! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module FilterShapes_m

  ! Read the filter shapes file.

  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error
  use MLSSignals_m, only: GetNameOfSignal, MaxSigLen, Signals, Signal_T
  
  implicit none

  ! More USEs below in each procedure, if they're only used therein.

  private
  ! Public procedures:
  public :: Open_Filter_Shapes_File
  public :: Read_DACS_Filter_Shapes_File, Read_Filter_Shapes_File
  public :: Close_Filter_Shapes_File
  public :: Destroy_DACS_Filter_Database, Destroy_Filter_Shapes_Database
  public :: Dump_DACS_Filter_Database, Dump_Filter_Shapes_Database

  public :: Dump
  interface Dump
    module procedure Dump_DACS_Filter
  end interface Dump

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

  type, public :: DACSFilterShape_T
    integer :: LogApod, LogFilter, LogNorm
    real(r8), dimension(:), pointer :: FilterGrid => NULL()  ! Abscissae, 2**logFilter+1
    real(r8), dimension(:), pointer :: FilterShape => NULL() ! Ordinates, ditto
    real(r8), dimension(:), pointer :: LO_Apod(:) => NULL()  ! Apodization, ditto
    real(r8), dimension(:), pointer :: CH_Norm(:) => NULL()  ! Normalization, 2**logNorm+1
    type (Signal_T) :: Signal
  end type DACSFilterShape_T

  type(DACSFilterShape_T), dimension(:), pointer, public, save :: &
    & DACSFilterShapes => NULL()

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter, private :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains

  ! ------------------------------------  Open_Filter_Shapes_File  -----
  subroutine Open_Filter_Shapes_File ( Filename, Lun )

    use IO_stuff, only: Get_Lun

    character(len=*), intent(in) :: Filename ! Name of the filter shape file
    integer, intent(out) :: Lun              ! Logical unit number to read it

    logical :: Exist, Opened
    integer :: Status

    call get_lun ( lun, msg=.false. )
    if ( lun < 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "No logical unit numbers available" )
    open ( unit=lun, file=filename, status='old', form='formatted', &
      & access='sequential', iostat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "Unable to open filter shapes file " // Filename )
  end subroutine Open_Filter_Shapes_File

  ! ------------------------------------  Read_Filter_Shapes_File  -----
  subroutine Read_Filter_Shapes_File ( Lun )
    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use Machine, only: IO_Error
    use Parse_Signal_m, only: Parse_Signal
    use Toggles, only: Gen, Switches, Toggle
    use Trace_M, only: Trace_begin, Trace_end

    integer, intent(in) :: Lun          ! Logical unit number to read it

    real(r8) :: DX                      ! To compute FilterGrid
    real(r8) :: LHS, RHS                ! For computing grid
    integer :: I, N                     ! Loop inductor, subscript
    character(80) :: Line               ! From the file
    integer :: Number_in_shape          ! How many points in each filter?
    integer :: NumFilterShapes          ! How many filter shapes in file?
    integer :: Offset                   ! From start of FilterShapes -- to extend it
    integer :: Sideband                 ! From parse signal
                                        ! shape array -- all the same
                                        ! for each signal.
    character(len=MaxSigLen) :: SigName ! Signal Name
    integer :: Status                   ! From read or allocate
    logical, pointer, dimension(:) :: Channels ! Result of parse signal
    integer, pointer, dimension(:) :: Signal_Indices   ! From Parse_Signal, q.v.
    type(filterShape_T), dimension(:), pointer :: TempFilterShapes

    namelist /Filter/ lhs, rhs, number_in_shape

    if ( toggle(gen) ) call trace_begin ( "Read_Filter_Shapes_File" )

    ! Determine the size of the created or expanded FilterShapes array
    offset = 0
    if ( associated(filterShapes) ) offset = size(filterShapes)
    numFilterShapes = offset
    do ! Loop over filter shapes
      call read_a_line ( lun, line, status ) ! Read the signal
      if ( status < 0 ) exit
      if ( status > 0 ) go to 99
      ! Read the lhs, rhs and number_in_shape
      read ( lun, filter, iostat=status )
      if ( status /= 0 ) go to 99
      ! Skip the shape, using LHS for a temp
      read ( lun, *, iostat=status ) (lhs, i = 1, number_in_shape)
      if ( status /= 0 ) go to 99
      numFilterShapes = numFilterShapes + 1
    end do
    rewind ( lun )

    ! Create or expand the FilterShapes array
    tempFilterShapes => filterShapes
    allocate ( filterShapes(numFilterShapes), stat=status )
    if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "FilterShapes")
    if ( associated(tempFilterShapes) ) then
      filterShapes(:offset) = tempFilterShapes
      deallocate ( tempFilterShapes, stat=status )
      if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate // "TempFilterShapes")
    end if

    ! Read and store the filter shapes
    n = offset
    do ! Loop over filter shapes
      call read_a_line ( lun, line, status ) ! Read the signal
      if ( status < 0 ) exit
      if ( status > 0 ) go to 99
      n = n + 1
      sigName = line
      nullify ( channels, signal_indices ) 
      call parse_signal ( sigName, signal_indices, sideband=sideband, &
        channels=channels )

      if ( .not. associated(signal_indices) ) &
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & trim(sigName) // " is not a valid signal." )
      ! Just take the first one.
      filterShapes(n)%signal = signals(signal_indices(1))
      filterShapes(n)%signal%sideband = sideband
      filterShapes(n)%signal%channels => channels
      call deallocate_test ( signal_indices, "Signal_Indices", moduleName )

      ! Read the lhs, rhs and number_in_shape
      read ( lun, filter, iostat=status )
      if ( status /= 0 ) go to 99

      call allocate_test ( filterShapes(n)%filterGrid,&
        & number_in_shape, 'filterShapes(n)%filterGrid', ModuleName )
      call allocate_test ( filterShapes(n)%filterShape,&
        & number_in_shape, 'filterShapes(n)%filterShape', ModuleName )

      ! Read the shape array and calculate the associated abscissae
      dx = ( rhs - lhs ) / (number_in_shape - 1)
      do i = 1, number_in_shape
        filterShapes(n)%filterGrid(i) = lhs + (i-1) * dx
      end do ! i
      read ( lun, *, iostat=status ) filterShapes(n)%filterShape
      if ( status /= 0 ) go to 99

    end do ! Loop over filter shapes

    if ( index(switches,'filt') /= 0 ) call dump_filter_shapes_database
    if ( toggle(gen) ) then
      call trace_end ( "Read_Filter_Shapes_File" )
    end if

    return
 98 call MLSMessage ( MLSMSG_Error, moduleName, "Unexpected end-of-file" )
 99 call io_error ( "While reading the filter shape file", status )
    call MLSMessage ( MLSMSG_Error, moduleName, "Input error" )

  end subroutine Read_Filter_Shapes_File

  ! -------------------------------  Read_DACS_Filter_Shapes_File  -----
  subroutine Read_DACS_Filter_Shapes_File ( Lun )
    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use Machine, only: IO_Error
    use Parse_Signal_m, only: Parse_Signal
    use Toggles, only: Gen, Switches, Toggle
    use Trace_M, only: Trace_begin, Trace_end

    integer, intent(in) :: Lun          ! Logical unit number to read it

    real(r8) :: DX                      ! To compute FilterGrid
    real(r8) :: LHS, RHS                ! For computing grid
    integer :: I, N                     ! Loop inductor, subscript
    integer :: LogApod                  ! Number_apod == 2 ** LogApod + 1
    integer :: LogFilt                  ! Number_shape == 2 ** LogSize + 1
    integer :: LogNorm                  ! Number_norm == 2 ** LogNorm + 1
                                        ! Default 7
    character(80) :: Line               ! From the file
    integer :: Number_apod              ! How many LO apodizations?
    integer :: Number_norm              ! How many channel normalizations?
    integer :: Number_shape             ! How many points in each filter?
    integer :: NumFilterShapes          ! How many filter shapes in file?
    integer :: Offset                   ! From start of FilterShapes -- to extend it
    integer :: Sideband                 ! From parse signal
                                        ! shape array -- all the same
                                        ! for each signal.
    character(len=MaxSigLen) :: SigName ! Signal Name
    integer :: Status                   ! From read or allocate
    logical, pointer, dimension(:) :: Channels ! Result of parse signal
    integer, pointer, dimension(:) :: Signal_Indices   ! From Parse_Signal, q.v.
    type(DACSfilterShape_T), dimension(:), pointer :: TempFilterShapes

    namelist /Filter/ lhs, rhs, logApod, logFilt, logNorm

    if ( toggle(gen) ) call trace_begin ( "Read_DACS_Filter_Shapes_File" )

    ! Determine the size of the created or expanded DACSfilterShapes array
    offset = 0
    if ( associated(DACSfilterShapes) ) offset = size(DACSfilterShapes)
    numFilterShapes = offset
    do ! Loop over filter shapes
      call read_a_line ( lun, line, status ) ! Read the signal
      if ( status < 0 ) exit
      if ( status > 0 ) go to 99
      ! Read the lhs, rhs and number_shape
      logApod = -1
      logFilt = -1
      logNorm = 7
      read ( lun, filter, iostat=status )
      if ( status /= 0 ) go to 99
      if ( logFilt < 0 ) &
        & call MLSMessage ( MLSMSG_Error, moduleName, "LogFilt < 0" )
      if ( logFilt /= logApod .and. logNorm /= logApod ) &
        & call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Require LogApod == LogFilt or LogApod == LogNorm" )
      ! Skip the shape, using LHS for a temp
      number_shape = 2 ** logFilt + 1
      read ( lun, *, iostat=status ) (lhs, i = 1, number_shape)
      if ( status /= 0 ) go to 99
      ! Skip the LO apodization
      number_apod = 2 ** logApod + 1
      read ( lun, *, iostat=status ) (lhs, i = 1, number_apod)
      if ( status /= 0 ) go to 99
      ! Skip the channel normalization
      number_norm = 2 ** logNorm + 1
      read ( lun, *, iostat=status ) (lhs, i = 1, number_norm)
      if ( status /= 0 ) go to 99
      numFilterShapes = numFilterShapes + 1
    end do
    rewind ( lun )

    ! Create or expand the DACSfilterShapes array
    tempFilterShapes => DACSfilterShapes
    allocate ( DACSfilterShapes(numFilterShapes), stat=status )
    if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "DACSfilterShapes")
    if ( associated(tempFilterShapes) ) then
      DACSfilterShapes(:offset) = tempFilterShapes
      deallocate ( tempFilterShapes, stat=status )
      if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate // "TempFilterShapes")
    end if

    ! Read and store the filter shapes.  Don't need error checks since
    ! if we get here we've been able to read the file once successfully
    n = offset
    do ! Loop over filter shapes
      call read_a_line ( lun, line, status ) ! Read the signal
      if ( status < 0 ) exit
      n = n + 1
      sigName = line
      nullify ( channels, signal_indices ) 
      call parse_signal ( sigName, signal_indices, sideband=sideband, &
        channels=channels )

      if ( .not. associated(signal_indices) ) &
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & trim(sigName) // " is not a valid signal." )
      ! Just take the first one.
      DACSfilterShapes(n)%signal = signals(signal_indices(1))
      DACSfilterShapes(n)%signal%sideband = sideband
      DACSfilterShapes(n)%signal%channels => channels
      call deallocate_test ( signal_indices, "Signal_Indices", moduleName )

      ! Read the lhs, rhs and number_shape
      logNorm = 7
      read ( lun, filter )
      DACSfilterShapes(n)%logApod = logApod
      DACSfilterShapes(n)%logFilter = logFilt
      DACSfilterShapes(n)%logNorm = logNorm
      number_shape = 2 ** logFilt + 1

      call allocate_test ( DACSfilterShapes(n)%filterGrid,&
        & number_shape, 'DACSfilterShapes(n)%filterGrid', ModuleName )
      call allocate_test ( DACSfilterShapes(n)%filterShape,&
        & number_shape, 'DACSfilterShapes(n)%filterShape', ModuleName )

      ! Read the shape array and calculate the associated abscissae
      dx = ( rhs - lhs ) / (number_shape - 1)
      do i = 1, number_shape
        DACSfilterShapes(n)%filterGrid(i) = lhs + (i-1) * dx
      end do ! i
      read ( lun, * ) DACSfilterShapes(n)%filterShape

      ! Read the LO_Apod and CH_Norm arrays
      number_apod = 2 ** logApod + 1
      call allocate_test ( DACSfilterShapes(n)%lo_apod,&
        & number_apod, 'DACSfilterShapes(n)%lo_apod', ModuleName )
      read ( lun, *, iostat=status ) DACSfilterShapes(n)%lo_apod
      if ( status /= 0 ) go to 99

      number_norm = 2 ** logNorm + 1
      call allocate_test ( DACSfilterShapes(n)%ch_norm,&
        & number_norm, 'DACSfilterShapes(n)%ch_norm', ModuleName )
      read ( lun, *, iostat=status ) DACSfilterShapes(n)%ch_norm
      if ( status /= 0 ) go to 99

    end do ! Loop over filter shapes

    if ( index(switches,'dacsfil') /= 0 ) call dump_DACS_filter_database
    if ( toggle(gen) ) then
      call trace_end ( "Read_DACS_Filter_Shapes_File" )
    end if

    return
 98 call MLSMessage ( MLSMSG_Error, moduleName, "Unexpected end-of-file" )
 99 call io_error ( "While reading the DACS filter shape file", status )
    call MLSMessage ( MLSMSG_Error, moduleName, "Input error" )

  end subroutine Read_DACS_Filter_Shapes_File

  ! -----------------------------------  Close_Filter_Shapes_File  -----
  subroutine Close_Filter_Shapes_File ( Lun )
    integer, intent(in) :: lun
    close ( lun )
  end subroutine Close_Filter_Shapes_File

  ! -----------------------------------  AddFilterShapeToDatabase  -----
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

  ! -------------------------------  AddDACSFilterShapeToDatabase  -----
  integer function AddDACSFilterShapeToDatabase ( database, item )

    ! Add a quantity template to a database, or create the database if it
    ! doesn't yet exist
    
    ! Dummy arguments
    type (DACSFilterShape_T), dimension(:), pointer :: database
    type (DACSFilterShape_T), intent(in) :: item

    ! Local variables
    type (DACSFilterShape_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddDACSFilterShapeToDatabase = newSize
  end function AddDACSFilterShapeToDatabase

  ! -----------------------------  Destroy_Filter_Shapes_Database  -----
  subroutine Destroy_Filter_Shapes_Database
    use Allocate_Deallocate, only: Deallocate_Test
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

  ! ------------------------  Destroy_DACS_Filter_Database  -----
  subroutine Destroy_DACS_Filter_Database
    use Allocate_Deallocate, only: Deallocate_Test
    integer :: I, Status
    if ( .not. associated(DACSfilterShapes) ) return
    do i = 1, size(filterShapes)
      call deallocate_test ( DACSfilterShapes(i)%filterGrid, &
        & "DACSFilterShapes(?)%filterGrid", moduleName )
      call deallocate_test ( DACSfilterShapes(i)%filterShape, &
        & "DACSFilterShapes(?)%filterShape", moduleName )
      call deallocate_test ( DACSfilterShapes(i)%signal%channels, &
        & "DACSFilterShapes(?)%signal%channels", moduleName )
      call deallocate_test ( DACSfilterShapes(i)%lo_apod, &
        & "DACSFilterShapes(?)%lo_apod", moduleName )
      call deallocate_test ( DACSfilterShapes(i)%ch_norm, &
        & "DACSFilterShapes(?)%ch_norm", moduleName )
    end do ! i
    deallocate ( DACSfilterShapes, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      MLSMSG_DeAllocate // "DACSFilterShapes" )
  end subroutine Destroy_DACS_Filter_Database

  ! --------------------------------  Dump_Filter_Shapes_Database  -----
  subroutine Dump_Filter_Shapes_Database ( where )
    use Dump_0, only: Dump
    use MoreTree, only: StartErrorMessage
    use Output_m, only: Output

    integer, intent(in), optional :: Where   ! Tree node index

    integer :: I                   ! Subscripts, loop inductors
    character(len=MaxSigLen) :: sigName
    if ( associated(filterShapes) ) then
      call output ( 'Filter Shapes: SIZE = ' )
      call output ( size(filterShapes), advance='yes' )
      do i = 1, size(filterShapes)
        call output ( i, 4 )
        call output ( ':    Signal = ' )
        call GetNameOfSignal ( filterShapes(i)%signal, sigName )
        call output ( sigName, advance='yes' )
        call dump ( filterShapes(i)%filterShape, name='FilterShape' )
        call dump ( filterShapes(i)%filterGrid, name='FilterGrid', &
          & width=4, format='(1x,1pg18.11)' )
      end do ! i
    else
      if ( present(where) ) call startErrorMessage ( where )
      call output ( 'No filter shapes database to dump.', advance='yes' )
    end if
  end subroutine Dump_Filter_Shapes_Database

  ! -------------------------------------------  Dump_DACS_Filter  -----
  subroutine Dump_DACS_Filter ( Filter )
    use Dump_0, only: Dump
    use Output_m, only: Output
    type (DACSFilterShape_T), intent(in) :: Filter
    character(len=MaxSigLen) :: sigName

    call output ( 'Signal = ' )
    call GetNameOfSignal ( filter%signal, sigName )
    call output ( trim(sigName), advance='yes' )
    call dump ( filter%filterShape, name='FilterShape' )
    call dump ( filter%filterGrid, name='FilterGrid', width=5, format='(f14.5)' )
    call dump ( filter%lo_apod, name='LO_Apod', width=7, format='(f10.6)' )
    call dump ( filter%ch_norm, name='CH_Norm' )
  end subroutine Dump_DACS_Filter

  ! ----------------------------------  Dump_DACS_Filter_Database  -----
  subroutine Dump_DACS_Filter_Database ( where )
    use MoreTree, only: StartErrorMessage
    use Output_m, only: Output

    integer, intent(in), optional :: Where   ! Tree node index

    integer :: I                   ! Subscripts, loop inductors
    if ( associated(DACSfilterShapes) ) then
      call output ( 'DACS Filter Shapes: SIZE = ' )
      call output ( size(DACSfilterShapes), advance='yes' )
      do i = 1, size(DACSfilterShapes)
        call output ( i, 4 )
        call output ( ':    ' )
        call dump_DACS_filter ( DACSfilterShapes(i) )
      end do ! i
    else
      if ( present(where) ) call startErrorMessage ( where )
      call output ( 'No DACS filter shapes database to dump.', advance='yes' )
    end if
  end subroutine Dump_DACS_Filter_Database

  ! ------------------------------------------------  Read_A_Line  -----
  subroutine Read_A_Line ( Lun, Line, Status )
  ! Read a line, skipping blank lines and lines that begin with !
    use Machine, only: IO_Error
    integer, intent(in) :: Lun
    character(len=*), intent(out) :: Line
    integer, intent(out) :: Status

    do
      read ( lun, '(a)', iostat=status ) line
      if ( status < 0 ) return
      if ( status > 0 ) go to 99
      if ( line == '' ) cycle           ! Skip blank lines
      line = adjustl(line)
      if ( line(1:1) /= '!' ) return
      ! Skip comments
    end do
 99 call io_error ( "While reading the filter shape file", status )
    call MLSMessage ( MLSMSG_Error, moduleName, "Input error" )
  end subroutine Read_A_Line

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module FilterShapes_m

! $Log$
! Revision 2.16  2004/05/22 02:27:50  vsnyder
! Use Get_Lun from io_stuff instead of repeating it here
!
! Revision 2.15  2004/03/30 00:45:31  vsnyder
! Remove USE for unreferenced symbol
!
! Revision 2.14  2004/02/14 00:23:48  vsnyder
! New DACS convolution algorithm
!
! Revision 2.13  2004/02/09 20:20:32  vsnyder
! Make an error message more precise and informative
!
! Revision 2.12  2004/02/06 00:45:27  vsnyder
! Allow more general relation between sizes of filter shape and
! apodization arrays.
!
! Revision 2.11  2004/01/21 22:01:37  vsnyder
! Trim signal name before printing, cosmetics
!
! Revision 2.10  2003/08/15 00:20:50  michael
! changed debugging switch for dump_DACS_filter_database
!
! Revision 2.9  2003/07/16 01:06:50  vsnyder
! Add DACS filter shapes
!
! Revision 2.8  2003/05/19 19:58:07  vsnyder
! Remove USEs for unreferenced symbols, remove unused local variables
!
! Revision 2.7  2003/05/10 22:20:57  livesey
! Tried to calm down -g1..
!
! Revision 2.6  2003/05/05 23:00:24  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.5.2.1  2003/04/21 20:07:08  vsnyder
! Count filter shapes, then allocate the right size
!
! Revision 2.5  2002/10/08 17:08:03  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.4  2002/09/13 22:06:22  vsnyder
! Move a few USEs from module scope to procedure scope
!
! Revision 2.3  2002/05/14 20:02:12  livesey
! Bug fix in handling of channels field in signal part of filter shape.
!
! Revision 2.2  2002/05/10 00:33:18  vsnyder
! Repair a botched comment
!
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
