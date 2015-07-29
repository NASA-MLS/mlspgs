! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Checkpoint

  ! Write and read some stuff on Fortran files instead of HDF

  use MLSKinds, only: R8

  implicit none

  private

  public :: CheckpointCommand, CloseCheckpointFile
  public :: MakeCheckpoint, OpenCheckpointFile, Restore
  public :: ReadVector, WriteVector

  type :: TemplateStuff_t ! Used to read, write and test a quantity template
    integer :: QuantityType
    integer :: NoInstances
    integer :: NoSurfs
    integer :: NoChans
    logical :: Coherent
    logical :: Stacked
    logical :: Regular
    logical :: MinorFrame
    logical :: MajorFrame
    logical :: LogBasis
    real(r8) :: MinValue
    integer :: noInstancesLowerOverlap
    integer :: noInstancesUpperOverlap
    real(r8) :: badValue
    integer :: unit
    integer :: instanceLen
    integer :: verticalCoordinate
    logical :: sharedVGrid
    integer :: vGridIndex
    integer :: hGridIndex              
    integer :: instanceOffset          
    integer :: grandTotalInstances     
    integer :: frequencyCoordinate
    logical :: sharedFGrid
    integer :: fGridIndex
    real(r8) :: lo
    integer :: signal
    integer :: sideband
    integer :: instrumentModule
    integer :: radiometer
    integer :: reflector
    integer :: molecule
    logical :: AssociatedSurfs
    logical :: AssociatedPhi
    logical :: AssociatedGeodLat
    logical :: AssociatedLon
    logical :: AssociatedTime
    logical :: AssociatedSolarTime
    logical :: AssociatedSolarZenith
    logical :: AssociatedLosAngle
    logical :: AssociatedFrequencies
    logical :: AssociatedSurfIndex
    logical :: AssociatedChanIndex
    integer :: BoundsSurfs(2,2)
    integer :: BoundsPhi(2,2)
    integer :: BoundsGeodLat(2,2)
    integer :: BoundsLon(2,2)
    integer :: BoundsTime(2,2)
    integer :: BoundsSolarTime(2,2)
    integer :: BoundsSolarZenith(2,2)
    integer :: BoundsLosAngle(2,2)
    integer :: BoundsFrequencies(2)
    integer :: BoundsSurfIndex(2,2)
    integer :: BoundsChanIndex(2,2)
  end type

  ! Used to read, write and test a vector
  type VectorValueStuff_t
    logical :: AssociatedValues
    logical :: AssociatedMask
    integer :: Bounds(2,2)
  end type VectorValueStuff_t

  ! From the FileName field of the most recent Checkpoint command
  character(len=255), save :: CheckpointFileName = ''

  ! Root of the most recent Checkpoint command, if any.
  integer, save :: CheckpointRoot = 0

  ! Root of the Vectors field, if any, from the most recent Checkpoint
  ! command, if any.
  integer, save :: CheckpointVectorsRoot = 0

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine CheckpointCommand ( Root )
    ! Save the root index of the Checkpoint command.
    ! Get the file name from the Checkpoint command.
    ! The vectors will be gotten later, when and if they're needed.
    use Init_Tables_Module, only: F_Filename, F_Vectors
    use MoreTree, only: Get_Field_ID
    use String_Table, only: Get_String
    use Tree, only: NSons, Sub_Rosa, Subtree
    integer, intent(in) :: Root
    integer :: I, Son
    checkpointRoot = root
    do i = 1, nsons(root)
      son = subtree(i,root)
      if ( get_field_id(son) == f_fileName ) then
        call get_string ( sub_rosa(subtree(2,son)), checkpointFileName )
      else if ( get_field_id(son) == f_vectors ) then
        checkpointVectorsRoot = son
      end if
    end do
  end subroutine CheckpointCommand

  logical function CheckQuantityTemplate ( Unit, QuantityTemplate )
    ! Read a quantity template from Unit.
    ! Return true if it appears to be the same as QuantityTemplate, else false.

    use Allocate_Deallocate, only: Test_Allocate
    use QuantityTemplates, only: QuantityTemplate_T

    integer, intent(in) :: Unit ! Fortran I/O unit, already open
    type (quantityTemplate_t), intent(in) :: QuantityTemplate

    integer :: Stat
    type (templateStuff_t) :: TemplateStuff

    real(r8), allocatable :: Surfs(:,:)
    real(r8), allocatable :: Phi(:,:)
    real(r8), allocatable :: GeodLat(:,:)
    real(r8), allocatable :: Lon(:,:)
    real(r8), allocatable :: Time(:,:)
    real(r8), allocatable :: SolarTime(:,:)
    real(r8), allocatable :: SolarZenith(:,:)
    real(r8), allocatable :: LosAngle(:,:)
    real(r8), allocatable :: Frequencies(:)
    integer, allocatable :: SurfIndex(:,:)
    integer, allocatable :: ChanIndex(:,:)

    checkQuantityTemplate = .false. ! Assume for now it won't be OK

    read ( unit, err=9, end=9 ) templateStuff

    if ( &
      templateStuff%QuantityType            /= QuantityTemplate%QuantityType     .or. &
      templateStuff%NoInstances             /= QuantityTemplate%NoInstances      .or. &
      templateStuff%NoSurfs                 /= QuantityTemplate%NoSurfs          .or. &
      templateStuff%NoChans                 /= QuantityTemplate%NoChans          .or. &
      templateStuff%Coherent            .neqv. QuantityTemplate%Coherent         .or. &
      templateStuff%Stacked             .neqv. QuantityTemplate%Stacked          .or. &
      templateStuff%Regular             .neqv. QuantityTemplate%Regular          .or. &
      templateStuff%MinorFrame          .neqv. QuantityTemplate%MinorFrame       .or. &
      templateStuff%MajorFrame          .neqv. QuantityTemplate%MajorFrame       .or. &
      templateStuff%LogBasis            .neqv. QuantityTemplate%LogBasis         .or. &
      templateStuff%MinValue                /= QuantityTemplate%MinValue         .or. &
      templateStuff%noInstancesLowerOverlap /= QuantityTemplate%noInstancesLowerOverlap .or. &
      templateStuff%noInstancesUpperOverlap /= QuantityTemplate%noInstancesUpperOverlap .or. &
      templateStuff%badValue                /= QuantityTemplate%badValue         .or. &
      templateStuff%unit                    /= QuantityTemplate%unit             .or. &
      templateStuff%instanceLen             /= QuantityTemplate%instanceLen      .or. &
      templateStuff%verticalCoordinate      /= QuantityTemplate%verticalCoordinate .or. &
      templateStuff%sharedVGrid         .neqv. QuantityTemplate%sharedVGrid      .or. &
      templateStuff%vGridIndex              /= QuantityTemplate%vGridIndex       .or. &
      templateStuff%hGridIndex              /= QuantityTemplate%hGridIndex       .or. &
      templateStuff%instanceOffset          /= QuantityTemplate%instanceOffset   .or. &
      templateStuff%grandTotalInstances     /= QuantityTemplate%grandTotalInstances .or. &
      templateStuff%frequencyCoordinate     /= QuantityTemplate%frequencyCoordinate .or. &
      templateStuff%sharedFGrid         .neqv. QuantityTemplate%sharedFGrid      .or. &
      templateStuff%fGridIndex              /= QuantityTemplate%fGridIndex       .or. &
      templateStuff%lo                      /= QuantityTemplate%lo               .or. &
      templateStuff%signal                  /= QuantityTemplate%signal           .or. &
      templateStuff%sideband                /= QuantityTemplate%sideband         .or. &
      templateStuff%instrumentModule        /= QuantityTemplate%instrumentModule .or. &
      templateStuff%radiometer              /= QuantityTemplate%radiometer       .or. &
      templateStuff%reflector               /= QuantityTemplate%reflector        .or. &
      templateStuff%molecule                /= QuantityTemplate%molecule ) return

    if ( &
      templateStuff%AssociatedSurfs       .neqv. allocated(QuantityTemplate%surfs) .or. &
      templateStuff%AssociatedPhi         .neqv. allocated(QuantityTemplate%phi) .or. &
      templateStuff%AssociatedGeodLat     .neqv. allocated(QuantityTemplate%geodLat) .or. &
      templateStuff%AssociatedLon         .neqv. allocated(QuantityTemplate%lon) .or. &
      templateStuff%AssociatedTime        .neqv. associated(QuantityTemplate%time) .or. &
      templateStuff%AssociatedSolarTime   .neqv. associated(QuantityTemplate%solarTime) .or. &
      templateStuff%AssociatedSolarZenith .neqv. associated(QuantityTemplate%solarZenith) .or. &
      templateStuff%AssociatedLosAngle    .neqv. associated(QuantityTemplate%losAngle) .or. &
      templateStuff%AssociatedFrequencies .neqv. associated(QuantityTemplate%frequencies) .or. &
      templateStuff%AssociatedSurfIndex   .neqv. associated(QuantityTemplate%surfIndex) .or. &
      templateStuff%AssociatedChanIndex   .neqv. associated(QuantityTemplate%chanIndex) ) return

    if ( templateStuff%AssociatedSurfs       ) then
      if ( any(templateStuff%boundsSurfs(:,1) /= lbound(QuantityTemplate%surfs)) .or. &
           any(templateStuff%boundsSurfs(:,2) /= ubound(QuantityTemplate%surfs)) ) return
      allocate ( surfs(templateStuff%boundsSurfs(1,1):templateStuff%boundsSurfs(1,2), &
                     & templateStuff%boundsSurfs(2,1):templateStuff%boundsSurfs(2,2)), stat=stat )
      call test_allocate ( stat, moduleName, 'Surfs', templateStuff%boundsSurfs(:,1), templateStuff%boundsSurfs(:,2) )
      read ( unit, err=9, end=9 ) surfs
      if ( any(surfs /= QuantityTemplate%surfs) ) return
    end if
    if ( templateStuff%AssociatedPhi         ) then
      if ( any(templateStuff%boundsPhi(:,1) /= lbound(QuantityTemplate%phi)) .or. &
           any(templateStuff%boundsPhi(:,2) /= ubound(QuantityTemplate%phi)) ) return
      allocate ( phi(templateStuff%boundsPhi(1,1):templateStuff%boundsPhi(1,2), &
                   & templateStuff%boundsPhi(2,1):templateStuff%boundsPhi(2,2)), stat=stat )
      call test_allocate ( stat, moduleName, 'phi', templateStuff%boundsPhi(:,1), templateStuff%boundsPhi(:,2) )
      read ( unit, err=9, end=9 ) phi
      if ( any(phi /= QuantityTemplate%phi) ) return
    end if
    if ( templateStuff%AssociatedGeodLat     ) then
      if ( any(templateStuff%boundsGeodLat(:,1) /= lbound(QuantityTemplate%geodLat)) .or. &
           any(templateStuff%boundsGeodLat(:,2) /= ubound(QuantityTemplate%geodLat)) ) return
      allocate ( geodLat(templateStuff%boundsGeodLat(1,1):templateStuff%boundsGeodLat(1,2), &
                       & templateStuff%boundsGeodLat(2,1):templateStuff%boundsGeodLat(2,2)), stat=stat )
      call test_allocate ( stat, moduleName, 'geodLat', templateStuff%boundsGeodLat(:,1), templateStuff%boundsGeodLat(:,2) )
      read ( unit, err=9, end=9 ) geodLat
      if ( any(geodLat /= QuantityTemplate%geodLat) ) return
    end if
    if ( templateStuff%AssociatedLon         ) then
      if ( any(templateStuff%boundsLon(:,1) /= lbound(QuantityTemplate%lon)) .or. &
           any(templateStuff%boundsLon(:,2) /= ubound(QuantityTemplate%lon)) ) return
      allocate ( lon(templateStuff%boundsLon(1,1):templateStuff%boundsLon(1,2), &
                   & templateStuff%boundsLon(2,1):templateStuff%boundsLon(2,2)), stat=stat )
      call test_allocate ( stat, moduleName, 'lon', templateStuff%boundsLon(:,1), templateStuff%boundsLon(:,2) )
      read ( unit, err=9, end=9 ) lon
      if ( any(lon /= QuantityTemplate%lon) ) return
    end if
    if ( templateStuff%AssociatedTime        ) then
      if ( any(templateStuff%boundsTime(:,1) /= lbound(QuantityTemplate%time)) .or. &
           any(templateStuff%boundsTime(:,2) /= ubound(QuantityTemplate%time)) ) return
      allocate ( time(templateStuff%boundsTime(1,1):templateStuff%boundsTime(1,2), & 
                    & templateStuff%boundsTime(2,1):templateStuff%boundsTime(2,2)), stat=stat )
      call test_allocate ( stat, moduleName, 'time', templateStuff%boundsTime(:,1), templateStuff%boundsTime(:,2) )
      read ( unit, err=9, end=9 ) time
      if ( any(time /= QuantityTemplate%time) ) return
    end if
    if ( templateStuff%AssociatedSolarTime   ) then
      if ( any(templateStuff%boundsSolarTime(:,1) /= lbound(QuantityTemplate%solarTime)) .or. &
           any(templateStuff%boundsSolarTime(:,2) /= ubound(QuantityTemplate%solarTime)) ) return
      allocate ( solarTime(templateStuff%boundsSolarTime(1,1):templateStuff%boundsSolarTime(1,2), &
                         & templateStuff%boundsSolarTime(2,1):templateStuff%boundsSolarTime(2,2)), stat=stat )
      call test_allocate ( stat, moduleName, 'solarTime', templateStuff%boundsSolarTime(:,1), templateStuff%boundsSolarTime(:,2) )
      read ( unit, err=9, end=9 ) solarTime
      if ( any(solarTime /= QuantityTemplate%solarTime) ) return
    end if
    if ( templateStuff%AssociatedSolarZenith ) then
       if ( any(templateStuff%boundsSolarZenith(:,1) /= lbound(QuantityTemplate%solarZenith)) .or. &
            any(templateStuff%boundsSolarZenith(:,2) /= ubound(QuantityTemplate%solarZenith)) ) return
      allocate ( solarZenith(templateStuff%boundsSolarZenith(1,1):templateStuff%boundsSolarZenith(1,2), &
                           & templateStuff%boundsSolarZenith(2,1):templateStuff%boundsSolarZenith(2,2)), stat=stat )
      call test_allocate ( stat, moduleName, 'solarZenith', templateStuff%boundsSolarZenith(:,1), &
                         & templateStuff%boundsSolarZenith(:,2) )
      read ( unit, err=9, end=9 ) solarZenith
      if ( any(solarZenith /= QuantityTemplate%solarZenith) ) return
    end if
    if ( templateStuff%AssociatedLosAngle    ) then
       if ( any(templateStuff%boundsLOSAngle(:,1) /= lbound(QuantityTemplate%LOSAngle)) .or. &
            any(templateStuff%boundsLOSAngle(:,2) /= ubound(QuantityTemplate%LOSAngle)) ) return
      allocate ( LOSAngle(templateStuff%boundsLOSAngle(1,1):templateStuff%boundsLOSAngle(1,2), &
                        & templateStuff%boundsLOSAngle(2,1):templateStuff%boundsLOSAngle(2,2)), stat=stat )
      call test_allocate ( stat, moduleName, 'LOSAngle', templateStuff%boundsLOSAngle(:,1), templateStuff%boundsLOSAngle(:,2) )
      read ( unit, err=9, end=9 ) LOSAngle
      if ( any(LOSAngle /= QuantityTemplate%LOSAngle) ) return
    end if
    if ( templateStuff%AssociatedFrequencies ) then
       if ( templateStuff%boundsFrequencies(1) /= lbound(QuantityTemplate%frequencies,1) .or. &
            templateStuff%boundsFrequencies(2) /= ubound(QuantityTemplate%frequencies,1) ) return
      allocate ( frequencies(templateStuff%boundsFrequencies(1):templateStuff%boundsFrequencies(2)),  stat=stat )
      call test_allocate ( stat, moduleName, 'frequencies', templateStuff%boundsFrequencies(1:1), &
                         & templateStuff%boundsFrequencies(2:2) )
      read ( unit, err=9, end=9 ) frequencies
      if ( any(frequencies /= QuantityTemplate%frequencies) ) return
    end if
    if ( templateStuff%AssociatedSurfIndex   ) then
       if ( any(templateStuff%boundsSurfIndex(:,1) /= lbound(QuantityTemplate%surfIndex)) .or. &
            any(templateStuff%boundsSurfIndex(:,2) /= ubound(QuantityTemplate%surfIndex)) ) return
      allocate ( surfIndex(templateStuff%boundsSurfIndex(1,1):templateStuff%boundsSurfIndex(1,2), &
                         & templateStuff%boundsSurfIndex(2,1):templateStuff%boundsSurfIndex(2,2)), stat=stat )
      call test_allocate ( stat, moduleName, 'surfIndex', templateStuff%boundsSurfIndex(:,1), templateStuff%boundsSurfIndex(:,2) )
      read ( unit, err=9, end=9 ) surfIndex
      if ( any(surfIndex /= QuantityTemplate%surfIndex) ) return
    end if
    if ( templateStuff%AssociatedChanIndex   ) then
       if ( any(templateStuff%boundsChanIndex(:,1) /= lbound(QuantityTemplate%chanIndex)) .or. &
            any(templateStuff%boundsChanIndex(:,2) /= ubound(QuantityTemplate%chanIndex)) ) return
      allocate ( chanIndex(templateStuff%boundsChanIndex(1,1):templateStuff%boundsChanIndex(1,2), &
                         & templateStuff%boundsChanIndex(2,1):templateStuff%boundsChanIndex(2,2)), stat=stat )
      call test_allocate ( stat, moduleName, 'chanIndex', templateStuff%boundsChanIndex(:,1), templateStuff%boundsChanIndex(:,2) )
      read ( unit, err=9, end=9 ) chanIndex
      if ( any(chanIndex /= QuantityTemplate%chanIndex) ) return
    end if

    checkQuantityTemplate = .true.
    ! Return deallocates all allocatable local variables

  9 return
  end function CheckQuantityTemplate

  subroutine CloseCheckpointFile ( Unit )
    integer, intent(in) :: Unit ! Fortran I/O unit, assumed to be open
    close ( unit )
  end subroutine CloseCheckpointFile

  subroutine MakeCheckpoint ( VectorDatabase, Name, Vectors )
    ! Write a checkpoint either of Vectors on Name, or using
    ! the specifications of a previous Checkpoint command
    use Machine, only: IO_Error
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use Tree, only: Decoration, NSons, Subtree
    use VectorsModule, only: Vector_t

    type (vector_t), intent(in) :: VectorDatabase(:)
    character(len=*), intent(in), optional :: Name ! of file to use, else
      ! use the one mentioned on a prior Checkpoint command, if any.
    type (vector_t), intent(in), optional :: Vectors(:) ! to save, else use
      ! the ones mentioned on a prior Checkpoint command, if any.

    integer :: I, Status, Unit

    if ( checkpointRoot <= 0 ) then
      if ( .not. present(name) .or. .not. present(vectors) ) return
    end if
    if ( present(vectors) ) then
      call openCheckpointFile ( unit, "UNKNOWN", name )
      do i = 1, size(vectors)
        call writeVector ( unit, vectors(i), status )
        if ( status /= 0 ) go to 9
      end do
    else
      if ( checkpointVectorsRoot <= 0 ) return
      call openCheckpointFile ( unit, "UNKNOWN", name )
      do i = 2, nsons(checkpointVectorsRoot)
        call writeVector ( unit, &
          & vectorDatabase(decoration(subtree(i,checkpointVectorsRoot))), &
          & status )
        if ( status /= 0 ) go to 9
      end do ! i
    end if

    call closeCheckpointFile ( unit )
    return
  9 continue
    if ( present(name) ) then
      call io_error ( "Trying to write checkpoint", status, name )
    else
      call io_error ( "Trying to write checkpoint", status, checkpointFileName )
    end if
    call MLSMessage ( MLSMSG_Error, moduleName, &
      & "Failed to write checkpoint" )

  end subroutine MakeCheckpoint

  subroutine OpenCheckpointFile ( Unit, Status, Name )
    ! Get a unit number.  Open it on the file NAME if NAME is present, else
    ! open it on the saved checkpoint file name
    use IO_Stuff, only: Get_LUN
    use Machine, only: IO_Error
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    integer, intent(out) :: Unit ! Fortran I/O unit
    character(len=*), intent(in) :: Status ! "OLD" or "UNKNOWN"
    character(len=*), intent(in), optional :: Name
    integer :: Stat
    call get_lun ( unit )
    if ( present(name) ) then
      open ( unit, file=name, form='unformatted', status=status, iostat=stat )
      if ( stat /= 0 ) &
        & call io_error ( 'Unable to open checkpoint file', stat, name )
    else
      open ( unit, file=checkpointFileName, form='unformatted', iostat=stat )
      if ( stat /= 0 ) &
        & call io_error ( 'Unable to open checkpoint file', stat, checkpointFileName )
    end if
    if ( stat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & 'Unable to open checkpoint file' )
  end subroutine OpenCheckpointFile

  subroutine ReadVector ( Unit, Vector, Status )
    ! Read a vector from Unit, while verifying its properties are the same
    ! as the one on Unit
    use VectorsModule, only: CreateMask, CreateVectorValue, Vector_T
    integer, intent(in) :: Unit    ! Fortran I/O unit, already open
    type (vector_t), intent(inout) :: Vector
    integer, intent(out) :: Status ! 0 = OK

    integer :: I
    integer :: NQ ! Number of vector quantities
    type (vectorValueStuff_t) :: VectorValueStuff

    status = 1 ! Assume failure

    read ( unit, err=9, end=9 ) nq
    if ( nq /= size(vector%quantities) ) return

    do i = 1, nq
      if ( .not. checkQuantityTemplate(unit,vector%quantities(i)%template) ) return
      read ( unit, err=9, end=9 ) vectorValueStuff
      if ( vectorValueStuff%associatedValues ) then
        call createVectorValue ( vector%quantities(i), what="ReadVector's vector%quantities(i)" )
        read ( unit, err=9, end=9 ) vector%quantities(i)%values
      end if
      if ( vectorValueStuff%associatedMask ) then
        call createMask ( vector%quantities(i), forWhom="ReadVector's vector%quantities(i)" )
        read ( unit, err=9, end=9 ) vector%quantities(i)%mask
      end if
    end do

    status = 0
  9 return
  end subroutine ReadVector

  subroutine Restore ( VectorDatabase, Name, Vectors )
    ! Restore a checkpoint either of Vectors from Name, or using
    ! the specifications of a previous Checkpoint command
    use Machine, only: IO_Error
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use Tree, only: Decoration, NSons, Subtree
    use VectorsModule, only: Vector_t

    type (vector_t), intent(inout) :: VectorDatabase(:)
    character(len=*), intent(in), optional :: Name ! of file to use, else
      ! use the one mentioned on a prior Checkpoint command, if any.
    type (vector_t), intent(inout), optional :: Vectors(:) ! to save, else use
      ! the ones mentioned on a prior Checkpoint command, if any.

    integer :: I, Status, Unit

    if ( checkpointRoot <= 0 ) then
      if ( .not. present(name) .or. .not. present(vectors) ) return
    end if
    if ( present(vectors) ) then
      call openCheckpointFile ( unit, "UNKNOWN", name )
      do i = 1, size(vectors)
        call readVector ( unit, vectors(i), status )
        if ( status /= 0 ) go to 9
      end do
    else
      if ( checkpointVectorsRoot <= 0 ) return
      call openCheckpointFile ( unit, "UNKNOWN", name )
      do i = 2, nsons(checkpointVectorsRoot)
        call readVector ( unit, &
          & vectorDatabase(decoration(subtree(i,checkpointVectorsRoot))), &
          & status )
        if ( status /= 0 ) go to 9
      end do ! i
    end if

    call closeCheckpointFile ( unit )
    return
  9 continue
    if ( present(name) ) then
      call io_error ( "Trying to write checkpoint", status, name )
    else
      call io_error ( "Trying to write checkpoint", status, checkpointFileName )
    end if
    call MLSMessage ( MLSMSG_Error, moduleName, &
      & "Failed to write checkpoint" )

  end subroutine Restore

  subroutine WriteVector ( Unit, Vector, Status )
    ! Write a vector on Unit so that ReadVector can read it
    use VectorsModule, only: Vector_T
    integer, intent(in) :: Unit    ! Fortran I/O unit, already open
    type (vector_t), intent(in) :: Vector
    integer, intent(out) :: Status ! = 0 OK

    integer :: I
    type (templateStuff_t) :: TemplateStuff
    type (vectorValueStuff_t) :: VectorValueStuff

    write ( unit, err=9, iostat=status ) size(vector%quantities)

    do i = 1, size(vector%quantities)
      templateStuff%QuantityType            = vector%quantities(i)%template%QuantityType
      templateStuff%NoInstances             = vector%quantities(i)%template%NoInstances
      templateStuff%NoSurfs                 = vector%quantities(i)%template%NoSurfs
      templateStuff%NoChans                 = vector%quantities(i)%template%NoChans
      templateStuff%Coherent                = vector%quantities(i)%template%Coherent
      templateStuff%Stacked                 = vector%quantities(i)%template%Stacked
      templateStuff%Regular                 = vector%quantities(i)%template%Regular
      templateStuff%MinorFrame              = vector%quantities(i)%template%MinorFrame
      templateStuff%MajorFrame              = vector%quantities(i)%template%MajorFrame
      templateStuff%LogBasis                = vector%quantities(i)%template%LogBasis
      templateStuff%MinValue                = vector%quantities(i)%template%MinValue
      templateStuff%noInstancesLowerOverlap = vector%quantities(i)%template%noInstancesLowerOverlap
      templateStuff%noInstancesUpperOverlap = vector%quantities(i)%template%noInstancesUpperOverlap
      templateStuff%badValue                = vector%quantities(i)%template%badValue
      templateStuff%unit                    = vector%quantities(i)%template%unit
      templateStuff%instanceLen             = vector%quantities(i)%template%instanceLen
      templateStuff%verticalCoordinate      = vector%quantities(i)%template%verticalCoordinate
      templateStuff%sharedVGrid             = vector%quantities(i)%template%sharedVGrid
      templateStuff%vGridIndex              = vector%quantities(i)%template%vGridIndex
      templateStuff%hGridIndex              = vector%quantities(i)%template%hGridIndex
      templateStuff%instanceOffset          = vector%quantities(i)%template%instanceOffset
      templateStuff%grandTotalInstances     = vector%quantities(i)%template%grandTotalInstances
      templateStuff%frequencyCoordinate     = vector%quantities(i)%template%frequencyCoordinate
      templateStuff%sharedFGrid             = vector%quantities(i)%template%sharedFGrid
      templateStuff%fGridIndex              = vector%quantities(i)%template%fGridIndex
      templateStuff%lo                      = vector%quantities(i)%template%lo
      templateStuff%signal                  = vector%quantities(i)%template%signal
      templateStuff%sideband                = vector%quantities(i)%template%sideband
      templateStuff%instrumentModule        = vector%quantities(i)%template%instrumentModule
      templateStuff%radiometer              = vector%quantities(i)%template%radiometer
      templateStuff%reflector               = vector%quantities(i)%template%reflector
      templateStuff%molecule                = vector%quantities(i)%template%molecule
      templateStuff%AssociatedSurfs         = allocated(vector%quantities(i)%template%surfs)
      templateStuff%AssociatedPhi           = allocated(vector%quantities(i)%template%phi)
      templateStuff%AssociatedGeodLat       = allocated(vector%quantities(i)%template%geodLat)
      templateStuff%AssociatedLon           = allocated(vector%quantities(i)%template%lon)
      templateStuff%AssociatedTime          = associated(vector%quantities(i)%template%time)
      templateStuff%AssociatedSolarTime     = associated(vector%quantities(i)%template%solarTime)
      templateStuff%AssociatedSolarZenith   = associated(vector%quantities(i)%template%solarZenith)
      templateStuff%AssociatedLosAngle      = associated(vector%quantities(i)%template%losAngle)
      templateStuff%AssociatedFrequencies   = associated(vector%quantities(i)%template%frequencies)
      templateStuff%AssociatedSurfIndex     = associated(vector%quantities(i)%template%surfIndex)
      templateStuff%AssociatedChanIndex     = associated(vector%quantities(i)%template%chanIndex)
      templateStuff%boundsSurfs = 0
      templateStuff%boundsPhi = 0
      templateStuff%boundsGeodLat = 0
      templateStuff%boundsLon = 0
      templateStuff%boundsTime = 0
      templateStuff%boundsSolarTime = 0
      templateStuff%boundsSolarZenith = 0
      templateStuff%boundsLosAngle = 0
      templateStuff%boundsFrequencies = 0
      templateStuff%boundsSurfIndex = 0
      templateStuff%boundsChanIndex = 0
      if ( templateStuff%AssociatedSurfs       ) then
        templateStuff%boundsSurfs(:,1) = lbound(vector%quantities(i)%template%surfs)
        templateStuff%boundsSurfs(:,2) = ubound(vector%quantities(i)%template%surfs)
      end if
      if ( templateStuff%AssociatedPhi         ) then
        templateStuff%boundsPhi(:,1) = lbound(vector%quantities(i)%template%phi)
        templateStuff%boundsPhi(:,2) = ubound(vector%quantities(i)%template%phi)
      end if
      if ( templateStuff%AssociatedGeodLat     ) then
        templateStuff%boundsGeodLat(:,1) = lbound(vector%quantities(i)%template%GeodLat)
      end if
      if ( templateStuff%AssociatedLon         ) then
        templateStuff%boundsLon(:,1) = lbound(vector%quantities(i)%template%Lon)
        templateStuff%boundsLon(:,2) = ubound(vector%quantities(i)%template%Lon)
      end if
      if( templateStuff%AssociatedTime        ) then
        templateStuff%boundsTime(:,1) = lbound(vector%quantities(i)%template%Time)
        templateStuff%boundsTime(:,2) = ubound(vector%quantities(i)%template%Time)
      end if
      if( templateStuff%AssociatedSolarTime   ) then
        templateStuff%boundsSolarTime(:,1) = lbound(vector%quantities(i)%template%SolarTime)
        templateStuff%boundsSolarTime(:,2) = ubound(vector%quantities(i)%template%SolarTime)
      end if
      if( templateStuff%AssociatedSolarZenith ) then
        templateStuff%boundsSolarZenith(:,1) = lbound(vector%quantities(i)%template%SolarZenith)
        templateStuff%boundsSolarZenith(:,2) = ubound(vector%quantities(i)%template%SolarZenith)
      end if
      if( templateStuff%AssociatedLosAngle    ) then
        templateStuff%boundsLosAngle(:,1) = lbound(vector%quantities(i)%template%LosAngle)
        templateStuff%boundsLosAngle(:,2) = ubound(vector%quantities(i)%template%LosAngle)
      end if
      if( templateStuff%AssociatedFrequencies ) then
        templateStuff%boundsFrequencies(1) = lbound(vector%quantities(i)%template%Frequencies,1)
        templateStuff%boundsFrequencies(2) = ubound(vector%quantities(i)%template%Frequencies,1)
      end if
      if ( templateStuff%AssociatedSurfIndex   ) then
        templateStuff%boundsSurfIndex(:,1) = lbound(vector%quantities(i)%template%SurfIndex)
        templateStuff%boundsSurfIndex(:,2) = ubound(vector%quantities(i)%template%SurfIndex)
      end if
      if ( templateStuff%AssociatedChanIndex   ) then
        templateStuff%boundsChanIndex(:,1) = lbound(vector%quantities(i)%template%ChanIndex)
        templateStuff%boundsChanIndex(:,2) = ubound(vector%quantities(i)%template%ChanIndex)
      end if
      write ( unit, err=9, iostat=status ) templateStuff

      if ( templateStuff%AssociatedSurfs       ) then
        write ( unit, err=9, iostat=status ) vector%quantities(i)%template%surfs
      end if
      if ( templateStuff%AssociatedPhi         ) then
        write ( unit, err=9, iostat=status ) vector%quantities(i)%template%Phi
      end if
      if ( templateStuff%AssociatedGeodLat     ) then
        write ( unit, err=9, iostat=status ) vector%quantities(i)%template%GeodLat
      end if
      if ( templateStuff%AssociatedLon         ) then
        write ( unit, err=9, iostat=status ) vector%quantities(i)%template%Lon
      end if
      if ( templateStuff%AssociatedTime        ) then
        write ( unit, err=9, iostat=status ) vector%quantities(i)%template%Time
      end if
      if ( templateStuff%AssociatedSolarTime   ) then
        write ( unit, err=9, iostat=status ) vector%quantities(i)%template%SolarTime
      end if
      if ( templateStuff%AssociatedSolarZenith ) then
        write ( unit, err=9, iostat=status ) vector%quantities(i)%template%SolarZenith
      end if
      if ( templateStuff%AssociatedLosAngle    ) then
        write ( unit, err=9, iostat=status ) vector%quantities(i)%template%LosAngle
      end if
      if ( templateStuff%AssociatedFrequencies ) then
        write ( unit, err=9, iostat=status ) vector%quantities(i)%template%Frequencies
      end if
      if ( templateStuff%AssociatedSurfIndex   ) then
        write ( unit, err=9, iostat=status ) vector%quantities(i)%template%SurfIndex
      end if
      if ( templateStuff%AssociatedChanIndex   ) then
        write ( unit, err=9, iostat=status ) vector%quantities(i)%template%ChanIndex
      end if

      vectorValueStuff%associatedValues = associated(vector%quantities(i)%values)
      vectorValueStuff%associatedMask = associated(vector%quantities(i)%mask)
      vectorValueStuff%bounds = 0
      if ( vectorValueStuff%associatedValues ) then
        vectorValueStuff%bounds(:,1) = lbound(vector%quantities(i)%values)
        vectorValueStuff%bounds(:,2) = ubound(vector%quantities(i)%values)
      end if
      write ( unit, err=9, iostat=status ) vectorValueStuff
      if ( vectorValueStuff%associatedValues ) &
        & write ( unit, err=9, iostat=status ) vector%quantities(i)%values
      if ( vectorValueStuff%associatedMask ) &
        & write ( unit, err=9, iostat=status ) vector%quantities(i)%mask

    end do

  9 return
  end subroutine WriteVector

!=======================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Checkpoint

! $Log$
! Revision 2.8  2015/07/29 00:29:54  vsnyder
! Convert Phi from pointer to allocatable
!
! Revision 2.7  2015/06/04 03:14:14  vsnyder
! Make Surfs component of quantity template allocatable
!
! Revision 2.6  2015/05/28 18:24:11  vsnyder
! Remove shared HGrid
!
! Revision 2.5  2014/01/11 01:44:18  vsnyder
! Decruftification
!
! Revision 2.4  2012/07/31 00:47:00  vsnyder
! Use CreateVectorValue and CreateMask abstractions in ReadVector
!
! Revision 2.3  2010/02/04 23:12:44  vsnyder
! Remove USE or declaration for unreferenced names
!
! Revision 2.2  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.1  2007/10/02 22:49:18  vsnyder
! Initial commit.  Might be useful someday
!
