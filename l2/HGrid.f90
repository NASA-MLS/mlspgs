! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module HGrid                    ! Horizontal grid information
!=============================================================================

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use EXPR_M, only: EXPR
  use Dump_0, only: DUMP
  use INIT_TABLES_MODULE, only: F_FRACTION, F_HEIGHT, F_INCLINATION, &
    & F_INTERPOLATIONFACTOR, F_MODULE, F_TYPE, F_VALUES, FIELD_FIRST,&
    & F_SOURCEL2GP, FIELD_LAST, L_EXPLICIT, L_FRACTIONAL, L_HEIGHT, L_L2GP,&
    & LIT_INDICES, PHYQ_DIMENSIONLESS, PHYQ_LENGTH
  use LEXER_CORE, only: PRINT_SOURCE
  use L1BData, only: DeallocateL1BData, L1BData_T, ReadL1BData
  use MLSCommon, only: L1BInfo_T, MLSChunk_T, NameLen, R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_allocate, &
    & MLSMSG_DeAllocate, MLSMSG_Error, MLSMSG_Info, MLSMSG_L1BRead
  use MLSNumerics, only: HUNT
  use OUTPUT_M, only: OUTPUT
  use STRING_TABLE, only: GET_STRING
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TOGGLES, only: GEN, TOGGLE
  use TREE, only: DECORATION, DUMP_TREE_NODE, NSONS, NULL_TREE, SOURCE_REF, &
                  SUB_ROSA, SUBTREE
  use UNITS, only: DEG2RAD, RAD2DEG
  use L2GPData, only: L2GPDATA_T

  implicit none

  private :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
  character (len=256) :: Id = &
       "$Id$"
  character (len=*), parameter :: ModuleName= "$RCSfile$"
  ! ---------------------------------------------------------------------------

  ! This module contains datatypes and routines for handling HGrid information
  ! HGrids are the horizontal gridding information that get into vector
  ! quantities.

  ! This is the main datatype, an HGrid.

  type HGrid_T
    integer :: NAME                 ! String index of name.
    integer :: noProfs              ! Number of profiles in this grid
    integer :: noProfsLowerOverlap  ! Number of profiles in the lower overlap
    integer :: noProfsUpperOverlap  ! Number of profiles in the upper overlap

    ! Now the various coordinates in the HGrid, all dimensioned (noProfs)
    real(r8), dimension(:), pointer :: phi => NULL()
    real(r8), dimension(:), pointer :: geodLat => NULL()
    real(r8), dimension(:), pointer :: lon => NULL()
    real(r8), dimension(:), pointer :: time => NULL()
    real(r8), dimension(:), pointer :: solarTime => NULL()
    real(r8), dimension(:), pointer :: solarZenith => NULL()
    real(r8), dimension(:), pointer :: losAngle => NULL()
  end type HGrid_T

! -----     Private declarations     ---------------------------------

  integer, private :: ERROR

! Error codes for "announce_error"
  integer, private, parameter :: LengthUnitMessage = 1
  integer, private, parameter :: NoFraction = LengthUnitMessage + 1
  integer, private, parameter :: NoHeight = NoFraction + 1
  integer, private, parameter :: UnitlessMessage = NoHeight + 1

contains ! =====     Public Procedures     =============================

  ! -----------------------------------------  AddHGridToDatabase  -----
  integer function AddHGridToDatabase ( database, item )

    ! Dummy arguments
    type (HGrid_T), dimension(:), pointer :: database
    type (HGrid_T), intent(in) :: item

    ! Local variables
    type (HGrid_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddHGridToDatabase = newSize
  end function AddHGridToDatabase

  ! -----------------------------------  CreateHGridFromMLSCFInfo  -----
  type(hGrid_T) function CreateHGridFromMLSCFInfo &
    & ( name, root, l1bInfo, l2gpDatabase, chunk ) result ( hGrid )

  ! This routine creates an hGrid based on the user requests.

    ! Dummy arguments
    integer, intent(in) :: NAME               ! String index of name
    integer, intent(in) :: ROOT               ! Root of hGrid subtree
    type (L1BInfo_T), intent(in) :: L1BINFO   ! File handles for l1b data
    type (L2GPData_T), pointer, dimension(:) :: L2GPDATABASE
    type (MLSChunk_T), intent(in) :: CHUNK    ! This chunk

    ! Local variables
    integer :: EXPR_UNITS(2)            ! Output from Expr subroutine
    double precision :: EXPR_VALUE(2)   ! Output from Expr subroutine
    integer :: hGridType
    real(r8) :: interpolationFactor
    integer :: instrumentModule
    real(r8) :: fraction, height
    type (L2GPData_T), pointer :: L2GP  ! The l2gp to use 

    integer :: keyNo            ! Entry in the mlscf line

    type (L1BData_T) :: l1bField ! L1B data

    real(r8) :: incline                 ! Orbital inclination / degrees
    real(r8) :: MINTIME, MAXTIME        ! Span for a chunk
    integer, dimension(2) :: PROFRANGE  ! Profile range
    integer :: A,B                      ! Elements of profile range

    integer :: FIELD              ! Subtree index of "field" node
    integer :: FIELD_INDEX        ! F_..., see Init_Tables_Module
    logical :: GOT_FIELD(field_first:field_last)
    integer :: L1BFLAG
    integer :: L1BITEM            ! Loop counter
    integer :: MAF                ! Loop counters
    integer :: NOMAFS             ! Number of MAFs of L1B data read
    integer :: PROF                     ! Loop counter
    integer :: SON                ! Son of Root
    integer :: STATUS             ! From Allocate, ReadL1B... etc.
    integer :: VALUESNODE               ! Node of tree for explicit

    character (len=NameLen) :: InstrumentModuleName

    ! Executable code

    hGrid%name = name
    if ( toggle(gen) ) call trace_begin ( "CreateHGridFromMLSCFInfo", root )

    got_field = .false.
    interpolationFactor = 1.0

    do keyNo = 2, nsons(root)
      son = subtree(keyNo,root)
      field = subtree(1,son)
      field_index = decoration(field)
      ! The tree checker prevents duplicate fields
      got_field(field_index) = .true.
      select case ( field_index )
      case ( f_type )
        hGridType = decoration(subtree(2,son))
      case ( f_module )
        instrumentModule = sub_rosa(subtree(2,son))
        call get_string ( instrumentModule , instrumentModuleName )
      case ( f_height )
        call expr ( subtree(2,son), expr_units, expr_value )
        height = expr_value(1)
        if ( expr_units(1) /= PHYQ_Length) &
          & call announce_error ( field, lengthUnitMessage )
      case ( f_fraction )
        call expr ( subtree(2,son), expr_units, expr_value )
        fraction = expr_value(1)
        if ( expr_units(1) /= PHYQ_Dimensionless) &
          & call announce_error ( field, unitlessMessage )
      case ( f_interpolationFactor )
        call expr ( subtree(2,son), expr_units, expr_value )
        interpolationFactor = expr_value(1)
        if ( expr_units(1) /= PHYQ_Dimensionless) &
          & call announce_error ( field, unitlessMessage )
      case ( f_values )
        valuesNode = son
      case ( f_sourceL2gp )
        l2gp => l2gpDatabase(decoration(decoration(subtree(2,son))))
      case ( f_inclination )
        call expr (subtree ( 2, son), expr_units, expr_value )
        incline = expr_value(1)
      case default ! Can't get here if tree_checker works correctly
      end select
    end do

    select case (hGridType)

    case ( l_height, l_fractional ) ! ----- Fractional or Height ------
      call CreateCommonHGrids ( l1bInfo, hGridType, chunk, &
        & got_field, root, height, fraction, interpolationFactor, &
        & instrumentModuleName, hGrid )

    case ( l_explicit ) ! ----------------- Explicit ------------------
      ! For explicit hGrids, do things differently
      hGrid%noProfs = nsons(valuesNode)-1
      hGrid%noProfsLowerOverlap = 0
      hGrid%noProfsUpperOverlap = 0
      call CreateEmptyHGrid(hGrid)
      hGrid%lon = 0.0_r8
      hGrid%time = 0.0_r8
      hGrid%solarTime = 0.0_r8
      hGrid%losAngle = 0.0_r8
      hGrid%solarZenith = 0.0_r8
      
      do prof = 1, hGrid%noProfs
        call expr ( subtree ( prof+1, valuesNode), expr_units, expr_value )
        hGrid%phi(prof) = expr_value(1)
        hGrid%geodLat(prof) = hGrid%phi(prof) !???? Sort this out later!
      end do

    case ( l_l2gp) ! -------------------- L2GP ------------------------
      
      ! Get the time from the l1b file
      call ReadL1BData ( l1bInfo%l1boaID, "MAFStartTimeTAI", l1bField, noMAFs, &
        & l1bFlag, firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex)
      if ( l1bFlag==-1) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_L1BRead//"MAFStartTimeTAI" )
      
      minTime = l1bField%dpField(1,1,1)
      maxTime = l1bField%dpField(1,1,noMAFs)

      call Hunt ( l2gp%time, (/minTime, maxTime/), profRange, allowTopValue=.true. )
      profRange(1) = min ( profRange(1)+1, l2gp%nTimes )
      
      a = profRange(1)
      b = profRange(2)
      hGrid%noProfs = b - a + 1
      call CreateEmptyHGrid(hGrid)
      hGrid%phi =         l2gp%geodAngle(a:b)
      hGrid%geodLat =     l2gp%latitude(a:b)
      hGrid%lon =         l2gp%longitude(a:b)
      hGrid%time =        l2gp%time(a:b)
      hGrid%solarTime =   l2gp%solarTime(a:b)
      hGrid%solarZenith = l2gp%solarZenith(a:b)
      hGrid%losAngle =    l2gp%losAngle(a:b)

    end select
    
    if ( toggle(gen) ) call trace_end ( "CreateHGridFromMLSCFInfo" )

  end function CreateHGridFromMLSCFInfo

  ! ------------------------------------- CreateCommonHGrids -----------
  subroutine CreateCommonHGrids ( l1bInfo, hGridType, &
    & chunk, got_field, root, height, fraction, interpolationFactor,&
    & instrumentModuleName, hGrid )
    ! This is part of ConstructHGridFromMLSCFInfo
    type (L1BInfo_T), intent(in)      :: L1BINFO
    integer, intent(in)               :: HGRIDTYPE
    type (MLSChunk_T), intent(in)     :: CHUNK
    logical, dimension(:), intent(in) :: GOT_FIELD
    integer, intent(in)               :: ROOT
    real(r8), intent(in)              :: INTERPOLATIONFACTOR
    real(r8), intent(in)              :: HEIGHT
    real(r8), intent(in)              :: FRACTION
    character (len=*)                 :: INSTRUMENTMODULENAME
    type (HGrid_T), intent(inout)     :: HGRID ! Needs inout as name set by caller

    ! Local variables / parameters
    integer, parameter :: L1B_MAFSTARTTIMETAI = 1
    integer, parameter :: L1B_TPGEODLAT       = l1b_MAFStartTimeTAI + 1
    integer, parameter :: L1B_TPLON           = l1b_tpGeodLat + 1
    integer, parameter :: L1B_TPGEODANGLE     = l1b_tpLon + 1
    integer, parameter :: L1B_TPSOLARZENITH   = l1b_tpGeodAngle + 1
    integer, parameter :: L1B_TPSOLARTIME     = l1b_tpSolarZenith + 1
    integer, parameter :: L1B_TPLOSANGLE      = l1b_tpSolarTime + 1
    integer, parameter :: NOL1BITEMSTOREAD=l1b_tpLosAngle
    integer, parameter :: FIRSTMODULARITEM=l1b_tpGeodLat

    real(r8), parameter :: SIXTH = 1.0_r8 / 6.0_r8

    character (len=15), DIMENSION(noL1BItemsToRead) :: l1bItemNames
    ! Entries in the above array follwing FirstModularItem are prefixed
    ! with either GHz or THz. 

    integer :: L1BITEM                  ! Loop counter
    integer :: L1BFLAG                  ! Flag
    integer :: NOMAFS                   ! Dimension
    integer :: STATUS                   ! Flag
    integer :: MAF                      ! Loop counter etc.
    real(r8) :: minAngle, maxAngle
    real(r8), dimension(:,:,:), pointer :: tpGeodAngle, tpGeodAlt

    ! MIFs it would choose in the non over/undersampled case
    real(r8), dimension(:), allocatable :: defaultField, interpolatedField

    type (L1BData_T) :: l1bField ! L1B data
    integer, dimension(:), allocatable :: defaultMIFs
    character (len=NameLen) :: L1BItemName

    ! Executable code
    l1bItemNames(l1b_mafstarttimetai ) = 'MAFStartTimeTAI'
    l1bItemNames(l1b_tpgeodlat       ) = 'tpGeodLat'
    l1bItemNames(l1b_tplon           ) = 'tpLon'
    l1bItemNames(l1b_tpgeodangle     ) = 'tpGeodAngle'
    l1bItemNames(l1b_tpsolarzenith   ) = 'tpSolarZenith'
    l1bItemNames(l1b_tpsolartime     ) = 'tpSolarTime'
    l1bItemNames(l1b_tplosangle      ) = 'tpLosAngle'

    select case ( hGridType )
    case ( l_Fractional )
      if ( .not. got_field(f_Fraction) ) &
        & call announce_error ( root, noFraction )
      l1bItemName = trim(instrumentModuleName)//"."//"tpGeodAngle"
    case ( l_Height )
      if ( .not. got_field(f_height) ) call announce_error ( root, noHeight )
      l1bItemName = trim(instrumentModuleName)//"."//"tpGeodAlt"
    end select

    call ReadL1BData ( l1bInfo%l1boaid, l1bItemName, l1bField, noMAFs, &
      & l1bFlag, firstMAF=chunk%firstMafIndex, lastMAF=chunk%lastMafIndex )
    if ( l1bFlag==-1) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_L1BRead//l1bItemName )
    
    ! allocate default MIFs
    allocate ( defaultMIFs(noMAFs), STAT=status )
    if ( status/=0) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_allocate//"defaultMIFs" )
    
    ! Work out which MIF should have the profile for each MAF.
    select case ( hGridType )
    case ( l_Fractional )
      ! A fractional hGrid, we need to read the tangent point phi
      tpGeodAngle => l1bField%dpField
      
      ! Loop over the MAFs
      do maf = 1, noMAFs
        ! ??? Think about missing data here! ***
        ! Probably need to do a pack on tpGeodAngle and then unpack on
        ! defaultMIFs
        
        minAngle=minval(tpGeodAngle(1,:,maf))
        maxAngle=maxval(tpGeodAngle(1,:,maf))
        
        call Hunt ( tpGeodAngle(1,:,maf), &
          & minAngle+fraction*(maxAngle-minAngle) , defaultMIFs(maf) )
      end do
    case ( l_Height )
      tpGeodAlt => l1bField%dpField
      
      ! Loop over the MAFs
      do maf = 1, noMAFs
        ! ??? Think about missing data here! ***
        ! Probably need to do a pack on tpGeodAngle and then unpack on
        ! defaultMIFs
        
        call Hunt ( tpGeodAlt(1,:,maf), height, defaultMIFs(maf) )
      end do
    end select
    
    ! Done with this piece of l1b data for the moment
    call DeallocateL1BData ( l1bField, l1bFlag )
    
    ! Now we have a default MIFs array; this is a list of the `standard'
    ! MIFs we would choose in the interpolationFactor=1 case.
    ! Work out how many profiles this is going to be.
    
    ! Create an empty hGrid
    hGrid%noProfs = NINT(noMAFs*interpolationFactor)
    call CreateEmptyHGrid(hGrid)
    
    hGrid%noProfsLowerOverlap = 0
    hGrid%noProfsUpperOverlap = 0
    
    ! Setup some arrays
    allocate ( defaultField(noMAFs), interpolatedField(hGrid%noProfs), &
      & STAT=status )
    if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_allocate//"defaultField and/or interpolatedField" )
    
    ! Now we go through all the important geolocation quantities, read them
    ! in, interpolate them if required and store the result in the hGrid
    
    do l1bItem = 1, NoL1BItemsToRead
      ! Get the name of the item to read
      l1bItemName = l1bItemNames(l1bItem)
      if ( l1bItem >= firstModularItem ) l1bItemName = &
        & trim(instrumentModuleName)//"."//l1bItemName
      
      ! Read it from the l1boa file
      call ReadL1BData ( l1bInfo%l1boaid, l1bItemName, l1bField,noMAFs, &
        & l1bFlag, firstMAF=chunk%firstMafIndex, lastMAF=chunk%lastMafIndex )
      if ( l1bFlag==-1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_L1BRead//l1bItemName )
      
      if ( l1bItem==1 ) then       ! do something special for time
        do maf = 1, noMAFs
          defaultField(maf) = l1bField%dpField(1,1,maf) + &
            & (defaultMIFs(maf)-1)*sixth
        end do
      else                         ! Otherwise this is fairly easy.
        do maf = 1, noMAFs
          defaultField(maf) = l1bField%dpField(1,defaultMIFs(maf),maf)
        end do
      end IF
      
      if ( interpolationFactor==1.0 ) then
        interpolatedField = defaultField
      else
        ! ??? Some interpolation is wanted.  I'm going to hold off writing
        ! this because we certaintly don't need it for 0.1 and probably
        ! won't till 1.0.  For the sake of getting things down I'll state
        ! here what I think would be implemented.  One would simply
        ! interpolate from the defaultField to the interpolatedField, using
        ! linear or spline I imagine.  However, there are issues with roll
        ! overs for quantities such as longitude and solarTime.  This is
        ! why I have chosen to defer this piece of code. NJL - 16 December
        ! 1999
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Sorry -- interpolation of hGrids is not yet supported" )
      end if
      
      select case ( l1bItem )
      case ( l1b_MAFStartTimeTAI )
        hGrid%time = interpolatedField
      case ( l1b_tpGeodLat )
        hGrid%geodLat = interpolatedField
      case ( l1b_tpLon )
        hGrid%lon = interpolatedField
      case ( l1b_tpGeodAngle )
        hGrid%phi = interpolatedField
      case ( l1b_tpSolarZenith )
        hGrid%solarZenith = interpolatedField
      case ( l1b_tpSolarTime )
        hGrid%solarTime = interpolatedField
      case ( l1b_tpLosAngle )
        hGrid%losAngle = interpolatedField
      end select
    end do
    
    deallocate ( defaultMIFs, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_DeAllocate // "defaultMIFs" )
    deallocate ( defaultField, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_DeAllocate // "defaultField" )
    deallocate ( interpolatedField, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_DeAllocate // "interpolatedField" )
    
    ! ??? This calculation may need attention! ***
    hGrid%noProfsLowerOverlap = &
      & NINT(chunk%noMAFsLowerOverlap*interpolationFactor)
    hGrid%noProfsUpperOverlap = &
      & NINT(chunk%noMAFsUpperOverlap*interpolationFactor)
    
  end subroutine CreateCommonHGrids

  ! ----------------------------------- CreateEmptyHGrid ---------------
  subroutine CreateEmptyHGrid ( hGrid )
    ! Just does allocates etc.
    type (HGrid_T), intent(inout) :: HGRID

    ! Executable code
    call Allocate_Test ( hGrid%phi, hGrid%noProfs, 'hGrid%phi', ModuleName)
    call Allocate_Test ( hGrid%geodLat, hGrid%noProfs, 'hGrid%geodLat', ModuleName)
    call Allocate_Test ( hGrid%lon, hGrid%noProfs, 'hGrid%lon', ModuleName)
    call Allocate_Test ( hGrid%time, hGrid%noProfs, 'hGrid%time', ModuleName)
    call Allocate_Test ( hGrid%solarTime, hGrid%noProfs, 'hGrid%solarTime', ModuleName)
    call Allocate_Test ( hGrid%solarZenith, hGrid%noProfs, 'hGrid%solarZenith', ModuleName)
    call Allocate_Test ( hGrid%losAngle, hGrid%noProfs, 'hGrid%losAngle', ModuleName)

  end subroutine CreateEmptyHGrid

  ! ---------------------------------------  DestroyHGridContents  -----
  subroutine DestroyHGridContents ( hGrid )

  ! This routine destroys the information associated with an hGrid

    ! Dummy arguments
    type (HGrid_T), intent(inout) :: hGrid

    ! Local Variables
    integer :: STATUS    ! From Deallocate

    ! Executable code
    
    call deallocate_test ( hGrid%phi, 'hGrid%phi', ModuleName )
    call deallocate_test ( hGrid%geodLat, 'hGrid%geodLat', ModuleName )
    call deallocate_test ( hGrid%lon, 'hGrid%lon', ModuleName )
    call deallocate_test ( hGrid%time, 'hGrid%time', ModuleName )
    call deallocate_test ( hGrid%solarTime, 'hGrid%solarTime', ModuleName )
    call deallocate_test ( hGrid%solarZenith, 'hGrid%solarZenith', ModuleName )
    call deallocate_test ( hGrid%losAngle, 'hGrid%losAngle', ModuleName )
    hGrid%noProfs = 0
  end subroutine DestroyHGridContents

  ! ---------------------------------------  DestroyHGridDatabase  -----
  subroutine DestroyHGridDatabase ( database )

  ! This subroutine destroys a quantity template database

    ! Dummy argument
    type (HGrid_T), dimension(:), pointer :: database

    ! Local variables
    integer :: hGridIndex, status
    if ( associated(database) ) then
      do hGridIndex=1,SIZE(database)
        call DestroyHGridContents ( database(hGridIndex) )
      end do
      deallocate ( database, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate // "database" )
    end if
  end subroutine DestroyHGridDatabase

! =====     Private Procedures     =====================================

  ! ---------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( WHERE, CODE )
    integer, intent(in) :: WHERE   ! Tree node where error was noticed
    integer, intent(in) :: CODE    ! Code for error message

    error = max(error,1)
    call output ( '***** At ' )
    call print_source ( source_ref(where) )
    call output ( ': ' )
    select case ( code )
    case ( lengthUnitMessage )
      call output ( "Value for the " )
      call dump_tree_node ( where, 0 )
      call output ( " field is required to be a length, e.g. km", &
        advance='yes' )
    case ( noFraction )
      call output ( "TYPE = FRACTIONAL but no fraction is specified", &
        & advance='yes' )
    case ( noHeight )
      call output ( "TYPE = HEIGHT but no height is specified", advance='yes' )
    case ( unitlessMessage )
      call output ( "Value for the " )
      call dump_tree_node ( where, 0 )
      call output ( " field is required to dimensionless", advance='yes' )
    end select
    end subroutine ANNOUNCE_ERROR

!=============================================================================
end module HGrid
!=============================================================================

!
! $Log$
! Revision 2.12  2001/04/24 22:21:05  livesey
! Gave up on latitude stuff
!
! Revision 2.11  2001/04/23 23:25:26  livesey
! Changed l2gpDatabase to pointer
!
! Revision 2.10  2001/04/21 01:24:55  livesey
! New version, tidied up, can create hGrid from L2GP now
!
! Revision 2.9  2001/04/20 23:11:48  livesey
! Added explicit functionality
!
! Revision 2.8  2001/03/02 01:27:06  livesey
! For new MLSSignals
!
! Revision 2.7  2001/02/22 23:44:29  livesey
! Typo
!
! Revision 2.6  2001/02/22 23:43:26  livesey
! Nullified pointer elements of HGrid_T
!
! Revision 2.5  2001/02/21 01:09:24  livesey
! Tidied stuff up a bit
!
! Revision 2.4  2001/02/09 19:30:16  vsnyder
! Move checking for required and duplicate fields to init_tables_module
!
! Revision 2.3  2001/02/09 00:38:22  livesey
! Various updates
!
! Revision 2.2  2001/02/08 01:50:11  vsnyder
! Move duplicate field checking to tree_checker, set by init_tables
!
! Revision 2.1  2000/12/04 23:34:38  vsnyder
! Move more of addItemToDatabase into the include.
!
! Revision 2.0  2000/09/11 19:18:01  ahanzel
! Changing revision to 2.0.
!
! Revision 1.1  2000/09/07 17:36:29  vsnyder
! Initial version 2.0
!
! Revision 1.9  2000/05/17 23:33:51  lungu
! Added dots between MLSInstrumentModuleName and l1bItemName so that is consistent with L1BOA file.
! Added check "if ( ASSOCIATED(database))deallocate(database)" so it doesn't chrash trying to dealocate
! an "empty" database.
!
! Revision 1.8  2000/05/17 18:15:23  livesey
! Finished off interaction with l2cf.
!

