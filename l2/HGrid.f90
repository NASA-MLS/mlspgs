! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module HGrid                    ! Horizontal grid information
!=============================================================================

  use EXPR_M, only: EXPR
  use INIT_TABLES_MODULE, only: F_FRACTION, F_HEIGHT, F_INTERPOLATIONFACTOR, &
                                F_MODULE, F_TYPE, FIELD_FIRST, FIELD_LAST, &
                                L_FRACTIONAL, L_HEIGHT, LIT_INDICES, &
                                PHYQ_DIMENSIONLESS, PHYQ_LENGTH
  use LEXER_CORE, only: PRINT_SOURCE
  use L1BData, only: DeallocateL1BData, L1BData_T, ReadL1BData
  use MLSCommon, only: L1BInfo_T, MLSChunk_T, NameLen, R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_allocate, &
    & MLSMSG_DeAllocate, MLSMSG_Error, MLSMSG_L1BRead
  use MLSNumerics, only: HUNT
  use OUTPUT_M, only: OUTPUT
  use STRING_TABLE, only: GET_STRING
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TOGGLES, only: GEN, TOGGLE
  use TREE, only: DECORATION, DUMP_TREE_NODE, NSONS, NULL_TREE, SOURCE_REF, &
                  SUB_ROSA, SUBTREE

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
    real(r8), dimension(:), pointer :: phi, geodLat, lon, time, &
      & solarTime, solarZenith, losAngle
  end type HGrid_T

! -----     Private declarations     ---------------------------------

  integer, private :: ERROR

! Error codes for "announce_error"
  integer, private, parameter :: LengthUnitMessage = 1
  integer, private, parameter :: NoFraction = LengthUnitMessage + 1
  integer, private, parameter :: NoHeight = NoFraction + 1
  integer, private, parameter :: NoInstrumentModule = NoHeight + 1
  integer, private, parameter :: NoType = NoInstrumentModule + 1
  integer, private, parameter :: UnitlessMessage = NoType + 1

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
    & ( name, root, l1bInfo, chunk ) result ( hGrid )

  ! This routine creates an hGrid based on the user requests.

    ! Dummy arguments
    integer, intent(in) :: NAME               ! String index of name
    integer, intent(in) :: ROOT               ! Root of hGrid subtree
    type (L1BInfo_T), intent(in) :: l1bInfo   ! File handles for l1b data
    type (MLSChunk_T), intent(in) :: chunk    ! This chunk

    ! Local parameters
    real(r8), parameter :: SIXTH = 1.0_r8 / 6.0_r8

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! IF MODIFYING THIS SECTION PLEASE TAKE CARE, SEE BELOW!
    integer, parameter :: NoL1BItemsToRead=7
    character (len=15), dimension(NoL1BItemsToRead), &
         &   parameter :: L1bItemsToRead = &
         & (/"MAFStartTimeTAI","tpGeodLat      ","tpLon          ",&
         &   "tpGeodAngle    ","tpSolarZenith  ","tpSolarTime    ",&
         &   "tpLosAngle     "/)
    integer, parameter :: TransitionToModularItems = 2
    ! Entries in the above array below TransitionToModularItems are prefixed
    ! with either GHz or THz.  The layout of the above array is critically
    ! bound to the "select case ( l1bItem )" code below.  So TAKE CARE! when
    ! modifing it.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! Local variables
    integer :: EXPR_UNITS(2)            ! Output from Expr subroutine
    double precision :: EXPR_VALUE(2)   ! Output from Expr subroutine
    integer :: hGridType
    real(r8) :: interpolationFactor
    integer :: instrumentModule
    real(r8) :: fraction, height

    integer :: keyNo            ! Entry in the mlscf line

    type (L1BData_T) :: l1bField ! L1B data
    real(r8), dimension(:,:,:), pointer :: tpGeodAngle, tpGeodAlt

    real(r8) :: minAngle, maxAngle

    integer, dimension(:), allocatable :: defaultMIFs
    integer :: FIELD              ! Subtree index of "field" node
    integer :: FIELD_INDEX        ! F_..., see Init_Tables_Module
    logical :: GOT_FIELD(field_first:field_last)
    character (len=NameLen) :: InstrumentModuleName
    integer :: L1BFLAG
    integer :: L1BITEM            ! Loop counter
    character (len=NameLen) :: L1BItemName
    integer :: MAF                ! Loop counters
    integer :: noMAFs             ! Number of MAFs of L1B data read
    integer :: SON                ! Son of Root
    integer :: STATUS             ! From Allocate, ReadL1B... etc.

    ! MIFs it would choose in the non over/undersampled case
    real(r8), dimension(:), allocatable :: defaultField, interpolatedField

    ! Executable code

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
        instrumentModule = decoration(subtree(2,son))
        call get_string ( lit_indices(instrumentModule), instrumentModuleName )
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
      case default ! Can't get here if tree_checker works correctly
      end select
    end do

    ! Now check the sanity of what we have

    if ( .not. got_field(f_type) ) call announce_error ( root, noType )
    if ( .not. got_field(f_module) ) &
      & call announce_error ( root, noInstrumentModule )

    ! This is where we will start reading the l1bdata to get the name to read

    select case ( hGridType )
    case ( l_Fractional )
      if ( .not. got_field(f_Fraction) ) &
        & call announce_error ( root, noFraction )
      l1bItemName = trim(instrumentModuleName)//"."//"tpGeodAngle"
    case ( l_Height )
      if ( .not. got_field(f_height) ) call announce_error ( root, noHeight )
      l1bItemName = trim(instrumentModuleName)//"."//"tpGeodAlt"
    end select

    ! Read the data

    call ReadL1BData ( l1bInfo%l1boaid, l1bItemName, l1bField, noMAFs, &
       & l1bFlag, firstMAF=chunk%firstMafIndex, lastMAF=chunk%lastMafIndex )
    if ( l1bFlag==-1) call MLSMessage ( MLSMSG_Error, ModuleName, &
       & MLSMSG_L1BRead//l1bItemName )

    ! allocate default MIFs

    allocate ( defaultMIFs(noMAFs), STAT=status )
    if ( status/=0) call MLSMessage ( MLSMSG_Error, ModuleName, &
       & MLSMSG_allocate//"defaultMIFs" )
   
    ! Work out which MIF should have the profile for each MAF.

    if ( hGridType==l_Fractional ) then
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
    else
      tpGeodAlt => l1bField%dpField

      ! Loop over the MAFs
      do maf = 1, noMAFs
        ! ??? Think about missing data here! ***
        ! Probably need to do a pack on tpGeodAngle and then unpack on
        ! defaultMIFs

        call Hunt ( tpGeodAlt(1,:,maf), height, defaultMIFs(maf) )
      end do
    end if

    ! Done with this piece of l1b data for the moment
    call DeallocateL1BData ( l1bField, l1bFlag )

    ! Now we have a default MIFs array; this is a list of the `standard'
    ! MIFs we would choose in the interpolationFactor=1 case.
    ! Work out how many profiles this is going to be.

    ! Create an empty hGrid
    hGrid%name = name
    hGrid%noProfs = NINT(noMAFs*interpolationFactor)
    allocate ( hGrid%phi(hGrid%noProfs), hGrid%geodLat(hGrid%noProfs), &
      & hGrid%lon(hGrid%noProfs), hGrid%time(hGrid%noProfs), &
      & hGrid%solarTime(hGrid%noProfs), hGrid%solarZenith(hGrid%noProfs), &
      & hGrid%losAngle(hGrid%noProfs), STAT=status )
    if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "hGrid information" )

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
      l1bItemName = l1bItemsToRead(l1bItem)
      if ( l1bItem >= TransitionToModularItems ) l1bItemName = &
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

      ! Now we have to save this field in the hGrid data.  This is rather a
      ! kludgy way of doing it but this worked out the least boring way to
      ! write the code.  See the definition of L1BItemsToRead above for
      ! reference.

      select case ( l1bItem )
      case ( 1 )
        hGrid%time = interpolatedField
      case ( 2 )
        hGrid%geodLat = interpolatedField
      case ( 3 )
        hGrid%lon = interpolatedField
      case ( 4 )
        hGrid%phi = interpolatedField
      case ( 5 )
        hGrid%solarZenith = interpolatedField
      case ( 6 )
        hGrid%solarTime = interpolatedField
      case ( 7 )
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

    if ( toggle(gen) ) call trace_end ( "CreateHGridFromMLSCFInfo" )
    
  end function CreateHGridFromMLSCFInfo

  ! ---------------------------------------  DestroyHGridContents  -----
  subroutine DestroyHGridContents ( hGrid )

  ! This routine destroys the information associated with an hGrid

    ! Dummy arguments
    type (HGrid_T), intent(out) :: hGrid

    ! Local Variables
    integer :: STATUS    ! From Deallocate

    ! Executable code

    deallocate ( hGrid%phi, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_DeAllocate // "hGrid%phi" )
    deallocate ( hGrid%geodLat, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_DeAllocate // "hGrid%geodLat" )
    deallocate ( hGrid%lon, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_DeAllocate // "hGrid%lon" )
    deallocate ( hGrid%time, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_DeAllocate // "hGrid%time" )
    deallocate ( hGrid%solarTime, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_DeAllocate // "hGrid%solarTime" )
    deallocate ( hGrid%solarZenith, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_DeAllocate // "hGrid%solarZenith" )
    deallocate ( hGrid%losAngle, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_DeAllocate // "hGrid%losAngle" )

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
    case ( noInstrumentModule )
      call output ( "No instrument module specified for the hGrid", &
        & advance='yes' )
    case ( noType )
      call output ( "No type specified for the hGrid", advance='yes' )
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

