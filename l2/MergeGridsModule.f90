! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MergeGridsModule

  ! This module contains code for merging operational gridded data with apriori
  ! information.

  implicit none
  private

  public :: MergeGrids

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains ! =================================== Public procedures

  ! ----------------------------------------- MergeGrid

  subroutine MergeGrids ( root, griddedDataBase )

    use GriddedData, only: GRIDDEDDATA_T, RGR, SETUPNEWGRIDDEDDATA, &
      & ADDGRIDDEDDATATODATABASE, NULLIFYGRIDDEDDATA, &
      & WRAPGRIDDEDDATA, CONCATENATEGRIDDEDDATA, COPYGRID
    use Init_tables_module, only: S_MERGE, S_CONCATENATE, S_DELETE
    use MLSCommon, only: R8
    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR
    use MoreTree, only: GET_SPEC_ID
    use Trace_M, only: TRACE_BEGIN, TRACE_END
    use Tree, only: NSONS, SUBTREE, DECORATE, DECORATION, NODE_ID, SUB_ROSA
    use Tree_Types, only: N_NAMED
    use Toggles, only: GEN, TOGGLE

    integer, intent(in) :: ROOT         ! Tree root
    type (griddedData_T), dimension(:), pointer :: griddedDataBase ! Database

    ! Local variables
    integer :: I                        ! Loop counter
    integer :: SON                      ! Tree node
    integer :: KEY                      ! Another node
    integer :: NAME                     ! Index into string table

    ! excutable code
    if ( toggle(gen) ) call trace_begin ( "MergeGrids", root )
    do i = 2, nsons(root) - 1           ! Skip the begin and end stuff
      son = subtree ( i, root )
      if ( node_id(son) == n_named ) then ! Is spec labed?
        key = subtree ( 2, son )
        name = sub_rosa ( subtree(1,son) )
      else
        key = son
      end if

      select case ( get_spec_id(key) )
      case ( s_merge )
        call decorate ( key, AddgriddedDataToDatabase ( griddedDataBase, &
          & MergeOneGrid ( key, griddedDataBase ) ) )
      case ( s_concatenate )
        call decorate ( key, AddgriddedDataToDatabase ( griddedDataBase, &
          & Concatenate ( key, griddedDataBase ) ) )
      case ( s_delete )
        call DeleteGriddedData ( key, griddedDatabase )
      case default
        ! Shouldn't get here is parser worked?
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Only merge and concatenate commands allowed in MergeGrids section' )
      end select
    end do
    if ( toggle(gen) ) call trace_end ( "MergeGrids" )

  end subroutine MergeGrids

  ! ---------------------------------------- Concatenate --
  type (griddedData_T) function Concatenate ( root, griddedDataBase ) &
    & result ( newGrid )
    use Tree, only: NSONS, SUBTREE, DECORATION
    use GriddedData, only: GRIDDEDDATA_T, NULLIFYGRIDDEDDATA, &
      & CONCATENATEGRIDDEDDATA, COPYGRID
    use Trace_M, only: TRACE_BEGIN, TRACE_END
    use Init_tables_module, only: F_A, F_B
    use Toggles, only: GEN, TOGGLE
    
    integer, intent(in) :: ROOT         ! Tree node
    type (griddedData_T), dimension(:), pointer :: griddedDataBase ! Database
    ! This routine parses the l2cf instructions that request
    ! a grid concatenation, and then performs the concatenation

    ! Local variables
    integer :: SON                    ! Tree node
    integer :: FIELD                  ! Another tree node
    integer :: FIELD_INDEX            ! Type of tree node
    integer :: VALUE                  ! Tree node
    integer :: I                      ! Loop counter

    type (griddedData_T), pointer :: A
    type (griddedData_T), pointer :: B

    ! Executable code
    call nullifyGriddedData ( newGrid ) ! for Sun's still useless compiler
    if ( toggle(gen) ) call trace_begin ( "Concatenate", root )

    ! Get the information from the l2cf    
    ! Note that init_tables_module has insisted that we have all
    ! arguments so we don't need a 'got' type arrangement
    do i = 2, nsons(root)
      son = subtree(i,root)
      field = subtree(1,son)
      value = subtree(2,son)
      field_index = decoration(field)
      select case ( field_index )
      case ( f_a ) 
        a => griddedDataBase ( decoration ( decoration ( value ) ) )
      case ( f_b )
        b => griddedDataBase ( decoration ( decoration ( value ) ) )
      end select
    end do

    ! Do the concatenation unless one or other is empty
    if ( .not. a%empty .and. .not. b%empty ) then
      call ConcatenateGriddedData ( A, B, newGrid )
    else if ( a%empty ) then
      ! Copy B into the result, of course, that may be empty too
      ! in which case the result is empty, no problem!
      call CopyGrid ( newGrid, b )
    else
      ! Otherwise a must be full, b empty
      call CopyGrid ( newGrid, a )
    end if

    if ( toggle(gen) ) call trace_end ( "Concatenate" )
  end function Concatenate

  ! ------------------------------------ DeleteGriddedData ---
  subroutine DeleteGriddedData ( root, griddedDataBase )
    use Tree, only: NSONS, SUBTREE, DECORATION
    use GriddedData, only: DESTROYGRIDDEDDATA, GRIDDEDDATA_T
    use Init_Tables_Module, only: F_GRID
    ! This routine deletes the grid indicated by the l2cf
    integer, intent(in) :: ROOT         ! Tree node
    type (griddedData_T), dimension(:), pointer :: GRIDDEDDATABASE ! Database
    ! Local variables
    type (griddedData_T), pointer :: GRID
    integer :: FIELD                    ! Tree node
    integer :: FIELD_INDEX              ! Tree node type
    integer :: I                        ! Counter
    integer :: SON                      ! Tree node
    integer :: VALUE                    ! Tree node

    ! Get the information from the l2cf    
    ! Note that init_tables_module has insisted that we have all
    ! arguments so we don't need a 'got' type arrangement
    ! In this case there is only one argument anyway
    do i = 2, nsons(root)
      son = subtree(i,root)
      field = subtree(1,son)
      value = subtree(2,son)
      field_index = decoration(field)
      select case ( field_index )
      case ( f_grid ) 
        grid => griddedDataBase ( decoration ( decoration ( value ) ) )
      end select
    end do
    call DestroyGriddedData ( grid )
  end subroutine DeleteGriddedData

  ! ----------------------------------------- MergeOneGrid
  type (griddedData_T) function MergeOneGrid ( root, griddedDataBase ) &
    & result ( newGrid )
    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use Tree, only: NSONS, SUBTREE, DECORATION
    use Units, only: PHYQ_Length, PHYQ_Pressure
    use GriddedData, only: GRIDDEDDATA_T, NULLIFYGRIDDEDDATA, COPYGRID, &
      & WRAPGRIDDEDDATA, SETUPNEWGRIDDEDDATA, RGR, V_IS_PRESSURE, SLICEGRIDDEDDATA
    use MLSNumerics, only: ESSENTIALLYEQUAL
    use L3ASCII, only: L3ASCII_INTERP_FIELD
    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR, MLSMSG_ALLOCATE
    use Trace_M, only: TRACE_BEGIN, TRACE_END
    use Init_tables_module, only: F_CLIMATOLOGY, F_HEIGHT, F_OPERATIONAL, &
      & F_SCALE
    use MLSCommon, only: R8
    use Toggles, only: GEN, TOGGLE
    use Expr_m, only: EXPR

    integer, intent(in) :: ROOT         ! Tree node
    type (griddedData_T), dimension(:), pointer :: griddedDataBase ! Database

    ! This routine creates a new grid being a merge of two others.
    ! The operational grid forms the bottom of the dataset
    ! The climatology grid the top.  The result has the horizontal
    ! coordinates to operational and the vertical coordiantes of climatology

    ! Note! This routine is far from efficient. Not least because it
    ! uses l3atascii_interp_field which isn't terribly efficient either.
    ! But hey, this isn't a key part of the software when it comes
    ! to a desire for speed (or at least it shouldn't be)

    ! I'll need to think about missing data at some point.

    ! Local parameters
    real (r8), parameter :: SCALEHEIGHT = 16.0e3_r8 ! Approximate scale height / m

    ! Local variables
    integer :: DAY                      ! Loop counter
    integer :: FIELD                    ! Another tree node
    integer :: FIELD_INDEX              ! Type of tree node
    integer :: I                        ! Loop inductor
    integer :: LAT                      ! Loop counter
    integer :: LON                      ! Loop counter
    integer :: LST                      ! Loop counter
    integer :: SON                      ! Tree node
    integer :: STATUS                   ! Flag from allocate
    integer :: SURF                     ! Loop counter
    integer :: SZA                      ! Loop counter
    integer :: EXPRUNITS(2)             ! Units for expr
    integer :: VALUE                    ! Value of tree node

    real (r8) :: CLIWEIGHT              ! Climatological 'weight'
    real (r8) :: HEIGHT                 ! Transition height
    real (r8) :: OPWEIGHT               ! Operational 'weight'
    real (r8) :: SCALE                  ! Transition scale
    real (r8) :: TOTALWEIGHT            ! Total weight
    real (r8) :: EXPRVALUES(2)          ! Value of expr
    real (r8) :: ZTRANS                 ! Transition 'height'
    real (r8) :: Z                      ! One 'height'
    real (r8) :: Z1, Z2                 ! Range of transition region

    real (rgr) :: CLIVAL                ! One interpolated value
    real (rgr) :: OPVAL                 ! One interpolated value
    real (rgr), pointer, dimension(:,:,:,:,:,:) :: CLIMAPPED
    real (rgr), pointer, dimension(:,:,:,:,:,:) :: OPERMAPPED
    real (r8), dimension(:), pointer :: MEANDATES ! Mean dates for new grid

    type (griddedData_T), pointer :: OPERATIONAL
    type (griddedData_T), pointer :: CLIMATOLOGY

    ! Executable code

    call nullifyGriddedData ( newGrid ) ! for Sun's still useless compiler
    if ( toggle(gen) ) call trace_begin ( "MergeOneGrid", root )

    ! Get the information from the l2cf    
    ! Note that init_tables_module has insisted that we have all
    ! arguments so we don't need a 'got' type arrangement
    do i = 2, nsons(root)
      son = subtree(i,root)
      field = subtree(1,son)
      value = subtree(2,son)
      field_index = decoration(field)
      select case ( field_index )
      case ( f_operational ) 
        operational => griddedDataBase ( decoration ( decoration ( value ) ) )
      case ( f_climatology )
        climatology => griddedDataBase ( decoration ( decoration ( value ) ) )
      case ( f_height )
        call expr ( value, exprUnits, exprValues )
        if ( exprUnits(1) /= phyq_pressure ) call MLSMessage ( &
          & MLSMSG_Error, ModuleName, &
          & 'Only pressure is allowed for the height field' )
        height = exprValues(1)
      case ( f_scale )
        call expr ( value, exprUnits, exprValues )
        if ( exprUnits(1) /= phyq_length ) call MLSMessage ( &
          & MLSMSG_Error, ModuleName, &
          & 'Only altitude is allowed for the scale field' )
        scale = exprValues(1)
      end select
    end do

    ! Think about cases where one or other grid is empty
    if ( climatology%empty ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'The climatology grid for the merge is empty' )
    if ( operational%empty ) then
      ! If no operational data, then just use climatology
      call CopyGrid ( newGrid, climatology )
      return
    end if

    ! Do some final sanity checks
    if ( operational%verticalCoordinate /= v_is_pressure ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Operational grid not on pressure surfaces' )
    if ( climatology%verticalCoordinate /= v_is_pressure ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Climatology grid not on pressure surfaces' )
    !     if ( climatology%units /= operational%units ) &
    !       & call MLSMessage ( MLSMSG_Error, ModuleName, &
    !       & 'The climatology and operational data describe different physical quantities' )
    if ( climatology%equivalentLatitude .neqv. operational%equivalentLatitude ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'The climatology and operational data are mixed latitude/equivalent latitude.' )

    ! OK, now we're ready to go.
    ! First we're going to 'wrap' the climatology to be sure that we can
    ! interpolate it in longitude.  The chances are that it has no longitudinal
    ! variation anyway, so this won't actually do anything
    call WrapGriddedData ( climatology )
    ! Create the result.  It has the same vertical coordinates as climatology
    ! But the same horizontal coordinates as operational.
    call SetupNewGriddedData ( newGrid, source=operational, &
      & noHeights=climatology%noHeights )
    ! Setup the rest of the quantity
    newGrid%sourceFileName = 'Result of merge'
    newGrid%quantityName   = 'Result of merge'
    newGrid%description    = 'Result of merge'
    newGrid%units          = climatology%units
    newGrid%verticalCoordinate = v_is_pressure
    newGrid%equivalentLatitude = climatology%equivalentLatitude
    newGrid%heights = climatology%heights
    newGrid%lats = operational%lats
    newGrid%lons = operational%lons
    newGrid%lsts = operational%lsts
    newGrid%szas = operational%szas
    newGrid%dateStarts = operational%dateStarts
    newGrid%dateEnds = operational%dateEnds

    ! Get the 'mean' dates for the result
    nullify ( meanDates )
    call Allocate_test ( meanDates, newGrid%noDates, 'meanDates', ModuleName )
    meanDates = ( newGrid%dateStarts + newGrid%dateEnds ) / 2.0

    ! Now create two fields the same shape as the new field that contain
    ! the operational and climatological data interpolated to our new locations.
    allocate ( operMapped ( &
      & newGrid%noHeights, newGrid%noLats, newGrid%noLons, &
      & newGrid%noLsts, newGrid%noSzas, newGrid%noDates ), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'operMapped' )
    allocate ( cliMapped ( &
      & newGrid%noHeights, newGrid%noLats, newGrid%noLons, &
      & newGrid%noLsts, newGrid%noSzas, newGrid%noDates ), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'operMapped' )

    call SliceGriddedData ( operational, operMapped, &
      & newGrid%heights, newGrid%lats, newGrid%lons, newGrid%lsts, &
      & newGrid%szas, meanDates, missingValue=newGrid%missingValue )
    call SliceGriddedData ( climatology, cliMapped, &
      & newGrid%heights, newGrid%lats, newGrid%lons, newGrid%lsts, &
      & newGrid%szas, meanDates, missingValue=newGrid%missingValue )

    zTrans = scaleHeight * ( 3.0 - log10 ( height ) )
    z1 = zTrans - scale/2.0
    z2 = zTrans + scale/2.0

    ! Now we're going to fill in the rest of the field
    do day = 1, newGrid%noDates
      do sza = 1, newGrid%noSzas
        do lst = 1, newGrid%noLsts
          do lon = 1, newGrid%noLons
            do lat = 1, newGrid%noLats
              do surf = 1, newGrid%noHeights
                ! Get the values
                cliVal = cliMapped ( surf, lat, lon, lst, sza, day )
                opVal = operMapped ( surf, lat, lon, lst, sza, day )
                ! Weight them by height
                z = scaleHeight * ( 3.0 - log10 ( newGrid%heights(surf) ) )
                if ( scale /= 0.0 ) then
                  cliWeight = ( z - z1 ) / ( z2-z1 )
                  opWeight = 1.0 - cliWeight
                end if

                cliWeight = min ( max ( cliWeight, 0.0_r8 ), 1.0_r8 )
                opWeight = min ( max ( opWeight, 0.0_r8 ), 1.0_r8 )

                ! Check for bad data in operational dataset
                if ( EssentiallyEqual ( opVal, newGrid%missingValue ) ) opWeight = 0.0

                ! Check for bad data in the climatology
                if ( EssentiallyEqual ( cliVal, newGrid%missingValue ) ) &
                  & call MLSMessage ( MLSMSG_Error, ModuleName, &
                  & 'There is a bad data point in the climatology field' )

                totalWeight = cliWeight + opWeight
                if ( totalWeight == 0.0 ) then
                  ! Presumably was in the region where operational was supposed
                  ! to dominate, but it's bad, so switch to a priori
                  cliWeight = 1.0
                  totalWeight = 1.0
                end if

                ! OK, store this value
                newGrid%field ( surf, lat, lon, lst, sza, day ) = &
                  & ( cliWeight*cliVal + opWeight*opVal ) / totalWeight
              end do
            end do
          end do
        end do
      end do
    end do

    ! Tidy up
    call Deallocate_test ( meanDates, 'meanDates', ModuleName )
    deallocate ( cliMapped, operMapped, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'operMapped or cliMapped' )

    if ( toggle(gen) ) call trace_end ( "MergeOneGrid" )
  end function MergeOneGrid

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MergeGridsModule

! $Log$
! Revision 2.13  2003/06/06 01:06:59  livesey
! Added DeleteGrids stuff
!
! Revision 2.12  2003/05/09 02:13:05  livesey
! Removed a dump
!
! Revision 2.11  2003/05/09 01:55:14  livesey
! Sped up Merge by using new SliceGriddedData routine.
!
! Revision 2.10  2003/04/04 00:08:26  livesey
! Added Concatenate capability, various reorganizations.
!
! Revision 2.9  2003/02/28 02:33:28  livesey
! Bug fix, careless with the old emacs.
!
! Revision 2.8  2003/02/28 02:25:50  livesey
! First working version.
!
! Revision 2.7  2003/02/19 19:15:13  pwagner
! Consistent with new GriddedData_T and rgr
!
! Revision 2.6  2002/11/22 12:21:14  mjf
! Added nullify routine(s) to get round Sun's WS6 compiler not
! initialising derived type function results.
!
! Revision 2.5  2002/10/08 17:36:21  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.4  2002/08/22 20:26:09  vsnyder
! Move another USE from module scope to procedure scope
!
! Revision 2.3  2002/08/21 02:23:39  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.2  2002/01/26 00:10:54  livesey
! It compiles at least
!
! Revision 2.1  2002/01/24 00:58:03  livesey
! First version, not much more than a stub
!
