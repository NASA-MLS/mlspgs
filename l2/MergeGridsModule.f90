! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MergeGridsModule

  ! This module contains code for merging operational gridded data with apriori
  ! information.

  use Expr_m, only: EXPR
  use Init_tables_module, only: F_CLIMATOLOGY, F_HEIGHT, F_OPERATIONAL, &
    & F_SCALE, S_MERGE
  use GriddedData, only: GRIDDEDDATA_T, SETUPNEWGRIDDEDDATA, &
    & ADDGRIDDEDDATATODATABASE, V_IS_PRESSURE
  use L3ASCII, only: L3ASCII_INTERP_FIELD
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR
  use MoreTree, only: GET_SPEC_ID
  use Trace_M, only: TRACE_BEGIN, TRACE_END
  use Tree, only: NSONS, SUBTREE, DECORATE, DECORATION, NODE_ID, SUB_ROSA
  use Tree_Types, only: N_NAMED
  use Toggles, only: GEN, TOGGLE
  use Units, only: PHYQ_Pressure, PHYQ_Length

  implicit none
  private

  public :: MergeGrids

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  !---------------------------------------------------------------------------

contains ! =================================== Public procedures

  ! ----------------------------------------- MergeGrid

  subroutine MergeGrids ( root, griddedData )
    integer, intent(in) :: ROOT         ! Tree root
    type (GriddedData_T), dimension(:), pointer :: griddedData ! Database
    
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
        ! Shouldn't get here if parser worked?
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Expecting only named specifiers in MergeGrids section' )
      end if

      if ( get_spec_id(key) == s_merge ) then
        call decorate ( key, AddGriddedDataToDatabase ( griddedData, &
          & MergeOneGrid ( key, griddedData ) ) )
      else
        ! Shouldn't get here is parser worked?
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Only merge commands allowed in MergeGrids section' )
      end if
    end do
    if ( toggle(gen) ) call trace_end ( "MergeGrids", root )
  end subroutine MergeGrids

  ! ----------------------------------------- MergeOneGrid
  type (GriddedData_T) function MergeOneGrid ( root, griddedData ) &
    & result ( newGrid )
    integer, intent(in) :: ROOT         ! Tree node
    type (GriddedData_T), dimension(:), pointer :: GRIDDEDDATA ! Database

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
    integer :: SURF                     ! Loop counter
    integer :: SZA                      ! Loop counter
    integer :: UNITS(2)                 ! Units for expr
    integer :: VALUE                    ! Value of tree node

    real (r8) :: CLIVAL                 ! Value from climatological grid
    real (r8) :: CLIWEIGHT              ! Climtaological 'weight'
    real (r8) :: HEIGHT                 ! Transition height
    real (r8) :: OPVAL                  ! Value from operational grid
    real (r8) :: OPWEIGHT               ! Operational 'weight'
    real (r8) :: SCALE                  ! Transition scale
    real (r8) :: TOTALWEIGHT            ! Total weight
    real (r8) :: VALUES(2)              ! Value of expr
    real (r8) :: ZTRANS                 ! Transition 'height'
    real (r8) :: Z                      ! One 'height'
    real (r8) :: Z1, Z2                 ! Range of transition region

    type (GriddedData_T), pointer :: OPERATIONAL
    type (GriddedData_T), pointer :: CLIMATOLOGY

    ! Executable code

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
        operational => griddedData ( decoration ( decoration ( value ) ) )
      case ( f_climatology )
        climatology => griddedData ( decoration ( decoration ( value ) ) )
      case ( f_height )
        call expr ( value, units, values )
        if ( units(1) /= phyq_pressure ) call MLSMessage ( &
          & MLSMSG_Error, ModuleName, &
          & 'Only pressure is allowed for the height field' )
        height = values(1)
      case ( f_scale )
        call expr ( value, units, values )
        if ( units(1) /= phyq_length ) call MLSMessage ( &
          & MLSMSG_Error, ModuleName, &
          & 'Only altitude is allowed for the scale field' )
        scale = values(1)
      end select
    end do
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
                call l3ascii_interp_field ( &
                  & climatology, &
                  & cliVal, &
                  & newGrid%heights(surf), &
                  & newGrid%lats(lat), &
                  & newGrid%lons(lon), &
                  & newGrid%lsts(lst), &
                  & newGrid%szas(sza), &
                  & 0.5 * ( newGrid%dateStarts(day)+newGrid%dateEnds(day) ) )
                call l3ascii_interp_field ( &
                  & operational, &
                  & opVal, &
                  & newGrid%heights(surf), &
                  & newGrid%lats(lat), &
                  & newGrid%lons(lon), &
                  & newGrid%lsts(lst), &
                  & newGrid%szas(sza), &
                  & 0.5 * ( newGrid%dateStarts(day)+newGrid%dateEnds(day) ) )
                ! Now work out the weighting of the two
                z = scaleHeight * ( 3.0 - log10 ( newGrid%heights(surf) ) )
                if ( scale /= 0.0 ) then
                  cliWeight = ( z - z1 ) / ( z2-z1 )
                  opWeight = 1.0 - cliWeight
                end if

                ! Would probably think about missing data here
                cliWeight = min ( max ( cliWeight, 0.0_r8 ), 1.0_r8 )
                opWeight = min ( max ( cliWeight, 0.0_r8 ), 1.0_r8 )
                totalWeight = cliWeight + opWeight

                ! OK, store this value
                if ( totalWeight > 0.0_r8 ) then
                  newGrid%field ( surf, lat, lon, lst, sza, day ) = &
                    & ( cliWeight*cliVal + opWeight*opVal ) / totalWeight
                else
!                   newGrid%field ( surf, lat, lon, lst, sza, day ) = &
!                     & newGrid%missing
                end if
              end do
            end do
          end do
        end do
      end do
    end do
    if ( toggle(gen) ) call trace_end ( "MergeOneGrid", root )
  end function MergeOneGrid

end module MergeGridsModule

! $Log$
! Revision 2.2  2002/01/26 00:10:54  livesey
! It compiles at least
!
! Revision 2.1  2002/01/24 00:58:03  livesey
! First version, not much more than a stub
!
