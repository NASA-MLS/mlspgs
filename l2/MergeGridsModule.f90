! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MergeGridsModule

  ! This module contains code for merging operational gridded data with apriori
  ! information.

  use Init_tables_module, only: S_MERGE
  use GriddedData, only: GRIDDEDDATA_T, SETUPNEWGRIDDEDDATA, &
    & ADDGRIDDEDDATATODATABASE
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR
  use MoreTree, only: GET_SPEC_ID
  use Tree, only: NSONS, SUBTREE, DECORATE, DECORATION, NODE_ID, SUB_ROSA
  use Tree_Types, only: N_NAMED

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
  end subroutine MergeGrids

  ! ----------------------------------------- MergeOneGrid
  type (GriddedData_T) function MergeOneGrid ( root, griddedData ) &
    & result ( newGrid )
    integer, intent(in) :: ROOT         ! Tree node
    type (GriddedData_T), dimension(:), intent(in) :: GRIDDEDDATA ! Database

    ! Executable code
  end function MergeOneGrid

end module MergeGridsModule

! $Log$
! Revision 2.1  2002/01/24 00:58:03  livesey
! First version, not much more than a stub
!
