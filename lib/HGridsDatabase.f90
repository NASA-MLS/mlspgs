! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module HGridsDatabase                   ! Horizontal grid information
  
  use MLSCommon, only: r8

  implicit none
  private

  public :: HGrid_T, AddHGridToDatabase, CreateEmptyHGrid, TrimHGrid, &
    & DestroyHGridContents, DestroyHGridDatabase, NullifyHGrid

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

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
    real(r8), dimension(:,:), pointer :: phi => NULL()
    real(r8), dimension(:,:), pointer :: geodLat => NULL()
    real(r8), dimension(:,:), pointer :: lon => NULL()
    real(r8), dimension(:,:), pointer :: time => NULL()
    real(r8), dimension(:,:), pointer :: solarTime => NULL()
    real(r8), dimension(:,:), pointer :: solarZenith => NULL()
    real(r8), dimension(:,:), pointer :: losAngle => NULL()
    integer, dimension(:), pointer :: mafIndex => NULL()
    integer, dimension(:), pointer :: mafCounter => NULL()
  end type HGrid_T

contains ! =========== Public procedures ===================================

 ! -----------------------------------------  AddHGridToDatabase  -----
  integer function AddHGridToDatabase ( database, item )
    
    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, &
      & MLSMSG_DeAllocate, MLSMSG_Error

    ! Dummy arguments
    type (HGrid_T), dimension(:), pointer :: database
    type (HGrid_T), intent(in) :: item

    ! Local variables
    type (HGrid_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddHGridToDatabase = newSize
  end function AddHGridToDatabase

  ! -------------------------------------------  CreateEmptyHGrid  -----
  subroutine CreateEmptyHGrid ( hGrid )
    ! Just does allocates etc.

    use Allocate_Deallocate, only: Allocate_Test

    type (HGrid_T), intent(inout) :: HGRID

    ! Executable code
    call Allocate_Test ( hGrid%phi, 1, hGrid%noProfs, &
      & 'hGrid%phi', ModuleName )
    call Allocate_Test ( hGrid%geodLat, 1, hGrid%noProfs, &
      & 'hGrid%geodLat', ModuleName )
    call Allocate_Test ( hGrid%lon, 1, hGrid%noProfs, &
      & 'hGrid%lon', ModuleName )
    call Allocate_Test ( hGrid%time, 1, hGrid%noProfs, &
      & 'hGrid%time', ModuleName )
    call Allocate_Test ( hGrid%solarTime, 1, hGrid%noProfs, &
      & 'hGrid%solarTime', ModuleName )
    call Allocate_Test ( hGrid%solarZenith, 1, hGrid%noProfs, &
      & 'hGrid%solarZenith', ModuleName )
    call Allocate_Test ( hGrid%losAngle, 1, hGrid%noProfs, &
      & 'hGrid%losAngle', ModuleName )

  end subroutine CreateEmptyHGrid

  ! -------------------------------------------------  TrimHGrid  ------
  subroutine TrimHGrid ( hGrid, side, NOTODELETE )

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use MLSCommon, only: R8
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error

    type (HGrid_T), intent(inout) :: HGrid
    integer, intent(in) :: SIDE         ! -1 = lower, 1 = upper
    integer, intent(in) :: NOTODELETE ! How many to delete, default all
    ! Local variables
    integer :: newNoProfs
    real(r8), dimension(:), pointer :: temp
    integer :: first, last

    ! Executable code
    select case ( side )
    case ( -1 )
      newNoProfs = hGrid%noProfs - noToDelete
      first = noToDelete + 1
      last = hGrid%noProfs
      hGrid%noProfsLowerOverlap = max ( hGrid%noProfsLowerOverlap - noToDelete, 0 )
    case ( 1 )
      newNoProfs = hGrid%noProfs - noToDelete
      first = 1
      last = newNoProfs
      hGrid%noProfsUpperOverlap = max ( hGrid%noProfsUpperOverlap - noToDelete, 0 )
    case default
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Invalid side argument to TrimHGrid' )
    end select
    if ( newNoProfs <= 0 ) call MLSMessage ( &
        & MLSMSG_Error, ModuleName, 'Too many profiles to delete' )

    nullify ( temp )
    call allocate_test ( temp, hGrid%noProfs, 'temp', ModuleName )
    hGrid%noProfs = newNoProfs

    ! Now allocate each entry and trim it
    temp = hGrid%phi(1,:)               ! ------------------------- Phi
    call Allocate_Test ( hGrid%phi, 1, newNoProfs, 'hGrid%phi', ModuleName)
    hGrid%phi(1,:) = temp ( first : last )
    temp = hGrid%geodLat(1,:)           ! ------------------------- GeodLat
    call Allocate_Test ( hGrid%geodLat, 1, newNoProfs, 'hGrid%geodLat', ModuleName)
    hGrid%geodLat(1,:) = temp ( first : last )
    temp = hGrid%lon(1,:)               ! ------------------------- Lon
    call Allocate_Test ( hGrid%lon, 1, newNoProfs, 'hGrid%lon', ModuleName)
    hGrid%lon(1,:) = temp ( first : last )
    temp = hGrid%time(1,:)              ! ------------------------- Time
    call Allocate_Test ( hGrid%time, 1, newNoProfs, 'hGrid%time', ModuleName)
    hGrid%time(1,:) = temp ( first : last )
    temp = hGrid%solarTime(1,:)         ! ------------------------- SolarTime
    call Allocate_Test ( hGrid%solarTime, 1, newNoProfs, 'hGrid%solarTime', ModuleName)
    hGrid%solarTime(1,:) = temp ( first : last )
    temp = hGrid%solarZenith(1,:)       ! ------------------------- SolarZenith
    call Allocate_Test ( hGrid%solarZenith, 1, newNoProfs, 'hGrid%solarZenith', ModuleName)
    hGrid%solarZenith(1,:) = temp ( first : last )
    temp = hGrid%losAngle(1,:)          ! ------------------------- LosAngle
    call Allocate_Test ( hGrid%losAngle, 1, newNoProfs, 'hGrid%losAngle', ModuleName)
    hGrid%losAngle(1,:) = temp ( first : last )

    call Deallocate_test ( temp, 'temp', ModuleName )
  end subroutine TrimHGrid
    
  ! ---------------------------------------  DestroyHGridContents  -----
  subroutine DestroyHGridContents ( hGrid )

    use Allocate_Deallocate, only: Deallocate_Test

  ! This routine destroys the information associated with an hGrid

    ! Dummy arguments
    type (HGrid_T), intent(inout) :: hGrid

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

    use MLSMessageModule, only: MLSMessage, MLSMSG_DeAllocate, MLSMSG_Error

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

  ! ----------------------------------------NullifyHGrid -----
  subroutine NullifyHGrid ( H )
    ! Given a hGrid, nullify all the pointers associated with it
    type ( HGrid_T ), intent(out) :: H

    ! Executable code
    nullify ( h%phi )
    nullify ( h%geodLat )
    nullify ( h%lon )
    nullify ( h%time )
    nullify ( h%solarTime )
    nullify ( h%solarZenith )
    nullify ( h%losAngle )
  end subroutine NullifyHGrid

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module HGridsDatabase

! $Log$
! Revision 2.1  2003/06/20 19:34:45  pwagner
! Quanities now share grids stored separately in databses
!
