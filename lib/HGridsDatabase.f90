! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module HGridsDatabase                   ! Horizontal grid information
  
  use MLSCommon, only: r8

  implicit none
  private

  public :: HGrid_T
  public :: AddHGridToDatabase, CreateEmptyHGrid, DestroyHGridContents, &
    & DestroyHGridDatabase, Dump, FindClosestMatch, NullifyHGrid, &
    & TrimHGrid

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
  end type HGrid_T

  interface DUMP
    module procedure DUMP_a_HGRID
    module procedure DUMP_HGRIDS
  end interface

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
    if ( newNoProfs < 0 ) call MLSMessage ( &
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

  ! ------------------------------------------------  DUMP_A_HGRID  -----
  subroutine DUMP_a_HGRID ( aHGRID )
    use OUTPUT_M, only: OUTPUT
    use STRING_TABLE, only: DISPLAY_STRING
    type(hGrid_T), intent(in) :: aHGRID
    integer :: J
      call output ( 'Name = ' )
      call display_string ( aHgrid%name )
      call output ( aHgrid%noProfs, before=' noProfs = ' )
      call output ( aHgrid%noProfsLowerOverlap, before=' lowerOverlap = ' )
      call output ( aHgrid%noProfsUpperOverlap, before=' upperOverlap = ', advance='yes' )
      call output ( ' prof       phi       geodLat           lon' )
      call output ( '          time     solarTime   solarZenith' )
      call output ( '      losAngle', advance='yes' )
      do j = 1, aHgrid%noProfs
        call output ( j, places=5 )
        call output ( aHgrid%phi(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%geodLat(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%lon(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%time(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%solarTime(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%solarZenith(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%losAngle(1,j), '(1x,1pg13.6)', advance='yes' )
      end do
  end subroutine DUMP_a_HGRID

  ! ------------------------------------------------  DUMP_HGRIDS  -----
  subroutine DUMP_HGRIDS ( HGRIDS )
    use OUTPUT_M, only: OUTPUT
    type(hGrid_T), intent(in) :: HGRIDS(:)
    integer :: I
    call output ( size(hgrids), before='HGRIDS: SIZE = ', advance='yes' )
    do i = 1, size(hgrids)
      call output ( i, 4, after=': ' )
      call dump ( hgrids(i) )
    end do
  end subroutine DUMP_HGRIDS

  ! ---------------------------------------- FindClosestMatch ---
  integer function FindClosestMatch ( reference, sought, instance )
    use MLSNumerics, only: HUNT
    ! This routine is best explained in context.  Given a 'sought' quantity
    ! e.g. ptan, radiance, the profile in reference (e.g. temperature) is found
    ! that is closest to profile 'instance' in reference
    real(r8), dimension(:), intent(in) :: REFERENCE ! e.g. temperature
    real(r8), dimension(:,:), intent(in) :: SOUGHT ! e.g. ptan, radiance
    integer, intent(in) :: INSTANCE
    ! Local variables
    integer :: BESTGUESS                ! A guessed index
    integer :: FIRSTGUESS               ! A guessed index
    integer :: HIGHGUESS                ! A guessed index
    integer :: LOWGUESS                 ! A guessed index
    real(r8) :: BESTCOST                ! A cost for a guess
    real(r8) :: COST                    ! A cost for a guess

    ! Executable code
    call Hunt ( reference, sought(1,instance), firstGuess, &
      & start=max(min(instance,size(reference)),1), &
      & allowTopValue=.true., nearest=.true. )

    ! Now check the ones either side
    lowGuess = max ( firstGuess-1, 1 )
    highGuess = min ( firstGuess+1, size(reference) )
    bestCost = 0.0
    do firstGuess = lowGuess, highGuess
      cost = sum ( abs ( reference(firstGuess) - sought(:,instance) ) )
      if ( ( firstGuess == lowGuess ) .or. ( cost < bestCost ) ) then
        bestGuess = firstGuess
        bestCost = cost
      end if
    end do

    FindClosestMatch = bestGuess
  end function FindClosestMatch

  ! ---------------------------------------- NullifyHGrid -----
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
! Revision 2.4  2004/05/20 19:48:25  vsnyder
! Move Dump*HGrid here from dumper
!
! Revision 2.3  2004/05/18 01:05:06  vsnyder
! Delete unused MAFIndex and MAFCounter fields from HGrid_T
!
! Revision 2.2  2003/07/07 20:20:28  livesey
! New FindClosestMatch routine
!
! Revision 2.1  2003/06/20 19:34:45  pwagner
! Quanities now share grids stored separately in databses
!
