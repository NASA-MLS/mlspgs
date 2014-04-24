! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module HGridsDatabase                   ! Horizontal grid information
  
  use MLSKINDS, only: R8

  implicit none
  private

  public :: HGRID_T
  public :: ADDHGRIDTODATABASE, CREATEEMPTYHGRID, DESTROYHGRIDCONTENTS, &
    & DESTROYHGRIDDATABASE, DUMP, FINDCLOSESTMATCH, NULLIFYHGRID, &
    & TRIMHGRID

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! This module contains datatypes and routines for handling HGrid information
  ! HGrids are the horizontal gridding information that get into vector
  ! quantities.

  ! This is the main datatype, an HGrid.

  type HGrid_T
    integer :: Name = 0                ! String index of name.            
    integer :: masterCoordinate = 0    ! Its string index; e.g. l_phiTan
    integer :: noProfs                 ! Number of profiles in this grid  
    integer :: noProfsLowerOverlap = 0 ! Number of profiles in the lower overlap
    integer :: noProfsUpperOverlap = 0 ! Number of profiles in the upper overlap

    ! This is the maf number passing nearest to the grid point
    integer, dimension(:), pointer    :: maf         => NULL()
    ! Now the various coordinates in the HGrid, all dimensioned (1,noProfs)
    real(r8), dimension(:,:), pointer :: phi         => NULL()
    real(r8), dimension(:,:), pointer :: geodLat     => NULL()
    real(r8), dimension(:,:), pointer :: lon         => NULL()
    real(r8), dimension(:,:), pointer :: time        => NULL()
    real(r8), dimension(:,:), pointer :: solarTime   => NULL()
    real(r8), dimension(:,:), pointer :: solarZenith => NULL()
    real(r8), dimension(:,:), pointer :: losAngle    => NULL()
  end type HGrid_T

  interface DUMP
    module procedure DUMP_a_HGRID
    module procedure DUMP_HGRIDS
  end interface

contains ! =========== Public procedures ===================================

 ! -----------------------------------------  AddHGridToDatabase  -----
  integer function AddHGridToDatabase ( database, item )
    
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ALLOCATE, &
      & MLSMSG_DEALLOCATE, MLSMSG_ERROR

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

    use ALLOCATE_DEALLOCATE, only: ALLOCATE_TEST

    type (HGrid_T), intent(inout) :: HGRID

    ! Executable code
    call Allocate_Test ( hGrid%maf, hGrid%noProfs, &
      & 'hGrid%maf', ModuleName )
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

    use ALLOCATE_DEALLOCATE, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ERROR

    type (HGrid_T), intent(inout) :: HGrid
    integer, intent(in) :: SIDE         ! -1 = lower, 1 = upper
    integer, intent(in) :: NOTODELETE ! How many to delete, default all
    ! Local variables
    integer :: newNoProfs
    integer, dimension(:), pointer  :: inttemp
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

    nullify ( inttemp, temp )
    call allocate_test ( inttemp, hGrid%noProfs, 'inttemp', ModuleName )
    call allocate_test ( temp, hGrid%noProfs, 'temp', ModuleName )
    hGrid%noProfs = newNoProfs

    ! Now allocate each entry and trim it
    inttemp = hGrid%maf(:)               ! ------------------------- Maf
    call Allocate_Test ( hGrid%maf, newNoProfs, 'hGrid%maf', ModuleName)
    hGrid%maf(:) = inttemp ( first : last )
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

    use ALLOCATE_DEALLOCATE, only: DEALLOCATE_TEST

  ! This routine destroys the information associated with an hGrid

    ! Dummy arguments
    type (HGrid_T), intent(inout) :: hGrid

    ! Executable code
    
    call deallocate_test ( hGrid%maf, 'hGrid%maf', ModuleName )
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

    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_DEALLOCATE, MLSMSG_ERROR

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
    use Intrinsic, only: lit_indices
    use output_m, only: newLine, output
    use string_table, only: display_string
    type(hGrid_T), intent(in) :: aHGRID
    integer :: IERR
    integer :: J
      IERR = 0
      call output ( 'Name = ', advance='no' )
      if ( aHgrid%name > 0 ) then
        call display_string ( aHgrid%name, ierr=ierr )
        if ( ierr /= 0 ) call output ( '(not found in string table)', advance='no' )
      else
        call output('(unknown)', advance='no' )
      endif
      call output ( aHgrid%noProfs, before=' noProfs = ', advance='no' )
      call output ( aHgrid%noProfsLowerOverlap, before=' lowerOverlap = ', advance='no' )
      call output ( aHgrid%noProfsUpperOverlap, before=' upperOverlap = ', advance='yes' )
      call output ( ' prof       phi       geodLat           lon', advance='no' )
      call output ( '          time     solarTime   solarZenith', advance='no' )
      call output ( '      losAngle  nearest maf', advance='yes' )
      do j = 1, aHgrid%noProfs
        call output ( j, places=5 )
        call output ( aHgrid%phi(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%geodLat(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%lon(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%time(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%solarTime(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%solarZenith(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%losAngle(1,j), '(1x,1pg13.6)', advance='no' )
        call output ( aHgrid%maf(j), places=6, advance='yes' )
      end do
      if ( aHgrid%masterCoordinate > 0 ) then
        call output( ' Master coordinate: ', advance='no' )
        call display_string ( lit_indices(aHgrid%masterCoordinate), ierr=ierr )
        call newLine
      endif
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
    use MLSNUMERICS, only: HUNT
    ! Given a sought quantity, 
    ! the profile in reference is found
    ! that is closest to profile 'instance' in sought quantity
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
    ! Get starting Guess from Hunt
    call Hunt ( reference, sought(1,instance), firstGuess, &
      & start=max(min(instance,size(reference)),1), &
      & allowTopValue=.true., nearest=.true. )

    ! Now check for better ones either side
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

    ! Executable code.  Nothing is needed since intent(out) causes
    ! default initialization, which nullifies the pointers because
    ! they all have default initialization => NULL()

  end subroutine NullifyHGrid

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module HGridsDatabase

! $Log$
! Revision 2.13  2014/04/24 23:50:25  pwagner
! Added masterCoordinate component
!
! Revision 2.12  2013/10/01 22:16:45  pwagner
! Added maf component to HGrid_T
!
! Revision 2.11  2013/08/12 23:47:07  vsnyder
! Default initialize more stuff in HGrid_t
!
! Revision 2.10  2013/03/01 01:04:51  pwagner
! Get R8 from MLSKinds
!
! Revision 2.9  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.8  2008/05/02 00:38:48  vsnyder
! Simplify NullifyHGrid
!
! Revision 2.7  2005/09/21 23:11:54  pwagner
! Unnecessary changes
!
! Revision 2.6  2005/08/25 20:18:05  pwagner
! Protect against crashing when dumping anonymous HGrids
!
! Revision 2.5  2005/06/22 17:25:49  pwagner
! Reworded Copyright statement, moved rcs id
!
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
