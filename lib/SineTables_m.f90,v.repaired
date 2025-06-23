! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module SineTables_m

  ! The ?FFT_M modules assume that only one sine table exists.  They keep
  ! track of the maximum it's been, and assume that any sine table they
  ! get is at least that size, unless they're told to reinitialize it.
  ! To avoid reinitializing sine tables that are smaller than the biggest
  ! one they've seen (which we can't know, because the variables they use
  ! for that purpose are private), the sine table and its size are here.

  ! It is intended that one will call CreateSineTable, and then use
  ! SineTable_R* as the S argument, and LogSize_SineTable_R* as the MS
  ! argument, in calls to the procedures in the ?FFT_M modules.

  use MLSKinds, only: R8

  implicit NONE
  private
  public CreateSineTable, DestroySineTable, LogSize_SineTable_R8, SineTable_R8

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  interface CreateSineTable
    module procedure CreateSineTable_r8
  end interface CreateSineTable

  interface DestroySineTable
    module procedure DestroySineTable_r8
  end interface DestroySineTable

  integer, save :: LogSize_SineTable_r8 = -1

  real(r8), allocatable, save :: SineTable_R8(:)

contains

!============================================  CreateSineTable_r8  =====

  subroutine CreateSineTable_r8 ( LogSize )
    ! Initialize SineTable_r8 to be of size 2**LogSize - 1, unless it's
    ! already allocated and at least that big.

    use Allocate_Deallocate, only: Byte_Size, Bytes, &
      & Test_Allocate, Test_Deallocate
    use DFFT_M, only: InitSineTable

    integer, intent(in) :: LogSize

    integer :: N, Stat

    if ( allocated(sineTable_r8) .and. logSize + 2 > logSize_SineTable_r8 ) then
      n = byte_size(sineTable_r8)
      deallocate ( sineTable_r8, stat=stat )
      call test_deallocate ( stat, moduleName, 'SineTable_r8', n )
    end if

    if ( .not. allocated(sineTable_r8) ) then
      n = 2**logSize - 1
      allocate ( sineTable_r8(n), stat=stat )
      call test_allocate ( stat, moduleName, 'SineTable_r8', (/1/), (/n/), &
        & bytes(sineTable_r8) )
      call initSineTable ( sineTable_r8, logSize )
      logSize_SineTable_r8 = logSize + 2
    end if

  end subroutine CreateSineTable_r8

!==========================================   DestroySineTable_r8  =====

  subroutine DestroySineTable_r8

    use Allocate_Deallocate, only: Byte_Size, Test_Deallocate

    integer :: N, Stat

    if ( allocated(sineTable_r8) ) then
      n = byte_size(sineTable_r8)
      deallocate ( sineTable_r8, stat=stat )
      call test_deallocate ( stat, moduleName, 'SineTable_r8', n )
    end if
    logSize_SineTable_r8 = -1

  end subroutine DestroySineTable_r8

!=================================================  not_used_here  =====

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module SineTables_m

! $Log$
! Revision 2.6  2014/09/05 00:18:06  vsnyder
! Get kinds from MLSKinds instead of MLSCommon.  Keep track of size in
! bytes instead of Memory_Units.
!
! Revision 2.5  2013/06/12 02:17:27  vsnyder
! UBYTES and BYTE_SIZE from Allocate_Deallocate
!
! Revision 2.4  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.3  2007/04/14 00:38:25  vsnyder
! Change SineTable from pointer to allocatable, improves FFT performance
!
! Revision 2.2  2005/06/22 17:25:50  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.1  2004/02/12 02:19:39  vsnyder
! Initial commit
!
