! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

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

  use MLSCommon, only: R8

  implicit NONE
  private
  public CreateSineTable, DestroySineTable, LogSize_SineTable_R8, SineTable_R8

  !---------------------------- RCS Ident Info -------------------------------
  character(len=*), private, parameter :: IdParm = &
    & "$Id$"
  character(len=len(idparm)), private :: Id = idParm
  character(len=*), private, parameter :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here
  !---------------------------------------------------------------------------

  interface CreateSineTable
    module procedure CreateSineTable_r8
  end interface CreateSineTable

  interface DestroySineTable
    module procedure DestroySineTable_r8
  end interface DestroySineTable

  integer, save :: LogSize_SineTable_r8 = -1

  real(r8), pointer, save :: SineTable_R8(:) => NULL()

contains

!============================================  CreateSineTable_r8  =====

  subroutine CreateSineTable_r8 ( LogSize )
    ! Initialize SineTable_r8 to be of size 2**LogSize - 1, unless it's
    ! already allocated and at least that big.

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use DFFT_M, only: InitSineTable

    integer, intent(in) :: LogSize

    if ( logSize + 2 > logSize_SineTable_r8 ) &
      & call deallocate_test ( sineTable_r8, 'SineTable_r8', moduleName )

    if ( .not. associated(sineTable_r8) ) then
      call allocate_test ( sineTable_r8, 2**logSize - 1, 'SineTable_r8', moduleName )
      call initSineTable ( sineTable_r8, logSize )
      logSize_SineTable_r8 = logSize + 2
    end if

  end subroutine CreateSineTable_r8

!==========================================   DestroySineTable_r8  =====

  subroutine DestroySineTable_r8

    use Allocate_Deallocate, only: Deallocate_Test

    call deallocate_test ( sineTable_r8, 'SineTable_r8', moduleName )
    logSize_SineTable_r8 = -1

  end subroutine DestroySineTable_r8

!=================================================  not_used_here  =====

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module SineTables_m

! $Log$
! Revision 2.1  2004/02/12 02:19:39  vsnyder
! Initial commit
!
