! Copyright (c) 1999, California Institute of Technology. ALL RIGHTS RESERVED.
! U.S. Government sponsorship under NASA Contract NAS7407 is acknowledged.

!=============================================================================
module Allocate_Deallocate
!=============================================================================

! This module contains procedures to allocate, test the allocation,
! and announce an error if it fails, and similarly for deallocation,
! for Double precision, Integer and Real arrays of one or two dimensions.

! **************************************************
! *****     Important Notice:                  *****
! **************************************************
! *****     All of the specific procedures     *****
! *****     of the Allocate_test generic       *****
! *****     deallocate the object to be        *****
! *****     allocated.  Therefore, it must     *****
! *****     be declared with => NULL()         *****
! **************************************************

  use MACHINE, only: MLS_GC_NOW
  use MLSCommon, only: r4, r8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error, MLSMSG_Warning

  implicit NONE
  private

  public :: ALLOCATE_TEST, DEALLOCATE_TEST, DEALLOC_STATUS, & 
    & SET_GARBAGE_COLLECTION

  integer, save :: DEALLOC_STATUS = 0
  logical, save :: COLLECT_GARBAGE_EACH_TIME = .false.

  interface ALLOCATE_TEST
    module procedure ALLOCATE_TEST_CHARACTER_1D
    module procedure ALLOCATE_TEST_CHARACTER_2D
    module procedure ALLOCATE_TEST_INTEGER_1D, ALLOCATE_TEST_INTEGER_2D
    module procedure ALLOCATE_TEST_INTEGER_3D
    module procedure ALLOCATE_TEST_LOGICAL_1D, ALLOCATE_TEST_LOGICAL_2D
    module procedure ALLOCATE_TEST_REALR4_1D, ALLOCATE_TEST_REALR4_2D
    module procedure ALLOCATE_TEST_REALR4_3D
    module procedure ALLOCATE_TEST_REALR8_1D, ALLOCATE_TEST_REALR8_2D
    module procedure ALLOCATE_TEST_REALR8_3D
  end interface

  interface DEALLOCATE_TEST
    module procedure DEALLOCATE_TEST_CHARACTER_1D
    module procedure DEALLOCATE_TEST_CHARACTER_2D
    module procedure DEALLOCATE_TEST_INTEGER_1D, DEALLOCATE_TEST_INTEGER_2D
    module procedure DEALLOCATE_TEST_INTEGER_3D
    module procedure DEALLOCATE_TEST_LOGICAL_1D, DEALLOCATE_TEST_LOGICAL_2D
    module procedure DEALLOCATE_TEST_REALR4_1D, DEALLOCATE_TEST_REALR4_2D
    module procedure DEALLOCATE_TEST_REALR4_3D
    module procedure DEALLOCATE_TEST_REALR8_1D, DEALLOCATE_TEST_REALR8_2D
    module procedure DEALLOCATE_TEST_REALR8_3D
  end interface

  !------------------------------- RCS Ident Info ------------------------------
  character(len=130), private :: id = & 
       "$Id$"
  character(len=*), parameter, private :: ModuleName = &
    & "$RCSfile$"
  !-----------------------------------------------------------------------------

contains ! =====     Private Procedures     ============================
  ! ---------------------------------  Allocate_Test_Character_1d  -----
  subroutine Allocate_Test_Character_1d ( To_Allocate, Dim1, Its_Name, &
    & ModuleName, LowBound )
    character(len=*), pointer, dimension(:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer, intent(in), optional :: LowBound     ! Lower bound, default 1
    integer :: MY_LOW, STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    my_low = 1
    if ( present(lowBound) ) my_low = lowBound
    allocate ( To_Allocate(my_low:dim1), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name )
  end subroutine Allocate_Test_Character_1d
  ! ---------------------------------  Allocate_Test_Character_2d  -----
  subroutine Allocate_Test_Character_2d ( To_Allocate, Dim1, Dim2, Its_Name, &
    & ModuleName )
    character(len=*), pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Second dimension of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer :: STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    allocate ( To_Allocate(dim1,dim2), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name )
  end subroutine Allocate_Test_Character_2d
  ! ------------------------------------  Allocate_Test_RealR8_1d  -----
  subroutine Allocate_Test_RealR8_1d ( To_Allocate, Dim1, Its_Name, &
    & ModuleName, LowBound )
    real (r8), pointer, dimension(:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer, intent(in), optional :: LowBound     ! Lower bound, default 1
    integer :: MY_LOW, STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    my_low = 1
    if ( present(lowBound) ) my_low = lowBound
    allocate ( To_Allocate(my_low:dim1), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name )
  end subroutine Allocate_Test_RealR8_1d
  ! ------------------------------------  Allocate_Test_RealR8_2d  -----
  subroutine Allocate_Test_RealR8_2d ( To_Allocate, Dim1, Dim2, Its_Name, &
    & ModuleName )
    real (r8), pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! First dimension of To_Allocate
    integer, intent(in) :: Dim2    ! Second dimension of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer :: STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    allocate ( To_Allocate(dim1,dim2), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name )
  end subroutine Allocate_Test_RealR8_2d
  ! ------------------------------------  Allocate_Test_RealR8_3d  -----
  subroutine Allocate_Test_RealR8_3d ( To_Allocate, Dim1, Dim2, Dim3, &
    & Its_Name, ModuleName )
    real (r8), pointer, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! First dimension of To_Allocate
    integer, intent(in) :: Dim2    ! Second dimension of To_Allocate
    integer, intent(in) :: Dim3    ! Third dimension of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer :: STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    allocate ( To_Allocate(dim1,dim2,dim3), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name )
  end subroutine Allocate_Test_RealR8_3d
  ! -----------------------------------  Allocate_Test_Integer_1d  -----
  subroutine Allocate_Test_Integer_1d ( To_Allocate, Dim1, Its_Name, &
    & ModuleName, LowBound )
    integer, pointer, dimension(:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer, intent(in), optional :: LowBound     ! Lower bound, default 1
    integer :: MY_LOW, STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    my_low = 1
    if ( present(lowBound) ) my_low = lowBound
    allocate ( To_Allocate(my_low:dim1), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name )
  end subroutine Allocate_Test_Integer_1d
  ! -----------------------------------  Allocate_Test_Integer_2d  -----
  subroutine Allocate_Test_Integer_2d ( To_Allocate, Dim1, Dim2, Its_Name, &
    & ModuleName )
    integer, pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! First dimension of To_Allocate
    integer, intent(in) :: Dim2    ! Second dimension of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer :: STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    allocate ( To_Allocate(dim1,dim2), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name )
  end subroutine Allocate_Test_Integer_2d
  ! -----------------------------------  Allocate_Test_Integer_3d  -----
  subroutine Allocate_Test_Integer_3d ( To_Allocate, Dim1, Dim2, Dim3, Its_Name, &
    & ModuleName )
    integer, pointer, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! First dimension of To_Allocate
    integer, intent(in) :: Dim2    ! Second dimension of To_Allocate
    integer, intent(in) :: Dim3    ! Third dimension of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer :: STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    allocate ( To_Allocate(dim1,dim2,dim3), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name )
  end subroutine Allocate_Test_Integer_3d
  ! -----------------------------------  Allocate_Test_Logical_1d  -----
  subroutine Allocate_Test_Logical_1d ( To_Allocate, Dim1, Its_Name, &
    & ModuleName, LowBound )
    logical, pointer, dimension(:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer, intent(in), optional :: LowBound     ! Lower bound, default 1
    integer :: MY_LOW, STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    my_low = 1
    if ( present(lowBound) ) my_low = lowBound
    allocate ( To_Allocate(my_low:dim1), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name )
  end subroutine Allocate_Test_Logical_1d
  ! --------------------------------------  Allocate_Test_RealR4_2d  -----
  subroutine Allocate_Test_Logical_2d ( To_Allocate, Dim1, Dim2, Its_Name, &
    & ModuleName )
    logical, pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! First dimension of To_Allocate
    integer, intent(in) :: Dim2    ! Second dimension of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer :: STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    allocate ( To_Allocate(dim1,dim2), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name )
  end subroutine Allocate_Test_Logical_2d
  ! --------------------------------------  Allocate_Test_RealR4_1d  -----
  subroutine Allocate_Test_RealR4_1d ( To_Allocate, Dim1, Its_Name, ModuleName, &
    & LowBound )
    real (r4), pointer, dimension(:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer, intent(in), optional :: LowBound     ! Lower bound, default 1
    integer :: MY_LOW, STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    my_low = 1
    if ( present(lowBound) ) my_low = lowBound
    allocate ( To_Allocate(my_low:dim1), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name )
  end subroutine Allocate_Test_RealR4_1d
  ! --------------------------------------  Allocate_Test_RealR4_2d  -----
  subroutine Allocate_Test_RealR4_2d ( To_Allocate, Dim1, Dim2, Its_Name, &
    & ModuleName )
    real (r4), pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! First dimension of To_Allocate
    integer, intent(in) :: Dim2    ! Second dimension of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer :: STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    allocate ( To_Allocate(dim1,dim2), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name )
  end subroutine Allocate_Test_RealR4_2d
  ! ------------------------------------  Allocate_Test_RealR4_3d  -----
  subroutine Allocate_Test_RealR4_3d ( To_Allocate, Dim1, Dim2, Dim3, &
    & Its_Name, ModuleName )
    real (r4), pointer, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! First dimension of To_Allocate
    integer, intent(in) :: Dim2    ! Second dimension of To_Allocate
    integer, intent(in) :: Dim3    ! Third dimension of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer :: STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    allocate ( To_Allocate(dim1,dim2,dim3), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name )
  end subroutine Allocate_Test_RealR4_3d
  ! -------------------------------  Deallocate_Test_Character_1d  -----
  subroutine Deallocate_Test_Character_1d ( To_Deallocate, Its_Name, ModuleName )
    character(len=*), pointer, dimension(:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    if ( associated(To_Deallocate) ) then
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      elseif ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
    end if
  end subroutine Deallocate_Test_Character_1d
  ! -------------------------------  Deallocate_Test_Character_2d  -----
  subroutine Deallocate_Test_Character_2d ( To_Deallocate, Its_Name, ModuleName )
    character(len=*), pointer, dimension(:,:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    if ( associated(To_Deallocate) ) then
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      elseif ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
    end if
  end subroutine Deallocate_Test_Character_2d
  ! ----------------------------------  Deallocate_Test_RealR8_1d  -----
  subroutine Deallocate_Test_RealR8_1d ( To_Deallocate, Its_Name, ModuleName )
    real (r8), pointer, dimension(:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    if ( associated(To_Deallocate) ) then
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      elseif ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
    end if
  end subroutine Deallocate_Test_RealR8_1d
  ! ----------------------------------  Deallocate_Test_RealR8_2d  -----
  subroutine Deallocate_Test_RealR8_2d ( To_Deallocate, Its_Name, ModuleName )
    real (r8), pointer, dimension(:,:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    if ( associated(To_Deallocate) ) then
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      elseif ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
    end if
  end subroutine Deallocate_Test_RealR8_2d
  ! ----------------------------------  Deallocate_Test_RealR8_3d  -----
  subroutine Deallocate_Test_RealR8_3d ( To_Deallocate, Its_Name, ModuleName )
    real (r8), pointer, dimension(:,:,:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    if ( associated(To_Deallocate) ) then
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      elseif ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
    end if
  end subroutine Deallocate_Test_RealR8_3d
  ! ---------------------------------  Deallocate_Test_Integer_1d  -----
  subroutine Deallocate_Test_Integer_1d ( To_Deallocate, Its_Name, ModuleName )
    integer, pointer, dimension(:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    if ( associated(To_Deallocate) ) then
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      elseif ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
    end if
  end subroutine Deallocate_Test_Integer_1d
  ! ---------------------------------  Deallocate_Test_Integer_2d  -----
  subroutine Deallocate_Test_Integer_2d ( To_Deallocate, Its_Name, ModuleName )
    integer, pointer, dimension(:,:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    if ( associated(To_Deallocate) ) then
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      elseif ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
    end if
  end subroutine Deallocate_Test_Integer_2d
  ! ---------------------------------  Deallocate_Test_Integer_3d  -----
  subroutine Deallocate_Test_Integer_3d ( To_Deallocate, Its_Name, ModuleName )
    integer, pointer, dimension(:,:,:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    if ( associated(To_Deallocate) ) then
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      elseif ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
    end if
  end subroutine Deallocate_Test_Integer_3d
  ! ---------------------------------  Deallocate_Test_Logical_1d  -----
  subroutine Deallocate_Test_Logical_1d ( To_Deallocate, Its_Name, ModuleName )
    logical, pointer, dimension(:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    if ( associated(To_Deallocate) ) then
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      elseif ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
    end if
  end subroutine Deallocate_Test_Logical_1d
  ! ---------------------------------  Deallocate_Test_Logical_1d  -----
  subroutine Deallocate_Test_Logical_2d ( To_Deallocate, Its_Name, ModuleName )
    logical, pointer, dimension(:,:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    if ( associated(To_Deallocate) ) then
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      elseif ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
    end if
  end subroutine Deallocate_Test_Logical_2d
  ! ------------------------------------  Deallocate_Test_RealR4_1d  -----
  subroutine Deallocate_Test_RealR4_1d ( To_Deallocate, Its_Name, ModuleName )
    real (r4), pointer, dimension(:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    if ( associated(To_Deallocate) ) then
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      elseif ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
    end if
  end subroutine Deallocate_Test_RealR4_1d
  ! ------------------------------------  Deallocate_Test_RealR4_2d  -----
  subroutine Deallocate_Test_RealR4_2d ( To_Deallocate, Its_Name, ModuleName )
    real (r4), pointer, dimension(:,:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    if ( associated(To_Deallocate) ) then
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      elseif ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
    end if
  end subroutine Deallocate_Test_RealR4_2d
  ! ----------------------------------  Deallocate_Test_RealR4_3d  -----
  subroutine Deallocate_Test_RealR4_3d ( To_Deallocate, Its_Name, ModuleName )
    real (r4), pointer, dimension(:,:,:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    if ( associated(To_Deallocate) ) then
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      elseif ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
    end if
  end subroutine Deallocate_Test_RealR4_3d

  ! ----------------------------------  Set_garbage_collection  -----
  subroutine Set_garbage_collection ( setting )
    logical :: setting
    collect_garbage_each_time = setting
  end subroutine Set_garbage_collection

end module Allocate_Deallocate

! $Log$
! Revision 2.9  2002/02/05 00:42:47  pwagner
! Optionally collects garbage after each deallocate_test
!
! Revision 2.8  2002/02/01 01:48:17  vsnyder
! Add [De]Allocate_Test_Character_2d
!
! Revision 2.7  2001/09/10 21:05:01  livesey
! Added logical 2d stuff
!
! Revision 2.6  2001/05/30 23:52:15  livesey
! Added 3D integer option.
!
! Revision 2.5  2001/04/02 20:50:27  vsnyder
! Add logical 1-D
!
! Revision 2.4  2001/03/24 23:25:57  livesey
! Added real(r4) 3d, and renamed routines
!
! Revision 2.3  2001/02/22 01:54:41  vsnyder
! Periodic commit
!
! Revision 2.2  2000/10/09 23:03:49  vsnyder
! Provide for lower bounds for allocated arrays
!
! Revision 2.0  2000/09/05 18:57:02  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:03  vsnyder
! Initial entry
!
