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

  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error, MLSMSG_Warning

  implicit NONE
  private

  public :: ALLOCATE_TEST, DEALLOCATE_TEST, DEALLOC_STATUS

  integer, save :: DEALLOC_STATUS = 0

  interface ALLOCATE_TEST
    module procedure ALLOCATE_TEST_CHARACTER_1D
    module procedure ALLOCATE_TEST_DOUBLE_1D, ALLOCATE_TEST_DOUBLE_2D
    module procedure ALLOCATE_TEST_DOUBLE_3D
    module procedure ALLOCATE_TEST_INTEGER_1D, ALLOCATE_TEST_INTEGER_2D
    module procedure ALLOCATE_TEST_REAL_1D, ALLOCATE_TEST_REAL_2D
  end interface

  interface DEALLOCATE_TEST
    module procedure DEALLOCATE_TEST_CHARACTER_1D
    module procedure DEALLOCATE_TEST_DOUBLE_1D, DEALLOCATE_TEST_DOUBLE_2D
    module procedure DEALLOCATE_TEST_DOUBLE_3D
    module procedure DEALLOCATE_TEST_INTEGER_1D, DEALLOCATE_TEST_INTEGER_2D
    module procedure DEALLOCATE_TEST_REAL_1D, DEALLOCATE_TEST_REAL_2D
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
  ! ------------------------------------  Allocate_Test_Double_1d  -----
  subroutine Allocate_Test_Double_1d ( To_Allocate, Dim1, Its_Name, &
    & ModuleName, LowBound )
    double precision, pointer, dimension(:) :: To_Allocate
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
  end subroutine Allocate_Test_Double_1d
  ! ------------------------------------  Allocate_Test_Double_2d  -----
  subroutine Allocate_Test_Double_2d ( To_Allocate, Dim1, Dim2, Its_Name, &
    & ModuleName )
    double precision, pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! First dimension of To_Allocate
    integer, intent(in) :: Dim2    ! Second dimension of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer :: STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    allocate ( To_Allocate(dim1,dim2), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name )
  end subroutine Allocate_Test_Double_2d
  ! ------------------------------------  Allocate_Test_Double_3d  -----
  subroutine Allocate_Test_Double_3d ( To_Allocate, Dim1, Dim2, Dim3, &
    & Its_Name, ModuleName )
    double precision, pointer, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! First dimension of To_Allocate
    integer, intent(in) :: Dim2    ! Second dimension of To_Allocate
    integer, intent(in) :: Dim3    ! Third dimension of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer :: STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    allocate ( To_Allocate(dim1,dim2,dim3), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name )
  end subroutine Allocate_Test_Double_3d
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
  ! --------------------------------------  Allocate_Test_Real_1d  -----
  subroutine Allocate_Test_Real_1d ( To_Allocate, Dim1, Its_Name, ModuleName, &
    & LowBound )
    real, pointer, dimension(:) :: To_Allocate
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
  end subroutine Allocate_Test_Real_1d
  ! --------------------------------------  Allocate_Test_Real_2d  -----
  subroutine Allocate_Test_Real_2d ( To_Allocate, Dim1, Dim2, Its_Name, &
    & ModuleName )
    real, pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! First dimension of To_Allocate
    integer, intent(in) :: Dim2    ! Second dimension of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer :: STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    allocate ( To_Allocate(dim1,dim2), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name )
  end subroutine Allocate_Test_Real_2d

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
      end if
    end if
  end subroutine Deallocate_Test_Character_1d
  ! ----------------------------------  Deallocate_Test_Double_1d  -----
  subroutine Deallocate_Test_Double_1d ( To_Deallocate, Its_Name, ModuleName )
    double precision, pointer, dimension(:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    if ( associated(To_Deallocate) ) then
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      end if
    end if
  end subroutine Deallocate_Test_Double_1d
  ! ----------------------------------  Deallocate_Test_Double_2d  -----
  subroutine Deallocate_Test_Double_2d ( To_Deallocate, Its_Name, ModuleName )
    double precision, pointer, dimension(:,:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    if ( associated(To_Deallocate) ) then
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      end if
    end if
  end subroutine Deallocate_Test_Double_2d
  ! ----------------------------------  Deallocate_Test_Double_3d  -----
  subroutine Deallocate_Test_Double_3d ( To_Deallocate, Its_Name, ModuleName )
    double precision, pointer, dimension(:,:,:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    if ( associated(To_Deallocate) ) then
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      end if
    end if
  end subroutine Deallocate_Test_Double_3d
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
      end if
    end if
  end subroutine Deallocate_Test_Integer_2d
  ! ------------------------------------  Deallocate_Test_Real_1d  -----
  subroutine Deallocate_Test_Real_1d ( To_Deallocate, Its_Name, ModuleName )
    real, pointer, dimension(:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    if ( associated(To_Deallocate) ) then
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      end if
    end if
  end subroutine Deallocate_Test_Real_1d
  ! ------------------------------------  Deallocate_Test_Real_2d  -----
  subroutine Deallocate_Test_Real_2d ( To_Deallocate, Its_Name, ModuleName )
    real, pointer, dimension(:,:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    if ( associated(To_Deallocate) ) then
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      end if
    end if
  end subroutine Deallocate_Test_Real_2d

end module Allocate_Deallocate

! $Log$
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
