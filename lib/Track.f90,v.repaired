! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Track_m

  implicit None

  private

  public :: TrackAllocate, TrackDeallocate, ReportLeaks

  interface TrackAllocate
    module procedure TrackAllocate_C1
    module procedure TrackAllocate_C2
    module procedure TrackAllocate_C3
    module procedure TrackAllocate_X1
    module procedure TrackAllocate_X2
    module procedure TrackAllocate_X3
    module procedure TrackAllocate_Z1
    module procedure TrackAllocate_Z2
    module procedure TrackAllocate_Z3
    module procedure TrackAllocate_I1
    module procedure TrackAllocate_I2
    module procedure TrackAllocate_I3
    module procedure TrackAllocate_I4
    module procedure TrackAllocate_R1
    module procedure TrackAllocate_R2
    module procedure TrackAllocate_R3
    module procedure TrackAllocate_R4
    module procedure TrackAllocate_D1
    module procedure TrackAllocate_D2
    module procedure TrackAllocate_D3
    module procedure TrackAllocate_D4
    module procedure TrackAllocate_L1
    module procedure TrackAllocate_L2
    module procedure TrackAllocate_L3
  ! For allocatable argument instead of pointer
    module procedure TrackAllocateA_C1
    module procedure TrackAllocateA_C2
    module procedure TrackAllocateA_I1
    module procedure TrackAllocateA_I2
    module procedure TrackAllocateA_I3
    module procedure TrackAllocateA_I4
    module procedure TrackAllocateA_L1
    module procedure TrackAllocateA_L2
    module procedure TrackAllocateA_L3
    module procedure TrackAllocateA_R1
    module procedure TrackAllocateA_R2
    module procedure TrackAllocateA_R3
    module procedure TrackAllocateA_R4
    module procedure TrackAllocateA_D1
    module procedure TrackAllocateA_D2
    module procedure TrackAllocateA_D3
    module procedure TrackAllocateA_D4
  end interface

  interface TrackDeallocate
    module procedure TrackDeallocate_C1
    module procedure TrackDeallocate_C2
    module procedure TrackDeallocate_C3
    module procedure TrackDeallocate_X1
    module procedure TrackDeallocate_X2
    module procedure TrackDeallocate_X3
    module procedure TrackDeallocate_Z1
    module procedure TrackDeallocate_Z2
    module procedure TrackDeallocate_Z3
    module procedure TrackDeallocate_I1
    module procedure TrackDeallocate_I2
    module procedure TrackDeallocate_I3
    module procedure TrackDeallocate_I4
    module procedure TrackDeallocate_R1
    module procedure TrackDeallocate_R2
    module procedure TrackDeallocate_R3
    module procedure TrackDeallocate_R4
    module procedure TrackDeallocate_D1
    module procedure TrackDeallocate_D2
    module procedure TrackDeallocate_D3
    module procedure TrackDeallocate_D4
    module procedure TrackDeallocate_L1
    module procedure TrackDeallocate_L2
    module procedure TrackDeallocate_L3
  ! For allocatable argument instead of pointer
    module procedure TrackDeallocateA_C1
    module procedure TrackDeallocateA_C2
    module procedure TrackDeallocateA_I1
    module procedure TrackDeallocateA_I2
    module procedure TrackDeallocateA_I3
    module procedure TrackDeallocateA_I4
    module procedure TrackDeallocateA_L1
    module procedure TrackDeallocateA_L2
    module procedure TrackDeallocateA_L3
    module procedure TrackDeallocateA_R1
    module procedure TrackDeallocateA_R2
    module procedure TrackDeallocateA_R3
    module procedure TrackDeallocateA_R4
    module procedure TrackDeallocateA_D1
    module procedure TrackDeallocateA_D2
    module procedure TrackDeallocateA_D3
    module procedure TrackDeallocateA_D4
  end interface

  integer, parameter :: WhereLen = 48

  type :: Track_C1_t
    character, pointer :: P(:) => NULL()
    character(len=whereLen) :: Where
  end type Track_C1_t
  type :: Track_C2_t
    character, pointer :: P(:,:) => NULL()
    character(len=whereLen) :: Where
  end type Track_C2_t
  type :: Track_C3_t
    character, pointer :: P(:,:,:) => NULL()
    character(len=whereLen) :: Where
  end type Track_C3_t
  type :: Track_X1_t
    complex, pointer :: P(:) => NULL()
    character(len=whereLen) :: Where
  end type Track_X1_t
  type :: Track_X2_t
    complex, pointer :: P(:,:) => NULL()
    character(len=whereLen) :: Where
  end type Track_X2_t
  type :: Track_X3_t
    complex, pointer :: P(:,:,:) => NULL()
    character(len=whereLen) :: Where
  end type Track_X3_t
  type :: Track_Z1_t
    complex(kind(0.0d0)), pointer :: P(:) => NULL()
    character(len=whereLen) :: Where
  end type Track_Z1_t
  type :: Track_Z2_t
    complex(kind(0.0d0)), pointer :: P(:,:) => NULL()
    character(len=whereLen) :: Where
  end type Track_Z2_t
  type :: Track_Z3_t
    complex(kind(0.0d0)), pointer :: P(:,:,:) => NULL()
    character(len=whereLen) :: Where
  end type Track_Z3_t
  type :: Track_I1_t
    integer, pointer :: P(:) => NULL()
    character(len=whereLen) :: Where
  end type Track_I1_t
  type :: Track_I2_t
    integer, pointer :: P(:,:) => NULL()
    character(len=whereLen) :: Where
  end type Track_I2_t
  type :: Track_I3_t
    integer, pointer :: P(:,:,:) => NULL()
    character(len=whereLen) :: Where
  end type Track_I3_t
  type :: Track_L1_t
    logical, pointer :: P(:) => NULL()
    character(len=whereLen) :: Where
  end type Track_L1_t
  type :: Track_L2_t
    logical, pointer :: P(:,:) => NULL()
    character(len=whereLen) :: Where
  end type Track_L2_t
  type :: Track_L3_t
    logical, pointer :: P(:,:,:) => NULL()
    character(len=whereLen) :: Where
  end type Track_L3_t
  type :: Track_I4_t
    integer, pointer :: P(:,:,:,:) => NULL()
    character(len=whereLen) :: Where
  end type Track_I4_t
  type :: Track_R1_t
    real, pointer :: P(:) => NULL()
    character(len=whereLen) :: Where
  end type Track_R1_t
  type :: Track_R2_t
    real, pointer :: P(:,:) => NULL()
    character(len=whereLen) :: Where
  end type Track_R2_t
  type :: Track_R3_t
    real, pointer :: P(:,:,:) => NULL()
    character(len=whereLen) :: Where
  end type Track_R3_t
  type :: Track_R4_t
    real, pointer :: P(:,:,:,:) => NULL()
    character(len=whereLen) :: Where
  end type Track_R4_t
  type :: Track_D1_t
    double precision, pointer :: P(:) => NULL()
    character(len=whereLen) :: Where
  end type Track_D1_t
  type :: Track_D2_t
    double precision, pointer :: P(:,:) => NULL()
    character(len=whereLen) :: Where
  end type Track_D2_t
  type :: Track_D3_t
    double precision, pointer :: P(:,:,:) => NULL()
    character(len=whereLen) :: Where
  end type Track_D3_t
  type :: Track_D4_t
    double precision, pointer :: P(:,:,:,:) => NULL()
    character(len=whereLen) :: Where
  end type Track_D4_t

  type(Track_C1_t), pointer, save :: Track_C1(:) => NULL()
  type(Track_C2_t), pointer, save :: Track_C2(:) => NULL()
  type(Track_C3_t), pointer, save :: Track_C3(:) => NULL()
  type(Track_X1_t), pointer, save :: Track_X1(:) => NULL()
  type(Track_X2_t), pointer, save :: Track_X2(:) => NULL()
  type(Track_X3_t), pointer, save :: Track_X3(:) => NULL()
  type(Track_Z1_t), pointer, save :: Track_Z1(:) => NULL()
  type(Track_Z2_t), pointer, save :: Track_Z2(:) => NULL()
  type(Track_Z3_t), pointer, save :: Track_Z3(:) => NULL()
  type(Track_I1_t), pointer, save :: Track_I1(:) => NULL()
  type(Track_I2_t), pointer, save :: Track_I2(:) => NULL()
  type(Track_I3_t), pointer, save :: Track_I3(:) => NULL()
  type(Track_I4_t), pointer, save :: Track_I4(:) => NULL()
  type(Track_R1_t), pointer, save :: Track_R1(:) => NULL()
  type(Track_R2_t), pointer, save :: Track_R2(:) => NULL()
  type(Track_R3_t), pointer, save :: Track_R3(:) => NULL()
  type(Track_R4_t), pointer, save :: Track_R4(:) => NULL()
  type(Track_D1_t), pointer, save :: Track_D1(:) => NULL()
  type(Track_D2_t), pointer, save :: Track_D2(:) => NULL()
  type(Track_D3_t), pointer, save :: Track_D3(:) => NULL()
  type(Track_D4_t), pointer, save :: Track_D4(:) => NULL()
  type(Track_L1_t), pointer, save :: Track_L1(:) => NULL()
  type(Track_L2_t), pointer, save :: Track_L2(:) => NULL()
  type(Track_L3_t), pointer, save :: Track_L3(:) => NULL()

  integer, save :: Num_C1 = 0
  integer, save :: Num_C2 = 0
  integer, save :: Num_C3 = 0
  integer, save :: Num_X1 = 0
  integer, save :: Num_X2 = 0
  integer, save :: Num_X3 = 0
  integer, save :: Num_Z1 = 0
  integer, save :: Num_Z2 = 0
  integer, save :: Num_Z3 = 0
  integer, save :: Num_I1 = 0
  integer, save :: Num_I2 = 0
  integer, save :: Num_I3 = 0
  integer, save :: Num_I4 = 0
  integer, save :: Num_R1 = 0
  integer, save :: Num_R2 = 0
  integer, save :: Num_R3 = 0
  integer, save :: Num_R4 = 0
  integer, save :: Num_D1 = 0
  integer, save :: Num_D2 = 0
  integer, save :: Num_D3 = 0
  integer, save :: Num_D4 = 0
  integer, save :: Num_L1 = 0
  integer, save :: Num_L2 = 0
  integer, save :: Num_L3 = 0

  integer, parameter :: InitSize = 100

  integer, parameter :: CharsPerInt = 4 ! Characters per integer

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine TrackAllocate_C1 ( What, Where, Module )
    character(len=*), pointer :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_c1_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_c1
      if ( .not. associated(track_c1(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_c1) ) then
      allocate(track_c1(initSize))
    else if ( i > size(track_c1) ) then
      temp => track_c1
      allocate ( track_c1(2*size(temp)) )
      track_c1(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_c1(i)%p => what(:)(1:1)
    track_c1(i)%where = trim(where) // "@" // trimModule(module)
    num_c1 = max(num_c1,i)
  end subroutine TrackAllocate_C1
  subroutine TrackAllocate_C2 ( What, Where, Module )
    character(len=*), pointer :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_c2_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_c2
      if ( .not. associated(track_c2(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_c2) ) then
      allocate(track_c2(initSize))
    else if ( i > size(track_c2) ) then
      temp => track_c2
      allocate ( track_c2(2*size(temp)) )
      track_c2(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_c2(i)%p => what(:,:)(1:1)
    track_c2(i)%where = trim(where) // "@" // trimModule(module)
    num_c2 = max(num_c2,i)
  end subroutine TrackAllocate_C2
  subroutine TrackAllocate_C3 ( What, Where, Module )
    character(len=*), pointer :: What(:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_c3_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_c3
      if ( .not. associated(track_c3(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_c3) ) then
      allocate(track_c3(initSize))
    else if ( i > size(track_c3) ) then
      temp => track_c3
      allocate ( track_c3(2*size(temp)) )
      track_c3(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_c3(i)%p => what(:,:,:)(1:1)
    track_c3(i)%where = trim(where) // "@" // trimModule(module)
    num_c3 = max(num_c3,i)
  end subroutine TrackAllocate_C3

  subroutine TrackAllocate_X1 ( What, Where, Module )
    complex, pointer :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_x1_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_x1
      if ( .not. associated(track_x1(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_x1) ) then
      allocate(track_x1(initSize))
    else if ( i > size(track_x1) ) then
      temp => track_x1
      allocate ( track_x1(2*size(temp)) )
      track_x1(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_x1(i)%p => what(:)
    track_x1(i)%where = trim(where) // "@" // trimModule(module)
    num_x1 = max(num_x1,i)
  end subroutine TrackAllocate_X1
  subroutine TrackAllocate_X2 ( What, Where, Module )
    complex, pointer :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_x2_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_x2
      if ( .not. associated(track_x2(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_x2) ) then
      allocate(track_x2(initSize))
    else if ( i > size(track_x2) ) then
      temp => track_x2
      allocate ( track_x2(2*size(temp)) )
      track_x2(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_x2(i)%p => what(:,:)
    track_x2(i)%where = trim(where) // "@" // trimModule(module)
    num_x2 = max(num_x2,i)
  end subroutine TrackAllocate_X2
  subroutine TrackAllocate_X3 ( What, Where, Module )
    complex, pointer :: What(:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_x3_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_x3
      if ( .not. associated(track_x3(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_x3) ) then
      allocate(track_x3(initSize))
    else if ( i > size(track_x3) ) then
      temp => track_x3
      allocate ( track_x3(2*size(temp)) )
      track_x3(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_x3(i)%p => what(:,:,:)
    track_x3(i)%where = trim(where) // "@" // trimModule(module)
    num_x3 = max(num_x3,i)
  end subroutine TrackAllocate_X3

  subroutine TrackAllocate_Z1 ( What, Where, Module )
    complex(kind(0.0d0)), pointer :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_z1_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_z1
      if ( .not. associated(track_z1(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_z1) ) then
      allocate(track_z1(initSize))
    else if ( i > size(track_z1) ) then
      temp => track_z1
      allocate ( track_z1(2*size(temp)) )
      track_z1(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_z1(i)%p => what(:)
    track_z1(i)%where = trim(where) // "@" // trimModule(module)
    num_z1 = max(num_z1,i)
  end subroutine TrackAllocate_Z1
  subroutine TrackAllocate_Z2 ( What, Where, Module )
    complex(kind(0.0d0)), pointer :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_z2_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_z2
      if ( .not. associated(track_z2(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_z2) ) then
      allocate(track_z2(initSize))
    else if ( i > size(track_z2) ) then
      temp => track_z2
      allocate ( track_z2(2*size(temp)) )
      track_z2(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_z2(i)%p => what(:,:)
    track_z2(i)%where = trim(where) // "@" // trimModule(module)
    num_z2 = max(num_z2,i)
  end subroutine TrackAllocate_Z2
  subroutine TrackAllocate_Z3 ( What, Where, Module )
    complex(kind(0.0d0)), pointer :: What(:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_z3_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_z3
      if ( .not. associated(track_z3(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_z3) ) then
      allocate(track_z3(initSize))
    else if ( i > size(track_z3) ) then
      temp => track_z3
      allocate ( track_z3(2*size(temp)) )
      track_z3(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_z3(i)%p => what(:,:,:)
    track_z3(i)%where = trim(where) // "@" // trimModule(module)
    num_z3 = max(num_z3,i)
  end subroutine TrackAllocate_Z3

  subroutine TrackAllocate_I1 ( What, Where, Module )
    integer, pointer :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_i1_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_i1
      if ( .not. associated(track_i1(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_i1) ) then
      allocate(track_i1(initSize))
    else if ( i > size(track_i1) ) then
      temp => track_i1
      allocate ( track_i1(2*size(temp)) )
      track_i1(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_i1(i)%p => what
    track_i1(i)%where = trim(where) // "@" // trimModule(module)
    num_i1 = max(num_i1,i)
  end subroutine TrackAllocate_I1
  subroutine TrackAllocate_I2 ( What, Where, Module )
    integer, pointer :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_i2_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_i2
      if ( .not. associated(track_i2(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_i2) ) then
      allocate(track_i2(initSize))
    else if ( i > size(track_i2) ) then
      temp => track_i2
      allocate ( track_i2(2*size(temp)) )
      track_i2(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_i2(i)%p => what
    track_i2(i)%where = trim(where) // "@" // trimModule(module)
    num_i2 = max(num_i2,i)
  end subroutine TrackAllocate_I2
  subroutine TrackAllocate_I3 ( What, Where, Module )
    integer, pointer :: What(:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_i3_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_i3
      if ( .not. associated(track_i3(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_i3) ) then
      allocate(track_i3(initSize))
    else if ( i > size(track_i3) ) then
      temp => track_i3
      allocate ( track_i3(2*size(temp)) )
      track_i3(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_i3(i)%p => what
    track_i3(i)%where = trim(where) // "@" // trimModule(module)
    num_i3 = max(num_i3,i)
  end subroutine TrackAllocate_I3
  subroutine TrackAllocate_I4 ( What, Where, Module )
    integer, pointer :: What(:,:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_i4_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_i4
      if ( .not. associated(track_i4(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_i4) ) then
      allocate(track_i4(initSize))
    else if ( i > size(track_i4) ) then
      temp => track_i4
      allocate ( track_i4(2*size(temp)) )
      track_i4(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_i4(i)%p => what
    track_i4(i)%where = trim(where) // "@" // trimModule(module)
    num_i4 = max(num_i4,i)
  end subroutine TrackAllocate_I4

  subroutine TrackAllocate_R1 ( What, Where, Module )
    real, pointer :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_r1_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r1
      if ( .not. associated(track_r1(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_r1) ) then
      allocate(track_r1(initSize))
    else if ( i > size(track_r1) ) then
      temp => track_r1
      allocate ( track_r1(2*size(temp)) )
      track_r1(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_r1(i)%p => what
    track_r1(i)%where = trim(where) // "@" // trimModule(module)
    num_r1 = max(num_r1,i)
  end subroutine TrackAllocate_R1
  subroutine TrackAllocate_R2 ( What, Where, Module )
    real, pointer :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_r2_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r2
      if ( .not. associated(track_r2(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_r2) ) then
      allocate(track_r2(initSize))
    else if ( i > size(track_r2) ) then
      temp => track_r2
      allocate ( track_r2(2*size(temp)) )
      track_r2(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_r2(i)%p => what
    track_r2(i)%where = trim(where) // "@" // trimModule(module)
    num_r2 = max(num_r2,i)
  end subroutine TrackAllocate_R2
  subroutine TrackAllocate_R3 ( What, Where, Module )
    real, pointer :: What(:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_r3_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r3
      if ( .not. associated(track_r3(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_r3) ) then
      allocate(track_r3(initSize))
    else if ( i > size(track_r3) ) then
      temp => track_r3
      allocate ( track_r3(2*size(temp)) )
      track_r3(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_r3(i)%p => what
    track_r3(i)%where = trim(where) // "@" // trimModule(module)
    num_r3 = max(num_r3,i)
  end subroutine TrackAllocate_R3
  subroutine TrackAllocate_R4 ( What, Where, Module )
    real, pointer :: What(:,:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_r4_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r4
      if ( .not. associated(track_r4(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_r4) ) then
      allocate(track_r4(initSize))
    else if ( i > size(track_r4) ) then
      temp => track_r4
      allocate ( track_r4(2*size(temp)) )
      track_r4(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_r4(i)%p => what
    track_r4(i)%where = trim(where) // "@" // trimModule(module)
    num_r4 = max(num_r4,i)
  end subroutine TrackAllocate_R4

  subroutine TrackAllocateA_C1 ( What, Where, Module )
    character(len=*), allocatable :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_c1_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_c1
      if ( .not. associated(track_c1(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_c1) ) then
      allocate(track_c1(initSize))
    else if ( i > size(track_c1) ) then
      temp => track_c1
      allocate ( track_c1(2*size(temp)) )
      track_c1(:size(temp)) = temp
      deallocate ( temp )
    end if
9   nullify ( track_c1(i)%p )
    track_c1(i)%where = trim(where) // "@" // trimModule(module)
    num_c1 = max(num_c1,i)
  end subroutine TrackAllocateA_C1
  subroutine TrackAllocateA_C2 ( What, Where, Module )
    character(len=*), allocatable :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_c2_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_c2
      if ( .not. associated(track_c2(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_c2) ) then
      allocate(track_c2(initSize))
    else if ( i > size(track_c2) ) then
      temp => track_c2
      allocate ( track_c2(2*size(temp)) )
      track_c2(:size(temp)) = temp
      deallocate ( temp )
    end if
9   nullify ( track_c2(i)%p )
    track_c2(i)%where = trim(where) // "@" // trimModule(module)
    num_c2 = max(num_c2,i)
  end subroutine TrackAllocateA_C2
  subroutine TrackAllocateA_I1 ( What, Where, Module )
    integer, allocatable :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_i1_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r1
      if ( .not. associated(track_i1(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_i1) ) then
      allocate(track_i1(initSize))
    else if ( i > size(track_i1) ) then
      temp => track_i1
      allocate ( track_i1(2*size(temp)) )
      track_i1(:size(temp)) = temp
      deallocate ( temp )
    end if
9   nullify ( track_i1(i)%p )
    track_i1(i)%where = trim(where) // "@" // trimModule(module)
    num_r1 = max(num_r1,i)
  end subroutine TrackAllocateA_I1
  subroutine TrackAllocateA_I2 ( What, Where, Module )
    integer, allocatable :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_i2_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r2
      if ( .not. associated(track_i2(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_i2) ) then
      allocate(track_i2(initSize))
    else if ( i > size(track_i2) ) then
      temp => track_i2
      allocate ( track_i2(2*size(temp)) )
      track_i2(:size(temp)) = temp
      deallocate ( temp )
    end if
9   nullify ( track_i2(i)%p )
    track_i2(i)%where = trim(where) // "@" // trimModule(module)
    num_r2 = max(num_r2,i)
  end subroutine TrackAllocateA_I2
  subroutine TrackAllocateA_I3 ( What, Where, Module )
    integer, allocatable :: What(:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_i3_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r3
      if ( .not. associated(track_i3(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_i3) ) then
      allocate(track_i3(initSize))
    else if ( i > size(track_i3) ) then
      temp => track_i3
      allocate ( track_i3(2*size(temp)) )
      track_i3(:size(temp)) = temp
      deallocate ( temp )
    end if
9   nullify ( track_i3(i)%p )
    track_i3(i)%where = trim(where) // "@" // trimModule(module)
    num_r3 = max(num_r3,i)
  end subroutine TrackAllocateA_I3
  subroutine TrackAllocateA_I4 ( What, Where, Module )
    integer, allocatable :: What(:,:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_i4_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r4
      if ( .not. associated(track_i4(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_i4) ) then
      allocate(track_i4(initSize))
    else if ( i > size(track_i4) ) then
      temp => track_i4
      allocate ( track_i4(2*size(temp)) )
      track_i4(:size(temp)) = temp
      deallocate ( temp )
    end if
9   nullify ( track_i4(i)%p )
    track_i4(i)%where = trim(where) // "@" // trimModule(module)
    num_r4 = max(num_r4,i)
  end subroutine TrackAllocateA_I4

  subroutine TrackAllocateA_L1 ( What, Where, Module )
    logical, allocatable :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_L1_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r1
      if ( .not. associated(track_L1(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_L1) ) then
      allocate(track_L1(initSize))
    else if ( i > size(track_L1) ) then
      temp => track_L1
      allocate ( track_L1(2*size(temp)) )
      track_L1(:size(temp)) = temp
      deallocate ( temp )
    end if
9   nullify ( track_L1(i)%p )
    track_L1(i)%where = trim(where) // "@" // trimModule(module)
    num_r1 = max(num_r1,i)
  end subroutine TrackAllocateA_L1
  subroutine TrackAllocateA_L2 ( What, Where, Module )
    logical, allocatable :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_L2_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r2
      if ( .not. associated(track_L2(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_L2) ) then
      allocate(track_L2(initSize))
    else if ( i > size(track_L2) ) then
      temp => track_L2
      allocate ( track_L2(2*size(temp)) )
      track_L2(:size(temp)) = temp
      deallocate ( temp )
    end if
9   nullify ( track_L2(i)%p )
    track_L2(i)%where = trim(where) // "@" // trimModule(module)
    num_r2 = max(num_r2,i)
  end subroutine TrackAllocateA_L2
  subroutine TrackAllocateA_L3 ( What, Where, Module )
    logical, allocatable :: What(:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_L3_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r3
      if ( .not. associated(track_L3(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_L3) ) then
      allocate(track_L3(initSize))
    else if ( i > size(track_L3) ) then
      temp => track_L3
      allocate ( track_L3(2*size(temp)) )
      track_L3(:size(temp)) = temp
      deallocate ( temp )
    end if
9   nullify ( track_L3(i)%p )
    track_L3(i)%where = trim(where) // "@" // trimModule(module)
    num_r3 = max(num_r3,i)
  end subroutine TrackAllocateA_L3

  subroutine TrackAllocateA_R1 ( What, Where, Module )
    real, allocatable :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_r1_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r1
      if ( .not. associated(track_r1(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_r1) ) then
      allocate(track_r1(initSize))
    else if ( i > size(track_r1) ) then
      temp => track_r1
      allocate ( track_r1(2*size(temp)) )
      track_r1(:size(temp)) = temp
      deallocate ( temp )
    end if
9   nullify ( track_r1(i)%p )
    track_r1(i)%where = trim(where) // "@" // trimModule(module)
    num_r1 = max(num_r1,i)
  end subroutine TrackAllocateA_R1
  subroutine TrackAllocateA_R2 ( What, Where, Module )
    real, allocatable :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_r2_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r2
      if ( .not. associated(track_r2(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_r2) ) then
      allocate(track_r2(initSize))
    else if ( i > size(track_r2) ) then
      temp => track_r2
      allocate ( track_r2(2*size(temp)) )
      track_r2(:size(temp)) = temp
      deallocate ( temp )
    end if
9   nullify ( track_r2(i)%p )
    track_r2(i)%where = trim(where) // "@" // trimModule(module)
    num_r2 = max(num_r2,i)
  end subroutine TrackAllocateA_R2
  subroutine TrackAllocateA_R3 ( What, Where, Module )
    real, allocatable :: What(:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_r3_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r3
      if ( .not. associated(track_r3(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_r3) ) then
      allocate(track_r3(initSize))
    else if ( i > size(track_r3) ) then
      temp => track_r3
      allocate ( track_r3(2*size(temp)) )
      track_r3(:size(temp)) = temp
      deallocate ( temp )
    end if
9   nullify ( track_r3(i)%p )
    track_r3(i)%where = trim(where) // "@" // trimModule(module)
    num_r3 = max(num_r3,i)
  end subroutine TrackAllocateA_R3
  subroutine TrackAllocateA_R4 ( What, Where, Module )
    real, allocatable :: What(:,:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_r4_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r4
      if ( .not. associated(track_r4(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_r4) ) then
      allocate(track_r4(initSize))
    else if ( i > size(track_r4) ) then
      temp => track_r4
      allocate ( track_r4(2*size(temp)) )
      track_r4(:size(temp)) = temp
      deallocate ( temp )
    end if
9   nullify ( track_r4(i)%p )
    track_r4(i)%where = trim(where) // "@" // trimModule(module)
    num_r4 = max(num_r4,i)
  end subroutine TrackAllocateA_R4

  subroutine TrackAllocate_D1 ( What, Where, Module )
    double precision, pointer :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_d1_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_d1
      if ( .not. associated(track_d1(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_d1) ) then
      allocate(track_d1(initSize))
    else if ( i > size(track_d1) ) then
      temp => track_d1
      allocate ( track_d1(2*size(temp)) )
      track_d1(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_d1(i)%p => what
    track_d1(i)%where = trim(where) // "@" // trimModule(module)
    num_d1 = max(num_d1,i)
  end subroutine TrackAllocate_D1
  subroutine TrackAllocate_D2 ( What, Where, Module )
    double precision, pointer :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_d2_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_d2
      if ( .not. associated(track_d2(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_d2) ) then
      allocate(track_d2(initSize))
    else if ( i > size(track_d2) ) then
      temp => track_d2
      allocate ( track_d2(2*size(temp)) )
      track_d2(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_d2(i)%p => what
    track_d2(i)%where = trim(where) // "@" // trimModule(module)
    num_d2 = max(num_d2,i)
  end subroutine TrackAllocate_D2
  subroutine TrackAllocate_D3 ( What, Where, Module )
    double precision, pointer :: What(:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_d3_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_d3
      if ( .not. associated(track_d3(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_d3) ) then
      allocate(track_d3(initSize))
    else if ( i > size(track_d3) ) then
      temp => track_d3
      allocate ( track_d3(2*size(temp)) )
      track_d3(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_d3(i)%p => what
    track_d3(i)%where = trim(where) // "@" // trimModule(module)
    num_d3 = max(num_d3,i)
  end subroutine TrackAllocate_D3
  subroutine TrackAllocate_D4 ( What, Where, Module )
    double precision, pointer :: What(:,:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_d4_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_d4
      if ( .not. associated(track_d4(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_d4) ) then
      allocate(track_d4(initSize))
    else if ( i > size(track_d4) ) then
      temp => track_d4
      allocate ( track_d4(2*size(temp)) )
      track_d4(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_d4(i)%p => what
    track_d4(i)%where = trim(where) // "@" // trimModule(module)
    num_d4 = max(num_d4,i)
  end subroutine TrackAllocate_D4

  subroutine TrackAllocateA_D1 ( What, Where, Module )
    double precision, allocatable :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_d1_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_d1
      if ( .not. associated(track_d1(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_d1) ) then
      allocate(track_d1(initSize))
    else if ( i > size(track_d1) ) then
      temp => track_d1
      allocate ( track_d1(2*size(temp)) )
      track_d1(:size(temp)) = temp
      deallocate ( temp )
    end if
9   nullify ( track_d1(i)%p )
    track_d1(i)%where = trim(where) // "@" // trimModule(module)
    num_d1 = max(num_d1,i)
  end subroutine TrackAllocateA_D1
  subroutine TrackAllocateA_D2 ( What, Where, Module )
    double precision, allocatable :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_d2_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_d2
      if ( .not. associated(track_d2(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_d2) ) then
      allocate(track_d2(initSize))
    else if ( i > size(track_d2) ) then
      temp => track_d2
      allocate ( track_d2(2*size(temp)) )
      track_d2(:size(temp)) = temp
      deallocate ( temp )
    end if
9   nullify ( track_d2(i)%p )
    track_d2(i)%where = trim(where) // "@" // trimModule(module)
    num_d2 = max(num_d2,i)
  end subroutine TrackAllocateA_D2
  subroutine TrackAllocateA_D3 ( What, Where, Module )
    double precision, allocatable :: What(:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_d3_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_d3
      if ( .not. associated(track_d3(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_d3) ) then
      allocate(track_d3(initSize))
    else if ( i > size(track_d3) ) then
      temp => track_d3
      allocate ( track_d3(2*size(temp)) )
      track_d3(:size(temp)) = temp
      deallocate ( temp )
    end if
9   nullify ( track_d3(i)%p )
    track_d3(i)%where = trim(where) // "@" // trimModule(module)
    num_d3 = max(num_d3,i)
  end subroutine TrackAllocateA_D3
  subroutine TrackAllocateA_D4 ( What, Where, Module )
    double precision, allocatable :: What(:,:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_d4_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_d4
      if ( .not. associated(track_d4(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_d4) ) then
      allocate(track_d4(initSize))
    else if ( i > size(track_d4) ) then
      temp => track_d4
      allocate ( track_d4(2*size(temp)) )
      track_d4(:size(temp)) = temp
      deallocate ( temp )
    end if
9   nullify ( track_d4(i)%p )
    track_d4(i)%where = trim(where) // "@" // trimModule(module)
    num_d4 = max(num_d4,i)
  end subroutine TrackAllocateA_D4

  subroutine TrackAllocate_L1 ( What, Where, Module )
    logical, pointer :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_l1_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_l1
      if ( .not. associated(track_l1(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_l1) ) then
      allocate(track_l1(initSize))
    else if ( i > size(track_l1) ) then
      temp => track_l1
      allocate ( track_l1(2*size(temp)) )
      track_l1(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_l1(i)%p => what(:)
    track_l1(i)%where = trim(where) // "@" // trimModule(module)
    num_l1 = max(num_l1,i)
  end subroutine TrackAllocate_L1
  subroutine TrackAllocate_L2 ( What, Where, Module )
    logical, pointer :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_l2_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_l2
      if ( .not. associated(track_l2(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_l2) ) then
      allocate(track_l2(initSize))
    else if ( i > size(track_l2) ) then
      temp => track_l2
      allocate ( track_l2(2*size(temp)) )
      track_l2(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_l2(i)%p => what(:,:)
    track_l2(i)%where = trim(where) // "@" // trimModule(module)
    num_l2 = max(num_l2,i)
  end subroutine TrackAllocate_L2
  subroutine TrackAllocate_L3 ( What, Where, Module )
    logical, pointer :: What(:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    type(track_l3_t), pointer :: temp(:)
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_l3
      if ( .not. associated(track_l3(i)%p) ) go to 9 ! Found empty one
    end do
    if ( .not. associated(track_l3) ) then
      allocate(track_l3(initSize))
    else if ( i > size(track_l3) ) then
      temp => track_l3
      allocate ( track_l3(2*size(temp)) )
      track_l3(:size(temp)) = temp
      deallocate ( temp )
    end if
9   track_l3(i)%p => what(:,:,:)
    track_l3(i)%where = trim(where) // "@" // trimModule(module)
    num_l3 = max(num_l3,i)
  end subroutine TrackAllocate_L3

  subroutine TrackDeallocate_C1 ( What, Where, Module )
    character(len=*), pointer :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_c1
      if ( associated(track_c1(i)%p, what(:)(1:1)) ) then
        nullify ( track_c1(i)%p )  ! Mark it free
        do num_c1 = num_c1, 1, -1  ! Reduce count to speed searches
          if ( associated(track_c1(num_c1)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_C1
  subroutine TrackDeallocate_C2 ( What, Where, Module )
    character(len=*), pointer :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_c2
      if ( associated(track_c2(i)%p, what(:,:)(1:1)) ) then
        nullify ( track_c2(i)%p )  ! Mark it free
        do num_c2 = num_c2, 1, -1  ! Reduce count to speed searches
          if ( associated(track_c2(num_c2)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_C2
  subroutine TrackDeallocate_C3 ( What, Where, Module )
    character(len=*), pointer :: What(:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_c3
      if ( associated(track_c3(i)%p, what(:,:,:)(1:1)) ) then
        nullify ( track_c3(i)%p )  ! Mark it free
        do num_c3 = num_c3, 1, -1  ! Reduce count to speed searches
          if ( associated(track_c3(num_c3)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_C3

  subroutine TrackDeallocate_X1 ( What, Where, Module )
    complex, pointer :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_x1
      if ( associated(track_x1(i)%p, what) ) then
        nullify ( track_x1(i)%p )  ! Mark it free
        do num_x1 = num_x1, 1, -1  ! Reduce count to speed searches
          if ( associated(track_x1(num_x1)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_X1
  subroutine TrackDeallocate_X2 ( What, Where, Module )
    complex, pointer :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_x2
      if ( associated(track_x2(i)%p, what) ) then
        nullify ( track_x2(i)%p )  ! Mark it free
        do num_x2 = num_x2, 1, -1  ! Reduce count to speed searches
          if ( associated(track_x2(num_x2)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_X2
  subroutine TrackDeallocate_X3 ( What, Where, Module )
    complex, pointer :: What(:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_x3
      if ( associated(track_x3(i)%p, what) ) then
        nullify ( track_x3(i)%p )  ! Mark it free
        do num_x3 = num_x3, 1, -1  ! Reduce count to speed searches
          if ( associated(track_x3(num_x3)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_X3

  subroutine TrackDeallocate_Z1 ( What, Where, Module )
    complex(kind(0.0d0)), pointer :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_z1
      if ( associated(track_z1(i)%p, what) ) then
        nullify ( track_z1(i)%p )  ! Mark it free
        do num_z1 = num_z1, 1, -1  ! Reduce count to speed searches
          if ( associated(track_z1(num_z1)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_Z1
  subroutine TrackDeallocate_Z2 ( What, Where, Module )
    complex(kind(0.0d0)), pointer :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_z2
      if ( associated(track_z2(i)%p, what) ) then
        nullify ( track_z2(i)%p )  ! Mark it free
        do num_z2 = num_z2, 1, -1  ! Reduce count to speed searches
          if ( associated(track_z2(num_z2)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_Z2
  subroutine TrackDeallocate_Z3 ( What, Where, Module )
    complex(kind(0.0d0)), pointer :: What(:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_z3
      if ( associated(track_z3(i)%p, what) ) then
        nullify ( track_z3(i)%p )  ! Mark it free
        do num_z3 = num_z3, 1, -1  ! Reduce count to speed searches
          if ( associated(track_z3(num_z3)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_Z3

  subroutine TrackDeallocate_I1 ( What, Where, Module )
    integer, pointer :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_i1
      if ( associated(track_i1(i)%p, what) ) then
        nullify ( track_i1(i)%p )  ! Mark it free
        do num_i1 = num_i1, 1, -1  ! Reduce count to speed searches
          if ( associated(track_i1(num_i1)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_I1
  subroutine TrackDeallocate_I2 ( What, Where, Module )
    integer, pointer :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_i2
      if ( associated(track_i2(i)%p, what) ) then
        nullify ( track_i2(i)%p )  ! Mark it free
        do num_i2 = num_i2, 1, -1  ! Reduce count to speed searches
          if ( associated(track_i2(num_i2)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_I2
  subroutine TrackDeallocate_I3 ( What, Where, Module )
    integer, pointer :: What(:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_i3
      if ( associated(track_i3(i)%p, what) ) then
        nullify ( track_i3(i)%p )  ! Mark it free
        do num_i3 = num_i3, 1, -1  ! Reduce count to speed searches
          if ( associated(track_i3(num_i3)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_I3
  subroutine TrackDeallocate_I4 ( What, Where, Module )
    integer, pointer :: What(:,:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_i4
      if ( associated(track_i4(i)%p, what) ) then
        nullify ( track_i4(i)%p )  ! Mark it free
        do num_i4 = num_i4, 1, -1  ! Reduce count to speed searches
          if ( associated(track_i4(num_i4)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_I4

  subroutine TrackDeallocate_R1 ( What, Where, Module )
    real, pointer :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r1
      if ( associated(track_r1(i)%p, what) ) then
        nullify ( track_r1(i)%p )  ! Mark it free
        do num_r1 = num_r1, 1, -1  ! Reduce count to speed searches
          if ( associated(track_r1(num_r1)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_R1
  subroutine TrackDeallocate_R2 ( What, Where, Module )
    real, pointer :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r2
      if ( associated(track_r2(i)%p, what) ) then
        nullify ( track_r2(i)%p )  ! Mark it free
        do num_r2 = num_r2, 1, -1  ! Reduce count to speed searches
          if ( associated(track_r2(num_r2)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_R2
  subroutine TrackDeallocate_R3 ( What, Where, Module )
    real, pointer :: What(:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r3
      if ( associated(track_r3(i)%p, what) ) then
        nullify ( track_r3(i)%p )  ! Mark it free
        do num_r3 = num_r3, 1, -1  ! Reduce count to speed searches
          if ( associated(track_r3(num_r3)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_R3
  subroutine TrackDeallocate_R4 ( What, Where, Module )
    real, pointer :: What(:,:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r4
      if ( associated(track_r4(i)%p, what) ) then
        nullify ( track_r4(i)%p )  ! Mark it free
        do num_r4 = num_r4, 1, -1  ! Reduce count to speed searches
          if ( associated(track_r4(num_r4)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_R4

  subroutine TrackDeallocateA_C1 ( What, Where, Module )
    character(len=*), allocatable, target :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_c1
      if ( associated(track_c1(i)%p, what) ) then
        nullify ( track_c1(i)%p )  ! Mark it free
        do num_c1 = num_c1, 1, -1  ! Reduce count to speed searches
          if ( associated(track_c1(num_c1)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocateA_C1
  subroutine TrackDeallocateA_C2 ( What, Where, Module )
    character(len=*), allocatable, target :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_c2
      if ( associated(track_c2(i)%p, what) ) then
        nullify ( track_c2(i)%p )  ! Mark it free
        do num_c2 = num_c2, 1, -1  ! Reduce count to speed searches
          if ( associated(track_c2(num_c2)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocateA_C2
  subroutine TrackDeallocateA_I1 ( What, Where, Module )
    integer, allocatable, target :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r1
      if ( associated(track_I1(i)%p, what) ) then
        nullify ( track_I1(i)%p )  ! Mark it free
        do num_r1 = num_r1, 1, -1  ! Reduce count to speed searches
          if ( associated(track_I1(num_r1)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocateA_I1
  subroutine TrackDeallocateA_I2 ( What, Where, Module )
    integer, allocatable, target :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r2
      if ( associated(track_I2(i)%p, what) ) then
        nullify ( track_I2(i)%p )  ! Mark it free
        do num_r2 = num_r2, 1, -1  ! Reduce count to speed searches
          if ( associated(track_I2(num_r2)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocateA_I2
  subroutine TrackDeallocateA_I3 ( What, Where, Module )
    integer, allocatable, target :: What(:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r3
      if ( associated(track_I3(i)%p, what) ) then
        nullify ( track_I3(i)%p )  ! Mark it free
        do num_r3 = num_r3, 1, -1  ! Reduce count to speed searches
          if ( associated(track_I3(num_r3)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocateA_I3
  subroutine TrackDeallocateA_I4 ( What, Where, Module )
    integer, allocatable, target :: What(:,:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r4
      if ( associated(track_I4(i)%p, what) ) then
        nullify ( track_I4(i)%p )  ! Mark it free
        do num_r4 = num_r4, 1, -1  ! Reduce count to speed searches
          if ( associated(track_I4(num_r4)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocateA_I4

  subroutine TrackDeallocateA_L1 ( What, Where, Module )
    logical, allocatable, target :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r1
      if ( associated(track_L1(i)%p, what) ) then
        nullify ( track_L1(i)%p )  ! Mark it free
        do num_r1 = num_r1, 1, -1  ! Reduce count to speed searches
          if ( associated(track_L1(num_r1)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocateA_L1
  subroutine TrackDeallocateA_L2 ( What, Where, Module )
    logical, allocatable, target :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r2
      if ( associated(track_L2(i)%p, what) ) then
        nullify ( track_L2(i)%p )  ! Mark it free
        do num_r2 = num_r2, 1, -1  ! Reduce count to speed searches
          if ( associated(track_L2(num_r2)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocateA_L2
  subroutine TrackDeallocateA_L3 ( What, Where, Module )
    logical, allocatable, target :: What(:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r3
      if ( associated(track_L3(i)%p, what) ) then
        nullify ( track_L3(i)%p )  ! Mark it free
        do num_r3 = num_r3, 1, -1  ! Reduce count to speed searches
          if ( associated(track_L3(num_r3)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocateA_L3

  subroutine TrackDeallocateA_R1 ( What, Where, Module )
    real, allocatable, target :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r1
      if ( associated(track_r1(i)%p, what) ) then
        nullify ( track_r1(i)%p )  ! Mark it free
        do num_r1 = num_r1, 1, -1  ! Reduce count to speed searches
          if ( associated(track_r1(num_r1)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocateA_R1
  subroutine TrackDeallocateA_R2 ( What, Where, Module )
    real, allocatable, target :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r2
      if ( associated(track_r2(i)%p, what) ) then
        nullify ( track_r2(i)%p )  ! Mark it free
        do num_r2 = num_r2, 1, -1  ! Reduce count to speed searches
          if ( associated(track_r2(num_r2)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocateA_R2
  subroutine TrackDeallocateA_R3 ( What, Where, Module )
    real, allocatable, target :: What(:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r3
      if ( associated(track_r3(i)%p, what) ) then
        nullify ( track_r3(i)%p )  ! Mark it free
        do num_r3 = num_r3, 1, -1  ! Reduce count to speed searches
          if ( associated(track_r3(num_r3)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocateA_R3
  subroutine TrackDeallocateA_R4 ( What, Where, Module )
    real, allocatable, target :: What(:,:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_r4
      if ( associated(track_r4(i)%p, what) ) then
        nullify ( track_r4(i)%p )  ! Mark it free
        do num_r4 = num_r4, 1, -1  ! Reduce count to speed searches
          if ( associated(track_r4(num_r4)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocateA_R4

  subroutine TrackDeallocate_D1 ( What, Where, Module )
    double precision, pointer :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_d1
      if ( associated(track_d1(i)%p, what) ) then
        nullify ( track_d1(i)%p )  ! Mark it free
        do num_d1 = num_d1, 1, -1  ! Reduce count to speed searches
          if ( associated(track_d1(num_d1)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_D1
  subroutine TrackDeallocate_D2 ( What, Where, Module )
    double precision, pointer :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_d2
      if ( associated(track_d2(i)%p, what) ) then
        nullify ( track_d2(i)%p )  ! Mark it free
        do num_d2 = num_d2, 1, -1  ! Reduce count to speed searches
          if ( associated(track_d2(num_d2)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_D2
  subroutine TrackDeallocate_D3 ( What, Where, Module )
    double precision, pointer :: What(:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_d3
      if ( associated(track_d3(i)%p, what) ) then
        nullify ( track_d3(i)%p )  ! Mark it free
        do num_d3 = num_d3, 1, -1  ! Reduce count to speed searches
          if ( associated(track_d3(num_d3)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_D3
  subroutine TrackDeallocate_D4 ( What, Where, Module )
    double precision, pointer :: What(:,:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_d4
      if ( associated(track_d4(i)%p, what) ) then
        nullify ( track_d4(i)%p )  ! Mark it free
        do num_d4 = num_d4, 1, -1  ! Reduce count to speed searches
          if ( associated(track_d4(num_d4)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_D4

  subroutine TrackDeallocateA_D1 ( What, Where, Module )
    double precision, allocatable, target :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_d1
      if ( associated(track_d1(i)%p, what) ) then
        nullify ( track_d1(i)%p )  ! Mark it free
        do num_d1 = num_d1, 1, -1  ! Reduce count to speed searches
          if ( associated(track_d1(num_d1)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocateA_D1
  subroutine TrackDeallocateA_D2 ( What, Where, Module )
    double precision, allocatable, target :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_d2
      if ( associated(track_d2(i)%p, what) ) then
        nullify ( track_d2(i)%p )  ! Mark it free
        do num_d2 = num_d2, 1, -1  ! Reduce count to speed searches
          if ( associated(track_d2(num_d2)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocateA_D2
  subroutine TrackDeallocateA_D3 ( What, Where, Module )
    double precision, allocatable, target :: What(:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_d3
      if ( associated(track_d3(i)%p, what) ) then
        nullify ( track_d3(i)%p )  ! Mark it free
        do num_d3 = num_d3, 1, -1  ! Reduce count to speed searches
          if ( associated(track_d3(num_d3)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocateA_D3
  subroutine TrackDeallocateA_D4 ( What, Where, Module )
    double precision, allocatable, target :: What(:,:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_d4
      if ( associated(track_d4(i)%p, what) ) then
        nullify ( track_d4(i)%p )  ! Mark it free
        do num_d4 = num_d4, 1, -1  ! Reduce count to speed searches
          if ( associated(track_d4(num_d4)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocateA_D4

  subroutine TrackDeallocate_L1 ( What, Where, Module )
    logical, pointer :: What(:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_l1
      if ( associated(track_l1(i)%p, what) ) then
        nullify ( track_l1(i)%p )  ! Mark it free
        do num_l1 = num_l1, 1, -1  ! Reduce count to speed searches
          if ( associated(track_l1(num_l1)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_L1
  subroutine TrackDeallocate_L2 ( What, Where, Module )
    logical, pointer :: What(:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_l2
      if ( associated(track_l2(i)%p, what) ) then
        nullify ( track_l2(i)%p )  ! Mark it free
        do num_l2 = num_l2, 1, -1  ! Reduce count to speed searches
          if ( associated(track_l2(num_l2)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_L2
  subroutine TrackDeallocate_L3 ( What, Where, Module )
    logical, pointer :: What(:,:,:)
    character(len=*), intent(in) :: Where, Module
    integer :: I
    if ( size(what) == 0 ) return ! Never associated with another, so we can't track them
    do i = 1, num_l3
      if ( associated(track_l3(i)%p, what) ) then
        nullify ( track_l3(i)%p )  ! Mark it free
        do num_l3 = num_l3, 1, -1  ! Reduce count to speed searches
          if ( associated(track_l3(num_l3)%p) ) return
        end do
        return
      end if
    end do
    write ( *, '("No allocation found for ", a, " in ", a)' ) trim(where), trim(trimModule(module))
  end subroutine TrackDeallocate_L3

  subroutine ReportLeaks ( Where )
    use HIGHOUTPUT, only: DUMPSIZE
    character(len=*), intent(in), optional :: Where
    double precision :: Total ! memory still allocated that we could find
    if ( present(where) )  write ( *, '(a)' ) trim(where)
    write ( *, '(a)' ) "Checking for leaks involving arrays of one, two, or three dimensions", &
                     & "of type character, complex, double precision complex, integer, real", &
                     & "double precision, and logical.  Other types and ranks not tracked."
    total = 0
    call reportLeaks_c1 ( total )
    call reportLeaks_c2 ( total )
    call reportLeaks_c3 ( total )
    call reportLeaks_x1 ( total )
    call reportLeaks_x2 ( total )
    call reportLeaks_x3 ( total )
    call reportLeaks_z1 ( total )
    call reportLeaks_z2 ( total )
    call reportLeaks_z3 ( total )
    call reportLeaks_i1 ( total )
    call reportLeaks_i2 ( total )
    call reportLeaks_i3 ( total )
    call reportLeaks_i4 ( total )
    call reportLeaks_r1 ( total )
    call reportLeaks_r2 ( total )
    call reportLeaks_r3 ( total )
    call reportLeaks_r4 ( total )
    call reportLeaks_d1 ( total )
    call reportLeaks_d2 ( total )
    call reportLeaks_d3 ( total )
    call reportLeaks_d4 ( total )
    call reportLeaks_l1 ( total )
    call reportLeaks_l2 ( total )
    call reportLeaks_l3 ( total )
    call dumpsize ( total, advance='yes', before='At least ', &
      & after=' total still allocated.' )
  end subroutine ReportLeaks

  subroutine ReportLeaks_c1 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_c1
      if ( associated(track_c1(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
         &  track_c1(i)%where, size(track_c1(i)%p)
        total = total + len(track_c1(i)%p) * size(track_c1(i)%p)
      end if
    end do
  end subroutine ReportLeaks_c1
  subroutine ReportLeaks_c2 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_c2
      if ( associated(track_c2(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
          & track_c2(i)%where, size(track_c2(i)%p)
        total = total + len(track_c2(i)%p) * size(track_c2(i)%p)
      end if
    end do
  end subroutine ReportLeaks_c2
  subroutine ReportLeaks_c3 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_c3
      if ( associated(track_c3(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
          & track_c3(i)%where, size(track_c3(i)%p)
        total = total + len(track_c3(i)%p) * size(track_c3(i)%p)
      end if
    end do
  end subroutine ReportLeaks_c3

  subroutine ReportLeaks_x1 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_x1
      if ( associated(track_x1(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
          & track_x1(i)%where, size(track_x1(i)%p)
        total = total + 2 * charsPerInt * size(track_x1(i)%p)
      end if
    end do
  end subroutine ReportLeaks_x1
  subroutine ReportLeaks_x2 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_x2
      if ( associated(track_x2(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
          & track_x2(i)%where, size(track_x2(i)%p)
        total = total + 2 * charsPerInt * size(track_x2(i)%p)
      end if
    end do
  end subroutine ReportLeaks_x2
  subroutine ReportLeaks_x3 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_x3
      if ( associated(track_x3(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
          & track_x3(i)%where, size(track_x3(i)%p)
        total = total + 2 * charsPerInt * size(track_x3(i)%p)
      end if
    end do
  end subroutine ReportLeaks_x3

  subroutine ReportLeaks_z1 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_z1
      if ( associated(track_z1(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
          & track_z1(i)%where, size(track_z1(i)%p)
        total = total + 4 * charsPerInt * size(track_z1(i)%p)
      end if
    end do
  end subroutine ReportLeaks_z1
  subroutine ReportLeaks_z2 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_z2
      if ( associated(track_z2(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
          & track_z2(i)%where, size(track_z2(i)%p)
        total = total + 4 * charsPerInt * size(track_z2(i)%p)
      end if
    end do
  end subroutine ReportLeaks_z2
  subroutine ReportLeaks_z3 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_z3
      if ( associated(track_z3(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
          & track_z3(i)%where, size(track_z3(i)%p)
        total = total + 4 * charsPerInt * size(track_z3(i)%p)
      end if
    end do
  end subroutine ReportLeaks_z3

  subroutine ReportLeaks_i1 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_i1
      if ( associated(track_i1(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
          & track_i1(i)%where, size(track_i1(i)%p)
        total = total + charsPerInt * size(track_i1(i)%p)
      end if
    end do
  end subroutine ReportLeaks_i1
  subroutine ReportLeaks_i2 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_i2
      if ( associated(track_i2(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
          & track_i2(i)%where, size(track_i2(i)%p)
        total = total + charsPerInt * size(track_i2(i)%p)
      end if
    end do
  end subroutine ReportLeaks_i2
  subroutine ReportLeaks_i3 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_i3
      if ( associated(track_i3(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
          & track_i3(i)%where, size(track_i3(i)%p)
        total = total + charsPerInt * size(track_i3(i)%p)
      end if
    end do
  end subroutine ReportLeaks_i3
  subroutine ReportLeaks_i4 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_i4
      if ( associated(track_i4(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
          & track_i4(i)%where, size(track_i4(i)%p)
        total = total + charsPerInt * size(track_i4(i)%p)
      end if
    end do
  end subroutine ReportLeaks_i4

  subroutine ReportLeaks_r1 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_r1
      if ( associated(track_r1(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
          & track_r1(i)%where, size(track_r1(i)%p)
        total = total + charsPerInt * size(track_r1(i)%p)
      end if
    end do
  end subroutine ReportLeaks_r1
  subroutine ReportLeaks_r2 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_r2
      if ( associated(track_r2(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
          & track_r2(i)%where, size(track_r2(i)%p)
        total = total + charsPerInt * size(track_r2(i)%p)
      end if
    end do
  end subroutine ReportLeaks_r2
  subroutine ReportLeaks_r3 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_r3
      if ( associated(track_r3(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
          & track_r3(i)%where, size(track_r3(i)%p)
        total = total + charsPerInt * size(track_r3(i)%p)
      end if
    end do
  end subroutine ReportLeaks_r3
  subroutine ReportLeaks_r4 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_r4
      if ( associated(track_r4(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
          & track_r4(i)%where, size(track_r4(i)%p)
        total = total + charsPerInt * size(track_r4(i)%p)
      end if
    end do
  end subroutine ReportLeaks_r4

  subroutine ReportLeaks_d1 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_d1
      if ( associated(track_d1(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
          & track_d1(i)%where, size(track_d1(i)%p)
        total = total + 2 * charsPerInt * size(track_d1(i)%p)
      end if
    end do
  end subroutine ReportLeaks_d1
  subroutine ReportLeaks_d2 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_d2
      if ( associated(track_d2(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
          & track_d2(i)%where, size(track_d2(i)%p)
        total = total + 2 * charsPerInt * size(track_d2(i)%p)
      end if
    end do
  end subroutine ReportLeaks_d2
  subroutine ReportLeaks_d3 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_d3
      if ( associated(track_d3(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
          & track_d3(i)%where, size(track_d3(i)%p)
        total = total + 2 * charsPerInt * size(track_d3(i)%p)
      end if
    end do
  end subroutine ReportLeaks_d3
  subroutine ReportLeaks_d4 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_d4
      if ( associated(track_d4(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
          & track_d4(i)%where, size(track_d4(i)%p)
        total = total + 2 * charsPerInt * size(track_d4(i)%p)
      end if
    end do
  end subroutine ReportLeaks_d4

  subroutine ReportLeaks_l1 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_l1
      if ( associated(track_l1(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
          & track_l1(i)%where, size(track_l1(i)%p)
        total = total + charsPerInt * size(track_l1(i)%p)
      end if
    end do
  end subroutine ReportLeaks_l1
  subroutine ReportLeaks_l2 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_l2
      if ( associated(track_l2(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
          & track_l2(i)%where, size(track_l2(i)%p)
        total = total + charsPerInt * size(track_l2(i)%p)
      end if
    end do
  end subroutine ReportLeaks_l2
  subroutine ReportLeaks_l3 ( Total )
    double precision, intent(inout) :: Total
    integer :: I
    do i = 1, num_l3
      if ( associated(track_l3(i)%p) ) then
        write ( *, '(a," Still allocated with ", i0, " elements")' ) &
          & track_l3(i)%where, size(track_l3(i)%p)
        total = total + charsPerInt * size(track_l3(i)%p)
      end if
    end do
  end subroutine ReportLeaks_l3

  function TrimModule ( Module )
    character(len=*), intent(in) :: Module
    character(len=len(module)) :: TrimModule
    if ( module(1:1) == '$' ) then
      ! Module is <dollar>RCSFile: <filename>,v <dollar>
      trimModule = module(11:(len_trim(module)-8))
    else
      trimModule = module
    end if
  end function TrimModule

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Track_m

! $Log$
! Revision 2.11  2018/05/15 03:17:36  vsnyder
! Add 1D and 2D character allocatable versions
!
! Revision 2.10  2017/10/31 23:42:57  vsnyder
! Add allocatable logical array support
!
! Revision 2.9  2016/09/07 23:54:37  vsnyder
! Add support for allocatable integer arrays
!
! Revision 2.8  2015/06/02 23:59:42  vsnyder
! Track allocation and deallocation of allocatable arrays
!
! Revision 2.7  2014/01/09 00:24:29  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.6  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.5  2008/05/20 01:59:29  vsnyder
! Add 4d integer, real, double
!
! Revision 2.4  2007/04/14 00:35:34  vsnyder
! Remove declarations for unused variables
!
! Revision 2.3  2006/08/04 18:14:15  vsnyder
! Add size tracking to ReportLeaks, simplify TrackDeallocate
!
! Revision 2.2  2006/07/29 03:42:34  vsnyder
! Can't track zero-size allocations
!
! Revision 2.1  2006/07/29 03:00:54  vsnyder
! Initial commit
!
