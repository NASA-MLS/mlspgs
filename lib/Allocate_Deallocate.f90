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
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error, MLSMSG_Warning

  implicit NONE
  private

  public :: ALLOCATE_TEST, DEALLOCATE_TEST, DEALLOC_STATUS, & 
    & SET_GARBAGE_COLLECTION, REPORTALLOCATEDEALLOCATE

  integer, save :: DEALLOC_STATUS = 0
  logical, save :: COLLECT_GARBAGE_EACH_TIME = .false.

  interface ALLOCATE_TEST
    module procedure ALLOCATE_TEST_CHARACTER_1D
    module procedure ALLOCATE_TEST_CHARACTER_2D
    module procedure ALLOCATE_TEST_CHARACTER_3D
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
    module procedure DEALLOCATE_TEST_CHARACTER_3D
    module procedure DEALLOCATE_TEST_INTEGER_1D, DEALLOCATE_TEST_INTEGER_2D
    module procedure DEALLOCATE_TEST_INTEGER_3D
    module procedure DEALLOCATE_TEST_LOGICAL_1D, DEALLOCATE_TEST_LOGICAL_2D
    module procedure DEALLOCATE_TEST_REALR4_1D, DEALLOCATE_TEST_REALR4_2D
    module procedure DEALLOCATE_TEST_REALR4_3D
    module procedure DEALLOCATE_TEST_REALR8_1D, DEALLOCATE_TEST_REALR8_2D
    module procedure DEALLOCATE_TEST_REALR8_3D
  end interface

  logical, public :: TRACKALLOCATES = .false. ! If true keep track of memory allocated
  logical, public :: CLEARONALLOCATE = .false. ! If true, zero all allocated stuff
  ! and report on it.

  integer, private, save :: NOWORDSALLOCATED=0 ! Number of 4 byte words allocated.

  !------------------------------- RCS Ident Info ------------------------------
  character(len=130), private :: id = & 
       "$Id$"
  character(len=*), parameter, private :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here 
  !-----------------------------------------------------------------------------

contains
  ! =====     Public Procedures      ============================
  subroutine ReportAllocateDeallocate ( name, moduleName, noWords )
    use Output_m, only: OUTPUT
    use Dump_0, only: DUMPSIZE
    ! Dummy arguments
    character (len=*), intent(in) :: NAME ! Name of thing allocated
    character (len=*), intent(in) :: MODULENAME ! Module that allocated it
    integer, intent(in) :: NOWORDS      ! No words allocated (or deallocated if -ve)
    ! Executable code
    if ( .not. trackAllocates ) return        ! Most probably will not be called anyway
    noWordsAllocated = noWordsAllocated + noWords
    call output ( 'Tracking: ' )
    if ( noWords < 0 ) then
      call output ( 'Dea' )
    else
      call output ( 'A' )
    end if
    call output ( 'llocated ' )
    call DumpSize ( abs ( noWords*4.0 ) )
    call output ( ' for ' // trim ( name ) // ' in ' )
    if ( moduleName(1:1) == '$' ) then
      ! The moduleNameIn is <dollar>RCSFile: <filename>,v <dollar>
      call output ( moduleName(11:(len_trim(moduleName)-8)) )
    else
      call output ( moduleName )
    end if
    call output ( ' total ' )
    call DumpSize ( noWordsAllocated*4.0, advance='yes' )
  end subroutine ReportAllocateDeallocate
    
  ! =====     Private Procedures     ============================
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
      & MLSMSG_Allocate // Its_Name // trim(bounds(1,dim1)) )
    if ( status == 0 .and. clearOnAllocate ) to_allocate = ' '
    if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, dim1/4+1 )
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
      & MLSMSG_Allocate // Its_Name // trim(bounds(1,dim1,1,dim2)) )
    if ( status == 0 .and. clearOnAllocate ) to_allocate = ' '
    if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, (dim1*dim2)/4+1 )
  end subroutine Allocate_Test_Character_2d
  ! ---------------------------------  Allocate_Test_Character_3d  -----
  subroutine Allocate_Test_Character_3d ( To_Allocate, Dim1, Dim2, Dim3, &
    & Its_Name, ModuleName )
    character(len=*), pointer, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Second dimension of To_Allocate
    integer, intent(in) :: Dim3    ! Third dimension of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer :: STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    allocate ( To_Allocate(dim1,dim2,dim3), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name // trim(bounds(1,dim1,1,dim2,1,dim3)) )
    if ( status == 0 .and. clearOnAllocate ) to_allocate = ' '
    if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, (dim1*dim2*dim3)/4+1 )
  end subroutine Allocate_Test_Character_3d
  ! ------------------------------------  Allocate_Test_RealR8_1d  -----
  subroutine Allocate_Test_RealR8_1d ( To_Allocate, Dim1, Its_Name, &
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
      & MLSMSG_Allocate // Its_Name // trim(bounds(my_low,dim1)) )
    if ( status == 0 .and. clearOnAllocate ) to_allocate = 0.0d0
    if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, 2*dim1 )
  end subroutine Allocate_Test_RealR8_1d
  ! ------------------------------------  Allocate_Test_RealR8_2d  -----
  subroutine Allocate_Test_RealR8_2d ( To_Allocate, Dim1, Dim2, Its_Name, &
    & ModuleName, Low1, Low2 )
    double precision, pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! First dimension of To_Allocate
    integer, intent(in) :: Dim2    ! Second dimension of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer, intent(in), optional :: Low1, Low2 ! Low bounds for dimensions
    integer :: MyLow1, MyLow2
    integer :: STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    myLow1 = 1; myLow2 = 1
    if ( present(low1) ) myLow1 = low1
    if ( present(low2) ) myLow2 = low2
    allocate ( To_Allocate(myLow1:dim1,myLow2:dim2), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name // trim(bounds(myLow1,dim1,myLow2,dim2)) )
    if ( status == 0 .and. clearOnAllocate ) to_allocate = 0.0d0
    if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, 2*dim1*dim2 )
  end subroutine Allocate_Test_RealR8_2d
  ! ------------------------------------  Allocate_Test_RealR8_3d  -----
  subroutine Allocate_Test_RealR8_3d ( To_Allocate, Dim1, Dim2, Dim3, &
    & Its_Name, ModuleName, Low1, Low2, Low3 )
    double precision, pointer, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! First dimension of To_Allocate
    integer, intent(in) :: Dim2    ! Second dimension of To_Allocate
    integer, intent(in) :: Dim3    ! Third dimension of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer, intent(in), optional :: Low1, Low2, Low3 ! Low bounds for dimensions
    integer :: MyLow1, MyLow2, MyLow3
    integer :: STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    myLow1 = 1; myLow2 = 1; myLow3 = 1
    if ( present(low1) ) myLow1 = low1
    if ( present(low2) ) myLow2 = low2
    if ( present(low3) ) myLow3 = low3
    allocate ( To_Allocate(myLow1:dim1,myLow2:dim2,myLow3:dim3), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name // trim(bounds(myLow1,dim1,myLow2,dim2,myLow3,dim3)) )
    if ( status == 0 .and. clearOnAllocate ) to_allocate = 0.0d0
    if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, 2*dim1*dim2*dim3 )
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
      & MLSMSG_Allocate // Its_Name // trim(bounds(lowBound,dim1)) )
    if ( status == 0 .and. clearOnAllocate ) to_allocate = 0
    if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, dim1 )
  end subroutine Allocate_Test_Integer_1d
  ! -----------------------------------  Allocate_Test_Integer_2d  -----
  subroutine Allocate_Test_Integer_2d ( To_Allocate, Dim1, Dim2, Its_Name, &
    & ModuleName, LowBound_1, LowBound_2 )
    integer, pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! First dimension of To_Allocate
    integer, intent(in) :: Dim2    ! Second dimension of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer, intent(in), optional :: LowBound_1, LowBound_2 ! default 1
    integer :: MY_LOW_1, MY_LOW_2, STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    my_low_1 = 1
    if ( present(lowBound_1) ) my_low_1 = lowBound_1
    my_low_2 = 1
    if ( present(lowBound_2) ) my_low_2 = lowBound_2
    allocate ( To_Allocate(my_low_1:dim1,my_low_2:dim2), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name // trim(bounds(my_low_1,dim1,my_low_2,dim2)) )
    if ( status == 0 .and. clearOnAllocate ) to_allocate = 0
    if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, dim1*dim2 )
  end subroutine Allocate_Test_Integer_2d
  ! -----------------------------------  Allocate_Test_Integer_3d  -----
  subroutine Allocate_Test_Integer_3d ( To_Allocate, Dim1, Dim2, Dim3, Its_Name, &
    & ModuleName, LowBound_1 )
    integer, pointer, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! First dimension of To_Allocate
    integer, intent(in) :: Dim2    ! Second dimension of To_Allocate
    integer, intent(in) :: Dim3    ! Third dimension of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer, intent(in), optional :: LowBound_1 ! default 1
    integer :: MY_LOW_1, STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    my_low_1 = 1
    if ( present(lowBound_1) ) my_low_1 = lowBound_1
    allocate ( To_Allocate(my_low_1:dim1,dim2,dim3), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name // trim(bounds(my_low_1,dim1,1,dim2,1,dim3)) )
    if ( status == 0 .and. clearOnAllocate ) to_allocate = 0
    if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, dim1*dim2*dim3 )
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
      & MLSMSG_Allocate // Its_Name // trim(bounds(lowBound,dim1)) )
    if ( status == 0 .and. clearOnAllocate ) to_allocate = .false.
    if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, dim1/4+1 )
  end subroutine Allocate_Test_Logical_1d
  ! -------------------------------------  Allocate_Test_Logical_2d  -----
  subroutine Allocate_Test_Logical_2d ( To_Allocate, Dim1, Dim2, Its_Name, &
    & ModuleName, Low1, Low2 )
    logical, pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! First dimension of To_Allocate
    integer, intent(in) :: Dim2    ! Second dimension of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer, intent(in), optional :: Low1, Low2 ! Low bounds for dimensions
    integer :: MyLow1, MyLow2
    integer :: STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    myLow1 = 1; myLow2 = 1
    if ( present(low1) ) myLow1 = low1
    if ( present(low2) ) myLow2 = low2
    allocate ( To_Allocate(myLow1:dim1,myLow2:dim2), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name // trim(bounds(myLow1,dim1,myLow2,dim2)) )
    if ( status == 0 .and. clearOnAllocate ) to_allocate = .false.
    if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, (dim1*dim2)/4+1 )
  end subroutine Allocate_Test_Logical_2d
  ! --------------------------------------  Allocate_Test_RealR4_1d  -----
  subroutine Allocate_Test_RealR4_1d ( To_Allocate, Dim1, Its_Name, ModuleName, &
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
      & MLSMSG_Allocate // Its_Name // trim(bounds(lowBound,dim1)) )
    if ( status == 0 .and. clearOnAllocate ) to_allocate = 0.0
    if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, dim1 )
  end subroutine Allocate_Test_RealR4_1d
  ! --------------------------------------  Allocate_Test_RealR4_2d  -----
  subroutine Allocate_Test_RealR4_2d ( To_Allocate, Dim1, Dim2, Its_Name, &
    & ModuleName, Low1, Low2 )
    real, pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! First dimension of To_Allocate
    integer, intent(in) :: Dim2    ! Second dimension of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer, intent(in), optional :: Low1, Low2 ! Low bounds for dimensions
    integer :: MyLow1, MyLow2
    integer :: STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    myLow1 = 1; myLow2 = 1
    if ( present(low1) ) myLow1 = low1
    if ( present(low2) ) myLow2 = low2
    allocate ( To_Allocate(myLow1:dim1,myLow2:dim2), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name // trim(bounds(myLow1,dim1,myLow2,dim2)) )
    if ( status == 0 .and. clearOnAllocate ) to_allocate = 0.0
    if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, dim1*dim2 )
  end subroutine Allocate_Test_RealR4_2d
  ! ------------------------------------  Allocate_Test_RealR4_3d  -----
  subroutine Allocate_Test_RealR4_3d ( To_Allocate, Dim1, Dim2, Dim3, &
    & Its_Name, ModuleName, Low1, Low2, Low3 )
    real, pointer, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! First dimension of To_Allocate
    integer, intent(in) :: Dim2    ! Second dimension of To_Allocate
    integer, intent(in) :: Dim3    ! Third dimension of To_Allocate
    character(len=*), intent(in) :: Its_Name, ModuleName
    integer, intent(in), optional :: Low1, Low2, Low3 ! Low bounds for dimensions
    integer :: MyLow1, MyLow2, MyLow3
    integer :: STATUS
    call deallocate_Test ( To_Allocate, Its_Name, ModuleName )
    myLow1 = 1; myLow2 = 1; myLow3 = 1
    if ( present(low1) ) myLow1 = low1
    if ( present(low2) ) myLow2 = low2
    if ( present(low3) ) myLow3 = low3
    allocate ( To_Allocate(myLow1:dim1,myLow2:dim2,myLow3:dim3), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // Its_Name // trim(bounds(myLow1,dim1,myLow2,dim2,myLow3,dim3)) )
    if ( status == 0 .and. clearOnAllocate ) to_allocate = 0.0
    if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, dim1*dim2*dim3 )
  end subroutine Allocate_Test_RealR4_3d
  ! -------------------------------  Deallocate_Test_Character_1d  -----
  subroutine Deallocate_Test_Character_1d ( To_Deallocate, Its_Name, ModuleName )
    character(len=*), pointer, dimension(:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    integer :: DIM1
    if ( associated(To_Deallocate) ) then
      if ( trackAllocates ) then
        dim1 = size ( to_deallocate, 1 )
      endif
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      else if ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
      nullify ( to_Deallocate )
      if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, -dim1/4-1 )
    end if
  end subroutine Deallocate_Test_Character_1d
  ! -------------------------------  Deallocate_Test_Character_2d  -----
  subroutine Deallocate_Test_Character_2d ( To_Deallocate, Its_Name, ModuleName )
    character(len=*), pointer, dimension(:,:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    integer :: DIM1, DIM2
    if ( associated(To_Deallocate) ) then
      if ( trackAllocates ) then
        dim1 = size ( to_deallocate, 1 )
        dim2 = size ( to_deallocate, 2 )
      endif
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      else if ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
      nullify ( to_Deallocate )
      if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, -(dim1*dim2)/4-1 )
    end if
  end subroutine Deallocate_Test_Character_2d
  ! -------------------------------  Deallocate_Test_Character_3d  -----
  subroutine Deallocate_Test_Character_3d ( To_Deallocate, Its_Name, ModuleName )
    character(len=*), pointer, dimension(:,:,:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    integer :: DIM1, DIM2, DIM3
    if ( associated(To_Deallocate) ) then
      if ( trackAllocates ) then
        dim1 = size ( to_deallocate, 1 )
        dim2 = size ( to_deallocate, 2 )
        dim3 = size ( to_deallocate, 3 )
      endif
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      else if ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
      nullify ( to_Deallocate )
      if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, -(dim1*dim2*dim3)/4-1 )
    end if
  end subroutine Deallocate_Test_Character_3d
  ! ----------------------------------  Deallocate_Test_RealR8_1d  -----
  subroutine Deallocate_Test_RealR8_1d ( To_Deallocate, Its_Name, ModuleName )
    double precision, pointer, dimension(:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    integer :: DIM1
    if ( associated(To_Deallocate) ) then
      if ( trackAllocates ) then
        dim1 = size ( to_deallocate, 1 )
      endif
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      else if ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
      nullify ( to_Deallocate )
      if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, -2*dim1 )
    end if
  end subroutine Deallocate_Test_RealR8_1d
  ! ----------------------------------  Deallocate_Test_RealR8_2d  -----
  subroutine Deallocate_Test_RealR8_2d ( To_Deallocate, Its_Name, ModuleName )
    double precision, pointer, dimension(:,:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    integer :: DIM1, DIM2
    if ( associated(To_Deallocate) ) then
      if ( trackAllocates ) then
        dim1 = size ( to_deallocate, 1 )
        dim2 = size ( to_deallocate, 2 )
      endif
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      else if ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
      nullify ( to_Deallocate )
      if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, -2*dim1*dim2 )
    end if
  end subroutine Deallocate_Test_RealR8_2d
  ! ----------------------------------  Deallocate_Test_RealR8_3d  -----
  subroutine Deallocate_Test_RealR8_3d ( To_Deallocate, Its_Name, ModuleName )
    double precision, pointer, dimension(:,:,:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    integer :: DIM1, DIM2, DIM3
    if ( associated(To_Deallocate) ) then
      if ( trackAllocates ) then
        dim1 = size ( to_deallocate, 1 )
        dim2 = size ( to_deallocate, 2 )
        dim3 = size ( to_deallocate, 3 )
      endif
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      else if ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
      nullify ( to_Deallocate )
      if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, -2*dim1*dim2*dim3 )
    end if
  end subroutine Deallocate_Test_RealR8_3d
  ! ---------------------------------  Deallocate_Test_Integer_1d  -----
  subroutine Deallocate_Test_Integer_1d ( To_Deallocate, Its_Name, ModuleName )
    integer, pointer, dimension(:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    integer :: DIM1
    if ( associated(To_Deallocate) ) then
      if ( trackAllocates ) then
        dim1 = size ( to_deallocate, 1 )
      endif
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      else if ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
      nullify ( to_Deallocate )
      if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, -dim1 )
    end if
  end subroutine Deallocate_Test_Integer_1d
  ! ---------------------------------  Deallocate_Test_Integer_2d  -----
  subroutine Deallocate_Test_Integer_2d ( To_Deallocate, Its_Name, ModuleName )
    integer, pointer, dimension(:,:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    integer :: DIM1, DIM2
    if ( associated(To_Deallocate) ) then
      if ( trackAllocates ) then
        dim1 = size ( to_deallocate, 1 )
        dim2 = size ( to_deallocate, 2 )
      endif
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      else if ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
      nullify ( to_Deallocate )
      if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, -dim1*dim2 )
    end if
  end subroutine Deallocate_Test_Integer_2d
  ! ---------------------------------  Deallocate_Test_Integer_3d  -----
  subroutine Deallocate_Test_Integer_3d ( To_Deallocate, Its_Name, ModuleName )
    integer, pointer, dimension(:,:,:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    integer :: DIM1, DIM2, DIM3
    if ( associated(To_Deallocate) ) then
      if ( trackAllocates ) then
        dim1 = size ( to_deallocate, 1 )
        dim2 = size ( to_deallocate, 2 )
        dim3 = size ( to_deallocate, 3 )
      endif
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      else if ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
      nullify ( to_Deallocate )
      if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, -dim1*dim2*dim3 )
    end if
  end subroutine Deallocate_Test_Integer_3d
  ! ---------------------------------  Deallocate_Test_Logical_1d  -----
  subroutine Deallocate_Test_Logical_1d ( To_Deallocate, Its_Name, ModuleName )
    logical, pointer, dimension(:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    integer :: DIM1
    if ( associated(To_Deallocate) ) then
      if ( trackAllocates ) then
        dim1 = size ( to_deallocate, 1 )
      endif
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      else if ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
      nullify ( to_Deallocate )
      if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, -dim1/4-1 )
    end if
  end subroutine Deallocate_Test_Logical_1d
  ! ---------------------------------  Deallocate_Test_Logical_1d  -----
  subroutine Deallocate_Test_Logical_2d ( To_Deallocate, Its_Name, ModuleName )
    logical, pointer, dimension(:,:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    integer :: DIM1, DIM2
    if ( associated(To_Deallocate) ) then
      if ( trackAllocates ) then
        dim1 = size ( to_deallocate, 1 )
        dim2 = size ( to_deallocate, 2 )
      endif
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      else if ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
      nullify ( to_Deallocate )
      if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, -(dim1*dim2)/4-1 )
    end if
  end subroutine Deallocate_Test_Logical_2d
  ! ------------------------------------  Deallocate_Test_RealR4_1d  -----
  subroutine Deallocate_Test_RealR4_1d ( To_Deallocate, Its_Name, ModuleName )
    real, pointer, dimension(:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    integer :: DIM1
    if ( associated(To_Deallocate) ) then
      if ( trackAllocates ) then
        dim1 = size ( to_deallocate, 1 )
      endif
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      else if ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
      nullify ( to_Deallocate )
      if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, -dim1 )
    end if
  end subroutine Deallocate_Test_RealR4_1d
  ! ------------------------------------  Deallocate_Test_RealR4_2d  -----
  subroutine Deallocate_Test_RealR4_2d ( To_Deallocate, Its_Name, ModuleName )
    real, pointer, dimension(:,:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    integer :: DIM1, DIM2
    if ( associated(To_Deallocate) ) then
      if ( trackAllocates ) then
        dim1 = size ( to_deallocate, 1 )
        dim2 = size ( to_deallocate, 2 )
      endif
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      else if ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
      nullify ( to_Deallocate )
      if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, -dim1*dim2 )
    end if
  end subroutine Deallocate_Test_RealR4_2d
  ! ----------------------------------  Deallocate_Test_RealR4_3d  -----
  subroutine Deallocate_Test_RealR4_3d ( To_Deallocate, Its_Name, ModuleName )
    real, pointer, dimension(:,:,:) :: To_Deallocate
    character(len=*) :: Its_Name, ModuleName
    integer :: STATUS
    integer :: DIM1, DIM2, DIM3
    if ( associated(To_Deallocate) ) then
      if ( trackAllocates ) then
        dim1 = size ( to_deallocate, 1 )
        dim2 = size ( to_deallocate, 2 )
        dim3 = size ( to_deallocate, 3 )
      endif
      deallocate ( To_Deallocate, stat=status )
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & MLSMSG_DeAllocate // Its_Name )
        dealloc_status = max(dealloc_status, status)
      else if ( collect_garbage_each_time ) then
        call mls_gc_now
      end if
      nullify ( to_Deallocate )
      if ( trackAllocates ) call ReportAllocateDeallocate ( its_name, moduleName, -dim1*dim2*dim3 )
    end if
  end subroutine Deallocate_Test_RealR4_3d

  ! ---------------------------------------  Set_garbage_collection  -----
  subroutine Set_garbage_collection ( setting )
    logical :: setting
    collect_garbage_each_time = setting
  end subroutine Set_garbage_collection

  ! =====  Private Procedures  ===========================================
  ! -------------------------------------------------------  Bounds  -----
  character(127) function Bounds ( Low1, High1, Low2, High2, Low3, High3 )
    integer, intent(in) :: Low1, High1
    integer, intent(in), optional :: Low2, High2, Low3, High3
    character(127) :: Temp
    write ( bounds, '("(",i0,":",i0)' ) low1, high1
    if ( present(low2) ) then
      write ( temp, '(",",i0,":",i0)' ) low2, high2
      bounds = trim(bounds) // temp
    end if
    if ( present(low3) ) then
      write ( temp, '(",",i0,":",i0)' ) low3, high3
      bounds = trim(bounds) // temp
    end if
    bounds = trim(bounds) // ')'
  end function Bounds

  ! ------------------------------------------------  Not_Used_Here  -----
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Allocate_Deallocate

! $Log$
! Revision 2.19  2004/12/28 00:24:39  vsnyder
! Use default REAL and DOUBLE PRECISION because R4 and R8 are not guaranteed
! to be different.
!
! Revision 2.18  2004/10/30 00:23:21  vsnyder
! Add low_bound1 to Allocate_Test_Integer_3d
!
! Revision 2.17  2004/10/02 02:42:35  vsnyder
! Add low bounds to Allocate_Test_Integer_2D
!
! Revision 2.16  2004/04/06 23:49:59  livesey
! Added ClearOnAllocate stuff
!
! Revision 2.15  2004/04/05 17:47:56  livesey
! Bug fix in the tracking stuff
!
! Revision 2.14  2004/04/03 05:43:53  livesey
! First attempt at memory tracking
!
! Revision 2.13  2004/01/30 21:20:35  vsnyder
! Add bounds to 'allocation failed' message
!
! Revision 2.12  2004/01/24 01:01:05  livesey
! Nullify all pointers after deallocating.
!
! Revision 2.11  2002/10/08 00:09:08  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.10  2002/07/01 23:47:08  vsnyder
! Add [De]AllocateTest_3d, cosmetic changes
!
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
