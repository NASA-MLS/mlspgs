! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module Allocate_Deallocate
!=============================================================================

! This module contains procedures to allocate, test the allocation,
! and announce an error if it fails, and similarly for deallocation,
! for Real, Double precision, Integer and Character arrays.

! **************************************************
! *****     Important Notice:                  *****
! **************************************************
! *****     All of the specific procedures     *****
! *****     of the Allocate_test generic       *****
! *****     deallocate the object to be        *****
! *****     allocated.  Therefore, it must     *****
! *****     be declared with => NULL() or      *****
! *****     nullified before its first use     *****
! **************************************************

  use MACHINE, only: MLS_GC_NOW
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error, MLSMSG_Warning

  implicit NONE
  private

  public :: ALLOCATE_TEST, DEALLOCATE_TEST, DEALLOC_STATUS, & 
    & SET_GARBAGE_COLLECTION, REPORTALLOCATEDEALLOCATE, &
    & Test_Allocate, Test_Deallocate

  interface ALLOCATE_TEST
    module procedure ALLOCATE_TEST_CHARACTER_1D
    module procedure ALLOCATE_TEST_CHARACTER_2D
    module procedure ALLOCATE_TEST_CHARACTER_3D
    module procedure ALLOCATE_TEST_COMPLEX_1D
    module procedure ALLOCATE_TEST_COMPLEX_2D
    module procedure ALLOCATE_TEST_COMPLEX_3D
    module procedure ALLOCATE_TEST_DCOMPLEX_1D
    module procedure ALLOCATE_TEST_DCOMPLEX_2D
    module procedure ALLOCATE_TEST_DCOMPLEX_3D
    module procedure ALLOCATE_TEST_INTEGER_1D, ALLOCATE_TEST_INTEGER_2D
    module procedure ALLOCATE_TEST_INTEGER_3D
    module procedure ALLOCATE_TEST_LOGICAL_1D, ALLOCATE_TEST_LOGICAL_2D
    module procedure ALLOCATE_TEST_LOGICAL_3D
    module procedure ALLOCATE_TEST_REALR4_1D, ALLOCATE_TEST_REALR4_2D
    module procedure ALLOCATE_TEST_REALR4_3D
    module procedure ALLOCATE_TEST_REALR8_1D, ALLOCATE_TEST_REALR8_2D
    module procedure ALLOCATE_TEST_REALR8_3D
  end interface

  interface DEALLOCATE_TEST
    module procedure DEALLOCATE_TEST_CHARACTER_1D
    module procedure DEALLOCATE_TEST_CHARACTER_2D
    module procedure DEALLOCATE_TEST_CHARACTER_3D
    module procedure DEALLOCATE_TEST_COMPLEX_1D
    module procedure DEALLOCATE_TEST_COMPLEX_2D
    module procedure DEALLOCATE_TEST_COMPLEX_3D
    module procedure DEALLOCATE_TEST_DCOMPLEX_1D
    module procedure DEALLOCATE_TEST_DCOMPLEX_2D
    module procedure DEALLOCATE_TEST_DCOMPLEX_3D
    module procedure DEALLOCATE_TEST_INTEGER_1D, DEALLOCATE_TEST_INTEGER_2D
    module procedure DEALLOCATE_TEST_INTEGER_3D
    module procedure DEALLOCATE_TEST_LOGICAL_1D, DEALLOCATE_TEST_LOGICAL_2D
    module procedure DEALLOCATE_TEST_LOGICAL_3D
    module procedure DEALLOCATE_TEST_REALR4_1D, DEALLOCATE_TEST_REALR4_2D
    module procedure DEALLOCATE_TEST_REALR4_3D
    module procedure DEALLOCATE_TEST_REALR8_1D, DEALLOCATE_TEST_REALR8_2D
    module procedure DEALLOCATE_TEST_REALR8_3D
  end interface

  integer, save :: DEALLOC_STATUS = 0

  logical, public :: CLEARONALLOCATE = .false. ! If true, zero all allocated stuff
  logical, public :: TRACKALLOCATES = .false. ! If true keep track of memory allocated
                                              ! and print every transaction
  ! and report on it.

  ! Element sizes (bytes)
  integer, parameter, public :: E_Ch = 1 ! Character
    ! Default integer, logical and real:
  integer, parameter, public :: E_Def = (bit_size(e_ch) + 7 ) / 8
  integer, parameter, public :: E_DP = 2 * e_def ! DP and Complex

  logical, save, private :: COLLECT_GARBAGE_EACH_TIME = .false.

  ! The next two used for tracking allocated memory
  ! (The 1st is public to enable reporting finer or coarser grains)
  real, save, public  :: MEMORY_UNITS = 1024. ! Report nothing smaller than KB
  double precision, save, public :: NoBytesAllocated=0.0d0 ! Number of MEMORY_UNITS allocated.


  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), parameter, private :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here 
  !-----------------------------------------------------------------------------

  interface ReportAllocateDeallocate
    module procedure ReportAllocateDeallocate_ints
    module procedure ReportAllocateDeallocate_real
  end interface

  interface Test_Deallocate
    module procedure Test_Deallocate_int_s
    module procedure Test_Deallocate_real_s
  end interface

contains
  ! =====     Public Procedures      ===================================

  !-----------------------------------   ReportAllocateDeallocate_ints  -----
  subroutine ReportAllocateDeallocate_ints ( name, moduleName, noBytes )
    ! Dummy arguments
    character (len=*), intent(in) :: NAME ! Name of thing allocated
    character (len=*), intent(in) :: MODULENAME ! Module that allocated it
    integer, intent(in) :: noBytes      ! No bytes allocated (or deallocated if -ve)
    !
    call ReportAllocateDeallocate( name, moduleName, real(noBytes) )
  end subroutine ReportAllocateDeallocate_ints

  !-----------------------------------   ReportAllocateDeallocate_real  -----
  subroutine ReportAllocateDeallocate_real ( name, moduleName, noBytes )
    use Output_m, only: OUTPUT, DUMPSIZE
    ! Dummy arguments
    character (len=*), intent(in) :: NAME ! Name of thing allocated
    character (len=*), intent(in) :: MODULENAME ! Module that allocated it
    real, intent(in) :: noBytes      ! No bytes allocated (or deallocated if -ve)

    ! Executable code
    if ( .not. trackAllocates ) return        ! Most probably will not be called anyway
    ! print *, 'noBytes: ', noBytes
    noBytesAllocated = noBytesAllocated + noBytes
    call output ( 'Tracking: ' )
    if ( noBytes < 0.0 ) then
      call output ( 'Dea' )
    else
      call output ( 'A' )
    end if
    call output ( 'llocated ' )
    call DumpSize ( abs ( noBytes ), units=MEMORY_UNITS )
    call output ( ' for ' // trim ( name ) // ' in ' )
    if ( moduleName(1:1) == '$' ) then
      ! The moduleNameIn is <dollar>RCSFile: <filename>,v <dollar>
      call output ( moduleName(11:(len_trim(moduleName)-8)) )
    else
      call output ( moduleName )
    end if
    call output ( ' total ' )
    call DumpSize ( noBytesAllocated, units=MEMORY_UNITS, advance='yes' )
  end subroutine ReportAllocateDeallocate_real

  ! -------------------------------------  Set_garbage_collection  -----
  subroutine Set_garbage_collection ( setting )
    logical :: setting
    collect_garbage_each_time = setting
  end subroutine Set_garbage_collection

  ! ----------------------------------------------  Test_Allocate  -----
  subroutine Test_Allocate ( Status, ModuleNameIn, ItsName, lBounds, uBounds, &
    & ElementSize )
  ! Test the status from an allocate.  If it's nonzero, issue a message.
  ! Track allocations if TrackAllocates is true and ElementSize is present
  ! and > 0.  
    integer, intent(in) :: Status
    character(len=*), intent(in) :: ModuleNameIn, ItsName
    integer, intent(in) :: Lbounds(:), Ubounds(:)
    integer, intent(in), optional :: ElementSize ! Bytes, <= 0 for no tracking
    real :: Amount
    character(127) :: Bounds
    integer :: I, L

    if ( status /= 0 ) then
      ! print *, 'status ', status
      write ( bounds, '("(",i0,":",i0, 2(:",",i0,":",i0))' ) &
        & ( lBounds(i), uBounds(i), i = 1, size(lBounds) )
      l = len_trim(bounds)+1
      bounds(l:l)= ')'
      call MLSMessage ( MLSMSG_Error, moduleNameIn, &
        & MLSMSG_Allocate // ItsName  // bounds(:l) )
    end if

    if ( .not. present(elementSize) ) return

    if ( elementSize > 0 ) then
      amount = memproduct(elementSize, ubounds-lbounds+1)
      if ( trackAllocates .and. present(elementSize) ) then
        call ReportAllocateDeallocate ( itsName, moduleNameIn, amount )
      else
        noBytesAllocated = noBytesAllocated + amount
      end if
    end if

  end subroutine Test_Allocate

  ! --------------------------------------------  Test_DeAllocate_int_s  -----
  subroutine Test_DeAllocate_int_s ( Status, ModuleNameIn, ItsName, Size )
  ! Test the status from a deallocate.  If it's nonzero, issue a message.
  ! Do garbage collection if Collect_garbage_each_time is true.
  ! Track deallocations if TrackAllocates is true and Size is present and > 0.
    integer, intent(in) :: Status
    character(len=*), intent(in) :: ModuleNameIn, ItsName
    integer, intent(in) :: Size ! in MEMORY_UNITS, <= 0 for no tracking
    call Test_Deallocate( Status, ModuleNameIn, ItsName, real(Size) )
  end subroutine Test_DeAllocate_int_s

  ! --------------------------------------------  Test_Deallocate_real_s  -----
  subroutine Test_Deallocate_real_s ( Status, ModuleNameIn, ItsName, Size )
  ! Test the status from a deallocate.  If it's nonzero, issue a message.
  ! Do garbage collection if Collect_garbage_each_time is true.
  ! Track deallocations if TrackAllocates is true and Size is present and > 0.
    integer, intent(in) :: Status
    character(len=*), intent(in) :: ModuleNameIn, ItsName
    real, intent(in), optional :: Size ! in MEMORY_UNITS, <= 0 for no tracking

    if ( status /= 0 ) then
      call MLSMessage ( MLSMSG_Warning, moduleNameIn, &
        & MLSMSG_DeAllocate // itsName )
      dealloc_status = max(dealloc_status, status)
    else if ( collect_garbage_each_time ) then
      call mls_gc_now
    end if
    if ( status == 0 .and. present(size) ) then
      if ( size > 0.0 ) then
        if ( trackAllocates ) then
          call ReportAllocateDeallocate ( itsName, moduleNameIn, -size )
        else
          noBytesAllocated = noBytesAllocated - size
        end if
      end if
    end if

  end subroutine Test_Deallocate_real_s

  ! =====     Private Procedures     ===================================
  ! ---------------------------------  Allocate_Test_Character_1d  -----
  subroutine Allocate_Test_Character_1d ( To_Allocate, Dim1, &
    & ItsName, ModuleName, LowBound, Fill )
    character(len=*), pointer, dimension(:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: LowBound     ! Lower bound, default 1
    character(len=*), intent(in), optional :: Fill
    character(len=1), parameter :: Default = ''
    integer, parameter :: ESize = E_ch
    include "Allocate_Test_1D.f9h"
  end subroutine Allocate_Test_Character_1d
  ! ---------------------------------  Allocate_Test_Character_2d  -----
  subroutine Allocate_Test_Character_2d ( To_Allocate, Dim1, Dim2, &
    & ItsName, ModuleName, Low1, Low2, Fill )
    character(len=*), pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2 ! default 1
    character(len=*), intent(in), optional :: Fill ! To fill allocated array
    character(len=1), parameter :: Default = ''
    integer, parameter :: ESize = E_ch
    include "Allocate_Test_2D.f9h"
  end subroutine Allocate_Test_Character_2d
  ! ---------------------------------  Allocate_Test_Character_3d  -----
  subroutine Allocate_Test_Character_3d ( To_Allocate, Dim1, Dim2, Dim3, &
    & ItsName, ModuleName, Low1, Low2, Low3, Fill )
    character(len=*), pointer, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    integer, intent(in) :: Dim3    ! Upper bound of third dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2, Low3 ! Low bounds for dimensions
    character(len=*), intent(in), optional :: Fill ! To fill allocated array
    character(len=1), parameter :: Default = ''
    integer, parameter :: ESize = E_ch
    include "Allocate_Test_3D.f9h"
  end subroutine Allocate_Test_Character_3d
  ! -----------------------------------  Allocate_Test_Complex_1d  -----
  subroutine Allocate_Test_Complex_1d ( To_Allocate, Dim1, &
    & ItsName, ModuleName, LowBound, Fill )
    Complex, pointer, dimension(:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: LowBound     ! Lower bound, default 1
    Complex, intent(in), optional :: Fill
    Complex, parameter :: Default = (0.0,0.0)
    integer, parameter :: ESize = 2 * e_def
    include "Allocate_Test_1D.f9h"
  end subroutine Allocate_Test_Complex_1d
  ! -----------------------------------  Allocate_Test_Complex_2d  -----
  subroutine Allocate_Test_Complex_2d ( To_Allocate, Dim1, Dim2, &
    & ItsName, ModuleName, Low1, Low2, Fill )
    Complex, pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2 ! default 1
    Complex, intent(in), optional :: Fill ! To fill allocated array
    Complex, parameter :: Default = (0.0,0.0)
    integer, parameter :: ESize = 2 * e_def
    include "Allocate_Test_2D.f9h"
  end subroutine Allocate_Test_Complex_2d
  ! -----------------------------------  Allocate_Test_Complex_3d  -----
  subroutine Allocate_Test_Complex_3d ( To_Allocate, Dim1, Dim2, Dim3, &
    & ItsName, ModuleName, Low1, Low2, Low3, Fill )
    Complex, pointer, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    integer, intent(in) :: Dim3    ! Upper bound of third dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2, Low3 ! Low bounds for dimensions
    Complex, intent(in), optional :: Fill ! To fill allocated array
    Complex, parameter :: Default = (0.0,0.0)
    integer, parameter :: ESize = 2 * e_def
    include "Allocate_Test_3D.f9h"
  end subroutine Allocate_Test_Complex_3d
  ! ----------------------------------  Allocate_Test_DComplex_1d  -----
  subroutine Allocate_Test_DComplex_1d ( To_Allocate, Dim1, &
    & ItsName, ModuleName, LowBound, Fill )
    Complex(kind(0.0d0)), pointer, dimension(:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: LowBound     ! Lower bound, default 1
    Complex(kind(0.0d0)), intent(in), optional :: Fill
    Complex(kind(0.0d0)), parameter :: Default = (0.0d0,0.0d0)
    integer, parameter :: ESize = 2 * e_dp
    include "Allocate_Test_1D.f9h"
  end subroutine Allocate_Test_DComplex_1d
  ! ----------------------------------  Allocate_Test_DComplex_2d  -----
  subroutine Allocate_Test_DComplex_2d ( To_Allocate, Dim1, Dim2, &
    & ItsName, ModuleName, Low1, Low2, Fill )
    Complex(kind(0.0d0)), pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2 ! default 1
    Complex(kind(0.0d0)), intent(in), optional :: Fill ! To fill allocated array
    Complex(kind(0.0d0)), parameter :: Default = (0.0d0,0.0d0)
    integer, parameter :: ESize = 2 * e_dp
    include "Allocate_Test_2D.f9h"
  end subroutine Allocate_Test_DComplex_2d
  ! ----------------------------------  Allocate_Test_DComplex_3d  -----
  subroutine Allocate_Test_DComplex_3d ( To_Allocate, Dim1, Dim2, Dim3, &
    & ItsName, ModuleName, Low1, Low2, Low3, Fill )
    Complex(kind(0.0d0)), pointer, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    integer, intent(in) :: Dim3    ! Upper bound of third dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2, Low3 ! Low bounds for dimensions
    Complex(kind(0.0d0)), intent(in), optional :: Fill ! To fill allocated array
    Complex(kind(0.0d0)), parameter :: Default = (0.0d0,0.0d0)
    integer, parameter :: ESize = 2 * e_dp
    include "Allocate_Test_3D.f9h"
  end subroutine Allocate_Test_DComplex_3d
  ! ------------------------------------  Allocate_Test_RealR8_1d  -----
  subroutine Allocate_Test_RealR8_1d ( To_Allocate, Dim1, &
    & ItsName, ModuleName, LowBound, Fill )
    double precision, pointer, dimension(:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: LowBound     ! Lower bound, default 1
    double precision, intent(in), optional :: Fill
    double precision, parameter :: Default = 0.0d0
    integer, parameter :: ESize = E_dp
    include "Allocate_Test_1D.f9h"
  end subroutine Allocate_Test_RealR8_1d
  ! ------------------------------------  Allocate_Test_RealR8_2d  -----
  subroutine Allocate_Test_RealR8_2d ( To_Allocate, Dim1, Dim2, &
    & ItsName, ModuleName, Low1, Low2, Fill )
    double precision, pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2 ! Low bounds for dimensions
    double precision, intent(in), optional :: Fill
    double precision, parameter :: Default = 0.0d0
    integer, parameter :: ESize = E_dp
    include "Allocate_Test_2D.f9h"
  end subroutine Allocate_Test_RealR8_2d
  ! ------------------------------------  Allocate_Test_RealR8_3d  -----
  subroutine Allocate_Test_RealR8_3d ( To_Allocate, Dim1, Dim2, Dim3, &
    & ItsName, ModuleName, Low1, Low2, Low3, Fill )
    double precision, pointer, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    integer, intent(in) :: Dim3    ! Upper bound of third dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2, Low3 ! Low bounds for dimensions
    double precision, intent(in), optional :: Fill
    double precision, parameter :: Default = 0.0d0
    integer, parameter :: ESize = E_dp
    include "Allocate_Test_3D.f9h"
  end subroutine Allocate_Test_RealR8_3d
  ! -----------------------------------  Allocate_Test_Integer_1d  -----
  subroutine Allocate_Test_Integer_1d ( To_Allocate, Dim1, &
    & ItsName, ModuleName, LowBound, Fill )
    integer, pointer, dimension(:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: LowBound     ! Lower bound, default 1
    integer, intent(in), optional :: Fill ! To fill allocated array
    integer, parameter :: Default = 0
    integer, parameter :: ESize = E_def
    include "Allocate_Test_1D.f9h"
  end subroutine Allocate_Test_Integer_1d
  ! -----------------------------------  Allocate_Test_Integer_2d  -----
  subroutine Allocate_Test_Integer_2d ( To_Allocate, Dim1, Dim2, &
    & ItsName, ModuleName, Low1, Low2, Fill )
    integer, pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2 ! default 1
    integer, intent(in), optional :: Fill ! To fill allocated array
    integer, parameter :: Default = 0
    integer, parameter :: ESize = E_def
    include "Allocate_Test_2D.f9h"
  end subroutine Allocate_Test_Integer_2d
  ! -----------------------------------  Allocate_Test_Integer_3d  -----
  subroutine Allocate_Test_Integer_3d ( To_Allocate, Dim1, Dim2, Dim3, &
    & ItsName, ModuleName, Low1, Low2, Low3, Fill )
    integer, pointer, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    integer, intent(in) :: Dim3    ! Upper bound of third dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2, Low3 ! default 1
    integer, intent(in), optional :: Fill ! To fill allocated array
    integer, parameter :: Default = 0
    integer, parameter :: ESize = E_def
    include "Allocate_Test_3D.f9h"
  end subroutine Allocate_Test_Integer_3d
  ! -----------------------------------  Allocate_Test_Logical_1d  -----
  subroutine Allocate_Test_Logical_1d ( To_Allocate, Dim1, &
    & ItsName, ModuleName, LowBound, Fill )
    logical, pointer, dimension(:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: LowBound     ! Lower bound, default 1
    logical, intent(in), optional :: Fill ! To fill allocated array
    logical, parameter :: Default = .false.
    integer, parameter :: ESize = E_def
    include "Allocate_Test_1D.f9h"
  end subroutine Allocate_Test_Logical_1d
  ! -------------------------------------  Allocate_Test_Logical_2d  -----
  subroutine Allocate_Test_Logical_2d ( To_Allocate, Dim1, Dim2, &
    & ItsName, ModuleName, Low1, Low2, Fill )
    logical, pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2 ! Low bounds for dimensions
    logical, intent(in), optional :: Fill
    logical, parameter :: Default = .false.
    integer, parameter :: ESize = e_def
    include "Allocate_Test_2D.f9h"
  end subroutine Allocate_Test_Logical_2d
  ! -------------------------------------  Allocate_Test_Logical_3d  -----
  subroutine Allocate_Test_Logical_3d ( To_Allocate, Dim1, Dim2, Dim3, &
    & ItsName, ModuleName, Low1, Low2, Low3, Fill )
    logical, pointer, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    integer, intent(in) :: Dim3    ! Upper bound of third dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2, Low3 ! Low bounds for dimensions
    logical, intent(in), optional :: Fill
    logical, parameter :: Default = .false.
    integer, parameter :: ESize = e_def
    include "Allocate_Test_3D.f9h"
  end subroutine Allocate_Test_Logical_3d
  ! --------------------------------------  Allocate_Test_RealR4_1d  -----
  subroutine Allocate_Test_RealR4_1d ( To_Allocate, Dim1, &
    & ItsName, ModuleName, LowBound, Fill )
    real, pointer, dimension(:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: LowBound     ! Lower bound, default 1
    real, intent(in), optional :: Fill
    real, parameter :: Default = 0.0
    integer, parameter :: ESize = e_def
    include "Allocate_Test_1D.f9h"
  end subroutine Allocate_Test_RealR4_1d
  ! --------------------------------------  Allocate_Test_RealR4_2d  -----
  subroutine Allocate_Test_RealR4_2d ( To_Allocate, Dim1, Dim2, &
    & ItsName, ModuleName, Low1, Low2, Fill )
    real, pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2 ! Low bounds for dimensions
    real, intent(in), optional :: Fill
    real, parameter :: Default = 0.0
    integer, parameter :: ESize = e_def
    include "Allocate_Test_2D.f9h"
  end subroutine Allocate_Test_RealR4_2d
  ! ------------------------------------  Allocate_Test_RealR4_3d  -----
  subroutine Allocate_Test_RealR4_3d ( To_Allocate, Dim1, Dim2, Dim3, &
    & ItsName, ModuleName, Low1, Low2, Low3, Fill )
    real, pointer, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    integer, intent(in) :: Dim3    ! Upper bound of third dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2, Low3 ! Low bounds for dimensions
    real, intent(in), optional :: Fill
    real, parameter :: Default = 0.0
    integer, parameter :: ESize = e_def
    include "Allocate_Test_3D.f9h"
  end subroutine Allocate_Test_RealR4_3d
  ! -------------------------------  Deallocate_Test_Character_1d  -----
  subroutine Deallocate_Test_Character_1d ( To_Deallocate, ItsName, ModuleName )
    character(len=*), pointer, dimension(:) :: To_Deallocate
    integer, parameter :: ESize = e_ch
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Character_1d
  ! -------------------------------  Deallocate_Test_Character_2d  -----
  subroutine Deallocate_Test_Character_2d ( To_Deallocate, ItsName, ModuleName )
    character(len=*), pointer, dimension(:,:) :: To_Deallocate
    integer, parameter :: ESize = e_ch
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Character_2d
  ! -------------------------------  Deallocate_Test_Character_3d  -----
  subroutine Deallocate_Test_Character_3d ( To_Deallocate, ItsName, ModuleName )
    character(len=*), pointer, dimension(:,:,:) :: To_Deallocate
    integer, parameter :: ESize = e_ch
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Character_3d
  ! ---------------------------------  Deallocate_Test_Complex_1d  -----
  subroutine Deallocate_Test_Complex_1d ( To_Deallocate, ItsName, ModuleName )
    complex, pointer, dimension(:) :: To_Deallocate
    integer, parameter :: ESize = 2 * e_def
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Complex_1d
  ! ---------------------------------  Deallocate_Test_Complex_2d  -----
  subroutine Deallocate_Test_Complex_2d ( To_Deallocate, ItsName, ModuleName )
    complex, pointer, dimension(:,:) :: To_Deallocate
    integer, parameter :: ESize = 2 * e_def
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Complex_2d
  ! ---------------------------------  Deallocate_Test_Complex_3d  -----
  subroutine Deallocate_Test_Complex_3d ( To_Deallocate, ItsName, ModuleName )
    complex, pointer, dimension(:,:,:) :: To_Deallocate
    integer, parameter :: ESize = 2 * e_def
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Complex_3d
  ! --------------------------------  Deallocate_Test_Dcomplex_1d  -----
  subroutine Deallocate_Test_Dcomplex_1d ( To_Deallocate, ItsName, ModuleName )
    complex(kind(0.0d0)), pointer, dimension(:) :: To_Deallocate
    integer, parameter :: ESize = 2 * e_dp
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Dcomplex_1d
  ! --------------------------------  Deallocate_Test_Dcomplex_2d  -----
  subroutine Deallocate_Test_Dcomplex_2d ( To_Deallocate, ItsName, ModuleName )
    complex(kind(0.0d0)), pointer, dimension(:,:) :: To_Deallocate
    integer, parameter :: ESize = 2 * e_dp
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Dcomplex_2d
  ! --------------------------------  Deallocate_Test_Dcomplex_3d  -----
  subroutine Deallocate_Test_Dcomplex_3d ( To_Deallocate, ItsName, ModuleName )
    complex(kind(0.0d0)), pointer, dimension(:,:,:) :: To_Deallocate
    integer, parameter :: ESize = 2 * e_dp
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Dcomplex_3d
  ! ----------------------------------  Deallocate_Test_RealR8_1d  -----
  subroutine Deallocate_Test_RealR8_1d ( To_Deallocate, ItsName, ModuleName )
    double precision, pointer, dimension(:) :: To_Deallocate
    integer, parameter :: ESize = e_dp
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_RealR8_1d
  ! ----------------------------------  Deallocate_Test_RealR8_2d  -----
  subroutine Deallocate_Test_RealR8_2d ( To_Deallocate, ItsName, ModuleName )
    double precision, pointer, dimension(:,:) :: To_Deallocate
    integer, parameter :: ESize = e_dp
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_RealR8_2d
  ! ----------------------------------  Deallocate_Test_RealR8_3d  -----
  subroutine Deallocate_Test_RealR8_3d ( To_Deallocate, ItsName, ModuleName )
    double precision, pointer, dimension(:,:,:) :: To_Deallocate
    integer, parameter :: ESize = e_dp
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_RealR8_3d
  ! ---------------------------------  Deallocate_Test_Integer_1d  -----
  subroutine Deallocate_Test_Integer_1d ( To_Deallocate, ItsName, ModuleName )
    integer, pointer, dimension(:) :: To_Deallocate
    integer, parameter :: ESize = e_def
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Integer_1d
  ! ---------------------------------  Deallocate_Test_Integer_2d  -----
  subroutine Deallocate_Test_Integer_2d ( To_Deallocate, ItsName, ModuleName )
    integer, pointer, dimension(:,:) :: To_Deallocate
    integer, parameter :: ESize = e_def
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Integer_2d
  ! ---------------------------------  Deallocate_Test_Integer_3d  -----
  subroutine Deallocate_Test_Integer_3d ( To_Deallocate, ItsName, ModuleName )
    integer, pointer, dimension(:,:,:) :: To_Deallocate
    integer, parameter :: ESize = e_def
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Integer_3d
  ! ---------------------------------  Deallocate_Test_Logical_1d  -----
  subroutine Deallocate_Test_Logical_1d ( To_Deallocate, ItsName, ModuleName )
    logical, pointer, dimension(:) :: To_Deallocate
    integer, parameter :: ESize = e_def
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Logical_1d
  ! ---------------------------------  Deallocate_Test_Logical_2d  -----
  subroutine Deallocate_Test_Logical_2d ( To_Deallocate, ItsName, ModuleName )
    logical, pointer, dimension(:,:) :: To_Deallocate
    integer, parameter :: ESize = e_def
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Logical_2d
  ! ---------------------------------  Deallocate_Test_Logical_3d  -----
  subroutine Deallocate_Test_Logical_3d ( To_Deallocate, ItsName, ModuleName )
    logical, pointer, dimension(:,:,:) :: To_Deallocate
    integer, parameter :: ESize = e_def
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Logical_3d
  ! ------------------------------------  Deallocate_Test_RealR4_1d  -----
  subroutine Deallocate_Test_RealR4_1d ( To_Deallocate, ItsName, ModuleName )
    real, pointer, dimension(:) :: To_Deallocate
    integer, parameter :: ESize = e_def
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_RealR4_1d
  ! ------------------------------------  Deallocate_Test_RealR4_2d  -----
  subroutine Deallocate_Test_RealR4_2d ( To_Deallocate, ItsName, ModuleName )
    real, pointer, dimension(:,:) :: To_Deallocate
    integer, parameter :: ESize = e_def
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_RealR4_2d
  ! ----------------------------------  Deallocate_Test_RealR4_3d  -----
  subroutine Deallocate_Test_RealR4_3d ( To_Deallocate, ItsName, ModuleName )
    real, pointer, dimension(:,:,:) :: To_Deallocate
    integer, parameter :: ESize = e_def
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_RealR4_3d

  ! ----------------------------------  memproduct  -----
  function memproduct ( elementSize, dimensions ) result( p )
    ! Find how many multiples of MEMORY_UNITS an array
    ! with elementSize and dimensions
    ! Note: we go through all this to forestall integer overflows
    integer, dimension(:), intent(in) :: dimensions
    integer, intent(in)               :: elementSize
    real                              :: p
    !
    p = ( elementSize / MEMORY_UNITS ) * product(real(dimensions))
    ! print *, 'p ', p
  end function memproduct

  ! ------------------------------------------------  Not_Used_Here  -----
  logical function not_used_here()
  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), parameter :: Idparm = &
    & "$Id$"
  character(len=len(idparm)) :: Id = idparm
  !-----------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Allocate_Deallocate

! $Log$
! Revision 2.29  2006/07/28 01:57:25  vsnyder
! Use real() instead of *1.0 to convert to real
!
! Revision 2.28  2006/07/19 22:24:18  vsnyder
! Track memory usage even if not reporting, for trace_m
!
! Revision 2.27  2005/12/16 23:25:37  pwagner
! dumpSize moved from dump0 to output_m
!
! Revision 2.26  2005/11/18 21:09:35  vsnyder
! Don't print trailing comma after last dimension in allocation failure message
!
! Revision 2.25  2005/10/03 18:04:51  pwagner
! Allocated memory now tracked in units of MEMORY_UNITS
!
! Revision 2.24  2005/09/16 23:38:20  vsnyder
! Add Complex allocators
!
! Revision 2.23  2005/07/20 16:22:52  pwagner
! Declared undeclared variable L; compiles successfully
!
! Revision 2.22  2005/07/20 01:32:55  vsnyder
! Add Test_Allocate and Test_Deallocate.  Regularize routines.  All now have an
! optional Fill argument.  For a given rank, calling sequences are the same
! (modulo types), including argument names.  Add [De]Allocate_Test_Logical_3d.
! Move guts of [De]Allocate_Test* into include files.
!
! Revision 2.21  2005/06/03 01:53:06  vsnyder
! New copyright notice, move Id to not_used_here to avoid cascades
!
! Revision 2.20  2005/06/01 02:30:30  vsnyder
! Add optional Fill argument to allocate_test_integer_...
!
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
