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

  use HighOutput, only: FinalMemoryReport
  use, Intrinsic :: Ieee_Arithmetic, only: Ieee_Signaling_Nan, Ieee_Value, &
    & Ieee_Support_Datatype
  use Iso_C_Binding, only: C_Intptr_T, C_Loc
  use Machine, only: MLS_Gc_Now
  use MLSCommon, only: NoBlocksAllocated, NoBlocksDeAllocated, &
    & NoBytesAllocated, TotalAllocated, TotalDeAllocated
  use Printit_M, only: MLSMSG_Allocate, MLSMSG_Deallocate, &
    & MLSMSG_Error, MLSMSG_Warning, Printitout, Sniprcsfrom
  use Track_M, only: Trackallocate, Trackdeallocate

  implicit none
!   private

  public :: Byte_Size, Bytes, Allocate_Test, Deallocate_Test, dealloc_status, & 
    & FinalMemoryReport, NoBytesAllocated, Reportallocatedeallocate, &
    & Same_Shape, Set_garbage_collection, Test_Allocate, Test_Deallocate

  interface Allocate_Test
    ! For separate scalar arguments for bounds
    module procedure Allocate_Test_Character_1D, Allocate_Test_Character_2D
    module procedure Allocate_Test_Character_3D
    module procedure Allocate_Test_Complex_1D, Allocate_Test_Complex_2D
    module procedure Allocate_Test_Complex_3D
    module procedure Allocate_Test_DComplex_1D, Allocate_Test_DComplex_2D
    module procedure Allocate_Test_DComplex_3D
    module procedure Allocate_Test_Integer_1D, Allocate_Test_Integer_2D
    module procedure Allocate_Test_Integer_3D, Allocate_Test_Integer_4D
    module procedure Allocate_Test_Logical_1D, Allocate_Test_Logical_2D
    module procedure Allocate_Test_Logical_3D
    module procedure Allocate_Test_RealR4_1D, Allocate_Test_RealR4_2D
    module procedure Allocate_Test_RealR4_3D, Allocate_Test_RealR4_4D
    module procedure Allocate_Test_RealR8_1D, Allocate_Test_RealR8_2D
    module procedure Allocate_Test_RealR8_3D, Allocate_Test_RealR8_4D
    ! For allocatable instead of pointer, with scalar arguments for bounds
    module procedure AllocateA_Test_Character_1D, AllocateA_Test_Character_2D
    module procedure AllocateA_Test_Integer_1D, AllocateA_Test_Integer_2D
    module procedure AllocateA_Test_Integer_3D, AllocateA_Test_Integer_4D
    module procedure AllocateA_Test_Logical_1D, AllocateA_Test_Logical_2D
    module procedure AllocateA_Test_Logical_3D
    module procedure AllocateA_Test_RealR4_1D, AllocateA_Test_RealR4_2D
    module procedure AllocateA_Test_RealR4_3D, AllocateA_Test_RealR4_4D
    module procedure AllocateA_Test_RealR8_1D, AllocateA_Test_RealR8_2D
    module procedure AllocateA_Test_RealR8_3D, AllocateA_Test_RealR8_4D
    ! For array arguments for bounds
    module procedure Allocate_Test_Character_2D_A, Allocate_Test_Character_3D_A
    module procedure Allocate_Test_Complex_2D_A, Allocate_Test_Complex_3D_A
    module procedure Allocate_Test_DComplex_2D_A, Allocate_Test_DComplex_3D_A
    module procedure Allocate_Test_Integer_2D_A
    module procedure Allocate_Test_Integer_3D_A, Allocate_Test_Integer_4D_A
    module procedure Allocate_Test_Logical_2D_A, Allocate_Test_Logical_3D_A
    module procedure Allocate_Test_RealR4_2D_A
    module procedure Allocate_Test_RealR4_3D_A, Allocate_Test_RealR4_4D_A
    module procedure Allocate_Test_RealR8_2D_A
    module procedure Allocate_Test_RealR8_3D_A, Allocate_Test_RealR8_4D_A
    ! For allocatable instead of pointer, with array arguments for bounds
    module procedure AllocateA_Test_Integer_2D_A
    module procedure AllocateA_Test_Integer_3D_A, AllocateA_Test_Integer_4D_A
    module procedure AllocateA_Test_RealR4_2D_A
    module procedure AllocateA_Test_RealR4_3D_A, AllocateA_Test_RealR4_4D_A
    module procedure AllocateA_Test_RealR8_2D_A
    module procedure AllocateA_Test_RealR8_3D_A, AllocateA_Test_RealR8_4D_A
  end interface

  interface Deallocate_Test
    module procedure Deallocate_Test_Character_1D, Deallocate_Test_Character_2D
    module procedure Deallocate_Test_Character_3D
    module procedure Deallocate_Test_Complex_1D, Deallocate_Test_Complex_2D
    module procedure Deallocate_Test_Complex_3D
    module procedure Deallocate_Test_DComplex_1D, Deallocate_Test_DComplex_2D
    module procedure Deallocate_Test_DComplex_3D
    module procedure Deallocate_Test_Integer_1D, Deallocate_Test_Integer_2D
    module procedure Deallocate_Test_Integer_3D, Deallocate_Test_Integer_4D
    module procedure Deallocate_Test_Logical_1D, Deallocate_Test_Logical_2D
    module procedure Deallocate_Test_Logical_3D
    module procedure Deallocate_Test_RealR4_1D, Deallocate_Test_RealR4_2D
    module procedure Deallocate_Test_RealR4_3D, Deallocate_Test_RealR4_4D
    module procedure Deallocate_Test_RealR8_1D, Deallocate_Test_RealR8_2D
    module procedure Deallocate_Test_RealR8_3D, Deallocate_Test_RealR8_4D
    ! For allocatable instead of pointer argument
    module procedure DeallocateA_Test_Character_1D, DeallocateA_Test_Character_2D
    module procedure DeallocateA_Test_Integer_1D, DeallocateA_Test_Integer_2D
    module procedure DeallocateA_Test_Integer_3D, DeallocateA_Test_Integer_4D
    module procedure DeallocateA_Test_Logical_1D, DeallocateA_Test_Logical_2D
    module procedure DeallocateA_Test_Logical_3D
    module procedure DeallocateA_Test_RealR4_1D, DeallocateA_Test_RealR4_2D
    module procedure DeallocateA_Test_RealR4_3D, DeallocateA_Test_RealR4_4D
    module procedure DeallocateA_Test_RealR8_1D, DeallocateA_Test_RealR8_2D
    module procedure DeallocateA_Test_RealR8_3D, DeallocateA_Test_RealR8_4D
  end interface

  ! Size in bytes of an array
  interface Byte_Size
    module procedure Byte_Size_Character_1D, Byte_Size_Character_2D
    module procedure Byte_Size_Character_3D
    module procedure Byte_Size_Complex_1D, Byte_Size_Complex_2D
    module procedure Byte_Size_Complex_3D
    module procedure Byte_Size_DComplex_1D, Byte_Size_DComplex_2D
    module procedure Byte_Size_DComplex_3D
    module procedure Byte_Size_Integer_1D, Byte_Size_Integer_2D
    module procedure Byte_Size_Integer_3D, Byte_Size_Integer_4D
    module procedure Byte_Size_Logical_1D, Byte_Size_Logical_2D
    module procedure Byte_Size_Logical_3D
    module procedure Byte_Size_RealR4_1D, Byte_Size_RealR4_2D
    module procedure Byte_Size_RealR4_3D, Byte_Size_RealR4_4D
    module procedure Byte_Size_RealR4_5D, Byte_Size_RealR4_6D
    module procedure Byte_Size_RealR8_1D, Byte_Size_RealR8_2D
    module procedure Byte_Size_RealR8_3D, Byte_Size_RealR8_4D
    module procedure Byte_Size_RealR8_5D, Byte_Size_RealR8_6D
  end interface

  ! Size in bytes of a scalar; includes character length
  interface BYTES
    module procedure Bytes_Character_1D, Bytes_Character_2D
    module procedure Bytes_Character_3D
    module procedure Bytes_Complex_1D, Bytes_Complex_2D, Bytes_Complex_3D
    module procedure Bytes_DComplex_1D, Bytes_DComplex_2D, Bytes_DComplex_3D
    module procedure Bytes_Integer_1D, Bytes_Integer_2D
    module procedure Bytes_Integer_3D, Bytes_Integer_4D
    module procedure Bytes_Logical_1D, Bytes_Logical_2D, Bytes_Logical_3D
    module procedure Bytes_RealR4_1D, Bytes_RealR4_2D
    module procedure Bytes_RealR4_3D, Bytes_RealR4_4D
    module procedure Bytes_RealR4_5D, Bytes_RealR4_6D
    module procedure Bytes_RealR8_1D, Bytes_RealR8_2D
    module procedure Bytes_RealR8_3D, Bytes_RealR8_4D
    module procedure Bytes_RealR8_5D, Bytes_RealR8_6D
  end interface

  interface ReportAllocateDeallocate
    module procedure ReportAllocateDeallocate_ints
    module procedure ReportAllocateDeallocate_Real
  end interface

  ! subroutine Same_Shape ( Ref, New, ItsName, ModuleName )
  ! If Ref is not associated, deallocate New.  Otherwise, if Ref and New
  ! have the same shape (but not necessarily the same bounds), do nothing.
  ! Otherwise, allocate New with the same bounds as Ref.
  ! subroutine Same_Shape ( New, ItsName, ModuleName, TheShape )
  ! If New is not associated, or if any of its dimensions are different
  ! from TheShape, allocate New with the shape given by TheShape.
  ! The reason to put TheShape at the end is so as not to have a generic
  ! conflict with the integer versions.
  interface Same_Shape
    module procedure Same_Shape_Character_1D, Same_Shape_Character_2D
    module procedure Same_Shape_Character_3D
    module procedure Same_Shape_Complex_1D, Same_Shape_Complex_2D
    module procedure Same_Shape_Complex_3D
    module procedure Same_Shape_DComplex_1D, Same_Shape_DComplex_2D
    module procedure Same_Shape_DComplex_3D
    module procedure Same_Shape_Integer_1D, Same_Shape_Integer_2D
    module procedure Same_Shape_Integer_3D, Same_Shape_Integer_4D
    module procedure Same_Shape_Logical_1D, Same_Shape_Logical_2D
    module procedure Same_Shape_Logical_3D
    module procedure Same_Shape_RealR4_1D, Same_Shape_RealR4_2D
    module procedure Same_Shape_RealR4_3D, Same_Shape_RealR4_4D
    module procedure Same_Shape_RealR8_1D, Same_Shape_RealR8_2D
    module procedure Same_Shape_RealR8_3D, Same_Shape_RealR8_4D
    module procedure Same_Shape_Character_1D_A, Same_Shape_Character_2D_A
    module procedure Same_Shape_Character_3D_A
    module procedure Same_Shape_Complex_1D_A, Same_Shape_Complex_2D_A
    module procedure Same_Shape_Complex_3D_A
    module procedure Same_Shape_DComplex_1D_A, Same_Shape_DComplex_2D_A
    module procedure Same_Shape_DComplex_3D_A
    module procedure Same_Shape_Integer_1D_A, Same_Shape_Integer_2D_A
    module procedure Same_Shape_Integer_3D_A, Same_Shape_Integer_4D_A
    module procedure Same_Shape_Logical_1D_A, Same_Shape_Logical_2D_A
    module procedure Same_Shape_Logical_3D_A
    module procedure Same_Shape_RealR4_1D_A, Same_Shape_RealR4_2D_A
    module procedure Same_Shape_RealR4_3D_A, Same_Shape_RealR4_4D_A
    module procedure Same_Shape_RealR8_1D_A, Same_Shape_RealR8_2D_A
    module procedure Same_Shape_RealR8_3D_A, Same_Shape_RealR8_4D_A
  end interface

  interface Test_Allocate
    module procedure Test_Allocate_nd, Test_Allocate_1d
  end interface

  interface Test_Deallocate
    module procedure Test_Deallocate_int_s
  end interface

  integer, save :: Dealloc_Status = 0

  logical, public :: ClearOnAllocate = .false. ! If true, zero all allocated stuff
  logical, public :: InitRealNaN = .false.     ! If true, and FILL is not present,
                                               ! and IEEE_SupportDataType is true,
                                               ! fill real with signaling NaN.
  integer, public :: TrackAllocates = 0 ! <= 0 => No tracking
                                        ! == 1 => Track using the Track_m module
                                        ! == 2 => Report all transactions
                                        ! == 3 => Both track and report
  integer, public :: AllocateLogUnit = -1 ! If > 0, log allocates to this unit
  real, public ::    AllocateLogLimit = 1. ! Don't log amounts smaller than this

  ! Element sizes (bytes)
  integer, parameter, public :: E_Ch = ( storage_size(' ') + 7 ) / 8 ! Character
    ! Default integer, logical and real:
  integer, parameter, public :: E_Def = (storage_size(0.0) + 7 ) / 8
  integer, parameter, public :: E_DP = 2 * e_def ! DP and Complex

  logical, save, private :: COLLECT_GARBAGE_EACH_TIME = .false.

  double precision, parameter  :: OUTPUT_UNITS = 1.d6 ! output in MB
  logical, parameter :: DEEBUG = .false.

  integer(selected_int_kind(18)) :: Serial_No = 0 ! for logging

  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), parameter, private :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here 
  !-----------------------------------------------------------------------------

contains
  ! =====     Public Procedures      ===================================

  !-----------------------------------   ReportAllocateDeallocate_ints  -----
  subroutine ReportAllocateDeallocate_ints ( name, moduleName, noBytes, bounds )
    ! Dummy arguments
    character (len=*), intent(in) :: NAME ! Name of thing allocated
    character (len=*), intent(in) :: MODULENAME ! Module that allocated it
    integer, intent(in) :: noBytes      ! No bytes allocated (or deallocated if -ve)
    character (len=*), intent(in), optional :: Bounds ! for allocation

    call ReportAllocateDeallocate( name, moduleName, real(noBytes), bounds )
  end subroutine ReportAllocateDeallocate_ints

  !-----------------------------------   ReportAllocateDeallocate_Real  -----
  subroutine ReportAllocateDeallocate_Real ( name, moduleName, noBytes, bounds )
    use Output_M, only: Output
    use HighOutput, only: Dumpsize
    ! Dummy arguments
    character (len=*), intent(in) :: NAME ! Name of thing allocated
    character (len=*), intent(in) :: MODULENAME ! Module that allocated it
    real, intent(in) :: noBytes      ! No bytes allocated (or deallocated if -ve)
    character (len=*), intent(in), optional :: Bounds ! for allocation

    ! Executable code
    if ( trackAllocates <= 0 ) return        ! Most probably will not be called anyway
    ! print *, 'noBytes: ', noBytes
    if ( noBytes > 0 ) then
      TotalAllocated = TotalAllocated + noBytes
      NoBlocksAllocated = NoBlocksAllocated + 1
    elseif ( noBytes < 0 ) then
      TotalDeAllocated = TotalDeAllocated - noBytes
      NoBlocksDeAllocated = NoBlocksDeAllocated + 1
    endif
    noBytesAllocated = noBytesAllocated + noBytes
    call output ( 'Tracking: ' )
    if ( noBytes < 0.0 ) then
      call output ( 'Dea' )
    else
      call output ( 'A' )
    end if
    call output ( 'llocated ' )
    call DumpSize ( abs ( noBytes ) )
    call output ( ' for ' // trim ( name ) )
    if ( present(bounds) ) call output ( trim(bounds) )
    call output ( ' in ' )
    if ( moduleName(1:1) == '$' ) then
      ! The moduleNameIn is <dollar>RCSFile: <filename>,v <dollar>
      call output ( moduleName(11:(len_trim(moduleName)-8)) )
    else
      call output ( moduleName )
    end if
    call output ( ' total ' )
    call DumpSize ( noBytesAllocated, advance='yes' )
  end subroutine ReportAllocateDeallocate_Real

  ! -------------------------------------  Set_garbage_collection  -----
  subroutine Set_garbage_collection ( setting )
    logical :: setting
    collect_garbage_each_time = setting
  end subroutine Set_garbage_collection

  ! -------------------------------------------  Test_Allocate_nd  -----
  subroutine Test_Allocate_nd ( Status, ModuleNameIn, ItsName, &
    & lBounds, uBounds, ElementSize, ERMSG, Address )
    use HIGHOUTPUT, only: outputNamedValue
  ! Test the status from an allocate.  If it's nonzero, issue a message.
  ! Track allocations if TrackAllocates is >= 2 and ElementSize is present
  ! and > 0.
    integer, intent(in) :: Status
    character(len=*), intent(in) :: ModuleNameIn, ItsName
    integer, intent(in), optional :: Lbounds(:), Ubounds(:)
    integer, intent(in), optional :: ElementSize ! Bytes, <= 0 for no tracking
    character(*), intent(in), optional :: ERMSG  ! from Allocate statement
    integer(c_intptr_t), intent(in), optional :: Address ! of allocated object
                                                         ! if status == 0
    real :: Amount
    character(127) :: Bounds
    integer :: I, L

    if ( status /= 0 .or. present(ubounds) .and. &
      & present(elementSize) .and. trackAllocates >= 2 ) then
      ! print *, 'status ', status
      if ( present(ubounds) ) then
        if ( present(lbounds) ) then
          write ( bounds, '("(",i0,":",i0, 6(:",",i0,":",i0))' ) &
            & ( lBounds(i), uBounds(i), i = 1, size(lBounds) )
        else
          write ( bounds, '("(",i0,6(:",",i0))' ) ubounds
        end if
        l = len_trim(bounds)+1
        bounds(l:l)= ')'
        write ( bounds(l+1:), '(", status = ", i0)' ) status
        l = len_trim(bounds)
      else
        l = 0
        bounds(1:1) = ''
      end if
      if ( status /= 0 ) then
        if ( present(ermsg) ) then
          bounds(l+2:) = ermsg
          l = len_trim(bounds)
        end if
        call myMessage ( MLSMSG_Error, moduleNameIn, &
          & MLSMSG_Allocate // ItsName  // bounds(:l) )
      end if
    end if
    if ( AllocateLogUnit > 0 ) call LogAllocate( ModuleNameIn, ItsName, &
      & lBounds, uBounds, ElementSize, Address )

    if ( .not. present(elementSize) ) return

    if ( present(ubounds) ) then
      if ( elementSize > 0 ) then
        amount = 0
        if ( present(lbounds) ) then
          amount = memproduct(elementSize, ubounds-lbounds+1)
        else
          amount = memproduct(elementSize, ubounds)
        end if
        if ( trackAllocates >= 2 ) then
          call ReportAllocateDeallocate ( itsName, moduleNameIn, amount, bounds )
        else
          if ( amount > 0 ) then
            TotalAllocated = TotalAllocated + amount
            NoBlocksAllocated = NoBlocksAllocated + 1
          elseif ( amount < 0 ) then
            TotalDeAllocated = TotalDeAllocated - amount
            NoBlocksDeAllocated = NoBlocksDeAllocated + 1
          endif
          noBytesAllocated = noBytesAllocated + amount
          if ( DEEBUG ) &
            & call outputNamedValue( itsname, &
            & (/noBytesAllocated, amount*1.d0/)/OUTPUT_UNITS )
        end if
      end if
    end if

  end subroutine Test_Allocate_nd

  ! -------------------------------------------  Test_Allocate_1d  -----
  subroutine Test_Allocate_1d ( Status, ModuleNameIn, ItsName, &
    & lBounds, uBounds, ElementSize, ERMSG, Address )
  ! Test the status from an allocate.  If it's nonzero, issue a message.
  ! Track allocations if TrackAllocates is >= 2 and ElementSize is present
  ! and > 0.
  ! This is just a handy wrapper so that you don't need to put brackets
  ! around the bounds of a 1-d array just to make the bounds arrays.
    integer, intent(in) :: Status
    character(len=*), intent(in) :: ModuleNameIn, ItsName
    integer, intent(in), optional :: Lbounds
    integer, intent(in)           :: Ubounds
    integer, intent(in), optional :: ElementSize ! Bytes, <= 0 for no tracking
    character(*), intent(in), optional :: ERMSG  ! from Allocate statement
    integer(c_intptr_t), intent(in), optional :: Address ! of allocated object
                                                         ! if status == 0
    integer :: MyLBounds(1)
    integer :: MyUBounds(1)
    myLBounds = 1
    if ( present(lBounds) ) myLBounds = lBounds
    myUBounds = uBounds
    call test_allocate ( status, moduleNameIn, itsName, &
    & myLBounds, myUBounds, elementSize, ERMSG, address )
  end subroutine Test_Allocate_1d

  ! --------------------------------------------  Test_DeAllocate_int_s  -----
  subroutine Test_DeAllocate_int_s ( Status, ModuleNameIn, ItsName, Size, &
    & ERMSG, Address )
  ! Test the status from a deallocate.  If it's nonzero, issue a message.
  ! Do garbage collection if Collect_garbage_each_time is true.
  ! Track deallocations if TrackAllocates >= 2 and Size is present and > 0.
    use HIGHOUTPUT, only: outputNamedValue
    integer, intent(in) :: Status
    character(len=*), intent(in) :: ModuleNameIn, ItsName
    integer, intent(in), optional :: Size ! in Bytes
    character(*), intent(in), optional :: ERMSG ! From ERRMSG= specifier
    integer(c_intptr_t), intent(in), optional :: Address ! of allocated object
                                                         ! before deallocation

    integer :: L
    character(255) :: Line

    if ( status /= 0 ) then
      write ( line, '(", status = ", i0)' ) status
      if ( present(ermsg) ) then
        l = len_trim(line)
        line(l+2:) = ermsg
      end if
      call myMessage ( MLSMSG_Warning, moduleNameIn, &
        & MLSMSG_DeAllocate // itsName // trim(line) )
      dealloc_status = max(dealloc_status, status)
    else if ( collect_garbage_each_time ) then
      call mls_gc_now
    end if
    if ( AllocateLogUnit > 0 ) &
      & call LogDeallocate( ModuleNameIn, ItsName, Size, Address )
    if ( status == 0 .and. present(size) ) then
      if ( size > 0.0 ) then
        if ( trackAllocates >= 2 ) then
          call ReportAllocateDeallocate ( itsName, moduleNameIn, -size )
        else
          TotalDeAllocated = TotalDeAllocated - size
          NoBlocksDeAllocated = NoBlocksDeAllocated + 1
          noBytesAllocated = noBytesAllocated - size
          if ( DEEBUG ) &
            & call outputNamedValue( itsname, &
            & (/noBytesAllocated, size*1.d0/)/OUTPUT_UNITS )
        end if
      end if
    end if

  end subroutine Test_Deallocate_int_s

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
    include "Allocate_Test_2D.f9h"
  end subroutine Allocate_Test_Character_2d
  ! -------------------------------  Allocate_Test_Character_2d_a  -----
  subroutine Allocate_Test_Character_2d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    character(len=*), pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim(2)    ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(2) ! default 1
    character(len=*), intent(in), optional :: Fill ! To fill allocated array
    character(len=1), parameter :: Default = ''
    include "Allocate_Test_2D_a.f9h"
  end subroutine Allocate_Test_Character_2d_a
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
    include "Allocate_Test_3D.f9h"
  end subroutine Allocate_Test_Character_3d
  ! -------------------------------  Allocate_Test_Character_3d_a  -----
  subroutine Allocate_Test_Character_3d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    character(len=*), pointer, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim(3)    ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(3) ! default 1
    character(len=*), intent(in), optional :: Fill ! To fill allocated array
    character(len=1), parameter :: Default = ''
    include "Allocate_Test_3D_a.f9h"
  end subroutine Allocate_Test_Character_3d_a
  ! -----------------------------------  Allocate_Test_Complex_1d  -----
  subroutine Allocate_Test_Complex_1d ( To_Allocate, Dim1, &
    & ItsName, ModuleName, LowBound, Fill )
    complex, pointer, dimension(:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: LowBound     ! Lower bound, default 1
    complex, intent(in), optional :: Fill
    complex, parameter :: Default = (0.0,0.0)
    include "Allocate_Test_1D.f9h"
  end subroutine Allocate_Test_Complex_1d
  ! -----------------------------------  Allocate_Test_Complex_2d  -----
  subroutine Allocate_Test_Complex_2d ( To_Allocate, Dim1, Dim2, &
    & ItsName, ModuleName, Low1, Low2, Fill )
    complex, pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2 ! default 1
    complex, intent(in), optional :: Fill ! To fill allocated array
    complex, parameter :: Default = (0.0,0.0)
    include "Allocate_Test_2D.f9h"
  end subroutine Allocate_Test_Complex_2d
  ! ---------------------------------  Allocate_Test_Complex_2d_a  -----
  subroutine Allocate_Test_Complex_2d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    complex, pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim(2)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(2) ! default 1
    complex, intent(in), optional :: Fill ! To fill allocated array
    complex, parameter :: Default = (0.0,0.0)
    include "Allocate_Test_2D_a.f9h"
  end subroutine Allocate_Test_Complex_2d_a
  ! -----------------------------------  Allocate_Test_Complex_3d  -----
  subroutine Allocate_Test_Complex_3d ( To_Allocate, Dim1, Dim2, Dim3, &
    & ItsName, ModuleName, Low1, Low2, Low3, Fill )
    complex, pointer, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    integer, intent(in) :: Dim3    ! Upper bound of third dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2, Low3 ! Low bounds for dimensions
    complex, intent(in), optional :: Fill ! To fill allocated array
    complex, parameter :: Default = (0.0,0.0)
    include "Allocate_Test_3D.f9h"
  end subroutine Allocate_Test_Complex_3d
  ! ---------------------------------  Allocate_Test_Complex_3d_a  -----
  subroutine Allocate_Test_Complex_3d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    complex, pointer, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim(3)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(3) ! default 1
    complex, intent(in), optional :: Fill ! To fill allocated array
    complex, parameter :: Default = (0.0,0.0)
    include "Allocate_Test_3D_a.f9h"
  end subroutine Allocate_Test_Complex_3d_a
  ! ----------------------------------  Allocate_Test_DComplex_1d  -----
  subroutine Allocate_Test_DComplex_1d ( To_Allocate, Dim1, &
    & ItsName, ModuleName, LowBound, Fill )
    complex(kind(0.0d0)), pointer, dimension(:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: LowBound     ! Lower bound, default 1
    complex(kind(0.0d0)), intent(in), optional :: Fill
    complex(kind(0.0d0)), parameter :: Default = (0.0d0,0.0d0)
    include "Allocate_Test_1D.f9h"
  end subroutine Allocate_Test_DComplex_1d
  ! ----------------------------------  Allocate_Test_DComplex_2d  -----
  subroutine Allocate_Test_DComplex_2d ( To_Allocate, Dim1, Dim2, &
    & ItsName, ModuleName, Low1, Low2, Fill )
    complex(kind(0.0d0)), pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2 ! default 1
    complex(kind(0.0d0)), intent(in), optional :: Fill ! To fill allocated array
    complex(kind(0.0d0)), parameter :: Default = (0.0d0,0.0d0)
    include "Allocate_Test_2D.f9h"
  end subroutine Allocate_Test_DComplex_2d
  ! --------------------------------  Allocate_Test_DComplex_2d_a  -----
  subroutine Allocate_Test_DComplex_2d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    complex(kind(0.0d0)), pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim(2)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(2) ! default 1
    complex(kind(0.0d0)), intent(in), optional :: Fill ! To fill allocated array
    complex(kind(0.0d0)), parameter :: Default = (0.0d0,0.0d0)
    include "Allocate_Test_2D_a.f9h"
  end subroutine Allocate_Test_DComplex_2d_a
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
    include "Allocate_Test_3D.f9h"
  end subroutine Allocate_Test_DComplex_3d
  ! --------------------------------  Allocate_Test_DComplex_3d_a  -----
  subroutine Allocate_Test_DComplex_3d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    complex(kind(0.0d0)), pointer, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim(3)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(3) ! default 1
    complex(kind(0.0d0)), intent(in), optional :: Fill ! To fill allocated array
    complex(kind(0.0d0)), parameter :: Default = (0.0d0,0.0d0)
    include "Allocate_Test_3D_a.f9h"
  end subroutine Allocate_Test_DComplex_3d_a
  ! ------------------------------------  Allocate_Test_RealR8_1d  -----
  subroutine Allocate_Test_RealR8_1d ( To_Allocate, Dim1, &
    & ItsName, ModuleName, LowBound, Fill )
    double precision, pointer, dimension(:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: LowBound     ! Lower bound, default 1
    double precision, intent(in), optional :: Fill
    double precision, parameter :: Default = 0.0d0
    include "Allocate_Test_1D.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
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
    include "Allocate_Test_2D.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine Allocate_Test_RealR8_2d
  ! ----------------------------------  Allocate_Test_RealR8_2d_a  -----
  subroutine Allocate_Test_RealR8_2d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    double precision, pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim(2)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(2) ! Low bounds for dimensions
    double precision, intent(in), optional :: Fill
    double precision, parameter :: Default = 0.0d0
    include "Allocate_Test_2D_a.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine Allocate_Test_RealR8_2d_a
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
    include "Allocate_Test_3D.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine Allocate_Test_RealR8_3d
  ! ----------------------------------  Allocate_Test_RealR8_3d_a  -----
  subroutine Allocate_Test_RealR8_3d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    double precision, pointer, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim(3)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(3) ! Low bounds for dimensions
    double precision, intent(in), optional :: Fill
    double precision, parameter :: Default = 0.0d0
    include "Allocate_Test_3D_a.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine Allocate_Test_RealR8_3d_a
  ! ------------------------------------  Allocate_Test_RealR8_4d  -----
  subroutine Allocate_Test_RealR8_4d ( To_Allocate, Dim1, Dim2, Dim3, Dim4, &
    & ItsName, ModuleName, Low1, Low2, Low3, Low4, Fill )
    double precision, pointer, dimension(:,:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    integer, intent(in) :: Dim3    ! Upper bound of third dim. of To_Allocate
    integer, intent(in) :: Dim4    ! Upper bound of fourth dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2, Low3, Low4 ! Low bounds for dimensions
    double precision, intent(in), optional :: Fill
    double precision, parameter :: Default = 0.0d0
    include "Allocate_Test_4D.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine Allocate_Test_RealR8_4d
  ! ----------------------------------  Allocate_Test_RealR8_4d_a  -----
  subroutine Allocate_Test_RealR8_4d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    double precision, pointer, dimension(:,:,:,:) :: To_Allocate
    integer, intent(in) :: Dim(4)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(4) ! Low bounds for dimensions
    double precision, intent(in), optional :: Fill
    double precision, parameter :: Default = 0.0d0
    include "Allocate_Test_4D_a.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine Allocate_Test_RealR8_4d_a
  ! -----------------------------------  Allocate_Test_Integer_1d  -----
  subroutine Allocate_Test_Integer_1d ( To_Allocate, Dim1, &
    & ItsName, ModuleName, LowBound, Fill )
    integer, pointer, dimension(:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: LowBound     ! Lower bound, default 1
    integer, intent(in), optional :: Fill ! To fill allocated array
    integer, parameter :: Default = 0
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
    include "Allocate_Test_2D.f9h"
  end subroutine Allocate_Test_Integer_2d
  ! ---------------------------------  Allocate_Test_Integer_2d_a  -----
  subroutine Allocate_Test_Integer_2d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    integer, pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim(2)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(2) ! default 1
    integer, intent(in), optional :: Fill ! To fill allocated array
    integer, parameter :: Default = 0
    include "Allocate_Test_2D_a.f9h"
  end subroutine Allocate_Test_Integer_2d_a
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
    include "Allocate_Test_3D.f9h"
  end subroutine Allocate_Test_Integer_3d
  ! ---------------------------------  Allocate_Test_Integer_3d_a  -----
  subroutine Allocate_Test_Integer_3d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    integer, pointer, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim(3)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(3) ! default 1
    integer, intent(in), optional :: Fill ! To fill allocated array
    integer, parameter :: Default = 0
    include "Allocate_Test_3D_a.f9h"
  end subroutine Allocate_Test_Integer_3d_a
  ! -----------------------------------  Allocate_Test_Integer_4d  -----
  subroutine Allocate_Test_Integer_4d ( To_Allocate, Dim1, Dim2, Dim3, Dim4, &
    & ItsName, ModuleName, Low1, Low2, Low3, Low4, Fill )
    integer, pointer, dimension(:,:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    integer, intent(in) :: Dim3    ! Upper bound of third dim. of To_Allocate
    integer, intent(in) :: Dim4    ! Upper bound of fourth dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2, Low3, Low4 ! default 1
    integer, intent(in), optional :: Fill ! To fill allocated array
    integer, parameter :: Default = 0
    include "Allocate_Test_4D.f9h"
  end subroutine Allocate_Test_Integer_4d
  ! ---------------------------------  Allocate_Test_Integer_4d_a  -----
  subroutine Allocate_Test_Integer_4d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    integer, pointer, dimension(:,:,:,:) :: To_Allocate
    integer, intent(in) :: Dim(4)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(4) ! default 1
    integer, intent(in), optional :: Fill ! To fill allocated array
    integer, parameter :: Default = 0
    include "Allocate_Test_4D_a.f9h"
  end subroutine Allocate_Test_Integer_4d_a
  ! -----------------------------------  Allocate_Test_Logical_1d  -----
  subroutine Allocate_Test_Logical_1d ( To_Allocate, Dim1, &
    & ItsName, ModuleName, LowBound, Fill )
    logical, pointer, dimension(:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: LowBound     ! Lower bound, default 1
    logical, intent(in), optional :: Fill ! To fill allocated array
    logical, parameter :: Default = .false.
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
    include "Allocate_Test_2D.f9h"
  end subroutine Allocate_Test_Logical_2d
  ! ---------------------------------  Allocate_Test_Logical_2d_a  -----
  subroutine Allocate_Test_Logical_2d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    logical, pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim(2)  ! Upper bounds To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(2) ! Lower bounds, default 1
    logical, intent(in), optional :: Fill ! To fill allocated array
    logical, parameter :: Default = .false.
    include "Allocate_Test_2D_a.f9h"
  end subroutine Allocate_Test_Logical_2d_a
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
    include "Allocate_Test_3D.f9h"
  end subroutine Allocate_Test_Logical_3d
  ! ---------------------------------  Allocate_Test_Logical_3d_a  -----
  subroutine Allocate_Test_Logical_3d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    logical, pointer, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim(3)  ! Upper bounds To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(3) ! Lower bounds, default 1
    logical, intent(in), optional :: Fill ! To fill allocated array
    logical, parameter :: Default = .false.
    include "Allocate_Test_3D_a.f9h"
  end subroutine Allocate_Test_Logical_3d_a
  ! --------------------------------------  Allocate_Test_RealR4_1d  -----
  subroutine Allocate_Test_RealR4_1d ( To_Allocate, Dim1, &
    & ItsName, ModuleName, LowBound, Fill )
    real, pointer, dimension(:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: LowBound     ! Lower bound, default 1
    real, intent(in), optional :: Fill
    real, parameter :: Default = 0.0
    include "Allocate_Test_1D.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
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
    include "Allocate_Test_2D.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine Allocate_Test_RealR4_2d
  ! ------------------------------------  Allocate_Test_RealR4_2d_a  -----
  subroutine Allocate_Test_RealR4_2d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    real, pointer, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim(2)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(2) ! Low bounds, default 1
    real, intent(in), optional :: Fill
    real, parameter :: Default = 0.0
    include "Allocate_Test_2D_a.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine Allocate_Test_RealR4_2d_a
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
    include "Allocate_Test_3D.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine Allocate_Test_RealR4_3d
  ! ------------------------------------  Allocate_Test_RealR4_3d_a  -----
  subroutine Allocate_Test_RealR4_3d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    real, pointer, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim(3)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(3) ! Low bounds, default 1
    real, intent(in), optional :: Fill
    real, parameter :: Default = 0.0
    include "Allocate_Test_3D_a.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine Allocate_Test_RealR4_3d_a
  ! ------------------------------------  Allocate_Test_RealR4_4d  -----
  subroutine Allocate_Test_RealR4_4d ( To_Allocate, Dim1, Dim2, Dim3, Dim4, &
    & ItsName, ModuleName, Low1, Low2, Low3, Low4, Fill )
    real, pointer, dimension(:,:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    integer, intent(in) :: Dim3    ! Upper bound of third dim. of To_Allocate
    integer, intent(in) :: Dim4    ! Upper bound of fourth dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2, Low3, Low4 ! Low bounds for dimensions
    real, intent(in), optional :: Fill
    real, parameter :: Default = 0.0
    include "Allocate_Test_4D.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine Allocate_Test_RealR4_4d
  ! ------------------------------------  Allocate_Test_RealR4_4d_a  -----
  subroutine Allocate_Test_RealR4_4d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    real, pointer, dimension(:,:,:,:) :: To_Allocate
    integer, intent(in) :: Dim(4)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(4) ! Low bounds, default 1
    real, intent(in), optional :: Fill
    real, parameter :: Default = 0.0
    include "Allocate_Test_4D_a.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine Allocate_Test_RealR4_4d_a

  ! For allocatable instead of pointer To_Allocate
  ! ---------------------------------  Allocate_Test_Character_1d  -----
  subroutine AllocateA_Test_Character_1d ( To_Allocate, Dim1, &
    & ItsName, ModuleName, LowBound, Fill )
    character(len=*), allocatable, dimension(:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: LowBound     ! Lower bound, default 1
    character(len=*), intent(in), optional :: Fill
    character(len=1), parameter :: Default = ''
    include "AllocateA_Test_1D.f9h"
  end subroutine AllocateA_Test_Character_1d
  ! ---------------------------------  Allocate_Test_Character_2d  -----
  subroutine AllocateA_Test_Character_2d ( To_Allocate, Dim1, Dim2, &
    & ItsName, ModuleName, Low1, Low2, Fill )
    character(len=*), allocatable, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2 ! default 1
    character(len=*), intent(in), optional :: Fill ! To fill allocated array
    character(len=1), parameter :: Default = ''
    include "AllocateA_Test_2D.f9h"
  end subroutine AllocateA_Test_Character_2d
  ! -----------------------------------  Allocate_Test_Integer_1d  -----
  subroutine AllocateA_Test_Integer_1d ( To_Allocate, Dim1, &
    & ItsName, ModuleName, LowBound, Fill )
    integer, allocatable, dimension(:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: LowBound     ! Lower bound, default 1
    integer, intent(in), optional :: Fill
    integer, parameter :: Default = 0
    include "AllocateA_Test_1D.f9h"
  end subroutine AllocateA_Test_Integer_1d
  ! ----------------------------------  AllocateA_Test_Integer_2d  -----
  subroutine AllocateA_Test_Integer_2d ( To_Allocate, Dim1, Dim2, &
    & ItsName, ModuleName, Low1, Low2, Fill )
    integer, allocatable, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2 ! Low bounds for dimensions
    integer, intent(in), optional :: Fill
    integer, parameter :: Default = 0
    include "AllocateA_Test_2D.f9h"
  end subroutine AllocateA_Test_Integer_2d
  ! --------------------------------  AllocateA_Test_Integer_2d_a  -----
  subroutine AllocateA_Test_Integer_2d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    integer, allocatable, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim(2)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(2) ! Low bounds for dimensions
    integer, intent(in), optional :: Fill
    integer, parameter :: Default = 0
    include "AllocateA_Test_2D_a.f9h"
  end subroutine AllocateA_Test_Integer_2d_a
  ! ----------------------------------  AllocateA_Test_Integer_3d  -----
  subroutine AllocateA_Test_Integer_3d ( To_Allocate, Dim1, Dim2, Dim3, &
    & ItsName, ModuleName, Low1, Low2, Low3, Fill )
    integer, allocatable, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    integer, intent(in) :: Dim3    ! Upper bound of third dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2, Low3 ! Low bounds for dimensions
    integer, intent(in), optional :: Fill
    integer, parameter :: Default = 0
    include "AllocateA_Test_3D.f9h"
  end subroutine AllocateA_Test_Integer_3d
  ! --------------------------------  AllocateA_Test_Integer_3d_a  -----
  subroutine AllocateA_Test_Integer_3d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    integer, allocatable, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim(3)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(3) ! Low bounds for dimensions
    integer, intent(in), optional :: Fill
    integer, parameter :: Default = 0
    include "AllocateA_Test_3D_a.f9h"
  end subroutine AllocateA_Test_Integer_3d_a
  ! ----------------------------------  AllocateA_Test_Integer_4d  -----
  subroutine AllocateA_Test_Integer_4d ( To_Allocate, Dim1, Dim2, Dim3, Dim4, &
    & ItsName, ModuleName, Low1, Low2, Low3, Low4, Fill )
    integer, allocatable, dimension(:,:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    integer, intent(in) :: Dim3    ! Upper bound of third dim. of To_Allocate
    integer, intent(in) :: Dim4    ! Upper bound of fourth dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2, Low3, Low4 ! Low bounds for dimensions
    integer, intent(in), optional :: Fill
    integer, parameter :: Default = 0
    include "AllocateA_Test_4D.f9h"
  end subroutine AllocateA_Test_Integer_4d
  ! --------------------------------  AllocateA_Test_Integer_4d_a  -----
  subroutine AllocateA_Test_Integer_4d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    integer, allocatable, dimension(:,:,:,:) :: To_Allocate
    integer, intent(in) :: Dim(4)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(4) ! Low bounds for dimensions
    integer, intent(in), optional :: Fill
    integer, parameter :: Default = 0
    include "AllocateA_Test_4D_a.f9h"
  end subroutine AllocateA_Test_Integer_4d_a
  ! -----------------------------------  Allocate_Test_Integer_1d  -----
  subroutine AllocateA_Test_Logical_1d ( To_Allocate, Dim1, &
    & ItsName, ModuleName, LowBound, Fill )
    logical, allocatable, dimension(:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: LowBound     ! Lower bound, default 1
    logical, intent(in), optional :: Fill
    logical, parameter :: Default = .false.
    include "AllocateA_Test_1D.f9h"
  end subroutine AllocateA_Test_Logical_1d
  ! ----------------------------------  AllocateA_Test_Logical_2d  -----
  subroutine AllocateA_Test_Logical_2d ( To_Allocate, Dim1, Dim2, &
    & ItsName, ModuleName, Low1, Low2, Fill )
    logical, allocatable, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2 ! Low bounds for dimensions
    logical, intent(in), optional :: Fill
    logical, parameter :: Default = .false.
    include "AllocateA_Test_2D.f9h"
  end subroutine AllocateA_Test_Logical_2d
  ! --------------------------------  AllocateA_Test_Logical_2d_a  -----
  subroutine AllocateA_Test_Logical_2d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    logical, allocatable, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim(2)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(2) ! Low bounds for dimensions
    logical, intent(in), optional :: Fill
    logical, parameter :: Default = .false.
    include "AllocateA_Test_2D_a.f9h"
  end subroutine AllocateA_Test_Logical_2d_a
  ! ----------------------------------  AllocateA_Test_Logical_3d  -----
  subroutine AllocateA_Test_Logical_3d ( To_Allocate, Dim1, Dim2, Dim3, &
    & ItsName, ModuleName, Low1, Low2, Low3, Fill )
    logical, allocatable, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    integer, intent(in) :: Dim3    ! Upper bound of third dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2, Low3 ! Low bounds for dimensions
    logical, intent(in), optional :: Fill
    logical, parameter :: Default = .false.
    include "AllocateA_Test_3D.f9h"
  end subroutine AllocateA_Test_Logical_3d
  ! --------------------------------  AllocateA_Test_Logical_3d_a  -----
  subroutine AllocateA_Test_Logical_3d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    logical, allocatable, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim(3)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(3) ! Low bounds for dimensions
    logical, intent(in), optional :: Fill
    logical, parameter :: Default = .false.
    include "AllocateA_Test_3D_a.f9h"
  end subroutine AllocateA_Test_Logical_3d_a
  ! -----------------------------------  AllocateA_Test_RealR8_1d  -----
  subroutine AllocateA_Test_RealR8_1d ( To_Allocate, Dim1, &
    & ItsName, ModuleName, LowBound, Fill )
    double precision, allocatable, dimension(:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: LowBound     ! Lower bound, default 1
    double precision, intent(in), optional :: Fill
    double precision, parameter :: Default = 0.0d0
    include "AllocateA_Test_1D.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine AllocateA_Test_RealR8_1d
  ! -----------------------------------  AllocateA_Test_RealR8_2d  -----
  subroutine AllocateA_Test_RealR8_2d ( To_Allocate, Dim1, Dim2, &
    & ItsName, ModuleName, Low1, Low2, Fill )
    double precision, allocatable, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2 ! Low bounds for dimensions
    double precision, intent(in), optional :: Fill
    double precision, parameter :: Default = 0.0d0
    include "AllocateA_Test_2D.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine AllocateA_Test_RealR8_2d
  ! ----------------------------------  AllocateA_Test_RealR8_2d_a  -----
  subroutine AllocateA_Test_RealR8_2d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    double precision, allocatable, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim(2)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(2) ! Low bounds for dimensions
    double precision, intent(in), optional :: Fill
    double precision, parameter :: Default = 0.0d0
    include "AllocateA_Test_2D_a.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine AllocateA_Test_RealR8_2d_a
  ! ------------------------------------  AllocateA_Test_RealR8_3d  -----
  subroutine AllocateA_Test_RealR8_3d ( To_Allocate, Dim1, Dim2, Dim3, &
    & ItsName, ModuleName, Low1, Low2, Low3, Fill )
    double precision, allocatable, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    integer, intent(in) :: Dim3    ! Upper bound of third dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2, Low3 ! Low bounds for dimensions
    double precision, intent(in), optional :: Fill
    double precision, parameter :: Default = 0.0d0
    include "AllocateA_Test_3D.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine AllocateA_Test_RealR8_3d
  ! ----------------------------------  AllocateA_Test_RealR8_3d_a  -----
  subroutine AllocateA_Test_RealR8_3d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    double precision, allocatable, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim(3)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(3) ! Low bounds for dimensions
    double precision, intent(in), optional :: Fill
    double precision, parameter :: Default = 0.0d0
    include "AllocateA_Test_3D_a.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine AllocateA_Test_RealR8_3d_a
  ! ------------------------------------  AllocateA_Test_RealR8_4d  -----
  subroutine AllocateA_Test_RealR8_4d ( To_Allocate, Dim1, Dim2, Dim3, Dim4, &
    & ItsName, ModuleName, Low1, Low2, Low3, Low4, Fill )
    double precision, allocatable, dimension(:,:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    integer, intent(in) :: Dim3    ! Upper bound of third dim. of To_Allocate
    integer, intent(in) :: Dim4    ! Upper bound of fourth dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2, Low3, Low4 ! Low bounds for dimensions
    double precision, intent(in), optional :: Fill
    double precision, parameter :: Default = 0.0d0
    include "AllocateA_Test_4D.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine AllocateA_Test_RealR8_4d
  ! ----------------------------------  AllocateA_Test_RealR8_4d_a  -----
  subroutine AllocateA_Test_RealR8_4d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    double precision, allocatable, dimension(:,:,:,:) :: To_Allocate
    integer, intent(in) :: Dim(4)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(4) ! Low bounds for dimensions
    double precision, intent(in), optional :: Fill
    double precision, parameter :: Default = 0.0d0
    include "AllocateA_Test_4D_a.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine AllocateA_Test_RealR8_4d_a
  ! -------------------------------------  AllocateA_Test_RealR4_1d  -----
  subroutine AllocateA_Test_RealR4_1d ( To_Allocate, Dim1, &
    & ItsName, ModuleName, LowBound, Fill )
    real, allocatable, dimension(:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: LowBound     ! Lower bound, default 1
    real, intent(in), optional :: Fill
    real, parameter :: Default = 0.0
    include "AllocateA_Test_1D.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine AllocateA_Test_RealR4_1d
  ! --------------------------------------  AllocateA_Test_RealR4_2d  -----
  subroutine AllocateA_Test_RealR4_2d ( To_Allocate, Dim1, Dim2, &
    & ItsName, ModuleName, Low1, Low2, Fill )
    real, allocatable, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2 ! Low bounds for dimensions
    real, intent(in), optional :: Fill
    real, parameter :: Default = 0.0
    include "AllocateA_Test_2D.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine AllocateA_Test_RealR4_2d
  ! ------------------------------------  AllocateA_Test_RealR4_2d_a  -----
  subroutine AllocateA_Test_RealR4_2d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    real, allocatable, dimension(:,:) :: To_Allocate
    integer, intent(in) :: Dim(2)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(2) ! Low bounds, default 1
    real, intent(in), optional :: Fill
    real, parameter :: Default = 0.0
    include "AllocateA_Test_2D_a.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine AllocateA_Test_RealR4_2d_a
  ! ------------------------------------  AllocateA_Test_RealR4_3d  -----
  subroutine AllocateA_Test_RealR4_3d ( To_Allocate, Dim1, Dim2, Dim3, &
    & ItsName, ModuleName, Low1, Low2, Low3, Fill )
    real, allocatable, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    integer, intent(in) :: Dim3    ! Upper bound of third dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2, Low3 ! Low bounds for dimensions
    real, intent(in), optional :: Fill
    real, parameter :: Default = 0.0
    include "AllocateA_Test_3D.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine AllocateA_Test_RealR4_3d
  ! ------------------------------------  AllocateA_Test_RealR4_3d_a  -----
  subroutine AllocateA_Test_RealR4_3d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    real, allocatable, dimension(:,:,:) :: To_Allocate
    integer, intent(in) :: Dim(3)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(3) ! Low bounds, default 1
    real, intent(in), optional :: Fill
    real, parameter :: Default = 0.0
    include "AllocateA_Test_3D_a.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine AllocateA_Test_RealR4_3d_a
  ! ------------------------------------  AllocateA_Test_RealR4_4d  -----
  subroutine AllocateA_Test_RealR4_4d ( To_Allocate, Dim1, Dim2, Dim3, Dim4, &
    & ItsName, ModuleName, Low1, Low2, Low3, Low4, Fill )
    real, allocatable, dimension(:,:,:,:) :: To_Allocate
    integer, intent(in) :: Dim1    ! Upper bound of first dim. of To_Allocate
    integer, intent(in) :: Dim2    ! Upper bound of second dim. of To_Allocate
    integer, intent(in) :: Dim3    ! Upper bound of third dim. of To_Allocate
    integer, intent(in) :: Dim4    ! Upper bound of fourth dim. of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low1, Low2, Low3, Low4 ! Low bounds for dimensions
    real, intent(in), optional :: Fill
    real, parameter :: Default = 0.0
    include "AllocateA_Test_4D.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine AllocateA_Test_RealR4_4d
  ! ------------------------------------  AllocateA_Test_RealR4_4d_a  -----
  subroutine AllocateA_Test_RealR4_4d_a ( To_Allocate, Dim, &
    & ItsName, ModuleName, Low, Fill )
    real, allocatable, dimension(:,:,:,:) :: To_Allocate
    integer, intent(in) :: Dim(4)  ! Upper bounds of To_Allocate
    character(len=*), intent(in) :: ItsName, ModuleName
    integer, intent(in), optional :: Low(4) ! Low bounds, default 1
    real, intent(in), optional :: Fill
    real, parameter :: Default = 0.0
    include "AllocateA_Test_4D_a.f9h"
    if ( InitRealNaN .and. IEEE_Support_DataType(to_allocate) .and. &
       & .not. present(fill) ) &
         & to_allocate = IEEE_Value ( to_allocate, IEEE_Signaling_NaN )
  end subroutine AllocateA_Test_RealR4_4d_a

  ! -------------------------------  Deallocate_Test_Character_1d  -----
  subroutine Deallocate_Test_Character_1d ( To_Deallocate, ItsName, ModuleName )
    character(len=*), pointer, dimension(:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Character_1d
  ! -------------------------------  Deallocate_Test_Character_2d  -----
  subroutine Deallocate_Test_Character_2d ( To_Deallocate, ItsName, ModuleName )
    character(len=*), pointer, dimension(:,:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Character_2d
  ! -------------------------------  Deallocate_Test_Character_3d  -----
  subroutine Deallocate_Test_Character_3d ( To_Deallocate, ItsName, ModuleName )
    character(len=*), pointer, dimension(:,:,:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Character_3d
  ! ---------------------------------  Deallocate_Test_Complex_1d  -----
  subroutine Deallocate_Test_Complex_1d ( To_Deallocate, ItsName, ModuleName )
    complex, pointer, dimension(:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Complex_1d
  ! ---------------------------------  Deallocate_Test_Complex_2d  -----
  subroutine Deallocate_Test_Complex_2d ( To_Deallocate, ItsName, ModuleName )
    complex, pointer, dimension(:,:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Complex_2d
  ! ---------------------------------  Deallocate_Test_Complex_3d  -----
  subroutine Deallocate_Test_Complex_3d ( To_Deallocate, ItsName, ModuleName )
    complex, pointer, dimension(:,:,:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Complex_3d
  ! --------------------------------  Deallocate_Test_Dcomplex_1d  -----
  subroutine Deallocate_Test_Dcomplex_1d ( To_Deallocate, ItsName, ModuleName )
    complex(kind(0.0d0)), pointer, dimension(:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Dcomplex_1d
  ! --------------------------------  Deallocate_Test_Dcomplex_2d  -----
  subroutine Deallocate_Test_Dcomplex_2d ( To_Deallocate, ItsName, ModuleName )
    complex(kind(0.0d0)), pointer, dimension(:,:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Dcomplex_2d
  ! --------------------------------  Deallocate_Test_Dcomplex_3d  -----
  subroutine Deallocate_Test_Dcomplex_3d ( To_Deallocate, ItsName, ModuleName )
    complex(kind(0.0d0)), pointer, dimension(:,:,:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Dcomplex_3d
  ! ----------------------------------  Deallocate_Test_RealR8_1d  -----
  subroutine Deallocate_Test_RealR8_1d ( To_Deallocate, ItsName, ModuleName )
    double precision, pointer, dimension(:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_RealR8_1d
  ! ----------------------------------  Deallocate_Test_RealR8_2d  -----
  subroutine Deallocate_Test_RealR8_2d ( To_Deallocate, ItsName, ModuleName )
    double precision, pointer, dimension(:,:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_RealR8_2d
  ! ----------------------------------  Deallocate_Test_RealR8_3d  -----
  subroutine Deallocate_Test_RealR8_3d ( To_Deallocate, ItsName, ModuleName )
    double precision, pointer, dimension(:,:,:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_RealR8_3d
  ! ----------------------------------  Deallocate_Test_RealR8_4d  -----
  subroutine Deallocate_Test_RealR8_4d ( To_Deallocate, ItsName, ModuleName )
    double precision, pointer, dimension(:,:,:,:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_RealR8_4d
  ! ---------------------------------  Deallocate_Test_Integer_1d  -----
  subroutine Deallocate_Test_Integer_1d ( To_Deallocate, ItsName, ModuleName )
    integer, pointer, dimension(:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Integer_1d
  ! ---------------------------------  Deallocate_Test_Integer_2d  -----
  subroutine Deallocate_Test_Integer_2d ( To_Deallocate, ItsName, ModuleName )
    integer, pointer, dimension(:,:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Integer_2d
  ! ---------------------------------  Deallocate_Test_Integer_3d  -----
  subroutine Deallocate_Test_Integer_3d ( To_Deallocate, ItsName, ModuleName )
    integer, pointer, dimension(:,:,:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Integer_3d
  ! ---------------------------------  Deallocate_Test_Integer_4d  -----
  subroutine Deallocate_Test_Integer_4d ( To_Deallocate, ItsName, ModuleName )
    integer, pointer, dimension(:,:,:,:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Integer_4d
  ! ---------------------------------  Deallocate_Test_Logical_1d  -----
  subroutine Deallocate_Test_Logical_1d ( To_Deallocate, ItsName, ModuleName )
    logical, pointer, dimension(:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Logical_1d
  ! ---------------------------------  Deallocate_Test_Logical_2d  -----
  subroutine Deallocate_Test_Logical_2d ( To_Deallocate, ItsName, ModuleName )
    logical, pointer, dimension(:,:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Logical_2d
  ! ---------------------------------  Deallocate_Test_Logical_3d  -----
  subroutine Deallocate_Test_Logical_3d ( To_Deallocate, ItsName, ModuleName )
    logical, pointer, dimension(:,:,:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_Logical_3d
  ! ------------------------------------  Deallocate_Test_RealR4_1d  -----
  subroutine Deallocate_Test_RealR4_1d ( To_Deallocate, ItsName, ModuleName )
    real, pointer, dimension(:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_RealR4_1d
  ! ------------------------------------  Deallocate_Test_RealR4_2d  -----
  subroutine Deallocate_Test_RealR4_2d ( To_Deallocate, ItsName, ModuleName )
    real, pointer, dimension(:,:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_RealR4_2d
  ! ----------------------------------  Deallocate_Test_RealR4_3d  -----
  subroutine Deallocate_Test_RealR4_3d ( To_Deallocate, ItsName, ModuleName )
    real, pointer, dimension(:,:,:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_RealR4_3d
  ! ----------------------------------  Deallocate_Test_RealR4_4d  -----
  subroutine Deallocate_Test_RealR4_4d ( To_Deallocate, ItsName, ModuleName )
    real, pointer, dimension(:,:,:,:) :: To_Deallocate
    include "Deallocate_Test.f9h"
  end subroutine Deallocate_Test_RealR4_4d

! For allocatable instead of pointer arguments

  ! -------------------------------  Deallocate_Test_Character_1d  -----
  subroutine DeallocateA_Test_Character_1d ( To_Deallocate, ItsName, ModuleName )
    character(len=*), allocatable, dimension(:) :: To_Deallocate
    include "DeallocateA_Test.f9h"
  end subroutine DeallocateA_Test_Character_1d
  ! -------------------------------  Deallocate_Test_Character_2d  -----
  subroutine DeallocateA_Test_Character_2d ( To_Deallocate, ItsName, ModuleName )
    character(len=*), allocatable, dimension(:,:) :: To_Deallocate
    include "DeallocateA_Test.f9h"
  end subroutine DeallocateA_Test_Character_2d
  ! --------------------------------  DeallocateA_Test_Integer_1d  -----
  subroutine DeallocateA_Test_Integer_1d ( To_Deallocate, ItsName, ModuleName )
    integer, allocatable, dimension(:) :: To_Deallocate
    include "DeallocateA_Test.f9h"
  end subroutine DeallocateA_Test_Integer_1d
  ! --------------------------------  DeallocateA_Test_Integer_2d  -----
  subroutine DeallocateA_Test_Integer_2d ( To_Deallocate, ItsName, ModuleName )
    integer, allocatable, dimension(:,:) :: To_Deallocate
    include "DeallocateA_Test.f9h"
  end subroutine DeallocateA_Test_Integer_2d
  ! --------------------------------  DeallocateA_Test_Integer_3d  -----
  subroutine DeallocateA_Test_Integer_3d ( To_Deallocate, ItsName, ModuleName )
    integer, allocatable, dimension(:,:,:) :: To_Deallocate
    include "DeallocateA_Test.f9h"
  end subroutine DeallocateA_Test_Integer_3d
  ! --------------------------------  DeallocateA_Test_Integer_4d  -----
  subroutine DeallocateA_Test_Integer_4d ( To_Deallocate, ItsName, ModuleName )
    integer, allocatable, dimension(:,:,:,:) :: To_Deallocate
    include "DeallocateA_Test.f9h"
  end subroutine DeallocateA_Test_Integer_4d
  ! --------------------------------  DeallocateA_Test_Logical_1d  -----
  subroutine DeallocateA_Test_Logical_1d ( To_Deallocate, ItsName, ModuleName )
    logical, allocatable, dimension(:) :: To_Deallocate
    include "DeallocateA_Test.f9h"
  end subroutine DeallocateA_Test_Logical_1d
  ! --------------------------------  DeallocateA_Test_Logical_2d  -----
  subroutine DeallocateA_Test_Logical_2d ( To_Deallocate, ItsName, ModuleName )
    logical, allocatable, dimension(:,:) :: To_Deallocate
    include "DeallocateA_Test.f9h"
  end subroutine DeallocateA_Test_Logical_2d
  ! --------------------------------  DeallocateA_Test_Logical_3d  -----
  subroutine DeallocateA_Test_Logical_3d ( To_Deallocate, ItsName, ModuleName )
    logical, allocatable, dimension(:,:,:) :: To_Deallocate
    include "DeallocateA_Test.f9h"
  end subroutine DeallocateA_Test_Logical_3d
  ! ----------------------------------  DeallocateA_Test_RealR8_1d  -----
  subroutine DeallocateA_Test_RealR8_1d ( To_Deallocate, ItsName, ModuleName )
    double precision, allocatable, dimension(:) :: To_Deallocate
    include "DeallocateA_Test.f9h"
  end subroutine DeallocateA_Test_RealR8_1d
  ! ----------------------------------  DeallocateA_Test_RealR8_2d  -----
  subroutine DeallocateA_Test_RealR8_2d ( To_Deallocate, ItsName, ModuleName )
    double precision, allocatable, dimension(:,:) :: To_Deallocate
    include "DeallocateA_Test.f9h"
  end subroutine DeallocateA_Test_RealR8_2d
  ! ----------------------------------  DeallocateA_Test_RealR8_3d  -----
  subroutine DeallocateA_Test_RealR8_3d ( To_Deallocate, ItsName, ModuleName )
    double precision, allocatable, dimension(:,:,:) :: To_Deallocate
    include "DeallocateA_Test.f9h"
  end subroutine DeallocateA_Test_RealR8_3d
  ! ----------------------------------  DeallocateA_Test_RealR8_4d  -----
  subroutine DeallocateA_Test_RealR8_4d ( To_Deallocate, ItsName, ModuleName )
    double precision, allocatable, dimension(:,:,:,:) :: To_Deallocate
    include "DeallocateA_Test.f9h"
  end subroutine DeallocateA_Test_RealR8_4d
  ! ------------------------------------  DeallocateA_Test_RealR4_1d  -----
  subroutine DeallocateA_Test_RealR4_1d ( To_Deallocate, ItsName, ModuleName )
    real, allocatable, dimension(:) :: To_Deallocate
    include "DeallocateA_Test.f9h"
  end subroutine DeallocateA_Test_RealR4_1d
  ! ------------------------------------  DeallocateA_Test_RealR4_2d  -----
  subroutine DeallocateA_Test_RealR4_2d ( To_Deallocate, ItsName, ModuleName )
    real, allocatable, dimension(:,:) :: To_Deallocate
    include "DeallocateA_Test.f9h"
  end subroutine DeallocateA_Test_RealR4_2d
  ! ----------------------------------  DeallocateA_Test_RealR4_3d  -----
  subroutine DeallocateA_Test_RealR4_3d ( To_Deallocate, ItsName, ModuleName )
    real, allocatable, dimension(:,:,:) :: To_Deallocate
    include "DeallocateA_Test.f9h"
  end subroutine DeallocateA_Test_RealR4_3d
  ! ----------------------------------  DeallocateA_Test_RealR4_4d  -----
  subroutine DeallocateA_Test_RealR4_4d ( To_Deallocate, ItsName, ModuleName )
    real, allocatable, dimension(:,:,:,:) :: To_Deallocate
    include "DeallocateA_Test.f9h"
  end subroutine DeallocateA_Test_RealR4_4d

  ! --------------------------------------------  Bytes functions  -----
  integer function Byte_Size_Character_1D ( A ) result ( Byte_Size )
    character(len=*), intent(in) :: A(:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_Character_1D
  integer function Byte_Size_Character_2D ( A ) result ( Byte_Size )
    character(len=*), intent(in) :: A(:,:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_Character_2D
  integer function Byte_Size_Character_3D ( A ) result ( Byte_Size )
    character(*), intent(in) :: A(:,:,:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_Character_3D
  integer function Byte_Size_Complex_1D ( A ) result ( Byte_Size )
    complex, intent(in) :: A(:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_Complex_1D
  integer function Byte_Size_Complex_2D ( A ) result ( Byte_Size )
    complex, intent(in) :: A(:,:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_Complex_2D
  integer function Byte_Size_Complex_3D ( A ) result ( Byte_Size )
    complex, intent(in) :: A(:,:,:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_Complex_3D
  integer function Byte_Size_DComplex_1D ( A ) result ( Byte_Size )
    complex(kind(0.0d0)), intent(in) :: A(:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_DComplex_1D
  integer function Byte_Size_DComplex_2D ( A ) result ( Byte_Size )
    complex(kind(0.0d0)), intent(in) :: A(:,:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_DComplex_2D
  integer function Byte_Size_DComplex_3D ( A ) result ( Byte_Size )
    complex(kind(0.0d0)), intent(in) :: A(:,:,:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_DComplex_3D
  integer function Byte_Size_Integer_1D ( A ) result ( Byte_Size )
    integer, intent(in) :: A(:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_Integer_1D
  integer function Byte_Size_Integer_2D ( A ) result ( Byte_Size )
    integer, intent(in) :: A(:,:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_Integer_2D
  integer function Byte_Size_Integer_3D ( A ) result ( Byte_Size )
    integer, intent(in) :: A(:,:,:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_Integer_3D
  integer function Byte_Size_Integer_4D ( A ) result ( Byte_Size )
    integer, intent(in) :: A(:,:,:,:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_Integer_4D
  integer function Byte_Size_Logical_1D ( A ) result ( Byte_Size )
    logical, intent(in) :: A(:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_Logical_1D
  integer function Byte_Size_Logical_2D ( A ) result ( Byte_Size )
    logical, intent(in) :: A(:,:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_Logical_2D
  integer function Byte_Size_Logical_3D ( A ) result ( Byte_Size )
    logical, intent(in) :: A(:,:,:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_Logical_3D
  integer function Byte_Size_RealR4_1D ( A ) result ( Byte_Size )
    real, intent(in) :: A(:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_RealR4_1D
  integer function Byte_Size_RealR4_2D ( A ) result ( Byte_Size )
    real, intent(in) :: A(:,:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_RealR4_2D
  integer function Byte_Size_RealR4_3D ( A ) result ( Byte_Size )
    real, intent(in) :: A(:,:,:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_RealR4_3D
  integer function Byte_Size_RealR4_4D ( A ) result ( Byte_Size )
    real, intent(in) :: A(:,:,:,:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_RealR4_4D
  integer function Byte_Size_RealR4_5D ( A ) result ( Byte_Size )
    real, intent(in) :: A(:,:,:,:,:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_RealR4_5D
  integer function Byte_Size_RealR4_6D ( A ) result ( Byte_Size )
    real, intent(in) :: A(:,:,:,:,:,:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_RealR4_6D
  integer function Byte_Size_RealR8_1D ( A ) result ( Byte_Size )
    double precision, intent(in) :: A(:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_RealR8_1D
  integer function Byte_Size_RealR8_2D ( A ) result ( Byte_Size )
    double precision, intent(in) :: A(:,:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_RealR8_2D
  integer function Byte_Size_RealR8_3D ( A ) result ( Byte_Size )
    double precision, intent(in) :: A(:,:,:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_RealR8_3D
  integer function Byte_Size_RealR8_4D ( A ) result ( Byte_Size )
    double precision, intent(in) :: A(:,:,:,:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_RealR8_4D
  integer function Byte_Size_RealR8_5D ( A ) result ( Byte_Size )
    double precision, intent(in) :: A(:,:,:,:,:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_RealR8_5D
  integer function Byte_Size_RealR8_6D ( A ) result ( Byte_Size )
    double precision, intent(in) :: A(:,:,:,:,:,:)
    byte_size = bytes(a) * size(a)
  end function Byte_Size_RealR8_6D

  !??? We'd prefer to use storage_size(a) but Intel 13.0.1 complains ???
  integer function Bytes_Character_1D ( A ) result ( BYTES )
    character(len=*), intent(in) :: A(:)
    character, parameter :: B = '' ! Intel 13.0.1 doesn't like storage_size(a)
    ! When Intel is happy with storage_size(a), DO NOT multiply by len(a)!
    ! Storage_Size includes character length.
    bytes = ( ( storage_size(b) + 7 ) / 8 ) * len(a)
!     bytes = ( ( storage_size(a) + 7 ) / 8 )
  end function Bytes_Character_1D
  integer function Bytes_Character_2D ( A ) result ( BYTES )
    character(len=*), intent(in) :: A(:,:)
    character, parameter :: B = '' ! Intel 13.0.1 doesn't like storage_size(a)
    ! When Intel is happy with storage_size(a), DO NOT multiply by len(a)!
    ! Storage_Size includes character length.
    bytes = ( ( storage_size(b) + 7 ) / 8 ) * len(a)
!     bytes = ( ( storage_size(a) + 7 ) / 8 )
  end function Bytes_Character_2D
  integer function Bytes_Character_3D ( A ) result ( BYTES )
    character(len=*), intent(in) :: A(:,:,:)
    character, parameter :: B = '' ! Intel 13.0.1 doesn't like storage_size(a)
    ! When Intel is happy with storage_size(a), DO NOT multiply by len(a)!
    ! Storage_Size includes character length.
    bytes = ( ( storage_size(b) + 7 ) / 8 ) * len(a)
!     bytes = ( ( storage_size(a) + 7 ) / 8 )
  end function Bytes_Character_3D
  integer function Bytes_Complex_1D ( A ) result ( BYTES )
    complex, intent(in) :: A(:)
    complex, parameter :: B = (0.0,0.0) ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_Complex_1D
  integer function Bytes_Complex_2D ( A ) result ( BYTES )
    complex, intent(in) :: A(:,:)
    complex, parameter :: B = (0.0,0.0) ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_Complex_2D
  integer function Bytes_Complex_3D ( A ) result ( BYTES )
    complex, intent(in) :: A(:,:,:)
    complex, parameter :: B = (0.0,0.0) ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_Complex_3D
  integer function Bytes_DComplex_1D ( A ) result ( BYTES )
    complex(kind(0.0d0)), intent(in) :: A(:)
    complex(kind(0.0d0)), parameter :: B = (0.0d0,0.0d0) ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_DComplex_1D
  integer function Bytes_DComplex_2D ( A ) result ( BYTES )
    complex(kind(0.0d0)), intent(in) :: A(:,:)
    complex(kind(0.0d0)), parameter :: B = (0.0d0,0.0d0) ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_DComplex_2D
  integer function Bytes_DComplex_3D ( A ) result ( BYTES )
    complex(kind(0.0d0)), intent(in) :: A(:,:,:)
    complex(kind(0.0d0)), parameter :: B = (0.0d0,0.0d0) ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_DComplex_3D
  integer function Bytes_Integer_1D ( A ) result ( BYTES )
    integer, intent(in) :: A(:)
    integer, parameter :: B = 0 ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_Integer_1D
  integer function Bytes_Integer_2D ( A ) result ( BYTES )
    integer, intent(in) :: A(:,:)
    integer, parameter :: B = 0 ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_Integer_2D
  integer function Bytes_Integer_3D ( A ) result ( BYTES )
    integer, intent(in) :: A(:,:,:)
    integer, parameter :: B = 0 ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_Integer_3D
  integer function Bytes_Integer_4D ( A ) result ( BYTES )
    integer, intent(in) :: A(:,:,:,:)
    integer, parameter :: B = 0 ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_Integer_4D
  integer function Bytes_Logical_1D ( A ) result ( BYTES )
    logical, intent(in) :: A(:)
    logical, parameter :: B = .false. ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_Logical_1D
  integer function Bytes_Logical_2D ( A ) result ( BYTES )
    logical, intent(in) :: A(:,:)
    logical, parameter :: B = .false. ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_Logical_2D
  integer function Bytes_Logical_3D ( A ) result ( BYTES )
    logical, intent(in) :: A(:,:,:)
    logical, parameter :: B = .false. ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_Logical_3D
  integer function Bytes_RealR4_1D ( A ) result ( BYTES )
    real, intent(in) :: A(:)
    real, parameter :: B = 0.0 ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_RealR4_1D
  integer function Bytes_RealR4_2D ( A ) result ( BYTES )
    real, intent(in) :: A(:,:)
    real, parameter :: B = 0.0 ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_RealR4_2D
  integer function Bytes_RealR4_3D ( A ) result ( BYTES )
    real, intent(in) :: A(:,:,:)
    real, parameter :: B = 0.0 ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_RealR4_3D
  integer function Bytes_RealR4_4D ( A ) result ( BYTES )
    real, intent(in) :: A(:,:,:,:)
    real, parameter :: B = 0.0 ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_RealR4_4D
  integer function Bytes_RealR4_5D ( A ) result ( BYTES )
    real, intent(in) :: A(:,:,:,:,:)
    real, parameter :: B = 0.0 ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_RealR4_5D
  integer function Bytes_RealR4_6D ( A ) result ( BYTES )
    real, intent(in) :: A(:,:,:,:,:,:)
    real, parameter :: B = 0.0 ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_RealR4_6D
  integer function Bytes_RealR8_1D ( A ) result ( BYTES )
    double precision, intent(in) :: A(:)
    double precision, parameter :: B = 0.0d0 ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_RealR8_1D
  integer function Bytes_RealR8_2D ( A ) result ( BYTES )
    double precision, intent(in) :: A(:,:)
    double precision, parameter :: B = 0.0d0 ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_RealR8_2D
  integer function Bytes_RealR8_3D ( A ) result ( BYTES )
    double precision, intent(in) :: A(:,:,:)
    double precision, parameter :: B = 0.0d0 ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_RealR8_3D
  integer function Bytes_RealR8_4D ( A ) result ( BYTES )
    double precision, intent(in) :: A(:,:,:,:)
    double precision, parameter :: B = 0.0d0 ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_RealR8_4D
  integer function Bytes_RealR8_5D ( A ) result ( BYTES )
    double precision, intent(in) :: A(:,:,:,:,:)
    double precision, parameter :: B = 0.0d0 ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_RealR8_5D
  integer function Bytes_RealR8_6D ( A ) result ( BYTES )
    double precision, intent(in) :: A(:,:,:,:,:,:)
    double precision, parameter :: B = 0.0d0 ! Intel 13.0.1 doesn't like storage_size(a)
    bytes = ( ( storage_size(b) + 7 ) / 8 )
  end function Bytes_RealR8_6D

  ! --------------------------------------------------  myMessage  -----
  subroutine myMessage ( severity, name, line )
    use MLSMessageModule, only: MLSMessage
    ! Args
    integer, intent(in)           :: severity
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: line
    call MLSMessage ( severity, ModuleName // '%' &
      & // trim(name), trim(line) )
  end subroutine myMessage

  subroutine myMessage_old ( severity, name, line )
    ! Args
    integer, intent(in)           :: severity
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: line

    ! Local variables
    integer :: nChars
    character(len=len(line) + len(name) + 3) :: thus
    ! Executable
    nChars = len(line)
    thus = line
    if ( len_trim(name) > 0 ) then
      nChars = len(line) + len_trim(snipRCSfrom(name)) + 3
      thus = '(' // trim(snipRCSfrom(name)) // ') ' // line
    end if
    if ( severity > MLSMSG_Warning ) then
      call PrintItOut( thus(1:nChars), SEVERITY, exitStatus = 1  )
    else
      call PrintItOut( thus(1:nChars), SEVERITY  )
    end if
  end subroutine myMessage_old

  ! ------------------------------------  Same_Shape_Character_1d  -----
  subroutine Same_Shape_Character_1d ( Ref, New, ItsName, ModuleName )
    character(len=*), pointer, dimension(:) :: Ref, New
    include "Same_Shape_1d.f9h"
  end subroutine Same_Shape_Character_1d
  ! ------------------------------------  Same_Shape_Character_2d  -----
  subroutine Same_Shape_Character_2d ( Ref, New, ItsName, ModuleName )
    character(len=*), pointer, dimension(:,:) :: Ref, New
    include "Same_Shape_md.f9h"
  end subroutine Same_Shape_Character_2d
  ! ------------------------------------  Same_Shape_Character_3d  -----
  subroutine Same_Shape_Character_3d ( Ref, New, ItsName, ModuleName )
    character(len=*), pointer, dimension(:,:,:) :: Ref, New
    include "Same_Shape_md.f9h"
  end subroutine Same_Shape_Character_3d
  ! --------------------------------------  Same_Shape_Complex_1d  -----
  subroutine Same_Shape_Complex_1d ( Ref, New, ItsName, ModuleName )
    complex, pointer, dimension(:) :: Ref, New
    include "Same_Shape_1d.f9h"
  end subroutine Same_Shape_Complex_1d
  ! --------------------------------------  Same_Shape_Complex_2d  -----
  subroutine Same_Shape_Complex_2d ( Ref, New, ItsName, ModuleName )
    complex, pointer, dimension(:,:) :: Ref, New
    include "Same_Shape_md.f9h"
  end subroutine Same_Shape_Complex_2d
  ! --------------------------------------  Same_Shape_Complex_3d  -----
  subroutine Same_Shape_Complex_3d ( Ref, New, ItsName, ModuleName )
    complex, pointer, dimension(:,:,:) :: Ref, New
    include "Same_Shape_md.f9h"
  end subroutine Same_Shape_Complex_3d
  ! -------------------------------------  Same_Shape_Dcomplex_1d  -----
  subroutine Same_Shape_Dcomplex_1d ( Ref, New, ItsName, ModuleName )
    complex(kind(0.0d0)), pointer, dimension(:) :: Ref, New
    include "Same_Shape_1d.f9h"
  end subroutine Same_Shape_Dcomplex_1d
  ! -------------------------------------  Same_Shape_Dcomplex_2d  -----
  subroutine Same_Shape_Dcomplex_2d ( Ref, New, ItsName, ModuleName )
    complex(kind(0.0d0)), pointer, dimension(:,:) :: Ref, New
    include "Same_Shape_md.f9h"
  end subroutine Same_Shape_Dcomplex_2d
  ! -------------------------------------  Same_Shape_Dcomplex_3d  -----
  subroutine Same_Shape_Dcomplex_3d ( Ref, New, ItsName, ModuleName )
    complex(kind(0.0d0)), pointer, dimension(:,:,:) :: Ref, New
    include "Same_Shape_md.f9h"
  end subroutine Same_Shape_Dcomplex_3d
  ! ---------------------------------------  Same_Shape_RealR8_1d  -----
  subroutine Same_Shape_RealR8_1d ( Ref, New, ItsName, ModuleName )
    double precision, pointer, dimension(:) :: Ref, New
    include "Same_Shape_1d.f9h"
  end subroutine Same_Shape_RealR8_1d
  ! ---------------------------------------  Same_Shape_RealR8_2d  -----
  subroutine Same_Shape_RealR8_2d ( Ref, New, ItsName, ModuleName )
    double precision, pointer, dimension(:,:) :: Ref, New
    include "Same_Shape_md.f9h"
  end subroutine Same_Shape_RealR8_2d
  ! ---------------------------------------  Same_Shape_RealR8_3d  -----
  subroutine Same_Shape_RealR8_3d ( Ref, New, ItsName, ModuleName )
    double precision, pointer, dimension(:,:,:) :: Ref, New
    include "Same_Shape_md.f9h"
  end subroutine Same_Shape_RealR8_3d
  ! ---------------------------------------  Same_Shape_RealR8_4d  -----
  subroutine Same_Shape_RealR8_4d ( Ref, New, ItsName, ModuleName )
    double precision, pointer, dimension(:,:,:,:) :: Ref, New
    include "Same_Shape_md.f9h"
  end subroutine Same_Shape_RealR8_4d
  ! --------------------------------------  Same_Shape_Integer_1d  -----
  subroutine Same_Shape_Integer_1d ( Ref, New, ItsName, ModuleName )
    integer, pointer, dimension(:) :: Ref, New
    include "Same_Shape_1d.f9h"
  end subroutine Same_Shape_Integer_1d
  ! --------------------------------------  Same_Shape_Integer_2d  -----
  subroutine Same_Shape_Integer_2d ( Ref, New, ItsName, ModuleName )
    integer, pointer, dimension(:,:) :: Ref, New
    include "Same_Shape_md.f9h"
  end subroutine Same_Shape_Integer_2d
  ! --------------------------------------  Same_Shape_Integer_3d  -----
  subroutine Same_Shape_Integer_3d ( Ref, New, ItsName, ModuleName )
    integer, pointer, dimension(:,:,:) :: Ref, New
    include "Same_Shape_md.f9h"
  end subroutine Same_Shape_Integer_3d
  ! --------------------------------------  Same_Shape_Integer_4d  -----
  subroutine Same_Shape_Integer_4d ( Ref, New, ItsName, ModuleName )
    integer, pointer, dimension(:,:,:,:) :: Ref, New
    include "Same_Shape_md.f9h"
  end subroutine Same_Shape_Integer_4d
  ! --------------------------------------  Same_Shape_Logical_1d  -----
  subroutine Same_Shape_Logical_1d ( Ref, New, ItsName, ModuleName )
    logical, pointer, dimension(:) :: Ref, New
    include "Same_Shape_1d.f9h"
  end subroutine Same_Shape_Logical_1d
  ! --------------------------------------  Same_Shape_Logical_2d  -----
  subroutine Same_Shape_Logical_2d ( Ref, New, ItsName, ModuleName )
    logical, pointer, dimension(:,:) :: Ref, New
    include "Same_Shape_md.f9h"
  end subroutine Same_Shape_Logical_2d
  ! --------------------------------------  Same_Shape_Logical_3d  -----
  subroutine Same_Shape_Logical_3d ( Ref, New, ItsName, ModuleName )
    logical, pointer, dimension(:,:,:) :: Ref, New
    include "Same_Shape_md.f9h"
  end subroutine Same_Shape_Logical_3d
  ! -----------------------------------------  Same_Shape_RealR4_1d  -----
  subroutine Same_Shape_RealR4_1d ( Ref, New, ItsName, ModuleName )
    real, pointer, dimension(:) :: Ref, New
    include "Same_Shape_1d.f9h"
  end subroutine Same_Shape_RealR4_1d
  ! -----------------------------------------  Same_Shape_RealR4_2d  -----
  subroutine Same_Shape_RealR4_2d ( Ref, New, ItsName, ModuleName )
    real, pointer, dimension(:,:) :: Ref, New
    include "Same_Shape_md.f9h"
  end subroutine Same_Shape_RealR4_2d
  ! ---------------------------------------  Same_Shape_RealR4_3d  -----
  subroutine Same_Shape_RealR4_3d ( Ref, New, ItsName, ModuleName )
    real, pointer, dimension(:,:,:) :: Ref, New
    include "Same_Shape_md.f9h"
  end subroutine Same_Shape_RealR4_3d
  ! ---------------------------------------  Same_Shape_RealR4_4d  -----
  subroutine Same_Shape_RealR4_4d ( Ref, New, ItsName, ModuleName )
    real, pointer, dimension(:,:,:,:) :: Ref, New
    include "Same_Shape_md.f9h"
  end subroutine Same_Shape_RealR4_4d
  ! ----------------------------------  Same_Shape_Character_1d_a  -----
  subroutine Same_Shape_Character_1d_a ( New, ItsName, ModuleName, TheShape )
    character(len=*), pointer, dimension(:) :: New
    include "Same_Shape_1d_a.f9h"
  end subroutine Same_Shape_Character_1d_a
  ! ----------------------------------  Same_Shape_Character_2d_a  -----
  subroutine Same_Shape_Character_2d_a ( New, ItsName, ModuleName, TheShape )
    character(len=*), pointer, dimension(:,:) :: New
    include "Same_Shape_md_a.f9h"
  end subroutine Same_Shape_Character_2d_a
  ! ----------------------------------  Same_Shape_Character_3d_a  -----
  subroutine Same_Shape_Character_3d_a ( New, ItsName, ModuleName, TheShape )
    character(len=*), pointer, dimension(:,:,:) :: New
    include "Same_Shape_md_a.f9h"
  end subroutine Same_Shape_Character_3d_a
  ! ------------------------------------  Same_Shape_Complex_1d_a  -----
  subroutine Same_Shape_Complex_1d_a ( New, ItsName, ModuleName, TheShape )
    complex, pointer, dimension(:) :: New
    include "Same_Shape_1d_a.f9h"
  end subroutine Same_Shape_Complex_1d_a
  ! ------------------------------------  Same_Shape_Complex_2d_a  -----
  subroutine Same_Shape_Complex_2d_a ( New, ItsName, ModuleName, TheShape )
    complex, pointer, dimension(:,:) :: New
    include "Same_Shape_md_a.f9h"
  end subroutine Same_Shape_Complex_2d_a
  ! ------------------------------------  Same_Shape_Complex_3d_a  -----
  subroutine Same_Shape_Complex_3d_a ( New, ItsName, ModuleName, TheShape )
    complex, pointer, dimension(:,:,:) :: New
    include "Same_Shape_md_a.f9h"
  end subroutine Same_Shape_Complex_3d_a
  ! -----------------------------------  Same_Shape_Dcomplex_1d_a  -----
  subroutine Same_Shape_Dcomplex_1d_a ( New, ItsName, ModuleName, TheShape )
    complex(kind(0.0d0)), pointer, dimension(:) :: New
    include "Same_Shape_1d_a.f9h"
  end subroutine Same_Shape_Dcomplex_1d_a
  ! -----------------------------------  Same_Shape_Dcomplex_2d_a  -----
  subroutine Same_Shape_Dcomplex_2d_a ( New, ItsName, ModuleName, TheShape )
    complex(kind(0.0d0)), pointer, dimension(:,:) :: New
    include "Same_Shape_md_a.f9h"
  end subroutine Same_Shape_Dcomplex_2d_a
  ! -----------------------------------  Same_Shape_Dcomplex_3d_a  -----
  subroutine Same_Shape_Dcomplex_3d_a ( New, ItsName, ModuleName, TheShape )
    complex(kind(0.0d0)), pointer, dimension(:,:,:) :: New
    include "Same_Shape_md_a.f9h"
  end subroutine Same_Shape_Dcomplex_3d_a
  ! -------------------------------------  Same_Shape_RealR8_1d_a  -----
  subroutine Same_Shape_RealR8_1d_a ( New, ItsName, ModuleName, TheShape )
    double precision, pointer, dimension(:) :: New
    include "Same_Shape_1d_a.f9h"
  end subroutine Same_Shape_RealR8_1d_a
  ! -------------------------------------  Same_Shape_RealR8_2d_a  -----
  subroutine Same_Shape_RealR8_2d_a ( New, ItsName, ModuleName, TheShape )
    double precision, pointer, dimension(:,:) :: New
    include "Same_Shape_md_a.f9h"
  end subroutine Same_Shape_RealR8_2d_a
  ! -------------------------------------  Same_Shape_RealR8_3d_a  -----
  subroutine Same_Shape_RealR8_3d_a ( New, ItsName, ModuleName, TheShape )
    double precision, pointer, dimension(:,:,:) :: New
    include "Same_Shape_md_a.f9h"
  end subroutine Same_Shape_RealR8_3d_a
  ! -------------------------------------  Same_Shape_RealR8_4d_a  -----
  subroutine Same_Shape_RealR8_4d_a ( New, ItsName, ModuleName, TheShape )
    double precision, pointer, dimension(:,:,:,:) :: New
    include "Same_Shape_md_a.f9h"
  end subroutine Same_Shape_RealR8_4d_a
  ! ------------------------------------  Same_Shape_Integer_1d_a  -----
  subroutine Same_Shape_Integer_1d_a ( New, ItsName, ModuleName, TheShape )
    integer, pointer, dimension(:) :: New
    include "Same_Shape_1d_a.f9h"
  end subroutine Same_Shape_Integer_1d_a
  ! ------------------------------------  Same_Shape_Integer_2d_a  -----
  subroutine Same_Shape_Integer_2d_a ( New, ItsName, ModuleName, TheShape )
    integer, pointer, dimension(:,:) :: New
    include "Same_Shape_md_a.f9h"
  end subroutine Same_Shape_Integer_2d_a
  ! ------------------------------------  Same_Shape_Integer_3d_a  -----
  subroutine Same_Shape_Integer_3d_a ( New, ItsName, ModuleName, TheShape )
    integer, pointer, dimension(:,:,:) :: New
    include "Same_Shape_md_a.f9h"
  end subroutine Same_Shape_Integer_3d_a
  ! ------------------------------------  Same_Shape_Integer_4d_a  -----
  subroutine Same_Shape_Integer_4d_a ( New, ItsName, ModuleName, TheShape )
    integer, pointer, dimension(:,:,:,:) :: New
    include "Same_Shape_md_a.f9h"
  end subroutine Same_Shape_Integer_4d_a
  ! ------------------------------------  Same_Shape_Logical_1d_a  -----
  subroutine Same_Shape_Logical_1d_a ( New, ItsName, ModuleName, TheShape )
    logical, pointer, dimension(:) :: New
    include "Same_Shape_1d_a.f9h"
  end subroutine Same_Shape_Logical_1d_a
  ! ------------------------------------  Same_Shape_Logical_2d_a  -----
  subroutine Same_Shape_Logical_2d_a ( New, ItsName, ModuleName, TheShape )
    logical, pointer, dimension(:,:) :: New
    include "Same_Shape_md_a.f9h"
  end subroutine Same_Shape_Logical_2d_a
  ! ------------------------------------  Same_Shape_Logical_3d_a  -----
  subroutine Same_Shape_Logical_3d_a ( New, ItsName, ModuleName, TheShape )
    logical, pointer, dimension(:,:,:) :: New
    include "Same_Shape_md_a.f9h"
  end subroutine Same_Shape_Logical_3d_a
  ! -------------------------------------  Same_Shape_RealR4_1d_a  -----
  subroutine Same_Shape_RealR4_1d_a ( New, ItsName, ModuleName, TheShape )
    real, pointer, dimension(:) :: New
    include "Same_Shape_1d_a.f9h"
  end subroutine Same_Shape_RealR4_1d_a
  ! -------------------------------------  Same_Shape_RealR4_2d_a  -----
  subroutine Same_Shape_RealR4_2d_a ( New, ItsName, ModuleName, TheShape )
    real, pointer, dimension(:,:) :: New
    include "Same_Shape_md_a.f9h"
  end subroutine Same_Shape_RealR4_2d_a
  ! -------------------------------------  Same_Shape_RealR4_3d_a  -----
  subroutine Same_Shape_RealR4_3d_a ( New, ItsName, ModuleName, TheShape )
    real, pointer, dimension(:,:,:) :: New
    include "Same_Shape_md_a.f9h"
  end subroutine Same_Shape_RealR4_3d_a
  ! -------------------------------------  Same_Shape_RealR4_4d_a  -----
  subroutine Same_Shape_RealR4_4d_a ( New, ItsName, ModuleName, TheShape )
    real, pointer, dimension(:,:,:,:) :: New
    include "Same_Shape_md_a.f9h"
  end subroutine Same_Shape_RealR4_4d_a

  ! ------------------------------=--- LogAllocate
  subroutine LogAllocate ( ModuleNameIn, ItsName, lBounds, uBounds, &
    & ElementSize, Address )
    character(len=*), intent(in)  :: ModuleNameIn, ItsName
    integer, intent(in), optional :: Lbounds(:), Ubounds(:)
    integer, intent(in), optional :: ElementSize
    integer(c_intptr_t), intent(in), optional :: Address ! of allocated object
                                                         ! if status == 0
    integer(c_intptr_t) :: Addr
    integer(c_intptr_t) :: Amount
    amount = 0
    addr = 0
    if ( present(address) ) addr = address
    if ( present(lbounds) .and. present(ubounds) .and. present(elementSize) ) then
      amount = memproduct(elementSize, ubounds-lbounds+1)
    end if
    serial_no = serial_no + 1
    if ( amount >= AllocateLogLimit ) then
      write( AllocateLogUnit, "(2(i0,1x),a,i0,2(1x,a))" ) addr, serial_no, 'A ', &
      & amount, trim(snipRCSfrom(ModuleNameIn)), trim(ItsName)
      flush ( AllocateLogUnit ) ! In case there's a crash soon after...
    end if
  end subroutine LogAllocate

  ! ------------------------------=--- LogDeallocate
  subroutine LogDeallocate ( ModuleNameIn, ItsName, Size, Address )
    character(len=*), intent(in)  :: ModuleNameIn, ItsName
    integer, intent(in), optional :: Size
    integer(c_intptr_t), intent(in), optional :: Address ! of deallocated object

    integer(c_intptr_t) :: Addr
    integer :: mySize
    ! Executable
    mySize = 0.
    addr = 0
    if ( present(address) ) addr = address
    if ( present(Size) ) mySize = Size
    serial_no = serial_no + 1
    if ( mySize >= AllocateLogLimit ) then
      write( AllocateLogUnit, "(2(i0,1x),a,i0,2(1x,a))" ) addr, serial_no, 'D ', &
        & mySize, trim(snipRCSfrom(ModuleNameIn)), trim(ItsName)
      flush ( AllocateLogUnit ) ! In case there's a crash soon after...
    end if
  end subroutine LogDeallocate

  ! ----------------------------------  memproduct  -----
  function memproduct ( elementSize, dimensions ) result( p )
    ! Find how many multiples of MEMORY_UNITS an array
    ! with elementSize and dimensions
    ! Note: we go through all this to forestall integer overflows
    integer, dimension(:), intent(in) :: dimensions
    integer, intent(in)               :: elementSize
    integer(c_intptr_t)               :: p
    !
    p = real(elementSize) * product(real(dimensions))
    ! print *, 'p ', p
  end function memproduct

  ! ------------------------------------------------  Not_Used_Here  -----
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Allocate_Deallocate

! $Log$
! Revision 2.58  2020/04/27 21:29:46  pwagner
! Separately track allocates/deallocates
!
! Revision 2.57  2019/08/19 21:58:18  pwagner
! Moved NoBytesAllocated to MLSCommon from Allocate_Deallocate; print warnings through MLSMessage
!
! Revision 2.56  2018/05/15 03:17:36  vsnyder
! Add 1D and 2D character allocatable versions
!
! Revision 2.55  2018/04/19 02:00:36  vsnyder
! Compute address for allocate/deallocate tracking.  Remove USE statements for
! unused names.
!
! Revision 2.54  2018/04/10 22:51:33  vsnyder
! Don't leave bounds(1:1) undefined for messages about scalars
!
! Revision 2.53  2017/10/31 23:42:57  vsnyder
! Add allocatable logical array support
!
! Revision 2.52  2016/09/07 23:54:37  vsnyder
! Add support for allocatable integer arrays
!
! Revision 2.51  2015/06/04 01:55:17  vsnyder
! Don't reallocate if bounds are the same
!
! Revision 2.50  2015/06/02 23:51:40  vsnyder
! Routines for allocatable arrays
!
! Revision 2.49  2015/03/27 01:10:38  vsnyder
! Log allocated/deallocated entity's address, if provided
!
! Revision 2.48  2014/09/04 23:37:01  vsnyder
! Change Test_Allocate to Test_Allocate_nd, Add Test_Allocate_1d, add
! generic for them.  Keep track of memory in bytes instead of Memory_Units
! (formerly 1024 bytes), delete Memory_Units.
!
! Revision 2.47  2014/08/06 23:18:04  vsnyder
! Add option controlled by private parameter to fill REAL arrays with sNaN
!
! Revision 2.46  2014/06/30 23:24:17  pwagner
! Can log allocations/deallocations to separate file
!
! Revision 2.45  2014/05/29 18:20:27  pwagner
! Extra debugging possibility
!
! Revision 2.44  2014/01/09 00:25:06  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.43  2013/10/09 01:02:05  vsnyder
! Add SnipRCSFromName to MyMessage
!
! Revision 2.42  2013/08/28 00:38:17  pwagner
! Added a local version of MyMessage to evade possible circular dependency
!
! Revision 2.41  2013/08/16 02:05:34  vsnyder
! Make Same_Shape public (oops), add Same_Shape..._a versions
!
! Revision 2.40  2013/08/16 01:06:22  vsnyder
! Add Same_Shape
!
! Revision 2.39  2013/06/12 02:50:45  vsnyder
! Repair a couple of harmless but ugly blunders
!
! Revision 2.38  2013/06/12 02:17:03  vsnyder
! Add and use BYTES and Byte_Size
!
! Revision 2.37  2011/03/02 01:59:07  vsnyder
! Make lbounds and ubounds arguments of Test_Allocate optional
!
! Revision 2.36  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.35  2008/05/21 21:49:26  vsnyder
! Add allocation routines with arrays for bounds
!
! Revision 2.34  2008/05/20 01:59:28  vsnyder
! Add 4d integer, real, double
!
! Revision 2.33  2007/07/27 00:19:30  vsnyder
! Better printing in Test_[De]Allocate, work on memory tracking
!
! Revision 2.32  2006/08/05 02:36:57  vsnyder
! Pass bounds argument of ReportAllocateDeallocate_ints forward
!
! Revision 2.31  2006/08/04 18:12:56  vsnyder
! Add 'bounds' to ReportAlllocateDeallocate, don't look at ElementSize if it's not present
!
! Revision 2.30  2006/07/29 03:01:07  vsnyder
! Use track_m stuff
!
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
