!=================================
PROGRAM blastest ! tests blas routines
!=================================

   USE MLSCommon , ONLY: r4, r8
   USE time_m , ONLY: wait, time_now, WAIT_LOOP_LIMITS, &
   & RETRY, INIT_RETRY, TRY_AGAIN, RETRY_SUCCESS, TOO_MANY_RETRIES   
   
   IMPLICIT NONE

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This program tests the blas subroutines.

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it

! Variables
   integer, parameter          :: ntests = 6
   real(r8), dimension(ntests) :: times
   integer, dimension(ntests)  :: lengths
   integer                     :: test, loops

   lengths = (/ 5000, 7500, 10000, 12500, 15000, 17500 /)
   print *, 'How many loops (e.g., 100)'
   read *, loops
   lengths = (/ 500, 750, 1000, 1250, 1500, 2000 /)
  ! Matrix - vector products
	do test=1, ntests
 	  times(test) = gemv_test(lengths(test), loops)
	enddo
   print *, 'Test results (gemv)'
   print *, 'Test      Length          Time          time/mflop'
   do test=1, ntests
     print*, test, lengths(test), times(test), 1d6*times(test)/lengths(test)**2
   enddo
   lengths = (/ 50, 75, 100, 125, 150, 200 /)
  ! Matrix - matrix products (r4)
	do test=1, ntests
 	  times(test) = sgemm_test(lengths(test), loops)
	enddo
   print *, 'Test results (sgemm)'
   print *, 'Test      Length          Time          time/mflop'
   do test=1, ntests
     print*, test, lengths(test), times(test), 1d6*times(test)/lengths(test)**3
   enddo
   lengths = (/ 50, 75, 100, 125, 150, 200 /)
  ! Matrix - matrix products (r8)
	do test=1, ntests
 	  times(test) = dgemm_test(lengths(test), loops)
	enddo
   print *, 'Test results (dgemm)'
   print *, 'Test      Length          Time          time/mflop'
   do test=1, ntests
     print*, test, lengths(test), times(test), 1d6*times(test)/lengths(test)**3
   enddo
contains
  subroutine fill_vector ( vector, startv, endv )
    ! Fill vector from vector(1) to vector (length) with values
    ! ranging from startv to endv linearly
    real(r8), dimension(:), intent(out)   :: vector
    real(r8), intent(in)                  :: startv
    real(r8), intent(in)                  :: endv
    
    ! Private
    real(r8) :: dv
    integer ::   i, length
    if (size(vector) <= 0) return
    length = size(vector)
    vector(1) = startv
    if ( length == 1) return
    dv = (endv - startv) / (length-1)
    do i = 2, length
      vector(i) = startv + (i-1)*dv
    enddo
  end subroutine fill_vector

  subroutine fill_svector ( vector, startv, endv )
    ! Fill vector from vector(1) to vector (length) with values
    ! ranging from startv to endv linearly
    real(r4), dimension(:), intent(out)   :: vector
    real(r4), intent(in)                  :: startv
    real(r4), intent(in)                  :: endv
    
    ! Private
    real(r4) :: dv
    integer ::   i, length
    if (size(vector) <= 0) return
    length = size(vector)
    vector(1) = startv
    if ( length == 1) return
    dv = (endv - startv) / (length-1)
    do i = 2, length
      vector(i) = startv + (i-1)*dv
    enddo
  end subroutine fill_svector

  function gemv_test ( length, loops )
    ! Time how long it takes to perform matrix-vector product of length length
    integer, intent(in)   :: length, loops
    real(r8)                  :: gemv_test
    real(r8) :: t0, t1
    integer :: i
    
    ! Private
    real(r8), parameter :: startv=-1.d0, endv = 1.d0
    real(r8), dimension(:), pointer   :: v1, v2
    real(r8), dimension(:,:), pointer :: m
    real(r8) :: v1sq, v2sq, v1dotv2, v2dotv1
    
    allocate(v1(length), v2(length), m(length, length))
    call fill_vector(v1, startv, endv)
    call fill_vector(v2, startv, endv)
    do i=1, length
      call fill_vector(m(:,i), startv, endv)
    enddo
    call time_now(t0)
    do i=1, loops
    call dgemv('n', length, length, 1.d0, m, length, v1, 1, &
      & 1.d0, v2, 1)
    enddo
    call time_now(t1)
    gemv_test = t1 - t0
    deallocate(v1, v2, m)
  end function gemv_test

  function dgemm_test ( length, loops )
    ! Time how long it takes to perform matrix-matrix product of length length
    ! (type r8)
    integer, intent(in)   :: length, loops
    real(r8)                  :: dgemm_test
    real(r8) :: t0, t1
    integer :: i
    
    ! Private
    real(r8), parameter :: startv=-1.d0, endv = 1.d0
    real(r8), dimension(:,:), pointer   :: m1, m2
    real(r8), dimension(:,:), pointer :: m
    real(r8) :: v1sq, v2sq, v1dotv2, v2dotv1
    
    allocate(m1(length, length), m2(length, length), m(length, length))
    do i=1, length
      call fill_vector(m(:,i), startv, endv)
      call fill_vector(m1(:,i), startv, endv)
      call fill_vector(m2(:,i), startv, endv)
    enddo
    call time_now(t0)
    do i=1, loops
    call dgemm('n', 'n', length, length, length, 1.d0, m1, length, m2, length, &
      & 1.d0, m, length)
    enddo
    call time_now(t1)
    dgemm_test = t1 - t0
    deallocate(m1, m2, m)
  end function dgemm_test

  function sgemm_test ( length, loops )
    ! Time how long it takes to perform matrix-vector product of length length
    ! (type r4)
    integer, intent(in)   :: length, loops
    real(r4)                  :: sgemm_test
    real(r8) :: t0, t1
    integer :: i
    
    ! Private
    real(r4), parameter :: startv=-1.0, endv = 1.0
    real(r4), dimension(:,:), pointer   :: m1, m2
    real(r4), dimension(:,:), pointer :: m
    real(r4) :: v1sq, v2sq, v1dotv2, v2dotv1
    
    allocate(m1(length, length), m2(length, length), m(length, length))
    do i=1, length
      call fill_svector(m(:,i), startv, endv)
      call fill_svector(m1(:,i), startv, endv)
      call fill_svector(m2(:,i), startv, endv)
    enddo
    call time_now(t0)
    do i=1, loops
    call sgemm('n', 'n', length, length, length, 1.0, m1, length, m2, length, &
      & 1.0, m, length)
    enddo
    call time_now(t1)
    sgemm_test = t1 - t0
    deallocate(m1, m2, m)
  end function sgemm_test

!==================
END PROGRAM blastest
!==================

!# $Log$
!#
