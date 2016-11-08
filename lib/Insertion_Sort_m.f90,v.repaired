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
module Insertion_Sort_m
!=============================================================================

  implicit NONE

  public

  interface Insertion_Sort
    module procedure Insertion_Sort_P_D, Insertion_Sort_P_S
    module procedure Insertion_Sort_D, Insertion_Sort_S
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine Insertion_Sort_P_D ( R, P )
  ! Sort R, producing the permutation P that sorted R
    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(inout) :: R(0:)
    integer, intent(out) :: P(:)    ! Assumed to be the same size as R(1:)

    integer :: I, IP, J
    real(rk) :: RZ

    r(0) = -0.5 * huge(1.0_rk)      ! Sentinel simplifies algorithm
    p = [ ( i, i = 1, ubound(r,1) ) ]

    do i = 2, ubound(r,1)
      rz = r(i)
      ip = p(i)
      j = i
      do while ( rz < r(j-1) )
        r(j) = r(j-1)
        p(j) = p(j-1)
        j = j - 1
      end do
      r(j) = rz
      p(j) = ip
    end do

  end subroutine Insertion_Sort_P_D

  subroutine Insertion_Sort_P_S ( R, P )
  ! Sort R, producing the permutation P that sorted R
    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(inout) :: R(0:)
    integer, intent(out) :: P(:)    ! Assumed to be the same size as R(1:)

    integer :: I, IP, J
    real(rk) :: RZ

    r(0) = -0.5 * huge(1.0_rk)      ! Sentinel simplifies algorithm
    p = [ ( i, i = 1, ubound(r,1) ) ]

    do i = 2, ubound(r,1)
      rz = r(i)
      ip = p(i)
      j = i
      do while ( rz < r(j-1) )
        r(j) = r(j-1)
        p(j) = p(j-1)
        j = j - 1
      end do
      r(j) = rz
      p(j) = ip
    end do

  end subroutine Insertion_Sort_P_S

  subroutine Insertion_Sort_D ( R )
  ! Sort R.
    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(inout) :: R(0:)

    integer :: I, J
    real(rk) :: RZ

    r(0) = -0.5 * huge(1.0_rk)      ! Sentinel simplifies algorithm

    do i = 2, ubound(r,1)
      rz = r(i)
      j = i
      do while ( rz < r(j-1) )
        r(j) = r(j-1)
        j = j - 1
      end do
      r(j) = rz
    end do

  end subroutine Insertion_Sort_D

  subroutine Insertion_Sort_S ( R )
  ! Sort R.
    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(inout) :: R(0:)

    integer :: I, J
    real(rk) :: RZ

    r(0) = -0.5 * huge(1.0_rk)      ! Sentinel simplifies algorithm

    do i = 2, ubound(r,1)
      rz = r(i)
      j = i
      do while ( rz < r(j-1) )
        r(j) = r(j-1)
        j = j - 1
      end do
      r(j) = rz
    end do

  end subroutine Insertion_Sort_S

!=============================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Insertion_Sort_m

! $Log$
! Revision 2.1  2016/11/08 21:15:45  vsnyder
! Initial commit
!
