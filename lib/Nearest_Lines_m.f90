! Copyright 2015, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module Nearest_Lines_m
!=============================================================================

  ! Compute the points on two lines that are nearest to each other.
  ! For derivation, see wvs-136.

  implicit NONE
  private

  public :: Nearest_Lines

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine Nearest_Lines ( L1, L2, S1, S2, Parallel )
    use Geolocation_0, only: ECR_t, RG
    type(ECR_t), intent(inout) :: L1(2), L2(2) ! Given lines, as C + s * U
                                               ! where C is a point on the
                                               ! line and U is an unit vector
                                               ! along the line.  Just in case,
                                               ! U are made unit vectors here.
    real(rg), intent(out) :: S1, S2            ! L1(1) + S1 * L1(2) and
                                               ! L2(1) + S2 * L2(2) are points
                                               ! on L1 and L2 that are nearest
                                               ! to each other, unless L1 and
                                               ! L2 are parallel, in which case
                                               ! S1 = 0 and
                                               ! |L1(1) - L2(1) - S2 * L2(2)|
                                               ! is the distance between them
    logical, intent(out) :: Parallel           ! "L1 and L2 are parallel"

    type(ECR_t) :: C
    real(rg) :: B, D, E ! Various dot products, see wvs-136
    real(rg) :: Z       ! 1 - b**2

    l1(2) = l1(2) / l1(2)%norm2()
    l2(2) = l2(2) / l2(2)%norm2()
    b = l1(2) .dot. l2(2)         ! U1 .dot. U2
    z = 1 - b*b
    c = l1(1) - l2(1)             ! Cw = P2 - P1
    e = l2(2) .dot. c             ! U2 .dot. Cw
    parallel = abs(z) < sqrt(tiny(1.0_rg))
    if ( parallel ) then
      s1 = 0
      s2 = e
    else
      d = l1(2) .dot. c           ! U1 .dot. Cw
      s1 = ( b * e - d ) / z
      s2 = ( e - b * d ) / z
    end if

  end subroutine Nearest_Lines

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

end module Nearest_Lines_m

! $Log$
! Revision 2.1  2016/11/03 01:54:30  vsnyder
! Initial commit
!
