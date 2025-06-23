! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Rotation_m

  ! Rotate a three-dimensional vector through a specified angle about
  ! another three-dimensional vector.

  implicit none
  private

  public :: Rotate_3d

  interface Rotate_3d
    module procedure Rotate_3d_D, Rotate_3d_s
  end interface

  ! Kronecker delta
  integer, parameter :: D(3,3) = reshape( [ 1, 0, 0, &
                                          & 0, 1, 0, &
                                          & 0, 0, 1 ], [ 3, 3 ] )
  ! Index of element of About_Vector to multiply by Levi_Civita and sin(phi)
  ! The diagonal element is arbitrary (but it must be in 1..3) because that
  ! element of About_Vector is multiplied by zero.
  integer, parameter :: K(3,3) = reshape( [ 1, 3, 2, &
                                          & 3, 1, 1, &
                                          & 2, 1, 1  ], [ 3, 3 ] )
  ! (i,j,k(i,j)) element of the Levi_Civita symbol (others are zero)
  integer, parameter :: L(3,3) = reshape( [ 0, -1, 1, &
                                          & 1, 0, -1, &
                                          & -1, 1, 0  ], [ 3, 3 ] )

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains ! ============= Public Procedures ==========================

  subroutine Rotate_3d_D ( Original_Vector, Phi, About_Vector, Rotated_Vector )

    ! Rotate Original_Vector by Phi radians about About_Vector giving
    ! Rotated_Vector.

    integer, parameter :: RK = kind(0.0d0)

    real(rk), intent(in) :: Original_Vector(3) ! in XYZ coordinates
    real(rk), intent(in) :: Phi                ! to rotate, Radians
    real(rk), intent(in) :: About_Vector(3)    ! in XYZ coordinates
                                               ! An unit vector parallel to it
                                               ! is constructed.
    real(rk), intent(out) :: Rotated_Vector(3) ! in XYZ coordinates

    include "Rotate_3d.f9h"

  end subroutine Rotate_3d_D

  subroutine Rotate_3d_S ( Original_Vector, Phi, About_Vector, Rotated_Vector )

    ! Rotate Original_Vector by Phi about About_Vector giving Rotated_Vector.

    integer, parameter :: RK = kind(0.0e0)

    real(rk), intent(in) :: Original_Vector(3) ! in XYZ coordinates
    real(rk), intent(in) :: Phi                ! to rotate, Radians
    real(rk), intent(in) :: About_Vector(3)    ! in XYZ coordinates
                                               ! An unit vector parallel to it
                                               ! is constructed.
    real(rk), intent(out) :: Rotated_Vector(3) ! in XYZ coordinates

    include "Rotate_3d.f9h"

  end subroutine Rotate_3d_S

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Rotation_m

! $Log$
! Revision 2.3  2015/02/05 21:41:14  vsnyder
! Move dummy argument declarations from the include file.
!
! Revision 2.2  2014/11/14 23:57:56  vsnyder
! Correct some comments, add comments about arguments
!
! Revision 2.1  2014/11/04 01:25:45  vsnyder
! Initial commit
!
