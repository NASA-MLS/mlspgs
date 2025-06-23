! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Cross_m

  ! Compute the cross product of two three-dimensional vectors.

  implicit none
  private

  public :: Cross

  interface Cross
    procedure Cross_r4, Cross_r8, Cross2_r4, Cross2_r8
  end interface

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains ! ============= Public Procedures ==========================

  pure function Cross_r4 ( XYZ, Norm ) result ( Cross )
    real, intent(in) :: XYZ(3,2)
    logical, intent(in), optional :: Norm ! Normalize a nonzero result
    real :: Cross(3)
    real :: N
    cross = (/  xyz(2,1)*xyz(3,2) - xyz(3,1)*xyz(2,2),  &
              -(xyz(1,1)*xyz(3,2) - xyz(3,1)*xyz(1,2)), &
                xyz(1,1)*xyz(2,2) - xyz(2,1)*xyz(1,2) /)
    if ( present(norm) ) then
      if ( norm ) then
        n = norm2(cross)
        if ( n /= 0.0 ) cross = cross / norm2(cross)
      end if
    end if
  end function Cross_r4

  pure function Cross_r8 ( XYZ, Norm ) result ( Cross )
    double precision, intent(in) :: XYZ(3,2)
    logical, intent(in), optional :: Norm ! Normalize a nonzero result
    double precision :: Cross(3)
    double precision :: N
    cross = (/  xyz(2,1)*xyz(3,2) - xyz(3,1)*xyz(2,2),  &
              -(xyz(1,1)*xyz(3,2) - xyz(3,1)*xyz(1,2)), &
                xyz(1,1)*xyz(2,2) - xyz(2,1)*xyz(1,2) /)
    if ( present(norm) ) then
      if ( norm ) then
        n = norm2(cross)
        if ( n /= 0.0 ) cross = cross / norm2(cross)
      end if
    end if
  end function Cross_r8

  pure function Cross2_r4 ( A, B, Norm ) result ( Cross )
    real, intent(in) :: A(3), B(3)
    logical, intent(in), optional :: Norm ! Normalize a nonzero result
    real :: Cross(3)
    real :: N
    cross = (/  a(2)*b(3) - a(3)*b(2),  &
              -(a(1)*b(3) - a(3)*b(1)), &
                a(1)*b(2) - a(2)*b(1) /)
    if ( present(norm) ) then
      if ( norm ) then
        n = norm2(cross)
        if ( n /= 0.0 ) cross = cross / norm2(cross)
      end if
    end if
  end function Cross2_r4

  pure function Cross2_r8 ( A, B, Norm ) result ( Cross )
    double precision, intent(in) :: A(3), B(3)
    logical, intent(in), optional :: Norm ! Normalize a nonzero result
    double precision :: Cross(3)
    double precision :: N
    cross = (/  a(2)*b(3) - a(3)*b(2),  &
              -(a(1)*b(3) - a(3)*b(1)), &
                a(1)*b(2) - a(2)*b(1) /)
    if ( present(norm) ) then
      if ( norm ) then
        n = norm2(cross)
        if ( n /= 0.0 ) cross = cross / norm2(cross)
      end if
    end if
  end function Cross2_r8

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Cross_m

! $Log$
! Revision 2.1  2014/12/12 01:16:43  vsnyder
! Initial commit
!
