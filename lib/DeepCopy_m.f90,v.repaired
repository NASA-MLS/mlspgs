! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module DeepCopy_m

  ! call deepCopy ( a, b ) does a deep copy of b to a by destroying a,
  ! cloning a to be like b, and then, if a is not disassociated,
  ! copying the value of b to a.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test

  implicit NONE
  private
  public :: DeepCopy

  interface DeepCopy
    module procedure DeepCopy_C1, DeepCopy_C2, DeepCopy_C3 ! Complex
    module procedure DeepCopy_D1, DeepCopy_D2, DeepCopy_D3 ! Double Precision
    module procedure DeepCopy_I1, DeepCopy_I2, DeepCopy_I3 ! Integer
    module procedure DeepCopy_R1, DeepCopy_R2, DeepCopy_R3 ! Real
  end interface

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains

  subroutine DeepCopy_C1 ( A, B )
    complex, pointer :: A(:), B(:)
    include 'deepCopy_1.f9h'
  end subroutine DeepCopy_C1

  subroutine DeepCopy_D1 ( A, B )
    double precision, pointer :: A(:), B(:)
    include 'deepCopy_1.f9h'
  end subroutine DeepCopy_D1

  subroutine DeepCopy_I1 ( A, B )
    integer, pointer :: A(:), B(:)
    include 'deepCopy_1.f9h'
  end subroutine DeepCopy_I1

  subroutine DeepCopy_R1 ( A, B )
    real, pointer :: A(:), B(:)
    include 'deepCopy_1.f9h'
  end subroutine DeepCopy_R1

  subroutine DeepCopy_C2 ( A, B )
    complex, pointer :: A(:,:), B(:,:)
    include 'deepCopy.f9h'
  end subroutine DeepCopy_C2

  subroutine DeepCopy_D2 ( A, B )
    double precision, pointer :: A(:,:), B(:,:)
    include 'deepCopy.f9h'
  end subroutine DeepCopy_D2

  subroutine DeepCopy_I2 ( A, B )
    integer, pointer :: A(:,:), B(:,:)
    include 'deepCopy.f9h'
  end subroutine DeepCopy_I2

  subroutine DeepCopy_R2 ( A, B )
    real, pointer :: A(:,:), B(:,:)
    include 'deepCopy.f9h'
  end subroutine DeepCopy_R2

  subroutine DeepCopy_C3 ( A, B )
    complex, pointer :: A(:,:,:), B(:,:,:)
    include 'deepCopy.f9h'
  end subroutine DeepCopy_C3

  subroutine DeepCopy_D3 ( A, B )
    double precision, pointer :: A(:,:,:), B(:,:,:)
    include 'deepCopy.f9h'
  end subroutine DeepCopy_D3

  subroutine DeepCopy_I3 ( A, B )
    integer, pointer :: A(:,:,:), B(:,:,:)
    include 'deepCopy.f9h'
  end subroutine DeepCopy_I3

  subroutine DeepCopy_R3 ( A, B )
    real, pointer :: A(:,:,:), B(:,:,:)
    include 'deepCopy.f9h'
  end subroutine DeepCopy_R3

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module DeepCopy_m

! $Log$
! Revision 2.2  2012/07/18 22:13:13  vsnyder
! Improve a comment
!
! Revision 2.1  2012/07/10 03:10:40  vsnyder
! Initial commit
!
