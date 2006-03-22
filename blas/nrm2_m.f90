
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module NRM2_M ! From Math77, converted to Fortran 90
!=============================================================================

  implicit NONE
  private

  public :: NRM2, SNRM2, DNRM2

  interface NRM2
    module procedure SNRM2, DNRM2
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName="$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  function SNRM2 ( X ) result ( NRM2 )

    integer, parameter :: RK = kind(0.0)
    real(rk), intent(in) :: X(:)
    real(rk) :: NRM2
    include 'nrm2.f9h'

  end function SNRM2

  function DNRM2 ( X ) result ( NRM2 )

    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: X(:)
    real(rk) :: NRM2
    include 'nrm2.f9h'

  end function DNRM2

!=======================================================================
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module NRM2_M

! $Log$
