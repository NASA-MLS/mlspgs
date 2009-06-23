! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module SINHZ_Z_M

  implicit NONE
  private
  public SINHZ_Z, D_SINHZ_Z3

  interface SINHZ_Z; module procedure D_SINHZ_Z, S_SINHZ_Z; end interface
  interface D_SINHZ_Z3; module procedure D_D_SINHZ_Z3, S_D_SINHZ_Z3; end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! --------------------------------------------------  D_SINHZ_Z  -----
  pure function D_SINHZ_Z ( Z ) result ( SINHZ_Z )
  ! Compute sinh(z)/z, taking care when Z gets near zero.
  ! Write sinh(z)/z as exp(-z) (exp(2*z)-1)/z.

    use CRREXP_M, only: CRREXP
    integer, parameter :: RK = kind(0.0d0)

    complex(rk), intent(in) :: Z
    complex(rk) :: SINHZ_Z

    sinhz_z = exp(-z) * crrexp(z+z)
  end function D_SINHZ_Z

  ! --------------------------------------------------  S_SINHZ_Z  -----
  pure function S_SINHZ_Z ( Z ) result ( SINHZ_Z )
  ! Compute sinh(z)/z, taking care when Z gets near zero.
  ! Write sinh(z)/z as exp(-z) (exp(2*z)-1)/z.

    use CRREXP_M, only: CRREXP
    integer, parameter :: RK = kind(0.0e0)

    complex(rk), intent(in) :: Z
    complex(rk) :: SINHZ_Z

    sinhz_z = exp(-z) * crrexp(z+z)
  end function S_SINHZ_Z

  ! -----------------------------------------------  D_D_SINHZ_Z3  -----
  pure function D_D_SINHZ_Z3 ( Z ) result ( D_SINHZ_Z3 )
  ! Compute the derivative of sinh(z)/z, divided by z, which is
  ! (z cosh(z) - sinh(z)) / z**3, taking care when |Z| is small.

    integer, parameter :: RK = kind(0.0d0)
    complex(rk), intent(in) :: Z
    complex(rk) :: D_SINHZ_Z3

    include 'd_sinhz_z3.f9h'

  end function D_D_SINHZ_Z3

  ! -----------------------------------------------  S_D_SINHZ_Z3  -----
  pure function S_D_SINHZ_Z3 ( Z ) result ( D_SINHZ_Z3 )
  ! Compute the derivative of sinh(z)/z, divided by z, which is
  ! (z cosh(z) - sinh(z)) / z**3, taking care when |Z| is small.

    integer, parameter :: RK = kind(0.0e0)
    complex(rk), intent(in) :: Z
    complex(rk) :: D_SINHZ_Z3

    include 'd_sinhz_z3.f9h'

  end function S_D_SINHZ_Z3

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module SINHZ_Z_M

! $Log$
! Revision 2.2  2005/06/22 18:08:20  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.1  2003/05/05 23:00:26  livesey
! Merged in feb03 newfwm branch
!
! Revision 1.1.2.1  2003/04/30 01:46:29  vsnyder
! Initial commit
!
