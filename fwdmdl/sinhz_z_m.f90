! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module SINHZ_Z_M

  implicit NONE
  private
  public SINHZ_Z, D_SINHZ_Z3

  interface SINHZ_Z; module procedure D_SINHZ_Z, S_SINHZ_Z; end interface
  interface D_SINHZ_Z3; module procedure D_D_SINHZ_Z3, S_D_SINHZ_Z3; end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    &  "$RCSfile$"
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

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module SINHZ_Z_M

! $Log$
! Revision 1.1.2.1  2003/04/30 01:46:29  vsnyder
! Initial commit
!
