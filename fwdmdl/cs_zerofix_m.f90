! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module CS_ZeroFix_M

!--------------------------------------------------------------------
!  Fix a complex number such that small parts are zeroed out

  implicit NONE
  private
  public :: CS_ZeroFix

  interface CS_ZeroFix
    module procedure CS_ZeroFix_D, CS_ZeroFix_S
  end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    &  "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  elemental function CS_ZeroFix_D ( W ) result ( V )

    integer, parameter :: RK = kind(0.0d0)

    complex(rk), intent(in) :: W
    complex(rk) :: V

    include "cs_zerofix.f9h"

  end function CS_ZeroFix_D

  elemental function CS_ZeroFix_S ( W ) result ( V )

    integer, parameter :: RK = kind(0.0e0)

    complex(rk), intent(in) :: W
    complex(rk) :: V

    include "cs_zerofix.f9h"

  end function CS_ZeroFix_S

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module CS_ZeroFix_M

! $Log$
! Revision 2.1  2003/02/04 01:41:06  vsnyder
! Initial commit
!

