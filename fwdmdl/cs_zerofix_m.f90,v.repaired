! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module CS_ZeroFix_M

!--------------------------------------------------------------------
!  Fix a complex number such that small parts are zeroed out

  implicit NONE
  private
  public :: CS_ZeroFix

  interface CS_ZeroFix
    module procedure CS_ZeroFix_D, CS_ZeroFix_S
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
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

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module CS_ZeroFix_M

! $Log$
! Revision 2.3  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.2  2005/06/22 18:08:18  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.1  2003/02/04 01:41:06  vsnyder
! Initial commit
!

