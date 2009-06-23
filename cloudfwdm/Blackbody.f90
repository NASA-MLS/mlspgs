! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Blackbody

! ------------------------------------------------------
! USE PLANCK FUNCTION TO COMPUTE BRIGHTNESS TEMPERATURE
! ------------------------------------------------------

        use MLSCommon, only: r8
        IMPLICIT NONE
        Private
        Public :: planck

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
      
contains

      subroutine planck(temp,freq,tb)
        real(r8) :: temp,tb,freq
        real(r8) :: h
        real(r8) :: k

        h = 6.6256
        k = 1.3805

        tb=h*freq*1.e-2_r8/(exp(h*freq*1.e-2_r8/k/temp)-1.)/k

      end subroutine planck

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Blackbody

! $Log$
! Revision 1.5  2007/10/04 18:20:32  vsnyder
! Remove tabs, which are not allowed by the Fortran standard
!
! Revision 1.4  2005/06/22 18:27:38  pwagner
! Cant have access declared outside module scope
!
! Revision 1.3  2002/10/08 17:08:07  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 1.2  2001/09/21 15:51:37  jonathan
! modified F95 version
!
