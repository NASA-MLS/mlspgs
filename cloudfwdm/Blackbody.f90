! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Blackbody

! ------------------------------------------------------
! USE PLANCK FUNCTION TO COMPUTE BRIGHTNESS TEMPERATURE
! ------------------------------------------------------

        use MLSCommon, only: r8
        IMPLICIT NONE
        Private
        Public :: planck

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
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

end module Blackbody

! $Log$
