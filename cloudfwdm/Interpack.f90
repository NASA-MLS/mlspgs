! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Interpack

      use MLSCommon, only: r8
      IMPLICIT NONE
      Private
      Public :: LOCATE

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
 !---------------------------------------------------------------------------
      
contains

!--------------------------------------------------------------------
!
      SUBROUTINE LOCATE(XX,N,NP,X,J)
      integer :: jl,ju,n,np,jm,j
      real(r8) :: XX(NP),x

	if(x .le. xx(1)) then
	j = 1
	return
	endif

	if(x .ge. xx(n) ) then
	j = n-1
	return
	endif

      JL=0
      JU=N+1
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GO TO 10
      ENDIF
      J=JL

      END SUBROUTINE LOCATE

end module Interpack

! $Log: Interpack.f90,v      





