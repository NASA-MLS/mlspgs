!--------------------------------------------------------------------
!
      SUBROUTINE LOCATE(XX,N,NP,X,J)
      use MLSCommon, only: r8
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
      RETURN
      END

! $Log: interp_pack.f90,v      





