c note:
c	np	- physical dimension
c	n	- old array dimension
c	m	- new array dimension
c
c--------------------------------------------------------------------
c
	subroutine lin_interp(x,y,np,n,x1,y1)
	integer n,j,np
	real*4 x(np), y(np)
	real*4 x1, y1

	call locate(x,n,np,x1,j)
	y1 = ((x(j+1)-x1)*y(j) + (x1-x(j))*y(j+1))/(x(j+1)-x(j))

	return
	end
c
c--------------------------------------------------------------------
c
      SUBROUTINE LOCATE(XX,N,NP,X,J)
	integer jl,ju,n,np,jm,j
      real*4 XX(NP),x

	if(x .lt. xx(1)) then
	j = 1
	return
	endif

	if(x .gt. xx(n)) then
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

! $Log: interp_pack.f,v      
