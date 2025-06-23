
! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!==============================================================================
Module NumRecipesModule
!==============================================================================
	
	PUBLIC


!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! Contents:

! Subroutines -- Init 
!                Loc2RS
!                cordTransform
!                findad
!                FFSM
!                FFSM1
!                Reconstruct
!                Diagnostics
!                DataGenerate
!                DataGenerate1
! Function -- DataField

! Remarks:  This is a module for the main Core processing.

! Parameters

	Save

Contains

! /---------------------------------------------------------------\
! |*** SPLINE
! \---------------------------------------------------------------/

      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
      PARAMETER(NMAX=1000)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
      IF (YP1.GT..99E30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))	&
           /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99E30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      End Subroutine SPLINE


! /---------------------------------------------------------------\
! |*** SPLINT
! \---------------------------------------------------------------/

      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      DIMENSION XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.0) PAUSE 'Bad XA input.'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+	&
            ((A-A**3)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.0
!    *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.0
      RETURN
      END Subroutine SPLINT


!===================
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here
End Module NumRecipesModule
!===================

! $Log$
! Revision 1.4  2005/06/23 19:07:39  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 1.3  2003/04/06 19:32:34  jdone
! minor header editing
!
! Revision 1.2  2001/08/06 17:04:05  ybj
! fix the array size
!
! Revision 1.1  2001/02/27 20:54:21  ybj
! Interpolation Routines From Numerical Recipes
!
! Revision 1.1  2000/10/05 18:17:41  nakamura
! Module split from synoptic.f90 and modified to be more like the standard template.
!
