! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

! $Id$

C === (start of toc) ===
C ddot      Returns dot product of two double-precision vectors
C === (end of toc) ===

C === (start of api) ===
C dble ddot ( int n, dble x(:), int incx, dble y(:), int incy )
C === (end of api) ===

      DOUBLE PRECISION FUNCTION DDOT(N,X,INCX,Y,INCY)
C>> 1994-11-11 DDOT  Krogh   Declared all vars.
c>> 1994-10-20 DDOT   Krogh  Changes to use M77CON
c>> 1994-04-19 DDOT   Krogh   Minor -- Made code versions line up.
C>> 1985-08-02 DDOT   Lawson  Initial code.
c--D replaces "?": ?DOT
C
C     RETURNS THE DOT PRODUCT OF X AND Y.
C     DDOT = SUM FOR I = 0 TO N-1 OF  X(LX+I*INCX) * Y(LY+I*INCY),
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
C

      INTEGER N, INCX, INCY, IX, IY, I, M, MP1, NS
      DOUBLE PRECISION X(*),Y(*)
      DDOT = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
         DDOT = DDOT + X(IX)*Y(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
         DDOT = DDOT + X(I)*Y(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
         DDOT = DDOT + X(I)*Y(I) + X(I+1)*Y(I+1) +
     $   X(I + 2)*Y(I + 2) + X(I + 3)*Y(I + 3) + X(I + 4)*Y(I + 4)
   50 CONTINUE
      RETURN
C
C         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DDOT = DDOT + X(I)*Y(I)
   70     CONTINUE
      RETURN
      END
! $Log$
! Revision 1.3  2002/10/10 23:48:57  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 1.2  2001/03/06 18:17:56  pwagner
! clean more amitious--rms * from MLSCONFGl[13]/MakeFC
!
