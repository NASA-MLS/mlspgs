C Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
C U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

C === (start of toc) ===
C mdsdot      Returns dot product of a double- by a single-precision vector
C msddot      Returns dot product of a single- by a double-precision vector
C === (end of toc) ===

C === (start of api) ===
C dble mdsdot ( int n, dble x(:), int incx, real y(:), int incy )
C dble msddot ( int n, real x(:), int incx, dble y(:), int incy )
C === (end of api) ===

C > >       DOUBLE PRECISION      FUNCTION MSDDOT(N,X,INCX,Y,INCY)
C > > C Dot product of mixed type input
C > >       INTEGER N, INCX, INCY, IX, IY, I, M, MP1, NS
C > >       REAL             X(*)
C > >       DOUBLE PRECISION             Y(*)
C > >       DOUBLE PRECISION             DDOT
C > >       EXTERNAL             DDOT
C > >       INTEGER N_MAX
C > >       PARAMETER(NMAX=1000)
C > >       DOUBLE PRECISION             X_DOUBLE(NMAX)
C > >       IF ( N .LE. 0 ) THEN
C > >          SDDOT = 0.D0
C > >       ELSEIF ( N .EQ. 1 ) THEN
C > >          SDDOT = X(1)*Y(1)
C > >       ELSEIF ( N .LE. NMAX ) THEN
C > >          DO I=1, N
C > >             X_DOUBLE(I) = X(1 + (I-1)*INCX)
C > >          ENDDO
C > >          SDDOT = DDOT(N,X_DOUBLE,1,Y,INCY)
C > >       ELSE
C > >          PRINT *, N, ' too large in function sddot'
C > >          STOP
C > >       ENDIF
C > >       RETURN
C > >       END
      DOUBLE PRECISION FUNCTION MSDDOT(N,X,INCX,Y,INCY)
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
      REAL X(*)
      DOUBLE PRECISION Y(*)
      MSDDOT = 0.0D0
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
         MSDDOT = MSDDOT + X(IX)*Y(IY)
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
         MSDDOT = MSDDOT + X(I)*Y(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
         MSDDOT = MSDDOT + X(I)*Y(I) + X(I+1)*Y(I+1) +
     $   X(I + 2)*Y(I + 2) + X(I + 3)*Y(I + 3) + X(I + 4)*Y(I + 4)
   50 CONTINUE
      RETURN
C
C         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          MSDDOT = MSDDOT + X(I)*Y(I)
   70     CONTINUE
      RETURN
      END
      DOUBLE PRECISION      FUNCTION MDSDOT(N,X,INCX,Y,INCY)
C Dot product of mixed type input
      INTEGER N, INCX, INCY, IX, IY, I, M, MP1, NS
      DOUBLE PRECISION             X(*)
      REAL             Y(*)
      DOUBLE PRECISION             MSDDOT
      EXTERNAL             MSDDOT
         MDSDOT = MSDDOT(N,Y,INCY,X,INCX)
      RETURN
      END

C $Log$
C Revision 1.2  2002/09/13 22:49:10  pwagner
C Change external names to MSDDOT and MDSDOT
C
C Revision 1.1  2002/09/13 18:02:07  pwagner
C First commit
C
