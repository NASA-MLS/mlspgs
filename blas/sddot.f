C Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
C U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

      DOUBLE PRECISION      FUNCTION SDDOT(N,X,INCX,Y,INCY)
C Dot product of mixed type input
      INTEGER N, INCX, INCY, IX, IY, I, M, MP1, NS
      REAL             X(*)
      DOUBLE PRECISION             Y(*)
      DOUBLE PRECISION             DDOT
      EXTERNAL             DDOT
      INTEGER N_MAX
      PARAMETER(NMAX=1000)
      DOUBLE PRECISION             X_DOUBLE(NMAX)
      IF ( N .LE. 0 ) THEN
         SDDOT = 0.D0
      ELSEIF ( N .EQ. 1 ) THEN
         SDDOT = X(1)*Y(1)
      ELSEIF ( N .LE. NMAX ) THEN
         DO I=1, N
            X_DOUBLE(I) = X(1 + (I-1)*INCX)
         ENDDO
         SDDOT = DDOT(N,X_DOUBLE,1,Y,INCY)
      ELSE
         PRINT *, N, ' too large in function sddot'
         STOP
      ENDIF
      RETURN
      END
      DOUBLE PRECISION      FUNCTION DSDOT(N,X,INCX,Y,INCY)
C Dot product of mixed type input
      INTEGER N, INCX, INCY, IX, IY, I, M, MP1, NS
      DOUBLE PRECISION             X(*)
      REAL             Y(*)
      DOUBLE PRECISION             SDDOT
      EXTERNAL             SDDOT
         DSDOT = SDDOT(N,Y,INCY,X,INCX)
      RETURN
      END

C $Log$
