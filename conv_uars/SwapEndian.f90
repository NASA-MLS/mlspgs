MODULE SwapEndian

! Contains functions to swap endianess of input data into the endianess of
!  the local machine

  IMPLICIT NONE

  INTERFACE SwapBytes
     MODULE PROCEDURE SwapReal, SwapInteger
  END INTERFACE

  INTERFACE SwapBig
     MODULE PROCEDURE SwapBigReal, SwapBigInteger
  END INTERFACE

  INTERFACE SwapLittle
     MODULE PROCEDURE SwapLittleReal, SwapLittleInteger
  END INTERFACE

  CONTAINS

    REAL FUNCTION SwapReal (r4)

      REAL :: r4
      CHARACTER(len=1) :: cbuf(4)

      cbuf = TRANSFER (r4, cbuf)
      SwapReal = TRANSFER (cbuf(4:1:-1), r4)

    END FUNCTION SwapReal

    INTEGER FUNCTION SwapInteger (i4)

      INTEGER :: i4
      CHARACTER(len=1) :: cbuf(4)

      cbuf = TRANSFER (i4, cbuf)
      SwapInteger = TRANSFER (cbuf(4:1:-1), i4)

    END FUNCTION SwapInteger

    REAL FUNCTION SwapBigReal (r4)

      REAL :: r4
      CHARACTER(LEN=1) :: cbuf

      IF (ICHAR(TRANSFER (1, cbuf)) == 1) THEN   ! machine is little endian
         SwapBigReal = SwapBytes (r4)
      ELSE
         SwapBigReal = r4
      ENDIF

    END FUNCTION SwapBigReal

    REAL FUNCTION SwapBigInteger (i4)

      integer :: i4
      CHARACTER(LEN=1) :: cbuf

      IF (ICHAR(TRANSFER (1, cbuf)) == 1) THEN   ! machine is little endian
         SwapBigInteger = SwapBytes (i4)
      ELSE
         SwapBigInteger = i4
      ENDIF

    END FUNCTION SwapBigInteger

    REAL FUNCTION SwapLittleReal (r4)

      REAL :: r4
      CHARACTER(LEN=1) :: cbuf

      IF (ICHAR(TRANSFER (1, cbuf)) == 0) THEN   ! machine is big endian
         SwapLittleReal = SwapBytes (r4)
      ELSE
         SwapLittleReal = r4
      ENDIF

    END FUNCTION SwapLittleReal

    REAL FUNCTION SwapLittleInteger (i4)

      integer :: i4
      CHARACTER(LEN=1) :: cbuf

      IF (ICHAR(TRANSFER (1, cbuf)) == 0) THEN   ! machine is big endian
         SwapLittleInteger = SwapBytes (i4)
      ELSE
         SwapLittleInteger = i4
      ENDIF

    END FUNCTION SwapLittleInteger

END MODULE SwapEndian
