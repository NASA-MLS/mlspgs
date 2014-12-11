! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

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

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

CONTAINS ! ============= Public procedures ===================================

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

  pure FUNCTION SwapShort (n2)

    use MLSKinds, only: I2

    integer(i2), intent(in) :: n2
    integer(i2) :: SwapShort

    CHARACTER(len=1) :: cbuf(2)

    cbuf = TRANSFER (n2, cbuf)
    SwapShort = TRANSFER (cbuf(2:1:-1), n2)

  END FUNCTION SwapShort

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

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

END MODULE SwapEndian

! $Log$
