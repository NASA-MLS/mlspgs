! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE MLSL1Utils              ! Common utilities for the MLSL1 program
!=============================================================================

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSFile: $"
  !---------------------------------------------------------------------------

  ! This module contains utility routines for the MLSL1 program

CONTAINS

  FUNCTION BigEndianStr (str) RESULT (number)

    !! Extracts BIG ENDIAN string into integer (MSB first)
    !! Note: will use only a maximum of the first 4 characters

    CHARACTER (LEN=*), INTENT(IN) :: str

    INTEGER :: number

    INTEGER, PARAMETER :: factor = 256
    INTEGER :: i, slen

    slen = MIN (LEN (str), 4)   ! allow only 4 chars MAX

    number = 0
    DO i = 1, slen
       number = number + factor**(slen-i) * ICHAR (str(i:i))
    ENDDO

  END FUNCTION BigEndianStr

  SUBROUTINE ExtractBigEndians (str, ibige)

    !! Returns array of integers converted from Big Endian string array

    CHARACTER (LEN=*), DIMENSION (:), INTENT(IN) :: str
    INTEGER, DIMENSION(:), POINTER  :: ibige

    INTEGER :: i

    IF (.NOT. ASSOCIATED (ibige)) THEN
       ALLOCATE (ibige(SIZE(str)))     ! Allocate the needed space
    ENDIF

    !! Consider reallocating receiving array to be the exact size?

    DO i = 1, (MIN (SIZE(ibige), SIZE(str)))
       ibige(i) = BigEndianStr (str(i))
    ENDDO

  END SUBROUTINE ExtractBigEndians

  SUBROUTINE SwapBytes (source, dest)

    !! Swaps consecutive bytes in a character string

    CHARACTER (len=*), INTENT (IN) :: source
    CHARACTER (len=*), INTENT (OUT) :: dest

    CHARACTER (len=1) :: temp
    INTEGER :: i, slen

    slen = MIN (LEN(source), LEN(dest))
    IF (slen == 1) THEN
       dest(1:1) = source(1:1)  !! What else can be done?!
       RETURN
    ENDIF
    DO i = 1, (slen-1), 2
       temp = source(i:i)   ! in case dest = source
       dest(i:i) = source(i+1:i+1)
       dest(i+1:i+1) = temp
    ENDDO
    dest(slen:slen) = source(slen:slen)  ! Last character (for odd slen)

  END SUBROUTINE SwapBytes

  FUNCTION QNan () RESULT (xnan)

    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_value, ieee_quiet_nan

    !! Return a quiet NaN for marking missing/undefined data

    REAL :: xnan

    xnan = ieee_value (0.0, ieee_quiet_nan)

  END FUNCTION QNan

  FUNCTION Finite (x) RESULT (is_finite)

    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite

    !! Return whether or not input number is finite

    REAL :: x
    LOGICAL :: is_finite

    is_finite = ieee_is_finite (x)

  END FUNCTION Finite

  FUNCTION GetIndexedAvg (value, indx) RESULT (avg)

    !! Return average of indexed value array

    REAL, DIMENSION(:) :: value
    INTEGER, DIMENSION (:) :: indx  ! Indexes for "value"

    REAL :: avg

    INTEGER :: i, navg, low, high
    REAL :: sum

    low = LBOUND (value, 1)
    high = UBOUND (value, 1)

    navg = 0
    sum = 0.0
    DO i = 1, SIZE (indx)
       IF (indx(i) >= low .AND. indx(i) <= high) THEN
          IF (Finite(value(indx(i)))) THEN
             navg = navg + 1
             sum = sum + value(indx(i))
          ENDIF
       ENDIF
    ENDDO

    IF (navg > 0) THEN
       avg = sum / navg
    ELSE
       avg = Qnan()
    ENDIF

  END FUNCTION GetIndexedAvg

!=============================================================================
END MODULE MLSL1Utils
!=============================================================================

! $Log$
! Revision 2.2  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.1  2001/02/23 20:52:08  perun
! Version 0.5 commit
!
