! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
MODULE MLSL1Utils              ! Common utilities for the MLSL1 program
!=============================================================================

  use MLSFillValues, only: isNaN
  implicit none

  PRIVATE

  PUBLIC :: BigEndianStr, ExtractBigEndians, SwapBytes, QNan, Finite, &
       GetIndexedAvg

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  INTERFACE Finite
     MODULE PROCEDURE Finite_S, Finite_D
  END INTERFACE

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

    number = ICHAR (str(1:1))
    DO i = 2, slen
       number = number * factor + ICHAR (str(i:i))
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

    USE ieee_arithmetic, ONLY: ieee_value, ieee_quiet_nan

    !! Return a quiet NaN for marking missing/undefined data

    REAL :: xnan

    xnan = ieee_value (0.0, ieee_quiet_nan)

  END FUNCTION QNan

  FUNCTION Finite_S (s) RESULT (is_finite)

    USE ieee_arithmetic, ONLY: ieee_is_finite

    !! Return whether or not input number is finite

    REAL :: s
    LOGICAL :: is_finite

    is_finite = ieee_is_finite (s)

  END FUNCTION Finite_S

  FUNCTION Finite_D (d) RESULT (is_finite)

    USE ieee_arithmetic, ONLY: ieee_is_finite

    !! Return whether or not input number is finite

    DOUBLE PRECISION :: d
    LOGICAL :: is_finite

    is_finite = ieee_is_finite (d)

  END FUNCTION Finite_D

  FUNCTION GetIndexedAvg (value, indx, debug) RESULT (avg)

    !! Return average of indexed value array

    REAL, DIMENSION(:) :: value
    INTEGER, DIMENSION (:) :: indx  ! Indexes for "value"
    logical, optional, intent(in) :: debug

    REAL :: avg

    INTEGER :: i, navg, low, high
    REAL :: sum
    logical :: myDebug

    low = LBOUND (value, 1)
    high = UBOUND (value, 1)
    myDebug = .false.
    if ( present(debug) ) myDebug = debug

    navg = 0
    sum = 0.0
    DO i = 1, SIZE (indx)
       IF (indx(i) >= low .AND. indx(i) <= high) THEN
          IF (Finite(value(indx(i)))) THEN
             navg = navg + 1
             sum = sum + value(indx(i))
             if ( myDebug .and. isNaN(value(indx(i))) ) then
               print *, 'i, indx(i), value(indx(i)) ', i, indx(i), value(indx(i))
             endif
          elseif ( myDebug ) then
             print *, 'i, indx(i), value(indx(i)) ', i, indx(i), value(indx(i))
          ENDIF
       elseif ( myDebug ) then
           print *, 'i, indx(i), low, high ', i, indx(i), low, high
       ENDIF
    ENDDO

    IF (navg > 0) THEN
       avg = sum / navg
    ELSE
       avg = Qnan()
       if ( myDebug ) print *, 'navg = 0, so returning NaN'
    ENDIF

  END FUNCTION GetIndexedAvg

!=============================================================================
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here
END MODULE MLSL1Utils
!=============================================================================

! $Log$
! Revision 2.8  2015/01/21 19:30:38  pwagner
! Gets isNaN from MLSFillValues
!
! Revision 2.7  2015/01/13 18:40:43  pwagner
! Can warn if producing a NaN
!
! Revision 2.6  2005/06/23 18:41:36  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.5  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
! Revision 2.4  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
! Revision 2.2  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.1  2001/02/23 20:52:08  perun
! Version 0.5 commit
!
