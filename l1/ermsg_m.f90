MODULE ERMSG_M
!     .  Copyright (C) 1989-1999, California Institute of Technology.
!     .  U. S. Government sponsorship under
!     .  NASA contract NAS7-918 is acknowledged.

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  PRIVATE

  PUBLIC ERFIN, ERM1, IERM1, SERM1, DERM1
  INTERFACE ERM1; MODULE PROCEDURE IERM1, SERM1, DERM1; END INTERFACE
  PUBLIC ERMN, IERMN, SERMN, DERMN
  INTERFACE ERMN; MODULE PROCEDURE IERMN, SERMN, DERMN; END INTERFACE
  PUBLIC ERMOR, ERMSET, ERMSG, ERV1, IERV1, SERV1, DERV1
  INTERFACE ERV1; MODULE PROCEDURE IERV1, SERV1, DERV1; END INTERFACE
  PUBLIC ERVN, IERVN, SERVN, DERVN
  INTERFACE ERVN; MODULE PROCEDURE IERVN, SERVN, DERVN; END INTERFACE

  INTEGER, SAVE, PUBLIC :: IDELTA=0 ! increment to LEVEL to calculate
                                    ! IALPHA in ERMSG, q.v.
  INTEGER, SAVE, PUBLIC :: OLDLEV=0 ! last value of LEVEL passed to ERMSG

  INTEGER, SAVE, PRIVATE :: IALPHA  ! = IDELTA + LEVEL; used to determine
                                    ! whether to print

CONTAINS
  SUBROUTINE ERFIN
!>> 1998-04-13 ERFIN  Snyder  Convert to F90
!>> 1994-11-11 ERFIN  CLL     Typing all variables.
!>> 1985-09-23 ERFIN  Lawson  Initial code.
 
    WRITE (*,"(1X,72('$')/' ')")
    IF (ialpha >= 2) STOP
    RETURN
  END SUBROUTINE ERFIN
!=======================================================================
  SUBROUTINE IERM1 (SUBNAM, INDIC, LEVEL, MSG, LABEL, VALUE, FLAG)
!>> 1998-04-13 IERM1  Snyder  Convert to F90
!>> 1990-01-18 IERM1  CLL     Typed all variables.
!>> 1985-08-02 IERM1  Lawson  Initial code.
    INTEGER, INTENT(in) :: INDIC, LEVEL, VALUE
    CHARACTER(LEN=*), INTENT(in) :: SUBNAM, MSG, LABEL
    CHARACTER(LEN=1), INTENT(in) :: FLAG
    CALL ERMSG (SUBNAM, INDIC, LEVEL, MSG, ',')
    CALL ERV1 (LABEL, VALUE, FLAG)
    RETURN
  END SUBROUTINE IERM1

  SUBROUTINE SERM1 (SUBNAM, INDIC, LEVEL, MSG, LABEL, VALUE, FLAG)
!>> 1998-04-13 SERM1  Snyder  Convert to F90
!>> 1990-01-18 SERM1  CLL     Typed all variables.
!>> 1985-08-02 SERM1  Lawson  Initial code.
    INTEGER, INTENT(in) :: INDIC, LEVEL
    REAL, INTENT(IN) :: VALUE
    CHARACTER(LEN=*), INTENT(in) :: SUBNAM, MSG, LABEL
    CHARACTER(LEN=1), INTENT(in) :: FLAG
    CALL ERMSG (SUBNAM, INDIC, LEVEL, MSG, ',')
    CALL ERV1 (LABEL, VALUE, FLAG)
    RETURN
  END SUBROUTINE SERM1

  SUBROUTINE DERM1 (SUBNAM, INDIC, LEVEL, MSG, LABEL, VALUE, FLAG)
!>> 1998-04-13 DERM1  Snyder  Convert to F90
!>> 1990-01-18 DERM1  CLL     Typed all variables.
!>> 1985-08-02 DERM1  Lawson  Initial code.
    INTEGER, INTENT(in) :: INDIC, LEVEL
    DOUBLE PRECISION, INTENT(IN) :: VALUE
    CHARACTER(LEN=*), INTENT(in) :: SUBNAM, MSG, LABEL
    CHARACTER(LEN=1), INTENT(in) :: FLAG
    CALL ERMSG (SUBNAM, INDIC, LEVEL, MSG, ',')
    CALL ERV1 (LABEL, VALUE, FLAG)
    RETURN
  END SUBROUTINE DERM1
!=======================================================================

! These subroutines call ERMSG to initiate an error message and then
! call ERVN to print an array of values with labels.
!
! ------------------------------------------------------------------
! SUBROUTINE ARGUMENTS
! --------------------
! SUBNAM   A name that identifies the subprogram in which
!          the error occurs.
!
! INDIC    An integer printed as part of the mininal error
!          message. It together with SUBNAM can be used to
!          uniquely identify an error.
!
! LEVEL    The user sets LEVEL=2,0,or -2 to specify the
!          nominal action to be taken.
!          See subr ERMSG for interpretation of LEVEL.
!
! MSG      Message to be printed to describe the error.
!
! LABELS() An array of character data to identify the individual
!          scalars.  There must be one character label for each scal
!          The character string length is used to determine how many
!          labelled scalars can be output on a single line (of 132).
!
! VALUES() An array of scalars to be printed with their identifiers.
!
! FLAG     A single character, which when set to '.' will
!          call the subroutine ERFIN and will just RETURN
!          when set to any other character.
!
! ------------------------------------------------------------------
!
! Kris Stewart, JPL, 1983 July.
! C.Lawson & S.Chan, JPL, 1983,Nov.
!
! ------------------------------------------------------------------

  SUBROUTINE IERMN (SUBNAM, INDIC, LEVEL, MSG, LABELS, VALUES, FLAG)
!>> 1998-04-13 IERMN  Snyder  Convert to F90
!>> 1994-10-20 IERMN  Krogh   Changes to use M77CON
!>> 1987-10-02 IERMN  Lawson  Initial code.
    INTEGER, INTENT(in) :: LEVEL, INDIC
    CHARACTER(len=*), INTENT(in) :: SUBNAM, MSG, LABELS(:)
    CHARACTER(len=1), INTENT(in) :: FLAG
    INTEGER, INTENT(in) :: VALUES(:)

    CALL ermsg (subnam, indic, level,msg, ',')
    CALL ervn (labels, values, flag)
    RETURN
  END SUBROUTINE IERMN

  SUBROUTINE SERMN (SUBNAM, INDIC, LEVEL, MSG, LABELS, VALUES, FLAG)
!>> 1998-04-13 SERMN  Snyder  Convert to F90
!>> 1994-10-20 SERMN  Krogh   Changes to use M77CON
!>> 1987-10-02 SERMN  Lawson  Initial code.
    INTEGER, INTENT(in) :: LEVEL, INDIC
    CHARACTER(len=*), INTENT(in) :: SUBNAM, MSG, LABELS(:)
    CHARACTER(len=1), INTENT(in) :: FLAG
    REAL, INTENT(in) :: VALUES(:)

    CALL ermsg (subnam, indic, level,msg, ',')
    CALL ervn (labels, values, flag)
    RETURN
  END SUBROUTINE SERMN

  SUBROUTINE DERMN (SUBNAM, INDIC, LEVEL, MSG, LABELS, VALUES, FLAG)
!>> 1998-04-13 DERMN  Snyder  Convert to F90
!>> 1994-10-20 DERMN  Krogh   Changes to use M77CON
!>> 1987-10-02 DERMN  Lawson  Initial code.
    INTEGER, INTENT(in) :: LEVEL, INDIC
    CHARACTER(len=*), INTENT(in) :: SUBNAM, MSG, LABELS(:)
    CHARACTER(len=1), INTENT(in) :: FLAG
    DOUBLE PRECISION, INTENT(in) :: VALUES(:)

    CALL ermsg (subnam, indic, level,msg, ',')
    CALL ervn (labels, values, flag)
    RETURN
  END SUBROUTINE DERMN
!=======================================================================
  SUBROUTINE ERMOR (MSG, FLAG)
!>> 1998-04-13 ERMOR  WVSnyder  Convert to F90
!>> 1985-09-20 ERMOR  Lawson  Initial code.

! --------------------------------------------------------------
! SUBROUTINE ARGUMENTS
! --------------------
! MSG      Message to be printed as part of the diagnostic.
!
! FLAG     A single character,which when set to '.' will
!          call the subroutine ERFIN and will just RETURN
!          when set to any other character.
!
! --------------------------------------------------------------
!
    CHARACTER(len=*), INTENT(in) :: MSG
    CHARACTER(len=1), INTENT(in) :: FLAG

    IF (ialpha >= -1) THEN
      WRITE (*,*) msg
      IF (flag == '.') CALL erfin
    END IF
    RETURN
  END SUBROUTINE ERMOR
!=======================================================================
! ERMSET resets IDELTA.
  SUBROUTINE ERMSET (IDEL)
    INTEGER, INTENT(in) :: IDEL
    idelta = idel
    RETURN
  END SUBROUTINE ERMSET
!=======================================================================
  SUBROUTINE ERMSG (SUBNAM, INDIC, LEVEL, MSG, FLAG)
!>> 1998-04-13 ERMSG  WVSnyder  Convert to F90
!>> 1994-11-11 ERMSG  Krogh     Declared all vars.
!>> 1992-10-20 ERMSG  WVSnyder  added ERLSET, ERLGET
!>> 1985-09-25 ERMSG  Lawson    Initial code.

!     --------------------------------------------------------------

! ERMSG initiates an error message, and manages the saved values IDELTA
! and IALPHA to control the level of action. This is intended to
! be the only subroutine that assigns a value to IALPHA.

!
!     --------------------------------------------------------------
!     SUBROUTINE ARGUMENTS
!     --------------------
!     SUBNAM   A name that identifies the subprogram in which
!              the error occurs.
!
!     INDIC    An integer printed as part of the mininal error
!              message. It together with SUBNAM can be used to
!              uniquely identify an error.
!
!     LEVEL    The user sets LEVEL=2,0,or -2 to specify the
!              nominal action to be taken by ERMSG. The
!              module ERMSG_M contains a private variable
!              IDELTA, whose nominal value is zero. The
!              subroutine will compute IALPHA = LEVEL + IDELTA
!              and proceed as follows:
!              If (IALPHA >= 2)        Print message and STOP.
!              If (IALPHA=-1,0,1)      Print message and return.
!              If (IALPHA <= -2)       Just RETURN.
!
!     MSG      Message to be printed as part of the diagnostic.
!
!     FLAG     A single character,which when set to '.' will
!              call the subroutine ERFIN and will just RETURN
!              when set to any other character.
!
!     --------------------------------------------------------------
!
!     C.Lawson & S.Chan, JPL, 1983 Nov
!
!     ------------------------------------------------------------------
    INTEGER, INTENT(in) :: LEVEL, INDIC
    CHARACTER(len=*), INTENT(in) :: SUBNAM, MSG
    CHARACTER(len=1), INTENT(in) :: FLAG

    oldlev = level
    ialpha = level + idelta
    IF (ialpha >= -1) THEN
      WRITE (*,"('0',72('$')/' SUBPROGRAM ',A,' REPORTS ERROR NO. ',I4)") &
        subnam, indic
      WRITE (*,*) msg
      IF (flag == '.') CALL erfin
    ENDIF
    RETURN
  END SUBROUTINE ERMSG
!=======================================================================
! ------------------------------------------------------------
! SUBROUTINE ARGUMENTS
! --------------------
! LABEL     An identifing name to be printed with VALUE.
!
! VALUE     A integer to be printed.
!
! FLAG      See write up for FLAG in ERMSG.
!
! ------------------------------------------------------------

  SUBROUTINE IERV1 (LABEL, VALUE, FLAG)
!>> 1998-04-13 IERV1  Snyder  Convert to F90
!>> 1985-09-20 IERV1  Lawson  Initial code.
    INTEGER, INTENT(in) :: VALUE
    CHARACTER(len=*), INTENT(in) ::  LABEL
    CHARACTER(len=1), INTENT(in) ::  FLAG

    IF (ialpha >= -1) THEN
      WRITE (*,*) label, ' = ', value
      IF (flag == '.') CALL erfin
    END IF
    RETURN
  END SUBROUTINE IERV1

  SUBROUTINE SERV1 (LABEL, VALUE, FLAG)
!>> 1998-04-13 SERV1  Snyder  Convert to F90
!>> 1985-09-20 SERV1  Lawson  Initial code.
    REAL, INTENT(in) :: VALUE
    CHARACTER(len=*), INTENT(in) ::  LABEL
    CHARACTER(len=1), INTENT(in) ::  FLAG

    IF (ialpha >= -1) THEN
      WRITE (*,*) label, ' = ', value
      IF (flag == '.') CALL erfin
    END IF
    RETURN
  END SUBROUTINE SERV1

  SUBROUTINE DERV1 (LABEL, VALUE, FLAG)
!>> 1998-04-13 DERV1  Snyder  Convert to F90
!>> 1985-09-20 DERV1  Lawson  Initial code.
    DOUBLE PRECISION, INTENT(in) :: VALUE
    CHARACTER(len=*), INTENT(in) ::  LABEL
    CHARACTER(len=1), INTENT(in) ::  FLAG

    IF (ialpha >= -1) THEN
      WRITE (*,*) label, ' = ', value
      IF (flag == '.') CALL erfin
    END IF
    RETURN
  END SUBROUTINE DERV1
!=======================================================================

! These subroutines print an array of values with labels as part of an
! error msg. Uses IALPHA that must have been previously set by ERMSG.

!     ------------------------------------------------------------------
!
! LABELS() An array of character data to identify the individual
!          scalars.  There must be one character label for each scal
!          The character string length is used to determine how many
!          labelled scalars can be output on a single line (of 132).
!
! VALUES() An array of single precision scalars to be printed with
!          their identifiers.
!
! FLAG     A single character, which when set to '.' will
!          call the subroutine ERFIN and will just RETURN
!          when set to any other character.

  SUBROUTINE IERVN (LABELS, VALUES, FLAG)
!>> 1998-04-13 IERVN  Snyder Convert to F90
!>> 1994-10-20 IERVN  Krogh  Changes to use M77CON
!>> 1989-11-10 IERVN  CLL
!>> 1987-10-02 IERVN  Lawson  Initial code.
! ----------------------------------------------------------------------
! 1989-11-10 CLL  Changed to keep length of printed line not
! greater than MAXCOL = 75 characters.
! ----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: LABELS(:)
    INTEGER, INTENT(in) :: VALUES(SIZE(LABELS))
    CHARACTER(LEN=1), INTENT(in) :: FLAG

    INTEGER I, K, LENIDV, NUMBER
    INTEGER, PARAMETER :: MAXCOL = 75

    IF (ialpha >= -1) THEN
      lenidv = LEN (labels)
      number = maxcol / (lenidv+17)
      DO i = 1, SIZE(labels), number
        WRITE(*,"(4(2x,a,'=',i14))") &
          (labels(k), values(k), k = i, MIN(SIZE(labels), i+number-1) )
      END DO
      IF (flag == '.') CALL erfin
    ENDIF
    RETURN
  END SUBROUTINE IERVN

  SUBROUTINE SERVN (LABELS, VALUES, FLAG)
!>> 1998-04-13 SERVN  Snyder Convert to F90
!>> 1994-10-20 SERVN  Krogh  Changes to use M77CON
!>> 1989-11-10 SERVN  CLL
!>> 1987-10-02 SERVN  Lawson  Initial code.
! ----------------------------------------------------------------------
! 1989-11-10 CLL  Changed to keep length of printed line not
! greater than MAXCOL = 75 characters.
! ----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: LABELS(:)
    REAL, INTENT(in) :: VALUES(SIZE(LABELS))
    CHARACTER(LEN=1), INTENT(in) :: FLAG

    INTEGER I, K, LENIDV, NUMBER
    INTEGER, PARAMETER :: MAXCOL = 75

    IF (ialpha >= -1) THEN
      lenidv = LEN (labels)
      number = maxcol / (lenidv+17)
      DO i = 1, SIZE(labels), number
        WRITE(*,"(4(2x,a,'=',g14.7))") &
          (labels(k), values(k), k = i, MIN(SIZE(labels), i+number-1) )
      END DO
      IF (flag == '.') CALL erfin
    ENDIF
    RETURN
  END SUBROUTINE SERVN

  SUBROUTINE DERVN (LABELS, VALUES, FLAG)
!>> 1998-04-13 DERVN  Snyder Convert to F90
!>> 1994-10-20 DERVN  Krogh  Changes to use M77CON
!>> 1989-11-10 DERVN  CLL
!>> 1987-10-02 DERVN  Lawson  Initial code.
! ----------------------------------------------------------------------
! 1989-11-10 CLL  Changed to keep length of printed line not
! greater than MAXCOL = 75 characters.
! ----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: LABELS(:)
    DOUBLE PRECISION, INTENT(in) :: VALUES(SIZE(LABELS))
    CHARACTER(LEN=1), INTENT(in) :: FLAG

    INTEGER I, K, LENIDV, NUMBER
    INTEGER, PARAMETER :: MAXCOL = 75

    IF (ialpha >= -1) THEN
      lenidv = LEN (labels)
      number = maxcol / (lenidv+17)
      DO i = 1, SIZE(labels), number
        WRITE(*,"(4(2x,a,'=',g14.7))") &
          (labels(k), values(k), k = i, MIN(SIZE(labels), i+number-1) )
      END DO
      IF (flag == '.') CALL erfin
    ENDIF
    RETURN
  END SUBROUTINE DERVN
END MODULE ERMSG_M

! $Log$
! Revision 2.1  2002/03/29 20:21:45  perun
! Version 1.0 commit
!
