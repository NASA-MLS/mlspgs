module ERMSG_M
!     .  Copyright (C) 1989-1999, California Institute of Technology.
!     .  U. S. Government sponsorship under
!     .  NASA contract NAS7-918 is acknowledged.
  implicit NONE
  private
  public ERFIN, ERM1
  interface ERM1; module procedure IERM1, SERM1, DERM1; end interface
  public ERMN
  interface ERMN; module procedure IERMN, SERMN, DERMN; end interface
  public ERMOR, ERMSET, ERMSG, ERV1
  interface ERV1; module procedure IERV1, SERV1, DERV1; end interface
  public ERVN
  interface ERVN; module procedure IERVN, SERVN, DERVN; end interface
  integer, save, public :: IDELTA=0 ! increment to LEVEL to calculate
                                    ! IALPHA in ERMSG, q.v.
  integer, save, public :: OLDLEV=0 ! last value of LEVEL passed to ERMSG
  integer, save, private :: IALPHA  ! = IDELTA + LEVEL; used to determine
                                    ! whether to print
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
  subroutine ERFIN
!>> 1998-04-13 ERFIN  Snyder  Convert to F90
!>> 1994-11-11 ERFIN  CLL     Typing all variables.
!>> 1985-09-23 ERFIN  Lawson  Initial code.
    write (*,"(1X,72('$')/' ')")
    if (ialpha >= 2) stop
    return
  end subroutine ERFIN
!=======================================================================
  subroutine IERM1 (SUBNAM, INDIC, LEVEL, MSG, LABEL, VALUE, FLAG)
!>> 1998-04-13 IERM1  Snyder  Convert to F90
!>> 1990-01-18 IERM1  CLL     Typed all variables.
!>> 1985-08-02 IERM1  Lawson  Initial code.
    integer, intent(in) :: INDIC, LEVEL, VALUE
    character(LEN=*), intent(in) :: SUBNAM, MSG, LABEL
    character(LEN=1), intent(in) :: FLAG
    call ERMSG (SUBNAM, INDIC, LEVEL, MSG, ',')
    call ERV1 (LABEL, VALUE, FLAG)
    return
  end subroutine IERM1
  subroutine SERM1 (SUBNAM, INDIC, LEVEL, MSG, LABEL, VALUE, FLAG)
!>> 1998-04-13 SERM1  Snyder  Convert to F90
!>> 1990-01-18 SERM1  CLL     Typed all variables.
!>> 1985-08-02 SERM1  Lawson  Initial code.
    integer, intent(in) :: INDIC, LEVEL
    real, intent(IN) :: VALUE
    character(LEN=*), intent(in) :: SUBNAM, MSG, LABEL
    character(LEN=1), intent(in) :: FLAG
    call ERMSG (SUBNAM, INDIC, LEVEL, MSG, ',')
    call ERV1 (LABEL, VALUE, FLAG)
    return
  end subroutine SERM1
  subroutine DERM1 (SUBNAM, INDIC, LEVEL, MSG, LABEL, VALUE, FLAG)
!>> 1998-04-13 DERM1  Snyder  Convert to F90
!>> 1990-01-18 DERM1  CLL     Typed all variables.
!>> 1985-08-02 DERM1  Lawson  Initial code.
    integer, intent(in) :: INDIC, LEVEL
    double precision, intent(IN) :: VALUE
    character(LEN=*), intent(in) :: SUBNAM, MSG, LABEL
    character(LEN=1), intent(in) :: FLAG
    call ERMSG (SUBNAM, INDIC, LEVEL, MSG, ',')
    call ERV1 (LABEL, VALUE, FLAG)
    return
  end subroutine DERM1
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
  subroutine IERMN (SUBNAM, INDIC, LEVEL, MSG, LABELS, VALUES, FLAG)
!>> 1998-04-13 IERMN  Snyder  Convert to F90
!>> 1994-10-20 IERMN  Krogh   Changes to use M77CON
!>> 1987-10-02 IERMN  Lawson  Initial code.
    integer, intent(in) :: LEVEL, INDIC
    character(len=*), intent(in) :: SUBNAM, MSG, LABELS(:)
    character(len=1), intent(in) :: FLAG
    integer, intent(in) :: VALUES(:)
    call ermsg (subnam, indic, level,msg, ',')
    call ervn (labels, values, flag)
    return
  end subroutine IERMN
  subroutine SERMN (SUBNAM, INDIC, LEVEL, MSG, LABELS, VALUES, FLAG)
!>> 1998-04-13 SERMN  Snyder  Convert to F90
!>> 1994-10-20 SERMN  Krogh   Changes to use M77CON
!>> 1987-10-02 SERMN  Lawson  Initial code.
    integer, intent(in) :: LEVEL, INDIC
    character(len=*), intent(in) :: SUBNAM, MSG, LABELS(:)
    character(len=1), intent(in) :: FLAG
    real, intent(in) :: VALUES(:)
    call ermsg (subnam, indic, level,msg, ',')
    call ervn (labels, values, flag)
    return
  end subroutine SERMN
  subroutine DERMN (SUBNAM, INDIC, LEVEL, MSG, LABELS, VALUES, FLAG)
!>> 1998-04-13 DERMN  Snyder  Convert to F90
!>> 1994-10-20 DERMN  Krogh   Changes to use M77CON
!>> 1987-10-02 DERMN  Lawson  Initial code.
    integer, intent(in) :: LEVEL, INDIC
    character(len=*), intent(in) :: SUBNAM, MSG, LABELS(:)
    character(len=1), intent(in) :: FLAG
    double precision, intent(in) :: VALUES(:)
    call ermsg (subnam, indic, level,msg, ',')
    call ervn (labels, values, flag)
    return
  end subroutine DERMN
!=======================================================================
  subroutine ERMOR (MSG, FLAG)
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
    character(len=*), intent(in) :: MSG
    character(len=1), intent(in) :: FLAG
    if (ialpha >= -1) then
      write (*,*) msg
      if (flag == '.') call erfin
    end if
    return
  end subroutine ERMOR
!=======================================================================
! ERMSET resets IDELTA.
  subroutine ERMSET (IDEL)
    integer, intent(in) :: IDEL
    idelta = idel
    return
  end subroutine ERMSET
!=======================================================================
  subroutine ERMSG (SUBNAM, INDIC, LEVEL, MSG, FLAG)
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
    integer, intent(in) :: LEVEL, INDIC
    character(len=*), intent(in) :: SUBNAM, MSG
    character(len=1), intent(in) :: FLAG
    oldlev = level
    ialpha = level + idelta
    if (ialpha >= -1) then
      write (*,"('0',72('$')/' SUBPROGRAM ',A,' REPORTS ERROR NO. ',I4)") &
        subnam, indic
      write (*,*) msg
      if (flag == '.') call erfin
    endif
    return
  end subroutine ERMSG
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
  subroutine IERV1 (LABEL, VALUE, FLAG)
!>> 1998-04-13 IERV1  Snyder  Convert to F90
!>> 1985-09-20 IERV1  Lawson  Initial code.
    integer, intent(in) :: VALUE
    character(len=*), intent(in) ::  LABEL
    character(len=1), intent(in) ::  FLAG
    if (ialpha >= -1) then
      write (*,*) label, ' = ', value
      if (flag == '.') call erfin
    end if
    return
  end subroutine IERV1
  subroutine SERV1 (LABEL, VALUE, FLAG)
!>> 1998-04-13 SERV1  Snyder  Convert to F90
!>> 1985-09-20 SERV1  Lawson  Initial code.
    real, intent(in) :: VALUE
    character(len=*), intent(in) ::  LABEL
    character(len=1), intent(in) ::  FLAG
    if (ialpha >= -1) then
      write (*,*) label, ' = ', value
      if (flag == '.') call erfin
    end if
    return
  end subroutine SERV1
  subroutine DERV1 (LABEL, VALUE, FLAG)
!>> 1998-04-13 DERV1  Snyder  Convert to F90
!>> 1985-09-20 DERV1  Lawson  Initial code.
    double precision, intent(in) :: VALUE
    character(len=*), intent(in) ::  LABEL
    character(len=1), intent(in) ::  FLAG
    if (ialpha >= -1) then
      write (*,*) label, ' = ', value
      if (flag == '.') call erfin
    end if
    return
  end subroutine DERV1
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
  subroutine IERVN (LABELS, VALUES, FLAG)
!>> 1998-04-13 IERVN  Snyder Convert to F90
!>> 1994-10-20 IERVN  Krogh  Changes to use M77CON
!>> 1989-11-10 IERVN  CLL
!>> 1987-10-02 IERVN  Lawson  Initial code.
! ----------------------------------------------------------------------
! 1989-11-10 CLL  Changed to keep length of printed line not
! greater than MAXCOL = 75 characters.
! ----------------------------------------------------------------------
    character(LEN=*), intent(in) :: LABELS(:)
    integer, intent(in) :: VALUES(SIZE(LABELS))
    character(LEN=1), intent(in) :: FLAG
    integer I, K, LENIDV, NUMBER
    integer, parameter :: MAXCOL = 75
    if (ialpha >= -1) then
      lenidv = len (labels)
      number = maxcol / (lenidv+17)
      do i = 1, size(labels), number
        write(*,"(4(2x,a,'=',i14))") &
          (labels(k), values(k), k = i, min(size(labels), i+number-1) )
      end do
      if (flag == '.') call erfin
    endif
    return
  end subroutine IERVN
  subroutine SERVN (LABELS, VALUES, FLAG)
!>> 1998-04-13 SERVN  Snyder Convert to F90
!>> 1994-10-20 SERVN  Krogh  Changes to use M77CON
!>> 1989-11-10 SERVN  CLL
!>> 1987-10-02 SERVN  Lawson  Initial code.
! ----------------------------------------------------------------------
! 1989-11-10 CLL  Changed to keep length of printed line not
! greater than MAXCOL = 75 characters.
! ----------------------------------------------------------------------
    character(LEN=*), intent(in) :: LABELS(:)
    real, intent(in) :: VALUES(SIZE(LABELS))
    character(LEN=1), intent(in) :: FLAG
    integer I, K, LENIDV, NUMBER
    integer, parameter :: MAXCOL = 75
    if (ialpha >= -1) then
      lenidv = len (labels)
      number = maxcol / (lenidv+17)
      do i = 1, size(labels), number
        write(*,"(4(2x,a,'=',g14.7))") &
          (labels(k), values(k), k = i, min(size(labels), i+number-1) )
      end do
      if (flag == '.') call erfin
    endif
    return
  end subroutine SERVN
  subroutine DERVN (LABELS, VALUES, FLAG)
!>> 1998-04-13 DERVN  Snyder Convert to F90
!>> 1994-10-20 DERVN  Krogh  Changes to use M77CON
!>> 1989-11-10 DERVN  CLL
!>> 1987-10-02 DERVN  Lawson  Initial code.
! ----------------------------------------------------------------------
! 1989-11-10 CLL  Changed to keep length of printed line not
! greater than MAXCOL = 75 characters.
! ----------------------------------------------------------------------
    character(LEN=*), intent(in) :: LABELS(:)
    double precision, intent(in) :: VALUES(SIZE(LABELS))
    character(LEN=1), intent(in) :: FLAG
    integer I, K, LENIDV, NUMBER
    integer, parameter :: MAXCOL = 75
    if (ialpha >= -1) then
      lenidv = len (labels)
      number = maxcol / (lenidv+17)
      do i = 1, size(labels), number
        write(*,"(4(2x,a,'=',g14.7))") &
          (labels(k), values(k), k = i, min(size(labels), i+number-1) )
      end do
      if (flag == '.') call erfin
    endif
    return
  end subroutine DERVN
end module ERMSG_M
! $Log$
! Revision 1.2  2000/05/04 23:39:24  vsnyder
! Initial conversion to Fortran 90 from Math 77 library
!
