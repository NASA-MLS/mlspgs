! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MACHINE
! This file must be preprocessed through a makefile
! containing sed commands to snip out offending lines
! (delimited by BAD_MATCH) depending upon whether to build 
! (1) a version allowing garbage collection (-gc on the link line); or
! (2) a version not allowing garbage collection
!
! The recipe for getting each of the versions above is
!     if  (1) set BAD_MATCH='no -gc'
! else if (2) set BAD_MATCH='-gc'
!
!     sed "/Start $BAD_MATCH/,/End $BAD_MATCH/ d" machine.f90 \
!       > my_machine.f90
! and then just compile my_machine.f90

!---------- Start -gc section
!           requires -gc among LDOPTS
!   (the following lines automatically deleted for version (2))
  use F90_GC, only: GCOLLECT, NCOLLECTIONS, &
   & DONT_EXPAND, DONT_GC, FULL_FREQUENCY, MAX_RETRIES, SILENT_GC
!---------- End -gc section
  use F90_IOSTAT				! everything; see iostat_msg_NAG
  use F90_UNIX_ENV, only: IARGC, NAG_GETARG => GETARG
  ! Exit and return an integer status to the invoking process
  use F90_UNIX_PROC, only: EXIT_WITH_STATUS => EXIT, SYSTEM
  implicit none

  character(LEN=2) :: END_LINE = ' ' // char(10)
  character(LEN=1) :: FILSEP = '/'      ! '/' for Unix, '\' for DOS or NT
  integer, parameter :: HP = 0          ! Offset for first argument for GETARG

  interface IO_ERROR; module procedure IO_ERROR_; end interface
  private IO_ERROR_
  public :: SHELL_COMMAND
  public :: MLS_DISABLE_AUTOGC, MLS_GC_NOW, MLS_HOWMANY_GC, MLS_CONFIG_GC

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

  subroutine IO_ERROR_ ( MESSAGE, IOSTAT, FILE )
  ! Print MESSAGE and FILE, and then do something reasonable with IOSTAT.
    character(LEN=*), intent(in) :: MESSAGE
    integer, intent(in) :: IOSTAT
    character(LEN=*), intent(in), optional :: FILE

    integer :: L
!   character(LEN=127) :: MSG

    write (*,*) message(:len_trim(message))
    if ( present(file) ) then
      l = len_trim(file)
      if ( l /= 0 ) write (*,*) file(:l)
    end if
    call print_iostat_msg_NAG (iostat) ! my version of a Lahey intrinsic
!   write (*,*) msg(:len_trim(msg))    ! Print the error message
    write (*,*) 'Error status code =', iostat
    return
  end subroutine IO_ERROR_

  subroutine GETARG ( ARGNUM, ARGVAL )
    integer, intent(in) :: ARGNUM  ! 0 = command name, 1 = first arg, etc.
    character(len=*), intent(out) :: ARGVAL   ! Blank if argnum out-of-range
    integer :: STATUS
    call nag_getarg ( argnum, argval, errno = status )
    if ( status /= 0 ) argval = ' '
  end subroutine GETARG

  subroutine Print_iostat_msg_NAG ( IOSTAT )
    integer, intent(in) :: IOSTAT

    ! This subroutine was made from NAG's own f90_iostat.f90
    ! found in /usr/local/lib/NAGWare
    ! then edited via the following 2 steps:

    ! sed -n -f machines/iostat.sed machines/f90_iostat.f90
    !                            > machines/iostat.f90

    ! cat -n machines/iostat.f90 | sort -n -r
    !               | sed -n 's/^.......//p' > machines/iostat_rev.f90

    ! where the file iostat.sed contains the following 5 lines
    ! (w/o the leading !):

    !s/integer, parameter :: \([A-Za-z_][A-Za-z_0-9]*\)/case(\1)/
    !s/'//g
    !s/^! [A-Z].*$/         print*, '&'/
    !s/=/!/
    !p

    ! after this lines from iostat_rev.f90 were copied and pasted below
	
    select case (iostat)
	! start of included file iostat_rev.f90
    case(IOERR_SCALE_FOLLOWED_BY_REPEAT)  ! 217
      print*, '! Scale factors cannot be followed by repeat count in FMT! specifier'
    case(IOERR_DIRECT_POSITION_INCOMPAT)  ! 216
      print*, '! Direct access is incompatible with the POSITION! specifier'
    case(IOERR_REAL_INPUT_OVERFLOW)       ! 215
      print*, '! Floating overflow on real number input'
    case(IOERR_END_OF_DIRECT_ACCESS)      ! 214
      print*, '! READ beyond end of direct access file on unit %d'
    case(IOERR_RECL_LE_ZERO)              ! 213
      print*, '! Invalid value for RECL! specifier (must be positive)'
    case(IOERR_REPLACE_OR_NEW_NEED_FILE)  ! 212
      print*, '! No FILE! specifier with STATUS=REPLACE or STATUS=NEW'
    case(IOERR_INTEGER64_TOO_BIG)         ! 211
      print*, '! Input value too large for 64-bit integer'
    case(IOERR_RW_AFTER_ENDFILE)          ! 210
      print*, '! READ/WRITE attempted after ENDFILE on unit %d'
    case(IOERR_TRUNCATE_FAILED)           ! 209
      print*, '! File truncation on unit %d failed'
    case(IOERR_CORRUPT_UNFORMATTED_FILE)  ! 208
      print*, '! Unformatted data file open on unit %d is corrupt'
    case(IOERR_RECORD_TOO_SHORT)          ! 207
      print*, '! Record too short for format requirement and PAD!NO on unit %d'
    case(IOERR_INPUT_LIST_TOO_BIG)        ! 206
      print*, '! Input list bigger than record length in unformatted READ on unit %d'
    case(IOERR_NO_NAMELIST_GROUP_NAME)    ! 205
      print*, '! Missing namelist group name in input of NAMELIST/%s/'
    case(IOERR_ZERO_STRIDE)               ! 204
      print*, '! Section stride has value zero in input for object %s of NAMELIST/%s/'
    case(IOERR_ZERO_SIZE_INPUT)           ! 203
      print*, '! Array section has zero size in input for object %s of NAMELIST/%s/'
    case(IOERR_SUBSCRIPT_OUT_OF_RANGE)    ! 202
      print*, '! Subscript (%d) out of range in input for object %s of NAMELIST/%s/'
    case(IOERR_EXPECTED_RPAREN)           ! 201
      print*, '! Expected ) but found %c in input for object %s of NAMELIST/%s/'
    case(IOERR_EXPECTED_COMMA)            ! 200
      print*, '! Expected , but found %c in input for object %s of NAMELIST/%s/'
    case(IOERR_SUBSTRING_OUT_OF_BOUNDS)   ! 199
      print*, '! Substring (%d:%d) out of bounds in input for object %s of NAMELIST/%s/'
    case(IOERR_ZERO_LENGTH_INPUT)         ! 198
      print*, '! Substring has zero length in input for object %s of NAMELIST/%s/'
    case(IOERR_EXPECTED_COLON)            ! 197
      print*, '! Expected : but found %c in input for object %s of NAMELIST/%s/'
    case(IOERR_BAD_INTEGER_LITERAL)       ! 196
      print*, '! Invalid integer literal in input for object %s of NAMELIST/%s/'
    case(IOERR_ARRAY_OF_ARRAY)            ! 195
      print*, '! Array component of array parent in input for object %s of NAMELIST/%s/'
    case(IOERR_UNKNOWN_COMPONENT)         ! 194
      print*, '! Unknown component \"%s\" in input for object %s of NAMELIST/%s/'
    case(IOERR_COMPONENT_NAME_TOO_LONG)   ! 193
      print*, '! Component name too long in input for object %s of NAMELIST/%s/'
    case(IOERR_UNEXPECTED_COMPONENT)      ! 192
      print*, '! Unexpected component specifier for object %s of NAMELIST/%s/'
    case(IOERR_UNEXPECTED_SUBSCRIPT)      ! 191
      print*, '! Unexpected subscript for object %s of NAMELIST/%s/'
    case(IOERR_UNKNOWN_OBJECT_NAME)       ! 190
      print*, '! Unknown group object \"%s\" in input for NAMELIST/%s/'
    case(IOERR_EXPECTED_EQUALS)           ! 189
      print*, '! Expected ! but found %c in NAMELIST input'
    case(IOERR_OBJECT_NAME_TOO_LONG)      ! 188
      print*, '! NAMELIST group name in input of NAMELIST/%s/ is too long'
    case(IOERR_NAMELIST_BAD_CHAR)         ! 187
      print*, '! Invalid character %c in NAMELIST input'
    case(IOERR_WRONG_NAMELIST)            ! 186
      print*, '! Expected NAMELIST group /%s/ but found /%s/'
    case(IOERR_GROUP_NAME_TOO_LONG)       ! 185
      print*, '! NAMELIST group name in input is too long'
    case(IOERR_NO_AMPERSAND)              ! 184
      print*, '! Expected & but found %c in NAMELIST input'
    case(IOERR_CANNOT_OPEN)               ! 183
      print*, '! Unknown OPEN failure on unit'
    case(IOERR_NO_INPUT_LOGICAL)          ! 182
      print*, '! No value found in LOGICAL input field'
    case(IOERR_BAD_INPUT_LOGICAL)         ! 181
      print*, '! Illegal character in LOGICAL input field'
    case(IOERR_BAD_REC)                   ! 180
      print*, '! Record number out of range'
    case(IOERR_NOT_DIRECT)                ! 179
      print*, '! Unit is not connected for DIRECT i/o'
    case(IOERR_BACKSPACE_FAILED)          ! 178
      print*, '! BACKSPACE failed to find the beginning of the previous record'
    case(IOERR_CANNOT_REWIND)             ! 177
      print*, '! File connected to unit is not capable of REWIND'
    case(IOERR_NEW_FILE_EXISTS)           ! 176
      print*, '! NeW file already exists'
    case(IOERR_NO_OLD_FILE)               ! 175
      print*, '! Cannot find OLD file'
    case(IOERR_NAME_TOO_LONG)             ! 174
      print*, '! File name too long'
    case(IOERR_ENDFILE_TWICE)             ! 173
      print*, '! ENDFILE applied twice to unit with no intervening file positioning'
    case(IOERR_CANNOT_KEEP)               ! 172
      print*, '! STATUS!KEEP is invalid for a SCRATCH file'
    case(IOERR_NO_RECL)                   ! 171
      print*, '! The RECL! specifier must be given for DIRECT access OPEN'
    case(IOERR_BAD_PAD)                   ! 170
      print*, '! Invalid value for PAD! specifier'
    case(IOERR_BAD_DELIM)                 ! 169
      print*, '! Invalid value for DELIM! specifier'
    case(IOERR_BAD_ACTION)                ! 168
      print*, '! INvalid value for ACTION! specifier'
    case(IOERR_BAD_POSITION)              ! 167
      print*, '! Invalid value for POSITION! specifier'
    case(IOERR_BAD_BLANKS)                ! 166
      print*, '! Invalid value for BLANKS! specifier'
    case(IOERR_BAD_FORM)                  ! 165
      print*, '! Invalid value for FORM! specifier'
    case(IOERR_BAD_ACCESS)                ! 164
      print*, '! Invalid value for ACCESS! specifier'
    case(IOERR_BAD_STATUS)                ! 163
      print*, '! Invalid value for STATUS! specifier'
    case(IOERR_DIFFERENT_POSITION)        ! 162
      print*, '! OPEN on connected unit has different POSITION! specifier'
    case(IOERR_DIFFERENT_ACTION)          ! 161
      print*, '! OPEN on connected unit has different ACTION! specifier'
    case(IOERR_DIFFERENT_RECL)            ! 160
      print*, '! OPEN on connected unit has different RECL! specifier'
    case(IOERR_DIFFERENT_FORM)            ! 159
      print*, '! OPEN on connected unit has different FORM! specifier'
    case(IOERR_DIFFERENT_ACCESS)          ! 158
      print*, '! OPEN on connected unit has different ACCESS! specifier'
    case(IOERR_DIFFERENT_STATUS)          ! 157
      print*, '! OPEN on connected unit with STATUS! specifier must have STATUS=OLD'
    case(IOERR_SCRATCH_NAMED)             ! 156
      print*, '! FILE! specifier on OPEN with STATUS=SCRATCH'
    case(IOERR_OLD_UNCONNECTED_NEED_FILE) ! 155
      print*, '! Unit is not connected on OPEN with STATUS!OLD and no FILE= specifier'
    case(IOERR_NOT_UNFORMATTED)           ! 154
      print*, '! Unit is not connected for UNFORMATTED i/o'
    case(IOERR_NOT_WRITE)                 ! 153
      print*, '! Unit is not connected for WRITE action'
    case(IOERR_NOT_FORMATTED)             ! 152
      print*, '! Unit is not connected for FORMATTED i/o'
    case(IOERR_NOT_READ)                  ! 151
      print*, '! Unit is not connected for READ action'
    case(IOERR_NOT_SEQUENTIAL)            ! 150
      print*, '! Unit is not connected for SEQUENTIAL i/o'
    case(IOERR_CANNOT_BACKSPACE)          ! 149
      print*, '! File connected to unit is not capable of BACKSPACE'
    case(IOERR_NOT_CONNECTED)             ! 148
      print*, '! Unit is not connected'
    case(IOERR_BAD_UNIT)                  ! 147
      print*, '! Unit number out of range'
    case(IOERR_READ_AFTER_WRITE)          ! 146
      print*, '! READ after WRITE with no intervening file positioning'
    case(IOERR_BAD_EDIT_FOR_CHARACTER)    ! 145
      print*, '! Invalid edit descriptor for character i/o-list item'
    case(IOERR_BAD_HEX)                   ! 144
      print*, '! Invalid character in hexadecimal integer input field'
    case(IOERR_BAD_OCTAL)                 ! 143
      print*, '! Invalid character in octal integer input field'
    case(IOERR_BAD_BINARY)                ! 142
      print*, '! Invalid character in binary integer input field'
    case(IOERR_BAD_INPUT_INTEGER)         ! 141
      print*, '! Invalid character in integer input field'
    case(IOERR_BAD_INPUT_REAL)            ! 140
      print*, '! Invalid character in real input field'
    case(IOERR_BAD_INPUT_EXPONENT)        ! 139
      print*, '! Invalid exponent in real input field'
    case(IOERR_ONLY_SIGN_FOUND)           ! 138
      print*, '! Sign in a numeric input field not followed by any digits'
    case(IOERR_CHAR_OVERLAPS_END)         ! 137
      print*, '! Character string edit descriptor does not terminate before format end'
    case(IOERR_BAD_EDIT_FOR_LOGICAL)      ! 136
      print*, '! Invalid edit descriptor for logical i/o-list item'
    case(IOERR_BAD_EDIT_FOR_INTEGER)      ! 135
      print*, '! Invalid edit descriptor for integer i/o-list item'
    case(IOERR_BAD_EDIT_FOR_REAL)         ! 134
      print*, '! Invalid edit descriptor for real i/o-list item'
    case(IOERR_CHAR_EDIT_IN_READ)         ! 133
      print*, '! Character string edit descriptor used on input'
    case(IOERR_REPEAT_FOR_POSITION)       ! 132
      print*, '! Repeat factor given for position edit descriptor'
    case(IOERR_NO_SPACING_FOR_X)          ! 131
      print*, '! No spacing specified for X edit descriptor'
    case(IOERR_REPEAT_FOR_CHAR_EDIT)      ! 130
      print*, '! Repeat factor given for character string edit descriptor'
    case(IOERR_MISSING_HOLLERITH_LENGTH)  ! 129
      print*, '! Missing length of H edit descriptor'
    case(IOERR_REPEAT_FOR_BLANK_INTERP)   ! 128
      print*, '! Repeat factor given for blank-interpretation edit descriptor'
    case(IOERR_REPEAT_FOR_SIGN)           ! 127
      print*, '! Repeat factor given for sign edit descriptor'
    case(IOERR_NO_EDIT_FOR_REPEAT)        ! 126
      print*, '! No edit descriptor following repeat factor'
    case(IOERR_BAD_EDIT)                  ! 125
      print*, '! Invalid edit descriptor'
    case(IOERR_BAD_BNBZ)                  ! 124
      print*, '! Expected N or Z following B in format specification'
    case(IOERR_EXPECTED_P)                ! 123
      print*, '! Expected P following signed integer constant in format specification'
    case(IOERR_EXPECTED_PERIOD)           ! 122
      print*, '! Expected decimal point in format specification'
    case(IOERR_FORMAT_MBNZ)               ! 121
      print*, '! Field/exponent width or repeat in format specification must be non-zero'
    case(IOERR_EXPECTED_INTEGER_VALUE)    ! 120
      print*, '! Expected integer literal constant in format specification'
    case(IOERR_UNEXPECTED_FORMAT_END)     ! 119
      print*, '! Unexpected end of format specification'
    case(IOERR_SUBFMT_TOO_DEEP)           ! 118
      print*, '! Sub-format groups nested too deeply'
    case(IOERR_NO_DATA_EDIT_IN_REVERSION) ! 117
      print*, '! No data edit descriptor in tail of format specification after reversion'
    case(IOERR_FORMAT_NO_ENDING_RPAREN)   ! 116
      print*, '! Format specification does not end with a right parenthesis'
    case(IOERR_FORMAT_NO_LPAREN)          ! 115
      print*, '! Format specification does not begin with a left parenthesis'
    case(IOERR_BAD_CHAR)                  ! 114
      print*, '! Invalid input for character editing'
    case(IOERR_BAD_COMPLEX)               ! 113
      print*, '! Invalid input for complex editing'
    case(IOERR_BAD_LOGICAL)               ! 112
      print*, '! Invalid input for logical editing'
    case(IOERR_BAD_REAL)                  ! 111
      print*, '! Invalid input for real editing'
    case(IOERR_INTEGER_TOO_BIG)           ! 110
      print*, '! Input value too large for default INTEGER type'
    case(IOERR_INTEGER_OVERFLOW_REPEAT)   ! 109
      print*, '! Repeat factor in list-directed input larger than HUGE(0)'
    case(IOERR_INTEGER2_TOO_BIG)          ! 108
      print*, '! Input value too large for INTEGER(KIND!2)'
    case(IOERR_INTEGER1_TOO_BIG)          ! 107
      print*, '! Input value too large for INTEGER(KIND!1)'
    case(IOERR_BAD_INTEGER)               ! 106
      print*, '! Invalid input for integer editing'
    case(IOERR_ZERO_REPEAT)               ! 105
      print*, '! Zero repeat factor in list-directed input'
    case(IOERR_INPUT_BUFFER_OVERFLOW)     ! 104
      print*, '! Record too long for input buffer'
    case(IOERR_BAD_EXPONENT)              ! 103
      print*, '! Exponent too large for w.d format'
    case(IOERR_BAD_SCALE)                 ! 102
      print*, '! Scale factor out of range'
    case(IOERR_INTERNAL_FILE_OVERFLOW)    ! 101
      print*, '! Internal file overflow'
    case(IOERR_BUFFER_OVERFLOW)           ! 100
      print*, '! Buffer overflow on output'

    !  Compiler run-time system
    !  The following IOSTAT values are defined for the NAGWare f90

    !  /usr/include/sys/errno.h.
    !  status returns.  These are usually documented in the file
    !  IOSTAT values between 1 and 99 are reserved for host system

    case(IOERR_OK)                        !   0
    case(IOERR_EOF)                       !  -1
    case(IOERR_EOR)                       !  -2

    !  The following IOSTAT values are always

    ! end of included file iostat_rev.f90
    case default
      print*, '! IOSTAT not defined in f90_iostat.f90'
    end select

  end subroutine print_iostat_msg_NAG

  subroutine SHELL_COMMAND ( Command, Status, Error )
  ! Submit a character variable to the system as a shell command.

    character(len=*), intent(in) :: Command  ! The command
    integer, intent(out), optional :: Status ! Its status, if the system
                                        !  has such a concept, else zero
    integer, intent(out), optional :: Error  ! Status of the routine to submit
                                        ! the command, if the system has
                                        ! such a concept, else zero

    integer :: MyError, MyStatus

    call system ( command, myStatus, myError)
    if ( present(error) ) error = myError
    if ( present(status) ) status = myStatus
  end subroutine SHELL_COMMAND

!---------- Start -gc section
!           requires -gc among LDOPTS
!   (the following lines automatically deleted for version (2))
  subroutine MLS_CONFIG_GC ( EXPAND, FREQUENCY, RETRIES, SILENT )
  ! Configures various parameters affecting garbage collection

    logical, optional, intent(in) :: expand    ! Autoexpand heap?
    integer, optional, intent(in) :: frequency ! How many incremental colls. betw. fulls
    integer, optional, intent(in) :: retries   ! How many attempts before giving up
    logical, optional, intent(in) :: silent    ! Quash report on each collection?
    if ( present(expand) ) dont_expand = .not. expand
    if ( present(frequency) ) full_frequency = frequency
    if ( present(retries) ) max_retries = retries
    if ( present(silent) ) silent_gc = silent
  end subroutine MLS_CONFIG_GC

  subroutine MLS_DISABLE_AUTOGC ( Which )
  ! Turns automatic garbage collection on/off

    character(len=*), intent(in) :: Which  ! 'On' or 'Off'
    if ( Which == 'On' .or. Which == 'ON' &
      & .or. Which == 'on' ) then
      DONT_GC = .false.
    else
      DONT_GC = .true.
    endif
  end subroutine MLS_DISABLE_AUTOGC

  subroutine MLS_GC_NOW
  ! Manually collects garbage when called

    CALL GCOLLECT
  end subroutine MLS_GC_NOW

  integer function MLS_HOWMANY_GC()
  ! Returns how many garbage collections have been performed

    MLS_HOWMANY_GC = NCOLLECTIONS()
  end function MLS_HOWMANY_GC
!---------- End -gc section
!---------- Start no -gc section
!           forbids -gc among LDOPTS
!   (the following lines automatically deleted for version (2))
  ! ----------------------------------------------
  ! The following are merely introduced to satisfy 
  ! NAG call interfaces to f90_gc
  ! Because we're assuming the link statement will lack -gc
  ! we have non-functional substitutes
  subroutine MLS_CONFIG_GC ( EXPAND, FREQUENCY, RETRIES, SILENT )
  ! Configures various parameters affecting garbage collection
   logical dont_EXPAND, SILENT_GC
   integer full_frequency, max_retries
    logical, optional, intent(in) :: expand    ! Autoexpand heap?
    integer, optional, intent(in) :: frequency ! How many incremental colls. betw. fulls
    integer, optional, intent(in) :: retries   ! How many attempts before giving up
    logical, optional, intent(in) :: silent    ! Quash report on each collection?
    if ( present(expand) ) dont_expand = .not. expand
    if ( present(frequency) ) full_frequency = frequency
    if ( present(retries) ) max_retries = retries
    if ( present(silent) ) silent_gc = silent
  end subroutine MLS_CONFIG_GC

  subroutine MLS_DISABLE_AUTOGC ( Which )
  ! Turns automatic garbage collection on/off
    character(len=*), intent(in) :: Which  ! 'On' or 'Off'
    logical dont_gc
    if ( Which == 'On' .or. Which == 'ON' &
      & .or. Which == 'on' ) then
      DONT_GC = .false.
    else
      DONT_GC = .true.
    endif
  end subroutine MLS_DISABLE_AUTOGC

  subroutine MLS_GC_NOW
  ! Manually collects garbage when called
  !    CALL GCOLLECT
  end subroutine MLS_GC_NOW

  integer function MLS_HOWMANY_GC()
  ! Returns how many garbage collections have been performed
    MLS_HOWMANY_GC = 0 ! NCOLLECTIONS()
  end function MLS_HOWMANY_GC
!---------- End no -gc section

end module MACHINE

! $Log$
! Revision 1.8  2002/02/05 00:37:19  pwagner
! Added garbage collection features
!
! Revision 1.7  2002/01/31 19:16:28  pwagner
! Brought up-to-date with shell_command using library as appropriate; untested
!
! Revision 1.6  2002/01/31 00:47:47  pwagner
! Removed interface in shell_command
!
! Revision 1.5  2002/01/30 19:50:30  vsnyder
! Added Shell_Command subroutine
!
! Revision 1.4  2001/05/15 20:34:37  pwagner
! Compatible with NAG f95-V4.1
!
! Revision 1.3  2001/05/04 23:25:17  vsnyder
! Added Exit_With_Status routine
!
! Revision 1.2  2001/03/21 00:42:20  pwagner
! Added print_iostat_msg_NAG
!
! Revision 1.1  2001/01/13 00:29:44  pwagner
! moved to lib/machines/MLSCONFG/machine.f90
!
! -- Revision history of earlier incarnations archived elsewhere--
!
! Revision 1.1  2000/10/19 17:40:52  pwagner
! first commit
!
! Revision 1.2  2000/10/12 22:54:12  vsnyder
! Correct a commented-out line that may get commented-in for another
! computer/os/compiler combination
!
! Revision 1.1  2000/10/12 22:21:11  vsnyder
! Change directory name from NAG to pclinuxNAG
!
! Revision 1.3  2000/10/09 22:15:55  vsnyder
! Moved machine.f90 from l2 to lib
!
! Revision 2.0  2000/09/05 18:58:04  ahanzel
! Changing file revision to 2.0.
