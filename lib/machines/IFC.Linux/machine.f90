! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module MACHINE
! Useful functions and procedures that should be part of the Fortran standard.
! In a few cases they have in fact been made part of a later standard, but
! languish here still for backward ciompatibility.

  use IFPORT, only: getenv, getpid, getgid, system, ierrno
  use iso_c_binding
  ! ----------------- uncomment the following line with v15 -----------
  ! use ifposix ! only with Intel v15
  ! -------------------------------------------------------------------
  implicit none

  character(len=2) :: END_LINE = ' ' // char(10)
  character(len=1) :: FILSEP = '/'      ! '/' for Unix, '\' for DOS or NT
  integer, parameter :: HP = 0          ! Offset for first argument for GETARG

  public :: crash_burn, create_script, execute, execute_command_line
  public :: getarg, getenv, getids, is_a_directory
  interface
    subroutine GETARG ( ARGNUM, ARGVAL )
      integer, intent(in) :: ARGNUM  ! 0 = command name, 1 = first arg, etc.
      character(len=*), intent(out) :: ARGVAL   ! Blank if argnum out-of-range
    end subroutine GETARG
  end interface

  interface IO_ERROR; module procedure IO_ERROR_; end interface
  private :: IO_ERROR_
  ! private :: USleep
  ! These supply sleep and usleep functions from the c library
  ! which we call as procedures, thowing away their return values
  interface
    subroutine usleep ( i ) bind (c)
      ! 4.3BSD, POSIX.1-2001.  POSIX.1-2001 declares this function obsolete;
      ! use nanosleep(2) instead.  POSIX.1-2008 removes the specification of
      !  usleep().    
      use, intrinsic          :: iso_c_binding
      integer, value      :: i
    end subroutine usleep
    subroutine sleep ( i ) bind (c)
    use, intrinsic          :: iso_c_binding
      integer, value      :: i
    end subroutine sleep
  end interface

  logical, public, save :: NEVERCRASH = .true. ! Change to false for testing

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! -------------------------------------------------  CRASH_BURN  -----
  subroutine CRASH_BURN ( CODE, MSG, OPT )
  ! Print CODE and MSG if they're present.
  ! Try to make the processor crash, so it produces a walkback if it
  ! supports such a thing.
  ! DON'T put an optional argument corresponding to OPT.  It may be
  ! referenced on some systems in an attempt to make the processor crash.
  ! STOP as a last resort.
  ! Compile everything with -gline to get a walkback -- even with -O.
  
  ! Crash is provoked by one of available methods:
  ! (a) rewind illegal unit number (default)
  ! (b) acos(2.)
  ! You may select which method to use by setting msg = '!x (rest of message)'
  ! where x is the method
    integer, intent(in), optional :: CODE
    character(len=*), intent(in), optional :: MSG
    integer, intent(in), optional :: OPT(:)
    integer :: I = -huge(0)
    real :: arg
    character(len=*), parameter :: DEFAULTMETHOD = 'a' ! 'a' or 'b'
    character(len=1) :: method
    if ( present(code) ) print *, 'In ifc MACHINE%CRASH_BURN, CODE =', code
    if ( present(msg) ) print *, 'In ifc MACHINE%CRASH_BURN, MSG = ', trim(msg)
    method = DEFAULTMETHOD
    if ( NEVERCRASH ) then
      method = 's' ! Just stop, don't crash nor burn
    elseif ( present(msg) ) then
      if ( msg(1:1) == '!' ) method = adjustl(msg(2:))
    endif
    select case (method)
    case (' ', 'a', 'A')
      rewind i
    case ('b', 'B')
      arg = 1.
      if ( .not. present(opt) ) arg = arg + 1.
      arg = acos(arg)
      if ( arg > 0. ) print *, 'arg ', arg
    case default
      ! When you just want to land safely; e.g., a sips production run
    end select
    stop
  end subroutine CRASH_BURN

  ! --------------------------------------------------  create_script  -----
  subroutine create_script( Name, lines, thenRun, status, cmd_output, delay )
  ! Create the shell script named Name composed of lines
  ! Optionally run and then return the status and any output
    character(len=*), intent(in)                 :: Name    ! its path and name
    character(len=*), dimension(:), intent(in)   :: lines
    logical, optional, intent(in)                :: thenRun
    integer, optional, intent(out)               :: status
    character(len=*), optional, intent(out)      :: cmd_output
    integer, optional, intent(in)                :: delay
    ! Internal variables
    logical                                      :: mustRun
    ! Executable
    mustRun = .false.
    if ( present(thenRun) ) mustRun = thenRun
    call write_TEXTFILE ( Name, lines )
    call Execute( 'chmod a+x ' // trim(name), status )
    if ( mustRun ) call Execute( trim(name), status, cmd_output, delay )
  end subroutine create_script

  ! --------------------------------------------------  Execute  -----
  subroutine Execute( Command, status, cmd_output, delay )
  ! Execute the Command
  ! Return the status and optionally any output
    character(len=*), intent(in)            :: Command
    integer, intent(out)                    :: status
    character(len=*), optional, intent(out) :: cmd_output
    integer, optional, intent(in)           :: delay
    ! Internal variables
    integer                                 :: pid, gid
    integer                                 :: myDelay
    integer                                 :: tempstatus
    character(len=64)                       :: tempfilename
    ! Executable
    myDelay = 2
    if ( present(delay) ) myDelay = delay
    if ( .not. present(cmd_output) ) then
      call execute_command_line( Command, exitstat=status )
    else
      ! Approach:
      ! (1) construct a hopefully unique temporary file name
      ! (2) delete that file (in case it already exists)
      ! (3) redirect the stdout from Command to it
      ! (4) read its contents into cmd_output
      ! (5) delete it
      cmd_output = ' '
      call getids( pid, gid )
      write( tempfilename, * ) pid
      tempfilename = '/tmp/Execute_temp.' // adjustl(tempfilename)
      call execute_command_line( '/bin/rm -f ' // trim(tempfilename), &
        & exitstat=tempstatus )
      call execute_command_line( trim(Command) // ' > ' // trim(tempfilename), &
        & exitstat=status )
      call usleep( myDelay ) ! To allow time for os to carry out our command
      call read_textfile ( tempfilename, cmd_output )
      call execute_command_line( '/bin/rm -f ' // trim(tempfilename), &
        & exitstat=tempstatus )
    endif
  end subroutine Execute

  subroutine EXIT_WITH_STATUS ( STATUS )
  ! Exit and return STATUS to the invoking process
    integer, intent(in) :: STATUS
    call exit(status)
  end subroutine EXIT_WITH_STATUS
  
  ! subroutine getenv( name, value )
  ! Find value of env variable name
  !  character(len=*), intent(in)  :: name
  !  character(len=*), intent(out) :: value
  !  ! Executable
  !  call get_environment_variable( name, value )
  ! end subroutine getenv

  subroutine getids( process_id, group_id )
  ! return process and group ids
    integer, intent(out) :: process_id
    integer, intent(out) :: group_id
    ! Executable
    process_id = getpid()
    group_id   = getgid()
  end subroutine getids

  subroutine IO_ERROR_ ( MESSAGE, IOSTAT, FILE )
  ! Print MESSAGE and FILE, and then do something reasonable with IOSTAT.

    character(len=*), intent(in) :: MESSAGE
    integer, intent(in) :: IOSTAT
    character(len=*), intent(in), optional :: FILE

    integer :: L
    integer :: io_errnumber, sys_errnumber

    write (*,*) message(:len_trim(message))
    if ( present(file) ) then
      l = len_trim(file)
      write (*,*) file(:l)
    end if
    write (*,*) 'Error status code =', iostat
    call errsns( io_errnumber, sys_errnumber )
    write (*,*) 'io error number =', io_errnumber
    write (*,*) 'sys error number =', sys_errnumber
    return
  end subroutine IO_ERROR_

  ! --------------------------------------------------  is_a_directory  -----
  ! Return TRUE if and only path is a directory
  ! Note: returns FALSE if path is an ordinary file
  logical  function is_a_directory ( path )
    character(len=*), intent(in) :: path
    inquire( directory=trim(path), exist=is_a_directory )
  end function is_a_directory

  ! --------------------------------------------------  SHELL_COMMAND  -----
  ! Request Shell to execute ommand
  ! May optionally return
  ! Status             an integer to indicate status
  ! Error              an integer to indicate an error
  ! cmd_output         a character string to hold the result (e.g. from 'cat file')
  subroutine SHELL_COMMAND ( Command, Status, Error, cmd_output )
  ! Submit a character variable to the system as a shell command.

    character(len=*), intent(in) :: Command  ! The command
    integer, intent(out), optional :: Status ! Its status, if the system
                                        !  has such a concept, else zero
    integer, intent(out), optional :: Error  ! Error flag of the routine to submit
                                        ! the command, if the system has
                                        ! such a concept, else zero

    character(len=*), intent(out), optional :: cmd_output  ! The command's output
    integer :: MyError, MyStatus, RmStatus
    character(len=*), parameter :: tempfile = '/tmp/cmd_output.txt'
    ! Executable
    if ( present(cmd_output) ) then
      cmd_output = ' '
      ! Does /tmp exist?
      if ( is_a_directory( '/tmp' ) ) then
        MyStatus = system( Command // '> ' // tempfile )
        call read_textfile( tempfile, cmd_output )
        RmStatus = system( 'rm -f ' // tempfile ) ! Housekeeping
      else
      MyStatus = system( Command )
      endif
    else
      MyStatus = system( Command )
    endif
    if ( MyStatus == -1 ) then
      MyStatus = 0
      MyError = ierrno( )
    else
      MyStatus = 0
      MyError = 0
    endif
    if ( present(error) ) error = MyError
    if ( present(status) ) status = MyStatus
  end subroutine SHELL_COMMAND

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

  subroutine READ_TEXTFILE ( File, string )
  ! read a textfile into a single string
  ! Stolen mustly from io_stuff in lib
    character(len=*), intent(in)  :: File ! its path and name
    character(len=*), intent(inout) :: string    ! its contents
    ! Internal variables
    character(len=1), dimension(len(string)) :: cArray
    integer :: i
    integer :: lun
    integer :: pos
    integer :: recrd
    integer :: status
    character(len=12) :: xfmt
    character(len=8) :: xlen
    ! Executable
    ! What format do we use for reading each line?
    xfmt = '(128a1)' ! This is the default; if lines are larger supply maxLineLen
    write( xlen, '(i8)' ) len(string)
    if ( index(xlen, '*') < 1 ) xfmt = '(' // trim(adjustl(xlen)) // 'a1)'
    ! Try to read the textfile
    open( newunit=lun, form='formatted', &
      & file=trim(File), status='old', iostat=status )
    if ( status /= 0 ) then
      write(*,*) 'machine%READ_TEXTFILE_ARR-E- Unable to open textfile ' // &
        & trim(File)
      return
    endif
    recrd = 0
    ! print *, 'xfmt: ', xfmt
    i = 0
    do
      status = 0
      carray = achar( 0 )
      read( UNIT=lun, fmt=xfmt, eor=50, end=500, err=50, advance='no' ) cArray
500   status = -1
50    if ( status /= 0 ) exit
      ! print *, cArray
      recrd = recrd + 1
      oneLine: do pos=1, len(string) - 1
        if ( any(carray(pos:pos+1) == achar(0)) ) exit oneLine
        i = min(i + 1, len(string))
        string(i:i) = carray(pos)
       enddo oneLine
      i = min(i + 1, len(string))
      string(i:i) = achar(13)
    enddo
    close ( lun )
  end subroutine READ_TEXTFILE
  
  subroutine write_TEXTFILE ( File, lines )
  ! write lines to a textfile
    character(len=*), intent(in)               :: File     ! its path and name
    character(len=*), intent(in), dimension(:) :: lines    ! its contents
    ! Internal variables
    integer :: i, lun, status
    ! Executable
    open( newunit=lun, form='formatted', &
      & file=trim(File), status='unknown', access='sequential', &
      & recl=size(lines)*len(lines(1)) + 1, iostat=status )
    if ( status /= 0 ) then
      write(*,*) 'machine%write_TEXTFILE_ARR-E- Unable to open textfile ' // &
        & trim(File)
      return
    endif
    do i=1, size(lines)
      write ( lun, '(a)', advance='yes' ) trim(lines(i))
    enddo
    close( UNIT=lun, iostat=status )
  end subroutine write_TEXTFILE

  ! ----------------- comment-out the following routine with v15 -----------
  ! Because v14 and earlier don't support this
  ! subroutine execute_command_line( Command, wait, exitstat )
  !   ! Dummy args
  !   character(len=*), intent(in)   :: Command
  !   logical, optional, intent(in)  :: wait
  !   integer, optional, intent(out) :: exitstat
  !   ! Internal variables
  !   integer                        :: status
  !   ! Executable
  !   status = system( Command )
  !   if ( present(exitstat) ) exitstat = status
  ! end subroutine execute_command_line
  
! ----------------------------------------------  not_used_here  -----
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MACHINE

! $Log$
! Revision 1.14  2016/02/29 19:45:52  pwagner
! Exploit c bindings for usleep, sleep
!
! Revision 1.13  2016/02/28 23:59:46  pwagner
! Added USleep needed here
!
! Revision 1.12  2016/02/26 19:37:27  pwagner
! Added create_script, execute, execute_command_line, getenv, getids
!
! Revision 1.11  2015/08/12 20:17:54  pwagner
! Added Is_A_Directory; Shell_Command now optionally returns stdout from cmd
!
! Revision 1.10  2009/09/15 20:01:12  pwagner
! Intel compiler crashes with walkback only for method 'a'
!
! Revision 1.9  2009/06/23 19:58:52  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 1.8  2008/01/16 17:48:10  pwagner
! Fleshed out ifort call to system
!
! Revision 1.7  2007/05/22 20:58:17  vsnyder
! Set ERROR and STATUS to avoid warning messages
!
! Revision 1.6  2007/04/24 22:30:27  pwagner
! Updated to include new ifc support for calls to exit, system
!
! Revision 1.5  2005/06/22 20:26:22  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 1.4  2004/08/19 00:13:12  pwagner
! Added crash_burn to provoke crash with walkback
!
! Revision 1.3  2003/04/01 22:53:02  pwagner
! Chucked stub getarg--the v7.0 ifc has getarg
!
! Revision 1.2  2002/05/02 09:49:55  hcp
! Just syncing up
!
! Revision 1.1  2001/09/17 10:45:18  pumphrey
! Added machine.f90 file for Intel ifc fortran
