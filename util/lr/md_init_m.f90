module MD_INIT_M

! Machine-dependent initialization.  This is a Unix version.
  use IO, only: OUNIT, PRUNIT
  use MI_INIT_M, only: OPEN_INPUT, OPEN_LISTING, OPEN_OUTPUT
  use TOGGLES, only : TOGGLE

  implicit NONE
  private

  public :: MD_INIT

  character, parameter :: FILSEP = '/'

contains

  subroutine MD_INIT

  ! Process the command line.

    character(len=127) :: ARG      ! Argument from the command line
    character(len=4) :: EXTS(3) = (/ '.grm', '.lr ', '.lls' /)
    character(len=127) :: FILES(3) ! 1 = input, 2 = output, 3 = listing
    integer :: I                   ! Subscript and loop inductor
    integer :: IARG                ! The argument number
    integer :: IFILE               ! Which file
    integer :: LDOT, LDOT1, LSEP   ! Positions in file names
    logical :: OPTS                ! Options noticed until option --
    integer :: TEST                ! Zero or 1, for showing toggle states

    files = ' '
    iarg = 1
    ifile = 0
    opts = .true.
    do
      call get_command_argument ( iarg, arg )
      iarg = iarg + 1
      if ( arg(1:3) == '-- ' ) opts = .false.
      if ( opts .and. arg(1:1) == '-' ) then  ! Process options
        do i = 2, len_trim(arg)
          if ( arg(i:i) == '?' .or. arg(i:i) == 'h' .or. arg(i:i) == 'H' ) then
            call usage
          else
            toggle(ichar(arg(i:i))) = 1 - toggle(ichar(arg(i:i)))
          end if
        end do
      else if ( arg /= ' ' ) then
        ifile = ifile + 1
        if ( ifile > size(files) ) then
          write ( *, * ) 'More than ', size(files), &
                         ' file names -- excess ignored.'
        else
          files(ifile) = arg
        end if
      else
        exit
      end if
    end do

    ! The command line has been read.  Now open the files.

    call places ( 1 )
    if ( ldot == 0 ) then
      ! There's no . after the last file separator.  Add '.grm' on the
      ! end of the input file name.
      ldot = len_trim(files(1)) + 1
      if ( files(1) /= ' ' ) files(1)(ldot:) = exts(1)
    end if
    ldot1 = ldot
    do ifile = 2, 3 ! output and listing files
      ! If there's no file, use the input file's base name
      if ( files(ifile) == ' ' ) files(ifile) = files(1)(:ldot1-1)
      i = len_trim(files(ifile))
      ! If the output file name ends with a file separator, use the base name
      ! of the input file without its path.
      if ( files(ifile)(i:i) == filsep ) &
        files(ifile)(i+1:) = files(1)(lsep+1:ldot1-1)
      call places ( ifile )
      i = len_trim(files(ifile))
      if ( ldot == 0 .and. files(ifile) /= ' ' ) &
        files(ifile)(i+1:) = exts(ifile)
    end do
    if ( toggle(ichar('d')) + toggle(ichar('D')) /= 0 ) then
    ! Debugging output:
      write ( *, '(a)' ) 'File names: '
      write ( *, '(1x,a)' ) ( trim(files(ifile)), ifile = 1, 3 )
      ldot = 0
      do i = lbound(toggle,1), ubound(toggle,1)
        test = merge(1,0,any([i==iachar('F'),i==iachar('L'),i==iachar('M'), &
                            & i==iachar('N'),i==iachar('R'),i==iachar('X')]))
        if ( toggle(i) /= test ) then
          if ( ldot == 0 ) write ( *, '(a)', advance='NO' ) 'Toggles: '
          ldot = 1
          write ( *, '(a)', advance='NO' ) achar(i)
        end if
      end do
      if ( ldot /= 0 ) write ( *, * )
    end if
    call open_input ( files(1) )
    call open_output ( files(2) )
    if ( toggle(ichar('s')) /= 0 ) then
      ounit = 6
    else
      ounit = prunit
      call open_listing ( files(3) )
    end if

  contains

    subroutine PLACES ( FINDEX )
    ! Look backward in the input file name for the last . after the last
    ! file separator (/ or \ depending on Unix or Dos-ish).
      integer, intent(in) :: FINDEX

      ldot = 0
      lsep = 0
      do i = len_trim(files(findex)), 1, -1
        if ( files(findex)(i:i) == filsep ) then
          lsep = i
          exit
        end if
        if ( files(findex)(i:i) == '.' .and. ldot == 0 ) ldot = i
      end do
    end subroutine PLACES

    subroutine Usage
      call get_command_argument ( 0, arg )
      print '(a)', 'Usage: ' // trim(arg) // ' [options] input output listing'
      print '(a)', '  ".grm" is appended onto input if there is no extension'
      print '(a)', '  ".lr" is appended onto output if there is no extension'
      print '(a)', '  ".lls" is appended onto listing if there is no extension'
      print '(a)', '  options:'
      print '(a)', '    -- => no more options'
      print '(a)', '    -?hH => this output'
      print '(a)', '    -<anything else>* inverts toggle:'
      print '(a)', '    -d or -D => debugging output'
      print '(a)', '    -F => output Fortran parameter statements (default on),'
      print '(a)', '          else output a formatted file'
      print '(a)', '    -L => output the tables (default on)'
      print '(a)', '    -M => print parsing machine (default on)'
      print '(a)', '    -N => print grammar neatly and quit (default off)'
      print '(a)', '    -R => print grammar neatly (default on)'
      print '(a)', '    -X => print cross reference (default on)'
      print '(a)', '    -2 => print some debugging information about the'
      print '(a)', '          configuration analysis phase (default off)'
      print '(a)', '    -3 => print more debugging information about the'
      print '(a)', '          configuration analysis phase (default off)'
      stop
    end subroutine Usage

  end subroutine MD_INIT

end module MD_INIT_M
