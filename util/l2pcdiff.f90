! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=================================
program l2pcdiff ! show diffs between swaths in two different files
!=================================

   use Declaration_table, only: allocate_decl, deallocate_decl, dump_decl
   use Dump_Options, only: dumpTableSide, rmsFormat
   use Hdf, only: dfacc_rdonly
   use HighOutput, only: OutputNamedValue
   use Init_tables_module, only: init_tables
   use Intrinsic, only: l_hdf
   use L2PC_m, only: Diff
   use Lexer_core, only: init_lexer
   use Machine, only: hp, getarg
   use MLSCommon, only: MLSFile_T
   use MLSFiles, only: HDFVersion_5, &
     & InitializeMLSFile, MLS_Exists, MLS_Inqswath
   use MLSHDF5, only: MLS_H5Open, MLS_H5Close
   use MLSMessageModule, only: MLSMessageConfig
   use MLSStringLists, only: catLists, ExpandStringRange
   use MLSStrings, only: WriteIntsToChars
   use Output_m, only: resumeOutput, suspendOutput, output
   use PrintIt_m, only: Set_Config
   use SDPToolkit, only: UseSDPToolkit
   use Time_M, only: Time_Now, time_config
   use Toggles, only: switches
   use Tree, only: allocate_tree, deallocate_tree, print_subtree
   
   implicit none

!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it
! LF95.Linux/test [options] [input files]
  integer, parameter ::          MAXFIELDSLENGTH = 256
  type options_T
    character(len=MAXFIELDSLENGTH) :: diffOptions = ''
    character(len=8) :: HEADSIDE = 'left' ! on which side stats headers printed
    integer     ::         Details         = 1
    logical     ::         timing          = .false.
    logical     ::         verbose         = .false.
    logical     ::         debug           = .false.
  end type options_T
  
  type ( options_T ) :: options

  integer, parameter                      :: MAXFILES = 100
  integer, parameter                      :: MAXNCHUNKS = 50
  character(len=*), parameter             :: SWATHPREFIX='R3:240.B25D:CO.S1.DACS-1 chisqBinned Core'

  character(len=255)                      :: filename          ! input filename
  character(len=255), dimension(MAXFILES) :: filenames
  integer                                 :: i, error ! Counting indices & Error flags
  type(MLSFile_T)                         :: L2PC1
  type(MLSFile_T)                         :: L2PC2
  integer                                 :: listSize
  integer                                 :: n_filenames
  integer                                 :: NUMSWATHSPERFILE
  integer                                 :: status
  character(len=16)                       :: string
  real                                    :: t1
  real                                    :: t2
  real                                    :: tFile
  logical, parameter                      :: USEALLINPUTSWATHS = .true.
  ! 
! Initialize the lexer, symbol table, and tree checker's tables:
!  ( Under some circumstances, you may need to increase these )
  call init_lexer ( n_chars=80000, n_symbols=4000, hash_table_size=611957 )
  call allocate_decl ( ndecls=8000 )
  call allocate_tree ( n_tree=2000000 )
  call init_tables
  call set_config ( useToolkit = .false., logFileUnit = -1 )
  time_config%use_wall_clock = .true.
  UseSDPToolkit = .false.
  SWITCHES = ' ' ! 'l2pc1,hess'
  CALL mls_h5open(error)
  n_filenames = 0
  do      ! Loop over filenames, options
     call get_filename(filename, options)
     if ( filename(1:1) == '-' ) cycle
     if ( filename == ' ' ) exit
     if ( mls_exists(trim(filename)) /= 0 ) then
       print *, 'Sorry--file not found: ', trim(filename)
       cycle
     endif
     n_filenames = n_filenames + 1
     filenames(n_filenames) = filename
  enddo
  dumpTableSide = options%headSide ! 'left'
  if ( options%debug ) call dumpProgramOptions
  if ( n_filenames < 2 ) then
    if ( options%verbose ) print *, 'Sorry -- need at least 2 input files to diff'
    stop
  endif
  rmsFormat = '(1pe8.1)'
  call time_now ( t1 )
  do i = 2, n_filenames, 2
    call time_now ( tFile )
    if ( options%verbose ) then
      print *, 'diff (1): ', trim(filenames(i-1))
      print *, '     (2): ', trim(filenames(i))
    endif
    status = InitializeMLSFile( L2PC1, type=l_hdf, access=DFACC_RDONLY, &
      & name=trim(filenames(i-1)), HDFVersion=HDFVERSION_5, &
      & shortname='l2pc1' )
    status = InitializeMLSFile( L2PC2, type=l_hdf, access=DFACC_RDONLY, &
      & name=trim(filenames(i)), HDFVersion=HDFVERSION_5, &
      & shortname='l2pc2' )
    if ( .true. ) call diff ( L2PC1, L2PC2, details=options%Details, &
      & options=trim(options%diffOptions) // 'fD' )
    if ( options%timing ) call sayTime('diff these files: ' // trim(filenames(i-1)) // &
      & ' and ' // trim(filenames(i)), tFile)
  enddo
  if ( options%timing ) call sayTime('diffing all files')
  call mls_h5close(error)
contains
!------------------------- get_filenames ---------------------
    subroutine get_filename(filename, options)
    ! Added for command-line processing
     character(len=255), intent(out) :: filename          ! filename
     type ( options_T ), intent(inout) :: options
     ! Local variables
     integer ::                         error = 1
     integer, save ::                   i = 1
     character(len=160)              :: Chars
  ! Get inputfile name, process command-line args
  ! (which always start with -)
    do
      call getarg ( i+hp, filename )
      print *, i, ' th Arg: ', trim(filename)
      error = 0
      if ( filename(1:1) /= '-' ) exit
      if ( filename(1:3) == '-h ' ) then
        call print_help
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
        exit
      else if ( filename(1:3) == '-f ' ) then
        call getarg ( i+1+hp, filename )
        i = i + 1
        exit
      else if ( filename(1:3) == '-d ' ) then
        call getarg ( i+1+hp, Chars )
        ! read(detailChars, '(i)') options%Details
        read(Chars, *) options%Details
        i = i + 1
        exit
      else if ( filename(1:4) == '-deb' ) then
        options%debug = .true.
        exit
      else if ( filename(1:5) == '-side' ) then
        call getarg ( i+1+hp, options%headSide )
        i = i + 1
        exit
      else if ( filename(1:4) == '-opt' ) then
        call getarg ( i+1+hp, options%diffOptions )
        i = i + 1
        exit
      else
        call print_help
      end if
      i = i + 1
    end do
    if ( error /= 0 ) then
      call print_help
    endif
    i = i + 1
    
  end subroutine get_filename
!------------------------- dumpProgramOptions ---------------------
  subroutine dumpProgramOptions
    ! dump options
    call outputNamedValue( 'headSide', trim(options%headSide)         )
    call outputNamedValue( 'Details        ', options%Details         )
    call outputNamedValue( 'diffOptions', trim(options%diffOptions)   )
    call outputNamedValue( 'timing         ', options%timing          )
    call outputNamedValue( 'verbose        ', options%verbose         )
    call outputNamedValue( 'debug          ', options%debug           )
  end subroutine dumpProgramOptions
!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'Usage:l2pcdiff [options] [filenames]'
      write (*,*) &
      & ' diffs l2pc between paired l2pc files [2k] and [2k - 1]'
      write (*,*) &
      & ' each file contains one or more l2pc bins'
      write (*,*) ' Options: -f filename => add filename to list of filenames'
      write (*,*) '                         (can do the same w/o the -f)'
      write (*,*) '          -d details  => level of details to show'
      write (*,*) '          -opt options  => pass options to diff'
      write (*,*) '          -v          => switch on verbose mode'
      write (*,*) '          -debug      => dump options, etc.'
      write (*,*) '          -side "s"   => print stat headers on one of'
      write (*,*) '                          {"top", "left", "right", "bottom"}'
      write (*,*) '          -h          => print brief help'
      write (*,*) '    (Notes)'
      write (*,*) ' (1) Use -opt to assemble list of usual diff options; e.g.'
      write (*,*) '    -rsb to show only rms, statistics, and tables'
      write (*,*) '  character  meaning                                       '
      write (*,*) '     ---     -------                                       '
      write (*,*) '      H       show rank, shape of array                    '
      write (*,*) '      h       diff hessian tensors                         '
      write (*,*) '      j       diff jacobian matrices                       '
      write (*,*) '      b       table of % vs. amount of differences (pdf)   '
      write (*,*) '      c       clean                                        '
      write (*,*) '      g       gaps                                         '
      write (*,*) '      m       mute (print only if different)               '
      write (*,*) '      r       rms -- min, max, etc. of differences         '
      write (*,*) '      s       stats -- number of differences               '
      write (*,*) '      p       transpose                                    '
      write (*,*) '      t       trim                                         '
      write (*,*) '      u       unique                                       '
      write (*,*) '      w       wholearray                                   '
      stop
  end subroutine print_help

!------------------------- print_string ---------------------
  subroutine print_string(string)
    character(len=*), intent(in) :: string
    write(*,'(a)') trim(string)
  end subroutine print_string

!------------------------- SayTime ---------------------
  subroutine SayTime ( What, startTime )
    character(len=*), intent(in) :: What
    real, intent(in), optional :: startTime
    real :: myt1
    if ( present(startTime) ) then
      myt1 = startTime
    else
      myt1 = t1
    endif
    call time_now ( t2 )
    call output ( "Timing for " // what // " = " )
    call output ( dble(t2 - myt1), advance = 'yes' )
  end subroutine SayTime

!==================
end program l2pcdiff
!==================

! $Log$
! Revision 1.4  2015/04/29 00:04:50  pwagner
! Improved -help page
!
! Revision 1.3  2014/01/09 00:31:26  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 1.2  2013/08/23 02:51:48  vsnyder
! Move PrintItOut to PrintIt_m
!
! Revision 1.1  2010/11/25 01:20:32  pwagner
! First commit
!
