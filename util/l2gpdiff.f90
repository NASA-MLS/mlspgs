! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=================================
program l2gpdiff ! show diffs between swaths in two different files
!=================================

   use Hdf, only: DFACC_CREATE, DFACC_RDWR, DFACC_READ
   use HDF5, only: h5fopen_f, h5fclose_f, h5fis_hdf5_f   
   use HDFEOS5, only: HE5T_NATIVE_CHAR
   use L2GPData, only: L2GPData_T, &
     & Diff, &
     & L2GPNameLen, MAXSWATHNAMESBUFSIZE
   use MACHINE, only: FILSEP, HP, IO_ERROR, GETARG
   use MLSCommon, only: R8
   use MLSFiles, only: mls_exists, MLS_IO_GEN_OPENF, MLS_IO_GEN_CLOSEF, &
     & HDFVERSION_4, HDFVERSION_5, MLS_INQSWATH
   use MLSHDF5, only: mls_h5open, mls_h5close
   use MLSMessageModule, only: MLSMessageConfig
   use MLSStrings, only: GetStringElement, NumStringElements
   use output_m, only: output
   use PCFHdr, only: GlobalAttributes
   use Time_M, only: Time_Now, USE_WALL_CLOCK
   
   implicit none

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it
! LF95.Linux/test [options] [input files]
  type options_T
    logical     :: verbose = .false.
    integer     :: Details = 1
    logical     :: stats = .false.
    logical     :: rms = .false.
    logical     :: ignoreBadChunks = .false.
    logical     :: showMissing = .false.
  end type options_T
  
  type ( options_T ) :: options

  integer, parameter ::          MAXFILES = 100
  ! logical ::          columnsOnly
  character(len=255) :: filename          ! input filename
  ! character(len=255) :: outputFile= 'default.he5'        ! output filename
  character(len=255), dimension(MAXFILES) :: filenames
  integer            :: n_filenames
  integer     ::  i, j, count, status, error ! Counting indices & Error flags
  integer     ::  hdfversion1
  integer     ::  hdfversion2
  logical     :: is_hdf5
  character(len=MAXSWATHNAMESBUFSIZE) :: swathList1
  character(len=MAXSWATHNAMESBUFSIZE) :: swathList2
  character(len=10) :: swath_name
  character(len=10) :: column_name
  character(len=*), parameter :: SWATHPREFIX='R3:240.B25D:CO.S1.DACS-1 chisqBinned Core'
  logical, parameter :: USEALLINPUTSWATHS = .true.
  integer :: NUMSWATHSPERFILE
  integer :: listSize
  ! logical     :: verbose = .false.
  real        :: t1
  real        :: t2
  real        :: tFile
  ! 
  MLSMessageConfig%useToolkit = .false.
  MLSMessageConfig%logFileUnit = -1
  USE_WALL_CLOCK = .true.
  CALL mls_h5open(error)
  n_filenames = 0
  do      ! Loop over filenames
     call get_filename(filename, n_filenames, options)
     if ( filename(1:1) == '-' ) cycle
     if ( filename == ' ' ) exit
     if ( mls_exists(trim(filename)) /= 0 ) then
       print *, 'Sorry--file not found: ', trim(filename)
       cycle
     endif
     n_filenames = n_filenames + 1
     filenames(n_filenames) = filename
  enddo
  if ( n_filenames < 2 ) then
    if ( options%verbose ) print *, 'Sorry -- need at least 2 input files to diff'
    stop
  endif
  if ( USEALLINPUTSWATHS ) then
    NUMSWATHSPERFILE = mls_InqSwath ( trim(filenames(1)), &
      & swathList1, listSize, hdfVersion=HDFVERSION_5)
    if ( NUMSWATHSPERFILE < 1 ) then
      print *, ' NUMSWATHSPERFILE: ', NUMSWATHSPERFILE
      print *, ' file1: ', trim(filenames(1))
      print *, ' swathList1: ', trim(swathList1)
      print *, ' listSize: ', listSize
      stop
    endif
  else
    NUMSWATHSPERFILE = 2
    swathList1 = ''
  endif
  call time_now ( t1 )
  count = 0
  do i = 2, n_filenames, 2
    call time_now ( tFile )
    if ( options%verbose ) then
      print *, 'diff (1): ', trim(filenames(i-1))
      print *, '     (2): ', trim(filenames(i))
    endif
    call diff(trim(filenames(i-1)), &
      & trim(filenames(i)), details=options%Details, stats=options%stats, &
      & rms=options%rms, ignoreBadChunks=options%ignoreBadChunks, &
      & showMissing=options%showMissing)
    call sayTime('diff these files: ' // trim(filenames(i-1)) // &
      & ' and ' // trim(filenames(i)), tFile)
  enddo
  call sayTime('diffing all files')
  call mls_h5close(error)
contains
!------------------------- get_filenames ---------------------
    subroutine get_filename(filename, n_filenames, options)
    ! Added for command-line processing
     character(LEN=255), intent(out) :: filename          ! filename
     integer, intent(in)             :: n_filenames
     type ( options_T ), intent(inout) :: options
     ! character(LEN=*), intent(inout) :: outputFile        ! output filename
     ! logical, intent(inout)          :: verbose
     ! Local variables
     integer ::                         error = 1
     integer, save ::                   i = 1
     character(LEN=16)              :: detailChars
  ! Get inputfile name, process command-line args
  ! (which always start with -)
    do
      call getarg ( i+hp, filename )
      ! print *, i, ' th Arg: ', trim(filename)
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
        call getarg ( i+1+hp, detailChars )
        ! read(detailChars, '(i)') options%Details
        read(detailChars, *) options%Details
        i = i + 1
        exit
      else if ( filename(1:4) == '-ign' ) then
        options%ignorebadchunks = .true.
        exit
      else if ( filename(1:5) == '-rms ' ) then
        options%rms = .true.
        exit
      else if ( filename(1:3) == '-s ' ) then
        options%stats = .true.
        exit
      else if ( filename(1:6) == '-miss ' ) then
        options%showMissing = .true.
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
!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'Usage:l2gpcat [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options: -f filename => add filename to list of filenames'
      write (*,*) '                         (can do the same w/o the -f)'
      write (*,*) '          -d details  => level of details to show'
      write (*,*) '          -v          => switch on verbose mode'
      write (*,*) '          -ignore     => ignore bad chunks'
      write (*,*) '          -rms        => just print mean, rms'
      write (*,*) '          -s          => just show statistics'
      write (*,*) '          -miss       => just show which swaths are missing'
      write (*,*) '          -h          => print brief help'
      stop
  end subroutine print_help
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
end program l2gpdiff
!==================

! $Log$
! Revision 1.1  2004/06/16 17:55:06  pwagner
! First commit
!
