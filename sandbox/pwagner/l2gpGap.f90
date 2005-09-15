! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=================================
PROGRAM L2GPGap ! finds geod. angle gap between L2GPData files
!=================================

   use dump_0, only: dump
   use Hdf, only: DFACC_CREATE, DFACC_RDWR, DFACC_READ
   use HDF5, only: h5fopen_f, h5fclose_f, h5fis_hdf5_f   
   use HDFEOS5, only: HE5T_NATIVE_CHAR
   use L2GPData, only: Dump, L2GPData_T, ReadL2GPData, DestroyL2GPContents, &
     & L2GPNameLen, MAXSWATHNAMESBUFSIZE
   use MACHINE, only: FILSEP, HP, IO_ERROR, GETARG
   use MLSCommon, only: R8
   use MLSFiles, only: HDFVERSION_4, HDFVERSION_5, MLS_INQSWATH
   use MLSHDF5, only: mls_h5open, mls_h5close
   use MLSMessageModule, only: MLSMessageConfig, MLSMSG_Warning, &
     & MLSMessage
   use MLSStringLists, only: GetStringElement, NumStringElements
   use OUTPUT_M, only: OUTPUT
   
   IMPLICIT NONE

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This program tests the L2GPData subroutines.

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it
! LF95.Linux/test [options] [filenames]

  type options_T
     logical ::          compare = .false.
     logical ::          verbose = .false.
     ! integer ::          details = 1
     logical ::          columnsOnly = .false.
     character(len=255) ::  fields = ''
     CHARACTER(LEN=255) :: filename1          ! filename
     CHARACTER(LEN=255) :: filename2          ! filename
  end type options_T

  type ( options_T ) :: options
     integer            :: n_filenames
     INTEGER     ::  i, count, status, error ! Counting indices & Error flags
     INTEGER, PARAMETER ::  max_nsds = 1000  ! Maximum number of datasets in file.
     LOGICAL     :: is_hdf5
     ! 
  MLSMessageConfig%useToolkit = .false.   
  MLSMessageConfig%logFileUnit = -1       
  CALL mls_h5open(error)
  do      ! Loop over filenames
     n_filenames = 0
     call get_filename(n_filenames, options)
     if ( options%filename1 == ' ' ) exit
     n_filenames = n_filenames + 1
     call get_filename(n_filenames, options)
     if ( options%filename2 == ' ' ) then
       call print_help
       exit
     endif
     n_filenames = n_filenames + 1
     call h5fis_hdf5_f(trim(options%filename1), is_hdf5, error)
     if ( .not. is_hdf5 ) then
       print *, 'Sorry--not recognized as hdf5 file: ', trim(options%filename1)
     endif
     ! if ( options%verbose ) print *, 'Dumping swaths in ', trim(options%filename1)
     call gap_two_files(options)
  enddo
  call mls_h5close(error)
contains
!------------------------- get_filename ---------------------
    subroutine get_filename(n_filenames, options)
    ! Added for command-line processing
     integer, intent(in) ::             n_filenames
     type ( options_T ) :: options
     integer ::                         error = 1
     integer, save ::                   i = 1
     CHARACTER(LEN=255)              :: filename
  ! Get inputfile name, process command-line args
  ! (which always start with -)
!    i = 1
!    error = 0
    do
      call getarg ( i+hp, filename )
      ! print *, i, ' th Arg: ', trim(filename)
      error = 0
      if ( filename(1:1) /= '-' ) exit
      if ( filename(1:3) == '-h ' ) then
        call print_help
      elseif ( filename(1:3) == '-c ' ) then
        options%compare = .true.
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
      else if ( filename(1:3) == '-l ' ) then
        call getarg ( i+1+hp, options%fields )
        i = i + 1
      else if ( filename(1:3) == '-f ' ) then
        call getarg ( i+1+hp, filename )
        error = 0
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
    ! if (trim(filename) == ' ' .and. n_filenames == 0) call print_help
    if ( n_filenames == 0 ) then
      options%filename1 = filename
    else
      options%filename2 = filename
    endif
    
  end subroutine get_filename
!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'Usage: l2gpGap [options] [filename1 filename2] [..]'
      write (*,*) &
      & 'Displays gaps between successive pairs of l2gp files'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply them'
      write (*,*) ' Options: -f filename => use filename'
      write (*,*) '          -h          => print brief help'
      write (*,*) '          -c          => show comparisons with typical in-file gap'
      write (*,*) '          -v          => verbose'
      write (*,*) ' (by default, gaps times, geod. angles)'
      stop
  end subroutine print_help

   subroutine gap_two_files(options)
    ! Args
    type ( options_T ) :: options
    ! Internal variables
    character (len=MAXSWATHNAMESBUFSIZE) :: SwathList
    integer :: listsize
    logical, parameter            :: countEmpty = .true.
    type (L2GPData_T) :: l2gp1
    type (L2GPData_T) :: l2gp2
    integer :: i
    integer :: noSwaths
    integer :: nTimes
    character (len=L2GPNameLen) :: swath
    integer :: record_length
    integer :: status
    real(r8) :: timeGap
    real(r8) :: angleGap
    ! Executable code
    ! Get swath list
    noSwaths = mls_InqSwath ( options%filename1, SwathList, listSize, &
           & hdfVersion=HDFVERSION_5)
    noSwaths = NumStringElements(trim(swathList), countEmpty)
    if ( noSwaths < 1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'No swaths to dump--unable to count swaths in ' // trim(swathList) )
       return
    endif
    if ( options%verbose ) then
      call output('Gap between files ')
      call output(trim(options%filename1))
      call output(' and ')
      call output(trim(options%filename2), advance='yes')
    endif
    ! Loop over swaths in file 1
    do i = 1, noSwaths
      call GetStringElement (trim(swathList), swath, i, countEmpty )
      ! Allocate and fill l2gp
      ! print *, 'Reading swath from file: ', trim(swath)
      call ReadL2GPData ( options%filename1, trim(swath), l2gp1, &
           & hdfVersion=HDFVERSION_5 )
      call ReadL2GPData ( options%filename2, trim(swath), l2gp2, &
           & hdfVersion=HDFVERSION_5 )
      ! Compute the gaps
      nTimes = l2gp1%nTimes
      timeGap = l2gp2%time(1) - l2gp1%time(nTimes)
      angleGap = l2gp2%geodAngle(1) - l2gp1%geodAngle(nTimes)
      if ( l2gp1%nLevels > 1 ) then
        call showgaps(timeGap, PlusDegrees(angleGap))
        if ( options%compare ) then
          timeGap = l2gp2%time(2) - l2gp2%time(1)
          angleGap = l2gp2%geodAngle(2) - l2gp2%geodAngle(1)
          call showgaps(timeGap, PlusDegrees(angleGap))
        endif
      endif
      call DestroyL2GPContents ( l2gp1 )
      call DestroyL2GPContents ( l2gp2 )
    enddo
   end subroutine gap_two_files

   subroutine showgaps(timeGap, angleGap)
    real(r8) :: timeGap
    real(r8) :: angleGap
     ! Show the time and angle gaps
     call output(timeGap)
     call output(' (s)   ')
     call output(angleGap)
     call output(' (deg)', advance='yes')
   end subroutine showgaps

   function PlusDegrees(arg) result (plus)
    real(r8) :: arg
    real(r8) :: plus
     ! Return 0 < (360*n + arg) < 360 
     ! plus = 360n + arg
     integer :: n
     !
     if ( arg > 0.d0 ) then
       n = -arg/360
     else
       n = -arg/360 + 1
     endif
     plus = 360*n + arg
   end function PlusDegrees

!==================
END PROGRAM L2GPGap
!==================

! $Log$
