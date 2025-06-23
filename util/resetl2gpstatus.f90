! Copyright 2018, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=================================
program resetL2GPStatus ! resets status bits of L2GPData files, e.g. nrt
!=================================

   use BitStuff, only: BitsToBooleans, BooleansToBits
   use Dump_0, only: Dump
   use HDF, only: Dfacc_Rdwr
   use Intrinsic, only: L_Swath
   use L2gpData, only: L2gpData_T, L2gpnamelen, Maxswathnamesbufsize, Rgp, &
     & Appendl2GPData, Destroyl2GPcontents, Dump, &
     & Readl2GPData
   use Machine, only: Hp, Getarg
   use MLSCommon, only: MLSFile_T
   use MLSFiles, only: HDFVersion_5, MLS_Exists, &
     & MLS_Inqswath, InitializeMLSFile
   use MLSHDF5, only: MLS_H5open, MLS_H5close
   use MLSMessagemodule, only: MLSMessageconfig
   use MLSStringlists, only: GetstringElement, &
     & Intersection, NumstringElements
   use Output_M, only: Output
   use Time_M, only: Time_Now, Time_Config
   
   implicit none

!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

! Brief description of program
! Reset Status bits of L2GPData from list of files
! May thus set bits in nrt Temperature file to "Do not use"

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it
! cp /testing/workspace/pwagner/nrttests/v4.23/DL-01/004/outputs/MLS-Aura_L2GP-Temperature_v04-23-nrt-11_2017d111.he5 ./; NAG.Linux-6.2/test -v -profiles 5 10 MLS-Aura_L2GP-Temperature_v04-23-nrt-11_2017d111.he5LF95.Linux/test [options] [input files] -o [output file]

! NoWayManStatus is the default value to which we will reset the Status
! It corresponds to setting these Bit names
! bit    name   value if set
! 0    dontuse     1
! 3  postprocd     8
! Thus                                     1 + 8 = 9
  integer, parameter      :: NoWayManStatus = 9
!
! We'll default to or-ing the current status bits with those set by newstatus
! You can override this by the commandline option -over
  type options_T
    logical               :: overwrite = .false.        ! Overwrite with new status
    logical               :: verbose = .false.
    character(len=255)    :: swathNames = ' '           ! which swaths to reset
    integer               :: newStatus = NoWayManStatus ! reset Status to this
    integer, dimension(2) :: profiles = 0               ! Reset range of profiles   
  end type options_T
  
  type ( options_T ) :: options

  integer, parameter ::          MAXFILES = 100
  logical, parameter ::          countEmpty = .true.
  logical, parameter ::          DEEBUG = .false.
  character(len=255) :: filename          ! input filename
  character(len=255), dimension(MAXFILES) :: filenames
  integer            :: n_filenames
  integer     ::  status, error ! Counting indices & Error flags
  real        :: t1
  real        :: t2
  real        :: tFile
  character(len=L2GPNameLen)          :: swath
  character(len=MAXSWATHNAMESBUFSIZE) :: swathList
  character(len=MAXSWATHNAMESBUFSIZE) :: swathList1
  integer :: listSize
  integer :: NUMSWATHSSOFAR
  ! 
  MLSMessageConfig%useToolkit = .false.
  MLSMessageConfig%logFileUnit = -1
  time_config%use_wall_clock = .true.
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
  if ( n_filenames == 0 ) then
    if ( options%verbose ) print *, 'Sorry no files to reset'
  else
    call time_now ( t1 )

    numswathssofar = 0
    call reset_swaths
  endif
  call mls_h5close(error)
contains
!------------------------- reset_swaths ---------------------
! Identify and reset Status bits

  subroutine reset_swaths
    ! Internal variables
    integer :: i
    integer :: jj
    type(L2GPData_T) :: l2gp
    type( MLSFile_T ) :: l2gpFile
    integer, parameter :: MAXNUMPROFS = 3500
    integer :: numProfs
    integer :: numTotProfs
    real(rgp), dimension(MAXNUMPROFS) :: a
    integer :: N ! size of a
    ! Executable
    a = 0.
    if ( options%verbose ) print *, 'Resetting l2gp file ' // trim(filenames(1))
    numswathssofar = mls_InqSwath ( filenames(1), SwathList, listSize, &
           & hdfVersion=HDFVERSION_5)
    if ( options%swathNames /= ' ' ) then
      swathList1 = swathList
      swathList = Intersection( options%swathNames, swathList1 )
    endif
    call GetStringElement( swathList, swath, 1, countEmpty )
    numTotProfs = 0
    do jj=1, NumStringElements( swathList, countEmpty )
      call GetStringElement( swathList, swath, jj, countEmpty )
      call time_now ( tFile )
      if ( options%verbose ) print *, 'Resetting Status of swath: ', trim(swath)
      N = 0
      do i=1, n_filenames
        status = InitializeMLSFile( l2gpFile, type=l_swath, access=DFACC_RDWR, &
            & content='l2gp', name=trim(filenames(i)), hdfVersion=HDFVERSION_5 )
        if ( options%verbose ) print *, 'Reading from: ', trim(filenames(i))
        call ReadL2GPData( l2gpFile, swath, l2gp, numProfs )
        if ( DEEBUG ) call Dump( l2gp )
        if ( options%overwrite ) then
          call resetL2GPDataByOverwrite ( l2gp )
        else
          call resetL2GPDataByOr ( l2gp )
        endif
        ! l2gpFile%name = trim(filenames(i)) ! options%outputFile
        ! l2gpFile%stillOpen = .true.
        ! call MLS_CloseFile( l2gpFile )
        ! call WriteL2GPData ( l2gp, l2gpFile, swath )
        call AppendL2GPData ( l2gp, l2gpFile, swath )
        N = l2gp%nTimes
        a(1:N) = l2gp%Status
        if ( options%verbose ) call dump( a(1:N), 'Status' )
        call DestroyL2GPContents( l2gp )
        numTotProfs = numTotProfs + N
      enddo
      call sayTime('resetting this swath', tFile)
    enddo
    call sayTime('resetting all swaths')
  end subroutine reset_swaths

  ! These two subroutines reset the Status bits
  ! Either by
  ! (1) overwriting the current value
  subroutine resetL2GPDataByOverwrite ( l2gp )
    ! Args
    type(L2GPData_T) :: l2gp
    ! Executable
    if ( options%profiles(1) > 0 ) then
      l2gp%Status(options%profiles(1):options%profiles(2)) = options%newStatus
    else
      l2gp%Status = options%newStatus
    endif
  end subroutine resetL2GPDataByOverwrite

  ! (2) or-ing the current value with newstatus
  subroutine resetL2GPDataByOr ( l2gp )
    ! Args
    type(L2GPData_T)                :: l2gp
    ! Local variables
    integer, parameter              :: MAXBITNUM = 30
    integer                         :: i
    integer                         :: n1
    integer                         :: n2
    logical, dimension(0:MAXBITNUM) :: newset
    logical, dimension(0:MAXBITNUM) :: set
    ! Executable
    call BitsToBooleans( options%newStatus, newset )
    if ( options%profiles(1) > 0 ) then
      n1 = options%profiles(1)
      n2 = options%profiles(2)
    else
      n1 = 1
      n2 = size(l2gp%Status)
    endif
    do i=n1, n2
      call BitsToBooleans( l2gp%Status(i), set )
      call BooleansToBits ( set .or. newset, l2gp%Status(i) )
    enddo
  end subroutine resetL2GPDataByOr

!------------------------- get_filename ---------------------
    subroutine get_filename(filename, n_filenames, options)
    ! Added for command-line processing
     character(len=255), intent(out) :: filename          ! filename
     integer, intent(in)             :: n_filenames
     type ( options_T ), intent(inout) :: options
     ! Local variables
     integer ::                         error = 1
     integer, save ::                   i = 1
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
      elseif ( filename(1:5) == '-over' ) then
        options%overwrite = .true.
        exit
      else if ( filename(1:3) == '-f ' ) then
        call getarg ( i+1+hp, filename )
        i = i + 1
        exit
      else if ( filename(1:3) == '-s ' ) then
        call getarg ( i+1+hp, options%swathNames )
        i = i + 1
        exit
      else if ( filename(1:3) == '-st' ) then
        call igetarg ( i+1+hp, options%newStatus )
        i = i + 1
        exit
      elseif ( filename(1:5) == '-prof' ) then
        call igetarg ( i+1+hp, options%profiles(1) )
        i = i + 1
        call igetarg ( i+1+hp, options%profiles(2) )
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
    if (trim(filename) == ' ' .and. n_filenames == 0) then

    ! Last chance to enter filename
      print *,  "Enter the name of the HDFEOS4 or 5 L2GP file. " // &
       &  "The default output file name will be used."
      read(*,'(a)') filename
    endif
    
  end subroutine get_filename
!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'Usage: resetl2gpstatus [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options:'
      write (*,*) ' -f filename   => add filename to list of filenames'
      write (*,*) '                  (can do the same w/o the -f)'
      write (*,*) ' -s swaths     => reset only swaths'
      write (*,*) ' -status k     => reset status with bits from k'
      write (*,*) '                   (default is 9)'
      write (*,*) ' -over         => overwrite current status with new'
      write (*,*) '                   (otherwise or bits of current and new)'
      write (*,*) ' -v            => switch on verbose mode'
      write (*,*) ' -profiles m n => reset only profiles in range m n'
      write (*,*) ' -h            => print brief help'
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
!------------------------- igetarg ---------------------
  subroutine igetarg ( pos, iarg )
   integer, intent(in) :: pos
   integer, intent(out) :: iarg
   character(len=16) :: arg
   call getarg ( pos, arg )
   read(arg, *) iarg
  end subroutine igetarg

!==================
end program resetL2GPStatus
!==================

! $Log$
