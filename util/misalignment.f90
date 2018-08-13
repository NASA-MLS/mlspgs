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
program misalignment
!=================================

   use Dump_Options, only: RmsFormat
   use HighOutput, only: OutputNamedValue
   use L2GPData, only: L2GPData_T, MaxSwathNamesBufsize, &
     & DestroyL2GPContents, ReadL2GPData
   use Machine, only: Hp, Getarg
   use MLSCommon, only: DefaultUndefinedValue
   use MLSFiles, only: MLS_Exists, HDFVersion_5, MLS_Inqswath
   use MLSFinds, only: FindFirst
   use MLSHDF5, only: MLS_H5Open, MLS_H5Close
   use MLSHDFEOS, only: MLS_Swath_In_File
   use MLSStrings, only: WriteIntsToChars
   use Output_M, only: Beebeep => Beep, ResumeOutput, SuspendOutput, Output
   use PrintIt_M, only: Set_Config
   use Time_M, only: Time_Now, Time_Config
   
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
! IFC.Linux/test [options] [input files]
  character(len=*), parameter :: swath1 = 'Temperature-InitPtan'
  character(len=*), parameter :: swath2 = 'O3-StdProd' 
  integer, parameter ::          MAXFIELDSLENGTH = 256

  type options_T
    character(len=MAXFIELDSLENGTH) :: fields = '*' ! wild card means 'all'
    character(len=MAXSWATHNAMESBUFSIZE) :: swath1 = swath1
    character(len=MAXSWATHNAMESBUFSIZE) :: swath2 = swath2
    logical     ::         silent          = .false.
    logical     ::         verbose         = .false.
    logical     ::         debug           = .false.
    integer     ::         numDiffs = 0
  end type options_T
  
  type ( options_T ) :: options

  integer, parameter                      :: MAXFILES = 100

  character(len=255)                      :: filename          ! input filename
  character(len=255), dimension(MAXFILES) :: filenames
  integer                                 :: i, error ! Counting indices & Error flags
  integer                                 :: listSize
  integer                                 :: n_filenames
  integer                                 :: numDiffs = 0
  integer                                 :: NUMSWATHSPERFILE
  character(len=16)                       :: string
  character(len=MAXSWATHNAMESBUFSIZE)     :: swathList1
  real                                    :: t1
  real                                    :: t2
  real                                    :: tFile
  logical, parameter                      :: USEALLINPUTSWATHS = .true.
  ! 
  call set_config ( useToolkit = .false., logFileUnit = -1 )
  time_config%use_wall_clock = .true.
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
  if ( options%debug ) call dumpProgramOptions
  if ( options%verbose .and. options%silent ) then
    print *, 'Sorry--either verbose or silent; cant be both'
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
  rmsFormat = '(1pe8.1)'
  if ( options%silent ) call suspendOutput
  call time_now ( t1 )
  do i = 1, n_filenames
    call time_now ( tFile )
    if ( options%verbose ) then
      print *, ' About to check: ', trim(filenames(i))
    endif
    call DiffThese ( filenames(i) )
    options%numDiffs = options%numDiffs + numDiffs
    call sayTime('checking '  // trim(filenames(i)), tFile)
  enddo
  call sayTime('checking all files')
  call resumeOutput
  if ( options%silent .and. options%numDiffs > 0 ) then
    call WriteIntsToChars ( options%numDiffs, string )
    call print_string(string)
  endif
  call mls_h5close(error)
contains
!------------------------- DiffThese ---------------------
    subroutine DiffThese( filename )
    ! Added for command-line processing
      character(len=*), intent(in) :: filename          ! filename
      ! Local variables
      type(L2GPData_T) :: l2gp1
      type(L2GPData_T) :: l2gp2
      logical          :: alignment
      ! Executable
      if ( .not. MLS_swath_in_file( filename, options%swath1, HDFVERSION_5 ) .or. &
        & .not. MLS_swath_in_file( filename, options%swath2, HDFVERSION_5 ) ) return
      call ReadL2GPData( filename, options%swath1, l2gp1 )
      call ReadL2GPData( filename, options%swath2, l2gp2 )
      call Zap (l2gp1, l2gp2 )
      if ( options%verbose ) &
        & print *, 'Checking ' // trim(options%swath1) // &
        &          ' and '     // trim(options%swath2)
      alignment = .true.
      ! Now check for misaligned geolocations
      if ( .not. all( l2gp1%latitude     == l2gp2%latitude    ) ) call beep( 'latitude    ', alignment )
      if ( .not. all( l2gp1%longitude    == l2gp2%longitude   ) ) call beep( 'longitude   ', alignment )
      if ( .not. all( l2gp1%solarTime    == l2gp2%solarTime   ) ) call beep( 'solarTime   ', alignment )
      if ( .not. all( l2gp1%solarZenith  == l2gp2%solarZenith ) ) call beep( 'solarZenith ', alignment )
      if ( .not. all( l2gp1%losAngle     == l2gp2%losAngle    ) ) call beep( 'losAngle    ', alignment )
      if ( .not. all( l2gp1%time         == l2gp2%time        ) ) call beep( 'time        ', alignment )
      ! if ( .not. all( l2gp1%geodAngle    == l2gp2%geodAngle   ) ) call beep( 'geodAngle   ', alignment )
      if ( .not. alignment .and. options%verbose ) then
        print *, 'latitude     ', count( l2gp1%latitude      /= l2gp2%latitude     )
        print *, 'longitude    ', count( l2gp1%longitude     /= l2gp2%longitude    )
        print *, 'solarTime    ', count( l2gp1%solarTime     /= l2gp2%solarTime    )
        print *, 'solarZenith  ', count( l2gp1%solarZenith   /= l2gp2%solarZenith  )
        print *, 'losAngle     ', count( l2gp1%losAngle      /= l2gp2%losAngle     )
        print *, 'geodAngle    ', count( l2gp1%geodAngle     /= l2gp2%geodAngle    )
        print *, 'time         ', count( l2gp1%time          /= l2gp2%time         )
        print *, '1st          ', FindFirst( l2gp1%latitude      /= l2gp2%latitude     )
      endif
      if ( options%silent .and. .not. alignment ) &
        & print *, trim(filename), 'is misaligned'
      if ( options%verbose .and. alignment ) &
        & print *, trim(filename), 'is aligned'
      if ( .not. alignment ) numDiffs = numDiffs + &
        & count( l2gp1%latitude      /= l2gp2%latitude     )
      call destroyL2GPContents ( l2gp1 )
      call destroyL2GPContents ( l2gp2 )
    end subroutine DiffThese

!------------------------- beep ---------------------
    subroutine beep ( geolocation, alignment )
      character(len=*), intent(in) :: geolocation
      logical, intent(out)         :: alignment
      alignment = .false.
      if ( options%silent ) return
      print *, 'filename ', trim(filename)
      print *, 'swath1, swath2 ', trim(swath1), trim(swath2)
      print *, 'geolocation ', trim(geolocation)
    end subroutine beep

!------------------------- Zap ---------------------
    subroutine Zap( l2gp1, l2gp2 )
      ! set to Fill values geolocations associated with invalid profiles
      ! Args
      type(L2GPData_T) :: l2gp1
      type(L2GPData_T) :: l2gp2
      ! Internal variables
      logical :: badProfile
      integer :: profile
      logical :: zapped
      ! Executable
      zapped = .false.
      do profile=1, l2gp1%nTimes
        badProfile = mod( l2gp1%Status(profile), 2 ) /= 0
        badProfile = badProfile .or. mod( l2gp2%Status(profile), 2 ) /= 0
        badProfile = badProfile .or. all( l2gp1%l2gpPrecision(:,:,profile) <= 0. )
        badProfile = badProfile .or. all( l2gp2%l2gpPrecision(:,:,profile) <= 0. )
        if ( .not. badProfile ) cycle
        zapped = .true.
        ! Bad profile, so must zap geolocations
        l2gp1%latitude   (profile)     = defaultUndefinedValue
        l2gp1%longitude  (profile)     = defaultUndefinedValue
        l2gp1%time       (profile)     = defaultUndefinedValue
        l2gp1%losAngle   (profile)     = defaultUndefinedValue
        l2gp1%solarTime  (profile)     = defaultUndefinedValue
        l2gp1%solarZenith(profile)     = defaultUndefinedValue

        l2gp2%latitude   (profile)     = defaultUndefinedValue
        l2gp2%longitude  (profile)     = defaultUndefinedValue
        l2gp2%time       (profile)     = defaultUndefinedValue
        l2gp2%losAngle   (profile)     = defaultUndefinedValue
        l2gp2%solarTime  (profile)     = defaultUndefinedValue
        l2gp2%solarZenith(profile)     = defaultUndefinedValue
      enddo
      ! if ( zapped ) call beebeep( 'Zapped!' )
    end subroutine Zap
!------------------------- get_filename ---------------------
    subroutine get_filename(filename, options)
    ! Added for command-line processing
     character(len=255), intent(out) :: filename          ! filename
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
      elseif ( filename(1:8) == '-silent ' ) then
        options%silent = .true.
        exit
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
        exit
      else if ( filename(1:6) == '-field' ) then
        call getarg ( i+1+hp, options%fields )
        i = i + 1
        exit
      else if ( filename(1:3) == '-f ' ) then
        call getarg ( i+1+hp, filename )
        i = i + 1
        exit
      else if ( filename(1:4) == '-deb' ) then
        options%debug = .true.
        exit
      else if ( filename(1:3) == '-s1' ) then
        call getarg ( i+1+hp, options%swath1 )
        i = i + 1
        exit
      else if ( filename(1:3) == '-s2' ) then
        call getarg ( i+1+hp, options%swath2 )
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
    call outputNamedValue( 'fields', trim(options%fields) )
    call outputNamedValue( 'swath1', trim(options%swath1) )
    call outputNamedValue( 'swath2', trim(options%swath2) )
    call outputNamedValue( 'silent         ', options%silent          )
    call outputNamedValue( 'verbose        ', options%verbose         )
    call outputNamedValue( 'debug          ', options%debug           )
  end subroutine dumpProgramOptions
!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'Usage:misalignment [options] [filenames]'
      write (*,*) &
      & ' Checks for misaligned swaths in l2gp files'
      write (*,*) &
      & ' each file contains one or more swaths, each swath many named fields'
      write (*,*) &
      & ' optionally restrict diffs to certain fields, chunks, etc.'
      write (*,*) ' Options: -f filename => add filename to list of filenames'
      write (*,*) '                  (can do the same w/o the -f)'
      write (*,*) '  -s1 swath1  =>   use swath1 as reference (default is Temperature-InitPtan)'
      write (*,*) '  -s2 swath2  =>   check swaths against ref (default is O3-StdProd)'
      write (*,*) '  -v          => switch on verbose mode'
      write (*,*) '  -silent     => switch on silent mode'
      write (*,*) '                    (printing only if diffs found)'
      write (*,*) '  -debug      => dump options, etc.'
      write (*,*) '  -ignore     => ignore bad chunks'
      write (*,*) '  -h          => print brief help'
      write (*,*) ' (Notes)'
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
end program misalignment
!==================

! $Log$
! Revision 1.3  2016/10/04 22:32:34  pwagner
! Builds properly with some Dumps moved to Dump_1
!
! Revision 1.2  2015/11/17 19:53:30  pwagner
! Fill values put in bad profiles to reduce number of false positives
!
! Revision 1.1  2015/11/03 17:31:00  pwagner
! First commit
!
