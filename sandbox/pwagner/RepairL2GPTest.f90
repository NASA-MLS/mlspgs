! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=================================
PROGRAM RepairL2GPTest ! dumps L2GPData files
!=================================

   use dump_0, only: dump
   use Hdf, only: DFACC_CREATE, DFACC_RDWR, DFACC_READ
   use HDF5, only: h5fopen_f, h5fclose_f, h5fis_hdf5_f   
   use HDFEOS5, only: HE5T_NATIVE_CHAR, HE5F_ACC_RDONLY, &
     & he5_swclose, he5_swopen
   use L2GPData, only: L2GPData_T, &
     & Diff, Dump, ReadL2GPData, DestroyL2GPContents, &
     & L2GPNameLen, MAXSWATHNAMESBUFSIZE, RepairL2GP, setupnewl2gprecord
   use MACHINE, only: FILSEP, HP, IO_ERROR, GETARG
   use MLSCommon, only: DEFAULTUNDEFINEDVALUE, R4, R8
   use MLSFiles, only: MLS_IO_GEN_OPENF, MLS_IO_GEN_CLOSEF, &
     & HDFVERSION_4, HDFVERSION_5, MLS_INQSWATH
   use MLSHDF5, only: mls_h5open, mls_h5close
   use MLSMessageModule, only: MLSMessageConfig, MLSMSG_Warning, &
     & MLSMessage
   use MLSNumerics, only: ReplaceFillValues
   use MLSStringLists, only: ExpandStringRange, GetStringElement, &
     & NumStringElements
   use OUTPUT_M, only: OUTPUT
   use PCFHdr, only: GlobalAttributes
   
   implicit none

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
  logical, parameter :: IGNORE = .false.
  integer, parameter :: MAXNCHUNKS = 400
  integer, parameter :: MAXNTIMES = 4000
  type options_T
     integer, dimension(MAXNCHUNKS) ::chunks = 0 ! which chunks crashed
     character(len=255) ::            chunkRange = '1' ! e.g., '1,4-90'
     logical, dimension(MAXNTIMES) :: crashed = .false. ! which profiles affected
     integer ::                       details = 1
     character(len=255) ::            fields = '*' ! which fields to munge
     ! logical ::                     ignorebadchunks = .false.
     logical ::                       rms = .false.
     logical ::                       stats = .false.
     logical ::                       verbose = .false.
     logical ::                       wholeArray = .false.
  end type options_T

  type ( options_T ) :: options
     CHARACTER(LEN=255) :: filename          ! filename
     integer            :: n_filenames
     INTEGER     ::  i, count, status, error ! Counting indices & Error flags
     INTEGER, PARAMETER ::  max_nsds = 1000  ! Maximum number of datasets in file.
     LOGICAL     :: is_hdf5
     ! 
  MLSMessageConfig%useToolkit = .false.   
  MLSMessageConfig%logFileUnit = -1       
  CALL mls_h5open(error)
  n_filenames = 0
  do      ! Loop over filenames
     call get_filename(filename, n_filenames, options)
     if ( filename == ' ' ) exit
     n_filenames = n_filenames + 1
     call h5fis_hdf5_f(trim(filename), is_hdf5, error)
     if ( .not. is_hdf5 ) then
       print *, 'Sorry--not recognized as hdf5 file: ', trim(filename)
     endif
     call ExpandStringRange(options%chunkRange, options%chunks)
     if ( options%verbose ) then
       call dump(options%chunks, 'crashed chunks')
       print *, 'Munging and repairing swaths in ', &
         & trim(filename)
     endif
     call dump_one_file(trim(filename), options)
  enddo
  call mls_h5close(error)
contains
!------------------------- get_filename ---------------------
    subroutine get_filename(filename, n_filenames, options)
    ! Added for command-line processing
     CHARACTER(LEN=255), intent(out) :: filename          ! filename
     integer, intent(in) ::             n_filenames
     type ( options_T ) :: options
     integer ::                         error = 1
     integer, save ::                   i = 1
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
      elseif ( filename(1:3) == '-0 ' ) then
        options%details = 0
      elseif ( filename(1:3) == '-1 ' ) then
        options%details = -1
      elseif ( filename(1:3) == '-2 ' ) then
        options%details = -2
      ! else if ( filename(1:4) == '-ign' ) then
        ! options%ignorebadchunks = .true.
      else if ( filename(1:5) == '-rms ' ) then
        options%rms = .true.
      else if ( filename(1:3) == '-s ' ) then
        options%stats = .true.
      else if ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
      else if ( filename(1:2) == '-c' ) then
        call getarg ( i+1+hp, options%chunkRange )
        i = i + 1
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
    if (trim(filename) == ' ' .and. n_filenames == 0) then

    ! Last chance to enter filename
      print *,  "Enter the name of the HDF5 file. " // &
       &  "Datasets in the file will be listed shortly."
      read(*,'(a)') filename
    endif
    
  end subroutine get_filename
!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'Usage: RepairL2GPTest [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options: -f filename => use filename'
      write (*,*) '          -h          => print brief help'
      write (*,*) '          -c chunks   => munge chunks'
      write (*,*) '                          e.g., "1,4,8-12"'
      write (*,*) '          -l list     => dump only fields named in list'
      write (*,*) '          -0          => dump only scalars, 1-d array'
      write (*,*) '          -1          => dump only scalars'
      write (*,*) '          -2          => dump only swath names'
      write (*,*) '          -ignore     => ignore bad chunks'
      write (*,*) '          -rms        => just print mean, rms'
      write (*,*) '          -s          => just show statistics'
      write (*,*) '          -v          => verbose'
      write (*,*) ' (by default, dumps all fields in allswaths, but not attributes)'
      stop
  end subroutine print_help

   subroutine dump_one_file(filename, options)
    character(len=*), intent(in) :: filename          ! filename
    type ( options_T ) :: options
    character (len=MAXSWATHNAMESBUFSIZE) :: SwathList
    integer :: File1
    integer :: itime
    integer :: listsize
    logical, parameter            :: countEmpty = .true.
    type (L2GPData_T) :: l2gp1
    type (L2GPData_T) :: l2gp2
    integer :: i
    integer :: noSwaths
    character (len=L2GPNameLen) :: swath
    integer :: record_length
    integer :: status
    ! Get swath list
    noSwaths = mls_InqSwath ( filename, SwathList, listSize, &
           & hdfVersion=HDFVERSION_5)
    if ( options%details < -1 ) then
      call output('swaths in ' // trim(filename), advance='yes')
      call dump(trim(swathList))
      return
    endif
    ! print *, 'Opening: ', trim(filename)
    file1 = he5_swopen(trim(fileName), HE5F_ACC_RDONLY)
!     File1 = mls_io_gen_openF('swopen', .TRUE., status, &
!        & record_length, DFACC_READ, FileName=trim(filename), &
!        & hdfVersion=HDFVERSION_5, debugOption=.false. )
    ! print *, 'Status: ', status
    ! Executable code
    noSwaths = NumStringElements(trim(swathList), countEmpty)
    if ( noSwaths < 1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'No swaths to dump--unable to count swaths in ' // trim(swathList) )
    endif
    ! Loop over swaths in file
    do i = 1, noSwaths
      call GetStringElement (trim(swathList), swath, i, countEmpty )
      ! Allocate and fill l2gp
      ! print *, 'Reading swath from file: ', trim(swath)
      call ReadL2GPData ( file1, trim(swath), l2gp1, &
           & hdfVersion=HDFVERSION_5 )
      ! Reinitialize which chunks crashed
      options%crashed = .false.
      do itime=1, l2gp1%nTimes
        options%crashed(itime) = any(l2gp1%chunkNumber(itime) == options%chunks)
      enddo
      call CloneL2GPData(l2gp1, l2gp2)
      call MungeL2GPData(l2gp1)
      ! Diff the two swaths
      print *, 'Before Repair'
      call diff(l2gp1, l2gp2, options%details, options%wholeArray, &
        & options%stats, options%rms, IGNORE, options%fields)
      ! Repair the munged swath
      call RepairL2GP ( L2GP1, L2GP2, options%fields )
      ! Diff again
      print *, 'After repair'
      call diff(l2gp1, l2gp2, options%details, options%wholeArray, &
        & options%stats, options%rms, IGNORE, options%fields)
      call DestroyL2GPContents ( l2gp1 )
      call DestroyL2GPContents ( l2gp2 )
    enddo
!     status = mls_io_gen_closeF('swclose', File1, FileName=Filename, &
!       & hdfVersion=HDFVERSION_5, debugOption=.false.)
     status = he5_swclose(File1)
   end subroutine dump_one_file

   subroutine CloneL2GPData(l2gp1, l2gp2)
     ! Args
     type(L2GPData_T), intent(in)  :: l2gp1
     type(L2GPData_T), intent(out) :: l2gp2
     ! Internal variables
     ! Executable
     call SetupNewL2GPRecord ( l2gp2, l2gp1%nFreqs, l2gp1%nLevels, &
       & l2gp1%nTimes, l2gp1%nTimesTotal, &
       & FillIn=.true.)
     l2gp2%pressures    = l2gp1%pressures
     l2gp2%latitude     = l2gp1%latitude 
     l2gp2%longitude    = l2gp1%longitude
     l2gp2%frequency    = l2gp1%frequency
     l2gp2%solarTime    = l2gp1%solarTime
     l2gp2%solarZenith  = l2gp1%solarZenith
     l2gp2%losAngle     = l2gp1%losAngle 
     l2gp2%geodAngle    = l2gp1%geodAngle
     l2gp2%time         = l2gp1%time     
     l2gp2%chunkNumber  = l2gp1%chunkNumber
     l2gp2%l2gpvalue    = l2gp1%l2gpvalue     
     l2gp2%l2gpprecision= l2gp1%l2gpprecision 
     l2gp2%status       = l2gp1%status        
     l2gp2%quality      = l2gp1%quality       
   end subroutine CloneL2GPData
   
   subroutine MungeL2GPData(l2gp)
     ! Args
     type(L2GPData_T), intent(inout)  :: l2gp
     ! Internal variables
     ! Executable
     call Munge1dr4(l2gp%pressures    )
     call Munge1dr4(l2gp%latitude     )
     call Munge1dr4(l2gp%longitude    )
     call Munge1dr4(l2gp%frequency    )
     call Munge1dr4(l2gp%solarTime    )
     call Munge1dr4(l2gp%solarZenith  )
     call Munge1dr4(l2gp%losAngle     )
     call Munge1dr4(l2gp%geodAngle    )
     call Munge1dr8(l2gp%time         )
     call Munge1dint(l2gp%chunkNumber  )
     call Munge3dr4(l2gp%l2gpvalue    )
     call Munge3dr4(l2gp%l2gpprecision)
     call Munge1dint(l2gp%status       )
     call Munge1dr4(l2gp%quality      )
   end subroutine MungeL2GPData
   
   subroutine Munge1dr4(array)
     ! Args
     real(r4), dimension(:), intent(inout) :: array
     ! Internal variables
     integer :: n
     ! Executable
     n = size(array)
     where (options%crashed(1:n))
       array = DEFAULTUNDEFINEDVALUE
     end where
   end subroutine Munge1dr4
   
   subroutine Munge1dr8(array)
     ! Args
     real(r8), dimension(:), intent(inout) :: array
     ! Internal variables
     integer :: n
     ! Executable
     n = size(array)
     where (options%crashed(1:n))
       array = DEFAULTUNDEFINEDVALUE
     end where
   end subroutine Munge1dr8
   
   subroutine Munge1dint(array)
     ! Args
     integer, dimension(:), intent(inout) :: array
     ! Internal variables
     integer :: n
     ! Executable
     n = size(array)
     where (options%crashed(1:n))
       array = 513 ! 512 + 1
     end where
   end subroutine Munge1dint
   
   subroutine Munge3dr4(array)
     ! Args
     real(r4), dimension(:,:,:), intent(inout) :: array
     ! Internal variables
     integer :: i
     ! Executable
     do i=1, size(array, 3)
       if ( options%crashed(i) ) array(:,:,i) = DEFAULTUNDEFINEDVALUE
     end do
   end subroutine Munge3dr4
   
!==================
END PROGRAM RepairL2GPTest
!==================

! $Log$
! Revision 1.6  2004/10/27 16:48:55  pwagner
! -l list option added to specify which fields to dump
!
! Revision 1.5  2004/09/28 23:10:15  pwagner
! Much moved to MLSStringLists
!
! Revision 1.4  2004/07/20 22:42:40  pwagner
! Works with newer toolkit 5.2.11
!
! Revision 1.3  2004/03/03 19:10:38  pwagner
! Knows to write error messages to stdout
!
! Revision 1.2  2004/02/25 00:07:49  pwagner
! Many options added
!
! Revision 1.1  2004/02/21 00:10:10  pwagner
! First commit
!

