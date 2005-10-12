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
PROGRAM L2GPGap ! finds geod. angle gap between L2GPData files
!=================================

   use HDFEOS5, only: he5_swclose, he5_swopen, HE5F_ACC_RDONLY, &
     & he5_EHinqglatts
   use L2GPData, only: L2GPData_T, L2GPNameLen, MAXSWATHNAMESBUFSIZE, &
     & DestroyL2GPContents, ReadL2GPData, SetupNewL2GPRecord
   use MACHINE, only: HP, GETARG
   use MLSCommon, only: R4, R8, MLSFile_T
   use MLSFiles, only: HDFVERSION_5, MLS_INQSWATH
   use MLSHDF5, only: mls_h5open, mls_h5close
   use MLSHDFEOS, only: HE5_EHRDGLATT, MLS_ISGLATT
   use MLSMessageModule, only: MLSMessageConfig, MLSMSG_Error, MLSMSG_Warning, &
     & MLSMessage
   use MLSStringLists, only: GetStringElement, NumStringElements
   use OUTPUT_M, only: OUTPUT
   use HDF5, only: h5fis_hdf5_f   
   
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
     logical ::          useHGrids = .false.
     character(len=255) ::  fields = ''
     CHARACTER(LEN=255) :: filename1          ! filename
     CHARACTER(LEN=255) :: filename2          ! filename
  end type options_T

  type ( options_T ) :: options
     integer            :: n_filenames
     INTEGER     ::  error ! Counting indices & Error flags
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
      elseif ( filename(1:3) == '-H ' ) then
        options%useHGrids = .true.
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
      write (*,*) '          -H          => use HGrids'
      write (*,*) '          -v          => verbose'
      write (*,*) ' (by default, gaps times, geod. angles)'
      stop
  end subroutine print_help

   subroutine gap_two_files(options)
    ! Args
    type ( options_T ) :: options
    ! Internal variables
    real(r8) :: angleGap
    character(len=MAXSWATHNAMESBUFSIZE) :: attrList
    logical, parameter            :: countEmpty = .true.
    integer :: i
    integer, dimension(1) :: ints
    type (L2GPData_T) :: l2gp1
    type (L2GPData_T) :: l2gp2
    type (MLSFile_T) :: L2GPFile1
    type (MLSFile_T) :: L2GPFile2
    integer :: listsize
    integer :: noSwaths
    integer :: nTimes
    integer :: nTimes2
    integer :: status
    character (len=L2GPNameLen) :: swath
    character (len=MAXSWATHNAMESBUFSIZE) :: SwathList
    real(r8) :: timeGap
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
    
    ! Use l2gp data or HGrids?
    if ( options%useHGrids ) then
      L2GPFile1%name = options%filename1
      L2GPFile2%name = options%filename2
      L2GPFile1%FileID%f_id = he5_SWopen( options%filename1, HE5F_ACC_RDONLY )
      L2GPFile2%FileID%f_id = he5_SWopen( options%filename2, HE5F_ACC_RDONLY )
      status = he5_EHinqglatts(L2GPFile1%FileID%f_id, attrList, listSize)
      ! print *, 'status: ', status
      ! print *, 'attrList: ', trim(attrList)
      ! print *, 'listsize: ', listsize
      if ( .not. MLS_ISGLATT( L2GPFile1%FileID%f_id, 'HGrid_noProfs' ) ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'HGrid_noProfs attributes not in the supplied files' )
      else
        ! print *, 'Attribute HGrid_noProfs is there now, before call to HE5_EHRDGLATT'
      endif
      status = he5_EHinqglatts(L2GPFile1%FileID%f_id, attrList, listSize)
      ! print *, 'status: ', status
      ! print *, 'attrList: ', trim(attrList)
      ! print *, 'listsize: ', listsize
      if ( .not. MLS_ISGLATT( L2GPFile1%FileID%f_id, 'HGrid_phi' ) ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'HGrid_phi attributes not in the supplied files' )
      else
        ! print *, 'Attribute HGrid_phi is there now, before call to HE5_EHRDGLATT'
      endif
      status = he5_EHinqglatts(L2GPFile1%FileID%f_id, attrList, listSize)
      ! print *, 'status: ', status
      ! print *, 'attrList: ', trim(attrList)
      ! print *, 'listsize: ', listsize
      status = HE5_EHRDGLATT( L2GPFile1%FileID%f_id, 'HGrid_noProfs', ints )
      nTimes = ints(1)
      status = HE5_EHRDGLATT( L2GPFile2%FileID%f_id, 'HGrid_noProfs', ints )
      nTimes2 = ints(1)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'unable to count profiles in the supplied files' )
      ! print *, 'nTimes : ', nTimes
      ! print *, 'nTimes2: ', nTimes2
      call SetupNewL2GPRecord( l2gp1, nTimes=nTimes  )
      call SetupNewL2GPRecord( l2gp2, nTimes=nTimes2 )
      ! print *, 'Set up both l2gp data for filling'
      if ( .not. MLS_ISGLATT( L2GPFile1%FileID%f_id, 'HGrid_phi' ) ) &
        &  call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'attributes not in the supplied files' )
      status = HE5_EHRDGLATT( L2GPFile1%FileID%f_id, 'HGrid_phi', l2gp1%time )
      l2gp1%geodAngle = real(l2gp1%time, r4)
      ! print *, 'phi  : ', l2gp1%geodAngle(1)
      status = HE5_EHRDGLATT( L2GPFile2%FileID%f_id, 'HGrid_phi', l2gp2%time )
      l2gp2%geodAngle = real(l2gp2%time, r4)
      ! print *, 'phi  : ', l2gp2%geodAngle(1)
      status = HE5_EHRDGLATT( L2GPFile1%FileID%f_id, 'HGrid_time', l2gp1%time )
      ! print *, 'time : ', l2gp1%time(1)
      status = HE5_EHRDGLATT( L2GPFile2%FileID%f_id, 'HGrid_time', l2gp2%time )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'unable to phi, time values in the supplied files' )
      ! print *, 'time : ', l2gp2%time(1)
      timeGap = l2gp2%time(1) - l2gp1%time(nTimes)
      angleGap = l2gp2%geodAngle(1) - l2gp1%geodAngle(nTimes)
      call showgaps(timeGap, PlusDegrees(angleGap))
      if ( options%compare ) then
        timeGap = l2gp2%time(2) - l2gp2%time(1)
        angleGap = l2gp2%geodAngle(2) - l2gp2%geodAngle(1)
        call showgaps(timeGap, PlusDegrees(angleGap))
      endif
      status = he5_swclose( L2GPFile1%FileID%f_id )
      status = he5_swclose( L2GPFile2%FileID%f_id )
      call DestroyL2GPContents ( l2gp1 )
      call DestroyL2GPContents ( l2gp2 )
      return
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
! Revision 1.1  2005/09/15 00:13:20  pwagner
! First commit
!
