! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=================================
PROGRAM L2GPDump ! dumps L2GPData files
!=================================

   use Hdf, only: DFACC_CREATE, DFACC_RDWR, DFACC_READ
   use HDF5, only: h5fopen_f, h5fclose_f, h5fis_hdf5_f   
   use HDFEOS5, only: HE5T_NATIVE_CHAR
   use L2GPData, only: Dump, L2GPData_T, ReadL2GPData, DestroyL2GPContents, &
     & L2GPNameLen, MAXSWATHNAMESBUFSIZE
   use MACHINE, only: FILSEP, HP, IO_ERROR, GETARG
   use MLSCommon, only: R8
   use MLSFiles, only: MLS_IO_GEN_OPENF, MLS_IO_GEN_CLOSEF, &
     & HDFVERSION_4, HDFVERSION_5, MLS_INQSWATH
   use MLSHDF5, only: mls_h5open, mls_h5close
   use MLSMessageModule, only: MLSMessageConfig, MLSMSG_Warning, &
     & MLSMessage
   use MLSStrings, only: GetStringElement, NumStringElements
   use PCFHdr, only: GlobalAttributes
   
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

     integer ::          details = 1
     logical ::          columnsOnly = .false.
     logical ::          attributesToo = .false.
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
     call get_filename(filename, n_filenames, details, columnsOnly, attributesToo)
     if ( filename == ' ' ) exit
     n_filenames = n_filenames + 1
     call h5fis_hdf5_f(trim(filename), is_hdf5, error)
     if ( .not. is_hdf5 ) then
       print *, 'Sorry--not recognized as hdf5 file: ', trim(filename)
     endif
     call dump_one_file(trim(filename), details, columnsOnly, attributesToo)
  enddo
  call mls_h5close(error)
contains
!------------------------- get_filename ---------------------
    subroutine get_filename(filename, n_filenames, details, columnsOnly, attributesToo)
    ! Added for command-line processing
     CHARACTER(LEN=255), intent(out) :: filename          ! filename
     integer, intent(in) ::             n_filenames
     integer, intent(inout) ::          details
     logical, intent(inout) ::          columnsOnly
     logical, intent(inout) ::          attributesToo
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
        details = 0
      elseif ( filename(1:3) == '-1 ' ) then
        details = -1
      elseif ( filename(1:3) == '-2 ' ) then
        details = -2
      elseif ( filename(1:3) == '-a ' ) then
        attributesToo = .true.
      elseif ( filename(1:3) == '-c ' ) then
        columnsOnly = .true.
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
      & 'Usage: l2gpdump [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options: -f filename => use filename'
      write (*,*) '          -h          => print brief help'
      write (*,*) '          -0          => dump only scalars, 1-d array'
      write (*,*) '          -1          => dump only scalars'
      write (*,*) '          -2          => dump only swath names'
      write (*,*) '          -c          => dump only column abundances'
      write (*,*) '          -a          => dump attributes, too'
      stop
  end subroutine print_help

   subroutine dump_one_file(filename, details, columnsOnly, attributesToo)
    character(len=*), intent(in) :: filename          ! filename
    integer, intent(inout) ::          details
    logical, intent(inout) ::          attributesToo
    logical, intent(inout) ::          columnsOnly
    character (len=MAXSWATHNAMESBUFSIZE) :: SwathList
    integer :: File1
    integer :: listsize
    logical, parameter            :: countEmpty = .true.
    type (L2GPData_T) :: l2gp
    integer :: i
    integer :: noSwaths
    character (len=L2GPNameLen) :: swath
    integer :: record_length
    integer :: status
    ! Get swath list
    noSwaths = mls_InqSwath ( filename, SwathList, listSize, &
           & hdfVersion=HDFVERSION_5)
    print *, 'Opening: ', trim(filename)
    File1 = mls_io_gen_openF('swopen', .TRUE., status, &
       & record_length, DFACC_READ, FileName=trim(filename), &
       & hdfVersion=HDFVERSION_5, debugOption=.false. )
    ! print *, 'Status: ', status
    ! Executable code
    noSwaths = NumStringElements(trim(swathList), countEmpty)
    if ( noSwaths < 1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'No swaths to dump--unable to count swaths in ' // trim(swathList) )
    endif
    ! Loop over swaths in file 1
    do i = 1, noSwaths
      call GetStringElement (trim(swathList), swath, i, countEmpty )
      ! Allocate and fill l2gp
      ! print *, 'Reading swath from file: ', trim(swath)
      call ReadL2GPData ( file1, trim(swath), l2gp, &
           & hdfVersion=HDFVERSION_5 )
      ! print *, 'Dumping swath: ', trim(swath)
      ! print *, 'l2gp%nFreqs:  ', l2gp%nFreqs
      ! print *, 'l2gp%nLevels: ', l2gp%nLevels
      ! print *, 'l2gp%nTimes:  ', l2gp%nTimes
      ! print *, 'shape(l2gp%l2gpvalue):  ', shape(l2gp%l2gpvalue)
      ! Dump the swath- and file-level attributes
      if ( attributesToo ) call dump(file1, l2gp)
      ! Dump the actual swath
      call dump(l2gp, columnsOnly, details)
      call DestroyL2GPContents ( l2gp )
    enddo
    status = mls_io_gen_closeF('swclose', File1, FileName=Filename, &
      & hdfVersion=HDFVERSION_5, debugOption=.false.)
   end subroutine dump_one_file
!==================
END PROGRAM L2GPDump
!==================

! $Log$
! Revision 1.3  2004/03/03 19:10:38  pwagner
! Knows to write error messages to stdout
!
! Revision 1.2  2004/02/25 00:07:49  pwagner
! Many options added
!
! Revision 1.1  2004/02/21 00:10:10  pwagner
! First commit
!

