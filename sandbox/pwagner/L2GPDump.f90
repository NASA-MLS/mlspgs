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
PROGRAM L2GPDump ! dumps L2GPData files
!=================================

   use dump_0, only: dump
   use Hdf, only: DFACC_CREATE, DFACC_RDWR, DFACC_READ
   use HDF5, only: h5fopen_f, h5fclose_f, h5gopen_f, h5gclose_f, h5fis_hdf5_f   
   use HDFEOS5, only: HE5T_NATIVE_CHAR
   use intrinsic, only: l_swath
   use L2GPData, only: Dump, L2GPData_T, ReadL2GPData, DestroyL2GPContents, &
     & L2GPNameLen, MAXSWATHNAMESBUFSIZE
   use MACHINE, only: FILSEP, HP, IO_ERROR, GETARG
   use MLSCommon, only: R8
   use MLSFiles, only: HDFVERSION_4, HDFVERSION_5, MLS_INQSWATH, &
     & mls_io_gen_closeF, mls_io_gen_openF, split_path_name
   use MLSHDF5, only: mls_h5open, mls_h5close
   use MLSMessageModule, only: MLSMessageConfig, MLSMSG_Error, MLSMSG_Warning, &
     & MLSMessage
   use MLSStringLists, only: ExpandStringRange, &
     & GetStringElement, NumStringElements, &
     & stringElementNum
   use OUTPUT_M, only: OUTPUT, resumeOutput, suspendOutput
   use PCFHdr, only: GlobalAttributes
   
   IMPLICIT NONE

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This program dumps L2GPData files

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it
! LF95.Linux/test [options] [filenames]

  integer, parameter ::  max_nsds = 1000  ! Maximum number of datasets in file.
  integer, parameter :: MAXNCHUNKS = 50

  type options_T
    character(len=255) ::   chunks = '*' ! wild card means 'all'
     logical ::             verbose = .false.
     logical ::             senseInquiry = .true.
     integer ::             details = 1
     logical ::             columnsOnly = .false.
     logical ::             attributesToo = .false.
     character(len=255) ::  dsInquiry = ''
     character(len=255) ::  attrInquiry = ''
     character(len=255) ::  fields = ''
     character(len=255) ::  swaths = '*' ! wildcard, meaning all swaths
  end type options_T

  type ( options_T ) :: options
  character(LEN=255) :: filename          ! filename
  integer            :: n_filenames
  integer     ::  i, count, status, error ! Counting indices & Error flags
  logical     :: is_hdf5
  logical     :: is_present
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
     if ( options%dsInquiry /= ' ' ) then
       call suspendOutput
       is_present = IsDSInFile( trim(filename), trim(options%dsInquiry) )
       call Respond( options%senseInquiry, is_present, &
         & trim(filename), trim(options%dsInquiry) )
     elseif ( options%attrInquiry /= ' ' ) then
       call suspendOutput
       is_present = IsAttributeInFile( trim(filename), trim(options%attrInquiry) )
       call Respond( options%senseInquiry, is_present, &
         & trim(filename), trim(options%attrInquiry) )
     else
       if ( options%verbose ) print *, 'Dumping swaths in ', trim(filename)
       call dump_one_file(trim(filename), options)
     endif
     call resumeOutput
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
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
      elseif ( filename(1:3) == '-0 ' ) then
        options%details = 0
      elseif ( filename(1:3) == '-1 ' ) then
        options%details = -1
      elseif ( filename(1:3) == '-2 ' ) then
        options%details = -2
      elseif ( filename(1:3) == '-a ' ) then
        options%attributesToo = .true.
      elseif ( filename(1:3) == '-c ' ) then
        options%columnsOnly = .true.
      else if ( filename(1:6) == '-chunk' ) then
        call getarg ( i+1+hp, options%chunks )
        i = i + 1
      else if ( filename(1:6) == '-inqat' ) then
        call getarg ( i+1+hp, options%attrInquiry )
        i = i + 1
      else if ( filename(1:6) == '-inqds' ) then
        call getarg ( i+1+hp, options%dsInquiry )
        i = i + 1
      else if ( filename(1:7) == '-ninqat' ) then
        options%senseInquiry = .false.
        call getarg ( i+1+hp, options%attrInquiry )
        i = i + 1
      else if ( filename(1:7) == '-ninqds' ) then
        options%senseInquiry = .false.
        call getarg ( i+1+hp, options%dsInquiry )
        i = i + 1
      else if ( filename(1:3) == '-l ' ) then
        call getarg ( i+1+hp, options%fields )
        i = i + 1
      else if ( filename(1:3) == '-s ' ) then
        call getarg ( i+1+hp, options%swaths )
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
      & 'Usage: l2gpdump [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) &
      & ' optionally restrict dumps to certain fields, chunks, etc.'
      write (*,*) ' Options: -f filename => use filename'
      write (*,*) '          -h          => print brief help'
      write (*,*) '          -chunks cl  => dump only chunks named in cl'
      write (*,*) '          -[n]inqattr attr'
      write (*,*) '                      => print only if attribute attr [not] present'
      write (*,*) '          -[n]inqds ds'
      write (*,*) '                      => print only if dataset ds [not] present'
      write (*,*) '          -l list     => dump only fields named in list'
      write (*,*) '          -s slist    => dump only swaths named in slist'
      write (*,*) '          -c          => dump only column abundances'
      write (*,*) '          -a          => dump attributes, too'
      write (*,*) '          -v          => verbose'
      write (*,*) '     (details level)'
      write (*,*) '          -0          => dump only scalars, 1-d array'
      write (*,*) '          -1          => dump only scalars'
      write (*,*) '          -2          => dump only swath names'
      write (*,*) '    (Notes)'
      write (*,*) ' (1) by default, dumps all fields in allswaths, but not attributes'
      write (*,*) ' (2) by default, detail level is -1'
      write (*,*) ' (2) details levels, -l options are all mutually exclusive'
      write (*,*) ' (4) the list of chunks may include the range operator "-"'
      stop
  end subroutine print_help
  
  function IsAttributeInFile( file, attribute ) result(sooDesu)
    use MLSHDF5, only: IsHDF5ItemPresent
    use HDF5, only: h5fopen_f, H5F_ACC_RDONLY_F
    ! Dummy args
    character(len=*), intent(in) :: file
    character(len=*), intent(in) :: attribute
    logical :: sooDesu
    ! Local variables
    integer :: fileID
    integer :: grpID
    integer :: status
    character(len=len(attribute)) :: path, name
    ! TRUE if attribute in file
    call h5fopen_f ( trim(file), H5F_ACC_RDONLY_F, fileID, status )
    if ( status /= 0 ) call defeat('Opening file')
    call split_path_name ( attribute, path, name )
    call h5gopen_f( fileID, trim(path), grpID, status )
    if ( status /= 0 ) call defeat('Opening group')
    sooDesu = IsHDF5ItemPresent ( grpID, name, '-a' )
    call h5gclose_f(grpID, status)
    if ( status /= 0 ) call defeat('Closing group')
    call h5fclose_f(fileID, status)
    if ( status /= 0 ) call defeat('Closing file')
  end function IsAttributeInFile

  function IsDSInFile( file, DS ) result(sooDesu)
    use MLSHDF5, only: IsHDF5ItemPresent
    use HDF5, only: h5fopen_f, H5F_ACC_RDONLY_F
    ! Dummy args
    character(len=*), intent(in) :: file
    character(len=*), intent(in) :: DS
    logical :: sooDesu
    ! Local variables
    integer :: fileID
    integer :: status
    integer :: grpID
    character(len=len(DS)) :: path, name
    ! TRUE if DS in file
    call h5fopen_f ( trim(file), H5F_ACC_RDONLY_F, fileID, status )
    call split_path_name ( DS, path, name )
    call h5gopen_f( fileID, trim(path), grpID, status )
    sooDesu = IsHDF5ItemPresent ( grpID, name, '-d' )
    call h5gclose_f(grpID, status)
    call h5fclose_f(fileID, status)
  end function IsDSInFile
  
  subroutine Defeat(msg)
    character(len=*), intent(in) :: msg
    call resumeOutput
    call output('Serious error: ' // msg, advance='yes')
    call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'stopping' )
  end subroutine Defeat
  
  subroutine Respond( sense, test, file, name )
    ! print only if sense matches test
    ! Dummy args
    logical, intent(in)          :: sense
    logical, intent(in)          :: test
    character(len=*), intent(in) :: file
    character(len=*), intent(in) :: name
    character(len=*), parameter  :: Found = 'found'
    character(len=*), parameter  :: notFound = 'not found'
    character(len=16)            :: answer
    ! Executable
    if ( sense .neqv. test ) return
    if ( sense ) then
      answer = found
    else
      answer = notfound
    endif
    call resumeOutput
    call output(trim(name) // ' ' // trim(answer) &
      & // ' in ' // trim(file), advance='yes' )
  end subroutine Respond

   subroutine dump_one_file(filename, options)
    ! Dummy args
    character(len=*), intent(in) :: filename          ! filename
    type ( options_T ) :: options
    ! Local variables
    integer, dimension(MAXNCHUNKS) :: chunks
    character (len=MAXSWATHNAMESBUFSIZE) :: SwathList
    integer :: File1
    integer :: listsize
    logical, parameter            :: countEmpty = .true.
    type (L2GPData_T) :: l2gp
    integer :: i
    integer :: nChunks
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
    ! Loop over swaths in file 1
    do i = 1, noSwaths
      call GetStringElement (trim(swathList), swath, i, countEmpty )
      ! Is this one of the swaths we wished to dump?
      if ( options%swaths /= '*' ) then
        status = stringElementNum(options%swaths, trim(swath), countEmpty)
        if ( status < 1 ) cycle
      endif
      ! Allocate and fill l2gp
      ! print *, 'Reading swath from file: ', trim(swath)
      call ReadL2GPData ( trim(filename), trim(swath), l2gp, &
           & hdfVersion=HDFVERSION_5 )
!       call ReadL2GPData ( file1, trim(swath), l2gp, &
!            & hdfVersion=HDFVERSION_5 )
      ! print *, 'Dumping swath: ', trim(swath)
      ! print *, 'l2gp%nFreqs:  ', l2gp%nFreqs
      ! print *, 'l2gp%nLevels: ', l2gp%nLevels
      ! print *, 'l2gp%nTimes:  ', l2gp%nTimes
      ! print *, 'shape(l2gp%l2gpvalue):  ', shape(l2gp%l2gpvalue)
      ! Dump the swath- and file-level attributes
      if ( options%attributesToo ) then
        File1 = mls_io_gen_openF(l_swath, .TRUE., status, &
         & record_length, DFACC_READ, FileName=trim(filename), &
         & hdfVersion=HDFVERSION_5, debugOption=.false. )
        call dump(file1, l2gp)
        status = mls_io_gen_closeF(l_swath, File1, FileName=Filename, &
        & hdfVersion=HDFVERSION_5, debugOption=.false.)
      endif
      ! Dump the actual swath
      if ( options%verbose ) print *, 'swath: ', trim(swath)
      if ( options%chunks == '*' ) then
        call dump(l2gp, options%columnsOnly, options%details, options%fields)
      else
        call ExpandStringRange(options%chunks, chunks, nchunks)
        if ( nchunks < 1 ) cycle
        call dump(l2gp, chunks(1:nChunks), &
          & options%columnsOnly, options%details, options%fields)
      endif
      call DestroyL2GPContents ( l2gp )
    enddo
!     status = mls_io_gen_closeF('swclose', File1, FileName=Filename, &
!       & hdfVersion=HDFVERSION_5, debugOption=.false.)
   end subroutine dump_one_file
!==================
END PROGRAM L2GPDump
!==================

! $Log$
! Revision 1.9  2006/01/24 23:57:43  pwagner
! May optionally restrict dumps to specific chunks
!
! Revision 1.8  2005/09/27 17:06:16  pwagner
! Added -s option; changed to conform with new MLSFiles module
!
! Revision 1.7  2005/02/11 21:13:59  pwagner
! Now -2 option works correctly
!
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

