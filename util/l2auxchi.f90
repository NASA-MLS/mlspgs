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
program l2auxchi ! dumps chi^sq read from L2AUX files
!=================================

   use Dump_1, only: dump
   use Hdf, only: dfacc_read
   use Hdf5, only: h5fis_hdf5_f, &
     & H5gclose_f, h5gopen_f
   use L1bdata, only: l1bdata_t, &
     & Deallocatel1bdata, readl1bdata
   use Machine, only: hp, getarg
   use MLSCommon, only: namelen
   use MLSKinds, only: r8
   use MLSFiles, only: filenotfound, wildcardhdfversion, &
     & MLS_exists, mls_hdf_version, mls_sfstart, mls_sfend, &
     & HDFVersion_5
   use MLSHdf5, only: getallhdf5dsnames, &
     & MLS_h5open, MLS_H5Close
   use MLSMessageModule, only: MLSMSG_Error, MLSMSG_Warning, &
     & MLSMessage
   use MLSStats1, only: fillValueRelation, stat_t, statistics
   use MLSStringLists, only: catlists, getstringelement, numstringelements
   use MLSStrings, only: lowercase, writeintstochars
   use Output_m, only: blanks, newline, output, resumeoutput, suspendoutput
   use Printit_m, only: set_config
   use Time_m, only: time_now, time_config
   
   implicit none

!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

! Brief description of program
! diffs chi^2 quantities from l2aux files w.r.t. reference value

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"
! Then run it
! LF95.Linux/test [options] [input files]

  type options_T
    logical     :: silent  = .false. ! Not used
    logical     :: tabulate = .false.        ! tabulate
    logical     :: verbose = .false.
    logical     :: list    = .false.
    logical     :: stats   = .true.
    logical     :: rms     = .false.    ! Not used
    integer     :: numDiffs = 0     ! Not used
    real(r8)    :: referenceValue = 1.
  end type options_T
  
  type ( options_T ) :: options


  integer, parameter ::          MAXDS = 500
  integer, parameter ::          MAXSDNAMESBUFSIZE = MAXDS*namelen
  integer, parameter ::          MAXFILES = 100
  character(len=255) :: filename          ! input filename
  character(len=255), dimension(MAXFILES) :: filenames
  integer            :: n_filenames
  integer     ::  i, error ! Counting indices & Error flags
  logical     :: is_hdf5
  character (len=MAXSDNAMESBUFSIZE) :: mySdList
  character (len=MAXSDNAMESBUFSIZE) :: names
  character(len=16) :: string
  real        :: t1
  real        :: t2
  real, dimension( MAXDS, MAXFILES ) :: table
  real        :: tFile
  ! 
  call set_config ( useToolkit = .false., logFileUnit = -1 )
  time_config%use_wall_clock = .true.
  FILLVALUERELATION = '<' ! we want to know % chi^2 are < 1 (or whatever)
  CALL mls_h5open(error)
  n_filenames = 0
  do      ! Loop over filenames
     call get_filename(filename, n_filenames, options)
     if ( filename(1:1) == '-' ) cycle
     if ( filename == ' ' ) exit
     call h5fis_hdf5_f(trim(filename), is_hdf5, error)
     if ( .not. is_hdf5 ) then
       print *, 'Sorry--not recognized as hdf5 file: ', trim(filename)
       cycle
     endif
     n_filenames = n_filenames + 1
     filenames(n_filenames) = filename
  enddo
  if ( n_filenames == 0 ) then
    if ( options%verbose ) print *, 'Sorry no input files supplied'
    stop
  elseif ( options%verbose .and. options%silent ) then
    print *, 'Sorry--either verbose or silent; cant be both'
    stop
  endif
  if ( options%verbose ) call dumpSettings ( options, n_filenames, filenames ) 
  if ( options%silent ) call suspendOutput
  call time_now ( t1 )
  if ( options%verbose .and. .not. options%list ) &
    & print *, 'Compare chi^2 values to: ', options%referenceValue
  do i=1, n_filenames
    call time_now ( tFile )
    if ( options%list ) then
      print *, 'DS Names in: ', trim(filenames(i))
      call GetAllHDF5DSNames (trim(filenames(i)), '/', mysdList)
      call dump(mysdList, 'DS names')
    else 
      if ( options%verbose ) then
        print *, 'get chi^2 from: ', trim(filenames(i))
      endif
      call dumpchisq( trim(filenames(i)), &
      & HDFVERSION_5, options, table(:, i), names )
      call sayTime('reading this file', tFile)
    endif
  enddo
  if ( .not. options%list) call sayTime('dreading all files')
  call resumeOutput
  if ( options%tabulate ) call tabulate ( table, names )
  if ( options%silent .and. options%numDiffs > 0 ) then
    call WriteIntsToChars ( options%numDiffs, string )
    call print_string(string)
  endif
  call mls_h5close(error)
contains
!------------------------- dumpSettings ---------------------
    subroutine dumpSettings( options, n_filenames, filenames )
    ! Added for command-line processing
     integer, intent(in)              :: n_filenames
     character(len=255), dimension(:) :: filenames
     type ( options_T ), intent(in)   :: options
     ! Local variables
     integer :: i
     ! print *, 'silent ?            ', options%silent 
     print *, 'tabulate?           ', options%tabulate
     print *, 'verbose?            ', options%verbose
     print *, 'list   ?            ', options%list   
     print *, 'stats  ?            ', options%stats  
     print *, 'rms    ?            ', options%rms    
     print *, 'referenceValue      ', options%referenceValue
     print *, 'num files           ', n_filenames
     do i=1, n_filenames
       print *, i, trim(filenames(i))
     enddo
    end subroutine dumpSettings

!------------------------- get_filename ---------------------
    subroutine get_filename(filename, n_filenames, options)
    ! Added for command-line processing
     character(LEN=255), intent(out)       :: filename          ! filename
     integer, intent(in)                   :: n_filenames
     type ( options_T ), intent(inout) :: options
     ! Local variables
     integer ::                         error = 1
     integer, save ::                   i = 1
     character(LEN=16)       :: number
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
      elseif ( filename(1:2) == '-t' ) then
        options%tabulate = .true.
        exit
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
        exit
      elseif ( filename(1:3) == '-l ' ) then
        options%list = .true.
        exit
      else if ( filename(1:5) == '-rms ' ) then
        options%rms = .true.
        exit
      else if ( filename(1:3) == '-s ' ) then
        options%stats = .true.
        exit
      else if ( filename(1:3) == '-r ' ) then
        call getarg ( i+1+hp, number )
        read( number, * ) options%referenceValue
        i = i + 1
        exit
      else if ( filename(1:3) == '-f ' ) then
        call getarg ( i+1+hp, filename )
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
      print *,  "Enter the name of the HDF5 l2aux file. "
      read(*,'(a)') filename
    endif
    
  end subroutine get_filename
!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'Usage:l2auxchi [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options: -f filename     => add filename to list of filenames'
      write (*,*) '                           (can do the same w/o the -f)'
      write (*,*) '          -v              => switch on verbose mode'
      write (*,*) '          -silent         => switch on silent mode'
      write (*,*) '                            (printing only if diffs found)'
      write (*,*) '          -l              => just list sd names in files'
      write (*,*) '          -r value        => compare chi^2 against value'
      write (*,*) '          -rms            => just print mean, rms'
      write (*,*) '          -s              => just show % statistics'
      write (*,*) '          -t[abulate]     => print data in tables (dont)'
      write (*,*) '          -h              => print brief help'
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

  ! ---------------------- dumpchisq  ---------------------------
  subroutine dumpchisq( file1, hdfVersion, options, table, names )
  !------------------------------------------------------------------------

    ! Given file name file1,
    ! This routine prints the diffs from referenceValue for chi^2

    ! Arguments

    character (len=*), intent(in) :: file1 ! Name of file
    integer, intent(in)           :: hdfVersion
    type ( options_T )            :: options
    real, dimension( : ) :: table
    character(len=*) :: names

    ! Local
    logical, parameter            :: countEmpty = .true.
    logical :: file_exists
    integer :: grpid
    integer :: i
    type(l1bdata_t) :: L1BDATA  ! Result
    type(Stat_T) :: L1BStat
    character (len=MAXSDNAMESBUFSIZE) :: mySdList
    integer :: NoMAFs
    integer :: noSds
    integer :: sdfid1
    character (len=80) :: sdName
    integer, dimension(3) :: shp
    integer :: status
    integer :: the_hdfVersion
    
    ! Executable code
    the_hdfVersion = HDFVERSION_5
    the_hdfVersion = hdfVersion
    file_exists = ( mls_exists(trim(File1)) == 0 )
    if ( .not. file_exists ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'File 1 not found; make sure the name and path are correct' &
        & // trim(file1) )
    endif
    if ( the_hdfVersion == WILDCARDHDFVERSION ) then
      the_hdfVersion = mls_hdf_version(File1, hdfVersion)
      if ( the_hdfVersion == FILENOTFOUND ) &
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'File 1 not found; make sure the name and path are correct' &
          & // trim(file1) )
    endif
    call GetAllHDF5DSNames (trim(File1), '/', mysdList)
    if ( options%verbose) then
      call output ( '============ DS names in ', advance='no' )
      call output ( trim(file1) //' ============', advance='yes' )
    endif
    if ( mysdList == ' ' ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'No way yet to find sdList in ' // trim(File1) )
      return
    else
      if ( options%verbose ) call dump(mysdList, 'DS names')
    endif

    sdfid1 = mls_sfstart(File1, DFACC_READ, hdfVersion=hdfVersion)
    if (sdfid1 == -1 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      &  'Failed to open l1b file ' // trim(File1) )
    end if
	 call h5gOpen_f (sdfid1,'/', grpID, status)
    if ( status /= 0 ) then
	   call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	& 'Unable to open group to read attribute in l2aux file' )
    endif
    noSds = NumStringElements(trim(mysdList), countEmpty)
    if ( noSds < 1 ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'No sdNames cp to file--unable to count sdNames in ' // trim(mysdList) )
    endif
    ! Loop over sdNames in file 1
    names = ''
    do i = 1, noSds
      call GetStringElement (trim(mysdList), sdName, i, countEmpty )
      if ( index( lowercase(trim(sdName)), 'chisq' ) < 1 ) cycle
      ! Allocate and fill l2aux
      if ( options%verbose ) print *, 'About to read ', trim(sdName)
      call ReadL1BData ( sdfid1, trim(sdName), L1bData, NoMAFs, status, &
        & hdfVersion=the_hdfVersion, NEVERFAIL=.true., L2AUX=.true. )
      if ( status /= 0 ) then
	     call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Unable to find ' // trim(sdName) // ' in ' // trim(File1) )
        cycle
      endif
      shp = shape(L1bData%DpField)
      call dump( L1bData%DpField, name=trim(sdName), &
        & options="-rs", FillValue=options%referenceValue )
      L1BStat%count = 0
      call statistics(&
        & reshape( L1bData%DpField, (/ product(shp) /) ), &
        & L1BStat, fillValue=options%referenceValue )
      call DeallocateL1BData ( l1bData )
      table( i ) = L1BStat%fillcount * 1.d0/L1BStat%count
      options%numDiffs = options%numDiffs + L1BStat%fillcount
      names = catlists(names, sdName)
    enddo
	 call h5gClose_f (grpID, status)
    if ( status /= 0 ) then
	   call MLSMessage ( MLSMSG_Warning, ModuleName, &
       & 'Unable to close group in l2aux file: ' // trim(File1) // ' after diffing' )
    endif
	 status = mls_sfend(sdfid1, hdfVersion=the_hdfVersion)
    if ( status /= 0 ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
       & "Unable to close L2aux file: " // trim(File1) // ' after diffing')
  end subroutine dumpchisq

!------------------------- print_string ---------------------
  subroutine print_string(string)
    character(len=*), intent(in) :: string
    write(*,'(a)') trim(string)
  end subroutine print_string

!------------------------- tabulate ---------------------
  subroutine tabulate ( table, phases )
    ! Args
    real, dimension(:,:), intent(in) :: table
    character(len=*), optional, intent(in) :: phases
    ! Internal variables
    integer :: chunk
    character(len=32) :: myPhase
    integer :: numPhases
    integer :: phase
    ! Executable
    numPhases = size(table, 1)
    call blanks(30, fillchar='*', advance='yes')
    call output('chunk', advance='no')
    call blanks(3)
    do phase=1, numPhases
      if ( present(phases) ) then
        call GetStringElement( phases, myPhase, phase, .FALSE.)
        call output(trim(myPhase), advance='no')
        call blanks(3)
      else
        call output('phase', advance='no')
        call output(phase, advance='no')
        call blanks(3)
      endif
    enddo
    call newline
    do chunk=1, size(table, 2)
      if ( any(table(:, chunk) /= options%referenceValue) ) then
        call output(chunk, advance='no')
        call blanks(1)
        call output(table(:, chunk), advance='no')
        call newline
      endif
    enddo
  end subroutine tabulate

!==================
end program l2auxchi
!==================

! $Log$
! Revision 1.7  2014/03/07 21:44:11  pwagner
! Name_Len changed to nameLen; got from MLSCommon
!
! Revision 1.6  2013/08/23 02:51:47  vsnyder
! Move PrintItOut to PrintIt_m
!
! Revision 1.5  2010/06/09 18:14:14  pwagner
! Made compatible with new dump api
!
! Revision 1.4  2009/04/13 20:43:17  pwagner
! Fixed a bug preventing macros file from using its own macros properly
!
! Revision 1.3  2006/08/23 00:00:37  pwagner
! Many untested changes; compiles successfully
!
! Revision 1.2  2006/05/24 22:22:05  pwagner
! Increased MAXDS to more realistic value
!
! Revision 1.1  2006/05/24 20:39:41  pwagner
! First commit
!
