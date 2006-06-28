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
program l2auxdump ! dumps datasets, attributes from L2AUX files
!=================================

   use Dump_0, only: DUMP, RELATIONFORPCTAGES
   use Hdf, only: DFACC_CREATE, DFACC_READ
   use HDF5, only: h5fis_hdf5_f, &
     & H5GCLOSE_F, H5GOPEN_F, H5DOPEN_F, H5DCLOSE_F, h5gcreate_f
   use L1BData, only: L1BData_T, NAME_LEN, &
     & DeallocateL1BData, Diff, ReadL1BData
   use MACHINE, only: FILSEP, HP, IO_ERROR, GETARG
   use MLSCommon, only: R8
   use MLSFiles, only: FILENOTFOUND, WILDCARDHDFVERSION, &
     & mls_exists, mls_hdf_version, mls_sfstart, mls_sfend, &
     & HDFVERSION_4, HDFVERSION_5
   use MLSHDF5, only: DumpHDF5Attributes, DumpHDF5DS, &
     & GetAllHDF5AttrNames, GetAllHDF5DSNames, &
     & IsHDF5AttributePresent, mls_h5open, mls_h5close
   use MLSMessageModule, only: MLSMessageConfig, MLSMSG_Error, MLSMSG_Warning, &
     & MLSMessage
   use MLSStringLists, only: GetStringElement, NumStringElements
   use MLSStrings, only: LOWERCASE, WriteIntsToChars
   use output_m, only: resumeOutput, suspendOutput, output
   use Time_M, only: Time_Now, time_config
   
   implicit none

!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

! Brief description of program
! dumps datasets or attributes from l2aux files

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"
! Then run it
! LF95.Linux/test [options] [input files]

  type options_T
    logical     :: verbose            = .false.
    logical     :: la                 = .false.
    logical     :: ls                 = .false.
    logical     :: stats              = .false.
    logical     :: rms                = .false.
    character(len=128) :: DSName      = '' ! Extra dataset if attributes under one
    character(len=128) :: root        = '/'
    character(len=1024) :: attributes = ''
    character(len=1024) :: datasets   = '*'
  end type options_T
  
  type ( options_T ) :: options


  integer, parameter ::          MAXDS = 500
  integer, parameter ::          MAXSDNAMESBUFSIZE = MAXDS*NAME_LEN
  integer, parameter ::          MAXFILES = 100
  integer, parameter ::          hdfVersion = HDFVERSION_5
  character(len=255) :: filename          ! input filename
  character(len=255), dimension(MAXFILES) :: filenames
  integer            :: n_filenames
  integer     ::  i, count, status, error ! Counting indices & Error flags
  logical     :: is_hdf5
  character (len=MAXSDNAMESBUFSIZE) :: mySdList
  integer     :: sdfid1
  character(len=16) :: string
  real        :: t1
  real        :: t2
  real        :: tFile
  ! 
  MLSMessageConfig%useToolkit = .false.
  MLSMessageConfig%logFileUnit = -1
  time_config%use_wall_clock = .true.
  relationforpctages = '<' ! we want to know % chi^2 are < 1 (or whatever)
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
  endif
  if ( options%verbose ) call dumpSettings ( options, n_filenames, filenames ) 
  call time_now ( t1 )
  do i=1, n_filenames
    call time_now ( tFile )
    if ( options%ls ) then
      print *, 'DS Names in: ', trim(filenames(i)), ' ', trim(options%root)
      call GetAllHDF5DSNames (trim(filenames(i)), trim(options%root), mysdList)
      call dump(mysdList, 'DS names')
    endif 
    if ( options%la ) then
      sdfid1 = mls_sfstart(filenames(i), DFACC_READ, hdfVersion=hdfVersion)
      print *, 'Attribute Names in: ', trim(filenames(i)), ' ', &
        & trim(options%root)//trim(options%DSName)
      if ( options%DSName /= ' ' ) then
        call GetAllHDF5AttrNames ( sdfid1, mysdList, &
          & DSName=trim(options%root)//trim(options%DSName) )
      else
        call GetAllHDF5AttrNames ( sdfid1, mysdList, &
          & groupName=trim(options%root) )
      endif
      call dump(mysdList, 'Attribute names')
      status = mls_sfend(sdfid1, hdfVersion=hdfVersion)
    endif
    if ( (options%attributes // options%datasets) == ' ' ) cycle
    if ( options%verbose ) then
      print *, 'Reading from: ', trim(filenames(i))
    endif
    sdfid1 = mls_sfstart(filenames(i), DFACC_READ, hdfVersion=hdfVersion)
    if (sdfid1 == -1 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      &  'Failed to open l2aux file ' // trim(filenames(i)) )
    end if
    if ( options%datasets /= ' ' ) &
      & call DumpHDF5DS ( sdfid1, trim(options%root), trim(options%datasets), &
      & rms=options%rms, stats=options%stats )
    if ( options%attributes /= ' ' ) then
      if ( options%DSName /= ' ' ) then
        call DumpHDF5Attributes ( sdfid1, trim(options%attributes), &
          & DSName=trim(options%root)//trim(options%DSName), &
          & rms=options%rms, stats=options%stats )
      else
        call DumpHDF5Attributes ( sdfid1, trim(options%attributes), &
          & groupName=trim(options%root), rms=options%rms, stats=options%stats )
      endif
    endif
	 status = mls_sfend(sdfid1, hdfVersion=hdfVersion)
    call sayTime('reading this file', tFile)
  enddo
  if ( .not. (options%la .or. options%ls) ) call sayTime('reading all files')
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
     print *, 'verbose?            ', options%verbose
     print *, 'list attributes  ?  ', options%la   
     print *, 'list datasets  ?    ', options%ls
     print *, 'stats  ?            ', options%stats  
     print *, 'rms    ?            ', options%rms    
     print *, 'root                ', options%root
     print *, 'attributes          ', options%attributes
     print *, 'datasets            ', options%datasets
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
      elseif ( filename(1:3) == '-A ' ) then
        options%attributes = '*'
        exit
      elseif ( filename(1:3) == '-D ' ) then
        options%datasets = '*'
        exit
      elseif ( filename(1:3) == '-nA ' ) then
        options%attributes = ''
        exit
      elseif ( filename(1:3) == '-nD ' ) then
        options%datasets = ''
        exit
      else if ( filename(1:3) == '-a ' ) then
        call getarg ( i+1+hp, options%attributes )
        i = i + 1
        exit
      else if ( filename(1:3) == '-d ' ) then
        call getarg ( i+1+hp, options%datasets )
        i = i + 1
        exit
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
        exit
      elseif ( filename(1:3) == '-la ' ) then
        options%la = .true.
        options%datasets = ''
        exit
      elseif ( filename(1:3) == '-ls ' ) then
        options%ls = .true.
        options%datasets = ''
        exit
      else if ( filename(1:3) == '-r ' ) then
        call getarg ( i+1+hp, options%root )
        i = i + 1
        exit
      else if ( filename(1:5) == '-rms ' ) then
        options%rms = .true.
        exit
      else if ( filename(1:3) == '-s ' ) then
        options%stats = .true.
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
      & 'Usage:l2auxdump [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options: -f filename     => add filename to list of filenames'
      write (*,*) '                           (can do the same w/o the -f)'
      write (*,*) '          -v              => switch on verbose mode'
      write (*,*) '          -A              => dump all attributes'
      write (*,*) '          -D              => dump all datasets (default)'
      write (*,*) '          -nA             => do not dump attributes (default)'
      write (*,*) '          -nD             => do not dump datasets'
      write (*,*) '          -r root         => limit to group based at root'
      write (*,*) '                             (default is "/")'
      write (*,*) '          -a a1,a2,..     => dump just attributes named a1,a2,..'
      write (*,*) '          -d d1,d2,..     => dump just datasets named a1,a2,..'
      write (*,*) '          -la             => just list attribute names in files'
      write (*,*) '          -ls             => just list sd names in files'
      write (*,*) '          -rms            => just print mean, rms'
      write (*,*) '          -s              => just show % statistics'
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

!==================
end program l2auxdump
!==================

! $Log$
