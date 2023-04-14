! Copyright 2023, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=================================
program createattributes ! creates file attributes
!=================================

   use Allocate_Deallocate, only: Allocate_Test, DeAllocate_Test
   use Dump_0, only: Dump
   use HDF, only: Dfacc_Rdwr
   use HDF5, only: H5DClose_f, H5DOpen_f, H5GClose_f, H5GOpen_f
   use HDFEOS5, only: MLS_CharType
   use Intrinsic, only: L_HDF, L_Swath
   use Io_Stuff, only: Get_NLines, Read_TextFile
   use Machine, only: Hp, Getarg
   use MLSCommon, only: FileNameLen, MLSFile_T
   use MLSFiles, only: HDFVersion_5, &
     & Dump, MLS_CloseFile, MLS_OpenFile, MLS_Exists, &
     & InitializeMLSFile
   use MLSHDFEOS, only: MLS_EhwrGlAtt
   use MLSHDF5, only: MLS_H5open, MLS_H5close
   use MLSHDF5, only: MakeHDF5Attribute
   use MLSMessageModule, only: MLSMessageConfig, &
     & MLSMessage, MLSMSG_Error, MLSMSG_Warning
   use MLSStringLists, only: GetStringElement
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
! Creates or resets values of file attributes in 
! (1) an L2GPData file
! (2) an L2AUXData file
! according to a text file containing lines of name=value

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"

! Then run it something like this
! Usage:
! (1)
!    createattributes \ 
!       -l2gp \
!       -Vf attr_file \
!       L2GP_File
! (2)
!    createattributes \ 
!       -l2aux \
!       -Vf attr_file \
!       L2AUX_File

! Here is an example of an attr_file
! #------env file for createattributes----
! newattr_a=newvalue_a
! newattr_b=newvalue_b
! newattr_c=newvalue_c
! newattr_d=newvalue_d
! GranuleYear=3023
! root=/
! DSName=chunk number
! attr_a=chunk_a
! attr_b=chunk_b
! attr_c=chunk_c
! attr_d=chunk_d
! 
! You can observe that the rhs, even in the case of chunk number with
! its embedded space, is not surrounded by quote marks. In fact,
! surrounding the rhs with quote marks will generate an hdf error, so
! don't do that.

! Notes and Limitations
! (1) The attribute values will be character-valued scalars.
!     If a prior attribute with the same name existed, it will
!     be silently deleted, then replaced with the new one.
!     If the prior attribute was, say integer-valued, the new one
!     will be character-valued, which may confuse whatever application
!     is intended to read it.
! (2) One possible improvement would be to allow for attribute values
!     that are numeric. A bit of syntactic trickery to do this
!     would be to allow an optional keyword preceding the attribute name.
!     E.g.,
!      integer FirstMAF=5089784
!      double TAI93At0zOfGranule=9.43142e+08
! (3) More ambitious plans would allow for numeric arrays using syntax like
!      integer(:) OrbitNumber=97608, 97609, 97610, 97611, 97612, 97613, 97614, 97615, 97616, 97617, 97618, 97619, 97620, 97621, 97622, -1
!      double(:) OrbitPeriod=5933.18, 5933.15, 5933.15, 5933.21, 5933.29, 5933.29, 5933.26, 5933.18, 5933.14, 5933.14, 5933.16, 5933.24, 5933.31, 5933.25, 5933.25, 0

  type options_T
    logical               :: debug       = .false.
    logical               :: silent      = .false.
    logical               :: verbose     = .false.
    logical               :: dryrun      = .false.
    character(len=255)    :: attr_file   = ' ' 
    character(len=8  )    :: filetype    = ' '       ! 'l2gp' or 'l2aux'
    character(len=255)    :: root        = '/' 
    character(len=255)    :: DSName      = ' ' 
  end type options_T
  
  type ( options_T ) :: options

  integer, parameter ::          MAXFILES = 10 ! Usually just 1
  logical, parameter ::          countEmpty = .true.
  character(len=255) :: filename          ! input filename
  character(len=255), dimension(MAXFILES) :: filenames
  integer            :: n_filenames
  integer     ::  status, error ! Counting indices & Error flags
  real        :: t1
  real        :: t2
  real        :: tFile
  type( MLSFile_T ) :: L2GPFile
  integer :: ACCESS
  integer :: k
  integer :: nLines
  character(len=FileNameLen), dimension(:), pointer  :: attrLines => null()
  ! 
  MLSMessageConfig%useToolkit = .false.
  MLSMessageConfig%logFileUnit = -1
  time_config%use_wall_clock = .true.
  CALL mls_h5open(error)
  ACCESS = Dfacc_Rdwr
  n_filenames = 0
  do      ! Loop over filenames
     call get_filename( filename, n_filenames, options )
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
    print *, 'Sorry no files to insert'
  elseif ( len_trim(options%filetype) == 0 ) then
    print *, 'Sorry, type must be either l2gp or l2aux'
  else
    call time_now ( t1 )
    if ( options%silent ) then
      options%debug   = .false.
      options%verbose = .false.
    endif

    if ( options%verbose ) &
      & call Dump ( filenames(1:n_filenames), 'files', width=1 )
    if ( len_trim(options%attr_file) < 1 ) then
      print *, 'Sorry, you must supply a file of attribute names=values'
      call MLSMessage( MLSMSG_Error, ModuleName, &
        & 'Sorry, you must supply a file of attribute names=values' )
    endif
    call get_nLines ( options%attr_file, nLines )
    if ( nLines < 0 ) then
      print *, 'Sorry, unable to open attribute file ' // trim(options%attr_file)
      call MLSMessage( MLSMSG_Error, ModuleName, &
        & 'Sorry, unable to open attribute file ' // trim(options%attr_file) )
    elseif ( nLines < 1 ) then
      print *, '0 lines in attribute file ' // trim(options%attr_file)
      call MLSMessage( MLSMSG_Warning, ModuleName, &
        & 'Unexpectedly empty attribute file ' // trim(options%attr_file) )
    endif
    if ( options%verbose ) &
      & print *, "attribute file: ", trim(options%attr_file)
    call allocate_test ( attrLines, nLines, 'attrLines', &
      & trim(ModuleName) // 'processLine' )
    attrLines = ' '
    call read_textFile( options%attr_file, attrLines )
    if ( options%verbose ) then
      print *, 'attributes read from ' // trim(options%attr_file)
      print *, 'nLines', nLines 
    endif
    if ( options%filetype == 'l2gp' ) then
      status = InitializeMLSFile( L2GPFile, type=l_swath, access=ACCESS, &
        & content='L2GP', name=trim(filenames(1)), hdfVersion=HDFVERSION_5 )
    else
      status = InitializeMLSFile( L2GPFile, type=l_hdf, access=ACCESS, &
        & content='L2GP', name=trim(filenames(1)), hdfVersion=HDFVERSION_5 )
    endif
    call mls_openFile( L2GPFile, Status )
    do k=1, nLines
      attrLines(k) = adjustl(attrLines(k))
      if ( options%verbose ) print *,  trim(attrLines(k))
      ! skip comment lines
      if ( index(attrLines(k), '#') == 1 ) cycle
      if ( options%filetype == 'l2gp' ) then
        call create_l2gpattr( L2GPFile, attrLines(k) )
      else
        call create_l2auxattr( L2GPFile, attrLines(k) )
      endif
    enddo
    call deallocate_test ( attrLines, 'attrLines', &
      & trim(ModuleName) // 'processLine' )
    call mls_CloseFile( L2GPFile, Status )
  endif
  call mls_h5close(error)
contains

!------------------------- create_l2gpattr ---------------------
! Create the file-level attribute and its value

  subroutine create_l2gpattr( L2GPFile, Line )
    ! Dummy args
    type( MLSFile_T ) :: L2GPFile
    character(len=FileNameLen) :: Line
    ! Internal variables
    character(len=fileNameLen)      :: name
    character(len=fileNameLen)      :: valu
    character(len=1), parameter     :: eqls = '='
    integer                         :: status
    
    ! --------------------------------------------------------------------
    ! --------------------------------------------------------------------
    ! Executable
    ! ------------------------------------------------------------------------
    call time_now ( tFile )
    call getStringElement( line, name, 1, countEmpty, eqls )
    call getStringElement( line, valu, 2, countEmpty, eqls )
    if ( options%dryrun .or. options%debug ) then
      print *, 'name: ', trim(name)
      print *, 'value: ', trim(valu)
    endif
    if ( .not. options%dryrun ) then
      status = mls_EHwrglatt( L2GPFile%fileID%f_id, &
        & trim(name), MLS_CHARTYPE, 1, &
        & trim(valu) )
      call sayTime('Creating this attribute')
    endif
  end subroutine create_l2gpattr

!------------------------- create_l2auxattr ---------------------
! Create the file-level attribute and its value

  subroutine create_l2auxattr( L2GPFile, Line )
    ! Dummy args
    type( MLSFile_T ) :: L2GPFile
    character(len=FileNameLen)  :: Line
    ! Internal variables
    character(len=fileNameLen)      :: name
    character(len=fileNameLen)      :: valu
    character(len=1), parameter     :: eqls = '='
    integer                         :: grp_id
    integer                         :: setid
    integer                         :: status
    
    ! --------------------------------------------------------------------
    ! --------------------------------------------------------------------
    ! Executable
    ! ------------------------------------------------------------------------
    call time_now ( tFile )
    call getStringElement( line, name, 1, countEmpty, eqls )
    call getStringElement( line, valu, 2, countEmpty, eqls )
    if ( options%dryrun .or. options%debug ) then
      print *, 'name: ', trim(name)
      print *, 'value: ', trim(valu)
    endif
    ! Are we using this line to set the root or the DSName?
    if ( trim(name) == 'root' ) then
      options%root = valu
      return
    elseif ( trim(name) == 'DSName' ) then
      options%DSName = valu
      return
    elseif ( .not. options%dryrun ) then
      call h5gopen_f( L2GPFile%fileID%f_id, trim(options%root), grp_id, status )
      if ( status /= 0 ) then
        call MLSMessage( MLSMSG_Error, ModuleName, &
          & 'Sorry, Cant open group: ' // trim(options%root) )
      endif
      ! Is our attribute a group-level attribute or specific to a dataset?
      if ( len_trim(options%DSName) < 1 ) then

        call MakeHDF5Attribute( grp_id, &
          & trim(name), trim(valu), skip_if_already_there=.false. )
      else
        call h5dOpen_f ( grp_id, trim(options%DSName), setID, status )
        call MakeHDF5Attribute( setid, &
          & trim(name), trim(valu), skip_if_already_there=.false. )
        call h5dclose_f ( setid, status )
      endif
      call h5gclose_f( grp_id, status )
      call sayTime('Creating this attribute')
    endif
  end subroutine create_l2auxattr

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
      elseif ( filename(1:4) == '-deb' ) then
        options%debug = .true.
        exit
      elseif ( filename(1:4) == '-dry' ) then
        options%dryrun = .true.
        exit
      elseif ( filename(1:7) == '-silent' ) then
        options%silent = .true.
        exit
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
        exit
      else if ( filename(1:3) == '-f ' ) then
        call getarg ( i+1+hp, filename )
        i = i + 1
        exit
      else if ( filename(1:5) == '-l2gp' ) then
        options%filetype = 'l2gp'
        exit
      else if ( filename(1:6) == '-l2aux' ) then
        options%filetype = 'l2aux'
        exit
      else if ( filename(1:3) == '-Vf' ) then
        call getarg ( i+1+hp, options%attr_file )
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
      print *,  "Enter the name of the L2GP file. " // &
       &  "The default output file name will be used."
      read(*,'(a)') filename
    endif
    
  end subroutine get_filename
!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'Creates or resets file attributes'
      write (*,*) &
      & '    of an L2GP or an L2AUX file'
      write (*,*) &
      & 'Usage: createattributes [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options:'
      write (*,*) ' -dryrun       => dont execute, just describe'
      write (*,*) ' -f filename   => add filename to list of l2gp filenames'
      write (*,*) '                  (can do the same w/o the -f)'
      write (*,*) ' -Vf attr_file => settings file; one line per attribute'
      write (*,*) ' -debug        => switch on debug mode'
      write (*,*) ' -v            => switch on verbose mode'
      write (*,*) ' -silent       => switch on silent mode'
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
    if ( options%silent ) return
    call output ( "Timing for " // what // " = " )
    call output ( dble(t2 - myt1), advance = 'yes' )
  end subroutine SayTime
!------------------------- dgetarg ---------------------
  subroutine dgetarg ( pos, darg )
   integer, intent(in) :: pos
   double precision, intent(out) :: darg
   character(len=16) :: arg
   call getarg ( pos, arg )
   read(arg, *) darg
  end subroutine dgetarg
!------------------------- igetarg ---------------------
  subroutine igetarg ( pos, iarg )
   integer, intent(in) :: pos
   integer, intent(out) :: iarg
   character(len=16) :: arg
   call getarg ( pos, arg )
   read(arg, *) iarg
  end subroutine igetarg

!==================
end program createattributes
!==================

! $Log$
