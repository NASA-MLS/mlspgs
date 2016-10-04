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
program l1bcat ! catenates l1b files, e.g. l1boa
!=================================

   use Dump_1, only: Dump
   use Hdf, only: DFACC_Create, DFACC_RDWR, DFACC_Read
   use HDF5, only: H5FIs_HDF5_F, H5GCreate_F, H5GClose_F
   use Intrinsic, only: l_hdf
   use io_stuff, only: get_lun
   use L1BData, only: cpL1BData
   use machine, only: hp, getarg
   use MLSFiles, only: HDFVersion_5, Dump, InitializeMLSFile, &
     & MLS_OpenFile, MLS_CloseFile
   use MLSCommon, only: MLSFile_T
   use MLSHDF5, only: GetAllHDF5DSNames, MLS_H5Open, MLS_H5Close
   use MLSStringLists, only: catLists, GetStringElement, &
     & Intersection, NumStringElements, &
     & RemoveElemFromList, StringElement, StringElementNum
   use output_m, only: output
   use PrintIt_m, only: Set_Config
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
! Catenate l1b dataset from list of input files to a single output file

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it
! LF95.Linux/test [options] [input files] -o [output file]

  integer, parameter ::          MAXFILES = 100
  type options_T
    logical               :: l2aux = .false. ! Don't complain if missing counterMAF?
    logical               :: verbose = .false.
    character(len=255) ::    outputFile= 'default.h5'        ! output filename
    logical ::               noDupDSNames = .false.            ! cp 1st, ignore rest
    logical     ::           list = .false.
    character(len=255) ::    DSNames = ' '              ! which datasets to copy
    character(len=255) ::    rename = ' '               ! how to rename them
    character(len=255), dimension(MAXFILES) :: filenames
  end type options_T

  type ( options_T ) :: options

  logical, parameter :: COUNTEMPTY = .false.
  integer, parameter :: MAXNUMSDPERFILE = 500
  integer, parameter :: MAXSDNAMESBUFSIZE = 80*MAXNUMSDPERFILE
  logical     :: createdYet
  logical, parameter :: DEEBUG = .true.
  character(len=255) :: filename          ! input filename
  integer            :: n_filenames
  integer     ::  i, j, error ! Counting indices & Error flags
  integer     :: elem
  logical     :: is_hdf5
  character (len=MAXSDNAMESBUFSIZE) :: mySdList
  integer :: numdsets
  integer :: numdsetssofar
  character(len=MAXSDNAMESBUFSIZE) :: rename = ' '               ! how to rename them
  character(len=MAXSDNAMESBUFSIZE) :: sdListAll
  character(len=255) :: sdName
  character (len=MAXSDNAMESBUFSIZE) :: tempSdList
  real        :: t1
  real        :: t2
  real        :: tFile
  !
  call set_config ( useToolkit = .false., logFileUnit = -1 )
  time_config%use_wall_clock = .true.
  CALL mls_h5open(error)
  n_filenames = 0
  createdYet = .false.
  do      ! Loop over filenames
     call get_filename( filename, n_filenames, options )
     if ( filename(1:1) == '-' ) cycle
     if ( filename == ' ' ) exit
     call h5fis_hdf5_f( trim(filename), is_hdf5, error )
     if ( .not. is_hdf5 ) then
       print *, 'Sorry--not recognized as hdf5 file: ', trim(filename)
       cycle
     endif
     n_filenames = n_filenames + 1
     options%filenames(n_filenames) = filename
  enddo
  call time_now ( t1 )
  if ( n_filenames == 0 ) then
    if ( options%verbose ) print *, 'Sorry no input files supplied'
  elseif ( options%list ) then
    do i=1, n_filenames
      print *, 'DS Names in: ', trim(options%filenames(i))
      call GetAllHDF5DSNames (trim(options%filenames(i)), '/', mysdList)
      call dump(mysdList, 'DS names')
    enddo
  else
    call catenate_all
  endif
  if ( .not. options%list ) call sayTime('copying all files')
  call mls_h5close(error)
contains

  ! Copy l1b data from one of the filenames to the outputFile
  ! Copy all the MAFs
  subroutine catenate_all
    if ( options%verbose ) print *, 'Copy l1b data to: ', trim(options%outputFile)
    numdsetssofar = 0
    sdListAll = ''
    do i=1, n_filenames
      call time_now ( tFile )
      if ( options%verbose ) then
        print *, 'Copying from: ', trim(options%filenames(i))
      endif
      if ( options%noDupDSNames.or. options%DSNames /= ' ' ) then
        call GetAllHDF5DSNames (trim(options%filenames(i)), '/', mysdList)
        numdsets = NumStringElements(trim(mysdList), COUNTEMPTY)
        if ( options%DSNames /= ' ' ) then
          tempSdList = mysdList
          mysdList = Intersection( options%DSNames, tempSdList )
          if ( mysdList == ' ' ) cycle
          rename = ' '
          do j=1, NumStringElements( mysdList, countEmpty )
            call GetStringElement( mysdList, sdName, j, countEmpty )
            elem = StringElementNum( options%DSNames, sdName, countEmpty )
            rename = catLists( rename, &
              & StringElement( options%rename, elem, countempty ) )
          enddo
        elseif ( numdsetssofar > 0 ) then
          ! Remove any duplicates
          do j=1, numdsetssofar
            call GetStringElement(sdListAll, sdName, j, countEmpty)
            tempSdList = mysdList
            call RemoveElemFromList (tempSdList, mysdList, trim(sdName))
          enddo
          if ( DEEBUG ) then
            print *, 'sds to cp'
            print *, trim(mysdList)
          endif
        endif
        call cpAllL1BData( trim(options%filenames(i)), &
        & trim(options%outputFile), &
        & hdfVersion=HDFVERSION_5, sdList=mysdList, rename=rename )
        tempSdList = sdListAll
        sdListAll = catlists(tempSdList, mysdList)
        numdsetssofar = NumStringElements(sdListAll, countEmpty)
      else
        call GetAllHDF5DSNames ( trim(options%filenames(i)), '/', mysdList )
        rename = mysdList
        call cpAllL1BData( trim(options%filenames(i)), &
        & trim(options%outputFile), &
        & hdfVersion=HDFVERSION_5, sdList=mysdList, rename=rename )
      endif
      call sayTime('copying this file', tFile)
      createdYet = .true.
    enddo
  end subroutine catenate_all

!------------------------- cpAllL1BData ---------------------
    ! copy all the l1bdata sets from filename1 to filename2
    subroutine cpAllL1BData( filename1, filename2, &
        & hdfVersion, sdList, rename )
      ! Args
      character(len=*), intent(in)            :: filename1
      character(len=*), intent(in)            :: filename2
      integer, intent(in)                     :: hdfversion
      character(len=*), intent(in)            :: sdList
      character(len=*), intent(in)            :: rename
      ! Internal variables
      integer                                 :: j
      type(MLSFile_t), pointer                :: L1BFile1
      type(MLSFile_t), pointer                :: L1BFile2
      type(MLSFile_t), target                 :: L1BFileT1
      type(MLSFile_t), target                 :: L1BFileT2
      character(len=128)                      :: sdNew
      character(len=128)                      :: sdOld
      integer                                 :: status
      ! Executable
      nullify ( L1BFile1, L1BFile2 )
      print *, 'sdList ' // trim(sdList)
      print *, 'rename ' // trim(rename)
      status = InitializeMLSFile ( L1BFileT1, type=l_hdf, access=DFACC_READ, &
        & name=filename1, HDFVersion=hdfVersion )
      status = InitializeMLSFile ( L1BFileT2, type=l_hdf, access=DFACC_CREATE, &
        & name=filename2, HDFVersion=hdfVersion )
      if ( createdyet ) L1BFileT2%access = DFACC_RDWR
      L1BFile1 => L1BFileT1
      L1BFile2 => L1BFileT2
      do j=1, NumStringElements( sdList, countEmpty=.true. )
        call GetStringElement ( sdList, sdOld, j, countEmpty=.true. )
        call GetStringElement ( rename, sdNew, j, countEmpty=.true. )
        print *, 'sdOld ', trim(sdOld)
        print *, 'sdNew ', trim(sdNew)
        call CpL1BData ( L1BFile1, L1BFile2, sdOld, sdNew, l2aux=options%l2aux )
        L1BFileT2%access = DFACC_RDWR
      enddo
    end subroutine cpAllL1BData

!------------------------- get_filename ---------------------
    subroutine get_filename(filename, n_filenames, options)
    ! Added for command-line processing
     character(LEN=255), intent(out)   :: filename          ! filename
     integer, intent(in)               :: n_filenames
     type ( options_T ), intent(inout) :: options
     ! Local variables
     type(MLSFile_t)                   :: L1BFile
     integer                           :: error = 1
     integer                           :: fileaccess
     integer, save                     :: i = 1
     character(LEN=255)                :: groupname
     integer                           :: grp_id
     integer                           :: status
  ! Get inputfile name, process command-line args
  ! (which always start with -)
    do
      call getarg ( i+hp, filename )
      ! print *, i, ' th Arg: ', trim(filename)
      error = 0
      if ( filename(1:1) /= '-' ) exit
      if ( filename(1:3) == '-h ' ) then
        call print_help
      elseif ( filename(1:3) == '-o ' ) then
        call getarg ( i+1+hp, options%outputFile )
        i = i + 1
        exit
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
        exit
      elseif ( filename(1:3) == '-l ' ) then
        options%list = .true.
        exit
      elseif ( filename(1:6) == '-l2aux' ) then
        options%l2aux = .true.
        exit
      elseif ( filename(1:3) == '-no' ) then
        options%noDupDSNames = .true.
        exit
      else if ( filename(1:3) == '-s ' ) then
        call getarg ( i+1+hp, options%DSNames )
        i = i + 1
        exit
      else if ( filename(1:3) == '-r ' ) then
        call getarg ( i+1+hp, options%rename )
        i = i + 1
        exit
      else if ( filename(1:3) == '-f ' ) then
        call getarg ( i+1+hp, filename )
        i = i + 1
        exit
      else if ( filename(1:3) == '-Sf' ) then
        call getarg ( i+1+hp, filename )
        call read_file ( filename, options%DSNames )
        i = i + 1
        exit
      else if ( filename(1:3) == '-Rf' ) then
        call getarg ( i+1+hp, filename )
        call read_file ( filename, options%rename )
        i = i + 1
        exit
      else if ( filename(1:7) == '-create' ) then
        call getarg ( i+1+hp, groupname )
        fileaccess = DFACC_Create
        if ( createdYet ) fileaccess = DFACC_RDWR
        print *, 'Attempting to create group ' // trim(groupname)
        status = InitializeMLSFile ( L1BFile, type=l_hdf, &
          & access=fileaccess, content='l2aux', name=options%outputFile, &
          & HDFVersion=HDFVERSION_5 )
        call mls_openFile ( L1BFile, status )
        if ( status /= 0 ) then
          print *, 'could not create file ' // trim(options%outputFile)
          stop
        endif
        call h5GCreate_f ( L1BFile%fileID%f_id, groupname, grp_id, status )
        if ( status /= 0 ) then
          print *, 'could not create group ' // trim(groupname)
          stop
        endif
        call h5GClose_f ( grp_id, status )
        call mls_CloseFile( L1BFile )
        createdYet = .true.
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
      print *,  "Enter the name of the HDF5 l1b file. " // &
       &  "The default output file name will be used."
      read(*,'(a)') filename
    endif

  end subroutine get_filename

!------------------------- read_file ---------------------
  subroutine read_file ( filename, string )
    ! Read contents of file into comma-separated string
    ! Args
    character(len=*), intent(in)        :: filename
    character(len=*), intent(out)       :: string
    ! Internal varaibles
    character(len=80)                   :: line
    integer                             :: lun
    integer                             :: status
    ! Executable
    string = ' '
    call get_lun( lun )
    open(UNIT=lun, form='formatted', &
      & file=trim(FileName), status='old', iostat=status )
    if ( status /= 0 ) then
      write(*,*) 'read_file error Unable to open textfile ' // &
        & trim(FileName)
      return
    endif
    do
      read( lun, *, end=500, err=50 ) line
      line = adjustl(line)
      if (len_trim(line) < 1 .or. line(1:1) == '#' ) cycle
      string = trim(string) // ',' // line
500   status = -1
50    if ( status /= 0 ) exit
    enddo
    ! Snip off leading ','
    if ( string(1:1) == ',' ) then
      string(1:1) = ' '
      string = adjustl(string)
    endif
    close(unit=lun)
    
  end subroutine read_file

!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'Usage:l1bcat [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options:'
      write (*,*) ' -f filename            => add filename to list of filenames'
      write (*,*) '                        (can do the same w/o the -f)'
      write (*,*) ' -o ofile               => copy data sets to ofile'
      write (*,*) ' -create gname          => create group named gname'
      write (*,*) '                        (may be repeated)'
      write (*,*) '                        (must appear after -o ofile)'
      write (*,*) ' -v                     => switch on verbose mode'
      write (*,*) ' -l                     => just list lib names in files'
      write (*,*) ' -l2aux                 => input files may be missing counterMAF'
      write (*,*) ' -nodup                 => if dup dataset names, cp 1st only'
      write (*,*) ' -h                     => print brief help'
      write (*,*) ' -s name1,name2,..      => copy only datasets so named; otherwise all'
      write (*,*) ' -r rename1,rename2,..  => if and how to rename the copied datasets'
      write (*,*) ' -Sf file               => copy only datasets named in file'
      write (*,*) ' -Rf file               => rename them according to names in file'
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
end program l1bcat
!==================

! $Log$
! Revision 1.2  2016/06/10 16:15:31  pwagner
! Corrected call to cpAllL1BData; clarified meaning of -l2aux option
!
! Revision 1.1  2016/04/20 00:21:13  pwagner
! First commit
!
