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
program L2AUXcat ! catenates split L2AUX files, e.g. dgm
!=================================

   use Dump_0, only: Dump
   use Dump_1, only: Dump
   use HDF, only: DFACC_Create, DFACC_RDWR, DFACC_RDOnly
   use HDF5, only: H5FIs_HDF5_F, H5GCreate_F, H5GClose_F
   use Intrinsic, only: L_HDF
   use Io_Stuff, only: Get_Lun, Read_TextFile
   use L2AUXData, only: L2AUXData_T, MaxSDNamesBufSize, &
     & CpL2AUXData, DestroyL2AUXContents, ReadL2AUXData, &
     & ResizeL2AUXData, WriteL2AUXData
   use Machine, only: Hp, Getarg
   use MLSFiles, only: HDFVersion_5, Dump, InitializeMLSFile, &
     & MLS_OpenFile, MLS_CloseFile, Split_Path_Name
   use MLSFinds, only: FindFirst, FindLast
   use MLSCommon, only: MLSFile_T, DefaultUndefinedValue
   use MLSHDF5, only: GetAllHDF5DSNames, MLS_H5open, MLS_H5close
   use MLSKinds, only: R8
   use MLSStrings, only: Asciify
   use MLSStringLists, only: CatLists, GetStringElement, &
     & Intersection, NumStringElements, &
     & RemoveElemFromList, StringElement, StringElementNum
   use Output_M, only: Output
   use PCFHdr, only: H5_ReadGlobalAttr, H5_WriteGlobalAttr
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

! Brief description of program
! Catenate L2AUX from list of input files to a single output file
! May thus "unsplit" the dgm files

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it
! LF95.Linux/test [options] [input files] -o [output file]

  integer, parameter ::   MAXFILES = 750
  type Options_T
    logical            :: verbose = .false.
    character(len=255) :: glAttrFile= ''                ! file with global attrs
    character(len=255) :: outputFile= 'default.h5'      ! output filename
    character(len=255) :: inputFile= ''                 ! file with list of inputs       
    logical            :: noDupDSNames = .false.        ! cp 1st, ignore rest
    logical            :: ignoreFills  = .false.
    logical            :: list = .false.
    logical            :: flatten = .false.
    character(len=255) :: DSNames = ' '                 ! which datasets to copy
    character(len=255) :: rename = ' '                  ! how to rename them
    character(len=255), dimension(MAXFILES) :: filenames
  end type Options_T

  type ( Options_T ) :: options

  logical, parameter :: COUNTEMPTY = .false.
  logical     :: createdYet
  logical, parameter :: DEEBUG = .false.
  character(len=255) :: filename          ! input filename
  character(LEN=255)                :: groupname
  integer            :: n_filenames
  integer     ::  i, j, status, error ! Counting indices & Error flags
  integer     :: elem
  logical     :: is_hdf5
  type(MLSFile_t)                   :: L2AUXFile1
  type(MLSFile_t)                   :: L2AUXFile2
  character (len=MaxSDNamesBufSize) :: mySdList
  integer :: numdsets
  integer :: numdsetssofar
  character(len=255) ::    rename = ' '               ! how to rename them
  character(len=MaxSDNamesBufSize) :: sdListAll
  character(len=255) :: sdName
  character (len=MaxSDNamesBufSize) :: tempSdList
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
  ! Did we use the -F commandline option to supply the input file list?
  if ( len_trim(options%inputFile) > 0 ) then
    print *, 'inputFile :', trim(options%inputFile)
    ! call read_file ( options%inputFile, allInputFiles )
    ! print *, 'File list :', trim(allInputFiles)
    ! n_filenames = NumStringElements ( allInputFiles, countEmpty )
    call read_textfile ( options%inputFile, options%filenames, nLines=n_filenames )
    ! Must replace all nulls with spaces
    do i=1, n_filenames
      options%filenames(i) = asciify( options%filenames(i), how='snip' )
    enddo
    print *, 'n :', n_filenames
    ! call List2Array( allInputFiles, options%filenames, countEmpty )
    if ( options%verbose ) then
      call dump( options%filenames(1:n_filenames), width=1, options='-t' )
    endif
  endif
  call time_now ( t1 )
  ! stop
  if ( n_filenames == 0 ) then
    if ( options%verbose ) print *, 'Sorry no input files supplied'
  elseif ( options%list ) then
    do i=1, n_filenames
      print *, 'DS Names in: ', trim(options%filenames(i))
      call GetAllHDF5DSNames (trim(options%filenames(i)), '/', mysdList)
      call dump(mysdList, 'DS names')
    enddo
  elseif ( options%ignoreFills ) then
    call catenate_non_Fills
  else
    call catenate_all
  endif
  ! Do we copy global attributes from a particular file
  if ( len_trim(options%glAttrFile) > 0 ) then
    status = InitializeMLSFile( l2auxFile1, type=l_hdf, access=DFACC_RDONLY, &
      & content='l2aux', name=options%glAttrFile, hdfVersion=HDFVERSION_5 )
    status = InitializeMLSFile( l2auxFile2, type=l_hdf, access=DFACC_RDWR, &
      & content='l2aux', name=options%outputFile, hdfVersion=HDFVERSION_5 )
    call h5_ReadGlobalAttr ( l2auxFile1 )
    call h5_WriteGlobalAttr ( l2auxFile2 )
  endif
  if ( .not. options%list ) call sayTime('copying all files')
  call mls_h5close(error)
contains

!------------------------- catenate_non_Fills ---------------------
! catenate the MAFs not containing Fills
! In terms of the MAF numbers, we will 
! write starting with N+1 and continue up to M
! So, if a represents values already in the file, and 
! b represents value we will overwrite, then after overwriting we will have

! a[1]   ..  a[N] b[k] b[k+1] .. b[M] a[M+1] ..
! <---   N   --->  <---   M-k+1  --->

  subroutine catenate_non_Fills
    ! Internal variables
    integer :: i1
    integer :: i2
    integer :: iBloc
    integer :: jj
    type(l2auxData_T) :: l2aux
    type( MLSFile_T ) :: l2auxFile
    type(MLSFile_T), dimension(MAXFILES)    :: L2AUXFiles
    type(l2auxData_T) :: ol2aux
    integer :: l2FileHandle
    integer :: N ! last profile number of a
    integer :: M ! size of b
    real(r8) :: MissingValue
    integer :: nBlocs
    integer, parameter :: numFilesPerBloc = 40
    integer :: status
    integer, parameter :: FORBDIMS = 10
    character(len=42), dimension(FORBDIMS), parameter :: forbiddens = (/ &
      & 'HDFEOS INFORMATION/coremetadata.0        ', &
      & 'HDFEOS INFORMATION/xmlmetadata           ', &
      & 'PCF                                      ', &
      & 'chunk number                             ', &
      & 'leap seconds                             ', &
      & 'master.ident                             ', &
      & 'phase timing                             ', & 
      & 'section timing                           ', & 
      & 'solar zenith                             ', & 
      & 'utc pole                                 ' /)
      
    ! Executable
    MissingValue = real(defaultUndefinedValue, r8)
    if ( options%verbose ) print *, 'Catenate l2aux data to: ', &
      & trim(options%outputFile)
    if ( DEEBUG ) then
      print *, 'files to concatenate'
      do i=1, n_filenames
        print *, trim(options%filenames(i))
      enddo
    endif
    call GetAllHDF5DSNames ( trim(options%filenames(1)), '/', mysdList )
    if ( DEEBUG ) then
      print *, 'datasets in file'
      print *, trim(mysdList)
    endif
    if ( options%DSNames /= ' ' ) then
      tempSdList = mysdList
      mysdList = Intersection( options%DSNames, tempSdList )
      if ( mysdList == ' ' ) return
    endif
    status = InitializeMLSFile( l2auxFile, type=l_hdf, access=DFACC_CREATE, &
      & content='l2aux', name='unknown', hdfVersion=HDFVERSION_5 )
    l2auxFile%name = options%outputFile
    l2auxFile%stillOpen = .false.
    call mls_openFile( l2auxFile, Status )
    if ( status /= 0 ) then
      print *, 'Unable to open trim(options%outputFile'
      stop
    endif
    l2FileHandle = l2auxFile%FileID%f_id
    
    nBlocs = (n_filenames - 1)/numFilesPerBloc + 1

    do jj=1, NumStringElements( mysdList, countEmpty ) ! Loop of datasets
      call GetStringElement( mysdList, sdName, jj, countEmpty )
      sdName = adjustl(sdName)
      call time_now ( tFile )
      ! Must skip forbidden datasets
      if ( any(trim(sdName) == forbiddens) ) cycle
      if ( trim(sdName) == 'HDFEOS INFORMATION/coremetadata.0' ) then
        print *, trim(sdName)
        print *, trim(forbiddens(1))
        print *, trim(sdName) == trim(forbiddens(1))
        stop
      endif
      if ( options%verbose ) print *, 'Catenating dataset: ', trim(sdName)
      N = 0
      M = 0
      i2 = 0
      do iBloc = 1, nBlocs ! Loop of blocs
        i1 = i2 + 1
        i2 = min( i2 + numFilesPerBloc, n_filenames )
        ! 1st open the files
        do i=i1, i2
          status = InitializeMLSFile( l2auxFiles(i), type=l_hdf, access=DFACC_RDONLY, &
            & content='l2aux', name=options%filenames(i), hdfVersion=HDFVERSION_5 )
          l2auxFiles(i)%name = options%filenames(i)
          l2auxFiles(i)%stillOpen = .false.
          call mls_openFile( l2auxFiles(i), Status )
          if ( status /= 0 ) then
            print *, 'Unable to open trim(options%filenames(i)'
            stop
          endif
          if ( options%verbose .and. jj == 1 ) call Dump( l2auxFiles(i), details=1 )
        enddo
        ! Read each file, putting its values into ol2aux
        do i=i1, i2
          if ( options%verbose ) print *, 'Reading from: ', trim(options%filenames(i))
          if ( i == 1 ) then
            call Readl2auxData( l2auxFiles(i)%FileID%f_id, trim(sdName), ol2aux, &
              & hdfVersion=HDFVERSION_5 )
            print *, 'shape l2auxvalues (after reading first): ', shape(ol2aux%values)
            cycle
          endif
          call Readl2auxData( l2auxFiles(i)%FileID%f_id, trim(sdName), l2aux, &
            & hdfVersion=HDFVERSION_5)
          ! We will use ChunkNumbers to determine N and M
          N = FindFirst( l2aux%values(1,1,:) /= MissingValue ) - 1
          M = FindLast( l2aux%values(1,1,:) /= MissingValue )
          call ResizeL2AUXData ( ol2aux, M )
          ol2aux%values      (:,:,N+1:M) = l2aux%values        (:,:,N+1:M)    

          call Destroyl2auxContents( l2aux )
        enddo
        ! Last, close each file
        do i=i1, i2
          call mls_closeFile( l2auxFiles(i), Status )
        enddo
      enddo ! End Loop of blocs
      call sayTime('Reading this dataset', tFile)
      ! print *, 'shape l2auxvalues (before writing): ', shape(ol2aux%values)
      call Writel2auxData ( ol2aux, l2FileHandle, status, sdName, &
        & hdfVersion=HDFVERSION_5 )
      call sayTime('Writing this dataset', tFile)
      call Destroyl2auxContents( ol2aux )
    enddo ! End Loop of datasets
    call mls_closeFile( l2auxFile, Status )
    call sayTime( 'catenating all datasets' )
  end subroutine catenate_non_Fills

  ! Copy l2aux data from one of the filenames to the outputFile
  ! Copy all the MAFs
  subroutine catenate_all
    character(len=128) :: grpname
    character(len=128) :: itemname
    character(len=16) :: myOptions
    myOptions = ' '
    if ( options%verbose ) myOptions = '-v'
    if ( options%verbose ) print *, 'Copy l2aux data to: ', trim(options%outputFile)
    numdsetssofar = 0
    sdListAll = ''
    do i=1, n_filenames
      call time_now ( tFile )
      if ( options%verbose ) then
        print *, 'Copying from: ', trim(options%filenames(i))
      endif
      call GetAllHDF5DSNames (trim(options%filenames(i)), '/', mysdList)
      numdsets = NumStringElements( trim(mysdList), COUNTEMPTY )
      ! Remove metadata, etc.
      call RemoveElemFromList ( mysdList, tempSdList, 'HDFEOS INFORMATION/coremetadata.0' )
      call RemoveElemFromList ( tempSdList, mysdList, 'HDFEOS INFORMATION/xmlmetadata' )
      call RemoveElemFromList ( mysdList, tempSdList, 'PCF' )
      call RemoveElemFromList ( tempSdList, mysdList, 'l2cf' )
      call RemoveElemFromList ( mysdList, tempSdList, 'leap seconds' )
      call RemoveElemFromList ( tempSdList, mysdList, 'utc pole' )
      if ( options%noDupDSNames .or. options%DSNames /= ' ' ) then
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
        call cpL2AUXData( trim(options%filenames(i)), &
        & trim(options%outputFile), create2=.not. createdYet, &
        & hdfVersion=HDFVERSION_5, sdList=mysdList, rename=rename, options=myOptions )
        tempSdList = sdListAll
        sdListAll = catlists(tempSdList, mysdList)
        numdsetssofar = NumStringElements(sdListAll, countEmpty)
      elseif ( options%flatten ) then
        sdListAll = ' '
        do j=1, numdsets
          call GetStringElement( mysdList, sdName, j, countEmpty )
          call split_path_name ( sdName, grpname, itemname )
          ! Did we create a group to put all the datsets under?
          if ( len_trim(groupname) > 0 ) &
            & itemname = trim(groupname) // '/' // itemname
          print *, 'sdname: ', trim(sdName)
          print *, 'grpname: ', trim(grpname)
          print *, 'itemname: ', trim(itemname)
          sdListAll = catLists ( sdListAll,  itemname )
        enddo
        print *, 'sdListAll: ', trim(sdListAll)
        ! stop
        call cpL2AUXData( trim(options%filenames(i)), &
        & trim(options%outputFile), create2=.not. createdYet, &
        & hdfVersion=HDFVERSION_5, sdList=mysdList, rename=sdListAll, options=myOptions )
        tempSdList = sdListAll
        sdListAll = catlists(tempSdList, mysdList)
        numdsetssofar = NumStringElements(sdListAll, countEmpty)
      else
        call cpL2AUXData( trim(options%filenames(i)), &
        & trim(options%outputFile), create2=.not. createdYet, &
        & hdfVersion=HDFVERSION_5, options=myOptions )
      endif
      call sayTime('copying this file', tFile)
      createdYet = .true.
    enddo
  end subroutine catenate_all
!------------------------- get_filename ---------------------
    subroutine get_filename(filename, n_filenames, options)
    ! Added for command-line processing
     character(LEN=255), intent(out)   :: filename          ! filename
     integer, intent(in)               :: n_filenames
     type ( Options_T ), intent(inout) :: options
     ! Local variables
     type(MLSFile_t)                   :: L2AUXFile
     integer                           :: error = 1
     integer                           :: fileaccess
     integer, save                     :: i = 1
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
      elseif ( filename(1:5) == '-flat' ) then
        options%flatten = .true.
        exit
      elseif ( filename(1:3) == '-l ' ) then
        options%list = .true.
        exit
      elseif ( filename(1:4) == '-ign' ) then
        options%ignoreFills = .true.
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
      else if ( filename(1:3) == '-F ' ) then
        call getarg ( i+1+hp, options%inputFile )
        i = i + 1
      else if ( filename(1:3) == '-f ' ) then
        call getarg ( i+1+hp, filename )
        i = i + 1
        exit
      else if ( filename(1:3) == '-g ' ) then
        call getarg ( i+1+hp, options%glAttrFile )
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
        status = InitializeMLSFile ( L2AUXFile, type=l_hdf, &
          & access=fileaccess, content='l2aux', name=options%outputFile, &
          & HDFVersion=HDFVERSION_5 )
        call mls_openFile ( L2AUXFile, status )
        call h5GCreate_f ( L2AUXFile%fileID%f_id, groupname, grp_id, status )
        call h5GClose_f ( grp_id, status )
        call mls_CloseFile( L2AUXFile )
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
    & 'Usage:l2auxcat [options] [filenames]'
    write (*,*) ' Options: '
    write (*,*) '-f filename   => add filename to list of filenames'
    write (*,*) '               (can do the same w/o the -f)'
    write (*,*) '-F infile     => read list of filenames from infile'
    write (*,*) '-g glattrfile => read global attrs from glattrfile'
    write (*,*) '-create gname => create group named gname'
    write (*,*) '              (may be repeated)'
    write (*,*) '-o ofile      => copy data sets to ofile'
    write (*,*) '-v            => switch on verbose mode'
    write (*,*) '-flatten      => copy the datasets to the root / of ofile'
    write (*,*) '                 or to gname if -create gname among options'
    write (*,*) '-l            => just list l2aux names in files'
    write (*,*) '-ign          => ignore MAFs with Fill Values'
    write (*,*) '-nodup        => if dup dataset names, cp 1st only'
    write (*,*) '-h            => print brief help'
    write (*,*) '-s name1,name2,..'
    write (*,*) '              => copy only datasets so named; otherwise all'
    write (*,*) '-r rename1,rename2,..'
    write (*,*) '              => if and how to rename the copied datasets'
    write (*,*) '-Sf file      => copy data only sets named in file'
    write (*,*) '-Rf file      => rename them according to names in file'
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
end program L2AUXcat
!==================

! $Log$
! Revision 1.10  2017/08/30 23:03:21  pwagner
! Added the -flatten commandline option
!
! Revision 1.9  2017/05/17 22:21:25  pwagner
! Works properly with dual-l2 nrt scripts
!
! Revision 1.8  2016/08/25 22:56:29  pwagner
! Can now successfully cat all 351 chunks of Pleiades-style dgm
!
! Revision 1.7  2016/06/13 23:26:41  pwagner
! Added commandline option -F; upped max num of input files tto 750
!
! Revision 1.6  2015/08/05 20:36:28  pwagner
! Option -ign can catenate non-Fills
!
! Revision 1.5  2013/08/23 02:51:47  vsnyder
! Move PrintItOut to PrintIt_m
!
! Revision 1.4  2006/05/19 22:48:30  pwagner
! May rename copied SDs
!
! Revision 1.3  2005/09/23 21:01:13  pwagner
! use_wall_clock now a component of time_config
!
! Revision 1.2  2005/06/22 19:27:33  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 1.1  2005/04/29 21:57:22  pwagner
! First commit
!
