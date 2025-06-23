! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=================================
program mangleL2AUX ! mangles L2AUX files by resetting range of maf values
!=================================

   use Dump_0, only: DUMP
   use Hdf, only: DFACC_CREATE, DFACC_RDWR, DFACC_READ, &
     & DFNT_CHAR8, DFNT_FLOAT32, DFNT_INT32, DFNT_FLOAT64
   use HDF5, only: h5fopen_f, h5fclose_f, h5fis_hdf5_f, &
     & H5GCLOSE_F, H5GOPEN_F, H5DOPEN_F, H5DCLOSE_F, h5gcreate_f
   use HDFEOS5, only: HE5T_NATIVE_CHAR
   use L1BData, only: L1BData_T, NAME_LEN, &
     & DeallocateL1BData, Diff, ReadL1BData
   use L2GPData, only: cpL2GPData, L2GPData_T, ReadL2GPData, DestroyL2GPContents, &
     & L2GPNameLen, MAXSWATHNAMESBUFSIZE
   use MACHINE, only: FILSEP, HP, IO_ERROR, GETARG
   use MLSCommon, only: R8
   use MLSFiles, only: FILENOTFOUND, WILDCARDHDFVERSION, &
     & mls_exists, mls_hdf_version, mls_sfstart, mls_sfend, &
     & MLS_IO_GEN_OPENF, MLS_IO_GEN_CLOSEF, &
     & HDFVERSION_4, HDFVERSION_5, MLS_INQSWATH
   use MLSHDF5, only: GetAllHDF5DSNames, saveAsHDF5DS, &
     & IsHDF5AttributePresent, mls_h5open, mls_h5close
   use MLSMessageModule, only: MLSMessageConfig, MLSMSG_Error, MLSMSG_Warning, &
     & MLSMessage
   use MLSStringLists, only: GetStringElement, NumStringElements
   use output_m, only: output
   use PCFHdr, only: GlobalAttributes
   use Time_M, only: Time_Now, USE_WALL_CLOCK
   
   implicit none

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! mangles L2AUX from list of input files to a single output file
! May thus "unsplit" the dgm files

! To use this, copy it into
! mlspgs/tests/l2
! then enter "make depends" followed by "make"

  type options_T
    integer, dimension(2)  :: mifs = 0
    integer, dimension(2)  :: mafs = 0
    logical     :: verbose = .false.
    logical     :: list = .false.
    logical     :: showdiff = .false.
    character(len=255) :: outputFile= 'default.h5'        ! output filename
  end type options_T
  
  type ( options_T ) :: options

! Then run it
! LF95.Linux/test [options] [input files] -o [output file]

  integer, parameter ::          MAXDS = 50
  integer, parameter ::          MAXSDNAMESBUFSIZE = MAXDS*NAME_LEN
  integer, parameter ::          MAXFILES = 100
  logical ::          columnsOnly
  character(len=255) :: filename          ! input filename
  character(len=255), dimension(MAXFILES) :: filenames
  integer            :: n_filenames
  integer     ::  i, count, status, error ! Counting indices & Error flags
  logical     :: is_hdf5
  character (len=MAXSDNAMESBUFSIZE) :: mySdList
  real        :: t1
  real        :: t2
  real        :: tFile
  ! 
  MLSMessageConfig%useToolkit = .false.
  MLSMessageConfig%logFileUnit = -1
  USE_WALL_CLOCK = .true.
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
  else
    call time_now ( t1 )
    if ( options%verbose .and. .not. options%list) &
      & print *, 'Copy l1b data to: ', trim(options%outputFile)
    do i=1, n_filenames
      call time_now ( tFile )
      if ( options%list ) then
        print *, 'DS Names in: ', trim(filenames(i))
        call GetAllHDF5DSNames (trim(filenames(i)), '/', mysdList)
        call dump(mysdList, 'DS names')
      else 
        if ( options%verbose ) then
          print *, 'Mangling from: ', trim(filenames(i))
        endif
        call mangleCopy(trim(filenames(i)), &
        & trim(options%outputFile), (i==1), &
        & HDFVERSION_5, options)
        call sayTime('Mangling this file', tFile)
      endif
    enddo
    if ( .not. options%list) call sayTime('Mangling all files')
  endif
  call mls_h5close(error)
contains
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
      elseif ( filename(1:6) == '-mafs ' ) then
        call getarg ( i+1+hp, number )
        if ( NumStringElements(trim(number), .true.) < 1 ) then
          call output('usage: -mafs must be followed by "maf1,maf2"', advance='yes')
        elseif ( NumStringElements(trim(number), .true.) < 2 ) then
          read(trim(number), *) options%mafs(1)
        else
          read(trim(number), *) options%mafs(1), options%mafs(2)
        endif
        i = i + 1
        exit
      elseif ( filename(1:6) == '-mifs ' ) then
        call getarg ( i+1+hp, number )
        if ( NumStringElements(trim(number), .true.) < 1 ) then
          call output('usage: -mifs must be followed by "mif1,mif2"', advance='yes')
        elseif ( NumStringElements(trim(number), .true.) < 2 ) then
          read(trim(number), *) options%mifs(1)
        else
          read(trim(number), *) options%mifs(1), options%mifs(2)
        endif
        i = i + 1
        exit
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
      elseif ( filename(1:3) == '-d ' ) then
        options%showdiff = .true.
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
      print *,  "Enter the name of the HDFEOS5 L2AUX file. " // &
       &  "The default output file name will be used."
      read(*,'(a)') filename
    endif
    
  end subroutine get_filename
!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'Usage:mangleL2AUX [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options: -f filename     => add filename to list of filenames'
      write (*,*) '                           (can do the same w/o the -f)'
      write (*,*) '          -o ofile        => copy sds to ofile'
      write (*,*) '          -d              => show resulting diffs'
      write (*,*) '          -v              => switch on verbose mode'
      write (*,*) '          -mafs maf1,maf2 => mangle range of mafs'
      write (*,*) '          -mifs mif1,mif2 => mangle range of mifs'
      write (*,*) '          -l              => just list sd names in files'
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

  ! ---------------------- mangleCopy  ---------------------------
  subroutine mangleCopy(file1, file2, create2, hdfVersion, options)
  !------------------------------------------------------------------------

    ! Given file names file1 and file2,
    ! This routine copies all the l1bdata from 1 to 2 after mangling it
    ! If file2 doesn't exist yet, or if create2 is TRUE, it'll create it

    ! Arguments

    character (len=*), intent(in) :: file1 ! Name of file 1
    character (len=*), intent(in) :: file2 ! Name of file 2
    logical, intent(in) :: create2
    integer, intent(in) :: hdfVersion
    type ( options_T ), intent(in) :: options

    ! Local
    integer :: QuantityType
    integer :: sdfid1
    integer :: sdfid2
    integer :: grpid
    integer :: sd_id
    integer :: status
    integer :: the_hdfVersion
    logical :: file_exists
    integer :: file_access
    integer :: listsize
    integer :: noSds
    character (len=MAXSDNAMESBUFSIZE) :: mySdList
    logical, parameter            :: countEmpty = .true.
    logical :: isl1boa
    ! type (L2AUXData_T) :: l2aux
    type(l1bdata_t) :: L1BDATA  ! Result
    type(l1bdata_t) :: L1BDATA2 ! for diff
    integer :: i
    integer :: maf1, maf2
    integer :: mif1, mif2
    integer :: NoMAFs
    character (len=80) :: sdName
    
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

    isl1boa = (index(trim(mysdList), '/GHz') > 0)
    file_exists = ( mls_exists(trim(File2)) == 0 )
    if ( file_exists ) then
      file_access = DFACC_RDWR
    else
      file_access = DFACC_CREATE
    endif
    if ( create2 ) file_access = DFACC_CREATE
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
    sdfId2 = mls_sfstart(trim(file2), file_access, &
              & hdfVersion=hdfVersion)
    if (sdfid2 == -1 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Failed to open l1b ' // trim(File2) )
    end if
    noSds = NumStringElements(trim(mysdList), countEmpty)
    if ( noSds < 1 ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'No sdNames cp to file--unable to count sdNames in ' // trim(mysdList) )
    endif
    ! May need to create groups '/sc', '/GHz', and '/THz'
    if ( isL1boa .and. file_access == DFACC_CREATE ) then
      call h5gcreate_f(sdfId2, '/sc', grpID, status)
      if ( status /= 0 ) then
	     call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	& 'Unable to create group /sc in ' // trim(File2) )
      endif
      call h5gcreate_f(sdfId2, '/GHz', grpID, status)
      if ( status /= 0 ) then
	     call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	& 'Unable to create group /GHz in ' // trim(File2) )
      endif
      call h5gcreate_f(sdfId2, '/THz', grpID, status)
      if ( status /= 0 ) then
	     call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	& 'Unable to create group /THz in ' // trim(File2) )
      endif
    endif
    ! Loop over sdNames in file 1
    do i = 1, noSds
      call GetStringElement (trim(mysdList), sdName, i, countEmpty )
      ! Allocate and fill l2aux
      if ( options%verbose ) print *, 'About to read ', trim(sdName)
      ! call ReadL2AUXData ( sdfid1, trim(sdName), QuantityType, l2aux, &
      !      & checkDimNames=.false., hdfVersion=hdfVersion )
      call ReadL1BData ( sdfid1, trim(sdName), L1bData, NoMAFs, status, &
        & hdfVersion=the_hdfVersion, NEVERFAIL=.true. )
      if ( status /= 0 ) then
        if ( options%verbose ) print *, 'Unknown type for ', trim(sdName)
        if ( options%verbose ) print *, 'status ', status
        cycle
      endif
      if ( options%verbose ) print *, 'Read ', NoMAFs, ' mafs ', L1BData%MaxMIFs, ' mifs'
      ! Mangle the data
      if ( trim(L1BData%data_type) == 'double' ) then
        maf1 = max(L1BData%FirstMAF+1, options%mafs(1))
        maf2 = max(options%mafs(2), maf1)
        mif1 = max(1, options%mifs(1))
        if ( options%mifs(2) < 1 ) then
          mif2 = L1BData%MaxMIFs
        else
          mif2 = min(max(options%mifs(2), mif1), L1BData%MaxMIFs)
        endif
        if ( options%verbose ) then
          print *, 'About to mangle ', trim(sdName)
          print *, 'mif range ', mif1, mif2
          print *, 'maf range ', maf1, maf2
        endif
        L1BData%DpField(:, mif1:mif2, maf1:maf2) = -999.99
      elseif ( options%verbose) then
        call output('sd: ', advance='no')
        call output(trim(sdName), advance='no')
        call output('  Data type: ', advance='no')
        call output(trim(L1BData%data_type), advance='no')
        call output('  not mangled ', advance='yes')
      endif
      ! Write the filled l1bdata to file2
      if ( options%verbose) print *, 'About to write ', trim(sdName)
      call WriteL1BData(l1bdata, sdfid2, status, trim(sdName))
      ! Deallocate memory used by the l2aux
      call DeallocateL1BData ( l1bData )
      ! Show diff between the resulting l1bdata sets
      if ( options%showdiff ) then
        call ReadL1BData ( sdfid1, trim(sdName), L1bData, NoMAFs, status, &
          & hdfVersion=the_hdfVersion, NEVERFAIL=.true. )
        call ReadL1BData ( sdfid2, trim(sdName), L1bData2, NoMAFs, status, &
          & hdfVersion=the_hdfVersion, NEVERFAIL=.true. )
        call diff(L1bData, L1bData2)
        call DeallocateL1BData ( l1bData )
        call DeallocateL1BData ( l1bData2 )
      endif
    enddo
	 call h5gClose_f (grpID, status)
    if ( status /= 0 ) then
	   call MLSMessage ( MLSMSG_Warning, ModuleName, &
       & 'Unable to close group in l2aux file: ' // trim(File1) // ' after cping' )
    endif
	 status = mls_sfend(sdfid1, hdfVersion=the_hdfVersion)
    if ( status /= 0 ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
       & "Unable to close L2aux file: " // trim(File1) // ' after cping')
	 status = mls_sfend(sdfid2, hdfVersion=the_hdfVersion)
    if ( status /= 0 ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
       & "Unable to close L2aux file: " // trim(File2) // ' after cping')
  end subroutine mangleCopy

  ! ---------------------- WriteL1BData  ---------------------------
  subroutine WriteL1BData(l1bdata, fileID, status, sdName)
  !------------------------------------------------------------------------
    ! Arguments
    type(l1bdata_t), intent(in) :: L1BDATA
    integer, intent(in)         :: fileID   ! file handle
    integer, intent(out)        :: status   ! /= 0 if trouble
    character(len=*), intent(in):: sdName
    
    ! Local variables
    integer :: total_DS_size
    ! Executable
    select case (trim(L1BData%data_type))
    case ('double')
      if ( L1BData%trueRank == 3 ) then
        call SaveAsHDF5DS( fileID, trim(sdName), &
          & real( &
          &   L1BData%DpField) )
      elseif ( L1BData%trueRank == 2 ) then
        call SaveAsHDF5DS( fileID, trim(sdName), &
          & real( &
          &   L1BData%DpField(1, :, :)) )
      else
        call SaveAsHDF5DS( fileID, trim(sdName), &
          & real( &
          &   L1BData%DpField(1, 1, :)) )
      endif
    case ('integer')
      if ( L1BData%trueRank == 3 ) then
        ! not implemented yet
        call output('3d integer hdf5 writes not implemented yet', advance='yes')
        ! call SaveAsHDF5DS( fileID, trim(sdName), &
        !   & L1BData%intField )
      elseif ( L1BData%trueRank == 2 ) then
        call SaveAsHDF5DS( fileID, trim(sdName), &
          & L1BData%intField(1, :, :) )
      else
        call SaveAsHDF5DS( fileID, trim(sdName), &
          & L1BData%intField(1, 1, :) )
      endif
    case ('character')
      if ( L1BData%trueRank == 3 ) then
        ! not implemented yet
        call output('3d char hdf5 writes not implemented yet', advance='yes')
        ! call SaveAsHDF5DS( fileID, trim(sdName), &
        !   & L1BData%charField )
      elseif ( L1BData%trueRank == 2 ) then
        call SaveAsHDF5DS( fileID, trim(sdName), &
          & L1BData%charField(1, :, :) )
      else
        call SaveAsHDF5DS( fileID, trim(sdName), &
          & L1BData%charField(1, 1, :) )
      endif
    case default
      call output('sd: ', advance='no')
      call output(trim(sdName), advance='no')
      call output(' Data type: ', advance='no')
      call output(trim(L1BData%data_type), advance='no')
      call output(' not a recognized type ', advance='yes')
    end select
  end subroutine WriteL1BData

!==================
end program mangleL2AUX
!==================

! $Log$
