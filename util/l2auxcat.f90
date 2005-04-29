! Copyright (c) 2005, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contracts NAS7-1407/NAS7-03001 is acknowledged.

!=================================
program L2AUXcat ! catenates split L2AUX files, e.g. dgm
!=================================

   use Dump_0, only: DUMP
   use Hdf, only: DFACC_CREATE, DFACC_RDWR, DFACC_READ
   use HDF5, only: h5fis_hdf5_f
   use HDFEOS5, only: HE5T_NATIVE_CHAR
   use L2AUXData, only: L2AUXDATA_T, MAXSDNAMESBUFSIZE, &
     & cpL2AUXData, WriteL2AUXData
   use MACHINE, only: HP, GETARG
   use MLSFiles, only: HDFVERSION_5
   use MLSHDF5, only: GetAllHDF5DSNames, mls_h5open, mls_h5close
   use MLSMessageModule, only: MLSMessageConfig
   use MLSStringLists, only: catLists, GetStringElement, NumStringElements, &
     & RemoveElemFromList
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
! Catenate L2AUX from list of input files to a single output file
! May thus "unsplit" the dgm files

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it
! LF95.Linux/test [options] [input files] -o [output file]

  integer, parameter ::          MAXFILES = 100
  type options_T
    logical     :: verbose = .false.
    character(len=255) :: outputFile= 'default.h5'        ! output filename
    logical ::          noDupDSNames = .false.            ! cp 1st, ignore rest
    logical     :: list = .false.
    character(len=255), dimension(MAXFILES) :: filenames
  end type options_T
  
  type ( options_T ) :: options

  logical, parameter :: COUNTEMPTY = .false.
  logical, parameter :: DEEBUG = .false.
  character(len=255) :: filename          ! input filename
  integer            :: n_filenames
  integer     ::  i, j, status, error ! Counting indices & Error flags
  logical     :: is_hdf5
  character (len=MAXSDNAMESBUFSIZE) :: mySdList
  integer :: numdsets
  integer :: numdsetssofar
  character(len=MAXSDNAMESBUFSIZE) :: sdListAll
  character(len=255) :: sdName
  character (len=MAXSDNAMESBUFSIZE) :: tempSdList
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
     options%filenames(n_filenames) = filename
  enddo
  if ( n_filenames == 0 ) then
    if ( options%verbose ) print *, 'Sorry no input files supplied'
  else
    call time_now ( t1 )
    if ( options%verbose .and. .not. options%list) print *, 'Copy l2aux data to: ', trim(options%outputFile)
    numdsetssofar = 0
    sdListAll = ''
    do i=1, n_filenames
      call time_now ( tFile )
      if ( options%list ) then
        print *, 'DS Names in: ', trim(options%filenames(i))
        call GetAllHDF5DSNames (trim(options%filenames(i)), '/', mysdList)
        call dump(mysdList, 'DS names')
      else 
        if ( options%verbose ) then
          print *, 'Copying from: ', trim(options%filenames(i))
        endif
        if ( options%noDupDSNames) then
          call GetAllHDF5DSNames (trim(options%filenames(i)), '/', mysdList)
          numdsets = NumStringElements(trim(mysdList), COUNTEMPTY)
          if ( numdsetssofar > 0 ) then
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
          call cpL2AUXData(trim(options%filenames(i)), &
          & trim(options%outputFile), create2=(i==1), &
          & hdfVersion=HDFVERSION_5, sdList=mysdList)
          tempSdList = sdListAll
          sdListAll = catlists(tempSdList, mysdList)
          numdsetssofar = NumStringElements(sdListAll, countEmpty)
        else
          call cpL2AUXData(trim(options%filenames(i)), &
          & trim(options%outputFile), create2=(i==1), &
          & hdfVersion=HDFVERSION_5)
        endif
        call sayTime('copying this file', tFile)
      endif
    enddo
    if ( .not. options%list ) call sayTime('copying all files')
  endif
  call mls_h5close(error)
contains
!------------------------- get_filename ---------------------
    subroutine get_filename(filename, n_filenames, options)
    ! Added for command-line processing
     character(LEN=255), intent(out)   :: filename          ! filename
     integer, intent(in)               :: n_filenames
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
      elseif ( filename(1:3) == '-no' ) then
        options%noDupDSNames = .true.
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
      & 'Usage:l2auxcat [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options: -f filename => add filename to list of filenames'
      write (*,*) '                         (can do the same w/o the -f)'
      write (*,*) '          -o ofile    => copy data sets to ofile'
      write (*,*) '          -v          => switch on verbose mode'
      write (*,*) '          -l          => just list l2aux names in files'
      write (*,*) '          -nodup      => if dup dataset names, cp 1st only'
      write (*,*) '          -h          => print brief help'
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
