! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=================================
program L2GPcat ! catenates split L2GPData files, e.g. dgg
!=================================

   use Hdf, only: DFACC_CREATE, DFACC_RDWR, DFACC_READ
   use HDF5, only: h5fopen_f, h5fclose_f, h5fis_hdf5_f   
   use HDFEOS5, only: HE5T_NATIVE_CHAR
   use L2GPData, only: cpL2GPData, L2GPData_T, ReadL2GPData, DestroyL2GPContents, &
     & L2GPNameLen, MAXSWATHNAMESBUFSIZE
   use MACHINE, only: FILSEP, HP, IO_ERROR, GETARG
   use MLSCommon, only: R8
   use MLSFiles, only: MLS_IO_GEN_OPENF, MLS_IO_GEN_CLOSEF, &
     & HDFVERSION_4, HDFVERSION_5, MLS_INQSWATH
   use MLSMessageModule, only: MLSMessageConfig
   use MLSStrings, only: GetStringElement, NumStringElements
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
! Catenate L2GPData from list of input files to a single output file
! May thus "unsplit" the dgg files

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it
! LF95.Linux/test [options] [input files] -o [output file]

  integer, parameter ::          MAXFILES = 100
  logical ::          columnsOnly
  character(len=255) :: filename          ! input filename
  character(len=255) :: outputFile= 'default.he5'        ! output filename
  character(len=255), dimension(MAXFILES) :: filenames
  integer            :: n_filenames
  integer     ::  i, count, status, error ! Counting indices & Error flags
  logical     :: is_hdf5
  logical     :: verbose = .false.
  real        :: t1
  real        :: t2
  real        :: tFile
  ! 
  MLSMessageConfig%useToolkit = .false.
  MLSMessageConfig%logFileUnit = -1
  USE_WALL_CLOCK = .true.
  CALL h5open_f(error)
  n_filenames = 0
  do      ! Loop over filenames
     call get_filename(filename, n_filenames, outputFile, verbose)
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
    if ( verbose ) print *, 'Sorry no input files to copy'
  else
    call time_now ( t1 )
    if ( verbose ) print *, 'Copy l2gp data to: ', trim(outputFile)
    do i=1, n_filenames
      call time_now ( tFile )
      if ( verbose ) print *, 'Copying from: ', trim(filenames(i))
      call cpL2GPData(trim(filenames(i)), &
        & trim(outputFile), create2=(i==1), &
        & hdfVersion=HDFVERSION_5, notUnlimited=.true., andGlAttributes=.true.)
      call sayTime('copying this file', tFile)
    enddo
    call sayTime('copying all files')
  endif
  call h5close_f(error)
contains
!------------------------- get_filename ---------------------
    subroutine get_filename(filename, n_filenames, outputFile, verbose)
    ! Added for command-line processing
     character(LEN=255), intent(out) :: filename          ! filename
     integer, intent(in)             :: n_filenames
     character(LEN=*), intent(inout) :: outputFile        ! output filename
     logical, intent(inout)          :: verbose
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
        call getarg ( i+1+hp, outputFile )
        i = i + 1
        exit
      elseif ( filename(1:3) == '-v ' ) then
        verbose = .true.
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
      print *,  "Enter the name of the HDFEOS5 L2GP file. " // &
       &  "The default output file name will be used."
      read(*,'(a)') filename
    endif
    
  end subroutine get_filename
!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'Usage:l2gpcat [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options: -f filename => add filename to list of filenames'
      write (*,*) '                         (can do the same w/o the -f)'
      write (*,*) '          -o ofile    => copy swaths to ofile'
      write (*,*) '          -v          => switch on verbose mode'
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
end program L2GPcat
!==================

! $Log$
! Revision 1.1  2004/02/25 00:03:34  pwagner
! First commit
!

