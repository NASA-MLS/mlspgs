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
program l2gpcat ! catenates split L2GPData files, e.g. dgg
!=================================

   use Hdf, only: DFACC_CREATE, DFACC_RDWR, DFACC_READ
   use HDF5, only: h5fopen_f, h5fclose_f, h5fis_hdf5_f   
   use HDFEOS5, only: HE5T_NATIVE_CHAR
   use L2GPData, only: cpL2GPData, L2GPData_T, ReadL2GPData, DestroyL2GPContents, &
     & L2GPNameLen, MAXSWATHNAMESBUFSIZE
   use MACHINE, only: FILSEP, HP, IO_ERROR, GETARG
   use MLSCommon, only: R8
   use MLSFiles, only: mls_exists, &
     & HDFVERSION_4, HDFVERSION_5, MLS_INQSWATH
   use MLSHDF5, only: mls_h5open, mls_h5close
   use MLSMessageModule, only: MLSMessageConfig
   use MLSStringLists, only: catLists, GetStringElement, GetUniqueList, &
     & Intersection, NumStringElements, RemoveElemFromList, &
     & StringElement, StringElementNum
   use output_m, only: output
   use PCFHdr, only: GlobalAttributes
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
! Catenate L2GPData from list of input files to a single output file
! May thus "unsplit" the dgg files

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it
! LF95.Linux/test [options] [input files] -o [output file]
  type options_T
    logical     :: verbose = .false.
    character(len=255) ::    outputFile= 'default.he5'  ! output filename       
    logical ::               columnsOnly = .false.
    logical ::               noDupSwaths = .false.      ! cp 1st, ignore rest   
    character(len=3) ::      convert= ' '               ! e.g., '425'
    character(len=255) ::    swathNames = ' '           ! which swaths to copy
    character(len=255) ::    rename = ' '               ! how to rename them
    integer, dimension(2) :: freqs = 0                  ! Keep range of freqs   
    integer, dimension(2) :: levels = 0                 ! Keep range of levels   
    integer, dimension(2) :: profiles = 0               ! Keep range of profiles   
  end type options_T
  
  type ( options_T ) :: options

  integer, parameter ::          MAXFILES = 100
  logical, parameter ::          countEmpty = .true.
  logical     :: createdYet
  logical, parameter ::          DEEBUG = .false.
  ! logical ::          columnsOnly
  character(len=255) :: filename          ! input filename
  ! character(len=255) :: outputFile= 'default.he5'        ! output filename
  character(len=255), dimension(MAXFILES) :: filenames
  integer            :: n_filenames
  integer     ::  i, j, status, error ! Counting indices & Error flags
  integer     :: elem
  integer     ::  hdfversion1
  integer     ::  hdfversion2
  logical     :: is_hdf5
  ! logical     :: verbose = .false.
  real        :: t1
  real        :: t2
  real        :: tFile
  character(len=255) ::    rename = ' '               ! how to rename them
  character(len=L2GPNameLen)          :: swath
  character(len=MAXSWATHNAMESBUFSIZE) :: swathList
  character(len=MAXSWATHNAMESBUFSIZE) :: swathList1
  character(len=MAXSWATHNAMESBUFSIZE) :: swathListAll
  integer :: listSize
  integer :: NUMSWATHSPERFILE
  integer :: NUMSWATHSSOFAR
  ! 
  MLSMessageConfig%useToolkit = .false.
  MLSMessageConfig%logFileUnit = -1
  time_config%use_wall_clock = .true.
  CALL mls_h5open(error)
  n_filenames = 0
!   do      ! Loop over input
!     read (*, '(a)') swathList
!     if ( swathList == 'stop' ) stop
!     read (*, '(a)') swathList1
!     swathListAll = catLists(swathList, swathList1)
!     swathList = swathListAll
!     call GetUniqueList(swathList, swathListAll, &
!       & i, countEmpty)
!     print *, 'no unique: ', i
!     print *, 'unique: ', trim(swathListAll)
!     read (*, '(a)') swath
!     call RemoveElemFromList (swathListAll, swathList, trim(swath))
!     print *, 'After removing: ', trim(swathList)
!   enddo
  do      ! Loop over filenames
     call get_filename(filename, n_filenames, options)
     if ( filename(1:1) == '-' ) cycle
     if ( filename == ' ' ) exit
     if ( mls_exists(trim(filename)) /= 0 ) then
       print *, 'Sorry--file not found: ', trim(filename)
       cycle
     endif
     n_filenames = n_filenames + 1
     filenames(n_filenames) = filename
  enddo
  createdYet = .false.
  if ( n_filenames == 0 ) then
    if ( options%verbose ) print *, 'Sorry no input files to copy'
  else
    ! Check that the hdfversions of the input files accord with convert mode
    status = 0
    do i=1, n_filenames
     call h5fis_hdf5_f(filenames(i), is_hdf5, error)
     select case (options%convert)
     case ('425')
       if ( is_hdf5 ) then
         print *, 'Sorry--not recognized as hdf4 file: ', trim(filenames(i))
         status = 1
         cycle
       else
         hdfVersion1 = HDFVERSION_4
         hdfVersion2 = HDFVERSION_5
       endif
     case ('524')
       if ( .not. is_hdf5 ) then
         print *, 'Sorry--not recognized as hdf5 file: ', trim(filenames(i))
         status = 1
         cycle
       else
         hdfVersion1 = HDFVERSION_5
         hdfVersion2 = HDFVERSION_4
       endif
     case default
       if ( .not. is_hdf5 ) then
         hdfVersion1 = HDFVERSION_4
       else
         hdfVersion1 = HDFVERSION_5
       endif
       hdfVersion2 = hdfVersion1
     end select
    enddo
    call time_now ( t1 )
    swathListAll = ''
    numswathssofar = 0
    if ( options%verbose ) print *, 'Copy l2gp data to: ', trim(options%outputFile)
    do i=1, n_filenames
      call time_now ( tFile )
      if ( options%verbose ) print *, 'Copying from: ', trim(filenames(i))
      if ( options%noDupSwaths .or. options%swathNames /= ' ' ) then
        numswathsperfile = mls_InqSwath ( trim(filenames(i)), &
          & swathList, listSize, hdfVersion=hdfVersion1)
        if ( DEEBUG ) then
          print *, 'swaths in file'
          print *, trim(swathList)
        endif
        if ( DEEBUG ) then
          print *, 'all swaths'
          print *, trim(swathListAll)
        endif
        if ( options%swathNames /= ' ' ) then
          swathList1 = swathList
          swathList = Intersection( options%swathNames, swathList1 )
          if ( swathList == ' ' ) cycle
          rename = ' '
          do j=1, NumStringElements( swathList, countEmpty )
            call GetStringElement( swathList, swath, j, countEmpty )
            elem = StringElementNum( options%swathNames, swath, countEmpty )
            rename = catLists( rename, &
              & StringElement( options%rename, elem, countempty ) )
          enddo
        elseif ( numswathssofar > 0 ) then
          ! Remove any duplicates
          do j=1, numswathssofar
            call GetStringElement(swathListAll, swath, j, countEmpty)
            swathList1 = swathList
            call RemoveElemFromList (swathList1, swathList, trim(swath))
            ! Crude hAck--really should fix removeElem procedure
            if ( swathList(1:1) == ',' ) then
              swathList1 = swathList(2:)
              swathList = swathList1
            endif
          enddo
          if ( DEEBUG ) then
            print *, 'swaths to cp'
            print *, trim(swathList)
          endif
        endif
        if ( any( (/options%freqs(2), options%levels(2), &
          & options%profiles(2)/) > 0 ) &
          & ) then
          call cpL2GPData(trim(filenames(i)), &
          & trim(options%outputFile), create2=.not. createdYet, &
          & hdfVersion1=hdfVersion1, hdfVersion2=hdfVersion2, &
          & swathList=trim(swathList), rename=rename, &
          & notUnlimited=.true., andGlAttributes=.true., &
          & rFreqs=options%freqs, rLevels=options%levels, rTimes=options%profiles)
        else
          call cpL2GPData(trim(filenames(i)), &
          & trim(options%outputFile), create2=.not. createdYet, &
          & hdfVersion1=hdfVersion1, hdfVersion2=hdfVersion2, &
          & swathList=trim(swathList), rename=rename, &
          & notUnlimited=.true., andGlAttributes=.true.)
        endif
        swathList1 = swathListAll
        swathListAll = catlists(swathList1, swathList)
        numswathssofar = NumStringElements(swathListAll, countEmpty)
      else
        if ( any( (/options%freqs(2), options%levels(2), options%profiles(2)/) &
          & > 0 ) ) then
          call cpL2GPData(trim(filenames(i)), &
          & trim(options%outputFile), create2=.not. createdYet, &
          & hdfVersion1=hdfVersion1, hdfVersion2=hdfVersion2, &
          & notUnlimited=.true., andGlAttributes=.true., &
          & rFreqs=options%freqs, rLevels=options%levels, rTimes=options%profiles)
        else
          call cpL2GPData(trim(filenames(i)), &
          & trim(options%outputFile), create2=.not. createdYet, &
          & hdfVersion1=hdfVersion1, hdfVersion2=hdfVersion2, &
          & notUnlimited=.true., andGlAttributes=.true.)
        endif
      endif
      call sayTime('copying this file', tFile)
      createdYet = .true.
    enddo
    call sayTime('copying all files')
  endif
  call mls_h5close(error)
contains
!------------------------- get_filename ---------------------
    subroutine get_filename(filename, n_filenames, options)
    ! Added for command-line processing
     character(LEN=255), intent(out) :: filename          ! filename
     integer, intent(in)             :: n_filenames
     type ( options_T ), intent(inout) :: options
     ! character(LEN=*), intent(inout) :: outputFile        ! output filename
     ! logical, intent(inout)          :: verbose
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
      elseif ( filename(1:5) == '-425 ' ) then
        options%convert = '425'
        exit
      elseif ( filename(1:5) == '-524 ' ) then
        options%convert = '524'
        exit
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
        exit
      elseif ( filename(1:3) == '-no' ) then
        options%noDupSwaths = .true.
        exit
      else if ( filename(1:3) == '-f ' ) then
        call getarg ( i+1+hp, filename )
        i = i + 1
        exit
      else if ( filename(1:3) == '-s ' ) then
        call getarg ( i+1+hp, options%swathNames )
        i = i + 1
        exit
      else if ( filename(1:3) == '-r ' ) then
        call getarg ( i+1+hp, options%rename )
        i = i + 1
        exit
      elseif ( filename(1:5) == '-freq' ) then
        call igetarg ( i+1+hp, options%freqs(1) )
        i = i + 1
        call igetarg ( i+1+hp, options%freqs(2) )
        i = i + 1
        exit
      elseif ( filename(1:6) == '-level' ) then
        call igetarg ( i+1+hp, options%levels(1) )
        i = i + 1
        call igetarg ( i+1+hp, options%levels(2) )
        i = i + 1
        exit
      elseif ( filename(1:5) == '-prof' ) then
        call igetarg ( i+1+hp, options%profiles(1) )
        i = i + 1
        call igetarg ( i+1+hp, options%profiles(2) )
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
      print *,  "Enter the name of the HDFEOS4 or 5 L2GP file. " // &
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
      write (*,*) ' Options: -f filename   => add filename to list of filenames'
      write (*,*) '                           (can do the same w/o the -f)'
      write (*,*) '          -o ofile      => copy swaths to ofile'
      write (*,*) '          -425          => convert from hdf4 to hdf5'
      write (*,*) '          -524          => convert from hdf5 to hdf4'
      write (*,*) '          -v            => switch on verbose mode'
      write (*,*) '          -nodup        => if dup swath names, cp 1st only'
      write (*,*) '          -freqs m n    => keep only freqs in range m n'
      write (*,*) '          -levels m n   => keep only levels in range m n'
      write (*,*) '          -profiles m n => keep only profiles in range m n'
      write (*,*) '          -s name1,name2,..'
      write (*,*) '             => copy only swaths so named; otherwise all'
      write (*,*) '          -r rename1,rename2,..'
      write (*,*) '             => if and how to rename the copied swaths'
      write (*,*) '          -h            => print brief help'
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
!------------------------- igetarg ---------------------
  subroutine igetarg ( pos, iarg )
   integer, intent(in) :: pos
   integer, intent(out) :: iarg
   character(len=16) :: arg
   call getarg ( pos, arg )
   read(arg, *) iarg
  end subroutine igetarg

!==================
end program L2GPcat
!==================

! $Log$
! Revision 1.9  2006/05/19 20:55:36  pwagner
! May rename copied swaths
!
! Revision 1.8  2006/04/06 23:04:21  pwagner
! Optionally cp only ranges of freq, level, profile
!
! Revision 1.7  2005/10/29 00:13:56  pwagner
! Removed unused procedures from use statements
!
! Revision 1.6  2005/09/23 21:01:13  pwagner
! use_wall_clock now a component of time_config
!
! Revision 1.5  2005/06/22 19:27:33  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 1.4  2004/12/06 19:13:12  pwagner
! With -nodup option ignores dup swath names after 1st
!
! Revision 1.3  2004/08/07 00:15:55  pwagner
! All stringlist stuff was moved from mlsstrings to mlsstringlists
!
! Revision 1.2  2004/05/06 21:50:48  pwagner
! Uses mls_h5open/close
!
! Revision 1.1  2004/04/30 18:54:22  pwagner
! First commit
!
