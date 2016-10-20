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
program l2pcdump ! dumps datasets, attributes from l2pc files
!=================================

   use Allocate_Deallocate, only: allocate_test, deallocate_test
   use DECLARATION_TABLE, only: ALLOCATE_DECL, DEALLOCATE_DECL
   use Dump_0, only: Dump
   use dump_options, only: Intplaces
   use Hdf, only: DFACC_READ
   use HDF5, only: h5fis_hdf5_f, h5gclose_f, h5gopen_f
   use HIGHOUTPUT, only: DUMP
   use INIT_TABLES_MODULE, only: Init_tables
   use Intrinsic, only: l_hdf
   use L2PC_m, only: DestroyL2PC, Dump, L2PCDatabase, &
     & ReadCompleteHDF5L2PCFile, PopulateL2PCBin
   use LEXER_CORE, only: INIT_LEXER
   use MACHINE, only: HP, GETARG
   use MLSCommon, only: R8, MLSFile_T
   use MLSFiles, only: FILENOTFOUND, &
     & InitializeMLSFile, mls_exists, mls_sfstart, mls_sfend, &
     & HDFVERSION_5, mls_hdf_version, WILDCARDHDFVERSION
   use MLSHDF5, only: DumpHDF5Attributes, DumpHDF5DS, &
     & GetAllHDF5AttrNames, GetAllHDF5DSNames, &
     & mls_h5open, mls_h5close
   use MLSMessageModule, only: MLSMessageConfig, MLSMSG_Error, MLSMSG_Warning, &
     & dumpConfig, MLSMessage
   use MLSStats1, only: FILLVALUERELATION, Stat_T, dump, STATISTICS
   use MLSStringLists, only: catLists, GetStringElement, List2Array, &
     & NumStringElements, StringElementNum
   use MLSStrings, only: lowercase, trim_safe
   use output_m, only: outputOptions, output
   use PrintIt_m, only: Set_Config
   use Time_M, only: Time_Now, time_config
   use TOGGLES, only: SWITCHES
   use TREE, only: ALLOCATE_TREE, DEALLOCATE_TREE
   
   implicit none

!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

! Brief description of program
! dumps datasets or attributes from l2pc files

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"
! Then run it
! LF95.Linux/test [options] [input files]

  type options_T
    integer             :: details            = 3       ! Print matrices and hessians
    character(len=32)  ::  dumpOptions        = ''
    logical             :: laconic            = .false. ! Print contents only
    logical             :: verbose            = .false. ! Print (lots) extra
    logical             :: la                 = .false.
    logical             :: ls                 = .false.
    logical             :: timereads          = .false. ! Just time how long to read
    character(len=32), dimension(3)  :: blockNames  = (/ '*', '*', '*' /) ! If not all
    character(len=128)  :: root        = '/'
    character(len=1024) :: attributes = ''
  end type options_T
  
  type ( options_T ) :: options
  integer, parameter ::          MAXFILES = 100
  integer, parameter ::          hdfVersion = HDFVERSION_5
  ! character(len=8)   ::          dumpOptions
  character(len=255) ::          filename          ! input filename
  character(len=255), dimension(MAXFILES) :: filenames
  integer     ::                 i, status, error ! Counting indices & Error flags
  logical     ::                 is_hdf5
  type(MLSFile_T), target ::     MLSFileT
  type(MLSFile_T), pointer ::    MLSFile
  integer            ::          n_filenames
  integer     ::                 j, j1, j2
  integer     ::                 sdfid1
  real        ::                 t1
  real        ::                 t2
  real        ::                 tFile
  ! 
  call set_config ( useToolkit = .false., logFileUnit = -1 )
  time_config%use_wall_clock = .true.
  INTPLACES = '8'
  CALL mls_h5open(error)
  call init_lexer ( n_chars=80000, n_symbols=4000, hash_table_size=611957 )
  call allocate_decl ( ndecls=8000 )
  call allocate_tree ( n_tree=2000000 )
  call init_tables
  switches = 'hess,nl2cf'
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
  ! dumpOptions = '-'
  ! if ( options%laconic ) dumpOptions = trim(dumpOptions) // 'l'
  call time_now ( t1 )
  j2 = 0
  do i=1, n_filenames
    j1 = j2 + 1
    call time_now ( tFile )
    if ( options%verbose ) then
      print *, 'Reading from: ', trim(filenames(i))
    endif
    status = InitializeMLSFile ( MLSFileT, type=l_hdf, access=DFACC_READ, &
     & name=filenames(i), HDFVersion=hdfVersion )
    MLSFile => MLSFileT
    call ReadCompleteHDF5L2PCFile ( MLSFile, 0, shallow=.true. )
    j2 = size( L2PCDatabase )
    if ( options%verbose ) print *, 'Size of l2pc database: ', j2
    do j = j1, j2
      print *, 'j: ', j
      call PopulateL2PCBin( j )
      ! Maybe it can't print name?
      L2PCDataBase(j)%j%name = 0
      L2PCDataBase(j)%h%name = 0
      print *, 'About to dump'
      call Dump( L2PCDataBase(j), &
        & details=options%details, onlyTheseBlocks=options%blockNames, &
        & options=options%dumpOptions )
      print *, 'About to destroy'
      call DestroyL2PC( L2PCDataBase(j) )
    enddo
    ! print *, 'About to mls_sfend on ', sdfid1
    if ( .not. options%laconic ) call sayTime('reading this file', tFile)
  enddo
  call deallocate_decl
  call deallocate_tree
  if ( .not. ( options%laconic .or. options%la .or. options%ls ) ) &
    &  call sayTime('reading all files')
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
     print *, 'details             ', options%details
     print *, 'blockNames          ', trim(options%blockNames(1))
     print *, '                    ', trim(options%blockNames(2))
     print *, '                    ', trim(options%blockNames(3))
     print *, 'dumpOptions         ', options%dumpOptions
     print *, 'laconic?            ', options%laconic
     print *, 'verbose?            ', options%verbose
     print *, 'list attributes  ?  ', options%la   
     print *, 'list datasets  ?    ', options%ls
     print *, 'just time reads?    ', options%timereads
     print *, 'num files           ', n_filenames
     print *, 'switches            ', trim_safe(switches)
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
     character(LEN=96)       :: blocknames
  ! Get inputfile name, process command-line args
  ! (which always start with -)
    do
      call getarg ( i+hp, filename )
      ! print *, i, ' th Arg: ', trim(filename)
      error = 0
      if ( filename(1:1) /= '-' ) exit
      if ( filename(1:3) == '-h ' ) then
        call print_help
      elseif ( filename(1:9) == '-details ' ) then
        call getarg ( i+1+hp, number )
        read( number, * ) options%details
        i = i + 1
        exit
      else if ( filename(1:3) == '-d ' ) then
        call getarg ( i+1+hp, options%dumpOptions )
        i = i + 1
      elseif ( filename(1:4) == '-lac' ) then
        options%laconic = .true.
        exit
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
        exit
      elseif ( filename(1:4) == '-la ' ) then
        options%la = .true.
        exit
      elseif ( filename(1:4) == '-ls ' ) then
        options%ls = .true.
        exit
      else if ( filename(1:3) == '-r ' ) then
        call getarg ( i+1+hp, options%root )
        i = i + 1
        exit
      else if ( filename(1:3) == '-bl' ) then
        call getarg ( i+1+hp, blockNames )
        call List2Array ( blockNames, options%blockNames, .true. )
        i = i + 1
        exit
      else if ( filename(1:3) == '-t ' ) then
        options%timereads = .true.
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
      print *,  "Enter the name of the HDF5 l2pc file. "
      read(*,'(a)') filename
    endif
    
  end subroutine get_filename
!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'Usage:l2pcdump [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options: -f filename     => add filename to list of filenames'
      write (*,*) '                           (can do the same w/o the -f)'
      write (*,*) '          -details number => set dump details to number'
      write (*,*) '                            1 dump matrices only'
      write (*,*) '                            2 dump hessians only'
      write (*,*) '                            3 dump both matrices and hessians'
      write (*,*) '                            (defaults to 3)'
      write (*,*) '          -d options      => pass options to dump routines'
      write (*,*) '          -lac            => switch on laconic mode'
      write (*,*) '          -v              => switch on verbose mode'
      write (*,*) '          -la             => just list attribute names in files'
      write (*,*) '          -ls             => just list sd names in files'
      write (*,*) '          -A              => dump all attributes'
      write (*,*) '          -D              => dump all datasets (default)'
      write (*,*) '          -nA             => do not dump attributes (default)'
      write (*,*) '          -nD             => do not dump datasets'
      write (*,*) '          -bl[ocks] b1,b2,.. => dump just blocks named b1,b2,..'
      write (*,*) '          -h              => print brief help'
      write (*,*) 'Note about -d options:'
      write (*,*) '  any of {mh*b[]d[]} where the "[]" mean themselves here'
      write (*,*) '  m: dump only matrices                   '
      write (*,*) '  h: dump only hessians                   '
      write (*,*) '  *: dump matrices and hessians           '
      write (*,*) '  b[HCN]: dump only HCN blocks            '
      write (*,*) '  d[dopts]: pass dopts when dumping arrays'
      write (*,*) '  default is *                            '
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
end program l2pcdump
!==================

! $Log$
! Revision 1.5  2014/01/09 00:31:26  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 1.4  2013/08/23 02:51:48  vsnyder
! Move PrintItOut to PrintIt_m
!
! Revision 1.3  2011/02/18 23:10:33  pwagner
! Passes -d opts to dump routines
!
! Revision 1.2  2010/09/03 22:13:18  pwagner
! May dump blocks by name
!
! Revision 1.1  2010/08/14 00:11:59  pwagner
! First commit
!
