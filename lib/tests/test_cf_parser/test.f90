! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

program TEST

! A test harness for the MLS CF parser.  It also illustrates how to use
! it, including how to cause the input to come from a Fortran unit instead
! of standard input.

  use GETCF_M, only: GetCF, InitGetCF
  use MACHINE ! At least HP, for command lines, IO_ERROR, and maybe GETARG
  use MLSCF, only: MLSCF_T
  use OUTPUT_M, only: PRUNIT
  use TOGGLES, only: CON, GEN, LEX, PAR, SYN, TAB, TOGGLE

  type(mlscf_t) :: CF_DATA
  logical :: DO_DUMP = .false.     ! Dump declaration table
  logical :: DO_DUMP_EARLY = .false.    ! Dump declaration table before check
  logical :: DO_LISTING = .false.  ! List input
  logical :: DUMP_TREE = .false.   ! Dump tree after parsing
  integer :: ERROR                 ! Error flag from GetCF
  integer :: I                     ! counter for command line arguments
  integer :: INUNIT = -1           ! Fortran unit number for input -1 = stdin
  integer :: J                     ! index within option
  character(len=80) :: LINE        ! to read command line arguments
  integer :: STATUS

!---------------------------- RCS Ident Info -------------------------------
  character (len=256) :: Id = &
       "$Id$"
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  call InitGetCF

  i = 1+hp
  do ! Process the command line options to set toggles
    call getarg ( i, line )
    if ( line(1:1) == '-' ) then
      do j = 2, len(line)
        if ( line(j:j) == 'c' ) then
          toggle(con) = .true.
        else if ( line(j:j) == 'g' ) then
          toggle(gen) = .true.
        else if ( line(j:j) == 'l' ) then
          toggle(lex) = .true.
        else if ( line(j:j) == 'p' ) then
          toggle(par) = .true.
        else if ( line(j:j) == 'a' ) then
          toggle(syn) = .true.
        else if ( line(j:j) == 'A' ) then
          dump_tree = .true.
        else if ( line(j:j) == 'd' ) then
          do_dump = .true.
        else if ( line(j:j) == 'D' ) then
          do_dump_early = .true.
        else if ( line(j:j) == 't' ) then
          toggle(tab) = .true.
        else if ( line(j:j) == 'v' ) then
          do_listing = .true.
        end if
      end do
    else
  exit
    end if
    i = i + 1
  end do

  if ( line /= ' ' ) then     ! Process input file name
    open ( 98, file=trim(line), form='formatted', iostat=status )
    if ( status /= 0 ) then
      call io_error ( 'While opening input file', status, trim(line) )
      stop
    end if
    inUnit = 98
    call getarg ( i+1, line )
    if ( line /= ' ' ) then   ! Process output file name
      open ( 99, file=trim(line), form='formatted', iostat=status )
      if ( status /= 0 ) then
        call io_error ( 'While opening output file', status, trim(line) )
        stop
      end if
      prUnit = 99
    end if
  end if

  call getCF ( cf_data, error, inUnit=inUnit, listing=do_listing, &
    & dump=do_dump, dumpEarly=do_dump_early, dumpTree=dump_tree, &
    & dumpTables=.true. )

end program TEST

! $Log$
! Revision 1.1  2000/10/11 20:54:40  vsnyder
! Initial addition
!
! Revision 2.4  2000/10/03 01:39:10  vsnyder
! Add the copyright notice.
!
! Revision 2.3  2000/10/03 01:09:19  vsnyder
! Add getting input and output files from command line, and telling the
! parser to use them.
!
! Revision 2.2  2000/10/03 00:54:44  vsnyder
! Revised to account for changing getL2CF_m.f90 to getCF_m.f90
!
! Revision 2.1  2000/09/29 23:30:09  vsnyder
! Revised to account for getL2CF_m
!
! Revision 2.0  2000/09/05 17:41:51  dcuddy
! Change revision to 2.0
!
! Revision 1.1  2000/07/06 01:43:12  vsnyder
! Initial check-in
