program TEST

  use GETCF_M, only: GetCF, InitGetCF
  use MACHINE ! At least HP, for command lines, and maybe GETARG
  use MLSCF, only: MLSCF_T
  use TOGGLES, only: CON, GEN, LEX, PAR, SYN, TAB, TOGGLE

  logical :: DO_DUMP = .false.     ! Dump declaration table
  logical :: DO_DUMP_EARLY = .false.    ! Dump declaration table before check
  logical :: DO_LISTING = .false.  ! List input
  logical :: DUMP_TREE = .false.   ! Dump tree after parsing
  integer :: ERROR                 ! Error flag from GetCF
  integer :: I                     ! counter for command line arguments
  integer :: J                     ! index within option
  type(mlscf_t) :: CF_DATA
  character(len=80) :: LINE        ! to read command line arguments

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

  call getCF ( cf_data, error, inUnit=-1, listing=do_listing, &
    & dump=do_dump, dumpEarly=do_dump_early, dumpTables=.true. )

end program TEST

! $Log$
! Revision 2.1  2000/09/29 23:30:09  vsnyder
! Revised to account for getL2CF_m
!
! Revision 2.0  2000/09/05 17:41:51  dcuddy
! Change revision to 2.0
!
! Revision 1.1  2000/07/06 01:43:12  vsnyder
! Initial check-in
!
