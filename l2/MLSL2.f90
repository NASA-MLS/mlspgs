! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

program MLSL2
  use DECLARATION_TABLE, only: ALLOCATE_DECL, DUMP_DECL
  use INIT_TABLES_MODULE, only: INIT_TABLES, LIT_INDICES
  use LEXER_CORE, only: INIT_LEXER
  use MACHINE ! At least HP for command lines, and maybe GETARG, too
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use OUTPUT_M, only: OUTPUT, PRUNIT
  use Open_Init, only: CloseMLSCF, OpenMLSCF
  use PARSER, only: CONFIGURATION
  use STRING_TABLE, only: DO_LISTING, INUNIT
  use TOGGLES, only: CON, GEN, LEVELS, LEX, PAR, SYN, TAB, TOGGLE
  use TREE, only: ALLOCATE_TREE, PRINT_SUBTREE
  use TREE_CHECKER, only: CHECK_TREE
  use TREE_WALKER, only: WALK_TREE_TO_DO_MLS_L2
  use UNITS, only: INIT_UNITS

  implicit NONE

  integer, parameter :: L2CF_UNIT = 20  ! Unit # if L2CF is read by Fortran

  logical :: DO_DUMP = .false.     ! Dump declaration table
  logical :: DUMP_TREE = .false.   ! Dump tree after parsing
  integer :: ERROR                 ! Error flag from check_tree
  integer :: FIRST_SECTION         ! Index of son of root of first n_cf node
  integer :: I                     ! counter for command line arguments
  integer :: J                     ! index within option
  character(len=80) :: LINE        ! Into which is read the command args
  logical :: PCF = .false.         ! Open L2CF using PCF
  integer :: ROOT                  ! of the abstract syntax tree
  integer :: STATUS                ! From OPEN

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: Id = & 
     "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

! Initialize the lexer, symbol table, and tree checker's tables:
  call init_lexer ( n_chars=10000, n_symbols=1000, hash_table_size=1003 )
  call allocate_decl ( ndecls=1000 )
  call allocate_tree ( n_tree=10000 )
  call init_tables
  call init_units ( lit_indices )

  i = 1+hp
  do ! Process the command line options to set toggles
    call getarg ( i, line )
    if ( line(1:2) == '--' ) then       ! "word" options
      if ( line(3:6) == 'pcf ' ) then
        pcf = .true.
      else if ( line(3:) == ' ' ) then  ! "--" means "no more options"
        i = i + 1
        call getarg ( i, line )
  exit
      else
        print *, 'unrecognized option ', trim(line), ' ignored.'
      end if
    else if ( line(1:1) == '-' ) then   ! "letter" options
      j = 1
      do while ( j < len(line) )
        j = j + 1
        select case ( line(j:j) )
        case ( ' ' )
      exit
        case ( 'A' )
          dump_tree = .true.
        case ( 'a' )
          toggle(syn) = .true.
        case ( 'c' )
          toggle(con) = .true.
        case ( 'd' )
          do_dump = .true.
        case ( 'g' )
          toggle(gen) = .true.
          if ( j < len(line) ) then
            if ( line(j+1:j+1) >= '0' .and. line(j+1:j+1) <= '9' ) then
              j = j + 1
              levels(gen) = ichar(line(j:j)) - ichar('0')
            end if
          end if
        case ( 'h', 'H', '?' )
          call getarg ( 0+hp, line )
          print *, 'Usage: ', trim(line), ' [options] [--] [L2CF-name]'
          print *, ' Options:'
          print *, '  -A: Dump the un-decorated abstract syntax tree.'
          print *, '  -a: Dump the decorated type-checked abstract syntax tree.'
          print *, '  -c: Trace expression evaluation and tree decoration.'
          print *, '  -d: Dump the declaration table after type checking'
          print *, '  -g[digit]: Trace "generation".  Bigger digit means ', &
          &                      'more output.'
          print *, '  -l: Trace lexical analysis.'
          print *, '  -M: Send output through MLSMessage.'
          print *, '  -p: Trace parsing.'
          print *, '  -t: Trace declaration table construction.'
          print *, '  -v: List the configuration file.'
          print *, '  The above options can be concatenated after one hyphen.'
          print *, '  Options a, c, g1, l, p and t can be toggled in the ', &
          &          'configuration file by'
          print *, '  @A, @C, @G, @L, @P and @S respectively.  @T in the ', &
          &          'configuration file'
          print *, '  dumps the string table at that instant.'
          stop
        case ( 'l' )
          toggle(lex) = .true.
        case ( 'M' )
          prunit = -2
        case ( 'p' )
          toggle(par) = .true.
        case ( 't' )
          toggle(tab) = .true.
        case ( 'v' )
          do_listing = .true.
        case default
          print *, 'Unrecognized option -', line(j:j), ' ignored.'
        end select
      end do
    else    
  exit
    end if
    i = i + 1
  end do

! Parse the L2CF, producing an abstract syntax tree
  if ( pcf ) then
    call openMLSCF
  else if ( line /= ' ' ) then
    open ( l2cf_unit, file=line, status='old', &
      & form='formatted', access='sequential', iostat=status )
    if ( status /= 0 ) then
      call io_error ( "While opening L2CF", status, line )
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Unable to open L2CF file " // trim(line) )
    end if
    inunit = l2cf_unit
  end if
  call configuration ( root )
  if ( pcf ) then
    call closeMLSCF
  else
    close ( l2cf_unit )  ! Don't worry about the status
  end if
  if ( root <= 0 ) then
    call output ( &
      'A syntax error occurred -- there is no abstract syntax tree', &
      advance='yes' )
  else
    if ( dump_tree ) call print_subtree ( root, 0 )

    ! Check that supra-syntactic conditions are met, e.g. correct
    ! types for fields of commands, correct command order, etc.
    call check_tree ( root, error, first_section )
    if ( do_dump ) call dump_decl
    if ( toggle(syn) ) then
      call output ( 'Begin abstract syntax tree:', advance='yes' )
      call print_subtree ( root, 0 )
      call output ( 'End abstract syntax tree', advance='yes' )
    end if
    if ( error == 0 .and. first_section /= 0 ) then
      ! Now do the L2 processing.
      call walk_tree_to_do_MLS_L2 ( root, error, first_section )
    end if
  end if
end program MLSL2

! $Log$
! Revision 2.6  2001/02/23 02:38:34  vsnyder
! Open L2CF either by PCF or by Fortran OPEN or expect it on stdin
!
! Revision 2.5  2001/02/22 23:51:00  vsnyder
! Improved usage messages
!
! Revision 2.4  2001/02/22 23:05:12  vsnyder
! Display usage if -h, -H or -? option is present.
!
! Revision 2.3  2000/10/12 00:33:47  vsnyder
! Insert CVS variables and copyright notice
!
