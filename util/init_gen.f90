program INIT_GEN

  ! Generate the parameter statements and references to add_ident
  ! for init_tables_module and others similar to it.

  ! An input file specifies the names.
  ! You can put one or two names on a line, separated by spaces.
  ! The first name is the parameter name WITHOUT ITS PREFIX.
  ! The second name is the entity's text -- e.g. a field name.  If
  ! it's not present, the first name is used for the entity's text.

  ! Blank lines don't count.

  ! The command line control everything else.  See the PRINT statements
  ! that explain the usage.

  use MACHINE

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = "$RCSfile$"
!---------------------------------------------------------------------------

  character(len=30) :: CapName               ! Capitalized? name
  integer, parameter :: Decl_Wid = len("  integer, parameter :: ")
  logical :: DoCap = .false.                 ! "Capitalize declaration"
  character(len=40) :: FirstName = '1', LastName = ' '
  integer :: I, J                            ! subscripts, loop inductors
  character(len=127) :: InFile
  character(len=30) :: Index_Name            ! Array for output from add_ident
  integer :: Index_Wid                       ! Width of Index_Name
  integer :: IOSTAT
  character(len=127) :: LINE
  integer :: Margin = 0                      ! Internal margin
  integer :: MaxWid                          ! of a name
  character(len=30), dimension(:), allocatable :: Names ! for the symbol table
  integer :: Nspaces
  integer :: NumNames
  character(len=127) :: Out_Add, Out_Parm    ! File names
  character(len=40), dimension(:), allocatable :: P_Names   ! Parameter names
  !                                            that go with Names
  character(len=10) :: Pfx = ' '             ! for declarations

  i = hp
  do
    i = i + 1
    call getarg ( i, line )
    if ( line(1:3) == '-c ' ) then
      doCap = .true.
    else if ( line(1:2) == '-f' ) then
      if ( line(3:3) == ' ' ) then
        i = i + 1
        call getarg ( i, line(3:) )
      end if
      firstName = line(3:)
    else if ( line(1:2) == '-l' ) then
      if ( line(3:3) == ' ' ) then
        i = i + 1
        call getarg ( i, line(3:) )
      end if
      lastName = line(3:)
    else if ( line(1:2) == '-m' ) then
      if ( line(3:3) == ' ' ) then
        i = i + 1
        call getarg ( i, line(3:) )
      end if
      read ( line(3:), *, iostat=iostat ) margin
      if ( iostat /= 0 ) then
        call io_error ( 'While reading "m" option', iostat )
        margin = 0
      end if
    else if ( line(1:2) == '-p' ) then
      if ( line(3:3) == ' ' ) then
        i = i + 1
        call getarg ( i, line(3:) )
      end if
      pfx = line(3:)
    else if ( line(1:1) == '-' ) then
      call getarg ( hp, line )
      print *, 'Usage: ', trim(line), ' [options] input_file out_parm out_add index_name'
      print *, ' Options:  -c: Capitalize entire name in declarations'
      print *, '           -f[ ]first_name -- The expression from which the'
      print *, '             first name gets its value.  Default 1.'
      print *, '           -l[ ]last_name -- A name made equal to the last name.'
      print *, '             Not emitted if not present.'
      print *, '           -m[ ]margin -- The "middle" margin, which is otherwise'
      print *, '             computed from the width of the widest name.'
      print *, '           -p[ ]prefix -- Prefix for variable names.'
      print *, '           -<anything else>: This output.'
      print *, ' input_file: The file of names'
      print *, ' out_parm:   The file to store the generated parameter declarations.'
      print *, ' out_add:    The file to store the generated references to add_ident.'
      print *, ' index_name: The name of the array in which to store the outputs'
      print *, '             of add_ident.'
      stop
    else
  exit
    end if
  end do

  call getarg ( i, inFile )
  if ( inFile == ' ' ) then
    print *, 'No input file name given.'
    stop
  end if
  call getarg ( i+1, out_parm )
  if ( out_parm == ' ' ) then
    print *, 'No parameter declarations file name given.'
    stop
  end if
  call getarg ( i+2, out_add )
  if ( out_add == ' ' ) then
    print *, 'No file name given for references to add_ident.'
    stop
  end if
  call getarg ( i+3, index_name )
  if ( index_name == ' ' ) then
    print *, 'No name given for the array to store references to add_ident.'
    stop
  end if
  index_Wid = len_trim(index_name)

  open ( 10, file=inFile, form='formatted', status='old', iostat=iostat )
  if ( iostat /= 0 ) then
    call io_error ( 'Opening input file', iostat, inFile )
    stop
  end if

  open ( 11, file=out_add, form='formatted', iostat=iostat )
  if ( iostat /= 0 ) then
    call io_error ( 'Opening file for results of add_ident', iostat, inFile )
    stop
  end if

  open ( 12, file=out_parm, form='formatted', iostat=iostat )
  if ( iostat /= 0 ) then
    call io_error ( 'Opening file for parameter declarations', iostat, inFile )
    stop
  end if

  ! Count the names
  numNames = 0
  do
    read ( 10, '(a)', iostat=iostat ) line
    if ( iostat < 0 ) exit
    if ( iostat > 0 ) then
      call io_error ( 'Reading input file', iostat, inFile )
      stop
    end if
    if ( line == ' ' ) cycle
    line = adjustl(line)
    if ( line(1:1) /= '#' ) numNames = numNames + 1
  end do
  rewind ( 10 )

  ! Allocate space for the names
  allocate ( names(numNames), stat=iostat )
  allocate ( P_names(numNames), stat=iostat )
  if ( iostat /= 0 ) then
    call io_error ( 'Allocating table to store the names', iostat )
    stop
  end if

  ! Read the names.  Make sure the file size hasn't changed.  Calculate
  ! the width of the widest name.
  maxWid = 0
  i = 0
  do
    read ( 10, '(a)', iostat=iostat ) line
    if ( iostat < 0 ) then
      if ( i == numNames ) exit
      print *, 'The file of names changed size!'
      stop
    end if
    if ( iostat > 0 ) then
      call io_error ( 'Reading input file', iostat, inFile )
      stop
    end if
    if ( line == ' ' ) cycle
    line = adjustl(line)
    if ( line(1:1) == '#' ) cycle
    i = i + 1
    if ( i > numNames ) then
      print *, 'The file of names changed size!'
      stop
    end if
    j = index(line,' ')
    p_names(i) = line(:j-1)
    line = adjustl(line(j+1:))
    if ( line == ' ' ) then
      names(i) = p_names(i)
    else
      names(i) = line
    end if
    p_names(i) = trim(pfx) // p_names(i)
    maxWid = max(maxWid, len_trim(p_names(i)))
  end do
  line = ' '
  close ( 10 )

  if ( lastName /= ' ' ) maxWid = max(maxWid, len_trim(lastName))

  ! Output the parameter declarations.

  line = ' ' ! used to output spaces

100 format ( "  integer, parameter :: ",a," = ",a,a)
  do i = 1, numNames
    capName = capitalize(p_names(i))
    capName(1:1) = cap(capName(1:1))
    if ( margin /= 0 ) then
      nspaces = margin - decl_wid - len_trim(p_names(i)) - 4
    else
      nspaces = maxwid - len_trim(capName)
    end if
    write ( 12, 100 ) trim(capName), line(:nspaces), &
                    & trim(firstName)
    firstName = trim(p_names(i)) // ' + 1'
  end do
  if ( lastName /= ' ' ) then
    lastName = capitalize ( lastName )
    lastName(1:1) = cap(lastName(1:1))
    if ( margin /= 0 ) then
      nspaces = margin - decl_wid - len_trim(lastName) - 4
    else
      nspaces = maxwid - len_trim(lastName)
    end if
    write ( 12, 100 ) trim(lastName), line(:nspaces), &
                    & trim(p_names(numNames))
  end if
  close ( 12 )

  ! Output the references to "add_ident"
  do i = 1, numNames
    if ( margin /= 0 ) then
      nspaces = margin - index_Wid - 10 - len_trim(names(i))
    else
      nspaces = maxwid - len_trim(p_names(i))
    end if
    write ( 11, '(4x,a,"(",a,") = ", a, "add_ident ( ''",a,"'' )")' ) &
      & index_name(:index_Wid), trim(p_names(i)), &
      & line(:nspaces), trim(names(i))
  end do
  close ( 11 )

contains

  character function Cap ( TheChar )
    character, intent(in) :: TheChar
    cap = theChar
    if ( cap >= 'a' .and. cap <= 'z' ) &
      & cap = achar(iachar(cap) + iachar('A') - iachar('a'))
  end function Cap

  function Capitalize ( TheName ) result ( CapName )
    character(len=*), intent(in) :: TheName
    character(len(theName)) :: CapName
    integer :: J
    capName = theName
    if ( doCap ) then
      do j = 1, len_trim(theName)
        capName(j:j) = cap(capName(j:j))
      end do
    end if
  end function Capitalize

end program INIT_GEN

! $Log$
