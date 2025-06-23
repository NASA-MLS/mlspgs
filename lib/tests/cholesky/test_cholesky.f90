! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

program TEST_CHOLESKY

! Test the CholeskyFactor subroutine in MatrixModule_0

  use Machine ! At least for HP, and maybe for GETARG
  use MatrixModule_0, only: CholeskyFactor, CreateBlock, Densify, Dump, &
    & M_Banded, M_Full, MatrixElement_T, operator(.TX.), operator(+), &
    & SolveCholesky
  use MLSCommon, only: R8
  implicit NONE

  integer :: BandHeight = 0
  logical :: BAND = .false.
  !          -f,           -i           -n or            -s option specified:
  integer :: DUMP_FAC = 0, DUMP_IN = 0, DUMP_NORMAL = 0, DUMP_SOL = 0
  logical :: EARLY = .false.            ! Print timing early (and at the end)
  integer :: I, IERR, NR, NC, OVERLAP
  character(len=127) :: LINE
  logical :: PROGRESS = .false.
  logical :: RANDOM = .false.
  type(MatrixElement_T) :: RHS          ! RHS of normal equations (random #)
  real(r8), allocatable :: RHSD(:,:)    ! RHS of normal equations (random #)
  real(r8), allocatable :: S(:)         ! Solution, from DCHOL
  logical :: Square = .false.           ! Set by the -q option
  type(MatrixElement_T) :: S1, S2       ! Steps in solution, from SolveCholesky
  real :: T1, T2, T3, T4, T5, T6, T7, T8 ! For timing
  real(r8) :: U = 0.0_r8                ! for DCHOL
  type(MatrixElement_T) :: Z            ! The input
  type(MatrixElement_T) :: ZF           ! Factored normal equations
  type(MatrixElement_T) :: ZTEST        ! ZF .TX. ZF
  type(MatrixElement_T) :: ZTZ          ! Z^T Z
  real(r8), pointer :: ZTZD(:,:)        ! Z^T Z densified

  !---------------------------- RCS Ident Info -------------------------------
  character (len=256) :: Id = &
    & "$Id$"
  !---------------------------------------------------------------------------

  i = 1
  do
    call getarg ( i+hp, line )
    if ( line(1:3) == "-b" ) then
      band = .true.
    else if ( line(1:2) == "-e" ) then
      early = .true.
    else if ( line(1:2) == "-f" ) then
      if ( line(3:) == " " ) then
        i = i + 1
        call getarg ( i+hp, line(3:) )
      end if
      read (line(3:),*) dump_fac
    else if ( line(1:2) == "-i" ) then
      if ( line(3:) == " " ) then
        i = i + 1
        call getarg ( i+hp, line(3:) )
      end if
      read (line(3:),*) dump_in
    else if ( line(1:2) == "-n" ) then
      if ( line(3:) == " " ) then
        i = i + 1
        call getarg ( i+hp, line(3:) )
      end if
      read (line(3:),*) dump_normal
    else if ( line(1:3) == "-p " ) then
      progress = .true.
    else if ( line(1:2) == "-q" ) then
      if ( line(3:) == " " ) then
        i = i + 1
        call getarg ( i+hp, line(3:) )
      end if
      read ( line(3:), * ) nr
      nc = nr
      square = .true.
    else if ( line(1:3) == "-r" ) then
      random = .true.
    else if ( line(1:3) == "-s" ) then
      dump_sol = 1
    else if ( line /= " " ) then
      call getarg ( 0+hp, line )
      print *, 'Usage: ', trim(line), ' [options]'
      print *, ' Options: -b => create a banded matrix'
      print *, '          -e => print timing early (and at the end too)'
      print *, '          -f<number> => dump the factor if number > 0'
      print *, '          -i<number> => dump the input if number > 0'
      print *, '          -n<number> => dump normal equations if number > 0'
      print *, '          -p => print "progress" lines'
      print *, '          -q[ ]<number> => Matrix is square with "number" rows'
      print *, '               Doesn''t work with -b'
      print *, '          -r => create a matrix full of random numbers'
      print *, '          -s => dump the solution'
      print *, '          If number <= 1, only the structure is dumped.'
      stop
    else
      exit
    end if
    i = i + 1
  end do

  if ( band ) then ! Want to fill matrix with random column stripes
    print *, 'Enter nonzeroes per column, number of columns, and number of'
    print *, 'rows that adjacent columns overlap: '
    read (*,*) NR, NC, OVERLAP
    overlap = max(min(overlap,nr),0)
    call createBlock ( z, nr + (nr-overlap)*(nc-1), nc, M_Banded, nr*nc )
    z%r2(0) = 0
    do i = 1, nc
      z%r1(i) = (i-1)*(nr-overlap)+1
      z%r2(i) = i*nr
    end do
    call random_number ( z%values )
  else
    if ( .not. square ) then
      print *, 'Enter the number of rows and columns in the matrix:'
      read (*,*) NR, NC
    end if
    call createBlock ( z, nr, nc, M_Full )
    if ( random ) then ! Fill matrix with random numbers
      call random_number ( z%values )
    else
      print *, "Enter the input matrix, column-by-column: "
      read (*,*) z%values
    end if
  end if

  if ( progress ) call sayTime ( "Forming RHS" )
  allocate ( rhsd(nc,1) )     ! Right-hand side of normal equations
  call random_number ( rhsd )
  call createBlock ( rhs, nc, 1, M_Full )
  rhs%values = rhsd
  allocate ( s(nc) )
  s = rhsd(:,1)

  if ( progress ) call sayTime ( "Forming normal equations" )
  call cpu_time ( t1 )
  ztz = z .tx. z              ! Normal equations matrix
  call cpu_time ( t2 )
  if ( early ) print *, "Normal equations = ", t2 - t1, " seconds for", &
    z%nrows, z%ncols, "; Band Height =", BandHeight
  if ( ztz%kind == M_Full ) then
    BandHeight = ztz%nrows
  else
    BandHeight = maxval(ztz%r2(1:ztz%ncols)-ztz%r2(0:ztz%ncols-1))
  end if
  allocate ( ztzd(nc,nc) )
  if ( progress ) call sayTime ( "Densifying normal equations" )
  call densify ( ztzd, ztz )  ! A copy of ZTZ, because DCHOL will clobber it

  ! Get a solution using the Cholesky decomposition routine from Math77
  if ( progress ) call sayTime ( "Beginning DCHOL" )
  call cpu_time ( t3 )
  call dchol ( ztzd, nc, nc, s, 0.0d0, 1.0d-12, ierr )
  call cpu_time ( t4 )
  if ( early ) print *, "DCHOL =", t4 - t3, " seconds for", shape(ztzd)
  if ( ierr /= 0 ) then
    print *, "DCHOL reports IERR =", ierr
    stop
  end if
  if ( dump_sol > 0 ) then
100 format (a, 1p, 4(1x,g13.7)/5(1x,g13.7))
    print 100, "RHS =         ", rhsd
    print 100, "Solution =    ", s
  end if

  if ( dump_in > 0 ) call dump ( z, "Input", dump_in - 1 )
  if ( dump_normal > 0 ) call dump ( ztz, "Normal equations", dump_normal - 1 )
  if ( progress ) call sayTime ( "Beginning CholeskyFactor" )
  call cpu_time ( t5 )
  call choleskyFactor ( zf, ztz )
  call cpu_time ( t6 )
  if ( early ) print *, "Cholesky factor =", t6 - t5, " seconds."
  if ( dump_fac > 0 ) call dump ( zf, "Factor of input", dump_fac - 1 )
  if ( progress ) call sayTime ( "Beginning SolveCholesky" )
  call solveCholesky ( zf, s1, rhs, .true. )
  if ( dump_sol > 0 ) call dump ( s1, "S1", dump_sol - 1 )
  call solveCholesky ( zf, s2, s1 )
  call cpu_time ( t7 )
  if ( early ) print *, "Back solve =", t7 - t6, " seconds."
  if ( dump_sol > 0 ) call dump ( s2, "S2", dump_sol - 1 )
  if ( progress ) call sayTime ( "Beginning zf .tx. zf" )
  ztest = zf .tx. zf
  call cpu_time ( t8 )
  if ( progress ) call sayTime ( "Finished zf .tx. zf" )
  if ( early ) print *, "Square factors=", t8 - t7, " seconds."
  ztest%values = -ztest%values
  ztest = ztz + ztest ! actually ztz - ztest
  print *, "Z^T Z - G^T G =", MAXVAL( ABS( ztest%values  ) )
  print *, "Error in solution =", MAXVAL( ABS (s2%values(:,1) - s ) )
  print *, "Normal equations = ", t2 - t1, " seconds for", &
    z%nrows, z%ncols, "; Band Height =", BandHeight
  print *, "DCHOL =", t4 - t3, " seconds for", shape(ztzd)
  print *, "Cholesky factor =", t6 - t5, " seconds."
  print *, "Back solve =", t7 - t6, " seconds."
  print *, "Square factors=", t8 - t7, " seconds."
  ! Now do it in place and see if it gets the same answer
! call choleskyFactor ( ztz )
! print *, "inPlace - separate = ", maxval( abs(zf%values-ztz%values) )

contains
  subroutine SayTime ( Prefix )
    character(len=*), intent(in) :: Prefix
    character(len=10) :: Time
    write ( *, '(1x,a)', advance='no' ) Prefix
    write ( *, '(a)', advance='no' ) ' at '
    call date_and_time ( time=time )
    write ( *, '(a)' ) time(1:2) // ':' // time(3:4) // ':' // time(5:)
  end subroutine SayTime
end program TEST_CHOLESKY

! $Log$
! Revision 1.3  2000/11/16 02:41:45  vsnyder
! Corrected data type of constant arguments of dchol.
!
! Revision 1.2  2000/11/15 23:31:39  vsnyder
! Correct use of a constant with unspecified intent dummy argument
!
! Revision 1.1  2000/11/09 17:29:47  pwagner
! First commit
!
! Revision 2.1  2000/10/31 23:32:32  vsnyder
! Initial commit
!
