! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

! $Id$

program TEST_MULTIPLY

  use MatrixModule_0, only: Densify, Dump, MatrixElement_T, &
    & operator(.XT.), Sparsify
  use MLSCommon, only: R8

  implicit NONE

  real(r8), dimension(:,:), allocatable :: A, B, Z, ZT
  type(MatrixElement_T) :: AB, BB, ZB
  logical :: DUMP_IT = .true., DUMP_SPARSE = .false., QUIET = .false.
  integer :: I
  character(len=127) :: LINE
  integer :: NDIFF

  integer :: NRA, NCA, NCB

  i = 0
  do
    i = i + 1
    call getarg ( i, line )
    if ( line == " " ) exit
    if ( line == "-s " ) then
      dump_sparse = .true.
    else if ( line == "-n " ) then
      dump_it = .false.
    else if ( line == "-d " ) then
      dump_it = .true.
    else if ( line == "-q " ) then
      quiet = .true.
    end if
  end do
  if ( .not. quiet ) &
  &  write(0,*) "Enter numbers of rows and columns of A, and columns of B: "
  read *, NRA, NCA, NCB
  allocate ( A(NRA,NCA), B(NRA,NCB) )
  if ( .not. quiet ) &
  &  write(0,*) "Enter A, a column at a time: "
  read *, A
  if ( .not. quiet ) &
  &  write(0,*) "Enter B, a column at a time: "
  read *, B
  call sparsify ( a, ab )
  call sparsify ( b, bb )
  if ( dump_it ) then
    if ( dump_sparse ) then
      call dump ( ab, name = "Input A =" )
      call dump ( bb, name = "Input B =" )
    else
      call dump ( a, name = "Input A =" )
      call dump ( b, name = "Input B =" )
    end if
  end if
  zb = ab .xt. bb
  allocate ( z(zb%nrows,zb%ncols), zt(zb%nrows,zb%ncols) )
  call densify ( z, zb )
  zt = matmul(transpose(a),b)
  ndiff = count(z-zt > 0.0)
  if ( dump_it .or. ndiff > 0 ) then
    if ( dump_sparse ) then
      call dump ( zb, name = "Output Z = A^T B =" )
    else
      call dump ( z, name = "Output Z = A^T B =" )
    end if
  end if
  print *, "Kinds =", ab%kind, bb%kind, &
    & "Number of differences to correct answer = ", ndiff
  if ( ndiff > 0 ) then
    call dump ( zt, "True answer =" )
    call dump ( z-zt, "Z - ZT =" )
  end if
end program TEST_MULTIPLY

! $Log$
! Revision 2.2  2000/10/12 20:10:36  vsnyder
! Get R8 from MLSCommon instead of MatrixMultiply_0
!
! Revision 2.1  2000/10/11 19:10:33  vsnyder
! Initial entry
!
