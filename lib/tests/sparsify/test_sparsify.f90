! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

! $Id$

program TEST_SPARSIFY

  use MatrixModule_0, only: Dump, MatrixElement_T, Sparsify
  use MLSCommon, only: R8

  real(r8), allocatable :: Z(:,:)
  type(MatrixElement_T) :: B

  integer :: NR, NC

  write(0,*) "Enter number of rows and columns in Z:"
  read *, NR, NC
  allocate ( Z(NR,NC) )
  write(0,*) "Enter Z, a column at a time"
  read *, Z
  call dump ( z, name="Input =" )
  call sparsify ( z, b )
  call dump ( b )
end program TEST_SPARSIFY

! $Log$
! Revision 2.3  2000/10/12 20:10:36  vsnyder
! Get R8 from MLSCommon instead of MatrixMultiply_0
!
! Revision 2.2  2000/10/10 23:14:22  vsnyder
! Added copyright, cvs Id and Log
!
