program TEST_SPARSIFY

  use MatrixModule_0, only: Dump, MatrixElement_T, R8, Sparsify

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
