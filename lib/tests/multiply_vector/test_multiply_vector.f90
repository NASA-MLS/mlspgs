program TEST_MULTIPLY_VECTOR

! Test the matrix-vector multiply routine in MatrixModule_0.

  use Machine ! At least for HP, and maybe for GETARG
  use MatrixModule_0, only: CreateBlock, Densify, Dump, M_Banded, &
    & M_Column_Sparse, M_Full, MatrixElement_T, MultiplyMatrixVector, &
    & Sparsify
  use MLSCommon, only: R8

  logical :: BAND = .false.        ! make a Banded matrix
  integer :: BW                    ! Band width of banded matrix
  logical :: CHECK = .true.        ! Check the result
  logical :: COL = .false.         ! make a Column-sparse matrix
  integer :: DUMP_M = 0            ! Dump level for matrix: 0 => nothing,
                                   !  1 => summary 2 => details
  integer :: DUMP_P = 0            ! Dump level for product: 0 => nothing,
                                   !  1 => summary 2 => details
  integer :: DUMP_V = 0            ! Dump level for vector: 0 => nothing,
                                   !  1 => summary 2 => details
  integer :: I, J, N               ! Subscripts and loop inductors
  character(len=127) :: LINE       ! Command line arguments
  type(matrixElement_T) :: M       ! Matrix
  real(r8), pointer :: MD(:,:)     ! Matrix, densified
  integer :: NCM, NRM              ! Size of Matrix
  integer :: NNZ                   ! Number of nonzeroes in M
  integer :: NV                    ! Size of P and V
  real(r8), allocatable :: P(:), V(:)    ! Product and Vector
  real(r8) :: SPARSITY             ! What fraction of M is nonzeroes?
  real :: T1, T2                   ! For timing

  !---------------------------- RCS Ident Info -------------------------------
  character (len=256) :: Id = &
    & "$Id$"
  !---------------------------------------------------------------------------

  i = 1
  do
    call getarg ( i+hp, line )
    if ( line(1:3) == "-b" ) then
      band = .true.
    else if ( line(1:3) == "-c" ) then
      col = .true.
    else if ( line(1:2) == "-m" ) then
      read (line(3:),*) dump_m
    else if ( line(1:2) == "-n " ) then
      check = .false.
    else if ( line(1:2) == "-p" ) then
      read (line(3:),*) dump_p
    else if ( line(1:2) == "-v" ) then
      read (line(3:),*) dump_v
    else if ( line /= " " ) then
      call getarg ( 0+hp, line )
      print *, 'Usage: ', trim(line), ' [options]'
      print *, ' Options: -b => create a banded matrix'
      print *, '          -c => create a column sparse matrix'
      print *, '          -m<number> => dump the matrix if number > 0'
      print *, '          -n => timing only -- don''t check the result'
      print *, '          -p<number> => dump the product if number > 0'
      print *, '          -v<number> => dump the vector if number > 0'
      print *, '          If number <= 1, only the structure is dumped.'
      stop
    else
      exit
    end if
    i = i + 1
  end do

  print *, 'Enter the size of the "vector": '
  read (*,*) nv
  allocate ( p(nv), v(nv) )
  ncm = nv
  nrm = nv

  if ( band ) then       ! Want a banded matrix
    print *, 'Enter the bandwidth (diagonal + superdiagonals; do not count &
      &subdiagonals): '
    read ( *,* ) bw
    bw = min( nrm, max( (nrm+ncm-1)/ncm, bw ) )
    nnz = ncm + ncm*(ncm-1) - (ncm-bw)*(ncm-bw+1)
    call createBlock ( m, nrm, ncm, M_banded, nnz )
    m%r1(1:bw-1) = 1
    do i = bw, ncm
      m%r1(i) = i-bw+1
    end do ! i = 1, ncm
    do i = 1, ncm
      m%r2(i) = m%r2(i-1) + min( 2*bw-1, i+bw-1, nrm+1-m%r1(i) )
    end do ! i = 1, ncm
    call random_number ( m%values )
    if ( check ) then
      allocate ( md(nrm,ncm) )
      call densify ( md, m )
    end if
  else if ( col ) then   ! Want a column-sparse matrix
    print *, 'Enter the faction of M that is to be nonzeroes, in (0,1]: '
    read ( *,* ) sparsity
    if ( sparsity <= 0.0 .or. sparsity > 1.0 ) then
      print *, 'Sparsity is not in (0,1]'
      stop
    end if
    n = nrm * ncm
    nnz = max(INT(n * sparsity), ncm)
    call createBlock ( m, nrm, ncm, M_Column_Sparse, nnz )
    call random_number ( m%values )
    m%r2 = m%values(:,1) * n
    j = 1
    do while ( j /= 0 )
      call isort ( m%r2, 1, nnz )
      j = 0 ! number of duplicates
      do i = 2, nnz
        if ( m%r2(i) == m%r2(i-1) ) then
          call random_number ( m%values(i,1) )
          m%r2(i) = m%values(i,1) * n
          j = j + 1
        end if
      end do
    end do
    call random_number ( m%values )
    j = 0
    do i = 1, nnz
      n = m%r2(i) / nrm + 1
      do while ( j < n )
        j = j + 1
        m%r1(j) = m%r1(j-1)
      end do
      m%r1(j) = i
      m%r2(i) = mod(m%r2(i), nrm) + 1
    end do ! i = 1, nnz
      do while ( j < ncm )
        j = j + 1
        m%r1(j) = m%r1(j-1)
      end do
    if ( check ) then
      allocate ( md(nrm,ncm) )
      call densify ( md, m )
    end if
  else                   ! Want a full matrix
    call createBlock ( m, nrm, ncm, M_Full )
    call random_number ( m%values )
    md => m%values
  end if
  if ( dump_m > 0 ) call dump ( m, "Matrix", dump_m > 1 )

  call random_number ( v )
  if ( dump_v == 1 ) then
    print *, 'Vector has', nv, ' elements'
  else if ( dump_v > 1 ) then
    call dump ( v, "Vector" )
  end if

  call cpu_time ( t1 )
  call multiplyMatrixVector ( m, v, p )
  call cpu_time ( t2 )
  print *, 'Time for multiply =', t2-t1

  if ( dump_p  == 1 ) then
    print *, 'Product has', nv, ' elements'
  else if ( dump_v > 1 ) then
    call dump ( p, "Product" )
  end if

  if ( check ) then
    if ( dump_m > 2 ) call dump ( md, "Dense matrix" )
    print *, 'Error =', &
    & maxval( abs( matmul(transpose(md),reshape(v,(/size(v)/))) - &
    &              reshape(p,(/size(p)/)) ) )
  end if
end program TEST_MULTIPLY_VECTOR

! $Log$
! Revision 2.3  2000/11/09 23:23:57  vsnyder
! Moved to tests/multiply
!
! Revision 2.2  2000/11/09 01:22:06  vsnyder
! Removed "private" in a parameter statement -- can't work in a main program
!
! Revision 2.1  2000/10/31 23:35:52  vsnyder
! Initial commit
!
