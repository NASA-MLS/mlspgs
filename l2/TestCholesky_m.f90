module TestCholesky_M

  use Biggify_M, only: Biggify
  use MatrixModule_0, only: DenseCholesky, Densify
  use MatrixModule_1, only: GetMatrixElement, Matrix_T, Matrix_Cholesky_T, &
    & Matrix_SPD_T
  use MLSCommon, only: R8

contains

  subroutine TestCholesky ( NormalEquations, SparseCholesky )

    type(matrix_SPD_T), intent(in) :: NormalEquations
    type(matrix_Cholesky_T), intent(in), optional :: SparseCholesky

    double precision, external :: DDOT
    real(r8), pointer, save :: BigCholesky(:,:), BigNeq(:,:)
    real(r8) :: U

    integer :: I, J           ! Subscripts, loop inductors
    integer :: IERR           ! Error signal from DCHOL
    integer :: II, JJ         ! Subscripts of block(i,j) in bigCholesky
    integer :: M, N           ! Dimensions of bigCholesky
    integer :: MM, NN         ! Dimensions of a block of NormalEquations

    if ( .not. present(sparseCholesky) ) then
      m = sum(normalEquations%m%row%nelts)
      n = sum(normalEquations%m%col%nelts)
      allocate ( bigCholesky(m,n) )
      call biggify ( normalEquations%m, bigNeq )
      u = 0.0
      call denseCholesky ( bigCholesky, bigNeq, ierr )
      if ( ierr /= 0 ) then
        print *, 'DCHOL reports error at ', ierr
        stop
      end if
      return
    end if
    do j = 1, n
      do i = 1, j
        if ( abs(getMatrixElement(sparseCholesky%m,i,j)-bigCholesky(i,j)) > 0.0_r8 ) then
          if ( abs( (getMatrixElement(sparseCholesky%m,i,j)-bigCholesky(i,j)) / &
                    (getMatrixElement(sparseCholesky%m,i,j)+bigCholesky(i,j)) ) > 1.0e-6_r8 ) then
            print *, 'Sparse and Dense differ at (', i, ',', j, ')'
            ierr = ierr + 1
          end if
        end if
      end do
    end do
    if ( ierr /= 0 ) then
      print *, 'Sparse and Dense differ at ', ierr, ' places'
      stop
    end if
    deallocate ( bigCholesky, bigNeq, stat=ierr )
  end subroutine TestCholesky 

end module TestCholesky_M
