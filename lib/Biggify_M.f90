! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Biggify_M

  use Allocate_Deallocate, only: Allocate_Test
  use MatrixModule_0, only: Densify
  use MatrixModule_1, only: Matrix_T, RC_Info
  use MLSCommon, only: R8
  use VectorsModule, only: Vector_T

  interface Biggify
    module procedure Biggify_Matrix, Biggify_Vector
  end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

  ! ---------------------------------------------  Biggify_Matrix  -----
  subroutine Biggify_Matrix ( Matrix, Array )

  ! Convert Matrix from a Matrix_T to Array, a Fortran array.
  ! The primary purpose of this routine is for debugging, but it
  ! could also be useful for output.

    type(matrix_T), intent(in) :: Matrix
    real(r8), dimension(:,:), pointer :: Array

    integer :: I, J           ! Subscripts, loop inductors
    integer :: II, JJ         ! Subscripts of block(i,j) in Array
    integer :: M, N           ! Dimensions of Array
    integer :: MM, NN         ! Dimensions of a block of Matrix

    m = sum(matrix%row%nelts)
    n = sum(matrix%col%nelts)
    call allocate_test ( array, m, n, 'Array', moduleName )
    jj = 1
    do j = 1, matrix%col%nb
      nn = matrix%block(1,j)%ncols
      ii = 1
      do i = 1, matrix%row%nb
        mm = matrix%block(i,1)%nrows
        call densify ( Array(ii:ii+mm-1,jj:jj+nn-1), &
          & matrix%block(i,j) )
        ii = ii + mm
      end do ! i
      jj = jj + nn
    end do ! j

  end subroutine Biggify_Matrix

  ! ---------------------------------------------  Biggify_Vector  -----
  subroutine Biggify_Vector ( RC, Vector, Array )

  ! Convert Vector from a Vector_T to Array, a Fortran array, in the
  ! order specified by RC.
  ! The primary purpose of this routine is for debugging, but it
  ! could also be useful for output.

    type(rc_Info), intent(in) :: RC
    type(vector_T), intent(in) :: Vector
    real(r8), dimension(:), pointer :: Array

    integer :: I              ! Subscript in RC%inst, RC%quant
    integer :: M              ! Subscript in Array of beginning of quantity
    integer :: N              ! Number of elements in a quantity
    integer :: NB             ! Number of blocks specified by RC, minus 1
    !                           if RC%Extra

    call allocate_test ( array, vector%template%totalElements, 'Array', &
      & moduleName )
    m = 1
    nb = rc%nb
    if ( rc%extra ) nb = nb - 1
    do i = 1, nb
      n = rc%nelts(i)
      array(m:m+n-1) = vector%quantities(rc%quant(i))%values(:,rc%inst(i))
      m = m + n
    end do
    if ( rc%extra ) array(m) = 1.0 ! So it at least has a value
  end subroutine Biggify_Vector

end module Biggify_M

! $Log$
! Revision 2.1  2001/06/01 01:14:27  vsnyder
! For debugging only
!
