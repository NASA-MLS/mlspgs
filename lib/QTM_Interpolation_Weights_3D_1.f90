! Copyright 2015, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module QTM_Interpolation_Weights_3D_1
!=============================================================================

  implicit NONE
  private

  ! Get interpolation weights from QTM_Interpolation_Weights_3D_m and put
  ! them into a column-sparse matrix.  This is in a different module from
  ! QTM_Interpolation_Weights_3D_m so that the latter does not depend upon
  ! either Allocate_Deallocate, MatrixModule_0, or MLSMessageModule.

  public :: Fill_Matrix_From_Weights, QTM_Interpolation_Weights

  interface QTM_Interpolation_Weights
    module procedure &
      & QTM_Interpolation_Weights_Geo_3D_Incoherent_Line_Matrix, &
      & QTM_Interpolation_Weights_Geo_3D_Line_Matrix
  end interface QTM_Interpolation_Weights

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine Fill_Matrix_From_Weights ( Weights, N_Weights, Matrix )
    ! Create a column-sparse matrix block from Weights and N_Weights

    use MatrixModule_0, only: CreateBlock, M_Column_Sparse, MatrixElement_t
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use QTM_Interpolation_Weights_3D_m, only: Weight_t

    type(weight_t), intent(in) :: Weights (:,:)
    integer, intent(in) :: N_Weights(:)
    type(matrixElement_t), intent(inout) :: Matrix

    integer :: I
    integer :: NCols                     ! from Matrix
    logical :: New                       ! Re-create the matrix block
    integer :: NNZ                       ! Number of nonzeroes in Weights
    integer :: NRows                     ! from Matrix
    integer :: NVal                      ! Number of values so far

    nCols = matrix%nCols
    nRows = matrix%nRows
    if ( size(n_weights) > nCols ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "Matrix has fewer columns than line has points" )

    nnz = sum(n_weights)
    ! Re-create the matrix block as column sparse with the appropriate number
    ! of nonzeroes.
    new = matrix%kind /= m_column_sparse
    if ( associated(matrix%r2) ) new = size(matrix%r2) < nnz
    if ( associated(matrix%value1) ) new = size(matrix%value1) < nnz
    if ( new ) then
      call createBlock ( matrix, nRows, nCols, m_column_sparse, nnz, &
        & forWhom="Fill_Matrix_From_Weights" )
    else
      matrix%r2(nnz+1:) = 0
    end if
    nVal = 0

    ! Fill the nonzeroes locators (matrix%r2 and matrix%r1) and the values
    do i = 1, nCols
      matrix%r2(nVal+1:nVal+n_weights(i)) = weights(1:n_weights(i),i)%which
      matrix%value1(nVal+1:nVal+n_weights(i)) = weights(1:n_weights(i),i)%weight
      nVal = nVal + n_weights(i)
      matrix%r1(i) = nVal
    end do

  end subroutine Fill_Matrix_From_Weights

  subroutine QTM_Interpolation_Weights_Geo_3D_Incoherent_Line_Matrix ( QTM_Tree, &
                         & Heights, Line, Points, Matrix )

    ! Construct interpolation weights and positions for a point in a 3D grid
    ! for which the horizontal grid is QTM.  The 3D grid is assumed to be
    ! stacked but not necessarily coherent.  The positions are array element
    ! order positions, which can be converted to 2D subscripts (QTM serial
    ! number, height) by the Subscripts function in the Array_Stuff module.

    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: ECR_t, RG
    use MatrixModule_0, only: MatrixElement_T
    use QTM_Interpolation_Weights_3D_m, only: QTM_Interpolation_Weights, S_QTM_T, &
      & Weight_t

    type(QTM_tree_t), intent(in) :: QTM_Tree
    real(rg), intent(in) :: Heights(:,:) ! Extents are (heights, QTM_Tree%N_In)
                                         ! Assume Heights and QTM_Tree came from
                                         ! the same Geolocation_t structure
                                         ! (see Geolocation_m).
    type(ECR_t), intent(inout) :: Line(2) ! Vector to a point on the line,
                                         ! vector along the line.  Line(2) is
                                         ! made an unit vector here, just in
                                         ! case it wasn't on input.
    type(S_QTM_t), intent(in) :: Points(:) ! Points along Line
    type(matrixElement_t), intent(inout), optional :: Matrix ! Weights

    integer :: N_Weights(size(points))
    type(weight_t) :: Weights(6,size(points))

    if ( size(heights,2) == 1 ) then ! Coherent heights
      call QTM_Interpolation_Weights ( QTM_Tree, Heights(:,1), &
                                     & Line, Points, Matrix )
      return
    end if

    call QTM_interpolation_weights ( QTM_Tree, &
                         & heights, line, points, weights, n_weights )

    call fill_matrix_from_weights ( weights, n_weights, matrix )

  end subroutine QTM_Interpolation_Weights_Geo_3D_Incoherent_Line_Matrix

  subroutine QTM_Interpolation_Weights_Geo_3D_Line_Matrix ( QTM_Tree, &
                         & Heights, Line, Points, Matrix )

    ! Construct interpolation weights and positions for a point in a 3D grid
    ! for which the horizontal grid is QTM.  The 3D grid is assumed to be
    ! stacked but not necessarily coherent.  The positions are array element
    ! order positions, which can be converted to 2D subscripts (QTM serial
    ! number, height) by the Subscripts func

    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: ECR_t, RG
    use MatrixModule_0, only: MatrixElement_T
    use QTM_Interpolation_Weights_3D_m, only: QTM_Interpolation_Weights, S_QTM_T, &
      & Weight_t
  
    type(QTM_tree_t), intent(in) :: QTM_Tree
    real(rg), intent(in) :: Heights(:)   ! Extents are (heights, QTM_Tree%N_In)
                                         ! Assume Heights and QTM_Tree came from
                                         ! the same Geolocation_t structure
                                         ! (see Geolocation_m).
    type(ECR_t), intent(inout) :: Line(2) ! Vector to a point on the line,
                                         ! vector along the line.  Line(2) is
                                         ! made an unit vector here, just in
                                         ! case it wasn't on input.
    type(S_QTM_t), intent(in) :: Points(:) ! Points along Line
    type(matrixElement_t), intent(inout), optional :: Matrix ! Weights

    integer :: N_Weights(size(points))
    type(weight_t) :: Weights(6,size(points))

    call QTM_interpolation_weights ( QTM_Tree, &
                         & heights, line, points, weights, n_weights )

    call fill_matrix_from_weights ( weights, n_weights, matrix )

  end subroutine QTM_Interpolation_Weights_Geo_3D_Line_Matrix

!=============================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module QTM_Interpolation_Weights_3D_1

! $Log$
! Revision 2.1  2016/04/16 02:04:59  vsnyder
! Initial commit
!
