! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Sparse_Eta_m

! Get interpolation coefficients in a Sparse_t matrix

  use MLSKinds, only: RP

  implicit NONE

  private
  public :: Sparse_Eta_1D, Sparse_Eta_nD

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Sparse_Eta_1D ( Basis, Grid, Eta, What, Row1, Rown, Create, &
                           & Sorted, Resize )

    ! Compute Eta for linear interpolation from 1D Basis to 1D Grid.

    use Sparse_m, only: Sparse_t, Create_Sparse, Add_Element, Resize_Sparse

    real(rp), intent(in) :: Basis(:)
    real(rp), intent(in) :: Grid(:)
    type(sparse_t), intent(inout) :: Eta    ! Might not be created here
    integer, intent(in), optional :: What   ! String index
    integer, intent(in), optional :: Row1, RowN ! Part of Grid to use, default all
    logical, intent(in), optional :: Create ! Create Eta, default true
    logical, intent(in), optional :: Sorted ! "Grid is sorted" -- default true
    logical, intent(in), optional :: Resize ! Re-size Eta%E to Eta%NE

    logical, parameter :: CheckSorted = .false. ! Check whether Grid is sorted
    real(rp) :: Del_Basis
    integer :: I, J, K
    logical :: MyCreate
    integer :: MyRow1, MyRowN
    logical :: MySorted
    integer :: N_Basis, N_Grid
    real(rp) :: V(2) ! Values of coefficient

    mySorted = .true.
    if ( present(sorted) ) mySorted = sorted

    n_basis = size(basis)
    n_grid = size(grid)
 
    myCreate = .true.
    if ( present(create) ) myCreate = create

    if ( myCreate ) call create_sparse ( eta, size(grid), size(basis), &
                                       & 2*size(grid), what=what )

    myRow1 = 1
    myRowN = n_grid
    if ( present(row1) ) myRow1 = min(row1,n_grid)
    if ( present(rowN) ) myRowN = min(rowN,n_grid)

    if ( mySorted ) then
      if ( grid(myRow1) > grid(myRowN) ) then ! Grid is in the opposite order
        i = myRow1
        myRow1 = myRowN
        myRowN = i
      end if

      if ( checkSorted ) then
        block
          use Dump_0, only: Dump
          use MLSMessageModule, only: MLSMessage, MLSMSG_Crash
          use Output_m, only: Output
          j = sign(1,myRowN-myRow1)
          do i = myRow1, myRowN-j, j
            if ( grid(i) > grid(i+j) ) then
              call dump ( grid(min(myRow1,myRowN):max(myRow1,myRowN)), 'Grid', &
                & lBound=min(myRow1,myRowN) )
              call output ( i, before='Grid(' )
              call output ( grid(i), before=' ) = ' )
              call output ( i+j, before=' < Grid(' )
              call output ( grid(i+j), before=' ) = ', advance='yes' )
              call MLSMessage ( MLSMSG_Crash, moduleName, &
                & "You lied!  Grid is claimed to be sorted, but it isn't" )
            end if
          end do
        end block
      end if

      if ( myRow1 <= myRowN ) then ! Process grid in increasing order
        ! Coefficients below Basis(1) are all 1.0
        do i = myRow1, myRowN
          if ( grid(i) > basis(1) ) exit
          call add_element ( eta, 1.0_rp, i, 1 )
        end do
        ! Coefficients between Basis(1) and Basis(n_basis) are "hat" functions
        do j = 2, n_basis
          del_basis = 1.0_rp / ( basis(j) - basis(j-1) )
          do while ( i <= n_grid )
            if ( grid(i) > basis(j) ) exit
            if ( basis(j-1) < grid(i) ) then
              v = [ (basis(j)-grid(i)) * del_basis, &
                &   (grid(i)-basis(j-1)) * del_basis ]
              do k = 1, 2
                if ( v(k) /= 0.0 ) call add_element ( eta, v(k), i, j+k-2 )
              end do
            end if
            i = i + 1
          end do
        end do
        ! Coefficients above Basis(n_basis) are all 1.0
        do i = myRowN, i, -1
          if ( basis(n_basis) >= grid(i) ) exit
          call add_element ( eta, 1.0_rp, i, n_basis )
        end do
      else ! Process grid in decreasing order
        ! Coefficients below Basis(1) are all 1.0
        do i = myRow1, myRowN, -1
          if ( grid(i) > basis(1) ) exit
          call add_element ( eta, 1.0_rp, i, 1 )
        end do
        ! Coefficients between Basis(1) and Basis(n_basis) are "hat" functions
        do j = 2, n_basis
          del_basis = 1.0_rp / ( basis(j) - basis(j-1) )
          do while ( i > myRowN )
            if ( grid(i) > basis(j) ) exit
            if ( basis(j-1) < grid(i) ) then
              v = [ (basis(j)-grid(i)) * del_basis, &
                &   (grid(i)-basis(j-1)) * del_basis ]
              do k = 1, 2
                if ( v(k) /= 0.0 ) call add_element ( eta, v(k), i, j+k-2 )
              end do
            end if
            i = i - 1
          end do
        end do
        ! Coefficients above Basis(n_basis) are all 1.0
        do i = myRowN, i
          if ( basis(n_basis) >= grid(i) ) exit
          call add_element ( eta, 1.0_rp, i, n_basis )
        end do
      end if
    else ! Grid needs to be sorted
      block
        use Sort_m, only: SortP
        integer :: P(n_grid)
        call sortp ( grid, 1, n_grid, p ) ! grid(p(:)) are now sorted
        ! Coefficients below Basis(1) are all 1.0
        do i = 1, n_grid
          if ( grid(p(i)) > basis(1) ) exit
          call add_element ( eta, 1.0_rp, p(i), 1 )
        end do
        ! Coefficients between basis(1) and Basis(n_basis) are "hat" functions
        do j = 2, n_basis
          del_basis = 1.0_rp / ( basis(j) - basis(j-1) )
          do while ( i <= n_grid )
            if ( grid(p(i)) > basis(j) ) exit
            if ( basis(j-1) < grid(p(i)) ) then
              v = [ (basis(j)-grid(p(i))) * del_basis, &
                &   (grid(p(i))-basis(j-1)) * del_basis ]
              do k = 1, 2
                if ( v(k) /= 0.0 ) then
              if ( v(k) /= 0.0 ) call add_element ( eta, v(k), p(i), j+k-2 )
                end if
              end do
            end if
            i = i + 1
          end do
        end do
        ! Coefficients above Basis(n_basis) are all 1.0
        do i = n_grid, i, -1
          if ( basis(n_basis) >= grid(p(i)) ) exit
        call add_element ( eta, 1.0_rp, p(i), n_basis )
        end do
      end block
    end if

    if ( present(resize) ) then
      if ( resize ) call resize_sparse ( eta )
    end if

  end subroutine Sparse_Eta_1D

  subroutine Sparse_Eta_nD ( L, R, P, What, Resize )

    ! Compute P as the product of L and R, e.g. Freq and Zeta x Phi, for
    ! linear interpolation from an nD Basis (e.g. Freq x Zeta x Phi, to 1D
    ! Grid, e.g. line-of-sight).

    use MLSMessageModule, only: MLSMessage, MLSMsg_Error
    use Sparse_m, only: Sparse_t, Create_Sparse, Add_Element, Resize_Sparse

    type(sparse_t), intent(inout) :: L      ! Assume it's been created
    type(sparse_t), intent(inout) :: R      ! Assume it's been created
    type(sparse_t), intent(out) :: P        ! Created here
    integer, intent(in), optional :: What   ! String index
    logical, intent(in), optional :: Resize ! Re-size Eta%E to Eta%NE

    integer :: I     ! Row index
    integer :: J, K  ! Column indices
    integer :: M     ! Column dimension of L
    real(rp) :: V    ! L(i,j) * R(i,k)

    if ( size(l%rows) /= size(r%rows) ) &
      & call MLSMessage ( MLSMsg_Error, moduleName, &
        & "L and R in Sparse_Eta_nD have different row dimension extents" )

    call create_sparse ( p, size(l%rows), size(l%cols)*size(r%cols), &
                       & 4*size(l%rows), &
                       & ubnd=[ l%ubnd, r%ubnd ], &
                       & lbnd=[ l%lbnd, r%lbnd ], what=what )

    m = size(l%cols)

    ! Compute the column elements of P in each row as the products of
    ! all elements of L with all elements of R in that row.
    do i = 1, size(l%rows)
      if ( l%rows(i) == 0 .or. r%rows(i) == 0 ) cycle ! Don't include empty row
      k = r%e(r%rows(i))%nr ! First element of R in row I
      do
        j = l%e(l%rows(i))%nr       ! First element of L in row I
        do
          v = l%e(j)%v * r%e(k)%v
          if ( v /= 0.0 ) &
            & call add_element ( p, v, i, m*(r%e(k)%c-1) + l%e(j)%c )
          j = l%e(j)%nr             ! Next element of L in row I
          if ( j == l%e(l%rows(i))%nr ) exit ! Back to the first element?
        end do
        k = r%e(k)%nr               ! Next element of R in row I
        if ( k == r%e(r%rows(i))%nr ) exit ! Back to the first element?
      end do
    end do ! rows

    if ( present(resize) ) then
      if ( resize ) call resize_sparse ( p )
    end if

  end subroutine Sparse_Eta_nD

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Sparse_Eta_m

! $Log$
! Revision 2.1  2017/11/01 18:52:11  vsnyder
! Initial commit
!
