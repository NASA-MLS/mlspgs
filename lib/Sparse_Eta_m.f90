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

! Compute interpolation coefficients in a Sparse_Eta_t matrix

  use MLSKinds, only: RP
  use Sparse_m, only: Sparse_t

  implicit NONE

  private
  public :: Sparse_Eta_t

  type, extends(sparse_t) :: Sparse_Eta_t
    ! No new components
  contains
    procedure, pass(eta) :: Eta_1D => Sparse_Eta_1D
    procedure, pass(p) :: Eta_nD => Sparse_Eta_nD
    generic :: Eta => Eta_1D, Eta_nD
  end type Sparse_Eta_t

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Sparse_Eta_1D ( Basis, Grid, Eta, What, Row1, Rown, Create, &
                           & Sorted, Resize )

    ! Compute Eta for linear interpolation from 1D Basis to 1D Grid.
    use Allocate_Deallocate, only: Test_Allocate

    real(rp), intent(in) :: Basis(:)
    real(rp), intent(in) :: Grid(:)
    class(sparse_eta_t), intent(inout) :: Eta ! Might not be created here
    integer, intent(in), optional :: What     ! String index for dumps
    integer, intent(in), optional :: Row1, RowN ! Part of Grid to use,
                                              ! default all
    logical, intent(in), optional :: Create   ! Force Eta to be created, default
                                              ! false.  Eta is created anyway if
                                              ! rows or cols are not allocated
                                              ! or the wrong sizes.
    logical, intent(in), optional :: Sorted   ! "Grid is sorted" -- default true
    logical, intent(in), optional :: Resize   ! Re-size Eta%E to Eta%NE --
                                              ! default false

    logical, parameter :: CheckSorted = .false. ! Check whether Grid is sorted

    real(rp) :: Del_Basis
    integer :: I, J
    character(127) :: Msg
    logical :: MyCreate
    integer :: MyRow1, MyRowN
    logical :: MySorted
    integer :: N_Basis, N_Grid
    integer :: PR   ! Previous row, to avoid creating duplicates
    integer :: Stat
    real(rp) :: V   ! Value of coefficient

    mySorted = .true.
    if ( present(sorted) ) mySorted = sorted

    n_basis = size(basis)
    n_grid = size(grid)
 
    myCreate = .false.
    if ( present(create) ) myCreate = create

    if ( allocated(eta%rows) .and. allocated(eta%cols) ) then
      if ( size(eta%rows) < n_grid .or. size(eta%cols) /= n_basis ) &
        & myCreate = .true.
    else
      myCreate = .true.
    end if

    if ( myCreate ) then
      call eta%create ( n_grid, n_basis, 2*size(grid), what=what )
    else if ( .not. allocated(eta%e) ) then
      allocate ( eta%e(2*n_grid), stat=stat, errmsg=msg )
      call test_allocate ( stat, moduleName, "Sparse%E", ermsg=msg )
    end if

    eta%nRows = n_grid
    if ( present(row1) ) eta%nRows = row1
    if ( present(rowN) ) eta%nRows = rowN

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
          call eta%add_element ( 1.0_rp, i, 1 )
          pr = i
        end do
        ! Coefficients between Basis(1) and Basis(n_basis) are "hat" functions
        do j = 2, n_basis
          del_basis = 1.0_rp / ( basis(j) - basis(j-1) )
          do while ( i <= myRowN )
            if ( grid(i) > basis(j) ) exit
            if ( basis(j-1) < grid(i) .and. i /= pr ) then
              v = ( basis(j)-grid(i) ) * del_basis
              if ( v /= 0 ) call eta%add_element ( v, i, j-1 )
              v = (grid(i)-basis(j-1)) * del_basis
              if ( v /= 0 ) call eta%add_element ( v, i, j )
            end if
            pr = i
            i = i + 1
          end do
        end do
        ! Coefficients above Basis(n_basis) are all 1.0
        do i = myRowN, i, -1
          ! I think ">" ought to be ">=" to avoid duplicates, but it somehow
          ! sometimes neglects to do one, so we check for duplicates.
          if ( basis(n_basis) > grid(i) ) exit
          if ( i /= pr ) call eta%add_element ( 1.0_rp, i, n_basis )
        end do
      else ! Process grid in decreasing order
        ! Coefficients below Basis(1) are all 1.0
        do i = myRow1, myRowN, -1
          if ( grid(i) > basis(1) ) exit
          call eta%add_element ( 1.0_rp, i, 1 )
          pr = i
        end do
        ! Coefficients between Basis(1) and Basis(n_basis) are "hat" functions
        do j = 2, n_basis
          del_basis = 1.0_rp / ( basis(j) - basis(j-1) )
          do while ( i > myRowN )
            if ( grid(i) > basis(j) ) exit
            if ( basis(j-1) < grid(i) .and. i /= pr ) then
              v = ( basis(j)-grid(i) ) * del_basis
              if ( v /= 0 ) call eta%add_element ( v, i, j-1 )
              v = (grid(i)-basis(j-1)) * del_basis
              if ( v /= 0 ) call eta%add_element ( v, i, j )
            end if
            pr = i
            i = i - 1
          end do
        end do
        ! Coefficients above Basis(n_basis) are all 1.0
        do i = myRowN, i
          ! I think ">" ought to be ">=" to avoid duplicates, but it somehow
          ! sometimes neglects to do one, so we check for duplicates.
          if ( basis(n_basis) > grid(i) ) exit
          if ( i /= pr ) call eta%add_element ( 1.0_rp, i, n_basis )
          pr = i
        end do
      end if
    else ! Grid needs to be sorted
      block
        use Sort_m, only: SortP
        integer :: P(n_grid)
        call sortp ( grid, 1, n_grid, p ) ! grid(p(:)) are now sorted
        ! Coefficients below Basis(1) are all 1.0
        do i = myRow1, myRowN
          if ( grid(p(i)) > basis(1) ) exit
          call eta%add_element ( 1.0_rp, p(i), 1 )
          pr = i
        end do
        ! Coefficients between basis(1) and Basis(n_basis) are "hat" functions
        do j = 2, n_basis
          del_basis = 1.0_rp / ( basis(j) - basis(j-1) )
          do while ( i <= myRowN )
            if ( grid(p(i)) > basis(j) ) exit
            if ( basis(j-1) < grid(p(i)) .and. i /= pr ) then
              v = ( basis(j)-grid(i) ) * del_basis
              if ( v /= 0 ) call eta%add_element ( v, i, j-1 )
              v = (grid(i)-basis(j-1)) * del_basis
              if ( v /= 0 ) call eta%add_element ( v, i, j )
            end if
            i = i + 1
          end do
        end do
        ! Coefficients above Basis(n_basis) are all 1.0
        do i = myRowN, i, -1
          ! I think ">" ought to be ">=" to avoid duplicates, but it somehow
          ! sometimes neglects to do one, so we check for duplicates.
          if ( basis(n_basis) > grid(p(i)) ) exit
          if ( i /= pr ) call eta%add_element ( 1.0_rp, p(i), n_basis )
          pr = i
        end do
      end block
    end if

    if ( present(resize) ) then
      if ( resize ) call eta%resize
    end if

  end subroutine Sparse_Eta_1D

  subroutine Sparse_Eta_nD ( L, R, P, What, Resize, One_Row_OK )

    ! Compute P as the product of L and R, e.g. Freq and Zeta x Phi, for
    ! linear interpolation from an nD Basis (e.g. Freq x Zeta x Phi, to 1D
    ! Grid, e.g. line-of-sight).

    use MLSMessageModule, only: MLSMessage, MLSMsg_Error
    use MoreMessage, only: MLSMessage

    class(sparse_eta_t), intent(in) :: L    ! Assume it's been created
    class(sparse_eta_t), intent(in) :: R    ! Assume it's been created
    class(sparse_eta_t), intent(out) :: P   ! Created here
    integer, intent(in), optional :: What   ! String index
    logical, intent(in), optional :: Resize ! Re-size Eta%E to Eta%NE
    logical, intent(in), optional :: One_Row_OK ! OK if L has only one row,
                                            ! e.g., for Frequency.  Apply that
                                            ! row of L to every row of R.

    integer :: I       ! Row index
    integer :: J, K    ! Column indices
    integer :: FL, FR  ! Index of final column in row of left, right factors
    integer :: M       ! Column dimension of L
    logical :: OK      ! One_Row_OK, default false
    integer :: LR, LRI ! L row index, increment
    real(rp) :: V      ! L(i,j) * R(i,k)

    ok = .false.
    if ( present(one_row_ok) ) ok = one_row_ok
    lri = 1            ! Increment row for L
    if ( ok .and. ( l%nRows == 1 .or. size(l%rows) == 1 ) ) then
      lri = 0          ! Don't increment row for L
    else if ( l%nRows /= 0 .and. r%nRows /= 0 ) then
      if ( l%nRows /= r%nRows ) &
      & call MLSMessage ( MLSMsg_Error, moduleName, &
        & "L(%d) and R(%d) in Sparse_Eta_nD have different numbers of useful rows", &
        & [ l%nRows, r%nRows ] )
      p%nRows = l%nRows
    else if ( size(l%rows) /= size(r%rows) ) then
      call MLSMessage ( MLSMsg_Error, moduleName, &
        & "L(%d) and R(%d) in Sparse_Eta_nD have different row dimension extents", &
        & [ size(l%rows), size(r%rows) ] )
    end if

    call p%create ( size(r%rows), size(l%cols)*size(r%cols), 4*size(r%rows), &
                  & ubnd=[ l%ubnd, r%ubnd ], &
                  & lbnd=[ l%lbnd, r%lbnd ], what=what )

    m = size(l%cols)

    ! Compute the column elements of P in each row as the outer product of
    ! all elements of L in that row (or its only row) with all elements of R
    ! in that row.
    lr = 1
    do i = 1, r%nRows
      fr = r%rows(i)
      if ( l%rows(lr) /= 0 .and. fr /= 0 ) then ! Don't include empty row
        k = r%e(fr)%nr ! First element of R in row I
        if ( k /= 0 ) then
          do
            fl = l%rows(lr)
            j = l%e(fl)%nr       ! First element of L in row I
            if ( j /= 0 ) then
              do
                v = l%e(j)%v * r%e(k)%v
                if ( v /= 0.0 ) &
                  & call p%add_element ( v, i, m*(r%e(k)%c-1) + l%e(j)%c )
                j = l%e(j)%nr             ! Next element of L in row I
                if ( j == l%e(fl)%nr ) exit ! Back to the first element?
              end do
            end if
            k = r%e(k)%nr               ! Next element of R in row I
            if ( k == r%e(fr)%nr ) exit ! Back to the first element?
          end do
        end if
      end if
      lr = lr + lri
    end do ! rows
    p%nRows = r%nRows

    if ( present(resize) ) then
      if ( resize ) call p%resize ( p%ne )
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
! Revision 2.3  2018/03/07 00:21:19  vsnyder
! Don't make names of procedures that are type-bound public.  Add Eta generic
! for Sparse_Eta_1D and Sparse_Eta_nD.  Allow the left factor in Sparse_Eta_nD
! to have only one row.
!
! Revision 2.2  2017/11/29 00:33:39  vsnyder
! Add Sparse_Eta_t as extension of Sparse_t.  Make Eta_1D and Eta_nD type-
! bound to Sparse_Eta_t.  Better criteria to create.  Ad hoc hand waving
! to avoid creating duplicates.  Correct bugs in Eta_nD.
!
! Revision 2.1  2017/11/01 18:52:11  vsnyder
! Initial commit
!
