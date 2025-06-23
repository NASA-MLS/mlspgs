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
    procedure, pass(eta) :: Eta_0D => Sparse_Eta_0D
    procedure, pass(eta) :: Eta_1D => Sparse_Eta_1D
    procedure, pass(p) :: Eta_nD => Sparse_Eta_nD
    procedure, pass(eta) :: Eta_QTM => Sparse_Eta_QTM
    generic :: Eta => Eta_0D, Eta_1D, Eta_nD, Eta_QTM
  end type Sparse_Eta_t

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Sparse_Eta_0D ( Basis, Grid, Eta, What, Create, Resize )

    ! Compute Eta for linear interpolation from 1D Basis to a single point.
    use Allocate_Deallocate, only: Test_Allocate
    use Pure_Hunt_m, only: PureHunt

    real(rp), intent(in) :: Basis(:)
    real(rp), intent(in) :: Grid
    class(sparse_eta_t), intent(inout) :: Eta ! Might not be created here
    integer, intent(in), optional :: What     ! String index for dumps
    logical, intent(in), optional :: Create   ! Force Eta to be created, default
                                              ! false.  Eta is created anyway if
                                              ! rows or cols are not allocated
                                              ! or the wrong sizes.
    logical, intent(in), optional :: Resize   ! Re-size Eta%E to Eta%NE --
                                              ! default false

    real(rp) :: Del_Basis
    integer :: JLO, JHI
    character(127) :: Msg
    logical :: MyCreate
    integer :: N_Basis
    integer :: Stat
    real(rp) :: V   ! Value of coefficient

    n_basis = size(basis)
 
    myCreate = .false.
    if ( present(create) ) myCreate = create

    if ( allocated(eta%rows) .and. allocated(eta%cols) ) then
      if ( size(eta%rows) < 1 .or. size(eta%cols) /= n_basis ) &
        & myCreate = .true.
    else
      myCreate = .true.
    end if

    if ( myCreate ) then
      call eta%create ( 1, n_basis, 2, what=what )
    else if ( .not. allocated(eta%e) ) then
      allocate ( eta%e(2), stat=stat, errmsg=msg )
      call test_allocate ( stat, moduleName, "Sparse%E", ermsg=msg )
    end if

    eta%nRows = 1

    if ( grid <= basis(1) ) then
      ! Coefficient below Basis(1) is 1.0
      call eta%add_element ( 1.0_rp, 1, 1 )
    else if ( grid >= basis(n_basis) ) then
      ! Coefficient above Basis(n_basis) is 1.0
      call eta%add_element ( 1.0_rp, 1, n_basis )
    else
      call purehunt ( grid, basis, n_basis, jlo, jhi )
      ! Assume Basis is increasing.  Basis(JLO) <= Grid <= Basis(JHI) here.
      jhi = jlo + 1 ! In case PureHunt returned JHI == JLO
      ! "Hat" function between Basis(JLO) and Basis(JHI)
      del_basis = 1.0_rp / ( basis(jhi) - basis(jlo) )
      v = ( basis(jhi)-grid ) * del_basis
      if ( v /= 0 ) call eta%add_element ( v, 1, jlo )
      v = (grid-basis(jlo)) * del_basis
      if ( v /= 0 ) call eta%add_element ( v, 1, jhi )
    end if

    if ( present(resize) ) then
      if ( resize ) call eta%resize
    end if

  end subroutine Sparse_Eta_0D

  subroutine Sparse_Eta_1D ( Basis, Grid, Eta, What, Row1, Rown, Create, &
                           & Empty, Sorted, Resize, Basis1, BasisN )

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
    logical, intent(in), optional :: Empty    ! Make Eta empty before starting,
                                              ! default false.  If absent or
                                              ! false, add new ones.  Not quite
                                              ! as traumatic as Create=.true.
    logical, intent(in), optional :: Sorted   ! "Grid is sorted" -- default true
    logical, intent(in), optional :: Resize   ! Re-size Eta%E to Eta%NE --
                                              ! default false
    integer, intent(in), optional :: Basis1, BasisN ! Part of Basis to use,
                                              ! default all

    logical, parameter :: CheckSorted = .false. ! Check whether Grid is sorted

    real(rp) :: Del_Basis
    integer :: I, J
    character(127) :: Msg
    integer :: MyBasis1, MyBasisN
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

    if ( present(empty) ) then
      if ( empty ) call eta%empty
    end if

    if ( allocated(eta%rows) .and. allocated(eta%cols) ) then
      if ( size(eta%rows) < n_grid .or. size(eta%cols) /= n_basis ) &
        & myCreate = .true.
    else
      myCreate = .true.
    end if

    if ( myCreate ) then
      call eta%create ( n_grid, n_basis, 2*n_grid, what=what )
    else if ( .not. allocated(eta%e) ) then
      allocate ( eta%e(2*n_grid), stat=stat, errmsg=msg )
      call test_allocate ( stat, moduleName, "Sparse%E", ermsg=msg )
      eta%ne = 0
    end if

    eta%nRows = n_grid
    if ( present(row1) ) eta%nRows = row1
    if ( present(rowN) ) eta%nRows = rowN

    myRow1 = 1
    myRowN = n_grid
    if ( present(row1) ) myRow1 = min(row1,n_grid)
    if ( present(rowN) ) myRowN = min(rowN,n_grid)

    myBasis1 = 1
    myBasisN = n_basis
    if ( present(basis1) ) myBasis1 = basis1
    if ( present(basisN) ) myBasisN = basisN

    pr = 0 ! Make sure first row is not believed to be a duplicate

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
        ! Coefficients below Basis(myBasis1) are all 1.0
        do i = myRow1, myRowN
          if ( grid(i) > basis(myBasis1) ) exit
          call eta%add_element ( 1.0_rp, i, 1 )
          pr = i
        end do
        ! Coefficients between Basis(1) and Basis(n_basis) are "hat" functions
        do j = myBasis1+1, myBasisN
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
        ! Coefficients above Basis(myBasisN) are all 1.0
        do i = myRowN, i, -1
          ! I think ">" ought to be ">=" to avoid duplicates, but it somehow
          ! sometimes neglects to do one, so we check for duplicates.
          if ( basis(myBasisN) > grid(i) ) exit
          if ( i /= pr ) call eta%add_element ( 1.0_rp, i, n_basis )
        end do
      else ! Process grid in decreasing order
        ! Coefficients below Basis(1) are all 1.0
        do i = myRow1, myRowN, -1
          if ( grid(i) > basis(myBasis1) ) exit
          call eta%add_element ( 1.0_rp, i, 1 )
          pr = i
        end do
        ! Coefficients between Basis(1) and Basis(n_basis) are "hat" functions
        do j = myBasis1+1, myBasisN
          del_basis = 1.0_rp / ( basis(j) - basis(j-1) )
          do while ( i >= myRowN )
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
          if ( basis(myBasisN) > grid(i) ) exit
          if ( i /= pr ) call eta%add_element ( 1.0_rp, i, n_basis )
          pr = i
        end do
      end if
    else ! Grid needs to be sorted
      block
        use Sort_m, only: SortP
        integer :: P(n_grid)
        integer :: Pi ! p(i)
        call sortp ( grid, 1, n_grid, p ) ! grid(p(:)) are now sorted
        ! Coefficients below Basis(1) are all 1.0
        do i = myRow1, myRowN
          pi = p(i)
          if ( grid(pi) > basis(myBasis1) ) exit
          call eta%add_element ( 1.0_rp, pi, 1 )
          pr = i
        end do
        ! Coefficients between basis(1) and Basis(n_basis) are "hat" functions
        do j = myBasis1+1, myBasisN
          del_basis = 1.0_rp / ( basis(j) - basis(j-1) )
          do while ( i <= myRowN )
            pi = p(i)
            if ( grid(pi) > basis(j) ) exit
            if ( basis(j-1) < grid(pi) .and. i /= pr ) then
              v = ( basis(j)-grid(pi) ) * del_basis
              if ( v /= 0 ) call eta%add_element ( v, pi, j-1 )
              v = (grid(pi)-basis(j-1)) * del_basis
              if ( v /= 0 ) call eta%add_element ( v, pi, j )
            end if
            pr = i
            i = i + 1
          end do
        end do
        ! Coefficients above Basis(n_basis) are all 1.0
        do i = myRowN, i, -1
          pi = p(i)
          ! I think ">" ought to be ">=" to avoid duplicates, but it somehow
          ! sometimes neglects to do one, so we check for duplicates.
          if ( basis(myBasisN) > grid(pi) ) exit
          if ( i /= pr ) call eta%add_element ( 1.0_rp, pi, n_basis )
          pr = i
        end do
      end block
    end if

    if ( present(resize) ) then
      if ( resize ) call eta%resize
    end if

  end subroutine Sparse_Eta_1D

  subroutine Sparse_Eta_QTM ( Basis, Lon, GeodLat, Eta, What, Row1, Rown, &
                            & Create, Empty, Resize )

    ! Compute interpolation coefficients from a QTM Basis to a 2D Grid
    ! represented by longitude and latitude.  The latitude can be geodetic or
    ! geocentric, which will get the wrong answer if Basis and Grid differ.

    use Allocate_Deallocate, only: Test_Allocate
    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: H_Geod
    use QTM_m, only: Stack_t
    use QuantityTemplates, only: RT
    use Triangle_Interpolate_m, only: Triangle_Interpolate

    type (QTM_Tree_t), intent(in) :: Basis
    real (rt), intent(in) :: Lon(:)     ! Degrees
    real (rt), intent(in) :: GeodLat(:) ! Degrees, same size as Lon
    class(sparse_eta_t), intent(inout) :: Eta ! Might not be created here
    integer, intent(in), optional :: What     ! String index for dumps
    integer, intent(in), optional :: Row1, RowN ! Part of Grid to use,
                                              ! default all
    logical, intent(in), optional :: Create   ! Force Eta to be created, default
                                              ! false.  Eta is created anyway if
                                              ! rows or cols are not allocated
                                              ! or the wrong sizes.
    logical, intent(in), optional :: Empty    ! Make Eta empty before starting,
                                              ! default false.  If absent or
                                              ! false, add new ones.  Not quite
                                              ! as traumatic as Create=.true.
    logical, intent(in), optional :: Resize   ! Re-size Eta%E to Eta%NE --
                                              ! default false

    integer :: Facet       ! Facet index in Basis
    type (h_geod) :: Geod  ! One Lon and GeodLat
    integer :: I           ! Index in Lon and GeodLat
    integer :: J           ! Indices of vertices within Facet
    character(127) :: Msg
    logical :: MyCreate
    integer :: MyRow1, MyRowN
    integer :: N_Basis, N_Grid
    type(stack_t) :: Stack ! To make QTM searches faster
    integer :: Stat
    real(rp) :: W(3)       ! Barycentric interpolation weights

    n_basis = basis%n_in   ! Vertices within or adjacent to the polygon
    n_grid = size(lon)

    myCreate = .false.
    if ( present(create) ) myCreate = create

    if ( present(empty) ) then
      if ( empty ) call eta%empty
    end if

    if ( allocated(eta%rows) .and. allocated(eta%cols) ) then
      if ( size(eta%rows) < n_grid .or. size(eta%cols) /= n_basis ) &
        & myCreate = .true.
    else
      myCreate = .true.
    end if

    if ( myCreate ) then
      call eta%create ( n_grid, n_basis, 2*n_grid, what=what )
    else if ( .not. allocated(eta%e) ) then
      allocate ( eta%e(2*n_grid), stat=stat, errmsg=msg )
      call test_allocate ( stat, moduleName, "Sparse%E", ermsg=msg )
    end if

    do i = 1, n_grid
      geod%lon%d = lon(i)   ! h_geod components have different kind, so the
      geod%lat = geodLat(i) !   h_geod() constructor cannot be used.
      facet = basis%find_facet ( geod, stack )
      if ( facet /= 0 ) then ! No interpolation coefficient outside the QTM
        call triangle_interpolate ( basis%Q(facet)%geo%lon%d, &
                                  & basis%Q(facet)%geo%lat, &
                                  & lon(i), geodLat(i), w )
        do j = 1, 3
          call eta%add_element ( w(j), i, basis%Q(facet)%ser(j) )
        end do
      end if
    end do

  end subroutine Sparse_Eta_QTM

  subroutine Sparse_Eta_nD ( L, R, P, What, Resize, One_Row_OK, Flags )

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
    logical, intent(in), optional :: Flags(:)   ! Don't produce elements of P
                                            ! in columns for which Flags is false

    integer :: Col     ! Column index to put into P%E(j)%C
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
            j = l%e(fl)%nr       ! First element of L in row LR
            if ( j /= 0 ) then
              do
                col = m*(r%e(k)%c-1) + l%e(j)%c
                v = 0
                if ( present(flags) ) then
                  if ( flags(col) ) v = l%e(j)%v * r%e(k)%v
                else
                  v = l%e(j)%v * r%e(k)%v
                end if
                if ( v /= 0.0 ) call p%add_element ( v, i, col )
                j = l%e(j)%nr             ! Next element of L in row LR
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
      if ( resize ) call p%resize
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
! Revision 2.14  2020/05/07 21:27:15  vsnyder
! Add Basis1, BasisN to control range of applicability
!
! Revision 2.13  2019/05/07 20:06:53  vsnyder
! Set NE=0 if Eta%E is allocated here when .not. myCreate
!
! Revision 2.12  2018/12/04 02:42:11  vsnyder
! Add Sparse_Eta_QTM
!
! Revision 2.11  2018/10/26 02:52:59  vsnyder
! Add Empty optional argument to Sparse_Eta_1D
!
! Revision 2.10  2018/10/23 20:44:38  vsnyder
! Make sure PR is defined everywhere in Sparse_Eta_1D
!
! Revision 2.9  2018/10/11 00:32:26  vsnyder
! Give PR a value in the loop wherein the grid and basis are within ranges.
! Otherwise, I might be compared to an unitialized value, and skip a row.
!
! Revision 2.8  2018/09/05 21:01:30  vsnyder
! Correct off-by-one error that only occurs if basis is in descending order
! and the grid does not extend beyond the bases.
!
! Revision 2.7  2018/08/20 23:40:03  vsnyder
! Correct error in the case that the grid needs to be sorted
!
! Revision 2.6  2018/05/24 03:21:43  vsnyder
! Add flags to allow saying 'don\'t bother with this column'
!
! Revision 2.5  2018/05/17 02:16:48  vsnyder
! Add Sparse_Eta_0D
!
! Revision 2.4  2018/04/11 19:30:29  vsnyder
! Call p%resize without argument. Repair some comments
!
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
