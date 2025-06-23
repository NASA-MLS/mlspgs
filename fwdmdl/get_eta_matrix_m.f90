! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Get_Eta_Matrix_m

  use MLSKinds, only: RP, IP

  implicit NONE

  private
  public :: Eta_D_T, Eta_S_T
  public :: Clean_Out_Nonzeros, Dump, Dump_Eta_Column_Sparse, Dump_Eta_Transpose
  public :: Eta_Func, Eta_Func_1d, Eta_Func_2d
  public :: Get_Column_Sparsity
  public :: Get_Eta_Column_Sparse, Get_Eta_Column_Sparse_fl_nz, Get_Eta_Sparse
  public :: Get_Eta_Sparse_1d, Get_Eta_Sparse_1d_fl, Get_Eta_Sparse_1d_nz
  public :: Get_Eta_Sparse_2d, Get_Eta_Sparse_2d_fl, Get_Eta_Sparse_2d_fl_nz
  public :: Get_Eta_Sparse_2d_nz
  public :: Get_Eta_Stru, Get_Eta_Stru_D_D, Get_Eta_Stru_D_S, Get_Eta_Stru_S_S
  public :: Get_Eta_1d_Hunt, Get_Eta_ZP
  public :: Interpolate_Stru
  public :: Interpolate_Stru_1D_D, Interpolate_Stru_1D_S
  public :: Interpolate_Stru_2D_D, Interpolate_Stru_2D_S
  public :: Multiply_Eta_Column_Sparse
  public :: Select_NZ_List

  interface Dump
    module procedure Dump_Eta_Column_Sparse, Dump_Eta_Transpose
  end interface

  interface Eta_Func
    module procedure Eta_Func_1d, Eta_Func_2d
  end interface Eta_Func

  interface Get_Eta_Sparse
    module procedure Get_Eta_Column_Sparse, Get_Eta_Column_Sparse_fl_nz
    module procedure Get_Eta_Sparse_1d, Get_Eta_Sparse_1d_fl, Get_Eta_Sparse_1d_nz
    module procedure Get_Eta_Sparse_2d, Get_Eta_Sparse_2d_fl, Get_Eta_Sparse_2d_nz
    module procedure Get_Eta_Sparse_2d_fl_nz
  end interface

  interface Get_Eta_Stru
    module procedure Get_Eta_Stru_D_D, Get_Eta_Stru_D_S, Get_Eta_Stru_S_S
  end interface

  interface Interpolate_Stru
    module procedure Interpolate_Stru_1D_D, Interpolate_Stru_1D_S
    module procedure Interpolate_Stru_2D_D, Interpolate_Stru_2D_S
  end interface

  type :: Eta_S_T
    integer :: L
    real :: Eta
  end type Eta_S_T

  type :: Eta_D_T
    integer :: L
    double precision :: Eta
  end type Eta_D_T

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains

! -------------------------------------------  Clean_Out_Nonzeros  -----
  subroutine Clean_Out_Nonzeros ( Eta_Array, Do_Calc, NZ, NNZ )
    ! Clean out old nonzeros in Eta_Array and old true values in Do_Calc
    use MLSKinds, only: RP
    real(rp), intent(inout) :: Eta_Array(:,:)
    logical, intent(inout), optional :: Do_Calc(:,:) ! Same shape as Eta_Array
    integer, intent(inout), optional :: NZ(:,:)      ! Row subscripts of nonzeros
    integer, intent(inout), optional :: NNZ(:)       ! How many nonzeros in a column
    integer :: I
    if ( present(nz) ) then
      do i = 1, size(eta_array,2)
        eta_array(nz(:nnz(i),i),i) = 0
        if ( present(do_calc) ) do_calc(nz(:nnz(i),i),i) = .false.
        nnz(i) = 0
      end do
    else
      eta_array = 0
      if ( present(do_calc) ) do_calc=.false.
    end if
  end subroutine Clean_Out_Nonzeros

! ---------------------------------------  Dump_Eta_Column_Sparse  -----
  subroutine Dump_Eta_Column_Sparse ( Eta, Nz, NNz, Grids_f, Name, Short, &
    & Show_Z, ZP_only, Small, Width, WhichSps )
    use Array_Stuff, only: Subscripts
    use Dump_0, only: Dump
    use intrinsic, only: Lit_Indices
    use Load_SPS_Data_m, only: Grids_t
    use Output_m, only: Blanks, NewLine, Output
    use String_Table, only: Display_String
    real(rp), intent(in) :: Eta(:,:)
    integer, intent(in) :: Nz(:,:) ! Nonzeroes in each column
    integer, intent(in) :: NNz(:)  ! Number of nonzeroes in each column
    type (grids_t), intent(in) :: Grids_f ! to compute subscripts
    character(len=*), intent(in), optional :: Name
    logical, intent(in), optional :: Short   ! "Don't print zero columns", default false
    logical, intent(in), optional :: Show_Z  ! "Print zeroes", default false
    logical, intent(in), optional :: ZP_Only ! No Freq coordinate for Eta cols
    real(rp), intent(in), optional :: Small  ! Don't print smaller values
    integer, intent(in), optional :: Width   ! How many nonzeros per line
    integer, intent(in), optional :: WhichSps ! Which species to dump
    integer :: Dims(4), I, J, J1, J2, L, MySps, MyWidth
    integer :: Places, Sps, Subs(4)
    logical :: MyShort, MyShow_Z, MyZP
    real(rp) :: MySmall
    logical :: Saw_NZ

    myShort = .false.
    if ( present(short) ) myShort = short
    myShow_Z = .false.
    if ( present(show_Z) ) myShow_Z = show_Z
    myZP = .false.
    if ( present(ZP_Only) ) myZP = ZP_Only
    mySmall = 0
    if ( present(small) ) mySmall = small
    myWidth = 5
    if ( present(width) ) myWidth = width
    mySps = -1
    if ( present(whichSps) ) mySps = whichSps

    j2 = 0
    do sps = 1, ubound(grids_f%l_v,1)
      j1 = j2
      if ( myZP ) then
        j2 = grids_f%l_zp(sps)
        dims = [ size(eta,1), &
               & grids_f%l_z(sps) - grids_f%l_z(sps-1), &
               & grids_f%l_p(sps) - grids_f%l_p(sps-1), 1 ]
      else
        j2 = grids_f%l_v(sps)
        dims = [ size(eta,1), &
               & grids_f%l_f(sps) - grids_f%l_f(sps-1), &
               & grids_f%l_z(sps) - grids_f%l_z(sps-1), &
               & grids_f%l_p(sps) - grids_f%l_p(sps-1) ]
      end if
      if ( mySps > 0 .and. sps /= mySps ) cycle
      if ( present(name) ) then
        call output ( name )
        call blanks ( 1 )
      end if
      if ( grids_f%mol(sps) /= 0 ) then
        call display_string ( lit_indices(grids_f%mol(sps)) )
        call blanks ( 1 )
      end if
      if ( grids_f%l_f(sps) - grids_f%l_f(sps-1) > 1 ) &
        call output ( 'frequency dependent ' )
      if ( grids_f%qty(sps) > 0 ) then
        call display_string ( lit_indices(grids_f%qty(sps)), &
          & before='quantity: ' )
        call blanks ( 1 )
      end if
      call showList ( dims, dims, '')
      dims = [dims(2:), 1]
      call output ( ' state-vector-index ( path-index eta )*', advance='yes' )
      if ( all(nnz(j1+1:j2) == 0) ) then
        call output ( j1+1, before='Columns ' )
        call output ( j2, before=':' )
        call output ( ' are zero', advance='yes' )
        cycle
      end if
      do j = j1+1, j2 ! Corresponds to elements of grids_f%l_v...
        if ( myShort .and. nnz(j) == 0 ) cycle
        subs = subscripts ( j-j1, dims )
        if ( nnz(j) == size(eta,1) ) then
          call output ( 'All rows in column ' )
          call showList ( dims, subs, ' are nonzero')
          call dump ( eta(nz(:nnz(j),j),j), name='Values' )
        else
          places = 5
          places = sum(int(log10(real(subs)))) + count(dims/=1) + 3
          l = 0
          saw_nz = .false.
          do i = 1, nnz(j)
            if ( .not. myShow_Z .and. abs(eta(nz(i,j),j)) <= mySmall ) cycle
            if ( abs(eta(nz(i,j),j)) > mySmall .or. myShow_Z ) then
              if ( .not. saw_nz ) call showList ( dims, subs, '')
              saw_nz = .true.
              l = l + 1
              if ( l > myWidth ) then
                l = 1
                call newLine
                call blanks ( places )
              end if
              call output ( nz(i,j), format='(i4)' )
              call output ( eta(nz(i,j),j), format='(1pg14.6)' )
            end if
          end do
          if ( saw_nz ) call newLine
        end if
      end do
    end do
  contains
    subroutine ShowList ( Dims, Subs, After )
    use Output_m, only: Output
    ! Print "(" subs where dims/=1 ")"
    integer, intent(in) :: Dims(:), Subs(:)
    character(*), intent(in) :: After
    integer :: S
    do s = 1, size(subs)-1
      if ( dims(s) /= 1 ) exit
    end do
    call output ( subs(s), before="(" )
    do while ( s < size(subs) )
      s = s + 1
      if ( dims(s) /= 1 ) call output ( subs(s), before="," )
    end do
    call output ( ")" )
    end subroutine ShowList
  end subroutine Dump_Eta_Column_Sparse

! -------------------------------------------  Dump_Eta_Transpose  -----
  subroutine Dump_Eta_Transpose ( Eta, NPF, Grids_f, Name, N, Width, ZP_Only, &
                                & DoBig, P_only )
    use intrinsic, only: Lit_Indices
    use Load_SPS_Data_m, only: Grids_t
    use Output_m, only: Blanks, NewLine, Output
    use String_Table, only: Display_String
    real(rp), intent(in), contiguous, target :: Eta(:,:)
    integer, intent(in) :: NPF               ! Fine path length = how much
                                             ! of first dimension of Eta to use
    type (grids_t), intent(in) :: Grids_f    ! to compute subscripts
    character(len=*), intent(in), optional :: Name
    integer, intent(in), optional :: N       ! Which species to dump
    integer, intent(in), optional :: Width   ! How many per line, default 5
    logical, intent(in), optional :: ZP_Only ! Eta has no frequency spread
    logical, intent(in), optional :: DoBig   ! Show subs for entire state vector
    logical, intent(in), optional :: P_only  ! Second dimension of Eta is P,
                                             ! not Z, overrides ZP_Only
    integer :: Dims(4,ubound(grids_f%l_z,1))
    real(rp), pointer :: Eta3(:,:,:)
    logical :: FrqDep(ubound(grids_f%l_v,1))
    integer :: I, J, J2(0:ubound(grids_f%l_z,1)), K, L, N1, N2, NC, Sps, W
    integer :: BigSubs(3) ! In entire state vector
    integer :: Subs(3)    ! Only within one species
    equivalence ( subs(1), j ), ( subs(2), k ), (subs(3), l )
    logical :: MyBig, MyP, MyZP, SawNZ
    integer :: MyWidth

    myWidth = 5
    if ( present(width) ) myWidth = width
    myZP = .false.
    if ( present(ZP_Only) ) myZP = ZP_Only
    myBig = .false.
    if ( present(doBig) ) myBig = .true.
    myP = .false.
    if ( present(p_only) ) myP = p_only
    n1 = 1
    n2 = ubound(grids_f%l_z,1)
    if ( myP ) n2 = ubound(grids_f%l_p,1)
    if ( present(n) ) then
      n1 = min(n,n2)
      n2 = n1
    end if
    j2(0) = 0
    do sps = 1, n2
      frqDep(sps) = grids_f%l_f(sps) - grids_f%l_f(sps-1) > 1
      if ( myP ) then
        dims(:,sps) = [ npf, &
                      & grids_f%l_p(sps) - grids_f%l_p(sps-1), 1, 1 ]
        j2(sps) = j2(sps-1) + dims(2,sps)
      else if ( myZP .or. .not. frqDep(sps) ) then
        dims(:,sps) = [ npf, &
                      & grids_f%l_z(sps) - grids_f%l_z(sps-1), &
                      & grids_f%l_p(sps) - grids_f%l_p(sps-1), 1 ]
        j2(sps) = j2(sps-1) + dims(2,sps) * dims(3,sps)
      else
        dims(:,sps) = [ npf, &
                      & grids_f%l_f(sps) - grids_f%l_f(sps-1), &
                      & grids_f%l_z(sps) - grids_f%l_z(sps-1), &
                      & grids_f%l_p(sps) - grids_f%l_p(sps-1) ]
        j2(sps) = j2(sps-1) + dims(2,sps) * dims(3,sps) * dims(4,sps)
      end if
    end do
    do sps = n1, n2
      nc = 3 ! Upper bound of useful part of dims(sps,:)
      if ( frqDep(sps) ) nc = 4
      if ( myP ) then
        nc = 2
      else if ( myZP ) then
        nc = 3
      end if
      if ( present(name) ) then
        call output ( name )
        call blanks ( 1 )
      end if
      if ( grids_f%mol(sps) /= 0 ) then
        call display_string ( lit_indices(grids_f%mol(sps)) )
        call blanks ( 1 )
      end if
      if ( frqDep(sps) ) call output ( 'frequency dependent ' )
      if ( grids_f%qty(sps) > 0 ) then
        call display_string ( lit_indices(grids_f%qty(sps)), &
          & before='quantity: ' )
        call blanks ( 1 )
      end if
      call showList ( dims(:nc,sps) )
      call output ( ' state-vector-index ( path-index eta )*', advance='yes' )
      do i = 1, npf ! path length
        sawNZ = .false.
        w = 0
        eta3(1:dims(2,sps),1:dims(3,sps),1:dims(4,sps)) => eta(i,j2(sps-1)+1:j2(sps))
        do l = 1, size(eta3,3)
          do k = 1, size(eta3,2)
            do j = 1, size(eta3,1)
              if ( eta3(j,k,l) /= 0 ) then
                if ( .not. sawNZ ) call output ( i, places=4, after='#' )
                sawNZ = .true.
                if ( w >= myWidth ) then
                  call newLine
                  call blanks ( 5 )
                  w = 0
                end if
                w = w + 1
                if ( myBig ) then ! Show Grids_f%v subscripts as if it were 4D
                  if ( myZP ) then
                    bigSubs(1) = subs(1) + grids_f%l_z(sps-1)
                    bigSubs(2) = subs(2) + grids_f%l_p(sps-1)
                  else if ( .not. frqDep(sps) ) then
                    bigSubs(1) = subs(2) + grids_f%l_z(sps-1)
                    bigSubs(2) = subs(3) + grids_f%l_p(sps-1)
                    nc = 3
                  else
                    bigSubs(1) = subs(1) + grids_f%l_f(sps-1)
                    bigSubs(2) = subs(2) + grids_f%l_z(sps-1)
                    bigSubs(3) = subs(3) + grids_f%l_p(sps-1)
                  end if
                  call showlist ( bigSubs(:nc-1) )
                else ! Show subscripts for one species as if in Vector_T%value3
                  call showlist ( subs(:nc-1) )
                end if
                call output ( eta3(j,k,l),  format='(1x,1pg14.6)' )
              end if
            end do
          end do
        end do
        if ( sawNZ ) call newLine
      end do
    end do
  contains
    subroutine ShowList ( A )
      integer, intent(in) :: A(:)
      integer :: S
      call output ( a(1), before=' (' )
      do s = 2, size(a)
        call output ( a(s), before=',' )
      end do
      call output ( ")" )
    end subroutine ShowList
  end subroutine Dump_Eta_Transpose

! ---------------------------------  ShowList (for dump routines)  -----
  subroutine ShowList ( J, Dims, Subs, After, Grids_f )
    use Load_SPS_Data_m, only: Grids_t
    use Output_m, only: Output
    ! Print "(" subs where dims/=1 ")"
    integer, intent(in) :: J, Dims(4), Subs(4)
    character(*), intent(in) :: After
    type (grids_t), intent(in), optional :: Grids_f ! to compute subscripts
    integer :: S
    if ( present(grids_f) ) then
      do s = 1, 3
        if ( dims(s) /= 1 ) exit
      end do
      call output ( subs(s), before="(" )
      do while ( s < 4 )
        s = s + 1
        if ( dims(s) /= 1 ) call output ( subs(s), before="," )
      end do
      call output ( ")" )
    else
      call output ( j, places=4, after=after )
    end if
  end subroutine ShowList

! --------------------------------------------------  Eta_Func_1d  -----
  pure function Eta_Func_1d ( Basis, Grid_Pt ) result ( Eta )

! Compute the eta matrix

! Inputs

    real(rp), intent(in) :: Basis(:) ! basis break points
    real(rp), intent(in) :: Grid_pt  ! grid point

! Outputs

    real(rp) :: Eta(size(basis))

! Internals

    integer(ip) :: J, N_coeffs
    real(rp) :: Del_basis

! The first coefficient is one for all values of grid_pt below basis(1)
! until grid_pt = basis(1),then it ramps down in the usual triangular sense.
! J is the coefficient index.

    n_coeffs = size(basis)
    eta = 0.0_rp

! first basis calculation

    if ( grid_pt <= basis(1) ) eta(1) = 1.0_rp

! Normal triangular function for j=2 to j=n_coeffs-1

    do j = 2 , n_coeffs
      del_basis = basis(j) - basis(j-1)
      if ( grid_pt <= basis(j) .and. &
        &  basis(j-1) <= grid_pt ) then
        eta(j-1) = (basis(j) - grid_pt) / del_basis
        eta(j) =   (grid_pt - basis(j-1)) / del_basis
      end if
    end do

! last basis calculation

    if ( grid_pt > basis(n_coeffs) ) eta(n_coeffs) = 1.0_rp

  end function Eta_Func_1d

! ------------------------------------------  Get_Column_Sparsity  -----
  subroutine Get_Column_Sparsity ( A, Nz, NNz )
    ! This isn't very efficient because it fills Nz and NNz with zeroes.
    ! If one knew that they had been filled with zeroes ab initio, and then
    ! kept properly, one could make Nz and NNz intent(inout) and use NNz
    ! to put zeroes into Nz.
    real(rp), intent(in) :: A(:,:)
    integer, intent(out) :: Nz(:,:) ! Nonzeroes in each column
    integer, intent(out) :: NNz(:)  ! Number of nonzeroes in each column
    integer :: I, J

    nnz = 0
    do j = 1, size(a,2)
      nz(:,j) = 0
      do i = 1, size(a,1)
        if ( a(i,j) /= 0.0 ) then
          nnz(j) = nnz(j) + 1
          nz(nnz(j),j) = i
        end if
      end do
    end do

  end subroutine Get_Column_Sparsity

! ----------------------------------------------  Get_Eta_1d_Hunt  -----
  subroutine Get_Eta_1d_Hunt ( Basis, Grid_Pt, Eta, IX )
  ! Find Grid_Pt in Basis, then compute Eta(1:2).  Use constant
  ! interpolation of Grid_Pt is below or above Basis.
    use MLSNumerics, only: Hunt

  ! Inputs

    real(rp), intent(in) :: Basis(:) ! basis break points
    real(rp), intent(in) :: Grid_pt  ! grid point

  ! Outputs

    real(rp), intent(out) :: Eta(2)
    integer, intent(out) :: IX ! if 1 <= IX < size(basis),
                               ! Basis(ix) <= Grid_Pt < Basis(ix+1)

  ! Internals

    integer(ip) :: N_coeffs

    if ( grid_pt < basis(1) ) then
      ix = 0
      eta = (/ 1.0, 0.0 /) ! constant extrapolation below range
      return
    end if

    n_coeffs = size(basis)
    if ( grid_pt >= basis(n_coeffs) ) then
      ix = n_coeffs
      eta = (/ 0.0, 1.0 /) ! constant extrapolation above range
      return
    end if

    ! basis(ix) <= grid_pt < basis(ix+1)
    call hunt ( basis, grid_pt, ix )
    eta(1) = ( grid_pt - basis(ix) ) / &
           & ( basis(ix+1) - basis(ix) )
    eta(2) = 1.0 - eta(1)

  end subroutine Get_Eta_1d_Hunt

!---------------------------------------------------  Eta_Func_2d  -----
  function Eta_Func_2d ( Basis, Grid, Sorted ) result ( Eta )

! Compute the eta matrix

    use Sort_m, only: SortP
! Inputs

    real(rp), intent(in) :: Basis(:) ! basis break points
    real(rp), intent(in) :: Grid(:)  ! grid values
    logical, optional, intent(in) :: Sorted ! "Grid is sorted"

! Outputs

    real(rp) :: Eta(size(grid),size(basis))

! Internals

    integer(ip) :: I, J, N_coeffs, N_Grid, P(size(grid)), PI
    real(rp) :: Del_basis
    logical :: MySorted

! Things below go more efficiently if Grid is sorted

    n_grid = size(grid)

    mySorted = n_grid == 1
    if ( .not. mySorted ) then
      if ( grid(1) <= grid(2) .and. present(sorted) ) mySorted = sorted
    end if

    if ( mySorted ) then
      do i = 1, n_grid
        p(i) = i
      end do
    else
      call sortp ( grid, 1, n_grid, p ) ! grid(p(:)) are now sorted
    end if

! The first coefficient is one for all values of grid_pt below basis(1)
! until grid_pt = basis(1),then it ramps down in the usual triangular sense.
! J is the coefficient index.

    n_coeffs = size(basis)
    eta = 0.0_rp

! first basis calculation

    do i = 1, n_grid
      pi = p(i)
      if ( grid(pi) > basis(1) ) exit
      eta(pi,1) = 1.0_rp
    end do

! Normal triangular function for j=2 to j=n_coeffs-1

    do j = 2 , n_coeffs
      del_basis = 1.0 / (basis(j) - basis(j-1))
      do while ( i <= n_grid )
        pi = p(i)
        if ( grid(pi) > basis(j) ) exit
        if ( basis(j-1) <= grid(pi) ) then
          eta(pi,j-1) = (basis(j) - grid(pi)) * del_basis
          eta(pi,j) =   (grid(pi) - basis(j-1)) * del_basis
        end if
        i = i + 1
      end do
    end do

! last basis calculation

    do i = n_grid, i, -1
      pi = p(i)
      if ( basis(n_coeffs) >= grid(pi) ) exit
      eta(pi,n_coeffs) = 1.0_rp
    end do

  end function Eta_Func_2d

!---------------------------------------------  Get_Eta_Sparse_1d  -----
  subroutine Get_Eta_Sparse_1d ( Basis, Grid_Pt, Eta )

! Compute the eta matrix

! Inputs

    real(rp), intent(in) :: Basis(:) ! basis break points
    real(rp), intent(in) :: Grid_pt  ! grid point

! Outputs

    real(rp), intent(out) :: Eta(:)  ! representation basis function

! Internals

    integer(ip) :: J, N_coeffs
    real(rp) :: Del_basis

! The first coefficient is one for all values of grid_pt below basis(1)
! until grid_pt = basis(1),then it ramps down in the usual triangular sense.
! J is the coefficient index.
    n_coeffs = size(basis)
    eta = 0.0_rp

! first basis calculation

    if ( grid_pt <= basis(1) ) eta(1) = 1.0_rp

! Normal triangular function for j=2 to j=n_coeffs-1

    do j = 2 , n_coeffs
      del_basis = 1.0 / (basis(j) - basis(j-1))
       if ( grid_pt <= basis(j) .and. &
         &  basis(j-1) <= grid_pt ) then
        eta(j-1) = (basis(j) - grid_pt) * del_basis
        eta(j) =   (grid_pt - basis(j-1)) * del_basis
      end if
    end do

! last basis calculation

    if ( basis(n_coeffs) < grid_pt ) eta(n_coeffs) = 1.0_rp

  end subroutine Get_Eta_Sparse_1d

!------------------------------------------  Get_Eta_Sparse_1d_fl  -----
  subroutine Get_Eta_Sparse_1d_fl ( Basis, Grid_Pt, Eta, First, Last )

! Compute the eta matrix

! Inputs

    real(rp), intent(in) :: Basis(:) ! basis break points
    real(rp), intent(in) :: Grid_pt  ! grid point

! Outputs

    real(rp), intent(out) :: Eta(:)  ! representation basis function
    integer, intent(out) :: First, Last ! First and last where Eta is not zero

! Internals

    integer(ip) :: J, N_coeffs
    real(rp) :: Del_basis

! The first coefficient is one for all values of grid_pt below basis(1)
! until grid_pt = basis(1),then it ramps down in the usual triangular sense.
! J is the coefficient index.
    n_coeffs = size(basis)
    eta = 0.0_rp

! first basis calculation

    if ( grid_pt <= basis(1) ) then
      eta(1) = 1.0
      first = 1
      last = 1
      return
    end if

! Normal triangular function for j=2 to j=n_coeffs-1.  Basis is sorted.

    do j = 2 , n_coeffs
      if ( grid_pt <= basis(j) .and. &
        &  basis(j-1) <= grid_pt ) then
        del_basis = 1.0 / (basis(j) - basis(j-1))
        eta(j-1) = (basis(j) - grid_pt) * del_basis
        eta(j) =   (grid_pt - basis(j-1)) * del_basis
        first = j - 1
        last = j
        return
      end if
    end do

! last basis calculation

    if ( basis(n_coeffs) < grid_pt ) then
      eta(n_coeffs) = 1.0_rp
      first = n_coeffs
      last = n_coeffs
    end if

  end subroutine Get_Eta_Sparse_1d_fl

!------------------------------------------  Get_Eta_Sparse_1d_nz  -----
  subroutine Get_Eta_Sparse_1d_nz ( Basis, Grid_Pt, Eta, Not_zero )

! Compute the eta matrix

! Inputs

    real(rp), intent(in) :: Basis(:)    ! basis break points
    real(rp), intent(in) :: Grid_pt     ! grid point

! Outputs

    real(rp), intent(out) :: Eta(:)     ! representation basis function
    logical, intent(out) :: Not_zero(:) ! where Eta is not zero

! Internals

    integer(ip) :: J, N_coeffs
    real(rp) :: Del_basis

! The first coefficient is one for all values of grid_pt below basis(1)
! until grid_pt = basis(1),then it ramps down in the usual triangular sense.
! J is the coefficient index.
    n_coeffs = size(basis)
    eta = 0.0_rp
    not_zero = .false.

! first basis calculation

    if ( grid_pt <= basis(1) ) then
      eta(1) = 1.0
      not_zero(1) = .true.
      return
    end if

! Normal triangular function for j=2 to j=n_coeffs-1.  Basis is sorted.

    do j = 2 , n_coeffs
      if ( grid_pt <= basis(j) .and. &
        &  basis(j-1) <= grid_pt ) then
        del_basis = 1.0 / (basis(j) - basis(j-1))
        eta(j-1) = (basis(j) - grid_pt) * del_basis
        eta(j) =   (grid_pt - basis(j-1)) * del_basis
        not_zero(j-1) = .true.
        not_zero(j) = .true.
        return
      end if
    end do

! last basis calculation

    if ( basis(n_coeffs) < grid_pt ) then
      eta(n_coeffs) = 1.0_rp
      not_zero(n_coeffs) = .true.
    end if

  end subroutine Get_Eta_Sparse_1d_nz

!---------------------------------------------  Get_Eta_Sparse_2d  -----
  subroutine Get_Eta_Sparse_2d ( Basis, Grid, Eta, Sorted )

! Compute the eta matrix.  Basis is assumed to be sorted.  Grid need not
! be sorted, but if it is, Sorted can be set .true. to suppress sorting it
! here.

    use Sort_m, only: SortP
! Inputs

    real(rp), intent(in) :: Basis(:) ! basis break points
    real(rp), intent(in) :: Grid(:)  ! grid values
    logical, optional, intent(in) :: Sorted ! "Grid is sorted"

! Outputs

    real(rp), intent(out) :: Eta(:,:) ! representation basis function

! Internals

    integer(ip) :: I, J, N_coeffs, N_Grid, P(size(grid)), PI
    real(rp) :: Del_basis
    logical :: MySorted

! The first coefficient is one for all values of grid below basis(1)
! until grid = basis(1),then it ramps down in the usual triangular sense,
! then the last coefficient is one for all values of grid above basis(n_coeffs).
! I is the independent variable grid index and J is the coefficient index

    n_coeffs = size(basis)
    eta = 0.0_rp
    ! eta(:,:n_coeffs) = 0.0_rp

! Things below go more efficiently if Grid is sorted

    n_grid = size(grid)

    mySorted = n_grid == 1
    if ( .not. mySorted ) then
      if ( grid(1) <= grid(2) .and. present(sorted) ) mySorted = sorted
    end if

    if ( mySorted ) then

! first basis calculation

      do i = 1, n_grid
        if ( grid(i) > basis(1) ) exit
        eta(i,1) = 1.0_rp
      end do

! Normal triangular function for j=2 to j=n_coeffs-1

      do j = 2 , n_coeffs
        del_basis = 1.0 / (basis(j) - basis(j-1))
        do while ( i <= n_grid )
          if ( grid(i) > basis(j) ) exit
          if ( basis(j-1) <= grid(i) ) then
            eta(i,j-1) = (basis(j) - grid(i)) * del_basis
            eta(i,j) =   (grid(i) - basis(j-1)) * del_basis
          end if
          i = i + 1
        end do
      end do

! last basis calculation

      do i = i, n_grid
        if ( basis(n_coeffs) >= grid(i) ) exit
        eta(i,n_coeffs) = 1.0_rp
      end do

    else

      call sortp ( grid, 1, n_grid, p ) ! grid(p(:)) are now sorted

! first basis calculation

      do i = 1, n_grid
        pi = p(i)
        if ( grid(pi) > basis(1) ) exit
        eta(pi,1) = 1.0_rp
      end do

! Normal triangular function for j=2 to j=n_coeffs-1

      do j = 2 , n_coeffs
        del_basis = 1.0 / (basis(j) - basis(j-1))
        do while ( i <= n_grid )
          pi = p(i)
          if ( grid(pi) > basis(j) ) exit
          eta(pi,j-1) = (basis(j) - grid(pi)) * del_basis
          eta(pi,j) =   (grid(pi) - basis(j-1)) * del_basis
          i = i + 1
        end do
      end do

! last basis calculation

      do i = i, n_grid
        pi = p(i)
        if ( basis(n_coeffs) >= grid(pi) ) exit
        eta(pi,n_coeffs) = 1.0_rp
      end do

    end if

  end subroutine Get_Eta_Sparse_2d

!------------------------------------------  Get_Eta_Sparse_2d_fl  -----
  subroutine Get_Eta_Sparse_2d_fl ( Basis, Grid, Eta, First, Last, Sorted )

! Compute the eta matrix.  Basis is assumed to be sorted.  Grid need not
! be sorted, but if it is, Sorted can be set .true. to suppress sorting it
! here.

!   use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use Sort_m, only: SortP

! Inputs

    real(rp), intent(in) :: Basis(:) ! basis break points
    real(rp), intent(in) :: Grid(:)  ! grid values
    logical, optional, intent(in) :: Sorted ! "Grid is sorted"

! Outputs

    real(rp), intent(out) :: Eta(:,:) ! representation basis function
    integer, intent(out) :: First(:)  ! First nonzero in each row of Eta
    integer, intent(out) :: Last(:)   ! Last nonzero in each row of Eta

! Internals

    integer(ip) :: I, J, N_coeffs, N_Grid, P(size(grid)), PI
    real(rp) :: Del_basis
    logical :: MySorted

! The first coefficient is one for all values of grid below basis(1)
! until grid = basis(1),then it ramps down in the usual triangular sense,
! then the last coefficient is one for all values of grid above basis(n_coeffs).
! I is the independent variable grid index and J is the coefficient index

    n_coeffs = size(basis)
    eta = 0.0_rp

! Things below go more efficiently if Grid is sorted

    n_grid = size(grid)

    mySorted = n_grid == 1
    if ( .not. mySorted ) then
      if ( grid(1) <= grid(2) .and. present(sorted) ) mySorted = sorted
    end if

    if ( mySorted ) then

! first basis calculation

      do i = 1, n_grid
        if ( grid(i) > basis(1) ) exit
        eta(i,1) = 1.0_rp
        first(i) = 1
        last(i) = 1
      end do

! Normal triangular function for j=2 to j=n_coeffs-1.  Both Basis and
! Grid(p(:)) are sorted, so we don't need to start from i=1.

      do j = 2, n_coeffs
        del_basis = 1.0_rp / ( basis(j) - basis(j-1) )
        do while ( i <= n_grid )
          if ( grid(i) > basis(j) ) exit
!         if ( basis(j-1) <= grid(i) ) then
            eta(i,j-1) = (basis(j) - grid(i)) * del_basis
            eta(i,j  ) = (grid(i) - basis(j-1)) * del_basis
            first(i) = j - 1
            last(i) = j
!         else
!           call MLSMessage ( MLSMSG_Error, ModuleName, &
!             & "Grid probably out of order" )
!         end if
          i = i + 1
        end do
      end do

! last basis calculation

      do i = i, n_grid
        eta(i,n_coeffs) = 1.0_rp
        first(i) = n_coeffs
        last(i) = n_coeffs
      end do

    else

      call sortp ( grid, 1, n_grid, p ) ! grid(p(:)) are now sorted

! first basis calculation

      do i = 1, n_grid
        pi = p(i)
        if ( grid(pi) > basis(1) ) exit
        eta(pi,1) = 1.0_rp
        first(pi) = 1
        last(pi) = 1
      end do

! Normal triangular function for j=2 to j=n_coeffs-1.  Both Basis and
! Grid(p(:)) are sorted, so we don't need to start from i=1.

      do j = 2, n_coeffs
        del_basis = 1.0_rp / ( basis(j) - basis(j-1) )
        do while ( i <= n_grid )
          pi = p(i)
          if ( grid(pi) > basis(j) ) exit
!         if ( basis(j-1) <= grid(pi) ) then
            eta(pi,j-1) = (basis(j) - grid(pi)) * del_basis
            eta(pi,j  ) = (grid(pi) - basis(j-1)) * del_basis
            first(pi) = j - 1
            last(pi) = j
!         else
!           call MLSMessage ( MLSMSG_Error, ModuleName, &
!             & "Why does the grid appear to be out of order?" )
!         end if
          i = i + 1
        end do
      end do

! last basis calculation

      do i = i, n_grid
        pi = p(i)
        eta(pi,n_coeffs) = 1.0_rp
        first(pi) = n_coeffs
        last(pi) = n_coeffs
      end do

    end if

  end subroutine Get_Eta_Sparse_2d_fl

!---------------------------------------  Get_Eta_Sparse_2d_fl_nz  -----
  subroutine Get_Eta_Sparse_2d_fl_nz ( Basis, Grid, Eta, Not_zero, First, Last, &
    & Sorted )

! Compute the eta matrix.  Basis is assumed to be sorted.  Grid need not
! be sorted, but if it is, Sorted can be set .true. to suppress sorting it
! here.

!   use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use Sort_m, only: SortP

! Inputs

    real(rp), intent(in) :: Basis(:) ! basis break points
    real(rp), intent(in) :: Grid(:)  ! grid values
    logical, optional, intent(in) :: Sorted ! "Grid is sorted"

! Outputs

    real(rp), intent(out) :: Eta(:,:) ! representation basis function
    logical, intent(out) :: Not_zero(:,:) ! where the above is not zero
    integer, intent(out) :: First(:)  ! First nonzero in each row of Eta
    integer, intent(out) :: Last(:)   ! Last nonzero in each row of Eta

! Internals

    integer(ip) :: I, J, N_coeffs, N_Grid, P(size(grid)), PI
    real(rp) :: Del_basis
    logical :: MySorted

! The first coefficient is one for all values of grid below basis(1)
! until grid = basis(1),then it ramps down in the usual triangular sense,
! then the last coefficient is one for all values of grid above basis(n_coeffs).
! I is the independent variable grid index and J is the coefficient index

    n_coeffs = size(basis)
    eta = 0.0_rp
    not_zero = .false.

! Things below go more efficiently if Grid is sorted

    n_grid = size(grid)

    mySorted = n_grid == 1
    if ( .not. mySorted ) then
      if ( grid(1) <= grid(2) .and. present(sorted) ) mySorted = sorted
    end if

    if ( mySorted ) then

! first basis calculation

      do i = 1, n_grid
        if ( grid(i) > basis(1) ) exit
        eta(i,1) = 1.0_rp
        not_zero(i,1) = .true.
        first(i) = 1
        last(i) = 1
      end do

! Normal triangular function for j=2 to j=n_coeffs-1.  Both Basis and
! Grid(p(:)) are sorted, so we don't need to start from i=1.

      do j = 2, n_coeffs
        del_basis = 1.0_rp / ( basis(j) - basis(j-1) )
        do while ( i <= n_grid )
          if ( grid(i) > basis(j) ) exit
!         if ( basis(j-1) <= grid(i) ) then
            eta(i,j-1) = (basis(j) - grid(i)) * del_basis
            eta(i,j  ) = (grid(i) - basis(j-1)) * del_basis
            not_zero(i,j-1) = .true.
            not_zero(i,j) = .true.
            first(i) = j - 1
            last(i) = j
!         else
!           call MLSMessage ( MLSMSG_Error, ModuleName, &
!             & "Grid probably out of order" )
!         end if
          i = i + 1
        end do
      end do

! last basis calculation

      do i = i, n_grid
        eta(i,n_coeffs) = 1.0_rp
        not_zero(i,n_coeffs) = .true.
        first(i) = n_coeffs
        last(i) = n_coeffs
      end do

    else

      call sortp ( grid, 1, n_grid, p ) ! grid(p(:)) are now sorted

! first basis calculation

      do i = 1, n_grid
        pi = p(i)
        if ( grid(pi) > basis(1) ) exit
        eta(pi,1) = 1.0_rp
        not_zero(pi,1) = .true.
        first(pi) = 1
        last(pi) = 1
      end do

! Normal triangular function for j=2 to j=n_coeffs-1.  Both Basis and
! Grid(p(:)) are sorted, so we don't need to start from i=1.

      do j = 2, n_coeffs
        del_basis = 1.0_rp / ( basis(j) - basis(j-1) )
        do while ( i <= n_grid )
          pi = p(i)
          if ( grid(pi) > basis(j) ) exit
!         if ( basis(j-1) <= grid(pi) ) then
            eta(pi,j-1) = (basis(j) - grid(pi)) * del_basis
            eta(pi,j  ) = (grid(pi) - basis(j-1)) * del_basis
            not_zero(pi,j-1) = .true.
            not_zero(pi,j) = .true.
            first(pi) = j - 1
            last(pi) = j
!         else
!           call MLSMessage ( MLSMSG_Error, ModuleName, &
!             & "Why does the grid appear to be out of order?" )
!         end if
          i = i + 1
        end do
      end do

! last basis calculation

      do i = i, n_grid
        pi = p(i)
        eta(pi,n_coeffs) = 1.0_rp
        not_zero(pi,n_coeffs) = .true.
        first(pi) = n_coeffs
        last(pi) = n_coeffs
      end do

    end if

  end subroutine Get_Eta_Sparse_2d_fl_nz

!-------------------------------------------  Get_Eta_Sparse_2d_nz  -----
  subroutine Get_Eta_Sparse_2d_nz ( Basis, Grid, Eta, Not_zero, Sorted )

! Compute the eta matrix.  Basis is assumed to be sorted.  Grid need not
! be sorted, but if it is, Sorted can be set .true. to suppress sorting it
! here.

    use Sort_m, only: SortP
! Inputs

    real(rp), intent(in) :: Basis(:) ! basis break points
    real(rp), intent(in) :: Grid(:)  ! grid values
    logical, optional, intent(in) :: Sorted ! "Grid is sorted"

! Outputs

    real(rp), intent(out) :: Eta(:,:) ! representation basis function
    logical, intent(out) :: Not_zero(:,:) ! where the above is not zero

! Internals

    integer(ip) :: I, J, N_coeffs, N_Grid, P(size(grid)), PI
    real(rp) :: Del_basis
    logical :: MySorted

! The first coefficient is one for all values of grid below basis(1)
! until grid = basis(1),then it ramps down in the usual triangular sense,
! then the last coefficient is one for all values of grid above basis(n_coeffs).
! I is the independent variable grid index and J is the coefficient index

    n_coeffs = size(basis)
    eta = 0.0_rp
    not_zero = .false.

! Things below go more efficiently if Grid is sorted

    n_grid = size(grid)

    mySorted = n_grid == 1
    if ( .not. mySorted ) then
      if ( grid(1) <= grid(2) .and. present(sorted) ) mySorted = sorted
    end if

    if ( mySorted ) then

! first basis calculation

      do i = 1, n_grid
        if ( grid(i) > basis(1) ) exit
        eta(i,1) = 1.0_rp
        not_zero(i,1) = .true.
      end do

! Normal triangular function for j=2 to j=n_coeffs-1.  Both Basis and
! Grid are sorted, so we don't need to start from i=1.

      do j = 2, n_coeffs
        del_basis = 1.0_rp / ( basis(j) - basis(j-1) )
        do while ( i <= n_grid )
          if ( grid(i) > basis(j) ) exit
!         if ( basis(j-1) <= grid(i) ) then
            eta(i,j-1) = (basis(j) - grid(i)) * del_basis
            eta(i,j  ) = (grid(i) - basis(j-1)) * del_basis
            not_zero(i,j-1) = .true.
            not_zero(i,j) = .true.
!         else
!           stop "grid not sorted"
!         end if
          i = i + 1
        end do
      end do

! last basis calculation

      do i = i, n_grid
        if ( basis(n_coeffs) >= grid(i) ) exit
        eta(i,n_coeffs) = 1.0_rp
        not_zero(i,n_coeffs) = .true.
      end do

    else

      call sortp ( grid, 1, n_grid, p ) ! grid(p(:)) are now sorted

! first basis calculation

      do i = 1, n_grid
        pi = p(i)
        if ( grid(pi) > basis(1) ) exit
        eta(pi,1) = 1.0_rp
        not_zero(pi,1) = .true.
      end do

! Normal triangular function for j=2 to j=n_coeffs-1.  Both Basis and
! Grid(p(:)) are sorted, so we don't need to start from i=1.

      do j = 2, n_coeffs
        del_basis = 1.0_rp / ( basis(j) - basis(j-1) )
        do while ( i <= n_grid )
          pi = p(i)
          if ( grid(pi) > basis(j) ) exit
          if ( basis(j-1) <= grid(pi) ) then
            eta(pi,j-1) = (basis(j) - grid(pi)) * del_basis
            eta(pi,j  ) = (grid(pi) - basis(j-1)) * del_basis
            not_zero(pi,j-1) = .true.
            not_zero(pi,j) = .true.
          end if
          i = i + 1
        end do
      end do

! last basis calculation

      do i = i, n_grid
        pi = p(i)
        if ( basis(n_coeffs) >= grid(pi) ) exit
        eta(pi,n_coeffs) = 1.0_rp
        not_zero(pi,n_coeffs) = .true.
      end do

    end if

  end subroutine Get_Eta_Sparse_2d_nz

! ----------------------------------------  Get_Eta_Column_Sparse  -----
  subroutine Get_Eta_Column_Sparse ( Basis, Grid, Eta, Row1, RowN, NZ, NNZ, Update )

  ! Compute rows Row1 to RowN of the eta matrix.  Basis is assumed to be
  ! sorted.  Rows Row1 to RowN of Grid are assumed to be sorted.  Row1 > RowN
  ! is allowed.  If Update is false, set nonzero elements of Eta to zero.
  ! It is assumed that Eta has initially been set to zero by the caller.

    use Sort_m, only: Sort

    real(rp), intent(in) :: Basis(:) ! basis break points
    real(rp), intent(in) :: Grid(:)  ! grid values
    real(rp), intent(inout) :: Eta(:,:)
    integer, intent(in) :: Row1, RowN
    integer, intent(inout) :: NZ(:,:) ! Locations of nonzeros in columns of Eta
    integer, intent(inout) :: NNZ(:)  ! Numbers of nonzeros in columns of Eta
    logical, intent(in) :: Update

    real(rp) :: Del_basis
    integer :: I, J, N_Coeffs, MyRow1, MyRowN

    if ( size(grid) == 0 ) return
    myRow1 = min(row1,ubound(grid,1))
    myRowN = min(rowN,ubound(grid,1))
    if ( grid(myRow1) > grid(myRowN) ) then ! grid is in opposite order.
      myRow1 = min(rowN,ubound(grid,1))
      myRowN = min(row1,ubound(grid,1))
    end if

    if ( .not. update ) then
      do i = 1, size(eta,2)
        eta(nz(:nnz(i),i),i) = 0.0
        nnz(i) = 0
      end do
    end if

    n_coeffs = size(basis)
    if ( n_coeffs <= 0 ) return

    if ( myRow1 <= myRowN ) then

      ! First basis calculation

      do i = myRow1, myRowN
        if ( grid(i) > basis(1) ) exit
        eta(i,1) = 1.0_rp
        nnz(1) = nnz(1) + 1
        nz(nnz(1),1) = i
      end do

      ! Normal triangular function for j=2 to j=n_coeffs-1.  Both Basis and
      ! Grid are sorted, so we don't need to start from i=myRow1.

      do j = 2, n_coeffs
        del_basis = 1.0_rp / ( basis(j) - basis(j-1) )
        do while ( i <= myRowN )
          if ( grid(i) > basis(j) ) exit
          eta(i,j-1) = (basis(j) - grid(i)) * del_basis
          if ( eta(i,j-1) /= 0.0 ) then
            nnz(j-1) = nnz(j-1) + 1
            nz(nnz(j-1),j-1) = i
          end if
          eta(i,j  ) = (grid(i) - basis(j-1)) * del_basis
          if ( eta(i,j) /= 0.0 ) then
            nnz(j) = nnz(j) + 1
            nz(nnz(j),j) = i
          end if
          i = i + 1
        end do
      end do

      ! last basis calculation

      do i = i, myRowN
        if ( basis(n_coeffs) >= grid(i) ) exit
        eta(i,n_coeffs) = 1.0_rp
        nnz(n_coeffs) = nnz(n_coeffs) + 1
        nz(nnz(n_coeffs),n_coeffs) = i
      end do

    else

      ! First basis calculation
      do i = myRow1, myRowN, -1
        if ( grid(i) > basis(1) ) exit
        eta(i,1) = 1.0_rp
        nnz(1) = nnz(1) + 1
        nz(nnz(1),1) = i
      end do

      ! Normal triangular function for j=2 to j=n_coeffs-1.  Both Basis and
      ! Grid are sorted, so we don't need to start from i=myRow1.

      do j = 2, n_coeffs
        del_basis = 1.0_rp / ( basis(j) - basis(j-1) )
        do while ( i >= myRowN )
          if ( grid(i) > basis(j) ) exit
          eta(i,j-1) = (basis(j) - grid(i)) * del_basis
          if ( eta(i,j-1) /= 0.0 ) then
            nnz(j-1) = nnz(j-1) + 1
            nz(nnz(j-1),j-1) = i
          end if
          eta(i,j  ) = (grid(i) - basis(j-1)) * del_basis
          if ( eta(i,j) /= 0.0 ) then
            nnz(j) = nnz(j) + 1
            nz(nnz(j),j) = i
          end if
          i = i - 1
        end do
      end do

      ! last basis calculation

      do i = i, myRowN, -1
        if ( basis(n_coeffs) >= grid(i) ) exit
        eta(i,n_coeffs) = 1.0_rp
        nnz(n_coeffs) = nnz(n_coeffs) + 1
        nz(nnz(n_coeffs),n_coeffs) = i
      end do

      ! Indices of nonzeros have come out in reverse order.
      ! If we are updating they're appended to those that were
      ! created in forward order.  It's simplest to sort them.
      do j = 1, n_coeffs
        call sort ( nz(:,j), 1, nnz(j) )
      end do

    end if

  end subroutine Get_Eta_Column_Sparse
  
! ----------------------------------  Get_Eta_Column_Sparse_fl_nz  -----
  subroutine Get_Eta_Column_Sparse_fl_nz ( Basis, Grid, Eta, Row1, RowN, &
    & NZ, NNZ, First, Last, Not_Zero )

  ! Compute rows Row1 to RowN of the eta matrix.  Basis is assumed to be
  ! sorted.  Rows Row1 to RowN of Grid are assumed to be sorted.  Row1 > RowN
  ! is allowed.   It is assumed that Eta has initially been set to zero by the
  ! caller.

    use Sort_m, only: Sort

    real(rp), intent(in) :: Basis(:) ! basis break points
    real(rp), intent(in) :: Grid(:)  ! grid values
    real(rp), intent(inout) :: Eta(:,:)
    integer, intent(in) :: Row1, RowN
    integer, intent(inout) :: NZ(:,:) ! Locations of nonzeros in columns of Eta
    integer, intent(inout) :: NNZ(:)  ! Numbers of nonzeros in columns of Eta
    integer, intent(out) :: First(:)  ! First nonzero in each row of Eta
    integer, intent(out) :: Last(:)   ! Last nonzero in each row of Eta
    logical, intent(inout) :: Not_zero(:,:) ! where Eta is not zero

    real(rp) :: Del_basis
    integer :: I, J, N_Coeffs, MyRow1, MyRowN

    myRow1 = row1
    myRowN = rowN
    if ( grid(row1) > grid(rowN) ) then ! grid is in opposite order.
      myRow1 = rowN
      myRowN = row1
    end if

    do i = 1, size(eta,2)
      eta(nz(:nnz(i),i),i) = 0.0
      not_zero(nz(:nnz(i),i),i) = .true.
      nnz(i) = 0
    end do

    n_coeffs = size(basis)
    if ( n_coeffs <= 0 ) return

    if ( myRow1 <= myRowN ) then

      ! First basis calculation

      do i = myRow1, myRowN
        if ( grid(i) > basis(1) ) exit
        eta(i,1) = 1.0_rp
        nnz(1) = nnz(1) + 1
        nz(nnz(1),1) = i
        not_zero(i,1) = .true.
        first(i) = 1
        last(i) = 1
      end do

      ! Normal triangular function for j=2 to j=n_coeffs-1.  Both Basis and
      ! Grid are sorted, so we don't need to start from i=myRow1.

      do j = 2, n_coeffs
        del_basis = 1.0_rp / ( basis(j) - basis(j-1) )
        do while ( i <= myRowN )
          if ( grid(i) > basis(j) ) exit
          eta(i,j-1) = (basis(j) - grid(i)) * del_basis
          eta(i,j  ) = (grid(i) - basis(j-1)) * del_basis
          nnz(j-1) = nnz(j-1) + 1
          nz(nnz(j-1),j-1) = i
          nnz(j) = nnz(j) + 1
          nz(nnz(j),j) = i
          not_zero(i,j-1) = .true.
          not_zero(i,j) = .true.
          first(i) = j - 1
          last(i) = j
          i = i + 1
        end do
      end do

      ! last basis calculation

      do i = i, myRowN
        eta(i,n_coeffs) = 1.0_rp
        nnz(n_coeffs) = nnz(n_coeffs) + 1
        nz(nnz(n_coeffs),n_coeffs) = i
        not_zero(i,n_coeffs) = .true.
        first(i) = n_coeffs
        last(i) = n_coeffs
      end do

    else

      ! First basis calculation
      do i = myRow1, myRowN, -1
        if ( grid(i) > basis(1) ) exit
        eta(i,1) = 1.0_rp
        nnz(1) = nnz(1) + 1
        nz(nnz(1),1) = i
        not_zero(i,1) = .true.
        first(i) = 1
        last(i) = 1
      end do

      ! Normal triangular function for j=2 to j=n_coeffs-1.  Both Basis and
      ! Grid are sorted, so we don't need to start from i=myRow1.

      do j = 2, n_coeffs
        del_basis = 1.0_rp / ( basis(j) - basis(j-1) )
        do while ( i >= myRowN )
          if ( grid(i) > basis(j) ) exit
          eta(i,j-1) = (basis(j) - grid(i)) * del_basis
          eta(i,j  ) = (grid(i) - basis(j-1)) * del_basis
          nnz(j-1) = nnz(j-1) + 1
          nz(nnz(j-1),j-1) = i
          nnz(j) = nnz(j) + 1
          nz(nnz(j),j) = i
          not_zero(i,j-1) = .true.
          not_zero(i,j) = .true.
          first(i) = j - 1
          last(i) = j
          i = i - 1
        end do
      end do

      ! last basis calculation

      do i = i, myRowN, -1
        eta(i,n_coeffs) = 1.0_rp
        nnz(n_coeffs) = nnz(n_coeffs) + 1
        nz(nnz(n_coeffs),n_coeffs) = i
        not_zero(i,n_coeffs) = .true.
        first(i) = n_coeffs
        last(i) = n_coeffs
      end do

      ! Indices of nonzeros have come out in reverse order.
      ! It's simplest to sort them.
      do j = 1, n_coeffs
        call sort ( nz(:,j), 1, nnz(j) )
      end do

    end if

  end subroutine Get_Eta_Column_Sparse_fl_nz

! ---------------------------------------------  Get_Eta_Stru_D_D  -----
  subroutine Get_Eta_Stru_D_D ( Basis, Grid, Eta )
    ! Compute linear interpolating coefficients from Basis to Grid. It is
    ! assumed that Basis(:) is sorted, but it is not assumed that Grid(:)
    ! is sorted.  It is assumed that Grid(:) changes slowly and smoothly.
    ! It is assumed that Size(Grid) == SIze(Eta).
    ! To interpolate from some T whose coordinates are Basis(:) to
    ! a value corresponding to Grid(k), compute
    ! Basis(eta(k)%l)*eta(k)%eta + Basis(eta(k)%l+1)*(1.0-eta(k)%eta).
    ! If Grid(k) < Basis(1) then Eta(k)%L = 1 and Eta(k)%Eta = 1.0.
    ! If Grid(k) >= Basis(size(basis)) then Eta(k)%L = size(basis-1) and
    ! Eta(k)%Eta = 0.0.  That is, outside the range of Basis(:), constant
    ! extrapolation is used.
    integer, parameter :: KB = kind(0.0d0)
    integer, parameter :: KG = kind(0.0d0)
    include "Get_Eta_Stru.f9h"
  end subroutine Get_Eta_Stru_D_D

! ---------------------------------------------  Get_Eta_Stru_D_S  -----
  subroutine Get_Eta_Stru_D_S ( Basis, Grid, Eta )
    ! Compute linear interpolating coefficients from Basis to Grid. It is
    ! assumed that Basis(:) is sorted, but it is not assumed that Grid(:)
    ! is sorted.  It is assumed that Grid(:) changes slowly and smoothly.
    ! It is assumed that Size(Grid) == SIze(Eta).
    ! To interpolate from some T whose coordinates are Basis(:) to
    ! a value corresponding to Grid(k), compute
    ! Basis(eta(k)%l)*eta(k)%eta + Basis(eta(k)%l+1)*(1.0-eta(k)%eta).
    ! If Grid(k) < Basis(1) then Eta(k)%L = 1 and Eta(k)%Eta = 1.0.
    ! If Grid(k) >= Basis(size(basis)) then Eta(k)%L = size(basis-1) and
    ! Eta(k)%Eta = 0.0.  That is, outside the range of Basis(:), constant
    ! extrapolation is used.
    integer, parameter :: KB = kind(0.0d0)
    integer, parameter :: KG = kind(0.0e0)
    include "Get_Eta_Stru.f9h"
  end subroutine Get_Eta_Stru_D_S

! ---------------------------------------------  Get_Eta_Stru_S_S  -----
  subroutine Get_Eta_Stru_S_S ( Basis, Grid, Eta )
    ! Compute linear interpolating coefficients from Basis to Grid. It is
    ! assumed that Basis(:) is sorted, but it is not assumed that Grid(:)
    ! is sorted.  It is assumed that Grid(:) changes slowly and smoothly.
    ! It is assumed that Size(Grid) == SIze(Eta).
    ! To interpolate from some T whose coordinates are Basis(:) to
    ! a value corresponding to Grid(k), compute
    ! Basis(eta(k)%l)*eta(k)%eta + Basis(eta(k)%l+1)*(1.0-eta(k)%eta).
    ! If Grid(k) < Basis(1) then Eta(k)%L = 1 and Eta(k)%Eta = 1.0.
    ! If Grid(k) >= Basis(size(basis)) then Eta(k)%L = size(basis-1) and
    ! Eta(k)%Eta = 0.0.  That is, outside the range of Basis(:), constant
    ! extrapolation is used.
    integer, parameter :: KB = kind(0.0e0)
    integer, parameter :: KG = kind(0.0e0)
    include "Get_Eta_Stru.f9h"
  end subroutine Get_Eta_Stru_S_S

! ---------------------------------------------------  Get_Eta_ZP  -----
  subroutine Get_Eta_ZP ( Z_Basis, Z_Path, Eta_Z, NZ_Z, NNZ_Z, Tan_Ind, &
                        & P_Basis, P_Path, Eta_P, NZ_P, NNZ_P, &
                        & Eta_ZP, NZ_ZP, NNZ_ZP, NonZero_ZP, Update )

    use GLNP, only: NG, NGP1
    real(rp), intent(in) :: Z_Basis(:)     ! Zeta basis break points
    real(rp), intent(in) :: Z_Path(:)      ! Zeta on the path
    real(rp), intent(inout) :: Eta_Z(:,:)  ! Interpolating coefficients from
                                           ! Z_Basis to Z_Path
    integer, intent(inout) :: NZ_Z(:,:)    ! Nonzeros in Eta_Z
    integer, intent(inout) :: NNZ_Z(:)     ! Numbers of nonzeros in columns of
                                           ! Eta_Z
    integer, intent(in) :: Tan_Ind         ! Index of tangent in Z_Path etc.
    real(rp), intent(in) :: P_Basis(:)     ! Phi basis break points
    real(rp), intent(in) :: P_Path(:)      ! Phi on the path
    real(rp), intent(inout) :: Eta_P(:,:)  ! Interpolating coefficients from
                                           ! P_Basis to P_Path
    integer, intent(inout) :: NZ_P(:,:)    ! Nonzeros in Eta_P
    integer, intent(inout) :: NNZ_P(:)     ! Numbers of nonzeros in columns of
                                           ! Eta_P
    real(rp), intent(inout) :: Eta_ZP(:,:) ! Eta_Z * Eta_P on the path for
                                           ! each state vector Z X P element
    integer, intent(inout) :: NZ_ZP(:,:)   ! Nonzeros in Eta_ZP
    integer, intent(inout) :: NNZ_ZP(:)    ! Numbers of nonzeros in columns of
                                           ! Eta_ZP
    logical, intent(out), optional :: NonZero_ZP(:,:) ! Where Eta_ZP is nonzero
    logical, intent(in), optional :: Update ! Clear nonzeros using NZ_, NNZ_
                                           ! for each Eta first, otherwise
                                           ! clear everything

    integer :: I
    logical :: MyUpdate

    myUpdate = .false.
    if ( present(update) ) myUpdate = update

    if ( myUpdate ) then
      ! Replace previously calculated nonzeros with zeros.
      do i = 1, size(nnz_p)
        eta_p(nz_p(:nnz_p(i),i),i) = 0
      end do
      do i = 1, size(nnz_z)
        eta_z(nz_z(:nnz_z(i),i),i) = 0
      end do
      do i = 1, size(nnz_zp)
        eta_zp(nz_zp(:nnz_zp(i),i),i) = 0
      end do
    else
      eta_p = 0
      eta_z = 0
      eta_zp = 0
    end if
    nnz_p = 0
    nnz_z = 0
    nnz_zp = 0

    call get_eta_sparse ( z_basis, z_path, eta_z, &
    &      tan_ind, 1, nz_z, nnz_z, .false. )
    ! Fine grid points between tangent points, if any, aren't used.
    eta_z(tan_ind+1:tan_ind+ng,:) = 0.0_rp
    call get_eta_sparse ( z_basis, z_path, eta_z, &
      &    tan_ind+ngp1, size(z_path), nz_z, nnz_z, .true. )
    call get_eta_sparse ( p_basis, p_path, eta_p, &
      &    1, size(p_path), nz_p, nnz_p, .false. )
    ! Fine grid points between tangent points, if any, aren't used.
    eta_p(tan_ind+1:tan_ind+ng,:) = 0.0_rp

    if ( present(nonZero_ZP) ) then
      call multiply_eta_column_sparse ( &
        & eta_z, nz_z, nnz_z, eta_p, nz_p, nnz_p, eta_zp, nz_zp, nnz_zp, &
        & nonZero_ZP )
    else
      call multiply_eta_column_sparse ( &
        & eta_z, nz_z, nnz_z, eta_p, nz_p, nnz_p, eta_zp, nz_zp, nnz_zp )
    end if

  end subroutine Get_Eta_ZP

! ----------------------------------------  Interpolate_Stru_1D_D  -----
  subroutine Interpolate_Stru_1D_D ( From, Eta, To )
    ! Interpolate from From to To using interpolating coefficients Eta
    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: From(:)
    type (eta_d_t), intent(in) :: Eta(:) ! Gotten from Get_Eta_Stru
    real(rk), intent(out) :: To(:)       ! Eta and To assumed same length
    include "Interpolate_Stru_1D.f9h"
  end subroutine Interpolate_Stru_1D_D

! ----------------------------------------  Interpolate_Stru_1D_S  -----
  subroutine Interpolate_Stru_1D_S ( From, Eta, To )
    ! Interpolate from From to To using interpolating coefficients Eta
    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: From(:)
    type (eta_s_t), intent(in) :: Eta(:) ! Gotten from Get_Eta_Stru
    real(rk), intent(out) :: To(:)       ! Eta and To assumed same length
    include "Interpolate_Stru_1D.f9h"
  end subroutine Interpolate_Stru_1D_S

! ----------------------------------------  Interpolate_Stru_2D_D  -----
  subroutine Interpolate_Stru_2D_D ( From, Eta_1, Eta_2, To )
    ! Interpolate from From to To using interpolating coefficients Eta
    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: From(:,:)
    type (eta_d_t), intent(in) :: Eta_1(:), Eta_2(:) ! Gotten from Get_Eta_Stru
    real(rk), intent(out) :: To(:)       ! Eta's and To assumed same length
    include "Interpolate_Stru_2D.f9h"
  end subroutine Interpolate_Stru_2D_D

! ----------------------------------------  Interpolate_Stru_2D_S  -----
  subroutine Interpolate_Stru_2D_S ( From, Eta_1, Eta_2, To )
    ! Interpolate from From to To using interpolating coefficients Eta
    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: From(:,:)
    type (eta_s_t), intent(in) :: Eta_1(:), Eta_2(:) ! Gotten from Get_Eta_Stru
    real(rk), intent(out) :: To(:)       ! Eta's and To assumed same length
    include "Interpolate_Stru_2D.f9h"
  end subroutine Interpolate_Stru_2D_S

! -----------------------------------  Multiply_Eta_Column_Sparse  -----
  subroutine Multiply_Eta_Column_Sparse (Eta_1, NZ_1, NNZ_1, Eta_2, NZ_2, NNZ_2, &
    & Eta_P, NZ_P, NNZ_P, Not_Zero_P )

  ! Multiply each column of Eta_1 by every column of Eta_2, each product
  ! giving a column of Eta_P, so that the number of columns of Eta_P is
  ! the product of the numbers of columns of Eta_1 and Eta_2.  Columns of
  ! Eta_1 are stepped more rapidly than columns of Eta_2.

  ! It is assumed that NNZ_P(i) and NZ_P(:NNZ_P(i),i) have defined
  ! values for each I from 1 to size(Eta_P,2).

    real(rp), intent(in) :: Eta_1(:,:), Eta_2(:,:) ! Values
    integer, intent(in) :: NZ_1(:,:), NZ_2(:,:)    ! Nonzero elements of Eta_*
    integer, intent(in) :: NNZ_1(:), NNZ_2(:)      ! Numbers of nonzeros
    real(rp), intent(inout) :: Eta_P(:,:)   ! Only the nonzeros are replaced
    integer, intent(inout) :: NZ_P(:,:)     ! Nonzero elements of Eta_P
    integer, intent(inout) :: NNZ_P(:)      ! Numbers of nonzeros
    logical, intent(inout), optional :: Not_Zero_P(:,:) ! Nonzero elements of Eta_P

    integer :: I1, I2, I3, J1, J2, K

    ! First remove nonzero elements from Eta_P
    if ( present(not_zero_p) ) then
      do k = 1, size(eta_p,2)
        not_zero_p(nz_p(:nnz_p(k),k),k) = .false.
      end do
    end if

    do k = 1, size(eta_p,2)
      eta_p(nz_p(:nnz_p(k),k),k) = 0.0
      nnz_p(k) = 0
    end do

    ! Multiply nonzero elements of Eta_1 and Eta_2
    ! Assume elements of each column of NZ_1 and NZ_2 are in order.
    k = 0

    do j2 = 1, size(eta_2,2)
      do j1 = 1, size(eta_1,2)
        k = k + 1
        i1 = 1
        i2 = 1
        do
          if ( i1 > nnz_1(j1) ) exit
          if ( i2 > nnz_2(j2) ) exit
          if ( nz_1(i1,j1) < nz_2(i2,j2) ) then
            i1 = i1 + 1
          else if ( nz_1(i1,j1) > nz_2(i2,j2) ) then
            i2 = i2 + 1
          else ! Both columns have nonzero in the same row
            i3 = nz_1(i1,j1)
            nnz_p(k) = nnz_p(k) + 1
            nz_p(nnz_p(k),k) = i3
            eta_p(i3,k) = eta_1(i3,j1) * eta_2(i3,j2)
            i1 = i1 + 1
            i2 = i2 + 1
          end if
        end do
      end do
    end do

    ! Note where nonzeros are
    if ( present(not_zero_p) ) then
      do k = 1, size(eta_p,2)
        not_zero_p(nz_p(:nnz_p(k),k),k) = .true.
      end do
    end if

  end subroutine Multiply_Eta_Column_Sparse

! -----------------------------------------------  Select_NZ_List  -----
  subroutine Select_NZ_List ( NZ_1, NNZ_1, List, I_start, NZ_2, NNZ_2 )
  ! Select the nonzeros from NZ_1 that are in List(I_start:) and put their
  ! indices in List into NZ_2.
    integer, intent(in) :: NZ_1(:,:)
    integer, intent(in) :: NNZ_1(:)
    integer, intent(in) :: List(:)
    integer, intent(in) :: I_Start
    integer, intent(out) :: NZ_2(:,:)
    integer, intent(out) :: NNZ_2(:)

    integer :: J
    integer :: I1, I2

    nnz_2 = 0
    ! This is just a merge
    do j = 1, size(nz_1,2)
      i1 = 1
      i2 = i_start
      do
        if ( i1 > nnz_1(j) .or. i2 > size(list) ) exit
        if ( nz_1(i1,j) < list(i2) ) then
          i1 = i1 + 1
        else if ( nz_1(i1,j) > list(i2) ) then
          i2 = i2 + 1
        else
          nnz_2(j) = nnz_2(j) + 1
          nz_2(nnz_2(j),j) = i2
          i1 = i1 + 1
          i2 = i2 + 1
        end if
      end do
    end do ! j

  end subroutine Select_NZ_List

!=========================================================================

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Get_Eta_Matrix_m
!---------------------------------------------------
! $Log$
! Revision 2.37  2018/08/07 18:10:15  vsnyder
! Make Dump_Eta_Transpose more like sparse dump
!
! Revision 2.36  2018/05/17 01:55:00  vsnyder
! Cannonball polishing
!
! Revision 2.35  2018/05/14 23:30:19  vsnyder
! Cannonball polishing
!
! Revision 2.34  2018/04/26 02:58:13  vsnyder
! size(Grids_f%l_z)==size(Grids_f%l_p) so N2 doesn't need to depend on myP
!
! Revision 2.33  2018/04/26 00:09:38  vsnyder
! Add Clean_Out_Nonzeroes, add P_Only to Dump_Eta_Transpose
!
! Revision 2.32  2018/03/21 01:58:10  vsnyder
! Restore NC for each sps in case it was changed by .not. frqDep
!
! Revision 2.31  2017/11/21 00:03:42  vsnyder
! Print the quantity name, print subscripts only for nonzer elements
!
! Revision 2.30  2017/09/20 01:06:04  vsnyder
! Embellish the dump
!
! Revision 2.29  2017/06/01 22:50:25  vsnyder
! Add Get_Column_Sparsity, more dump spiffing
!
! Revision 2.28  2017/05/24 20:29:49  vsnyder
! Spiff a dump
!
! Revision 2.27  2016/10/18 00:29:13  vsnyder
! Add Get_Eta_ZP
!
! Revision 2.26  2013/05/13 23:40:52  vsnyder
! Check for zero size grid
!
! Revision 2.25  2013/02/28 21:06:42  vsnyder
! Try not to access array out of bounds
!
! Revision 2.24  2011/06/23 01:13:31  vsnyder
! Add Get_Eta_1D_Hunt
!
! Revision 2.23  2010/06/07 23:17:20  vsnyder
! Add another representation for and interpolation using etas
!
! Revision 2.22  2009/06/23 18:26:11  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.21  2009/04/17 22:41:26  vsnyder
! Repair bug when n_grid == 1
!
! Revision 2.20  2009/01/16 23:44:15  vsnyder
! Add Dump_Eta_Column_Sparse.  Double-check whether things are sorted.
! Corrections to sparse calculations.
!
! Revision 2.19  2008/10/14 20:59:16  vsnyder
! Undo changes to Select_NZ_List for now
!
! Revision 2.18  2008/10/13 23:42:28  vsnyder
! Compute first,last correctly for 2d cases
!
! Revision 2.17  2007/06/25 20:36:08  vsnyder
! Add column sparse representation
!
! Revision 2.16  2007/06/06 01:15:43  vsnyder
! Use the PI variable instead of P(I)
!
! Revision 2.15  2007/01/20 01:06:39  vsnyder
! Add Get_Eta_Sparse_2d_fl_nz
!
! Revision 2.14  2007/01/19 02:38:29  vsnyder
! Separate paths for sorted and unsorted grid
!
! Revision 2.13  2006/12/13 21:53:11  vsnyder
! Delete declaration of unused variable
!
! Revision 2.12  2006/12/09 02:23:45  vsnyder
! Use generic instead of optional, add First, Last
!
! Revision 2.11  2006/06/29 02:51:30  vsnyder
! Correct a recently-introduced bug
!
! Revision 2.10  2006/06/29 01:44:02  vsnyder
! Make the 1d stuff work.  It wasn't used until metrics was revised
!
! Revision 2.9  2005/12/23 03:12:42  vsnyder
! Simplify last basis calculation
!
! Revision 2.8  2005/12/22 20:55:16  vsnyder
! Improve handling of sorted grid, some cannonball polishing
!
! Revision 2.7  2005/12/10 01:52:44  vsnyder
! Added Get_Eta_Sparse_1d, made Get_Eta_Sparse generic
!
! Revision 2.6  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.5  2003/06/20 02:02:11  vsnyder
! Sort GRID and merge with Basis instead of using WHERE
!
! Revision 2.4  2003/03/21 17:02:19  pwagner
! Initalizating whole eta, not_zero arrays instead of sections
!
! Revision 2.3  2002/10/08 17:08:04  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.2  2002/09/06 20:50:28  vsnyder
! Insert copyright notice, bound some array assignments
!
! Revision 2.1  2002/09/06 18:17:19  vsnyder
! Cosmetic changes, move USEs from module scope to procedure scope
!
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model

! Revision 1.1.2.3  2001/09/12 21:38:49  zvi
! Added CVS stuff


