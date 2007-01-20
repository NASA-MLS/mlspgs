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

  use MLSCommon, only: RP, IP

  implicit NONE

  private
  public :: Eta_Func, Eta_Func_1d, Eta_Func_2d
  public :: Get_Eta_Sparse
  public :: Get_Eta_Sparse_1d, Get_Eta_Sparse_1d_fl, Get_Eta_Sparse_1d_nz
  public :: Get_Eta_Sparse_2d, Get_Eta_Sparse_2d_fl_nz, Get_Eta_Sparse_2d_nz

  interface Eta_Func
    module procedure Eta_Func_1d, Eta_Func_2d
  end interface Eta_Func

  interface Get_Eta_Sparse
    module procedure Get_Eta_Sparse_1d, Get_Eta_Sparse_1d_fl, Get_Eta_Sparse_1d_nz
    module procedure Get_Eta_Sparse_2d, Get_Eta_Sparse_2d_fl_nz, Get_Eta_Sparse_2d_nz
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains

!---------------------------------------------------  Eta_Func_1d  -----
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

    integer(ip) :: I, J, N_coeffs, N_Grid, P(size(grid))
    real(rp) :: Del_basis
    logical :: MySorted

! Things below go more efficiently if Grid is sorted

    mySorted = .false.
    if ( present(sorted) ) mySorted = sorted

    n_grid = size(grid)
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
      if ( grid(p(i)) > basis(1) ) exit
      eta(p(i),1) = 1.0_rp
    end do

! Normal triangular function for j=2 to j=n_coeffs-1

    do j = 2 , n_coeffs
      del_basis = 1.0 / (basis(j) - basis(j-1))
      do while ( i <= n_grid )
        if ( grid(p(i)) > basis(j) ) exit
        if ( basis(j-1) <= grid(p(i)) ) then
          eta(p(i),j-1) = (basis(j) - grid(p(i))) * del_basis
          eta(p(i),j) =   (grid(p(i)) - basis(j-1)) * del_basis
        end if
        i = i + 1
      end do
    end do

! last basis calculation

    do j = n_grid, i, -1
      if ( basis(n_coeffs) >= grid(p(j)) ) exit
      eta(p(j),n_coeffs) = 1.0_rp
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
    integer, intent(out) :: First, Last ! First and last where Etae is not zero

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

    integer(ip) :: I, J, N_coeffs, N_Grid, P(size(grid))
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

    mySorted = .false.
    if ( present(sorted) ) mySorted = sorted

    n_grid = size(grid)
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
        if ( grid(p(i)) > basis(1) ) exit
        eta(p(i),1) = 1.0_rp
      end do

! Normal triangular function for j=2 to j=n_coeffs-1

      do j = 2 , n_coeffs
        del_basis = 1.0 / (basis(j) - basis(j-1))
        do while ( i <= n_grid )
          if ( grid(p(i)) > basis(j) ) exit
          if ( basis(j-1) <= grid(p(i)) ) then
            eta(p(i),j-1) = (basis(j) - grid(p(i))) * del_basis
            eta(p(i),j) =   (grid(p(i)) - basis(j-1)) * del_basis
          end if
          i = i + 1
        end do
      end do

! last basis calculation

      do i = i, n_grid
        if ( basis(n_coeffs) >= grid(p(i)) ) exit
        eta(p(i),n_coeffs) = 1.0_rp
      end do

    end if

  end subroutine Get_Eta_Sparse_2d

!---------------------------------------  Get_Eta_Sparse_2d_fl_nz  -----
  subroutine Get_Eta_Sparse_2d_fl_nz ( Basis, Grid, Eta, Not_zero, First, Last, &
    & Sorted )

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

    mySorted = .false.
    if ( present(sorted) ) mySorted = sorted

    n_grid = size(grid)
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
          if ( basis(j-1) <= grid(i) ) then
            eta(i,j-1) = (basis(j) - grid(i)) * del_basis
            eta(i,j  ) = (grid(i) - basis(j-1)) * del_basis
            not_zero(i,j-1) = .true.
            not_zero(i,j) = .true.
            first(i) = j - 1
            last(i) = j
          end if
          i = i + 1
        end do
      end do

! last basis calculation

      do i = i, n_grid
        if ( basis(n_coeffs) >= grid(i) ) exit
        eta(i,n_coeffs) = 1.0_rp
        not_zero(i,n_coeffs) = .true.
        first(i) = n_coeffs
        last(i) = n_coeffs
      end do

    else

      call sortp ( grid, 1, n_grid, p ) ! grid(p(:)) are now sorted

! first basis calculation

      do i = 1, n_grid
        if ( grid(p(i)) > basis(1) ) exit
        eta(p(i),1) = 1.0_rp
        not_zero(p(i),1) = .true.
        first(p(i)) = 1
        last(p(i)) = 1
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
            first(p(i)) = j - 1
            last(p(i)) = j
          end if
          i = i + 1
        end do
      end do

! last basis calculation

      do i = i, n_grid
        if ( basis(n_coeffs) >= grid(p(i)) ) exit
        eta(p(i),n_coeffs) = 1.0_rp
        not_zero(p(i),n_coeffs) = .true.
        first(p(i)) = n_coeffs
        last(p(i)) = n_coeffs
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

    mySorted = .false.
    if ( present(sorted) ) mySorted = sorted

    n_grid = size(grid)
    if ( mySorted ) then

! first basis calculation

      do i = 1, n_grid
        if ( grid(i) > basis(1) ) exit
        eta(i,1) = 1.0_rp
        not_zero(i,1) = .true.
      end do

! Normal triangular function for j=2 to j=n_coeffs-1.  Both Basis and
! Grid(p(:)) are sorted, so we don't need to start from i=1.

      do j = 2, n_coeffs
        del_basis = 1.0_rp / ( basis(j) - basis(j-1) )
        do while ( i <= n_grid )
          if ( grid(i) > basis(j) ) exit
          if ( basis(j-1) <= grid(i) ) then
            eta(i,j-1) = (basis(j) - grid(i)) * del_basis
            eta(i,j  ) = (grid(i) - basis(j-1)) * del_basis
            not_zero(i,j-1) = .true.
            not_zero(i,j) = .true.
          end if
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
        if ( grid(p(i)) > basis(1) ) exit
        eta(p(i),1) = 1.0_rp
        not_zero(p(i),1) = .true.
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
        if ( basis(n_coeffs) >= grid(p(i)) ) exit
        eta(p(i),n_coeffs) = 1.0_rp
        not_zero(p(i),n_coeffs) = .true.
      end do

    end if

  end subroutine Get_Eta_Sparse_2d_nz

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Get_Eta_Matrix_m
!---------------------------------------------------
! $Log$
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


