! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Get_Eta_Matrix_m

  use MLSCommon, only: RP, IP

  implicit NONE

  private
  public :: Get_Eta_Sparse

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    & "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------
  subroutine Get_Eta_Sparse ( Basis, Grid, Eta, Not_zero )

! Compute the eta matrix

    use Sort_m, only: SortP
! Inputs

    real(rp), intent(in) :: Basis(:) ! basis break points
    real(rp), intent(in) :: Grid(:)  ! grid values

! Outputs

    real(rp), intent(out) :: Eta(:,:) ! representation basis function
    logical, optional, intent(out) :: Not_zero(:,:) ! where the above is not zero

! Internals

    integer(ip) :: I, J, N_coeffs, N_Grid, P(size(grid))
    real(rp) :: Del_basis

! Things below go more efficiently if Grid is sorted

    n_grid = size(grid)
    call sortp ( grid, 1, n_grid, p ) ! grid(p(i)) are now sorted

! The first coefficient is one for all values of grid below basis(1)
! until grid = basis(1),then it ramps down in the usual triangular sense.
! I is the independent variable grid index and J is the coefficient index

    n_coeffs = size(basis)
    eta = 0.0_rp
    ! eta(:,:n_coeffs) = 0.0_rp
    if ( present(not_zero) ) then

      ! not_zero(:,:n_coeffs) = .false.
      not_zero = .false.

! first basis calculation

      do i = 1, n_grid
        if ( grid(p(i)) > basis(1) ) exit
        eta(p(i),1) = 1.0_rp
        not_zero(p(i),1) = .true.
      end do

! Normal triangular function for j=2 to j=n_coeffs-1.  Both Basis and
! Grid are sorted, so we don't need to start from i=1.

      do j = 2 , n_coeffs
        del_basis = basis(j) - basis(j-1)
        do while ( i <= n_grid )
          if ( grid(p(i)) > basis(j) ) exit
          if ( basis(j-1) <= grid(p(i)) ) then
            eta(p(i),j-1) = (basis(j) - grid(p(i))) / del_basis
            eta(p(i),j) =   (grid(p(i)) - basis(j-1)) / del_basis
            not_zero(p(i),j-1) = .true.
            not_zero(p(i),j) = .true.
          end if
          i = i + 1
        end do
      end do

! last basis calculation

      do while ( i <= n_grid )
        if ( basis(n_coeffs) < grid(p(i)) ) then
          eta(p(i),n_coeffs) = 1.0_rp
          not_zero(p(i),n_coeffs) = .true.
        end if
        i = i + 1
      end do

    else ! not_zero is not present

! first basis calculation

      do i = 1, n_grid
        if ( grid(p(i)) > basis(1) ) exit
        eta(p(i),1) = 1.0_rp
      end do

! Normal triangular function for j=2 to j=n_coeffs-1

      do j = 2 , n_coeffs
        del_basis = basis(j) - basis(j-1)
        do while ( i <= n_grid )
          if ( grid(p(i)) > basis(j) ) exit
          if ( basis(j-1) <= grid(p(i)) ) then
            eta(p(i),j-1) = (basis(j) - grid(p(i))) / del_basis
            eta(p(i),j) =   (grid(p(i)) - basis(j-1)) / del_basis
          end if
          i = i + 1
        end do
      end do

! last basis calculation

      do while ( i <= n_grid )
        if ( basis(n_coeffs) < grid(p(i)) ) then
          eta(p(i),n_coeffs) = 1.0_rp
        end if
        i = i + 1
      end do

    end if

  end subroutine Get_Eta_Sparse

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Get_Eta_Matrix_m
!---------------------------------------------------
! $Log$
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


