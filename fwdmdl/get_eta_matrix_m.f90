
! This module computes the eta matrix

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
!---------------------------------------------------------------------------
contains
!---------------------------------------------------
  subroutine Get_Eta_Sparse ( Basis, Grid, Eta, Not_zero )

! Inputs

    real(rp), intent(in) :: Basis(:) ! basis break points
    real(rp), intent(in) :: Grid(:)  ! grid values

! Outputs

    real(rp), intent(out) :: Eta(:,:) ! representation basis function
    logical, optional, intent(out) :: Not_zero(:,:) ! where the above is not zero

! Internals

    integer(ip) :: J, N_coeffs
    real(rp) :: Del_basis

! The first coefficient is one for all values of grid below basis(1)
! until grid = basis(1),then it ramps down in the usual triangular sense.
! I is the independent variable grid index and J is the coefficient index

    n_coeffs = size(basis)
    eta = 0.0_rp
    if ( present(not_zero) ) then

      not_zero = .false.

! The wheres could be replaced with a search routine which would speed
! this up some more but you have to worry about sorting grid.
! first basis calculation

      where ( grid <= basis(1) )
        eta(:,1) = 1.0_rp
        not_zero(:,1) = .true.
      end where

! Normal triangular function for j=2 to j=n_coeffs-1

      do j = 2 , n_coeffs
        del_basis = basis(j) - basis(j-1)
        where ( basis(j-1) <= grid .and. grid <= basis(j) )
          eta(:,j-1) = (basis(j) - grid)/ del_basis
          eta(:,j) =   (grid - basis(j-1))/ del_basis
          not_zero(:,j-1) = .true.
          not_zero(:,j) = .true.
        end where
      end do

! last basis calculation

      where(basis(n_coeffs) < grid)
        eta(:,n_coeffs) = 1.0_rp
        not_zero(:,n_coeffs) = .true.
      end where

    else

! first basis calculation

      where ( grid <= basis(1) ) eta(:,1) = 1.0_rp

! Normal triangular function for j=2 to j=n_coeffs-1

      do j = 2 , n_coeffs
        del_basis = basis(j) - basis(j-1)
        where ( basis(j-1) <= grid .and. grid <= basis(j) )
          eta(:,j-1) = (basis(j) - grid)/ del_basis
          eta(:,j) =   (grid - basis(j-1))/ del_basis
        end where
      end do

! last basis calculation

      where ( basis(n_coeffs) < grid ) eta(:,n_coeffs) = 1.0_rp

    end if

  end subroutine Get_Eta_Sparse

end module Get_Eta_Matrix_m
!---------------------------------------------------
! $Log$
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model

! Revision 1.1.2.3  2001/09/12 21:38:49  zvi
! Added CVS stuff


