!
! This module computes the eta matrix
!
MODULE get_eta_matrix_m
!
 use MLSCommon, only: RP, IP
!
 IMPLICIT none

 Private
 Public :: get_eta_sparse

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
 "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
 "$RCSfile$"
!---------------------------------------------------------------------------
CONTAINS
!---------------------------------------------------
 SUBROUTINE get_eta_sparse(basis,grid,eta,not_zero)
!
! Inputs
!
  REAL(rp), INTENT(in) :: basis(:) ! basis break points
  REAL(rp), INTENT(in) :: grid(:)  ! grid values
!
! Outputs
!
  REAL(rp), INTENT(out) :: eta(:,:) ! representation basis function
  LOGICAL, OPTIONAL, INTENT(out) :: not_zero(:,:) ! where the above is not zero
!
! Internals
!
  INTEGER(ip) :: j, n_coeffs
  REAL(rp) :: del_basis
!
! The first coefficient is one for all values of grid below basis(1)
! until grid = basis(1),then it ramps down in the usual triangular sense
! i is the independent variable grid index and j is the coefficient index
!
  n_coeffs = SIZE(basis)
  eta = 0.0_rp
  IF(PRESENT(not_zero)) THEN

    not_zero = .false.
!
! The wheres could be replaced with a search routine which would speed
! this up some more but you have to worry about sorting grid.
! first basis calculation
!
    WHERE(grid <= basis(1))
      eta(:,1) = 1.0_rp
      not_zero(:,1) = .true.
    ENDWHERE
!
! Normal triangular function for j=2 to j=n_coeffs-1
!
    DO j = 2 , n_coeffs
      del_basis = basis(j) - basis(j-1)
      WHERE(basis(j-1) <= grid .AND. grid <= basis(j))
        eta(:,j-1) = (basis(j) - grid)/ del_basis
        eta(:,j) =   (grid - basis(j-1))/ del_basis
        not_zero(:,j-1) = .true.
        not_zero(:,j) = .true.
      ENDWHERE
    ENDDO
!
! last basis calculation
!
    WHERE(basis(n_coeffs) < grid)
      eta(:,n_coeffs) = 1.0_rp
      not_zero(:,n_coeffs) = .TRUE.
    ENDWHERE
!
  ELSE
!
! first basis calculation
!
    WHERE(grid <= basis(1)) eta(:,1) = 1.0_rp
!
! Normal triangular function for j=2 to j=n_coeffs-1
!
    DO j = 2 , n_coeffs
      del_basis = basis(j) - basis(j-1)
      WHERE(basis(j-1) <= grid .AND. grid <= basis(j))
        eta(:,j-1) = (basis(j) - grid)/ del_basis
        eta(:,j) =   (grid - basis(j-1))/ del_basis
      ENDWHERE
    ENDDO
!
! last basis calculation
!
    WHERE(basis(n_coeffs) < grid) eta(:,n_coeffs) = 1.0_rp
!
  ENDIF
!
 END SUBROUTINE get_eta_sparse
!
END MODULE get_eta_matrix_m
!---------------------------------------------------
! $Log$
! Revision 1.1.2.3  2001/09/12 21:38:49  zvi
! Added CVS stuff
!
!
