MODULE piq_int_m
  use MLSCommon, only: RP, IP
!
! This computes the piq (sans mass) used in the L2PC
! hydrostatic function
! This assumes the mean molecular mass is invariant
! between z_grid levels
! uses L2PC style triangular temperature representation basis
! argument z_ref allows a reference pressure that is /= z_grid(1)
! however adding this feature will cause program to run twice as slow
!
  IMPLICIT none
!
  Private
  Public :: piq_int
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
 "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
 "$RCSfile$"
!---------------------------------------------------------------------------
  CONTAINS
!---------------------------------------------------------------------------

SUBROUTINE piq_int(z_grid,t_basis,z_ref,piq)
!
  REAL(rp), INTENT(in) :: z_grid(:), t_basis(:), z_ref
  REAL(rp), INTENT(out) :: piq(:,:)
!
! inside code variables
!
  REAL(rp) :: a,c,aa,cc
  REAL(rp), DIMENSION(:), ALLOCATABLE :: b,d,bb,dd
  INTEGER(ip) :: n_lvls,n_coeffs,i,ind(1)
!
! begin code
! Establish dimensions
!
  n_lvls = SIZE(z_grid)
  n_coeffs = SIZE(t_basis)
!
! allocate limit arrays
!
  ALLOCATE(b(1:n_lvls))
  ALLOCATE(d(1:n_lvls))
  ALLOCATE(bb(1:n_lvls))
  ALLOCATE(dd(1:n_lvls))
!
! locate z_ref relative to t_basis
! I don't know if this is the fastest way to do this but this is a
! method derived from idl's ind_scl routine
!
  ind = MIN(MAX(PACK((/(i,i=0,n_coeffs)/), &
           (/MINVAL((/z_ref,(t_basis(i),i=1,n_coeffs)/))-1.0_rp, &
           (t_basis(i),i=1,n_coeffs)/) <= z_ref .AND. &
           z_ref < (/(t_basis(i),i=1,n_coeffs), &
           MAXVAL((/z_ref,(t_basis(i),i=1,n_coeffs)/))+1.0_rp/)),2),n_coeffs-2)
!
! NOTE it seems like MIN/MAX need to be done in pairs when comparing
! against arrays and scalers.
! for all coeffients below ind use
! initial coefficient
!
  a = MAX(t_basis(1),z_ref)
  b = MAX(MIN(MAX(z_grid,t_basis(1)),t_basis(2)),z_ref)
  c = MIN(z_ref,t_basis(2))
  d = MIN(MIN(MAX(z_grid,t_basis(1)),z_ref),t_basis(2))
  piq(:,1) = MAX(MIN(z_grid,t_basis(1)),z_ref) - z_ref &
           + MIN(MIN(z_grid,t_basis(1)),z_ref) - MIN(t_basis(1),z_ref) &
           + ((t_basis(2)-0.5_rp*(a+b))*(b-a) &
           -  (t_basis(2)-0.5_rp*(c+d))*(c-d))/(t_basis(2)-t_basis(1))
!
! these are wholly negative contributions only
!
  DO i = 2,ind(1) - 1
    b = MIN(MAX(z_grid,t_basis(i-1)),t_basis(i))
    d = MIN(MAX(z_grid,t_basis(i)),t_basis(i+1))
    piq(:,i) = (t_basis(i) - b)*(0.5_rp*(t_basis(i) - b) &
             / (t_basis(i) - t_basis(i-1)) - 1.0_rp) &
             - 0.5_rp*(t_basis(i+1) - d)**2 / (t_basis(i+1)-t_basis(i))
  END DO
!
! coefficients where z_ref is amongst t_basis
!
  DO i = ind(1),ind(1)+1
    a = MAX(t_basis(i-1),z_ref)
    b = MAX(MIN(MAX(z_grid,t_basis(i-1)),t_basis(i)),z_ref)
    c = MIN(z_ref,t_basis(i))
    d = MIN(MIN(MAX(z_grid,t_basis(i-1)),z_ref),t_basis(i))
    aa = MAX(t_basis(i),z_ref)
    bb = MAX(MIN(MAX(z_grid,t_basis(i)),t_basis(i+1)),z_ref)
    cc = MIN(z_ref,t_basis(i+1))
    dd = MIN(MIN(MAX(z_grid,t_basis(i)),z_ref),t_basis(i+1))
    piq(:,i) = ((0.5_rp*(a+b)-t_basis(i-1))*(b-a) &
             -  (0.5_rp*(c+d)-t_basis(i-1))*(c-d)) &
             / (t_basis(i)-t_basis(i-1)) &
             + ((t_basis(i+1)-0.5_rp*(aa+bb))*(bb-aa) &
             -  (t_basis(i+1)-0.5_rp*(cc+dd))*(cc-dd)) &
             / (t_basis(i+1)-t_basis(i))
  END DO
!
! for all coeffients above ind use
!
  DO i = ind(1)+2,n_coeffs-1
    b = MIN(MAX(z_grid,t_basis(i-1)),t_basis(i))
    d = MIN(MAX(z_grid,t_basis(i)),t_basis(i+1))
    piq(:,i) = 0.5_rp*(b - t_basis(i-1))**2 / (t_basis(i)-t_basis(i-1)) &
             + (d - t_basis(i))*(1.0 - 0.5_rp*(d - t_basis(i)) &
             / (t_basis(i+1)-t_basis(i)))
  END DO
!
! upper coefficient
!
  a = MAX(t_basis(n_coeffs-1),z_ref)
  b = MAX(MIN(MAX(z_grid,t_basis(n_coeffs-1)),t_basis(n_coeffs)),z_ref)
  c = MIN(z_ref,t_basis(n_coeffs))
  d = MIN(MIN(MAX(z_grid,t_basis(n_coeffs-1)),z_ref),t_basis(n_coeffs))
  piq(:,n_coeffs) = ((0.5_rp*(a+b)-t_basis(n_coeffs-1))*(b-a) &
                  -  (0.5_rp*(c+d)-t_basis(n_coeffs-1))*(c-d)) &
                  / (t_basis(n_coeffs)-t_basis(n_coeffs-1)) &
                  + MAX(MAX(z_grid,t_basis(n_coeffs)),z_ref) &
                  - MAX(t_basis(n_coeffs),z_ref) &
                  - MIN(z_ref,MAX(z_grid,t_basis(n_coeffs))) &
                  + MIN(t_basis(n_coeffs),z_ref)
!
! deallocate limit arrays
!
  DEALLOCATE(b)
  DEALLOCATE(d)
  DEALLOCATE(bb)
  DEALLOCATE(dd)
!
  END SUBROUTINE piq_int
!
END MODULE piq_int_m
!---------------------------------------------------
! $Log$
! Revision 1.1.2.3  2001/09/13 22:51:23  zvi
! Separating allocation stmts
!
! Revision 1.1.2.2  2001/09/12 21:38:52  zvi
! Added CVS stuff
!
