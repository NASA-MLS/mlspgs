!
! This routine automatically makes an appropriate z_grid for the
! users specified input
!
MODULE make_z_grid_m
!
  USE MLSCommon, ONLY: rp, r4, r8, ip
  USE Sort_m
  USE Allocate_deallocate, only: Allocate_test, Deallocate_test
!
  IMPLICIT NONE
!
!---------------------------- RCS Ident Info -------------------------------
character (len=*), parameter, private :: IdParm = &
    & "$Id$"
character (len=len(idParm)) :: Id = IdParm
character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
!---------------------------------------------------------------------------
!
  CONTAINS
!---------------------------------------------------------------------------
  SUBROUTINE make_z_grid(zetas,z_grid,tan_inds,logp_eq_th,mult)
!
! inputs:
!
  REAL(rp), INTENT(in) :: zetas(:) ! This contains output pointing,
!                        vertical bases for all considered species,
!                        and a surface zeta (if desired) all concatenated
!                        into one vector.
! elements in zetas do not need to be ordered.
!
! outputs:
!
  REAL(rp), DIMENSION(:), POINTER :: z_grid(:) ! a suitable
!                preselected integration grid
  INTEGER(ip), DIMENSION(:), POINTER:: tan_inds(:) ! a
!                suitable set of pointing indicies.
!
! Keywords (Optional variables):
!
  REAL(rp), OPTIONAL, INTENT(in) :: logp_eq_th ! threshold value for
!                              recognizing equality
  INTEGER(ip), OPTIONAL, INTENT(in) :: mult ! A multiplicative factor for
!                      increasing the z_grid density. ie 2 means have
!                      a grid point between the minimal necessary z_grid
! Local variables:
!
  LOGICAL,  DIMENSION(:), POINTER :: mask
!
  INTEGER(ip) :: i ! generic counter
  INTEGER(ip) :: n ! number of elements in minimal z_grid
  INTEGER(ip) :: n_ele ! number of elements in zetas
  INTEGER(ip) :: n_grid ! number of elements in z_grid
  INTEGER(ip) :: m ! grid resolution factor
!
  REAL(rp) :: thresh ! Equality threshold
  REAL(rp), DIMENSION(:), POINTER :: z
  REAL(rp), DIMENSION(:), POINTER :: frac
  REAL(rp), DIMENSION(:), POINTER :: z1
  REAL(rp), DIMENSION(:), POINTER :: del_z
!
! BEGIN CODE

  nullify ( mask, z, frac, z1, del_z )
  nullify ( z_grid, tan_inds )
! Sort the input
!
  n_ele = SIZE(zetas)
  CALL ALLOCATE_TEST(z,n_ele,'z',ModuleName)
  CALL ALLOCATE_TEST(mask,n_ele,'mask',ModuleName)
  z = zetas
  CALL SORT(z,1,n_ele)
  IF (PRESENT(logp_eq_th)) THEN
    thresh = logp_eq_th
  ELSE
    thresh = 0.0001_rp
  ENDIF
!
  mask = (CSHIFT(z,1) - z) > thresh
  mask(n_ele) = .true.
  n = COUNT(mask)
  CALL ALLOCATE_TEST(z1,n,'z1',ModuleName)
  CALL ALLOCATE_TEST(del_z,n,'del_z',ModuleName)
  z1 = PACK(z,mask)
  del_z = CSHIFT(z1,1) - z1
!
! compute total grid points including multiplicative factors
!
  IF (PRESENT(mult)) THEN
    m = mult
  ELSE
    m = 1
  ENDIF
  n_grid = m * (n - 1) + 1
!
! Fill in sub interval points
!
  CALL ALLOCATE_TEST(frac,m,'frac',ModuleName)
  frac = (/(REAL(i,KIND=rp) / m,i = 0,m-1)/)
!
! Create and compute z_grid
!
  CALL ALLOCATE_TEST(z_grid,n_grid,'z_grid',ModuleName)
  z_grid(1:n_grid-1) = RESHAPE(SPREAD(z1(1:n-1),1,m) &
                     +  SPREAD(del_z(1:n-1),1,m) &
                     *  SPREAD(frac,2,n-1), (/n_grid - 1/))
  z_grid(n_grid) = z1(n)
!
! Create a tangent pointing array index
!
  CALL ALLOCATE_TEST(tan_inds,n_grid,'tan_inds',ModuleName)
  tan_inds = (/(i,i = 1, n_grid)/)
  CALL DEALLOCATE_TEST(z,'z',ModuleName)
  CALL DEALLOCATE_TEST(mask,'mask',ModuleName)
  CALL DEALLOCATE_TEST(z1,'z1',ModuleName)
  CALL DEALLOCATE_TEST(del_z,'del_z',ModuleName)
  CALL DEALLOCATE_TEST(frac,'frac',ModuleName)
!
  END SUBROUTINE make_z_grid
!
END MODULE make_z_grid_m
