! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Make_Z_Grid_m

  implicit NONE

  private
  public :: Make_Z_Grid

!---------------------------- RCS Ident Info -------------------------------
character (len=*), parameter :: IdParm = &
    & "$Id$"
character (len=len(idParm)) :: Id = IdParm
character (len=*), parameter :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains
!---------------------------------------------------------------------------
  subroutine Make_Z_Grid ( zetas, z_grid, tan_inds, logp_eq_th, mult )

! This routine automatically makes an appropriate z_grid for the
! users specified input

    use Allocate_deallocate, only: Allocate_test, Deallocate_test
    use MLSCommon, only: rp, ip
    use Sort_m, only: Sort

! inputs:

    real(rp), intent(in) :: zetas(:) ! this contains output pointing,
!                          vertical bases for all considered species,
!                          and a surface zeta (if desired) all concatenated
!                          into one vector.
!   elements in zetas do not need to be ordered.

! outputs:

    real(rp), dimension(:), pointer :: z_grid(:) ! a suitable
!                  preselected integration grid
    integer(ip), dimension(:), pointer:: tan_inds(:) ! a
!                  suitable set of pointing indicies.

! Keywords (Optional variables):

    real(rp), optional, intent(in) :: logp_eq_th ! threshold value for
!                                recognizing equality
    integer(ip), optional, intent(in) :: mult ! A multiplicative factor for
!                        increasing the z_grid density. ie 2 means have
!                        a grid point between the minimal necessary z_grid
! Local variables:

    logical :: mask(size(zetas))

    integer(ip) :: i ! generic counter
    integer(ip) :: n ! number of elements in minimal z_grid
    integer(ip) :: n_ele ! number of elements in zetas
    integer(ip) :: n_grid ! number of elements in z_grid
    integer(ip) :: m ! grid resolution factor

    real(rp) :: thresh ! equality threshold
    real(rp) :: z(size(zetas))
    real(rp), dimension(:), pointer :: frac
    real(rp), dimension(:), pointer :: z1
    real(rp), dimension(:), pointer :: del_z

! BEGIN CODE

! Sort the input

    n_ele = size(zetas)
    z = zetas
    call sort ( z, 1, n_ele )
    if ( present(logp_eq_th) ) then
      thresh = logp_eq_th
    else
      thresh = 0.0001_rp
    end if

    mask = (cshift(z,1) - z) > thresh
    mask(n_ele) = .true.
    n = count(mask)

! compute total grid points including multiplicative factors

    if ( present(mult) ) then
      m = mult
    else
      m = 1
    end if
    n_grid = m * (n - 1) + 1
    nullify ( z_grid, tan_inds )
    call allocate_test ( z_grid, n_grid, 'z_grid', modulename )

! Create a tangent pointing array index

    call allocate_test ( tan_inds, n_grid, 'tan_inds', modulename )
    tan_inds = (/(i,i = 1, n_grid)/)

! Fill the grid

    if ( m == 1 ) then ! not subgridding -- use the simple algorithm
      z_grid = pack(z,mask)
      return
    end if

! Fill in sub interval scale factors

    nullify ( frac, z1, del_z )
    call allocate_test ( z1, n, 'z1', modulename )
    call allocate_test ( del_z, n, 'del_z', modulename )
    call allocate_test ( frac, m, 'frac', ModuleName )

    z1 = pack(z,mask)
    del_z = cshift(z1,1) - z1
    frac = (/(real(i,kind=rp) / m, i = 0,m-1)/)

! Create and compute z_grid

    z_grid(1:n_grid-1) = reshape(spread(z1(1:n-1),1,m) &
                       +  spread(del_z(1:n-1),1,m) &
                       *  spread(frac,2,n-1), (/n_grid - 1/))
    z_grid(n_grid) = z1(n)

! Clean up

    call deallocate_test ( frac, 'frac', modulename )
    call deallocate_test ( del_z, 'del_z', modulename )
    call deallocate_test ( z1, 'z1', modulename )

  end subroutine Make_Z_Grid

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Make_Z_Grid_m

! $Log$
! Revision 2.5  2002/10/08 17:08:05  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.4  2002/10/07 19:27:50  vsnyder
! Use a simpler algorithm when not subgridding
!
! Revision 2.3  2002/10/04 00:04:09  vsnyder
! Move USE statements from module scope to procedure scope.  Insert copyright
! notice.  Cosmetic changes.
!
