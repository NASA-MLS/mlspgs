! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Make_Z_Grid_m

  implicit NONE

  private
  public :: Make_Z_Grid, Default_Thresh

  real, parameter :: Default_Thresh = 0.0001 ! Threshold for discarding dups

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains
!---------------------------------------------------------------------------
  subroutine Make_Z_Grid ( zetas, z_grid, logp_eq_th, mult )

! This routine automatically makes an appropriate z_grid for the
! users specified input

    use Allocate_deallocate, only: Allocate_test, Bytes, Test_Allocate
    use MLSCommon, only: rp, ip
    use Sort_m, only: Sort

! inputs:

    real(rp), intent(in) :: zetas(:) ! this contains output pointing,
!                          vertical bases for all considered species,
!                          and a surface zeta (if desired) all concatenated
!                          into one vector.
!                          Elements in zetas do not need to be ordered.

! outputs:

    real(rp), allocatable :: z_grid(:) ! a suitable preselected integration
!                          grid

! Keywords (Optional variables):

    real(rp), optional, intent(in) :: logp_eq_th ! threshold value for
!                          eliminating duplicates
    integer(ip), optional, intent(in) :: mult ! A multiplicative factor for
!                          increasing the z_grid density. N means have N
!                          additional grid points between the minimal
!                          necessary z_grid

! Local variables:

    logical :: mask(size(zetas))

    integer(ip) :: i      ! generic counter                      
    integer(ip) :: n      ! number of elements in minimal z_grid 
    integer(ip) :: n_ele  ! number of elements in zetas
    integer(ip) :: n_grid ! number of elements in z_grid
    integer(ip) :: m      ! grid resolution factor

    real(rp) :: thresh    ! threshold for eliminating duplicates
    real(rp) :: z(size(zetas))

! BEGIN CODE

! Sort the input

    n_ele = size(zetas)
    z = zetas
    call sort ( z, 1, n_ele )
    if ( present(logp_eq_th) ) then
      thresh = logp_eq_th
    else
      thresh = Default_Thresh
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
    allocate ( z_grid(n_grid), stat=i )
    call test_allocate ( i, moduleName, 'Z_Grid', [1], [n_grid], &
      & bytes(z_grid) )

! Fill the grid

    z_grid(1:n) = pack(z,mask)
    if ( m == 1 ) return ! not subgridding

! Create z_grid with subgridding: z_grid = z_grid + del_z * frac

    z_grid(n_grid) = z_grid(n)
    z_grid(1:n_grid-1) = &
      & reshape(spread(z_grid(1:n-1),1,m) + &
      &         spread(z_grid(2:n) - z_grid(1:n-1),1,m) * &
      &         spread((/(real(i,kind=rp) / m, i = 0,m-1)/),2,n-1), &
      &         (/n_grid - 1/))

  end subroutine Make_Z_Grid

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Make_Z_Grid_m

! $Log$
! Revision 2.11  2013/06/14 20:23:02  vsnyder
! Eliminate tan_press -- nobody used it
!
! Revision 2.10  2013/06/12 02:32:06  vsnyder
! Make z_grid allocatable instead of pointer
!
! Revision 2.9  2009/06/23 18:26:11  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.8  2007/01/12 21:45:03  vsnyder
! Make tolerance for discarding duplicates a parameter, publish it
!
! Revision 2.7  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.6  2003/09/19 18:10:51  vsnyder
! Make tangent point indices optional
!
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
