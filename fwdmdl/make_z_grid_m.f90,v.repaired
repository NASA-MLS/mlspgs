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
  subroutine Make_Z_Grid ( zetas, z_grid, logp_eq_th, mult, last )

! This routine automatically makes an appropriate z_grid for the
! users specified input

    use Allocate_deallocate, only: Allocate_Test
    use MLSCommon, only: rp, ip
    use Sort_m, only: Sort

    ! inputs:

    real(rp), intent(in) :: zetas(:) ! this contains output pointing,
                        !  vertical bases for all considered species,
                        !  and a surface zeta (if desired) all concatenated
                        !  into one vector.
                        !  Elements in zetas do not need to be ordered.

    ! outputs:

    real(rp), allocatable, intent(out) :: z_grid(:) ! a suitable preselected
                        !  integration grid, with duplicates or near-duplicates
                        !  from zetas eliminated, and intermediate zetas
                        !  inserted if Mult is present and > 1.

    ! Keywords (Optional variables):

    real(rp), optional, intent(in) :: logp_eq_th ! threshold value for
                        !  eliminating duplicates
    integer(ip), optional, intent(in) :: mult ! A multiplicative factor for
                        !  increasing the z_grid density. N means have N
                        !  additional grid points between the minimal
                        !  necessary z_grid
    logical, optional, intent(in) :: Last ! Keep the last of a set of (nearly)
                        ! duplicate values.  Default: Keep the first one.

    ! Local variables:

    integer(ip) :: I      ! generic counter
    integer(ip) :: Keep(size(zetas)) ! indices of non-duplicates
    integer(ip) :: N      ! number of elements in minimal z_grid
    integer(ip) :: N_Ele  ! number of elements in zetas
    integer(ip) :: N_Grid ! number of elements in z_grid
    integer(ip) :: M      ! grid resolution factor

    logical :: First

    real(rp), allocatable :: Frac(:) ! fractions [1,2,...m-1]/m, for subgridding
    real(rp) :: Thresh    ! threshold for eliminating duplicates
    real(rp) :: Z(size(zetas))

    ! Sort the input

    n_ele = size(zetas)
    z = zetas
    call sort ( z, 1, n_ele )

    thresh = Default_Thresh
    if ( present(logp_eq_th) ) thresh = logp_eq_th

    ! Find non-duplicates.  This eliminates consecutive elements of Z that
    ! are closer together than the threshold.  If your input grid spacing is
    ! less than the threshold, it eliminates the entire grid, except for the
    ! first (last) one, even if the difference between the first and last is
    ! greater than the threshold. A more complicated method is required to
    ! get a minimum spacing of the threshold even if consecutive elements are
    ! all closer together.

    first = .true.
    if ( present(last) ) first = .not. last
    if ( first ) then
      n = 1
      keep(1) = 1
      do i = 2, n_ele
        if ( z(i) - z(i-1) > thresh ) then
          n = n + 1
          keep(n) = i
        end if
      end do
    else
      n = 0
      do i = 2, n_ele
        if ( z(i) - z(i-1) > thresh ) then
          n = n + 1
          keep(n) = i - 1
        end if
      end do
      if ( n == 0 ) then
        n = 1
        keep(1) = z(1)
      else if ( z(n_ele) - z(keep(n)) > thresh ) then
        n = n + 1
        keep(n) = n_ele
      end if
    end if

    ! Compute total grid points including subgrid points

    m = 1
    if ( present(mult) ) m = max(1,mult)
    n_grid = m * (n - 1) + 1
    call allocate_test ( z_grid, n_grid, 'Z_Grid', moduleName )

    ! Fill the grid

    z_grid(1:n_grid:m) = z(keep(1:n))

    if ( m > 1 ) then

      ! Insert subgridding.

      allocate ( frac(m-1) )
      frac = [ ( real(i,kind=rp), i = 1, m-1 ) ] / m
      do i = 1, n_grid-m, m
        z_grid(i+1:i+m-1) = z_grid(i) + ( z_grid(i+m) - z_grid(i) ) * frac
      end do

    end if

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
! Revision 2.14  2016/09/08 20:52:22  vsnyder
! Add 'last' argument
!
! Revision 2.13  2016/09/03 00:29:14  vsnyder
! Eliminate several array temps
!
! Revision 2.12  2014/01/11 01:28:53  vsnyder
! Decruftification
!
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
