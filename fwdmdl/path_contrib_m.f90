module Path_Contrib_M

  implicit NONE
  private
  public :: Get_GL_Inds, PATH_CONTRIB

  interface PATH_CONTRIB
    module procedure PATH_CONTRIB_SCALAR, PATH_CONTRIB_POLARIZED
  end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

!-------------------------------------------  Path_Contrib_Scalar  -----
! Estimate the contributions (along the path) of each interval of the
! (coarse) pre-selected integration grid.  Use that estimate to select
! where to do Gauss-Legendre quadrature.  Then allocate and fill
! arrays that control the Gauss-Legendre quadratures.

  subroutine Path_Contrib_Scalar ( incoptdepth, e_rflty, tol, do_gl )

    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use MLSCommon, only: RP, IP

  ! inputs

    real(rp), intent(in) :: incoptdepth(:) ! layer optical depth
    real(rp), intent(in) :: e_rflty        ! earth reflectivity
    real(rp), intent(in) :: tol            ! accuracy target in K

  ! outputs

    logical(ip), intent(inout) :: do_gl(:) ! TRUEs added for indicies to do
  !                                          gl computation

  ! Internal stuff

    real(rp) :: dtaudn(size(incoptdepth))  ! path derivative of the
                                           ! transmission function
    integer(ip) :: i, i_tan, n_path

    real(rp), parameter :: temp = 250.0_rp
    real(rp), parameter :: TolScale = 2.0_rp / temp ! 2.0 comes from centered
                                           ! difference used to compute dtaudn

  ! start code

    n_path = size(incoptdepth)
    i_tan = n_path / 2

  ! Compute the indefinite sum of (-incoptdepth).

    dtaudn(1) = 0.0_rp
    do i = 2 , i_tan
      dtaudn(i) = dtaudn(i-1) - incoptdepth(i)
    end do

    dtaudn(i_tan+1) = dtaudn(i_tan)

    do i = i_tan+2, n_path
      dtaudn(i) = dtaudn(i-1) - incoptdepth(i-1)
    end do

  ! compute the tau path derivative ~ exp(Tau) dTau/ds.

    dtaudn = (eoshift(dtaudn,1,dtaudn(n_path)) -             &
              eoshift(dtaudn,-1,dtaudn(1))) * exp(dtaudn)

    dtaudn(i_tan+1:n_path) = dtaudn(i_tan+1:n_path) * e_rflty

  ! find where the tau derivative is large.  Remember, we've
  ! been subtracting, so "large" means "large and negative."

    where ( dtaudn < -tol * tolscale ) do_gl = .true.

  end subroutine Path_Contrib_Scalar

!----------------------------------------  Path_Contrib_Polarized  -----
! Estimate the contributions (along the path) of each interval of the
! (coarse) pre-selected integration grid.  Use that estimate to select
! where to do Gauss-Legendre quadrature.  Then allocate and fill
! arrays that control the Gauss-Legendre quadratures.

  subroutine Path_Contrib_Polarized ( incoptdepth, e_rflty, tol, do_gl )

    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use MLSCommon, only: RP, IP

  ! inputs

    complex(rp), intent(in) :: incoptdepth(:,:) ! layer optical depth
                                           ! First dimension is 3
    real(rp), intent(in) :: e_rflty        ! earth reflectivity
    real(rp), intent(in) :: tol            ! accuracy target in K

  ! outputs

    logical(ip), intent(out) :: do_gl(:)       ! set true for indicies to do
  !                                              gl computation

  ! Internal stuff

    complex(rp) :: dtaudn(3,size(incoptdepth,2))    ! path derivative of the
                                           ! transmission function
    integer(ip) :: i, i_tan, n_path

    real(rp), parameter :: temp = 250.0_rp
    real(rp), parameter :: TolScale = 2.0_rp / temp ! 2.0 comes from centered
                                           ! difference used to compute dtaudn

  ! start code

    n_path = size(incoptdepth,2)
    i_tan = n_path / 2

  ! Compute the indefinite sum of (-incoptdepth).

    dtaudn(:,1) = 0.0_rp
    do i = 2 , i_tan
      dtaudn(:,i) = dtaudn(:,i-1) - incoptdepth(:,i)
    end do

    dtaudn(:,i_tan+1) = dtaudn(:,i_tan)

    do i = i_tan+2, n_path
      dtaudn(:,i) = dtaudn(:,i-1) - incoptdepth(:,i-1)
    end do

  ! compute the tau path derivative ~ exp(Tau) dTau/ds.

    dtaudn = (eoshift(dtaudn,+1,dim=2) -             &
      &       eoshift(dtaudn,-1,dim=2) ) * exp(dtaudn)

    dtaudn(:,i_tan+1:n_path) = dtaudn(:,i_tan+1:n_path) * e_rflty

  ! find where the tau derivative is large.

    do i = 2, n_path - 1
      if ( max(maxval(abs( real(dtaudn(:,i)))), &
        &      maxval(abs(aimag(dtaudn(:,i))))) > tol * tolscale ) &
        & do_gl(i) = .true.
    end do

  end subroutine Path_Contrib_Polarized

  ! ------------------------------------------------  Get_GL_inds  -----
  subroutine Get_GL_inds ( Do_GL, GL_Inds, GL_Ndx, NGL )
  ! Create and fill the arrays that control application of GL

    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use GLnp, only: NG
    use MLSCommon, only: IP

    logical(ip), intent(inout) :: DO_GL(:)     ! Set true for indicies to do
                                               ! gl computation.  First and
                                               ! last are set false here.
    integer, dimension(:), pointer :: GL_INDS  ! Index of GL indices
    integer, dimension(:,:), pointer :: GL_NDX ! Packed Index array of GL intervals
    integer, intent(out) :: NGL                ! How much of GL_INDS to use

    integer :: I, NO_GL_NDX, N_PATH

    integer, parameter :: NGP1 = NG + 1
    integer, parameter :: GLIR(ng) = (/ (i, i = 1, ng ) /)  ! for > n_path/2
    integer, parameter :: GLIL(ng) = (/ (i ,i = -ng, -1) /) ! for <= n_path/2

  ! Allocate the output index arrays (and a temp array)

    nullify ( gl_inds, gl_ndx )

    n_path = size(do_gl)

  ! The first and last index must be false

    do_gl((/1,n_path/)) = .FALSE.

    no_gl_ndx = count(do_gl)
    ngl = Ng * no_gl_ndx

    call allocate_test ( gl_inds, ngl, 'gl_inds', moduleName )
    call allocate_test ( gl_ndx, no_gl_ndx, 2, 'gl_ndx', moduleName )

  ! Indices where GL is needed are ones where do_gl is true.

    gl_ndx(:,1) = pack((/(i,i=1,n_path)/),do_gl)

  ! Compute the gl indicies

    do i = 1 , no_gl_ndx
      if ( gl_ndx(i,1) > n_path/2 ) then
        gl_ndx(i,2) = 1 - Ng + Ngp1 * (gl_ndx(i,1) - 1)
        gl_inds(ng*(i-1)+1:ng*i) = gl_ndx(i,2) + glir
      else
        gl_ndx(i,2) = 1 +      Ngp1 * (gl_ndx(i,1) - 1)
        gl_inds(ng*(i-1)+1:ng*i) = gl_ndx(i,2) + glil
      end if
    end do

  end subroutine Get_GL_inds

!-----------------------------------------------------------------------
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Path_Contrib_M

! $Log$
! Revision 2.4  2003/01/30 19:31:18  vsnyder
! Undo change that didn't work -- tried to compute gl_inds without array temp
!
! Revision 2.3  2003/01/18 02:22:58  vsnyder
! IMAG should have been AIMAG
!
! Revision 2.2  2003/01/18 01:42:42  vsnyder
! Added complex 3-quantity path_contrib subroutine.
! Separated Get_GL_Inds into a separate subroutine.
!
! Revision 2.1  2003/01/08 00:09:56  vsnyder
! Moved from rad_tran, where it didn't seem to belong
!
