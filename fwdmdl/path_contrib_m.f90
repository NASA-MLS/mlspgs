module Path_Contrib_M

  implicit NONE
  private
  public :: PATH_CONTRIB

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

!--------------------------------------------------  Path_Contrib  -----
! Estimate the contributions (along the path) of each interval of the
! (coarse) pre-selected integration grid.  Use that estimate to select
! where to do Gauss-Legendre quadrature.  Then allocate and fill
! arrays that control the Gauss-Legendre quadratures.

  subroutine Path_Contrib ( incoptdepth, e_rflty, tol, &
    &                       do_gl, gl_inds, gl_ndx )

    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use GLnp, only: NG
    use MLSCommon, only: RP, IP

  ! inputs

    real(rp), intent(in) :: incoptdepth(:) ! layer optical depth
    real(rp), intent(in) :: e_rflty        ! earth reflectivity
    real(rp), intent(in) :: tol            ! accuracy target in K

  ! outputs

    logical(ip), intent(out) :: do_gl(:)   ! set true for indicies to do
  !                                          gl computation
    integer, dimension(:), pointer :: GL_INDS ! Index of GL indices
    integer, dimension(:,:), pointer :: GL_NDX ! Packed Index array of GL intervals

  ! Internal stuff

    real(rp) :: dtaudn(size(incoptdepth))  ! path derivative of the
                                           ! transmission function
    integer(ip), dimension(:,:), pointer :: GL_INDGEN ! Temp. array of indeces
    integer(ip) :: i, i_tan, n_path, no_gl_ndx

    real(rp), parameter :: temp = 250.0_rp
    real(rp), parameter :: TolScale = 2.0_rp / temp ! 2.0 comes from centered
                                           ! difference used to compute dtaudn

    integer, parameter :: NGP1 = NG + 1
    integer, parameter :: GLIR(ng) = (/ (i, i = 1, ng ) /)  ! for > n_path/2
    integer, parameter :: GLIL(ng) = (/ (i ,i = -ng, -1) /) ! for <= n_path/2

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

    do_gl = dtaudn < -tol * tolscale

  ! The first and last index must be false

    do_gl((/1,n_path/)) = .FALSE.

  ! Allocate the output index arrays (and a temp array)

    nullify ( gl_inds, gl_ndx, gl_indgen )

    no_gl_ndx = count(do_gl(1:n_path))

    call allocate_test ( gl_inds, Ng * no_gl_ndx, 'gl_inds', moduleName )
    call allocate_test ( gl_ndx, no_gl_ndx, 2, 'gl_ndx', moduleName )
    call allocate_test ( gl_indgen, Ng, no_gl_ndx, 'gl_indgen', moduleName )

  ! Indices where GL is needed are ones where do_gl is true.

    gl_ndx(:,1) = pack((/(i,i=1,n_path)/),do_gl(1:n_path))

  ! Compute the gl indicies

    do i = 1 , no_gl_ndx
      if ( gl_ndx(i,1) > n_path/2 ) then
        gl_ndx(i,2) = 1 - Ng + Ngp1 * (gl_ndx(i,1) - 1)
        gl_indgen(:,i) = glir
      else
        gl_ndx(i,2) = 1 +      Ngp1 * (gl_ndx(i,1) - 1)
        gl_indgen(:,i) = glil
      end if
    end do

    gl_inds = reshape(spread(gl_ndx(:,2),1,Ng) + gl_indgen, (/Ng*no_gl_ndx/))

  ! Dispose of the temp array

    call deallocate_test ( gl_indgen, 'gl_indgen', moduleName )

  end subroutine Path_Contrib

!----------------------------------------------------------------------
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Path_Contrib_M

! $Log$
