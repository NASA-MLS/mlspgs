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
    use MLSCommon, only: RK => RP, IP

  ! inputs

    real(rk), intent(in) :: incoptdepth(:) ! layer optical depth
    real(rk), intent(in) :: e_rflty        ! earth reflectivity
    real(rk), intent(in) :: tol            ! accuracy target in K

  ! outputs

    logical, intent(inout) :: do_gl(:)     ! TRUEs added for indicies to do
  !                                          gl computation

  ! Internal stuff

    real(rk) :: dtaudn(size(incoptdepth))  ! path derivative of the
                                           ! transmission function
    integer(ip) :: i, i_tan, n_path

    real(rk), parameter :: temp = 250.0_rk
    real(rk), parameter :: TolScale = 2.0_rk / temp ! 2.0 comes from centered
                                           ! difference used to compute dtaudn

  ! start code

    n_path = size(incoptdepth)
    i_tan = n_path / 2

  ! Compute the indefinite sum of (-incoptdepth).

    dtaudn(1) = 0.0_rk
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
    use CS_Expmat_M, only: CS_Expmat
    use MLSCommon, only: RK => RP, IP

  ! inputs

    complex(rk), intent(in) :: incoptdepth(:,:,:) ! layer optical depth
                                           ! (2,2,:)
    real(rk), intent(in) :: e_rflty        ! earth reflectivity
    real(rk), intent(in) :: tol            ! accuracy target in K

  ! outputs

    logical, intent(inout) :: do_gl(:)     ! set true for indicies to do
  !                                          gl computation

  ! Internal stuff

    complex(rk) :: deltau(2,2,size(incoptdepth,3))    ! exp(incoptdepth)
    complex(rk) :: dtaudn(2,2,size(incoptdepth,3))    ! path derivative of the
                                           ! transmission function
    complex(rk), parameter :: Ident(2,2) = reshape( (/ 1.0, 0.0, &
      &                                                0.0, 1.0 /), (/ 2,2 /) )

    complex(rk) :: tmp(2,2,size(incoptdepth,3)), tmp1(2,2,size(incoptdepth,3))

    integer(ip) :: i, i_tan, n_path

    real(rk), parameter :: temp = 250.0_rk
    real(rk), parameter :: TolScale = 2.0_rk / temp ! 2.0 comes from centered
                                           ! difference used to compute dtaudn

  ! start code

    n_path = size(incoptdepth,3)
    i_tan = n_path / 2

  ! Compute exp(incoptdepth) for all but the first and last levels

    do i = 1, n_path - 1
      call cs_expmat ( incoptdepth(:,:,i), deltau(:,:,i) )
    end do

    tmp(:,:,1) = ident
    tmp1(:,:,1) = ident

  ! Multiply the exp(incoptdepth) matrices together.

    do i = 2, i_tan
      tmp(:,:,i) =  matmul ( tmp(:,:,i-1),  deltau(:,:,i-1) )
      tmp1(:,:,i) = matmul ( tmp1(:,:,i-1), conjg(tmp(:,:,i)) )
    end do

    tmp(:,:,i_tan+1) = tmp(:,:,i_tan) * sqrt(e_rflty)
    tmp1(:,:,i_tan+1) = tmp1(:,:,i_tan) * e_rflty

    do i = i_tan+2, n_path
      tmp(:,:,i) =  matmul ( tmp(:,:,i-1),  deltau(:,:,i-1) )
      tmp1(:,:,i) = matmul ( tmp1(:,:,i-1), conjg(tmp(:,:,i)) )
    end do

  ! Where is the derivative of that product large?

    dtaudn = 0.5_rk * ( eoshift(tmp1(:,:,1:n_path),+1,dim=3) - &
                      & eoshift(tmp1(:,:,1:n_path),-1,dim=3) )

    do i = 2, n_path-1
      if ( any(abs(real(dtaudn(:,:,i)))  > tolScale*tol ) .or. &
           any(abs(aimag(dtaudn(:,:,i))) > tolScale*tol ) ) do_gl(i) = .true.
    end do

  end subroutine Path_Contrib_Polarized

  ! ------------------------------------------------  Get_GL_inds  -----
  subroutine Get_GL_inds ( Do_GL, GL_Inds, NGL )
  ! Fill the array that controls application of GL

    use GLnp, only: NG
    use MLSCommon, only: IP

    logical, intent(inout) :: DO_GL(:)         ! Set true for indicies to do
                                               ! gl computation.  First and
                                               ! last are set false here.
    integer(ip), intent(out) :: GL_INDS(:)     ! Indices of where to do GL
    integer(ip), intent(out) :: NGL            ! How much of GL_INDS to use

    integer :: I, N_PATH

    integer, parameter :: NGP1 = NG + 1
    integer, parameter :: GLIR(ng) = (/ (i, i = 1, ng ) /)  ! for > n_path/2
    integer, parameter :: GLIL(ng) = (/ (i ,i = -ng, -1) /) ! for <= n_path/2

    n_path = size(do_gl)

  ! The first and last index must be false

    do_gl((/1,n_path/)) = .FALSE.

    ngl = 0
    do i = 2, n_path-1 ! first and last elements of do_gl are false
      if ( do_gl(i) ) then
        ngl = ngl + ng
        if ( i > n_path / 2 ) then
          gl_inds(ngl-ng+1:ngl) = 1 - Ng + Ngp1 * (i - 1) + glir
        else
          gl_inds(ngl-ng+1:ngl) = 1 +      Ngp1 * (i - 1) + glil
        end if
      end if
    end do

  end subroutine Get_GL_inds

!-----------------------------------------------------------------------
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Path_Contrib_M

! $Log$
! Revision 2.6  2003/02/03 22:56:20  vsnyder
! Get rid of an array temp
!
! Revision 2.5  2003/01/31 01:53:28  vsnyder
! Calculate where to do GL with one less array temp
!
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
