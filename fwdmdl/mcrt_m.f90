! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MCRT_m

!       Magnetic    Condensed    Radiative    Transfer

  implicit NONE
  private
  public :: MCRT

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    &  "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

! ---------------------------------------------------------  Mcrt  -----
  subroutine Mcrt ( t_script, sqrt_earth_ref, del_tau, tau, radiance )

!       Magnetic    Condensed    Radiative    Transfer

! Compute the polarized radiative transfer using a condensed algorithm.
! Enter t_script(*) in a quantity that is linear in power.
! where T is intensity of a black body radiator having temperature t
! at frequency f (usually entered in units of kelvin being derived from
! T_{i} = h v / [ k ( exp{h v / (k t_{i}) } - 1)]).
! sqrt_earth_ref is the square root of Earth's reflectivity.

    use CS_Expmat_M, only: CS_Expmat
    use MLSCommon, only: Rk => Rp

    real(rk), intent(in) :: t_script(:)
    real(rk), intent(in) :: sqrt_earth_ref
    complex(rk), intent(inout) :: del_tau(:,:,:) ! 2 x 2 x : = exp(del_opcty)
                                                 ! Called P in the notes
                                                 ! The middle one is changed
    complex(rk), intent(out) :: tau(:,:,:)       ! 2 x 2 x :
    complex(rk), intent(out) :: radiance(2,2)

    integer :: i, i_tan, n_path

    complex(rk), parameter :: Ident(2,2) = reshape( (/ 1.0, 0.0, &
                                                     & 0.0, 1.0 /), (/2,2/) )

    n_path = size(t_script)

! Initialize segment 1 calculation

    tau(:,:,1) = ident
    ! radiance = t_script(1) * ident
    radiance(1,1) = t_script(1)
    radiance(1,2) = 0.0
    radiance(2,1) = 0.0
    radiance(2,2) = t_script(1)

! Proceed with first segment integration

    i_tan = n_path / 2
    do i = 2, i_tan
      tau(:,:,i) = matmul ( tau(:,:,i-1), del_tau(:,:,i) )
      call updaterad ( radiance, t_script(i), tau(:,:,i) )
    end do

! Tangent point (or Earth intersecting) layer.  If it's not
! an Earth intersecting layer, sqrt_earth_ref will be 1.0.

    del_tau(:,:,i_tan+1) = del_tau(:,:,i_tan)
    tau(:,:,i_tan+1) = sqrt_earth_ref * tau(:,:,i_tan)
    call updaterad ( radiance, t_script(i_tan+1), tau(:,:,i_tan+1) )

! Proceed with third segment integration, which includes the
! space radiance contribution.

    do i = i_tan+2, n_path
      tau(:,:,i) = matmul ( tau(:,:,i-1) , del_tau(:,:,i-1) )
      call updaterad ( radiance, t_script(i), tau(:,:,i) )  
    end do                                                       

  contains
  ! ..................................................  Updaterad  .....

    subroutine Updaterad ( Radiance, Scalar, Tau )

! Update the radiance from the tau:
! Radiance = Radiance + Scalar * Tau * conjg(transpose(Tau))
! We know that Radiance is Hermitian, so the diagonal elements are
! real and the off-diagonal elements are conjugates.  Exploiting
! this saves roughly half of the work.

      complex(rk), intent(inout) :: Radiance(2,2)
      real(rk), intent(in) :: Scalar
      complex(rk), intent(in) :: Tau(2,2)

      real(rk) t11r, t11i, t12r, t12i, t21r, t21i, t22r, t22i

      t11r = real(tau(1,1))
      t11i = aimag(tau(1,1))
      t12r = real(tau(1,2))
      t12i = aimag(tau(1,2))
      t21r = real(tau(2,1))
      t21i = aimag(tau(2,1))
      t22r = real(tau(2,2))
      t22i = aimag(tau(2,2))

!     radiance(1,1) = radiance(1,1) + scalar * (tau(1,1)*conjg(tau(1,1)) + &
!                                            &  tau(1,2)*conjg(tau(1,2)) )
      radiance(1,1) = radiance(1,1) + scalar * &
        & ( t11r * t11r + t11i * t11i + t12r * t12r + t12i * t12i )

!     radiance(1,2) = radiance(1,2) + scalar * (tau(1,1)*conjg(tau(2,1)) + &
!                                            &  tau(1,2)*conjg(tau(2,2)) )
      radiance(1,2) = radiance(1,2) + scalar * &
        & cmplx( t11r * t21r + t11i * t21i + t12r * t22r + t12i * t22i, &
        &       -t11r * t21i + t11i * t21r - t12r * t22i + t12i * t22r )

      radiance(2,1) = conjg(radiance(1,2))

!     radiance(2,2) = radiance(2,2) + scalar * (tau(2,1)*conjg(tau(2,1)) + &
!                                            &  tau(2,2)*conjg(tau(2,2)) )
      radiance(2,2) = radiance(2,2) + scalar * &
        & ( t21r * t21r + t21i * t21i + t22r * t22r + t22i * t22i )

    end subroutine Updaterad

  end subroutine Mcrt

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MCRT_m

! $Log$
! Revision 1.1.2.2  2003/02/15 00:28:42  vsnyder
! Don't exp(del_opcty) here, it's done by the caller
!
! Revision 1.1.2.1  2003/02/14 03:54:10  vsnyder
! Initial commit
!
