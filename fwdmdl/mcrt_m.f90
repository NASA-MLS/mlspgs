! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MCRT_m

!       Magnetic    Condensed    Radiative    Transfer

  implicit NONE
  private
  public :: MCRT, MCRT_Der

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
  subroutine Mcrt ( T_script, Sqrt_earth_ref, Del_tau, Prod, Tau, Radiance )

!       Magnetic    Condensed    Radiative    Transfer

! Compute the polarized radiative transfer using a condensed algorithm.
! Enter t_script(*) in a quantity that is linear in power.
! where T is intensity of a black body radiator having temperature t
! at frequency f (usually entered in units of kelvin being derived from
! T_{i} = h v / [ k ( exp{h v / (k t_{i}) } - 1)]).
! sqrt_earth_ref is the square root of Earth's reflectivity.

!{\begin{equation*}
! \begin{split}
!  {\bf I} =& \sum_{i=1}^n \Delta B_i {\bf \tau}_i,
!  \text{ where {\bf I} is {\tt Radiance} (output), }
!        \Delta B \text{ is {\tt T\_script} (input)},\\
!  {\bf \tau}_i =& {\bf P}_i {\bf P}_i^\dagger,~
!    {\bf \tau} \text{ is {\tt Tau} (output), {\bf P} is {\tt Prod} (output),}\\
!  {\bf P}_i = & \prod_{k=1}^{i-1} {\bf E}_k \text{ and } {\bf E} \text{ is
!    {\tt Del\_tau} (input), which is}\\
!  {\bf E}_i = & \exp \left( - \int_{s_i}^{s_{i-1}} {\bf G} (s^\prime)
!                \text{d} s^\prime \right )
! \end{split}
! \end{equation*}

    use MLSCommon, only: Rk => Rp

    real(rk), intent(in) :: T_script(:)          ! Called Delta B in some notes
    real(rk), intent(in) :: Sqrt_earth_ref
    complex(rk), intent(in) :: Del_tau(:,:,:)    ! 2 x 2 x : = exp(del_opcty).
                                                 ! Called P in some notes and
                                                 ! E in other notes.
    complex(rk), intent(out) :: Prod(:,:,:)      ! 2 x 2 x :.  Called P in some
                                                 ! notes.  Prod(:,:,I) is the
                                                 ! product of del_tau(:,:,:I).
    complex(rk), intent(out) :: Tau(:,:,:)       ! 2 x 2 x :.  
                                                 ! Matmul(Prod,conjg(Prod)).
    complex(rk), intent(out) :: Radiance(2,2)    ! Sum(Delta B_i * matmul (
                                                 !  Prod(:,:,i), conjg(
                                                 !   transpose(Prod(:,:,i)))))

    integer :: i, i_tan, n_path

    complex(rk), parameter :: Ident(2,2) = reshape( (/ 1.0, 0.0, &
                                                     & 0.0, 1.0 /), (/2,2/) )

    n_path = size(t_script)

! Initialize segment 1 calculation

    prod(:,:,1) = ident
    ! radiance = t_script(1) * ident
    radiance(1,1) = t_script(1)
    radiance(1,2) = 0.0
    radiance(2,1) = 0.0
    radiance(2,2) = t_script(1)

! Proceed with first segment integration

    i_tan = n_path / 2
    do i = 2, i_tan
      prod(:,:,i) = matmul ( prod(:,:,i-1), del_tau(:,:,i) )
      call updaterad ( radiance, t_script(i), prod(:,:,i), tau(:,:,i) )
    end do

! Tangent point (or Earth intersecting) layer.  If it's not
! an Earth intersecting layer, sqrt_earth_ref will be 1.0.

    prod(:,:,i_tan+1) = sqrt_earth_ref * prod(:,:,i_tan)
    call updaterad ( &
      & radiance, t_script(i_tan+1), prod(:,:,i_tan+1), tau(:,:,i_tan+1) )

! Proceed with third segment integration, which includes the
! space radiance contribution.

    do i = i_tan+2, n_path
      prod(:,:,i) = matmul ( prod(:,:,i-1) , del_tau(:,:,i-1) )
      call updaterad ( radiance, t_script(i), prod(:,:,i), tau(:,:,i) )  
    end do                                                       

  contains
  ! ..................................................  Updaterad  .....

    subroutine Updaterad ( Radiance, Scalar, Prod, Tau )

! Update the radiance from the prod:
! Radiance = Radiance + Scalar * Prod * conjg(transpose(Prod))
! We know that Radiance is Hermitian, so the diagonal elements are
! real and the off-diagonal elements are conjugates.  Exploiting
! this saves roughly half of the work.

      complex(rk), intent(inout) :: Radiance(2,2)
      real(rk), intent(in) :: Scalar
      complex(rk), intent(in) :: Prod(:,:) ! Really (2,2), but this may
                                           ! avoid a copy-in/copy-out
      complex(rk), intent(out) :: Tau(:,:) ! Really (2,2), but this may
                                           ! avoid a copy-in/copy-out

      real(rk) t11r, t11i, t12r, t12i, t21r, t21i, t22r, t22i

!     tau = matmul(prod,conjg(transpose(prod)))

      t11r = real(prod(1,1))
      t11i = aimag(prod(1,1))
      t12r = real(prod(1,2))
      t12i = aimag(prod(1,2))
      t21r = real(prod(2,1))
      t21i = aimag(prod(2,1))
      t22r = real(prod(2,2))
      t22i = aimag(prod(2,2))

!     tau(1,1) = (prod(1,1)*conjg(prod(1,1)) + prod(1,2)*conjg(prod(1,2)) )
      tau(1,1) =  ( t11r * t11r + t11i * t11i + t12r * t12r + t12i * t12i )

!     tau(1,2) = (prod(1,1)*conjg(prod(2,1)) + prod(1,2)*conjg(prod(2,2)) )
      tau(1,2) = cmplx( t11r * t21r + t11i * t21i + t12r * t22r + t12i * t22i, &
               &       -t11r * t21i + t11i * t21r - t12r * t22i + t12i * t22r )

      tau(2,1) = conjg(tau(1,2))

!     tau(2,2) = (prod(2,1)*conjg(prod(2,1)) + prod(2,2)*conjg(prod(2,2)) )
      tau(2,2) = ( t21r * t21r + t21i * t21i + t22r * t22r + t22i * t22i )

!     Specifying the rangees hopefully prevents accessing the dope vectors
      radiance(1:2,1:2) = radiance(1:2,1:2) + scalar * tau(1:2,1:2)

    end subroutine Updaterad

  end subroutine Mcrt

! -----------------------------------------------------  Mcrt_Der  -----
  subroutine Mcrt_Der ( T_script, D_T_script, D_E, Prod, Tau, D_Radiance )

!       Magnetic    Condensed    Radiative    Transfer    Derivative

! Compute the derivative of Radiance given T_script (aka Delta_B),
! D_T_script (aka D_Delta_B), D_E, Prod (product of E's), and Tau (P * P^*).

  !{\begin{equation*}
  ! \begin{split}
  !  \frac{\partial \bf I}{\partial x} = &\sum_{i=1}^n
  !   \left [ \frac{\partial {\bf \tau}_i}{\partial x} \Delta B_i +
  !           {\bf \tau}_i \frac{\partial \Delta B_i}{\partial x}
  !   \right ],\text{ where}
  !\\
  !  \frac{\partial {\bf \tau}_i}{\partial x} = &
  !   \frac{\partial {\bf P}_i {\bf P}_i^\dagger}{\partial x}
  !\\
  !  = & \sum_{k=1}^{i-1} \left [ {\bf E}_i \dots E_{k-1}
  !       \frac{\partial {\bf E}_k}{\partial x} {\bf E}_{k+1} \dots
  !       {\bf E}_{i-1} {\bf P}_i^\dagger +
  !       {\bf P}_i {\bf E}_{i-1}^\dagger \dots {\bf E}_{k+1}
  !       \frac{\partial {\bf E}_k^\dagger}{\partial x}
  !       {\bf E}_{k-1}^\dagger \dots {\bf E}_1^\dagger \right ]
  !\\
  !  = & \sum_{k=1}^{i-1} \left [ {\bf P}_k
  !       \frac{\partial {\bf E}_k}{\partial x} {\bf P}_{k+1}^{-1} {\bf \tau}_i
  !     + {\bf \tau}_i^\dagger {\bf P}_{k+1}^{-\dagger}
  !       \frac{\partial {\bf E}_k^\dagger}{\partial x} {\bf P}_k^\dagger \right ]
  ! \end{split}
  ! \end{equation*}
  ! Notice that $\tau = \tau^\dagger$ by construction (see {\tt MCRT}).  Then
  ! \begin{equation*}
  ! \begin{split}
  !  \frac{\partial {\bf I}}{\partial x} =& \sum_{i=1}^n
  !    \mathcal{Q}_i {\bf \tau}_i + {\bf \tau}_i^\dagger \mathcal{Q}_i^\dagger
  !    \text{, where}
  !\\
  !  \mathcal{Q}_i =& \frac12 \frac{\partial \Delta B_i}{\partial x} {\bf 1} +
  !    \Delta B_i \mathcal{W}_i \text{ and}
  !\\
  !  \mathcal{W}_i =& \sum_{k=1}^{i-1} {\bf P}_k
  !                     \frac{\partial {\bf E}_k}{\partial x} {\bf P}_{k+1}^{-1}
  ! \end{split}
  ! \end{equation*}
  ! {\bf 1} is the identity matrix (we're using ${\bf I}$ for the radiance
  ! matrix).

    use dExDt_M, only: dExDt               ! D (exp(2x2 matrix))
    use MLSCommon, only: Rk => Rp

    real(rk), intent(in) :: T_script(:)    ! Called Delta B above
    real(rk), intent(in) :: D_T_script(:)  ! 0.5 * D Delta B
    complex(rk), intent(in) :: D_E(:,:,:)  ! 2 x 2 x :.  D (deltau)
    complex(rk), intent(in) :: Prod(:,:,:) ! 2 x 2 x :.  Called P above.
    complex(rk), intent(in) :: Tau(:,:,:)  ! 2 x 2 x :. Matmul(Prod,conjg(Prod)).
    complex(rk), intent(out) :: D_Radiance(2,2)

    integer :: i, k
    complex(rk) :: PINV(2,2)               ! P_{k+1}^{-1}
    complex(rk) :: Q(2,2)
    complex(rk) :: Q_Tau(2,2)              ! Q x Tau
    complex(rk) :: W(2,2)

    d_radiance = 0.0_rk
    do i = 1, size(t_script)
      w = 0.0_rk
      do k = 1, i-1
        ! w = w + P_k DE_k P_{k+1}^{-1}
        pinv = reshape( (/ prod(2,2,k+1), -prod(2,1,k+1), &
          &             -prod(1,2,k+1), prod(1,1,k+1) /), (/2,2/) ) / &
          &   ( prod(1,1,k+1)*prod(2,2,k+1) - prod(1,2,k+1)*prod(2,1,k+1) )
        w = w + matmul ( matmul ( prod(1:2,1:2,k), d_e(1:2,1:2,k) ), pinv )
          &  
      end do
      q(1,1) = d_t_script(k) + t_script(k) * w(1,1)
      q(1,2) =                 t_script(k) * w(1,2)
      q(2,1) =                 t_script(k) * w(2,2)
      q(2,2) = d_t_script(k) + t_script(k) * w(2,2)
      q_tau = matmul(q,tau(1:2,1:2,i))
      d_radiance = d_radiance + q_tau + conjg(transpose(q_tau))
    end do

  end subroutine Mcrt_Der

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MCRT_m

! $Log$
! Revision 2.1  2003/05/05 23:00:26  livesey
! Merged in feb03 newfwm branch
!
! Revision 1.1.2.2  2003/02/15 00:28:42  vsnyder
! Don't exp(del_opcty) here, it's done by the caller
!
! Revision 1.1.2.1  2003/02/14 03:54:10  vsnyder
! Initial commit
!
