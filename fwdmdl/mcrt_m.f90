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
  subroutine Mcrt ( T_script, Sqrt_earth_ref, Del_tau, P_Stop, Prod, Tau, Radiance )

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
!  {\bf P}_i = & \prod_{k=2}^{i} {\bf E}_k \text{ and } {\bf E} \text{ is
!    {\tt Del\_tau} (input), which is}\\
!  {\bf E}_i = & \exp \left( - \int_{s_i}^{s_{i-1}} {\bf G} (s^\prime)
!                \text{d} s^\prime \right )
! \end{split}
! \end{equation*}
!
! $\mathbf{E}_1$ is the incremental transmissivity from the spacecraft to
! the top of the atmosphere, which is obviously zero.  That's why we start
! with $\mathbf{E}_2$.

    use MLSCommon, only: Rk => Rp

    real(rk), intent(in) :: T_script(:)          ! Called Delta B in some notes
    real(rk), intent(in) :: Sqrt_earth_ref
    complex(rk), intent(in) :: Del_tau(:,:,:)    ! 2 x 2 x path = exp(del_opcty).
                                                 ! Called P in some notes and
                                                 ! E in other notes.  It's E here.
    integer, intent(in) :: P_Stop                ! Stop here
    complex(rk), intent(out) :: Prod(:,:,:)      ! 2 x 2 x path.  Called P in some
                                                 ! notes.  Prod(:,:,I) is the
                                                 ! product of del_tau(:,:,:I).
    complex(rk), intent(out) :: Tau(:,:,:)       ! 2 x 2 x path.  
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
    tau(:,:,1) = ident
    ! radiance = t_script(1) * tau(1)
    radiance(1,1) = t_script(1)
    radiance(1,2) = 0.0_rk
    radiance(2,1) = 0.0_rk
    radiance(2,2) = t_script(1)

! Proceed with first segment integration

! We use 1:2 instead of : for the first two dimensions of the arrays of
! 2x2 matrices, in the hope that some compiler someday may be able to
! use the information for optimization.

    i_tan = n_path / 2
    do i = 2, min(i_tan,p_stop)
      prod(1:2,1:2,i) = matmul ( prod(1:2,1:2,i-1), del_tau(1:2,1:2,i) )
      call updaterad ( radiance, t_script(i), prod(1:2,1:2,i), tau(1:2,1:2,i) )
    end do

! Tangent point (or Earth intersecting) layer.  If it's not
! an Earth intersecting layer, sqrt_earth_ref will be 1.0.

    if ( p_stop <= i_tan ) return
    prod(1:2,1:2,i_tan+1) = sqrt_earth_ref * prod(1:2,1:2,i_tan)
    call updaterad ( &
      & radiance, t_script(i_tan+1), prod(1:2,1:2,i_tan+1), tau(1:2,1:2,i_tan+1) )

! Proceed with third segment integration, which includes the
! space radiance contribution.

    do i = i_tan+2, min(n_path,p_stop)
      prod(1:2,1:2,i) = matmul ( prod(1:2,1:2,i-1) , del_tau(1:2,1:2,i-1) )
      call updaterad ( radiance, t_script(i), prod(1:2,1:2,i), tau(1:2,1:2,i) )  
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

!     Specifying the ranges hopefully prevents accessing the dope vectors
      radiance(1:2,1:2) = radiance(1:2,1:2) + scalar * tau(1:2,1:2)

    end subroutine Updaterad

  end subroutine Mcrt

! -----------------------------------------------------  Mcrt_Der  -----
  subroutine Mcrt_Der ( T_script, Sqrt_earth_ref, E, D_E, Prod, Tau, P_Stop, &
    & D_Radiance, D_T_script )

!       Magnetic    Condensed    Radiative    Transfer    Derivative

! Compute d(Radiance)/d(Whatever) given T_script (aka Delta_B),
! D_T_script (aka D_Delta_B), D_E (d(Deltau)/d(Whatever)), Prod
! (product of E's), and Tau (P * P^*).

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
  !  = & \sum_{k=2}^{i} \left [ {\bf E}_2 \dots {\bf E}_{k-1}
  !       \frac{\partial {\bf E}_k}{\partial x} {\bf E}_{k+1} \dots
  !       {\bf E}_i {\bf P}_i^\dagger +
  !       {\bf P}_i {\bf E}_i^\dagger \dots {\bf E}_{k+1}
  !       \frac{\partial {\bf E}_k^\dagger}{\partial x}
  !       {\bf E}_{k-1}^\dagger \dots {\bf E}_2^\dagger \right ]
  !\\
  !  = & \sum_{k=2}^{i} \left [ {\bf P}_{k-1}
  !       \frac{\partial {\bf E}_k}{\partial x} {\bf P}_k^{-1} {\bf \tau}_i
  !     + {\bf \tau}_i^\dagger {\bf P}_k^{-\dagger}
  !       \frac{\partial {\bf E}_k^\dagger}{\partial x} {\bf P}_{k-1}^\dagger \right ]
  ! \end{split}
  ! \end{equation*}
  ! Notice that $\tau = \tau^\dagger$ by construction (see {\tt MCRT}).  Then
  ! \begin{equation*}
  ! \begin{split}
  !  \frac{\partial {\bf I}}{\partial x} =& \sum_{i=1}^n
  !    \mathcal{Q}_i {\bf \tau}_i + {\bf \tau}_i^\dagger \mathcal{Q}_i^\dagger
  !    \text{, where}
  !\\
  !  \mathcal{Q}_i =& \frac12 \frac{\partial \Delta B_i}{\partial x} \mathbf{1} +
  !    \Delta B_i \mathcal{W}_i \text{ and}
  !\\
  !  \mathcal{W}_i =& \sum_{k=2}^i {\bf P}_{k-1}
  !                     \frac{\partial {\bf E}_k}{\partial x} {\bf P}_k^{-1}
  ! \end{split}
  ! \end{equation*}
  ! {\bf 1} is the identity matrix (we're using ${\bf I}$ for the radiance
  ! matrix).

    use MLSCommon, only: Rk => Rp

    ! SVE == state_vector_elements
    real(rk), intent(in) :: T_script(:)     ! sve.  Called Delta B above
    real(rk), intent(in) :: Sqrt_earth_ref
    complex(rk), intent(in) :: E(:,:,:)     ! 2 x 2 x path.  Deltau
    complex(rk), intent(in) :: D_E(:,:,:,:) ! 2 x 2 x path x sve.
                                            ! D (deltau) / D (Whatever)
    complex(rk), intent(in) :: Prod(:,:,:)  ! 2 x 2 x path.  Called P above.
    complex(rk), intent(in) :: Tau(:,:,:)   ! 2 x 2 x path. Matmul(Prod,conjg(Prod)).
    integer, intent(in) :: P_Stop           ! Where to stop on the path
    complex(rk), intent(out) :: D_Radiance(:,:,:) ! 2 x 2 x sve
    real(rk), intent(in),optional :: D_T_script(:,:) ! path x sve. a.k.a D Delta B
    ! T script or Delta B depends only on temperature and frequency, so it's
    ! only needed for temperature derivatives.

    integer :: I_P, I_pp, I_Sv              ! Path, State vector indices
    integer :: I_Tan                        ! Tangent point
    complex(rk) :: PDET                     ! Det(P_{k+1})
    complex(rk) :: PINV(2,2,size(t_script)) ! P_{k+1}^{-1}
!   real(rk), parameter :: PTOL = (radix(0.0_rk)+0.0_rk) ** minexponent(0.0_rk)
      ! PTOL = sqrt(tiny(0.0_rk))
    real(rk), parameter :: PTOL = epsilon(0.0_rk) ! defines "PDET is too small"
    complex(rk) :: Q(2,2)
    complex(rk) :: Q_Tau(2,2)               ! Q x Tau
    integer :: NP                           ! Index of first PINV not computed
    complex(rk) :: W(2,2)

    i_tan = size(t_script) / 2

    ! Compute the P_{k+1}^{-1} matrices, but only up to where det(P_{k+1})
    ! becomes small.  We can get away with this for two reasons.  First, we
    ! use P_{k+1}^{-1} \tau_i to get \prod{j=k+1}^{i-1} E_j.  Second,
    ! because |P| <= 1.0, a product of P's is smaller than any of them.

    do np = 2, p_stop
      pdet = prod(1,1,np)*prod(2,2,np) - prod(1,2,np)*prod(2,1,np)
      if ( abs(pdet) <= ptol ) & ! getting too small
    exit
      pinv(:,:,np) = reshape( (/ prod(2,2,np), -prod(2,1,np), &
                &               -prod(1,2,np),  prod(1,1,np) /), (/2,2/) ) / &
                &    pdet
    end do
    ! NP-1 is the index of the last computed PINV

    d_radiance = 0.0_rk

    if ( present(d_t_script) ) then
      do i_sv = 1, size(d_radiance,3)      ! state vector elements
        d_radiance(1,1,i_sv) = d_t_script(1,i_sv) ! W(1)=0 and Tau(1)=Ident
        d_radiance(2,2,i_sv) = d_t_script(1,i_sv)
        w = 0.0_rk
        i_pp = 2
        do i_p = 2, np - 1                 ! path elements
          ! w = w + P_k DE_k P_{k+1}^{-1}
          if ( i_p /= i_tan + 1 ) then ! Skip adding to W at the tangent point
            w = w + matmul ( matmul ( &
            &                prod(1:2,1:2,i_p-1), d_e(1:2,1:2,i_pp,i_sv) ), &
            &       pinv(:,:,i_p) )
          else
            i_pp = i_p - 1
          end if
          q(1,1) = 0.5 * d_t_script(i_p,i_sv) + t_script(i_p) * w(1,1)
          q(1,2) =                              t_script(i_p) * w(1,2)
          q(2,1) =                              t_script(i_p) * w(2,1)
          q(2,2) = 0.5 * d_t_script(i_p,i_sv) + t_script(i_p) * w(2,2)
          q_tau = matmul(q,tau(1:2,1:2,i_p))
          d_radiance(1:2,1:2,i_sv) = d_radiance(1:2,1:2,i_sv) + &
                                   &   q_tau + conjg(transpose(q_tau))
          i_pp = i_pp + 1
        end do ! i_p
      end do ! i_sv
    else
      do i_sv = 1, size(d_radiance,3)      ! state vector elements
        w = 0.0_rk
        i_pp = 2
        do i_p = 2, np - 1                 ! path elements
          ! w = w + P_k DE_k P_{k+1}^{-1}
          if ( i_p /= i_tan + 1 ) then ! Skip adding to W at the tangent point
            w = w + matmul ( matmul ( &
            &                prod(1:2,1:2,i_p-1), d_e(1:2,1:2,i_pp,i_sv) ), &
            &       pinv(:,:,i_p) )
          else
            i_pp = i_p - 1
          end if
          q(1,1) = t_script(i_p) * w(1,1)
          q(1,2) = t_script(i_p) * w(1,2)
          q(2,1) = t_script(i_p) * w(2,1)
          q(2,2) = t_script(i_p) * w(2,2)
          q_tau = matmul(q,tau(1:2,1:2,i_p))
          d_radiance(1:2,1:2,i_sv) = d_radiance(1:2,1:2,i_sv) + &
                                   &   q_tau + conjg(transpose(q_tau))
          i_pp = i_pp + 1
        end do ! i_p
      end do ! i_sv
    end if

  end subroutine Mcrt_Der

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MCRT_m

! $Log$
! Revision 2.10  2003/08/15 20:29:26  vsnyder
! Implement polarized VMR derivatives
!
! Revision 2.9  2003/08/14 19:37:55  vsnyder
! Futzing with comments
!
! Revision 2.8  2003/06/09 20:52:37  vsnyder
! More work on polarized derivatives
!
! Revision 2.7  2003/05/28 01:25:02  vsnyder
! Hopefully squashed some more bugs in polarized derivatives
!
! Revision 2.6  2003/05/27 22:31:12  vsnyder
! More work on polarized derivatives
!
! Revision 2.5  2003/05/24 02:26:18  vsnyder
! More work on polarized temperature derivatives
!
! Revision 2.4  2003/05/15 03:29:44  vsnyder
! Implement polarized model's temperature derivatives
!
! Revision 2.3  2003/05/10 01:05:07  livesey
! Typo!
!
! Revision 2.2  2003/05/09 19:27:02  vsnyder
! Initial stuff for temperature derivatives
!
! Revision 2.1  2003/05/05 23:00:26  livesey
! Merged in feb03 newfwm branch
!
! Revision 1.1.2.2  2003/02/15 00:28:42  vsnyder
! Don't exp(del_opcty) here, it's done by the caller
!
! Revision 1.1.2.1  2003/02/14 03:54:10  vsnyder
! Initial commit
!
