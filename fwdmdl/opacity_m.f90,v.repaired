! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Opacity_M

!---------------------------------------------------------------------
! This module computes the incremental opacity matrices.


  implicit none
  private
  public :: Opacity

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Opacity ( CT, STCP, STSP, Alpha_path, Del_opcty )

    use MLSCommon, only: Rk => Rp

! Arguments

! Theta is the angle between the line of sight and magnetic field vectors.

! Phi is the angle between the plane defined by the line of sight and the
! magnetic field vector, and the "instrument field of view plane polarized"
! (IFOVPP) X axis.

    real(rk), intent(in) ::    CT(:)     ! Cos(Theta)
    real(rk), intent(in) ::    STCP(:)   ! Sin(Theta) Cos(Phi)
    real(rk), intent(in) ::    STSP(:)   ! Sin(Theta) Sin(Phi)

    ! Alpha_path(-1,:) is Alpha_sigma_m, Alpha_path(0,:) is Alpha_pi,
    ! Alpha_path(+1,:) is Alpha_sigma_p
    complex(rk), intent(in) ::  Alpha_path(-1:,:)

    complex(rk), intent(out) :: Del_opcty(:,:,:) ! (2,2,:)

! Local variables

    integer :: h_i
    real(rk) ::    rhopi1, rhopi2, rhosi1, rhosi2, rhorrx
    complex(rk) :: del_m, del_p
!   complex(rk), parameter :: im_one = (0.0, 1.0)

!{ Compute elements of the $\rho$ matrices.
!  Let $k = \cos\theta$, $s = \sin\theta \sin\phi$ and $c = \sin\theta \cos\phi$.
!  \begin{equation*}\begin{split}
!  {\bf\rho}_{-1} =& \left [ \begin{array}{cc}
!                      \cos^2 \phi + \sin^2 \phi \cos^2 \theta &
!                       -\sin\phi \cos\phi \sin^2 \theta - i \cos \theta \\
!                       -\sin\phi \cos\phi \sin^2 \theta + i \cos \theta &
!                       \sin^2 \phi + \cos^2 \phi \cos^2 \theta
!                     \end{array} \right ] \\
!                 =& \left [ \begin{array}{cc}
!                       1 - s^2   & -sc  - i k \\
!                       -sc + i k & 1 - c^2
!                     \end{array} \right ]
!                 =  \left [ \begin{array}{cc}
!                       \rho_{\sigma_1} & -\rho_r - i k \\
!                       -\rho_r + i k   & \rho_{\sigma_2}
!                     \end{array} \right ] = {\bf\rho}^H_{-1}\\
!  {\bf\rho}_0 =& \sin^2 \theta \left [ \begin{array}{cc}
!                      \sin^2 \phi & \sin\phi \cos\phi \\
!                      \sin\phi \cos\phi & \cos^2 \phi
!                 \end{array} \right ] =
!                 \left [ \begin{array}{cc}
!                      s^2 & sc \\
!                      sc  & c^2
!                 \end{array} \right ] = {\bf\rho}^T_0 =
!                 \left [ \begin{array}{cc}
!                      \rho_{\pi_1} & \rho_r \\
!                      \rho_r       & \rho_{\pi_2}
!                 \end{array} \right ] \\
!  {\bf\rho}_{+1} =& \left [ \begin{array}{cc}
!                      \cos^2 \phi + \sin^2 \phi \cos^2 \theta &
!                       -\sin\phi \cos\phi \sin^2 \theta + i \cos \theta \\
!                       -\sin\phi \cos\phi \sin^2 \theta - i \cos \theta &
!                       \sin^2 \phi + \cos^2 \phi \cos^2 \theta
!                     \end{array} \right ] = {\bf\rho}^T_{-1}
!%                    = \bar{\bf \rho}_{-1}}
!                     = {\bf \rho}^*_{-1}
!  \end{split}\end{equation*}

! cover all heights

    do h_i = 1, size(alpha_path,2)

      if ( all(alpha_path(:,h_i) == (0.0_rk,0.0_rk) ) ) then
        ! This is the most common case for derivatives.
        del_opcty(:,:,h_i) = 0
        cycle
      end if

! Compute elements of the rho matrices
      rhopi1 = stsp(h_i)**2            ! \rho_0 (1,1)
      rhopi2 = stcp(h_i)**2            ! \rho_0 (2,2)
      rhosi1 = 1.0_rk - rhopi1         ! \rho_{+-1} (1,1)
      rhosi2 = 1.0_rk - rhopi2         ! \rho_{+-1} (2,2)
      rhorrx = stsp(h_i) * stcp(h_i)   ! used in real part of all off-diagonal elements

! Compute the total incremental opacity for the segment
! Sum over the coefficients C and the species S

!{Let {\tt Alpha\_path(-1:+1,:)} be $\Delta_{\sigma_-}$, $\Delta_{\pi}$,
! and $\Delta_{\sigma_+}$ respectively.
!
!  Let {\tt del\_opcty(:,:,h\_i)} be $\delta_{h_i}$.  
!  Let $\Delta_+ = \Delta_{\sigma_+} + \Delta_{\sigma_-}$ and
!      $\Delta_- = \Delta_{\sigma_+} - \Delta_{\sigma_-}$.  Then
!  \begin{equation*}
!  \delta_{h_i} = \left[ \begin{array}{cc}
!                  \rho_{\pi_1} \Delta_\pi + \rho_{\sigma_1} \Delta_+ &
!                  \rho_r ( \Delta_\pi - \Delta_+ ) + i \cos \theta\, \Delta_- \\
!                  \rho_r ( \Delta_\pi - \Delta_+ ) - i \cos \theta\, \Delta_- &
!                  \rho_{\pi_2} \Delta_\pi + \rho_{\sigma_2} \Delta_+
!                 \end{array} \right]
!               = \rho_{-1} \Delta_{\sigma_-} + \rho_0 \Delta_\pi +
!                 \rho_{+1} \Delta_{\sigma_+}
!  \end{equation*}

      del_m = ct(h_i) * ( alpha_path(+1,h_i) - alpha_path(-1,h_i) )
      del_p = alpha_path(+1,h_i) + alpha_path(-1,h_i)

      del_opcty(1,1,h_i) = rhopi1 * alpha_path(0,h_i) + rhosi1 * del_p
!     del_opcty(1,2,h_i) = rhorrx * (alpha_path(0,h_i) - del_p) &
!       &                  - im_one * del_m
!     del_opcty(2,1,h_i) = rhorrx * (alpha_path(0,h_i) - del_p) &
!       &                  + im_one * del_m
      del_opcty(1,2,h_i) = rhorrx * (alpha_path(0,h_i) - del_p)
      del_opcty(2,1,h_i) = del_opcty(1,2,h_i) - cmplx(-aimag(del_m),real(del_m), kind=rk)
      del_opcty(1,2,h_i) = del_opcty(1,2,h_i) + cmplx(-aimag(del_m),real(del_m), kind=rk)
      del_opcty(2,2,h_i) = rhopi2 * alpha_path(0,h_i) + rhosi2 * del_p

    end do

  end subroutine Opacity

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Opacity_M

! $Log$
! Revision 2.6  2017/08/09 20:20:12  vsnyder
! Set del_opcty to zero where alpha_path is zero because testing is cheaper
! than evaluating everything to get zero.  This is the usual case for
! derivatives along the path, which depend upon only two state-vector
! elements.
!
! Revision 2.5  2010/02/04 23:09:53  vsnyder
! Use kind= in CMPLX
!
! Revision 2.4  2009/06/23 18:26:11  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.3  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.2  2003/05/05 23:00:26  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.1.2.1  2003/03/05 03:30:32  vsnyder
! Calculate the off-diagonal elements in the correct order
!
! Revision 2.1  2003/02/06 22:38:22  vsnyder
! Initial commit
!
