! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module dExDt_M

!---------------------------------------------------------------------
! This module computes the derivative of an exponential of a matrix.
! This isn't just the exponential of the matrix, as it would be for a
! scalar, because the argument and its derivative don't necessarily
! commute.  For example, the third term in the Taylor series for
! exp(z(x))) is z^2/2, so its derivative is ( z z' + z' z ) / 2, not z z'.

  implicit NONE
  private
  public :: dExDt

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    &  "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine dExDt ( A, dA, dEx )

    use CS_GetEv_M, only: CS_GetEvSD
    use MLSCommon, only: RK => Rp
    use SINHZ_Z_M, only: D_SINHZ_Z3, SINHZ_Z

    complex(rk), intent(in) :: A(2,2)
    complex(rk), intent(in) :: dA(2,2)
    complex(rk), intent(out) :: dEx(2,2)

    complex(rk) :: ACOF       ! Coefficients for A, A'
    complex(rk) :: dEv(2)     ! Sum = s', then difference**2 = h'
    complex(rk) :: Ev(2)      ! Sum = s,  then difference**2 = h
    complex(rk) :: D          ! 0.5(Difference between eigenvalues, z1 - z2)
    complex(rk) :: D_SINHD_D  ! exp(s) [ (sinh(d)/d)' / d ] =
    !                         ! exp(s) [ (d cosh d - sinh d) / d**3 ]
    complex(rk) :: DH2        ! h' / 2
    complex(rk) :: DIAG       ! Multiple of I to add to diagonal
    complex(rk) :: ES         ! exp(s)
    complex(rk) :: SINHD_D    ! exp(s) sinh(d)/d

!{\parskip 2pt
!  By applying the Hamilton-Cayley theorem (``A matrix satisfies its
!  characteristic polynomial'') to the Taylor series for a function of a
!  matrix, one can develop Sylvester's identity.  In the case of a function
!  of a $2 \times 2$ matrix A whose eigenvalues are $z_1$ and $z_2$ we have
!
!    $F({\bf A}) = \frac{F(z_1)}
!                       {z_1-z_2} ({\bf A} - z_2 {\bf I}) +
!                  \frac{F(z_2)}
!                       {z_2-z_1} ({\bf A} - z_1 {\bf I})$
!%
!  where {\bf I} is the $2 \times 2$ identity matrix,
!
!                              Or:
!
!    $F({\bf A}) = \frac{F(z_1) ({\bf A} - z_2 {\bf I}) - 
!                        F(z_2) ({\bf A} - z_1 {\bf I})}
!                       {z_1 - z_2}$
!
!  For any function $F()$ (in this case: $\exp()$ )
!
!  If you assume that {\bf A} is a function of some parameter $p$, and work
!  through $\frac{\text{d} F({\bf A})}{\text{d}p}$ you will eventually
!  find $\frac{\text{d}z_1}{\text{d}p}$ and $\frac{\text{d}z_1}{\text{d}p}$.
!  These derivatives approach infinity as the eigenvalues approach each other.
!  Using
!  $s = \frac12 ( z_1 + z_2 ) = \frac12 \text{tr}(\mathbf{A})$
!  and
!  $d = \frac12 ( z_1 - z_2 ) = \sqrt { s^2 - \det(\mathbf{A})}
!   = \sqrt h$,
!  this can be rewritten as
!%
!  $\exp({\bf A}) =
!    e^s \left[ \frac{\sinh d}d ( \mathbf{A} - s\, \mathbf{I}) +
!               \cosh d\, \mathbf{I} \right ]$.  After much work you will find
!
!  \begin{equation}\begin{split}\label{dexp}
!  \frac{\text{d} \exp(\mathbf{A})}{\text{d}p} =\;
!   e^s & \left\{ \frac{\sinh d}d
!         \left[ s^\prime \mathbf{A} + \mathbf{A}^\prime +
!                (\frac{h^\prime}2 - s^\prime s) \mathbf{I} \right]
!                 + \right.\\
!       & \left. \frac{d \cosh d - \sinh d}{d^3}
!           \left[ \frac{h^\prime}{2} \mathbf{A} +
!                  ( h s^\prime - \frac{h^\prime}2 s )\,\mathbf{I} ) \right]
!       \right \} \;.
!  \end{split}\end{equation}
!
!  Collecting terms, we have
! \begin{equation}\begin{split}
! \frac{\text{d} \exp(\mathbf{A})}{\text{d}p} =\;&
!   e^s \left [ s^\prime \frac{\sinh d}d +
!           \frac{h^\prime}{2} \frac{d \cosh d - \sinh d}{d^3} \right ] \mathbf{A}
!   + e^s \frac{\sinh d}d \mathbf{A}^\prime + \\
! &   e^s \left[ \frac{\sinh d}d
!                \left( \frac{h^\prime}{2} - s^\prime s \right) +
!                \frac{d \cosh d - \sinh d}{d^3}
!                \left( h s^\prime - \frac{h^\prime}{2} s \right) \right]
!         \mathbf{I}
! \end{split}\end{equation}
!
!  As the eigenvalues coalesce, no cancellations occur, and no infinities
!  arise if  the elements of $\mathbf{A}$ and $\mathbf{A}^\prime$ are
!  finite.  The $h^\prime$ and $s^\prime$ terms are clearly well behaved as
!  the eigenvalues coalesce.  To see that the functions of $d$ are well
!  behaved, write
!
!  $\frac{\sinh\,d}d = \sum_{k=0}^{\infty}
!  \frac{d^{2k}}{(2k+1)!}$ and
!  %
!  $\frac{d \cosh d - \sinh d}{d^3} =
!  \sum_{k=0}^{\infty} \frac{d^{2k}}{2^k k! (2k+3)!!}$, where $(2k+3)!! =
!  3 \cdot 5 \cdot 7 \cdot \cdot \cdot 2k+3$.
!
!  We can write $\frac{\sinh d}d$ as $e^{-d}\,\frac{e^{2d}-1}{2d}$ and use
!  the same software to evaluate this expression as in {\tt cs\_expmat\_m}.
!  The series for $\frac{d \cosh d - \sinh d}{d^3}$ converges extremely
!  rapidly for small $d$ (see {\tt sinhz\_z\_m}).

    call CS_GetEvSD ( A, Ev, dA, dEv )   ! Get the eigenvalues and their
                                         ! derivatives

    es = exp(Ev(1))                      ! exp(s)
    d = sqrt(Ev(2))                      ! z1 - z2
    sinhd_d = es * sinhz_z(d)            ! exp(s) sinh(d)/d
    d_sinhd_d = es * d_sinhz_z3(d)       ! exp(s) [ (sinh(d)/d)' / d ]
    dh2 = 0.5 * dEv(2)                   ! h' / 2
    diag = sinhd_d * ( dh2 - Ev(1)*dEv(1) ) + &
      &    d_sinhd_d * ( Ev(2)*dEv(1) - dh2*Ev(1) )
    acof = sinhd_d * dEv(1) + d_sinhd_d * dh2

    dEx(1,1) = acof * A(1,1) + sinhd_d * dA(1,1) + diag
    dEx(1,2) = acof * A(1,2) + sinhd_d * dA(1,2)
    dEx(2,1) = acof * A(2,1) + sinhd_d * dA(2,1)
    dEx(2,2) = acof * A(2,2) + sinhd_d * dA(2,2) + diag

  end subroutine dExDt

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module dExDt_M

! $Log$
! Revision 2.1  2003/05/05 23:00:25  livesey
! Merged in feb03 newfwm branch
!
! Revision 1.1.2.2  2003/04/30 02:07:41  vsnyder
! Fix a bogus comment
!
! Revision 1.1.2.1  2003/04/30 01:47:02  vsnyder
! Initial commit
!
