! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module CS_ExpMat_M

  implicit NONE
  private
  public :: CS_ExpMat

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

! ----------------------------------------------------  CS_ExpMat  -----
  subroutine CS_ExpMat ( A, Ex, Status )

!  Compute the exponential function of a complex, 2x2 matrix, A

    use CRREXP_M, only: CRREXP
    use CS_GetEv_M, only: CS_GetEv
    use CS_ZeroFix_M, only: CS_ZeroFix
    use MLSCommon, only: RK => Rp

    complex(rk), intent(in) :: A(2,2)
    complex(rk), intent(out) :: Ex(2,2)
    integer, optional, intent(out) :: Status ! zero => OK, else overflow in EXP

    real(rk), save :: Big = -1.0_rk          ! Log(Huge(0.0_rk)), eventually
    complex(rk) :: EH    ! (exp(h)-1) / h
    complex(rk) :: Ev(2) ! Eigenvalues of A
    complex(rk) :: EZ2   ! exp(z2)
    complex(rk) :: H     ! Ev(1) - Ev(2) = z1 - z2
    complex(rk) :: W

!{ Use Sylvester's identity for a function of a 2x2 matrix {\bf A}, whose
!  eigenvalues are: $z_1$ and $z_2$  ({\bf I} is the 2x2 Identity matrix)
!
!    $F({\bf A}) = \frac{F(z_1)}
!                       {z_1-z_2} ({\bf A} - z_2 {\bf I}) +
!                  \frac{F(z_2)}
!                       {z_2-z_1} ({\bf A} - z_1 {\bf I})$
!
!                              Or:
!
!    $F({\bf A}) = \frac{F(z_1) ({\bf A} - z_2 {\bf I}) - 
!                        F(z_2) ({\bf A} - z_1 {\bf I})}
!                       {z_1 - z_2}$
!
!  For any function $F()$ (in this case: $\exp()$ )
!
! Now, rearrange this to get
!%
! $e^{z_2} \left [ \frac{e^h-1}h ( {\bf A} - z_2 {\bf I} ) + {\bf I} \right ]$.
!%
! where $h = z_1 - z_2$.  This is well behaved as $h \rightarrow 0$, so we
! don't need to futz with L'H\^opital's rule.  Since
! $\lim_{h \rightarrow 0}\frac{\text{d}}{\text{d}h} \frac{e^h-1}h = 1$,
! the relative error in the first term no larger than the relative error
! in $h$ when $h$ is small.

!{This could also be written as
! $e^s \left[ \frac{\sinh d}d ( \mathbf{A} - s\, \mathbf{I}) +
!             \cosh d\, \mathbf{I} \right ]$, where $s = \frac12 (z_1+z_2)$
! and $d = \frac12 (z_1-z_2)$.  This is useful to compute the derivative
! without getting into trouble as $d \rightarrow 0$, but the present form
! is OK here, and more efficient.

    if ( abs(A(1,1)) + abs(A(1,2)) + abs(A(2,1)) + abs(A(2,2)) < &
      & epsilon(1.0_rk) ) then
    ! Use the first two terms of Taylor's series for exp(A) = I + A
      Ex(1,1) = A(1,1) + 1.0_rk
      Ex(1,2) = A(1,2)
      Ex(2,1) = A(2,1)
      Ex(2,2) = A(2,2) + 1.0_rk
      if ( present(status) ) status = 0
      return
    end if

!  First, Get the eigenvalues of the input matrix

    call CS_GetEv ( A, Ev )

!  Now compute Sylvester's identity (as rearranged to avoid trouble
!  as h -> 0).

    if ( present(status) ) then
      if ( big < 0.0_rk ) big = log(huge(0.0_rk))
      if ( real(ev(2)) >= big ) then
        status = 1
        return
      end if
      status = 0
    end if

    ez2 = exp(Ev(2))
    h = ev(1) - ev(2) ! z1 - z2

    if ( present(status) ) then
      if ( real(h) >= big ) then
        status = 1
        return
      end if
    end if

    eh = crrexp ( h ) ! (exp(h)-1)/h
    w = ez2 * eh      ! exp(z2) * (exp(h)-1)/h

    Ex(1,1) = cs_zerofix( w * (A(1,1)-ev(2)) + ez2 )

    Ex(1,2) = cs_zerofix( w * A(1,2) )

    Ex(2,1) = cs_zerofix( w * A(2,1) )

    Ex(2,2) = cs_zerofix( w * (A(2,2)-ev(2)) + ez2 )

  end subroutine CS_ExpMat

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module CS_ExpMat_M

! $Log$
! Revision 2.8  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.7  2005/06/22 18:08:18  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.6  2004/04/21 18:23:56  vsnyder
! Correct a TeXnicality
!
! Revision 2.5  2003/07/02 17:19:38  vsnyder
! Don't set status if it's not present!
!
! Revision 2.4  2003/07/02 00:43:16  vsnyder
! Make sure STATUS gets set
!
! Revision 2.3  2003/06/27 22:05:18  vsnyder
! Add 'status' argument
!
! Revision 2.2  2003/05/05 23:00:25  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.1.2.1  2003/03/06 21:53:23  vsnyder
! Correct a sign error
!
! Revision 2.1  2003/02/04 01:41:33  vsnyder
! Initial commit
!
