! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module CS_GetEv_M

!--------------------------------------------------------------------
!  Compute the eigenvalues and the derivatives of the eigenvalues a
!  complex 2x2 matrix.

  implicit NONE
  private
  public :: CS_GetEv, CS_GetEvSD, dEdt, dEdtSD

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    &  "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

!------------------------------------------------------  CS_GetEv  -----
  pure subroutine CS_GetEv ( A, Ev, dA, dEv )
!  Compute the eigenvalues of a complex 2x2 matrix A.  Optionallly
!  compute the derivatives of the eigenvalues.

    use CS_ZeroFix_M, only: CS_ZeroFix
    use MLSCOmmon, only: RK => Rp

    complex(rk), intent(in) :: A(2,2)
    complex(rk), intent(out) :: Ev(2)
    complex(rk), intent(in), optional :: dA(2,2)
    complex(rk), intent(out), optional :: dEv(2)
    complex(rk)  W, DW, DU, Q, DQ
    real(rk) :: QA

!  Compute the eigenvalues of the input matrix

!  w = Tr(a)
    w = A(2,2) + A(1,1)

!  The quadratic equation is:  z*z - w*z + det(a)
!  This solution doesn't check for pathologies.  See solve_quad.f9h

    q = w * w - 4.0_rk * ( A(2,2) * A(1,1) - A(1,2) * A(2,1) )
    qa = abs(q)
    q = sqrt(q)

    ev(1) = CS_ZeroFix ( 0.5_rk * (w - q) ) ! first  eigenvalue

    ev(2) = CS_ZeroFix ( 0.5_rk * (w + q) ) ! second eigenvalue

    if ( present(dA) .and. present(dEv) ) then
      ! Compute the derivatives of the eigenvalues

!{ Let $w = \text{Tr} {\bf A}$ and $u = \text{Det} {\bf A}$.  The
!  eigenvalues of ${\bf A}$ are solutions of $\lambda^2 - w \lambda + u =
!  0$; $\lambda = \frac12 ( w \pm q)$, so $\lambda^\prime = \frac12
!  (w^\prime \pm q^\prime)$, where $q^2 = w^2-4u$, $q^\prime = (w w^\prime -
!  2 u^\prime)/q$, $w^\prime = {\bf A}^\prime_{11} + {\bf A}^\prime_{22}$,
!  and $u^\prime = {\bf A}^\prime_{11} {\bf A}_{22} + {\bf A}_{11} {\bf
!  A}^\prime_{22} - {\bf A}^\prime_{12} {\bf A}_{21} - {\bf A}_{12} {\bf
!  A}^\prime_{21}$
!
!  In the code below, Q is $q^2$ from above, but DQ is still $q^\prime$.

      dw = dA(1,1) + dA(2,2)
      if ( qa == 0.0_rk ) then              ! Double Eigenvalue
        dEv = CS_ZeroFix( 0.5_rk * dw )
        return
      end if
      du = dA(1,1) * A(2,2) + dA(2,2) * A(1,1) - &
        &  dA(1,2) * A(2,1) + dA(2,1) * A(1,2)
      dq = (w * dw - 2.0_rk * du) * conjg(q)
      if ( qa >= 1.0_rk ) then
        dq = dq / qa      ! No overflow since qa >= 1
        dev(1) = CS_ZeroFix( 0.5_rk * (dw - dq) )
        dev(2) = CS_ZeroFix( 0.5_rk * (dw + dq) )
      else                ! Don't overflow while computing dq/da
        if ( abs(dq) < qa * huge(qa) ) then ! qa < 1 here, so no overflow
          dq = dq / qa
          dev(1) = CS_ZeroFix( 0.5_rk * (dw - dq) )
          dev(2) = CS_ZeroFix( 0.5_rk * (dw + dq) )
        else ! enormous derivatives due to really close eigenvalues
          dev(1) = cmplx(huge(qa),huge(qa),rk)
          dev(2) = -dev(1)
        end if
      end if
    end if

  end subroutine CS_GetEv

!----------------------------------------------------  CS_GetEvSD  -----
  pure subroutine CS_GetEvSD ( A, Ev, dA, dEv )
!  Compute the sum and difference of the eigenvalues of a complex 2x2 matrix
!  A.  Optionally compute their derivatives.

    use CS_ZeroFix_M, only: CS_ZeroFix
    use MLSCOmmon, only: RK => Rp

    complex(rk), intent(in) :: A(2,2)
    complex(rk), intent(out) :: Ev(2)   ! 1/2 ( sum of eigenvalues ),
                                        ! [ 1/2 ( difference of eigenvalues ) ]**2
    complex(rk), intent(in), optional :: dA(2,2) ! Derivative of A
    complex(rk), intent(out), optional :: dEv(2) ! Derivative of Ev
    complex(rk)  S, DU, Q

!  Compute the sum and the square of the difference of the eigenvalues of the
!  input matrix.  1/2 the sum of the eigenvalues is 1/2 the trace of the matrix.

!   s = Tr(a) =
    s = 0.5_rk * ( A(1,1) + A(2,2) )

!  The quadratic equation is:  z*z - 2 s*z + Det(A)
!  The square of 1/2 (difference of the eigenvalues) is s**2 - Det(A) =

    q = s * s - ( A(2,2) * A(1,1) - A(1,2) * A(2,1) )

    ev(1) = CS_ZeroFix ( s ) ! 1/2 (sum of eigenvalues)

    ev(2) = CS_ZeroFix ( q ) ! square of 1/2 (difference of eigenvalues)

    if ( present(dA) .and. present(dEv) ) then
      ! Compute the derivatives of ev(1) and ev(2)

      dEv(1) = dA(1,1) + dA(2,2)                 ! d tr(A)                    

      du = dA(1,1) * A(2,2) + dA(2,2) * A(1,1) & ! d det(A)                   
       & - dA(1,2) * A(2,1) - dA(2,1) * A(1,2)

      dEv(2) = Ev(1) * dEv(1) - du               ! d (tr(A)/4)**2 - d det(A)  

      dEv(1) = 0.5_rk * dEv(1)                   ! 1/2 d tr(A)                
    end if

  end subroutine CS_GetEvSD

!----------------------------------------------------------  dEdt  -----
!  Compute the derivatives of the eigenvalues of the 2X2 matrix A.

  subroutine dEdt ( A, Ev, dA, dEv )

    use CS_ZeroFix_M, only: CS_ZeroFix
    use MLSCommon, only: RK => Rp

    complex(rk), intent(in) :: A(2,2)   ! The matrix
    complex(rk), intent(in) :: Ev(2)    ! Its eigenvalues
    complex(rk), intent(in) :: dA(2,2)  ! Derivative of the matrix
    complex(rk), intent(out) :: dEv(2)  ! Derivative of its eigenvalues
    complex(rk)  W, DW, DU, Q, DQ
    real(rk) :: QA

!{ Let $w = \text{Tr} {\bf A}$ and $u = \text{Det} {\bf A}$.  The
!  eigenvalues of ${\bf A}$ are solutions of $\lambda^2 - w \lambda + u =
!  0$; $\lambda = \frac12 ( w \pm q)$, so $\lambda^\prime = \frac12
!  (w^\prime \pm q^\prime)$, where $q^2 = w^2-4u$, $q^\prime = (w w^\prime -
!  2 u^\prime)/q$, $w^\prime = {\bf A}^\prime_{11} + {\bf A}^\prime_{22}$,
!  and $u^\prime = {\bf A}^\prime_{11} {\bf A}_{22} + {\bf A}_{11} {\bf
!  A}^\prime_{22} - {\bf A}^\prime_{12} {\bf A}_{21} - {\bf A}_{12} {\bf
!  A}^\prime_{21}$
!
!  In the code below, Q is $q^2$ from above, but DQ is still $q^\prime$.

    w = ev(1) + ev(2)
    dw = dA(1,1) + dA(2,2)
    q = ev(2) - ev(1)
    qa = abs(q)
    if ( qa == 0.0_rk ) then              ! Double Eigenvalue
      dEv = CS_ZeroFix( 0.5_rk * dw )
      return
    end if
    du = dA(1,1) * A(2,2) + dA(2,2) * A(1,1) - &
      &  dA(1,2) * A(2,1) + dA(2,1) * A(1,2)
    dq = (w * dw - 2.0_rk * du) * conjg(q)
    if ( qa >= 1.0_rk ) then
      dq = dq / qa      ! No overflow since qa >= 1
      dev(1) = CS_ZeroFix( 0.5_rk * (dw - dq) )
      dev(2) = CS_ZeroFix( 0.5_rk * (dw + dq) )
    else                ! Don't overflow while computing dq/da
      if ( abs(dq) < qa * huge(qa) ) then ! qa < 1 here, so no overflow
        dq = dq / qa
        dev(1) = CS_ZeroFix( 0.5_rk * (dw - dq) )
        dev(2) = CS_ZeroFix( 0.5_rk * (dw + dq) )
      else ! enormous derivatives due to really close eigenvalues
        dev(1) = cmplx(huge(qa),huge(qa),rk)
        dev(2) = -dev(1)
      end if
    end if

  end subroutine dEdt

!--------------------------------------------------------  dEdtSD  -----
!  Compute the derivatives of 1/2 ( the sum of the eigenvalues ) and the
!  square of 1/2 ( the difference of the eigenvalues ) of the 2X2 matrix A.

  subroutine dEdtSD ( A, Ev, dA, dEv )

    use MLSCommon, only: RK => Rp

    complex(rk), intent(in) :: A(2,2)   ! The matrix
    complex(rk), intent(in) :: Ev(2)    ! 1/2 ( sum of eigenvalues ),
                                        ! [ 1/2 ( difference of eigenvalues ) ]**2
    complex(rk), intent(in) :: dA(2,2)  ! Derivative of the matrix
    complex(rk), intent(out) :: dEv(2)  ! Derivative of 1/2 ( sum of eigenvalues ),
                                        ! Derivative of Ev
    complex(rk)  DU

    ! Compute the derivatives of 1/2 ( sum of the eigenvalues ) and
    ! the square of 1/2 ( difference of the eigenvalues ).

    ! Compute the derivatives of ev(1) and ev(2)

    dEv(1) = dA(1,1) + dA(2,2)                 ! d tr(A)                    

    du = dA(1,1) * A(2,2) + dA(2,2) * A(1,1) & ! d det(A)                   
     & - dA(1,2) * A(2,1) - dA(2,1) * A(1,2)

    dEv(2) = Ev(1) * dEv(1) - du               ! d (tr(A)/4)**2 - d det(A)  

    dEv(1) = 0.5_rk * dEv(1)                   ! 1/2 d tr(A)                

  end subroutine dEdtSD

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module CS_GetEv_M

! $Log$
! Revision 2.5  2003/12/19 20:46:18  vsnyder
! Beautify a comment
!
! Revision 2.4  2003/12/18 00:38:47  vsnyder
! Cosmetic changes
!
! Revision 2.3  2003/02/06 18:58:58  vsnyder
! Put an explicit kind on CMPLX references, otherwise NAG gets an overflow
!
! Revision 2.2  2003/02/05 21:47:41  vsnyder
! Remove a USE that wasn't used
!
! Revision 2.1  2003/02/04 01:41:33  vsnyder
! Initial commit
!
