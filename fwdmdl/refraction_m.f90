module REFRACTION_M
  use MLSCommon, only: R8, RP, IP
  use GLNP, only: NG, GX, GW

  Implicit None

  Private
  Public :: refractive_index, comp_refcor, path_ds_dh

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
 "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
 "$RCSfile$"
!---------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------

SUBROUTINE refractive_index(p_path,t_path,n_path,h2o_path)

! This routine computes the refractive index as a function of altitude
! and phi. The returned value has one subtracted from it
! We could easily make this elemental but it might run slower due
! to multiple executions of IF(PRESENT(...))
!  ===============================================================
!  Declaration of variables for sub-program: refractive_index
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
! inputs
  REAL(rp), INTENT(IN) :: p_path(:) ! pressure(hPa) vector
  REAL(rp), INTENT(IN) :: t_path(:) ! temperature vector(K)
! output
  REAL(rp), INTENT(OUT) :: n_path(:) ! refractive indicies - 1
! Keywords
  REAL(rp), OPTIONAL, INTENT(IN) :: h2o_path(:) ! H2O vmr(ppv)
!  ----------------
!  Local variables:
!  ----------------
  REAL(rp), PARAMETER :: const1 = 0.0000776_rp
  REAL(rp), PARAMETER :: const2 = 4810.0_rp
! begin code
  n_path = const1 * p_path / t_path
  IF(PRESENT(h2o_path)) n_path = n_path*(1.0_rp + const2*h2o_path/t_path)
  RETURN

END SUBROUTINE refractive_index

!---------------------------------------------------------------
! This routine computes the integral described in Eqn. 8.11 of the
! MLS ATBD, pg. 46,  using the Gauss-Legendre method.
!
! For derivation of the code below, please see: "FWD Model" paper,
! Page 16, Eqn. 26 & 27
!
Subroutine comp_refcor(h_path,n_path,ht,del_s,ref_corr)
!
    REAL(rp), intent(in) :: H_PATH(:)
    REAL(rp), intent(in) :: N_PATH(:)

    Real(rp), intent(in) :: ht
!
    REAL(rp), INTENT(out) :: Del_s(:)
    Real(rp), intent(out) :: REF_CORR(:)

    Integer(ip) :: i,j,k,no_ele,mid

    Real(rp) :: INTEGRAND_GL(Ng)

    Real(rp) :: q, htan2, Nt2Ht2
    Real(rp) :: H,N,dndh,x1,x2,h1,h2,n1,n2,xm,ym,NH,eps
!
    no_ele = Size(n_path)
    mid = (no_ele + 1) / 2
!
!  Initialize the ref_corr array:
!
    ref_corr(1:no_ele) = 1.0_rp
!
    Del_s = 0.0
!
    htan2 = ht * ht
    Nt2Ht2 = (n_path(mid)*ht)**2
!
    i = 2
    j = mid
    Del_s(i:j) = &
       &    abs(Sqrt(abs(h_path(i-1:j-1)**2-htan2)) -  &
       &        Sqrt(abs(h_path( i : j )**2-htan2)))

    i = j+1
    j = no_ele-1
    Del_s(i:j) = &
       &    abs(Sqrt(abs(h_path(i+1:j+1)**2-htan2)) -  &
       &        Sqrt(abs(h_path( i : j )**2-htan2)))

! First, do the right hand side of the ray path:
!
    h2 = h_path(1)
    n2 = n_path(1)-1.0_rp
!
    q = (h_path(1)*n_path(1))**2-Nt2Ht2
    if(q < 1.0e-9_rp) q = 0.0_rp
    x2 = Sqrt(q)

    do j = 2, mid
!
      x1 = x2
      h1 = h2
      n1 = n2
      h2 = h_path(j)
      n2 = n_path(j)-1.0_rp

      q = (h_path(j)*n_path(j))**2 - Nt2Ht2
      if(q < 1.0e-9_rp) q = 0.0_rp
      x2 = Sqrt(q)

      eps = Log(n2/n1)/(h2-h1)
      xm = 0.5_rp *(x1 + x2)
      ym = 0.5_rp *(x1 - x2)
      do k = 1, Ng
        q = xm + ym * Gx(k)
        NH = Sqrt(q*q + Nt2Ht2)
        Call Solve_Hn(NH)
        Integrand_GL(k) = 1.0_rp/(N+H*dndh)
      end do
!
      q = SUM(integrand_GL(1:)*Gw(1:))
!
! And Finally - define the refraction correction:
!
      ref_corr(j) = q * ym / Del_s(j)
!
    end do
!
! Now, do the left hand side of the ray path:
!
    j = mid+1
    h2 = h_path(j)
    n2 = n_path(j)-1.0_rp

    q = (h_path(j)*n_path(j))**2 - Nt2Ht2
    if(q < 1.0e-9_rp) q = 0.0_rp
    x2 = Sqrt(q)

    do j = mid+1, no_ele-1
!
      x1 = x2
      h1 = h2
      n1 = n2
      h2 = h_path(j+1)
      n2 = n_path(j+1)-1.0_rp

      q = (h_path(j+1)*n_path(j+1))**2 - Nt2Ht2
      if(q < 1.0e-9_rp) q = 0.0_rp
      x2 = Sqrt(q)

      eps = Log(n2/n1)/(h2-h1)
      xm = 0.5_rp *(x2 + x1)
      ym = 0.5_rp *(x2 - x1)
      do k = 1, Ng
        q = xm + ym * Gx(k)
        NH = Sqrt(q*q + Nt2Ht2)
        Call Solve_Hn(NH)
        Integrand_GL(k) = 1.0_rp/(N+H*dndh)
      end do
!
      q = SUM(integrand_GL(1:)*Gw(1:))
!
! And Finally - define the refraction correction:
!
      ref_corr(j) = q * ym / Del_s(j)
!
    end do
!
    Return
!
Contains
!------------------------------------------------------------------
! Solve the equation h*(1.0+N(h)) = N*H, where N(h) is an exponential:
!    N(h) = n1*Exp(eps*(h-h1))

  Subroutine Solve_Hn(NH)

    Real(rp), intent(in) :: NH

    Integer :: iter
    Logical :: bound
    Real(rp) :: v1,v2,f1,f2,df,fpos,fneg,hpos,hneg

    Integer,  PARAMETER :: Max_Iter = 20
    Real(rp), PARAMETER :: Tiny = 1.0e-8_rp

     bound = .FALSE.
     f1 = h1 * (1.0 + n1) - NH
     f2 = h2 * (1.0 + n2) - NH
     df = (f2 - f1) / (h2 - h1)

     if(f1*f2 < 0.0_rp) then
       bound = .TRUE.
       if(f1 < 0.0_rp) then
         fneg = f1
         hneg = h1
         fpos = f2
         hpos = h2
       else
         fpos = f1
         hpos = h1
         fneg = f2
         hneg = h2
       endif
     else
       Print *,'** Warning from Solve_Hn: ROOT is NOT bound ..'
     endif

     iter = 1
     v2 = (h1*abs(f2)+h2*abs(f1))/(abs(f1)+abs(f2))
     f2 = v2*(1.0+n1*Exp(eps*(v2-h1)))-NH

     DO

       v1 = v2
       f1 = f2

       v2 = v1 - f1 / df

       if(bound) then
         if(v2 < min(hpos,hneg) .OR. v2 > max(hpos,hneg)) v2 = 0.5_rp*(hneg + hpos)
       endif

       f2 = v2*(1.0+n1*Exp(eps*(v2-h1)))-NH

       if(abs(f2) < Tiny .OR. abs(v2-v1) < Tiny) EXIT

       if(Iter == Max_Iter) EXIT

       iter = iter + 1
       df = (f2 - f1) / (v2 - v1)

       if(bound) then
         if(f2 < 0.0_rp) then
           fneg = f2
           hneg = v2
         else
           fpos = f2
           hpos = v2
         endif
         df = (fpos - fneg) / (hpos - hneg)
       endif

     END DO

     if(abs(f2) >= Tiny .AND. abs(v2-v1) >= Tiny) then
       Print *,'** Warning from Solve_Hn: DID NOT converged within ',Max_Iter, &
             & ' iterations ..'
     endif

     H = v2
     df = n1 * Exp(eps*(H-h1))
     dndh = eps * df
     N = 1.0_rp + df

  End Subroutine Solve_Hn

End Subroutine comp_refcor

!------------------------------------------------------------------
  ELEMENTAL REAL(rp) FUNCTION path_ds_dh(r_path,r_tan)
!
! inputs:
!
  REAL(rp), INTENT(in) :: r_path ! heights + req (km).
  REAL(rp), INTENT(in) :: r_tan ! tangent height + req (km).
!
! output:
!  REAL(rp), INTENT(out) :: path_ds_dh ! path length derivative wrt height.
!
! calculation
!
  path_ds_dh = r_path / SQRT(r_path**2 - r_tan**2)
!
  END FUNCTION path_ds_dh
!
!------------------------------------------------------------------

END module REFRACTION_M
! $Log$
! Revision 2.3  2002/02/16 10:32:18  zvi
! Make sure iteration in Solve_HN do not diverge
!
! Revision 2.2  2002/02/14 21:36:13  zvi
! Fix Sqrt() problem..
!
! Revision 2.1  2001/12/01 01:35:22  zvi
! Clerifying code.. easier to follow..
!
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model
!
! Revision 1.9.2.2  2001/09/12 21:38:53  zvi
! Added CVS stuff
!
! Revision 1.9.2.1  2001/09/10 10:02:32  zvi
! Cleanup..comp_path_entities_m.f90
!
! Revision 1.1 2000/06/09 00:08:14  Z.Shippony
! Initial conversion to Fortran 90
