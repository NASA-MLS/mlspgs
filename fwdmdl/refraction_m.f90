! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module REFRACTION_M

  use MLSCommon, only: RP
    
  implicit none

  private
  public :: Refractive_index, Comp_refcor, Path_ds_dh

  real(rp), parameter, public :: RefrAterm = 0.0000776_rp
  real(rp), parameter, public :: RefrBterm = 4810.0_rp

  interface Refractive_index
    module procedure Refractive_index_1, Refractive_index_0
  end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains

!--------------------------------------------  Refractive_index_1  -----
  subroutine Refractive_index_1 ( p_path, t_path, n_path, h2o_path )

  ! This routine computes the refractive index as a function of altitude
  ! and phi. The returned value has one subtracted from it
  ! We could easily make this elemental but it might run slower due
  ! to multiple executions of if ( PRESENT(...))
  !  ===============================================================
  !  Declaration of variables for sub-program: refractive_index
  !  ===============================================================
  !  ---------------------------
  !  Calling sequence variables:
  !  ---------------------------
  ! inputs
    real(rp), intent(in) :: p_path(:) ! pressure(hPa) vector
    real(rp), intent(in) :: t_path(:) ! temperature vector(K)
  ! output
    real(rp), intent(out) :: n_path(:) ! refractive indicies - 1
  ! Keywords
    real(rp), optional, intent(in) :: h2o_path(:) ! H2O vmr(ppv)

  ! begin code
    n_path = refrAterm * p_path / t_path
    if ( present(h2o_path) ) &
      & n_path = n_path *( 1.0_rp + refrBterm*h2o_path/t_path)
    return

  end subroutine Refractive_index_1

!--------------------------------------------  Refractive_index_0  -----
  subroutine Refractive_index_0 ( p_path, t_path, n_path, h2o_path )

  ! This routine computes the refractive index as a function of altitude
  ! and phi. The returned value has one subtracted from it
  ! We could easily make this elemental but it might run slower due
  ! to multiple executions of if ( PRESENT(...))
  !  ===============================================================
  !  Declaration of variables for sub-program: refractive_index
  !  ===============================================================
  !  ---------------------------
  !  Calling sequence variables:
  !  ---------------------------
  ! inputs
    real(rp), intent(in) :: p_path ! pressure(hPa) vector
    real(rp), intent(in) :: t_path ! temperature vector(K)
  ! output
    real(rp), intent(out) :: n_path ! refractive indicies - 1
  ! Keywords
    real(rp), optional, intent(in) :: h2o_path ! H2O vmr(ppv)

  ! begin code
    n_path = refrAterm * p_path / t_path
    if ( present(h2o_path) ) &
      & n_path = n_path *( 1.0_rp + refrBterm*h2o_path/t_path)
    return

  end subroutine Refractive_index_0

! --------------------------------------------------  Comp_refcor  -----

  subroutine Comp_refcor ( h_path, n_path, ht, del_s, ref_corr )

  ! This routine computes the integral described in Eqn. 8.11 of the
  ! MLS ATBD, pg. 44,  using the Gauss-Legendre method.

  ! For derivation of the code below, please see: "FWD Model" paper,
  ! Page 16, Eqn. 26 & 27

    use MLSCommon, only: RP, IP
    use GLNP, only: NG, GX, GW
    use MLSMessageModule, only: MLSMessage, MLSMSG_Warning

    real(rp), intent(in) :: H_PATH(:)
    real(rp), intent(in) :: N_PATH(:)

    real(rp), intent(in) :: ht

    real(rp), intent(out) :: Del_s(:)
    real(rp), intent(out) :: REF_CORR(:)

    integer(ip) :: i, j, k, no_ele, mid

    real(rp) :: INTEGRAND_GL(Ng)

    real(rp) :: q, htan2, Nt2Ht2
    real(rp) :: H, N, dndh, x1, x2, h1, h2, n1, n2, xm, ym, NH, eps

    real(rp), parameter :: Tiny = 1.0e-8_rp

    no_ele = size(n_path)
    mid = (no_ele + 1) / 2

  !  Initialize the ref_corr array:

    ref_corr(1:no_ele) = 1.0_rp

    Del_s = 0.0

    htan2 = ht * ht
    Nt2Ht2 = (n_path(mid)*ht)**2

    i = 2
    j = mid
    Del_s(i:j) = &
       &    abs(sqrt(abs(h_path(i-1:j-1)**2-htan2)) -  &
       &        sqrt(abs(h_path( i : j )**2-htan2)))

    i = j+1
    j = no_ele-1
    Del_s(i:j) = &
       &    abs(sqrt(abs(h_path(i+1:j+1)**2-htan2)) -  &
       &        sqrt(abs(h_path( i : j )**2-htan2)))

  ! First, do the right hand side of the ray path:

    h2 = h_path(1)
    n2 = n_path(1)-1.0_rp

    q = (h_path(1)*n_path(1))**2 - nt2ht2
    if ( abs(q) < tiny ) q = 0.0_rp
    x2 = sqrt(abs(q))

o1: do j = 2, mid

      x1 = x2
      h1 = h2
      n1 = n2
      h2 = h_path(j)
      n2 = n_path(j)-1.0_rp

      q = (h_path(j)*n_path(j))**2 - nt2ht2
      if ( abs(q) < tiny) q = 0.0_rp

      if ( q < 0.0_rp .or. n1*n2 <= 0.0_rp ) then
        ref_corr(j) = ref_corr(j-1)
        cycle
      end if

      x2 = sqrt(q)

      eps = log(n2/n1)/(h2-h1)
      xm = 0.5_rp *(x1 + x2)
      ym = 0.5_rp *(x1 - x2)
      do k = 1, ng
        q = xm + ym * gx(k)
        nh = sqrt(q*q + nt2ht2)
        call solve_hn ( nh )
        if ( h < 0.0 ) then
          ref_corr(j) = ref_corr(j-1)
          cycle o1
        end if
        integrand_gl(k) = 1.0_rp/(n+h*dndh)
      end do

  ! And Finally - define the refraction correction:

      ref_corr(j) = dot_product(integrand_gl,gw) * ym / Del_s(j)

    end do o1

  ! Now, do the left hand side of the ray path:

    j = mid+1
    h2 = h_path(j)
    n2 = n_path(j)-1.0_rp

    q = (h_path(j)*n_path(j))**2 - nt2ht2
    if ( abs(q) < tiny) q = 0.0_rp
    x2 = sqrt(abs(q))

o2: do j = mid+1, no_ele-1

      x1 = x2
      h1 = h2
      n1 = n2
      h2 = h_path(j+1)
      n2 = n_path(j+1)-1.0_rp

      q = (h_path(j+1)*n_path(j+1))**2 - nt2ht2
      if ( abs(q) < tiny) q = 0.0_rp

      if ( q < 0.0_rp .or. n1*n2 <= 0.0_rp ) then
        ref_corr(j) = ref_corr(j-1)
        cycle
      end if

      x2 = sqrt(q)

      eps = log(n2/n1)/(h2-h1)
      xm = 0.5_rp *(x2 + x1)
      ym = 0.5_rp *(x2 - x1)
      do k = 1, ng
        q = xm + ym * gx(k)
        nh = sqrt(q*q + nt2ht2)
        call solve_hn ( nh )
        if ( h < 0.0 ) then
          ref_corr(j) = ref_corr(j-1)
          cycle o2
        end if
        integrand_gl(k) = 1.0_rp/(n+h*dndh)
      end do

  ! And Finally - define the refraction correction:

      ref_corr(j) = dot_product(integrand_gl,gw) * ym / Del_s(j)

    end do o2

    return

  contains
  !------------------------------------------------------------------
  ! Solve the equation h*(1.0+N(h)) = N*H, where N(h) is an exponential:
  !    N(h) = n1*Exp(eps*(h-h1))

    subroutine Solve_Hn ( NH )

      real(rp), intent(in) :: NH

      integer :: iter
      real(rp) :: E ! exp(eps * (v2-h1))
      real(rp) :: v1, v2, f1, f2, df, hpos, hneg

      integer,  parameter :: Max_Iter = 20

      character(LEN=*), parameter :: Msg1 = &
        & 'From Solve_Hn routine: Could not bracket the root'
      character(LEN=*), parameter :: Msg2 = &
        & 'From Solve_Hn routine: Did not converge within 20 iterations'

       f1 = h1 * (1.0_rp + n1) - NH
       f2 = h2 * (1.0_rp + n2) - NH

       if ( f1*f2 > 0.0_rp ) then
         H = -1.0_rp
         call MLSMessage ( MLSMSG_Warning, ModuleName, Msg1)
         return 
       end if

       if ( f1 <= 0.0_rp ) then
         hneg = h1
         hpos = h2
       else
         hpos = h1
         hneg = h2
       end if

       iter = 1
       v2 = (h1 * abs(f2) + h2 * abs(f1)) / (abs(f1) + abs(f2))
       e = n1 * exp(eps*(v2-h1))
       f2 = v2 * (1.0_rp + e ) - NH
       df = 1.0_rp + e * ( 1.0_rp + eps * v2 )

       do

         v1 = v2
         f1 = f2

         v2 = v1 - f1 / df

         if ( v2 < min(hpos,hneg) .OR. v2 > max(hpos,hneg) ) &
             &  v2 = 0.5_rp * (hneg + hpos)

         e = n1 * exp(eps*(v2-h1))
         f2 = v2 * (1.0_rp + e ) - NH

         if ( abs(f2) < tiny .or. abs(v2-v1) < tiny ) exit

         if ( iter >= max_iter ) exit

         if ( f2 < 0.0_rp ) then
           hneg = v2
         else
           hpos = v2
         end if

         iter = iter + 1
         df = 1.0_rp + e * ( 1.0_rp + eps * v2 )

       end do

       if ( abs(f2) >= tiny .and. abs(v2-v1) >= tiny ) &
           & call MLSMessage ( MLSMSG_Warning, ModuleName, Msg2 )

       H = v2
       df = e
       dndh = eps * df
       N = 1.0_rp + df

    end subroutine Solve_Hn

  end subroutine Comp_refcor

! ---------------------------------------------------  Path_ds_dh  -----
  elemental real(rp) function Path_ds_dh ( r_path, r_tan )

  ! inputs:

    real(rp), intent(in) :: r_path ! heights + req (km).
    real(rp), intent(in) :: r_tan ! tangent height + req (km).

  ! output:
  !  REAL(rp), INTENT(out) :: path_ds_dh ! path length derivative wrt height.

  ! calculation

    path_ds_dh = r_path / sqrt(r_path**2 - r_tan**2)

  end function Path_ds_dh

!------------------------------------------------------------------

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

END module REFRACTION_M
! $Log$
! Revision 2.13  2002/10/08 17:08:06  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.12  2002/09/26 21:03:02  vsnyder
! Publish two constants, make Refractive_Index generic
!
! Revision 2.11  2002/09/26 18:02:36  livesey
! Bug fix (wouldn't compile)
!
! Revision 2.10  2002/09/26 00:27:55  vsnyder
! Insert copyright notice, move USEs from module scope to procedure scope,
! cosmetic changes.
!
! Revision 2.9  2002/03/15 06:53:02  zvi
! Some cosmetic changes
!
! Revision 2.8  2002/03/14 22:33:30  zvi
! Add protection against Log() blowout
!
! Revision 2.7  2002/03/14 20:31:14  zvi
! Make comp_refcor more robust
!
! Revision 2.6  2002/02/18 06:58:04  zvi
! Trimming some unused code..
!
! Revision 2.5  2002/02/18 01:01:58  zvi
! Let the program crash & burn for LARGE negative Sqrt Arg.
!
! Revision 2.4  2002/02/17 03:23:40  zvi
! Better code for convergance in Solve_Hn
!
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
