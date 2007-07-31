! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module REFRACTION_M

  use MLSCommon, only: RP
    
  implicit none

  private
  public :: Refractive_index, Refractive_Index_H2O_Update, Comp_refcor

  real(rp), parameter, public :: RefrAterm = 0.0000776_rp
  real(rp), parameter, public :: RefrBterm = 4810.0_rp

  interface Refractive_index
    module procedure Refractive_index_0,   Refractive_index_0_h2o
    module procedure Refractive_index_1,   Refractive_index_1_h2o
    module procedure Refractive_index_1_2, Refractive_index_1_2_h2o
    module procedure Refractive_index_2,   Refractive_index_2_h2o
  end interface

  interface Refractive_Index_H2O_Update
    module procedure Refractive_Index_H2O_Update_0
    module procedure Refractive_Index_H2O_Update_1
    module procedure Refractive_Index_H2O_Update_2
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains

!--------------------------------------------  Refractive_index_0  -----
  subroutine Refractive_index_0 ( p_path, t_path, n_path )

  ! This routine computes the refractive index as a function of altitude
  ! and phi. The returned value has one subtracted from it
  ! We could easily make this elemental.

  ! inputs
    real(rp), intent(in) :: p_path ! pressure(hPa)
    real(rp), intent(in) :: t_path ! temperature(K)
  ! output
    real(rp), intent(out) :: n_path ! refractive indicies - 1

  ! begin code
    n_path = refrAterm * p_path / t_path

  end subroutine Refractive_index_0

!----------------------------------------  Refractive_index_0_h2o  -----
  subroutine Refractive_index_0_h2o ( p_path, t_path, n_path, h2o_path )

  ! This routine computes the refractive index as a function of altitude
  ! and phi. The returned value has one subtracted from it

  ! inputs
    real(rp), intent(in) :: p_path ! pressure(hPa)
    real(rp), intent(in) :: t_path ! temperature(K)
  ! output
    real(rp), intent(out) :: n_path ! refractive indicies - 1
  ! Keywords
    real(rp), intent(in) :: h2o_path ! H2O vmr(ppv)

  ! begin code
    n_path = refrAterm * p_path / t_path * ( 1.0_rp + refrBterm*h2o_path/t_path)

  end subroutine Refractive_index_0_h2o

! --------------------------------  Refractive_Index_H2O_Update_0  -----
  subroutine Refractive_Index_H2O_Update_0 ( t_path, h2o_path, n_path )

  ! Update the refractive index N_Path previously computed using only
  ! temperature and pressure to account for H2O.

  ! inputs
    real(rp), intent(in) :: t_path    ! temperature(K)
    real(rp), intent(in) :: h2o_path  ! H2O vmr(ppv)
  ! inout
    real(rp), intent(inout) :: n_path ! refractive indicies - 1

    n_path = n_path * ( 1.0_rp + refrBterm*h2o_path/t_path)

  end subroutine Refractive_Index_H2O_Update_0

! -------------------------------------------  Refractive_index_1  -----
  subroutine Refractive_index_1 ( p_path, t_path, n_path )

  ! This routine computes the refractive index as a function of altitude
  ! and phi. The returned value has one subtracted from it.

  ! inputs
    real(rp), intent(in) :: p_path(:) ! pressure(hPa)
    real(rp), intent(in) :: t_path(:) ! temperature(K)
  ! output
    real(rp), intent(out) :: n_path(:) ! refractive indicies - 1

  ! begin code
    n_path = refrAterm * p_path / t_path

  end subroutine Refractive_index_1

! ---------------------------------------  Refractive_index_1_h2o  -----
  subroutine Refractive_index_1_h2o ( p_path, t_path, n_path, h2o_path )

  ! This routine computes the refractive index as a function of altitude
  ! and phi. The returned value has one subtracted from it.

  ! inputs
    real(rp), intent(in) :: p_path(:) ! pressure(hPa)
    real(rp), intent(in) :: t_path(:) ! temperature(K)
  ! output
    real(rp), intent(out) :: n_path(:) ! refractive indicies - 1
  ! Keywords
    real(rp), intent(in) :: h2o_path(:) ! H2O vmr(ppv)

  ! begin code
    n_path = refrAterm * p_path / t_path * ( 1.0_rp + refrBterm*h2o_path/t_path)

  end subroutine Refractive_index_1_h2o

! --------------------------------  Refractive_Index_H2O_Update_1  -----
  subroutine Refractive_Index_H2O_Update_1 ( t_path, h2o_path, n_path )

  ! Update the refractive index N_Path previously computed using only
  ! temperature and pressure to account for H2O.

  ! inputs
    real(rp), intent(in) :: t_path(:)    ! temperature(K)
    real(rp), intent(in) :: h2o_path(:)  ! H2O vmr(ppv)
  ! inout
    real(rp), intent(inout) :: n_path(:) ! refractive indicies - 1

    n_path = n_path * ( 1.0_rp + refrBterm*h2o_path/t_path)

  end subroutine Refractive_Index_H2O_Update_1

! -----------------------------------------  Refractive_index_1_2  -----
  subroutine Refractive_index_1_2 ( p_path, t_path, n_path )

  ! This routine computes the refractive index as a function of altitude
  ! and phi. The returned value has one subtracted from it.

  ! inputs
    real(rp), intent(in) :: p_path(:)   ! pressure(hPa)
    real(rp), intent(in) :: t_path(:,:) ! temperature(K)
  ! output
    real(rp), intent(out) :: n_path(:,:) ! refractive indicies - 1

    integer :: I

  ! begin code
    do i = 1, size(t_path,2)
      n_path(:,i) = refrAterm * p_path(:) / t_path(:,i)
    end do

  end subroutine Refractive_index_1_2

! -------------------------------------  Refractive_index_1_2_h2o  -----
  subroutine Refractive_index_1_2_h2o ( p_path, t_path, n_path, h2o_path )

  ! This routine computes the refractive index as a function of altitude
  ! and phi. The returned value has one subtracted from it.

  ! inputs
    real(rp), intent(in) :: p_path(:)   ! pressure(hPa)
    real(rp), intent(in) :: t_path(:,:) ! temperature(K)
  ! output
    real(rp), intent(out) :: n_path(:,:) ! refractive indicies - 1
  ! Keywords
    real(rp), intent(in) :: h2o_path(:,:) ! H2O vmr(ppv)

    integer :: I

  ! begin code
    do i = 1, size(t_path,2)
      n_path(:,i) = refrAterm * p_path(:) / t_path(:,i) * &
        & ( 1.0_rp + refrBterm*h2o_path(:,i)/t_path(:,i))
    end do

  end subroutine Refractive_index_1_2_h2o

! -------------------------------------------  Refractive_index_2  -----
  subroutine Refractive_index_2 ( p_path, t_path, n_path )

  ! This routine computes the refractive index as a function of altitude
  ! and phi. The returned value has one subtracted from it.

  ! inputs
    real(rp), intent(in) :: p_path(:,:) ! pressure(hPa)
    real(rp), intent(in) :: t_path(:,:) ! temperature(K)
  ! output
    real(rp), intent(out) :: n_path(:,:) ! refractive indicies - 1

  ! begin code
    n_path = refrAterm * p_path / t_path

  end subroutine Refractive_index_2

! ---------------------------------------  Refractive_index_2_h2o  -----
  subroutine Refractive_index_2_h2o ( p_path, t_path, n_path, h2o_path )

  ! This routine computes the refractive index as a function of altitude
  ! and phi. The returned value has one subtracted from it.

  ! inputs
    real(rp), intent(in) :: p_path(:,:) ! pressure(hPa)
    real(rp), intent(in) :: t_path(:,:) ! temperature(K)
  ! output
    real(rp), intent(out) :: n_path(:,:) ! refractive indicies - 1
  ! Keywords
    real(rp), intent(in) :: h2o_path(:,:) ! H2O vmr(ppv)

  ! begin code
    n_path = refrAterm * p_path / t_path * ( 1.0_rp + refrBterm*h2o_path/t_path)

  end subroutine Refractive_index_2_h2o

! --------------------------------  Refractive_Index_H2O_Update_2  -----
  subroutine Refractive_Index_H2O_Update_2 ( t_path, h2o_path, n_path )

  ! Update the refractive index N_Path previously computed using only
  ! temperature and pressure to account for H2O.

  ! inputs
    real(rp), intent(in) :: t_path(:,:)    ! temperature(K)
    real(rp), intent(in) :: h2o_path(:,:)  ! H2O vmr(ppv)
  ! inout
    real(rp), intent(inout) :: n_path(:,:) ! refractive indicies - 1

    n_path = n_path * ( 1.0_rp + refrBterm*h2o_path/t_path)

  end subroutine Refractive_Index_H2O_Update_2

! --------------------------------------------------  Comp_refcor  -----

  subroutine Comp_refcor ( tan_pt, h_path, n_path, ht, del_s, ref_corr, status )

  !{ This routine computes the integral
  !
  ! \begin{equation*}
  ! I = \Delta s^{\text{refr}}_{i\rightarrow i-1} =
  ! \int_{H_{i-1}}^{H_i} \frac{\mathcal{N}(H) H}
  !                           {\sqrt{\mathcal{N}(H)^2 H^2 - \mathcal{N}_t^2 H_t^2}}
  !                      \,\text{d} H\,,
  ! \end{equation*}
  !
  ! which is described in Eqn. 10.12 of the 19 August 2004 MLS ATBD, pg. 45,
  ! using the Gauss-Legendre method.
  !
  ! Substitute $x = \sqrt{\mathcal{N}(H)^2 H^2 - \mathcal{N}_t^2 H_t^2}$
  ! to remove the singularity, giving
  ! $I = \int_{x_{i-1}}^{x_i}
  !       \frac{\text{d}x}
  !            {N(h)+h\frac{\text{d}}{\text{d}h}N(h)}
  !    = \int_{x_{i-1}}^{x_i}
  !       \frac{\text{d}x}{\frac{\text{d}}{\text{d}h} h N(h)}$,
  ! where $h$ is the solution of $h N(h) = \sqrt{H_t^2 \mathcal{N}_t^2 + x^2}$.
  !
  ! Assume $\mathcal{N}(h) = 1+E$, where $E = n_{i-1} \exp(\epsilon (h-H_{i-1}))$,
  ! $\epsilon = \frac1{H_i-H_{i-1}} \log \frac{n_i}{n{_i-1}}$, and
  ! $n_i = \mathcal{N}(H_i)-1$.
  !
  ! If $\mathcal{N}(H_{i-1}) = \mathcal{N}(H_i)$ the integral is evaluable as
  ! $\left.\frac1{\mathcal{N}(H_i)} \sqrt{\mathcal{N}(H_i) H - \mathcal{N}_t H_t}
  ! \right|_{H=H_{i-1}}^{H_i}$.
  !
  ! If the refractive correction is not in the interval [1.0,1.3], or if
  ! $\frac{\text{d}}{\text{d}h} h N(h)$ changes sign, use the trapezoidal
  ! rule instead.  On the ends of the interval, $h=H_i$ and $E=n_i$, so
  ! $\mathcal{N}(H_i) = \mathcal{N}_i$, giving $I=\frac{x_i-x_{i-1}}2
  ! (f_{i-1}+f_i)$ where the integrand $f_i = \frac1{1+n_i(1+\epsilon H_i)}$.

  ! For derivation of the code below, please see: "FWD Model" paper,
  ! Page 16, Eqn. 26 & 27

    use Dump_0, only: Dump
    use GLNP, only: NG, GX=>gx_all, GW=>gw_all
    use MLSKinds, only: RP, IP
    use MLSMessageModule, only: MLSMessage, MLSMSG_Warning
    use Output_m, only: Output
    use Toggles, only: Switches

    integer, intent(in) :: Tan_pt      ! Tangent point index in H_Path etc.
    real(rp), intent(in) :: H_PATH(:)
    real(rp), intent(in) :: N_PATH(:)

    real(rp), intent(in) :: ht

    real(rp), intent(out) :: Del_s(:)
    real(rp), intent(out) :: REF_CORR(:)

    integer, intent(out) :: Status ! 0 = OK, 1 = failed to bracket root,
                                   ! 2 = too many iterations

    logical :: Bad ! Ref_Corr < 1 or Ref_corr > 1.3 somewhere
    integer(ip) :: HNDP ! "Solve_Hn Detail Printing"
    integer(ip) :: j, j1, j2, k, m, no_ele, stat

    real(rp) :: INTEGRAND_GL(Ng)

    real(rp) :: q, htan2, NH, Nt2Ht2
    real(rp) :: dndh, eps, H, h1, h2, N, n1, n2, t1, t2, x1, x2, xm, ym

    real(rp), parameter :: Htol = 1.0e-3_rp
    real(rp), parameter :: Tiny = 1.0e-8_rp

    status = 0
    hndp = index(switches,'hndp')

    no_ele = size(n_path)

  !  Initialize the ref_corr array:

    ref_corr(1:no_ele) = 1.0_rp

    Del_s = 0.0

    htan2 = ht * ht
    Nt2Ht2 = (n_path(tan_pt)*ht)**2

    bad = .false.

    j1 = tan_pt
    j2 = 2
    do m = -1, 1, 2
      h2 = h_path(j1)
      n2 = n_path(j1)

      q = (h2 * n2)**2 - nt2ht2
      if ( abs(q) < tiny ) q = 0.0_rp
      n2 = n2 - 1.0_rp
      x2 = sqrt(abs(q))

jl:   do j = j1, j2, m

        x1 = x2
        n1 = n2
        h1 = h2
        n2 = n_path(j+m)
        h2 = h_path(j+m)

        del_s(j) = abs(sqrt(abs(h1**2-htan2)) - sqrt(abs(h2**2-htan2)))

        q = (h2 * n2)**2 - nt2ht2
        if ( abs(q) < tiny) q = 0.0_rp
        n2 = n2 - 1.0_rp

        if ( q < 0.0_rp ) then
          if ( hndp > 0 ) then
            write ( *, '("q < 0 for Ref_Corr(",i0,"), q = ",1pg15.7)' ) j, q
            write ( *, * ) "n2, h2, nt2ht2 = ", n2+1.0, h2, nt2ht2
          end if
          ref_corr(j) = 100.0
          bad = .true.
          cycle
        end if

        x2 = sqrt(q)

        if ( abs(n2-n1) <= 1.0e-3*(n2+n1) ) then
          ! Where N is essentially constant, the integral is easy
          if ( hndp > 0 ) write ( *, '("N essentially constant for Ref_Corr(",i0,")")' ) j
          ref_corr(j) = abs( sqrt((n_path(j)*h2)**2-nt2ht2) - &
            &                sqrt((n_path(j)*h1)**2-nt2ht2) ) / &
            &           ( n_path(j) * del_s(j) )
        else if ( n1 <= 0.0 .or. n2 <= 0.0 ) then
          ! Force trapezoidal rule on untransformed integral because
          ! we can't compute eps to do the transformation
          if ( hndp > 0 ) write ( *, '("Trapezoidal for Ref_Corr(",i0,")")' ) j
          t1 = (1.0 + n1) * h1
          t2 = (1.0 + n2) * h2
          ref_corr(j) = 0.5 * abs( ( h2 - h1 ) * &
            & (t1/sqrt(t1**2-nt2ht2) + t2/sqrt(t2**2-nt2ht2)) ) / Del_s(j)
        else
          if ( hndp > 0 ) write ( *, '("Solving for H for Ref_Corr(",i0,")")' ) j
          eps = log(n2/n1)/(h2-h1)
          xm = 0.5_rp *(x2 + x1)      ! Midpoint of the interval
          ym = 0.5_rp *(x2 - x1)      ! Half of the interval length
          do k = 1, ng
            q = xm + ym * gx(k)       ! Gauss abscissa
            NH = sqrt(q*q + nt2ht2)
            ! Solve h*(1+n1*exp(eps*(h-h1))) = sqrt(q*q + nt2ht2) for h
            call solve_hn
            status = max(status,stat)
            if ( stat == 1 ) then
!               ref_corr(j) = ref_corr(j-1+2*m) ! Why did I do this?
!               ref_corr(j) = ref_corr(j-1)
!               cycle jl
              ! Force a different approximation
              integrand_gl = 100.0 * Del_s(j) / ym
              exit
            end if
            integrand_gl(k) = 1.0_rp/(n+h*dndh) ! = 1 / d(nh)/dh
          end do ! k
          ! And Finally - define the refraction correction:

          ref_corr(j) = dot_product(integrand_gl,gw) * ym / Del_s(j)
        end if

        ! If ref_corr(j) is unphysical, use the trapezoidal rule instead

        if ( ref_corr(j) < 1.0 .or. ref_corr(j) > 1.3 ) then
          t1 = 1.0 / (1.0 + n1*(1.0 + eps*h1))
          t2 = 1.0 / (1.0 + n2*(1.0 + eps*h2))
          if ( t1 * t2 > 0.0 ) then
            ! Use trapezoidal rule on transformed integral
            ref_corr(j) = ym * ( t1 + t2 ) / Del_s(j)
          else
            ! Derivative has a discontinuity in the interval, so
            ! the transformation doesn't work.  Use the trapezoidal
            ! rule on the untransformed integral.
            t1 = (1.0 + n1) * h1
            t2 = (1.0 + n2) * h2
            ref_corr(j) = 0.5 * abs( ( h2 - h1 ) * &
              & (t1/sqrt(t1**2-nt2ht2) + t2/sqrt(t2**2-nt2ht2)) ) / Del_s(j)
          end if
        end if

        bad = bad .or. ref_corr(j) < 1.0 .or. ref_corr(j) > 1.3

      end do jl ! j

      j1 = tan_pt + 1
      j2 = no_ele - 1
    end do ! m

    if ( bad ) then
      ! Things are still haywire.  For unphysical ref_corr replace by the
      ! average of adjacent panels, bounded by 1.0 ... 1.3.
      call MLSMessage ( MLSMSG_Warning, moduleName, 'Drastic Ref_Corr fixup needed' )
      if ( max(index(switches,'drcx'),index(switches,'DRCX')) /= 0 ) then
        call output ( tan_pt, before='Tan_Pt = ' )
        call output ( ht, before=', Ht = ', advance='yes' )
        call dump ( h_path, name='H_Path', format='(f14.4)' )
        call dump ( n_path-1.0, name='N_Path' )
        call dump ( ref_corr, name='Ref_Corr' )
      end if
      do j = 2, no_ele-1
        if ( ref_corr(j) < 1.0 .or. ref_corr(j) > 1.3 ) &
          & ref_corr(j) = 0.5 * ( max(1.0_rp,min(1.3_rp,ref_corr(j-1))) + &
                                & max(1.0_rp,min(1.3_rp,ref_corr(j+1))) )
      end do
      if ( max(index(switches,'drcx'),index(switches,'DRCX')) /= 0 ) then
        call dump ( ref_corr, name='Ref_Corr after fixup' )
        if ( index(switches,'DRCX') /= 0 ) stop
      end if
    end if

  contains

  ! .................................................... Solve_Hn  .....

  ! Solve the equation h*(1.0+N(h)) = N*H, where N(h) is an exponential:
  !    N(h) = n1*Exp(eps*(h-h1)), for h,
  ! using a Newton iteration

    subroutine Solve_Hn

      integer :: iter
      real(rp) :: E ! N1 * exp(eps * (h-h1))
      real(rp) :: h_old, f1, f2, dh, hpos, hneg
      logical :: Head ! Need debug stuff heading

      real(rp), parameter :: H_tol = 0.001 ! km

      character(LEN=*), parameter :: Msg1 = &
        & 'From Solve_Hn routine: Could not bracket the root'

      integer,  parameter :: Max_Iter = 20
      character(LEN=*), parameter :: Msg2 = &
        & 'From Solve_Hn routine: Did not converge within 20 iterations'

      character(len=*), parameter :: Msg3 = &
        & 'Solution is within range even though not initially bracketed'

      character(max(len(msg1),len(msg2))+10) :: Msg

      stat = 0 ! Assume it will work
      f1 = h1 * (1.0_rp + n1) - NH ! residual at H = H1
      f2 = h2 * (1.0_rp + n2) - NH ! residual at H = H2

      if ( f1*f2 > 0.0_rp ) then ! start in the middle
        stat = 1
        msg = msg1
        write ( msg(len(msg1)+1:), '(" at ", i0)' ) j
        call MLSMessage ( MLSMSG_Warning, ModuleName, trim(msg) )
        call dumpDiags ( NH, f1, f2 )
        h = h1 + 0.5*(1.0+gx(k)) * ( h2 - h1 )
!        h = 0.5 * ( h1 + h2 )
      else ! start at the mean of the ends weighted by the other's residuals
        h = (h1 * abs(f2) + h2 * abs(f1)) / (abs(f1) + abs(f2))
      end if
      hneg = min(h1,h2)
      hpos = max(h1,h2)

      iter = 1

      if ( hndp > 0 ) then
        write ( *, '(7(3x,a2,9x)/3g14.7,1p,3g14.6,g20.12)' ) &
          & 'h1', 'h2', 'h ', 'f1', 'f2', 'Q ', 'NH', &
          &  h1,   h2,   h,    f1,   f2,  q,     NH
        head = .true. ! Need the heading for each iteration
      end if
      
      do

        h_old = h

        e = n1 * exp(eps*(h-h1))
        f2 = h * (1.0_rp + e ) - NH
        if ( abs(f2) < tiny ) exit ! Are we near a zero?

!         This is probably a bad idea in that it prevents correcting from
!         overshooting
!         if ( f2 < 0.0_rp ) then
!           hneg = h
!         else
!           hpos = h
!         end if

        dh = f2 / ( 1.0_rp + e * ( 1.0_rp + eps * h ) ) ! f2 / (d f2 / dh)

        h = h_old - dh ! Take the Newton move

        if ( hndp > 0 ) then
          if ( head ) &
            & write ( *, '(6(3x,a2,9x))' ) 'H-','H+','H ','dH','e ','f2'
          head = .false.
          write ( *, '(3g14.7,1p,3g14.6)' ) hneg, hpos, h, dh, e, f2
        end if

        if ( abs(dh) < htol ) exit ! Is the Newton move tiny?

        if ( (h-hneg)*(h-hpos) > 0.0 ) then
!         if ( h < min(hpos,hneg) .OR. h > max(hpos,hneg) ) &
!           h = 0.5_rp * (hneg + hpos) ! Keep H in bounds
          if ( h > hpos ) then
            if ( f2 < 0.0 ) then
              h = hneg
            else
              h = hpos
            end if
          else
            if ( f2 > 0.0 ) then
              h = hpos
            else
              h = hneg
            end if
          end if
!           h = min(hpos,max(hneg,h))
          if ( hndp > 0 ) write ( *, '(28x,1pg14.7,44x,"H was out of bounds")' ) h
          if ( abs(hpos-hneg) <= h_tol ) exit
        end if

        if ( abs(hpos-hneg) <= h_tol ) exit

        if ( iter >= max_iter ) then
          stat = 2
          msg = msg2
          write ( msg(len(msg2)+1:), '(" at ", i0)' ) j
          call MLSMessage ( MLSMSG_Warning, ModuleName, trim(msg) )
          call dumpDiags ( NH, f1, f2 )
          exit
        end if

        iter = iter + 1

      end do

      if ( stat == 1 ) then
        ! Maybe the solution is good even though we couldn't initially
        ! bracket it
        if ( (h-h1)*(h-h2) <= 0.0 ) then
          stat = 0
          call MLSMessage ( MLSMSG_Warning, ModuleName, msg3 )
        end if
      end if

      N = 1.0_rp + e
      dndh = eps * e

    end subroutine Solve_Hn

    subroutine DumpDiags ( NH, F1, F2 )

      real(rp), intent(in), optional :: NH, F1, F2

      if ( max(index(switches,'drfc'),index(switches,'DRFC')) == 0 ) return
      call output ( tan_pt, before='Tan_PT = ' )
      call output ( ht, before=', Ht = ' )
      if ( present(NH) ) then
        call output ( stat, before=', stat = ', advance='yes' )
        call output ( n1, before='n1 = ' )
        call output ( n2, before=', n2 = ' )
        call output ( h1, format='(f8.3)', before=', h1 = ' )
        call output ( h2, format='(f8.3)', before=', h2 = ', advance='yes' )
        call output ( nh, format='(f8.3)', before='NH = ' )
        call output ( f1, before=', f1 = ' )
        call output ( f2, before=', f2 = ', advance='yes' )
      else
        call output ( bad, before=', bad = ', advance='yes' )
        call dump ( h_path, name='H_Path' )
        call dump ( n_path-1.0, name='N_Path' )
      end if
      if ( index(switches,'DRFC') > 0 ) stop
    end subroutine DumpDiags

  end subroutine Comp_refcor

!------------------------------------------------------------------

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module REFRACTION_M

! $Log$
! Revision 2.30  2007/07/27 00:17:41  vsnyder
! Print enough to run Comp_Refcor off line if "Drastic correction" message
! is produced, and switch drcx or DRCX is set.  Stop after printing if
! DRCX is set.
!
! Revision 2.29  2007/07/11 22:27:39  vsnyder
! More robust integration in Comp_Refcor
!
! Revision 2.28  2007/02/02 00:22:44  vsnyder
! Don't bracket the Newton move so tightly
!
! Revision 2.27  2007/02/01 02:51:07  vsnyder
! Improve Newton iteration in Solve_HN
!
! Revision 2.26  2006/12/13 02:32:03  vsnyder
! Drag the tangent point around instead of assuming it's the middle one
!
! Revision 2.25  2006/06/29 19:31:59  vsnyder
! Use entire integration formula, not just interior points, in case of Lobatto
!
! Revision 2.24  2005/12/22 20:58:22  vsnyder
! Added more ranks and H2O update
!
! Revision 2.23  2005/12/07 00:33:48  vsnyder
! Update references to ATBD in comments
!
! Revision 2.22  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.21  2004/09/01 01:47:56  vsnyder
! Add status argument
!
! Revision 2.20  2004/06/17 00:08:24  vsnyder
! Removed two unused variables
!
! Revision 2.19  2003/11/04 01:55:08  vsnyder
! simplify nonconverged case
!
! Revision 2.18  2003/11/03 23:15:15  vsnyder
! Get rid of path_ds_dh procedure -- a one-liner used in one place
!
! Revision 2.17  2003/09/26 18:23:34  vsnyder
! Reinstate a lost CVS comment
!
! Revision 2.16  2003/09/17 23:33:26  vsnyder
! Major revision
!
! Revision 2.15  2003/09/17 23:32:56  vsnyder
! Clean up a few loose ends before major revision
!
! Revision 2.14  2003/09/13 02:02:00  vsnyder
! Converges faster with derivatives instead of differences
!
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
