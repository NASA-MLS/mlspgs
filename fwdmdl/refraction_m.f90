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

  !{ Compute the refractive index
  !  $n = A \frac{P}{T} \left( 1 + B
  !           \frac{\mathbf{f}_{{\text{H}_2\text{O}}}}{T} \right)$
  !%
  !  and its derivatives w.r.t.\ temp and H2O.  Don't bother with
  !  derivatives w.r.t.\ $P$ or $\zeta$, which can be gotten easily without
  !  needing $A$ or $B$, and would require derivatives of $T$ and
  !  $\mathbf{f}_{{\text{H}_2\text{O}}}$ with respect to $P$ or $\zeta$.\\
  !
  !  $\frac{\text{d} n}{\text{d} T} = 
  !    \frac1T \left( \frac{A P}T - 2 n \right)$ or
  !    $- \frac{n}T$ if $\mathbf{f}_{{\text{H}_2\text{O}}}$ is not used.
  !\\[5pt]
  !  $\frac{\text{d} n}{\text{d} \mathbf{f}_{{\text{H}_2\text{O}}}} =
  !    \frac{A B P}{T^2}$
  !\\[5pt]
  !  $\frac{\text{d} n}{\text{d} P} = \frac{n}{P} +
  !    \frac{\text{d} n}{\text{d} T} \frac{\text{d} T}{\text{d} P} +
  !    \frac{\text{d} n}{\text{d} \mathbf{f}_{{\text{H}_2\text{O}}}}
  !    \frac{\text{d} \mathbf{f}_{{\text{H}_2\text{O}}}}{\text{d} P}$
  !\\[5pt]
  !  $\frac{\text{d} n}{\text{d} \zeta} = -n \ln 10 +
  !    \frac{\text{d} n}{\text{d} T} \frac{\text{d} T}{\text{d} \zeta} +
  !    \frac{\text{d} \mathbf{f}_{{\text{H}_2\text{O}}}}{\text{d} \zeta} =
  !  \frac{\text{d} n}{\text{d} P} \frac{\text{d} P}{\text{d} \zeta}$

  use MLSKINDS, only: RP, R8
    
  implicit none

  private
  public :: Refractive_Index, Refractive_Index_Deriv, Refractive_Index_f, &
          & Refractive_Index_H2O_Update, Comp_Refcor
  public :: MaxRefraction

  real(rp), parameter, public :: RefrAterm = 0.0000776_rp
  real(rp), parameter, public :: RefrBterm = 4810.0_rp

  interface Refractive_Index
    module procedure Refractive_Index_0,   Refractive_Index_0_h2o
    module procedure Refractive_Index_1,   Refractive_Index_1_h2o
    module procedure Refractive_Index_1_2, Refractive_Index_1_2_h2o
    module procedure Refractive_Index_2,   Refractive_Index_2_h2o
  end interface

  interface Refractive_Index_f
    module procedure Refractive_Index_0_f, Refractive_Index_0_h2o_f
  end interface

  interface Refractive_Index_deriv
    module procedure Refractive_Index_0_deriv
    module procedure Refractive_Index_0_h2o_deriv
    module procedure Refractive_Index_1_deriv
    module procedure Refractive_Index_1_h2o_deriv
  end interface

  interface Refractive_Index_H2O_Update
    module procedure Refractive_Index_H2O_Update_0
    module procedure Refractive_Index_H2O_Update_1
    module procedure Refractive_Index_H2O_Update_2
  end interface

  ! This is the maximum amount of refraction allowed
  real(r8), parameter :: MaxRefraction = 0.0004 ! Add one to get refractive index

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

!------------------------------------------  Refractive_Index_0_f  -----
  pure elemental function Refractive_Index_0_f ( p_path, t_path ) result ( n_path )

  ! This routine computes the refractive index as a function of pressure
  ! and temperature. The returned value has one subtracted from it.
  ! We could easily make this elemental.

  ! inputs
    real(rp), intent(in) :: p_path  ! pressure(hPa)
    real(rp), intent(in) :: t_path  ! temperature(K)
  ! output
    real(rp) :: n_path ! refractive indicies - 1

  ! begin code
    n_path = refrAterm * p_path / t_path

  end function Refractive_Index_0_f

!--------------------------------------------  Refractive_Index_0  -----
  subroutine Refractive_Index_0 ( p_path, t_path, n_path )

  ! This routine computes the refractive index as a function of pressure
  ! and temperature. The returned value has one subtracted from it.
  ! We could easily make this elemental.

  ! inputs
    real(rp), intent(in) :: p_path  ! pressure(hPa)
    real(rp), intent(in) :: t_path  ! temperature(K)
  ! output
    real(rp), intent(out) :: n_path ! refractive indicies - 1

  ! begin code
    n_path = refrAterm * p_path / t_path

  end subroutine Refractive_Index_0

!--------------------------------------  Refractive_Index_0_deriv  -----
  subroutine Refractive_Index_0_deriv ( p_path, t_path, n_path, dNdT )

  ! This routine computes the refractive index as a function of pressure
  ! and temperature. The returned value has one subtracted from it.
  ! We could easily make this elemental.

  ! inputs
    real(rp), intent(in) :: p_path  ! pressure(hPa)
    real(rp), intent(in) :: t_path  ! temperature(K)
  ! output
    real(rp), intent(out) :: n_path ! refractive indicies - 1
    real(rp), intent(out) :: dNdT   ! d(n_path)/d(t_path)

  ! begin code
    n_path = refrAterm * p_path / t_path

    dNdT = - n_path / t_path

  end subroutine Refractive_Index_0_deriv

!----------------------------------------  Refractive_Index_0_h2o_f  -----
  pure elemental function Refractive_Index_0_h2o_f ( p_path, t_path, h2o_path ) &
                                                   & result ( n_path )
  ! This routine computes the refractive index as a function of pressure,
  ! temperature, and H2O. The returned value has one subtracted from it.
  ! We could easily make this elemental.

  ! inputs
    real(rp), intent(in) :: p_path   ! pressure(hPa)
    real(rp), intent(in) :: t_path   ! temperature(K)
  ! Keywords
    real(rp), intent(in) :: h2o_path ! H2O vmr(ppv)
  ! output
    real(rp) :: n_path ! refractive indicies - 1

  ! begin code
    n_path = refrAterm * p_path / t_path * ( 1.0_rp + refrBterm*h2o_path/t_path)

  end function Refractive_Index_0_h2o_f

!----------------------------------------  Refractive_Index_0_h2o  -----
  subroutine Refractive_Index_0_h2o ( p_path, t_path, n_path, h2o_path )

  ! This routine computes the refractive index as a function of pressure,
  ! temperature, and H2O. The returned value has one subtracted from it.
  ! We could easily make this elemental.

  ! inputs
    real(rp), intent(in) :: p_path   ! pressure(hPa)
    real(rp), intent(in) :: t_path   ! temperature(K)
  ! output
    real(rp), intent(out) :: n_path  ! refractive indicies - 1
  ! Keywords
    real(rp), intent(in) :: h2o_path ! H2O vmr(ppv)

  ! begin code
    n_path = refrAterm * p_path / t_path * ( 1.0_rp + refrBterm*h2o_path/t_path)

  end subroutine Refractive_Index_0_h2o

!-----------------------------------  Refractive_Index_0_h2o_deriv  -----
  subroutine Refractive_Index_0_h2o_deriv ( p_path, t_path, h2o_path, &
                                          & n_path, dNdT, dNdH )

  ! This routine computes the refractive index as a function of pressure,
  ! temperature, and H2O. The returned value has one subtracted from it.
  ! We could easily make this elemental.  Its derivatives w.r.t. Temperature
  ! and H2O are also computed.

  ! inputs
    real(rp), intent(in) :: p_path   ! pressure(hPa)
    real(rp), intent(in) :: t_path   ! temperature(K)
    real(rp), intent(in) :: h2o_path ! H2O vmr(ppv)
  ! output
    real(rp), intent(out) :: n_path  ! refractive indicies - 1
    real(rp), intent(out) :: dNdT    ! d(n_path)/d(t_path)
    real(rp), intent(out) :: dNdH    ! d(n_path)/d(h2o_path)

    real(rp) :: APT                  ! A * P / T

  ! begin code
    apt = refrAterm * p_path / t_path
    n_path = apt * ( 1.0_rp + refrBterm*h2o_path/t_path)

    dNdT = ( apt - 2.0*n_path ) / t_path
    dNdH = apt * refrBterm / t_path  ! A B P /T^2

  end subroutine Refractive_Index_0_h2o_deriv

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

! -------------------------------------------  Refractive_Index_1  -----
  subroutine Refractive_Index_1 ( p_path, t_path, n_path )

  ! This routine computes the refractive index as a function of pressure
  ! and temperature. The returned value has one subtracted from it.

  ! inputs
    real(rp), intent(in) :: p_path(:)  ! pressure(hPa)
    real(rp), intent(in) :: t_path(:)  ! temperature(K)
  ! output
    real(rp), intent(out) :: n_path(:) ! refractive indicies - 1

  ! begin code
    n_path = refrAterm * p_path / t_path

  end subroutine Refractive_Index_1

!--------------------------------------  Refractive_Index_1_deriv  -----
  subroutine Refractive_Index_1_deriv ( p_path, t_path, n_path, dNdT )

  ! This routine computes the refractive index as a function of pressure
  ! and temperature. The returned value has one subtracted from it.
  ! Its derivative w.r.t. Temperature is also computed.

  ! inputs
    real(rp), intent(in) :: p_path(:)  ! pressure(hPa)
    real(rp), intent(in) :: t_path(:)  ! temperature(K)
  ! output
    real(rp), intent(out) :: n_path(:) ! refractive indicies - 1
    real(rp), intent(out) :: dNdT(:)   ! d(n_path)/d(t_path)
  ! optional

  ! begin code
    n_path = refrAterm * p_path / t_path

    dNdT = - n_path / t_path

  end subroutine Refractive_Index_1_deriv

! ---------------------------------------  Refractive_Index_1_h2o  -----
  subroutine Refractive_Index_1_h2o ( p_path, t_path, n_path, h2o_path )

  ! This routine computes the refractive index as a function of pressure,
  ! temperature, and H2O. The returned value has one subtracted from it.
  ! Its derivatives w.r.t. Temperature and H2O are also computed.

  ! inputs
    real(rp), intent(in) :: p_path(:)   ! pressure(hPa)
    real(rp), intent(in) :: t_path(:)   ! temperature(K)
  ! output
    real(rp), intent(out) :: n_path(:)  ! refractive indicies - 1
  ! Keywords
    real(rp), intent(in) :: h2o_path(:) ! H2O vmr(ppv)

  ! begin code
    n_path = refrAterm * p_path / t_path * ( 1.0_rp + refrBterm*h2o_path/t_path)

  end subroutine Refractive_Index_1_h2o

!-----------------------------------  Refractive_Index_1_h2o_deriv  -----
  subroutine Refractive_Index_1_h2o_deriv ( p_path, t_path, h2o_path, &
                                          & n_path, dNdT, dNdH )

  ! This routine computes the refractive index as a function of altitude
  ! and phi. The returned value has one subtracted from it.
  ! We could easily make this elemental.

  ! inputs
    real(rp), intent(in) :: p_path(:)    ! pressure(hPa)
    real(rp), intent(in) :: t_path(:)    ! temperature(K)
    real(rp), intent(in) :: h2o_path(:)  ! H2O vmr(ppv)
  ! output
    real(rp), intent(out) :: n_path(:)   ! refractive indicies - 1
    real(rp), intent(out) :: dNdT(:)     ! d(n_path)/d(t_path)
    real(rp), intent(out) :: dNdH(:)     ! d(n_path)/d(h2o_path)

    real(rp) :: APT(size(p_path))        ! A * P / T

  ! begin code
    apt = refrAterm * p_path / t_path
    n_path = apt * ( 1.0_rp + refrBterm*h2o_path/t_path)

    dNdT = ( apt - 2.0*n_path ) / t_path
    dNdH = apt * refrBterm / t_path      ! A B P /T^2

  end subroutine Refractive_Index_1_h2o_deriv

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

! -----------------------------------------  Refractive_Index_1_2  -----
  subroutine Refractive_Index_1_2 ( p_path, t_path, n_path )

  ! This routine computes the refractive index as a function of altitude
  ! and phi. The returned value has one subtracted from it.

  ! inputs
    real(rp), intent(in) :: p_path(:)    ! pressure(hPa)
    real(rp), intent(in) :: t_path(:,:)  ! temperature(K)
  ! output
    real(rp), intent(out) :: n_path(:,:) ! refractive indicies - 1

    integer :: I

  ! begin code
    do i = 1, size(t_path,2)
      n_path(:,i) = refrAterm * p_path(:) / t_path(:,i)
    end do

  end subroutine Refractive_Index_1_2

! -------------------------------------  Refractive_Index_1_2_h2o  -----
  subroutine Refractive_Index_1_2_h2o ( p_path, t_path, n_path, h2o_path )

  ! This routine computes the refractive index as a function of altitude
  ! and phi. The returned value has one subtracted from it.

  ! inputs
    real(rp), intent(in) :: p_path(:)     ! pressure(hPa)
    real(rp), intent(in) :: t_path(:,:)   ! temperature(K)
  ! output
    real(rp), intent(out) :: n_path(:,:)  ! refractive indicies - 1
  ! Keywords
    real(rp), intent(in) :: h2o_path(:,:) ! H2O vmr(ppv)

    integer :: I

  ! begin code
    do i = 1, size(t_path,2)
      n_path(:,i) = refrAterm * p_path(:) / t_path(:,i) * &
        & ( 1.0_rp + refrBterm*h2o_path(:,i)/t_path(:,i))
    end do

  end subroutine Refractive_Index_1_2_h2o

! -------------------------------------------  Refractive_Index_2  -----
  subroutine Refractive_Index_2 ( p_path, t_path, n_path )

  ! This routine computes the refractive index as a function of altitude
  ! and phi. The returned value has one subtracted from it.

  ! inputs
    real(rp), intent(in) :: p_path(:,:)  ! pressure(hPa)
    real(rp), intent(in) :: t_path(:,:)  ! temperature(K)
  ! output
    real(rp), intent(out) :: n_path(:,:) ! refractive indicies - 1

  ! begin code
    n_path = refrAterm * p_path / t_path

  end subroutine Refractive_Index_2

! ---------------------------------------  Refractive_Index_2_h2o  -----
  subroutine Refractive_Index_2_h2o ( p_path, t_path, n_path, h2o_path )

  ! This routine computes the refractive index as a function of altitude
  ! and phi. The returned value has one subtracted from it.

  ! inputs
    real(rp), intent(in) :: p_path(:,:)   ! pressure(hPa)
    real(rp), intent(in) :: t_path(:,:)   ! temperature(K)
  ! output
    real(rp), intent(out) :: n_path(:,:)  ! refractive indicies - 1
  ! Keywords
    real(rp), intent(in) :: h2o_path(:,:) ! H2O vmr(ppv)

  ! begin code
    n_path = refrAterm * p_path / t_path * ( 1.0_rp + refrBterm*h2o_path/t_path)

  end subroutine Refractive_Index_2_h2o

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

! --------------------------------------------------  Comp_Refcor  -----

  subroutine Comp_Refcor ( tan_pt, h_path, n_path, ht, del_s, ref_corr, &
                         & status )

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
  ! If the refractive correction is not in the interval [1.0,rmax], or if
  ! $\frac{\text{d}}{\text{d}h} h N(h)$ changes sign, use the trapezoidal
  ! rule instead.  On the ends of the interval, $h=H_i$ and $E=n_i$, so
  ! $\mathcal{N}(H_i) = \mathcal{N}_i$, giving $I=\frac{x_i-x_{i-1}}2
  ! (f_{i-1}+f_i)$ where the integrand $f_i = \frac1{1+n_i(1+\epsilon H_i)}$.

  ! For derivation of the code below, please see: "FWD Model" paper,
  ! Page 16, Eqn. 26 & 27

    use DUMP_0, only: DUMP
    use GLNP, only: NG, GX=>GX_ALL, GW=>GW_ALL
    use MLSKINDS, only: RP, IP
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ERROR, MLSMSG_WARNING
    use MLSSTRINGLISTS, only: SWITCHDETAIL
    use OUTPUT_M, only: OUTPUT
    use TOGGLES, only: SWITCHES

    integer, intent(in) :: Tan_pt  ! Tangent point index in H_Path etc.
    real(rp), intent(in) :: H_PATH(:)
    real(rp), intent(in) :: N_PATH(:)

    real(rp), intent(in) :: ht     ! Geometric tangent height from Earth center, km

    real(rp), intent(out) :: Del_s(:) ! Del_s(j) = path length from j to j+m,
                                   ! where m is the direction from the tangent
                                   ! point, except Del_s(0) = Del_s(no_ele) = 0.
    real(rp), intent(out) :: REF_CORR(:)

    integer, intent(out) :: Status ! 0 = OK, 1 = too many iterations,
                                   ! 4 = failed to bracket root,
                                   ! 5 = discontinuity, 6 = improper usage

    logical :: Bad      ! Ref_Corr < 1 somewhere
    integer(ip) :: HNDP ! "Solve_Hn Detail Printing"
    integer(ip) :: RCFX ! "Ref_Cor FiXup printing"
    integer(ip) :: j, j1, j2, k, m, no_ele, panel_stat, stat, My_Tan

    real(rp) :: INTEGRAND_GL(Ng)

    real(rp) :: q, htan2, NH, Nt2Ht2
    real(rp) :: dndh, eps, H, h1, h2, N, n1, n2, t1, t2, x1, x2, xm, ym

    real(rp), parameter :: Tiny = 1.0e-8_rp
!   real(rp), parameter :: Rmax = 1.3_rp

    status = 0
    hndp = switchDetail(switches,'hndp') ! bit 1 => Dump the arrays if trouble
                                         ! bit 2 => Dump the arrays and stop
                                         ! bit 4 => Dump the iterates
    rcfx = switchDetail(switches,'rcfx')

    no_ele = size(n_path)
    my_tan = min(tan_pt, no_ele)

  !  Initialize the ref_corr array:

    ref_corr(1:no_ele) = 1.0_rp

    htan2 = ht * ht
    Nt2Ht2 = (n_path(my_tan)*ht)**2

    bad = .false.

    j1 = my_tan
    j2 = 2
    del_s(1) = 0.0
    do m = -1, 1, 2
      h2 = h_path(j1)
      n2 = n_path(j1)

      q = (h2 * n2)**2 - nt2ht2
      if ( abs(q) < tiny ) q = 0.0_rp
      n2 = n2 - 1.0_rp
      x2 = sqrt(abs(q))

      do j = j1, j2, m

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
          ! Where N is essentially constant, the integral is easy, and probably
          ! almost exactly 1.0.
          if ( hndp > 0 ) write ( *, '("N essentially constant for Ref_Corr(",i0,")")' ) j
          ref_corr(j) = abs( sqrt((n_path(j+m)*h2)**2-nt2ht2) - &
            &                sqrt((n_path(j  )*h1)**2-nt2ht2) ) / &
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
          if ( hndp > 3 ) write ( *, '("Solving for H for Ref_Corr(",i0,")")' ) j
          panel_stat = 0
          eps = log(n2/n1)/(h2-h1)
          xm = 0.5_rp *(x2 + x1)      ! Midpoint of the interval
          ym = 0.5_rp *(x2 - x1)      ! Half of the interval length
          do k = 1, ng
            q = xm + ym * gx(k)       ! Gauss abscissa
            NH = sqrt(q*q + nt2ht2)
            ! Solve h*(1+n1*exp(eps*(h-h1))) = sqrt(q*q + nt2ht2) for h
            call solve_hn
            panel_stat = max(panel_stat,stat)
            if ( stat > 1 ) exit
            integrand_gl(k) = 1.0_rp/(n+h*dndh) ! = 1 / d(nh)/dh
          end do ! k
          status = max(status,panel_stat)

          ! And Finally - define the refraction correction:

          if ( panel_stat < 2 ) then
            ref_corr(j) = dot_product(integrand_gl,gw) * ym / Del_s(j)
          else
            ref_corr(j) = -1.0 ! To be corrected later
          end if
        end if

        ! If ref_corr(j) is unphysical, use the trapezoidal rule instead

        bad = bad .or. ref_corr(j) < 1.0

      end do ! j

      if ( my_tan < tan_pt ) exit
      j1 = tan_pt + 1
      j2 = no_ele - 1
      del_s(no_ele) = 0.0
    end do ! m

    if ( bad ) then
      ! Things are unphysical somewhere.  Replace ref_corr by 1.0
      ! wherever it's less than 1.0.  Alternatively, the average
      ! of neighbors, or a value interpolated from neighbors could be used.
      if ( rcfx >= 0 ) call output( 'Ref_Corr fixup needed.', advance='yes' )
      call dumpArrays
      if ( rcfx >= 0 ) call output ( 'Ref_Corr fixed at' )
      do j = 2, no_ele-1
        if ( ref_corr(j) < 1.0 ) then
          ref_corr(j) = 1.0
!           ref_corr(j) = 0.5 * ( max(1.0_rp,min(rmax,ref_corr(j-1))) + &
!                      & max(1.0_rp,min(rmax,ref_corr(j+1))) )
          if ( rcfx >= 0 ) call output ( j, before=' ' )
        end if
      end do
      if ( rcfx >= 0 ) then
        call output ( '', advance='yes' )
        call dump ( ref_corr, name='Ref_Corr after fixup' )
        if ( rcfx > 0 ) stop
      end if
    end if

    if ( status > 0 ) then ! Solve_Hn didn't converge somewhere
      call dumpArrays
      if ( iand(hndp,3) /= 0 ) call dumpArrays
      if ( iand(hndp,2) /= 0 ) stop
    end if

  contains

  ! .................................................... Solve_Hn  .....

  ! Solve the equation h*(1.0+N(h)) = N*H, where N(h) is an exponential:
  !    N(h) = n1*Exp(eps*(h-h1)), for h,

    subroutine Solve_Hn

      ! Stat = 0 => all OK, = 1 => too many iterations,
      !      = 4 => Could not bracket root (no solution),
      !      = 5 => Discontinuity (no solution),
      !      = 6 => improper usage (ought to stop).
      use Zero_m, only: ZERO

      integer :: iter
      real(rp) :: E ! N1 * exp(eps * (h-h1))
      real(rp) :: f1, f2, x1, x2
      logical :: Head ! Need debug stuff heading

      real(rp), parameter :: H_tol = 0.001_rp ! km

      integer,  parameter :: Max_Iter = 20
      character(LEN=*), parameter :: Msg1 = &
        & 'From Solve_Hn routine: Did not converge within 20 iterations'

      character(LEN=*), parameter :: Msg2 = &
        & 'From Solve_Hn routine: Zero finder thinks F has a discontinuity'

      character(LEN=*), parameter :: Msg3 = &
        & 'From Solve_Hn routine: Could not bracket the root'

      character(LEN=*), parameter :: Msg4 = &
        & 'From Solve_Hn routine: Improper usage of zero finder'

      character(max(len(msg1),len(msg2),len(msg3))+10) :: Msg

      iter = 1

      x1 = h1
      x2 = h2
      e = n1 * exp(eps*(x2-h1))
      f2 = x2 * (1.0_rp + e ) - NH

      head = .true.
      stat = 0
      do while ( iter < max_iter )

        iter = iter + 1
        e = n1 * exp(eps*(x1-h1))
        f1 = x1 * (1.0_rp + e ) - NH
        if ( hndp > 3 ) then
          if ( head ) then
            write ( *, '(5(3x,a2,9x)/2g14.7,1p,2g14.6,g20.12)' ) &
            & 'x1', 'x2', 'f1', 'f2', 'NH', &
            &  h1,   h2,   f1,   f2,   NH
          else
            write ( *, '(g14.7,14x,1pg14.6)') x1, f1
          end if
          head = .false.
        end if
        call zero ( x1, f1, x2, f2, stat, h_tol )
        if ( stat /= 1 ) exit

      end do

      select case ( stat )
      case ( 1 ) ! Too many iterations
        msg = msg1
        write ( msg(len(msg1)+1:), '(" at ", i0)' ) j
        call MLSMessage ( MLSMSG_Warning, ModuleName, trim(msg) )
        call dumpDiags ( NH, f1, f2 )
      case ( 2, 3 ) ! Normal termination, tolerance is/is not satisfied
        stat = 0 ! Worked OK
        e = n1 * exp(eps*(x1-h1))
        h = x1
      case ( 4 ) ! F has a discontinuity -- no solution
        msg = msg2
        write ( msg(len(msg2)+1:), '(" at ", i0)' ) j
        call MLSMessage ( MLSMSG_Warning, ModuleName, trim(msg) )
        call dumpDiags ( NH, f1, f2 )
      case ( 5 ) ! Couldn't bracket root -- no solution
        msg = msg3
        write ( msg(len(msg3)+1:), '(" at ", i0)' ) j
        call MLSMessage ( MLSMSG_Warning, ModuleName, trim(msg) )
        call dumpDiags ( NH, f1, f2 )
      case ( 6 ) ! Improper usage
        msg = msg4
        write ( msg(len(msg4)+1:), '(" at ", i0)' ) j
        call MLSMessage ( MLSMSG_Error, ModuleName, trim(msg) )
      end select

      N = 1.0_rp + e
      dndh = eps * e

    end subroutine Solve_Hn

    subroutine DumpArrays

      if ( rcfx >= 0 ) then
        call output ( tan_pt, before='Tan_Pt = ' )
        call output ( ht, before=' Ht = ' )
        call output ( no_ele, before=' No_Ele = ' )
        call dump ( h_path, name=' H_Path', format='(f14.4)' )
        call dump ( n_path-1.0, name='N_Path' )
        call dump ( ref_corr, name='Ref_Corr' )
      end if

    end subroutine DumpArrays

    subroutine DumpDiags ( NH, F1, F2 )

      real(rp), intent(in), optional :: NH, F1, F2

      integer :: DFRC
      dfrc = switchDetail(switches,'drfc')
      if ( dfrc < 0 ) return
      call output ( tan_pt, before='Tan_PT = ' )
      call output ( ht, before=', Ht = ' )
      call output ( no_ele, before=' No_Ele = ' )
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
      if ( dfrc > 0 ) stop

    end subroutine DumpDiags

  end subroutine Comp_Refcor

!------------------------------------------------------------------

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module REFRACTION_M

! $Log$
! Revision 2.46  2016/01/23 02:52:55  vsnyder
! Add LaTeX, procedures to compute derivatives
!
! Revision 2.45  2013/07/26 22:19:05  vsnyder
! Fiddle with dump switches
!
! Revision 2.44  2013/06/12 02:33:19  vsnyder
! Cruft removal
!
! Revision 2.43  2013/02/28 21:05:48  vsnyder
! Try to cope with short paths
!
! Revision 2.42  2012/02/16 22:44:36  pwagner
! Skip printing Ref_Corr msg unless switch rcfx set
!
! Revision 2.41  2011/05/09 18:03:33  pwagner
! Converted to using switchDetail
!
! Revision 2.40  2011/01/05 00:22:08  vsnyder
! TeXnicalities
!
! Revision 2.39  2009/07/09 23:53:53  vsnyder
! Replace MLSCommon with MLSKinds
!
! Revision 2.38  2009/06/23 18:26:11  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.37  2009/06/13 01:08:51  vsnyder
! Improve comments about dummy arguments
!
! Revision 2.36  2009/01/16 23:42:04  vsnyder
! Don't need to set all of Del_s to zero initially because all but 1 and
! no_ele get set anyway.  Add PRINT statement to not_used_here.
!
! Revision 2.35  2007/10/30 20:52:26  vsnyder
! Change diagnostic output criteria
!
! Revision 2.34  2007/10/02 22:35:17  vsnyder
! Cannonball polishing
!
! Revision 2.33  2007/09/07 22:13:00  vsnyder
! Don't put an upper bound on the refractive correction
!
! Revision 2.32  2007/09/07 03:00:32  vsnyder
! Use Math 77 zero finder
!
! Revision 2.31  2007/07/31 23:48:10  vsnyder
! Integrate away from tangent point to handle singularity correctly
!
! Revision 2.30  2007/07/27 00:17:41  vsnyder
! Print enough to run Comp_refcor off line if "Drastic correction" message
! is produced, and switch drcx or DRCX is set.  Stop after printing if
! DRCX is set.
!
! Revision 2.29  2007/07/11 22:27:39  vsnyder
! More robust integration in Comp_refcor
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
! Publish two constants, make Refractive_index generic
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
! Make Comp_refcor more robust
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
