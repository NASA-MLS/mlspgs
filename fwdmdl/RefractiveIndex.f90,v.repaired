! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module RefractiveIndex

!==================================================================
! COMPUTE COMPLEX REFRACTIVE INDEX OF WATER DROPS AND ICE PARTICLES
!
!* T  ---- TEMPERATURE IN  K                          ++ INPUT    *
!* F  ---- FREQUENCY  IN  GHZ                         ++ INPUT    *
!* E  ---- REFRACTIVE INDEX (E=M**2)                  ++ OUTPUT   *
!==================================================================

  implicit none

  public :: UKSUB     ! For Water
  public :: UKISUB    ! For Ice
  public :: UKISUB_dT ! For Ice, with temperature derivative

  interface UKISUB_dT
    module procedure UKISUB_dT_, UKISUB_dT_2
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  pure subroutine UKSUB ( F, T, M )  ! For Water
    !-----------------------
    use MLSKinds, only: r8
    real(r8), intent(in) :: F, T         ! GHz, K
    complex(r8), intent(out) :: M
    real(r8) :: th, fp, fs, e0, x, y     ! Liebe
    real(r8), parameter :: e1=5.48_r8    ! Liebe
    real(r8), parameter :: e2=3.51_r8    ! Liebe
!    real(r8) :: A, BV, C, E1_Lu, O0, O1, O2, R, S, Wl, X1, X2, XS, XC ! Lu

    th = 300.0/t - 1.0

!... from Liebe (1989)

    fp = 20.09 - (142.4 - 294*th)*th
    fs = 590.0 - 1500.0*th
    e0 = 77.66 + 103.3*th
    x = (e0-e1)/(1+(f/fp)**2) + (e1-e2)/(1+(f/fs)**2) + e2
    y = -( (e0-e1)*f/(fp*(1+(f/fp)**2)) + (e1-e2)*f/(fs*(1+(f/fs)**2)) )

!... Lu
!    wl=30./f
!    C=max(max(T,233.15)-273.16,-50.0)
!    E0=5.27137+(0.0216474-0.00131198*C)*C
!    A=-16.8129/(C+273.)+0.0609265
!    S=12.5664E+08
!    R=0.00033836*EXP(2513.98/(C+273.))
!    e1_Lu=78.54*(1.0-4.579E-3*(C-25.)+1.19E-5*(C-25.)**2-2.8E-8*(C-25.)**3)
!    X1=(R/WL)**(1.0-A)
!    XX=X1*X1
!    X2=A*PAI/2
!    XS=ASIN(X2)
!    XC=ACOS(X2)
!    BV=18.8496E+10
!    O1=(e1_Lu-E0)*(1.0+X1*XS)
!    O2=(e1_Lu-E0)*X1*XC
!    O0=1.0+2.0*X1*XS+XX
!    X=E0+O1/O0
!    Y=-( S*WL/BV+O2/O0 )
    m=sqrt(cmplx(x,y,kind=r8))

  end subroutine UKSUB

  pure subroutine UKISUB ( F, T, M )  ! For Ice
    ! ----------------------
    use Constants, only: Pi
    use MLSKinds, only: r8
    real(r8), intent(in) :: F, T         ! GHz, K
    complex(r8), intent(out) :: M
    real(r8) :: a, b, th, x, y           ! Hufford
!    real(r8), parameter :: c=1.16e-11_r8 ! Mishima
!    real(r8), parameter :: c = 6.543036549e-12 ! 1.11e-6/(30.0**3*2*pi) ! Dong Wu
    real(r8), parameter :: c = 1.11e-6/(30.0**3*2*pi) ! Dong Wu

    th = 300.0/t

    x = 3.15_r8

!{ from George Hufford, \emph{A model for the complex permittivity if
!  ice at frequencies below 1 THz}, {\bf International Journal of
!  Infrared and Millimeter Waves 12}, 7 (1991):\\
!  $\beta(T) = (0.445 \times 10^{-4} + 0.00211 \times 10^{-4} T )
!            + 0.585 \times 10^{-4} / ( 1 - T / 29.1 )^2$ GHz$^{-1}$
!  (T in degrees Celsius) =
!             $0.00211 \times 10^{-4} T - 0.131 \times 10^{-4}
!            + \frac{0.04947229919}{(302.2061046 - T)^2}$ (T in degrees Kelvin).
!  In terms of $\theta = 300/T$
!  (T in degrees Kelvin):
!  $\beta(T) = \frac{0.633 \times 10^{-4}} \theta - 0.131 \times 10^{-4}
!            + \left(\frac{7.36 \times 10^{-4} \theta}
!                         {\theta - 0.9927} \right)^2$.

    a = (50.4e-4 + 62.0e-4*(th-1)) * exp(-22.1*(th-1.))
!    b = 0.633e-4/th-0.131e-4 + (7.36e-4*th/(th-0.9927))**2
    b = 0.00211e-4 * t - 0.131e-4 + 0.04947229919/(302.2061046 - t)**2
!    y = -( a/f+b*f )

!... additional term from Mishima, as derived by Jonathan Jiang

!    y = -( a/f + b*f + c*f**3 )

!... and then modified by Dong Wu (2007-12-05)

    y = -( a/f + f*(b + c*f*f*sqrt(x)) )

    m = sqrt(cmplx(x,y,kind=r8))

  end subroutine UKISUB

  pure subroutine UKISUB_dT_ ( F, T, M, dM_dT )  ! For Ice
    ! ----------------------
    use Constants, only: Pi
    use MLSKinds, only: r8
    real(r8), intent(in) :: F, T         ! GHz, K
    complex(r8), intent(out) :: M, dM_dT
    real(r8) :: a, b, th, x, y           ! Hufford
!    real(r8), parameter :: c=1.16e-11_r8 ! Mishima
!    real(r8), parameter :: c = 6.543036549e-12 ! 1.11e-6/(30.0**3*2*pi) ! Dong Wu
    real(r8), parameter :: c = 1.11e-6/(30.0**3*2*pi) ! Dong Wu

    real(r8) :: da_dT, db_dt, de_dt, eth ! exp(-22.1*(th-1.))

    th = 300.0/t

    x = 3.15_r8

!  from George Hufford, \emph{A model for the complex permittivity if
!  ice at frequencies below 1 THz}, {\bf International Journal of
!  Infrared and Millimeter Waves 12}, 7 (1991):\\
!  $\beta(T) = (0.445 \times 10^{-4} + 0.00211 \times 10^{-4} T )
!            + 0.585 \times 10^{-4} / ( 1 - T / 29.1 )^2$ GHz$^{-1}$
!  (T in degrees Celsius) =
!             $0.00211 \times 10^{-4} T - 0.131 \times 10^{-4}
!            + \frac{0.04947229919}{(302.2061046 - T)^2}$ (T in degrees Kelvin).
!  In terms of $\theta = 300/T$
!  (T in degrees Kelvin):
!  $\beta(T) = \frac{0.633 \times 10^{-4}} \theta - 0.131 \times 10^{-4}
!            + \left(\frac{7.36 \times 10^{-4} \theta}
!                         {\theta - 0.9927} \right)^2$.

!{ From this
! \begin{equation*}\begin{split}
! \frac{\partial \epsilon}{\partial T} =\,&
!  -i \left( \frac1f \frac{\partial a}{\partial T} +
!   f \frac{\partial b}{\partial T} \right) \text{ where}\\
! \frac{\partial a}{\partial T} =\,&
!  \left( 62\times 10^{-4} \exp(-22.1(\theta-1)) -
!  22.1(50.4 + 62(\theta-1)) 10^{-4} \exp(-22.1(\theta-1)) \right)
!  \frac{\partial \theta}{\partial T} \\
!  =\,& (0.031836-0.13702\, \theta)\, \exp(-22.1(\theta-1))
!    \frac{\partial \theta}{\partial T} \\
! \frac{\partial b}{\partial T} =\,&
!  2.11 \times 10^{-7} + 2 \frac{0.04947229919}{(302.2061046 - T)^3} \\
! \frac{\partial \theta}{\partial T} =\,& -\frac{300}{T^2} = -\frac\theta{T}
! \end{split}\end{equation*}
!%
! Since $m = \sqrt{\epsilon}$ and $\frac{\partial\epsilon}{\partial T}$ is
! pure imaginary,
! \begin{equation*}
!  \frac{\partial m}{\partial T} =
!   \frac1{2m}\, \frac{\partial \epsilon}{\partial T} =
!   {m^*}\, \frac1{2|m|^2}\, \frac{\partial \epsilon}{\partial T} =
!   -\frac1{2|m|^2} \Im \frac{\partial\epsilon}{\partial T}
!    ( \Im m + i \Re m)
! \end{equation*}

    eth = exp(-22.1*(th-1.))
    a = (50.4e-4 + 62.0e-4*(th-1)) * eth
!    b = 0.633e-4/th-0.131e-4 + (7.36e-4*th/(th-0.9927))**2
    b = 2.11e-7 * t - 0.131e-4 + 0.04947229919/(302.2061046 - t)**2
!    y = -( a/f+b*f )

!... additional term from Mishima, as derived by Jonathan Jiang

!    y = -( a/f + b*f + c*f**3 )

!... and then modified by Dong Wu (2007-12-05)

    y = -( a/f + f*(b + c*f*f*sqrt(x)) )

    m = sqrt(cmplx(x,y,kind=r8))

    da_dT = -( 0.031836 - 0.13702 * th ) * th/T * eth
    db_dT = 2.11e-7 + 0.09894459838 / ( 302.2061046 - t )**3
    de_dt = -( da_dt / f + f * db_dt ) ! pure imaginary
    dm_dt = ( 0.5 / m ) * cmplx(0.0_r8,de_dt,kind=r8)

  end subroutine UKISUB_dT_

!  pure &
  subroutine UKISUB_dT_2 ( F, T, M, dM_dT, Y, dE_dT )  ! For Ice
    ! ----------------------
    use Constants, only: Pi
    use MLSKinds, only: r8
    real(r8), intent(in) :: F, T         ! GHz, K
    complex(r8), intent(out) :: M, dM_dT
    real(r8), intent(out) :: Y, dE_dT
    real(r8) :: a, b, th, x              ! Hufford
!    real(r8), parameter :: c=1.16e-11_r8 ! Mishima
!    real(r8), parameter :: c = 6.543036549e-12 ! 1.11e-6/(30.0**3*2*pi) ! Dong Wu
    real(r8), parameter :: c = 1.11e-6/(30.0**3*2*pi) ! Dong Wu

    real(r8) :: da_dT, db_dt, eth ! exp(-22.1*(th-1.))

    th = 300.0/t

    x = 3.15_r8

    eth = exp(-22.1*(th-1.))
    a = (50.4e-4 + 62.0e-4*(th-1)) * eth
!    b = 0.633e-4/th-0.131e-4 + (7.36e-4*th/(th-0.9927))**2
    b = 2.11e-7 * t - 0.131e-4 + 0.04947229919/(302.2061046 - t)**2
!    y = -( a/f+b*f )

!... additional term from Mishima, as derived by Jonathan Jiang

!    y = -( a/f + b*f + c*f**3 )

!... and then modified by Dong Wu (2007-12-05)

    y = -( a/f + f*(b + c*f*f*sqrt(x)) )

    m = sqrt(cmplx(x,y,kind=r8))
    da_dT = -( 0.031836 - 0.13702 * th ) * th/T * eth
    db_dT = 2.11e-7 + 0.09894459838 / ( 302.2061046 - t )**3
    de_dt = -( da_dt / f + f * db_dt ) ! pure imaginary
    dm_dt = ( 0.5 / m ) * cmplx(0.0_r8,de_dt,kind=r8)

  end subroutine UKISUB_dT_2

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module RefractiveIndex

! $Log$
! Revision 2.8  2010/02/04 23:09:28  vsnyder
! Use kind= in CMPLX
!
! Revision 2.7  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.6  2008/04/19 01:07:24  vsnyder
! Split UKISUB into state and derivatives, plus other minor revisions
!
! Revision 2.5  2008/02/29 01:59:11  vsnyder
! Added a sqrt(x) term
!
! Revision 2.4  2007/10/05 23:59:54  vsnyder
! Move water part of private COMX into UKSUB, ice part into UKISUB. then
! delete COMX.  Polish up formatting.
!
! Revision 2.3  2007/06/21 00:52:54  vsnyder
! Remove tabs, which are not part of the Fortran standard
!
! Revision 2.2  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.1  2003/01/31 18:35:55  jonathan
! moved from cldfwm
!
! Revision 1.6  2002/10/08 17:08:08  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 1.5  2002/04/16 03:42:25  jonathan
! fix a bug
!
! Revision 1.4  2002/04/15 22:22:32  jonathan
! check bug
!
! Revision 1.3  2002/04/15 17:13:28  jonathan
! add mishima term
!
! Revision 1.2  2001/09/21 15:51:37  jonathan
! modified F95 version
!
