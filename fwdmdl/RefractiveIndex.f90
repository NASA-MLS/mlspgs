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

  public :: UKSUB  ! For Water
  public :: UKISUB ! For Ice

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine UKSUB ( F, T, M )  ! For Water
    !-----------------------
    use MLSKinds, only: r8
    real(r8), intent(in) :: F, T
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
    m=sqrt(cmplx(x,y))

  end subroutine UKSUB

  subroutine UKISUB ( F, T, M )  ! For Ice
    ! ----------------------
    use MLSKinds, only: r8
    real(r8), intent(in) :: F, T
    complex(r8), intent(out) :: M
    real(r8) :: a, b, th, x, y           ! Hufford
    real(r8), parameter :: c=1.16e-11_r8 ! Mishima

    th = 300.0/t

!... from Hufford (1991) model

    a = (50.4 + 62.*(th-1)) * 1.0e-4*exp(-22.1*(th-1.))
    b = (0.633/th-0.131)*1e-4 + (7.36e-4*th/(th-0.9927))**2
    x = 3.15
!    y = -( a/f+b*f )

!... additional term from Mishima

    y = -( a/f + b*f + c*f**3 )
!    write(*,*)y,f,t

    m=sqrt(cmplx(x,y))

  end subroutine UKISUB

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module RefractiveIndex

! $Log$
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
