! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module MieTheory

! -------------------------------------------------------------------------
! COMPUTE MIE EFFICIENCIES
! -------------------------------------------------------------------------

   use MLSKinds, only: r8
   USE MLSMessageModule, only: MLSMessage, MLSMSG_Error
   use RefractiveIndex, only: UKSUB, UKISUB
   implicit none

   private
   public :: MieCoeff

   ! These are legal values of ispi
   integer, parameter, public :: ice = 1
   integer, parameter, public :: water = ice + 1

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine MieCoeff ( ISPI, f, t, nr, r, a, b, nab, nabr, bc )

    use constants, only: Pi

   ! Arguments
    integer, intent(in) :: ISPI            ! cloud type (1=ICE, 2=WATER)
    real(r8), intent(in) :: f              ! frequency in GHz
    real(r8), intent(in) :: t              ! Temperature in K
    integer, intent(in) :: nr              ! no of particle size
    integer, intent(in) :: nab             ! no of a/b terms
    real(r8), intent(in) :: r(nr)          ! particle radius
    integer, intent(out) :: nabr(nr)       ! truncation number for a and b
                                           ! 10um ---> 5; 2000 um ---> 20.
    complex(r8), intent(out) :: a(nr,nab),b(nr,nab) ! Mie coefficients
    real(r8), intent(out) :: bc(3,nr)      ! single particle Mie
   ! Internal variables                          efficiencies (abs,scat,ext)
    real(r8) :: wl                         ! wavelength in meters
    complex(r8) :: m                       ! refractive index
    complex(r8) :: a0,a1,w9,w0,w1,p1,p2,mx,mx1
    real(r8) :: x, x1
    real(r8) :: ab_err                     ! relative error in Mie efficiency
    parameter(ab_err = 1.e-3)
    real(r8) :: dab1, dab2                 ! a/b contribution to bc at i term
    integer :: i,j

!... initialization
    wl=0.3_r8/f
    a = 0.0
    b = 0.0

  part_size_loop: do j = 1, nr

      if ( ISPI == ice ) then
        call ukisub ( f, t, m )
      else if ( ISPI == water ) then
        call uksub ( f, t, m )
      else
        CALL MLSMessage ( MLSMSG_Error, ModuleName, &
          & ' Unrecognized ispi parameter--ice==1, water==2 ' )
      end if

      x=real(2.*pi*r(j)/wl)*1.d-6

      bc(2,j) =0.
      bc(3,j) =0.

      !... default of nabr is nab, the largest possible
      nabr(j)=nab

      !... x1 is 1/x
      x1=real(0.5*wl/pi/r(j))*1.d6
      mx=cmplx(m)*x
      mx1=1.0/mx

      ! ....... for CGI     !JJ
      !        w9=dcmplx(dcos(x),-dsin(x))
      !        w0=dcmplx(dsin(x),dcos(x))
      !        a0=cdcos(mx)/cdsin(mx)

      ! ....   for f95
      w9=cmplx(cos(x),-sin(x),kind=r8)
      w0=cmplx(sin(x),cos(x),kind=r8)
      a0=cos(mx)/sin(mx)

      do i=1, nab

        w1=(2.0*i-1.)*x1*w0-w9
        a1=-i*mx1+1.0/(i*mx1-a0)
        p1=a1/m+i*x1
        p2=m*a1+i*x1

        a(j,i) = cmplx( (p1*real(w1)-real(w0))/(p1*w1-w0) )
        b(j,i) = cmplx( (p2*real(w1)-real(w0))/(p2*w1-w0) )

        ! ... the factor of 2./x^2 is left out in dab since nabr is usually
        !  important for large x. For x<0.1, usually, nabr=2 is enough.
             dab1=(2.*i+1.)*((abs(a(j,i)))**2+(abs(b(j,i)))**2)
             dab2=(2.*i+1.)*real(a(j,i)+b(j,i))
             bc(2,j)=bc(2,j)+dab1
             bc(3,j)=bc(3,j)+dab2

        ! ... determine cutoff no. for higher order terms
             if ( i >= 2 ) then
               dab1=dab1/bc(2,j)
               dab2=dab2/bc(3,j)
               if ( abs(dab1) < ab_err .and. abs(dab2) < ab_err ) then
                  nabr(j) = i
                  exit ! goto 11     ! stop looping and use i as nabr
               end if
             end if

        a0=a1
        w9=w0
        w0=w1

      end do
      ! 11   bc(2,j)=bc(2,j)*2/x/x
      bc(2,j)=bc(2,j)*2/x/x
      bc(3,j)=bc(3,j)*2/x/x

    end do part_size_loop

  end subroutine MieCoeff

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MieTheory

! $Log$
! Revision 2.9  2013/06/12 02:20:35  vsnyder
! Cruft removal
!
! Revision 2.8  2010/02/04 23:09:28  vsnyder
! Use kind= in CMPLX
!
! Revision 2.7  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.6  2007/10/06 00:00:17  vsnyder
! Modernization, polishing
!
! Revision 2.5  2007/06/21 00:52:54  vsnyder
! Remove tabs, which are not part of the Fortran standard
!
! Revision 2.4  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.3  2003/10/09 16:24:05  jonathan
! some cleaning
!
! Revision 2.2  2003/10/09 16:21:48  jonathan
! some changes
!
! Revision 2.1  2003/01/31 18:35:55  jonathan
! moved from cldfwm
!
! Revision 1.5  2003/01/30 22:01:14  pwagner
! Cosmetic changes
!
! Revision 1.4  2002/10/08 17:08:07  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 1.3  2001/09/21 15:51:37  jonathan
! modified F95 version
!
