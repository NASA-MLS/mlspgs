! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module SizeDistribution

USE MLSSpecialFunctions, ONLY: gamma

! -------------------------------------------------------------------------  
! THIS MODULE CONTAINS VARIOUS PARTICLE SIZE DISTRIBUTIONS
! -------------------------------------------------------------------------

  implicit NONE
  private
  public :: DRP_SIZE

  ! These are legal values of ISPI argument for DRP_SIZE
  integer, public, parameter :: ice = 1
  integer, public, parameter :: water = ice + 1

  ! These are legal values of IPSD argument for DRP_SIZE
  integer, public, parameter :: mh_sd = 1000
  integer, public, parameter :: gamma_sd = 2000
  integer, public, parameter :: liu_curry_sd = 1100
  integer, public, parameter :: knollenberg_sd = 4000

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
  private :: not_used_here 
 !---------------------------------------------------------------------------
      
contains

      SUBROUTINE DRP_SIZE(ISPI, R, RN, NR, CWC, T, IPSD, Dm)

!========================================================================C
!     PRODUCE SIZE DISTRIBUTIONS FOR ICE PARTICLES AND WATER DROPLETS    C
!     J.JIANG, MAY 1, 2001                                               C
!========================================================================C

      use MLSCommon, only: r8
      use MLSMessageModule, only: MLSMessage, MLSMSG_Error

      ! Arguments
      INTEGER, intent(in) :: ISPI            ! ICE or WATER?

      INTEGER, intent(in) :: NR              ! NUMBER OF SIZE BINS

      INTEGER, intent(in) :: IPSD            ! SIZE-DISTRIBUTION TYPE
                                             ! See parameters MH_SD etc. above.

      REAL(r8), intent(out) :: RN(NR)        ! INTEGRATED NUMBER DENSITY OF EACH 
                                             ! SIZE BIN
      REAL(r8), intent(in) :: R(NR)          ! PARTICLE RADIUS (micron)
      REAL(r8), intent(in) :: CWC            ! WATER CONTENT (g/m3)
      REAL(r8), intent(in) :: T              ! TEMPERATURE (K)
      REAL(r8), intent(out) :: Dm            ! MASS-MEAN-DIAMETER

      ! Internal variables
      REAL(r8) :: DR0                        ! INTEGRATION STEP FOR SIZE BINS
      REAL(r8) :: RN0                        ! INTEGRATION NUMBER DENSITY 
      REAL(r8) :: DIAM                       ! PARTICLE DIAMETER (micron)
      REAL(r8) :: TEMPC                      ! TEMPERATURE (C)
      REAL(r8) :: IWC0                       ! IWC OF SIZE < 100 microns
      REAL(r8) :: IWC1                       ! IWC OF SIZE > 100 microns
 
      REAL(r8), parameter :: RAOI = 0.91_r8  ! ICE DENSITY (g/m3)

      REAL(r8) :: IWC00                      ! NORMALIZATION IWC (g/m3)
      REAL(r8) :: DIAM0                      ! NORMALIZATION DIAMETER (microns)
      INTEGER :: ITP
      REAL(r8) :: ALPHA0,MU1,AMU,BMU,RAO1,ARAO,BRAO,RC,C1,C2,A,B
      REAL(r8) :: SUM,sum1,sum2,al,dme,dmm

      INTEGER :: J, I, J1, J2

      REAL(r8) :: V1, B1

!------------------------------------------------------------------------

      IWC00=1.
      DIAM0=1.
      dm=0._r8
      TEMPC=t-273.15_r8

   select case (ISPI)
!	if(ISPI .eq. 1) then    
!            goto 1000
!        else if(ISPI .eq. 2) then
!            goto 2000
!        end if

! 1000   CONTINUE

   case (ICE)
     !========================
     !     ICE PARTICLES
     !========================

	  dr0   = 1._r8
	  sum   = 0._r8
     bin_loop: do i=1, nr
       !  do 1300 i=1,nr

	    j1=1
	    if(i.gt.1) j1=int((r(i-1)+r(i))*.5/dr0)
	    if(i.lt.nr)j2=int((r(i+1)+r(i))*.5/dr0)
	    if(i.eq.nr)j2=int((3*r(i)-r(i-1))*.5/dr0)
	       
	    rn(i)=0.

       diameter_loop: do j=j1, j2
         !  do 1200 j=j1,j2
	      diam=2.*j*dr0     

	      rn0 = 0.          

         select case (ipsd)
         case (mh_sd)
          ! IF(IPSD .EQ. 1000) THEN     ! MH DISTRIBUTION

             iwc0 = min(CWC, 0.252*(CWC/iwc00)**0.837)
             iwc1 = CWC - iwc0

             amu = 5.2 + 0.0013*tempc
             bmu = 0.026 - 1.2e-3*tempc
             arao   = 0.47 + 2.1e-3*tempc
             brao   = 0.018 - 2.1e-4*tempc

             alpha0= -4.99e-3 - 0.0494*log10(iwc0/iwc00) 
             ! if(alpha0 .lt. 0.) alpha0=0. ! NEED TO BE POSITIVE
             ! if(iwc1 .gt. 0) mu1 = amu + bmu*log10(iwc1/iwc00)
             ! if(iwc1 .gt. 0) rao1   = arao + brao*log10(iwc1/iwc00)
             alpha0 = max(0._r8, alpha0)
             if(iwc1 .gt. 0) then
               mu1  = amu + bmu*log10(iwc1/iwc00)
               rao1 = arao + brao*log10(iwc1/iwc00)
             endif

             !... MH size dist has two parts separated at 100 micron

             !           if(diam .le. 100) 

             rn0=iwc0*alpha0**5/3.14/raoi/4    &
             &        *diam*exp(-alpha0*diam)         ! diam <=100

             if(diam .gt. 0. .and. iwc1 .gt. 0.) then
    
             rn0=rn0+6*iwc1/raoi/sqrt(3.14**3*2)  &
             &         /exp(3*mu1+9./2*rao1*rao1)/diam/rao1/diam0**3    &
             &         *exp(-0.5*((log(diam/diam0)-mu1)/rao1)**2) ! diam >100

             else
               rn0=rn0
             endif
 
         case (gamma_sd)
          ! ELSE IF(IPSD .EQ. 2000) THEN   ! GAMMA DISTRIBUTION 

             dme = 150.
             al=1.0
             rn0 = diam**(-al)*exp(-(al+3.67)*diam/dme)

         case (liu_curry_sd)
          ! ELSE IF(IPSD .EQ. 1100) THEN   ! LIU-CURRY DISTRIBUTION 

            dmm = 750. + 10*tempc
	         rn0 = exp(-5.*(diam-100.)/(dmm+500.)**0.75)

         case (knollenberg_sd)
          ! ELSE IF(IPSD .EQ. 4000) THEN   ! KNOLLENBERG DISTRIBUTION 
            v1=-1.9
            b1=0.05
            rn0=(1.e-5+diam**v1)*exp(-b1*diam)
         case default
          ! ELSE 
             ! WRITE(*,*) 'SIZE-DISTRIBUTION NOT DEFINED!'
             ! STOP
             CALL MLSMessage(MLSMSG_Error, ModuleName, &
             & ' Unrecognized ipsd parameter ')

         end select
          ! ENDIF

	      rn(i)=rn(i)+rn0*2*dr0
	      sum=sum+raoi*rn0*diam**3*3.1416/6*1e-12*2*dr0

       end do diameter_loop
     end do bin_loop
     !1200 continue
     !1300 continue

     !... converted to meters
	  do i=1,nr
	    rn(i)=rn(i)/1.e-12
	  enddo
	  sum = sum/1.e-12

     !... normalized by iwc
     do i=1,nr                             
     rn(i)=rn(i)*CWC/max(1.e-29_r8, sum)   
     enddo                                 

     !... COMPUTE MASS-MEAN-DIAMETER
     SUM1=0.
     SUM2=0.
     do i=1,nr
        SUM1=SUM1+rn(i)*1.e-12*(2.*r(I))**4
        SUM2=SUM2+rn(i)*1.e-12*(2.*r(I))**3
     enddo
     Dm=sum1/max(1.e-29_r8,sum2)

     RETURN

   case (water)
      !==============================
      !     WATER DROPLETS
      !==============================

      ! 2000   continue

        itp=1

      !... cumulus
	   if(itp .eq. 1) then
	      rc = 20.
	      c1 = 5.
	      c2 = 1.
	   endif
      !... stratus
	   if(itp .eq. 2) then
	      rc = 10.
	      c1 = 6.
	      c2 = 1.
	   endif

		b = c1/(c2*rc**c2)
		A=1.e12*(3*CWC*c2*b**((c1+4)/c2))/((4*3.14*gamma((c1+4)/c2)*1.e6))

      dr0     = 0.1  !JJ
      sum     = 0.

      bin_loop_water: do i=1, nr
       ! do 2300 i=1,nr
         

        j1=1

        if(i.gt.1) j1=int((r(i-1)+r(i))*.5/dr0)
        if(i.lt.nr)j2=int((r(i+1)+r(i))*.5/dr0)
        if(i.eq.nr)j2=int((3*r(i)-r(i-1))*.5/dr0)

        rn(i)=0.
        
        do j=j1,j2
          diam=2.*j*dr0
          rn0=A*diam**c1*exp(-b*diam**c2)
	       rn(i)=rn(i)+rn0*2*dr0
               sum=sum+rn0*diam**3*3.1416/6*1e-12*2*dr0
	 enddo

      enddo bin_loop_water
      !2300	continue

      !... normalized by LWC
      do i=1,nr
         rn(i)=rn(i)*CWC/max(1.e-29_r8, sum)
      enddo

      !... COMPUTE MASS-MEAN-DIAMETER
      SUM1=0.
      SUM2=0.
      do i=1,nr
         SUM1=SUM1+rn(i)*(2.*r(I))**4
         SUM2=SUM2+rn(i)*(2.*r(I))**3
      enddo
      Dm=sum1/max(1.e-29_r8,sum2)

    case default
      CALL MLSMessage(MLSMSG_Error, ModuleName, &
        &' Unrecognized ispi parameter--1 is ice, 2 is water ')
    end select
!      RETURN

  END SUBROUTINE DRP_SIZE

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module SizeDistribution

! $Log$
! Revision 2.4  2004/06/03 18:07:28  jonathan
! add gama function to water size distribution
!
! Revision 2.3  2004/04/15 23:11:56  jonathan
! fix a bug..Note make affect DTcir at 100mb by a factor of ~2
!
! Revision 2.2  2003/05/05 23:00:25  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.1.2.1  2003/04/16 00:57:24  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.1  2003/01/31 18:35:55  jonathan
! moved from cldfwm
!
! Revision 1.11  2003/01/30 22:01:14  pwagner
! Cosmetic changes
!
! Revision 1.10  2002/12/01 19:51:02  dwu
! fix a bug
!
! Revision 1.9  2002/10/08 17:08:08  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 1.8  2002/04/16 03:42:39  jonathan
! fix a bug
!
! Revision 1.7  2001/11/19 19:27:17  jonathan
! fix error in MH distribution
!
! Revision 1.6  2001/10/24 21:53:41  jonathan
! changes in Dm formula
!
! Revision 1.5  2001/10/23 23:57:40  jonathan
! fix mass-mean-diameter
!
! Revision 1.4  2001/10/19 18:50:21  jonathan
! fixed a bug
!
! Revision 1.3  2001/10/18 17:01:07  jonathan
! pretect divided by zero
!
! Revision 1.2  2001/09/21 15:51:37  jonathan
! modified F95 version
!

