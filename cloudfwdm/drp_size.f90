
      SUBROUTINE DRP_SIZE(ISPI,R,RN,NR,CWC,T,IPSD,Dm)

!========================================================================C
!     PRODUCE SIZE DISTRIBUTIONS FOR ICE PARTICLES AND WATER DROPLETS    C
!     J.JIANG, MAY 1, 2001
!========================================================================C

      use MLSCommon, only: r8
      IMPLICIT NONE
      INTEGER :: ISPI,ITP

      INTEGER :: NR                          ! NUMBER OF SIZE BINS

      INTEGER :: IPSD                        ! SIZE-DISTRIBUTION TYPE
                                             ! 0 = USER DEFINED
                                             ! 1 = MH
                                             ! 2 = GAMMA
                                             ! 3 = LIU-CURRY 
                                             ! 4 = KNOLLENGBERG

      REAL(r8) :: DR0                        ! INTEGRATION STEP FOR SIZE BINS
      REAL(r8) :: RN0                        ! INTEGRATION NUMBER DENSITY
      REAL(r8) :: RN(NR)                     ! INTEGRATED NUMBER DENSITY OF EACH 
                                             ! SIZE BIN
      REAL(r8) :: R(NR)                      ! PARTICLE RADIUS (micron)
      REAL(r8) :: DIAM                       ! PARTICLE DIAMETER (micron)
      REAL(r8) :: TEMPC                      ! TEMPERATURE (C)
      REAL(r8) :: T                          ! TEMPERATURE (K)
      REAL(r8) :: CWC                        ! WATER CONTENT (g/m3)
      REAL(r8) :: IWC0                       ! IWC OF SIZE < 100 microns
      REAL(r8) :: IWC1                       ! IWC OF SIZE > 100 microns
 
      REAL(r8) :: RAOI                       ! ICE DENSITY (g/m3)
      PARAMETER (RAOI=0.91_r8)

      REAL(r8) :: IWC00                      ! NORMALIZATION IWC (g/m3)
      REAL(r8) :: DIAM0                      ! NORMALIZATION DIAMETER (microns)
      REAL(r8) :: Dm                         ! MASS-MEAN-DIAMETER

      REAL(r8) :: ALPHA0,MU1,AMU,BMU,RAO1,ARAO,BRAO,RC,C1,C2,A,B
      REAL(r8) :: SUM,sum1,sum2,JJ,JJ1,JJ2,ddr0,al,dme,dmm

      INTEGER :: J, I, J1, J2

      REAL(r8) :: V1, B1

!------------------------------------------------------------------------

      IWC00=1.
      DIAM0=1.
      TEMPC=t-273.15

	if(ISPI .eq. 1) then    
            goto 1000
        else if(ISPI .eq. 2) then
            goto 2000
        end if

 1000   CONTINUE

!========================
!     ICE PARTICLES
!========================

	dr0	= 1.
	sum	= 0.
	do 1300 i=1,nr

	j1=1
	if(i.gt.1) j1=int((r(i-1)+r(i))*.5/dr0)
	if(i.lt.nr)j2=int((r(i+1)+r(i))*.5/dr0)
	if(i.eq.nr)j2=int((3*r(i)-r(i-1))*.5/dr0)
		
	rn(i)=0.

	do 1200 j=j1,j2
	diam=2.*j*dr0

	rn0 = 0.

        IF(IPSD .EQ. 1000) THEN     ! MH DISTRIBUTION

           iwc0 = min(CWC, 0.252*(CWC/iwc00)**0.837)
           iwc1 = CWC - iwc0

           amu	= 5.2 + 0.0013*tempc
           bmu	= 0.026 - 1.2e-3*tempc
           arao	= 0.47 + 2.1e-3*tempc
           brao	= 0.018 - 2.1e-4*tempc

           alpha0= -4.99e-3 - 0.0494*log10(iwc0/iwc00) 
           if(alpha0 .lt. 0.) alpha0=0. ! NEED TO BE POSITIVE
           if(iwc1 .gt. 0) mu1	= amu + bmu*log10(iwc1/iwc00)
           if(iwc1 .gt. 0) rao1	= arao + brao*log10(iwc1/iwc00)

!... MH size dist has two parts separated at 100 micron

           if(diam .le. 100) rn0=iwc0*alpha0**5/3.14/raoi/4    &
     &        *diam*exp(-alpha0*diam)
           if(diam .gt. 100 .and. iwc1 .gt. 0) &
     &        rn0=6*iwc1/raoi/sqrt(3.14**3*2)  &
     &	      /exp(3*mu1+9./2*rao1*rao1)/diam/rao1/diam0**3    &
     &	      *exp(-0.5*((log(diam/diam0)-mu1)/rao1)**2)

        ELSE IF(IPSD .EQ. 2000) THEN   ! GAMMA DISTRIBUTION 

           dme = 150.
           al=1.0
           rn0 = diam**(-al)*exp(-(al+3.67)*diam/dme)

        ELSE IF(IPSD .EQ. 1100) THEN   ! LIU-CURRY DISTRIBUTION 

          dmm = 750. + 10*tempc
	  rn0 = exp(-5.*(diam-100.)/(dmm+500.)**0.75)

        ELSE IF(IPSD .EQ. 4000) THEN   ! KNOLLENBERG DISTRIBUTION 
          v1=-1.9
          b1=0.05
          rn0=(1.e-5+diam**v1)*exp(-b1*diam)
        ELSE 
!           WRITE(*,*) 'SIZE-DISTRIBUTION NOT DEFINED!'
!           STOP
! default using MH

           iwc0 = min(CWC, 0.252*(CWC/iwc00)**0.837)
           iwc1 = CWC - iwc0

           amu	= 5.2 + 0.0013*tempc
           bmu	= 0.026 - 1.2e-3*tempc
           arao	= 0.47 + 2.1e-3*tempc
           brao	= 0.018 - 2.1e-4*tempc

           alpha0= -4.99e-3 - 0.0494*log10(iwc0/iwc00) 
           if(alpha0 .lt. 0.) alpha0=0. ! NEED TO BE POSITIVE
           if(iwc1 .gt. 0) mu1	= amu + bmu*log10(iwc1/iwc00)
           if(iwc1 .gt. 0) rao1	= arao + brao*log10(iwc1/iwc00)

!... MH size dist has two parts separated at 100 micron

           if(diam .le. 100) rn0=iwc0*alpha0**5/3.14/raoi/4  &
     &        *diam*exp(-alpha0*diam)
           if(diam .gt. 100 .and. iwc1 .gt. 0)   &
     &        rn0=6*iwc1/raoi/sqrt(3.14**3*2)    &
     &	      /exp(3*mu1+9./2*rao1*rao1)/diam/rao1/diam0**3  &
     &	      *exp(-0.5*((log(diam/diam0)-mu1)/rao1)**2)

        ENDIF

	rn(i)=rn(i)+rn0
	sum=sum+raoi*rn0*diam**3*3.1416/6*1e-12*2*dr0

1200	continue
1300	continue

!... converted to meters
	do i=1,nr
	rn(i)=rn(i)/1e-12
	enddo
	sum = sum/1e-12

!... normalized by iwc
        do i=1,nr
        rn(i)=rn(i)*CWC/sum
        enddo

!... COMPUTE MASS-MEAN-DIAMETER
        SUM1=0.
        SUM2=0.
        do i=1,nr
           SUM1=SUM1+rn(i)*r(I)**4*2
           SUM2=SUM2+rn(i)*r(I)**3
        enddo
        Dm=sum1/sum2

!	write(*,*)iwc1,iwc0,cwc,sum
!	write(*,*)(rn(i),i=1,nr)

      RETURN

!==============================
!     WATER DROPLETS
!==============================

2000	continue

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
		A=1.e12*(3*CWC*c2*b*((c1+4)/c2)/(4*3.14*3.6e6))

        dr0     = 0.1  !JJ
        sum     = 0.

        do 2300 i=1,nr


        j1=1

        if(i.gt.1) j1=int((r(i-1)+r(i))*.5/dr0)
        if(i.lt.nr)j2=int((r(i+1)+r(i))*.5/dr0)
        if(i.eq.nr)j2=int((3*r(i)-r(i-1))*.5/dr0)

        rn(i)=0.
        
        do j=j1,j2
           diam=2.*j*dr0
           rn0=a*diam**c1*exp(-b*diam**c2)
	   rn(i)=rn(i)+rn0
	enddo

	sum=sum+rn0*diam**3*3.1416/6*1e-12*2*dr0

2300	continue

!... normalized by LWC
        do i=1,nr
           rn(i)=rn(i)*CWC/sum
        enddo

!... COMPUTE MASS-MEAN-DIAMETER
        SUM1=0.
        SUM2=0.
        do i=1,nr
           SUM1=SUM1+rn(i)*r(I)**4*2
           SUM2=SUM2+rn(i)*r(I)**3
        enddo
        Dm=sum1/sum2

	return

      END

! $Log: drp_size.f90,v      

