
! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!==============================================================================
Module DailyMapModule
!==============================================================================
	
  USE MLSCommon, ONLY: r8 
  USE L3SPData, ONLY: L3SPData_T
  USE global_data
  USE SDPToolkit, ONLY: PI
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, & 
       & MLSMSG_DEALLOCATE, MLSMSG_ALLOCATE
  Implicit None
  private
  PUBLIC :: Init, ClearMemory, CordTransform, FindAD, &
    & FFSM_Opt, FFSM, FFSMA, FFSMD, Reconstruct, Diagnostics, DataGenerate, &
    & DataGeneratePrec, CopyPrec2Data
  
  PRIVATE :: ID, ModuleName

  !------------------- RCS Ident Info -----------------------
  CHARACTER(LEN=130) :: Id = &
       
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !----------------------------------------------------------

  ! Contents:

  ! Subroutines -- Init 
  !                Loc2RS
  !                CordTransform
  !                findad
  !                FFSM
  !                Reconstruct
  !                Diagnostics
  !                DataGenerate
  !                DataGeneratePrec
  !                CopyPrec2Data
  ! Function -- 
  
  ! Remarks:  This is a prototype module for the main Core processing.

  ! Parameters

  integer FFTW_FORWARD,FFTW_BACKWARD
  parameter (FFTW_FORWARD=-1,FFTW_BACKWARD=1)
  
  integer FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL
  parameter (FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1)
  
  integer FFTW_ESTIMATE,FFTW_MEASURE
  parameter (FFTW_ESTIMATE=0,FFTW_MEASURE=1)
  
  integer FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
  parameter (FFTW_OUT_OF_PLACE=0)
  parameter (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)
  
  integer FFTW_THREADSAFE
  parameter (FFTW_THREADSAFE=128)
  
  !     Constants for the MPI wrappers:
  integer FFTW_TRANSPOSED_ORDER, FFTW_NORMAL_ORDER
  integer FFTW_SCRAMBLED_INPUT, FFTW_SCRAMBLED_OUTPUT
  parameter(FFTW_TRANSPOSED_ORDER=1, FFTW_NORMAL_ORDER=0)
  parameter(FFTW_SCRAMBLED_INPUT=8192)
  parameter(FFTW_SCRAMBLED_OUTPUT=16384)
  
  Save

Contains

  ! /---------------------------------------------------------------\
  ! |*** Initilization
  ! \---------------------------------------------------------------/
  
  Subroutine Init(nt_a_i, nt_d_i, tau0_i, delTad_i, 	&
       & c0_i, lonD0_i, tD0_i, lonA0_i, tA0_i, lat_i)
    
    ! Arguments

!    CHARACTER (LEN=3) :: mode                  ! asc/des/com/all
    
    CHARACTER (LEN=480) :: msr

    REAL (r8) :: lonD0_i, tD0_i, lonA0_i, tA0_i, lat_i, delTad_i, & 
         lonD0_i_temp, lonA0_i_temp
    REAL :: c0_i, tau0_i
          
    INTEGER :: nt_a_i, nt_d_i, err
    
    lat	 = lat_i
    nt_a = nt_a_i
    nt_d = nt_d_i
    
    if(nt_a > nt_d) then
       nt = nt_d
    else
       nt = nt_a
    end if
    
    nt_a = nt
    nt_d = nt
    
    c0 	 = c0_i
    tau0 = tau0_i 
    
    sina = 1.0/sqrt(1.0+c0*c0)
    cosa = c0/sqrt(1.0+c0*c0)
    ds   = tau0/sina

    krmax = 1.0/sina
    !! jpd
    if (abs(ds*cosa).gt.0.000001) then 
       nwave = PI/(ds*cosa)
    else
       nwave = 100000
    endif
    !! jpd
    dtad = delTad_i 
    
    if(lonD0_i < 0) then
       lonD0_i_temp = 2.0*PI+lonD0_i
    else 
       lonD0_i_temp = lonD0_i
    end if
    if(lonA0_i < 0) then
       lonA0_i_temp = 2.0*PI+lonA0_i
    else 	
       lonA0_i_temp = lonA0_i
    end if
    
    dlonad    = abs( lonD0_i-lonA0_i)
    if(lonD0_i < 0 .and. lonA0_i > 0 .or. lonD0_i > 0 .and. lonA0_i < 0) then
       dlonad    = abs( lonD0_i-lonA0_i)
    else if(lonD0_i < 0 .and. lonA0_i < 0 ) then
       dlonad    = abs( lonD0_i-lonA0_i)
    else
       dlonad    = abs( lonD0_i+lonA0_i)
    end if
    
    if(lonA0_i_temp < lonD0_i_temp) then
       dlonad  = abs( lonD0_i_temp-lonA0_i_temp )
    else 
       dlonad  = 2.0*PI-abs( lonD0_i_temp-lonA0_i_temp )
    end if
    
    d1lonad   = dlonad + c0*dtad
    
    lonD0   = lonD0_i
    tD0     = tD0_i
    lonA0   = lonA0_i
    tA0     = tA0_i
    
    !	  if(mode == 'com') then
    Allocate(lonD(nt), STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Allocate // ' lonD array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    Allocate(tD(nt), STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Allocate // ' tD array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    Allocate(sD(nt), STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Allocate // ' sD array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    Allocate(rD(nt), STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Allocate // ' rD array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    Allocate(lonA(nt), STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Allocate // ' lonA array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    Allocate(tA(nt), STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Allocate // ' tA array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    Allocate(sA(nt), STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Allocate // ' sA array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    Allocate(rA(nt), STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Allocate // ' rA array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    Allocate(Dscend(nt_d), STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Allocate // ' Dscend array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    Allocate(Ascend(nt_a), STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Allocate // ' Ascend array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    Allocate(DPrec(nt_d), STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Allocate // ' DPrec array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    Allocate(APrec(nt_a), STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Allocate // ' APrec array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    Allocate(wn(nt*2), STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Allocate // ' wn array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    Allocate(sigma(nt*2), STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Allocate // ' sigma array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    Allocate(wna(nt_a*2), STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Allocate // ' wna array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    Allocate(sigmaa(nt_a*2), STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Allocate // ' sigmaa array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    Allocate(wnd(nt_d*2), STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Allocate // ' wnd array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    Allocate(sigmad(nt_d*2), STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Allocate // ' sigmad array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    Allocate(phikr(nt*2), STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Allocate // ' phikr array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    Allocate(phikra(nt_a*2), STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Allocate // ' phikra array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    Allocate(phikrd(nt_d*2), STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Allocate // ' phikrd array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    !          else if(mode == 'asc') then
    !	    Allocate(lonDA(nt), tDA(nt), sDA(nt), rDA(nt))
    !	    Allocate(DAcend(nt))
    !	    Allocate(wn(nt*2), sigma(nt*2))
    !	    Allocate(phikrda(nt*2))
    !	  end if
  End Subroutine Init

  ! /---------------------------------------------------------------\
  ! |*** Clear Memory
  ! \---------------------------------------------------------------/

  Subroutine ClearMemory()

    CHARACTER (LEN=480) :: msr
    
    INTEGER :: err

    DeAllocate(lonD, STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Deallocate // ' lonD array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    DeAllocate(tD, STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Deallocate // ' tD array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    DeAllocate(sD, STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Deallocate // ' sD array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    DeAllocate(rD, STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Deallocate // ' rD array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    DeAllocate(lonA, STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Deallocate // ' lonA array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    DeAllocate(tA, STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Deallocate // ' tA array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    DeAllocate(sA, STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Deallocate // ' sA array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    DeAllocate(rA, STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Deallocate // ' rA array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    DeAllocate(Dscend, STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Deallocate // ' Dscend array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    DeAllocate(Ascend, STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Deallocate // ' Ascend array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    DeAllocate(DPrec, STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Deallocate // ' DPrec array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    DeAllocate(APrec, STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Deallocate // ' APrec array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    DeAllocate(wn, STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Deallocate // ' wn array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    DeAllocate(sigma, STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Deallocate // ' sigma array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    DeAllocate(wna, STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Deallocate // ' wna array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    DeAllocate(sigmaa, STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Deallocate // ' sigmaa array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    DeAllocate(wnd, STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Deallocate // ' wnd array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    DeAllocate(sigmad, STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Deallocate // ' sigmad array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    DeAllocate(phikr, STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Deallocate // ' phikr array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    DeAllocate(phikra, STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Deallocate // ' phikra array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    DeAllocate(phikrd, STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Deallocate // ' phikrd array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF
    
  End Subroutine ClearMemory
  
  ! /---------------------------------------------------------------\
  ! |*** Coordinate transformation
  ! \---------------------------------------------------------------/
  
  Subroutine CordTransform(mode)
    
    ! Arguments
    
    CHARACTER (LEN=3) :: mode                  ! asc/des/com/all
    
    ! Variables
    
    Integer :: i
    
    if(mode == 'com') then
       sD0 = lonD0*cosa - tD0*sina - tau0*(nt-1.0)/sina
       sA0 = sD0 + (-dlonad*c0 + dtad)*sina
       rD0 = (lonD0 + c0*tD0)*sina
       rA0 = (lonD0 - dlonad + c0*tD0 - c0*dtad)*sina
    else if(mode == 'asc') then
       sA0 = lonA0*cosa - tA0*sina - tau0*(nt-1.0)/sina
       rA0 = (lonA0 + c0*tA0)*sina
    else if(mode == 'des') then
       sD0 = lonD0*cosa - tD0*sina - tau0*(nt-1.0)/sina
       rD0 = (lonD0 + c0*tD0)*sina
    end if
    
    if(mode == 'com') then
       Do i = 1, nt
          rD(i) = rD0
          rA(i) = rA0
          sD(i) = sD0 + i*tau0/sina
          sA(i) = sA0 + i*tau0/sina
       End Do
    else if(mode == 'asc') then
       Do i = 1, nt
          rA(i) = rA0
          sA(i) = sA0 + i*tau0/sina
       End Do
    else if(mode == 'des') then
       Do i = 1, nt
          rD(i) = rD0
          sD(i) = sD0 + i*tau0/sina
       End Do
    end if
    
  End Subroutine CordTransform
  
  
  ! /---------------------------------------------------------------\
  ! |*** findad 
  ! \---------------------------------------------------------------/
  
  Subroutine FindAD(sina, cosa, t, r0, ep, em, frp, frm, rpi, rmi)
    
    ! Arguments
    
    Real :: sina, cosa, t, r0, ep, em, frp, frm, rpi, rmi
    
    ! Variables
    
    Real :: e, r, s
    Integer :: inc
    
    e = ep
    
    r = e*sina + t*cosa
    s = e*cosa - t*sina
    
    if (abs(r-r0) >= 1.e-6) then 
       
       if ( r > r0 ) then 
          inc = -1
          ep = e
       else 
          inc = 1
          em = e
       end if
       
       e = e + inc*2.0*PI
       r = e*sina + t*cosa
       
       if ( inc < 0) then 
          do
             if (r <= r0) exit
             ep = e
             e  = e + inc*2.0*PI
             r = e*sina + t*cosa
          end do
          em = e
       endif
       
       if ( inc > 0) then 
          do
             if (r >= r0) exit
             em = e
             e  = e + inc*2.0*PI
             r = e*sina + t*cosa
          end do
          ep = e
       endif
       
       rpi = ep*sina + t*cosa
       rmi = em*sina + t*cosa
       
       if( abs(r0-rpi) <= 1.e-8 .and. abs(r0-rmi) <= 1.e-8) then 
          frm = 1.0
          frp = 0.0
       else 
          frm = (r0-rpi)**2/((r0-rpi)**2+(r0-rmi)**2)
          frp = (r0-rmi)**2/((r0-rpi)**2+(r0-rmi)**2)
       end if
       
    else 
       frm = 1.0
       frp = 0.0
       rpi = ep*sina + t*cosa
       rmi = em*sina + t*cosa
    end if
    
  End Subroutine findad
  
  ! /---------------------------------------------------------------\
  ! |*** FFSM_Opt --- FFSM Mode Passing Routine
  ! \---------------------------------------------------------------/
  
  Subroutine FFSM_Opt(mode, l3sp, iLv, iLt)
    
    ! Variables
    
    CHARACTER (LEN=3) :: mode                  ! asc/des/com/all
    TYPE( L3SPData_T ), POINTER :: l3sp(:)
    
    integer iLv, iLt
    
    
    if(mode == 'com') then
       CALL FFSM(l3sp(1), iLv, iLt)
    else if(mode == 'asc') then
       CALL FFSMA(l3sp(1), iLv, iLt)
    else if(mode == 'des') then
       CALL FFSMD(l3sp(1), iLv, iLt)
    else if(mode == 'all') then
       CALL FFSM(l3sp(1), iLv, iLt)
       CALL FFSMA(l3sp(2), iLv, iLt)
       CALL FFSMD(l3sp(3), iLv, iLt)
    end if
    
  End Subroutine FFSM_Opt
  
  ! /---------------------------------------------------------------\
  ! |*** FFSM --- combined mode
  ! \---------------------------------------------------------------/
  
  Subroutine FFSM(l3sp, iLv, iLt)
    
    ! Variables
    
    TYPE( L3SPData_T ) :: l3sp
    !TYPE( L3SPData_T ), POINTER :: l3sp
    
    double complex phia, phid
    
    REAL(r8), Dimension((nt+1)/2) :: imgPtemp
    REAL(r8), Dimension(nt) :: & 
         & dscend_fft, ascend_fft, drealP, dimgP, arealP, aimgP 
    
    real :: ks, kr, expa, expd, krmax_g
    
    INTEGER :: plan, m, m1, i, j, flag, iLv, iLt, wIndex
    
    ! Fourier Transform of both descending & ascending series
    
    call rfftw_f77_create_plan(plan,nt,FFTW_FORWARD,FFTW_ESTIMATE)
    call rfftw_f77_one(plan,ascend,ascend_fft)
    call rfftw_f77_one(plan,dscend,dscend_fft)
    call rfftw_f77_destroy_plan(plan)
    
    ! Ascending Part
    
    Do i = nt/2, nt
       arealP(i) = ascend_fft(i+1-nt/2)/float(nt)
    End Do
    Do i = 1, nt/2-1
       arealP(i) = ascend_fft(nt/2-i+1)/float(nt)
    End Do
    
    Do i = 1, (nt+1)/2-1
       imgPTemp(i) = ascend_fft(nt+1-i)/float(nt)
    End Do
    
    IF (mod(nt, 2) == 0) THEN
       imgPTemp((nt+1)/2) = 0.0
    ENDIF
    
    aimgP(nt/2) = 0.0 
    Do i = nt/2+1, nt-1
       aimgP(i) = imgPtemp(i-nt/2)
    End Do
    Do i = 1, nt/2-1
       aimgP(i) = -imgPtemp(nt/2-i)
    End Do
    aimgP(nt) = 0.0
    
    ! Dscending Part
    
    Do i = 1, (nt+1)/2
       imgPTemp(i) = 0.0 
    End Do
    
    Do i = nt/2, nt
       drealP(i) = dscend_fft(i+1-nt/2)/float(nt)
    End Do
    Do i = 1, nt/2-1
       drealP(i) = dscend_fft(nt/2-i+1)/float(nt)
    End Do
    
    Do i = 1, (nt+1)/2-1
       imgPTemp(i) = dscend_fft(nt+1-i)/float(nt)
    End Do
    
    IF (mod(nt, 2) == 0) THEN
       imgPTemp((nt+1)/2) = 0.0
    ENDIF
    
    dimgP(nt/2) = 0.0 
    Do i = nt/2+1, nt-1
       dimgP(i) = imgPtemp(i-nt/2)
    End Do
    Do i = 1, nt/2-1
       dimgP(i) = -imgPtemp(nt/2-i)
    End Do
    dimgP(nt) = 0.0
    
    mtotal = 0
    
    krmax_g = krmax

    Do m = 1, nwave 
       m1 = m-1
       Do i = 1, nt
          
          mtotal = mtotal + 1
          
          wn(mtotal) = float(m1)
          
          if (abs(float(nt)*ds).gt.0.00001) ks = 2.*PI*(-nt/2+i)/(float(nt)*ds)
          !ks = 2.*PI*(-nt/2+i)/(float(nt)*ds)
          if (abs(sina) .gt.0.00001) kr = -ks*c0+float(m1)/sina
          !kr = -ks*c0+float(m1)/sina
          
          sigma(mtotal) = -ks*sina + kr*cosa
          
          if (kr < -krmax_g .or. kr > 0.0 .or. ks < 0.0) then 
             flag = 0
             phikr(mtotal) = CMPLX(0.0, 0.0)
             Goto 101
          endif
          
          expa = ks*sa0 + kr*ra0
          expd = ks*sd0 + kr*rd0
          
          flag = 0
          
          phia   = CMPLX(arealP(i), aimgP(i))*CMPLX(cos(expa), -sin(expa))
          phid   = CMPLX(drealP(i), dimgP(i))*CMPLX(cos(expd), -sin(expd))
          
          phikr(mtotal) = & 
               & (phia-phid*CMPLX(cos(d1lonad), -sin(d1lonad)))/& 
               & (1.0-CMPLX(cos(d1lonad), -sin(d1lonad)))

          if( phikr(mtotal) /= CMPLX(0.0, 0.0) ) flag = 1
          
          
101       continue 
          
          ! *** kr+
          
          if (kr > krmax_g .or. kr < 0.0 .or. ks < 0.0 ) then 
             if (flag == 0 .and. mtotal > 0) mtotal = mtotal-1
             Goto 102
          endif
          
          kr  = -ks*c0+float(m1-1)/sina
          
          expa = ks*sa0 + kr*ra0
          expd = ks*sd0 + kr*rd0
          
          phia   = CMPLX(arealP(i), aimgP(i))*CMPLX(cos(expa), -sin(expa))
          phid   = CMPLX(drealP(i), dimgP(i))*CMPLX(cos(expd), -sin(expd))
          
          phikr(mtotal) = & 
               (phid-phia)*CMPLX( cos(rd0/sina), -sin(rd0/sina) )/& 
               (1.0-CMPLX( cos(d1lonad), -sin(d1lonad) ))
          
102       continue 
          
       End Do
    End Do

    
    do j = 1, mtotal
       l3sp%waveNumber(iLv, iLt, j) =  int(wn(j))
       l3sp%frequency(iLv, iLt, j) = sigma(j)
       l3sp%l3spRelValue(iLv, iLt, j) = real(phikr(j))
       l3sp%l3spRelPrecision(iLv, iLt, j) = real(phikr(j))
       l3sp%l3spImgValue(iLv, iLt, j) = aimag(phikr(j))
       l3sp%l3spImgPrecision(iLv, iLt, j) = aimag(phikr(j))
    end do
    
    
  End Subroutine FFSM

  ! /---------------------------------------------------------------\
  ! |*** FFSMA --- single mode / ascending
  ! \---------------------------------------------------------------/
  
  Subroutine FFSMA(l3sp, iLv, iLt)
    
    Implicit None

    ! Variables
    
    !TYPE( L3SPData_T ), POINTER :: l3sp
    TYPE( L3SPData_T ) :: l3sp
    
    double complex phida

    REAL(r8), Dimension((nt_a+1)/2) :: imgPtemp
    REAL(r8), Dimension(nt_a) :: ascend_fft, arealP, aimgP
    
    REAL :: ks, kr, expda, krmax_g
    
    INTEGER :: plan, m, m1, i, j, flag, iLv, iLt, wIndex
    
    ! Fourier Transform of one series
    
    call rfftw_f77_create_plan(plan,nt_a,FFTW_FORWARD,FFTW_ESTIMATE)
    call rfftw_f77_one(plan,ascend,ascend_fft)
    call rfftw_f77_destroy_plan(plan)

    ! re-order the fft transform result

    Do i = 1, (nt_a+1)/2
       imgPTemp(i) = 0.0 
    End Do
    
    Do i = nt_a/2, nt_a
       arealP(i) = ascend_fft(i+1-nt_a/2)/float(nt_a)
    End Do
    Do i = 1, nt_a/2-1
       arealP(i) = ascend_fft(nt_a/2-i+1)/float(nt_a)
    End Do
    
    Do i = 1, (nt_a+1)/2-1
       imgPTemp(i) = ascend_fft(nt_a+1-i)/float(nt_a)
    End Do
    
    IF (mod(nt_a, 2) == 0) THEN
       imgPTemp((nt_a+1)/2) = 0.0
    ENDIF
    
    aimgP(nt_a/2) = 0.0 
    Do i = nt_a/2+1, nt_a
       aimgP(i) = imgPtemp(i-nt_a/2)
    End Do
    Do i = 1, nt_a/2-1
       aimgP(i) = -imgPtemp(nt_a/2-i)
    End Do
    
    mtotala = 0
    
    krmax_g = krmax/2.0
    
    Do m = 1, nwave 
       m1 = m-1
       Do i = 1, nt_a
          mtotala = mtotala + 1
          
          wna(mtotala) = float(m1)
          
          !! floating point error

          if (abs(float(nt_a)*ds) .gt. 0.0000001) ks = 2.0*PI*(-nt_a/2+i)/(float(nt_a)*ds)
               !ks = 2.0*PI*(-nt_a/2+i)/(float(nt_a)*ds)
          if (abs(sina) .gt. 0.0000001) kr = -ks*c0+float(m1)/sina
          !kr = -ks*c0+float(m1)/sina
          
          !! floating point error

          sigmaa(mtotala) = -ks*sina + kr*cosa
          
          if (kr < -krmax_g .or. kr > 0.0 .or. ks < 0.0) then 
             flag = 0
             Goto 101
          endif
          
          expda = ks*sa0 + kr*ra0
          
          flag = 0
          
          phida   = CMPLX(arealP(i), aimgP(i))*CMPLX(cos(expda), -sin(expda))
          
          phikra(mtotala) = phida 
          
 
          if( phikra(mtotala) /= CMPLX(0.0, 0.0) ) flag = 1
          
          
101       continue 
          
          ! *** kr+
          
          if (kr > krmax_g .or. kr < 0.0 .or. ks < 0.0 ) then 
             if (flag == 0 .and. mtotala > 0) mtotala = mtotala-1
             Goto 102
          endif
          
          kr  = -ks*c0+float(m1)/sina
          
          expda = ks*sa0 + kr*ra0
          
          phida   = CMPLX(arealP(i), aimgP(i))*CMPLX(cos(expda), -sin(expda))

          phikra(mtotala) = phida 
          
102       continue 
          
       End Do
    End Do
    
    
    do j = 1, mtotala
       l3sp%waveNumber(iLv, iLt, j) = int(wna(j)) 
       l3sp%frequency(iLv, iLt, j) = sigmaa(j)
       l3sp%l3spRelValue(iLv, iLt, j) = real(phikra(j))
       l3sp%l3spRelPrecision(iLv, iLt, j) = real(phikra(j))
       l3sp%l3spImgValue(iLv, iLt, j) = aimag(phikra(j))
       l3sp%l3spImgPrecision(iLv, iLt, j) = aimag(phikra(j))
    end do
    
    
  End Subroutine FFSMA
  
  
  ! /---------------------------------------------------------------\
  ! |*** FFSMD --- single mode / descending
  ! \---------------------------------------------------------------/
  
  Subroutine FFSMD(l3sp, iLv, iLt)

    Implicit None
    
    ! Variables
    
    !TYPE( L3SPData_T ), POINTER :: l3sp
    TYPE( L3SPData_T ) :: l3sp
    
    double complex phida
    
    REAL(r8), Dimension((nt_d+1)/2) :: imgPtemp
    REAL(r8), Dimension(nt_d) :: dscend_fft, drealP, dimgP
    
    REAL :: ks, kr, expda, krmax_g
    
    INTEGER ::  plan, m, m1, i, j, flag, iLv, iLt, wIndex
    
    ! Fourier Transform of one series

    call rfftw_f77_create_plan(plan,nt_d,FFTW_FORWARD,FFTW_ESTIMATE)
    call rfftw_f77_one(plan,dscend,dscend_fft)
    call rfftw_f77_destroy_plan(plan)

    ! re-order the fft transform result

    Do i = 1, (nt_d+1)/2
       imgPTemp(i) = 0.0 
    End Do
    
    Do i = nt_d/2, nt_d
       drealP(i) = dscend_fft(i+1-nt_d/2)/float(nt_d)
    End Do
    Do i = 1, nt_d/2-1
       drealP(i) = dscend_fft(nt_d/2-i+1)/float(nt_d)
    End Do
    
    Do i = 1, (nt_d+1)/2-1
       imgPTemp(i) = dscend_fft(nt_d+1-i)/float(nt_d)
    End Do
    
    IF (mod(nt_d, 2) == 0) THEN
       imgPTemp((nt_d+1)/2) = 0.0
    ENDIF
    
    dimgP(nt_d/2) = 0.0 
    
    Do i = nt_d/2+1, nt_d
       dimgP(i) = imgPtemp(i-nt_d/2)
    End Do
    Do i = 1, nt_d/2-1
       dimgP(i) = -imgPtemp(nt_d/2-i)
    End Do
    
    mtotald = 0
    
    krmax_g = krmax/2.0
    
    Do m = 1, nwave 
       m1 = m-1
       Do i = 1, nt_d
          
          mtotald = mtotald + 1
          
          wnd(mtotald) = float(m1)
          !! floating point 
          if (abs(float(nt_d)*ds).gt.0.000001) ks = 2.0*PI*(-nt_d/2+i)/(float(nt_d)*ds)
               !ks = 2.0*PI*(-nt_d/2+i)/(float(nt_d)*ds)
          if (abs(sina).gt.0.000001) kr = -ks*c0+float(m1)/sina
          !kr = -ks*c0+float(m1)/sina
          !! floating point
          sigmad(mtotald) = -ks*sina + kr*cosa

          if (kr < -krmax_g .or. kr > 0.0 .or. ks < 0.0) then 
             flag = 0
             Goto 101
          endif
          
          expda = ks*sd0 + kr*rd0
          
          flag = 0

          phida   = CMPLX(drealP(i), dimgP(i))*CMPLX(cos(expda), -sin(expda))
          
          phikrd(mtotald) = phida 
          
          
          if( phikrd(mtotald) /= CMPLX(0.0, 0.0) ) flag = 1
          
          
101       continue 
          
          ! *** kr+

          if (kr > krmax_g .or. kr < 0.0 .or. ks < 0.0 ) then 
             if (flag == 0 .and. mtotald > 0) mtotald = mtotald-1
             Goto 102
          endif
          
          kr  = -ks*c0+float(m1)/sina
          
          expda = ks*sd0 + kr*rd0
          
          phida   = CMPLX(drealP(i), dimgP(i))*CMPLX(cos(expda), -sin(expda))
          
          phikrd(mtotald) = phida 
          
102       continue 
          
       End Do
    End Do
    
    
    do j = 1, mtotald
       l3sp%waveNumber(iLv, iLt, j) =  int(wnd(j))
       l3sp%frequency(iLv, iLt, j) = sigmad(j)
       l3sp%l3spRelValue(iLv, iLt, j) = real(phikrd(j))
       l3sp%l3spRelPrecision(iLv, iLt, j) = real(phikrd(j))
       l3sp%l3spImgValue(iLv, iLt, j) = aimag(phikrd(j))
       l3sp%l3spImgPrecision(iLv, iLt, j) = aimag(phikrd(j))
    end do
    
    
  End Subroutine FFSMD

  
  ! /---------------------------------------------------------------\
  ! |*** Reconstruct 
  ! \---------------------------------------------------------------/
  
  Subroutine Reconstruct(mode, xtime, nlons, xlon_orig, result)
    
    ! Arguments
    
    CHARACTER (LEN=3) :: mode                  ! asc/des/com/all
    
    INTEGER :: nlons
    
    REAL :: xtime
    
    REAL(r8), Dimension(nlons) :: xlon, xlon_orig, result

    ! Variables

    REAL :: kr, ks, xloni, eoff, rp, rm, ep, em, sum, & 
         & argpa, argpd, argma, argmd, argp, argm, & 
         & epa, epd, ema, emd, rma, rmd, rpa, rpd , &
         & rdp, ran, rdn, frpa, frma, frpd, frmd, frp, frm, rpi, rmi

    INTEGER :: j, k

    !*** Convert to gradient 
    
    do k = 1, nlons
       xlon(k) = xlon_orig(k)*PI/180.0
    end do
    
    !*** find offset in longitude necessary to get r value right
    
    eoff = c0*xtime - mod( (c0*xtime), real(2.0*PI) )
    if (eoff < 0.0) eoff = 2.0*PI
    
    !*** Loop over longitude
    
    do k = 1, nlons
       
       result(k) = 0.0
       
       xloni = xlon(k) - eoff
       
       if (mode == 'com') then 
          
          !***find r value corresponding to this xlon & xtime
          
          epa = xloni
          ema = xloni
          epd = xloni
          emd = xloni
          
          rdp = xloni*sina + xtime*cosa

          call findad(sina, cosa, xtime, ra0, epa, ema, frpa, frma, rpi, rmi)
          call findad(sina, cosa, xtime, rd0, epd, emd, frpd, frmd, rpi, rmi)

          rpa = epa*sina + xtime*cosa
          rma = ema*sina + xtime*cosa
          rpd = epd*sina + xtime*cosa
          rmd = emd*sina + xtime*cosa
          
          if( (rpa-ra0)**2 >= (rma-ra0)**2 ) ran = (rma-ra0)**2
          if( (rpa-ra0)**2 <  (rma-ra0)**2 ) ran = (rpa-ra0)**2
          
          if( (rpd-rd0)**2 >= (rmd-rd0)**2 ) rdn = (rmd-rd0)**2
          if( (rpd-rd0)**2 <  (rmd-rd0)**2 ) rdn = (rpd-rd0)**2

          frpa = frpa*rdn/(rdn+ran)
          frma = frma*rdn/(rdn+ran)
          frpd = frpd*ran/(rdn+ran)
          frmd = frmd*ran/(rdn+ran)


          do j = 1, mtotal 

             kr     = wn(j)*sina+sigma(j)*cosa
             ks     = wn(j)*cosa-sigma(j)*sina

             argpa = (wn(j)*epa + sigma(j)*xtime) + kr*(ra0-rpa)/(2.0*PI)
             argma = (wn(j)*ema + sigma(j)*xtime) + kr*(ra0-rma)/(2.0*PI)
             argpd = (wn(j)*epd + sigma(j)*xtime) + kr*(rd0-rpd)/(2.0*PI)
             argmd = (wn(j)*emd + sigma(j)*xtime) + kr*(rd0-rmd)/(2.0*PI)
             !argpa = (wn(j)*epa + sigma(j)*xtime) + kr*(ra0-rpa)
             !argma = (wn(j)*ema + sigma(j)*xtime) + kr*(ra0-rma)
             !argpd = (wn(j)*epd + sigma(j)*xtime) + kr*(rd0-rpd)
             !argmd = (wn(j)*emd + sigma(j)*xtime) + kr*(rd0-rmd)
             
             sum =  real( frpa*phikr(j)*CMPLX(cos(argpa), sin(argpa)) ) + &
                  & real( frma*phikr(j)*CMPLX(cos(argma), sin(argma)) ) + &
                  & real( frpd*phikr(j)*CMPLX(cos(argpd), sin(argpd)) ) + &
                  & real( frmd*phikr(j)*CMPLX(cos(argmd), sin(argmd)) ) 
             
             if( abs(ks) > 1.e-4) then 
                sum = sum + & 
                     real( conjg(frpa*phikr(j))*CMPLX(cos(argpa), & 
                     -sin(argpa)) ) + &
 
                     real( conjg(frma*phikr(j))*CMPLX(cos(argma), & 
                     -sin(argma)) ) + &

                     real( conjg(frpd*phikr(j))*CMPLX(cos(argpd), & 
                     -sin(argpd)) ) + &

                     real( conjg(frmd*phikr(j))*CMPLX(cos(argmd), & 
                     -sin(argmd)) ) 
             endif
             
             result(k) = result(k) + sum

          end do

          
       else if (mode == 'des') then 

          ep = xloni
          em = xloni
          
          call findad( sina, cosa, xtime, rd0, ep, em, frp, frm, rpi, rmi)
          
          rp = rpi 
          rm = rmi 

          do j = 1, mtotald 

             kr     = wnd(j)*sina+sigmad(j)*cosa
             ks     = wnd(j)*cosa-sigmad(j)*sina

             argp = (wnd(j)*ep + sigmad(j)*xtime) + kr*(rd0-rp)/(2.0*PI)
             argm = (wnd(j)*em + sigmad(j)*xtime) + kr*(rd0-rm)/(2.0*PI)

             sum = real( frp*phikrd(j)*CMPLX(cos(argp), sin(argp)) ) + &
                  real( frm*phikrd(j)*CMPLX(cos(argm), sin(argm)) ) 
             
             if( abs(ks) > 1.e-4) then 
                sum = sum + & 
                     & real( conjg(frp*phikrd(j))*CMPLX(cos(argp),-sin(argp)))&
                     & + &
                     & real( conjg(frm*phikrd(j))*CMPLX(cos(argm), -sin(argm)))
             endif
             
             result(k) = result(k) + sum
             
          end do


       else if (mode == 'asc') then 
          
          ep = xloni
          em = xloni

          call findad( sina, cosa, xtime, ra0, ep, em, frp, frm, rpi, rmi)

          rp = rpi 
          rm = rmi 
                
          do j = 1, mtotala 

             kr     = wna(j)*sina+sigmaa(j)*cosa
             ks     = wna(j)*cosa-sigmaa(j)*sina

             argp = (wna(j)*ep + sigmaa(j)*xtime) + kr*(ra0-rp)/(2.0*PI)
             argm = (wna(j)*em + sigmaa(j)*xtime) + kr*(ra0-rm)/(2.0*PI)

             sum = real( frp*phikra(j)*CMPLX(cos(argp), sin(argp)) ) + &
                  & real( frm*phikra(j)*CMPLX(cos(argm), sin(argm)) ) 

             if( abs(ks) > 1.e-4) then 
                sum = sum + & 
                     & real(conjg(frp*phikra(j))*CMPLX(cos(argp), -sin(argp)))&
                     & + &
                     & real(conjg(frm*phikra(j))*CMPLX(cos(argm), -sin(argm))) 
             endif
             
             result(k) = result(k) + sum

          end do

       end if

    end do

  End Subroutine Reconstruct

  
  ! /---------------------------------------------------------------\
  ! |*** Diagnostics 
  ! \---------------------------------------------------------------/

  Subroutine Diagnostics(mode, xtime_in, xlon_in, result)

    ! Arguments

    CHARACTER (LEN=3) :: mode                  ! asc/des/com/all
    
    REAL(r8) :: xtime_in, xlon_in
    REAL :: xtime, xlon, result
    
    ! Variables
    
    REAL :: kr, ks, xloni, eoff, rp, rm, ep, em, sum, &
         & argpa, argpd, argma, argmd, argp, argm, &
         & epa, epd, ema, emd, rma, rmd, rpa, rpd , & 
         & rdp, ran, rdn, frpa, frma, frpd, frmd, frp, frm, rpi, rmi

    INTEGER :: j

    !*** find offset in longitude necessary to get r value right

    xtime = xtime_in
    xlon  = xlon_in

    eoff = c0*xtime - mod( (c0*xtime), real(2.0*PI) )
    if (eoff < 0.0) eoff = 2.0*PI
    
    result = 0.0
    
    xloni = xlon - eoff
    
    if (mode == 'com') then 

       !***find r value corresponding to this xlon & xtime

       epa = xloni
       ema = xloni
       epd = xloni
       emd = xloni

       rdp = xloni*sina + xtime*cosa
       
       call findad(sina, cosa, xtime, ra0, epa, ema, frpa, frma, rpi, rmi)
       call findad(sina, cosa, xtime, rd0, epd, emd, frpd, frmd, rpi, rmi)

       rpa = epa*sina + xtime*cosa
       rma = ema*sina + xtime*cosa
       rpd = epd*sina + xtime*cosa
       rmd = emd*sina + xtime*cosa
       
       if( (rpa-ra0)**2 >= (rma-ra0)**2 ) ran = (rma-ra0)**2
       if( (rpa-ra0)**2 <  (rma-ra0)**2 ) ran = (rpa-ra0)**2
       
       if( (rpd-rd0)**2 >= (rmd-rd0)**2 ) rdn = (rmd-rd0)**2
       if( (rpd-rd0)**2 <  (rmd-rd0)**2 ) rdn = (rpd-rd0)**2

       frpa = frpa*rdn/(rdn+ran)
       frma = frma*rdn/(rdn+ran)
       frpd = frpd*ran/(rdn+ran)
       frmd = frmd*ran/(rdn+ran)

       do j = 1, mtotal 

          kr     = wn(j)*sina+sigma(j)*cosa
          ks     = wn(j)*cosa-sigma(j)*sina

          argpa = (wn(j)*epa + sigma(j)*xtime) + kr*(ra0-rpa)/(2.0*PI)
          argma = (wn(j)*ema + sigma(j)*xtime) + kr*(ra0-rma)/(2.0*PI)
          argpd = (wn(j)*epd + sigma(j)*xtime) + kr*(rd0-rpd)/(2.0*PI)
          argmd = (wn(j)*emd + sigma(j)*xtime) + kr*(rd0-rmd)/(2.0*PI)
          
          sum = real( frpa*phikr(j)*CMPLX(cos(argpa), sin(argpa))  ) + &
               & real( frma*phikr(j)*CMPLX(cos(argma), sin(argma)) ) + &
               & real( frpd*phikr(j)*CMPLX(cos(argpd), sin(argpd)) ) + &
               & real( frmd*phikr(j)*CMPLX(cos(argmd), sin(argmd)) ) 
          
          if( abs(ks) > 1.e-4) then 
             sum = sum + & 
                  & real(conjg(frpa*phikr(j))*CMPLX(cos(argpa), -sin(argpa))) &
                  & + &
                  & real(conjg(frma*phikr(j))*CMPLX(cos(argma), -sin(argma))) &
                  & + &
                  & real(conjg(frpd*phikr(j))*CMPLX(cos(argpd), -sin(argpd))) &
                  & + &
                  & real(conjg(frmd*phikr(j))*CMPLX(cos(argmd), -sin(argmd))) 
          endif

          result = result + sum
          
       end do

    else if (mode == 'des') then 

       ep = xloni
       em = xloni

       call findad( sina, cosa, xtime, rd0, ep, em, frp, frm, rpi, rmi)

       rp = rpi 
       rm = rmi 

       do j = 1, mtotald 
          
          kr     = wnd(j)*sina+sigmad(j)*cosa
          ks     = wnd(j)*cosa-sigmad(j)*sina
          
          argp = (wnd(j)*ep + sigmad(j)*xtime) + kr*(rd0-rp)/(2.0*PI)
          argm = (wnd(j)*em + sigmad(j)*xtime) + kr*(rd0-rm)/(2.0*PI)
          
          sum = real( frp*phikrd(j)*CMPLX(cos(argp), sin(argp)) ) + &
               & real( frm*phikrd(j)*CMPLX(cos(argm), sin(argm)) ) 
          
          if( abs(ks) > 1.e-4) then 
             sum = sum + & 
                  & real( conjg(frp*phikrd(j))*CMPLX(cos(argp), -sin(argp)))+ &
                  & real( conjg(frm*phikrd(j))*CMPLX(cos(argm), -sin(argm))) 
          endif

          result = result + sum
          
       end do
       
    else if (mode == 'asc') then 

       ep = xloni
       em = xloni

       call findad( sina, cosa, xtime, ra0, ep, em, frp, frm, rpi, rmi)

       rp = rpi 
       rm = rmi 
                
       do j = 1, mtotala 

          kr     = wna(j)*sina+sigmaa(j)*cosa
          ks     = wna(j)*cosa-sigmaa(j)*sina

          argp = (wna(j)*ep + sigmaa(j)*xtime) + kr*(ra0-rp)/(2.0*PI)
          argm = (wna(j)*em + sigmaa(j)*xtime) + kr*(ra0-rm)/(2.0*PI)

          sum = real( frp*phikra(j)*CMPLX(cos(argp), sin(argp)) ) + &
               & real(frm*phikra(j)*CMPLX(cos(argm), sin(argm))) 

          if( abs(ks) > 1.e-4) then 
             sum = sum + & 
                 & real(conjg(frp*phikra(j))*CMPLX(cos(argp), -sin(argp))) + &
                 & real(conjg(frm*phikra(j))*CMPLX(cos(argm), -sin(argm))) 
          endif

          result = result + sum

       end do

    end if
    
    
  End Subroutine Diagnostics

  
  ! /---------------------------------------------------------------\
  ! |*** Generate Data
  ! \---------------------------------------------------------------/

  Subroutine DataGenerate(aField, dField)

    ! Arguments

    REAL (r8), DIMENSION(:) ::  aField, dField

    ! Variables

	 INTEGER :: i
         
	 Do i = 1, nt
            Dscend(nt+1-i) = dField(i) 
            Ascend(nt+1-i) = aField(i) 
         End Do

  End Subroutine DataGenerate


  Subroutine DataGeneratePrec(aField, dField)

    ! Arguments

    REAL (r8), DIMENSION(:) ::  aField, dField

    ! Variables

	 INTEGER :: i
         
	 Do i = 1, nt
            DPrec(nt+1-i) = dField(i) 
            APrec(nt+1-i) = aField(i) 
         End Do

  End Subroutine DataGeneratePrec


  Subroutine CopyPrec2Data()
         
	Dscend = DPrec
	Ascend = APrec

  End Subroutine CopyPrec2Data
       
!===================
End Module DailyMapModule
!===================

! $Log$
! Revision 1.11  2003/04/30 21:44:32  pwagner
! Work-around for LF95 infinite compile-time bug
!
! Revision 1.10  2003/03/22 03:42:55  jdone
! zero denominator check, individual allocate/deallocate, use only and indentation added
!
! Revision 1.9  2002/02/20 21:18:13  ybj
! *** empty log message ***
!
! Revision 1.8  2001/08/17 20:27:51  nakamura
! Fixed bug & core dump problems.
!
! Revision 1.7  2001/08/14 16:13:24  nakamura
! Adjusted initialization of fNum in FFS routines.
!
! Revision 1.6  2001/08/13 16:42:06  ybj
! *** empty log message ***
!
! Revision 1.5  2001/04/13 22:07:05  ybj
! reasonable values
!
! Revision 1.4  2001/04/11 18:29:30  ybj
! reasonable values
!
! Revision 1.3  2001/03/06 18:52:41  ybj
! *** empty log message ***
!
! Revision 1.2  2001/03/03 01:38:12  ybj
! with selected pressure levels
!
! Revision 1.1  2001/02/27 20:52:39  ybj
! FFSM Process
!
! Revision 1.1  2000/10/05 18:17:41  nakamura
! Module split from synoptic.f90 and modified to be more like the standard template.
!
