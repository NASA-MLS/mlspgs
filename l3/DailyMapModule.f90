
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
Module DailyMapModule
!===============================================================================
	
   	Use MLSCommon
   	USE MLSCF
   	USE MLSPCF3
   	USE L3CF
   	USE L3DMData
   	USE L2GPData
   	USE L2Interface
   	USE L3DMData
   	USE L3SPData
	Use global_data
	Implicit None
	PUBLIC

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
!                DataGenerate1
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

	Subroutine Init(mode, nt_a_i, nt_d_i, tau0_i, delTad_i, 	&
			c0_i, lonD0_i, tD0_i, lonA0_i, tA0_i, lat_i)

! Arguments

	  CHARACTER (LEN=3) :: mode                  ! asc/des/com/all

          Integer :: nt_a_i, nt_d_i
	  Real    :: c0_i, tau0_i
	  Real (r8)   :: lonD0_i, tD0_i, lonA0_i, tA0_i, lat_i, delTad_i
	 
          lat	    = lat_i
	  nt_a      = nt_a_i
	  nt_d      = nt_d_i

	  if(nt_a > nt_d) then
	    nt = nt_d
	  else
	    nt = nt_a
	  end if

	  nt_a = nt
	  nt_d = nt

	  c0 	    = c0_i
 	  tau0 	    = tau0_i 

 	  sina      = 1.0/sqrt(1.0+c0*c0)
 	  cosa      = c0/sqrt(1.0+c0*c0)
 	  ds        = tau0/sina

	  krmax     = 1.0/sina
 	  nwave     = PI/(ds*cosa)

 	  dtad	    = delTad_i 
 	  !dlonad    = PI-dtad*c0 
 	  dlonad    = abs(lonD0_i-lonA0_i) 

	  d1lonad   = dlonad + c0*dtad

	  lonD0   = lonD0_i
    	  tD0     = tD0_i
	  lonA0   = lonA0_i
    	  tA0     = tA0_i

!	  if(mode == 'com') then
	    Allocate(lonD(nt), tD(nt), sD(nt), rD(nt))
	    Allocate(lonA(nt), tA(nt), sA(nt), rA(nt))
	    Allocate(Dscend(nt_d), Ascend(nt_a))
	    Allocate(wn(nt*2), sigma(nt*2))
	    Allocate(wna(nt_a*2), sigmaa(nt_a*2))
	    Allocate(wnd(nt_d*2), sigmad(nt_d*2))
	    Allocate(phikr(nt*2))
	    Allocate(phikra(nt_a*2))
	    Allocate(phikrd(nt_d*2))
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

	    DeAllocate(lonD, tD, sD, rD)
	    DeAllocate(lonA, tA, sA, rA)
	    DeAllocate(Dscend, Ascend)
	    DeAllocate(wn, sigma)
	    DeAllocate(wna, sigmaa)
	    DeAllocate(wnd, sigmad)
	    DeAllocate(phikr)
	    DeAllocate(phikra)
	    DeAllocate(phikrd)
 	 
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

	  Integer :: inc
	  Real :: e, r, s

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

   	  !TYPE( L3SPData_T ), POINTER :: l3sp
   	  TYPE( L3SPData_T ) :: l3sp

	  integer plan, m, m1, i, j, flag, iLv, iLt

	  integer fNum(0:nt), wIndex

          real(r8), Dimension(nt) :: dscend_fft, ascend_fft
          real(r8), Dimension(nt) :: drealP, dimgP, arealP, aimgP 
          real(r8), Dimension((nt+1)/2) :: imgPtemp
	  double complex phia, phid

	  real	ks, kr, expa, expd, krmax_g

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
	  Do i = nt/2+1, nt
   	     aimgP(i) = imgPtemp(i-nt/2)
          End Do
	  Do i = 1, nt/2-1
   	     aimgP(i) = -imgPtemp(nt/2-i)
          End Do

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
	  Do i = nt/2+1, nt
   	     dimgP(i) = imgPtemp(i-nt/2)
          End Do
	  Do i = 1, nt/2-1
   	     dimgP(i) = -imgPtemp(nt/2-i)
          End Do

 	  mtotal = 0

	  krmax_g = krmax

 	  Do m = 1, nwave 
            m1 = m-1
   	    Do i = 1, nt

     	       mtotal = mtotal + 1

      	       wn(mtotal) = float(m1)

               ks = 2.0*PI*(-nt/2+i)/(nt*ds)
               kr = -ks*c0+float(m1)/sina

               sigma(mtotal) = -ks*sina + kr*cosa

               if (kr < -krmax_g .or. kr > 0.0 .or. ks < 0.0) then 
                  flag = 0
                  Goto 101
               endif

               expa = ks*sa0 + kr*ra0
               expd = ks*sd0 + kr*rd0

               flag = 0

     	       phia   = CMPLX(arealP(i), aimgP(i))*CMPLX(cos(expa), -sin(expa))
     	       phid   = CMPLX(drealP(i), dimgP(i))*CMPLX(cos(expd), -sin(expd))

     	       phikr(mtotal) = (phia-phid*CMPLX(cos(d1lonad), -sin(d1lonad)))/(1.0-CMPLX(cos(d1lonad), -sin(d1lonad)))

 
      	       if( phikr(mtotal) /= CMPLX(0.0, 0.0) ) flag = 1


101            continue 

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

     	       phikr(mtotal) = (phid-phia)*CMPLX( cos(rd0/sina), -sin(rd0/sina) )/(1.0-CMPLX( cos(d1lonad), -sin(d1lonad) ))

102            continue 

           End Do
         End Do

	 do j = 1, nt
	   fNum(j) = 0 
	 end do

	 do j = 1, mtotal
	   wIndex = int(wn(j))
	   fNum(wIndex) = fNum(wIndex) + 1
	   l3sp%waveNumber(iLv, iLt, wIndex+1) = wIndex 
	   l3sp%frequency(iLv, iLt, fNum(wIndex)) = sigma(j)
	   l3sp%l3spRelValue(iLv, iLt, wIndex+1, fNum(wIndex)) = real(phikr(j))
	   l3sp%l3spImgValue(iLv, iLt, wIndex+1, fNum(wIndex)) = aimag(phikr(j))
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

	  integer plan, m, m1, i, j, flag, iLv, iLt

	  integer fNum(0:nt_a), wIndex

          real(r8), Dimension(nt_a) :: ascend_fft
          real(r8), Dimension(nt_a) :: arealP, aimgP
          real(r8), Dimension((nt_a+1)/2) :: imgPtemp
	  double complex phida

	  real	ks, kr, expda, krmax_g


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

               ks = 2.0*PI*(-nt_a/2+i)/(nt_a*ds)
               kr = -ks*c0+float(m1)/sina

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


101            continue 

! *** kr+

      	       if (kr > krmax_g .or. kr < 0.0 .or. ks < 0.0 ) then 
                  if (flag == 0 .and. mtotala > 0) mtotala = mtotala-1
                  Goto 102
     	       endif

               kr  = -ks*c0+float(m1)/sina

     	       expda = ks*sa0 + kr*ra0

     	       phida   = CMPLX(arealP(i), aimgP(i))*CMPLX(cos(expda), -sin(expda))

     	       phikra(mtotala) = phida 

102            continue 

           End Do
         End Do

	 do j = 1, nt_a
	   fNum(j) = 0 
	 end do

	 do j = 1, mtotala
	   wIndex = int(wna(j))
	   fNum(wIndex) = fNum(wIndex) + 1
	   l3sp%waveNumber(iLv, iLt, wIndex+1) = wIndex 
	   l3sp%frequency(iLv, iLt, fNum(wIndex)) = sigmaa(j)
	   l3sp%l3spRelValue(iLv, iLt, wIndex+1, fNum(wIndex)) = real(phikra(j))
	   l3sp%l3spImgValue(iLv, iLt, wIndex+1, fNum(wIndex)) = aimag(phikra(j))
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

	  integer plan, m, m1, i, j, flag, iLv, iLt

	  integer fNum(0:nt_d), wIndex

          real(r8), Dimension(nt_d) :: dscend_fft
          real(r8), Dimension(nt_d) :: drealP, dimgP
          real(r8), Dimension((nt_d+1)/2) :: imgPtemp
	  double complex phida

	  real	ks, kr, expda, krmax_g


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

               ks = 2.0*PI*(-nt_d/2+i)/(nt_d*ds)
               kr = -ks*c0+float(m1)/sina

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


101            continue 

! *** kr+

      	       if (kr > krmax_g .or. kr < 0.0 .or. ks < 0.0 ) then 
                  if (flag == 0 .and. mtotald > 0) mtotald = mtotald-1
                  Goto 102
     	       endif

               kr  = -ks*c0+float(m1)/sina

     	       expda = ks*sd0 + kr*rd0

     	       phida   = CMPLX(drealP(i), dimgP(i))*CMPLX(cos(expda), -sin(expda))

     	       phikrd(mtotald) = phida 

102            continue 

           End Do
         End Do

	 do j = 1, nt_d
	   fNum(j) = 0 
	 end do

	 do j = 1, mtotald
	   wIndex = int(wnd(j))
	   fNum(wIndex) = fNum(wIndex) + 1
	   l3sp%waveNumber(iLv, iLt, wIndex+1) = wIndex 
	   l3sp%frequency(iLv, iLt, fNum(wIndex)) = sigmad(j)
	   l3sp%l3spRelValue(iLv, iLt, wIndex+1, fNum(wIndex)) = real(phikrd(j))
	   l3sp%l3spImgValue(iLv, iLt, wIndex+1, fNum(wIndex)) = aimag(phikrd(j))
         end do


	End Subroutine FFSMD


! /---------------------------------------------------------------\
! |*** Reconstruct 
! \---------------------------------------------------------------/

        Subroutine Reconstruct(mode, xtime, nlons, xlon_orig, result)

! Arguments

	  CHARACTER (LEN=3) :: mode                  ! asc/des/com/all

	  integer nlons

	  real xtime

	  real(r8), Dimension(nlons) :: xlon, xlon_orig, result

! Variables

	  integer j, k

	  real  kr, ks

	  real  xloni, eoff, rp, rm, ep, em, sum 
	  real  argpa, argpd, argma, argmd, argp, argm
 	  real  epa, epd, ema, emd
 	  real  rma, rmd, rpa, rpd , rdp, ran, rdn
 	  real  frpa, frma, frpd, frmd, frp, frm, rpi, rmi

!*** Convert to gradient 

 	  do k = 1, nlons
	    xlon(k) = xlon_orig(k)*PI/180.0
  	  end do

!*** find offset in longitude necessary to get r value right

 	 eoff = c0*xtime - mod( (c0*xtime), (2.0*PI) )
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

        	  argpa = (wn(j)*epa + sigma(j)*xtime) + kr*(ra0-rpa)
        	  argma = (wn(j)*ema + sigma(j)*xtime) + kr*(ra0-rma)
        	  argpd = (wn(j)*epd + sigma(j)*xtime) + kr*(rd0-rpd)
        	  argmd = (wn(j)*emd + sigma(j)*xtime) + kr*(rd0-rmd)

		  sum = real( frpa*phikr(j)*CMPLX(cos(argpa), sin(argpa)) ) + &
	      		real( frma*phikr(j)*CMPLX(cos(argma), sin(argma)) ) + &
	      		real( frpd*phikr(j)*CMPLX(cos(argpd), sin(argpd)) ) + &
	      		real( frmd*phikr(j)*CMPLX(cos(argmd), sin(argmd)) ) 

        	  if( abs(ks) > 1.e-4) then 
	   		sum = sum + real( conjg(frpa*phikr(j))*CMPLX(cos(argpa), -sin(argpa)) ) + &
	               		    real( conjg(frma*phikr(j))*CMPLX(cos(argma), -sin(argma)) ) + &
	               		    real( conjg(frpd*phikr(j))*CMPLX(cos(argpd), -sin(argpd)) ) + &
	               		    real( conjg(frmd*phikr(j))*CMPLX(cos(argmd), -sin(argmd)) ) 
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

        	  argp = (wnd(j)*ep + sigmad(j)*xtime) + kr*(rd0-rp)
        	  argm = (wnd(j)*em + sigmad(j)*xtime) + kr*(rd0-rm)

		  sum = real( frp*phikrd(j)*CMPLX(cos(argp), sin(argp)) ) + &
	      	        real( frm*phikrd(j)*CMPLX(cos(argm), sin(argm)) ) 

        	  if( abs(ks) > 1.e-4) then 
	   		sum = sum + real( conjg(frp*phikrd(j))*CMPLX(cos(argp), -sin(argp)) ) + &
	               		    real( conjg(frm*phikrd(j))*CMPLX(cos(argm), -sin(argm)) ) 
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

        	  argp = (wna(j)*ep + sigmaa(j)*xtime) + kr*(ra0-rp)
        	  argm = (wna(j)*em + sigmaa(j)*xtime) + kr*(ra0-rm)

		  sum = real( frp*phikra(j)*CMPLX(cos(argp), sin(argp)) ) + &
	      	        real( frm*phikra(j)*CMPLX(cos(argm), sin(argm)) ) 

        	  if( abs(ks) > 1.e-4) then 
	   		sum = sum + real( conjg(frp*phikra(j))*CMPLX(cos(argp), -sin(argp)) ) + &
	               		    real( conjg(frm*phikra(j))*CMPLX(cos(argm), -sin(argm)) ) 
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

	  real (r8) :: xtime_in, xlon_in
	  real xtime, xlon, result

! Variables

	  integer j

	  real  kr, ks

	  real  xloni, eoff, rp, rm, ep, em, sum 
	  real  argpa, argpd, argma, argmd, argp, argm
 	  real  epa, epd, ema, emd
 	  real  rma, rmd, rpa, rpd , rdp, ran, rdn
 	  real  frpa, frma, frpd, frmd, frp, frm, rpi, rmi

!*** find offset in longitude necessary to get r value right

	 xtime = xtime_in
	 xlon = xlon_in

 	 eoff = c0*xtime - mod( (c0*xtime), (2.0*PI) )
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

        	  argpa = (wn(j)*epa + sigma(j)*xtime) + kr*(ra0-rpa)
        	  argma = (wn(j)*ema + sigma(j)*xtime) + kr*(ra0-rma)
        	  argpd = (wn(j)*epd + sigma(j)*xtime) + kr*(rd0-rpd)
        	  argmd = (wn(j)*emd + sigma(j)*xtime) + kr*(rd0-rmd)

		  sum = real( frpa*phikr(j)*CMPLX(cos(argpa), sin(argpa)) ) + &
	      		real( frma*phikr(j)*CMPLX(cos(argma), sin(argma)) ) + &
	      		real( frpd*phikr(j)*CMPLX(cos(argpd), sin(argpd)) ) + &
	      		real( frmd*phikr(j)*CMPLX(cos(argmd), sin(argmd)) ) 

        	  if( abs(ks) > 1.e-4) then 
	   		sum = sum + real( conjg(frpa*phikr(j))*CMPLX(cos(argpa), -sin(argpa)) ) + &
	               		    real( conjg(frma*phikr(j))*CMPLX(cos(argma), -sin(argma)) ) + &
	               		    real( conjg(frpd*phikr(j))*CMPLX(cos(argpd), -sin(argpd)) ) + &
	               		    real( conjg(frmd*phikr(j))*CMPLX(cos(argmd), -sin(argmd)) ) 
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

        	  argp = (wnd(j)*ep + sigmad(j)*xtime) + kr*(rd0-rp)
        	  argm = (wnd(j)*em + sigmad(j)*xtime) + kr*(rd0-rm)

		  sum = real( frp*phikrd(j)*CMPLX(cos(argp), sin(argp)) ) + &
	      	        real( frm*phikrd(j)*CMPLX(cos(argm), sin(argm)) ) 

        	  if( abs(ks) > 1.e-4) then 
	   		sum = sum + real( conjg(frp*phikrd(j))*CMPLX(cos(argp), -sin(argp)) ) + &
	               		    real( conjg(frm*phikrd(j))*CMPLX(cos(argm), -sin(argm)) ) 
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

        	  argp = (wna(j)*ep + sigmaa(j)*xtime) + kr*(ra0-rp)
        	  argm = (wna(j)*em + sigmaa(j)*xtime) + kr*(ra0-rm)

		  sum = real( frp*phikra(j)*CMPLX(cos(argp), sin(argp)) ) + &
	      	        real( frm*phikra(j)*CMPLX(cos(argm), sin(argm)) ) 

        	  if( abs(ks) > 1.e-4) then 
	   		sum = sum + real( conjg(frp*phikra(j))*CMPLX(cos(argp), -sin(argp)) ) + &
	               		    real( conjg(frm*phikra(j))*CMPLX(cos(argm), -sin(argm)) ) 
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

         Real (r8), DIMENSION(:) ::  aField, dField


! Variables

	 Integer i

	 Do i = 1, nt
            Dscend(nt+1-i) = dField(i) 
            Ascend(nt+1-i) = aField(i) 
         End Do


	End Subroutine DataGenerate

!===================
End Module DailyMapModule
!===================

! $Log$
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
