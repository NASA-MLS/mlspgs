
! /---------------------------------------------------------------\
! /---------------------------------------------------------------\
! |*** Module global_data 
! \---------------------------------------------------------------/
! \---------------------------------------------------------------/

Module global_data

	Implicit None
	Save

	Real, Parameter :: PI=3.14159265

	Integer	:: nt, nwave, mtotal
	Real    :: orbitfreq, c0, lonD0, lonDA0, t0
	Real    :: tau0, tD0, tA0, tDA0, dtad, dlonad, d1lonad, sina, cosa, ds
	Real    :: sD0, sA0, rD0, rA0, sDA0, rDA0
	Real    :: krmax 
	Real    :: lat 

        Real, Allocatable, Dimension(:) :: lonD, tD, sD, rD
        Real, Allocatable, Dimension(:) :: lonA, tA, sA, rA
        !Real, Allocatable, Dimension(:) :: Dscend, Ascend
        double precision, Allocatable, Dimension(:) :: Dscend, Ascend, wn, sigma
        double complex, Allocatable, Dimension(:) :: phikr

        Real, Allocatable, Dimension(:) :: lonDA, tDA, sDA, rDA
        double precision, Allocatable, Dimension(:) :: DAcend
        double complex, Allocatable, Dimension(:) :: phikrda

	Type Loc_T
	   Real :: lon
	   Real :: t
	End Type Loc_T

	Type RS_T
	   Real :: r
	   Real :: s
	End Type RS_T

End Module global_data


! /---------------------------------------------------------------\
! /---------------------------------------------------------------\
! |*** Module main_func 
! \---------------------------------------------------------------/
! \---------------------------------------------------------------/

Module main_func
	
	Use global_data
	Implicit None

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

	Subroutine Init(ctype, nt_i, orbitfreq_i, c0_i, lonDA0_i, tDA0_i, lat_i)

	  Implicit None

          Integer :: nt_i, ctype
	  Real    :: orbitfreq_i, c0_i, lonDA0_i, tDA0_i, lat_i
	 
          lat	    = lat_i
	  nt        = nt_i
 	  orbitfreq = orbitfreq_i
	  c0 	    = c0_i
 	  tau0 	    = 1/orbitfreq

 	  sina      = 1.0/sqrt(1.0+c0*c0)
 	  cosa      = c0/sqrt(1.0+c0*c0)
 	  ds        = tau0/sina

	  krmax     = 1.0/sina
 	  nwave     = PI/(ds*cosa)

 	  dtad	    = tau0/2.0
 	  dlonad    = PI-dtad*c0 

	  d1lonad   = dlonad + c0*dtad

	  if(ctype == 0) then
	    lonD0   = lonDA0_i
    	    tD0     = tDA0_i
    	    tA0     = tD0-dtad
          else
	    lonDA0   = lonDA0_i
    	    tDA0     = tDA0_i
	  end if

	  if(ctype == 0) then
	    Allocate(lonD(nt), tD(nt), sD(nt), rD(nt))
	    Allocate(lonA(nt), tA(nt), sA(nt), rA(nt))
	    Allocate(Dscend(nt), Ascend(nt))
	    Allocate(wn(nt*2), sigma(nt*2))
	    Allocate(phikr(nt*2))
          else
	    Allocate(lonDA(nt), tDA(nt), sDA(nt), rDA(nt))
	    Allocate(DAcend(nt))
	    Allocate(wn(nt*2), sigma(nt*2))
	    Allocate(phikrda(nt*2))
	  end if
 	 
	End Subroutine Init

! /---------------------------------------------------------------\
! |*** Coordinate transformation
! \---------------------------------------------------------------/

        Subroutine Loc2RS(aLoc, aRS)

	 TYPE(Loc_T), INTENT(IN) :: aLoc
	 TYPE(RS_T), INTENT(OUT) :: aRS

	End Subroutine Loc2RS


        Subroutine cordTransform(ctype)

	  Implicit None

	  Integer :: i, ctype

          if(ctype == 0) then
	    sD0 = lonD0*cosa - tD0*sina - tau0*(nt-1.0)/sina
	    sA0 = sD0 + (-dlonad*c0 + dtad)*sina
	    rD0 = (lonD0 + c0*tD0)*sina
	    rA0 = (lonD0 - dlonad + c0*tD0 - c0*dtad)*sina
	  else
	    sDA0 = lonDA0*cosa - tDA0*sina - tau0*(nt-1.0)/sina
	    rDA0 = (lonDA0 + c0*tDA0)*sina
          end if 

          if(ctype == 0) then
	    Do i = 1, nt
   	       rD(i) = rD0
               rA(i) = rA0
               sD(i) = sD0 + i*tau0/sina
               sA(i) = sA0 + i*tau0/sina
            End Do
	  else
	    Do i = 1, nt
   	       rDA(i) = rDA0
               sDA(i) = sDA0 + i*tau0/sina
            End Do
          end if 

	End Subroutine cordTransform


! /---------------------------------------------------------------\
! |*** Data Field 
! \---------------------------------------------------------------/

	Real Function DataField(lon, lat, t)

	  Implicit None

	  Real :: lon, lat, t, m0, sigma0, k0, DataField

   	  m0      = 4.0
   	  sigma0  = 0.1
   	  !sigma0 = 0.0
   	  k0      = 2.25

   	  !DataField = cos(m0*lon + sigma0*t)
   	  !DataField = cos(2.0*lon + sigma0*t)
   	  !DataField = 2.0*cos(m0*lon - sigma0*t + lat) + cos(3.0*lon - 2.0*sigma0*t + lat) - 1.5*cos(1.0*lon - 3.0*sigma0*t)
   	  !DataField = 2.0*cos(4*lon - 2.*sigma0*t + lat) - cos(3.0*lon - 2.0*sigma0*t + lat) - 1.5*cos(1.0*lon + 1.0*sigma0*t)
   	  !DataField = cos(m0*lon + sigma0*t)*cos(k0*lat)
   	  DataField = 1.5 + 0.5*cos(lon) + cos(2.0*lon + 0.5*t)

	End Function DataField


! /---------------------------------------------------------------\
! |*** findad 
! \---------------------------------------------------------------/

	Subroutine findad(sina, cosa, t, r0, ep, em, frp, frm, rpi, rmi)

	  Implicit None

	  Integer :: inc
	  Real :: sina, cosa, t, r0, ep, em, frp, frm, rpi, rmi
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
! |*** FFSM --- combined mode
! \---------------------------------------------------------------/

        Subroutine FFSM()

	  Implicit None

	  integer plan, m, m1, i, j, flag

          double precision, Dimension(nt) :: dscend_fft, ascend_fft
          double precision, Dimension(nt) :: drealP, dimgP, arealP, aimgP 
          double precision, Dimension((nt+1)/2-1) :: imgPtemp
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

   	  aimgP(1) = 0.0 
	  Do i = nt/2+1, nt
   	     aimgP(i) = imgPtemp(i-nt/2)
          End Do
	  Do i = 2, nt/2
   	     aimgP(i) = -imgPtemp(nt/2-i)
          End Do

! Dscending Part

	  Do i = nt/2, nt
   	     drealP(i) = dscend_fft(i+1-nt/2)/float(nt)
          End Do
	  Do i = 1, nt/2-1
   	     drealP(i) = dscend_fft(nt/2-i+1)/float(nt)
          End Do

	  Do i = 1, (nt+1)/2-1
   	     imgPTemp(i) = dscend_fft(nt+1-i)/float(nt)
          End Do

   	  dimgP(1) = 0.0 
	  Do i = nt/2+1, nt
   	     dimgP(i) = imgPtemp(i-nt/2)
          End Do
	  Do i = 2, nt/2
   	     dimgP(i) = -imgPtemp(nt/2-i)
          End Do

          open(2, file="test.dat", status="replace")
          write(2, *) ds
          write(2, *) c0
          write(2, *) sina
          write(2, *) cosa

   	  Do i = 1, nt
	 	!write(2, *) i, dscend(i), dscend_fft(i)
	  End do

 	  mtotal = 0

	  print *, 'nwave=', nwave

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

	 write(2, *) mtotal
	 do j = 1, mtotal
	   write(2, *) j, wn(j), sigma(j), real(phikr(j)), aimag(phikr(j))
	   !write(2, *) j, wn(j), sigma(j), phikr(j)
         end do

	 close(2)


	End Subroutine FFSM


! /---------------------------------------------------------------\
! |*** FFSM1 --- single mode
! \---------------------------------------------------------------/

        Subroutine FFSM1()

	  Implicit None

	  integer plan, m, m1, i, j, flag

          double precision, Dimension(nt) :: dacend_fft
          double precision, Dimension(nt) :: darealP, daimgP
          double precision, Dimension((nt+1)/2-1) :: imgPtemp
	  double complex phida

	  real	ks, kr, expda, krmax_g


! Fourier Transform of one series

	  call rfftw_f77_create_plan(plan,nt,FFTW_FORWARD,FFTW_ESTIMATE)
	  call rfftw_f77_one(plan,dacend,dacend_fft)
	  call rfftw_f77_destroy_plan(plan)

! re-order the fft transform result

	  Do i = nt/2, nt
   	     darealP(i) = dacend_fft(i+1-nt/2)/float(nt)
          End Do
	  Do i = 1, nt/2-1
   	     darealP(i) = dacend_fft(nt/2-i+1)/float(nt)
          End Do

	  Do i = 1, (nt+1)/2-1
   	     imgPTemp(i) = dacend_fft(nt+1-i)/float(nt)
          End Do

   	  daimgP(1) = 0.0 
	  Do i = nt/2+1, nt
   	     daimgP(i) = imgPtemp(i-nt/2)
          End Do
	  Do i = 2, nt/2
   	     daimgP(i) = -imgPtemp(nt/2-i)
          End Do

          open(2, file="test.dat", status="replace")
          write(2, *) ds
          write(2, *) c0
          write(2, *) sina
          write(2, *) cosa

   	  Do i = 1, nt
	 	!write(2, *) i, dscend(i), dscend_fft(i)
	  End do

 	  mtotal = 0

	  print *, 'nwave=', nwave

	  krmax_g = krmax/2.0

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

               expda = ks*sda0 + kr*rda0

               flag = 0

     	       phida   = CMPLX(darealP(i), daimgP(i))*CMPLX(cos(expda), -sin(expda))

     	       phikrda(mtotal) = phida 

 
      	       if( phikrda(mtotal) /= CMPLX(0.0, 0.0) ) flag = 1


101            continue 

! *** kr+

      	       if (kr > krmax_g .or. kr < 0.0 .or. ks < 0.0 ) then 
                  if (flag == 0 .and. mtotal > 0) mtotal = mtotal-1
                  Goto 102
     	       endif

               kr  = -ks*c0+float(m1)/sina

     	       expda = ks*sda0 + kr*rda0

     	       phida   = CMPLX(darealP(i), daimgP(i))*CMPLX(cos(expda), -sin(expda))

     	       phikrda(mtotal) = phida 

102            continue 

           End Do
         End Do

	 write(2, *) mtotal
	 do j = 1, mtotal
	   write(2, *) j, wn(j), sigma(j), real(phikrda(j)), aimag(phikrda(j))
         end do

	 close(2)


	End Subroutine FFSM1


! /---------------------------------------------------------------\
! |*** Reconstruct 
! \---------------------------------------------------------------/

        Subroutine Reconstruct(ctype, xtime, nlons, xlon, result)

	  Implicit None


	  integer plan, m, m1, i, j, k, ctype, nlons

	  real	kr, ks

	  real  xtime, xloni, eoff, rp, rm, argp, argm, ep, em, sum 
	  real  argpa, argpd, argma, argmd, argp, argm
 	  real  epa, epd, ema, emd
 	  real  rma, rmd, rpa, rpd , rdp, ran, rdn
 	  real  frpa, frma, frpd, frmd, frp, frm, rpi, rmi
          double precision, Dimension(nlons) :: xlon, result 


!*** find offset in longitude necessary to get r value right

 	 eoff = c0*xtime - mod( (c0*xtime), (2.0*PI) )
 	 if (eoff < 0.0) eoff = 2.0*PI

!*** Loop over longitude

 	 do k = 1, nlons

  	    result(k) = 0.0

            xloni = xlon(k) - eoff

            if (ctype == 0) then 

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


 	    else 

  		ep = xloni
  		em = xloni

  		call findad( sina, cosa, xtime, rda0, ep, em, frp, frm, rpi, rmi)

  		rp = rpi 
  		rm = rmi 

  		do j = 1, mtotal 

        	  kr     = wn(j)*sina+sigma(j)*cosa
        	  ks     = wn(j)*cosa-sigma(j)*sina

        	  argp = (wn(j)*ep + sigma(j)*xtime) + kr*(rda0-rp)
        	  argm = (wn(j)*em + sigma(j)*xtime) + kr*(rda0-rm)

		  sum = real( frp*phikrda(j)*CMPLX(cos(argp), sin(argp)) ) + &
	      	        real( frm*phikrda(j)*CMPLX(cos(argm), sin(argm)) ) 

        	  if( abs(ks) > 1.e-4) then 
	   		sum = sum + real( conjg(frp*phikrda(j))*CMPLX(cos(argp), -sin(argp)) ) + &
	               		    real( conjg(frm*phikrda(j))*CMPLX(cos(argm), -sin(argm)) ) 
        	  endif

  		  result(k) = result(k) + sum

  		end do


 	    end if


 	 end do


 	close(3)
 	close(4)



	End Subroutine Reconstruct


! /---------------------------------------------------------------\
! |*** Diagnostics 
! \---------------------------------------------------------------/

        Subroutine Diagnostics(ctype, xtime, nlons, xlon, result)

	  Implicit None


	  integer plan, m, m1, i, j, k, ctype, nlons

	  real	kr, ks

	  real  xtime, xloni, eoff, rp, rm, argp, argm, ep, em, sum 
	  real  argpa, argpd, argma, argmd
 	  real  epa, epd, ema, emd
 	  real  rma, rmd, rpa, rpd , rdp, ran, rdn
 	  real  frpa, frma, frpd, frmd, frp, frm, rpi, rmi
          double precision, Dimension(nlons) :: xlon, result 


	End Subroutine Diagnostics


! /---------------------------------------------------------------\
! |*** Generate Data
! \---------------------------------------------------------------/

        Subroutine DataGenerate(filename1, filename2)

	 Implicit None

	 Integer i
	 Character(*) :: filename1, filename2
         Real	      :: temp1, temp2, temp3

         lonD(1) = lonD0
         lonA(1) = lonD0 - dlonad
         tD(1)   = tD0
         tA(1)   = tA0

         If (lonD(1) < -PI) lonD(1) = 2.0*PI + lonD(1)
   	 print *, 'lonD0=', lonD0

         If (lonA(1) < -PI) lonA(1) = 2.0*PI + lonA(1)

	 Do i = 2, nt
            lonD(i) = lonD(i-1) - c0*tau0
            If (lonD(i) < -PI) lonD(i) = 2.0*PI + lonD(i)
            lonA(i) = lonD(i) - dlonad
            If (lonA(i) < -PI) lonA(i) = 2.0*PI + lonA(i)
            tD(i) = tD0 + i*tau0
            tA(i) = tA0 + i*tau0
            temp1 = lonD(i)
            temp2 = tD(i)
            Dscend(nt+1-i) = DataField(temp1, lat, temp2) 
            temp1 = lonA(i)
            temp2 = tA(i)
            Ascend(nt+1-i) = DataField(temp1, lat, temp2) 
         End Do
         temp1 = lonD(1)
         temp2 = tD(1)
         Dscend(nt) = DataField(temp1, lat, temp2) 
         temp1 = lonA(1)
         temp2 = tA(1)
         Ascend(nt) = DataField(temp1, lat, temp2) 

         open(2, file=filename1, status="replace")
   	 write(2, '(4A10)') 'lonA', 'lonD', 'tA', 'tD'
   	 write(2, *) nt 
	 Do i = 1, nt
   	   write(2, '(4F10.5)') lonA(i), lonD(i), tA(i), tD(i)
         End Do
         close(2) 

         open(2, file=filename2, status="replace")
   	 write(2, '(6A10)') 'sA', 'sD', 'rA', 'rD', 'Ascend', 'Dscend'  
   	 write(2, *) nt 
	 Do i = 1, nt
   	   write(2, '(6F10.5)') sA(i), sD(i), rA(i), rD(i), Ascend(i), Dscend(i)  
         End Do
         close(2) 

         open(2, file="ascend.dat", status="replace")
   	 write(2, *) nt 
	 Do i = 1, nt
   	   write(2, '(F10.5)') Ascend(i) 
         End Do
         close(2) 

         open(2, file="dscend.dat", status="replace")
   	 write(2, *) nt 
	 Do i = 1, nt
   	   write(2, '(F10.5)') Dscend(i) 
         End Do
         close(2) 

	End Subroutine DataGenerate


! /---------------------------------------------------------------\
! |*** Generate Data for only one series
! \---------------------------------------------------------------/

        Subroutine DataGenerate1(filename1, filename2)

	 Implicit None

	 Integer i
	 Character(*) :: filename1, filename2
         Real	      :: temp1, temp2, temp3

         lonDA(1) = lonDA0
         tDA(1)   = tDA0

         If (lonDA(1) < -PI) lonDA(1) = 2.0*PI + lonDA(1)
   	 print *, 'lonDA0=', lonDA0

	 Do i = 2, nt
            lonDA(i) = lonDA(i-1) - c0*tau0
            If (lonDA(i) < -PI) lonDA(i) = 2.0*PI + lonDA(i)
            tDA(i) = tDA0 + i*tau0
            temp1 = lonDA(i)
            temp2 = tDA(i)
            DAcend(nt+1-i) = DataField(temp1, lat, temp2) 
         End Do
         temp1 = lonDA(1)
         temp2 = tDA(1)
         DAcend(nt) = DataField(temp1, lat, temp2) 

         open(2, file=filename1, status="replace")
   	 write(2, '(4A10)') 'lonDA',  'tDA'
   	 write(2, *) nt 
	 Do i = 1, nt
   	   write(2, '(2F10.5)') lonDA(i), tDA(i)
         End Do
         close(2) 

         open(2, file=filename2, status="replace")
   	 write(2, '(4A10)') 'sDA', 'rDA', 'DAcend'  
   	 write(2, *) nt 
	 Do i = 1, nt
   	   write(2, '(4F10.5)') sDA(i), rDA(i), DAcend(i)  
         End Do
         close(2) 

	End Subroutine DataGenerate1



End Module main_func



! /---------------------------------------------------------------\
! /---------------------------------------------------------------\
! |*** Main Program 
! \---------------------------------------------------------------/
! \---------------------------------------------------------------/


Program Synoptic

	Use global_data
	Use main_func

	implicit none
	
	! Variable definitions

	integer, Parameter :: nlons = 180 
        integer i, ctype
        double precision, Dimension(nlons) :: xlon, result
	real xloni, xtemp, xtime
 
        ctype = 1
	xtime = 5.5

 	do i = 1, nlons 
  		xlon(i) = ((i-1)*360.0/float(nlons)*PI/180.0 - PI)
 	end do

	!call Init(nt_i=180, orbitfreq_i=32.0, c0_i=2.0*PI, lonD0_i=PI, tD0_i=0.0, lat_i=0.0)
	!call Init(nt_i=128, orbitfreq_i=20.0, c0_i=2.0*PI, lonD0_i=0.8, tD0_i=0.0, lat_i=.0)
	!call Init(nt_i=128, orbitfreq_i=15.0, c0_i=2.0*PI, lonD0_i=0.8, tD0_i=0.0, lat_i=.0)
	call Init(ctype, nt_i=128, orbitfreq_i=20.0, c0_i=2.0*PI, lonDA0_i=0.8, tDA0_i=0.0, lat_i=0.0)

        print *, orbitfreq, c0, lonD0, tD0

        call cordTransform(ctype)
        if(ctype == 0) then
          call DataGenerate("cord.dat", "cord_rs.dat")
 	else
          call DataGenerate1("cord.dat", "cord_rs.dat")
 	end if

        if(ctype == 0) then
	   call FFSM()
 	else
	   call FFSM1()
 	end if

        call Reconstruct(ctype, xtime, nlons, xlon, result)

        open(3, file="synoptic.dat", status="replace")
	write(3, *) nlons
 	do i = 1, nlons
           xloni = xlon(i)
	   xtemp = DataField(xloni, 0.0, xtime) 
	   write(3, '(I5, 4F10.5)') i, xtemp, xlon(i), result(i)
        end do
        close(3)

End Program Synoptic
!=============================================================================

!
! $Log$
! Revision 1.1  2000/10/05 16:37:19  ybj
! L3 Core File (test) 

