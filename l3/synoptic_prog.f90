
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============
Program Synoptic
!===============

	Use global_data
	Use main_func
	Use MLSCommon

	implicit none

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &

   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This is the main program to run the Core processing.

! Parameters

	integer, Parameter :: nlons = 180

	! Variable definitions

        integer i, ctype
        real(r8), Dimension(nlons) :: xlon, result
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

!===================
End Program Synoptic
!===================

! $Log$
