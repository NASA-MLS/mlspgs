
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
Module global_data
!===============================================================================

	USE MLSCommon
	USE SDPTOOLKIT
	Implicit None
	Save
	PUBLIC

	PRIVATE :: ID, ModuleName

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &

   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!----------------------------------------------------------

! Remarks:  This is a prototype module which defines quantities used in the
!           Core processing.

! Parameters

	Integer	:: nt, nt_a, nt_d, nwave, mtotal
	Real    :: orbitfreq, c0, lonD0, lonA0, lonDA0, t0
	Real    :: tau0, tD0, tA0, tDA0, dtad, dlonad, d1lonad, sina, cosa, ds
	Real    :: sD0, sA0, rD0, rA0, sDA0, rDA0
	Real    :: krmax 
	Real    :: lat 

        Real, Allocatable, Dimension(:) :: lonD, tD, sD, rD
        Real, Allocatable, Dimension(:) :: lonA, tA, sA, rA
        !Real, Allocatable, Dimension(:) :: Dscend, Ascend
        Real(r8), Allocatable, Dimension(:) :: Dscend, Ascend, wn, sigma
        double complex, Allocatable, Dimension(:) :: phikr

        Real, Allocatable, Dimension(:) :: lonDA, tDA, sDA, rDA
        Real(r8), Allocatable, Dimension(:) :: DAcend
        double complex, Allocatable, Dimension(:) :: phikra, phikrd

	Type Loc_T
	   Real :: lon
	   Real :: t
	End Type Loc_T

	Type RS_T
	   Real :: r
	   Real :: s
	End Type RS_T

!=====================
End Module global_data
!=====================

! $Log$
! Revision 1.1  2000/10/05 18:17:01  nakamura
! Module split from synoptic.f90 and modified to be more like the standard template.
!
