
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!==============================================================================
Module global_data
!==============================================================================

  USE MLSCommon, ONLY: r8
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

  Type Loc_T
     Real :: lon
     Real :: t
  End Type Loc_T

  Type RS_T
     Real :: r
     Real :: s
  End Type RS_T

  double complex, Allocatable, Dimension(:) :: phikra, phikrd, phikr
  
  Real(r8), Allocatable, Dimension(:) :: DAcend, Dscend, Ascend, wn, sigma, & 
       & wna, sigmaa, wnd, sigmad, DPrec, APrec
  
  Real, Allocatable, Dimension(:) :: lonD, tD, sD, rD, lonA, tA, sA, rA, &
       & lonDA, tDA, sDA, rDA
  
  Real    :: orbitfreq, c0, lonD0, lonA0, lonDA0, t0, tau0, tD0, tA0, & 
       & tDA0, dtad, dlonad, d1lonad, sina, cosa, ds, sD0, sA0, rD0, rA0, &
       & sDA0, rDA0, krmax, lat 

  Integer	:: nt, nt_a, nt_d, nwave, mtotal, mtotala, mtotald
  
!=====================
End Module global_data
!=====================

! $Log$
! Revision 1.4  2003/03/22 02:56:31  jdone
! use only and indentation added
!
! Revision 1.3  2001/03/03 01:41:56  ybj
! *** empty log message ***
!
! Revision 1.2  2001/02/27 20:53:07  ybj
! global data
!
! Revision 1.1  2000/10/05 18:17:01  nakamura
! Module split from synoptic.f90 and modified to be more like the standard template.
!
