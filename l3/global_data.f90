
! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!==============================================================================
Module global_data
!==============================================================================

  USE MLSCommon, ONLY: r8
  Implicit None
  Save
  PUBLIC


!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! Description:  This module defines quantities used in the Core processing.
  
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
  
contains
!=====================
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here
End Module global_data
!=====================

! $Log$
! Revision 1.6  2005/06/23 19:07:38  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 1.5  2004/05/19 18:08:16  cvuu
! Add parameters DPrec and APrec
!
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
