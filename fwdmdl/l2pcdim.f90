module L2PCDIM
  use MLSCommon, only: I4
  implicit NONE
  public
!---------------------------- RCS Ident Info -------------------------------
  private ID, ModuleName
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
!  This file contains all the parameters defining the dimensions of
!  various entities in the L2PC program
!
!  Updated by Z. Shippony on May/5/97.  (For version 5.1)
!  (Increasing Ncf from 22 to: 43)
!  (Increasing Nsps from 13  to: 20)
!  (Commenting out all magnetic stuff ..)
!
!  Updated by Z. Shippony on May/5/98.  (For EOS version)
!  (Adding maximum number of Phi dimension)
!
!  Updated by Z. Shippony on Jan/8/01.  (For EOS version)
!  (Adding maximum number of MMAF per pass)
!
! Integer(i4) :: Maxmagch, Maxmagsps, Maxbfields
  Integer(i4), parameter :: NBAND=6      ! Max. # of bands
  Integer(i4), parameter :: NCH=90       ! Max. # of channels
  Integer(i4), parameter :: NLVL=100     ! Max. # of pre-selected major
!                                          grid points
  Integer(i4), parameter :: N2LVL=2*NLVL ! Twice Nlvl
  Integer(i4), parameter :: NCF=43       ! Max. # of profile coefficients
  Integer(i4), parameter :: NRAD=32      ! Max. # of radiances per channel
  Integer(i4), parameter :: MAXFFT=1024  ! Max. # of FFT points
  Integer(i4), parameter :: MAX_NO_MMAF=10 ! Max. # of MMAF per chunk
!
  Integer(i4), parameter :: NSPS=3       ! Max. # of species
  Integer(i4), parameter :: NPTG=54      ! Max. # of convolution hights
  Integer(i4), parameter :: MAX_NO_PHI=5
!
! Integer(i4), parameter :: MAX_NO_MMAF=70 ! Max. # of MMAF per chunk
! Integer(i4), parameter :: MAX_NO_PHI=13
! Integer(i4), parameter :: NSPS=20      ! Max. # of species
! Integer(i4), parameter :: NSPS=14      ! Max. # of species
! Integer(i4), parameter :: NPTG=67      ! Max. # of convolution hights

end module L2PCDIM
! $Log$
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
