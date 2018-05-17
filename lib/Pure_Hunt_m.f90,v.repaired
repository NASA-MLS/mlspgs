! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module Pure_Hunt_m
!=============================================================================

  ! A binary search routine with a hunt procedure, to start from last known
  ! location (if 0 < JLO < N) or from the begining otherwise.

  ! It's here instead of in Hunt_m so as not to suck in all the modules
  ! used therein.

  implicit none

  private

  public :: PureHunt

  interface purehunt
    module procedure purehunt_r4, purehunt_r8
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ---------------------------------------------------  PureHunt  -----
  ! A binary search routine with a hunt procedure, to start from last known
  ! location (if 0 < JLO < N) or from the begining otherwise.  Upon return,
  ! Array(JLO) <= Element <= Array(JLO+INC) where 0 <= INC <= 1 if
  ! Array(N) >= Array(1) and -1 <= INC <= 0 otherwise.  INC = 0 if N = 1,
  ! JLO = N and Array(N) >= Array(1), or JLO = 1 and Array(N) < Array(1).
  pure subroutine purehunt_r4 ( ELEMENT, ARRAY, N, JLO, JHI )
    integer, parameter :: RK = kind(0.0e0)
    include 'hunt.f9h'
  end subroutine purehunt_r4

  pure subroutine purehunt_r8 ( ELEMENT, ARRAY, N, JLO, JHI )
    integer, parameter :: RK = kind(0.0d0)
    include 'hunt.f9h'
  end subroutine purehunt_r8

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, id ! .mod files sometimes change if PRINT is added
  end function not_used_here

end module Pure_Hunt_m
!=============================================================================

! $Log$
! Revision 2.2  2018/05/17 01:30:42  vsnyder
! Start hunt with a secant estimate instead of 1:n
!
! Revision 2.1  2016/09/09 00:12:07  vsnyder
! Initial commit
!
