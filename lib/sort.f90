!     .  Copyright (C) 1989, California Institute of Technology.
!     .  All rights reserved.  U. S. Government sponsorship under
!     .  NASA contract NAS7-918 is acknowledged.
module Sort_M

!     Sort the M:N-vector A into ascending order.
!
!     To sort an array A' into descending order, let A = -A'
!     To sort an array A' into ascending order according to the
!     absolute value of the elements let A = ABS(A').
!     To sort an array A' into decending order according to the
!     absolute value of the elements let A = -ABS(A').
!
!     To keep track of the original elements, use SORTP.

!>> 2002-01-30 SORT   Snyder Convert to Fortran 90
!>> 1995-11-17 SSORT  Krogh  Converted SFTRAN to Fortran 77
!>> 1994-10-19 SSORT  Krogh  Changes to use M77CON
!>> 1988-11-22 SSORT  Snyder Initial code.

  interface SORT
    module procedure DSORT, ISORT, SSORT
  end interface

  interface SORTP ! Using a permutation vector
    module procedure DSORTP, ISORTP, SSORTP
  end interface

  interface SORTQ ! Using a pre-specified permutation vector
    module procedure DSORTQ, ISORTQ, SSORTQ
  end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine DSORT (A, M, N)
    double precision, intent(inout) :: A(*)
    double precision :: PARTN, TEMP
    include "sort.f9h"
  end subroutine DSORT

  subroutine ISORT (A, M, N)
    integer, intent(inout) :: A(*)
    integer :: PARTN, TEMP
    include "sort.f9h"
  end subroutine ISORT

  subroutine SSORT (A, M, N)
    real, intent(inout) :: A(*)
    real :: PARTN, TEMP
    include "sort.f9h"
  end subroutine SSORT

  subroutine DSORTP (A, M, N, P)
    double precision, intent(inout) :: A(*)
    integer, intent(in) :: M, N
    integer, intent(out) :: P(*)
    integer :: I
    do i = m, n
      p(i) = i
    end do
    call SORTQ ( A, M, N, P )
  end subroutine DSORTP

  subroutine ISORTP (A, M, N, P)
    integer, intent(inout) :: A(*)
    integer, intent(in) :: M, N
    integer, intent(out) :: P(*)
    integer :: I
    do i = m, n
      p(i) = i
    end do
    call SORTQ ( A, M, N, P )
  end subroutine ISORTP

  subroutine SSORTP (A, M, N, P)
    real, intent(inout) :: A(*)
    integer, intent(in) :: M, N
    integer, intent(out) :: P(*)
    integer :: I
    do i = m, n
      p(i) = i
    end do
    call SORTQ ( A, M, N, P )
  end subroutine SSORTP

  subroutine DSORTQ (A, M, N, P)
    double precision, intent(inout) :: A(*)
    double precision :: PARTN
    include "sortq.f9h"
  end subroutine DSORTQ

  subroutine ISORTQ (A, M, N, P)
    integer, intent(inout) :: A(*)
    integer :: PARTN
    include "sortq.f9h"
  end subroutine ISORTQ

  subroutine SSORTQ (A, M, N, P)
    real, intent(inout) :: A(*)
    real :: PARTN
    include "sortq.f9h"
  end subroutine SSORTQ

!---------------------------------------------------------------------------
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Sort_M

! $Log$
! Revision 2.3  2003/06/20 00:18:23  vsnyder
! Add SORTP and SORTQ
!
! Revision 2.2  2002/10/08 00:09:14  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.1  2002/01/30 19:36:34  vsnyder
! Initial commit
!
