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

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

  subroutine DSORT (A, M, N)
    double precision, intent(inout) :: A(*)
    double precision :: PARTN,TEMP
    include "sort.f9h"
  end subroutine DSORT

  subroutine ISORT (A, M, N)
    integer, intent(inout) :: A(*)
    integer :: PARTN,TEMP
    include "sort.f9h"
  end subroutine ISORT

  subroutine SSORT (A, M, N)
    real, intent(inout) :: A(*)
    real :: PARTN,TEMP
    include "sort.f9h"
  end subroutine SSORT
end module Sort_M

! $Log$
! Revision 2.1  2002/01/30 19:36:34  vsnyder
! Initial commit
!
