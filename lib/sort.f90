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
    module procedure ASORT, DSORT, ISORT, SSORT
  end interface

  interface SORTP ! Using a permutation vector
    module procedure ASORTP, DSORTP, ISORTP, SSORTP
  end interface

  interface SORTQ ! Using a pre-specified permutation vector
    module procedure ASORTQ, DSORTQ, ISORTQ, SSORTQ
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

  subroutine ASORT (A, M, N)
    character(len=*), intent(inout) :: A(*)
    character(len=len(a)) :: PARTN, TEMP
    include "sort.f9h"
  end subroutine ASORT

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

  subroutine ASORTP (A, M, N, P)
    character(len=*), intent(in) :: A(*)
    integer, intent(in) :: M, N
    integer, intent(out) :: P(*)
    integer :: I
    do i = m, n
      p(i) = i
    end do
    call SORTQ ( A, M, N, P )
  end subroutine ASORTP

  subroutine DSORTP (A, M, N, P)
    double precision, intent(in) :: A(*)
    integer, intent(in) :: M, N
    integer, intent(out) :: P(*)
    integer :: I
    do i = m, n
      p(i) = i
    end do
    call SORTQ ( A, M, N, P )
  end subroutine DSORTP

  subroutine ISORTP (A, M, N, P)
    integer, intent(in) :: A(*)
    integer, intent(in) :: M, N
    integer, intent(out) :: P(*)
    integer :: I
    do i = m, n
      p(i) = i
    end do
    call SORTQ ( A, M, N, P )
  end subroutine ISORTP

  subroutine SSORTP (A, M, N, P)
    real, intent(in) :: A(*)
    integer, intent(in) :: M, N
    integer, intent(out) :: P(*)
    integer :: I
    do i = m, n
      p(i) = i
    end do
    call SORTQ ( A, M, N, P )
  end subroutine SSORTP

  subroutine ASORTQ (A, M, N, P)
    character(len=*), intent(in) :: A(*)
    character(len=len(a)) :: PARTN
    include "sortq.f9h"
  end subroutine ASORTQ

  subroutine DSORTQ (A, M, N, P)
    double precision, intent(in) :: A(*)
    double precision :: PARTN
    include "sortq.f9h"
  end subroutine DSORTQ

  subroutine ISORTQ (A, M, N, P)
    integer, intent(in) :: A(*)
    integer :: PARTN
    include "sortq.f9h"
  end subroutine ISORTQ

  subroutine SSORTQ (A, M, N, P)
    real, intent(in) :: A(*)
    real :: PARTN
    include "sortq.f9h"
  end subroutine SSORTQ

      subroutine GSORTP (COMPAR, N, P)
!     .  Copyright (C) 1989, California Institute of Technology.
!     .  All rights reserved.  U. S. Government sponsorship under
!     .  NASA contract NAS7-918 is acknowledged.
!>> 1998-01-20  GSORTP  Snyder  Allow not initializing P.
!>> 1996-05-01  GSORTP  Krogh   Changes to use .C. and C%%.
!>> 1995-11-17  GSORTP  Krogh   Converted SFTRAN to Fortran 77.
!>> 1991-04-02  GSORTP  Snyder  Repair no permutation vector if m-n < 10
!>> 1988-11-22  GSORTP  Snyder  Initial code.
!
!     Sort an N-vector of objects of unknown type and organization.
!     P is set so that the P(J)'th element of the original sequence is
!     the J'th element of the sorted sequence.  The order is defined by
!     the integer function COMPAR.  An invocation COMPAR(I,J) should
!     return -1 if the I'th element of the data is to preceed the J'th
!     element in the sorted sequence, +1 if the J'th element is to
!     preceed the I'th element in the sorted sequence, and 0 if the I'th
!     and J'th elements are the same.
!
!     This subprogram is unaware of the data, and cannot manipulate it.
!     It is the caller's responsibility to make the data known to the
!     COMPAR function.
!
      integer COMPAR, N, P(*)
      external COMPAR
!
!     *****     Local Variables     ************************************
!
! BL      Left bound of the sub-array to be sorted at the next step.
! BR      Right bound of the sub array to be sorted at the next step.
! CL      Current left bound of the unsorted sub-array.
! CR      Current right bound of the unsorted sub-array.
! MYN     My N.
! PARTN   is the partition element.
! PTEMP   holds elements of P during exchanges.
! STACKL  keeps track of the left bounds of sub-arrays waiting to be
!         sorted.
! STACKR  keeps track of the right bounds of sub-arrays waiting to be
!         sorted.
! STKTOP  keeps track of the top of the stacks.
!
      integer BL,BR,CL,CR,MYN,PARTN,PTEMP
      integer STACKL(32),STACKR(32),STKTOP
!
!     *****     Executable Statements     ******************************
!
      do 20 cl = 1, n
         p(cl)=cl
   20 continue
      myn = abs(n)
      if (myn.ge.10) then
         stktop=1
         stackl(1)=1
         stackr(1)=myn
!           Start until loop
   40    continue
            bl=stackl(stktop)
            br=stackr(stktop)
            stktop=stktop-1
!           Choose a partitioning element.  Use the median of the first,
!           middle and last elements.  Sort them so the extreme elements
!           can serve as sentinels during partitioning.
            cl=(bl+br)/2
            partn=p(cl)
!%%         if ((*compar)( P[bl], partn ) > 0 ){
            if (compar(p(bl),partn).gt.0) then
               p(cl)=p(bl)
               p(bl)=partn
               partn=p(cl)
!%%         }
            end if
!%%         if ((*compar)( P[bl], P[br] ) > 0 ){
            if (compar(p(bl),p(br)).gt.0) then
               ptemp=p(bl)
               p(bl)=p(br)
               p(br)=ptemp
!%%         }
            end if
!%%         if ((*compar)( partn, P[br] ) > 0 ){
            if (compar(partn,p(br)).gt.0) then
               p(cl)=p(br)
               p(br)=partn
               partn=p(cl)
!%%         }
            end if
            p(cl)=p(br-1)
            p(br-1)=partn
!           Partition the sub-array around PARTN.  Exclude the above
!           considered elements from partitioning because they're al-
!           ready in the correct subfiles.  Stop scanning on equality to
!           prevent files containing equal values from causing a loop.
            cl=bl
            cr=br-1
!              Start forever block
   60       continue
   80          cl=cl+1
!%%               if ((*compar)( P[cl], partn ) < 0) goto L_80;
                  if (compar(p(cl),partn) .lt. 0) go to 80
  100          cr=cr-1
!%%               if ((*compar)( P[cr], partn ) > 0) goto L_100;
                  if (compar(p(cr),partn) .gt. 0) go to 100
               if (cl.gt.cr) go to 120
               ptemp=p(cl)
               p(cl)=p(cr)
               p(cr)=ptemp
               go to 60
!              End forever block
  120       continue
!           Put sub-arrays on the stack if they're big enough.  Put the
!           larger under the smaller, so the smaller will be done next.
!           This makes the upper bound of the stack depth log2 (myn).
!           (The "Hibbard" modification of quicksort).
            if (cl-bl .gt. br-cr) then
               if (cl-bl.gt.10) then
                  stktop=stktop+1
                  stackl(stktop)=bl
                  stackr(stktop)=cr
               end if
               if (br-cr.gt.10) then
                  stktop=stktop+1
                  stackl(stktop)=cl
                  stackr(stktop)=br
               end if
            else
               if (br-cr.gt.10) then
                  stktop=stktop+1
                  stackl(stktop)=cl
                  stackr(stktop)=br
               end if
               if (cl-bl.gt.10) then
                  stktop=stktop+1
                  stackl(stktop)=bl
                  stackr(stktop)=cr
               end if
            end if
!           End until loop
         if (stktop .ne. 0) go to 40
      end if
!     Clean up small subfiles by using insertion sort on everything.
      do 160 cr = 2, myn
         ptemp=p(cr)
         cl=cr
  140    continue
!%%      if ((*compar)( P[cl - 1], ptemp ) > 0) {
         if (compar(p(cl-1),ptemp).gt.0) then
            p(cl)=p(cl-1)
            cl=cl-1
            if (cl .gt. 1) go to 140
!%%      }
         end if
         p(cl)=ptemp
  160 continue
      return
      end subroutine GSORTP

!---------------------------------------------------------------------------
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Sort_M

! $Log$
! Revision 2.6  2005/01/12 03:01:59  vsnyder
! Added character sorts
!
! Revision 2.5  2004/10/01 23:09:18  vsnyder
! Add GSORTP
!
! Revision 2.4  2003/06/20 01:34:30  vsnyder
! Change intent for A and P
!
! Revision 2.3  2003/06/20 00:18:23  vsnyder
! Add SORTP and SORTQ
!
! Revision 2.2  2002/10/08 00:09:14  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.1  2002/01/30 19:36:34  vsnyder
! Initial commit
!
