! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Complete

  implicit NONE
  private
  public :: COMPLT

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine COMPLT (ISTATE, MAXSET)

    use Delete, only: DELCS
    use Immediate_Transitions, only: IMTRCS
    use LISTS, only: ITEM
    use SCRCOM
    use S3, only: FRSPRD, PRDIND, PRODCN
    use S5, only: BASIS
    use TABCOM, only: NTERMS, NUMPRD
    use Union, only: CSUN

    implicit NONE

    ! Complete the state ISTATE.

    integer ISTATE, MAXSET

    ! *****     External References     ************************************

    ! CSUN    constructs the union of two context sets.
    ! DELCS   deletes a context set.
    ! ERROR   prints error messages.
    ! IMTRCS  constructs the immediate transition context set.

    ! *****     Local Variables     ****************************************

    ! I       is a loop inductor and subscript for BASIS.
    ! ICH     is a value of ICHANG returned by CSUN.
    ! ICHANG  indicates whether a change to the configuration occurred
    !         during a loop iteration.
    ! IEND    is the upper limit for I.
    ! IPTR    is a pointer to a context set list.
    ! ISTART  is the lower limit for I.
    ! J       is a loop inductor and subscript for SCRTCH.
    ! K       is a loop inductor and subscript for MARK.
    ! LHS     is the next symbol scanned in a production, which will be the
    !         left hand side of the next production added to the closure if
    !         it is a nonterminal symbol.
    ! MARK    is an array used to indicate whether a production has already
    !         been added to the closure.

    integer I,ICH,ICHANG,IEND,IPTR,ISTART,J,K,LHS
    integer MARK(NUMPRD)

    !     *****     Procedures     *****************************************

    ! First move the basis to scratch.

    istart = basis(istate)
    iend = basis(istate+5) + 3
    if (iend - istart + 1 > scrsiz) call scratch_overflow
    j = 0
    do i = istart, iend, -3
      scrtch(j+1) = basis(i)
      scrtch(j+2) = basis(i+1)
      iptr = basis(i+2)
      scrtch(j+3) = iptr
      item(iptr) = item(iptr) + 1
      j = j + 3
    end do

    ! Unmark all productions.

    mark(1:numprd) = 0

    ! Complete the configuration set by adding all the immediate
    ! transitions to scratch.

    do ! until (ichang == 0)
      ichang = 0
      i = 1
      do ! until (i >= j)
      if (scrtch(i+1) < prdind(scrtch(i)+1)-prdind(scrtch(i))) then
         lhs = prodcn(prdind(scrtch(i)) + scrtch(i+1))

         ! If the dot is before a nonterminal symbol then there is an
         ! immediate transition from the production.

           if (lhs > nterms) then

           ! Generate the context set which is the same for all the
           ! immediate transitions from this production.

             call imtrcs (i, iptr)

             ! Add all the un-marked productions with the left hand side
             ! equal to the nonterminal LHS to the right of the dot.
             ! Union in the new context set if the production has already
             ! been included.

             k = frsprd(lhs)
             if ( k > 0 ) then   ! Otherwise it an unused symbol defined by &V
               do while (prodcn(prdind(k)) == lhs)
                 if (mark(k).eq.0) then
                   if (j+3 .gt. scrsiz) call scratch_overflow
                   mark(k) = j + 3
                   scrtch(j+1) = k
                   scrtch(j+2) = 1
                   scrtch(j+3) = iptr
                   item(iptr) = item(iptr) + 1
                   j = j + 3
                   ichang = 1
                 else
                   call csun (iptr, scrtch(mark(k)), ich)
                   if (ich /= 0) ichang = 1
                 end if
                 k = k + 1
               end do
             end if

             ! The call to delete the context set deletes the "extra"
             ! reference to the set and deletes the set completely from the
             ! list space if it was never referenced.  As a result of the
             ! call to IMTRCS, IPTR had its ref count incremented.

             call delcs (iptr)
           end if
         end if
         i = i + 3
         if (i >= j) exit
       end do
       if (ichang == 0) exit
    end do
    maxset = j

  contains

    subroutine scratch_overflow
      use Error_Handler, only: Error
      call error ('Configuration set too large',2)
    end subroutine scratch_overflow

  end subroutine COMPLT

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Complete

! $Log$
