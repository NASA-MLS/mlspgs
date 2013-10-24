! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Merge_Sets

  implicit NONE
  private
  public :: MERGE

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine MERGE (IBASIS, IRES, ICHNG)

    use ANACOM, only: INDBAS
    use Delete, only: DELCS
    use LISTS, only: LINT, NEXT
    use S3, only: HEADEN, PRDIND, PRODCN
    use S5, only: BASIS
    use Union, only: CSUN

    implicit NONE

    ! Merge the configuration set at BASIS(IBASIS) into the set of
    ! configuration sets.  Set IRES to the resulting configuration set
    ! corresponding to IBASIS.  Set ICHNG zero if a compatible set was
    ! found.  Set ICHNG 1 if no compatible sets were found.

    integer, intent(in) :: IBASIS
    integer, intent(out) :: IRES, ICHNG

    ! *****     External References     ************************************

    ! CSUN    unions two context sets.
    ! DELCS   deletes a reference to a context set.
    ! LINT    inquires whether two context sets have an intersection.

    ! *****     Local Variables     ****************************************

    ! BPR     is the pointer to productions for the state at IBASIS.
    ! I       is a loop induction variable and subscript.
    ! ICH     is the value of ICHNGE reported by CSUN.
    ! IEND    is the upper limit for I.
    ! IENT    is the entry symbol for IBASIS.
    ! II      is a temporary variable.
    ! IPTR    points to a basis being tested for compatibility with IBASIS.
    ! J       is a loop induction variable and subscript.
    ! JJ      is a temporary variable.
    ! KK      is a temporary variable.
    ! LL      is a temporary variable.
    ! NIPTR   is the next value for IPTR
    ! PPR     is the pointer to productions for the state at IPTR.
    ! SAME    indicates two configuration sets are identical.

    integer BPR, I, ICH, IEND, IENT, II, IPTR, J, JJ, KK, LL, NIPTR, PPR
    logical SAME

    ! *****     Procedures     *********************************************

    ! Search the list of configuration sets having the same entrance
    ! symbol as IBASIS.

    bpr = basis(ibasis)
    ient = prodcn(prdind(basis(bpr))+basis(bpr+1)-1)
    niptr = headen(ient)
  o:do
      iptr = niptr
      if ( iptr == 0 ) exit
      niptr = basis(iptr+2)
      ppr = basis(iptr)

      ! Does config IPTR have the same number of basis productions as
      ! config IBASIS?

      iend = ppr - basis(iptr+5)
      if (iend == bpr - basis(ibasis+5)) then

        ! Compare the basis configurations.

        same = .true.
        j = bpr
        do i = ppr, basis(iptr+5)+3, -3
          if (basis(i) /= basis(j)) cycle o
          if (basis(i+1) /= basis(j+1)) cycle o
          if (basis(i+2) /= basis(j+2)) same = .false.
          j = j - 3
        end do

        ! Are the config sets compatible?  They are not if there would be
        ! two intersecting context sets created by the merge where there
        ! was no intersection before.

        if (.not. same) then
          do i = 1, iend-3, 3
            do j = i+3, iend, 3
              ii = basis(ppr-i+3)
              jj = basis(bpr-j+3)
              kk = basis(ppr-j+3)
              ll = basis(bpr-i+3)
              if (  .not. lint(next(ii),next(jj)) &
              .and. .not. lint(next(kk),next(ll)) ) cycle
              if ( lint(next(ii),next(kk))  ) cycle
              if ( .not. lint(next(ll),next(jj)) ) cycle o
            end do
          end do
        end if

        ! The configuration sets are compatible.  Merge them by unioning
        ! the context set lists.  If the context set lists are equal then
        ! no change will occur when they are merged.  Return the old basis
        ! and delete the trial basis.  Then delete the context sets in the
        ! trial basis

        ichng = 0
        do i = 1, iend, 3
          call csun (basis(bpr-i+3), basis(ppr-i+3), ich)
          if (ich /= 0) ichng = 1
          call delcs (basis(bpr-i+3))
        end do

        ! The trial basis is always constructed as the last basis in the
        ! BASIS array.  Since we have merged it into another basis we can
        ! now delete it.  The context sets for the trial basis have already
        ! been released.  All that remains is to reset the pointer into
        ! the BASIS array.

        indbas = ibasis
        ires = iptr
        return
      end if

      ! Try the next config set.

    end do o

    ! The trial set is not compatible with any existing config set.
    ! Link the trial set into the entrance symbol chain and return it
    ! as the resulting set.

    ichng = 1
    ires = ibasis
    basis(ibasis+2) = headen(ient)
    headen(ient) = ibasis

  end subroutine MERGE

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Merge_Sets

! $Log$
