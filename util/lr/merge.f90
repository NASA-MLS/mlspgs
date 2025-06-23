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

  subroutine MERGE ( IBASIS, IRES, CHANGE )

    use Basis_m, only: BASIS, INDBAS, Items
    use Delete, only: DELCS
    use LISTS, only: LINT, LIST
    use Tables, only: HEADEN, PRDIND => Prod_Ind, PRODCN => Productions
    use Union, only: CSUN

    implicit NONE

    ! Merge the configuration set at BASIS(IBASIS) into the set of
    ! configuration sets.  Set IRES to the resulting configuration set
    ! corresponding to IBASIS.  Set CHANGE false if a compatible set was
    ! found.  Set CHANGE true if no compatible sets were found.

    integer, intent(in) :: IBASIS
    integer, intent(out) :: IRES
    logical, intent(out) :: CHANGE

    ! *****     External References     ************************************

    ! CSUN    unions two context sets.
    ! DELCS   deletes a reference to a context set.
    ! LINT    inquires whether two context sets have an intersection.

    ! *****     Local Variables     ****************************************

    ! BPR     is the pointer to productions for the state at IBASIS.
    ! CHU     "CSUN(A,B) says A added something to B."
    ! I       is a loop induction variable and subscript.
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

    logical :: CHU
    integer BPR, I, IEND, IENT, II, IPTR, J, JJ, KK, LL, NIPTR, PPR
    logical SAME

    ! *****     Procedures     *********************************************

    ! Search the list of configuration sets having the same entrance
    ! symbol as IBASIS.

    bpr = basis(ibasis)%item     ! Context set
    ient = prodcn(prdind(items(bpr)%prod)+items(bpr)%dot-1)
    niptr = headen(ient)
  o:do
      iptr = niptr
      if ( iptr == 0 ) exit
      niptr = basis(iptr)%same ! Same entrance symbol
      ppr = basis(iptr)%item   ! Context set

      ! Does config IPTR have the same number of basis productions as
      ! config IBASIS?

      iend = basis(iptr+1)%item - ppr
      if ( iend == basis(ibasis+1)%item - bpr ) then

        ! Compare the basis configurations.

        same = .true.
        j = bpr
        do i = ppr, basis(iptr+1)%item - 1
          if ( items(i)%prod /= items(j)%prod ) cycle o
          if ( items(i)%dot /= items(j)%dot ) cycle o
          if ( items(i)%set /= items(j)%set ) same = .false.
          j = j + 1
        end do

        ! Are the config sets compatible?  They are not if there would be
        ! two intersecting context sets created by the merge where there
        ! was no intersection before.

        if ( .not. same ) then
          do i = 1, iend-1
            do j = i+1, iend
              ii = items(ppr+i-1)%set
              jj = items(bpr+j-1)%set
              kk = items(ppr+j-1)%set
              ll = items(bpr+i-1)%set
              if (  .not. lint(list(ii)%next,list(jj)%next) &
              .and. .not. lint(list(kk)%next,list(ll)%next) ) cycle
              if ( lint(list(ii)%next,list(kk)%next)  ) cycle
              if ( .not. lint(list(ll)%next,list(jj)%next) ) cycle o
            end do
          end do
        end if

        ! The configuration sets are compatible.  Merge them by unioning
        ! the context set lists.  If the context set lists are equal then
        ! no change will occur when they are merged.  Return the old basis
        ! and delete the trial basis.  Then delete the context sets in the
        ! trial basis

        change = .false.
        do i = 1, iend
          call csun ( items(bpr+i-1)%set, items(ppr+i-1)%set, chu )
          if ( chu ) change = .true.
          call delcs (items(bpr+i-1)%set)
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

    change = .true.
    ires = ibasis
    basis(ibasis)%same = headen(ient)
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
! Revision 1.1  2013/10/24 22:41:14  vsnyder
! Initial commit
!
