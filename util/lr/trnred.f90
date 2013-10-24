! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Transitions_And_Reductions

  implicit NONE
  private
  public :: TRNRED

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine TRNRED (IBASIS, JMAX)

    use ANACOM, only: NXTRED, NXTTRN
    use Error_Handler, only: Error
    use Merge_Sets, only: Merge
    use SCRCOM
    use S3, only: PRDIND, PRODCN
    use S5, only: ADDBAS, BASIS, ENQUE, NEWBAS, RED, TRAN
    implicit NONE

  ! Attach lists of transitions and reductions to the state in SCRTCH.

  ! IBASIS  is the point in the BASIS array corresponding to the state in
  !         SCRTCH.
  ! JMAX    is the largest position in SCRTCH used by the state.

    integer, intent(in) :: IBASIS, JMAX

  ! *****     External References     ********************************

  ! ADDBAS  adds an element to the BASIS array.
  ! ENDBAS  cleans up after adding a new element to the BASIS array.
  ! ENQUE   adds a basis to the queue for processing by ANALYZ.
  ! ERROR   prints error messages.
  ! NEWBAS  gets ready to add new elements to the BASIS array.

  ! *****     Local Variables     ************************************

  ! I       is a loop induction variable and subscript.
  ! ICH     indicates whether MERGE made a change.
  ! IPATH   indicates whether a path is found to a new configuration.
  ! IRED    is the pointer to the next place a reduction may be stored in
  !         RED.
  ! ITRAN   is the pointer to the next place a transition may be stored in
  !         TRAN.
  ! LHS     is the left side of a production.
  ! MAXR    indicates the amount of space remaining in RED.
  ! MAXT    indicates the amount of space remaining in TRAN.
  ! NB      is the position in BASIS of a new state added by MERGE.
  ! NBASIS  is the position in BASIS of a newly constructed basis, to be
  !         merged into or added onto the existing set of configurations
  !         by MERGE.

    integer I, ICH, IPATH, IRED, ITRAN, LHS, MAXR, MAXT, NB, NBASIS

  ! *****     Procedures     *****************************************

  ! Calculate how much space remains in TRAN.
  ! If space was previously allocated for the transitions from IBASIS
  ! then reuse that space.  The amount of space required will never
  ! change since the transitions depend only on the completed basis.
  ! The completed basis is never changed, but its contexts are.

    maxt = basis(ibasis+8)
    itran = basis(ibasis+3)

    i = 1
    if (scrtch(2) < prdind(scrtch(1)+1) - prdind(scrtch(1))) then
      ipath = 1
      do while (ipath /= 0)
        lhs = prodcn(prdind(scrtch(i))+scrtch(i+1))
        call newbas (nbasis)
        do while (lhs == prodcn(prdind(scrtch(i))+scrtch(i+1)))
          call addbas (nbasis,scrtch(i),scrtch(i+1)+1,scrtch(i+2))
          i = i + 3
          ipath = 0
          if (i > jmax) go to 10
          if (scrtch(i+1) >= prdind(scrtch(i)+1) - prdind(scrtch(i))) go to 10
        end do
        ipath = 1
  10    continue
        call merge (nbasis, nb, ich)
        if (ich /= 0) call enque (nb)
        ! Add a transition to NB to the basis at IBASIS.
        if (itran >= maxt) call error ('Transition area (TRAN) overflow',2)
        tran(itran) = nb
        itran = itran + 1
      end do
    end if
    nxttrn = max(nxttrn, itran)
    basis(ibasis+8) = itran

    ! Calculate how much space remains in RED.
    ! If space was previously allocated for the reductions from IBASIS
    ! then reuse that space.  The amount of space required will never
    ! change since the reductions depend only on the completed basis.
    ! The completed basis is never changed, but its contexts are.

    maxr = basis(ibasis+9)
    ired = basis(ibasis+4)

    ! Construct or reconstruct reductions.

    do while (i <= jmax)
    ! Add the reduction of production SCRTCH(I) for context
    ! SCRTCH(I+2) to the state at IBASIS.
      if (ired >= maxr) call error ('Reduction area (RED) overflow',2)
      red(ired) = scrtch(i)
      red(ired+1) = scrtch(i+2)
      ired = ired + 2
      i = i + 3
    end do
    nxtred = max(nxtred, ired)
    basis(ibasis+9) = ired

  end subroutine TRNRED

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Transitions_And_Reductions

! $Log$
