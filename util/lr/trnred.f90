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

  ! Red_t is the type of object used to keep track of reductions
  type :: Red_t
    integer :: Prod   ! Production number
    integer :: Set    ! Context set index
  end type Red_t

  ! Red(state) keeps track of reductions for "state"
  type(red_t), allocatable, save, public :: Red(:)
  integer, parameter :: Red_Init = 1000
  integer, save :: Red_Size = 0 

  ! Tran(state) keeps track of transitions for "state".  The values are
  ! other state numbers.
  integer, allocatable, save, public :: Tran(:)
  integer, parameter :: Tran_Init = 1000
  integer, save :: Tran_Size = 0 

  integer, public :: NXTRED = 1, NXTTRN = 1

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine TRNRED ( IBASIS, JMAX )

    use Basis_m, only: ADDBAS, BASIS, ENQUE, NEWBAS
    use Complete, only: Scrtch
    use Tables, only: PRDIND => Prod_Ind, PRODCN => Productions
    use Merge_Sets, only: Merge
    use Toggles, only: Gen, Levels
    use Trace, only: Trace_Begin, Trace_End

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
  ! NEWBAS  gets ready to add new elements to the BASIS array.

  ! *****     Local Variables     ************************************

  ! CHANGE  "MERGE made a change."
  ! I       is a loop induction variable and subscript.
  ! IPATH   indicates whether a path is found to a new configuration.
  ! IRED    is the pointer to the next place a reduction may be stored in
  !         RED.
  ! ITRAN   is the pointer to the next place a transition may be stored in
  !         TRAN.
  ! LHS     is the left side of a production.
  ! NB      is the position in BASIS of a new state added by MERGE.
  ! NBASIS  is the position in BASIS of a newly constructed basis, to be
  !         merged into or added onto the existing set of configurations
  !         by MERGE.
  ! Temp_Red  is used to allocate a new Red array
  ! Temp_Tran is used to allocater a new Tran array

    logical :: CHANGE
    integer I, IPATH, IRED, ITRAN, LHS, NB, NBASIS
    type(red_t), allocatable :: Temp_Red(:)
    integer, allocatable :: Temp_Tran(:)

  ! *****     Procedures     *****************************************

  ! Calculate how much space remains in TRAN.
  ! If space was previously allocated for the transitions from IBASIS
  ! then reuse that space.  The amount of space required will never
  ! change since the transitions depend only on the completed basis.
  ! The completed basis is never changed, but its contexts are.

    if ( levels(gen) > 1 ) call trace_begin ( 'TRNRED' )

    itran = basis(ibasis)%tran   ! Start of transitions from this basis

    i = 1
    if (scrtch(1)%dot < prdind(scrtch(1)%prod+1) - prdind(scrtch(1)%prod)) then
      ipath = 1
      do while (ipath /= 0)
        lhs = prodcn(prdind(scrtch(i)%prod)+scrtch(i)%dot)
        call newbas ( nbasis )
        do while ( lhs == prodcn(prdind(scrtch(i)%prod)+scrtch(i)%dot) )
          call addbas ( nbasis, scrtch(i)%prod, scrtch(i)%dot+1, scrtch(i)%set )
          i = i + 1
          ipath = 0
          if ( i > jmax ) go to 10
          if ( scrtch(i)%dot >= &
               prdind(scrtch(i)%prod+1) - prdind(scrtch(i)%prod) ) go to 10
        end do
        ipath = 1
  10    continue
        call merge ( nbasis, nb, change )
        if ( change ) call enque (nb)
        ! Add a transition to NB to the basis at IBASIS.
        if ( itran > tran_size ) then
          if ( .not. allocated(tran) ) then
            allocate ( tran(tran_init) )
          else
            allocate ( temp_tran(2*tran_size) )
            temp_tran(:tran_size) = tran
            call move_alloc ( temp_tran, tran )
          end if
          tran_size = size(tran)
        end if
        tran(itran) = nb
        itran = itran + 1
      end do
    end if
    nxttrn = max(nxttrn, itran)
    basis(ibasis+1)%tran = itran

    ! Calculate how much space remains in RED.
    ! If space was previously allocated for the reductions from IBASIS
    ! then reuse that space.  The amount of space required will never
    ! change since the reductions depend only on the completed basis.
    ! The completed basis is never changed, but its contexts are.

    ired = basis(ibasis)%red   ! Start of reductions for this basis

    ! Construct or reconstruct reductions.

    do while (i <= jmax)
    ! Add the reduction of production SCRTCH(I)%prod for context
    ! SCRTCH(I)%set to the state at IBASIS.
      if ( ired >= red_size ) then
        if ( .not. allocated(red) ) then
          allocate ( red(red_init) )
        else
          allocate ( temp_red(2*red_size) )
          temp_red(:red_size) = red
          call move_alloc ( temp_red, red )
        end if
        red_size = size(red)
      end if
      red(ired)%prod = scrtch(i)%prod
      red(ired)%set = scrtch(i)%set
      ired = ired + 1
      i = i + 1
    end do
    nxtred = max(nxtred, ired)
    basis(ibasis+1)%red = ired

    if ( levels(gen) > 1 ) call trace_end ( 'TRNRED' )

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
! Revision 1.1  2013/10/24 22:41:14  vsnyder
! Initial commit
!
