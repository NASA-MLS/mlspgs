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

  use Basis_m, only: Item_t

  implicit NONE
  private
  public :: Complt

  integer, public, allocatable, save :: Closures_Index(:) ! 0 : Number of states
    ! Closures_Index(i) is index of last closure for state I in Closures.

  type(item_t), public, allocatable, save :: Closures(:)  ! Closure for state I
    ! is in Closures ( closures_index(i-1)+1 : closures_index(i) )

  type(item_t), private, allocatable, save :: SCRTCH(:)

  integer, parameter, private :: Init_Closures_Index = 1000 ! Initial size

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Complt ( ISTATE )

    use Basis_m, only: Basis, Increase_Items, Item_T, Items
    use Delete, only: DELCS
    use LISTS, only: LIST
    use Tables, only: FRSPRD => First_Production, NTERMS, NUMPRD, &
      & PRDIND => Prod_Ind, PRODCN => Productions
    use Union, only: CSUN

    implicit NONE

    ! Complete the state ISTATE.

    integer, intent(in) :: ISTATE

    ! *****     External References     ************************************

    ! CSUN    constructs the union of two context sets.
    ! DELCS   deletes a context set.
    ! ERROR   prints error messages.
    ! IMTRCS  constructs the immediate transition context set.

    ! *****     Local Variables     ****************************************

    ! CHANGE  indicates whether a change to the configuration occurred
    ! CHU     is a value of CHANGE returned by CSUN.
    ! CLOSIZ  is the size of the Closures array
    ! I       is a loop inductor and subscript for Basis.
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
    ! SCRSIZ  size of scratch array.

    logical :: Change, CHU
    integer :: CLOSIZ = 0  ! Closures array starts out not allocated
    integer I, IEND, IPTR, ISTART, J, K, LHS
    integer MARK(NUMPRD)
    integer :: SCRSIZ = 0  ! SRCTCH starts out not allocated

    !     *****     Procedures     *****************************************

    ! Make room for more closures if necessary

    if ( .not. allocated(closures_index) ) then
      call increase_closures_index
    else if ( istate > ubound(closures_index,1) ) then
      call increase_closures_index
    end if

    ! First move the set's configuration items to scratch.

    istart = basis(istate)%item      ! First item of the configuration set
    iend = basis(istate+1)%item - 1  ! Last item of the configuration set
    if ( iend - istart + 1 > scrsiz ) call increase_items ( scrtch, scrsiz )
    j = 0
    do i = istart, iend
      j = j + 1
      scrtch(j) = items(i)
      iptr = items(i)%set
      list(iptr)%item = list(iptr)%item + 1
    end do

    ! Unmark all productions.

    mark(1:numprd) = 0

    ! Complete the configuration set by adding all the immediate
    ! transitions to scratch.

    change = .true.
    do while ( change )
      change = .false.
      i = 1
      do ! until ( i > j )
         if ( scrtch(i)%dot < prdind(scrtch(i)%prod+1)-prdind(scrtch(i)%prod)) then
           lhs = prodcn(prdind(scrtch(i)%prod) + scrtch(i)%dot)

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
             if ( k > 0 ) then   ! Otherwise it is an unused symbol defined by
                                 ! symbol = name in the grammar
               do while (prodcn(prdind(k)) == lhs)
                 if ( mark(k) == 0 ) then
                   if ( j >= scrsiz ) call increase_items ( scrtch, scrsiz )
                   j = j + 1
                   mark(k) = j
                                    ! prod dot set
                   scrtch(j) = item_t(k, 1, iptr)
                   list(iptr)%item = list(iptr)%item + 1
                   change = .true.
                 else
                   call csun ( iptr, scrtch(mark(k))%set, chu )
                   if ( chu ) change = .true.
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
         i = i + 1
         if ( i > j ) exit
       end do ! until ( i > j )
    end do ! while ( change )

    ! Copy items from scrtch to closures
    closures_index(istate) = closures_index(istate-1) + j
    if ( closures_index(istate) > closiz ) &
      & call increase_items ( closures, closiz )
    closures(closures_index(istate-1)+1:closures_index(istate)) = &
      & scrtch(1:j)

  end subroutine Complt

  ! IMTRCS is private, called only by Complt.

  subroutine IMTRCS (Ibasis, IPTR)

    use FIRST_SETS, only: FIRST_PT        ! heads of first sets
    use LISTS, only: ADDLTL, COPYL, LIST
    use New_Context_Set, only: NEWCS
    use NULLABLE_M, only: NULLABLE        ! nullable nonterminals
    use Tables, only: PRDIND => Prod_Ind, PRODCN => Productions

  ! Construct the immediate transition context set for the
  ! configuration at SCRTCH(Ibasis).  Return the pointer to it in
  ! IPTR.

  ! The immediate transition context set is the set consisting of the
  ! THEADS of the substring of symbols in the production following the
  ! symbol after the dot, unioned with the context set of the
  ! production if that string is null or potentially null.

    integer Ibasis, IPTR

  ! *****     External References     ********************************

  ! ADDLTL  adds a list to a list (Set Union operation).
  ! COPYL   copies a list.
  ! NEWCS   prepares a new context set.

  ! *****     Local Variables     ************************************

  ! CHANGE  "ADDLTL changed a list."
  ! I       is a loop induction variable and subscript.
  ! IEND    is the upper limit for I.
  ! IP      is the pointer to the list being constructed.
  ! ISTART  is the lower limit for I.

    logical :: CHANGE
    integer :: I, IEND, IP, ISTART

  ! *****     Procedures     *****************************************

    i = scrtch(ibasis)%prod

  ! If there is no substring just return the context set for the
  ! production.

    if (scrtch(ibasis)%dot >= prdind(i+1)-prdind(i)-1) then
      iptr = scrtch(ibasis)%set
      list(iptr)%item = list(iptr)%item + 1
      return
    end if

  ! The set consists of the THEAD set of the symbol unioned with the
  ! THEAD set of the following substring if the symbol is nullable.
  ! If the rest of the RHS is nullable, the set includes the context
  ! set of the production.

    istart = scrtch(ibasis)%dot + prdind(i) + 1
    iend = prdind(i+1) - 1
    call copyl (first_pt(prodcn(istart)), ip)
    if ( nullable(prodcn(istart)) ) then
      do i = istart+1, iend
        call addltl ( first_pt(prodcn(i)), ip, change )
        if ( .not. nullable(prodcn(i)) ) go to 9
      end do

      ! Add the context set of the production.

      call addltl ( list(scrtch(ibasis)%set)%next, ip, change )
    end if

  ! Create a new reference to the context set with IP elements.

  9 continue
    call newcs ( IP, IPTR )

  end subroutine IMTRCS

  subroutine Increase_Closures_Index

    ! If Closures_Index is not allocated, allocate it with bounds
    ! 0 : Init_Closures_Index.

    ! If it is allocated, double its size.

    integer :: N ! Upper bound of Closures_Index if it's allocated
    integer, allocatable :: Temp(:)

    if ( .not. allocated(closures_index) ) then
      allocate ( closures_index(0:init_closures_index) )
      closures_index(0) = 0
      return
    end if

    n = ubound(closures_index,1)
    allocate ( temp(0:n) )
    temp(0:n) = closures_index
    call move_alloc ( temp, closures_index )

  end subroutine Increase_Closures_Index

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
! Revision 1.2  2014/01/14 00:11:42  vsnyder
! Revised LR completely
!
! Revision 1.1  2013/10/24 22:41:13  vsnyder
! Initial commit
!
