module HASH
  implicit NONE

  private
  public :: LOOKUP_AND_INSERT
  public :: DELETE
  public :: KEY                      ! Get the original key of an entry

  ! Parameters for value of STATUS returned by LOOKUP_AND_INSERT or DELETE
  ! (q.v.):
  integer, public, parameter :: FOUND = 1    ! Found or deleted G at LOC
  integer, public, parameter :: INSERTED = 2 ! Inserted G at LOC
  integer, public, parameter :: FULL = 3     ! Table is full
  integer, public, parameter :: NOT_KEY = 4  ! Key of H(1,LOC) /= key of G
  integer, public, parameter :: BAD_LOC = 5  ! LOC < 2 or LOC > SIZE(H,2)
  integer, public, parameter :: EMPTY = 6    ! Trying to delete empty cell

  ! Variables used for performance measurement:
  integer, public :: NLINK   ! number of links examined to find the
    ! predecessors of entries that must be moved.
  integer, public :: NMOVE   ! number of times entries are moved.
  integer, public :: NSPROB  ! number of probes (key comparisons) to find
    ! entries that are present.
  integer, public :: NUPROB  ! number of probes (key comparisons) to
    ! discover entries are not present.

  integer, private, parameter :: KZERO=262143 ! the key value used if G is
    ! zero in LOOKUP_AND_INSERT.  A reasonable value is something close to
    ! the square root of the largest integer.

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

! =====================================================     DELETE     =====

  subroutine DELETE ( H, LOC, STATUS )
  ! Delete the object at LOC from the hash table H.

    integer, intent(inout) :: H(:,:) ! See LOOKUP_AND_INSERT
    integer, intent(in) :: LOC       ! Location of object to delete.  Should
    ! be the result of the most recent call to LOOKUP_AND_INSERT -- remember,
    ! LOOKUP_AND_INSERT can move things around.
    integer, intent(out) :: STATUS   ! Result of operation.  Either FOUND
    ! or BAD_LOC or EMPTY (see parameters above).

    ! *****     Local Entities     *************************************

    integer :: NEXTP                 ! Next value to use for P, q.v.
    integer :: NT                    ! Number of entries, SIZE(H,2) - 1
    integer :: P                     ! Pointer to an element of H.

    ! *****     Procedures     *****************************************

    nt = size(h,2) - 1
    if ( loc < 2 .or. loc > nt+1 ) then
      status = bad_loc
      return
    end if
    p = loc - 2
    if ( h(1,p+2) == 0 ) then
      status = empty
      return
    end if
    status = found

!   Delete the object by moving its successor into its place and
!   making the successor's place empty.  This must be repeated for
!   the successor if the successor was a head of chain.  If the chain
!   is of length 2 and we delete the one that's not the head, we must
!   point the head at itself.

!   The following loop is executed either once or twice.  It is
!   executed twice if the predecessor of the head of a chain longer
!   than one entry is deleted.

    do
      nextp=mod( abs(h(1,p+2)), nt )
      if ( nextp == p ) exit ! chain of length 1
      if ( nextp == loc - 2 ) then ! deleting non-head of chain of length 2
        h(1,p+2) = h(1,p+2) - nextp + p
        exit
       end if
    ! Move object at NEXTP to P.
      h(1,p+2) = sign( h(1,nextp+2), h(1,p+2) )
      h(2:,p+2) = h(2:,nextp+2)
      p = nextp
      if ( h(1,nextp+2) < 0 ) exit
    end do
  ! Make the entry at NEXTP empty.
    h(1,nextp+2)=0

    return
  end subroutine DELETE

! ========================================================     KEY     =====

  integer function KEY ( H, LOC )
  ! Return ABS(key value used to insert in hash table).  < 0 means
  ! "LOC refers to an empty cell."  A result of zero is ambiguous.  It
  ! could result from a key initially zero, or a key initially equal to the
  ! value, KZERO, used to represent a zero key.

    integer, intent(in) :: H(:,:)    ! Hash table.  See LOOKUP_AND_INSERT.
    integer, intent(in) :: LOC       ! Location of object in H, usually output
      ! from LOOKUP_AND_INSERT.

    integer :: NT   ! Number of entries possible in H = SIZE(H,2) - 1.
    integer :: P    ! Index to an element of H.  Used during searching in H.
    integer :: PN   ! Next P
    integer :: W    ! ABS(H(1,P)).
    integer :: WL   ! ABS(H(1,LOC)).

    if ( h(1,loc) == 0 ) then
      key = -1
      return
    end if
    nt = size( H, 2 ) - 1
    wl = abs( h(1,loc) )
    p = loc - 2
    do
      w = abs( h(1, p+2) )
      pn = mod(w, nt)
      if ( h(1,p+2) > 0 ) then
        key = w - pn + p
        if ( key == kzero ) key = 0
        return
      end if
      p = pn
    end do
  end function KEY

! ==========================================     LOOKUP_AND_INSERT     =====

  subroutine LOOKUP_AND_INSERT ( G, H, INSERT, LOC, STATUS )

  ! Lookup or lookup and insert into a hash table that uses separate chains.

    integer, intent(in) :: G         ! Element to be found or inserted.
    integer, intent(inout) :: H(:,:) ! Table in which G is to be found,
      ! or into which G is to be inserted.  Column 1 of H is used
      ! internally for control information.  G is stored in row 1 of H, in
      ! some column between 1 and NT inclusive.  The initial values of
      ! H(1,:) must be zero.  H(1,1) must always be zero.  Rows 2: contain
      ! user data, which are not stored here, but may be moved around
      ! here.  H(2,1) contains the location of the last empty entry
      ! allocated.  If H(2,1) < 1 or H(2,1) > SIZE(H,2) the subroutine
      ! find_next_E_link_P_to_E resets it to SIZE(H,2) before searching
      ! for an empty space.  The best initial value is SIZE(H,2).
      ! SIZE(H,2) should be 1+prime, or at least 1+(a number with no small
      ! factors).
    logical, intent(in) :: INSERT    ! .TRUE. if G is to be inserted if it
      ! is not found, and .FALSE. otherwise.
    integer, intent(inout) :: LOC    ! Input: defines how a search is to
      ! start or continue. If LOC is zero the search starts at the probe
      ! location defined by applying the hash function to G.  If LOC is
      ! non-zero the search continues at the successor of LOC, where the
      ! successor is the next object having the same value of the hash
      ! function as the object at LOC.  LOC should be zero or the value
      ! returned by an immediately previous call.  This allows duplicate
      ! entries to be disambiguated externally.
      ! Output: the location where an object was found or inserted.
    integer, intent(out) :: STATUS   ! indicates the reason for returning.
      ! Values of STATUS are
      !   FOUND - The object in G was found at H(1,LOC).
      !   INSERTED - The object in G was not found.  If INSERT it was
      !     inserted in H(1,LOC).
      !   FULL - There was no more space to insert an object.
      !   NOT_KEY - LOC not zero, and G was not the key of the object at
      !     H(1,LOC).
      !   BAD_LOC - LOC < 2 or LOC > SIZE(H,2).

    ! *****     Method     *********************************************

    ! The method of hashing using separate chains, with the modification
    ! suggested by Butler Lampson to store the link and the key in the
    ! same space, is used.  The method is described in

    !     Donald E. Knuth, "The Art of Computer Programming, Volume 3:
    !     Sorting and Searching," Addison-Wesley (1973)

    ! on pages 517 - 518, and in exercise 13 on pages 543 - 544, for
    ! which the answer appears on page 689.

    ! Denote the key of an entry A by k(A), and the hash function by
    ! h(k(A)).  We define h(k(A)) = mod(k(A), NT).  Let q = [k(A)/NT].
    ! Then k(A) = q*NT + h(k(A)).  Thus k(A) can be reconstructed from
    ! q, NT and h(k(A)).  Instead of storing k(A) in each entry, we
    ! store W = q*NT + link, where "link" is a pointer to another
    ! entry A' such that h(A) = h(A').  By restricting k(A) > 0, we do
    ! not require extra space to store "link."

    ! We store W for an entry A at P in H(1,P), and use W = 0 to
    ! indicate an empty entry, W > 0 to indicate an entry that starts
    ! a chain (i.e. h(k(A)) = P), and W < 0 to indicate an entry that
    ! does not start a chain.  When W > 0, k(A) = W - mod(W, NT) + P.
    ! When W < 0, follow link = mod(-W, NT) until an entry A' is found
    ! at P' such that W' > 0.  Then k(A) = W - mod(W, NT) + P'.

    ! *****     Local Entities     *************************************

    integer :: E    ! location of the last empty element of H to have been
      ! assigned.  E is local storage for H(2,1).
    integer :: I    ! loop induction variable.
    integer :: K    ! a key computed from G.  K is abs(G) if G /= 0, and
      ! KZERO if G is zero.
    integer :: NE   ! Number of words per entry = SIZE(H,1).
    integer :: NEXTP ! the next value to be used for P.  See P.
    integer :: NPROBE ! the number of probes (key comparisons) used to
      ! discover whether the given entry G is in H.  NPROBE is added onto
      ! NSPROB or NUPROB in common /HCOUNT/ depending on whether the search
      ! was successful.
    integer :: NT   ! Number of entries possible in H = SIZE(H,2) - 1.
    integer :: P    ! index to an element of H.  Used during searching in H.
    integer :: P1   ! the first value of P.
    integer :: W    ! ABS(H(1,P)).

    ! *****     Procedures     *****************************************

    ! Get sizes of H
      ne = size(H,1)
      nt = size(H,2) - 1

    ! Calculate the key to use.  We don't allow non-positive keys in H.

      k = abs(g)
      if ( k == 0 ) k = kzero

    ! Calculate the initial probe location.

      p1 = mod(k,nt)
      status = inserted

    ! Search for the object.

      if ( loc == 0 ) then ! Start searching at P1
         if ( h(1,p1+2) > 0 ) then ! Head of chain.
            p = p1
            call search_chain
            if ( status == found ) return
            if ( insert ) call insert_at_end_of_chain
            return
         end if

         if ( .not.insert ) return
         nuprob=nuprob+1

         if ( h(1,p1+2) < 0) then ! Not head of chain, not empty.
         !  The entry at H(:,P1+2) must be Displaced.
         !  Find P, the predecessor of P1.
            nextp = p1
            do
               p = nextp
               w = abs(h(1,p+2))
               nextp = mod(w,nt)
               nlink = nlink + 1
               if ( nextp == p1 ) exit
            end do
            call find_next_E_link_P_to_E
            if ( status == full ) return
         !  Move the entry at P1 to E.
            nmove = nmove + 1
            h(:,e+2) = h(:,p1+2)
         end if

      !  At this point P1 is empty, either because it was originally
      !  empty or because it was not the head of a chain and was
      !  therefore made empty.

         h(1,p1+2) = k
         loc = p1 + 2
      else ! Continue searching at successor of LOC.
      !  Check whether LOC is within bounds.
         if ( loc < 2 .or. loc > nt+1 ) then
            status = bad_loc
            return
         end if
         p = loc - 2
         w = abs(h(1,p+2))
         nextp = mod(w,nt)
         if ( k-p1 /= w-nextp ) then
            status = not_key
            return
         end if
         if ( nextp /= p1 ) then ! not at the end of the chain
            p = nextp
            call search_chain
            if ( status == found ) return
         end if
         if ( insert ) call insert_at_end_of_chain
      end if
      return

   contains ! in LOOKUP_AND_INSERT

      subroutine find_next_E_link_P_to_E
         e = max(0,min(h(2,1),nt-1))
         h(1,1) = 0 ! so a bounds check isn't needed in the loop below
         do i = 0, 1
            do while ( h(1,e+2) /= 0 )
               e = e - 1
            end do
            if ( e > 0 ) exit
            if ( i /= 0 ) then
               status = full
               return
            end if
            e = nt - 1
         end do
         h(2,1) = e
      !  Since P points to the predecessor of P1, mod(W,NT) = P1.
         h(1,p+2) = sign(w-p1+e,h(1,p+2))
      end subroutine find_next_E_link_P_to_E


      subroutine insert_at_end_of_chain
      !  Reached end of chain without finding K/NT = H(1,P+2)/NT for any
      !  P.  At this point, P points to the last entry in the chain, and
      !  therefore mod(W,NT) = P1.  Insert K in an available space and
      !  link the chain to it.
         call find_next_E_link_P_to_E
         if ( status == full ) return
         h(1,e+2) = -k
         loc = e + 2
         nuprob = nuprob + nprobe
      end subroutine insert_at_end_of_chain


      subroutine search_chain
         nprobe = 0
         do
            w = abs(h(1,p+2))
            nextp = mod(w,nt)
            nprobe = nprobe + 1
            if ( k-p1 == w-nextp ) then
            !  Found K/NT = H(1,P+2)/NT.
               loc = p + 2
               status = found
               nsprob = nsprob + nprobe
               exit
            end if
            if ( nextp == p1 ) exit
            p = nextp
         end do
         ! Here P points to the end of the chain, and therefore the
         ! predecessor of P1
      end subroutine search_chain

   end subroutine LOOKUP_AND_INSERT
end module HASH

! $Log$
