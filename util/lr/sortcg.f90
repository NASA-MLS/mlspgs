! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Sort_Configurations

  implicit NONE
  private
  public :: Sortcg

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Sortcg ( Items )

    use Basis_m, only: Item_t
    use Tables, only: PRDIND => Prod_Ind, PRODCN => Productions
    use Toggles, only: Gen, Levels
    use Trace, only: Trace_Begin, Trace_End

    ! Sort the NSETS configurations in the scratch array.  Put reducing
    ! configurations after transitioning configurations.  Order transitioning
    ! configurations according to the symbol after the dot.  For
    ! reducing configurations, or transitioning configurations having the
    ! same symbol after the dot, order according to production number.

    type(item_t), intent(inout) :: Items(:)

    !     *****     Local Variables     ************************************

    ! H       is the increment separating sorted elements.
    ! I       is a loop induction variable and subscript.
    ! IP      is the pointer to the symbol the I iterator is examining.
    !         IP is also used as a temporary variable while exchanging
    !         configurations.
    ! ISYM    is the symbol after the dot in the production at Items(I).
    ! J       is a loop induction variable and subscript.
    ! JP      is the pointer to the symbol the J iterator is examining.
    ! JSYM    is the symbol after the dot in the production at Items(J).
    ! NSets   is the number of elements in Items
    ! V       is the I'th item during sorting
    ! W       is the (J-H)'th item during sorting

    integer H, I, IP, ISym, J, JP, JSym, NSets
    type(item_t) :: V, W

    ! *****     Procedures     *********************************************

    ! Use a shellsort algorithm.

    if ( levels(gen) > 1 ) call trace_begin ( 'SORTCG' )

    nsets = size(items)
    if ( nsets <= 1 ) go to 9
    h = 1
    do while ( h <= nsets )
      h = 3*h + 1
    end do
    do while (h > 1)
      h = h/3
      do i = h+1, nsets
        v = items(i)
        ip = prdind(v%prod)
        if ( v%dot >= prdind(v%prod+1) - ip) then
        ! The I'th configuration is a reducing configuration.
        ! Pretend the symbol after the dot is HUGE.
          isym = huge(1)
        else
          isym = prodcn(ip+v%dot)
        end if
        j = i
        do while ( j > h )
          w = items(j-h)
          jp = prdind(w%prod)
          if ( w%dot >= prdind(w%prod+1) - jp ) then
          ! The J'th configuration is a reducing configuration.
          ! Pretend the symbol after the dot is HUGE.
            jsym = huge(1)
          else
            jsym = prodcn(jp+w%dot)
          end if
          if ( jsym < isym ) exit
          if ( jsym == isym ) then
          ! The symbols after the dot are equal.  Order the
          ! configurations according to the production number.
            if ( w%prod <= v%prod ) exit
          end if
          items(j) = w
          j = j - h
        end do
        items(j) = v
      end do
    end do

  9 if ( levels(gen) > 1 ) call trace_end ( 'SORTCG' )

  end subroutine Sortcg

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Sort_Configurations

! $Log$
! Revision 1.2  2014/01/14 00:11:42  vsnyder
! Revised LR completely
!
! Revision 1.1  2013/10/24 22:41:14  vsnyder
! Initial commit
!
