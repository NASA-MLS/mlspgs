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
  public :: SORTCG

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine SORTCG (NSETS)
    use SIZES, only: MAXSHD
    use SCRCOM
    use S3, only: PRDIND, PRODCN
    implicit NONE

    ! Sort the NSETS configurations in the scratch array.  Put reducing
    ! configurations after transiting configurations.  Order transiting
    ! configurations according to the symbol after the dot.  For
    ! reducing configurations, or transiting configurations having the
    ! same symbol after the dot, order according to production number.

    integer, intent(in) :: NSETS

    !     *****     Local Variables     ************************************

    ! H       is the increment separating sorted elements.
    ! H3      is 3*H.
    ! I       is a loop induction variable and subscript.
    ! IP      is the pointer to the symbol the I iterator is examining.
    !         IP is also used as a temporary variable while exchanging
    !         configurations.
    ! IS      is SCRTCH(I) (the production number).
    ! ISYM    is the symbol after the dot in the production at SCRTCH(I).
    ! IS1     is SCRTCH(I+1) (the dot position).
    ! IS2     is SCRTCH(I+2) (the context list pointer).
    ! J       is a loop induction variable and subscript.
    ! JP      is the pointer to the symbol the J iterator is examining.
    ! JS      is SCRTCH(J) (the production number).
    ! JSYM    is the symbol after the dot in the production at SCRTCH(J).
    ! JS1     is SCRTCH(J+1) (the dot position).

    integer H, H3, I, IP, IS, ISYM, IS1, IS2, J, JP, JS, JSYM, JS1

    ! *****     Procedures     *********************************************

    ! Use a shellsort algorithm.

    if (nsets <= 3) return
    h = 1
    do while (h <= nsets)
      h = 3*h + 1
    end do
    do while (h > 1)
      h = h/3
      h3 = 3*h
      do i = h3+1, nsets, 3
        is = scrtch(i)
        is1 = scrtch(i+1)
        is2 = scrtch(i+2)
        ip = prdind(is)
        if (is1 .ge. prdind(is+1) - ip) then
        ! The I'th configuration is a reducing configuration.
        ! Pretend the symbol after the dot is MAXSHD+1.
          isym = maxshd + 1
        else
          isym = prodcn(ip+is1)
        end if
        j = i
        do while (j > h3)
          js = scrtch(j-h3)
          js1 = scrtch(j-h3+1)
          jp = prdind(js)
          if (js1 .ge. prdind(js+1) - jp) then
          ! The J'th configuration is a reducing configuration.
          ! Pretend the symbol after the dot is MAXSHD+1.
            jsym = maxshd + 1
          else
            jsym = prodcn(jp+js1)
          end if
          if (jsym .lt. isym) exit
          if (jsym .eq. isym) then
          ! The symbols after the dot are equal.  Order the
          ! configurations according to the production number.
            if (js .le. is) exit
          end if
          scrtch(j) = scrtch(j-h3)
          scrtch(j+1) = scrtch(j-h3+1)
          scrtch(j+2) = scrtch(j-h3+2)
          j = j - h3
        end do
        scrtch(j) = is
        scrtch(j+1) = is1
        scrtch(j+2) = is2
      end do
    end do

  end subroutine SORTCG

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
