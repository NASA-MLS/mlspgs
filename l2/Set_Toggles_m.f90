! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Set_Toggles_m

  implicit NONE
  public

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Set_Toggles ( String )

    ! Set toggles from String
    ! a => syn -- dump abstract syntax tree
    ! c => con -- trace type checker
    ! f => emit -- forward model
    ! g => gen -- most of everything after type checking
    ! l => lex -- lexer
    ! p => par -- parser
    ! t => tab -- declaration table activity
    ! Invalid toggles are ignored
    ! Any toggle can be followed by a single digit, which sets Levels

    use Toggles, only: Con, Emit, Gen, Lex, Par, Syn, Tab, Toggle, Levels
    character(*), intent(in) :: String

    character(*), parameter :: Good = 'cfglpat'
    integer, parameter :: Index(len(good)) = &
      & (/ con, emit, gen, lex, par, syn, tab /)

    integer :: i, j

    i = 0
    do
      if ( i >= len_trim(string) ) exit
      i = i + 1
      j = scan(good,string(i:i))
      if ( j /= 0 ) then
        toggle(j) = .true.
        if ( string(i+1:i+1) >= '0' .and. string(i+1:i+1) <= '9' ) then
          i = i + 1
          levels(j) = ichar(string(i:i)) - ichar('0')
        end if
      end if
    end do

  end subroutine Set_Toggles

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Set_Toggles_m

! $Log$
! Revision 2.1  2012/06/06 20:36:02  vsnyder
! Initial commit
!
