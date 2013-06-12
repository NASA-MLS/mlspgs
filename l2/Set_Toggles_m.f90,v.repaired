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

  implicit none
  private

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  public :: SET_TOGGLES

contains

  subroutine Set_Toggles ( String )

    ! Set toggles and toggle levels from String
    ! a => syn -- dump abstract syntax tree
    ! c => con -- trace type checker
    ! f => emit -- forward model
    ! g => gen -- most of everything after type checking
    ! l => lex -- lexer
    ! p => par -- parser
    ! t => tab -- declaration table activity
    ! Invalid toggles are ignored
    ! Any toggle can be followed by a single digit, which sets Levels

    use TOGGLES, only: CON, EMIT, GEN, LEX, PAR, SYN, TAB, TOGGLE, LEVELS
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
        toggle(index(j)) = .true.
        levels(index(j)) = 0
        if ( i+1 > len_trim(string) ) exit
        if ( string(i+1:i+1) >= '0' .and. string(i+1:i+1) <= '9' ) then
          i = i + 1
          levels(index(j)) = ichar(string(i:i)) - ichar('0')
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
! Revision 2.4  2013/06/12 02:39:21  vsnyder
! Use INDEX array to ensure correct subscripts
!
! Revision 2.3  2012/10/20 00:00:58  pwagner
! Fixed bug added with last bugfix
!
! Revision 2.2  2012/06/21 00:31:29  pwagner
! Protect against exceeding substring end
!
! Revision 2.1  2012/06/06 20:36:02  vsnyder
! Initial commit
!
