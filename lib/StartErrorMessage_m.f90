! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module StartErrorMessage_m

  implicit NONE
  private
  public StartErrorMessage

  interface StartErrorMessage
    module procedure StartErrorMessage_I, StartErrorMessage_TX
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================

  ! ----------------------------------------  StartErrorMessage_I  -----
  subroutine StartErrorMessage_I ( where )
    use LEXER_CORE, only: PRINT_SOURCE
    use OUTPUT_M, only: OUTPUT
    use TREE, only: WHERE_AT => WHERE
    integer, intent(in) :: Where             ! Tree node index
    call output ( '***** At ' , advance='no' )
    if ( where > 0 ) then
      call print_source ( where_at(where), advance='no', after=': ' )
    else
      call output ( '(no tree available): ', advance='no' )
    end if
  end subroutine StartErrorMessage_I

  ! ---------------------------------------  StartErrorMessage_TX  -----
  subroutine StartErrorMessage_TX ( where )
    use LEXER_CORE, only: PRINT_SOURCE
    use OUTPUT_M, only: OUTPUT
    use TREE, only: TX, WHERE_AT => WHERE
    type(tx), intent(in) :: Where            ! Tree node index
    call output ( '***** At ' , advance='no' )
    if ( where%i > 0 ) then
      call print_source ( where_at(where), advance='no', after=': ' )
    else
      call output ( '(no tree available): ', advance='no' )
    end if
  end subroutine StartErrorMessage_TX

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module StartErrorMessage_m

! $Log$
! Revision 2.1  2013/09/30 23:58:14  vsnyder
! Initial commit after moving from include
!
