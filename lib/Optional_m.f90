! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Optional_m

  implicit none
  private

  public :: Default, Default_Char, Default_Double, Default_Integer
  public :: Default_Logical, Default_Single

  interface Default
    module procedure Default_Char, Default_Double, Default_Integer
    module procedure Default_Logical, Default_Single
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  function Default_Char ( Opt, Default ) result ( R )
    character(len=*), intent(in), optional :: Opt
    character(len=*), intent(in) :: Default
    character(len=len(default)) :: R
    r = default
    if ( present(opt) ) r = opt
  end function Default_Char

  double precision function Default_Double ( Opt, Default ) result ( R )
    double precision, intent(in), optional :: Opt
    double precision, intent(in) :: Default
    r = default
    if ( present(opt) ) r = opt
  end function Default_Double

  integer function Default_Integer ( Opt, Default ) result ( R )
    integer, intent(in), optional :: Opt
    integer, intent(in) :: Default
    r = default
    if ( present(opt) ) r = opt
  end function Default_Integer

  logical function Default_Logical ( Opt, Default ) result ( R )
    logical, intent(in), optional :: Opt
    logical, intent(in) :: Default
    r = default
    if ( present(opt) ) r = opt
  end function Default_Logical

  real function Default_Single ( Opt, Default ) result ( R )
    real, intent(in), optional :: Opt
    real, intent(in) :: Default
    r = default
    if ( present(opt) ) r = opt
  end function Default_Single

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Optional_m

! $Log$
! Revision 2.1  2016/09/23 02:55:58  vsnyder
! Initial commit
!
