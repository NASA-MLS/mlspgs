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

! === (start of toc) ===                                                 
!     c o n t e n t s
!     - - - - - - - -
!         Functions, operations, routines
! Default            Returns either the value of an optional arg
!                       if it is present, otherwise a default value
! HowWeHandle        Returns one of two values depending on whether an 
!                       optional arg is present; if present sets it to TRUE
! === (end of toc) ===

  implicit none
  private

  public :: Default, Default_Char, Default_Double, Default_Integer
  public :: Default_Logical, Default_Single
  public :: HowWeHandle_int, HowWeHandle_log

  interface Default
    module procedure Default_Char, Default_Double, Default_Integer
    module procedure Default_Logical, Default_Single
  end interface

  ! interface HowWeHandle
  !  module procedure HowWeHandle_log, HowWeHandle_int
  ! end interface

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

  ! --------------------------------------------------------------------
  ! Function used to determine how to handle exceptions (crash or just set flag)
  ! depending on whether an optional flag named 'opt' is present
  ! If opt is not present                     => AbsentValue
  ! If opt is present                         => PresentValue
  !                                               and set opt => 1 or TRUE or ..
  ! E.g., whether to quit with just a warning message or with an error
  ! Assume the optional arg 'fail' is meant to allow for a soft failure:
  ! if an exception occurs, fail is seet to a non-zero value,
  ! a message is printed, and the procedure returns
  ! If instead fail is not present, we quit with a message
  !   severity = HowWeHandle( fail, MLSMSG_Error, MLSMSG_Warning )
  integer function HowWeHandle_int ( Opt, AbsentValue, PresentValue, OptValue ) &
    & result ( R )
    integer, intent(out), optional :: Opt
    integer, intent(in)            :: AbsentValue
    integer, intent(in)            :: PresentValue
    integer, intent(in), optional  :: OptValue    ! If present, set opt to me
    if ( present(opt) ) then
      r = PresentValue
      opt = 1
      if ( present(OptValue) ) opt = OptValue
    else
      r = AbsentValue
    endif
  end function HowWeHandle_int

  integer function HowWeHandle_log ( Opt, AbsentValue, PresentValue ) result ( R )
    logical, intent(out), optional :: Opt
    integer, intent(in)            :: AbsentValue
    integer, intent(in)            :: PresentValue
    if ( present(opt) )then
      r = PresentValue
      opt = .true.
    else
      r = AbsentValue
    endif
  end function HowWeHandle_log

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
! Revision 2.2  2016/11/03 21:00:16  pwagner
! Added HowWeHandle to determine how to handle an exception
!
! Revision 2.1  2016/09/23 02:55:58  vsnyder
! Initial commit
!
