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
! Data types and procedures to 
! (1) either return a default or override it depending
!
! === (start of toc) ===                                                 
!     c o n t e n t s
!     - - - - - - - -
!         Datatypes
!
!         Functions, operations, routines
! Default            Returns either the value of an optional arg
!                       if it is present, otherwise a default value
! Print_Default      Print an int array , or a default mesg if not present
! Raise              Construct an Exception_t datatype
! === (end of toc) ===

! === (start of api) ===
! === (end of api) ===

  implicit none
  private

  public :: Default, Default_Char, Default_Double, Default_Integer
  public :: Default_Logical, Default_Single
  public :: Print_Default

  interface Default
    module procedure Default_Char, Default_Double, Default_Integer
    module procedure Default_Logical, Default_Single
  end interface

  interface Print_Default
    module procedure PrintDefault_Intarr
  end interface

  integer, public, parameter :: DefaultCharLen = 32

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  function Default_Char ( Opt, Default, additional ) result ( R )
    ! The optional arg additional (if true) adds on the characters in
    ! Default even if Opt is present
    character(len=*), intent(in), optional           :: Opt
    character(len=*), intent(in)                     :: Default
    logical, intent(in), optional                    :: Additional
    character(len=max(DefaultCharLen, len(default))) :: R
    r = default
    if ( present(opt) ) r = opt
    if ( .not. present(additional) .or. .not. present(opt) ) return
    if ( additional ) r = trim(opt) // trim(Default)
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

  ! ------------------- PrintDefault_intarr --------------------
  ! The PrintDefault subroutine prints either the value
  ! if present, or else a default character string
  subroutine PrintDefault_intarr ( name, ints, mesg )
    ! Args
    character(len=*), intent(in)                     :: name
    integer, dimension (:), optional, intent(in)     :: ints
    character(len=*), intent(in)                     :: mesg
    ! Executable
    if ( present(ints) ) then
      print *, trim(name) // ': ', ints
    else
      print *, trim(name) // ' ' // trim(mesg)
    endif
  end subroutine PrintDefault_intarr

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
! Revision 2.7  2021/09/10 17:42:52  pwagner
! Moved most stuff into new Exception_m
!
! Revision 2.6  2021/02/05 05:10:33  pwagner
! Added a new Print_Default
!
! Revision 2.5  2019/10/31 22:56:51  pwagner
! Moved DumpException here from dump_1
!
! Revision 2.4  2017/10/27 23:13:09  pwagner
! Default_Char can now taked optional arg additional
!
! Revision 2.3  2017/01/19 23:31:03  pwagner
! Add the Exception_t data type and Raise and Pass_or_Catch functions
!
! Revision 2.2  2016/11/03 21:00:16  pwagner
! Added HowWeHandle to determine how to handle an exception
!
! Revision 2.1  2016/09/23 02:55:58  vsnyder
! Initial commit
!
