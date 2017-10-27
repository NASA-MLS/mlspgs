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
!     on the presence of an optional variable holding the overrding value
! (2) Handle exceptions in a flexible manner
!
! Should we move the Exception-handling stuff to a separate module?
  use MLSFinds, only: FindFirst

! === (start of toc) ===                                                 
!     c o n t e n t s
!     - - - - - - - -
!         Datatypes
! Exception_t        The container for kind, cause, times, etc.
!
!         Functions, operations, routines
! Default            Returns either the value of an optional arg
!                       if it is present, otherwise a default value
! HowWeHandle        Returns one of two values depending on whether an 
!                       optional arg is present; if present sets it to TRUE
! Pass_or_catch      Do we Process an Exception_t datatype, or Pass it on up
!                       the chain of calls?
! Raise              Construct an Exception_t datatype
! === (end of toc) ===

  implicit none
  private

  public :: Default, Default_Char, Default_Double, Default_Integer
  public :: Default_Logical, Default_Single
  public :: HowWeHandle_int, HowWeHandle_log
  public :: Raise, Pass_or_Catch

  interface Default
    module procedure Default_Char, Default_Double, Default_Integer
    module procedure Default_Logical, Default_Single
  end interface

  interface Pass_or_Catch
    module procedure Pass_or_Catch_Int
    module procedure Pass_or_Catch_char
  end interface

  ! The compiler can't disambiguate these; ideas?
  ! interface HowWeHandle
  !  module procedure HowWeHandle_log, HowWeHandle_int
  ! end interface

  ! This is the Exception type
  !   The Raise function creates one
  !   The Pass_or_Catch function processes one
  ! Ideally, what should happen is this:
  ! Say A_1 calls A_2 which calls .. which calls A_n. 
  ! A certain exception occurs in A_n.
  ! A_n Raises an exception by Creating that certain Exception datatype 
  ! and returns it to A_(n-1). If A_(n-1) doesn't handle it itself, it passes
  ! it back up to A_(n-2) which passes it .. finally landing back on 
  ! A_1's doorstep.
  ! So, the operations supplied here are actually Raise and Pass_or_Catch.
  ! The caller is obliged to respond appropriately when the callee returns
  ! an Exception.
  !
  ! Here's some  sample code
  ! use Optional_m, only: Exception_t, Pass_or_Catch
  ! . . .
  ! type(Exception_t) :: Exception
  ! call fragile ( args, Exception )
  ! if ( Exception%condition ) then ! You may not need this extra 'if'
  !    n = Pass_or_catch( Exception, names_we_catch )
  !    select case ( n )
  !    case (1)
  !      call handler_1
  !    case (2)
  !      call handler_2
  !    .   .   .
  !    case default
  !      call clean_up
  !      return ! Pass Exception to Caller
  !    end select
  ! endif
  ! ! Go on our merry way unless the handler cried 'stop'
  ! 
  ! While the lack of elegance abaove is obvious even to the coarsest observer
  ! we've made a start at least.
  ! And, dammit Jim, this is Fortran, not Java.
  public :: Exception_t
  integer, public, parameter :: DefaultCharLen = 32
  type Exception_t
    logical           :: condition   ! If conditional
    character(len=DefaultCharLen) :: name        ! So we can match on names
    character(len=DefaultCharLen) :: where       ! it was Raised originally
    integer           :: id          ! So we can match on id
    integer           :: try_number  ! If not the 1st time
    integer           :: severity    ! E.g., MLSMSG_Error
    integer           :: verboseness ! 0, 1, 2, .. Do we print; how much
    integer           :: int_result  ! Immediate cause
    real              :: sngl_result ! Immediate cause
    character(len=16) :: char_result ! Immediate cause
  end type Exception_t
  
  ! Pass_or_Catch may return these reserved values
  integer, public, parameter :: Pass = 0
  integer, public, parameter :: InternalError = Pass - 1

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

  ! --------------------------------------------------------------------
  ! One way to handle exceptions: does not use the Exception_t datatype
  ! Instead we proceed on a path determined by the presence or absence
  ! of an optional argument
  !

  ! --------------------------- HowWeHandle -----------------------------------
  ! Function used to determine how to handle exceptions (crash or just set flag)
  ! depending on whether an optional flag named 'opt' is present
  ! If opt is not present                     => AbsentValue
  ! If opt is present                         => PresentValue
  !                                               and set opt => 1 or TRUE or ..
  ! E.g., whether to quit with just a warning message or with an error
  ! Assume the optional arg 'fail' is meant to allow for a soft failure:
  ! if an exception occurs, fail is seet to a non-zero value,
  ! a message is printed, and the procedure returns
  ! If instead fail is not present, we quit with a message, as in:
  !   severity = HowWeHandle( fail, MLSMSG_Error, MLSMSG_Warning )
  !   call MLSMessage ( severity, ModuleName, &
  !     & 'Something happened' )
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
  
  ! --------------------------------------------------------------------
  ! The other way to handle exceptions uses the Exception_t datatype
  ! The procedure that detects the exception Raises an instance and
  ! then either Catches it or else Passes it back up the chain of
  ! calls

  ! ------------------- Pass_or_Catch --------------------
  ! The Pass_or_Catch subroutine decides how to process an Exception.
  ! Pass means send it back up one more step in the calling chain
  ! Catch means handle it right away, either dying or recovering somehow.
  ! For now, the decision whether to Catch or Pass is based solely
  ! on matching names or ids.
  ! Can't DumpException here because Dump_1 module uses us, which would
  ! create a circular dependence
  function Pass_or_Catch_int ( exception, names ) result ( which )
    ! Args
    type (Exception_t), intent(in)                            :: Exception
    integer, dimension(:), intent(in)                         :: names
    integer                                                   :: which
    ! Executable
    which = FindFirst( names, Exception%id )
  ! if ( Exception%verboseness > 2 .or. &
  !   & ( Exception%verboseness > 0 .and. which > 0 ) ) &
  !   & call DumpException( exception )
  end function Pass_or_Catch_int

  function Pass_or_Catch_char ( exception, names ) result ( which )
    ! Args
    type (Exception_t), intent(in)                            :: Exception
    character(len=*), dimension(:), intent(in)                :: names
    integer                                                   :: which
    ! Executable
    which = FindFirst( names, Exception%name )
  ! if ( Exception%verboseness > 2 .or. &
  !   & ( Exception%verboseness > 0 .and. which > 0 ) ) &
  !   & call DumpException( exception )
  end function Pass_or_Catch_char

  ! ------------------- Raise --------------------
  ! The Raise function constructs an Exception
  ! The caller can Catch it himself, or Pass it back up the Chain
  function Raise ( name, where, id, condition, try_number, severity, &
    & verboseness, int_result, sngl_result, char_result ) result ( Exception )
    ! Args
    character(len=*)   , optional, intent(in)     :: name
    character(len=*)   , optional, intent(in)     :: where
    integer            , optional, intent(in)     :: id 
    logical            , optional, intent(in)     :: condition  
    integer            , optional, intent(in)     :: try_number 
    integer            , optional, intent(in)     :: severity   
    integer            , optional, intent(in)     :: verboseness   
    integer            , optional, intent(in)     :: int_result 
    real               , optional, intent(in)     :: sngl_result
    character(len=*)   , optional, intent(in)     :: char_result
    type (Exception_t)                            :: Exception
    ! Executable
    Exception%name         = Default( name       , ' '     ) 
    Exception%where        = Default( where      , ' '     ) 
    Exception%id           = Default( id         , 0       ) 
    Exception%condition    = Default( condition  , .false. )
    Exception%try_number   = Default( try_number , 0       )
    Exception%severity     = Default( severity   , 5       )
    Exception%verboseness  = Default( verboseness, 0       )
    Exception%int_result   = Default( int_result , 0       )
    Exception%sngl_result  = Default( sngl_result, 0.e0    )
    Exception%char_result  = Default( char_result, ' '     )
  end function Raise
   
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
