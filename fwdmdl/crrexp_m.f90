module CRREXP_M

  implicit NONE
  private
  public CRREXP, CRREXP_1, SREXP

  interface CRREXP;   module procedure CRREXP_,   ZRREXP_;   end interface
  interface CRREXP_1; module procedure CRREXP_1_, ZRREXP_1_; end interface
  interface REXP;     module procedure SREXP,     DREXP;     end interface


!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    &  "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  pure function CRREXP_ ( Z ) result ( CRREXP )

!{ Evaluate $\frac{e^z-1}z$ as
!  $\frac{\overline{z}}{|z|^2} \left[ e^{i y} ( e^x - 1 ) + e^{i y} - 1 \right ] =
!   \frac{\overline{z}}{|z|^2}
!    \left[ e^{i y} ( e^x - 1 ) + \cos y - 1 + i \sin y \right ]$.  The
!  reason to arrange the computation in this way is that both $e^x-1$ and
!  $\cos y -1$ can be computed without cancellation if care is taken (write
!  a few terms of their Taylor series and you'll agree).

    integer, parameter :: RK = kind(0.0e0)
    complex(rk) :: CRREXP
    complex(rk), intent(in) :: Z
    include 'crrexp.f9h'
  end function CRREXP_

  pure function ZRREXP_ ( Z ) result ( CRREXP )

!{ Evaluate $\frac{e^z-1}z$ as
!  $\frac{\overline{z}}{|z|^2} \left[ e^{i y} ( e^x - 1 ) + e^{i y} - 1 \right ] =
!   \frac{\overline{z}}{|z|^2}
!    \left[ e^{i y} ( e^x - 1 ) + \cos y - 1 + i \sin y \right ]$.  The
!  reason to arrange the computation in this way is that both $e^x-1$ and
!  $\cos y -1$ can be computed without cancellation if care is taken (write
!  a few terms of their Taylor series and you'll agree).

    integer, parameter :: RK = kind(0.0d0)
    complex(rk) :: CRREXP
    complex(rk), intent(in) :: Z
    include 'crrexp.f9h'
  end function ZRREXP_

  pure function CRREXP_1_ (X) result ( CRREXP_1 )
!>> 2002-12-27 WVSnyder  Adapted from SRREXP
!--S replaces "?": ?REXP
! ----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION ( EXP(X) - 1 ) / X
!-----------------------------------------------------------------------
    integer, parameter :: RK = kind(0.0e0)
    complex(rk) :: CRREXP_1
    complex(rk), intent(in) :: X
    include 'crrexp_1.f9h'
  end function CRREXP_1_

  pure function ZRREXP_1_ (X) result ( CRREXP_1 )
!>> 2002-12-27 WVSnyder  Adapted from SRREXP
!--S replaces "?": ?REXP
! ----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION ( EXP(X) - 1 ) / X
!-----------------------------------------------------------------------
    integer, parameter :: RK = kind(0.0d0)
    complex(rk) :: CRREXP_1
    complex(rk), intent(in) :: X
    include 'crrexp_1.f9h'
  end function ZRREXP_1_

  pure function SREXP (X) result ( REXP )
!>> 2002-12-27 SREXP WVS JPL Adapt from Fortran 77
!>> 1994-10-20 SREXP Krogh   Changes to use M77CON
!>> 1994-05-20 SREXP WVS JPL Make SP and DP alike using CHGTYP
!>> 1993-05-06 SREXP WVS JPL Convert from NSWC to Math 77
!--S replaces "?": ?REXP
! ----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION EXP(X) - 1
!-----------------------------------------------------------------------
    integer, parameter :: RK = kind(0.0e0)
    real(rk) :: REXP
    real(rk), intent(in) :: X
    include 'rexp.f9h'
  end function SREXP

  pure function DREXP (X) result ( REXP )
!>> 2002-12-27 SREXP WVS JPL Adapt from Fortran 77
!>> 1994-10-20 SREXP Krogh   Changes to use M77CON
!>> 1994-05-20 SREXP WVS JPL Make SP and DP alike using CHGTYP
!>> 1993-05-06 SREXP WVS JPL Convert from NSWC to Math 77
!--S replaces "?": ?REXP
! ----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION EXP(X) - 1
!-----------------------------------------------------------------------
    integer, parameter :: RK = kind(0.0d0)
    real(rk) :: REXP
    real(rk), intent(in) :: X
    include 'rexp.f9h'
  end function DREXP

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module CRREXP_M

! $Log$
! Revision 2.1  2003/02/04 01:27:32  vsnyder
! Initial commit
!
