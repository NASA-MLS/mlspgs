! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module DNWT_CLONE

! Really simple Newton-method with DNWT's interface.

  use DNWT_TYPE, only: RK
  use DNWT_MODULE, only: NF_Best, NF_DX, NF_EvalF, NF_EvalJ, NF_NewX, &
    & NF_Solve, NF_Start, NF_TolX, NF_TolF, NWT_Options, NWT_t

  implicit NONE

  private

  public Alt_NWT, Alt_NWTA, Alt_NWTDB

  interface Alt_NWT; module procedure DNWT; end interface
  interface Alt_NWTA; module procedure DNWTA; end interface
  interface Alt_NWTDB; module procedure DNWTDB; end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! *****     Private data     *******************************************

  character(len=*), parameter :: ME = 'DNWT'

contains

! ***************************************************     DNWT     *****
  subroutine DNWT ( NFlag, AJ, Opt )

    integer, intent(out) :: NFlag
    class(nwt_t), intent(inout) :: AJ
    type(nwt_options), intent(in), optional :: Opt
    nflag = nf_start
    aj%dxnl = huge(aj%dxnl)
    if ( present(opt) ) aj%o = opt
    return
  end subroutine DNWT

! **************************************************     DNWTA     *****

  subroutine DNWTA ( NFlag, AJ )

    integer, intent(inout) :: NFlag
    class(nwt_t), intent(inout) :: AJ

    select case ( nflag )
    case ( nf_start )
      nflag = nf_evalf
    case ( nf_evalf )
      if ( aj%fnorm < (1.0 + aj%o%relsf) * aj%fnmin ) then
        nflag = nf_tolf
      else
        nflag = nf_evalj
      end if
    case ( nf_evalj )
      nflag = nf_solve
    case ( nf_solve )
      if ( aj%dxn <= aj%o%tolxa .or. aj%dxn <= aj%o%tolxr * aj%axmax ) then
        nflag = nf_tolx
      else if ( aj%dxn < aj%dxnl ) then
        aj%dxnl = aj%dxn
        nflag = nf_best
      else
        nflag = nf_dx
      end if
    case ( nf_best )
      nflag = nf_dx
    case ( nf_dx )
      nflag = nf_newx
    case ( nf_newx )
      nflag = nf_evalf
    end select
  end subroutine DNWTA 

! *************************************************     DNWTDB     *****

  subroutine DNWTDB ( AJ, WIDTH, WHY )
    type (NWT_T), intent(in), optional :: AJ
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: WHY ! printed if present
  end subroutine DNWTDB

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module DNWT_CLONE

! $Log$
! Revision 2.11  2019/06/24 23:29:26  pwagner
! Updated to reflect TA-01-143
!
! Revision 2.10  2018/07/23 22:20:21  vsnyder
! Add AJ to all argument lists for nwt* so the options can be (eventually)
! stored there instead of as saved module variables.
!
! Revision 2.9  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.8  2006/09/19 20:33:57  vsnyder
! Add NF_BEST detector so caller's BEST actions get done
!
! Revision 2.7  2006/08/11 20:58:38  vsnyder
! Add 'simple' method to use alternate Newton solver
!
! Revision 2.6  2006/05/27 02:59:36  vsnyder
! Add convergence test
!
! Revision 2.5  2006/05/20 02:17:23  vsnyder
! Make DNWTDB compatible with dnwt_module
!
! Revision 2.4  2005/06/22 18:57:01  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.3  2002/10/07 23:43:11  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.2  2001/06/01 21:28:21  livesey
! Added some more USE's to make it compile
!
! Revision 2.1  2001/06/01 01:34:52  vsnyder
! Initial commit.  For debugging only
!
