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

! Really simple Newton-method with DNWT's interface.  It doesn't even
! have a convergence test, so you need to limit the number of iterations.

  use DNWT_TYPE, only: RK
  use DNWT_MODULE, only: NF_DX, NF_EVALF, NF_EVALJ, NF_NEWX, &
    & NF_SOLVE, NF_START, NWT_T, FLAGNAME, NF_GMOVE, NF_BEST, NF_AITKEN,&
    & NF_DX_AITKEN, NF_TOLX, NF_TOLF, NF_TOLX_BEST, NF_TOO_SMALL, NF_FANDJ

  interface NWT; module procedure DNWT; end interface
  interface NWTA; module procedure DNWTA; end interface
  interface NWTDB; module procedure DNWTDB; end interface

    
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains
! ***************************************************     DNWT     *****
  subroutine DNWT ( NFLAG, XOPT, NOPT )

    integer, intent(out) :: NFLAG
    real(rk), intent(in) :: XOPT(*)
    integer, intent(in), optional :: NOPT(*)
    nflag = nf_start
    return
  end subroutine DNWT

! **************************************************     DNWTA     *****

  subroutine DNWTA ( NFLAG, AJ )

    integer, intent(inout) :: NFLAG
    type(nwt_t), intent(inout) :: AJ

    select case ( nflag )
    case ( nf_start )
      nflag = nf_evalf
    case ( nf_evalf )
      nflag = nf_evalj
    case ( nf_evalj )
      nflag = nf_solve
    case ( nf_solve )
      nflag = nf_dx
    case ( nf_dx )
      nflag = nf_newx
    case ( nf_newx )
      nflag = nf_evalf
    end select
  end subroutine DNWTA 

! *************************************************     DNWTDB     *****

  subroutine DNWTDB ( AJ, WIDTH, LEVEL, WHY )
    type (NWT_T), intent(in), optional :: AJ
    integer, intent(in), optional :: WIDTH
    integer, intent(in), optional :: LEVEL ! Absent, do everything
    ! Present and 1 do everything
    ! Present and 0 do CDXDXL, CGDX, DXN, FN, FNMIN, INC, SP, SQ, SQRT(FNXE)
    character(len=*), intent(in), optional :: WHY ! printed if present
  end subroutine DNWTDB
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module DNWT_CLONE

! $Log$
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
