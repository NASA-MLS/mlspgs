! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module DNWT_CLONE

! Really simple Newton-method with DNWT's interface.  It doesn't even
! have a convergence test, so you need to limit the number of iterations.

  use DNWT_TYPE, only: RK
  use DNWT_MODULE, only: NF_DX, NF_EVALF, NF_EVALJ, NF_NEWX, &
    & NF_SOLVE, NF_START, NWT_T

  interface NWT; module procedure DNWT; end interface
  interface NWTA; module procedure DNWTA; end interface
  interface NWTDB; module procedure DNWTDB; end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
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

  subroutine DNWTDB
  end subroutine DNWTDB
end module DNWT_CLONE

! $Log$
! Revision 2.1  2001/06/01 01:34:52  vsnyder
! Initial commit.  For debugging only
!
