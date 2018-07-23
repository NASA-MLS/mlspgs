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
  use DNWT_MODULE, only: NF_BEST, NF_DX, NF_EVALF, NF_EVALJ, NF_NEWX, &
    & NF_SOLVE, NF_START, NWT_T, FLAGNAME, NF_GMOVE, NF_BEST, NF_AITKEN,&
    & NF_DX_AITKEN, NF_TOLX, NF_TOLF, NF_TOLX_BEST, NF_TOO_SMALL, NF_FANDJ

  implicit NONE

  private

  public ALT_NWT, ALT_NWTA, ALT_NWTDB, ALT_NWTOP

  interface ALT_NWT; module procedure DNWT; end interface
  interface ALT_NWTA; module procedure DNWTA; end interface
  interface ALT_NWTDB; module procedure DNWTDB; end interface
  interface ALT_NWTOP; module procedure DNWTOP; end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! *****     Private data     *******************************************
  save

  ! Some of these aren't used.  They're here so that the DNWTOP from the
  ! "real" DNWT can be used without modification.

  real(rk) :: AJSCAL, DXMAXI, DXN_PREV, DXNOIS, RELSF, SPMINI, SPSTRT, TOLXA, TOLXR
  integer :: K1IT

  character(len=*), parameter :: ME = 'DNWT'

contains
! ***************************************************     DNWT     *****
  subroutine DNWT ( NFLAG, AJ, XOPT, NOPT )

    integer, intent(out) :: NFLAG
    type(nwt_t), intent(inout) :: AJ
    real(rk), intent(in) :: XOPT(*)
    integer, intent(in), optional :: NOPT(*)
    nflag = nf_start
    dxn_prev = huge(dxn_prev)
    if ( present(nopt) ) call alt_nwtop ( aj, nopt, xopt )
    call alt_nwtop ( aj ) ! default initialization
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
      if ( aj%fnorm < (1.0 + relsf) * aj%fnmin ) then
        nflag = nf_tolf
      else
        nflag = nf_evalj
      end if
    case ( nf_evalj )
      nflag = nf_solve
    case ( nf_solve )
      if ( aj%dxn <= tolxa .or. aj%dxn <= tolxr * aj%axmax ) then
        nflag = nf_tolx
      else if ( aj%dxn < dxn_prev ) then
        dxn_prev = aj%dxn
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

  subroutine DNWTDB ( AJ, WIDTH, LEVEL, WHY )
    type (NWT_T), intent(in), optional :: AJ
    integer, intent(in), optional :: WIDTH
    integer, intent(in), optional :: LEVEL ! Absent, do everything
    ! Present and 1 do everything
    ! Present and 0 do CDXDXL, CGDX, DXN, FN, FNMIN, INC, SP, SQ, SQRT(FNXE)
    character(len=*), intent(in), optional :: WHY ! printed if present
  end subroutine DNWTDB

! *************************************************     DNWTOP     *****
  subroutine DNWTOP ( AJ, NOPT, XOPT )
    ! Process option vector for DNWT.  With no arguments, does default
    ! initialization.
    use ERMSG_M, only: ERMSG, ERVN
    type (NWT_T), intent(in) :: AJ
    integer, intent(in), optional :: NOPT(*)
    real(rk), intent(in), optional :: XOPT(*)

    logical, save :: FIRST = .true.
    integer I, INDIC, K, KA, NACT
    integer, save :: IOPTS(8) = (/ 0, 0, 0, 0, 0, 0, 0, 150 /)
    character(len=4), parameter :: LABL(2) = (/ '(I) ', 'NOPT' /)
    ! ??? Defaults aren't what comments say
    real(rk),save :: VALUES(9) = (/ 0.1_rk, epsilon(values), epsilon(values), &
                                    huge(values), 0.01_rk, 0.0_rk, 0.0_rk, &
                                    0.0_rk, 0.0_rk /)
    real(rk),save :: VALNOM(9) = (/ 0.1_rk, epsilon(values), epsilon(values), &
                                    huge(values), 0.01_rk, 0.0_rk, 0.0_rk, &
                                    0.0_rk, 0.0_rk /)
    integer :: IVAL(2) ! for error messages

!*************** Start of executable code ******************

    if ( first ) then
      valnom(7) = sqrt(sqrt(tiny(values))) ** 3
      valnom(8) = sqrt(sqrt(epsilon(values))) ** 3
      values(7:8) = valnom(7:8)
      first = .false.
    end if
    if ( present(nopt) ) then
      if ( .not. present(xopt) ) then
        call ermsg ( me, 99, 2, 'NOPT present but XOPT absent', '.' )
        return
      end if
      i = 1
      do while ( nopt(i) /= 0 )
        k = nopt(i)
        ka = abs(k)
        select case ( ka )
!****************** Change *DATA* values *******************
        case ( 1, 2 ) ! Set K1IT
          if ( k > 0) iopts(ka) = nopt(i+1)
          k1it = max(iopts(1),iopts(2))
        case ( 3, 4, 7, 8 ) ! Change XSCAL and FSCAL indexes in data
                            ! or set up bounds option in data
          if ( k > 0) iopts(ka) = nopt(i+1)
        case ( 5, 6 ) ! Set flags in data for reverse communi-
                      ! cation and special matrix operations
          if ( k > 0) iopts(ka) = 1
          i = i - 1
        case ( 9, 10 ) ! Set to default values
          iopts = 0
          iopts(8) = 150
          values = valnom
          i = i - 1
        case ( 11:19 ) ! Set in data SPSTRT, SPMINI, AJSCAL,
                       ! DXMAXI, RELSF, and DXNOIS
          values(ka-10) = valnom(ka-10) ! Reset to nominal value
          if ( k > 0) & ! If indicated, set to user input value
            values(ka-10) = xopt(nopt(i+1))
        case default
          indic = 1
          nact = 2
          call ermsg ( me, indic, nact, 'INVALID NOPT', ',' )
          ival(1) = i
          ival(2) = nopt(i)
          call ervn ( labl, ival, '.' )
          return
        end select
        i = i + 2
      end do
    else
      k1it = max(iopts(1),iopts(2))
    end if
!                     Reset every time one of them changes
    spmini = values(2)
    ajscal = values(3)
    dxmaxi = values(4)
    relsf = max(values(5),epsilon(relsf))
    tolxa = values(7)
    tolxr = values(8)
    if ( present(nopt) ) return
!                     Only set during initialization
    spstrt = values(1)
    dxnois = values(6)
    return
  end subroutine DNWTOP

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
