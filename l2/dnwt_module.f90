! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module DNWT_MODULE
!>> 2000-12-19 W. V. Snyder Removed fwd communication and linear algebra
!>> 2000-03-21 DNWT_MODULE W. V. Snyder Converted to Fortran 90
!--D replaces "?": ?NWT_MODULE, ?NWT_TYPE, ?NWT, ?NWTA, ?NWTDB, ?NWTOP

! All versions use ERMSG.
! ERMSG and ERVN need ERFIN.

! ***** THIS PACKAGE IS STILL UNDER DEVELOPMENT ******

!****************** Program description ********************************

! These subroutines solve f(x)=0 (or find the least square solution)
! where f is a vector with NF components, and x is a vector with NX
! components.

! Reverse communication is used.  Reverse communication means that a
! return to the user's program is made whenever additional information
! is required. When this return is made NFLAG will be set to indicate
! the necessary action.  Negative values indicate that the calling
! program unit should carry out a computation and call NWTA again.
! Positive values indicate that NWTA believes the computation is
! finished.

! Assign values to XOPT( ) and NOPT( ).

! *************************************************
! **  CALL NWT ( NFLAG, XOPT, NOPT)                 **
! *************************************************

!     You will need a Matrix J, and seven vectors:  F, X, "Best X", DX,
!     "Candidate DX" Gradient and "Best Gradient".

!     Set X( ) = initial guess for the solution,
!         AJ%AXMAX = MAXVAL(ABS(X(1:NX)))
!     DO
! *************************************************
! **    CALL NWTA ( NFLAG, AJ )                  **
! *************************************************
!       SELECT CASE ( NFLAG ) ! NFLAG > 0 means "Done", <0 means "Continue"
!       CASE ( NF_EVALF )
!         IF ( too many function values ) EXIT
!         Compute f(x)
!         Compute the Jacobian matrix J if you feel like it
!         Set AJ%FNORM = L2 norm of f(x)
!         IF ( AJ%FNORM is small enough ) EXIT
!       CASE ( NF_EVALJ )
!         IF ( too many Jacobian values ) EXIT
!         Compute the Jacobian matrix J if you didn't do it when NFLAG
!         was NF_EVALF:
!           J(K,L) = Partial of F(K) / W.R.T. X(L), K = 1, NF, L = 1, NX
!         Triangularize J, and compute (negative of the) gradient =
!         -(Jacobian)**T * F.  This is the RHS of the normal equations
!         J**T * J * "Candidate DX" = -J**T * F.            
!         Set
!           AJ%DIAG = element on diagonal with smallest absolute value,
!                   after triangularization,
!           AJ%AJN = maximum L1 norm of column in upper triangle
!                   after triangularization,
!           AJ%FNMIN = L2 norm of residual, ||F + J * "Candidate DX"||
!                   (which can be gotten without solving for
!                   "Candidate DX": if -F is put as the last column of
!                   J before triangularization, either by Householder
!                   or by Cholesky factoring the normal equations,
!                   this is J(N+1,N+1)),
!           AJ%GRADN = L2 norm of Gradient.
!       CASE ( NF_SOLVE )
!         Apply Levenberg-Marquardt stabilization with parameter = AJ%SQ,
!         and solve (Jacobian) * "candidate DX" ~ -F for "candidate DX".
!         Set AJ%FNMIN as for NWT_FLAG = NF_EVALJ, but taking account
!         of Levenberg-Marquardt stabilization.
!         Set
!           AJ%DXN = L2 norm of "candidate DX",
!           AJ%GDX = (Gradient) .dot. ("candidate DX")
!       CASE ( NF_NEWX )
!         Set X = X + DX
!             AJ%AXMAX = MAXVAL(ABS(X)),
!             AJ%BIG = ANY ( DX > 10.0 * epsilon(X) * X )
!       CASE ( NF_GMOVE )
!         Set X = "Best X"
!             DX = AJ%GFAC * "Best Gradient"
!       CASE ( NF_BEST )
!         Set "Best X" = X, "Best Gradient" = Gradient
!       CASE ( NF_AITKEN )
!         Set DX = DX - "Candidate DX",
!             AJ%DXDX = dot_product( DX, DX )
!         IF ( AJ%DXDX /= 0.0 ) &
!           Set AJ%DXDXL = dot_product( DX, "Candidate DX" )
!       CASE ( NF_DX )
!         Set DX = "Candidate DX"
!       CASE ( NF_DX_AITKEN )
!         Set DX = AJ%CAIT * "Candidate DX"
!       CASE ( NF_TOLX, NF_TOLX_BEST, NF_TOLF, NF_TOO_SMALL )
!         IF ( NFLAG == NF_TOO_SMALL ) THEN
!           Take special action if requested accuracy is critical
!         END IF
!         If ( NFLAG == NF_TOLX_BEST ) copy "best X" to "X"
!         Convergence to desired solution.  Do whatever you want to
!         with the solution.
!         EXIT ! unless you have a really good reason to continue
!       CASE ( NF_FANDJ )
!         There is probably an error in the way F or J is computed.
!         A warning has been printed by the error processor.
!         IF ( you have confidence in F and J ) CYCLE
!         STOP
!       END SELECT
!       IF ( you want to return to a previous best X ) NFLAG = NF_START
!     END DO

  use ERMSG_M, only: ERMSG, ERVN
  use DNWT_TYPE, only: RK

  implicit NONE

  private
  public :: RK, DNWT, DNWTA, DNWTDB, DNWTOP
  public :: NWT, NWT_T, NWTA, NWTDB, NWTOP
  public :: FlagName

  type NWT_T               ! Stuff about the problem, neatly packaged.  This
                           ! is the type of the AJ argument of NWTA.  Unless
                           ! otherwise noted, components are inputs to ?NWTA.
                           ! See usage instructions above.
    real(rk) :: AJN        ! Largest L1 norm of column in upper triangle
                           ! of factored Jacobian matrix
    real(rk) :: AXMAX      ! MAXVAL(ABS(X))
    real(rk) :: CAIT = 0.0 ! Aitken parameter -- intent(out)
    real(rk) :: DIAG       ! Smallest | diagonal element | after factoring
    real(rk) :: DXDX       ! dot_product( DX, DX )
    real(rk) :: DXDXL      ! dot_product( "candidate DX", DX )
    real(rk) :: DXN        ! L2 Norm of candidate DX
    real(rk) :: DXNL = 0.0 ! L2 Norm of last candidate DX -- intent(out)
    real(rk) :: FNMIN      ! L2 Norm of F not in column space of the Jacobian
    real(rk) :: FNORM      ! L2 Norm of F at current X
    real(rk) :: GDX        ! dot_product( Gradient, "Candidate DX" )
    real(rk) :: GFAC       ! Amount of best gradient to use for DX
    real(rk) :: GRADN      ! L2 norm of Gradient
    real(rk) :: SQ = 0.0   ! Levenberg-Marquardt parameter -- intent(out)
    real(rk) :: SQT        ! Total Levenberg-Marquardt stabilization -- intent(out)
    logical :: BIG         ! ANY( DX > 10.0 * epsilon(X) * X )
    logical :: STARTING = .true. ! NWTA is still in "starting up" phase -- intent(out)
  end type NWT_T

  ! Start or restart:
  integer, parameter, public :: NF_START = 0
  ! Reasons for returning to user.  See description of usage above.
  ! Reasons to continue
  !      Evaluate F and AJ%FNORM:
  integer, parameter, public :: NF_EVALF = -1
  !      Evaluate J and do other things:
  integer, parameter, public :: NF_EVALJ = nf_evalf-1
  !      Solve for candidate DX:
  integer, parameter, public :: NF_SOLVE = nf_evalj-1
  !      Compute X = X + DX, set AJ%AXMAX = MAXVAL(ABS(X)),
  !      AJ%BIG = ANY ( ABS(DX) > 10.0 * epsilon(X) * ABS(X) );
  !      IF ( .not. AJ%STARTING ) set
  !      aj%dxdxl = dot_product( DX, "Candidate DX" ):
  integer, parameter, public :: NF_NEWX = nf_solve-1
  !      Set X = "Best X", DX = AJ%GFAC * "Best Gradient":
  integer, parameter, public :: NF_GMOVE = nf_newx-1
  !      Set "Best X" = X, "Best Gradient" = Gradient
  integer, parameter, public :: NF_BEST = nf_gmove-1
  !      Set DX = DX - "Candidate DX",
  !          AJ%DXDX = dot_product( DX, DX )
  !      IF ( AJ%DXDX /= 0.0 ) Set AJ%DXDXL = dot_product( DX, "Candidate DX" ):
  integer, parameter, public :: NF_AITKEN = nf_best-1
  !      Set DX = "Candidate DX"
  integer, parameter, public :: NF_DX = nf_aitken-1
  !      Set DX = AJ%CAIT * "Candidate DX"
  integer, parameter, public :: NF_DX_AITKEN = nf_dx-1
  ! Reasons to stop
  !      Convergence due to TOLXA and TOLXR tests:
  integer, parameter, public :: NF_TOLX = 1
  !      Convergence due to TOLXA and TOLXR tests, but solution
  !      is saved "best X":
  integer, parameter, public :: NF_TOLX_BEST = nf_tolx+1
  !      Convergence due to RELSF test:
  integer, parameter, public :: NF_TOLF = nf_tolx_best+1
  !      Convergence due to too small a move:
  integer, parameter, public :: NF_TOO_SMALL = nf_tolf+1
  !      F and J appear not to be consistent:
  integer, parameter, public :: NF_FANDJ = nf_too_small+1

  interface NWT; module procedure DNWT; end interface
  interface NWTA; module procedure DNWTA; end interface
  interface NWTDB; module procedure DNWTDB; end interface
  interface NWTOP; module procedure DNWTOP; end interface

! *****     Private data     *******************************************
  save

  real(rk) :: AJN, AJSCAL, CAIT, CDXDXL, CONDAI, DIAG
  real(rk) :: DXI, DXINC, DXMAXI, DXN, DXNBIG, DXNL, DXNOIS
  real(rk) :: FN, FNB, FNL, FNMIN, FNXE, FRZ, FRZB
  real(rk) :: FRZL, GFAC, GRADN, GRADNB, GRADNL
  integer :: IFL, INC, ITER, ITKEN
  integer :: K1IT, K2IT, KB
  integer :: NFL
  real(rk) :: RELSF, SPACT, SPB, SPFAC, SPG, SPINC, SPL, SPMINI
  real(rk) :: SPSTRT, SQB, SQL, SQMIN
  real(rk) :: TOLXA, TOLXR      ! WVS added 2000-04-05

  character(len=*), parameter :: ME = 'DNWT'

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

!>> 2000-03-21 W. V. Snyder Converted to Fortran 90.
!>> 1988-07-21 C. L. Lawson adapting code for the MATH77 library.
!-----------------------------------------------------------

! Original design and code by FRED T. KROGH
! Modifications for new usage coded by S. SINGLETARY
! Modifications to add bounds, JULY 1977
! JET PROPULSION LAB, PASADENA    JULY, 1974

! Variables in the calling sequence are of the following type

    integer, intent(out) :: NFLAG
    real(rk), intent(in) :: XOPT(*)
    integer, intent(in), optional :: NOPT(*)

! The above parameters are used as follows:

! NFLAG  = Integer used for communication with the user.  NWT sets NFLAG =
!     NF_START.  NFLAG is used by NWTA.  When NWTA returns control to the
!     calling program unit, NFLAG is used as shown above.

! XOPT( ) = Array used for values needed to set options; see NOPT below.

! NOPT( ) = Array of options.
!           If no options are being used, omit NOPT or set NOPT(1) = 0.
!           Otherwise set NOPT(1) /= 0; then NOPT( ) must be a vector,
!           with dimension that depends on the options being specified.
!           Any option below which uses NOPT(K) as a pointer to a loca-
!           tion in XOPT( ) requires NOPT(K) > 7 * NX.  In such a case,
!           IDIMX should be increased as required by the pointer and
!           the option.

!           Define I=1 initially, and after the specification of any
!           option let I=index of next location in NOPT.  Then
!    NOPT(I)
!    = 0    No more options.  The last location in NOPT( ) must always
!           equal 0.
!    = 1    Give the basic internal output for NOPT(I+1) iterations.
!    = 2    Same as 1 above, except AJ is also printed.
!    = 3    Not used.  Do not use.
!    = 4    Not used.  Do not use.
!    = 5    Not used.  Do not use.
!    = 6    Not used.  Do not use.
!    = 7    Not used.  Do not use.
!    = 8    Not used.  Do not use.
!    = 9    Reserved -- currently sets all options to their default values.
!    = 10   Set all options to their default values.
!    = 11   Set SPSTRT = XOPT(NOPT(I+1)), where SPSTRT is the starting
!           value of the normalized Levenberg-Marquardt parameter.
!           NOPT(I) = -11 gives SPSTRT = 0.1  .
!    = 12   Set SPMINI = XOPT(NOPT(I+1)), where SPMINI is the approximate
!           relative error in computing the Jacobian. NOPT(I)= -12 gives
!           SPMINI= 10**(1-D) where D is the number of significant
!           digits in the floating point numbers used.
!    = 13   Set AJSCAL = XOPT(NOPT(I+1)), where AJSCAL is the approximate
!           absolute error in computing the Jacobian. NOPT(I)= -13 gives
!           AJSCAL = 0  .
!    = 14   Set DXMAXI = XOPT(NOPT(I+1)), where DXMAXI is the largest value
!           permitted for DXINC at any time. NOPT(I)= -14 gives
!           DXMAXI = 10**30  .
!    = 15   Set RELSF = XOPT(NOPT(I+1)), where convergence is indicated
!           with NFLAG= 2 if FN - FNMIN < RELSF * FNMIN .
!           NOPT(I) = -15 gives RELSF= 0.01  .
!    = 16   Set DXNIOS = XOPT(NOPT(I+1)), where DXNOIS is the smallest move
!           considered to be significant when the norm of f increases.
!           The initial default is 0, but DXNOIS if updated later.
!    = 17   Set TOLXA = XOPT(NOPT(I+1)), where convergence is indicated
!           with NFLAG = NF_TOLX* if the (possibly scaled) norm of the
!           difference between the current x and the solution is
!           estimated to be .le. TOLXA.  The default value is
!           (underflow limit) ** .75.
!    = 18   Set TOLXR = XOPT(NOPT(I+1)), where convergence is indicated
!           with NFLAG = NF_TOLX* if the (possibly scaled) norm of the
!           difference between the current x and the solution is
!           estimated to be .le. TOLXR * max (abs(X(i))).  The default
!           value is (relative precision) ** .75.
!    = 19   Not used.  Do not use.

!           Any option once set remains set until turned off. In order
!           to return the option with NOPT(I)= K to the nominal state,
!           set NOPT(I)= -K.  The same number of cells in NOPT are requi-
!           red for the -K case, even though any additional cells requi-
!           red by the K case are not referenced.
!???        This appears not to be the way DNWTOP works

!****************** END OF INITIAL COMMENTS ****************

    ifl = nf_start
    nflag = ifl
    if ( present(nopt) ) call nwtop ( nopt, xopt )
    call nwtop ( ) ! default initialization
    return
  end subroutine DNWT

! **************************************************     DNWTA     *****

  subroutine DNWTA ( NFLAG, AJ )

! **************************************************************

! This subroutine is referenced by the user.
! Variables in the calling sequence are defined below.

    integer, intent(inout) :: NFLAG
    type(nwt_t), intent(inout) :: AJ

!                V A R I A B L E   D E F I N I T I O N S


! AJ       About the Jacobian.  See usage instructions above.
! AJN      Norm of the Jacobian
! AJSCAL   Approximate absolute error in computing the Jacobian
! AXMAX    MAXVAL(ABS(X))
! AXMAXB   AXMAX at best X.

! C0       Constant 0
! C1       Constant 1
! C10      Constant 10
! C1P025   Constant 1.025
! C1PXM4   Constant 1.E-4
! C4       Constant 4
! CAIT     Candidate factor for Aitken acceleration
! CBIG     Constant 'BIG' = HUGE()
! CDXDXL   COS of angle between DX and DXL ! ??? DX and DXL aren't defined
! CGDX     COS of angle between gradient and DX
! CONDAI   Crude estimate of reciprocal of condition number of the Jacobian
! CP01     constant .01
! CP125    constant .125
! CP25     constant .25
! CP5      constant .5
! CP76     constant .76
! CP9      constant .9

! DIAG     Minimum diagonal element of triangularized Jacobian
! DXI      Largest factor by which stepsize is allowed to increase
! DXINC    Largest value for DXN currently allowed
! DXMAXI   Largest value permitted for DXINC at any time
! DXN      Norm of current DX (=stepsize)
! DXNBIG   Value of DXN that resulted in a larger FN
! DXNL     Last value of DXN
! DXNOIS   C1PXM4 * largest move that gives a new best X (used to
!          see if increase in FN might be due to too small a step)

! ERMSG    Error message routine

! FN       Norm of F at X from which next step is being taken
! FNB      Value of FN at best X (smallest FN obtained up to this point)
! FNL      Last value of FN
! FNMIN    Norm of component of F orthogonal to the column space of
!          the Jacobian after stabilization
! FNXE     Used to test if F is behaving with near linearity. If
!          FNORM on the next iteration satisfies FNORM**2 <= FNXE
!          then linear behavior is assumed.
! FRZ      Norm of projection of F into the column space of the Jacobian
! FRZB     FRZ at best X
! FRZL     Last value of FRZ

! GFAC     Factor by which the gradient is multiplied when taking a
!          move from the best X
! GRADN    Norm of the gradient
! GRADNB   Value of GRADN at the best X
! GRADNL   Last value of GRADN

! IFL      Internal flag (NFLAG is usually set same as IFL)
!          See negative values of NF_... parameters above.
! INC      Flag set to indicate some past history
!          =-1  Starting
!          =0   Last X == best X
!          =J>0 FN > FNL  J times since last best X
! ITER     Number of current iteration
! ITKEN    Index used to keep track of Aitken accelerations

! K1IT     Counter on number of iterations to give basic internal
!          output
! K2IT     Counter on number of iterations to give extended internal
!          output
! KB       Flag set to indicate some past history
!          =2    Starting
!          =-J   (J>0) Making gradient move from best X (-KB serves as
!                a counter)
!          =1    Last F is bigger than the best F
!          =0    (INC >0) trying a move in wrong direction from the
!                best X
!          =0    (INC =0) new best X

! NFL      Internal flag to indicate reason for returning to user.
!          See positive values of NF_... parameters above.
! NFLAG    Set same as NFL or same as IFL.

! RELSF    Convergence is indicated with NFLAG=NF_TOLF if
!          FN-FNMIN < RELSF*FNMIN
! RND      Used to decide if answer has been determined to machine
!          accuracy (RND is 10 times the round-off level)

! SP       Value of current normalized Marquardt parameter
! SPACT    If SP<SPACT, SQ is set = SQMIN
! SPB      Normalized Marquardt parameter used to get best X
! SPFAC    Factor multiplying last normalized Marquardt parameter
!          to get one estimate to use for current normalized
!          Marquardt parameter
! SPG      Used to initialize SPL after taking a move from best X
!          after a previous failure
! SPINC    Lower bound permitted for normalized Marquardt parameter
!          based on past behavior
! SPL      Used to store normalized Marquardt parameter from last
!          iteration and also the current one
! SPMINI   Nominal minimum value for normalized Marquardt parameter
! SPSTRT   Starting value for normalized Marquardt parameter
! SQ       Actual value of Marquardt parameter. If SQ < 0, modified
!          Aitken acceleration is being used and ABS(SQ) gives the
!          acceleration factor
! SQB      Value of SQL from last best X
! SQL      CP9*AJN*SQL (frequently = CP9* last Marquardt parameter)
! SQMIN    Current minimum value of SQ

! TOLXA    absolute tolerance on X
! TOLXC    tolerance on X
! TP
! TP1      Temporary storage

!-----     Local Variables     ------------------------------------

    real(rk), parameter :: C0 = 0.0_rk
    real(rk), parameter :: C1 = 1.0_rk
    real(rk), parameter :: C10 = 10.0_rk
    real(rk), parameter :: C1P025 = 1.025_rk
    real(rk), parameter :: C1PXM4 = 1.0E-4_rk
    real(rk), parameter :: C4 = 4.0_rk
    real(rk), parameter :: CBIG = huge(c0)
    real(rk), parameter :: CP01 = 0.01_rk
    real(rk), parameter :: CP125 = 0.125_rk
    real(rk), parameter :: CP25 = 0.25_rk
    real(rk), parameter :: CP5 = 0.5_rk
    real(rk), parameter :: CP76 = 0.76_rk
    real(rk), parameter :: CP9 = 0.9_rk

    real(rk) :: AXMAX, AXMAXB, CGDX, SP, SQ, TP, TP1

    real(rk), parameter :: RND = 10.0_rk * epsilon(c0)

!-----     Executable Statements     ------------------------------

! Continuing after return that was expected to be final:
      if ( nflag > 0 ) &
!                nf_tolx nf_tolx_best nf_tolf nf_too_small nf_fandj
        & go to (735,    735,         470,    230,         228), nfl

! Continuing after return for reverse communication:
      if ( nflag < 0 ) &
!                nf_evalf nf_evalj nf_solve nf_newx nf_gmove nf_best
        & go to (160,     260,     540,     850,    240,     740, &
!                nf_aitken nf_dx nf_dx_aitken
        &        750,      770,  770), -ifl

! Initialization
   10 if (ifl /= nf_start) go to 222 ! retreat to best X
      ajn = c0
      condai = c0
      dxnl = c1
      aj%dxnl = c1
      fnl = c0
      inc = -1
      iter = 0
      itken = 0
      kb = 2
      spl = c0
      sq = c0
      axmax = aj%axmax
      aj%dxdxl = c0 ! So it's not undefined, because the user isn't expected
      ! to set it when starting.
   20 ifl = nf_evalf
      nflag = ifl
      return

! Re-enter after evaluating F.

  160 iter = iter + 1
      fn = aj%fnorm

! Test if a retreat to a previous -- best -- X should be made

      ifl = nf_evalj
      tp = spl*ajn*cp9
      if (fn < fnl) then
        if (fn >= fnb) then
          if ( (max(tp,sql) <= sqmin) .and. &
               ((fn*(fn/fnl)**2) > fnb) .and. (cdxdxl >= cp25)) go to 222
        end if
        go to 219
      end if

! Test for convergence

  200 if ( fnl > c0 ) then
        if ( dxn <= tolxa .or. &
           ((sq == c0).or.(spl <= spmini)) .and. dxn <= tolxr * axmax ) then
          nfl = nf_tolx
          if ( kb /= 0 ) then
            aj%fnorm = fnb
            nfl = nf_tolx_best
          end if
          nflag = nfl
          return
        end if
      end if
      if ( fn < fnl ) go to 735
      if ( inc >= 0 ) then
        if ( inc == 0 ) then ! Last X == Best X
          dxnbig = max(dxnl,dxnois)
        else      ! INC > 0 -- Last X not Best X
          if ( kb == 0 ) then
            if (fn < fnb) then
               gfac = -gfac
               if (gfac >=  0) go to 927
            end if
            go to 228
          end if
          if (dxnl <= dxnois) go to 219
          if ( kb < 0 ) go to 224
          if ((tp >= min(sql,sqb+sqb)) .or. (inc >= huge(inc))) go to 222
        end if
        sqmin = min(sqb,max(spl,spact)*ajn*c4)
        dxinc = max(dxnl*cp5,dxnois)
        inc = inc + 1
      end if

  219 sql = tp
      nflag = ifl
      return

! Retreat to a previous -- best -- X

  222 kb = -1
  224 kb = kb - 1
      if (kb == (-4)) then

! Test if Jacobian matrix is being computed properly

        kb = 0
        fnl = c0
        gfac = -gfac/cp01
      end if
      go to 230
  228 gfac = -gfac*cp01
      kb = -1
! Return to best X and step in gradient direction
  230 aj%gfac = gfac
      ifl = nf_gmove
      nflag = ifl
      return

! Re-enter after computing gradient step from best X

  240 continue
      aj%dxn = gradnb*gfac ! In case the caller wants to compute cosines
      if (gfac <= c0) go to 780
      axmax = axmaxb
      dxn = aj%dxn
      gfac = gfac*cp125
      spl = spg
      spg = spl/cp125
      fn = fnb
      fnxe = cp5 * fn **2
      frz = frzb
      gradn = gradnb
      go to 775

! Re-enter after computing Jacobian

  260 diag = aj%diag     ! Smallest | diagonal element | after factoring
      ajn = aj%ajn       ! Largest L1 norm of column in upper triangle
      fnmin = aj%fnmin   ! L2 Norm of F not in column space of the Jacobian
      gradn = aj%gradn   ! L2 norm of gradient
      if (gradn <= c0) then
         gradn = c1
         if (ajn <= c0) ajn=c1
      end if
      condai = max(min(condai,diag/ajn),(diag/ajn)**2)
      spact = cp25*condai
      frz = sqrt(abs((fn-fnmin)*(fn+fnmin)))

      if (iter == 1) then ! First iteration
        sp = spstrt
        spinc = sp
        if (frz <= rnd*fnmin) go to 865
        dxi = cp125
      else
        if (fn >= fnl) then

! Select new stabilizing parameter for case when F has increased

          kb = 1
          if (dxnl < dxnois) then
            sp = c1pxm4*sp
            sqmin = min (sqmin, sp*ajn)
!           dxnl = dxnl / cp01
            go to 520
          end if
          spfac = max((fn/fnl)**2,c10)
          spinc = spl+spl+spact
          dxi = c1
          go to 515
        end if
        if (fn >= fnb) go to 460
      end if

! New best X

      dxinc = dxmaxi
      spinc = min(spinc,max(spmini,abs(ajscal)/ajn))
      sqmin = c0
      spb = spl
      sqb = sql
      if (inc < 0) go to 520
      if (inc /= 0) then
        if (kb < 0) then
          dxi = cp25
        else if (kb == 1) then
          dxi = cp25*dxnl/dxnbig
        end if
        kb = 0
        inc = 0
      end if

! Test for convergence in sense of FN < (RELSF*FNMIN)

  460 if ( inc == 0 ) then
        tp = fn - (c1+relsf)*fnmin
        if ( (tp < c0) .and. &
           & ( (sq == c0) .or. (spl <= spmini) .or. &
           &   ( tp <= -(c1+relsf)*spl*fn) ) ) go to 865 ! May eventually
           ! come back to next statement -- see computed GO TO on NFL near top
      end if

! Select new stabilizing parameter

  470 tp = min((frz/frzl),(fn/fnl))

! Test if F appears almost linear over last step

      if (fn**2 < fnxe) then ! F appears almost linear.
        spfac = cp125*(c1p025-cdxdxl)*tp**2
        if (spl <= spinc) spinc = cp25*spinc
        dxi = min(dxi,cp25)
      else ! F not linear over last step
        spfac = tp*(fn**2)/fnxe
        if (cdxdxl >= cp9) spfac = min(spfac,cp5)
        dxi = c1
      end if
  515 if (gradn > gradnl) spfac = spfac*gradnl/gradn
      sp = max(spinc,min(spb,spl*spfac))
  520 sp = max(sp,sqmin/ajn)
      spl = sp
      sq = sp*ajn
      aj%sqt = sq
      if (sp < spact) sq = sqmin
      gradnl = gradn

! Do Levenberg-Marquardt stabilization, solve for "Candidate DX"

  530 aj%sq = sq
      aj%starting = inc < 0
      ifl = nf_solve
      nflag = ifl
      return

! Re-enter after doing Levenberg-Marquardt stabilization and solving for
! candidate DX

  540 fnmin = aj%fnmin
      dxn = aj%dxn
      cgdx = aj%gdx/(gradn*dxn) ! Cosine ( gradient, dx )
      condai = min(cgdx,gradn/(dxn*ajn**2),diag/ajn)
      fnxe = fnmin**2
      if (sq /= c0) fnxe = fnxe-(spl*ajn*dxn)**2
      tp = dxinc
      if (inc < 0) then
        if (dxn <= dxinc) then
          inc = huge(inc)
          cdxdxl = c0
          go to 735
        end if
      else
        cdxdxl = aj%dxdxl/(dxn*dxnl)
        tp1 = min(cp5,dxi*((c1-cdxdxl)**2))
        if (tp*tp1 > dxnl) tp = dxnl/tp1
        if (dxn <= tp .or. sp >= 1.0e12_rk) then
          if ( inc == 0 ) go to 200
          cait = cbig
          go to 755
        end if
      end if

! Step length is too large

      if (k1it /= 0) then
         write(*,3001) dxn,cdxdxl,spl,sq,cgdx,condai,inc,kb
 3001    format(' DXN=',1PG10.3,'  CDXDXL=',G10.3,'  SPL=', &
     &          G10.3,'  SQ=',G10.3,'  CGDX=',G10.3,        &
     &          '  CI=',G10.3,'  I,K=',I2,',',I2)
      end if
      spinc = sp
      sp = c4*spl+min(condai+condai,min(cp125,condai)*(dxn-tp)/tp)
      sq = sp*ajn
      spl = sqrt(sp**2+spl**2)
      aj%sqt = spl*ajn
      dxi = dxi*cp5
      go to 530

! Store best X, the gradient, and other constants used if a
! return needs to be made to the best X

  735 ifl = nf_best
      nflag = ifl
      return

! Re-enter here after saving X as "best X" and Gradient as "Best gradient"

  740 axmaxb = axmax
      dxnois = max(dxnois,c1pxm4*dxnl)
      gfac = min(cp125*dxn/gradn,ajn**(-2))
      spg = spl+cp01
      fnb = fn
      frzb = frz
!     fnminb = fnmin
      gradnb = gradn
      if ( abs(cdxdxl) < cp9 .or. sq /= c0 .or. dxn >= dxnl ) go to 755

! Compute scalar factor for (modified) Aitken acceleration

      ifl = nf_aitken
      nflag = ifl
      return

! Re-enter here after computing some numbers needed for (modified)
! Aitken acceleration

  750 tp = aj%dxdx
      if (tp /= c0) then
         tp = c1+aj%dxdxl/tp
         tp = min(c10, max(tp,cp01))
         if (abs(cait-tp) < abs(cp25*(cait+tp)-cp5)) then
            cait = tp
            itken = iter+2
            go to 760
         end if
         if (itken < iter) cait = tp
      end if

! End of logic for Aitken acceleration

! Store DX

  755 ifl = nf_dx
      nflag = ifl
      return
! Store the Aitken-modified DX
  760 sq = -cait
      dxn = cait*dxn
      aj%dxn = dxn ! In case the caller wants to compute cosines
      aj%cait = cait
      ifl = nf_dx_aitken
      nflag = ifl
      return

! Re-enter here after DX = "Candidate DX" or DX = CAIT * "Candidate DX"

  770 dxnl = dxn
      aj%dxnl = dxnl
      ! Come here after returning from a gradient move
  775 fnl = fn
      frzl = frz
      fnxe = cp25*fn**2+cp76*fnxe
  780 if (k1it /= 0) then
         k1it = k1it-1
         write(*,3002) iter,fn,frz,fnmin,fnxe,spl,         &
     &                 sq,dxn,gradn,cgdx,cdxdxl,ajn,diag,condai, &
     &                 inc,kb
 3002    format(' ITER=',I4,'  FN=',G10.3,  &
     &          '  FRZ=',G10.3,'  FNMIN=',G10.3,'  FNXE=',    &
     &          G10.3,'  SPL=',G10.3,'  SQ=',G10.3/' DXN=',   &
     &          G10.3,'  GRADN=',G10.3,'  CGDX=',G10.3,       &
     &          '  CDXDXL=',G10.3,'  AJN=',G10.3,'  DIAG=',   &
     &          G10.3,'  CI=',G10.3,'  I,K=',I2,',',I2)
      end if

! Compute new X and test for too small a correction

      aj%starting = inc < 0
      ifl = nf_newx
      nflag = ifl
      return

! Re-enter after computing new X

  850 if ( .not. aj%big ) then
         if ((spl <= spmini).or.(sq == c0))  go to 870
         dxnois = max(dxnois, c10*dxn)
      end if
      axmax = aj%axmax
      go to 20
  865 nfl = nf_tolf
      nflag = nfl
      return

  870 if ((inc > 0) .and. (kb /= 0)) go to 222

! Convergence -- move too small

      nfl = nf_too_small
      nflag = nfl
      return

! Error processing

  927 call ermsg ( me, 2, 0, 'J or F may be in error', ',' )
      nfl = nf_fandj
      nflag = nfl
      return

  end subroutine DNWTA

! *************************************************     DNWTOP     *****
  subroutine DNWTOP ( NOPT, XOPT )
    ! Process option vector for DNWT.  With no arguments, does default
    ! initialization.
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
      valnom(7) = tiny(values) / sqrt(sqrt(tiny(values)))
      valnom(8) = epsilon(values) / sqrt(sqrt(epsilon(values)))
      values(7:8) = valnom(7:8)
      first = .false.
    end if      
    if ( present(nopt) ) then
      if ( .not. present(xopt) ) then
        call ermsg ( me, 99, 2, 'NOPT present but XOPT absent', ',' )
        return
      end if
      i = 1
      do while ( nopt(i) /= 0 )
        k = nopt(i)
        ka = abs(k)
        select case ( ka )
!****************** Change *DATA* values *******************
        case ( 1, 2 ) ! Set K1IT and K2IT (data and common)
          if (k > 0) iopts(ka) = nopt(i+1)
          k1it = max(iopts(1),iopts(2))
          k2it = iopts(2)
        case ( 3, 4, 7, 8 ) ! Change XSCAL and FSCAL indexes in data
                            ! or set up bounds option in data
          if (k > 0) iopts(ka) = nopt(i+1)
        case ( 5, 6 ) ! Set flags in data for reverse communi-
                      ! cation and special matrix operations
          if (k > 0) iopts(ka) = 1
          i = i - 1
        case ( 9, 10 ) ! Set to default values
          iopts = 0
          iopts(8) = 150
          values = valnom
          i = i - 1
        case ( 11:19 ) ! Set in data SPSTRT, SPMINI, AJSCAL,
                       ! DXMAXI, RELSF, and DXNOIS
          values(ka-10) = valnom(ka-10) ! Reset to nominal value
          if (k > 0) & ! If indicated, set to user input value
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
      k2it = iopts(2)
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

! *************************************************     DNWTDB     *****

  subroutine DNWTDB

!   Print the scalars in the module

    namelist /DNWTDB_OUT/ AJN, AJSCAL, CAIT, CDXDXL, CONDAI, DIAG
    namelist /DNWTDB_OUT/ DXI, DXINC, DXMAXI, DXN, DXNBIG, DXNL, DXNOIS
    namelist /DNWTDB_OUT/ FN, FNB, FNL, FNMIN, FNXE, FRZ, FRZB
    namelist /DNWTDB_OUT/ FRZL, GFAC, GRADN, GRADNB, GRADNL
    namelist /DNWTDB_OUT/ IFL, INC, ITER, ITKEN
    namelist /DNWTDB_OUT/ K1IT, K2IT, KB
    namelist /DNWTDB_OUT/ NFL
    namelist /DNWTDB_OUT/ RELSF, SPACT, SPB, SPFAC, SPG, SPINC, SPL, SPMINI
    namelist /DNWTDB_OUT/ SPSTRT, SQB, SQL, SQMIN

    character(len=9) :: IFLname, NFLname

!   write (*,dnwtdb_out)
    write ( *, '(a)' ) &
      & '       AJN        AJSCAL          CAIT        CDXDXL        CONDAI'
    write ( *, '(5es14.7)' ) AJN, AJSCAL, CAIT, CDXDXL, CONDAI
    write ( *, '(a)' ) &
      & '      DIAG           DXI         DXINC        DXMAXI           DXN'
    write ( *, '(5es14.7)' ) DIAG, DXI, DXINC, DXMAXI, DXN
    write ( *, '(a)' ) &
      & '    DXNBIG          DXNL        DXNOIS            FN           FNB'
    write ( *, '(5es14.7)' ) DXNBIG, DXNL, DXNOIS, FN, FNB
    write ( *, '(a)' ) &
      & '       FNL         FNMIN          FNXE           FRZ          FRZB'
    write ( *, '(5es14.7)' ) FNL, FNMIN, FNXE, FRZ, FRZB
    write ( *, '(a)' ) &
      & '      FRZL          GFAC         GRADN        GRADNB        GRADNL'
    write ( *, '(5es14.7)' ) FRZL, GFAC, GRADN, GRADNB, GRADNL
    call flagName ( ifl, iflName ); call flagName ( nfl, nflName )
    write ( *, '(a)' ) &
      & '       IFL        INC    ITER   ITKEN    K1IT    K2IT      KB       NFL'
    write ( *, '(1x,a9,i11,5i8,1x,a9)' ) adjustr(iflName), INC, ITER, ITKEN, &
      & K1IT, K2IT, KB, adjustr(nflName)
    write ( *, '(a)' ) &
      & '     RELSF         SPACT           SPB         SPFAC           SPG'
    write ( *, '(5es14.7)' ) RELSF, SPACT, SPB, SPFAC, SPG
    write ( *, '(a)' ) &
      & '     SPINC           SPL        SPMINI        SPSTRT           SQB'
    write ( *, '(5es14.7)' ) SPINC, SPL, SPMINI, SPSTRT, SQB
    write ( *, '(a)' ) &
      & '       SQL         SQMIN         TOLXA         TOLXR'
    write ( *, '(5es14.7)' ) SQL, SQMIN, TOLXA, TOLXR
  end subroutine DNWTDB

! ***********************************************     FlagName     *****

  subroutine FlagName ( NFlag, ItsName )
  ! Return the name of *NWTA's flag
    integer, intent(in) :: NFlag
    character(len=*), intent(out) :: ItsName
    select case ( nflag )
    case ( NF_EVALF )
      itsName = 'EVALF'
    case ( NF_EVALJ )
      itsName = 'EVALJ'
    case ( NF_SOLVE )
      itsName = 'SOLVE'
    case ( NF_NEWX )
      itsName = 'NEWX'
    case ( NF_GMOVE )
      itsName = 'GMOVE'
    case ( NF_BEST )
      itsName = 'BEST'
    case ( NF_AITKEN )
      itsName = 'AITKEN'
    case ( NF_DX )
      itsName = 'DX'
    case ( NF_DX_AITKEN )
      itsName = 'AITKEN'
    case ( NF_START )
      itsName = 'START'
    case ( NF_TOLX )
      itsName = 'TOLX'
    case ( NF_TOLX_BEST )
      itsName = 'TOLX_BEST'
    case ( NF_TOLF )
      itsName = 'TOLF'
    case ( NF_TOO_SMALL )
      itsName = 'TOO_SMALL'
    case ( NF_FANDJ )
      itsName = 'FANDJ'
    case default
      itsName = 'What???'
    end select
  end subroutine FlagName

end module DNWT_MODULE

! $Log$
! Revision 2.19  2001/06/13 23:57:12  vsnyder
! Use NF_START instead of zero to start IFL
!
! Revision 2.18  2001/06/13 23:50:57  vsnyder
! Correct not-restarting-correctly in DNWT
!
! Revision 2.17  2001/06/01 23:04:11  vsnyder
! Corrected two bugs having to do with Aitken acceleration
!
! Revision 2.16  2001/06/01 01:38:27  vsnyder
! Correct some comments, and some cosmetic changes
!
! Revision 2.15  2001/05/25 20:14:01  vsnyder
! Replace 1p5e format that NAG doesn't like with 5es format
!
! Revision 2.14  2001/05/25 04:58:50  livesey
! Had to comment out some of Van's code to make it compile on NAG.
!
! Revision 2.13  2001/05/24 23:28:06  vsnyder
! Improve internal output; initial value for dxnl=1, not zero
!
! Revision 2.12  2001/05/24 20:23:41  vsnyder
! Keep aj%dxn up to date after Gradient or Aitken moves
!
! Revision 2.11  2001/05/24 18:16:11  vsnyder
! Added NF_START parameter; initially define aj%dxdxl
!
! Revision 2.10  2001/05/22 19:10:33  vsnyder
! Compute aj%starting; improve some comments
!
! Revision 2.9  2001/05/18 01:02:43  vsnyder
! Checking for 'big' move incorrectly
!
! Revision 2.8  2001/05/17 20:15:28  vsnyder
! Make sure aj%dxdxl has an initial value
!
! Revision 2.7  2001/05/12 01:10:24  vsnyder
! Correct 'iter' calculation, add 'FlagName' subroutine
!
! Revision 2.6  2001/05/03 02:00:39  vsnyder
! Insert copyright notice
!
! Revision 2.5  2001/04/28 01:47:39  vsnyder
! Get correct initial value for nflag
!
! Revision 2.4  2001/04/12 00:03:21  vsnyder
! Put 'Total Levenberg-Marquart' in AJ
!
! Revision 2.3  2001/02/16 00:19:34  vsnyder
! Corrected some comments; changed a dummy argument name.
!
! Revision 2.2  2001/02/06 23:33:54  vsnyder
! Initially botched putting in the CVS variables; got it right this time.
!
! Revision 2.1  2001/02/06 23:23:53  vsnyder
! Initial commit
!
