!*==dnwt_module.f    processed by SPAG 7.20RA at 19:07 on 16 Jul 2018
! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module DNWT_Module
!>> 2000-12-19 W. V. Snyder Removed fwd communication and linear algebra
!>> 2000-03-21 DNWT_Module W. V. Snyder Converted to Fortran 90
!>> 2018-07-17 W. V. Snyder Restructured
!--D replaces "?": ?NWT_MODULE, ?NWT_TYPE, ?NWT, ?NWTA, ?NWTDB, ?NWTOP

! All versions use ERMSG.
! ERMSG and ERVN need ERFIN.

! ***** THIS PACKAGE IS STILL UNDER DEVELOPMENT ******

!****************** Program description ********************************

! These subroutines solve F(X)=0 (or find the least square solution)
! where F is a vector with NF components, and X is a vector with NX
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
! **  CALL NWT ( NFLAG, AJ, XOPT, NOPT)                 **
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
!         Compute F(X)
!         Compute the Jacobian matrix J if you feel like it
!         Set AJ%FNORM = L2 norm of F(X)
!         IF ( AJ%FNORM is small enough ) EXIT
!       CASE ( NF_EVALJ )
!         IF ( too many Jacobian values ) EXIT
!         Compute the Jacobian matrix J if you didn't do it when NFLAG
!         was NF_EVALF:
!           J(K,L) = Partial of F(K) / W.R.T. X(L), K = 1, NF, L = 1, NX
!         Compute (negative of the) gradient =
!         -(Jacobian)**T * F.  This is the RHS of the normal equations
!         J**T * J * "Candidate DX" = -J**T * F.
!         Triangularize J.
!         Set
!           AJ%DIAG = element on diagonal with smallest absolute value,
!                   after triangularization,
!           AJ%AJN = maximum L1 norm of column in upper triangle
!                   after triangularization,
!           AJ%FNMIN = L2 Norm of F not in column space of the Jacobian after
!                   stabilization, ||F + J * "Candidate DX"||, (which can be
!                   gotten without solving for "Candidate DX": if -F is put as
!                   the last column of J before triangularization, either by
!                   Householder or by Cholesky factoring the normal equations,
!                   this is J(N+1,N+1)),
!           AJ%GRADN = L2 norm of Gradient.
!       CASE ( NF_LEV )
!         Compute quantities necessary to determine the Levenberg-Marquardt
!         stabilization parameter: Solve U^T q = dX where U is the Cholesky
!         factor of J^T J, compute |q|**2 and store it in AJ%QNSQ.  This is
!         for the Moré-Sorensen algorithm.
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
!         Best solution so far.
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
!         IF ( you don't have confidence in F and J ) STOP
!         If you don't stop, a gradient move is forced
!           (next NFLAG == NF_Gmove).
!       END SELECT
!       IF ( you want to return to a previous best X ) NFLAG = NF_START
!     END DO

  use ERMSG_M, only: Ermsg, Ervn
  use DNWT_Type, only: RK
  use Output_M, only: NewLine, Output

  implicit NONE

  private
  public :: RK, DNWT, DNWTA, DNWTDB, DNWTOP, DNWT_GUTS
  public ::     NWT,  NWTA,  NWTDB,  NWTOP,  NWT_GUTS, NWT_T
  public :: FlagName

  ! Start or restart:
  integer, parameter, public :: NF_START = 0
  ! Reasons for returning to user.  See description of usage above.
  ! Reasons to continue
  !      Evaluate F and AJ%FNORM:
  integer, parameter, public :: NF_EVALF = NF_START - 1
  !      Evaluate J and do other things:
  integer, parameter, public :: NF_EVALJ = nf_evalf-1
  !      Calculate quantities necessary to determine the Levenberg-Marquardt
  !      stabilization parameter
  integer, parameter, public :: NF_LEV = nf_evalj-1
  !      Solve for candidate DX:
  integer, parameter, public :: NF_SOLVE = nf_lev-1
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
  integer, parameter, public :: NF_TOLX = NF_START + 1
  !      Convergence due to TOLXA and TOLXR tests, but solution
  !      is saved "best X":
  integer, parameter, public :: NF_TOLX_BEST = nf_tolx+1
  !      Convergence due to RELSF test:
  integer, parameter, public :: NF_TOLF = nf_tolx_best+1
  !      Convergence due to too small a move:
  integer, parameter, public :: NF_TOO_SMALL = nf_tolf+1
  !      F and J appear not to be consistent:
  integer, parameter, public :: NF_FANDJ = nf_too_small+1
  !      Indicate not ready to return:
  integer, parameter, private :: NF_NOTHING = nf_fandj+1
  !      In case you want to add some flags of your own...
  integer, parameter, public :: NF_BIGGEST_FLAG = NF_NOTHING
  integer, parameter, public :: NF_SMALLEST_FLAG = NF_DX_AITKEN

  real(rk), parameter, private :: EPS = epsilon(1.0_rk), CBIG = huge(1.0_rk)

  type NWT_Options    ! Options that can be set by DNWTOP
    real(rk) :: AJSCAL = eps         ! Approximate absolute error in computing the
                                     ! Jacobian
    real(rk) :: DXMAXI = 0.1*cbig    ! Largest value permitted for DXINC at any time
    real(rk) :: DXNOIS = 0.0         ! 1.0E-4_rk * largest move that gives a new best X
                                     ! (used to see if increase in AJ%FNORM might be
                                     ! due to too small a step)
    integer :: K1IT = 0              ! Number of iterations to give basic internal output
    real(rk) :: RELSF = max(0.01_rk,eps) ! Convergence is indicated with NFLAG = NF_TOLF if
                                     ! AJ%FNORM - AJ%FNMIN < RELSF*AJ%FNMIN
    real(rk) :: SPMINI = eps         ! Nominal minimum value for normalized Marquardt
                                     ! parameter
    real(rk) :: SPSTRT = 0.1_rk      ! Starting value for normalized Marquardt parameter
    real(rk) :: TOLXA = 0.0          ! Absolute tolerance on X.  Convergence is indicated
                                     ! with NFLAG = NF_TOLX* if the (possibly scaled)
                                     ! norm of the difference between the current x and
                                     ! the solution is estimated to be <= TOLXA.
    real(rk) :: TOLXR = 0.0          ! Relative tolerance on X.  Convergence is indicated
                                     ! with NFLAG = NF_TOLX* if the (possibly scaled)
                                     ! norm of the difference between the current x and
                                     ! the solution is estimated to be <= TOLXR * max
                                     ! (abs(X(i))).
  end type NWT_Options

  type(nwt_options), parameter :: Default_Options = NWT_options()

  type NWT_T               ! Stuff about the problem, neatly packaged.  This
                           ! is the type of the AJ argument of NWTA.
    ! Inputs to ?NWTA.  See usage instructions above.
    real(rk) :: AJN        ! Largest L1 norm of column in upper triangle
                           ! of factored Jacobian matrix
    real(rk) :: AXMAX      ! MAXVAL(ABS(X))
    real(rk) :: DIAG       ! Smallest | diagonal element | after factoring
    real(rk) :: DXDX       ! dot_product( DX, DX )
    real(rk) :: DXDXL      ! dot_product( "candidate DX", DX )
    real(rk) :: DXN        ! L2 Norm of candidate DX (Newton step length)
    real(rk) :: FNMIN      ! L2 Norm of F not in column space of the Jacobian
    real(rk) :: FNORM      ! L2 Norm of F at current X
    real(rk) :: GDX        ! dot_product( Gradient, "Candidate DX" )
    real(rk) :: GFAC       ! Factor by which the gradient is to be multiplied
                           ! when taking a gradient move from the best X
    real(rk) :: GRADN      ! L2 norm of Gradient
    real(rk) :: QNSQ       ! Square of norm of solution of U^T q = x
    real(rk) :: SQMIN = 0.0 ! Current minimum value of SQ ; recommend to use
                           ! this as a floor for Levenberg-Marquardt
                           ! stabilization in the caller
    logical :: BIG         ! ANY( DX > 10.0 * epsilon(X) * X )
    ! Outputs from ?NWTA
    real(rk) :: CAIT = 0.0 ! Candidate factor for Aitken acceleration
    real(rk) :: DXINC      ! Largest value for DXN currently allowed
    real(rk) :: SQ = 0.0   ! Actual value of Marquardt parameter. If SQ < 0,
                           ! modified Aitken acceleration is being used and
                           ! ABS(SQ) gives the acceleration factor
    real(rk) :: SQT = 0.0  ! Total Levenberg-Marquardt stabilization
    logical :: STARTING = .true. ! NWTA is still in "starting up" phase
    ! Values retained from return until next call
    real(rk) :: AXMAXB     ! AXMAX at best X (not actually used)
    real(rk) :: CDXDXL     ! COS of angle between DX (current move) and
                           ! DXL (previous move)
    real(rk) :: CONDAI     ! Crude estimate of reciprocal of condition number
                           ! of the Jacobian
    real(rk) :: DXBAD      ! This size for DX is likely a disaster
    real(rk) :: DXFAIL     ! Like DXINC, but more so.
    real(rk) :: DXI        ! Largest factor by which stepsize is allowed to
                           ! increase
    real(rk) :: DXNBIG     ! A size of move that should be safe.
    real(rk) :: DXNL       ! L2 Norm of last candidate DX
    real(rk) :: FNMINB     ! FNMIN at best X, not best FNMIN
    real(rk) :: FNORMB     ! L2 Norm of F at best X
    real(rk) :: FNORML     ! Last value of AJ%FNORM
    real(rk) :: FNXE       ! Used to test if F is behaving with near linearity.
                           ! If ( FNORM**2 on the next iteration ) <= FNXE
                           ! then linear behavior is assumed.
    real(rk) :: FRZ        ! Norm of projection of F into the column space of
                           ! the Jacobian
    real(rk) :: FRZB       ! FRZ at best X (not actually used)
    real(rk) :: FRZL       ! Last value of FRZ
    real(rk) :: GRADNB     ! L2 norm of Gradient at best F
    real(rk) :: GRADNL     ! Last value of GRADN
    integer :: KFAIL       ! Failure count: 0 no failures; k > 0 we have
                           ! successfully taken k - 1 steps of size DXBAD.
    real(rk) :: SP         ! Current normalized Marquardt parameter
    real(rk) :: SPACT      ! Set to 0.25 * CONDAI after evaluating J (but not
                           ! after CONDAI is calculated after solving for DX).
                           ! If SP <SPACT, SQ is set = SQMIN
    real(rk) :: SPB        ! Normalized Marquardt parameter used to get best X
    real(rk) :: SPFAC      ! Factor multiplying last normalized Marquardt
                           ! parameter to get one estimate to use for current
                           ! normalized Marquardt parameter
    real(rk) :: SPG        ! Used to initialize SPL after taking a gradient
                           ! move from best X after a previous failure
    real(rk) :: SPINC      ! Lower bound permitted for normalized Marquardt
                           ! parameter based on past behavior
    real(rk) :: SPL        ! Normalized Marquardt parameter from last
                           ! iteration, and then the current one
    real(rk) :: SQB        ! Value of SQL from last best X
    real(rk) :: SQL        ! 0.9 * AJN * SPL 
                           ! (frequently = 0.9 * last Marquardt parameter)
    integer :: WHERE_TO = NF_Start   ! Internal value of NFLAG, to remember
                           ! reason for return, and therefore where to go next.

    ! Temporary only, but here for printing by ?NWTDB
    real(rk) :: CGDX       ! COS of angle between gradient and DX
  end type NWT_T

  interface NWT
    module procedure DNWT
  end interface
  interface NWTA
    module procedure DNWTA
  end interface
  interface NWTDB
    module procedure DNWTDB
  end interface
  interface NWTOP
    module procedure DNWTOP
  end interface

! *****     Private data     *******************************************
  save

  ! Variables that can be set from the option vector:

  real(rk) :: AJSCAL
  real(rk) :: DXMAXI
  real(rk) :: DXNOIS
  real(rk) :: RELSF
  real(rk) :: SPMINI, SPSTRT
  real(rk) :: TOLXA, TOLXR      ! WVS added 2000-04-05
  integer :: K1IT

  integer :: INC, ITER, ITKEN
  integer :: KB

  type NWT_GUTS    ! Public type returned by DNWT_GUTS
    real(rk) :: AJSCAL
    real(rk) :: DXMAXI, DXNOIS
    integer :: INC, ITER, ITKEN
    integer :: K1IT, KB
    real(rk) :: RELSF, SPMINI
    real(rk) :: SPSTRT
    real(rk) :: TOLXA, TOLXR      ! WVS added 2000-04-05
  end type NWT_GUTS

  character(len=*), parameter :: ME = 'DNWT'

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains
! ***************************************************     DNWT     *****
  subroutine DNWT ( NFLAG, AJ, XOPT, NOPT )

!>> 2000-03-21 W. V. Snyder Converted to Fortran 90.
!>> 1988-07-21 C. L. Lawson adapting code for the MATH77 library.
!-----------------------------------------------------------

! Original design and code by FRED T. KROGH
! Modifications for new usage coded by S. SINGLETARY
! Modifications to add bounds, JULY 1977
! JET PROPULSION LAB, PASADENA    JULY, 1974
! Re-structured by VAN SNYDER, JPL, 17 July 2018

! Variables in the calling sequence are of the following types

    integer, intent(out) :: NFLAG
    type(nwt_t), intent(inout) :: AJ
    real(rk), intent(in) :: XOPT(*)
    integer, intent(in), optional :: NOPT(*)

! The above arguments are used as follows:

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
!           with NFLAG = NF_TOLF if AJ%FNORM - AJ%FNMIN < RELSF * AJ%FNMIN .
!           NOPT(I) = -15 gives RELSF = 0.01  .
!    = 16   Set DXNIOS = XOPT(NOPT(I+1)), where DXNOIS is the smallest move
!           considered to be significant when the norm of f increases.
!           The initial default is 0, but DXNOIS is updated later.
!    = 17   Set TOLXA = XOPT(NOPT(I+1)), where convergence is indicated
!           with NFLAG = NF_TOLX* if the (possibly scaled) norm of the
!           difference between the current x and the solution is
!           estimated to be <= TOLXA.  The default value is
!           (underflow limit) ** .75.
!    = 18   Set TOLXR = XOPT(NOPT(I+1)), where convergence is indicated
!           with NFLAG = NF_TOLX* if the (possibly scaled) norm of the
!           difference between the current x and the solution is
!           estimated to be <= TOLXR * max (abs(X(i))).  The default
!           value is (relative precision) ** .75.
!    = 19   Not used.  Do not use.

!           Any option once set remains set until turned off. In order
!           to return the option with NOPT(I)= K to the nominal state,
!           set NOPT(I)= -K.  The same number of cells in NOPT are requi-
!           red for the -K case, even though any additional cells requi-
!           red by the K case are not referenced.
!???        This appears not to be the way DNWTOP works

!****************** END OF INITIAL COMMENTS ****************

    nflag = nf_start
    aj%where_to = nf_start
    if ( present(nopt) ) call nwtop ( aj, nopt, xopt )
    call nwtop ( aj ) ! default initialization

  end subroutine DNWT

! **************************************************     DNWTA     *****

  subroutine DNWTA ( NFLAG, AJ )

! **********************************************************************

! This subroutine is referenced by the user.
! Arguments in the calling sequence are defined below.

    integer, intent(inout) :: NFLAG ! See usage instructions above
    type(nwt_t), intent(inout) :: AJ

!                V A R I A B L E   D E F I N I T I O N S


! AJ       About the Jacobian.  See usage instructions above.
! AJSCAL   Approximate absolute error in computing the Jacobian

! CBIG     Constant 'BIG' = HUGE()

! DXMAXI   Largest value permitted for DXINC at any time
! DXNOIS   1.0E-4_rk * largest move that gives a new best X (used to
!          see if increase in AJ%FNORM might be due to too small a step)

! INC      Flag set to indicate some past history
!          =-1  Starting
!          =0   Last X == best X
!          =J>0 AJ%FNORM > AJ%FNORML  J times since last best X
! ITER     Number of current iteration
! ITKEN    Index used to keep track of Aitken accelerations

! K1IT     Number of iterations to give basic internal output
! KB       Flag set to indicate some past history
!          =2    Starting
!          =-J   (J>0) Making gradient move from best X (-KB serves as
!                a counter)
!          =1    Last F is bigger than the best F
!          =0    (INC >0) trying a move in wrong direction from the
!                best X
!          =0    (INC =0) new best X

! RELSF    Convergence is indicated with NFLAG = NF_TOLF if
!          AJ%FNORM - AJ%FNMIN < RELSF*AJ%FNMIN
! RND      Used to decide if answer has been determined to machine
!          accuracy (RND is 10 times the round-off level)

! SPMINI   Nominal minimum value for normalized Marquardt parameter
! SPSTRT   Starting value for normalized Marquardt parameter

! TOLXA    absolute tolerance on X
! TOLXC    tolerance on X
! TP       Temporary storage

!-----     Parameters and Unsaved Local Variables     ------------------

    integer, parameter :: INCBIG = huge(inc) / 2

    real(rk) :: TP       ! Widely used temporary variable

    real(rk), parameter :: RND = 10.0_rk * epsilon(0.0)

    ! Values for Before_Solve argument of Do_NF_Solve
    integer, parameter :: Solve_Only = 0
    integer, parameter :: Set_New_Marquardt = 1 + 2 * Solve_Only
    integer, parameter :: Test_GradN_Increase = 1 + 2 * Set_New_Marquardt
    integer, parameter :: Test_F_Linear = 1 + Test_GradN_Increase

!-----     Executable Statements     -----------------------------------

    if ( nflag == nf_start ) then
      ! Either initialization, or user is forcing retreat to best X
      ! followed by gradient move

      if ( aj%where_to == nf_start ) then
      ! Initialization
        aj%ajn = 0.0
        aj%condai = 0.0
        aj%dxnl = 1.0
        aj%fnormb = sqrt(huge(aj%fnormb))
        aj%dxfail = huge(aj%dxfail)
        aj%dxbad = aj%dxfail
        aj%kfail = 0
        aj%fnorml = 0.0
        inc = -1 ! Starting
      ! inc = 0  ! Pretend we're at a "best"
        iter = 0
        itken = 0
        kb = 2
        aj%spl = 0.0
        aj%sq = 0.0
        aj%dxdxl = 0.0 ! So it's not undefined, because the user isn't expected
                       ! to set it when starting.
        call do_NF_Evalf
      else
        ! User forcing retreat to best X followed by gradient move
        call Do_NF_GMove_First
      end if

    else if ( nflag /= aj%where_to ) then

      call ermsg ( me, 99, 2, &
        & 'NFLAG was changed to be inconsistent with AJ%WHERE_TO', '.' )

    else

      ! Continuing after return for reverse communication for the caller
      ! to do some computation:

      select case ( aj%where_to )        ! Reason for returning
      case ( nf_evalf )                  ! Resumed after evaluating F

        ! Re-enter after evaluating F.

        iter = iter + 1

        ! Test if a retreat to a previous -- best -- X should be made

        tp = 0.9_rk * aj%spl * aj%ajn     ! .9 times sql on previous iteration
        if ( aj%fnorm >= aj%fnorml ) then ! Things are not getting better

        ! AJ%FNORML will have been set to zero if there are four consecutive
        ! gradient moves.
        ! Test for X convergence

          if ( x_converge() ) return
          if ( inc >= 0 ) then   ! We are not starting
            if ( inc == 0 ) then ! Last X == Best X
              aj%dxnbig = max(aj%dxnl,dxnois)

            else  ! inc > 0 -- Not starting, Last X not Best X
              if ( kb == 0 ) then
                if ( aj%fnorm < aj%fnormb ) then ! New best |F|
                  aj%gfac = -aj%gfac
                  if ( aj%gfac >= 0 ) then
                    ! Error processing
                    call ermsg ( me, 2, 0, 'J or F may be in error', '.' )
                    nflag = nf_fandj
                    aj%where_to = nflag
                    return
                  end if
                end if
                call Do_NF_GMove_Small
                return
              end if
              if ( aj%dxnl <= dxnois ) then ! Last Newton move tiny?
                call do_NF_EvalJ ( tp )
                return
              end if
              if ( kb < 0 ) then ! Gradient move last time
                call do_NF_Gmove_second
                return
              end if
! WVS revised this 2012-08-27 based upon advice from FTK
!             if ( (tp >= min(aj%sql,aj%sqb+aj%sqb)) .or. (inc >= incbig) ) then
              if ( aj%fnmin > (1.0+0.5*relsf)*aj%fnminb .or. &
                 & inc >= incbig ) then
                call do_NF_GMove_First
                return
              end if
            end if
            aj%sqmin = min(aj%sqb, max(aj%spl, aj%spact) * 4.0 * aj%ajn)
            aj%dxinc = max(aj%dxnl * 0.5,dxnois)
            inc = inc + 1
          end if
        else if ( aj%fnorm >= aj%fnormb ) then ! Things are getting better, but
                                               ! no better than at the best X
          ! If our Marquardt parameter is small and AJ%FNORM does not appear to
          ! be getting close to or smaller than AJ%FNORMB, and we are not
          ! making a move in a quite different direction we give up and go back
          ! to the previous best X and take a gradient move.

          if ( (max(tp,aj%sql) <= aj%sqmin) .and. &
             & ((aj%fnorm*(aj%fnorm/aj%fnorml)**2) > aj%fnormb) .and. &
             & (aj%cdxdxl >= 0.25_rk) ) then
            call Do_NF_GMove_First
            return
          end if
        end if

        call do_NF_EvalJ ( tp )

      case ( nf_evalj )              ! Resumed after evaluating J

        ! Re-enter after evaluating Jacobian

        ! diag =  Smallest | diagonal element | after factoring
        ! ajn =   Largest L1 norm of column in upper triangle
        ! fnmin = L2 Norm of F not in column space of the Jacobian
        ! gradn = L2 norm of gradient
        if ( aj%gradn <= 0.0 ) then
          aj%gradn = 1.0
          if ( aj%ajn <= 0.0 ) aj%ajn = 1.0
        end if
        tp = aj%diag/aj%ajn
        aj%condai = max(min(aj%condai,tp),tp**2)
        aj%spact = 0.25 * aj%condai
        aj%frz = sqrt(abs((aj%fnorm-aj%fnmin)*(aj%fnorm+aj%fnmin)))

        if ( iter == 1 ) then ! First iteration
          aj%sp = spstrt
          aj%spinc = aj%sp
          if ( aj%frz <= rnd*aj%fnmin ) then
            call do_nf_tolf
            return
          end if
          aj%dxi = 0.125_rk
        else
          if ( aj%fnorm >= aj%fnorml ) then  ! AJ%FNORM increased

            ! Select new stabilizing parameter for case when AJ%FNORM has
            ! increased

            kb = 1
            if ( aj%dxnl < dxnois ) then
              aj%sp = 1.0E-4_rk * aj%sp
              aj%sqmin = min(aj%sqmin, aj%sp * aj%ajn)
            ! aj%dxnl = aj%dxnl * 100.0
              call do_NF_Solve ( before_Solve = set_new_marquardt )
            end if
            aj%spfac = max( (aj%fnorm / aj%fnorml)**2, 10.0_rk )
            aj%spinc = aj%spl + aj%spl + aj%spact
            aj%dxi = 1.0
            call do_NF_Solve ( before_Solve = test_gradn_increase )
          end if
          if ( aj%fnorm >= aj%fnormb ) then ! Is function norm larger than best?
            if ( .not. f_converge() ) &
              & call do_NF_Solve ( before_Solve = test_f_linear )
            return
          end if
        end if

        ! New best X

        if ( aj%kfail /= 0 ) then
          if ( kb /=0 ) then ! Got new best immediately
            aj%dxfail = min(1.0625_rk*aj%dxfail,aj%dxbad)
          else if ( aj%dxn >= 0.9_rk*aj%dxbad ) then
            aj%dxbad = aj%dxbad*(1.25_rk**aj%kfail)
            aj%dxfail = aj%dxbad
            aj%kfail = aj%kfail + 1
          else
            aj%dxfail = 0.75_rk*aj%dxfail + 0.25_rk*aj%dxbad
            aj%kfail = 1
          end if
          aj%dxinc = min(dxmaxi,aj%dxfail)
        else
          aj%dxinc = dxmaxi
        end if

        aj%spinc = min(aj%spinc, max(spmini, abs(ajscal) / aj%ajn))
        aj%sqmin = 0.0
        aj%spb = aj%spl
        aj%sqb = aj%sql
        if ( inc < 0 ) then ! Starting
          call do_NF_Solve ( before_Solve = set_new_marquardt )
          return
        end if
        if ( inc /= 0 ) then
          if ( kb < 0 ) then
            aj%dxi = 0.25
          else if ( kb == 1 ) then
          ! We had a function increase prior to getting a new best.
            aj%dxi = 0.25 * aj%dxnl / aj%dxnbig
          end if
          kb = 0
          inc = 0
        end if
        if ( .not. f_converge() ) &
          & call do_NF_Solve ( before_Solve = test_f_linear )

      case ( nf_solve )              ! Resumed after solving for DX

      ! Re-enter after doing Levenberg-Marquardt stabilization and solving for
      ! candidate DX

        ! fnmin = L2 Norm of F not in column space of the Jacobian
        aj%cgdx = aj%gdx/(aj%gradn*aj%dxn) ! Cosine ( gradient, dx )
        aj%condai = min(aj%cgdx,aj%gradn/(aj%dxn*aj%ajn**2),aj%diag/aj%ajn)
        aj%fnxe = aj%fnmin**2 - (aj%sq*aj%dxn)**2
        if ( inc >= 0 ) then
          aj%cdxdxl = aj%dxdxl/(aj%dxn*aj%dxnl) ! cosine between DX and last DX
          if ( aj%fnxe > aj%fnormb**2 ) then
            if ( k1it /= 0 ) call nwtdb ( width=9, level=0, why='Give up' )
            call Do_NF_GMove_First
            return
          end if
          if ( aj%dxn < 1.1_rk*aj%dxinc ) then
            if ( aj%sq <= 0.0 .or. aj%dxn > 0.90_rk*aj%dxinc ) then
              if ( inc==0 ) then
                if ( .not. x_converge() ) call do_nf_best
                return
              end if
              aj%cait = cbig
              call do_NF_DX
              return
            end if
!! WVS: was          dxinc = aj%dxi * aj%dxn * (1.0 - aj%cdxdxl)**2
!! This was clearly wrong, since it says "If consecutive Newton moves are in
!! the same direction, take shorter steps"
            aj%dxinc = aj%dxi * aj%dxn / (1.025 - aj%cdxdxl)**2
          end if
        else if ( aj%dxn <= aj%dxinc ) then ! Move length OK?
          inc = incbig
          aj%cdxdxl = 0.0 ! Don't think about Aitken
          call do_nf_best
          return
        end if

        ! Newton step length is too large or too small

        if ( k1it /= 0 ) call nwtdb ( width=9, level=0, why='Step length' )

        !{ Calculate {\tt qn}$^2 = {\bf q}^T{\bf q} = ||{\bf q}||^2$, used
        !  to determine the Levenberg-Marquardt parameter $\lambda$, where
        !  ${\bf q} = {\bf U}^{-T} {\delta\bf x}$ and ${\delta\bf x}$ is the
        !  Newton step computed with the current value of $\lambda$.  The
        !  method to determine $\lambda$ and the reason for being interested
        !  in ${\bf q}^T{\bf q}$ are explained below.

        nflag = nf_lev
        aj%where_to = nflag

      case ( nf_lev )              ! Resumed after doing the calculations for
                                   ! the Moré-Sorensen algorithm

!{ Determine the Levenberg-Marquardt parameter $\lambda$, known here as AJ%SQ.
!  The last Marquardt parameter used was AJ%SQ, which means that AJ%SQ**2 was
!  added to the diagonal of the normal equations.  A little bit of
!  derivation so you have some idea of what is going on.  Start with
!
!  $({\bf H} + \lambda {\bf I}) {\delta\bf x} = -{\bf g}$,
!
!  where ${\bf H}$ is the Hessian matrix.  In our case, ${\bf J}^T {\bf J}$,
!  where ${\bf J}$ is the Jacobian matrix, is an approximate Hessian.
!  Differentiating both sides with respect to $\lambda$ gives
!
!  $({\bf H} + \lambda {\bf I}) \frac{\text{d}{\delta\bf x}}{\text{d}\lambda}
!  + {\bf I}\, {\delta\bf x} = 0$ or
!  $({\bf H} + \lambda {\bf I}) \frac{\text{d}{\delta\bf x}}{\text{d}\lambda}
!  = -{\delta\bf x}$.
!
!  Cholesky factoring $({\bf H} + \lambda {\bf I})$, we have
!  ${\bf U}^T {\bf U} \frac{\text{d}{\delta\bf x}}{\text{d}\lambda} =
!  -{\delta\bf x}$ or ${\bf U} \frac{\text{d}{\delta\bf x}}{\text{d}\lambda} =
!  -{\bf U}^{-T} {\delta\bf x} = -{\bf q}$.
!
!  Multiplying ${\bf U}^T {\bf U} \frac{\text{d}{\delta\bf x}}{\text{d}\lambda} =
!  -{\delta\bf x}$ on the left by $\frac{\text{d}{\delta\bf x}}{\text{d}\lambda}^T$
!  and substituting
!  ${\bf U} \frac{\text{d}{\delta\bf x}}{\text{d}\lambda} = - {\bf q}$ we have
!  $\frac{\text{d}{\delta\bf x}}{\text{d}\lambda}^T {\bf U}^T {\bf U}
!   \frac{\text{d}{\delta\bf x}}{\text{d}\lambda}
!  = -\frac{\text{d}{\delta\bf x}}{\text{d}\lambda}^T {\delta\bf x} =
!  {\bf q}^T {\bf q}$.
!
!  We intend to find a zero of
!  $\phi = \frac1{|| {\delta\bf x} ||} - \frac1{\rho}$, where $\rho$ is the
!  desired step length for ${\delta\bf x}$,
!  treating ${\delta\bf x}$ as a function of $\lambda$.
!
!  $\frac{\text{d}\phi}{\text{d}\lambda} =
!  \frac{\partial ( (|| {\delta\bf x} ||^2)^{-1/2}
!   - \frac1{\rho})}{\partial{\delta\bf x}}
!  \frac{\text{d}{\delta\bf x}}{\text{d}\lambda}
!  = \frac{-1}{|| {\delta\bf x} ||^{3}} {\delta\bf x}^T
!    \frac{\text{d}{\delta\bf x}}{\text{d}\lambda}
!  = \frac{-1}{|| {\delta\bf x} ||^{3}} (- {\bf q}^T {\bf q})$
!
!  The Newton method gives us
!
!  $\lambda_{\text{new}} = \lambda_{\text{old}} -
!    \phi / \frac{\text{d}\phi}{\text{d}\lambda_{\text{old}}}
!    = \lambda_{\text{old}} - \frac{|| {\delta\bf x} ||^3}{{\bf q}^T {\bf q}}
!      \left( \frac1{|| {\delta\bf x} ||} - \frac1{\rho} \right)
!    = \lambda_{\text{old}} - \frac{|| {\delta\bf x} ||^2}{|| {\bf q} ||^2}
!      \left(1 - \frac{|| {\bf \delta x} ||}{\rho}\right)$
!
! If you have been looking at Mor\'e and Sorensen, this update is the
! negative of theirs since they define $\lambda$ as the negative of the
! $\lambda$ above.  Also, Mor\'e and Sorensen write $\lambda$ as the quantity
! added to the diagonal of the normal equations, which we call $\lambda^2$, so
! the iteration above and the code below are really iterating on what we call
! $\lambda^2$.  I.e., from a least-squares point of view we write
!%
! \begin{equation*}
! \left[ \begin{array}{c} {\bf J} \\
!                 \lambda {\bf I} \\
! \end{array} \right] {\delta \bf x} \simeq
! \left[ \begin{array}{c} -{\bf f} \\
!                          {\bf 0} \\
! \end{array} \right]
! \end{equation*}
! which, in normal-equations form is $({\bf J}^T {\bf J} + \lambda^2 {\bf I})
! {\delta\bf x} = -{\bf J}^T {\bf f} = -{\bf g}$.

        aj%sq = sqrt(max(0.0_rk, &
                         aj%sq**2-(aj%dxn**2/aj%qnsq)*(1.0 - aj%dxn/aj%dxinc)))
        aj%sp = aj%sq / aj%ajn
        ! Solve again with new Marquardt parameter
        call do_NF_Solve ( before_Solve = solve_only )

      case ( nf_newx )               ! Resumed after calculating new X

        ! Re-enter after computing new X

        if ( .not. aj%big ) then ! All dx's are very small
          if ( (aj%spl <= spmini) .or. (aj%sq == 0.0) ) then

            if ( (inc > 0) .and. (kb /= 0) ) then ! Go do gradient move
              call Do_NF_GMove_First
              return
            end if

            ! Convergence -- move too small

            nflag = nf_too_small
            aj%where_to = nflag
            return
          else
            dxnois = max(dxnois,10.0_rk*aj%dxn) ! Set so steps don't get too small.
          end if
        end if

        call do_NF_Evalf

      case ( nf_gmove )              ! Resumed after gradient move from best X

        ! Re-enter after computing gradient step from best X

        aj%dxn = aj%gradnb * aj%gfac ! In case the caller wants to compute cosines
        if ( aj%gfac <= 0.0 ) then
          call do_NF_NewX ( .true. )
          return
        end if
        aj%axmax = aj%axmaxb
        aj%gfac = 0.125 * aj%gfac
        aj%spl = aj%spg
        aj%spg = 8.0 * aj%spl
        aj%fnxe = 0.5 * aj%fnorm**2
        aj%frz = aj%frzb
        aj%gradn = aj%gradnb
        call do_NF_NewX

      case ( nf_best )               ! Resumed after saving X as best X

        ! Re-enter here after saving X as "best X" and Gradient as
        ! "Best gradient"

        aj%axmaxb = aj%axmax
        dxnois = max(dxnois,1.0E-4_rk*aj%dxnl)
        aj%gfac = min(0.125*aj%dxn/aj%gradn,aj%ajn**(-2))
        aj%spg = aj%spl + 0.01
        aj%fnormb = aj%fnorm
        aj%frzb = aj%frz
        aj%fnminb = aj%fnmin
        aj%gradnb = aj%gradn
        if ( abs(aj%cdxdxl) < 0.9_rk .or. aj%sq /= 0.0 .or. &
           & aj%dxn >= aj%dxnl ) then
          call do_NF_DX
        else

          ! Compute scalar factor for (modified) Aitken acceleration

          nflag = nf_aitken
          aj%where_to = nflag
        end if

      case ( nf_aitken )             ! Resumed after computing some numbers
                                     ! needed for (modified) Aitken
                                     ! acceleration

        tp = aj%dxdx     ! | dX | **2
        if ( tp /= 0.0 ) then
          tp = 1.0 + aj%dxdxl/tp
          tp = min(10.0_rk,max(tp,0.01_rk))
          if ( abs(aj%cait-tp) < abs(0.25_rk*(aj%cait+tp)-0.5_rk) ) then
          ! Set DX = Aitken-modified DX = AJ%CAIT * "Candidate DX"
            aj%cait = tp
            itken = iter + 2
            aj%sq = -aj%cait
            aj%dxn = aj%cait*aj%dxn
            nflag = nf_dx_aitken
            aj%where_to = nflag
            return
          end if
          if ( itken < iter ) aj%cait = tp
        end if
        call do_NF_DX

        ! End of logic for Aitken acceleration

      case ( nf_dx, nf_dx_aitken )   ! Resumed after DX = "Candidate DX" or
                                     ! DX = AJ%CAIT * "Candidate DX"

        aj%dxnl = aj%dxn
        call do_NF_NewX

      ! Continuing after return that was expected to be final:
      case ( nf_tolx, nf_tolx_best ) ! Converged because distance from C to
                                     ! solution was estimated to be small and
                                     ! DX is small
        call do_nf_best
      case ( nf_tolf )               ! Converged because |F| is small
        call do_NF_Solve ( before_Solve = test_f_linear )
      case ( nf_too_small )          ! Converged because DX was too small
        call do_NF_GMove
      case ( nf_fandj )              ! F and J appeared not to be consistent
        call Do_NF_GMove_Small
      end select

    end if

  contains

    ! ...............................................  Do_NF_Best  .....
    subroutine Do_NF_Best

    ! Store best X, the gradient, and other constants used if a
    ! return needs to be made to the best X

      aj%fnormb = aj%fnorm
      nflag = nf_best
      aj%where_to = nflag 
    end subroutine Do_NF_Best

    ! .................................................  Do_NF_DX  .....
    subroutine Do_NF_DX

    ! Store Candidate DX as DX

      nflag = nf_dx
      aj%where_to = nflag
    end subroutine Do_NF_DX

    ! ..............................................  Do_NF_EvalF  .....
    subroutine Do_NF_EvalF

    ! Evaluate F

      nflag = nf_evalf
      aj%where_to = nflag
    end subroutine Do_NF_EvalF

    ! ..............................................  Do_NF_EvalJ  .....
    subroutine Do_NF_EvalJ ( New_SQL )
      real(rk), intent(in) :: New_SQL

    ! Evaluate J

      aj%sql = new_sql
      nflag = nf_evalj
      aj%where_to = nflag
    end subroutine Do_NF_EvalJ

    ! ..............................................  Do_NF_GMove  .....
    subroutine Do_NF_GMove

    ! Retreat to a previous -- best -- X and step in gradient direction.  An
    ! alternative is to use the gradients at the best X and the current X to
    ! compute the directional derivatives along the line from the best X to the
    ! current X, fit a cubic polynomial to F using Hermite interpolation, and
    ! go to the minimum of that polynomial if it's between the current X and
    ! the best X.

      nflag = nf_gmove
      aj%where_to = nflag
    end subroutine Do_NF_GMove

    ! ........................................  Do_NF_GMove_First  .....
    subroutine Do_NF_GMove_First

    ! Set up for the first gradient move after a successful Newton move

      kb = -1 ! Set consecutive gradient move counter
      aj%dxbad = 0.75_rk*min(aj%dxn,aj%dxbad)
      aj%dxfail = 0.125_rk*aj%dxbad
      aj%kfail = 1
      call do_NF_GMove_Second
    end subroutine Do_NF_GMove_First

    ! .......................................  Do_NF_GMove_Second  .....
    subroutine Do_NF_GMove_Second

    ! Set up for a gradient move after a previous gradient move, or
    ! finish setting up for a first gradient move.

      kb = kb - 1
      if ( kb <= -4 ) then   ! Four consecutive gradient moves !

        ! Setup to test if Jacobian matrix is being computed properly

        kb = 0
        aj%fnorml = 0.0
        aj%gfac = -100.0 * aj%gfac
      end if
      call do_NF_GMove
    end subroutine Do_NF_GMove_Second

    ! ........................................  Do_NF_GMove_Small  .....
    subroutine Do_NF_GMove_Small

    ! Set up for a really small gradient move.

      kb = -1 ! Restart consecutive gradient move counter
      aj%gfac = -0.01 * aj%gfac
      call do_NF_GMove
    end subroutine Do_NF_GMove_Small

    ! ...............................................  Do_NF_NewX  .....
    subroutine Do_NF_NewX ( Gfac_Negative )
      logical, intent(in), optional :: Gfac_Negative ! Present means true

    ! Compute new X and test for too small a correction

      if ( .not. present(gfac_negative) ) then
        aj%fnorml = aj%fnorm
        aj%frzl = aj%frz
        aj%fnxe = 0.25_rk*aj%fnorm**2 + 0.76_rk*aj%fnxe
      end if
      if ( k1it /= 0 ) then ! Do requested debugging output
        k1it = k1it - 1
        call nwtdb ( width=9, level=0, why='Before new X' )
      end if

      aj%starting = inc < 0
      nflag = nf_newx
      aj%where_to = nflag
    end subroutine Do_NF_NewX

    ! ..............................................  Do_NF_Solve  .....

    subroutine Do_NF_Solve ( Before_Solve )
      integer, intent(in) :: Before_Solve

      ! Select new stabilizing parameter.  Also comes back here if user
      ! continues after TOLF convergence.

      if ( before_solve >= test_F_Linear ) then

        tp = min((aj%frz/aj%frzl),(aj%fnorm/aj%fnorml))

        ! Test if F appears almost linear over last step

        if ( aj%fnorm**2 < 1.125 * aj%fnxe ) then ! F appears almost linear.
          aj%spfac = 0.125 * (1.025 - aj%cdxdxl) * tp**2
          if ( aj%spl <= aj%spinc ) aj%spinc = 0.25 * aj%spinc
          aj%dxi = min(aj%dxi,0.25_rk)
        else ! F not linear over last step
        ! aj%spfac = tp*(aj%fnorm**2)/aj%fnxe
        ! if ( aj%cdxdxl >= 0.9_rk) aj%spfac = min(aj%spfac,0.5_rk)
        ! On 2002/07/25, FTK recommended:
          aj%spfac = min( tp*(aj%fnorm**2)/aj%fnxe, &
                        & 32.0 * (1.025 - aj%cdxdxl)**2 )
          aj%dxi = 1.0
        end if
      end if
      if ( before_solve >= test_GradN_Increase ) then
        if ( aj%gradn > aj%gradnl ) aj%spfac = aj%spfac * aj%gradnl / aj%gradn
        aj%sp = max(aj%spinc, min(aj%spb, aj%spl * aj%spfac))
      end if
      if ( before_solve >= set_New_Marquardt ) then
        aj%sp = max(aj%sp,aj%sqmin / aj%ajn)
    !   if ( inc == 0 ) aj%sp = min(aj%sp, 0.5 * aj%spl)
        aj%spl = aj%sp
        aj%sq = aj%sp * aj%ajn
        aj%sqt = aj%sq
        if ( aj%sp < aj%spact .and. inc >= 0 ) aj%sq = aj%sqmin
      end if

    ! Do Levenberg-Marquardt stabilization, solve for "Candidate DX"

      aj%starting = inc<0
      nflag = nf_solve
      aj%where_to = nflag
    end subroutine Do_NF_Solve

    ! ...............................................  Do_NF_Tolf  .....

    subroutine Do_NF_Tolf

    ! Announce function (not DX) convergence

      nflag = nf_tolf
      aj%where_to = nflag
    end subroutine Do_NF_Tolf

    ! ...............................................  F_Converge  .....
    logical function F_Converge ( )

    ! Test for convergence in sense of AJ%FNORM < (1 + RELSF)*AJ%FNMIN

      f_converge = .false.
      if ( inc==0 ) then
        tp = aj%fnorm - (1.0+relsf)*aj%fnmin
        if ( tp < 0.0 ) then
                    ! WVS changed this 2012-09-10 to allow convergence
                    ! with nonzero Levenberg-Marquardt parameter
!           if ( (aj%sq <= aj%sqmin) .or. (aj%spl <= spmini) .or. &
!            &   (tp <= -(1.0+relsf)*aj%spl*aj%fnorm) ) then
          call do_nf_tolf
          f_converge = .true.
        end if
      end if

    end function F_Converge

    ! ...............................................  X_Converge  .....
    logical function X_Converge ( )

    ! Test whether convergence has occurred due to very small \delta X

      x_converge = .false.
      if ( aj%fnorml > 0.0 ) then
      ! Didn't just do four gradient moves from the best X
        if ( aj%dxn <= tolxa .or. &
           & ( (aj%sq == 0.0) .or. (aj%spl <= spmini) ) .and. &
           &   aj%dxn <= tolxr*aj%axmax ) then
        ! Convergence if aj%dxn is small enough or we have a suffiently small
        ! Levenberg-Marquardt parameter and we had a very small move.
          nflag = nf_tolx
          if ( kb /= 0 ) then  ! Made gradient move, or F increased
            aj%fnorm = aj%fnormb
            nflag = nf_tolx_best ! Solution is at best X, not current X
          end if
          aj%where_to = nflag
          x_converge = .true.
        end if
      end if

    end function X_Converge

  end subroutine DNWTA

! *************************************************     DNWTOP     *****
  subroutine DNWTOP ( AJ, NOPT, XOPT )
    ! Process option vector for DNWT.  With no arguments, does default
    ! initialization.
    type(nwt_t), intent(inout) :: AJ
    integer, intent(in), optional :: NOPT(*)
    real(rk), intent(in), optional :: XOPT(*)

    logical, save :: FIRST = .true.
    integer I, INDIC, K, KA, NACT
    integer, save :: IOPTS(8) = (/ 0, 0, 0, 0, 0, 0, 0, 150 /)
    character(len=4), parameter :: LABL(2) = (/ '(I) ', 'NOPT' /)
    ! ??? Defaults aren't what comments say
    real(rk),save :: VALUES(9) = (/ 0.1_rk, epsilon(values), epsilon(values), &
                                  & 0.1*huge(values), 0.01_rk, 0.0_rk, 0.0_rk, &
                                  & 0.0_rk, 0.0_rk /)
    real(rk),save :: VALNOM(9) = (/ 0.1_rk, epsilon(values), epsilon(values), &
                                  & 0.1*huge(values), 0.01_rk, 0.0_rk, 0.0_rk, &
                                  & 0.0_rk, 0.0_rk /)
    integer :: IVAL(2) ! for error messages

!*************** Start of executable code ******************

    if ( first ) then
      valnom(7) = sqrt(sqrt(tiny(values)))**3
      valnom(8) = sqrt(sqrt(epsilon(values)))**3
      values(7:8) = valnom(7:8)
      first = .false.
    end if
    if ( present(nopt) ) then
      if ( .not.present(xopt) ) then
        call ermsg ( me, 99, 2, 'NOPT present but XOPT absent', '.' )
        return
      end if
      i = 1
      do while ( nopt(i) /= 0 )
        k = nopt(i)
        ka = abs(k)
        select case (ka)
!****************** Change *DATA* values *******************
        case ( 1, 2 )      ! Set K1IT
          if ( k > 0 ) iopts(ka) = nopt(i+1)
          k1it = max(iopts(1),iopts(2))
        case ( 3, 4, 7, 8) ! Change XSCAL and FSCAL indices in data
                           ! or set up bounds option in data
          if ( k > 0 ) iopts(ka) = nopt(i+1)
        case ( 5, 6 )      ! Set flags in data for reverse communi-
                           ! cation and special matrix operations
          if ( k > 0 ) iopts(ka) = 1
          i = i - 1
        case ( 9, 10 )     ! Set to default values
          iopts = 0
          iopts(8) = 150
          values = valnom
          i = i - 1
        case ( 11: 19 )    ! Set in data SPSTRT, SPMINI, AJSCAL,
                           ! DXMAXI, RELSF, and DXNOIS
          values(ka-10) = valnom(ka-10) ! Reset to nominal value
          if ( k > 0 ) & ! If indicated, set to user input value
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
      aj%where_to = nf_start
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

  subroutine DNWTDB ( AJ, WIDTH, LEVEL, WHY )

!   Print the scalars in the module.  Print the stuff in AJ if it's
!   present.  Print the integer scalars first.  Then print
!   max(5,min(9,WIDTH)) real scalars per line if width is present, else
!   print five per line.

    type (NWT_T), intent(in), optional :: AJ
    integer, intent(in), optional :: WIDTH
    integer, intent(in), optional :: LEVEL ! Absent, do everything
    ! Present and 1 (default) do everything
    ! Present and 0 don't do internal variables,
    ! and if AJ is present skip AXMAX, AXMAXB, CAIT, DXI, DXNL, FNMINB,
    !   FNORMB, FNORML, FRZB, FRZL, GRADNB, GRADNL, from AJ
    character(len=*), intent(in), optional :: WHY ! printed if present

!   namelist /DNWTDB_OUT/ AJSCAL
!   namelist /DNWTDB_OUT/ DXMAXI, DXNOIS
!   namelist /DNWTDB_OUT/ INC, ITER, ITKEN
!   namelist /DNWTDB_OUT/ K1IT, KB
!   namelist /DNWTDB_OUT/ RELSF, SPMINI
!   namelist /DNWTDB_OUT/ SPSTRT

    integer :: I                        ! Which one in the line is being worked
    character(len=9) :: Where_To_Name
    integer :: MyLevel, MyWidth
    character(len=132) :: Name_Line     ! For names
    character(len=132) :: Output_Line   ! For values

    interface Add_To_Line
      procedure Add_To_Line_I, Add_To_Line_L, Add_To_Line_R
    end interface

    myLevel = 1
    if ( present(level) ) myLevel = level
    myWidth = 5
    if ( present(width) ) myWidth = max(5,min(9,width))
    if ( myLevel == 0 ) myWidth = min(8,myWidth)

    output_line = '-----     DNWT internal variables     ----------------------&
    &------------------------------------------------------------------'
    call output ( output_line(1:14*myWidth), advance='yes' )

    call flagName ( aj%where_to, where_to_name )

    call output (  &
      & '  WHERE_TO        INC    ITER   ITKEN    K1IT      KB' )
    if ( present(why) ) call output ( '  WHY' )
    call newLine
    write ( output_line, '(1x,a9,i11,4i8,1x,a9)' ) adjustr(where_to_name), INC, &
         & ITER, ITKEN, K1IT, KB
    call output ( trim(output_line) )
    if ( present(why) ) call output ( '  '//trim(why) )
    call newLine

    i = 1
    name_line = ''
    output_line = ''
    if ( myLevel > 0 ) call add_to_line ( ajscal, 'AJSCAL' )
    if ( myLevel > 0 ) call add_to_line ( dxmaxi, 'DXMAXI' )
    if ( myLevel > 0 ) call add_to_line ( dxnois, 'DXNOIS' )
    if ( myLevel > 0 ) call add_to_line ( relsf,  'RELSF' )
    if ( myLevel > 0 ) call add_to_line ( spmini, 'SPMINI' )
    if ( myLevel > 0 ) call add_to_line ( spstrt, 'SPSTRT' )
    if ( myLevel > 0 ) call add_to_line ( tolxa,  'TOLXA' )
    if ( myLevel > 0 ) call add_to_line ( tolxr,  'TOLXR' )
    if ( i/=1 ) call print_lines

    if ( present(aj) ) then
      output_line = '-----     DNWT external variables     ----------------------&
        &------------------------------------------------------------------'
      call output ( output_line(1:14*myWidth), advance='yes' )
      output_line = ''

      call add_to_line ( aj%ajn,        'AJN' )
      if ( myLevel > 0 ) call add_to_line ( aj%axmax,      'AXMAX' )
      if ( myLevel > 0 ) call add_to_line ( aj%axmaxb,     'AXMAXB' )
      if ( myLevel > 0 ) call add_to_line ( aj%cait,       'CAIT' )
      call add_to_line ( aj%cdxdxl,     'CDXDXL' )
      call add_to_line ( aj%cgdx,       'CGDX' )
      call add_to_line ( aj%condai,     'CONDAI' )
      call add_to_line ( aj%diag,       'DIAG' )
      call add_to_line ( aj%dxbad,      'DXBAD' )
      if ( myLevel > 0 ) call add_to_line ( aj%dxi,        'DXI' )
      call add_to_line ( aj%dxinc,      'DXINC' )
      if ( myLevel > 0 ) call add_to_line ( aj%dxnbig,     'DXNBIG' )
      call add_to_line ( aj%dxfail,     'DXFAIL' )
      call add_to_line ( aj%dxdx,       'DXDX' )
      call add_to_line ( aj%dxdxl,      'DXDXL' )
      call add_to_line ( aj%dxn,        'DXN' )
      if ( myLevel > 0 ) call add_to_line ( aj%dxnl,       'DXNL' )
      call add_to_line ( aj%fnmin,      'FNMIN' )
      if ( myLevel > 0 ) call add_to_line ( aj%fnminb,     'FNMINB' )
      call add_to_line ( aj%fnorm,      'FNORM' )
      if ( myLevel > 0 ) call add_to_line ( aj%fnormb,     'FNORMB' )
      if ( myLevel > 0 ) call add_to_line ( aj%fnorml,     'FNORML' )
      call add_to_line ( sqrt(aj%fnxe), 'FNXE**.5')
      call add_to_line ( aj%frz,        'FRZ' )
      if ( myLevel > 0 ) call add_to_line ( aj%frzb,       'FRZB' )
      if ( myLevel > 0 ) call add_to_line ( aj%frzl,       'FRZL' )
      call add_to_line ( aj%gdx,        'GDX' )
      if ( myLevel > 0 ) call add_to_line ( aj%gfac,       'GFAC' )
      call add_to_line ( aj%gradn,      'GRADN' )
      if ( myLevel > 0 ) call add_to_line ( aj%gradnb,     'GRADNB' )
      if ( myLevel > 0 ) call add_to_line ( aj%gradnl,     'GRADNL' )
      call add_to_line ( aj%kfail,      'KFAIL' )
      call add_to_line ( aj%qnsq,       'QNSQ' )
      call add_to_line ( aj%sp,         'SP ' )
      if ( myLevel > 0 ) call add_to_line ( aj%spact,      'SPACT' )
      if ( myLevel > 0 ) call add_to_line ( aj%spb,        'SPB' )
      if ( myLevel > 0 ) call add_to_line ( aj%spfac,      'SPFAC' )
      if ( myLevel > 0 ) call add_to_line ( aj%spg,        'SPG' )
      if ( myLevel > 0 ) call add_to_line ( aj%spinc,      'SPINC' )
      call add_to_line ( aj%spl,        'SPL' )
      call add_to_line ( aj%sq,         'SQ' )
      if ( myLevel > 0 ) call add_to_line ( aj%sqb,        'SQB' )
      if ( myLevel > 0 ) call add_to_line ( aj%sql,        'SQL' )
      call add_to_line ( aj%sqmin,      'SQMIN' )
      call add_to_line ( aj%sqt,        'SQT' )
      call add_to_line ( aj%kfail,      'KFAIL' )
      call add_to_line ( aj%big,        'BIG' )
      call add_to_line ( aj%starting,   'STARTING' )
      if ( i /=1 ) call print_lines

    end if
    output_line = '------------------------------------------------------------&
      &------------------------------------------------------------------'
    call output ( output_line(1:14*myWidth), advance='yes' )

  contains
    subroutine Add_To_Line_I ( Value, Name )
      integer, intent(in) :: Value
      character(len=*), intent(in) :: Name
      name_line(1+14*(i-1):10+14*(i-1)) = name
      name_line(1+14*(i-1):10+14*(i-1)) &
        & = adjustr(name_line(1+14*(i-1):10+14*(i-1)))
      write ( output_line(1+14*(i-1):10+14*(i-1)), '(I10)' ) value
      i = i + 1
      if ( i > myWidth ) call print_lines
    end subroutine Add_To_Line_I
    subroutine Add_To_Line_L ( Value, Name )
      logical, intent(in) :: Value
      character(len=*), intent(in) :: Name
      name_line(1+14*(i-1):10+14*(i-1)) = name
      name_line(1+14*(i-1):10+14*(i-1)) &
        & = adjustr(name_line(1+14*(i-1):10+14*(i-1)))
      write ( output_line(1+14*(i-1):1+14*(i-1)), '(L1)' ) value
      output_line(1+14*(i-1):10+14*(i-1)) &
        & = adjustr(output_line(1+14*(i-1):10+14*(i-1)))
      i = i + 1
      if ( i > myWidth ) call print_lines
    end subroutine Add_To_Line_L
    subroutine Add_To_Line_R ( Value, Name )
      real(rk), intent(in) :: Value
      character(len=*), intent(in) :: Name
      name_line(1+14*(i-1):10+14*(i-1)) = name
      name_line(1+14*(i-1):10+14*(i-1)) &
        & = adjustr(name_line(1+14*(i-1):10+14*(i-1)))
      write ( output_line(1+14*(i-1):14*i), '(es14.7)' ) value
      i = i + 1
      if ( i > myWidth ) call print_lines
    end subroutine Add_To_Line_R
    subroutine Print_Lines
      call output(trim(name_line),advance='yes')
      call output(trim(output_line),advance='yes')
      i = 1
      name_line = ''
      output_line = ''
    end subroutine Print_Lines
  end subroutine DNWTDB

! **********************************************     DNWT_GUTS     *****

  subroutine DNWT_GUTS ( GUTS )
  ! Get the private saved variables of DNWT, for debugging purposes
    type(NWT_GUTS), intent(out) :: GUTS

    guts%ajscal = ajscal
    guts%dxmaxi = dxmaxi
    guts%dxnois = dxnois
    guts%inc = inc
    guts%iter = iter
    guts%itken = itken
    guts%k1it = k1it
    guts%kb = kb
    guts%relsf = relsf
    guts%spmini = spmini
    guts%spstrt = spstrt
    guts%tolxa = tolxa
    guts%tolxr = tolxr
  end subroutine DNWT_GUTS

! ***********************************************     FlagName     *****

  subroutine FlagName ( NFlag, ItsName )
  ! Return the name of *NWTA's flag
    integer, intent(in) :: NFlag
    character(len=*), intent(out) :: ItsName
    select case (nflag)
    case (NF_EVALF)
      itsName = 'EVALF'
    case (NF_EVALJ)
      itsName = 'EVALJ'
    case (NF_LEV)
      itsName = 'LEV'
    case (NF_SOLVE)
      itsName = 'SOLVE'
    case (NF_NEWX)
      itsName = 'NEWX'
    case (NF_GMOVE)
      itsName = 'GMOVE'
    case (NF_BEST)
      itsName = 'BEST'
    case (NF_AITKEN)
      itsName = 'AITKEN'
    case (NF_DX)
      itsName = 'DX'
    case (NF_DX_AITKEN)
      itsName = 'DX_AITKEN'
    case (NF_START)
      itsName = 'START'
    case (NF_TOLX)
      itsName = 'TOLX'
    case (NF_TOLX_BEST)
      itsName = 'TOLX_BEST'
    case (NF_TOLF)
      itsName = 'TOLF'
    case (NF_TOO_SMALL)
      itsName = 'TOO_SMALL'
    case (NF_FANDJ)
      itsName = 'FANDJ'
    case default
      itsName = 'What???'
    end select
  end subroutine FLAGNAME

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
    character (len=*), parameter :: IdParm = &
       "$Id$"
    character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1)==ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function NOT_USED_HERE
!---------------------------------------------------------------------------

end module DNWT_Module

! $Log$
! Revision 2.56  2018/07/23 23:29:56  vsnyder
! Remove ChiSqMinNorm and ChiSqNorm from DNWT AJ structure because DNWT
! doesn't compute them or use them.
!
! Revision 2.53  2015/05/15 23:41:01  vsnyder
! Add FNORMB to AJ
!
! Revision 2.52  2012/09/13 18:06:56  vsnyder
! Add SQMIN to AJ structure, allow convergence with nonzero lambda
!
! Revision 2.51  2012/08/30 23:04:02  vsnyder
! Use FNMIN at best X to decide on gradient move
!
! Revision 2.50  2012/04/20 01:29:21  vsnyder
! Change dxinc calculation, add QNSQ output
!
! Revision 2.49  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.48  2007/09/22 00:25:16  vsnyder
! Initialize DXMAXI to 0.1*huge() so as not to get overflow
!
! Revision 2.47  2007/01/11 20:46:33  vsnyder
! Correct a comment
!
! Revision 2.46  2006/06/08 23:55:59  vsnyder
! Tighten More' Sorensen convergence criteria, TeXnicalities
!
! Revision 2.45  2006/06/06 15:22:51  vsnyder
! Some code restructuring and commenting
!
! Revision 2.44  2006/06/03 00:15:59  vsnyder
! Respect initial Levenberg-Marquardt parameter
!
! Revision 2.43  2006/05/30 22:44:36  vsnyder
! Handle Newton moves after gradient moves better
!
! Revision 2.42  2005/06/22 18:57:01  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.41  2003/02/03 23:11:23  vsnyder
! Implement Mor\'e and Sorensen calculation of the Lavenberg parameter
!
! Revision 2.40  2003/01/18 01:39:57  vsnyder
! More output, prepare for More and Sorensen
!
! Revision 2.39  2003/01/16 21:48:37  vsnyder
! More work on internal output
!
! Revision 2.38  2003/01/16 04:06:58  vsnyder
! Change some output in DNWT
!
! Revision 2.37  2003/01/15 01:49:51  vsnyder
! Add an '12-interesting-variables' dump
!
! Revision 2.36  2002/10/25 22:24:33  livesey
! Changed order of nwt_t
!
! Revision 2.35  2002/10/23 01:15:06  livesey
! Added chiSq stuff
!
! Revision 2.34  2002/10/07 23:43:11  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.33  2002/09/23 22:04:11  vsnyder
! Correct an error message, better computation of defaults for two options
!
! Revision 2.32  2002/09/21 00:33:47  vsnyder
! Repair an error message
!
! Revision 2.31  2002/09/21 00:21:44  vsnyder
! Correct a FlagName result
!
! Revision 2.30  2002/09/19 01:25:53  vsnyder
! More on when to stop going uphill
!
! Revision 2.29  2002/09/14 02:46:00  vsnyder
! Move test for retreat-to-best, some housecleaning
!
! Revision 2.27  2002/09/11 23:41:53  vsnyder
! Correct improved test for retreating to best X
!
! Revision 2.26  2002/09/10 23:52:56  vsnyder
! Improved test for retreating to best X
!
! Revision 2.25  2002/07/26 22:46:40  vsnyder
! More better output in DNWTDB
!
! Revision 2.24  2002/07/26 01:19:18  vsnyder
! Better output in DNWTDB, changed how Levenberg_Marquardt changes, per Fred
!
! Revision 2.23  2002/07/24 20:32:15  vsnyder
! Moved AXMAX, AXMAXB, CGDX, SP and SQ from DNWTA to module scope.
! Added a "get DNWT's guts" routine.
!
! Revision 2.22  2002/07/24 01:08:11  vsnyder
! Made most local variables SAVE.  Mark Filipiak noticed this problem.
!
! Revision 2.21  2002/02/14 21:53:10  vsnyder
! Add parameters for largest and smallest nwt_flag values
!
! Revision 2.20  2002/01/09 00:00:04  pwagner
! Replaced write or print statements with calls to output
!
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
