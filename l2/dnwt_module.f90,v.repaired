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
!>> 2000-03-21 DNWT_MODULE W. V. Snyder Converted to Fortran 90
!>> 2018-07-17 W. V. Snyder Restructured
!--D replaces "?": ?NWT_MODULE, ?NWT_TYPE, ?NWT, ?NWTA, ?NWTDB, ?NWTOP

!****************** Program description ********************************

! These subroutines solve f(x)=0 (or find the least square solution)
! where f is a vector with NF components, and x is a vector with NX
! components.

! Reverse communication is used.  Reverse communication means that a
! return to the user's program is made whenever additional information
! is required. When this return is made NFlag will be set to indicate
! the necessary action.  Negative values indicate that the calling
! program unit should carry out a computation and call NWTA again.
! Positive values indicate that NWTA believes the computation is
! finished.

! Assign values to OPT.

! *************************************************
! **  CALL NWT ( NFlag, AJ, Opt )                **
! *************************************************

!     You will need a Matrix J, and seven vectors:  F, X, "Best X", DX,
!     "Candidate DX" Gradient and "Best Gradient".

!     Set options using OPT as explained below, or set values
!     in AJ%O (see type NWT_Options below).

!     Set X( ) = initial guess for the solution, and
!         AJ%AXMAX = MAXVAL(ABS(X(1:NX)))
!     DO
! *************************************************
! **    CALL NWTA ( NFlag, AJ )                  **
! *************************************************
!       SELECT CASE ( NFlag ) ! NFlag > 0 means "Done", <0 means "Continue"
!       CASE ( NF_EVALF )
!         IF ( too many function values ) EXIT
!         Compute f(x)
!         Compute the Jacobian matrix J if you feel like doing it now (and
!         maybe do the other stuff required for NF_EvalJ now too)
!         Set AJ%FNORM = L2 norm of f(x)
!         IF ( AJ%FNORM is small enough ) EXIT
!       CASE ( NF_EVALJ )
!         IF ( too many Jacobian values ) EXIT
!         Compute the Jacobian matrix J if you didn't do it when NFlag
!         was NF_EVALF:
!           J(K,L) = Partial of F(K) / W.R.T. X(L), K = 1:NF, L = 1:NX
!         Compute (negative of the) gradient =
!         -(Jacobian)**T * F.  This is the RHS of the normal equations
!         J^T * J * "Candidate DX" = -J^T * F.
!         Triangularize J using QR, or J^T J using Cholesky.
!         Set
!           AJ%DIAG = element on diagonal of triangular factor with smallest
!                   absolute value,
!           AJ%AJN = maximum L1 norm of column in upper triangular factor,
!           AJ%FNMIN = L2 Norm of F not in column space of the Jacobian after
!                   stabilization, ||F + J * "Candidate DX"||, (which can be
!                   gotten without solving for "Candidate DX": if -F is put as
!                   the last column of J before triangularization, either by QR
!                   or by Cholesky factoring the normal equations, this is
!                   R(N+1,N+1) (or U(N+1,N+1) after factorization),
!           AJ%GRADN = L2 norm of Gradient.
!           The Q factor of the Jacobian factored as QR is the same as the
!           Cholesky factor U of U^T U = J^T J, except signs on the diagonal
!           might be different.
!       CASE ( NF_LEV )
!         Compute quantities necessary to determine the Levenberg-Marquardt
!         stabilization parameter: Solve U^T q = dX where U is the Cholesky
!         factor of J^T J, compute |q|**2 and store it in AJ%QNSQ.  This is
!         for the Moré-Sorensen algorithm.
!       CASE ( NF_SOLVE )
!         Apply Levenberg-Marquardt stabilization with parameter = AJ%SQ,
!         and solve the least-squares problem
!           (Jacobian) * "candidate DX" ~ -F
!         for "candidate DX".
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
!         IF ( NFlag == NF_TOO_SMALL ) THEN
!           Take special action if requested accuracy is critical
!         END IF
!         If ( NFlag == NF_TOLX_BEST ) copy "best X" to "X"
!         Convergence to desired solution.  Do whatever you want to
!         with the solution.
!         EXIT ! unless you have a really good reason to continue
!       CASE ( NF_FANDJ )
!         There is probably an error in the way F or J is computed.
!         A warning has been printed by the error message processor.
!         IF ( you don't have confidence in F and J ) STOP
!         If you don't stop, a gradient move is forced
!           (next NFlag == NF_Gmove).
!       END SELECT
!       IF ( you want to return to a previous best X ) NFlag = NF_START
!     END DO

  use Ermsg_m, only: Ermsg
  use DNWT_Type, only: RK
  use Output_M, only: NewLine, Output

  implicit NONE

  private
  public :: RK, DNWT, DNWTA, DNWTDB
  public ::      NWT,  NWTA,  NWTDB, NWT_Options, NWT_t
  public :: FlagName

  ! Start, or restart with gradient move from "Best":
  integer, parameter, public :: NF_START = 0

  ! Reasons for returning to user.  See description of usage above.

  ! Reasons to continue

  !      Evaluate F and AJ%FNORM:
  integer, parameter, public :: NF_EVALF = NF_START - 1

  !      Evaluate J and do other things:
  integer, parameter, public :: NF_EVALJ = nf_evalf-1

  !      Calculate quantities necessary to determine the Levenberg-Marquardt
  !      stabilization parameter:
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

  !      Set "Best X" = X, "Best Gradient" = Gradient:
  integer, parameter, public :: NF_BEST = nf_gmove-1

  !      Set DX = DX - "Candidate DX",
  !          AJ%DXDX = dot_product( DX, DX )
  !      IF ( AJ%DXDX /= 0.0 ) Set AJ%DXDXL = dot_product( DX, "Candidate DX" ):
  integer, parameter, public :: NF_AITKEN = nf_best-1

  !      Set DX = "Candidate DX":
  integer, parameter, public :: NF_DX = nf_aitken-1

  !      Set DX = AJ%CAIT * "Candidate DX":
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

  !      In case you want to add some flags for your own use...  If ?NWTA sees
  !      flags that are out of range, it prints a message but otherwise does
  !      nothing.
  integer, parameter, public :: NF_BIGGEST_FLAG = NF_FANDJ
  integer, parameter, public :: NF_SMALLEST_FLAG = NF_DX_AITKEN

  real(rk), parameter, private :: EPS = epsilon(1.0_rk), CBIG = huge(1.0_rk)

  ! Options for DNWT:

  type NWT_Options
    real(rk) :: AJSCAL = eps      ! Approximate absolute error in computing the
                                  ! Jacobian
    real(rk) :: DXMAXI = 0.1*cbig ! Largest value permitted for DXINC at any
                                  ! time
    real(rk) :: DXNOIS = 0.0      ! 1.0E-4_rk * largest move that gives a new
                                  ! best X (used to see if increase in AJ%FNORM
                                  ! might be due to too small a step)
    integer :: K1IT = 0           ! Number of iterations for which to give
                                  ! internal output
    real(rk) :: RELSF = max(0.01_rk,eps) ! Convergence is indicated with
                                  ! NFlag = NF_TOLF if AJ%FNORM - AJ%FNMIN <
                                  ! RELSF * AJ%FNMIN
    real(rk) :: SPMINI = eps      ! Nominal minimum value for normalized
                                  ! Marquardt parameter
    real(rk) :: SPSTRT = 0.1_rk   ! Starting value for normalized Marquardt
                                  ! parameter
    real(rk) :: TOLXA = tiny(1.0_rk)**0.75 ! Absolute tolerance on X. 
                                  ! Convergence is indicated with NFlag =
                                  ! NF_TOLX or NF_TOLX_BEST if the (possibly
                                  ! scaled) norm of the difference between the
                                  ! current X and the solution is estimated to
                                  ! be <= TOLXA.
    real(rk) :: TOLXR = eps**0.75 ! Relative tolerance on X.  Convergence is
                                  ! indicated with NFlag = NF_TOLX or
                                  ! NF_TOLX_BEST if the (possibly scaled) norm
                                  ! of the difference between the current X and
                                  ! the solution is estimated to be <= TOLXR *
                                  ! maxval (abs(X)).
  end type NWT_Options

  type(nwt_options), parameter :: Default_Options = NWT_options()

  type NWT_t               ! Stuff about the problem, neatly packaged.  This
                           ! is the type of the AJ argument of NWT*.
    type(nwt_options) :: O = Default_Options

    ! Inputs to ?NWTA.  Values to be set depend upon NFlag.
    ! See usage instructions above.
    real(rk) :: AJN        ! Largest L1 norm of column in upper triangle
                           ! of factored Jacobian matrix.  See NF_EvalJ
                           ! description above.
    real(rk) :: AXMAX      ! MAXVAL(ABS(X)).  See NF_NewX description above.
    real(rk) :: DIAG       ! Smallest | diagonal element | after factoring the
                           ! Jacobian matrix.  See NF_EvalF description above.
    real(rk) :: DXDX       ! dot_product( DX, DX ) (Newton step length)^2.
                           ! See NF_Aitken description above.
    real(rk) :: DXDXL      ! dot_product( "candidate DX", DX ).  See NF_Aitken
                           ! description above.
    real(rk) :: DXN        ! L2 Norm of candidate DX (Newton step length).  See
                           ! NF_Solve description above.
    real(rk) :: FNMIN      ! L2 Norm of F not in column space of the Jacobian.
                           ! See NF_EvalJ description above.
    real(rk) :: FNORM      ! L2 Norm of F at current X.  See NF_EvalF
                           ! description above.
    real(rk) :: GDX        ! dot_product( Gradient, "Candidate DX" ).  See
                           ! NF_Solve description above.
    real(rk) :: GRADN      ! L2 norm of Gradient = || J^T F ||.  See NF_EvalJ
                           ! description above.
    real(rk) :: QNSQ       ! Square of norm of solution of U^T q = x.  See
                           ! NF_Lev description above.
    logical :: BIG         ! ANY( DX > 10.0 * epsilon(X) * X ).  See NF_NewX
                           ! description above.

    ! Outputs from ?NWTA.  DO NOT CHANGE THESE!
    real(rk) :: CAIT = 0.0 ! Candidate factor for Aitken acceleration.  See
                           ! NF_DX_Aitken description above.
    real(rk) :: CGDX = 0.0 ! COS of angle between gradient and DX.   Not needed
                           ! from return until next call.  Here for printing by
                           ! ?NWTDB or the user.
    real(rk) :: GFAC = 1.0 ! Factor by which the gradient is to be multiplied
                           ! when taking a gradient move from the best X.  See
                           ! NF_Gmove description above.
    real(rk) :: SQ = 0.0   ! Actual value of Marquardt parameter. If SQ < 0,
                           ! modified Aitken acceleration is being used and
                           ! ABS(SQ) gives the acceleration factor.  See
                           ! NF_Solve description above.
    real(rk) :: SQMIN = 0.0 ! Current minimum value of SQ; recommend to use
                           ! this as a floor for Levenberg-Marquardt
                           ! stabilization in the caller.
    logical :: STARTING = .true. ! NWTA is still in "starting up" phase

    ! Values retained from return until next call.  DO NOT CHANGE THESE!
    real(rk) :: AXMAXB     ! AXMAX at best X (not actually used)
    real(rk) :: CDXDXL     ! COS of angle between DX (current move) and
                           ! DXL (previous move)
    real(rk) :: CONDAI     ! Crude estimate of reciprocal of condition number
                           ! of the Jacobian
    real(rk) :: DXBAD      ! This size for DX is likely a disaster
    real(rk) :: DXFAIL     ! Like DXINC, but more so.
    real(rk) :: DXI        ! Largest factor by which stepsize is allowed to
                           ! increase
    real(rk) :: DXINC      ! Largest value for DXN currently allowed
    real(rk) :: DXNBIG = 0.1*cbig ! A size of move that should be safe.
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
    integer :: INC         ! Flag set to indicate some past history
                           ! =-1  Starting
                           ! =0   Last X == best X
                           ! =J>0 AJ%FNORM > AJ%FNORML J times since last best X
    integer :: ITER        ! Number of current iteration
    integer :: ITKEN       ! Index used to keep track of Aitken accelerations
    integer :: KB          ! Flag set to indicate some past history
                           ! =2    Starting
                           ! =-J   (J>0) Making gradient move from best X
                           !       (-KB serves as a counter)
                           ! =1    Last F is bigger than the best F
                           ! =0    (INC >0) trying a move in wrong direction
                           !       from the best X
                           ! =0    (INC =0) new best X
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
    real(rk) :: SQT = 0.0  ! Total Levenberg-Marquardt stabilization
    integer :: WHERE_TO = NF_Start   ! Internal value of NFlag, to remember
                           ! reason for return, and therefore where to go next.

  contains
    procedure, pass(AJ) :: DNWT
    generic :: NWT => DNWT
    procedure, pass(AJ) :: DNWTA
    generic :: NWTA => DNWTA
    procedure, pass(AJ) :: DNWTDB
    generic :: NWTDB => DNWTDB
  end type NWT_t

  interface NWT; module procedure DNWT; end interface
  interface NWTA; module procedure DNWTA; end interface
  interface NWTDB; module procedure DNWTDB; end interface

  ! This should be internal to DNWTDB, but as of 2 July 2019 some processors
  ! still do not support internal procedures as specific procedures for a
  ! generic identifier.

  interface Add_To_Line
    procedure Add_To_Line_I, Add_To_Line_L, Add_To_Line_R
  end interface

! *****     Private data     *******************************************

  character(len=*), parameter :: ME = 'DNWT'

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

! ***************************************************     DNWT     *****

  subroutine DNWT ( NFlag, AJ, Opt )

! **********************************************************************

!>> 2000-03-21 W. V. Snyder Converted to Fortran 90.
!>> 1988-07-21 C. L. Lawson adapting code for the MATH77 library.
!-----------------------------------------------------------

! Original design and code by FRED T. KROGH
! Modifications for new usage coded by S. SINGLETARY
! Modifications to add bounds, JULY 1977
! JET PROPULSION LAB, PASADENA    JULY, 1974
! Re-structured by VAN SNYDER, JPL, 17 July 2018

! Arguments in the calling sequence are of the following types:

    integer, intent(out) :: NFlag
    class (NWT_t), intent(inout) :: AJ
    type (NWT_Options), intent(in), optional :: Opt

! The above arguments are used as follows:

! NFlag  = Integer used for communication with the user.  NWT sets NFlag =
!     NF_START.  NFlag is used by NWTA.  When NWTA returns control to the
!     calling program unit, NFlag is used as shown above.

! AJ =  NWTA internal state.

! Opt = Options; see type NWT_Options above.  Options can also be set by
!       the caller using AJ%O = Opt.

    aj%where_to = nf_start
    nflag = aj%where_to
    if ( present(opt) ) aj%o = opt

  end subroutine DNWT

! **************************************************     DNWTA     *****

  subroutine DNWTA ( NFlag, AJ )

! **********************************************************************

! Arguments in the calling sequence are of the following types:

    integer, intent(inout) :: NFlag    ! See usage instructions above
    class (nwt_t), intent(inout) :: AJ ! "About the Jacobian."  See usage
                                       ! instructions above.

!-----     Parameters and Unsaved Local Variables     ------------------

    integer, parameter :: IncBig = huge(0) / 2
    real(rk) :: TP       ! Widely used temporary variable.  Not needed from
                         ! one call to the next.

! RND      Used to decide if answer has been determined to machine
!          accuracy.

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

      if ( aj%where_to == nf_start ) then ! Initialization
      ! Initialization
        aj%ajn = 0.0
        aj%condai = 0.0
        aj%dxnl = 1.0
        aj%fnormb = sqrt(huge(aj%fnormb))
        aj%dxfail = huge(aj%dxfail)
        aj%dxbad = aj%dxfail
        aj%kfail = 0
        aj%fnorml = 0.0
        aj%inc = -1 ! Starting
      ! aj%inc = 0  ! Pretend we're at a "best"
        aj%iter = 0
        aj%itken = 0
        aj%kb = 2
        aj%spl = 0.0
        aj%sql = 0.0
        aj%sq = 0.0
        aj%dxdxl = 0.0 ! So it's not undefined, because the user isn't expected
                       ! to set it when starting.
        call do_NF_Evalf
      else ! User is forcing retreat to best X followed by gradient move
        call Do_NF_GMove_First
      end if

    else if ( nflag < nf_smallest_flag .or. nflag > nf_biggest_flag ) then

      call ermsg ( me // 'A', 98, 2, &
        & 'Invalid value of NFLAG', '.' )

    else if ( nflag /= aj%where_to ) then

      call ermsg ( me // 'A', 98, 2, &
        & 'NFLAG was changed to be inconsistent with AJ%WHERE_TO', '.' )

    else

! Continuing after return for reverse communication:

      select case ( aj%where_to )
      case ( nf_evalf ) ! Resumed here after evaluating F.

        aj%iter = aj%iter + 1

        ! Test if a retreat to a previous -- best -- X should be made

        tp = 0.9 * aj%spl * aj%ajn ! .9 times aj%sql on previous iteration
        if ( aj%fnorm < aj%fnorml ) then ! Things are getting better
          !{ \hspace*{20mm} $|F| < |F|_{\text{last}}$.
          if ( aj%fnorm >= aj%fnormb ) then ! But no better than the best X.
            !{ \hspace*{25mm} $|F| \geq |F|_{\text{best}}$.

            ! If our Marquardt parameter is small and aj%fnorm does not appear
            ! to be getting close to smaller than aj%fnormb, and we are not
            ! making a move in a quite different direction we give up and go
            ! back to previous best.

            if ( (max(tp, aj%sql) <= aj%sqmin) .and. &
               & ((aj%fnorm * (aj%fnorm/aj%fnorml)**2) > aj%fnormb) .and. &
               & (aj%cdxdxl >= 0.25) ) then
              call Do_NF_GMove_First
              return ! -----  GMove  ---------------------------------->
            end if
          end if
        else ! ( aj%fnorm >= aj%fnorml ) -- F norm increased

          !{ \hspace*{20mm} $|F| \geq |F|_{\text{last}}$.

          ! The residual expected if the problem is locally linear is
          ! aj%fnmin. Hopefully, the definition of "locally linear" is
          ! consistent with the Levenberg-Marquardt parameter. If aj%fnmin is
          ! smaller than the best one -- aj%fnminb, maybe it's OK to take
          ! another Newton move after the norm increases. But we don't want to
          ! do that if there have been a lot of "uphill" Newton moves.

          ! AJ%FNORML will have been set to zero if there are four consecutive
          ! gradient moves.

          ! Test for X convergence

          if ( x_converge() ) return ! -----  TolX or Tolx_Best  ------>
          if ( aj%inc >= 0 ) then ! We are not starting
            if ( aj%inc == 0 ) then ! Last X == Best X
              aj%dxnbig = max(aj%dxnl, aj%o%dxnois)

            else      ! aj%inc > 0 -- Not starting, Last X not Best X
              if ( aj%kb == 0 ) then ! Trying to move in wrong direction from
                                     ! best X.
                if ( aj%fnorm < aj%fnormb ) then
                  !{  \hspace*{35mm} $|F| < |F|_{\text{best}}$.
                  aj%gfac = -aj%gfac
                  if ( aj%gfac >=  0 ) then
                    ! Error processing -- F and J appear not to be consistent.
                    call ermsg ( me, 2, 0, 'J or F might be in error', '.' )
                    aj%where_to = nf_fandj
                    nflag = aj%where_to
                    return ! -----  FandJ  ---------------------------->
                  end if
                end if
                call Do_NF_GMove_Small
                return ! -----  GMove  -------------------------------->
              end if
              if ( aj%dxnl <= aj%o%dxnois ) then ! Last Newton move tiny?
                call do_nf_evalj ( tp )
                return ! ---------------------------------------------->
              end if
              if ( aj%kb < 0 ) then ! Gradient move last time?
                call Do_NF_GMove_Second
                return ! -----  GMove  -------------------------------->
              end if

! WVS revised this 2012-08-27 based upon advice from FTK
!             if ( (tp >= min(aj%sql,aj%sqb+aj%sqb)) .or. (aj%inc >= incbig) ) then
              if ( aj%fnmin > (1.0+0.5*aj%o%relsf)*aj%fnminb .or. &
                 & aj%inc >= incbig ) then
                ! Residual of locally-linear problem is too much bigger than
                ! the best estimate of the residual of the locally-linear
                ! problem, or we took an "uphill" Newton move after a
                ! gradient move after starting.
                call Do_NF_GMove_First
                return ! -----  GMove  -------------------------------->
              end if
            end if
            aj%sqmin = min ( aj%sqb, 4.0 * max( aj%spl, aj%spact ) * aj%ajn )
            aj%dxinc = max ( 0.5 * aj%dxnl, aj%o%dxnois )
            aj%inc = aj%inc + 1 ! Increase count of "uphill" Newton moves
          end if
        end if

        call do_NF_EvalJ ( tp ) ! -----  EvalJ  ----------------------->

      case ( nf_evalj ) ! Resumed here after computing Jacobian

        if ( aj%gradn <= 0.0 ) then
           aj%gradn = 1.0
           if ( aj%ajn <= 0.0 ) aj%ajn = 1.0
        end if
        tp = aj%diag / aj%ajn
        aj%condai = max ( min ( aj%condai,tp ), tp**2 )
        aj%spact = 0.25 * aj%condai

        !{ \hspace*{15mm}FRZ $= \sqrt{ | ( |F|^2 - |F|^2_{\text{min} } |) }$
        aj%frz = sqrt ( abs ( (aj%fnorm-aj%fnmin) * (aj%fnorm+aj%fnmin) ) )

        if ( aj%iter == 1 ) then ! First iteration
          aj%sp = aj%o%spstrt
          aj%spinc = aj%sp
          if ( aj%frz <= rnd * aj%fnmin ) then
            call Do_NF_Tolf
            return ! -----  TolF  ------------------------------------->
          end if
          aj%dxi = 0.125
        else
          if ( aj%fnorm >= aj%fnorml ) then

            !{ \hspace*{24mm} $|F| \geq |F|_{\text{last}}$.

            ! Select new Levenberg-Marquardt stabilizing parameter for the
            ! case when AJ%FNORM has increased -- an "uphill" Newton move.

            aj%kb = 1
            if ( aj%dxnl < aj%o%dxnois ) then
              aj%sp = 1.0e-4 * aj%sp
              aj%sqmin = min ( aj%sqmin, aj%sp * aj%ajn )
!             aj%dxnl = aj%dxnl * 100.0
              call do_NF_Solve ( before_Solve = set_new_marquardt )
              return ! -----  Solve  ---------------------------------->
            end if
            aj%spfac = max ( (aj%fnorm/aj%fnorml)**2, 10.0_rk )
            aj%spinc = aj%spl + aj%spl + aj%spact
            aj%dxi = 1.0
            call do_NF_Solve ( before_Solve = test_GradN_Increase )
            return ! -----  Solve  ------------------------------------>
          end if
          if ( aj%fnorm >= aj%fnormb ) then ! Is function norm larger than best?
            !{ \hspace*{24mm} $|F| \geq |F|_{\text{best}}$.
            if ( .not. f_converge() ) &
              & call do_NF_Solve ( before_Solve = test_f_linear )
            return ! -----  Solve  ------------------------------------>
          end if
        end if

        ! New best X

        if ( aj%kfail /= 0 ) then
           if ( aj%kb == 0 ) then ! Got new best immediately
              if ( aj%dxn >= 0.9_rk * aj%dxbad ) then
                 aj%dxbad = aj%dxbad * (1.25_rk ** aj%kfail)
                 aj%dxfail = aj%dxbad
                 aj%kfail = aj%kfail + 1
              else
                 aj%dxfail = 0.75_rk * aj%dxfail + 0.25_rk * aj%dxbad
                 aj%kfail = 1
              end if
           else
              aj%dxfail = min ( 1.0625_rk * aj%dxfail, aj%dxbad )
           end if
           aj%dxinc = min ( aj%o%dxmaxi, aj%dxfail )
        else
           aj%dxinc = aj%o%dxmaxi
        end if

        aj%spinc = min ( aj%spinc, &
                       & max ( aj%o%spmini, abs(aj%o%ajscal) / aj%ajn ) )
        aj%sqmin = 0.0
        aj%spb = aj%spl
        aj%sqb = aj%sql
        if ( aj%inc < 0 ) then ! Starting
          call do_NF_Solve ( before_Solve = set_new_marquardt )
          return ! -----  Solve  -------------------------------------->
        end if
        if ( aj%inc /= 0 ) then
          if ( aj%kb < 0 ) then
            aj%dxi = 0.25
          else if ( aj%kb == 1 ) then
          ! We had a function increase prior to getting a new best.
            aj%dxi = 0.25 * aj%dxnl / aj%dxnbig ! prior to getting a new best.
          end if
          aj%kb = 0
          aj%inc = 0
        end if

        ! Test for convergence in sense of AJ%FNORM < (aj%o%RELSF * aj%FNMIN)

        if ( .not. f_converge() ) &
           & call do_NF_Solve ( before_Solve = test_F_Linear )

        ! -----  TolF or Solve  --------------------------------------->

      case ( nf_solve ) ! Resumed here after doing Levenberg-Marquardt
                        ! stabilization and solving for candidate DX

        !{ \hspace*{14mm} {\tt CGDX}
        !  $= \cos( \nabla \mathbf{f} - \delta \mathbf{x} ) =
        !  \frac{\nabla \mathbf{f} \cdot \delta \mathbf{x}}
        !       {|\nabla \mathbf{f}|\,|\delta \mathbf{x}|}$
        aj%cgdx = aj%gdx / (aj%gradn * aj%dxn) ! Cosine ( gradient, dx )
        aj%condai = min ( aj%cgdx, &
                        & aj%gradn / ( aj%dxn * aj%ajn**2 ), &
                        & aj%diag / aj%ajn )
        aj%fnxe = aj%fnmin**2 - (aj%sq * aj%dxn)**2
        if ( aj%inc < 0 ) then ! Starting
          if ( aj%dxn <= aj%dxinc ) then ! Move length OK?
            aj%inc = incbig
            aj%cdxdxl = 0.0 ! Don't think about Aitken
            call do_nf_best ! Go save best X, then either DX or Aitken
            return ! -----  Best  ------------------------------------->
          end if
        else
          !{ \hspace*{18mm} {\tt CDXDXL }$= \cos( \delta \mathbf{x} -
          !                        \delta \mathbf{x}_{\text{last}}) =
          ! \frac{\delta \mathbf{x} \cdot \delta \mathbf{x}_{\text{last}}}
          !      {|\delta \mathbf{x}| \, |\delta \mathbf{x}_{\text{last}}|}$
          aj%cdxdxl = aj%dxdxl / ( aj%dxn * aj%dxnl )
          if ( aj%fnxe > aj%fnormb**2 ) then
            if ( aj%o%k1it /= 0 ) call nwtdb ( aj, width=9, why='Give up' )
            call Do_NF_GMove_First
            return ! -----  GMove  ------------------------------------>
          end if
          if ( aj%dxn < 1.1_rk * aj%dxinc )  then
            if ( aj%sq <= 0.0 .or. aj%dxn > 0.90_rk * aj%dxinc ) then
              if ( aj%inc == 0 ) then ! Previous move not gradient move
                if ( .not. x_converge() ) call do_nf_best
                return ! -----  TolX or TolX_Best or Best  ------------>
              end if
              aj%cait = cbig
              call do_nf_dx  ! Store Candidate DX as DX
              return ! -----  DX  ------------------------------------->
            end if
! WVS: was  aj%dxinc = aj%dxi * aj%dxn * (1.0 - aj%cdxdxl)**2
! This was clearly wrong, since it says "If consecutive Newton moves are in
! the same direction, take shorter steps"
            aj%dxinc = aj%dxi * aj%dxn / (1.025 - aj%cdxdxl)**2
          end if
        end if

        ! Newton step length is too large or too small

        if ( aj%o%k1it /= 0 ) call nwtdb ( aj, width=9, why='Step length' )

        aj%where_to = nf_lev
        nflag = aj%where_to

        !{ Calculate {\tt qn}$^2 = {\bf q}^T{\bf q} = ||{\bf q}||^2$, used to
        !  determine the Levenberg-Marquardt parameter $\lambda$, where ${\bf q}
        !  = {\bf U}^{-T} {\delta\bf x}$ and ${\delta\bf x}$ is the Newton step
        !  computed with the current value of $\lambda$.  The method to
        !  determine $\lambda$ and the reason for being interested in ${\bf
        !  q}^T{\bf q}$ are explained below.

        ! -----  Lev  ------------------------------------------------->

      case ( nf_lev ) ! Resumed here after calculating QN**2.
      
!{ Determine the Levenberg-Marquardt parameter $\lambda$, known here as SQ,
!  using the Mor\'e-Sorensen method.
!  The last parameter used was SQ, which means that {\tt SQ}$^2 = \lambda^2$
!  was added to the diagonal of the normal equations.  A little bit of 
!  derivation so you have some idea of what is going on.  Start with
!
!  $({\bf H} + \lambda {\bf I})\, {\delta\bf x} = -{\bf g}$,
!
!  where ${\bf H}$ is the Hessian matrix.  In our case, ${\bf J}^T {\bf J}$,
!  where ${\bf J}$ is the Jacobian matrix, is an approximate Hessian. 
!  Differentiating both sides with respect to $\lambda$ gives
!
!  $({\bf H} + \lambda {\bf I})\, \frac{\partial{\delta\bf x}}{\partial\lambda}
!  + {\bf I}\, {\delta\bf x} = 0$ or
!  $({\bf H} + \lambda {\bf I})\, \frac{\partial{\delta\bf x}}{\partial\lambda}
!  = -{\delta\bf x}$.
!
!  Cholesky factoring $({\bf H} + \lambda {\bf I})$, we have
!  ${\bf U}^T {\bf U} \frac{\partial{\delta\bf x}}{\partial\lambda} =
!  -{\delta\bf x}$ or ${\bf U} \frac{\partial{\delta\bf x}}{\partial\lambda} =
!  -{\bf U}^{-T} {\delta\bf x} = -{\bf q}$.
!
!  Multiplying ${\bf U}^T {\bf U} \frac{\partial{\delta\bf x}}{\partial\lambda} =
!  -{\delta\bf x}$ on the left by $\frac{\partial{\delta\bf x}}{\partial\lambda}^T$
!  and substituting
!  ${\bf U} \frac{\partial{\delta\bf x}}{\partial\lambda} = - {\bf q}$ we have
!  $\frac{\partial{\delta\bf x}}{\partial\lambda}^T {\bf U}^T {\bf U}
!   \frac{\partial{\delta\bf x}}{\partial\lambda}
!  = -\frac{\partial{\delta\bf x}}{\partial\lambda}^T {\delta\bf x} =
!  {\bf q}^T {\bf q}$.
!
!  We intend to find a zero of
!  $\phi = \frac1{|| {\delta\bf x} ||} - \frac1{\rho}$, where $\rho$ is the
!  desired step length for ${\delta\bf x}$,
!  treating ${\delta\bf x}$ as a function of $\lambda$.
!
!  $\frac{\partial\phi}{\partial\lambda} =
!  \frac{\partial ( (|| {\delta\bf x} ||^2)^{-1/2}
!   - \frac1{\rho})}{\partial{\delta\bf x}}
!  \frac{\partial{\delta\bf x}}{\partial\lambda}
!  = \frac{-1}{|| {\delta\bf x} ||^{3}} {\delta\bf x}^T
!    \frac{\partial{\delta\bf x}}{\partial\lambda}
!  = \frac{-1}{|| {\delta\bf x} ||^{3}} (- {\bf q}^T {\bf q})$
!
!  The Newton method gives us
!
!  $\lambda_{\text{new}} = \lambda_{\text{old}} -
!    \phi / \frac{\partial\phi}{\partial\lambda_{\text{old}}}
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
! {\,\delta\bf x} = -{\bf J}^T {\bf f} = -{\bf g}$.

        aj%sq = sqrt ( max( 0.0_rk, &
                          & aj%sq**2 - (aj%dxn**2 / aj%qnsq) * &
                                     & (1.0_rk - aj%dxn / aj%dxinc) ) )
        aj%sp = aj%sq / aj%ajn
        call do_NF_Solve ( before_Solve = solve_only )

        ! -----  Solve  ----------------------------------------------->

      case ( nf_newx ) ! Resumed here after computing new X

        if ( .not. aj%big ) then ! All dx components are very small
          if ( (aj%spl <= aj%o%spmini) .or. (aj%sq == 0.0) ) then
            if ( (aj%inc > 0) .and. (aj%kb /= 0) ) then ! Go do gradient move
              call Do_NF_GMove_First
            else ! Convergence -- move too small
              nflag = nf_too_small
              aj%where_to = nflag
            end if
            return ! -----  GMove or TooSmall  ------------------------>
          end if
          ! Set so steps don't get too small:
          aj%o%dxnois = max ( aj%o%dxnois, 10.0 * aj%dxn )
        end if
        call do_NF_EvalF

      case ( nf_gmove ) ! Resumed here after computing gradient step from best X

        ! In case the caller wants to compute cosines:
        aj%dxn = aj%gradnb * aj%gfac
        if ( aj%gfac <= 0.0 ) then
          call do_NF_NewX
        else
          aj%axmax = aj%axmaxb
          aj%gfac = 0.125 * aj%gfac
          aj%spl = aj%spg
          aj%spg = 8.0 * aj%spl
          aj%fnorm = aj%fnormb
          aj%fnxe = 0.5 * aj%fnorm **2
          aj%frz = aj%frzb
          aj%gradn = aj%gradnb
          call do_NF_NewX_First
        end if

        ! -----  NewX  ------------------------------------------------>

      case ( nf_best ) ! Resumed here after saving X as "best X" and
                       ! Gradient as "Best gradient"

        aj%axmaxb = aj%axmax
        aj%o%dxnois = max ( aj%o%dxnois, 1.0e-4 * aj%dxnl )
        aj%gfac = min ( 0.125 * aj%dxn / aj%gradn, aj%ajn**(-2) )
        aj%spg = aj%spl + 0.01
        aj%fnormb = aj%fnorm
        aj%frzb = aj%frz
        aj%fnminb = aj%fnmin
        aj%gradnb = aj%gradn
        if ( abs(aj%cdxdxl) < 0.9 .or. aj%sq /= 0.0 .or. aj%dxn >= aj%dxnl ) then
          call do_nf_dx ! Store Candidate DX as DX
        else

          ! Compute scalar factor for (modified) Aitken acceleration

          aj%where_to = nf_aitken
          nflag = aj%where_to
        end if

        ! -----  DX or Aitken  ---------------------------------------->

      case ( nf_aitken ) ! Resumed here after computing some numbers needed
                         ! for (modified) Aitken acceleration

        tp = aj%dxdx
        if ( tp /= 0.0 ) then
           tp = 1.0 + aj%dxdxl/tp
           tp = min ( 10.0_rk, max ( tp, 0.01_rk ) )
           if ( abs(aj%cait-tp) < abs ( 0.25 * (aj%cait+tp)-0.5 ) ) then
              ! Set DX = Aitken-modified DX = aj%CAIT * "Candidate DX"
              aj%cait = tp
              aj%itken = aj%iter + 2
              aj%sq = -aj%cait
              aj%dxn = aj%cait*aj%dxn
              aj%where_to = nf_dx_aitken
              nflag = aj%where_to
              return ! -----  Dx_Aitken  ------------------------------>
           end if
           if ( aj%itken < aj%iter ) aj%cait = tp
        end if

        ! End of logic for Aitken acceleration

        call do_nf_dx ! Store Candidate DX as DX

        ! -----  DX  -------------------------------------------------->

      case ( nf_dx, nf_dx_aitken ) ! Resumed here after DX = "Candidate DX"
                                   ! or DX = aj%CAIT * "Candidate DX"

        aj%dxnl = aj%dxn
        call do_NF_NewX_First

        !  -----  NewX  ----------------------------------------------->
! Continuing after return that was expected to be final:

      case ( nf_tolx, nf_tolx_best )
        call do_NF_Best ! -----  Best  -------------------------------->
      case ( nf_tolf )
        call do_NF_Solve ( before_Solve = test_F_Linear )
        ! -----  Solve  ----------------------------------------------->
      case ( nf_too_small )
        call do_nf_gmove ! -----  GMove  ------------------------------>
      case ( nf_fandj )
        call do_NF_GMove_Small ! -----  GMove  ------------------------>
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

      aj%kb = -1 ! Set consecutive gradient move counter
      aj%dxbad = 0.75 * min ( aj%dxn, aj%dxbad )
      aj%dxfail = 0.125 * aj%dxbad
      aj%kfail = 1
      call do_NF_GMove_Second
    end subroutine Do_NF_GMove_First

    ! .......................................  Do_NF_GMove_Second  .....
    subroutine Do_NF_GMove_Second

    ! Set up for a gradient move after a previous gradient move, or
    ! finish setting up for a first gradient move.

      aj%kb = aj%kb - 1
      if ( aj%kb <= -4 ) then   ! Four consecutive gradient moves !

        ! Setup to test if Jacobian matrix is being computed properly

        aj%kb = 0
        aj%fnorml = 0.0
        aj%gfac = -100.0 * aj%gfac
      end if
      call do_NF_GMove
    end subroutine Do_NF_GMove_Second

    ! ........................................  Do_NF_GMove_Small  .....
    subroutine Do_NF_GMove_Small

    ! Set up for a really small gradient move.

      aj%kb = -1 ! Restart consecutive gradient move counter
      aj%gfac = -0.01_rk * aj%gfac
      call do_NF_GMove
    end subroutine Do_NF_GMove_Small

    ! .........................................  Do_NF_NewX_First  .....
    subroutine Do_NF_NewX_First

    ! Compute new X and test for too small a correction

      aj%fnorml = aj%fnorm
      aj%frzl = aj%frz
      aj%fnxe = 0.25 * aj%fnorm**2 + 0.76_rk * aj%fnxe
      call do_NF_NewX
    end subroutine Do_NF_NewX_First


    ! ...............................................  Do_NF_NewX  .....
    subroutine Do_NF_NewX

    ! Compute new X and test for too small a correction

      if ( aj%o%k1it /= 0 ) then ! Do requested debugging output
        aj%o%k1it = aj%o%k1it - 1
        call nwtdb ( aj, width=9, why='Before new X' )
      end if

      aj%starting = aj%inc < 0
      nflag = nf_newx
      aj%where_to = nflag
    end subroutine Do_NF_NewX

    ! ..............................................  Do_NF_Solve  .....

    subroutine Do_NF_Solve ( Before_Solve )
      integer, intent(in) :: Before_Solve

      ! Select new stabilizing parameter.  Also comes back here if user
      ! continues after TOLF convergence.

      if ( before_solve >= test_F_Linear ) then

        tp = min ( aj%frz / aj%frzl, aj%fnorm / aj%fnorml )

        ! Test if F appears almost linear over last step

        if ( aj%fnorm**2 < 1.125 * aj%fnxe ) then ! F appears almost linear.
          aj%spfac = 0.125 * (1.025 - aj%cdxdxl) * tp**2
          if ( aj%spl <= aj%spinc ) aj%spinc = 0.25 * aj%spinc
          aj%dxi = min ( aj%dxi, 0.25_rk )
        else ! F not linear over last step
        ! aj%spfac = tp*(aj%fnorm**2)/aj%fnxe
        ! if ( aj%cdxdxl >= 0.9_rk) aj%spfac = min(aj%spfac,0.5_rk)
        ! On 2002/07/25, FTK recommended:
          aj%spfac = min ( tp * (aj%fnorm**2) / aj%fnxe, &
                         & 32.0 * (1.025 - aj%cdxdxl)**2 )
          aj%dxi = 1.0
        end if
      end if
      if ( before_solve >= test_GradN_Increase ) then
        if ( aj%gradn > aj%gradnl ) aj%spfac = aj%spfac * aj%gradnl / aj%gradn
        aj%sp = max ( aj%spinc, min ( aj%spb, aj%spl * aj%spfac ) )
      end if
      if ( before_solve >= set_New_Marquardt ) then
        aj%sp = max ( aj%sp,aj%sqmin / aj%ajn )
    !   if ( aj%inc == 0 ) aj%sp = min(aj%sp, 0.5 * aj%spl)
        aj%spl = aj%sp
        aj%sq = aj%sp * aj%ajn
        aj%sqt = aj%sq
        if ( aj%sp < aj%spact .and. aj%inc >= 0 ) aj%sq = aj%sqmin
        aj%gradnl = aj%gradn
      end if

    ! Do Levenberg-Marquardt stabilization, solve for "Candidate DX"

      aj%starting = aj%inc < 0
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

    ! Test for convergence in sense of AJ%FNORM < (1 + AJ%O%RELSF)*AJ%FNMIN

      f_converge = .false.
      if ( aj%inc == 0 ) then
        tp = aj%fnorm - ( 1.0 + aj%o%relsf ) * aj%fnmin
        if ( tp < 0.0 ) then
                    ! WVS changed this 2012-09-10 to allow convergence
                    ! with nonzero Levenberg-Marquardt parameter
!           if ( (aj%sq <= aj%sqmin) .or. (aj%spl <= aj%o%spmini) .or. &
!            &   (tp <= -(1.0+aj%o%relsf)*aj%spl*aj%fnorm) ) then
          call do_nf_tolf
          f_converge = .true.
        end if
      end if

    end function F_Converge

    ! ...............................................  X_Converge  .....
    logical function X_Converge ( )
    ! Test whether convergence has occurred due to very small \delta X
      x_converge = .false.
      if ( aj%fnorml > 0.0 ) then ! Not doing gradient moves from the best X
        if ( aj%dxn <= aj%o%tolxa .or. &
           & ((aj%sq == 0.0).or.(aj%spl <= aj%o%spmini)) .and. &
             & aj%dxn <= aj%o%tolxr * aj%axmax ) then
          ! Convergence if aj%dxn is small enough or we have a suffiently small
          ! Levenberg-Marquardt parameter and we had a very small move.
          aj%where_to = nf_tolx
          if ( aj%kb /= 0 ) then
            aj%fnorm = aj%fnormb
            aj%where_to = nf_tolx_best
          end if
          nflag = aj%where_to
          x_converge = .true.
        end if
      end if
    end function X_Converge

  end subroutine DNWTA

! *************************************************     DNWTDB     *****

  subroutine DNWTDB ( AJ, WIDTH, WHY )

!   Print the scalars in the module.  Print the stuff in AJ if it's
!   present.  Print the integer scalars first.  Then print
!   max(5,min(9,WIDTH)) real scalars per line if width is present, else
!   print five per line.

    class (NWT_T), intent(in) :: AJ
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: WHY ! printed if present

    integer :: I                        ! Which one in the line is being worked
    character(len=9) :: IFLname, NFLname
    integer :: MyWidth
    character(len=132) :: Name_Line     ! For names
    character(len=132) :: Output_Line   ! For values

    myWidth = 5
    if ( present(width) ) myWidth = max(5,min(9,width))

    output_line = '-----     DNWT options     ---------------------------------&
      &------------------------------------------------------------------'
    call output ( output_line(1:14*myWidth), advance='yes' )

    call flagName ( aj%where_to, iflName ); call flagName ( aj%where_to, nflName )

!   write ( *, '(a)' ) &
    call output( '  WHERE_TO        INC    ITER   ITKEN    K1IT      KB' )
    if ( present(why) ) call output ( '  WHY' )
    call newLine
    write ( output_line, '(1x,a9,i11,4i8)' ) adjustr(iflName), aj%inc, &
      & aj%iter, aj%itken, aj%o%k1it, aj%kb
    call output ( trim(output_line) )
    if ( present(why) ) call output ( '  ' // trim(why) )
    call newLine

    i = 1
    name_line = ''
    output_line = ''
    call add_to_line ( aj%o%ajscal, 'AJSCAL', name_line, i, myWidth, output_Line )
    call add_to_line ( aj%o%dxmaxi, 'DXMAXI', name_line, i, myWidth, output_Line )
    call add_to_line ( aj%o%dxnois, 'DXNOIS', name_line, i, myWidth, output_Line )
    call add_to_line ( aj%o%relsf,  'RELSF',  name_line, i, myWidth, output_Line )
    call add_to_line ( aj%o%spmini, 'SPMINI', name_line, i, myWidth, output_Line )
    call add_to_line ( aj%o%spstrt, 'SPSTRT', name_line, i, myWidth, output_Line )
    call add_to_line ( aj%o%tolxa,  'TOLXA',  name_line, i, myWidth, output_Line )
    call add_to_line ( aj%o%tolxr,  'TOLXR',  name_line, i, myWidth, output_Line )
    if ( i /= 1 ) call print_lines ( name_line, output_line, i )

    output_line = '-----     DNWT variables     -------------------------------&
      &------------------------------------------------------------------'
    call output ( output_line(1:14*myWidth), advance='yes' )
    output_line = ''

    call add_to_line ( aj%ajn,        'AJN',      name_line, i, myWidth, output_Line )
    call add_to_line ( aj%axmax,      'AXMAX',    name_line, i, myWidth, output_Line )
    call add_to_line ( aj%axmaxb,     'AXMAXB',   name_line, i, myWidth, output_Line )
    call add_to_line ( aj%cait,       'CAIT',     name_line, i, myWidth, output_Line )
    call add_to_line ( aj%cdxdxl,     'CDXDXL',   name_line, i, myWidth, output_Line )
    call add_to_line ( aj%cgdx,       'CGDX',     name_line, i, myWidth, output_Line )
    call add_to_line ( aj%condai,     'CONDAI',   name_line, i, myWidth, output_Line )
    call add_to_line ( aj%diag,       'DIAG',     name_line, i, myWidth, output_Line )
    call add_to_line ( aj%dxbad,      'DXBAD',    name_line, i, myWidth, output_Line )
    call add_to_line ( aj%dxi,        'DXI',      name_line, i, myWidth, output_Line )
    call add_to_line ( aj%dxinc,      'DXINC',    name_line, i, myWidth, output_Line )
    call add_to_line ( aj%dxnbig,     'DXNBIG',   name_line, i, myWidth, output_Line )
    call add_to_line ( aj%dxfail,     'DXFAIL',   name_line, i, myWidth, output_Line )
    call add_to_line ( aj%dxdx,       'DXDX',     name_line, i, myWidth, output_Line )
    call add_to_line ( aj%dxdxl,      'DXDXL',    name_line, i, myWidth, output_Line )
    call add_to_line ( aj%dxn,        'DXN',      name_line, i, myWidth, output_Line )
    call add_to_line ( aj%dxnl,       'DXNL',     name_line, i, myWidth, output_Line )
    call add_to_line ( aj%fnmin,      'FNMIN',    name_line, i, myWidth, output_Line )
    call add_to_line ( aj%fnminb,     'FNMINB',   name_line, i, myWidth, output_Line )
    call add_to_line ( aj%fnorm,      'FNORM',    name_line, i, myWidth, output_Line )
    call add_to_line ( aj%fnormb,     'FNORMB',   name_line, i, myWidth, output_Line )
    call add_to_line ( aj%fnorml,     'FNORML',   name_line, i, myWidth, output_Line )
    call add_to_line ( sqrt(aj%fnxe), 'FNXE**.5', name_line, i, myWidth, output_Line )
    call add_to_line ( aj%frz,        'FRZ',      name_line, i, myWidth, output_Line )
    call add_to_line ( aj%frzb,       'FRZB',     name_line, i, myWidth, output_Line )
    call add_to_line ( aj%frzl,       'FRZL',     name_line, i, myWidth, output_Line )
    call add_to_line ( aj%gdx,        'GDX',      name_line, i, myWidth, output_Line )
    call add_to_line ( aj%gfac,       'GFAC',     name_line, i, myWidth, output_Line )
    call add_to_line ( aj%gradn,      'GRADN',    name_line, i, myWidth, output_Line )
    call add_to_line ( aj%gradnb,     'GRADNB',   name_line, i, myWidth, output_Line )
    call add_to_line ( aj%gradnl,     'GRADNL',   name_line, i, myWidth, output_Line )
    call add_to_line ( aj%kfail,      'KFAIL',    name_line, i, myWidth, output_Line )
    call add_to_line ( aj%qnsq,       'QNSQ',     name_line, i, myWidth, output_Line )
    call add_to_line ( aj%sp,         'SP ',      name_line, i, myWidth, output_Line )
    call add_to_line ( aj%spact,      'SPACT',    name_line, i, myWidth, output_Line )
    call add_to_line ( aj%spb,        'SPB',      name_line, i, myWidth, output_Line )
    call add_to_line ( aj%spfac,      'SPFAC',    name_line, i, myWidth, output_Line )
    call add_to_line ( aj%spg,        'SPG',      name_line, i, myWidth, output_Line )
    call add_to_line ( aj%spinc,      'SPINC',    name_line, i, myWidth, output_Line )
    call add_to_line ( aj%spl,        'SPL',      name_line, i, myWidth, output_Line )
    call add_to_line ( aj%sq,         'SQ',       name_line, i, myWidth, output_Line )
    call add_to_line ( aj%sqb,        'SQB',      name_line, i, myWidth, output_Line )
    call add_to_line ( aj%sql,        'SQL',      name_line, i, myWidth, output_Line )
    call add_to_line ( aj%sqmin,      'SQMIN',    name_line, i, myWidth, output_Line )
    call add_to_line ( aj%sqt,        'SQT',      name_line, i, myWidth, output_Line )
    call add_to_line ( aj%kfail,      'KFAIL',    name_line, i, myWidth, output_Line )
    call add_to_line ( aj%big,        'BIG',      name_line, i, myWidth, output_Line )
    call add_to_line ( aj%starting,   'STARTING', name_line, i, myWidth, output_Line )
    if ( i /= 1 ) call print_lines ( name_line, output_line, i )

    output_line = '------------------------------------------------------------&
      &------------------------------------------------------------------'
    call output ( output_line(1:14*myWidth), advance='yes' )

  contains
  end subroutine DNWTDB

  ! The next three subroutines ought to be internal to DNWTDB, but as of
  ! 2 July 2019, some processors still don't support Fortran 2008.  In
  ! Fortran 2003, internal procedures were prohibited to be specific
  ! procedures for a generic identifier.

  subroutine Add_To_Line_R ( Value, Name, Name_Line, I, MyWidth, Output_Line )
    real(rk), intent(in) :: Value
    character(len=*), intent(in) :: Name
    character(len=*), intent(inout) :: Name_Line
    integer, intent(inout) :: I
    integer, intent(inout) :: MyWidth
    character(len=*), intent(inout) :: Output_Line
    name_line(1+14*(i-1):10+14*(i-1)) = name
    name_line(1+14*(i-1):10+14*(i-1)) = adjustr(name_line(1+14*(i-1):10+14*(i-1)))
    write ( output_line(1+14*(i-1):14*i), '(es13.6)' ) value
    i = i + 1
    if ( i > myWidth ) call print_lines ( name_line, output_line, i )
  end subroutine Add_To_Line_R

  subroutine Add_To_Line_I ( Value, Name, Name_Line, I, MyWidth, Output_Line )
    integer, intent(in) :: Value
    character(len=*), intent(in) :: Name
    character(len=*), intent(inout) :: Name_Line
    integer, intent(inout) :: I
    integer, intent(inout) :: MyWidth
    character(len=*), intent(inout) :: Output_Line
    name_line(1+14*(i-1):10+14*(i-1)) = name
    name_line(1+14*(i-1):10+14*(i-1)) = adjustr(name_line(1+14*(i-1):10+14*(i-1)))
    write ( output_line(1+14*(i-1):10+14*(i-1)), '(I10)' ) value
    i = i + 1
    if ( i > myWidth ) call print_lines ( name_line, output_line, i )
  end subroutine Add_To_Line_I

  subroutine Add_To_Line_L ( Value, Name, Name_Line, I, MyWidth, Output_Line )
    logical, intent(in) :: Value
    character(len=*), intent(in) :: Name
    character(len=*), intent(inout) :: Name_Line
    integer, intent(inout) :: I
    integer, intent(inout) :: MyWidth
    character(len=*), intent(inout) :: Output_Line
    name_line(1+14*(i-1):10+14*(i-1)) = name
    name_line(1+14*(i-1):10+14*(i-1)) = adjustr(name_line(1+14*(i-1):10+14*(i-1)))
    write ( output_line(1+14*(i-1):1+14*(i-1)), '(L1)' ) value
    output_line(1+14*(i-1):10+14*(i-1)) = adjustr(output_line(1+14*(i-1):10+14*(i-1)))
    i = i + 1
    if ( i > myWidth ) call print_lines ( name_line, output_line, i )
  end subroutine Add_To_Line_L

  ! The next subroutine is referenced by Add_To_Line_*, and ought to be
  ! internal to DNWTDB.

  subroutine Print_Lines ( Name_Line, Output_Line, I )
    character, intent(inout) :: Name_Line, Output_Line
    integer, intent(out) :: I
    call output ( trim(name_line), advance='yes' )
    call output ( trim(output_line), advance='yes' )
    i = 1
    name_line = ''
    output_line = ''
  end subroutine Print_Lines

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
    case ( NF_LEV )
      itsName = 'LEV'
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
      itsName = 'DX_AITKEN'
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

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module DNWT_MODULE

! $Log$
! Revision 2.60  2019/09/23 20:04:38  vsnyder
! Add comments about Newton moves after 'uphill' moves
!
! Revision 2.59  2019/08/20 23:52:00  vsnyder
! Move Add_To_Line* from internal to module scope because some processors
! (NAG in particular) do not yet support internal procedures as specific
! procedures for a generic identifier.
!
! Revision 2.58  2019/06/24 23:29:26  pwagner
! Updated to reflect TA-01-143
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
