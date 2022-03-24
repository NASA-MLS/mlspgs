! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module MLSNumerics              ! Some low level numerical stuff
!=============================================================================

 use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
 use HyperSlabs, only: Rerank
 use Dump_0, Only : Dump
 use HighOutput, only: OutputNamedValue
 use Hunt_M, only: Hunt, Huntbox, Huntrange
 use MatrixModule_0, only: Createblock, M_Absent, MatrixElement_T, Sparsify
 use MLSFillValues, only: IsfillValue
 use MLSKinds, only: I4, R4, R8, Rm
 use MLSMessageModule, only: MLSMSG_Error, MLSMSG_Warning, &
   & MLSMessage
 use MLSFinds, only: Findfirst, Findlast
 use MLSStrings, only: Capitalize, Trim_Safe
 use Optional_M, only: Default
 use Output_M, only: Blanks, Output
 use Pure_Hunt_M, only: Purehunt
 use Symm_Tri, only: Factor_Symm_Tri, Solve_Factored_Symm_Tri

  implicit none

  private
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! This module contains some low level numerical stuff, hunting, interpolating,
  ! approximating roots, derivatives, integrals, and some stray functions
  ! An area of focus is approximating functions that are expensive to evaluate
  ! by assembling a lookup table of precomputed values
  ! (a) Using an array of its values at regularly spaced arguments
  ! (b) Using an array of its values at an array of specified arguments
  ! (c) Like (a) but using a user-defined datatype, the LookUpTable_0
  !
! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (parameters and datatypes)
! Coefficients_nprec       Coefficients to speed up interpolation
!                            (See InterpolateArraySetup)
! LookUpTable_0_nprec    Look Up Table
!                            y(x[i]) where x[i+1] - x[i] = constant

!         Functions, operations, routines
! AGM                      Return the arithmetic-geometric mean of two numbers
! Average                  Compute the average of a rank-one real array
! Battleship               By ever-widening evaluations find integer root
!                             or floating root to within a given tolerance
!                             (Generalized Binary Search)
! ClosestElement           Find index(es) in array closest to test value
!                           (array may be multidimensional, non-monotonic)
! Destroy                  Deallocate y values in LookUpTable_0
! d2Fdx2Approximate        Compute 2nd derivative using LookUpTable_0
! dFdxApproximate          Compute derivative using LookUpTable_0
! Dump                     Dump Coefficients structure
! Dump                     Dump uniform discrete function structure
! F_Of_X                   Use any recognized method to approximate f
! FApproximate             Use Uniformly Discretized Fn as an approximation
! FInvApproximate          Use it to invert a function (may not be unique)
! FillLookUpTable          Fill table with evaluations at regularly-spaced args
!                            to be used in place of later, frequent evaluations;
!                            reversing role of (table, xtable) => function^(-1)
! FindInRange              Finds indices of array elements for which values
!                            lie within a range
!                            (list need not be monotonic)
! Hunt                     Finds index of item(s) in list closest to prey
!                           (list must be monotonic)
! HuntBox                  Finds indices of n-dimensional box enclosing coords
! HuntRange                Finds index range between which
!                           all item(s) in list lie within range of values
! IFApproximate            Compute integral using LookUpTable_0
! InterpolateArraySetup    Compute coefficients for InterpolateUsingSetup
! InterpolateArrayTeardown Deallocate tables created by InterpolateArraySetup
! InterpolateExtrapolate   Like InterpolateValues, but extrapolate using
!                            average slope instead of slope at the end.
! Interpolate_Regular_To_Irregular 2D - to - 2D interpolation using composite
!                          1D interpolation
! InterpolateValues        Interpolate for new y value(s):
!                            given old (x,y), new (x), method
! Interpolate_2D_Composite 2D - to - 2D interpolation using composite 1D
!                          interpolation
! LinearInterpolate        Do a single linear interpolation in n dimensions
! PureHunt                 Like Hunt, but may be quicker due to optimization
! ReadLookUpTable          Read a LookUpTable from a text file
! Setup                    Fill y values in LookUpTable
! Simpsons                 Apply Simpson's rule to integrate--function form
! SimpsonsSub              Apply Simpson's rule to integrate--a subroutine
! UseLookUpTable           Use LookUpTable to approximate function
!                            (or its derivatives or its integral)
! WriteLookUpTable         Write a LookUpTable to a text file
! === (end of toc) ===

! === (start of api) ===
! nprec  AGM ( nprec x, nprec y )
! nprec  Average ( nprec A(:) )
! Battleship( int extern fun, int root, [int n1], [int maxPhase1], [int ns(:)], &
!    [int b], [char* options], [int status] )
! Battleship( log extern fun, int root, [int n1], [int maxPhase1], [int ns(:)], &
!    [log b], [char* options], [int status] )
! Battleship( nprec extern fun, nprec root, nprec arg1, nprec delta, &
!    [int maxPhase1], [int ns(:)], [int status] )
! BivariateLinearInterpolation ( real X_Basis(:), real Y_Basis(:),
!    real Table_2d(:,:), real X_Grid(:), real Y_Grid(:), real Out(:) )
! ClosestElement ( nprec test, nprec array, int indices, [char* options] )
! Destroy ( LookUpTable_0_nprec LkUpTable )
! nprec  dFdxApproximate ( nprec x, LookUpTable_0_nprec LkUpTable )
! nprec  d2Fdx2Approximate ( nprec x, LookUpTable_0_nprec LkUpTable )
! Dump ( coefficients_nprec Coeffs )
! Dump ( LookUpTable_0_nprec LkUpTable, [int Details] )
! nprec  F_Of_X
!         (see FApproximate and UseLookupTable)
! nprec  FApproximate ( nprec x, LookUpTable_0_nprec LkUpTable )
! nprec  FInvApproximate ( nprec y, LookUpTable_0_nprec LkUpTable, &
!     [nprec xS], [nprec xE] )
! FillLookUpTable ( nprec extern fun, nprec table(:), nprec x1, nprec x2, &
!   [int N], [nprec xtable(:)] )
! FindInRange ( num list(:), num vrange(2), int which(:), [how_many], [options] )
! FindInRange_2d ( num list(:,:), num vrange(2), int which(:,:), [how_many], [options] )
! Hunt ( nprec list, nprec values, int indices(:), &
!   [int start], [log allowTopValue], [log allowBelowValue], &
!   [log nearest], [log logSpace], [log fail] )
! HuntBox ( nprec gridPoints, int MGridPoints(:), nprec coords, int indices(:), &
!   [nprec vertices] )
! HuntRange ( num list(:), num vrange(2), int irange(2), options )
! nprec  IFApproximate ( LookUpTable_0_nprec LkUpTable, &
!     [nprec xS], [nprec xE] )
! InterpolateValues ( nprec oldX(:), nprec oldY(:), &
!   nprec newX(:), nprec newY(:), char* method, [char* extrapolate], &
!   [nprec badValue], [nprec rangeofperiod(2)], [log missingRegions], &
!   [nprec dyByDx(:), log skipNewY], [nprec IntYdX(:)] )
! InterpolateValues ( nprec oldX(:), nprec oldY(:,:), &
!   nprec newX(:), nprec newY(:,:), char* method, [char* extrapolate], &
!   [nprec badValue], [log missingRegions], [nprec dyByDx(:,:)], &
!   [MatrixElement_T dNewByDOld], [log skipNewY], [nprec IntYdX(:,:)] )
! Interpolate_2d_Composite ( XOld, YOld, ZOld, XNew, YNew, ZNew, &
!   XMethod, YMethod, [XExtrapolate], [YExtrapolate] )
! Interpolate_Regular_To_Irregular ( XOld, YOld, ZOld, &
!   XYNew, ZNew, XMethod, YMethod, [XExtrapolate], [YExtrapolate] )
! InterpolateArraySetup ( nprec OldX(:), nprec NewX(:), char* Method, &
!   coefficients_nprec Coeffs, &[log Extrapolate], [int Width], [log DyByDx],
!   [matrixElement_T dNewByDOld], [log IntYdX] )
! InterpolateArrayTeardown ( Coefficients_nprec Coeffs )
! nprec LinearInterpolate ( nprec values(:,:,..,:), nprec coords(:), &
!    nprec verts(:,:,..,:) )
! PureHunt ( nprec element, nprec array, int n, int jlo, int jhi )
! ReadLookUpTable ( char* filename, int n, nprec x(:), nprec y(:) )
! Setup ( LookUpTable_0_nprec LkUpTable, int N, nprec x1, nprec x2, [ nprec y(:)], &
!    [char* BC], [nprec yLeft], [nprec yRight], [extern nprec fun] )
! nprec Simpsons ( int n, nprec h, nprec y(:) )
! SimpsonsSub ( nprec y(:), nprec h, int n, nprec r )
! nprec UseLookUpTable ( nprec x, nprec table(:), [nprec x1], [nprec x2], &
!    [nprec xtable(:), [nprec missingValue], [char* options], &
!    [nprec xS], [nprec xE] )
! WriteLookUpTable ( char* filename, int n, nprec x(:), nprec y(:) )

! In the above types, "nprec" can be either r4 or r8. num can be any of
! int, r4, or r8. A, B, Precision, and array can be any 
! multidimensional arrays up to rank 3. In a scalar context A and B may be scalars
! === (end of api) ===

  public :: AGM, Average, Battleship, BivariateLinearInterpolation
  public :: ClosestElement
  public :: Destroy, dFdXApproximate, d2FdX2Approximate, Dump
  public :: F_Of_X, FApproximate, FindInRange, FinvApproximate
  public :: Hunt, HuntBox, HuntRange, IfApproximate
  public :: InterpolateArraySetup, InterpolateArrayTeardown
  public :: InterpolateExtrapolate
  public :: InterpolateExtrapolate_d, InterpolateExtrapolate_s
  public :: InterpolateScalar_r4, InterpolateScalar_r8
  public :: InterpolateValues
  public :: Interpolate_Regular_To_Irregular, Interpolate_2D_Composite
  public :: LinearInterpolate
  public :: PureHunt
  public :: Setup, Simpsons, SimpsonsSub, SolveQuadratic
  public :: FillLookupTable, ReadLookupTable, UseLookupTable, WriteLookupTable

  type, public :: Coefficients ( RK )
    integer, kind :: RK
    private
    !{ Coefficients for linear interpolation:  Let $\{x\} =$ {\tt oldX} and
    !  $\{\chi\} =$ {\tt newX}.  Then
    !  {\tt Gap(j)} $= g_j = x_{i+1}-x_i$, $A_j = \frac{x_{i+1}-\chi_j}{g_j}$,    
    !                                  and $B_j = 1-A_j = \frac{\chi_j-x_i}{g_j}$,
    !  where $i$ is such that $x_i \leq \chi_j < x_{i+1}$ for $1 \leq j \leq$
    !  {\tt size(newX)}.
    !  {\tt lowerInds} = $\{i\, |\, x_i \leq \chi_j < x_{i+1},\,1 \leq j \leq$
    !  {\tt size(newX)}$\}$.
    !  $|\{g\}| = |\{A\}| = |\{B\}|$, and
    !  others that depend upon them, $=|\{\chi\}|$ = {\tt size(newX)}.
    !  Coefficients for differentiation in the linear case (and of the linear
    !  terms in the spline case) are just $-1$ and $+1$, so they're not
    !  computed here.
    integer, allocatable :: LowerInds(:)
    real(rk), allocatable :: A(:), B(:), Gap(:)
    !{ Coefficients for spline interpolation:
    !  $C = (A^3-A) \frac{g^2}6$.  $D = (B^3-B) \frac{g^2}6$.
    real(rk), allocatable :: C(:), D(:)
    !{ {\tt dX(i)} $= \delta_i = x_i-x_{i-1}$ for $1 < i \leq$ {\tt size(oldX)} and
    !  $\Delta_i = x_{i+1}-x_{i-1} = \delta_i + \delta_{i+1}$ for $1 < i <$
    !  {\tt size(oldX)}.
    !  {\tt P(i)} $= p_i = \frac1{\delta_i o_{i-1} + 2 \Delta_i}$ and
    !  {\tt O(i)} $= o_i = -\delta_{i+1} p_i$
    !  for $1 < i <$ {\tt size(oldX)} and zero at the ends.
    !  $p_i$ is the inverse of the diagonal of $L$ in the $LU$ factorization
    !  of the symmetric tridiagonal linear system for splines, $\delta_i$ is
    !  the subdiagonal of $L$, and $o_i$ is the negative of the
    !  superdiagonal of $U$.  The nonzero part of each row of the original
    !  system is of the form $[\delta_{i-1},\, 2 \Delta_i ,\, \delta_i]$.
    !  See wvs-086.
    real(rk), allocatable :: dX(:), P(:), O(:)
    !{ In the periodic continuity case, {\tt Col} is the rightmost column of
    !  $U$ and {\tt Row} is the bottom row of $L$.  During the backsolve, the
    !  last element of {\tt Col} is used instead of the last element of $o$,
    !  and the last element of {\tt Row} is used instead of the last element
    !  of $\delta$.
    real(rk), allocatable :: Col(:), Row(:)
    !{ Coefficients for spline derivatives:
    !  $E = \frac{\text{d}C}{\text{d}A} = \frac6g \frac{\text{d}C}{\text{d}x}
    !     = 3 A^2 - 1$.
    !  $F = \frac{\text{d}D}{\text{d}A} = \frac6g \frac{\text{d}D}{\text{d}x}
    !     = 3 B^2 - 1$.
    real(rk), allocatable :: E(:), F(:)
    !{ Coefficients for integration:
    !  $\int A \text{d}x =  \frac{x(x_{i+1}-\frac{x}2)}g$.
    !  $\int B \text{d}x = -\frac{x(x_i    -\frac{x}2)}g
    !                    = x - \int A \text{d}x$.\\
    !  $\int C \text{d}x = -\frac{g^2}6
    !                      \left( \int A \text{d}x + \frac{g A^4}4 \right)$.
    !  $\int D \text{d}x = -\frac{g^2}6
    !                      \left( \int B \text{d}x - \frac{g B^4}4 \right)$.
    real(rk), allocatable :: AI(:), BI(:), &
      &                  CI(:), DI(:)
    ! Stuff for extrapolation == "B"ad
    logical, allocatable :: BadValue(:)
  end type Coefficients

  ! This is a family of datatypes, the Look Up Table:
  ! a list of function values y[i]
  ! evaluated over uniformly-spaced x-values x[i] ( = x1 + (i-1) dx )
  ! which some call a look up table
  ! (although that name unfairly limits its suggested use)
  ! In our usage it is an efficient alternative to the *LookupTable procedures
  ! implemented in this module, for they are not limited by requiring
  ! that the x-values be uniformly spaced (or even monotonic)
  
  ! By default when we use one of these datatypes to approximate the
  ! actual function at some x we 
  ! choose the corresponding y at the xi closest to x
  ! This behavior may be modified by the "method" field which may be one of
  ! ' '   choose whichever xi is closest to x          (default)
  ! 'l'   choose the lower of two xis closest to x
  ! 'u'   choose the upper of two xis closest to x
  ! 'i'   interpolate between two xis closest to x
  ! 'q'   quadratic interpolation (Simpson's rule)
  ! 's'   (cubic) splines (not efficiently implemented)

  ! Note that points outside the range of xi (i.e. less than the smallest
  ! or greater than the largest) are treated according to "BC" which may be
  ! 'clamped'   use y at the nearest xi                (default)                        
  ! 'cyclic'    assume y is cyclic with period (x2-x1) (how would we know?)          
  ! 'scyclic'   assume y is signed cyclic 
  !                    with half period (x2-x1)        (same question)          
  ! 'express'   use yLeft or yRight                    (of doubtful utility)                

  ! In cases where a conflict seems to arise between the actions
  ! dictated by "method" and "BC", "BC" will prevail

  ! Notes and limitations:
  ! Another family of datatypes may someday be needed,
  !    which would utilize non-uniformly spaced x values
  ! Another BC type may someday be needed,
  !    which would permit extrapolating x values outside [x1, x2] range

  type, public :: LookUpTable_0 ( RK )
    integer, kind :: RK
    integer :: N = 0                               ! The number of values xi
    character(len=8) :: BC = 'clamped'             ! boundary conditions
    character(len=1) :: method = ' '               ! which of closest xi to use
    real(rk) :: x1                                 ! x1 <= xi <= x2
    real(rk) :: x2
    real(rk) :: yLeft  = 0.                        ! Assume all x < x1 are this
    real(rk) :: yRight = 0.                        ! Assume all x > x2 are this
    real(rk), dimension(:), allocatable :: y       ! y(xi)
    ! type(Coefficients(rk)) :: Coeffs             ! in case we'll use splines
  contains
    procedure :: DestroyLookUpTable_0_r4
    procedure :: DestroyLookUpTable_0_r8
    generic   :: Destroy => destroyLookUpTable_0_r4, destroyLookUpTable_0_r8
    procedure :: DumpLookUpTable_0_r4
    procedure :: DumpLookUpTable_0_r8
    generic   :: Dump => DumpLookUpTable_0_r4, DumpLookUpTable_0_r8
    procedure :: setUpLookUpTable_0_r4
    procedure :: setUpLookUpTable_0_r8
    generic   :: Setup => setUpLookUpTable_0_r4, setUpLookUpTable_0_r8
  end type LookUpTable_0

  ! No need yet for a class member of this datatype yet
  ! type(LookUpTable_0(r8)), save :: MLSLkUpTable

  ! This data type
  ! (1) Relaxes the requirement that the x values are uniformly spaced
  !     although they must be sorted from smallest to largest
  ! (2) Includes optional polynomial interpolation via
  ! y[x] = a_0 + a_1 x + a_2 x^2 + .. + a_P x^P
  ! which is used only if P > 1
  
  ! However, we must decide, if N > P, among options (a) or (b)
  ! (a) choose only the P closest x values to uniquely determine a
  ! each time we are givin a new x
  ! (b) do a least squares over all y and x to cakculate a
  ! While similar in spirit to how we might imagine generalizing linear
  ! interpolation, (a) would be costly and sems to offer little value
  ! which (b) says we can calculate a just once and then go to town
  type, public, extends(LookUpTable_0) :: LookUpTable_1
    integer :: P = 0                               ! The degree of the polynomial
    real(rk), dimension(:), allocatable :: x       ! the x values
    real(rk), dimension(:), allocatable :: a       ! the a values
  end type LookUpTable_1


  interface AGM
    module procedure AGM_D, AGM_S
  end interface

  interface Average
    module procedure Average_D, Average_S
  end interface

  interface Battleship
    module procedure Battleship_int, Battleship_log, Battleship_r4, Battleship_r8
  end interface

  interface BivariateLinearInterpolation
    module procedure BivariateLinearInterp_D_D, BivariateLinearInterp_D_S
    module procedure BivariateLinearInterp_S_S
  end interface

  interface ClosestElement
    module procedure ClosestElement_r4_1d, ClosestElement_r8_1d
    module procedure ClosestElement_r4_2d, ClosestElement_r8_2d
    module procedure ClosestElement_r4_3d, ClosestElement_r8_3d
  end interface

  interface CreateXArray
    module procedure CreateXArray_r4, CreateXArray_r8
  end interface

  interface CSpline
    module procedure D_CSpline, S_CSpline
  end interface

  interface Destroy
    module procedure destroyLookUpTable_0_r4, destroyLookUpTable_0_r8
  end interface

  interface d2Fdx2Approximate
    module procedure d2Fdx2Approximate_r4, d2Fdx2Approximate_r8
  end interface

  interface dFdxApproximate
    module procedure dFdxApproximate_r4, dFdxApproximate_r8
  end interface

  interface Dump
    module procedure DumpCoefficients_r4, DumpCoefficients_r8
    module procedure DumpLookUpTable_0_r4, DumpLookUpTable_0_r8
  end interface

  interface F_Of_X
    module procedure FApproximate_r4, FApproximate_r8
    module procedure FLookup_r4, FLookup_r8
    module procedure FXYLookup_r4, FXYLookup_r8
  end interface

  interface FApproximate
    module procedure FApproximate_r4, FApproximate_r8
  end interface

  interface FInvApproximate
    module procedure FInvApproximate_r4, FInvApproximate_r8
  end interface

  interface FillLookUpTable
    module procedure FillLookUpTable_r4, FillLookUpTable_r8
  end interface

  interface FindInRange
    module procedure FindInRange_int, FindInRange_r4, FindInRange_r8
    module procedure FindInRange_2d_int, FindInRange_2d_r4, FindInRange_2d_r8
  end interface

  interface IFApproximate
    module procedure IFApproximate_r4, IFApproximate_r8
  end interface

  interface InterpolateArraySetup
    module procedure InterpolateArraySetup_r4, InterpolateArraySetup_r8
  end interface

  interface InterpolateExtrapolate
    module procedure InterpolateExtrapolate_d, InterpolateExtrapolate_d_1
    module procedure InterpolateExtrapolate_s, InterpolateExtrapolate_s_1
  end interface

  interface Interpolate_Regular_To_Irregular
    module procedure Interpolate_Regular_To_Irregular_r4, &
                     Interpolate_Regular_To_Irregular_r8
  end interface

  interface InterpolateArrayTeardown
    module procedure InterpolateArrayTeardown_r4, InterpolateArrayTeardown_r8
  end interface

  interface InterpolateValues
    module procedure InterpolateArray_r4, InterpolateArray_r8
    module procedure InterpolateScalar_r4, InterpolateScalar_r8
    module procedure InterpolateUsingSetup_r4, InterpolateUsingSetup_r8
    module procedure InterpolateScalarUsingSetup_r4, InterpolateScalarUsingSetup_r8
    module procedure Interp_Bilinear_2d_1d_r4, Interp_Bilinear_2d_1d_r8
  end interface

  interface Interpolate_2d_Composite
    module procedure Interpolate_2d_Composite_r4, Interpolate_2d_Composite_r8
  end interface

  interface LinearInterpolate
    module procedure LinearInterpolate_1d_r4, LinearInterpolate_1d_r8
    module procedure LinearInterpolate_2d_r4, LinearInterpolate_2d_r8
    module procedure LinearInterpolate_3d_r4, LinearInterpolate_3d_r8
    module procedure LinearInterpolate_4d_r4, LinearInterpolate_4d_r8
  end interface

  interface pcspl
    module procedure D_PCSPL, S_PCSPL
  end interface

  interface psimpsons
    module procedure psimpsons_r4, psimpsons_r8
  end interface

  interface reposit
    module procedure reposit_r4, reposit_r8
  end interface

  interface ReadLookUpTable
    module procedure ReadLookUpTable_r4, ReadLookUpTable_r8
  end interface
  
  interface SetUp
    module procedure InterpolateArraySetup_r4, InterpolateArraySetup_r8
    module procedure setUpLookUpTable_0_r4, setUpLookUpTable_0_r8
  end interface

  interface Simpsons
    module procedure Simpsons_r4, Simpsons_r8
  end interface

  interface SimpsonsSub
    module procedure Simps_r4, Simps_r8
  end interface

  interface SolveQuadratic
    module procedure SolveQuadratic_r4, SolveQuadratic_r8
  end interface

  interface UseLookUpTable
    module procedure UseLookUpTable_r4, UseLookUpTable_r8
  end interface
  
  interface WriteLookUpTable
    module procedure WriteLookUpTable_r4, WriteLookUpTable_r8
  end interface

  ! These are arrays in name only used when implementing cubic splines
  real, private, dimension(1) :: newYr4, newdYr4
  double precision, private, dimension(1) :: newYr8, newdYr8
  
  ! Parameters used in the AGM functions
  integer, parameter                      :: MAXITER = 20

contains

! -------------------------------------------------  AGM  -----
! Return the arithmetic-geometric mean
! Given two numbers, form the sequences
! a_n, g_n
! a[0] = x
! g[0] = y
! a[i+1] = (a[i]+g[i])/2
! g[i+1] = sqrt(a[i]*g[i])
! Find the limit as i -> Inf
  double precision function AGM_D ( x, y ) result ( R )
    double precision, intent(in) :: x, y
    ! Internal variables
    integer                      :: i
    double precision             :: a0, g0, ai, gi, eps
    include 'AGM.f9h'
  end function AGM_D

  real function AGM_S ( x, y ) result ( R )
    real, intent(in)             :: x, y
    ! Internal variables
    integer                      :: i
    real                         :: a0, g0, ai, gi, eps
    include 'AGM.f9h'
  end function AGM_S

! -------------------------------------------------  Average  -----
! Return the arithmetic mean
! If you want the geometric mean instead, you can use this formula
!   geom_mean = exp ( Average( Log (a) ) )
! Assuming all a > 0
  double precision function Average_D ( A ) result ( R )
    double precision, intent(in) :: A(:)
    r = sum(a) / size(a)
  end function Average_D

  real function Average_S ( A ) result ( R )
    real, intent(in) :: A(:)
    r = sum(a) / size(a)
  end function Average_S

! -------------------------------------------------  Battleship  -----

  ! This family of routines finds a root of a function
  ! by repeatedly evaluating it. Each evaluation is a "shot". What we consider
  ! a "hit" depends on the options parameter (see below).
  ! Warning--the default behavior is not the usual root-finder's
  ! (+ve on one side, -ve on the other)
  ! Instead what we do is this:
  !
  ! A returned value
  ! of "0" (or the optional parameter b) is short. Any other value is long.
  ! The root is the longest argument that is still short.
  ! (If you instead wish the shortest argument still long just add 1
  ! or utilize the options string)

  ! Example: a direct-access read of n chars from a file where the iostatus
  ! is 0 if we don't try to read too many chars, but non-zero if we do
  ! Create your own function that takes the number of chars to be read
  ! as its sole argument and that returns the iostat as its value
  ! Battleship will then calculate the exact number characters in the file
  ! (which NAG cares about; Lahey doesn't care)

  ! If you wish to do the usual root-finder where the returned values change
  ! sign, set the options appropriately to "-x"

  ! Method:
  ! We take shots during 2 phases:
  ! (1) outbound: ever-widening circles of radius 1 2 4 8 .. (n) (2n) ..
  !     or else prescribed shots in array ns[:]
  ! (2) inbound: once root is crossed, ever narowing circles around it
  !     (until "You sank my battleship!")

  ! The options string (if supplied) modifies how this search operates
  !  options            search goal
  !  -------            -----------
  !    -s (default)     largest root for which fun(root) = 0 (or b)
  !                      (useful for io status)
  !    -r               reverse of "-s"
  !                       i.e., all tests return non-zero below root
  !    -x               root where f(root) crosses 0 (or b)
  !                      (assumes (fun(n)-b) changes sign at n=root)

  ! It can also be used with a logical-valued function
  ! in this case we shoot until we encounter TRUE
  ! similar to integer-value version if we make the mapping
  ! 0 -> FALSE
  ! 1 -> TRUE
  ! -r option or b can be used to reverse the sense

  ! A version for use with real-valued functions can be
  ! operated as a root-finder where the root is
  ! desired to be found within a tolerance of delta
  ! The real-valued version won't be as precise or as efficient
  ! as, say, the Zero subroutine in the Zero_m module
  
  ! A generalization of the Binary Search Algorithm to
  ! two phases (outbound to first bracket the target, then inbound)

  subroutine Battleship_int( fun, root, n1, maxPhase1, ns, b, options, status )
    ! Args
    integer, external                          :: fun
    integer, optional, intent(in)              :: n1 ! 1st circle
    integer, optional, intent(in)              :: maxPhase1 ! max phase1 shots
    integer, optional, dimension(:), intent(in):: ns ! array of phase1 shots
    integer, intent(out)                       :: root ! root
    integer, optional, intent(in)              :: b ! is short
    character(len=*), optional, intent(in)     :: options
    integer, optional, intent(out)             :: status ! /= 0 if failed
    ! Internal variables
    integer :: flast
    integer :: fnext
    integer :: isShort
    character(len=8) :: myOptions
    integer :: shot
    integer :: x0
    integer :: x1
    integer :: x2
    ! Executable
    isShort = 0
    if ( present(b) ) isShort = b
    myOptions = '-s'
    if ( present(options) ) myOptions = options
    root = -1 ! in case we can't find root
    if ( present(status) ) status = 1
    ! Phase 1
    ! Some error checks
    if ( present(maxPhase1) ) then
      if ( (index(myOptions, 's') > 0 .and. fun(n1) /= isShort) ) return
      if ( (index(myOptions, 'r') > 0 .and. fun(n1) == isShort) ) return
      if ( maxPhase1 < 1 ) return
      if ( .not. present(n1) ) return
      x2 = n1 ! Initialize things
      fnext = fun(x2)
      do shot = 1, maxPhase1
        x1 = x2
        x2 = 2*x1
        flast = fnext
        fnext = fun(x2)
        if ( index(myOptions, 's') > 0 ) then
          if ( fnext /= isShort ) exit
        elseif ( index(myOptions, 'r') > 0 ) then
          if ( fnext == isShort ) exit
        else
          if ( (fnext-isShort)*(flast-isShort) <= 0 ) exit
        end if
      enddo
      if ( shot > maxPhase1 ) return ! No shot was long enough
    else
      if ( .not. present(ns) ) return
      x2 = ns(1) ! Initialize things
      fnext = fun(x2)
      if ( (index(myOptions, 's') > 0 .and. fnext /= isShort) ) return
      if ( (index(myOptions, 'r') > 0 .and. fnext == isShort) ) return
      do shot = 2, size(ns)
        x1 = x2
        x2 = ns(shot)
        flast = fnext
        fnext = fun(x2)
        if ( index(myOptions, 's') > 0 ) then
          if ( fnext /= isShort ) exit
        elseif ( index(myOptions, 'r') > 0 ) then
          if ( fnext == isShort ) exit
        else
          if ( (fnext-isShort)*(flast-isShort) <= 0 ) exit
        end if
      enddo
      if ( shot > size(ns) ) return ! No shot was long enough
    end if
    ! Phase 2
    if ( present(status) ) status = 0
    ! Narrow the spashes, always keeping root between x0 and x2
    x0 = x1
    do
      x1 = (x0 + x2) / 2
      ! This test should prevent us from looping endlessly
      if ( x1 == x0 ) then
        ! apparently x0 = x2 - 1, so we've found our root
        if ( index(myOptions, 's') > 0 .or. index(myOptions, 'r') > 0 ) then
          root = x1
        else
          if ( fun(x1) == isShort ) then
            root = x1
          elseif ( fun(x2) == isShort ) then
            root = x2
          else
            root = -1
          end if
        end if
        return
      end if
      if ( index(myOptions, 's') > 0 ) then
        if ( fun(x1) == isShort ) then
          x0 = x1
          ! x2 = x2
        else
          ! x0 = x0
          x2 = x1
        end if
      elseif ( index(myOptions, 'r') > 0 ) then
        if ( fun(x1) == isShort ) then
          ! x0 = x0
          x2 = x1
        else
          x0 = x1
          ! x2 = x2
        end if
      else
        if ( (fun(x2)-isShort)*(fun(x1)-isShort) == 0 ) then
          if ( fun(x2) == isShort ) then
            root = x2
          else
            root = x1
          end if
          return
        elseif ( (fun(x2)-isShort)*(fun(x1)-isShort) < 0 ) then
          x0 = x1
        else
          x2 = x1
        end if
      end if
    enddo
  end subroutine Battleship_int

  subroutine Battleship_log( fun, root, n1, maxPhase1, ns, b, options, status )
    ! Args
    logical, external                          :: fun
    integer, optional, intent(in)              :: n1 ! 1st circle
    integer, optional, intent(in)              :: maxPhase1 ! max phase1 shots
    integer, optional, dimension(:), intent(in):: ns ! array of phase1 shots
    integer, intent(out)                       :: root ! root
    logical, optional, intent(in)              :: b ! is short
    character(len=*), optional, intent(in)     :: options
    integer, optional, intent(out)             :: status ! /= 0 if failed
    ! Internal variables
    logical :: flast
    logical :: fnext
    logical :: isShort
    character(len=8) :: myOptions
    integer :: shot
    integer :: x0
    integer :: x1
    integer :: x2
    ! Executable
    isShort = .FALSE.
    if ( present(b) ) isShort = b
    myOptions = '-s'
    if ( present(options) ) myOptions = options
    root = -1 ! in case we can't find root
    if ( present(status) ) status = 1
    ! Phase 1
    ! Some error checks
    if ( present(maxPhase1) ) then
      if ( (index(myOptions, 's') > 0 .and. ( fun(n1) .neqv. isShort) ) ) return
      if ( (index(myOptions, 'r') > 0 .and. ( fun(n1) .eqv. isShort) ) ) return
      if ( maxPhase1 < 1 ) return
      if ( .not. present(n1) ) return
      x2 = n1 ! Initialize things
      fnext = fun(x2)
      do shot = 1, maxPhase1
        x1 = x2
        x2 = 2*x1
        flast = fnext
        fnext = fun(x2)
        if ( index(myOptions, 's') > 0 ) then
          if ( fnext .neqv. isShort ) exit
        elseif ( index(myOptions, 'r') > 0 ) then
          if ( fnext .eqv. isShort ) exit
        else
          if ( fnext .neqv. flast ) exit
        end if
      enddo
      if ( shot > maxPhase1 ) return ! No shot was long enough
    else
      if ( .not. present(ns) ) return
      x2 = ns(1) ! Initialize things
      fnext = fun(x2)
      if ( (index(myOptions, 's') > 0 .and. ( fnext .neqv. isShort ) ) ) return
      if ( (index(myOptions, 'r') > 0 .and. ( fnext .eqv. isShort) ) ) return
      do shot = 2, size(ns)
        x1 = x2
        x2 = ns(shot)
        flast = fnext
        fnext = fun(x2)
        if ( index(myOptions, 's') > 0 ) then
          if ( fnext .neqv. isShort ) exit
        elseif ( index(myOptions, 'r') > 0 ) then
          if ( fnext .eqv. isShort ) exit
        else
          if ( fnext .neqv. flast ) exit
        end if
      end do
      if ( shot > size(ns) ) return ! No shot was long enough
    end if
    ! Phase 2
    if ( present(status) ) status = 0
    ! Narrow the spashes, always keeping root between x0 and x2
    x0 = x1
    do
      x1 = (x0 + x2) / 2
      ! This test should prevent us from looping endlessly
      if ( x1 == x0 ) then
        ! apparently x0 = x2 - 1, so we've found our root
        if ( index(myOptions, 's') > 0 .or. index(myOptions, 'r') > 0 ) then
          root = x1
        else
          if ( fun(x1) .eqv. isShort ) then
            root = x1
          elseif ( fun(x2) .eqv. isShort ) then
            root = x2
          else
            root = -1
          end if
        end if
        return
      end if
      if ( index(myOptions, 's') > 0 ) then
        if ( fun(x1) .eqv. isShort ) then
          x0 = x1
          ! x2 = x2
        else
          ! x0 = x0
          x2 = x1
        end if
      elseif ( index(myOptions, 'r') > 0 ) then
        if ( fun(x1) .eqv. isShort ) then
          ! x0 = x0
          x2 = x1
        else
          x0 = x1
          ! x2 = x2
        end if
      else
        if ( (fun(x2) .eqv. isShort) .or. (fun(x1) .eqv. isShort) ) then
          if ( fun(x2) .eqv. isShort ) then
            root = x2
          else
            root = x1
          end if
          return
        elseif ( fun(x2) .neqv. fun(x1) ) then
          x0 = x1
        else
          x2 = x1
        end if
      end if
    enddo
  end subroutine Battleship_log

  subroutine Battleship_r4( fun, root, arg1, delta, maxPhase1, ns, status )
    integer, parameter :: RK = kind(0.0e0)
      include "Battleship.f9h"
  end subroutine Battleship_r4

  subroutine Battleship_r8( fun, root, arg1, delta, maxPhase1, ns, status )
    integer, parameter :: RK = kind(0.0d0)
      include "Battleship.f9h"
  end subroutine Battleship_r8

! ------------------------------------  BivariateLinearInterp_D_D  -----
  subroutine BivariateLinearInterp_D_D ( X_Basis, Y_Basis, Table_2d, &
      & X_Grid, Y_Grid, Out )

      ! Interpolate linearly in Table_2d whose coordinates are (X_Basis,Y_Basis),
      ! which are assumed to be sorted, to Out at each (X_Grid,Y_Grid), which
      ! are not necessarily sorted.  The Hunt routine assumes X_Grid and Y_Grid
      ! change slowly and smoothly.

      integer, parameter :: KT = kind(0.0d0) ! Kind for table and basis
      integer, parameter :: KO = kind(0.0d0) ! Kind for out and grids

      ! Extents for Table_2d are (size(X_basis,1),size(Y_Basis,1))
      real(kt), intent(in) :: X_Basis(:), Y_Basis(:), Table_2D(:,:)

      ! Extents for X_Grid, Y_Grid, Out are all the same.
      real(ko), intent(in) :: X_Grid(:), Y_Grid(:)
      real(ko), intent(out) :: Out(:)

      include "BivariateLinearInterpolation.f9h"

    end subroutine BivariateLinearInterp_D_D

! ------------------------------------  BivariateLinearInterp_D_S  -----
  subroutine BivariateLinearInterp_D_S ( X_Basis, Y_Basis, Table_2d, &
      & X_Grid, Y_Grid, Out )

      ! Interpolate linearly in Table_2d whose coordinates are (X_Basis,Y_Basis),
      ! which are assumed to be sorted, to Out at each (X_Grid,Y_Grid), which
      ! are not necessarily sorted.  The Hunt routine assumes X_Grid and Y_Grid
      ! change slowly and smoothly.

      integer, parameter :: KT = kind(0.0d0) ! Kind for table and basis
      integer, parameter :: KO = kind(0.0e0) ! Kind for out and grids

      ! Extents for Table_2d are (size(X_basis,1),size(Y_Basis,1))
      real(kt), intent(in) :: X_Basis(:), Y_Basis(:), Table_2D(:,:)

      ! Extents for X_Grid, Y_Grid, Out are all the same.
      real(ko), intent(in) :: X_Grid(:), Y_Grid(:)
      real(ko), intent(out) :: Out(:)

      include "BivariateLinearInterpolation.f9h"

    end subroutine BivariateLinearInterp_D_S

! ------------------------------------  BivariateLinearInterp_S_S  -----
  subroutine BivariateLinearInterp_S_S ( X_Basis, Y_Basis, Table_2d, &
      & X_Grid, Y_Grid, Out )

      ! Interpolate linearly in Table_2d whose coordinates are (X_Basis,Y_Basis),
      ! which are assumed to be sorted, to Out at each (X_Grid,Y_Grid), which
      ! are not necessarily sorted.  The Hunt routine assumes X_Grid and Y_Grid
      ! change slowly and smoothly.

      integer, parameter :: KT = kind(0.0e0) ! Kind for table and basis
      integer, parameter :: KO = kind(0.0e0) ! Kind for out and grids

      ! Extents for Table_2d are (size(X_basis,1),size(Y_Basis,1))
      real(kt), intent(in) :: X_Basis(:), Y_Basis(:), Table_2D(:,:)

      ! Extents for X_Grid, Y_Grid, Out are all the same.
      real(ko), intent(in) :: X_Grid(:), Y_Grid(:)
      real(ko), intent(out) :: Out(:)

      include "BivariateLinearInterpolation.f9h"

    end subroutine BivariateLinearInterp_S_S

! -----------------------------------------------  ClosestElement  -----

  ! This family of routines finds the element within a multidimensional
  ! array nearest a test value
  ! The array of indices locate that nearest element
  ! The following options really make sense only for 1-d searches
  ! and so are ignored for arrays of rank 2 or higher

  ! options  none, one, or more of the following:
  ! (default)   choose pt in array closest to x
  !   l         always choose lower of two closest x's in array
  !   u         always choose upper of two closest x's in array
  !   p         assume array is presorted so array[i] < array[i+1]
  !              (greatly speeds up search)

  subroutine ClosestElement_r4_1d ( test, array, indices, options )
    integer, parameter :: RK = kind(0.0e0)

    ! Dummy arguments
    real(rk), intent(in)               :: test
    real(rk), dimension(:), intent(in) :: array
    integer, dimension(:), intent(out) :: indices ! Result
    character(len=*), optional, intent(in)      :: options
    include "ClosestElement.f9h"

  end subroutine ClosestElement_r4_1d

  subroutine ClosestElement_r8_1d ( test, array, indices, options )
    integer, parameter :: RK = R8

    ! Dummy arguments
    real(rk), intent(in)               :: test
    real(rk), dimension(:), intent(in) :: array
    integer, dimension(:), intent(out) :: indices ! Result
    character(len=*), optional, intent(in)      :: options
    include "ClosestElement.f9h"

  end subroutine ClosestElement_r8_1d

  subroutine ClosestElement_r4_2d ( test, array, indices, options )
    integer, parameter :: RK = R4

    ! Dummy arguments
    real(rk), intent(in)               :: test
    real(rk), dimension(:,:), intent(in) :: array
    integer, dimension(:), intent(out) :: indices ! Result
    character(len=*), optional, intent(in)      :: options
    integer, dimension(1)              :: indices_1d ! Result
    call ClosestElement( test, &
      & reshape(array, (/ size(array,1)*size(array,2) /) ), &
      & indices_1d )
    call rerank( indices_1d(1), shape(array), indices )
  end subroutine ClosestElement_r4_2d

  subroutine ClosestElement_r8_2d ( test, array, indices, options )
    integer, parameter :: RK = kind(0.0d0)

    ! Dummy arguments
    real(rk), intent(in)               :: test
    real(rk), dimension(:,:), intent(in) :: array
    integer, dimension(:), intent(out) :: indices ! Result
    character(len=*), optional, intent(in)      :: options
    integer, dimension(1)              :: indices_1d ! Result
    call ClosestElement( test, &
      & reshape(array, (/ size(array,1)*size(array,2) /) ), &
      & indices_1d )
    call rerank( indices_1d(1), shape(array), indices )
  end subroutine ClosestElement_r8_2d

  subroutine ClosestElement_r4_3d ( test, array, indices, options )
    integer, parameter :: RK = kind(0.0e0)

    ! Dummy arguments
    real(rk), intent(in)               :: test
    real(rk), dimension(:,:,:), intent(in) :: array
    integer, dimension(:), intent(out) :: indices ! Result
    character(len=*), optional, intent(in)      :: options
    integer, dimension(1)              :: indices_1d ! Result
    call ClosestElement( test, &
      & reshape(array, (/ size(array,1)*size(array,2)*size(array,3) /) ), &
      & indices_1d )
    call rerank( indices_1d(1), shape(array), indices )
  end subroutine ClosestElement_r4_3d

  subroutine ClosestElement_r8_3d ( test, array, indices, options )
    integer, parameter :: RK = kind(0.0d0)

    ! Dummy arguments
    real(rk), intent(in)               :: test
    real(rk), dimension(:,:,:), intent(in) :: array
    integer, dimension(:), intent(out) :: indices ! Result
    character(len=*), optional, intent(in)      :: options
    integer, dimension(1)              :: indices_1d ! Result
    call ClosestElement( test, &
      & reshape(array, (/ size(array,1)*size(array,2)*size(array,3) /) ), &
      & indices_1d )
    call rerank( indices_1d(1), shape(array), indices )
  end subroutine ClosestElement_r8_3d

! ------------------------------------------------------  Destroy  -----

  ! This family of routines deallocates a LookUpTable's arrays
  ! We ought to do this for both _0 and _1 levels
  subroutine destroyLookUpTable_0_r4 ( LkUpTable )
    ! Args
    class(LookUpTable_0(r4)) :: LkUpTable ! Intent(out) would clobber retainable values
    ! Executable
    LkUpTable%N       = 0
    call deallocate_test ( LkUpTable%y, "LkUpTable%y", ModuleName )
  end subroutine destroyLookUpTable_0_r4

  subroutine destroyLookUpTable_0_r8 ( LkUpTable )
    ! Args
    class(LookUpTable_0(r8)) :: LkUpTable ! Intent(out) would clobber retainable values
    ! Executable
    LkUpTable%N       = 0
    call deallocate_test ( LkUpTable%y, "LkUpTable%y", ModuleName )
  end subroutine destroyLookUpTable_0_r8

! --------------------------------------------  d2Fdx2Approximate  -----

  ! This family of routines use a LookUpTable to approximate a 
  ! function's 2nd derivative

  ! Args: (* means optional)
  ! x        pt at which to evaluate
  ! LkUpTable      the Look Up Table type

  function d2Fdx2Approximate_r4 ( x, LkUpTable ) result(value)
    integer, parameter :: RK = kind(0.0e0)
    ! Args
    real(rk), intent(in)                 :: x
    type(LookUpTable_0(r4)), intent(in)       :: LkUpTable
    real(rk)                             :: value
    ! Internal variables
    real(rk) :: arg
    integer  :: itsSign
    character :: options
    real(rk), dimension(:), pointer      :: xArray => null()
    ! Executable
    call reposit( x, LkUpTable, arg, itsSign )
    options = LkUpTable%method
    if ( options == 'q' ) options = '2'
    if ( arg < LkUpTable%x1 ) then
      value = 0.
    elseif ( arg > LkUpTable%x2 ) then
      value = 0.
    elseif ( options == '/s/' ) then ! How to calculate spline's 2nd der?
      call createXArray( xArray, LkUpTable )
      call InterpolateValues( xArray, LkUpTable%y, (/arg/), newYr4, dyByDx=newdYr4, &
        & method='S' )
      value = newdYr4(1)
      call deallocate_test( xArray, 'xArray (r4)', ModuleName )
    else
      value = UseLookUpTable ( arg, LkUpTable%y, LkUpTable%x1, LkUpTable%x2, options=options // '2' )
    end if
    value = itsSign*value
  end function d2Fdx2Approximate_r4

  function d2Fdx2Approximate_r8 ( x, LkUpTable ) result(value)
    integer, parameter :: RK = kind(0.0d0)
    ! Args
    real(rk), intent(in)                 :: x
    type(LookUpTable_0(r8)), intent(in)       :: LkUpTable
    real(rk)                             :: value
    ! Internal variables
    real(rk) :: arg
    integer  :: itsSign
    character :: options
    real(rk), dimension(:), pointer      :: xArray => null()
    ! Executable
    call reposit( x, LkUpTable, arg, itsSign )
    options = LkUpTable%method
    if ( options == 'q' ) options = '2'
    if ( arg < LkUpTable%x1 ) then
      value = 0.
    elseif ( arg > LkUpTable%x2 ) then
      value = 0.
    elseif ( options == '/s/' ) then ! How to calculate spline's 2nd der?
      call createXArray( xArray, LkUpTable )
      call InterpolateValues( xArray, LkUpTable%y, (/arg/), newYr8, dyByDx=newdYr8, &
        & method='S' )
      value = newdYr8(1)
      call deallocate_test( xArray, 'xArray (r8)', ModuleName )
    else
      value = UseLookUpTable ( arg, LkUpTable%y, LkUpTable%x1, LkUpTable%x2, options=options // '2' )
    end if
    value = itsSign*value
  end function d2Fdx2Approximate_r8

! ----------------------------------------------  dFdxApproximate  -----

  ! This family of routines use a LookUpTable to approximate a 
  ! function's derivative

  ! Args: (* means optional)
  ! x        pt at which to evaluate
  ! LkUpTable      the Look Up Table type

  function dFdxApproximate_r4 ( x, LkUpTable ) result(value)
    integer, parameter :: RK = kind(0.0e0)
    ! Args
    real(rk), intent(in)                 :: x
    type(LookUpTable_0(r4)), intent(in)       :: LkUpTable
    real(rk)                             :: value
    ! Internal variables
    real(rk) :: arg
    integer  :: itsSign
    character :: options
    real(rk), dimension(:), pointer      :: xArray => null()
    ! Executable
    call reposit( x, LkUpTable, arg, itsSign )
    options = LkUpTable%method
    if ( options == 'q' ) options = '1'
    if ( arg < LkUpTable%x1 ) then
      value = 0.
    elseif ( arg > LkUpTable%x2 ) then
      value = 0.
    elseif ( options == 's' ) then
      call createXArray( xArray, LkUpTable )
      call InterpolateValues( xArray, LkUpTable%y, (/arg/), newYr4, dyByDx=newdYr4, &
        & method='S' )
      value = newdYr4(1)
      call deallocate_test( xArray, 'xArray (r4)', ModuleName )
    else
      value = UseLookUpTable ( arg, LkUpTable%y, LkUpTable%x1, LkUpTable%x2, options=options // '1' )
    end if
    value = itsSign*value
  end function dFdxApproximate_r4

  function dFdxApproximate_r8 ( x, LkUpTable ) result(value)
    integer, parameter :: RK = kind(0.0d0)
    ! Args
    real(rk), intent(in)                 :: x
    type(LookUpTable_0(r8)), intent(in)       :: LkUpTable
    real(rk)                             :: value
    ! Internal variables
    real(rk) :: arg
    integer  :: itsSign
    character :: options
    real(rk), dimension(:), pointer      :: xArray => null()
    ! Executable
    call reposit( x, LkUpTable, arg, itsSign )
    options = LkUpTable%method
    if ( options == 'q' ) options = '1'
    if ( arg < LkUpTable%x1 ) then
      value = 0.
    elseif ( arg > LkUpTable%x2 ) then
      value = 0.
    elseif ( options == 's' ) then
      call createXArray( xArray, LkUpTable )
      call InterpolateValues( xArray, LkUpTable%y, (/arg/), newYr8, dyByDx=newdYr8, &
        & method='S' )
      value = newdYr8(1)
      call deallocate_test( xArray, 'xArray (r8)', ModuleName )
    else
      value = UseLookUpTable ( arg, LkUpTable%y, LkUpTable%x1, LkUpTable%x2, options=options // '1' )
    end if
    value = itsSign*value
  end function dFdxApproximate_r8

! ---------------------------------------------------------  Dump  -----
  subroutine DumpCoefficients_r4 ( Coeffs, Name )
    type(coefficients(r4)), intent(in) :: Coeffs
    character(len=*), optional, intent(in) :: Name
    include 'DumpCoefficients.f9h'
  end subroutine DumpCoefficients_r4

  subroutine DumpCoefficients_r8 ( Coeffs, Name )
    type(coefficients(r8)), intent(in) :: Coeffs
    character(len=*), optional, intent(in) :: Name
    include 'DumpCoefficients.f9h'
  end subroutine DumpCoefficients_r8

  ! We ought to do this for both _0 and _1 levels
  subroutine DumpLookUpTable_0_r4 ( LkUpTable, name, details )
    ! Args
    class(LookUpTable_0(r4)) :: LkUpTable ! Intent(out) would clobber retainable values
    character(len=*), optional, intent(in) :: name
    integer, optional, intent(in) :: details
    ! Internal variables
    integer :: myDetails
    ! Executable
    myDetails = 0
    if ( present(details) ) myDetails = details
    if ( present(name) ) call output( trim(name) // ' ', advance='no' )
    call output( 'uniform discrete function (r4)', advance='yes' )
    call blanks( 32, fillchar='-', advance='yes' )
    call outputNamedValue( 'N           ', LkUpTable%N             )
    call outputNamedValue( 'x1          ', LkUpTable%x1            )
    call outputNamedValue( 'x2          ', LkUpTable%x2            )
    call outputNamedValue( 'BC          ', LkUpTable%BC            )
    call outputNamedValue( 'method      ', LkUpTable%method        )
    call outputNamedValue( 'yLeft       ', LkUpTable%yLeft         )
    call outputNamedValue( 'yRight      ', LkUpTable%yRight        )
    if ( .not. allocated(LkUpTable%y) ) then
      call output( '(y values not associated)', advance='yes' )
      return
    end if
    if ( myDetails < 1 ) then
      call outputNamedValue( 'min(y)           ', minval(LkUpTable%y) )
      call outputNamedValue( 'max(y)           ', maxval(LkUpTable%y) )
      return
    end if
    call dump ( LkUpTable%y, name='y' )
  end subroutine DumpLookUpTable_0_r4

  subroutine DumpLookUpTable_0_r8 ( LkUpTable, name, details )
    ! Args
    class(LookUpTable_0(r8)) :: LkUpTable ! Intent(out) would clobber retainable values
    character(len=*), optional, intent(in) :: name
    integer, optional, intent(in) :: details
    ! Internal variables
    integer :: myDetails
    ! Executable
    myDetails = 0
    if ( present(details) ) myDetails = details
    if ( present(name) ) call output( trim(name) // ' ', advance='no' )
    call output( 'uniform discrete function (r8)', advance='yes' )
    call blanks( 32, fillchar='-', advance='yes' )
    call outputNamedValue( 'N           ', LkUpTable%N             )
    call outputNamedValue( 'x1          ', LkUpTable%x1            )
    call outputNamedValue( 'x2          ', LkUpTable%x2            )
    call outputNamedValue( 'BC          ', LkUpTable%BC            )
    call outputNamedValue( 'method      ', LkUpTable%method        )
    call outputNamedValue( 'yLeft       ', LkUpTable%yLeft         )
    call outputNamedValue( 'yRight      ', LkUpTable%yRight        )
    if ( .not. allocated(LkUpTable%y) ) then
      call output( '(y values not associated)', advance='yes' )
      return
    end if
    if ( myDetails < 1 ) then
      call outputNamedValue( 'min(y)           ', minval(LkUpTable%y) )
      call outputNamedValue( 'max(y)           ', maxval(LkUpTable%y) )
      return
    end if
    call dump ( LkUpTable%y, name='y' )
  end subroutine DumpLookUpTable_0_r8
! -------------------------------------------------  FLookup_r4  -----

  ! This family of routines uses lookup tables to approximate F

  ! Args: (* means optional)
  ! x        pt at which to evaluate
  ! x1, x2   (if present) lower, upper bounds of x range
  ! xtable   (if present) x values (should be monotonic)
  ! ytable   y values

  function FLookup_r4 ( x, x1, x2, ytable, options ) result(value)
    integer, parameter :: RK = kind(0.0e0)
    ! Args
    real(rk), intent(in)                        :: x
    real(rk), dimension(:), intent(in)          :: ytable
    real(rk), intent(in)                        :: x1
    real(rk), intent(in)                        :: x2
    character(len=*), optional, intent(in)      :: options
    real(rk)                                    :: value    ! x are presorted
    value = UseLookUpTable_r4 ( x, ytable, x1, x2, options=options )
  end function FLookup_r4

  function FLookup_r8 ( x, x1, x2, ytable, options ) result(value)
    integer, parameter :: RK = kind(0.0d0)
    ! Args
    real(rk), intent(in)                        :: x
    real(rk), dimension(:), intent(in)          :: ytable
    real(rk), intent(in)                        :: x1
    real(rk), intent(in)                        :: x2
    character(len=*), optional, intent(in)      :: options
    real(rk)                                    :: value    ! x are presorted
    value = UseLookUpTable_r8 ( x, ytable, x1, x2, options=options )
  end function FLookup_r8

  function FXYLookup_r4 ( x, xtable, ytable, options ) result(value)
    integer, parameter :: RK = kind(0.0e0)
    ! Args
    real(rk), intent(in)                        :: x
    real(rk), dimension(:), intent(in)          :: xtable
    real(rk), dimension(:), intent(in)          :: ytable
    character(len=*), optional, intent(in)      :: options
    real(rk)                                    :: value    ! x are presorted
    value = UseLookUpTable_r4 ( x, ytable, xtable=xtable, &
      & options=Default(options, '-p', additional=.true.) )
  end function FXYLookup_r4

  function FXYLookup_r8 ( x, xtable, ytable, options ) result(value)
    integer, parameter :: RK = kind(0.0d0)
    ! Args
    real(rk), intent(in)                        :: x
    real(rk), dimension(:), intent(in)          :: xtable
    real(rk), dimension(:), intent(in)          :: ytable
    character(len=*), optional, intent(in)      :: options
    real(rk)                                    :: value    ! x are presorted
    value = UseLookUpTable_r8 ( x, ytable, xtable=xtable, &
      & options=Default(options, '-p', additional=.true.) )
  end function FXYLookup_r8
    
! -------------------------------------------------  FApproximate  -----

  ! This family of routines use a LookUpTable to approximate a 
  ! (costly-to-evaluate) function based on its values at a set of points

  ! Args: (* means optional)
  ! x        pt at which to evaluate
  ! LkUpTable      the Look Up Table type

  function FApproximate_r4 ( x, LkUpTable ) result(value)
    integer, parameter :: RK = kind(0.0e0)
    ! Args
    real(rk), intent(in)                 :: x
    type(LookUpTable_0(r4)), intent(in)       :: LkUpTable
    real(rk)                             :: value
    ! Internal variables
    real(rk) :: arg
    integer  :: itsSign
    character :: options
    real(rk), dimension(:), pointer      :: xArray => null()
    ! Executable
    call reposit( x, LkUpTable, arg, itsSign )
    ! call outputNamedValue( 'x', x )
    ! call outputNamedValue( 'arg', arg )
    ! call outputNamedValue( 'itsSign', itsSign )
    options = LkUpTable%method
    if ( options == 'q' ) options = '0'
    if ( arg < LkUpTable%x1 ) then
      value = LkUpTable%yLeft
    elseif ( arg == LkUpTable%x1 ) then
      value = LkUpTable%y(1)
    elseif ( arg == LkUpTable%x2 ) then
      value = LkUpTable%y(LkUpTable%N)
    elseif ( arg > LkUpTable%x2 ) then
      value = LkUpTable%yRight
    elseif ( options == 's' ) then
      call createXArray( xArray, LkUpTable )
      call InterpolateValues( xArray, LkUpTable%y, (/arg/), newYr4, method='S' )
      value = newYr4(1)
      call deallocate_test( xArray, 'xArray (r4)', ModuleName )
    else
      ! call outputNamedValue( 'options', options )
      value = UseLookUpTable ( arg, LkUpTable%y, LkUpTable%x1, LkUpTable%x2, options=options )
    end if
    value = itsSign*value
  end function FApproximate_r4

  function FApproximate_r8 ( x, LkUpTable ) result(value)
    integer, parameter :: RK = kind(0.0d0)
    ! Args
    real(rk), intent(in)                 :: x
    type(LookUpTable_0(r8)), intent(in)       :: LkUpTable
    real(rk)                             :: value
    ! Internal variables
    real(rk) :: arg
    integer  :: itsSign
    character :: options
    real(rk), dimension(:), pointer      :: xArray => null()
    ! Executable
    call reposit( x, LkUpTable, arg, itsSign )
    ! call outputNamedValue( 'x', x )
    ! call outputNamedValue( 'arg', arg )
    ! call outputNamedValue( 'itsSign', itsSign )
    options = LkUpTable%method
    if ( options == 'q' ) options = '0'
    if ( arg < LkUpTable%x1 ) then
      value = LkUpTable%yLeft
    elseif ( arg == LkUpTable%x1 ) then
      value = LkUpTable%y(1)
    elseif ( arg == LkUpTable%x2 ) then
      value = LkUpTable%y(LkUpTable%N)
    elseif ( arg > LkUpTable%x2 ) then
      value = LkUpTable%yRight
    elseif ( options == 's' ) then
      call createXArray( xArray, LkUpTable )
      call InterpolateValues( xArray, LkUpTable%y, (/arg/), newYr8, method='S' )
      value = newYr8(1)
      call deallocate_test( xArray, 'xArray (r8)', ModuleName )
    else
      ! call outputNamedValue( 'options', options )
      value = UseLookUpTable ( arg, LkUpTable%y, LkUpTable%x1, LkUpTable%x2, options=options )
    end if
    value = itsSign*value
  end function FApproximate_r8

! ----------------------------------------------  FInvApproximate  -----

  ! This family of routines use a LookUpTable to approximately invert 
  ! a function possibly restricting the search to a range [xS, xE]

  ! Args: (* means optional)
  ! y        y value to invert
  ! LkUpTable      the Look Up Table type
  ! [xS,xE]  range in which to search (otherwise [LkUpTable%x1, LkUpTable%x2])

  ! Will return x such that y[x] is approximately y

  function FInvApproximate_r4 ( y, LkUpTable, xS, xE ) result(x)
    integer, parameter :: RK = kind(0.0e0)
    ! Args
    real(rk), intent(in)                 :: y
    type(LookUpTable_0(r4)), intent(in)  :: LkUpTable
    real(rk), optional, intent(in)       :: xS
    real(rk), optional, intent(in)       :: xE
    real(rk)                             :: x
    ! Internal variables
    real(rk), dimension(:), pointer      :: xArray => null()
    ! Executable
    call createXArray( xArray, LkUpTable )
    x = UseLookUpTable( y, xArray, xtable=LkUpTable%y, yBottom=xS, yTop=xE )
    call deallocate_test( xArray, 'xArray (r4)', ModuleName )
  end function FInvApproximate_r4

  function FInvApproximate_r8 ( y, LkUpTable, xS, xE ) result(x)
    integer, parameter :: RK = kind(0.0d0)
    ! Args
    real(rk), intent(in)                 :: y
    type(LookUpTable_0(r8)), intent(in)  :: LkUpTable
    real(rk), optional, intent(in)       :: xS
    real(rk), optional, intent(in)       :: xE
    real(rk)                             :: x
    ! Internal variables
    real(rk), dimension(:), pointer      :: xArray => null()
    ! Executable
    call createXArray( xArray, LkUpTable )
    x = UseLookUpTable( y, xArray, xtable=LkUpTable%y, yBottom=xS, yTop=xE )
    call deallocate_test( xArray, 'xArray (r8)', ModuleName )
  end function FInvApproximate_r8

! ----------------------------------------------  FillLookUpTable  -----

  ! This family of routines fills a table with evaluations of a function
  ! at regularly-spaced points
  ! Subsequently, instead of evaluating the function you can address the
  ! array at the index of its closest element

  ! This only makes sense if you're going to evaluate the function many more
  ! times than takes to fill the table and you're willing 
  ! to accept whatever error may result from using the ClosestValue
  ! Of course that error could be large if you will evaluate the function
  ! outside the range [x1, x2]

  subroutine FillLookUpTable_r4 ( fun, table, x1, x2, N, xtable )
    integer, parameter :: RK = kind(0.0e0)
    include 'FillLookUpTable.f9h'
  end subroutine FillLookUpTable_r4 

  subroutine FillLookUpTable_r8 ( fun, table, x1, x2, N, xtable )
    integer, parameter :: RK = kind(0.0d0)
    include 'FillLookUpTable.f9h'
  end subroutine FillLookUpTable_r8

! --------------------------------------------------  FindInRange  -----
! This family of subroutines search not for a single index
! but for all at which the corresponding elements
! lie within a range of values, inclusive
! If none, return (/ 0, .., 0 /)
! Special interpretation: inclusive means
! if vrange(1) == vrange(2), any values of list also == vrange
! are within that range
! Unlike Hunt, list need not be monotonic
! options may include any of the following characters
!    character            meaning
!    ---------            -------
!       a                 ignore sign
!       r                 reverse sense (i.e. find outside range)
!       c                 modulo 360 degress
  subroutine FindInRange_int ( list, vrange, which, how_many, options )
    integer, parameter :: RK = kind(0.0e0)
    ! Dummy args
    integer, dimension(:) :: list
    integer, dimension(2) :: vrange
    include 'FindInRange.f9h'
  end subroutine FindInRange_int

  subroutine FindInRange_r4 ( list, vrange, which, how_many, options )
    integer, parameter :: RK = kind(0.0e0)
    ! Dummy args
    real(rk), dimension(:) :: list
    real(rk), dimension(2) :: vrange
    include 'FindInRange.f9h'
  end subroutine FindInRange_r4

  subroutine FindInRange_r8 ( list, vrange, which, how_many, options )
    integer, parameter :: RK = kind(0.0d0)
    ! Dummy args
    real(rk), dimension(:) :: list
    real(rk), dimension(2) :: vrange
    include 'FindInRange.f9h'
  end subroutine FindInRange_r8

  subroutine FindInRange_2d_int ( list, vrange, which, how_many, options )
    ! Dummy args
    integer, dimension(:,:) :: list
    integer, dimension(2) :: vrange
    include 'FindInRange_2d.f9h'
  end subroutine FindInRange_2d_int

  subroutine FindInRange_2d_r4 ( list, vrange, which, how_many, options )
    integer, parameter :: RK = kind(0.0e0)
    ! Dummy args
    real(rk), dimension(:,:) :: list
    real(rk), dimension(2) :: vrange
    include 'FindInRange_2d.f9h'
  end subroutine FindInRange_2d_r4

  subroutine FindInRange_2d_r8 ( list, vrange, which, how_many, options )
    integer, parameter :: RK = kind(0.0d0)
    ! Dummy args
    real(rk), dimension(:,:) :: list
    real(rk), dimension(2) :: vrange
    include 'FindInRange_2d.f9h'
  end subroutine FindInRange_2d_r8

! ------------------------------------------------  IFApproximate  -----

  ! This family of routines use a LookUpTable to approximate a 
  ! function's integral

  ! Args: (* means optional)
  ! LkUpTable      the Look Up Table type
  ! [xS,xE]  range over which to integrate (otherwise [LkUpTable%x1, LkUpTable%x2])

  function IFApproximate_r4 ( LkUpTable, xS, xE ) result(value)
    integer, parameter :: RK = kind(0.0e0)
    ! Args
    type(LookUpTable_0(r4)), intent(in)  :: LkUpTable
    real(rk), optional, intent(in)       :: xS
    real(rk), optional, intent(in)       :: xE
    real(rk)                             :: value
    ! Internal variables
    character :: options
    real(rk), dimension(:), pointer      :: xArray => null()
    real(rk) :: x1, x2
    ! Executable
    x1 = LkUpTable%x1
    x2 = LkUpTable%x2
    if ( present(xS) ) x1 = xS
    if ( present(xE) ) x2 = xE
    options = LkUpTable%method
    if ( options == 'q' ) options = 'S'
    if ( x2 < LkUpTable%x1 ) then
      value = 0.
    elseif ( x1 > LkUpTable%x2 ) then
      value = 0.
    elseif ( options == '/s/' ) then ! We'll just use our standard integration
      call createXArray( xArray, LkUpTable )
      call MLSMessage &
      & ( MLSMSG_Warning, ModuleName, "Unable to integrate with cubic spline")
      call deallocate_test( xArray, 'xArray (r4)', ModuleName )
    else
      value = UseLookUpTable ( x2, LkUpTable%y, LkUpTable%x1, LkUpTable%x2, options=options // 'S' ) &
           & - &
           &  UseLookUpTable ( x1, LkUpTable%y, LkUpTable%x1, LkUpTable%x2, options=options // 'S' )
    end if
  end function IFApproximate_r4

  function IFApproximate_r8 ( LkUpTable, xS, xE ) result(value)
    integer, parameter :: RK = kind(0.0d0)
    ! Args
    type(LookUpTable_0(r8)), intent(in)  :: LkUpTable
    real(rk), optional, intent(in)       :: xS
    real(rk), optional, intent(in)       :: xE
    real(rk)                             :: value
    ! Internal variables
    character :: options
    real(rk), dimension(:), pointer      :: xArray => null()
    real(rk) :: x1, x2
    ! Executable
    x1 = LkUpTable%x1
    x2 = LkUpTable%x2
    if ( present(xS) ) x1 = xS
    if ( present(xE) ) x2 = xE
    options = LkUpTable%method
    if ( options == 'q' ) options = 'S'
    if ( x2 < LkUpTable%x1 ) then
      value = 0.
    elseif ( x1 > LkUpTable%x2 ) then
      value = 0.
    elseif ( options == '/s/' ) then ! We'll just use our standard integration
      call createXArray( xArray, LkUpTable )
      call MLSMessage &
      & ( MLSMSG_Warning, ModuleName, "Unable to integrate with cubic spline")
      call deallocate_test( xArray, 'xArray (r8)', ModuleName )
    else
      value = UseLookUpTable ( x2, LkUpTable%y, LkUpTable%x1, LkUpTable%x2, options=options // 'S' ) &
           & - &
           &  UseLookUpTable ( x1, LkUpTable%y, LkUpTable%x1, LkUpTable%x2, options=options // 'S' )
    end if
  end function IFApproximate_r8

! ------------------------------------------  InterpolateArray_r4  -----

  ! This next subroutine is a workhorse interpolation routine, loosely based on
  ! my (Nathaniel) IDL routine of the same name.

  ! Method is one of 'L'inear, 'C'spline, or 'S'pline
  !                                (Numerical Recipes, more later no doubt)
  ! Extrapolate is one of 'A'llow, 'C'onstant, 'B'ad or 'P'eriodic (Spline only)

  ! The 'C' spline was simply moved here from the fwdmdl directory
  ! it appears similar to the 'S'pline except constraining the newY values
  ! that exceed either
  ! (a) 125% of the value that would result from linear interpolation; or
  ! (b) the hard bounds set by the optional parameters [yMin, ymax]

  ! Notes:
  !   oldX must be monotonically increasing or decreasing
  !   newX can be in any order
  !   one can't ask for spline interpolation with missing regions.
  !   missingRegions will probably slow the code down, as will extrapolate=B

  subroutine InterpolateArray_r4 ( oldX, oldY, newX, newY, method, extrapolate, &
    & badValue, missingRegions, dyByDx, dNewByDOld, skipNewY, IntYdX, Second )
    integer, parameter :: RK = kind(0.0e0)

    ! Dummy arguments
    real(rk), dimension(:), intent(IN) :: oldX
    real(rk), dimension(:,:), intent(IN) :: oldY   ! See Second argument below
    real(rk), dimension(:), intent(IN) :: newX
    real(rk), dimension(:,:), intent(OUT) :: newY  ! See Second argument below

    character (len=*), intent(in) :: method ! See comments above
    character (len=*), optional, intent(in) :: extrapolate ! See comments above
    real(rk), optional, intent(in) :: badvalue
    logical, optional, intent(in) :: missingRegions ! Allow missing regions
    real(rk), dimension(:,:), optional, intent(out) :: dyByDx
    type (matrixElement_T), intent(out), optional :: dNewByDOld ! Derivatives
    logical, optional, intent(in) :: SKIPNEWY ! Don't compute newY
    real(rk), dimension(:,:), optional, intent(out) :: IntYdX ! Antiderivative
                                              ! of Y at X
    logical, optional, intent(in) :: Second   ! Interpolate on the second
                                              ! dimension of OldY, NewY, default
                                              ! false.

    type(coefficients(r4)) :: Coeffs

    include "InterpolateArray.f9h"

  end subroutine InterpolateArray_r4

! ------------------------------------------  InterpolateArray_r8  -----

  subroutine InterpolateArray_r8 ( oldX, oldY, newX, newY, method, extrapolate, &
    & badValue, missingRegions, dyByDx, dNewByDOld, skipNewY, IntYdX, Second )
    integer, parameter :: RK = kind(0.0d0)

    ! Dummy arguments
    real(rk), dimension(:), intent(IN) :: oldX
    real(rk), dimension(:,:), intent(IN) :: oldY   ! See Second argument below
    real(rk), dimension(:), intent(IN) :: newX
    real(rk), dimension(:,:), intent(OUT) :: newY  ! See Second argument below

    character (len=*), intent(in) :: method ! See comments above
    character (len=*), optional, intent(in) :: extrapolate ! See comments above
    real(rk), optional, intent(in) :: badvalue
    logical, optional, intent(in) :: missingRegions ! Allow missing regions
    real(rk), dimension(:,:), optional, intent(out) :: dyByDx
    type (matrixElement_T), intent(OUT), optional :: dNewByDOld ! Derivatives
    logical, optional, intent(in) :: SKIPNEWY ! Don't compute newY
    real(rk), dimension(:,:), optional, intent(out) :: IntYdX ! Antiderivative
                                              ! of Y at X
    logical, optional, intent(in) :: Second   ! Interpolate on the second
                                              ! dimension of OldY, NewY, default
                                              ! false.

    type(coefficients(r8)) :: Coeffs

    include "InterpolateArray.f9h"

  end subroutine InterpolateArray_r8

! -------------------------------------  InterpolateArraySetup_r4  -----

  subroutine InterpolateArraySetup_r4 ( OldX, NewX, Method, Coeffs, &
    & Extrapolate, Width, DyByDx, dNewByDOld, IntYdX, fail )

    integer, parameter :: RK = kind(0.0e0)

    real(rk), intent(in) :: OldX(:), NewX(:)
    character(len=*), intent(in) :: Method
    type(coefficients(r4)), intent(out) :: Coeffs
    character(len=*), intent(in), optional :: Extrapolate ! See comments above
    integer, intent(in), optional :: Width ! Second dimension for OldY when
                                           ! interpolations get done
    logical, optional, intent(in) :: DyByDx ! just a signal to
                                           ! compute more coeffs for splines
    type (matrixElement_T), intent(inout), optional :: dNewByDOld ! Derivatives
                                           ! inout so createBlock can clean up
                                           ! after an old one
    logical, optional, intent(in) :: IntYdX ! just a signal to
                                           ! compute more coeffs for splines
    logical, optional, intent(out) :: Fail    ! True for failure

    include "InterpolateArraySetup.f9h"

  end subroutine InterpolateArraySetup_r4

! -------------------------------------  InterpolateArraySetup_r8  -----

  subroutine InterpolateArraySetup_r8 ( OldX, NewX, Method, Coeffs, &
    & Extrapolate, Width, DyByDx, dNewByDOld, IntYdX, fail )

    integer, parameter :: RK = kind(0.0d0)

    real(rk), intent(in) :: OldX(:), NewX(:)
    character(len=*), intent(in) :: Method
    type(coefficients(r8)), intent(out) :: Coeffs
    character(len=*), intent(in), optional :: Extrapolate ! See comments above
    integer, intent(in), optional :: Width ! Second dimension for OldY when
                                           ! interpolations get done
    logical, optional, intent(in) :: DyByDx ! just a signal to
                                           ! compute more coeffs for splines
    type (matrixElement_T), intent(inout), optional :: dNewByDOld ! Derivatives
                                           ! inout so createBlock can clean up
                                           ! after an old one
    logical, optional, intent(in) :: IntYdX ! just a signal to
                                           ! compute more coeffs for splines
    logical, optional, intent(out) :: Fail    ! True for failure

    include "InterpolateArraySetup.f9h"

  end subroutine InterpolateArraySetup_r8

! ----------------------------------  InterpolateArrayTeardown_r4  -----

  subroutine InterpolateArrayTeardown_r4 ( Coeffs )

    type(coefficients(r4)), intent(inout) :: Coeffs

    include "InterpolateArrayTeardown.f9h"

  end subroutine InterpolateArrayTeardown_r4

! ----------------------------------  InterpolateArrayTeardown_r8  -----

  subroutine InterpolateArrayTeardown_r8 ( Coeffs )

    type(coefficients(r8)), intent(inout) :: Coeffs

    include "InterpolateArrayTeardown.f9h"

  end subroutine InterpolateArrayTeardown_r8

! -------------------------------------  InterpolateExtrapolate_d  -----

  subroutine InterpolateExtrapolate_d ( OldX, OldY, NewX, NewY, Second )

    ! Interpolate ( OldX(:), OldY(:,:) ) to ( NewX(:), NewY(:,:) ), where
    ! the dimension upon which to interpolate is merge(2,1,second).
    ! Extrapolate outside the range of OldX using the average slope, not
    ! the slope nearest the end where NewX is outside the range of OldX.

    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: OldX(:), OldY(:,:), NewX(:)
    real(rk), intent(out) :: NewY(:,:)
    logical, intent(in) :: Second
    type(coefficients(r8)) :: Coeffs

    include 'InterpolateExtrapolate.f9h'

  end subroutine InterpolateExtrapolate_d

! -----------------------------------  InterpolateExtrapolate_d_1  -----

  subroutine InterpolateExtrapolate_d_1 ( OldX, OldY, NewX, NewY, Second )

    ! Interpolate ( OldX(:), OldY(:) ) to ( NewX(:), NewY(:) ).
    ! Extrapolate outside the range of OldX using the average slope, not
    ! the slope nearest the end where NewX is outside the range of OldX.

    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: OldX(:), OldY(:), NewX(:)
    real(rk), intent(out) :: NewY(:)
    logical, intent(in) :: Second
    type(coefficients(r8)) :: Coeffs

    include 'InterpolateExtrapolate_1.f9h'

  end subroutine InterpolateExtrapolate_d_1

! -------------------------------------  InterpolateExtrapolate_s  -----

  subroutine InterpolateExtrapolate_s ( OldX, OldY, NewX, NewY, Second )

    ! Interpolate ( OldX(:), OldY(:,:) ) to ( NewX(:), NewY(:,:) ), where
    ! the dimension upon which to interpolate is merge(2,1,second).
    ! Extrapolate outside the range of OldX using the average slope, not
    ! the slope nearest the end where NewX is outside the range of OldX.

    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: OldX(:), OldY(:,:), NewX(:)
    real(rk), intent(out) :: NewY(:,:)
    logical, intent(in) :: Second
    type(coefficients(r4)) :: Coeffs

    include 'InterpolateExtrapolate.f9h'

  end subroutine InterpolateExtrapolate_s

! -----------------------------------  InterpolateExtrapolate_s_1  -----

  subroutine InterpolateExtrapolate_s_1 ( OldX, OldY, NewX, NewY, Second )

    ! Interpolate ( OldX(:), OldY(:) ) to ( NewX(:), NewY(:) ).
    ! Extrapolate outside the range of OldX using the average slope, not
    ! the slope nearest the end where NewX is outside the range of OldX.

    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: OldX(:), OldY(:), NewX(:)
    real(rk), intent(out) :: NewY(:)
    logical, intent(in) :: Second
    type(coefficients(r4)) :: Coeffs

    include 'InterpolateExtrapolate_1.f9h'

  end subroutine InterpolateExtrapolate_s_1

! -----------------------------------------  InterpolateScalar_r4  -----

! This subroutine is a scalar wrapper for the array one.

  subroutine InterpolateScalar_r4 ( oldX, oldY, newX, newY, method, extrapolate, &
    & badValue, missingRegions, dyByDx, RangeOfPeriod, skipNewY, IntYdX, YMIN, YMAX )
    integer, parameter :: RK = kind(0.0e0)

    ! Dummy arguments
    real(rk), dimension(:), intent(in) :: oldX
    real(rk), dimension(:), intent(in) :: oldY
    real(rk), dimension(:), intent(in) :: newX
    real(rk), dimension(:), intent(out) :: newY

    character (len=*), intent(in) :: method ! See comments above
    character (len=*), optional, intent(in) :: extrapolate ! See comments above
    real(rk), optional, intent(in) :: badValue
    real(rk), dimension(:), optional, intent(out) :: dyByDx
    real(rk), dimension(2), optional, intent(in) :: rangeofperiod   ! for periodic data
    logical, optional, intent(in) :: missingRegions ! Allow missing regions
    logical, optional, intent(in) :: SKIPNEWY ! Don't compute newY
    real(rk), dimension(:), optional, intent(out) :: IntYdX ! Antiderivative
                                              ! of Y at X
    real(rk), optional, intent(in) :: YMIN, YMAX

    include "InterpolateScalar.f9h"
  end subroutine InterpolateScalar_r4

! -----------------------------------------  InterpolateScalar_r8  -----

  subroutine InterpolateScalar_r8 ( oldX, oldY, newX, newY, method, extrapolate, &
    & badValue, missingRegions, dyByDx, RangeOfPeriod, skipNewY, IntYdX, YMIN, YMAX )
    integer, parameter :: RK = kind(0.0d0)

    ! Dummy arguments
    real(rk), dimension(:), intent(in) :: oldX
    real(rk), dimension(:), intent(in) :: oldY
    real(rk), dimension(:), intent(in) :: newX
    real(rk), dimension(:), intent(out) :: newY

    character (len=*), intent(in) :: method ! See comments above
    character (len=*), optional, intent(in) :: extrapolate ! See comments above
    real(rk), optional, intent(in) :: badValue
    logical, optional, intent(in) :: missingRegions ! Allow missing regions
    real(rk), dimension(:), optional, intent(out) :: dyByDx
    real(rk), dimension(2), optional, intent(in) :: rangeofperiod   ! for periodic data
    logical, optional, intent(in) :: SKIPNEWY ! Don't compute newY
    real(rk), dimension(:), optional, intent(out) :: IntYdX ! Antiderivative
                                              ! of Y at X
    real(rk), optional, intent(in) :: YMIN, YMAX

    include "InterpolateScalar.f9h"
  end subroutine InterpolateScalar_r8

! -------------------------------  InterpolateScalarUsingSetup_r4  -----

  subroutine InterpolateScalarUsingSetup_r4 ( coeffs, oldX, oldY, newX, newY, &
    & method, extrapolate, badValue, missingRegions, dyByDx, skipNewY, IntYdX )
    integer, parameter :: RK = kind(0.0e0)

    ! Dummy arguments
    type(coefficients(r4)), intent(in) :: Coeffs
    real(rk), dimension(:), intent(in) :: oldX
    real(rk), dimension(:), intent(in) :: oldY
    real(rk), dimension(:), intent(in) :: newX
    real(rk), dimension(:), intent(out) :: newY

    character (len=*), intent(in) :: method ! See comments above
    character (len=*), optional, intent(in) :: extrapolate ! See comments above
    real(rk), optional, intent(in) :: badvalue
    logical, optional, intent(in) :: missingRegions ! Allow missing regions
    real(rk), dimension(:), optional, intent(out) :: dyByDx
    logical, optional, intent(in) :: SKIPNEWY ! Don't compute newY
    real(rk), dimension(:), optional, intent(out) :: IntYdX ! Antiderivative
                                              ! of Y at X

    include "InterpolateScalarUsingSetup.f9h"

  end subroutine InterpolateScalarUsingSetup_r4

! -------------------------------  InterpolateScalarUsingSetup_r8  -----

  subroutine InterpolateScalarUsingSetup_r8 ( coeffs, oldX, oldY, newX, newY, &
    & method, extrapolate, badValue, missingRegions, dyByDx, skipNewY, IntYdX )
    integer, parameter :: RK = kind(0.0d0)

    ! Dummy arguments
    type(coefficients(r8)), intent(in) :: Coeffs
    real(rk), dimension(:), intent(in) :: oldX
    real(rk), dimension(:), intent(in) :: oldY
    real(rk), dimension(:), intent(in) :: newX
    real(rk), dimension(:), intent(out) :: newY

    character (len=*), intent(in) :: method ! See comments above
    character (len=*), optional, intent(in) :: extrapolate ! See comments above
    real(rk), optional, intent(in) :: badvalue
    logical, optional, intent(in) :: missingRegions ! Allow missing regions
    real(rk), dimension(:), optional, intent(out) :: dyByDx
    logical, optional, intent(in) :: SKIPNEWY ! Don't compute newY
    real(rk), dimension(:), optional, intent(out) :: IntYdX ! Antiderivative
                                              ! of Y at X

    include "InterpolateScalarUsingSetup.f9h"

  end subroutine InterpolateScalarUsingSetup_r8

! -------------------------------------  InterpolateUsingSetup_r4  -----

  subroutine InterpolateUsingSetup_r4 ( coeffs, oldX, oldY, newX, newY, &
    & method, extrapolate, badValue, missingRegions, dyByDx, skipNewY, IntYdX, &
    & Second )
    integer, parameter :: RK = kind(0.0e0)

    ! Dummy arguments
    type(coefficients(r4)), intent(in) :: Coeffs
    real(rk), dimension(:), intent(in) :: oldX
    real(rk), dimension(:,:), intent(in) :: oldY   ! See Second argument below
    real(rk), dimension(:), intent(in) :: newX
    real(rk), dimension(:,:), intent(out) :: newY  ! See Second argument below

    character (len=*), intent(in) :: method ! See comments above
    character (len=*), optional, intent(in) :: extrapolate ! See comments above
    real(rk), optional, intent(in) :: badvalue
    logical, optional, intent(in) :: missingRegions ! Allow missing regions
    real(rk), dimension(:,:), optional, intent(out) :: dyByDx
    logical, optional, intent(in) :: SKIPNEWY ! Don't compute newY
    real(rk), dimension(:,:), optional, intent(out) :: IntYdX ! Antiderivative
                                              ! of Y at X
    logical, optional, intent(in) :: Second   ! Interpolate on the second
                                              ! dimension of OldY, NewY, default
                                              ! false.

    include "InterpolateUsingSetup.f9h"

  end subroutine InterpolateUsingSetup_r4

! -------------------------------------  InterpolateUsingSetup_r8  -----

  subroutine InterpolateUsingSetup_r8 ( coeffs, oldX, oldY, newX, newY, &
    & method, extrapolate, badValue, missingRegions, dyByDx, skipNewY, IntYdX, &
    & Second )
    integer, parameter :: RK = kind(0.0d0)

    ! Dummy arguments
    type(coefficients(r8)), intent(in) :: Coeffs
    real(rk), dimension(:), intent(in) :: oldX
    real(rk), dimension(:,:), intent(in) :: oldY   ! See Second argument below
    real(rk), dimension(:), intent(in) :: newX
    real(rk), dimension(:,:), intent(out) :: newY  ! See Second argument below

    character (len=*), intent(in) :: method ! See comments above
    character (len=*), optional, intent(in) :: extrapolate ! See comments above
    real(rk), optional, intent(in) :: badvalue
    logical, optional, intent(in) :: missingRegions ! Allow missing regions
    real(rk), dimension(:,:), optional, intent(out) :: dyByDx
    logical, optional, intent(in) :: SKIPNEWY ! Don't compute newY
    real(rk), dimension(:,:), optional, intent(out) :: IntYdX ! Antiderivative
                                              ! of Y at X
    logical, optional, intent(in) :: Second   ! Interpolate on the second
                                              ! dimension of OldY, NewY, default
                                              ! false.

    include "InterpolateUsingSetup.f9h"

  end subroutine InterpolateUsingSetup_r8

  ! -----------------------------------  Interp_Bilinear_2d_1d_r4  -----
  subroutine Interp_Bilinear_2d_1d_r4 ( XOld, XNew, YOld, YNew, Zold, Znew, &
    & Update )

    ! Given ZOld on coordinates (XOld x YOld), interpolate to (XNew,YNew)
    ! to give ZNew.  ZOld must have shape (size(xOld),size(yOld)), while
    ! XNew, YNew and ZNew must have the same shapes.

    integer, parameter :: RK = kind(1.0e0)
    real(rk), intent(in) :: XOld(:)
    real(rk), intent(in) :: XNew(:)
    real(rk), intent(in) :: YOld(:)
    real(rk), intent(in) :: YNew(:)
    real(rk), intent(in) :: ZOld(:,:)
    real(rk), intent(inout) :: ZNew(:)
    logical, intent(in), optional :: Update ! Add interpolate to Znew
    include 'Interp_Bilinear_2d_1d.f9h'

  end subroutine Interp_Bilinear_2d_1d_r4

  ! -----------------------------------  Interp_Bilinear_2d_1d_r8  -----
  subroutine Interp_Bilinear_2d_1d_r8 ( XOld, XNew, YOld, YNew, Zold, Znew, &
    & Update )

    ! Given ZOld on coordinates (XOld x YOld), interpolate to (XNew,YNew)
    ! to give ZNew.  ZOld must have shape (size(xOld),size(yOld)), while
    ! XNew, YNew and ZNew must have the same shapes.

    integer, parameter :: RK = kind(1.0d0)
    real(rk), intent(in) :: XOld(:)
    real(rk), intent(in) :: XNew(:)
    real(rk), intent(in) :: YOld(:)
    real(rk), intent(in) :: YNew(:)
    real(rk), intent(in) :: ZOld(:,:)
    real(rk), intent(inout) :: ZNew(:)
    logical, intent(in), optional :: Update ! Add interpolate to Znew
    include 'Interp_Bilinear_2d_1d.f9h'

  end subroutine Interp_Bilinear_2d_1d_r8

  ! ------------------------  Interpolate_Regular_To_Irregular_r4  -----
  subroutine Interpolate_Regular_To_Irregular_r4 ( XOld, YOld, ZOld, &
    & XNew, YNew, ZNew )
    ! Given XOld, YOld, ZOld, XYNew, interpolate to ZNew.
    ! The shape of ZOld must be (size(xOld),size(yOld)).
    ! The shape of ZNew must be the same as Xnew and YNew.
    ! The only method supported is linear with constant extrapolation.
    integer, parameter :: RK = kind(1.0e0)
    real(rk), intent(in) :: XOld(:)
    real(rk), intent(in) :: YOld(:)
    real(rk), intent(in) :: ZOld(:,:)
    real(rk), intent(in) :: XNew(:,:), YNew(:,:)
    real(rk), intent(out) :: ZNew(:,:)
    include 'Interpolate_Regular_To_Irregular.f9h'
   end subroutine Interpolate_Regular_To_Irregular_r4

  ! ------------------------  Interpolate_Regular_To_Irregular_r8  -----
  subroutine Interpolate_Regular_To_Irregular_r8 ( XOld, YOld, ZOld, &
    & XNew, YNew, ZNew )
    ! Given XOld, YOld, ZOld, XYNew, interpolate to ZNew.
    ! The shape of ZOld must be (size(xOld),size(yOld)).
    ! The shape of ZNew must be the same as Xnew and YNew.
    ! The only method supported is linear with constant extrapolation.
    integer, parameter :: RK = kind(1.0d0)
    real(rk), intent(in) :: XOld(:)
    real(rk), intent(in) :: YOld(:)
    real(rk), intent(in) :: ZOld(:,:)
    real(rk), intent(in) :: XNew(:,:), YNew(:,:)
    real(rk), intent(out) :: ZNew(:,:)
    include 'Interpolate_Regular_To_Irregular.f9h'
   end subroutine Interpolate_Regular_To_Irregular_r8

  ! --------------------------------  Interpolate_2d_Composite_r4  -----
  subroutine Interpolate_2d_Composite_r4 ( XOld, YOld, ZOld, XNew, YNew, ZNew, &
    & XMethod, YMethod, XExtrapolate, YExtrapolate )
    ! Given XOld, YOld, ZOld, XNew and YNew, interpolate to ZNew.
    ! The shape of ZOld must be (size(xOld),size(yOld)).
    ! The shape of ZNew must be (size(xNew),size(yNew)).
    ! Interpolation is first done in the X direction using XMethod to produce
    ! a temporary variable with shape (size(xNew),size(yOld)).  Then
    ! interpolation is done in the Y direction, from that temporary variable,
    ! using YMethod, to produce ZNew.
    integer, parameter :: RK = kind(1.0e0)
    real(rk), intent(in) :: XOld(:)
    real(rk), intent(in) :: YOld(:)
    real(rk), intent(in) :: ZOld(:,:)
    real(rk), intent(in) :: XNew(:)
    real(rk), intent(in) :: YNew(:)
    real(rk), intent(out) :: ZNew(:,:)
    character (len=*), intent(in) :: XMethod, YMethod ! See comments above
    character (len=*), optional, intent(in) :: XExtrapolate, YExtrapolate ! See comments above
    type(coefficients(r4)) :: XCoeffs, YCoeffs
    include 'Interpolate_2d_Composite.f9h'
  end subroutine Interpolate_2d_Composite_r4

  ! --------------------------------  Interpolate_2d_Composite_r8  -----
  subroutine Interpolate_2d_Composite_r8 ( XOld, YOld, ZOld, XNew, YNew, ZNew, &
    & XMethod, YMethod, XExtrapolate, YExtrapolate )
    ! Given XOld, YOld, ZOld, XNew and YNew, interpolate to ZNew.
    ! The shape of ZOld must be (size(xOld),size(yOld)).
    ! The shape of ZNew must be (size(xNew),size(yNew)).
    ! Interpolation is first done in the X direction using XMethod to produce
    ! a temporary variable with shape (size(xNew),size(yOld)).  Then
    ! interpolation is done in the Y direction, from that temporary variable,
    ! using YMethod, to produce ZNew.
    integer, parameter :: RK = kind(1.0d0)
    real(rk), intent(in) :: XOld(:)
    real(rk), intent(in) :: YOld(:)
    real(rk), intent(in) :: ZOld(:,:)
    real(rk), intent(in) :: XNew(:)
    real(rk), intent(in) :: YNew(:)
    real(rk), intent(out) :: ZNew(:,:)
    character (len=*), intent(in) :: XMethod, YMethod ! See comments above
    character (len=*), optional, intent(in) :: XExtrapolate, YExtrapolate ! See comments above
    type(coefficients(r8)) :: XCoeffs, YCoeffs
    include 'Interpolate_2d_Composite.f9h'
  end subroutine Interpolate_2d_Composite_r8

! --------------------------------------------  LinearInterpolate  -----

  ! This family of functions return a value interpolated across
  ! multiple dimensions

  ! Args: (* means optional)
  ! values   the 2^n values at the hypercube's vertices
  ! coords   2*n vertices coordinates

  function LinearInterpolate_1d_r4 ( values, coords, vertices ) result(value)
    integer, parameter :: RK = kind(0.0e0)
    include "LinearInterpolate_1d.f9h"
  end function LinearInterpolate_1d_r4

  function LinearInterpolate_1d_r8 ( values, coords, vertices ) result(value)
    integer, parameter :: RK = kind(0.0d0)
    include "LinearInterpolate_1d.f9h"
  end function LinearInterpolate_1d_r8

  function LinearInterpolate_2d_r4 ( values, coords, vertices ) result(value)
    integer, parameter :: RK = kind(0.0e0)
    include "LinearInterpolate_2d.f9h"
  end function LinearInterpolate_2d_r4

  function LinearInterpolate_2d_r8 ( values, coords, vertices ) result(value)
    integer, parameter :: RK = kind(0.0d0)
    include "LinearInterpolate_2d.f9h"
  end function LinearInterpolate_2d_r8

  function LinearInterpolate_3d_r4 ( values, coords, vertices ) result(value)
    integer, parameter :: RK = kind(0.0e0)
    include "LinearInterpolate_3d.f9h"
  end function LinearInterpolate_3d_r4

  function LinearInterpolate_3d_r8 ( values, coords, vertices ) result(value)
    integer, parameter :: RK = kind(0.0d0)
    include "LinearInterpolate_3d.f9h"
  end function LinearInterpolate_3d_r8

  function LinearInterpolate_4d_r4 ( values, coords, vertices ) result(value)
    integer, parameter :: RK = kind(0.0e0)
    include "LinearInterpolate_4d.f9h"
  end function LinearInterpolate_4d_r4

  function LinearInterpolate_4d_r8 ( values, coords, vertices ) result(value)
    integer, parameter :: RK = kind(0.0d0)
    include "LinearInterpolate_4d.f9h"
  end function LinearInterpolate_4d_r8

! ----------------------------------------------  ReadLookUpTable  -----

  ! This family of routines Reads a table from a file where it is stored
  ! Something like this (although without the "!")
  ! Celsius2Fahrenheit.txt
  ! file of Celsius Fahrenheit conversions
  ! C        F
  ! -40.   -40.
  !  0.     32.
  ! 1.e2   212
  !
  ! Returning
  ! N = 3
  ! x = (-40., 0., 1.e2)
  ! y = (-40., 32., 212.)
  ! Note:
  ! We ignore lines with non-numerical values
  ! You may use the resulting x and y arrays as args to UseLookUpTable
  ! where
  !       x -> xtable
  !       y -> table
  subroutine ReadLookUpTable_r4 ( filename, N, x, y )
    integer, parameter :: RK = kind(0.0e0)
    include 'ReadLookUpTable.f9h'
  end subroutine ReadLookUpTable_r4

  subroutine ReadLookUpTable_r8 ( filename, N, x, y )
    integer, parameter :: RK = kind(0.0d0)
    include 'ReadLookUpTable.f9h'
  end subroutine ReadLookUpTable_r8

! --------------------------------------------------------  SetUp  -----

  ! This family of routines sets up a LookUpTable of the appropriate type
  ! We ought to do this for both _0 and _1 levels
  subroutine setUpLookUpTable_0_r4 ( LkUpTable, N, x1, x2, &
    & y, BC, method, yLeft, yRight, fun )
    integer, parameter :: RK = kind(0.0e0)
    ! Args
    class(LookUpTable_0(r4)) :: LkUpTable ! Intent(out) would clobber defaults
    integer, intent(in)  :: N
    real(rk), intent(in) :: x1
    real(rk), intent(in) :: x2
    real(rk), dimension(:), optional, intent(in) :: y
    character(len=*), optional, intent(in) :: BC
    character(len=*), optional, intent(in) :: method
    real(rk), optional, intent(in) :: yLeft
    real(rk), optional, intent(in) :: yRight
    real(rk), optional, external :: fun
    ! Internal args
    integer :: i
    ! Executable
    LkUpTable%N       = N
    LkUpTable%x1      = x1
    LkUpTable%x2      = x2
    call allocate_test ( LkUpTable%y, N, "LkUpTable%y", ModuleName )
    if ( present(y        ) ) LkUpTable%y           = y
    if ( present(BC       ) ) LkUpTable%BC          = BC
    if ( present(method   ) ) LkUpTable%method      = method
    if ( present(yLeft    ) ) LkUpTable%yLeft       = yLeft
    if ( present(yRight   ) ) LkUpTable%yRight      = yRight
    if ( .not. present(fun) ) return
    do i=1, N
      LkUpTable%y(i) = fun( x1 + (i-1)*(x2-x1)/(N-1) )
    enddo
  end subroutine setUpLookUpTable_0_r4

  subroutine setUpLookUpTable_0_r8 ( LkUpTable, N, x1, x2, &
    & y, BC, method, yLeft, yRight, fun )
    integer, parameter :: RK = kind(0.0d0)
    ! Args
    class(LookUpTable_0(r8)) :: LkUpTable ! Intent(out) would clobber defaults
    integer, intent(in)  :: N
    real(rk), intent(in) :: x1
    real(rk), intent(in) :: x2
    real(rk), dimension(:), optional, intent(in) :: y
    character(len=*), optional, intent(in) :: BC
    character(len=*), optional, intent(in) :: method
    real(rk), optional, intent(in) :: yLeft
    real(rk), optional, intent(in) :: yRight
    real(rk), optional, external :: fun
    ! Internal args
    integer :: i
    ! Executable
    LkUpTable%N       = N
    LkUpTable%x1      = x1
    LkUpTable%x2      = x2
    call allocate_test ( LkUpTable%y, N, "LkUpTable%y", ModuleName )
    if ( present(y        ) ) LkUpTable%y           = y
    if ( present(BC       ) ) LkUpTable%BC          = BC
    if ( present(method   ) ) LkUpTable%method      = method
    if ( present(yLeft    ) ) LkUpTable%yLeft       = yLeft
    if ( present(yRight   ) ) LkUpTable%yRight      = yRight
    if ( .not. present(fun) ) return
    do i=1, N
      LkUpTable%y(i) = fun( x1 + (i-1)*(x2-x1)/(N-1) )
    enddo
  end subroutine setUpLookUpTable_0_r8

! ------------------------------------------------  SimSubroutine  -----
  ! This family provides subroutine apis to integration by Simpson's rule
  subroutine Simps_r4 ( F, DX, N, R )
!  Simpson's Integration of discrete equal spacing
    integer, parameter :: RK = kind(0.0e0)
    include 'simpson.f9h'
  end subroutine Simps_r4

  subroutine Simps_r8 ( F, DX, N, R )
!  Simpson's Integration of discrete equal spacing
    integer, parameter :: RK = kind(0.0d0)
    include 'simpson.f9h'
  end subroutine Simps_r8

! -----------------------------------------------  SolveQuadratic  -----
   ! Solve

   ! a x^2 + b x + c = 0

   ! for x = r1 + i imPart, and x = r2 - i imPart
   ! where i^2 = -1

   ! Of course because a, b, and c are all purely real
   ! if imPart /= 0, r1 = r2

   ! We bother with this to avoid truncation that would result
   ! from taking a difference between like-signed quantities
   ! b and + or - sqrt(disc)

   ! Special cases (which you may prefer to intercept yourself)
   ! If a = 0, we return the same root in both r1 and r2
   ! If a = b = 0, we divide c by zero

  subroutine SolveQuadratic_r4 ( A, B, C, R1, R2, ImPart )
   integer, parameter :: RK = kind(0.0e0)
   real(rk), intent(in)  :: A, B, C
   real(rk), intent(out) :: R1, R2, ImPart
   include 'SolveQuadratic.f9h'
 end subroutine SolveQuadratic_r4

  subroutine SolveQuadratic_r8 ( A, B, C, R1, R2, ImPart )
   integer, parameter :: RK = kind(0.0d0)
   real(rk), intent(in)  :: A, B, C
   real(rk), intent(out) :: R1, R2, ImPart
   include 'SolveQuadratic.f9h'
 end subroutine SolveQuadratic_r8

! -----------------------------------------------  UseLookUpTable  -----

  ! This family of routines use a LookUpTable to approximate a costly-to-evaluate
  ! function based on its values at a set of points
  ! They are workhorse functions used by the more specialized functions
  ! that approximate derivatives, integrals, etc.

  ! Args:       (* means optional)
  ! x              pt at which to evaluate
  ! table          table of function values
  ! x1,x2          * range of equally-spaced argument values represented in table
  ! xtable         * table of argument values
  ! xS,xE          * range of x-values in which to search
  ! yBottom, yTop  * range of y-values in which to search
  !                 (ignored if xtable is sorted)

  ! options  * none, one, or more of the following:
  ! (default)   choose pt in xtable closest to x
  !    p        xtable is sorted (necessary for the other options)
  !    l        always choose lower of two closest x's in xtable
  !    u        always choose upper of two closest x's in xtable
  !    i        interpolate among two closest, but never extrapolate
  !    0        use quadratic interpolation
  !    1        return 1st derivative instead
  !    2        return 2nd derivative (at nearest pt)
  !    S        return definite integral from x1 to x
  !    C        return definite integral from x to x2

  ! notes:
  ! (1) You should supply either xtable or the range (x1,x2)
  !      (but not both)
  ! (2) If you supply xtable, and that xtable is sorted so that
  !     xtable[i] < xtable[i+1] include 'p' among the options
  ! (3) If you supply (x1,x2) the resulting x values will be automatically
  !     sorted, so you may omit 'p' from among the options
  ! (4) An unsorted xtable is incompatible with any of the options
  !      i.e., you may not integrate, differentiate, interpolate, etc.
  ! (5) The range yBottom, yTop is incompatible with a sorted xtable

  function UseLookUpTable_r4 ( x, table, x1, x2, xtable, &
    & missingValue, options, xS, xE, yBottom, yTop ) result(value)
    integer, parameter :: RK = kind(0.0e0)
    include 'UseLookUpTable.f9h'
  end function UseLookUpTable_r4 

  function UseLookUpTable_r8 ( x, table, x1, x2, xtable, &
    & missingValue, options, xS, xE, yBottom, yTop ) result(value)
    integer, parameter :: RK = kind(0.0d0)
    include 'UseLookUpTable.f9h'
  end function UseLookUpTable_r8 

! ----------------------------------------------  WriteLookUpTable  -----

  ! This family of routines Writes a table to a file where it is stored
  ! Something like this (although without the "!")
  ! Celsius2Fahrenheit.txt
  ! file of Celsius Fahrenheit conversions
  ! C        F
  ! -40.   -40.
  !  0.     32.
  ! 1.e2   212
  !
  ! Given the input
  ! N = 3
  ! x = (-40., 0., 1.e2)
  ! y = (-40., 32., 212.)
  subroutine WriteLookUpTable_r4 ( filename, N, x, y )
    integer, parameter :: RK = kind(0.0e0)
    include 'WriteLookUpTable.f9h'
  end subroutine WriteLookUpTable_r4

  subroutine WriteLookUpTable_r8 ( filename, N, x, y )
    integer, parameter :: RK = kind(0.0d0)
    include 'WriteLookUpTable.f9h'
  end subroutine WriteLookUpTable_r8

!==================== Private Procedures ===============================
! ...........................................  BattleRealToInt_r4  .....
  ! This is a utility function to convert a real-valued function into the 
  ! integer-valued function Battleship expects
  function BattleRealToInt_r4 ( arg, arg1, delta, setUp ) result( iarg )
    integer, parameter :: RK = R4
    ! Args
    real(rk), intent(in)           :: arg
    real(rk), optional, intent(in) :: arg1, delta
    logical, optional, intent(in)  :: setUp
    integer                        :: iarg
    ! Internal variables
    logical                        :: mySetup
    real(rk), save                 :: x1, eps
    ! Executable
    mySetup = .false.
    if ( present(setUp) ) mySetup = setUp
    if ( mySetUp ) then
      x1 = arg1
      eps = delta
    else
      iarg = int( (arg - x1) / eps ) + 1
    endif
  end function BattleRealToInt_r4

! ...........................................  BattleRealToInt_r8  .....
  function BattleRealToInt_r8 ( arg, arg1, delta, setUp ) result( iarg )
    integer, parameter :: RK = r8
    ! Args
    real(rk), intent(in)           :: arg
    real(rk), optional, intent(in) :: arg1, delta
    logical, optional, intent(in)  :: setUp
    integer                        :: iarg
    ! Internal variables
    logical                        :: mySetup
    real(rk), save                 :: x1, eps
    ! Executable
    mySetup = .false.
    if ( present(setUp) ) mySetup = setUp
    if ( mySetUp ) then
      x1 = arg1
      eps = delta
    else
      iarg = int( (arg - x1) / eps ) + 1
    endif
  end function BattleRealToInt_r8

! ..............................................  createXArray_r4  .....
  ! This family creates an array of x values appropriate
  ! for the LkUpTable
  subroutine CreateXArray_r4( xArray, LkUpTable )
    integer, parameter :: RK = R4
    ! Args
    real(rk), dimension(:), pointer      :: xArray
    type(LookUpTable_0(r4)), intent(in)       :: LkUpTable
    ! Internal variables
    integer :: i
    ! Executable
    call allocate_test( xArray, LkUpTable%N, 'xArray (r4)', ModuleName )
    do i=1, LkUpTable%N
      xArray(i) = LkUpTable%x1 + (i-1)*(LkUpTable%x2-LkUpTable%x1)/(LkUpTable%N-1)
    enddo
  end subroutine CreateXArray_r4

! ..............................................  createXArray_r8  .....
  subroutine CreateXArray_r8( xArray, LkUpTable )
    integer, parameter :: RK = R8
    ! Args
    real(rk), dimension(:), pointer      :: xArray
    type(LookUpTable_0(r8)), intent(in)       :: LkUpTable
    ! Internal variables
    integer :: i
    ! Executable
    call allocate_test( xArray, LkUpTable%N, 'xArray (r8)', ModuleName )
    do i=1, LkUpTable%N
      xArray(i) = LkUpTable%x1 + (i-1)*(LkUpTable%x2-LkUpTable%x1)/(LkUpTable%N-1)
    enddo
  end subroutine CreateXArray_r8

! .................................................  psimpsons_r4  .....
  ! This family of functions performs a part of asimpson's rule integration
  ! from x1 to x
  function Psimpsons_r4 ( x, x1, x2, h, y ) result (sum)
    integer, parameter :: RK = R4
    ! Args
    real(rk), intent(in)               :: x
    real(rk), intent(in)               :: x1
    real(rk), intent(in)               :: x2
    real(rk), intent(in)               :: h
    real(rk), dimension(:), intent(in) :: y
    real(rk)                           :: sum
    ! Internal variables
    integer :: i
    integer, dimension(2) :: indices
    integer :: n
    real(rk), dimension(size(y))      :: xArray
    real(rk)                           :: yofx
    ! Executable
    sum = 0.
    if ( size(y) < 2 ) return
    n = (x2 - x1) / h + 1.5
    do i=1, n
      xArray(i) = x1 + (i-1)*h
    enddo
    ! 1st--we need to find how many intervals, n, are to the left of x
    call ClosestElement( x, xArray, indices, 'l' )
    n = indices(1)
    sum = simpsons( n, h, y )
    yofx = UseLookUpTable ( x, y, x1, x2, xArray )
    ! Now integrate over the distance from x[n] to x
    ! using, less accurately, a trapezoidal quadrature
    sum = sum + ( (x - xArray(n))/2 * ( y(n) + yofx ) )
  end function Psimpsons_r4

! .................................................  psimpsons_r8  .....
  function Psimpsons_r8 ( x, x1, x2, h, y ) result (sum)
    integer, parameter :: RK = R8
    ! Args
    real(rk), intent(in)               :: x
    real(rk), intent(in)               :: x1
    real(rk), intent(in)               :: x2
    real(rk), intent(in)               :: h
    real(rk), dimension(:), intent(in) :: y
    real(rk)                           :: sum
    ! Internal variables
    integer :: i
    integer, dimension(2) :: indices
    integer :: n
    real(rk), dimension(size(y))      :: xArray
    real(rk)                           :: yofx
    ! Executable
    sum = 0.
    if ( size(y) < 2 ) return
    n = (x2 - x1) / h + 1.5
    do i=1, n
      xArray(i) = x1 + (i-1)*h
    enddo
    ! 1st--we need to find how many intervals, n, are to the left of x
    call ClosestElement( x, xArray, indices, 'l' )
    n = indices(1)
    sum = simpsons( n, h, y )
    yofx = UseLookUpTable ( x, y, x1, x2, xArray )
    ! Now integrate over the distance from x[n] to x
    ! using, less accurately, a trapezoidal quadrature
    sum = sum + ( (x - xArray(n))/2 * ( y(n) + yofx ) )
  end function Psimpsons_r8

! ...................................................  reposit_r4  .....
  ! This family repositions x if it is outside [x1, x2]
  ! considering the type of BC determined by the LkUpTable
  subroutine Reposit_r4 ( x, LkUpTable, p, itsSign )
    integer, parameter :: RK = R4
    ! Args
    real(rk), intent(in)                 :: x ! Given this x
    type(LookUpTable_0(r4)), intent(in)       :: LkUpTable
    real(rk), intent(out)                :: p ! reposition it here
    integer, intent(out)                 :: itsSign
    ! Internal variables
    real(rk) :: d
    real(rk) :: Period ! The "period"
    integer :: k ! How many periods between home and x
    ! Executable
    p = x
    itsSign = 1
    if (  LkUpTable%x1 <= x .and. x <= LkUpTable%x2 ) return
    select case ( LkUpTable%BC )
    case ( 'cyclic', 'scyclic' )
      ! How far from "home" are we?
      ! We'll attempt to write this as k*Period + p
      d = max( LkUpTable%x1 - x, x - LkUpTable%x2 )
      Period = LkUpTable%x2 - LkUpTable%x1
      k = d / Period + 0.9999
      if ( x < LkUpTable%x1 ) then
        p = x + k*Period
      else
        p = x - k*Period
      end if
      ! Now if signed cyclic (any non-trig examples of such a beast?)
      ! we also need to know if k is even or odd
      ! if even, we're still the same sign as home
      ! if odd, need to reverse sign
      k = abs(k)
      if ( mod(k, 2) /= 0 .and. LkUpTable%BC == 'scyclic' ) itsSign = -1
    case ( 'express' )
      ! The caller will utilize LkUpTable%y(Left)(Right)
    case default
    ! case ( 'clamped' )
      p = max( p, LkUpTable%x1 )
      p = min( p, LkUpTable%x2 )
    end select
  end subroutine Reposit_r4

! ...................................................  reposit_r8  .....
  subroutine Reposit_r8( x, LkUpTable, p, itsSign )
    integer, parameter :: RK = R8
    ! Args
    real(rk), intent(in)                 :: x ! Given this x
    type(LookUpTable_0(r8)), intent(in)       :: LkUpTable
    real(rk), intent(out)                :: p ! reposition it here
    integer, intent(out)                 :: itsSign
    ! Internal variables
    real(rk) :: d
    real(rk) :: Period ! The "period"
    integer :: k ! How many periods between home and x
    ! Executable
    p = x
    itsSign = 1
    if (  LkUpTable%x1 <= x .and. x <= LkUpTable%x2 ) return
    select case ( LkUpTable%BC )
    case ( 'cyclic', 'scyclic' )
      ! How far from "home" are we?
      ! We'll attempt to write this as k*Period + p
      d = max( LkUpTable%x1 - x, x - LkUpTable%x2 )
      Period = LkUpTable%x2 - LkUpTable%x1
      k = d / Period + 0.9999
      if ( x < LkUpTable%x1 ) then
        p = x + k*Period
      else
        p = x - k*Period
      end if
      ! Now if signed cyclic (any non-trig examples of such a beast?)
      ! we also need to know if k is even or odd
      ! if even, we're still the same sign as home
      ! if odd, need to reverse sign
      k = abs(k)
      if ( mod(k, 2) /= 0 .and. LkUpTable%BC == 'scyclic' ) itsSign = -1
    case ( 'express' )
      ! The caller will utilize LkUpTable%y(Left)(Right)
    case default
    ! case ( 'clamped' )
      p = max( p, LkUpTable%x1 )
      p = min( p, LkUpTable%x2 )
    end select
  end subroutine Reposit_r8

! ..................................................  simpsons_r4  .....
  ! This family of functions performs simpson's rule integrations
  ! for either even or odd n
  function Simpsons_r4 ( n, h, y ) result (sum)
    integer, parameter :: RK = R4
    ! Args
    integer, intent(in)                :: n
    real(rk), intent(in)               :: h
    real(rk), dimension(:), intent(in) :: y
    real(rk)                           :: sum
    call SimpsonsSub( y, h, n, sum )
  end function Simpsons_r4

! ..................................................  simpsons_r8  .....
  function Simpsons_r8 ( n, h, y ) result (sum)
    integer, parameter :: RK = R8
    ! Args
    integer, intent(in)                :: n
    real(rk), intent(in)               :: h
    real(rk), dimension(:), intent(in) :: y
    real(rk)                           :: sum
    call SimpsonsSub( y, h, n, sum )
  end function Simpsons_r8

! ....................................................  D_CSpline  .....
  ! This family was moved here from fwdmdl
  subroutine D_CSpline (xin, xout, yin, yout, nin, nout, ymin, ymax)
    integer, parameter :: rk = kind(0.0d0)
    include 'cspline.f9h'
  end subroutine D_CSpline

! ....................................................  s_cspline  .....
  subroutine S_CSpline (xin, xout, yin, yout, nin, nout, ymin, ymax)
    integer, parameter :: rk = kind(0.0e0)
    include 'cspline.f9h'
  end subroutine S_CSpline

! ......................................................  D_PCSPL  .....
  subroutine D_PCSpl ( tau, c, n, ibcbeg, ibcend )
    integer, parameter :: RK = kind(0.0d0)
    include 'pcspl.f9h'
  end subroutine D_PCSPL

! ......................................................  S_PCSPL  .....
  subroutine S_PCSPL ( TAU, C, N, IBCBEG, IBCEND )
    integer, parameter :: RK = kind(0.0e0)
    include 'pcspl.f9h'
  end subroutine S_PCSPL

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, id ! .mod files sometimes change if PRINT is added
  end function not_used_here

end module MLSNumerics
!=============================================================================

!
! $Log$
! Revision 2.101  2022/03/24 21:43:07  pwagner
! Added AGM function; renamed UDF dummy arg LkUpTable
!
! Revision 2.100  2018/12/05 01:00:46  pwagner
! Made Destroy, Dump, and SetUp typebound procedures
!
! Revision 2.99  2018/12/03 23:19:58  pwagner
! Changed name of datatype to more natural LookUpTable
!
! Revision 2.98  2017/12/07 02:22:11  vsnyder
! Remove unused parameter declaration
!
! Revision 2.97  2017/11/03 23:34:06  pwagner
! Fixed error where we omitted module procedure FindInRange_2d versions
!
! Revision 2.96  2017/11/03 19:56:36  pwagner
! Most array gymnastics moved from MLSFillValues to HyperSlabs module
!
! Revision 2.95  2017/11/02 00:09:38  pwagner
! Added Read,WriteLookupTable, 2d versions of FindInRange, and F_Of_X
!
! Revision 2.94  2017/10/31 23:46:29  vsnyder
! Make Coefficients and LookUpTable_0 parameterized types
!
! Revision 2.93  2017/10/17 23:41:31  pwagner
! Removed unused stuff
!
! Revision 2.92  2016/09/14 20:13:19  vsnyder
! Get PureHunt from Pure_Hunt_m
!
! Revision 2.91  2016/08/23 20:27:40  pwagner
! InterpolateArraySetup may return after failure if optional arg fail present
!
! Revision 2.90  2016/07/28 01:40:24  vsnyder
! Remove unused USE
!
! Revision 2.89  2016/06/02 02:11:51  vsnyder
! Add Average function
!
! Revision 2.88  2015/07/29 00:24:05  vsnyder
! Add InterpolateExtrapolate_[sd]_1
!
! Revision 2.87  2015/05/27 22:42:15  vsnyder
! Move Hunt-related stuff to Hunt_m to eliminate circular dependence
!
! Revision 2.86  2015/04/29 00:53:01  vsnyder
! Make some specific procedures public
!
! Revision 2.85  2015/04/11 01:28:25  vsnyder
! Add 'Second' argument to several routines
!
! Revision 2.84  2015/04/07 02:47:30  vsnyder
! Add InterpolateExtrapolate
!
! Revision 2.83  2015/03/28 01:49:22  vsnyder
! Moved Cross to Cross_m module
!
! Revision 2.82  2014/01/09 00:24:29  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.81  2013/08/13 00:58:43  vsnyder
! Move SolveQuadratic into MLSNumerics
!
! Revision 2.80  2013/08/12 23:47:25  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.79  2013/05/31 23:30:37  vsnyder
! Add two-argument cross product routines
!
! Revision 2.78  2013/05/31 02:37:11  vsnyder
! Add cross product
!
! Revision 2.77  2013/04/12 00:35:56  vsnyder
! Describe 'index' argument of Hunt more precisely
!
! Revision 2.76  2013/02/11 17:19:04  pwagner
! Battleship returns status if cant find root
!
! Revision 2.75  2012/12/20 01:06:15  vsnyder
! Add Interpolate_Regular_To_Irregular
!
! Revision 2.74  2012/06/12 18:10:00  pwagner
! Battleship can now find non-integer roots
!
! Revision 2.73  2012/05/25 20:53:21  pwagner
! Added multidimensional LinearInterpolate and HuntBox
!
! Revision 2.72  2012/04/20 23:55:22  pwagner
! Remove unused code, misleading comments
!
! Revision 2.71  2011/11/18 02:42:35  vsnyder
! Add Interpolate_2d_Composite
!
! Revision 2.70  2011/08/26 17:52:49  pwagner
! purehunt recovers optimized functionality of fwdmdls own hunt
!
! Revision 2.69  2011/08/26 00:23:56  pwagner
! Moved Simpson and CSpline functionality here from fwdmdl
!
! Revision 2.68  2011/08/20 00:47:14  vsnyder
! use IEEE_Arithmetic to get IEEE_Is_NaN
!
! Revision 2.67  2011/08/17 00:48:57  pwagner
! Fixed bug in Simpsons; made it public
!
! Revision 2.66  2011/03/22 23:37:56  pwagner
! Rerank now taken from MLSFillValues module
!
! Revision 2.65  2010/06/07 23:31:49  vsnyder
! Add BivariateLinearInterpolation
!
! Revision 2.64  2009/12/08 21:43:14  vsnyder
! Get Symm_Tri
!
! Revision 2.63  2009/11/17 23:34:38  vsnyder
! Alphabetize public procedure names, add components for periodic continuity
! to Coefficients_R* types, add Name argument to DumpCoefficients_r*
!
! Revision 2.62  2009/06/20 02:32:58  vsnyder
! Precompute more stuff, handle identical abscissae, in spline case
!
! Revision 2.61  2009/06/13 02:27:39  vsnyder
! Several coefficients_R8 should have been coefficients_R4
!
! Revision 2.60  2008/09/03 20:43:48  pwagner
! Added FindInRange
!
! Revision 2.59  2008/06/06 22:52:21  pwagner
! EssentiallyEqual moved to MLSFillValues
!
! Revision 2.58  2008/05/02 00:41:42  vsnyder
! Delete unused symbol
!
! Revision 2.57  2008/01/07 21:36:33  pwagner
! Replace DEFAULTUNDEFINEDVALUE with user-settable undefinedValue
!
! Revision 2.56  2007/08/20 22:03:31  pwagner
! Fixed another bug in HuntRange
!
! Revision 2.55  2007/08/13 17:28:46  pwagner
! Fixed obvious bugs in HuntRange
!
! Revision 2.54  2007/08/07 23:55:02  pwagner
! Added new data type, functions for approximating
!
! Revision 2.53  2007/07/31 22:48:27  pwagner
! UseLookUpTable can now differentiate, integrate
!
! Revision 2.52  2007/07/25 20:09:25  vsnyder
! Delete USE for unused entity
!
! Revision 2.51  2007/07/23 23:18:26  pwagner
! Battleship may be used with logical-valued function
!
! Revision 2.50  2007/06/21 00:49:52  vsnyder
! Remove tabs, which are not part of the Fortran standard
!
! Revision 2.49  2007/04/02 22:53:26  pwagner
! Added Battleship integer rootfinder
!
! Revision 2.48  2007/03/14 23:58:05  pwagner
! Improved precision when interpolating
!
! Revision 2.47  2007/03/02 18:21:14  pwagner
! Added LookUpTable routines
!
! Revision 2.46  2006/10/04 03:20:08  vsnyder
! Better comments for HUNT results
!
! Revision 2.45  2006/08/05 02:36:36  vsnyder
! Delete unused symbols
!
! Revision 2.44  2006/08/03 01:57:22  vsnyder
! Return an optional status flag from Hunt
!
! Revision 2.43  2006/01/14 00:53:05  pwagner
! Added HuntRange
!
! Revision 2.42  2006/01/05 03:46:47  vsnyder
! Add Interp_Bilinear_2d_1d_r*
!
! Revision 2.41  2006/01/05 00:56:03  pwagner
! Added ClosestElement for multidimensional, non-monotonic Hunting
!
! Revision 2.40  2005/12/16 00:02:05  pwagner
! FillValue-related stuff moved to new MLSFillValues module
!
! Revision 2.39  2005/11/04 18:47:34  pwagner
! Added IsMonotonic; Monotonize made public
!
! Revision 2.38  2005/08/15 20:35:37  pwagner
! Added Monotonize interfaces
!
! Revision 2.37  2005/08/06 01:36:30  vsnyder
! Add Dump for interpolation coefficients
!
! Revision 2.36  2005/08/05 20:34:47  pwagner
! ReplaceFillValues can now to interpolate to bridge across MissingValues
!
! Revision 2.35  2005/08/03 16:36:46  pwagner
! antiderivatives; added replaceFillValues
!
! Revision 2.34  2005/06/22 17:25:50  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.33  2005/05/12 20:47:56  pwagner
! Added new versions of EssentaillyEqual discounting Fills, NaNs
!
! Revision 2.32  2005/01/19 17:16:48  pwagner
! Now gets dump from dump_0
!
! Revision 2.31  2004/09/28 23:14:24  pwagner
! Added isFillValue function
!
! Revision 2.30  2004/09/10 23:52:29  livesey
! Added the logSpace options for Hunt (not actually needed but never
! mind).
!
! Revision 2.29  2003/09/11 23:09:18  livesey
! Added skipNewY argument
!
! Revision 2.28  2003/04/04 00:10:08  livesey
! Added EssentiallyEqual
!
! Revision 2.27  2002/11/25 18:51:19  vsnyder
! More interfaces
!
! Revision 2.26  2002/11/23 00:01:00  vsnyder
! Modify interpolation to separate setup and teardown
!
! Revision 2.25  2002/10/08 00:09:12  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.24  2002/10/04 16:40:30  pwagner
! Fixed missing close quotes on include lines
!
! Revision 2.23  2002/10/04 01:48:27  vsnyder
! Get rid of a local variable
!
! Revision 2.22  2002/10/04 00:48:05  vsnyder
! Move declarations of local variables to includes.  Use generic
! InterpolateValues in InterpolateScalar routines.  Move guts of
! InterpolateScalar routines to an include file.
!
! Revision 2.21  2002/09/13 18:08:12  pwagner
! May change matrix precision rm from r8
!
! Revision 2.20  2002/09/11 17:43:38  pwagner
! Began changes needed to conform with matrix%values type move to rm from r8
!
! Revision 2.19  2002/05/24 17:00:55  livesey
! Fixed bug with interpolation with only one oldX value
!
! Revision 2.18  2001/12/01 01:03:07  livesey
! Bug fix with derivatives.
!
! Revision 2.17  2001/11/14 01:47:40  livesey
! Added nearest option to Hunt
!
! Revision 2.16  2001/10/19 23:41:43  livesey
! Fixed bug with dyByDx extrapolation
!
! Revision 2.15  2001/09/24 17:27:50  pwagner
! Removed random number things
!
! Revision 2.14  2001/09/20 20:54:55  pwagner
! Added random number stuff from MATH77
!
! Revision 2.13  2001/07/10 17:55:23  livesey
! Cosmetic changes only
!
! Revision 2.12  2001/07/06 18:44:40  dwu
! Add a feature in InterpolateScaler to handle periodic data
!
! Revision 2.11  2001/06/07 21:59:41  pwagner
! Added Copyright statement
!
! Revision 2.10  2001/05/08 23:25:58  livesey
! Added a nullify for oldSecond (how did I miss that?)
!
! Revision 2.9  2001/05/03 23:13:28  livesey
! Made Van's changes compile with NAG
!
! Revision 2.8  2001/05/03 21:54:49  vsnyder
! Added some nullify's, did some cosmetic changes
!
! Revision 2.7  2001/04/28 19:42:48  livesey
! Hunt now correctly handles cases where list is of size one.
! Changed to lower case keywords and reformatted.
!
! Revision 2.6  2001/04/28 07:05:05  livesey
! Minor bug fix in spline
!
! Revision 2.5  2001/04/11 22:43:19  vsnyder
! Let sparsify do the deallocate_test
!
! Revision 2.4  2001/03/06 00:35:23  livesey
! Missed one pointer nullification
!
! Revision 2.3  2001/03/05 01:20:36  livesey
! Nullified pointers
!
! Revision 2.2  2001/02/22 01:59:52  vsnyder
! Remove declarations for unused variables and parameters
!
! Revision 2.1  2001/02/09 00:38:55  livesey
! Various changes
!
! Revision 2.0  2000/09/05 17:41:06  dcuddy
! Change revision to 2.0
!
! Revision 1.13  2000/06/23 01:08:48  vsnyder
! Delete unused variables (except ID) to keep NAG f95 happy
!
