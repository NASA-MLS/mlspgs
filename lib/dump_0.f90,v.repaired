! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Dump_0

! Low-level dump routines -- for some arrays of intrinsic type.
! It should handle most combinations of rank and type
! Its behavior depends on optional parameters
! The actual output device determined by output_m module
!
! In general, dumping a whole array of values will be presented as a matrix
! up to 10 values across
! Instead of a whole array, or in addition, one may dump a condensed summary
! showing min, max, percentages of non-zero values, etc.

! We take special care when dumping arrays whose values are
! logical-valued
!   bit-valued
!   character strings
!   rank > 1 (only for integer, single, double, or complex)

! See also dump_1.f90

  use HyperSlabs, only: Bandwidth, Collapse
  use BitStuff, only: MaxBitNumber, WhichBitsAreSet
  use Dump_Options, only: AfterSub, AuBrick, MyBandwidth=>Bandwidth, &
    & Clean, CollapseIt, CollapseOptions, &
    & DefaultDumpOptions, DefaultMaxLon, DefaultPCTFormat, &
    & DefaultWidth, Dopt_Collapse, Dopt_Transpose, Dot, DiffRMSMeansRMS, &
    & DontDumpIfAllEqual, Dopts, Gaps, IntPlaces, Laconic, MaxNumNANs, &
    & NameHasBeenPrinted, NameOnEachLine, NaNs, &
    & PCTFormat, PrintFillValue, PrintName, PrintNameAtLineEnd, &
    & Ratios, RMS, RMSFormat, &
    & SDFormatDefault, SDFormatDefaultCmplx, Stats, StatsOnOneLine, &
    & TheDumpBegins, TheDumpEnds, ItsShape, MyTranspose=>Transpose, &
    & TrimIt, Unique, WholeArray
  use HighOutput, only: BlanksToColumn, NumNeedsFormat, OutputNamedValue
  use MLSFillValues, only: InfFunction, IsFinite, IsInfinite, IsNaN, NaNFunction, &
    & WhereAreTheInfs, WhereAreTheNaNs
  use MLSFinds, only: FindAll, FindUnique
  use MLSStats1, only: FillValueRelation, &
    & MLSMax, MLSMean, MLSMin, MLSStdDev
  use MLSStringLists, only: CatLists, NumStringElements, OptionDetail
  use MLSStrings, only: Delete, Indexes, ReadIntsFromChars
  use Optional_m, only: Default
  use Output_m, only: OutputOptions, &
    & Blanks, GetOutputStatus, Newline, Output

  implicit none
  private

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (subroutines and functions)
! AsExplicitValues         Dump an array in a form suitable for pasting into
!                            an l2cf as values for an explicit Fill; e.g.
!                            [ 258.0K, 312.1K, .. ]
! Dump                     Dump an array to output
! Dump_2x2xN               Dump polarized incremental optical depth
! Empty                    Print mesg that array is empty
! FinishLine               Print a newLine, optionally echoing the item name
! ILog10                   Returns integer part of base-10 log of an int
! Name_And_Size            Print array name and its shape

! PrintName                Print the item name unless already done so
! PrintRMSetc              Prints a nicely-formatted list of min, max, etc
! === (end of toc) ===

! === (start of api) ===
! AsExplicitValues  ( Array,
!      [char* unit,], [char separator], [char brackets[2] )
! Dump ( Array, char* name,
!      [fillvalue], [int width], [char* format],
!      [int lbound], [char* options] ) 
!       where array can be a 1, 2, or 3d array of
!       chars, ints, reals, or doubles,
!       and fillValue is a scalar of the same type, if present
! Dump_2x2xN ( cmplx array(2, 2, :) &
!      [, char* name] [, char* Format] [, char* Options] )
! Empty ( [ char* Name ] )
! FinishLine
! int ILog10 ( int int )
! PrintName ( char* Name [, log nameHasBeenPrintedAlready] )
! PrintRMSetc ( char* Name, num min, num max, num rms, num mean )
!
! Options are described in Dump_Options
! By their judicious use you might choose to print an entire array or
! instead just a statistical sampling.
!
! Note that the appearance of Name may be affected by the use
! of any of 3 possible flags in the Name field of Dump. The method
! of appending the flags has an air of hackery-quackery about it. E.g.,
! If Name = 'MyName\h'
! then Dump might print this
! ---------------------------MyName---------------------
!    (:) all 3500 values are 0
! So the flags are separated from the final character of Name by a '\'
! flag             effect
!  b               print Name in a Banner
!  h               print Name as a Headline
!  n               omit Newline after printing name
! === (end of api) ===

  public :: &
    & AsExplicitValues, &
    & Dump, Dump_Sparse, Dump_2x2xN, Empty, FinishLine, ILog10, Name_And_Size, &
    & PrintName, PrintRMSEtc

  ! =====     Public Generics     ======================================

  interface Dump        ! Dump n-d arrays of homogeneous type
    module procedure Dump_1D_BIT, Dump_1D_Char, Dump_1D_Complex, Dump_1D_DComplex
    module procedure Dump_1D_Double, Dump_1D_Integer, Dump_1D_Integer_2B
    module procedure Dump_1D_Logical, Dump_1D_Real
    module procedure Dump_2D_Char, Dump_2D_Complex, Dump_2D_DComplex
    module procedure Dump_2D_Double, Dump_2D_Integer, Dump_2D_Integer_2B
    module procedure Dump_2D_Logical, Dump_2D_Real
    module procedure Dump_3D_Char, Dump_3D_Double, Dump_3D_Integer
    module procedure Dump_3D_Real, Dump_3D_Complex, Dump_3D_DComplex
    module procedure Dump_4D_Double, Dump_4D_Real
  end interface

  interface Dump_Sparse ! Dump only the nonzeroes in 1D arrays
    module procedure Dump_1D_Sparse_Double, Dump_1D_Sparse_Real
    module procedure Dump_2D_Row_Sparse_Double, Dump_2D_Row_Sparse_Real
  end interface

  interface Dump_2x2xN  ! For polarized incremental optical depth
    module procedure Dump_2x2xN_Complex, Dump_2x2xN_DComplex
  end interface

  ! =====     Private Generics     =====================================
  interface DumpCollapsedArray
    module procedure DumpCollapsedArray_1D_Double, DumpCollapsedArray_1D_Real
    module procedure DumpCollapsedArray_2D_Double, DumpCollapsedArray_2D_Real
    module procedure DumpCollapsedArray_3D_Double, DumpCollapsedArray_3D_Real
    module procedure DumpCollapsedArray_1D_Integer
    module procedure DumpCollapsedArray_2D_Integer
    module procedure DumpCollapsedArray_3D_Integer
  end interface

  Interface PrintIt
    module procedure PrintIt_Char, PrintIt_Double, PrintIt_Int, PrintIt_Real
    module procedure PrintIt_Complex, PrintIt_DComplex
  end Interface

  Interface PrintRMSEtc
    module procedure PrintRMSEtc_Double, PrintRMSEtc_int, PrintRMSEtc_Real
  end Interface

  Interface Say_Fill
    module procedure Say_Fill_Char, Say_Fill_Double, Say_Fill_Int
    module procedure Say_Fill_Real, Say_Fill_Complex, Say_Fill_DComplex
  end Interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! These public parameters can't be reconfigured outside the module
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------

  ! These are private variables declared module-wide purely for convenience
  integer, parameter          :: MaxNumElements    = 2000 
  logical                     :: DumpTheseZeros           
  character(len=16)           :: MyOptions                
  character(len=16)           :: MyPCTFormat              
  integer                     :: Bwidth                   
  integer                     :: myRank                   
  integer                     :: numFill               
  integer                     :: numNonFill               
  integer                     :: indx2BSliced             
  real                        :: Pctnzero                 
  logical, save               :: ThisIsADiff = .false.    
  integer                     :: How_many                 
  integer, dimension(1024)    :: which                    
  complex, parameter          :: one_c4 = (1., 0.)        
  character(len=*), parameter :: OldNameOnEachLine = ' ' ! Should be in dump_options

contains

  ! -----------------------------------------------  AsExplicitValues  -----
  ! Print lines that can be pasted into an l2cf as part of an
  ! explicit Fill using values of a 1d array.
  ! The defaults are chosen aptly. Change them only if you know
  ! what you're doing.
  subroutine AsExplicitValues ( Array, Unit, Separator, Brackets )
    ! Args
    real, intent(in)                       :: Array(:)
    character(len=*), intent(in), optional :: Unit      ! E.g., 'K'
    character(len=1), intent(in), optional :: Separator ! E.g., ','
    character(len=*), intent(in), optional :: Brackets  ! E.g., '[]'
    ! Internal variables
    integer                                :: Column
    integer                                :: i
    character(len=8)                       :: myUnit
    character(len=1)                       :: mySeparator
    character(len=2)                       :: myBrackets 
    ! Executable
    myUnit = Default( Unit, ' ' )
    mySeparator = Default( Separator, ',' )
    myBrackets = Default( Brackets, '[]' )
    if ( size(array) < 1 ) return
    call output ( myBrackets(1:1), advance='no' )
    call output ( ' $', advance='yes' )
    do i=1, size(array)
      call output( array(i), advance='no' )
      if ( len_trim(myUnit) > 0 ) call output( trim(myUnit), advance='no' )
      if ( i < size(array) ) then
        call output ( mySeparator, advance='no' )
        if ( len_trim(mySeparator) > 0 ) call Blanks( 1 )
      endif
      Column = getOutputStatus( 'column' )
      if ( Column > 72 ) then
        call output ( ' $', advance='yes' )
      endif

    enddo
    if ( getOutputStatus( 'start' ) /= 1 ) call output ( ' $', advance='yes' )
    call output ( myBrackets(2:2), advance='yes' )
    
  end subroutine AsExplicitValues

  ! -----------------------------------------------  Dump_1D_Bit  -----
  subroutine Dump_1D_Bit ( Array, Name, BITNameS, FillValue, Options )
    integer, intent(in) :: Array(:)
    character(len=*), intent(in) :: Name
    character(len=*), intent(in) :: BITNameS
    integer, intent(in), optional :: FillValue
    character(len=*), optional, intent(in) :: options

    integer :: howMany, J, K
    integer, dimension(MAXBITNUMBER) :: ints
    integer :: myFillValue
    integer :: NumBitNames
    integer :: NumZeroRows
    integer, dimension(MAXBITNUMBER+1) :: set
    ! Executable
    call theDumpBegins ( options )
    myFillValue = 0
    if ( present(FillValue) ) myFillValue = FillValue

    NumBitNames = NumStringElements( bitNames, countEmpty=.true. )
    if ( NumBitNames < 1 ) NumBitNames = MAXBITNUMBER+1
    NumBitNames = min( numBitNames, MAXBITNUMBER+1 )
    numZeroRows = 0
    if ( any(shape(array) == 0) ) then
      call empty ( name )
    else
      call name_and_size ( name, dopts(clean)%v, size(array) )
      call newline
      call output( trim(BitNames), advance='yes' )
      do j=1, size(array)
        if ( array(j) == myFillValue ) then
          numZeroRows = numZeroRows + 1
        else
          call WhichBitsAreSet( array(j), set, howMany )
          ints = 0
          do k=1, howMany
            ints( 1 + set(k) ) = 1
          end do
          if ( numZeroRows > 0 ) then
            call output ( ' ' , advance='no' )
            call output ( numZeroRows , advance='no' )
            call output ( ' lines of ', advance='no' )
            call output ( myFillValue , advance='no' )
            call output ( ' not printed', advance='no' )
            call finishLine
          end if
          call output( ints(1:numBitNames), format='(i3)', advance='no' )
          call finishLine
          numZeroRows = 0
        end if
      end do
    end if
    call theDumpEnds
  end subroutine Dump_1D_Bit

  ! -----------------------------------------------  Dump_1D_Char  -----
  subroutine Dump_1D_Char ( Array, Name, FillValue, Width, Options, Maxlon, &
    & TheShape )
    character(len=*), intent(in) :: Array(:)
    character(len=*), intent(in), optional :: Name
    character(len=*), intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), optional, intent(in) :: Options
    integer, intent(in), optional :: Maxlon
    character(len=*), intent(in), optional :: TheShape

    integer :: J, K
    integer :: LON
    integer :: NumZeroRows
    character(len=len(array)) :: MyFillValue
    integer :: MyWidth

    call theDumpBegins ( options )
    myFillValue = ' '
    if ( present(FillValue) ) myFillValue = FillValue
    MyWidth = DefaultWidth ! 10
    if ( present(width) ) MyWidth = width

    lon = len(array(1))
    if ( dopts(trimIt)%v ) lon = maxval(len_trim(array))
    if ( present(maxlon) ) then
      lon = min(lon, maxlon)
    else if ( .not. dopts(trimIt)%v ) then
      lon = min(lon, DEFAULTMaxlon)
    end if

    numZeroRows = 0
    if ( any(shape(array) == 0) ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, dopts(clean)%v, 1, TheShape )
      call output ( array(1)(1:lon), advance='no' )
      call finishLine
    else
      call name_and_size ( name, dopts(clean)%v, size(array) )
      if ( getOutputStatus( 'start' ) /= 1 ) call newLine
      do j = 1, size(array), MyWidth
        DumpTheseZeros = dopts(clean)%v .or. &
          & any(array(j:min(j+2*myWidth-1, size(array))) /= myFillValue) &
          & .or. &
          & size(array) <= myWidth
        if (.not. dopts(clean)%v) then
          if ( DumpTheseZeros ) then
            call say_fill ( (/ j-1, size(array) /), numZeroRows, &
              & myFillValue, inc=1 )
          else
            numZeroRows = numZeroRows + 1
          end if
        end if
        if ( DumpTheseZeros ) then
          do k = j, min(j+MyWidth-1, size(array))
              call output ( ' ' // array(k)(1:lon) // ' ' , advance='no' )
          end do
          call newLine
          numZeroRows = 0
        end if
      end do ! j
      call say_fill ( (/ j-MyWidth, size(array) /), numZeroRows, &
        & myFillValue )
    end if
    call theDumpEnds
  end subroutine Dump_1D_Char

  ! --------------------------------------------  Dump_1D_Complex  -----
  subroutine Dump_1D_Complex ( Array, Name, Width, Format, &
    & FillValue, Lbound, Options, TheShape )
    integer, parameter :: RK = kind(0.0e0)
    complex(rk), intent(in) :: Array(:)
    character(len=*), intent(in), optional :: Name
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    real(rk), intent(in), optional :: FillValue
    integer, intent(in), optional :: Lbound ! Low bound for Array             
    character(len=*), optional, intent(in) :: Options
    character(len=*), intent(in), optional :: TheShape

    integer :: J, K, MyWidth
    character(len=64) :: MyFormat
    complex(rk) :: myFillValue
    integer :: Base
    integer :: NumZeroRows
    ! Executable
    myFormat = sdFormatDefaultCmplx
    call theDumpBegins ( options )
    myFillValue = 0.
    myWidth = 3
    if ( present(width) ) myWidth = width
    if ( present(format) ) myFormat = format
    if ( index(myFormat, '*') > 0 ) &
      & myFormat = numNeedsFormat( one_c4*maxval(abs(array)), format )
    base = 0
    if ( present(lbound) ) base = lbound - 1

    include 'dump1db.f9h'
  end subroutine Dump_1D_Complex

  ! --------------------------------------------  Dump_1D_DComplex  -----
  subroutine Dump_1D_DComplex ( Array, Name, Width, Format, &
    & FillValue, Lbound, Options, TheShape )
    integer, parameter :: RK = kind(0.0d0)
    complex(rk), intent(in) :: Array(:)
    character(len=*), intent(in), optional :: Name
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    real(rk), intent(in), optional :: FillValue
    integer, intent(in), optional :: Lbound ! Low bound for Array             
    character(len=*), optional, intent(in) :: Options
    character(len=*), intent(in), optional :: TheShape

    integer :: J, K, MyWidth
    character(len=64) :: MyFormat
    complex(rk) :: myFillValue
    integer :: Base
    integer :: NumZeroRows
    ! Executable
    myFormat = sdFormatDefaultCmplx
    call theDumpBegins ( options )
    myFillValue = 0.
    myWidth = 3
    if ( present(width) ) myWidth = width
    if ( present(format) ) myFormat = format
    if ( index(myFormat, '*') > 0 ) &
      & myFormat = numNeedsFormat( one_c4*maxval(abs(array)), format )
    base = 0
    if ( present(lbound) ) base = lbound - 1

    include 'dump1db.f9h'
  end subroutine Dump_1D_DComplex

 ! ---------------------------------------------  Dump_1D_Double  -----
  subroutine Dump_1D_Double ( Array, Name, &
    & FillValue, Width, Format, Lbound, Options, TheShape, Unit )
    double precision, intent(in) :: Array(:)
    character(len=*), intent(in), optional :: Name
    double precision, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Lbound ! Low bound for Array
    character(len=*), optional, intent(in) :: Options
    character(len=*), intent(in), optional :: TheShape
    integer, intent(in), optional :: Unit

    integer, dimension(MaxNumElements) :: Counts
    double precision, dimension(MaxNumElements) :: Elements
    double precision :: MyFillValue
    integer :: Base, J, K
    character(len=64) :: MyFormat
    integer :: MyWidth
    integer :: NumZeroRows
    integer :: nUnique
    integer :: SU                ! Save unit
    su = outputOptions%prUnit
    if ( present(unit) ) outputOptions%prUnit = unit
    myFillValue = 0.
    if ( present(FillValue) ) myFillValue=FillValue
    myFormat = sdFormatDefault
    include 'dump1d.f9h'
    include 'dump1db.f9h'
    outputOptions%prUnit = su
  end subroutine Dump_1D_Double

  ! --------------------------------------------  Dump_1D_Integer  -----
  subroutine Dump_1D_Integer ( Array, Name, &
    & FillValue, Format, Width, Lbound, Options, TheShape, Unit )
    integer, intent(in) :: Array(:)
    character(len=*), intent(in), optional :: Name
    integer, intent(in), optional :: FillValue
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Width ! How many numbers per line (10)?
    integer, intent(in), optional :: Lbound ! Low bound for Array
    character(len=*), optional, intent(in) :: Options
    character(len=*), intent(in), optional :: TheShape
    integer, intent(in), optional :: Unit

    integer, dimension(MaxNumElements) :: counts
    integer, dimension(MaxNumElements) :: elements
    integer :: MyFillValue
    integer :: Base, J, K
    character(len=64) :: MyFormat
    integer :: MyWidth
    integer :: NumZeroRows
    integer :: nUnique
    integer :: SU                ! Save unit
    if ( present( Format ) ) call outputNamedValue ( 'format', format )
    su = outputOptions%prUnit
    if ( present(unit) ) outputOptions%prUnit = unit
    myFillValue = 0.
    if ( present(FillValue) ) myFillValue=FillValue
    myFormat = 'places=' // INTPLACES ! To sneak places arg into call to output
    include 'dump1d.f9h'
    include 'dump1db.f9h'
    outputOptions%prUnit = su
  end subroutine Dump_1D_Integer

  ! --------------------------------------------  Dump_1D_Integer_2B  -----
  subroutine Dump_1D_Integer_2B ( Array, Name, &
    & FillValue, Format, Width, Lbound, Options, TheShape, Unit )
    use ISO_C_BINDING, only: C_int16_t
    integer(C_int16_t), intent(in) :: Array(:)
    character(len=*), intent(in), optional :: Name
    integer, intent(in), optional :: FillValue
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Width ! How many numbers per line (10)?
    integer, intent(in), optional :: Lbound ! Low bound for Array
    character(len=*), optional, intent(in) :: Options
    character(len=*), intent(in), optional :: TheShape
    integer, intent(in), optional :: Unit

    call dump( int(array), Name, &
    & FillValue, Format, Width, Lbound, Options, TheShape, Unit )
  end subroutine Dump_1D_Integer_2B

  ! ----------------------------------------------  Dump_1D_Logical ----
  subroutine Dump_1D_Logical ( Array, Name, Lbound, Options, TheShape )
    logical, intent(in) :: Array(:)
    character(len=*), intent(in), optional :: Name
    integer, intent(in), optional :: Lbound ! Low bound of Array
    character(len=*), optional, intent(in) :: Options
    character(len=*), intent(in), optional :: TheShape

    integer :: Base, J, K, N, Ntrue, Nfalse
    integer, dimension(size(Array))        :: whereTrue, whereFalse

    ! Executable
    call theDumpBegins ( options )
    base = 0
    if ( present(lbound) ) base = lbound - 1

    if ( any(shape(array) == 0) ) then
      call empty ( name ) ! The array appears to have no elements at all
    else if ( size(array) == 1 .and. base == 0 ) then
      ! The array appears to have only a single element
      call name_and_size ( name, dopts(clean)%v, 1 )
      call output ( array(1), advance='no' )
      call finishLine
    else if ( all(array .eqv. .true.) ) then
      ! The array appears to have only true elements
      call name_and_size ( name, dopts(clean)%v, size(array) )
      call output ( '(' )
      if ( present(lbound) ) call output ( lbound )
      call output ( ':) all ', advance='no' )
      call output( trim(arrayShapeToString(shape(array))), advance='no' )
      call output( ' values are T', advance='no' )
      call finishLine
    else if ( all(array .eqv. .false.) ) then
      ! The array appears to have only false elements
      call name_and_size ( name, dopts(clean)%v, size(array) )
      call output ( '(' )
      if ( present(lbound) ) call output ( lbound )
      call output ( ':) all ', advance='no' )
      call output( trim(arrayShapeToString(shape(array))), advance='no' )
      call output( ' values are F', advance='no' )
      call finishLine
    else if ( dopts(NaNs)%v ) then
      ! Display the indices where true and the indices where false
      call PrintName ( Name, NameHasBeenPrinted )
      call FindAll( array, whereTrue, Ntrue, whereFalse )
      Nfalse = size(array) - Ntrue
      NameHasBeenPrinted = .false.
      call dump( whereTrue(1:Ntrue), 'Indices where TRUE' )
      call dump( whereFalse(1:Nfalse), 'Indices where FALSE' )
    else if ( dopts(gaps)%v ) then
      ! Display gaps, i.e. only the exceptional elements;
      ! T among a host of otherwise F's
      ! or F among a host of otherwise T's
      call name_and_size ( name, dopts(clean)%v, size(array) )
      if ( getOutputStatus( 'start' ) /= 1 ) call newLine
      k = 0
      if ( size(array)/100 > 100 )  call showColumnNums( 100 )
      do j=1, size(array), 100
        N = min(k+100, size(array)) - k
        call output( k+1 , advance='no' )
        call output ( ' through ', advance='no' )
        call output( k+N, advance='no' )
        call finishLine
        ntrue = count( array(k+1 : min(k+100, size(array))) )
        nfalse = n - ntrue
        ! Special cases if all true or all false
        if ( nfalse < 1 ) then
          call output( 'all values are T', advance='no' )
        elseif ( ntrue < 1 ) then
          call output( 'all values are F', advance='no' )
        elseif ( ntrue > 0 .and. ( ntrue < nfalse .or. nfalse < 1 ) ) then
          call output( array(k+1 : k+N), onlyif=.true., advance='no' )
        else
          call output( array(k+1 : k+N), onlyif=.false., advance='no' )
        end if
        call finishLine
        k = k + 100
      end do
      call showColumnNums( 100 )
    else
      ! Run-of-the-mill matrix dumps of 'T' and 'F'
      call name_and_size ( name, dopts(clean)%v, size(array) )
      if ( getOutputStatus( 'start' ) /= 1 ) call newLine
      do j = 1, size(array), 34
        if (.not. dopts(clean)%v) then
          call output ( j+base, max(4,ilog10(size(array))+1) , advance='no' )
          call output ( afterSub , advance='no' )
        end if
        if ( dopts(dot)%v ) then ! print .false. as a dot
          do k = j, min(j+33, size(array))
            call output ( merge('T','.',array(k)) , advance='no' )
          end do
        else
          do k = j, min(j+33, size(array))
            call output ( array(k) , advance='no' )
          end do
        end if
        call newLine
      end do
    end if
    call theDumpEnds
  end subroutine Dump_1D_Logical

  ! -----------------------------------------------  Dump_1D_Real  -----
  subroutine Dump_1D_Real ( Array, Name, &
    & FillValue, Width, Format, Lbound, Options, TheShape, Unit )
    real, intent(in) :: Array(:)
    character(len=*), intent(in), optional :: Name
    real, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Lbound ! Low bound for Array
    character(len=*), optional, intent(in) :: Options
    character(len=*), intent(in), optional :: TheShape
    integer, intent(in), optional :: Unit

    integer, dimension(MaxNumElements) :: Counts
    real, dimension(MaxNumElements) :: Elements
    real :: MyFillValue
    integer :: Base, J, K
    character(len=64) :: MyFormat
    integer :: MyWidth
    integer :: NumZeroRows
    integer :: nUnique
    integer :: SU                ! Save unit
    su = outputOptions%prUnit
    if ( present(unit) ) outputOptions%prUnit = unit
    myFillValue = 0.
    if ( present(FillValue) ) myFillValue=FillValue
    myFormat = sdFormatDefault
    include 'dump1d.f9h'
    include 'dump1db.f9h'
    outputOptions%prUnit = su
  end subroutine Dump_1D_Real

  ! -----------------------------------------------  Dump_2D_Char  -----
  recursive subroutine Dump_2D_Char ( Array, Name, FillValue, Width, Maxlon, &
    & Options, TheShape )
    character(len=*), intent(in) :: Array(:,:)
    character(len=*), intent(in), optional :: Name
    character(len=*), intent(in), optional :: FillValue
    integer, intent(in), optional :: Maxlon
    character(len=*), optional, intent(in) :: options
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: TheShape

    integer :: I, J, K
    integer :: Lon
    integer :: NumZeroRows
    character(len=len(array)) :: MyFillValue
    integer :: MyWidth

    call theDumpBegins ( options )
    if ( dopts(myTranspose)%v ) then
      call dump ( transpose(array), name, &
        & fillvalue, width, &
        & options=snipoption(options, dopt_transpose) )
      return
    end if
    myFillValue = ' '
    if ( present(FillValue) ) myFillValue = FillValue
    MyWidth = DefaultWidth ! 10
    if ( present(width) ) MyWidth = width
    lon = len(array(1,1))
    if ( dopts(trimIt)%v ) lon = maxval(len_trim(array))
    if ( present(maxlon) ) lon = min(lon, maxlon)

    numZeroRows = 0
    if ( any(shape(array) == 0) ) then
      call empty ( name )
    else if ( size( Array,1) == 1 ) then
      call name_and_size ( name, dopts(clean)%v, 1, TheShape )
      call output ( array(1,:), advance='no' )
      call finishLine
    else if ( size( Array,2) == 1 ) then
      call dump ( array(:,1), name, fillValue=fillValue, maxlon=maxlon, options=options )
    else
      call name_and_size ( name, dopts(clean)%v, size(array) )
      if ( getOutputStatus( 'start' ) /= 1 ) call newLine
      do i = 1, size( Array,1)
        do j = 1, size( Array,2), MyWidth
          DumpTheseZeros = dopts(clean)%v .or. &
            & any(array(i,j:min(j+2*MyWidth-1, size( Array,2))) /= myFillValue) &
            & .or. &
            & ( j+MyWidth >= size( Array,2) .and. &
            & any(array(min(i+1, size( Array,1)),1:min(1+MyWidth-1, size( Array,2))) &
            & /= myFillValue) &
            & )
          if (.not. dopts(clean)%v) then
            if ( DumpTheseZeros ) then
              call say_fill ( (/ i, size( Array,1), j-1, size( Array,2) /), &
                & numZeroRows, myFillValue, inc=3 )
            else
              numZeroRows = numZeroRows + 1
            end if
          end if
          if ( DumpTheseZeros ) then
            do k = j, min(j+myWidth-1, size( Array,2))
                call output ( ' ' // array(i,k)(1:lon) // ' ' , advance='no' )
            end do
            call newLine
            numZeroRows = 0
          end if
        end do ! j
      end do ! i
      call say_fill ( (/ i-1, size( Array,1), j-MyWidth, size( Array,2) /), &
        & numZeroRows, myFillValue )
    end if
    call theDumpEnds
  end subroutine Dump_2D_Char

  ! --------------------------------------------  Dump_2D_Complex  -----
  recursive subroutine Dump_2D_Complex ( Array, Name, Width, Format, &
    & FillValue, Options, Lbound, TheShape )
    integer, parameter :: RK = kind(0.0e0)
    complex(rk), intent(in) :: Array(:,:)
    character(len=*), intent(in), optional :: Name
    integer, intent(in), optional :: Width ! How many per line?
    character(len=*), optional :: Format
    real, intent(in), optional :: FillValue
    character(len=*), intent(in), optional :: options
    integer, intent(in), optional :: LBound ! to print for first dimension
    character(len=*), intent(in), optional :: TheShape

    integer :: Base, I, J, K
    integer :: MyWidth
    integer :: NumZeroRows
    complex(rk) :: MyFillValue
    character(len=64) :: MyFormat
    ! Executable

    myFormat = sdFormatDefaultCmplx
    call theDumpBegins ( options )
    if ( dopts(myTranspose)%v ) then
      call dump ( transpose(array), name, &
        & width, format, fillvalue, &
        & options=snipoption(options, dopt_transpose) )
      return
    end if

    myWidth = 3
    if ( present(width) ) myWidth = width

    if ( present(format) ) myFormat = format
    if ( index(myFormat, '*') > 0 ) &
      & myFormat = numNeedsFormat( one_c4*maxval(abs(array)), format )

    myFillValue = 0.0_rk
    if ( present(fillValue) ) myFillValue = fillValue

    base = 1
    if ( present(lbound) ) base = lbound

    include 'dump2db.f9h'
  end subroutine Dump_2D_Complex

  ! --------------------------------------------  Dump_2D_DComplex  -----
  recursive subroutine Dump_2D_DComplex ( Array, Name, Width, Format, &
    & FillValue, Options, Lbound, TheShape )
    integer, parameter :: RK = kind(0.0d0)
    complex(rk), intent(in) :: Array(:,:)
    character(len=*), intent(in), optional :: Name
    integer, intent(in), optional :: Width ! How many per line?
    character(len=*), optional :: Format
    real(rk), intent(in), optional :: FillValue
    character(len=*), intent(in), optional :: options
    integer, intent(in), optional :: LBound ! to print for first dimension
    character(len=*), intent(in), optional :: TheShape

    integer :: Base, I, J, K
    integer :: myWidth
    integer :: NumZeroRows
    complex(rk) :: myFillValue
    character(len=64) :: MyFormat

    ! Executable

    myFormat = sdFormatDefaultCmplx
    call theDumpBegins ( options )
    if ( dopts(myTranspose)%v ) then
      call dump ( transpose(array), name, &
        & width, format, fillvalue, &
        & options=snipoption(options, dopt_transpose) )
      return
    end if

    myWidth = 3
    if ( present(width) ) myWidth = width

    if ( present(format) ) myFormat = format
    if ( index(myFormat, '*') > 0 ) &
      & myFormat = numNeedsFormat( one_c4*maxval(abs(array)), format )

    myFillValue = 0.0_rk
    if ( present(fillValue) ) myFillValue = fillValue

    base = 1
    if ( present(lbound) ) base = lbound

    include 'dump2db.f9h'
  end subroutine Dump_2D_DComplex

  ! ---------------------------------------------  Dump_2D_Double  -----
 recursive subroutine Dump_2D_Double ( Array, Name, &
    & FillValue, Width, Format, Lbound, Options, TheShape, Unit )
    double precision, intent(in) :: Array(:,:)
    character(len=*), intent(in), optional :: Name
    double precision, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Lbound ! to print for first dimension
    character(len=*), intent(in), optional :: Options
    character(len=*), intent(in), optional :: TheShape
    integer, intent(in), optional :: Unit

    integer :: Base, I, J, K
    integer :: NumZeroRows
    double precision :: MyFillValue
    character(len=64) :: MyFormat
    integer :: nUnique
    integer :: MyWidth
    integer, dimension(MaxNumElements) :: Counts
    double precision, dimension(MaxNumElements) :: Elements
    integer :: SU                ! Save unit

    myFormat = sdFormatDefault
    myFillValue = 0.
    if ( present(FillValue) ) myFillValue=FillValue

    su = outputOptions%prUnit
    if ( present(unit) ) outputOptions%prUnit = unit
    base = 1
    if ( present(lbound) ) base = lbound
    include 'dump2d.f9h'
    include 'dump2db.f9h'
    outputOptions%prUnit = su
  end subroutine Dump_2D_Double

  ! --------------------------------------------  Dump_2D_Integer  -----
  recursive subroutine Dump_2D_Integer ( Array, Name, &
    & FillValue, Width, Format, Lbound, Options, TheShape, Unit )
    integer, intent(in) :: Array(:,:)
    character(len=*), intent(in), optional :: Name
    integer, intent(in), optional :: FillValue
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Width ! How many numbers per line (10)?
    integer, intent(in), optional :: Lbound ! to print for first dimension
    character(len=*), intent(in), optional :: Options
    character(len=*), intent(in), optional :: TheShape
    integer, intent(in), optional :: Unit

    integer :: Base, I, J, K
    integer :: MyWidth
    integer :: NumZeroRows
    integer :: myFillValue
    character(len=64) :: MyFormat
    integer :: nUnique
    integer, dimension(MaxNumElements) :: Counts
    integer, dimension(MaxNumElements) :: Elements
    integer :: SU                ! Save unit

    myFormat = 'places=' // INTPLACES ! To sneak places arg into call to output
    myFillValue = 0.
    if ( present(FillValue) ) myFillValue=FillValue

    su = outputOptions%prUnit
    if ( present(unit) ) outputOptions%prUnit = unit
    base = 1
    if ( present(lbound) ) base = lbound

    include 'dump2d.f9h'
    include 'dump2db.f9h'
    outputOptions%prUnit = su
  end subroutine Dump_2D_Integer

  ! --------------------------------------------  Dump_2D_Integer_2B  -----
  recursive subroutine Dump_2D_Integer_2B ( Array, Name, &
    & FillValue, Width, Format, Lbound, Options, TheShape, Unit )
    use ISO_C_BINDING, only: C_int16_t
    integer(C_int16_t), intent(in) :: Array(:,:)
    character(len=*), intent(in), optional :: Name
    integer, intent(in), optional :: FillValue
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Width ! How many numbers per line (10)?
    integer, intent(in), optional :: Lbound ! to print for first dimension
    character(len=*), intent(in), optional :: Options
    character(len=*), intent(in), optional :: TheShape
    integer, intent(in), optional :: Unit

    call dump( int(array), Name, &
    & FillValue, Width, Format, Lbound, Options, TheShape, Unit )
  end subroutine Dump_2D_Integer_2B

  ! --------------------------------------------  Dump_2D_Logical  -----
  recursive subroutine Dump_2D_Logical ( Array, Name, Options, TheShape )
    logical, intent(in) :: Array(:,:)
    character(len=*), intent(in), optional :: Name
    character(len=*), intent(in), optional :: Options
    character(len=*), intent(in), optional :: TheShape

    integer :: I, J, K
    integer, parameter :: MyWidth = 34

    call theDumpBegins ( options )
    if ( dopts(myTranspose)%v ) then
      call dump ( transpose(array), name, &
        & options=snipoption(options, dopt_transpose) )
      return
    end if

    if ( any(shape(array) == 0) ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, dopts(clean)%v, 1, TheShape )
      call output ( array(1,1), advance='no' )
      call finishLine
    else if ( size( Array,2) == 1 ) then
      call dump ( array(:,1), name, options=options )
    else if ( all(array .eqv. .true.) ) then
      call name_and_size ( name, dopts(clean)%v, size(array), TheShape )
      call output ( ': all ', advance='no' )
      call output( trim(arrayShapeToString(shape(array))), advance='no' )
      call output( ' values are T', advance='no' )
      call finishLine
    else if ( all(array .eqv. .false.) ) then
      call name_and_size ( name, dopts(clean)%v, size(array), TheShape )
      call output ( ': all ', advance='no' )
      call output( trim(arrayShapeToString(shape(array))), advance='no' )
      call output( ' values are F', advance='no' )
      call finishLine
    else
      call name_and_size ( name, dopts(clean)%v, size(array) )
      if ( getOutputStatus( 'start' ) /= 1 ) call newLine
      do i = 1, size( Array,1)
        do j = 1, size( Array,2), myWidth
          if (.not. dopts(clean)%v) then
            call output ( i, places=max(4,ilog10(size( Array,1))+1) , advance='no' )
            call output ( j, places=max(4,ilog10(size( Array,2))+1) , advance='no' )
            call output ( afterSub , advance='no' )
          end if
          if ( dopts(dot)%v ) then
            do k = j, min(j+33, size(array,2))
              call output ( merge('T','.',array(i,k)) , advance='no' )
            end do
          else
            do k = j, min(j+33, size(array,2))
              call output ( array(i,k) , advance='no' )
            end do
          end if
          call newLine
        end do ! j
      end do ! i
    end if
    call theDumpEnds
  end subroutine Dump_2D_Logical

  ! -----------------------------------------------  Dump_2D_Real  -----
  recursive subroutine Dump_2D_Real ( Array, Name, &
    & FillValue, Width, Format, Lbound, Options, TheShape, Unit )
    real, intent(in) :: Array(:,:)
    character(len=*), intent(in), optional :: Name
    real, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Lbound
    character(len=*), intent(in), optional :: Options
    character(len=*), intent(in), optional :: TheShape
    integer, intent(in), optional :: Unit

    integer :: Base, I, J, K
    integer :: NumZeroRows
    real :: MyFillValue
    character(len=64) :: MyFormat
    integer :: nUnique
    integer :: MyWidth
    integer, dimension(MaxNumElements) :: Counts
    real, dimension(MaxNumElements) :: Elements
    integer :: SU                ! Save unit

    myFormat = sdFormatDefault
    myFillValue = 0.
    if ( present(FillValue) ) myFillValue=FillValue

    su = outputOptions%prUnit
    if ( present(unit) ) outputOptions%prUnit = unit
    base = 1
    if ( present(lbound) ) base = lbound

    include 'dump2d.f9h'
    include 'dump2db.f9h'
    outputOptions%prUnit = su
  end subroutine Dump_2D_Real

  ! --------------------------------------  Dump_1D_Sparse_Double  -----
  subroutine Dump_1D_Sparse_Double ( Array, Name, Format, Width, Parens, &
                                   & RowNum, Exclude, DidOne )
    double precision, intent(in) :: Array(:)
    character(len=*), intent(in), optional :: Name
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Width   ! Items per line, default 5
    logical, intent(in), optional :: Parens
    integer, intent(in), optional :: RowNum  ! Print this first
    double precision, intent(in), optional :: Exclude ! Don't print these
    logical, intent(out), optional :: DidOne ! Printed a value

    logical :: NeedRow ! Need to print RowNum
    integer :: J, N, W
    character(len=64) :: MyFormat
    logical :: MyParens ! Enclose element number in parens instead of colon after

    myFormat = sdFormatDefault
    if ( present(format) ) myFormat = format
    w = 5
    if ( present(width) ) w = width
    myParens = .false.
    if ( present(parens) ) myParens = parens
    if ( present(name) ) call output ( trim(name), advance='yes' )
    needRow = present(rowNum)

    n = 0 ! Number on the line so far
    if ( present(didOne) ) didOne = .false.
    do j = 1, size(array)
      if ( array(j) == 0 ) cycle
      if ( present(exclude) ) then
        if ( array(j) == exclude ) cycle
      end if
      if ( needRow ) call output ( rowNum, places=4, after='#' )
      needRow = .false.
      if ( n >= w ) then
        call newLine
        call blanks ( 5 )
        n = 0
      end if
      n = n + 1
      if ( myParens ) then
        call output ( j, before=' (', after=') ' )
      else
        call output ( j, places=4, after=': ' )
      end if
      call output ( array(j), format=myFormat )
      if ( present(didOne) ) didOne = .true.
    end do
    if ( n /= 0 ) call newLine

  end subroutine Dump_1D_Sparse_Double

  ! ----------------------------------------  Dump_1D_Sparse_Real  -----
  subroutine Dump_1D_Sparse_Real ( Array, Name, Format, Width, Parens, &
                                 & RowNum, Exclude, DidOne )
    real, intent(in) :: Array(:)
    character(len=*), intent(in), optional :: Name
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Width  ! Items per line, default 5
    logical, intent(in), optional :: Parens
    integer, intent(in), optional :: RowNum ! Print this first
    real, intent(in), optional :: Exclude   ! Don't print these
    logical, intent(out), optional :: DidOne ! Printed a value

    logical :: NeedRow ! Need to print RowNum
    integer :: J, N, W
    character(len=64) :: MyFormat
    logical :: MyParens ! Enclose element number in parens instead of colon after
    needRow = present(rowNum)

    myFormat = sdFormatDefault
    if ( present(format) ) myFormat = format
    w = 5
    if ( present(width) ) w = width
    myParens = .false.
    if ( present(parens) ) myParens = parens
    if ( present(name) ) call output ( trim(name), advance='yes' )

    n = 0 ! Number on the line so far
    if ( present(didOne) ) didOne = .false.
    do j = 1, size(array)
      if ( array(j) == 0 ) cycle
      if ( present(exclude) ) then
        if ( array(j) == exclude ) cycle
      end if
      if ( needRow ) call output ( rowNum, places=4, after='#' )
      needRow = .false.
      if ( n >= w ) then
        call newLine
        call blanks ( 5 )
        n = 0
      end if
      n = n + 1
      if ( myParens ) then
        call output ( j, before=' (', after=') ' )
      else
        call output ( j, places=4, after=': ' )
      end if
      call output ( array(j), format=myFormat )
      if ( present(didOne) ) didOne = .true.
    end do
    if ( n /= 0 ) call newLine

  end subroutine Dump_1D_Sparse_Real

  ! ----------------------------------  Dump_2D_Row_Sparse_Double  -----
  subroutine Dump_2D_Row_Sparse_Double ( Array, Name, Format, Width, Exclude )
    double precision, intent(in) :: Array(:,:)
    character(len=*), intent(in), optional :: Name
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Width ! Items per line, default 5
    double precision, intent(in), optional :: Exclude ! Don't print these

    logical :: DidAnyRow, DidOne
    integer :: I ! Row number

    if ( present(name) ) call output ( trim(name), advance='yes' )
    didAnyRow = .false.
    do i = 1, size(array,1)
      call dump_sparse ( array(i,:), format=format, width=width, parens=.true., &
        & rowNum=i, exclude=exclude, didOne = didOne )
      didAnyRow = didAnyRow .or. didOne
    end do
    if ( .not. didAnyRow ) then
      if ( present(exclude) ) then
        call output ( exclude, &
          & before="Either there are no nonzero elements, or they're all = ", &
          & advance="yes" )
      else
        call output ( "There are no nonzero elements", advance="yes" )
      end if
    end if

  end subroutine Dump_2D_Row_Sparse_Double

  ! ------------------------------------  Dump_2D_Row_Sparse_Real  -----
  subroutine Dump_2D_Row_Sparse_Real ( Array, Name, Format, Width, Exclude )
    real, intent(in) :: Array(:,:)
    character(len=*), intent(in), optional :: Name
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Width ! Items per line, default 5
    real, intent(in), optional :: Exclude  ! Don't print these

    logical :: DidAnyRow, DidOne
    integer :: I ! Row number

    if ( present(name) ) call output ( trim(name), advance='yes' )
    didAnyRow = .false.
    do i = 1, size(array,1)
      call dump_sparse ( array(i,:), format=format, width=width, parens=.true., &
        & rowNum=i, exclude=exclude, didOne = didOne )
      didAnyRow = didAnyRow .or. didOne
    end do
    if ( .not. didAnyRow ) then
      if ( present(exclude) ) then
        call output ( exclude, &
          & before="Either there are no nonzero elements, or they're all = ", &
          & advance="yes" )
      else
        call output ( "There are no nonzero elements", advance="yes" )
      end if
    end if

  end subroutine Dump_2D_Row_Sparse_Real

  ! -----------------------------------------  Dump_2x2xN_Complex  -----
  subroutine Dump_2x2xN_Complex ( Array, Name, Format, Options )
  ! This is for dumping polarized incremental optical depth
    integer, parameter :: RK = kind(0.0e0)
    complex(rk), intent(in) :: Array(:,:,:) ! Better be 2x2xn
    character(len=*), intent(in), optional :: Name
    character(len=*), intent(in), optional :: Format
    character(len=*), intent(in), optional :: Options

    integer :: J
    character(len=64) :: MyFormat

    call theDumpBegins ( options )
    myFormat = sdFormatDefaultCmplx
    if ( present(format) ) myFormat = format
    if ( index(myFormat, '*') > 0 ) &
      & myFormat = numNeedsFormat( one_c4*maxval(abs(array)), format )

    if ( any(shape(array) == 0) ) then
      call empty ( name )
    else
      call name_and_size ( name, dopts(clean)%v, size(array) )
      if ( getOutputStatus( 'start' ) /= 1 ) call newLine
      do j = 1, size( Array,3)
        if (.not. dopts(clean)%v) then
          call output ( j, max(3,ilog10(size(array))+1) , advance='no' )
          call output ( afterSub , advance='no' )
        end if
        call output ( array(1,1,j), myFormat , advance='no' )
        call output ( array(1,2,j), myFormat, advance='no' )
        call finishLine
        call blanks ( max(3,ilog10(size(array))+1) + len(afterSub) )
        call output ( array(2,1,j), myFormat , advance='no' )
        call output ( array(2,2,j), myFormat, advance='no' )
        call finishLine
      end do
    end if
    call theDumpEnds
  end subroutine Dump_2x2xN_Complex

  ! ----------------------------------------  Dump_2x2xN_DComplex  -----
  subroutine Dump_2x2xN_DComplex ( Array, Name, Format, Options )
  ! This is for dumping polarized incremental optical depth
    integer, parameter :: RK = kind(0.0d0)
    complex(rk), intent(in) :: Array(:,:,:) ! Better be 2x2xn
    character(len=*), intent(in), optional :: Name
    character(len=*), intent(in), optional :: Format
    character(len=*), intent(in), optional :: Options

    integer :: J
    character(len=64) :: MyFormat

    call theDumpBegins ( options )
    myFormat = sdFormatDefaultCmplx
    if ( present(format) ) myFormat = format
    if ( index(myFormat, '*') > 0 ) &
      & myFormat = numNeedsFormat( one_c4*maxval(abs(array)), format )

    if ( any(shape(array) == 0) ) then
      call empty ( name )
    else
      call name_and_size ( name, dopts(clean)%v, size(array) )
      if ( getOutputStatus( 'start' ) /= 1 ) call newLine
      do j = 1, size( Array,3)
        if (.not. dopts(clean)%v) then
          call output ( j, max(3,ilog10(size(array))+1) , advance='no' )
          call output ( afterSub , advance='no' )
        end if
        call output ( array(1,1,j), myFormat , advance='no' )
        call output ( array(1,2,j), myFormat, advance='no' )
        call finishLine
        call blanks ( max(3,ilog10(size(array))+1) + len(afterSub) )
        call output ( array(2,1,j), myFormat , advance='no' )
        call output ( array(2,2,j), myFormat, advance='no' )
        call finishLine
      end do
    end if
    call theDumpEnds
  end subroutine Dump_2x2xN_DComplex

  ! -----------------------------------------------  Dump_3D_Char  -----
  subroutine Dump_3D_Char ( Array, Name, FillValue, Width, &
    & Maxlon, Options, TheShape )
    character(len=*), intent(in) :: Array(:,:,:)
    character(len=*), intent(in), optional :: Name
    character(len=*), intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    integer, intent(in), optional :: Maxlon
    character(len=*), intent(in), optional :: Options
    character(len=*), intent(in), optional :: TheShape

    integer :: LON
    integer :: I, J, K, L
    integer :: NumZeroRows
    character(len=len(array)) :: MyFillValue
    integer :: MyWidth

    call theDumpBegins ( options )
    myFillValue = ' '
    if ( present(FillValue) ) myFillValue = FillValue

    lon = len(array(1,1,1))
    if ( dopts(trimIt)%v ) lon = maxval(len_trim(array))
    if ( present(maxlon) ) lon = min(lon, maxlon)
    MyWidth = DefaultWidth ! 10
    if ( present(width) ) MyWidth = width

    numZeroRows = 0
    if ( any(shape(array) == 0) ) then
      call empty ( name )
    else if ( size( Array,1) == 1 ) then
      call dump ( array(1,:,:), name, fillValue=fillValue, &
        & maxlon=maxlon, options=options, width=width )
    else if ( size( Array,2) == 1 ) then
      call dump ( array(:,1,:), name, fillValue=fillValue, &
        & maxlon=maxlon, options=options, width=width )
    else if ( size( Array,3) == 1 ) then
      call dump ( array(:,:,1), name, fillValue=fillValue, &
        & maxlon=maxlon, options=options, width=width )
    else
      call name_and_size ( name, dopts(clean)%v, size(array) )
      if ( getOutputStatus( 'start' ) /= 1 ) call newLine
      do i = 1, size( Array,1)
        do j = 1, size( Array,2)
          do k = 1, size( Array,3), MyWidth
            DumpTheseZeros = dopts(clean)%v .or. &
              & any(array(i,j,k:min(k+2*MyWidth-1, size( Array,3))) /= myFillValue) &
              & .or. &
              & ( k+MyWidth >= size( Array,3) .and. &
              & any(array(i,min(j+1, size( Array,2)),1:min(1+MyWidth-1, size( Array,3))) &
              & /= myFillValue) &
              & )
            if (.not. dopts(clean)%v) then
              if ( DumpTheseZeros ) then
                call say_fill ( (/ i, size( Array,1), j-1, size( Array,2), &
                  & k, size( Array,3) /), numZeroRows, myFillValue, inc=3 )
              else
                numZeroRows = numZeroRows + 1
              end if
            end if
          if ( DumpTheseZeros ) then
              do l = k, min(k+MyWidth-1, size( Array,3))
                  call output ( ' ' // array(i,j,l)(1:lon) // ' ' , advance='no' )
              end do
              call newLine
              numZeroRows = 0
            end if
          end do
        end do
      end do
      call say_fill ( (/ i-1, size( Array,1), j-1, size( Array,2), &
        & k-MyWidth, size( Array,3) /), numZeroRows, myFillValue )
    end if
    call theDumpEnds
  end subroutine Dump_3D_Char

  ! --------------------------------------------  Dump_3D_Complex  -----
  subroutine Dump_3D_Complex ( Array, Name, &
    & FillValue, Width, Format, Lbound, Options, TheShape )
    integer, parameter :: RK = kind(0.0e0)
    complex(rk), intent(in) :: Array(:,:,:)
    character(len=*), intent(in), optional :: Name
    real, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Lbound
    character(len=*), intent(in), optional :: options
    character(len=*), intent(in), optional :: TheShape

    integer :: Base, I, ISlice, J, K, L
    integer :: NumZeroRows
    complex(rk) :: myFillValue
    character(len=64) :: MyFormat
    integer :: MyWidth

    ! Executable
    myFormat = sdFormatDefaultCmplx
    call theDumpBegins ( options )
    myFillValue = 0.
    if ( present(FillValue) ) myFillValue=FillValue

    if ( present(format) ) myFormat = format
    if ( index(myFormat, '*') > 0 ) &
      & myFormat = numNeedsFormat( one_c4*maxval(abs(array)), format )
    
    MyWidth = 3
    if ( present(width) ) MyWidth = width

    include 'dump3db.f9h'
  end subroutine Dump_3D_Complex

  ! -------------------------------------------  Dump_3D_DComplex  -----
  subroutine Dump_3D_DComplex ( Array, Name, &
    & FillValue, Width, Format, Lbound, Options, TheShape )
    integer, parameter :: RK = kind(0.0d0)
    complex(rk), intent(in) :: Array(:,:,:)
    character(len=*), intent(in), optional :: Name
    real(rk), intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Lbound
    character(len=*), intent(in), optional :: options
    character(len=*), intent(in), optional :: TheShape

    integer :: Base, I, ISlice, J, K, L
    integer :: NumZeroRows
    complex(rk) :: myFillValue
    character(len=64) :: MyFormat
    integer :: myWidth

    ! Executable
    myFormat = sdFormatDefaultCmplx
    call theDumpBegins ( options )
    myFillValue = 0.
    if ( present(FillValue) ) myFillValue=FillValue

    if ( present(format) ) myFormat = format
    if ( index(myFormat, '*') > 0 ) &
      & myFormat = numNeedsFormat( one_c4*maxval(abs(array)), format )
    
    MyWidth = 3
    if ( present(width) ) MyWidth = width

    include 'dump3db.f9h'
  end subroutine Dump_3D_DComplex

  ! ---------------------------------------------  Dump_3D_Double  -----
  subroutine Dump_3D_Double ( Array, Name, &
    & FillValue, Width, Format, Lbound, Options, TheShape )
    double precision, intent(in) :: Array(:,:,:)
    character(len=*), intent(in), optional :: Name
    double precision, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Lbound
    character(len=*), intent(in), optional :: options
    character(len=*), intent(in), optional :: TheShape

    integer :: Base, I, ISlice, J, K, L
    integer :: NumZeroRows
    double precision :: myFillValue
    character(len=64) :: myFormat
    integer :: nUnique
    integer :: MyWidth
    integer, dimension(MaxNumElements) :: counts
    double precision, dimension(MaxNumElements) :: elements
    myFormat = sdFormatDefault
    myFillValue = 0.
    if ( present(FillValue) ) myFillValue=FillValue
    include 'dump3d.f9h'
    include 'dump3db.f9h'
  end subroutine Dump_3D_Double

  ! --------------------------------------------  Dump_3D_Integer  -----
  subroutine Dump_3D_Integer ( Array, Name, &
    & FillValue, Format, Width, Lbound, Options, TheShape )
    integer, intent(in) :: Array(:,:,:)
    character(len=*), intent(in), optional :: Name
    integer, intent(in), optional :: FillValue
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Width ! How many numbers per line (10)?
    integer, intent(in), optional :: Lbound ! Low bound for Array
    character(len=*), intent(in), optional :: options
    character(len=*), intent(in), optional :: TheShape

    integer :: Base, I, ISlice, J, K, L
    integer :: NumZeroRows
    integer, dimension(3) :: which
    integer :: how_many

    character(len=64) :: MyFormat
    integer :: myFillValue
    integer :: nUnique
    integer :: MyWidth
    integer, dimension(MaxNumElements) :: counts
    integer, dimension(MaxNumElements) :: elements
    myFormat = 'places=' // INTPLACES ! To sneak places arg into call to output
    myFillValue = 0.
    if ( present(FillValue) ) myFillValue=FillValue
    include 'dump3d.f9h'
    include 'dump3db.f9h'
  end subroutine Dump_3D_Integer

  ! ---------------------------------------------  Dump_3D_Real  -----
  subroutine Dump_3D_Real ( Array, Name, &
    & FillValue, Width, Format, Lbound, Options, TheShape )
    real, intent(in) :: Array(:,:,:)
    character(len=*), intent(in), optional :: Name
    real, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Lbound
    character(len=*), intent(in), optional :: options
    character(len=*), intent(in), optional :: TheShape

    integer :: Base, I, ISlice, J, K, L
    integer :: NumZeroRows
    real    :: myFillValue
    character(len=64) :: MyFormat
    integer :: myWidth
    integer :: nUnique
    integer, dimension(MaxNumElements) :: counts
    real, dimension(MaxNumElements) :: elements
    myFormat = sdFormatDefault
    myFillValue = 0.
    if ( present(FillValue) ) myFillValue=FillValue
    include 'dump3d.f9h'
    include 'dump3db.f9h'
  end subroutine Dump_3D_Real

  ! ---------------------------------------------  Dump_4D_Double  -----
  subroutine Dump_4D_Double ( Array, Name, &
    & FillValue, Width, Format, Lbound, Options )
    double precision, intent(in) :: Array(:,:,:,:)
    character(len=*), intent(in), optional :: Name
    double precision, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Lbound
    character(len=*), intent(in), optional :: Options

    integer :: ISlice, J
    character(len=64) :: MyFormat
    integer :: nUnique
    integer, dimension(MaxNumElements) :: Counts
    double precision, dimension(MaxNumElements) :: Elements
    myFormat = sdFormatDefault
    include 'dump4d.f9h'
    include 'dump4db.f9h'
  end subroutine Dump_4D_Double

  ! ---------------------------------------------  Dump_4D_Real  -----
  subroutine Dump_4D_Real ( Array, Name, &
    & FillValue, Width, Format, Lbound, Options )
    real, intent(in) :: Array(:,:,:,:)
    character(len=*), intent(in), optional :: Name
    real, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Lbound
    character(len=*), intent(in), optional :: Options

    integer :: ISlice, J
    character(len=64) :: MyFormat
    integer :: nUnique
    integer, dimension(MaxNumElements) :: Counts
    real, dimension(MaxNumElements) :: Elements
    myFormat = sdFormatDefault
    include 'dump4d.f9h'
    include 'dump4db.f9h'
  end subroutine Dump_4D_Real

  ! ------------------------------------------------------  Empty  -----
  subroutine Empty ( Name )
    character(len=*), intent(in), optional :: Name

    if ( present(name) .and. .not. NameHasBeenPrinted ) then
      call output ( name , advance='no' )
      call output ( ' is ' , advance='no' )
      nameHasBeenPrinted = .true.
    end if
    call output ( 'empty', advance='yes' )

  end subroutine Empty

  ! -------------------------------------------------  FinishLine  -----
  ! Print a newLine, optionally echoing the item name
  subroutine FinishLine
    ! Executable
    if ( len_trim(nameOnEachLine) < 1 ) then
    ! No op
    elseif( PrintNameAtLineEnd ) then
      call blanksToColumn ( 80-len_trim(nameOnEachLine) )
      call output ( trim(nameOnEachLine), advance='no' )
    else
      call blanks (2)
      call output ( trim(nameOnEachLine), advance='no' )
    end if
    call newLine ( dont_make_blank_line=.true. )
  end subroutine FinishLine

  ! -----------------------------------------------------  ILOG10  -----
  integer function ILOG10(int)
    integer, intent(in) :: int
    ilog10=nint(log10(real(int)))
  end function ILOG10

  ! ----------------------------------------------  Name_And_Size  -----
  subroutine Name_And_Size ( Name, Clean, Size, TheShape )
    character(len=*), intent(in), optional :: Name
    logical, intent(in) :: Clean
    integer, intent(in) :: Size
    character(len=*), intent(in), optional :: TheShape

    if ( present(name) .and. .not. dopts(laconic)%v ) then
      if ( len_trim(name) < 1 ) return
      if ( .not. nameHasBeenPrinted ) then
        if ( .not. clean .and. .not. present(theShape) ) then
          call PrintName ( Name, NameHasBeenPrinted )
          return
        else
          call output ( name , advance='no' )
        endif
        if ( present(theShape) ) call output ( theShape , advance='no' )
      end if
      if ( clean ) then 
        call output ( trim(" \ ") ) ! This goofiness is to outwit an incorrect
                                    ! Intel compiler.
        call output ( size , advance='no' )
      end if
      if ( size == 1 ) call output ( ' ' , advance='no' )
      nameHasBeenPrinted = .true.
    end if

  end subroutine Name_And_Size

  ! =====     Private Procedures     ===================================

  ! -------------------------------------------------  Dump_Bogus  -----
  ! Never used--just here to tell Makefiles that dumpstats.f9h is
  ! a prerequisite for Dump_0 because perl script f90makedep.pl
  ! won't follow .f9h files to look for more uses and includes.
  ! When will we repair the perl script?
  !
  subroutine Dump_Bogus ( Array, Name, &
    & FillValue, Format, Width, Lbound, Options )
    integer, intent(in) :: Array(:)
    character(len=*), intent(in), optional :: Name
    integer, intent(in), optional :: FillValue
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Width ! How many numbers per line (10)?
    integer, intent(in), optional :: Lbound ! Low bound for Array
    character(len=*), optional, intent(in) :: options
    integer :: myFillValue
    myFillValue = 0
    myRank = 0
    include 'dumpstats.f9h'
  end subroutine Dump_Bogus

  ! -----------------------------------------  DumpCollapsedArray  -----
  ! This family of subroutines dumps a lower-rank representation of
  ! an array
  subroutine DumpCollapsedArray_1D_Integer ( Array, Name, FillValue, Options )
    INTEGER, intent(in) :: Array(:)
    character(len=*), intent(in) :: Name
    INTEGER, intent(in), optional :: FillValue
    character(len=*), intent(in), optional :: options
    call output( 'Did not expect to dump collapsed 1d array', advance='yes' )
  end subroutine DumpCollapsedArray_1D_Integer

  subroutine DumpCollapsedArray_1D_Double ( Array, Name, FillValue, Options )
    double precision, intent(in) :: Array(:)
    character(len=*), intent(in) :: Name
    double precision, intent(in), optional :: FillValue
    character(len=*), intent(in), optional :: options
    call output( 'Did not expect to dump collapsed 1d array', advance='yes' )
  end subroutine DumpCollapsedArray_1D_Double

  subroutine DumpCollapsedArray_1D_Real ( Array, Name, FillValue, Options )
    real, intent(in) :: Array(:)
    character(len=*), intent(in) :: Name
    real, intent(in), optional :: FillValue
    character(len=*), intent(in), optional :: options
    call output( 'Did not expect to dump collapsed 1d array', advance='yes' )
  end subroutine DumpCollapsedArray_1D_Real

  subroutine DumpCollapsedArray_2D_Double ( Array, Name, FillValue, Options )
    double precision, intent(in) :: Array(:,:)
    character(len=*), intent(in) :: Name
    double precision, intent(in), optional :: FillValue
    character(len=*), intent(in), optional :: options
    ! For dumping lower-rank collapsed representations
    double precision, dimension(size( Array, 1)) :: nums
    logical, dimension(size( Array, 1))          :: logs

    ! dump numerical representation
    if ( index(collapseOptions, 'num') > 0 ) then
      call collapse( Array, nums, options=collapseOptions )
      call dump( nums, name, fillvalue, Options=options )
    end if
    if ( index(collapseOptions, 'any') > 0 .or. &
      &  index(collapseOptions, 'all') > 0 ) then
      call collapse( Array, logs=logs, options=collapseOptions )
      call dump( logs, name, options=options )
    end if
  end subroutine DumpCollapsedArray_2D_Double

  subroutine DumpCollapsedArray_2D_Real ( Array, Name, FillValue, Options )
    real, intent(in) :: Array(:,:)
    character(len=*), intent(in) :: Name
    real, intent(in), optional :: FillValue
    character(len=*), intent(in), optional :: options
    ! For dumping lower-rank collapsed representations
    real, dimension(size( Array, 1)) :: nums
    logical, dimension(size( Array, 1))          :: logs

    ! dump numerical representation
    if ( index(collapseOptions, 'num') > 0 ) then
      call collapse( Array, nums, options=collapseOptions )
      call dump( nums, name, fillvalue, Options=options )
    end if
    if ( index(collapseOptions, 'any') > 0 .or. &
      &  index(collapseOptions, 'all') > 0 ) then
      call collapse( Array, logs=logs, options=collapseOptions )
      call dump( logs, name, options=options )
    end if
  end subroutine DumpCollapsedArray_2D_Real

  subroutine DumpCollapsedArray_2D_Integer ( Array, Name, FillValue, Options )
    INTEGER, intent(in) :: Array(:,:)
    character(len=*), intent(in) :: Name
    INTEGER, intent(in), optional :: FillValue
    character(len=*), intent(in), optional :: options
    ! For dumping lower-rank collapsed representations
    INTEGER, dimension(size( Array, 1)) :: nums
    logical, dimension(size( Array, 1))          :: logs

    ! dump numerical representation
    if ( index(collapseOptions, 'num') > 0 ) then
      call collapse( Array, nums, options=collapseOptions )
      call dump( nums, name, fillvalue, Options=options )
    end if
    if ( index(collapseOptions, 'any') > 0 .or. &
      &  index(collapseOptions, 'all') > 0 ) then
      call collapse( Array, logs=logs, options=collapseOptions )
      call dump( logs, name, options=options )
    end if
  end subroutine DumpCollapsedArray_2D_Integer

  subroutine DumpCollapsedArray_3D_Double ( Array, Name, FillValue, Options )
    double precision, intent(in) :: Array(:,:,:)
    character(len=*), intent(in) :: Name
    double precision, intent(in), optional :: FillValue
    character(len=*), intent(in), optional :: options
    ! For dumping lower-rank collapsed representations
    double precision, dimension(size( Array, 1),size( Array, 2)) :: nums
    logical, dimension(size( Array, 1),size( Array, 2))          :: logs

    ! dump numerical representation
    if ( index(collapseOptions, 'num') > 0 ) then
      call collapse( Array, nums, options=collapseOptions )
      call dump( nums, name, fillvalue, Options=options )
    end if
    if ( index(collapseOptions, 'any') > 0 .or. &
      &  index(collapseOptions, 'all') > 0 ) then
      call collapse( Array, logs=logs, options=collapseOptions )
      call dump( logs, name, options=options )
    end if
  end subroutine DumpCollapsedArray_3D_Double

  subroutine DumpCollapsedArray_3D_Real ( Array, Name, FillValue, Options )
    real, intent(in) :: Array(:,:,:)
    character(len=*), intent(in) :: Name
    real, intent(in), optional :: FillValue
    character(len=*), intent(in), optional :: options
    ! For dumping lower-rank collapsed representations
    real, dimension(size( Array, 1),size( Array, 2)) :: nums
    logical, dimension(size( Array, 1),size( Array, 2))          :: logs

    ! dump numerical representation
    if ( index(collapseOptions, 'num') > 0 ) then
      call collapse( Array, nums, options=collapseOptions )
      call dump( nums, name, fillvalue, Options=options )
    end if
    if ( index(collapseOptions, 'any') > 0 .or. &
      &  index(collapseOptions, 'all') > 0 ) then
      call collapse( Array, logs=logs, options=collapseOptions )
      call dump( logs, name, options=options )
    end if
  end subroutine DumpCollapsedArray_3D_Real

  subroutine DumpCollapsedArray_3D_Integer ( Array, Name, FillValue, Options )
    INTEGER, intent(in) :: Array(:,:,:)
    character(len=*), intent(in) :: Name
    INTEGER, intent(in), optional :: FillValue
    character(len=*), intent(in), optional :: options
    ! For dumping lower-rank collapsed representations
    INTEGER, dimension(size( Array, 1),size( Array, 2)) :: nums
    logical, dimension(size( Array, 1),size( Array, 2))          :: logs

    ! dump numerical representation
    if ( index(collapseOptions, 'num') > 0 ) then
      call collapse( Array, nums, options=collapseOptions )
      call dump( nums, name, fillvalue, Options=options )
    end if
    if ( index(collapseOptions, 'any') > 0 .or. &
      &  index(collapseOptions, 'all') > 0 ) then
      call collapse( Array, logs=logs, options=collapseOptions )
      call dump( logs, name, options=options )
    end if
  end subroutine DumpCollapsedArray_3D_Integer

  ! -----------------------------------------  ArrayShapeToString  -----
  function ArrayShapeToString ( ArrayShape ) result ( String )
    ! Given an array of integers return the shape as a string
    ! E.g., given (/4,2,6/) return '4*2*6'
    integer, dimension(:), intent(in) :: ArrayShape
    character(len=16) :: String
    ! Internal variables
    integer :: i
    ! Executable
    string = ' '
    if ( size(arrayShape) < 1 ) return
    do i=1, size( arrayshape )
      string = catLists( string, arrayShape(i), inseparator='*' )
    end do
    
  end function ArrayShapeToString
  
  ! ----------------------------------------------------  PrintIt  -----
  ! This family of subroutines exists only so that we can generically call
  ! output with either a numeric arg or a character string, 
  ! trimming if the latter
  subroutine PrintIt_Char ( it, format )
    character(len=*) :: it
    character(len=*), intent(in), optional :: Format
    call output ( trim(it), advance='no' )
  end subroutine PrintIt_Char

  subroutine PrintIt_int ( it, format )
    integer :: it
    character(len=*), intent(in), optional :: Format
    call output ( it, format=format, advance='no' )
  end subroutine PrintIt_int

  subroutine PrintIt_Real ( it, format )
    real :: it
    character(len=*), intent(in), optional :: Format
    call output ( it, format=format, advance='no' )
  end subroutine PrintIt_Real

  subroutine PrintIt_Double ( it, format )
    double precision :: it
    character(len=*), intent(in), optional :: Format
    call output ( it, format=format, advance='no' )
  end subroutine PrintIt_Double

  subroutine PrintIt_Complex ( it, format )
    integer, parameter :: RK = kind(0.0e0)
    complex(rk) :: it
    character(len=*), intent(in), optional :: Format
    call output ( it, format=format, advance='no' )
  end subroutine PrintIt_Complex

  subroutine PrintIt_DComplex ( it, format )
    integer, parameter :: RK = kind(0.0d0)
    complex(rk) :: it
    character(len=*), intent(in), optional :: Format
    call output ( it, format=format, advance='no' )
  end subroutine PrintIt_DComplex

  ! -------------------------------------------  PrintPercentages  -----
  ! Prints a nicely-formatted summary of equal, unequal, etc.
  ! using output
  subroutine PrintPercentages ( name, equal, unequal )
    character(len=*), intent(in), optional :: Name
    integer, intent(in) :: equal
    integer, intent(in) :: unequal
    if ( equal+unequal < 1 ) return
    myPCTFormat  = DefaultPCTFormat
    if ( PCTFormat /= '*' ) myPCTFormat = PCTFormat
    if ( present(name) ) call output ( trim(name), advance='no' )
    call blanks( 1, advance='no' )
    call output( fillvaluerelation, advance='no' )
    call output ( ', !', advance='no' )
    call output( fillvaluerelation, advance='no' )
    call output ( ' (%) ', advance='no' )
    call output ( equal, advance='no' )
    call output ( ': ', advance='no' )
    call output ( unequal, advance='no' )
    if ( .not. statsOnOneLine ) then
      call newline
      call blanks(10)
    end if
    call output ( '( ', advance='no' )
    call output ( 100*equal/(equal+unequal+0.), format = myPCTFormat, advance='no' )
    call output ( ': ', advance='no' )
    call output ( 100*unequal/(equal+unequal+0.), format = myPCTFormat, advance='no' )
    call output ( ' )', advance='no' )
    call finishLine
  end subroutine PrintPercentages

  ! ------------------------------------------------  PrintRMSetc  -----
  ! This family of routines prints a nicely-formatted list of min, max, etc.
  ! using output
  subroutine PrintRMSetc_Double ( Name, min, max, rms, mean )
    character(len=*), intent(in), optional :: Name
    double precision, intent(in) :: min
    double precision, intent(in) :: max
    double precision, intent(in) :: rms
    double precision, intent(in), optional :: mean
    !
    character(len=16) :: originalSDFormat
    !
    originalSDFormat = outputOptions%sdFormatDefault
    outputOptions%sdFormatDefault = rmsFormat
    include 'printRMSetc.f9h'
    outputOptions%sdFormatDefault = originalSDFormat
  end subroutine PrintRMSetc_Double

  subroutine PrintRMSetc_Real ( Name, min, max, rms, mean )
    character(len=*), intent(in), optional :: Name
    real, intent(in) :: min
    real, intent(in) :: max
    real, intent(in) :: rms
    real, intent(in), optional :: mean
    character(len=16) :: originalSDFormat
    !
    originalSDFormat = outputOptions%sdFormatDefault
    outputOptions%sdFormatDefault = rmsFormat
    include 'printRMSetc.f9h'
    outputOptions%sdFormatDefault = originalSDFormat
  end subroutine PrintRMSetc_Real

  subroutine PrintRMSetc_int ( Name, in_min, in_max, rms, mean )
    character(len=*), intent(in), optional :: Name
    integer, intent(in)        :: in_min
    integer, intent(in)        :: in_max
    real, intent(in)           :: rms
    real, intent(in), optional :: mean
    ! Internal variables
    real                       :: min
    real                       :: max
    ! Executable
    min = in_min
    max = in_max
    include 'printRMSetc.f9h'
  end subroutine PrintRMSetc_int

  ! This family of subroutines print subscripts to the left
  ! of each dumped row, sometimes noting that repeated lines
  ! that have been omitted for brevity
  ! ----------------------------------------------  Say_Fill_Char  -----
  subroutine Say_Fill_Char ( Subs, NumZeroRows, Fill, Inc, Format  )
    character(len=*), intent(in) :: Fill
    include 'Say_Fill.f9h'
  end subroutine Say_Fill_Char

  ! -------------------------------------------  Say_Fill_Complex  -----
  subroutine Say_Fill_Complex ( Subs, NumZeroRows, Fill, Inc, Format )
    integer, parameter :: RK = kind(0.0e0)
    complex(rk), intent(in) :: Fill
    include 'Say_Fill.f9h'
  end subroutine Say_Fill_Complex

  ! ------------------------------------------  Say_Fill_DComplex  -----
  subroutine Say_Fill_DComplex ( Subs, NumZeroRows, Fill, Inc, Format )
    integer, parameter :: RK = kind(0.0d0)
    complex(rk), intent(in) :: Fill
    include 'Say_Fill.f9h'
  end subroutine Say_Fill_DComplex

  ! --------------------------------------------  Say_Fill_Double  -----
  subroutine Say_Fill_Double ( Subs, NumZeroRows, Fill, Inc, Format )
    double precision, intent(in) :: Fill
    include 'Say_Fill.f9h'
  end subroutine Say_Fill_Double

  ! -----------------------------------------------  Say_Fill_Int  -----
  subroutine Say_Fill_Int ( Subs, NumZeroRows, Fill, Inc, Format )
    integer, intent(in) :: Fill
    include 'Say_Fill.f9h'
  end subroutine Say_Fill_Int

  ! ----------------------------------------------  Say_Fill_Real  -----
  subroutine Say_Fill_Real ( Subs, NumZeroRows, Fill, Inc, Format )
    real, intent(in) :: Fill
    include 'Say_Fill.f9h'
  end subroutine Say_Fill_Real

  ! ----------------------------------------------------  Say_Subs -----
  subroutine Say_Subs ( Subs, NumZeroRows )
    integer, intent(in) :: Subs(:)
    integer, intent(in) :: NumZeroRows
    call say_subs_only ( subs )
    call output ( ' ' , advance='no' )
    call output ( numZeroRows , advance='no' )
    call output ( ' lines of ', advance='no' )
  end subroutine Say_Subs

  ! -----------------------------------------------  Say_Subs_Only -----
  subroutine Say_Subs_Only ( Subs )
    integer, intent(in) :: Subs(:)
    integer :: I
    do i = 1, size(subs), 2
      call output ( subs(i), places=max(4,ilog10(subs(i+1))+2) , advance='no' )
    end do
    call output ( afterSub , advance='no' )
  end subroutine Say_Subs_Only
  
  ! ----------------------------------------------  ShowColumnNums -----
  subroutine ShowColumnNums( lineLength, skip )
    ! show column numbers
    ! usu. printed below other data to help guide the eye as to where
    ! stuff actually is
    ! E.g.,
    !   TTTT  TT  T
    !       FF  FF FFFF
    !   12345678901234567890          <-- These are the 
    !            1         2          <-- two lines we print here
    ! Args:
    integer, intent(in)           :: lineLength ! How many per line
    integer, intent(in), optional :: skip       ! start numbering after skip
    ! Internal variables
    integer                          :: bloc
    character(len=10), dimension(2)  :: line
    integer                          :: multiple
    integer                          :: mySkip
    integer                          :: numBlocs
    ! Executable
    numBlocs = 1 + (lineLength-1)/10
    mySkip = 0
    if ( present(skip) ) mySkip = skip
    line(1) = '1234567890'
!   line(2) = '         1'
    do multiple=1, 2
      if ( mySkip > 0 ) call blanks(mySkip)
      do bloc = 1, numBlocs
        write(line(2), '(i10)' ) bloc
        call output( line(multiple), advance='no' )
      end do
      call newline
    end do
  end subroutine ShowColumnNums

  ! --------------------------------------------------  SnipOption -----
  function SnipOption ( options, particular ) result ( snipped )
    ! Snip from options a particular option commanding dump or diff
    ! Args:
    character(len=*), optional, intent(in) :: options
    character(len=1), optional, intent(in) :: particular
    character(len=16)                      :: snipped
    ! Executable
    snipped = ' '
    if ( .not. present(options) .or. .not. present(particular) ) return
    snipped = delete( options, particular )
  end function SnipOption

  ! --------------------------------------------------  UniqueOnly -----
  logical function UniqueOnly ( options )
    character(len=*), intent(in), optional :: options
    ! Executable
    uniqueonly = .true.
    if ( .not. present(options) ) return
    uniqueonly = uniqueonly .and. &
      & .not. any( indexes(options, (/'w','s','r'/)) > 0 )

  end function UniqueOnly

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Dump_0

! $Log$
! Revision 2.163  2021/09/09 22:58:02  pwagner
! Added subroutine AsExplicitValues
!
! Revision 2.162  2021/04/01 23:45:18  pwagner
! DefaultWidth allplied consistently to all character-values
!
! Revision 2.161  2019/10/25 20:54:58  pwagner
! Dump full row of zeros if narrower than my width
!
! Revision 2.160  2019/07/22 22:21:38  pwagner
! -N opt now will display the indices where true and the indices where false if array is logical-valued
!
! Revision 2.159  2019/04/24 19:17:01  vsnyder
! Add Unit argument to several dumps
!
! Revision 2.158  2018/10/27 01:37:20  vsnyder
! Spiff sparse dumps
!
! Revision 2.157  2018/08/21 01:52:14  vsnyder
! Add Exclude argument to sparse dumps
!
! Revision 2.156  2018/08/16 02:17:43  vsnyder
! Add row sparse dump
!
! Revision 2.155  2018/04/25 01:47:43  vsnyder
! Add 1D sparse real and double precision dumps
!
! Revision 2.154  2018/04/13 00:22:49  pwagner
! Improved comments; explain use of appearance flags like 'MyName\h'
!
! Revision 2.153  2018/02/28 19:51:35  pwagner
! Moved PrintName to dump_options
!
! Revision 2.152  2017/12/07 02:40:10  vsnyder
! Remove some unreferenced use named.  Delete some unreferenced variable
! declarations.  Don't use host-associated DO index variables; make them
! local.
!
! Revision 2.151  2017/11/30 20:46:19  pwagner
! Avoid leaving OldNameOnEachLine undefined
!
! Revision 2.150  2017/11/03 20:01:40  pwagner
! Most array gymnastics moved from MLSFillValues to HyperSlabs module
!
! Revision 2.149  2017/10/11 20:56:23  pwagner
! Fixed bug in Dump_2D_Logical
!
! Revision 2.148  2017/09/15 22:37:55  pwagner
! Dont repeat name if empty
!
! Revision 2.147  2017/09/14 18:31:14  pwagner
! Take care to skip a space when dumping chars; aso not to print an unwanted blank line
!
! Revision 2.146  2017/09/07 23:44:30  pwagner
! Added PrintNameAsHeadline and PrintNameInBanner options to dump
!
! Revision 2.145  2017/09/07 20:59:35  pwagner
! Dont print 2 spaces before name
!
! Revision 2.144  2017/08/03 20:37:31  pwagner
! Try harder to avoid messing up dumps called with 'c(lean)'
!
! Revision 2.143  2017/07/31 23:05:47  pwagner
! Moved TheDumpEnds to dump_options; try to avoid blank lines
!
! Revision 2.142  2017/07/31 22:18:22  vsnyder
! Option to print FALSE as dot to make it easier to see TRUE
!
! Revision 2.141  2017/07/19 22:42:53  pwagner
! Added PrintName; may PrintNameAtLineEnd; PrintRMSetc now public
!
! Revision 2.140  2016/10/20 23:05:41  pwagner
! Separated row, column indexes in 2d arrays
!
! Revision 2.139  2016/07/28 03:28:17  vsnyder
! Moved comments to Dump_1 and Dump_Options
!
! Revision 2.138  2016/07/28 01:42:27  vsnyder
! Refactoring dump and diff
!
! Revision 2.137  2016/04/05 23:54:57  pwagner
! -v verbose option added; usu. will print name on each line
!
! Revision 2.136  2016/03/31 22:59:03  pwagner
! Added dumpTextfile
!
! Revision 2.135  2016/03/23 00:22:17  pwagner
! Diff now able to print name on each line; repaired error in printRMSEtc when array is integer
!
! Revision 2.134  2016/01/12 00:46:51  pwagner
! May override DEFAULTWidth when dumping char array
!
! Revision 2.133  2015/08/25 18:38:27  vsnyder
! Include Lbound in 'all values are the same' dumps
!
! Revision 2.132  2015/01/29 01:23:29  vsnyder
! Make sure MyFillValue has a value before references
!
! Revision 2.131  2014/08/06 23:02:21  vsnyder
! Use kind C_Int16_t from ISO_C_Bindings instead of INTEGER*2
!
! Revision 2.130  2014/01/09 00:24:29  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.129  2013/09/26 15:25:41  pwagner
! Added 2-byte integer dumps
!
! Revision 2.128  2013/08/12 23:47:25  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.127  2013/06/28 18:08:38  pwagner
! Note if all logical elements equal
!
! Revision 2.126  2013/06/19 23:14:39  pwagner
! Remove more unused stuff; give default value to myDirect
!
! Revision 2.125  2013/01/10 00:18:43  pwagner
! Dumps with -N just show where NaNs are
!
! Revision 2.124  2012/09/11 21:10:24  pwagner
! Requires 'N' option to show where NaNs, Infs are located
!
! Revision 2.123  2012/07/05 23:47:31  pwagner
! Added restoreDumpConfig
!
! Revision 2.122  2012/06/22 20:25:52  pwagner
! Specify advance arg because we may now set default to 'yes'
!
! Revision 2.121  2012/06/13 23:59:37  pwagner
! dumpDumpOptions optionally dumps available diff/dump options
!
! Revision 2.120  2012/01/09 22:25:55  pwagner
! Distinguish 'r' option to print rms of ratios and 'R' option for rms of values
!
! Revision 2.119  2011/12/15 01:47:45  pwagner
! Accepts W[i] option; deletes more prolix notices of rank reduction
!
! Revision 2.118  2011/12/07 01:14:44  pwagner
! Added option to show bandwidth of banded arrays
!
! Revision 2.117  2011/11/11 00:30:22  vsnyder
! Simplify notice of rank reduction
!
! Revision 2.116  2011/07/26 20:40:24  pwagner
! Added 4d diffs, too
!
! Revision 2.115  2011/07/23 00:16:25  vsnyder
! Use (1.,0.) instead of CMPLX(1.,0.) to avoid Intel gripe
!
! Revision 2.114  2011/07/15 23:23:43  pwagner
! Can now dump 4d arrays, with certain restrictions
!
! Revision 2.113  2011/07/12 00:15:01  pwagner
! Improved dumps; format option now more flexible; complex arrays parallel real ones
!
! Revision 2.112  2011/06/23 17:27:41  pwagner
! Added function to difference args with option to supply fillvalue or period
!
! Revision 2.111  2011/04/26 20:54:12  pwagner
! Can now dump dates
!
! Revision 2.110  2011/04/20 22:27:09  pwagner
! Fixed long-standing bug in dumping 2d char array
!
! Revision 2.109  2011/04/18 19:12:54  vsnyder
! Turn on nameHasBeenPrinted in more places
!
! Revision 2.108  2011/04/04 23:08:27  pwagner
! May diff 2d integer arrays
!
! Revision 2.107  2011/03/08 00:04:40  pwagner
! Skip printing row of zeros only if multiple
!
! Revision 2.106  2011/02/25 18:50:26  pwagner
! Added optional width arg to char array dumps
!
! Revision 2.105  2011/02/05 01:01:41  pwagner
! Added new dopt_ dump options; transpose and collapse work better
!
! Revision 2.104  2011/01/20 01:16:01  pwagner
! New '-l' option dumps collapsed representation of higher ranked array
!
! Revision 2.103  2011/01/04 00:48:26  pwagner
! DEFAULTMaxlon can now set max width of 1d char array dumps
!
! Revision 2.102  2010/10/14 18:44:01  pwagner
! Can now dump lists
!
! Revision 2.101  2010/02/04 23:05:39  vsnyder
! Remove USE and declaration for unreferenced names.
! Declare a kind parameter in unfiltereddiff*
!
! Revision 2.100  2010/01/29 21:08:21  pwagner
! gave myFillValue a value to stop lf95 complaints
!
! Revision 2.99  2009/11/20 01:15:11  pwagner
! Decided against using allFinite
!
! Revision 2.98  2009/11/20 01:12:50  pwagner
! Added new option H or 'shape' to just show array rank, Shape
!
! Revision 2.97  2009/10/30 23:02:41  pwagner
! Should not double-print name if only whole array diff
!
! Revision 2.96  2009/10/26 18:53:33  pwagner
! Fixed bug only NAG found
!
! Revision 2.95  2009/10/19 17:33:26  pwagner
! Trying to prevent double-printing of name
!
! Revision 2.94  2009/10/13 00:09:04  pwagner
! Percentages printed with better format; dumpstats.f9h now a direct prerequisite
!
! Revision 2.93  2009/09/10 20:58:00  pwagner
! 3 ways to summarize diffs: 'b' (table), 'r' (rms), 's' (number of differences)
!
! Revision 2.92  2009/08/18 20:41:17  pwagner
! making INTPLACES public can change appearance of dumped ints
!
! Revision 2.91  2009/08/17 16:55:41  pwagner
! Among options string 'd' means 'direct'
!
! Revision 2.90  2009/06/26 00:15:18  pwagner
! Added dumpDumpOptions
!
! Revision 2.89  2009/06/24 22:36:32  pwagner
! Make use of dump1-3.f9h include files
!
! Revision 2.88  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.87  2009/06/16 17:12:55  pwagner
! Changed api for dump, diff routines; now rely on options for most optional behavior
!
! Revision 2.86  2009/05/14 21:59:04  pwagner
! New optional arg AS_GAPS dor dumping 1d bit logicals
!
! Revision 2.85  2009/05/08 00:40:09  pwagner
! Avoid printing blank name when dumped with empty name string
!
! Revision 2.84  2009/04/01 23:30:49  pwagner
! Improved appearance when printing 'nnn lines of xxx not printed.'
!
! Revision 2.83  2008/11/24 19:34:47  pwagner
! Less wasteful of memory; should not segment dault so often
!
! Revision 2.82  2008/10/24 23:21:49  pwagner
! Limits 3d diffs to prevent segment faults
!
! Revision 2.81  2008/08/27 16:23:41  pwagner
! Added dumpSums
!
! Revision 2.80  2008/07/10 00:13:30  pwagner
! SelfDiff can now print half wave lengths
!
! Revision 2.79  2008/07/09 16:30:21  pwagner
! selfdiff can take optional arg unique
!
! Revision 2.78  2008/06/18 20:56:18  pwagner
! New optional arg 'unique' dumps print unique elements, counts only
!
! Revision 2.77  2008/01/09 20:53:22  pwagner
! When NaNs short-circuit dump or diff, print how many and where
!
! Revision 2.76  2008/01/07 21:37:57  pwagner
! Replace DEFAULTUNDEFINEDVALUE with user-settable undefinedValue
!
! Revision 2.75  2007/10/18 23:37:41  pwagner
! Added dumpTable
!
! Revision 2.74  2007/10/12 23:34:49  pwagner
! Using new howfar procedure to summarize diffs after peeling away outliers
!
! Revision 2.73  2007/09/20 18:39:37  pwagner
! Dont interrupt tables of dumped numbers with stamping
!
! Revision 2.72  2007/09/13 21:09:57  pwagner
! -rms and -s combined add new info about & how near
!
! Revision 2.71  2007/07/17 00:21:58  pwagner
! diff stats when rms concentrate on showing ratios
!
! Revision 2.70  2007/07/11 22:29:23  vsnyder
! Add more dumps for complex, add 2x2xN dumps
!
! Revision 2.69  2007/06/14 18:38:36  pwagner
! Allow separate rmsFormat to be set at class level
!
! Revision 2.68  2007/04/14 00:32:47  vsnyder
! Remove OPTIONAL attribute from Name1 and Name args of DIFF_...
!
! Revision 2.67  2007/04/03 02:49:23  vsnyder
! Don't look at non-present dummy argument
!
! Revision 2.66  2007/03/23 00:14:30  pwagner
! new optional maxlon arg to dumping n-d chars; sdFormatDefault for numeric dumps public and changeable
!
! Revision 2.65  2007/03/07 21:01:45  pwagner
! Some small changes unrelated to real bugs elsewhere
!
! Revision 2.64  2007/01/31 00:05:43  pwagner
! Added interface for dumping bit arrays
!
! Revision 2.63  2006/08/23 20:06:25  pwagner
! Restored more backward-compatible arglist to Dump_STRLIST
!
! Revision 2.62  2006/08/22 20:40:04  pwagner
! Added required arg=separator to Dump_STRLIST
!
! Revision 2.61  2006/07/11 00:24:04  pwagner
! use fillValue properly when computing rms etc.
!
! Revision 2.60  2006/06/24 23:07:04  pwagner
! Changes to reduce memory footprint computing statistics
!
! Revision 2.59  2006/06/09 18:50:12  pwagner
! Avoid dumping an entire array if all elements the same
!
! Revision 2.58  2006/05/24 20:38:14  pwagner
! Allow any of 3 ordering relations for dumping pct
!
! Revision 2.57  2006/04/20 01:09:30  vsnyder
! Don't print lines of zeroes in complex dumps
!
! Revision 2.56  2006/03/22 23:48:52  vsnyder
! Print 'lines ... not printed' instead of 'rows...' to avoid confusion with matrices
!
! Revision 2.55  2006/03/15 17:34:28  pwagner
! Fixed bug causing incorrect rms when diffing with fill values
!
! Revision 2.54  2006/03/03 23:04:55  pwagner
! May dump logical-valued hashes
!
! Revision 2.53  2006/02/28 21:42:29  pwagner
! Replace fillvalues with 0 before computing rms (should actually remove)
!
! Revision 2.52  2006/01/27 01:01:37  pwagner
! Can now dump hashes
!
! Revision 2.51  2005/12/17 00:58:56  pwagner
! dumps of rms, etc. should appear uniform
!
! Revision 2.50  2005/12/16 23:25:13  pwagner
! dumpSize moved from dump0 to output_m
!
! Revision 2.49  2005/12/16 00:04:06  pwagner
! Changes to reflect new MLSFillValues module
!
! Revision 2.48  2005/11/04 18:49:02  pwagner
! Added SelfDiff procedures to dump diffs among consecutive 1d array elems
!
! Revision 2.47  2005/10/03 18:05:52  pwagner
! Allocated memory now dumped in units of MEMORY_UNITS
!
! Revision 2.46  2005/07/20 01:33:47  vsnyder
! Simplify DumpSize routines
!
! Revision 2.45  2005/06/22 17:25:48  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.44  2005/05/12 20:43:54  pwagner
! Added diff routines
!
! Revision 2.43  2004/12/14 21:31:59  pwagner
! Optional trim arg added to multi-dim char dumps
!
! Revision 2.42  2004/08/16 17:09:14  pwagner
! 3d integer dumps interface redone like others
!
! Revision 2.41  2004/08/04 23:19:01  pwagner
! Much moved from MLSStrings to MLSStringLists
!
! Revision 2.40  2004/07/23 19:47:20  vsnyder
! Add Lbound to Dump_1d_[double,integer,real]
!
! Revision 2.39  2004/07/23 18:34:59  vsnyder
! Add Lbound argument to Dump_1d_Logical
!
! Revision 2.38  2004/07/21 19:57:21  pwagner
! New rms, wholeArray options to some dumps
!
! Revision 2.37  2004/06/11 19:04:29  pwagner
! which_ints_are_it renamed findAll, moved MLSSets.f90
!
! Revision 2.36  2004/05/27 23:25:33  pwagner
! Added stats parameter; also more 1-d get fillvalue
!
! Revision 2.35  2004/04/05 17:47:43  livesey
! Split dumpsize into real and integer versions
!
! Revision 2.34  2004/04/03 05:43:23  livesey
! Added DumpSize
!
! Revision 2.33  2004/03/30 00:44:10  vsnyder
! Remove unused variable declaration
!
! Revision 2.32  2004/02/26 21:53:31  pwagner
! Can dump ,-separated string list
!
! Revision 2.31  2004/01/21 22:02:19  vsnyder
! Don't use number of entities per line as default width for counter
!
! Revision 2.30  2003/09/19 02:00:14  vsnyder
! More about the goofy Intel compiler
!
! Revision 2.29  2003/09/15 17:43:41  livesey
! Cosmetic change for fussy (and wrong) intel compiler
!
! Revision 2.28  2003/09/06 00:48:40  vsnyder
! Specify default formats with a module parameter instead of literals.
! Change default (1x,1pg13.6) to (1pg14.6) to avoid problems with length
! calculation in output_m.
!
! Revision 2.27  2003/08/08 20:45:42  vsnyder
! Made say_fill_* generic, made them test for numZeroRows, and made them
! optionally do say_subs_only.  This simplified several dump routines.
! Added optional Format arguments in several more routines.
!
! Revision 2.26  2003/07/04 02:41:33  vsnyder
! Substantial simplification by putting little things into subroutines
!
! Revision 2.25  2003/07/02 01:07:27  vsnyder
! Add complex output
!
! Revision 2.24  2003/05/21 19:20:40  vsnyder
! Start a new line after \ 1 if size==1 and clean
!
! Revision 2.23  2003/05/06 00:15:03  pwagner
! Fixed incompatibility with FilterShapes
!
! Revision 2.20.2.4  2003/05/05 23:00:05  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.20.2.3  2003/04/18 20:26:05  vsnyder
! Add Width and Format arguments to 1D_Real and 1D_Double
!
! Revision 2.20.2.2  2003/03/27 23:18:33  vsnyder
! Put new-lines in better places
!
! Revision 2.20.2.1  2003/03/14 00:25:47  vsnyder
! Add Dump_2D_Logical, cosmetic changes
!
! Revision 2.20  2002/12/02 23:34:14  pwagner
! Now can dump name/value pairs
!
! Revision 2.19  2002/10/08 00:09:08  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.18  2002/09/13 18:08:12  pwagner
! May change matrix precision rm from r8
!
! Revision 2.17  2002/02/14 23:21:18  vsnyder
! Work on dumping masks
!
! Revision 2.16  2001/12/08 00:47:51  pwagner
! Added Dump_1d_Real for s.p. arrays
!
! Revision 2.15  2001/11/29 23:50:53  pwagner
! Added optional blase arg to Dump_nd_Char; fixed bug where optional
! format not passed from Dump_3d_int
!
! Revision 2.14  2001/11/28 23:32:01  livesey
! Fixed bug where Dump_2d_Integer didn't pass format to 1d dump.
!
! Revision 2.13  2001/10/25 23:30:39  pwagner
! Improved Dump_nd_Double to skip rows (e.g., of zeros)
!
! Revision 2.12  2001/10/24 18:11:14  pwagner
! which_ints_are_it now works properly
!
! Revision 2.11  2001/10/23 22:40:37  pwagner
! Now dumps 1d,2d,3d char arrays and 3d ints
!
! Revision 2.10  2001/09/28 22:43:20  vsnyder
! Don't print rows of zeroes
!
! Revision 2.9  2001/09/11 22:52:32  livesey
! Added printing of sizes
!
! Revision 2.8  2001/05/11 22:44:54  vsnyder
! Print transpose of 2d-double if it would take fewer lines.  Get rid of
! double printing of "without mask"
!
! Revision 2.7  2001/05/08 20:27:24  vsnyder
! Added an optional 'format' argument in a few more places
!
! Revision 2.6  2001/05/08 17:21:02  livesey
! Added a `clean' option to the array dumps.  This omits the indices at
! the start, making it easier for other programs to read output.
!
! Revision 2.5  2001/05/03 02:12:34  vsnyder
! Insert copyright notice, clean up CVS stuff, cosmetics
!
! Revision 2.4  2001/03/10 03:39:58  vsnyder
! Improve handling of "name" if size==1 or size==0
!
! Revision 2.3  2001/03/02 01:32:08  livesey
! Handles larger arrays better
!
! Revision 2.2  2001/02/28 21:35:27  livesey
! Added dump logical 1d
!
! Revision 2.1  2000/09/13 20:38:50  vsnyder
! Initial code
!
