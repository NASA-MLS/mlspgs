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
module MLSStats1                 ! Calculate statistics of rank n arrays
!=============================================================================
  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use HighOutput, only: OutputNamedValue
  use MLSCommon, only: Interval_T
  use MLSKinds, only: R4, R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use MLSFinds, only: FindAll, FindFirst, FindLast
  use MLSStringLists, only: Catlists
  use MLSStrings, only: Lowercase
  use Output_M, only: Blanks, Newline, Output
  use Sort_M, only: Sort, Sortp

  implicit none
! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (data types and parameters)
! Stat_t                          Basic user-defined data type containing
!                                  all standard statistics
! FillValueRelation               Whether to use '=' (default) or '<', '>'

!     (subroutines and functions)
! AllStats                        Computes some or all standard statistics
! ConfidInterval                  Computes the confidence interval
! Dump                            Prints a STAT_T with all standard statistics
! SetUp                           Initialize and set up arrays for the Stat_T datatype
! MLSMin                          Finds min of an array
!                                   (excluding FillValues or negative precisions)
! MLSMax, MLSMean, MLSMedian,
!    MLSStddev, MLSRms            Similar to MLSMin for other statistics
! Hownear                         Finds how near 2 arrays are in %
! Howfar                          Inverse of hownear--given a % finds stats of 
!                                   diffs of nearest %
!                                 If mode='median' finds stats of array1
!                                   differences from its median value
! Pdf                             Finds the pdf for a sample x given a Stat_T
!                                   or with x1, x2 its integral over [x1, x2]
! Ratios                          Statistics of ratio between 2 arrays;
!                                    can be used to track changes to a reference
!                                    goldbrick standard
! Reset                           Resets or Destroys a Stat_T
! Statfunction                    Given a set of values it returns a STAT_T
! Statistics                      Similar to STATFUNCTION, but may accumulate
!                                   statistics in STAT_T over multiple calls
! === (end of toc) ===

! === (start of api) ===
!     (user-defined types)
! Stat_T 
!             (
!             int   count, 
!             int   fillcount, 
!             r8    max   ,
!             r8    mean  ,
!             r8    median,
!             r8    min   ,
!             r8    stddev,
!             r8    rms   ,
!             int   indexing(3),
!             int   nbins, 
!             r8    bounds(2)  ,
!             *int  bincount(:)
!             ) 

!     (subroutines and functions)
! allstats( values, &
!      & [int nbins], [*bounds(2)], [log addedData], [fillValue], [precision], &
!      & [int count], [int fillcount], &
!      & [min], [max], [mean], [stddev], [rms], [median], [int indexing(3)], &
!      & [int bincount], [log doDump] )
!      values may go up to rank 4, typed either r4 or r8
! confidInterval( values, pct(:), intervals(:), &
!      & [char* mode], [type(values) estimate] )
! dump( stat_t statistic, [char which] )
! type(values) mls$fun( values, [fillvalue] )
! howfar( array1, array2, pct(:), stat_t gaps, char* mode, &
!      & [array1AtN(3,:)],  [array2AtN(3,:)] )
! hownear( array1, array2, pct, [type(array1) gaps(:)], [type(array1) gapratios(:)] )
! r8 pdf( r8 x, stat_t statistic )
! r8 pdf( r8 x1, r8 x2, stat_t statistic )
! ratios( array1, array2, type(array1) exvalues, type(array1) exratios, &
!      & type(array1) minratio, type(array1) maxratio, type(array1) meanratio, &
!      & type(array1) stddevratio, type(array1) rmsratio, &
!      & type(array1) medianratio, &
!      & [fillValue], [char op] )
! reset( stat_t statistic )
! SetUp( stat_t statistic, [int NBins], [r8 Bounds(2)], [log verbose] )
! stat_t statFunction( r8 values(:), [r8 fillValue], [r8 precision(:) )
! statistics( r8 values(:), stat_t statistic, [r8 fillValue], [r8 precision(:) )
!
! Notes:
! (1) Unless specified explicitly, 'values', 'precision', 'pct',,
! 'array1', 'array2', etc. may be any numerically typed array with
! 0 < rank < 4 except complex
! (2) 'fillvalue', 'min', 'max', 'mean', 'stddev', 'rms', and 'median'
! are scalars with the same numerical type as 'values'
! (3) The $fun in mls$fun can be any of the names in (2) beginning with 'min'
! === (end of api) ===
  private
  
  public :: Stat_t             ! the data type
  public :: Allstats, ConfidInterval, Dump, Howfar, Hownear, Ratios, Statistics
  public :: MLScount
  public :: MLSMin, MLSMax, MLSMean, MLSMedian, MLSStddev, MLSRms
  public :: Pdf
  public :: Reset, SetUp
  public :: Statfunction
  
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! The functions and subroutines are based on
  ! (1) The MATH77 library subroutine; and
  ! (2) optionally filtering out fill values from array entries
  ! (3) Optionally, instead of filtering fill values, an array of precisions,
  ! with the same shape and numeric type as the original values,
  ! can be supplied. In this case, values for which the corresponding
  ! precisions < 0 will be filtered out.
  ! If both fillValue and precision are supplied, 
  ! the precision array will be ignored.
  
  ! Interfaces have been supplied so that most procedures accept
  ! arrays of rank 1, 2, or 3 and either single or double precision
  
  ! Note that multidimensional arrays are reshaped into longer
  ! 1-d arrays, rather than returning reduced rank results
  
  ! In addition to min, max, mean, stddev functions of MATH77
  ! median and rms functions have been created.
  ! Missing so far is a mode function, chi^2, 
  ! confidence level, etc.
  
  ! If called via a subroutine,
  ! an optional returned array is bincount. If supplied, nbin and bounds must
  ! also be supplied and satisfy the conditions: 
  ! nbins > 2, X2 > X1, size(bincount) >= nbins
  ! where bounds = (X1, X2)
  ! bincount will be suitable for plotting a histogram of the data. 
  ! Let h = (X2-X1) / (nbins-2)
  ! Then the bins are the intervals between the points in the following set:
  ! {-Inf  X1  X1+h  X1+2h  ..  X2-2h  X2-h  X2  Inf}
  
  ! Also, if called via a subroutine, successive calls may "merge" data
  ! allowing the calculation of statistics for a larger cumulative data set
  ! Except, in this case, medians are not found correctly
  ! (hard to fix this bug; perhaps why median not part of original MATH77)

  ! This is the main datatype, a statistic type that both holds
  ! and can progressively accumulate info
  
  ! However, accumulating data is inconsistent with obtaining correct
  ! values for the median and for the indexes of the max, mean, and medians
  ! stored in the indexing array
  
  ! (Would making fillValue a component makes sense?)

  type Stat_T
    integer  :: count = 0     ! If > 0, merging data from prior call(s)
    integer  :: fillcount = 0 ! Number of times fillValues ignored
    real(r8) :: max
    real(r8) :: mean
    real(r8) :: median        ! incorrect if merging data
    real(r8) :: min
    real(r8) :: stddev
    real(r8) :: rms
    integer, dimension(3) :: indexing = 0  !index of  (/ max, median, min /)
    ! The next 3 deal with histogramming data
    integer  :: nbins = 0  !NCELLS; put > 2 if histogramming
    real(r8), dimension(2) :: bounds = 1.d0 ! X1, X2; put X2 > X1
    integer, dimension(:), pointer :: bincount => null() ! IHIST
  end type Stat_T
  
  type(Stat_T), save :: MLSStat
  
  ! consider whether one of {"=" (default), "<", ">"}
  ! when calculating %ages (and possibly rms, etc.)
  character(len=1), public, save :: fillValueRelation = '='
  logical, public, save          :: showIndexing = .false. ! show where max, min
  logical, public, save          :: statsOnOneLine = .false.

  interface allstats
    module procedure allstats_d1r4, allstats_d2r4, allstats_d3r4, allstats_d4r4
    module procedure allstats_d1r8, allstats_d2r8, allstats_d3r8, allstats_d4r8
  end interface
  
  interface dump
    module procedure dump_all
    module procedure dump_selected
  end interface
  
  interface dump_if_selected
    module procedure dump_if_selected_int_array
    module procedure dump_if_selected_int_scalar
    module procedure dump_if_selected_r8_array
    module procedure dump_if_selected_r8_scalar
  end interface
  
  interface filterValues
    module procedure filterValues_r4
    module procedure filterValues_r8
  end interface
  
  interface ConfidInterval
    module procedure ConfidInterval_int, ConfidInterval_r4, ConfidInterval_r8
  end interface

  interface howfar
    module procedure howfar_d1int, howfar_d2int, howfar_d3int
    module procedure howfar_d1r4, howfar_d2r4, howfar_d3r4, howfar_d4r4
    module procedure howfar_d1r8, howfar_d2r8, howfar_d3r8, howfar_d4r8
  end interface
  
  interface hownear
    module procedure hownear_d1int, hownear_d2int, hownear_d3int
    module procedure hownear_d1r4, hownear_d2r4, hownear_d3r4, hownear_d4r4
    module procedure hownear_d1r8, hownear_d2r8, hownear_d3r8, hownear_d4r8
  end interface
  
  interface MLSMin
    module procedure MLSMin_d1int, MLSMin_d2int, MLSMin_d3int
    module procedure MLSMin_d1r4, MLSMin_d2r4, MLSMin_d3r4
    module procedure MLSMin_d1r8, MLSMin_d2r8, MLSMin_d3r8
  end interface
  
  interface MLSMax
    module procedure MLSMax_d1int, MLSMax_d2int, MLSMax_d3int
    module procedure MLSMax_d1r4, MLSMax_d2r4, MLSMax_d3r4
    module procedure MLSMax_d1r8, MLSMax_d2r8, MLSMax_d3r8
  end interface
  
  interface MLSMean
    module procedure MLSMean_d1int, MLSMean_d2int, MLSMean_d3int
    module procedure MLSMean_d1r4, MLSMean_d2r4, MLSMean_d3r4
    module procedure MLSMean_d1r8, MLSMean_d2r8, MLSMean_d3r8
  end interface
  
  interface MLSMedian
    module procedure MLSMedian_d1int, MLSMedian_d2int, MLSMedian_d3int
    module procedure MLSMedian_d1r4, MLSMedian_d2r4, MLSMedian_d3r4
    module procedure MLSMedian_d1r8, MLSMedian_d2r8, MLSMedian_d3r8
  end interface
  
  interface MLSStddev
    module procedure MLSStddev_d1int, MLSStddev_d2int, MLSStddev_d3int
    module procedure MLSStddev_d1r4, MLSStddev_d2r4, MLSStddev_d3r4, MLSStddev_d4r4
    module procedure MLSStddev_d1r8, MLSStddev_d2r8, MLSStddev_d3r8, MLSStddev_d4r8
  end interface
  
  interface MLSRms
    module procedure MLSRms_d1int, MLSRms_d2int, MLSRms_d3int
    module procedure MLSRms_d1r4, MLSRms_d2r4, MLSRms_d3r4
    module procedure MLSRms_d1r8, MLSRms_d2r8, MLSRms_d3r8
  end interface
  
  interface mlscount
    module procedure mlscount_d1int
    module procedure mlscount_d1r4
    module procedure mlscount_d1r8
  end interface
  
  interface pdf
    module procedure pdf_1, pdf_2
  end interface
  
  interface ratios
    module procedure ratios_d1int, ratios_d2int, ratios_d3int
    module procedure ratios_d1r4, ratios_d2r4, ratios_d3r4, ratios_d4r4
    module procedure ratios_d1r8, ratios_d2r8, ratios_d3r8, ratios_d4r8
  end interface
  
  interface shrinkarray
    module procedure shrinkarray_int, shrinkarray_r4, shrinkarray_r8
  end interface
  
  interface stat1
    module procedure STAT1_r4, STAT1_r8
  end interface
  
  real, parameter :: FAC = 64.0E0
  logical, parameter :: DEEBUG = .false.
  ! Function names
  ! Stored using sequence of integers
  integer, parameter              :: FN_COUNT     = 1
  integer, parameter              :: FN_MIN       = FN_COUNT + 1
  integer, parameter              :: FN_MAX       = FN_MIN + 1
  integer, parameter              :: FN_MEAN      = FN_MAX + 1
  integer, parameter              :: FN_STDDEV    = FN_MEAN + 1
  integer, parameter              :: FN_RMS       = FN_STDDEV + 1
  integer, parameter              :: FN_MEDIAN    = FN_RMS + 1
  integer, parameter              :: FN_FILLCOUNT = FN_MEDIAN + 1
  integer, parameter              :: FN_COUNT_TOT = FN_FILLCOUNT + 1

contains
      ! ------------------- allstats_d1r4 -----------------------
      subroutine allstats_d1r4( values, &
        & nbins, bounds, addedData, fillValue, precision, &
        & count, fillcount, min, max, mean, stddev, rms, median, &
        & bincount, indexing, doDump )
        integer, parameter :: rk = r4
        include 'allstats_d1.f9h'
      end subroutine allstats_d1r4

      ! ------------------- allstats_d1r8 -----------------------
      subroutine allstats_d1r8( values, &
        & nbins, bounds, addedData, fillValue, precision, &
        & count, fillcount, min, max, mean, stddev, rms, median, &
        & bincount, indexing, doDump )
        integer, parameter :: rk = r8
        include 'allstats_d1.f9h'
      end subroutine allstats_d1r8

      ! ------------------- allstats_d2r4 -----------------------
      subroutine allstats_d2r4( values, &
        & nbins, bounds, addedData, fillValue, precision, &
        & count, fillcount, min, max, mean, stddev, rms, median, &
        & bincount, indexing, doDump )
        ! Args
        real(r4), dimension(:,:), intent(in)           :: values
        integer, optional, intent(in)                  :: nbins
        real(r4), dimension(2), optional, intent(in)   :: bounds
        logical, optional, intent(in)                  :: addeddata
        real(r4), optional, intent(in)                 :: fillValue
        real(r4), dimension(:,:), optional, intent(in) :: precision
        integer, optional, intent(inout)               :: count
        integer, optional, intent(inout)               :: fillcount
        real(r4), optional, intent(out)                :: min
        real(r4), optional, intent(out)                :: max
        real(r4), optional, intent(out)                :: mean
        real(r4), optional, intent(out)                :: stddev
        real(r4), optional, intent(out)                :: rms
        real(r4), optional, intent(out)                :: median
        integer, dimension(:), optional, intent(out)   :: bincount
        integer, dimension(3), optional, intent(out)   :: indexing
        logical , optional, intent(in)                 :: doDump
        ! Internal variables
        integer, dimension(2)                          :: shp
        ! Executable
        shp = shape(values)
        if ( .not. present(precision) ) then
          call allstats_d1r4(reshape(values, (/shp(1)*shp(2)/)), &
            & nbins, bounds, addedData, fillValue, &
            & count=count, fillcount=fillcount, min=min, max=max, mean=mean, &
            & stddev=stddev, rms=rms, median=median, bincount=bincount, &
            & indexing=indexing, doDump=doDump)
        else
          call allstats_d1r4(reshape(values, (/shp(1)*shp(2)/)), &
            & nbins, bounds, addedData, fillValue, &
            & reshape(precision, (/shp(1)*shp(2)/)),&
            & count=count, fillcount=fillcount, min=min, max=max, mean=mean, &
            & stddev=stddev, rms=rms, median=median, bincount=bincount, &
            & indexing=indexing, doDump=doDump)
        endif
      end subroutine allstats_d2r4

      ! ------------------- allstats_d2r8 -----------------------
      subroutine allstats_d2r8( values, &
        & nbins, bounds, addedData, fillValue, precision, &
        & count, fillcount, min, max, mean, stddev, rms, median, &
        & bincount, indexing, doDump )
        ! Args
        real(r8), dimension(:,:), intent(in)           :: values
        integer, optional, intent(in)                  :: nbins
        real(r8), dimension(2), optional, intent(in)   :: bounds
        logical, optional, intent(in)                  :: addeddata
        real(r8), optional, intent(in)                 :: fillValue
        real(r8), dimension(:,:), optional, intent(in) :: precision
        integer, optional, intent(inout)               :: count
        integer, optional, intent(inout)               :: fillcount
        real(r8), optional, intent(out)                :: min
        real(r8), optional, intent(out)                :: max
        real(r8), optional, intent(out)                :: mean
        real(r8), optional, intent(out)                :: stddev
        real(r8), optional, intent(out)                :: rms
        real(r8), optional, intent(out)                :: median
        integer, dimension(:), optional, intent(out)   :: bincount
        integer, dimension(3), optional, intent(out)   :: indexing
        logical , optional, intent(in)                 :: doDump
        ! Internal variables
        integer, dimension(2)                          :: shp
        ! Executable
        shp =shape(values)
        if ( .not. present(precision) ) then
          call allstats_d1r8(reshape(values, (/shp(1)*shp(2)/)), &
            & nbins, bounds, addedData, fillValue, &
            & count=count, fillcount=fillcount, min=min, max=max, mean=mean, &
            & stddev=stddev, rms=rms, median=median, bincount=bincount, &
            & indexing=indexing, doDump=doDump)
        else
          call allstats_d1r8(reshape(values, (/shp(1)*shp(2)/)), &
            & nbins, bounds, addedData, fillValue, &
            & reshape(precision, (/shp(1)*shp(2)/)),&
            & count=count, fillcount=fillcount, min=min, max=max, mean=mean, &
            & stddev=stddev, rms=rms, median=median, bincount=bincount, &
            & indexing=indexing, doDump=doDump)
        endif
      end subroutine allstats_d2r8

      ! ------------------- allstats_d3r4 -----------------------
      subroutine allstats_d3r4( values, &
        & nbins, bounds, addedData, fillValue, precision, &
        & count, fillcount, min, max, mean, stddev, rms, median, &
        & bincount, indexing, doDump )
        ! Args
        real(r4), dimension(:,:,:), intent(in)         :: values
        integer, optional, intent(in)                  :: nbins
        real(r4), dimension(2), optional, intent(in)   :: bounds
        logical, optional, intent(in)                  :: addeddata
        real(r4), optional, intent(in)                 :: fillValue
        real(r4), dimension(:,:,:), optional, intent(in)   :: precision
        integer, optional, intent(inout)               :: count
        integer, optional, intent(inout)               :: fillcount
        real(r4), optional, intent(out)                :: min
        real(r4), optional, intent(out)                :: max
        real(r4), optional, intent(out)                :: mean
        real(r4), optional, intent(out)                :: stddev
        real(r4), optional, intent(out)                :: rms
        real(r4), optional, intent(out)                :: median
        integer, dimension(:), optional, intent(out)   :: bincount
        integer, dimension(3), optional, intent(out)   :: indexing
        logical , optional, intent(in)                 :: doDump
        ! Internal variables
        integer, dimension(3)                          :: shp
        ! Executable
        shp =shape(values)
        if ( .not. present(precision) ) then
          call allstats_d1r4(reshape(values, (/shp(1)*shp(2)*shp(3)/)), &
            & nbins, bounds, addedData, fillValue, &
            & count=count, fillcount=fillcount, min=min, max=max, mean=mean, &
            & stddev=stddev, rms=rms, median=median, bincount=bincount, &
            & indexing=indexing, doDump=doDump)
        else
          call allstats_d1r4(reshape(values, (/shp(1)*shp(2)*shp(3)/)), &
            & nbins, bounds, addedData, fillValue, &
            & reshape(precision, (/shp(1)*shp(2)/)),&
            & count=count, fillcount=fillcount, min=min, max=max, mean=mean, &
            & stddev=stddev, rms=rms, median=median, bincount=bincount, &
            & indexing=indexing, doDump=doDump)
        endif
      end subroutine allstats_d3r4

      ! ------------------- allstats_d3r8 -----------------------
      subroutine allstats_d3r8( values, &
        & nbins, bounds, addedData, fillValue, precision, &
        & count, fillcount, min, max, mean, stddev, rms, median, &
        & bincount, indexing, doDump )
        ! Args
        real(r8), dimension(:,:,:), intent(in)         :: values
        integer, optional, intent(in)                  :: nbins
        real(r8), dimension(2), optional, intent(in)   :: bounds
        logical, optional, intent(in)                  :: addeddata
        real(r8), optional, intent(in)                 :: fillValue
        real(r8), dimension(:,:,:), optional, intent(in)   :: precision
        integer, optional, intent(inout)               :: count
        integer, optional, intent(inout)               :: fillcount
        real(r8), optional, intent(out)                :: min
        real(r8), optional, intent(out)                :: max
        real(r8), optional, intent(out)                :: mean
        real(r8), optional, intent(out)                :: stddev
        real(r8), optional, intent(out)                :: rms
        real(r8), optional, intent(out)                :: median
        integer, dimension(:), optional, intent(out)   :: bincount
        integer, dimension(3), optional, intent(out)   :: indexing
        logical , optional, intent(in)                 :: doDump
        ! Internal variables
        integer, dimension(3)                          :: shp
        ! Executable
        shp =shape(values)
        if ( .not. present(precision) ) then
          call allstats_d1r8(reshape(values, (/shp(1)*shp(2)*shp(3)/)), &
            & nbins, bounds, addedData, fillValue, &
            & count=count, fillcount=fillcount, min=min, max=max, mean=mean, &
            & stddev=stddev, rms=rms, median=median, bincount=bincount, &
            & indexing=indexing, doDump=doDump)
        else
          call allstats_d1r8(reshape(values, (/shp(1)*shp(2)*shp(3)/)), &
            & nbins, bounds, addedData, fillValue, &
            & reshape(precision, (/shp(1)*shp(2)/)),&
            & count=count, fillcount=fillcount, min=min, max=max, mean=mean, &
            & stddev=stddev, rms=rms, median=median, bincount=bincount, &
            & indexing=indexing, doDump=doDump)
        endif
      end subroutine allstats_d3r8

      ! ------------------- allstats_d4r4 -----------------------
      subroutine allstats_d4r4( values, &
        & nbins, bounds, addedData, fillValue, precision, &
        & count, fillcount, min, max, mean, stddev, rms, median, &
        & bincount, indexing, doDump )
        ! Args
        real(r4), dimension(:,:,:,:), intent(in)       :: values
        integer, optional, intent(in)                  :: nbins
        real(r4), dimension(2), optional, intent(in)   :: bounds
        logical, optional, intent(in)                  :: addeddata
        real(r4), optional, intent(in)                 :: fillValue
        real(r4), dimension(:,:,:,:), optional, intent(in)   :: precision
        integer, optional, intent(inout)               :: count
        integer, optional, intent(inout)               :: fillcount
        real(r4), optional, intent(out)                :: min
        real(r4), optional, intent(out)                :: max
        real(r4), optional, intent(out)                :: mean
        real(r4), optional, intent(out)                :: stddev
        real(r4), optional, intent(out)                :: rms
        real(r4), optional, intent(out)                :: median
        integer, dimension(:), optional, intent(out)   :: bincount
        integer, dimension(3), optional, intent(out)   :: indexing
        logical , optional, intent(in)                 :: doDump
        ! Internal variables
        integer, dimension(4)                          :: shp
        ! Executable
        shp =shape(values)
        if ( .not. present(precision) ) then
          call allstats_d1r4(reshape(values, (/shp(1)*shp(2)*shp(3)*shp(4)/)), &
            & nbins, bounds, addedData, fillValue, &
            & count=count, fillcount=fillcount, min=min, max=max, mean=mean, &
            & stddev=stddev, rms=rms, median=median, bincount=bincount, &
            & indexing=indexing, doDump=doDump)
        else
          call allstats_d1r4(reshape(values, (/shp(1)*shp(2)*shp(3)*shp(4)/)), &
            & nbins, bounds, addedData, fillValue, &
            & reshape(precision, (/shp(1)*shp(2)/)),&
            & count=count, fillcount=fillcount, min=min, max=max, mean=mean, &
            & stddev=stddev, rms=rms, median=median, bincount=bincount, &
            & indexing=indexing, doDump=doDump)
        endif
      end subroutine allstats_d4r4

      ! ------------------- allstats_d4r8 -----------------------
      subroutine allstats_d4r8( values, &
        & nbins, bounds, addedData, fillValue, precision, &
        & count, fillcount, min, max, mean, stddev, rms, median, &
        & bincount, indexing, doDump )
        ! Args
        real(r8), dimension(:,:,:,:), intent(in)       :: values
        integer, optional, intent(in)                  :: nbins
        real(r8), dimension(2), optional, intent(in)   :: bounds
        logical, optional, intent(in)                  :: addeddata
        real(r8), optional, intent(in)                 :: fillValue
        real(r8), dimension(:,:,:,:), optional, intent(in)   :: precision
        integer, optional, intent(inout)               :: count
        integer, optional, intent(inout)               :: fillcount
        real(r8), optional, intent(out)                :: min
        real(r8), optional, intent(out)                :: max
        real(r8), optional, intent(out)                :: mean
        real(r8), optional, intent(out)                :: stddev
        real(r8), optional, intent(out)                :: rms
        real(r8), optional, intent(out)                :: median
        integer, dimension(:), optional, intent(out)   :: bincount
        integer, dimension(3), optional, intent(out)   :: indexing
        logical , optional, intent(in)                 :: doDump
        ! Internal variables
        integer, dimension(4)                          :: shp
        ! Executable
        shp =shape(values)
        if ( .not. present(precision) ) then
          call allstats_d1r8(reshape(values, (/shp(1)*shp(2)*shp(3)*shp(4)/)), &
            & nbins, bounds, addedData, fillValue, &
            & count=count, fillcount=fillcount, min=min, max=max, mean=mean, &
            & stddev=stddev, rms=rms, median=median, bincount=bincount, &
            & indexing=indexing, doDump=doDump)
        else
          call allstats_d1r8(reshape(values, (/shp(1)*shp(2)*shp(3)*shp(4)/)), &
            & nbins, bounds, addedData, fillValue, &
            & reshape(precision, (/shp(1)*shp(2)/)),&
            & count=count, fillcount=fillcount, min=min, max=max, mean=mean, &
            & stddev=stddev, rms=rms, median=median, bincount=bincount, &
            & indexing=indexing, doDump=doDump)
        endif
      end subroutine allstats_d4r8

      ! ------------------- mlscount -----------------------
      ! This family of functions returns the count of values that are
      ! (depending on method)
      ! 'non-fill'   /= fillValue (default)
      ! 'fill'        = fillValue
      ! 'total'       independent of FillValue (same as size(values))
      ! By intention and by design 'total' = 'non-fill' + 'fill'
      ! If fillValue is absent, then 'non-fill' and 'total' will return
      ! identical values while 'fill' will return 0
      function mlscount_d1int(values, fillValue, method) result(iValue)
        integer, parameter                             :: KINDVALUE = r4
        character(len=*), optional, intent(in)         :: method
        ! Args
        integer, dimension(:), intent(in)      :: values
        integer                                :: iValue
        integer, optional                      :: FillValue
        if ( present(fillValue) ) then
          iValue = mlscount( real(values, KINDVALUE), &
            & fillValue=real(fillValue, KINDVALUE), &
            & method=method )
        else
          iValue = mlscount( real(values, KINDVALUE), &
            & method=method )
        endif
      end function mlscount_d1int

      function mlscount_d1r4(values, fillValue, method) result(iValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_COUNT ! unused
        character(len=*), optional, intent(in)         :: method
        character(len=8)                               :: myMethod
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        integer                                :: iValue
        include 'stats0.f9h'
        myMethod = 'non-fill'
        if ( present(method) ) myMethod = method
        select case( lowercase(myMethod(1:1)) )
        case ('n')
          iValue = count
        case ('f')
          iValue = fillcount
        case ('t')
          iValue = count + fillcount
        case default
          iValue = count
        end select
      end function mlscount_d1r4

      function mlscount_d1r8(values, fillValue, method) result(iValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_COUNT ! unused
        character(len=*), optional, intent(in)         :: method
        character(len=8)                               :: myMethod
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        integer                                :: iValue
        include 'stats0.f9h'
        myMethod = 'non-fill'
        if ( present(method) ) myMethod = method
        select case( lowercase(myMethod(1:1)) )
        case ('n')
          iValue = count
        case ('f')
          iValue = fillcount
        case ('t')
          iValue = count + fillcount
        case default
          iValue = count
        end select
      end function mlscount_d1r8

      ! ------------------- MLSMin_d1int -----------------------
      function MLSMin_d1int(values, fillValue) result(iValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_MIN
        ! Args
        integer, dimension(:), intent(in)      :: values
        integer                                :: iValue
        integer, optional                      :: FillValue
        if ( present(fillValue) ) then
          iValue = MLSMin( real(values, KINDVALUE), real(fillValue, KINDVALUE) )
        else
          iValue = MLSMin( real(values, KINDVALUE) )
        endif
      end function MLSMin_d1int

      ! ------------------- MLSMin_d2int -----------------------
      function MLSMin_d2int(values, fillValue) result(iValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_MIN
        ! Args
        integer, dimension(:,:), intent(in)      :: values
        integer                                :: iValue
        integer, optional                      :: FillValue
        if ( present(fillValue) ) then
          iValue = MLSMin( real(values, KINDVALUE), real(fillValue, KINDVALUE) )
        else
          iValue = MLSMin( real(values, KINDVALUE) )
        endif
      end function MLSMin_d2int

      ! ------------------- MLSMin_d3int -----------------------
      function MLSMin_d3int(values, fillValue) result(iValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_MIN
        ! Args
        integer, dimension(:,:,:), intent(in)      :: values
        integer                                :: iValue
        integer, optional                      :: FillValue
        if ( present(fillValue) ) then
          iValue = MLSMin( real(values, KINDVALUE), real(fillValue, KINDVALUE) )
        else
          iValue = MLSMin( real(values, KINDVALUE) )
        endif
      end function MLSMin_d3int

      ! ------------------- MLSMax_d1int -----------------------
      function MLSMax_d1int(values, fillValue) result(iValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_max
        ! Args
        integer, dimension(:), intent(in)      :: values
        integer                                :: iValue
        integer, optional                      :: FillValue
        if ( present(fillValue) ) then
          iValue = MLSMax( real(values, KINDVALUE), real(fillValue, KINDVALUE) )
        else
          iValue = MLSMax( real(values, KINDVALUE) )
        endif
      end function MLSMax_d1int

      ! ------------------- MLSMax_d2int -----------------------
      function MLSMax_d2int(values, fillValue) result(iValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_max
        ! Args
        integer, dimension(:,:), intent(in)      :: values
        integer                                :: iValue
        integer, optional                      :: FillValue
        if ( present(fillValue) ) then
          iValue = MLSMax( real(values, KINDVALUE), real(fillValue, KINDVALUE) )
        else
          iValue = MLSMax( real(values, KINDVALUE) )
        endif
      end function MLSMax_d2int

      ! ------------------- MLSMax_d3int -----------------------
      function MLSMax_d3int(values, fillValue) result(iValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_max
        ! Args
        integer, dimension(:,:,:), intent(in)      :: values
        integer                                :: iValue
        integer, optional                      :: FillValue
        if ( present(fillValue) ) then
          iValue = MLSMax( real(values, KINDVALUE), real(fillValue, KINDVALUE) )
        else
          iValue = MLSMax( real(values, KINDVALUE) )
        endif
      end function MLSMax_d3int

      ! ------------------- MLSMean_d1int -----------------------
      function MLSMean_d1int(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_mean
        ! Args
        integer, dimension(:), intent(in)      :: values
        real(KINDVALUE)                        :: rValue
        integer, optional                      :: FillValue
        if ( present(fillValue) ) then
          rValue = MLSMean( real(values, KINDVALUE), real(fillValue, KINDVALUE) )
        else
          rValue = MLSMean( real(values, KINDVALUE) )
        endif
      end function MLSMean_d1int

      ! ------------------- MLSMean_d2int -----------------------
      function MLSMean_d2int(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_mean
        ! Args
        integer, dimension(:,:), intent(in)      :: values
        real(KINDVALUE)                        :: rValue
        integer, optional                      :: FillValue
        if ( present(fillValue) ) then
          rValue = MLSMean( real(values, KINDVALUE), real(fillValue, KINDVALUE) )
        else
          rValue = MLSMean( real(values, KINDVALUE) )
        endif
      end function MLSMean_d2int

      ! ------------------- MLSMean_d3int -----------------------
      function MLSMean_d3int(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_mean
        ! Args
        integer, dimension(:,:,:), intent(in)      :: values
        real(KINDVALUE)                        :: rValue
        integer, optional                      :: FillValue
        if ( present(fillValue) ) then
          rValue = MLSMean( real(values, KINDVALUE), real(fillValue, KINDVALUE) )
        else
          rValue = MLSMean( real(values, KINDVALUE) )
        endif
      end function MLSMean_d3int

      ! ------------------- MLSStddev_d1int -----------------------
      function MLSStddev_d1int(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_stddev
        ! Args
        integer, dimension(:), intent(in)      :: values
        real(KINDVALUE)                        :: rValue
        integer, optional                      :: FillValue
        if ( present(fillValue) ) then
          rValue = MLSStddev( real(values, KINDVALUE), real(fillValue, KINDVALUE) )
        else
          rValue = MLSStddev( real(values, KINDVALUE) )
        endif
      end function MLSStddev_d1int

      ! ------------------- MLSStddev_d2int -----------------------
      function MLSStddev_d2int(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_stddev
        ! Args
        integer, dimension(:,:), intent(in)      :: values
        real(KINDVALUE)                        :: rValue
        integer, optional                      :: FillValue
        if ( present(fillValue) ) then
          rValue = MLSStddev( real(values, KINDVALUE), real(fillValue, KINDVALUE) )
        else
          rValue = MLSStddev( real(values, KINDVALUE) )
        endif
      end function MLSStddev_d2int

      ! ------------------- MLSStddev_d3int -----------------------
      function MLSStddev_d3int(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_stddev
        ! Args
        integer, dimension(:,:,:), intent(in)      :: values
        real(KINDVALUE)                        :: rValue
        integer, optional                      :: FillValue
        if ( present(fillValue) ) then
          rValue = MLSStddev( real(values, KINDVALUE), real(fillValue, KINDVALUE) )
        else
          rValue = MLSStddev( real(values, KINDVALUE) )
        endif
      end function MLSStddev_d3int

      ! ------------------- MLSMedian_d1int -----------------------
      function MLSMedian_d1int(values, fillValue) result(iValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_median
        ! Args
        integer, dimension(:), intent(in)      :: values
        integer                                :: iValue
        integer, optional                      :: FillValue
        if ( present(fillValue) ) then
          iValue = MLSMedian( real(values, KINDVALUE), real(fillValue, KINDVALUE) )
        else
          iValue = MLSMedian( real(values, KINDVALUE) )
        endif
      end function MLSMedian_d1int

      ! ------------------- MLSMedian_d2int -----------------------
      function MLSMedian_d2int(values, fillValue) result(iValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_median
        ! Args
        integer, dimension(:,:), intent(in)      :: values
        integer                                :: iValue
        integer, optional                      :: FillValue
        if ( present(fillValue) ) then
          iValue = MLSMedian( real(values, KINDVALUE), real(fillValue, KINDVALUE) )
        else
          iValue = MLSMedian( real(values, KINDVALUE) )
        endif
      end function MLSMedian_d2int

      ! ------------------- MLSMedian_d3int -----------------------
      function MLSMedian_d3int(values, fillValue) result(iValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_median
        ! Args
        integer, dimension(:,:,:), intent(in)      :: values
        integer                                :: iValue
        integer, optional                      :: FillValue
        if ( present(fillValue) ) then
          iValue = MLSMedian( real(values, KINDVALUE), real(fillValue, KINDVALUE) )
        else
          iValue = MLSMedian( real(values, KINDVALUE) )
        endif
      end function MLSMedian_d3int

      ! ------------------- MLSRms_d1int -----------------------
      function MLSRms_d1int(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_rms
        ! Args
        integer, dimension(:), intent(in)      :: values
        real(KINDVALUE)                        :: rValue
        integer, optional                      :: FillValue
        if ( present(fillValue) ) then
          rValue = MLSRms( real(values, KINDVALUE), real(fillValue, KINDVALUE) )
        else
          rValue = MLSRms( real(values, KINDVALUE) )
        endif
      end function MLSRms_d1int

      ! ------------------- MLSRms_d2int -----------------------
      function MLSRms_d2int(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_rms
        ! Args
        integer, dimension(:,:), intent(in)      :: values
        real(KINDVALUE)                        :: rValue
        integer, optional                      :: FillValue
        if ( present(fillValue) ) then
          rValue = MLSRms( real(values, KINDVALUE), real(fillValue, KINDVALUE) )
        else
          rValue = MLSRms( real(values, KINDVALUE) )
        endif
      end function MLSRms_d2int

      ! ------------------- MLSRms_d3int -----------------------
      function MLSRms_d3int(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_rms
        ! Args
        integer, dimension(:,:,:), intent(in)      :: values
        real(KINDVALUE)                        :: rValue
        integer, optional                      :: FillValue
        if ( present(fillValue) ) then
          rValue = MLSRms( real(values, KINDVALUE), real(fillValue, KINDVALUE) )
        else
          rValue = MLSRms( real(values, KINDVALUE) )
        endif
      end function MLSRms_d3int

      ! ------------------- MLSMin_d1r4 -----------------------
      function MLSMin_d1r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_MIN
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function MLSMin_d1r4

      ! ------------------- MLSMin_d1r8 -----------------------
      function MLSMin_d1r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_MIN
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function MLSMin_d1r8

      ! ------------------- MLSMin_d2r4 -----------------------
      function MLSMin_d2r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_MIN
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function MLSMin_d2r4

      ! ------------------- MLSMin_d2r8 -----------------------
      function MLSMin_d2r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_MIN
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function MLSMin_d2r8

      ! ------------------- MLSMin_d3r4 -----------------------
      function MLSMin_d3r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_MIN
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function MLSMin_d3r4

      ! ------------------- MLSMin_d3r8 -----------------------
      function MLSMin_d3r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_MIN
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function MLSMin_d3r8

      ! ------------------- MLSMax_d1r4 -----------------------
      function MLSMax_d1r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_max
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function MLSMax_d1r4

      ! ------------------- MLSMax_d1r8 -----------------------
      function MLSMax_d1r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_max
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function MLSMax_d1r8

      ! ------------------- MLSMax_d2r4 -----------------------
      function MLSMax_d2r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_max
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function MLSMax_d2r4

      ! ------------------- MLSMax_d2r8 -----------------------
      function MLSMax_d2r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_max
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function MLSMax_d2r8

      ! ------------------- MLSMax_d3r4 -----------------------
      function MLSMax_d3r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_max
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function MLSMax_d3r4

      ! ------------------- MLSMax_d3r8 -----------------------
      function MLSMax_d3r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_max
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function MLSMax_d3r8

      ! ------------------- MLSMean_d1r4 -----------------------
      function MLSMean_d1r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_mean
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function MLSMean_d1r4

      ! ------------------- MLSMean_d1r8 -----------------------
      function MLSMean_d1r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_mean
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function MLSMean_d1r8

      ! ------------------- MLSMean_d2r4 -----------------------
      function MLSMean_d2r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_mean
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function MLSMean_d2r4

      ! ------------------- MLSMean_d2r8 -----------------------
      function MLSMean_d2r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_mean
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function MLSMean_d2r8

      ! ------------------- MLSMean_d3r4 -----------------------
      function MLSMean_d3r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_mean
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function MLSMean_d3r4

      ! ------------------- MLSMean_d3r8 -----------------------
      function MLSMean_d3r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_mean
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function MLSMean_d3r8

      ! ------------------- MLSMedian_d1r4 -----------------------
      function MLSMedian_d1r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_median
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function MLSMedian_d1r4

      ! ------------------- MLSMedian_d1r8 -----------------------
      function MLSMedian_d1r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_median
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function MLSMedian_d1r8

      ! ------------------- MLSMedian_d2r4 -----------------------
      function MLSMedian_d2r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_median
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function MLSMedian_d2r4

      ! ------------------- MLSMedian_d2r8 -----------------------
      function MLSMedian_d2r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_median
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function MLSMedian_d2r8

      ! ------------------- MLSMedian_d3r4 -----------------------
      function MLSMedian_d3r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_median
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function MLSMedian_d3r4

      ! ------------------- MLSMedian_d3r8 -----------------------
      function MLSMedian_d3r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_median
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function MLSMedian_d3r8

      ! ------------------- MLSStddev_d1r4 -----------------------
      function MLSStddev_d1r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_stddev
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function MLSStddev_d1r4

      ! ------------------- MLSStddev_d1r8 -----------------------
      function MLSStddev_d1r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_stddev
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function MLSStddev_d1r8

      ! ------------------- MLSStddev_d2r4 -----------------------
      function MLSStddev_d2r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_stddev
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function MLSStddev_d2r4

      ! ------------------- MLSStddev_d2r8 -----------------------
      function MLSStddev_d2r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_stddev
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function MLSStddev_d2r8

      ! ------------------- MLSStddev_d3r4 -----------------------
      function MLSStddev_d3r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_stddev
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function MLSStddev_d3r4

      ! ------------------- MLSStddev_d3r8 -----------------------
      function MLSStddev_d3r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_stddev
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function MLSStddev_d3r8

      ! ------------------- MLSStddev_d4r4 -----------------------
      function MLSStddev_d4r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_stddev
        ! Args
        real(KINDVALUE), dimension(:,:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function MLSStddev_d4r4

      ! ------------------- MLSStddev_d4r8 -----------------------
      function MLSStddev_d4r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_stddev
        ! Args
        real(KINDVALUE), dimension(:,:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function MLSStddev_d4r8

      ! ------------------- MLSRms_d1r4 -----------------------
      function MLSRms_d1r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_rms
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function MLSRms_d1r4

      ! ------------------- MLSRms_d1r8 -----------------------
      function MLSRms_d1r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_rms
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function MLSRms_d1r8

      ! ------------------- MLSRms_d2r4 -----------------------
      function MLSRms_d2r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_rms
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function MLSRms_d2r4

      ! ------------------- MLSRms_d2r8 -----------------------
      function MLSRms_d2r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_rms
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function MLSRms_d2r8

      ! ------------------- MLSRms_d3r4 -----------------------
      function MLSRms_d3r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_rms
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function MLSRms_d3r4

      ! ------------------- MLSRms_d3r8 -----------------------
      function MLSRms_d3r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_rms
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function MLSRms_d3r8
      
      ! ------------------- ConfidInterval -----------------------
      ! This family of routines, given a confidence measured as a percentage,
      ! returns the lower and upper bounds of a conffidence interval
      ! within which the data match an estimate
      ! By default the estimate will be the mean, but you may specify
      ! another estimate
      subroutine ConfidInterval_int( values, &
        & pct, intervals, mode, estimate )
        integer, parameter                             :: KINDVALUE = r4
        integer, dimension(:), intent(in)                   :: values
        real(KINDVALUE), dimension(:), intent(in)           :: pct
        type(interval_t), dimension(:), intent(out)         :: intervals
        character(len=*), intent(in), optional              :: mode
        integer, intent(in), optional                       :: estimate
        ! Internal variables
        ! Executable
        if ( present(estimate) ) then
          call ConfidInterval_r4( real(values, kind=kindvalue), pct, intervals, &
            & mode, real(estimate, kind=kindvalue) )
        else
          call ConfidInterval_r4( real(values, kind=kindvalue), pct, intervals, &
            & mode )
        endif
      end subroutine ConfidInterval_int

      subroutine ConfidInterval_r4( values, &
        & pct, intervals, mode, estimate )
        integer, parameter                             :: KINDVALUE = r4
        real(KINDVALUE), dimension(:), intent(in)           :: values
        real(KINDVALUE), dimension(:), intent(in)           :: pct
        type(interval_t), dimension(:), intent(out)         :: intervals
        character(len=*), intent(in), optional              :: mode
        real(KINDVALUE), intent(in), optional               :: estimate
        ! Internal variables
        character(len=8)                                    :: myMode
        real(KINDVALUE), dimension(1)                       :: array2
        type(Stat_T), dimension(size(pct))                  :: gaps
        ! Executable
        myMode = 'mean'
        if ( present(mode) ) myMode = mode
        if ( lowercase(myMode(1:4)) == 'esti' .and. .not. present(estimate) ) &
          & call MLSMessage ( MLSMSG_Error, moduleName,  &
          & "This mode of ConfidInterval requires an explicit estimate " )
        if ( present(estimate) ) then
          myMode = 'esti'
          array2(1) = estimate
        endif
        call howfar_d1r4( values, array2, pct, gaps, mode )
        intervals%Bottom = gaps%min
        intervals%Top    = gaps%max
      end subroutine ConfidInterval_r4

      subroutine ConfidInterval_r8( values, &
        & pct, intervals, mode, estimate )
        integer, parameter                             :: KINDVALUE = r8
        real(KINDVALUE), dimension(:), intent(in)           :: values
        real(KINDVALUE), dimension(:), intent(in)           :: pct
        type(interval_t), dimension(:), intent(out)         :: intervals
        character(len=*), intent(in), optional              :: mode
        real(KINDVALUE), intent(in), optional               :: estimate
        ! Internal variables
        character(len=8)                                    :: myMode
        real(KINDVALUE), dimension(1)                       :: array2
        type(Stat_T), dimension(size(pct))                  :: gaps
        ! Executable
        myMode = 'mean'
        if ( present(mode) ) myMode = mode
        if ( lowercase(myMode(1:4)) == 'esti' .and. .not. present(estimate) ) &
          & call MLSMessage ( MLSMSG_Error, moduleName,  &
          & "This mode of ConfidInterval requires an explicit estimate " )
        if ( present(estimate) ) then
          myMode = 'esti'
          array2(1) = estimate
        endif
        call howfar_d1r8( values, array2, pct, gaps, mode )
        intervals%Bottom = gaps%min
        intervals%Top    = gaps%max
      end subroutine ConfidInterval_r8

      ! ------------------- howfar -----------------------
      ! This family of routines, keeping a given percentage of points
      ! in two arrays "near" each other in value,
      ! Returns the statistics of the remaining differences
      ! where nearness is defined as within a gap according to mode; either
      !  mode      gap
      !  'abs'  absolute value
      !  'med'  absolute value of difference of array1 from its median
      !           (ignoring array2)
      !  'rel'  relative value
      
      ! i.e., if absolute gap
      ! | array1(i) - array2(i) | < gap
      ! if relative gap
      ! | array1(i) - array2(i) | < gap * max( abs(array1(i)), abs(array2(i)) )
      
      ! ----------- Update ----------------------
      ! There are 5 possible values for mode
      !  value          meaning
      !  -----          -------
      !   abs           how far apart |array1 - array2|
      !   rel           how far apart |array1 - array2| / max(|array1|,|array2|)
      !   esti          set estimate = array2(1)
      !   mean          set estimate = mean(array1)
      !   medi          set estimate = median(array1)
      !
      ! Where an estimate is set, we find 2 values alpha, beta such that
      ! alpha < pct of array1 < beta
      !  gaps%min = alpha
      !  gaps%max = beta
      
      ! An inverse, in a sense, of the hownear procedures
      ! In another sense, this a different view of the histogram's
      ! bincount of our usual Stat_T type
      ! because if pct=90% => min, max of 90% of our date
      ! because if pct=10% => min, max of 10% of our date
      
      ! Some of the updated usage is for convenience in
      ! calculating confidence intervals (see ConfidInterval)
      
      ! Note that the statistic returned is inout--otherwise
      ! we would clobber nbins, bounds, bincount, etc.
      subroutine howfar_d1int( array1, array2, pct, gaps, mode, array1AtN, array2AtN )
        integer, parameter                             :: KINDVALUE = r4
        integer, dimension(:), intent(in)                   :: array1, array2
        real(KINDVALUE), dimension(:), intent(in)           :: pct
        type(Stat_T), dimension(:), intent(inout)           :: gaps
        character(len=*), intent(in)              :: mode ! 'rel' or 'abs'          
        integer, dimension(:,:), optional, intent(out)      :: array1AtN
        integer, dimension(:,:), optional, intent(out)      :: array2AtN
        real(KINDVALUE), dimension(3,size(pct))             :: rArray1AtN, rArray2AtN
        call howfar_d1r4( real(array1, KINDVALUE), real(array2, KINDVALUE), &
          & pct, gaps, mode, rArray1AtN, rArray2AtN )
        if ( present(array1AtN) ) array1AtN = rArray1AtN
        if ( present(array2AtN) ) array2AtN = rArray2AtN
      end subroutine howfar_d1int

      subroutine howfar_d2int( array1, array2, pct, gaps, mode, array1AtN, array2AtN )
        integer, parameter                             :: KINDVALUE = r4
        integer, dimension(:,:), intent(in)                 :: array1, array2
        real(KINDVALUE), dimension(:), intent(in)           :: pct
        type(Stat_T), dimension(:), intent(inout)           :: gaps
        character(len=*), intent(in)              :: mode ! 'rel' or 'abs'          
        integer, dimension(:,:), optional, intent(out)      :: array1AtN
        integer, dimension(:,:), optional, intent(out)      :: array2AtN
        real(KINDVALUE), dimension(3,size(pct))             :: rArray1AtN, rArray2AtN
        call howfar_d2r4( real(array1, KINDVALUE), real(array2, KINDVALUE), &
          & pct, gaps, mode, rArray1AtN, rArray2AtN )
        if ( present(array1AtN) ) array1AtN = rArray1AtN
        if ( present(array2AtN) ) array2AtN = rArray2AtN
      end subroutine howfar_d2int

      subroutine howfar_d3int( array1, array2, pct, gaps, mode, array1AtN, array2AtN )
        integer, parameter                             :: KINDVALUE = r4
        integer, dimension(:,:,:), intent(in)               :: array1, array2
        real(KINDVALUE), dimension(:), intent(in)           :: pct
        type(Stat_T), dimension(:), intent(inout)           :: gaps
        character(len=*), intent(in)              :: mode ! 'rel' or 'abs'          
        integer, dimension(:,:), optional, intent(out)      :: array1AtN
        integer, dimension(:,:), optional, intent(out)      :: array2AtN
        real(KINDVALUE), dimension(3,size(pct))             :: rArray1AtN, rArray2AtN
        call howfar_d3r4( real(array1, KINDVALUE), real(array2, KINDVALUE), &
          & pct, gaps, mode, rArray1AtN, rArray2AtN )
        if ( present(array1AtN) ) array1AtN = rArray1AtN
        if ( present(array2AtN) ) array2AtN = rArray2AtN
      end subroutine howfar_d3int

      subroutine howfar_d1r4( array1, array2, pct, gaps, mode, array1AtN, array2AtN )
        integer, parameter                             :: KINDVALUE = r4
        include 'howfar.f9h'
      end subroutine howfar_d1r4

      subroutine howfar_d2r4( array1, array2, pct, gaps, mode, array1AtN, array2AtN )
        integer, parameter                             :: KINDVALUE = r4
        real(KINDVALUE), dimension(:,:), intent(in)         :: array1, array2
        real(KINDVALUE), dimension(:), intent(in)           :: pct
        type(Stat_T), dimension(:), intent(inout)           :: gaps
        character(len=*), intent(in)              :: mode ! 'rel' or 'abs'          
        real(KINDVALUE), dimension(:, :), &
        & optional, intent(out)                 :: array1AtN, array2AtN
        integer, dimension(2)                          :: shp
        ! Executable
        shp = shape(array1)
        call howfar( reshape(array1, (/shp(1)*shp(2)/)), &
          & reshape(array2, (/shp(1)*shp(2)/)), &
          & pct, gaps, mode, array1AtN, array2AtN )
      end subroutine howfar_d2r4

      subroutine howfar_d3r4( array1, array2, pct, gaps, mode, array1AtN, array2AtN )
        integer, parameter                             :: KINDVALUE = r4
        real(KINDVALUE), dimension(:,:,:), intent(in)       :: array1, array2
        real(KINDVALUE), dimension(:), intent(in)           :: pct
        type(Stat_T), dimension(:), intent(inout)           :: gaps
        character(len=*), intent(in)              :: mode ! 'rel' or 'abs'          
        real(KINDVALUE), dimension(:, :), &
          & optional, intent(out)                 :: array1AtN, array2AtN
        integer, dimension(3)                          :: shp
        ! Executable
        shp = shape(array1)
        call howfar( reshape(array1, (/shp(1)*shp(2)*shp(3)/)), &
          & reshape(array2, (/shp(1)*shp(2)*shp(3)/)), &
          & pct, gaps, mode, array1AtN, array2AtN )
      end subroutine howfar_d3r4

      subroutine howfar_d4r4( array1, array2, pct, gaps, mode, array1AtN, array2AtN )
        integer, parameter                             :: KINDVALUE = r4
        real(KINDVALUE), dimension(:,:,:,:), intent(in)       :: array1, array2
        real(KINDVALUE), dimension(:), intent(in)           :: pct
        type(Stat_T), dimension(:), intent(inout)           :: gaps
        character(len=*), intent(in)              :: mode ! 'rel' or 'abs'          
        real(KINDVALUE), dimension(:, :), &
          & optional, intent(out)                 :: array1AtN, array2AtN
        ! Executable
        call howfar( reshape(array1, (/product(shape(array1))/)), &
          & reshape(array2, (/product(shape(array1))/)), &
          & pct, gaps, mode, array1AtN, array2AtN )
      end subroutine howfar_d4r4

      subroutine howfar_d1r8( array1, array2, pct, gaps, mode, array1AtN, array2AtN )
        integer, parameter                             :: KINDVALUE = r8
        include 'howfar.f9h'
      end subroutine howfar_d1r8

      subroutine howfar_d2r8( array1, array2, pct, gaps, mode, array1AtN, array2AtN )
        integer, parameter                             :: KINDVALUE = r8
        real(KINDVALUE), dimension(:,:), intent(in)         :: array1, array2
        real(KINDVALUE), dimension(:), intent(in)           :: pct
        type(Stat_T), dimension(:), intent(inout)           :: gaps
        character(len=*), intent(in)              :: mode ! 'rel' or 'abs'          
        real(KINDVALUE), dimension(:, :), &
          & optional, intent(out)                 :: array1AtN, array2AtN
        integer, dimension(2)                          :: shp
        ! Executable
        shp = shape(array1)
        call howfar( reshape(array1, (/shp(1)*shp(2)/)), &
          & reshape(array2, (/shp(1)*shp(2)/)), &
          & pct, gaps, mode, array1AtN, array2AtN )
      end subroutine howfar_d2r8

      subroutine howfar_d3r8( array1, array2, pct, gaps, mode, array1AtN, array2AtN )
        integer, parameter                             :: KINDVALUE = r8
        real(KINDVALUE), dimension(:,:,:), intent(in)       :: array1, array2
        real(KINDVALUE), dimension(:), intent(in)           :: pct
        type(Stat_T), dimension(:), intent(inout)           :: gaps
        character(len=*), intent(in)              :: mode ! 'rel' or 'abs'          
        real(KINDVALUE), dimension(:, :), &
          & optional, intent(out)                 :: array1AtN, array2AtN
        integer, dimension(3)                          :: shp
        ! Executable
        shp = shape(array1)
        call howfar( reshape(array1, (/shp(1)*shp(2)*shp(3)/)), &
          & reshape(array2, (/shp(1)*shp(2)*shp(3)/)), &
          & pct, gaps, mode, array1AtN, array2AtN )
      end subroutine howfar_d3r8

      subroutine howfar_d4r8( array1, array2, pct, gaps, mode, array1AtN, array2AtN )
        integer, parameter                             :: KINDVALUE = r8
        real(KINDVALUE), dimension(:,:,:,:), intent(in)       :: array1, array2
        real(KINDVALUE), dimension(:), intent(in)           :: pct
        type(Stat_T), dimension(:), intent(inout)           :: gaps
        character(len=*), intent(in)              :: mode ! 'rel' or 'abs'          
        real(KINDVALUE), dimension(:, :), &
          & optional, intent(out)                 :: array1AtN, array2AtN
        ! Executable
        call howfar( reshape(array1, (/product(shape(array1))/)), &
          & reshape(array2, (/product(shape(array2))/)), &
          & pct, gaps, mode, array1AtN, array2AtN )
      end subroutine howfar_d4r8

      ! ------------------- hownear -----------------------
      ! This family of routines finds the percentages of points
      ! in two arrays "near" each other in value
      ! where nearness is defined as within a gap either
      !  mode      gap
      !  'abs'  absolute value
      !  'rel'  relative value
      ! i.e., if absolute gap
      ! | array1(i) - array2(i) | < gap
      ! if relative gap
      ! | array1(i) - array2(i) | < gap * max( abs(array1(i)), abs(array2(i)) )
      subroutine hownear_d1int( array1, array2, pct, gaps, gapratios )
        integer, parameter                             :: KINDVALUE = r4
        integer, dimension(:), intent(in)                   :: array1, array2
        real(KINDVALUE), dimension(:), intent(out)          :: pct
        real(KINDVALUE), dimension(:), optional, intent(in) :: gaps
        real(KINDVALUE), dimension(:), optional, intent(in) :: gapratios
        call hownear_d1r4( real(array1, KINDVALUE), real(array2, KINDVALUE), &
          & pct, gaps, gapratios )
      end subroutine hownear_d1int

      subroutine hownear_d2int( array1, array2, pct, gaps, gapratios )
        integer, parameter                             :: KINDVALUE = r4
        integer, dimension(:,:), intent(in)                 :: array1, array2
        real(KINDVALUE), dimension(:), intent(out)          :: pct
        real(KINDVALUE), dimension(:), optional, intent(in) :: gaps
        real(KINDVALUE), dimension(:), optional, intent(in) :: gapratios
        call hownear_d2r4( real(array1, KINDVALUE), real(array2, KINDVALUE), &
          & pct, gaps, gapratios )
      end subroutine hownear_d2int

      subroutine hownear_d3int( array1, array2, pct, gaps, gapratios )
        integer, parameter                             :: KINDVALUE = r4
        integer, dimension(:,:,:), intent(in)               :: array1, array2
        real(KINDVALUE), dimension(:), intent(out)          :: pct
        real(KINDVALUE), dimension(:), optional, intent(in) :: gaps
        real(KINDVALUE), dimension(:), optional, intent(in) :: gapratios
        call hownear_d3r4( real(array1, KINDVALUE), real(array2, KINDVALUE), &
          & pct, gaps, gapratios )
      end subroutine hownear_d3int

      subroutine hownear_d1r4( array1, array2, pct, gaps, gapratios )
        integer, parameter                             :: KINDVALUE = r4
        include 'hownear.f9h'
      end subroutine hownear_d1r4

      subroutine hownear_d2r4( array1, array2, pct, gaps, gapratios )
        integer, parameter                             :: KINDVALUE = r4
        real(KINDVALUE), dimension(:,:), intent(in)         :: array1, array2
        real(KINDVALUE), dimension(:), intent(out)          :: pct
        real(KINDVALUE), dimension(:), optional, intent(in) :: gaps
        real(KINDVALUE), dimension(:), optional, intent(in) :: gapratios
        integer, dimension(2)                          :: shp
        ! Executable
        shp = shape(array1)
        call hownear( reshape(array1, (/shp(1)*shp(2)/)), &
          & reshape(array2, (/shp(1)*shp(2)/)), &
          & pct, gaps, gapratios )
      end subroutine hownear_d2r4

      subroutine hownear_d3r4( array1, array2, pct, gaps, gapratios )
        integer, parameter                             :: KINDVALUE = r4
        real(KINDVALUE), dimension(:,:,:), intent(in)       :: array1, array2
        real(KINDVALUE), dimension(:), intent(out)          :: pct
        real(KINDVALUE), dimension(:), optional, intent(in) :: gaps
        real(KINDVALUE), dimension(:), optional, intent(in) :: gapratios
        integer, dimension(3)                          :: shp
        ! Executable
        shp = shape(array1)
        call hownear( reshape(array1, (/shp(1)*shp(2)*shp(3)/)), &
          & reshape(array2, (/shp(1)*shp(2)*shp(3)/)), &
          & pct, gaps, gapratios )
      end subroutine hownear_d3r4

      subroutine hownear_d4r4( array1, array2, pct, gaps, gapratios )
        integer, parameter                             :: KINDVALUE = r4
        real(KINDVALUE), dimension(:,:,:,:), intent(in)       :: array1, array2
        real(KINDVALUE), dimension(:), intent(out)          :: pct
        real(KINDVALUE), dimension(:), optional, intent(in) :: gaps
        real(KINDVALUE), dimension(:), optional, intent(in) :: gapratios
        ! Executable
        call hownear( reshape(array1, (/product(shape(array1))/)), &
          & reshape(array2, (/product(shape(array2))/)), &
          & pct, gaps, gapratios )
      end subroutine hownear_d4r4

      subroutine hownear_d1r8( array1, array2, pct, gaps, gapratios )
        integer, parameter                             :: KINDVALUE = r8
        include 'hownear.f9h'
      end subroutine hownear_d1r8

      subroutine hownear_d2r8( array1, array2, pct, gaps, gapratios )
        integer, parameter                             :: KINDVALUE = r8
        real(KINDVALUE), dimension(:,:), intent(in)         :: array1, array2
        real(KINDVALUE), dimension(:), intent(out)          :: pct
        real(KINDVALUE), dimension(:), optional, intent(in) :: gaps
        real(KINDVALUE), dimension(:), optional, intent(in) :: gapratios
        integer, dimension(2)                          :: shp
        ! Executable
        shp = shape(array1)
        call hownear( reshape(array1, (/shp(1)*shp(2)/)), &
          & reshape(array2, (/shp(1)*shp(2)/)), &
          & pct, gaps, gapratios )
      end subroutine hownear_d2r8

      subroutine hownear_d3r8( array1, array2, pct, gaps, gapratios )
        integer, parameter                             :: KINDVALUE = r8
        real(KINDVALUE), dimension(:,:,:), intent(in)       :: array1, array2
        real(KINDVALUE), dimension(:), intent(out)          :: pct
        real(KINDVALUE), dimension(:), optional, intent(in) :: gaps
        real(KINDVALUE), dimension(:), optional, intent(in) :: gapratios
        integer, dimension(3)                          :: shp
        ! Executable
        shp = shape(array1)
        call hownear( reshape(array1, (/shp(1)*shp(2)*shp(3)/)), &
          & reshape(array2, (/shp(1)*shp(2)*shp(3)/)), &
          & pct, gaps, gapratios )
      end subroutine hownear_d3r8

      subroutine hownear_d4r8( array1, array2, pct, gaps, gapratios )
        integer, parameter                             :: KINDVALUE = r8
        real(KINDVALUE), dimension(:,:,:,:), intent(in)       :: array1, array2
        real(KINDVALUE), dimension(:), intent(out)          :: pct
        real(KINDVALUE), dimension(:), optional, intent(in) :: gaps
        real(KINDVALUE), dimension(:), optional, intent(in) :: gapratios
        ! Executable
        call hownear( reshape(array1, (/product(shape(array1))/)), &
          & reshape(array2, (/product(shape(array2))/)), &
          & pct, gaps, gapratios )
      end subroutine hownear_d4r8

      ! ------------------- pdf -----------------------
      ! This family of functions finds the probability density function
      ! for distribution of values saved in the data type statistic
      ! There are two modes of operation:
      ! 1-arg
      ! -----
      ! returns pdf(x, statistic) = f(x)

      ! 2-arg
      ! -----
      ! returns pdf(x1, x2, statistic) = Integral( f(x) dx, [x1, x2] )
      ! (which is actually cdf(x2) - cdf(x1))s

      ! Method:
      ! Let the bincount of the cell containing x be bincount(x)
      ! the width of that cell be h
      ! and the total counts over all bins be count
      ! Then
      ! pdf(x) = bincount(x) / (h count)

      ! Notes:
      ! (1) You must have first called the statistics subroutine
      ! to create the statistics data type
      ! (2) Statistics%nbins must be > 0
      function pdf_1(x, statistic) result(pdf)
        ! Given a value x and a statistic data type 
        ! (for which see statistics subroutine below)
        ! returns the probablility density function
        !
        ! Args
        real(r8), intent(in)     :: x
        type(Stat_T), intent(in) :: statistic
        real(r8)                 :: pdf
        ! Internal variables
        integer                         :: i
        real(r8)                        :: eta
        real(r8)                        :: h
        real(r8), dimension(:), pointer :: xbins
        ! Executable
        if ( statistic%nbins < 3 .or. .not. associated(statistic%bincount) ) then
          call MLSMessage ( MLSMSG_Error, moduleName,  &
            & "You must have binned the data before calling pdf" )
        endif
        nullify (xbins)
        call allocate_test(xbins, statistic%nbins-1, 'xbins', moduleName)
        xbins(1) = statistic%bounds(1)
        h = ( statistic%bounds(2) - statistic%bounds(1) ) / (statistic%nbins-2)
        do i=2, statistic%nbins-1
          eta = (i-1._r8) / (statistic%nbins-2)
          xbins(i) = (1._r8 - eta)*statistic%bounds(1) + eta*statistic%bounds(2)
        enddo
        i = FindFirst( x < xbins )
        i = max(i, 1)
        pdf = statistic%bincount(i) / (h*statistic%count)
        call deallocate_test(xbins, 'xbins', moduleName)
      end function pdf_1
      
      function pdf_2(x1, x2, statistic) result(pdf)
        ! Given values x1, x2,  and a statistic data type 
        ! (for which see statistics subroutine below)
        ! returns an integrated probablility density function
        !
        ! We want 
        !   Sum ( pdf(x) h, [x1 < x < x2] )
        !
        ! Args
        real(r8), intent(in)     :: x1
        real(r8), intent(in)     :: x2
        type(Stat_T), intent(in) :: statistic
        real(r8)                 :: pdf
        ! Internal variables
        integer                         :: i
        integer                         :: i1
        integer                         :: i2
        real(r8)                        :: eta
        real(r8)                        :: h
        real(r8)                        :: hi
        real(r8)                        :: pdfi
        real(r8), dimension(:), pointer :: xbins
        real(r8)                        :: xL
        real(r8)                        :: xR
        ! Executable
        if ( statistic%nbins < 3 .or. .not. associated(statistic%bincount) ) then
          call MLSMessage ( MLSMSG_Error, moduleName,  &
            & "You must have binned the data before calling pdf" )
        endif
        nullify (xbins)
        call allocate_test(xbins, statistic%nbins-1, 'xbins', moduleName)
        xbins(1) = statistic%bounds(1)
        h = ( statistic%bounds(2) - statistic%bounds(1) ) / (statistic%nbins-2)
        do i=2, statistic%nbins-1
          eta = (i-1._r8) / (statistic%nbins-2)
          xbins(i) = (1._r8 - eta)*statistic%bounds(1) + eta*statistic%bounds(2)
        enddo
        i1 = FindFirst( x1 < xbins )
        i2 = FindLast( x1 < xbins )
        i1 = max(i1, 1)
        i2 = max(i2, i1)
        pdf = 0.
        xR = min(x1, x2) - 1000*max(abs(x1), abs(x2)) ! This means -Infinity
        do i=i1, i2
          xL = xR
          if ( i < statistic%nbins ) then
            xL = xbins(i)
          else
            xR = max(x1, x2) + 1000*max(abs(x1), abs(x2)) ! This means +Infinity
          endif
          ! There are 4 possible cases to consider:
          ! (1) xL < x1 < xR < x2
          ! (2) x1 < xL < x2 < xR
          ! (3) x1 < xL < xR < x2
          ! (4) xL < x1 < x2 < xR
          hi = xR - xL
          pdfi = statistic%bincount(i) / (hi*statistic%count)
          if ( xL < x1 .and. xR < x2 ) then
            pdf = pdf + (xR-x1)*pdfi
          elseif ( x1 < xL .and. x2 < xR ) then
            pdf = pdf + (x2-xL)*pdfi
          elseif ( x1 < xL .and. xR < x2 ) then
            pdf = pdf + (xR-xL)*pdfi
          elseif ( xL < x1 .and. x2 < xR ) then
            pdf = pdf + (x2-x1)*pdfi
          else
            call MLSMessage ( MLSMSG_Error, moduleName,  &
            & "Unanticipated relation among xL, xR, x1, x2" )
          endif
        enddo
        call deallocate_test(xbins, 'xbins', moduleName)
      end function pdf_2
      
      ! ------------------- RATIOS -----------------------
      ! This family of routines finds both the absolute and relative RATIOS
      ! of one array w.r.t. another 'weighting' array (ignoring signs in both)
      
      ! The result, thinking in classic goldbrick style, would be

      ! **max absolute is exvalues(1) (corresponding ratio is exratios(1))**
      ! **max ratio is exratios(2) (corresponding absolute is exvalues(2))**

      ! E.g. let
      ! array1 = (/ 1,  -0.5, 2, +100 /)
      ! array2 = (/ 10, 0.5, 10, 1000 /)
      ! Then
      ! max value is 100 (corresponding ratio is 0.100) which is the 4th element
      ! max ratio is 1.0 (corresponding value is 0.5)   which is the 2nd element
      ! Thus we return
      ! exvalues = (/ 100,   0.5 /)
      ! exratios = (/ 0.100, 1.0 /)
      
      ! This may be useful in checking for non-negligible differences
      ! between two arrays representing results from two different analyses
      ! or noting the magnitude of changes to the goldbrick
      !
      ! If array2 is all 0, ratios, which would be undefined, are
      ! all set to -1.0
      
      ! optionally, first operates on array1 and array2, replacing them with
      ! op     replace array1 with               replace array2 with
      ! ---    -------------------               -------------------
      ! '-'      array1 - array2             max( abs(array1), abs(array2) )
      ! '+'      array1 + array2             max( abs(array1), abs(array2) )
      ! '*' min( abs(array1), abs(array2) )  max( abs(array1), abs(array2) )
      
      ! '-' is precisely what we would want for the goldbrick if 
      ! array1 is the reference standard
      ! array2 is the new result
      !
      ! Couldn't we have thought of a better name than "RATIOS"?
      subroutine ratios_d1int( array1, array2, exvalues, exratios, &
        & minratio, maxratio, meanratio, stddevratio, rmsratio, medianratio, &
        & fillValue, op )
        integer, parameter                             :: KINDVALUE = r4
        ! Args
        integer, dimension(:), intent(in)      :: array1, array2
        integer, dimension(2), intent(out)     :: exvalues
        integer, dimension(2), intent(out)     :: exratios
        integer, optional, intent(out)         :: minratio
        integer, optional, intent(out)         :: maxratio
        integer, optional, intent(out)         :: meanratio
        integer, optional, intent(out)         :: stddevratio
        integer, optional, intent(out)         :: rmsratio
        integer, optional, intent(out)         :: medianratio
        integer, optional                      :: FillValue
        character, optional, intent(in)        :: op
        ! Internal variables
        real(r4), dimension(2) :: rvalues
        real(r4), dimension(2) :: rratios
        real(KINDVALUE)        :: rmin, rmax, rmean, rstddev, rrms, rmedian
        call ratios_d1r4( real(array1, KINDVALUE), real(array2, KINDVALUE), &
          & rvalues, rratios, &
          & rmin, rmax, rmean, rstddev, rrms, rmedian, &
          & real(fillValue, KINDVALUE) )
        exvalues = rvalues
        exratios = rratios
        if ( present(minratio   ) ) minratio    = rmin
        if ( present(maxratio   ) ) maxratio    = rmax
        if ( present(meanratio  ) ) meanratio   = rmean
        if ( present(stddevratio) ) stddevratio = rstddev
        if ( present(rmsratio   ) ) rmsratio    = rrms
        if ( present(medianratio) ) medianratio = rmedian
      end subroutine ratios_d1int

      subroutine ratios_d2int( array1, array2, exvalues, exratios, &
        & minratio, maxratio, meanratio, stddevratio, rmsratio, medianratio, &
        & fillValue, op )
        integer, parameter                             :: KINDVALUE = r4
        ! Args
        integer, dimension(:,:), intent(in)    :: array1, array2
        integer, dimension(2), intent(out)     :: exvalues
        integer, dimension(2), intent(out)     :: exratios
        integer, optional, intent(out)         :: minratio
        integer, optional, intent(out)         :: maxratio
        integer, optional, intent(out)         :: meanratio
        integer, optional, intent(out)         :: stddevratio
        integer, optional, intent(out)         :: rmsratio
        integer, optional, intent(out)         :: medianratio
        integer, optional                      :: FillValue
        character, optional, intent(in)        :: op
        ! Internal variables
        real(r4), dimension(2) :: rvalues
        real(r4), dimension(2) :: rratios
        real(KINDVALUE)        :: rmin, rmax, rmean, rstddev, rrms, rmedian
        call ratios_d2r4( real(array1, KINDVALUE), real(array2, KINDVALUE), &
          & rvalues, rratios, &
          & rmin, rmax, rmean, rstddev, rrms, rmedian, &
          & real(fillValue, KINDVALUE) )
        exvalues = rvalues
        exratios = rratios
        if ( present(minratio   ) ) minratio    = rmin
        if ( present(maxratio   ) ) maxratio    = rmax
        if ( present(meanratio  ) ) meanratio   = rmean
        if ( present(stddevratio) ) stddevratio = rstddev
        if ( present(rmsratio   ) ) rmsratio    = rrms
        if ( present(medianratio) ) medianratio = rmedian
      end subroutine ratios_d2int

      subroutine ratios_d3int( array1, array2, exvalues, exratios, &
        & minratio, maxratio, meanratio, stddevratio, rmsratio, medianratio, &
        & fillValue, op )
        integer, parameter                             :: KINDVALUE = r4
        ! Args
        integer, dimension(:,:,:), intent(in)  :: array1, array2
        integer, dimension(2), intent(out)     :: exvalues
        integer, dimension(2), intent(out)     :: exratios
        integer, optional, intent(out)         :: minratio
        integer, optional, intent(out)         :: maxratio
        integer, optional, intent(out)         :: meanratio
        integer, optional, intent(out)         :: stddevratio
        integer, optional, intent(out)         :: rmsratio
        integer, optional, intent(out)         :: medianratio
        integer, optional                      :: FillValue
        character, optional, intent(in)        :: op
        ! Internal variables
        real(r4), dimension(2) :: rvalues
        real(r4), dimension(2) :: rratios
        real(KINDVALUE)        :: rmin, rmax, rmean, rstddev, rrms, rmedian
        call ratios_d3r4( real(array1, KINDVALUE), real(array2, KINDVALUE), &
          & rvalues, rratios, &
          & rmin, rmax, rmean, rstddev, rrms, rmedian, &
          & real(fillValue, KINDVALUE) )
        exvalues = rvalues
        exratios = rratios
        if ( present(minratio   ) ) minratio    = rmin
        if ( present(maxratio   ) ) maxratio    = rmax
        if ( present(meanratio  ) ) meanratio   = rmean
        if ( present(stddevratio) ) stddevratio = rstddev
        if ( present(rmsratio   ) ) rmsratio    = rrms
        if ( present(medianratio) ) medianratio = rmedian
      end subroutine ratios_d3int

      subroutine ratios_d1r4( array1, array2, exvalues, exratios, &
        & minratio, maxratio, meanratio, stddevratio, rmsratio, medianratio, &
        & fillValue, op )
        integer, parameter                             :: KINDVALUE = r4
        ! Args
        include 'ratios.f9h'
      end subroutine ratios_d1r4

      subroutine ratios_d1r8( array1, array2, exvalues, exratios, &
        & minratio, maxratio, meanratio, stddevratio, rmsratio, medianratio, &
        & fillValue, op )
        integer, parameter                             :: KINDVALUE = r8
        ! Args
        include 'ratios.f9h'
      end subroutine ratios_d1r8

      subroutine ratios_d2r4( array1, array2, exvalues, exratios, &
        & minratio, maxratio, meanratio, stddevratio, rmsratio, medianratio, &
        & fillValue, op )
        integer, parameter                             :: KINDVALUE = r4
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: array1, array2
        real(KINDVALUE), dimension(2), intent(out)     :: exvalues
        real(KINDVALUE), dimension(2), intent(out)     :: exratios
        real(KINDVALUE), optional, intent(out)         :: minratio
        real(KINDVALUE), optional, intent(out)         :: maxratio
        real(KINDVALUE), optional, intent(out)         :: meanratio
        real(KINDVALUE), optional, intent(out)         :: stddevratio
        real(KINDVALUE), optional, intent(out)         :: rmsratio
        real(KINDVALUE), optional, intent(out)         :: medianratio
        real(KINDVALUE), optional                      :: FillValue
        character, optional, intent(in)        :: op
        ! Internal variables
        integer, dimension(2)                          :: shp
        ! Executable
        shp = shape(array1)
        call ratios( reshape(array1, (/shp(1)*shp(2)/)), &
          & reshape(array2, (/shp(1)*shp(2)/)), &
          & exvalues, exratios, &
          & minratio, maxratio, meanratio, stddevratio, rmsratio, medianratio, &
          & fillValue, op )
      end subroutine ratios_d2r4

      subroutine ratios_d2r8( array1, array2, exvalues, exratios, &
        & minratio, maxratio, meanratio, stddevratio, rmsratio, medianratio, &
        & fillValue, op )
        integer, parameter                             :: KINDVALUE = r8
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: array1, array2
        real(KINDVALUE), dimension(2), intent(out)     :: exvalues
        real(KINDVALUE), dimension(2), intent(out)     :: exratios
        real(KINDVALUE), optional, intent(out)         :: minratio
        real(KINDVALUE), optional, intent(out)         :: maxratio
        real(KINDVALUE), optional, intent(out)         :: meanratio
        real(KINDVALUE), optional, intent(out)         :: stddevratio
        real(KINDVALUE), optional, intent(out)         :: rmsratio
        real(KINDVALUE), optional, intent(out)         :: medianratio
        real(KINDVALUE), optional                      :: FillValue
        character, optional, intent(in)        :: op
        ! Internal variables
        integer, dimension(2)                          :: shp
        ! Executable
        shp = shape(array1)
        call ratios( reshape(array1, (/shp(1)*shp(2)/)), &
          & reshape(array2, (/shp(1)*shp(2)/)), &
          & exvalues, exratios, &
          & minratio, maxratio, meanratio, stddevratio, rmsratio, medianratio, &
          & fillValue, op )
      end subroutine ratios_d2r8

      subroutine ratios_d3r4( array1, array2, exvalues, exratios, &
        & minratio, maxratio, meanratio, stddevratio, rmsratio, medianratio, &
        & fillValue, op )
        integer, parameter                             :: KINDVALUE = r4
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: array1, array2
        real(KINDVALUE), dimension(2), intent(out)     :: exvalues
        real(KINDVALUE), dimension(2), intent(out)     :: exratios
        real(KINDVALUE), optional, intent(out)         :: minratio
        real(KINDVALUE), optional, intent(out)         :: maxratio
        real(KINDVALUE), optional, intent(out)         :: meanratio
        real(KINDVALUE), optional, intent(out)         :: stddevratio
        real(KINDVALUE), optional, intent(out)         :: rmsratio
        real(KINDVALUE), optional, intent(out)         :: medianratio
        real(KINDVALUE), optional                      :: FillValue
        character, optional, intent(in)        :: op
        ! Internal variables
        integer, dimension(3)                          :: shp
        ! Executable
        shp = shape(array1)
        call ratios( reshape(array1, (/shp(1)*shp(2)*shp(3)/)), &
          & reshape(array2, (/shp(1)*shp(2)*shp(3)/)), &
          & exvalues, exratios, &
          & minratio, maxratio, meanratio, stddevratio, rmsratio, medianratio, &
          & fillValue, op )
      end subroutine ratios_d3r4

      subroutine ratios_d4r4( array1, array2, exvalues, exratios, &
        & minratio, maxratio, meanratio, stddevratio, rmsratio, medianratio, &
        & fillValue, op )
        integer, parameter                             :: KINDVALUE = r4
        ! Args
        real(KINDVALUE), dimension(:,:,:,:), intent(in)  :: array1, array2
        real(KINDVALUE), dimension(2), intent(out)     :: exvalues
        real(KINDVALUE), dimension(2), intent(out)     :: exratios
        real(KINDVALUE), optional, intent(out)         :: minratio
        real(KINDVALUE), optional, intent(out)         :: maxratio
        real(KINDVALUE), optional, intent(out)         :: meanratio
        real(KINDVALUE), optional, intent(out)         :: stddevratio
        real(KINDVALUE), optional, intent(out)         :: rmsratio
        real(KINDVALUE), optional, intent(out)         :: medianratio
        real(KINDVALUE), optional                      :: FillValue
        character, optional, intent(in)        :: op
        ! Executable
        call ratios( reshape(array1, (/product(shape(array1))/)), &
          & reshape(array2, (/product(shape(array2))/)), &
          & exvalues, exratios, &
          & minratio, maxratio, meanratio, stddevratio, rmsratio, medianratio, &
          & fillValue, op )
      end subroutine ratios_d4r4

      subroutine ratios_d3r8( array1, array2, exvalues, exratios, &
        & minratio, maxratio, meanratio, stddevratio, rmsratio, medianratio, &
        & fillValue, op )
        integer, parameter                             :: KINDVALUE = r8
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: array1, array2
        real(KINDVALUE), dimension(2), intent(out)     :: exvalues
        real(KINDVALUE), dimension(2), intent(out)     :: exratios
        real(KINDVALUE), optional, intent(out)         :: minratio
        real(KINDVALUE), optional, intent(out)         :: maxratio
        real(KINDVALUE), optional, intent(out)         :: meanratio
        real(KINDVALUE), optional, intent(out)         :: stddevratio
        real(KINDVALUE), optional, intent(out)         :: rmsratio
        real(KINDVALUE), optional, intent(out)         :: medianratio
        real(KINDVALUE), optional                      :: FillValue
        character, optional, intent(in)        :: op
        ! Internal variables
        integer, dimension(3)                          :: shp
        ! Executable
        shp = shape(array1)
        call ratios( reshape(array1, (/shp(1)*shp(2)*shp(3)/)), &
          & reshape(array2, (/shp(1)*shp(2)*shp(3)/)), &
          & exvalues, exratios, &
          & minratio, maxratio, meanratio, stddevratio, rmsratio, medianratio, &
          & fillValue, op )
      end subroutine ratios_d3r8
      
      subroutine ratios_d4r8( array1, array2, exvalues, exratios, &
        & minratio, maxratio, meanratio, stddevratio, rmsratio, medianratio, &
        & fillValue, op )
        integer, parameter                             :: KINDVALUE = r8
        ! Args
        real(KINDVALUE), dimension(:,:,:,:), intent(in)  :: array1, array2
        real(KINDVALUE), dimension(2), intent(out)     :: exvalues
        real(KINDVALUE), dimension(2), intent(out)     :: exratios
        real(KINDVALUE), optional, intent(out)         :: minratio
        real(KINDVALUE), optional, intent(out)         :: maxratio
        real(KINDVALUE), optional, intent(out)         :: meanratio
        real(KINDVALUE), optional, intent(out)         :: stddevratio
        real(KINDVALUE), optional, intent(out)         :: rmsratio
        real(KINDVALUE), optional, intent(out)         :: medianratio
        real(KINDVALUE), optional                      :: FillValue
        character, optional, intent(in)        :: op
        ! Executable
        call ratios( reshape(array1, (/product(shape(array1))/)), &
          & reshape(array2, (/product(shape(array2))/)), &
          & exvalues, exratios, &
          & minratio, maxratio, meanratio, stddevratio, rmsratio, medianratio, &
          & fillValue, op )
      end subroutine ratios_d4r8
      
      ! ------------------ SetUp -------------------------
      ! Sets up a stat_t by setting counters to zero 
      ! and optionally allocating a pointer
      subroutine SetUp( statistic, NBins, Bounds, verbose )
        ! Arg
        type(stat_t)                                 :: statistic
        integer, optional, intent(in)                :: NBins
        real(r8), dimension(2), optional, intent(in) :: Bounds
        logical, optional, intent(in)                :: verbose
        ! Internal variables
        integer                                      :: bin
        real(r8)                                     :: dx
        logical                                      :: myVerbose
        ! Executable
        myVerbose = .false.
        if ( present(verbose) ) myVerbose = verbose
        if ( myVerbose ) &
          & call Output ( 'Setting up a Stat_T data type', advance='yes' )
        statistic%count = 0
        statistic%fillcount = 0
        statistic%nbins = 0
        if ( associated(statistic%bincount) ) &
          & call deallocate_test( statistic%bincount, 'xbins', moduleName )
        if ( present(NBins) .and. present(Bounds) ) then
          if ( NBins < 1 ) return
          call allocate_test( statistic%bincount, NBins, 'xbins', moduleName )
          if ( .not. myVerbose ) return
          dx = (Bounds(2) - Bounds(1)) / max( 1, (NBins - 1) )
          call OutputNamedValue ( 'Containing bins', NBins )
          do bin=1, NBins
            call OutputNamedValue ( 'bin', (/ (bin-1)*dx, bin*dx /) )
          enddo
        endif
      end subroutine SetUp

      ! ------------------ reset -------------------------
      ! resets a stat_t by setting counters to zero and deallocating a pointer
      subroutine reset( statistic )
        ! Arg
        type(stat_t) :: statistic
        ! Executable
        statistic%count = 0
        statistic%fillcount = 0
        statistic%nbins = 0
        if ( associated(statistic%bincount) ) &
          & call deallocate_test( statistic%bincount, 'xbins', moduleName )
      end subroutine reset

      ! ------------------- statFunction -----------------------
      function statFunction(values, fillValue, precision) result(statistic)
        ! Given a 1-d set of dbl prec values, returns a Stat_t typed result
        ! Does not bother with bincount array
        ! Thus does not permit histogramming;
        ! for that see statistics subroutine below
        ! Args
        real(r8), dimension(:), intent(in)             :: values
        real(r8), optional, intent(in)                 :: fillValue
        real(r8), dimension(:), optional,  intent(in)  :: precision
        type(stat_T)                                   :: statistic
        call allstats(values, fillValue=fillValue, precision=precision, &
          & min=statistic%min, max=statistic%max, mean=statistic%mean, &
          & stddev=statistic%stddev, rms=statistic%rms)
      end function statFunction
      
      ! ------------------- statistics -----------------------
      subroutine statistics(values, statistic, fillValue, precision)
        ! Given a 1-d set of dbl prec values, calculates a Stat_t typed result
        ! Note that the caller has responsibility for 
        ! (1) setting nbins and bounds
        ! (2) Resetting count (or else all calls get merged together)
        ! (3) allocating bincount array (and deallocating later)
        ! Args
        real(r8), dimension(:), intent(in)            :: values
        real(r8), optional, intent(in)                :: fillValue
        type(stat_T), intent(inout)                   :: statistic
        real(r8), dimension(:), optional, intent(in)  :: precision
        ! Internal variables
        logical                                       :: addedData
        integer                                       :: nbins
        real(r8), dimension(2)                        :: bounds
        ! Executable
        ! default values
        addedData = statistic%count > 0
        nbins = 1
        if ( statistic%nbins > 0 ) nbins = statistic%nbins
        bounds = 1.d0
        if ( any(statistic%bounds /= 1.0d0) ) bounds = statistic%bounds
        MLSStat%count = statistic%count
        if ( associated(statistic%bincount) ) then
          call allstats(values, &
            & nbins, bounds, addeddata, fillValue, precision, &
            & statistic%count, statistic%fillcount, &
            & statistic%min, statistic%max, statistic%mean, &
            & statistic%stddev, statistic%rms, statistic%median, &
            & statistic%bincount, statistic%indexing)
        else
          call allstats( values, &
            & nbins, bounds, addeddata, fillValue, precision, &
            & count=statistic%count, min=statistic%min, max=statistic%max, &
            & mean=statistic%mean, &
            & stddev=statistic%stddev, rms=statistic%rms, &
            & median=statistic%median, indexing=statistic%indexing )
        endif
      end subroutine statistics
      
      ! ------------------- dump_all -----------------------
      subroutine dump_all( statistic, oneLine )
        ! Dumps all details of statistic
        type(stat_T), intent(in)         :: statistic
        logical, optional, intent(in)    :: oneLine
        ! 
        logical :: myOneLine
        !
        myOneLine = statsOnOneLine
        if ( present(oneLine) ) myOneLine = oneLine
        call output('count:  ')
        call output(statistic%count)
        call blanks(4)
        call output('max:    ')
        call output(statistic%max)
        call blanks(4)
        call output('min:    ')
        call output(statistic%min)
        call blanks(4)
        call output('median:    ')
        call output(statistic%median)
        if ( .not. myOneLine ) then
          call newline
        else
          call blanks(4)
        endif
        call output('mean:   ')
        call output(statistic%mean)
        call blanks(4)
        call output('stddev: ')
        call output(statistic%stddev)
        call blanks(4)
        call output('rms:    ')
        call output(statistic%rms)
        if ( myOneLine ) return
        call newline
        if ( any(statistic%indexing /= 0) .and. showIndexing ) then
          call output('(locations)  ')
          call output('max:    ')
          call output(statistic%indexing(1))
          call blanks(4)
          call output('min:    ')
          call output(statistic%indexing(2))
          call blanks(4)
          call output('median:    ')
          call output(statistic%indexing(3))
          call newline
        endif
        if ( statistic%nbins > 2 ) then
          call output('x1,x2: ')
          call output(statistic%bounds)
          call newline
          call output('bincounts: ')
          call newline
          call output(statistic%bincount)
          call newline
        endif
        if ( statistic%fillcount > 0 ) then
          call output('bad or filtered counts: ')
          call output(statistic%fillcount)
          call blanks(4)
          call output('( ')
          call output( (100.*statistic%fillcount) / &
            & (statistic%fillcount + statistic%count) )
          call output(' % )')
          call newline
        endif
      end subroutine dump_all
      
      ! ------------------- dump_selected -----------------------
      subroutine dump_selected( statistic, which, oneLine )
        ! Dumps selected details of statistic
        type(stat_T), intent(in)         :: statistic
        character(len=*), intent(in)     :: which ! E.g., 'max,min'; '*' means all
        logical, optional, intent(in)    :: oneLine
        ! 
        logical :: myOneLine
        !
        myOneLine = statsOnOneLine
        if ( present(oneLine) ) myOneLine = oneLine
        call dump_if_selected( statistic%count, which, 'count', 'no' )
        call dump_if_selected( statistic%max, which, 'max', 'no' )
        call dump_if_selected( statistic%min, which, 'min', 'no' )
        call dump_if_selected( statistic%mean, which, 'mean', 'no' )
        call dump_if_selected( statistic%median, which, 'median', 'no' )
        if ( .not. myOneLine ) then
          call newline
        else
          call blanks(4)
        endif
        call dump_if_selected( statistic%stddev, which, 'stddev', 'no' )
        call dump_if_selected( statistic%rms, which, 'rms', 'no' )
        if ( myOneLine ) return
        if ( showIndexing ) then
          call newline
          call dump_if_selected( statistic%indexing, which, 'indexing' )
        endif
        if ( statistic%nbins < 3 .or. index(which, 'bin') < 1 ) return
        call output('x1,x2: ')
        call output(statistic%bounds)
        call newline
        call output('bincounts: ')
        call newline
        call output(statistic%bincount)
        call newline
      end subroutine dump_selected
      
      ! ------------------- Private Procedures -----------------------
      ! ------------------- dump_if_selected -----------------------
      subroutine dump_if_selected_r8_array( what, fields, name )
        ! Dumps selected details of statistic
        real(r8), dimension(:), intent(in) :: what
        character(len=*), intent(in)       :: fields ! E.g., 'max,min'; '*' means all
        character(len=*), intent(in)       :: name
        ! 
        if ( .not. showMe( name, fields ) ) return
        call output(trim(name) // ': ')
        call newline
        call output(what)
        call newline
      end subroutine dump_if_selected_r8_array

      subroutine dump_if_selected_r8_scalar( what, fields, name, advance )
        ! Dumps selected details of statistic
        real(r8), intent(in)               :: what
        character(len=*), intent(in)       :: fields ! E.g., 'max,min'; '*' means all
        character(len=*), intent(in)       :: advance
        character(len=*), intent(in)       :: name
        ! 
        if ( .not. showMe( name, fields ) ) return
        call output(trim(name) // ': ')
        if ( advance == 'yes' ) call newline
        call output(what)
        if ( advance == 'yes' ) then
          call newline
        else
          call blanks(3)
        endif
      end subroutine dump_if_selected_r8_scalar

      subroutine dump_if_selected_int_array( what, fields, name )
        ! Dumps selected details of statistic
        integer, dimension(:), intent(in)  :: what
        character(len=*), intent(in)       :: fields ! E.g., 'max,min'; '*' means all
        character(len=*), intent(in)       :: name
        ! 
        if ( .not. showMe( name, fields ) ) return
        call output(trim(name) // ': ')
        call newline
        call output(what)
        call newline
      end subroutine dump_if_selected_int_array

      subroutine dump_if_selected_int_scalar( what, fields, name, advance )
        ! Dumps selected details of statistic
        integer, intent(in)                :: what
        character(len=*), intent(in)       :: fields ! E.g., 'max,min'; '*' means all
        character(len=*), intent(in)       :: advance
        character(len=*), intent(in)       :: name
        ! 
        if ( .not. showMe( name, fields ) ) return
        call output(trim(name) // ': ')
        if ( advance == 'yes' ) call newline
        call output(what)
        if ( advance == 'yes' ) then
          call newline
        else
          call blanks(3)
        endif
      end subroutine dump_if_selected_int_scalar

      logical function showMe( field, fields )                      
        ! Determine whether this field should be dumped or not               
        character(len=*), intent(in) :: field
        character(len=*), intent(in) :: fields
        !                                                                    
        if ( fields == '*' ) then                                     
          showMe = .true.                                                 
        else                                                                 
          showMe = ( index(LowerCase(fields), LowerCase(trim(field))) > 0 )  
        endif                                                                
      end function showMe                                                    

      ! ------------------- filterValues_r4 -----------------------
      subroutine filterValues_r4(values, XTAB, NX, fillValue, precision)
      ! Return an array filtered of any fillValues
      ! or where corresponding precision array < 0
      ! If neither fillValue, precision supplied, return all values
      ! If only fillValues in array, return length 1 array containing fillValue
      ! If all precisions < 0, return length 1 array containing -999.99
      ! In all cases, allocate Array of appropriate size
      ! Args
      real(r4), dimension(:), intent(in)             :: values
      real(r4), dimension(:), pointer                :: xtab
      integer, intent(out)                           :: NX
      real(r4), optional, intent(in)                 :: fillValue
      real(r4), dimension(:), optional, intent(in)   :: precision
      ! Internal variables
      integer, dimension(size(values))               :: which
      integer                                        :: i
      ! Executable
      nullify(xtab)
      if ( present(fillValue) ) then
        ! call findAll(values /= fillvalue, which, how_many=NX)
        select case ( fillvaluerelation )
        case ( '=' )
          ! call findAll(.not. isFillValue(values, fillvalue), which, how_many=NX)
          call findAll(values /= fillvalue, which, how_many=NX)
        case ( '<' )
          call findAll(values >= fillValue, which, how_many=NX)
        case ( '>' )
          call findAll(values <= fillValue, which, how_many=NX)
        case default
          ! call findAll(.not. isFillValue(values, fillvalue), which, how_many=NX)
          call findAll(values /= fillvalue, which, how_many=NX)
        end select
        if ( NX < 1 ) then
          call allocate_test(XTAB, 1, 'XTAB', moduleName)
          XTAB = fillValue
        else
          call allocate_test(XTAB, NX, 'XTAB', moduleName)
          do i=1, NX
            XTAB(i) = values(which(i))
          enddo
        endif
      elseif ( present(precision) ) then
        call findAll(precision >= 0._r4, which, how_many=NX)
        if ( NX < 1 ) then
          call allocate_test(XTAB, 1, 'XTAB', moduleName)
          XTAB = -999.99_r4
        else
          call allocate_test(XTAB, NX, 'XTAB', moduleName)
          do i=1, NX
            XTAB(i) = values(which(i))
          enddo
        endif
      else
        NX = size(values)
        call allocate_test(XTAB, NX, 'XTAB', moduleName)
        XTAB = values
      endif
      end subroutine filterValues_r4

      ! ------------------- filterValues_r8 -----------------------
      subroutine filterValues_r8(values, XTAB, NX, fillValue, precision)
      ! Return an array filtered of any fillValues
      ! or where corresponding precision array < 0
      ! If fillValue is absent, return all values
      ! If neither fillValue, precision supplied, return all values
      ! If only fillValues in array, return length 1 array containing fillValue
      ! If all precisions < 0, return length 1 array containing -999.99
      ! In all cases, allocate Array of appropriate size
      ! Args
      real(r8), dimension(:), intent(in)             :: values
      real(r8), dimension(:), pointer                :: xtab
      integer, intent(out)                           :: NX
      real(r8), optional, intent(in)                 :: fillValue
      real(r8), dimension(:), optional, intent(in)   :: precision
      ! Internal variables
      integer, dimension(size(values))               :: which
      integer                                        :: i
      ! Executable
      nullify(xtab)
      if ( present(fillValue) ) then
        ! call findAll(values /= fillvalue, which, how_many=NX)
        select case ( fillvaluerelation )
        case ( '=' )
          ! call findAll(.not. isFillValue(values, fillvalue), which, how_many=NX)
          call findAll(values /= fillvalue, which, how_many=NX)
        case ( '<' )
          call findAll(values >= fillValue, which, how_many=NX)
        case ( '>' )
          call findAll(values <= fillValue, which, how_many=NX)
        case default
          ! call findAll(.not. isFillValue(values, fillvalue), which, how_many=NX)
          call findAll(values /= fillvalue, which, how_many=NX)
        end select
        if ( NX < 1 ) then
          call allocate_test(XTAB, 1, 'XTAB', moduleName)
          XTAB = fillValue
        else
          call allocate_test(XTAB, NX, 'XTAB', moduleName)
          do i=1, NX
            XTAB(i) = values(which(i))
          enddo
        endif
      elseif ( present(precision) ) then
        call findAll(precision >= 0._r8, which, how_many=NX)
        if ( NX < 1 ) then
          call allocate_test(XTAB, 1, 'XTAB', moduleName)
          XTAB = -999.99_r8
        else
          call allocate_test(XTAB, NX, 'XTAB', moduleName)
          do i=1, NX
            XTAB(i) = values(which(i))
          enddo
        endif
      else
        NX = size(values)
        call allocate_test(XTAB, NX, 'XTAB', moduleName)
        XTAB = values
      endif
      end subroutine filterValues_r8

      ! ------------------- STAT1_r4 -----------------------
      subroutine STAT1_r4(XTAB, NX, STATS, IHIST, NCELLS, X1, X2)
      real(r4) ::             COUNT, DELTA, PREV
      real(r4) ::             SCALE, RSCALE, SCLNEW, STATS(5), SUMSCL
      real(r4) ::             TEMP, TEST, X, X1, X2, XMAX, XMEAN, XMIN
      real(r4), intent(in) :: XTAB(:)
      include 'stats1.f9h'
      end subroutine STAT1_r4

      ! ------------------- STAT1_r8 -----------------------
      subroutine STAT1_r8(XTAB, NX, STATS, IHIST, NCELLS, X1, X2)
      real(r8) ::             COUNT, DELTA, PREV
      real(r8) ::             SCALE, RSCALE, SCLNEW, STATS(5), SUMSCL
      real(r8) ::             TEMP, TEST, X, X1, X2, XMAX, XMEAN, XMIN
      real(r8), intent(in) :: XTAB(:)
      include 'stats1.f9h'
      end subroutine STAT1_r8
      
      ! This family shrinks an array, deleting indices from set delete
      function shrinkarray_int( array, delete, shrunken ) result(newsize)
        integer, dimension(:), intent(in)    :: array
        ! Local variables
        integer, dimension(:), intent(out)   :: shrunken
        ! Executable
        include 'shrinkArray.f9h'
      end function shrinkarray_int
      
      function shrinkarray_r4( array, delete, shrunken ) result(newsize)
        real(r4), dimension(:), intent(in)    :: array
        ! Local variables
        real(r4), dimension(:), intent(out)   :: shrunken
        ! Executable
        include 'shrinkArray.f9h'
      end function shrinkarray_r4
      
      function shrinkarray_r8( array, delete, shrunken ) result(newsize)
        real(r8), dimension(:), intent(in)    :: array
        ! Local variables
        real(r8), dimension(:), intent(out)   :: shrunken
        ! Executable
        include 'shrinkArray.f9h'
      end function shrinkarray_r8

      function shcount ( condition ) result ( howmany )
        ! args
        logical, dimension(:), intent(in)  :: condition
        integer                            :: howmany
        howmany = count(condition)
        
      end function shcount

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

!=============================================================================
end module MLSStats1
!=============================================================================

!
! $Log$
! Revision 2.29  2022/02/03 18:50:45  pwagner
! Fixed spelling error; capitalized some function names
!
! Revision 2.28  2019/04/04 23:50:20  pwagner
! Added ConfidInterval and SetUp routines
!
! Revision 2.27  2017/12/01 01:22:08  pwagner
! CamelCase use statements; improved comments
!
! Revision 2.26  2016/09/09 20:37:32  pwagner
! Improved explanation of ratios, including args and usage
!
! Revision 2.25  2014/01/09 00:24:29  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.24  2013/08/12 23:47:25  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.23  2011/07/26 20:43:51  pwagner
! Added some 4d interfaces
!
! Revision 2.22  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.21  2008/11/24 19:31:49  pwagner
! Less wasteful of memory; should not segment dault so often
!
! Revision 2.20  2008/05/23 01:15:29  pwagner
! New public variables to control dumps; e.g. statsOnOneLine
!
! Revision 2.19  2007/11/01 23:28:42  pwagner
! Added mlscount function
!
! Revision 2.18  2007/10/24 23:56:58  pwagner
! Changed name of component to 'indexing'
!
! Revision 2.17  2007/10/24 00:19:35  pwagner
! Added mindexes component to MLSStat to hold index of max, min, median
!
! Revision 2.16  2007/10/12 23:36:16  pwagner
! Added howfar procedures for comparing two arrays
!
! Revision 2.15  2007/09/13 21:07:44  pwagner
! Added hownear
!
! Revision 2.14  2007/07/17 00:25:07  pwagner
! Added ratios, reset
!
! Revision 2.13  2007/03/07 21:03:45  pwagner
! Avoiding isFillValue in filterValues (did not fix bug in isFillValue yet)
!
! Revision 2.12  2007/02/06 17:54:13  pwagner
! Correctly tracks fillcount; dumps as fillcount and as %
!
! Revision 2.11  2006/08/21 23:38:41  pwagner
! Added pdf function; speedier median algorithm
!
! Revision 2.10  2006/08/12 00:08:21  pwagner
! Corrected bug in filling median
!
! Revision 2.9  2006/07/11 00:22:16  pwagner
! Most mls.. functions can take integer arrays
!
! Revision 2.8  2006/03/08 01:13:38  pwagner
! Added median as statistical component
!
! Revision 2.7  2006/02/01 23:44:37  pwagner
! Added doDump option to allStats
!
! Revision 2.6  2005/12/16 00:04:29  pwagner
! Changes to reflect new MLSFillValues module
!
! Revision 2.5  2005/06/22 17:25:50  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.4  2005/03/24 21:16:40  pwagner
! Avoid assigning to undefined values
!
! Revision 2.3  2004/09/28 23:15:35  pwagner
! Uses isFillValue to allow slight tolerance
!
! Revision 2.2  2004/09/15 18:03:46  pwagner
! Added optional precision array
!
! Revision 2.1  2004/09/13 20:40:38  pwagner
! First commit
!
