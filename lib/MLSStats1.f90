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
module MLSStats1                 ! Calculate Min, Max, Mean, rms, std deviation
!=============================================================================
  use Allocate_Deallocate, only: allocate_Test, Deallocate_Test
  use MLSCommon, only: r4, r8
  use MLSFillValues, only: isFillValue
  use MLSSets, only: findAll
  use MLSStringLists, only: catLists
  use MLSStrings, only: lowerCase
  use OUTPUT_M, only: BLANKS, NEWLINE, OUTPUT

  implicit none
  private
  
  public :: STAT_T             ! The data type
  public :: ALLSTATS, DUMP, STATISTICS  ! subroutines
  public :: MLSMIN, MLSMAX, MLSMEAN, MLSMEDIAN, MLSSTDDEV, MLSRMS ! functions
  public :: STATFUNCTION
  
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! The functions and subroutines are based on
  ! (1) The MATH77 library subroutine; and
  ! (2) optionally filtering out fill values from array entries
  ! Optionally, instead of filtering fill values, an array of precisions,
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

  ! This is the main datatype, a stat.
  ! (Would making fillValue a component makes sense?)

  type Stat_T
    integer :: count = 0      ! If > 0, merging data from prior call(s)
    integer :: fillcount = 0  ! Number of times fillValues ignored
    real(r8) :: max
    real(r8) :: mean
    real(r8) :: median        ! incorrect if merging data
    real(r8) :: min
    real(r8) :: stddev
    real(r8) :: rms
    ! The next 3 deal with histogramming data
    integer :: nbins = 0  !NCELLS; put > 2 if histogramming
    real(r8), dimension(2) :: bounds = 1.d0 ! X1, X2; put X2 > X1
    integer, dimension(:), pointer :: bincount => null() ! IHIST
  end type Stat_T
  
  type(Stat_T), save :: MLSStat
  
  interface ALLSTATS
    module procedure allstats_d1r4, allstats_d2r4, allstats_d3r4
    module procedure allstats_d1r8, allstats_d2r8, allstats_d3r8
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
  
  interface mlsmin
    module procedure mlsmin_d1r4, mlsmin_d2r4, mlsmin_d3r4
    module procedure mlsmin_d1r8, mlsmin_d2r8, mlsmin_d3r8
  end interface
  
  interface mlsmax
    module procedure mlsmax_d1r4, mlsmax_d2r4, mlsmax_d3r4
    module procedure mlsmax_d1r8, mlsmax_d2r8, mlsmax_d3r8
  end interface
  
  interface mlsmean
    module procedure mlsmean_d1r4, mlsmean_d2r4, mlsmean_d3r4
    module procedure mlsmean_d1r8, mlsmean_d2r8, mlsmean_d3r8
  end interface
  
  interface mlsmedian
    module procedure mlsmedian_d1r4, mlsmedian_d2r4, mlsmedian_d3r4
    module procedure mlsmedian_d1r8, mlsmedian_d2r8, mlsmedian_d3r8
  end interface
  
  interface mlsstddev
    module procedure mlsstddev_d1r4, mlsstddev_d2r4, mlsstddev_d3r4
    module procedure mlsstddev_d1r8, mlsstddev_d2r8, mlsstddev_d3r8
  end interface
  
  interface mlsrms
    module procedure mlsrms_d1r4, mlsrms_d2r4, mlsrms_d3r4
    module procedure mlsrms_d1r8, mlsrms_d2r8, mlsrms_d3r8
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
  integer, parameter              :: FN_MIN = 1
  integer, parameter              :: FN_MAX = FN_MIN + 1
  integer, parameter              :: FN_MEAN = FN_MAX + 1
  integer, parameter              :: FN_STDDEV = FN_MEAN + 1
  integer, parameter              :: FN_RMS = FN_STDDEV + 1
  integer, parameter              :: FN_MEDIAN = FN_RMS + 1

contains
      ! ------------------- allstats_d1r4 -----------------------
      subroutine allstats_d1r4( values, &
        & nbins, bounds, addedData, fillValue, precision, &
        & count, min, max, mean, stddev, rms, median, bincount, doDump )
        integer, parameter :: rk = r4
        include 'allstats_d1.f9h'
      end subroutine allstats_d1r4

      ! ------------------- allstats_d1r8 -----------------------
      subroutine allstats_d1r8( values, &
        & nbins, bounds, addedData, fillValue, precision, &
        & count, min, max, mean, stddev, rms, median, bincount, doDump )
        integer, parameter :: rk = r8
        include 'allstats_d1.f9h'
      end subroutine allstats_d1r8

      ! ------------------- allstats_d2r4 -----------------------
      subroutine allstats_d2r4( values, &
        & nbins, bounds, addedData, fillValue, precision, &
        & count, min, max, mean, stddev, rms, median, bincount, doDump )
        ! Args
        real(r4), dimension(:,:), intent(in)           :: values
        integer, optional, intent(in)                  :: nbins
        real(r4), dimension(2), optional, intent(in)   :: bounds
        logical, optional, intent(in)                  :: addeddata
        real(r4), optional, intent(in)                 :: fillValue
        real(r4), dimension(:,:), optional, intent(in) :: precision
        integer, optional, intent(inout)               :: count
        real(r4), optional, intent(out)                :: min
        real(r4), optional, intent(out)                :: max
        real(r4), optional, intent(out)                :: mean
        real(r4), optional, intent(out)                :: stddev
        real(r4), optional, intent(out)                :: rms
        real(r4), optional, intent(out)                :: median
        integer, dimension(:), optional, intent(out)   :: bincount
        logical , optional, intent(in)                 :: doDump
        ! Internal variables
        integer, dimension(2)                          :: shp
        ! Executable
        shp =shape(values)
        if ( .not. present(precision) ) then
          call allstats_d1r4(reshape(values, (/shp(1)*shp(2)/)), &
            & nbins, bounds, addedData, fillValue, &
            & count=count, min=min, max=max, mean=mean, &
            & stddev=stddev, rms=rms, median=median, bincount=bincount, &
            & doDump=doDump)
        else
          call allstats_d1r4(reshape(values, (/shp(1)*shp(2)/)), &
            & nbins, bounds, addedData, fillValue, &
            & reshape(precision, (/shp(1)*shp(2)/)),&
            & count=count, min=min, max=max, mean=mean, &
            & stddev=stddev, rms=rms, median=median, bincount=bincount, &
            & doDump=doDump)
        endif
      end subroutine allstats_d2r4

      ! ------------------- allstats_d2r8 -----------------------
      subroutine allstats_d2r8( values, &
        & nbins, bounds, addedData, fillValue, precision, &
        & count, min, max, mean, stddev, rms, median, bincount, doDump )
        ! Args
        real(r8), dimension(:,:), intent(in)           :: values
        integer, optional, intent(in)                  :: nbins
        real(r8), dimension(2), optional, intent(in)   :: bounds
        logical, optional, intent(in)                  :: addeddata
        real(r8), optional, intent(in)                 :: fillValue
        real(r8), dimension(:,:), optional, intent(in) :: precision
        integer, optional, intent(inout)               :: count
        real(r8), optional, intent(out)                :: min
        real(r8), optional, intent(out)                :: max
        real(r8), optional, intent(out)                :: mean
        real(r8), optional, intent(out)                :: stddev
        real(r8), optional, intent(out)                :: rms
        real(r8), optional, intent(out)                :: median
        integer, dimension(:), optional, intent(out)   :: bincount
        logical , optional, intent(in)                 :: doDump
        ! Internal variables
        integer, dimension(2)                          :: shp
        ! Executable
        shp =shape(values)
        if ( .not. present(precision) ) then
          call allstats_d1r8(reshape(values, (/shp(1)*shp(2)/)), &
            & nbins, bounds, addedData, fillValue, &
            & count=count, min=min, max=max, mean=mean, &
            & stddev=stddev, rms=rms, median=median, bincount=bincount, &
            & doDump=doDump)
        else
          call allstats_d1r8(reshape(values, (/shp(1)*shp(2)/)), &
            & nbins, bounds, addedData, fillValue, &
            & reshape(precision, (/shp(1)*shp(2)/)),&
            & count=count, min=min, max=max, mean=mean, &
            & stddev=stddev, rms=rms, median=median, bincount=bincount, &
            & doDump=doDump)
        endif
      end subroutine allstats_d2r8

      ! ------------------- allstats_d3r4 -----------------------
      subroutine allstats_d3r4( values, &
        & nbins, bounds, addedData, fillValue, precision, &
        & count, min, max, mean, stddev, rms, median, bincount, doDump )
        ! Args
        real(r4), dimension(:,:,:), intent(in)         :: values
        integer, optional, intent(in)                  :: nbins
        real(r4), dimension(2), optional, intent(in)   :: bounds
        logical, optional, intent(in)                  :: addeddata
        real(r4), optional, intent(in)                 :: fillValue
        real(r4), dimension(:,:,:), optional, intent(in)   :: precision
        integer, optional, intent(inout)               :: count
        real(r4), optional, intent(out)                :: min
        real(r4), optional, intent(out)                :: max
        real(r4), optional, intent(out)                :: mean
        real(r4), optional, intent(out)                :: stddev
        real(r4), optional, intent(out)                :: rms
        real(r4), optional, intent(out)                :: median
        integer, dimension(:), optional, intent(out)   :: bincount
        logical , optional, intent(in)                 :: doDump
        ! Internal variables
        integer, dimension(3)                          :: shp
        ! Executable
        shp =shape(values)
        if ( .not. present(precision) ) then
          call allstats_d1r4(reshape(values, (/shp(1)*shp(2)*shp(3)/)), &
            & nbins, bounds, addedData, fillValue, &
            & count=count, min=min, max=max, mean=mean, &
            & stddev=stddev, rms=rms, median=median, bincount=bincount, &
            & doDump=doDump)
        else
          call allstats_d1r4(reshape(values, (/shp(1)*shp(2)*shp(3)/)), &
            & nbins, bounds, addedData, fillValue, &
            & reshape(precision, (/shp(1)*shp(2)/)),&
            & count=count, min=min, max=max, mean=mean, &
            & stddev=stddev, rms=rms, median=median, bincount=bincount, &
            & doDump=doDump)
        endif
      end subroutine allstats_d3r4

      ! ------------------- allstats_d3r8 -----------------------
      subroutine allstats_d3r8( values, &
        & nbins, bounds, addedData, fillValue, precision, &
        & count, min, max, mean, stddev, rms, median, bincount, doDump )
        ! Args
        real(r8), dimension(:,:,:), intent(in)         :: values
        integer, optional, intent(in)                  :: nbins
        real(r8), dimension(2), optional, intent(in)   :: bounds
        logical, optional, intent(in)                  :: addeddata
        real(r8), optional, intent(in)                 :: fillValue
        real(r8), dimension(:,:,:), optional, intent(in)   :: precision
        integer, optional, intent(inout)               :: count
        real(r8), optional, intent(out)                :: min
        real(r8), optional, intent(out)                :: max
        real(r8), optional, intent(out)                :: mean
        real(r8), optional, intent(out)                :: stddev
        real(r8), optional, intent(out)                :: rms
        real(r8), optional, intent(out)                :: median
        integer, dimension(:), optional, intent(out)   :: bincount
        logical , optional, intent(in)                 :: doDump
        ! Internal variables
        integer, dimension(3)                          :: shp
        ! Executable
        shp =shape(values)
        if ( .not. present(precision) ) then
          call allstats_d1r8(reshape(values, (/shp(1)*shp(2)*shp(3)/)), &
            & nbins, bounds, addedData, fillValue, &
            & count=count, min=min, max=max, mean=mean, &
            & stddev=stddev, rms=rms, median=median, bincount=bincount, &
            & doDump=doDump)
        else
          call allstats_d1r8(reshape(values, (/shp(1)*shp(2)*shp(3)/)), &
            & nbins, bounds, addedData, fillValue, &
            & reshape(precision, (/shp(1)*shp(2)/)),&
            & count=count, min=min, max=max, mean=mean, &
            & stddev=stddev, rms=rms, median=median, bincount=bincount, &
            & doDump=doDump)
        endif
      end subroutine allstats_d3r8

      ! ------------------- mlsmin_d1r4 -----------------------
      function mlsmin_d1r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_MIN
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function mlsmin_d1r4

      ! ------------------- mlsmin_d1r8 -----------------------
      function mlsmin_d1r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_MIN
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function mlsmin_d1r8

      ! ------------------- mlsmin_d2r4 -----------------------
      function mlsmin_d2r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_MIN
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function mlsmin_d2r4

      ! ------------------- mlsmin_d2r8 -----------------------
      function mlsmin_d2r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_MIN
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function mlsmin_d2r8

      ! ------------------- mlsmin_d3r4 -----------------------
      function mlsmin_d3r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_MIN
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function mlsmin_d3r4

      ! ------------------- mlsmin_d3r8 -----------------------
      function mlsmin_d3r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_MIN
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function mlsmin_d3r8

      ! ------------------- mlsmax_d1r4 -----------------------
      function mlsmax_d1r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_max
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function mlsmax_d1r4

      ! ------------------- mlsmax_d1r8 -----------------------
      function mlsmax_d1r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_max
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function mlsmax_d1r8

      ! ------------------- mlsmax_d2r4 -----------------------
      function mlsmax_d2r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_max
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function mlsmax_d2r4

      ! ------------------- mlsmax_d2r8 -----------------------
      function mlsmax_d2r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_max
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function mlsmax_d2r8

      ! ------------------- mlsmax_d3r4 -----------------------
      function mlsmax_d3r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_max
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function mlsmax_d3r4

      ! ------------------- mlsmax_d3r8 -----------------------
      function mlsmax_d3r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_max
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function mlsmax_d3r8

      ! ------------------- mlsmean_d1r4 -----------------------
      function mlsmean_d1r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_mean
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function mlsmean_d1r4

      ! ------------------- mlsmean_d1r8 -----------------------
      function mlsmean_d1r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_mean
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function mlsmean_d1r8

      ! ------------------- mlsmean_d2r4 -----------------------
      function mlsmean_d2r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_mean
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function mlsmean_d2r4

      ! ------------------- mlsmean_d2r8 -----------------------
      function mlsmean_d2r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_mean
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function mlsmean_d2r8

      ! ------------------- mlsmean_d3r4 -----------------------
      function mlsmean_d3r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_mean
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function mlsmean_d3r4

      ! ------------------- mlsmean_d3r8 -----------------------
      function mlsmean_d3r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_mean
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function mlsmean_d3r8

      ! ------------------- mlsmedian_d1r4 -----------------------
      function mlsmedian_d1r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_median
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function mlsmedian_d1r4

      ! ------------------- mlsmedian_d1r8 -----------------------
      function mlsmedian_d1r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_median
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function mlsmedian_d1r8

      ! ------------------- mlsmedian_d2r4 -----------------------
      function mlsmedian_d2r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_median
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function mlsmedian_d2r4

      ! ------------------- mlsmedian_d2r8 -----------------------
      function mlsmedian_d2r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_median
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function mlsmedian_d2r8

      ! ------------------- mlsmedian_d3r4 -----------------------
      function mlsmedian_d3r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_median
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function mlsmedian_d3r4

      ! ------------------- mlsmedian_d3r8 -----------------------
      function mlsmedian_d3r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_median
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function mlsmedian_d3r8

      ! ------------------- mlsstddev_d1r4 -----------------------
      function mlsstddev_d1r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_stddev
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function mlsstddev_d1r4

      ! ------------------- mlsstddev_d1r8 -----------------------
      function mlsstddev_d1r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_stddev
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function mlsstddev_d1r8

      ! ------------------- mlsstddev_d2r4 -----------------------
      function mlsstddev_d2r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_stddev
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function mlsstddev_d2r4

      ! ------------------- mlsstddev_d2r8 -----------------------
      function mlsstddev_d2r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_stddev
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function mlsstddev_d2r8

      ! ------------------- mlsstddev_d3r4 -----------------------
      function mlsstddev_d3r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_stddev
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function mlsstddev_d3r4

      ! ------------------- mlsstddev_d3r8 -----------------------
      function mlsstddev_d3r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_stddev
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function mlsstddev_d3r8

      ! ------------------- mlsrms_d1r4 -----------------------
      function mlsrms_d1r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_rms
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function mlsrms_d1r4

      ! ------------------- mlsrms_d1r8 -----------------------
      function mlsrms_d1r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_rms
        ! Args
        real(KINDVALUE), dimension(:), intent(in)      :: values
        include 'stats0.f9h'
      end function mlsrms_d1r8

      ! ------------------- mlsrms_d2r4 -----------------------
      function mlsrms_d2r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_rms
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function mlsrms_d2r4

      ! ------------------- mlsrms_d2r8 -----------------------
      function mlsrms_d2r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_rms
        ! Args
        real(KINDVALUE), dimension(:,:), intent(in)    :: values
        include 'stats0.f9h'
      end function mlsrms_d2r8

      ! ------------------- mlsrms_d3r4 -----------------------
      function mlsrms_d3r4(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r4
        integer, parameter                             :: FN = FN_rms
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function mlsrms_d3r4

      ! ------------------- mlsrms_d3r8 -----------------------
      function mlsrms_d3r8(values, fillValue) result(rValue)
        integer, parameter                             :: KINDVALUE = r8
        integer, parameter                             :: FN = FN_rms
        ! Args
        real(KINDVALUE), dimension(:,:,:), intent(in)  :: values
        include 'stats0.f9h'
      end function mlsrms_d3r8
      
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
            & statistic%count, statistic%min, statistic%max, statistic%mean, &
            & statistic%stddev, statistic%rms, statistic%median, &
            & statistic%bincount)
        else
          call allstats(values, &
            & nbins, bounds, addeddata, fillValue, precision, &
            & count=statistic%count, min=statistic%min, max=statistic%max, &
            & mean=statistic%mean, &
            & stddev=statistic%stddev, rms=statistic%rms, median=statistic%rms )
        endif
      end subroutine statistics
      
      ! ------------------- dump_all -----------------------
      subroutine dump_all(statistic)
        ! Dumps all details of statistic
        type(stat_T), intent(in)         :: statistic
        ! 
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
        call newline
        call output('mean:   ')
        call output(statistic%mean)
        call blanks(4)
        call output('stddev: ')
        call output(statistic%stddev)
        call blanks(4)
        call output('rms:    ')
        call output(statistic%rms)
        call newline
        if ( statistic%nbins < 3 ) return
        call output('x1,x2: ')
        call output(statistic%bounds)
        call newline
        call output('bincounts: ')
        call newline
        call output(statistic%bincount)
        call newline
      end subroutine dump_all
      
      ! ------------------- dump_selected -----------------------
      subroutine dump_selected( statistic, which )
        ! Dumps selected details of statistic
        type(stat_T), intent(in)         :: statistic
        character(len=*), intent(in)     :: which ! E.g., 'max,min'; '*' means all
        ! 
        call dump_if_selected( statistic%count, which, 'count', 'no' )
        call dump_if_selected( statistic%max, which, 'max', 'no' )
        call dump_if_selected( statistic%min, which, 'min', 'no' )
        call dump_if_selected( statistic%mean, which, 'mean', 'no' )
        call dump_if_selected( statistic%median, which, 'median', 'no' )
        call newline
        call dump_if_selected( statistic%stddev, which, 'stddev', 'no' )
        call dump_if_selected( statistic%rms, which, 'rms', 'no' )
        call newline
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
        call findAll(.not. isFillValue(values, fillvalue), which, how_many=NX)
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
        call findAll(.not. isFillValue(values, fillvalue), which, how_many=NX)
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

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

!=============================================================================
end module MLSStats1
!=============================================================================

!
! $Log$
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
