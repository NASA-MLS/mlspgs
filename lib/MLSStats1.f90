! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module MLSStats1                 ! Calculate Min, Max, Mean, rms, std deviation
!=============================================================================
  use Allocate_Deallocate, only: allocate_Test, Deallocate_Test
  use MLSCommon, only: r4, r8
  use MLSSets, only: findAll
  use OUTPUT_M, only: BLANKS, NEWLINE, OUTPUT

  implicit none
  private
  
  public :: STAT_T, &
    & ALLSTATS, DUMP, MLSMIN, MLSMAX, MLSMEAN, MLSSTDDEV, MLSRMS, &
    & STATFUNCTION, STATISTICS
  
  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       & "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character(len=*), parameter, private :: ModuleName = &
       & "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

  ! The functions and subroutines are based on
  ! (1) The MATH77 library subroutine; and
  ! (2) optionally filtering out fill values from array entries
  ! Interfaces have been supplied so that most procedures accept
  ! arrays of rank 1, 2, or 3 and either sngl or dbl precision
  ! In addition to the usual min, max, mean, stddev functions
  ! an rms function has been created
  ! Missing so far is a median function, mode function, chi^2, 
  ! confidence level, etc.
  ! This is the main datatype, a stat.
  ! (Would making fillValue a component makes sense?)
  ! (How about a separate count of times fillValue had to be ignored?)

  type Stat_T
    integer :: count = 0    ! If > 0, merging data from prior run(s)
    real(r8) :: min
    real(r8) :: max
    real(r8) :: mean
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
  
  interface mlsstddev
    module procedure mlsstddev_d1r4, mlsstddev_d2r4, mlsstddev_d3r4
    module procedure mlsstddev_d1r8, mlsstddev_d2r8, mlsstddev_d3r8
  end interface
  
  interface mlsrms
    module procedure mlsrms_d1r4, mlsrms_d2r4, mlsrms_d3r4
    module procedure mlsrms_d1r8, mlsrms_d2r8, mlsrms_d3r8
  end interface
  
  real, parameter :: FAC = 64.0E0

  ! Function names
  ! Stored using sequence of integers
  integer, parameter              :: FN_MIN = 1
  integer, parameter              :: FN_MAX = FN_MIN + 1
  integer, parameter              :: FN_MEAN = FN_MAX + 1
  integer, parameter              :: FN_STDDEV = FN_MEAN + 1
  integer, parameter              :: FN_RMS = FN_STDDEV + 1

contains
      ! ------------------- allstats_d1r4 -----------------------
      subroutine allstats_d1r4(values, nbins, bounds, addedData, fillValue, &
        & count, min, max, mean, stddev, rms, bincount)
        ! Args
        real(r4), dimension(:), intent(in)             :: values
        integer, optional, intent(in)                  :: nbins
        real(r4), dimension(2), optional, intent(in)   :: bounds
        logical, optional, intent(in)                  :: addeddata
        real(r4), optional, intent(in)                 :: fillValue
        integer, optional, intent(inout)               :: count
        real(r4), optional, intent(out)                :: min
        real(r4), optional, intent(out)                :: max
        real(r4), optional, intent(out)                :: mean
        real(r4), optional, intent(out)                :: stddev
        real(r4), optional, intent(out)                :: rms
        integer, dimension(:), optional, intent(out)   :: bincount
        ! Internal variables
        real(r4), dimension(5)                         :: stats
        logical                                        :: myAddedData
        integer                                        :: ncells
        integer                                        :: n
        integer                                        :: nx
        real(r4)                                       :: x1, x2
        real(r4), dimension(:), pointer                :: xtab => null()
        real(r4)                                       :: absMu, sigma
        ! Executable
        ncells = 1
        if ( present(nbins) .and. present(bincount) ) ncells = nbins
        x1 = 1.
        x2 = 1.
        if ( present(bounds) ) then
          x1 = bounds(1)
          x2 = bounds(2)
        endif
        stats(1) = MLSStat%count
        stats(2) = MLSStat%min
        stats(3) = MLSStat%max
        stats(4) = MLSStat%mean
        stats(5) = MLSStat%stddev
        myAddedData = .false.
        if ( present(addedData) ) myAddedData = addedData
        if ( .not. myAddedData ) stats(1) = 0  ! Reset count to start again
        if ( present(fillValue) ) then
          call filterValues_r4(values, XTAB, NX, fillValue)
          call STAT1_r4(XTAB, NX, STATS, bincount, NCELLS, X1, X2)
          call Deallocate_test ( XTAB, 'XTAB', ModuleName )
        else
          call STAT1_r4(values, size(values), STATS, bincount, NCELLS, X1, X2)
        endif
        MLSStat%count  = stats(1)
        MLSStat%min    = stats(2)
        MLSStat%max    = stats(3)
        MLSStat%mean   = stats(4)
        MLSStat%stddev = stats(5)
        if ( present(count ) ) count  = stats(1)
        if ( present(min   ) ) min    = stats(2)
        if ( present(max   ) ) max    = stats(3)
        if ( present(mean  ) ) mean   = stats(4)
        if ( present(stddev) ) stddev = stats(5)
        if ( present(rms   ) ) then
          n = stats(1)
          absMu = ABS(stats(4))
          sigma = stats(5)
          if ( n == 0 .or. (absMu == 0. .and. sigma == 0.) ) then
            rms = 0.
          elseif ( sigma < absMu ) then
            rms = absMu*sqrt( ((n-1.)/n)*(sigma/absMu)**2 + 1. )
          else
            rms = sigma*sqrt( ((n-1.)/n) + (absMu/sigma)**2 )
          endif
        endif
        
      end subroutine allstats_d1r4

      ! ------------------- allstats_d1r8 -----------------------
      subroutine allstats_d1r8(values, nbins, bounds, addedData, fillValue, &
        & count, min, max, mean, stddev, rms, bincount)
        ! Args
        real(r8), dimension(:), intent(in)             :: values
        integer, optional, intent(in)                  :: nbins
        real(r8), dimension(2), optional, intent(in)   :: bounds
        logical, optional, intent(in)                  :: addeddata
        real(r8), optional, intent(in)                 :: fillValue
        integer, optional, intent(inout)               :: count
        real(r8), optional, intent(out)                :: min
        real(r8), optional, intent(out)                :: max
        real(r8), optional, intent(out)                :: mean
        real(r8), optional, intent(out)                :: stddev
        real(r8), optional, intent(out)                :: rms
        integer, dimension(:), optional, intent(out)   :: bincount
        ! Internal variables
        real(r8), dimension(5)                         :: stats
        logical                                        :: myAddedData
        integer                                        :: ncells
        integer                                        :: n
        integer                                        :: nx
        real(r8)                                       :: x1, x2
        real(r8), dimension(:), pointer                :: xtab => null()
        real(r8)                                       :: absMu, sigma
        ! Executable
        ncells = 1
        if ( present(nbins) .and. present(bincount) ) ncells = nbins
        x1 = 1.
        x2 = 1.
        if ( present(bounds) ) then
          x1 = bounds(1)
          x2 = bounds(2)
        endif
        stats(1) = MLSStat%count
        stats(2) = MLSStat%min
        stats(3) = MLSStat%max
        stats(4) = MLSStat%mean
        stats(5) = MLSStat%stddev
        myAddedData = .false.
        if ( present(addedData) ) myAddedData = addedData
        if ( .not. myAddedData ) stats(1) = 0  ! Reset count to start again
        if ( present(fillValue) ) then
          call filterValues_r8(values, XTAB, NX, fillValue)
          call STAT1_r8(XTAB, NX, STATS, bincount, NCELLS, X1, X2)
          call Deallocate_test ( XTAB, 'XTAB', ModuleName )
        else
          call STAT1_r8(values, size(values), STATS, bincount, NCELLS, X1, X2)
        endif
        MLSStat%count  = stats(1)
        MLSStat%min    = stats(2)
        MLSStat%max    = stats(3)
        MLSStat%mean   = stats(4)
        MLSStat%stddev = stats(5)
        if ( present(count ) ) count  = stats(1)
        if ( present(min   ) ) min    = stats(2)
        if ( present(max   ) ) max    = stats(3)
        if ( present(mean  ) ) mean   = stats(4)
        if ( present(stddev) ) stddev = stats(5)
        if ( present(rms   ) ) then
          n = stats(1)
          absMu = ABS(stats(4))
          sigma = stats(5)
          if ( n == 0 .or. (absMu == 0. .and. sigma == 0.) ) then
            rms = 0.
          elseif ( sigma < absMu ) then
            rms = absMu*sqrt( ((n-1.)/n)*(sigma/absMu)**2 + 1. )
          else
            rms = sigma*sqrt( ((n-1.)/n) + (absMu/sigma)**2 )
          endif
        endif
        
      end subroutine allstats_d1r8

      ! ------------------- allstats_d2r4 -----------------------
      subroutine allstats_d2r4(values, nbins, bounds, addedData, fillValue, &
        & count, min, max, mean, stddev, rms, bincount)
        ! Args
        real(r4), dimension(:,:), intent(in)           :: values
        integer, optional, intent(in)                  :: nbins
        real(r4), dimension(2), optional, intent(in)   :: bounds
        logical, optional, intent(in)                  :: addeddata
        real(r4), optional, intent(in)                 :: fillValue
        integer, optional, intent(inout)               :: count
        real(r4), optional, intent(out)                :: min
        real(r4), optional, intent(out)                :: max
        real(r4), optional, intent(out)                :: mean
        real(r4), optional, intent(out)                :: stddev
        real(r4), optional, intent(out)                :: rms
        integer, dimension(:), optional, intent(out)   :: bincount
        ! Internal variables
        integer, dimension(2)                          :: shp
        ! Executable
        shp =shape(values)
        call allstats_d1r4(reshape(values, (/shp(1)*shp(2)/)), &
          & nbins, bounds, addedData, fillValue, &
          & count, min, max, mean, stddev, rms, bincount)
      end subroutine allstats_d2r4

      ! ------------------- allstats_d2r8 -----------------------
      subroutine allstats_d2r8(values, nbins, bounds, addedData, fillValue, &
        & count, min, max, mean, stddev, rms, bincount)
        ! Args
        real(r8), dimension(:,:), intent(in)           :: values
        integer, optional, intent(in)                  :: nbins
        real(r8), dimension(2), optional, intent(in)   :: bounds
        logical, optional, intent(in)                  :: addeddata
        real(r8), optional, intent(in)                 :: fillValue
        integer, optional, intent(inout)               :: count
        real(r8), optional, intent(out)                :: min
        real(r8), optional, intent(out)                :: max
        real(r8), optional, intent(out)                :: mean
        real(r8), optional, intent(out)                :: stddev
        real(r8), optional, intent(out)                :: rms
        integer, dimension(:), optional, intent(out)   :: bincount
        ! Internal variables
        integer, dimension(2)                          :: shp
        ! Executable
        shp =shape(values)
        call allstats_d1r8(reshape(values, (/shp(1)*shp(2)/)), &
          & nbins, bounds, addedData, fillValue, &
          & count, min, max, mean, stddev, rms, bincount)
      end subroutine allstats_d2r8

      ! ------------------- allstats_d3r4 -----------------------
      subroutine allstats_d3r4(values, nbins, bounds, addedData, fillValue, &
        & count, min, max, mean, stddev, rms, bincount)
        ! Args
        real(r4), dimension(:,:,:), intent(in)         :: values
        integer, optional, intent(in)                  :: nbins
        real(r4), dimension(2), optional, intent(in)   :: bounds
        logical, optional, intent(in)                  :: addeddata
        real(r4), optional, intent(in)                 :: fillValue
        integer, optional, intent(inout)               :: count
        real(r4), optional, intent(out)                :: min
        real(r4), optional, intent(out)                :: max
        real(r4), optional, intent(out)                :: mean
        real(r4), optional, intent(out)                :: stddev
        real(r4), optional, intent(out)                :: rms
        integer, dimension(:), optional, intent(out)   :: bincount
        ! Internal variables
        integer, dimension(3)                          :: shp
        ! Executable
        shp =shape(values)
        call allstats_d1r4(reshape(values, (/shp(1)*shp(2)*shp(3)/)), &
          & nbins, bounds, addedData, fillValue, &
          & count, min, max, mean, stddev, rms, bincount)
      end subroutine allstats_d3r4

      ! ------------------- allstats_d3r8 -----------------------
      subroutine allstats_d3r8(values, nbins, bounds, addedData, fillValue, &
        & count, min, max, mean, stddev, rms, bincount)
        ! Args
        real(r8), dimension(:,:,:), intent(in)         :: values
        integer, optional, intent(in)                  :: nbins
        real(r8), dimension(2), optional, intent(in)   :: bounds
        logical, optional, intent(in)                  :: addeddata
        real(r8), optional, intent(in)                 :: fillValue
        integer, optional, intent(inout)               :: count
        real(r8), optional, intent(out)                :: min
        real(r8), optional, intent(out)                :: max
        real(r8), optional, intent(out)                :: mean
        real(r8), optional, intent(out)                :: stddev
        real(r8), optional, intent(out)                :: rms
        integer, dimension(:), optional, intent(out)   :: bincount
        ! Internal variables
        integer, dimension(3)                          :: shp
        ! Executable
        shp =shape(values)
        call allstats_d1r8(reshape(values, (/shp(1)*shp(2)*shp(3)/)), &
          & nbins, bounds, addedData, fillValue, &
          & count, min, max, mean, stddev, rms, bincount)
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
      function statFunction(values, fillValue) result(statistic)
        ! Given a 1-d set of dbl prec values, returns a Stat_t typed result
        ! Does not bother with bincount array
        ! Thus does not permit histogramming;
        ! for that see statistics subroutine below
        ! Args
        real(r8), dimension(:), intent(in)  :: values
        real(r8), optional, intent(in)      :: fillValue
        type(stat_T)                        :: statistic
        call allstats(values, fillValue=fillValue, &
          & min=statistic%min, max=statistic%max, mean=statistic%mean, &
          & stddev=statistic%stddev, rms=statistic%rms)
      end function statFunction
      
      ! ------------------- statistics -----------------------
      subroutine statistics(values, statistic, fillValue)
        ! Given a 1-d set of dbl prec values, calculates a Stat_t typed result
        ! Note that the caller has responsibility for 
        ! (1) setting nbins and bounds
        ! (2) Resetting count (or else all calls get merged together)
        ! (3) allocating bincount array (and deallocating later)
        ! Args
        real(r8), dimension(:), intent(in)  :: values
        real(r8), optional, intent(in)      :: fillValue
        type(stat_T), intent(inout)         :: statistic
        ! Internal variables
        logical                             :: addedData
        integer                             :: nbins
        real(r8), dimension(2)              :: bounds
        ! Executable
        ! default values
        addedData = statistic%count > 0
        nbins = 1
        if ( statistic%nbins > 0 ) nbins = statistic%nbins
        bounds = 1.d0
        if ( any(statistic%bounds /= 1.0d0) ) bounds = statistic%bounds
        MLSStat%count = statistic%count
        if ( associated(statistic%bincount) ) then
          call allstats(values, nbins, bounds, addeddata, fillValue, &
            & statistic%count, statistic%min, statistic%max, statistic%mean, &
            & statistic%stddev, statistic%rms, statistic%bincount)
        else
          call allstats(values, nbins, bounds, addeddata, fillValue, &
            & count=statistic%count, min=statistic%min, max=statistic%max, &
            & mean=statistic%mean, &
            & stddev=statistic%stddev, rms=statistic%rms)
        endif
      end subroutine statistics
      
      ! ------------------- dump -----------------------
      subroutine dump(statistic)
        ! Dumps all details of statistic
        type(stat_T), intent(in)         :: statistic
        call output('count:  ')
        call output(statistic%count)
        call blanks(4)
        call output('max:    ')
        call output(statistic%max)
        call blanks(4)
        call output('min:    ')
        call output(statistic%min)
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
      end subroutine dump
      
      ! ------------------- Private Procedures -----------------------
      ! ------------------- filterValues_r4 -----------------------
      subroutine filterValues_r4(values, XTAB, NX, fillValue)
      ! Return an array filtered of any fillValues
      ! If fillValue is absent, return all values
      ! If only fillValues in array, return length 1 array containing fillValue
      ! In all cases, allocate Array of appropriate size
      ! Args
      real(r4), dimension(:), intent(in)             :: values
      real(r4), dimension(:), pointer                :: xtab
      integer, intent(out)                           :: NX
      real(r4), optional, intent(in)                 :: fillValue
      ! Internal variables
      integer, dimension(size(values))               :: which
      integer                                        :: i
      ! Executable
      nullify(xtab)
      if ( present(fillValue) ) then
        call findAll(values /= fillvalue, which, how_many=NX)
        if ( NX < 1 ) then
          call allocate_test(XTAB, 1, 'XTAB', moduleName)
          XTAB = fillValue
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
      subroutine filterValues_r8(values, XTAB, NX, fillValue)
      ! Return an array filtered of any fillValues
      ! If fillValue is absent, return all values
      ! If only fillValues in array, return length 1 array containing fillValue
      ! In all cases, allocate Array of appropriate size
      ! Args
      real(r8), dimension(:), intent(in)             :: values
      real(r8), dimension(:), pointer                :: xtab
      integer, intent(out)                           :: NX
      real(r8), optional, intent(in)                 :: fillValue
      ! Internal variables
      integer, dimension(size(values))               :: which
      integer                                        :: i
      ! Executable
      nullify(xtab)
      if ( present(fillValue) ) then
        call findAll(values /= fillvalue, which, how_many=NX)
        if ( NX < 1 ) then
          call allocate_test(XTAB, 1, 'XTAB', moduleName)
          XTAB = fillValue
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

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

!=============================================================================
end module MLSStats1
!=============================================================================

!
! $Log$
! Revision 2.1  2004/09/13 20:40:38  pwagner
! First commit
!
