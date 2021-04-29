! Copyright 2017, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module HyperSlabs              ! .. and array gymnastics
!=============================================================================

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use MLSFillValues, only: FilterValues
  use MLSCommon, only: MLS_HyperStart
  use MLSFinds, only: FindAll, FindFirst, FindLast
  use MLSKinds ! Everything
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use MLSStringLists, only: ExtractSubstring, WriteIntsToList
  use MLSStrings, only: Lowercase
  use Output_M, only: Blanks, Output

  implicit none

  private

  public :: Bandwidth, Collapse, Depopulate, Repopulate
  public :: EmbedArray, EssentiallyEqual, ExtractArray
  public :: HyperTrace
  public :: Extremum
  public :: GatherBloc
  public :: HalfWaves
  public :: Rerank

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! This module contains stuff related to:
! (1) embedding/extracting blocs  or hyperslabs of elements
! (2) operations that reduce, expand, or compare arrays

! === (start of toc) ===                                                 
!     c o n t e n t s
!     - - - - - - - -

!         Functions, operations, routines
! Bandwidth         Calculate the Bandwidth of a banded array
! Collapse          Collapse an array into a lower rank array
!                    e.g., by summing over elements of dropped index
! Depopulate        Finds locations of non-zero values in presumably
!                     sparse array
! EmbedArray        Replace a hyperslab of elements in the larger array
!                     with the smaller
! EssentiallyEqual  Returns true if two real arguments 'close enough'
!                    (See comments below for interpretation
!                     of array versions)
! ExtractArray      Extract a hyperslab of elements from the larger array
!                     (optionally allocates slab first)
! Extremum          Returns the value farther from 0
! GatherBloc        Gather a subset of elements from a larger array
! HalfWaves         Count consecutive values with same sign
! HyperTrace        Show where a hyperslab of elements fit in the larger array
!                     as a string list of indexes
! Repopulate        Restores non-zero values in presumably
!                     sparse array to their proper locations
! === (end of toc) ===                                                   

! === (start of api) ===
! Bandwidth ( array, int bwidth, [real pctnzero], [int indxColpse] )
! Collapse ( arrayIn, nums, logs, [testvalue], [char* options] )
! Depopulate ( array[:], int i[:], int n, &
!   & [nprec testvalue], [nprec values[:], [char* options] )
! Depopulate ( array[:,:], int i[:], int j[:], int n, &
!   & [nprec testvalue], [nprec values[:], [char* options] )
! Depopulate ( array[:,:,:], int i[:], int j[:],  int k[:], int n, &
!   & [nprec testvalue], [nprec values[:], [char* options] )
! EmbedArray ( slab, array, &
!   int start[2], int count[2], int stride[2], int block[2], [char* options] )
! log EssentiallyEqual ( nprec A, nprec B, &
!   [nprec FillValue], [nprec Precision] )
! ExtractArray ( slab, array, &
!   int start[2], int count[2], int stride[2], int block[2], [char* options] )
! nprec Extremum( nprec arg1, nprec arg2 )
! GatherBloc ( bloc, array, int which1[:], &
!   [int which2(:)], [int which3(:), [char* options] )
! halfWaves( nprec array[:], int lengths[:], [int nWaves] )
! HyperTrace ( slab, &
!   int start[2], int count[2], int stride[2], int block[2], [char* options] )
! Repopulate ( array[:], int i[:], int n, &
!   & [nprec values[:], [char* options] )
! Repopulate ( array[:,:], int i[:], int j[:], int n, &
!   & [nprec values[:], [char* options] )
! Repopulate ( array[:,:,:], int i[:], int j[:],  int k[:], int n, &
!   & [nprec values[:], [char* options] )
!
! More Notes:
! (1) For the meaning of start, count, stride, and block, see
!     the wiki page https://mls.jpl.nasa.gov/team/wiki/index.php/Hyperslab
!     or the comments accompanying subroutine Gather in the l2/FillUtils module
! (2) In the case of GatherBloc of this module, which1, which2, ..
!     are the indices of components 1, 2, .. of the larger array to be
!     collected and returned in the smaller bloc
! === (end of api) ===                                                 
  interface BandWidth
    module procedure BandWidth_1dr4, BandWidth_1dr8, BandWidth_1dint
    module procedure BandWidth_2dr4, BandWidth_2dr8, BandWidth_2dint
    module procedure BandWidth_3dr4, BandWidth_3dr8, BandWidth_3dint
  end interface

  interface Collapse
    module procedure Collapse_2d_r4, Collapse_2d_r8, Collapse_2d_int
    module procedure Collapse_3d_r4, Collapse_3d_r8, Collapse_3d_int
  end interface

  interface Depopulate
    module procedure Depopulate_1d_r4, Depopulate_1d_r8
    module procedure Depopulate_2d_r4, Depopulate_2d_r8
    module procedure Depopulate_3d_r4, Depopulate_3d_r8
    module procedure Depopulate_1d_int
  end interface
  
  interface EmbedArray
    module procedure EmbedArray_1d_r4, EmbedArray_1d_r8
    module procedure EmbedArray_2d_r4, EmbedArray_2d_r8
    module procedure EmbedArray_3d_r4, EmbedArray_3d_r8
    module procedure EmbedArray_1d_int
  end interface

  interface HyperTrace
    module procedure HyperTrace_1d
    module procedure HyperTrace_2d
    module procedure HyperTrace_3d
  end interface

  interface EssentiallyEqual
    module procedure EssentiallyEqual_r4, EssentiallyEqual_r8
    module procedure EssentiallyEqual_r4_1d, EssentiallyEqual_r8_1d
    module procedure EssentiallyEqual_r4_2d, EssentiallyEqual_r8_2d
    module procedure EssentiallyEqual_r4_3d, EssentiallyEqual_r8_3d
  end interface

  interface ExtractArray
    module procedure ExtractArray_1d_r4, ExtractArray_1d_r8
    module procedure ExtractArray_2d_r4, ExtractArray_2d_r8
    module procedure ExtractArray_3d_r4, ExtractArray_3d_r8
    module procedure ExtractArray_1d_int
  end interface

  interface Extremum
    module procedure Extremum_REAL, Extremum_DOUBLE, Extremum_int
  end interface

  interface GatherBloc
    module procedure GatherBloc_1d_r4, GatherBloc_1d_r8
    module procedure GatherBloc_2d_r4, GatherBloc_2d_r8
    module procedure GatherBloc_3d_r4, GatherBloc_3d_r8
    module procedure GatherBloc_1d_int
  end interface

  interface halfWaves
    module procedure halfWaves_REAL, halfWaves_DOUBLE
  end interface

  interface Repopulate
    module procedure Repopulate_1d_r4, Repopulate_1d_r8
    module procedure Repopulate_2d_r4, Repopulate_2d_r8
    module procedure Repopulate_3d_r4, Repopulate_3d_r8
    module procedure Repopulate_1d_int
  end interface

  integer, parameter       :: MAXDEPOPULATED     = 1000000
  logical, parameter       :: DEEBUG             = .false.

contains

  ! ---------------------------------------------  Bandwidth  -----
  ! This family of routines calculates the Bandwidth of a presumably
  ! banded matrix
  ! Optionally it also returns the %age of elements that are non-zero
  ! Passed a higher-rank object, that object is first collapsed to a matrix
  ! Passed a rank-one object, it complains
  !
  ! Is there any advantage to be gained from generalizing to the case of
  ! a "banded" matrix where the elements outside the band are Fill values
  ! instead of zero?
  subroutine BandWidth_1dint( iarray, bwidth, pctnzero )
    ! Args
    integer, dimension(:), intent(in )    :: iarray
    integer, intent(out)                  :: bwidth
    real, optional, intent(out)           :: pctnzero
    call output( 'Sorry, Bandwidth inappropriate for 1-d array', advance='yes' )
    bwidth = 0
    if ( present(pctnzero) ) pctnzero = 0.
  end subroutine BandWidth_1dint

  subroutine BandWidth_1dr4( array, bwidth, pctnzero )
    ! Args
    integer, parameter :: RK = R4
    real(rk), dimension(:), intent(in )   :: array
    integer, intent(out)                  :: bwidth
    real, optional, intent(out)           :: pctnzero
    ! Internal variables
    call output( 'Sorry, Bandwidth inappropriate for 1-d array', advance='yes' )
    bwidth = 0
    if ( present(pctnzero) ) pctnzero = 0.
  end subroutine BandWidth_1dr4

  subroutine BandWidth_1dr8( array, bwidth, pctnzero )
    ! Args
    integer, parameter :: RK = R8
    real(rk), dimension(:), intent(in )   :: array
    integer, intent(out)                  :: bwidth
    real, optional, intent(out)           :: pctnzero
    ! Internal variables
    call output( 'Sorry, Bandwidth inappropriate for 1-d array', advance='yes' )
    bwidth = 0
    if ( present(pctnzero) ) pctnzero = 0.
  end subroutine BandWidth_1dr8

  subroutine BandWidth_2dint( iarray, bwidth, pctnzero )
    ! Args
    integer, dimension(:,:), intent(in )  :: iarray
    integer, intent(out)                  :: bwidth
    real, optional, intent(out)           :: pctnzero
    integer, parameter :: RK = R4
    ! Internal variables
    real(rk), dimension(size(iarray, 1), size(iarray, 2)) :: array
    array = iarray
    call BandWidth( array, bwidth, pctnzero )
  end subroutine BandWidth_2dint

  subroutine BandWidth_2dr4( array, bwidth, pctnzero )
    ! Args
    integer, parameter :: RK = R4
    real(rk), dimension(:,:), intent(in ) :: array
    integer, intent(out)                  :: bwidth
    real, optional, intent(out)           :: pctnzero
    ! Internal variables
    include 'BandWidth.f9h'
  end subroutine BandWidth_2dr4

  subroutine BandWidth_2dr8( array, bwidth, pctnzero )
    ! Args
    integer, parameter :: RK = R8
    real(rk), dimension(:,:), intent(in ) :: array
    integer, intent(out)                  :: bwidth
    real, optional, intent(out)           :: pctnzero
    ! Internal variables
    include 'BandWidth.f9h'
  end subroutine BandWidth_2dr8

  subroutine BandWidth_3dint( array, bwidth, pctnzero, indxColpse )
    ! Args
    integer, dimension(:,:,:), intent(in)    :: array
    integer, intent(out)                  :: bwidth
    real, optional, intent(out)           :: pctnzero
    integer, optional, intent(in)         :: indxColpse ! which index to collapse
    integer, dimension(:,:), pointer      :: ColpsedArray => null()
    include 'BandWidth3d.f9h'
  end subroutine BandWidth_3dint

  subroutine BandWidth_3dr4( array, bwidth, pctnzero, indxColpse )
    ! Args
    integer, parameter :: RK = R4
    real(rk), dimension(:,:,:), intent(in)   :: array
    integer, intent(out)                  :: bwidth
    real, optional, intent(out)           :: pctnzero
    integer, optional, intent(in)         :: indxColpse ! which index to collapse
    ! Internal variables
    real(rk), dimension(:,:), pointer     :: ColpsedArray => null()
    include 'BandWidth3d.f9h'
    ! call output( 'Sorry, Bandwidth inappropriate for 3-d array', advance='yes' )
  end subroutine BandWidth_3dr4

  subroutine BandWidth_3dr8( array, bwidth, pctnzero, indxColpse )
    ! Args
    integer, parameter :: RK = R8
    real(rk), dimension(:,:,:), intent(in)   :: array
    integer, intent(out)                  :: bwidth
    real, optional, intent(out)           :: pctnzero
    integer, optional, intent(in)         :: indxColpse ! which index to collapse
    real(rk), dimension(:,:), pointer     :: ColpsedArray => null()
    ! Internal variables
    include 'BandWidth3d.f9h'
    ! call output( 'Sorry, Bandwidth inappropriate for 3-d array', advance='yes' )
  end subroutine BandWidth_3dr8

  ! ---------------------------------------------  Collapse  -----
  ! This family of routines collapses an array[i1,i2,..in] to lower-rank
  ! arrays nums or logs[i1,i2,..,i(n-1)] by operating on all the elements of 
  ! the last index in depending upon options:
  ! options             nums                      logs
  ! 'num[op]'         use 'op'
  ! 'any[op]'                                 use 'any' with 'op'
  ! 'all[op]'                                 use 'all' with 'op' (default)
  
  !   op                          meaning
  !   '+'               sum                any or all even (default)
  !   '>'               max                any or all > testvalue
  !   '<'               min                any or all < testvalue
  !   '='                                  any or all = testvalue
  !   '!'                                  reverse sense of test
  !   ('!' may be added to any of ops above, e.g.: 'all[!+]' would
  !    return TRUE if all of the elements were odd)
  ! If testvalue is not supplied, it defaults to 0
  ! Options may contain separate fields for the nums and logs; e.g.,
  ! 'num[>]all[!=]' would return 
  !    nums    max values of elements over last index
  !    logs    true where all the elements were non-zero over last index
  ! 
  subroutine Collapse_2d_r4 ( array, nums, logs, testvalue, options )
    integer, parameter :: RK = R4
    include 'Collapse_2d.f9h'
  end subroutine Collapse_2d_r4

  subroutine Collapse_3d_r4 ( array, nums, logs, testvalue, options )
    integer, parameter :: RK = R4
    include 'Collapse_3d.f9h'
  end subroutine Collapse_3d_r4

  subroutine Collapse_2d_r8 ( array, nums, logs, testvalue, options )
    integer, parameter :: RK = r8
    include 'Collapse_2d.f9h'
  end subroutine Collapse_2d_r8

  subroutine Collapse_3d_r8 ( array, nums, logs, testvalue, options )
    integer, parameter :: RK = r8
    include 'Collapse_3d.f9h'
  end subroutine Collapse_3d_r8

  subroutine Collapse_2d_int ( iarray, inums, logs, itestvalue, options )
    integer, dimension(:,:), intent(in)           :: iarray ! The larger array
    integer , dimension(:), optional, intent(out) :: inums
    logical , dimension(:), optional, intent(out) :: logs
    integer, intent(in), optional                 :: itestvalue
    character(len=*), intent(in), optional        :: options
    ! Local variables
    real(r4), dimension(size(iarray, 1), size(iarray, 2)) :: array
    real(r4), dimension(size(iarray, 1))                  :: nums
    real(r4)                                              :: testvalue
    ! Executable
    array = iarray
    testvalue = 0
    if ( present(itestvalue) ) testvalue = itestvalue
    call Collapse( array, nums, logs, testvalue, options )
    if ( present(inums) ) inums = nums
  end subroutine Collapse_2d_int

  subroutine Collapse_3d_int ( iarray, inums, logs, itestvalue, options )
    integer, dimension(:,:,:), intent(in)           :: iarray ! The larger array
    integer , dimension(:,:), optional, intent(out) :: inums
    logical , dimension(:,:), optional, intent(out) :: logs
    integer, intent(in), optional                   :: itestvalue
    character(len=*), intent(in), optional          :: options
    ! Local variables
    real(r4), dimension(size(iarray, 1), size(iarray, 2), size(iarray, 3)) &
      &                                                   :: array
    real(r4), dimension(size(iarray, 1), size(iarray, 2)) :: nums
    real(r4)                                              :: testvalue
    ! Executable
    array = iarray
    testvalue = 0
    if ( present(itestvalue) ) testvalue = itestvalue
    call Collapse( array, nums, logs, testvalue, options )
    if ( present(inums) ) inums = nums
  end subroutine Collapse_3d_int

  ! ---------------------------------------------  Depopulate  -----
  ! This family of routines returns only non-zero values of an array[i,j,..]
  ! or, depending upon options:
  ! options             which to return
  !   '+'              values  >  testvalue
  !   '/'              values  /= testvalue
  !   '%'              compare values with testvalue*maxval(values)
  !  '||'              test according to |array[i]|
  !  'r '              reverse sense; i.e. return just the "zeros"
  ! If testvalue is not supplied, it defaults to 0
  ! Optionally we return just the indices i for 1-d arrays, (i,j) for 2-d,
  ! (i,j,k) for 3-d, etc.
  ! This is similar in spirit to "sparifying" a full block Hessian
  ! See also Repopulate
  subroutine Depopulate_1d_int ( iarray, i, n, itestvalue, ivalues, options )
    integer, dimension(:), intent(out) :: i
    integer, dimension(:), intent(in) :: iarray ! The larger array
    integer, intent(out) :: n
    integer, parameter :: RK = R4
    integer, intent(in), optional     :: itestvalue
    integer, dimension(:), intent(out), optional :: ivalues
    character(len=*), intent(in), optional :: options
    ! Local variables
    real(rk), dimension(size(iarray)) :: array ! The sparse array
    real(rk) :: testvalue
    real(rk), dimension(MAXDepopulateD) :: values
    ! Executable
    testvalue = 0._rk
    if ( present(itestvalue) ) testvalue = real( itestvalue, rk )
    array = iarray
    call Depopulate ( array, i, n, testvalue, values, options )
    if ( present(ivalues) ) ivalues = values(1:size(ivalues))
  end subroutine Depopulate_1d_int

  subroutine Depopulate_1d_r4 ( array, i, n, testvalue, values, options )
    integer, parameter :: RK = R4
    include 'Decimate_1d.f9h'
  end subroutine Depopulate_1d_r4

  subroutine Depopulate_2d_r4 ( array, i, j, n, testvalue, values, options )
    integer, parameter :: RK = R4
    include 'Decimate_2d.f9h'
  end subroutine Depopulate_2d_r4

  subroutine Depopulate_3d_r4 ( array, i, j, k, n, testvalue, values, options )
    integer, parameter :: RK = R4
    include 'Decimate_3d.f9h'
  end subroutine Depopulate_3d_r4

  subroutine Depopulate_1d_r8 ( array, i, n, testvalue, values, options )
    integer, parameter :: RK = r8
    include 'Decimate_1d.f9h'
  end subroutine Depopulate_1d_r8

  subroutine Depopulate_2d_r8 ( array, i, j, n, testvalue, values, options )
    integer, parameter :: RK = r8
    include 'Decimate_2d.f9h'
  end subroutine Depopulate_2d_r8

  subroutine Depopulate_3d_r8 ( array, i, j, k, n, testvalue, values, options )
    integer, parameter :: RK = r8
    include 'Decimate_3d.f9h'
  end subroutine Depopulate_3d_r8

  ! ---------------------------------------------  EmbedArray  -----
  ! This family of routines replace a bloc of elements in a larger
  ! array with corresponding elements from a smaller
  subroutine EmbedArray_1d_int ( ibloc, iarray, start, count, stride, block, options )
    integer, dimension(:), pointer :: ibloc
    integer, dimension(:), pointer :: iarray ! The larger array
    integer, parameter :: RK = R4
    integer, dimension(:), intent(in)     :: start
    integer, dimension(:), intent(in)     :: count
    integer, dimension(:), intent(in)     :: stride
    integer, dimension(:), intent(in)     :: block
    character(len=*), intent(in), optional :: options
    ! Local variables
    real(rk), dimension(:), pointer :: bloc
    real(rk), dimension(:), pointer :: array ! The larger array
    ! Executable
    nullify ( array, bloc )
    call allocate_test ( bloc, size(ibloc), ModuleName, &
      & "1-d bloc for int embedding" )
    call allocate_test ( array, size(iarray), ModuleName, &
      & "1-d array for int embedding" )
    array = iarray
    bloc = ibloc
    call EmbedArray ( bloc, array, start, count, stride, block, options )
    iarray = array
    call deallocate_test ( array, ModuleName, "1-d array for int embedding" )
    call deallocate_test ( bloc, ModuleName, "1-d bloc for int embedding" )
  end subroutine EmbedArray_1d_int

  subroutine EmbedArray_1d_r4 ( slab, array, start, count, stride, block, options )
    integer, parameter :: RK = R4
    character, parameter :: DIRECTION = 'm'
    include 'EmbedExtract_1d.f9h'
  end subroutine EmbedArray_1d_r4

  subroutine EmbedArray_1d_r8 ( slab, array, start, count, stride, block, options )
    integer, parameter :: RK = R8
    character, parameter :: DIRECTION = 'm'
    include 'EmbedExtract_1d.f9h'
  end subroutine EmbedArray_1d_r8

  subroutine ExtractArray_1d_int ( ibloc, iarray, start, count, stride, block, options )
    integer, dimension(:), pointer :: ibloc
    integer, dimension(:) :: iarray ! The larger array
    integer, parameter :: RK = R4
    integer, dimension(:), intent(in)     :: start
    integer, dimension(:), intent(in)     :: count
    integer, dimension(:), intent(in)     :: stride
    integer, dimension(:), intent(in)     :: block
    character(len=*), intent(in), optional :: options
    ! Local variables
    real(rk), dimension(:), pointer :: bloc
    real(rk), dimension(:), pointer :: array ! The larger array
    ! Executable
    nullify ( array, bloc )
    if ( associated(ibloc) ) then
      call allocate_test ( bloc, size(ibloc), ModuleName, &
        & "1-d bloc for int extracting" )
    end if
    call allocate_test ( array, size(iarray), ModuleName, &
        & "1-d array for int extracting" )
    array = iarray
    call ExtractArray ( bloc, array, start, count, stride, block, options )
    if ( .not. associated(ibloc) ) then
      call allocate_test ( ibloc, size(bloc), ModuleName, &
        & "1-d bloc for int extracting" )
    end if
    ibloc = bloc
    call deallocate_test ( array, ModuleName, "1-d array for int extracting" )
    call deallocate_test ( bloc, ModuleName, "1-d bloc for int extracting" )
  end subroutine ExtractArray_1d_int

  subroutine ExtractArray_1d_r4 ( slab, array, start, count, stride, block, options )
    integer, parameter :: RK = R4
    character, parameter :: DIRECTION = 'x'
    include 'EmbedExtract_1d.f9h'
  end subroutine ExtractArray_1d_r4

  subroutine ExtractArray_1d_r8 ( slab, array, start, count, stride, block, options )
    integer, parameter :: RK = R8
    character, parameter :: DIRECTION = 'x'
    include 'EmbedExtract_1d.f9h'
  end subroutine ExtractArray_1d_r8

  subroutine EmbedArray_2d_r4 ( slab, array, start, count, stride, block, options )
    integer, parameter :: RK = R4
    character, parameter :: DIRECTION = 'm'
    include 'EmbedExtract_2d.f9h'
  end subroutine EmbedArray_2d_r4

  subroutine EmbedArray_2d_r8 ( slab, array, start, count, stride, block, options )
    integer, parameter :: RK = R8
    character, parameter :: DIRECTION = 'm'
    include 'EmbedExtract_2d.f9h'
  end subroutine EmbedArray_2d_r8

  subroutine ExtractArray_2d_r4 ( slab, array, start, count, stride, block, options )
    integer, parameter :: RK = R4
    character, parameter :: DIRECTION = 'x'
    include 'EmbedExtract_2d.f9h'
  end subroutine ExtractArray_2d_r4

  subroutine ExtractArray_2d_r8 ( slab, array, start, count, stride, block, options )
    integer, parameter :: RK = R8
    character, parameter :: DIRECTION = 'x'
    include 'EmbedExtract_2d.f9h'
  end subroutine ExtractArray_2d_r8

  subroutine EmbedArray_3d_r4 ( slab, array, start, count, stride, block, options )
    integer, parameter :: RK = R4
    character, parameter :: DIRECTION = 'm'
    include 'EmbedExtract_3d.f9h'
  end subroutine EmbedArray_3d_r4

  subroutine EmbedArray_3d_r8 ( slab, array, start, count, stride, block, options )
    integer, parameter :: RK = R8
    character, parameter :: DIRECTION = 'm'
    include 'EmbedExtract_3d.f9h'
  end subroutine EmbedArray_3d_r8

  subroutine ExtractArray_3d_r4 ( slab, array, start, count, stride, block, options )
    integer, parameter :: RK = R4
    character, parameter :: DIRECTION = 'x'
    include 'EmbedExtract_3d.f9h'
  end subroutine ExtractArray_3d_r4

  subroutine ExtractArray_3d_r8 ( slab, array, start, count, stride, block, options )
    integer, parameter :: RK = R8
    character, parameter :: DIRECTION = 'x'
    include 'EmbedExtract_3d.f9h'
  end subroutine ExtractArray_3d_r8

  ! ---------------------------------------------  HyperTrace  -----
  ! Write as a string list the indexes inside the larger array
  ! corresponding to each element in the hyperslab
  ! if performing an Embed or Extract operation.
  !
  ! E.g.,
  ! Say the larger array is rank 3, with shape (5,5,4)
  ! (although those particular values are not relevant in what follows)
  !
  ! Now if a hyperlab of array is selected by setting
  ! 
  ! start  = (/ 0,0,0 /)
  ! count  = (/ 2,2,2 /)
  ! stride = (/ 2,2,1 /)
  ! block  = (/ 1,1,1 /)
  !
  ! then slab would be returned as the 2x2x2 3d matrix of stringlists
  !     1   1   1# 1,1,1             1,1,2                        
  !     1   2   1# 1,3,1             1,3,2                       
  !     2   1   1# 3,1,1             3,1,2            
  !     2   2   1# 3,3,1             3,3,2            
  ! (which as been flattened by dump into the 2 2x2 matrices shown above)
  
  subroutine HyperTrace_1d ( slab, start, count, stride, block )
    character(len=*), dimension(:)        :: slab
    integer, dimension(:), intent(in)     :: start
    integer, dimension(:), intent(in)     :: count
    integer, dimension(:), intent(in)     :: stride
    integer, dimension(:), intent(in)     :: block
    ! Local variables
    integer                                :: II
    integer                                :: k
    integer                                :: n1
    ! integer                                :: n1b
    integer                                :: m
    integer                                :: m1
    ! integer                                :: m1b
    ! Executable
    do II=1, count(1)
      n1 = start(1) - MLS_HyperStart + (II-1)*stride(1) + 1
      ! n1b = start(1) - MLS_HyperStart + (II-1)*stride(1) + block(1)
      m1 = (II-1)*block(1) + 1
      ! m1b = (II-1)*block(1) + block(1)
      do m=0, block(1) - 1
        call WriteIntsToList ( &
          & (/ n1 + m /), &
          & slab( m1 + m ) )
      enddo
    enddo
  end subroutine HyperTrace_1d

  subroutine HyperTrace_2d ( slab, start, count, stride, block )
    character(len=*), dimension(:,:)      :: slab
    integer, dimension(:), intent(in)     :: start
    integer, dimension(:), intent(in)     :: count
    integer, dimension(:), intent(in)     :: stride
    integer, dimension(:), intent(in)     :: block
    ! Local variables
    integer                                :: II
    integer                                :: JJ
    integer                                :: j
    integer                                :: k
    integer                                :: n1
    ! integer                                :: n1b
    integer                                :: n2
    integer                                :: m
    integer                                :: m1
    ! integer                                :: m1b
    integer                                :: m2
    ! Executable
    do JJ=1, count(2)
      do II=1, count(1)
        do j=1, block(2)
          n2 = start(2) - MLS_HyperStart + (JJ-1)*stride(2) + j
          m2 = (JJ-1)*block(2) + j
          n1 = start(1) - MLS_HyperStart + (II-1)*stride(1) + 1
          ! n1b = start(1) - MLS_HyperStart + (II-1)*stride(1) + block(1)
          m1 = (II-1)*block(1) + 1
          ! m1b = (II-1)*block(1) + block(1)
          do k=0, block(1) - 1
            call WriteIntsToList ( &
              & (/ n1 + k, n2 /), &
              & slab(m1 + k, m2 ) )
          enddo
        enddo
      enddo
    enddo
  end subroutine HyperTrace_2d

  subroutine HyperTrace_3d ( slab, start, count, stride, block )
    character(len=*), dimension(:,:,:)    :: slab
    integer, dimension(:), intent(in)     :: start
    integer, dimension(:), intent(in)     :: count
    integer, dimension(:), intent(in)     :: stride
    integer, dimension(:), intent(in)     :: block
    ! Local variables
    integer                                :: II
    integer                                :: JJ
    integer                                :: KK
    integer                                :: j
    integer                                :: k
    integer                                :: p
    integer                                :: n1
    ! integer                                :: n1b
    integer                                :: n2
    integer                                :: n3
    integer                                :: m
    integer                                :: m1
    ! integer                                :: m1b
    integer                                :: m2
    integer                                :: m3
    ! Executable
    do KK=1, count(3)
      do JJ=1, count(2)
        do II=1, count(1)
          do k=1, block(3)
            n3 = start(3) - MLS_HyperStart + (KK-1)*stride(3) + k
            m3 = (KK-1)*block(3) + k
            do j=1, block(2)
              n2 = start(2) - MLS_HyperStart + (JJ-1)*stride(2) + j
              m2 = (JJ-1)*block(2) + j
              n1 = start(1) - MLS_HyperStart + (II-1)*stride(1) + 1
              ! n1b = start(1) - MLS_HyperStart + (II-1)*stride(1) + block(1)
              m1 = (II-1)*block(1) + 1
              ! m1b = (II-1)*block(1) + block(1)
              ! slab(m1:m1b, m2, m3) = array(n1:n1b, n2, n3)
              do p=0, block(1) - 1
                call WriteIntsToList ( &
                  & (/ n1 + p, n2, n3 /), &
                  & slab(m1 + p, m2, m3 ) )
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine HyperTrace_3d

! ---------------------------------------------  EssentiallyEqual  -----

  ! This family of routines checks to see if two reals are essentially
  ! the same.
  elemental logical function EssentiallyEqual_r4 ( A, B )
    real(r4), intent(in) :: A
    real(r4) ,intent(in) :: B
    EssentiallyEqual_r4 = &
      & a >= nearest ( b, -1.0_r4 ) .and. a <= nearest ( b, 1.0_r4 )
  end function EssentiallyEqual_r4

  elemental logical function EssentiallyEqual_r8 ( A, B )
    real(r8), intent(in) :: A
    real(r8) ,intent(in) :: B
    EssentiallyEqual_r8 = &
      & a >= nearest ( b, -1.0_r8 ) .and. a <= nearest ( b, 1.0_r8 )
  end function EssentiallyEqual_r8

  ! The following functions are slightly different:
  ! A scalar logical testing that every element of two arrays are equal
  ! You may filter out values in either array equal to the FillValue
  ! or for which the corresponding Precision array is negative
  ! or where either element is not finite
  ! Warn if an element of one array is finite while the other is not
  function EssentiallyEqual_r4_1d ( A, B, FillValue, Precision ) &
    & result(equal)
    real(r4), dimension(:), intent(in)             :: A
    real(r4), dimension(:), intent(in)             :: B
    real(r4), intent(in)                           :: fillValue
    real(r4), dimension(:), optional, intent(in)   :: precision
    logical                                        :: equal
    real(r4), dimension(size(A))                   :: atab
    real(r4), dimension(size(B))                   :: btab
    logical                                        :: warn
    equal = .false.
    call filterValues(A, ATAB, B, BTAB, warn, fillValue, precision)
    if ( .not. warn ) equal = all( &
      & a >= nearest ( b, -1.0_r4 ) .and. a <= nearest ( b, 1.0_r4 ) &
      & )
  end function EssentiallyEqual_r4_1d

  function EssentiallyEqual_r8_1d ( A, B, FillValue, Precision ) &
    & result(equal)
    real(r8), dimension(:), intent(in)             :: A
    real(r8), dimension(:), intent(in)             :: B
    real(r8), intent(in)                           :: fillValue
    real(r8), dimension(:), optional, intent(in)   :: precision
    logical                                        :: equal
    real(r8), dimension(size(A))                   :: atab
    real(r8), dimension(size(B))                   :: btab
    logical                                        :: warn
    equal = .false.
    call filterValues(A, ATAB, B, BTAB, warn, fillValue, precision)
    if ( .not. warn ) equal = all( &
      & a >= nearest ( b, -1.0_r8 ) .and. a <= nearest ( b, 1.0_r8 ) &
      & )
  end function EssentiallyEqual_r8_1d

  function EssentiallyEqual_r4_2d ( A, B, FillValue, Precision ) &
    & result(equal)
    real(r4), dimension(:,:), intent(in)             :: A
    real(r4), dimension(:,:), intent(in)             :: B
    real(r4), intent(in)                             :: fillValue
    real(r4), dimension(:,:), optional, intent(in)   :: precision
    logical                                        :: equal
    ! Internal variables
    integer, dimension(2)                          :: shp
    shp =shape(a)
    if ( present(Precision) ) then
      equal = EssentiallyEqual_r4_1d(reshape(a, (/shp(1)*shp(2)/)), &
        & reshape(b, (/shp(1)*shp(2)/)), &
        & FillValue, reshape(Precision, (/shp(1)*shp(2)/)) )
    else
      equal = EssentiallyEqual_r4_1d(reshape(a, (/shp(1)*shp(2)/)), &
        & reshape(b, (/shp(1)*shp(2)/)), FillValue )
    endif
      
  end function EssentiallyEqual_r4_2d

  function EssentiallyEqual_r8_2d ( A, B, FillValue, Precision ) &
    & result(equal)
    real(r8), dimension(:,:), intent(in)             :: A
    real(r8), dimension(:,:), intent(in)             :: B
    real(r8), intent(in)                             :: fillValue
    real(r8), dimension(:,:), optional, intent(in)   :: precision
    logical                                        :: equal
    ! Internal variables
    integer, dimension(2)                          :: shp
    shp =shape(a)
    if ( present(Precision) ) then
      equal = EssentiallyEqual_r8_1d(reshape(a, (/shp(1)*shp(2)/)), &
        & reshape(b, (/shp(1)*shp(2)/)), &
        & FillValue, reshape(Precision, (/shp(1)*shp(2)/)) )
    else
      equal = EssentiallyEqual_r8_1d(reshape(a, (/shp(1)*shp(2)/)), &
        & reshape(b, (/shp(1)*shp(2)/)), FillValue )
    endif
      
  end function EssentiallyEqual_r8_2d

  function EssentiallyEqual_r4_3d ( A, B, FillValue, Precision ) &
    & result(equal)
    real(r4), dimension(:,:,:), intent(in)             :: A
    real(r4), dimension(:,:,:), intent(in)             :: B
    real(r4), intent(in)                               :: fillValue
    real(r4), dimension(:,:,:), optional, intent(in)   :: precision
    logical                                        :: equal
    ! Internal variables
    integer, dimension(3)                          :: shp
    shp =shape(a)
    if ( present(Precision) ) then
      equal = EssentiallyEqual_r4_1d(reshape(a, (/shp(1)*shp(2)*shp(3)/)), &
        & reshape(b, (/shp(1)*shp(2)*shp(3)/)), &
        & FillValue, reshape(Precision, (/shp(1)*shp(2)*shp(3)/)) )
    else
      equal = EssentiallyEqual_r4_1d(reshape(a, (/shp(1)*shp(2)*shp(3)/)), &
        & reshape(b, (/shp(1)*shp(2)*shp(3)/)), FillValue )
    endif
      
  end function EssentiallyEqual_r4_3d

  function EssentiallyEqual_r8_3d ( A, B, FillValue, Precision ) &
    & result(equal)
    real(r8), dimension(:,:,:), intent(in)             :: A
    real(r8), dimension(:,:,:), intent(in)             :: B
    real(r8), intent(in)                               :: fillValue
    real(r8), dimension(:,:,:), optional, intent(in)   :: precision
    logical                                        :: equal
    ! Internal variables
    integer, dimension(3)                          :: shp
    shp =shape(a)
    if ( present(Precision) ) then
      equal = EssentiallyEqual_r8_1d(reshape(a, (/shp(1)*shp(2)*shp(3)/)), &
        & reshape(b, (/shp(1)*shp(2)*shp(3)/)), &
        & FillValue, reshape(Precision, (/shp(1)*shp(2)*shp(3)/)) )
    else
      equal = EssentiallyEqual_r8_1d(reshape(a, (/shp(1)*shp(2)*shp(3)/)), &
        & reshape(b, (/shp(1)*shp(2)*shp(3)/)), FillValue )
    endif
      
  end function EssentiallyEqual_r8_3d

! -------------------------------------------------  Extremum  -----
! Returns the arg with the large absolute value
  elemental function Extremum_DOUBLE(arg1, arg2) result(value)
    double precision, intent(in) :: arg1, arg2
    double precision             :: value
    if ( abs(arg2) > abs(arg1) ) then
      value = arg2
    else
      value = arg1
    endif
  end function Extremum_DOUBLE

  elemental function Extremum_REAL(arg1, arg2) result(value)
    real, intent(in)        :: arg1, arg2
    real                    :: value
    if ( abs(arg2) > abs(arg1) ) then
      value = arg2
    else
      value = arg1
    endif
  end function Extremum_REAL

  elemental function Extremum_int(arg1, arg2) result(value)
    integer, intent(in)        :: arg1, arg2
    integer                    :: value
    if ( abs(arg2) > abs(arg1) ) then
      value = arg2
    else
      value = arg1
    endif
  end function Extremum_int

! -------------------------------------------------  GatherBloc  -----
! This family of subroutines gathers {array[i], i in which1} into bloc}
! if 'a' is among options, we allocate bloc
! We could also code ScatterArray, but choose not to do so
  subroutine GatherBloc_1d_int ( ibloc, iarray, which1, options )
    integer, dimension(:), pointer :: ibloc
    integer, dimension(:), pointer :: iarray ! The larger array
    integer, parameter :: RK = R4
    integer, dimension(:), intent(in)     :: which1
    character(len=*), intent(in), optional :: options
    ! Local variables
    real(rk), dimension(:), pointer :: bloc => null()
    real(rk), dimension(:), pointer :: array => null() ! The larger array
    integer :: status
    ! Executable
    if ( associated(ibloc) ) then
      allocate( bloc(size(ibloc)), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, "unable to allocate 1-d bloc for int gathering" )
    endif
    allocate( array(size(iarray)), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, "unable to allocate 1-d array for int gathering" )
    array = iarray
    call GatherBloc ( bloc, array, which1, options )
    if ( .not. associated(ibloc) ) then
      allocate( ibloc(size(bloc)), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, "unable to allocate 1-d bloc for int extracting" )
    endif
    ibloc = bloc
    deallocate( bloc, array, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, "unable to deallocate 1-d bloc for int extracting" )
  end subroutine GatherBloc_1d_int

  subroutine GatherBloc_1d_r4 ( bloc, array, which1, options )
    integer, parameter :: RK = r4
    character, parameter :: DIRECTION = 'g'
    include 'GatherScatter_1d.f9h'
  end subroutine GatherBloc_1d_r4

  subroutine GatherBloc_2d_r4 ( bloc, array, which1, which2, options )
    integer, parameter :: RK = r4
    character, parameter :: DIRECTION = 'g'
    include 'GatherScatter_2d.f9h'
  end subroutine GatherBloc_2d_r4

  subroutine GatherBloc_3d_r4 ( bloc, array, which1, which2, which3, options )
    integer, parameter :: RK = r4
    character, parameter :: DIRECTION = 'g'
    include 'GatherScatter_3d.f9h'
  end subroutine GatherBloc_3d_r4

  subroutine GatherBloc_1d_r8 ( bloc, array, which1, options )
    integer, parameter :: RK = R8
    character, parameter :: DIRECTION = 'g'
    include 'GatherScatter_1d.f9h'
  end subroutine GatherBloc_1d_r8

  subroutine GatherBloc_2d_r8 ( bloc, array, which1, which2, options )
    integer, parameter :: RK = R8
    character, parameter :: DIRECTION = 'g'
    include 'GatherScatter_2d.f9h'
  end subroutine GatherBloc_2d_r8

  subroutine GatherBloc_3d_r8 ( bloc, array, which1, which2, which3, options )
    integer, parameter :: RK = R8
    character, parameter :: DIRECTION = 'g'
    include 'GatherScatter_3d.f9h'
  end subroutine GatherBloc_3d_r8

! -------------------------------------------------  HalfWaves  -----
! Calculates the half-wave lengths of an array
! i.e., the number of consecutive pos. or neg. values
  subroutine halfWaves_DOUBLE( array, lengths, nWaves )
    real(r8), dimension(:), intent(in)         :: array
    integer, dimension(:), intent(out)         :: lengths
    integer, optional, intent(out)             :: nWaves
    integer, parameter :: RK = R8
    include 'halfWaves.f9h'
  end subroutine halfWaves_DOUBLE

  subroutine halfWaves_REAL( array, lengths, nWaves )
    real(r4), dimension(:), intent(in)         :: array
    integer, dimension(:), intent(out)         :: lengths
    integer, optional, intent(out)             :: nWaves
    integer, parameter :: RK = R4
    include 'halfWaves.f9h'
  end subroutine halfWaves_REAL

  ! ---------------------------------------------  Repopulate  -----
  ! This family of routines restores non-zero values of a presumably
  ! sparse array to their proper locations, depending on options
  ! options             does which 
  !   ' '              replaces value at (i, j, ..)
  !   '+'              added to value at (i, j, ..)
  !   '-'              subtracted from value at (i, j, ..)
  !   '>'              replaces value at (i, j, ..) only if bigger
  !   '<'              replaces value at (i, j, ..) only if smaller
  ! See Depopulate
  subroutine Repopulate_1d_int ( iarray, i, n, ivalues, options )
    integer, dimension(:), intent(in) :: i
    integer, dimension(:), intent(inout) :: iarray ! The larger array
    integer, intent(in) :: n
    integer, parameter :: RK = R4
    integer, dimension(:), intent(in) :: ivalues
    character(len=*), intent(in), optional :: options
    ! Local variables
    real(rk), dimension(size(iarray)) :: array ! The sparse array
    real(rk), dimension(size(ivalues)) :: values
    ! Executable
    values = ivalues
    array = iarray
    call Repopulate ( array, i, n, values, options )
    iarray = array
  end subroutine Repopulate_1d_int

  subroutine Repopulate_1d_r4 ( array, i, n, values, options )
    integer, parameter :: RK = R4
    include 'Repopulate_1d.f9h'
  end subroutine Repopulate_1d_r4

  subroutine Repopulate_2d_r4 ( array, i, j, n, values, options )
    integer, parameter :: RK = R4
    include 'Repopulate_2d.f9h'
  end subroutine Repopulate_2d_r4

  subroutine Repopulate_3d_r4 ( array, i, j, k, n, values, options )
    integer, parameter :: RK = R4
    include 'Repopulate_3d.f9h'
  end subroutine Repopulate_3d_r4

  subroutine Repopulate_1d_r8 ( array, i, n, values, options )
    integer, parameter :: RK = r8
    include 'Repopulate_1d.f9h'
  end subroutine Repopulate_1d_r8

  subroutine Repopulate_2d_r8 ( array, i, j, n, values, options )
    integer, parameter :: RK = r8
    include 'Repopulate_2d.f9h'
  end subroutine Repopulate_2d_r8

  subroutine Repopulate_3d_r8 ( array, i, j, k, n, values, options )
    integer, parameter :: RK = r8
    include 'Repopulate_3d.f9h'
  end subroutine Repopulate_3d_r8

  ! ----------------------------------  private procedures  -----
  subroutine announce_error(message, int1, int2, dontstop)
    character(len=*), intent(in) :: message
    integer, optional, intent(in) :: int1
    integer, optional, intent(in) :: int2
    logical, optional, intent(in) :: dontstop
    logical :: keepgoing
    !
    keepgoing = .false.
    if ( present(dontstop) ) keepgoing=dontstop
    if ( .not. keepgoing ) then
      call output('*** Error in MLSFillValues module ***', advance='yes')
    endif
    call output(trim(message), advance='no')
    call blanks(3)
    if ( present(int1) ) write(*,'(i4)',advance='no') int1
    call blanks(3)
    if ( present(int2) ) write(*,'(i4)', advance='no') int2
    if ( .not. keepgoing ) stop
  end subroutine announce_error

  ! Find multidimensional set of indices in an array
  ! with shape shp corresponding to 1-d address
  !
  ! We shall assume that the first index is the fastest, then the 2nd, ..
  ! Our method is the following:
  ! Let the size of the kth index be s[k]
  ! Then we seek the array i[k] such that
  ! address = i[1] + s[1] ( i[2] + s[2] ( i[3] + .. + i[N] ) .. )
  ! (Where we assume 0-based indexing, like c, 
  !   rather than 1-based, as Fortran uses)
  ! We can build this by parts as follows
  ! a[N]   = i[N]
  ! a[N-1] = i[N-1] + s[N-1] i[N]
  ! a[N-2] = i[N-2] + s[N-2] ( i[N-1] + s[N-1] i[N] )
  ! .   .   .
  ! a[1] = address
  ! Where N is the rank of the array
  ! Note then that the recurrences hold
  ! a[N-1] - a[N] s[N-1]   = i[N-1]
  ! a[N-2] - a[N-1] s[N-2] = i[N-2]
  ! .   .   .
  ! a[1] - a[2] s[1]       = i[1]
  !
  ! From this last we realize that
  ! i[1] = a[1] mod(s[1])
  ! Solve it for i[1], then a[2] = ( a[1] - i[1] ) / s[1]
  ! Then for succeeding values of k
  ! i[k] = a[k] mod(s[k])

  ! Remember to modify each of these if we wish to use
  ! Fortran-style indexes which start at 1, not 0, as follows
  !
  ! i'[k] = i[k] + 1, k > 1
  ! i'[1] = i[1]
  ! address = i'[1] + s[1] ( i'[2] - 1 + s[2] ( i[3] - 1 + .. + i[N] ) .. )
  ! Actually, let's just use Subscripts from Array_Stuff
  subroutine rerank( address, shp, indices )
    use Array_Stuff, only: Subscripts
    integer, intent(in)                :: address
    integer, dimension(:), intent(in)  :: shp
    integer, dimension(:), intent(out) :: indices
    indices = Subscripts( address, shp )
  end subroutine rerank

  subroutine rerank_old( address, shp, indices )
    integer, intent(in)                :: address
    integer, dimension(:), intent(in)  :: shp
    integer, dimension(:), intent(out) :: indices
    ! Local variables
    integer :: aofk
    integer :: k
    integer :: N
    integer, parameter :: OFFSET = 1 ! at what index do arrays start?
    !
    N = size(shp)
    if ( N < 2 ) then
      indices(1) = address
      return
    end if
    aofk = address - OFFSET
    do k=1, N
      indices(k) = MOD(aofk, shp(k))
      aofk = ( aofk - indices(k) ) / shp(k)
    enddo
    indices = indices + OFFSET ! Converting to Fortran-style, beginning with 1
  end subroutine rerank_old
  
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module HyperSlabs
!=============================================================================

!
! $Log$
! Revision 2.3  2021/04/29 22:52:55  pwagner
! Added HyperTrace
!
! Revision 2.2  2021/04/15 22:44:27  pwagner
! Now uses MLS_HyperStart
!
! Revision 2.1  2017/11/03 19:55:09  pwagner
! First commit
!
