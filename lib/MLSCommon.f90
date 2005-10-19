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
module MLSCommon                ! Common definitions for the MLS software
!=============================================================================

  use ieee_arithmetic, only: ieee_is_finite
  use MLSKinds ! Everything
  use SDPTOOLKIT, only: PGSd_PC_FILE_PATH_MAX

  implicit none
  private

! === (start of toc) ===                                                 
!     c o n t e n t s                                                    
!     - - - - - - - -                                                    

!     (data types and parameters)

! i1, i2, i4    integer types
! r4, r8        floating point types
! ip, rp        integer, floating point types used in forward model
! rv            floating point type used in vector quantity values
! rm            floating point type used in matrix values
! NameLen       character-length of quantity names
! LineLen       character-length of most input
! FileNameLen   character-length of path/filenames
! BareFNLen     character-length of filenames
! L1BInfo_T     L1B data file names, etc.

!     (subroutines and functions)
! === (end of toc) ===                                                   
! === (start of api) ===
! === (end of api) ===

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! This module contains simple definitions that are common to all the MLS PGS
  ! f90 software.

  public :: FileIDs_T
  public :: FilterValues
  public :: InRange
  public :: IsFinite
  public :: MLSFile_T
  public :: L1BInfo_T
  public :: Range_T
  public :: TAI93_Range_T

  ! Make parameters gotten from MLSKinds public

  public :: i1
  public :: i2
  public :: i4
  public :: r4
  public :: r8
  public :: rm
  public :: rp
  public :: ip
  public :: rv

  ! Now we have the lengths for various strings

  integer, public, parameter :: NameLen=32
  integer, public, parameter :: LineLen=132
  integer, public, parameter :: FileNameLen=max(PGSd_PC_FILE_PATH_MAX, 132) ! was 132
  integer, public, parameter :: BareFNLen=64      ! Bare file name length (w/o path)

  real, public, parameter ::    DEFAULTUNDEFINEDVALUE = -999.99 ! Try to use in lib, l2
  real, parameter, private :: FILLVALUETOLERANCE = 0.2 ! Poss. could make it 1
  ! --------------------------------------------------------------------------
  
  ! A type to hold the hdf ids

  type FileIds_T
    integer :: f_id     = 0 ! File id, handle, or io unit
    integer :: grp_id   = 0 ! group id
    integer :: sd_id    = 0 ! sd or swath id
  end type Fileids_T
  ! --------------------------------------------------------------------------

  ! A PCFid range

  type Range_T
    integer :: Bottom   = 0
    integer :: Top      = 0
  end type Range_T
  ! --------------------------------------------------------------------------

  ! Moved here from MLSFiles module
  ! Information describing the files used by the mls software
  ! Stop passing file handles back & forth between routines
  ! -- pass one of these instead
  ! (Not used yet; maybe someday)
  type MLSFile_T
    character (LEN=16) :: content=""  ! e.g., 'l1brad', 'l2gp', 'l2aux', ..
    character (LEN=8) :: lastOperation=""  ! 'open','close','read','write'
    character (LEN=FileNameLen) :: Name=""  ! its name (usu. w/path)
    character (LEN=NameLen) :: ShortName=""  ! its short name; e.g. 'H2O'
    character (LEN=8) :: typeStr=""  ! one of {'swath', 'hdf', ..}
    integer :: type=0  ! one of {l_swath, l_hdf, ..}
    integer :: access=0  ! one of {DFACC_RDONLY, DFACC_CREATE, ..}
    integer :: HDFVersion=0  ! its hdf version if hdf(eos)
    integer :: PCFId=0      ! its PCF ID (ref), if any
    integer :: recordLength=0! its max record_length, if any
    integer :: errorCode=0  ! non-zero usu. means trouble
    logical :: StillOpen=.false.
    type(Range_T) :: PCFidRange
    type(Fileids_T) :: FileID
  end type MLSFile_T

  ! The next datatype describes the information on the L1B data files in use

  type L1BInfo_T
    integer :: L1BOAId=0     ! The HDF ID (handle) for the L1BOA file
    ! Id(s) for the L1BRAD file(s)
    integer, dimension(:), pointer :: L1BRADIds=>NULL()
    character (LEN=FileNameLen) :: L1BOAFileName=""  ! L1BOA file name
    character (LEN=FileNameLen), dimension(:), pointer :: &
         & L1BRADFileNames=>NULL()
  end type L1BInfo_T

  ! --------------------------------------------------------------------------

  ! The TAI93 time range

  type TAI93_Range_T
    real(r8) :: startTime ! TAI93 format
    real(r8) :: endTime   ! TAI93 format
  end type TAI93_Range_T
  ! --------------------------------------------------------------------------

  interface FilterValues
    module procedure FilterValues_REAL, FilterValues_DOUBLE
    module procedure FilterValues_REAL_2d, FilterValues_DOUBLE_2d
    module procedure FilterValues_REAL_3d, FilterValues_DOUBLE_3d
  end interface
  
  interface IsFillValue
    module procedure IsFillValue_REAL, IsFillValue_DOUBLE
  end interface
  
  interface IsFinite
    module procedure IsFinite_REAL, IsFinite_DOUBLE, IsFinite_INTEGER
  end interface
  
  logical, parameter ::   DEEBUG = .false.

  contains

! ------------------------------------------------- FilterValues ---
  subroutine filterValues_REAL(a, ATAB, b, BTAB, warn, fillValue, precision)
      ! Return arrays filtered of any fillValues
      ! or where corresponding precision array < 0
      ! or whose values are not finite
      ! "Filter" means offending elements set to 0.
      ! Returned arrays are allocated and assigned values as appropriate
      ! Args
      real, dimension(:), intent(in)             :: a
      real, dimension(:), intent(in)             :: b
      real, dimension(:), intent(out)            :: atab
      real, dimension(:), intent(out)            :: btab
      logical, intent(out)                           :: warn
      real, optional, intent(in)                 :: fillValue
      real, dimension(:), optional, intent(in)   :: precision
      ! Internal variables
      integer                                        :: i
      real                                       :: myFillValue
      integer                                        :: n
      ! Executable
      myFillValue = 0.
      if ( present(FillValue) ) myFillValue = FillValue
      warn = .false.
      n=size(a)
      if ( n /= size(b) ) then
        call announce_error('a and b different sizes', n, size(b))
      elseif ( n /= size(atab) ) then
        call announce_error('a and atab different sizes', n, size(b))
      elseif ( n /= size(btab) ) then
        call announce_error('a and btab different sizes', n, size(b))
      elseif ( DEEBUG ) then
        call output('Filtering 1-d reals ', advance='yes')
        !call dump_name_v_pairs( (/n, size(b), size(atab), size(btab) /), &
        !  & 'size(a), size(b), size(atab), size(btab)', width=4)
      endif
      atab = a
      btab = b
      do i=1, N
        if ( .not. ieee_is_finite(a(i)) .or. .not. ieee_is_finite(b(i)) ) then
          atab(i) = myFillValue
          btab(i) = myFillValue
          warn = warn .or. ieee_is_finite(a(i)) .or. ieee_is_finite(b(i))
        endif
      enddo
      if ( present(fillValue) ) then
        do i=1, N
          if ( isFillValue(a(i), FillValue) .or. isFillValue(b(i), FillValue) ) then
            atab(i) = myFillValue
            btab(i) = myFillValue
          endif
        enddo
      endif
      if ( present(precision) ) then
        do i=1, N
          if ( (precision(i) < 0.) ) then
            atab(i) = myFillValue
            btab(i) = myFillValue
          endif
        enddo
      endif
  end subroutine filterValues_REAL

  subroutine filterValues_DOUBLE(a, ATAB, b, BTAB, warn, fillValue, precision)
      ! Return arrays filtered etc.
      ! Args
      double precision, dimension(:), intent(in)             :: a
      double precision, dimension(:), intent(in)             :: b
      double precision, dimension(:), intent(out)            :: atab
      double precision, dimension(:), intent(out)            :: btab
      logical, intent(out)                           :: warn
      double precision, optional, intent(in)                 :: fillValue
      double precision, dimension(:), optional, intent(in)   :: precision
      ! Internal variables
      integer                                        :: i
      double precision                                       :: myFillValue
      integer                                        :: n
      ! Executable
      myFillValue = 0.d0
      if ( present(FillValue) ) myFillValue = FillValue
      n=size(a)
      if ( n /= size(b) ) then
        call announce_error('a and b different sizes', n, size(b))
      elseif ( n /= size(atab) ) then
        call announce_error('a and atab different sizes', n, size(b))
      elseif ( n /= size(btab) ) then
        call announce_error('a and btab different sizes', n, size(b))
      elseif ( DEEBUG ) then
        call output('Filtering 1-d double precision ', advance='yes')
        !call dump_name_v_pairs( (/n, size(b), size(atab), size(btab)/) , &
        !  & 'size(a), size(b), size(atab), size(btab)', width=4)
      endif
      atab = a
      btab = b
      do i=1, N
        if ( .not. ieee_is_finite(a(i)) .or. .not. ieee_is_finite(b(i)) ) then
          atab(i) = myFillValue
          btab(i) = myFillValue
          warn = warn .or. ieee_is_finite(a(i)) .or. ieee_is_finite(b(i))
        endif
      enddo
      if ( present(fillValue) ) then
        do i=1, N
          if ( isFillValue(a(i), FillValue) .or. isFillValue(b(i), FillValue) ) then
            atab(i) = myFillValue
            btab(i) = myFillValue
          endif
        enddo
      endif
      if ( present(precision) ) then
        do i=1, N
          if ( (precision(i) < 0.d0) ) then
            atab(i) = myFillValue
            btab(i) = myFillValue
          endif
        enddo
      endif
  end subroutine filterValues_DOUBLE

  subroutine filterValues_REAL_2d(a, ATAB, b, BTAB, warn, fillValue, precision)
      ! Return arrays filtered etc.
      ! Args
      real, dimension(:,:), intent(in)             :: a
      real, dimension(:,:), intent(in)             :: b
      real, dimension(:,:), intent(out)            :: atab
      real, dimension(:,:), intent(out)            :: btab
      logical, intent(out)                           :: warn
      real, optional, intent(in)                 :: fillValue
      real, dimension(:,:), optional, intent(in)   :: precision
      ! Internal variables
      integer, dimension(2)                          :: shp
      real, dimension(size(a,1)*size(a,2))       :: a1
      real, dimension(size(b,1)*size(b,2))       :: b1
      ! Executable
      shp = shape(a)
      if ( DEEBUG ) then
        call dump_name_v_pairs(shp, width=4)
        call dump_name_v_pairs(shape(b), width=4)
      endif
      if ( present(precision) ) then
        call filterValues(reshape(a, (/shp(1)*shp(2)/)), &
        & a1, &
        & reshape(b, (/shp(1)*shp(2)/)), &
        & b1, &
        & warn, fillValue, reshape(precision, (/shp(1)*shp(2)/)) )
      else
        call filterValues(reshape(a, (/shp(1)*shp(2)/)), &
        & a1, &
        & reshape(b, (/shp(1)*shp(2)/)), &
        & b1, &
        & warn, fillValue)
      endif
      atab = reshape(a1, shp)
      btab = reshape(b1, shp)
  end subroutine filterValues_REAL_2d

  subroutine filterValues_DOUBLE_2d(a, ATAB, b, BTAB, warn, fillValue, precision)
      ! Return arrays filtered etc.
      ! Args
      double precision, dimension(:,:), intent(in)             :: a
      double precision, dimension(:,:), intent(in)             :: b
      double precision, dimension(:,:), intent(out)            :: atab
      double precision, dimension(:,:), intent(out)            :: btab
      logical, intent(out)                           :: warn
      double precision, optional, intent(in)                 :: fillValue
      double precision, dimension(:,:), optional, intent(in)   :: precision
      ! Internal variables
      integer, dimension(2)                          :: shp
      double precision, dimension(size(a,1)*size(a,2))       :: a1
      double precision, dimension(size(b,1)*size(b,2))       :: b1
      ! Executable
      shp = shape(a)
      if ( DEEBUG ) then
        call dump_name_v_pairs(shp, width=4)
        call dump_name_v_pairs(shape(b), width=4)
      endif
      if ( present(precision) ) then
        call filterValues(reshape(a, (/shp(1)*shp(2)/)), &
        & a1, &
        & reshape(b, (/shp(1)*shp(2)/)), &
        & b1, &
        & warn, fillValue, reshape(precision, (/shp(1)*shp(2)/)) )
      else
        call filterValues(reshape(a, (/shp(1)*shp(2)/)), &
        & a1, &
        & reshape(b, (/shp(1)*shp(2)/)), &
        & b1, &
        & warn, fillValue)
      endif
      atab = reshape(a1, shp)
      btab = reshape(b1, shp)
  end subroutine filterValues_DOUBLE_2d

  subroutine filterValues_REAL_3d(a, ATAB, b, BTAB, warn, fillValue, precision)
      ! Return arrays filtered etc.
      ! Args
      real, dimension(:,:,:), intent(in)             :: a
      real, dimension(:,:,:), intent(in)             :: b
      real, dimension(:,:,:), intent(out)            :: atab
      real, dimension(:,:,:), intent(out)            :: btab
      logical, intent(out)                           :: warn
      real, optional, intent(in)                 :: fillValue
      real, dimension(:,:,:), optional, intent(in)   :: precision
      ! Internal variables
      integer, dimension(3)                          :: shp
      real, dimension(size(a,1)*size(a,2)*size(a,3))       :: a1
      real, dimension(size(b,1)*size(b,2)*size(b,3))       :: b1
      ! Executable
      shp = shape(a)
      if ( DEEBUG ) then
        call dump_name_v_pairs(shp, width=4)
        call dump_name_v_pairs(shape(b), width=4)
      endif
      if ( present(precision) ) then
        call filterValues(reshape(a, (/shp(1)*shp(2)*shp(3)/)), &
        & a1, &
        & reshape(b, (/shp(1)*shp(2)*shp(3)/)), &
        & b1, &
        & warn, fillValue, reshape(precision, (/shp(1)*shp(2)*shp(3)/)) )
      else
        call filterValues(reshape(a, (/shp(1)*shp(2)*shp(3)/)), &
        & a1, &
        & reshape(b, (/shp(1)*shp(2)*shp(3)/)), &
        & b1, &
        & warn, fillValue)
      endif
      atab = reshape(a1, shp)
      btab = reshape(b1, shp)
  end subroutine filterValues_REAL_3d

  subroutine filterValues_DOUBLE_3d(a, ATAB, b, BTAB, warn, fillValue, precision)
      ! Return arrays filtered etc.
      ! Args
      double precision, dimension(:,:,:), intent(in)             :: a
      double precision, dimension(:,:,:), intent(in)             :: b
      double precision, dimension(:,:,:), intent(out)            :: atab
      double precision, dimension(:,:,:), intent(out)            :: btab
      logical, intent(out)                           :: warn
      double precision, optional, intent(in)                 :: fillValue
      double precision, dimension(:,:,:), optional, intent(in)   :: precision
      ! Internal variables
      integer, dimension(3)                          :: shp
      double precision, dimension(size(a,1)*size(a,2)*size(a,3))       :: a1
      double precision, dimension(size(b,1)*size(b,2)*size(b,3))       :: b1
      ! Executable
      shp = shape(a)
      if ( DEEBUG ) then
        call dump_name_v_pairs(shp, width=4)
        call dump_name_v_pairs(shape(b), width=4)
      endif
      if ( present(precision) ) then
        call filterValues(reshape(a, (/shp(1)*shp(2)*shp(3)/)), &
        & a1, &
        & reshape(b, (/shp(1)*shp(2)*shp(3)/)), &
        & b1, &
        & warn, fillValue, reshape(precision, (/shp(1)*shp(2)*shp(3)/)) )
      else
        call filterValues(reshape(a, (/shp(1)*shp(2)*shp(3)/)), &
        & a1, &
        & reshape(b, (/shp(1)*shp(2)*shp(3)/)), &
        & b1, &
        & warn, fillValue)
      endif
      atab = reshape(a1, shp)
      btab = reshape(b1, shp)
  end subroutine filterValues_DOUBLE_3d

  elemental function inRange(arg, range) result(relation)
    ! Is arg in range?
    integer, intent(in)       :: arg
    type(Range_T), intent(in) :: range
    logical                   :: relation
    relation = (arg < (range%top + 1)) .and. (arg > (range%bottom - 1))
  end function inRange
!=============================================================================
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
      call output('*** Error in MLSCommon module ***', advance='yes')
    endif
    call output(trim(message), advance='no')
    call blanks(3)
    if ( present(int1) ) write(*,'(i4)',advance='no') int1
    call blanks(3)
    if ( present(int2) ) write(*,'(i4)', advance='no') int2
    if ( .not. keepgoing ) stop
  end subroutine announce_error

  subroutine blanks(num, advance)
    integer, intent(in) :: num
    character(len=*), optional, intent(in) :: advance
    integer :: i
    character(len=3) :: myAdvance
    myAdvance = 'no'
    if ( present(advance) ) myAdvance = advance
    do i=1, num
      write(*, '(a1)', advance='no') ' '
    enddo
    if ( myAdvance == 'yes' ) write(*, '(a1)', advance='yes') ''
    
  end subroutine blanks

  subroutine output(str, advance)
    character(len=*), intent(in) :: str
    character(len=*), optional, intent(in) :: advance
    write(*, '(a1)', advance=advance) trim(str)
    
  end subroutine output

  ! ----------------------------------  DUMP_NAME_V_PAIRS  -----
  subroutine DUMP_NAME_V_PAIRS ( VALUES, WIDTH )
    integer, intent(in)                         :: values(:)
    integer, intent(in), optional :: WIDTH ! How many pairs per line (1)?
    
    integer :: J, K, L
    integer :: MyWidth
    character(len=24) :: myName
    MyWidth = 1
    if ( present(width) ) myWidth = max(width, 1)
    if ( size(values) < 1 ) return
    l = 0
    do j=1, size(values), MyWidth
      do k=1, MyWidth
        call blanks(3, advance='no')
        l = l + 1
        if ( l <= size(values) ) then
          write(myName, *) 'integer # ', l, ': '
          call output(myName,  advance='no')
          call blanks(3, advance='no')
          write(*,'(i4)', advance='no') values(l)
        end if
      end do
      call output(' ', advance='yes')
    end do

  end subroutine DUMP_NAME_V_PAIRS
! ------------------------------------------------- IsFillValue ---

  ! This family of routines checks to see if an arg is a fillValue
  elemental logical function IsFillValue_REAL ( A, FILLVALUE )
    real, intent(in) :: A
    real ,intent(in), optional :: FILLVALUE
    real  :: MYFILLVALUE
    myFillValue = DEFAULTUNDEFINEDVALUE
    if ( present(fillValue) ) myFillValue = fillValue
    IsFillValue_REAL = &
      & abs(a - myFillValue) < FILLVALUETOLERANCE
  end function IsFillValue_REAL

  elemental logical function IsFillValue_DOUBLE ( A, FILLVALUE )
    double precision, intent(in) :: A
    double precision ,intent(in), optional :: FILLVALUE
    double precision  :: MYFILLVALUE
    myFillValue = DEFAULTUNDEFINEDVALUE
    if ( present(fillValue) ) myFillValue = fillValue
    IsFillValue_DOUBLE = &
      & abs(a - myFillValue) < Real(FILLVALUETOLERANCE, kind(A))
  end function IsFillValue_DOUBLE

! ------------------------------------------------- IsFinite ---

  ! This family of routines checks to see if an arg is finite
  elemental logical function IsFinite_REAL ( A ) result( finite )
    real, intent(in) :: A
    finite = ieee_is_finite(a)
  end function isfinite_real
  elemental logical function IsFinite_DOUBLE ( A ) result( finite )
    double precision, intent(in) :: A
    finite = ieee_is_finite(a)
  end function isfinite_DOUBLE
  elemental logical function IsFinite_INTEGER ( A ) result( finite )
    integer, intent(in) :: A
    finite = .true.
  end function isfinite_INTEGER

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MLSCommon
!=============================================================================

!
! $Log$
! Revision 2.26  2005/10/19 22:53:01  vsnyder
! Move kinds to MLSKinds
!
! Revision 2.25  2005/06/22 17:25:49  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.24  2005/06/14 20:31:10  pwagner
! Added and changed some fields of MLSFile_T
!
! Revision 2.23  2005/05/31 17:49:15  pwagner
! Added new fields to MLSFile_T
!
! Revision 2.22  2005/05/12 20:46:46  pwagner
! Added filterValues and isFinite procedures (Should they be elsewhere?)
!
! Revision 2.21  2004/08/03 17:58:25  pwagner
! Now holds DEFAULTUNDEFINEDVALUE to be used elsewhere
!
! Revision 2.20  2004/06/10 01:00:50  vsnyder
! Move FindFirst, FindNext from MLSCommon to MLSSets
!
! Revision 2.19  2004/05/19 19:16:40  vsnyder
! Move MLSChunks_t to Chunks_m
!
! Revision 2.18  2004/01/09 00:38:04  pwagner
! Added FindNext function
!
! Revision 2.17  2003/06/20 19:31:39  pwagner
! Changes to allow direct writing of products
!
! Revision 2.16  2003/02/17 03:52:49  livesey
! Bit the bullet and changed rm to r4.
!
! Revision 2.15  2002/12/05 19:44:24  pwagner
! Moved MLSFile_T from MLSFiles to MLSCommon
!
! Revision 2.14  2002/11/06 00:16:48  pwagner
! Added toc/api blocks
!
! Revision 2.13  2002/10/08 00:09:11  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.12  2002/09/13 18:08:12  pwagner
! May change matrix precision rm from r8
!
! Revision 2.11  2002/08/28 22:16:18  pwagner
! Added rm, rv types
!
! Revision 2.10  2002/02/19 23:10:54  pwagner
! Added BareFNLen
!
! Revision 2.9  2002/01/09 23:51:27  pwagner
! Connected FileNameLen with PGSd_PC_FILE_PATH_MAX
!
! Revision 2.8  2001/11/14 18:03:32  livesey
! Changed FindFirst to return 0 not -1 if not found
!
! Revision 2.7  2001/09/09 02:47:58  livesey
! Moved FindFirst into MLSCommon
!
! Revision 2.6.2.3  2001/09/09 01:53:27  livesey
! Bug fix
!
! Revision 2.6.2.2  2001/09/09 01:35:46  livesey
! Moved FindFirst in from MLSL2Common
!
! Revision 2.6.2.1  2001/09/08 22:32:24  livesey
! Added RP and IP
!
! Revision 2.6  2001/04/20 23:10:53  livesey
! Initialised parameters in L1BINFO
!
! Revision 2.5  2001/03/10 18:48:17  livesey
! Really nullified the pointer!
!
! Revision 2.4  2001/03/10 07:06:46  livesey
! Nullified L1BRadfileNames in L1BInfo
!
! Revision 2.3  2001/02/09 00:38:55  livesey
! Various changes
!
! Revision 2.2  2001/01/26 23:46:35  pwagner
! Restored L1BInfo from l1/MLSL1Common back to lib/MLSCommon
!
! Revision 2.0  2000/09/05 17:41:06  dcuddy
! Change revision to 2.0
!
! Revision 1.14  2000/09/02 01:58:30  vsnyder
! Cosmetic changes
!
