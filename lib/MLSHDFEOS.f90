! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MLSHDFEOS

  ! This module contains MLS specific routines to do common HDFEOS
  ! tasks not specifically found in he5_swapi or he5_gdapi.

  use HDFEOS, only: swdefdfld, swdefgfld, swdefdim, swdiminfo
  use HDFEOS5, only: HE5_SWdefchunk, HE5_swdefdfld, HE5_SWdefgfld, &
    & HE5_swdefdim, HE5_swdiminfo
  use HE5_SWAPI, only: HE5_SWSETFILL, HE5_SWWRATTR, HE5_SWWRLATTR
  use HE5_SWAPI_CHARACTER_ARRAY, only: HE5_EHWRGLATT_CHARACTER_ARRAY
  use HE5_SWAPI_CHARACTER_SCALAR, only: HE5_EHWRGLATT_CHARACTER_SCALAR
  use HE5_SWAPI_DOUBLE, only: HE5_EHWRGLATT_DOUBLE, &
    & HE5_SWRDFLD_DOUBLE, HE5_SWRDFLD_DOUBLE_2D, HE5_SWRDFLD_DOUBLE_3D, &
    & HE5_SWWRFLD_DOUBLE, HE5_SWWRFLD_DOUBLE_2D, HE5_SWWRFLD_DOUBLE_3D
  use HE5_SWAPI_INTEGER, only: HE5_EHWRGLATT_INTEGER, &
    & HE5_SWRDFLD_INTEGER, HE5_SWWRFLD_INTEGER
  use HE5_SWAPI_REAL, only: HE5_EHWRGLATT_REAL, &
    & HE5_SWRDFLD_REAL, HE5_SWRDFLD_REAL_2D, HE5_SWRDFLD_REAL_3D, &
    & HE5_SWWRFLD_REAL, HE5_SWWRFLD_REAL_2D, HE5_SWWRFLD_REAL_3D
  use MLSFiles, only: HDFVERSION_4, HDFVERSION_5, WILDCARDHDFVERSION, &
    & mls_hdf_version
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
  use MLSStrings, only: GetStringElement, NumStringElements
  use SWAPI_DOUBLE, only: SWRDFLD_DOUBLE, SWRDFLD_DOUBLE_2D, SWRDFLD_DOUBLE_3D, &
    &                     SWWRFLD_DOUBLE, SWWRFLD_DOUBLE_2D, SWWRFLD_DOUBLE_3D
  use SWAPI_INTEGER, only: SWRDFLD_INTEGER, SWRDFLD_INTEGER_2D, SWRDFLD_INTEGER_3D, &
    &                      SWWRFLD_INTEGER, SWWRFLD_INTEGER_2D, SWWRFLD_INTEGER_3D
  use SWAPI_REAL, only: SWRDFLD_REAL, SWRDFLD_REAL_2D, SWRDFLD_REAL_3D, &
    &                   SWWRFLD_REAL, SWWRFLD_REAL_2D, SWWRFLD_REAL_3D

  implicit NONE
  private

  public :: HE5_EHWRGLATT, MLS_DFLDSETUP, MLS_GFLDSETUP, &
    & MLS_SWDEFDIM, MLS_SWDIMINFO, MLS_SWRDFLD, MLS_SWSETFILL, MLS_SWWRFLD
  logical, parameter :: HE5_SWSETFILL_BROKEN = .true.
  character(len=*), parameter :: SETFILLTITLE = '_FillValue'

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

! HE5_EHWRGLATT     Sets the global attributes at file level
! MLS_GFLDSETUP     Sets up a geolocation field in a swath
! MLS_SWRDFLD       Reads a field from a swath, data or geolocation
! MLS_SWSETFILL     Sets the fill value for one or more fields
! MLS_SWWRFLD       Writes a field to a swath, data or geolocation
! === (end of toc) ===

! === (start of api) ===
! int HE5_EHWRGLATT (int fileID, char* attrName, int datatype, int count, value)
!     value can be one of:
!    {char* value, char* value(:), int value(:), r4 value(:), r8 value(:)}
! int MLS_SWSETFILL (int swathID, char* names, int datatype, value) 
!     value can be one of:
!    {char* value, int value, r4 value, r8 value}
! === (end of api) ===
  interface HE5_EHWRGLATT   ! From its name, might better be in he5_ehapi.f90
    module procedure HE5_EHWRGLATT_CHARACTER_SCALAR, HE5_EHWRGLATT_DOUBLE, &
    HE5_EHWRGLATT_INTEGER, HE5_EHWRGLATT_REAL, HE5_EHWRGLATT_CHARACTER_ARRAY
  end interface

  interface MLS_SWSETFILL
    module procedure MLS_SWSETFILL_DOUBLE, &
      & MLS_SWSETFILL_INTEGER, MLS_SWSETFILL_REAL
  end interface

  interface MLS_SWRDFLD
    module procedure &
      & MLS_SWRDFLD_DOUBLE_1D, MLS_SWRDFLD_DOUBLE_2D, MLS_SWRDFLD_DOUBLE_3D, &
      & MLS_SWRDFLD_INTEGER, &
      & MLS_SWRDFLD_REAL_1D, MLS_SWRDFLD_REAL_2D, MLS_SWRDFLD_REAL_3D
  end interface

  interface MLS_SWWRFLD
    module procedure &
      & MLS_SWWRFLD_DOUBLE_1D, MLS_SWWRFLD_DOUBLE_2D, MLS_SWWRFLD_DOUBLE_3D, &
      & MLS_SWWRFLD_INTEGER, &
      & MLS_SWWRFLD_REAL_1D, MLS_SWWRFLD_REAL_2D, MLS_SWWRFLD_REAL_3D
  end interface

contains ! ======================= Public Procedures =========================

  integer function MLS_SWdefdim ( SWATHID, DIMNAME, DIMSIZE, FILENAME, &
    & hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 7
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: DIMNAME       ! Dimension name
    integer, intent(in)          :: DIMSIZE
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    ! Internal variables
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    MLS_SWdefdim = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_SWdefdim = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWdefdim' )
      endif
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = mls_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_SWdefdim = swdefdim(swathid, DIMNAME, DIMSIZE)
    case (HDFVERSION_5)
        MLS_SWdefdim = HE5_SWdefdim(swathid, DIMNAME, DIMSIZE)
    case default
      MLS_SWdefdim = -1
    end select
    if ( myDontFail .or. MLS_SWdefdim /= -1 ) return
    CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to define dimension ' // trim(dimname) )

  end function MLS_SWdefdim

  integer function MLS_SWdiminfo ( SWATHID, DIMNAME, FILENAME, &
    & hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 7
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: DIMNAME       ! Dimension name
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    ! Internal variables
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    MLS_SWdiminfo = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_SWdiminfo = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWdiminfo' )
      endif
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = mls_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_SWdiminfo = SWdiminfo(swathid, DIMNAME)
    case (HDFVERSION_5)
        MLS_SWdiminfo = HE5_SWdiminfo(swathid, DIMNAME)
    case default
      MLS_SWdiminfo = -1
    end select
    if ( myDontFail .or. MLS_SWdiminfo /= -1 ) return
    CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to get info on dimension ' // trim(dimname) )

  end function MLS_SWdiminfo

  integer function MLS_dfldsetUP ( SWATHID, FIELDNAME, DIMNAME, MAXDIMList, &
    & DATATYPE, MERGE, CHUNK_RANK, CHUNK_DIMS, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 7
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    character(len=*), intent(in) :: DIMNAME       ! Dimension name
    character(len=*), intent(in) :: MAXDIMLIST
    integer, intent(in)          :: DATATYPE
    integer, intent(in)          :: MERGE
    integer, intent(in)          :: CHUNK_RANK
    integer, dimension(RANK), intent(in) :: CHUNK_DIMS
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    ! Internal variables
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    mls_dfldsetup = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        mls_dfldsetup = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to mls_dfldsetup' )
      endif
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = mls_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      mls_dfldsetup = swdefdfld(swathid, FIELDName, DIMNAME, Datatype, &
         & MERGE)
    case (HDFVERSION_5)
      if ( chunk_rank /= 0 ) &
        & mls_dfldsetup = HE5_SWdefchunk(swathid, chunk_rank, chunk_dims)
      if ( mls_dfldsetup == 0 ) &
        & mls_dfldsetup = HE5_SWdefdfld(swathid, FIELDName, DIMNAME, MAXDIMLIST, &
        & Datatype, MERGE)
    case default
      mls_dfldsetup = -1
    end select
    if ( myDontFail .or. mls_dfldsetup /= -1 ) return
    CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to set up data field ' // trim(fieldname) )

  end function MLS_dfldsetUP

  integer function MLS_GFLDSETUP ( SWATHID, FIELDNAME, DIMNAME, MAXDIMList, &
    & DATATYPE, MERGE, CHUNK_RANK, CHUNK_DIMS, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 7
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    character(len=*), intent(in) :: DIMNAME       ! Dimension name
    character(len=*), intent(in) :: MAXDIMLIST
    integer, intent(in)          :: DATATYPE
    integer, intent(in)          :: MERGE
    integer, intent(in)          :: CHUNK_RANK
    integer, dimension(RANK), intent(in) :: CHUNK_DIMS
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    ! Internal variables
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    mls_gfldsetup = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        mls_gfldsetup = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to mls_gfldsetup' )
      endif
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = mls_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      mls_gfldsetup = swdefgfld(swathid, FIELDName, DIMNAME, Datatype, &
         & MERGE)
    case (HDFVERSION_5)
      if ( chunk_rank /= 0 ) &
        & mls_gfldsetup = HE5_SWdefchunk(swathid, chunk_rank, chunk_dims)
      if ( mls_gfldsetup == 0 ) &
        & mls_gfldsetup = HE5_SWdefgfld(swathid, FIELDName, DIMNAME, MAXDIMLIST, &
        & Datatype, MERGE)
    case default
      mls_gfldsetup = -1
    end select
    if ( myDontFail .or. mls_gfldsetup /= -1 ) return
    CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to set up geoloc field ' // trim(fieldname) )

  end function MLS_GFLDSETUP

  integer function MLS_SWSETFILL_DOUBLE ( SWATHID, FIELDNAMES, DATATYPE, &
    & FILLVALUE )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAMES     ! Field names
    integer, intent(in) :: DATATYPE
    double precision, intent(in) :: FILLVALUE

    integer :: Field
    integer :: numFields
    character(len=len(FIELDNAMES)) :: FIELDNAME
    mls_swsetfill_double = 0
    numFields = NumStringElements(fieldnames, .false.)
    if ( numFields < 1 ) return
    do Field=1, numFields
      call GetStringElement(fieldnames, fieldname, Field, .true.)
      if ( HE5_SWSETFILL_BROKEN ) then
        mls_swsetfill_double = he5_swwrlattr(swathid, trim(fieldname), &
          & trim(SETFILLTITLE), &
          & DATATYPE, 1, (/ FILLVALUE /) )
      else
        mls_swsetfill_double = HE5_SWsetfill(swathid, trim(fieldname), &
        & datatype, fillvalue)
      endif
      if ( mls_swsetfill_double == -1 ) return
    enddo

  end function MLS_SWSETFILL_DOUBLE

  integer function MLS_SWSETFILL_INTEGER ( SWATHID, FIELDNAMES, DATATYPE, &
    & FILLVALUE )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAMES     ! Field names
    integer, intent(in) :: DATATYPE
    integer, intent(in) :: FILLVALUE

    integer :: Field
    integer :: numFields
    character(len=len(FIELDNAMES)) :: FIELDNAME
    mls_swsetfill_integer = 0
    numFields = NumStringElements(fieldnames, .false.)
    if ( numFields < 1 ) return
    do Field=1, numFields
      call GetStringElement(fieldnames, fieldname, Field, .true.)
      if ( HE5_SWSETFILL_BROKEN ) then
        mls_swsetfill_integer = he5_swwrlattr(swathid, trim(fieldname), &
          & trim(SETFILLTITLE), &
          & DATATYPE, 1, (/ FILLVALUE /) )
      else
        mls_swsetfill_integer = HE5_SWsetfill(swathid, trim(fieldname), &
        & datatype, fillvalue)
      endif
      if ( mls_swsetfill_integer == -1 ) return
    enddo

  end function MLS_SWSETFILL_INTEGER

  integer function MLS_SWSETFILL_REAL ( SWATHID, FIELDNAMES, DATATYPE, &
    & FILLVALUE )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAMES     ! Field names
    integer, intent(in) :: DATATYPE
    real, intent(in) :: FILLVALUE
    integer, external :: HE5_SWsetfill
    integer :: Field
    integer :: numFields
    character(len=len(FIELDNAMES)) :: FIELDNAME
    mls_swsetfill_real = 0
    numFields = NumStringElements(fieldnames, .false.)
    if ( numFields == -1 ) return
    do Field=1, numFields
      call GetStringElement(fieldnames, fieldname, Field, .true.)
      if ( HE5_SWSETFILL_BROKEN ) then
        mls_swsetfill_real = he5_swwrlattr(swathid, trim(fieldname), &
          & trim(SETFILLTITLE), &
          & DATATYPE, 1, (/ FILLVALUE /) )
      else
        mls_swsetfill_real = HE5_SWsetfill(swathid, trim(fieldname), &
        & datatype, fillvalue)
      endif
      if ( mls_swsetfill_real == -1 ) return
    enddo

  end function MLS_SWSETFILL_REAL

  integer function MLS_SWRDFLD_DOUBLE_1D ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 1
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    double precision, dimension(:), intent(out) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    ! Internal variables
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    mls_swrdfld_double_1d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        mls_swrdfld_double_1d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWRDFLD_DOUBLE_1D' )
      endif
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = mls_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      mls_swrdfld_double_1d = SWRDFLD_DOUBLE(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      mls_swrdfld_double_1d = HE5_SWRDFLD_DOUBLE(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case default
      mls_swrdfld_double_1d = -1
    end select
    if ( myDontFail .or. mls_swrdfld_double_1d /= -1 ) return
    CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to read 1-d double field ' // trim(fieldname) )

  end function MLS_SWRDFLD_DOUBLE_1D

  integer function MLS_SWRDFLD_DOUBLE_2d ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 2
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    double precision, dimension(:,:), intent(out) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    ! Internal variables
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    mls_swrdfld_double_2d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        mls_swrdfld_double_2d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWRDFLD_DOUBLE_2d' )
      endif
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = mls_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      mls_swrdfld_double_2d = SWRDFLD_DOUBLE_2D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      mls_swrdfld_double_2d = HE5_SWRDFLD_DOUBLE_2D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case default
      mls_swrdfld_double_2d = -1
    end select
    if ( myDontFail .or. mls_swrdfld_double_2d /= -1 ) return
    CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to read 2-d double field ' // trim(fieldname) )

  end function MLS_SWRDFLD_DOUBLE_2d

  integer function MLS_SWRDFLD_DOUBLE_3d ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 3
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    double precision, dimension(:,:,:), intent(out) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    ! Internal variables
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    mls_swrdfld_double_3d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        mls_swrdfld_double_3d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWRDFLD_DOUBLE_3d' )
      endif
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = mls_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      mls_swrdfld_double_3d = SWRDFLD_DOUBLE_3D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      mls_swrdfld_double_3d = HE5_SWRDFLD_DOUBLE_3D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case default
      mls_swrdfld_double_3d = -1
    end select
    if ( myDontFail .or. mls_swrdfld_double_3d /= -1 ) return
    CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to read 1-d double field ' // trim(fieldname) )

  end function MLS_SWRDFLD_DOUBLE_3d

  integer function MLS_SWRDFLD_INTEGER ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 1
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    integer, dimension(:), intent(out) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    ! Internal variables
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    mls_swrdfld_integer = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        mls_swrdfld_integer = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWRDFLD_integer' )
      endif
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = mls_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      mls_swrdfld_integer = SWRDFLD_INTEGER(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      mls_swrdfld_integer = HE5_SWRDFLD_INTEGER(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case default
      mls_swrdfld_integer = -1
    end select
    if ( myDontFail .or. mls_swrdfld_integer /= -1 ) return
    CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to read 1-d double field ' // trim(fieldname) )

  end function MLS_SWRDFLD_INTEGER

  integer function MLS_SWRDFLD_REAL_1D ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 1
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    real, dimension(:), intent(out) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    ! Internal variables
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    mls_swrdfld_REAL_1d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        mls_swrdfld_REAL_1d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWRDFLD_REAL_1D' )
      endif
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = mls_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      mls_swrdfld_REAL_1d = SWRDFLD_REAL(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      mls_swrdfld_REAL_1d = HE5_SWRDFLD_REAL(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case default
      mls_swrdfld_REAL_1d = -1
    end select
    if ( myDontFail .or. mls_swrdfld_REAL_1d /= -1 ) return
    CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to read 1-d REAL field ' // trim(fieldname) )

  end function MLS_SWRDFLD_REAL_1D

  integer function MLS_SWRDFLD_REAL_2d ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 2
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    real, dimension(:,:), intent(out) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    ! Internal variables
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    mls_swrdfld_REAL_2d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        mls_swrdfld_REAL_2d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWRDFLD_REAL_2d' )
      endif
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = mls_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      mls_swrdfld_REAL_2d = SWRDFLD_REAL_2D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      mls_swrdfld_REAL_2d = HE5_SWRDFLD_REAL_2D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case default
      mls_swrdfld_REAL_2d = -1
    end select
    if ( myDontFail .or. mls_swrdfld_REAL_2d /= -1 ) return
    CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to read 2-d REAL field ' // trim(fieldname) )

  end function MLS_SWRDFLD_REAL_2d

  integer function MLS_SWRDFLD_REAL_3d ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 3
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    real, dimension(:,:,:), intent(out) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    ! Internal variables
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    mls_swrdfld_REAL_3d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        mls_swrdfld_REAL_3d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWRDFLD_REAL_3d' )
      endif
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = mls_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      mls_swrdfld_REAL_3d = SWRDFLD_REAL_3D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      mls_swrdfld_REAL_3d = HE5_SWRDFLD_REAL_3D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case default
      mls_swrdfld_REAL_3d = -1
    end select
    if ( myDontFail .or. mls_swrdfld_REAL_3d /= -1 ) return
    CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to read 1-d REAL field ' // trim(fieldname) )

  end function MLS_SWRDFLD_REAL_3d

  integer function MLS_SWWRFLD_DOUBLE_1D ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 1
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    double precision, dimension(:), intent(in) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    ! Internal variables
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    mls_SWWRFLD_double_1d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        mls_SWWRFLD_double_1d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWWRFLD_DOUBLE_1D' )
      endif
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = mls_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      mls_SWWRFLD_double_1d = SWWRFLD_DOUBLE(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      mls_SWWRFLD_double_1d = HE5_SWWRFLD_DOUBLE(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case default
      mls_SWWRFLD_double_1d = -1
    end select
    if ( myDontFail .or. mls_SWWRFLD_double_1d /= -1 ) return
    CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to read 1-d double field ' // trim(fieldname) )

  end function MLS_SWWRFLD_DOUBLE_1D

  integer function MLS_SWWRFLD_DOUBLE_2d ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 2
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    double precision, dimension(:,:), intent(in) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    ! Internal variables
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    mls_SWWRFLD_double_2d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        mls_SWWRFLD_double_2d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWWRFLD_DOUBLE_2d' )
      endif
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = mls_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      mls_SWWRFLD_double_2d = SWWRFLD_DOUBLE_2D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      mls_SWWRFLD_double_2d = HE5_SWWRFLD_DOUBLE_2D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case default
      mls_SWWRFLD_double_2d = -1
    end select
    if ( myDontFail .or. mls_SWWRFLD_double_2d /= -1 ) return
    CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to read 2-d double field ' // trim(fieldname) )

  end function MLS_SWWRFLD_DOUBLE_2d

  integer function MLS_SWWRFLD_DOUBLE_3d ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 3
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    double precision, dimension(:,:,:), intent(in) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    ! Internal variables
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    mls_SWWRFLD_double_3d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        mls_SWWRFLD_double_3d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWWRFLD_DOUBLE_3d' )
      endif
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = mls_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      mls_SWWRFLD_double_3d = SWWRFLD_DOUBLE_3D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      mls_SWWRFLD_double_3d = HE5_SWWRFLD_DOUBLE_3D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case default
      mls_SWWRFLD_double_3d = -1
    end select
    if ( myDontFail .or. mls_SWWRFLD_double_3d /= -1 ) return
    CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to read 1-d double field ' // trim(fieldname) )

  end function MLS_SWWRFLD_DOUBLE_3d

  integer function MLS_SWWRFLD_INTEGER ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 1
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    integer, dimension(:), intent(in) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    ! Internal variables
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    mls_SWWRFLD_integer = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        mls_SWWRFLD_integer = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWWRFLD_integer' )
      endif
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = mls_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      mls_SWWRFLD_integer = SWWRFLD_INTEGER(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      mls_SWWRFLD_integer = HE5_SWWRFLD_INTEGER(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case default
      mls_SWWRFLD_integer = -1
    end select
    if ( myDontFail .or. mls_SWWRFLD_integer /= -1 ) return
    CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to read 1-d double field ' // trim(fieldname) )

  end function MLS_SWWRFLD_INTEGER

  integer function MLS_SWWRFLD_REAL_1D ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 1
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    real, dimension(:), intent(in) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    ! Internal variables
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    mls_SWWRFLD_REAL_1d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        mls_SWWRFLD_REAL_1d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWWRFLD_REAL_1D' )
      endif
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = mls_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      mls_SWWRFLD_REAL_1d = SWWRFLD_REAL(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      mls_SWWRFLD_REAL_1d = HE5_SWWRFLD_REAL(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case default
      mls_SWWRFLD_REAL_1d = -1
    end select
    if ( myDontFail .or. mls_SWWRFLD_REAL_1d /= -1 ) return
    CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to read 1-d REAL field ' // trim(fieldname) )

  end function MLS_SWWRFLD_REAL_1D

  integer function MLS_SWWRFLD_REAL_2d ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 2
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    real, dimension(:,:), intent(in) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    ! Internal variables
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    mls_SWWRFLD_REAL_2d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        mls_SWWRFLD_REAL_2d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWWRFLD_REAL_2d' )
      endif
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = mls_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      mls_SWWRFLD_REAL_2d = SWWRFLD_REAL_2D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      mls_SWWRFLD_REAL_2d = HE5_SWWRFLD_REAL_2D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case default
      mls_SWWRFLD_REAL_2d = -1
    end select
    if ( myDontFail .or. mls_SWWRFLD_REAL_2d /= -1 ) return
    CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to read 2-d REAL field ' // trim(fieldname) )

  end function MLS_SWWRFLD_REAL_2d

  integer function MLS_SWWRFLD_REAL_3d ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 3
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    real, dimension(:,:,:), intent(in) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    ! Internal variables
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    mls_SWWRFLD_REAL_3d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        mls_SWWRFLD_REAL_3d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWWRFLD_REAL_3d' )
      endif
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = mls_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      mls_SWWRFLD_REAL_3d = SWWRFLD_REAL_3D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      mls_SWWRFLD_REAL_3d = HE5_SWWRFLD_REAL_3D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case default
      mls_SWWRFLD_REAL_3d = -1
    end select
    if ( myDontFail .or. mls_SWWRFLD_REAL_3d /= -1 ) return
    CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to read 1-d REAL field ' // trim(fieldname) )

  end function MLS_SWWRFLD_REAL_3d

! ======================= Private Procedures =========================  

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MLSHDFEOS

! $Log$
! Revision 2.4  2003/04/17 23:04:51  pwagner
! Now sits between L2GPData and HDFEOS(5)
!
! Revision 2.3  2003/04/15 23:16:46  pwagner
! Removed prints; begun MLS_SWRDFLD stuff
!
! Revision 2.2  2003/04/15 21:58:54  pwagner
! Now sets _FillValue attribute because swsetfill seems broken
!
! Revision 2.1  2003/04/11 23:28:44  pwagner
! FIrst commit
!
