! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MLSHDFEOS

  ! This module contains MLS specific routines to do common HDFEOS
  ! tasks not specifically found in he5_swapi or he5_gdapi.

  use Hdf, only: DFNT_CHAR8, DFNT_FLOAT32, DFNT_FLOAT64, &
    & DFNT_INT8, DFNT_INT16, DFNT_INT32, DFNT_INT64
  use HDFEOS, only: gdattach, gdcreate, &
    & swattach, swcreate, swdefdfld, swdefgfld, swdefdim, swdetach, &
    & swdiminfo, swinqdflds, swinqswath
  use HDFEOS5, only: HE5T_NATIVE_FLOAT, HE5T_NATIVE_DOUBLE, HE5T_NATIVE_SCHAR, &
    & HE5T_NATIVE_INT, HE5T_NATIVE_INT8, HE5T_NATIVE_INT16, HE5T_NATIVE_INT64
  use HDFEOS5, only: HE5_GDattach, HE5_GDcreate, &
    & HE5_SWattach, HE5_SWcreate, HE5_SWdefchunk, HE5_SWdetach, HE5_swdefdfld, &
    & HE5_SWdefgfld, HE5_swdefdim, HE5_swdiminfo, HE5_swinqdflds, &
    & HE5_swinqswath
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
  use MLSStrings, only: GetStringElement, NumStringElements, StringElementNum
  use SWAPI_DOUBLE, only: SWRDFLD_DOUBLE, SWRDFLD_DOUBLE_2D, SWRDFLD_DOUBLE_3D, &
    &                     SWWRFLD_DOUBLE, SWWRFLD_DOUBLE_2D, SWWRFLD_DOUBLE_3D
  use SWAPI_INTEGER, only: SWRDFLD_INTEGER, SWRDFLD_INTEGER_2D, SWRDFLD_INTEGER_3D, &
    &                      SWWRFLD_INTEGER, SWWRFLD_INTEGER_2D, SWWRFLD_INTEGER_3D
  use SWAPI_REAL, only: SWRDFLD_REAL, SWRDFLD_REAL_2D, SWRDFLD_REAL_3D, &
    &                   SWWRFLD_REAL, SWWRFLD_REAL_2D, SWWRFLD_REAL_3D

  implicit NONE
  private

  public :: HE5_EHWRGLATT, MLS_DFLDSETUP, MLS_GFLDSETUP, &
    & MLS_SWDEFDIM, MLS_SWDIMINFO, MLS_SWRDFLD, MLS_SWSETFILL, MLS_SWWRFLD, &
    & MLS_SWATTACH, MLS_SWCREATE, MLS_SWDETACH, MLS_GDCREATE, &
    & mls_swath_in_file
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

  interface MLS_SWATH_IN_FILE
    module procedure MLS_SWATH_IN_FILE_ARR, &
      & MLS_SWATH_IN_FILE_SCA
  end interface

  interface MLS_SWSETFILL
    module procedure MLS_SWSETFILL_DOUBLE, &
      & MLS_SWSETFILL_INTEGER, MLS_SWSETFILL_REAL
  end interface

  interface MLS_SWRDFLD
    module procedure &
      & MLS_SWRDFLD_CHAR_1D, &
      & MLS_SWRDFLD_DOUBLE_1D, MLS_SWRDFLD_DOUBLE_2D, MLS_SWRDFLD_DOUBLE_3D, &
      & MLS_SWRDFLD_INTEGER, &
      & MLS_SWRDFLD_REAL_1D, MLS_SWRDFLD_REAL_2D, MLS_SWRDFLD_REAL_3D
  end interface

  interface MLS_SWWRFLD
    module procedure &
      & MLS_SWWRFLD_CHAR_1D, &
      & MLS_SWWRFLD_DOUBLE_1D, MLS_SWWRFLD_DOUBLE_2D, MLS_SWWRFLD_DOUBLE_3D, &
      & MLS_SWWRFLD_INTEGER, &
      & MLS_SWWRFLD_REAL_1D, MLS_SWWRFLD_REAL_2D, MLS_SWWRFLD_REAL_3D
  end interface

  integer, public, parameter :: MAXNODFIELDS = 1000
  integer, public, parameter :: MAXDLISTLENGTH = 4096

  ! Print debugging stuff?
  logical, parameter :: DEEBUG = .true.  

contains ! ======================= Public Procedures =========================

  integer function MLS_GDCREATE ( FILEID, GRIDNAME, &
   &  xdimsize, ydimsize, upleft, lowright, FileName, hdfVersion )
    integer, intent(in) :: FILEID      ! ID returned by mls_swopen
    character(len=*), intent(in) :: GRIDNAME       ! Swath name
    integer, intent(in) :: xdimsize
    integer, intent(in) :: ydimsize
    double precision, dimension(2), intent(in) :: upleft
    double precision, dimension(2), intent(in) :: lowright
    character(len=*), optional, intent(in) :: FILENAME  ! File name
    integer, optional, intent(in) :: hdfVersion
    ! Internal variables
    logical :: alreadyThere
    integer :: myHdfVersion
    logical :: needsFileName
    logical, parameter :: MUSTCREATE = .true.
    MLS_GDCREATE = 0
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_GDCREATE' )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = mls_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      if ( MUSTCREATE ) then
        alreadyThere = .false.
      else
        alreadyThere = (gdattach(FileID, trim(GRIDNAME)) >= 0)
      endif
      if ( alreadyThere ) then
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'GRIDNAME call to MLS_GDCREATE already there: ' // trim(GRIDNAME))
      endif
      MLS_GDCREATE = gdcreate(Fileid, trim(GRIDNAME), &
        & xdimsize, ydimsize, upleft, lowright)
    case (HDFVERSION_5)
      if ( MUSTCREATE ) then
        alreadyThere = .false.
      else
        alreadyThere = (he5_gdattach(FileID, trim(GRIDNAME)) >= 0)
      endif
      if ( alreadyThere ) then
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'GRIDNAME call to MLS_GDCREATE already there: ' // trim(GRIDNAME))
      endif
      MLS_GDCREATE = he5_gdcreate(Fileid, trim(GRIDNAME), &
        & xdimsize, ydimsize, upleft, lowright)
    case default
      MLS_GDCREATE = -1
    end select
    if ( MLS_GDCREATE /= -1 ) return
    CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to create grid name ' // trim(GRIDNAME) )

  end function MLS_GDCREATE

  integer function MLS_SWATTACH ( FILEID, SWATHNAME, FileName, &
    &  hdfVersion, DONTFAIL )
    integer, intent(in) :: FILEID      ! ID returned by mls_swopen
    character(len=*), intent(in) :: SWATHNAME       ! Swath name
    character(len=*), optional, intent(in) :: FILENAME  ! File name
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    ! Internal variables
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    if (Deebug) print *, 'swattaching ', trim(SWATHNAME)
    MLS_SWATTACH = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_SWATTACH = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWATTACH' )
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
      MLS_SWATTACH = swattach(Fileid, trim(swathName))
    case (HDFVERSION_5)
      MLS_SWATTACH = he5_swattach(Fileid, trim(swathName))
    case default
      MLS_SWATTACH = -1
    end select
    if (Deebug) print *, ' (swath id is ', MLS_SWATTACH, ')'
    if ( MLS_SWATTACH /= -1 ) then
      return
    elseif ( myDontFail ) then
      CALL MLSMessage ( MLSMSG_Warning, moduleName,  &
          & 'Failed to attach swath name ' // trim(swathname) )
    else
      CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to attach swath name ' // trim(swathname) )
    endif
  end function MLS_SWATTACH

  integer function MLS_SWCREATE ( FILEID, SWATHNAME, FileName, hdfVersion )
    use hdf5, only: h5eSet_auto_f
    integer, intent(in) :: FILEID      ! ID returned by mls_swopen
    character(len=*), intent(in) :: SWATHNAME       ! Swath name
    character(len=*), optional, intent(in) :: FILENAME  ! File name
    integer, optional, intent(in) :: hdfVersion
    ! Internal variables
    logical :: alreadyThere
    integer :: myHdfVersion
    logical :: needsFileName
    integer :: status
    logical, parameter :: ALWAYSTRYSWATTACH = .true.
    logical, parameter :: MUSTCREATE = .true.
    MLS_SWCREATE = 0
    ! All necessary input supplied?
    print *, 'Now in mls_swcreate'
    print *, 'SWATHNAME: ', trim(SWATHNAME)
    if ( present(filename) ) print *, 'filename: ', trim(filename)
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FileName)) then
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWCREATE' )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = mls_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      if ( MUSTCREATE ) then
        alreadyThere = .false.
      elseif ( .not. present(filename) .or. ALWAYSTRYSWATTACH ) then
        if(DEEBUG) print *, 'About to call swattach with FileID: ', Fileid, ' swathname ', &
          & trim(swathName)
        alreadyThere = (swattach(FileID, trim(swathName)) >= 0)
      else
        if(DEEBUG) print *, 'About to call mls_swath_in_file'
        alreadyThere = mls_swath_in_file(FileName, swathName, HDFVERSION_4)
      endif
      if ( alreadyThere ) then
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Swathname call to MLS_SWCREATE already there: ' // trim(swathName))
      endif
      if(DEEBUG) print *, 'About to call swcreate with FileID: ', Fileid, ' swathname ', &
        & trim(swathName)
      MLS_SWCREATE = swcreate(Fileid, trim(swathName))
    case (HDFVERSION_5)
      call h5eSet_auto_f ( 0, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to turn error messages off in MLS_SWCREATE')
      if ( MUSTCREATE ) then
        alreadyThere = .false.
      elseif ( .not. present(filename) .or. ALWAYSTRYSWATTACH  ) then
        if(DEEBUG) print *, 'About to call he5_swattach with FileID: ', Fileid, ' swathname ', &
          & trim(swathName)
        alreadyThere = (he5_swattach(FileID, trim(swathName)) >= 0)
      else
        if(DEEBUG) print *, 'About to call mls_swath_in_file'
        alreadyThere = mls_swath_in_file(FileName, swathName, HDFVERSION_5)
      endif
      call h5eSet_auto_f ( 1, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to turn error messages back on in MLS_SWCREATE')
      if ( alreadyThere ) then
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Swathname call to MLS_SWCREATE already there: ' // trim(swathName))
      endif
      MLS_SWCREATE = he5_swcreate(Fileid, trim(swathName))
    case default
      MLS_SWCREATE = -1
    end select
    if ( MLS_SWCREATE /= -1 ) return
    CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to create swath name ' // trim(swathname) )

  end function MLS_SWCREATE

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

  integer function MLS_SWDETACH ( SWATHID, FileName, hdfVersion )
    integer, intent(in) :: SWATHID      ! ID returned by mls_swattach
    character(len=*), optional, intent(in) :: FILENAME  ! File name
    integer, optional, intent(in) :: hdfVersion
    ! Internal variables
    integer :: myHdfVersion
    logical :: needsFileName
    if (Deebug) print *, 'swattaching ', swathid
    MLS_SWDETACH = 0
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWDETACH' )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = mls_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_SWDETACH = swdetach(swathid)
    case (HDFVERSION_5)
      MLS_SWDETACH = he5_swdetach(swathid)
    case default
      MLS_SWDETACH = -1
    end select
    if ( MLS_SWDETACH /= -1 ) return
    CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to detach swath id ' )

  end function MLS_SWDETACH

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
        & he2he5_DataType(Datatype), MERGE)
!>       if ( Datatype == DFNT_CHAR8 ) then
!>         if(DEEBUG) print *,' Hey, we just tried to set up a char-valued hdfeos5 field'
!>         if(DEEBUG) print *,' Data type: ', Datatype
!>         if(DEEBUG) print *,' he2he5_DataType(Datatype): ', he2he5_DataType(Datatype)
!>         if(DEEBUG) print *,' HE5T_NATIVE_SCHAR: ', HE5T_NATIVE_SCHAR
!>       endif
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
    print *, 'mls_gfldsetup'
    print *, 'FIELDName: ', trim(FIELDName)
    print *, 'myHdfVersion: ', myHdfVersion
    print *, 'swathid: ', swathid
    print *, 'Datatype: ', Datatype
    print *, 'he5_Datatype: ', he2he5_DataType(Datatype)
    print *, 'chunk_rank: ', chunk_rank
    print *, 'chunk_dims: ', chunk_dims
    print *, 'DIMNAME: ', trim(DIMNAME)
    print *, 'MAXDIMLIST: ', trim(MAXDIMLIST)
    print *, 'MERGE: ', MERGE
    select case (myHdfVersion)
    case (HDFVERSION_4)
      mls_gfldsetup = swdefgfld(swathid, FIELDName, DIMNAME, Datatype, &
         & MERGE)
    case (HDFVERSION_5)
      if ( chunk_rank /= 0 ) &
        & mls_gfldsetup = HE5_SWdefchunk(swathid, chunk_rank, chunk_dims)
      if ( mls_gfldsetup == 0 ) &
        & mls_gfldsetup = HE5_SWdefgfld(swathid, FIELDName, DIMNAME, MAXDIMLIST, &
        & he2he5_DataType(Datatype), MERGE)
    case default
      mls_gfldsetup = -1
    end select
    print *, 'mls_gfldsetup returns: ', mls_gfldsetup
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

  integer function MLS_SWRDFLD_CHAR_1D ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 1
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    character(len=*), dimension(:), intent(out) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    ! Internal variables
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    ! Declare these as externals to try to fool hdfeos(5)
    integer, external :: swrdfld
    integer, external :: he5_swrdfld
    mls_swrdfld_CHAR_1D = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        mls_swrdfld_CHAR_1D = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWRDFLD_CHAR_1D' )
      endif
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = mls_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    if (myDontFail) then
      if (.not. is_datafield_in_swath(swathid, trim(fieldname), myHdfVersion) )&
        & return
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      if (myDontFail) then
        if (.not. is_swath_datatype_right(swathid, trim(fieldname), &
          & DFNT_CHAR8, myHdfVersion) ) return
      endif
      mls_swrdfld_CHAR_1D = SWRDFLD(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      mls_swrdfld_CHAR_1D = HE5_SWRDFLD(swathid, trim(fieldname), &
        & start, stride, edge, values)
      if (myDontFail) then
        if (.not. is_swath_datatype_right(swathid, trim(fieldname), &
          & HE5T_NATIVE_SCHAR, myHdfVersion) ) return
      endif
    case default
      mls_swrdfld_CHAR_1D = -1
    end select
    if ( myDontFail .or. mls_swrdfld_CHAR_1D /= -1 ) return
    CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to read 1-d char field ' // trim(fieldname) )

  end function MLS_SWRDFLD_CHAR_1D

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

  integer function MLS_SWWRFLD_CHAR_1D ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 1
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    character(len=*), dimension(:), intent(in) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    ! Internal variables
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    ! Declare these as externals to try to fool hdfeos(5)
    integer, external :: swwrfld
    integer, external :: he5_swwrfld
    integer, dimension(12) :: dfrank
    integer, dimension(12) :: numbertype
    character(len=80)      :: fieldlist
    integer                :: nflds
    ! begin execution
    mls_SWWRFLD_CHAR_1D = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        mls_SWWRFLD_CHAR_1D = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWWRFLD_CHAR_1D' )
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
      mls_SWWRFLD_CHAR_1D = SWWRFLD(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      mls_SWWRFLD_CHAR_1D = HE5_SWWRFLD(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case default
      mls_SWWRFLD_CHAR_1D = -1
    end select
    needsFileName = is_swath_datatype_right(swathid, fieldname, DFNT_CHAR8, &
      & myHdfVersion, rank_out=dfrank, numbertype_out=numbertype, &
      & fieldlist_out=fieldlist, nflds_out=nflds)
!>     if(DEEBUG) print *, 'num data fields  ', nflds
!>     if(DEEBUG) print *, 'data field ranks ', dfrank
!>     if(DEEBUG) print *, 'data field types ', numbertype
!>     if(DEEBUG) print *, 'data field list  ', trim(fieldlist)
    if ( myDontFail .or. mls_SWWRFLD_CHAR_1D /= -1 ) return
    CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to write 1-d char field ' // trim(fieldname) )

  end function MLS_SWWRFLD_CHAR_1D

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
          & 'Failed to write 1-d double field ' // trim(fieldname) )

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
          & 'Failed to write 2-d double field ' // trim(fieldname) )

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
          & 'Failed to write 1-d double field ' // trim(fieldname) )

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
          & 'Failed to write 1-d integer field ' // trim(fieldname) )

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
          & 'Failed to write 1-d REAL field ' // trim(fieldname) )

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
    if(DEEBUG) print *, 'swathid: ', swathid
    if(DEEBUG) print *, 'myHdfVersion: ', myHdfVersion
    if(DEEBUG) print *, 'fieldname: ', trim(fieldname)
    if(DEEBUG) print *, 'start: ', start
    if(DEEBUG) print *, 'stride: ', stride
    if(DEEBUG) print *, 'edge: ', edge
    if(DEEBUG) print *, 'shape(values): ', shape(values)
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
          & 'Failed to write 2-d REAL field ' // trim(fieldname) )

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
          & 'Failed to write 1-d REAL field ' // trim(fieldname) )

  end function MLS_SWWRFLD_REAL_3d

  logical function mls_swath_in_file_sca(filename, swath, HdfVersion)
    ! Returns .true. if swath found in file, .false. otherwise
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: swath
    integer, intent(in) :: HdfVersion
    ! Internal variables
    integer                          :: listsize
    character(len=MAXDLISTLENGTH)    :: fieldlist
    integer                          :: nswaths
    ! Begin execution
    nswaths = 0
    mls_swath_in_file_sca = .false.
    select case (HdfVersion)
    case (HDFVERSION_4)
      nswaths = swinqswath(trim(filename), fieldlist, listsize)
    case (HDFVERSION_5)
      nswaths = HE5_swinqswath(trim(filename), fieldlist, listsize)
    end select
    if ( nswaths == 0 ) return
    mls_swath_in_file_sca = &
      & ( StringElementNum(fieldlist, trim(swath), .true.) > 0 )
  end function mls_swath_in_file_sca

  logical function mls_swath_in_file_arr(filename, swaths, HdfVersion, &
    & present )
    ! Array version of the above
    character(len=*), intent(in) :: filename
    character(len=*), dimension(:), intent(in) :: swaths
    integer, intent(in) :: HdfVersion
    logical, dimension(:), intent(out) :: present
    ! Internal variables
    integer                          :: listsize
    character(len=MAXDLISTLENGTH)    :: fieldlist
    integer                          :: nswaths
    integer                          :: i
    ! Begin execution
    nswaths = 0
    mls_swath_in_file_arr = .false.
    present = .false.
    if ( size(swaths) > size(present) ) &
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'array to small to hold values in mls_swath_in_file_arr' )
    select case (HdfVersion)
    case (HDFVERSION_4)
      nswaths = swinqswath(trim(filename), fieldlist, listsize)
    case (HDFVERSION_5)
      nswaths = HE5_swinqswath(trim(filename), fieldlist, listsize)
    end select
    if ( nswaths == 0 ) return
    do i=1, size(swaths)
      present(i) = &
      & ( StringElementNum(fieldlist, trim(swaths(i)), .true.) > 0 )
    enddo
    mls_swath_in_file_arr = any( present )
  end function mls_swath_in_file_arr

! ======================= Private Procedures =========================  

  logical function is_datafield_in_swath(swathid, field, HdfVersion)
    ! Returns .true. if datafield found in swath, .false. otherwise
    integer, intent(in) :: swathid
    character(len=*), intent(in) :: field
    integer, intent(in) :: HdfVersion
    ! Internal variables
    integer, dimension(MAXNODFIELDS) :: rank
    integer, dimension(MAXNODFIELDS) :: numbertype
    character(len=MAXDLISTLENGTH)    :: fieldlist
    integer                         :: nflds
    ! Begin execution
    nflds = 0
    is_datafield_in_swath = .false.
    select case (HdfVersion)
    case (HDFVERSION_4)
      nflds = swinqdflds(swathid, fieldlist, rank, numbertype)
    case (HDFVERSION_5)
      nflds = HE5_swinqdflds(swathid, fieldlist, rank, numbertype)
    end select
    if ( nflds == 0 ) return
    is_datafield_in_swath = &
      & ( StringElementNum(fieldlist, trim(field), .true.) > 0 )
  end function is_datafield_in_swath

  logical function is_swath_datatype_right(swathid, field, datatype, &
    & HdfVersion, rank_out, numbertype_out, fieldlist_out, nflds_out)
    integer, intent(in)                          :: swathid
    character(len=*), intent(in)                 :: field
    integer, intent(in)                          :: datatype
    integer, intent(in)                          :: HdfVersion
    integer, dimension(:), intent(out), optional :: rank_out
    integer, dimension(:), intent(out), optional :: numbertype_out
    character(len=*), intent(out), optional      :: fieldlist_out
    integer, intent(out), optional               :: nflds_out
    ! Internal variables
    integer, dimension(MAXNODFIELDS) :: rank
    integer, dimension(MAXNODFIELDS) :: numbertype
    character(len=MAXDLISTLENGTH)    :: fieldlist
    integer                         :: nflds
    ! Begin execution
    nflds = 0
    is_swath_datatype_right = .false.
    select case (HdfVersion)
    case (HDFVERSION_4)
      nflds = swinqdflds(swathid, fieldlist, rank, numbertype)
    case (HDFVERSION_5)
      nflds = HE5_swinqdflds(swathid, fieldlist, rank, numbertype)
    end select
    if ( present(nflds_out) ) nflds_out = nflds
    if ( present(rank_out) ) rank_out = 0
    if ( present(numbertype_out) ) numbertype_out = 0
    if ( present(fieldlist_out) ) fieldlist_out = ' '
    if ( nflds == 0 ) return
    if ( present(rank_out) ) rank_out = rank(1:size(rank_out))
    if ( present(numbertype_out) ) &
      & numbertype_out = numbertype(1:size(numbertype_out))
    if ( present(fieldlist_out) ) fieldlist_out = fieldlist
    nflds = StringElementNum(fieldlist, trim(field), .true.)
    is_swath_datatype_right = ( nflds > 0 .and. nflds <= MAXNODFIELDS )
    if ( .not. is_swath_datatype_right ) return
    is_swath_datatype_right = &
      & ( numbertype(nflds) == datatype )
  end function is_swath_datatype_right

  ! ---------------------------------------------  he2he5_DataType  -----

  ! This function converts hdfeos2 datatypes to
  ! corresponding hdfeos5 numbers
  ! (unless they are already HE5 types in which it returns them unchanged)

  function he2he5_DataType(dataType) result (HE5_dataType)

    ! Arguments
    integer, intent(IN)       :: dataType
    integer                   :: HE5_dataType
    ! begin
    HE5_dataType = dataType
    select case (dataType)
    case(DFNT_FLOAT32)
      HE5_dataType = HE5T_NATIVE_FLOAT
    case(DFNT_FLOAT64)
      HE5_dataType = HE5T_NATIVE_DOUBLE
    case(DFNT_CHAR8)
      HE5_dataType = HE5T_NATIVE_SCHAR
    case(DFNT_INT8)
      HE5_dataType = HE5T_NATIVE_INT8
    case(DFNT_INT16)
      HE5_dataType = HE5T_NATIVE_INT16
    case(DFNT_INT32)
      HE5_dataType = HE5T_NATIVE_INT
    case(DFNT_INT64)
      HE5_dataType = HE5T_NATIVE_INT64
    !case default
      !HE5_dataType = HE5T_NATIVE_INT
    end select
  end function he2he5_DataType

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MLSHDFEOS

! $Log$
! Revision 2.10  2003/07/09 21:49:07  pwagner
! Wont try swattaching just to see if swath already there
!
! Revision 2.9  2003/07/02 00:54:25  pwagner
! Failed attempt to turn off hdfeos5 warnings by turning off hdf5 ones
!
! Revision 2.8  2003/06/26 00:05:40  pwagner
! Added optional DONTFAIL arg to MLS_SWATTACH
!
! Revision 2.7  2003/06/20 19:31:39  pwagner
! Changes to allow direct writing of products
!
! Revision 2.6  2003/06/06 22:49:12  pwagner
! Added mls_sw(gd)create
!
! Revision 2.5  2003/04/21 19:32:27  pwagner
! Can read/write 1-d char fields
!
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
