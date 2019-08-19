! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module readGriddedUtils ! Collection of subroutines to read TYPE GriddedData_T

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test, Bytes, &
    & Test_Allocate
  use, Intrinsic :: ISO_C_Binding, only: C_Intptr_T, C_Loc
  use Dump_0, Only : Dump
  use GriddedData, only: GriddedData_T, RGR, V_Is_Altitude, V_Is_GPH, &
    & V_Is_Pressure, V_Is_Theta, &
    & AddGriddedDataToDatabase, Dump, SetupNewGriddedData
  use HDFeos, only: HDFe_NentDim, &
    & GDOpen, GDAttach, GDDetach, GDClose, GDFldinfo, &
    & GDInqGrid, GDNentries, GDInqDims, GDInqFlds
  use HDF, only: Dfacc_Create, Dfacc_Rdonly, Dfacc_Rdwr, &
    & Dfnt_Float32, Dfnt_Float64
  use HighOutput, only: OutputNamedValue
  use L3ASCII, only: L3ASCII_Read_Field
  use Lexer_Core, only: Print_Source
  use MLSCommon, only: LineLen, NameLen, FileNameLen, &
    & UndefinedValue, MLSFile_T
  use MLSFiles, only: FileNotFound, Dump, &
    & GetPCFromRef, MLS_Exists, MLS_OpenFile, MLS_CloseFile, &
    & Split_Path_Name, MLS_OpenFile, MLS_CloseFile
  use MLSKinds, only: R4, R8
  use MLSMessageModule, only: MLSMsg_Error, MLSMsg_Info, MLSMsg_Warning, &
    & MLSMessage
  use MLSStrings, only: Capitalize, HHMMSS_Value, Lowercase
  use MLSStringLists, only: GetStringElement, NumStringElements, &
    & List2array, ReplaceSubstring, StringElementNum, SwitchDetail
  use Output_M, only: Output
  use SDPtoolkit, only: PGS_S_Success, &
    & PGS_Io_Gen_Closef, PGS_Io_Gen_Openf, PGSd_Io_Gen_Rseqfrm, &
    & PGSd_Gct_Inverse, &
    & UseSDPtoolkit
  use Toggles, only: Switches
  use Tree, only: Dump_Tree_Node, Where

  implicit none
  private

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (subroutines and functions)
! Read_Climatology       read l3ascii-formatted climatology file
! Read_dao               read meteorology file in one of supported desciptions
! Read_geos5_or_merra    read a geos5 or merra meteorology file
! Read_geos5_7           read netCDF4-formatted meteorology file
! Read_merra_2           read netCDF4-formatted meteorology merra file
! Read_ncep_gdas         read meteorology file in one of supported desciptions
! Read_ncep_strat        read meteorology file in one of supported desciptions
! ReadGloriaFile         read binary-formatted file designed by G. Manney
! Write_Merra            write this kind of meteorology file
!
! === (end of toc) ===

  public:: Announce_error   
  public:: Read_Climatology   
  public:: Read_dao           
  public:: Read_geos5_or_merra
  public:: Read_geos5_7       
  public:: Read_merra_2     
  public:: Read_ncep_gdas     
  public:: Read_ncep_strat    
  public:: ReadGloriaFile     
  public:: Write_Merra     

  integer :: ERROR
  character(len=8), public :: LIT_DESCRIPTION = 'none'
  ! First we'll define some global parameters and data types.
  logical, parameter :: COUNTEMPTY=.true.
  character (len=*), parameter :: DEFAULTDAODIMLIST = 'XDim,YDim,Height,Time'
  character (len=*), parameter :: DEFAULTDAOFIELDNAME = 'TMPU'
  character (len=*), parameter :: DEFAULTGEOS5FIELDNAME = 'T'
  character (len=*), parameter :: DEFAULTNCEPGRIDNAME = 'TMP_3'
  character (len=*), parameter :: DEFAULTNCEPSTRATFIELDNAME = 'Temperature'

  character (len=*), parameter :: lit_dao = 'dao'
  character (len=*), parameter :: lit_ncep = 'ncep'
  character (len=*), parameter :: lit_strat = 'strat'
  character (len=*), parameter :: lit_geos5 = 'geos5'
  integer, parameter :: MAXLISTLENGTH=Linelen ! Max length list of grid names
  integer, parameter :: NENTRIESMAX=200 ! Max num of entries
  integer, parameter :: MAXDS = 1024 ! 500
  integer, parameter :: MAXSDNAMESBUFSIZE = MAXDS*NAMELEN
contains

  ! ------------------------------------------------  announce_error  -----
  subroutine announce_error ( lcf_where, full_message, &
    & extra_message, extra_number, use_toolkit, non_fatal )

    ! Arguments
    integer, intent(in)    :: lcf_where
    character(len=*), intent(in)    :: full_message
    logical, intent(in), optional :: use_toolkit
    character(len=*), intent(in), optional    :: extra_message
    integer, intent(in), optional    :: extra_number
    logical, intent(in), optional :: non_fatal

    ! Local variables
    logical :: just_print_it
    logical :: its_non_fatal
    logical, parameter :: default_output_by_toolkit = .true.
    logical :: verbose
    ! Executable
    its_non_fatal = .false.
    if(present(use_toolkit)) then
      just_print_it = use_toolkit
    else if(default_output_by_toolkit) then
      just_print_it = .false.
    else
      just_print_it = .true.
    end if
    verbose = ( SwitchDetail(switches, 'apr') > -1 )
    its_non_fatal = just_print_it
    if ( present(non_fatal) ) its_non_fatal = non_fatal

    if( .not. its_non_fatal ) then
      error = max(error,1)
      call output ( '***** At ' )

      if(lcf_where > 0) then
        call print_source ( where(lcf_where) )
      else
        call output ( '(no lcf node available)' )
      end if

      if ( verbose ) then
        call output ( ': ' )
        call output ( "The " );
        if ( lcf_where > 0 ) then
          call dump_tree_node ( lcf_where, 0 )
        else
          call output ( '(no lcf tree available)' )
        end if
      end if

      call output(" Caused the following error:", advance='yes', &
        & from_where=ModuleName)
      call output(trim(full_message), advance='yes', &
        & from_where=ModuleName)
      if(present(extra_message)) then
        call output('error number ', advance='no')
        call output(extra_message, advance='yes')
      end if
      if(present(extra_number)) then
        call output('error number ', advance='no')
        call output(extra_number, places=9, advance='yes')
      end if
    else
      call output ( '***Non-fatal error in module ' )
      call output ( ModuleName, advance='yes' )
      call output ( trim(full_message), advance='yes' )
      if ( present(extra_message) ) then
        call output ( 'Error number ' )
        call output ( extra_message, advance='yes' )
      end if
      if ( present(extra_number) ) then
        call output ( 'Error number ' )
        call output ( extra_number, advance='yes' )
      end if
    end if

  end subroutine announce_error

  ! ----------------------------------------------- Read_GEOS5_7
  subroutine Read_GEOS5_7( GEOS5File, lcf_where, v_type, &
    & the_g_data, GeoDimList, fieldName, date, sumDelp )

    use Dates_module, only: UTC2TAI93s
    use Dump_1, only: Dump
    use HDF5, only: Hsize_t
    use MLSHDF5, only: GetAllHDF5DSNames, &
      & GetHDF5DSRank, gethDF5DSDims, LoadFromHDF5DS, &
      & ReadHDF5Attribute

    ! This routine reads a gmao geos5_7 file, named something like
    ! DAS.ops.asm.avg3_3d_Nv.GEOS571.20110930_0300.V01.nc4 (pressure with

    ! This file is formatted in the following way:
    ! each gridded quantity, not an hdfeos grrid, by the way, is
    ! a rank 4 field, e.g. T, is given with dimensions
    ! 'Time, Height, YDim, XDim'
    ! We'll simply copy it into a single gridded data type

    ! Arguments
    type(MLSFile_T)                :: GEOS5file
    integer, intent(IN) :: lcf_where    ! node of the lcf that provoked me
    integer, intent(IN) :: v_type       ! vertical coordinate; an 'enumerated' type
    type( GriddedData_T ) :: the_g_data ! Result
    character (len=*), optional, intent(IN) :: GeoDimList ! Comma-delimited dim names
    character (len=*), optional, intent(IN) :: fieldName ! Name of gridded field
    character (len=*), optional, intent(IN) :: date ! date (offset)
    logical, optional, intent(IN)           :: sumDelp ! sum the DELP to make PL
    ! Local variables
    character (len=NAMELEN) :: actual_field_name
    integer(c_intptr_t) :: Addr         ! For tracing
    character(len=19) :: datestring ! will be in the form of yyyy-MM-ddTHH:MM:ss
    logical :: DEEBUG
    integer(kind=hsize_t) :: dim1(1), dim4(4)
    integer :: error, rank
    character(len=256) :: errormsg
    integer, dimension(4) :: idim4
    character (len=MAXSDNAMESBUFSIZE) :: mySdList
    integer :: S                     ! Size in bytes of an object to deallocate
    real(r8), dimension(:), pointer :: temp1d => null()
    real(r8), dimension(:,:,:,:), pointer :: temp4d => null()
    character(len=16) :: the_units
    real(r4), parameter :: FILLVALUE = 1.e15 !this value might be wrong
    integer :: mydate, mytime, year, month, day, hour, minute, second
    logical :: verbose

    ! Executable
    if(present(fieldName)) then
      actual_field_name=fieldName
    else
      actual_field_name=DEFAULTGEOS5FIELDNAME
    end if
    nullify( temp1d, temp4d )
    DEEBUG = ( index(lowercase(actual_field_name), 'inq') > 0 )
    verbose = ( switchDetail(switches, 'geos5') > -1 ) .or. DEEBUG
    call GetAllHDF5DSNames ( GEOS5File%Name, '/', mysdList )
    if ( verbose ) call dump( mysdList, 'DS names' )

    ! Fill in value
    the_g_data%quantityName = actual_field_name
    the_g_data%description = lit_geos5
    the_g_data%verticalCoordinate = v_type
    the_g_data%nodates = 1
    the_g_data%empty = .true.
    the_g_data%missingvalue = FILLVALUE ! this value might be wrong
    call allocate_test ( the_g_data%lsts, 1, "the_g_data%lsts", moduleName )
    the_g_data%lsts = the_g_data%missingvalue
    the_g_data%nolsts = 1
    the_g_data%noszas = 1
    the_g_data%lsts = the_g_data%missingValue ! Know how to read this yet?
    call allocate_test ( the_g_data%szas, 1, "the_g_data%szas", moduleName )
    the_g_data%szas = the_g_data%missingvalue

    call mls_openfile(geos5File, error)
    if (error .gt. 0) call announce_error(lcf_where, "Could not open "// GEOS5File%Name)

    ! Get lons
    call GetHDF5DSRank (geos5file%fileid%f_id, 'lon', rank)
    if (rank /= 1) call announce_error(lcf_where, "lon must be a 1-dim array: "// GEOS5File%Name)

    call GetHDF5DSDims(geos5file%fileid%f_id, 'lon', dim1)
    the_g_data%noLons = dim1(1)

    call allocate_test ( temp1d, the_g_data%noLons, "temp1d", moduleName )

    call LoadFromHDF5DS (geos5file%fileid%f_id, 'lon', temp1d)
    call allocate_test ( the_g_data%lons, the_g_data%nolons, "the_g_data%nolons", &
      & moduleName )
    the_g_data%lons = temp1d
    s = size(temp1d) * storage_size(temp1d) / 8
    call deallocate_test ( temp1d, "temp1d", moduleName )

    ! Fill dateStart and dateEnd
    mytime = 0. ! In case it's not found
    if (.not. ReadHDF5Attribute(geos5file%fileid%f_id, 'time', &
        & 'begin_time', mytime, error=errormsg)) then
        call announce_error (lcf_where, errormsg // ' in file ' &
        & // geos5file%name)
    end if

    if ( mytime < 0 ) then
      call announce_error (lcf_where, "Invalid 'begin_time' value in " // geos5file%name)
    end if

    mydate = 10000*2001 + 1*100 + 1 ! In case it's not found
    if (.not. ReadHDF5Attribute(geos5file%fileid%f_id, 'time', &
        & 'begin_date', mydate, error=errormsg)) then
        call announce_error (lcf_where, errormsg // ' in file ' &
        & // geos5file%name)
    end if

    if (mydate < 0) then
        call announce_error (lcf_where, "Invalid 'begin_date' value in " // geos5file%name)
    end if

    year = mydate / 10000
    month = mod(mydate, 10000) / 100
    day = mod(mydate, 100)
    if (year > 9999 .or. month > 12 .or. month < 1 .or. day > 31 .or. day < 1) then
        call announce_error (lcf_where, "Invalid 'begin_date' value in " // geos5file%name)
    end if    
    
    hour = mytime / 10000
    minute = mod (mytime, 10000) / 100
    second = mod(mytime, 100)
    if (hour > 24 .or. minute > 60 .or. second > 60) then
        call announce_error (lcf_where, "Invalid 'begin_time' value in " // geos5file%name)
    end if

    write ( datestring, &
      & '(I4.4, A1,   I2.2,  A1, I2.2,  A1, I2.2, A1,  I2.2,    A1, I2.2)') &
      &   year, '-',  month, '-', day, 'T', hour, ':', minute, ':', second

    call allocate_test ( the_g_data%datestarts, the_g_data%nodates, &
      & "the_g_data%datestarts", moduleName )
    call allocate_test ( the_g_data%dateends, the_g_data%nodates, &
      & "the_g_data%dateends", moduleName )

    the_g_data%datestarts(1) = utc2tai93s(datestring)

    the_g_data%dateends(1) = the_g_data%datestarts(1) ! start and end are the same

    ! Get lats
    call GetHDF5DSRank (geos5file%fileid%f_id, 'lat', rank)
    if (rank /= 1) call announce_error(lcf_where, "lat must be a 1-dim array: " // geos5file%name)

    call GetHDF5DSDims(geos5file%fileid%f_id, 'lat', dim1)
    the_g_data%noLats = dim1(1)

    call allocate_test ( temp1d, the_g_data%nolats, "temp1d", moduleName )

    call LoadFromHDF5DS (geos5file%fileid%f_id, 'lat', temp1d)
    call allocate_test ( the_g_data%lats, the_g_data%nolats, "the_g_data%lats", &
      & moduleName )
    the_g_data%lats = temp1d
    call deallocate_test ( temp1d, "temp1d", moduleName )

    ! Get heights
    call GetHDF5DSRank (geos5file%fileid%f_id, 'lev', rank)
    if (rank /= 1) call announce_error(lcf_where, "lev must be a 1-dim array: " // geos5file%name)
    
    call GetHDF5DSDims (geos5file%fileid%f_id, 'lev', dim1)
    the_g_data%noheights = dim1(1)

    call allocate_test ( temp1d, the_g_data%noheights, "temp1d", moduleName )
    
    call LoadFromHDF5DS (geos5file%fileid%f_id, 'lev', temp1d)
    call allocate_test ( the_g_data%heights, the_g_data%noheights, &
      & "the_g_data%heights", moduleName )
    
    ! We cannot load directly into the the_g_data's array because
    ! this is an array of 32-bit float, while the data from file is
    ! 64-bit float
    the_g_data%heights = temp1d
    call deallocate_test ( temp1d, "temp1d", moduleName )

    the_g_data%heightsunits = 'hPa' ! the units is according the file_specification

    ! The following is according to GEOS-5.7.2 file specification
    select case ( lowercase(actual_field_name) )
    case ( 'h' )
        the_units = 'm'
    case ( 'pl', 'ps', 'delp' )
        the_units = 'Pa'
    case ( 't' )
        the_units = 'K'
    case default
      call announce_error( lcf_where, &
        & "Unexpected field name: " // actual_field_name, non_fatal=.true. )
    end select
    the_g_data%units = the_units

    ! read the field, field should be either PL or T
    call GetHDF5DSRank (geos5file%fileid%f_id, capitalize(actual_field_name), rank)
    if (rank /= 4) then
        call announce_error ( lcf_where, &
        trim(capitalize(actual_field_name)) // " must be a 4-dim array: (geos57)" // geos5file%name)
    end if

    call GetHDF5DSDims (geos5file%fileid%f_id, capitalize(actual_field_name), dim4)
    if (dim4(4) /= 1) then
        call announce_error(geos5file%fileid%f_id, &
        "The fourth dimension of " // actual_field_name // " is not 1 in: " // geos5file%name)
    end if
    idim4 = dim4

    call allocate_test ( temp4d, idim4(1), idim4(2), idim4(3), idim4(4), &
        & 'temp4d', ModuleName // 'Read_geos57' )

    call LoadFromHDF5DS (geos5file%fileid%f_id, capitalize(actual_field_name), temp4d)

    ! call allocate_test ( the_g_data%field(:,:,:,1,1,1), idim4(1:3), &
    !     & 'the_g_data%field', ModuleName // 'Read_geos57' )
    allocate ( the_g_data%field(idim4(3), idim4(2), idim4(1), 1, 1, 1), STAT=error)
    addr = 0
    if ( error == 0 ) then
      if ( size(the_g_data%field) > 0 ) addr = transfer(c_loc( &
        & the_g_data%field(1,1,1,1,1,1)), addr)
    end if
    call test_allocate ( error, moduleName, 'the_g_data%field', (/1,1,1,1,1,1/), &
      & (/ idim4(3), idim4(2), idim4(1), 1, 1, 1 /), bytes(the_g_data%field), &
      & address=addr )
    the_g_data%field(:,:,:,1,1,1) = reshape(temp4d(:,:,:,1), order=(/3,2,1/), &
    shape=(/the_g_data%noheights, the_g_data%nolats, the_g_data%nolons/))
    !  do i1 = 1, idim4(1)
    !    do i2 = 1, idim4(2)
    !      do i3 = 1, idim4(3)
    !        the_g_data%field(i3,i2,i1,1,1,1) = temp4d(i1,i2,i3,1)
    !      end do
    !    end do
    !  end do
    call deallocate_test ( temp4d, 'temp4d', ModuleName // 'Read_geos57' )
    
    ! Read file successful
    the_g_data%empty = .false.
    call mls_closefile(geos5File, error)
    if (error .gt. 0) call announce_error(lcf_where, "Could not close "// GEOS5File%Name)
    if ( verbose ) call dump( the_g_data, details = 1 )
  end subroutine Read_GEOS5_7

  ! ----------------------------------------------- Read_GEOS5_or_Merra
  subroutine Read_GEOS5_or_Merra( GEOS5File, lcf_where, v_type, &
    & the_g_data, GeoDimList, fieldName, date, sumDelp )

    ! This routine reads a gmao GEOS5 or MERRA file, named something like
    ! DAS.ops.asm.tavg3d_prs_v.GEOS500.20060320_0000.V01 (pressure with
    ! fieldname='PL') or 
    ! DAS.ops.asm.tavg3d_dyn_v.GEOS500.20060320_0000.V01 (temperature
    ! with fieldname='T') or
    ! MERRA300.prod.assim.inst6_3d_ana_Nv.20050131.hdf (pressure with
    ! fieldname='DELP' or 'T')
    ! returning a filled data
    ! structure appropriate for older geos5 or newer style gmao merra

    ! FileName and the_g_data are required args
    ! GeoDimList, if present, should be the Dimensions' short names
    ! as a comma-delimited character string in the order:
    ! longitude, latitude, vertical level, time

    ! fieldName, if present, should be the rank 4
    ! like temperature

    ! This file is formatted in the following way:
    ! At each gridded quantity, e.g. EOSGRID,
    ! a rank 4 field, e.g. T, is given with dimensions
    ! 'Time, Height, YDim, XDim'
    ! We'll simply copy it into a single gridded data type
    
    ! Arguments
    type(MLSFile_T)                :: GEOS5file
    integer, intent(IN) :: lcf_where    ! node of the lcf that provoked me
    integer, intent(IN) :: v_type       ! vertical coordinate; an 'enumerated' type
    type( GriddedData_T ) :: the_g_data ! Result
    character (len=*), optional, intent(IN) :: GeoDimList ! Comma-delimited dim names
    character (len=*), optional, intent(IN) :: fieldName ! Name of gridded field
    character (len=*), optional, intent(IN) :: date ! date (offset)
    logical, optional, intent(IN)           :: sumDelp ! sum the DELP to make PL

    ! Local Variables
    integer :: itime, ilat, ilon, ilev
    integer :: file_id, gd_id
    integer :: inq_success
    integer :: nentries, ngrids, ndims, nfields
    integer :: strbufsize
    logical,  parameter       :: CASESENSITIVE = .false.
    integer, parameter :: GRIDORDER=1   ! What order grid written to file
    character (len=MAXLISTLENGTH) :: gridlist
    character (len=MAXLISTLENGTH) :: dimlist, actual_dim_list
    character (len=MAXLISTLENGTH), dimension(1) :: dimlists
    character (len=16), dimension(NENTRIESMAX) :: dimNames
    character (len=MAXLISTLENGTH) :: fieldlist
    integer, parameter :: MAXNAMELENGTH=NameLen         ! Max length of grid name
    character (len=MAXNAMELENGTH) :: gridname, actual_field_name
    integer, dimension(NENTRIESMAX) :: dims, rank, numberTypes
    integer                        :: our_rank, numberType
    logical                        :: mySum

    integer :: start(4), stride(4), edge(4)
    integer :: status
    character(len=16) :: the_units
    integer                        :: timeIndex
    !                                  These start out initialized to one
    integer                        :: nlon=1, nlat=1, nlev=1, ntime=1
    integer, external :: GDRDFLD
    real(r4), parameter :: FILLVALUE = 1.e15
    real(r4), dimension(:,:,:,:), pointer :: all_the_fields => null()
    real(r8), dimension(:), pointer :: dim_field, pb
    real(r8) :: dateValue
    logical :: DEEBUG
    ! Executable code
    dateValue = 0.
    if ( present(date) ) dateValue = HHMMSS_value( date, status )
    if(present(fieldName)) then
      actual_field_name=fieldName
    else
      actual_field_name=DEFAULTGEOS5FIELDNAME
    end if
    
    DEEBUG = ( index(lowercase(actual_field_name), 'inq') > 0 )
    
    ! We need the timeIndex as an offset for merra files
    timeIndex = int(1.5 + (dateValue - 1.)/(6*60*60)) ! How many quarter-days?
    if ( DEEBUG ) then
      print *, 'date (s) ', dateValue
      print *, 'time index ', timeIndex
    end if
    mySum = .false.
    if ( present(sumDelp) ) mySum = sumDelp
    ! Find list of grid names on this file (This has been core dumping on me)
    if(DEEBUG) print *, 'About to find grid list of file ', trim(GEOS5File%Name)
    gridlist = ''
    inq_success = gdinqgrid(GEOS5File%Name, gridlist, strbufsize)
    if (inq_success < 0) then
      call announce_error(lcf_where, "Could not inquire gridlist "// &
        & trim(GEOS5File%Name))
    else if ( strbufsize > MAXLISTLENGTH ) then
       CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'list size too big in Read_GEOS5_or_Merra ' // trim(GEOS5File%Name), MLSFile=GEOS5File )
    else if ( strbufsize < MAXLISTLENGTH .and. strbufsize > 0 ) then
      gridlist = gridlist(1:strbufsize) // ' '
    end if
    if(DEEBUG) print *, 'grid list ', trim(gridlist)

    error = 0
    file_id = gdopen(GEOS5File%Name, DFACC_RDONLY)

    if (file_id < 0) then
      call announce_error(lcf_where, "Could not open "// GEOS5File%Name)
    end if

    ! Find grid name corresponding to the GRIDORDER'th one
    ngrids = NumStringElements(gridlist, COUNTEMPTY)

    if(ngrids <= 0) then
      call announce_error(lcf_where, "NumStringElements of gridlist <= 0")
    else if(ngrids /= inq_success) then
      call announce_error(lcf_where, "NumStringElements of gridlist /= inq_success")
    else if(ngrids < GRIDORDER) then
      call announce_error(lcf_where, "NumStringElements of gridlist < GRIDORDER")
    end if

    call GetStringElement(gridlist, gridname, GRIDORDER, COUNTEMPTY)

    gd_id = gdattach(file_id, gridname)
    if (gd_id < 0) then
      call announce_error(lcf_where, "Could not attach "//trim(gridname))
    end if

    ! Now find dimsize(), dimname(), etc.
    nentries = gdnentries(gd_id, HDFE_NENTDIM, strbufsize)

    if(nentries <= 0) then
      call announce_error(lcf_where, "nentries of gd_id <= 0")
    else if(nentries > NENTRIESMAX) then
      call announce_error(lcf_where, "nentries of gd_id > NENTRIESMAX")
    end if

    dimlist = ''
    ndims = gdinqdims(gd_id, dimlist, dims)

    if(ndims <= 0) then
      call announce_error(lcf_where, "ndims of gd_id <= 0")
    else if(ndims > NENTRIESMAX) then
      call announce_error(lcf_where, "ndims of gd_id > NENTRIESMAX")
    end if

    fieldlist = ''
    nfields = gdinqflds(gd_id, fieldlist, rank, numberTypes)

    if(nfields <= 0) then
      call announce_error(lcf_where, "nfields of gd_id <= 0")
    else if(nfields > NENTRIESMAX) then
      call announce_error(lcf_where, "nfields of gd_id > NENTRIESMAX")
    end if

    if(.not. CASESENSITIVE) then
      fieldlist = Capitalize(fieldlist)
    end if

    if(DEEBUG) print *, 'nentries ', nentries
    if(DEEBUG) print *, 'ndims ', ndims
    if(DEEBUG) print *, 'dimlist ', dimlist
    if(DEEBUG) print *, 'nfields ', nfields
    if(DEEBUG) print *, 'fieldlist ', fieldlist
    if(DEEBUG) print *, 'actual_field_name ', actual_field_name

    actual_dim_list = ' '
    if(present(GeoDimList)) then
      actual_dim_list=GeoDimList
    end if
    if ( actual_dim_list == ' ' ) then
      actual_dim_list = DEFAULTDAODIMLIST
    end if
    call List2Array (actual_dim_list, dimNames, countEmpty)

    ! Check that our requested field is present
    ! We might, after all, have been doing an 'inq'
    if ( &
      & StringElementNum( fieldlist, trim(actual_field_name), COUNTEMPTY )&
      &  < 1 ) then
      the_g_data%empty = .true.
      call output( trim(actual_field_name) // ' not found in ' // &
        & trim(GEOS5File%Name), advance='yes' )
      call outputNamedValue( 'fieldList', trim(fieldlist) )
      status = gddetach(gd_id)
      if(status /= 0) &
        & call announce_error(lcf_where, "failed to detach from grid " &
        & //trim(gridname))
      status = gdclose(file_id)
      if(status /= 0) &
        & call announce_error(lcf_where, "failed to close file " //trim(GEOS5File%Name))
        return
    end if
    ! Now find the rank of our field
    inq_success = gdfldinfo(gd_id, trim(actual_field_name), our_rank, dims, &
      & numbertype, dimlists(1))

    if ( inq_success /= PGS_S_SUCCESS ) then
      the_g_data%empty = .true.
      CALL MLSMessage ( MLSMSG_Warning, moduleName,  &
        & trim(actual_field_name) // ' was inq_fail ' // trim(GEOS5File%Name) )
      return
    end if
    dimlist = trim(dimlists(1))
    if(DEEBUG) print *, 'our_rank ', our_rank
    if(DEEBUG) print *, 'dims ', dims(1:our_rank)
    if(DEEBUG) print *, 'dimlist ', dimlist

    nlon = dims(1)
    nlat = dims(2)
    nlev = dims(3)
    ntime = dims(4)

    the_g_data%quantityName = actual_field_name
    the_g_data%description = lit_geos5
    the_g_data%verticalCoordinate = v_type

    the_g_data%noLons = nlon
    the_g_data%noLats = nlat
    the_g_data%noHeights = nlev
    ! the_g_data%noLsts = ntime
    if ( LIT_DESCRIPTION == lit_geos5 ) then
      the_g_data%noDates = ntime
      if ( ntime > 1 ) then
        the_g_data%empty = .true.
        CALL MLSMessage ( MLSMSG_Info, moduleName,  &
          & 'Was apparently not geos5 ' // trim(GEOS5File%Name) )
        return
      end if
    else
      the_g_data%noDates = 1
    end if
    ! The following is an awful hack
    ! to prevent me from the having to read the units attribute
    select case ( lowercase(actual_field_name) )
    case ( 'pl', 'delp' )
      the_units = 'Pa'
    case ( 't' )
      the_units = 'K'
    case default
      the_units = 'Pa'
    end select
    if(DEEBUG) print *, 'our quantity name ', the_g_data%quantityName
    if(DEEBUG) print *, 'our description ', the_g_data%description
    if(DEEBUG) print *, 'our units ', the_units
    if(DEEBUG) print *, 'our vertical coord ', the_g_data%verticalCoordinate
    if(DEEBUG) print *, 'v_type ', v_type

    ! Setup the grid
    call SetupNewGriddedData ( the_g_data, noHeights=nlev, noLats=nlat, &
      & noLons=nlon, noLsts=1, noSzas=1, noDates=the_g_data%noDates, &
      & missingValue=FILLVALUE, units=the_units, verticalCoordinate=v_type, &
      & heightsunits='hPa' )
    if(DEEBUG) print *, '(Again) our quantity name ', the_g_data%quantityName
    if(DEEBUG) print *, 'our description ', the_g_data%description
    if(DEEBUG) print *, 'our units ', the_g_data%units
    if(DEEBUG) print *, 'our vertical coord ', the_g_data%verticalCoordinate
    if(DEEBUG) print *, 'v_type ', v_type
    if(DEEBUG) print *, 'About to allocate ', dims(1:4)
    call allocate_test( all_the_fields, dims(1), dims(2), dims(3), dims(4), &
        & 'all_the_fields', ModuleName // 'Read_GEOS5_or_Merra' )
    all_the_fields = the_g_data%missingValue
    start = 0                                                             
    stride = 1                                                               
    edge = dims(1:4)                                                        
    if(DEEBUG) print *, 'About to read ' // trim(actual_field_name)
    if(DEEBUG) print *, 'Start ', Start
    if(DEEBUG) print *, 'Stride ', Stride
    if(DEEBUG) print *, 'Edge ', Edge
    status = gdrdfld(gd_id, trim(actual_field_name), start, stride, edge, &  
      & all_the_fields)                                                          
    if(status /= 0) &
      & call announce_error(lcf_where, "failed to read field " &
      & //trim(actual_field_name))
    if ( mySum ) then
      ! If we are so asked
      ! Sum up all the delps, starting from the top of the atmosphere
      nullify(pb)
      call allocate_test( pb, nlev+1, 'pb', &
        & ModuleName // 'Read_GEOS5_or_Merra' )
      do itime=1, ntime
        do ilat=1, nlat
          do ilon=1, nlon
            pb = 0.0
            ! pb(nlev+1) = 1.0 ! We must assume the topmost level is 1 hPa
            ! do ilev=nlev, 1, -1
            !  pb(ilev) = pb(ilev+1) + all_the_fields(ilon, ilat, ilev, itime)
            pb(1) = 1.0 ! We must assume the topmost level is 1 hPa
            do ilev=1, nlev
              pb(ilev+1) = pb(ilev) + all_the_fields(ilon, ilat, ilev, itime)
            end do
            do ilev=1, nlev
              all_the_fields(ilon, ilat, ilev, itime) = &
                & 0.5*( pb(ilev) + pb(ilev+1) )
            end do
          end do
        end do
      end do
      call deallocate_test( pb, 'pb', &
        & ModuleName // 'Read_GEOS5_or_Merra' )
    end if
    ! The actual dimlist is this                    XDim,YDim,Height,TIME
    ! Need to reshape it so that the order becomes: Height,YDim,XDim,TIME
    if ( DEEBUG) then
      print *, LIT_DESCRIPTION // ' Before reshaping'
      call dump(all_the_fields(:,1,1,1), 'x-slice')
      call dump(all_the_fields(1,:,1,1), 'y-slice')
      call dump(all_the_fields(1,1,:,1), 'p-slice')
      call dump(all_the_fields(1,1,1,:), 't-slice')
    end if
    if ( LIT_DESCRIPTION == lit_geos5 ) then
      ! GEOS5 format
      the_g_data%field(:,:,:,1,1,:) = reshape( all_the_fields, &
        & shape=(/nLev, nlat, nlon, ntime/), order=(/3,2,1,4/) &
        & )
    else if ( size(all_the_fields, 4) < 2 ) then
      ! Claimed to be MERRA (but array size says otherwise)
      the_g_data%empty = .true.
      call output( trim(actual_field_name) // ' had the wrong size (field, 4) ' // &
        & trim(GEOS5File%Name), advance='yes' )
      call deallocate_test( all_the_fields, 'all_the_fields', &
        & ModuleName // 'Read_GEOS5_or_Merra' )
      deallocate(the_g_data%field, STAT=status)
      if ( status /= 0 ) &
        & call announce_error( lcf_where, "failed to deallocate field" )
      return
    else
      ! MERRA format
      the_g_data%field(:,:,:,1,1,1) = reshape( all_the_fields(:,:,:,timeIndex), &
        & shape=(/nLev, nlat, nlon/), order=(/3,2,1/) &
        & )
    end if
    if ( DEEBUG ) then
      print *, LIT_DESCRIPTION // ' After reshaping'
      call dump(the_g_data%field(1,1,:,1,1,1), 'x-slice')
      call dump(the_g_data%field(1,:,1,1,1,1), 'y-slice')
      call dump(the_g_data%field(:,1,1,1,1,1), 'p-slice')
      call dump(the_g_data%field(1,1,1,:,1,1), 't-slice')
    end if
    call deallocate_test( all_the_fields, 'all_the_fields', &
        & ModuleName // 'Read_GEOS5_or_Merra' )
    ! Now read the dims
    nullify(dim_field)
    ! call read_the_dim(gd_id, 'XDim', dims(1), dim_field)
    call read_the_dim(gd_id, trim(dimNames(1)), dims(1), dim_field)
    the_g_data%lons = dim_field
    ! call read_the_dim(gd_id, 'YDim', dims(2), dim_field)
    call read_the_dim(gd_id, trim(dimNames(2)), dims(2), dim_field)
    the_g_data%lats = dim_field
    ! call read_the_dim(gd_id, 'Height', dims(3), dim_field)
    call read_the_dim(gd_id, trim(dimNames(3)), dims(3), dim_field)
    the_g_data%Heights = dim_field
    ! call read_the_dim(gd_id, 'Time', dims(4), dim_field)
    call read_the_dim(gd_id, trim(dimNames(4)), dims(4), dim_field)
    ! the_g_data%lsts = dim_field
    ! Because GMAO background files have no date field
    if ( all(dim_field == 0._r8) ) then
      the_g_data%dateStarts = dateValue
      the_g_data%dateEnds = dateValue
    else if (LIT_DESCRIPTION == lit_geos5 ) then
      the_g_data%dateStarts = dim_field
      the_g_data%dateEnds = dim_field
    else
      the_g_data%dateStarts = dim_field(timeIndex)
      the_g_data%dateEnds = dim_field(timeIndex)
    end if
    if ( DEEBUG ) then
      call dump(the_g_data%dateStarts, 'dateStarts')
      call dump(the_g_data%dateEnds, 'dateEnds')
    end if
    call deallocate_test ( dim_field, "dim_field", moduleName )
    ! Have not yet figured out how to assign these
    ! Probably will have to read metadata
    the_g_data%Szas = the_g_data%missingValue
    ! the_g_data%DateStarts = the_g_data%missingValue
    ! the_g_data%DateEnds = the_g_data%missingValue
    the_g_data%lsts = the_g_data%missingValue
    ! Close grid
    status = gddetach(gd_id)
    if(status /= 0) &
      & call announce_error(lcf_where, "failed to detach from grid " &
      & //trim(gridname))
    status = gdclose(file_id)
    if(status /= 0) &
      & call announce_error(lcf_where, "failed to close file " //trim(GEOS5File%Name))
  contains 
    subroutine read_the_dim(gd_id, field_name, field_size, values)
      ! Arguments
      integer, intent(in) :: gd_id
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: field_size
      real(r8), dimension(:), pointer :: values
      ! Local variables
      integer :: status
      integer, dimension(1) :: start, stride, edge
      ! Executable
      call allocate_test ( values, field_size, "values", moduleName )
      start = 0                                                             
      stride = 1                                                             
      edge = field_size                                                       
      values = the_g_data%missingValue
      status = gdrdfld(gd_id, trim(field_name), start, stride, edge, &
        & values)                                                     
      if ( status /= 0 ) &
        & call announce_error(lcf_where, "failed to read dim field " &
        & // trim(field_name))
    end subroutine read_the_dim
  end subroutine Read_GEOS5_or_Merra

  ! ---------------------------------------------------  Read_DAO  -----
  subroutine Read_DAO(DAOFile, lcf_where, v_type, &
    & the_g_data, GeoDimList, fieldName)


    ! What was once called 'dao' now better known as 'gmao'
    ! This routine reads a gmao file, named something like
    ! DAS.flk.asm.tsyn3d_mis_p.GEOS401.2003011800.2003011818.V01
    ! returning a filled data
    ! structure appropriate for newer style dao
    ! (And just how does it differ from the older style? 
    !  For one thing it holds 4-d, not 3-d fields)

    ! FileName and the_g_data are required args
    ! GeoDimList, if present, should be the Dimensions' short names
    ! as a comma-delimited character string in the order:
    ! longitude, latitude, vertical level, time

    ! fieldName, if present, should be the rank 4
    ! like temperature

    ! This file is formatted in the following way:
    ! At each gridded quantity, e.g. EOSGRID,
    ! a rank 4 field, e.g. TMPU, is given with dimensions
    ! 'Time, Height, YDim, XDim'
    ! We'll simply copy it into a single gridded data type

    ! Arguments
    type(MLSFile_T)                :: DAOFile
    integer, intent(IN) :: lcf_where    ! node of the lcf that provoked me
    integer, intent(IN) :: v_type       ! vertical coordinate; an 'enumerated' type
    type( GriddedData_T ) :: the_g_data ! Result
    character (len=*), optional, intent(IN) :: GeoDimList ! Comma-delimited dim names
    character (len=*), optional, intent(IN) :: fieldName ! Name of gridded field

    ! Local Variables
    integer :: file_id, gd_id
    integer :: inq_success
    integer :: nentries, ngrids, ndims, nfields
    integer :: strbufsize
    logical,  parameter       :: CASESENSITIVE = .false.
    integer, parameter :: GRIDORDER=1   ! What order grid written to file
    character (len=MAXLISTLENGTH) :: gridlist
    character (len=MAXLISTLENGTH) :: dimlist, actual_dim_list
    character (len=MAXLISTLENGTH), dimension(1) :: dimlists
    character (len=16), dimension(NENTRIESMAX) :: dimNames
    character (len=MAXLISTLENGTH) :: fieldlist
    integer, parameter :: MAXNAMELENGTH=NameLen         ! Max length of grid name
    character (len=MAXNAMELENGTH) :: gridname, actual_field_name
    integer, dimension(NENTRIESMAX) :: dims, rank, numberTypes
    integer                        :: our_rank, numberType

    integer :: start(4), stride(4), edge(4)
    integer :: status
    !                                  These start out initialized to one
    integer                        :: nlon=1, nlat=1, nlev=1, ntime=1
    integer, external :: GDRDFLD
    real(r4), parameter :: FILLVALUE = 1.e15
    real(r4), dimension(:,:,:,:), pointer :: all_the_fields
    real(r8), dimension(:), pointer :: dim_field
    logical, parameter :: DEEBUG = .false.
    ! Executable code
    ! Find list of grid names on this file (This has been core dumping on me)
    if(DEEBUG) print *, 'About to find grid list of file ', trim(DAOFile%Name)
    gridlist = ''
    inq_success = gdinqgrid(DAOFile%Name, gridlist, strbufsize)
    if (inq_success < 0) then
      call announce_error(lcf_where, "Could not inquire gridlist "// &
        & trim(DAOFile%Name))
    else if ( strbufsize > MAXLISTLENGTH ) then
       CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'list size too big in Read_DAO ' // trim(DAOFile%Name), MLSFile=DAOFile )
    else if ( strbufsize < MAXLISTLENGTH .and. strbufsize > 0 ) then
      gridlist = gridlist(1:strbufsize) // ' '
    end if
    if(DEEBUG) print *, 'grid list ', trim(gridlist)

    error = 0
    file_id = gdopen(DAOFile%Name, DFACC_RDONLY)

    if (file_id < 0) then
      call announce_error(lcf_where, "Could not open "// DAOFile%Name)
    end if

    ! Find grid name corresponding to the GRIDORDER'th one
    ngrids = NumStringElements(gridlist, COUNTEMPTY)

    if(ngrids <= 0) then
      call announce_error(lcf_where, "NumStringElements of gridlist <= 0")
    else if(ngrids /= inq_success) then
      call announce_error(lcf_where, "NumStringElements of gridlist /= inq_success")
    else if(ngrids < GRIDORDER) then
      call announce_error(lcf_where, "NumStringElements of gridlist < GRIDORDER")
    end if

    call GetStringElement(gridlist, gridname, GRIDORDER, COUNTEMPTY)

    gd_id = gdattach(file_id, gridname)
    if (gd_id < 0) then
      call announce_error(lcf_where, "Could not attach "//trim(gridname))
    end if

    ! Now find dimsize(), dimname(), etc.
    nentries = gdnentries(gd_id, HDFE_NENTDIM, strbufsize)

    if(nentries <= 0) then
      call announce_error(lcf_where, "nentries of gd_id <= 0")
    else if(nentries > NENTRIESMAX) then
      call announce_error(lcf_where, "nentries of gd_id > NENTRIESMAX")
    end if

    dimlist = ''
    ndims = gdinqdims(gd_id, dimlist, dims)

    if(ndims <= 0) then
      call announce_error(lcf_where, "ndims of gd_id <= 0")
    else if(ndims > NENTRIESMAX) then
      call announce_error(lcf_where, "ndims of gd_id > NENTRIESMAX")
    end if

    fieldlist = ''
    nfields = gdinqflds(gd_id, fieldlist, rank, numberTypes)

    if(nfields <= 0) then
      call announce_error(lcf_where, "nfields of gd_id <= 0")
    else if(nfields > NENTRIESMAX) then
      call announce_error(lcf_where, "nfields of gd_id > NENTRIESMAX")
    end if

    if(.not. CASESENSITIVE) then
      fieldlist = Capitalize(fieldlist)
    end if

    if(present(fieldName)) then
      actual_field_name=fieldName
    else
      actual_field_name=DEFAULTDAOFIELDNAME
    end if

    actual_dim_list = ' '
    if(present(GeoDimList)) then
      actual_dim_list=GeoDimList
    end if
    if ( actual_dim_list == ' ' ) then
      actual_dim_list = DEFAULTDAODIMLIST
    end if
    call List2Array (actual_dim_list, dimNames, countEmpty)

    ! Now find the rank of our field
    inq_success = gdfldinfo(gd_id, trim(actual_field_name), our_rank, dims, &
      & numbertype, dimlists(1))

    dimlist = trim(dimlists(1))
    if(DEEBUG) print *, 'our_rank ', our_rank
    if(DEEBUG) print *, 'dims ', dims(1:our_rank)
    if(DEEBUG) print *, 'dimlist ', dimlist

    nlon = dims(1)
    nlat = dims(2)
    nlev = dims(3)
    ntime = dims(4)

    the_g_data%quantityName = actual_field_name
    the_g_data%description = lit_dao
    the_g_data%verticalCoordinate = v_type

    the_g_data%noLons = nlon
    the_g_data%noLats = nlat
    the_g_data%noHeights = nlev
    ! the_g_data%noLsts = ntime
    the_g_data%noDates = ntime
    the_g_data%units = 'K'
    if(DEEBUG) print *, 'our quantity name ', the_g_data%quantityName
    if(DEEBUG) print *, 'our description ', the_g_data%description
    if(DEEBUG) print *, 'our units ', the_g_data%units
    if(DEEBUG) print *, 'our vertical coord ', the_g_data%verticalCoordinate
    if(DEEBUG) print *, 'v_type ', v_type

    ! Setup the grid
    call SetupNewGriddedData ( the_g_data, noHeights=nlev, noLats=nlat, &
      & noLons=nlon, noLsts=1, noSzas=1, noDates=ntime, &
      & missingValue=FILLVALUE, units='K', verticalCoordinate=v_type, &
      & heightsunits='hPa' )
      ! & noLons=nlon, noLsts=ntime, noSzas=1, noDates=1, missingValue=FILLVALUE )
    if(DEEBUG) print *, '(Again) our quantity name ', the_g_data%quantityName
    if(DEEBUG) print *, 'our description ', the_g_data%description
    if(DEEBUG) print *, 'our units ', the_g_data%units
    if(DEEBUG) print *, 'our vertical coord ', the_g_data%verticalCoordinate
    if(DEEBUG) print *, 'v_type ', v_type
    nullify ( all_the_fields )
    call allocate_test ( all_the_fields, dims(1), dims(2), dims(3), dims(4), &
      & "all_the_fields", moduleName )
    all_the_fields = the_g_data%missingValue
    start = 0                                                             
    stride = 1                                                               
    edge = dims(1:4)                                                        
    if(DEEBUG) print *, 'About to read ' // trim(actual_field_name)
    if(DEEBUG) print *, 'Start ', Start
    if(DEEBUG) print *, 'Stride ', Stride
    if(DEEBUG) print *, 'Edge ', Edge
    status = gdrdfld(gd_id, trim(actual_field_name), start, stride, edge, &  
      & all_the_fields)                                                          
    if(status /= 0) &
      & call announce_error(lcf_where, "failed to read field " &
      & //trim(actual_field_name))
    ! The actual dimlist is this                    XDim,YDim,Height,TIME
    ! Need to reshape it so that the order becomes: Height,YDim,XDim,TIME
    if ( DEEBUG) then
      print *, 'dao Before reshaping'
      call dump(all_the_fields(:,1,1,1), 'x-slice')
      call dump(all_the_fields(1,:,1,1), 'y-slice')
      call dump(all_the_fields(1,1,:,1), 'p-slice')
      call dump(all_the_fields(1,1,1,:), 't-slice')
    end if
    ! the_g_data%field(:,:,:,:,1,1) = reshape( all_the_fields, &
    the_g_data%field(:,:,:,1,1,:) = reshape( all_the_fields, &
      & shape=(/nLev, nlat, nlon, ntime/), order=(/3,2,1,4/) &
      & )
    if ( DEEBUG ) then
      print *, 'dao After reshaping'
      call dump(the_g_data%field(1,1,:,1,1,1), 'x-slice')
      call dump(the_g_data%field(1,:,1,1,1,1), 'y-slice')
      call dump(the_g_data%field(:,1,1,1,1,1), 'p-slice')
      call dump(the_g_data%field(1,1,1,:,1,1), 't-slice')
    end if
    call deallocate_test ( all_the_fields, "all_the_fields", moduleName )
    ! Now read the dims
    nullify(dim_field)
    ! call read_the_dim(gd_id, 'XDim', dims(1), dim_field)
    call read_the_dim(gd_id, trim(dimNames(1)), dims(1), dim_field)
    the_g_data%lons = dim_field
    ! call read_the_dim(gd_id, 'YDim', dims(2), dim_field)
    call read_the_dim(gd_id, trim(dimNames(2)), dims(2), dim_field)
    the_g_data%lats = dim_field
    ! call read_the_dim(gd_id, 'Height', dims(3), dim_field)
    call read_the_dim(gd_id, trim(dimNames(3)), dims(3), dim_field)
    the_g_data%Heights = dim_field
    ! call read_the_dim(gd_id, 'Time', dims(4), dim_field)
    call read_the_dim(gd_id, trim(dimNames(4)), dims(4), dim_field)
    ! the_g_data%lsts = dim_field
    the_g_data%dateStarts = dim_field
    the_g_data%dateEnds = dim_field
    call deallocate_test ( dim_field, "dim_field", moduleName )
    ! Have not yet figured out how to assign these
    ! Probably will have to read metadata
    the_g_data%Szas = the_g_data%missingValue
    ! the_g_data%DateStarts = the_g_data%missingValue
    ! the_g_data%DateEnds = the_g_data%missingValue
    the_g_data%lsts = the_g_data%missingValue
    ! Close grid
    status = gddetach(gd_id)
    if(status /= 0) &
      & call announce_error(lcf_where, "failed to detach from grid " &
      & //trim(gridname))
    status = gdclose(file_id)
    if(status /= 0) &
      & call announce_error(lcf_where, "failed to close file " //trim(DAOFile%Name))
  contains 
     subroutine read_the_dim(gd_id, field_name, field_size, values)
       ! Arguments
       integer, intent(in) :: gd_id
       character(len=*), intent(in) :: field_name
       integer, intent(in) :: field_size
       real(r8), dimension(:), pointer :: values
       ! Local variables
       integer :: status
       integer, dimension(1) :: start, stride, edge
       ! Executable
       call allocate_test ( values, field_size, "values", moduleName )
       start = 0                                                             
       stride = 1                                                             
       edge = field_size                                                       
       values = the_g_data%missingValue
       status = gdrdfld(gd_id, trim(field_name), start, stride, edge, &
         & values)                                                     
       if ( status /= 0 ) &
         & call announce_error(lcf_where, "failed to read dim field " &
         & // trim(field_name))
     end subroutine read_the_dim
  end subroutine Read_DAO

  ! ----------------------------------------------- Read_MERRA_2
  subroutine Read_MERRA_2( GEOS5File, lcf_where, v_type, &
    & the_g_data, GeoDimList, fieldName, date, sumDelp )

  use Dates_Module, only: GetStartingDate, UTC2TAI93s
  use Dump_1, only: Dump
  use HDF5, only: Hsize_T
  use MLSHDF5, only: GetAllHDF5DSNames, &
    & GetHDF5DSRank, GethDF5DSDims, LoadFromHDF5DS, &
    & ReadHDF5Attribute

    ! This routine reads a gmao merra2 file, named something like
    ! MERRA2_200.inst6_3d_ana_Np.19980724.nc4

    ! This file is formatted in the following way:
    ! each gridded quantity, not an hdfeos grid, by the way, is
    ! a rank 4 field, e.g. T or H, is given with dimensions
    ! 'Time, Height, YDim, XDim'
    ! We'll simply copy it into a single gridded data type

    ! Arguments
    type(MLSFile_T)                :: GEOS5file
    integer, intent(IN) :: lcf_where    ! node of the lcf that provoked me
    integer, intent(IN) :: v_type       ! vertical coordinate; an 'enumerated' type
    type( GriddedData_T ) :: the_g_data ! Result
    character (len=*), optional, intent(IN) :: GeoDimList ! Comma-delimited dim names
    character (len=*), optional, intent(IN) :: fieldName ! Name of gridded field
    character (len=*), optional, intent(IN) :: date ! date (offset)
    logical, optional, intent(IN)           :: sumDelp ! sum the DELP to make PL
    ! Local variables
    character (len=NAMELEN) :: actual_field_name
    integer(c_intptr_t) :: Addr         ! For tracing
    character(len=19) :: datestring ! will be in the form of yyyy-MM-ddTHH:MM:ss
    logical :: DEEBUG
    integer(kind=hsize_t) :: dim1(1), dim4(4)
    integer :: error, rank
    character(len=256) :: errormsg
    integer :: i
    integer :: itime, ilat, ilon, ilev
    integer, dimension(4) :: idim4
    character (len=MAXSDNAMESBUFSIZE) :: mySdList
    logical                        :: mySum
    integer :: S                     ! Size in bytes of an object to deallocate
    real(r8), dimension(:), pointer :: temp1d => null()
    real(r8), dimension(:,:,:,:), pointer :: temp4d => null()
    character(len=16) :: the_units
    real(r4), parameter :: FILLVALUE = 1.e15 !this value might be wrong
    integer :: mydate, mytime, year, month, day, hour, minute, second
    logical :: verbose

    ! Executable
    if(present(fieldName)) then
      actual_field_name=fieldName
    else
      actual_field_name=DEFAULTGEOS5FIELDNAME
    end if
    mySum = .false.
    if ( present(sumDelp) ) mySum = sumDelp
    nullify( temp1d, temp4d )
    DEEBUG = ( index(lowercase(actual_field_name), 'inq') > 0 )
    verbose = ( switchDetail(switches, 'merra_2') > -1 ) .or. DEEBUG
    call GetAllHDF5DSNames ( GEOS5File%Name, '/', mysdList )
    if ( verbose ) call dump( mysdList, 'DS names' )

    ! Fill in value
    ! Note that we could read the 'long_name' attribute here
    the_g_data%quantityName = actual_field_name
    the_g_data%description = lit_geos5
    the_g_data%verticalCoordinate = v_type
    the_g_data%nodates = 1
    the_g_data%empty = .true.
    ! Note that we could read the 'missing_value' attribute here
    the_g_data%missingvalue = FILLVALUE ! this value might be wrong
    call allocate_test ( the_g_data%lsts, 1, "the_g_data%lsts", moduleName )
    the_g_data%lsts = the_g_data%missingvalue
    the_g_data%nolsts = 1
    the_g_data%noszas = 1
    the_g_data%lsts = the_g_data%missingValue ! Know how to read this yet?
    call allocate_test ( the_g_data%szas, 1, "the_g_data%szas", moduleName )
    the_g_data%szas = the_g_data%missingvalue

    call mls_openfile(geos5File, error)
    if (error .gt. 0) call announce_error(lcf_where, "Could not open "// GEOS5File%Name)

    ! Read the begin_date attribute
    if (.not. ReadHDF5Attribute(geos5file%fileid%f_id, 'time', &
        & 'begin_date', mydate, error=errormsg)) then
        call announce_error (lcf_where, errormsg // ' in file ' &
        & // geos5file%name)
    end if

    if (mydate < 0) then
        call announce_error (lcf_where, "Invalid 'begin_date' value in " // geos5file%name)
    end if

    year = mydate / 10000
    month = mod(mydate, 10000) / 100
    day = mod(mydate, 100)
    if (year > 9999 .or. month > 12 .or. month < 1 .or. day > 31 .or. day < 1) then
        call announce_error (lcf_where, "Invalid 'begin_date' value in " // geos5file%name)
    end if
    ! print *, 'year, month, day: ', year, month, day
    
    ! Get lons
    call GetHDF5DSRank (geos5file%fileid%f_id, 'lon', rank)
    if (rank /= 1) call announce_error(lcf_where, "lon must be a 1-dim array: "// GEOS5File%Name)

    call GetHDF5DSDims(geos5file%fileid%f_id, 'lon', dim1)
    the_g_data%noLons = dim1(1)

    call allocate_test ( temp1d, the_g_data%noLons, "temp1d", moduleName )

    call LoadFromHDF5DS (geos5file%fileid%f_id, 'lon', temp1d)
    call allocate_test ( the_g_data%lons, the_g_data%nolons, "the_g_data%nolons", &
      & moduleName )
    the_g_data%lons = temp1d
    s = size(temp1d) * storage_size(temp1d) / 8
    call deallocate_test ( temp1d, "temp1d", moduleName )

    ! Defer Filling dateStart and dateEnd

    ! Get lats
    call GetHDF5DSRank (geos5file%fileid%f_id, 'lat', rank)
    if (rank /= 1) &
      & call announce_error(lcf_where, "lat must be a 1-dim array: " // geos5file%name)

    call GetHDF5DSDims(geos5file%fileid%f_id, 'lat', dim1)
    the_g_data%noLats = dim1(1)

    call allocate_test ( temp1d, the_g_data%nolats, "temp1d", moduleName )

    call LoadFromHDF5DS (geos5file%fileid%f_id, 'lat', temp1d)
    call allocate_test ( the_g_data%lats, the_g_data%nolats, "the_g_data%lats", &
      & moduleName )
    the_g_data%lats = temp1d
    call deallocate_test ( temp1d, "temp1d", moduleName )

    ! Get heights
    call GetHDF5DSRank (geos5file%fileid%f_id, 'lev', rank)
    if (rank /= 1) &
      & call announce_error(lcf_where, "lev must be a 1-dim array: " // geos5file%name)
    
    call GetHDF5DSDims (geos5file%fileid%f_id, 'lev', dim1)
    the_g_data%noheights = dim1(1)

    call allocate_test ( temp1d, the_g_data%noheights, "temp1d", moduleName )
    
    call LoadFromHDF5DS (geos5file%fileid%f_id, 'lev', temp1d)
    call allocate_test ( the_g_data%heights, the_g_data%noheights, &
      & "the_g_data%heights", moduleName )
    
    ! We cannot load directly into the the_g_data's array because
    ! this is an array of 32-bit float, while the data from file is
    ! 64-bit float
    the_g_data%heights = temp1d
    call deallocate_test ( temp1d, "temp1d", moduleName )

    the_g_data%heightsunits = 'hPa' ! the units is according the file_specification

    ! The following is according to Merra_2 file specification
    ! First: try to read the 'units' attribute
    if (.not. ReadHDF5Attribute(geos5file%fileid%f_id, actual_field_name, &
        & 'units', the_g_data%units, error=errormsg) ) then
      call announce_error ( lcf_where, errormsg // ' in file ' &
      & // geos5file%name, non_fatal=.true. )
      select case ( lowercase(actual_field_name) )
      case ( 'h' )
          the_units = 'm'
      case ( 'pl', 'ps', 'delp' )
          the_units = 'Pa'
      case ( 't' )
          the_units = 'K'
      case default
        call announce_error( lcf_where, &
          & "Unexpected field name: " // actual_field_name, non_fatal=.true. )
      end select
      the_g_data%units = the_units
    else
      call outputNamedValue( 'Name of units attribute for ' &
        & // trim(actual_field_name), &
        & the_g_data%units )
    endif

    ! read the field, field should be either PL or T
    call GetHDF5DSRank (geos5file%fileid%f_id, capitalize(actual_field_name), rank)
    if (rank /= 4) then
        call announce_error (lcf_where, &
        trim(capitalize(actual_field_name)) // " must be a 4-dim array: (merra2)" // geos5file%name)
    end if

    call GetHDF5DSDims (geos5file%fileid%f_id, capitalize(actual_field_name), dim4)
    idim4 = dim4

    the_g_data%nodates = idim4(4)
    if ( verbose ) call outputnamedValue ( 'nodates', the_g_data%nodates )
    ! call Dump( geos5file )

    ! Now assign values to DateStarts and Ends
    call allocate_test ( the_g_data%datestarts, the_g_data%nodates, &
      & "the_g_data%datestarts", moduleName )
    call allocate_test ( the_g_data%dateends, the_g_data%nodates, &
      & "the_g_data%dateends", moduleName )
    ! Let's just assume we can read a data field called 'time'
    ! If not, we're dead anyway
    call allocate_test ( temp1d, the_g_data%nodates, "temp1d", moduleName )

    call LoadFromHDF5DS ( geos5file%fileid%f_id, 'time', temp1d )
    if ( verbose ) call Dump( temp1d, 'time field' )
    ! This time field turns out to be minutes since the start of the day 
    ! (Of course! How else would anyone code a time field?)
    do i=1, idim4(4)
      mytime = 60*temp1d(i) ! to convert to second
      second = mod( mytime, 60 )
      minute = mod( mytime, 3600 ) - second
      hour = mytime / 3600
      write ( datestring, &
        & '(I4.4, A1,   I2.2,  A1, I2.2,  A1, I2.2, A1,  I2.2,    A1, I2.2)') &
        &   year, '-',  month, '-', day, 'T', hour, ':', minute, ':', second
      the_g_data%datestarts(i) = utc2tai93s(datestring)
    enddo
    if ( verbose ) then
      call outputNamedValue( 'datestring', trim(datestring) )
      call GetStartingDate( datestring )
      call outputNamedValue( 'Starting date', trim(datestring) )
      call Dump( the_g_data%datestarts, 'datestarts field' )
    endif
    ! Next we will read an attribute called "begin_date"
    the_g_data%dateends = the_g_data%datestarts
    call deallocate_test ( temp1d, "temp1d", moduleName )

    call allocate_test ( temp4d, idim4(1), idim4(2), idim4(3), idim4(4), &
        & 'temp4d', ModuleName // 'Read_geos57' )

    call LoadFromHDF5DS ( geos5file%fileid%f_id, &
      & capitalize(actual_field_name), temp4d )

    if ( mySum ) then
      ! If we are so asked
      ! Sum up all the delps, starting from the top of the atmosphere
      nullify(temp1d)
      call allocate_test( temp1d, idim4(2)+1, 'temp1d', &
        & ModuleName // 'Read__Merra_2' )
      do itime=1, idim4(4) ! ntime
        do ilat=1, idim4(2) ! nlat
          do ilon=1, idim4(1) ! nlon
            temp1d = 0.0
            ! temp1d(nlev+1) = 1.0 ! We must assume the topmost level is 1 hPa
            ! do ilev=nlev, 1, -1
            !  temp1d(ilev) = temp1d(ilev+1) + all_the_fields(ilon, ilat, ilev, itime)
            temp1d(1) = 1.0 ! We must assume the topmost level is 1 hPa
            do ilev=1, idim4(3)
              temp1d(ilev+1) = temp1d(ilev) + temp4d(ilon, ilat, ilev, itime)
            end do
            do ilev=1, idim4(3)
              temp4d(ilon, ilat, ilev, itime) = &
                & 0.5*( temp1d(ilev) + temp1d(ilev+1) )
            end do
          end do
        end do
      end do
      call deallocate_test( temp1d, 'temp1d', &
        & ModuleName // 'Read__Merra_2' )
    end if
    allocate ( the_g_data%field &
      & (idim4(3), idim4(2), idim4(1), 1, 1, idim4(4)), &
      & STAT=error )
    addr = 0
    if ( error == 0 ) then
      if ( size(the_g_data%field) > 0 ) addr = transfer( c_loc( &
        & the_g_data%field(1,1,1,1,1,1)), addr )
    end if
    call test_allocate ( error, moduleName, 'the_g_data%field', &
      & (/1,1,1,1,1,1/), &
      & (/ idim4(3), idim4(2), idim4(1), 1, 1, idim4(4) /), &
      & bytes(the_g_data%field), &
      & address=addr )
    the_g_data%field(:,:,:,1,1,:) = &
      & reshape( temp4d(:,:,:,:), order=(/3,2,1,4/), &
      & shape=(/ the_g_data%noheights, the_g_data%nolats, &
      & the_g_data%nolons, the_g_data%noDates /) )
    call deallocate_test ( temp4d, 'temp4d', ModuleName // 'Read_geos57' )
    
    ! Read file successful
    the_g_data%empty = .false.
    call mls_closefile(geos5File, error)
    if (error .gt. 0) call announce_error(lcf_where, "Could not close "// GEOS5File%Name)
    if ( verbose ) call dump( the_g_data, details = 1 )
  end subroutine Read_Merra_2

  ! ---------------------------------------------  Read_NCEP_GDAs  -----
  subroutine Read_NCEP_GDAs ( NCEPFile, lcf_where, v_type, &
    & the_g_data, GeoDimList, gridName, missingValue )

    ! This routine reads a ncep global assimilation model product
    ! (GDAS) file, named something like
    ! 6207ea2e-1dd2-11b2-b61e-ae069991db9a.hdfeos 
    ! returning a filled data structure appropriate for it

    ! FileName and the_g_data are required args
    ! GeoDimList, if present, should be the Dimensions' short names
    ! as a comma-delimited character string in the order:
    ! longitude, latitude, vertical level, time

    ! fieldName, if present, should be the rank 2 object
    ! like temperature
    
    ! This file is formatted in the following way:
    ! At each gridded quantity, e.g. TMP_3,
    ! a series of 2d slices at a fixed pressure surface are
    ! given as fields, with names like "ISOBARIC LEVEL AT 975 (hPa)"
    ! We will read them, stacking them up in a single gridded data type
    
    ! ncep is peculiar in 3 other ways:
    ! (1) It has 0 entries
    ! (2) It has 0 dimensions
    ! (3) (Related to 2) The lats and lons have to be inferred, not read

    ! Arguments
    type(MLSFile_T)                :: NCEPFile
    integer, intent(IN) :: lcf_where    ! node of the lcf that provoked me
    integer, intent(IN) :: v_type       ! vertical coordinate
    type( GriddedData_T ) :: the_g_data ! Result
    character (len=*), optional, intent(IN) :: GeoDimList ! E.g., 'X.Y,..'
    character (len=*), optional, intent(IN) :: gridName ! Name of grid
    real(rgr), optional, intent(IN) :: missingValue

    ! Local Variables
    integer :: file_id, gd_id
    integer :: inq_success
    integer :: nentries, ndims, nfields
    integer :: strbufsize
    character(len=NameLen) :: my_gridname
    integer, parameter :: MAXLISTLENGTH=10*Linelen ! Max length list grid names
    character (len=MAXLISTLENGTH), dimension(1) :: dimlists
    character (len=MAXLISTLENGTH) :: fieldlist
    character (len=MAXLISTLENGTH) :: gridlist
    character (len=MAXLISTLENGTH) :: dimlist
    integer, parameter :: MAXNAMELENGTH=NameLen         ! Max length of grid name
    character (len=MAXNAMELENGTH) :: actual_field_name
    integer, dimension(NENTRIESMAX) :: dims, rank, numberTypes
    integer                        :: our_rank, numberType

    integer :: start(2), stride(2), edge(2)
    integer :: status
    integer :: field
    !                                  These start out initialized to one
    integer                        :: nlon=1, nlat=1, nlev=1
    integer, external :: GDRDFLD
    logical, parameter :: USEPROJECTFORDIMS = .false. ! This doesn't work, yet
    real(r4), dimension(:), pointer     :: pressures
    real(r4), dimension(:,:), pointer   :: field_data
    real(r4), dimension(:,:,:), pointer :: all_the_fields, now_the_fields
    real(r8) :: pressure
    integer :: xdimsize, ydimsize
    real(r8), dimension(2) :: upLeft, lowRight
    real(r8), dimension(15) :: projParams
    integer :: projcode, zonecode, spherecode, directioncode
    integer :: i
    integer, external :: gdgridinfo, GDprojinfo, PGS_GCT_Init
    integer, parameter :: Lon_offset = -180  ! Longitude start in degrees
    integer, parameter :: Lat_offset = -90   ! Latitude start in degrees
    logical, parameter :: CHECKFORDIMSANYWAY = .false.
    logical, parameter :: DEEBUG = .false.
    ! Executable
    projParams = 0.d0
    my_gridname = DEFAULTNCEPGRIDNAME
    if ( present(gridname) ) my_gridname = gridname
    ! Don't find list of grid names on this file --because it dumped core on me
    gridlist = ''
    inq_success = 0 ! inq_success = gdinqgrid(FileName, gridlist, strbufsize)
    error = 0
    file_id = gdopen(NCEPFile%name, DFACC_RDONLY)
    if(DEEBUG) print *, 'fileID: ', file_id

    if (file_id < 0) then
      call announce_error(lcf_where, "Could not open "// NCEPFile%name)
    end if

    gd_id = gdattach(file_id, my_gridname)
    if(DEEBUG) print *, 'gridID: ', gd_id
    if (gd_id < 0) then
      call announce_error(lcf_where, "Could not attach "//trim(my_gridname))
    end if

    if ( USEPROJECTFORDIMS ) then
      ! Projection info
      status = GDprojinfo(gd_id, projcode, zonecode, spherecode, projParams)
      if (status /= 0) then
        call announce_error(lcf_where, "Could not get info on projection")
      end if
      if(DEEBUG) print *, 'projcode: ', projcode
      if(DEEBUG) print *, 'zonecode: ', zonecode
      if(DEEBUG) print *, 'spherecode: ', spherecode
      if(DEEBUG) print *, 'projParams: ', projParams
      directioncode = PGSd_GCT_INVERSE
      status = PGS_GCT_Init(projcode, projParams, directioncode)
      if (status /= 0) then
        call announce_error(lcf_where, "Could not initialize gct")
      end if
    end if
    ! Now find dimsize(), dimname(), etc.
    ! ncep file don't have entries, nor dims
    if ( CHECKFORDIMSANYWAY ) then
      nentries = gdnentries(gd_id, HDFE_NENTDIM, strbufsize)

      if(nentries <= 0) then
        call announce_error(lcf_where, "nentries of gd_id <= 0", &
          & 'nentries:', nentries)
      else if(nentries > NENTRIESMAX) then
        call announce_error(lcf_where, "nentries of gd_id > NENTRIESMAX", &
          & 'nentries:', nentries)
      end if

      dimlist = ''
      ndims = gdinqdims(gd_id, dimlist, dims)

      if(ndims <= 0) then
        call announce_error(lcf_where, "ndims of gd_id <= 0", &
          & 'ndims:', ndims)
      else if(ndims > NENTRIESMAX) then
        call announce_error(lcf_where, "ndims of gd_id > NENTRIESMAX", &
          & 'ndims:', ndims)
      end if
    end if

    status = gdgridinfo(gd_id, xdimsize, ydimsize, &
      & upLeft, lowRight )

    if(status /= 0) then
      call announce_error(lcf_where, "failed to obtain gridinfo", &
        & 'status:', status)
    end if
    if(DEEBUG) print *, 'xdimsize: ', xdimsize
    if(DEEBUG) print *, 'ydimsize: ', ydimsize
    if(DEEBUG) print *, 'upLeft  : ', upLeft
    if(DEEBUG) print *, 'lowRight: ', lowRight

    fieldlist = ''
    nfields = gdinqflds(gd_id, fieldlist, rank, numberTypes)
    if(DEEBUG) print *, 'nfields: ', nfields

    if(nfields <= 0) then
      call announce_error(lcf_where, "nfields of gd_id <= 0", &
        & 'nfields:', nfields)
      return
    else if(nfields > NENTRIESMAX) then
      call announce_error(lcf_where, "nfields of gd_id > NENTRIESMAX", &
        & 'nfields:', nfields)
      return
    end if

    the_g_data%noHeights = 0
    nullify ( pressures )
    call allocate_test ( pressures, nfields, "pressures", moduleName )
    ! Loop over data fields
    nullify ( all_the_fields, now_the_fields )
    do field=1, nfields
      call getStringElement(fieldlist, actual_field_name, field, .false.)
      if(DEEBUG) print *, 'actual_field_name: ', trim(actual_field_name)
      if ( actual_field_name == ' ' ) cycle
      call ncepFieldNameTohPa(trim(actual_field_name), pressure)
      if(DEEBUG) print *, 'inferred pressure: ', pressure
      if ( pressure < 0._r8 ) cycle
      ! Now find the rank of our field
      inq_success = gdfldinfo(gd_id, trim(actual_field_name), our_rank, dims, &
        & numbertype, dimlists(1))

      if ( our_rank /= 2 ) then
        call announce_error(lcf_where, "rank /= 2")
        cycle
      end if
      call allocate_test ( field_data, dims(1), dims(2), "field_data", moduleName )
      
      field_data = the_g_data%missingValue
      dimlist = trim(dimlists(1))

      nlon = dims(1)
      nlat = dims(2)

      the_g_data%quantityName = my_gridname
      the_g_data%description = lit_ncep
      the_g_data%verticalCoordinate = v_type

      the_g_data%noLons = nlon
      the_g_data%noLats = nlat
      the_g_data%noHeights = the_g_data%noHeights + 1
      pressures(the_g_data%noHeights) = pressure
      start(1) = 0
      start(2) = 0
      stride = 1
      edge(1) = dims(1)
      edge(2) = dims(2)
      status = gdrdfld(gd_id, trim(actual_field_name), start, stride, edge, &
        & field_data)
      if(status /= 0) &
        & call announce_error(lcf_where, "failed to read field " &
        & //trim(actual_field_name))
      ! call dump(field_data(1,:), trim(actual_field_name))
      ! Assemble the_fields into 3-d array all_the_fields
      if (the_g_data%noHeights == 1) then
        call allocate_test ( all_the_fields, 1, dims(1), dims(2), "all_the_fields", &
          & moduleName )
        all_the_fields(1, :, :) = field_data
      else
        call allocate_test ( now_the_fields, the_g_data%noHeights, dims(1), dims(2), &
          & "now_the_fields", moduleName )
        now_the_fields(1:the_g_data%noHeights - 1, :, :) = all_the_fields
        now_the_fields(the_g_data%noHeights, :, :) = field_data
        call deallocate_test ( all_the_fields, "all_the_fields", moduleName )
        all_the_fields => now_the_fields
      end if
      call deallocate_test ( field_data, "field_data", moduleName )
    end do
    the_g_data%noLsts = 0
    the_g_data%units = 'K'
    if(DEEBUG) print *, 'our quantity name ', the_g_data%quantityName
    if(DEEBUG) print *, 'our description ', the_g_data%description
    if(DEEBUG) print *, 'our units ', the_g_data%units

    ! Setup the grid
    nLev = the_g_data%noHeights
    if(DEEBUG) print *, 'NoHeights: ', nLev
    if(DEEBUG) print *, 'NoLons   : ', nLon
    if(DEEBUG) print *, 'NoLats   : ', nLat
    call SetupNewGriddedData ( the_g_data, noHeights=nlev, noLats=nlat, &
      & noLons=nlon, noLsts=1, noSzas=1, noDates=1 )
    if(DEEBUG) print *, '(Again) our quantity name ', the_g_data%quantityName
    if(DEEBUG) print *, 'our description ', the_g_data%description
    if(DEEBUG) print *, 'our units ', the_g_data%units
    ! The dimlist as stacked up is this             Height,XDim,YDim
    ! Need to reshape it so that the order becomes: Height,YDim,XDim
    if(DEEBUG) then
      print *, 'Before reshaping'
      call dump(all_the_fields(1,:,1), 'x-slice')
      call dump(all_the_fields(1,1,:), 'y-slice')
      call dump(all_the_fields(:,1,1), 'p-slice')
    end if
    the_g_data%field(:,:,:,1,1,1) = reshape( all_the_fields, &
      & shape=(/nLev, nlat, nlon/), order=(/1, 3, 2/) &
      & )
    if(DEEBUG) then
      print *, 'After reshaping'
      call dump(the_g_data%field(1,1,:,1,1,1), 'x-slice')
      call dump(the_g_data%field(1,:,1,1,1,1), 'y-slice')
      call dump(the_g_data%field(:,1,1,1,1,1), 'p-slice')
    end if
    the_g_data%heights = pressures(1:the_g_data%noHeights)
    call deallocate_test ( all_the_fields, "all_the_fields", moduleName )
    call deallocate_test ( pressures, "pressures", moduleName )

    ! Insert code to transform upLeft and LowRight
    ! info into lats and lons
    if ( USEPROJECTFORDIMS ) then
      call announce_error(lcf_where, "Read_NCEP_GDAs unable to use projection")
    else
    ! Now read the dims
    !  But .. they aren't there!
      do i=1, xdimsize
        the_g_data%lons(i) = i - 1 + Lon_offset
      end do
      do i=1, ydimsize
        the_g_data%lats(i) = i - 1 + Lat_offset
      end do
    end if
    ! Have not yet figured out how to assign these
    ! Probably will have to read metadata
    the_g_data%lsts = the_g_data%missingValue
    the_g_data%Szas = the_g_data%missingValue
    the_g_data%DateStarts = the_g_data%missingValue
    the_g_data%DateEnds = the_g_data%missingValue

    ! Close grid
    status = gddetach(gd_id)
    if(status /= 0) &
      & call announce_error(lcf_where, "failed to detach from grid " // &
      & trim(my_gridname))
    status = gdclose(file_id)
    if(status /= 0) &
      & call announce_error(lcf_where, "failed to close file " //trim(NCEPFile%name))
  contains 
     subroutine read_the_dim ( gd_id, field_name, field_size, values )
       ! Arguments
       integer, intent(in) :: gd_id
       character(len=*), intent(in) :: field_name
       integer, intent(in) :: field_size
       real(r4), dimension(:), pointer :: values
       ! Local variables
       integer :: status
       integer, dimension(1) :: start, stride, edge
       ! Executable
       call allocate_test ( values, field_size, "values", moduleName )
       start = 0                                                             
       stride = 1                                                             
       edge = field_size
       values = the_g_data%missingValue
       status = gdrdfld(gd_id, trim(field_name), start, stride, edge, &
         & values)                                                     
       if ( status /= 0 ) &
         & call announce_error(lcf_where, "failed to read dim field " &
         & // trim(field_name))
     end subroutine read_the_dim
  end subroutine Read_NCEP_GDAs

  ! --------------------------------------------  Read_NCEP_Strat  -----
  subroutine Read_NCEP_Strat(NCEPFile, lcf_where, v_type, &
    & the_g_data, GeoDimList, fieldName)
    use HDFeos5, only: He5_HDFe_Nentdim, He5f_Acc_Rdonly, &
      & He5_Gdopen, He5_Gdattach, He5_Gddetach, He5_Gdclose, &
      & He5_Gdnentries, He5_Gdinqgrid, He5_Gdinqdims, He5_Gdinqflds, &
      & He5_Gdfldinfo, He5_Gdgridinfo
    use MLSHDFeos, only: Hsizes
    use HDF5, only: Size_T

    ! This routine reads a ncep stratospheric combined product file,
    ! named something like nmct_030126.he5
    ! returning a filled data structure appropriate for it

    ! FileName and the_g_data are required args
    ! GeoDimList, if present, should be the Dimensions' short names
    ! as a comma-delimited character string in the order:
    ! longitude, latitude, vertical level, time

    ! The file itself determines which data field, e.g., Temperature,
    ! is being defined
    
    ! fieldName, if present, should be a character string matching one of the
    ! fields in the grid: {'Temperature', 'Moisture', 'Height'}
    
    ! gridName, if present, should be a character string matching one of the
    ! grids in the file: {'NorthernHemisphere', 'SouthernHemisphere'}
    
    ! This file is formatted in the following way:
    ! At each gridded quantity, e.g. /NorthernHemisphere/Temperature,
    ! a Polar stereographic projection is defined

    ! Arguments
    type(MLSFile_T)                :: NCEPFile
    integer, intent(IN) :: lcf_where    ! node of the lcf that provoked me
    integer, intent(IN) :: v_type       ! vertical coordinate
    type( GriddedData_T ) :: the_g_data ! Result
    character (len=*), optional, intent(IN) :: GeoDimList ! E.g., 'X.Y,..'
    character (len=*), optional, intent(IN) :: fieldName  ! Name of grid quant.

    ! Local Variables
    integer :: file_id
    integer :: gd_id(2)
    integer :: inq_success
    integer :: nentries, ngrids, ndims, nfields
    integer(kind=size_t) :: strbufsize
    logical,  parameter       :: CASESENSITIVE = .false.
    integer, parameter :: GRIDORDER=1   ! What order grid written to file
    character (len=MAXLISTLENGTH) :: gridlist
    character (len=MAXLISTLENGTH) :: dimlist, actual_dim_list
    character (len=MAXLISTLENGTH), dimension(1) :: dimlists
    character (len=MAXLISTLENGTH), dimension(1) :: maxdimlists
    character (len=16), dimension(NENTRIESMAX) :: dimNames
    character (len=MAXLISTLENGTH), dimension(2) :: names
    character (len=MAXLISTLENGTH) :: fieldlist
    integer, parameter :: MAXNAMELENGTH=NameLen         ! Max length of grid name
    character (len=MAXNAMELENGTH) :: actual_field_name
    integer, dimension(NENTRIESMAX) :: rank, numberTypes
    integer(kind=size_t), dimension(NENTRIESMAX) :: dims
    integer                        :: our_rank, numberType
    integer, dimension(NENTRIESMAX,NENTRIESMAX) :: dims_temp

    integer :: start(3), stride(3), edge(3)
    integer :: status, j
    !                                  These start out initialized to one
    integer                        :: nlon=0, nlat=0, nlev=0, ntime=1
    real(r4), parameter :: FILLVALUE = 1.e15
    !real(r4), dimension(:,:,:,:), pointer :: all_the_fields
    real(r4), dimension(:,:,:), pointer :: all_the_fields
    real(r4), dimension(:,:,:), pointer :: t_all_fields
    logical, parameter :: DEEBUG = .false.
    real(r8):: upleftpt(2), lowrightpt(2)
    integer :: i
    integer(kind=size_t) :: xdimsize=1, ydimsize=2
    integer, parameter :: Lon_offset = -180  ! Longitude start in degrees
    integer, parameter :: Lat_offset = -90   ! Latitude start in degrees
    integer :: fdims3, sdims3

    ! Some day we'll get all these from MLSHDFEOS
    integer, external :: he5_gdrdfld
    ! Executable code
    ! Find list of grid names on this file (This has been core dumping on me)
    if(DEEBUG) print *, 'About to find grid list of file ', trim(NCEPFile%name)
    gridlist = ''
    inq_success = he5_gdinqgrid(NCEPFile%name, gridlist, strbufsize)
    if (inq_success < 0) then
      call announce_error(lcf_where, &
        & "Could not inquire gridlist "// trim(NCEPFile%name))
    else if ( strbufsize > MAXLISTLENGTH ) then
       CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'list size too big in Read_ncep ' // trim(NCEPFile%name), MLSFile=NCEPFile )
    else if ( strbufsize < MAXLISTLENGTH .and. strbufsize > 0 ) then
      gridlist = gridlist(1:strbufsize) // ' '
    end if

    if(DEEBUG) print *, 'grid list ', trim(gridlist)
    error = 0

    file_id = he5_gdopen(NCEPFile%name, HE5F_ACC_RDONLY)
    if (file_id < 0) then
      call announce_error(lcf_where, "Could not open "// NCEPFile%name)
    end if

    ! Find grid name corresponding to the GRIDORDER'th one
    ngrids = NumStringElements(gridlist, COUNTEMPTY)
    if(ngrids <= 0) then
      call announce_error(lcf_where, "NumStringElements of gridlist <= 0")
    else if(ngrids /= inq_success) then
      call announce_error(lcf_where, &
        & "NumStringElements of gridlist /= inq_success")
    else if(ngrids < GRIDORDER) then
      call announce_error(lcf_where, &
        & "NumStringElements of gridlist < GRIDORDER")
    end if

    !if ( present(gridName) ) then
    !    mygridname = gridName
    !else
      !call GetStringElement(gridlist, mygridname, GRIDORDER, COUNTEMPTY)
    !end if
    do i=1, ngrids
      call GetStringElement(gridlist, names(i), i, COUNTEMPTY)
      if (DEEBUG) print *,'name = ', trim(names(i))
        gd_id(i) = he5_gdattach(file_id, trim(names(i)))
      if (gd_id(i) < 0) then
          !call announce_error(lcf_where,"Could not attach "//trim(mygridname))
        call MLSMessage (MLSMSG_Warning, ModuleName, & 
          & "Could not attach "//trim(names(i)))
           exit 
      end if
        !Now find dimsize(), dimname(), etc.
        nentries = he5_gdnentries(gd_id(i), HE5_HDFE_NENTDIM, strbufsize)

        if(nentries <= 0) then
          call announce_error(lcf_where, "nentries of gd_id <= 0")
        else if(nentries > NENTRIESMAX) then
           call announce_error(lcf_where, "nentries of gd_id > NENTRIESMAX")
        end if

        dimlist = ''
        ndims = he5_gdinqdims(gd_id(i), dimlist, dims)
        if(ndims <= 0) then
           call announce_error(lcf_where, "ndims of gd_id <= 0")
        else if(ndims > NENTRIESMAX) then
           call announce_error(lcf_where, "ndims of gd_id > NENTRIESMAX")
        end if

        fieldlist = ''
        nfields = he5_gdinqflds(gd_id(i), fieldlist, rank, numberTypes)
        if(nfields <= 0) then
           call announce_error(lcf_where, "nfields of gd_id <= 0")
        else if(nfields > NENTRIESMAX) then
           call announce_error(lcf_where, "nfields of gd_id > NENTRIESMAX")
        end if

        if(.not. CASESENSITIVE) then
           fieldlist = Capitalize(fieldlist)
        end if

        if(present(fieldName)) then
           actual_field_name=fieldName
        else
           actual_field_name=DEFAULTNCEPSTRATFIELDNAME
        end if

        actual_dim_list = ' '
        if(present(GeoDimList)) then
           actual_dim_list=GeoDimList
        end if
        if ( actual_dim_list == ' ' ) then
           actual_dim_list = DEFAULTDAODIMLIST
        end if

        call List2Array (actual_dim_list, dimNames, countEmpty)

        ! Now find the rank of our field
        inq_success = he5_gdfldinfo(gd_id(i), trim(actual_field_name), & 
          & our_rank, dims, numbertype, dimlists(1), maxdimlists(1))

        dimlist = trim(dimlists(1))
        dims_temp(1:our_rank,i) = dims(1:our_rank)

        if(DEEBUG) print *, 'our_rank ', our_rank
        if(DEEBUG) print *, 'dims ', dims(1:our_rank)
        if(DEEBUG) print *, 'dimlist ', dimlist

        nlon = dims(1) 
        nlat = dims(2) 
        nlev = nlev + dims(3)
        if(DEEBUG) print *, 'nlon, nlat, nlev: ',nlon, nlat, nlev
    end do

    the_g_data%quantityName = actual_field_name
    the_g_data%description = lit_strat
    the_g_data%verticalCoordinate = v_type

    the_g_data%noLons = nlon
    the_g_data%noLats = nlat
    the_g_data%noHeights = nlev
    ! the_g_data%noLsts = ntime
    the_g_data%noDates = ntime
    the_g_data%units = 'K'
    if(DEEBUG) print *, 'our quantity name ', the_g_data%quantityName
    if(DEEBUG) print *, 'our description ', the_g_data%description
    if(DEEBUG) print *, 'our units ', the_g_data%units

    ! Setup the grid
    call SetupNewGriddedData ( the_g_data, noHeights=nlev, noLats=nlat, &
      & noLons=nlon, noLsts=1, noSzas=1, noDates=ntime, missingValue=FILLVALUE )
      ! & noLons=nlon, noLsts=ntime, noSzas=1, noDates=1, missingValue=FILLVALUE )
    if(DEEBUG) print *, '(Again) our quantity name ', the_g_data%quantityName
    if(DEEBUG) print *, 'our description ', the_g_data%description
    if(DEEBUG) print *, 'our units ', the_g_data%units
    nullify ( t_all_fields )
    call allocate_test ( t_all_fields, nlon, nlat, nlev, "t_all_fields", &
      & moduleName )
    t_all_fields = the_g_data%missingValue

    nullify ( all_the_fields )
    do i=1, ngrids
       if (gd_id(i) < 0) then
         call MLSMessage (MLSMSG_Warning, ModuleName, & 
           & "Could not attach "//trim(names(i)))
         exit 
       end if

       dims(1) = dims_temp(1,i)
       dims(2) = dims_temp(2,i)
       dims(3) = dims_temp(3,i)
       if (DEEBUG) print *, 'dims1, dims2, dims3: ', dims(1), dims(2), dims(3) 

       call allocate_test ( all_the_fields, int(dims(1)), int(dims(2)), int(dims(3)), &
         & "all_the_fields", moduleName )
       all_the_fields = the_g_data%missingValue

       status = he5_gdgridinfo(gd_id(i), xdimsize, ydimsize, upleftpt, &
          & lowrightpt)
       start = 0
       stride = 1
       !edge = dims(1:3)
       edge(1) = xdimsize 
       edge(2) = ydimsize 
       edge(3) = dims(3)
       if(DEEBUG) print *, 'About to read ' // trim(actual_field_name)
       if(DEEBUG) print *, 'Start ', Start
       if(DEEBUG) print *, 'Stride ', Stride
       if(DEEBUG) print *, 'Edge ', Edge
       if(DEEBUG) print *, 'xdimsize ', xdimsize
       if(DEEBUG) print *, 'ydimsize ', ydimsize

       status = he5_gdrdfld(gd_id(i), trim(actual_field_name), &
         & hsizes(start), hsizes(stride), hsizes(edge), all_the_fields)
       if(status /= 0) &
          & call announce_error(lcf_where, "failed to read field " &
          & //trim(actual_field_name))

       ! The actual dimlist is this                    XDim,YDim,Height
       ! Need to reshape it so that the order becomes: Height,YDim,XDim
       if ( DEEBUG) then
         print *, 'Before reshaping'
          call dump(all_the_fields(:,1,1), 'x-slice')
          call dump(all_the_fields(1,:,1), 'y-slice')
          call dump(all_the_fields(1,1,:), 'p-slice')
       end if

       fdims3 = (dims(3))*(i-1)+1
       sdims3 = (dims(3))*i
       if (DEEBUG) print *, 'fsdims3: ',fdims3, sdims3

       t_all_fields(1:nlon,1:nlat,fdims3:sdims3) = all_the_fields(:,:,:)
       if ( DEEBUG) then
          print *, 't_all_fields '
          call dump(t_all_fields(1:65,1,fdims3), 'x-slice')
          call dump(t_all_fields(1,1:65,fdims3), 'y-slice')
          call dump(t_all_fields(1,1,fdims3:nlev), 'p-slice')
       end if

        ! Close grid
        status = he5_gddetach(gd_id(i))
        if(status /= 0) &
           & call announce_error(lcf_where, "failed to detach from grid " &
           & //trim(names(i)))
        call deallocate_test ( all_the_fields, "all_the_fields", moduleName )
    end do

    the_g_data%field(:,:,:,1,1,1) = reshape( t_all_fields, &
      & shape=(/nlev, nlat, nlon/), order=(/3,2,1/))

    if ( DEEBUG ) then
      print *, 'After reshaping'
      call dump(the_g_data%field(1,1,:,1,1,1), 'x-slice')
      call dump(the_g_data%field(1,:,1,1,1,1), 'y-slice')
      call dump(the_g_data%field(:,1,1,1,1,1), 'p-slice')
    end if

    call deallocate_test ( t_all_fields, "t_all_fields", moduleName )

    ! Now read the dims
    do j=1, nlon
        the_g_data%lons(j) = j - 1 + Lon_offset
    end do
    do j=1, nlat
        the_g_data%lats(j) = j - 1 + Lat_offset
    end do

    ! Have not yet figured out how to assign these
    ! Probably will have to read metadata
    the_g_data%lsts = the_g_data%missingValue
    the_g_data%Szas = the_g_data%missingValue
    the_g_data%DateStarts = the_g_data%missingValue
    the_g_data%DateEnds = the_g_data%missingValue

    status = he5_gdclose(file_id)
    if(status /= 0) &
      & call announce_error(lcf_where, "failed to close file " //trim(NCEPFile%name))

  contains 
     subroutine read_the_dim(gd_id, field_name, field_size, values)
       ! Arguments
       integer, intent(in) :: gd_id
       character(len=*), intent(in) :: field_name
       integer, intent(in) :: field_size
       real(r8), dimension(:), pointer :: values
       ! Local variables
       integer :: status
       integer, dimension(1) :: start, stride, edge
       ! Executable
       call allocate_test ( values, field_size, "values", moduleName )
       start = 0                                                             
       stride = 1                                                             
       edge = field_size                                                       
       values = the_g_data%missingValue
       status = he5_gdrdfld(gd_id, trim(field_name), &
         & hsizes(start), hsizes(stride), hsizes(edge), values)                                                     
       if ( status /= 0 ) &
         & call announce_error(lcf_where, "failed to read dim field " &
         & // trim(field_name))
     end subroutine read_the_dim
  end subroutine Read_NCEP_Strat

  ! ---------------------------------------------  ReadGloriaFile  -----
  type ( GriddedData_T) function ReadGloriaFile ( gloriaFile ) result ( grid )
    type(MLSFile_T)                :: GloriaFile
    ! This function reads data in Gloria's gridded data format into memory
    ! Note that Gloria has a very specific format for her gridded data

    ! Local parameters
    integer, parameter :: NOLATS = 45   ! Number of latitudes
    integer, parameter :: NOLONS = 72   ! Number of longitudes
    integer, parameter :: NOHEIGHTS = 22  ! Number of surfaces
    real(r8), parameter :: DLAT = 4.0   ! Latitude spacing
    integer, parameter :: HEADERLEN = 61
    integer, parameter :: OFFSET = 4    ! Extra 4 byte word in file
    integer, parameter :: BUFFERLEN = &
      & headerLen + offset + noHeights*noLats*noLons*4 ! Size of file

    ! Local variables
    integer :: UARSDAY                  ! UARS day for data
    integer :: I                        ! Loop counter
    integer :: LAT                      ! Loop counter
    integer :: LON                      ! Loop counter
    integer :: SURF                     ! Loop counter
    integer :: LUN                      ! Unit number
    integer :: STATUS                   ! Flag
    character (len=headerLen) :: COMMENT
    logical :: EXIST                    ! Flag
    logical :: OPENED                   ! Flag
    real (r4) :: ONEVALUE               ! One element of the grid
    character (len=1), dimension(4) :: TESTCHAR ! For finding our Endian
    integer :: TESTINT              ! For finding our Endian
    character (len=1), dimension(bufferLen) :: BUFFER
    logical :: NEEDSWAP                 ! Flag for endian
        
    ! Executable code

    ! Setup the grid
    call SetupNewGriddedData ( grid, noHeights=noHeights, noLats=noLats, &
      & noLons=noLons, noLsts=1, noSzas=1, noDates=1 )
    grid%equivalentLatitude = .false.
    grid%description = 'gloria'

    ! Now fill the coordinates
    grid%verticalCoordinate = v_is_pressure
    do i = 1, noHeights
      grid%heights(i) = 10.0**(3.0-(i-1)/6.0)
    end do
    do i = 1, noLats
      grid%lats(i) = dLat * ( i - 2.0*noLats/dLat ) - dLat/2.0
    end do
    do i = 1, noLons
      grid%lons(i) = (i-1)*360.0 / noLons
    end do
    grid%lsts(1) = 0.0
    grid%szas(1) = 0.0

    ! Now read Gloria's data file
    ! Perhaps Paul could make this use MLSFiles later? !?????????
    do lun = 20, 99
      inquire ( unit=lun, exist=exist, opened=opened )
      if ( exist .and. .not. opened ) exit
    end do
    open ( unit=lun, file=gloriaFile%name, status='old', form='unformatted', &
      & access='direct', iostat=status, recl=bufferLen )
    if ( status /= 0 ) then
      call output ( 'IOSTAT value is:' )
      call output ( status, advance='yes' )
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Problem opening Gloria format datafile', MLSFile=gloriaFile )
    end if
    ! Read the data
    read ( lun, rec=1 ) buffer
    ! Done with the file
    close ( lun ) 

    ! Get the comment field to work out the day
    write ( comment, '(61A1)' ) buffer ( 1 : headerLen )
    ! Now work out the uars day
    i = index ( comment, 'DAY ' )
    read ( comment(i+4:i+8), * ) uarsDay

    ! Now convert this into TAI time, choose midday for the day
    ! Also ignore leapseconds. uarsDay 1 was 12 Sep 1991.  From
    ! 1 Jan 1993 to 11 Sep 1991 is 478 days, 477.5 puts it at midday
    grid%dateStarts = (uarsDay-477.5) * 60.0 * 60.0 * 24.0
    grid%dateEnds = grid%dateStarts
    ! Now read the field, into a character array first to deal with endian
    ! issues.

    ! Now, what's our endian
    testInt = 1
    testChar = transfer ( testInt, testChar )
    if ( ( ichar(testChar(1)) == 1 ) .and. &
      &  ( ichar(testChar(2)) == 0 ) .and. &
      &  ( ichar(testChar(3)) == 0 ) .and. &
      &  ( ichar(testChar(4)) == 0 ) ) then
      needSwap = .true.
    else if ( ( ichar(testChar(1)) == 0 ) .and. &
      &       ( ichar(testChar(2)) == 0 ) .and. &
      &       ( ichar(testChar(3)) == 0 ) .and. &
      &       ( ichar(testChar(4)) == 1 ) ) then
      needSwap = .false.
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "I'm very confused about the endian of this machine" )
    end if

    ! Put the if statement on the outside to avoid if's in loops
    i = headerLen + offset + 1
    if ( needSwap ) then
      do surf = 1, noHeights
        do lon = 1, noLons
          do lat = 1, noLats
            grid%field ( surf, lat, lon, 1, 1, 1 ) = &
              & transfer ( buffer(i+3:i:-1), oneValue )
            i = i + 4
          end do
        end do
      end do
    else
      ! No byte swapping required, just store it
      do surf = 1, noHeights
        do lon = 1, noLons
          do lat = 1, noLats
            grid%field ( surf, lat, lon, 1, 1, 1 ) = &
              & transfer ( buffer(i:i+3), oneValue )
            i = i + 4
          end do
        end do
      end do
    end if

    ! Patch the missing data
    where ( grid%field(:,:,:,:,:,:) > 0.95e12 )
      grid%field(:,:,:,:,:,:) = grid%missingValue
    end where

  end function ReadGloriaFile

  ! -------------------------------------------  Read_Climatology  -----
  subroutine Read_Climatology ( climFile, root, aprioriData, returnStatus, &
    & mlspcf_l2clim_start, mlspcf_l2clim_end, missingValue, echo_data, &
    & dump_data )
    ! Brief description of program
    ! This subroutine reads a l3ascii file and returns
    ! the data_array to the caller
    use Toggles, only: Gen, Levels, Toggle
    use Trace_m, only: Trace_Begin, Trace_End

    ! Arguments
    type(MLSFile_T)                :: ClimFile
    type (GriddedData_T), dimension(:), pointer :: aprioriData 
    integer, intent(in) :: ROOT        ! Root of the L2CF abstract syntax tree
    integer, intent(out) :: returnStatus ! E.g., FILENOTFOUND
    integer, optional, intent(IN) :: mlspcf_l2clim_start, mlspcf_l2clim_end
    real, optional, intent(IN) :: missingValue
    logical, optional, intent(in) :: echo_data        ! echo climatology quantity name
    logical, optional, intent(in) :: dump_data        ! dump climatology data

    ! Local
    integer, parameter :: version=1
    character (len=FileNameLen)            :: fname   ! Physical file name
    character (len=FileNameLen)            :: path    ! Physical path
    character (len=FileNameLen)            :: ExactName    ! Full and exact path

    ! These determine how much extra to output
    logical, parameter :: debug=.false.
    logical, parameter :: ECHO_GRIDDED_QUANTITIES=.false. ! echo_data overrides
    logical, parameter :: DUMP_GRIDDED_QUANTITIES=.false. ! dump_data overrides

    logical :: end_of_file
    type (GriddedData_T)        :: gddata 
    integer :: ErrType
    logical :: echo
    logical :: dump_now
    integer :: Me = -1             ! String index for trace
    integer:: processCli, CliUnit
    logical :: use_PCF

    ! begin
    call trace_begin ( me, "Read_Climatology", root, &
      & cond=toggle(gen) .and. levels(gen) > 1 )
    ErrType= 0
    returnStatus = FILENOTFOUND
    end_of_file=.false.
    if(present(echo_data)) then
      echo = echo_data
    else
      echo = ECHO_GRIDDED_QUANTITIES
    end if

    if(present(dump_data)) then
      dump_now = dump_data
    else
      dump_now = DUMP_GRIDDED_QUANTITIES
    end if

    use_PCF = present(mlspcf_l2clim_start) &
      & .and. present(mlspcf_l2clim_end) &
      & .and. UseSDPToolkit .and. .false.

    ! use PCF

    if ( use_PCF ) then
      call split_path_name(ClimFile%name, path, fname)
      processCli = GetPCFromRef( fname, mlspcf_l2clim_start, mlspcf_l2clim_end, &
        & .true., ErrType, version, debugOption=debug, ExactName=ExactName )
      if(ErrType /= 0) then
        call announce_error ( root, &
          &"Climatology file name " // trim(fname) // " unmatched in PCF", &
          & 'error number: ', extra_number=ErrType)
        go to 9
      end if
      if ( mls_exists(trim(ExactName)) /= 0 ) then
        ErrType = FILENOTFOUND
        go to 9
      endif
      ErrType = Pgs_io_gen_openF ( processCli, PGSd_IO_Gen_RSeqFrm, 0, &
        cliUnit, version )
    else
      fname = climFile%name
      if ( mls_exists(trim(fname)) /= 0 ) then
        ErrType = FILENOTFOUND
        go to 9
      endif
      ! use Fortran open
      if(debug) call output('opening ' // fname, advance = 'yes')
      call MLS_OpenFile( ClimFile )
      CliUnit = ClimFile%FileID%f_id
      
    end if

    if(debug) then
      if(.not. end_of_file) then
        call output('Not yet eof on io unit', advance = 'yes')
      else
        call output('Starting at eof on io unit', advance = 'yes')
      end if
      call dump( ClimFile )
    end if


    if ( ErrType == PGS_S_SUCCESS ) then
      do while (.not. end_of_file)
        if(debug) call output('reading l3ascii file', advance = 'yes')
        call l3ascii_read_field ( CliUnit, gddata, end_of_file, ErrType)
        select case ( gddata%verticalCoordinate )
        case ( v_is_pressure )
          gddata%heightsUnits = 'hPa'
        case ( v_is_altitude )
          gddata%heightsUnits = 'km'
        case ( v_is_gph )
          gddata%heightsUnits = 'km'
        case ( v_is_theta )
          gddata%heightsUnits = 'K'
        case default
          gddata%heightsUnits = 'hPa'
        end select
        if(ErrType == 0) then
          if(debug) then
            call output('adding to grid database', advance='yes')
            call output('adding grid template to database ', advance='yes')
          end if
          if(echo .or. debug) then
            call output('quantity name ' // gddata%quantityName, advance='yes')
            call output('description ' // gddata%description, advance='yes')
            call output('units ' // gddata%units, advance='yes')
          end if

          if(dump_now) then
            call Dump(gddata, root)
          end if

          ErrType = AddGriddedDataToDatabase(aprioriData, gddata)

          nullify (gddata%heights)
          nullify (gddata%lats)
          nullify (gddata%lons)
          nullify (gddata%lsts)
          nullify (gddata%szas)
          nullify (gddata%dateStarts)
          nullify (gddata%dateEnds)
          nullify (gddata%field)
          if ( present(missingValue) ) then
            if (missingValue /= 0._rgr) &
              & aprioriData(size(aprioriData))%missingValue = missingValue
          end if
          if(debug) call output('Destroying our grid template', advance='yes')
        end if
      end do !(.not. end_of_file)
      ! ok, done with this file and unit number
      if( use_PCF ) then
        ErrType = Pgs_io_gen_CloseF ( CliUnit )
        ! use Fortran close
      else
        if(debug) call output('closing ' // fname, advance = 'yes')
        call MLS_CloseFile( ClimFile )
        ErrType = 0
      end if
      if(ErrType /= 0) then
        call announce_error (ROOT, &
          &"Error closing " // fname, &
          &'error number: ', extra_number=ErrType)
      end if
      returnStatus = ErrType
    else
      call announce_error (ROOT, &
        &"Error opening " // fname, &
        &'error number: ', extra_number=ErrType)
    end if

9   continue
    call trace_end ( "Read_Climatology", root, &
      & cond=toggle(gen) .and. levels(gen) > 1 )

  end subroutine Read_Climatology

  ! ------------------------------------------------  Write_Merra  -----
  subroutine Write_Merra( createGrid, GEOS5File, lcf_where, &
    & the_g_data, GeoDimList, fieldName )

    ! This routine Writes a gmao MERRA file, named something like
    ! MERRA300.prod.assim.inst6_3d_ana_Nv.20050131.hdf (pressure with
    ! fieldname='DELP' or 'T')
    ! based on a filled data
    ! structure appropriate for newer style gmao merra
    ! (For one thing it holds 4-d, not 3-d fields)

    ! FileName and the_g_data are required args
    ! GeoDimList, if present, should be the Dimensions' short names
    ! as a comma-delimited character string in the order:
    ! longitude, latitude, vertical level, time

    ! fieldName, if present, should be the rank 4
    ! like temperature

    ! This file is formatted in the following way:
    ! At each gridded quantity, e.g. EOSGRID,
    ! a rank 4 field, e.g. T, is given with dimensions
    ! 'Time, Height, YDim, XDim'
    ! We'll simply copy it into a single gridded data type
    
    ! This is so similar to Write_geos5 that we ought to combine the two into
    ! a single routine, except the extra timeIndex business would need
    ! careful handling

    ! Arguments
    logical, intent(in)                     :: createGrid
    type(MLSFile_T)                         :: GEOS5file
    integer, intent(IN)                     :: lcf_where    ! node of the lcf that provoked me
    type( GriddedData_T )                   :: the_g_data ! Result
    character (len=*), optional, intent(IN) :: GeoDimList ! Comma-delimited dim names
    character (len=*), optional, intent(IN) :: fieldName ! Name of gridded field

    ! Local Variables
    integer :: file_id, gd_id
    integer, dimension(NENTRIESMAX) :: dims

    integer :: start(4), stride(4), edge(4)
    integer :: status
    character(len=*), parameter    :: GridName = "EOSGRID"
    integer, parameter :: gctp_geo = 0
    real(r8) :: uplft(2), lowrgt(2)
    real(r8) :: projparm(13)
    logical, parameter :: DEEBUG = .false.
    integer, external :: gdopen, gdclose, gdcreate, gdattach, &
      & gddefdim, gddeffld, gddefproj, gddetach, gdwrfld
    ! Executable code
    uplft = (/ -180000000.000000 ,90000000.000000 /)
    lowrgt = (/ 180000000.000000,-90000000.000000 /)
    projparm = 0.0
    dims(1) = the_g_data%noLons
    dims(2) = the_g_data%noLats
    dims(3) = the_g_data%noHeights
    dims(4) = the_g_data%noDates
    edge = dims(1:4)                                                        
    start = 0
    stride = 1
    error = 0
    if ( createGrid ) then
      file_id = gdopen(GEOS5File%Name, DFACC_CREATE)
      gd_id = gdcreate(file_id, gridName, the_g_data%noLons, &
           & the_g_data%noLats, uplft, lowrgt)
      status = gddefdim( gd_Id, 'TDim', the_g_data%noDates )
      status = gddefdim( gd_Id, 'ZDim', the_g_data%noHeights )
      status = gddefdim( gd_Id, 'YDim', the_g_data%noLats )
      status = gddefdim( gd_Id, 'XDim', the_g_data%noLons )
      status = gddefproj( gd_Id, GCTP_GEO, 0, 0, projparm )
      status = gddeffld( gd_Id, 'TIME', 'TDim', DFNT_FLOAT64, 0 )
      status = gddeffld( gd_Id, 'Height', 'ZDim', DFNT_FLOAT32, 0 )
      status = gddeffld( gd_Id, 'YDim', 'YDim', DFNT_FLOAT32, 0 )
      status = gddeffld( gd_Id, 'XDim', 'XDim', DFNT_FLOAT32, 0 )
      status = gdwrfld( gd_Id, 'TIME', start(4), stride(4), edge(4), the_g_data%dateStarts )
      status = gdwrfld( gd_Id, 'Height', start(3), stride(3), edge(3), the_g_data%heights )
      status = gdwrfld( gd_Id, 'YDim', start(2), stride(2), edge(2), the_g_data%Lats )
      status = gdwrfld( gd_Id, 'XDim', start(1), stride(1), edge(1), the_g_data%Lons )
    else
      file_id = gdopen(GEOS5File%Name, DFACC_RDWR)
      gd_id = gdattach(file_id, gridname)
    end if

    if (file_id < 0) then
      call announce_error(lcf_where, "Could not open "// GEOS5File%Name)
    end if
    status = gddeffld( gd_Id, fieldName, 'TIME,Height,YDim,XDim', DFNT_FLOAT32, 0 )
    status = gdwrfld( gd_Id, 'XDim', start(1), stride(1), edge(1), the_g_data%Lons )
    status = gddetach( gd_Id )
    status = gdclose( file_id )

  end subroutine Write_Merra

  !    -- Private --
  ! ------------------------------------------------  ncepFieldNameTohPa  -----
  subroutine ncepFieldNameTohPa ( field_name, pressure )

    ! Snip INTRO and TAIL off field_name
    ! If remainder can be interpreted as a number return that
    ! Otherwise return -999.99
    ! Arguments
    character(len=*), intent(in)    :: field_name
    real(r8), intent(out)           :: pressure
    ! Local variables
    character(len=*), parameter     :: INTRO = 'ISOBARIC LEVEL AT'
    character(len=*), parameter     :: TAIL = '(hPa)'
    character(len=len(field_name))  :: beheaded
    character(len=len(field_name))  :: remainder
    ! Local variables
    integer :: status
    ! Executable
    pressure = undefinedValue !-999.99
    call ReplaceSubString(field_name, beheaded, INTRO, ' ' )
    call ReplaceSubString(beheaded, remainder, TAIL, ' ' )
    ! Now attempt to interpret remainder as a number
    ! print *, 'Attempting to interpret: ', trim(remainder)
    read(remainder, *, iostat=status) pressure
    if ( status /= 0 ) pressure = undefinedValue !-999.99
  end subroutine ncepFieldNameTohPa

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module readGriddedUtils

! $Log$
! Revision 2.6  2019/08/19 22:04:10  pwagner
! We may now freely read the 'H' field from geos59 files
!
! Revision 2.5  2019/01/29 21:48:23  pwagner
! Removed more unused stuff
!
! Revision 2.4  2018/03/14 17:14:25  pwagner
! Unless verbose, omit useless part of error mesg; correct toc desc of subroutines
!
! Revision 2.3  2017/07/10 18:21:33  pwagner
! Print more if verbose
!
! Revision 2.2  2017/03/17 00:08:40  pwagner
! May read DELP fron merra_2 native; quit with apt mesg if cant read climatology
!
! Revision 2.1  2017/03/07 21:18:10  pwagner
! First commit
!
