! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module ncep_dao ! Collections of subroutines to handle TYPE GriddedData_T

  use GriddedData, only: GriddedData_T, rgr, v_is_pressure, &
    & AddGriddedDataToDatabase, Dump, SetupNewGriddedData, NullifyGriddedData
  use HDFEOS, only: HDFE_NENTDIM, &
    & gdopen, gdattach, gddetach, gdclose, gdfldinfo, &
    & gdinqgrid, gdnentries, gdinqdims, gdinqflds
  use Hdf, only: DFACC_RDONLY
  use l3ascii, only: l3ascii_read_field
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: LineLen, NameLen, FileNameLen, R8, R4, I4
  use MLSFiles, only: FILENOTFOUND, &
    & GetPCFromRef, MLS_HDF_VERSION, mls_io_gen_closeF, mls_io_gen_openF, &
    & split_path_name
  use MLSStrings, only: GetStringElement, NumStringElements, Capitalize, &
    & List2Array, LowerCase, ReplaceSubString
  use OUTPUT_M, only: BLANKS, OUTPUT
  use SDPToolkit, only: PGS_S_SUCCESS, PGS_PC_GETREFERENCE, &
    & PGS_IO_GEN_CLOSEF, PGS_IO_GEN_OPENF, PGSD_IO_GEN_RSEQFRM, &
    & PGSd_GCT_INVERSE, &
    & UseSDPToolkit
  use TREE, only: DUMP_TREE_NODE, SOURCE_REF
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error

  implicit none
  private

  private :: Id,ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  character(LEN=130) :: id = & 
    "$Id$"
  character(LEN=*), parameter :: ModuleName="$RCSfile$"
  private :: not_used_here 
  !-----------------------------------------------------------------------------

  public:: source_file_already_read
  public::OBTAIN_CLIM, READ_CLIMATOLOGY, OBTAIN_DAO, Obtain_NCEP
  public::ReadGriddedData, ReadGloriaFile

  integer :: ERROR

  ! First we'll define some global parameters and data types.
  character (len=*), parameter :: DEFAULTDAODIMLIST = 'XDim,YDim,Height,TIME'
  character (len=*), parameter :: DEFAULTDAOFIELDNAME = 'TMPU'
  character (len=*), parameter :: DEFAULTNCEPDIMLIST = 'YDim,XDim'
  character (len=*), parameter :: DEFAULTNCEPGRIDNAME = 'TMP_3'
  character (len=*), parameter :: GEO_FIELD1 = 'Latitude'
  character (len=*), parameter :: GEO_FIELD2 = 'Longitude'
  character (len=*), parameter :: GEO_FIELD3 = 'Height'
  character (len=*), parameter :: GEO_FIELD4 = 'Time'

  character (len=*), parameter :: lit_dao = 'dao'
  character (len=*), parameter :: lit_ncep = 'ncep'
  character (len=*), parameter :: lit_clim = 'clim'
  integer, parameter :: MAXLISTLENGTH=Linelen ! Max length list of grid names
  integer, parameter :: NENTRIESMAX=200 ! Max num of entries

contains

  ! ----------------------------------------------- ReadGriddedData
  subroutine ReadGriddedData(FileName, lcf_where, description, v_type, &
    & the_g_data, returnStatus, &
    & GeoDimList, fieldName, missingValue)

    ! This routine reads a Gridded Data file, returning a filled data
    ! structure and the  appropriate for 'ncep' or 'dao'

    ! FileName and the_g_data are required args
    ! GeoDimList, if present, should be the Dimensions' short names
    ! as a comma-delimited character string in the order:
    ! longitude, latitude, vertical level, time

    ! fieldName, if present, should be the rank 3 or higher object
    ! like temperature

    ! Arguments
    character (LEN=*), intent(IN) :: FileName ! Name of the file containing the grid(s)
    integer, intent(IN) :: lcf_where    ! node of the lcf that provoked me
    integer, intent(IN) :: v_type       ! vertical coordinate; an 'enumerated' type
    type( GriddedData_T ) :: the_g_data ! Result
    character (LEN=*), intent(IN) :: description ! e.g., 'dao'
    integer, intent(out) :: returnStatus ! E.g., FILENOTFOUND
    character (LEN=*), optional, intent(IN) :: GeoDimList ! Comma-delimited dim names
    character (LEN=*), optional, intent(IN) :: fieldName ! Name of gridded field
    real(rgr), optional, intent(IN) :: missingValue

    ! Local Variables
    character ( len=NameLen) :: my_description   ! In case mixed case
    logical, parameter :: DEEBUG = .false.
    ! Executable code
    my_description = lowercase(description)
    if ( DEEBUG ) print *, 'Reading ' // trim(my_description) // ' data'

    returnStatus = mls_hdf_version(FileName)
    if ( returnStatus == FILENOTFOUND ) then
      call SetupNewGriddedData ( the_g_data, empty=.true. )
      return
    else
      returnStatus = 0
    endif
    ! According to the kinds of gridded data files we can read
    select case ( trim(my_description) )
    case ('clim')
      call announce_error(lcf_where, 'READGriddedData called with climatology' &
        & // ' description; appropriate for ncep/dao only')
    case ('dao')
      call Read_dao(FileName, lcf_where, v_type, &
        & the_g_data, GeoDimList, fieldName, missingValue)
      if ( DEEBUG ) then
        print *, '(Returned from read_dao)'
        print *, 'Quantity Name ' // trim(the_g_data%QuantityName)
        print *, 'Description   ' // trim(the_g_data%description)
        print *, 'Units         ' // trim(the_g_data%units)
      endif
    case ('ncep')
      call Read_ncep(FileName, lcf_where, v_type, &
        & the_g_data, GeoDimList, fieldName, missingValue)
      if ( DEEBUG ) then
        print *, '(Returned from read_ncep)'
        print *, 'Quantity Name ' // trim(the_g_data%QuantityName)
        print *, 'Description   ' // trim(the_g_data%description)
        print *, 'Units         ' // trim(the_g_data%units)
      endif
    case ('olddao', 'oldncep')
      call Read_old(FileName, lcf_where, my_description(4:), v_type, &
        & the_g_data, GeoDimList, fieldName)
    case default
      call announce_error(lcf_where, 'READGriddedData called with unknown' &
        & // ' description: ' // trim(my_description))
    end select

  end subroutine ReadGriddedData

  ! ----------------------------------------------- Read_dao
  subroutine Read_dao(FileName, lcf_where, v_type, &
    & the_g_data, GeoDimList, fieldName, missingValue)
    use Dump_0, only: Dump

    ! This routine reads a dao file, named something like
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
    character (LEN=*), intent(IN) :: FileName ! Name of the file containing the grid(s)
    integer, intent(IN) :: lcf_where    ! node of the lcf that provoked me
    integer, intent(IN) :: v_type       ! vertical coordinate; an 'enumerated' type
    type( GriddedData_T ) :: the_g_data ! Result
    character (LEN=*), optional, intent(IN) :: GeoDimList ! Comma-delimited dim names
    character (LEN=*), optional, intent(IN) :: fieldName ! Name of gridded field
    real(rgr), optional, intent(IN) :: missingValue

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
    integer, parameter :: MAXNAMELENGTH=NameLen		! Max length of grid name
    character (len=MAXNAMELENGTH) :: gridname, actual_field_name
    integer, dimension(NENTRIESMAX) :: dims, rank, numberTypes
    integer                        :: our_rank, numberType

    integer :: start(4), stride(4), edge(4)
    integer :: status
    !                                  These start out initialized to one
    integer                        :: nlon=1, nlat=1, nlev=1, ntime=1
    integer, parameter             :: i_longitude=1
    integer, parameter             :: i_latitude=i_longitude+1
    integer, parameter             :: i_vertical=i_latitude+1
    integer, parameter             :: i_time=i_vertical+1
    integer, external :: GDRDFLD
    logical, parameter :: COUNTEMPTY=.true.
    real(r4), parameter :: FILLVALUE = 1.e15
    real(r4), dimension(:,:,:,:), pointer :: all_the_fields
    real(r8), dimension(:), pointer :: dim_field
    logical, parameter :: DEEBUG = .false.
    ! Executable code
    ! Find list of grid names on this file (This has been core dumping on me)
    if(DEEBUG) print *, 'About to find grid list of file ', trim(FileName)
    inq_success = gdinqgrid(FileName, gridlist, strbufsize)
    if (inq_success < 0) then
      call announce_error(lcf_where, "Could not inquire gridlist "// trim(FileName))
    end if
    if(DEEBUG) print *, 'grid list ', trim(gridlist)

    error = 0
    file_id = gdopen(FileName, DFACC_RDONLY)

    if (file_id < 0) then
      call announce_error(lcf_where, "Could not open "// FileName)
    end if

    ! Find grid name corresponding to the GRIDORDER'th one
    ngrids = NumStringElements(gridlist, COUNTEMPTY)

    if(ngrids <= 0) then
      call announce_error(lcf_where, "NumStringElements of gridlist <= 0")
    elseif(ngrids /= inq_success) then
      call announce_error(lcf_where, "NumStringElements of gridlist /= inq_success")
    elseif(ngrids < GRIDORDER) then
      call announce_error(lcf_where, "NumStringElements of gridlist < GRIDORDER")
    endif

    call GetStringElement(gridlist, gridname, GRIDORDER, COUNTEMPTY)

    gd_id = gdattach(file_id, gridname)
    if (gd_id < 0) then
      call announce_error(lcf_where, "Could not attach "//trim(gridname))
    end if

    ! Now find dimsize(), dimname(), etc.
    nentries = gdnentries(gd_id, HDFE_NENTDIM, strbufsize)

    if(nentries <= 0) then
      call announce_error(lcf_where, "nentries of gd_id <= 0")
    elseif(nentries > NENTRIESMAX) then
      call announce_error(lcf_where, "nentries of gd_id > NENTRIESMAX")
    endif

    ndims = gdinqdims(gd_id, dimlist, dims)

    if(ndims <= 0) then
      call announce_error(lcf_where, "ndims of gd_id <= 0")
    elseif(ndims > NENTRIESMAX) then
      call announce_error(lcf_where, "ndims of gd_id > NENTRIESMAX")
    endif

    nfields = gdinqflds(gd_id, fieldlist, rank, numberTypes)

    if(nfields <= 0) then
      call announce_error(lcf_where, "nfields of gd_id <= 0")
    elseif(nfields > NENTRIESMAX) then
      call announce_error(lcf_where, "nfields of gd_id > NENTRIESMAX")
    endif

    if(.not. CASESENSITIVE) then
      fieldlist = Capitalize(fieldlist)
    endif

    if(present(fieldName)) then
      actual_field_name=fieldName
    else
      actual_field_name=DEFAULTDAOFIELDNAME
    endif

    actual_dim_list = ' '
    if(present(GeoDimList)) then
      actual_dim_list=GeoDimList
    endif
    if ( actual_dim_list == ' ' ) then
      actual_dim_list = DEFAULTDAODIMLIST
    endif
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

    call nullifyGriddedData ( the_g_data ) ! for Sun's still useless compiler
    ! Setup the grid
    call SetupNewGriddedData ( the_g_data, noHeights=nlev, noLats=nlat, &
      & noLons=nlon, noLsts=1, noSzas=1, noDates=ntime, missingValue=FILLVALUE )
      ! & noLons=nlon, noLsts=ntime, noSzas=1, noDates=1, missingValue=FILLVALUE )
    if(DEEBUG) print *, '(Again) our quantity name ', the_g_data%quantityName
    if(DEEBUG) print *, 'our description ', the_g_data%description
    if(DEEBUG) print *, 'our units ', the_g_data%units
    allocate(all_the_fields(dims(1), dims(2), dims(3), dims(4)), stat=status)
    all_the_fields = the_g_data%missingValue
    if ( status /= 0 ) &
        & call announce_error(lcf_where, "failed to allocate field_data")
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
    endif
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
    endif
    deallocate(all_the_fields)
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
    deallocate(dim_field)        ! Before leaving, some light housekeeping
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
      & call announce_error(lcf_where, "failed to close file " //trim(FileName))
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
         if ( associated(values) ) then
           deallocate(values, stat=status)
           if ( status /= 0 ) &
             & call announce_error(lcf_where, "failed to deallocate dim field")
         endif
         allocate(values(field_size), stat=status)
         if ( status /= 0 ) &
           & call announce_error(lcf_where, "failed to allocate dim field")
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
  end subroutine Read_dao

  ! ----------------------------------------------- Read_ncep
  subroutine Read_ncep(FileName, lcf_where, v_type, &
    & the_g_data, GeoDimList, gridName, missingValue)
    use Dump_0, only: Dump

    ! This routine reads a ncep file, named something like
    ! 6207ea2e-1dd2-11b2-b61e-ae069991db9a.hdfeos 
    ! returning a filled data
    ! structure appropriate for newer style ncep

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
    character (LEN=*), intent(IN) :: FileName ! Name of the file of the grid(s)
    integer, intent(IN) :: lcf_where    ! node of the lcf that provoked me
    integer, intent(IN) :: v_type       ! vertical coordinate
    type( GriddedData_T ) :: the_g_data ! Result
    character (LEN=*), optional, intent(IN) :: GeoDimList ! E.g., 'X.Y,..'
    character (LEN=*), optional, intent(IN) :: gridName ! Name of grid
    real(rgr), optional, intent(IN) :: missingValue

    ! Local Variables
    integer :: file_id, gd_id
    integer :: inq_success
    integer :: nentries, ngrids, ndims, nfields
    integer :: strbufsize
    character(len=NameLen) :: my_gridname
    integer, parameter :: MAXLISTLENGTH=10*Linelen ! Max length list grid names
    character (len=MAXLISTLENGTH), dimension(1) :: dimlists
    character (len=MAXLISTLENGTH) :: fieldlist
    character (len=MAXLISTLENGTH) :: gridlist
    character (len=MAXLISTLENGTH) :: dimlist
    integer, parameter :: MAXNAMELENGTH=NameLen		! Max length of grid name
    character (len=MAXNAMELENGTH) :: actual_field_name
    integer, dimension(NENTRIESMAX) :: dims, rank, numberTypes
    integer                        :: our_rank, numberType

    integer :: start(2), stride(2), edge(2)
    integer :: status
    integer :: field
    !                                  These start out initialized to one
    integer                        :: nlon=1, nlat=1, nlev=1, ntime=1
    integer, parameter             :: i_longitude=1
    integer, parameter             :: i_latitude=i_longitude+1
    integer, parameter             :: i_vertical=i_latitude+1
    integer, parameter             :: i_time=i_vertical+1
    integer, external :: GDRDFLD
    logical, parameter :: COUNTEMPTY = .true.
    logical, parameter :: USEPROJECTFORDIMS = .false. ! This doesn't work, yet
    real(r4), dimension(:), pointer     :: pressures
    real(r4), dimension(:,:), pointer   :: field_data
    real(r4), dimension(:,:,:), pointer :: all_the_fields, now_the_fields
    real(r8) :: pressure
    real(r4), dimension(:), pointer :: dim_field
    integer :: xdimsize, ydimsize
    real(r8), dimension(2) :: upLeft, lowRight
    real(r8), dimension(15) :: projParams
    integer :: projcode, zonecode, spherecode, directioncode
    integer :: i
    integer, external :: gdgridinfo, GDprojinfo, PGS_GCT_Init
    integer, parameter :: Lon_offset = -180  ! Longitude start in degrees
    integer, parameter :: Lat_offset = -90   ! Latitude start in degrees
    real(r4), parameter :: FILLVALUE = 0.e0
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
    file_id = gdopen(FileName, DFACC_RDONLY)
    if(DEEBUG) print *, 'fileID: ', file_id

    if (file_id < 0) then
      call announce_error(lcf_where, "Could not open "// FileName)
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
      elseif(nentries > NENTRIESMAX) then
        call announce_error(lcf_where, "nentries of gd_id > NENTRIESMAX", &
          & 'nentries:', nentries)
      end if

      ndims = gdinqdims(gd_id, dimlist, dims)

      if(ndims <= 0) then
        call announce_error(lcf_where, "ndims of gd_id <= 0", &
          & 'ndims:', ndims)
      elseif(ndims > NENTRIESMAX) then
        call announce_error(lcf_where, "ndims of gd_id > NENTRIESMAX", &
          & 'ndims:', ndims)
      endif
    endif

    status = gdgridinfo(gd_id, xdimsize, ydimsize, &
      & upLeft, lowRight )

    if(status /= 0) then
      call announce_error(lcf_where, "failed to obtain gridinfo", &
        & 'status:', status)
    endif
    if(DEEBUG) print *, 'xdimsize: ', xdimsize
    if(DEEBUG) print *, 'ydimsize: ', ydimsize
    if(DEEBUG) print *, 'upLeft  : ', upLeft
    if(DEEBUG) print *, 'lowRight: ', lowRight

    nfields = gdinqflds(gd_id, fieldlist, rank, numberTypes)
    if(DEEBUG) print *, 'nfields: ', nfields

    if(nfields <= 0) then
      call announce_error(lcf_where, "nfields of gd_id <= 0", &
        & 'nfields:', nfields)
      return
    elseif(nfields > NENTRIESMAX) then
      call announce_error(lcf_where, "nfields of gd_id > NENTRIESMAX", &
        & 'nfields:', nfields)
      return
    endif

    the_g_data%noHeights = 0
    allocate(pressures(nfields), stat=status)
    if ( status /= 0 ) &
        & call announce_error(lcf_where, "failed to allocate pressures")
    ! Loop over data fields
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
      endif
      allocate(field_data(dims(1), dims(2)), stat=status)
      if ( status /= 0 ) &
        & call announce_error(lcf_where, "failed to allocate field_data")
      
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
        allocate(all_the_fields(1, dims(1), dims(2)))
        all_the_fields(1, :, :) = field_data
      else
        allocate(now_the_fields(the_g_data%noHeights, dims(1), dims(2)))
        now_the_fields(1:the_g_data%noHeights - 1, :, :) = all_the_fields
        now_the_fields(the_g_data%noHeights, :, :) = field_data
        deallocate(all_the_fields)
        all_the_fields => now_the_fields
      endif
      deallocate(field_data)
    enddo
    the_g_data%noLsts = 0
    the_g_data%units = 'K'
    if(DEEBUG) print *, 'our quantity name ', the_g_data%quantityName
    if(DEEBUG) print *, 'our description ', the_g_data%description
    if(DEEBUG) print *, 'our units ', the_g_data%units

    call nullifyGriddedData ( the_g_data ) ! for Sun's still useless compiler
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
    endif
    the_g_data%field(:,:,:,1,1,1) = reshape( all_the_fields, &
      & shape=(/nLev, nlat, nlon/), order=(/1, 3, 2/) &
      & )
    if(DEEBUG) then
      print *, 'After reshaping'
      call dump(the_g_data%field(1,1,:,1,1,1), 'x-slice')
      call dump(the_g_data%field(1,:,1,1,1,1), 'y-slice')
      call dump(the_g_data%field(:,1,1,1,1,1), 'p-slice')
    endif
    the_g_data%heights = pressures(1:the_g_data%noHeights)
    deallocate(all_the_fields)
    deallocate(pressures)

    ! Insert code to transform upLeft and LowRight
    ! info into lats and lons
    if ( USEPROJECTFORDIMS ) then
      call announce_error(lcf_where, "read_ncep unable to use projection")
    else
    ! Now read the dims
    !  But .. they aren't there!
      !nullify(dim_field)
      !call read_the_dim(gd_id, 'XDim', dims(1), dim_field)
      !the_g_data%lons = dim_field
      !call read_the_dim(gd_id, 'YDim', dims(2), dim_field)
      !the_g_data%lats = dim_field
      !deallocate(dim_field)        ! Before leaving, some light housekeeping
      do i=1, xdimsize
        the_g_data%lons(i) = i - 1 + Lon_offset
      enddo
      do i=1, ydimsize
        the_g_data%lats(i) = i - 1 + Lat_offset
      enddo
    endif
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
      & call announce_error(lcf_where, "failed to close file " //trim(FileName))
    contains 
       subroutine read_the_dim(gd_id, field_name, field_size, values)
         ! Arguments
         integer, intent(in) :: gd_id
         character(len=*), intent(in) :: field_name
         integer, intent(in) :: field_size
         real(r4), dimension(:), pointer :: values
         ! Local variables
         integer :: status
         integer, dimension(1) :: start, stride, edge
         ! Executable
         if ( associated(values) ) then
           deallocate(values, stat=status)
           if ( status /= 0 ) &
             & call announce_error(lcf_where, "failed to deallocate dim field")
         endif
         allocate(values(field_size), stat=status)
         if ( status /= 0 ) &
           & call announce_error(lcf_where, "failed to allocate dim field")
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
  end subroutine Read_ncep

  ! ----------------------------------------------- Read_old
  subroutine Read_old(FileName, lcf_where, description, v_type, &
    & the_g_data, GeoDimList, fieldName)

    ! This routine reads a Gridded Data file, returning a filled data
    ! structure appropriate for older style 'ncep' or 'dao'

    ! FileName and the_g_data are required args
    ! GeoDimList, if present, should be the Dimensions' short names
    ! as a comma-delimited character string in the order:
    ! longitude, latitude, vertical level, time

    ! fieldName, if present, should be the rank 3 or higher object
    ! like temperature

    ! Arguments
    character (LEN=*), intent(IN) :: FileName ! Name of the file containing the grid(s)
    integer, intent(IN) :: lcf_where    ! node of the lcf that provoked me
    integer, intent(IN) :: v_type       ! vertical coordinate; an 'enumerated' type
    type( GriddedData_T ) :: the_g_data ! Result
    character (LEN=*), intent(IN) :: description ! e.g., 'dao'
    character (LEN=*), optional, intent(IN) :: GeoDimList ! Comma-delimited dim names
    character (LEN=*), optional, intent(IN) :: fieldName ! Name of gridded field

    ! Local Variables
    integer :: file_id, gd_id
    integer :: inq_success
    integer :: nentries, ngrids, ndims, nfields
    integer :: strbufsize
    logical,  parameter       :: CASESENSITIVE = .false.
    integer, parameter :: GRIDORDER=1   ! What order grid written to file
    integer, parameter :: MAXLISTLENGTH=Linelen ! Max length list of grid names
    integer, parameter :: NENTRIESMAX=20 ! Max num of entries
    character (len=MAXLISTLENGTH) :: gridlist
    character (len=MAXLISTLENGTH) :: dimlist, actual_dim_list
    character (len=MAXLISTLENGTH), dimension(1) :: dimlists
    character (len=MAXLISTLENGTH) :: fieldlist
    integer, parameter :: MAXNAMELENGTH=NameLen		! Max length of grid name
    character (len=MAXNAMELENGTH) :: gridname, actual_field_name
    integer, dimension(NENTRIESMAX) :: dims, rank, numberTypes
    integer                        :: our_rank, numberType

    !                                  These start out initialized to one
    integer                        :: nlon=1, nlat=1, nlev=1, ntime=1
    integer, parameter             :: i_longitude=1
    integer, parameter             :: i_latitude=i_longitude+1
    integer, parameter             :: i_vertical=i_latitude+1
    integer, parameter             :: i_time=i_vertical+1
    integer, external :: GDRDFLD
    logical, parameter :: COUNTEMPTY=.true.
    integer :: status
    logical, parameter :: DEEBUG = .false.
    ! Executable code
    ! Find list of grid names on this file (This has been core dumping on me)
    if(DEEBUG) print *, 'About to find grid list'
    if(DEEBUG) print *, "proceed (yes) or (no)"
    ! read *, gridlist
    ! if ( index(gridlist, 'y') < 1 ) stop
    inq_success = gdinqgrid(FileName, gridlist, strbufsize)
    if (inq_success < 0) then
      call announce_error(lcf_where, "Could not inquire gridlist "// FileName)
    end if
    if(DEEBUG) print *, 'grid list ', trim(gridlist)

    error = 0
    file_id = gdopen(FileName, DFACC_RDONLY)

    if (file_id < 0) then
      call announce_error(lcf_where, "Could not open "// FileName)
    end if

    ! Find grid name corresponding to the GRIDORDER'th one
    ngrids = NumStringElements(gridlist, COUNTEMPTY)

    if(ngrids <= 0) then
      call announce_error(lcf_where, "NumStringElements of gridlist <= 0")
    elseif(ngrids /= inq_success) then
      call announce_error(lcf_where, "NumStringElements of gridlist /= inq_success")
    elseif(ngrids < GRIDORDER) then
      call announce_error(lcf_where, "NumStringElements of gridlist < GRIDORDER")
    endif

    call GetStringElement(gridlist, gridname, GRIDORDER, COUNTEMPTY)

    gd_id = gdattach(file_id, gridname)
    if (gd_id < 0) then
      call announce_error(lcf_where, "Could not attach "//FileName)
    end if

    ! Now find dimsize(), dimname(), etc.
    nentries = gdnentries(gd_id, HDFE_NENTDIM, strbufsize)

    if(nentries <= 0) then
      call announce_error(lcf_where, "nentries of gd_id <= 0")
    elseif(nentries > NENTRIESMAX) then
      call announce_error(lcf_where, "nentries of gd_id > NENTRIESMAX")
    endif

    ndims = gdinqdims(gd_id, dimlist, dims)

    if(ndims <= 0) then
      call announce_error(lcf_where, "ndims of gd_id <= 0")
    elseif(ndims > NENTRIESMAX) then
      call announce_error(lcf_where, "ndims of gd_id > NENTRIESMAX")
    endif

    nfields = gdinqflds(gd_id, fieldlist, rank, numberTypes)

    if(nfields <= 0) then
      call announce_error(lcf_where, "nfields of gd_id <= 0")
    elseif(nfields > NENTRIESMAX) then
      call announce_error(lcf_where, "nfields of gd_id > NENTRIESMAX")
    endif

    if(.not. CASESENSITIVE) then
      fieldlist = Capitalize(fieldlist)
    endif

    if(present(fieldName)) then
      actual_field_name=fieldName
    else
      actual_field_name=DEFAULTDAOFIELDNAME
    endif

    if(present(GeoDimList)) then
      actual_dim_list=GeoDimList
    else
      actual_dim_list=GEO_FIELD1 // ',' // &
        & GEO_FIELD2 // ',' // &
        & GEO_FIELD3 // ',' // &
        & GEO_FIELD4
    endif

    ! Now find the rank of our field
    inq_success = gdfldinfo(gd_id, trim(actual_field_name), our_rank, dims, &
      & numbertype, dimlists(1))

    dimlist = trim(dimlists(1))

    nlon = dims(1)
    nlat = dims(2)
    nlev = dims(3)
    ntime = dims(4)

    the_g_data%quantityName = actual_field_name
    the_g_data%description = description
    the_g_data%verticalCoordinate = v_type

    the_g_data%noLons = nlon
    the_g_data%noLats = nlat
    the_g_data%noHeights = nlev
    the_g_data%noLsts = ntime

    ! Close grid
    status = gddetach(gd_id)
    if(status /= 0) &
      & call announce_error(lcf_where, "failed to detach from grid " // &
      & trim(gridname))
    status = gdclose(file_id)
    if(status /= 0) &
      & call announce_error(lcf_where, "failed to close file " // &
      & trim(FileName))
  end subroutine Read_old

  ! --------------------------------------------- ReadGloriaFile -------
  type ( GriddedData_T) function ReadGloriaFile ( filename ) result ( grid )
    character (len=*), intent(in) :: FILENAME
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
    integer(i4) :: DUMMY                ! For record length words
    character (len=headerLen) :: COMMENT
    logical :: EXIST                    ! Flag
    logical :: OPENED                   ! Flag
    real (r4) :: ONEVALUE               ! One element of the grid
    character (len=1), dimension(4) :: TESTCHAR ! For finding our Endian
    integer(i4) :: TESTINT              ! For finding our Endian
    character (len=1), dimension(bufferLen) :: BUFFER
    logical :: NEEDSWAP                 ! Flag for endian
        
    ! Executable code

    call nullifyGriddedData ( grid ) ! for Sun's still useless compiler
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
    open ( unit=lun, file=filename, status='old', form='unformatted', &
      & access='direct', iostat=status, recl=bufferLen )
    if ( status /= 0 ) then
      call output ( 'IOSTAT value is:' )
      call output ( status, advance='yes' )
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Problem opening Gloria format datafile' )
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

  ! --------------------------------------------- Obtain_Clim ----------
  subroutine OBTAIN_CLIM ( aprioriData, root, &
    & mlspcf_l2clim_start, mlspcf_l2clim_end )

    ! An atavism--
    ! a throwback to when ncep files were opened
    ! independently of being required by the lcf

    !Arguments 
    type (GriddedData_T), dimension(:), pointer :: aprioriData 
    ! Input a priori database
    integer, intent(in) :: root        ! Root of the L2CF abstract syntax tree
    integer, intent(in) :: mlspcf_l2clim_start, mlspcf_l2clim_end

    !Local Variables

    type (GriddedData_T):: qty
    integer:: CliUnit, processCli, returnStatus, version

    logical :: end_of_file = .false.

   error = 0
   
   if( .not. UseSDPToolkit ) then
      call announce_error(root, &
      & 'Detached from toolkit--climatology files must be opened via lcf')
      return
   endif

    do CliUnit = mlspcf_l2clim_start, mlspcf_l2clim_end

      !     Open one Climatology file as a generic file for reading
      version = 1
      returnStatus = Pgs_io_gen_openF ( CliUnit, PGSd_IO_Gen_RSeqFrm, 0, &
        processCli, version )
      if ( returnStatus == PGS_S_SUCCESS ) then

        do while (.not. end_of_file)

          call l3ascii_read_field ( processCli, qty, end_of_file)
          returnStatus = AddGriddedDataToDatabase(aprioriData, qty)

          nullify (qty%lats)
          nullify (qty%lons)
          nullify (qty%lsts)
          nullify (qty%szas)
          nullify (qty%dateStarts)
          nullify (qty%dateEnds)
          nullify (qty%field)
          ! No, this is a bad idea (according to njl)
          !        call DestroyGridTemplateContents ( qty )

        end do !(.not. end_of_file)

        end_of_file = .false.

      end if

    end do ! CliUnit = mlspcf_l2clim_start, mlspcf_l2clim_end

    return

  end subroutine OBTAIN_CLIM

  ! --------------------------------------------------  Read_Climatology
  subroutine READ_CLIMATOLOGY ( input_fname, root, aprioriData, &
    & mlspcf_l2clim_start, mlspcf_l2clim_end, missingValue, echo_data, &
    & dump_data )
    ! Brief description of program
    ! This subroutine reads a l3ascii file and returns
    ! the data_array to the caller

    ! Arguments
    character*(*), intent(in) :: input_fname			! Physical file name
    type (GriddedData_T), dimension(:), pointer :: aprioriData 
    integer, intent(in) :: ROOT        ! Root of the L2CF abstract syntax tree
    integer, optional, intent(IN) :: mlspcf_l2clim_start, mlspcf_l2clim_end
    real, optional, intent(IN) :: missingValue
    logical, optional, intent(in) :: echo_data        ! echo climatology quantity name
    logical, optional, intent(in) :: dump_data        ! dump climatology data

    ! Local
    integer, parameter :: version=1
    character (LEN=FileNameLen)            :: fname   ! Physical file name
    character (LEN=FileNameLen)            :: path	   ! Physical path

    ! These determine how much extra to output
    logical, parameter :: debug=.false.
    logical, parameter :: ECHO_GRIDDED_QUANTITIES=.false.	! echo_data overrides
    logical, parameter :: DUMP_GRIDDED_QUANTITIES=.false.	! dump_data overrides

    logical :: end_of_file
    type (GriddedData_T)        :: gddata 
    integer :: ErrType
    logical :: echo
    logical :: dump_now
    integer:: processCli, CliUnit, record_length
    logical :: use_PCF

    ! begin
    end_of_file=.false.
    if(present(echo_data)) then
      echo = echo_data
    else
      echo = ECHO_GRIDDED_QUANTITIES
    endif

    if(present(dump_data)) then
      dump_now = dump_data
    else
      dump_now = DUMP_GRIDDED_QUANTITIES
    endif

    use_PCF = present(mlspcf_l2clim_start) &
      & .and. present(mlspcf_l2clim_end) &
      & .and. UseSDPToolkit

    ! use PCF

    if ( use_PCF ) then
      call split_path_name(input_fname, path, fname)
      processCli = GetPCFromRef(fname, mlspcf_l2clim_start, mlspcf_l2clim_end, &
        & .true., ErrType, version, debugOption=debug)
      if(ErrType /= 0) then
        !    CALL MLSMessage (MLSMSG_Error, ModuleName, &
        !              &"Climatology file name unmatched in PCF")
        call announce_error (ROOT, &
          &"Climatology file name " // trim(fname) // " unmatched in PCF", &
          & 'error number: ', extra_number=ErrType)
        return
      endif
      ErrType = Pgs_io_gen_openF ( processCli, PGSd_IO_Gen_RSeqFrm, 0, &
        cliUnit, version )
    else
      fname = input_fname
      ! use Fortran open
      if(debug) call output('opening ' // fname, advance = 'yes')
      CliUnit = mls_io_gen_openF ( 'open', .true., ErrType, &
	& record_length, PGSd_IO_Gen_RSeqFrm, FileName=fname)
    endif

    if(debug) then
      if(.not. end_of_file) then
        call output('Not yet eof on io unit', advance = 'yes')
      else
        call output('Starting at eof on io unit', advance = 'yes')
      endif
    endif


    if ( ErrType == PGS_S_SUCCESS ) then
      do while (.not. end_of_file)
        if(debug) call output('reading l3ascii file', advance = 'yes')
        call l3ascii_read_field ( CliUnit, gddata, end_of_file, ErrType)
        if(ErrType == 0) then
          if(debug) then
            call output('adding to grid database', advance='yes')
            call output('adding grid template to database ', advance='yes')
          endif
          if(echo .or. debug) then
            call output('quantity name ' // gddata%quantityName, advance='yes')
            call output('description ' // gddata%description, advance='yes')
            call output('units ' // gddata%units, advance='yes')
          endif

          if(dump_now) then
            call Dump(gddata, root)
          endif

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
          endif
          if(debug) call output('Destroying our grid template', advance='yes')
        endif
      end do !(.not. end_of_file)
      ! ok, done with this file and unit number
      if( use_PCF ) then
        ErrType = Pgs_io_gen_CloseF ( CliUnit )
	! use Fortran close
      else
        if(debug) call output('closing ' // fname, advance = 'yes')
        ErrType = mls_io_gen_CloseF ('close', CliUnit )
      endif
      if(ErrType /= 0) then
        call announce_error (ROOT, &
          &"Error closing " // fname, &
          &'error number: ', extra_number=ErrType)
      endif
    else
      call announce_error (ROOT, &
        &"Error opening " // fname, &
        &'error number: ', extra_number=ErrType)
    endif

  end subroutine READ_CLIMATOLOGY

  ! -------------------------------------------------  Obtain_DAO
  subroutine OBTAIN_DAO ( aprioriData, root, &
    & mlspcf_l2dao_start, mlspcf_l2dao_end )

    ! An atavism--
    ! a throwback to when ncep files were opened
    ! independently of being required by the lcf

    type (GriddedData_T), dimension(:), pointer :: aprioriData 
    ! Input a priori database
    integer, intent(in) :: ROOT        ! Root of L2CF abstract syntax tree
    integer, intent(IN) :: mlspcf_l2dao_start, mlspcf_l2dao_end

    ! Local Variables

    !    real(R8) :: data_array(XDIM, YDIM, ZDIM)
    integer :: DAOFileHandle, DAO_Version
    character (LEN=132) :: DAOphysicalFilename
    type (GriddedData_T):: qty
    integer :: returnStatus
    !   integer :: sd_id
    character (LEN=80) :: vname

    !    ALLOCATE (data_array(XDIM, YDIM, ZDIM), stat=returnStatus)

   error = 0
   
   if( .not. UseSDPToolkit ) then
      call announce_error(root, &
      & 'Detached from toolkit--DAO files must be opened via lcf')
      return
   endif

    DAO_Version = 1
    vname = "TMPU" ! for now

    ! Get the DAO file name from the PCF
    do DAOFileHandle = mlspcf_l2dao_start, mlspcf_l2dao_end
      returnStatus = Pgs_pc_getReference ( DAOFileHandle, DAO_Version, &
        DAOphysicalFilename )
      if ( returnStatus == PGS_S_SUCCESS ) then
        ! Open the HDF-EOS file and read gridded data
        returnStatus = AddGriddedDataToDatabase(aprioriData, qty)
        !        call read_dao ( DAOphysicalFilename, vname, data_array )
        if(returnStatus > 0) then
          call ReadGriddedData ( DAOphysicalFilename, root, &
            & 'dao', v_is_pressure, aprioriData(returnStatus), returnStatus )
        endif
      end if
    end do ! DAOFileHandle = mlspcf_l2_dao_start, mlspcf_l2_dao_end
  end subroutine Obtain_DAO

  ! ------------------------------------------------  Obtain_NCEP  -----
  subroutine Obtain_NCEP ( aprioriData, root, &
    & mlspcf_l2ncep_start, mlspcf_l2ncep_end )
    ! An atavism--
    ! a throwback to when ncep files were opened
    ! independently of being required by the lcf
    ! Arguments
    type (GriddedData_T), dimension(:), pointer :: aprioriData 
    ! Input a priori database
    integer, intent(in) :: ROOT        ! Root of the L2CF abstract syntax tree
    integer, intent(IN) :: mlspcf_l2ncep_start, mlspcf_l2ncep_end

    ! Local Variables
    integer :: NCEPFileHandle, NCEP_Version
    character (len=132) :: NCEPphysicalFilename
    type (GriddedData_T):: QTY
    integer :: RETURNSTATUS

   error = 0
   
   if( .not. UseSDPToolkit ) then
      call announce_error(root, &
      & 'Detached from toolkit--ncep files must be opened via lcf')
      return
   endif

    NCEP_Version = 1
    !    vname = "TMP_3" ! for now
    ! Get the NCEP file name from the PCF

    do NCEPFileHandle = mlspcf_l2ncep_start, mlspcf_l2ncep_end

      returnStatus = Pgs_pc_getReference ( NCEPFileHandle, NCEP_Version, &
        NCEPphysicalFilename )

      if ( returnStatus == PGS_S_SUCCESS ) then

        ! Open the HDF-EOS file and read gridded data

        returnStatus = AddGriddedDataToDatabase(aprioriData, qty)

        !        call read_ncep ( NCEPphysicalFilename, data_array )
        if(returnStatus > 0) then
          call ReadGriddedData ( NCEPphysicalFilename, root, &
            & 'ncep', v_is_pressure, aprioriData(returnStatus), returnStatus )
        endif

      end if

    end do ! NCEPFileHandle = mlspcf_l2_ncep_start, mlspcf_l2_ncep_end
  end subroutine Obtain_NCEP

  ! --------------------------------  source_file_already_read  -----
  function source_file_already_read(GriddedDataBase, source_file, &
    & field_name ) result (already)
    ! check if source file among those already read to form database
    ! returns .TRUE. if already read, .FALSE. if not or if database is empty

    ! optionally checks that field name is also matched for that partcilar
    ! source file

    ! Arguments
    !   
    type (GriddedData_T), dimension(:), pointer :: GriddedDataBase
    character (LEN=*), intent(in) :: source_file
    character (LEN=*), optional, intent(in) :: field_name
    logical :: already

    ! Local
    integer :: i

    ! Begin
    already = .false.

    if(.not. associated(GriddedDataBase)) then
      return
    elseif(len(source_file) == 0) then
      return
    elseif(size(GriddedDataBase) == 0) then
      return
    endif
    
    call output('Database, filenames:', advance='no')
    call blanks(3)
    call output(GriddedDatabase%sourceFilename, advance='yes')
    
    call output('Database, fieldNames:', advance='no')
    call blanks(3)
    call output(GriddedDatabase%quantityName, advance='yes')
    
    call output('Source file:', advance='no')
    call blanks(3)
    call output(source_file, advance='yes')
    
    call output('field_name:', advance='no')
    call blanks(3)
    call output(field_name, advance='yes')
    
    do i=1, size(GriddedDataBase)
      if(trim(adjustl(source_file)) == trim(adjustl(GriddedDataBase(i)%sourceFileName))) then
        already = .true.
        exit
      endif
    enddo

    if(present(field_name) .and. already) then
      if(trim(adjustl(field_name)) == trim(adjustl(GriddedDataBase(i)%quantityName))) then
        already = .true.
      else
        already = .false.
      endif
    endif
  end function source_file_already_read

  ! ----- utility procedures ----
  
  ! ------------------------------------------------  ncepFieldNameTohPa  -----
  subroutine ncepFieldNameTohPa ( field_name, pressure )

    ! Snip INTRO and TAIL off filed_name
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
    pressure = -999.99
    call ReplaceSubString(field_name, beheaded, INTRO, ' ' )
    call ReplaceSubString(beheaded, remainder, TAIL, ' ' )
    ! Now attempt to interpret remainder as a number
    ! print *, 'Attempting to interpret: ', trim(remainder)
    read(remainder, *, iostat=status) pressure
    if ( status /= 0 ) pressure = -999.99
  end subroutine ncepFieldNameTohPa

  ! ------------------------------------------------  announce_error  -----
  subroutine announce_error ( lcf_where, full_message, &
    & extra_message, extra_number, use_toolkit )

    ! Arguments
    integer, intent(in)    :: lcf_where
    character(LEN=*), intent(in)    :: full_message
    logical, intent(in), optional :: use_toolkit
    character(LEN=*), intent(in), optional    :: extra_message
    integer, intent(in), optional    :: extra_number

    ! Local variables
    logical :: just_print_it
    logical, parameter :: default_output_by_toolkit = .true.

    if(present(use_toolkit)) then
      just_print_it = use_toolkit
    elseif(default_output_by_toolkit) then
      just_print_it = .false.
    else
      just_print_it = .true.
    endif

    if(.not. just_print_it) then
      error = max(error,1)
      call output ( '***** At ' )

      if(lcf_where > 0) then
        call print_source ( source_ref(lcf_where) )
      else
        call output ( '(no lcf node available)' )
      endif

      call output ( ': ' )
      call output ( "The " );
      if(lcf_where > 0) then
        call dump_tree_node ( lcf_where, 0 )
      else
        call output ( '(no lcf tree available)' )
      endif

      call output(" Caused the following error:", advance='yes', &
        & from_where=ModuleName)
      call output(trim(full_message), advance='yes', &
        & from_where=ModuleName)
      if(present(extra_message)) then
        call output('error number ', advance='no')
        call output(extra_message, advance='yes')
      endif
      if(present(extra_number)) then
        call output('error number ', advance='no')
        call output(extra_number, places=9, advance='yes')
      endif
    else
      call output ( '***Error in module ' )
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

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module ncep_dao

! $Log$
! Revision 2.31  2003/05/06 00:31:45  vsnyder
! Delete trailing blanks from two too-long-to-be-standard lines
!
! Revision 2.30  2003/04/04 18:34:22  pwagner
! Sets empty field if FILENOTFOUND
!
! Revision 2.29  2003/04/02 00:14:25  pwagner
! The 4th index is now noDates instead of ntime
!
! Revision 2.28  2003/03/01 00:23:27  pwagner
! missingValue an optional arg to reading GriddedData and climatology
!
! Revision 2.27  2003/02/28 02:27:12  livesey
! Now uses missingValue stuff
!
! Revision 2.26  2003/02/27 21:51:02  pwagner
! Commented out the last prints
!
! Revision 2.25  2003/02/27 18:38:49  pwagner
! Removed some intent(out); Lahey takes perverse delight in resetting such to undefined
!
! Revision 2.24  2003/02/21 21:01:06  pwagner
! Actually uses GeoDimList; filters Fill values
!
! Revision 2.23  2003/02/20 21:23:40  pwagner
! More successful; cant read metadata yet
!
! Revision 2.22  2003/02/19 19:14:16  pwagner
! Many changes; not perfect yet
!
! Revision 2.21  2002/11/22 12:57:59  mjf
! Added nullify routine(s) to get round Sun's WS6 compiler not
! initialising derived type function results.
!
! Revision 2.20  2002/10/08 00:09:13  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.19  2002/02/05 04:13:43  livesey
! Minor bug fix
!
! Revision 2.18  2002/01/23 22:34:37  livesey
! Added ReadGloriaFile functionality
!
! Revision 2.17  2002/01/09 23:48:40  pwagner
! Added toc; each print became call output
!
! Revision 2.16  2001/10/26 23:14:37  pwagner
! Complies with Gridded data dump
!
! Revision 2.15  2001/09/10 23:37:10  livesey
! Tidied up a bit, moved much stuff to GriddedData.f90
!
! Revision 2.14  2001/07/12 22:03:55  livesey
! Some minor ish changes.  Needs an overhaul at some point.
!
! Revision 2.13  2001/06/04 23:57:40  pwagner
! Splits path from l2cf-defined file name before getPCfromRef
!
! Revision 2.12  2001/05/09 23:30:13  pwagner
! Detachable from toolkit
!
! Revision 2.11  2001/04/12 22:04:47  vsnyder
! Improve an error message
!
! Revision 2.10  2001/04/10 20:05:30  livesey
! Tidied up
!
! Revision 2.9  2001/03/30 00:26:19  pwagner
! Added source_file_already_read
!
! Revision 2.8  2001/03/29 00:51:03  pwagner
! AddGridTemplatetoDatabase now works
!
! Revision 2.7  2001/03/28 00:25:14  pwagner
! More changes, but not perfect yet
!
! Revision 2.6  2001/03/27 17:28:31  pwagner
! Can dump gridded database
!
! Revision 2.5  2001/03/24 00:29:32  pwagner
! Now seems to read climatology files better
!
! Revision 2.4  2001/03/21 00:47:29  pwagner
! Changes to READ_CLIMATOLOGY, announce_error
!
! Revision 2.3  2001/03/20 00:42:11  pwagner
! Improved Read_Climatology
!
! Revision 2.2  2001/03/15 21:40:30  pwagner
! Eliminated unused routines from USE statements
!
! Revision 2.1  2001/03/15 21:26:57  pwagner
! Moved non-l3ascii methods from GriddedData here
!
