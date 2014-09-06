! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module ncep_dao ! Collections of subroutines to handle TYPE GriddedData_T

  use allocate_deallocate, only: allocate_test, deallocate_test, bytes, &
    & test_allocate
  use dump_0, only : dump
  use griddeddata, only: griddeddata_t, rgr, v_is_altitude, v_is_gph, &
    & v_is_pressure, v_is_theta, &
    & addgriddeddatatodatabase, dump, setupnewgriddeddata, nullifygriddeddata
  use HDFeos, only: HDFe_nentdim, &
    & gdopen, gdattach, gddetach, gdclose, gdfldinfo, &
    & gdinqgrid, gdnentries, gdinqdims, gdinqflds
  use HDF, only: dfacc_create, dfacc_rdonly, dfacc_rdwr, &
    & dfnt_float32, dfnt_float64
  use highoutput, only: outputnamedvalue
  use l3ascii, only: l3ascii_read_field
  use lexer_core, only: print_source
  use MLScommon, only: linelen, namelen, filenamelen, &
    & undefinedvalue, MLSfile_t
  use MLSfiles, only: filenotfound, HDFversion_5, &
    & dump, getpcfromref, MLS_HDF_version, MLS_openfile, MLS_closefile, &
    & split_path_name, MLS_openfile, MLS_closefile
  use MLSkinds, only: r4, r8
  use MLSmessagemodule, only: MLSmsg_error, MLSmsg_info, MLSmsg_warning, &
    & MLSmessage
  use MLSstrings, only: capitalize, hhmmss_value, lowercase
  use MLSstringlists, only: getstringelement, numstringelements, &
    & list2array, replacesubstring, stringelementnum
  use output_m, only: output
  use SDPtoolkit, only: PGS_s_success, &
    & PGS_io_gen_closef, PGS_io_gen_openf, PGSd_io_gen_rseqfrm, &
    & PGSd_gct_inverse, &
    & useSDPtoolkit
  use tree, only: dump_tree_node, where

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
! read_climatology     read l3ascii-formatted climatology file
! ReadGriddedData      read meteorology file in one of supported desciptions
! WriteGriddedData     write meteorology file in one of supported desciptions
! ReadGloriaFile       read binary-formatted file designed by G. Manney
!
! The supported descriptions are
! geos5_7         Geos5.7.2  gmao files in netcdf4/hdf5 format
! geos5           Geos5x.x   gmao files in hdfeos2/hdf4 format
! dao or gmao     Geos4x.x   gmao files in hdfeos2/hdf4 format
! merra           Geos5x.x   gmao reanalysis files in hdfeos2/hdf4 format
! ncep            GDAS       ncep files in hdfeos2/hdf4 format
! strat           STRAT      ncep combined in hdfeos5/hdf5 format
! === (end of toc) ===

  public:: read_climatology
  public:: ReadGriddedData, ReadGloriaFile
  public:: WriteGriddedData

  integer :: ERROR
  character(len=8) :: LIT_DESCRIPTION

  ! First we'll define some global parameters and data types.
  logical, parameter :: COUNTEMPTY=.true.
  character (len=*), parameter :: DEFAULTDAODIMLIST = 'XDim,YDim,Height,Time'
  character (len=*), parameter :: DEFAULTDAOFIELDNAME = 'TMPU'
  character (len=*), parameter :: DEFAULTGEOS5FIELDNAME = 'T'
  character (len=*), parameter :: DEFAULTNCEPDIMLIST = 'YDim,XDim'
  character (len=*), parameter :: DEFAULTNCEPGRIDNAME = 'TMP_3'
  character (len=*), parameter :: DEFAULTNCEPSTRATFIELDNAME = 'Temperature'
  character (len=*), parameter :: GEO_FIELD1 = 'Latitude'
  character (len=*), parameter :: GEO_FIELD2 = 'Longitude'
  character (len=*), parameter :: GEO_FIELD3 = 'Height'
  character (len=*), parameter :: GEO_FIELD4 = 'Time'

  logical, parameter :: GEOS5MAYBEMERRATOO = .false. ! Causes confusion if true

  character (len=*), parameter :: lit_dao = 'dao'
  character (len=*), parameter :: lit_ncep = 'ncep'
  character (len=*), parameter :: lit_strat = 'strat'
  character (len=*), parameter :: lit_clim = 'clim'
  character (len=*), parameter :: lit_geos5 = 'geos5'
  integer, parameter :: MAXLISTLENGTH=Linelen ! Max length list of grid names
  integer, parameter :: NENTRIESMAX=200 ! Max num of entries
  integer, parameter ::          MAXDS = 1024 ! 500
  integer, parameter ::          MAXSDNAMESBUFSIZE = MAXDS*NAMELEN

  interface ReadGRIDDEDDATA
    module procedure ReadGriddedData_MLSFile
    module procedure ReadGriddedData_Name
  end interface

contains

  ! ----------------------------------------------- ReadGriddedData
  ! This family of routines reads a Gridded Data file, returning a filled data
  ! structure and the  appropriate for the description
  ! which may be one of 'geos5_7', 'geos5', 'gmao', 'dao', 'merra', 'ncep', 'strat'

  ! FileName and the_g_data are required args
  ! GeoDimList should be the Dimensions' short names
  ! as a comma-delimited character string in the order:
  ! longitude, latitude, vertical level, time

  ! fieldName should name the rank 3 or higher object
  ! like temperature

  ! date is needed only for background files
  ! or any other case in which different files come with the date geolocations
  ! unset
  subroutine ReadGriddedData_MLSFile( GriddedFile, lcf_where, description, v_type, &
    & the_g_data, returnStatus, &
    & GeoDimList, fieldNames, missingValue, &
    & date, sumDelp, deferReading, litDescription, verbose )

    use MLSStats1, only: MLSMIN, MLSMAX, MLSMEAN
    use Toggles, only: Gen, Levels, Toggle
    use Trace_m, only: TRACE_BEGIN, TRACE_END

    ! Arguments
    type(MLSFile_T)                         :: GriddedFile
    integer, intent(IN)                     :: lcf_where    ! node of the lcf that provoked me         
    integer, intent(IN)                     :: v_type       ! vertical coordinate; an 'enumerated' type
    type( GriddedData_T )                   :: the_g_data ! Result
    character (LEN=*), intent(IN)           :: description ! e.g., 'dao'
    integer, intent(out)                    :: returnStatus ! E.g., FILENOTFOUND
    character (LEN=*), intent(IN)           :: GeoDimList ! Comma-delimited dim names
    character (LEN=*), intent(IN)           :: fieldNames ! Name of gridded field
    real(rgr), optional, intent(IN)         :: missingValue
    character (LEN=*), optional, intent(IN) :: date ! offset
    logical, optional, intent(IN)           :: sumDelp ! sum the DELP to make PL
    logical, optional, intent(IN)           :: deferReading ! don't read yet
    character(len=*), optional, intent(out) :: litDescription
    logical, optional, intent(IN)           :: verbose

    ! Local Variables
    character ( len=len(fieldNames)) :: fieldName   ! In case we supply two
    logical, parameter :: DEEBUG = .false.
    integer            :: Me = -1                   ! String index for trace
    logical            :: myDefer
    logical            :: myVerbose
    ! Executable code
    call trace_begin ( me, "ReadGriddedData_MLSFile", lcf_where, &
      & cond=toggle(gen) .and. levels(gen) > 0 )
    myDefer = .false.
    if ( present(deferReading) ) myDefer = deferReading
    myVerbose = deebug
    if ( present(verbose) ) myVerbose = verbose .or. deebug
    
    LIT_DESCRIPTION = lowercase(description)
    ! The following once allowed us to describe both geos5-classic and merra
    ! reanalysis with the same griddedorigin field='geos5'
    ! However, we have since seen the light: that trick was a later
    ! source of confusion--we will eventually simply remove litDescription
    if ( LIT_DESCRIPTION == 'geos5' .and. GEOS5MAYBEMERRATOO ) &
      & LIT_DESCRIPTION = GEOS5orMERRA( GriddedFile )
    if ( myVerbose ) &
      & call output( 'Reading ' // trim(LIT_DESCRIPTION) // ' data', advance='yes' )

    call nullifyGriddedData ( the_g_data ) ! for Sun's still useless compiler
    the_g_data%empty = .true.
    the_g_data%QuantityName = '(none)'
    the_g_data%description  = '(none)'
    the_g_data%units        = '(none)'
    the_g_data%equivalentLatitude = .false.
    returnStatus = mls_hdf_version(GriddedFile%Name)
    if ( returnStatus == FILENOTFOUND ) then
      call SetupNewGriddedData ( the_g_data, empty=.true. )
      go to 9
    else
      returnStatus = 0
    end if
    ! Are we deferring the actual reading until later?
    ! If so, we'll store the dummy args for use later
    if ( myDefer ) then
      the_g_data%sourceFileName     = GriddedFile%Name
      the_g_data%description        = description
      the_g_data%dimList            = GeoDimList
      the_g_data%fieldNames         = fieldNames
      the_g_data%verticalCoordinate = v_type
      the_g_data%empty              = .true.
      the_g_data%equivalentLatitude = .false.
      go to 9
    end if
    
    ! Nope, let's read it!
    call getStringElement(fieldNames, fieldName, 1, COUNTEMPTY)
    ! According to the kinds of gridded data files we can read
    select case ( trim(LIT_DESCRIPTION) )
    case ('geos5_7')
      ! Check that hdf version is OK
      if ( mls_hdf_version( GriddedFile%Name ) /= HDFVERSION_5 ) then
        the_g_data%empty = .true.
        returnStatus = FILENOTFOUND
        call MLSMessage( MLSMSG_Warning, ModuleName, &
          & 'Not an hdf5 file so could not be geos 5.7.x' )
        if ( present(litDescription) ) litDescription = LIT_DESCRIPTION
        go to 9
      end if
      call Read_geos5_7( GriddedFile, lcf_where, v_type, &
        & the_g_data, GeoDimList, fieldName, date, sumDelp )
      if ( .not. myVerbose ) then
        ! Do nothing -- but some compilers may complain
      else if ( the_g_data%empty .or. .not. associated(the_g_data%field) ) then
          call output( 'File appears not to be ' // trim(LIT_DESCRIPTION), advance='yes' )
      else
        call output( '(Returned from read_geos5_7)', advance='yes' )
        call output( 'Quantity Name   ' // trim(the_g_data%QuantityName), advance='yes' )
        call output( 'Description     ' // trim(the_g_data%description), advance='yes' )
        call output( 'Units           ' // trim(the_g_data%units), advance='yes' )
        call outputNamedValue( 'Vertical Coord, type, type(P)  ', &
          & (/ the_g_data%verticalCoordinate, v_type, v_is_pressure /) )
        call outputNamedValue( 'max val  ', mlsmax( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'min val  ', mlsmin( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'meanval  ', mlsmean( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue('associated(the_g_data%field)', associated(the_g_data%field))
        call outputNamedValue('NumStringElements', NumStringElements(fieldNames, COUNTEMPTY))
      end if
    case ('geos5')
      call Read_geos5_or_merra( GriddedFile, lcf_where, v_type, &
        & the_g_data, GeoDimList, fieldName, date )
      if ( .not. myVerbose ) then
        ! Do nothing -- but some compilers may complain
      else if ( the_g_data%empty .or. .not. associated(the_g_data%field) ) then
          call output( 'File appears not to be ' // trim(LIT_DESCRIPTION), advance='yes' )
      else
        call output( '(Returned from read_geos5)', advance='yes' )
        call output( 'Quantity Name   ' // trim(the_g_data%QuantityName), advance='yes' )
        call output( 'Description     ' // trim(the_g_data%description), advance='yes' )
        call output( 'Units           ' // trim(the_g_data%units), advance='yes' )
        call outputNamedValue( 'Vertical Coord, type, type(P)  ', &
          & (/ the_g_data%verticalCoordinate, v_type, v_is_pressure /) )
        call outputNamedValue( 'max val  ', mlsmax( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'min val  ', mlsmin( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'meanval  ', mlsmean( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue('associated(the_g_data%field)', associated(the_g_data%field))
        call outputNamedValue('NumStringElements', NumStringElements(fieldNames, COUNTEMPTY))
      end if
    case ('dao', 'gmao')
      call Read_dao(GriddedFile, lcf_where, v_type, &
        & the_g_data, GeoDimList, fieldName)
      if ( .not. myVerbose ) then
        ! Do nothing -- but some compilers may complain
      else if ( the_g_data%empty .or. .not. associated(the_g_data%field) ) then
          call output( 'File appears not to be ' // trim(LIT_DESCRIPTION), advance='yes' )
      else
        call output( '(Returned from read_dao)', advance='yes' )
        call output( 'Quantity Name   ' // trim(the_g_data%QuantityName), advance='yes' )
        call output( 'Description     ' // trim(the_g_data%description), advance='yes' )
        call output( 'Units           ' // trim(the_g_data%units), advance='yes' )
        call outputNamedValue( 'Vertical Coord, type, type(P)  ', &
          & (/ the_g_data%verticalCoordinate, v_type, v_is_pressure /) )
        call outputNamedValue( 'max val  ', mlsmax( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'min val  ', mlsmin( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'meanval  ', mlsmean( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue('associated(the_g_data%field)', associated(the_g_data%field))
        call outputNamedValue('NumStringElements', NumStringElements(fieldNames, COUNTEMPTY))
      end if
    case ('merra')
      ! In case we were forced to acknowledge that out GMAO files
      ! are MERRA and not GEOS5
      if ( NumStringElements(fieldNames, COUNTEMPTY) > 1 ) &
        & call getStringElement(fieldNames, fieldName, 2, COUNTEMPTY)
      call Read_geos5_or_merra( GriddedFile, lcf_where, v_type, &
        & the_g_data, GeoDimList, fieldName, date, sumDelp )
      if ( .not. myVerbose ) then
        ! Do nothing -- but some compilers may complain
      else if ( the_g_data%empty .or. .not. associated(the_g_data%field) ) then
          call output( 'File appears not to be ' // trim(LIT_DESCRIPTION), advance='yes' )
      else
        call output( '(Returned from read merra)', advance='yes' )
        call output( 'Quantity Name   ' // trim(the_g_data%QuantityName), advance='yes' )
        call output( 'Description     ' // trim(the_g_data%description), advance='yes' )
        call output( 'Units           ' // trim(the_g_data%units), advance='yes' )
        call outputNamedValue( 'Vertical Coord, type, type(P)  ', &
          & (/ the_g_data%verticalCoordinate, v_type, v_is_pressure /) )
        call outputNamedValue( 'max val  ', mlsmax( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'min val  ', mlsmin( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'meanval  ', mlsmean( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue('associated(the_g_data%field)', associated(the_g_data%field))
        call outputNamedValue('NumStringElements', NumStringElements(fieldNames, COUNTEMPTY))
      end if
    case ('ncep')
      ! These are ncep global assimilation model data
      ! in hdfeos format
      call Read_ncep_gdas(GriddedFile, lcf_where, v_type, &
        & the_g_data, GeoDimList, fieldName, missingValue)
      if ( .not. myVerbose ) then
        ! Do nothing -- but some compilers may complain
      else if ( the_g_data%empty .or. .not. associated(the_g_data%field) ) then
          call output( 'File appears not to be ' // trim(LIT_DESCRIPTION), advance='yes' )
      else
        call output( '(Returned from read_ncep_gdas)', advance='yes' )
        call output( 'Quantity Name   ' // trim(the_g_data%QuantityName), advance='yes' )
        call output( 'Description     ' // trim(the_g_data%description), advance='yes' )
        call output( 'Units           ' // trim(the_g_data%units), advance='yes' )
        call outputNamedValue( 'Vertical Coord, type, type(P)  ', &
          & (/ the_g_data%verticalCoordinate, v_type, v_is_pressure /) )
        call outputNamedValue( 'max val  ', mlsmax( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'min val  ', mlsmin( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'meanval  ', mlsmean( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue('associated(the_g_data%field)', associated(the_g_data%field))
        call outputNamedValue('NumStringElements', NumStringElements(fieldNames, COUNTEMPTY))
      end if
    case ('strat')
      ! These are ncep stratospheric analysis combined data
      ! in hdfeos5 format
      call Read_ncep_strat(GriddedFile, lcf_where, v_type, &
        & the_g_data, GeoDimList, fieldName)
      if ( .not. myVerbose ) then
        ! Do nothing -- but some compilers may complain
      else if ( the_g_data%empty .or. .not. associated(the_g_data%field) ) then
          call output( 'File appears not to be ' // trim(LIT_DESCRIPTION), advance='yes' )
      else
        call output( '(Returned from read_ncep_strat)', advance='yes' )
        call output( 'Quantity Name   ' // trim(the_g_data%QuantityName), advance='yes' )
        call output( 'Description     ' // trim(the_g_data%description), advance='yes' )
        call output( 'Units           ' // trim(the_g_data%units), advance='yes' )
        call outputNamedValue( 'Vertical Coord, type, type(P)  ', &
          & (/ the_g_data%verticalCoordinate, v_type, v_is_pressure /) )
        call outputNamedValue( 'max val  ', mlsmax( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'min val  ', mlsmin( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'meanval  ', mlsmean( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue('associated(the_g_data%field)', associated(the_g_data%field))
        call outputNamedValue('NumStringElements', NumStringElements(fieldNames, COUNTEMPTY))
      end if
    case default
      call announce_error(lcf_where, 'ReadGriddedData_MLSFile called with unknown' &
        & // ' description: ' // trim(LIT_DESCRIPTION))
    end select
    if ( present(litDescription) ) litDescription = LIT_DESCRIPTION

9   continue
    call trace_end ( "rprocessOneAprioriFile", &
      & cond=toggle(gen) .and. levels(gen) > 0 )

  end subroutine ReadGriddedData_MLSFile

  subroutine ReadGriddedData_Name( FileName, lcf_where, description, v_type, &
    & the_g_data, returnStatus, &
    & GeoDimList, fieldNames, missingValue, &
    & date, sumDelp, deferReading, litDescription, verbose )
    ! Arguments
    character (LEN=*), intent(IN)           :: FileName
    integer, intent(IN)                     :: lcf_where    ! node of the lcf that provoked me         
    integer, intent(IN)                     :: v_type       ! vertical coordinate; an 'enumerated' type
    type( GriddedData_T )                   :: the_g_data ! Result
    character (LEN=*), intent(IN)           :: description ! e.g., 'dao'
    integer, intent(out)                    :: returnStatus ! E.g., FILENOTFOUND
    character (LEN=*), intent(IN)           :: GeoDimList ! Comma-delimited dim names
    character (LEN=*), intent(IN)           :: fieldNames ! Name of gridded field
    real(rgr), optional, intent(IN)         :: missingValue
    character (LEN=*), optional, intent(IN) :: date ! offset
    logical, optional, intent(IN)           :: sumDelp ! sum the DELP to make PL
    logical, optional, intent(IN)           :: deferReading ! don't read yet
    character(len=*), optional, intent(out) :: litDescription
    logical, optional, intent(IN)           :: verbose
  end subroutine ReadGriddedData_Name

  ! ----------------------------------------------- Read_geos5_7
  subroutine Read_geos5_7( GEOS5File, lcf_where, v_type, &
    & the_g_data, GeoDimList, fieldName, date, sumDelp )

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use DATES_MODULE, only: UTC2TAI93S
    use HDF5, only: HSIZE_T
    use MLSHDF5, only: GETALLHDF5DSNAMES, &
      & GETHDF5DSRANK, GETHDF5DSDIMS, LOADFROMHDF5DS, &
      & READHDF5ATTRIBUTE
    use MLSSTRINGLISTS, only: SWITCHDETAIL
    use TOGGLES, only: SWITCHES

    ! This routine reads a gmao geos5_7 file, named something like
    ! DAS.ops.asm.avg3_3d_Nv.GEOS571.20110930_0300.V01.nc4 (pressure with

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
    character (LEN=*), optional, intent(IN) :: GeoDimList ! Comma-delimited dim names
    character (LEN=*), optional, intent(IN) :: fieldName ! Name of gridded field
    character (LEN=*), optional, intent(IN) :: date ! date (offset)
    logical, optional, intent(IN)           :: sumDelp ! sum the DELP to make PL
    ! Local variables
    character (len=NAMELEN) :: actual_field_name
    character(len=19) :: datestring ! will be in the form of yyyy-MM-ddTHH:MM:ss
    logical :: DEEBUG
    integer(kind=hsize_t) :: dim1(1), dim4(4)
    integer :: error, rank
    character(len=256) :: errormsg
    integer :: i1, i2, i3            ! looping indexes
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
    if ( verbose ) call dump(mysdList, 'DS names')

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

    !if (.not. ReadHDF5Attribute(geos5file%fileid%f_id, 'time', &
    !    & 'time_increment', timeinc, error=errormsg)) then
    !    call announce_error (lcf_where, errormsg // ' in file ' &
    !    & // geos5file%name)
    !end if

    !if (timeinc <= 0) then
    !    call announce_error (lcf_where, "Invalid 'time_increment' value in " // geos5file%name)
    !end if

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

    !hour = timeinc / 10000
    !minute = mod(timeinc, 10000) / 100
    !second = mod(timeinc, 100)
    !if (hour > 24 .or. minute > 60 .or. second > 60) then
    !    call announce_error (lcf_where, "Invalid 'time_increment' value in " // geos5file%name)
    !end if

    !the_g_data%dateends(1) = the_g_data%datestarts(1) + hour * 3600 + minute * 60 + second
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
    case ( 'pl' )
        the_units = 'Pa'
    case ( 't' )
        the_units = 'K'
    case default
        call announce_error(lcf_where, "Unexpected data: " // actual_field_name)
    end select
    
    the_g_data%units = the_units

    ! read the field, field should be either PL or T
    call GetHDF5DSRank (geos5file%fileid%f_id, capitalize(actual_field_name), rank)
    if (rank /= 4) then
        call announce_error (lcf_where, &
        capitalize(actual_field_name) // " must be a 3-dim array: " // geos5file%name)
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
    call test_allocate ( error, moduleName, 'the_g_data%field', (/1,1,1,1,1,1/), &
      & (/ idim4(3), idim4(2), idim4(1), 1, 1, 1 /), bytes(the_g_data%field) )
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
  end subroutine Read_geos5_7

  ! ----------------------------------------------- Read_geos5_or_merra
  subroutine Read_geos5_or_merra( GEOS5File, lcf_where, v_type, &
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
    character (LEN=*), optional, intent(IN) :: GeoDimList ! Comma-delimited dim names
    character (LEN=*), optional, intent(IN) :: fieldName ! Name of gridded field
    character (LEN=*), optional, intent(IN) :: date ! date (offset)
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
    integer, parameter             :: i_longitude=1
    integer, parameter             :: i_latitude=i_longitude+1
    integer, parameter             :: i_vertical=i_latitude+1
    integer, parameter             :: i_time=i_vertical+1
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
          & 'list size too big in Read_geos5_or_merra ' // trim(GEOS5File%Name), MLSFile=GEOS5File )
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
        & 'all_the_fields', ModuleName // 'Read_geos5_or_merra' )
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
        & ModuleName // 'Read_geos5_or_merra' )
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
        & ModuleName // 'Read_geos5_or_merra' )
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
        & ModuleName // 'Read_geos5_or_merra' )
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
        & ModuleName // 'Read_geos5_or_merra' )
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
  end subroutine Read_geos5_or_merra

  ! ----------------------------------------------- Read_dao
  subroutine Read_dao(DAOFile, lcf_where, v_type, &
    & the_g_data, GeoDimList, fieldName)

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test

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
    character (LEN=*), optional, intent(IN) :: GeoDimList ! Comma-delimited dim names
    character (LEN=*), optional, intent(IN) :: fieldName ! Name of gridded field

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
    integer, parameter             :: i_longitude=1
    integer, parameter             :: i_latitude=i_longitude+1
    integer, parameter             :: i_vertical=i_latitude+1
    integer, parameter             :: i_time=i_vertical+1
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
          & 'list size too big in Read_dao ' // trim(DAOFile%Name), MLSFile=DAOFile )
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
  end subroutine Read_dao

  ! ----------------------------------------------- Read_ncep_gdas
  subroutine Read_ncep_gdas ( NCEPFile, lcf_where, v_type, &
    & the_g_data, GeoDimList, gridName, missingValue )

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test

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
    character (LEN=*), optional, intent(IN) :: GeoDimList ! E.g., 'X.Y,..'
    character (LEN=*), optional, intent(IN) :: gridName ! Name of grid
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
    integer, parameter             :: i_longitude=1
    integer, parameter             :: i_latitude=i_longitude+1
    integer, parameter             :: i_vertical=i_latitude+1
    integer, parameter             :: i_time=i_vertical+1
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
      call announce_error(lcf_where, "Read_ncep_gdas unable to use projection")
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
  end subroutine Read_ncep_gdas

  ! ----------------------------------------------- Read_ncep_strat
  subroutine Read_ncep_strat(NCEPFile, lcf_where, v_type, &
    & the_g_data, GeoDimList, fieldName)
    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use HDFEOS5, only: HE5_HDFE_NENTDIM, HE5F_ACC_RDONLY, &
      & HE5_GDOPEN, HE5_GDATTACH, HE5_GDDETACH, HE5_GDCLOSE, &
      & HE5_GDNENTRIES, HE5_GDINQGRID, HE5_GDINQDIMS, HE5_GDINQFLDS, &
      & HE5_GDFLDINFO, HE5_GDGRIDINFO
    use MLSHDFEOS, only: HSIZES
    use hdf5, only: SIZE_T

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
    character (LEN=*), optional, intent(IN) :: GeoDimList ! E.g., 'X.Y,..'
    character (LEN=*), optional, intent(IN) :: fieldName  ! Name of grid quant.

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
    integer, parameter             :: i_longitude=1
    integer, parameter             :: i_latitude=i_longitude+1
    integer, parameter             :: i_vertical=i_latitude+1
    integer, parameter             :: i_time=i_vertical+1
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
  end subroutine Read_ncep_strat

  ! --------------------------------------------- ReadGloriaFile -------
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

  ! --------------------------------------------------  Read_Climatology
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
    character (LEN=FileNameLen)            :: fname   ! Physical file name
    character (LEN=FileNameLen)            :: path    ! Physical path

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
      processCli = GetPCFromRef(fname, mlspcf_l2clim_start, mlspcf_l2clim_end, &
        & .true., ErrType, version, debugOption=debug)
      if(ErrType /= 0) then
        call announce_error ( root, &
          &"Climatology file name " // trim(fname) // " unmatched in PCF", &
          & 'error number: ', extra_number=ErrType)
        go to 9
      end if
      ErrType = Pgs_io_gen_openF ( processCli, PGSd_IO_Gen_RSeqFrm, 0, &
        cliUnit, version )
    else
      fname = climFile%name
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

  ! ----------------------------------------------- WriteGriddedData
  subroutine WriteGriddedData( GriddedFile, lcf_where, description, v_type, &
    & the_g_data, returnStatus, &
    & GeoDimList, fieldNames, missingValue, date, sumDelp )

    ! This routine Writes a Gridded Data file, based on an array of filled data
    ! structures and the  appropriate for the description
    ! which may be one of 'geos5', 'gmao', 'dao', 'merra', 'ncep', 'strat'
    ! (But see shortc. and lim. below)

    ! FileName and the_g_data are required args
    ! GeoDimList should be the Dimensions' short names
    ! as a comma-delimited character string in the order:
    ! longitude, latitude, vertical level, time

    ! fieldNames should be the data fields' short names
    ! as a comma-delimited character string in the order
    ! of the elements of the_g_data array

    ! fieldName should name the rank 3 or higher object
    ! like temperature
    
    ! ---- Shortcomings and limitations: ----
    ! We have only coded it for description='merra' so far
    ! Any other file description will generate an error
    ! Obviously we should only code what we plan to use
    ! Would it be useful to be able to write a geos5_7 file, perhaps
    ! degrading its resolution thereby?
    
    ! Arguments
    type(MLSFile_T)                         :: GriddedFile
    integer, intent(IN)                     :: lcf_where    ! node of the lcf that provoked me         
    integer, intent(IN)                     :: v_type       ! vertical coordinate; an 'enumerated' type
    type( GriddedData_T ), dimension(:)     :: the_g_data   ! Result
    character (LEN=*), intent(IN)           :: description ! e.g., 'dao'
    integer, intent(out)                    :: returnStatus ! E.g., FILENOTFOUND
    character (LEN=*), intent(IN)           :: GeoDimList ! Comma-delimited dim names
    character (LEN=*), intent(IN)           :: fieldNames ! Names of gridded field
    real(rgr), optional, intent(IN)         :: missingValue
    character (LEN=*), optional, intent(IN) :: date ! offset
    logical, optional, intent(IN)           :: sumDelp ! sum the DELP to make PL

    ! Local Variables
    integer :: fieldIndex
    character ( len=len(fieldNames)) :: fieldName   ! In case we supply two
    character ( len=NameLen) :: my_description   ! In case mixed case
    logical, parameter :: DEEBUG = .false.
    ! Executable code
    returnStatus = 0
    my_description = lowercase(description)
    if ( my_description == 'geos5' ) &
      & my_description = GEOS5orMERRA( GriddedFile )
    if ( DEEBUG ) print *, 'Writeing ' // trim(my_description) // ' data'

    do fieldIndex=1, size(the_g_data)
      call getStringElement(fieldNames, fieldName, fieldIndex, COUNTEMPTY)
      ! According to the kinds of gridded data files we can Write
      select case ( trim(my_description) )
      case ('geos5')
        call announce_error(lcf_where, 'WriteGriddedData called with illegal' &
          & // ' description: ' // trim(my_description))
      case ('dao', 'gmao')
        call announce_error(lcf_where, 'WriteGriddedData called with illegal' &
          & // ' description: ' // trim(my_description))
      case ('merra')
        call Write_merra( (fieldIndex==1), GriddedFile, lcf_where, &
          & the_g_data(fieldIndex), GeoDimList, fieldName )
      case ('ncep')
        call announce_error(lcf_where, 'WriteGriddedData called with illegal' &
          & // ' description: ' // trim(my_description))
      case ('strat')
        call announce_error(lcf_where, 'WriteGriddedData called with illegal' &
          & // ' description: ' // trim(my_description))
      case default
        call announce_error(lcf_where, 'WriteGriddedData called with unknown' &
          & // ' description: ' // trim(my_description))
      end select
    end do
  end subroutine WriteGriddedData

  ! ----------------------------------------------- Write_merra
  subroutine Write_merra( createGrid, GEOS5File, lcf_where, &
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
    character (LEN=*), optional, intent(IN) :: GeoDimList ! Comma-delimited dim names
    character (LEN=*), optional, intent(IN) :: fieldName ! Name of gridded field

    ! Local Variables
    integer :: file_id, gd_id
    logical,  parameter       :: CASESENSITIVE = .false.
    integer, parameter :: GRIDORDER=1   ! What order grid written to file
    integer, parameter :: MAXNAMELENGTH=NameLen         ! Max length of grid name
    integer, dimension(NENTRIESMAX) :: dims

    integer :: start(4), stride(4), edge(4)
    integer :: status
    integer, parameter             :: i_longitude=1
    integer, parameter             :: i_latitude=i_longitude+1
    integer, parameter             :: i_vertical=i_latitude+1
    integer, parameter             :: i_time=i_vertical+1
    character(len=*), parameter    :: GridName = "EOSGRID"
    integer, parameter :: gctp_geo = 0
    real(r4), parameter :: FILLVALUE = 1.e15
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

  end subroutine Write_merra

  ! ----- utility procedures ----
  function GEOS5orMERRA( File ) result( fType )
    ! Attempt to identify file as
    ! (1) GEOS5, or
    ! (2) MERRA
    ! based on size of dims(4) of field "T"
    ! Arguments
    type(MLSFile_T)                :: File
    character(len=8)               :: fType
    ! Local Variables
    integer :: file_id, gd_id
    integer :: inq_success
    character (len=MAXLISTLENGTH) :: gridlist
    character (len=MAXLISTLENGTH), dimension(1) :: dimlists
    integer :: ngrids
    integer                        :: our_rank, numberType
    integer :: strbufsize
    logical,  parameter       :: CASESENSITIVE = .false.
    integer, parameter :: GRIDORDER=1   ! What order grid written to file
    integer, parameter :: MAXNAMELENGTH=NameLen         ! Max length of grid name
    character (len=MAXNAMELENGTH) :: gridname
    integer, dimension(NENTRIESMAX) :: dims
    integer :: status
    logical, parameter :: DEEBUG = .false.
    ! Executable
    fType = 'geos5'
    gridlist = ''
    inq_success = gdinqgrid(File%Name, gridlist, strbufsize)
    if (inq_success < 0) then
      call announce_error(0, "Could not inquire gridlist "// &
        & trim(File%Name))
    else if ( strbufsize > MAXLISTLENGTH ) then
       CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'list size too big in Read_merra ' // trim(File%Name), MLSFile=File )
    else if ( strbufsize < MAXLISTLENGTH .and. strbufsize > 0 ) then
      gridlist = gridlist(1:strbufsize) // ' '
    end if
    file_id = gdopen(File%Name, DFACC_RDONLY)

    if (file_id < 0) then
      call announce_error(0, "Could not open "// File%Name)
    end if
    ! Find grid name corresponding to the GRIDORDER'th one
    ngrids = NumStringElements(gridlist, COUNTEMPTY)

    if(ngrids <= 0) then
      call announce_error(0, "NumStringElements of gridlist <= 0")
    else if(ngrids /= inq_success) then
      call announce_error(0, "NumStringElements of gridlist /= inq_success")
    else if(ngrids < GRIDORDER) then
      call announce_error(0, "NumStringElements of gridlist < GRIDORDER")
    end if

    call GetStringElement(gridlist, gridname, GRIDORDER, COUNTEMPTY)

    gd_id = gdattach(file_id, gridname)
    if (gd_id < 0) then
      call announce_error(0, "Could not attach "//trim(gridname))
    end if

    ! Now find the rank of our field
    inq_success = gdfldinfo(gd_id, 'T', our_rank, dims, &
      & numbertype, dimlists(1))

    if ( DEEBUG ) then
      print *, 'our_rank ', our_rank
      print *, 'dims ', dims(1:our_rank)
    end if

    if ( dims(4) > 1 ) then
      fType = 'merra'
    else
      fType = 'geos5'
    end if
    status = gddetach(gd_id)
    status = gdclose(file_id)
  end function GEOS5orMERRA
  
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
    else if(default_output_by_toolkit) then
      just_print_it = .false.
    else
      just_print_it = .true.
    end if

    if(.not. just_print_it) then
      error = max(error,1)
      call output ( '***** At ' )

      if(lcf_where > 0) then
        call print_source ( where(lcf_where) )
      else
        call output ( '(no lcf node available)' )
      end if

      call output ( ': ' )
      call output ( "The " );
      if(lcf_where > 0) then
        call dump_tree_node ( lcf_where, 0 )
      else
        call output ( '(no lcf tree available)' )
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

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module ncep_dao

! $Log$
! Revision 2.79  2014/09/06 00:05:02  pwagner
! Stores geost_7 lats in the_g_data after reading them
!
! Revision 2.78  2014/09/05 00:27:11  vsnyder
! More complete and accurate allocate/deallocate size tracking.  Add some
! tracing.
!
! Revision 2.77  2014/06/04 18:29:45  pwagner
! Account for memory usage accurately; allow deferReading
!
! Revision 2.76  2014/04/02 23:02:52  pwagner
! Removed redundant open_ and close_MLSFile
!
! Revision 2.75  2014/01/09 00:24:29  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.74  2013/09/24 23:27:14  vsnyder
! Use Get_Where or Print_Source to start error messages
!
! Revision 2.73  2012/07/11 20:01:06  pwagner
! Dont use DEEBUG until after its defined
!
! Revision 2.72  2012/07/02 20:17:42  pwagner
! Protect against bounds error when debugging with NAG
!
! Revision 2.71  2012/05/15 16:53:06  pwagner
! Allocate first, then give values
!
! Revision 2.70  2012/05/08 17:46:05  pwagner
! Fixed bug that caused confusion when geos5 not ffound
!
! Revision 2.69  2012/03/06 19:33:11  pwagner
! Remove more unused things
!
! Revision 2.68  2012/03/01 20:00:35  pwagner
! When verbose, note when file does not fit description
!
! Revision 2.67  2011/11/30 21:28:23  pwagner
! Converted most debugging prints to verbose ones; more robust
!
! Revision 2.66  2011/08/30 22:23:46  pwagner
! Should not try to read as geos5_7 if not hdf5
!
! Revision 2.65  2011/08/03 22:50:03  pwagner
! Repaired syntax of write; removed unused variables; updated toc
!
! Revision 2.64  2011/08/03 15:36:23  honghanh
! Move one line of code around for better readability
!
! Revision 2.63  2011/08/02 17:26:37  honghanh
! Remove code that read time_increment attribute
!
! Revision 2.62  2011/08/02 17:23:59  honghanh
! Make datestarts and dateends the same.
!
! Revision 2.61  2011/08/02 16:48:51  honghanh
! Implement readgeos5_7
!
! Revision 2.60  2011/06/29 21:38:23  pwagner
! Uses output api
!
! Revision 2.59  2011/05/05 23:09:22  pwagner
! verbose when requested via optional arg
!
! Revision 2.58  2011/04/27 17:36:12  pwagner
! Consistently sets grid%empty when not found or wrong description
!
! Revision 2.57  2011/04/20 16:56:16  pwagner
! First steps toward being able to read next GMAO format
!
! Revision 2.56  2010/03/31 18:12:10  pwagner
! Added WriteGriddedData
!
! Revision 2.55  2010/02/04 23:08:00  vsnyder
! Remove USE or declaration for unused names
!
! Revision 2.54  2009/11/04 23:16:49  pwagner
! Fixed bugs in recognizing and reading merra files
!
! Revision 2.53  2009/09/29 23:35:09  pwagner
! Changes needed by 64-bit build
!
! Revision 2.52  2009/08/17 16:52:19  pwagner
! From MERRA file may read DELP field and perform sum
!
! Revision 2.51  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.50  2009/06/16 17:25:58  pwagner
! Removed most unused stuff
!
! Revision 2.49  2008/12/02 23:11:13  pwagner
! mls_io_gen_[openF,closeF] functions now private; use MLSFile_T interfaces instead
!
! Revision 2.48  2008/09/17 23:21:19  pwagner
! Allow date string in gridded data to offset gmao background files
!
! Revision 2.47  2008/01/07 21:38:12  pwagner
! Replace DEFAULTUNDEFINEDVALUE with user-settable undefinedValue
!
! Revision 2.46  2007/06/21 00:49:52  vsnyder
! Remove tabs, which are not part of the Fortran standard
!
! Revision 2.45  2006/11/01 20:30:11  pwagner
! More unused debugging prints
!
! Revision 2.44  2006/06/13 22:11:45  pwagner
! Correctly sets units, heightsUnits
!
! Revision 2.43  2006/05/19 19:55:07  pwagner
! Corrected a misspelling Lahey missed
!
! Revision 2.42  2006/05/18 18:39:32  cvuu
! Add subroutine Read_geos5
!
! Revision 2.41  2006/05/12 21:25:56  pwagner
! verticalCoordinate now correctly set by read_dao
!
! Revision 2.40  2006/01/25 00:58:44  pwagner
! Cleared out some commented-out stuff
!
! Revision 2.39  2005/06/22 17:24:59  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.38  2005/06/14 20:37:06  pwagner
! Interfaces changed to accept MLSFile_T args
!
! Revision 2.37  2004/08/04 23:19:01  pwagner
! Much moved from MLSStrings to MLSStringLists
!
! Revision 2.36  2004/08/03 17:59:35  pwagner
! Gets DEFAULTUNDEFINEDVALUE from MLSCommon
!
! Revision 2.35  2004/06/23 17:11:50  pwagner
! read_climatology now returns status, e.g. FILENOTFOUND
!
! Revision 2.34  2004/01/23 01:12:27  pwagner
! Some care taken in handling ...inq.. functions
!
! Revision 2.33  2003/10/06 13:30:40  cvuu
! Modified to handle reading the ncep data for origin=strat
!
! Revision 2.32  2003/06/03 20:42:25  pwagner
! Can read strat (ncep) files; untested
!
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
