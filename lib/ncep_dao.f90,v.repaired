! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module ncep_dao ! Collections of subroutines to handle DataType GriddedData_T

  ! This module's name is in the running for worst-chosen of all
  ! A better choice might be Read_Gridded_Data_m
  ! use, Intrinsic :: ISO_C_Binding, only: C_Intptr_T, C_Loc
  use GriddedData, only: GriddedData_T, V_Is_Pressure, RGR, &
    & SetupNewGriddedData, NullifyGriddedData
  use HDFeos, only: GDOpen, GDAttach, GDDetach, GDClose, GDFldinfo, &
    & GDInqGrid
  use HDF, only: Dfacc_Rdonly
  use HighOutput, only: OutputNamedValue
  use Intrinsic, only: L_HDFeos
  use MLSCommon, only: LineLen, NameLen, &
    & MLSFile_T
  use MLSFiles, only: FileNotFound, HDFversion_5, InitializeMLSFile, &
    & MLS_HDF_Version
  use MLSMessagemodule, only: MLSMsg_Error, MLSMsg_Warning, &
    & MLSMessage
  use MLSStrings, only: Lowercase
  use MLSStringlists, only: GetStringElement, NumStringElements
  use Output_M, only: Output
  use ReadGriddedUtils, only: Lit_Description, Announce_Error, &
    & Read_Dao, Read_Geos5_Or_Merra, Read_Merra_2, &
    & Read_Geos5_7, Read_Ncep_Gdas, Read_Ncep_Strat, ReadGloriaFile, &
    & Write_Merra

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
! ReadGriddedData      read meteorology file in one of supported desciptions
! WriteGriddedData     write meteorology file in one of supported desciptions
!
! The supported descriptions are
! geos5_7         Geos5.7 +  gmao files in netcdf4/hdf5 format
! geos5           Geos5x.x   gmao files in hdfeos2/hdf4 format
! dao or gmao     Geos4x.x   gmao files in hdfeos2/hdf4 format
! merra           Geos5x.x   gmao reanalysis files in hdfeos2/hdf4 format
! merra_2         merra2     merra2 files in netcdf4/hdf5 format
! ncep            GDAS       ncep files in hdfeos2/hdf4 format
! strat           STRAT      ncep combined in hdfeos5/hdf5 format
! === (end of toc) ===

  public:: ReadGriddedData, ReadGloriaFile
  public:: WriteGriddedData

  ! First we'll define some global parameters and data types.
  logical, parameter :: COUNTEMPTY=.true.
!  character (len=*), parameter :: DEFAULTDAODIMLIST = 'XDim,YDim,Height,Time'
!  character (len=*), parameter :: DEFAULTDAOFIELDNAME = 'TMPU'
!  character (len=*), parameter :: DEFAULTGEOS5FIELDNAME = 'T'
!  character (len=*), parameter :: DEFAULTNCEPDIMLIST = 'YDim,XDim'
!  character (len=*), parameter :: DEFAULTNCEPGRIDNAME = 'TMP_3'
!  character (len=*), parameter :: DEFAULTNCEPSTRATFIELDNAME = 'Temperature'
!  character (len=*), parameter :: GEO_FIELD1 = 'Latitude'
!  character (len=*), parameter :: GEO_FIELD2 = 'Longitude'
!  character (len=*), parameter :: GEO_FIELD3 = 'Height'
!  character (len=*), parameter :: GEO_FIELD4 = 'Time'
!  logical, parameter :: GEOS5MAYBEMERRATOO = .false. ! Causes confusion if true
! 
!  character (len=*), parameter :: lit_dao = 'dao'
!  character (len=*), parameter :: lit_ncep = 'ncep'
!  character (len=*), parameter :: lit_strat = 'strat'
!  character (len=*), parameter :: lit_clim = 'clim'
!  character (len=*), parameter :: lit_geos5 = 'geos5'
! 
! 
  integer, parameter :: MAXLISTLENGTH=Linelen ! Max length list of grid names
  integer, parameter :: NENTRIESMAX=200 ! Max num of entries
! integer, parameter ::          MAXDS = 1024 ! 500
! integer, parameter ::          MAXSDNAMESBUFSIZE = MAXDS*NAMELEN

  interface ReadGRIDDEDDATA
    module procedure ReadGriddedData_MLSFile
    module procedure ReadGriddedData_Name
  end interface

contains

  ! ----------------------------------------------- ReadGriddedData_MLSFile
  ! This family of routines reads a Gridded Data file, returning a filled data
  ! structure and the  appropriate for the description
  ! which may be one of 
  ! 'geos5_7', 'geos5', 'gmao', 'dao', 'merra', 'merra_2', 'ncep', 'strat'

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

    use MLSStats1, only: MLSMin, MLSMax, MLSMean
    use Toggles, only: Gen, Levels, Toggle
    use Trace_M, only: Trace_Begin, Trace_End

    ! Arguments
    type(MLSFile_T)                         :: GriddedFile
    integer, intent(in)                     :: lcf_where    ! node of the lcf that provoked me         
    integer, intent(in)                     :: v_type       ! vertical coordinate; an 'enumerated' type
    type( GriddedData_T )                   :: the_g_data ! Result
    character (len=*), intent(in)           :: description ! e.g., 'dao'
    integer, intent(out)                    :: returnStatus ! E.g., FILENOTFOUND
    character (len=*), intent(in)           :: GeoDimList ! Comma-delimited dim names
    character (len=*), intent(in)           :: fieldNames ! Name of gridded field
    real(rgr), optional, intent(in)         :: missingValue
    character (len=*), optional, intent(in) :: date ! offset
    logical, optional, intent(in)           :: sumDelp ! sum the DELP to make PL
    logical, optional, intent(in)           :: deferReading ! don't read yet
    character(len=*), optional, intent(out) :: litDescription
    logical, optional, intent(in)           :: verbose

    ! Local Variables
    character ( len=len(fieldNames)) :: fieldName   ! In case we supply two
    logical, parameter :: DEEBUG = .false.
    integer            :: Me = -1                   ! String index for trace
    logical            :: myDefer
    logical            :: myVerbose
    ! Executable code
    call trace_begin ( me, "ReadGriddedData_MLSFile", lcf_where, &
      & cond=toggle(gen) .and. levels(gen) > 0 )
    if ( present(litDescription) ) litDescription = 'none'
    myDefer = .false.
    if ( present(deferReading) ) myDefer = deferReading
    myVerbose = deebug
    if ( present(verbose) ) myVerbose = verbose .or. deebug
    
    LIT_DESCRIPTION = lowercase(description)
    if ( myVerbose ) &
      & call output( 'Reading ' // trim(LIT_DESCRIPTION) // ' data', advance='yes' )

    call nullifyGriddedData ( the_g_data ) ! for Sun's still useless compiler
    the_g_data%empty = .true.
    the_g_data%QuantityName = '(none)'
    the_g_data%description  = '(none)'
    the_g_data%units        = '(none)'
    the_g_data%equivalentLatitude = .false.
    the_g_data%noYear = .false.
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
    case ( 'geos5_7', 'merra_2' ) ! Actually 5.7 or merra_2; netCDF4 format
      ! Check that hdf version is OK
      if ( mls_hdf_version( GriddedFile%Name ) /= HDFVERSION_5 ) then
        the_g_data%empty = .true.
        returnStatus = FILENOTFOUND
        call MLSMessage( MLSMSG_Warning, ModuleName, &
          & 'Not an hdf5 file so could not be geos 5.7.x or merra2' )
        if ( present(litDescription) ) litDescription = LIT_DESCRIPTION
        go to 9
      end if
      if ( trim(LIT_DESCRIPTION) == 'geos5_7' ) then
        call Read_geos5_7( GriddedFile, lcf_where, v_type, &
          & the_g_data, GeoDimList, fieldName, date, sumDelp )
      else
        call Read_merra_2( GriddedFile, lcf_where, v_type, &
          & the_g_data, GeoDimList, fieldName, date, sumDelp )
      endif
      if ( .not. myVerbose ) then
        ! Do nothing -- but some compilers may complain
      else if ( the_g_data%empty .or. .not. associated(the_g_data%field) ) then
          call output( 'File appears not to be ' // trim(LIT_DESCRIPTION), advance='yes' )
      else
        call output( '(Returned from read_geos5_7 or merra_2)', advance='yes' )
        call output( 'Quantity Name   ' // trim(the_g_data%QuantityName), advance='yes' )
        call output( 'Description     ' // trim(the_g_data%description), advance='yes' )
        call output( 'Units           ' // trim(the_g_data%units), advance='yes' )
        call outputNamedValue( 'Vertical Coord, type, type(P)  ', &
          & (/ the_g_data%verticalCoordinate, v_type, v_is_pressure /) )
        call outputNamedValue( 'max val  ', MLSMax( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'min val  ', MLSMin( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'meanval  ', MLSMean( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue('associated(the_g_data%field)', associated(the_g_data%field))
        call outputNamedValue('NumStringElements', NumStringElements(fieldNames, COUNTEMPTY))
      end if
    case ( 'geos5' )
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
        call outputNamedValue( 'max val  ', MLSMax( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'min val  ', MLSMin( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'meanval  ', MLSMean( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
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
        call outputNamedValue( 'max val  ', MLSMax( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'min val  ', MLSMin( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'meanval  ', MLSMean( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
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
        call outputNamedValue( 'max val  ', MLSMax( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'min val  ', MLSMin( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'meanval  ', MLSMean( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
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
        call outputNamedValue( 'max val  ', MLSMax( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'min val  ', MLSMin( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'meanval  ', MLSMean( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
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
        call outputNamedValue( 'max val  ', MLSMax( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'min val  ', MLSMin( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue( 'meanval  ', MLSMean( the_g_data%field(:,:,:,1,1,1), the_g_data%missingValue ) )
        call outputNamedValue('associated(the_g_data%field)', associated(the_g_data%field))
        call outputNamedValue('NumStringElements', NumStringElements(fieldNames, COUNTEMPTY))
      end if
    case default
      call announce_error( lcf_where, 'ReadGriddedData_MLSFile called with unknown' &
        & // ' description: ' // trim(LIT_DESCRIPTION) )
    end select
    if ( present(litDescription) ) litDescription = LIT_DESCRIPTION

9   continue
    call trace_end ( "ReadGriddedData_MLSFile", &
      & cond=toggle(gen) .and. levels(gen) > 0 )

  end subroutine ReadGriddedData_MLSFile

  subroutine ReadGriddedData_Name( FileName, lcf_where, description, v_type, &
    & the_g_data, returnStatus, &
    & GeoDimList, fieldNames, missingValue, &
    & date, sumDelp, deferReading, litDescription, verbose )
    ! Arguments
    character (len=*), intent(in)           :: FileName
    integer, intent(in)                     :: lcf_where    ! node of the lcf that provoked me         
    integer, intent(in)                     :: v_type       ! vertical coordinate; an 'enumerated' type
    type( GriddedData_T )                   :: the_g_data ! Result
    character (len=*), intent(in)           :: description ! e.g., 'dao'
    integer, intent(out)                    :: returnStatus ! E.g., FILENOTFOUND
    character (len=*), intent(in)           :: GeoDimList ! Comma-delimited dim names
    character (len=*), intent(in)           :: fieldNames ! Name of gridded field
    real(rgr), optional, intent(in)         :: missingValue
    character (len=*), optional, intent(in) :: date ! offset
    logical, optional, intent(in)           :: sumDelp ! sum the DELP to make PL
    logical, optional, intent(in)           :: deferReading ! don't read yet
    character(len=*), optional, intent(out) :: litDescription
    logical, optional, intent(in)           :: verbose
    ! Internal variables
    type(MLSFile_T)                         :: GriddedFile

    ! Executable
    returnStatus = InitializeMLSFile( GriddedFile, content = 'gridded', &
      & name=FileName, &
      & type=l_hdfeos, access=DFACC_RDONLY )
    call ReadGriddedData_MLSFile( GriddedFile, lcf_where, description, v_type, &
      & the_g_data, returnStatus, &
      & GeoDimList, fieldNames, missingValue, &
      & date, sumDelp, deferReading, litDescription, verbose )

  end subroutine ReadGriddedData_Name

  ! -------------------------------------------  WriteGriddedData  -----
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
    integer, intent(in)                     :: lcf_where    ! node of the lcf that provoked me         
    integer, intent(in)                     :: v_type       ! vertical coordinate; an 'enumerated' type
    type( GriddedData_T ), dimension(:)     :: the_g_data   ! Result
    character (len=*), intent(in)           :: description ! e.g., 'dao'
    integer, intent(out)                    :: returnStatus ! E.g., FILENOTFOUND
    character (len=*), intent(in)           :: GeoDimList ! Comma-delimited dim names
    character (len=*), intent(in)           :: fieldNames ! Names of gridded field
    real(rgr), optional, intent(in)         :: missingValue
    character (len=*), optional, intent(in) :: date ! offset
    logical, optional, intent(in)           :: sumDelp ! sum the DELP to make PL

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
    if ( DEEBUG ) print *, 'Writing ' // trim(my_description) // ' data'

    do fieldIndex=1, size(the_g_data)
      call getStringElement(fieldNames, fieldName, fieldIndex, COUNTEMPTY)
      ! According to the kinds of gridded data files we can Write
      select case ( trim(my_description) )
      case ('geos5')
        call announce_error( lcf_where, 'WriteGriddedData called with illegal' &
          & // ' description: ' // trim(my_description) )
      case ('dao', 'gmao')
        call announce_error( lcf_where, 'WriteGriddedData called with illegal' &
          & // ' description: ' // trim(my_description) )
      case ('merra')
        call Write_merra( (fieldIndex==1), GriddedFile, lcf_where, &
          & the_g_data(fieldIndex), GeoDimList, fieldName )
      case ('ncep')
        call announce_error( lcf_where, 'WriteGriddedData called with illegal' &
          & // ' description: ' // trim(my_description) )
      case ('strat')
        call announce_error( lcf_where, 'WriteGriddedData called with illegal' &
          & // ' description: ' // trim(my_description) )
      case default
        call announce_error( lcf_where, 'WriteGriddedData called with unknown' &
          & // ' description: ' // trim(my_description) )
      end select
    end do
  end subroutine WriteGriddedData

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
      call announce_error( 0, "Could not inquire gridlist "// &
        & trim(File%Name) )
    else if ( strbufsize > MAXLISTLENGTH ) then
       CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'list size too big in Read_merra ' // trim(File%Name), MLSFile=File )
    else if ( strbufsize < MAXLISTLENGTH .and. strbufsize > 0 ) then
      gridlist = gridlist(1:strbufsize) // ' '
    end if
    file_id = gdopen(File%Name, DFACC_RDONLY)

    if (file_id < 0) then
      call announce_error( 0, "Could not open "// File%Name )
    end if
    ! Find grid name corresponding to the GRIDORDER'th one
    ngrids = NumStringElements(gridlist, COUNTEMPTY)

    if(ngrids <= 0) then
      call announce_error( 0, "NumStringElements of gridlist <= 0" )
    else if(ngrids /= inq_success) then 
      call announce_error( 0, "NumStringElements of gridlist /= inq_success" )
    else if(ngrids < GRIDORDER) then
      call announce_error( 0, "NumStringElements of gridlist < GRIDORDER" )
    end if

    call GetStringElement(gridlist, gridname, GRIDORDER, COUNTEMPTY)

    gd_id = gdattach(file_id, gridname)
    if (gd_id < 0) then
      call announce_error( 0, "Could not attach "//trim(gridname) )
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
! Revision 2.87  2017/10/17 23:39:48  pwagner
! Removed unused stuff
!
! Revision 2.86  2017/07/10 18:22:46  pwagner
! Correct arg to trace_end; preset noYear
!
! Revision 2.85  2017/03/07 21:17:35  pwagner
! The specific procedures for each origin type moved to a Utils module
!
! Revision 2.84  2017/01/13 01:29:47  pwagner
! Pre-initialize litDescription = 'none' in ReadGriddedData
!
! Revision 2.83  2016/09/22 22:22:47  pwagner
! Change message to acknowledge version may be later than geos5_7
!
! Revision 2.82  2016/07/28 01:42:27  vsnyder
! Refactoring dump and diff
!
! Revision 2.81  2015/03/28 00:56:14  vsnyder
! Stuff to trace allocate/deallocate addresses -- mostly commented out
! because NAG build 1017 doesn't yet allow arrays as arguments to C_LOC.
!
! Revision 2.80  2015/02/06 01:10:36  pwagner
! Added body missing from ReadGriddedData with filename api
!
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
