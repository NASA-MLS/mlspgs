

! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ncep_dao ! Collections of subroutines to handle TYPE GriddedData_T
  !=============================================================================


  use Allocate_Deallocate, only: Allocate_test, Deallocate_test
  use GriddedData, only: GriddedData_T, v_is_pressure, v_is_altitude, &
    & v_is_gph, v_is_theta
  use HDFEOS, only: HDFE_NENTDIM, HDFE_NENTDFLD, &
    & gdopen, gdattach, gddetach, gdclose, gdfldinfo, &
    & gdinqgrid, gdnentries, gdinqdims, gdinqflds, gddiminfo
  use Hdf, only: SUCCEED, DFACC_RDONLY
  use l3ascii, only: l3ascii_read_field
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: R8, LineLen, NameLen, FileNameLen
  use MLSFiles, only: GetPCFromRef, mls_io_gen_closeF, mls_io_gen_openF, &
    &                split_path_name
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Allocate, &
    & MLSMSG_Deallocate, MLSMSG_Warning
  use MLSStrings, only: GetStringElement, NumStringElements, Capitalize, &
    & GetIntHashElement, LowerCase
  use OUTPUT_M, only: OUTPUT
  use SDPToolkit, only: PGS_S_SUCCESS, PGS_PC_GETREFERENCE, &
    & PGS_IO_GEN_CLOSEF, PGS_IO_GEN_OPENF, PGSD_IO_GEN_RSEQFRM, &
    & UseSDPToolkit
  use Toggles, only: Gen, Levels, Switches, Toggle
  use Trace_M, only: Trace_begin, Trace_end
  use TREE, only: DUMP_TREE_NODE, SOURCE_REF

  implicit none
  public

  private :: Id,ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  character(LEN=130) :: id = & 
    "$Id$"
  character(LEN=*), parameter :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------



  public::SetupNewGridTemplate, DestroyGridTemplateContents, &
    &   AddGridTemplateToDatabase, DestroyGridTemplateDatabase, &
    & Dump_Gridded_Database, source_file_already_read
  public::OBTAIN_CLIM, READ_CLIMATOLOGY, OBTAIN_DAO, Obtain_NCEP
  public::ReadGriddedData
  private::announce_error
  private::DEFAULTFIELDNAME, GEO_FIELD1, GEO_FIELD2, GEO_FIELD3, GEO_FIELD4
  private::lit_dao, lit_ncep, lit_clim
  integer, private :: ERROR

  ! First we'll define some global parameters and data types.

  character (len=*), parameter :: DEFAULTFIELDNAME = 'TMPU'
  character (len=*), parameter :: GEO_FIELD1 = 'Latitude'
  character (len=*), parameter :: GEO_FIELD2 = 'Longitude'
  character (len=*), parameter :: GEO_FIELD3 = 'Height'
  character (len=*), parameter :: GEO_FIELD4 = 'Time'

  character (len=*), parameter :: lit_dao = 'dao'
  character (len=*), parameter :: lit_ncep = 'ncep'
  character (len=*), parameter :: lit_clim = 'clim'

  ! This datatype stores a single gridded atmospheric quantity.  For example
  ! temperature, if an uncertainty field is also required, this is stored in a
  ! separate quantity.


  ! --------------------------------------------------------------------------

contains

  ! Now we have some subroutines to deal with these quantitites

  ! This first routine sets up a new quantity template according to the user
  ! input.  This may be based on a previously supplied template (with possible
  ! modifications), or created from scratch.

  subroutine SetupNewGridTemplate(qty, source, noHeights, noLats, noLons, noLsts, noSzas, noDates)

    ! Dummy arguments
    type (GriddedData_T), intent(OUT) :: qty ! Result

    type (GriddedData_T), optional, intent(IN) :: source ! Template

    integer, optional, intent(IN) :: noHeights, noLats, noLons, noLsts, noSzas, noDates

    ! Local variables
    integer :: status           ! Status from allocates etc.

    ! Executable code

    ! First, if we have a template setup according to that
    if (present(source)) then
      qty%noHeights=source%noHeights
      qty%noLats=source%noLats
      qty%noLons=source%noLons
      qty%noLsts=source%noLsts
      qty%noSzas=source%noSzas
      qty%noDates=source%noDates


    else ! We have no template, setup a very bare quantity
      qty%noHeights=1
      qty%noLats=1
      qty%noLons=1
      qty%noLsts=1
      qty%noSzas=1
      qty%noDates=1

    endif

    ! Now, see if the user asked for modifications to this
    if (present(noHeights)) qty%noHeights=noHeights
    if (present(noLats)) qty%noLats=noLats
    if (present(noLons)) qty%noLons=noLons
    if (present(noLsts)) qty%noLsts=noLsts
    if (present(noSzas)) qty%noSzas=noSzas
    if (present(noDates)) qty%noDates=noDates

    ! First the vertical/horizontal coordinates
    call Allocate_test ( qty%heights, qty%noHeights, "qty%heights", ModuleName )
    call Allocate_test ( qty%lats, qty%noLats, "qty%lats", ModuleName )
    call Allocate_test ( qty%lons, qty%noLons, "qty%lons", ModuleName )
    call Allocate_test ( qty%lsts, qty%noLsts, "qty%lsts", ModuleName )
    call Allocate_test ( qty%szas, qty%noSzas, "qty%szas", ModuleName )

    ! Now the temporal coordinates
    call Allocate_test ( qty%dateStarts, qty%noDates, "qty%dateStarts", ModuleName )
    call Allocate_test ( qty%dateEnds, qty%noDates, "qty%dateEnds", ModuleName )

    ! Now the data itself
    allocate(qty%field(qty%noHeights, qty%noLats, qty%noLons,  &
      qty%noLsts, qty%noSzas, qty%noDates), STAT=status)

    if (status /= 0) call announce_error(0,  &
      & MLSMSG_Allocate//"field")

  end subroutine SetupNewGridTemplate

  ! --------------------------------------------------------------------------

  ! This subroutine destroys a quantity template

  subroutine DestroyGridTemplateContents(qty)

    ! Dummy argument
    type (GriddedData_T), intent(INOUT) :: qty

    ! Local variables
    integer :: STATUS

    ! Executable code
    call Deallocate_test ( qty%heights, "qty%heights", ModuleName )
    call Deallocate_test ( qty%lats, "qty%lats", ModuleName )
    call Deallocate_test ( qty%lons, "qty%lons", ModuleName )
    call Deallocate_test ( qty%lsts, "qty%lsts", ModuleName )
    call Deallocate_test ( qty%szas, "qty%szas", ModuleName )

    ! Now the temporal coordinates
    call Deallocate_test ( qty%dateStarts, "qty%dateStarts", ModuleName )
    call Deallocate_test ( qty%dateEnds, "qty%dateEnds", ModuleName )

    ! Now the data itself
    deallocate(qty%field, STAT=status)

    if (status /= 0) call announce_error(0,  &
      & MLSMSG_Allocate//"field")

  end subroutine DestroyGridTemplateContents

  ! --------------------------------------------------------------------------

  ! This subroutine adds a quantity template to a database, or creates the
  ! database if it doesn't yet exist

  !  SUBROUTINE AddGridTemplateToDatabase(database,qty)
  integer function AddGridTemplateToDatabase(database,item)

    ! Dummy arguments
    type (GriddedData_T), dimension(:), pointer :: database
    type (GriddedData_T), intent(IN) :: item

    ! Local variables
    type (GriddedData_T), dimension(:), pointer :: tempDatabase

    ! Executable code

    include "addItemToDatabase.f9h"
    AddGridTemplateToDatabase = newSize

  end function AddGridTemplateToDatabase

  ! --------------------------------------------------------------------------

  ! This subroutine destroys a quantity template database

  subroutine DestroyGridTemplateDatabase(database)

    ! Dummy argument
    type (GriddedData_T), dimension(:), pointer :: database

    ! Local variables
    integer :: qtyIndex, status

    if ( toggle(gen) ) call trace_begin ( "DestroyGridTemplateDatabase" )

    if (associated(database)) then
      do qtyIndex=1,size(database)
        call DestroyGridTemplateContents(database(qtyIndex))
      enddo
      deallocate(database, stat=status)
      if (status /= 0) call announce_error(0,  &
        & MLSMSG_DeAllocate//"database")
    endif
    if ( toggle(gen) ) then
      call trace_end ( "DestroyGridTemplateDatabase" )
    end if
  end subroutine DestroyGridTemplateDatabase


  !----------------- Beginning of Paul's code ------------------

  !---------------------------- ReadGriddedData ---------------------
  subroutine ReadGriddedData(FileName, lcf_where, description, v_type, &
    & the_g_data, GeoDimList, fieldName)
    !------------------------------------------------------------------------

    ! This routine reads a Gridded Data file, returning a filled data structure and the !
    ! appropriate for 'ncep' or 'dao'

    ! FileName and the_g_data are required args
    ! GeoDimList, if present, should be the Dimensions' short names
    ! as a comma-delimited character string in the order:
    ! longitude, latitude, vertical level, time

    ! fieldName, if present, should be the rank 3 or higher object
    ! like temperature

    ! Arguments

    character (LEN=*), intent(IN) :: FileName ! Name of the file containing the grid(s)
    integer, intent(IN) :: lcf_where			! node of the lcf that provoked me
    integer, intent(IN) :: v_type       ! vertical coordinate; an 'enumerated' type
    type( GriddedData_T ), intent(OUT) :: the_g_data ! Result
    character (LEN=*), intent(IN) :: description ! e.g., 'dao'
    character (LEN=*), optional, intent(IN) :: GeoDimList ! Comma-delimited dim names
    character (LEN=*), optional, intent(IN) :: fieldName ! Name of gridded field

    ! Local Variables

    integer :: edges(4)
    integer :: file_id, gd_id
    integer :: inq_success
    integer :: i
    integer :: nentries, ngrids, ndims, nfields
    integer :: strbufsize
    !  character (len=80) :: msg, mnemonic
    !  integer :: status

    logical,  parameter       :: CASESENSITIVE = .false.
    integer, parameter :: GRIDORDER=1   ! What order grid written to file
    integer, parameter :: MAXLISTLENGTH=Linelen ! Max length list of grid names
    integer, parameter :: NENTRIESMAX=20 ! Max num of entries
    character (len=MAXLISTLENGTH) :: gridlist
    character (len=MAXLISTLENGTH) :: dimlist, actual_dim_list
    character (len=MAXLISTLENGTH), dimension(1) :: dimlists
    character (len=MAXLISTLENGTH) :: fieldlist
    integer, parameter :: MAXNAMELENGTH=NameLen		! Max length of grid name
    character (len=MAXNAMELENGTH) :: gridname, actual_field_name, the_dim
    integer, dimension(NENTRIESMAX) :: dims, rank, numberTypes, start, stride
    integer                        :: our_rank, size, numberType
    !                                  These start out initialized to one
    integer                        :: nlon=1, nlat=1, nlev=1, ntime=1
    integer, parameter             :: i_longitude=1
    integer, parameter             :: i_latitude=i_longitude+1
    integer, parameter             :: i_vertical=i_latitude+1
    integer, parameter             :: i_time=i_vertical+1
    ! External functions
    !  integer, external :: gdopen, gdattach, gdrdfld, gddetach, gdclose
    !  integer, external :: gdinqgrid, gdnentries, gdinqdims, gdinqflds
    integer, external :: GDRDFLD
    logical, parameter :: COUNTEMPTY=.true.
    logical            :: descrpt_is_legal
    logical            :: descrpt_is_misplcd

    ! - - - begin - - -

    ! Check if description is legal
    descrpt_is_legal = (lowercase(description(:len(lit_dao))) == lit_dao) &
      & .or. &
      & (lowercase(description(:len(lit_ncep))) == lit_ncep) &
      & .or. &
      & (lowercase(description(:len(lit_clim))) == lit_clim)

    descrpt_is_misplcd = lowercase(description(:len(lit_clim))) == lit_clim

    if(descrpt_is_misplcd) then
      call announce_error(lcf_where, 'READGriddedData called with climatology' &
        & // ' description')
      return
    elseif(.not. descrpt_is_legal) then
      call announce_error(lcf_where, 'READGriddedData called with unknown' &
        & // ' description: ' // description)
      return
    endif

    error = 0
    file_id = gdopen(FileName, DFACC_RDONLY)

    if (file_id < 0) then
      call announce_error(lcf_where, "Could not open "// FileName)
    end if

    ! Find list of grid names on this file
    inq_success = gdinqgrid(FileName, gridlist, strbufsize)
    if (inq_success < 0) then
      call announce_error(lcf_where, "Could not inquire gridlist "// FileName)
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
      actual_field_name=DEFAULTFIELDNAME
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

    !-----------------------------
  end subroutine ReadGriddedData
  !-----------------------------

  ! ------------------------------------------------  OBTAIN_CLIM  -----
  !=====================================================================
  subroutine OBTAIN_CLIM ( aprioriData, root, &
    & mlspcf_l2clim_start, mlspcf_l2clim_end )
    !=====================================================================

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
    character (LEN=256) :: msg, mnemonic
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
          returnStatus = AddGridTemplateToDatabase(aprioriData, qty)

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
    !============================
  end subroutine OBTAIN_CLIM
  !============================


  ! --------------------------------------------------  READ_CLIMATOLOGY  -----
  subroutine READ_CLIMATOLOGY ( input_fname, root, aprioriData, &
    & mlspcf_l2clim_start, mlspcf_l2clim_end, echo_data, dump_data )
    ! --------------------------------------------------
    ! Brief description of program
    ! This subroutine reads a l3ascii file and returns
    ! the data_array to the caller

    ! Arguments

    character*(*), intent(in) :: input_fname			! Physical file name
    type (GriddedData_T), dimension(:), pointer :: aprioriData 
    integer, intent(in) :: ROOT        ! Root of the L2CF abstract syntax tree
    integer, optional, intent(IN) :: mlspcf_l2clim_start, mlspcf_l2clim_end
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
    logical :: dump
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
      dump = dump_data
    else
      dump = DUMP_GRIDDED_QUANTITIES
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
          & error_number=ErrType)
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

          if(dump) then
            call Dump_Gridded_Data(gddata, root)
          endif

          ErrType = AddGridTemplateToDatabase(aprioriData, gddata)

          nullify (gddata%heights)
          nullify (gddata%lats)
          nullify (gddata%lons)
          nullify (gddata%lsts)
          nullify (gddata%szas)
          nullify (gddata%dateStarts)
          nullify (gddata%dateEnds)
          nullify (gddata%field)

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
          &"Error closing " // fname, error_number=ErrType)
      endif

    else

      call announce_error (ROOT, &
        &"Error opening " // fname, error_number=ErrType)
    endif

  end subroutine READ_CLIMATOLOGY

  ! -------------------------------------------------  OBTAIN_DAO  -----
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
    character (len=256) :: mnemonic, msg
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
        returnStatus = AddGridTemplateToDatabase(aprioriData, qty)

        !        call read_dao ( DAOphysicalFilename, vname, data_array )
        if(returnStatus > 0) then
          call ReadGriddedData ( DAOphysicalFilename, root, &
            & 'dao', v_is_pressure, aprioriData(returnStatus) )
        endif

      end if

    end do ! DAOFileHandle = mlspcf_l2_dao_start, mlspcf_l2_dao_end

    !===========================
  end subroutine Obtain_DAO
  !===========================

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
    !    character (len=80) :: MSG, MNEMONIC
    integer :: NCEPFileHandle, NCEP_Version
    character (len=132) :: NCEPphysicalFilename
    type (GriddedData_T):: QTY
    integer :: RETURNSTATUS
    !    character (len=80) :: VNAME

    !   allocate ( data_array(XDIM, YDIM, ZDIM), stat=returnStatus )

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

        returnStatus = AddGridTemplateToDatabase(aprioriData, qty)

        !        call read_ncep ( NCEPphysicalFilename, data_array )
        if(returnStatus > 0) then
          call ReadGriddedData ( NCEPphysicalFilename, root, &
            & 'ncep', v_is_pressure, aprioriData(returnStatus) )
        endif

      end if

    end do ! NCEPFileHandle = mlspcf_l2_ncep_start, mlspcf_l2_ncep_end
    !===========================
  end subroutine Obtain_NCEP
  !===========================

  ! --------------------------------  Dump_Gridded_Database  -----
  subroutine Dump_Gridded_Database(GriddedData, root)
    use Dump_0, only: Dump

    ! Imitating what dump_pointing_grid_database does, but for gridded data
    ! which may come from climatology, ncep, dao

    type (GriddedData_T), dimension(:), pointer :: GriddedData 

    integer, intent(in) :: ROOT        ! Root of the L2CF abstract syntax tree

    ! Local Variables
    logical, parameter :: MAYDUMPFIELDVALUES = .false.
    integer            :: i

    if ( .not. associated(GriddedData)) then
      call announce_error(ROOT, 'Gridded database still null')
      return
    endif

    call output ( 'database: a priori grids: SIZE = ' )
    call output ( size(GriddedData), advance='yes' )
    do i = 1, size(GriddedData)

      call output ( 'item number ' )
      call output ( i, advance='yes' )

      call Dump_Gridded_Data(GriddedData(i), root)
    end do ! i
  end subroutine Dump_Gridded_Database

  ! --------------------------------  Dump_Gridded_Data  -----
  subroutine Dump_Gridded_Data(GriddedData, root)
    use Dump_0, only: Dump

    ! Imitating what dump_pointing_grid_database does, but for gridded data
    ! which may come from climatology, ncep, dao

    type (GriddedData_T) :: GriddedData 

    integer, intent(in) :: ROOT        ! Root of the L2CF abstract syntax tree

    ! Local Variables
    logical, parameter :: MAYDUMPFIELDVALUES = .false.

    call output('quantity name ' // GriddedData%quantityName, advance='yes')
    call output('description ' // GriddedData%description, advance='yes')
    call output('units ' // GriddedData%units, advance='yes')

    call output ( ' ************ Geometry ********** ' ,advance='yes')

    call output ( ' Vertical coordinate = ' )
    call output ( GriddedData%verticalCoordinate, advance='yes' )
    call output ( ' No. of heights = ' )
    call output ( GriddedData%noHeights, advance='yes' )
    call dump ( GriddedData%heights, &
      & '    Heights =' )

    call output ( ' Equivalent latitude = ' )
    call output ( GriddedData%equivalentLatitude, advance='yes' )
    call output ( ' No. of latitudes = ' )
    call output ( GriddedData%noLats, advance='yes' )
    call dump ( GriddedData%lats, &
      & '    latitudes =' )

    call output ( ' No. of longitudes = ' )
    call output ( GriddedData%noLons, advance='yes' )
    call dump ( GriddedData%lons, &
      & '    longitudes =' )

    call output ( ' No. of local times = ' )
    call output ( GriddedData%noLsts, advance='yes' )
    call dump ( GriddedData%lsts, &
      & '    local times =' )

    call output ( ' No. of solar zenith angles = ' )
    call output ( GriddedData%noSzas, advance='yes' )
    call dump ( GriddedData%szas, &
      & '    solar zenith angles =' )

    call output ( ' No. of dates = ' )
    call output ( GriddedData%noDates, advance='yes' )
    call dump ( GriddedData%dateStarts, &
      & '    starting dates =' )
    call dump ( GriddedData%dateEnds, &
      & '    ending dates =' )

    if(MAYDUMPFIELDVALUES) then
      call output ( ' ************ tabulated field values ********** ' ,advance='yes')

      ! No dump for 6-dimensional double arrays yet, anyway
      !     call dump ( GriddedData%field, &
      !      & '    gridded field values =' )
    endif

  end subroutine Dump_Gridded_Data

  ! --------------------------------  source_file_already_read  -----
  function source_file_already_read(GriddedDataBase, source_file, field_name &
    & ) result (already)
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

    print*,'Database, filenames:',GriddedDatabase%sourceFilename
    print*,'Database, fieldNames:',GriddedDatabase%quantityName
    print*,'Source file:',source_file
    print*,'field_name:',field_name

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

  ! ------------------------------------------------  announce_error  -----
  subroutine announce_error ( lcf_where, full_message, use_toolkit, &
    & error_number )

    ! Arguments

    integer, intent(in)    :: lcf_where
    character(LEN=*), intent(in)    :: full_message
    logical, intent(in), optional :: use_toolkit
    integer, intent(in), optional    :: error_number
    ! Local
    !  character (len=80) :: msg, mnemonic
    !  integer :: status
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
      !    CALL Pgs_smf_getMsg(status, mnemonic, msg)
      !    CALL MLSMessage (level, ModuleName, &
      !              &trim(full_message)//" "//mnemonic//" "//msg)
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
      if(present(error_number)) then
        call output('error number ', advance='no')
        call output(error_number, places=9, advance='yes')
      endif
    else
      print*, '***Error in module ', ModuleName
      print*, trim(full_message)
      if(present(error_number)) then
        print*, 'error number ', error_number
      endif
    endif

    !===========================
  end subroutine announce_error
  !===========================

  !=============================================================================
end module ncep_dao
!=============================================================================

!
! $Log$
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
