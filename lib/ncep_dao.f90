

! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE ncep_dao ! Collections of subroutines to handle TYPE GriddedData_T
!=============================================================================


  use GriddedData, only: GriddedData_T, v_is_pressure, v_is_altitude, &
  & v_is_gph, v_is_theta
  use HDFEOS, only: HDFE_NENTDIM, HDFE_NENTDFLD, &
  & gdopen, gdattach, gddetach, gdclose, gdfldinfo, &
  & gdinqgrid, gdnentries, gdinqdims, gdinqflds, gddiminfo
  use Hdf, only: SUCCEED, DFACC_RDONLY
  use l3ascii, only: l3ascii_read_field
  use LEXER_CORE, only: PRINT_SOURCE
  USE MLSCommon, only: R8, LineLen, NameLen
  USE MLSFiles, only: GetPCFromRef, mls_io_gen_closeF, mls_io_gen_openF
  USE MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Allocate, &
  & MLSMSG_Deallocate, MLSMSG_Warning
  USE MLSStrings, only: GetStringElement, NumStringElements, Capitalize, &
  & GetIntHashElement, LowerCase
  use OUTPUT_M, only: OUTPUT
  USE SDPToolkit, only: PGS_S_SUCCESS, PGS_PC_GETREFERENCE, &
  & PGS_IO_GEN_CLOSEF, PGS_IO_GEN_OPENF, PGSD_IO_GEN_RSEQFRM
  use TREE, only: DUMP_TREE_NODE, SOURCE_REF

  IMPLICIT NONE
  PUBLIC

  PRIVATE :: Id,ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = & 
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

 

  public::SetupNewGridTemplate, DestroyGridTemplateContents, &
	&   AddGridTemplateToDatabase, DestroyGridTemplateDatabase, &
	& Dump_Gridded_Database
  public::OBTAIN_CLIM, READ_CLIMATOLOGY, OBTAIN_DAO, Obtain_NCEP
  public::ReadGriddedData
  private::announce_error
  private::DEFAULTFIELDNAME, GEO_FIELD1, GEO_FIELD2, GEO_FIELD3, GEO_FIELD4
  private::lit_dao, lit_ncep, lit_clim
  integer, private :: ERROR

  ! First we'll define some global parameters and data types.

   CHARACTER (len=*), PARAMETER :: DEFAULTFIELDNAME = 'TMPU'
   CHARACTER (len=*), PARAMETER :: GEO_FIELD1 = 'Latitude'
   CHARACTER (len=*), PARAMETER :: GEO_FIELD2 = 'Longitude'
   CHARACTER (len=*), PARAMETER :: GEO_FIELD3 = 'Height'
   CHARACTER (len=*), PARAMETER :: GEO_FIELD4 = 'Time'

   CHARACTER (len=*), PARAMETER :: lit_dao = 'dao'
   CHARACTER (len=*), PARAMETER :: lit_ncep = 'ncep'
   CHARACTER (len=*), PARAMETER :: lit_clim = 'clim'

! This datatype stores a single gridded atmospheric quantity.  For example
! temperature, if an uncertainty field is also required, this is stored in a
! separate quantity.


  ! --------------------------------------------------------------------------

  CONTAINS

  ! Now we have some subroutines to deal with these quantitites

  ! This first routine sets up a new quantity template according to the user
  ! input.  This may be based on a previously supplied template (with possible
  ! modifications), or created from scratch.

  SUBROUTINE SetupNewGridTemplate(qty, source, noHeights, noLats, noLons, noLsts, noSzas, noDates)

    ! Dummy arguments
    TYPE (GriddedData_T), INTENT(OUT) :: qty ! Result

    TYPE (GriddedData_T), OPTIONAL, INTENT(IN) :: source ! Template

    INTEGER, OPTIONAL, INTENT(IN) :: noHeights, noLats, noLons, noLsts, noSzas, noDates

    ! Local variables
    INTEGER :: status           ! Status from allocates etc.

    ! Executable code

    ! First, if we have a template setup according to that
    IF (PRESENT(source)) THEN
       qty%noHeights=source%noHeights
       qty%noLats=source%noLats
       qty%noLons=source%noLons
       qty%noLsts=source%noLsts
       qty%noSzas=source%noSzas
       qty%noDates=source%noDates

      
    ELSE ! We have no template, setup a very bare quantity
       qty%noHeights=1
       qty%noLats=1
       qty%noLons=1
       qty%noLsts=1
       qty%noSzas=1
       qty%noDates=1

    ENDIF

    ! Now, see if the user asked for modifications to this
    IF (PRESENT(noHeights)) qty%noHeights=noHeights
    IF (PRESENT(noLats)) qty%noLats=noLats
    IF (PRESENT(noLons)) qty%noLons=noLons
    IF (PRESENT(noLsts)) qty%noLsts=noLsts
    IF (PRESENT(noSzas)) qty%noSzas=noSzas
    IF (PRESENT(noDates)) qty%noDates=noDates
    ! First the vertical coordinates

    ALLOCATE (qty%heights(qty%noHeights),STAT=status)
    IF (status /= 0) call announce_error(0,  &
         & MLSMSG_Allocate//"heights")

    ! Now the geolocation coordinates
    ALLOCATE (qty%lats(qty%noLats),STAT=status)
    IF (status /= 0) call announce_error(0,  &
         & MLSMSG_Allocate//"lats")

    ALLOCATE (qty%lons(qty%noLons),STAT=status)
    IF (status /= 0) call announce_error(0,  &
         & MLSMSG_Allocate//"lons")

    ALLOCATE (qty%lsts(qty%noLsts),STAT=status)
    IF (status /= 0) call announce_error(0,  &
         & MLSMSG_Allocate//"lsts")

    ALLOCATE (qty%szas(qty%noSzas),STAT=status)
    IF (status /= 0) call announce_error(0,  &
         & MLSMSG_Allocate//"szas")

    !Now the temporal coordinates
    ALLOCATE (qty%DateStarts(qty%noDates),STAT=status)
    IF (status /= 0) call announce_error(0,  &
         & MLSMSG_Allocate//"DateStarts")

    ALLOCATE (qty%DateEnds(qty%noDates),STAT=status)
    IF (status /= 0) call announce_error(0,  &
         & MLSMSG_Allocate//"DateEnds")

    !Now the data itself
    ALLOCATE(qty%field(qty%noHeights, qty%noLats, qty%noLons,  &
             qty%noLsts, qty%noSzas, qty%noDates), STAT=status)

    IF (status /= 0) call announce_error(0,  &
         & MLSMSG_Allocate//"field")


  END SUBROUTINE SetupNewGridTemplate

  ! --------------------------------------------------------------------------

  ! This subroutine destroys a quantity template

  SUBROUTINE DestroyGridTemplateContents(qty)

    ! Dummy argument
    TYPE (GriddedData_T), INTENT(INOUT) :: qty
    ! Local variables
    INTEGER status

    ! Executable code

    DEALLOCATE (qty%heights, STAT=status)

    IF (status /= 0) call announce_error(0,  &
         & MLSMSG_DeAllocate//"heights")

    DEALLOCATE (qty%lats, STAT=status)

    IF (status /= 0) call announce_error(0,  &
         & MLSMSG_DeAllocate//"lats")

    DEALLOCATE (qty%lons, STAT=status)

    IF (status /= 0) call announce_error(0,  &
         & MLSMSG_DeAllocate//"lons")

    DEALLOCATE (qty%lsts, STAT=status)

    IF (status /= 0) call announce_error(0,  &
         & MLSMSG_DeAllocate//"lsts")

    DEALLOCATE (qty%szas, STAT=status)

    IF (status /= 0) call announce_error(0,  &
         & MLSMSG_DeAllocate//"szas")

    DEALLOCATE (qty%DateStarts, STAT=status)

    IF (status /= 0) call announce_error(0,  &
         & MLSMSG_DeAllocate//"DateStarts")

    DEALLOCATE (qty%DateEnds, STAT=status)

    IF (status /= 0) call announce_error(0,  &
         & MLSMSG_DeAllocate//"DateEnds")

    DEALLOCATE (qty%field, STAT=status)

    IF (status /= 0) call announce_error(0,  &
         & MLSMSG_DeAllocate//"field")    


  END SUBROUTINE DestroyGridTemplateContents

  ! --------------------------------------------------------------------------

  ! This subroutine adds a quantity template to a database, or creates the
  ! database if it doesn't yet exist

!  SUBROUTINE AddGridTemplateToDatabase(database,qty)
  INTEGER FUNCTION AddGridTemplateToDatabase(database,item)

    ! Dummy arguments
    TYPE (GriddedData_T), DIMENSION(:), POINTER :: database
    TYPE (GriddedData_T), INTENT(IN) :: item
!    TYPE (GriddedData_T), INTENT(IN) :: qty

    ! Local variables
    TYPE (GriddedData_T), DIMENSION(:), POINTER :: tempDatabase
!    INTEGER :: newSize,status

    ! Executable code

!    IF (ASSOCIATED(database)) THEN
       ! Check we don't already have one of this name
!       IF (LinearSearchStringArray(database%quantityName, qty%quantityName, &
!            & caseInsensitive=.TRUE.)/=0) CALL MLSMessage(MLSMSG_Error,&
!            & ModuleName,MLSMSG_Duplicate//qty%quantityName)
!       newSize=SIZE(database)+1
!    ELSE
!       newSize=1
!    ENDIF
!    ALLOCATE(tempDatabase(newSize),STAT=status)
!    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
!         & "Allocation failed for tempDatabase")

!    IF (newSize>1) tempDatabase(1:newSize-1)=database
!    tempDatabase(newSize)=qty
!    IF (ASSOCIATED(database)) THEN
!       DEALLOCATE(database, STAT=status)
!       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
!         & MLSMSG_DeAllocate//"database")
!    end if
!    database=>tempDatabase

    include "addItemToDatabase.f9h"
    AddGridTemplateToDatabase = newSize

  END FUNCTION AddGridTemplateToDatabase

  ! --------------------------------------------------------------------------

  ! This subroutine destroys a quantity template database

  SUBROUTINE DestroyGridTemplateDatabase(database)

    ! Dummy argument
    TYPE (GriddedData_T), DIMENSION(:), POINTER :: database

    ! Local variables
    INTEGER :: qtyIndex, status

    IF (ASSOCIATED(database)) THEN
       DO qtyIndex=1,SIZE(database)
          CALL DestroyGridTemplateContents(database(qtyIndex))
       ENDDO
       DEALLOCATE(database, stat=status)
       IF (status /= 0) call announce_error(0,  &
         & MLSMSG_DeAllocate//"database")
    ENDIF
  END SUBROUTINE DestroyGridTemplateDatabase


!----------------- Beginning of Paul's code ------------------

    !---------------------------- ReadGriddedData ---------------------
  SUBROUTINE ReadGriddedData(FileName, lcf_where, description, v_type, &
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

    CHARACTER (LEN=*), INTENT(IN) :: FileName ! Name of the file containing the grid(s)
    INTEGER, INTENT(IN) :: lcf_where			! node of the lcf that provoked me
    INTEGER, INTENT(IN) :: v_type			! vertical coordinate; an 'enumerated' type
     TYPE( GriddedData_T ), INTENT(OUT) :: the_g_data ! Result
    CHARACTER (LEN=*), INTENT(IN) :: description ! e.g., 'dao'
    CHARACTER (LEN=*), OPTIONAL, INTENT(IN) :: GeoDimList ! Comma-delimited dim names
	 CHARACTER (LEN=*), OPTIONAL, INTENT(IN) :: fieldName ! Name of gridded field

    ! Local Variables

  integer :: edges(4)
  integer :: file_id, gd_id
  integer :: inq_success
  integer :: i
  integer :: nentries, ngrids, ndims, nfields
  integer :: strbufsize
!  character (len=80) :: msg, mnemonic
!  integer :: status

    LOGICAL,  PARAMETER       :: CASESENSITIVE = .FALSE.
  integer, parameter :: GRIDORDER=1				! What order grid written to file
  integer, parameter :: MAXLISTLENGTH=LineLen		! Max length list of grid names
  integer, parameter :: NENTRIESMAX=20		   ! Max num of entries
  character (len=MAXLISTLENGTH) :: gridlist
  character (len=MAXLISTLENGTH) :: dimlist, actual_dim_list
  character (len=MAXLISTLENGTH), DIMENSION(1) :: dimlists
  character (len=MAXLISTLENGTH) :: fieldlist
  integer, parameter :: MAXNAMELENGTH=NameLen		! Max length of grid name
  character (len=MAXNAMELENGTH) :: gridname, actual_field_name, the_dim
  INTEGER, DIMENSION(NENTRIESMAX) :: dims, rank, numberTypes, start, stride
  INTEGER                        :: our_rank, size, numberType
	!                                  These start out initialized to one
  INTEGER                        :: nlon=1, nlat=1, nlev=1, ntime=1
  INTEGER, PARAMETER             :: i_longitude=1
  INTEGER, PARAMETER             :: i_latitude=i_longitude+1
  INTEGER, PARAMETER             :: i_vertical=i_latitude+1
  INTEGER, PARAMETER             :: i_time=i_vertical+1
  ! External functions
!  integer, external :: gdopen, gdattach, gdrdfld, gddetach, gdclose
!  integer, external :: gdinqgrid, gdnentries, gdinqdims, gdinqflds
  INTEGER, EXTERNAL :: GDRDFLD
  logical, parameter :: COUNTEMPTY=.TRUE.
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
	
	if(descrpt_is_misplcd) THEN
		call announce_error(lcf_where, 'READGriddedData called with climatology' &
		& // ' description')
		return
	elseif(.NOT. descrpt_is_legal) then
		call announce_error(lcf_where, 'READGriddedData called with unknown' &
		& // ' description: ' // description)
		return
	endif
	
    error = 0
  file_id = gdopen(FileName, DFACC_RDONLY)

  IF (file_id < 0) THEN
	CALL announce_error(lcf_where, "Could not open "// FileName)
  END IF

! Find list of grid names on this file
  inq_success = gdinqgrid(FileName, gridlist, strbufsize)
  IF (inq_success < 0) THEN
	CALL announce_error(lcf_where, "Could not inquire gridlist "// FileName)
  END IF

! Find grid name corresponding to the GRIDORDER'th one
	ngrids = NumStringElements(gridlist, COUNTEMPTY)
	
	IF(ngrids <= 0) THEN
		CALL announce_error(lcf_where, "NumStringElements of gridlist <= 0")
	ELSEIF(ngrids /= inq_success) THEN
		CALL announce_error(lcf_where, "NumStringElements of gridlist /= inq_success")
	ELSEIF(ngrids < GRIDORDER) THEN
		CALL announce_error(lcf_where, "NumStringElements of gridlist < GRIDORDER")
	ENDIF
	
	CALL GetStringElement(gridlist, gridname, GRIDORDER, COUNTEMPTY)

  gd_id = gdattach(file_id, gridname)
  IF (gd_id < 0) THEN
		CALL announce_error(lcf_where, "Could not attach "//FileName)
  END IF

! Now find dimsize(), dimname(), etc.
	nentries = gdnentries(gd_id, HDFE_NENTDIM, strbufsize)

	IF(nentries <= 0) THEN
		CALL announce_error(lcf_where, "nentries of gd_id <= 0")
	ELSEIF(nentries > NENTRIESMAX) THEN
		CALL announce_error(lcf_where, "nentries of gd_id > NENTRIESMAX")
	ENDIF

	ndims = gdinqdims(gd_id, dimlist, dims)

	IF(ndims <= 0) THEN
		CALL announce_error(lcf_where, "ndims of gd_id <= 0")
	ELSEIF(ndims > NENTRIESMAX) THEN
		CALL announce_error(lcf_where, "ndims of gd_id > NENTRIESMAX")
	ENDIF

	nfields = gdinqflds(gd_id, fieldlist, rank, numberTypes)

	IF(nfields <= 0) THEN
		CALL announce_error(lcf_where, "nfields of gd_id <= 0")
	ELSEIF(nfields > NENTRIESMAX) THEN
		CALL announce_error(lcf_where, "nfields of gd_id > NENTRIESMAX")
	ENDIF
	
	IF(.NOT. CASESENSITIVE) THEN
		fieldlist = Capitalize(fieldlist)
	ENDIF

	IF(PRESENT(fieldName)) THEN
		actual_field_name=fieldName
	ELSE
		actual_field_name=DEFAULTFIELDNAME
	ENDIF

	IF(PRESENT(GeoDimList)) THEN
		actual_dim_list=GeoDimList
	ELSE
		actual_dim_list=GEO_FIELD1 // ',' // &
		& GEO_FIELD2 // ',' // &
		& GEO_FIELD3 // ',' // &
		& GEO_FIELD4
	ENDIF

	! Now find the rank of our field
	
	inq_success = gdfldinfo(gd_id, TRIM(actual_field_name), our_rank, dims, &
	& numbertype, dimlists(1))

	dimlist = TRIM(dimlists(1))

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
  END SUBROUTINE ReadGriddedData
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

    logical :: end_of_file = .FALSE.

    do CliUnit = mlspcf_l2clim_start, mlspcf_l2clim_end

!     Open one Climatology file as a generic file for reading
      version = 1
      returnStatus = Pgs_io_gen_openF ( CliUnit, PGSd_IO_Gen_RSeqFrm, 0, &
                                        processCli, version )
      if ( returnStatus == PGS_S_SUCCESS ) then

      do while (.NOT. end_of_file)

        call l3ascii_read_field ( processCli, qty, end_of_file)
        returnStatus = AddGridTemplateToDatabase(aprioriData, qty)
        call DestroyGridTemplateContents ( qty )

      end do !(.not. end_of_file)
		
		end_of_file = .FALSE.

      end if

    end do ! CliUnit = mlspcf_l2clim_start, mlspcf_l2clim_end

  return
  !============================
  end subroutine OBTAIN_CLIM
  !============================


  ! --------------------------------------------------  READ_CLIMATOLOGY  -----
  SUBROUTINE READ_CLIMATOLOGY ( fname, root, aprioriData, &
  & mlspcf_l2clim_start, mlspcf_l2clim_end, echo_data, dump_data )
  ! --------------------------------------------------
  ! Brief description of program
  ! This subroutine reads a l3ascii file and returns
  ! the data_array to the caller

  ! Arguments

  character*(*), intent(in) :: fname			! Physical file name
    type (GriddedData_T), dimension(:), pointer :: aprioriData 
    integer, intent(in) :: ROOT        ! Root of the L2CF abstract syntax tree
	 INTEGER, OPTIONAL, INTENT(IN) :: mlspcf_l2clim_start, mlspcf_l2clim_end
    LOGICAL, OPTIONAL, intent(in) :: echo_data        ! echo climatology quantity name
    LOGICAL, OPTIONAL, intent(in) :: dump_data        ! dump climatology data
	 
	 ! Local
    type (GriddedData_T)        :: gddata 
	 INTEGER :: ErrType
	 INTEGER, PARAMETER :: version=1
	 LOGICAL :: end_of_file
	 LOGICAL, PARAMETER :: debug=.FALSE.
	 LOGICAL, PARAMETER :: ECHO_GRIDDED_QUANTITIES=.TRUE.
	 LOGICAL, PARAMETER :: DUMP_GRIDDED_QUANTITIES=.TRUE.
	 LOGICAL :: echo
	 LOGICAL :: dump
    integer:: processCli, CliUnit, record_length
	 
	! begin
	end_of_file=.FALSE.
	if(present(echo_data)) then
		echo = echo_data
	ELSE
		echo = ECHO_GRIDDED_QUANTITIES
	ENDIF
	
	if(present(dump_data)) then
		dump = dump_data
	ELSE
		dump = DUMP_GRIDDED_QUANTITIES
	ENDIF
	
	! use PCF

	if(present(mlspcf_l2clim_start) .and. present(mlspcf_l2clim_end)) then

	 CliUnit = GetPCFromRef(fname, mlspcf_l2clim_start, mlspcf_l2clim_end, &
  & .TRUE., ErrType, version)
  
  IF(ErrType /= 0) THEN
!    CALL MLSMessage (MLSMSG_Error, ModuleName, &
!              &"Climatology file name unmatched in PCF")
    CALL announce_error (ROOT, &
              &"Climatology file name " // fname // " unmatched in PCF")
	RETURN
  ENDIF

      ErrType = Pgs_io_gen_openF ( CliUnit, PGSd_IO_Gen_RSeqFrm, 0, &
                                        processCli, version )

	! use Fortran open
	else
	
		if(debug) call output('opening ' // fname, advance = 'yes')

		CliUnit = mls_io_gen_openF ( 'open', .true., ErrType, &
	& record_length, PGSd_IO_Gen_RSeqFrm, FileName=fname)
	
	endif

	if(debug) then
		if(.NOT. end_of_file) then
			call output('Not yet eof on io unit', advance = 'yes')
		else
			call output('Starting at eof on io unit', advance = 'yes')
		endif
	endif


      if ( ErrType == PGS_S_SUCCESS ) then

      do while (.NOT. end_of_file)

			if(debug) call output('reading l3ascii file', advance = 'yes')

        call l3ascii_read_field ( CliUnit, gddata, end_of_file, ErrType)

			if(ErrType == 0) then
		  		if(debug) then
					call output('adding to grid database', advance='yes')
					call output('adding grid template to database ', advance='yes')
				endif
		 	 if(echo .OR. debug) then
				call output('quantity name ' // gddata%quantityName, advance='yes')
				call output('description ' // gddata%description, advance='yes')
				call output('units ' // gddata%units, advance='yes')
			endif

		  if(dump) then
			call Dump_Gridded_Data(gddata, root)
			endif

        ErrType = AddGridTemplateToDatabase(aprioriData, gddata)

			if(debug) call output('Destroying our grid template', advance='yes')
			
		endif
        call DestroyGridTemplateContents ( gddata )

      end do !(.not. end_of_file)
		
	! ok, done with this file and unit number
	if(present(mlspcf_l2clim_start) .and. present(mlspcf_l2clim_end)) then
      ErrType = Pgs_io_gen_CloseF ( CliUnit )

	! use Fortran close
	else
	
		if(debug) call output('closing ' // fname, advance = 'yes')
		ErrType = mls_io_gen_CloseF ('close', CliUnit )
		
	endif

	if(ErrType /= 0) then
    		CALL announce_error (ROOT, &
              &"Error closing " // fname, error_number=ErrType)
	endif
	
		else

    		CALL announce_error (ROOT, &
              &"Error opening " // fname, error_number=ErrType)
		endif

	END SUBROUTINE READ_CLIMATOLOGY

  ! -------------------------------------------------  OBTAIN_DAO  -----
  subroutine OBTAIN_DAO ( aprioriData, root, &
  & mlspcf_l2dao_start, mlspcf_l2dao_end )
 
	! An atavism--
	! a throwback to when ncep files were opened
	! independently of being required by the lcf

    type (GriddedData_T), dimension(:), pointer :: aprioriData 
    ! Input a priori database
    integer, intent(in) :: ROOT        ! Root of L2CF abstract syntax tree
	 INTEGER, INTENT(IN) :: mlspcf_l2dao_start, mlspcf_l2dao_end

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
			IF(returnStatus > 0) THEN
				call ReadGriddedData ( DAOphysicalFilename, root, &
				& 'dao', v_is_pressure, aprioriData(returnStatus) )
			ENDIF

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
	 INTEGER, INTENT(IN) :: mlspcf_l2ncep_start, mlspcf_l2ncep_end

    ! Local Variables
!    character (len=80) :: MSG, MNEMONIC
    integer :: NCEPFileHandle, NCEP_Version
    character (len=132) :: NCEPphysicalFilename
    type (GriddedData_T):: QTY
    integer :: RETURNSTATUS
!    character (len=80) :: VNAME

!   allocate ( data_array(XDIM, YDIM, ZDIM), stat=returnStatus )

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
			IF(returnStatus > 0) THEN
				call ReadGriddedData ( NCEPphysicalFilename, root, &
				& 'ncep', v_is_pressure, aprioriData(returnStatus) )
			ENDIF

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
    logical, parameter :: MAYDUMPFIELDVALUES = .FALSE.
	 integer            :: i

	if ( .NOT. associated(GriddedData)) then
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
    logical, parameter :: MAYDUMPFIELDVALUES = .FALSE.

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

		CALL output("Caused the following error:", advance='yes', &
		& from_where=ModuleName)
		CALL output(trim(full_message), advance='yes', &
		& from_where=ModuleName)
		if(present(error_number)) then
			CALL output('error number ', advance='no')
			CALL output(error_number, places=9, advance='yes')
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
END MODULE ncep_dao
!=============================================================================

!
! $Log$
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
