! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

program test
  use GriddedData, only: GriddedData_T, v_is_pressure
  use ncep_dao, only: AddGridTemplateToDatabase, &
  & READ_CLIMATOLOGY, ReadGriddedData
  use Hdf, only: DFACC_READ, DFACC_RDONLY
  use HDFEOS, only: HDFE_NENTDIM, HDFE_NENTDFLD, &
  & gdopen, gdattach, gddetach, gdclose, gdinqattrs, &
  & gdinqgrid, gdnentries, gdinqdims, gdinqflds, gddiminfo, &
  & swopen, swclose
  USE MLSCommon, only: R8, LineLen, NameLen
  USE MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Allocate, &
  & MLSMSG_Deallocate, MLSMSG_Warning
  USE MLSStrings, only: GetStringElement, NumStringElements, Capitalize, &
  & Count_words, ReadCompleteLineWithoutComments, GetIntHashElement

  implicit NONE
  
  integer, external :: gdfldinfo

! Local
	CHARACTER(LEN=*), PARAMETER :: FileName = 'test_ncep.data'
	CHARACTER(LEN=*), PARAMETER :: description = 'ncep'
!	CHARACTER(LEN=*), PARAMETER :: FileName = 'MLS-Aura_L3DM_ClO_V0-5-C01_1996-044.dat'
	CHARACTER(LEN=*), PARAMETER :: path = '/nas/user2/pwagner/mlspgs/tests/l2/'
    type (GriddedData_T) :: gddata 
	 INTEGER, parameter :: unit=1
	 INTEGER :: iostat, i, file_id
	logical :: include_path
!	integer, external :: swopen

do i=1, 1

	include_path = i==2
	
	if(include_path) THEN
		print *, ' ' 
		print *, 'Including path: ', path
		open ( unit=unit, file=path//filename, status="old", action="read", iostat=iostat )
	else
		print *, ' ' 
		print *, 'Excluding path'
		open ( unit=unit, file=filename, status="old", action="read", iostat=iostat )
	endif

	if(iostat /= 0) then
		print *, 'Problem opening ' //FileName
	else
		close (unit)
	endif

	if(include_path) THEN
          CALL ReadGriddedData(path//FileName, 0, description, v_is_pressure, &
            & gddata, 'XDim,YDim,Height,TIME', 'Some_field_name')
	else
          CALL ReadGriddedData(FileName, 0, description, v_is_pressure, &
            & gddata, 'XDim,YDim,Height,TIME', 'Some_field_name')
	endif
	
	if(include_path) THEN
          file_id = swopen(path//FileName, DFACC_READ)
	else
          file_id = swopen(FileName, DFACC_READ)
	endif

	if(file_id == -1) then
		print *, 'Problem opening swath' //FileName
	else
		iostat = swclose(file_id)
	endif
	
enddo

          CALL LookAtGriddedData(FileName, &
            & gddata, 'XDim,YDim,Height,TIME', 'TMPU')

	CONTAINS

    !---------------------------- LookAtGriddedData ---------------------
  SUBROUTINE LookAtGriddedData(FileName, the_g_data, GeoDimList, fieldName)
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
     TYPE( GriddedData_T ), INTENT(OUT) :: the_g_data ! Result
    CHARACTER (LEN=*), OPTIONAL, INTENT(IN) :: GeoDimList ! Comma-delimited dim names
	 CHARACTER (LEN=*), OPTIONAL, INTENT(IN) :: fieldName ! Name of gridded field

    ! Local Variables

  integer :: edges(4)
  integer :: file_id, gd_id
  integer :: inq_success
  integer :: i
  integer :: nentries, ngrids, ndims, nfields, nattrs
  integer :: strbufsize

    LOGICAL,  PARAMETER       :: CASESENSITIVE = .FALSE.
  integer, parameter :: GRIDORDER=1				! What order grid written to file
  integer, parameter :: MAXLISTLENGTH=1024		! Max length list of grid names
  integer, parameter :: NENTRIESMAX=40		   ! Max num of entries
  character (len=MAXLISTLENGTH) :: gridlist
  character (len=MAXLISTLENGTH) :: dimlist, actual_dim_list
  character (len=MAXLISTLENGTH), DIMENSION(NENTRIESMAX) :: dimlists
  character (len=MAXLISTLENGTH) :: attrlist
  character (len=MAXLISTLENGTH) :: fieldlist
  integer, parameter :: MAXNAMELENGTH=NameLen		! Max length of grid name
  character (len=MAXNAMELENGTH) :: gridname, actual_field_name, the_dim
  INTEGER, DIMENSION(NENTRIESMAX) :: dims, rank, start, stride, numberTypes
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

   CHARACTER (len=*), PARAMETER :: DEFAULTFIELDNAME = 'TMPU'
   CHARACTER (len=*), PARAMETER :: GEO_FIELD1 = 'Latitude'
   CHARACTER (len=*), PARAMETER :: GEO_FIELD2 = 'Longitude'
   CHARACTER (len=*), PARAMETER :: GEO_FIELD3 = 'Height'
   CHARACTER (len=*), PARAMETER :: GEO_FIELD4 = 'Time'

  ! - - - begin - - -

  file_id = gdopen(FileName, DFACC_RDONLY)

	CALL TellMe('FileName', FileName)

! Find list of grid names on this file
  inq_success = gdinqgrid(FileName, gridlist, strbufsize)
	CALL TellMe('gridlist', gridlist)

! Find grid name corresponding to the GRIDORDER'th one
	ngrids = NumStringElements(gridlist, COUNTEMPTY)
	CALL TellMe('ngrids', intvalue=ngrids)
		
	CALL GetStringElement(gridlist, gridname, GRIDORDER, COUNTEMPTY)
	CALL TellMe('gridname', gridname)

  gd_id = gdattach(file_id, gridname)

! Now find dimsize(), dimname(), etc.
	nentries = gdnentries(gd_id, HDFE_NENTDIM, strbufsize)

	ndims = gdinqdims(gd_id, dimlist, dims)
	CALL TellMe('ndims', intvalue=ndims)
	CALL TellMe('dimlist', dimlist)

	nattrs = gdinqattrs(gd_id, attrlist, strbufsize)
	CALL TellMe('nattrs', intvalue=nattrs)
	CALL TellMe('attrlist', attrlist)

	nfields = gdinqflds(gd_id, fieldlist, rank, numberTypes)
	CALL TellMe('nfields', intvalue=nfields)
	CALL TellMe('fieldlist', fieldlist)

	CALL TellMe('The ranks of each field')
	CALL TellMe('field                             rank')
	DO i=1, nfields
		CALL GetStringElement(fieldlist, actual_field_name, i, COUNTEMPTY)
		CALL TellMe('rank(' // TRIM(actual_field_name) // ')', intValue=rank(i))
	ENDDO

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
	
!	our_rank = GetIntHashElement(fieldlist, rank, actual_field_name, inq_success, &
!  & countEmpty)

	CALL TellMe('Calling gdfldinfo for field ', 'no_such_field')

	inq_success = gdfldinfo(gd_id, 'no_such_field', our_rank, dims, &
	& numbertype, dimlist)

	CALL TellMe('gdfield info', intValue=inq_success)

	CALL TellMe('Calling gdfldinfo using dimlists(1) for field ', actual_field_name)

	inq_success = gdfldinfo(gd_id, TRIM(actual_field_name), our_rank, &
	& dims, &
	& numbertype, dimlists(1))

	dimlist = TRIM(dimlists(1))

	CALL TellMe('gdfield info', intValue=inq_success)

	CALL TellMe('our_rank', intValue=our_rank)
	CALL TellMe('numbertype', intValue=numbertype)
	CALL TellMe('dimlist', dimlist)

	CALL TellMe('The size of each dim of field', actual_field_name)
	CALL TellMe('dim                             size')
	DO i=1, our_rank
		CALL GetStringElement(dimlist, the_dim, i, COUNTEMPTY)
		CALL TellMe('size(' // TRIM(the_dim) // ')', intValue=dims(i))
	ENDDO

	nlon = dims(1)
	nlat = dims(2)
	nlev = dims(3)
	ntime = dims(4)
		  
	the_g_data%quantityName = actual_field_name
	  
	the_g_data%noLons = nlon
	the_g_data%noLats = nlat
	the_g_data%noHeights = nlev
	the_g_data%noLsts = ntime

  !-----------------------------
  END SUBROUTINE LookAtGriddedData
  !-----------------------------

	SUBROUTINE TellMe(What, StringValue, IntValue)

	! Args
	CHARACTER (LEN=*) , INTENT(IN) :: What
	CHARACTER (LEN=*) , OPTIONAL, INTENT(IN) :: StringValue
	INTEGER , OPTIONAL, INTENT(IN) :: IntValue
	
	! Private
	CHARACTER (LEN=LINELEN) :: theLine
	
	IF(PRESENT(StringValue)) THEN
		theLine = TRIM(adjustl(What)) // ' = ' // TRIM(adjustl(StringValue))
		print *, theLine
	ELSEIF(PRESENT(IntValue)) THEN
		print *, TRIM(adjustl(What)), ' = ', IntValue
	ELSE
		print *, TRIM(adjustl(What))
	ENDIF

	END SUBROUTINE TellMe

end program test

! $Log$
!
! MLS-Aura_L3DM_ClO_V0-5-C01_1996-044.dat
! /nas/user2/pwagner/mlspgs/tests/l2
