! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===========================================================
module ObtainNCEP !provides subroutines to access NCEP files
!===========================================================

! use DATES_MODULE
  use GriddedData, only: GriddedData_T, AddGridTemplateToDatabase
  use Hdf, only: DFACC_RDONLY, FAIL, SUCCEED
  use HDFEOS, only: HDFE_NENTDIM, HDFE_NENTDFLD
  use MLSCommon, only: R8
  use MLSPCF2, only: MLSPCF_L2NCEP_END, MLSPCF_L2NCEP_START
  use MLSStrings, only: GetStringElement, NumStringElements
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use SDPToolkit, only: Pgs_pc_getReference, PGS_S_SUCCESS
! use VerticalCoordinate
! use ???, only: Pgs_smf_getMsg

  implicit none
  private
  public :: Obtain_NCEP, READ_NCEP

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = & 
     "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  !Parameters
  integer, parameter :: XDIM=15
  integer, parameter :: YDIM=181
  integer, parameter :: ZDIM=360
  character (len=4), parameter :: LevelhPa(15)=(/'1000', '850 ','700 ','500 ', '400 ', &
         '300 ', '250 ', '200 ', '150 ', '100 ', '70  ', '50  ', '30  ','20  ','10  '/)

  integer, parameter :: LevelhPaI(15)=(/1000, 850, 700, 500, 400 , &
         300, 250, 200, 150, 100, 70, 50, 30, 20, 10/)

contains ! =====     Public Procedures     =============================

  ! ------------------------------------------------  Obtain_NCEP  -----
  subroutine Obtain_NCEP ( aprioriData, root )

    ! Arguments
    type (GriddedData_T), dimension(:), pointer :: aprioriData 
    ! Input a priori database
    integer, intent(in) :: ROOT        ! Root of the L2CF abstract syntax tree

    ! Local Variables
    real(R8) :: DATA_ARRAY(xdim, ydim, zdim)
!   integer :: I
    character (len=80) :: MSG, MNEMONIC
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

        call read_ncep ( NCEPphysicalFilename, data_array )

!       Create a GriddedDtata Template and copy data_array into it

!        CALL SetupNewGridTemplate(qty, noHeights=XDIM, noLats=YDIM, noLons=ZDIM,&
!                                  noLsts=1,noSzas=1, noDates=1)
!        do i = 1,zdim
!          qty%heights(i) = LevelhPaI(i)
!        end do

!        do i=1, ydim
!          qty%lats(i) = -91.0 + i
!        end do

!        do i = 1, ZDIM
!          qty%lons(i) = -181.0 +i
!        end do

!        qty%field(:,:,:,1,1,1) = data_array(:,:,:)
!       DEALLOCATE (data_array, stat=returnStatus)
 
        returnStatus = AddGridTemplateToDatabase(aprioriData, qty)
      else

! This isn't necessarily an error--we just don't need 30 different ncep files
! reserved by the MLSPCF2 file
!        call Pgs_smf_getMsg ( returnStatus, mnemonic, msg )
!        call MLSMessage (MLSMSG_Error, ModuleName, &
!                         "Error opening NCEP file:  "//mnemonic//" "//msg)
      end if

    end do ! NCEPFileHandle = mlspcf_l2_ncep_start, mlspcf_l2_ncep_end
!===========================
  end subroutine Obtain_NCEP
!===========================

  ! --------------------------------------------------  READ_NCEP  -----
  SUBROUTINE READ_NCEP ( fname, data_array )
  ! --------------------------------------------------
  ! Brief description of program
  ! This subroutine reads a NCEP correlative file and returns
  ! the data_array to the caller

  ! Arguments

  character*(*), intent(in) :: fname			! Physical file name
  real(R8) ::  data_array(:,:,:)

  ! - - - local declarations - - -
  integer :: edges(4)
  integer :: file_id, gd_id
  integer :: inq_success
  integer :: i
  integer :: nentries, ngrids, ndims, nfields
  integer :: strbufsize
  character (len=80) :: msg, mnemonic
  integer :: start(4)
  integer :: status
  integer :: stride(4)

  integer, parameter :: GRIDORDER=1				! What order grid written to file
  integer, parameter :: MAXLISTLENGTH=80		! Max length list of grid names
  integer, parameter :: NENTRIESMAX=20		   ! Max num of entries
  character (len=MAXLISTLENGTH) :: gridlist
  character (len=MAXLISTLENGTH) :: dimlist
  character (len=MAXLISTLENGTH) :: fieldlist
  integer, parameter :: MAXNAMELENGTH=16		! Max length of grid name
  character (len=MAXNAMELENGTH) :: gridname
  INTEGER, DIMENSION(NENTRIESMAX) :: dims, rank, numberType
  ! External functions
  integer, external :: gdopen, gdattach, gdrdfld, gddetach, gdclose
  integer, external :: gdinqgrid, gdnentries, gdinqdims, gdinqflds
  logical, parameter :: COUNTEMPTY=.TRUE.
 
  ! - - - begin - - -

  file_id = gdopen(fname, DFACC_RDONLY)

  IF (file_id /= SUCCEED) THEN
    CALL Pgs_smf_getMsg(status, mnemonic, msg)
    CALL MLSMessage (MLSMSG_Error, ModuleName, &
              &"Could not open "// fname//" "//mnemonic//" "//msg)
  END IF

! Find list of grid names on this file
  inq_success = gdinqgrid(fname, gridlist, strbufsize)
  IF (inq_success /= SUCCEED) THEN
    CALL Pgs_smf_getMsg(status, mnemonic, msg)
    CALL MLSMessage (MLSMSG_Error, ModuleName, &
              &"Could not inquire gridlist "// fname//" "//mnemonic//" "//msg)
  END IF

! Find grid name corresponding to the GRIDORDER'th one
	ngrids = NumStringElements(gridlist, COUNTEMPTY)
	
	IF(ngrids <= 0) THEN
    CALL MLSMessage (MLSMSG_Error, ModuleName, &
              &"NumStringElements of gridlist <= 0")
	ELSEIF(ngrids < GRIDORDER) THEN
    CALL MLSMessage (MLSMSG_Error, ModuleName, &
              &"NumStringElements of gridlist < GRIDORDER")
	ENDIF
	
	CALL GetStringElement(gridlist, gridname, GRIDORDER, COUNTEMPTY)

  gd_id = gdattach(file_id, gridname)
  IF (gd_id /= SUCCEED) THEN
    CALL Pgs_smf_getMsg(status, mnemonic, msg)
    CALL MLSMessage (MLSMSG_Error, ModuleName, &
               "Could not attach "//fname//" "//mnemonic//" "//msg)
  END IF

! Now find dimsize(), dimname(), etc.
	nentries = gdnentries(gd_id, HDFE_NENTDIM, strbufsize)

	IF(nentries <= 0) THEN
    CALL MLSMessage (MLSMSG_Error, ModuleName, &
              &"nentries of gd_id <= 0")
	ELSEIF(nentries > NENTRIESMAX) THEN
    CALL MLSMessage (MLSMSG_Error, ModuleName, &
              &"nentries of gd_id > NENTRIESMAX")
	ENDIF

	ndims = gdinqdims(gd_id, dimlist, dims)

	IF(ndims <= 0) THEN
    CALL MLSMessage (MLSMSG_Error, ModuleName, &
              &"ndims of gd_id <= 0")
	ELSEIF(ndims > NENTRIESMAX) THEN
    CALL MLSMessage (MLSMSG_Error, ModuleName, &
              &"ndims of gd_id > NENTRIESMAX")
	ENDIF

	nfields = gdinqflds(gd_id, fieldlist, rank, numberType)

	IF(nfields <= 0) THEN
    CALL MLSMessage (MLSMSG_Error, ModuleName, &
              &"nfields of gd_id <= 0")
	ELSEIF(nfields > NENTRIESMAX) THEN
    CALL MLSMessage (MLSMSG_Error, ModuleName, &
              &"nfields of gd_id > NENTRIESMAX")
	ENDIF

  start(1) = 0
  start(2) = 0

  stride(1) = 1
  stride(2) = 1

  edges(1) = 360
  edges(2) = 181


!   In this subroutine, we read the entire field.  By manipulating the start 
!   and edges arrays, it is possible to read a subset of the entire array.  
!   For example, to read a 3D section defined by x=100,224 y=50,149 
!   z=15,16 you would set the start and edges arrays to the following:

!   start(1) = 0    time start location
!   start(2) = 15   z-dim start location
!   start(3) = 50   y-dim start location
!   start(4) = 100  x-dim start location

!   edges(1) = 1    time length
!   edges(2) = 2    z-dim length
!   edges(3) = 100  y-dim length
!   edges(4) = 125  x-dim length

! NCEP has 15 pressure levels of temperature, see LevelhPa above

  do i = 1, 15
     status = gdrdfld(gd_id, 'ISOBARIC LEVEL AT '//TRIM(LevelhPa(i))//' (hPa)', start,&
                      stride,edges, data_array(i,:,:))

     IF(status /= PGS_S_SUCCESS) THEN
        CALL Pgs_smf_getMsg(status, mnemonic, msg)
        CALL MLSMessage (MLSMSG_Error, ModuleName, &
                      "Error reading grid field:  "//mnemonic//" "//msg)
     END IF
  end do

  status = gddetach(gd_id)
  IF(status /= PGS_S_SUCCESS) THEN
     CALL Pgs_smf_getMsg(status, mnemonic, msg)
     CALL MLSMessage (MLSMSG_Error, ModuleName, &
                      "Error detaching  grid:  "//mnemonic//" "//msg)
  END IF

  status = gdclose(file_id)
  IF(status /= PGS_S_SUCCESS) THEN
     CALL Pgs_smf_getMsg(status, mnemonic, msg)
     CALL MLSMessage (MLSMSG_Error, ModuleName, &
                      "Error closing  NCEP File:  "//mnemonic//" "//msg)
  END IF

  RETURN
  END SUBROUTINE read_ncep

! =====     Private Procedures     =====================================

END MODULE ObtainNCEP

! $Log$
! Revision 2.6  2001/03/03 00:11:29  pwagner
! Began transformations to act like L2GPData module for Gridded data
!
! Revision 2.5  2001/02/23 17:44:50  pwagner
! Corrected num of args to GetStringElement
!
! Revision 2.4  2001/02/23 00:06:52  pwagner
! Using some MLSStrings.f90 functions
!
! Revision 2.3  2001/02/21 00:37:51  pwagner
! Uses more of GriddedData
!
! Revision 2.2  2001/02/13 22:59:36  pwagner
! l2 modules can only use MLSPCF2
!
! Revision 2.1  2000/10/12 00:34:56  vsnyder
! Comment-out apparently unnecessary USEs; add "only" to the others
!
! Revision 2.0  2000/09/05 18:57:06  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry
!

