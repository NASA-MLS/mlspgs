! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===========================================================
module ObtainNCEP !provides subroutines to access NCEP files
!===========================================================

  use dates_module
  use GriddedData
  use Hdf
  use MLSCommon
  use MLSPCF, only: MLSPCF_L2NCEP_END, MLSPCF_L2NCEP_START
  use MLSStrings
  use MLSMessageModule
  use VerticalCoordinate
! use ???, only: Pgs_smf_getMsg

  implicit none
  private
  public :: Obtain_NCEP

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
    type (griddeddata_t), dimension(:), pointer :: aprioriData 
    ! Input a priori database
    integer, intent(in) :: ROOT        ! Root of the L2CF abstract syntax tree

    ! Local Variables
    real :: DATA_ARRAY(xdim, ydim, zdim)
!   integer :: I
    character (len=80) :: MSG, MNEMONIC
    integer :: NCEPFileHandle, NCEP_Version
    character (len=132) :: NCEPphysicalFilename
!   type (GriddedData_T):: QTY
    integer :: RETURNSTATUS
    character (len=80) :: VNAME

!   allocate ( data_array(XDIM, YDIM, ZDIM), stat=returnStatus )

    NCEP_Version = 1
    vname = "TMP_3" ! for now
! Get the NCEP file name from the PCF

    do NCEPFileHandle = mlspcf_l2ncep_start, mlspcf_l2ncep_end

      returnStatus = Pgs_pc_getReference ( NCEPFileHandle, NCEP_Version, &
                                           NCEPphysicalFilename )

      if ( returnStatus == PGS_S_SUCCESS ) then

! Open the HDF-EOS file and read gridded data

        call read_ncep ( NCEPphysicalFilename, vname, data_array )

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
 
!        CALL AddGridTemplateToDatabase(aprioriData, qty)
      else

        call Pgs_smf_getMsg ( returnStatus, mnemonic, msg )
        call MLSMessage (MLSMSG_Error, ModuleName, &
                         "Error opening NCEP file:  "//mnemonic//" "//msg)
      end if

    end do ! NCEPFileHandle = mlspcf_l2_ncep_start, mlspcf_l2_ncep_end
!===========================
  end subroutine Obtain_NCEP
!===========================

! =====     Private Procedures     =====================================

  ! --------------------------------------------------  READ_NCEP  -----
  SUBROUTINE READ_NCEP ( fname, vname, data_array )
  ! --------------------------------------------------
  ! Brief description of program
  ! This subroutine reads a NCEP correlative file and returns
  ! the data_array to the caller

  ! Arguments

  character*(*), intent(in) :: fname, vname
  real ::  data_array(:,:,:)

  ! - - - local declarations - - -
  integer :: edges(4)
  integer :: file_id, gd_id
  integer :: i
  character (len=80) :: msg, mnemonic
  integer :: start(4)
  integer :: status
  integer :: stride(4)

  ! External functions
  integer :: gdopen, gdattach, gdrdfld, gddetach, gdclose
 
  ! - - - begin - - -

  file_id = gdopen(fname, DFACC_RDONLY)

  IF (file_id < 0) THEN
    CALL Pgs_smf_getMsg(status, mnemonic, msg)
    CALL MLSMessage (MLSMSG_Error, ModuleName, &
              &"Could not open "// fname//" "//mnemonic//" "//msg)
 
  END IF

  gd_id = gdattach(file_id, vname)
  IF (gd_id < 0) THEN
    CALL Pgs_smf_getMsg(status, mnemonic, msg)
    CALL MLSMessage (MLSMSG_Error, ModuleName, &
               "Could not attach "//fname//" "//mnemonic//" "//msg)

  END IF

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

END MODULE ObtainNCEP

! $Log$
! Revision 2.0  2000/09/05 18:57:06  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry
!

