! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===========================================================

MODULE ObtainNCEP !provides subroutines to access NCEP files

!===========================================================
USE MLSCommon
USE MLSStrings 
USE MLSMessageModule
USE GriddedData
USE VerticalCoordinate
USE MLSPCF
USE dates_module    
USE Hdf
USE MLSCF

IMPLICIT NONE

PUBLIC

PRIVATE :: Id, ModuleName
!------------------------------- RCS Ident Info ------------------------------
CHARACTER(LEN=130) :: id = & 
   "$Id$"
CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!-----------------------------------------------------------------------------

!Parameters
  INTEGER, PARAMETER :: XDIM=15
  INTEGER, PARAMETER :: YDIM=181
  INTEGER, PARAMETER :: ZDIM=360
  CHARACTER (LEN=4), PARAMETER :: LevelhPa(15)=(/'1000', '850 ','700 ','500 ', '400 ', &
         '300 ', '250 ', '200 ', '150 ', '100 ', '70  ', '50  ', '30  ','20  ','10  '/)

  INTEGER, PARAMETER :: LevelhPaI(15)=(/1000, 850, 700, 500, 400 , &
         300, 250, 200, 150, 100, 70, 50, 30, 20, 10/)
CONTAINS

!=======================================================
  SUBROUTINE read_ncep (fname, vname, data_array)
! --------------------------------------------------
! Brief description of program
! This subroutine reads a NCEP correlative file and returns
! the data_array to the caller


! Arguments

  CHARACTER*(*), INTENT(IN) :: fname, vname
  REAL*4 ::  data_array(:,:,:)
! - - - local declarations - - -
  
  INTEGER*4 sd_id, sds_id, status
  INTEGER*4 sds_index
  INTEGER*4 start(4), edges(4), stride(4)
  INTEGER i,j,k
  INTEGER file_id, gd_id
  CHARACTER (LEN=80) :: msg, mnemonic
  INTEGER :: sdstart, gdopen, gdattach, gdrdfld, gddetach, &
                       & gdclose, gdinqdims, gdinqflds
  integer strbufsize, ndim, nfld, dim(27),numbertype(27), rank(27)
  CHARACTER (LEN=1000) :: dimname, fieldlist
 
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

  SUBROUTINE Obtain_NCEP(aprioriData, l2cf_info)
 
    TYPE (GriddedData_T), DIMENSION(:), POINTER :: aprioriData 
    ! Input a priori database
    TYPE (MLSCF_T) :: l2cf_info        ! The information from the l2cf file


! Local Variables
 
    TYPE (GriddedData_T):: qty
    CHARACTER (LEN=80) :: msg, mnemonic
    INTEGER ::  returnStatus
    INTEGER :: NCEPFileHandle, NCEP_Version, sd_id, i
    CHARACTER (LEN=132) :: NCEPphysicalFilename
    CHARACTER (LEN=80) :: vname

    REAL*4, ALLOCATABLE :: data_array(:, :, :)
    ALLOCATE (data_array(XDIM, YDIM, ZDIM), stat=returnStatus)

    NCEP_Version = 1
    vname = "TMP_3" ! for now
! Get the NCEP file name from the PCF

    DO NCEPFileHandle = mlspcf_l2ncep_start, mlspcf_l2ncep_end

      returnStatus = Pgs_pc_getReference(NCEPFileHandle, NCEP_Version, NCEPphysicalFilename)

      IF(returnStatus == PGS_S_SUCCESS) THEN

! Open the HDF-EOS file and read gridded data

        CALL read_ncep (NCEPphysicalFilename, vname, data_array)

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
        DEALLOCATE (data_array, stat=returnStatus)
 
!        CALL AddGridTemplateToDatabase(aprioriData, qty)
      ELSE

        CALL Pgs_smf_getMsg(returnStatus, mnemonic, msg)
        CALL MLSMessage (MLSMSG_Error, ModuleName, &
                         "Error opening NCEP file:  "//mnemonic//" "//msg)
      END IF

    END DO ! NCEPFileHandle = mlspcf_l2_ncep_start, mlspcf_l2_ncep_end
!===========================
  END SUBROUTINE Obtain_NCEP
!===========================

END MODULE ObtainNCEP

! $Log$
! Revision 2.0  2000/09/05 18:57:06  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.3  2000/06/19 22:46:27  lungu
! Added code so that it reads one pressure level at a time.
! Changed vname so that it corresponds to the one in file.
!

