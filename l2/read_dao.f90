! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===========================================================

MODULE ObtainDAO !provides subroutines to access DAO files

!===========================================================
USE MLSCommon
USE MLSStrings 
USE MLSMessageModule
USE GriddedData
USE VerticalCoordinate
USE MLSPCF
USE MLSCF
USE dates_module    

IMPLICIT NONE

PUBLIC

PRIVATE :: Id, ModuleName
!------------------------------- RCS Ident Info ------------------------------
CHARACTER(LEN=130) :: id = & 
   "$Id$"
CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile $"
!-----------------------------------------------------------------------------

!Parameters
  INTEGER, PARAMETER :: XDIM=360
  INTEGER, PARAMETER :: YDIM=181
  INTEGER, PARAMETER :: ZDIM=36

CONTAINS

!===============================================
  SUBROUTINE read_dao (fname, vname, data_array)
!===============================================

! Brief description of program
! This subroutine reads a DAO correlative file and returns
! the data_array to the caller


! Arguments

  CHARACTER*(*), INTENT(IN) :: fname, vname
  REAL*4 data_array(XDIM,YDIM, ZDIM)
! - - - local declarations - - -

  INTEGER*4 sd_id, sds_id, status
  INTEGER*4 sds_index
  INTEGER*4 start(4), edges(4), stride(4)
  INTEGER i,j,k
  INTEGER file_id, gd_id
  CHARACTER (LEN=256) ::  mnemonic, msg
  INTEGER :: sdstart, gdopen, gdattach, gdrdfld, gddetach, &
                       & gdclose, gdinqdims, gdinqflds
  integer strbufsize, ndim, nfld, dim(27),numbertype(27), rank(27)
  CHARACTER (LEN=1000) :: dimname, fieldlist


! - - - begin - - -
!  sd_id = sdstart(fname, DFACC_RDONLY)

  file_id = gdopen(fname, DFACC_RDONLY)

  IF (file_id < 0) THEN
    CALL Pgs_smf_getMsg(status, mnemonic, msg)
    CALL MLSMessage (MLSMSG_Error, ModuleName, &
              "Could not open "// fname//" "//mnemonic//" "//msg)
 
  END IF

  gd_id = gdattach(file_id,"EOSGRID")
  IF (gd_id < 0) THEN
    CALL Pgs_smf_getMsg(status, mnemonic, msg)
    CALL MLSMessage (MLSMSG_Error, ModuleName, &
               "Could not attach "//fname//" "//mnemonic//" "//msg)

  END IF
  ndim = gdinqdims(gd_id , dimname, dim)
  start(1) = 0
  start(2) = 0
  start(3) = 0
  start(4) = 0

  stride(1) = 1
  stride(2) = 1
  stride(3) = 1
  stride(4) = 1

  edges(1) = 1
  do i = 1, ndim-1
     edges(1+i) = dim(ndim+1-i) 
  end do
  edges(2) = ZDIM
 ! edges(2) = YDIM
 ! edges(3) = XDIM


!   In this subroutine, we read the entire field.  By manipulating the start 
!   and edges arrays, it is possible to read a subset of the entire array.  
!   For example, to read a 3D section defined by x=100,224 y=50,149 
!   z=15,16 you would set the start and edges arrays to the following:

!   start(0) = 0    time start location
!   start(1) = 15   z-dim start location
!   start(2) = 50   y-dim start location
!   start(3) = 100  x-dim start location

!   edges(0) = 1    time length
!   edges(1) = 2    z-dim length
!   edges(2) = 100  y-dim length
!   edges(3) = 125  x-dim length



!  status = gdrdfld(gd_id, vname, start, stride, edges, data_array)
   status = PGS_S_SUCCESS
  IF(status /= PGS_S_SUCCESS) THEN
     CALL Pgs_smf_getMsg(status, mnemonic, msg)
     CALL MLSMessage (MLSMSG_Error, ModuleName, &
                      "Error reading grid field:  "//mnemonic//" "//msg)
  END IF

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
                      "Error closing  DAO File:  "//mnemonic//" "//msg)
  END IF

  RETURN
  END SUBROUTINE read_dao

  SUBROUTINE Obtain_DAO(aprioriData, l2cf_info)
 
 
    TYPE (GriddedData_T), DIMENSION(:), POINTER :: aprioriData 
    ! Input a priori database
    TYPE (MLSCF_T) :: l2cf_info        ! The information from the l2cf file

! Local Variables

    TYPE (GriddedData_T):: qty
    CHARACTER (LEN=256) :: msg, mnemonic
    INTEGER ::  returnStatus
    INTEGER :: DAOFileHandle, DAO_Version, sd_id
    CHARACTER (LEN=132) :: DAOphysicalFilename
    CHARACTER (LEN=80) :: vname 
    REAL*4 :: data_array(XDIM, YDIM, ZDIM)
!    ALLOCATE (data_array(XDIM, YDIM, ZDIM), stat=returnStatus)
    DAO_Version = 1
    vname = "TMPU" ! for now


! Get the DAO file name from the PCF

    DO DAOFileHandle = mlspcf_l2dao_start, mlspcf_l2dao_end

      returnStatus = Pgs_pc_getReference(DAOFileHandle, DAO_Version, DAOphysicalFilename)
      IF(returnStatus == PGS_S_SUCCESS) THEN

! Open the HDF-EOS file and read gridded data
        CALL read_dao (DAOphysicalFilename, vname, data_array)

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
!        DEALLOCATE (data_array, stat=returnStatus)
!        CALL AddGridTemplateToDatabase(aprioriData, qty)

      ELSE

        CALL Pgs_smf_getMsg(returnStatus, mnemonic, msg)
        CALL MLSMessage (MLSMSG_Error, ModuleName, &
                         "Error opening DAO file:  "//mnemonic//" "//msg)
      END IF

    END DO ! DAOFileHandle = mlspcf_l2_dao_start, mlspcf_l2_dao_end

  END SUBROUTINE Obtain_DAO


END MODULE ObtainDAO

! $Log$
! Revision 2.0  2000/09/05 18:57:06  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.3  2000/06/30 00:14:08  lungu
! Commented out call to gdrdfld, which doesn't
! work with the NAG compiler and Toolkit-5.2.6v1.00.NAG.
!

