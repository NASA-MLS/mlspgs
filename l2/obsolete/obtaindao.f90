! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===========================================================
module ObtainDAO !provides subroutines to access DAO files
!===========================================================

! use DATES_MODULE
  use HDF, only: DFACC_RDONLY
  use GriddedData, only: GriddedData_T, AddGridTemplateToDatabase
 use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use MLSPCF2, only: MLSPCF_L2DAO_END, MLSPCF_L2DAO_START
! use MLSStrings
  use SDPToolkit, only: Pgs_pc_getReference, PGS_S_SUCCESS
! use ???, only: Pgs_smf_getMsg
! use VerticalCoordinate

  implicit none
  private
  public :: OBTAIN_DAO, READ_DAO

  !------------------------------- RCS Ident Info ------------------------------
  character(len=130) :: id = & 
     "$Id$"
  character(len=*), parameter :: ModuleName="obtain_dao.f90,v"
  !-----------------------------------------------------------------------------

  !Parameters
  integer, parameter :: XDIM=360
  integer, parameter :: YDIM=181
  integer, parameter :: ZDIM=36

contains ! =====     Public Procedures     =============================

  ! -------------------------------------------------  OBTAIN_DAO  -----
  subroutine OBTAIN_DAO ( aprioriData, root )
 
    type (GriddedData_T), dimension(:), pointer :: aprioriData 
    ! Input a priori database
    integer, intent(in) :: ROOT        ! Root of L2CF abstract syntax tree

! Local Variables

    real(R8) :: data_array(XDIM, YDIM, ZDIM)
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
        call read_dao ( DAOphysicalFilename, vname, data_array )

!       Create a GriddedDtata Template and copy data_array into it

!        CALL SetupNewGridTemplate(qty, noHeights=XDIM, noLats=YDIM, noLons=ZDIM,&
!                                  noLsts=1,noSzas=1, noDates=1)
!        do i = 1, zdim
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
        returnStatus = AddGridTemplateToDatabase(aprioriData, qty)

      else

        call Pgs_smf_getMsg ( returnStatus, mnemonic, msg )
        call MLSMessage ( MLSMSG_Error, ModuleName, &
                         "Error opening DAO file:  "//mnemonic//" "//msg)
      end if

    end do ! DAOFileHandle = mlspcf_l2_dao_start, mlspcf_l2_dao_end

  end subroutine Obtain_DAO

  ! ---------------------------------------------------  READ_DAO  -----

  !=====================================================================
  subroutine READ_DAO ( FNAME, VNAME, DATA_ARRAY )
  !=====================================================================

  ! Brief description of program
  ! This subroutine reads a DAO correlative file and returns
  ! the data_array to the caller

  ! Arguments
    character*(*), intent(in) :: FNAME, VNAME
    real(R8) :: DATA_ARRAY(xdim,ydim, zdim)

  ! - - - local declarations - - -

    integer :: DIM(27)
    character (len=1000) :: DIMNAME
    integer :: EDGES(4)
    integer :: FILE_ID, GD_ID
    integer :: I
    character (len=256) :: MSG, MNEMONIC
    integer :: NDIM
!   integer :: SD_ID
    integer :: START(4)
    integer :: STATUS
    integer :: STRIDE(4)

  ! - - - external functions - -- 
    integer :: GDOPEN, GDATTACH, GDRDFLD, GDDETACH, &
                         & GDCLOSE, GDINQDIMS

  ! - - - begin - - -
  !  sd_id = sdstart(fname, DFACC_RDONLY)

    file_id = gdopen(fname, DFACC_RDONLY)

    if ( file_id < 0 ) then
      call Pgs_smf_getMsg ( status, mnemonic, msg )
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Could not open "// fname//" "//mnemonic//" "//msg )
    end if

    gd_id = gdattach(file_id,"EOSGRID")

    if (gd_id < 0) then
      call Pgs_smf_getMsg ( status, mnemonic, msg )
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Could not attach "//fname//" "//mnemonic//" "//msg )
    end if
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
  !  edges(2) = ZDIM
  !  edges(3) = YDIM
  !  edges(4) = XDIM


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



    status = gdrdfld(gd_id, vname, start, stride, edges, data_array)
    if ( status /= PGS_S_SUCCESS ) then
      call Pgs_smf_getMsg ( status, mnemonic, msg )
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Error reading grid field:  "//mnemonic//" "//msg)
    end if

    status = gddetach(gd_id)
    if ( status /= PGS_S_SUCCESS ) then
      call Pgs_smf_getMsg ( status, mnemonic, msg )
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Error detaching  grid:  "//mnemonic//" "//msg)
    end if

    status = gdclose(file_id)
    if ( status /= PGS_S_SUCCESS ) then
      call Pgs_smf_getMsg ( status, mnemonic, msg )
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Error closing  DAO File:  "//mnemonic//" "//msg)
    end if

    return
  end subroutine READ_DAO

! =====     Private Procedures     =====================================

END MODULE ObtainDAO

! $Log$
! Revision 2.4  2001/03/03 00:11:29  pwagner
! Began transformations to act like L2GPData module for Gridded data
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

