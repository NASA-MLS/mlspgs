
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
module OutputL2GP
!===============================================================================

   use Allocate_Deallocate, only: DEALLOCATE_TEST, DEALLOC_STATUS
   use Hdf, only: DFNT_CHAR8, DFNT_FLOAT32, DFNT_INT32
   use HDFEOS, only: SWATTACH, SWCREATE, SWDEFDFLD, SWDEFDIM, SWDEFGFLD, &
     & SWDETACH
   use L2GPData, only: L2GPData_T
   use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
   use STRING_TABLE, only: GET_STRING
   use SWAPI, only: SWWRFLD
   implicit none
   public

!------------------- RCS Ident Info -----------------------
   character(len=130), private :: Id = &
     & "$Id$"
   character (len=*), parameter, private :: ModuleName= &
     & "$RCSfile$"
!----------------------------------------------------------

! Contents:

! Subroutines -- OutputL2GP_createFile
!                OutputL2GP_writeGeo
!                OutputL2GP_writeData
!                DeallocateL2GP

! Remarks:  This module contains parameters and subroutines used in
!           producing an L2GP output file.

! Parameters

   character (len=*), parameter :: DATA_FIELD1 = 'l2gpValue'
   character (len=*), parameter :: DATA_FIELD2 = 'l2gpPrecision'
   character (len=*), parameter :: DATA_FIELD3 = 'l2gpStatus'
   character (len=*), parameter :: DATA_FIELD4 = 'quality'

   character (len=*), parameter :: GEO_FIELD1 = 'latitude'
   character (len=*), parameter :: GEO_FIELD2 = 'longitude'
   character (len=*), parameter :: GEO_FIELD3 = 'time'
   character (len=*), parameter :: GEO_FIELD4 = 'ccsdsTime'
   character (len=*), parameter :: GEO_FIELD5 = 'solarTime'
   character (len=*), parameter :: GEO_FIELD6 = 'solarZenith'
   character (len=*), parameter :: GEO_FIELD7 = 'losAngle'
   character (len=*), parameter :: GEO_FIELD8 = 'geodAngle'
   character (len=*), parameter :: GEO_FIELD9 = 'chunkNumber'
   character (len=*), parameter :: GEO_FIELD10 = 'pressures'
   character (len=*), parameter :: GEO_FIELD11 = 'frequency'

   integer, parameter :: CCSDS_LEN = 27         ! len of CCSDS time string
   integer, parameter :: HDFE_AUTOMERGE = 1     ! merge fields with share dim
   integer, parameter :: HDFE_NOMERGE = 0       ! don't merge

contains ! =====     Public Procedures     =============================

  ! --------------------------------------  OutputL2GP_createFile  -----
  subroutine OutputL2GP_createFile (L2FileHandle, l2gp, flag)

  ! Brief description of subroutine
  ! This subroutine sets up the structural definitions in an empty L2GP file.

  ! Arguments

    integer, intent(in) :: L2FileHandle

    type( L2GPData_T ), intent(inout) :: l2gp

    integer, intent(out) :: flag

  ! Parameters

    character (len=*), parameter :: DIM_NAME1 = 'noInstances'
    character (len=*), parameter :: DIM_NAME2 = 'noSurfs'
    character (len=*), parameter :: DIM_NAME3 = 'noFreqs'
    character (len=*), parameter :: DIM_NAME4 = 'CCSDSLen,noInstances'
    character (len=*), parameter :: DIM_NAME12 = 'noSurfs,noInstances'
    character (len=*), parameter :: DIM_NAME123 = 'noFreqs,noSurfs,noInstances'
    character (len=*), parameter :: DIM_NAME41 = 'CCSDSLen,noInstances'

    character (len=*), parameter :: DIM_ERR = 'Failed to define dimension '
    character (len=*), parameter :: GEO_ERR = &
      & 'Failed to define geolocation field '
    character (len=*), parameter :: DAT_ERR = 'Failed to define data field '

  ! Variables

    character (len=480) :: MSR
    character (len=132) :: NAME   ! From l2gp%name

    integer :: SWID, STATUS

    flag = 0

  ! Create the swath within the file

    call get_string ( l2gp%name, name )
    swid = swcreate(L2FileHandle, name)
    if ( swid == -1 ) then
      msr = 'Failed to create swath ' // name
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

 ! Define dimensions

    status = swdefdim(swid, DIM_NAME1, l2gp%noInstances)
    if ( status == -1 ) then
      msr = DIM_ERR // DIM_NAME1
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    if ( l2gp%noSurfs > 0 ) then
      status = swdefdim(swid, DIM_NAME2, l2gp%noSurfs)
      if ( status == -1 ) then
        msr = DIM_ERR // DIM_NAME2
        call MLSMessage ( MLSMSG_Error, ModuleName, msr )
      end if
    end if

    if ( l2gp%noFreqs > 0 ) then
      status = swdefdim(swid, DIM_NAME3, l2gp%noFreqs)
      if ( status == -1 ) then
        msr = DIM_ERR // DIM_NAME3
        call MLSMessage ( MLSMSG_Error, ModuleName, msr )
      end if
    end if

    status = swdefdim(swid, DIM_NAME4, l2gp%noInstances)
    if ( status == -1 ) then
      msr = DIM_ERR // DIM_NAME4
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

  ! Define horizontal geolocation fields using above dimensions

    status = swdefgfld(swid, GEO_FIELD1, DIM_NAME1, DFNT_FLOAT32, &
                       HDFE_NOMERGE)
    if ( status == -1 ) then
      msr = GEO_ERR // GEO_FIELD1
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swdefgfld(swid, GEO_FIELD2, DIM_NAME1, DFNT_FLOAT32, &
                       HDFE_NOMERGE)
    if ( status == -1 ) then
      msr = GEO_ERR // GEO_FIELD2
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swdefgfld(swid, GEO_FIELD3, DIM_NAME1, DFNT_FLOAT32, &
                       HDFE_NOMERGE)
    if ( status == -1 ) then
      msr = GEO_ERR // GEO_FIELD3
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swdefgfld(swid, GEO_FIELD4, DIM_NAME1, DFNT_CHAR8, &
                       HDFE_NOMERGE)
    if ( status == -1 ) then
      msr = GEO_ERR // GEO_FIELD4
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swdefgfld(swid, GEO_FIELD5, DIM_NAME1, DFNT_FLOAT32, &
                       HDFE_NOMERGE)
    if ( status == -1 ) then
      msr = GEO_ERR // GEO_FIELD5
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swdefgfld(swid, GEO_FIELD6, DIM_NAME1, DFNT_FLOAT32, &
                       HDFE_NOMERGE)
    if ( status == -1 ) then
      msr = GEO_ERR // GEO_FIELD6
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swdefgfld(swid, GEO_FIELD7, DIM_NAME1, DFNT_FLOAT32, &
                       HDFE_NOMERGE)
    if ( status == -1 ) then
      msr = GEO_ERR // GEO_FIELD7
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swdefgfld(swid, GEO_FIELD8, DIM_NAME1, DFNT_FLOAT32, &
                       HDFE_NOMERGE)
    if ( status == -1 ) then
      msr = GEO_ERR // GEO_FIELD8
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swdefgfld(swid, GEO_FIELD9, DIM_NAME1, DFNT_INT32, HDFE_NOMERGE)
    if ( status == -1 ) then
      msr = GEO_ERR // GEO_FIELD9
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    if ( l2gp%noSurfs > 0 ) then
      status = swdefgfld(swid, GEO_FIELD10, DIM_NAME2, DFNT_FLOAT32, &
                         HDFE_NOMERGE)
      if ( status == -1 ) then
        msr = GEO_ERR // GEO_FIELD10
        call MLSMessage ( MLSMSG_Error, ModuleName, msr )
      end if
    end if

    if ( l2gp%noFreqs > 0 ) then
      status = swdefgfld(swid, GEO_FIELD11, DIM_NAME3, DFNT_FLOAT32, &
                         HDFE_NOMERGE)
      if ( status == -1 ) then
        msr = GEO_ERR // GEO_FIELD11
        call MLSMessage ( MLSMSG_Error, ModuleName, msr )
      end if
    end if

  ! Define data fields using above dimensions

    if ( (l2gp%noFreqs > 0) .and. (l2gp%noSurfs > 0) ) then

      status = swdefdfld(swid, DATA_FIELD1, DIM_NAME123, DFNT_FLOAT32, &
                         HDFE_NOMERGE)

      if ( status == -1 ) then
        msr = DAT_ERR // DATA_FIELD1 // ' for 3D quantity.'
        call MLSMessage ( MLSMSG_Error, ModuleName, msr )
      end if


      status = swdefdfld(swid, DATA_FIELD2, DIM_NAME123, DFNT_FLOAT32, &
                         HDFE_NOMERGE)

      if ( status == -1 ) then
        msr = DAT_ERR // DATA_FIELD2 // ' for 3D quantity.'
        call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if


    else if ( l2gp%noSurfs > 0 ) then

       status = swdefdfld(swid, DATA_FIELD1, DIM_NAME12, DFNT_FLOAT32, &
                          HDFE_NOMERGE)

       if ( status == -1 ) then
         msr = DAT_ERR // DATA_FIELD1 //  ' for 2D quantity.'
         call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if

       status = swdefdfld(swid, DATA_FIELD2, DIM_NAME12, DFNT_FLOAT32, &
                          HDFE_NOMERGE)

       if ( status == -1 ) then
         msr = DAT_ERR // DATA_FIELD2 //  ' for 2D quantity.'
         call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if

    else

       status = swdefdfld(swid, DATA_FIELD1, DIM_NAME1, DFNT_FLOAT32, &
                          HDFE_NOMERGE)

       if ( status == -1 ) then
         msr = DAT_ERR // DATA_FIELD1 // ' for 1D quantity.'
         call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if

       status = swdefdfld(swid, DATA_FIELD2, DIM_NAME1, DFNT_FLOAT32, &
                          HDFE_NOMERGE)

       if ( status == -1 ) then
         msr = DAT_ERR // DATA_FIELD2 // ' for 1D quantity.'
         call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if

    end if

    status = swdefdfld(swid, DATA_FIELD3, DIM_NAME1, DFNT_CHAR8, &
                       HDFE_NOMERGE)
    if ( status == -1 ) then
      msr = DAT_ERR // DATA_FIELD3
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swdefdfld(swid, DATA_FIELD4, DIM_NAME1, DFNT_FLOAT32, &
                       HDFE_NOMERGE)
    if ( status == -1 ) then
      msr = DAT_ERR // DATA_FIELD4
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

  ! Detach from the swath interface.  This stores the swath info within the
  ! file and must be done before writing or reading data to or from the
  ! swath.

    status = swdetach(swid)
    if ( status == -1 ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Failed to detach from swath interface after definition.' )
      flag = -1

    end if

  !--------------------------------------
  end subroutine OutputL2GP_createFile
  !--------------------------------------

  !-----------------------------------------  OutputL2GP_writeGeo  -----
  subroutine OutputL2GP_writeGeo (l2gpGeo, swfid)

  ! Brief description of subroutine
  ! This subroutine writes the geolocation fields to an L2GP output file.

  ! Arguments

    type( L2GPData_T ), intent(inout) :: l2gpGeo

    integer, intent(in) :: swfid

  ! Parameters

    character (len=*), parameter :: WR_ERR = &
      & 'Failed to write geolocation field '

  ! Variables

    character (len=480) :: msr
    character (len=132) :: NAME    ! From l2gpGeo%name

    integer :: status, swid
    integer :: start(2), stride(2), edge(2)

    call get_string ( l2gpGeo%name, name )
    swid = swattach (swfid, name)

  ! Write data to the fields
      
    stride(1) = 1
    start(1) = 0
    edge(1) = l2gpGeo%noInstances

    status = swwrfld(swid, GEO_FIELD1, start, stride, edge, &
                     REAL(l2gpGeo%latitude))
    if ( status == -1 ) then
      msr = WR_ERR // GEO_FIELD1
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if
    status = swwrfld(swid, GEO_FIELD2, start, stride, edge, &
                     REAL(l2gpGeo%longitude))

    if ( status == -1 ) then
      msr = WR_ERR // GEO_FIELD2
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if
    status = swwrfld(swid, GEO_FIELD3, start, stride, edge, &
                     REAL(l2gpGeo%time))

    if ( status == -1 ) then
      msr = WR_ERR // GEO_FIELD3
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if
    status = swwrfld( swid, GEO_FIELD4, start, stride, edge, &
                      l2gpGeo%ccsdsTime)

    if ( status == -1 ) then
      msr = WR_ERR // GEO_FIELD4
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if
    status = swwrfld(swid, GEO_FIELD5, start, stride, edge, &
                     REAL(l2gpGeo%solarTime))

    if ( status == -1 ) then
      msr = WR_ERR // GEO_FIELD5
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if
    status = swwrfld(swid, GEO_FIELD6, start, stride, edge, &
                     REAL(l2gpGeo%solarZenith))

    if ( status == -1 ) then
      msr = WR_ERR // GEO_FIELD6
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if
    status = swwrfld(swid, GEO_FIELD7, start, stride, edge, &
                     REAL(l2gpGeo%losAngle))

    if ( status == -1 ) then
      msr = WR_ERR // GEO_FIELD7
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if
    status = swwrfld(swid, GEO_FIELD8, start, stride, edge, &
                     REAL(l2gpGeo%geodAngle))

    if ( status == -1 ) then
      msr = WR_ERR // GEO_FIELD8
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if
    status = swwrfld(swid, GEO_FIELD9, start, stride, edge, &
                     l2gpGeo%chunkNumber)

    if ( status == -1 ) then
      msr = WR_ERR // GEO_FIELD9
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if
    if ( l2gpGeo%noSurfs > 0 ) then

      edge(1) = l2gpGeo%noSurfs

      status = swwrfld(swid, GEO_FIELD10, start, stride, edge, &
                       REAL(l2gpGeo%pressures))

      if ( status == -1 ) then
        msr = WR_ERR // GEO_FIELD10
        call MLSMessage ( MLSMSG_Error, ModuleName, msr )
      end if
    end if
    if ( l2gpGeo%noFreqs > 0 ) then

      edge(1) = l2gpGeo%noFreqs
      l2gpGeo%frequency = 0
      status = swwrfld(swid, GEO_FIELD11, start, stride, edge, &
                       REAL(l2gpGeo%frequency))

      if ( status == -1 ) then
        msr = WR_ERR // GEO_FIELD11
        call MLSMessage ( MLSMSG_Error, ModuleName, msr )
      end if
    end if

  ! Detach from the swath interface.  

    status = swdetach(swid)

    if ( status == -1 ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Failed to detach from swath interface' )
    end if
         
  !------------------------------------
  end subroutine OutputL2GP_writeGeo
  !------------------------------------

  !----------------------------------------  OutputL2GP_writeData  -----
   subroutine OutputL2GP_writeData(l2gpData, swfid)

  ! Brief description of subroutine
  ! This subroutine writes the data fields to an L2GP output file.

  ! Arguments

    type( L2GPData_T ), intent(inout) :: l2gpData

    integer, intent(in) :: swfid

  ! Parameters

    character (len=*), parameter :: WR_ERR = 'Failed to write data field '

  ! Variables

    character (len=480) :: msr
   character (len=132) :: NAME     ! From l2gpData%name

    integer :: status
    integer :: start(3), stride(3), edge(3)
    integer :: swid

  ! Write data to the fields

    start = 0
    stride = 1
    edge(1) = l2gpData%noFreqs
    edge(2) = l2gpData%noSurfs
    edge(3) = l2gpData%noInstances
    call get_string ( l2gpData%name, name )
    swid = swattach (swfid, name)
    if ( l2gpData%noFreqs > 0 ) then

  ! Value and Precision are 3-D fields

      status = swwrfld(swid, DATA_FIELD1, start, stride, edge, &
        & reshape(l2gpData%l2gpValue, (/size(l2gpData%l2gpValue)/)) )
      if ( status == -1 ) then
        msr = WR_ERR // DATA_FIELD1
        call MLSMessage ( MLSMSG_Error, ModuleName, msr )
      end if
      status = swwrfld(swid, DATA_FIELD2, start, stride, edge, &
        & reshape(REAL(l2gpData%l2gpPrecision), (/size(l2gpData%l2gpPrecision)/)) )
      if ( status == -1 ) then
        msr = WR_ERR // DATA_FIELD2
        call MLSMessage ( MLSMSG_Error, ModuleName, msr )
      end if

    else if ( l2gpData%noSurfs > 0 ) then

  ! Value and Precision are 2-D fields

      status = swwrfld( swid, DATA_FIELD1, start(2:3), stride(2:3), &
                        edge(2:3), REAL(l2gpData%l2gpValue(1,:,:) ))
      if ( status == -1 ) then
        msr = WR_ERR // DATA_FIELD1
        call MLSMessage ( MLSMSG_Error, ModuleName, msr )
      end if
      status = swwrfld( swid, DATA_FIELD2, start(2:3), stride(2:3), &
                        edge(2:3), REAL(l2gpData%l2gpPrecision(1,:,:) ))
      if ( status == -1 ) then
        msr = WR_ERR // DATA_FIELD2
        call MLSMessage ( MLSMSG_Error, ModuleName, msr )
      end if
    else

  ! Value and Precision are 1-D fields

      status = swwrfld( swid, DATA_FIELD1, start(3:3), stride(3:3), edge(3:3), &
                        REAL(l2gpData%l2gpValue(1,1,:) ))
      if ( status == -1 ) then
        msr = WR_ERR // DATA_FIELD1
        call MLSMessage ( MLSMSG_Error, ModuleName, msr )
      end if
      status = swwrfld( swid, DATA_FIELD2, start(3:3), stride(3:3), edge(3:3), &
                        REAL(l2gpData%l2gpPrecision(1,1,:) ))
      if ( status == -1 ) then
        msr = WR_ERR // DATA_FIELD2
        call MLSMessage ( MLSMSG_Error, ModuleName, msr )
      end if
    end if

! 1-D status & quality fields

    status = swwrfld(swid, DATA_FIELD3, start(3:3), stride(3:3), edge(3:3), &
                     l2gpData%l2gpStatus)
    if ( status == -1 ) then
      msr = WR_ERR // DATA_FIELD3
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if
    l2gpData%quality = 0
    status = swwrfld(swid, DATA_FIELD4, start(3:3), stride(3:3), edge(3:3), &
                     REAL(l2gpData%quality))
    if ( status == -1 ) then
      msr = WR_ERR // DATA_FIELD4
      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

  !     Detach from the swath interface.

    status = swdetach(swid)
    if ( status == -1 ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Failed to detach  from swath interface' )
    end if


  !-------------------------------------
  end subroutine OutputL2GP_writeData
  !-------------------------------------

  !----------------------------------------------  DeallocateL2GP  -----
  subroutine DeallocateL2GP ( l2gp, flag )

  ! Brief description of subroutine
  ! This subroutine deallocates the internal field pointers of the L2GP_T
  ! derived type, after the calling program has finished with the data.

  ! Arguments

    type( L2GPData_T ), intent(inout) :: l2gp

    integer, intent(out) :: flag

  ! Parameters

  ! Functions

  ! Variables

    flag = 0
    dealloc_status = 0

! Horizontal geolocation fields

    if ( associated(l2gp%time) ) &
      & call deallocate_test ( l2gp%time, "l2gp%time", ModuleName )
    if ( associated(l2gp%chunkNumber) ) &
      & call deallocate_test ( l2gp%chunkNumber, "l2gp%chunkNumber", ModuleName )
    if ( associated(l2gp%latitude) ) &
      & call deallocate_test ( l2gp%latitude, "l2gp%latitude", ModuleName )
    if ( associated(l2gp%longitude) ) &
      & call deallocate_test ( l2gp%longitude, "l2gp%longitude", ModuleName )
    if ( associated(l2gp%solarTime) ) &
      & call deallocate_test ( l2gp%solarTime, "l2gp%solarTime", ModuleName )
    if ( associated(l2gp%solarZenith) ) &
      & call deallocate_test ( l2gp%solarZenith, "l2gp%solarZenith", ModuleName )
    if ( associated(l2gp%losAngle) ) &
      & call deallocate_test ( l2gp%losAngle, "l2gp%losAngle", ModuleName )
    if ( associated(l2gp%geodAngle) ) &
      & call deallocate_test ( l2gp%geodAngle, "l2gp%geodAngle", ModuleName )
    if ( associated(l2gp%ccsdsTime) ) &
      & call deallocate_test ( l2gp%ccsdsTime, "l2gp%ccsdsTime", ModuleName )

  ! Vertical geolocation field

    if ( associated(l2gp%pressures) ) &
      & call deallocate_test ( l2gp%pressures, "l2gp%pressures", ModuleName )

  ! Frequency "geolocation field"

    if ( associated(l2gp%frequency) ) &
      & call deallocate_test ( l2gp%frequency, "l2gp%frequency", ModuleName )

  ! Data fields

    if ( associated(l2gp%l2gpValue) ) &
      & call deallocate_test ( l2gp%l2gpValue, "l2gp%l2gpValue", ModuleName )
    if ( associated(l2gp%l2gpPrecision) ) &
      & call deallocate_test ( l2gp%l2gpPrecision, "l2gp%l2gpPrecision", &
      & ModuleName )
    if ( associated(l2gp%l2gpStatus) ) &
      & call deallocate_test ( l2gp%l2gpStatus, "l2gp%l2gpStatus", ModuleName )
    if ( associated(l2gp%quality) ) &
      & call deallocate_test ( l2gp%quality, "l2gp%quality", ModuleName )

    if ( dealloc_status /= 0 ) flag = -1

  !-------------------------------
  end subroutine DeallocateL2GP
!-------------------------------

!====================
end module OutputL2GP
!====================

!# $Log$
!# Revision 2.2  2000/09/13 22:44:47  ahanzel
!# Removed old log entries in file.
!#
!# Revision 2.1  2000/09/12 21:24:52  vsnyder
!# Revised to use cf parser output directly.
!#

