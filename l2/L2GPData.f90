! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module L2GPData                 ! Data types for storing L2GP data internally
!=============================================================================

  use Allocate_Deallocate, only: Allocate_test, Deallocate_test
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error, MLSMSG_Warning
  use MLSCommon, only: R8

  implicit none

  !---------------------------- RCS Ident Info -------------------------------
  character(len=256), private :: Id = &
    & "$Id$"
  character(len=*), parameter, private :: ModuleName = &
    & "$RCSfile$"
  !---------------------------------------------------------------------------


  ! This module defines datatypes and gives basic routines for storing and
  ! manipulating L2GP data.

  integer, parameter :: L2GPNameLen = 80
  integer, parameter, private :: CCSDSLen = 32

  ! This datatype is the main one, it simply defines one l2gp swath

  type L2GPData_T

    integer :: Name              ! String index of name for quantity to
                                 ! be output
 
    integer :: noInstances       ! Total number of horizontal instances
                                 ! (either scans or profiles)
    integer :: noSurfs           ! Total number of surfaces
    integer :: noFreqs           ! Number of frequencies in breakdown

    ! Now we store the geolocation fields, first the vertical one:
    real (r8), pointer, dimension(:) :: pressures ! Vertical coords (noSurfs)

    ! Now the horizontal geolocation information. Dimensioned (noInstances)
    real (r8), pointer, dimension(:) :: latitude, longitude, solarTime, &
      & solarZenith, losAngle, geodAngle
    real (r8), pointer, dimension(:) :: time
    integer, pointer, dimension(:) :: chunkNumber
    character (len=CCSDSLen), pointer, dimension(:) :: ccsdsTime

    ! Now we store the `frequency' geolocation field

    real (r8), pointer, dimension(:) :: frequency
    !        dimensioned (noFreqs)

    ! Finally we store the data fields

    real (r8), pointer, dimension(:,:,:) :: l2gpValue
    real (r8), pointer, dimension(:,:,:) :: l2gpPrecision
    ! dimensioned (noFreqs, noSurfs, noInstances)

    character (len=1), pointer, dimension(:) :: l2gpStatus
    !                (status is a reserved word in F90)
    real (r8), pointer, dimension(:) :: quality
    ! Both the above dimensioned (noInstances)

  end type L2GPData_T

contains ! =====     Public Procedures     =============================

  !------------------------------------------  SetupNewL2GPRecord  -----
  subroutine SetupNewL2GPRecord ( noSurfs, noFreqs, l2gp )

  ! This routine sets up the arrays for an l2gp datatype.

    ! Dummy arguments
    integer, intent(in) :: noSurfs, noFreqs ! Dimensions
    type (L2GPData_T), intent(out)  :: l2gp

    ! Local variables
    integer :: freqsArrayLen, noInstances, status, surfsArrayLen

    call allocate_test ( l2gp%pressures, noSurfs, "l2gp%pressures", &
      & ModuleName )
    surfsArrayLen = max(noSurfs,1)

    call allocate_test ( l2gp%frequency, noFreqs, "l2gp%frequency", ModuleName )
    freqsArrayLen = max(1,noFreqs)

    noInstances = 0            ! Initially empty, so the following
                               ! allocates don't actually do anything
    allocate ( l2gp%latitude(noInstances), l2gp%longitude(noInstances), &
      & l2gp%solarTime(noInstances), l2gp%solarZenith(noInstances), &
      & l2gp%losAngle(noInstances), l2gp%geodAngle(noInstances), &
      & l2gp%chunkNumber(noInstances), l2gp%ccsdsTime(noInstances), &
      & l2gp%time(noInstances), STAT=status )
    if ( status /=0 ) call MLSMessage ( MLSMSG_error, ModuleName, &
      & MLSMSG_Allocate // "l2gp horizontal coordinates" )

    allocate ( l2gp%l2gpValue(freqsArrayLen,surfsArrayLen,noInstances), &
      & l2gp%l2gpPrecision(freqsArrayLen,surfsArrayLen,noInstances), &
      & STAT=status )
    if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "l2gp data/precision fields" )

    allocate ( l2gp%l2gpStatus(noInstances),l2gp%quality(noInstances), STAT=status )
    if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "l2gp status/quality fields" )

    l2gp%noInstances = noInstances
    l2gp%noSurfs = noSurfs
    l2gp%noFreqs = noFreqs
  end subroutine SetupNewL2GPRecord
    
  !-----------------------------------------  DestroyL2GPContents  -----
  SUBROUTINE DestroyL2GPContents ( L2GP )

  ! This routine deallocates all the arrays allocated above.

    ! Dummy arguments
    type (L2GPData_T), intent(inout) :: L2GP
    ! Local variables

    integer status

    ! Executable code

    call deallocate_test ( l2gp%pressures, "l2gp%pressures", ModuleName )
    call deallocate_test ( l2gp%latitude, "l2gp%latitude", ModuleName )
    call deallocate_test ( l2gp%longitude, "l2gp%longitude", ModuleName )
    call deallocate_test ( l2gp%solarTime, "l2gp%solarTime", ModuleName )
    call deallocate_test ( l2gp%solarZenith, "l2gp%solarZenith", ModuleName )
    call deallocate_test ( l2gp%losAngle, "l2gp%losAngle", ModuleName )
    call deallocate_test ( l2gp%losAngle, "l2gp%losAngle", ModuleName )
    call deallocate_test ( l2gp%geodAngle, "l2gp%geodAngle", ModuleName )
    call deallocate_test ( l2gp%chunkNumber, "l2gp%chunkNumber", ModuleName )
    call deallocate_test ( l2gp%ccsdsTime, "l2gp%ccsdsTime", ModuleName )
    call deallocate_test ( l2gp%time, "l2gp%time", ModuleName )
    call deallocate_test ( l2gp%frequency, "l2gp%frequency", ModuleName )
    call deallocate_test ( l2gp%l2gpValue, "l2gp%l2gpValue", ModuleName )
    call deallocate_test ( l2gp%l2gpPrecision, "l2gp%l2gpPrecision", ModuleName )
    call deallocate_test ( l2gp%l2gpStatus, "l2gp%l2gpStatus", ModuleName )
    call deallocate_test ( l2gp%quality, "l2gp%quality", ModuleName )
    l2gp%noInstances = 0
    l2gp%noSurfs = 0
    l2gp%noFreqs = 0
  end subroutine DestroyL2GPContents

  !---------------------------------------  ExpandL2GPDataInPlace  -----
  subroutine ExpandL2GPDataInPlace ( l2gp, newNoInstances )

  ! This subroutine expands an L2GPData_T in place allowing the user to add
  ! more profiles to it.

    ! Dummy arguments
    type (L2GPData_T), intent(inout) :: l2gp
    integer, intent(in) :: newNoInstances

    ! Local variables
    integer :: status                   ! From ALLOCATE
    type (L2GPData_T) :: tempL2gp       ! For copying data around

    ! Executable code
    ! First do a sanity check

    if ( newNoInstances<l2gp%noInstances ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "The number of profiles requested is fewer than those already present" )

    tempL2gp = l2gp ! to get copies of all of the pointers

    ! Now go through the parameters one by one, allocate new space, copy
    ! the previous contents, and remove the old information.

    call allocate_test ( l2gp%latitude, newNoInstances, "l2gp%latitude", ModuleName )
    l2gp%latitude(1:l2gp%noInstances) = templ2gp%latitude
    call deallocate_test ( templ2gp%latitude, "templ2gp%latitude", ModuleName )
    call allocate_test ( l2gp%longitude, newNoInstances, "l2gp%longitude", ModuleName )
    l2gp%longitude(1:l2gp%noInstances) = templ2gp%longitude
    call deallocate_test ( templ2gp%longitude, "templ2gp%longitude", ModuleName )
    call allocate_test ( l2gp%solarTime, newNoInstances, "l2gp%solarTime", ModuleName )
    l2gp%solarTime(1:l2gp%noInstances) = templ2gp%solarTime
    call deallocate_test ( templ2gp%solarTime, "templ2gp%solarTime", ModuleName )
    call allocate_test ( l2gp%solarZenith, newNoInstances, "l2gp%solarZenith", ModuleName )
    l2gp%solarZenith(1:l2gp%noInstances) = templ2gp%solarZenith
    call deallocate_test ( templ2gp%solarZenith, "templ2gp%solarZenith", ModuleName )
    call allocate_test ( l2gp%losAngle, newNoInstances, "l2gp%losAngle", ModuleName )
    l2gp%losAngle(1:l2gp%noInstances) = templ2gp%losAngle
    call deallocate_test ( templ2gp%losAngle, "templ2gp%losAngle", ModuleName )
    call allocate_test ( l2gp%geodAngle, newNoInstances, "l2gp%geodAngle", ModuleName )
    l2gp%geodAngle(1:l2gp%noInstances) = templ2gp%geodAngle
    call deallocate_test ( templ2gp%geodAngle, "templ2gp%geodAngle", ModuleName )
    call allocate_test ( l2gp%time, newNoInstances, "l2gp%time", ModuleName )
    l2gp%time(1:l2gp%noInstances) = templ2gp%time
    call deallocate_test ( templ2gp%time, "templ2gp%time", ModuleName )

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    call allocate_test ( l2gp%chunkNumber, newNoInstances, "l2gp%chunkNumber", ModuleName )
    l2gp%chunkNumber(1:l2gp%noInstances) = templ2gp%chunkNumber
    call deallocate_test ( templ2gp%chunkNumber, "templ2gp%chunkNumber", ModuleName )
    call allocate_test ( l2gp%ccsdsTime, newNoInstances, "l2gp%ccsdsTime", ModuleName )
    l2gp%ccsdsTime(1:l2gp%noInstances) = templ2gp%ccsdsTime
    call deallocate_test ( templ2gp%ccsdsTime, "templ2gp%ccsdsTime", ModuleName )

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    call allocate_test ( l2gp%l2gpValue, &
      & l2gp%noFreqs, l2gp%noSurfs, newNoInstances, "l2gp%l2gpValue", ModuleName )
    l2gp%l2gpValue(:,:,1:l2gp%noInstances) = templ2gp%l2gpValue
    call deallocate_test ( templ2gp%l2gpValue, "templ2gp%l2gpValue", ModuleName )
    call allocate_test ( l2gp%l2gpPrecision,&
      & l2gp%noFreqs, l2gp%noSurfs, newNoInstances, "l2gp%l2gpValue", ModuleName )
    l2gp%l2gpPrecision(:,:,1:l2gp%noInstances) = templ2gp%l2gpPrecision
    call deallocate_test ( templ2gp%l2gpPrecision, "templ2gp%l2gpPrecision", ModuleName )

    call allocate_test ( l2gp%l2gpStatus, newNoInstances, "l2gp%l2gpStatus", &
      & ModuleName )
    l2gp%l2gpStatus(1:l2gp%noInstances) = templ2gp%l2gpStatus
    call deallocate_test ( templ2gp%l2gpStatus, "templ2gp%l2gpStatus", ModuleName )

    call allocate_test ( l2gp%quality, newNoInstances, "l2gp%quality", ModuleName )
    l2gp%quality(1:l2gp%noInstances) = templ2gp%quality
    call deallocate_test ( templ2gp%quality, "templ2gp%quality", Modulename )
    l2gp%noInstances=newNoInstances

  end subroutine ExpandL2GPDataInPlace

  !-------------------------------------------  AddL2GPToDatabase  -----
  integer function AddL2GPToDatabase( DATABASE, ITEM )

  ! This function adds an l2gp data type to a database of said types,
  ! creating a new database if it doesn't exist.  The result value is
  ! the size -- where L2gp is put.

    ! Dummy arguments
    type (l2gpdata_t), dimension(:), pointer :: DATABASE
    type (l2gpdata_t), intent(in) :: ITEM

    ! Local variables
    type (L2GPData_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    database(newSize) = item
    AddL2GPToDatabase = newSize
  end function AddL2GPToDatabase

  ! --------------------------------------------------------------------------

  ! This subroutine destroys a quantity template database

  subroutine DestroyL2GPDatabase ( DATABASE )

    ! Dummy argument
    type (L2GPData_T), dimension(:), pointer :: DATABASE

    ! Local variables
    integer :: l2gpIndex, status

    if ( associated(database)) then
      do l2gpIndex = 1, SIZE(database)
        call DestroyL2GPContents ( database(l2gpIndex) )
      end do
      deallocate ( database, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
         & MLSMSG_deallocate // "database" )
    end if
  end subroutine DestroyL2GPDatabase


!=============================================================================
end module L2GPData
!=============================================================================

!
! $Log$
! Revision 2.0  2000/09/05 18:57:03  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry
!

