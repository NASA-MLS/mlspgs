! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module L2AUXData                 ! Data types for storing L2AUX data internally
!=============================================================================

  use Allocate_Deallocate, only: Allocate_test, Deallocate_test
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error, MLSMSG_Warning
  use MLSCommon, only: R8

  implicit none

  !---------------------------- RCS Ident Info -------------------------------
  character(len=256), private :: Id = &
    & "$Id$"
  character(len=*), private, parameter :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------


  ! This module defines datatypes and gives basic routines for storing and
  ! manipulating L2AUX data.

  integer, parameter :: L2AUXNameLen = 80
  integer, parameter :: CCSDSLen = 27

  ! This is a set of possible values for dimension%dimensionFamily

  integer, parameter :: NoL2AUXDimTypes = 7
  character(len=25), dimension(NoL2AUXDimTypes), parameter :: &
    & L2AUXDimNames= (/ &
       & "MLS Channel             ", &
       & "Intermediate Frequency  ", &
       & "Upper sideband Frequency", &
       & "Lower sideband Frequency", &
       & "Minor Frame             ", &
       & "Major Frame             ", &
       & "Geodetic Angle          " /)
  character (len=10), dimension(NoL2AUXDimTypes), parameter :: &
    & L2AUXDimUnits= (/ &
       & "          ", &
       & "MHz       ", &
       & "MHz       ", &
       & "MHz       ", &
       & "          ", &
       & "          ", &
       & "Degrees   "/)

  integer, parameter :: L2AUXDim_None = 0
  integer, parameter :: L2AUXDim_Channel = 1
  integer, parameter :: L2AUXDim_IntermediateFrequency = 2
  integer, parameter :: L2AUXDim_USBFrequency = 3
  integer, parameter :: L2AUXDim_LSBFrequency = 4
  integer, parameter :: L2AUXDim_MIF = 5
  integer, parameter :: L2AUXDim_MAF = 6
  integer, parameter :: L2AUXDim_GeodAngle = 7

  ! This datatype describes a dimension for an L2AUX quantity

  type L2AUX_Dimension_T
     integer :: noValues        ! Length of this dimension
     integer :: dimensionFamily ! What is this dimension
     real(r8), dimension(:), pointer :: values ! (noValues)
  end type L2AUX_Dimension_T

  ! This datatype describes an l2aux quantity itself.
  ! The dimensions will typically be ordered as follows:
  ! [Channel or frequency], MIF, [MAF or time or geodAngle]

  type L2AUXData_T

    ! A name for the L2AUX quantity, goes into SD name
    integer :: Name              ! String index of name for quantity to
                                 ! be output
    integer :: noDimensionsUsed  ! Number of dimensions used in quantity (max 3)

    ! The dimensions for the quantity
    type (L2AUX_Dimension_T), dimension(3) :: dimensions

    real(r8), pointer, dimension(:,:,:) :: values
  end type L2AUXData_T

contains ! =====     Public Procedures     =============================
  ! ---------------------------------------  SetupNewL2AUXRecord   -----
  subroutine SetupNewL2AUXRecord ( dimensionFamilies, dimSizes, l2aux )

  ! This first routine sets up the arrays for an l2aux datatype.
  ! The user supplies a set of three dimensionFamilies (e.g. L2AUXDim_MAF)
  ! Quantities can have upto three valid dimensions.  L2AUXDim_None can be used
  ! to indicate later dimensions are invalid.

    ! Dummy arguments
    integer, dimension(3), intent(in) :: dimensionFamilies
    integer, dimension(3), intent(in) :: dimSizes
    type (L2AUXData_T), intent(out) :: l2aux

    ! Local variables
    integer :: dimIndex
    integer :: status

    ! Fill the dimensions data structure

    l2aux%noDimensionsUsed = COUNT(dimensionFamilies /= L2AUXDim_None)
    l2aux%dimensions%dimensionFamily = dimensionFamilies
    l2aux%dimensions%noValues = dimSizes

    ! Allocate the values for each dimension
    do dimIndex = 1, 3
      if ( dimensionFamilies(dimIndex)/=L2AUXDim_None ) then
        allocate ( l2aux%dimensions(dimIndex)%values(dimSizes(dimIndex)), &
          & STAT=status)
        if ( status/=0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
          & MLSMSG_Allocate // "l2aux dimension values" )
      else
        l2aux%dimensions(dimIndex)%noValues=1
      end if
    end do

    ! Allocate the values for the data itself

    allocate ( l2aux%values( &
         & l2aux%dimensions(1)%noValues, &
         & l2aux%dimensions(2)%noValues, &
         & l2aux%dimensions(3)%noValues), STAT=status)
    if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate// "l2aux values" )

  end subroutine SetupNewL2AUXRecord
    
  !----------------------------------------  DestroyL2AUXContents  -----
  subroutine DestroyL2AUXContents ( l2aux )

  ! This routine deallocates all the arrays allocated above.

    ! Dummy arguments
    type (L2AUXData_T), intent(inout) :: l2aux
    ! Local variables
    integer :: status
    ! Executable code
    deallocate ( l2aux%dimensions(1)%values, &
         & l2aux%dimensions(2)%values, &
         & l2aux%dimensions(3)%values, stat=status)
    if (status /= 0) call MLSMessage ( MLSMSG_Warning, ModuleName, &
         & MLSMSG_DeAllocate // "l2aux%dimensions" ) 
    deallocate ( l2aux%values, stat=status)
    if (status /= 0) call MLSMessage ( MLSMSG_Warning, ModuleName, &
         & MLSMSG_DeAllocate // "l2aux%values" ) 
  end subroutine DestroyL2AUXContents

  !--------------------------------------  ExpandL2AUXDataInPlace  -----
  subroutine ExpandL2AUXDataInPlace ( l2aux, newSize )

  ! This subroutine expands an L2AUXData_T in place, allowing the user to
  ! add more `profiles' to it.  Note that the `profile' dimension is the last
  ! one.

    ! Dummy arguments
    type (L2AUXData_T), intent(inout) :: l2aux
    integer, intent(in) :: newSize

    ! Local variables
    integer :: status           ! From ALLOCATE
    ! The following are temporary arrays for copying data around
    real (r8), dimension(:), pointer :: temp1D
    real (r8), dimension(:,:,:), pointer :: temp3D
    integer :: expandingDimension
    integer :: oldSize

    ! Executable code

    ! First identity which is the `last' dimension, this is the one to expand

    expandingDimension = l2aux%noDimensionsUsed

    ! Now see how long this is
    oldSize = l2aux%dimensions(expandingDimension)%noValues
    ! Do a sanity check
    if ( newSize<oldSize ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "This l2aux is getting smaller not bigger" )

    ! Now expand this dimension
    temp1D => l2aux%dimensions(expandingDimension)%values
    l2aux%dimensions(expandingDimension)%noValues = newSize
    allocate ( l2aux%dimensions(expandingDimension)%values(newSize), &
      & STAT=status )
    if (status/=0) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "New dimension information")
    l2aux%dimensions(expandingDimension)%values(1:oldSize)=temp1D
    call deallocate_test ( temp1D, "temp1D", ModuleName )

    ! Now expand the data in this dimension
    temp3d => l2aux%values
    select case (expandingDimension)
    case (1)
      call allocate_test ( l2aux%values, newSize, 1, 1, "l2aux%values", &
        & ModuleName )
      l2aux%values(1:oldSize,:,:) = temp3d
    case(2)
      call allocate_test ( l2aux%values, &
        & l2aux%dimensions(1)%noValues, newSize, 1, "l2aux%values", ModuleName )
      l2aux%values(:,1:oldSize,:) = temp3d
    case(3)
      call allocate_test ( l2aux%values, &
        & l2aux%dimensions(1)%noValues, l2aux%dimensions(2)%noValues, &
        & newSize, "l2aux%values", ModuleName )
      l2aux%values(:,:,1:oldSize) = temp3d
    end select
    call deallocate_test ( temp3d, "temp3d", ModuleName )

    ! That's it.

  end subroutine ExpandL2AUXDataInPlace

  !------------------------------------------  AddL2AUXToDatabase  -----
  integer function AddL2AUXToDatabase ( DATABASE, ITEM )

  ! This subroutine adds an l2aux data type to a database of said types,
  ! creating a new database if it doesn't exist.  The result value is
  ! the length of the database -- where L2aux is put.

    ! Dummy arguments
    type (L2AUXData_T), dimension(:), pointer :: DATABASE
    type (L2AUXData_T), intent(in) :: ITEM

    ! local variables
    type (L2AUXData_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    database(newSize) = item
    AddL2AUXToDatabase = newSize
  end function AddL2AUXToDatabase

  ! ---------------------------------------  DestroyL2AUXDatabase  -----
  SUBROUTINE DestroyL2AUXDatabase ( DATABASE )

  ! This subroutine destroys a quantity template database

    ! Dummy argument
    type (L2AUXData_T), dimension(:), pointer :: DATABASE

    ! Local variables
    integer :: l2auxIndex, status

    if ( associated(database) ) then
      do l2auxIndex = 1, SIZE(database)
        call DestroyL2AUXContents ( database(l2auxIndex) )
      end do
      deallocate ( database, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Warning,ModuleName, &
        & MLSMSG_DeAllocate // "database")
    end if
  end subroutine DestroyL2AUXDatabase


!=============================================================================
end module L2AUXData
!=============================================================================

!
! $Log$
! Revision 2.0  2000/09/05 18:57:02  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry
!

