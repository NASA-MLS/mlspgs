! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module L2AUXData                 ! Data types for storing L2AUX data internally
!=============================================================================

  use Allocate_Deallocate, only: Allocate_test, Deallocate_test
  use Hdf, only: DFACC_READ, DFNT_FLOAT64, SFCREATE, SFDIMID, SFEND, &
    & SFENDACC, SFSTART, SFRDATA, sfn2index, sfselect, sfgetinfo, &
    & sfgdinfo
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error, MLSMSG_Warning
  use MLSCommon, only: R8
  use MLSStrings, only: LinearSearchStringArray

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
  integer, parameter :: L2AUXDimNameLen = 25
  character(len=L2AUXDimNameLen), dimension(NoL2AUXDimTypes), parameter :: &
    & L2AUXDimNames= (/ &
       & "MLS Channel             ", &
       & "Intermediate Frequency  ", &
       & "Upper sideband Frequency", &
       & "Lower sideband Frequency", &
       & "Minor Frame             ", &
       & "Major Frame             ", &
       & "Geodetic Angle          " /)
  integer, parameter :: L2AUXDimUnitLen = 10
  character (len=L2AUXDimUnitLen), dimension(NoL2AUXDimTypes), parameter :: &
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


  ! -------------------------------------------------------------------------

  SUBROUTINE ReadL2AUXData(sd_id, quantityname, l2aux, numProfs, &
       firstProf, lastProf)
    !------------------------------------------------------------------------

    ! This routine reads an l2aux file, returning a filled data structure and the !
    ! number of profiles read.

    ! Arguments

    CHARACTER (LEN=*), INTENT(IN) :: quantityname ! Name of L2AUX quantity
    INTEGER, INTENT(IN) :: sd_id ! Returned by sfstart before calling us
    INTEGER, INTENT(IN), OPTIONAL :: firstProf, lastProf ! Defaults to first and last
    TYPE( L2AUXData_T ), INTENT(OUT) :: l2aux ! Result
    INTEGER, INTENT(OUT) :: numProfs ! Number actually read

    ! Parameters

    CHARACTER (LEN=*), PARAMETER :: SZ_ERR = 'Failed to get size of &
         &dimension '
    CHARACTER (LEN=*), PARAMETER :: MLSMSG_INPUT = 'Error in input argument '
    CHARACTER (LEN=*), PARAMETER :: MLSMSG_l2auxRead = 'Unable to read l2aux &
                                                     &field:'
    INTEGER, PARAMETER :: MAXRANK = 3
    INTEGER, PARAMETER :: MAXDIMSIZES = 300
    LOGICAL, PARAMETER :: CHECKDIMSIZES = .TRUE.	! .TRUE. only while debugging

    ! Functions

    !INTEGER, EXTERNAL :: swattach, swdetach, swdiminfo, swinqdims, swrdfld

    ! Variables

    CHARACTER (LEN=80) :: list
    CHARACTER (LEN=480) :: msr

    INTEGER :: sds_index, sds_id, rank, data_type, num_attrs, dim, dim_id
    INTEGER :: dim_sizes(MAXRANK)
    INTEGER :: dim_families(MAXRANK)
    CHARACTER (LEN=LEN(quantityname)) :: sds_name
    CHARACTER (LEN=L2AUXDimNameLen) :: dim_name
    CHARACTER (LEN=1)                  :: dim_char

    INTEGER :: alloc_err, first, freq, lev, nDims, size, status
    INTEGER :: start(3), stride(3), edge(3), dims(3)
    INTEGER :: nFreqs, nLevels, nTimes, nFreqsOr1, nLevelsOr1

    LOGICAL :: firstCheck, lastCheck

    REAL, ALLOCATABLE :: realFreq(:), realSurf(:), realProf(:), real3(:,:,:)

    ! Attach to the file for reading

! found below in sfgetinfo--where we will check it for self consistency
!    l2aux%Name = quantityname

    ! find SD data set identifier
    sds_index = sfn2index(sd_id, quantityname)
    IF (sds_index == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         &get sds_index.')
         
    sds_id = sfselect(sd_id, sds_index)
    IF (sds_id == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         &get sds_id.')
         
    status = sfgetinfo(sds_id, sds_name, rank, dim_sizes, data_type, &
    & num_attrs)
    IF (status == -1) THEN
       CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         & get sf info.')
    ELSEIF (sds_name /= quantityname) THEN
       CALL MLSMessage(MLSMSG_Error, ModuleName, 'quantityname &
         & fails to match sf info.')
    ENDIF

    ! Check optional input arguments

    firstCheck = PRESENT(firstProf)
    lastCheck = PRESENT(lastProf)


    ! Uncertain what to do with those just yet
    ! Now find dimension family of dimension; e.g., MAF
    DO dim=1, rank
    	WRITE(dim_char, '(I1)') dim
    	dim_id = sfdimid(sds_id, dim)
        if(dim_id == -1) THEN
           msr = 'Failed to &
           & get dim_id for dim index number ' // dim_char
           CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
        ELSE
            status = sfgdinfo(dim_id, dim_name, dim_sizes(dim), data_type, &
            & num_attrs)
            IF(status == -1) THEN
                  msr = 'Failed to &
                  & get dim_info for dim index number ' // dim_char
                  CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ELSE
                dim_families(dim) = LinearSearchStringArray(L2AUXDimNames, dim_name)
                IF(dim_families(dim) == 0) THEN
                     msr = 'Failed to &
                     & find ' //dim_name // ' among L2AuxDimNames'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                ENDIF
            ENDIF
         ENDIF
    ENDDO
    ! Allocate result

    CALL SetupNewl2auxRecord ( dim_families, dim_sizes, l2aux )

    ! Allocate temporary arrays

    ! Read the SD
    start = 0
    stride = 1
    status = sfrdata(sds_id, start, stride, dim_sizes, l2aux%values)
    IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         & write SD.')

    ! Deallocate local variables


    ! Terminate access to the data set

    status = sfendacc(sds_id)
    IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         &end access to sds_id after reading.')

    !  After reading, detach from hdf interface

    status = sfend(sd_id)
    IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         &detach from SD file after reading.')

    !-----------------------------
  END SUBROUTINE ReadL2AUXData
  !-----------------------------

!=============================================================================
end module L2AUXData
!=============================================================================

!
! $Log$
! Revision 2.2  2000/12/04 21:48:29  pwagner
! ReadL2AUXData completed
!
! Revision 2.1  2000/12/02 01:12:00  pwagner
! Added ReadL2AUXData
!
! Revision 2.0  2000/09/05 18:57:02  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry
!

