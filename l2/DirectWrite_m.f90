! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=======================================================================================

module DirectWrite_m  ! alternative to Join/OutputAndClose methods
                      ! appropriate L2 Files

!=======================================================================================

    ! 
    ! Purpose:
    ! Write plain hdf-formatted files ala l2aux 
    ! or hdfeos-formatted l2gp for datasets that
    ! are too big to keep all chunks stored in memory
    ! or simply take too much time doing i/o
    ! so instead write them out chunk-by-chunk

  use Allocate_Deallocate, only: Allocate_test, DeAllocate_test
  use INIT_TABLES_MODULE, only: L_PRESSURE, L_ZETA, L_L2GP, L_L2AUX
  use MLSCommon, only: FindFirst, RV
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error, MLSMSG_Warning
  use OUTPUT_M, only: blanks, OUTPUT
  use STRING_TABLE, only: GET_STRING
  use VectorsModule, only: VectorValue_T

  implicit none
  private
  public :: DirectData_T, &
    & AddDirectToDatabase, &
    & DestroyDirectDatabase, DirectWrite_L2Aux, DirectWrite_L2GP, Dump, &
    & ExpandDirectDB, ExpandSDNames, &
    & SetupNewDirect

  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), parameter :: IdParm = &
    & "$Id$"
  character(len=len(idParm)) :: Id = idParm
  character(len=*), parameter :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here 
  !-----------------------------------------------------------------------------

  interface DUMP
    module procedure DumpDirectWrite
    module procedure DumpDirectDB
  end interface

  type DirectData_T
    integer :: type ! l_l2aux or l_l2gp  ! should be at least L2GPNameLen
    integer :: Handle     ! Some bit of toolkit foolishness
    integer :: NSDNames
    character(len=80), dimension(:), pointer :: sdNames => null()
    character(len=1024) :: fileNameBase ! E.g., 'H2O'
    character(len=1024) :: fileName ! E.g., '/data/../MLS..H2O..'
  end type DirectData_T
  logical, parameter :: DEEBUG = .false.
  ! For Announce_Error
  integer :: ERROR
  integer, save :: lastProfTooBigWarns = 0
  integer, parameter :: MAXNUMWARNS = 40

contains ! ======================= Public Procedures =========================

  !-------------------------------------------  AddDirectToDatabase  -----
  integer function AddDirectToDatabase( DATABASE, ITEM )

    ! This function adds a directly writeabledata type to a database
    ! creating a new database if it doesn't exist.  The result value is
    ! the size -- where it is put.

    ! Dummy arguments
    type (DirectData_T), dimension(:), pointer :: DATABASE
    type (DirectData_T), intent(in) :: ITEM

    ! Local variables
    type (DirectData_T), dimension(:), pointer :: tempDatabase
    !This include causes real trouble if you are compiling in a different 
    !directory.
    include "addItemToDatabase.f9h" 

    AddDirectToDatabase = newSize
  end function AddDirectToDatabase

  ! --------------------------------------------------------------------------

  ! This subroutine destroys a directly writeable database

  subroutine DestroyDirectDatabase ( DATABASE )

    ! Dummy argument
    type (DirectData_T), dimension(:), pointer :: DATABASE

    ! Local variables
    integer :: directIndex, status

    if ( associated(database) ) then
      do directIndex=1, size(DATABASE)
        status = 0
        if ( associated(database(directIndex)%sdNames) ) &
          & deallocate ( database(directIndex)%sdNames, stat=status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_deallocate // "sdNames" )
      enddo
       deallocate ( database, stat=status )
       if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & MLSMSG_deallocate // "database" )
    end if
  end subroutine DestroyDirectDatabase

  ! ------------------------------------------ DirectWrite_L2GP --------
  subroutine DirectWrite_L2GP ( L2gpFileHandle, &
    & quantity, quantity_precision, sdName, chunkNo, &
    & hdfVersion, fileName, createSwath )

    ! Purpose:
    ! Write standard hdfeos-formatted files ala l2gp for datasets that
    ! are too big to keep all chunks stored in memory
    ! so instead write them out profile-by-profile
    use L2GPData, only: L2GPData_T, L2GPNameLen, &
      & AppendL2GPData, DestroyL2GPContents

    integer, intent(in) :: L2gpFileHandle
    type (VectorValue_T), intent(in) :: QUANTITY
    type (VectorValue_T), pointer :: QUANTITY_PRECISION
    ! integer, intent(in) :: SDNAME       ! Name of sd in output file
    character(len=*), intent(in) :: SDNAME       ! Name of sd in output file
    integer, intent(in) :: HDFVERSION   ! Version of HDF file to write out
    integer, intent(in)              :: chunkNo
    character(len=*), intent(in), optional :: FILENAME
    logical, intent(in), optional :: createSwath
    ! Local variables
    type (L2GPData_T) :: l2gp
    integer :: OFFSET
    integer :: FIRSTINSTANCE
    integer :: LASTINSTANCE
    integer :: TOTALPROFS
    integer :: NOTOWRITE

    ! Executable code
    ! Size the problem
    firstInstance = quantity%template%noInstancesLowerOverlap + 1
    lastInstance = quantity%template%noInstances - &
      & quantity%template%noInstancesUpperOverlap
    noToWrite = lastInstance - firstInstance + 1
!     if ( noToWrite <= 0 ) then
!       call MLSMessage ( MLSMSG_Warning, 'No profiles in chunk to write' )
!       return
!     end if
    offset = quantity%template%instanceOffset
    totalProfs = quantity%template%grandTotalInstances

    ! Check sanity
    if ( offset + noToWrite - 1 > totalProfs .and. totalProfs > 0 ) &
      & call MLSMessage ( MLSMSG_Warning, ModuleName, &
      'last profile > grandTotalInstances for ' // trim(sdName) )

    ! Convert vector quantity to l2gp
    call vectorValue_to_l2gp(quantity, Quantity_precision, l2gp, &
      & sdname, chunkNo, offset=0, &
      & firstInstance=firstInstance, lastInstance=lastInstance)
    ! Output the l2gp into the file
    call AppendL2GPData(l2gp, l2gpFileHandle, &
      & sdName, filename, offset, TotalProfs, hdfVersion, createSwath)
    ! Clear up our temporary l2gp
    call DestroyL2GPContents(l2gp)
  end subroutine DirectWrite_L2GP

  ! ------------------------------------------- DirectWrite_L2Aux_FileName --------
  subroutine DirectWrite_L2Aux ( fileID, quantity, precision, sdName, &
    & hdfVersion, chunkNo, chunks )

    ! Purpose:
    ! Write plain hdf-formatted files ala l2aux for datasets that
    ! are too big to keep all chunks stored in memory
    ! so instead write them out chunk-by-chunk
    use MLSCommon, only: MLSCHUNK_T
    use MLSFiles, only: HDFVERSION_4, HDFVERSION_5
    use VectorsModule, only: VectorValue_T

    type (VectorValue_T), intent(in) :: QUANTITY
    type (VectorValue_T), pointer :: PRECISION
    ! integer, intent(in) :: SDNAME       ! Name of sd in output file
    character(len=*), intent(in) :: SDNAME       ! Name of sd in output file
    integer, intent(in) :: FILEID       ! ID of output file
    integer, intent(in) :: HDFVERSION   ! Version of HDF file to write out
    integer, intent(in) :: CHUNKNO      ! Index into chunks
    type (MLSChunk_T), dimension(:), intent(in) :: CHUNKS
    ! Local parameters
    logical, parameter :: DEEBUG = .FALSE.
    integer, parameter :: MAXFILES = 100             ! Set for an internal array
    ! executable code

    if ( (lastProfTooBigWarns <= MAXNUMWARNS) &
      & .and. &
      & (quantity%template%instanceOffset+quantity%template%noInstances - &
      & quantity%template%noInstancesLowerOverlap - &
      & quantity%template%noInstancesUpperOverlap) &
      & > &
      & quantity%template%grandTotalInstances &
      & .and. &
      & quantity%template%grandTotalInstances > 0 ) then
      lastProfTooBigWarns = lastProfTooBigWarns + 1
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'last profile > grandTotalInstances for ' // trim(sdName))
      call output('instanceOffset: ', advance='no')
      call output(quantity%template%instanceOffset, advance='yes')
      call output('noInstances: ', advance='no')
      call output(quantity%template%noInstances, advance='yes')
      call output('noInstancesLowerOverlap: ', advance='no')
      call output(quantity%template%noInstancesLowerOverlap, advance='yes')
      call output('noInstancesUpperOverlap: ', advance='no')
      call output(quantity%template%noInstancesUpperOverlap, advance='yes')
      call output('last profile: ', advance='no')
      call output((quantity%template%instanceOffset+quantity%template%noInstances - &
        & quantity%template%noInstancesLowerOverlap - &
        & quantity%template%noInstancesUpperOverlap), advance='yes')
      call output('grandTotalInstances: ', advance='no')
      call output(quantity%template%grandTotalInstances, advance='yes')
      if ( lastProfTooBigWarns > MAXNUMWARNS ) &
          & call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Max no. of warnings reached--suppressing further ones')
    endif
    select case (hdfversion)
    case (HDFVERSION_4)
      call DirectWrite_L2Aux_hdf4 ( quantity, sdName, fileID, &
        & chunkNo, chunks )
      if ( associated(precision) ) & 
        & call DirectWrite_L2Aux_hdf4 ( precision, &
        & trim(sdName) // 'precision', fileID, chunkNo, chunks )
    case (HDFVERSION_5)
      call DirectWrite_L2Aux_hdf5 ( quantity, sdName, fileID, &
        & chunkNo, chunks )
      if ( associated(precision) ) & 
        & call DirectWrite_L2Aux_hdf5 ( precision, &
        & trim(sdName) // 'precision', fileID, chunkNo, chunks )
    case default
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unsupported hdfVersion for DirectWrite_L2Aux (currently only 4 or 5)' )
    end select
  end subroutine DirectWrite_L2Aux

  ! ------------------------------------------ DirectWrite_L2Aux_hdf4 --------
  subroutine DirectWrite_L2Aux_hdf4 ( quantity, sdName, fileID, &
    & chunkNo, chunks )

    use Hdf, only: SFN2INDEX, SFSELECT, SFCREATE, &
      & SFENDACC, DFNT_FLOAT32, SFWDATA_F90
    use Intrinsic, only: L_None
    use MLSCommon, only: MLSCHUNK_T, R4, R8
    use MLSFiles, only: HDFVERSION_4
    use VectorsModule, only: VectorValue_T

    type (VectorValue_T), intent(in) :: QUANTITY
    ! integer, intent(in) :: SDNAME       ! Name of sd in output file
    character(len=*), intent(in) :: SDNAME       ! Name of sd in output file
    integer, intent(in) :: FILEID       ! ID of output file
    integer, intent(in) :: CHUNKNO      ! Index into chunks
    type (MLSChunk_T), dimension(:), intent(in) :: CHUNKS

    ! Local parameters
    integer, parameter :: MAXFILES = 100             ! Set for an internal array
    integer, parameter :: HDFVERSION = HDFVERSION_4

    ! Local variables
    integer :: SDINDEX                  ! Index of sd
    integer :: SDID                     ! Handle for sd
    integer :: STATUS                   ! Status flag
    integer :: START(3)                 ! HDF array starting position
    integer :: STRIDE(3)                ! HDF array stride
    integer :: SIZES(3)                 ! HDF array sizes
    integer :: NODIMS                   ! Also index of maf dimension
    type ( MLSChunk_T ) :: LASTCHUNK    ! The last chunk in the file
    real (r8) :: HUGER4

    ! executable code

    hugeR4 = real ( huge(0.0_r4), r8 )

    if ( quantity%template%frequencyCoordinate == L_None ) then
      noDims = 2
    else
      noDims = 3
    end if

    ! Create or access the SD
    sdIndex = sfn2index ( fileID, trim(sdName) )
    if ( sdIndex == -1 ) then
      lastChunk = chunks(size(chunks))
      sizes(noDims) = lastChunk%lastMAFIndex - lastChunk%noMAFSUpperOverlap + 1
      sizes(noDims-1) = quantity%template%noSurfs
      if ( noDims == 3 ) sizes(1) = quantity%template%noChans
      sdId  = sfCreate ( fileID, trim(sdName), DFNT_FLOAT32, &
        & noDims, sizes )
    else
      sdId = sfSelect ( fileID, sdIndex )
    end if
    if ( sdId == -1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
     & 'Error accessing SD '//trim(sdName) // ' (hdf4)')

    ! What exactly will be our contribution
    stride = 1
    start = 0
    sizes = 1
    sizes(noDims) = quantity%template%noInstances - &
      & quantity%template%noInstancesLowerOverlap - &
      & quantity%template%noInstancesUpperOverlap
    sizes(noDims-1) = quantity%template%noSurfs
    if ( noDims == 3 ) sizes(1) = quantity%template%noChans
    ! start(noDims) = quantity%template%mafIndex ( &
    !  & 1+quantity%template%noInstancesLowerOverlap )
    start(noDims) = quantity%template%instanceOffset
    ! if ( quantity%template%minorFrame ) &
    !  & start(noDims) = quantity%template%instanceOffset + 1

    ! Now write it out
    status = 0
    if ( sizes(noDims) > 0 ) &
      & status = SFWDATA_F90(sdId, start(1:noDims), &
      & stride(1:noDims), sizes(1:noDims), &
      & real ( max ( -hugeR4, min ( hugeR4, &
      &   quantity%values ( :, &
      &   1+quantity%template%noInstancesLowerOverlap : &
      &    quantity%template%noInstances - quantity%template%noInstancesUpperOverlap &
      &  ) ) ) ) )
    if ( status /= 0 ) then
      call announce_error (0,&
        & "Error writing SDS data " // trim(sdName) // " to l2aux file:  " )
    end if
    if ( DEEBUG ) then
      call output('noDims: ', advance='no')
      call output(noDims, advance='yes')
      call output('start: ', advance='no')
      call output(start, advance='yes')
      call output('stride: ', advance='no')
      call output(stride, advance='yes')
      call output('sizes: ', advance='no')
      call output(sizes, advance='yes')
      call output('shape(values): ', advance='no')
      call output(shape(quantity%values), advance='yes')
      call output('first instance: ', advance='no')
      call output(1+quantity%template%noInstancesLowerOverlap, advance='yes')
      call output('last instance: ', advance='no')
      call output(quantity%template%noInstances &
        & - quantity%template%noInstancesUpperOverlap, advance='yes')
    end if

    ! End access to the SD and close the file
    status = sfEndAcc ( sdId )
    if ( status == -1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Error ending access to direct write sd (hdf4)' )

  end subroutine DirectWrite_L2Aux_hdf4

  ! ------------------------------------------ DirectWrite_L2Aux_hdf5 --------
  subroutine DirectWrite_L2Aux_hdf5 ( quantity, sdName, fileID, &
    & chunkNo, chunks )

    use Intrinsic, only: L_None
    use MLSCommon, only: MLSCHUNK_T
    use MLSFiles, only: HDFVERSION_5
    use MLSHDF5, only: ISHDF5DSPRESENT, SaveAsHDF5DS
    use PCFHdr, only: h5_writeglobalattr
    use VectorsModule, only: VectorValue_T

    type (VectorValue_T), intent(in) :: QUANTITY
    ! integer, intent(in) :: SDNAME       ! Name of sd in output file
    character(len=*), intent(in) :: SDNAME       ! Name of sd in output file
    integer, intent(in) :: FILEID       ! ID of output file
    integer, intent(in) :: CHUNKNO      ! Index into chunks
    type (MLSChunk_T), dimension(:), intent(in) :: CHUNKS

    ! Local parameters
    integer, parameter :: MAXFILES = 100             ! Set for an internal array
    integer, parameter :: HDFVERSION = HDFVERSION_5

    ! Local variables
    logical :: already_there
    integer :: first_maf
    integer :: last_maf
    type ( MLSChunk_T ) :: LASTCHUNK    ! The last chunk in the file
    integer :: NODIMS                   ! Also index of maf dimension
    integer :: Num_qty_values
    integer :: SIZES(3)                 ! HDF array sizes
    integer :: START(3)                 ! HDF array starting position
    integer :: STRIDE(3)                ! HDF array stride
    integer :: total_DS_size
    logical, parameter :: MAYCOLLAPSEDIMS = .false.

    ! executable code
    Num_qty_values = size(quantity%values, 1)*size(quantity%values, 2)

    if ( quantity%template%frequencyCoordinate == L_None &
      & .and. MAYCOLLAPSEDIMS) then
      noDims = 2
    else
      noDims = 3
    end if

    ! Create or access the SD
    already_there = IsHDF5DSPresent(fileID, trim(sdName))
    if ( .not. already_there ) then
      lastChunk = chunks(size(chunks))
      sizes(noDims) = lastChunk%lastMAFIndex - lastChunk%noMAFSUpperOverlap + 1
      sizes(noDims-1) = quantity%template%noSurfs
      if ( noDims == 3 ) sizes(1) = quantity%template%noChans
    end if

    ! What exactly will be our contribution
    stride = 1
    start = 0
    sizes = 1
    sizes(noDims) = quantity%template%noInstances - &
      & quantity%template%noInstancesLowerOverlap - &
      & quantity%template%noInstancesUpperOverlap
    sizes(noDims-1) = quantity%template%noSurfs
    if ( noDims == 3 ) sizes(1) = quantity%template%noChans
    ! start(noDims) = quantity%template%mafIndex ( &
    !  & 1+quantity%template%noInstancesLowerOverlap )
    start(noDims) = quantity%template%instanceOffset
    ! if ( quantity%template%minorFrame ) &
    !  & start(noDims) = quantity%template%instanceOffset + 1
    first_maf = 1+quantity%template%noInstancesLowerOverlap
    last_maf = quantity%template%noInstances &
      &       - quantity%template%noInstancesUpperOverlap

    if ( DEEBUG ) then
      print *, 'sdname ', trim(sdName)
      print *, 'already_there ', already_there
      print *, 'noDims ', noDims
      print *, 'instanceOffset ', quantity%template%instanceOffset
      print *, 'minorFrame ', quantity%template%minorFrame
      print *, 'start ', start
      print *, 'sizes ', sizes
      print *, 'shape(quantity%values) ', shape(quantity%values)
      print *, 'first_maf ', first_maf
      print *, 'last_maf ', last_maf
    endif
    ! Make certain things will fit
    if ( noDims == 3 ) then
      total_DS_size = sizes(1)*sizes(2)*sizes(3)
      if ( DEEBUG ) then
        print *, 'total_DS_size ', total_DS_size
        print *, 'Num_qty_values ', Num_qty_values
      endif
      if ( total_DS_size > Num_qty_values ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Number of 3d array elements to write > number stored in qty values' )
      call SaveAsHDF5DS( fileID, trim(sdName), &
        & real( &
        &   reshape(quantity%values(:,first_maf:last_maf), sizes(1:3)) &
        & ), start, sizes, may_add_to=.true., adding_to=already_there )
    else
      total_DS_size = sizes(1)*sizes(2)
      if ( DEEBUG ) then
        print *, 'total_DS_size ', total_DS_size
        print *, 'Num_qty_values ', Num_qty_values
      endif
      if ( total_DS_size > Num_qty_values ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Number of 2d array elements to write > number stored in qty values' )
      call SaveAsHDF5DS( fileID, trim(sdName), &
        & real( &
        &   reshape(quantity%values(:,first_maf:last_maf), sizes(1:2)) &
        & ), start, sizes, may_add_to=.true., adding_to=already_there)
    endif

    ! Now some attribute stuff
    ! This first call to write dataset-specific stuff needs work
    ! basically repeat the steps you go through in SetupNewl2auxRecord;
    ! see e.g. Join
    ! call WriteL2AUXAttributes(fileID, l2aux, trim(dataProduct%name))
    ! However I'm too lazy--these files won't be archived at DAAC
    ! so the attribute writing can wait
    call h5_writeglobalattr(fileID, skip_if_already_there=.true.)

  end subroutine DirectWrite_L2Aux_hdf5

  !------------------------------------------  DumpDirectDB  -----
  subroutine DumpDirectDB ( directDB, Details )

    ! This routine dumps the DirectWrite DB

    ! Dummy arguments
    type (DirectData_T), dimension(:), pointer  :: directDB
    integer, intent(in), optional :: DETAILS

    ! Local variables
    integer :: i
    call output ( '========== DirectWrite Data Base ==========', advance='yes' )
    call output ( ' ', advance='yes' )
    if ( size(directDB) < 1 ) then
      call output ( '**** directWrite Database empty ****', advance='yes' )
      return
    endif
    do i = 1, size(directDB)
      call DumpDirectWrite(directDB(i), Details)
    end do
  end subroutine DumpDirectDB

  !------------------------------------------  DumpDirectWrite  -----
  subroutine DumpDirectWrite ( directWrite, Details )

    ! This routine dumps the DirectWrite DB

    ! Dummy arguments
    type (DirectData_T)  :: directWrite
    integer, intent(in), optional :: DETAILS

    ! Local variables
    integer :: i
    integer :: myDetails
    ! Executable code
    myDetails = 1
    if ( present(details) ) myDetails = details

    call output ( 'File Name (base): ')
    call output ( trim(directWrite%fileNameBase), advance='yes' )
    call output ( 'Full File Name  : ')
    call output ( trim(directWrite%fileName), advance='yes' )
    call output ( 'Type  : ')
    call output ( directWrite%type, advance='no' )
    call blanks ( 3 )
    if ( directWrite%type == l_l2aux ) then
      call output ( '(l2aux)', advance='yes')
    elseif ( directWrite%type == l_l2gp ) then
      call output ( '(l2gp)', advance='yes')
    else
      call output ( '(unknown)', advance='yes')
    endif
    call output ( 'Handle  : ')
    call output ( directWrite%Handle, advance='yes' )
    call output ( 'Num sci. datasets  : ')
    call output ( directWrite%NSDNames, advance='yes' )
    if ( .not. associated(directWrite%SDNames)  .or. myDetails < 1 ) then
      call output ( 's.d. names not associated', advance='yes' )
      return
    endif
    do i = 1, size(directWrite%sdNames)
      call blanks ( 3 )
      call output(trim(directWrite%SDNames(i)), advance='yes')
    end do
  end subroutine DumpDirectWrite

  !------------------------------------------  ExpandDirectDB  -----
  subroutine ExpandDirectDB ( directDB, fileNameBase, directData, isNew )

    ! This routine expands the DirectWrite DB if necessary
    ! Returning either a match or a new one

    ! Dummy arguments
    type (DirectData_T), dimension(:), pointer  :: directDB
    type (DirectData_T), pointer :: directData
    character(len=*), intent(in) :: fileNameBase
    logical, intent(out) :: isNew

    ! Local variables
    integer :: dbID
    character(len=80), dimension(:), pointer :: sdNames => null()
    type (DirectData_T)                      :: tempDirectData
    ! Begin executable
    isNew = .true.
    ! Do we have any database yet?
    if ( .not. associated(directDB) ) then
      ! print *, 'Setting up initial database with ', trim(fileNameBase)
      call SetupNewDirect(tempDirectData, 0)
      dbID = AddDirectToDatabase( directDB, tempDirectData )
      directData => directDB(dbID)
    else
      ! Check if the fileNameBase already there
      dbID = FindFirst(fileNameBase == directDB%fileNameBase)
      if ( dbID > 0 ) then
        directData => directDB(dbID)
        isNew = .false.
        ! print *, 'Recognized ', trim(fileNameBase), ' in ', trim(directData%FileName)
      else
        call SetupNewDirect(tempDirectData, 0)
        dbID = AddDirectToDatabase( directDB, tempDirectData )
        directData => directDB(dbID)
        ! print *, 'Adding to database ', trim(fileNameBase)
      endif
    endif
  end subroutine ExpandDirectDB

  !------------------------------------------  ExpandSDNames  -----
  subroutine ExpandSDNames ( directData, sdName )

    ! This routine adds to the sdNames arrays for a DirectWrite if necessary

    ! Dummy arguments
    type (DirectData_T), intent(inout)  :: directData
    character(len=*), intent(in) :: sdName

    ! Local variables
    logical :: alreadyThere
    character(len=80), dimension(:), pointer :: sdNames => null()
    integer :: newsize
    integer :: status
    ! Check if the sdName already there
    ! print *, 'Check if the sdName already there'
    if ( associated(directData%sdNames) ) then
      newSize=SIZE(directData%sdNames)+1
      alreadyThere = any(sdName == directData%sdNames)
      ! if ( alreadyThere ) print *, trim(sdName), ' already there'
      if ( alreadyThere ) return
    else
      ! print *, 'Allocating for ', trim(sdName)
      newSize=1
    end if
    directData%NsdNames = newsize
    ! print *, ' allocating ', newsize
    allocate(sdNames(newSize),STAT=status)
    if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "temp sdnames")
    if ( newSize>1 ) sdNames(1:newSize-1) = directData%sdNames
    if ( ASSOCIATED(directData%sdNames) ) &
      & deallocate ( directData%sdNames, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_DeAllocate // "sdNames" )
    directData%sdNames => sdNames
    directData%sdNames(newSize) = sdName
  end subroutine ExpandSDNames

  !------------------------------------------  SetupNewDirect  -----
  subroutine SetupNewDirect ( directData, NsdNames )

    ! This routine sets up the arrays for a directWrite

    ! Dummy arguments
    type (DirectData_T), intent(inout)  :: directData
    integer, intent(in) :: NsdNames            ! Dimensions

    ! Local variables
    ! Allocate the sdNames
    directData%NSDNames = NSDNames
    nullify(directData%sdNames)
    if ( NSDNames <= 0 ) return
    call allocate_test ( directData%sdNames, NsdNames, "directData%sdNames", &
         & ModuleName )
  end subroutine SetupNewDirect

! =====     Private Procedures     =====================================
  ! ---------------------------------------------  vectorValue_to_l2gp  -----
  subroutine vectorValue_to_l2gp (QUANTITY, Quantity_precision, l2gp, &
    & name, chunkNo, &
    & offset, firstInstance, lastInstance)
    use Intrinsic, only: L_None
    use L2GPData, only: L2GPData_T, RGP, &
      & SetupNewl2gpRecord, &
      & ExpandL2GPDataInPlace
    type (VectorValue_T), intent(in) :: QUANTITY
    type (VectorValue_T), pointer :: QUANTITY_PRECISION
    type (L2GPData_T)                :: l2gp
    character(len=*), intent(in)     :: name
    ! integer, intent(in)            :: nameIndex
    integer, intent(in)              :: chunkNo
    integer, intent(in), optional    :: offset
    integer, intent(in), optional    :: firstInstance
    integer, intent(in), optional    :: lastInstance
    ! Local variables
    integer :: noSurfsInL2GP
    integer :: noFreqsInL2GP
    integer :: useFirstInstance
    integer :: useLastInstance
    integer :: firstProfile
    integer :: lastProfile
    integer :: nProfiles

    real(rv) :: HUGERGP
    
    hugeRgp = real ( huge(0.0_rgp), rv )

    ! Work out what to do with the first and last instance information
    
    if ( present(firstInstance) ) then
      useFirstInstance = firstInstance
    else
      useFirstInstance = quantity%template%noInstancesLowerOverlap+1
    end if

    if ( present(lastInstance) ) then
      useLastInstance = lastInstance
    else
      useLastInstance = quantity%template%noInstances- &
        & quantity%template%noInstancesUpperOverlap
    end if

    ! noOutputInstances = useLastInstance-useFirstInstance+1
    ! Now create an empty L2GP record with this dimension

    if (any(quantity%template%verticalCoordinate == (/l_Pressure, l_Zeta /) )) then
      noSurfsInL2GP = quantity%template%noSurfs
    else
      noSurfsInL2GP = 0
    end if

    if ( quantity%template%frequencyCoordinate == l_None) then
       noFreqsInL2GP=0
    else
       noFreqsInL2GP=quantity%template%noChans
    end if
    
    if ( present(offset) ) then
      firstProfile=offset+1
    else
      firstProfile=1
    endif
    nProfiles = useLastInstance - useFirstInstance + 1
    lastProfile = firstProfile - 1 + nProfiles
    call SetupNewl2gpRecord ( l2gp, noFreqsInL2GP, noSurfsInL2GP, lastProfile )
    ! Setup the standard stuff, only pressure as it turns out.
    if ( quantity%template%verticalCoordinate == l_Pressure ) &
      & l2gp%pressures = quantity%template%surfs(:,1)
    if ( quantity%template%verticalCoordinate == l_Zeta ) &
      & l2gp%pressures = 10.0**(-quantity%template%surfs(:,1))
    ! It inherits its quantity type from the quantity template
    l2gp%quantityType=quantity%template%quantityType
    ! Do something about frequency
    if ( associated ( quantity%template%frequencies ) ) then
      l2gp%frequency = quantity%template%frequencies
    else
      l2gp%frequency = 0.0
    end if
    ! call ExpandL2GPDataInPlace ( l2gp, &
    !  & lastProfile-firstProfile+1 )

    ! Now copy the information from the quantity to the l2gpData

    ! name is an integer, but L2GP%name is Character data
    l2gp%nameIndex = 0
    l2gp%name = name

    ! Now fill the data, first the geolocation
    l2gp%latitude(firstProfile:lastProfile) = &
      & quantity%template%geodLat(1,useFirstInstance:useLastInstance)
    l2gp%longitude(firstProfile:lastProfile) = &
      & quantity%template%lon(1,useFirstInstance:useLastInstance)
    l2gp%solarTime(firstProfile:lastProfile) = &
      & quantity%template%solarTime(1,useFirstInstance:useLastInstance)
    l2gp%solarZenith(firstProfile:lastProfile) = &
      & quantity%template%solarZenith(1,useFirstInstance:useLastInstance)
    l2gp%losAngle(firstProfile:lastProfile) = &
      & quantity%template%losAngle(1,useFirstInstance:useLastInstance)
    l2gp%geodAngle(firstProfile:lastProfile) = &
      & quantity%template%phi(1,useFirstInstance:useLastInstance)
    l2gp%time(firstProfile:lastProfile) = &
      & quantity%template%time(1,useFirstInstance:useLastInstance)
    l2gp%chunkNumber(firstProfile:lastProfile)=chunkNo

    ! Now the various data quantities.

    l2gp%l2gpValue(:,:,firstProfile:lastProfile) = &
      & reshape ( max ( -hugeRgp, min ( hugeRgp, &
      &   quantity%values(:,useFirstInstance:useLastInstance) ) ), &
      &  (/max(l2gp%nFreqs,1),max(l2gp%nLevels,1),lastProfile-firstProfile+1/))
    if (associated(quantity_precision)) then
      l2gp%l2gpPrecision(:,:,firstProfile:lastProfile) = &
        & reshape ( max ( -hugeRgp, min ( hugeRgp, &
        &   quantity_precision%values(:,useFirstInstance:useLastInstance) ) ), &
        &  (/max(l2gp%nFreqs,1),max(l2gp%nLevels,1),lastProfile-firstProfile+1/))
    else
      l2gp%l2gpPrecision(:,:,firstProfile:lastProfile) = 0.0
    end if
    l2gp%status(firstProfile:lastProfile) = 'G'
    l2gp%quality(firstProfile:lastProfile) = 0.0
  end subroutine vectorValue_to_l2gp

  ! ---------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( Where, Full_message, Code, Penalty )

    use LEXER_CORE, only: PRINT_SOURCE
    use TREE, only: SOURCE_REF

    integer, intent(in) :: Where   ! Tree node where error was noticed
    character(LEN=*), intent(in) :: Full_Message
    integer, intent(in), optional :: Code    ! Code for error message
    integer, intent(in), optional :: Penalty
    integer :: myPenalty

    myPenalty = 1
    if ( present(penalty) ) myPenalty = penalty
    error = max(error,myPenalty)

    call output ( '***** At ' )
    if ( where > 0 ) then
      call print_source ( source_ref(where) )
    else
      call output ( '(no lcf node available)' )
    end if
    call output ( ModuleName // ' complained: ' )

    call output ( trim(full_message), advance='yes', &
      & from_where=ModuleName )
    if ( present(code) ) then
      select case ( code )
      end select
    end  if
  end subroutine ANNOUNCE_ERROR

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module DirectWrite_m

! $Log$
! Revision 2.13  2003/11/14 23:38:45  pwagner
! Uses DirectWrite databse in preference to repeated calls to toolkit
!
! Revision 2.12  2003/08/28 23:51:39  livesey
! Various bug fixes and simplifications
!
! Revision 2.11  2003/08/01 20:39:34  pwagner
! Skip warning mesg about grandTotalInstances when 0
!
! Revision 2.10  2003/07/18 16:05:26  pwagner
! Stops printing warnings after 40 times
!
! Revision 2.9  2003/07/15 23:41:47  pwagner
! l2aux always rank 3 unless MAYCOLLAPSEDIMS; disabled most printing
!
! Revision 2.8  2003/07/09 21:49:53  pwagner
! Tries to figure out in advance whether to create swath
!
! Revision 2.7  2003/07/07 21:03:43  pwagner
! Removed DirectWrite_L2Aux_FileName; tries to deal sensibly with all-overlap chunks
!
! Revision 2.6  2003/07/07 17:32:30  livesey
! New approach to DirectWrite
!
! Revision 2.5  2003/07/02 00:55:27  pwagner
! Some improvements in DirectWrites of l2aux, l2gp
!
! Revision 2.4  2003/06/25 18:24:57  pwagner
! A fix to vectorValue_to_l2gp
!
! Revision 2.3  2003/06/24 23:53:27  pwagner
! Allows unassociated precisions as args to direct write
!
! Revision 2.2  2003/06/23 23:55:17  pwagner
! Added DirectData_T to keep track of data written directly
!
! Revision 2.1  2003/06/20 19:43:16  pwagner
! First commit
!
