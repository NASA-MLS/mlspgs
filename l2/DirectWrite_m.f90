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

  use Allocate_Deallocate, only: Allocate_test
  use INIT_TABLES_MODULE, only: L_PRESSURE, L_ZETA
  use MLSCommon, only: RV
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error
  use OUTPUT_M, only: OUTPUT
  use STRING_TABLE, only: GET_STRING
  use VectorsModule, only: VectorValue_T

  implicit none
  private
  public :: DirectData_T, AddDirectToDatabase, DestroyDirectDatabase, &
    & DirectWrite_L2Aux, DirectWrite_L2GP, SetupNewDirect

  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), parameter :: IdParm = &
    & "$Id$"
  character(len=len(idParm)) :: Id = idParm
  character(len=*), parameter :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here 
  !-----------------------------------------------------------------------------

  interface DirectWrite_L2GP
    module procedure DirectWrite_L2GP_fileID
    module procedure DirectWrite_L2GP_fileName
  end interface

  interface DirectWrite_L2Aux
    module procedure DirectWrite_L2Aux_fileID
    module procedure DirectWrite_L2Aux_fileName
  end interface

  type DirectData_T
    integer :: type ! l_l2aux or l_l2gp  ! should be at least L2GPNameLen
    character(len=80), dimension(:), pointer :: sdNames => null()
    character(len=1024) :: fileNameBase ! E.g., 'H2O'
  end type DirectData_T
  ! For Announce_Error
  integer :: ERROR

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

  ! ----------------------------------------------- DirectWrite_L2GP --------
  subroutine DirectWrite_L2GP_filename ( quantity, quantity_precision, sdName, &
    & file, chunkNo, hdfVersion, firstProf, lastProf )

    ! Purpose:
    ! Write standard hdfeos-formatted files ala l2gp for datasets that
    ! are too big to keep all chunks stored in memory
    ! so instead write them out profile-by-profile
    use Hdf, only: DFACC_CREATE, DFACC_RDWR
    use L2GPData, only: L2GPNameLen
    use MLSFiles, only:  &
      & GetPCFromRef, split_path_name, mls_exists, &
      & mls_io_gen_openF, mls_io_gen_closeF
    use MLSL2Options, only: TOOLKIT
    use MLSPCF2, only: mlspcf_l2gp_start, mlspcf_l2gp_end

    type (VectorValue_T), intent(in) :: QUANTITY
    type (VectorValue_T), pointer :: QUANTITY_PRECISION
    ! integer, intent(in) :: SDNAME       ! Name of sd in output file
    character(len=*), intent(in) :: SDNAME       ! Name of sd in output file
    integer, intent(in) :: FILE         ! Name of output file
    integer, intent(in) :: HDFVERSION   ! Version of HDF file to write out
    integer, intent(in)              :: chunkNo
    integer, intent(in), optional :: firstProf, lastProf ! Defaults to first and last
    ! Local parameters
    logical, parameter :: DEBUG = .FALSE.
    character (len=132) :: FILE_BASE    ! From the FILE location in string table
    character (len=1024) :: FILENAME    ! The actual filename
    integer :: L2gpFileHandle, L2gp_Version
    integer, parameter :: MAXFILES = 100             ! Set for an internal array
    character (len=132) :: path
    character (len=L2GPNameLen) :: swathname    ! From the FILE location in string table
    integer :: fileaccess               ! DFACC_CREATE or DFACC_RDWR
    integer :: record_length
    integer :: ReturnStatus
    ! executable code

    ! Setup information, sanity checks etc.
    call get_string ( file, file_base, strip=.true. )
    ! call get_string ( sdname, swathname, strip=.true. )
    swathname=sdName
    L2gp_Version = 1
    if ( TOOLKIT ) then
      call split_path_name(file_base, path, file_base)
      if ( DEBUG ) call output('file_base after split: ', advance='no')
      if ( DEBUG ) call output(trim(file_base), advance='yes')

      L2gpFileHandle = GetPCFromRef(file_base, &
      & mlspcf_l2gp_start, mlspcf_l2gp_end, &
      & TOOLKIT, returnStatus, L2gp_Version, DEBUG, &
      & exactName=Filename)
    else
      Filename = file_base
      returnStatus = 0
    end if
    if ( returnStatus /= 0 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
         &  "Error finding l2gp file matching:  "// trim(file_base))
    endif
    if ( mls_exists(trim(Filename)) == 0 ) then
      fileaccess = DFACC_RDWR
    else
      fileaccess = DFACC_CREATE
    endif
    L2gpFileHandle = mls_io_gen_openF('sw', .true., returnStatus, &
        & record_length, FileAccess, FileName, hdfVersion=hdfVersion)
    call DirectWrite_L2GP_fileID( L2gpFileHandle, &
      & quantity, quantity_precision, sdName, chunkNo, &
      & hdfVersion, firstProf, lastProf )
    ReturnStatus = mls_io_gen_closeF('swclose', L2gpFileHandle, FileName=FileName, &
      & hdfVersion=hdfVersion)
    if ( ReturnStatus /= 0 ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
       & "Unable to close L2gp file: " // trim(FileName) // ' after reading')
  end subroutine DirectWrite_L2GP_filename

  ! ------------------------------------------ DirectWrite_L2GP_fileID --------
  subroutine DirectWrite_L2GP_fileID ( L2gpFileHandle, &
    & quantity, quantity_precision, sdName, chunkNo, &
    & hdfVersion, firstProf, lastProf )

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
    integer, intent(in), optional :: firstProf, lastProf ! Defaults to first and last
    ! Local parameters
    type (L2GPData_T) :: l2gp
    logical, parameter :: DEBUG = .FALSE.
    character (len=L2GPNameLen) :: swathname    ! From the FILE location in string table
    integer :: offset
    integer :: myFirstProf
    integer :: myLastProf
    integer :: FirstInstance
    integer :: LastInstance
    integer :: TotNumProfs
    ! executable code qty%template%instanceOffset + 1
    
    ! Non-overlapped portion of quantity:
    FirstInstance = quantity%template%noInstancesLowerOverlap+1
    LastInstance = quantity%template%noInstances- &
        & quantity%template%noInstancesUpperOverlap
        
    myFirstProf = 1
    if ( present(firstprof) ) myFirstProf = firstprof
    myLastProf = 1
    if ( present(lastprof) ) myLastProf = lastprof

    TotNumProfs = myLastProf
    offset = myFirstProf - 1
    if ( .not. present(firstprof) ) offset=quantity%template%instanceOffset
    call vectorValue_to_l2gp(quantity, Quantity_precision, l2gp, &
      & sdname, chunkNo, offset=0, &
      & firstInstance=firstInstance, lastInstance=lastInstance)
    call AppendL2GPData(l2gp, l2gpFileHandle, &
      & swathname, offset, TotNumProfs, hdfVersion)
    call DestroyL2GPContents(l2gp)
  end subroutine DirectWrite_L2GP_fileID

  ! ------------------------------------------- DirectWrite_L2Aux_FileName --------
  subroutine DirectWrite_L2Aux_FileName ( quantity, precision, sdName, file, &
    & hdfVersion, chunkNo, chunks )

    ! If slave, requests permission from master to access file
    ! Opens File (if it exists already) or else creates it
    ! Calls for direct write by file handle
    ! Closes file
    ! If slave, tells master it has finished
    use Hdf, only: DFACC_CREATE, DFACC_RDWR
    use L2ParInfo, only: PARALLEL, REQUESTDIRECTWRITEPERMISSION, &
      & FINISHEDDIRECTWRITE
    use MLSCommon, only: MLSCHUNK_T
    use MLSFiles, only: &
      & GetPCFromRef, split_path_name, mls_sfstart, mls_sfend
    use MLSL2Options, only: TOOLKIT
    use MLSPCF2, only: mlspcf_l2fwm_full_start, mlspcf_l2fwm_full_end
    use VectorsModule, only: VectorValue_T

    type (VectorValue_T), intent(in) :: QUANTITY
    type (VectorValue_T), pointer :: PRECISION
    ! integer, intent(in) :: SDNAME       ! Name of sd in output file
    character(len=*), intent(in) :: SDNAME       ! Name of sd in output file
    integer, intent(in) :: FILE         ! Name of output file
    integer, intent(in) :: HDFVERSION   ! Version of HDF file to write out
    integer, intent(in) :: CHUNKNO      ! Index into chunks
    type (MLSChunk_T), dimension(:), intent(in) :: CHUNKS
    ! Local parameters
    logical, parameter :: DEBUG = .FALSE.
    character (len=132) :: FILE_BASE    ! From the FILE location in string table
    character (len=1024) :: FILENAME    ! The actual filename
    integer :: fileID, L2fwm_Version
    integer, parameter :: MAXFILES = 100             ! Set for an internal array
    integer, save :: NOCREATEDFILES=0   ! Number of files created

    character (len=132) :: path
    integer :: status
    integer :: ReturnStatus
    logical :: createFile
    logical :: you_may
    ! Saved variable - used to work out file information.
    integer, dimension(maxFiles), save :: CREATEDFILENAMES = 0
    ! executable code

    ! Setup information, sanity checks etc.
    call get_string ( file, file_base, strip=.true. )
    l2fwm_Version = 1
    if ( TOOLKIT ) then
      call split_path_name(file_base, path, file_base)
      if ( DEBUG ) call output('file_base after split: ', advance='no')
      if ( DEBUG ) call output(trim(file_base), advance='yes')

      fileID = GetPCFromRef(file_base, mlspcf_l2fwm_full_start, &
      & mlspcf_l2fwm_full_end, &
      & TOOLKIT, returnStatus, L2fwm_Version, DEBUG, &
      & exactName=Filename)
    else
      Filename = file_base
      returnStatus = 0
    end if
    if ( returnStatus /= 0 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
         &  "Error finding l2gp file matching:  "// trim(file_base))
    endif
    ! Setup information, sanity checks etc.
    ! If we're a slave, we need to request permission from the master.
    if ( parallel%slave ) then
      call RequestDirectWritePermission ( file, createFile, .true., you_may )
    else
      createFile = .not. any ( createdFilenames == file )
      if ( createFile ) then
        noCreatedFiles = noCreatedFiles + 1
        if ( noCreatedFiles > maxFiles ) call MLSMessage ( &
          & MLSMSG_Error, ModuleName, 'Too many direct write files (hdf4)' )
        createdFilenames ( noCreatedFiles ) = file
      end if
    end if

    ! Create or open the file
    if ( createFile ) then
      fileID = mls_sfstart ( trim(filename), DFACC_CREATE, &
       & hdfVersion=hdfVersion )
    else
      fileID = mls_sfstart ( trim(filename), DFACC_RDWR, hdfVersion=hdfVersion )
    end if
    call DirectWrite_L2Aux_FileID ( fileID, quantity, precision, sdName, &
      & hdfVersion, chunkNo, chunks )

    status = mls_sfend( fileID, hdfVersion=hdfVersion)
    if ( status == -1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Error ending closing direct write l2aux file' )

    ! Tell master we're done
    if ( parallel%slave ) call FinishedDirectWrite
  end subroutine DirectWrite_L2Aux_FileName

  ! ------------------------------------------- DirectWrite_L2Aux_FileName --------
  subroutine DirectWrite_L2Aux_FileID ( fileID, quantity, precision, sdName, &
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
    logical, parameter :: DEBUG = .FALSE.
    integer, parameter :: MAXFILES = 100             ! Set for an internal array
    ! executable code

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
  end subroutine DirectWrite_L2Aux_FileID

  ! ----------------------------------------------- DirectWrite_L2Aux_hdf4 --------
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
    sizes(noDims) = quantity%template%noInstances - &
      & quantity%template%noInstancesLowerOverlap - &
      & quantity%template%noInstancesUpperOverlap
    sizes(noDims-1) = quantity%template%noSurfs
    if ( noDims == 3 ) sizes(1) = quantity%template%noChans
    ! start(noDims) = quantity%template%mafIndex ( &
    !  & 1+quantity%template%noInstancesLowerOverlap )
    start(noDims) = quantity%template%instanceOffset
    if ( quantity%template%minorFrame ) &
      & start(noDims) = quantity%template%instanceOffset + 1

    ! Now write it out
    status = SFWDATA_F90(sdId, start(1:noDims), &
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

    ! End access to the SD and close the file
    status = sfEndAcc ( sdId )
    if ( status == -1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Error ending access to direct write sd (hdf4)' )

  end subroutine DirectWrite_L2Aux_hdf4

  ! ----------------------------------------------- DirectWrite_L2Aux_hdf5 --------
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
    logical, parameter :: DEEBUG = .false.
    integer :: first_maf
    integer :: last_maf
    type ( MLSChunk_T ) :: LASTCHUNK    ! The last chunk in the file
    integer :: NODIMS                   ! Also index of maf dimension
    integer :: Num_qty_values
    integer :: SIZES(3)                 ! HDF array sizes
    integer :: START(3)                 ! HDF array starting position
    integer :: STRIDE(3)                ! HDF array stride
    integer :: total_DS_size

    ! executable code
    Num_qty_values = size(quantity%values, 1)*size(quantity%values, 2)

    ! Create or access the SD
    already_there = IsHDF5DSPresent(fileID, trim(sdName))
    if ( .not. already_there ) then
      lastChunk = chunks(size(chunks))
      sizes(noDims) = lastChunk%lastMAFIndex - lastChunk%noMAFSUpperOverlap + 1
      sizes(noDims-1) = quantity%template%noSurfs
      if ( noDims == 3 ) sizes(1) = quantity%template%noChans
    end if

    if ( quantity%template%frequencyCoordinate == L_None ) then
      noDims = 2
    else
      noDims = 3
    end if

    ! What exactly will be our contribution
    stride = 1
    start = 0
    sizes(noDims) = quantity%template%noInstances - &
      & quantity%template%noInstancesLowerOverlap - &
      & quantity%template%noInstancesUpperOverlap
    sizes(noDims-1) = quantity%template%noSurfs
    if ( noDims == 3 ) sizes(1) = quantity%template%noChans
    ! start(noDims) = quantity%template%mafIndex ( &
    !  & 1+quantity%template%noInstancesLowerOverlap )
    start(noDims) = quantity%template%instanceOffset
    if ( quantity%template%minorFrame ) &
      & start(noDims) = quantity%template%instanceOffset + 1
    first_maf = 1+quantity%template%noInstancesLowerOverlap
    last_maf = quantity%template%noInstances &
      &       - quantity%template%noInstancesUpperOverlap

    if ( DEEBUG ) then
      print *, 'sdname ', trim(sdName)
      print *, 'already_there ', already_there
      print *, 'noDims ', noDims
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
        & ), start, sizes, may_add_to=.true., adding_to=already_there)
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
        & real(quantity%values(:,first_maf:last_maf)), &
        & start, sizes, may_add_to=.true., adding_to=already_there)
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

  !------------------------------------------  SetupNewDirect  -----
  subroutine SetupNewDirect ( directData, NsdNames )

    ! This routine sets up the arrays for an l2gp datatype.

    ! Dummy arguments
    type (DirectData_T), intent(inout)  :: directData
    integer, intent(in) :: NsdNames            ! Dimensions

    ! Local variables
    ! Allocate the sdNames

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
    call output ( ' OutputAndClose complained: ' )

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
