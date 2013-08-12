! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

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

  use ALLOCATE_DEALLOCATE, only: ALLOCATE_TEST
  use INIT_TABLES_MODULE, only: L_PRESSURE, L_ZETA, &
    & L_L2GP, L_L2AUX, L_L2DGG, L_L2FWM
  use MLSCOMMON, only: DEFAULTUNDEFINEDVALUE, MLSFILE_T
  use MLSKINDS, only: RV
  use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ALLOCATE, MLSMSG_DEALLOCATE, &
    & MLSMSG_ERROR, MLSMSG_WARNING
  use MLSFINDS, only: FINDFIRST
  use MLSSTRINGLISTS, only: SWITCHDETAIL
  use OUTPUT_M, only: BLANKS, OUTPUT, OUTPUTNAMEDVALUE
  use STRING_TABLE, only: GET_STRING
  use TOGGLES, only: SWITCHES
  use VECTORSMODULE, only: VECTOR_T, VECTORVALUE_T, DUMP

  implicit none
  private
  public :: DirectData_T, &
    & AddDirectToDatabase, &
    & DestroyDirectDatabase, DirectWrite, Dump, &
    & ExpandDirectDB, ExpandSDNames, FileNameToID, &
    & SetupNewDirect

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  interface DirectWrite
    module procedure DirectWrite_L2Aux_MF
    module procedure DirectWrite_L2GP_MF
    module procedure DirectWriteVector_L2Aux_MF
    module procedure DirectWriteVector_L2GP_MF
  end interface

  interface DUMP
    module procedure DumpDirectWrite
    module procedure DumpDirectDB
  end interface

  type DirectData_T
    integer :: type       ! l_l2aux or l_l2gp
    integer :: autoType   ! 0 if not to be filled automatically, else type
    integer :: fileIndex  ! Index into character tables
    integer :: Handle     ! Set by toolkit
    integer :: NSDNames   ! How many sources
    character(len=80), dimension(:), pointer :: sdNames => null() ! source names
    character(len=1024) :: fileNameBase ! E.g., 'H2O'
    character(len=1024) :: fileName ! E.g., '/data/../MLS..H2O...he5'
  end type DirectData_T
  integer, parameter :: S2US  = 1000000 ! How many microseconds in a s
  integer, parameter :: DELAY = 1*S2US  ! How long to sleep in microseconds
  logical, parameter :: DEEBUG = .false.
  logical, parameter :: MAYWRITEPOSTOVERLAPS = .true.
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

  ! ------------------------------------------ DirectWriteVector_L2GP_MF --------
  subroutine DirectWriteVector_L2GP_MF ( L2gpFile, &
    & vector, &
    & chunkNo, HGrids, createSwath, lowerOverlap, upperOverlap )

    ! Purpose:
    ! Write standard hdfeos-formatted files ala l2gp for datasets that
    ! are too big to keep all chunks stored in memory
    ! so instead write them out profile-by-profile
    use HGRIDSDATABASE, only: HGRID_T
    ! Args:
    type(MLSFile_T)               :: L2GPFile
    type (Vector_T), intent(in)   :: VECTOR
    integer, intent(in)              :: chunkNo
    type (HGrid_T), dimension(:), pointer ::     HGrids
    logical, intent(in), optional :: createSwath
    logical, intent(in), optional :: lowerOverlap
    logical, intent(in), optional :: upperOverlap
    ! Local variables
    integer                       :: j
    type (VectorValue_T), pointer :: QUANTITY
    type (VectorValue_T), pointer :: precision
    type (VectorValue_T), pointer :: quality
    type (VectorValue_T), pointer :: status
    type (VectorValue_T), pointer :: Convergence
    character(len=32)             :: SDNAME       ! Name of sd in output file
    ! Executable
    nullify(precision, quality, status, convergence)
    do j = 1, size(vector%quantities)
      quantity => vector%quantities(j)
      call get_string( quantity%template%name, sdname )
      call DirectWrite_L2GP_MF ( L2gpFile, &
        & quantity, precision, quality, status, Convergence, &
        & sdName, chunkNo, HGrids, createSwath, lowerOverlap, upperOverlap )
    enddo
  end subroutine DirectWriteVector_L2GP_MF

  ! ------------------------------------------ DirectWrite_L2GP_MF --------
  subroutine DirectWrite_L2GP_MF ( L2gpFile, &
    & quantity, precision, quality, status, Convergence, &
    & sdName, chunkNo, HGrids, createSwath, lowerOverlap, upperOverlap )

    ! Purpose:
    ! Write standard hdfeos-formatted files ala l2gp for datasets that
    ! are too big to keep all chunks stored in memory
    ! so instead write them out profile-by-profile
    use HDF, only: DFACC_CREATE, DFACC_RDONLY, DFACC_RDWR
    use HGRIDSDATABASE, only: HGRID_T
    use L2GPDATA, only: L2GPDATA_T, &
      & APPENDL2GPDATA, DESTROYL2GPCONTENTS, DUMP
    use READAPRIORI, only: WRITEAPRIORIATTRIBUTES
    ! Args:
    type(MLSFile_T)               :: L2GPFile
    type (VectorValue_T), intent(in) :: QUANTITY
    type (VectorValue_T), pointer :: precision
    type (VectorValue_T), pointer :: quality
    type (VectorValue_T), pointer :: status
    type (VectorValue_T), pointer :: Convergence
    character(len=*), intent(in) :: SDNAME       ! Name of sd in output file
    integer, intent(in)              :: chunkNo
    type (HGrid_T), dimension(:), pointer ::     HGrids
    logical, intent(in), optional :: createSwath
    logical, intent(in), optional :: lowerOverlap
    logical, intent(in), optional :: upperOverlap
    ! Local variables
    ! logical :: DeeBug
    integer :: FIRSTINSTANCE
    integer :: GRANDTOTALINSTANCES
    type (L2GPData_T) :: l2gp
    integer :: LASTINSTANCE
    integer :: OFFSET
    character(len=8) :: overlaps ! 'lower', 'upper', or 'none'
    integer :: NOTOWRITE
    integer :: TOTALPROFS

    logical :: verbose
    ! Executable
    verbose = ( switchDetail(switches, 'direct') > -1 )
    if ( present(createSwath) ) then
      if ( createSwath ) then
        call outputNamedValue( '  creating swath', trim(sdname) )
      elseif ( verbose ) then
        call outputNamedValue( '  adding to swath', trim(sdname) )
      endif
    endif
    ! Size the problem
    overlaps = 'none'
    if ( present(lowerOverlap) ) then
      if ( lowerOverlap ) overlaps = 'lower'
    endif
    if ( present(upperOverlap) ) then
      if ( upperOverlap ) overlaps = 'upper'
    endif
    select case ( overlaps )
    case ( 'lower' )
      firstInstance =  1
      lastInstance = quantity%template%noInstancesLowerOverlap
      ! DeeBug = .true.
    case ( 'upper' )
      firstInstance =  quantity%template%noInstances - &
        & quantity%template%noInstancesUpperOverlap + 1
      lastInstance = quantity%template%noInstances
      ! DeeBug = .true.
    case ( 'none' )
      firstInstance = quantity%template%noInstancesLowerOverlap + 1
      lastInstance = quantity%template%noInstances - &
        & quantity%template%noInstancesUpperOverlap
    end select
    noToWrite = lastInstance - firstInstance + 1
    offset = quantity%template%instanceOffset
    grandtotalinstances = quantity%template%grandTotalInstances
    totalProfs = grandtotalinstances
    ! Check sanity
    if ( offset + noToWrite - 1 > grandTotalInstances .and. grandTotalInstances > 0 ) then
      call output('offset: ', advance='no')
      call output(offset, advance='no')
      call output('   noToWrite: ', advance='no')
      call output(noToWrite, advance='no')
      call output('   grandTotalInstances: ', advance='no')
      call output(grandTotalInstances, advance='yes')
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'last profile > grandTotalInstances for ' // trim(sdName), &
        & MLSFile=L2GPFile )
      ! Last-ditch effort resets offset
      offset = max(0, min(offset, grandTotalInstances - noToWrite))
    endif
    if ( L2GPFile%access == DFACC_RDONLY )  &
      & call MLSMessage(MLSMSG_Error, ModuleName, &
      & 'l2gp file is rdonly', MLSFile=L2GPFile)
    ! Convert vector quantity to l2gp
    call vectorValue_to_l2gp( quantity, &
      & precision, quality, status, convergence, l2gp, &
      & sdname, chunkNo, HGrids, offset=0, &
      & firstInstance=firstInstance, lastInstance=lastInstance)
    ! Output the l2gp into the file
    ! if ( l2gp%name == 'Temperature-InitPtan' .and. deebug ) then
    if ( deebug ) then
      call output('firstInstance: ', advance='no')
      call output(firstInstance, advance='no')
      call output('  lastInstance: ', advance='no')
      call output(lastInstance, advance='no')
      call output('  TotalProfs: ', advance='no')
      call output(TotalProfs, advance='yes')
      call dump(l2gp, Details=-1)
    endif
    if ( verbose ) call outputNamedValue( 'DW L2GP qty name', trim(sdName) )
    call usleep ( delay ) ! Should we make this parallel%delay?
    call AppendL2GPData( l2gp, l2gpFile, &
      & sdName, offset, lastprofile=lastInstance, &
      & TotNumProfs=TotalProfs, createSwath=createSwath )
    if ( L2GPFile%access == DFACC_CREATE ) L2GPFile%access = DFACC_RDWR
    call writeAPrioriAttributes( l2gpFile, dontreplace=.true. )
    if ( switchDetail(switches, 'l2gp') > -1 ) call dump(l2gp)
    ! Clear up our temporary l2gp
    call DestroyL2GPContents(l2gp)
  end subroutine DirectWrite_L2GP_MF

  ! ------------------------------------------- DirectWriteVector_L2Aux_MF --------
  subroutine DirectWriteVector_L2Aux_MF ( L2AUXFile, Vector, &
    & chunkNo, chunks, FWModelConfig, &
    & lowerOverlap, upperOverlap, single, options )

    ! Purpose:
    ! Write plain hdf-formatted files ala l2aux for datasets that
    ! are too big to keep all chunks stored in memory
    ! so instead write them out chunk-by-chunk
    
    ! Despite the name the routine takes vector quantities, not l2aux ones
    ! It dooes so an entrire vector's worth of vector quantities
    use CHUNKS_M, only: MLSCHUNK_T
    use FORWARDMODELCONFIG, only: FORWARDMODELCONFIG_T
    use MLSSTRINGS, only: WRITEINTSTOCHARS
    ! Args:
    type(ForwardModelConfig_T), dimension(:), pointer :: FWModelConfig
    type (Vector_T), intent(in)   :: VECTOR
    type(MLSFile_T)               :: L2AUXFile
    integer, intent(in)           :: CHUNKNO      ! Index into chunks
    type (MLSChunk_T), dimension(:), intent(in) :: CHUNKS
    logical, intent(in), optional :: lowerOverlap
    logical, intent(in), optional :: upperOverlap
    logical, intent(in), optional :: single       ! Write only the 1st instance
    character(len=*), intent(in), optional :: options
    ! Local parameters
    integer                       :: j
    logical                       :: nameQtyByTemplate
    type (VectorValue_T), pointer :: QUANTITY
    type (VectorValue_T), pointer :: PRECISION
    character(len=32)             :: SDNAME       ! Name of sd in output file
    logical :: verbose
    ! Executable
    verbose = ( switchDetail(switches, 'direct') > -1 )
    nameQtyByTemplate = .true.
    if ( present(options) ) nameQtyByTemplate = &
      & .not. ( index(options, 'num') > 0 )
    nullify(precision)
    if ( verbose ) then
      if ( vector%name > 0 ) then
        call get_string( vector%name, sdName )
      else
        sdname = '(unknown)'
      endif
      call outputNamedValue( 'DW L2AUX vector name', trim(sdName) )
    endif
    do j = 1, size(vector%quantities)
      quantity => vector%quantities(j)
      if ( nameQtyByTemplate ) then
        call get_string( quantity%template%name, sdname )
      else
        call writeIntsToChars ( j, sdName )
        sdName = 'Quantity ' // trim(sdName)
      endif
      call DirectWrite_L2Aux_MF ( L2AUXFile, quantity, precision, sdName, &
        & chunkNo, chunks, FWModelConfig, &
        & lowerOverlap, upperOverlap, single, options )
    enddo
  end subroutine DirectWriteVector_L2Aux_MF


  ! ------------------------------------------- DirectWrite_L2Aux_MF --------
  subroutine DirectWrite_L2Aux_MF ( L2AUXFile, quantity, precision, sdName, &
    & chunkNo, chunks, FWModelConfig, &
    & lowerOverlap, upperOverlap, single, options )

    ! Purpose:
    ! Write plain hdf-formatted files ala l2aux for datasets that
    ! are too big to keep all chunks stored in memory
    ! so instead write them out chunk-by-chunk
    
    ! Despite the name the routine takes vector quantities, not l2aux ones
    use CHUNKS_M, only: MLSCHUNK_T
    use CHUNKDIVIDE_M, only: CHUNKDIVIDECONFIG
    use FORWARDMODELCONFIG, only: FORWARDMODELCONFIG_T
    use HDF, only: DFACC_RDonly
    use MLSFILES, only: HDFVERSION_4, HDFVERSION_5, &
      & MLS_CLOSEFILE, MLS_OPENFILE

    ! Args:
    type(ForwardModelConfig_T), dimension(:), pointer :: FWModelConfig
    type (VectorValue_T), intent(in) :: QUANTITY
    type (VectorValue_T), pointer :: PRECISION
    character(len=*), intent(in) :: SDNAME       ! Name of sd in output file
    type(MLSFile_T)                :: L2AUXFile
    integer, intent(in) :: CHUNKNO      ! Index into chunks
    type (MLSChunk_T), dimension(:), intent(in) :: CHUNKS
    logical, intent(in), optional :: lowerOverlap
    logical, intent(in), optional :: upperOverlap
    logical, intent(in), optional :: single       ! Write only the 1st instance
    character(len=*), intent(in), optional :: options
    ! Local parameters
    logical :: alreadyOpen
    logical, parameter :: DEEBUG = .false.
    logical :: deebughere
    integer :: lastMAF
    integer, parameter :: MAXFILES = 100             ! Set for an internal array
    integer :: returnStatus
    character(len=*), parameter :: sdDebug = "R1A:118.B1F:PT.S0.FB25-1 Core"
    logical :: verbose
    ! Executable
    verbose = ( switchDetail(switches, 'direct') > -1 )

    alreadyOpen = L2AUXFile%stillOpen
    deebughere = ( deebug .or. sdname == sdDebug ) .and. .false.
    lastMAF = (quantity%template%instanceOffset+quantity%template%noInstances - &
        & quantity%template%noInstancesLowerOverlap - &
        & quantity%template%noInstancesUpperOverlap)
    if ( chunkNo == size(chunks) .and. MAYWRITEPOSTOVERLAPS .and. &
      & ChunkDivideConfig%allowPostOverlaps ) &
      & lastMAF = quantity%template%instanceOffset + &
      & quantity%template%noInstances - quantity%template%noInstancesLowerOverlap
    if ( .not. alreadyOpen ) then
      call mls_openFile(L2AUXFile, returnStatus)
      if ( returnStatus /= 0 ) &
        call MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to open l2aux file', MLSFile=L2AUXFile)
    endif
    if ( L2AUXFile%access == DFACC_RDONLY )  &
      & call MLSMessage(MLSMSG_Error, ModuleName, &
      & 'l2aux file is rdonly', MLSFile=L2AUXFile)
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
          & 'last profile > grandTotalInstances for ' // trim(sdName), &
          & MLSFile=L2AUXFile)
      call output('instanceOffset: ', advance='no')
      call output(quantity%template%instanceOffset, advance='yes')
      call output('noInstances: ', advance='no')
      call output(quantity%template%noInstances, advance='yes')
      call output('noInstancesLowerOverlap: ', advance='no')
      call output(quantity%template%noInstancesLowerOverlap, advance='yes')
      call output('noInstancesUpperOverlap: ', advance='no')
      call output(quantity%template%noInstancesUpperOverlap, advance='yes')
      call output('last instance: ', advance='no')
      call output(lastMAF, advance='yes')
      call output('grandTotalInstances: ', advance='no')
      call output(quantity%template%grandTotalInstances, advance='yes')
      if ( lastProfTooBigWarns > MAXNUMWARNS ) &
          & call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Max no. of warnings reached--suppressing further ones')
    endif
    if ( deebughere ) then
      print *, 'Direct Writing to ', trim(L2AUXFile%name)
      print *, 'hdfVersion ', L2AUXFile%hdfVersion
      call output('instanceOffset: ', advance='no')
      call output(quantity%template%instanceOffset, advance='yes')
      call output('noInstances: ', advance='no')
      call output(quantity%template%noInstances, advance='yes')
      call output('noInstancesLowerOverlap: ', advance='no')
      call output(quantity%template%noInstancesLowerOverlap, advance='yes')
      call output('noInstancesUpperOverlap: ', advance='no')
      call output(quantity%template%noInstancesUpperOverlap, advance='yes')
      call output('last profile: ', advance='no')
      call output(lastMAF, advance='yes')
      call output('grandTotalInstances: ', advance='no')
      call output(quantity%template%grandTotalInstances, advance='yes')
    endif
    select case (l2AUXFile%hdfversion)
    case (HDFVERSION_4)
      call DirectWrite_L2Aux_MF_hdf4 ( quantity, sdName, L2AUXFile, &
        & chunkNo, chunks, &
        & lowerOverlap=lowerOverlap, upperOverlap=upperOverlap )
      if ( associated(precision) ) & 
        & call DirectWrite_L2Aux_MF_hdf4 ( precision, &
        & trim(sdName) // 'precision', L2AUXFile, chunkNo, chunks, &
        & lowerOverlap=lowerOverlap, upperOverlap=upperOverlap )
    case (HDFVERSION_5)
      call DirectWrite_L2Aux_MF_hdf5 ( quantity, sdName, L2AUXFile, &
        & chunkNo, chunks, &
        & FWModelConfig, lowerOverlap=lowerOverlap, upperOverlap=upperOverlap, &
        & single=single, options=options )
      if ( associated(precision) ) & 
        & call DirectWrite_L2Aux_MF_hdf5 ( precision, &
        & trim(sdName) // 'precision', L2AUXFile, chunkNo, chunks, &
        & FWModelConfig, lowerOverlap=lowerOverlap, upperOverlap=upperOverlap, &
        & single=single, options=options )
    case default
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unsupported hdfVersion for DirectWrite_L2Aux (currently only 4 or 5)' )
    end select
    if ( verbose ) call outputNamedValue( 'DW L2AUX qty name', trim(sdName) )
    if ( .not. alreadyOpen )  call mls_closeFile(L2AUXFile, returnStatus)
    L2AUXFile%errorCode = returnStatus
    L2AUXFile%lastOperation = 'write'
    if ( switchDetail(switches, 'l2aux') < 0 ) return
    call dump(quantity)
    if ( associated(precision) ) call dump(precision)
  end subroutine DirectWrite_L2Aux_MF

  ! ------------------------------------------ DirectWrite_L2Aux_MF_hdf4 --------
  subroutine DirectWrite_L2Aux_MF_hdf4 ( quantity, sdName, L2AUXFile, &
    & chunkNo, chunks, lowerOverlap, upperOverlap )

    use CHUNKS_M, only: MLSCHUNK_T
    use HDF, only: SFN2INDEX, SFSELECT, SFCREATE, &
      & SFENDACC, DFNT_FLOAT32, SFWDATA_F90
    use INTRINSIC, only: L_NONE
    use MLSKINDS, only: R4, R8
    use MLSFILES, only: HDFVERSION_4

    ! Args:
    type (VectorValue_T), intent(in) :: QUANTITY
    ! integer, intent(in) :: SDNAME       ! Name of sd in output file
    character(len=*), intent(in) :: SDNAME       ! Name of sd in output file
    type(MLSFile_T)                :: L2AUXFile
    integer, intent(in) :: CHUNKNO      ! Index into chunks
    type (MLSChunk_T), dimension(:), intent(in) :: CHUNKS
    logical, intent(in), optional :: lowerOverlap
    logical, intent(in), optional :: upperOverlap

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
    sdIndex = sfn2index ( L2AUXFile%fileID%f_id, trim(sdName) )
    if ( sdIndex == -1 ) then
      lastChunk = chunks(size(chunks))
      sizes(noDims) = lastChunk%lastMAFIndex - lastChunk%noMAFSUpperOverlap + 1
      sizes(noDims-1) = quantity%template%noSurfs
      if ( noDims == 3 ) sizes(1) = quantity%template%noChans
      sdId  = sfCreate ( L2AUXFile%fileID%f_id, trim(sdName), DFNT_FLOAT32, &
        & noDims, sizes )
    else
      sdId = sfSelect ( L2AUXFile%fileID%f_id, sdIndex )
    end if
    if ( sdId == -1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
     & 'Error accessing SD '//trim(sdName) // ' (hdf4)', MLSFile=L2AUXFile)

    ! What exactly will be our contribution
    stride = 1
    start = 0
    sizes = 1
    sizes(noDims) = quantity%template%noInstances - &
      & quantity%template%noInstancesLowerOverlap - &
      & quantity%template%noInstancesUpperOverlap
    sizes(noDims-1) = quantity%template%noSurfs
    if ( noDims == 3 ) sizes(1) = quantity%template%noChans
    start(noDims) = quantity%template%instanceOffset

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
        & 'Error ending access to direct write sd (hdf4)', MLSFile=L2AUXFile )

  end subroutine DirectWrite_L2Aux_MF_hdf4

  ! ------------------------------------------ DirectWrite_L2Aux_MF_hdf5 --------
  subroutine DirectWrite_L2Aux_MF_hdf5 ( quantity, sdName, L2AUXFile, &
    & chunkNo, chunks, FWModelConfig, &
    & lowerOverlap, upperOverlap, single, options )

    use CHUNKS_M, only: MLSCHUNK_T
    use CHUNKDIVIDE_M, only: CHUNKDIVIDECONFIG
    use FORWARDMODELCONFIG, only: FORWARDMODELCONFIG_T
    use FORWARDMODELSUPPORT, only: SHOWFWDMODELNAMES
    use HDF5, only: H5GCLOSE_F, H5GOPEN_F
    use INTRINSIC, only: L_NONE
    use L2AUXDATA, only:  L2AUXDATA_T, PHASENAMEATTRIBUTES, &
      & DESTROYL2AUXCONTENTS, &
      & SETUPNEWL2AUXRECORD, WRITEL2AUXATTRIBUTES
    use MLSFILES, only: HDFVERSION_5, DUMP
    use MLSHDF5, only: ISHDF5ATTRIBUTEPRESENT, ISHDF5DSPRESENT, &
      & MAKEHDF5ATTRIBUTE, SAVEASHDF5DS
    use MLSL2TIMINGS, only: SHOWTIMINGNAMES
    use PCFHDR, only: H5_WRITEMLSFILEATTR, H5_WRITEGLOBALATTR
    use QUANTITYTEMPLATES, only: WRITEATTRIBUTES

    ! Args:
    type (VectorValue_T), intent(in) :: QUANTITY
    character(len=*), intent(in) :: SDNAME       ! Name of sd in output file
    type(MLSFile_T)                :: L2AUXFile
    integer, intent(in) :: CHUNKNO      ! Index into chunks
    type (MLSChunk_T), dimension(:), intent(in) :: CHUNKS
    type(ForwardModelConfig_T), dimension(:), pointer :: FWModelConfig
    logical, intent(in), optional :: lowerOverlap
    logical, intent(in), optional :: upperOverlap
    logical, intent(in), optional :: single       ! Write 1st instance only
    character(len=*), intent(in), optional :: options

    ! Local parameters
    integer, parameter :: MAXFILES = 100             ! Set for an internal array
    integer, parameter :: HDFVERSION = HDFVERSION_5

    ! Local variables
    logical :: addQtyAttributes
    logical :: already_there
    ! logical :: attributes_there
    integer :: first_maf
    integer :: grp_id
    type (L2AUXData_T) :: l2aux
    integer :: last_maf
    type ( MLSChunk_T ) :: LASTCHUNK    ! The last chunk in the file
    logical :: mySingle
    integer :: NODIMS                   ! Also index of maf dimension
    integer :: Num_qty_values
    character(len=8) :: overlaps        ! 'lower', 'upper', or 'none'
    integer :: returnStatus
    integer :: SIZES(3)                 ! HDF array sizes
    integer :: START(3)                 ! HDF array starting position
    integer :: STRIDE(3)                ! HDF array stride
    integer :: total_DS_size
    logical, parameter :: MAYCOLLAPSEDIMS = .false.
    ! logical, parameter :: DEEBUG = .true.

    ! executable code
    Num_qty_values = size(quantity%values, 1)*size(quantity%values, 2)
    addQtyAttributes = .false.
    if ( present(options) ) addQtyAttributes = ( index(options, 'A') > 0 )
    if ( quantity%template%frequencyCoordinate == L_None &
      & .and. MAYCOLLAPSEDIMS) then
      noDims = 2
    else
      noDims = 3
    end if
    
    overlaps = 'none'
    if ( present(lowerOverlap) ) then
      if ( lowerOverlap ) overlaps = 'lower'
    endif
    if ( present(upperOverlap) ) then
      if ( upperOverlap ) overlaps = 'upper'
    endif
    
    mySingle = .false.
    if ( present(single) ) mySingle = single

    ! Create or access the SD
    already_there = IsHDF5DSPresent(L2AUXFile%fileID%f_id, trim(sdName))
    if ( .not. already_there ) then
      lastChunk = chunks(size(chunks))
      sizes(noDims) = lastChunk%lastMAFIndex - lastChunk%noMAFSUpperOverlap + 1
      if ( MAYWRITEPOSTOVERLAPS .and. ChunkDivideConfig%allowPostOverlaps ) &
        & sizes(noDims) = lastChunk%lastMAFIndex + 1
      sizes(noDims-1) = quantity%template%noSurfs
      if ( noDims == 3 ) sizes(1) = quantity%template%noChans
    end if

    ! What exactly will be our contribution
    stride = 1
    start = 0
    sizes = 1
    sizes(noDims-1) = quantity%template%noSurfs
    if ( noDims == 3 ) sizes(1) = quantity%template%noChans
    start(noDims) = quantity%template%instanceOffset
    select case ( overlaps )
    case ( 'lower' )
      sizes(noDims) = quantity%template%noInstancesLowerOverlap
      first_maf = 1
      last_maf = quantity%template%noInstances
    case ( 'upper' )
      sizes(noDims) = quantity%template%noInstancesUpperOverlap
      first_maf = quantity%template%noInstances &
        &       - quantity%template%noInstancesUpperOverlap + 1
      last_maf = quantity%template%noInstances
    case ( 'none' )
      sizes(noDims) = quantity%template%noInstances - &
        & quantity%template%noInstancesLowerOverlap - &
        & quantity%template%noInstancesUpperOverlap
      if ( MAYWRITEPOSTOVERLAPS .and. ChunkDivideConfig%allowPostOverlaps .and. &
        & chunkNo == size(chunks) ) &
        & sizes(noDims) = quantity%template%noInstances - &
        & quantity%template%noInstancesLowerOverlap
      first_maf = 1+quantity%template%noInstancesLowerOverlap
      last_maf = quantity%template%noInstances &
        &       - quantity%template%noInstancesUpperOverlap
      if ( MAYWRITEPOSTOVERLAPS .and. ChunkDivideConfig%allowPostOverlaps .and. &
        & chunkNo == size(chunks) ) &
        & last_maf = quantity%template%noInstances
    end select
    
    if ( mySingle ) then
      ! In case we write only the 1st instance
      last_maf = first_maf
      sizes(noDims) = 1
    endif

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
        & 'Number of 3d array elements to write > number stored in qty values', &
        & MLSFile=L2AUXFile )
      call SaveAsHDF5DS( L2AUXFile%fileID%f_id, trim(sdName), &
        & real( &
        &   reshape(quantity%values(:,first_maf:last_maf), sizes(1:3)) &
        & ), start, sizes, may_add_to=.true., adding_to=already_there, &
        & fillValue=DEFAULTUNDEFINEDVALUE )
    else
      total_DS_size = sizes(1)*sizes(2)
      if ( DEEBUG ) then
        print *, 'total_DS_size ', total_DS_size
        print *, 'Num_qty_values ', Num_qty_values
      endif
      if ( total_DS_size > Num_qty_values ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Number of 2d array elements to write > number stored in qty values', &
        & MLSFile=L2AUXFile )
      call SaveAsHDF5DS( L2AUXFile%fileID%f_id, trim(sdName), &
        & real( &
        &   reshape(quantity%values(:,first_maf:last_maf), sizes(1:2)) &
        & ), start, sizes, may_add_to=.true., adding_to=already_there, &
        & fillValue=DEFAULTUNDEFINEDVALUE)
    endif

    ! call mls_CloseFile(L2AUXFILE)
    ! attributes_there = IsHDF5AttributeInFile( L2AUXFile%name, 'Phase Names' )
    ! call mls_OpenFile(L2AUXFILE)
    ! attributes_there = .false.
    ! Now some attribute stuff
    if ( PHASENAMEATTRIBUTES ) then
      if ( DeeBUG ) call dump(L2AUXFile, details=1)
      call h5gopen_f(L2AUXFile%fileID%f_id, '/', grp_id, returnstatus)
      if ( returnstatus /= 0 ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open group "/" to write attribute', &
        & MLSFile=L2AUXFile )
      if ( associated ( FWModelConfig ) ) then
        call MakeHDF5Attribute(grp_id, &
          & 'ForwardModel Names', ShowFwdModelNames(FWModelConfig), &
          & skip_if_already_there=.false., dont_trim=.true.)
      end if
      call MakeHDF5Attribute(grp_id, &
        & 'Phase Names', showTimingNames('phases', .true.), &
        & skip_if_already_there=.false., dont_trim=.true.)
      call h5gclose_f(grp_id, returnstatus)
    endif

    if ( already_there ) return
    ! sd-level attributes
    call SetupNewL2AUXRecord ( L2AUX, quantity%template, &
      & first_MAF, last_MAF-first_MAF+1 )
    if ( DEEBUG ) then
      call output('Writing attributes to: ', advance='no')
      call output(trim(sdName), advance='yes')
    endif
    call WriteL2AUXAttributes(L2AUXFile%fileID%f_id, l2aux, trim(sdName))
    if ( addQtyAttributes ) &
      & call WriteAttributes ( L2AUXFile%fileID%f_id, trim(sdName), quantity%template )
    ! Deallocate memory used by the l2aux
    call DestroyL2AUXContents ( l2aux )
    ! file-level attributes
    call h5_writeglobalattr(L2AUXFile%fileID%f_id, skip_if_already_there=.true.)
    if ( PHASENAMEATTRIBUTES ) then
      call h5gopen_f(L2AUXFile%fileID%f_id, '/', grp_id, returnstatus)
      if ( .not. &
        & IsHDF5AttributePresent('/', L2AUXFile%fileID%f_id, 'Section Names') ) &
        & call MakeHDF5Attribute(grp_id, &
        & 'Section Names', trim(showTimingNames('sections', .true.)), .true.)
      call h5gclose_f(grp_id, returnstatus)
    endif
    call h5_writeMLSFileAttr(L2AUXFile, skip_if_already_there=.true.)

  end subroutine DirectWrite_L2Aux_MF_hdf5

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
    if ( .not. associated(directDB) ) then
      call output ( '**** directWrite Database empty ****', advance='yes' )
      return
    elseif ( size(directDB) < 1 ) then
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

    call output ( 'File Name index: ')
    call output ( directWrite%fileIndex )
    call output ( '   File Name (base): ')
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
    elseif ( directWrite%type == l_l2dgg ) then
      call output ( '(l2dgg)', advance='yes')
    elseif ( directWrite%type == l_l2fwm ) then
      call output ( '(l2fwm)', advance='yes')
    else
      call output ( '(unknown)', advance='yes')
    endif
    if ( directWrite%autoType < 1 ) then
      call output ( '(Is not', advance='no')
    else
      call output ( '(Is', advance='no')
    endif
    call output (' eligible to be auto-filled)', advance='yes' )
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
    type (DirectData_T)                      :: tempDirectData
    ! Begin executable
    isNew = .true.
    ! Do we have any database yet?
    if ( size(directDB) < 1 ) then
      ! print *, 'Setting up initial database with ', trim(fileNameBase)
      call SetupNewDirect(tempDirectData, 0)
      dbID = AddDirectToDatabase( directDB, tempDirectData )
      directData => directDB(dbID)
    else
      ! Check if the fileNameBase already there
      ! print *, 'Checking if ', trim(fileNameBase), ' is in database'
      do dbID = 1, size(directDB)
        ! print *, dbID, trim(directDB(dbID)%fileNameBase)
      enddo
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
    integer :: i
    integer :: newsize
    character(len=80), dimension(:), pointer :: sdNames => null()
    integer :: status
    ! Check if the sdName already there
    ! print *, 'Check if the sdName already there'
    if ( associated(directData%sdNames) ) then
      newSize=SIZE(directData%sdNames)+1
      ! alreadyThere = any(sdName == directData%sdNames)
      alreadyThere = .false.
      do i=1, size(directData%sdNames)
        alreadyThere = alreadyThere .or. (sdName == directData%sdNames(i))
      enddo
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

  !------------------------------------------  FileNameToID  -----
  function FileNameToID ( fileName, DataBase )  result(ID)

    ! Given filename, returns index; if name not found in db, returns 0

    ! Dummy arguments
    type (DirectData_T), dimension(:), pointer :: DATABASE
    character(len=*), intent(in) :: FileName
    integer                      :: ID

    ! Local variables
    ! Executable
    id = 0
    if ( .not. associated(dataBase) .or. len_trim(filename) < 1 ) return
    do id =1, size(database)
      if ( fileName == dataBase(id)%fileName ) exit
      if ( fileName == dataBase(id)%fileNameBase ) exit
    enddo
    if ( id > size(database) ) id = 0
  end function FileNameToID

  !------------------------------------------  SetupNewDirect  -----
  subroutine SetupNewDirect ( directData, NsdNames )

    ! This routine sets up the arrays for a directWrite

    ! Dummy arguments
    type (DirectData_T), intent(inout)  :: directData
    integer, intent(in) :: NsdNames            ! Dimensions

    ! Allocate the sdNames
    directData%autoType = 0
    directData%Handle = 0
    directData%Type = 0
    directData%NSDNames = NSDNames
    directData%FileNameBase = ' '
    directData%FileName = ' '
    nullify(directData%sdNames)
    if ( NSDNames <= 0 ) return
    call allocate_test ( directData%sdNames, NsdNames, "directData%sdNames", &
         & ModuleName )
  end subroutine SetupNewDirect

! =====     Private Procedures     =====================================
  ! ---------------------------------------------  vectorValue_to_l2gp  -----
  subroutine vectorValue_to_l2gp (QUANTITY, &
    & precision, quality, status, convergence,  l2gp, &
    & name, chunkNo, HGrids, offset, firstInstance, lastInstance)
    use HGRIDSDATABASE, only: HGRID_T
    use INTRINSIC, only: L_NONE
    use L2GPDATA, only: L2GPDATA_T, RGP, &
      & SETUPNEWL2GPRECORD
    ! Args:
    type (VectorValue_T), intent(in) :: QUANTITY
    type (VectorValue_T), pointer :: precision
    type (VectorValue_T), pointer :: quality
    type (VectorValue_T), pointer :: status
    type (VectorValue_T), pointer :: convergence
    type (L2GPData_T)                :: l2gp
    character(len=*), intent(in)     :: name
    integer, intent(in)              :: chunkNo
    type (HGrid_T), dimension(:), pointer ::     HGrids
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
    if ( DEEBUG ) print *, 'noFreqsInL2GP, noSurfsInL2GP, lastProfile: ', noFreqsInL2GP, noSurfsInL2GP, lastProfile
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
    L2GP%nTimesTotal = quantity%template%grandTotalInstances

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
    if (associated(precision)) then
      l2gp%l2gpPrecision(:,:,firstProfile:lastProfile) = &
        & reshape ( max ( -hugeRgp, min ( hugeRgp, &
        &   precision%values(:,useFirstInstance:useLastInstance) ) ), &
        &  (/max(l2gp%nFreqs,1),max(l2gp%nLevels,1),lastProfile-firstProfile+1/))
    else
      l2gp%l2gpPrecision(:,:,firstProfile:lastProfile) = 0.0
    end if
    if (associated(quality)) then
      l2gp%quality(firstProfile:lastProfile) = &
        & quality%values(1,useFirstInstance:useLastInstance)
    else
      l2gp%quality(firstProfile:lastProfile) = 0.0
    endif
    if (associated(status)) then
      l2gp%status(firstProfile:lastProfile) = &
        & status%values(1,useFirstInstance:useLastInstance)
    else
      l2gp%status(firstProfile:lastProfile) = 0
    endif
    if (associated(convergence)) then
      l2gp%convergence(firstProfile:lastProfile) = &
        & convergence%values(1,useFirstInstance:useLastInstance)
    else
      l2gp%convergence(firstProfile:lastProfile) = 0.0
    endif
    if ( DEEBUG ) print *, 'Vector converted to l2gp; name: ', trim(name)
    if ( DEEBUG ) print *, 'firstProfile, lastProfile: ', firstProfile, lastProfile
    if ( DEEBUG ) print *, 'useFirstInstance, useLastInstance: ', useFirstInstance, useLastInstance
  end subroutine vectorValue_to_l2gp

  ! ---------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( Where, Full_message, Code, Penalty )

    use LEXER_CORE, only: PRINT_SOURCE
    use TREE, only: SOURCE_REF

    ! Args:
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

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module DirectWrite_m

! $Log$
! Revision 2.56  2013/08/12 23:49:41  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.55  2013/05/08 20:17:35  pwagner
! Changed sleep to usleep to link with NAG
!
! Revision 2.54  2013/05/07 22:23:30  pwagner
! A workaround to prevent unexpected Missing/Fill values
!
! Revision 2.53  2012/05/10 00:48:19  pwagner
! A little less clumsy at creating files, swaths
!
! Revision 2.52  2012/03/12 17:09:44  pwagner
! Use new writeAPrioriAttributes api
!
! Revision 2.51  2012/02/24 21:19:15  pwagner
! May DirectWrite a /single instance only
!
! Revision 2.50  2011/12/15 01:52:55  pwagner
! 'num' option names quantity in entireVector DirectWrites numerically
!
! Revision 2.49  2011/11/01 21:44:31  pwagner
! May optionally add quantity Attributes via options='A' field
!
! Revision 2.48  2011/05/09 18:06:59  pwagner
! May directwrite an entire vector
!
! Revision 2.47  2010/02/04 23:12:44  vsnyder
! Remove USE or declaration for unreferenced names
!
! Revision 2.46  2010/02/04 19:10:26  pwagner
! Change condition on directwrite db from associated to zero size
!
! Revision 2.45  2010/01/08 00:11:19  pwagner
! Added ability to write MLSFile_T fields as attributes
!
! Revision 2.44  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.43  2009/04/23 23:02:52  pwagner
! May specify upperOverlap or lowerOverlap in DirectWrites
!
! Revision 2.42  2008/09/20 00:03:00  pwagner
! Added print statement to not_used_here
!
! Revision 2.41  2006/10/11 22:58:00  pwagner
! Will write convergence ratio as another quality-like field of l2gp
!
! Revision 2.40  2006/06/20 00:13:34  pwagner
! Fwm files should haves MAFs matching radiance files
!
! Revision 2.39  2006/05/11 19:39:19  pwagner
! Should not segment fault if dumping empty DB
!
! Revision 2.38  2006/04/20 23:24:12  pwagner
! More bugs squashed related to extra-range MAFs; one crashed final chunk
!
! Revision 2.37  2006/04/19 20:48:13  pwagner
! Undid most of the changes regarding extra MAFs; perhaps fixed bugs
!
! Revision 2.36  2006/04/12 22:20:42  pwagner
! Attempted workaround for hdf5-1.6.5 bugs rewriting string attributes
!
! Revision 2.35  2006/04/11 23:33:14  pwagner
! Fixed bug which added excess profiles
!
! Revision 2.34  2006/01/26 00:34:50  pwagner
! demoted more use statements from module level to speed Lahey compiles
!
! Revision 2.33  2005/08/19 23:28:02  pwagner
! Trying to avoid possibility of Lahey-causes memory leak
!
! Revision 2.32  2005/07/12 17:38:57  pwagner
! Writes APriori File names as an attribute to every l2gp file
!
! Revision 2.31  2005/06/22 18:57:01  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.30  2005/06/14 20:41:17  pwagner
! Interfaces changed to accept MLSFile_T args
!
! Revision 2.29  2004/11/29 21:52:41  livesey
! Bug fix for handling cases where no forward models are defined.
!
! Revision 2.28  2004/08/03 18:02:01  pwagner
! Sets fillValue for l2aux type
!
! Revision 2.27  2004/07/22 20:42:57  cvuu
! May write ForwardModel Names as file-level attributes and fix the phase names
!
! Revision 2.26  2004/06/29 18:05:26  pwagner
! May write phase, section names as file-level attributes
!
! Revision 2.25  2004/06/10 00:58:45  vsnyder
! Move FindFirst, FindNext from MLSCommon to MLSSets
!
! Revision 2.24  2004/05/19 19:16:09  vsnyder
! Move MLSChunk_t to Chunks_m
!
! Revision 2.23  2004/05/05 21:31:48  pwagner
! More debug printing
!
! Revision 2.22  2004/03/03 19:25:45  pwagner
! Fixed poorly understood prob with single chunks; more initializing in setup
!
! Revision 2.21  2004/02/19 23:54:26  pwagner
! Dumps during directWrites if l2gp or l2aux switch set
!
! Revision 2.20  2004/02/11 17:23:25  pwagner
! l2gp status an integer, not a char
!
! Revision 2.19  2004/02/10 19:30:55  pwagner
! Cures serial directWrites from writing more than one chunk at a time
!
! Revision 2.18  2004/02/05 23:38:41  pwagner
! Writes attributes to directwrite l2aux file
!
! Revision 2.17  2004/01/23 01:09:48  pwagner
! Only directwrite files entered in global settings eligible to be auto-sourced
!
! Revision 2.16  2004/01/22 06:38:26  livesey
! Typo fix
!
! Revision 2.15  2004/01/22 00:56:35  pwagner
! Fixed many bugs in auto-distribution of DirectWrites
!
! Revision 2.14  2003/12/03 17:50:54  pwagner
! L2GP tracks both nTimes (for this slave) and nTimesTotal (done by all)
!
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
