! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

! -------------------------------------------------------
module WriteMetadata ! Populate metadata and write it out
! -------------------------------------------------------

  use Dump_0, only: Dump
  use HDF, only: Dfacc_Rdwr
  use HighOutput, only: Beverbose, OutputnamedValue
  use Init_Tables_Module, only: L_L2dgg, L_L2gp, L_HDF, L_Swath
  use Io_Stuff, only: Read_TextFile
  use Lexer_Core, only: Print_Source
  use MLSCommon, only: Filenamelen, Namelen, L2metaData_T, MLSFile_T
  use MLSKinds, only: R8
  use MLSFiles, only: Dump, Getpcfromref, MLS_CloseFile, MLS_OpenFile, &
    & MLS_Sfstart, MLS_Sfend, Split_Path_Name
  use MLSL2options, only: Toolkit, Sharedpcf
  use MLSMessagemodule, only: MLSMessage, MLSMSG_Warning
  use MLSPcf2, only: MLSPcf_Mcf_L2gp_End, MLSPcf_Mcf_L2gp_Start, &
    & MLSPcf_Mcf_L2log_Start
  use MLSStrings, only: Iscomment, Lowercase, Streq
  use MLSStringlists, only: Extractsubstring, &
    & GethashElement, GetstringElement
  use Output_M, only: Output, Blanks
  use Pcfhdr, only: Writeinputpointer, Writepcf2hdr, Globalattributes
  use Sdptoolkit, only: Pgsd_Met_Group_Name_L, &
    & Pgsd_Met_Num_Of_Groups, Pgsd_Pc_File_Path_Max, Pgs_Pc_Getreference, &
    & Pgspc_W_No_Reference_Found, Pgs_S_Success, Pgsmet_W_MetaData_Not_Set
  use Tree, only: Where

  implicit none

  private

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (parameters and datatypes)
! MCGROUP_T          metadata control group
! MCPARAM_T          metadata control object; e.g. PGEVersion
! PCFData_T          metadata-relevant data taken from PCF

!     (subroutines and functions)
! addMCGroupToDatabase  Adds a new MCGroup to existing database of them
! addMCParamToDatabase  Adds a new MC Parameter to existing database of them
! addToMetadata         Adds a new metadata field onto an existing file
! destroyMCGroup        Deallocates any arrays in a MCGroup
! destroyMCGroupDB      Destroys each MCGroup in a database of them, 
!                        then deallocates it, too
! dumpMCGroup           Dumps all the parameters in a MCGroup
! dumpMCGroupDB         dumps each MCGroup in a database of them
! get_l2gp_mcf          Finds the mcf relevant for a given species based on
!                        its file name 
! mctree                parse list of mc odl-like character args
! NullifyPCFData        Nullify all pointers associated with a PCFData_T
! Populate_metadata_std Write metadata to a standard l2 product file, e.g. H2O
! Populate_metadata_oth Write metadata to other l2 product files, e.g. DGG, DGM
! readMCF               Read a metadata control file (mcf) into a database of
!                         metadata control groups
! writemetalog          Write metadata for the log file
! writeMCF              Write a database of metadata control groups to
!                         a metadata control file (mcf)
! === (end of toc) ===

! === (start of api) ===
! addMCGroupToDatabase (type MCGroup_T dataBase(:), type MCGroup_T item)
! addMCParamToDatabase (type MCParam_T dataBase(:), type MCParam_T item)
! addToMetadata (char* mcf_grp, char* mcf_grp, value)
!      where value can be a single, double, integer, or string
! DestroyMCGroup (type MCGroup_T MCGroup)
! DestroyMCGroupDB (type MCGroup_T MCGroupDB(:))
! dumpMCGroup (type MCGroup_T MCGroup)
! dumpMCGroupDB (type MCGroup_T MCGroupDB(:))
! Get_l2gp_mcf (char* file_base, char* meta_name, int mcf, int version)
! mctree ( char* list(:), type MCGroup mcdb*(:) )
! NullifyPCFData (type PCFData_T P)
! Populate_metadata_std ( int hdf_file, int mcf_file, &
!   & char* Field_name, type(L2Metadata_T) l2metaData, &
!   & [int hdfVersion], [int Metadata_error], &
!   & [int filetype] )
! Populate_metadata_oth ( int hdf_file, int mcf_file, &
!   & int NumQuantitiesPerFile, char* QuantityNames, &
!   & type(L2Metadata_T) l2metaData, &
!   & [int hdfVersion], [int Metadata_error], &
!   & [int filetype] )
! readMCF (type MLSFile_T MLSFile, type MCGroup_T MCGroups(:))
! writeMCF (type MLSFile_T MLSFile, type MCGroup_T MCGroups(:))
! WriteMetaLog ([int metadata_error])
! === (end of api) ===

  public :: AddmcgrouptoDatabase, AddmcparamtoDatabase, AddtometaData, &
    & Destroy, Destroymcgroup, Destroymcgroupdb, &
    & Dump, Dumpmcgroup, Dumpmcgroupdb, &
    & Get_L2gp_Mcf, Mctree, NullifypcfData, &
    & Populate_MetaData_Std, Populate_MetaData_Oth, &
    & Readmcf, Writemcf, Writemetalog


! This data type is used to store User-defined Runtime Parameters and other
! information taken from the PCF.
! Why are these len=80 instead of, say, len=NameLen?
! If they need to be 80, then shouldn't NameLen be 80, too?
  type, public :: MCGROUP_T
    character (len=80) :: name = ' '
    character (len=80) :: class         = ' '
    character (len=80) :: type = ' '
    type(MCGROUP_T), dimension(:), pointer :: groups => null()
    type(MCPARAM_T), dimension(:), pointer :: params => null()
  end type MCGROUP_T

  type, public :: MCPARAM_T
    character (len=80) :: name          = ' '
    character (len=80) :: mandatory     = ' '
    character (len=80) :: data_location = ' '
    character (len=80) :: class         = ' '
    character (len=80) :: num_val       = ' '
    character (len=80) :: type          = ' '
    character (len=256) :: value        = ' '
    type(MCGROUP_T), dimension(:), pointer :: groups => null()
    type(MCPARAM_T), dimension(:), pointer :: params => null()
  end type MCPARAM_T

  type, public :: PCFData_T

    ! cycle # of processing run
    character (len=4) :: cycle         ! add to output files

    ! id string processing run
    character (len=16) :: RunID = ' '

    ! version string in PCF output file names
    character (len=15) :: PGEVersion   ! add to output files

    character(len=27) :: StartUTC
    character(len=27) :: EndUTC

    ! The annotation text to be written to the header of every scientific
    ! data file for which we need metadata
    ! In practice, this annotation will be the contents of the PCF file itself

    character (len=1), pointer :: AnText(:) => null()

    !     How to choose the mcf file to be used in writing metadata
    !                          (3)
    ! The SPECIES of the corresponding
    ! MCF file names are chosen from an associative array structure
    ! or hash table, where the keys are the possible species and
    ! the hash
    ! how to associate l2gp species names with mcf files

    ! if an associative array, then l2gp file names contain the "keys"
    ! the following is a comma-delimited list of possible species
    ! names that may be such keys

    ! Note that the text matching in options (2) and (3)
    ! will not be case sensitive unless you set MCFCASESENSITIVE

    !                          (*)
    ! Finally, you can avoid any of this trickery by using the metaName=
    ! field for each file. For example, if you set metaName='o3', the
    ! mcf matching o3 will be selected, no matter what you named the l2gp file.
    character (len=fileNameLen) :: spec_keys

    ! the following is a comma-delimited list of possible
    ! mcf file name parts that may be the corresponding hash values
    character (len=fileNameLen) :: spec_mcfnames      

    ! DOIs
    character (len=20*NameLen) :: spec_doinames   

    ! name of the log file (without the path)
    character (len=fileNameLen) :: logGranID

    ! identifier_product_doi
    character (len=32) :: DOI     ! unique for each type

  end type PCFData_T

  interface addToMetadata
    module procedure addToMetadata_dbl
    module procedure addToMetadata_int
    module procedure addToMetadata_sngl
    module procedure addToMetadata_str
  end interface
  interface DESTROY
    module procedure DESTROYMCGROUP
    module procedure DESTROYMCGROUPDB
  end interface
  interface DUMP
    module procedure DUMPMCPARAM
    module procedure DUMPMCGROUP
    module procedure DUMPMCGROUPDB
  end interface
  integer, public, parameter :: INVENTORYMETADATA = 2
  logical, public, parameter :: MCFCASESENSITIVE = .FALSE.
  logical, public, parameter :: ANNOTATEWITHPCF = .TRUE.
  logical, public, parameter :: SETINPUTPOINTER = .not. ANNOTATEWITHPCF
  logical, public, parameter :: SFINBETWEENSTARTEND = .FALSE.
  integer, public, parameter :: MCFFORL2GPOPTION = 3     ! 1, public, 2 or 3
  integer, private :: Module_error
  integer, private :: nest_degree
  type(PCFData_T), public, save :: L2PCF
  character(len=*), parameter :: spaces = '                                  ' &
    & // '                                                                    '
  logical, private, parameter :: countEmpty = .true.

contains

  !-----------------------------------------  AddMCGroupToDatabase  -----
  integer function AddMCGroupToDatabase ( DATABASE, ITEM )

  ! This routine adds a vector to a database of such vectors, 
  ! creating the database if necessary.

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

    ! Dummy arguments
    type (MCGroup_T), dimension(:), pointer :: DATABASE
    type (MCGroup_T), intent(in) ::            ITEM

    ! Local variables
    type (MCGroup_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddMCGroupToDatabase = newSize
  end function AddMCGroupToDatabase

  !-----------------------------------------  AddMCParamToDatabase  -----
  integer function AddMCParamToDatabase ( DATABASE, ITEM )

  ! This routine adds a vector to a database of such vectors, 
  ! creating the database if necessary.

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

    ! Dummy arguments
    type (MCParam_T), dimension(:), pointer :: DATABASE
    type (MCParam_T), intent(in) ::            ITEM

    ! Local variables
    type (MCParam_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddMCParamToDatabase = newSize
  end function AddMCParamToDatabase

  ! ------------ addToMetadata
  ! This family of routines adds on to existing metadata in a file
  subroutine addToMetadata_dbl( MCF_GRP, ATTRNAME, VALUE )
    ! Args
    character (len = PGSd_MET_GROUP_NAME_L), intent(in) :: MCF_GRP
    character(len=*), intent(in)     :: ATTRNAME
    double precision, intent(in)     :: VALUE
    integer, external :: PGS_MET_SETATTR_D
    ! Local variables
    integer :: status
    ! Executable
    status = pgs_met_setAttr_d ( mcf_grp, trim(attrName), value )
    if ( status /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
        & "Error in writing " // trim(attrName) )
    end if
  end subroutine addToMetadata_dbl

  subroutine addToMetadata_int( MCF_GRP, ATTRNAME, VALUE )
    ! Args
    character (len = PGSd_MET_GROUP_NAME_L), intent(in) :: MCF_GRP
    character(len=*), intent(in)     :: ATTRNAME
    integer, intent(in)     :: VALUE
    integer, external :: PGS_MET_SETATTR_I
    ! Local variables
    integer :: status
    ! Executable
    status = pgs_met_setAttr_i ( mcf_grp, trim(attrName), value )
    if ( status /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
        & "Error in writing " // trim(attrName) )
    end if
  end subroutine addToMetadata_int

  subroutine addToMetadata_sngl( MCF_GRP, ATTRNAME, VALUE )
    ! Args
    character (len = PGSd_MET_GROUP_NAME_L), intent(in) :: MCF_GRP
    character(len=*), intent(in)     :: ATTRNAME
    real, intent(in)                 :: VALUE
    integer, external :: PGS_MET_SETATTR_D
    ! Local variables
    integer :: status
    ! Executable
    status = pgs_met_setAttr_d ( mcf_grp, trim(attrName), value*1.0d0 )
    if ( status /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
        & "Error in writing " // trim(attrName) )
    end if
  end subroutine addToMetadata_sngl

  subroutine addToMetadata_str( MCF_GRP, ATTRNAME, VALUE )
    ! Args
    character (len = PGSd_MET_GROUP_NAME_L), intent(in) :: MCF_GRP
    character(len=*), intent(in)     :: ATTRNAME
    character(len=*), intent(in)     :: VALUE
    integer, external :: PGS_MET_SETATTR_S
    ! Local variables
    integer :: status
    ! Executable
    status = pgs_met_setAttr_s ( mcf_grp, trim(attrName), value )
    if ( status /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
        & "Error in writing " // trim(attrName) )
    end if
  end subroutine addToMetadata_str

  ! ----------------------------------------DestroyMCGroup -----
  recursive subroutine DestroyMCGroup ( MCGROUP )
    ! Destroy a group of inventory metadata control parameters 
    ! If the group contains other groups, destroy them, too
    use Allocate_Deallocate, only: Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type ( MCGROUP_T ), intent(inout)               :: MCGROUP
    ! Internal variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: error, i, s

    ! Executable code
    if ( associated(MCGroup%groups) ) then
      do i=1, size(MCGroup%groups)
        call DestroyMCGroup( MCGroup%groups(i) )
      end do
      s = size(MCGroup%groups) * storage_size(MCGroup%groups) / 8
      addr = 0
      if ( s > 0 ) addr=transfer(c_loc(MCGroup%groups(1)), addr)
      deallocate ( MCGroup%groups, stat=error )
      call test_deallocate ( error, moduleName, 'MCGroup%groups', s, address=addr )
    end if
    if ( associated(MCGroup%params) ) then
      do i=1, size(MCGroup%params)
        call DestroyMCParam( MCGroup%params(i) )
      end do
      s = size(MCGroup%params) * storage_size(MCGroup%params) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(MCGroup%params(1)), addr)
      deallocate ( MCGroup%params, stat=error )
      call test_deallocate ( error, moduleName, 'MCGroup%params', s, address=addr )
    end if
  end subroutine DestroyMCGroup

  ! ----------------------------------------DestroyMCParam -----
  recursive subroutine DestroyMCParam ( MCParam )
    ! Destroy a metadata control parameter
    ! If it contains other params and groups, destroy them, too
    use Allocate_Deallocate, only: Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type ( MCParam_T ), intent(inout)               :: MCParam
    ! Internal variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: error, i, s

    ! Executable code
    if ( associated(MCParam%groups) ) then
      do i=1, size(MCParam%groups)
        call DestroyMCGroup( MCParam%groups(i) )
      end do
      s = size(MCParam%groups) * storage_size(MCParam%groups) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(MCParam%groups(1)), addr)
      deallocate ( MCParam%groups, stat=error )
      call test_deallocate ( error, moduleName, 'MCParam%groups', s, address=addr )
    end if
    if ( associated(MCParam%params) ) then
      do i=1, size(MCParam%params)
        call DestroyMCParam( MCParam%params(i) )
      end do
      s = size(MCParam%params) * storage_size(MCParam%params) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(MCParam%params(1)), addr)
      deallocate ( MCParam%params, stat=error )
      call test_deallocate ( error, moduleName, 'MCParam%params', s, address=addr )
    end if
  end subroutine DestroyMCParam

  ! ----------------------------------------DestroyMCGroupDB -----
  subroutine DestroyMCGroupDB ( MCGROUPDB )
    ! Destroy a database of metadata control groups
    use Allocate_Deallocate, only: Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type ( MCGROUP_T ), dimension(:), pointer         :: MCGROUPDB
    ! Internal variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: error, i, s

    ! Executable code
    if ( .not. associated(MCGROUPDB) ) return
    do i=1, size(MCGroupDB)
      call DestroyMCGroup( MCGroupDB(i) )
    end do
    s = size(MCGroupDB) * storage_size(MCGroupDB) / 8
    addr = 0
    if ( s > 0 ) addr = transfer(c_loc(MCGroupDB(1)), addr)
    deallocate ( MCGroupDB, stat=error )
    call test_deallocate ( error, moduleName, 'MCGroupDB', s, address=addr )
  end subroutine DestroyMCGroupDB

  ! ----------------------------------------dumpMCGroup -----
  recursive subroutine dumpMCGroup ( MCGROUP )
    ! dump a group of inventory metadata control parameters 
    ! If the group contains other groups, dump them, too
    type ( MCGROUP_T ), intent(in)               :: MCGROUP
    ! Internal variables
    integer :: i
    integer :: indent

    ! Executable code
    indent = nest_degree*3
    if ( len_trim(MCGroup%name) > 0 ) &
      & call outputnamedValue( spaces(:indent)//'group name', &
      & MCGroup%name )
    if ( len_trim(MCGroup%class) > 0 ) &
      & call outputnamedValue( spaces(:indent)//'class', &
      & MCGroup%class )
    if ( len_trim(MCGroup%type) > 0 ) &
      & call outputnamedValue( spaces(:indent)//'type', &
      & MCGroup%type )
    nest_degree = nest_degree + 1
    indent = nest_degree*3
    if ( .not. associated(MCGroup%groups) ) then
      call blanks( indent )
      call output( '(No subgroups of this group) ', advance='yes' )
    else
      do i=1, size(MCGroup%groups)
        call blanks( indent )
        call output( 'Sub-group ', advance='no' )
        call output( i, advance='yes' )
        call dumpMCGroup( MCGroup%groups(i) )
      enddo
    endif
    if ( .not. associated(MCGroup%params) ) then
      call blanks( indent )
      call output( '(No parameters in this group) ', advance='yes' )
    else
      do i=1, size(MCGroup%params)
        call dumpMCParam( MCGroup%params(i) )
      enddo
    endif
    nest_degree = nest_degree - 1
  end subroutine dumpMCGroup

  ! ----------------------------------------dumpMCGroupDB -----
  subroutine dumpMCGroupDB ( MCGROUPDB )
    ! dump a database of metadata control groups
    type ( MCGROUP_T ), dimension(:), pointer         :: MCGROUPDB
    ! Internal variables
    integer :: i

    ! Executable code
    if ( .not. associated(MCGROUPDB) ) then
      call output( ' (MCGroup database not associated) ', advance='yes' )
      return
    endif
    call output( 'Dumping MCGroup database with ', advance='no' )
    call output( size(MCGroupDB), advance='no' )
    call output( ' groups ', advance='yes' )
    do i=1, size(MCGroupDB)
      nest_degree = 1
      call dumpMCGroup( MCGroupDB(i) )
    enddo
  end subroutine dumpMCGroupDB

  ! ----------------------------------------dumpMCParam -----
  recursive subroutine dumpMCParam ( MCParam )
    ! dump a group of inventory metadata control parameters 
    ! If the group contains other groups, dump them, too
    type ( MCParam_T ), intent(in)               :: MCParam
    ! Internal variables
    integer :: i
    integer :: indent

    ! Executable code
    indent = nest_degree*3
    if ( len_trim(MCParam%name) > 0 ) &
      & call outputnamedValue( spaces(:indent)//'parameter name', &
      & MCParam%name )
    nest_degree = nest_degree + 1
    indent = nest_degree*3
    if ( .not. associated(MCParam%groups) ) then
      call blanks( indent )
      call output( '(No subgroups of this parameter) ', advance='yes' )
    else
      do i=1, size(MCParam%groups)
        call blanks( indent )
        call output( 'Sub-group ', advance='no' )
        call output( i, advance='yes' )
        call dumpMCGroup( MCParam%groups(i) )
      enddo
    endif
    if ( .not. associated(MCParam%params) ) then
      call blanks( indent )
      call output( '(No parameters in this parameter) ', advance='yes' )
    else
      do i=1, size(MCParam%params)
        call dumpMCParam( MCParam%params(i) )
      enddo
    endif
    indent = nest_degree*3
    indent = indent + 4
    if ( len_trim(MCParam%name) > 0 ) &
      & call outputnamedValue( spaces(:indent)//'name', &
      & MCParam%name )
    if ( len_trim(MCParam%mandatory) > 0 ) &
      & call outputnamedValue( spaces(:indent)//'mandatory', &
      & MCParam%mandatory )
    if ( len_trim(MCParam%data_location) > 0 ) &
      & call outputnamedValue( spaces(:indent)//'data_location', &
      & MCParam%data_location )
    if ( len_trim(MCParam%class) > 0 ) &
      & call outputnamedValue( spaces(:indent)//'class', &
      & MCParam%class )
    if ( len_trim(MCParam%num_val) > 0 ) &
      & call outputnamedValue( spaces(:indent)//'num_val', &
      & MCParam%num_val )
    if ( len_trim(MCParam%type) > 0 ) &
      & call outputnamedValue( spaces(:indent)//'type', &
      & MCParam%type )
    if ( len_trim(MCParam%value) > 0 ) &
      & call outputnamedValue( spaces(:indent)//'value', &
      & MCParam%value )
    nest_degree = nest_degree - 1
  end subroutine dumpMCParam

  ! -----------------------------------------------  Get_l2gp_mcf  -----

  subroutine Get_l2gp_mcf ( File_base, meta_name, Mcf, doiIdentifier )

  ! metadata configuration file (mcf) PCF number corresponding to l2gp number
  ! sdid
  ! Arguments
    character(len=*), intent(in) ::  File_base
    character(len=*), intent(in) ::  meta_name
    integer, intent(inout) ::        Mcf
    ! integer, intent(in), optional :: Version
    character(len=*), intent(out) ::  doiIdentifier

    ! Local
    character (len=PGSd_PC_FILE_PATH_MAX) :: Sd_full
    !character (len=NameLen) :: Sd_path
    character (len=NameLen) :: Sd_name

    character (len=PGSd_PC_FILE_PATH_MAX) :: Mcf_full
    character (len=NameLen) :: Mcf_path
    character (len=NameLen) :: Mcf_name
    integer :: ReturnStatus, MyVersion, I

    character (len=1), parameter :: COMMA = ','

    ! Find species name
    ! assume sd_name is "*l2gp-species_"
    ! hence enclosed between "_" chars after an l2gp

    logical, parameter :: DEBUG = .false.

    ! Begin

    if(.NOT. TOOLKIT) then
      mcf=0
      return
   endif

    ! if ( MCFFORL2GPOPTION == 1 ) then
    !  mcf = mcf+1
    !  return
    ! end if

    if ( DEBUG ) then
      call output('file_base: ', advance='no')
      call output(trim(file_base), advance='yes')
      call output('meta_name: ', advance='no')
      call output(trim(meta_name), advance='yes')
    end if

    if ( len(TRIM(file_base)) <= 0 ) then
      mcf=0
      return
    end if

    if ( meta_name == ' ' ) then
      ! * The following complication is attempting to infer
      ! * the species name based on the file name
      ! * It is necessary only if you disdained to specify the metaName
      ! * field in the output command
      ! Get full file name for typical MCF file
      do i=mlspcf_mcf_l2gp_start, mlspcf_mcf_l2gp_end

        myVersion = 1

        returnStatus = PGS_PC_GetReference(i, myVersion, mcf_full)

        if ( returnStatus == PGS_S_SUCCESS ) then 
          exit
        end if

      end do

      if ( DEBUG ) then
        call output('returnStatus: ', advance='no')
        call output(returnStatus, advance='yes')
      end if

      if ( returnStatus /= PGS_S_SUCCESS ) then 
        mcf = 0
        return
      end if

      ! Split full_file_names into path+name

      call split_path_name ( mcf_full, mcf_path, mcf_name )

      ! If   text matching not case sensitive, shift to lower case
      if ( .NOT. MCFCASESENSITIVE ) then
        sd_full = LowerCase(file_base)
        mcf_name = LowerCase(mcf_name)
      else
        sd_full = file_base
      end if

      if ( DEBUG ) then
        call output('mcf_full: ', advance='no')
        call output(trim(mcf_full), advance='yes')
        call output('mcf_path: ', advance='no')
        call output(trim(mcf_path), advance='yes')
        call output('mcf_name: ', advance='no')
        call output(trim(mcf_name), advance='yes')
      end if

      ! Either we were given a short form of the sd_full
      ! e.g., 'h2o', or else a much longer one like 'mls-aura_...'
      if ( index(sd_full, 'mls-aura_') > 0 ) then
        ! Get species name assuming e.g. 'mls-aura_l2gp-h2O_'
        call ExtractSubString(sd_full, sd_name, 'mls-aura_l2gp-', '_')
      else
        sd_name = sd_full
      endif

      if ( DEBUG ) then
        call output('sd_full: ', advance='no')
        call output(trim(sd_full), advance='yes')
        call output('sd_name: ', advance='no')
        call output(trim(sd_name), advance='yes')
      end if

      if ( len(trim(sd_name)) <= 0 ) then
        mcf=0
        return
      end if

      ! if ( MCFFORL2GPOPTION == 3 ) then

        ! get mcfspecies name from associative array
        ! if the species name not found in spec_keys, it will return ','
        call GetHashElement ( l2pcf%spec_keys, l2pcf%spec_mcfnames      , &
          & trim(sd_name), sd_full, .TRUE. )
        call GetHashElement ( l2pcf%spec_keys, l2pcf%spec_doinames      , &
          & trim(sd_name), doiIdentifier, .TRUE. )

        if ( DEBUG ) then
          call output('keys: ', advance='no')
          call output(trim(l2pcf%spec_keys), advance='yes')
          call output('hash: ', advance='no')
          call output(trim(l2pcf%spec_mcfnames      ), advance='yes')
          call output('hash for species name: ', advance='no')
          call output(trim(sd_full), advance='yes')
        end if

        if ( trim(sd_full) == COMMA ) then
          mcf = 0
          return
        end if

        sd_name = trim(sd_full)

      ! end if
   elseif ( .NOT. MCFCASESENSITIVE ) then
     sd_name = Lowercase(meta_name)
   else
     sd_name = meta_name
   endif
   ! Now try to find mcf file corresponding to species name   
   ! assuming, e.g. '*h2o.*'                                  
   if ( DEBUG ) then                                       
     call output('Trying for a pcf name matching: ', advance='no')                   
     call output(trim(sd_name), advance='yes')     
   end if                                                  
   mcf = GetPCFromRef(trim(sd_name), &
    & mlspcf_mcf_l2gp_start, mlspcf_mcf_l2gp_end, &
    & MCFCASESENSITIVE, returnStatus)

   if ( returnStatus /= 0 ) then
     mcf = 0
   endif

  end subroutine Get_l2gp_mcf

   
  ! -------------------------------------- mctree ----------------------
  ! parse list of mc odl-like commands entered as character strings
  ! each of one of the following form s
  ! (a) key = value
  ! (b) end_group [= value]
  ! (c) end_parameter [= value]
  ! To start a new group, choose the special key
  !   group_name
  ! To start a new parameter, choose the special key
  !   parameter_name
  subroutine mctree ( list, mcdb )
     character(len=*), dimension(:), intent(in)   :: list
     type (MCGroup_T), dimension(:), pointer      :: mcdb
     ! Internal variables
     integer :: i
     character(len=64)  :: key
     type (MCGROUP_T)   :: group
     integer            :: newsize
     type (MCPARAM_T)   :: param
     character(len=128) :: value
     ! Executable
     if( size(list) < 1 ) return
     i = 0
     do
       i = i+1
       if ( i > size(list) ) exit
       
       ! Which form is our list_line?
       if ( index( trim(list(i)), 'end_group' ) > 0 ) then
         newsize = addmcgrouptodatabase( mcdb, group )
       elseif ( any( streq( list(i), &
         & (/ 'end_parameter', 'end_object   ' /) , options='-cfh' ) ) )  then
         newsize = addmcparamtodatabase( group%params, param )
       elseif( isComment(list(i)) ) then
         ! Just a comment
       elseif( len_trim(list(i)) < 1 ) then
         ! Just a blank line
       elseif ( index( trim(list(i)), '=' ) > 0 ) then
         ! separate into key, value
         call GetStringElement ( list(i), key, &
           & 1, countEmpty, '=' )
         call GetStringElement ( list(i), value, &
           & 2, countEmpty, '=' )
         ! Are we beginning a new group or parameter
         if ( any( streq( key, &
           & (/ 'group_name', 'group     ' /) , options='-cfh' ) ) )  then
           call insertGroup ( list, mcdb, i, value )
         elseif ( any( streq( key, &
           & (/ 'parameter_name', 'object        ' /) , options='-cfh' ) ) )  then
           call insertParameter ( list, group%params, i, value )
         else
           call setParam ( param, key, value )
         endif
       else
         ! unrecognized line
         call output( trim(list(i)) // ' is unrecognized', advance='yes' )
       endif
     enddo
   end subroutine mctree

  ! ----------------------------------------NullifyPCFData -----
  subroutine NullifyPCFData ( P )
    ! Given a PCFData, nullify all the pointers associated with it
    type ( PCFData_T ), intent(out) :: P

    ! Executable code
    nullify ( p%AnText )
  end subroutine NullifyPCFData

  ! --------------------------------------  Populate_metadata_std  -----

  subroutine Populate_metadata_std ( HDF_FILE, MCF_FILE, &
    & Field_name, l2metaData, hdfVersion, Metadata_error, &
    & filetype )

    ! This is the standard way to write meta data
    ! It should work unchanged for the standard l2gp files (e.g. BrO)
    ! and, with minor changes, for the l2gp file marked "other"
    !
    ! the l2aux files, also called dgm
    ! and the dgg files will probably require special treatment
    ! the log file should be able to use the one stolen from level 3
    !

    !    USE InitPCFs, ONLY: L2PCF

    !Arguments

    integer                        :: HDF_FILE, MCF_FILE
    character (len=*)              :: Field_name
    type (L2Metadata_T) :: l2metaData
    integer, optional, intent(in)  :: hdfVersion
    integer, optional, intent(out) :: Metadata_error
    integer, optional, intent(in)  :: filetype  ! 'sw' or 'hdf'

    !Local Variables

    integer :: hdfReturn
    integer :: returnStatus
    integer :: sdid, hdf_sdid

    character (len=PGSd_PC_FILE_PATH_MAX) :: physical_filename
    character (len=PGSd_PC_FILE_PATH_MAX) :: mcf_filename
    integer :: version

    ! the group have to be defined as 49 characters long. The C interface is 50.
    ! The cfortran.h mallocs an extra 1 byte for the null character '\0/1, 
    ! therefore making the actual length of a 
    ! string pass of 50.

    character (len = PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)

    ! Externals

    integer, external :: PGS_MET_remove

    !Executable code

    module_error = 0
    if ( present(metadata_error)) metadata_error=1
    if(.NOT. TOOLKIT) return

    version = 1
    if ( mcf_file > 0 ) then
      returnStatus = PGS_PC_GetReference (MCF_FILE, version , mcf_filename)
    else
      returnStatus = PGSPC_W_NO_REFERENCE_FOUND
    end if

    if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
      & "Error: failed to find PCF ref for MCF_FILE in populate_metadata_std.") 
      return
    end if

    version = 1
    returnStatus = PGS_PC_GetReference (HDF_FILE, version , physical_filename)

    if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
      & "Error: failed to find PCF ref for HDF_FILE in populate_metadata_std.") 
      return
    end if

    call first_grouping(HDF_FILE, MCF_FILE, groups)
    call measured_parameter (HDF_FILE, field_name, groups, 1)
    call third_grouping (HDF_FILE, hdf_sdid, groups, &
      & hdfVersion, filetype, l2metaData)

    if ( SFINBETWEENSTARTEND ) then
      sdid = mls_sfstart (physical_fileName, DFACC_RDWR, &
        & hdfVersion=hdfVersion, addingMetaData=.true.)
    else
      sdid = hdf_sdid
    endif 

    if ( sdid == -1 ) then
      call announce_error ( 0, &
      & "Error: failed to open the hdf file in populate_metadata_std: "&
      & //TRIM(physical_fileName)) 
      return
    end if

    hdfReturn = mls_sfend( sdid, hdfVersion=hdfVersion, addingMetadata=.true. )
    if ( hdfReturn /= 0 ) then
        call announce_error ( 0, &
        & "Error: metadata mls_sfend in populate_metadata_std.", &
        & error_number=hdfReturn) 
    end if

    ! Annotate the file with the PCF

    if ( ANNOTATEWITHPCF ) then
      call PCF2Hdr( physical_filename, hdfVersion, filetype=filetype )
    end if

    returnStatus = pgs_met_remove() 

    if ( present(metadata_error)) metadata_error=module_error

   if(BeVerbose( 'pro', -1 )) then
       call announce_success(physical_filename, mcf_filename, 'standard')
   end if

  end subroutine Populate_metadata_std

  ! --------------------------------------  Populate_metadata_oth  -----

  subroutine Populate_metadata_oth ( HDF_FILE, MCF_FILE, &
    & NumQuantitiesPerFile, QuantityNames, l2metaData, hdfVersion, Metadata_error, &
    & filetype )

    ! This is specially to write meta data for heterogeneous files
    ! It should work unchanged for the 'OTH' l2gp files (e.g. ML2OTH.001.MCF)
    ! and, with minor changes, for the l2aux files 

    !Arguments

    integer :: HDF_FILE, MCF_FILE, NumQuantitiesPerFile
    character (len=*), dimension(:) :: QuantityNames
    type (L2Metadata_T) :: l2metaData
    integer, optional, intent(out) :: Metadata_error
    integer, optional, intent(in) :: hdfVersion
    integer, optional, intent(in)  :: filetype  ! 'sw' or 'hdf'

    !Local Variables

    integer :: HdfReturn
    integer :: ReturnStatus
    integer :: Sdid, hdf_sdid

    character (len=PGSd_PC_FILE_PATH_MAX) :: physical_filename
    character (len=PGSd_PC_FILE_PATH_MAX) :: mcf_filename
    integer :: Version, Indx

    ! the group have to be defined as 49 characters long. The C interface is 50.
    ! The cfortran.h mallocs an extra 1 byte for the null character '\0/1, 
    ! therefore making the actual length of a 
    ! string pass of 50.

    character (len = PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)

    ! Externals

    integer, external :: PGS_MET_remove

    !Executable code

    module_error = 0
    if ( present(metadata_error)) metadata_error=1
    if(.NOT. TOOLKIT) return

    version = 1
    if ( MCF_FILE > 0 ) then
      returnStatus = PGS_PC_GetReference (MCF_FILE, version , mcf_filename)
    else
      returnStatus = PGSPC_W_NO_REFERENCE_FOUND
    end if

    if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
      & "Error: failed to find PCF ref for MCF_FILE in populate_metadata_oth.") 
      return
    end if

    version = 1
    returnStatus = PGS_PC_GetReference (HDF_FILE, version , physical_filename)

    if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
      & "Error: failed to find PCF ref for HDF_FILE in populate_metadata_oth.") 
      return
    end if

    call first_grouping(HDF_FILE, MCF_FILE, groups)

    do indx=1, numquantitiesperfile

      call measured_parameter (HDF_FILE, &
        & QuantityNames(indx), groups, indx)

    end do

    call third_grouping (HDF_FILE, hdf_sdid, groups, &
      & hdfVersion, filetype, l2metaData)

    if ( SFINBETWEENSTARTEND ) then
      sdid = mls_sfstart (physical_fileName, DFACC_RDWR, &
        & hdfVersion=hdfVersion, addingMetaData=.true.)
    else
      sdid = hdf_sdid
    endif 

    if ( sdid == -1 ) then
      call announce_error ( 0, &
      & "Error: failed to open the hdf file: "//TRIM(physical_fileName))
      return
    end if

    hdfReturn = mls_sfend( sdid, hdfVersion=hdfVersion, addingMetadata=.true. )
    if ( hdfReturn /= 0 ) then
        call announce_error ( 0, &
        & "Error: metadata mls_sfend in populate_metadata_oth.", &
        & error_number=hdfReturn) 
    end if

! Annotate the file with the PCF

    if ( ANNOTATEWITHPCF ) then
      call PCF2Hdr( physical_filename, hdfVersion, filetype=filetype )
      ! & hdfVersion, isHDFEOS)
    end if

    returnStatus = pgs_met_remove() 

    if ( present(metadata_error)) metadata_error=module_error

   if(BeVerbose( 'pro', -1 )) then
       call announce_success(physical_filename, mcf_filename, 'others')
   end if

  end subroutine Populate_metadata_oth

  ! ----------------------------------------readMCF -----
  subroutine readMCF ( MLSFILE, MCGROUPS )
    ! Given a metadata control file, read a group of inventory metadata control 
    ! parameters from it
    ! Note that the file must be an MLSFile_T
    type ( MLSFile_T ), intent(inout)            :: MLSFILE
    type ( MCGROUP_T ), pointer, dimension(:) :: MCGROUPS
    ! Internal variables
    integer :: list_len
    integer, parameter :: MAXLINELENGTH  = 256
    integer, parameter :: MAXARRAYLENGTH = 512
    character(len=MAXLINELENGTH), dimension(MAXARRAYLENGTH) :: list
    logical :: verbose

    ! Executable code
    verbose = BeVerbose( 'mcf', -1 )
    list = ' '
    call read_textfile( MLSFILE%Name, list, nLines=list_len )
    if ( verbose ) then
      call outputnamedValue( 'list_len', list_len )
      call dump( list(:list_len), 'list', width=1 )
    endif
    call mctree( list(:list_len), MCGroups )
  end subroutine readMCF

  ! ----------------------------------------writeMCF -----
  subroutine writeMCF ( MLSFILE, MCGROUPS )
    ! Given an MLSFile_T, write a group of inventory metadata control 
    ! parameters to it
    ! Note that the group may contain other groups, i.e. be nested
    type ( MLSFile_T ), intent(inout)            :: MLSFILE
    type ( MCGROUP_T ), intent(in), dimension(:) :: MCGROUPS
    ! Internal variables
    integer :: error
    integer :: i
    logical :: needMasterGroup

    ! Executable code
    if ( .not. MLSfile%stillOpen ) call mls_OpenFile ( MLSFile, error )
    call dump( MLSFile, details=1 )
    if ( error /= 0 ) call announce_error ( 0, &
      & 'error opening MLSFile in writeMCF' )
    needMasterGroup = .not. streq(MCGroups(1)%type, 'mastergroup', options='-cf' )
    call outputnamedvalue ( 'MCGroups(1)%type', trim(MCGroups(1)%type) )
    if ( needMasterGroup ) then
      write( MLSFile%fileID%f_id, * ) 'GROUP = INVENTORYMETADATA'
      write( MLSFile%fileID%f_id, * ) 'GROUPTYPE = MASTERGROUP'
      nest_degree = 1
    else
      nest_degree = 0
    endif
    do i=1, size(MCGroups)
      call writeMCGroup( MLSFile, MCGroups(i) )
    enddo
    if ( needMasterGroup ) write( MLSFile%fileID%f_id, * ) 'END_GROUP = INVENTORYMETADATA'
    write( MLSFile%fileID%f_id, * ) 'END'
    call mls_CloseFile ( MLSFile )
  end subroutine writeMCF

! Lori's routines

  ! -----------------------------------------------  WriteMetaLog  -----
  subroutine WriteMetaLog ( Metadata_error )

    ! Brief description of subroutine
    ! This subroutine writes metadata for the log file to a separate ASCII file.

    ! Arguments

    ! type( PCFData_T ), intent(in) :: Pcf
    integer, optional, intent(out) :: Metadata_error

    ! Parameters
    ! These are PCF numbers that *must* be in the PCF file

    integer, parameter :: ASCII_FILE = 101
    integer, parameter :: THE_LOG_FILE = 100

    ! Functions

    integer, external :: PGS_MET_init, PGS_MET_remove
    integer, external :: PGS_MET_setattr_s, PGS_MET_write
    integer, external :: PGS_PC_getconfigdata

    ! Internal

    character (len=1) :: NullStr
    character (len=45) :: Sval
    character (len=32) :: Mnemonic
    character (len=FileNameLen) :: Physical_filename
    character (len=PGSd_PC_FILE_PATH_MAX) :: mcf_filename
    character (len=480) :: Msg, Msr
    character (len=PGSd_MET_GROUP_NAME_L) :: Groups(PGSd_MET_NUM_OF_GROUPS)

    integer :: Result, Indx, Version
    logical, parameter :: DEBUG=.FALSE.

    ! Begin
    module_error = 0
    if ( present(metadata_error)) metadata_error=1
    if(.NOT. TOOLKIT) return

    nullStr = ''

    version = 1
    result = PGS_PC_GetReference (mlspcf_mcf_l2log_start, version, &
      & mcf_filename)

    if ( result /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
        & "Error: failed to find PCF ref for Log mcf in WriteMetaLog.") 
      return
    end if

    ! Check that the PCF file contains the entries for the ASCII file
    ! and for the log file itself
    result = pgs_pc_getconfigdata (ASCII_FILE, physical_filename)

    if ( result /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
        & "Error: failed to find PCF ref for the ASCII in WriteMetaLog.") 
      return
    end if

    version = 1
    result = PGS_PC_GetReference (THE_LOG_FILE, version , physical_filename)

    if ( result /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
        & "Error: failed to find PCF ref for the Log in WriteMetaLog.") 
      return
    end if

! Initialize the MCF file

    if ( DEBUG ) then
      call output('Initialize the MCF file', advance='yes')
    end if

    result = pgs_met_init(mlspcf_mcf_l2log_start, groups)
    if ( result /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
      & "metadata initialization error in WriteMetaLog.") 
      return
    end if

! Set PGE values

    result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), "LocalGranuleID", &
                               l2pcf%logGranID)

    sval = 'c' // l2pcf%cycle
    result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), "LocalVersionID", &
                               sval)

    indx = INDEX (l2pcf%startUTC, "T")

    result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                               "RangeBeginningDate", l2pcf%startUTC(1:indx-1))
!    sval= '00:00:00.000000'
    result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                               "RangeBeginningTime", l2pcf%startUTC(indx+1:))
    result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                               "RangeEndingDate", l2pcf%endUTC(1:indx-1))
!    sval= '23:59:59.999999'
    result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                               "RangeEndingTime", l2pcf%endUTC(indx+1:))

    result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), "PGEVersion", &
                               l2pcf%PGEVersion)

    if ( result /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
      & "Error in setting PGEVersion attribute in WriteMetaLog.", &
      &  error_number=result)
        return
    end if

    ! Write the metadata and their values to an ASCII file

    if ( DEBUG ) then
      call output('Write the metadata and their values to an ASCII file', &
      &  advance='yes')
    end if

    result = pgs_met_write(groups(1), nullStr, ASCII_FILE)

    if ( result /= PGS_S_SUCCESS ) then
       call Pgs_smf_getMsg(result, mnemonic, msg)
       msr = "Error: failed to write metadata in WriteMetaLog." &
       & // trim(mnemonic) // ':  ' // trim(msg)
       CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
    end if

    result = pgs_met_remove()

    if ( present(metadata_error)) metadata_error=module_error

   if(BeVerbose( 'pro', -1 )) then
       call announce_success(physical_filename, mcf_filename, 'Log')
   end if

   end subroutine WriteMetaLog

  ! ---------------------- Private procedures ------------

  ! ---------------------------------------------  announce_success  -----
  subroutine announce_success ( Name, mcf, l2_type )
    character(LEN=*), intent(in) :: Name
    character(LEN=*), intent(in) :: mcf
    character(LEN=*), intent(in) :: l2_type

    call output ( 'Level 2 output product metadata type : ' )
    call output ( trim(l2_type), advance='yes')
    call blanks(15)
    call output ( 'name : ' )
    call blanks(8)
    call output ( trim(Name), advance='yes')
    call blanks(15)
    call output ( 'mcf name : ' )
    call blanks(8)
    call output ( trim(mcf), advance='yes')

  end subroutine announce_success

  ! ------------------------------------------------  Announce_Error  -----
  subroutine Announce_Error ( lcf_where, full_message, use_toolkit, &
    & error_number )

    ! Arguments

    integer, intent(in) :: Lcf_where
    character(LEN=*), intent(in) :: Full_message
    logical, intent(in), optional :: Use_toolkit
    integer, intent(in), optional :: Error_number

    ! Local
    logical :: Just_print_it
    logical, parameter :: Default_output_by_toolkit = .true.
    character (LEN=132) :: attrname, errmsg
    integer :: Toolbox_error_num

    just_print_it = .not. default_output_by_toolkit
    if ( present(use_toolkit) ) just_print_it = .not. use_toolkit

    if ( .not. just_print_it ) then
      if ( TOOLKIT ) module_error = 1
      call output ( '***** At ' )

      if ( lcf_where > 0 ) then
        call print_source ( where(lcf_where) )
      else
        call output ( '(no lcf node available)' )
      end if

      call output ( " Caused the following error:", advance='yes', &
       & from_where=ModuleName)
      call output ( trim(full_message), advance='yes', &
        & from_where=ModuleName)
      if ( present(error_number) ) then
        call output ( 'Error number ', advance='no' )
        call output ( error_number, places=9, advance='yes' )
        Toolbox_error_num = error_number
        call Pgs_smf_getMsg (Toolbox_error_num, attrname, errmsg)
        call MLSMessage (MLSMSG_WARNING, ModuleName, &
            &  TRIM(attrname) //' ' // TRIM(errmsg) )
      end if
    else
      call output ( '***Error in module ' )
      call output ( ModuleName, advance='yes' )
      call output ( trim(full_message), advance='yes' )
      if ( present(error_number) ) then
        call output ( 'Error number ' )
        call output ( error_number, advance='yes' )
      end if
    end if

!===========================
  end subroutine Announce_Error
!===========================

  ! ---------------------------------------------  First_grouping  -----

  subroutine First_grouping ( HDF_FILE, MCF_FILE, Groups )

    ! This writes the metadata for the following attributes:
    ! (attributes marked automatic are not explicitly written, however)

    ! SizeMBECSDataGranule (automatic)
    ! ReProcessingPlanned
    ! LocalGranuleID
    ! DayNightFlag
    ! ProductionDateTime (automatic)
    ! LocalVersionID
    !

    !    USE InitPCFs, ONLY: L2PCF

    !Arguments

    integer :: HDF_FILE, MCF_FILE
    ! type(PCFData_T) :: l2pcf

    ! the group have to be defined as 49 characters long. The C interface is 50.
    ! The cfortran.h mallocs an extra 1 byte for the null character '\0/1, 
    ! therefore making the actual length of a 
    ! string pass of 50.

    character (len = PGSd_MET_GROUP_NAME_L) :: Groups(PGSd_MET_NUM_OF_GROUPS)

    !Local Variables

    integer, parameter :: INVENTORY=2 !, ARCHIVE=1
    character (len=132) :: Attrname
    integer :: Indx, Version
    character (len=PGSd_PC_FILE_PATH_MAX) :: Physical_filename
    character (len=PGSd_PC_FILE_PATH_MAX) :: Sval
    integer :: ReturnStatus

    ! Externals

    integer, external :: PGS_MET_init, &
      &  PGS_MET_setAttr_s
    ! Executable code

    version = 1
    returnStatus = PGS_PC_GetReference (HDF_FILE, version , physical_filename)

    if ( returnStatus /= PGS_S_SUCCESS ) then 
      call announce_error ( 0, &
      & "Error in getting ref for PCF number in 1st grouping." ) 
    end if

    returnStatus = pgs_met_init (MCF_FILE, groups)

    if ( returnStatus /= PGS_S_SUCCESS ) then 
      call announce_error ( 0, &
      & "Metadata initialization error" )
    end if

    ! Set PGE values 

    ! ECSDataGranule

    attrName = 'ReprocessingPlanned'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'further update anticipated using enhanced PGE')
    if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
        & "Error in writing ReprocessingPlanned attribute." )
    end if

    attrName = 'LocalGranuleID'
    sval = physical_filename
    indx = INDEX (sval, "/", .TRUE.) + 1  ! Begin after last "/"
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, sval(indx:))
    if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
       & "Error in writing LocalGranuleID attribute." ) 
    end if

    attrName = 'DayNightFlag'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, 'Both')
    if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
        & "Error in writing DayNightFlag attribute." )
    end if

    attrName = 'LocalVersionID'
    sval = 'c' // l2pcf%cycle
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, sval)
    if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
        & "Error in writing LocalVersionID attribute." ) 
    end if

  end subroutine First_grouping

  recursive subroutine insertGroup ( list, groups, i, name )
    ! Loop over list until end_group
    character(len=*), dimension(:), intent(in)   :: list
    type (MCGroup_T), dimension(:), pointer      :: groups
    integer, intent(inout)                       :: i
    character(len=*), intent(in)                 :: name
    ! Internal variables
    type (MCGroup_T)   :: group
    character(len=64)  :: key
    integer            :: newsize
    type (MCPARAM_T)   :: param
    character(len=128) :: value
    logical :: verbose
    ! Executable
    verbose = BeVerbose( 'mcf', -1 )
    group%name = name
    if ( verbose ) call output( 'Inserting group named '// trim(name), advance='yes' )
    do
      i = i+1
      if ( i > size(list) ) exit
      ! Which form is our list_line?
      if ( streq( list(i), 'end_group', options='-cfh' ) ) then
        newsize = addmcgrouptodatabase( groups, group )
        exit
      elseif ( any( streq( list(i), &
        & (/ 'end_parameter', 'end_object   ' /) , options='-cfh' ) ) )  then
        newsize = addmcparamtodatabase( group%params, param )
      elseif( isComment(list(i)) ) then
        ! Just a comment
      elseif( len_trim(list(i)) < 1 ) then
        ! Just a blank line
      elseif ( index( trim(list(i)), '=' ) > 0 ) then
        ! separate into key, value
        call GetStringElement ( list(i), key, &
          & 1, countEmpty, '=' )
        call GetStringElement ( list(i), value, &
          & 2, countEmpty, '=' )
        ! Are we beginning a new group or parameter? Maybe just grouptype?
        if ( streq( key, &
          & 'grouptype' , options='-cf' ) )  then
          group%type = value
          ! print *, "setting named group " // trim(name) // ' type to ' // trim(value)
        elseif ( streq( key, &
          & 'class' , options='-cf' ) )  then
          group%class = value
        elseif ( any( streq( key, &
          & (/ 'group_name', 'group     ' /) , options='-cfh' ) ) )  then
          call insertGroup ( list, group%groups, i, value )
        elseif ( any( streq( key, &
          & (/ 'parameter_name', 'object        ' /) , options='-cfh' ) ) )  then
          call insertParameter ( list, group%params, i, value )
        else
          call setParam ( param, key, value )
        endif
      else
        ! unrecognized line
        call output( trim(list(i)) // ' is unrecognized', advance='yes' )
      endif
    enddo
  end subroutine insertGroup

  recursive subroutine insertParameter ( list, parameters, i, name )
    ! Loop over list until end_parameter
    character(len=*), dimension(:), intent(in)   :: list
    type (MCParam_T), dimension(:), pointer      :: parameters
    integer, intent(inout)                       :: i
    character(len=*), intent(in)                 :: name
    ! Internal variables
    integer :: newsize
    character(len=64)  :: key
    type (MCParam_T)   :: param
    character(len=128) :: value
    logical :: verbose
    ! Executable
    verbose = BeVerbose( 'mcf', -1 )
    if ( verbose ) call output( 'Inserting parameter named '// trim(name), advance='yes' )
    param%name = name
    do
      i = i+1
      if ( i > size(list) ) exit
      ! Which form is our list_line?
      ! print *, trim(list(i)), streq( list(i), 'end_object', options='-cfh' )
      if ( any( streq( list(i), &
        & (/ 'end_parameter', 'end_object   ' /) , options='-cfh' ) ) )  then
        newsize = addmcparamtodatabase( parameters, param )
        exit
      elseif( isComment(list(i)) ) then
        ! Just a comment
      elseif( len_trim(list(i)) < 1 ) then
        ! Just a blank line
      elseif ( index( trim(list(i)), '=' ) > 0 ) then
        ! separate into key, value
        call GetStringElement ( list(i), key, &
          & 1, countEmpty, '=' )
        call GetStringElement ( list(i), value, &
          & 2, countEmpty, '=' )
        if ( any( streq( key, &
          & (/ 'group_name', 'group     ' /) , options='-cfh' ) ) )  then
          call insertGroup ( list, param%groups, i, value )
        elseif ( any( streq( key, &
          & (/ 'parameter_name', 'object        ' /) , options='-cfh' ) ) )  then
          call insertParameter ( list, param%params, i, value )
        else
          call setParam ( param, key, value )
        endif
      else
        ! unrecognized line
        call output( trim(list(i)) // ' is unrecognized', advance='yes' )
      endif
    enddo
  end subroutine insertParameter

  subroutine setParam ( param, key, value )
    ! set key component of param to value
    type (MCParam_T), intent(inout):: param
    character(len=*), intent(in)   :: key, value
    ! Internal variables
    ! Executable
    ! print *, 'Setting parameter key named ', trim(key), ' to ', trim(value)
    select case (adjustl(lowercase(key)))
    case ('name')
      param%name          = value
    case ('mandatory')
      param%mandatory     = value
    case ('data_location')
      param%data_location = value
    case ('class')
      param%class         = value
    case ('num_val')
      param%num_val       = value
    case ('type')
      param%type          = value
    case ('value')
      param%value         = value
    ! case default
    end select
  end subroutine setParam

  ! -----------------------------------------  Measured_parameter  -----

  subroutine Measured_parameter ( HDF_FILE, Field_name, Groups, Class_num)

    ! This writes the attributes corresponding to the measured parameter container:
    !
    ! ParameterName
    ! AutomaticQualityFlag
    ! AutomaticQualityFlagExplanation
    ! OperationalQualityFlag
    ! OperationalQualityFlagExplanation
    ! ScienceQualityFlag
    ! ScienceQualityFlagExplanation
    ! QAPercentInterpolatedData
    ! QAPercentMissingData
    ! QAPercentOutOfBoundsData
    !

    !    USE InitPCFs, ONLY: L2PCF

    !Arguments

    integer                       :: HDF_FILE
    character(len=*)              :: Field_name
    integer, intent(in)           :: Class_num

    ! the group have to be defined as 49 characters long. The C interface is 50.
    ! The cfortran.h mallocs an extra 1 byte for the null character '\0/1, 
    ! therefore making the actual length of a 
    ! string pass of 50.

    character (len = PGSd_MET_GROUP_NAME_L) :: Groups(PGSd_MET_NUM_OF_GROUPS)

    ! Local Variables

    integer, parameter :: INVENTORY=2 !, ARCHIVE=1
    character (len=132) :: Attrname
    character (len=2) :: Class
    character (len=PGSd_PC_FILE_PATH_MAX) :: Physical_filename
    integer :: ReturnStatus
    character (len=PGSd_PC_FILE_PATH_MAX) :: Sval
    integer :: Version

    ! Externals

    integer, external :: &
      &  PGS_MET_setAttr_s, PGS_MET_SETATTR_I

    !Executable code
    version = 1
    returnStatus = PGS_PC_GetReference (HDF_FILE, version , physical_filename)


    if ( returnStatus /= PGS_S_SUCCESS ) then 
      call announce_error ( 0, &
        & "Error in getting ref for PCF id in measured_parameter." )
    end if

    ! A moment of misplaced priorities: we have decided to create
    ! links between valid data field names and the name L2gpValue
    ! Instead of doing this in a routine created for that purpose,
    ! we have violated good design principles and stuck it here for now
    !                (moved to WriteL2GPData)
    ! MeasuredParameterContainer
    if ( class_num <= 0 ) then
      class(:2) = '0 '
    else if ( class_num < 10 ) then
      write ( class(:1), '(I1)' ) class_num
      class(2:2) = ' '
    else
      write ( class(:2), '(I2)' ) class_num
    end if

    if ( field_name /= ' ' ) then
      sval = adjustl(field_name)
    else
      sval = 'Miscellaneous'
    end if
    attrName = 'ParameterName' // '.' // class
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, sval)
    if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
        & "Error in writing ParameterName attribute." )
    end if

    ! QAFlags Group
    ! These have been moved to the MCF files
    ! QAStats Group

    attrName = 'QAPercentInterpolatedData' // '.' // class
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, 0)
    if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
        & "Error in writing QAPercentInterpolatedData attribute." )
    end if

    attrName = 'QAPercentMissingData' // '.' // class
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, 0)
    if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
        & "Error in writing QAPercentMissingData attribute." )
    end if

    attrName = 'QAPercentOutofBoundsData' // '.' // class
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, 0)
    if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
        & "Error in writing QAPercentOutofBoundsData attribute." )
    end if

  end subroutine Measured_parameter

  ! ---------------------------------------- PCF2Hdr -----
  ! Make sure file type will be recognized
  ! before calling writePCF2Hdr
  ! (which only knows l_hdf, l_hdfeos, l_swath)
  subroutine PCF2Hdr( filename, hdfVersion, filetype )
    ! Args
    character(len=*), intent(in)  :: filename
    integer, optional, intent(in) :: hdfVersion
    integer, optional, intent(in) :: fileType
    ! Internal variables
    integer :: myType
    ! Executable
    if ( .not. associated(l2pcf%anText) ) return
    myType = l_hdf
    if ( present(filetype) ) then
      myType = filetype
      select case(filetype)
      case (l_l2dgg,l_l2gp)
        myType = l_swath
      ! case default
        ! myType = filetype
      end select
    endif
    call writePCF2Hdr( filename, l2pcf%anText, hdfVersion, myType )

  end subroutine PCF2Hdr

  ! ---------------------------------------------  Third_grouping  -----

  subroutine Third_grouping ( HDF_FILE, hdf_sdid, Groups, &
    & hdfVersion, filetype, l2metaData )

    ! This writes the following metadata attributes:

    ! OrbitNumber
    ! StartOrbitNumber
    ! StopOrbitNumber
    ! EquatorCrossingLongitude
    ! EquatorCrossingTime
    ! EquatorCrossingDate
    ! ShortName
    ! VersionID
    ! InputPointer
    ! LocalityValue
    ! VerticalSpatialDomainType
    ! VerticalSpatialDomainValue
    ! ZOneIdentifier
    ! WestBoundingCoordinate
    ! NortBoundingCoordinate
    ! EastBoundingCoordinate
    ! SouthBoundingCoordinate
    ! RangeBeginningDate
    ! RangeBeginningTime
    ! RangeEndingDate
    ! RangeEndingTime
    ! PGEVersion
    !


    ! Arguments

    integer :: HDF_FILE
    integer :: hdf_sdid

    ! the group have to be defined as 49 characters long. The C interface is 50.
    ! The cfortran.h mallocs an extra 1 byte for the null character '\0/1, 
    ! therefore making the actual length of a 
    ! string pass of 50.

    character (len = PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)
    integer, optional, intent(in) :: hdfVersion
    integer, optional, intent(in)  :: filetype  ! 'sw' or 'hdf'
    type (L2Metadata_T) :: l2metaData

    !Local Variables

    integer :: HdfReturn
    integer :: ReturnStatus

    real(r8) :: Dval
    integer, parameter :: INVENTORY=2 !, ARCHIVE=1
    character (len=PGSd_PC_FILE_PATH_MAX) :: physical_filename
    character (len=PGSd_PC_FILE_PATH_MAX) :: sval
    character (len=132) :: attrname
    integer :: version, indx
    integer :: inptptr_filetype
    integer :: startOrbit, stopOrbit


    ! Externals

    integer, external :: PGS_MET_setattr_d, &
         PGS_MET_setAttr_s, PGS_MET_SETATTR_I, &
         PGS_MET_write, PGS_MET_remove
    !Executable code

    inptptr_filetype = l_swath
    if ( present(filetype) ) inptptr_filetype = filetype
    version = 1
    returnStatus = PGS_PC_GetReference (HDF_FILE, version , physical_filename)

    if ( returnStatus /= PGS_S_SUCCESS ) then 
      call announce_error ( 0, &
      & "Error in getting ref for PCF id in third_grouping.") 
    end if

    ! Orbit Calculated Spatial Domain Container
    ! Per James Johnson's email 6/12/03, use either OrbitNumber OR 
    ! (StartOrbitNumber and StopOrbitNumber)

    !attrName = 'OrbitNumber' // '.1'
    !returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, -1)
    ! This change to confirm James Johnson suggestion on 6/12/03
    !returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, 999)
    !if ( returnStatus /= PGS_S_SUCCESS ) then
    !  call announce_error ( 0, &
    !  & "Error in writing OrbitNumber attribute.") 
    !end if

    ! Start, Stop orbit numbers: level one has actual calculated numbers
    ! but, for now at least, we'll not trouble
    ! For now, use 99999 for invalid value
    ! Write start, stop orbit numbers from l1 to l2

    attrName = 'StartOrbitNumber' // '.1'
    if (GlobalAttributes%OrbNum(1) == -1) then
        startOrbit = 99999
    else
        startOrbit = GlobalAttributes%OrbNum(1)
    end if
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, &
        & startOrbit)
    !returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, 99999)
    !returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, -1)
    if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
      & "Error in writing StartOrbitNumber attribute.")
    end if

    attrName = 'StopOrbitNumber' // '.1'
    if (maxval(GlobalAttributes%OrbNum) == -1) then
        stopOrbit = 99999
    else
        stopOrbit = maxval(GlobalAttributes%OrbNum)
    end if
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, &
        & stopOrbit)
    !returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, 99999)
    !returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, -1)
    if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
      & "Error in writing StopOrbitNumber attribute.")
    end if

    attrName = 'EquatorCrossingLongitude' // '.1'
    dval = 0.0
    returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
    if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
      & "Error in writing EquatorCrossingLongitude attribute.") 
    end if

    attrName = 'EquatorCrossingTime' // '.1'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         '00:00:00')
    if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
      & "Error in writing EquatorCrossingTime attribute.") 
    end if

    indx = INDEX (L2PCF%startUTC, "T")
    attrName = 'EquatorCrossingDate' // '.1'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         L2PCF%startUTC(1:indx-1))
    if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
      & "Error in writing EquatorCrossingDate attribute.") 
    end if

    ! InputPointer

    attrName = 'InputPointer'
    if ( SETINPUTPOINTER ) then
      call announce_error ( 0, &
      & "Invalid method of writing InputPointer attribute.") 
    else
      returnStatus = WriteInputPointer(groups(INVENTORY), attrName, &
        & filetype=inptptr_filetype)
    endif
    if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
      & "Error in writing InputPointer attribute.") 
    end if

    ! Locality Value

    attrName = 'LocalityValue'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'Limb')
    if ( returnStatus /= PGS_S_SUCCESS ) then
       call announce_error ( 0, &
       & "Error in writing LocalityValue attribute.")
     end if

     ! VerticalSpatialDomain Product-Specific Attribute

     attrName = 'VerticalSpatialDomainType' // '.1'
     returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
          'Atmosphere Layer')
     if ( returnStatus /= PGS_S_SUCCESS ) then
       call announce_error ( 0, &
       & "Error in writing VerticalSpatialDomainType attribute.")
     end if

     attrName = 'VerticalSpatialDomainValue' // '.1'
     !returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
     !     'Brightness Temperature')
     ! The following change bring level 2 into agreement level 3
     ! Further change to confirm James Johnson suggestion 6/12/03
     ! will await in the next delivery which it also unify all 3 levels
     returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
          'Atmosphere Profile')
     if ( returnStatus /= PGS_S_SUCCESS ) then
       call announce_error ( 0, &
       & "Error in writing VerticalSpatialDomainValue attribute.")
     end if

     ! HorizontalSpatialDomainContainer

     attrName = 'ZoneIdentifier'
     returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
          'Other Grid System')
     if ( returnStatus /= PGS_S_SUCCESS ) then
       call announce_error ( 0, &
       & "Error in writing ZoneIdentifier attribute.")
     end if

     if ( .not. sharedPCF ) then
        l2metaData%minLon = -180.0
        l2metaData%maxLon = 180.0
        l2metaData%minLat = -90.0
        l2metaData%maxLat = 90.0
     end if

     attrName = 'WestBoundingCoordinate'
     dval = l2metaData%minLon 
     returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
     if ( returnStatus /= PGS_S_SUCCESS ) then
       call announce_error ( 0, &
       & "Error in writing WestBoundingCoordinate attribute.")
     end if

     attrName = 'NorthBoundingCoordinate'
     dval = l2metaData%maxLat 
     returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
     if ( returnStatus /= PGS_S_SUCCESS ) then
       call announce_error ( 0, &
       & "Error in writing NorthBoundingCoordinate attribute.")
     end if

     attrName = 'EastBoundingCoordinate'
     dval = l2metaData%maxLon 
     returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
     if ( returnStatus /= PGS_S_SUCCESS ) then
       call announce_error ( 0, &
       & "Error in writing EastBoundingCoordinate attribute.")
     end if

     attrName = 'SouthBoundingCoordinate'
     dval = l2metaData%minLat 
     returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
     if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
      & "Error in writing SouthBoundingCoordinate attribute.")
    end if

    indx = INDEX (L2PCF%startUTC, "T")

    ! RangeDateTime Group

    attrName = 'RangeBeginningDate'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), &
         attrName, L2PCF%startUTC(1:indx-1))
    if ( returnStatus /= PGS_S_SUCCESS ) then 
      call announce_error ( 0, &
      & "Error in setting RangeBeginningDate attribute.") 
    end if

    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), &
         "RangeBeginningTime", L2PCF%startUTC(indx+1:))
    if ( returnStatus /= PGS_S_SUCCESS ) then 
      call announce_error ( 0, &
      & "Error in setting RangeBeginningTime attribute.") 
    end if

    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), "RangeEndingDate", &
         L2PCF%endUTC(1:indx-1))
    if ( returnStatus /= PGS_S_SUCCESS ) then 
      call announce_error ( 0, &
      & "Error in setting RangeEndingDate attribute.") 
    end if

    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), "RangeEndingTime", &
         L2PCF%endUTC(indx+1:))
    if ( returnStatus /= PGS_S_SUCCESS ) then 
      call announce_error ( 0, &
      & "Error in setting RangeEndingTime attribute.") 
    end if

    ! PGEVersion

    attrName = 'PGEVersion'
    sval = L2PCF%PGEVersion
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, sval)
    if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
      & "Error in setting PGEVersion attribute.") 
    end if

    ! Production Location

    attrName = 'ProductionLocation'
    sval = GlobalAttributes%productionLoc
    returnStatus = pgs_met_setAttr_s ( groups(INVENTORY), attrName, sval )
    if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
      & "Error in setting ProductionLocation attribute.") 
    end if

    attrName = 'DataProducer'
    sval = 'MLS_SIPS'
    returnStatus = pgs_met_setAttr_s ( groups(INVENTORY), attrName, sval )
    if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
      & "Error in setting DataProducer attribute.") 
    end if

    ! DOI
    ! This should be in the MCF
    attrName = 'identifier_product_doi'
    sval = l2metaData%doiIdentifier
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, sval)
    if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( 0, &
      & "Error in setting DOI attribute.") 
    end if

    ! All done
    hdf_sdid = mls_sfstart (physical_fileName, DFACC_RDWR, &
      & hdfVersion=hdfVersion, addingMetaData=.true.) 

    if ( hdf_sdid == -1 ) then
      call announce_error ( 0, &
      & "Error: failed to open hdf file: "//physical_fileName) 
    end if

    returnStatus = pgs_met_write (groups(INVENTORY), "coremetadata.0", hdf_sdid)

    if ( returnStatus /= PGS_S_SUCCESS .AND. &
      &  returnStatus /= PGSMET_W_METADATA_NOT_SET ) then 
      if ( returnStatus == PGSMET_W_METADATA_NOT_SET ) then 
        call announce_error ( 0, &
        & "Error: some of the mandatory parameters not set.") 
      else
        call output('hdf_sdid: ', advance='no')
        call output(hdf_sdid, advance='yes')
        call output('hdfVersion: ', advance='no')
        call output(hdfVersion, advance='yes')
        call announce_error ( 0, &
        & "Error: metadata write failed in third_grouping." // trim(physical_fileName), &
        & error_number=returnStatus) 
      end if
    end if

    if ( SFINBETWEENSTARTEND ) then
      hdfReturn = mls_sfend( hdf_sdid, hdfVersion=hdfVersion, addingMetadata=.true. )
      if ( hdfReturn /= 0 ) then
          call announce_error ( 0, &
          & "Error: metadata mls_sfend in third_grouping.", &
          & error_number=hdfReturn) 
      end if

      returnStatus = pgs_met_remove() 
    endif

  end subroutine Third_grouping

  ! ----------------------------------------writeMCGroup -----
  recursive subroutine writeMCGroup ( MLSFILE, MCGROUP )
    ! Given an MLSFile_T, write a set of inventory metadata control parameters 
    ! to it
    ! If the group contains other groups, do what seems proper
    type ( MLSFile_T ), intent(inout)            :: MLSFILE
    type ( MCGROUP_T ), intent(in)               :: MCGROUP
    ! Internal variables
    integer :: i
    integer :: indent

    ! Executable code
    call output( 'Writing group ' // trim(MCGroup%name), advance='yes' )
    if ( len_trim(MCGroup%name) < 1 ) return
    nest_degree = nest_degree + 1
    indent = nest_degree*3
    ! write( MLSFile%fileID%f_id, * ) spaces(:indent), '#GROUP_namelen = ', len_trim(MCGroup%name)
    write( MLSFile%fileID%f_id, * ) spaces(:indent), 'GROUP = ' // trim(adjustl(MCGroup%name))
    indent = indent + 3
    if ( len_trim(MCGroup%class) > 0 ) &
      & write( MLSFile%fileID%f_id, * ) spaces(:indent), &
      & 'CLASS = ' // trim(adjustl(MCGroup%class))
    if ( len_trim(MCGroup%type) > 0 ) &
      & write( MLSFile%fileID%f_id, * ) spaces(:indent), &
      & 'GROUPTYPE = ' // trim(adjustl(MCGroup%type))
    if ( associated(MCGroup%groups) ) then
      do i=1, size(MCGroup%groups)
        call writeMCGroup( MLSFile, MCGroup%groups(i) )
      enddo
    endif
    if ( associated(MCGroup%params) ) then
      do i=1, size(MCGroup%params)
        call writeMCParam( MLSFile, MCGroup%params(i) )
      enddo
    endif
    indent = nest_degree*3
    write( MLSFile%fileID%f_id, * )  spaces(:indent), 'END_GROUP = ' // trim(adjustl(MCGroup%name))
    nest_degree = nest_degree - 1
  end subroutine writeMCGroup

  recursive subroutine writeMCParam ( MLSFILE, MCPARAM )
    type ( MLSFile_T ), intent(inout)            :: MLSFILE
    type ( MCParam_T ), intent(in)               :: MCPARAM
    ! Internal variables
    integer :: i
    integer :: indent
    ! Executable
    call output( 'Writing parameter ' // trim(MCParam%name), advance='yes' )
    if ( len_trim(MCParam%name) < 1 ) return
    nest_degree = nest_degree + 1
    indent = nest_degree*3
    write( MLSFile%fileID%f_id, * ) spaces(:indent), &
      & 'OBJECT = ' // trim(adjustl(MCParam%name))
    indent = indent + 3
    write( MLSFile%fileID%f_id, * ) spaces(:indent), &
      & 'Mandatory = ' // trim(adjustl(MCParam%mandatory))
    write( MLSFile%fileID%f_id, * ) spaces(:indent), &
      & 'Data_Location = ' // trim(adjustl(MCParam%data_location))
    if ( len_trim(MCParam%class) > 0 ) &
      & write( MLSFile%fileID%f_id, * ) spaces(:indent), &
      & 'CLASS = ' // trim(adjustl(MCParam%class))
    if ( len_trim(MCParam%num_val) > 0 ) &
      & write( MLSFile%fileID%f_id, * ) spaces(:indent), &
      & 'NUM_VAL = ' // trim(adjustl(MCParam%num_val))
    if ( len_trim(MCParam%type) > 0 ) &
      & write( MLSFile%fileID%f_id, * ) spaces(:indent), &
      & 'TYPE = ' // trim(adjustl(MCParam%type))
    if ( len_trim(MCParam%value) > 0 ) &
      & write( MLSFile%fileID%f_id, * ) spaces(:indent), &
      & 'Value = ' // trim(adjustl(MCParam%value))

    if ( associated(MCParam%groups) ) then
      do i=1, size(MCParam%groups)
        call writeMCGroup( MLSFile, MCParam%groups(i) )
      enddo
    endif
    if ( associated(MCParam%params) ) then
      do i=1, size(MCParam%params)
        call writeMCParam( MLSFile, MCParam%params(i) )
      enddo
    endif
    indent = nest_degree*3
    write( MLSFile%fileID%f_id, * ) spaces(:indent), &
      & 'END_OBJECT = ' // trim(adjustl(MCParam%name))
    nest_degree = nest_degree - 1
  end subroutine writeMCParam

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module WriteMetadata 
! $Log$
! Revision 2.84  2018/05/31 22:48:40  pwagner
! Changed name of attribute to identifier_product_doi to please DAAC
!
! Revision 2.83  2018/04/19 01:14:16  vsnyder
! Remove USE statements for unused names
!
! Revision 2.82  2015/03/28 02:56:12  vsnyder
! Added stuff to trace allocate/deallocate addresses
!
! Revision 2.81  2014/09/05 01:28:53  vsnyder
! More complete and accurate allocate/deallocate size tracking
!
! Revision 2.80  2014/04/07 18:10:09  pwagner
! Added new metadata attribute: DataProducer
!
! Revision 2.79  2014/03/26 17:47:44  pwagner
! Added ProductionLocation, identifier_product_doi metadata
!
! Revision 2.78  2014/03/07 19:30:30  pwagner
! Housekeeping; insert comments suggesting use of NameLen instead of hard-coded charlens
!
! Revision 2.77  2014/01/09 00:30:24  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.76  2013/12/05 01:42:03  pwagner
! Added RunID component to pcf; started using BeVerbose
!
! Revision 2.75  2013/09/24 23:47:23  vsnyder
! Use Where instead of Source_Ref for messages
!
! Revision 2.74  2013/07/26 00:07:35  pwagner
! Added routines to read, write, and parse MCFs
!
! Revision 2.73  2013/07/18 22:39:37  pwagner
! Added datatypes, functions for metadata control
!
! Revision 2.72  2013/04/19 20:07:28  pwagner
! Fixed long-standing bug in closing file after adding metadata
!
! Revision 2.71  2013/04/12 00:06:35  pwagner
! Added addToMetadata
!
! Revision 2.70  2011/07/12 22:35:03  honghanh
! Change l_grid to l_hdfeos
!
! Revision 2.69  2011/06/29 21:52:53  pwagner
! no metadata is always an error now
!
! Revision 2.68  2011/05/09 18:28:30  pwagner
! Converted to using switchDetail
!
! Revision 2.67  2009/06/23 18:46:19  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.66  2009/06/02 17:53:15  cvuu
! Add NRT Lat and Lon bounding to metadata
!
! Revision 2.65  2008/05/02 00:33:53  vsnyder
! Delete unused symbol
!
! Revision 2.64  2008/02/22 21:37:05  pwagner
! DGG file was hybrid by default; now it will be pure HDFEOS
!
! Revision 2.63  2007/12/07 01:50:32  pwagner
! Removed unused dummy variables, etc.
!
! Revision 2.62  2007/06/21 00:54:08  vsnyder
! Remove tabs, which are not part of the Fortran standard
!
! Revision 2.61  2006/03/15 23:52:24  pwagner
! Removed InputVersion component from PCF, l2cf
!
! Revision 2.60  2006/02/21 19:12:00  pwagner
! GetHashElement is now a generic
!
! Revision 2.59  2005/07/21 23:45:03  pwagner
! Removed unused l1b fileinfo fields from l2pcf
!
! Revision 2.58  2005/06/22 18:57:02  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.57  2005/01/27 00:35:06  pwagner
! ReprocessingActual field dropped from product metadata
!
! Revision 2.56  2005/01/14 21:37:30  pwagner
! ch3cn now has mcf, so special treatment lifted
!
! Revision 2.55  2004/12/15 23:34:32  pwagner
! Change REPROCESSINGACTUAL to unknown; light housecleaning
!
! Revision 2.54  2004/12/14 21:41:53  pwagner
! Repaired double-0 cycle; automaticQ.. and OperationalQ.. shifted to mcf
!
! Revision 2.53  2004/08/04 23:19:58  pwagner
! Much moved from MLSStrings to MLSStringLists
!
! Revision 2.52  2004/01/30 00:31:10  pwagner
! Stops useless warnings about pgs_met_remove return value
!
! Revision 2.51  2003/09/12 16:32:00  cvuu
! Output the right orbit numbers for StartOrbitNumber and StopOrbitNumber
!
! Revision 2.50  2003/09/04 22:41:09  pwagner
! Added some extra printing if metadata write fails in third_grouping
!
! Revision 2.49  2003/09/03 23:56:24  pwagner
! Can get species name from trivial file name fragment
!
! Revision 2.48  2003/08/11 17:29:36  cvuu
! Change to output attributes StartOrbitNumber and StopOrbitNumber instead of OrbitNumber
!
! Revision 2.47  2003/08/01 20:25:08  pwagner
! Removed ref to undefined sd_path
!
! Revision 2.46  2003/07/23 18:27:44  cvuu
! brought closer to James Johnson want to; quick and dirty fixed for CH3CN
!
! Revision 2.45  2003/07/07 23:48:18  pwagner
! Changed in interfaces to make filetype a lit_name; l2pcf a saved variable
!
! Revision 2.44  2003/06/09 22:49:35  pwagner
! Reduced everything (PCF, PUNISH.., etc.) to TOOLKIT
!
! Revision 2.43  2003/05/30 23:50:32  pwagner
! Relies on lib/PCFHdr to know what mesg goes in InputPointer
!
! Revision 2.42  2003/04/03 22:58:40  pwagner
! Alias now set in lib/L2GPData instead of l2/write_meta
!
! Revision 2.41  2003/03/15 00:15:52  pwagner
! Wont quit if pgs_met_remove returns non-zero value
!
! Revision 2.40  2003/03/11 00:21:36  pwagner
! Interfaces fit new WritePCF2Hdr flixibility
!
! Revision 2.39  2003/03/07 00:55:38  pwagner
! Set INPUTPOINTERMESG when annotating with PCF as a private param /HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/PCF
!
! Revision 2.38  2003/02/27 21:54:28  pwagner
! Now annotates with PCF; does not set inputPointer
!
! Revision 2.37  2003/02/10 22:02:41  pwagner
! Passes isHDFEOS to PCFHdr
!
! Revision 2.36  2003/02/01 00:29:40  pwagner
! Passes hdfVersion to writePCF2Hdr
!
! Revision 2.35  2003/01/30 01:01:05  pwagner
! Lets PCFHdr prepare and write input pointer
!
! Revision 2.34  2003/01/14 00:42:24  pwagner
! Creates soft link to type2precisionname, too
!
! Revision 2.33  2002/12/11 22:21:05  pwagner
! Makes soft link to data field name from L2gpValue field in hdf5 l2gp
!
! Revision 2.32  2002/12/06 23:37:22  pwagner
! addingMetaData now optional arg to mls_sfstart
!
! Revision 2.31  2002/11/22 12:31:16  mjf
! Added nullify routine(s) to get round Sun's WS6 compiler not
! initialising derived type function results.
!
! Revision 2.30  2002/10/08 17:36:23  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.29  2002/08/29 16:55:35  pwagner
! Moved writeInputPointer to PCFHdr in lib
!
! Revision 2.28  2002/08/28 22:27:26  pwagner
! Now writes input pointer metadata like L3MMData
!
! Revision 2.27  2002/03/13 18:31:01  pwagner
! Opens, writes metadata, closes each file only once
!
! Revision 2.26  2002/02/22 19:20:12  pwagner
! Added mcf file name to announce_success
!
! Revision 2.25  2002/02/22 01:18:06  pwagner
! get_l2gp_mcf accepts meta_name arg; now uses getPCFFromRef
!
! Revision 2.24  2002/01/31 00:37:45  pwagner
! Checks return value of mls_sfend
!
! Revision 2.23  2002/01/26 00:08:20  pwagner
! Housekeeping; changed proclaim to announce_success
!
! Revision 2.22  2002/01/23 21:50:04  pwagner
! Uses hdfVersion optional parameter
!
! Revision 2.21  2002/01/22 17:45:26  pwagner
! Removed bogus declaration of pgs_met_startf
!
! Revision 2.20  2002/01/18 23:11:09  pwagner
! Uses mls_sfstart; writePCF2Hdr controlled by ANNOTATEWITHPCF; dumps core (hooray)
!
! Revision 2.19  2002/01/09 23:54:24  pwagner
! Now gets FileNameLen from MLSCommon, not SDPToolkit
!
! Revision 2.18  2001/05/17 22:33:27  pwagner
! Prints info if pro switch set
!
! Revision 2.17  2001/05/08 23:28:46  pwagner
! Sorry-bad metadata functions in SDPToolkit
!
! Revision 2.16  2001/05/07 23:30:05  pwagner
! Gets all pgs_ from SDPToolkit
!
! Revision 2.15  2001/05/07 18:04:51  pwagner
! Always check for CREATEMETADATA first
!
! Revision 2.14  2001/05/03 20:34:20  vsnyder
! Cosmetic changes
!
! Revision 2.13  2001/04/20 23:52:15  vsnyder
! Consider the penalty from MLSL2Options in Announce_Error.  Numerous
! cosmetic changes.  Remove dump of tree node from Announce_Error.
!
! Revision 2.12  2001/04/19 23:51:40  pwagner
! Moved anText to become component of PCFData_T
!
! Revision 2.11  2001/04/16 23:49:11  pwagner
! Tiny change to announce_error
!
! Revision 2.10  2001/04/16 17:43:35  pwagner
! mcf text matching case sensitive if MCFCASESENSITIVE
!
! Revision 2.9  2001/04/13 23:47:26  pwagner
! Writes multiple measuredcontainers properly
!
! Revision 2.8  2001/04/13 00:28:03  pwagner
! Turned get_l2gp_mcf into a subroutine
!
! Revision 2.7  2001/04/12 00:22:17  pwagner
! Added announce_error
!
! Revision 2.6  2001/04/11 20:20:37  pwagner
! Checks correctly on PCF file before attempting LOG.met
!
! Revision 2.5  2001/04/10 23:03:57  pwagner
! Finally seems to work
!
! Revision 2.4  2001/04/09 23:45:47  pwagner
! Deleted unused old populate_metadata; some fixes
!
! Revision 2.3  2001/04/03 23:51:28  pwagner
! Many changes; some may be right
!
! Revision 2.2  2001/04/02 23:42:18  pwagner
! Added populate_metadata_oth
!
! Revision 2.1  2001/02/13 22:59:36  pwagner
! l2 modules can only use MLSPCF2
!
! Revision 2.0  2000/09/05 18:57:07  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.4  2000/06/30 00:16:41  lungu
! Made dval REAL(r8).
!
