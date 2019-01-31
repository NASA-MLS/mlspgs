! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module L2AUXData                 ! Data types for storing L2AUX data internally

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test, &
    & Test_Allocate, Test_Deallocate
  use Dump_0, only: Dump
  use HDF, only: Dfacc_Create, Dfacc_Rdonly, Dfacc_Rdwr, Dfnt_Float32, &
    & Dfnt_Int32, Sfcreate, Sfdimid, Sfendacc, Sfgdinfo, Sfginfo, &
    & Sfn2index, SfrData_F90, Sfsdmname, Sfsdscale, Sfselect, SfwData_F90
  use Init_Tables_Module, only: L_Baseline, L_Channel, L_Chisqchan, &
    & L_ChisqmMAF, L_ChisqmMIF, L_Chunk, L_Cloudextinction, &
    & L_CloudinducedRadiance, L_Cloudradsensitivity, L_Cloudwater, &
    & L_Dnwt_Ajn, L_Dnwt_Axmax, L_Dnwt_Cait, L_Dnwt_Chisqminnorm, &
    & L_Dnwt_Chisqnorm, L_Dnwt_Diag, L_Dnwt_Dxdx, L_Dnwt_Dxdxl, &
    & L_Dnwt_Dxn, L_Dnwt_Dxnl, L_Dnwt_Flag, L_Dnwt_Fnmin, L_Dnwt_Fnorm, &
    & L_Dnwt_Gdx, L_Dnwt_Gfac, L_Dnwt_Gradn, L_Dnwt_Sq, L_Dnwt_Sqt, &
    & L_Effectiveopticaldepth, L_Elevoffset, L_Frequency, L_Geodangle, &
    & L_Height, L_Heightoffset, L_Intermediatefrequency, L_Iteration, &
    & L_Jacobian_Cols, L_Jacobian_Rows, L_Limbsidebandfraction, &
    & L_Lostransfunc, L_Losvel, L_Lsbfrequency, L_Maf, &
    & L_Massmeandiameterice, L_Massmeandiameterwater, L_Mif, &
    & L_Mifextinction, L_Mifextinctionv2, L_Noisebandwidth, L_None, &
    & L_NoradsperMIF, L_Numj, L_Opticaldepth, L_Orbitinclination, L_Ascdescmode, &
    & L_Phasetiming, L_Phitan, L_Pressure, L_Ptan, L_Radiance, &
    & L_Reflspill, L_Refltemp, L_Scanresidual, L_Sceci, L_Scgeocalt, &
    & L_Scveleci, L_Scvelecr, L_SingleChannelRadiance, L_Sizedistribution, &
    & L_SpaceRadiance, L_StrayRadiance, L_Surfacetype, &
    & L_Systemtemperature, L_Tngteci, L_Tngtgeocalt, L_Tngtgeodalt, &
    & L_Totalextinction, L_Usbfrequency, L_Vmr, L_Xyz
  use Intrinsic, only: L_HDF, Lit_Indices
  use Lexer_Core, only: Print_Source
  use MLSCommon, only: DefaultUndefinedValue, MLSFile_T
  use MLSKinds, only: R8, R4
  use MLSL2Options, only: Default_HDFversion_Read, Default_HDFversion_Write
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Deallocate, &
    & MLSMSG_Error, MLSMSG_Warning
  use MLSSignals_M, only: GetModuleName, Modules
  use MLSStrings, only: Lowercase
  use MLSStringLists, only: Array2List, GetStringElement, List2Array, &
    & NumStringElements, StringElement
  use Output_M, only: Output
  use QuantityTemplates, only: QuantityTemplate_T
  use String_Table, only: Get_String, Display_String

  implicit none

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (data types and parameters)
! L2AUX_Dimension_T               Dimension for an L2AUX quantity
! L2AUXData_T                     An L2AUX quantity
! L2AUXRANK                       Rank (num of dims) of L2AUXData_T%values

!     (subroutines and functions)
! AddL2AUXToDatabase              Adds an l2aux data type to a database of that type
! cpL2AUXData                     Copies an l2aux quantity from file1 to file2
! DestroyL2AUXContents            Deallocates all the arrays for one l2aux
! DestroyL2AUXDatabase            Deallocates all the arrays for entire database
! Dump                            Prints info on one quantity or entire database
! ResizeL2AUXData                 Expands an l2aux quantity to take more profiles
! ReadL2AUXData                   Reads an l2aux quantity from a file
! SetupNewL2AUXRecord             Allocates the arrays for an l2aux quantity
! WriteHDF5Data                   Writes a named array to a file as an hdf dataset
! WriteL2AUXData                  Writes an l2aux quantity to a file
! WriteL2AUXAttributes            Writes l2aux sttributes to a file
! === (end of toc) ===

! === (start of api) ===
! int L2AUXRANK
!     (user-defined types)
! L2AUX_Dimension_T  ( int NoValues, int DIMENSIONFAMILY, *r8 values(:) )
! L2AUXData_T  ( int name, int INSTRUMENTMODULE,
!     log minorframe, log majorframe, L2AUX_Dimension_T dimensions(:),
!     *r8 values(:,:,:) )

!     (subroutines and functions)
! SetupNewL2AUXRecord ( int dimensionFamilies(L2AUXRANK), 
!    int dimSizes(L2AUXRANK), int dimStarts(L2AUXRANK), L2AUXData_T l2aux )
! CpL2AUXData ( MLSFile_t L2AUXFile1,  MLSFile_t L2AUXFile1, &
!    [log create], [char* sdList], [char* rename], [char* options] )
! DestroyL2AUXContents ( L2AUXData_T l2aux )
! ResizeL2AUXData ( L2AUXData_T l2aux, int newSize )
! int AddL2AUXToDatabase ( *L2AUXData_T DATABASE(:), L2AUXData_T ITEM )
! DestroyL2AUXDatabase ( *L2AUXData_T DATABASE(:) )
! Dump ( l2auxData_T L2aux(:), [char* Name], [int Details], [char* options] )
!    or Dump ( l2auxData_T L2aux, [int Details], [char* options] )
! ReadL2AUXData ( int sd_id, char* quantityname, l2auxData_T l2aux, 
!    [int firstProf], [int lastProf] )
! WriteHDF5Data ( real(r8) array(:,:,:), int l2FileHandle, int returnStatus, 
!    char* sdName )
! WriteL2AUXData ( l2auxData_T l2aux, int l2FileHandle, int returnStatus, 
!    [char* sdName], [int NoMAFS], [log WriteCounterMAF], [char* DimNames] )
! === (end of api) ===

  private
  public :: L2AUX_Dimension_T, L2AUXData_T, L2AUXRANK
  public :: AddL2AUXToDatabase, cpL2AUXData, DestroyL2AUXDatabase, Dump
  public :: SetupNewL2AUXRecord, DestroyL2AUXContents, ResizeL2AUXData
  public :: ReadL2AUXData, WriteHDF5Data, WriteL2AUXData, WriteL2AUXAttributes

  interface cpL2AUXData
    module procedure cpL2AUXData_Name
    module procedure cpL2AUXData_MLSFile
  end interface

  interface Dump
    module procedure Dump_L2AUX
    module procedure Dump_L2AUX_Database
  end interface

  interface ReadL2AUXData
    module procedure ReadL2AUXData_FileHandle
    module procedure ReadL2AUXData_MLSFile
  end interface

  interface WriteHDF5Data
    module procedure WriteHDF5Data_2d
    module procedure WriteHDF5Data_3d
  end interface

  interface WriteL2AUXData
    module procedure WriteL2AUXData_FileHandle
    module procedure WriteL2AUXData_MLSFile
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! This module defines datatypes and gives basic routines for storing and
  ! manipulating L2AUX data.
  ! L2AUX data is a rank three array, with the first index being a channel like
  ! index, the second a height like index, and the last an instance index.  It
  ! is directly conformable to a vector quantity.
  ! If a dimension is absent, it is flagged appropriately but the array is still
  ! rank 3.  However, when an HDF file is written or read missing dimensions
  ! are swallowed.

  ! For Announce_Error
    integer :: ERROR

  real, parameter    :: UNDEFINED_VALUE = DEFAULTUNDEFINEDVALUE ! -999.99 
  integer, parameter :: L2AUXRANK=3     ! Dimensionality of L2AUXData_T%values
  logical, parameter :: DEEBUG = .false.
  ! Write phase and section names as file-level attributes? (only for hdf5)
  logical, public, save :: PHASENAMEATTRIBUTES = .true.

  ! This datatype describes a dimension for an L2AUX quantity
  type L2AUX_Dimension_T
     integer :: NOVALUES = 0    ! Length of this dimension
     integer :: DIMENSIONFAMILY ! What is this dimension
     real(r8), dimension(:), pointer :: VALUES=>NULL() ! (noValues)
  end type L2AUX_Dimension_T

  ! This datatype describes an l2aux quantity itself.
  ! The dimensions will typically be ordered as follows:
  ! [Channel or frequency], MIF, [MAF or time or geodAngle]

  type L2AUXData_T
    integer :: NAME                     ! String index of name to be output
    integer :: INSTRUMENTMODULE = 0     ! From source vector
    integer :: QUANTITYTYPE = 0         ! From source vector
    logical :: MINORFRAME               ! Is this a minor frame quantity
    logical :: MAJORFRAME               ! Is this a major frame quantity
    ! type (QuantityTemplate_T) :: TEMPLATE = 0 ! Template for this l2aux quantity.
    ! The dimensions for the quantity
    type (L2AUX_Dimension_T), dimension(L2AUXRank) :: DIMENSIONS
    character(len=48)                              :: DIM_Names ! ','-separated
    character(len=48)                              :: DIM_Units ! ','-separated
    ! The values of the quantity
    real(r8), pointer, dimension(:,:,:)            :: VALUES=>NULL()
    character(len=24)                              :: VALUE_Units
  end type L2AUXData_T

  ! How long may the list of sd names grow (~80 x max num. of sds/file)
  integer, public, parameter :: MAXNUMSDPERFILE = 500
  integer, public, parameter :: MAXSDNAMESBUFSIZE = 80*MAXNUMSDPERFILE
  
  ! One quantity is recognized, but not copied, the other uncrecognized
  integer, private, parameter :: UNRECOGNIZEDQUANTITYTYPE = -999
  integer, private, parameter :: RECOGNIZEDBUTNOTCOPIEDQT = UNRECOGNIZEDQUANTITYTYPE + 1

contains ! =====     Public Procedures     =============================

  ! ------------------------------------------------- cpL2AUXData_MLSFile  -----

  subroutine cpL2AUXData_MLSFile( L2AUXFile1, L2AUXFile2, &
    & create2, sdList, rename, options )
    use Dump_1, only: Dump
    use MLSFiles, only: AreTheSameFile
    use MLSHDF5, only: GetAllHDF5DSNames
    !-------------------------------------------------------------------

    ! Given file names file1 and file2,
    ! This routine copies all the l2auxdata from 1 to 2
    ! If file2 doesn't exist yet, or if create2 is TRUE, it'll create it

    ! Arguments

    type(MLSFile_T), pointer      :: L2AUXfile1 ! file 1
    type(MLSFile_T), pointer      :: L2AUXfile2 ! file 2
    logical, optional, intent(in) :: create2 ! Force creation of new file2
    character (len=*), optional, intent(in) :: sdList ! Copy only these, unless
    character (len=*), optional, intent(in) :: rename
    character (len=*), optional, intent(in) :: options ! E.g., '-v'

    ! Local
    logical :: allSDs
    logical, parameter            :: countEmpty = .true.
    integer :: i
    type (L2AUXData_T) :: l2aux
    logical :: myCreate2
    character (len=8) :: myOptions
    character (len=MAXSDNAMESBUFSIZE) :: mySdList
    integer :: noSds
    integer :: originalAccess
    integer :: QuantityType
    logical :: renameSwaths
    character (len=80) :: sdName
    integer :: status
    logical :: verbose

    myOptions = ' '
    if ( present(options) ) myOptions = options
    verbose = ( index(myOptions, 'v') > 0 )
    renameSwaths  = present(rename)
    if ( renameSwaths ) renameSwaths = ( rename /= ' ' )
    myCreate2 = .false.
    if ( present(create2) ) myCreate2 = create2
    originalAccess = L2AUXFile2%access
    if ( myCreate2 ) L2AUXFile2%access = DFACC_CREATE
    allSDs = .not. present(sdList)
    if ( present(sdList) ) allSDs = (sdList == '*')
    if ( .not. allSDs ) then
      mysdList = sdList
      if ( verbose ) call dump(mysdList, 'DS names')
    else
      call GetAllHDF5DSNames (trim(L2AUXFile1%Name), '/', mysdList)
      if ( verbose ) call output ( '============ DS names in ', advance='no' )
      if ( verbose ) call output ( trim(L2AUXfile1%Name) //' ============', &
        & advance='yes' )
      if ( mysdList == ' ' ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'No way yet to find sdList ', MLSFile=L2AUXFile1 )
        return
      else
        if ( verbose ) call dump(mysdList, 'DS names')
      end if
    end if

    if ( AreTheSameFile( L2AUXfile1, L2AUXfile2 ) ) &
      & call MLSMessage ( MLSMSG_Error, trim(ModuleName) // ' cpL2AUXData', &
      & 'input and output files are the same', MLSFile=L2AUXFile1 )
    noSds = NumStringElements(trim(mysdList), countEmpty)
    if ( noSds < 1 ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'No sdNames cp to file--unable to count sdNames in ' // trim(mysdList) )
    end if
    ! Loop over sdNames in file 1
    do i = 1, noSds
      call GetStringElement (trim(mysdList), sdName, i, countEmpty )
        QuantityType = GetQuantityTypeFromName(trim(sdName)) ! l_radiance
      if ( QuantityType == RECOGNIZEDBUTNOTCOPIEDQT ) then
        cycle
      elseif ( QuantityType < 1 ) then
        call output('Quantity type: ', advance='no')
        call output(QuantityType, advance='yes')
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
              & 'Unrecognized quantity type for sd:' // trim(sdName) )
        cycle
      end if
      call ReadL2AUXData ( L2AUXFile1, trim(sdName), l2aux, QuantityType, &
           & checkDimNames=.false. )
      ! Write the filled l2aux to file2
      ! print *, 'writing ', trim(sdName)
      if ( renameSwaths ) sdName = StringElement ( rename, i, countEmpty )
      call WriteL2AUXData( l2aux, L2AUXFile2, status, trim(sdName) )
      ! Deallocate memory used by the l2aux
      call DestroyL2AUXContents ( l2aux )
      L2AUXFile2%access = originalAccess
    end do
  end subroutine cpL2AUXData_MLSFile

  ! ------------------------------------------------- cpL2AUXData_Name  -----

  subroutine cpL2AUXData_Name(file1, file2, create2, hdfVersion, sdList, rename, &
    & options)
    use Dump_1, only: Dump
    use HDF, only: Dfacc_Read, Dfacc_Rdwr
    use HDF5, only: H5gclose_F, H5gopen_F, H5dopen_F, H5dclose_F
    use MLSFiles, only: FileNotFound, WildCardHDFVersion, &
      & MLS_Exists, MLS_HDF_Version, MLS_Sfstart, MLS_Sfend
    use MLSHDF5, only: GetAllHDF5DSNames, GetHDF5Attribute, &
      & IsHDF5AttributePresent
    !-------------------------------------------------------------------

    ! Given file names file1 and file2,
    ! This routine copies all the l2auxdata from 1 to 2
    ! If file2 doesn't exist yet, or if create2 is TRUE, it'll create it

    ! Arguments

    character (len=*), intent(in) :: file1 ! Name of file 1
    character (len=*), intent(in) :: file2 ! Name of file 2
    logical, optional, intent(in) :: create2 ! Force creation of new file2
    integer, optional, intent(in) :: hdfVersion       !                      '*'
    character (len=*), optional, intent(in) :: sdList ! Copy only these, unless
    character (len=*), optional, intent(in) :: rename
    character (len=*), optional, intent(in) :: options ! E.g., '-v'

    ! Local
    logical :: allSDs
    logical, parameter            :: countEmpty = .true.
    logical :: file_exists
    integer :: file_access
    integer :: grpid
    integer :: i
    type (L2AUXData_T) :: l2aux
    character (len=8) :: myOptions
    character (len=MAXSDNAMESBUFSIZE) :: mySdList
    integer :: noSds
    integer :: QuantityType
    logical :: renameSwaths
    integer :: sd_id
    integer :: sdfid1
    integer :: sdfid2
    character (len=80) :: sdName
    integer :: status
    integer :: the_hdfVersion
    logical :: verbose
    
    ! Executable code
    the_hdfVersion = DEFAULT_HDFVERSION_WRITE
    if ( present(hdfVersion) ) the_hdfVersion = hdfVersion
    myOptions = ' '
    if ( present(options) ) myOptions = options
    verbose = ( index(myOptions, 'v') > 0 )
    renameSwaths  = present(rename)
    if ( renameSwaths ) renameSwaths = ( rename /= ' ' )
    file_exists = ( mls_exists(trim(File1)) == 0 )
    if ( .not. file_exists ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'File 1 not found; make sure the name and path are correct' &
        & // trim(file1) )
    end if
    if ( the_hdfVersion == WILDCARDHDFVERSION ) then
      the_hdfVersion = mls_hdf_version(File1, hdfVersion)
      if ( the_hdfVersion == FILENOTFOUND ) &
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'File 1 not found; make sure the name and path are correct' &
          & // trim(file1) )
    end if
    allSDs = .not. present(sdList)
    if ( present(sdList) ) allSDs = (sdList == '*')
    if ( .not. allSDs ) then
      mysdList = sdList
      if ( verbose ) call dump(mysdList, 'DS names')
    else
      call GetAllHDF5DSNames (trim(File1), '/', mysdList)
      if ( verbose ) call output ( '============ DS names in ', advance='no' )
      if ( verbose ) call output ( trim(file1) //' ============', advance='yes' )
      if ( mysdList == ' ' ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'No way yet to find sdList in ' // trim(File1) )
        return
      else
        if ( verbose ) call dump(mysdList, 'DS names')
      end if
    end if

    file_exists = ( mls_exists(trim(File2)) == 0 )
    if ( file_exists ) then
      file_access = DFACC_RDWR
    else
      file_access = DFACC_CREATE
    end if
    if ( present(create2) ) then
      if ( create2 ) file_access = DFACC_CREATE
    end if
    sdfid1 = mls_sfstart(File1, DFACC_READ, hdfVersion=hdfVersion)
    if (sdfid1 == -1 ) then
      call announce_error ( 0, 'Failed to open l2aux ' // &
      &  trim(File1) )
    end if
      call h5gOpen_f (sdfid1,'/', grpID, status)
    if ( status /= 0 ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
           & 'Unable to open group to read attribute in l2aux file' )
    end if
    sdfId2 = mls_sfstart(trim(file2), file_access, &
              & hdfVersion=hdfVersion)
    if (sdfid2 == -1 ) then
      call announce_error ( 0, 'Failed to open l2aux ' // &
      &  trim(File2) )
    end if
    noSds = NumStringElements(trim(mysdList), countEmpty)
    if ( noSds < 1 ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'No sdNames cp to file--unable to count sdNames in ' // trim(mysdList) )
    end if
    ! Loop over sdNames in file 1
    do i = 1, noSds
      call GetStringElement (trim(mysdList), sdName, i, countEmpty )
      ! Allocate and fill l2aux
      call h5dOpen_f (grpid,trim(sdName), sd_ID, status)
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
              & 'Unable to open sd to read attribute in l2aux file' )
      end if
      ! Get QuantityType attribute--unfortunately they're all 0; what gives?
      if ( .not. IsHDF5AttributePresent(sd_id, 'QuantityType') .or. .true.) then
        QuantityType = GetQuantityTypeFromName(trim(sdName)) ! l_radiance
      else
        call GetHDF5Attribute ( sd_id, 'QuantityType', QuantityType )
      end if
      call h5dClose_f (sd_ID, status)
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
              & 'Unable to close sd to read attribute in l2aux file' )
      end if
      if ( QuantityType < 1 ) then
        call output('Quantity type: ', advance='no')
        call output(QuantityType, advance='yes')
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
              & 'Unrecognized quantity type for sd:' // trim(sdName) )
        cycle
      end if
      call ReadL2AUXData ( sdfid1, trim(sdName), l2aux, QuantityType, &
           & checkDimNames=.false., hdfVersion=hdfVersion )
      if ( renameSwaths ) sdName = StringElement ( rename, i, countEmpty )
      ! Write the filled l2aux to file2
      call WriteL2AUXData(l2aux, sdfid2, status, trim(sdName), &
        & hdfVersion=hdfVersion)
      ! Deallocate memory used by the l2aux
      call DestroyL2AUXContents ( l2aux )
    end do
    call h5gClose_f (grpID, status)
    if ( status /= 0 ) then
    call MLSMessage ( MLSMSG_Warning, ModuleName, &
       & 'Unable to close group in l2aux file: ' // trim(File1) // ' after cping' )
    end if
    status = mls_sfend(sdfid1, hdfVersion=the_hdfVersion)
    if ( status /= 0 ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
       & "Unable to close L2aux file: " // trim(File1) // ' after cping')
    status = mls_sfend(sdfid2, hdfVersion=the_hdfVersion)
    if ( status /= 0 ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
       & "Unable to close L2aux file: " // trim(File2) // ' after cping')
  end subroutine cpL2AUXData_Name

  ! ---------------------------------------  SetupNewL2AUXRecord   -----
  !   (option 1)
  ! subroutine SetupNewL2AUXRecord ( l2aux, &
  !  & quantityTemplate, firstMAF, noMAFs )
  !   (option 2)
  ! subroutine SetupNewL2AUXRecord ( l2aux, dimFamilies, dimSizes, dimStarts, &
  !  & quantityType )
  subroutine SetupNewL2AUXRecord ( l2aux, &
   & quantityTemplate, firstMAF, noMAFs, &
   & inputDimFamilies, inputDimSizes, inputDimStarts, inputQuantityType )

    ! This first routine sets up the arrays for an l2aux datatype.
    ! (Option 1)
    ! The user supplies a Quantity template, 1st MAF and number of MAFs
    ! Then we deduce the dimension families and their sizes
    ! (Option 2)
    ! The user supplies a set of three dimensionFamilies (e.g. l_maf)
    ! and their sizes and starts
    ! Quantities can have up to three valid dimensions.  l_none can be used
    ! to indicate later dimensions are invalid.

    ! Dummy arguments
    type (L2AUXData_T), intent(out)       :: l2aux
    !    ( option 1 )
    type (QuantityTemplate_T), intent(in), optional :: quantityTemplate
    integer, intent(in), optional                   :: firstMAF
    integer, intent(in), optional                   :: noMAFs
    !    ( option 2 )
    integer, dimension(L2AUXRank), optional, intent(in) :: &
                                                    & inputDimFamilies
    integer, dimension(L2AUXRank), optional, intent(in) :: inputDimSizes
    integer, dimension(L2AUXRank), optional, intent(in) :: inputDimStarts
    integer, optional, intent(in) :: inputQuantityType
    ! Local variables
    integer, dimension(L2AUXRank) :: dimSizes
    integer, dimension(L2AUXRank) :: dimStarts
    integer             :: quantityType
    integer :: dimIndex
    integer :: status
    integer, dimension(L2AUXRank)   :: dimEnds
    integer, dimension(3        )   :: dim_names
    character(len=16), dimension(3) :: extra_name
    character(len=16)               :: framing
    ! integer                         :: option_number   ! (1 or 2; see above)

    ! Executable
    if ( present(quantityTemplate) ) then
      ! option_number = 1
      quantityType = quantityTemplate%quantityType
    else if ( present(inputQuantityType) ) then
      ! option_number = 1
      quantityType = inputQuantityType
      dimSizes = inputDimSizes
      dimStarts = inputdimStarts
    else
      call MLSMessage(MLSMSG_Error, ModuleName, &
      & 'args to SetupL2AUXData incompatible with options 1 or 2')
    end if
    ! Fill the dimensions data structure
    call GetQuantityAttributes( quantityType, &
       & framing, l2aux%VALUE_Units, dim_names )
    l2aux%dimensions%dimensionFamily = dim_names
    ! end if
    l2aux%minorFrame = (framing == 'minor')
    l2aux%majorFrame = (framing == 'major')
    ! Name the dimensions (e.g. 'frequency')
    do dimIndex=1, L2AUXRank
      call GetDimString( l2aux%dimensions(dimIndex)%dimensionFamily, &
        & extra_name(dimIndex) )
      if ( present (firstMAF) ) then
        call GetDimStart( l2aux%dimensions(dimIndex)%dimensionFamily, &
        & quantityTemplate, firstMAF, dimStarts(dimIndex) )
      end if
      if ( present (noMAFs) ) then
        call GetDimSize( l2aux%dimensions(dimIndex)%dimensionFamily, &
        & quantityTemplate, noMAFs, dimSizes(dimIndex) )
      end if
    end do
    call Array2List(extra_name, l2aux%DIM_Names)
    ! Name the dimensions' units (e.g. 'K')
    do dimIndex=1, L2AUXRank
      call GetQuantityAttributes( l2aux%dimensions(dimIndex)%dimensionFamily, &
       & units_name=extra_name(dimIndex) )
    end do
    call Array2List(extra_name, l2aux%DIM_Units)
    l2aux%dimensions%noValues = dimSizes

    dimEnds = dimStarts + max(1,dimSizes) - 1

    ! Allocate the values for each dimension
    do dimIndex = 1, L2AUXRank
      if ( l2aux%dimensions(dimIndex)%dimensionFamily /= L_None ) then
        allocate (l2aux%dimensions(dimIndex)%values( &
          & dimStarts(dimIndex):dimEnds(dimIndex)), &
          & STAT=status)
        if ( status/=0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
          & MLSMSG_Allocate // "l2aux dimension values" )
        if ( present (quantityTemplate) ) then
          call GetDimValues( l2aux%dimensions(dimIndex)%dimensionFamily, &
          & quantityTemplate, l2aux%dimensions(dimIndex)%values )
        end if
      else
        l2aux%dimensions(dimIndex)%noValues=1
      end if
    end do

    ! Allocate the values for the data itself

    nullify ( l2aux%values )
    call allocate_test( l2aux% values, &
      & dimEnds(1), &
      & dimEnds(2), &
      & dimEnds(3), &
      & 'l2aux%values', ModuleName, &
      & dimStarts(1), &
      & dimStarts(2), &
      & dimStarts(3) &
      )
  end subroutine SetupNewL2AUXRecord
    
  !----------------------------------------  DestroyL2AUXContents  -----
  subroutine DestroyL2AUXContents ( l2aux )

  ! This routine deallocates all the arrays allocated above.

    ! Dummy arguments
    type (L2AUXData_T), intent(inout) :: l2aux
    ! Local variables
    integer :: dim
    ! Executable code
    
    do dim=1,L2AUXRank
      if (l2aux%dimensions(dim)%dimensionFamily /= L_None) &
        call deallocate_test(l2aux%dimensions(dim)%values, &
        "l2aux%dimensions", ModuleName)
    end do
    call deallocate_test(l2aux%values,"l2aux%values",ModuleName ) 
  end subroutine DestroyL2AUXContents

  !--------------------------------------  ResizeL2AUXData  -----
  subroutine ResizeL2AUXData ( l2aux, newSize )

  ! This subroutine resizes an L2AUXData_T in place, allowing the user to
  ! add more `profiles' to it.  Note that the `profile' dimension is the last
  ! one.
  
  ! New feature: if newSize is < old size
  ! we effectively contract the l2AUXData_T

    ! Dummy arguments
    type (L2AUXData_T), intent(inout) :: L2aux
    integer, intent(in) :: NewSize

    ! Local variables
    integer :: OldSize
    ! The following are temporary arrays for copying data around
    real (r8), dimension(:), pointer :: Temp1D
    real (r8), dimension(:,:,:), pointer :: Temp3D

    ! Executable code

    nullify ( temp1D, temp3D )

    if ( l2aux%dimensions(3)%dimensionFamily == L_None ) &
      call MLSMessage (MLSMSG_Error, ModuleName, &
        & "This l2aux is not expandable") ! Why not?

    ! Now see how long this is
    oldSize = l2aux%dimensions(3)%noValues
    ! Do a usage check
    if ( newSize < oldSize ) call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & "This l2aux is getting smaller not bigger" )

    ! Now expand this dimension
    temp1D => l2aux%dimensions(3)%values
    ! Nullify old one so allocate_test doesn't clobber temp1D
    nullify ( l2aux%dimensions(3)%values )
    call allocate_test(l2aux%dimensions(3)%values,newSize, &
      & 'New l2aux%dimensions(3)%values', ModuleName)
    l2aux%dimensions(3)%noValues = newSize
    if ( oldSize < newSize ) then
      l2aux%dimensions(3)%values(1:oldSize) = temp1D
    else
      l2aux%dimensions(3)%values = temp1D(1:newSize)
    end if
    ! Now we can loose the old values
    call deallocate_test ( temp1D, "temp1D", ModuleName )

    ! Now expand the data in this dimension, save old field
    temp3d => l2aux%values
    ! Need to nullify old one, else next call will clobber it
    ! loosing us temp3d
    nullify ( l2aux%values )
    call allocate_test ( l2aux%values, &
      & max(1, l2aux%dimensions(1)%noValues), &
      & max(1, l2aux%dimensions(2)%noValues), &
      & newSize, "l2aux%values", ModuleName )
    if ( oldSize < newSize ) then
      l2aux%values(:,:,1:oldSize) = temp3d
    else
      l2aux%values = temp3d(:,:,1:newSize)
    end if

    ! Now we can set temp3d loose
    call deallocate_test ( temp3d, "temp3d", ModuleName )

  end subroutine ResizeL2AUXData

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

    AddL2AUXToDatabase = newSize
  end function AddL2AUXToDatabase

  ! ---------------------------------------  DestroyL2AUXDatabase  -----
  subroutine DestroyL2AUXDatabase ( DATABASE )

  ! This subroutine destroys the l2aux database

    use Toggles, only: Gen, Toggle
    use Trace_m, only: Trace_Begin, Trace_End

    ! Dummy argument
    type (L2AUXData_T), dimension(:), pointer :: DATABASE

    ! Local variables
    integer :: l2auxIndex, status
    integer :: Me = -1       ! String index for trace

    call trace_begin ( me, "DestroyL2AUXDatabase", cond=toggle(gen) )

    if ( associated(database) ) then
      do l2auxIndex = 1, size(database)
        call DestroyL2AUXContents ( database(l2auxIndex) )
      end do
      deallocate ( database, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Warning,ModuleName, &
        & MLSMSG_DeAllocate // "database")
    end if
    call trace_end ( "DestroyL2AUXDatabase", cond=toggle(gen) )
  end subroutine DestroyL2AUXDatabase

  ! ----------------------------------------  Dump_L2AUX_DataBase  -----

  subroutine Dump_L2AUX_DataBase ( L2aux, Name, Details, Options )
    use HighOutput, only: StyledOutput

    ! Dummy arguments
    type (l2auxData_T), intent(in) ::          L2AUX(:)
    character(len=*), intent(in), optional :: Name
    integer, intent(in), optional :: DETAILS
    character(len=*), intent(in), optional :: options

    ! Local variables
    integer :: i
    
    call output ( '============ L2AUX Data Base ============', advance='yes' )
    call output ( ' ', advance='yes' )
    if ( present(name) ) then
      call StyledOutput ( name, options )
    end if
    if ( size(l2aux) < 1 ) then
      call output ( '**** L2AUX Database empty ****', advance='yes' )
      return
    end if
    do i = 1, size(l2aux)
      call dump( l2aux(i), Details, options )
    end do
      
  end subroutine Dump_L2AUX_DATABASE

  ! -------------------------------------------------  Dump_L2AUX  -----

  subroutine Dump_L2AUX ( L2aux, Details, Options )

    ! Dummy arguments
    type (l2auxData_T), intent(in) ::          L2AUX
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump multidim arrays
    !                                        ! -1 Skip even 1-d arrays
    !                                        ! -2 Skip all but name
    !                                        ! >0 Dump even multi-dim arrays
    !                                        ! Default 1
    character(len=*), intent(in), optional :: options

    ! Local variables
    integer :: dim, ierr
    integer :: MYDETAILS

    ! Executable code
    myDetails = 1
    if ( present(details) ) myDetails = details
    
    call output ( 'L2AUX Data: ')
    call display_string ( l2aux%name, ierr=ierr )
    if ( ierr /= 0 ) call output ( '(not found in string table)')
    if ( myDetails < -1 ) return
    if ( associated(modules) ) then
      call output ( '    instrumentmodule: ')
      call display_string ( modules(l2aux%instrumentmodule)%name, &
        & advance='no', ierr=ierr ) 
      if ( ierr /= 0 ) then
        call output ( ' (not found in string table)', advance='no')
      else
        call output ( l2aux%instrumentmodule, before=' = ', advance='no')
      end if
    end if
    call output ( '', advance='yes')
    call output ( '  Minor Frame? (t/f): ')
    call output ( l2aux%minorframe, advance='no')
    call output ( '  Major Frame? (t/f): ')
    call output ( l2aux%majorframe, advance='yes')
    if ( myDetails < 0 ) return
    do dim=1, l2auxrank
      call output ( '  dimension: ')
      call output ( dim )
      call output ( '           ')
      if ( associated(l2aux%dimensions(dim)%values) ) then
        call output ( '  nValues: ')
        call output ( l2aux%dimensions(dim)%novalues, 3, advance='no')
        call output ( '           ')
        call output ( '  dimension family: ')
        call output ( l2aux%dimensions(dim)%dimensionfamily, 3, advance='yes')
        call dump ( l2aux%dimensions(dim)%values, &
          & 'dim values:', options=options )
       else
        call output ( ' is not associated', advance='yes')
       end if
    end do
    if ( myDetails < 1 ) return
    call dump ( l2aux%values, 'values:', options=options )
 
  end subroutine Dump_L2AUX
    
  !-----------------------------------------------  ReadL2AUXData_FileHandle  -----
  subroutine ReadL2AUXData_FileHandle( sd_id, quantityname, l2aux, inQuantityType, &
    & firstProf, lastProf, checkDimNames, hdfVersion )

  use MLSFiles, only: InitializeMLSFile

    ! This routine reads an l2aux file, returning a filled data structure
    ! and the number of profiles read.

    ! Arguments

    character (len=*), intent(in) :: quantityname ! Name of L2AUX quantity = sdname in writing routine
    integer, intent(in) :: sd_id ! Returned by sfstart before calling us
    integer, intent(in), optional :: inQuantityType ! Lit index
    integer, intent(in), optional :: firstProf, lastProf ! Defaults to first and last
    type( L2AUXData_T ), intent(OUT) :: l2aux ! Result
    logical, optional, intent(in) :: checkDimNames
    integer, intent(in), optional :: hdfVersion
    ! Local variables
    integer :: myhdfVersion
    integer :: QuantityType
    type( MLSFile_T ), target  :: l2auxFile
    type( MLSFile_T ), pointer :: l2auxPointer
    integer :: status
    ! Executable code
    myhdfVersion = default_hdfversion_read
    if ( present(hdfVersion) ) myhdfVersion = hdfVersion
    if ( present(inQuantityType) ) then
      QuantityType = inQuantityType
    else
      QuantityType = GetQuantityTypeFromName(trim(quantityname))
    endif
    status = InitializeMLSFile(l2auxFile, type=l_hdf, access=DFACC_RDONLY, &
      & content='l2aux', name='unknown', hdfVersion=myhdfVersion)
    l2auxFile%FileID%f_id = sd_id
    l2auxFile%stillOpen = .true.
    l2auxPointer => l2auxFile
    if ( deebug ) print *, 'About to call ReadL2AUXData_MLSFile'
    call ReadL2AUXData( l2auxPointer, quantityname, l2aux, &
      & quantityType, firstProf, lastProf, checkDimNames )

  end subroutine ReadL2AUXData_FileHandle

  !-----------------------------------------------  ReadL2AUXData_MLSFile  -----
  subroutine ReadL2AUXData_MLSFile( L2AUXFile, quantityname, l2aux, inQuantityType, &
    & firstProf, lastProf, checkDimNames )

 use MLSFiles, only: Dump, HDFVersion_4, HDFVersion_5, &
   & MLS_CloseFile, MLS_OpenFile
 use Trace_M, only: Trace_Begin, Trace_End

    ! This routine reads an l2aux file, returning a filled data structure
    ! and the number of profiles read.

    ! Arguments

    character (len=*), intent(in) :: quantityname ! Name of L2AUX quantity = sdname in writing routine
    type(MLSFile_T), pointer      :: L2AUXFile
    integer, intent(in), optional :: inQuantityType ! Lit index
    integer, intent(in), optional :: firstProf, lastProf ! Defaults to first and last
    type( L2AUXData_T ), intent(OUT) :: l2aux ! Result
    logical, optional, intent(in) :: checkDimNames

    ! Local variables
    logical :: alreadyOpen
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: returnStatus
    integer :: QuantityType

    ! Executable code
    call trace_begin ( me, 'ReadL2AUXData_MLSFile', cond=.false. )
    returnStatus = 0
    if ( present(inQuantityType) ) then
      QuantityType = inQuantityType
    else
      QuantityType = GetQuantityTypeFromName(trim(quantityname))
    endif
    alreadyOpen = L2AUXFile%stillOpen
    if ( L2AUXFile%access == DFACC_CREATE ) then
      call MLSMessage(MLSMSG_Error, trim(ModuleName) // ' ReadL2AUXData_MLSFile', &
        & 'Attempt to open l2aux file for reading with create access', MLSFile=L2AUXFile)
    end if
    if ( .not. alreadyOpen ) then
      call mls_openFile(L2AUXFile, returnStatus)
      if ( returnStatus /= 0 ) &
        call MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to open l2aux file for reading', MLSFile=L2AUXFile)
    end if
    if ( deebug ) print *, 'Must choose to read based on ', L2AUXFile%hdfVersion
    if ( deebug ) call Dump ( L2AUXFILE )
    select case (L2AUXFile%hdfVersion)
    case (HDFVERSION_4)
      call ReadL2AUXData_MF_hdf4(L2AUXFile, quantityname, quantityType, l2aux, &
       & firstProf, lastProf, checkDimNames)
    case (HDFVERSION_5)
      call ReadL2AUXData_MF_hdf5(L2AUXFile, quantityname, quantityType, l2aux, &
        & firstProf, lastProf, checkDimNames)
    case default
    end select
    if ( .not. alreadyOpen )  call mls_closeFile(L2AUXFile, returnStatus)

    L2AUXFile%errorCode = returnStatus
    L2AUXFile%lastOperation = 'read'
    call trace_end ( 'ReadL2AUXData_MLSFile', cond=.false. )
  end subroutine ReadL2AUXData_MLSFile

  ! -----------------------------------------  ReadL2AUXData_MF_hdf4  -----
  subroutine ReadL2AUXData_MF_hdf4( L2AUXFile, quantityname, quantityType, l2aux, &
    & firstProf, lastProf, checkDimNames )

    ! This routine reads an l2aux file, returning a filled data structure and the
    ! number of profiles read.

    ! Arguments

    ! Name of L2AUX quantity = sdname in writing routine
    character (len=*), intent(in) :: quantityname 
    type(MLSFile_T), pointer      :: L2AUXFile
    integer, intent(in) :: QuantityType
    integer, intent(in), optional :: firstProf, lastProf ! Defaults to first and last
    type( L2AUXData_T ), intent(OUT) :: l2aux ! Result
    logical, optional, intent(in) :: checkDimNames

    ! Parameters

    integer, parameter :: MAXRANK = 3
    logical            :: myCHECKDIMNAMES ! .TRUE. only for actual l2auxfiles

    ! Variables

    character (LEN=480) :: msr

    integer :: sds_index, sds_id, rank, data_type, num_attrs, dim, dim_id
    integer :: data_dim_sizes(MAXRANK), file_dim_sizes(MAXRANK), dim_size1
    integer :: dim_families(MAXRANK)
    character (LEN=len(quantityname)) :: sds_name
    character (LEN=132) :: dim_name
    character (LEN=1)                  :: dim_char

    integer :: status
    integer :: start(3), stride(3)

    ! logical :: firstCheck, lastCheck
    real (r4), dimension(:,:,:), pointer :: TMPVALUES

    myCHECKDIMNAMES = .false.
    if ( present(checkDimNames) ) myCHECKDIMNAMES = checkDimNames
    ! Attach to the file for reading

    ! find SD data set identifier
    sds_index = sfn2index(L2AUXFile%FileID%f_id, quantityname)
    if (sds_index == -1) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Failed to get sds_index for '//trim(quantityName), MLSFile=L2AUXFile )

    sds_id = sfselect(L2AUXFile%FileID%f_id, sds_index)
    if (sds_id == -1) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Failed to get sds_id.', MLSFile=L2AUXFile )

    status = sfginfo(sds_id, sds_name, rank, file_dim_sizes, data_type, &
      & num_attrs)

    if (status == -1) then
      call MLSMessage(MLSMSG_Error, ModuleName, 'Failed to get sf info.', MLSFile=L2AUXFile)
    else if (sds_name /= quantityname) then
      call MLSMessage(MLSMSG_Error, ModuleName, &
        & 'quantityname  fails to match sf info.', MLSFile=L2AUXFile)
    end if

    ! Check optional input arguments

    ! firstCheck = present(firstProf)
    ! lastCheck = present(lastProf)

    ! Uncertain what to do with those just yet
    ! Now find dimension family of dimension; e.g., MAF
    dim_families = l_none
    data_dim_sizes = 1
    file_dim_sizes = 1

    do dim=1, rank
      write(dim_char, '(I1)') dim
      dim_id = sfdimid(sds_id, dim-1)  ! dim starts from 0
      if(dim_id == -1) then
        msr = 'Failed to get dim_id for dim index number ' // dim_char
        call MLSMessage(MLSMSG_Error, ModuleName, msr, MLSFile=L2AUXFile)
      else
        status = sfgdinfo(dim_id, dim_name, dim_size1, data_type, &
          & num_attrs)
        if(status == -1) then
          msr = 'Failed to get dim_info for dim index number ' // &
            & dim_char
          call MLSMessage(MLSMSG_Error, ModuleName, msr, MLSFile=L2AUXFile)
        else
          file_dim_sizes(dim) = dim_size1
          select case ( trim(dim_name) )
          case ( 'channel' )
            dim_families(1) = l_channel
            data_dim_sizes(1) = dim_size1
          case ( 'GHz.mif' )
            dim_families(2) = l_mif
            data_dim_sizes(2) = dim_size1
          case ( 'GHz.maf' )
            dim_families(3) = l_maf
            data_dim_sizes(3) = dim_size1
          case default
            if ( myCHECKDIMNAMES ) then
              call MLSMessage ( MLSMSG_Error, ModuleName, &
                & 'Unrecognized dimension in l2aux:'//trim(dim_name), MLSFile=L2AUXFile )
            else
              dim_families(dim) = l_channel
              data_dim_sizes(dim) = dim_size1
            end if
          end select
        end if
      end if
    end do

    ! Allocate result
!   call SetupNewl2auxRecord ( dim_families, data_dim_sizes, (/1,1,1/), l2aux )
    call SetupNewl2auxRecord ( l2aux, inputDimFamilies=dim_families, &
      & inputDimSizes=data_dim_sizes, inputDimStarts=(/1,1,1/), inputQuantityType=quantityType )

    ! Read the SD
    start = 0
    stride = 1
    nullify ( tmpValues )
    call Allocate_test ( tmpValues, max(file_dim_sizes(1),1), &
      & max(file_dim_sizes(2),1), max(file_dim_sizes(3),1), &
      & 'tmpValues', ModuleName )

    ! Not sure this isn't cheating with the array but we'll see
    status = sfrdata_f90(sds_id, start(1:rank), stride(1:rank), file_dim_sizes(1:rank), &
      & tmpValues )
    if (status == -1) call MLSMessage(MLSMSG_Error, ModuleName, &
      & 'Failed to read SD.', MLSFile=L2AUXFile)
    l2aux%values = reshape ( tmpValues, &
      & (/ data_dim_sizes(1), data_dim_sizes(2), data_dim_sizes(3) /) )

    call Deallocate_test ( tmpValues, 'tmpValues', ModuleName )

    ! Deallocate local variables


    ! Terminate access to the data set

    status = sfendacc(sds_id)
    if (status == -1) call MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
      &end access to sds_id after reading.', MLSFile=L2AUXFile)

  end subroutine ReadL2AUXData_MF_hdf4


  ! -----------------------------------------  ReadL2AUXData_MF_hdf5  -----
  subroutine ReadL2AUXData_MF_hdf5( L2AUXFile, quantityname, quantityType, &
    & L2AUX, firstProf, lastProf, checkDimNames )
    use HighOutput, only: BeVerbose
    use L1BData, only: L1BData_T, ReadL1BData

    ! This routine reads an l2aux file, returning a filled data structure and the
    ! number of profiles read.
    ! Assumptions
    ! The data format is radiance-like, single- or double-precision
    ! You don't really care which names the original dimensions had
    ! You want the new dimension families to be 'channel', 'MIF', 'MAF'
    ! and in that order

    ! Arguments

    ! Name of L2AUX quantity = sdname in writing routine
    character (len=*), intent(in) :: quantityname 
    type(MLSFile_T), pointer      :: L2AUXFile
    integer, intent(in) :: QUANTITYTYPE ! Lit index
    integer, intent(in), optional :: firstProf, lastProf ! Defaults to first and last
    type( L2AUXData_T ), intent(OUT) :: l2aux ! Result
    logical, optional, intent(in) :: checkDimNames

    ! Parameters
    type(l1bdata_t)               :: L1BDATA1 ! Intermediate Result
    integer                       :: NoMAFs
    logical, parameter            :: NEVERFAIL = .TRUE.
    integer, dimension(L2AUXRank) :: dim_families
    integer, dimension(L2AUXRank) :: data_dim_sizes
    integer                       :: status
    ! Executable
    if ( BeVerbose( 'l2aux', 0 ) ) &
      & call output( 'Attempting to read ' // trim(quantityname), advance='yes' )
    CALL ReadL1BData( L2AUXFile, QuantityName, L1BDATA1, NoMAFs, status, &
      & FirstMAF=firstProf, LastMAF=lastProf, NEVERFAIL=NEVERFAIL, &
      & dontPad=.true., L2AUX=.true. )
    if ( status /= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read ' &
      & // trim(QuantityName) // ' (perhaps too unlike a radiance)', MLSFile=L2AUXFile )

    dim_families(1) = l_channel
    select case ( L1BDATA1%data_type(1:1) )
    case ( 'd' )
      data_dim_sizes = shape(L1BDATA1%DpField)
      if ( deebug ) print *, 'shape(L1BDATA1%DpField): ', shape(L1BDATA1%DpField)
      dim_families(2) = l_mif                          
      dim_families(3) = l_maf                          
      call SetupNewl2auxRecord ( l2aux, inputDimFamilies=dim_families, &
       & inputDimSizes=data_dim_sizes, inputDimStarts=(/1,1,1/), &
       & inputQuantityType=quantityType )
      l2aux%values = L1BDATA1%DpField
      deallocate( L1BDATA1%DpField, stat=status )
    case ( 'i' )
      data_dim_sizes = shape(L1BDATA1%IntField)
      if ( deebug ) print *, 'shape(L1BDATA1%IntField): ', shape(L1BDATA1%IntField)
      dim_families(2) = l_mif                          
      dim_families(3) = l_maf                          
      call SetupNewl2auxRecord ( l2aux, inputDimFamilies=dim_families, &
       & inputDimSizes=data_dim_sizes, inputDimStarts=(/1,1,1/), &
       & inputQuantityType=quantityType )
      l2aux%values = L1BDATA1%IntField
      deallocate( L1BDATA1%IntField, stat=status )
    case ( 'c' )
      data_dim_sizes = shape(L1BDATA1%CharField)
      if ( deebug ) print *, 'shape(L1BDATA1%CharField): ', shape(L1BDATA1%CharField)
      dim_families(2) = l_mif                          
      dim_families(3) = l_maf                          
      call SetupNewl2auxRecord ( l2aux, inputDimFamilies=dim_families, &
       & inputDimSizes=data_dim_sizes, inputDimStarts=(/1,1,1/), &
       & inputQuantityType=quantityType )
      l2aux%values = iachar(L1BDATA1%CharField(:,:,:)(1:1))
      deallocate( L1BDATA1%CharField, stat=status )
    case default
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to cope with type l1bdata%' // trim(L1BDATA1%data_type) // &
      &  'while reading' // &
      & trim(QuantityName), MLSFile=L2AUXFile )
    end select
    if ( status /= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to deallocate l1bdata%' // trim(L1BDATA1%data_type) // &
      &  'while reading' // &
      & trim(QuantityName), MLSFile=L2AUXFile )
  end subroutine ReadL2AUXData_MF_hdf5

  ! ---------------------------------------------  WriteHDF5Data  -----

  subroutine WriteHDF5Data_2d( array, sd_id, returnStatus, sdName, &
    & already_there, start, sizes )

  use MLSFiles, only: Split_Path_Name
  use MLSHDF5, only: IsHDF5GroupPresent, MakeNestedGroups, SaveasHDF5ds

  ! Write an array to the file with l2FileHandle
  ! Optionally, write a bogus CounterMAF sd so the
  ! resulting file can masquerade as an l1BRad
  ! (Note that this bogus sd should only be written once for each file)
    real, dimension(:,:), intent(in)      :: array
    integer, intent(in)                     :: sd_id  ! From h5fopen or sfstart
    character (len=*), intent(in)           :: SDNAME ! may contain "/"
    integer, intent(out)                    :: returnStatus ! 0 unless error
    logical, intent(in), optional           :: already_there
    integer, dimension(:), intent(in), optional :: start
    integer, dimension(:), intent(in), optional :: sizes

    ! Local variables
    character(len=128)                      :: barename
    logical, parameter                      :: countEmpty = .true.
    character(len=128), dimension(25)       :: groupNames
    integer                                 :: grp_id ! Perhaps just root
    integer                                 :: n
    character(len=1024)                     :: path
    ! Executable
    grp_id = sd_id
    ! Is here a '/' in sdname?
    if ( index( sdname, '/' ) > 0 ) then
      ! Do the containing groups exist yet?
      call split_path_name ( sdname, path, bareName )
      if ( .not. IsHDF5GroupPresent( sd_id, trim(path)) ) then
        call List2Array ( path, groupnames, countEmpty, inseparator='/' )
        n = NumStringElements ( path, countEmpty, inseparator='/' )
        n = max( n, 2 )
        call MakeNestedGroups( grp_id, groupNames(1:n-1) )
      endif
    endif
    call SaveAsHDF5DS( grp_id, trim(sdName), &
      & array, start, sizes, may_add_to=.true., adding_to=already_there, &
      & fillValue=DEFAULTUNDEFINEDVALUE )
    returnStatus = 0

  end subroutine WriteHDF5Data_2d

  subroutine WriteHDF5Data_3d( array, sd_id, returnStatus, sdName, &
    & already_there, start, sizes )

  use MLSFiles, only: Split_Path_Name
  use MLSHDF5, only: IsHDF5GroupPresent, MakeNestedGroups, SaveasHDF5ds

  ! Write an array to the file with l2FileHandle
  ! Optionally, write a bogus CounterMAF sd so the
  ! resulting file can masquerade as an l1BRad
  ! (Note that this bogus sd should only be written once for each file)
    real, dimension(:,:,:), intent(in)      :: array
    integer, intent(in)                     :: sd_id  ! From h5fopen or sfstart
    character (len=*), intent(in)           :: SDNAME ! may contain "/"
    integer, intent(out)                    :: returnStatus ! 0 unless error
    logical, intent(in), optional           :: already_there
    integer, dimension(:), intent(in), optional :: start
    integer, dimension(:), intent(in), optional :: sizes

    ! Local variables
    character(len=128)                      :: barename
    logical, parameter                      :: countEmpty = .true.
    character(len=128), dimension(25)       :: groupNames
    integer                                 :: grp_id ! Perhaps just root
    integer                                 :: n
    character(len=1024)                     :: path
    ! Executable
    grp_id = sd_id
    ! Is here a '/' in sdname?
    if ( index( sdname, '/' ) > 0 ) then
      ! Do the containing groups exist yet?
      call split_path_name ( sdname, path, bareName )
      if ( .not. IsHDF5GroupPresent( sd_id, trim(path)) ) then
        call List2Array ( path, groupnames, countEmpty, inseparator='/' )
        n = NumStringElements ( path, countEmpty, inseparator='/' )
        n = max( n, 2 )
        call MakeNestedGroups( grp_id, groupNames(1:n-1) )
      endif
    endif
    call SaveAsHDF5DS( grp_id, trim(sdName), &
      & array, start, sizes, may_add_to=.true., adding_to=already_there, &
      & fillValue=DEFAULTUNDEFINEDVALUE )
    returnStatus = 0

  end subroutine WriteHDF5Data_3d

  ! ---------------------------------------------  WriteL2AUXData_FileHandle  -----

  subroutine WriteL2AUXData_FileHandle(l2aux, sd_id, returnStatus, sdName, &
    & NoMAFS, WriteCounterMAF, DimNames, Reuse_dimNames, hdfVersion)

  use MLSFiles, only: InitializeMLSFile

  ! Write l2aux to the file with l2FileHandle
  ! Optionally, write a bogus CounterMAF sd so the
  ! resulting file can masquerade as an l1BRad
  ! (Note that this bogus sd should only be written once for each file)
    type (L2AUXData_T), intent(inout) :: L2AUX
    integer, intent(in) :: sd_id                      ! From h5fopen or sfstart
    character (len=*), optional, intent(in) :: SDNAME ! Defaults to l2aux%name
    character (len=*), optional, intent(in) :: DimNames ! Comma-separated list
                                                        ! Otherwise automatic
                                                        ! (Requiring l2cf)
    integer, intent(in), optional :: NoMAFS
    logical, intent(in), optional :: WriteCounterMAF  ! Write bogus CounterMAF
    logical, intent(in), optional :: Reuse_dimNames   ! We already wrote them
    integer, intent(in), optional :: hdfVersion
    integer, intent(out) :: returnStatus           ! 0 unless error

    ! Local variables
    integer :: myhdfVersion
    type( MLSFile_T ) :: l2auxFile
    ! Executable code
    myhdfVersion = default_hdfversion_write
    if ( present(hdfVersion) ) myhdfVersion = hdfVersion
    returnStatus = InitializeMLSFile(l2auxFile, type=l_hdf, access=DFACC_RDWR, &
      & content='l2aux', name='unknown', hdfVersion=myhdfVersion)
    l2auxFile%FileID%f_id = sd_id
    l2auxFile%stillOpen = .true.
    call WriteL2AUXData(l2aux, l2auxFile, returnStatus, sdName, &
    & NoMAFS, WriteCounterMAF, DimNames, Reuse_dimNames)
  end subroutine WriteL2AUXData_FileHandle

  ! ---------------------------------------------  WriteL2AUXData_MLSFile  -----

  subroutine WriteL2AUXData_MLSFile(l2aux, L2AUXFile, returnStatus, sdName, &
    & NoMAFS, WriteCounterMAF, DimNames, Reuse_dimNames)

    use MLSFiles, only: HDFversion_4, HDFversion_5, &
      & MLS_CloseFile, MLS_OpenFile
    use Trace_M, only: Trace_Begin, Trace_End

  ! Write l2aux to the file with l2FileHandle
  ! Optionally, write a bogus CounterMAF sd so the
  ! resulting file can masquerade as an l1BRad
  ! (Note that this bogus sd should only be written once for each file)
    type (L2AUXData_T), intent(inout) :: L2AUX
    type(MLSFile_T)                :: L2AUXFile
    character (len=*), optional, intent(in) :: SDNAME ! Defaults to l2aux%name
    character (len=*), optional, intent(in) :: DimNames ! Comma-separated list
                                                        ! Otherwise automatic
                                                        ! (Requiring l2cf)
    integer, intent(in), optional :: NoMAFS
    logical, intent(in), optional :: WriteCounterMAF  ! Write bogus CounterMAF
    logical, intent(in), optional :: Reuse_dimNames   ! We already wrote them
    integer, intent(out) :: returnStatus           ! 0 unless error

    ! Local variables
    logical :: alreadyOpen
    integer :: Me = -1                  ! String index for trace cacheing

    ! Executable code
    call trace_begin ( me, 'WriteL2AUXData_MLSFile', cond=.false. )
    alreadyOpen = L2AUXFile%stillOpen
    if ( .not. alreadyOpen ) then
      call mls_openFile(L2AUXFile, returnStatus)
      if ( returnStatus /= 0 ) &
        call MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to open l2aux file for writing', MLSFile=L2AUXFile)
    end if
    if ( L2AUXFile%access == DFACC_RDONLY )  &
      & call MLSMessage(MLSMSG_Error, ModuleName, &
      & 'l2aux file is rdonly', MLSFile=L2AUXFile)
    select case (L2AUXFile%hdfVersion)
    case (HDFVERSION_4)
      call WriteL2AUXData_MF_hdf4(l2aux, L2AUXFile, returnStatus, sdName, &
    & NoMAFS, WriteCounterMAF, DimNames, Reuse_dimNames)
    case (HDFVERSION_5)
      call WriteL2AUXData_MF_hdf5( l2aux, L2AUXFile, returnStatus, sdName, &
    & NoMAFS, WriteCounterMAF )
    case default
      call MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unrecognized hdfVersion for l2aux file', MLSFile=L2AUXFile)
    end select
    if ( .not. alreadyOpen )  call mls_closeFile(L2AUXFile, returnStatus)
    L2AUXFile%errorCode = returnStatus
    L2AUXFile%lastOperation = 'write'
    call trace_end ( 'WriteL2AUXData_MLSFile', cond=.false. )
  end subroutine WriteL2AUXData_MLSFile

  ! ----------------------------------------  WriteL2AUXData_MF_hdf5  -----
  subroutine WriteL2AUXData_MF_hdf5(l2aux, L2AUXFile, returnStatus, sdName, &
    & NoMAFS, WriteCounterMAF )
  ! Write l2aux to the L2AUXFile
  ! Optionally, write a bogus CounterMAF sd so the
  ! resulting file can masquerade as an l1BRad
  ! (Note that this bogus sd should only be written once for each file;
  !  also note the attempt to convert l2aux%values to KIND of l1b radiances)
    use HDF5, only: H5gclose_F, H5gopen_F
    use MLS_DataProducts, only: DataProducts_T
    use MLSAUXData, only: Build_MLSAUXData
    use MLSHDF5, only: IsHDF5attributepresent, MakeHDF5Attribute
    use MLSL2Timings, only: ShowTimingNames
    use PCFHDR, only: H5_WriteGlobalAttr

    type (L2AUXData_T), intent(inout) :: L2AUX
    type(MLSFile_T)                :: L2AUXFile
    character (len=*), optional, intent(in) :: SDNAME ! Defaults to l2aux%name
    integer, intent(in), optional :: NoMAFS
    logical, intent(in), optional :: WriteCounterMAF  ! Write bogus CounterMAF
    integer, intent(out) :: returnStatus           ! 0 unless error

    ! Local variables
    integer :: grp_id
    integer :: myNoMAFS, MAF
    type(DataProducts_T) :: dataProduct
    logical :: myWriteCounterMAF
    integer, dimension(:), pointer :: CounterMAF ! bogus array
    logical, parameter             :: ALWAYSWRITEAS32BITS = .true.

    ! Executable code
    returnStatus = 0
    myWriteCounterMAF = .false.
    if ( present(WriteCounterMAF) ) myWriteCounterMAF = WriteCounterMAF
    myNoMAFS = 1
    if ( any(l2aux%dimensions%dimensionFamily == L_MAF) ) &
     & myNoMAFS = l2aux%dimensions(3)%noValues
    if ( present(NoMAFS) ) myNoMAFS = NoMAFS
    
    if ( .not. L2AUXFile%stillOpen ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'File not opened' , MLSFile=L2AUXFile )
    
    if ( .not. associated ( l2aux%values ) ) then
      call announce_error (0,&
        & "l2aux values not associated yet ", L2AUXFile=L2AUXFile )
      returnStatus = 1
    else
      if ( present(sdName) ) then
        dataProduct%name = sdName
      else
        call get_string ( l2aux%name, dataProduct%name, strip=.true. )
      end if
      if ( myWriteCounterMAF .or. ALWAYSWRITEAS32BITS ) then
        dataProduct%data_type = 'real'    ! same type as l1bradiances
      else
        dataProduct%data_type = 'double'  ! type of L2AUXData_T%values
      end if
      ! dims(1) = size(l2aux%values, 1)
      ! dims(2) = size(l2aux%values, 2)
      ! dims(3) = size(l2aux%values, 3)
      ! call Dump_L2AUX(l2AUX)
      if ( myWriteCounterMAF .or. ALWAYSWRITEAS32BITS ) then
        call WriteHDF5Data (  real(l2aux%values, r4), L2AUXFile%FileID%f_id, &
          & returnStatus, trim(dataProduct%name) )
        call WriteL2AUXAttributes(L2AUXFile%FileID%f_id, l2aux, trim(dataProduct%name))
      else
        call MLSMessage ( MLSMSG_Error, ModuleName // '%WriteL2AUXData_MF_hdf5', &
          & 'Unexpected execution path' , MLSFile=L2AUXFile )
      end if
      call h5_writeglobalattr(L2AUXFile%FileID%f_id, skip_if_already_there=.false.)
      ! Write phase and section names as file-level attributes?
      if ( PHASENAMEATTRIBUTES ) then
        call h5gopen_f(L2AUXFile%FileID%f_id, '/', grp_id, returnstatus)
        if ( .not. &
          & IsHDF5AttributePresent('/', L2AUXFile%FileID%f_id, 'Phase Names') ) &
          & call MakeHDF5Attribute(grp_id, &
          & 'Phase Names', trim(showTimingNames('phases', .true.)), .true.)
        if ( .not. &
          & IsHDF5AttributePresent('/', L2AUXFile%FileID%f_id, 'Section Names') ) &
          & call MakeHDF5Attribute(grp_id, &
          & 'Section Names', trim(showTimingNames('sections', .true.)), .true.)
        call h5gclose_f(grp_id, returnstatus)
      end if
      if ( .not. myWriteCounterMAF ) return
    
      ! Now create and write bogus counterMAF array
      if ( myNoMAFS < 1 ) then
        call announce_error(0, &
        & "Too few MAFs to fake CounterMAFs in l2aux file:  ", L2AUXFile=L2AUXFile )
        return
      end if
      nullify (CounterMAF)
      call allocate_test(CounterMAF,myNoMAFS,'counterMAF',ModuleName)
      ! dims(1) = myNoMAFS
      dataProduct%name = 'counterMAF'
      dataProduct%data_type = 'integer'
      do MAF=0, myNoMAFS-1
        counterMAF(MAF+1) = MAF
      end do
      call Build_MLSAuxData(L2AUXFile%FileID%f_id, dataProduct, counterMAF, &
      & myNoMAFS )
      call Deallocate_Test(CounterMAF,"CounterMAF",ModuleName)
    end if
  end subroutine WriteL2AUXData_MF_hdf5

  ! ----------------------------------------  WriteL2AUXData_MF_hdf4  -----
  subroutine WriteL2AUXData_MF_hdf4( l2aux, L2AUXFile, returnStatus, sdName, &
    & NoMAFS, WriteCounterMAF, DimNames, Reuse_dimNames )
  ! Write l2aux to the file with l2FileHandle
  ! Optionally, write a bogus CounterMAF sd so the
  ! resulting file can masquerade as an l1BRad
  ! (Note that this bogus sd should only be written once for each file)
    type (L2AUXData_T), intent(inout) :: L2AUX
    type(MLSFile_T)                :: L2AUXFile
    character (len=*), optional, intent(in) :: SDNAME ! Defaults to l2aux%name
    character (len=*), optional, intent(in) :: DimNames ! Comma-separated list
                                                        ! Otherwise automatic
                                                        ! (Requiring l2cf)
    integer, intent(in), optional :: NoMAFS
    logical, intent(in), optional :: WriteCounterMAF  ! Write bogus CounterMAF
    logical, intent(in), optional :: Reuse_dimNames   ! We already wrote them
    integer, intent(out) :: returnStatus           ! 0 unless error

    ! Local variables
    integer :: NODIMENSIONSUSED         ! No. real dimensions
    integer, dimension(:), pointer :: DIMSIZES ! Size of each dimension
    integer :: sdId                   ! sd id from sfcreate
    logical, dimension(L2AUXRank) :: GOODDIM ! Not L_None
    integer :: dimensionInData     ! Index
    integer :: dimensionInFile          ! Index
    integer :: dimID                    ! ID for HDF
    character (len=132) :: NAMESTRING   ! Name for sd
    character (len=132) :: DIMNAME      ! Name for dimension
    integer :: status                   ! Flag
    integer, parameter, dimension(L2AUXRank) :: stride = (/ 1, 1, 1/)
    integer, parameter, dimension(L2AUXRank) :: start = (/ 0, 0, 0/)
    logical :: myWriteCounterMAF
    logical :: myReuse_dimNames
    integer :: myNoMAFS, MAF
    integer :: myL2AUXRANK
    integer, dimension(:), pointer :: CounterMAF ! bogus array
    real(r8) :: HUGER4

    ! Executable code
    hugeR4 = real ( huge ( 0.0_r4), r8 )
    myWriteCounterMAF = .false.
    if ( present(WriteCounterMAF) ) myWriteCounterMAF = WriteCounterMAF
    myReuse_dimNames = .false.
    if ( present(Reuse_dimNames) ) myReuse_dimNames = Reuse_dimNames
    myNoMAFS = 1
    if ( any(l2aux%dimensions%dimensionFamily == L_MAF) ) &
     & myNoMAFS = l2aux%dimensions(3)%noValues
    if ( present(NoMAFS) ) myNoMAFS = NoMAFS
    myL2AUXRANK = L2AUXRANK
    if ( present(DimNames) ) myL2AUXRANK = NumStringElements(DimNames, countEmpty=.true.)
    nullify ( dimSizes )
    error = 0

    if (present(sdName)) then
      nameString=sdName
    else
      call get_string ( l2aux%name, nameString, strip=.true. )
    end if

    goodDim=l2aux%dimensions%dimensionFamily /= L_None
    noDimensionsUsed = min(COUNT(goodDim), myL2AUXRank)
    call allocate_test(dimSizes,noDimensionsUsed,'dimSizes',ModuleName)
    if ( present(DimNames) ) then
      dimSizes = PACK(l2aux%dimensions(1:noDimensionsUsed)%noValues, &
        & goodDim(1:noDimensionsUsed))
    else
      dimSizes = PACK(l2aux%dimensions%noValues, &
        & goodDim)
    end if

    ! Create the sd within the file
    sdId= SFcreate ( L2AUXFile%FileID%f_id, nameString, DFNT_FLOAT32, &
      & noDimensionsUsed, dimSizes)

    if ( .not. myReuse_dimNames ) then
    ! Now define the dimensions
    dimensionInFile=0
    do dimensionInData=1, myL2AUXRank
      if (l2aux%dimensions(dimensionInData)%dimensionFamily &
        &    /= L_None) then
        dimID=SFDIMID(sdId,dimensionInFile)
        ! Construct dimension name. For minor frame quantities MIF and MAF
        ! names are global, otherwise specific
        if ( present(DimNames) ) then
          call GetStringElement(trim(DimNames), dimName, dimensionInData, &
            & countEmpty=.TRUE.)
        else if ( (dimensionInData > 1) .and. (l2aux%minorFrame) ) then
          call GetModuleName( l2aux%instrumentModule, dimName)
          if (len_trim(dimName) < len(dimName)) dimName=TRIM(dimName)//'.'
          call get_string (lit_indices(l2aux%dimensions(dimensionInData)%dimensionFamily), &
            & dimName(LEN_TRIM(dimName)+1:))
        else
          call get_string (l2aux%name, dimName, strip=.true. )
          if (len_trim(dimName) < len(dimName)) dimName=TRIM(dimName)//'.'
          call get_string (lit_indices(l2aux%dimensions(dimensionInData)%dimensionFamily), &
            & dimName(LEN_TRIM(dimName)+1:))
        end if
        ! Write dimension name
        status=SFSDMName(dimID,TRIM(dimName))
        if ( status /= 0 ) then
          call output("dim name: ")
          call output(TRIM(dimName), advance='yes')
          call announce_error (  0, &
          & "Error setting dimension name to SDS l2aux file:", L2AUXFile=L2AUXFile)
        end if
        ! Write dimension scale
        status=SFSDScale(dimID, dimSizes(dimensionInFile+1), DFNT_FLOAT32, &
          & l2aux%dimensions(dimensionInData)%values)
        if ( status /= 0 ) then
          call output("dimID: ")
          call output(dimID, advance='yes')
          call output("dim name: ")
          call output(TRIM(dimName), advance='yes')
          call announce_error ( 0, &
          & "Error writing dimension scale in l2auxFile:", L2AUXFile=L2AUXFile )
        end if
        dimensionInFile=dimensionInFile+1
      end if
    end do
    end if

    ! Now write the data
    ! Make sure the data can be placed in a real
    status= SFWDATA_F90(sdId, start(1:noDimensionsUsed), &
      & stride(1:noDimensionsUsed), dimSizes, real ( &
      & max ( -hugeR4, min ( hugeR4, l2aux%values ) ) ) )
    if ( status /= 0 ) then
      call announce_error (0, &
      & "Error writing SDS data to  l2aux file:  ", L2AUXFile=L2AUXFile )
    end if

    call Deallocate_Test(dimSizes,"dimSizes",ModuleName)
    
    ! Terminate access to sd
    status = sfendacc(sdId)
    if ( status /= 0 ) then
      call announce_error (0,&
      & "Error ending access to the sd  ", L2AUXFile=L2AUXFile )
    end if
    returnStatus = error
    if ( .not. myWriteCounterMAF ) return
    
    ! Now create and write bogus counterMAF array
    if ( myNoMAFS < 1 ) then
      call announce_error(0, &
      & "Too few MAFs to fake CounterMAFs in l2aux file:  ", L2AUXFile=L2AUXFile )
      return
    end if
    nullify (CounterMAF, dimSizes)
    call allocate_test(CounterMAF,myNoMAFS,'counterMAF',ModuleName)
    call allocate_test(dimSizes,1,'dimSizes',ModuleName)
    dimSizes(1) = myNoMAFS
    sdId= SFcreate ( L2AUXFile%FileID%f_id, 'counterMAF', DFNT_INT32, &
      & 1, dimSizes)
    do MAF=0, myNoMAFS-1
      counterMAF(MAF+1) = MAF
    end do
    status= SFWDATA_F90(sdId, (/ 0 /), &
      & (/ 1 /) , dimSizes, CounterMAF)
    if ( status /= 0 ) then
      call announce_error (0,&
      & "Error writing counterMAF data to  l2aux file:  ", L2AUXFile=L2AUXFile )
    end if
    call Deallocate_Test(dimSizes,"dimSizes",ModuleName)
    call Deallocate_Test(CounterMAF,"CounterMAF",ModuleName)
    
    ! Terminate access to sd
    status = sfendacc(sdId)
    if ( status /= 0 ) then
      call announce_error (0,&
      & "Error ending access to the sd  ", L2AUXFile=L2AUXFile )
    end if
    returnStatus = error
    ! call Dump_L2AUX(l2AUX)

  end subroutine WriteL2AUXData_MF_hdf4

  ! ---------------------------------------  WriteL2AUXAttributes  -----
  subroutine WriteL2AUXAttributes ( L2FileHandle, l2aux, name)
    use MLSHDF5, only: MakeHDF5Attribute
    ! Writes the pertinent attributes for an l2aux
    ! Arguments
    integer, intent(in) :: L2FileHandle
    type (L2AUXData_T), intent(inout) :: L2AUX
    character(len=*) :: name
    ! Internal variables
    integer :: dim
    character(len=16), dimension(L2AUXRank) :: dim_name
    character(len=16), dimension(L2AUXRank) :: dim_unit
    character(len=16) :: dim_of_i
    character(len=16) :: framing
    character(len=2) :: i_char
    logical :: is_timing
    integer :: ndims
    character(len=*), parameter :: ottff = '1,2,3,4,5'
    ! Executable
    is_timing = ( index( lowercase(name), 'timing') > 0 )
    if ( is_timing ) then
      l2aux%majorframe = .false.
      l2aux%minorframe = .false.
      l2aux%minorframe = .false.
      l2aux%DIM_Names  = 'chunk,' // name(1:5) // ',none'
      l2aux%DIM_Units  =  'none,none,none'
      l2aux%VALUE_Units=  's'
    end if
    if ( DEEBUG ) then
      call output('Writing attributes to: ', advance='no')
      call output(trim(Name), advance='yes')
    end if
    call MakeHDF5Attribute(L2FileHandle, name, 'Title', name)
    call MakeHDF5Attribute(L2FileHandle, name, 'Units', &
      & trim(l2aux%VALUE_Units))
    call MakeHDF5Attribute(L2FileHandle, name, 'DimensionNames', &
      & trim(l2aux%DIM_Names))
    if ( l2aux%majorframe ) then
      framing = 'major'
    else if ( l2aux%minorframe ) then
      framing = 'minor'
    else
      framing = 'neither'
    end if
    call MakeHDF5Attribute(L2FileHandle, name, 'Framing', trim(framing))
    call MakeHDF5Attribute(L2FileHandle, name, 'InstrumentModule', &
      & l2aux%instrumentmodule)
    call MakeHDF5Attribute(L2FileHandle, name, 'QuantityType', &
      & l2aux%quantitytype)
    call MakeHDF5Attribute(L2FileHandle, name, 'MissingValue', &
      & (/ real(UNDEFINED_VALUE, r8) /) )
    dim_name = ' '
    dim_unit = ' '
    ndims = min( NumStringElements(trim(l2aux%DIM_Names), .true.), L2AUXRank )
    if ( ndims < 1 ) return
    call List2Array(trim(l2aux%DIM_Names), dim_name, .true.)
    call List2Array(trim(l2aux%DIM_Units), dim_unit, .true.)
    ! loop of dimensions
    do dim=1, ndims
      call GetStringElement (ottff, i_char, dim, .true.)
      dim_of_i = 'dim ' // trim(i_char)
      if ( trim(dim_unit(dim)) == ' ' ) dim_unit(dim) = 'none'
      call MakeHDF5Attribute(L2FileHandle, name, trim(dim_of_i), &
        & trim(dim_name(dim)))
      call MakeHDF5Attribute(L2FileHandle, name, trim(dim_of_i) // ' units', &
        & trim(dim_unit(dim)))
      if ( l2aux%dimensions(dim)%noValues > 0 ) then
        if ( AreDimValuesNonTrivial(l2aux%dimensions(dim)%DimensionFamily) ) then
          call MakeHDF5Attribute(L2FileHandle, name, trim(dim_of_i)// ' values', &
          & real(l2aux%dimensions(dim)%values))
        end if
      end if
    end do
  end subroutine WriteL2AUXAttributes

  ! -------------------------------------------------  GetDimSize  -----
  subroutine GetDimSize ( nameType, quantityTemplate, noMAFs, dim_size)

  ! Given a dim name type, e.g. l_MIF,
  ! returns corresponding dimension size, e.g. quantityTemplate%noSurfs

    ! Dummy arguments
    integer, intent(in)                :: noMAFs
    integer, intent(in)                :: nameType
    integer, intent(out)               :: dim_size
    type (QuantityTemplate_T), intent(in) :: quantityTemplate

    ! Executable code
    dim_size = 1
    select case (nameType)                                       
    case ( l_channel )  
      dim_size = quantityTemplate%noChans
    case ( l_frequency )  
      dim_size = quantityTemplate%noChans
    case ( l_geodAngle )  
      dim_size = quantityTemplate%noInstances
    case ( l_height )  
      dim_size = quantityTemplate%noSurfs
    case ( l_MAF )  
      dim_size = noMAFs
    case ( l_MIF )  
      dim_size = quantityTemplate%noSurfs
    case ( l_pressure )  
      dim_size = quantityTemplate%noSurfs
    case ( l_xyz )  
      dim_size = 3

    end select                                                       

  end subroutine GetDimSize

  ! ------------------------------------------------  GetDimStart  -----
  subroutine GetDimStart ( nameType, quantityTemplate, firstMaf, dim_start)

  ! Given a dim name type, e.g. l_MAF,
  ! returns corresponding dimension start, e.g. firstMaf

    ! Dummy arguments
    integer, intent(in)                :: firstMaf
    integer, intent(in)                :: nameType
    integer, intent(out)               :: dim_start
    type (QuantityTemplate_T), intent(in) :: quantityTemplate

    ! Executable code
    dim_start = 1
    select case (nameType)                                       
    case ( l_MAF )  
      dim_start = firstMaf

    end select                                                       

  end subroutine GetDimStart

  ! -----------------------------------------------  GetDimString  -----
  subroutine GetDimString ( nameType, dim_string)

  ! Given a dim name type, e.g. l_vmr,
  ! returns corresponding character string

    ! Dummy arguments
    integer, intent(in)                :: nameType
    character(len=*), intent(out)      :: dim_string

    ! Executable code
    dim_string = 'none'
    select case (nameType)                                       
    case ( l_channel )  
      dim_string = 'channel'
    case ( l_chunk )  
      dim_string = 'chunk'
    case ( l_frequency )  
      dim_string = 'frequency'
    case ( l_geodAngle )  
      dim_string = 'geodAngle'
    case ( l_height )  
      dim_string = 'height'
    case ( l_iteration )  
      dim_string = 'iteration'
    case ( l_MAF )  
      dim_string = 'MAF'
    case ( l_MIF )  
      dim_string = 'MIF'
    case ( l_pressure )  
      dim_string = 'pressure'
    case ( l_xyz )  
      dim_string = 'xyz'

    end select                                                       

  end subroutine GetDimString

  ! -------------------------------------  AreDimValuesNonTrivial  -----
  function AreDimValuesNonTrivial ( nameType )

  ! Given a dim name type, e.g. l_channel,
  ! tell whether the corresponding values are non-trivial
  ! Trivial means, e.g. {1 2 3 4} or {0 0 0 0}
    integer, intent(in)                   :: nameType
    logical                               :: AreDimValuesNonTrivial
    AreDimValuesNonTrivial = any( &
      & nameType == (/ &
      & l_frequency, &
      & L_IntermediateFrequency, l_USBFrequency, L_LSBFrequency, &
      & l_geodAngle, l_height, l_pressure &
      & /) &
      & )
  end function AreDimValuesNonTrivial

  ! -----------------------------------------------  GetDimValues  -----
  subroutine GetDimValues ( nameType, quantityTemplate, dim_values )

  ! Given a dim name type, e.g. l_channel,
  ! fills dimension with appropriate values, e.g. channels
  ! (Assumes dim_values already allocated)

    ! Dummy arguments
    integer, intent(in)                   :: nameType
    real(r8), dimension(:), intent(out)   :: dim_values
    type (QuantityTemplate_T), intent(in) :: quantityTemplate

    ! Internal variables
    integer :: i
    ! Executable code
    dim_values = 0._r8
    select case (nameType)                                       
    case ( l_channel, l_MAF, l_MIF, l_xyz )  
      do i=1, SIZE(dim_values)
       dim_values(i) = i
      end do
    case ( l_frequency, &
         & L_IntermediateFrequency, l_USBFrequency, L_LSBFrequency )  
      dim_values = quantityTemplate%frequencies
    case ( l_geodAngle )  
      dim_values = quantityTemplate%phi(1,:)   ! If not stacked, phi(surf, :)
    case ( l_height, l_pressure )  
      dim_values = quantityTemplate%surfs(:,1) ! If incoherent, surfs(:, inst)

    end select

  end subroutine GetDimValues

  ! --------------------------------------  GetQuantityAttributes  -----
  subroutine GetQuantityAttributes ( quantityType, &
   & framing, units_name, dim_names)

  ! Given a quantity type, e.g. l_vmr,
  ! returns major/minor/neither framing distinction, default unit name,
  ! and the 3 dimension names, e.g. (/l_channel, l_MIF, l_MAF/)

    use Init_Tables_Module, only: First_Lit, Last_Lit

    ! Dummy arguments
    integer, intent(in)                :: quantityType
    character(len=*), intent(out), optional      :: framing
    character(len=*), intent(out), optional      :: units_name
    integer, dimension(3), intent(out), optional :: dim_names

    type :: Attrib_t
      character(7) :: Framing = 'neither'
      character(10) :: Units_Name = 'NoUnits'
      integer :: dim_names(3) = (/ l_channel, l_MIF, l_MAF /)
    end type

    type(attrib_t), save :: Attrib(first_lit:last_lit)
    logical, save :: First = .true.

    if ( first ) then ! Can't do this with a DATA statement because Attrib_t has default initialization
      first = .false.
    ! Units should be consistent with InitQuantityTemplates in
    ! ConstructQuantitytypes. The ones that are commented out are in the list
    ! of members of the "quantity" type in init_tables_module.  If they become
    ! eligible to be written in L2AUX, they need to be in this list.  But
    ! check the fields first!
    !                                                  framing    units           dim_names
    ! Default for those not listed here is:           'neither', 'NoUnits   ', (/ l_channel, l_MIF, l_MAF /)
!     attrib(l_azimuth)                    = attrib_t('neither', 'deg       ', (/ l_none, l_none, l_none /) )
      attrib(l_baseline)                   = attrib_t('minor  ', 'K         ', (/ l_frequency, l_MIF, l_MAF /) )
!     attrib(l_boundaryPressure)           = attrib_t('major  ', 'hPa       ', (/ l_none, l_none, l_MAF /) )
      attrib(l_channel)                    = attrib_t('major  ', 'channel   ', (/ l_none, l_none, l_none /) )
      attrib(l_chisqchan)                  = attrib_t('major  ', 'channel   ', (/ l_channel, l_none, l_MAF /) )
      attrib(l_chisqmmaf)                  = attrib_t('major  ', 'NoUnits   ', (/ l_none, l_none, l_MAF /) )
      attrib(l_chisqmmif)                  = attrib_t('minor  ', 'NoUnits   ', (/ l_none, l_MIF, l_MAF /) )
      attrib(l_chunk)                      = attrib_t('major  ', 'chunk     ', (/ l_none, l_none, l_none /) )
      attrib(l_cloudExtinction)            = attrib_t('neither', 'NoUnits   ', (/ l_channel, l_none, l_MAF /) )
!     attrib(l_cloudIce)                   = attrib_t('minor  ', 'gm/m^3    ', (/ l_none, l_MIF, l_MAF /) )
      attrib(l_cloudInducedRadiance)       = attrib_t('minor  ', 'K         ', (/ l_channel, l_MIF, l_MAF /) )
!     attrib(l_cloudMinMax)                = attrib_t('minor  ', 'K         ', (/ l_channel, l_MIF, l_MAF /) )
      attrib(l_cloudRadSensitivity)        = attrib_t('minor  ', 'K         ', (/ l_channel, l_none, l_MAF /) )
!     attrib(l_cloudTemperature)           = attrib_t('minor  ', 'K         ', (/ l_none, l_MIF, l_MAF /) )
      attrib(l_cloudWater)                 = attrib_t('neither', 'NoUnits   ', (/ l_channel, l_none, l_MAF /) )
!     attrib(l_columnAbundance)            = attrib_t('major  ', 'mol/cm^2  ', (/ l_none, l_none, l_MAF /) )
      attrib(l_dnwt_ajn)                   = attrib_t('neither', 'NoUnits   ', (/ l_none, l_iteration, l_chunk /) )
      attrib(l_dnwt_axmax)                 = attrib_t('neither', 'NoUnits   ', (/ l_none, l_iteration, l_chunk /) )
      attrib(l_dnwt_cait)                  = attrib_t('neither', 'NoUnits   ', (/ l_none, l_iteration, l_chunk /) )
      attrib(l_dnwt_chiSqMinNorm)          = attrib_t('neither', 'NoUnits   ', (/ l_none, l_iteration, l_chunk /) )
      attrib(l_dnwt_chiSqNorm)             = attrib_t('neither', 'NoUnits   ', (/ l_none, l_iteration, l_chunk /) )
!     attrib(l_dnwt_chiSqRatio)            = attrib_t('neither', 'NoUnits   ', (/ l_none, l_iteration, l_chunk /) )
!     attrib(l_dnwt_count)                 = attrib_t('neither', 'NoUnits   ', (/ l_none, l_iteration, l_chunk /) )
      attrib(l_dnwt_diag)                  = attrib_t('neither', 'NoUnits   ', (/ l_none, l_iteration, l_chunk /) )
      attrib(l_dnwt_dxdx)                  = attrib_t('neither', 'NoUnits   ', (/ l_none, l_iteration, l_chunk /) )
      attrib(l_dnwt_dxdxl)                 = attrib_t('neither', 'NoUnits   ', (/ l_none, l_iteration, l_chunk /) )
      attrib(l_dnwt_dxn)                   = attrib_t('neither', 'NoUnits   ', (/ l_none, l_iteration, l_chunk /) )
      attrib(l_dnwt_dxnl)                  = attrib_t('neither', 'NoUnits   ', (/ l_none, l_iteration, l_chunk /) )
      attrib(l_dnwt_flag)                  = attrib_t('neither', 'NoUnits   ', (/ l_none, l_iteration, l_chunk /) )
      attrib(l_dnwt_fnmin)                 = attrib_t('neither', 'NoUnits   ', (/ l_none, l_iteration, l_chunk /) )
      attrib(l_dnwt_fnorm)                 = attrib_t('neither', 'NoUnits   ', (/ l_none, l_iteration, l_chunk /) )
      attrib(l_dnwt_gdx)                   = attrib_t('neither', 'NoUnits   ', (/ l_none, l_iteration, l_chunk /) )
      attrib(l_dnwt_gfac)                  = attrib_t('neither', 'NoUnits   ', (/ l_none, l_iteration, l_chunk /) )
      attrib(l_dnwt_gradn)                 = attrib_t('neither', 'NoUnits   ', (/ l_none, l_iteration, l_chunk /) )
      attrib(l_dnwt_sq)                    = attrib_t('neither', 'NoUnits   ', (/ l_none, l_iteration, l_chunk /) )
      attrib(l_dnwt_sqt)                   = attrib_t('neither', 'NoUnits   ', (/ l_none, l_iteration, l_chunk /) )
!     attrib(l_earthRadius)                = attrib_t('major  ', 'm         ', (/ l_none, l_none, l_MAF /) )
!     attrib(l_earthRefl)                  = attrib_t('neither', 'NoUnits   ', (/ l_none, l_none, l_none /) )
!     attrib(l_ecrToFOV)                   = attrib_t('minor  ', 'NoUnits   ', (/ l_none, l_MIF, l_MAF /) )
      attrib(l_effectiveOpticalDepth)      = attrib_t('minor  ', 'NoUnits   ', (/ l_channel, l_MIF, l_MAF /) )
      attrib(l_elevOffset)                 = attrib_t('neither', 'deg       ', (/ l_channel, l_MIF, l_MAF /) )
!     attrib(l_extinction)                 = attrib_t('neither', '1/km      ', (/ l_channel, l_none, l_MAF /) )
!     attrib(l_extinctionV2)               = attrib_t('neither', '1/km      ', (/ l_channel, l_none, l_MAF /) )
!     attrib(l_fieldAzimuth)               = attrib_t('minor  ', 'deg       ', (/ l_none, l_MIF, l_MAF /) )
!     attrib(l_fieldElevation)             = attrib_t('minor  ', 'deg       ', (/ l_none, l_MIF, l_MAF /) )
!     attrib(l_fieldStrength)              = attrib_t('minor  ', 'gauss     ', (/ l_none, l_MIF, l_MAF /) )
      attrib(l_frequency)                  = attrib_t('major  ', 'MHz       ', (/ l_none, l_none, l_none /) )
!     attrib(l_geolocation)                = attrib_t('major  ', 'NoUnits   ', (/ l_none, l_none, l_none /) )
!     attrib(l_gph)                        = attrib_t('major  ', 'm         ', (/ l_none, l_none, l_MAF /) )
      attrib(l_heightOffset)               = attrib_t('minor  ', 'm         ', (/ l_channel, l_MIF, l_MAF /) )
!     attrib(l_isotopeRatio)               = attrib_t('neither', 'NoUnits   ', (/ l_none, l_none, l_none /) )
!     attrib(l_iwc)                        = attrib_t('minor  ', 'gm/m^3    ', (/ l_none, l_MIF, l_MAF /) )
      attrib(l_jacobian_cols)              = attrib_t('neither', 'NoUnits   ', (/ l_none, l_iteration, l_chunk /) )
      attrib(l_jacobian_rows)              = attrib_t('neither', 'NoUnits   ', (/ l_none, l_iteration, l_chunk /) )
!     attrib(l_l1bMAFBaseline)             = attrib_t('major  ', 'K         ', (/ l_none, l_none, l_MAF /) )
!     attrib(l_l1bMIF_TAI)                 = attrib_t('minor  ', 's         ', (/ l_none, l_MIF, l_MAF /) )
      attrib(l_limbSidebandFraction)       = attrib_t('neither', 'NoUnits   ', (/ l_channel, l_none, l_none /) )
!     attrib(l_lineCenter)                 = attrib_t('neither', 'MHz       ', (/ l_none, l_none, l_none /) )
!     attrib(l_lineWidth)                  = attrib_t('neither', 'MHz       ', (/ l_none, l_none, l_none /) )
!     attrib(l_lineWidth_TDep)             = attrib_t('neither', 'MHz       ', (/ l_none, l_none, l_none /) )
      attrib(l_losTransFunc)               = attrib_t('neither', 'NoUnits   ', (/ l_frequency, l_MIF, l_MAF /) )
      attrib(l_losVel)                     = attrib_t('minor  ', 'm/s       ', (/ l_xyz, l_MIF, l_MAF /) )
!     attrib(l_lowestRetrievedPressure)    = attrib_t('neither', 'log10(hPa)', (/ l_none, l_none, l_none /) )
!     attrib(l_magneticField)              = attrib_t('minor  ', 'gauss     ', (/ l_xyz, l_MIF, l_MAF /) )
      attrib(l_massMeanDiameterIce)        = attrib_t('minor  ', 'NoUnits   ', (/ l_none, l_MIF, l_MAF /) )
      attrib(l_massMeanDiameterWater)      = attrib_t('minor  ', 'NoUnits   ', (/ l_none, l_MIF, l_MAF /) )
!     attrib(l_mifDeadTime)                = attrib_t('neither', 's         ', (/ l_none, l_none, l_none /) )
      attrib(l_mifExtinction)              = attrib_t('minor  ', '1/km      ', (/ l_channel, l_MIF, l_MAF /) )
!     attrib(l_mifExtinctionExtrapolation) = attrib_t('neither', 'NoUnits   ', (/ l_none, l_none, l_none /) )
!     attrib(l_mifExtinctionForm)          = attrib_t('neither', 'NoUnits   ', (/ l_none, l_none, l_none /) )
      attrib(l_mifExtinctionV2)            = attrib_t('minor  ', '1/km      ', (/ l_channel, l_MIF, l_MAF /) )
      attrib(l_noiseBandwidth)             = attrib_t('neither', 'MHz       ', (/ l_channel, l_none, l_none /) )
!     attrib(l_noRadsBinned)               = attrib_t('minor  ', 'NoUnits   ', (/ l_channel, l_MIF, l_MAF /) )
      attrib(l_noRadsPerMIF)               = attrib_t('minor  ', 'NoUnits   ', (/ l_none, l_MIF, l_MAF /) )
!     attrib(l_numGrad)                    = attrib_t('neither', 'iteration ', (/ l_none, l_iteration, l_chunk /) )
      attrib(l_numJ)                       = attrib_t('neither', 'iteration ', (/ l_none, l_iteration, l_chunk /) )
      attrib(l_opticalDepth)               = attrib_t('minor  ', 'NoUnits   ', (/ l_channel, l_MIF, l_MAF /) )
      attrib(l_orbitInclination)           = attrib_t('minor  ', 'deg       ', (/ l_none, l_none, l_none /) )
      attrib(l_AscDescMode)                  = attrib_t('neither', 'NoUnits   ', (/ l_none, l_none, l_none /) )
      attrib(l_phaseTiming)                = attrib_t('minor  ', 's         ', (/ l_none, l_MIF, l_MAF /) )
      attrib(l_phiTan)                     = attrib_t('minor  ', 'deg       ', (/ l_none, l_MIF, l_MAF /) )
      attrib(l_ptan)                       = attrib_t('minor  ', 'log10(hPa)', (/ l_none, l_MIF, l_MAF /) )
!     attrib(l_quality)                    = attrib_t('major  ', 'NoUnits   ', (/ l_none, l_MIF, l_MAF /) )
      attrib(l_radiance)                   = attrib_t('minor  ', 'K         ', (/ l_channel, l_MIF, l_MAF /) )
!     attrib(l_refGPH)                     = attrib_t('major  ', 'm         ', (/ l_none, l_MIF, l_MAF /) )
      attrib(l_reflTemp)                   = attrib_t('major  ', 'K         ', (/ l_none, l_none, l_MAF /) )
      attrib(l_reflSpill)                  = attrib_t('major  ', 'NoUnits   ', (/ l_channel, l_none, l_MAF /) )
!     attrib(l_RHI)                        = attrib_t('minor  ', 'NoUnits   ', (/ l_none, l_MIF, l_MAF /) )
      attrib(l_scanResidual)               = attrib_t('minor  ', 'm         ', (/ l_none, l_MIF, l_MAF /) )
!     attrib(l_scatteringAngle)            = attrib_t('minor  ', 'deg       ', (/ l_none, l_MIF, l_MAF /) )
      attrib(l_scECI)                      = attrib_t('minor  ', 'm         ', (/ l_xyz, l_MIF, l_MAF /) )
      attrib(l_scGeocAlt)                  = attrib_t('minor  ', 'm         ', (/ l_xyz, l_MIF, l_MAF /) )
      attrib(l_scVelECI)                   = attrib_t('minor  ', 'm/s       ', (/ l_xyz, l_MIF, l_MAF /) )
      attrib(l_scVelECR)                   = attrib_t('minor  ', 'm/s       ', (/ l_xyz, l_MIF, l_MAF /) )
      attrib(l_singleChannelRadiance)      = attrib_t('minor  ', 'K         ', (/ l_none, l_MIF, l_MAF /) )
      attrib(l_sizeDistribution)           = attrib_t('minor  ', 'NoUnits   ', (/ l_none, l_MIF, l_MAF /) )
      attrib(l_spaceRadiance)              = attrib_t('neither', 'K         ', (/ l_none, l_none, l_none /) )
!     attrib(l_status)                     = attrib_t('major  ', 'NoUnits   ', (/ l_none, l_none, l_MAF /) )
      attrib(l_strayRadiance)              = attrib_t('major  ', 'K         ', (/ l_channel, l_none, l_MAF /) )
!     attrib(l_surfaceHeight)              = attrib_t('major  ', 'm         ', (/ l_channel, l_none, l_MAF /) )
      attrib(l_surfacetype)                = attrib_t('neither', 'NoUnits   ', (/ l_none, l_none, l_none /) )
      attrib(l_systemTemperature)          = attrib_t('neither', 'K         ', (/ l_channel, l_none, l_none /) )
!     attrib(l_totalPowerWeight)           = attrib_t('neither', 'NoUnits   ', (/ l_channel, l_none, l_MAF /) )
      attrib(l_tngtECI)                    = attrib_t('minor  ', 'm         ', (/ l_xyz, l_MIF, l_MAF /) )
      attrib(l_tngtGeocAlt)                = attrib_t('minor  ', 'm         ', (/ l_xyz, l_MIF, l_MAF /) )
      attrib(l_tngtGeodAlt)                = attrib_t('minor  ', 'm         ', (/ l_xyz, l_MIF, l_MAF /) )
!     attrib(l_TScat)                      = attrib_t('major  ', 'K         ', (/ l_channel, l_none, l_MAF /) )
      attrib(l_vmr)                        = attrib_t('neither', 'vmr       ', (/ l_channel, l_none, l_MAF /) )
      ! The following are not quantity types listed under t_quantityType in
      ! init_tables_module, or in InitQuantityTemplates in
      ! ConstructQuantityTemplates.  At least some of them are dimension
      ! types, for which a second call to this subroutine fetches the units
      ! attribute of the dimension (instead of the vector quantity).
      ! Others might be fossils.
      attrib(l_geodAngle)                  = attrib_t('neither', 'deg       ', (/ l_none, l_none, l_none /) )
      attrib(l_iteration)                  = attrib_t('major  ', 'iteration ', (/ l_none, l_none, l_none /) )
      attrib(l_MAF)                        = attrib_t('major  ', 'MAF       ', (/ l_none, l_none, l_none /) )
      attrib(l_MIF)                        = attrib_t('major  ', 'MIF       ', (/ l_none, l_none, l_none /) )
      attrib(l_totalExtinction)            = attrib_t('neither', 'NoUnits   ', (/ l_channel, l_none, l_MAF /) )
      attrib(l_xyz)                        = attrib_t('major  ', 'xyz       ', (/ l_none, l_none, l_none /) )
    end if

    ! Executable code
    if ( present(framing) )    framing    = attrib(quantityType)%framing
    if ( present(units_name) ) units_name = attrib(quantityType)%units_name
    if ( present(dim_names) )  dim_names  = attrib(quantityType)%dim_names

  end subroutine GetQuantityAttributes

  ! ------------------------------------  GetQuantityTypeFromName  -----
  function GetQuantityTypeFromName (name)  result(quantityType)

  ! Given quantity name, e.g. '/R4:640.B29M:HOCL.S0.MB11-3 chisqMMIF CorePlusR4'
  ! return a quantity type, e.g. l_chisqmmif

    ! Dummy arguments
    integer                :: quantityType
    character(len=*), intent(in)       :: name
    ! Local variables
    character(len=len(name)) :: myName
    ! Executable code
    quantityType = -999
    if ( name == '' ) return
    myName = lowerCase(name)
    if ( index(trim(myName), 'chisq') > 0 ) then
      if ( index(trim(myName), 'mmaf') > 0 ) then
        quantityType = l_chisqmmaf
      else if ( index(trim(myName), 'mmif') > 0 ) then
        quantityType = l_chisqmmif
      else if ( index(trim(myName), 'chan') > 0 ) then
        quantityType = l_chisqchan
      else
        ! This is binned chi^2
        ! unlike the others, it is neither major nor minor frame
        ! we'll say it is most like dnwt_chisqnorm
        quantityType = l_dnwt_chiSqNorm
      end if
    else if ( index(trim(myName), 'noradspermif') > 0 ) then
      ! Just like chisqmmif
      quantityType = l_chisqmmif
    else if ( index(trim(myName), 'pcf') > 0 ) then
      quantityType = RECOGNIZEDBUTNOTCOPIEDQT ! was -999
    else if ( index(trim(myName), 'coremetadata') > 0 ) then
      quantityType = RECOGNIZEDBUTNOTCOPIEDQT ! was -999
    else
      quantityType = l_radiance
    end if                      
                      
  end function GetQuantityTypeFromName

  ! ---------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( WHERE, full_message, CODE, L2AUXFile )
    use MLSFiles, only: Dump
    use Tree, only: Where_At => Where
    integer, intent(in) :: WHERE   ! Tree node where error was noticed
    character(len=*), intent(in)    :: full_message
    integer, intent(in), optional :: CODE    ! Code for error message
    type(MLSFile_T), optional :: L2AUXFile

    error = max(error,1)
    call output ( '***** At ' )
    if ( where > 0 ) then
      call print_source ( where_at(where) )
    else
      call output ( '(no lcf node available)' )
    end if
    call output ( ' L2AUXData complained: ' )


    call output ( " Caused the following error: ", advance='yes', &
      & from_where=ModuleName )
    call output ( trim(full_message), advance='yes', &
      & from_where=ModuleName )
    if ( present(code) ) then
      select case ( code )
      end select
    end if
    if ( present(L2AUXFile) ) call dump(L2AUXFile)
  end subroutine ANNOUNCE_ERROR

!=============================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module L2AUXData
!=============================================================================


! $Log$
! Revision 2.105  2019/01/31 19:50:55  pwagner
! Removed unused stuff
!
! Revision 2.104  2018/04/19 01:14:16  vsnyder
! Remove USE statements for unused names
!
! Revision 2.103  2017/10/20 20:05:07  pwagner
! Reduce default printing
!
! Revision 2.102  2017/10/11 23:55:55  pwagner
! Take pains to read int- and char-valued datasets, too
!
! Revision 2.101  2017/08/10 22:46:13  pwagner
! Add WriteHDF5Data
!
! Revision 2.100  2016/08/26 00:17:52  pwagner
! Removed two unused dummy args
!
! Revision 2.99  2016/07/28 01:45:07  vsnyder
! Refactor dump and diff
!
! Revision 2.98  2016/05/19 23:19:02  pwagner
! Added api for cpL2AUXData
!
! Revision 2.97  2015/08/03 21:43:03  pwagner
! Made quantityType optional in call to ReadL2AUXData
!
! Revision 2.96  2015/03/28 02:47:33  vsnyder
! Added stuff to trace allocate/deallocate addresses
!
! Revision 2.95  2014/09/05 01:03:35  vsnyder
! More complete and accurate allocate/deallocate size tracking
!
! Revision 2.94  2014/09/05 00:49:06  vsnyder
! EmpiricalGeometry.f90 -- Wrong comment
!
! Revision 2.93  2014/07/18 23:17:11  pwagner
! Aimed for consistency in names passed to allocate_test
!
! Revision 2.92  2014/04/07 18:06:58  pwagner
! May specify AscDescMode when DirectWrite-ing swaths
!
! Revision 2.91  2014/03/31 23:43:29  pwagner
! Commented-out unused stuff; renamed procedure ResizeL2AUXData, generalizing it to expand or contract
!
! Revision 2.90  2013/09/24 23:47:22  vsnyder
! Use Where instead of Source_Ref for messages
!
! Revision 2.89  2013/08/31 02:29:12  vsnyder
! Replace MLSMessageCalls with trace_begin and trace_end
!
! Revision 2.88  2013/07/19 01:24:29  vsnyder
! Sort some stuff, revise GetQuantityAttributes a bit more
!
! Revision 2.87  2013/07/18 01:10:57  vsnyder
! Remove scVel since it's ambiguous whether it's ECI or ECR, and nobody
! uses it anyway.

! Revision 2.86  2012/01/25 01:16:41  pwagner
! Improved error msg; snipped commented-out lines

! Revision 2.85  2011/07/07 00:39:15  pwagner
! Accepts options as arg for dumps

! Revision 2.84  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away

! Revision 2.83  2007/10/24 00:15:53  pwagner
! Removed unused declarations

! Revision 2.82  2007/08/13 17:39:42  pwagner
! Push some procedures onto new MLSCallStack

! Revision 2.81  2007/06/21 00:54:08  vsnyder
! Remove tabs, which are not part of the Fortran standard

! Revision 2.80  2006/05/19 22:49:15  pwagner
! May rename copied SDs

! Revision 2.79  2006/01/26 00:34:50  pwagner
! demoted more use statements from module level to speed Lahey compiles

! Revision 2.78  2005/12/21 18:45:29  pwagner
! Should recognize but not copy coremetadata, pcf

! Revision 2.77  2005/12/14 01:45:21  pwagner
! Attribute values to phase, section timing more reasonable

! Revision 2.76  2005/10/11 17:39:58  pwagner
! Added MLSFile interface to cpL2AUXData

! Revision 2.75  2005/09/21 23:17:34  pwagner
! Unnecessary changes

! Revision 2.74  2005/08/25 20:21:41  pwagner
! Ensure returnStatus defined in ReadL2AUXData_MLSFile

! Revision 2.73  2005/08/19 23:27:02  pwagner
! Uses '*' as wildcard sdList string in cpL2AUXData

! Revision 2.72  2005/08/05 20:38:31  pwagner
! L2AUXFile arg to ReadL2AUXFile now a pointer

! Revision 2.71  2005/07/06 00:29:26  pwagner
! optional arg options determines whether cpL2AUXData dumps DS names

! Revision 2.70  2005/06/22 18:57:01  pwagner
! Reworded Copyright statement, moved rcs id

! Revision 2.69  2005/06/14 20:41:02  pwagner
! Interfaces changed to accept MLSFile_T args

! Revision 2.68  2005/03/03 02:10:51  vsnyder
! Remove unused symbols, spiff up some dumps

! Revision 2.67  2004/08/19 00:19:28  pwagner
! Tells ReadL1BData to skip warnings about missing counterMAFs

! Revision 2.66  2004/08/17 17:09:45  pwagner
! L2AUX files shouldn't need padding when read as l1b

! Revision 2.65  2004/08/16 23:43:22  livesey
! Added ability to output minor frame baselines

! Revision 2.64  2004/08/04 23:19:57  pwagner
! Much moved from MLSStrings to MLSStringLists

! Revision 2.63  2004/06/29 18:05:26  pwagner
! May write phase, section names as file-level attributes

! Revision 2.62  2004/04/16 00:48:13  livesey
! Added singleChannelRadiance output

! Revision 2.61  2004/03/08 22:33:29  pwagner
! Bypass reading QuantityType attribute (why always 0)

! Revision 2.60  2004/02/26 22:05:06  pwagner
! Can copy l2aux file w/o knowing ds names; acts more gracefully if no attributes

! Revision 2.59  2004/02/05 23:36:41  pwagner
! WriteL2AUXAttributes now public

! Revision 2.58  2004/01/27 21:38:09  pwagner
! Added cpL2AUXData

! Revision 2.57  2003/09/03 05:25:49  livesey
! Bug fix in hdf5 readl2auxdata.

! Revision 2.56  2003/07/15 23:39:47  pwagner
! Disabled most printing

! Revision 2.55  2003/05/30 00:10:02  livesey
! Bug fix with reflTemp

! Revision 2.54  2003/05/30 00:08:54  livesey
! Added antenna loss terms

! Revision 2.53  2003/05/29 16:43:02  livesey
! Renamed sideband fraction

! Revision 2.52  2003/05/12 02:06:32  livesey
! Bound r8->r4 conversion

! Revision 2.51  2003/04/25 19:55:09  livesey
! Added more useful error message

! Revision 2.50  2003/03/07 00:42:13  pwagner
! Abbreviated Units names; removed spaces from attribute names

! Revision 2.49  2003/02/21 23:42:21  pwagner
! Also writes Fill Value attribute

! Revision 2.48  2003/02/12 21:52:34  pwagner
! Renames blank dim units to none

! Revision 2.47  2003/02/07 21:44:56  pwagner
! Capitalized 1st letter of each attribute name

! Revision 2.46  2003/01/30 01:02:28  pwagner
! Writing attributes for hdf5 files; global and data set

! Revision 2.45  2003/01/18 02:37:03  livesey
! Made the readl2aux data stuff work from the l2cf by adding the
! quantityType argument.

! Revision 2.44  2003/01/17 23:11:26  pwagner
! Moved most ops out of LoinL2AUXData to SetupL2AUXData

! Revision 2.43  2003/01/14 00:41:43  pwagner
! Added GetQuantityAttributes and getDimString; new fields in L2AUXData_T

! Revision 2.42  2002/12/10 00:41:28  pwagner
! In principle can now read hdf5-formatted l2aux files; untested; njl has other plans

! Revision 2.41  2002/12/07 00:25:42  pwagner
! Using SaveAsHDF5DS to write l2aux%values; it works

! Revision 2.40  2002/12/06 01:06:13  pwagner
! Finally writes radaiance-like l2aux as hdf5 files

! Revision 2.39  2002/12/05 19:46:23  pwagner
! Changes to speed up compiling tree-walker

! Revision 2.38  2002/12/03 18:04:02  pwagner
! Repaired bug that caused WriteL2AUXData files to be tiny

! Revision 2.37  2002/12/02 23:42:12  pwagner
! Optional param checkDimNames to ReadL2AUXData; defaults to FALSE

! Revision 2.36  2002/12/02 19:11:13  pwagner
! Corrected data types of counterMAF and dimensions

! Revision 2.35  2002/11/29 22:46:28  livesey
! Various bug fixes / improvements.

! Revision 2.34  2002/11/29 18:50:07  livesey
! Initialised a variable

! Revision 2.33  2002/11/26 22:16:41  jonathan
! Comment-out dump_l2aux diagnostics

! Revision 2.32  2002/11/25 18:04:52  pwagner
! Consistent with latest changes to MLSAuxData

! Revision 2.31  2002/11/22 21:48:02  pwagner
! Fleshed out WriteL2AUXData_hdf5; untested yet

! Revision 2.30  2002/11/13 01:09:47  pwagner
! Beginnings of attempt to write hdf5 L2AUX; incomplete

! Revision 2.29  2002/11/08 23:14:41  pwagner
! Should work again with mlsl2

! Revision 2.28  2002/11/08 18:25:33  jonathan
! Changes to allow writing rank 2; also reuse_DimNames

! Revision 2.27  2002/11/06 02:01:06  livesey
! Changes to fill from l2aux

! Revision 2.26  2002/11/06 00:18:37  pwagner
! Can WriteL2AUXData w/o l2cf: useable by small utility programs

! Revision 2.25  2002/10/08 17:36:21  pwagner
! Added idents to survive zealous Lahey optimizer

! Revision 2.24  2002/08/21 01:04:53  livesey
! Changed to single precision for data

! Revision 2.23  2002/08/15 21:47:04  pwagner
! WriteL2AuxData now returns non-zero status if it fails

! Revision 2.22  2001/11/01 21:03:59  pwagner
! Uses new sfwdata_f90 generic; added toc

! Revision 2.21  2001/10/26 23:13:18  pwagner
! Provides a single dump module interface and details

! Revision 2.20  2001/10/08 23:41:27  pwagner
! Improved dump routines

! Revision 2.19  2001/10/05 23:32:27  pwagner
! Added majorframe to data type; trimmed unused stuff

! Revision 2.18  2001/08/06 18:35:24  pwagner
! Added dump_l2aux

! Revision 2.17  2001/07/11 20:50:46  dwu
! fix problem in readl2auxdata

! Revision 2.16  2001/05/30 23:53:31  livesey
! Changed for new version of L1BData

! Revision 2.15  2001/05/12 00:18:40  livesey
! Tidied up issues with array bounds etc.

! Revision 2.14  2001/05/03 20:32:19  vsnyder
! Cosmetic changes

! Revision 2.13  2001/05/02 22:24:20  pwagner
! Removed SDPToolkit use

! Revision 2.12  2001/04/26 02:44:17  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic

! Revision 2.11  2001/04/12 22:19:33  vsnyder
! Improved an error message

! Revision 2.10  2001/04/10 22:27:47  vsnyder
! Nullify explicitly instead of with <initialization> so as not to give
! pointers the SAVE attribute.  <initialization> is NOT executed on each
! entry to a procedure.

! Revision 2.9  2001/04/07 00:14:27  pwagner
! Added announce_error

! Revision 2.8  2001/03/15 18:42:29  livesey
! Removed quotes from dimension name prefixes

! Revision 2.7  2001/03/08 02:20:12  livesey
! Added strip argument to a call to get_string

! Revision 2.6  2001/03/06 22:40:47  livesey
! Working version

! Revision 2.5  2001/02/14 23:41:33  livesey
! Removed irrelevant numProfs argument

! Revision 2.4  2001/01/03 00:46:19  pwagner
! Changed sfgetinfo to sfginfo

! Revision 2.3  2000/12/04 23:34:38  vsnyder
! Move more of addItemToDatabase into the include.

! Revision 2.2  2000/12/04 21:48:29  pwagner
! ReadL2AUXData completed

! Revision 2.1  2000/12/02 01:12:00  pwagner
! Added ReadL2AUXData

! Revision 2.0  2000/09/05 18:57:02  ahanzel
! Changing file revision to 2.0.

! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry


