! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module L2AUXData                 ! Data types for storing L2AUX data internally

  use Allocate_Deallocate, only: Allocate_test, Deallocate_test
  use Dump_0, only: DUMP
  use Hdf, only: &
    & DFNT_FLOAT32, DFNT_INT32, &
    & SFCREATE, SFDIMID, SFSDSCALE, &
    & SFENDACC, SFRDATA_F90, SFN2INDEX, SFSELECT, SFGINFO, &
    & SFGDINFO, SFSDMNAME, SFWDATA_F90
  use INIT_TABLES_MODULE, only: &
    L_BASELINE, L_BOUNDARYPRESSURE, &
    L_CHANNEL, L_CHISQCHAN, L_CHISQMMAF, L_CHISQMMIF, L_CHUNK, L_CLOUDICE, &
    L_CLOUDEXTINCTION, L_CLOUDWATER, &
    L_CLOUDINDUCEDRADIANCE, L_CLOUDRADSENSITIVITY, &
    L_COLUMNABUNDANCE, L_DNWT_AJN, L_DNWT_AXMAX, &
    L_DNWT_CAIT, L_DNWT_CHISQMINNORM, L_DNWT_CHISQNORM, L_DNWT_DIAG, &
    L_DNWT_DXDX, L_DNWT_DXDXL, L_DNWT_DXN, L_DNWT_DXNL, L_DNWT_FLAG, &
    L_DNWT_FNMIN, L_DNWT_FNORM, L_DNWT_GDX, L_DNWT_GFAC, L_DNWT_GRADN, &
    L_DNWT_SQ, L_DNWT_SQT, &
    L_EARTHREFL, L_EARTHRADIUS, L_EFFECTIVEOPTICALDEPTH, L_ELEVOFFSET, &
    L_EXTINCTION, L_FREQUENCY, L_GEODALTITUDE, L_GEODANGLE, &
    L_HEIGHT, L_HEIGHTOFFSET, L_INTERMEDIATEFREQUENCY, &
    L_ITERATION, L_JACOBIAN_COLS, L_JACOBIAN_ROWS, &
    L_LOSTRANSFUNC, L_LOSVEL, L_LSBFREQUENCY, L_MAGNETICFIELD, &
    L_MAF, L_MASSMEANDIAMETERICE, L_MASSMEANDIAMETERWATER, L_MIF, &
    L_NOISEBANDWIDTH, L_NONE, L_NUMJ, L_ORBITINCLINATION, L_OPTICALDEPTH, &
    L_PHITAN, L_PRESSURE, L_PTAN, L_RADIANCE, &
    L_SCANRESIDUAL, L_SCECI, L_SCGEOCALT, L_SCVEL, &
    L_SCVELECI, L_SCVELECR, L_LIMBSIDEBANDFRACTION, L_SIZEDISTRIBUTION, &
    L_SPACERADIANCE, L_SURFACETYPE, L_SYSTEMTEMPERATURE, &
    L_TNGTECI, L_TNGTGEOCALT, L_TNGTGEODALT, &
    L_TOTALEXTINCTION, L_USBFREQUENCY, L_VMR, L_XYZ
  use intrinsic, only: LIT_INDICES !, L_CHANNEL, &
!    & L_MAF, L_MIF, L_NONE
  use L1BData, only: L1BDATA_T, READL1BDATA
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: R8, R4
  use MLSL2Options, only: DEFAULT_HDFVERSION_READ, DEFAULT_HDFVERSION_WRITE
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ALLOCATE, MLSMSG_DEALLOCATE, &
    & MLSMSG_ERROR, MLSMSG_WARNING
  use MLSSignals_m, only: GETMODULENAME, MODULES
  use MLSStrings, only: Array2List, GetStringElement, List2Array, &
    & NumStringElements
  use Output_M, only: OUTPUT
  use QuantityTemplates, only: QuantityTemplate_T
  use STRING_TABLE, only: GET_STRING, DISPLAY_STRING
  use Tree, only: SOURCE_REF

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
! DestroyL2AUXContents            Deallocates all the arrays for one l2aux
! DestroyL2AUXDatabase            Deallocates all the arrays for entire database
! Dump                            Prints info on one quantity or entire database
! ExpandL2AUXDataInPlace          Expands an l2aux quantity to take more profiles
! ReadL2AUXData                   Reads an l2aux quantity from a file
! SetupNewL2AUXRecord             Allocates the arrays for an l2aux quantity
! WriteL2AUXData                  Writes an l2aux quantity to a file
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
! DestroyL2AUXContents ( L2AUXData_T l2aux )
! ExpandL2AUXDataInPlace ( L2AUXData_T l2aux, int newSize )
! int AddL2AUXToDatabase ( *L2AUXData_T DATABASE(:), L2AUXData_T ITEM )
! DestroyL2AUXDatabase ( *L2AUXData_T DATABASE(:) )
! Dump ( l2auxData_T L2aux(:), [char* Name], [int Details] )
!    or Dump ( l2auxData_T L2aux, [int Details] )
! ReadL2AUXData (int sd_id, char* quantityname, l2auxData_T l2aux, 
!    [int firstProf], [int lastProf])
! WriteL2AUXData(l2auxData_T l2aux, int l2FileHandle, int returnStatus, 
!    [char* sdName], [int NoMAFS], [log WriteCounterMAF], [char* DimNames])
! === (end of api) ===

  private
  public :: L2AUX_Dimension_T, L2AUXData_T, L2AUXRANK
  public :: AddL2AUXToDatabase, DestroyL2AUXDatabase, Dump
  public :: SetupNewL2AUXRecord, DestroyL2AUXContents, ExpandL2AUXDataInPlace
  public :: ReadL2AUXData, WriteL2AUXData
! public :: GetDimString, GetQuantityAttributes

  interface DUMP
    module procedure Dump_L2AUX
    module procedure Dump_L2AUX_Database
  end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
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

  real, parameter    :: UNDEFINED_VALUE = -999.99 ! Same as %template%badvalue
  integer, parameter :: L2AUXRANK=3     ! Dimensionality of L2AUXData_T%values

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
    ! The dimensions for the quantity
    type (L2AUX_Dimension_T), dimension(L2AUXRank) :: DIMENSIONS
    character(len=48)                              :: DIM_Names ! ','-separated
    character(len=48)                              :: DIM_Units ! ','-separated
    ! The values of the quantity
    real(r8), pointer, dimension(:,:,:) :: VALUES=>NULL()
    character(len=24)                              :: VALUE_Units
  end type L2AUXData_T

contains ! =====     Public Procedures     =============================

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
    ! Quantities can have upto three valid dimensions.  l_none can be used
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
    ! integer, dimension(L2AUXRank) :: dimensionFamilies
    integer, dimension(L2AUXRank) :: dimSizes
    integer, dimension(L2AUXRank) :: dimStarts
    integer             :: quantityType
    integer :: dimIndex
    integer :: status
    integer, dimension(L2AUXRank)   :: dimEnds
    integer, dimension(3        )   :: dim_names
    character(len=16), dimension(3) :: extra_name
    character(len=16)               :: framing
    integer                         :: option_number   ! (1 or 2; see above)

    ! Executable
    if ( present(quantityTemplate) ) then
      option_number = 1
      quantityType = quantityTemplate%quantityType
    elseif ( present(inputQuantityType) ) then
      option_number = 1
      quantityType = inputQuantityType
      dimSizes = inputDimSizes
      dimStarts = inputdimStarts
    else
      call MLSMessage(MLSMSG_Error, ModuleName, &
      & 'args to SetupL2AUXData incompatible with options 1 or 2')
    endif
    ! Fill the dimensions data structure
    ! l2aux%dimensions%dimensionFamily = dimensionFamilies
    ! if ( present(quantityType) ) then
    call GetQuantityAttributes( quantityType, &
       & framing, l2aux%VALUE_Units, dim_names )
    l2aux%dimensions%dimensionFamily = dim_names
    ! endif
    l2aux%minorFrame = (framing == 'minor')
    l2aux%majorFrame = (framing == 'major')
    ! Name the dimensions (e.g. 'frequency')
    do dimIndex=1, L2AUXRank
      call GetDimString( l2aux%dimensions(dimIndex)%dimensionFamily, &
        & extra_name(dimIndex) )
      if ( present (firstMAF) ) then
        call GetDimStart( l2aux%dimensions(dimIndex)%dimensionFamily, &
        & quantityTemplate, firstMAF, dimStarts(dimIndex) )
      endif
      if ( present (noMAFs) ) then
        call GetDimSize( l2aux%dimensions(dimIndex)%dimensionFamily, &
        & quantityTemplate, noMAFs, dimSizes(dimIndex) )
      endif
    enddo
    call Array2List(extra_name, l2aux%DIM_Names)
    ! Name the dimensions' units (e.g. 'K')
    do dimIndex=1, L2AUXRank
      call GetQuantityAttributes( l2aux%dimensions(dimIndex)%dimensionFamily, &
       & framing, extra_name(dimIndex), dim_names )
    enddo
    call Array2List(extra_name, l2aux%DIM_Units)
    l2aux%dimensions%noValues = dimSizes

    dimEnds = dimStarts + max(1,dimSizes) - 1

    ! Allocate the values for each dimension
    do dimIndex = 1, L2AUXRank
      ! if ( dimensionFamilies(dimIndex)/=L_None ) then
      if ( l2aux%dimensions(dimIndex)%dimensionFamily /= L_None ) then
        allocate (l2aux%dimensions(dimIndex)%values( &
          & dimStarts(dimIndex):dimEnds(dimIndex)), &
          & STAT=status)
        if ( status/=0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
          & MLSMSG_Allocate // "l2aux dimension values" )
        if ( present (quantityTemplate) ) then
          call GetDimValues( l2aux%dimensions(dimIndex)%dimensionFamily, &
          & quantityTemplate, l2aux%dimensions(dimIndex)%values )
        endif
      else
        l2aux%dimensions(dimIndex)%noValues=1
      end if
    end do

    ! Allocate the values for the data itself

    allocate ( l2aux%values( &
      & dimStarts(1):dimEnds(1), &
      & dimStarts(2):dimEnds(2), &
      & dimStarts(3):dimEnds(3)), STAT=status )
    if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate// "l2aux values" )

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

  !--------------------------------------  ExpandL2AUXDataInPlace  -----
  subroutine ExpandL2AUXDataInPlace ( l2aux, newSize )

  ! This subroutine expands an L2AUXData_T in place, allowing the user to
  ! add more `profiles' to it.  Note that the `profile' dimension is the last
  ! one.

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

    if (l2aux%dimensions(3)%dimensionFamily==L_None) &
      call MLSMessage (MLSMSG_Error, ModuleName, &
        & "This l2aux is not expandable")

    ! Now see how long this is
    oldSize = l2aux%dimensions(3)%noValues
    ! Do a sanity check
    if ( newSize<oldSize ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "This l2aux is getting smaller not bigger" )

    ! Now expand this dimension
    temp1D => l2aux%dimensions(3)%values
    ! Nullify old one so allocate_test doesn't clobber temp1D
    nullify ( l2aux%dimensions(3)%values )
    call allocate_test(l2aux%dimensions(3)%values,newSize, &
      & 'New l2aux%dimensions(3)%values', ModuleName)
    l2aux%dimensions(3)%noValues = newSize
    l2aux%dimensions(3)%values(1:oldSize)=temp1D
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
    l2aux%values(:,:,1:oldSize) = temp3d

    ! Now we can loose temp3d
    call deallocate_test ( temp3d, "temp3d", ModuleName )

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

    AddL2AUXToDatabase = newSize
  end function AddL2AUXToDatabase

  ! ---------------------------------------  DestroyL2AUXDatabase  -----
  subroutine DestroyL2AUXDatabase ( DATABASE )

  ! This subroutine destroys the l2aux database

    ! Dummy argument
    type (L2AUXData_T), dimension(:), pointer :: DATABASE

    ! Local variables
    integer :: l2auxIndex, status

    if ( associated(database) ) then
      do l2auxIndex = 1, size(database)
        call DestroyL2AUXContents ( database(l2auxIndex) )
      end do
      deallocate ( database, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Warning,ModuleName, &
        & MLSMSG_DeAllocate // "database")
    end if
  end subroutine DestroyL2AUXDatabase

  ! ------------------------------------------ Dump_L2AUX_DataBase ------------

  subroutine Dump_L2AUX_DataBase ( L2aux, Name, Details )

    ! Dummy arguments
    type (l2auxData_T), intent(in) ::          L2AUX(:)
    character(len=*), intent(in), optional :: Name
    integer, intent(in), optional :: DETAILS

    ! Local variables
    integer :: i
    
    call output ( '============ L2AUX Data Base ============', advance='yes' )
    call output ( ' ', advance='yes' )
    if ( present(name) ) then
      call output ( 'L2AUX Database name: ', advance='no' )
      call output ( name, advance='yes' )
    endif
    if ( size(l2aux) < 1 ) then
      call output ( '**** L2AUX Database empty ****', advance='yes' )
      return
    endif
    do i = 1, size(l2aux)
      call dump(l2aux(i), Details)
    end do
      
  end subroutine Dump_L2AUX_DATABASE

  ! ------------------------------------------ Dump_L2AUX ------------

  subroutine Dump_L2AUX ( L2aux, Details )

    ! Dummy arguments
    type (l2auxData_T), intent(in) ::          L2AUX
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump multidim arrays
    !                                        ! -1 Skip even 1-d arrays
    !                                        ! -2 Skip all but name
    !                                        ! >0 Dump even multi-dim arrays
    !                                        ! Default 1

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
          & advance='yes', ierr=ierr ) 
        if ( ierr /= 0 ) call output ( '(not found in string table)', &
          & advance='yes')
        call output ( '    (its index): ')
        call output ( l2aux%instrumentmodule, advance='no')
      endif
      call output ( ' ', advance='yes')
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
          call dump ( l2aux%dimensions(dim)%values, 'dim values:' )
         else
          call output ( ' is not associated', advance='yes')
         endif
      enddo
      if ( myDetails < 1 ) return
      call dump ( l2aux%values, 'values:' )
 
  end subroutine Dump_L2AUX
    
  !------------------------------------------------ ReadL2AUXData ------------
  subroutine ReadL2AUXData(sd_id, quantityname, quantityType, l2aux, firstProf, lastProf, &
    & checkDimNames, hdfVersion)

  use MLSFiles, only: HDFVERSION_4, HDFVERSION_5

    ! This routine reads an l2aux file, returning a filled data structure
    ! and the number of profiles read.

    ! Arguments

    character (LEN=*), intent(IN) :: quantityname ! Name of L2AUX quantity = sdname in writing routine
    integer, intent(IN) :: sd_id ! Returned by sfstart before calling us
    integer, intent(in) :: QuantityType ! Lit index
    integer, intent(IN), optional :: firstProf, lastProf ! Defaults to first and last
    type( L2AUXData_T ), intent(OUT) :: l2aux ! Result
    logical, optional, intent(in) :: checkDimNames
    integer, intent(in), optional :: hdfVersion
    ! Local variables
    integer :: myhdfVersion
    ! Executable code
    myhdfVersion = default_hdfversion_read
    if ( present(hdfVersion) ) myhdfVersion = hdfVersion
    select case (myhdfVersion)
    case (HDFVERSION_4)
      call ReadL2AUXData_hdf4(sd_id, quantityname, quantityType, l2aux, firstProf, lastProf, &
    & checkDimNames)
    case (HDFVERSION_5)
      call ReadL2AUXData_hdf5(sd_id, quantityname, quantityType, l2aux, firstProf, lastProf, &
    & checkDimNames)
    case default
    end select

  end subroutine ReadL2AUXData

  !------------------------------------------------ ReadL2AUXData_hdf4 ------------
  subroutine ReadL2AUXData_hdf4(sd_id, quantityname, quantityType, l2aux, firstProf, lastProf, &
    & checkDimNames )

    ! This routine reads an l2aux file, returning a filled data structure and the !
    ! number of profiles read.

    ! Arguments

    character (LEN=*), intent(IN) :: quantityname ! Name of L2AUX quantity = sdname in writing routine
    integer, intent(IN) :: sd_id ! Returned by sfstart before calling us
    integer, intent(in) :: QuantityType
    integer, intent(IN), optional :: firstProf, lastProf ! Defaults to first and last
    type( L2AUXData_T ), intent(OUT) :: l2aux ! Result
    logical, optional, intent(in) :: checkDimNames

    ! Parameters

    character (LEN=*), parameter :: SZ_ERR = 'Failed to get size of &
      &dimension '
    character (LEN=*), parameter :: MLSMSG_INPUT = 'Error in input argument '
    character (LEN=*), parameter :: MLSMSG_l2auxRead = 'Unable to read l2aux &
      &field:'
    integer, parameter :: MAXRANK = 3
    integer, parameter :: MAXDIMSIZES = 300
    logical            :: myCHECKDIMNAMES	! .TRUE. only for actual l2auxfiles

    ! Functions

    !INTEGER, EXTERNAL :: swattach, swdetach, swdiminfo, swinqdims, swrdfld

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

    logical :: firstCheck, lastCheck
    real (r4), dimension(:,:,:), pointer :: TMPVALUES

    myCHECKDIMNAMES = .false.
    if ( present(checkDimNames) ) myCHECKDIMNAMES = checkDimNames
    ! Attach to the file for reading

    ! find SD data set identifier
    sds_index = sfn2index(sd_id, quantityname)
    if (sds_index == -1) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Failed to get sds_index for '//trim(quantityName) )

    sds_id = sfselect(sd_id, sds_index)
    if (sds_id == -1) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Failed to get sds_id.' )

    status = sfginfo(sds_id, sds_name, rank, file_dim_sizes, data_type, &
      & num_attrs)

    if (status == -1) then
      call MLSMessage(MLSMSG_Error, ModuleName, 'Failed to get sf info.')
    elseif (sds_name /= quantityname) then
      call MLSMessage(MLSMSG_Error, ModuleName, &
        & 'quantityname  fails to match sf info.')
    endif

    ! Check optional input arguments

    firstCheck = present(firstProf)
    lastCheck = present(lastProf)

    ! Uncertain what to do with those just yet
    ! Now find dimension family of dimension; e.g., MAF
    dim_families = l_none
    data_dim_sizes = 1
    file_dim_sizes = 1

    do dim=1, rank
      write(dim_char, '(I1)') dim
      dim_id = sfdimid(sds_id, dim-1)		! dim starts from 0
      if(dim_id == -1) then
        msr = 'Failed to get dim_id for dim index number ' // dim_char
        call MLSMessage(MLSMSG_Error, ModuleName, msr)
      else
        status = sfgdinfo(dim_id, dim_name, dim_size1, data_type, &
          & num_attrs)
        if(status == -1) then
          msr = 'Failed to get dim_info for dim index number ' // &
            & dim_char
          call MLSMessage(MLSMSG_Error, ModuleName, msr)
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
                & 'Unrecognized dimension in l2aux:'//trim(dim_name) )
            else
              dim_families(dim) = l_channel
              data_dim_sizes(dim) = dim_size1
            endif
          end select
        endif
      endif
    enddo

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
      & 'Failed to read SD.')
    l2aux%values = reshape ( tmpValues, &
      & (/ data_dim_sizes(1), data_dim_sizes(2), data_dim_sizes(3) /) )

    call Deallocate_test ( tmpValues, 'tmpValues', ModuleName )

    ! Deallocate local variables


    ! Terminate access to the data set

    status = sfendacc(sds_id)
    if (status == -1) call MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
      &end access to sds_id after reading.')

    !  After reading, detach from hdf interface

    !     status = sfend(sd_id)
    !     if (status == -1) call MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
    !          &detach from SD file after reading.')

  end subroutine ReadL2AUXData_hdf4


  !------------------------------------------------ ReadL2AUXData_hdf5 ------------
  subroutine ReadL2AUXData_hdf5(sd_id, quantityname, quantityType, l2aux, firstProf, lastProf, &
    & checkDimNames)

    ! This routine reads an l2aux file, returning a filled data structure and the !
    ! number of profiles read.
    ! Assumptions
    ! The data format is radiance-like, single- or double-precision
    ! You don't really care which names the original dimensions had
    ! You want the new dimension families to be 'channel', 'MIF', 'MAF'
    ! and in that order

    ! Arguments

    character (LEN=*), intent(IN) :: quantityname ! Name of L2AUX quantity = sdname in writing routine
    integer, intent(IN) :: sd_id ! Returned by sfstart before calling us
    integer, intent(in) :: QUANTITYTYPE ! Lit index
    integer, intent(IN), optional :: firstProf, lastProf ! Defaults to first and last
    type( L2AUXData_T ), intent(OUT) :: l2aux ! Result
    logical, optional, intent(in) :: checkDimNames

    ! Parameters
    type(l1bdata_t)               :: L1BDATA ! Intermediate Result
    integer                       :: NoMAFs
    logical, parameter            :: NEVERFAIL = .TRUE.
    integer, dimension(L2AUXRank) :: dim_families
    integer, dimension(L2AUXRank) :: data_dim_sizes
    integer                       :: status
    ! Executable
    ! call MLSMessage ( MLSMSG_Error,ModuleName, &
    !      & 'Sorry--unable to read hdf5-formatted l2aux files yet' )
    CALL ReadL1BData(sd_id, QuantityName, L1bData, NoMAFs, status, &
    & FirstMAF=firstProf, LastMAF=lastProf, NEVERFAIL=NEVERFAIL)
    if ( status /= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read ' &
      & // trim(QuantityName) // ' (perhaps too unlike a radiance)' )

    dim_families(1) = l_channel                      
    data_dim_sizes = shape(L1BDATA%DpField)          
    dim_families(2) = l_mif                          
    dim_families(3) = l_maf                          
!   call SetupNewl2auxRecord ( dim_families, data_dim_sizes, (/1,1,1/), l2aux )
    call SetupNewl2auxRecord ( l2aux, inputDimFamilies=dim_families, &
     & inputDimSizes=data_dim_sizes, inputDimStarts=(/1,1,1/), inputQuantityType=quantityType )
    l2aux%values = L1BDATA%DpField
    deallocate(L1BDATA%DpField, stat=status)
    if ( status /= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to deallocate l1bdata%dpField while reading' &
      & // trim(QuantityName) )
  end subroutine ReadL2AUXData_hdf5

  !----------------------------------------------------- WriteL2AUXData ------

  subroutine WriteL2AUXData(l2aux, l2FileHandle, returnStatus, sdName, &
    & NoMAFS, WriteCounterMAF, DimNames, Reuse_dimNames, hdfVersion)

  use MLSFiles, only: HDFVERSION_4, HDFVERSION_5

  ! Write l2aux to the file with l2FileHandle
  ! Optionally, write a bogus CounterMAF sd so the
  ! resulting file can masquerade as an l1BRad
  ! (Note that this bogus sd should only be written once for each file)
    type (L2AUXData_T), intent(in) :: L2AUX
    integer, intent(in) :: L2FILEHANDLE                 ! From h5fopen or sfstart
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
    ! Executable code
    myhdfVersion = default_hdfversion_write
    if ( present(hdfVersion) ) myhdfVersion = hdfVersion
    select case (myhdfVersion)
    case (HDFVERSION_4)
      call WriteL2AUXData_hdf4(l2aux, l2FileHandle, returnStatus, sdName, &
    & NoMAFS, WriteCounterMAF, DimNames, Reuse_dimNames)
    case (HDFVERSION_5)
      call WriteL2AUXData_hdf5(l2aux, l2FileHandle, returnStatus, sdName, &
    & NoMAFS, WriteCounterMAF, DimNames, Reuse_dimNames)
    case default
    end select
  end subroutine WriteL2AUXData

  !----------------------------------------------------- WriteL2AUXData_hdf5 ------

  subroutine WriteL2AUXData_hdf5(l2aux, l2FileHandle, returnStatus, sdName, &
    & NoMAFS, WriteCounterMAF, DimNames, Reuse_dimNames)
  ! Write l2aux to the file with l2FileHandle
  ! Optionally, write a bogus CounterMAF sd so the
  ! resulting file can masquerade as an l1BRad
  ! (Note that this bogus sd should only be written once for each file;
  !  also note the attempt to convert l2aux%values to KIND of l1b radiances)
  use MLS_DataProducts, only: DATAPRODUCTS_T
  use MLSAuxData, only: BUILD_MLSAUXDATA
  use MLSHDF5, only: SaveAsHDF5DS
  use PCFHdr, only: h5_writeglobalattr

    type (L2AUXData_T), intent(in) :: L2AUX
    integer, intent(in) :: L2FILEHANDLE                 ! From h5fopen
    character (len=*), optional, intent(in) :: SDNAME ! Defaults to l2aux%name
    character (len=*), optional, intent(in) :: DimNames ! Comma-separated list
                                                        ! Otherwise automatic
                                                        ! (Requiring l2cf)
    integer, intent(in), optional :: NoMAFS
    logical, intent(in), optional :: WriteCounterMAF  ! Write bogus CounterMAF
    logical, intent(in), optional :: Reuse_dimNames   ! We already wrote them
    integer, intent(out) :: returnStatus           ! 0 unless error

    ! Local variables
    integer :: myNoMAFS, MAF
    integer, dimension(3) :: dims
    type(DataProducts_T) :: dataProduct
    logical :: myWriteCounterMAF
    logical :: myReuse_dimNames
    integer, dimension(:), pointer :: CounterMAF ! bogus array
    logical, parameter             :: ALWAYSWRITEAS32BITS = .true.

    ! Executable code
    returnStatus = 0
    myWriteCounterMAF = .false.
    if ( present(WriteCounterMAF) ) myWriteCounterMAF = WriteCounterMAF
    myReuse_dimNames = .false.
    if ( present(Reuse_dimNames) ) myReuse_dimNames = Reuse_dimNames
    myNoMAFS = 1
    if ( any(l2aux%dimensions%dimensionFamily == L_MAF) ) &
     & myNoMAFS = l2aux%dimensions(3)%noValues
    if ( present(NoMAFS) ) myNoMAFS = NoMAFS
    
	 !  call announce_error (0,&
    !  & "hdf5 version of WriteL2AUXData_hdf5 not ready yet " )
    !  returnStatus = 1
    if ( .not. associated ( l2aux%values ) ) then
	   call announce_error (0,&
        & "l2aux values not associated yet " )
      returnStatus = 1
    else
      if ( present(sdName) ) then
        dataProduct%name = sdName
      else
        call get_string ( l2aux%name, dataProduct%name, strip=.true. )
      endif
      if ( myWriteCounterMAF .or. ALWAYSWRITEAS32BITS ) then
        dataProduct%data_type = 'real'    ! same type as l1bradiances
      else
        dataProduct%data_type = 'double'  ! type of L2AUXData_T%values
      endif
      dims(1) = size(l2aux%values, 1)
      dims(2) = size(l2aux%values, 2)
      dims(3) = size(l2aux%values, 3)
      ! call Dump_L2AUX(l2AUX)
      if ( myWriteCounterMAF .or. ALWAYSWRITEAS32BITS ) then
        ! call Build_MLSAuxData(l2FileHandle, dataProduct, real(l2aux%values, r4))
        call SaveAsHDF5DS (l2FileHandle, trim(dataProduct%name), &
         & real(l2aux%values, r4))
        call WriteL2AUXAttributes(l2FileHandle, l2aux, trim(dataProduct%name))
      else
        ! call Build_MLSAuxData(l2FileHandle, dataProduct, l2aux%values)
        ! call SaveAsHDF5DS (l2FileHandle, trim(dataProduct%name), l2aux%values)
      endif
      call h5_writeglobalattr(l2FileHandle, skip_if_already_there=.false.)
      if ( .not. myWriteCounterMAF ) return
    
      ! Now create and write bogus counterMAF array
      if ( myNoMAFS < 1 ) then
        call announce_error(0, &
        & "Too few MAFs to fake CounterMAFs in l2aux file:  " )
        return
      endif
      nullify (CounterMAF)
      call allocate_test(CounterMAF,myNoMAFS,'counterMAF',ModuleName)
      dims(1) = myNoMAFS
      dataProduct%name = 'counterMAF'
      dataProduct%data_type = 'integer'
      ! sdId= SFcreate ( l2FileHandle, 'counterMAF', DFNT_INT8, &
      !  & 1, dimSizes)
      do MAF=0, myNoMAFS-1
        counterMAF(MAF+1) = MAF
      enddo
      ! status= SFWDATA_F90(sdId, (/ 0 /), &
      !  & (/ 1 /) , dimSizes, CounterMAF)
      call Build_MLSAuxData(l2FileHandle, dataProduct, counterMAF, &
      & myNoMAFS )
      call Deallocate_Test(CounterMAF,"CounterMAF",ModuleName)
    endif
  end subroutine WriteL2AUXData_hdf5

  !----------------------------------------------------- WriteL2AUXData_hdf4 ------

  subroutine WriteL2AUXData_hdf4(l2aux, l2FileHandle, returnStatus, sdName, &
    & NoMAFS, WriteCounterMAF, DimNames, Reuse_dimNames)
  ! Write l2aux to the file with l2FileHandle
  ! Optionally, write a bogus CounterMAF sd so the
  ! resulting file can masquerade as an l1BRad
  ! (Note that this bogus sd should only be written once for each file)
    type (L2AUXData_T), intent(in) :: L2AUX
    integer, intent(in) :: L2FILEHANDLE
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
    endif

    goodDim=l2aux%dimensions%dimensionFamily /= L_None
    noDimensionsUsed = min(COUNT(goodDim), myL2AUXRank)
    call allocate_test(dimSizes,noDimensionsUsed,'dimSizes',ModuleName)
    if ( present(DimNames) ) then
      dimSizes = PACK(l2aux%dimensions(1:noDimensionsUsed)%noValues, &
        & goodDim(1:noDimensionsUsed))
    else
      dimSizes = PACK(l2aux%dimensions%noValues, &
        & goodDim)
    endif

    ! Create the sd within the file
    sdId= SFcreate ( l2FileHandle, nameString, DFNT_FLOAT32, &
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
        elseif ( (dimensionInData > 1) .and. (l2aux%minorFrame) ) then
          call GetModuleName( l2aux%instrumentModule, dimName)
          if (len_trim(dimName) < len(dimName)) dimName=TRIM(dimName)//'.'
          call get_string (lit_indices(l2aux%dimensions(dimensionInData)%dimensionFamily), &
            & dimName(LEN_TRIM(dimName)+1:))
        else
          call get_string (l2aux%name, dimName, strip=.true. )
          if (len_trim(dimName) < len(dimName)) dimName=TRIM(dimName)//'.'
          call get_string (lit_indices(l2aux%dimensions(dimensionInData)%dimensionFamily), &
            & dimName(LEN_TRIM(dimName)+1:))
        endif
        ! Write dimension name
        status=SFSDMName(dimID,TRIM(dimName))
        if ( status /= 0 ) then
		  		call output("dim name: ")
		  		call output(TRIM(dimName), advance='yes')
		  		call announce_error (  0, &
          & "Error setting dimension name to SDS l2aux file:")
!		  		call MLSMessage ( MLSMSG_Error, ModuleName, &
!          & "Error setting dimension name to SDS l2aux file:")
	     endif
        ! Write dimension scale
        status=SFSDScale(dimID, dimSizes(dimensionInFile+1), DFNT_FLOAT32, &
          & l2aux%dimensions(dimensionInData)%values)
        if ( status /= 0 ) then
		  		call output("dimID: ")
		  		call output(dimID, advance='yes')
		  		call output("dim name: ")
		  		call output(TRIM(dimName), advance='yes')
		      call announce_error ( 0, &
          & "Error writing dimension scale in l2auxFile:" )
!		      call MLSMessage ( MLSMSG_Error, ModuleName, &
!          & "Error writing dimension scale in l2auxFile:" )
	     endif
        dimensionInFile=dimensionInFile+1
      endif
    end do
    endif

    ! Now write the data
    ! Make sure the data can be placed in a real
    status= SFWDATA_F90(sdId, start(1:noDimensionsUsed), &
      & stride(1:noDimensionsUsed), dimSizes, real ( &
      & max ( -hugeR4, min ( hugeR4, l2aux%values ) ) ) )
    if ( status /= 0 ) then
	   call announce_error (0, &
      & "Error writing SDS data to  l2aux file:  " )
    endif

    call Deallocate_Test(dimSizes,"dimSizes",ModuleName)
    
    ! Terminate access to sd
    status = sfendacc(sdId)      ! 
    if ( status /= 0 ) then
	   call announce_error (0,&
      & "Error ending access to the sd  " )
    endif
    returnStatus = error
    if ( .not. myWriteCounterMAF ) return
    
    ! Now create and write bogus counterMAF array
    if ( myNoMAFS < 1 ) then
      call announce_error(0, &
      & "Too few MAFs to fake CounterMAFs in l2aux file:  " )
      return
    endif
    nullify (CounterMAF, dimSizes)
    call allocate_test(CounterMAF,myNoMAFS,'counterMAF',ModuleName)
    call allocate_test(dimSizes,1,'dimSizes',ModuleName)
    dimSizes(1) = myNoMAFS
    sdId= SFcreate ( l2FileHandle, 'counterMAF', DFNT_INT32, &
      & 1, dimSizes)
    do MAF=0, myNoMAFS-1
      counterMAF(MAF+1) = MAF
    enddo
    status= SFWDATA_F90(sdId, (/ 0 /), &
      & (/ 1 /) , dimSizes, CounterMAF)
    if ( status /= 0 ) then
	   call announce_error (0,&
      & "Error writing counterMAF data to  l2aux file:  " )
    endif
    call Deallocate_Test(dimSizes,"dimSizes",ModuleName)
    call Deallocate_Test(CounterMAF,"CounterMAF",ModuleName)
    
    ! Terminate access to sd
    status = sfendacc(sdId)
    if ( status /= 0 ) then
	   call announce_error (0,&
      & "Error ending access to the sd  " )
    endif
    returnStatus = error
    ! call Dump_L2AUX(l2AUX)

  end subroutine WriteL2AUXData_hdf4

  ! ----------------------------------  WriteL2AUXAttributes  -----
  subroutine WriteL2AUXAttributes ( L2FileHandle, l2aux, name)
  use MLSHDF5, only: MakeHDF5Attribute
  ! Writes the pertinent attributes for an l2aux
  ! Arguments
  integer, intent(in) :: L2FileHandle
  type (L2AUXData_T), intent(in) :: L2AUX
  character(len=*) :: name
  ! Internal variables
  integer :: dim
  integer :: ndims
  character(len=16), dimension(L2AUXRank) :: dim_name
  character(len=16), dimension(L2AUXRank) :: dim_unit
  character(len=16) :: dim_of_i
  character(len=16) :: framing
  character(len=2) :: i_char
  character(len=*), parameter :: ottff = '1,2,3,4,5'
  ! Executable
  call MakeHDF5Attribute(L2FileHandle, name, 'Title', name)
  call MakeHDF5Attribute(L2FileHandle, name, 'Units', &
    & trim(l2aux%VALUE_Units))
  call MakeHDF5Attribute(L2FileHandle, name, 'DimensionNames', &
    & trim(l2aux%DIM_Names))
  if ( l2aux%majorframe) then
    framing = 'major'
  elseif ( l2aux%minorframe) then
    framing = 'minor'
  else
    framing = 'neither'
  endif
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
      endif
    endif
  enddo
  end subroutine WriteL2AUXAttributes

  ! ----------------------------------  GetDimSize  -----
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

  ! ----------------------------------  GetDimStart  -----
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

  ! ----------------------------------  GetDimString  -----
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

  ! ----------------------------------  AreDimValuesNonTrivial  -----
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

  ! ----------------------------------  GetDimValues  -----
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
      enddo
    case ( l_frequency, &
         & L_IntermediateFrequency, l_USBFrequency, L_LSBFrequency )  
      dim_values = quantityTemplate%frequencies
    case ( l_geodAngle )  
      dim_values = quantityTemplate%phi(1,:)   ! If not stacked, phi(surf, :)
    case ( l_height, l_pressure )  
      dim_values = quantityTemplate%surfs(:,1) ! If incoherent, surfs(:, inst)

    end select

  end subroutine GetDimValues

  ! ----------------------------------  GetQuantityAttributes  -----
  subroutine GetQuantityAttributes ( quantityType, &
   & framing, units_name, dim_names)

  ! Given a quantity type, e.g. l_vmr,
  ! returns major/minor/neither framing distinction, default unit name,
  ! and the 3 dimension names, e.g. (/l_channel, l_MIF, l_MAF/)

    ! Dummy arguments
    integer, intent(in)                :: quantityType
    character(len=*), intent(out)      :: framing
    character(len=*), intent(out)      :: units_name
    integer, dimension(3), intent(out) :: dim_names

    ! Executable code
    units_name = ' '
    select case (quantityType)                                       
    case ( l_channel )  
      framing = 'major'
      units_name = 'channel'
      dim_names = (/ l_none, l_none, l_none /)                  
    case ( l_chisqchan )  
      framing = 'major'
      dim_names = (/ l_channel, l_none, l_MAF /)                  
    case ( l_chisqmmaf )  
      framing = 'major'
      dim_names = (/ l_none, l_none, l_MAF /)                  
    case ( l_chisqmmif )  
      framing = 'minor'
      dim_names = (/ l_none, l_MIF, l_MAF /)                  
    case ( l_chunk )  
      framing = 'major'
      units_name = 'chunk'
      dim_names = (/ l_none, l_none, l_none /)                  
    case ( l_cloudInducedRadiance )  
      framing = 'minor'
      units_name = 'K'
      dim_names = (/ l_channel, l_MIF, l_MAF /)                  
    case ( l_cloudExtinction )  
      framing = 'neither'
      dim_names = (/ l_channel, l_none, l_MAF /)                  
    case ( l_cloudRadSensitivity )  
      framing = 'minor'
      units_name = 'K'
      dim_names = (/ l_channel, l_none, l_MAF /)                  
    case ( l_cloudWater )  
      framing = 'neither'
      dim_names = (/ l_channel, l_none, l_MAF /)                  
    case ( l_dnwt_ajn )  
      framing = 'neither'
      dim_names = (/ l_none, l_iteration, l_chunk /)                  
    case ( l_dnwt_axmax )  
      framing = 'neither'
      dim_names = (/ l_none, l_iteration, l_chunk /)                  
    case ( l_dnwt_cait, l_dnwt_chiSqMinNorm, l_dnwt_chiSqNorm, l_dnwt_diag, &
      & l_dnwt_dxdx, l_dnwt_dxdxl, l_dnwt_dxn, l_dnwt_dxnl, l_dnwt_flag, &
      & l_dnwt_fnmin, l_dnwt_fnorm, l_dnwt_gdx, l_dnwt_gfac, l_dnwt_gradn, &
      & l_dnwt_sq, l_dnwt_sqt )  
      framing = 'neither'
      dim_names = (/ l_none, l_iteration, l_chunk /)                  
    case ( l_effectiveOpticalDepth )  
      framing = 'minor'
      dim_names = (/ l_channel, l_MIF, l_MAF /)                  
    case ( l_elevOffset )  
      framing = 'neither'
      units_name = 'deg'
      dim_names = (/ l_channel, l_MIF, l_MAF /)                  
    case ( l_frequency )  
      framing = 'major'
      units_name = 'frequency'
      dim_names = (/ l_none, l_none, l_none /)                  
    case ( l_geodAngle )  
      framing = 'neither'
      units_name = 'deg'
      dim_names = (/ l_none, l_none, l_none /)                  
    case ( l_heightOffset )  
      framing = 'minor'
      units_name = 'deg'
      dim_names = (/ l_channel, l_MIF, l_MAF /)                  
    case ( l_iteration )  
      framing = 'major'
      units_name = 'iteration'
      dim_names = (/ l_none, l_none, l_none /)                  
    case ( l_jacobian_cols )  
      framing = 'neither'
      dim_names = (/ l_none, l_iteration, l_chunk /)                  
    case ( l_jacobian_rows )  
      framing = 'neither'
      dim_names = (/ l_none, l_iteration, l_chunk /)                  
    case ( l_losTransFunc )  
      framing = 'neither'
      dim_names = (/ l_frequency, l_MIF, l_MAF /)                  
    case ( l_losVel )  
      framing = 'minor'
      dim_names = (/ l_xyz, l_MIF, l_MAF /)                  
    case ( l_MAF )  
      framing = 'major'
      units_name = 'MAF'
      dim_names = (/ l_none, l_none, l_none /)                  
    case ( l_massMeanDiameterIce )  
      framing = 'neither'
      dim_names = (/ l_channel, l_MIF, l_MAF /)                  
    case ( l_massMeanDiameterWater )  
      framing = 'neither'
      dim_names = (/ l_channel, l_MIF, l_MAF /)                  
    case ( l_MIF )  
      framing = 'major'
      units_name = 'MIF'
      dim_names = (/ l_none, l_none, l_none /)                  
    case ( l_noiseBandwidth )  
      framing = 'neither'
      units_name = 'MHz'
      dim_names = (/ l_channel, l_none, l_none /)                  
    case ( l_numJ )  
      framing = 'neither'
      dim_names = (/ l_none, l_iteration, l_chunk /)                  
    case ( l_opticalDepth )  
      framing = 'minor'
      dim_names = (/ l_channel, l_MIF, l_MAF /)                  
    case ( l_orbitInclination )  
      framing = 'minor'
      units_name = 'deg'
      dim_names = (/ l_none, l_none, l_none /)                  
    case ( l_phiTan )  
      framing = 'minor'
      units_name = 'deg'
      dim_names = (/ l_none, l_MIF, l_MAF /)                  
    case ( l_ptan )  
      framing = 'minor'
      units_name = 'log10(hPa)'
      dim_names = (/ l_none, l_MIF, l_MAF /)                  
    case ( l_radiance )  
      framing = 'minor'
      units_name = 'K'
      dim_names = (/ l_channel, l_MIF, l_MAF /)                  
    case ( l_sizedistribution )  
      framing = 'neither'
      dim_names = (/ l_channel, l_MIF, l_MAF /)                  
    case ( l_scanResidual )  
      framing = 'minor'
      units_name = 'm'
      dim_names = (/ l_none, l_MIF, l_MAF /)                  
    case ( l_scECI )  
      framing = 'minor'
      units_name = 'm'
      dim_names = (/ l_xyz, l_MIF, l_MAF /)                  
    case ( l_scVel )  
      framing = 'minor'
      units_name = 'm/s'
      dim_names = (/ l_xyz, l_MIF, l_MAF /)                  
    case ( l_scVelECI )  
      framing = 'minor'
      units_name = 'm/s'
      dim_names = (/ l_xyz, l_MIF, l_MAF /)                  
    case ( l_scVelECR )  
      framing = 'minor'
      units_name = 'm/s'
      dim_names = (/ l_xyz, l_MIF, l_MAF /)                  
    case ( l_scGeocAlt )  
      framing = 'minor'
      units_name = 'm'
      dim_names = (/ l_xyz, l_MIF, l_MAF /)                  
    case ( l_limbSidebandFraction )  
      framing = 'neither'
      dim_names = (/ l_channel, l_none, l_none /)                  
    case ( l_spaceRadiance )  
      framing = 'neither'
      units_name = 'K'
      dim_names = (/ l_none, l_none, l_none /)                  
    case ( l_surfacetype )  
      framing = 'neither'
      dim_names = (/ l_none, l_none, l_none /)                  
    case ( l_systemTemperature )  
      framing = 'neither'
      units_name = 'K'
      dim_names = (/ l_channel, l_none, l_none /)                  
    case ( l_tngtECI )  
      framing = 'minor'
      units_name = 'm'
      dim_names = (/ l_xyz, l_MIF, l_MAF /)                  
    case ( l_tngtGeodAlt )  
      framing = 'minor'
      units_name = 'm'
      dim_names = (/ l_xyz, l_MIF, l_MAF /)                  
    case ( l_tngtGeocAlt )  
      framing = 'minor'
      units_name = 'm'
      dim_names = (/ l_xyz, l_MIF, l_MAF /)                  
    case ( l_totalExtinction )  
      framing = 'neither'
      dim_names = (/ l_channel, l_none, l_MAF /)                  
    case ( l_vmr )  
      framing = 'neither'
      dim_names = (/ l_channel, l_none, l_MAF /)                  
      units_name = 'vmr'
    case ( l_xyz )  
      framing = 'major'
      units_name = 'xyz'
      dim_names = (/ l_none, l_none, l_none /)                  
    case default                                                     
      framing = 'neither'
      dim_names = (/ l_channel, l_MIF, l_MAF /)                  
    end select                                                       

  end subroutine GetQuantityAttributes

  ! ---------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( WHERE, full_message, CODE )
    integer, intent(in) :: WHERE   ! Tree node where error was noticed
	character(LEN=*), intent(in)    :: full_message
    integer, intent(in), optional :: CODE    ! Code for error message

    error = max(error,1)
    call output ( '***** At ' )
!    call print_source ( source_ref(where) )
    if ( where > 0 ) then
      call print_source ( source_ref(where) )
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
  end subroutine ANNOUNCE_ERROR

!=============================================================================
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module L2AUXData
!=============================================================================

!
! $Log$
! Revision 2.53  2003/05/29 16:43:02  livesey
! Renamed sideband fraction
!
! Revision 2.52  2003/05/12 02:06:32  livesey
! Bound r8->r4 conversion
!
! Revision 2.51  2003/04/25 19:55:09  livesey
! Added more useful error message
!
! Revision 2.50  2003/03/07 00:42:13  pwagner
! Abbreviated Units names; removed spaces from attribute names
!
! Revision 2.49  2003/02/21 23:42:21  pwagner
! Also writes Fill Value attribute
!
! Revision 2.48  2003/02/12 21:52:34  pwagner
! Renames blank dim units to none
!
! Revision 2.47  2003/02/07 21:44:56  pwagner
! Capitalized 1st letter of each attribute name
!
! Revision 2.46  2003/01/30 01:02:28  pwagner
! Writing attributes for hdf5 files; global and data set
!
! Revision 2.45  2003/01/18 02:37:03  livesey
! Made the readl2aux data stuff work from the l2cf by adding the
! quantityType argument.
!
! Revision 2.44  2003/01/17 23:11:26  pwagner
! Moved most ops out of LoinL2AUXData to SetupL2AUXData
!
! Revision 2.43  2003/01/14 00:41:43  pwagner
! Added GetQuantityAttributes and getDimString; new fields in L2AUXData_T
!
! Revision 2.42  2002/12/10 00:41:28  pwagner
! In principle can now read hdf5-formatted l2aux files; untested; njl has other plans
!
! Revision 2.41  2002/12/07 00:25:42  pwagner
! Using SaveAsHDF5DS to write l2aux%values; it works
!
! Revision 2.40  2002/12/06 01:06:13  pwagner
! Finally writes radaiance-like l2aux as hdf5 files
!
! Revision 2.39  2002/12/05 19:46:23  pwagner
! Changes to speed up compiling tree-walker
!
! Revision 2.38  2002/12/03 18:04:02  pwagner
! Repaired bug that caused WriteL2AUXData files to be tiny
!
! Revision 2.37  2002/12/02 23:42:12  pwagner
! Optional param checkDimNames to ReadL2AUXData; defaults to FALSE
!
! Revision 2.36  2002/12/02 19:11:13  pwagner
! Corrected data types of counterMAF and dimensions
!
! Revision 2.35  2002/11/29 22:46:28  livesey
! Various bug fixes / improvements.
!
! Revision 2.34  2002/11/29 18:50:07  livesey
! Initialised a variable
!
! Revision 2.33  2002/11/26 22:16:41  jonathan
! Comment-out dump_l2aux diagnostics
!
! Revision 2.32  2002/11/25 18:04:52  pwagner
! Consistent with latest changes to MLSAuxData
!
! Revision 2.31  2002/11/22 21:48:02  pwagner
! Fleshed out WriteL2AUXData_hdf5; untested yet
!
! Revision 2.30  2002/11/13 01:09:47  pwagner
! Beginnings of attempt to write hdf5 L2AUX; incomplete
!
! Revision 2.29  2002/11/08 23:14:41  pwagner
! Should work again with mlsl2
!
! Revision 2.28  2002/11/08 18:25:33  jonathan
! Changes to allow writing rank 2; also reuse_DimNames
!
! Revision 2.27  2002/11/06 02:01:06  livesey
! Changes to fill from l2aux
!
! Revision 2.26  2002/11/06 00:18:37  pwagner
! Can WriteL2AUXData w/o l2cf: useable by small utility programs
!
! Revision 2.25  2002/10/08 17:36:21  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.24  2002/08/21 01:04:53  livesey
! Changed to single precision for data
!
! Revision 2.23  2002/08/15 21:47:04  pwagner
! WriteL2AuxData now returns non-zero status if it fails
!
! Revision 2.22  2001/11/01 21:03:59  pwagner
! Uses new sfwdata_f90 generic; added toc
!
! Revision 2.21  2001/10/26 23:13:18  pwagner
! Provides a single dump module interface and details
!
! Revision 2.20  2001/10/08 23:41:27  pwagner
! Improved dump routines
!
! Revision 2.19  2001/10/05 23:32:27  pwagner
! Added majorframe to data type; trimmed unused stuff
!
! Revision 2.18  2001/08/06 18:35:24  pwagner
! Added dump_l2aux
!
! Revision 2.17  2001/07/11 20:50:46  dwu
! fix problem in readl2auxdata
!
! Revision 2.16  2001/05/30 23:53:31  livesey
! Changed for new version of L1BData
!
! Revision 2.15  2001/05/12 00:18:40  livesey
! Tidied up issues with array bounds etc.
!
! Revision 2.14  2001/05/03 20:32:19  vsnyder
! Cosmetic changes
!
! Revision 2.13  2001/05/02 22:24:20  pwagner
! Removed SDPToolkit use
!
! Revision 2.12  2001/04/26 02:44:17  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.11  2001/04/12 22:19:33  vsnyder
! Improved an error message
!
! Revision 2.10  2001/04/10 22:27:47  vsnyder
! Nullify explicitly instead of with <initialization> so as not to give
! pointers the SAVE attribute.  <initialization> is NOT executed on each
! entry to a procedure.
!
! Revision 2.9  2001/04/07 00:14:27  pwagner
! Added announce_error
!
! Revision 2.8  2001/03/15 18:42:29  livesey
! Removed quotes from dimension name prefixes
!
! Revision 2.7  2001/03/08 02:20:12  livesey
! Added strip argument to a call to get_string
!
! Revision 2.6  2001/03/06 22:40:47  livesey
! Working version
!
! Revision 2.5  2001/02/14 23:41:33  livesey
! Removed irrelevant numProfs argument
!
! Revision 2.4  2001/01/03 00:46:19  pwagner
! Changed sfgetinfo to sfginfo
!
! Revision 2.3  2000/12/04 23:34:38  vsnyder
! Move more of addItemToDatabase into the include.
!
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

