! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module L2PC_m
  !=============================================================================

  ! This module contains data types etc. for dealing with the new EMLS L2PC
  ! files.  The first version dealt with ascii files, but later versions
  ! must be HDF5.

  use Allocate_Deallocate, only: Allocate_test, Deallocate_test
  use dump_0, only: dump
  use Intrinsic, only: L_CHANNEL, L_GEODALTITUDE, L_ZETA, L_NONE, L_VMR, &
    & L_RADIANCE, L_NONE, L_INTERMEDIATEFREQUENCY, L_LATITUDE, L_FIELDAZIMUTH, &
    & L_ROWS, L_COLUMNS, L_ADOPTED, L_TEMPERATURE, L_TSCAT, Lit_Indices, &
    & PHYQ_DIMENSIONLESS, PHYQ_TEMPERATURE, PHYQ_VMR
  use HessianModule_0, only: CreateBlock, HessianElement_T, &
    & H_Absent, H_Sparse, H_Full, H_Unknown, DestroyBlock
  use HessianModule_1, only: Hessian_T, CreateBlock, DestroyHessian, CreateEmptyHessian
  use machine, only: io_error
  use ManipulateVectorQuantities, only: DOVECTORSMATCH
  use MatrixModule_0, only: M_ABSENT, M_BANDED, M_COLUMN_SPARSE, M_FULL, &
    & MATRIXELEMENT_T, M_UNKNOWN, DESTROYBLOCK
  use MatrixModule_1, only: COPYMATRIXVALUE, CREATEBLOCK, CREATEEMPTYMATRIX, &
    & DESTROYMATRIX, DUMP, DUMP_STRUCT, FINDBLOCK, GETACTUALMATRIXFROMDATABASE, &
    & MATRIX_T, MATRIX_DATABASE_T
  use MLSCommon, only: R8, R4, MLSFile_T
  use MLSFiles, only: DumpMLSFile => Dump
  use MLSMessageModule, only: MLSMESSAGE, &
    & MLSMSG_ALLOCATE, MLSMSG_DEALLOCATE, MLSMSG_ERROR
  use MLSSets, only: FindFirst
  use MLSSignals_m, only: GETSIGNALNAME
  use MLSStringLists, only: switchDetail
  use MLSStrings, only: writeIntsToChars
  use Molecules, only: L_EXTINCTION
  use MoreTree, only: GetStringIndexFromString, GetLitIndexFromString
  use Output_m, only: newLine, output, outputNamedValue
  use Parse_Signal_m, only: Parse_Signal
  use QuantityTemplates, only: ADDQUANTITYTEMPLATETODATABASE, QUANTITYTEMPLATE_T, &
    & SETUPNEWQUANTITYTEMPLATE, INFLATEQUANTITYTEMPLATEDATABASE, &
    & DESTROYQUANTITYTEMPLATECONTENTS, COPYQUANTITYTEMPLATE, NULLIFYQUANTITYTEMPLATE, &
    & COPYQUANTITYTEMPLATE
  use String_Table, only: DISPLAY_STRING, GET_STRING
  use TOGGLES, only: SWITCHES
  use Tree, only: DECORATION, NSONS, SUBTREE
  use VectorsModule, only: assignment(=), ADDVECTORTEMPLATETODATABASE, &
    & ADDVECTORTODATABASE, CONSTRUCTVECTORTEMPLATE, COPYVECTOR, CREATEVECTOR, &
    & DESTROYVECTORINFO, DUMP, NULLIFYVECTORTEMPLATE, VECTORTEMPLATE_T, VECTOR_T

  implicit NONE
  private

  public :: AddBinSelectorToDatabase, AddL2PCToDatabase, AdoptVectorTemplate, &
    & binSelector_T, BinSelectors, CreateDefaultBinSelectors, &
    & DefaultSelector_FieldAzimuth, DefaultSelector_Latitude, &
    & DestroyL2PC, DestroyL2PCDatabase, DestroyBinSelectorDatabase, &
    & Dump, FlushL2PCBins, L2PC_T, &
    & LoadMatrix, LoadVector, OutputHDF5L2PC, &
    & PopulateL2PCBin, PopulateL2PCBinByName, &
    & ReadCompleteHDF5L2PCFile

  interface DUMP
    module procedure DUMPONEL2PC, DumpL2PCDatabase, DumpL2PCFile
  end interface

  interface DUMP_PRIVATE
    module procedure DUMPONEL2PC, DumpL2PCDatabase, DumpL2PCFile, DumpL2PCInfo
  end interface
  ! This is an update to the L2PCs where we can store both Jacobians and Hessians
  ! Previously the L2PC's were just Matrix_Ts, now they're more diversified
  ! As this contains pointers to vector_T's and so on, I maintain a private
  ! set of databases of these in this module.  We can't use the main databases,
  ! as these would get destroyed at the end of each chunk.

  ! The supporting databases
  integer, dimension(:), pointer, public, save :: FileIDDatabase => NULL()
  type(QuantityTemplate_T), dimension(:), pointer, save :: L2PCQTS => NULL()
  type(VectorTemplate_T), dimension(:), pointer, save :: L2PCVTS => NULL()
  type(Vector_T), dimension(:), pointer, save :: L2PCVS => NULL()

  integer :: counterStart
  parameter ( counterStart = huge (0) / 4 )

  ! This type holds an l2pc
  type L2PC_T
    integer :: NAME                     ! The name of the L2PC bin
    type(Matrix_T) :: J                 ! The Jacobian
    logical :: GOTH                     ! Set true if also have Hessian Info
    type(Hessian_T) :: H                ! The Hessian
  end type L2PC_T

  ! The L2PC database
  type(L2PC_T), dimension(:), pointer, public, save :: L2PCDatabase => NULL()

  ! This datatype describes a selection rule for l2pc bins.
  type BinSelector_T
    integer :: selectorType             ! What quantity type does this apply to
    integer :: molecule                 ! What molecule does it apply to
    integer :: nameFragment             ! A possible name fragment
    logical :: exact                    ! We require an exact numeric match for this selector
    real(r8), dimension(2) :: heightRange ! The height range for this selector
    real(r8) :: cost                    ! The cost for that range
  end type BinSelector_T

  type(BinSelector_T), dimension(:), pointer, save :: BINSELECTORS => NULL()

  type L2PCInfo_T
    integer :: fileID     ! What is the HDF5 file ID
    integer :: binID      ! What is the groupID for the bin
    integer :: blocksID   ! What is the groupID for the blocks
    integer :: hBlocksID  ! What is the groupID for the hessian blocks
    integer, dimension(:,:), pointer :: BLOCKID => NULL()
    character(len=64) :: matrixName
  end type L2PCINFO_T
  type(L2PCInfo_T), dimension(:), pointer, save :: L2PCINFO => NULL()

  ! Default bin selectors (see CreateDefaultBinSelectors below)
  integer, parameter :: DEFAULTSELECTOR_LATITUDE = 1
  integer, parameter :: DEFAULTSELECTOR_FIELDAZIMUTH = 2
  
  logical, parameter :: DIEIFDESTROYFAILS = .false.

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ============= Public Procedures ==========================

  ! ------------------------------------ AddBinSelectorToDatabase --
  integer function AddBinSelectorToDatabase ( database, item )
    type (BinSelector_T), dimension(:), pointer :: DATABASE
    type (BinSelector_T) :: item
    ! Local variables
    type (BinSelector_T), dimension(:), pointer :: TEMPDATABASE

    include "addItemToDatabase.f9h"
    AddBinSelectorToDatabase = newSize
  end function AddBinSelectorToDatabase

  ! ------------------------------------  Add fileID to database ----
  integer function AddFileIDToDatabase ( Database, Item )

    ! This function simply adds a fileID  to a database of this type

    integer, dimension(:), pointer :: Database
    integer :: Item

    integer, dimension(:), pointer :: TempDatabase

    include "addItemToDatabase.f9h"

    AddFileIDToDatabase = newSize
  end function AddFileIDToDatabase

  ! ------------------------------------  Add l2pc to database ----
  integer function AddL2PCToDatabase ( Database, Item )

    ! This function simply adds an l2pc  to a database of said l2pc s.

    type(L2PC_T), dimension(:), pointer :: Database
    type(L2PC_T) :: Item

    type(L2PC_T), dimension(:), pointer :: TempDatabase

    include "addItemToDatabase.f9h"

    AddL2PCToDatabase = newSize
  end function AddL2PCToDatabase

  ! --------------------------------------- DumpL2PCDatabase ---------------
  subroutine DumpL2PCDatabase ( L2pcDB, details )
    ! This subroutine dumps an l2pc to stdout

    ! Dummy arguments
    type (l2pc_t), dimension(:), intent(in), target :: L2pcDB
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump multidim arrays
    !                                        ! -1 Skip even 1-d arrays
    !                                        ! -2 Skip all but size
    !                                        ! >0 Dump even multi-dim arrays
    !                                        ! Default 0
    ! Local variables
    integer :: i
    ! Executable
    do i=1, size(L2PCDB)
      call DumpOneL2PC( L2PCDB(i), details )
    enddo
  end subroutine DumpL2PCDatabase

  ! --------------------------------------- DumpL2PCFile ---------------
  subroutine DumpL2PCFile ( L2PCFile, details )
    ! This subroutine dumps an l2pc to stdout

    ! Dummy arguments
    type (MLSFile_T), intent(in) :: L2PCFile
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump multidim arrays
    !                                        ! -1 Skip even 1-d arrays
    !                                        ! -2 Skip all but size
    !                                        ! >0 Dump even multi-dim arrays
    !                                        ! Default 0
    ! Local variables
    integer :: i
    ! Executable
    call dumpMLSFile( L2PCFile, details=1 )
    call dump ( fileIDDataBase, 'file id database', format='(i10)' )
    do i=1, size(L2PCDataBase)
      if ( L2PCFile%FileID%f_id /= fileIDDataBase(i) ) cycle
      call newline
      call outputNamedValue( 'db index', i )
      call DumpOneL2PC( L2PCDataBase(i), details )
    enddo
  end subroutine DumpL2PCFile

  ! --------------------------------------- DumpOneL2PC ---------------
  subroutine DumpOneL2PC ( L2pc, details )
    ! This subroutine dumps an l2pc to stdout

    use HessianModule_1, only: Dump

    ! Dummy arguments
    type (l2pc_t), intent(in), target :: L2pc
    integer, intent(in), optional :: DETAILS ! passed to Vector and Matrix dumps

    ! Executable code
    call output( '- Dump of L2PC -', advance='yes' )
    if ( l2pc%name > 0 ) then
      call display_string ( l2pc%name, before='name: ' )
      call output ( l2pc%name, before=' (', after=')', advance='yes' )
    else
      call output( '*** Uh-oh, name not found in string table', advance='yes' )
    endif
    ! First dump the xStar and yStar
    call dump ( l2pc%j%col%vec, details=details, name='xStar' )
    call dump ( l2pc%j%row%vec, details=details, name='yStar' )

    ! Now dump kStar
    call dump ( l2pc%j, 'kStar', details )

    ! Now dump the Hessian
    if ( l2pc%goth ) call dump ( l2pc%h, 'hStar', details )

  end subroutine DumpOneL2PC

  ! --------------------------------------- DumpL2PCInfo ---------------
  subroutine DumpL2PCInfo ( L2PCInfo, details )
    ! This subroutine dumps an l2pc to stdout

    ! Dummy arguments
    type (L2PCInfo_T), intent(in) :: L2PCInfo
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump multidim arrays
    !                                        ! -1 Skip even 1-d arrays
    !                                        ! -2 Skip all but size
    !                                        ! >0 Dump even multi-dim arrays
    !                                        ! Default 0
    ! Local variables
    integer :: i
    ! Executable
    call outputNamedValue ( 'fileID', L2PCInfo%fileID )
    call outputNamedValue ( 'binID', L2PCInfo%binID )
    call outputNamedValue ( 'blocksID', L2PCInfo%blocksID )
    call outputNamedValue ( 'hblocksID', L2PCInfo%hBlocksID )
    call outputNamedValue ( 'matrixName', trim(L2PCInfo%matrixName) )
  end subroutine DumpL2PCInfo

  ! ---------------------------------- LoadMatrix ----------
  subroutine LoadMatrix ( matrix, name, message, IgnoreHessian )
    ! This function copies the matrix from the l2pc database into the
    ! given matrix (presumably from the main dabase)

    ! Dummy arguments
    type(Matrix_T), intent(inout) :: MATRIX
    integer, intent(in) :: NAME
    character (len=*), intent(out) :: MESSAGE
    logical, intent(in), optional :: IgnoreHessian

    ! Local variables
    integer :: I                        ! Index
    type(Matrix_T), pointer :: SOURCE

    ! Executable code
    message = ''
    i = FindFirst ( l2pcDatabase%name, name )
    if ( i == 0 ) then
      message = 'No such l2pc matrix'
      return
    end if
    call PopulateL2PCBin ( i, ignoreHessian )
    source => l2pcDatabase(i)%j
    if ( .not. DoVectorsMatch ( matrix%row%vec, source%row%vec ) ) then
      message = 'Rows do not match for loading'
      return
    endif
    if ( .not. DoVectorsMatch ( matrix%col%vec, source%col%vec ) ) then
      message = 'Columns do not match for loading'
      return
    endif
    call CopyMatrixValue ( matrix, source, allowNameMismatch=.true. )
  end subroutine LoadMatrix

  ! ---------------------------------- LoadVector ----------
  subroutine LoadVector ( vector, name, source, message, IgnoreHessian )
    ! This function copies the matrix from the l2pc database into the
    ! given matrix (presumably from the main dabase)

    ! Dummy arguments
    type(Vector_T), intent(inout) :: VECTOR
    integer, intent(in) :: NAME
    integer, intent(in) :: SOURCE       ! l_rows or l_columns
    character (len=*), intent(out) :: MESSAGE
    logical, intent(in), optional :: IgnoreHessian

    ! Local variables
    integer :: I                        ! Index
    type(Matrix_T), pointer :: MATRIX
    type(Vector_T), pointer :: SOURCEVECTOR

    ! Executable code
    message = ''
    i = FindFirst ( l2pcDatabase%name, name )
    if ( i == 0 ) then
      message = 'No such l2pc matrix'
      return
    end if
    call PopulateL2PCBin ( i, ignoreHessian )
    matrix => l2pcDatabase(i)%j
    select case ( source )
    case ( l_rows )
      sourceVector => matrix%row%vec
    case ( l_columns )
      sourceVector => matrix%col%vec
    end select
    if ( .not. DoVectorsMatch ( vector, sourceVector ) ) then
      message = 'Vector does not match for loading'
      return
    endif
    ! Use the clone option to avoid fussyness with vector templates
    call CopyVector ( vector, sourceVector, allowNameMismatch=.true. )
  end subroutine LoadVector

  ! ----------------------------------- AdoptVectorTemplate --
  type (VectorTemplate_T) function AdoptVectorTemplate ( bin, name, quantityTemplates, source, message ) &
    & result ( vectorTemplate )
    ! This function is used to import a vector template from our private l2pc database
    ! into the mainstream database.  It also adopts the individual quantities into
    ! the mainstream quantity template database.

    ! Dummy arguments
    integer, intent(in) :: BIN          ! Name of MATRIX
    integer, intent(in) :: NAME         ! Name for new vector template
    type(QuantityTemplate_T), dimension(:), pointer :: QUANTITYTEMPLATES ! Mainstream database
    integer, intent(in) :: SOURCE       ! L_COLUMNS, L_ROWS
    character (len=*), intent(out) :: MESSAGE ! Error message

    ! Local variables
    integer :: I                        ! Array index
    type(VectorTemplate_T), pointer :: SOURCETEMPLATE ! Sought vector template
    type(Matrix_T), pointer :: M        ! The matrix
    integer :: QTY                      ! Loop counter
    type(QuantityTemplate_T) :: QT      ! A template

    ! Executable code
    message = ''
    call NullifyVectorTemplate ( vectorTemplate )

    i = FindFirst ( l2pcDatabase%name, bin )
    if ( i == 0 ) then
      message = 'No such l2pc bin found'
      return
    end if

    m => l2pcDatabase(i)%j
    if ( source == l_rows ) then
      sourceTemplate => m%row%vec%template
    else
      sourceTemplate => m%col%vec%template
    end if

    vectorTemplate = sourceTemplate
    vectorTemplate%name = name
    ! We need our own quantities array to point to the mainstream database
    nullify ( vectorTemplate%quantities )
    call Allocate_test ( vectorTemplate%quantities, vectorTemplate%noQuantities, &
      & 'adopted vectorTemplate%quantities', ModuleName )
    ! Copy each quantity from the l2pc database into the mainstream one
    do qty = 1, vectorTemplate%noQuantities
      call CopyQuantityTemplate ( qt, l2pcQTs ( sourceTemplate%quantities(qty) ) )
      i = FindFirst ( quantityTemplates%name, qt%name )
      if ( i /= 0 ) then
        ! Found this quantity already. Is it marked as one to be adopted
        if ( quantityTemplates(i)%quantityType == l_adopted ) then
          vectorTemplate%quantities(qty) = i
          ! Shallow copy this into the database over the 'adopted' placeholder,
          ! Destroy old placeholder first.
          call DestroyQuantityTemplateContents ( quantityTemplates(i) )
          quantityTemplates(i) = qt
        else
          call get_string ( qt%name, message, strip=.true. )
          message = 'Quantity ' // trim(message) // ' already exist, cannot adopt'
        end if
      else
        vectorTemplate%quantities(qty) = AddQuantityTemplateToDatabase ( quantityTemplates, qt )
      end if
      call NullifyQuantityTemplate ( qt ) ! To avoid treading on one in database
    end do
  end function AdoptVectorTemplate

  ! -------------------------------------- CreateDefaultBinSelectors --
  subroutine CreateDefaultBinSelectors
    ! This routine creates two bin selectors we know we need even
    ! if the user doesn't supply them in the l2cf.  If you add more
    ! here please make sure to also add DefaultSelector_... parmeters
    ! above.
    ! Local variables
    type (BinSelector_T) :: SEL
    integer :: DUMMY
    ! Executable code
    ! Create a latitude based one.  This is the one used by default if no
    ! other is supplied
    sel%selectorType = l_latitude
    sel%molecule = 0
    sel%nameFragment = 0
    sel%exact = .false.
    sel%heightRange = 0.0_r8
    sel%cost = 1.0_r8
    dummy = AddBinSelectorToDatabase ( binSelectors, sel )

    ! Create a template for a field_azimith one for the polarLinear model
    sel%selectorType = l_fieldAzimuth
    sel%molecule = 0
    sel%nameFragment = 0
    sel%exact = .true.
    sel%heightRange = 0.0_r8
    sel%cost = 1.0
    dummy = AddBinSelectorToDatabase ( binSelectors, sel )

  end subroutine CreateDefaultBinSelectors

  ! -------------------------------------- DestroyBinSelectorDatabase
  subroutine DestroyBinSelectorDatabase 
    ! Local variables
    integer :: STATUS                   ! Flag from deallocate
    ! Executable code
    if ( .not. associated ( binSelectors ) ) return
    deallocate ( binSelectors, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//"bin selectors binSelectors" )
  end subroutine DestroyBinSelectorDatabase

  ! ----------------------------------------------- DestroyL2PC ----
  subroutine DestroyL2PC ( l2pc )
    ! Dummy arguments
    type (L2pc_t), intent(inout), target :: L2PC

    integer :: QUANTITY                 ! Loop index
    integer :: VECTOR                   ! Loop index
    integer :: JH                       ! Loop index

    type (Vector_T), pointer :: V       ! Temporary pointer

    ! Exectuable code
    do jh = 1, 2
      if ( jh == 2 .and. .not. l2pc%gotH ) cycle
      do vector = 1, 2
        if ( vector == 1 ) then
          if ( jh == 1 ) then
            v => l2pc%j%col%vec
          else
            v => l2pc%h%col%vec
          end if
        else
          if ( jh == 1 ) then
            v => l2pc%j%row%vec
          else
            v => l2pc%h%row%vec
          end if
        end if

        do quantity = 1, size(v%quantities)
          call DestroyQuantityTemplateContents (v%quantities(quantity)%template )
        end do
        call DestroyVectorInfo ( v )
      end do
    end do

    ! Destory kStar
    call DestroyMatrix ( l2pc%j )
    if ( l2pc%goth ) call DestroyHessian ( l2pc%h )

  end subroutine DestroyL2PC

  ! ------------------------------------------- DestroyL2PCDatabase ---
  subroutine DestroyL2PCDatabase

    ! Local variables
    integer :: I, Status

    if (associated(l2pcDatabase)) then
      do i = 1, size(l2pcDatabase)
        call DestroyL2PC ( l2pcDatabase(i) )
      end do
      deallocate ( l2pcDatabase, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_deallocate // "l2pcDatabase" )
    end if

    ! Also destroy the info database (i.e. close files)
    call DestroyL2PCInfoDatabase
  end subroutine DestroyL2PCDatabase

  ! ------------------------------------ FlushL2PCBins -------------
  subroutine FlushL2PCBins
    ! Local variables
    integer :: BIN                      ! Loop counter
    integer :: I,J,K                    ! Loop counters

    type ( L2pc_t ), pointer :: L2PC
    type ( MatrixElement_T), pointer :: M0

    ! Executable code
    if ( .not. associated ( l2pcDatabase ) ) return
    do bin = 1, size ( l2pcDatabase )
      l2pc => l2pcDatabase ( bin )
      do i = 1, l2pc%j%row%NB
        do j = 1, l2pc%j%col%NB
          call DestroyBlock ( l2pc%j%block(i,j) )
        end do
      end do
      if ( l2pc%goth ) then
        do i = 1, l2pc%h%row%NB
          do j = 1, l2pc%h%col%NB
            do k = 1, l2pc%h%col%NB
              call DestroyBlock ( l2pc%h%block(i,j,k) )
            end do
          end do
        end do
      end if
    end do
  end subroutine FlushL2PCBins

  ! --------------------------------------------- OutputHDF5L2PC
  subroutine OutputHDF5L2PC ( filename, matrices, hessians, quantitiesNode, packed, dontPack )
  use HDF5, only: H5FCREATE_F, H5FClose_F, H5F_ACC_TRUNC_F
    character (len=*), intent(in) :: FILENAME
    type (Matrix_Database_T), dimension(:), pointer :: MATRICES
    type (Hessian_T), dimension(:), pointer :: HESSIANS
    integer, intent(in) :: QUANTITIESNODE
    logical, intent(in) :: PACKED
    integer, dimension(:), pointer :: DONTPACK

    ! Local variables
    integer :: FILEID                   ! ID of file
    integer :: STATUS                   ! From HDF
    integer :: FIELD                    ! Node index
    integer :: DB_INDEX                 ! Index of matrix
    integer :: NXT_INDEX                ! Index of next matrix (Hessian?)
    logical :: GOTH                     ! True if there is a Hessian
    type (Matrix_T), pointer :: tmpMatrix
    type (Hessian_T), pointer :: tmpHessian

    ! Executable code
    call H5FCreate_F ( trim(filename), H5F_ACC_TRUNC_F, fileID, &
      & status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open hdf5 l2pc file for output.' )
    field = 2
    do
      if ( field > nsons(quantitiesNode) ) exit
      db_index = decoration(decoration(subtree(field, quantitiesNode )))
      call GetActualMatrixFromDatabase ( matrices(db_index), tmpMatrix )
      goth = .false.
      if ( field < nsons(quantitiesNode) ) then
        nxt_index = decoration(decoration(subtree(field+1, quantitiesNode )))
        if ( nxt_index < 0 ) then
          tmpHessian => hessians ( -nxt_index )
          goth = .true.
          field = field + 1
        end if
      end if
      if ( goth ) then
        call writeOneHDF5L2PC ( tmpMatrix, fileID, packed, dontPack, hessian=tmpHessian )
      else
        call writeOneHDF5L2PC ( tmpMatrix, fileID, packed, dontPack )
      end if
      field = field + 1
    end do ! Loop over fields
    call H5FClose_F ( fileID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName,&
      & 'Unable to close hdf5 l2pc file.' )

  end subroutine OutputHDF5L2PC

  ! --------------------------------------- WriteOneHDF5L2PC -----------
  subroutine WriteOneHDF5L2PC ( JACOBIAN, fileID, packed, dontPack, hessian )
    use HessianModule_0, only: OptimizeBlock
    use HDF5, only: H5GCLOSE_F, H5GCREATE_F
    use MLSHDF5, only: MakeHDF5Attribute, SaveAsHDF5DS
    ! This subroutine writes an l2pc to a file in hdf5 format

    ! Dummy arguments
    type (matrix_T), intent(in), target :: JACOBIAN
    integer, intent(in) :: fileID
    logical, intent(in) :: PACKED
    integer, dimension(:), pointer :: DONTPACK
    type (Hessian_T), intent(in), optional :: HESSIAN

    ! Local variables
    integer :: I, J, K                  ! Loop index (row, col, col)
    integer :: BlocksGID                ! ID of group containing all blocks
    integer :: BlockGID                 ! ID of this block group
    integer :: matrixID                 ! ID of hdf5 group containing matrix
    integer :: STATUS                   ! Error flag

    integer, dimension(jacobian%row%nb) :: ROWBLOCKMAP
    integer, dimension(jacobian%col%nb) :: COLBLOCKMAP
    logical, dimension(jacobian%row%vec%template%noQuantities), target :: ROWPACK
    logical, dimension(jacobian%col%vec%template%noQuantities), target :: COLPACK

    character ( len=32 ) :: NAME        ! A name for output

    type (MatrixElement_T), pointer :: M0 ! A Matrix_0 within kStar
    type (HessianElement_T), pointer :: H0 ! A Hessian_0 within hStar

    ! Executable code
    ! If we have a Hessian, make sure it matches the jacobian
    if ( present ( hessian ) ) then
      if ( jacobian%row%vec%name /= hessian%row%vec%name ) call MLSMessage ( &
        & MLSMSG_Error, ModuleName, 'Rows of Jacobian and Hessian do not match' )
      if ( jacobian%col%vec%name /= hessian%col%vec%name ) call MLSMessage ( &
        & MLSMSG_Error, ModuleName, 'Columns of Jacobian and Hessian do not match' )
      if ( jacobian%row%instFirst .neqv. hessian%row%instFirst ) call MLSMessage ( &
        & MLSMSG_Error, ModuleName, 'Row order of Jacobian and Hessian do not match' )
      if ( jacobian%col%instFirst .neqv. hessian%col%instFirst ) call MLSMessage ( &
        & MLSMSG_Error, ModuleName, 'Column order of Jacobian and Hessian do not match' )
    end if

    ! Work out which quantities we can skip
    if ( packed ) then
      call MakeMatrixPackMap ( jacobian, rowPack, colPack, rowBlockMap, colBlockMap, dontPack, hessian )
    else
      rowPack = .true.
      colPack = .true.
      do i = 1, jacobian%row%nb
        rowBlockMap(i) = i
      end do
      do j = 1, jacobian%col%nb
        colBlockMap(j) = j
      end do
    end if

    ! First create the group for this.
    call get_string ( jacobian%name, name, strip=.true., cap=.true. )
    call h5gCreate_f ( fileID, trim(name), matrixID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create group for l2pc matrix' )

    ! Now dump the vectors
    call WriteVectorAsHDF5L2PC ( matrixID, jacobian%col%vec, 'Columns', colPack )
    call WriteVectorAsHDF5L2PC ( matrixID, jacobian%row%vec, 'Rows', rowPack )

    ! Now create a group for the blocks
    call h5gCreate_f ( matrixID, 'Blocks', blocksGID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create group for all l2pc matrix blocks' )

    ! Make the flags
    call MakeHDF5Attribute ( blocksGID, 'rowInstanceFirst', jacobian%row%instFirst )
    call MakeHDF5Attribute ( blocksGID, 'colInstanceFirst', jacobian%col%instFirst )

    ! Now loop over the blocks and write them.
    do i = 1, jacobian%row%NB
      do j = 1, jacobian%col%NB
        ! Do we write this block?
        if ( rowPack(jacobian%row%quant(i)) .and. &
          &  colPack(jacobian%col%quant(j)) ) then
          ! Identify the block
          m0 => jacobian%block(i, j)
          ! Get a name for this group for the block
          ! write ( name, * ) 'Block', rowBlockMap(i), colBlockMap(j)
          call CreateBlockName ( rowBlockMap(i), colBlockMap(j), 0, name )
          ! Create a grop for this block
          call h5gCreate_f ( blocksGID, trim(name), blockGID, status )
          if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Unable to create group for jacobian matrix block' )
          ! Stick some attributes on it
          call get_string ( &
            & jacobian%row%vec%quantities(&
            &    jacobian%row%quant(i))%template%name, name )
          call MakeHDF5Attribute ( blockGID, 'rowQuantity', trim(name) )
          call get_string ( &
            & jacobian%col%vec%quantities(&
            &    jacobian%col%quant(j))%template%name, name )
          call MakeHDF5Attribute ( blockGID, 'colQuantity', trim(name) )
          call MakeHDF5Attribute ( blockGID, 'rowInstance', jacobian%row%inst(i) )
          call MakeHDF5Attribute ( blockGID, 'colInstance', jacobian%col%inst(j) )
          call MakeHDF5Attribute ( blockGID, 'kind', m0%kind )
          ! Write the datasets
          if ( m0%kind /= m_absent ) then
            call SaveAsHDF5DS ( blockGID, 'values', real ( m0%values, r4 ) )
            if ( m0%kind /= m_full ) then
              call MakeHDF5Attribute ( blockGID, 'noValues', size(m0%values) )
              call SaveAsHDF5DS ( blockGID, 'r1', m0%r1 )
              call SaveAsHDF5DS ( blockGID, 'r2', m0%r2 )
            end if                      ! Not sparse or banded
          end if                        ! Not absent
          ! Close group for block
          call h5gClose_f ( blockGID, status )
          if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Unable to close group for l2pc matrix block' )
        end if                          ! Do this block
      end do
    end do

    ! Now close blocks group
    call h5gClose_f ( blocksGID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close qroup for all l2pc matrix blocks' )

    ! Create a group for the Hessian blocks if any
    if ( present ( hessian ) ) then
      call h5gCreate_f ( matrixID, 'HessianBlocks', blocksGID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create group for all l2pc matrix hessian blocks' )

      ! Make the flags
      call MakeHDF5Attribute ( blocksGID, 'rowInstanceFirst', hessian%row%instFirst )
      call MakeHDF5Attribute ( blocksGID, 'colInstanceFirst', hessian%col%instFirst )

      ! Now loop over the blocks and write them.
      do i = 1, hessian%row%NB
        do j = 1, hessian%col%NB
          do k = 1, hessian%col%NB
            ! Do we write this block?
            if ( rowPack(hessian%row%quant(i)) .and. &
              &  colPack(hessian%col%quant(j)) .and. &
              &  colPack(hessian%col%quant(k)) ) then
              ! Identify the block
              h0 => hessian%block ( i, j, k )
              call optimizeBlock ( h0 )
              ! Skip the absent ones in order to make life simpler
!             if ( h0%kind == h_absent ) cycle
              ! Get a name for this group for the block
              call CreateBlockName ( rowBlockMap(i), colBlockMap(j), colBlockMap(k), name )
              ! Create a group for this block
              call h5gCreate_f ( blocksGID, trim(name), blockGID, status )
              if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
                & 'Unable to create group for hessian matrix block' )
              ! Stick some attributes on it
              ! Row name
              call get_string ( &
                & hessian%row%vec%quantities(&
                &    hessian%row%quant(i))%template%name, name )
              call MakeHDF5Attribute ( blockGID, 'rowQuantity', trim(name) )
              call MakeHDF5Attribute ( blockGID, 'rowInstance', hessian%row%inst(i) )
              ! First column name
              call get_string ( &
                & hessian%col%vec%quantities(&
                &    hessian%col%quant(j))%template%name, name )
              call MakeHDF5Attribute ( blockGID, 'col1Quantity', trim(name) )
              call MakeHDF5Attribute ( blockGID, 'col1Instance', hessian%col%inst(j) )
              ! Second column name
              call get_string ( &
                & hessian%col%vec%quantities(&
                &    hessian%col%quant(k))%template%name, name )
              call MakeHDF5Attribute ( blockGID, 'col2Quantity', trim(name) )
              call MakeHDF5Attribute ( blockGID, 'col2Instance', hessian%col%inst(j) )

              call MakeHDF5Attribute ( blockGID, 'kind', h0%kind )
              ! Write the datasets

              if ( h0%kind == h_full ) then
                call SaveAsHDF5DS ( blockGID, 'values', real ( h0%values, r4 ) )
              end if
              if ( h0%kind == h_sparse ) then
                call MakeHDF5Attribute ( blockGID, 'noValues', h0%tuplesFilled )
                call SaveAsHDF5DS ( blockGID, 'i', h0%tuples(1:h0%tuplesFilled)%i )
                call SaveAsHDF5DS ( blockGID, 'j', h0%tuples(1:h0%tuplesFilled)%j )
                call SaveAsHDF5DS ( blockGID, 'k', h0%tuples(1:h0%tuplesFilled)%k )
                call SaveAsHDF5DS ( blockGID, 'h', h0%tuples(1:h0%tuplesFilled)%h )
              end if
              ! Close group for block
              call h5gClose_f ( blockGID, status )
              if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
                & 'Unable to close group for l2pc matrix block' )
            end if                          ! Do this block
          end do
        end do
      end do

      ! Now close blocks group
      call h5gClose_f ( blocksGID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close qroup for all l2pc matrix blocks' )
    end if

    ! Close matrix group
    call h5gClose_f ( matrixID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close group for l2pc matrix' ) 
 end subroutine WriteOneHDF5L2PC

  ! ======================================= PRIVATE PROCEDURES ====================

  ! ------------------------------------ AddL2PCInfoToDatabase --
  integer function AddL2PCInfoToDatabase ( database, item )
    type (L2PCInfo_T), dimension(:), pointer :: DATABASE
    type (L2PCInfo_T) :: item
    ! Local variables
    type (L2PCInfo_T), dimension(:), pointer :: TEMPDATABASE

    include "addItemToDatabase.f9h"
    AddL2PCInfoToDatabase = newSize
  end function AddL2PCInfoToDatabase

  ! ------------------------------------ DestroyL2PCInfoDatabase ----
  subroutine DestroyL2PCInfoDatabase
  use HDF5, only: H5FClose_F, H5GCLOSE_F, H5ESET_AUTO_F
    ! Local variables
    integer :: I                ! Loop counter
    integer :: STATUS           ! Flag from HDF

    ! Executable code

    if ( .not. associated(l2pcInfo) ) return
    call h5eSet_auto_f ( 0, status )
    do i = 1, size ( l2pcInfo )
      if ( index ( switches, 'l2pc' ) > 0 ) &
        & call outputNamedValue( 'Destroying l2pc db entry number ', i )
      call h5gClose_f ( l2pcInfo(i)%blocksID, status )
      if ( status /= 0 .and. DIEIFDESTROYFAILS ) then
        call dump_private( l2pcInfo(i) )
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close Blocks group for preserved input l2pc' )
      endif
      if ( l2pcInfo(i)%hblocksID /= 0 ) then
        call h5gClose_f ( l2pcInfo(i)%hblocksID, status )
        if ( status /= 0 .and. DIEIFDESTROYFAILS ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to close Hessian Blocks group for preserved input l2pc' )
      end if
      call h5gClose_f ( l2pcInfo(i)%binID, status )
      if ( status /= 0 .and. DIEIFDESTROYFAILS ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close matrix group for preserved input l2pc' )
      ! Close the file?
      if ( count ( l2pcInfo%fileID == l2pcInfo(i)%fileID ) == 1 ) then
        ! We're the only one (left?) with this file, close it.
        call h5fClose_f ( l2pcInfo(i)%fileID, status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to close hdf5 preserved input l2pc file' )
      end if
      l2pcInfo(i)%fileID = 0
    end do
    call h5eSet_auto_f ( 1, status )
    deallocate ( l2pcInfo, stat=i )
    if ( i /= 0 .and. DIEIFDESTROYFAILS ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//'l2pcInfo' )
  end subroutine DestroyL2PCInfoDatabase

  ! ----------------------------------- MakeMatrixPackMap -----------
  subroutine MakeMatrixPackMap ( j, rowPack, colPack, rowBlockMap, colBlockMap, dontpack, h )
    ! This subroutine fills the boolean arrays rowPack, colPack
    ! (each length row/col%noQuantities) with a flag set true
    ! if the quantity has any derivatives at all

    ! Dummy arguments
    type (Matrix_t), intent(in) :: J
    logical, intent(out), dimension(J%row%vec%template%noQuantities) :: ROWPACK
    logical, intent(out), dimension(J%col%vec%template%noQuantities) :: COLPACK
    integer, intent(out), dimension(J%row%NB) :: ROWBLOCKMAP
    integer, intent(out), dimension(J%col%NB) :: COLBLOCKMAP
    integer, dimension(:), pointer :: DONTPACK ! Quantities not to pack
    type (Hessian_T), intent(in), optional :: H ! A hessian to consider also

    ! Local variables
    integer :: ROWQ                     ! Loop counter
    integer :: COLQ                     ! Loop counter
    integer :: COLQ1                    ! Loop counter
    integer :: COLQ2                    ! Loop counter
    integer :: ROWI                     ! Loop counter
    integer :: COLI                     ! Loop counter
    integer :: COLI1                    ! Loop counter
    integer :: COLI2                    ! Loop counter
    integer :: ROWBLOCK                 ! Block index
    integer :: COLBLOCK                 ! Block index
    integer :: COLBLOCK1                ! Block index
    integer :: COLBLOCK2                ! Block index
    logical, dimension(J%row%NB) :: ROWBLOCKFLAG ! Flags per row block
    logical, dimension(J%col%NB) :: COLBLOCKFLAG ! Flags per row block

    ! Executable code
    rowPack = .false.
    colPack = .false.

    ! Do a nested loop over cols/rows
    ! I tried to be fancy with cycles etc. but the code got really messy.
    ! This simple approach is probably the clearest.
    do rowQ = 1, j%row%vec%template%noQuantities
      do rowI = 1, j%row%vec%quantities(rowQ)%template%noInstances
        rowBlock = FindBlock ( j%row, rowQ, rowI )
        do colQ = 1, j%col%vec%template%noQuantities
          do colI = 1, j%col%vec%quantities(colQ)%template%noInstances
            colBlock = FindBlock ( j%col, colQ, colI )
            if ( j%block ( rowBlock, colBlock ) % kind /= M_Absent ) then
              rowPack ( rowQ ) = .true.
              colPack ( colQ ) = .true.
            end if
          end do
        end do
      end do
    end do

    ! Do the same thing for any hessian supplied
    if ( present ( h ) ) then
      do rowQ = 1, h%row%vec%template%noQuantities
        do rowI = 1, h%row%vec%quantities(rowQ)%template%noInstances
          rowBlock = FindBlock ( h%row, rowQ, rowI )
          do colQ1 = 1, h%col%vec%template%noQuantities
            do colI1 = 1, h%col%vec%quantities(colQ1)%template%noInstances
              colBlock1 = FindBlock ( h%col, colQ1, colI1 )
              do colQ2 = 1, h%col%vec%template%noQuantities
                do colI2 = 1, h%col%vec%quantities(colQ2)%template%noInstances
                  colBlock2 = FindBlock ( h%col, colQ2, colI2 )
                  if ( h%block ( rowBlock, colBlock1, colBlock2 ) % kind /= H_Absent ) then
                    rowPack ( rowQ ) = .true.
                    colPack ( colQ1 ) = .true.
                    colPack ( colQ2 ) = .true.
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
    end if

    ! Go back and set to true any we know we want to keep
    if ( associated ( dontPack ) ) then
      do colQ = 1, j%col%vec%template%noQuantities
        if ( any ( dontPack == j%col%vec%template%quantities(colQ) ) ) &
          & colPack ( colQ ) = .true.
      end do
      do rowQ = 1, j%row%vec%template%noQuantities
        if ( any ( dontPack == j%row%vec%template%quantities(rowQ) ) ) &
          & rowPack ( rowQ ) = .true.
      end do
    end if

    ! Now work out the block mappings
    rowBlockFlag = rowPack ( j%row%quant )
    do rowBlock = 1, j%row%nb
      rowBlockMap ( rowBlock ) = count ( rowBlockFlag ( 1:rowBlock ) )
    end do
    colBlockFlag = colPack ( j%col%quant )
    do colBlock = 1, j%col%nb
      colBlockMap ( colBlock ) = count ( colBlockFlag ( 1:colBlock ) )
    end do

  end subroutine MakeMatrixPackMap

  ! --------------------------------------- PopulateL2PCBinByName ---
  subroutine PopulateL2PCBinByName ( name, IgnoreHessian )
    integer, intent(in) :: NAME         ! Name index
    logical, intent(in), optional :: IgnoreHessian
    integer :: INDEX
    ! Executable code
    index = FindFirst ( l2pcDatabase%name, name )
    if ( index /= 0 ) call PopulateL2PCBin ( index, ignoreHessian )
  end subroutine PopulateL2PCBinByName

  ! --------------------------------------- Populate L2PCBin --------
  subroutine PopulateL2PCBin ( bin, IgnoreHessian )
    use HDF5, only: H5GCLOSE_F, H5GOPEN_F, H5GGET_OBJ_INFO_IDX_F, &
        & H5GN_MEMBERS_F
    use MLSHDF5, only: GetHDF5Attribute, LoadFromHDF5DS

    integer, intent(in) :: BIN ! The bin index to populate
    logical, intent(in), optional :: IgnoreHessian

    ! Local variables
    type(L2pc_t), pointer :: L2PC  ! This l2pc
    type(L2PCInfo_T), pointer :: INFO ! Info for this l2pc
    type(MatrixElement_T), pointer :: M0 
    integer :: BLOCKROW  ! Loop counter
    integer :: BLOCKCOL  ! Loop counter
    character ( len=64 ) :: NAME ! Name of a block
    integer :: BLOCKID   ! Group ID for a block
    integer :: BLOCKSID  ! Group ID for all the non-Hessian blocks
    integer :: KIND      ! Kind for this block
    logical :: NoHessian ! Don't try to get the Hessian
    integer :: NOVALUES  ! Number of values for this block
    integer :: STATUS    ! Flag from HDF5
    integer :: i, nmembers, objtype
    ! Executable code

    l2pc => l2pcDatabase ( bin )
    if ( .not. any ( l2pc%j%block%kind == m_unknown ) ) return
    info => l2pcInfo ( bin )

    call h5gn_members_f( info%binID, 'Blocks', nmembers, status )
    if ( index ( switches, 'l2pc' ) > 0 ) then
      call outputNamedValue ( 'nmembers under ' // 'Blocks', nmembers )
      do i=0, nmembers-1
        call h5gget_obj_info_idx_f( info%binID, 'Blocks', i, &
                                           name, objtype, status )
        call outputNamedValue ( 'member name: ', trim(name) )
      enddo
    endif

    call h5gOpen_f ( info%binId, 'Blocks', blocksId, status )
    if ( status /= 0 ) then
      call outputNamedValue( 'binID', info%binID )
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open Blocks group for l2pc matrix block '//trim(info%matrixName) )
    endif
    do blockRow = 1, l2pc%j%row%NB
      do blockCol = 1, l2pc%j%col%NB
        ! Skip blocks we know about or are absent
        m0 => l2pc%j%block ( blockRow, blockCol )
        if ( m0%kind /= m_unknown ) cycle
        ! Access this block
        ! write ( name, * ) 'Block', blockRow, blockCol
        call CreateBlockName ( blockRow, blockCol, 0, name )
        ! This doesn't work any more (Why not?)
        ! call h5gOpen_f ( info%blocksId, trim(name), blockId, status )
        info%blocksId = blocksID ! When did these diverge?
        call h5gOpen_f ( blocksId, trim(name), blockId, status )
        if ( status /= 0 ) then
          call outputNamedValue( 'Block name', trim(name) )
          call outputNamedValue( 'Blocks ID', info%blocksId )
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Unable to open group for l2pc matrix block '//trim(name) )
        elseif ( switchDetail( switches, 'l2pc' ) > 0 ) then
          call outputNamedValue( 'Block name', trim(name) )
          call outputNamedValue( 'Blocks ID', info%blocksId )
          call outputNamedValue( 'Block ID', blockId )
        endif
        ! Get kind of block
        call GetHDF5Attribute ( blockID, 'kind', kind )

        ! Read the block
        if ( kind == m_banded .or. kind == m_column_sparse ) then
          call GetHDF5Attribute ( blockID, 'noValues', noValues )
          call CreateBlock ( l2pc%j, blockRow, blockCol, kind, noValues )
          m0 => l2pc%j%block ( blockRow, blockCol )
          call LoadFromHDF5DS ( blockId, 'r1', m0%r1 )
          call LoadFromHDF5DS ( blockId, 'r2', m0%r2 )
        else
          call CreateBlock ( l2pc%j, blockRow, blockCol, kind )
          m0 => l2pc%j%block ( blockRow, blockCol )
        end if
        if ( kind /= m_absent ) &
          call LoadFromHDF5DS ( blockID, 'values', m0%values )
        call h5gClose_f ( blockId, status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to close group for l2pc matrix block '//trim(name) )
      end do
    end do
    call h5gClose_f ( blocksId, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close Blocks group for l2pc matrix '//trim(info%matrixName) )

    noHessian = .false. ! Assume we want the Hessian
    if ( present(ignoreHessian) ) noHessian = ignoreHessian
    if ( l2pc%goth .and. .not. noHessian ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to fully populate L2PC Bin when we have a Hessian matrix' )
      !!! stop !!! CODE HERE
    endif

    if ( index ( switches, 'spa' ) /= 0 ) call dump_struct ( l2pc%j, 'Populated l2pc bin' )

  end subroutine PopulateL2PCBin

  ! --------------------------------------- ReadCompleteHDF5L2PC -------
  subroutine ReadCompleteHDF5L2PCFile ( MLSFile, Where )
    use HDF5, only: H5F_ACC_RDONLY_F, h5fopen_f, H5GN_MEMBERS_F
    use Trace_M, only: Trace_begin, Trace_end
    use Toggles, only: Toggle, gen
    type (MLSFile_T), pointer   :: MLSFile
    integer, intent(in) :: Where ! In the L2CF tree, for tracing
    ! character (len=*), intent(in) :: FILENAME

    ! Local variables
    integer :: FILEID          ! From hdf5
    type (L2PCInfo_T) :: INFO  ! Info for one bin
    type (L2pc_t) :: L2PC    ! The l2pc read from one bin
    integer :: STATUS          ! Flag from HDF5
    integer :: NOBINS          ! Number of bins
    integer :: BIN             ! Loop counter
    integer :: DUMMY           ! Ignored return from AddToDatabase

    ! Executable code

    if ( toggle (gen) ) call trace_begin ( "ReadCompleteHDF5L2PC", where )
    call h5fopen_f ( MLSFile%name, H5F_ACC_RDONLY_F, fileID, status )
    MLSFile%FileID%f_id = fileID
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open hdf5 l2pc file for input:'//trim(MLSFile%name), &
      & MLSFile=MLSFile )

    ! Get the number of bins
    call h5gn_members_f ( fileID, '/', noBins, status ) 
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get number of bins from input l2pc file:'//trim(MLSFile%name), &
      & MLSFile=MLSFile )

    if ( index ( switches, 'l2pc' ) /= 0 ) then
      call output ( 'Reading l2pc ' )
      call output ( trim(MLSFile%name), advance='yes' )
      call output ( 'Number of bins: ' )
      call output ( noBins, advance='yes' )
    endif

    ! Don't forget HDF5 numbers things from zero
    do bin = 0, noBins-1
      call ReadOneHDF5L2PCRecord ( L2PC, MLSFile, bin, &
        & shallow=.true., info=Info )
      dummy = AddL2PCToDatabase ( l2pcDatabase, L2PC )
      dummy = AddFileIDToDatabase ( fileIDDatabase, MLSFile%fileID%f_id )
      if ( index ( switches, 'spa' ) /= 0 ) call Dump_struct ( l2pc%j, 'One l2pc bin' ) 

      ! Now nullify the pointers in l2pc so we don't clobber the one we've written
      nullify ( l2pc%j%block )
      nullify ( l2pc%j%row%nelts, l2pc%j%row%inst, l2pc%j%row%quant )
      nullify ( l2pc%j%col%nelts, l2pc%j%col%inst, l2pc%j%col%quant )
      nullify ( l2pc%j%row%vec%template%quantities, l2pc%j%col%vec%template%quantities )
      nullify ( l2pc%j%row%vec%quantities, l2pc%j%col%vec%quantities )

      nullify ( l2pc%h%block )
      nullify ( l2pc%h%row%nelts, l2pc%h%row%inst, l2pc%h%row%quant )
      nullify ( l2pc%h%col%nelts, l2pc%h%col%inst, l2pc%h%col%quant )
      nullify ( l2pc%h%row%vec%template%quantities, l2pc%h%col%vec%template%quantities )
      nullify ( l2pc%h%row%vec%quantities, l2pc%h%col%vec%quantities )
      dummy = AddL2PCInfoToDatabase ( l2pcInfo, Info )
    end do

    ! Don't close the file, we're keeping it open to read blocks from it later
    if ( toggle (gen) ) call trace_end ( "ReadCompleteHDF5L2PC" )
  end subroutine ReadCompleteHDF5L2PCFile

  ! --------------------------------------- ReadOneHDF5L2PCRecord ------------
  subroutine ReadOneHDF5L2PCRecord ( l2pc, MLSFile, l2pcIndex, shallow, info )
    use HDF5, only: H5GCLOSE_F, H5GOPEN_F, H5GGET_OBJ_INFO_IDX_F, &
      & H5GN_MEMBERS_F
    use MLSHDF5, only: GetHDF5Attribute, LoadFromHDF5DS
    type ( L2pc_t ), intent(out), target :: L2PC
    type (MLSFile_T), pointer   :: MLSFile
    ! integer, intent(in) :: FILEID       ! HDF5 ID of input file
    integer, intent(in) :: L2PCINDEX        ! Index of l2pc entry to read
    logical, optional, intent(in) :: SHALLOW ! Don't read blocks
    type(L2PCInfo_T), intent(out), optional :: INFO ! Information output

    ! Local variables
    integer :: BLOCKID                  ! Id of block group (one block at a time)
    integer :: BLOCKSID                 ! Id of blocks group
    integer :: HBLOCKSID                ! Id of HessianBlocks group
    integer :: I, J, K                  ! Block loop indices
    integer :: KIND                     ! Kind of block (absent, etc.)
    integer :: MATRIXID                 ! HDF5 for matrix group
    integer :: NMEMBERS
    integer :: NOVALUES                 ! For banded or column sparse cases
    integer :: OBJTYPE                  ! From HDF5
    integer :: STATUS                   ! Flag from HDF5
    integer :: STRINGINDEX              ! Index of string
    integer :: XSTAR                    ! Linearisation state vector index
    integer :: YSTAR                    ! Radiances for xStar vector index

    logical :: COLINSTANCEFIRST         ! Flag for matrix
    logical :: ROWINSTANCEFIRST         ! Flag for matrix
    logical :: MYSHALLOW                ! Value of shallow

    type (MatrixElement_T), pointer :: M0 ! A Matrix block within kStar
    type (HessianElement_T), pointer :: H0 ! A Hessian block within kStar
    character ( len=64 ) :: MATRIXNAME  ! Name for matrix
    character ( len=64 ) :: NAME        ! Name for block group


    ! Executable code
    myShallow = .false.
    if ( present ( shallow ) ) myShallow = shallow

    if ( index ( switches, 'l2pc' ) /= 0 ) &
      & call output ( 'Reading bin from l2pc file', advance='yes' )

    call h5gGet_obj_info_idx_f ( MLSFile%fileID%f_id, '/', l2pcIndex, matrixName, &
      & objType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get information on matrix in input l2pc file', &
      & MLSFile=MLSFile )
    if ( index ( switches, 'l2pc' ) /= 0 ) &
      & call outputNamedValue ( 'Reading matrix', matrixName )
    call h5gOpen_f ( MLSFile%fileID%f_id, trim(matrixName), matrixId, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open matrix in input l2pc file' , &
      & MLSFile=MLSFile )

    ! Read the row and column vectors
    MLSFile%fileID%grp_id = matrixId
    call ReadOneVectorFromHDF5 ( MLSFile, 'Columns', xStar )
    call ReadOneVectorFromHDF5 ( MLSFile, 'Rows', yStar )

    ! Get the instance first information
    call h5gOpen_f ( matrixID, 'Blocks', blocksID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to access Blocks group of input l2pc matrix' , &
      & MLSFile=MLSFile )
    MLSFile%fileID%sd_id = blocksID
    call GetHDF5Attribute ( MLSFile, 'rowInstanceFirst', rowInstanceFirst )
    call GetHDF5Attribute ( MLSFile, 'colInstanceFirst', colInstanceFirst )

    stringIndex = GetStringIndexFromString ( "'"//trim(matrixName)//"'" )

    ! Create the matrix
    call CreateEmptyMatrix ( l2pc%j, stringIndex, l2pcVs(yStar), l2pcVs(xStar), &
      & row_quan_first = .not. rowInstanceFirst,&
      & col_quan_first = .not. colInstanceFirst )

    ! Must *not* forget to name l2pc 
    ! (l2pc%name will itself never be used beyond this point)
    l2pc%name = l2pc%j%name

    ! Fill up the information
    if ( present ( info ) ) then
      info%fileID = MLSFile%fileID%f_id
      info%binID = matrixID
      info%blocksID = blocksID
      info%matrixName = matrixName
    end if

    if ( switchDetail( switches, 'l2pc' ) > 0 ) then
      do i = 1, l2pc%j%row%NB
        do j = 1, l2pc%j%col%NB
          call CreateBlockName ( i, j, 0, name )
          call h5gOpen_f ( blocksId, trim(name), blockId, status )
          if ( status /= 0 ) then
            call outputNamedValue( 'Block name', trim(name) )
            call outputNamedValue( 'Blocks ID', blocksId )
            call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Unable to open group for l2pc matrix block '//trim(name) , &
              & MLSFile=MLSFile )
          elseif ( switchDetail( switches, 'l2pc' ) > 0 ) then
            call outputNamedValue( 'Block name', trim(name) )
            call outputNamedValue( 'Blocks ID', blocksId )
            call outputNamedValue( 'Block ID', blockId )
          endif
          call h5gClose_f ( blockId, status )
        enddo
      enddo
    endif
    ! Loop over blocks and read them
    if ( .not. myShallow ) then
      do i = 1, l2pc%j%row%NB
        do j = 1, l2pc%j%col%NB
          ! Access this block
          call CreateBlockName ( i, j, 0, name )
          call h5gOpen_f ( blocksId, trim(name), blockId, status )
          if ( status /= 0 ) then
            call outputNamedValue( 'Block name', trim(name) )
            call outputNamedValue( 'Blocks ID', blocksId )
            call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Unable to open group for l2pc matrix block '//trim(name) , &
              & MLSFile=MLSFile )
          elseif ( switchDetail( switches, 'l2pc' ) > 0 ) then
            call outputNamedValue( 'Block name', trim(name) )
            call outputNamedValue( 'Blocks ID', blocksId )
            call outputNamedValue( 'Block ID', blockId )
          endif
          ! Could check it's the block we're expecting but I think I'll be lazy
          MLSFile%fileID%sd_id = blockID
          call GetHDF5Attribute ( MLSFile, 'kind', kind )
          if ( kind == m_banded .or. kind == m_column_sparse ) then
            call GetHDF5Attribute ( MLSFile, 'noValues', noValues )
            call CreateBlock ( l2pc%j, i, j, kind, noValues )
            m0 => l2pc%j%block ( i, j )
            call LoadFromHDF5DS ( MLSFile, 'r1', m0%r1 )
            call LoadFromHDF5DS ( MLSFile, 'r2', m0%r2 )
          else
            call CreateBlock ( l2pc%j, i, j, kind )
            m0 => l2pc%j%block ( i, j )
          end if
          if ( kind /= m_absent ) &
            call LoadFromHDF5DS ( MLSFile, 'values', m0%values )
          call h5gClose_f ( blockId, status )
          if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Unable to close group for input l2pc matrix block '//trim(name) , &
            & MLSFile=MLSFile )
        end do
      end do
    else
      ! Otherwise, flag the whole matrix as unknown
      do i = 1, l2pc%j%row%NB
        do j = 1, l2pc%j%col%NB
          call CreateBlock ( l2pc%j, i, j, m_unknown )
        end do
      end do
    end if

    call h5gClose_f ( blocksId, status )
    ! Look for any Hessian blocks
    ! We begin by assuming there are none
    call h5gn_members_f( MLSFile%fileID%f_id, trim(matrixName), nmembers, status )
    if ( index ( switches, 'l2pc' ) > 0 ) then
      call outputNamedValue ( 'nmembers under ' // trim(matrixName), nmembers )
      do i=0, nmembers-1
        call h5gget_obj_info_idx_f( MLSFile%fileID%f_id, trim(matrixName), i, &
                                           name, objtype, status )
        call outputNamedValue ( 'member name: ', trim(name) )
      enddo
    endif
    status = 1
    if ( nmembers > 3 ) call h5gOpen_f ( matrixID, 'HessianBlocks', hBlocksID, status )
    l2pc%goth = status == 0 .and. nmembers > 3 ! Got Hessian?
    if ( l2pc%goth ) then
      MLSFile%fileID%sd_id = blocksID

      l2pc%h = CreateEmptyHessian ( stringIndex, l2pcVs(yStar), l2pcVs(xStar) )
      if ( present ( info ) ) info%hBlocksID = hBlocksID

      ! Loop over blocks and read them
      if ( .not. myShallow ) then
        do i = 1, l2pc%h%row%NB
          do j = 1, l2pc%h%col%NB
            do k = 1, l2pc%h%col%NB
              ! Access this block
              call CreateBlockName ( i, j, k, name )
              call h5gOpen_f ( blocksId, trim(name), blockId, status )
              if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
                & 'Unable to open group for l2pc matrix block '//trim(name) , &
                & MLSFile=MLSFile )
              ! Could check it's the block we're expecting but I think I'll be lazy
              MLSFile%fileID%sd_id = blockID
              call GetHDF5Attribute ( MLSFile, 'kind', kind )

              select case ( kind )
              case ( h_absent )
                call CreateBlock ( l2pc%h, i, j, k, kind, noValues )
              case ( h_sparse )
                call GetHDF5Attribute ( MLSFile, 'noValues', noValues )
                call CreateBlock ( l2pc%h, i, j, k, kind, noValues )
                h0 => l2pc%h%block ( i, j, k )
                call LoadFromHDF5DS ( MLSFile, 'i', h0%tuples%i )
                call LoadFromHDF5DS ( MLSFile, 'j', h0%tuples%j )
                call LoadFromHDF5DS ( MLSFile, 'k', h0%tuples%k )
                call LoadFromHDF5DS ( MLSFile, 'h', h0%tuples%h )
              case ( h_full )
                call CreateBlock ( l2pc%h, i, j, k, kind )
                h0 => l2pc%h%block ( i, j, k )                
                call LoadFromHDF5DS ( MLSFile, 'values', h0%values )
              end select
              call h5gClose_f ( blockId, status )
              if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
                & 'Unable to close group for input l2pc matrix block '//trim(name) , &
                & MLSFile=MLSFile )
            end do
          end do
        end do
      else
        ! Otherwise, flag the whole matrix as unknown
        do i = 1, l2pc%h%row%NB
          do j = 1, l2pc%h%col%NB
            do k = 1, l2pc%h%col%NB
              call CreateBlock ( l2pc%h, i, j, k, h_unknown )
            end do
          end do
        end do
      end if
    end if                              ! Looking for Hessians

    ! finish up, though if in shallow mode, then keep the groups open
    if ( .not. myShallow ) then
      call h5gClose_f ( blocksID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close Blocks group for input l2pc' , &
        & MLSFile=MLSFile )
      call h5gClose_f ( matrixID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close matrix group for input l2pc' , &
        & MLSFile=MLSFile )
    endif
    i = switchDetail( switches, 'l2pc' )
    if ( i >= 0 ) then
      call dump ( l2pc, details=i-1 )
      call output ( 'Done reading bin from l2pc file', advance='yes' )
    end if

  end subroutine ReadOneHDF5L2PCRecord

  ! --------------------------------------- ReadOneVectorFromHDF5 ------
  subroutine ReadOneVectorFromHDF5 ( MLSFile, name, vector )
    use HDF5, only: H5GCLOSE_F, H5GOPEN_F, H5GGET_OBJ_INFO_IDX_F
    use MLSHDF5, only: GetHDF5Attribute, IsHDF5AttributePresent, &
      & IsHDF5DSPresent, LoadFromHDF5DS
    use MLSSignals_m, only: Radiometers
    ! Read a vector from an l2pc HDF5 and adds it to internal databases.
    ! Dummy arguments
    type (MLSFile_T) :: MLSFile
    character (len=*), intent(in) :: NAME ! Name of vector
    integer, intent(out) :: VECTOR      ! Index of vector read in L2PCVectors

    ! Local variables
    integer :: FREQUENCYCOORDINATE      ! Enumeration
    integer :: I                        ! Index into various arrays
    integer :: MOLECULE                 ! Enumeration
    integer :: NAMEINDEX                ! Quantity name
    integer :: NOCHANS                  ! Dimension
    integer :: NOINSTANCES              ! Dimension
    integer :: NOINSTANCESOR1           ! Dimension
    integer :: NOQUANTITIES             ! Number of quantities in vector
    integer :: NOSURFS                  ! Dimension
    integer :: NOSURFSOR1               ! Dimension
    integer :: OBJTYPE                  ! Irrelevant argument to HDF5
    integer :: QID                      ! HDF5 ID of quantity group
    integer :: QT0                      ! Offset into l2pc quantity templates database
    integer :: QTINDEXOFFSET            ! First free index in quantity template database
    integer :: QUANTITY                 ! Loop inductor
    integer :: QUANTITYTYPE             ! Enumerated
    integer :: RADIOMETER               ! Enumerated
    integer :: SIDEBAND                 ! Sideband -1,0,1
    integer :: SIGNAL                   ! Index
    integer :: STATUS                   ! Flag from HDF
    integer :: STRINGINDEX              ! A string index
    integer :: UNIT                     ! Units for quantity
    integer :: VERTICALCOORDINATE       ! Enumeratire
    integer :: VID                      ! HDF5 ID of vector group
    integer :: VTINDEX                  ! Index of vector template
    character (len=64) :: THISNAME      ! Name of this quantity
    character (len=64) :: WORD          ! General string

    logical :: COHERENT                 ! Flag
    logical :: STACKED                  ! Flag
    logical :: LOGBASIS                 ! Flag

    integer, dimension(:), pointer :: SIGINDS  ! Index into signals database
    integer, dimension(:), pointer :: QTINDS   ! Quantity indices

    character (len=64), pointer, dimension(:) :: QUANTITYNAMES ! Names of quantities
    type ( QuantityTemplate_T), pointer :: QT    ! Template for the quantity
    type ( VectorTemplate_T) :: VT    ! Template for the vector
    type ( Vector_T ) :: V              ! The vector

    ! Executable code
    nullify ( sigInds, qtInds )

    if ( index ( switches, 'l2pc' ) /= 0 ) then
      call output ( 'Reading ' )
      call output ( trim(name) )
      call output ( ' vector from l2pc', advance='yes' )
    end if
    ! Open the vector group
    call h5gOpen_f ( MLSFile%fileID%grp_id, name, vId, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open group for vector '//trim(name), &
      & MLSFile=MLSFile )

    ! Get the number of quantities
    MLSFile%fileID%sd_id = vID
    call GetHDF5Attribute ( MLSFile, 'noQuantities', noQuantities )
    call allocate_test ( qtInds, noQuantities, 'qtInds', ModuleName )
    ! Work out the order of the quantities
    nullify ( quantityNames )
    call Allocate_test ( quantityNames, noQuantities, 'quantityNames', ModuleName )
    do quantity = 1, noQuantities
      if ( index ( switches, 'l2pc' ) /= 0) &
        & call output ( quantity, before='Identifying quantity ', advance='yes' )
      call h5gget_obj_info_idx_f ( MLSFile%fileID%grp_id, name, quantity-1, thisName, &
        & objType, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to access a quantity within '//trim(name) , &
      & MLSFile=MLSFile )
      call h5gOpen_f ( vId, trim(thisName), qId, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open quantity '//trim(thisName)//' in vector '//trim(name), &
        & MLSFile=MLSFile )
      MLSFile%fileID%sd_id = qID
      call GetHDF5Attribute ( MLSFile, 'index', i )
      quantityNames ( i ) = thisName
      call h5gClose_f ( qId, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close quantity '//trim(thisName)//' in vector '//trim(name), &
        & MLSFile=MLSFile )
    end do

    qtIndexOffset = InflateQuantityTemplateDatabase ( l2pcQTs, noQuantities )
    qt0 = size ( l2pcQTs ) - noQuantities

    ! Now go through quantities in order
    do quantity = 1, noQuantities
      qt => l2pcQTs ( qtIndexOffset + quantity - 1 )
      if ( index ( switches, 'l2pc' ) /= 0 ) then
        call output ( quantity, before='Reading quantity ' )
        call output ( ': ' )
        call output ( trim(quantityNames(quantity)), advance='yes' )
      end if

      ! Point to this quanity
      call h5gOpen_f ( vId, trim(quantityNames(quantity)), qId, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open quantity '//trim(quantityNames(quantity))// &
        & ' in vector '//trim(name) , &
        & MLSFile=MLSFile )

      ! Store the name
      stringIndex = GetStringIndexFromString ( quantityNames(quantity) )
      nameIndex = stringIndex

      ! Get the quantity type
      MLSFile%fileID%sd_id = qID
      call GetHDF5Attribute ( MLSFile, 'type', word )
      quantityType = GetLitIndexFromString ( word )

      ! Get other info as appropriate
      signal = 0
      sideband = 0
      molecule = 0
      radiometer = 0
      frequencyCoordinate = l_none
      unit = phyq_dimensionless
      select case ( quantityType )
      case ( l_vmr )
        call GetHDF5Attribute ( MLSFile, 'molecule', word )
        molecule = GetLitIndexFromString ( word )
        if ( molecule == l_extinction ) then
          call GetHDF5Attribute ( MLSFile, 'radiometer', word )
          stringIndex = GetStringIndexFromString ( word )
          radiometer = FindFirst ( Radiometers%prefix, stringIndex )
          frequencyCoordinate = l_intermediateFrequency
        else
          frequencyCoordinate = l_none
        end if
        unit = phyq_vmr
      case ( l_radiance, l_TScat )
        call GetHDF5Attribute ( MLSFile, 'signal', word )
        if ( quantityType == l_radiance ) then
          call Parse_Signal ( word, sigInds, sideband=sideband )
        else
          nullify ( qt%channels )
          call Parse_Signal ( word, sigInds, sideband=sideband, channels=qt%channels )
        end if
        signal = sigInds(1)
        verticalCoordinate = l_geodAltitude
        frequencyCoordinate = l_channel
        call deallocate_test(sigInds,'sigInds',ModuleName)
        unit = phyq_temperature
      case ( l_temperature )
        unit = phyq_temperature
      case default
      end select

      ! Now read the dimensions for the quantity
      call GetHDF5Attribute ( MLSFile, 'noChans', noChans )
      call GetHDF5Attribute ( MLSFile, 'noInstances', noInstances )
      call GetHDF5Attribute ( MLSFile, 'noSurfs', noSurfs )
      call GetHDF5Attribute ( MLSFile, 'verticalCoordinate', word )
      verticalCoordinate = GetLitIndexFromString ( word )
      ! Look for frequency coordinate (optional for backwards compatability with
      ! older l2pc files.
      if ( IsHDF5AttributePresent ( qID, 'frequencyCoordinate' ) ) then
        call GetHDF5Attribute ( MLSFile, 'frequencyCoordinate', word )
        frequencyCoordinate = GetLitIndexFromString ( trim(word) )
      end if
      call GetHDF5Attribute ( MLSFile, 'logBasis', logBasis )
      call GetHDF5Attribute ( MLSFile, 'coherent', coherent )
      call GetHDF5Attribute ( MLSFile, 'stacked', stacked )

      qt%name = nameIndex
      call SetupNewQuantityTemplate ( qt, &
        & noInstances=noInstances, noSurfs=noSurfs, noChans=noChans, &
        & coherent=coherent, stacked=stacked )

      qt%quantityType = quantityType
      qt%molecule = molecule
      qt%radiometer = radiometer
      qt%signal = signal
      qt%sideband = sideband
      qt%verticalCoordinate = verticalCoordinate
      qt%frequencyCoordinate = frequencyCoordinate
      qt%unit = unit

      if (qt%coherent) then
        noInstancesOr1 = 1
      else
        noInstancesOr1 = qt%noInstances
      endif
      if (qt%stacked) then
        noSurfsOr1 = 1
      else
        noSurfsOr1 = qt%noSurfs
      endif

      ! Get the surfaces and phis
      call LoadFromHDF5DS ( MLSFile, 'surfs', qt%surfs )
      call LoadFromHDF5DS ( MLSFile, 'phi', qt%phi )
      if ( IsHDF5DSPresent ( qId, 'solarZenith' ) ) then
        call LoadFromHDF5DS ( MLSFile, 'solarZenith', qt%solarZenith )
      else
        qt%solarZenith = 0.0
      endif
      ! Try to get the frequencies
      if ( IsHDF5DSPresent ( qId, 'frequencies' ) ) then
        call Allocate_test( qt%frequencies, qt%noChans, &
          & 'qt%frequencies', ModuleName )
        call LoadFromHDF5DS ( MLSFile, 'frequencies', qt%frequencies )
      end if

      ! Now record the index for this quantity template
      qtInds ( quantity ) = quantity + qt0

      ! For the moment, close the quantity. We'll come back to it later to fill
      ! up the values for the vector
      call h5gClose_f ( qId, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close quantity '//trim(thisName)//' in vector '//trim(name), &
        & MLSFile=MLSFile )
    end do

    ! Now create a vector template with these quantities
    call ConstructVectorTemplate ( 0, l2pcQTs, qtInds, vt, forWhom=moduleName )
    vtIndex = AddVectorTemplateToDatabase ( l2pcVTs, vt )

    call deallocate_test ( qtInds, 'qtInds', ModuleName )

    ! Now create a vector for this vector template
    if ( index ( switches, 'l2pc' ) /= 0 ) &
      & call output ( 'Creating vector', advance='yes' )
    v = CreateVector ( 0, l2pcVTs(vtIndex), l2pcQTs, vectorNameText='_v_'//name )
    if ( index ( switches, 'l2pc' ) /= 0 ) &
      & call output ( 'Adding vector to database', advance='yes' )
    vector = AddVectorToDatabase ( l2pcVs, v )

    ! Now go through the quantities again and read the values
    do quantity = 1, noQuantities
      if ( index ( switches, 'l2pc' ) /= 0 ) then
        call output ( 'Reading values for ' )
        call output ( trim(quantityNames(quantity)), advance='yes' )
      end if
      call h5gOpen_f ( vId, trim(quantityNames(quantity)), qId, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open quantity '//trim(quantityNames(quantity))// &
        & ' in vector '//trim(name), &
        & MLSFile=MLSFile )
      MLSFile%fileID%sd_id = qID
      call LoadFromHDF5DS ( MLSFile, 'values', &
        & l2pcVs(vector)%quantities(quantity)%values )
      call h5gClose_f ( qId, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close quantity '//trim(thisName)//' in vector '//trim(name), &
        & MLSFile=MLSFile )
    end do

    call Deallocate_test ( quantityNames, 'quantityNames', ModuleName )

  end subroutine ReadOneVectorFromHDF5

  ! --------------------------------------- WriteVectorAsHDF5L2PC ----------
  subroutine WriteVectorAsHDF5L2PC ( location, vector, name, packInfo )
  use HDF5, only: H5GCLOSE_F, H5GCREATE_F
  use MLSHDF5, only: MakeHDF5Attribute, SaveAsHDF5DS
    use MLSSignals_m, only: Radiometers
    integer, intent(in) :: LOCATION     ! The HDF5 location for the vector
    type (Vector_T), intent(in), target :: VECTOR ! The vector to write
    character(len=*), intent(in) :: NAME ! Name for vector
    logical, dimension(:), intent(in) :: PACKINFO ! Pack vector if set

    ! Local variables
    character ( len=132 ) :: LINE       ! A line of text
    character ( len=32 ) :: QNAME       ! Name of quantity
    integer :: QID                      ! HDF5 ID for quantity group
    integer :: QINDEX                   ! Quantity index
    integer :: QUANTITY                 ! Index
    integer :: STATUS                   ! Flag from HDF5
    integer :: VID                      ! HDF5 ID for vector group
    type (QuantityTemplate_T), pointer :: QT ! The template

    ! Executable code
    ! Create a group for the vector
    call h5gCreate_f ( location, trim(name), vID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create group for vector ' // trim(name) )
    call MakeHDF5Attribute ( vID, 'noQuantities', count(packInfo) )
    qIndex = 0
    do quantity = 1, size ( vector%quantities )
      if ( packInfo(quantity) ) then
        qIndex = qIndex + 1
        qt => vector%quantities(quantity)%template
        ! Create a group for the quantity
        call get_string ( qt%name, qName )
        call h5gCreate_f ( vID, trim(qName), qID, status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to create group for quantity ' // trim(qName) )

        ! Write quantity type and index
        call MakeHDF5Attribute ( qID, 'index', qIndex )
        call get_string ( lit_indices(qt%quantityType), line )
        call MakeHDF5Attribute ( qID, 'type', trim(line) )
        ! Write other info as appropriate
        select case ( qt%quantityType )
        case (l_vmr)
          call get_string ( lit_indices(qt%molecule), line )
          call MakeHDF5Attribute ( qID, 'molecule', trim(line) )
          if ( qt%molecule == l_extinction ) then
            call get_string ( radiometers(qt%radiometer)%prefix, line )
            call MakeHDF5Attribute ( qID, 'radiometer', trim(line) )
          end if
        case (l_radiance)
          call GetSignalName ( qt%signal, line, sideband=qt%sideband )
          call MakeHDF5Attribute ( qID, 'signal', trim(line) )
        case (l_TScat)
          call GetSignalName ( qt%signal, line, sideband=qt%sideband, &
            & otherChannels=qt%channels )
          call MakeHDF5Attribute ( qID, 'signal', trim(line) )
        end select
        ! Write out the dimensions for the quantity
        call MakeHDF5Attribute ( qID, 'noChans', qt%noChans )
        call MakeHDF5Attribute ( qID, 'noInstances', qt%noInstances )
        call MakeHDF5Attribute ( qID, 'noSurfs', qt%noSurfs )
        call get_string ( lit_indices(qt%verticalCoordinate), line )
        call MakeHDF5Attribute ( qID, 'verticalCoordinate', trim(line) )
        call get_string ( lit_indices(qt%frequencyCoordinate), line )
        call MakeHDF5Attribute ( qID, 'frequencyCoordinate', trim(line) )
        call MakeHDF5Attribute ( qID, 'logBasis', qt%logBasis )
        call MakeHDF5Attribute ( qID, 'coherent', qt%coherent )
        call MakeHDF5Attribute ( qID, 'stacked', qt%stacked )
        ! Write out important coordinates
        call SaveAsHDF5DS ( qID, 'phi', qt%phi )
        call SaveAsHDF5DS ( qID, 'surfs', qt%surfs )
        call SaveAsHDF5DS ( qID, 'solarZenith', qt%solarZenith )
        if ( associated ( qt%frequencies ) ) &
          & call SaveAsHDF5DS ( qID, 'frequencies', qt%frequencies )
        ! Write out the values
        call SaveAsHDF5DS ( qID, 'values', vector%quantities(quantity)%values )
        ! Close the group
        call h5gClose_f ( qID, status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to close quantity group for ' // trim(qName) )
      end if
    end do

    ! Finish off
    call h5gClose_f ( vID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close vector group for ' // trim(name) )

  end subroutine WriteVectorAsHDF5L2PC

  subroutine CreateBlockName ( row, column1, column2, name )
    ! Because the Intel compiler produces a different result from
    ! write(str, *) int
    ! Args
    integer, intent(in)           :: row, column1, column2
    character(len=*), intent(out) :: name
    ! Internal variables
    character(len=16), dimension(3) :: strs
    ! Executable
    call writeIntsToChars( (/row, column1, column2/), strs )
    name = ' Block' // ' ' // trim(adjustl(strs(1))) // ' ' // trim(adjustl(strs(2)))
    if ( column2 > 0 ) name = trim ( name ) // ' ' // trim(adjustl(strs(3)))
  end subroutine CreateBlockName

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module L2PC_m

! $Log$
! Revision 2.95  2010/05/19 00:32:06  vsnyder
! Pass IgnoreHessian into PopulateL2PCBin
!
! Revision 2.94  2010/05/13 23:45:14  pwagner
! Temporary expedients for l2pc files with Hessians; needs more code
!
! Revision 2.93  2010/04/30 22:55:32  vsnyder
! Give more meaningful names to vectors
!
! Revision 2.92  2010/04/17 01:43:35  vsnyder
! Simplify dump.  Set l2pc%goth if there is a Hessian block.
!
! Revision 2.91  2010/03/24 20:48:55  vsnyder
! Replace 'continue' by 'cycle'.  Call OptimizeBlock from WriteOneHDF5L2PC
! when writing a Hessian, to avoid trying to write empty arrays.  Add some
! comments.  Fix a wrong-subscript bug.
!
! Revision 2.90  2010/03/19 20:18:47  pwagner
! Moved DumpL2PCInfo out of public scope to appease ifort v10 which crashed the goldbrick
!
! Revision 2.89  2010/03/17 20:59:10  pwagner
! Code around bombs in DestroyL2PCInfoDatabase
!
! Revision 2.88  2010/02/25 18:04:45  pwagner
! l2pc type can now hold matrix and Hessian types
!
! Revision 2.87  2010/02/04 23:08:00  vsnyder
! Remove USE or declaration for unused names
!
! Revision 2.86  2009/11/17 23:43:51  vsnyder
! Add ability to output TScat
!
! Revision 2.85  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.84  2008/12/18 21:31:18  pwagner
! May now dump an l2pc or allL2PCs (use with caution)
!
! Revision 2.83  2008/09/30 22:28:03  vsnyder
! Remove AuxGrids -- didn't need them after all
!
! Revision 2.82  2008/06/06 01:55:29  vsnyder
! Process AuxGrids field of vector quantities
!
! Revision 2.81  2008/05/03 01:49:55  vsnyder
! Repair a comment
!
! Revision 2.80  2007/10/03 23:58:46  vsnyder
! Add 'where' for tracing
!
! Revision 2.79  2007/04/26 20:30:32  pwagner
! Bugfix for way ifc writes ints to strings
!
! Revision 2.78  2006/08/05 02:11:58  vsnyder
! Add ForWhom argument to ConstructVectorTemplate
!
! Revision 2.77  2006/08/04 01:53:43  vsnyder
! Need to set the name, not quantityType before SetupNewQuantityTemplate
!
! Revision 2.76  2006/08/03 01:10:06  vsnyder
! Put l2cf names in leak track database
!
! Revision 2.75  2005/06/29 00:42:35  pwagner
! Passes MLSFiles to GetHDF5Attribute, LoadFromHDF5DS
!
! Revision 2.74  2005/06/22 17:25:49  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.73  2004/06/10 00:57:47  vsnyder
! Move FindFirst, FindNext from MLSCommon to MLSSets
!
! Revision 2.72  2004/02/11 01:08:44  livesey
! Removed print statement.
!
! Revision 2.71  2004/01/30 23:27:32  livesey
! Capitalized the names of the bins on output
!
! Revision 2.70  2004/01/24 01:44:36  livesey
! Removed print statement
!
! Revision 2.69  2004/01/24 01:01:22  livesey
! Improvements to the adoption stuff
!
! Revision 2.68  2004/01/23 19:07:49  livesey
! More on adoption / loading
!
! Revision 2.67  2004/01/23 05:37:34  livesey
! Added the adoption stuff
!
! Revision 2.66  2003/09/15 17:45:03  livesey
! Added target declaration for fussy intel compiler
!
! Revision 2.65  2003/08/14 20:24:23  livesey
! Added the exact bin selector stuff
!
! Revision 2.64  2003/08/13 00:47:34  livesey
! Added the default bin selectors stuff
!
! Revision 2.63  2003/08/08 23:04:45  livesey
! Added the dontPack stuff for output.
!
! Revision 2.62  2003/06/20 19:31:39  pwagner
! Changes to allow direct writing of products
!
! Revision 2.61  2003/05/13 04:46:24  livesey
! Renamed a routine to avoid conflict with VectorHDF5
!
! Revision 2.60  2003/02/06 01:37:12  livesey
! Set solarZenith to zero if reading from hdf5 and not present
!
! Revision 2.59  2003/02/06 01:34:21  livesey
! Added writing of solarZenith, and reading if it is present.
!
! Revision 2.58  2003/02/06 00:45:37  livesey
! Added name fragment to binSelector_t
!
! Revision 2.57  2003/02/05 21:55:45  livesey
! Bin selectors don't contain signal information anymore.
!
! Revision 2.56  2003/01/28 02:41:10  livesey
! Bug fix in writing matrix, now gets right element from database in all
! cases.
!
! Revision 2.55  2003/01/13 19:28:20  pwagner
! Moved several uses to speed Lahey buggs compiler
!
! Revision 2.54  2002/11/25 11:47:34  mjf
! Put SaveAsHDF5DS back to writing r4 reals.
!
! Revision 2.53  2002/11/22 12:47:47  mjf
! Added nullify routine(s) to get round Sun's WS6 compiler not
! initialising derived type function results.
!
! Revision 2.52  2002/11/20 21:06:25  livesey
! Various bug fixes to make the packed storage work.
!
! Revision 2.51  2002/10/08 00:09:10  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.50  2002/10/05 00:42:31  livesey
! Modified to use new pack/unpack literals/strings
!
! Revision 2.49  2002/10/02 23:20:19  livesey
! Added radiometer and saving as single precision
!
! Revision 2.48  2002/09/11 17:43:38  pwagner
! Began changes needed to conform with matrix%values type move to rm from r8
!
! Revision 2.47  2002/08/28 21:43:38  livesey
! Cosmetic changes
!
! Revision 2.46  2002/08/28 20:42:39  livesey
! Now uses inflateQuantityTemplateDatabase, much much faster on reading
! l2pcs.
!
! Revision 2.45  2002/08/28 01:12:20  livesey
! Bug fix in frequency coordinate.
!
! Revision 2.44  2002/08/23 01:23:44  livesey
! Now optionally reads frequency coordinate and frequencies
!
! Revision 2.43  2002/08/21 23:10:11  livesey
! No longer scans each bin for blocks
!
! Revision 2.42  2002/08/20 21:01:40  livesey
! Bug fix in DestroyL2PCInfoDatabase
!
! Revision 2.41  2002/08/07 00:05:14  livesey
! Moved H5Open/Close_F into tree walker.  Made some error messages more
! informative.
!
! Revision 2.40  2002/07/25 15:53:40  mjf
! Initialised some elements of qt in ReadOneVectorFromASCII (already
! done in ReadOneVectorFromHDF5)
!
! Revision 2.39  2002/07/22 03:29:23  livesey
! Bug fix
!
! Revision 2.38  2002/07/22 03:25:23  livesey
! Bug fix and rework of loading vectors.
!
! Revision 2.37  2002/07/17 06:00:21  livesey
! Got hdf5 l2pc reading stuff working
!
! Revision 2.36  2002/07/11 22:21:00  pwagner
! These hdf5-savvy versions transferred from he5lib
!
! Revision 1.3  2002/07/09 17:38:17  livesey
! Added logBasis attribute
!
! Revision 1.2  2002/06/22 23:12:13  livesey
! Added sparsity dump
!
! Revision 1.1  2002/06/12 17:59:53  livesey
! Work on HDF5 stuff
!
! Revision 2.32  2002/05/21 01:13:24  livesey
! New file format includes name for reading into IDL
!
! Revision 2.31  2002/03/15 21:22:42  livesey
! Slight modification to binSelectors stuff
!
! Revision 2.30  2002/03/13 22:00:53  livesey
! Added output of vertical coordinate for xStar/yStar.
! Note reading routine deduces value I think.
!
! Revision 2.29  2002/03/13 01:28:48  livesey
! Commented out an error message I think I don't need
!
! Revision 2.28  2002/02/08 22:51:28  livesey
! Minor changes and tidy up
!
! Revision 2.27  2002/02/05 04:14:26  livesey
! Bug fix, still problems with packed write
!
! Revision 2.26  2002/01/23 00:50:50  pwagner
! Initialize binselectors to null
!
! Revision 2.25  2002/01/22 18:14:29  livesey
! Fixed typo.
!
! Revision 2.24  2002/01/21 23:11:49  livesey
! Completed BinSelectors support
!
! Revision 2.23  2002/01/21 21:13:42  livesey
! Added BinSelector definitions and support
!
! Revision 2.22  2002/01/18 00:34:23  livesey
! Added packed option to writeonel2pc, with supporting code.
!
! Revision 2.21  2001/06/21 22:45:43  livesey
! Found another one of those `clobbering the latest database item' gotchas in the code.
!
! Revision 2.20  2001/06/06 17:27:29  pwagner
! DEBUG parameter; checks on some more reads
!
! Revision 2.19  2001/05/23 22:00:59  livesey
! Changed counter start to avoid conflict with L2Parallel
!
! Revision 2.18  2001/05/02 20:24:03  livesey
! Removed some unused variables.
!
! Revision 2.17  2001/05/02 20:23:10  vsnyder
! Specify a name for a created vector
!
! Revision 2.16  2001/04/28 01:38:36  livesey
! Now maintains proper IDs for its own quantity and vector templates
!
! Revision 2.15  2001/04/28 01:26:50  livesey
! This one seems to work!
!
! Revision 2.14  2001/04/27 17:36:33  livesey
! An interim version
!
! Revision 2.13  2001/04/27 07:24:50  livesey
! Interim version.  Still loosing parts of kStar on add to database.
! Might this be a compiler bug?
!
! Revision 2.12  2001/04/27 07:05:28  livesey
! Not sure what happened, needed to restore the use statements.
!
! Revision 2.11  2001/04/26 23:55:17  livesey
! Interim version
!
! Revision 2.10  2001/04/26 22:32:04  vsnyder
! Alphabetize USEs and declarations
!
! Revision 2.9  2001/04/26 22:12:21  livesey
! Fixed, gets l_zeta, l_none
!
! Revision 2.8  2001/04/26 22:08:39  livesey
! Add check on vertical coordinates
!
! Revision 2.7  2001/04/26 20:02:26  livesey
! Made l2pc database a saved array in L2PC_m
!
! Revision 2.6  2001/04/26 19:33:03  livesey
! Working version, reads and writes, (but no arithmetic :-) )
!
! Revision 2.5  2001/04/26 02:48:08  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.4  2001/04/26 00:06:33  livesey
! Many changes. Working towards a working read routine
!
! Revision 2.3  2001/04/25 20:32:42  livesey
! Interim version, tidied up write
!
! Revision 2.2  2001/04/24 20:20:48  livesey
! Word bin dropped from various places e.g. type
!
! Revision 2.1  2001/04/24 20:07:44  livesey
! Moved in from l2
!

