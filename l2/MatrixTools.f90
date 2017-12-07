! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module MatrixTools                      ! Various tools for matrices

  ! This module provides some tools for dealing matrices not already present in
  ! MatrixModule_0 and MatrixModule_1.  In particular the DumpBlocks subroutine.

  use MatrixModule_0, only: &
    & M_ABSENT, M_BANDED, M_COLUMN_SPARSE, M_FULL, MATRIXELEMENT_T
  use MatrixModule_1, only: FINDBLOCK, MATRIX_T, RC_INFO
  use MLSKinds, only: R8, RM
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, PVMERRORMESSAGE
  use PVM, only: PVMDATADEFAULT, PVMFINITSEND, PVMFSEND
  use PVMIDL, only: PVMIDLPACK

  implicit none
  private

  public :: DumpBlocks, CombineChannelsInMatrix, PVMSendMatrix

  ! Local paramters
  integer, parameter :: MTXMSGTAG = 202

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! =====  Public procedures  ===================================

  ! --------------------------------------------------  DumpBlocks  ----
  subroutine DumpBlocks ( key, matrices, hessians )
    ! This routine can be called whenever a DumpBlocks command is issued in the
    ! l2cf.  It can be used to dump requested blocks from the l2cf.  It dumps
    ! the blocks specified by the Cartesian product of the rowQuantity and
    ! colQuantity fields.  If the rowChannels, colChannels, rowSurfaces,
    ! colSurfaces, rowInstances or colInstances fields are specified, they
    ! are used for all blocks.  If this is not what is desired, use a separate
    ! DumpBlocks command for each block.  If the noAbsent field is set, it does
    ! not dump absent blocks.

    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use Dump_0, only: Dump
    use Expr_m, only: Expr
    use HessianModule_0, only: H_Absent, HessianElement_T
    use HessianModule_1, only: Dump, Hessian_T
    use Init_Tables_Module, only: F_AllHessians, F_AllMatrices, &
      & F_ColChannels, F_ColInstances, F_ColQuantity, F_ColSurfaces, &
      & F_Details, F_Diagonal, F_Hessian, F_Matrix, F_NoAbsent, &
      & F_RowChannels, F_RowInstances, F_RowQuantity, F_RowSurfaces, &
      & F_Structure
    use Lexer_Core, only: Print_Source
    use MatrixModule_1, only: Dump, Dump_Struct, GetFromMatrixDatabase, &
      & Matrix_Database_t
    use MLSStringLists, only: SwitchDetail
    use MoreTree, only: Get_Boolean, Get_Field_Id
    use Output_M, only: NewLine, Output
    use String_Table, only: Display_String
    use Toggles, only: Switches
    use Tree, only: Decoration, NSons, Subtree, Where
    use VectorsModule, only: GetVectorQtyByTemplateIndex, &
      & VectorValue_T

    ! Dummy arguments
    integer, intent(in) :: Key          ! L2CF node
    type (Matrix_Database_T), dimension(:), pointer :: Matrices ! Matrix database
    type (Hessian_T), dimension(:), pointer :: Hessians

    ! Local variables
    integer :: Col                      ! Matrix or Hessian block column
    integer :: Col2                     ! Hessian block second column
    integer :: ColChannelsNode          ! Tree node
    integer :: ColInstance              ! Loop counter
    integer :: ColInstance2             ! Loop counter
    integer :: ColInstancesNode         ! Tree node
    integer :: ColQI, Colqi2            ! Index of column quantity within vector
    integer :: ColQuantityIx            ! Index in ColQIs array
    integer :: ColQuantityIx2           ! Index in ColQIs array for dumping Hessians
    integer :: ColQuantityNode          ! Tree node
    integer :: ColSurfacesNode          ! Tree node
    integer :: DetailReduction
    integer :: Details                  ! 0 => just shapes, >0 => values, default 1
    logical :: Diagonal                 ! Dump only the diagonal
    logical :: DoAny                    ! Any non-absent blocks?
    integer :: FieldIndex               ! Type for tree node
    integer :: HessianIndex             ! Hessian database index
    integer :: Hessian1, HessianEnd     ! Range for HessianIndex
    integer :: MatrixIndex              ! Matrix database index
    integer :: Matrix1, MatrixEnd       ! Range for MatrixIndex
    integer :: NColQ                    ! How many column quantities?
    logical :: NoAbsent                 ! Don't dump absent blocks
    integer :: NODE                     ! Loop counter
    integer :: NRowQ                    ! How many row quantities?
    integer :: ROW                      ! Matrix block row
    integer :: RowChannelsNode          ! Tree node
    integer :: RowInstance              ! Loop counter
    integer :: RowInstancesNode         ! Tree node
    integer :: RowQI                    ! Index of row quantity within vector
    integer :: RowQuantityIx            ! Index in RowQIs array
    integer :: RowQuantityNode          ! Tree node
    integer :: RowSurfacesNode          ! Tree node
    integer :: Son                      ! Tree node
    integer :: Units(2) ! of the Details expr -- known to be phyq_dimensionless
    double precision :: Values(2) ! of the Details expr

    integer, dimension(:), pointer :: ColInds  ! Which column instances?
    integer, dimension(:), pointer :: ColInds2 ! Which column instances for Hessians?
    integer, dimension(:), pointer :: ColQIs   ! Which column quantities?
    integer, dimension(:), pointer :: RowInds  ! Which row instances?
    integer, dimension(:), pointer :: RowQIs   ! Which row quantities?

    logical :: AllHessians   ! Dump all hessians
    logical :: AllMatrices   ! Dump all matrices
    logical :: Fail          ! No matrix to get from database
    logical :: Stru          ! Dump sparsity structure instead of values

    type (VectorValue_T), pointer :: ColQ   ! Col quantity
    type (VectorValue_T), pointer :: RowQ   ! Row quantity
    type (Matrix_T), pointer :: Matrix      ! The matrix to dump
    type (MatrixElement_T), pointer :: MB   ! A block from the matrix
    type (HessianElement_T), pointer :: HB  ! A block from the Hessian

    ! Error codes for Announce_Error
    integer, parameter :: Duplicate = 1     ! Duplicate quantity name specified
    integer, parameter :: NeedHessianDatabase = duplicate + 1 ! Need some matrix!
    integer, parameter :: NeedMatrixDatabase = NeedHessianDatabase + 1 ! Need some matrix!
    integer, parameter :: NeedSomething = NeedMatrixDatabase + 1  ! Need /all or matrix
    integer, parameter :: NoSuchQuantity = NeedSomething + 1  ! for row or column selection
    integer, parameter :: OutOfRange = NoSuchQuantity + 1     ! Index out of range
    integer, parameter :: QuantAndMol = outOfRange + 1 ! Both quantity and molecule
    integer, parameter :: Redundant = quantAndMol + 1  ! Both /all and matrix

    ! Executable code

    ! Don't do it if the "nodb" switch is set.
    if ( switchDetail(switches, 'nodb') > -1 ) return

    ! Nor if we reduce the details level sufficiently
    DetailReduction = switchDetail(switches, 'red')
    if ( DetailReduction < 0 ) then ! The 'red' switch is absent
      DetailReduction = 0
    else if ( DetailReduction == 0 ) then ! By default, reduce details level by 2
      DetailReduction = 2
    end if


    ! Set defaults

    allHessians = .false.
    allMatrices = .false.
    colChannelsNode = 0
    colInstancesNode = 0
    colQuantityNode = 0
    colSurfacesNode = 0
    details = 1 - detailReduction
    diagonal = .false.
    hessianIndex = -1
    hessian1 = 1
    hessianEnd = -1 ! In case no Hessians are to be dumped
    matrixIndex = -1
    matrix1 = 1
    matrixEnd = -1  ! In case no Matrices are to be dumped
    noAbsent = .false.
    rowChannelsNode = 0
    rowInstancesNode = 0
    rowQuantityNode = 0
    rowSurfacesNode = 0
    stru = .false.

    nullify ( colInds, colInds2, colQIs, rowInds, rowQIs )

    ! First go through the parsed information.
    do node = 2, nsons(key)                ! Skip the DumpBlocks son
      son = subtree(node,key)              ! This argument
      fieldIndex = get_field_id(son)       ! ID for this field
      select case ( fieldIndex )
      case ( f_allHessians )
        if ( .not. associated(hessians) ) then
          call announce_error ( needHessianDatabase, son )
          return
        end if
        allHessians = get_Boolean ( son )
      case ( f_allMatrices )
        if ( .not. associated(matrices) ) then
          call announce_error ( needMatrixDatabase, son )
          return
        end if
        allMatrices = get_Boolean ( son )
      case ( f_details )
        call expr ( subtree(2,son), units, values )
        details = nint(values(1)) - detailReduction
      case ( f_diagonal )
        diagonal = get_Boolean ( son )
      case ( f_hessian )
        ! The decoration is the negative of the index; see Fill, where
        ! the Hessian spec is processed.
        hessianIndex = -decoration(decoration(subtree(2,son)))
      case ( f_matrix )
        ! And yet the corresponding index for the matrix is not the negative
        ! Whatever justification was invoked for the Hessian's sign change 
        ! was cheerfully ignored in the matrix's case.
        ! Long term--let's forget the sign changes.
        matrixIndex = decoration(decoration(subtree(2,son)))
      case ( f_rowQuantity )
        rowQuantityNode = son
      case ( f_colQuantity )
        colQuantityNode = son
      case ( f_rowChannels )
        rowChannelsNode = son
      case ( f_colChannels )
        colChannelsNode = son
      case ( f_rowSurfaces )
        rowSurfacesNode = son
      case ( f_colSurfaces )
        colSurfacesNode = son
      case ( f_rowInstances )
        rowInstancesNode = son
      case ( f_colInstances )
        colInstancesNode = son
      case ( f_noAbsent )
        noAbsent = get_Boolean ( son )
      case ( f_structure )
        stru = get_Boolean ( son )
      case default ! shouldn't get here if the type checker worked
      end select
    end do

    ! Was a matrix or hessian specified?
    if ( .not. allMatrices .and. matrixIndex < 0 .and. &
      &  .not. allHessians .and. hessianIndex < 0 ) then
      call announce_error ( needSomething, key )
      return
    end if
    
    if ( details < -1 ) return ! Don't do it if details too small

    if ( allMatrices ) then
      if ( matrixIndex > 0 ) call announce_error ( redundant, key )
      if ( details < 0 ) then
        call dump ( matrices, details=details )
        return
      end if
      matrix1 = 1
      matrixEnd = size(matrices)
    else if ( matrixIndex > 0 ) then
      matrix1 = matrixIndex
      matrixEnd = matrixIndex
    end if

    if ( allHessians ) then
      if ( hessianIndex > 0 ) call announce_error ( redundant, key )
      if ( details < 0 ) then
        call dump ( hessians, details=details )
        return
      end if
      hessian1 = 1
      matrixEnd = size(hessians)
    else if ( hessianIndex > 0 ) then
      hessian1 = hessianIndex
      hessianEnd = hessianIndex
    end if

    do matrixIndex = matrix1, matrixEnd
      ! Identify the matrix
      call GetFromMatrixDatabase ( matrices(matrixIndex), matrix, fail )
      if ( fail ) then
        call output ( matrixIndex, &
          & before='No matrix to get from database at index ', advance='yes' )
        cycle
      end if
      call output ( 'Dump of ' )
      if ( diagonal ) call output ( 'diagonal of ' )
      call display_string ( matrix%name )
      call print_source ( where(key), before=' at ' )
      if ( matrix%row%vec%name /= 0 ) &
        & call display_string ( matrix%row%vec%name, before=', row vector ' )
      if ( matrix%col%vec%name /= 0 ) &
        & call display_string ( matrix%col%vec%name, before=', column vector ' )
      if ( noAbsent ) call output ( ', /noAbsent' )
      call newLine
      if ( stru ) then
        call dump_struct ( matrix )
      else

        ! Get the row and column quantities
        call allocate_test ( colQIs, matrix%col%nb, 'colQIs', moduleName )
        call allocate_test ( rowQIs, matrix%row%nb, 'rowQIs', moduleName )

        call getQuantities ( matrix%col, colQIs, nColQ, colQuantityNode, 'column' )
        call getQuantities ( matrix%row, rowQIs, nRowQ, rowQuantityNode, 'row' )

        ! Dump the specified blocks
        do rowQuantityIx = 1, nRowQ
          rowQI = rowQIs(rowQuantityIx)
          rowQ => matrix%row%vec%quantities(rowQI)
          ! Fill some flags arrays
          call FillIndicesArray ( rowInstancesNode, rowQ%template%noInstances, &
            & rowInds )
          do colQuantityIx = 1, nColQ
            colQI = colQIs(colQuantityIx)
            colQ => matrix%col%vec%quantities(colQI)

            ! Fill some flags arrays
            call FillIndicesArray ( colInstancesNode, colQ%template%noInstances, &
              & colInds )

            doAny = .not. noAbsent
            if ( noAbsent ) then
      o:      do colInstance = 1, size(colInds)
                do rowInstance = 1, size(rowInds)
                  row = FindBlock ( matrix%row, rowQI, rowInds(rowInstance) )
                  col = FindBlock ( matrix%col, colQI, colInds(colInstance) )
                  if ( .not. diagonal .or. row == col ) then
                    mb => matrix%block ( row, col )
                    doAny = mb%kind /= m_absent
                    if ( doAny ) exit o
                  end if
                end do
              end do o
            end if

            if ( doAny ) call DumpOneMatrixBlock
            call deallocate_test ( colInds, 'colInds', ModuleName )

          end do ! colQuantityIx

          call deallocate_test ( rowInds, 'rowInds', ModuleName )
        end do ! rowQuantityIx
      end if

    end do ! matrixIndex

    do hessianIndex = hessian1, hessianEnd
      call output ( 'Dump of ' )
      call display_string ( hessians(hessianIndex)%name )
      call print_source ( where(key), before=' at ' )
      if ( hessians(hessianIndex)%row%vec%name /= 0 ) &
        & call display_string ( hessians(hessianIndex)%row%vec%name, before=', row vector ' )
      if ( hessians(hessianIndex)%col%vec%name /= 0 ) &
        & call display_string ( hessians(hessianIndex)%col%vec%name, before=', column vectors ' )
      if ( noAbsent ) call output ( ', /noAbsent' )
      call newLine

      ! Get the row and column quantities
      call allocate_test ( colQIs, hessians(hessianIndex)%col%nb, 'colQIs', moduleName )
      call allocate_test ( rowQIs, hessians(hessianIndex)%row%nb, 'rowQIs', moduleName )

      call getQuantities ( hessians(hessianIndex)%col, colQIs, nColQ, colQuantityNode, 'column' )
      call getQuantities ( hessians(hessianIndex)%row, rowQIs, nRowQ, rowQuantityNode, 'row' )

      ! Dump the specified blocks
      do rowQuantityIx = 1, nRowQ
        rowQI = rowQIs(rowQuantityIx)
        rowQ => hessians(hessianIndex)%row%vec%quantities(rowQI)
        ! Fill some flags arrays
        call FillIndicesArray ( rowInstancesNode, rowQ%template%noInstances, &
          & rowInds )
        do colQuantityIx = 1, nColQ
          colQI = colQIs(colQuantityIx)
          colQ => hessians(hessianIndex)%col%vec%quantities(colQI)

          ! Fill some flags arrays
          call FillIndicesArray ( colInstancesNode, colQ%template%noInstances, &
            & colInds )
          do colQuantityIx2 = 1, nColQ
            colQI2 = colQIs(colQuantityIx2)

            ! Fill some flags arrays
            call FillIndicesArray ( colInstancesNode, colQ%template%noInstances, &
              & colInds2 )

            ! Find the block to dump.
            nullify ( hb )
      o2:   do rowInstance = 1, size(rowInds)
              do colInstance = 1, size(colInds)
                do colInstance2 = 1, size(colInds2)
                  row = FindBlock ( hessians(hessianIndex)%row, rowQI, rowInds(rowInstance) )
                  col = FindBlock ( hessians(hessianIndex)%col, colQI, colInds(colInstance) )
                  col2 = FindBlock ( hessians(hessianIndex)%col, colQI2, colInds2(colInstance2) )
                  hb => hessians(hessianIndex)%block ( row, col, col2 )
                  if ( associated(hb) ) then
                    if ( hb%kind /= h_absent .or. .not. noAbsent ) then
                      if ( details >= 0 ) then
                        call display_string ( &
                          & hessians(hessianIndex)%row%vec%quantities(row)%template%name, before='[' )
                        call display_string ( &
                          & hessians(hessianIndex)%col%vec%quantities(col)%template%name, before=',' )
                        call display_string ( &
                          & hessians(hessianIndex)%col%vec%quantities(col2)%template%name, before=',' )
                        call output ( '] ' )
                      end if
                      call DumpOneHessianBlock ( hb, (/row,col,col2/) )
                      exit O2
                    end if
                  end if
                end do
              end do
            end do o2

          end do ! colQuantityIx2
        end do ! colQuantityIx
      end do ! rowQuantityIx

      call deallocate_test ( rowInds, 'rowInds', ModuleName )
      call deallocate_test ( colInds, 'colInds', ModuleName )
      call deallocate_test ( colInds2, 'colInds2', ModuleName )

    end do ! hessianIndex

    call deallocate_test ( colQIs, 'colQIs', moduleName )
    call deallocate_test ( rowQIs, 'rowQIs', moduleName )

  contains
    ! ...........................................  Announce_Error  .....
    subroutine Announce_Error ( What, Where, Number, Text )
      use MoreTree, only: StartErrorMessage
      use TREE, only: SUB_ROSA
      integer, intent(in) :: What             ! Error code 
      integer, intent(in) :: Where            ! Tree node  
      integer, intent(in), optional :: Number ! to stick into message
      character(len=*), intent(in), optional :: Text ! to stick into message

      call startErrorMessage ( where )
      select case ( what )
      case ( duplicate )
        call display_string ( sub_rosa(where), before=': Duplicate quantity ' )
        call output ( ' not used.', advance='yes' )
      case ( needHessianDatabase )
        call output ( 'There is no Hessian database', advance='yes' )
      case ( needMatrixDatabase )
        call output ( 'There is no matrix database', advance='yes' )
      case ( needSomething )
        call output ( &
          & 'Need at least one of /allMatrices, matrix, /allHessians, hessian', &
          & advance='yes' )
      case ( noSuchQuantity )
        call display_string ( sub_rosa(where) )
        call output ( ' is not a '//text//' quantity.', advance='yes' )
      case ( outOfRange )
        call output ( number, before=': Index ', after=' is out of range.', &
          & advance='yes' )
      case ( quantAndMol )
        call output ( 'Both quantities and molecules specified; molecules used.', &
          & advance='yes' )
      case ( redundant )
        call output ( ': Both /allMatrices and a Matrix specified, ' // &
          & 'or /allMatrices and a Hessian specified; /all used.', &
          & advance='yes' )
      end select
    end subroutine Announce_Error

    ! .......................................  DumpOneMatrixBlock  .....
    subroutine DumpOneMatrixBlock

      use MatrixModule_0, only: Densify, GetDiagonal

      ! Local variables
      integer :: CC                       ! Loop counter
      integer :: CS                       ! Loop counter
      integer :: NoColChannels            ! Number selected
      integer :: NoColSurfaces            ! Number selected
      integer :: NoRowChannels            ! Number selected
      integer :: NoRowSurfaces            ! Number selected
      integer :: RC                       ! Loop counter
      integer :: RS                       ! Loop counter

      integer, dimension(:), pointer :: RowChanInds ! Indices
      integer, dimension(:), pointer :: ColChanInds ! Indices
      integer, dimension(:), pointer :: RowSurfInds ! Indices
      integer, dimension(:), pointer :: ColSurfInds ! Indices

      integer :: ColInstance, RowInstance         ! Do index variables
      real(rm), dimension(:,:), pointer :: Val    ! The values from the block
      real(rm), dimension(:), pointer :: Val_1D   ! The diagonal
      real(r8), dimension(:,:), pointer :: ToDump ! The 2D matrix to dump

      nullify ( rowChanInds, colChanInds )
      nullify ( rowSurfInds, colSurfInds )
      nullify ( toDump, val_1d )

      ! Set up the index arrays
      call getSurfOrChanInds ( rowQ%template%noChans, rowChannelsNode, &
        & rowChanInds, 'Row Channels: ' )
      call getSurfOrChanInds ( rowQ%template%noSurfs, rowSurfacesNode, &
        & rowSurfInds, 'Row Surfaces: ' )
      call getSurfOrChanInds ( colQ%template%noChans, colChannelsNode, &
        & colChanInds, 'Column Channels: ' )
      call getSurfOrChanInds ( colQ%template%noSurfs, colSurfacesNode, &
        & colSurfInds, 'Column Surfaces: ' )

      noRowChannels = size(rowChanInds)
      noRowSurfaces = size(rowSurfInds)
      noColChannels = size(colChanInds)
      noColSurfaces = size(colSurfInds)

      if ( .not. diagonal ) &
        & call allocate_test ( toDump, &
          & noRowChannels*noRowSurfaces, &
          & noColChannels*noColSurfaces, &
          & 'toDump', ModuleName )

      ! Loop over the row and column instances
      do colInstance = 1, size(colInds)
        do rowInstance = 1, size(rowInds)

          row = FindBlock ( matrix%row, rowQI, rowInds(rowInstance) )
          col = FindBlock ( matrix%col, colQI, colInds(colInstance) )

          mb => matrix%block ( row, col )
          if ( noAbsent .and. mb%kind == m_absent ) cycle

          ! Dump a header
          call display_string ( &
            & matrix%row%vec%quantities(matrix%row%quant(row))%template%name )
          call output ( matrix%row%inst(row), before=':', after=', ' )
          call display_string ( &
            & matrix%col%vec%quantities(matrix%col%quant(col))%template%name )
          call output ( matrix%col%inst(col), before=':', after=' (' )
          call output ( row, after=',' )
          call output ( col, after=') is ' )
          nullify ( val )
          select case ( mb%kind )
          case ( m_absent )
            call output ( 'absent' )
          case ( m_column_sparse, m_banded )
            if ( mb%kind == m_column_sparse ) then
              call output ( 'column sparse' )
            else
              call output ( 'banded' )
            end if
            if ( diagonal ) then
              call allocate_test ( val_1d, min(mb%nRows, mb%nCols), 'val_1d', &
                & moduleName )
              call getDiagonal ( mb, val_1d )
            else
              call allocate_test ( val, mb%nRows, mb%nCols, 'val', &
                & ModuleName )
              call densify ( val , mb )
            end if
          case ( m_full )
            call output ( 'full' )
            if ( diagonal ) then
              call allocate_test ( val_1d, min(mb%nRows, mb%nCols), 'val_1d', &
                & moduleName )
              call getDiagonal ( mb, val_1d )
            else
              val => mb%values
            end if
          case default
          end select
          if ( .not. diagonal ) then
            call output ( noRowChannels*noRowSurfaces, before=' ' )
            call output ( noColChannels*noColSurfaces, before='x' )
          end if

          if ( mb%kind /= m_absent .and. details > 0 ) then
            if ( diagonal ) then
              call dump ( val_1d, name=', number dumped:', options=what_options(clean=.true.) )
              call deallocate_test ( val_1d, 'val_1d', moduleName )
            else
              do cs = 1, noColSurfaces
                do cc = 1, noColChannels
                  do rs = 1, noRowSurfaces
                    do rc = 1, noRowChannels
                      todump ( rc + (rs-1)*noRowChannels, &
                        &      cc + (cs-1)*noColChannels ) = &
                        & val ( rowChanInds(rc) + &
                        &      (rowSurfInds(rs)-1)*rowQ%template%noChans, &
                        &       colChanInds(cc) + &
                        &      (colSurfInds(cs)-1)*colQ%template%noChans )
                    end do
                  end do
                end do
              end do
              if ( mb%kind /= m_full ) &
                & call deallocate_test ( val, 'val', ModuleName )
              call dump ( toDump, name=', number dumped:', options=what_options(clean=.true.) )
            end if
          else
            call newLine
          end if

        end do
      end do

      call deallocate_test ( toDump,       'toDump',       ModuleName )
      call deallocate_test ( rowChanInds,  'rowChanInds',  ModuleName )
      call deallocate_test ( colChanInds,  'colChanInds',  ModuleName )
      call deallocate_test ( rowSurfInds,  'rowSurfInds',  ModuleName )
      call deallocate_test ( colSurfInds,  'colSurfInds',  ModuleName )

    end subroutine DumpOneMatrixBlock

    ! ......................................  DumpOneHessianBlock  .....
    subroutine DumpOneHessianBlock ( HB, Indices )
      use HessianModule_0, only: Dump, HessianElement_T
      type(HessianElement_T), intent(in) :: HB
      integer, intent(in) :: Indices(:)
      ! For now ignore the row and column instances
      call dump ( hb, details=details, indices=indices )
    end subroutine DumpOneHessianBlock

    ! .........................................  FillIndicesArray  .....
    subroutine FillIndicesArray ( Node, Num, Inds )
    ! Fill Inds with sons 2..n of Node, or 1..Num if Node == 0

      use Declaration_Table, only: NUM_VALUE, RANGE
      use Expr_M, only: EXPR

      integer, intent(in) :: Node         ! Tree node
      integer, intent(in) :: Num          ! Maximum size of Inds
      integer, pointer :: Inds(:)         ! Sons of Node, or 1...n

      logical :: Error
      logical :: Flags(num)               ! Flags(son of Node) = .true.
      integer :: I
      integer :: SonIx                    ! Son Index for Node
      integer :: TYPE                     ! From expr
      integer, dimension(2) :: UNITS      ! Units from expr
      real(r8), dimension(2) :: VALUE     ! Value from expr

      nullify ( inds )

      if ( node /= 0 ) then
        error = .false.
        flags = .false.
        do sonIx = 2, nsons(node)
          call expr (subtree(sonIx,node), units, value, type)
          if ( nint(value(1)) < 1 .or. nint(value(1)) > num ) then
            error = .true.
            call announce_error ( outOfRange, subtree(sonIx,node), nint(value(1)) )
          end if
          select case (type)
          case (num_value)
            if ( .not. error ) flags(nint(value(1))) = .true.
          case (range)
            if ( nint(value(2)) < 1 .or. nint(value(2)) > num ) then
              error = .true.
              call announce_error ( outOfRange, subtree(sonIx,node), nint(value(2)) )
            end if
            if ( .not. error ) flags(nint(value(1)):nint(value(2))) = .true.
          case default
          end select
        end do
        if ( error ) &
          & call MLSMessage ( MLSMSG_Error, moduleName, 'Index out of range' )
        call allocate_test ( inds, count(flags), 'Inds', moduleName )
        inds = pack ( (/ (i, i=1, num) /), flags )
      else
        call allocate_test ( inds, num, 'Inds', moduleName )
        inds = (/ (i, i=1, num) /)
      end if

    end subroutine FillIndicesArray

    ! ............................................  GetQuantities  .....
    subroutine GetQuantities ( RC, QIs, NQIs, QInode, Text )
    ! Get quantities from the tree or the matrix.
      use MatrixModule_1, only: RC_Info
      use TREE, only: SUBTREE
      type(RC_Info), intent(in) :: RC  ! Row or Column info for Matrix
      integer, intent(out) :: QIs(:)   ! Quantity indices
      integer, intent(out) :: NQIs     ! Number of QIs actually used
      integer, intent(in) :: QInode    ! Tree node for quantities -- zero if not specified
      character(len=*), intent(in) :: Text  ! 'row' or 'column' for error message

      integer :: I                     ! Subscript, loop inductor
      integer :: QI                    ! A quantity index -- may go into QIs
      integer :: Son                   ! of Node
      type(vectorValue_t), pointer :: Q  ! A vector quantity

      nQIs = 0
      if ( qInode /= 0 ) then
        do i = 2, nsons(qInode)
          son = subtree(i,qInode)
          ! Identify the quantity.  The decoration is an index into the
          ! quantity templates database.  We need to get it as an index into
          ! the vector that describes the row or column.
          Q => GetVectorQtyByTemplateIndex ( rc%vec, &
            & decoration(decoration(son)), qi ) ! qi is output too
          if ( .not. associated (Q) ) then
            call announce_error ( noSuchQuantity, son, text=text )
          else if ( any(QIs(:nQIs) == qi ) ) then
            call announce_error ( duplicate, son )
          else
            nQIs = nQIs + 1
            QIs(nQIs) = qi
          end if
        end do
      else
        do i = 1, rc%nb
          qi = rc%quant(i)
          if ( all(QIs(:nQIs) /= qi ) ) then
            nQIs = nQIs + 1
            QIs(nQIs) = qi
          end if
        end do
      end if
    end subroutine GetQuantities

    ! ........................................  GetSurfOrChanInds  .....
    subroutine GetSurfOrChanInds ( Num, Node, Inds, Text )
      integer, intent(in) :: Num              ! vec%template%no[Chans,Surfs]     
      integer, intent(in) :: Node             ! Tree node                        
      integer, pointer :: Inds(:)             ! Indices                          
      character(len=*), intent(in) :: Text    ! For the dump                     

      call fillIndicesArray ( node, num, inds )

      call output ( text )
      if ( size(inds) == 0 ) then
        call output ( '(none)', advance='yes' )
      else if ( size(inds) == 1 ) then
        call output ( inds(1), advance='yes' )
      else if ( any(inds(2:)-inds(1:size(inds)-1) /= 1) ) then
        call newLine
        call dump ( inds )
      else
        call output ( inds(1), after=':' )
        call output ( inds(size(inds)), advance='yes' )
      end if

    end subroutine GetSurfOrChanInds

  end subroutine DumpBlocks

  ! ------------------------------------  CombineChannelsInMatrix  -----
  subroutine CombineChannelsInMatrix ( mOut, mIn )
    ! This subroutine might belong better in L2PC_m, as it's really only for
    ! working with that kind of matrix.  It takes a matrix with a
    ! single quantity for the row vector (i.e. radiance) and multiple
    ! quantiites for the column vector and remaps it to a 'downsampled'
    ! set of channels
    use ManipulateVectorQuantities, only: FILLWITHCOMBINEDCHANNELS
    use MatrixModule_1, only: CLEARMATRIX
    use MatrixModule_0, only: MULTIPLYMATRIX_XY
    use Output_m, only: OUTPUT
    use String_Table, only: DISPLAY_STRING
    type (Matrix_T), intent(inout) :: MOUT ! Result matrix
    type (Matrix_T), intent(in) :: MIN ! Source matrix
    
    ! Local variables
    type (MatrixElement_T) :: MAPPING   ! A mapping block
    character (len=80) :: MESSAGE       ! A possible error message
    integer :: ROW, COL                 ! Loop counters

    ! First some sanity checks
    ! Check that the column vectors are compatible
    if ( mOut%col%vec%template%name /= mIn%col%vec%template%name ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Matrices do not share a column vector in CombineChannels' )
    ! Check that the row vectors contain only a single quantity.
    if ( mOut%row%vec%template%noQuantities /= 1 .or. &
      &  mIn%row%vec%template%noQuantities /= 1 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Matrices must have single-quantity row vector in CombineChannels' )
    
    ! Now get the quantity filling command to do the first part of the work for us.
    ! This also checks the sanity of mapping mIn%row%vec to mOut%row%vec
    call FillWithCombinedChannels ( &
      & mOut%row%vec%quantities(1), mIn%row%vec%quantities(1), &
      & message, mapping )
    if ( message /= '' ) call MLSMessage ( MLSMSG_Error, ModuleName, message )

    ! OK, that's filled the vector now do the matrix
    call ClearMatrix ( mOut )

    do row = 1, mOut%row%nb
      do col = 1, mOut%col%nb
        ! Put out a useful progress message
        call output ( 'Block [ ' )
        call output ( row )
        call output ( ', ' )
        call output ( col )
        call output ( ' ] -- [ ' )
        call display_string ( mIn%row%vec%quantities(mIn%row%quant(row))%template%name, strip=.true. )
        call output ( '(' )
        call output ( mIn%row%inst(row) )
        call output ( '), ' )
        call display_string ( mIn%col%vec%quantities(mIn%col%quant(col))%template%name, strip=.true. )
        call output ( '(' )
        call output ( mIn%col%inst(col) )
        call output ( ') ]: ' )

        if ( mIn%block(row,col)%kind /= m_absent ) then
          ! This does all the creation etc.
          call output ( 'present', advance='yes' )
          call MultiplyMatrix_XY ( mapping, mIn%block(row,col), mOut%block(row,col) )
        else
          call output ( 'absent', advance='yes' )
        end if
      end do
    end do

  end subroutine CombineChannelsInMatrix

  ! -----------------------------------------------  PVMSendBlock  -----
  subroutine PVMPackBlock ( BLOCK )
    ! Dummy arguments
    type (MatrixElement_T), intent(in) :: BLOCK ! The block of the matrix

    ! Local variables
    integer :: INFO                     ! Flag

    ! Executable code
    call PVMIDLPack ( (/ block%kind, block%nRows, block%nCols /), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing block info" )

    if ( any (block%kind == (/M_Banded, M_Column_Sparse /) ) ) then
      call PVMIDLPack ( (/ size(block%values) /), info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing number of values" )
      call PVMIDLPack ( block%R1, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing R1" )
      call PVMIDLPack ( block%R2, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing R2" )
    end if

    if ( block%kind /= M_Absent ) then
      call PVMIDLPack ( real(block%values,kind=r8), info )
    end if

  end subroutine PVMPackBlock

  ! --------------------------------------------------  PVMSendRC  -----
  subroutine PVMPackRC ( RC )
    ! Dummy argument
    use QuantityPVM, only: PVMSENDQUANTITY
    type ( RC_Info ), intent(in) :: RC  ! the rcinfo to send

    ! Local variables
    integer :: INFO                     ! Flag
    integer :: QTY                      ! Loop counter

    ! Executable code
    ! Pack the number of blocks
    call PVMIDLPack ( rc%nb, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Packing NB' )

    ! Pack instFirst
    call PVMIDLPack(rc%instFirst, info)
    if (info /= 0) call PVMErrorMessage (info, 'Packing instFirst')

    ! Pack the indices
    call PVMIDLPack ( rc%nelts, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Packing Nelts' )
    call PVMIDLPack ( rc%inst, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Packing Inst' )
    call PVMIDLPack ( rc%quant, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Packing Quant' )

    ! Pack the size of the vector
    call PVMIDLPack ( size (rc%vec%quantities), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Packing no quantities' )
    do qty = 1, size (rc%vec%quantities)
      call PVMSendQuantity ( rc%vec%quantities(qty), 0, &
      & justPack=.true., noValues=.true., noMask=.true. )
    end do

  end subroutine PVMPackRC

  ! ----------------------------------------------  PVMSendMatrix  -----
  subroutine PVMSendMatrix ( MATRIX, TID, JUSTPACK )
    ! Dummy arguments
    type (Matrix_T), intent(in) :: MATRIX ! The matrix to send
    integer, intent(in) :: TID          ! The task to send it to
    logical, intent(in), optional :: JUSTPACK ! 

    ! Local variables
    logical :: MYJUSTPACK               ! Copy of justPack
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! Flag from PVM
    integer :: I,J                      ! Loop counters

    ! Executable code
    myJustPack = .false.
    if (present(justPack)) myJustPack = justPack

    if (.not. myJustPack) then
        call PVMFInitSend ( PvmDataDefault, bufferID )
    end if
    call PVMPackRC ( matrix%col )
    call PVMPackRC ( matrix%row )

    do j = 1, size(matrix%block,2)
      do i = 1, size(matrix%block,1)
        call PVMPackBlock ( matrix%block(i,j) )
      end do
    end do

    if (.not. myJustPack) then
      call PVMFSend ( tid, MtxMsgTag, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "sending vector values" )
    end if

  end subroutine PVMSendMatrix

  function What_Options( Clean, Transpose, Trim ) result( Options )
    use Dump_Options, only: Dopt_Clean, Dopt_Transpose, Dopt_Trim
    use MLSStrings, only: trim_safe
    logical, optional, intent(in) :: Clean
    logical, optional, intent(in) :: Transpose
    logical, optional, intent(in) :: Trim
    character(len=8) :: options
    options = ' '
    if ( present(clean) ) then
      if ( clean ) options = trim_safe(options) // Dopt_Clean
    end if
    if ( present(transpose) ) then
      if ( transpose ) options = trim_safe(options) // Dopt_Transpose
    end if
    if ( present(trim) ) then
      if ( trim ) options = trim_safe(options) // Dopt_Trim
    end if
  end function What_Options

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MatrixTools

! $Log$
! Revision 1.48  2017/07/27 01:40:43  vsnyder
! Print the correct sizes of rows and columns
!
! Revision 1.47  2016/07/28 03:27:34  vsnyder
! Cannonball polishing
!
! Revision 1.46  2016/07/28 00:40:34  vsnyder
! Remove unreferenced USE
!
! Revision 1.45  2015/03/28 02:49:25  vsnyder
! Added stuff to trace allocate/deallocate addresses
!
! Revision 1.44  2014/09/05 01:13:10  vsnyder
! More complete and accurate allocate/deallocate size tracking.  Convert
! some local pointer temps to allocatable.  Get kinds from MLSKinds instead
! of from MLSCommon.
!
! Revision 1.43  2014/09/05 00:49:07  vsnyder
! EmpiricalGeometry.f90 -- Wrong comment
!
! Revision 1.42  2014/08/06 23:33:45  vsnyder
! Remove USE for Num_Value, which is not referenced
!
! Revision 1.41  2014/02/28 01:08:20  vsnyder
! Remove unused names
!
! Revision 1.40  2014/02/28 00:21:12  vsnyder
! Move type and units checking to type checker
!
! Revision 1.39  2013/09/24 23:47:22  vsnyder
! Use Where instead of Source_Ref for messages
!
! Revision 1.38  2013/06/12 02:38:02  vsnyder
! Cruft removal
!
! Revision 1.37  2012/01/27 02:57:54  vsnyder
! Remove extraneous ' at ' in DumpBlocks
!
! Revision 1.36  2011/10/04 20:25:33  honghanh
! Fixed bug in PVM_Pack_RC
!
! Revision 1.35  2011/09/01 20:37:08  honghanh
! Fix the bug in PVMSendMatrix for sending NB
!
! Revision 1.34  2011/08/20 00:49:37  vsnyder
! Remove unused use names and variable declarations
!
! Revision 1.33  2010/08/06 23:02:36  pwagner
! Moved to using only switchdetail; negative index deplored
!
! Revision 1.32  2010/03/24 20:51:43  vsnyder
! Add code to dump Hessians.  Spiff up some error messages.
!
! Revision 1.31  2010/02/10 20:00:25  vsnyder
! More output from the dump
!
! Revision 1.30  2009/12/17 23:56:41  vsnyder
! Correct a formatting blunder
!
! Revision 1.29  2009/12/17 23:44:43  vsnyder
! Print block coordinates in dump
!
! Revision 1.28  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 1.27  2009/06/16 17:40:19  pwagner
! Changed api for dump, diff routines; now rely on options for most optional behavior
!
! Revision 1.26  2007/10/05 23:41:06  vsnyder
! Don't reference element 1 of a zero-size array
!
! Revision 1.25  2007/10/02 22:41:27  vsnyder
! Don't crash if a matrix can't be dumped
!
! Revision 1.24  2006/09/21 18:48:07  pwagner
! Reduce level of dumps in SIDS version
!
! Revision 1.23  2006/09/20 00:43:21  vsnyder
! Cannonball polishing
!
! Revision 1.22  2006/09/19 20:33:13  vsnyder
! Add /diagonal field
!
! Revision 1.21  2006/07/27 03:53:28  vsnyder
! Handle details field correctly
!
! Revision 1.20  2006/07/19 22:27:24  vsnyder
! Add /allMatrices, details= and /structure fields
!
! Revision 1.19  2005/06/22 18:57:02  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 1.18  2005/03/15 23:50:16  pwagner
! PVMERRORMESSAGE now part of MLSMessageModule
!
! Revision 1.17  2004/10/29 20:54:11  vsnyder
! Remove USE for DUMP -- not referenced, some cosmetics
!
! Revision 1.16  2004/09/28 00:40:06  livesey
! More informative log in CombineChannelsInMatrix
!
! Revision 1.15  2004/09/25 00:16:53  livesey
! Added CombineChannelsInMatrix
!
! Revision 1.14  2004/01/21 22:00:46  vsnyder
! Remove unused variable declarations
!
! Revision 1.13  2003/10/10 23:28:33  vsnyder
! Substantial reorganization
!
! Revision 1.12  2003/10/09 21:03:16  vsnyder
! Don't do anything if the 'nodb' switch is set
!
! Revision 1.11  2003/10/07 01:17:36  vsnyder
! DumpBlocks now dumps the blocks specified by the Cartesian product of the
! rowQuantity and colQuantity fields.  If the rowChannels, colChannels,
! rowSurfaces, colSurfaces, rowInstances or colInstances fields are specified,
! they are used for all blocks.  If this is not what is desired, use a separate
! DumpBlocks command for each block.  If the noAbsent field is set, it does
! not dump absent blocks.
!
! Revision 1.10  2003/08/15 20:28:05  vsnyder
! Put a new line after 'with Y rows and X columns' for absent blocks
!
! Revision 1.9  2003/05/21 19:19:06  vsnyder
! Plug some memory leaks.  Use "name" argument of "dump" routines.  Dump
! numbers of rows and columns.
!
! Revision 1.8  2002/10/08 17:36:21  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 1.7  2002/09/13 18:10:10  pwagner
! May change matrix precision rm from r8
!
! Revision 1.6  2001/09/13 00:54:06  livesey
! Fixed a bug with dump blocks dumping values when kind==m_absent
!
! Revision 1.5  2001/07/17 17:32:15  livesey
! Added PVM pack stuff
!
! Revision 1.4  2001/05/08 21:33:23  livesey
! Added the CVS log stuff, whoops!
!
