! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MatrixTools                      ! Various tools for matrices

  ! This module provides some tools for dealing matrices not already present in
  ! MatrixModule_0 and MatrixModule_1.  In particular the DumpBlocks subroutine.

  use PVM, only: PVMDATADEFAULT, PVMFINITSEND, PVMFSEND, PVMERRORMESSAGE
  use PVMIDL, only: PVMIDLPACK
  use MatrixModule_1, only: MATRIX_T, MATRIX_DATABASE_T, &
    & FINDBLOCK, GETFROMMATRIXDATABASE, RC_INFO
  use MatrixModule_0, only: MATRIXELEMENT_T, DENSIFY, &
    & M_ABSENT, M_BANDED, M_COLUMN_SPARSE, M_FULL
  use MLSCommon, only: R8, RM
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error

  implicit none
  private

  public :: DumpBlocks, CombineChannelsInMatrix

  ! Local paramters
  integer, parameter :: MTXMSGTAG = 202

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), save :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains ! =====  Public procedures  ===================================

  ! --------------------------------------------------  DumpBlocks  ----
  subroutine DumpBlocks ( key, matrices )
    ! This routine can be called whenever a DumpBlocks command is issued in the
    ! l2cf.  It can be used to dump requested blocks from the l2cf.  It dumps
    ! the blocks specified by the Cartesian product of the rowQuantity and
    ! colQuantity fields.  If the rowChannels, colChannels, rowSurfaces,
    ! colSurfaces, rowInstances or colInstances fields are specified, they
    ! are used for all blocks.  If this is not what is desired, use a separate
    ! DumpBlocks command for each block.  If the noAbsent field is set, it does
    ! not dump absent blocks.

    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use DUMP_0, only: DUMP
    use Init_Tables_Module, only: F_COLCHANNELS, F_COLINSTANCES, F_COLQUANTITY, &
      & F_COLSURFACES, F_MATRIX, F_NOABSENT, &
      & F_ROWCHANNELS, F_ROWINSTANCES, F_ROWQUANTITY, F_ROWSURFACES
    use MoreTree, only: GET_BOOLEAN, GET_FIELD_ID
    use Output_M, only: NewLine, OUTPUT
    use String_Table, only: DISPLAY_STRING
    use Toggles, only: Switches
    use Tree, only: NSONS, SUBTREE, DECORATION
    use VectorsModule, only: GETVECTORQTYBYTEMPLATEINDEX, VECTORVALUE_T

    ! Dummy arguments
    integer, intent(in) :: KEY          ! L2CF node
    type (Matrix_Database_T), dimension(:), pointer :: MATRICES ! Matrix database

    ! Local variables
    integer :: COL                      ! Matrix block column
    integer :: COLCHANNELSNODE          ! Tree node
    integer :: COLINSTANCE              ! Loop counter
    integer :: COLINSTANCESNODE         ! Tree node
    integer :: COLQI                    ! Index of column quantity within vector
    integer :: COLQuantityIx            ! Index in ColQIs array
    integer :: COLQuantityNode          ! Tree node
    integer :: COLSURFACESNODE          ! Tree node
    logical :: DoAny                    ! Any non-absent blocks?
    integer :: FIELDINDEX               ! Type for tree node
    integer :: MATRIXINDEX              ! Matrix database index
    integer :: NColQ                    ! How many column quantities?
    logical :: NoAbsent                 ! Don't dump absent blocks
    integer :: NODE                     ! Loop counter
    integer :: NRowQ                    ! How many row quantities?
    integer :: ROW                      ! Matrix block row
    integer :: ROWCHANNELSNODE          ! Tree node
    integer :: ROWINSTANCE              ! Loop counter
    integer :: ROWINSTANCESNODE         ! Tree node
    integer :: ROWQI                    ! Index of row quantity within vector
    integer :: ROWQuantityIx            ! Index in RowQIs array
    integer :: ROWQuantityNode          ! Tree node
    integer :: ROWSURFACESNODE          ! Tree node
    integer :: SON                      ! Tree node

    integer, dimension(:), pointer :: ColInds ! Which column instances?
    integer, dimension(:), pointer :: ColQIs  ! Which column quantities?
    integer, dimension(:), pointer :: RowInds ! Which row instances?
    integer, dimension(:), pointer :: RowQIs  ! Which row quantities?

    type (VectorValue_T), pointer :: COLQ ! Row quantity
    type (VectorValue_T), pointer :: ROWQ ! Row quantity
    type (Matrix_T), pointer :: MATRIX  ! The matrix to dump
    type (MatrixElement_T), pointer :: MB ! A block from the matrix

    ! Error codes for Announce_Error
    integer, parameter :: Duplicate = 1   ! Duplicate quantity name specified
    integer, parameter :: OutOfRange = duplicate + 1 ! Index out of range

    ! Executable code

    ! Don't do it if the "nodb" switch is set.
    if ( index(switches, 'nodb') /= 0 ) return

    ! Set defaults
    rowChannelsNode = 0
    colChannelsNode = 0
    rowQuantityNode = 0
    colQuantityNode = 0
    rowSurfacesNode = 0
    colSurfacesNode = 0
    rowInstancesNode = 0
    colInstancesNode = 0
    noAbsent = .false.

    nullify ( colInds, colQIs, rowInds, rowQIs )

    ! First go through the parsed information.
    do node = 2, nsons(key)                ! Skip the DumpBlocks son
      son = subtree(node,key)              ! This argument
      fieldIndex = get_field_id(son)       ! ID for this field
      select case ( fieldIndex )
      case ( f_matrix )
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
      case default ! shouldn't get here if the type checker worked
      end select
    end do

    ! Identify the matrix
    call GetFromMatrixDatabase ( matrices(matrixIndex), matrix )
    call output ( 'Dump of ' )
    call display_string ( matrix%name, advance='yes' )

    ! Get the row and column quantities
    call allocate_test ( colQIs, matrix%col%nb, 'colQIs', moduleName )
    call allocate_test ( rowQIs, matrix%row%nb, 'rowQIs', moduleName )

    call getQuantities ( matrix%col, colQIs, nColQ, colQuantityNode, 'column' )
    call getQuantities ( matrix%row, rowQIs, nRowQ, rowQuantityNode, 'row' )

    ! Dump the specified blocks
    do rowQuantityIx = 1, nRowQ
      rowQI = rowQIs(rowQuantityIx)
      rowQ => matrix%row%vec%quantities(rowQI)
      do colQuantityIx = 1, nColQ
        colQI = colQIs(colQuantityIx)
        colQ => matrix%col%vec%quantities(colQI)

        ! Fill some flags arrays
        call FillIndicesArray ( rowInstancesNode, rowQ%template%noInstances, &
          & rowInds )
        call FillIndicesArray ( colInstancesNode, colQ%template%noInstances, &
          & colInds )

        doAny = .not. noAbsent
        if ( noAbsent ) then
  o:      do colInstance = 1, size(colInds)
            do rowInstance = 1, size(rowInds)
              row = FindBlock ( matrix%row, rowQI, rowInds(rowInstance) )
              col = FindBlock ( matrix%col, colQI, colInds(colInstance) )
              mb => matrix%block ( row, col )
              doAny = mb%kind /= m_absent
              if ( doAny ) exit o
            end do
          end do o
        end if

        if ( doAny ) call DumpOneBlock

        call deallocate_test ( rowInds, 'rowInds', ModuleName )
        call deallocate_test ( colInds, 'colInds', ModuleName )

      end do
    end do

    call deallocate_test ( colQIs, 'colQIs', moduleName )
    call deallocate_test ( rowQIs, 'rowQIs', moduleName )

  contains
    ! ...........................................  Announce_Error  .....
    subroutine Announce_Error ( What, Where, Number )
      use LEXER_CORE, only: PRINT_SOURCE
      use TREE, only: SOURCE_REF, SUB_ROSA
      integer, intent(in) :: What      ! Error code
      integer, intent(in) :: Where     ! Tree node
      integer, intent(in) :: Number    ! Stuff to stick into message

      call output ( '***** At ' )
      call print_source ( source_ref(where) )
      select case ( what )
      case ( duplicate )
        call output ( ': Duplicate quantity ' )
        call display_string ( sub_rosa(where) )
        call output ( ' not used.', advance='yes' )
      case ( outOfRange )
        call output ( number, before=': Index ', after=' is out of range.', &
          & advance='yes' )
      end select
    end subroutine Announce_Error

    ! .............................................  DumpOneBlock  .....
    subroutine DumpOneBlock

      ! Local variables
      integer :: CC                       ! Loop counter
      integer :: CS                       ! Loop counter
      integer :: NOCOLCHANNELS            ! Number selected
      integer :: NOCOLSURFACES            ! Number selected
      integer :: NOROWCHANNELS            ! Number selected
      integer :: NOROWSURFACES            ! Number selected
      integer :: RC                       ! Loop counter
      integer :: RS                       ! Loop counter

      integer, dimension(:), pointer :: ROWCHANINDS ! Indices
      integer, dimension(:), pointer :: COLCHANINDS ! Indices
      integer, dimension(:), pointer :: ROWSURFINDS ! Indices
      integer, dimension(:), pointer :: COLSURFINDS ! Indices

      real(rm), dimension(:,:), pointer :: VAL    ! The values from the block
      real(r8), dimension(:,:), pointer :: TODUMP ! The 2D matrix to dump

      nullify ( rowChanInds, colChanInds )
      nullify ( rowSurfInds, colSurfInds )
      nullify ( toDump )

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

      call allocate_test ( toDump, &
        & noRowChannels*noRowSurfaces, &
        & noColChannels*noColSurfaces, &
        & 'toDump', ModuleName )

      call newLine
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
          call output ( matrix%col%inst(col), before=':' )
          nullify ( val )
          select case ( mb%kind )
          case ( m_absent )
            call output ( ' is absent ' )
          case ( m_column_sparse, m_banded )
            if ( mb%kind == m_column_sparse ) then
              call output ( ' is column sparse ' )
            else
              call output ( ' is banded ' )
            end if
            call allocate_test ( val, mb%nRows, mb%nCols, &
              & 'val', ModuleName )
            call densify ( val , mb )
          case ( m_full )
            call output ( ' is full ' )
            val => mb%values
          case default
          end select
          call output ( mb%nRows )
          call output ( mb%nCols, before='x', after='.  ' )

          if ( mb%kind /= m_absent ) then
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
            call dump ( toDump, name='Number dumped:', clean=.true. )
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

    end subroutine DumpOneBlock

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
      call allocate_test ( inds, num, 'Inds', moduleName )

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
      integer, intent(in) :: QInode    ! Tree node -- zero if not specified
      character(len=*), intent(in) :: Text  ! 'row' or 'column' for error message

      integer :: I                     ! Subscript, loop inductor
      integer :: QI                    ! A quantity index -- may go into QIs
      type(vectorValue_t), pointer :: Q  ! A vector quantity

      nQIs = 0
      if ( QInode /= 0 ) then
        do i = 2, nsons(QInode)
          ! Identify the row quantity.  The decoration is an index into the
          ! quantity templates database.  We need to get it as an index into
          ! the vector that describes the row.
          Q => GetVectorQtyByTemplateIndex ( rc%vec, &
            & decoration(decoration(subtree(i,QInode))), qi ) ! qi is output too
          if ( .not. associated (Q) ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
            & text // ' quantity was not found.' )
          if ( any(QIs(:nQIs) == qi ) ) then
            call announce_error ( duplicate, subtree(i,QInode), subtree(i,QInode) )
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
      if ( size(inds) == 1 ) then
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
    real(r8), dimension(:,:), pointer :: values
    integer :: n_rows, n_cols

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
!      call PVMIDLPack ( block%values, info )
        n_rows = size(block%values, 1)
        n_cols = size(block%values, 2)
        allocate(values(n_rows, n_cols))
        values = block%values
        call PVMIDLPack ( values, info )
        deallocate(values)
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing block values" )
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
    call PVMIDLPack ( (/ rc%nelts /), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'Packing NB' )

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

    if (.not. myJustPack) call PVMFInitSend ( PvmDataDefault, bufferID )
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

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MatrixTools

! $Log$
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
