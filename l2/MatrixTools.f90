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

  public :: DumpBlocks

  ! Local paramters
  integer, parameter :: MTXMSGTAG = 202

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
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

    ! Dummy arguments
    integer, intent(in) :: KEY          ! L2CF node
    type (Matrix_Database_T), dimension(:), pointer :: MATRICES ! Matrix database

    ! Local variables
    integer :: COLCHANNELSNODE          ! Tree node
    integer :: COLINSTANCESNODE         ! Tree node
    integer :: COLQUANTITYINDEX         ! Index for column quantity
    integer :: COLQuantityNode          ! Tree node
    integer :: COLQuantitySon           ! Index for son of COLQuantityNode
    integer :: COLSURFACESNODE          ! Tree node
    integer :: FIELDINDEX               ! Type for tree node
    integer :: MATRIXINDEX              ! Matrix database index
    logical :: NoAbsent                 ! Don't dump absent blocks
    integer :: NODE                     ! Loop counter
    integer :: ROWCHANNELSNODE          ! Tree node
    integer :: ROWINSTANCESNODE         ! Tree node
    integer :: ROWQUANTITYINDEX         ! Index for row quantity
    integer :: ROWQuantityNode          ! Tree node
    integer :: ROWQuantitySon           ! Index for son of ROWQuantityNode
    integer :: ROWSURFACESNODE          ! Tree node
    integer :: SON                      ! Tree node

    type (Matrix_T), pointer :: MATRIX  ! The matrix to dump

    ! Executable code

    ! Don't do it if the "nodb" switch is set.
    if ( index(switches, 'nodb') /= 0 ) return

    ! Set defaults
    rowChannelsNode = 0
    colChannelsNode = 0
    rowSurfacesNode = 0
    colSurfacesNode = 0
    rowInstancesNode = 0
    colInstancesNode = 0
    noAbsent = .false.

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
      case default
      end select
    end do

    ! Identify the matrix
    call GetFromMatrixDatabase ( matrices(matrixIndex), matrix )
    call output ( 'Dump of ' )
    call display_string ( matrix%name, advance='yes' )

    ! Dump the specified blocks
    do rowQuantitySon = 2, nsons(rowQuantityNode)
      rowQuantityIndex = decoration(decoration(subtree(rowQuantitySon,rowQuantityNode)))
      do colQuantitySon = 2, nsons(colQuantityNode)
        colQuantityIndex = decoration(decoration(subtree(colQuantitySon,colQuantityNode)))
        call DumpOneBlock
      end do
    end do

  contains
    ! .............................................  DumpOneBlock  .....
    subroutine DumpOneBlock

      use VectorsModule, only: GETVECTORQTYBYTEMPLATEINDEX, VECTORVALUE_T

      ! Local variables
      integer :: CC                       ! Loop counter
      integer :: COL                      ! Matrix block column
      integer :: COLINSTANCE              ! Loop counter
      integer :: COLQI                    ! Index of column quantity within vector
      integer :: CS                       ! Loop counter
      logical :: DoAny                    ! Any non-absent blocks?
      integer :: I                        ! Loop counter
      integer :: MAXIND                   ! Maximum dimension
      integer :: NOCOLCHANNELS            ! Number selected
      integer :: NOCOLSURFACES            ! Number selected
      integer :: NOROWCHANNELS            ! Number selected
      integer :: NOROWSURFACES            ! Number selected
      integer :: RC                       ! Loop counter
      integer :: ROW                      ! Matrix block row
      integer :: ROWINSTANCE              ! Loop counter
      integer :: ROWQI                    ! Index of row quantity within vector
      integer :: RS                       ! Loop counter

      integer, dimension(:), pointer :: INDGEN ! 1,2,3,4...
      integer, dimension(:), pointer :: ROWCHANINDS ! Indices
      integer, dimension(:), pointer :: COLCHANINDS ! Indices
      integer, dimension(:), pointer :: ROWSURFINDS ! Indices
      integer, dimension(:), pointer :: COLSURFINDS ! Indices

      logical, dimension(:), pointer :: ROWCHANNELS ! Do we want this channel?
      logical, dimension(:), pointer :: COLCHANNELS ! Do we want this channel?
      logical, dimension(:), pointer :: ROWSURFACES ! Do we want this surface?
      logical, dimension(:), pointer :: COLSURFACES ! Do we want this surface?
      logical, dimension(:), pointer :: ROWINSTANCES ! Do we want this surface?
      logical, dimension(:), pointer :: COLINSTANCES ! Do we want this surface?

      real(rm), dimension(:,:), pointer :: VAL ! The values from the block
      real(r8), dimension(:,:), pointer :: TODUMP ! The 2D matrix to dump

      type (VectorValue_T), pointer :: ROWQ ! Row quantity
      type (VectorValue_T), pointer :: COLQ ! Row quantity
      type (MatrixElement_T), pointer :: MB ! A block from the matrix

      nullify ( rowChannels, colChannels )
      nullify ( rowSurfaces, colSurfaces )
      nullify ( rowInstances, colInstances )
      nullify ( rowChanInds, colChanInds )
      nullify ( rowSurfInds, colSurfInds )
      nullify ( indgen )
      nullify ( toDump )

      ! Now identify the row and column quantity rowQuantityIndex and
      ! colQuantityIndex are indices into the quantity templates database.  Now,
      ! we need to get them as indices into the vectors describing rows and
      ! columns.
      rowQ => GetVectorQtyByTemplateIndex ( matrix%row%vec, &
        & rowQuantityIndex, rowQI )
      colQ => GetVectorQtyByTemplateIndex ( matrix%col%vec, &
        & colQuantityIndex, colQI )

      if ( .not. associated (rowQ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'row quantity was not found.' )
      if ( .not. associated (colQ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'col quantity was not found.' )

      ! Now setup flags arrays
      call allocate_test ( rowChannels, rowQ%template%noChans, &
        & 'rowChannels', ModuleName )
      call allocate_test ( colChannels, colQ%template%noChans, &
        & 'colChannels', ModuleName )
      call allocate_test ( rowSurfaces, rowQ%template%noSurfs, &
        & 'rowSurfaces', ModuleName )
      call allocate_test ( colSurfaces, colQ%template%noSurfs, &
        & 'colSurfaces', ModuleName )
      call allocate_test ( rowInstances, rowQ%template%noInstances, &
        & 'rowInstances', ModuleName )
      call allocate_test ( colInstances, colQ%template%noInstances, &
        & 'colInstances', ModuleName )

      ! Now fill flags arrays
      call FillFlagsArray ( rowChannelsNode, rowChannels, noRowChannels )
      call FillFlagsArray ( colChannelsNode, colChannels, noColChannels )
      call FillFlagsArray ( rowSurfacesNode, rowSurfaces, noRowSurfaces )
      call FillFlagsArray ( colSurfacesNode, colSurfaces, noColSurfaces )
      call FillFlagsArray ( rowInstancesNode, rowInstances )
      call FillFlagsArray ( colInstancesNode, colInstances )

      ! Now set indices arrays
      maxInd = max ( rowQ%template%noChans, colQ%template%noChans,&
        & rowQ%template%noSurfs, colQ%template%noSurfs )
      call allocate_test ( indgen, maxInd, 'indgen', ModuleName )
      indgen = (/ ( i, i = 1, maxInd ) /)

      call allocate_test ( rowChanInds, noRowChannels, &
        & 'rowChanInds', ModuleName )
      call allocate_test ( colChanInds, noColChannels, &
        & 'colChanInds', ModuleName )
      call allocate_test ( rowSurfInds, noRowSurfaces, &
        & 'rowSurfInds', ModuleName )
      call allocate_test ( colSurfInds, noColSurfaces, &
        & 'colSurfInds', ModuleName )

      rowChanInds = pack ( indgen(1:rowQ%template%noChans), rowChannels )
      colChanInds = pack ( indgen(1:colQ%template%noChans), colChannels )
      rowSurfInds = pack ( indgen(1:rowQ%template%noSurfs), rowSurfaces )
      colSurfInds = pack ( indgen(1:colQ%template%noSurfs), colSurfaces )

      doAny = .not. noAbsent
      if ( noAbsent ) then
o:      do colInstance = 1, colQ%template%noInstances
          if ( colInstances(colInstance) ) then
            do rowInstance = 1, rowQ%template%noInstances
              if ( rowInstances(rowInstance) ) then
                row = FindBlock ( matrix%row, rowQI, rowInstance )
                col = FindBlock ( matrix%col, colQI, colInstance )
                mb => matrix%block ( row, col )
                doAny = mb%kind /= m_absent
                if ( doAny ) exit o
              end if
            end do
          end if
        end do o
      end if

      if ( doAny ) then
        call NewLine
        call dumpIndex ( rowChannels, rowChanInds, 'Row channels: ' )
        call dumpIndex ( rowSurfaces, rowSurfInds, 'Row Surfaces: ' )
        call dumpIndex ( colChannels, colChanInds, 'Column channels: ' )
        call dumpIndex ( colSurfaces, colSurfInds, 'Column surfaces: ' )
      end if

      call allocate_test ( toDump, &
        & noRowChannels*noRowSurfaces, &
        & noColChannels*noColSurfaces, &
        & 'toDump', ModuleName )

      ! Now loop over the row and column instances
      do colInstance = 1, colQ%template%noInstances
        if ( colInstances(colInstance) ) then
          do rowInstance = 1, rowQ%template%noInstances
            if ( rowInstances(rowInstance) ) then

              row = FindBlock ( matrix%row, rowQI, rowInstance )
              col = FindBlock ( matrix%col, colQI, colInstance )

              mb => matrix%block ( row, col )
              if ( noAbsent .and. mb%kind == m_absent ) cycle

              ! Dump a header
              call output ( 'Block for ' )
              call display_string ( rowQ%template%name )
              call output ( rowInstance, before=' instance ', after=', ' )
              call display_string ( colQ%template%name )
              call output ( colInstance, before=' instance ' )
              nullify ( val )
              select case ( mb%kind )
              case ( m_absent )
                call output ( ' is absent,', advance='yes' )
              case ( m_column_sparse, m_banded )
                if ( mb%kind == m_column_sparse ) then
                  call output ( ' is column sparse,', advance='yes' )
                else
                  call output ( ' is banded,', advance='yes' )
                end if
                call allocate_test ( val, mb%nRows, mb%nCols, &
                  & 'val', ModuleName )
                call densify ( val , mb )
              case ( m_full )
                call output ( ' is full,', advance='yes' )
                val => mb%values
              case default
              end select
              call output ( mb%nRows, before='with ' )
              call output ( mb%nCols, before=' rows and ', after=' columns.  ' )

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
                call dump ( toDump, name='Number of elements dumped:', clean=.true. )
              else
                call NewLine
              end if

            end if
          end do
        end if
      end do

      call deallocate_test ( toDump,       'toDump',       ModuleName )
      call deallocate_test ( indgen,       'indgen',       ModuleName )
      call deallocate_test ( rowChanInds,  'rowChanInds',  ModuleName )
      call deallocate_test ( colChanInds,  'colChanInds',  ModuleName )
      call deallocate_test ( rowSurfInds,  'rowSurfInds',  ModuleName )
      call deallocate_test ( colSurfInds,  'colSurfInds',  ModuleName )
      call deallocate_test ( rowChannels,  'rowChannels',  ModuleName )
      call deallocate_test ( colChannels,  'colChannels',  ModuleName )
      call deallocate_test ( rowSurfaces,  'rowSurfaces',  ModuleName )
      call deallocate_test ( colSurfaces,  'colSurfaces',  ModuleName )
      call deallocate_test ( rowInstances, 'rowInstances', ModuleName )
      call deallocate_test ( colInstances, 'colInstances', ModuleName )

    end subroutine DumpOneBlock

    ! ................................................  DumpIndex  .....
    subroutine DumpIndex ( Flags, Indices, Name )
    ! Dump an index array as 1:n if every index is there, else as-is.
      logical, intent(in) :: Flags(:)
      integer, intent(in) :: Indices(:)
      character(len=*), intent(in) :: Name
      call output ( Name )
      if ( .not. all(flags) ) then
        if ( size(indices) /= 1 ) call newLine
        call dump ( indices )
      else
        if ( size(flags) /= 1 ) call output ( '1:' )
        call output ( size(flags), advance='yes' )
      end if
    end subroutine DumpIndex

    ! ...........................................  FillFlagsArray  .....
    subroutine FillFlagsArray ( Node, Array, NFlags )
    ! Set Array(I) = .true. for every element that is a son of Node.

      use Declaration_Table, only: NUM_VALUE, RANGE
      use Expr_M, only: EXPR

      integer, intent(in) :: Node              ! Tree node
      logical, intent(out) :: Array(:)         ! Array to fill
      integer, intent(out), optional :: NFlags ! Count(Array)

      integer :: SonIx                    ! Son Index for Node
      integer :: TYPE                     ! From expr
      integer, dimension(2) :: UNITS      ! Units from expr
      real(r8), dimension(2) :: VALUE     ! Value from expr

      if ( node /= 0 ) then
        array = .false.
        do sonIx = 2, nsons(node)
          call expr (subtree(sonIx,node), units, value, type)
          select case (type)
          case (num_value)
            array(nint(value(1))) = .true.
          case (range)
            array(nint(value(1)):nint(value(2))) = .true.
          case default
          end select
        end do
      else
        array = .true.
      end if
      if ( present(nFlags) ) nFlags = count(array)

    end subroutine FillFlagsArray

  end subroutine DumpBlocks

  ! ----------------------------------------- PVMSendBlock
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

  ! ----------------------------------------- PVMSendRC
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

  ! ----------------------------------------- PVMSendMatrix
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
