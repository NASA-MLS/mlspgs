! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MatrixTools                      ! Various tools for matrices

  ! This module provides some tools for dealing matrices not already present in
  ! MatrixModule_0 and MatrixModule_1.  In particular the DumpBlock subroutine.

  use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
  use DUMP_0, only: DUMP
  use MatrixModule_1, only: MATRIX_T, MATRIX_DATABASE_T, &
    & FINDBLOCK, GETFROMMATRIXDATABASE
  use MatrixModule_0, only: MATRIXELEMENT_T, DENSIFY, &
    & M_ABSENT, M_BANDED, M_COLUMN_SPARSE, M_FULL
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use Init_Tables_Module, only: F_COLCHANNELS, F_COLQUANTITY, F_COLSURFACES, &
    & F_MATRIX, F_ROWCHANNELS, F_ROWQUANTITY, F_ROWSURFACES, F_ROWINSTANCES, &
    & F_COLINSTANCES
  use Tree, only: NSONS, SUBTREE, DECORATION
  use MoreTree, only: GET_FIELD_ID
  use VectorsModule, only: GETVECTORQTYBYTEMPLATEINDEX, VECTORVALUE_T
  use Declaration_Table, only: NUM_VALUE, RANGE
  use Output_M, only: OUTPUT
  use String_Table, only: DISPLAY_STRING
  use Expr_M, only: EXPR

  implicit none
  private

  public :: DumpBlock

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  !---------------------------------------------------------------------------

contains ! ================ Public procedures ================================

  subroutine DumpBlock ( key, matrices )
    ! This routine can be called whenever a DumpBlock command is issued in the
    ! l2cf.  It can be used to dump a requested block from the l2cf

    ! Dummy arguments
    integer, intent(in) :: KEY          ! L2CF node
    type (Matrix_Database_T), dimension(:), pointer :: MATRICES ! Matrix database

    ! Local variables
    integer :: CC                       ! Loop counter
    integer :: COL                      ! Matrix block column
    integer :: COLCHANNELSNODE          ! Tree node
    integer :: COLINSTANCE              ! Loop counter
    integer :: COLINSTANCESNODE         ! Tree node
    integer :: COLQI                    ! Index of column quantity within vector
    integer :: COLQUANTITYINDEX         ! Index for column quantity
    integer :: COLSURFACESNODE          ! Tree node
    integer :: CS                       ! Loop counter
    integer :: FIELDINDEX               ! Type for tree node
    integer :: GSON                     ! Tree node
    integer :: I                        ! Loop counter
    integer :: MATRIXINDEX              ! Matrix database index
    integer :: MAXIND                   ! Maximum dimension
    integer :: NOCOLCHANNELS            ! Number selected
    integer :: NOCOLSURFACES            ! Number selected
    integer :: NODE                     ! Loop counter
    integer :: NOROWCHANNELS            ! Number selected
    integer :: NOROWSURFACES            ! Number selected
    integer :: RC                       ! Loop counter
    integer :: ROW                      ! Matrix block row
    integer :: ROWCHANNELSNODE          ! Tree node
    integer :: ROWINSTANCE              ! Loop counter
    integer :: ROWINSTANCESNODE         ! Tree node
    integer :: ROWQI                    ! Index of row quantity within vector
    integer :: ROWQUANTITYINDEX         ! Index for row quantity
    integer :: ROWSURFACESNODE          ! Tree node
    integer :: RS                       ! Loop counter
    integer :: SON                      ! Tree node
    integer :: TYPE                     ! From expr

    integer, dimension(2) :: UNITS      ! Units from expr
    integer, dimension(:), pointer :: INDGEN ! 1,2,3,4...
    integer, dimension(:), pointer :: ROWCHANINDS ! Indices
    integer, dimension(:), pointer :: COLCHANINDS ! Indices
    integer, dimension(:), pointer :: ROWSURFINDS ! Indices
    integer, dimension(:), pointer :: COLSURFINDS ! Indices

    real(r8), dimension(2) :: VALUE     ! Value from expr

    logical, dimension(:), pointer :: ROWCHANNELS ! Do we want this channel?
    logical, dimension(:), pointer :: COLCHANNELS ! Do we want this channel?
    logical, dimension(:), pointer :: ROWSURFACES ! Do we want this surface?
    logical, dimension(:), pointer :: COLSURFACES ! Do we want this surface?
    logical, dimension(:), pointer :: ROWINSTANCES ! Do we want this surface?
    logical, dimension(:), pointer :: COLINSTANCES ! Do we want this surface?

    real(r8), dimension(:,:), pointer :: VAL ! The values from the block
    real(r8), dimension(:,:), pointer :: TODUMP ! The 2D matrix to dump

    type (Matrix_T), pointer :: MATRIX  ! The matrix to dump
    type (VectorValue_T), pointer :: ROWQ ! Row quantity
    type (VectorValue_T), pointer :: COLQ ! Row quantity
    type (MatrixElement_T), pointer :: MB ! A block from the matrix

    ! Executable code

    nullify ( rowChannels, colChannels )
    nullify ( rowSurfaces, colSurfaces )
    nullify ( rowInstances, colInstances )
    nullify ( rowChanInds, colChanInds )
    nullify ( rowSurfInds, colSurfInds )
    nullify ( indgen )
    nullify ( toDump )

    ! Set defaults
    rowChannelsNode = 0
    colChannelsNode = 0
    rowSurfacesNode = 0
    colSurfacesNode = 0
    rowInstancesNode = 0
    colInstancesNode = 0

    ! First go through the parsed information.
    do node = 2, nsons(key)                ! Skip the DumpBlock son
      son = subtree(node,key)              ! This argument
      fieldIndex = get_field_id(son)    ! ID for this field
      if (nsons(son) > 1 ) gson = subtree(2,son)
      select case ( fieldIndex )
      case ( f_matrix )
        matrixIndex = decoration(decoration(gson))
      case ( f_rowQuantity )
        rowQuantityIndex = decoration(decoration(gson))
      case ( f_colQuantity )
        colQuantityIndex = decoration(decoration(gson))
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
      case default
      end select
    end do

    ! Identify the matrix
    call GetFromMatrixDatabase ( matrices(matrixIndex), matrix )

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

    ! Now setup flag arrays
    call allocate_test (rowChannels, rowQ%template%noChans, &
      &'rowChannels', ModuleName )
    call allocate_test (colChannels, colQ%template%noChans, &
      & 'colChannels', ModuleName )
    call allocate_test (rowSurfaces, rowQ%template%noSurfs, &
      & 'rowSurfaces', ModuleName )
    call allocate_test (colSurfaces, colQ%template%noSurfs, &
      & 'colSurfaces', ModuleName )
    call allocate_test (rowInstances, rowQ%template%noInstances, &
      & 'rowInstances', ModuleName )
    call allocate_test (colInstances, colQ%template%noInstances, &
      & 'colInstances', ModuleName )

    ! Now fill these arrays
    if ( rowChannelsNode /= 0 ) then
      rowChannels = .false.
      do node = 2, nsons(rowChannelsNode)
        call expr (subtree(node,rowChannelsNode), units, value, type)
        select case (type)
        case (num_value)
          rowChannels(int(value(1))) = .true.
        case (range)
          rowChannels(int(value(1)):int(value(2))) = .true.
        case default
        end select
      end do
    else
      rowChannels = .true.
    end if

    if ( colChannelsNode /= 0 ) then
      colChannels = .false.
      do node = 2, nsons(colChannelsNode)
        call expr (subtree(node,colChannelsNode), units, value, type)
        select case (type)
        case (num_value)
          colChannels(int(value(1))) = .true.
        case (range)
          colChannels(int(value(1)):int(value(2))) = .true.
        case default
        end select
      end do
    else
      colChannels = .true.
    end if

    if ( rowSurfacesNode /= 0 ) then
      rowSurfaces = .false.
      do node = 2, nsons(rowSurfacesNode)
        call expr (subtree(node,rowSurfacesNode), units, value, type)
        select case (type)
        case (num_value)
          rowSurfaces(int(value(1))) = .true.
        case (range)
          rowSurfaces(int(value(1)):int(value(2))) = .true.
        case default
        end select
      end do
    else
      rowSurfaces = .true.
    end if

    if ( colSurfacesNode /= 0 ) then
      colSurfaces = .false.
      do node = 2, nsons(colSurfacesNode)
        call expr (subtree(node,colSurfacesNode), units, value, type)
        select case (type)
        case (num_value)
          colSurfaces(int(value(1))) = .true.
        case (range)
          colSurfaces(int(value(1)):int(value(2))) = .true.
        case default
        end select
      end do
    else
      colSurfaces = .true.
    end if

    if ( rowInstancesNode /= 0 ) then
      rowInstances = .false.
      do node = 2, nsons(rowInstancesNode)
        call expr (subtree(node,rowInstancesNode), units, value, type)
        select case (type)
        case (num_value)
          rowInstances(int(value(1))) = .true.
        case (range)
          rowInstances(int(value(1)):int(value(2))) = .true.
        case default
        end select
      end do
    else
      rowInstances = .true.
    end if

    if ( colInstancesNode /= 0 ) then
      colInstances = .false.
      do node = 2, nsons(colInstancesNode)
        call expr (subtree(node,colInstancesNode), units, value, type)
        select case (type)
        case (num_value)
          colInstances(int(value(1))) = .true.
        case (range)
          colInstances(int(value(1)):int(value(2))) = .true.
        case default
        end select
      end do
    else
      colInstances = .true.
    end if

    ! Now set indices arrays
    noRowChannels = count ( rowChannels )
    noColChannels = count ( colChannels )
    noRowSurfaces = count ( rowSurfaces )
    noColSurfaces = count ( colSurfaces )

    maxInd = max ( rowQ%template%noChans, colQ%template%noChans,&
      & rowQ%template%noSurfs, colQ%template%noSurfs )
    call allocate_test ( indgen, maxInd, 'indgen', ModuleName )
    do i = 1, maxInd
      indgen(i) = i
    end do

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

    call allocate_test ( toDump, &
      & noRowChannels*noRowSurfaces, &
      & noColChannels*noColSurfaces, &
      & 'toDump', ModuleName )

    ! Now loop over the row and column instances
    call output ( 'Dump of ' )
    call display_string ( matrix%name, advance='yes' )
    call output ( 'Row channels:', advance='yes' )
    call dump ( rowChanInds )
    call output ( 'Row surfaces:', advance='yes' )
    call dump ( rowSurfInds )
    call output ( 'Column channels:', advance='yes' )
    call dump ( colChanInds )
    call output ( 'Column surfaces:', advance='yes' )
    call dump ( colSurfInds )
    do colInstance = 1, colQ%template%noInstances
      if ( colInstances(colInstance) ) then
        do rowInstance = 1, rowQ%template%noInstances
          if ( rowInstances(rowInstance) ) then
            ! Dump a header
            call output ( 'Block for ' )
            call display_string ( rowQ%template%name )
            call output ( ' instance ' )            
            call output ( rowInstance )
            call output ( ', ' )
            call display_string ( colQ%template%name )
            call output ( ' instance ' )            
            call output ( colInstance )

            row = FindBlock ( matrix%row, rowQI, rowInstance )
            col = FindBlock ( matrix%col, colQI, colInstance )

            mb => matrix%block ( row, col )
            nullify ( val )
            select case ( mb%kind )
            case ( m_absent )
              call output ( ' is absent.', advance='yes' )
            case ( m_column_sparse, m_banded )
              if ( mb%kind == m_column_sparse ) then
                call output ( ' is column sparse.', advance='yes' )
              else
                call output ( ' is banded.', advance='yes' )
              end if
              call allocate_test ( val, mb%nRows, mb%nCols, &
                & 'val', ModuleName )
              call densify ( val , mb )
            case ( m_full )
              call output ( ' is full.', advance='yes' )
              val => mb%values
            case default
            end select

            if ( associated(val) ) then
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
            endif

            call output ( 'Number of elements dumped:\ ' )
            call output ( size(toDump), advance='yes' )
            call dump ( toDump, clean=.true. )

          end if
        end do
      end if
    end do

    call deallocate_test ( toDump, 'toDump', ModuleName )
    call deallocate_test ( indgen, 'indgen', ModuleName )
    
  end subroutine DumpBlock

end module MatrixTools

! $Log$
