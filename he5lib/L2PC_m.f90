! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module L2PC_m
  !=============================================================================

  ! This module contains data types etc. for dealing with the new EMLS L2PC
  ! files.  The first version will deal with ascii files, but later versions
  ! will probably be HDF.

  use Allocate_Deallocate, only: Allocate_test, Deallocate_test
  use Declaration_Table, only: DECLS, ENUM_VALUE, GET_DECL, DUMP_DECL
  use Intrinsic, only: Lit_Indices, L_CHANNEL, L_GEODALTITUDE, L_ZETA, L_NONE, L_VMR, &
    & L_RADIANCE, L_PTAN
  use machine, only: io_error
  use MLSCommon, only: R8
  use VectorsModule, only: assignment(=), DESTROYVECTORINFO, &
    & VECTORTEMPLATE_T, VECTOR_T, VECTORVALUE_T, CREATEVECTOR, ADDVECTORTODATABASE,&
    & ADDVECTORTEMPLATETODATABASE, CONSTRUCTVECTORTEMPLATE
  use MatrixModule_1, only: CREATEBLOCK, CREATEEMPTYMATRIX, &
    & DESTROYMATRIX, MATRIX_T, DUMP, FINDBLOCK, MATRIX_DATABASE_T, GETFROMMATRIXDATABASE, &
    & DUMP_STRUCT
  use MatrixModule_0, only: M_ABSENT, M_BANDED, M_COLUMN_SPARSE, M_FULL, &
    & MATRIXELEMENT_T
  use MLSHDF5, only: MakeHDF5Attribute, GetHDF5Attribute, SaveAsHDF5DS, LoadFromHDF5DS
  use HDF5, only: H5FCREATE_F, H5FCLOSE_F, H5F_ACC_TRUNC_F, &
    & H5GCLOSE_F, H5GCREATE_F, H5GOPEN_F, H5GGET_OBJ_INFO_IDX_F
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR, &
    & MLSMSG_ALLOCATE, MLSMSG_DEALLOCATE
  use MLSSignals_m, only: GETSIGNALNAME
  use Output_m, only: output
  use Parse_Signal_m, only: Parse_Signal
  use QuantityTemplates, only: ADDQUANTITYTEMPLATETODATABASE, QUANTITYTEMPLATE_T
  use String_Table, only: GET_STRING
  use Symbol_Table, only: ENTER_TERMINAL, DUMP_SYMBOL_CLASS
  use Symbol_Types, only: T_IDENTIFIER
  use TOGGLES, only: TAB, TOGGLE, SWITCHES
  use Tree, only: DECORATION, NSONS, SUBTREE

  implicit NONE
  private
  
  public :: AddL2PCToDatabase, DestroyL2PC, DestroyL2PCDatabase, WriteOneL2PC
  public :: Open_l2pc_file, read_l2pc_file, close_l2pc_file, binSelector_T
  public :: BinSelectors, DestroyBinSelectorDatabase,  AddBinSelectorToDatabase
  public :: OutputHDF5L2PC

  ! This is the third attempt to do this.  An l2pc is simply a Matrix_T.
  ! As this contains pointers to vector_T's and so on, I maintain a private
  ! set of databases of these in this module.  We can't use the main databases,
  ! as these would get destroyed at the end of each chunk.

  ! The l2pc database and supporting databases
  type(QuantityTemplate_T), dimension(:), pointer, save :: L2PCQTS => NULL()
  type(VectorTemplate_T), dimension(:), pointer, save :: L2PCVTS => NULL()
  type(Vector_T), dimension(:), pointer, save :: L2PCVS => NULL()
  type(Matrix_T), dimension(:), pointer, public, save :: L2PCDatabase => NULL()

  integer :: counterStart
  parameter ( counterStart = huge (0) / 4 )

  ! This datatype describes a selection rule for l2pc bins.
  type BinSelector_T
    integer :: selectorType             ! What quantity type does this apply to
    integer :: molecule                 ! What molecule does it apply to
    integer, dimension(:), pointer :: signals => NULL() ! What signals does this apply to
    integer, dimension(:), pointer :: sidebands => NULL() ! What sidebands
    real(r8), dimension(2) :: heightRange ! The height range for this selector
    real(r8) :: cost                    ! The cost for that range
  end type BinSelector_T

  type(BinSelector_T), dimension(:), pointer, save :: BINSELECTORS => NULL()

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  ! Local saved variables - these keep track of the l2pc vector/quantity databases
  integer, save :: L2PCQTCOUNTER = CounterStart ! To place in qt%id
  integer, save :: L2PCVTCOUNTER = CounterStart ! To place in vt%id
    
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

  ! ------------------------------------  Add l2pc  to database ----
  integer function AddL2PCToDatabase ( Database, Item )
    
    ! This function simply adds an l2pc  to a database of said l2pc s.
    
    type(Matrix_T), dimension(:), pointer :: Database
    type(Matrix_T) :: Item
    
    type(Matrix_T), dimension(:), pointer :: TempDatabase

    include "addItemToDatabase.f9h"

    AddL2PCToDatabase = newSize
  end function AddL2PCToDatabase

  ! -----------------------------------  Close_L2PC_File  -----
  subroutine Close_L2PC_File ( Lun )
    integer, intent(in) :: lun
    close ( lun )
  end subroutine Close_L2PC_File

  ! -------------------------------------- DestroyBinSelectorDatabase
  subroutine DestroyBinSelectorDatabase 
    ! Local variables
    integer :: I                        ! Loop counter
    integer :: STATUS                   ! Flag from deallocate
    ! Executable code
    if ( .not. associated ( binSelectors ) ) return
    do i = 1, size(binSelectors)
      call Deallocate_test ( binSelectors(i)%signals, &
        & 'binSelectors%signals', ModuleName )
      call Deallocate_test ( binSelectors(i)%sidebands, &
        & 'binSelectors%sidebands', ModuleName )
    end do
    deallocate ( binSelectors, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//"bin selectors binSelectors" )
  end subroutine DestroyBinSelectorDatabase

  ! ----------------------------------------------- DestroyL2PC ----
  subroutine DestroyL2PC ( l2pc )
    ! Dummy arguments
    type (Matrix_T), intent(inout), target :: L2PC

    integer :: QUANTITY                 ! Loop index
    integer :: VECTOR                   ! Loop index

    type (QuantityTemplate_T), pointer :: Qt ! Temporary pointer
    type (Vector_T), pointer :: V       ! Temporary pointer

    ! Exectuable code
    do vector = 1, 2
      if ( vector == 1 ) then
        v => l2pc%col%vec
      else
        v => l2pc%row%vec
      end if
      
      do quantity = 1, size(v%quantities)
        qt => v%quantities(quantity)%template
        call deallocate_test (qt%surfs, 'qt%surfs', ModuleName)
        call deallocate_test (qt%phi, 'qt%phi', ModuleName)
        call deallocate_test (v%quantities(quantity)%values, 'q%values',&
          & ModuleName)
      end do
      deallocate (v%quantities)
    end do
    
    ! Destory kStar
    call DestroyMatrix ( l2pc )
    
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
  end subroutine DestroyL2PCDatabase

  ! ------------------------------------ open_l2pc_file ------------
  subroutine Open_L2PC_File ( Filename, Lun )

    character(len=*), intent(in) :: Filename ! Name of the antenna pattern file
    integer, intent(out) :: Lun              ! Logical unit number to read it

    logical :: Exist, Opened
    integer :: Status

    do lun = 20, 99
      inquire ( unit=lun, exist=exist, opened=opened )
      if ( exist .and. .not. opened ) exit
    end do
    if ( opened .or. .not. exist ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "No logical unit numbers available" )
    open ( unit=lun, file=filename, status='old', form='formatted', &
      & access='sequential', iostat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "Unable to open l2pc file " // Filename )
  end subroutine Open_L2PC_File

  ! --------------------------------------------- OutputHDF5L2PC
  subroutine OutputHDF5L2PC ( filename, matrices, quantitiesNode, packed )
    character (len=*), intent(in) :: FILENAME
    type (Matrix_Database_T), dimension(:), pointer :: MATRICES
    integer, intent(in) :: QUANTITIESNODE
    logical, intent(in) :: PACKED

    ! Local variables
    integer :: FILEID                   ! ID of file
    integer :: STATUS                   ! From HDF
    integer :: FIELD                    ! Node index
    integer :: DB_INDEX                 ! Index of matrix
    type (Matrix_T), pointer :: tmpMatrix

    ! Executable code
    call H5Open_F ( status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open hdf5 system' )
    call H5FCreate_F ( trim(filename), H5F_ACC_TRUNC_F, fileID, &
      & status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open hdf5 l2pc file for output.' )
    do field = 2, nsons(quantitiesNode)
      db_index = decoration(decoration(subtree(field, quantitiesNode )))
      call GetFromMatrixDatabase ( matrices(db_index), tmpMatrix )
      call writeOneHDF5L2PC ( tmpMatrix, fileID, packed )
    end do ! in_field_no = 2, nsons(gson)
    call H5FClose_F ( fileID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName,&
      & 'Unable to close hdf5 l2pc file.' )
    call H5Close_F ( status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close hdf5 system' )

  end subroutine OutputHDF5L2PC

  ! ------------------------------------- Read_l2pc_file ------
  subroutine Read_l2pc_file ( Lun )
    use Trace_M, only: Trace_begin, Trace_end
    use Toggles, only: Toggle, gen
    ! Read all the bins in an l2pc file
    integer, intent(in) :: lun

    ! Local variables
    type (Matrix_T) :: L2pc
    integer :: Dummy
    logical :: Eof

    ! Executable code
    if ( toggle (gen) ) call trace_begin ( "Read_l2pc_file" )
    eof = .false.
    do while (.not. eof )
      call ReadOneASCIIL2PC ( l2pc, lun, eof )
      if (.not. eof) dummy = AddL2PCToDatabase ( l2pcDatabase, l2pc )
      if ( index ( switches, 'spa' ) /= 0 ) call Dump_struct ( l2pc, 'One l2pc bin' ) 

      ! Now nullify the pointers in l2pc so we don't clobber the one we've written
      nullify ( l2pc%block )
      nullify ( l2pc%row%nelts, l2pc%row%inst, l2pc%row%quant )
      nullify ( l2pc%col%nelts, l2pc%col%inst, l2pc%col%quant )
      nullify ( l2pc%row%vec%template%quantities, l2pc%col%vec%template%quantities )
      nullify ( l2pc%row%vec%quantities, l2pc%col%vec%quantities )
      
    end do

    if ( toggle (gen) ) call trace_end ( "Read_l2pc_file" )
  end subroutine Read_l2pc_file

  ! --------------------------------------- WriteOneHDF5L2PC -----------
  subroutine WriteOneHDF5L2PC ( L2pc, fileID, packed )
    ! This subroutine writes an l2pc to a file in hdf5 format

    ! Dummy arguments
    type (matrix_T), intent(in), target :: L2PC
    integer, intent(in) :: fileID
    logical, intent(in), optional :: PACKED

    ! Local variables
    integer :: BlockCol                 ! Index
    integer :: BlockRow                 ! Index
    integer :: BlocksGroupID            ! ID of group containing all blocks
    integer :: BlockGroupID             ! ID of this block group
    integer :: matrixID                 ! ID of hdf5 group containing matrix
    integer :: Quantity                 ! Loop counter
    integer :: Vector                   ! Loop counter
    integer :: AdjustedIndex            ! Index into quantities not skipped
    integer :: STATUS                   ! Error flag

    logical, dimension(l2pc%row%vec%template%noQuantities), target :: ROWPACK
    logical, dimension(l2pc%col%vec%template%noQuantities), target :: COLPACK
    logical, dimension(:), pointer :: THISPACK
    logical :: MYPACKED

    character ( len=32 ) :: NAME        ! A name for output

    type (MatrixElement_T), pointer :: M0 ! A Matrix0 within kStar
    type (QuantityTemplate_T), pointer :: Qt  ! Temporary pointers
    type (Vector_T), pointer :: V       ! Temporary pointer

    ! Executable code

    myPacked = .false.
    if ( present ( packed ) ) myPacked = packed

    ! Work out which quantities we can skip
    if ( myPacked ) then
      call MakeMatrixPackMap ( l2pc, rowPack, colPack )
    else
      rowPack = .true.
      colPack = .true.
    end if

    ! First create the group for this.
    call get_string ( l2pc%name, name, strip=.true. )
    call h5gCreate_f ( fileID, trim(name), matrixID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create group for l2pc matrix' )

    ! Now dump the vectors
    call WriteVectorAsHDF5 ( matrixID, l2pc%col%vec, 'Columns', colPack )
    call WriteVectorAsHDF5 ( matrixID, l2pc%row%vec, 'Rows', rowPack )

    ! Now create a group for the blocks
    call h5gCreate_f ( matrixID, 'Blocks', blocksGroupID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create group for all l2pc matrix blocks' )

    ! Make the flags
    call MakeHDF5Attribute ( blocksGroupID, 'rowInstanceFirst', l2pc%row%instFirst )
    call MakeHDF5Attribute ( blocksGroupID, 'colInstanceFirst', l2pc%col%instFirst )

    ! Now loop over the blocks and write them.
    do blockRow = 1, l2pc%row%NB
      do blockCol = 1, l2pc%col%NB
        ! Do we write this block?
        if ( rowPack(l2pc%row%quant(blockRow)) .and. &
          &  colPack(l2pc%col%quant(blockCol)) ) then
          ! Identify the block
          m0 => l2pc%block(blockRow, blockCol)
          ! Get a name for this group for the block
          write ( name, * ) 'Block', blockRow, blockCol
          ! Create a grop for this block
          call h5gCreate_f ( blocksGroupID, trim(name), blockGroupID, status )
          if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Unable to create group for l2pc matrix block' )
          ! Stick some attributes on it
          call get_string ( &
            & l2pc%row%vec%quantities(&
            &    l2pc%row%quant(blockRow))%template%name, name )
          call MakeHDF5Attribute ( blockGroupID, 'rowQuantity', trim(name) )
          call get_string ( &
            & l2pc%col%vec%quantities(&
            &    l2pc%col%quant(blockCol))%template%name, name )
          call MakeHDF5Attribute ( blockGroupID, 'colQuantity', trim(name) )
          call MakeHDF5Attribute ( blockGroupID, 'rowInstance', l2pc%row%inst(blockRow) )
          call MakeHDF5Attribute ( blockGroupID, 'colInstance', l2pc%col%inst(blockCol) )
          call MakeHDF5Attribute ( blockGroupID, 'kind', m0%kind )
          ! Write the datasets
          if ( m0%kind /= m_absent ) then
            call SaveAsHDF5DS ( blockGroupID, 'values', m0%values )
            if ( m0%kind /= m_full ) then
              call MakeHDF5Attribute ( blockGroupID, 'noValues', size(m0%values) )
              call SaveAsHDF5DS ( blockGroupID, 'r1', m0%r1 )
              call SaveAsHDF5DS ( blockGroupID, 'r2', m0%r2 )
            end if                      ! Not sparse or banded
          end if                        ! Not absent
          ! Close group for block
          call h5gClose_f ( blockGroupID, status )
          if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Unable to close group for l2pc matrix block' )
        end if                          ! Do this block
      end do
    end do

    ! Now close blocks group
    call h5gClose_f ( blocksGroupID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close qroup for all l2pc matrix blocks' )

    ! Close matrix group
    call h5gClose_f ( matrixID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close group for l2pc matrix' ) 
 end subroutine WriteOneHDF5L2PC
  
  ! --------------------------------------- WriteOneL2PC ---------------
  subroutine WriteOneL2PC ( L2pc, Unit, packed )
    ! This subroutine writes an l2pc to a file
    ! Currently this file is ascii, later it will be
    ! some kind of HDF file

    ! Dummy arguments
    type (matrix_T), intent(in), target :: L2pc
    integer, intent(in) :: Unit
    logical, intent(in), optional :: PACKED

    ! Local parameters
    character (len=*), parameter :: rFmt = "(4(2x,1pg15.8))"
    character (len=*), parameter :: iFmt = "(8(2x,i6))"

    ! Local variables
    integer :: BlockCol                 ! Index
    integer :: BlockRow                 ! Index
    integer :: Quantity                 ! Loop counter
    integer :: Vector                   ! Loop counter
    integer :: AdjustedIndex            ! Index into quantities not skipped

    logical, dimension(l2pc%row%vec%template%noQuantities), target :: ROWPACK
    logical, dimension(l2pc%col%vec%template%noQuantities), target :: COLPACK
    logical, dimension(:), pointer :: THISPACK
    logical :: MYPACKED

    character (len=132) :: Line, Word1, Word2 ! Line of text

    type (MatrixElement_T), pointer :: M0 ! A Matrix0 within kStar
    type (QuantityTemplate_T), pointer :: Qt  ! Temporary pointers
    type (Vector_T), pointer :: V       ! Temporary pointer

    ! Executable code

    myPacked = .false.
    if ( present ( packed ) ) myPacked = packed

    ! Work out which quantities we can skip
    if ( myPacked ) then
      call MakeMatrixPackMap ( l2pc, rowPack, colPack )
    else
      rowPack = .true.
      colPack = .true.
    end if

    ! First dump the xStar and yStar
    do vector = 1, 2
      ! Identify vector
      if ( vector == 1 ) then
        write (unit,*) 'xStar'
        v => l2pc%col%vec
        thisPack => colPack
      else
        write (unit,*) 'yStar'
        v => l2pc%row%vec
        thisPack => rowPack
      end if

      write (unit,*) count ( thisPack )
      ! Loop over quantities
      do quantity = 1, size(v%quantities)
        if ( thisPack(quantity) ) then
          qt => v%quantities(quantity)%template

          ! Write quantity name - will be ignored on the read (at least by
          ! fortran, IDL pays attention
          call get_string ( qt%name, line )
          write (unit,*) trim(line)
          
          ! Write quantity type
          call get_string ( lit_indices(qt%quantityType), line )
          write (unit,*) trim(line)
          
          ! Write other info associated with type
          select case ( qt%quantityType )
          case (l_vmr)
            call get_string ( lit_indices(qt%molecule), line )
            write (unit,*) trim(line)
          case (l_radiance)
            call GetSignalName ( qt%signal, line, sideband=qt%sideband )
            write (unit,*) trim(line)
          end select
          
          ! Write out the dimensions for the quantity and the edges
          write (unit,*) qt%noChans, qt%noSurfs, qt%noInstances,&
            &  'noChans, noSurfs, noInstances'
          call get_string ( lit_indices(qt%verticalCoordinate), line )
          write (unit,*) qt%coherent, qt%stacked, trim(line), &
            &  ' coherent, stacked, verticalCoordinate'
!           if ( all (qt%verticalCoordinate /= (/ l_none, l_zeta /)) &
!             & .and. (vector==1) .and. (qt%quantityType /= l_ptan) ) &
!             &   call MLSMessage(MLSMSG_Error,ModuleName, &
!             &     "Only zeta coordinates allowed (or none) for xStar.")
          write (unit,*) 'surfs'
          write (unit, rFmt) qt%surfs
          write (unit,*) 'phi'
          write (unit, rFmt) qt%phi
        end if
      end do                            ! First loop over quantities

      ! Now do a second loop and write the values
      adjustedIndex = 1
      do quantity = 1, count ( thisPack )
        if ( thisPack(quantity) ) then
          write (unit,*) 'values', adjustedIndex
          write (unit,rFmt) v%quantities(quantity)%values
          adjustedIndex = adjustedIndex + 1
        end if
      end do                            ! Second loop over quantities

    end do                              ! Loop over xStar/yStar
    
    ! Now dump kStar
    write (unit,*) 'kStar'
    write (unit,*) l2pc%row%instFirst, l2pc%col%instFirst, 'Instances first'
    do blockRow = 1, l2pc%row%NB
      do blockCol = 1, l2pc%col%NB
        ! Print the type of the matrix
        if ( rowPack(l2pc%row%quant(blockRow)) .and. &
          &  colPack(l2pc%col%quant(blockCol)) ) then
          m0 => l2pc%block(blockRow, blockCol)
          write (unit,*) blockRow, blockCol, m0%kind,&
            & 'row, col, kind'
          call get_string ( &
            & l2pc%row%vec%quantities(&
            &    l2pc%row%quant(blockRow))%template%name, word1 )
          call get_string ( &
            & l2pc%col%vec%quantities(&
            &    l2pc%col%quant(blockCol))%template%name, word2 )
          write (unit,*) trim(word1), l2pc%row%inst(blockRow), ' , ',&
            &            trim(word2), l2pc%col%inst(blockCol)
          select case (m0%kind)
          case (M_Absent)
          case (M_Banded, M_Column_sparse)
            write (unit,*) size(m0%values), ' no values'
            write (unit,*) 'R1'
            write (unit,iFmt) m0%R1
            write (unit,*) 'R2'
            write (unit,iFmt) m0%R2
            write (unit,*) 'values'
            write (unit,rFmt) m0%values
          case (M_Full)
            write (unit,rFmt) m0%values
          end select
        end if
      end do
    end do

  end subroutine WriteOneL2PC

  ! ======================================= PRIVATE PROCEDURES ====================

  ! ----------------------------------- MakeMatrixPackMap -----------
  subroutine MakeMatrixPackMap ( m, rowPack, colPack )
    ! This subroutine fills the boolean arrays rowPack, colPack
    ! (each length row/col%noQuantities) with a flag set true
    ! if the quantity has any derivatives at all
    
    ! Dummy arguments
    type (Matrix_T), intent(in) :: M
    logical, intent(out), dimension(M%row%vec%template%noQuantities) :: ROWPACK
    logical, intent(out), dimension(M%col%vec%template%noQuantities) :: COLPACK

    ! Local variables
    integer :: ROWQ                     ! Loop counter
    integer :: COLQ                     ! Loop counter
    integer :: ROWI                     ! Loop counter
    integer :: COLI                     ! Loop counter
    integer :: ROWBLOCK                 ! Block index
    integer :: COLBLOCK                 ! Block index

    ! Executable code

    rowPack = .false.
    colPack = .false.

    ! Do a nested loop over cols/rows
    ! I tried to be fancy with cycles etc. but the code got really messy.
    ! This simple approach is probably the clearest.
    do colQ = 1, m%col%vec%template%noQuantities
      do colI = 1, m%col%vec%quantities(colQ)%template%noInstances
        colBlock = FindBlock ( m%col, colQ, colI )
        do rowQ = 1, m%row%vec%template%noQuantities
          do rowI = 1, m%row%vec%quantities(rowQ)%template%noInstances
            rowBlock = FindBlock ( m%row, rowQ, rowI )
            if ( m%block ( rowBlock, colBlock ) % kind /= M_Absent ) then
              rowPack ( rowQ ) = .true.
              colPack ( colQ ) = .true.
            end if
          end do
        end do
      end do
    end do
  end subroutine MakeMatrixPackMap

  ! --------------------------------------- WriteL2PC ---------------
  subroutine ReadOneASCIIL2PC ( L2pc, Unit, Eof )
    ! This subroutine writes an l2pc to a file
    ! Currently this file is ascii, later it will be
    ! some kind of HDF file

    ! Dummy arguments
    type (Matrix_T), intent(out), target :: L2pc
    integer, intent(in) :: Unit
    logical, intent(inout) :: Eof

    ! Local variables
    integer :: BLOCKCOL                 ! Index
    integer :: BLOCKKIND                ! Kind of matrix block
    integer :: BLOCKROW                 ! Index
    integer :: NOVALUES                 ! For banded/sparse matrices
    integer :: TESTBLOCKROW             ! Test index
    integer :: TESTBLOCKCOL             ! Test index

    logical :: COLINSTFIRST             ! Matrix order
    logical :: ROWINSTFIRST             ! Matrix order

    character (len=132) :: Line         ! Line of text

    type (MatrixElement_T), pointer :: M0 ! A Matrix0 within kStar
    integer :: XSTAR                    ! Linearisation state vector index
    integer :: YSTAR                    ! Radiances for xStar vector index

    ! executable code

    ! First read the xStar and yStar
    call ReadOneVectorFromASCII ( unit, xStar, eof )
    if ( eof ) return

    call ReadOneVectorFromASCII ( unit, yStar, eof )
    if ( eof ) return

    ! Now read kStar
    read (unit,*) line                  ! Line saying kStar
    read (unit,*) rowInstFirst, colInstFirst ! Flags
    call CreateEmptyMatrix ( l2pc, 0, l2pcVs(yStar), l2pcVs(xStar), &
      & row_quan_first = .not. rowInstFirst,&
      & col_quan_first = .not. colInstFirst )

    ! Loop over blocks and read them
    do blockRow = 1, l2pc%row%NB
      do blockCol = 1, l2pc%col%NB
        ! Read the type of the matrix and a set of test indices
        read (unit,*) testBlockRow, testBlockCol, blockKind
        if (testBlockRow /= blockRow) call MLSMessage(MLSMSG_Error,ModuleName,&
          & 'Bad row number for kStar')
        if (testBlockCol /= blockCol) call MLSMessage(MLSMSG_Error,ModuleName,&
          & 'Bad col number for kStar')
        read (unit,*) line              ! String giving info, we can ignore.
        select case (blockKind)
        case (M_Absent)
          call CreateBlock ( l2pc, blockRow, blockCol, blockKind )
        case (M_Banded, M_Column_sparse)
          read (unit,*) noValues
          call CreateBlock ( l2pc, blockRow, blockCol, blockKind, noValues )
          m0 => l2pc%block ( blockRow, blockCol )
          read (unit,*) line ! 'R1'
          read (unit,*) m0%R1
          read (unit,*) line ! 'R2'
          read (unit,*) m0%R2
          read (unit,*) line ! 'values'
          read (unit,*) m0%values
        case (M_Full)
          call CreateBlock ( l2pc, blockRow, blockCol, blockKind )
          m0 => l2pc%block ( blockRow, blockCol )
          read (unit,*) m0%values
        end select
      end do
    end do
  end subroutine ReadOneASCIIL2PC

  ! ------------------------------------------ ReadOneVectorFromASCII ----
  subroutine ReadOneVectorFromASCII ( unit, vector, eof )
    ! Reads a vector from l2pc file and adds it to internal databases. This
    ! is internal as having it inside the above routine screws up databases.
    
    ! Dummy arguments
    integer, intent(in) :: UNIT         ! File unit
    integer, intent(out) :: VECTOR      ! Index of Vector read in L2PCVs
    logical, intent(out) :: EOF         ! Flag

    ! Local variables
    integer :: QUANTITY                 ! Loop counter
    integer :: SIDEBAND                 ! From parse signal
    integer :: STATUS                   ! Flag
    integer :: STRINGINDEX              ! Index of string
    integer :: NOQUANTITIES             ! Number of quantities in a vector
    integer :: NOINSTANCESOR1           ! For allocates
    integer :: NOSURFSOR1               ! For allocates
    integer :: VTINDEX                  ! Index for this vector template
    
    integer, dimension(:), pointer :: SIGINDS ! Result of parse signal
    integer, dimension(:), pointer :: QTINDS  ! Quantity indices
    
    character (len=132) :: Line       ! Line of text
    
    type (QuantityTemplate_T) :: QT     ! Temporary quantity template
    type (VectorTemplate_T) :: VT       ! Temporary template for vector
    type (Vector_T) :: V                ! Temporary vector
    type (Decls) :: Decl                ! From tree
    
    ! Executable code
    
    eof = .false. 
    nullify ( sigInds, qtInds )
    read (unit,*, IOSTAT=status) line
    if (status == -1 ) then
      eof = .true.
      return
    end if
    
    ! Note we fill this later
    read (unit,*) noQuantities
    call allocate_test ( qtInds, noQuantities, 'qtInds', ModuleName )
    
    ! Loop over quantities
    do quantity = 1, noQuantities
      
      ! Nullify stuff so we don't clobber arrays now in databases
      nullify ( qt%surfs, qt%phi )

      ! Read and ignore the name, we're going on the type, at least for the
      ! moment.
      read (unit,*, IOSTAT=status) line
      
      ! Read quantity type
      read (unit,*, IOSTAT=status) line
      if (status /= 0 ) then
        call io_error('io error in L2PC_m: ReadOneVector' // &
          & ' Fortran read of line for quantity type', status)
        call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'An io-error occured' )
      end if
      stringIndex = enter_terminal ( trim(line), t_identifier )
      decl = get_decl ( stringIndex, type=enum_value )
      qt%quantityType = decl%units
      qt%name = stringIndex

      ! Set defaults for coordinates, radiance is the later exception
      qt%verticalCoordinate = l_zeta
      qt%frequencyCoordinate = l_none

      ! Read other info associated with type
      select case ( qt%quantityType )
      case (l_vmr)
        read (unit,*, IOSTAT=status) line
        if (status /= 0 ) then
        call io_error('io error in L2PC_m: ReadOneVector' // &
          & ' Fortran read of line for l_vmr', status)
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'An io-error occured' )
        end if
        stringIndex = enter_terminal ( trim(line), t_identifier )
        decl = get_decl ( stringIndex, type=enum_value )
        qt%molecule = decl%units
        qt%name = stringIndex
      case (l_radiance)
        read (unit,*, IOSTAT=status) line
        if (status /= 0 ) then
        call io_error('io error in L2PC_m: ReadOneVector' // &
          & ' Fortran read of line for radiance', status)
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'An io-error occured' )
        end if
        call Parse_Signal (line, sigInds, sideband=sideband)
        qt%signal = sigInds(1)
        qt%sideband = sideband
        qt%frequencyCoordinate = l_channel
        qt%verticalCoordinate = l_geodAltitude
        call deallocate_test(sigInds,'sigInds',ModuleName)
      case default
      end select
      
      ! Next read the dimensions for the quantity
      read (unit,*, IOSTAT=status) qt%noChans, qt%noSurfs, qt%noInstances
      if (status /= 0 ) then
        call io_error('io error in L2PC_m: ReadOneVector' // &
          & ' Fortran read of NoChans, noSurfs, noInstances', status)
          print *, 'Quantity: ', quantity
          print *, 'QuantityName: ', qt%name
          print *, 'QuantityType: ', qt%quantityType
          print *, 'l_vmr: ', l_vmr
          print *, 'last literal read: ', trim(line)
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'An io-error occured' )
      end if
      qt%instanceLen = qt%noChans* qt%noSurfs
      read (unit,*, IOSTAT=status) qt%coherent, qt%stacked
      if (status /= 0 ) then
        call io_error('io error in L2PC_m: ReadOneVector' // &
          & ' Fortran read of coherent, stacked', status)
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'An io-error occured' )
      end if
      
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
      
      call Allocate_test( qt%surfs, qt%noSurfs, noInstancesOr1, 'qt%surfs', ModuleName)
      call Allocate_test( qt%phi, noSurfsOr1, qt%noInstances, 'qt%phi', ModuleName)
      
      read (unit,*, IOSTAT=status) line              ! Line saying surfs
      if (status /= 0 ) then
        call io_error('io error in L2PC_m: ReadOneVector' // &
          & ' Fortran read of line saying surfs', status)
        call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'An io-error occured' )
      end if
      read (unit,*, IOSTAT=status) qt%surfs
      if (status /= 0 ) then
        call io_error('io error in L2PC_m: ReadOneVector' // &
          & ' Fortran read of surfs', status)
        call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'An io-error occured' )
      end if
      read (unit,*, IOSTAT=status) line              ! Line saying phi
      if (status /= 0 ) then
        call io_error('io error in L2PC_m: ReadOneVector' // &
          & ' Fortran read of line saying phi', status)
        call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'An io-error occured' )
      end if
      read (unit,*, IOSTAT=status) qt%phi
      if (status /= 0 ) then
        call io_error('io error in L2PC_m: ReadOneVector' // &
          & ' Fortran read of phi', status)
        call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'An io-error occured' )
      end if
      
      ! Now add this template to our private database 
      qt%id = l2pcQTCounter
      l2pcQTCounter = l2pcQtCounter + 1
      qtInds(quantity) = AddQuantityTemplateToDatabase ( l2pcQTs, qt )
      
    end do                            ! First Loop over quantities
    
    ! Now create a vector template with these quantities
    call ConstructVectorTemplate ( 0, l2pcQTs, qtInds, vt )
    vt%id = l2pcVTCounter
    l2pcVtCounter = l2pcVtCounter + 1
    vtIndex = AddVectorTemplateToDatabase ( l2pcVTs, vt )
    
    call deallocate_test ( qtInds, 'qtInds', ModuleName )
    
    ! Now create a vector for this vector template
    v = CreateVector ( 0, l2pcVTs(vtIndex), l2pcQTs, vectorNameText='_v' )
    vector = AddVectorToDatabase ( l2pcVs, v )
    
    ! Now go through the quantities again and read the values
    do quantity = 1, noQuantities
      read (unit,*, IOSTAT=status) line              ! Just a comment line
      if (status /= 0 ) then
        call io_error('io error in L2PC_m: ReadOneVector' // &
          & ' Fortran read of comment in loop of quantities', status)
        call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'An io-error occured' )
      end if
      read (unit,*, IOSTAT=status) l2pcVs(vector)%quantities(quantity)%values
      if (status /= 0 ) then
        call io_error('io error in L2PC_m: ReadOneVector' // &
          & ' Fortran read of quantities', status)
        call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'An io-error occured' )
      end if
    end do
  end subroutine ReadOneVectorFromASCII

  ! --------------------------------------- ReadOneHDF5L2PC ------------
  subroutine ReadOneHDF5L2PCRecord ( l2pc, fileID, index )
    type ( Matrix_T ), intent(out), target :: L2PC
    integer, intent(in) :: FILEID       ! HDF5 ID of input file
    integer, intent(in) :: INDEX        ! Index of l2pc entry to read

    ! Local variables
    integer :: BLOCKCOL                 ! Loop counter
    integer :: BLOCKID                  ! Id of block group (one block at a time)
    integer :: BLOCKROW                 ! Loop counter
    integer :: BLOCKSID                 ! Id of blocks group
    integer :: KIND                     ! Kind of block (absent, etc.)
    integer :: MATRIXID                 ! HDF5 for matrix group
    integer :: NOVALUES                 ! For banded or column sparse cases
    integer :: OBJTYPE                  ! From HDF5
    integer :: STATUS                   ! Flag from HDF5
    integer :: XSTAR                    ! Linearisation state vector index
    integer :: YSTAR                    ! Radiances for xStar vector index

    logical :: COLINSTANCEFIRST         ! Flag for matrix
    logical :: ROWINSTANCEFIRST         ! Flag for matrix

    type (MatrixElement_T), pointer :: M0 ! A Matrix0 within kStar
    character ( len=64 ) :: MATRIXNAME  ! Name for matrix
    character ( len=64 ) :: NAME        ! Name for block group

    ! Executable code
    call h5gGet_obj_info_idx_f ( fileID, '/', index, matrixName, &
      & objType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get information on matrix in l2pc file' )
    call h5gOpen_f ( fileID, trim(matrixName), matrixId, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open matrix in l2pc file' )

    ! Read the row and column vectors
    call ReadOneVectorFromHDF5 ( matrixId, 'Columns', xStar )
    call ReadOneVectorFromHDF5 ( matrixId, 'Rows', yStar )

    ! Get the instance first information
    call h5gOpen_f ( matrixID, 'Blocks', blocksID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to access Blocks group of matrix' )
    call GetHDF5Attribute ( blocksID, 'rowInstanceFirst', rowInstanceFirst )
    call GetHDF5Attribute ( blocksID, 'colInstanceFirst', colInstanceFirst )

    ! Create the matrix
    call CreateEmptyMatrix ( l2pc, 0, l2pcVs(yStar), l2pcVs(xStar), &
      & row_quan_first = .not. rowInstanceFirst,&
      & col_quan_first = .not. colInstanceFirst )

    ! Loop over blocks and read them
    do blockRow = 1, l2pc%row%NB
      do blockCol = 1, l2pc%col%NB
        ! Access this block
        write ( name, * ) 'Block', blockRow, blockCol
        call h5gOpen_f ( blocksId, trim(name), blockId, status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to open group for l2pc matrix block '//trim(name) )
        ! Could check it's the block we're expecting but I think I'll be lazy
        call GetHDF5Attribute ( blockID, 'kind', kind )
        if ( kind == m_banded .or. kind == m_column_sparse ) then
          call GetHDF5Attribute ( blockID, 'noValues', noValues )
          call CreateBlock ( l2pc, blockRow, blockCol, kind, noValues )
          m0 => l2pc%block ( blockRow, blockCol )
          call LoadFromHDF5DS ( blockId, 'r1', m0%r1 )
          call LoadFromHDF5DS ( blockId, 'r1', m0%r2 )

        else
          call CreateBlock ( l2pc, blockRow, blockCol, kind )
          m0 => l2pc%block ( blockRow, blockCol )
        end if
        if ( kind /= m_absent ) &
          call LoadFromHDF5DS ( blockID, 'values', m0%values )
        call h5gClose_f ( blockId, status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to close group for l2pc matrix block '//trim(name) )
      end do
    end do
    
    ! Finish up
    call h5gClose_f ( blocksID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close Blocks group' )
    call h5gClose_f ( matrixID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close matrix group' )

  end subroutine ReadOneHDF5L2PCRecord

  ! --------------------------------------- ReadOneVectorFromHDF5 ------
  subroutine ReadOneVectorFromHDF5 ( location, name, vector )
    ! Read a vector from an l2pc HDF5 and adds it to internal databases.
    ! Dummy arguments
    integer, intent(in) :: LOCATION     ! Node in HDF5
    character (len=*), intent(in) :: NAME ! Name of vector
    integer, intent(out) :: VECTOR      ! Index of vector read in L2PCVectors

    ! Local variables
    integer :: QUANTITY                 ! Loop inductor
    integer :: NOQUANTITIES             ! Number of quantities in vector
    integer :: VID                      ! HDF5 ID of vector group
    integer :: QID                      ! HDF5 ID of quantity group
    integer :: OBJTYPE                  ! Irrelevant argument to HDF5
    integer :: INDEX                    ! Index into various arrays
    integer :: STATUS                   ! Flag from HDF
    integer :: STRINGINDEX              ! A string index
    integer :: SIDEBAND                 ! Sideband -1,0,1
    integer :: NOSURFSOR1               ! Dimension
    integer :: NOINSTANCESOR1           ! Dimension
    integer :: VTINDEX                  ! Index of vector template
    character (len=64) :: THISNAME      ! Name of this quantity
    character (len=64) :: WORD          ! General string

    integer, dimension(:), pointer :: SIGINDS ! Index into signals database
    integer, dimension(:), pointer :: QTINDS ! Quantity indices

    character (len=64), pointer, dimension(:) :: QUANTITYNAMES ! Names of quantities
    type ( QuantityTemplate_T) :: QT    ! Template for the quantity
    type ( VectorTemplate_T) :: VT    ! Template for the vector
    type ( Vector_T ) :: V              ! The vector
    type ( Decls ) :: DECL              ! From tree


    ! Executable code
    nullify ( sigInds, qtInds )

    ! Open the vector group
    call h5gOpen_f ( location, name, vId, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open group for vector '//trim(name) )

    ! Get the number of quantities
    call GetHDF5Attribute ( vId, 'noQuantities', noQuantities )
    ! Work out the order of the quantities
    nullify ( quantityNames )
    call Allocate_test ( quantityNames, noQuantities, 'quantityNames', ModuleName )
    do quantity = 1, noQuantities
      call h5gget_obj_info_idx_f ( location, name, quantity, thisName, &
        & objType, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to acces a quantity within '//trim(name) )
      call h5gOpen_f ( vId, trim(thisName), qId, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open quantity '//trim(thisName)//' in vector '//trim(name) )
      call GetHDF5Attribute ( qId, 'index', index )
      quantityNames ( index ) = thisName
      call h5gClose_f ( qId, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close quantity '//trim(thisName)//' in vector '//trim(name) )
    end do

    ! Now go through quantities in order
    do quantity = 1, noQuantities
      ! Nullify stuff so we don't clobber our databases
      nullify ( qt%surfs, qt%phi )

      ! Point to this quanity
      call h5gOpen_f ( vId, trim(quantityNames(quantity)), qId, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open quantity '//trim(quantityNames(quantity))// &
        & ' in vector '//trim(name) )

      ! Store the name
      stringIndex = enter_terminal ( trim(quantityNames(quantity)), t_identifier ) 
      qt%name = stringIndex

      ! Get the quantity type
      call GetHDF5Attribute ( qId, 'type', word )
      stringIndex = enter_terminal ( trim(word), t_identifier ) 
      decl = get_decl ( stringIndex, type=enum_value )
      qt%quantityType = decl%units

      ! Get other info as appropriate
      select case ( qt%quantityType )
      case ( l_vmr )
        call GetHDF5Attribute ( qID, 'molecule', word )
        stringIndex = enter_terminal ( trim(word), t_identifier ) 
        decl = get_decl ( stringIndex, type=enum_value )
        qt%molecule = decl%units
      case ( l_radiance )
        call GetHDF5Attribute ( qID, 'signal', word )
        call Parse_Signal ( word, sigInds, sideband=sideband)
        qt%signal = sigInds(1)
        qt%sideband = sideband
        qt%frequencyCoordinate = l_channel
        qt%verticalCoordinate = l_geodAltitude
        call deallocate_test(sigInds,'sigInds',ModuleName)
      case default
      end select

      ! Now read the dimensions fo rhte quantity
      call GetHDF5Attribute ( qID, 'noInstances', qt%noInstances )
      call GetHDF5Attribute ( qID, 'noSurfs', qt%noSurfs )
      call GetHDF5Attribute ( qID, 'noChans', qt%noChans )
      call GetHDF5Attribute ( qId, 'verticalCoordinate', word )
      stringIndex = enter_terminal ( trim(word), t_identifier ) 
      decl = get_decl ( stringIndex, type=enum_value )
      qt%verticalCoordinate = decl%units
      call GetHDF5Attribute ( qId, 'coherent', qt%coherent )
      call GetHDF5Attribute ( qId, 'stacked', qt%stacked )

      ! This creates a rather minimal quantity template, it seems to be, but
      ! I'm modeling it after the ASCII version.  I wonder why I didn't use
      ! ConstructQuantityTemplte? 
      
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
      
      call Allocate_test( qt%surfs, qt%noSurfs, noInstancesOr1, 'qt%surfs', ModuleName)
      call Allocate_test( qt%phi, noSurfsOr1, qt%noInstances, 'qt%phi', ModuleName)

      ! Get the surfaces and phis
      call LoadFromHDF5DS ( qId, 'surfs', qt%surfs )
      call LoadFromHDF5DS ( qId, 'phi', qt%phi )

      ! Now add this template to our private database 
      qt%id = l2pcQTCounter
      l2pcQTCounter = l2pcQtCounter + 1
      qtInds(quantity) = AddQuantityTemplateToDatabase ( l2pcQTs, qt )

      ! For the moment, close the quantity. We'll come back to it later to fill
      ! up the values for the vector
      call h5gClose_f ( qId, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close quantity '//trim(thisName)//' in vector '//trim(name) )
    end do

    ! Now create a vector template with these quantities
    call ConstructVectorTemplate ( 0, l2pcQTs, qtInds, vt )
    vt%id = l2pcVTCounter
    l2pcVtCounter = l2pcVtCounter + 1
    vtIndex = AddVectorTemplateToDatabase ( l2pcVTs, vt )
    
    call deallocate_test ( qtInds, 'qtInds', ModuleName )
    
    ! Now create a vector for this vector template
    v = CreateVector ( 0, l2pcVTs(vtIndex), l2pcQTs, vectorNameText='_v' )
    vector = AddVectorToDatabase ( l2pcVs, v )
    
    ! Now go through the quantities again and read the values
    do quantity = 1, noQuantities
      call h5gOpen_f ( vId, trim(quantityNames(quantity)), qId, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open quantity '//trim(quantityNames(quantity))// &
        & ' in vector '//trim(name) )
      call LoadFromHDF5DS ( qId, 'values', &
        & l2pcVs(vector)%quantities(quantity)%values )
      call h5gClose_f ( qId, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close quantity '//trim(thisName)//' in vector '//trim(name) )
    end do

    call Deallocate_test ( quantityNames, 'quantityNames', ModuleName )
      
  end subroutine ReadOneVectorFromHDF5

  ! --------------------------------------- WriteVectorAsHDF5 ----------
  subroutine WriteVectorAsHDF5 ( location, vector, name, packInfo )
    integer, intent(in) :: LOCATION     ! The HDF5 location for the vector
    type (Vector_T), intent(in) :: VECTOR ! The vector to write
    character(len=*), intent(in) :: NAME ! Name for vector
    logical, dimension(:), intent(in) :: PACKINFO ! Pack vector if set

    ! Local variables
    character ( len=132 ) :: LINE       ! A line of text
    character ( len=32 ) :: QNAME       ! Name of quantity
    integer :: STATUS                   ! Flag from HDF5
    integer :: QUANTITY                 ! Index
    integer :: QID                      ! HDF5 ID for quantity group
    integer :: VID                      ! HDF5 ID for vector group
    type (QuantityTemplate_T), pointer :: QT ! The template

    ! Executable code
    ! Create a group for the vector
    call h5gCreate_f ( location, trim(name), vID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create group for vector ' // trim(name) )
    call MakeHDF5Attribute ( vID, 'noQuantities', count(packInfo) )
    do quantity = 1, size(vector%quantities )
      if ( packInfo(quantity) ) then
        qt => vector%quantities(quantity)%template
        ! Create a group for the quantity
        call get_string ( qt%name, qName )
        call h5gCreate_f ( vID, trim(qName), qID, status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to create group for quantity ' // trim(qName) )

        ! Write quantity type and index
        call MakeHDF5Attribute ( qID, 'index', quantity )
        call get_string ( lit_indices(qt%quantityType), line )
        call MakeHDF5Attribute ( qID, 'type', trim(line) )
        ! Write other info as appropriate
        select case ( qt%quantityType )
        case (l_vmr)
          call get_string ( lit_indices(qt%molecule), line )
          call MakeHDF5Attribute ( qID, 'molecule', trim(line) )
        case (l_radiance)
          call GetSignalName ( qt%signal, line, sideband=qt%sideband )
          call MakeHDF5Attribute ( qID, 'signal', trim(line) )
        end select
        ! Write out the dimensions for the quantity
        call MakeHDF5Attribute ( qID, 'noChans', qt%noChans )
        call MakeHDF5Attribute ( qID, 'noSurfs', qt%noSurfs )
        call MakeHDF5Attribute ( qID, 'noInstances', qt%noInstances )
        call get_string ( lit_indices(qt%verticalCoordinate), line )
        call MakeHDF5Attribute ( qID, 'verticalCoordinate', trim(line) )
        call MakeHDF5Attribute ( qID, 'coherent', qt%coherent )
        call MakeHDF5Attribute ( qID, 'stacked', qt%stacked )
        ! Write out important coordinates
        call SaveAsHDF5DS ( qID, 'phi', qt%phi )
        call SaveAsHDF5DS ( qID, 'surfs', qt%surfs )
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

  end subroutine WriteVectorAsHDF5


end module L2PC_m

! $Log$
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

