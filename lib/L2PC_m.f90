! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module L2PC_m
  !=============================================================================

  ! This module contains data types etc. for dealing with the new EMLS L2PC
  ! files.  The first version will deal with ascii files, but later versions
  ! will probably be HDF.

  use Allocate_Deallocate, only: Allocate_test, Deallocate_test
  use Declaration_Table, only: DECLS, ENUM_VALUE, GET_DECL
  use Intrinsic, only: Lit_Indices, l_zeta, l_none
  use MLSCommon, only: R8
  use VectorsModule, only: DESTROYVECTORINFO, VECTOR_T, VECTORVALUE_T
  use MatrixModule_1, only: CREATEBLOCK, CREATEEMPTYMATRIX, DESTROYMATRIX, MATRIX_T
  use MatrixModule_0, only: M_ABSENT, M_BANDED, M_COLUMN_SPARSE, M_FULL, &
    & MATRIXELEMENT_T
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR, &
    & MLSMSG_ALLOCATE, MLSMSG_DEALLOCATE
  use Parse_Signal_m, only: Parse_Signal
  use QuantityTemplates, only: QUANTITYTEMPLATE_T
  use Intrinsic, only: L_RADIANCE, L_TEMPERATURE, L_VMR
  use String_Table, only: GET_STRING, DISPLAY_STRING
  use MLSSignals_m, only: GETSIGNALNAME
  use SYMBOL_TABLE, only: ENTER_TERMINAL
  use SYMBOL_TYPES, only: T_IDENTIFIER
  use Tree, only: DECORATION

  implicit NONE
  private
  
  public :: AddL2PCToDatabase, DestroyL2PC, DestroyL2PCDatabase, WriteOneL2PC
  public :: Open_l2pc_file, read_l2pc_file, close_l2pc_file

  ! Public types
  type, public :: l2pc_T
    logical :: readFromFile=.false.     ! See below
    type (Vector_T) :: xStar            ! The linearisation x
    type (Vector_T) :: yStar            ! The corresponding y
    type (Matrix_T) :: kStar            ! The jacobian matrix
  end type l2pc_T
  ! The read from file stuff is needed, because if it was created by hand the
  ! software will deallocate all the arrays associated with xStar, yStar and
  ! kStar when destory vectors etc. are called.  If on the other hand it was
  ! read from a file, we'll have to destroy these ourselves

  ! The l2pc database
  type(l2pc_T), dimension(:), pointer, public :: L2PCDatabase => NULL()

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains ! ============= Public Procedures ==========================

  ! --------------------------------------- WriteOneL2PC ---------------
  subroutine WriteOneL2PC ( l2pc, unit )
    ! This subroutine writes an l2pc to a file
    ! Currently this file is ascii, later it will be
    ! some kind of HDF file

    ! Dummy arguments
    type (l2pc_T), intent(in), target :: l2pc
    integer, intent(in) :: unit

    ! Local parameters
    character (len=*), parameter :: rFmt = "(4(2x,1pg15.8))"
    character (len=*), parameter :: iFmt = "(8(2x,i6))"


    ! Local variables
    integer :: blockRow                 ! Index
    integer :: blockCol                 ! Index
    integer :: quantity                 ! Loop counter
    integer :: instance                 ! Loop counter
    integer :: vector                   ! Loop counter

    type (QuantityTemplate_T), pointer :: qt  ! Temporary pointers
    type (Vector_T), pointer :: v       ! Temporary pointer
    type (MatrixElement_T), pointer :: m0 ! A Matrix0 within kStar

    character (len=132) :: line,word1,word2 ! Line of text

    ! executable code

    ! First dump the xStar and yStar
    do vector = 1, 2
      ! Identify vector
      if ( vector == 1 ) then
        write (unit,*) 'xStar'
        v => l2pc%xStar
      else
        write (unit,*) 'yStar'
        v => l2pc%yStar
      end if

      write (unit,*) size(v%quantities)
      ! Loop over quantities
      do quantity = 1, size(v%quantities)
        qt => v%quantities(quantity)%template

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
        write (unit,*) qt%coherent, qt%stacked, &
          &  'coherent, stacked'
        if ( any (qt%verticalCoordinate /= (/ l_none, l_zeta /)) &
          & .and. (vector==1) ) &
          & call MLSMessage(MLSMSG_Error,ModuleName, &
          & "Only zeta coordinates allowed (or none) for xStar.")
        write (unit,*) 'surfs'
        write (unit, rFmt) qt%surfs
        write (unit,*) 'phi'
        write (unit, rFmt) qt%phi

        ! Write the values
        write (unit,*) 'values'
        write (unit,rFmt) v%quantities(quantity)%values

      end do                            ! Loop over quantities
    end do                              ! Loop over xStar/yStar

    ! Now dump kStar
    write (unit,*) 'kStar'
    write (unit,*) l2pc%kStar%row%instFirst, l2pc%kStar%col%instFirst, 'Instances first'
    do blockRow = 1, l2pc%kStar%row%NB
      do blockCol = 1, l2pc%kStar%col%NB
        ! Print the type of the matrix
        m0 => l2pc%kStar%block(blockRow, blockCol)
        write (unit,*) blockRow, blockCol, m0%kind,&
          & 'row, col, kind'
        call get_string ( &
          & l2pc%kStar%row%vec%quantities(&
          &    l2pc%kStar%row%quant(blockRow))%template%name, word1 )
        call get_string ( &
          & l2pc%kStar%col%vec%quantities(&
          &    l2pc%kStar%col%quant(blockCol))%template%name, word2 )
        write (unit,*) trim(word1), l2pc%kStar%row%inst(blockRow), ' , ',&
          &            trim(word2), l2pc%kStar%col%inst(blockCol)
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
      end do
    end do

  end subroutine WriteOneL2PC
  
  ! --------------------------------------- WriteL2PC ---------------
  subroutine ReadOneL2PC ( l2pc, unit, eof )
    ! This subroutine writes an l2pc to a file
    ! Currently this file is ascii, later it will be
    ! some kind of HDF file

    ! Dummy arguments
    type (l2pc_T), intent(out), target :: l2pc
    integer, intent(in) :: unit
    logical, intent(inout) :: eof

    ! Local variables
    integer :: BLOCKCOL                 ! Index
    integer :: BLOCKROW                 ! Index
    integer :: BLOCKKIND                ! Kind of matrix block
    integer :: INSTANCE                 ! Loop counter
    integer :: NOINSTANCESOR1           ! For allocates
    integer :: NOSURFSOR1               ! For allocates
    integer :: NOQUANTITIES             ! Number of quantities in a vector
    integer :: NOVALUES                 ! For banded/sparse matrices
    integer :: QUANTITY                 ! Loop counter
    integer :: SIDEBAND                 ! From parse signal
    integer :: STATUS                   ! Flag
    integer :: STRINGINDEX              ! Index of string
    integer :: TESTBLOCKROW             ! Test index
    integer :: TESTBLOCKCOL             ! Test index
    integer :: VECTOR                   ! Loop counter

    integer, dimension(:), pointer :: SIGINDS ! Result of parse signal

    logical :: ROWINSTFIRST             ! Matrix order
    logical :: COLINSTFIRST             ! Matrix order

    type (QuantityTemplate_T), pointer :: qt  ! Temporary pointers
    type (VectorValue_T), pointer :: vv ! Temporary pointer
    type (Vector_T), pointer :: v       ! Temporary pointer
    type (MatrixElement_T), pointer :: m0 ! A Matrix0 within kStar
    type (Decls) :: decl                ! From the declaration table

    character (len=132) :: line         ! Line of text

    ! executable code

    nullify ( sigInds )
    l2pc%readFromFile = .true.

    ! First read the xStar and yStar
    do vector = 1, 2
      ! Identify vector
      if ( vector == 1 ) then
        read (unit,*, IOSTAT=status) line              ! Comment line
        if (status == -1 ) then
          eof = .true.
          return
        end if
        v => l2pc%xStar
      else
        read (unit,*) line              ! Comment line
        v => l2pc%yStar
      end if

      ! Note we fill this later

      read (unit,*) noQuantities
      allocate (v%quantities(noQuantities), STAT=status)
      if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName,&
        & MLSMSG_Allocate//"v%quantities")

      ! Loop over quantities
      do quantity = 1, noQuantities

        ! Create a template for this quantity
        allocate (v%quantities(quantity)%template, STAT=status)
        if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName,&
          & MLSMSG_Allocate//"v%quantities(:)%template")
        vv => v%quantities(quantity)
        qt => vv%template

        ! Read quantity type
        read (unit,*) line
        stringIndex = enter_terminal ( trim(line), t_identifier )
        decl = get_decl ( stringIndex, type=enum_value )
        qt%quantityType = decl%units

        ! Write other info associated with type
        select case ( qt%quantityType )
        case (l_vmr)
          read (unit,*) line
          stringIndex = enter_terminal ( trim(line), t_identifier )
          decl = get_decl ( stringIndex, type=enum_value )
          qt%molecule = decl%units
        case (l_radiance)
          read (unit,*) line
          call Parse_Signal (line, sigInds, sideband=sideband)
          qt%signal = sigInds(1)
          qt%sideband = sideband
          call deallocate_test(sigInds,'sigInds',ModuleName)
        case default
        end select

        ! Next read the dimensions for the quantity
        read (unit,*) qt%noChans, qt%noSurfs, qt%noInstances
        qt%instanceLen = qt%noChans* qt%noSurfs
        read (unit,*) qt%coherent, qt%stacked

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
        call Allocate_test( vv%values, qt%instanceLen, qt%noInstances,&
          & 'vv%values', ModuleName )

        read (unit,*) line              ! Line saying surfs
        read (unit,*) qt%surfs
        read (unit,*) line              ! Line saying phi
        read (unit,*) qt%phi

        ! Read the values
        read (unit,*) line              ! Line saying values
        read (unit,*) vv%values

      end do                            ! Loop over quantities

      ! Create a dummy template for this vector
      allocate (v%template, STAT=status)
      if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName,&
        & MLSMSG_Allocate//"v%template")
      v%template%noQuantities = noQuantities
      v%template%totalInstances = 0
      v%template%totalElements = 0
      do quantity = 1, noQuantities
        v%template%totalInstances = v%template%totalInstances + &
          & v%quantities(quantity)%template%noInstances
        v%template%totalElements = v%template%totalElements + &
          & v%quantities(quantity)%template%noInstances * &
          & v%quantities(quantity)%template%instanceLen
      end do
      ! Ignore the quantities stuff for the moment.
      ! If later versions need it we'll put something here
    end do                              ! Loop over xStar/yStar

    ! Now read kStar
    read (unit,*) line                  ! Line saying kStar
    read (unit,*) rowInstFirst, colInstFirst ! Flags
    call CreateEmptyMatrix ( l2pc%kStar, 0, l2pc%yStar, l2pc%xStar,&
      & row_quan_first = .not. rowInstFirst,&
      & col_quan_first = .not. colInstFirst )
    ! Loop over blocks and read them
    do blockRow = 1, l2pc%kStar%row%NB
      do blockCol = 1, l2pc%kStar%col%NB
        ! Read the type of the matrix and a set of test indices
        read (unit,*) testBlockRow, testBlockCol, blockKind
        if (testBlockRow /= blockRow) call MLSMessage(MLSMSG_Error,ModuleName,&
          & 'Bad row number for kStar')
        if (testBlockCol /= blockCol) call MLSMessage(MLSMSG_Error,ModuleName,&
          & 'Bad col number for kStar')
        read (unit,*) line              ! String giving info, we can ignore.
        select case (blockKind)
        case (M_Absent)
        case (M_Banded, M_Column_sparse)
          read (unit,*) noValues
          call CreateBlock ( l2pc%kStar, blockRow, blockCol, blockKind, noValues )
          m0 => l2pc%kStar%block ( blockRow, blockCol )
          read (unit,*) line ! 'R1'
          read (unit,*) m0%R1
          read (unit,*) line ! 'R2'
          read (unit,*) m0%R2
          read (unit,*) line ! 'values'
          read (unit,*) m0%values
        case (M_Full)
          call CreateBlock ( l2pc%kStar, blockRow, blockCol, blockKind )
          m0 => l2pc%kStar%block ( blockRow, blockCol )
          read (unit,*) m0%values
        end select
      end do
    end do

  end subroutine ReadOneL2PC

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

  ! -----------------------------------  Close_L2PC_File  -----
  subroutine Close_L2PC_File ( Lun )
    integer, intent(in) :: lun
    close ( lun )
  end subroutine Close_L2PC_File

  ! ------------------------------------- Read_l2pc_file ------
  subroutine Read_l2pc_file ( lun )
    use Trace_M, only: Trace_begin, Trace_end
    use Toggles, only: Toggle, gen
    ! Read all the bins in an l2pc file
    integer, intent(in) :: lun

    ! Local variables
    type (l2pc_T) :: l2pc
    integer :: dummy
    logical :: eof

    ! Executable code
    if ( toggle (gen) ) call trace_begin ( "Read_l2pc_file" )
    eof = .false.
    do while (.not. eof )
      call ReadOneL2PC ( l2pc, lun, eof )
      dummy = AddL2PCToDatabase ( l2pcDatabase, l2pc )
    end do
    if ( toggle (gen) ) call trace_end ( "Read_l2pc_file" )
  end subroutine Read_l2pc_file

  ! ------------------------------------  Add l2pc  to database ----
  integer function AddL2PCToDatabase ( Database, Item )
    
    ! This function simply adds an l2pc  to a database of said l2pc s.
    
    type(l2pc_T), dimension(:), pointer :: Database
    type(l2pc_T) :: Item
    
    type(l2pc_T), dimension(:), pointer :: TempDatabase
    
    include "addItemToDatabase.f9h"

    AddL2PCToDatabase = newSize
  end function AddL2PCToDatabase

  ! ----------------------------------------------- DestroyL2PC ----
  subroutine DestroyL2PC ( l2pc )
    ! Dummy arguments
    type (l2pc_T), intent(inout), target :: L2PC

    integer :: QUANTITY                 ! Loop index
    integer :: VECTOR                   ! Loop index

    type (Vector_T), pointer :: v       ! Temporary pointer
    type (QuantityTemplate_T), pointer :: qt ! Temporary pointer

    ! Exectuable code
    ! If this wasn't created from a file we don't have to do anything,
    ! as the main program will call destroy vector info etc.
    ! if we did read it from a file, we need to do a bit of work but not much
    if ( l2pc%readFromFile ) then
      do vector = 1, 2
        if ( vector == 1 ) then
          v => l2pc%xStar
        else
          v => l2pc%yStar
        end if
        
        do quantity = 1, size(v%quantities)
          qt => v%quantities(quantity)%template
          call deallocate_test (qt%surfs, 'qt%surfs', ModuleName)
          call deallocate_test (qt%phi, 'qt%phi', ModuleName)
          call deallocate_test (v%quantities(quantity)%values, 'q%values',&
            & ModuleName)
        end do
        deallocate (v%template)
        deallocate (v%quantities)
      end do
      
      ! Destory kStar
      call DestroyMatrix ( l2pc%kStar)

    end if
  end subroutine DestroyL2PC

  ! ------------------------------------------- DestroyL2PCDatabase ---
  subroutine DestroyL2PCDatabase

    ! Local variables
    integer :: i, status

    if (associated(l2pcDatabase)) then
      do i = 1, size(l2pcDatabase)
        call DestroyL2PC ( l2pcDatabase(i) )
      end do
      deallocate ( l2pcDatabase, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_deallocate // "l2pcDatabase" )
    end if
  end subroutine DestroyL2PCDatabase

end module L2PC_m

! $Log$
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

