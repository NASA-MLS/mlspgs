! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module L2PC_m
  !=============================================================================

  ! This module contains data types etc. for dealing with the new EMLS L2PC
  ! files.  The first version will deal with ascii files, but later versions
  ! will probably be HDF.

  use Declaration_Table, only: DECLS, ENUM_VALUE, GET_DECL
  use Intrinsic, only: Lit_Indices
  use MLSCommon, only: R8
  use VectorsModule, only: DESTROYVECTORINFO, VECTOR_T, VECTORVALUE_T
  use MatrixModule_1, only: DESTROYMATRIX, MATRIX_T
  use MatrixModule_0, only: M_ABSENT, M_BANDED, M_COLUMN_SPARSE, M_FULL, &
    & MATRIXELEMENT_T
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR, &
    & MLSMSG_ALLOCATE, MLSMSG_DEALLOCATE
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
    type (Vector_T) :: xStar            ! The linearisation x
    type (Vector_T) :: yStar            ! The corresponding y
    type (Matrix_T) :: kStar            ! The jacobian matrix
  end type l2pc_T

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
          call GetSignalName ( qt%signal, line )
          write (unit,*) trim(line)
        end select

        ! Write out the dimensions for the quantity and the edges
        write (unit,*) qt%noSurfs, qt%noInstances, qt%noChans,&
          &  'noSurfs, noInstances, noChans'
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
    integer :: INSTANCE                 ! Loop counter
    integer :: NOQUANTITIES             ! Number of quantities in a vector
    integer :: QUANTITY                 ! Loop counter
    integer :: STATUS                   ! Flag
    integer :: STRINGINDEX              ! Index of string
    integer :: VECTOR                   ! Loop counter

    type (QuantityTemplate_T), pointer :: qt  ! Temporary pointers
    type (Vector_T), pointer :: v       ! Temporary pointer
    type (MatrixElement_T), pointer :: m0 ! A Matrix0 within kStar
    type (Decls) :: decl                ! From the declaration table

    character (len=132) :: line         ! Line of text

    ! executable code

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
        qt => v%quantities(quantity)%template

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
!          call Parse_Signal (line, sigInds, 
        end select

!         ! Write out the dimensions for the quantity and the edges
!         write (unit,*) qt%noSurfs, qt%noInstances, qt%noChans,&
!           &  'noSurfs, noInstances, noChans'
!         write (unit,*) 'surfs'
!         write (unit, rFmt) qt%surfs
!         write (unit,*) 'phi'
!         write (unit, rFmt) qt%phi

!         ! Write the values
!         write (unit,*) 'values'
!         write (unit,rFmt) v%quantities(quantity)%values

      end do                            ! Loop over quantities

      ! Create a dummy template for this vector
      allocate (v%template, STAT=status)
      if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName,&
        & MLSMSG_Allocate//"v%template")
      v%template%noQuantities = noQuantities
!      v%template%totalInstances = sum(v%quantities%template%noInstances)
!      v%template%totalElements = sum(v%quantities%template%noInstances *&
!        & v%quantities%template%instanceLen)
      ! Ignore the quantities stuff for the moment.
      ! If later versions need it we'll put something here
    end do                              ! Loop over xStar/yStar

!     ! Now dump kStar
!     write (unit,*) 'kStar'
!     do blockRow = 1, l2pc%kStar%row%NB
!       do blockCol = 1, l2pc%kStar%col%NB
!         ! Print the type of the matrix
!         m0 => l2pc%kStar%block(blockRow, blockCol)
!         write (unit,*) blockRow, blockCol, m0%kind,&
!           & 'row, col, kind'
!         call get_string ( &
!           & l2pc%kStar%row%vec%quantities(&
!           &    l2pc%kStar%row%quant(blockRow))%template%name, word1 )
!         call get_string ( &
!           & l2pc%kStar%col%vec%quantities(&
!           &    l2pc%kStar%col%quant(blockCol))%template%name, word2 )
!         write (unit,*) trim(word1), l2pc%kStar%row%inst(blockRow), ' , ',&
!           &            trim(word2), l2pc%kStar%col%inst(blockCol)
!         select case (m0%kind)
!         case (M_Absent)
!         case (M_Banded, M_Column_sparse)
!           write (unit,*) 'R1'
!           write (unit,iFmt) m0%R1
!           write (unit,*) 'R2'
!           write (unit,iFmt) m0%R2
!           write (unit,*) 'values'
!           write (unit,rFmt) m0%values
!         case (M_Full)
!           write (unit,rFmt) m0%values
!         end select
!       end do
!     end do

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
  subroutine Read_l2pc_file ( lun, l2pcDatabase )
    ! Read all the bins in an l2pc file
    integer, intent(in) :: lun
    type (l2pc_T), dimension(:), pointer :: l2pcDatabase

    ! Local variables
    type (l2pc_T) :: l2pc
    integer :: dummy
    logical :: eof

    ! Executable code
    eof = .false.
    do while (.not. eof )
      call ReadOneL2PC ( l2pc, lun, eof )
      dummy = AddL2PCToDatabase ( l2pcDatabase, l2pc )
    end do
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
    type (l2pc_T), intent(inout) :: L2PC

    ! Exectuable code
    call DestroyVectorInfo ( l2pc%xStar )
    call DestroyVectorInfo ( l2pc%yStar )
    call DestroyMatrix ( l2pc%kStar )

  end subroutine DestroyL2PC

  ! ------------------------------------------- DestroyL2PCDatabase ---
  subroutine DestroyL2PCDatabase (l2pcDatabase )
    ! Dummy arguments
    type (l2pc_T), dimension(:), pointer :: l2pcDatabase

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

