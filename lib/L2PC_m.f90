! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module L2PC_m
  !=============================================================================

  ! This module contains data types etc. for dealing with the new EMLS L2PC
  ! files.  The first version will deal with ascii files, but later versions
  ! will probably be HDF.

  use MLSCommon, only: R8
  use VectorsModule, only: DESTROYVECTORINFO, VECTOR_T, VECTORVALUE_T
  use MatrixModule_1, only: DESTROYMATRIX, MATRIX_T
  use MatrixModule_0, only: M_ABSENT, M_BANDED, M_COLUMN_SPARSE, M_FULL, &
    & MATRIXELEMENT_T
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR, &
    & MLSMSG_ALLOCATE, MLSMSG_DEALLOCATE
  use QuantityTemplates, only: QUANTITYTEMPLATE_T
  use Intrinsic, only: L_RADIANCE, L_TEMPERATURE, L_VMR
  use String_Table, only: GET_STRING
  use MLSSignals_m, only: GETSIGNALNAME

  implicit NONE
  private
  
  public :: AddL2PCToDatabase, DestroyL2PC, DestroyL2PCDatabase

  ! Public types
  type, public :: l2pc_T
    type (Vector_T) :: xStar            ! The linearisation x
    type (Vector_T) :: yStar            ! The corresponding y
    type (Matrix_T) :: kStar            ! The jacobian matrix
  end type l2pc_T

  !---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
    & "$Id$"
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
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

    ! Local variables
    integer :: blockRow                 ! Index
    integer :: blockCol                 ! Index
    integer :: quantity                 ! Loop counter
    integer :: instance                 ! Loop counter
    integer :: vector                   ! Loop counter

    type (QuantityTemplate_T), pointer :: qt  ! Temporary pointers
    type (Vector_T), pointer :: v       ! Temporary pointer
    type (MatrixElement_T), pointer :: m0 ! A Matrix0 within kStar

    character (len=132) :: line         ! Line of text

    ! executable code

    ! First dump the xStar and yStar
    do vector = 1, 2
      ! Identify vector
      if ( vector == 1 ) then
        v => l2pc%xStar
      else
        v => l2pc%yStar
      end if

      ! Loop over quantities
      do quantity = 1, size(v%quantities)
        qt => v%quantities(quantity)%template

        ! Write quantity type
        call get_string ( qt%quantityType, line )
        write (unit,*) trim(line)
        
        ! Write other info associated with type
        select case ( qt%quantityType )
        case (l_vmr)
          call get_string ( qt%molecule, line )
          write (unit,*) trim(line)
        case (l_radiance)
          call GetSignalName ( qt%signal, line )
          write (unit,*) trim(line)
        end select

        ! Write out the dimensions for the quantity and the edges
        write (unit,*) qt%noSurfs, qt%noInstances, qt%noChans
        write (unit,*) qt%surfs
        write (unit,*) qt%phi

        ! Write the values
        write (unit,*) v%quantities(quantity)%values

      end do                            ! Loop over quantities
    end do                              ! Loop over xStar/yStar

    ! Now dump kStar
    do blockRow = 1, l2pc%kStar%row%NB
      do blockCol = 1, l2pc%kStar%row%NB
        ! Print the type of the matrix
        m0 => l2pc%kStar%block(blockRow, blockCol)
        write (unit,*) m0%kind
        select case (m0%kind)
        case (M_Absent)
        case (M_Banded, M_Column_sparse)
          write (unit,*) m0%R1
          write (unit,*) m0%R2
          write (unit,*) m0%values
        case (M_Full)
          write (unit,*) m0%values
        end select
      end do
    end do

  end subroutine WriteOneL2PC
  
  ! --------------------------------------- WriteL2PC ---------------
  subroutine ReadOneL2PC ( l2pc, unit )
    ! This subroutine writes an l2pc to a file
    ! Currently this file is ascii, later it will be
    ! some kind of HDF file

    ! Dummy arguments
    type (l2pc_T), intent(out), target :: l2pc
    integer, intent(in) :: unit

    ! Local variables
    integer :: blockRow                 ! Index
    integer :: blockCol                 ! Index
    integer :: quantity                 ! Loop counter
    integer :: instance                 ! Loop counter
    integer :: vector                   ! Loop counter

    type (QuantityTemplate_T), pointer :: qt  ! Temporary pointers
    type (Vector_T), pointer :: v       ! Temporary pointer
    type (MatrixElement_T), pointer :: m0 ! A Matrix0 within kStar

    character (len=132) :: line         ! Line of text

    ! executable code

    ! First read the xStar and yStar
!     do vector = 1, 2
!       ! Identify vector
!       if ( vector == 1 ) then
!         v => l2pc%xStar
!       else
!         v => l2pc%yStar
!       end if

!       ! Loop over quantities
!       do quantity = 1, size(v%quantities)
!         qt => v%quantities(quantity)%template

!         ! Write quantity type
!         call get_string ( qt%quantityType, line )
!         write (unit,*) trim(line)
        
!         ! Write other info associated with type
!         select case ( qt%quantityType )
!         case (l_vmr)
!           call get_string ( qt%molecule, line )
!           write (unit,*) trim(line)
!         case (l_radiance)
!           call GetSignalName ( qt%signal, line )
!           write (unit,*) trim(line)
!         end select

!         ! Write out the dimensions for the quantity and the edges
!         write (unit,*) qt%noSurfs, qt%noInstances, qt%noChans
!         write (unit,*) qt%surfs
!         write (unit,*) qt%phi

!         ! Write the values
!         write (unit,*) v%quantities(quantity)%values

!       end do                            ! Loop over quantities
!     end do                              ! Loop over xStar/yStar

!     ! Now dump kStar
!     do blockRow = 1, l2pc%kStar%row%NB
!       do blockCol = 1, l2pc%kStar%row%NB
!         ! Print the type of the matrix
!         m0 => l2pc%kStar%block(blockRow, blockCol)
!         write (unit,*) m0%kind
!         select case (m0%kind)
!         case (M_Absent)
!         case (M_Banded, M_Column_sparse)
!           write (unit,*) m0%R1
!           write (unit,*) m0%R2
!           write (unit,*) m0%values
!         case (M_Full)
!           write (unit,*) m0%values
!         end select
!       end do
!     end do

  end subroutine ReadOneL2PC
  
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
! Revision 2.2  2001/04/24 20:20:48  livesey
! Word bin dropped from various places e.g. type
!
! Revision 2.1  2001/04/24 20:07:44  livesey
! Moved in from l2
!

