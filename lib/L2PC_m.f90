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
  use Init_Tables_Module, only: L_RADIANCE, L_TEMPERATURE, L_VMR
  use String_Table, only: GET_STRING
  use MLSSignals_m, only: GETSIGNALNAME

  implicit NONE
  private
  
  public :: AddL2PCBinToDatabase, DestroyL2PCBin

  ! Public types
  type, public :: l2pcBin_T
    type (Vector_T) :: xStar            ! The linearisation x
    type (Vector_T) :: yStar            ! The corresponding y
    type (Matrix_T) :: kStar            ! The jacobian matrix
  end type l2pcBin_T

  !---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
    & "$Id$"
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains ! ============= Public Procedures ==========================

  ! --------------------------------------- WriteL2PCBin ---------------
  subroutine WriteL2PCBin ( l2pcBin, unit )
    ! This subroutine writes an l2pc bin to a file
    ! Currently this file is ascii, later it will be
    ! some kind of HDF file

    ! Dummy arguments
    type (l2pcBin_T), intent(in), target :: l2pcBin
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
        v => l2pcBin%xStar
      else
        v => l2pcBin%yStar
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
    do blockRow = 1, l2pcBin%kStar%row%NB
      do blockCol = 1, l2pcBin%kStar%row%NB
        ! Print the type of the matrix
        m0 => l2pcBin%kStar%block(blockRow, blockCol)
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

  end subroutine WriteL2PCBin
  
  ! --------------------------------------- WriteL2PCBin ---------------
  subroutine ReadL2PCBin ( l2pcBin, unit )
    ! This subroutine writes an l2pc bin to a file
    ! Currently this file is ascii, later it will be
    ! some kind of HDF file

    ! Dummy arguments
    type (l2pcBin_T), intent(out), target :: l2pcBin
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
!         v => l2pcBin%xStar
!       else
!         v => l2pcBin%yStar
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
!     do blockRow = 1, l2pcBin%kStar%row%NB
!       do blockCol = 1, l2pcBin%kStar%row%NB
!         ! Print the type of the matrix
!         m0 => l2pcBin%kStar%block(blockRow, blockCol)
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

  end subroutine ReadL2PCBin
  
  ! ------------------------------------  Add l2pc bin to database ----
  integer function AddL2PCBinToDatabase ( Database, Item )
    
    ! This function simply adds an l2pc bin to a database of said l2pc bins.
    
    type(l2pcBin_T), dimension(:), pointer :: Database
    type(l2pcBin_T) :: Item
    
    type(l2pcBin_T), dimension(:), pointer :: TempDatabase
    
    include "addItemToDatabase.f9h"

    AddL2PCBinToDatabase = newSize
  end function AddL2PCBinToDatabase

  ! ----------------------------------------------- DestroyL2PCBin ----
  subroutine DestroyL2PCBin ( l2pcBin )
    ! Dummy arguments
    type (l2pcBin_T), intent(inout) :: L2PCBIN

    ! Exectuable code
    call DestroyVectorInfo ( l2pcBin%xStar )
    call DestroyVectorInfo ( l2pcBin%yStar )
    call DestroyMatrix ( l2pcBin%kStar )

  end subroutine DestroyL2PCBin

end module L2PC_m

! $Log$
! Revision 2.1  2001/04/24 20:07:44  livesey
! Moved in from l2
!

