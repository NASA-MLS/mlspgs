! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module L2PC_PFA_STRUCTURES
  use MLSCommon, only: I4, R4, R8
  use L2PC_File_Parameters, only: MAX_NO_ELMNTS_PER_SV_COMPONENT, &
                                  MAX_NO_POINTINGS
  use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR, MLSMSG_DEALLOCATE
  implicit none
!  This has the standard structure constructions used for
!    computing the data for the l2pc_xx.dat file.
!  Originally written by W. G. Read on 8/13/1990.
!  Modified, Mar/10/2000, Z. Shippony, to fit for EOS concept
!---------------------------- RCS Ident Info -------------------------------
  private :: Id, ModuleName
  character (LEN=256) :: Id = &
       "$Id$"
  character (LEN=*), parameter :: ModuleName = &
       "$RCSfile$"
!---------------------------------------------------------------------------
!_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
! Parameter declarations
  integer(i4), parameter :: MAXSPS = 20
  integer(i4), parameter :: MAXSUBSPS = 15
  integer(i4), parameter :: MAXGEOM = 13
  integer(i4), parameter :: MAXGEOPHYS = 2
! maxsps + maxgeom + maxgeophys + 4 (for pointing and 3 X field strength
! computations) should equal max_no_sv_components.
! More structures :
! These are for catagories of state vector types. This catagorization is
! because the method of differentiation differs according to catagory.

! This structure contains the a priori limb point pressure
  type LIMB_PRESS
    character(len=8) :: NAME
    integer(i4) :: NO_LIN_VALUES
    logical DER_CALC(6)
    real(r4) :: LIN_VAL(max_no_pointings)
  end type LIMB_PRESS

! This structure contains geometric and any other parameters characterized
! with a single value
  type GEOM_PARAM
    character(len=8) :: NAME
    logical DER_CALC(6)
    real(r4) :: LIN_VAL
  end type GEOM_PARAM

! This structure contains spectroscopic information about the species
  type SPECTRO_PARAM
    character(len=1) :: type
    character(len=8) :: NAME
    integer(i4) :: SPECTAG
    integer(i4) :: NO_PHI_VALUES
    integer(i4) :: NO_ZETA_VALUES
    logical DER_CALC(6)
    real(r4) :: PHI_BASIS(12)
    real(r4) :: ZETA_BASIS(52)
  end type SPECTRO_PARAM

! This structure contains geophysical information about the atmosphere
  type GEOPHYS_PARAM
    character(len=8) :: NAME
    integer(i4) :: NO_LIN_VALUES
    logical DER_CALC(6)
    real(r4) :: LIN_VAL(50)
    real(r4) :: BASIS_PEAKS(52)
  end type GEOPHYS_PARAM

! This structure contains atmospheric composition of the atmosphere
  type ATMOS_COMP
    character(len=8) :: NAME
    integer(i4) :: SPECTAG
    integer(i4) :: NO_LIN_VALUES
    logical FWD_CALC(6)
    logical DER_CALC(6)
    real(r4) :: LIN_VAL(50)
    real(r4) :: BASIS_PEAKS(52)
  end type ATMOS_COMP

! This NEW structure contains info about the K matrix
  type K_MATRIX_INFO
    character(len=8) :: NAME
    integer(i4) :: FIRST_DIM_INDEX
    integer(i4) :: NO_ZETA_BASIS
    integer(i4) :: NO_PHI_BASIS
    real(r4)    :: PHI_BASIS(12)
    real(r4)    :: ZETA_BASIS(52)
  end type K_MATRIX_INFO

  integer(i4), parameter :: MAXPFALINES = 6
  integer(i4), parameter :: MAXPFACH = 45
  integer(i4), parameter :: MAXFILTPTS = 161
  integer(i4), parameter :: MAXLINES = 35

  type :: PFA_SLAB
    integer(i4) :: NO_SPS
    integer(i4) :: NO_LINES
    integer(i4) :: SPECTAG
    character(len=8) :: NAME
    real(r8) :: QLOG(3)
    real(r8) :: V0(maxlines)
    real(r8) :: EL(maxlines)
    real(r8) :: STR(maxlines)
    real(r8) :: W(maxlines)
    real(r8) :: PS(maxlines)
    real(r8) :: N(maxlines)
    real(r8) :: N1(maxlines)
    real(r8) :: N2(maxlines)
    real(r8) :: GAMMA(maxlines)
    real(r8) :: DELTA(maxlines)
  end type PFA_SLAB

!------------------------------------------------------------
! This structure contains the "slabs preps arrays"
  type SLABS_STRUCT
    integer(i4) :: no_lines
    real(r8), dimension(:), pointer :: v0s => NULL()
    real(r8), dimension(:), pointer :: x1 => NULL()
    real(r8), dimension(:), pointer :: y => NULL()
    real(r8), dimension(:), pointer :: yi => NULL()
    real(r8), dimension(:), pointer :: slabs1 => NULL()
    real(r8), dimension(:), pointer :: dx1_dv0 => NULL()
    real(r8), dimension(:), pointer :: dy_dv0 => NULL()
    real(r8), dimension(:), pointer :: dslabs1_dv0 => NULL()
  end type SLABS_STRUCT

contains
  
  ! -------------------------------------------- AllocateOneSlabs ---------
  subroutine AllocateOneSlabs ( slabs, nl )
    ! Allocates the commonly used items in a slabs, or all if the optional
    ! Full parameter is set
    type (slabs_struct), intent(inout) :: slabs ! Slabs to allocate
    integer, intent(in) :: nl         ! Number of lines
    
    ! Local variables
    integer :: myl

    ! Executable code
    myl = MAX(1,nl)
    call Allocate_test ( slabs%v0s, myl,         'v0s',         ModuleName )
    call Allocate_test ( slabs%x1, myl,          'x1',          ModuleName )
    call Allocate_test ( slabs%y, myl,           'y',           ModuleName )
    call Allocate_test ( slabs%yi, myl,          'yi',          ModuleName )
    call Allocate_test ( slabs%slabs1, myl,      'slabs1',      ModuleName )
    call Allocate_test ( slabs%dx1_dv0, myl,     'dx1_dv0',     ModuleName )
    call Allocate_test ( slabs%dy_dv0, myl,      'dy_dv0',      ModuleName )
    call Allocate_test ( slabs%dslabs1_dv0, myl, 'dslabs1_dv0', ModuleName )
  end subroutine AllocateOneSlabs
  
  ! ------------------------------------------ DeallocateOneSlabs ---------
  subroutine DeallocateOneSlabs ( slabs, inName )
    ! Allocates the commonly used items in a slabs, or all if the optional
    ! Full parameter is set
    type (slabs_struct), intent(inout) :: slabs ! Slabs to allocate
    character (len=*), intent(in) :: inName ! ModuleName of caller
    
    ! Executable code
    call Deallocate_test ( slabs%v0s,         'v0s',         ModuleName )
    call Deallocate_test ( slabs%x1,          'x1',          ModuleName )
    call Deallocate_test ( slabs%y,           'y',           ModuleName )
    call Deallocate_test ( slabs%yi,          'yi',          ModuleName )
    call Deallocate_test ( slabs%slabs1,      'slabs1',      ModuleName )
    call Deallocate_test ( slabs%dx1_dv0,     'dx1_dv0',     ModuleName )
    call Deallocate_test ( slabs%dy_dv0,      'dy_dv0',      ModuleName )
    call Deallocate_test ( slabs%dslabs1_dv0, 'dslabs1_dv0', ModuleName )
  end subroutine DeallocateOneSlabs
 
  ! ------------------------------------------- DestroyCompleteSlabs -----
  subroutine DestroyCompleteSlabs ( slabs )
    ! Destroys all the components of a slabs
    type (slabs_struct), dimension(:,:), pointer :: slabs

    integer :: I
    integer :: J
    ! Executable code
    do i = 1, size(slabs,2)
      do j = 1, size(slabs,1)
        call DeallocateOneSlabs ( slabs(j,i), ModuleName )
      end do
    end do
    deallocate ( slabs, stat=i )
    if ( i /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//'slabs' )
  end subroutine DestroyCompleteSlabs
  
end module L2PC_PFA_STRUCTURES
! $Log$
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model
!
! Revision 1.11.2.6  2001/09/13 01:50:26  livesey
! Added DestroyCompleteSlabs
!
! Revision 1.11.2.5  2001/09/10 23:48:22  livesey
! Added some use statements
!
! Revision 1.11.2.4  2001/09/10 20:46:56  livesey
! Trimmed stuff out
!
! Revision 1.11.2.3  2001/09/10 20:03:57  livesey
! Tidied up a bit
!
! Revision 1.11.2.2  2001/09/10 20:03:29  livesey
! Added DeallocateOneSlabs
!
! Revision 1.11.2.1  2001/09/10 19:56:52  livesey
! Added AllocateOneSlabs
!
! Revision 1.11  2001/06/21 13:07:08  zvi
! Speed enhancement MAJOR update
!
! Revision 1.10  2001/06/07 23:39:31  pwagner
! Added Copyright statement
!
! Revision 1.9  2001/04/03 07:32:45  zvi
! Modify the spectral structure - eliminating sps_ from the names
!
! Revision 1.8  2001/03/31 23:40:55  zvi
! Eliminate l2pcdim (dimension parameters) move to allocatable ..
!
! Revision 1.7  2001/02/19 22:20:40  zvi
! Latest modification: Conv/NoConv
!
! Revision 1.6  2001/02/19 22:14:21  zvi
!
! Revision 1.1  2000/06/21 21:56:15  zvi
! First version D.P.
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
