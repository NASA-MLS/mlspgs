! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module L2PC_PFA_STRUCTURES
  use MLSCommon, only: I4, R4, R8
  use L2PC_File_Parameters, only: MAX_NO_POINTINGS
  implicit none
  public

  interface DUMP
    module procedure Dump_Slabs_Struct, Dump_Slabs_Struct_2D
  end interface

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
  private :: not_used_here 
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
    ! For temperature derivatives.  Most are logarithmic derivatives,
    ! so dz_dT really means 1/z dz_dT.
    real(r8), dimension(:), pointer :: dv0s_dT => NULL()    ! not * 1 / v0s
    real(r8), dimension(:), pointer :: dx1_dT => NULL()     ! / x1
    real(r8), dimension(:), pointer :: dy_dT => NULL()      ! / y
    real(r8), dimension(:), pointer :: dyi_dT => NULL()     ! / yi
    real(r8), dimension(:), pointer :: dslabs1_dT => NULL() ! / slabs1

  end type SLABS_STRUCT

contains

  ! -------------------------------------------- AllocateOneSlabs ---------
  subroutine AllocateOneSlabs ( slabs, nl, TempDer )
    ! Allocates the items in a slabs structure
    use Allocate_Deallocate, only: ALLOCATE_TEST
    type (slabs_struct), intent(inout) :: slabs ! Slabs to allocate
    integer, intent(in) :: nl                   ! Number of lines
    logical, intent(in), optional :: TempDer    ! "Allocate temperature
                                                !  derivative fields"

    ! Local variables
    integer :: myl
    logical :: MyDer

    ! Executable code
    myl = MAX(1,nl)
    myDer = .false.
    if ( present(tempDer) ) myDer = tempDer

    call Allocate_test ( slabs%v0s,         myl, 'v0s',         ModuleName )
    call Allocate_test ( slabs%x1,          myl, 'x1',          ModuleName )
    call Allocate_test ( slabs%y,           myl, 'y',           ModuleName )
    call Allocate_test ( slabs%yi,          myl, 'yi',          ModuleName )
    call Allocate_test ( slabs%slabs1,      myl, 'slabs1',      ModuleName )
    call Allocate_test ( slabs%dx1_dv0,     myl, 'dx1_dv0',     ModuleName )
    call Allocate_test ( slabs%dy_dv0,      myl, 'dy_dv0',      ModuleName )
    call Allocate_test ( slabs%dslabs1_dv0, myl, 'dslabs1_dv0', ModuleName )
    if ( myDer ) then
      call Allocate_test ( slabs%dv0s_dT,    myl, 'dv0s_dT',    ModuleName )
      call Allocate_test ( slabs%dx1_dT,     myl, 'dx1_dT',     ModuleName )
      call Allocate_test ( slabs%dy_dT,      myl, 'dy_dT',      ModuleName )
      call Allocate_test ( slabs%dyi_dT,     myl, 'dyi_dT',     ModuleName )
      call Allocate_test ( slabs%dslabs1_dT, myl, 'dslabs1_dT', ModuleName )
    end if
    slabs%no_lines = nl
    if ( nl == 0 ) then
      slabs%v0s = 0.0_r8
      slabs%x1 = 0.0_r8
      slabs%y = 0.0_r8
      slabs%yi = 0.0_r8
      slabs%slabs1 = 0.0_r8
      slabs%dx1_dv0 = 0.0_r8
      slabs%dy_dv0 = 0.0_r8
      slabs%dslabs1_dv0 = 0.0_r8
      if ( myDer ) then
        slabs%dv0s_dT = 0.0_r8
        slabs%dx1_dT = 0.0_r8
        slabs%dy_dT = 0.0_r8
        slabs%dyi_dT = 0.0_r8
        slabs%dslabs1_dT = 0.0_r8
      end if
    end if
  end subroutine AllocateOneSlabs

  ! --------------------------------------------  AllocateSlabs  ----------
  subroutine AllocateSlabs ( Slabs, No_Ele, Catalog, Caller, TempDer )
  ! Allocate an array of slabs structures, and then the items in each one

    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Error
    use SpectroscopyCatalog_m, only: Catalog_T

    type (slabs_struct), dimension(:,:), pointer :: Slabs
    integer, intent(in) :: No_Ele
    type (catalog_t), dimension(:), intent(in) :: Catalog
    character(len=*), intent(in) :: Caller
    logical, intent(in), optional :: TempDer    ! "Allocate temperature
                                                !  derivative fields"

    integer :: I, J

    allocate ( slabs(no_ele, size(catalog)), stat=i )
    if ( i /= 0 ) call MLSMessage ( MLSMSG_Error, Caller, &
      & MLSMSG_Allocate//"slabs" )

    do i = 1, size(catalog)
      do j = 1, no_ele
        call AllocateOneSlabs ( slabs(j,i), size(catalog(i)%lines), TempDer )
      end do
    end do

  end subroutine AllocateSlabs

  ! ------------------------------------------ DeallocateAllSlabs ---------
  subroutine DeallocateAllSlabs ( Slabs, inName )
    ! Allocates the items in a slabs
    type (slabs_struct), intent(inout), dimension(:,:) :: Slabs
    character(len=*), intent(in) :: InName

    integer :: I
    integer :: J
    ! Executable code
    do i = 1, size(slabs,2)
      do j = 1, size(slabs,1)
        call DeallocateOneSlabs ( slabs(j,i), inName )
      end do
    end do
  end subroutine DeallocateAllSlabs

  ! ------------------------------------------ DeallocateOneSlabs ---------
  subroutine DeallocateOneSlabs ( slabs, inName )
    ! Allocates the items in a slabs
    use Allocate_Deallocate, only: DEALLOCATE_TEST
    type (slabs_struct), intent(inout) :: slabs ! Slabs to deallocate
    character (len=*), intent(in) :: inName ! ModuleName of caller

    ! Executable code
    call Deallocate_test ( slabs%v0s,         'v0s',         inName )
    call Deallocate_test ( slabs%x1,          'x1',          inName )
    call Deallocate_test ( slabs%y,           'y',           inName )
    call Deallocate_test ( slabs%yi,          'yi',          inName )
    call Deallocate_test ( slabs%slabs1,      'slabs1',      inName )
    call Deallocate_test ( slabs%dx1_dv0,     'dx1_dv0',     inName )
    call Deallocate_test ( slabs%dy_dv0,      'dy_dv0',      inName )
    call Deallocate_test ( slabs%dslabs1_dv0, 'dslabs1_dv0', inName )
    call Deallocate_test ( slabs%dv0s_dT,     'dv0s_dT',     inName )
    call Deallocate_test ( slabs%dx1_dT,      'dx1_dT',      inName )
    call Deallocate_test ( slabs%dy_dT,       'dy_dT',       inName )
    call Deallocate_test ( slabs%dyi_dT,      'dyi_dT',      inName )
    call Deallocate_test ( slabs%dslabs1_dT,  'dslabs1_dT',  inName )
  end subroutine DeallocateOneSlabs

  ! ------------------------------------------- DestroyCompleteSlabs -----
  subroutine DestroyCompleteSlabs ( Slabs )
    ! Destroys all the components of a slabs
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Deallocate
    type (slabs_struct), dimension(:,:), pointer :: Slabs

    integer :: I
    ! Executable code
    call deallocateAllSlabs ( slabs, moduleName )
    deallocate ( slabs, stat=i )
    if ( i /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//'slabs' )
  end subroutine DestroyCompleteSlabs

  ! ------------------------------------------  Dump_Slabs_Struct  -----
  subroutine Dump_Slabs_Struct ( The_Slabs_Struct, Name )

    use Dump_0, only: Dump
    use Output_m, only: Output

    type(slabs_struct), intent(in) :: The_Slabs_Struct
    character(len=*), intent(in), optional :: Name

    integer :: NL

    call output ( 'Slabs_Struct ' )
    if ( present(name) ) call output ( trim(name) )
    nl = the_slabs_struct%no_lines
    if ( nl == 0 ) then
      call output ( ' is empty', advance='yes' )
    else
      call output ( '', advance='yes' )
      call dump ( the_slabs_struct%v0s(:nl), name='v0s' )
      call dump ( the_slabs_struct%x1(:nl), name='x1' )
      call dump ( the_slabs_struct%y(:nl), name='y' )
      call dump ( the_slabs_struct%yi(:nl), name='yi' )
      call dump ( the_slabs_struct%slabs1(:nl), name='slabs1' )
      call dump ( the_slabs_struct%dx1_dv0(:nl), name='dx1_dv0' )
      call dump ( the_slabs_struct%dy_dv0(:nl), name='dy_dv0' )
      call dump ( the_slabs_struct%dslabs1_dv0(:nl), name='dslabs1_dv0' )
      if ( associated (the_slabs_struct%dslabs1_dT) ) then
        call dump ( the_slabs_struct%dv0s_dT(:nl), name='dv0s_dT' )
        call dump ( the_slabs_struct%dx1_dT(:nl), name='dx1_dT' )
        call dump ( the_slabs_struct%dy_dT(:nl), name='dy_dT' )
        call dump ( the_slabs_struct%dyi_dT(:nl), name='dyi_dT' )
        call dump ( the_slabs_struct%dslabs1_dT(:nl), name='dslabs1_dT' )
      end if
    end if

  end subroutine Dump_Slabs_Struct


  ! ---------------------------------------  Dump_Slabs_Struct_2D  -----
  subroutine Dump_Slabs_Struct_2D ( The_Slabs_Struct, Name )

    use Output_m, only: Output

    type(slabs_struct), intent(in) :: The_Slabs_Struct(:,:)
    character(len=*), intent(in), optional :: Name

    integer :: I, J

    call output ( 'Slabs Struct' )
    if ( present(name) ) call output ( ' ' // trim(name) )
    call output ( ', SIZE = ' )
    call output ( size(the_slabs_struct,1) )
    call output ( ' X ' )
    call output ( size(the_slabs_struct,2), advance='yes' )
    do j = 1, size(the_slabs_struct,2)
      do i = 1, size(the_slabs_struct,1)
        call output ( 'Item ' )
        call output ( i )
        call output ( ', ' )
        call output ( j, advance='yes' )
        call dump ( the_slabs_struct(i,j) )
      end do
    end do

  end subroutine Dump_Slabs_Struct_2D

  ! ----------------------------------------------  not_used_here  -----
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module L2PC_PFA_STRUCTURES
! $Log$
! Revision 2.10  2004/03/20 04:08:55  vsnyder
! Steps along the way to analytic temperature derivatives
!
! Revision 2.9  2003/07/09 23:39:56  vsnyder
! Add AllocateSlabs
!
! Revision 2.8  2003/07/04 02:46:33  vsnyder
! Create DeallocateAllSlabs subroutine, futzing
!
! Revision 2.7  2003/05/17 01:20:52  vsnyder
! Futzing
!
! Revision 2.6  2003/05/16 23:52:08  livesey
! Removed reference to spectags.  What does this module do now anyway?
!
! Revision 2.5  2003/05/05 23:00:25  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.4.2.2  2003/03/13 00:26:50  vsnyder
! Don't dump an empty SLABS structure
!
! Revision 2.4.2.1  2003/03/12 21:49:39  vsnyder
! Add dumpers for scalar and 2-d Slabs_Struct
!
! Revision 2.4  2002/10/08 17:08:05  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.3  2002/10/03 22:16:50  vsnyder
! Move a USE from module scope to procedure scope
!
! Revision 2.2  2002/09/26 23:58:35  livesey
! Clear arrays in zero lines case (think about whether we need to do the
! max(nl,1) stuff later
!
! Revision 2.1  2002/05/23 22:03:35  zvi
! Prevention of zero allocation size
!
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
