! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module EmpiricalGeometry                ! For empirically obtaining orbit information

  use Allocate_Deallocate, only: Allocate_test
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error

  implicit none
  private

  public :: EmpiricalLongitude, InitEmpiricalGeometry, ForgetOptimumLon0, &
    & ChooseOptimumLon0

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName = &
    "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

  ! ---- Local private declations

  real(r8), pointer, dimension(:), private, save :: &
    & EMPIRICALTERMS => NULL() ! Fourier terms
  real(r8), private, save :: LON0       ! Longitude of first equator crossing
  logical, private, save :: LON0VALID = .false.
  integer, private, save :: NOITERATIONS ! When finding lon0 
  
contains ! ========================= Public Procedures ====================

  ! ------------------------------------------------ EmpiricalLongitude ---
  subroutine EmpiricalLongitude ( geodAngle, lon, tryLon0  ) 
    ! This function returns an empirical longitude for a given
    ! geodetic angle.
    ! Argument

    use Units, only: Deg2Rad
    real(r8), dimension(:), intent(in) :: GEODANGLE
    real(r8), dimension(size(geodAngle)), intent(out) :: LON
    real(r8), optional, intent(in) :: TRYLON0

    ! Local variables
    integer :: term
    real(r8) :: useLon0

    ! Exectuable code

    if ( .not. associated ( empiricalTerms ) ) call MLSMessage ( &
      & MLSMSG_Error, ModuleName, 'EmpiricalGeomtry information not given in l2cf' )

    if ( present(tryLon0) ) then
      useLon0 = tryLon0
    else
      if ( .not. lon0Valid ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Optimum lon0 has not been established" )
      useLon0 = lon0
    end if

    lon = useLon0
    lon = lon + empiricalTerms(1) * geodAngle

    do term = 2, size(empiricalTerms)
      lon = lon + empiricalTerms(term) * &
        & sin ( deg2Rad * geodAngle * 2*(term-1) )
    end do
    lon = NormalizeLongitude ( lon )
  end subroutine EmpiricalLongitude

  ! ----------------------------------------------- InitEmpiricalGeomtry --
  subroutine InitEmpiricalGeometry ( root )
    ! This subroutine sets up the empirical geometry from l2cf information

    use Expr_M, only: EXPR
    use Init_Tables_Module, only: F_ITERATIONS, F_TERMS
    use Intrinsic, only: PHYQ_DimensionLess
    use MoreTree, only: Get_Field_ID
    use Tree, only: NSONS, SUBTREE

    integer, intent(in) :: ROOT         ! Root of tree

    ! Local parameters
    integer, parameter :: INITNOITERATIONS = 5

    ! Local variables
    real(r8), dimension(2) :: VALUE     ! From EXPR
    integer, dimension(2) :: TheUnits   ! From EXPR
    integer :: I,J                      ! Loop inductors
    integer :: SON                      ! Son of root
    integer :: NOTERMS                  ! Number of terms
    integer :: FIELDINDEX               ! From parser

    ! Executable code

     noIterations = initNoIterations

    ! First get the information from the l2cf
    do i = 2, nsons(root)
      son = subtree( i, root )
      fieldIndex = get_field_id(son)
      select case ( fieldIndex )
      case ( f_terms )
        noTerms = nsons ( son ) - 1
        call Allocate_Test ( empiricalTerms, noTerms, 'empiricalTerms', &
          & ModuleName )
        do j = 2, noTerms + 1
          call expr ( subtree(j,son), theUnits, value )
          if ( theUnits(1) /= PHYQ_Dimensionless ) call MLSMessage ( MLSMSG_Error, &
            & ModuleName, "No units expected for empirical terms" )
          empiricalTerms(j-1) = value(1)
        end do
      case ( f_iterations )
        call expr ( subtree(2,son), theUnits, value )
        if ( theUnits(1) /= PHYQ_Dimensionless ) call MLSMessage ( MLSMSG_Error, &
          & ModuleName, "No units expected for iterations" )
        noIterations = value(1)
      end select
    end do

  end subroutine InitEmpiricalGeometry

  ! -------------------------------------------------- ChooseOptimumLon0 -----
  subroutine ChooseOptimumLon0 ( l1bInfo, chunk )

    use Allocate_Deallocate, only: Deallocate_test
    use Chunks_m, only: MLSChunk_T
    use L1BData, only: L1BData_T, ReadL1BData, DeallocateL1BData, Name_Len, &
      & AssembleL1BQtyName
    use MLSCommon, only: L1BInfo_T
    use MLSFiles, only: mls_hdf_version
    use MLSL2Options, only: LEVEL1_HDFVERSION

    type (L1BInfo_T), intent(in) :: L1BINFO ! Where to find L1 files
    type (MLSChunk_T), intent(in) :: CHUNK ! This chunk

    ! Local parameters
    integer, parameter :: NOLON0OPTIONS = 18

    ! Local variables
    integer :: BESTOPTION(1)            ! With lowest cost
    real(r8), dimension(noLon0Options) :: options
    real(r8), dimension(noLon0Options) :: cost
    type (L1BData_T) :: tpGeodAngle     ! From L1B
    type (L1BData_T) :: tpLon           ! From L1B
    integer :: NOMAFS                   ! From ReadL1B
    integer :: FLAG                     ! From ReadL1B
    integer :: I,J                      ! Loop counters
    real(r8), dimension(:), pointer :: TESTPHI ! Angle
    real(r8), dimension(:), pointer :: TESTLON ! Longitude to match
    real(r8), dimension(:), pointer :: GUESSLON ! Attempt at match
    real(r8) :: LOWLIMIT, HILIMIT       ! For new options
    real(r8) :: DELTA                   ! For new options
    integer ::  hdfVersion
    character(len=Name_Len) :: l1bItemName

    ! Executable code

    hdfVersion = mls_hdf_version(trim(l1bInfo%L1BOAFileName), LEVEL1_HDFVERSION)
    if ( hdfversion <= 0 ) &                                                
      & call MLSMessage ( MLSMSG_Error, ModuleName, &                      
      & 'Illegal hdf version for l1boa file (file missing or non-hdf?)' )    
    ! Now we want to establish the value of lon0
    nullify ( testPhi, testLon, guessLon )
    l1bItemName = AssembleL1BQtyName ( "GHz.tpGeodAngle", hdfVersion, .false. )
    call ReadL1BData ( l1bInfo%l1boaID, trim(l1bItemName), tpGeodAngle, noMAFs, flag, &
      & firstMAF = chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
      & hdfVersion=hdfVersion )
    call Allocate_test ( testPhi, noMAFs, 'testPhi', ModuleName )
    testPhi = tpGeodAngle%DpField(1,1,:)
    call DeallocateL1BData ( tpGeodAngle )

    l1bItemName = AssembleL1BQtyName ( "GHz.tpLon", hdfVersion, .false. )
    call ReadL1BData ( l1bInfo%l1boaID, trim(l1bItemName), tpLon, noMAFs, flag, &
      & firstMAF = chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
      & hdfVersion=hdfVersion )
    call Allocate_test ( testLon, noMAFs, 'testLon', ModuleName )
    testLon = tpLon%DpField(1,1,:)
    call DeallocateL1BData ( tpLon )

    call Allocate_test ( guessLon, size(testLon), 'guessLon', ModuleName )

    call EmpiricalLongitude ( testPhi, guessLon, tryLon0=0.0_r8 )
    lon0 = - sum ( NormalizeLongitude ( guessLon - testLon ) ) / size(testLon)
    lon0Valid = .true.

    call Deallocate_test ( testPhi, 'testPhi', ModuleName )
    call Deallocate_test ( testLon, 'testLon', ModuleName )
    call Deallocate_test ( guessLon, 'guessLon', ModuleName )

  end subroutine ChooseOptimumLon0

  ! ------------------------------------------------ ForgetOptimumLon0 -----
  subroutine ForgetOptimumLon0
    lon0Valid = .false.
  end subroutine ForgetOptimumLon0

  ! ========================================= PRIVATE PROCEUDRES ==========

  ! ----------------------------------------------- NormalizeLongitude ----
  elemental function NormalizeLongitude ( lon )
    real(r8), intent(in) :: LON
    real(r8) :: NormalizeLongitude
    ! Executable code
    NormalizeLongitude = modulo ( lon, 360.0_r8 )
    if (NormalizeLongitude > 180.0) NormalizeLongitude = NormalizeLongitude - 360.0
  end function NormalizeLongitude

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module EmpiricalGeometry

! $Log$
! Revision 2.10  2004/05/19 19:16:09  vsnyder
! Move MLSChunk_t to Chunks_m
!
! Revision 2.9  2003/08/15 23:58:20  vsnyder
! Get PHYQ_... directly from Intrinsic instead of indirectly via Units
!
! Revision 2.8  2002/12/11 22:17:05  pwagner
! Added error checks on hdf version
!
! Revision 2.7  2002/11/13 01:05:28  pwagner
! Actually reads hdf5 radiances
!
! Revision 2.6  2002/10/08 17:36:20  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.5  2002/08/20 22:43:37  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.4  2001/12/16 00:58:06  livesey
! New method for computing lon0 (much more efficient)
!
! Revision 2.3  2001/12/14 01:43:02  livesey
! Various bug fixes
!
! Revision 2.2  2001/12/13 23:05:08  vsnyder
! Add a kind parameter in MODULO; Add a CVS Log comment
!
