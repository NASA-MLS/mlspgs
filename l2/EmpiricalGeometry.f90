! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module EmpiricalGeometry                ! For empirically obtaining orbit information

  use Allocate_Deallocate, only: Allocate_test, Deallocate_test
  use Expr_M, only: EXPR
  use MLSCommon, only: L1BInfo_T, R8, MLSChunk_T
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use MoreTree, only: Get_Field_ID
  use Tree, only: NSONS, SUBTREE, NODE_ID
  use Units, only: Deg2Rad, PHYQ_Angle, PHYQ_DimensionLess
  use L1BData, only: L1BData_T, ReadL1BData, DeallocateL1BData
  use Init_Tables_Module, only: F_ITERATIONS, F_TERMS

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
    integer, intent(in) :: ROOT         ! Root of tree

    ! Local parameters
    integer, parameter :: INITNOITERATIONS = 5

    ! Local variables
    real(r8), dimension(2) :: VALUE     ! From EXPR
    integer, dimension(2) :: UNITS      ! From EXPR
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
          call expr ( subtree(j,son), units, value )
          if ( units(1) /= PHYQ_Dimensionless ) call MLSMessage ( MLSMSG_Error, &
            & ModuleName, "No units expected for empirical terms" )
          empiricalTerms(j-1) = value(1)
        end do
      case ( f_iterations )
        call expr ( subtree(2,son), units, value )
        if ( units(1) /= PHYQ_Dimensionless ) call MLSMessage ( MLSMSG_Error, &
          & ModuleName, "No units expected for iterations" )
        noIterations = value(1)
      end select
    end do

  end subroutine InitEmpiricalGeometry

  ! -------------------------------------------------- ChooseOptimumLon0 -----
  subroutine ChooseOptimumLon0 ( l1bInfo, chunk )
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

    ! Executable code

    ! Now we want to establish the value of lon0
    nullify ( testPhi, testLon, guessLon )
    call ReadL1BData ( l1bInfo%l1boaID, "GHz.tpGeodAngle", tpGeodAngle, noMAFs, flag, &
      & firstMAF = chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex )
    call Allocate_test ( testPhi, noMAFs, 'testPhi', ModuleName )
    testPhi = tpGeodAngle%DpField(1,1,:)
    call DeallocateL1BData ( tpGeodAngle )

    call ReadL1BData ( l1bInfo%l1boaID, "GHz.tpLon", tpLon, noMAFs, flag, &
      & firstMAF = chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex )
    call Allocate_test ( testLon, noMAFs, 'testLon', ModuleName )
    testLon = tpLon%DpField(1,1,:)
    call DeallocateL1BData ( tpLon )

    call Allocate_test ( guessLon, size(testLon), 'guessLon', ModuleName )

    ! First get a dataset to compare against
    do i = 1, noLon0Options
      options(i) = -180.0 + (i-1) * 360.0 / noLon0Options
    end do

    do i = 1, noIterations
      do j = 1, noLon0Options
        call EmpiricalLongitude ( testPhi, guessLon, tryLon0=options(j) )
        cost(j) = sum ( abs( NormalizeLongitude( testLon - guessLon ) ) )
      end do
      bestOption = minloc ( cost, 1 )
      if ( i < noIterations ) then
        lowLimit = options ( max ( bestOption(1)-1, 1 ) )
        hiLimit = options ( min ( bestOption(1)+1, noLon0Options ) )
        delta = (hiLimit-lowLimit)/(noLon0Options-1)
        do j = 1, noLon0Options
          options(j) = lowLimit + (j-1)*delta
        end do
      end if
    end do
    lon0 = options(bestOption(1))
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

end module EmpiricalGeometry

! $Log$
! Revision 2.3  2001/12/14 01:43:02  livesey
! Various bug fixes
!
! Revision 2.2  2001/12/13 23:05:08  vsnyder
! Add a kind parameter in MODULO; Add a CVS Log comment
!
