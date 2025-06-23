! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module EmpiricalGeometry                ! For empirically obtaining orbit information

  use HGridsdatabase, only: L1BGeolocation, L1BSubsample
  use MLSKinds, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning

  implicit none
  private

  public :: EmpiricalLongitude, InitEmpiricalGeometry, &
            DestroyEmpiricalGeometry, &
            ForgetOptimumLon0, CFM_InitEmpiricalGeometry, &
            ChooseOptimumLon0, CFM_ResetEmpiricalGeometry

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! ---- Local private declations

  real(r8), pointer, dimension(:), private, save :: &
    & EMPIRICALTERMS => NULL() ! Fourier terms
  real(r8), private, save :: LON0       ! Longitude of first equator crossing
  logical, private, save :: LON0VALID = .false.
  
contains ! ========================= Public Procedures ====================

  ! -----------------------------------------  EmpiricalLongitude  -----
  subroutine EmpiricalLongitude ( geodAngle, lon, tryLon0  ) 
    ! This function returns an empirical longitude for a given
    ! geodetic angle.
    use Constants, only: Deg2Rad

    ! Arguments
    real(r8), dimension(:), intent(in) :: GEODANGLE
    real(r8), dimension(size(geodAngle)), intent(out) :: LON
    real(r8), optional, intent(in) :: TRYLON0

    ! Local variables
    integer :: term
    real(r8) :: useLon0

    ! Executable code
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

  ! ---------------------------------------  InitEmpiricalGeomtry  -----
  subroutine InitEmpiricalGeometry ( root )
    ! This subroutine sets up the empirical geometry from l2cf information

    use Allocate_Deallocate, only: Allocate_test
    use Expr_M, only: EXPR
    use Init_Tables_Module, only: F_ITERATIONS, F_TERMS
    use MoreTree, only: Get_Field_ID
    use Toggles, only: gen, levels, toggle
    use Trace_m, only: trace_begin, trace_end
    use Tree, only: nsons, subtree

    integer, intent(in) :: ROOT         ! Root of tree

    ! Local variables
    real(r8), dimension(2) :: VALUE     ! From EXPR
    integer, dimension(2) :: TheUnits   ! From EXPR
    integer :: I,J                      ! Loop inductors
    integer :: Me = -1                  ! String index for trace
    integer :: SON                      ! Son of root
    integer :: NOTERMS                  ! Number of terms
    integer :: FIELDINDEX               ! From parser

    ! Executable code
    call trace_begin ( me, "InitEmpiricalGeometry", root, &
      & cond=toggle(gen) .and. levels(gen) > 0 )

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
          empiricalTerms(j-1) = value(1)
        end do
      case ( f_iterations )
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & "The Iterations field no longer does anything" )
      end select
    end do
    call trace_end ( "InitEmpiricalGeometry", &
      & cond=toggle(gen) .and. levels(gen) > 0 )

  end subroutine InitEmpiricalGeometry

  ! -----------------------------------  DestroyEmpiricalGeometry  -----

  subroutine DestroyEmpiricalGeometry

    use Allocate_Deallocate, only: Deallocate_test

    call deallocate_test ( empiricalTerms, 'empiricalTerms', ModuleName )

  end subroutine DestroyEmpiricalGeometry

  ! -----------------------------------  CFM_InitEmpiricalGeomtry  -----
  ! This module is for the callable forward model (CFM), in which no tree
  ! is available, to set up EmpiricalGeometry
  subroutine CFM_InitEmpiricalGeometry ( numberIterations, terms )
     use Allocate_Deallocate, only: Allocate_test

     integer, intent(in) :: numberIterations
     real(r8), dimension(:), intent(in) :: terms

     ! Executables
     call Allocate_Test ( empiricalTerms, size(terms), 'empiricalTerms', &
          & ModuleName )
     empiricalTerms = terms

  end subroutine

  ! ----------------------------------  CFM_ResetEmpiricalGeomtry  -----
  ! This module is for the callable forward model (CFM), in which
  ! EmpiricalGeometry is expected to be called multiple times.
  subroutine CFM_ResetEmpiricalGeometry

     ! Excutables
     call DestroyEmpiricalGeometry
     call ForgetOptimumLon0

  end subroutine

  ! ------------------------------------------  ChooseOptimumLon0  -----
  subroutine ChooseOptimumLon0 ( filedatabase, chunk, moduleStr )

    use Allocate_Deallocate, only: allocate_test, deallocate_test
    use Chunks_m, only: MLSChunk_t
    use MLSCommon, only: MLSFile_t

    type (MLSFile_T), dimension(:), pointer :: FILEDATABASE          
    type (MLSChunk_T), intent(in)           :: CHUNK ! This chunk    
    character(len=*), intent(in)            :: moduleStr             

    ! Local variables
    double precision, dimension(:), pointer :: FullArray
    real(r8), dimension(:), pointer         :: TESTPHI ! Angle
    real(r8), dimension(:), pointer         :: TESTLON ! Longitude to match
    real(r8), dimension(:), pointer         :: GUESSLON ! Attempt at match

    ! Executable code

    ! Now we want to establish the value of lon0
    nullify ( testPhi, testLon, guessLon )

    call L1BGeoLocation ( filedatabase, "tpGeodAngle", moduleStr, fullArray )
    call L1BSubsample ( chunk, FullArray, values=testPhi )
    call Deallocate_test ( fullArray, 'fullArray', ModuleName )

    call L1BGeoLocation ( filedatabase, "tpLon", moduleStr, fullArray )
    call L1BSubsample ( chunk, FullArray, values=testLon )
    call Deallocate_test ( fullArray, 'fullArray', ModuleName )

    call Allocate_test ( guessLon, size(testLon), 'guessLon', ModuleName )

    call EmpiricalLongitude ( testPhi, guessLon, tryLon0=0.0_r8 )
    lon0 = - sum ( NormalizeLongitude ( guessLon - testLon ) ) / size(testLon)
    lon0Valid = .true.

    call Deallocate_test ( testPhi, 'testPhi', ModuleName )
    call Deallocate_test ( testLon, 'testLon', ModuleName )
    call Deallocate_test ( guessLon, 'guessLon', ModuleName )

  end subroutine ChooseOptimumLon0

  ! ------------------------------------------  ForgetOptimumLon0  -----
  subroutine ForgetOptimumLon0
    lon0Valid = .false.
  end subroutine ForgetOptimumLon0

  ! ========================================= PRIVATE PROCEUDRES ==========

  ! -----------------------------------------  NormalizeLongitude  -----
  elemental function NormalizeLongitude ( lon )
    real(r8), intent(in) :: LON
    real(r8) :: NormalizeLongitude
    ! Executable code
    NormalizeLongitude = modulo ( lon, 360.0_r8 )
    if (NormalizeLongitude > 180.0) NormalizeLongitude = NormalizeLongitude - 360.0
  end function NormalizeLongitude

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module EmpiricalGeometry

! $Log$
! Revision 2.25  2016/08/09 21:07:18  pwagner
! Survives encounter with non-satellite data
!
! Revision 2.24  2015/06/19 00:37:56  pwagner
! Changed to avoid reading L1BOA
!
! Revision 2.23  2014/09/05 00:51:39  vsnyder
! Add DestroyEmpiricalGeometry, some cannonball polishing
!
! Revision 2.22  2014/09/05 00:49:06  vsnyder
! EmpiricalGeometry.f90 -- wrong comment.
!
! Revision 2.21  2014/08/15 02:56:25  vsnyder
! Remove noIterations, which was only used in code that was removed in 2001
!
! Revision 2.20  2014/03/07 19:21:44  pwagner
! Name_Len changed to nameLen; got from MLSCommon
!
! Revision 2.19  2014/03/04 17:37:48  pwagner
! Must not truncate f.p. terms
!
! Revision 2.18  2014/03/01 03:10:56  vsnyder
! Move units checking to init_tables_module
!
! Revision 2.17  2010/09/17 16:47:08  honghanh
! Add subroutine to deallocate empiricalTerms
!
! Revision 2.15  2010/02/04 23:12:44  vsnyder
! Remove USE or declaration for unreferenced names
!
! Revision 2.14  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.13  2009/05/13 20:41:55  vsnyder
! Get constants from Constants, kinds from MLSKinds
!
! Revision 2.12  2005/06/22 18:57:01  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.11  2005/05/31 17:51:17  pwagner
! Began switch from passing file handles to passing MLSFiles
!
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
