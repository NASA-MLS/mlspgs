! Copyright 2015, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module Path_Representation_m
!=============================================================================

  ! Types to represent information about the path of integration of the
  ! radiative-transfer equation, and information along it.

  ! Procedures to calculate the path.

  use Geolocation_0, only: ECR_t, RG

  implicit NONE
  private

  type, public :: Path_t
    ! A path is represented by two lines, one from the instrument to the
    ! tangent or intersection, and one from there onward. The given line is the
    ! one from the instrument to the tangent or intersection, represented by
    ! Lines(1,1) + s * Lines(2,1). A line defined by
    ! Lines(1,2) + s * Lines(2,2) is produced, which is the continuation
    ! of Lines(:,1) after the tangent point if Lines(:,1) does not
    ! intersect the Earth reference ellipsoid, or the reflection of
    ! Lines(:,1) if it does intersect the Earth reference ellipsoid. 
    ! Lines(2,1) and Lines(2,2) are made to be unit vectors.
    type(ECR_t) :: Lines(2,2)
    real(rg) :: SMax(2), SMin(2) ! Interesting intervals of S along Lines
  contains
    procedure :: Get_Path_Ready
  end type Path_t

  type, public :: Facets_and_Vertices_t
    integer :: NF=0                     ! Number of used elements in Facets
    integer :: NV=0                     ! Number of used elements in Vertices
    integer, allocatable :: Facets(:)   ! Indices in QTM_Tree%Q
    integer, allocatable :: Vertices(:) ! Indices in QTM_Tree%Geo_In
  end type Facets_and_Vertices_t

  public :: Geometric_Path, Path_Continuation, Tangent_Geodetic_Coordinates
  public :: Union_Paths

  interface Union_Paths
    procedure Add_One_Path_To_Union, Add_Paths_To_Union
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine Geometric_Path ( FwdModelIn, FwdModelExtra, MAF, Ref_MIF, &
                            & Tngt, Tngt_Geod, H, Path )

    ! Use the instrument position, tangent position, and line-of-sight, for
    ! a specified MAF and reference MIF, to calculate segments of the path
    ! of integration.  See wvs-143.

    use Geolocation_0, only: ECR_t, H_V_Geod, Norm2
    use Intrinsic, only: L_ECRtoFOV, L_ScECR
    use MLSKinds, only: RP
    use VectorsModule, only: GetVectorQuantityByType, Vector_T, VectorValue_T

    type(vector_t), intent(in) :: FwdModelIn, FwdModelExtra
    integer, intent(in) :: MAF
    integer, intent(in) :: Ref_MIF
    type(ECR_t), intent(in) :: Tngt           ! Tangent position at MAF and
                                              ! Ref_MIF
    type(h_v_geod), intent(in) :: Tngt_Geod   ! Geodetic surface coordinates of
                                              ! the tangent at MAF and Ref_MIF,
                                              ! from Tangent_Geodetic_Coordinates
    real(rp), intent(in) :: H                 ! Geodetic height at the tangent
                                              ! point on the integration path
    type(path_t), intent(out) :: Path

    real(rp) :: Dist                          ! | tngt - inst |
    type(vectorValue_t), pointer :: ECRtoFOV
    type(ECR_t) :: Fwm_Tngt                   ! Forward model's tangent position
    type(ECR_t) :: Inst                       ! Instrument position
    type(ECR_t) :: LOS                        ! Line-of-sight from third column
                                              ! of ECRtoFOV
    type(ECR_t) :: Normal                     ! Normal to viewing plane
    type(vectorValue_t), pointer :: ScECR     ! Instrument position
    type(ECR_t) :: Tngt_surf                  ! Vector to Tngt_Geod
    type(ECR_t) :: V                          ! Unit vector in tangent direction

    ECRtoFOV => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_ECRtoFOV )
    scECR => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_scECR )

    inst%xyz = scECR%value3(1:3,ref_MIF,MAF)
    dist = norm2 ( tngt - inst )

    los%xyz = ECRtoFOV%value3(7:9,ref_MIF,MAF) ! Third column is LOS

    normal = los .cross. tngt            ! Normal to viewing plane
    tngt_surf = tngt_geod%ECR()          ! Surface ECR of tangent
    v = tngt - tngt_surf                 ! V_r = T_r - T_g in wvs-143
    v = v / v%norm2()                    ! Unit vector \hat{V} in wvs-143
    fwm_tngt = tngt_surf + h * v         ! \mathbf{T}_i in wvs-143.
    los = normal%cross_norm ( fwm_tngt ) ! U_i in wvs-143
    path%lines(2,1) = normal .cross. fwm_tngt ! LOS for FWM
    path%lines(1,1) = fwm_tngt - dist * path%lines(2,1) ! Put the reference near
                                         ! the instrument.  Using this instead
                                         ! of Fwm_Tngt puts it outside the
                                         ! reference ellipsoid, an assumption
                                         ! by Path%Get_Path_Ready
    call path%get_path_ready ( inst )    ! Compute path%lines(:,2), either
                                         ! a continuation after the tangent,
                                         ! or a reflection from the Earth's
                                         ! surface.  Fill path%sMin and
                                         ! path%sMax.

  end subroutine Geometric_Path

  subroutine Get_Path_Ready ( Path, Instrument )

    ! Given a line defined by a point in ECR, and a vector in ECR parallel
    ! to that line, compute the position of the tangent on that line, and
    ! the extents of the interesting parts of the line.

    ! The given line is Path%Lines(1,1) + s * Path%Lines(2,1).
    ! Path%Lines(2,1) is made an unit vector here.  A line defined by
    ! Path%Lines(1,2) + s * Path%Lines(2,2) is produced, which is the
    ! continuation of Path%Lines(:,1) after the tangent point if
    ! Path%Lines(:,1) does not intersect the Earth reference ellipsoid, or
    ! the reflection of Path%Lines(:,1) if it does intersect the Earth
    ! reference ellipsoid.  Path%Lines(2,2) is an unit vector.

    ! The tangent or surface-intersection point is Path%Lines(1,2).

    class(path_t), intent(inout) :: Path
    type(ECR_t), intent(in), optional :: Instrument ! Position in ECR

    type(ECR_t) :: D_Inst ! Instrument - path%lines(1,1)
    real(rg) :: S_Inst    ! S-value of instrument position
    real(rg) :: Tangent   ! S-value of tangent point
    real(rg) :: Tan_Dir   ! +/- 1; Tan_Dir * Lines(2,i) is directed from
                          ! Lines(1,i) toward the tangent point.

    ! Make Path%Lines(1:2) an unit vector to simplify later calculations.
    ! Get the tangent or intersection point of Path%Lines(:,1) with the Earth
    ! reference ellipsoid.
    ! Path%Lines(1,1) + s * Path%Lines(2,1) describes the line before the
    ! tangent point.
    ! Path%Lines(1,2) + s * Path%Lines(2,2) describes the line after the
    ! tangent point.
    ! We could, in principle, divide the line into an arbitrary number of
    ! segments, to account (approximately) for refraction.
    call path_continuation ( path%lines, tangent )
    tan_dir = sign(1.0_rg,tangent)
    if ( present(instrument) ) then
      d_inst = instrument - path%lines(1,1)
      s_inst = d_inst%norm2() - tangent
    end if
    if ( tan_dir >= 0 ) then
      ! Direction of Path%Lines(2,1) is from Path%Lines(1,1) toward tangent,
      ! therefore direction of Path%Lines(2,2) is from tangent toward +infinity
      path%sMin(1) = -sqrt(huge(0.0_rg)) ! Allow intersections before Lines(1,1)
      if ( present(instrument) ) path%sMin(1) = s_inst
      path%sMax(1) = tangent             ! Disallow intersections after tangent
      path%sMin(2) = 0                   ! Disallow intersections before tangent
      path%sMax(2) = sqrt(huge(0.0_rg))  ! Allow intersections arbitrarily far away
    else ! sMax(1) < 0
      ! Direction of Lines(2,1) is from Lines(1,1) away from tangent, therefore
      ! direction of Lines(2,2) is from -infinity toward tangent
      path%sMin(1) = tangent             ! Disallow intersections before tangent
      path%sMax(1) = sqrt(huge(0.0_rg))  ! Allow intersections arbitrarily far away
      path%sMin(2) = -sqrt(huge(0.0_rg)) ! Allow intersections arbitrarily far away
      if ( present(instrument) ) path%sMin(2) = s_inst
      path%sMax(2) = 0                   ! Disallow intersections after tangent
    end if

  end subroutine Get_Path_Ready

  subroutine Path_Continuation ( Lines, Tangent )

    ! Calculate the path continuation after the tangent point or intersection.
    ! Lines(1,1) + s * Lines(2,1) describes the line before the tangent point.
    ! Lines(1,2) + s * Lines(2,2) describes the line after the tangent point.
    ! We could, in principle, divide the line into an arbitrary number of
    ! segments, to account (approximately) for refraction.

    use Geolocation_0, only: Norm2
    use Line_And_Ellipsoid_m, only: Line_And_Ellipsoid, &
      & Exact_Line_Nearest_Ellipsoid
    use Line_And_Plane_m, only: Line_Reflection
    type(ECR_t), intent(inout) :: Lines(2,2)
    real(rg), intent(out) :: Tangent   ! S-value of tangent point

    type(ECR_t) :: Grad           ! Gradient to Earth reference ellipsoid
                                  ! at tangent point
    real(rg) :: H                 ! Tangent height above Earth surface,
                                  ! Lines is an intersection if H < 0.
    integer :: I
    type(ECR_t), allocatable :: Ints(:) ! Intersections with Earth reference ellipsoid
    real(rg), allocatable :: W(:) ! Where line and ellipsoid intersect

    ! Make Lines(2,1) an unit vector, to simplify later calculations.
    lines(2,1) = lines(2,1) / lines(2,1)%norm2()

    ! Get the tangent or intersection point of Lines(:,1) with the Earth
    ! reference ellipsoid.
    call exact_line_nearest_ellipsoid ( lines(:,1), tangent, h=h )
    if ( h >= 0 ) then ! No intersection, Lines(:,2) is colinear with Lines(:,1)
      lines(1,2) = lines(1,1) + tangent * lines(2,1)
      lines(2,2) = lines(2,1) ! Parallel to incident line
    else               ! Compute reflection direction
      ! Assume lines(1,1) is not inside the ellipsoid
      call line_and_ellipsoid ( lines(:,1), ints, w )
      ! If lines(:,1) is not tangent to the Earth's surface, there are
      ! two intersections.  Choose the one closest to lines(1,1).
      i = minloc(norm2(ints - lines(1,1)),1)
      lines(1,2) = ints(i) ! The continuation starts at that intersection
      tangent = w(i)
      grad = lines(1,2)%grad() ! Gradient to Earth reference ellipsoid at
      ! the intersection Compute Lines(2,2) such that Lines(2,2) is at the
      ! same angle from Grad as Lines(2,1), but on the opposite side of the
      ! tangent from Lines(2,1).
      call line_reflection ( lines(2,1), grad, lines(2,2) )
      deallocate ( ints, w )
    end if
  end subroutine Path_Continuation

  function Tangent_Geodetic_Coordinates ( FwdModelIn, FwdModelExtra, MAF, &
                                        & Ref_MIF ) result ( Tngt_Geod )

    !{ Compute the surface geodetic coordinates of the tangent at the
    !  specified reference MIF and MAF.  The result is $\mathbf{T}_g$
    !  in wvs-143.

    use Geolocation_0, only: ECR_t, H_V_Geod
    use Intrinsic, only: L_TngtECR
    use VectorsModule, only: GetVectorQuantityByType, Vector_T, VectorValue_T

    type(vector_t), intent(in) :: FwdModelIn, FwdModelExtra
    integer, intent(in) :: MAF
    integer, intent(in) :: Ref_MIF

    type(h_v_geod) :: Tngt_Geod  ! Tangent geodetic coordinates at surface

    type(ECR_t) :: Tngt                       ! Tangent position
    type(vectorValue_t), pointer :: TngtECR   ! Tangent position

    tngtECR => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_tngtECR )
    tngt%xyz = tngtECR%value3(1:3,ref_MIF,MAF) ! Tangent ECR at ref MIF and MAF
    tngt_geod = tngt%geod ( )            ! Geodetic coordinates of ref tangent
    tngt_geod%v = 0                      ! Surface geod. coordinates of tangent

  end function Tangent_Geodetic_Coordinates

  ! Specific procedures for Union_Paths generic:

  subroutine Add_One_Path_To_Union ( Unions, Path )
    ! Construct the union of union%facets with path%facets, and the union of
    ! union%vertices with path%vertices.
    use MLSSets, only: Union
    type(facets_and_vertices_t), intent(inout) :: Unions
    type(facets_and_vertices_t), intent(in) :: Path

!   call move_alloc ( union(union%facets,path%facets), union%facets ) ! doesn't work
    if ( allocated(unions%facets) ) then
      unions%facets = union(unions%facets,path%facets)
    else
      unions%facets = path%facets
    end if
    unions%nf = size(unions%facets)
!   call move_alloc ( union(unions%vertices,path%vertices), unions%vertices ) ! doesn't work
    if ( allocated(unions%vertices) ) then
      unions%vertices = union(unions%vertices,path%vertices)
    else
      unions%vertices = path%vertices
    end if
    unions%nv = size(unions%vertices)
  end subroutine Add_One_Path_To_Union

  subroutine Add_Paths_To_Union ( Unions, Paths )
    type(facets_and_vertices_t), intent(inout) :: Unions
    type(facets_and_vertices_t), intent(in) :: Paths(:)
    integer :: I
    do i = 1, size(paths)
      call union_paths ( unions, paths(i) )
    end do
  end subroutine Add_Paths_To_Union

!=============================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Path_Representation_m

! $Log$
! Revision 2.11  2019/06/24 23:28:50  pwagner
! Updated to reflect TA-01-143
!
! Revision 2.10  2017/10/31 23:45:31  vsnyder
! Add Geometric_Path and Tangent_Geodetic_Coordinates
!
! Revision 2.9  2017/03/11 00:47:39  vsnyder
! Add Union_paths
!
! Revision 2.8  2016/11/23 00:06:55  vsnyder
! Remove Value_t in favor of types in Indexed_Values_m
!
! Revision 2.7  2016/11/12 01:32:23  vsnyder
! Add NP and NZ components to Value_t and Flag_t, and default initialize
!
! Revision 2.6  2016/11/11 01:46:41  vsnyder
! Add Facets_and_Vertices_t
!
! Revision 2.5  2016/11/08 02:04:56  vsnyder
! Use Exact_Line_Nearest_Ellipsoid
!
! Revision 2.4  2016/11/07 23:50:55  vsnyder
! Remove unused variable declaration
!
! Revision 2.3  2016/11/07 23:48:52  vsnyder
! Make Path_Continuation public
!
! Revision 2.2  2016/11/04 01:26:32  vsnyder
! Spiff some comments
!
! Revision 2.1  2016/10/26 01:20:19  vsnyder
! Initial commit
!
