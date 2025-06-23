! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module QTM_Tangent_Metrics_m

  implicit NONE
  private

  public :: QTM_Tangent_Metrics

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  ! ----------------------------------------  QTM_Tangent_Metrics  -----
                                   ! Inputs
  subroutine QTM_Tangent_Metrics ( Tan_Pt, QTM, H_Ref, Tan_Ind_F, &
                                   ! Outputs
    &                              H_Surf, H_Tan, &
                                   ! Optional inputs
    &                              Z_Ref, Tan_press, Surf_temp, Surf_height )

  ! Compute the surface height and the tangent height at Tan_Pt.

    use Geolocation_0, only: H_v_geod !, RG
    use Generate_QTM_m, only: QTM_Tree_t
    use MLSKinds, only: RP
    use QTM_Interpolation_Weights_m, only: QTM_Interpolation_Weights
    use QTM_Interpolation_Weights_m, only: Value_QTM_2D_list_t

    ! Inputs:
    class(h_v_geod), intent(in) :: Tan_Pt ! Tangent point coordinates
                                       ! (lon,lat, ht) degrees, m
    type(QTM_Tree_t), intent(in) :: QTM   ! Horizontal grid
    real(rp), intent(in) :: H_Ref(:,:) ! Geodetic heights (above Earth Reference
                                       ! Ellipsoid) by fine zeta and index of
                                       ! vertex near the path, km
    integer, intent(in) :: Tan_Ind_f   ! First (vertical) index of tangent
                                       ! point in H_Ref, which is on the fine
                                       ! zeta grid.

    ! Outputs:
    real(rp), intent(out) :: H_Surf    ! Geodetic height of the pressure
                                       ! reference surface at Tan_Pt, km --
                                       ! interpolated in Surf_Height if
                                       ! present(Surf_Height) else interpolated
                                       ! in row 1 of H_REF, km
    real(rp), intent(out) :: H_Tan     ! Tangent height above H_Surf (negative
                                       ! for Earth-intersecting rays), km.
                                       ! Interpolated in Tan_Ind_F row of H_Ref.

    ! optional inputs
    ! We have Z_Ref, Tan_Press and Surf_Temp if the pointing is below the
    ! surface index (where zeta Z_Ref) and we don't have Surf_Height.
    ! Otherwise we have Surf_Height.
    real(rp), optional, intent(in) :: Z_Ref        ! Zeta corresponding to H_Ref(1,:)
    real(rp), optional, intent(in) :: Tan_press    ! Tangent pressure
    real(rp), optional, intent(in) :: Surf_temp(:) ! Surface temperature on QTM grid.
    real(rp), optional, intent(in) :: Surf_height(:) ! Surface height in km
                                       ! above mean sea level (whatever that
                                       ! means) on QTM grid.

!     type(value_QTM_2D_list_t(rg)) :: Eta ! Interpolation coefficients
    type(value_QTM_2D_list_t) :: Eta   ! Interpolation coefficients

    ! Get interpolation coefficients and QTM serial numbers from QTM to Tan_Pt.
    call QTM_Interpolation_Weights ( QTM, Tan_Pt%h_t, Eta )

    if ( present(surf_height) ) then
      ! We set the surface reference at the actual surface height if we
      ! have it, and adjust r_eq and h_tan relative to this, and adjust
      ! h_ref accordingly.
      h_surf = dot_product ( surf_height(eta%v%j), eta%v%v )
    else
      ! If we don't have the actual surface height, we set the surface
      ! reference at the input z_ref and adjust r_eq and h_tan relative to
      ! this, and adjust h_ref accordingly.  H_Ref is adjacent to the path,
      ! so index it using the inverse of f_and_v%vertices, which would be
      ! qtm%path_vertices(eta%n) if we hadn't already computed it and put it
      ! into eta%np..
      h_surf = dot_product ( h_ref(1,eta%v%jp), eta%v%v )
    end if

    if ( present(tan_press) .and. present(surf_temp) .and. present(z_ref) .and. &
      &  .not. present(surf_height) ) then
      ! Earth intersecting ray. Compute GP height (km) of tangent pressure
      ! below surface. This will be negative because tan_press < z_ref.
      ! present(tan_press) requires present(surf_temp) and present(z_ref). 
      ! We don't need to subtract h_surf here because this gives km from
      ! the z_ref surface.
      h_tan = dot_product ( surf_temp(eta%v%j),eta%v%v ) * &
            & (tan_press-z_ref)/14.8
    else
      ! H_Ref is adjacent to the path, so index it using the inverse of
      ! f_and_v%vertices, which would be qtm%path_vertices(eta%n) if we
      ! hadn't already computed it and put it into eta%np.
      h_tan = dot_product ( h_ref(tan_ind_f,eta%v%jp), eta%v%v ) - h_surf
    end if

  end subroutine QTM_Tangent_Metrics

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module QTM_Tangent_Metrics_m

! $Log$
! Revision 2.9  2018/05/14 23:45:19  vsnyder
! Move Hessians stuff to Hessians_m
!
! Revision 2.8  2017/08/28 20:28:08  livesey
! Changed the n,nf,np,nz elements to j,jf,...
!
! Revision 2.7  2016/11/23 00:12:28  vsnyder
! Use types from Indexed_Values_m.
!
! Revision 2.6  2016/11/17 01:48:49  vsnyder
! Use eta%np instead of qtm%path_vertices(eta%n)
!
! Revision 2.5  2016/11/12 01:38:18  vsnyder
! Use subscript of vertex near path instead of QTM serial number for H_Ref
!
! Revision 2.4  2016/11/11 01:58:36  vsnyder
! Make computation involving Z_Ref contingent upon Z_Ref being present.
!
! Revision 2.3  2016/11/02 22:50:36  vsnyder
! Use Value_t from Path_Representation instead of Weight_t from
! QTM_Interpolation_Weights_m
!
! Revision 2.2  2016/09/02 00:22:19  vsnyder
! Calculate H_Surf and H_Tan if they're requested
!
! Revision 2.1  2016/08/24 22:59:43  vsnyder
! Initial commit
!
