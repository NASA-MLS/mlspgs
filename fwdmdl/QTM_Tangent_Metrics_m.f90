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

    use Geolocation_0, only: H_t
    use Generate_QTM_m, only: QTM_Tree_t
    use MLSKinds, only: RP
    use QTM_Interpolation_Weights_m, only: Weight_t, QTM_Interpolation_Weights

    ! Inputs:
    class(h_t), intent(in) :: Tan_Pt   ! Tangent point horizontal coordinates
    type(QTM_Tree_t), intent(in) :: QTM   ! Horizontal grid
    real(rp), intent(in) :: H_Ref(:,:) ! Heights by fine zeta and QTM grid
                                       ! serial number, km
    integer, intent(in) :: Tan_Ind_f   ! First (vertical) index of tangent
                                       ! point in H_Ref, which is on the fine
                                       ! zeta grid.

    ! Outputs:
    real(rp), intent(out) :: H_Surf    ! Height of the pressure reference
                                       ! surface at Tan_Pt, km -- interpolated
                                       ! in Surf_Height if present(Surf_Height)
                                       ! else interpolated in row 1 of H_REF, km
    real(rp), intent(out) :: H_Tan     ! Tangent height above H_Surf (negative
                                       ! for Earth-intersecting rays).
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

    type(Weight_t) :: Eta(3)           ! Interpolation coefficients

    ! Get interpolation coefficients and QTM serial numbers from QTM to Tan_Pt
    call QTM_Interpolation_Weights ( QTM, Tan_Pt, Eta )

    if ( present(surf_height) ) then
      ! We set the surface reference at the actual surface height if we
      ! have it, and adjust r_eq and h_tan relative to this, and adjust
      ! h_ref accordingly.
      h_surf = dot_product ( surf_height(eta%which), eta%weight )
    else
      ! If we don't have the actual surface height, we set the surface
      ! reference at the input z_ref and adjust r_eq and h_tan relative
      ! to this, and adjust h_ref accordingly.
      h_surf = dot_product ( h_ref(1,eta%which), eta%weight )
    end if

    if ( present(tan_press) .and. present(surf_temp) .and. &
      &  .not. present(surf_height) ) then
      ! Earth intersecting ray. Compute GP height (km) of tangent pressure
      ! below surface. This will be negative because tan_press < z_ref.
      ! present(tan_press) requires present(surf_temp).  We don't need to
      ! subtract h_surf here because this gives km from the z_ref surface.
      h_tan = dot_product ( surf_temp(eta%which),eta%weight ) * &
            & (tan_press-z_ref)/14.8
    else
      h_tan = dot_product ( h_ref(tan_ind_f,eta%which), eta%weight ) - h_surf
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
! Revision 2.2  2016/09/02 00:22:19  vsnyder
! Calculate H_Surf and H_Tan if they're requested
!
! Revision 2.1  2016/08/24 22:59:43  vsnyder
! Initial commit
!
