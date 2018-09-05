! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Magnetic_Field_On_Path_m

  implicit NONE

  private
  public :: Magnetic_Field_On_Path

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------- -----------------

contains

  subroutine Magnetic_Field_On_Path ( Grids_Mag, Tan_Pt_f, H_Path, R_Eq, &
                                    & Z_Path, Phi_Path, Rot, Mag_Path )

    use Comp_Eta_DoCalc_Sparse_m, only: Comp_Eta => Comp_Eta_DoCalc_Sparse
    use Comp_Sps_Path_Sparse_m, only: Comp_1_Sps_Path_Sparse_No_Frq_2D
    use GLNP, only: NGP1
    use Intrinsic, only: L_GeocAltitude, L_GeodAltitude,  L_Zeta
    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use Sparse_Eta_m, only: Sparse_Eta_t

    type(grids_t), intent(in) :: Grids_Mag
    integer, intent(in) :: Tan_Pt_f
    real(rp), intent(in) :: H_Path(:)  ! Heights on path, km above R_Eq
    real(rp), intent(in) :: R_Eq       ! Equivalent Earth radius
    real(rp), intent(in) :: Z_Path(:), Phi_Path(:)
    real(rp), intent(in) :: Rot(3,3)   ! Rotation from ECR to FOV
    real(rp), intent(out) :: Mag_Path(:,:)

    type(sparse_eta_t) :: Eta_Mag_zp
    integer :: J

    select case ( grids_mag%z_coord(1) ) ! grids_mag is a one-quantity grid
      case ( l_geocAltitude, l_geodAltitude )
        ! load_sps_data_m%Fill_Grids_2 has converted geocentric height
        ! to geodetic altitude above the geoid, and converted to km.
        ! H_Path is km above the equivalent Earth center, so subtract R_Eq
        ! to get it on the same basis as Grids_Mag%zet_basis, which is
        ! altitude in km here, not zeta.
        call comp_eta ( grids_mag, tan_pt_f, h_path-r_eq, phi_path, &
                      & eta_mag_zp, skip=ngp1 )
      case ( l_zeta )
        call comp_eta ( grids_mag, tan_pt_f, z_path, phi_path, &
                      & eta_mag_zp, skip=ngp1 )
    end select

    call Comp_1_Sps_Path_Sparse_No_Frq_2D ( grids_mag, eta_mag_zp, &
      & mag_path(:,1:3) )

    do concurrent ( j = 1: size(mag_path,1) )
      ! Rotate mag_path from ECR to IFOVPP (R1A) coordinates.
      mag_path(j,1:3) = matmul ( rot, mag_path(j,1:3) )
      ! Put the magnitude of mag_path(j,1:3) in mag_path(j,4)
      mag_path(j,4) = norm2(mag_path(j,1:3))
      ! Normalize mag_path(j,1:3).
      if ( mag_path(j,4) /= 0.0_rp ) then
        mag_path(j,1:3) = mag_path(j,1:3) / mag_path(j,4)
      else
        mag_path(j,1:3) = 0.0_rp
        mag_path(j,3) = 1.0_rp ! arbitrarily, theta=0 for zero field
      end if
    end do

  end subroutine Magnetic_Field_On_Path

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Magnetic_Field_On_Path_m

! $Log$
! Revision 2.1  2018/09/05 20:53:33  vsnyder
! Initial commit
!
