! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Two_D_Hydrostatic_m

  implicit none

  private
  public :: Two_D_Hydrostatic
  public :: Two_D_Hydrostatic_Coherent, Two_D_Hydrostatic_General

  interface Two_D_Hydrostatic
    module procedure Two_D_Hydrostatic_Coherent, Two_D_Hydrostatic_General
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

  ! ---------------------------------  Two_D_Hydrostatic_Coherent  -----
  subroutine Two_D_Hydrostatic_Coherent ( Grids_tmp, Z_Ref, H_Refs, Z_Grid, &
                               & T_Grid, H_Grid, dHidZij, &
                               & Eta_ZZ, Vertices, dHidTlm, ddHdHdTl0 )

    ! Compute the 2 dimensional hydrostatic stuff.  Assumes all the
    ! RefGPH are at the same Zeta.

    use Constants, only: Deg2Rad
    use Geometry, only: GeodToGeocLat
    use Hydrostatic_m, only: Hydrostatic
    use Intrinsic, only: L_Geodetic
    use Load_sps_data_m, only: Grids_T
    use MLSKinds, only: RP, IP
    use Sparse_m, only: Sparse_t

    ! Inputs:

    type (Grids_T), intent(in) :: Grids_tmp   ! All Temperature's coordinates
                                       ! and values.
    real(rp), intent(in) :: Z_Ref      ! Reference pressure for all profiles.
    real(rp), intent(in) :: H_Refs(:)  ! Reference geopotential heights (m) at
                                       ! z_ref for all profiles.
    real(rp), intent(in) :: Z_Grid(:)  ! pressures for which heights and
                                       ! temperatures are needed.

    ! Outputs:

    real(rp), intent(out):: T_Grid(:,:)    ! Computed temperatures, interpolated
                                           ! from Grids_Tmp to Z_Grid.
    real(rp), intent(out):: H_Grid(:,:)    ! Computed heights (km).
    real(rp), intent(out):: dHidZij(:,:)   ! Derivative of height wrt zeta.

    ! Optional inputs

    class(sparse_t), optional, intent(in) :: Eta_ZZ ! Interpolation
                                           ! coefficients from
                                           ! Grids_tmp%Zet_Basis to Z_Grid
    integer, optional, intent(in) :: Vertices(:) ! to select a subset of
                                           ! the horizontal basis for Grids_Tmp,
                                           ! usually for QTM.

    ! Optional outputs

    real(rp), optional, intent(out):: dHidTlm(:,:,:) ! Derivative of height wrt
                                           ! temperatures on output phi grid.
    real(rp), optional, intent(out):: ddHdHdTl0(:,:,:) ! second order derivative
                                           ! at the tangent only---used for
                                           ! antenna effects.

    ! Internal stuff

    integer(ip) :: I, J, N
    integer(ip) :: P_Coeffs ! Size of interesting part of Grids_tmp%phi_basis
    logical :: QTM          ! ...%the_Hgrid%type == l_QTM
    integer(ip) :: Z_Coeffs ! Size of interesting part of Grids_tmp%zet_basis
    real(rp), pointer :: T(:,:) ! Rank-2 view of Grids_tmp%values

    real(rp) :: Lat         ! Geocentric latitude in Radians

    ! Begin execution

    p_coeffs = Grids_tmp%l_p(1) ! - Grids_tmp%l_p(0), which is always zero
    z_coeffs = Grids_tmp%l_z(1) ! - Grids_tmp%l_z(0), which is always zero
    QTM = grids_tmp%isQTM(1)

    ! Compute the 2 d hydrostatic by computing the 1 d hydrostatic at
    ! each profile.

    t(1:z_coeffs,1:p_coeffs) => Grids_tmp%values ! Get rank-2 view

    n = p_coeffs
    if ( present(vertices) ) n = size(vertices)

    associate ( template => Grids_tmp%qtyStuff(1)%qty%template )
      do i = 1, n

        j = i
        if ( present(vertices) ) j = vertices(i)

        ! Compute the geocentric latitude in radians.

        if ( .not. QTM ) then
          lat = template%geodLat(1,j)
        else
          lat = template%the_HGrid%QTM_tree%geo_in(j)%lat
        end if

        if ( template%latitudeCoordinate == l_geodetic ) then
          ! Result is in Radians even though argument is in Degrees
          lat = GeodToGeocLat ( lat )
        else ! Latitude is geocentric, but in degrees, and we want radians
          lat = lat * deg2rad
        end if

        if ( present(eta_zz) ) then
          if ( present(ddhdhdtl0) ) then ! needs dhidtlm
            call hydrostatic ( lat, Grids_tmp%zet_basis, t(:,j), &
               & z_grid, z_ref, 0.001_rp*h_refs(j), eta_zz, &
               & t_grid(:,i), h_grid(:,i), dhidzij(:,i), &
               & dhidtlm(:,:,i), ddhdhdtl0(:,:,i) )
          else if ( present(dhidtlm) ) then
            call hydrostatic ( lat, Grids_tmp%zet_basis, t(:,j), &
               & z_grid, z_ref, 0.001_rp*h_refs(j), eta_zz, &
               & t_grid(:,i), h_grid(:,i), dhidzij(:,i), &
               & dhidtlm(:,:,i) )
          else
            call hydrostatic ( lat, Grids_tmp%zet_basis, t(:,j), &
               & z_grid, z_ref, 0.001_rp*h_refs(j), eta_zz, &
               & t_grid(:,i), h_grid(:,i), dhidzij(:,i) )
          end if
        else
          if ( present(ddhdhdtl0) ) then ! needs dhidtlm
            call hydrostatic ( lat, Grids_tmp%zet_basis, t(:,j), &
               & z_grid, z_ref, 0.001_rp*h_refs(j), t_grid(:,i), h_grid(:,i), &
               & dhidzij(:,i), dhidtlm(:,:,i), ddhdhdtl0(:,:,i) )
          else if ( present(dhidtlm) ) then
            call hydrostatic ( lat, Grids_tmp%zet_basis, t(:,j), &
               & z_grid, z_ref, 0.001_rp*h_refs(j), t_grid(:,i), h_grid(:,i), &
               & dhidzij(:,i), dhidtlm(:,:,i) )
          else
            call hydrostatic ( lat, Grids_tmp%zet_basis, t(:,j), &
               & z_grid, z_ref, 0.001_rp*h_refs(j), t_grid(:,i), h_grid(:,i), &
               & dhidzij(:,i) )
          end if
        end if
      end do
    end associate

  end subroutine Two_D_Hydrostatic_Coherent

  ! ----------------------------------  Two_D_Hydrostatic_General  -----
  subroutine Two_D_Hydrostatic_General ( Grids_tmp, Z_Refs, H_Refs, Z_Grid, &
                               & T_Grid, H_Grid, dHidZij, &
                               & Eta_ZZ, Vertices, dHidTlm, ddHdHdTl0 )

    ! Compute the 2 dimensional hydrostatic stuff.  Doesn't assume all the
    ! RefGPH are at the same Zeta, unless size(Z_Refs) == 1.

    use Constants, only: Deg2Rad
    use Geometry, only: GeodToGeocLat
    use Hydrostatic_m, only: Hydrostatic
    use Intrinsic, only: L_Geodetic
    use Load_sps_data_m, only: Grids_T
    use MLSKinds, only: RP, IP
    use Sparse_m, only: Sparse_t

    ! Inputs:

    type (Grids_T), intent(in) :: Grids_tmp   ! All Temperature's coordinates
                                       ! and values.
    real(rp), intent(in) :: Z_Refs(:)  ! Reference pressures for all profiles.
    real(rp), intent(in) :: H_Refs(:)  ! Reference geopotential heights (m) at
                                       ! z_refs for all profiles.
    real(rp), intent(in) :: Z_Grid(:)  ! pressures for which heights and
                                       ! temperatures are needed.

    ! Outputs:

    real(rp), intent(out):: T_Grid(:,:)    ! Computed temperatures, interpolated
                                           ! from Grids_Tmp to Z_Grid.
    real(rp), intent(out):: H_Grid(:,:)    ! Computed heights (km).
    real(rp), intent(out):: dHidZij(:,:)   ! Derivative of height wrt zeta.

    ! Optional inputs

    class(sparse_t), optional, intent(in) :: Eta_ZZ ! Interpolation
                                           ! coefficients from
                                           ! Grids_tmp%Zet_Basis to Z_Grid
    integer, optional, intent(in) :: Vertices(:) ! to select a subset of
                                           ! the horizontal basis for Grids_Tmp,
                                           ! usually for QTM.

    ! Optional outputs

    real(rp), optional, intent(out):: dHidTlm(:,:,:) ! Derivative of height wrt
                                           ! temperatures on output phi grid.
    real(rp), optional, intent(out):: ddHdHdTl0(:,:,:) ! second order derivative
                                           ! at the tangent only---used for
                                           ! antenna effects.

    ! Internal stuff

    integer(ip) :: I, J, K, N
    integer(ip) :: P_Coeffs ! Size of interesting part of Grids_tmp%phi_basis
    logical :: QTM          ! ...%the_Hgrid%type == l_QTM
    integer(ip) :: Z_Coeffs ! Size of interesting part of Grids_tmp%zet_basis
    real(rp), pointer :: T(:,:) ! Rank-2 view of Grids_tmp%values

    real(rp) :: Lat         ! Geocentric latitude in Radians

    ! Begin execution

    p_coeffs = Grids_tmp%l_p(1) ! - Grids_tmp%l_p(0), which is always zero
    z_coeffs = Grids_tmp%l_z(1) ! - Grids_tmp%l_z(0), which is always zero
    QTM = grids_tmp%isQTM(1)

    ! Compute the 2 d hydrostatic by computing the 1 d hydrostatic at
    ! each profile.

    t(1:z_coeffs,1:p_coeffs) => Grids_tmp%values ! Get rank-2 view

    n = p_coeffs
    if ( present(vertices) ) n = size(vertices)

    associate ( template => Grids_tmp%qtyStuff(1)%qty%template )
      do i = 1, n

        j = i
        if ( present(vertices) ) j = vertices(i)

        ! Compute the geocentric latitude in radians.

        if ( .not. QTM ) then
          lat = template%geodLat(1,j)
        else
          lat = template%the_HGrid%QTM_tree%geo_in(j)%lat
        end if

        if ( template%latitudeCoordinate == l_geodetic ) then
          ! Result is in Radians even though argument is in Degrees
          lat = GeodToGeocLat ( lat )
        else ! Latitude is geocentric, but in degrees, and we want radians
          lat = lat * deg2rad
        end if

        k = min(i,ubound(z_refs,1))
        if ( present(eta_zz) ) then
          if ( present(ddhdhdtl0) ) then ! needs dhidtlm
            call hydrostatic ( lat, Grids_tmp%zet_basis, t(:,j), &
               & z_grid, z_refs(k), 0.001_rp*h_refs(j), eta_zz, &
               & t_grid(:,i), h_grid(:,i), dhidzij(:,i), &
               & dhidtlm(:,:,i), ddhdhdtl0(:,:,i) )
          else if ( present(dhidtlm) ) then
            call hydrostatic ( lat, Grids_tmp%zet_basis, t(:,j), &
               & z_grid, z_refs(k), 0.001_rp*h_refs(j), eta_zz, &
               & t_grid(:,i), h_grid(:,i), dhidzij(:,i), &
               & dhidtlm(:,:,i) )
          else
            call hydrostatic ( lat, Grids_tmp%zet_basis, t(:,j), &
               & z_grid, z_refs(k), 0.001_rp*h_refs(j), eta_zz, &
               & t_grid(:,i), h_grid(:,i), dhidzij(:,i) )
          end if
        else
          if ( present(ddhdhdtl0) ) then ! needs dhidtlm
            call hydrostatic ( lat, Grids_tmp%zet_basis, t(:,j), &
               & z_grid, z_refs(k), 0.001_rp*h_refs(j), t_grid(:,i), h_grid(:,i), &
               & dhidzij(:,i), dhidtlm(:,:,i), ddhdhdtl0(:,:,i) )
          else if ( present(dhidtlm) ) then
            call hydrostatic ( lat, Grids_tmp%zet_basis, t(:,j), &
               & z_grid, z_refs(k), 0.001_rp*h_refs(j), t_grid(:,i), h_grid(:,i), &
               & dhidzij(:,i), dhidtlm(:,:,i) )
          else
            call hydrostatic ( lat, Grids_tmp%zet_basis, t(:,j), &
               & z_grid, z_refs(k), 0.001_rp*h_refs(j), t_grid(:,i), h_grid(:,i), &
               & dhidzij(:,i) )
          end if
        end if
      end do
    end associate

  end subroutine Two_D_Hydrostatic_General

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Two_D_Hydrostatic_m
!---------------------------------------------------
! $Log$
! Revision 2.32  2018/05/14 23:32:40  vsnyder
! Change to sparse eta representation
!
! Revision 2.31  2017/01/14 02:57:11  vsnyder
! Make Eta_ZZ polymorphic
!
! Revision 2.30  2016/12/02 02:04:50  vsnyder
! Use 'P' Eta list for Eta_ZZ
!
! Revision 2.29  2016/11/23 20:11:35  vsnyder
! Use Hydrostatic_All_ZZ if Eta_ZZ is present
!
! Revision 2.28  2016/11/11 01:50:02  vsnyder
! Add Vertices argument.  Change H_Refs units from km to m, and do the
! conversion to km for the 1-d model here instead of in callers, so as
! not to reqire an array temp.
!
! Revision 2.27  2016/10/24 22:15:11  vsnyder
! Use geodLat component instead of computing it from phi & orbIncline
!
! Revision 2.26  2016/08/23 00:43:11  vsnyder
! Components within or adjacent to the polygon are now within the QTM_Tree_t
! structure instead of the HGrid_t structure.
!
! Revision 2.25  2016/06/03 23:40:55  vsnyder
! Convert geodetic to geocentric latitude if needed
!
! Revision 2.24  2016/05/12 15:21:13  pwagner
! Avoid referring to unassociated the_HGrid
!
! Revision 2.23  2016/05/10 00:09:02  vsnyder
! Compute hydrostatic equilibrium on profiles at QTM vertices
!
! Revision 2.22  2015/04/11 00:45:03  vsnyder
! Add units (km) in h_grid comment
!
! Revision 2.21  2015/03/28 02:12:13  vsnyder
! Use Orbit_Plane_Minor_Axis_sq from Geometry
!
! Revision 2.20  2014/09/05 21:28:07  vsnyder
! Cannonball polishing
!
! Revision 2.19  2013/06/12 02:33:37  vsnyder
! Cruft removal
!
! Revision 2.18  2009/06/23 18:26:11  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.17  2009/05/13 20:03:02  vsnyder
! Get constants from Constants, kinds from MLSKinds
!
! Revision 2.16  2007/01/17 23:51:00  vsnyder
! Make dhidtlm optional
!
! Revision 2.15  2006/09/28 21:54:33  vsnyder
! Remove unused symbols
!
! Revision 2.14  2006/09/28 21:00:47  vsnyder
! Improved computation of csq again
!
! Revision 2.13  2005/12/22 20:59:18  vsnyder
! Improved computation of csq
!
! Revision 2.12  2005/06/22 18:08:20  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.11  2003/05/05 23:00:26  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.10.2.1  2003/03/20 01:42:26  vsnyder
! Revise Grids_T structure
!
! Revision 2.10  2002/10/10 19:52:45  vsnyder
! Get rid of several array temps.  Cosmetic changes.
!
! Revision 2.9  2002/10/08 17:08:06  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.8  2002/09/26 20:14:01  vsnyder
! Get PI from Units module
!
! Revision 2.7  2002/09/25 22:55:12  vsnyder
! Move USE statements from module scope to procedure scope.  Convert
! allocatable arrays to automatic arrays.  Cosmetic changes.
!
! Revision 2.6  2002/07/05 07:52:53  zvi
! Coor. switch (phi,z) -> (z,phi)
!
! Revision 2.5  2002/06/24 21:11:25  zvi
! Adding Grids_tmp stracture and modifying calling sequences
!
! Revision 2.2  2002/06/07 14:59:00  bill
! fixed latitude calculation--wgr
!
! Revision 2.1  2002/02/02 11:20:08  zvi
! Some cosmetic changes
!
! Revision 2.0  2001/09/17 20:26:28  livesey
! New forward model
!
! Revision 1.1.2.3  2001/09/13 22:51:25  zvi
! Separating allocation stmts
!
! Revision 1.1.2.2  2001/09/12 21:38:55  zvi
! Added CVS stuff
