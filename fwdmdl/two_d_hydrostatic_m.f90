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
  public Two_D_Hydrostatic

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

  ! ------------------------------------------  Two_D_Hydrostatic  -----
  subroutine Two_D_Hydrostatic ( Grids_tmp, z_refs, h_refs, z_grid, beta, &
                              &  t_grid, h_grid, dhidzij, dhidtlm, ddhdhdtl0 )

  ! Compute the 2 dimensional hydrostatic stuff if the horizontal grid is
  ! PhiTan based, or the 3 dimensional hydrostatic stuff if the horizontal
  ! grid is QTM based.

    use Constants, only: Deg2Rad
    use Earth_Constants, only: F2 => Earth_Axis_Ratio_Squared, & ! b^2/a^2
      & EarthRadA
    use Geometry, only: Orbit_Plane_Minor_Axis_sq
    use Hydrostatic_m, only: Hydrostatic
    use Intrinsic, only: L_Geodetic, L_QTM
    use Load_sps_data_m, ONLY: Grids_T
    use MLSKinds, only: RP, IP

  ! Inputs:

    type (Grids_T), intent(in) :: Grids_tmp   ! All Temperature's coordinates
    real(rp), intent(in) :: z_refs(:)  ! Reference pressures
    real(rp), intent(in) :: h_refs(:)  ! Reference geopotential heights (km) at
                                       ! z_refs. The horizontal basis for these
                                       ! is aligned with Grids_tmp%phi_basis.
    real(rp), intent(in) :: z_grid(:)  ! pressures for which heights and
                                       ! temperatures are needed
    real(rp), intent(in) :: beta       ! spacecraft beta angle (Radians).  Not
                                       ! used if grids_tmp%qtyStuff(1)%qty%...
                                       ! template%the_Hgrid%type == L_QTM

  ! Outputs:

    real(rp), intent(out):: t_grid(:,:)    ! Computed temperatures
    real(rp), intent(out):: h_grid(:,:)    ! Computed heights (km)
    real(rp), intent(out):: dhidzij(:,:)   ! Derivative of height wrt zeta
    real(rp), optional, intent(out):: dhidtlm(:,:,:) ! Derivative of height wrt
                                           ! temperatures on output phi grid
    real(rp), optional, intent(out):: ddhdhdtl0(:,:,:) ! second order derivative
                                           ! at the tangent only---used for
                                           ! antenna effects
  ! Internal stuff

    integer(ip) :: I, J1, J2
    integer(ip) :: P_coeffs ! Size of interesting part of Grids_tmp%phi_basis
                            ! if .not. QTM, else the size of the interesting
                            ! part of Grids_tmp%QTM_Geo
    logical :: QTM          ! ...%the_Hgrid%type == l_QTM
    integer(ip) :: Z_coeffs ! Size of interesting part of Grids_tmp%zet_basis

    real(rp) :: CSQ ! C**2
    real(rp) :: Lat         ! Geocentric latitude in Radians
    real(rp) :: SinBeta, SinPhi, SinPhiSQ ! Used to compute Lat if not QTM

  ! Begin execution

    p_coeffs = Grids_tmp%l_p(1) ! - Grids_tmp%l_p(0), which is always zero
    z_coeffs = Grids_tmp%l_z(1) ! - Grids_tmp%l_z(0), which is always zero

    QTM = associated(grids_tmp%qtyStuff(1)%qty%template%the_Hgrid)
    if ( QTM ) QTM = grids_tmp%qtyStuff(1)%qty%template%the_Hgrid%type == L_QTM
    if ( .not. QTM ) then
      ! Phitan-based two-dimensional grid

      !{ Compute the orbit-plane projected minor axis $c$, where
      !  $c^2 = \frac{a^2\,b^2}{a^2 \sin^2 \beta + b^2 \cos^2 \beta} =
      !         \frac{a^2}{\left(\frac{a^2}{b^2}-1\right) \sin^2 \beta + 1}$

      sinBeta = sin(beta)
      csq = orbit_plane_minor_axis_sq ( beta )

    end if

    ! compute the 2 d hydrostatic

    j2 = 0
    do i = 1, p_coeffs

      !{ Compute the geocentric latitude $\lambda\, = \, 
      !  \sin^{-1} ( \sin \gamma \, \sin \beta )$, where
      !  $\sin^2 \gamma = 
      !   \frac{c^4 \sin^2 \phi}{a^4 \cos^2 \phi + c^4 \sin^2 \phi}$
      !  and $\gamma$ is the geocentric angle in the orbit plane ellipse
      !  between $\mathbf{R}^{\oplus}$ and the $x$ axis.

      if ( .not. QTM ) then
        sinPhi = sin(Grids_tmp%phi_basis(i))
        sinPhiSQ = sinPhi**2
        lat = asin(csq * sinPhi * sinBeta &
          & / sqrt(earthrada**4*(1.0_rp-sinPhiSQ) + csq**2*sinPhiSQ))
      else ! QTM-based grid
        associate ( template => Grids_tmp%qtyStuff(i)%qty%template )
          ! If QTM, all of the QTM grid is used, so we don't need to get
          ! a subset of the indices that depend upon the phi window from
          ! grids_tmp.
          lat = deg2rad * template%the_HGrid%QTM_geo(i)%lat
          if ( template%latitudeCoordinate == l_geodetic ) then
            ! Similar to GeocToGeodLat from Geometry, but result is radians
            lat = atan2 ( f2 * sin(lat), cos(lat) ) ! Convert geodetic to geocentric
          end if
        end associate
      end if

      j1 = j2
      j2 = j1 + z_coeffs
      if ( present(ddhdhdtl0) ) then ! needs dhidtlm
        call hydrostatic ( lat, Grids_tmp%zet_basis, Grids_tmp%values(j1+1:j2), &
           & z_grid, z_refs(i), h_refs(i), t_grid(:,i), h_grid(:,i), &
           & dhidzij(:,i), dhidtlm(:,:,i), ddhdhdtl0(:,:,i) )
      else if ( present(dhidtlm) ) then
        call hydrostatic ( lat, Grids_tmp%zet_basis, Grids_tmp%values(j1+1:j2), &
           & z_grid, z_refs(i), h_refs(i), t_grid(:,i), h_grid(:,i), &
           & dhidzij(:,i), dhidtlm(:,:,i) )
      else
        call hydrostatic ( lat, Grids_tmp%zet_basis, Grids_tmp%values(j1+1:j2), &
           & z_grid, z_refs(i), h_refs(i), t_grid(:,i), h_grid(:,i), &
           & dhidzij(:,i) )
      end if
    end do

  end subroutine Two_D_Hydrostatic

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
