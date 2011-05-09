! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Add_Points_m

  implicit NONE

  private

  public :: Add_Points

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Add_Points ( more_h, more_phi, more_zeta, min_z, &
    &                     z_glgrid, nz_if, z_coarse,          &
    &                     h_path, phi_path, Vert_Inds, &
    &                     npc, npf, tan_pt_c, tan_pt_f )

    ! Add points to z_coarse, h_path and phi_path.

    use GLNP, only: GX, NGNEW, NGP1
    use MLSKINDS, only: RP
    use MLSSTRINGLISTS, only: SWITCHDETAIL
    use OUTPUT_M, only: OUTPUT
    use TOGGLES, only: SWITCHES

    real(rp), intent(in) :: More_H(:)       ! Heights to add
    real(rp), intent(in) :: More_Phi(:)     ! Phis to add
    real(rp), intent(in) :: More_Zeta(:)    ! Zetas to add
    integer, intent(in) :: Min_Z            ! Index of new minimum zeta in
                                            ! More_Zeta if there is one, else <= 0.
    real(rp), intent(in) :: Z_GLgrid(:)     ! Reference grid
    integer, intent(in) :: NZ_IF            ! Size(Z_GLGrid)
    real(rp), intent(inout) :: Z_Coarse(:)  ! Zeta on the coarse grid
    real(rp), intent(inout) :: H_Path(:)    ! Height on the fine grid
    real(rp), intent(inout) :: Phi_Path(:)  ! Phi on the fine grid
    integer, intent(inout) :: Vert_Inds(:)  ! Indices of fine grid in H_Glgrid etc.
    integer, intent(inout) :: NPC, NPF      ! Number of points in coarse, fine grids
    integer, intent(inout) :: Tan_Pt_C, Tan_Pt_F ! Tangent in coarse, fine grids

    real(rp), parameter :: Eta(ngnew) = 0.5*(1.0+gx) ! Gauss points on 0..1
    real(rp), parameter :: Phi_Tol = 0.25 * gx(1)**2 ! First GL point

    integer :: I, J, K
    integer :: New   ! Point in the fine path
    integer :: New_C ! Coarse point in the fine path
    logical :: print_more_points
    integer :: S     ! Direction of change in vert_inds at insertion point, +/-1

    print_more_points = switchDetail(switches, 'ZMOR' ) > -1

    do i = size(more_h), 1, -1 ! min zeta will be last if it's here at all
      if ( npc >= size(z_coarse) ) exit ! No room for any more
      new = minloc(abs(more_phi(i)-phi_path),1)
      if ( more_phi(i) > phi_path(new) ) new = new + 1
      ! Here phi_path(new) <= more_phi(i) <= phi_path(new+1)
      ! Find the nearest previous coarse grid point.
      new_c = new - mod(new-1,ngp1) + 1
      if ( new_c > tan_pt_f ) new_c = new_c + 1 ! don't add points between the
                                                ! tangent points
      if ( min(abs(more_phi(i)-phi_path(new)), &
        &      abs(more_phi(i)-phi_path(new+1))) <= &
        &  phi_tol * (phi_path(new_c+ngp1)-phi_path(new_c)) ) then
        ! New point is very close to an existing point
        if ( abs(more_phi(i)-phi_path(new)) > abs(more_phi(i)-phi_path(new+1)) ) &
          & new = new + 1
        ! "New" is now the index of the point
        if ( i == min_z .and. (new == new_c .or. new == new_c + ngp1) ) then
          ! "New" is at a coarse grid point
          new_c = new
          if ( i == min_z ) then
            ! "New" is for a new minimum zeta.
            ! All we need to do is change the zeta.
            j = new / ngp1 + 1
            z_coarse(j) = more_zeta(i)
            if ( print_more_points ) then
              call output ( j, before='Changed z_coarse(' )
              call output ( z_coarse(j), before=') to new minimum zeta ', advance='yes' )
            end if
          end if
          new_c = 0 ! Indicate nothing more to do
        end if
      end if
      if ( new_c > 0 .and. new_c < npf ) then
        ! New point is not near an existing point: add one to the coarse path
        npc = npc + 1
        npf = npf + ngp1
        j = new_c / ngp1 + 1 ! Point in coarse path
        z_coarse(j+1:npc) = z_coarse(j:npc-1)               ! Make room
        h_path(new_c+ngp1:npf) = h_path(new_c:npf-ngp1)
        phi_path(new_c+ngp1:npf) = phi_path(new_c:npf-ngp1)
        vert_inds(new_c+ngp1:npf) = vert_inds(new_c:npf-ngp1)
        z_coarse(j) = more_zeta(i)     ! Insert the new coarse point
        h_path(new_c) = more_h(i)
        phi_path(new_c) = more_phi(i)
        h_path(new_c+1:new_c+ngnew) = h_path(new_c) + & ! and fine point
          & eta * (h_path(new_c+ngp1)-h_path(new_c))
        phi_path(new_c+1:new_c+ngnew) = phi_path(new_c) + &
          & eta * (phi_path(new_c+ngp1)-phi_path(new_c))
        h_path(new_c+ngp1+1:new_c+ngp1+ngnew) = h_path(new_c+ngp1) + &
          & eta * (h_path(new_c+2*ngp1)-h_path(new_c+ngp1))
        phi_path(new_c+ngp1+1:new_c+ngp1+ngnew) = phi_path(new_c+ngp1) + &
          & eta * (phi_path(new_c+2*ngp1)-phi_path(new_c+ngp1))
!         if ( print_more_points ) then
          call output ( 'Added new ' )
          if ( i == min_z ) call output ( 'minimum ' )
          call output ( j, before='z_coarse(' )
          call output ( z_coarse(j), before=') = ' )
          call output ( new_c, before=', new_c = ' )
          call output ( more_phi(i), before=' at Phi = ' )
          call output ( more_h(i), before=', H = ', advance='yes' )
!         end if
        j = minval(abs(more_zeta(i)-z_glgrid(:nz_if)))
        ! For now, assume more_zeta(i) == z_glgrid(j)
        vert_inds(new_c) = j
        if ( new_c > 1 ) then
          s = sign(1,vert_inds(new_c-1) - vert_inds(new_c))
          vert_inds(new_c-ngnew:new_c-1) = (/(vert_inds(new_c-ngp1)-s*k,k=1,ngnew)/)
        end if
        if ( new_c < npf-ngp1 ) then
          s = sign(1,vert_inds(new_c) - vert_inds(new_c+1))
          vert_inds(new_c+1:new_c+ngnew) = (/(vert_inds(new_c)-s*k,k=1,ngnew)/)
        end if
        if ( new_c < tan_pt_f ) then
          tan_pt_c = tan_pt_c + 1
          tan_pt_f = tan_pt_f + ngp1
        end if
      end if
    end do ! i

  end subroutine Add_Points

!------------------------------------------------------------------

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Add_Points_m

! $Log$
! Revision 2.4  2011/05/09 17:26:50  pwagner
! Converted to using switchDetail
!
! Revision 2.3  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.2  2007/06/08 22:05:57  vsnyder
! More work on min zeta
!
! Revision 2.1  2007/02/01 02:44:29  vsnyder
! Initial commit
!
