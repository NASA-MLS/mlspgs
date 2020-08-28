! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Get_Chi_Out_m

  use Sparse_m, only: Sparse_t
  implicit none
  private
  public :: Get_Chi_Out

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Get_Chi_Out ( ZetaTan, PhiTan, SCgeocAlt, Grids_tmp, Ref_zeta, &
           & Ref_gph, Elev_offset, Req, Grids_f, H2O_ind, &
           & Tan_Chi_Out, dhdz_out, dx_dh_out, dxdt_tan, d2xdxdt_tan, &
           & InstRefr, MAF, ScECR_MIF, ECRtoFOV )

! Compute the output angles to interpolate to

    use Get_Chi_Angles_m, only: Get_Chi_Angles
    use Get_Chi_Angles_3D_m, only: Get_Chi_Angles
    use Load_sps_data_m, only: Grids_T
    use MLSKinds, only: RP
    use Refraction_m, only: MaxRefraction, Refractive_Index
    use Sparse_Eta_m, only: Sparse_Eta_t
    use Two_D_Hydrostatic_m, only: Two_D_Hydrostatic
    use VectorsModule, only: VectorValue_t

! Inputs

    real(rp), intent(in) :: ZetaTan(:)   ! Tangent zeta profle for input mmaf
    real(rp), intent(in) :: PhiTan(:)    ! Tangent phi profle for input mmaf
                                         ! (in radians)
    real(rp), intent(in) :: SCgeocAlt(:) ! Spacecraft geocentric altitude
                                         ! profile in METERS
    type (grids_t), intent(in) :: Grids_tmp  ! All Temperature's coordinates

    real(rp), intent(in) :: Ref_zeta     ! Zeta for input gph's
    real(rp), intent(in) :: Ref_gph(:)   ! Reference geopotential heights along the
                           ! temperature horizontal basis in METERS
    real(rp), intent(in) :: Elev_offset  ! Elevation offset in radians
    real(rp), intent(in) :: Req(:)       ! Earth radius (M) in equivalent circle
    type(grids_t), intent(in) :: Grids_f ! Grids other than temperature
    integer, intent(in) :: H2O_ind       ! Where is H2O in Grids_F?

! Outputs

    real(rp), intent(out) :: Tan_Chi_Out(:) ! Computed tangent pointing angles
                                ! in radians, corrected for refraction
    real(rp), intent(out) :: dhdz_out(:)    ! dh/dz on the output grid
    real(rp), intent(out) :: dx_dh_out(:)   ! d(chi)/dh on the output grid

! Outputs not computed if size(dxdt_tan) is zero
    real(rp), intent(out), optional :: dxdt_tan(:,:)      ! Computed dchi dt.
    real(rp), intent(out), optional :: d2xdxdt_tan(:,:)   ! Computed d2chi dxdt.
                                   ! at each phi value representation (km)

! Optional inputs for QTM
    real(rp), intent(in), optional :: InstRefr ! Index of refraction -1 at inst.
    integer, intent(in), optional :: MAF ! Used to index ScECR_MIF and ECRtoFOV
    type(VectorValue_t), intent(in), optional :: ScECR_MIF
    type(VectorValue_t), intent(in), optional :: ECRtoFOV

! The matrix hierarchy is n_out,n_t_zeta,n_t_phi

! Internal variables

    integer :: ht_i, n_out, n_t_phi, n_t_zeta

    logical :: Derivs ! Compute dxdt_tan and d2xdxdt_tan

    real(rp) :: d2hdhdt_tan(size(zetatan), grids_tmp%l_z(1), grids_tmp%l_p(1))
    real(rp) :: dhdt_tan(size(zetatan), grids_tmp%l_z(1), grids_tmp%l_p(1))
    real(rp) :: dhdz_tan(size(zetatan), grids_tmp%l_p(1))
! ifort 17 complains about this, but gets different answers if this is deleted
real(rp) :: h2o_tan_out(merge(size(zetatan),0,h2o_ind>0))
! or if the declaration is changed to
! real(rp) :: h2o_tan_out(0)
    real(rp) :: Height_Tan(size(zetatan), grids_tmp%l_p(1))
    real(rp) :: H_Tan_Out(size(zetatan))
    real(rp) :: N_Tan_Out(size(zetatan))
    real(rp) :: Temp_Tan(size(zetatan), grids_tmp%l_p(1))
    real(rp) :: T_Tan_Out(size(zetatan))

    type(sparse_eta_t) :: Eta_T ! Interpolation coefficients from temperature
                                ! phi basis to PhiTan

! Dimensions we use

    n_out = size(zetatan)
    n_t_phi = grids_tmp%l_p(1)  ! Size of temperature phi basis
    n_t_zeta = grids_tmp%l_z(1) ! Size of temperature zeta basis

    call two_d_hydrostatic ( grids_tmp, ref_zeta, ref_gph, zetatan, &
             & temp_tan, height_tan, dhdz_tan, dHidTlm=dhdt_tan, &
             & ddHdHdTl0=d2hdhdt_tan )

! Tangent heights, temperatures, and dhdz, for tangent pressures along phi

    call eta_t%eta ( grids_tmp%phi_basis, phitan, sorted=.false. )

    do ht_i = 1, n_out
      h_tan_out(ht_i) = eta_t%row_dot_vec ( ht_i, height_tan(ht_i,:) )
      t_tan_out(ht_i) = eta_t%row_dot_vec ( ht_i, temp_tan(ht_i,:) )
      dhdz_out(ht_i) = eta_t%row_dot_vec ( ht_i, dhdz_tan(ht_i,:) )
    end do

! Compute refractive index

    if ( h2o_ind > 0 ) then
      block
        integer :: beg_p, beg_z  ! Begin of h2o info in grids_f
        integer :: end_p, end_z  ! End of h2o info in grids_f
        type(sparse_eta_t) :: Eta_P, Eta_Z, Eta_ZP
        real(rp) :: h2o_tan_out(n_out)
        ! Compute the tangent water vapor profile along input phi_tan
        beg_z = grids_f%l_z(h2o_ind-1) + 1
        end_z = grids_f%l_z(h2o_ind)
        beg_p = grids_f%l_p(h2o_ind-1) + 1
        end_p = grids_f%l_p(h2o_ind)

        call eta_p%eta ( grids_f%phi_basis(beg_p:end_p), phitan, sorted=.false. )
        call eta_z%eta ( grids_f%zet_basis(beg_z:end_z), zetatan, sorted=.false. )

        call eta_zp%eta_nd ( eta_z, eta_p )
        call eta_zp%sparse_dot_vec ( grids_f%c(h2o_ind)%v1, h2o_tan_out )

        if ( grids_f%lin_log(h2o_ind) ) h2o_tan_out = exp(h2o_tan_out)

! Compute refractive index depending on h2o

        call refractive_index ( 10.0**(-zetatan), t_tan_out, n_tan_out, &
                              & h2o_path=h2o_tan_out )
      end block

    else

! Compute refractive index not depending on h2o

      call refractive_index ( 10.0**(-zetatan), t_tan_out, n_tan_out )

    end if
    n_tan_out = min ( n_tan_out, maxRefraction )

! Compute output angles to interpolate to

    derivs = present(dxdt_tan)
    if ( derivs ) derivs = size(dxdt_tan) > 0
    if ( derivs ) then

      block
        ! Variables for one tangent phi's derivatives w.r.t.
        ! temperature zeta and phi.  Get_Chi_Angles wants derivatives at the
        ! tangent as one-D zeta*phi instead of two-D zeta X phi.
        real(rp), target :: dhdt_1_tan_1(n_t_zeta*n_t_phi)
        real(rp), pointer :: dhdt_1_tan(:,:)
        real(rp), target :: d2hdhdt_1_tan_1(n_t_zeta*n_t_phi)
        real(rp), pointer :: d2hdhdt_1_tan(:,:)
        integer :: J, K       ! Subscripts

        dhdt_1_tan_1 = 0
        dhdt_1_tan(1:n_t_zeta,1:n_t_phi) => dhdt_1_tan_1
        d2hdhdt_1_tan_1 = 0
        d2hdhdt_1_tan(1:n_t_zeta,1:n_t_phi) => d2hdhdt_1_tan_1

        do ht_i = 1, n_out
          ! Interpolate derivatives to phi for Ht_i
          ! An iterator to traverse a row of a Sparse_t would be useful here
          k = eta_t%rows(ht_i)                  ! Last element in the row
          if ( k /= 0 ) then                    ! Row isn't empty
            do
              j = eta_t%e(k)%c                  ! Element's column subscript
              dhdt_1_tan(:,j)    = dhdt_tan(ht_i,:,j)    * eta_t%e(k)%v
              d2hdhdt_1_tan(:,j) = d2hdhdt_tan(ht_i,:,j) * eta_t%e(k)%v
              k = eta_t%e(k)%nr                 ! Next element in the row
              if ( k == eta_t%rows(ht_i) ) exit ! Back to the last one?
            end do
          end if

          if ( .not. present(instRefr) ) then ! not QTM
            ! We can't compute the index of refraction at the instrument
            ! (actually index - 1.0), so use zero.
            ! SCgeocAlt and RefGPH are in meters, but Get_Chi_Angles wants them
            ! in km.
            call get_chi_angles ( 0.001_rp*scgeocalt(ht_i), n_tan_out(ht_i), 0.0_rp, &
              & h_tan_out(ht_i), phitan(ht_i), req(ht_i), elev_offset, &
              & tan_chi_out(ht_i), dx_dh_out(ht_i), dhdz_out(ht_i), &
              & dhdt_1_tan_1, d2hdhdt_1_tan_1, &
              & dxdt_tan(ht_i,:), d2xdxdt_tan(ht_i,:) )
            ! Put zeroes back where we put nonzeroes above
            ! An iterator to traverse a row of a Sparse_t would be useful here
          else ! QTM
            call get_chi_angles ( n_tan_out(ht_i), instRefr, &
              & scECR_MIF%value3(1:3,ht_i,MAF), ECRtoFOV%value3(1:3,ht_i,MAF), &
              & h_tan_out(ht_i), req(ht_i), elev_offset, &
              & tan_chi_out(ht_i), &
              & dx_dh_out(ht_i), dhdz_out(ht_i), dhdt_1_tan_1, &
              & d2hdhdt_1_tan_1, dxdt_tan(ht_i,:), d2xdxdt_tan(ht_i,:) )
          end if
          k = eta_t%rows(ht_i)                  ! Last element in the row
          if ( k /= 0 ) then                    ! Row isn't empty
            do
              j = eta_t%e(k)%c                  ! Element's column subscript
              dhdt_1_tan(:,j) = 0
              d2hdhdt_1_tan(:,j) = 0
              k = eta_t%e(k)%nr                 ! Next element in the row
              if ( k == eta_t%rows(ht_i) ) exit ! Back to the last one?
            end do
          end if
        end do

      end block

    else

      if ( .not. present(instRefr) ) then ! not QTM
        do ht_i = 1, n_out
           ! We can't compute the index of refraction at the instrument, so
           ! pretend it's 1.0 by setting the instrument index to zero.
          call get_chi_angles( 0.001_rp*scgeocalt(ht_i), n_tan_out(ht_i), 0.0_rp, &
            & h_tan_out(ht_i), phitan(ht_i), req(ht_i), elev_offset, &
            & tan_chi_out(ht_i), dx_dh_out(ht_i), dhdz_out(ht_i) )
        end do
      else ! QTM
        do ht_i = 1, n_out
          call get_chi_angles ( n_tan_out(ht_i), instRefr, &
            & scECR_MIF%value3(1:3,ht_i,MAF), ECRtoFOV%value3(1:3,ht_i,MAF), &
            & h_tan_out(ht_i), req(ht_i), elev_offset, &
            & tan_chi_out(ht_i) )
        end do
      end if

    end if

  end subroutine Get_Chi_Out

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Get_Chi_Out_m

! $Log$
! Revision 2.29  2020/08/28 21:41:58  vsnyder
! Set up to calculate chi angles for QTM
!
! Revision 2.28  2019/06/24 23:28:17  pwagner
! Updated to reflect TA-01-143
!
! Revision 2.27  2018/09/07 17:15:20  pwagner
! Code around a NAG-6.1 compiler bug
!
! Revision 2.26  2018/08/22 01:40:27  vsnyder
! Convert to sparse interpolators
!
! Revision 2.25  2016/11/11 01:52:23  vsnyder
! Change units for ScGeocAlt and Ref_GPH from km to m, and do the conversion
! from m to km here instead of in callers, to avoid the need for an array
! temp.
!
! Revision 2.24  2016/10/24 22:14:24  vsnyder
! Eliminate orbit inclination because two_d_hydrostatic no longer needs it
!
! Revision 2.23  2016/05/10 00:03:55  vsnyder
! Cannonball polishing
!
! Revision 2.22  2016/01/23 02:53:49  vsnyder
! Get MaxRefraction from refraction_m instead of geometry
!
! Revision 2.21  2009/06/23 18:26:11  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.20  2008/07/31 17:57:54  vsnyder
! Make some arguments optional
!
! Revision 2.19  2006/11/30 02:06:49  vsnyder
! Clean up some surds of ASSOCIATED->SIZE transition
!
! Revision 2.18  2006/11/30 01:29:28  vsnyder
! Clean up some surds of ASSOCIATED->SIZE transition
!
! Revision 2.17  2006/11/11 03:44:13  vsnyder
! Depointerize some variables
!
! Revision 2.16  2005/12/22 20:53:29  vsnyder
! Simplify calcuation of water vapor
!
! Revision 2.15  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.14  2004/05/17 23:25:04  livesey
! More guards
!
! Revision 2.13  2004/05/17 22:06:07  livesey
! Guard against large refractive indices
!
! Revision 2.12  2003/05/05 23:00:25  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.11.2.2  2003/03/20 01:42:26  vsnyder
! Revise Grids_T structure
!
! Revision 2.11.2.1  2003/03/07 23:17:06  vsnyder
! Use ASSOCIATED instead of PRESENT of GET_CHI_OUT
!
! Revision 2.11  2002/10/08 17:08:04  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.10  2002/10/03 22:25:13  vsnyder
! Move USE statements from module scope to procedure scope.  Convert some
! allocated arrays to automatic arrays.  Cosmetic changes.
!
! Revision 2.9  2002/07/11 20:51:49  bill
! made req an mmif quantity
!
! Revision 2.8  2002/07/05 07:52:48  zvi
! Coor. switch (phi,z) -> (z,phi)
!
! Revision 2.7  2002/06/28 11:06:51  zvi
! computes dx_dh & dh_dz on output grid as well
!
! Revision 2.6  2002/06/24 21:11:24  zvi
! Adding Grids_tmp structure and modifying calling sequences
!
! Revision 1.0  2002/06/24 22:20:45  bill
! Intitial release ...
