module I_AND_K_STAR_M
  use FOV_CONVOLVE_M, only: FOV_CONVOLVE
  use HYDROSTATIC_INTRP, only: GET_PRESSURES
  use L2PC_FILE_PARAMETERS, only: DEG2RAD, MAX_NO_SV_COMPONENTS, &
                                  MSVD => max_no_sv_derivatives, &
                                  MXCO => max_no_elmnts_per_sv_component, &
                                  RAD2DEG
  use L2PC_FILE_STRUCTURES, only: L2PC_HEADER_ONE
  use L2PC_PFA_STRUCTURES, only: ATMOS_COMP, GEOM_PARAM, GEOPHYS_PARAM, &
                                 LIMB_PRESS
  use L2PCDIM, only: MNP => max_no_phi, NLVL
  use MLSCommon, only: I4, R4
  use S_CSINTERP_M, only: CSINTERP
  use S_CSPLINE_M, only: CSPLINE
  use S_LINTRP_M, only: LINTRP
  use SCSPLINE_DER_M, only: CSPLINE_DER
  implicit NONE
  private
  public :: I_AND_K_STAR

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

! This subroutine transfers the derivatives over from the internal
! convolution grid to the users specified points. This module uses
! cubic spline interpolation to do the job.
!
  Subroutine I_AND_K_STAR ( header1, ptg_press, geometric, geophysic,     &
 &           atmospheric, no_geom, no_geophys, no_sps_tbl, sps_tbl,       &
 &           conv_press, conv_hts, conv_temp, ptg_hts, ptg_angles, dx_dt, &
 &           d2x_dxdt, band, center_angle, fft_pts, fft_press,            &
 &           i_raw, In_Dir, Aaap, k_star_geophys, k_star_atmos, a_grid,   &
 &           t_grid, z_grid, N_lvls, no_conv_hts, c_yaw, s_yaw, c_roll,   &
 &           s_roll, c_pitch, s_pitch, elev_183, elev_205, azim_183,      &
 &           azim_205, azim_ref, geocsrad, i_star_all, k_star_all,        &
 &           fft_index, no_phi_g, no_phi_f, Ier )
!
    integer(i4), intent(in) :: NO_GEOM
    integer(i4), intent(in) :: NO_GEOPHYS
    integer(i4), intent(in) :: N_LVLS
    integer(i4), intent(in) :: NO_CONV_HTS

    type(l2pc_header_one), intent(in) :: HEADER1
    type(limb_press), intent(in) :: PTG_PRESS
    type(geom_param), intent(in) :: GEOMETRIC(no_geom)
    type(geophys_param), intent(in) :: GEOPHYSIC(no_geophys)
    type(atmos_comp), intent(in) :: ATMOSPHERIC(*)
    integer(i4), intent(in) :: NO_SPS_TBL
    integer(i4), intent(in) :: SPS_TBL(no_sps_tbl)
    real(r4), intent(in) :: CONV_PRESS(no_conv_hts)
    real(r4), intent(in) :: CONV_HTS(no_conv_hts)
    real(r4), intent(in) :: CONV_TEMP(*) ! not used
    real(r4), intent(in) :: PTG_HTS(*)
    real(r4), intent(in) :: PTG_ANGLES(no_conv_hts)
    real(r4), intent(in) :: DX_DT(Nlvl,*)
    real(r4), intent(in) :: D2X_DXDT(Nlvl,*)
    integer(i4), intent(in) :: BAND
    real(r4), intent(in) :: CENTER_ANGLE
    integer(i4), intent(in) :: FFT_PTS
    real(r4), intent(out) :: FFT_PRESS(2**fft_pts)
    real(r4), intent(in) :: I_RAW(no_conv_hts)
    character(len=*), intent(in) :: IN_DIR, Aaap
    real(r4), intent(in) :: K_STAR_GEOPHYS(Nlvl,mxco,mnp,no_geophys)
    real(r4), intent(in) :: K_STAR_ATMOS(Nlvl,mxco,mnp,*)
    real(r4), intent(in) :: A_GRID(N_lvls)
    real(r4), intent(in) :: T_GRID(N_lvls)
    real(r4), intent(in) :: Z_GRID(N_lvls)
    real(r4), intent(in) :: C_YAW, S_YAW
    real(r4), intent(in) :: C_ROLL, S_ROLL
    real(r4), intent(in) :: C_PITCH, S_PITCH
    real(r4), intent(in) :: ELEV_183, ELEV_205
    real(r4), intent(in) :: AZIM_183, AZIM_205, AZIM_REF
    real(r4), intent(in) :: GEOCSRAD
    real(r4), intent(out) :: I_STAR_ALL(ptg_press%no_lin_values)
    real(r4), intent(out) :: K_STAR_ALL(msvd,mnp,ptg_press%no_lin_values)
    integer(i4), intent(out) :: FFT_INDEX(*)
    integer(i4), intent(in) :: NO_PHI_G(no_geophys)
    integer(i4), intent(in) :: NO_PHI_F(no_sps_tbl)
!   integer(i4), intent(inout) :: IER    ! For Debugging
    integer(i4), intent(out) :: IER
!
! -----     Local Variables     ----------------------------------------
!
    real(r4) :: A, A2B2, A63B2, A_63, AZIM_ANGLE, B, B2, C_D, C_D_ELEV
    real(r4) :: C_E
    integer(i4) :: COMP_NDX
    real(r4) :: DE_DX, DE_DX_63, ELEV_OFFST, E_Z
    real(r4) :: FFT_ANGLES(max(2**fft_pts,no_conv_hts))
    integer(i4) :: GEOM_I, GEOPHYS_I
!   integer(i4) :: ICH                   ! FOR DEBUGGING
    integer(i4) :: IS, J, K, M, NF, NTR, PTG_I
    real(r4) :: Q, R, RAD(max(2**fft_pts,no_conv_hts)), ROOT
    real(r4) :: SC1(Nlvl), SC2(Nlvl), SC3(Nlvl), SC4(Nlvl)
    real(r4) :: S_D, S_D_ELEV, S_E
    integer(i4) :: SPS_I, SV_ELMNT, SV_I
    real(r4) :: TERM(Nlvl)
    real(r4) :: TX(1), TY(1)
!
! -----     Statement functions     ------------------------------------
!
! These statement functions compensate for the use of ACOSD, COSD and SIND,
! which are intrinsic in some compilers, but are not part of any standard.
!
    real(r4) :: ACOSD, COSD, SIND, X
    acosd(x) = rad2deg * acos(x)
    cosd(x) = cos(x*deg2rad)
    sind(x) = sin(x*deg2rad)
!
! -----     Begin the code     -----------------------------------------
!
! Compute the ratio of the strengths
!
! This subroutine is called by channel
!
!   ich = Ier                       ! DEBUG
    Ier = 0                         ! DEBUG
!
    Rad(1:no_conv_hts) = i_raw(1:no_conv_hts)
!
! Compute the convolution of the mixed radiances
!
    fft_angles(1:no_conv_hts) = ptg_angles
    Call fov_convolve(fft_angles,Rad,center_angle,                     &
   &                  1,no_conv_hts,band,fft_pts,In_Dir,Aaap,Ier)
    if (Ier /= 0) Return
!
! Determine radiometer 1=63, 2=205, 3=183
!
    if (band  >=  5) then                  ! 183
      azim_angle = azim_ref + azim_183
      elev_offst = elev_183
    else if (band  >=  2) then             ! 205
      azim_angle = azim_ref + azim_205
      elev_offst = elev_205
    else                                   ! 63  / pointing reference
      azim_angle = azim_ref
      elev_offst = 0.0
    end if
!
! At this stage we have I vs X, now we want I vs X^63 our common radiometer
! grid.
!
    e_z = s_pitch * c_yaw
    b = c_pitch * c_roll - s_pitch * s_yaw * s_roll
    q = s_pitch * s_yaw * c_roll + c_pitch * s_roll
    a = -cosd(azim_angle) * e_z + sind(azim_angle) * q
!
    b2 = b * b
    a2b2 = a * a + b2
!
! These are the coresponding reference values
!
    a_63 = -cosd(azim_ref) * e_z + sind(azim_ref) * q
    a63b2 = a_63 * a_63 + b2
!
    ntr = 2**fft_pts
    do ptg_i = 1, ntr
!
! Compute equivalent epsilon for the band
!
      e_z = Cos(fft_angles(ptg_i))
!
      q = a2b2 - e_z * e_z
      if (q < 0.0) then
        Ier = 1
        Print *,'** Error in subroutine: i_and_k_star'
        Print *,'   Taking Sqrt(Arg) for Arg < 0'
        Return
      end if
!
      root = sqrt(q)
      c_e = (e_z * a + b * root) / a2b2
      s_e = (e_z * b - a * root) / a2b2
!
! Adjust the fft_angles accordingly
!
      c_d_elev = cosd(elev_offst)
      s_d_elev = sind(elev_offst)
      q = c_e * (a_63 * c_d_elev - b * s_d_elev) +                     &
   &      s_e * (a_63 * s_d_elev + b * c_d_elev)
      if (abs(q) > 1.0) then
        Ier = 1
        Print *,'** Error in subroutine: i_and_k_star'
        Print *,'   Taking ACos(Arg) for abs(Arg) > 1'
        Return
      end if
!
      fft_angles(ptg_i) = ACos(q)
!
      end do
!
!  Get 'ntr' pressures associated with the fft_angles:
!
    Call get_pressures('a',a_grid,t_grid,z_grid,N_lvls,fft_angles,     &
   &                   fft_press,ntr,Ier)
    if (Ier /= 0) Return
!
! Make sure the fft_press array is MONOTONICALY increasing:
!
    is = 1
    do while (is < ntr.and.fft_press(is) >= fft_press(is+1))
      is = is + 1
    end do
!
    k = 1
    Rad(k) = Rad(is)
    fft_index(k) = is
    fft_press(k) = fft_press(is)
!
    do ptg_i = is+1, ntr
      q = fft_press(ptg_i)
      if (q > fft_press(k)) then
        k = k + 1
        fft_press(k) = q
        Rad(k) = Rad(ptg_i)
        fft_index(k) = ptg_i
      end if
    end do
!
    if (k == ntr) fft_index(1) = -2
!
! Interpolate the output values and store the radiances in: i_star_all
!
    j = ptg_press%no_lin_values
    Call Cspline_der(fft_press,ptg_press%lin_val,Rad,i_star_all,       &
   &                  term,k,j)
!
! Find out if user wants pointing derivatives
!
    if (ptg_press%der_calc(band)) then
!
! Derivatives wanted,find index location k_star_all and write the derivative
!
      is = 1
      sv_i = max_no_sv_components
      do while (ptg_press%name  /=  header1%sv_components(is) .and.    &
   &           is  <  sv_i)
        is = is + 1
      end do
!
      sv_elmnt = header1%sv_component_first_elmnt_index(is)
!
      do nf = 1, mnp
        k_star_all(sv_elmnt,nf,1:j) = term(1:j)
      end do
!
    end if
!
! Now transfer the other fwd_mdl derivatives to the output pointing
! values
!
! ****************** atmospheric derivatives ******************
!
    do sps_i = 1, no_sps_tbl
!
! no_sps_tbl is the number of species included for this channel
! check to determine if derivative is desired for that species
!
      is = sps_tbl(sps_i)
      if (atmospheric(is)%der_calc(band)) then
!
! Derivatives needed continue to process
! Find index that matches name in header
!
        comp_ndx = 1
        sv_i = max_no_sv_components
        do while (atmospheric(is)%name  /=                             &
   &      header1%sv_components(comp_ndx).and.comp_ndx < sv_i)
          comp_ndx = comp_ndx + 1
        end do
!
        m = header1%sv_component_first_elmnt_index(comp_ndx)
        sv_elmnt = m + atmospheric(is)%no_lin_values - 1
        if (sv_elmnt > msvd) then
          Ier = 1
          print 900,' variable: sv_elmnt',msvd,sv_elmnt,msvd
          Return
        end if
!
        do nf = 1, no_phi_f(sps_i)
!
          do sv_i = 1, atmospheric(is)%no_lin_values
!
! run through representation basis coefficients
!
            sv_elmnt = m + sv_i - 1
!
            Rad(1:no_conv_hts) = k_star_atmos(1:no_conv_hts,sv_i,nf,is)
!
! Now Convolve the derivative
!
            fft_angles(1:no_conv_hts) = ptg_angles
            Call fov_convolve(fft_angles,Rad,center_angle,             &
   &             1,no_conv_hts,band,fft_pts,In_Dir,Aaap,Ier)
            if (Ier /= 0) Return
!
            if (fft_index(1) > 0) then
              do ptg_i = 1, k
                Rad(ptg_i) = Rad(fft_index(ptg_i))
              end do
            end if
!
! Interpolate onto the output grid, and store in k_star_all ..
!
            Call Lintrp(fft_press,ptg_press%lin_val,Rad, &
              k_star_all(sv_elmnt,nf,1:j),k,j)
!
          end do
!
        end do
!
! *** DEBUG
!
!         if (ich == 61.and.atmospheric(is)%name == 'H2O') then
!           sv_i = 9
!           r = -1.666667
!           sv_elmnt = m + sv_i - 1
!           write(*,911) sv_i,r
!           do nf = 1, no_phi_f(sps_i)
!             write(*,918) nf,char(92),j
!             write(*,908) (k_star_all(sv_elmnt,nf,ptg_i),ptg_i=1,j)
!           end do
!         end if
!
!908  Format(6(1x,1pe12.5))
!918  Format('di_dh2o_phi_',i1,a1,i2,'n')
!
!911  Format(/,21x,'CONVOLVED k_star_atmos_dump',/,4x,
!    &'Derivatrives of Radiance with respect to H2O',i2.2,
!    &' (Zeta=',f7.4,')',/)
!
! *** END DEBUG
!

      end if
!
    end do
!
! ********************* geophysical derivatives ******************
!
    do geophys_i = 1, no_geophys
!
! check to determine if derivative is desired for this parameter
!
      if (geophysic(geophys_i)%der_calc(band)) then
!
! Derivatives needed continue to process
! Find index that matches name in header
!
        comp_ndx = 1
        sv_i = max_no_sv_components
        do while (geophysic(geophys_i)%name  /=                        &
   &      header1%sv_components(comp_ndx).and.comp_ndx < sv_i)
          comp_ndx = comp_ndx + 1
        end do
!
        m = header1%sv_component_first_elmnt_index(comp_ndx)
!
! Check for temperature, as this requires special treatment
!
        if (geophysic(geophys_i)%name == 'TEMP') then
!
          sv_elmnt = m + geophysic(geophys_i)%no_lin_values - 1
          if (sv_elmnt > msvd) then
            Ier = 1
            print 900,' variable: sv_elmnt',msvd,sv_elmnt,msvd
            Return
          end if
!
          do nf = 1, no_phi_g(geophys_i)
!
            do sv_i = 1, geophysic(geophys_i)%no_lin_values
!
! run through representation basis coefficients
!
              sv_elmnt = m + sv_i - 1
!
! Integrand over temperature derivative plus pointing differential
!
              do ptg_i = 1, no_conv_hts
                q = d2x_dxdt(ptg_i,sv_i)
                Rad(ptg_i) = i_raw(ptg_i) * q +                        &
   &                           k_star_geophys(ptg_i,sv_i,nf,geophys_i)
              end do
!
! Now, Convolve:
!
              fft_angles(1:no_conv_hts) = ptg_angles
              Call fov_convolve(fft_angles,Rad,center_angle,           &
   &               1,no_conv_hts,band,fft_pts,In_Dir,Aaap,Ier)
              if (Ier /= 0) Return
!
              if (fft_index(1) > 0) then
                do ptg_i = 1, k
                  Rad(ptg_i) = Rad(fft_index(ptg_i))
                end do
              end if
!
              Call Cspline(fft_press,ptg_press%lin_val,Rad, &
                k_star_all(sv_elmnt,nf,1:j),k,j)
!
! Now the convolution of radiance with the derivative antenna field
!
              Rad(1:no_conv_hts) = i_raw(1:no_conv_hts)
!
! Now, Convolve:
!
              fft_angles(1:no_conv_hts) = ptg_angles
              Call fov_convolve(fft_angles,Rad,center_angle,           &
   &               2,no_conv_hts,band,fft_pts,In_Dir,Aaap,Ier)
              if (Ier /= 0) Return
!
              if (fft_index(1) > 0) then
                do ptg_i = 1, k
                  Rad(ptg_i) = Rad(fft_index(ptg_i))
                end do
              end if
!
              Call Cspline(fft_press,ptg_press%lin_val,Rad,term,k,j)
!
! Transfer dx_dt from convolution grid onto the output grid
!
              Call Lintrp(conv_press,ptg_press%lin_val,                &
   &                      dx_dt(1:,sv_i),sc3,no_conv_hts,j)
!
              do ptg_i = 1, j
                r = sc3(ptg_i) * term(ptg_i)
                q = k_star_all(sv_elmnt,nf,ptg_i)
                k_star_all(sv_elmnt,nf,ptg_i) = q + r
              end do
!
! the convolution of the radiance weighted hydrostatic derivative
! with the antenna derivative
!
              Rad(1:no_conv_hts) = &
                dx_dt(1:no_conv_hts,sv_i) * i_raw(1:no_conv_hts)
!
! Now, convolve:
!
              fft_angles(1:no_conv_hts) = ptg_angles
              Call fov_convolve(fft_angles,Rad,center_angle,           &
   &               2,no_conv_hts,band,fft_pts,In_Dir,Aaap,Ier)
              if (Ier /= 0) Return
!
              if (fft_index(1) > 0) then
                do ptg_i = 1, k
                  Rad(ptg_i) = Rad(fft_index(ptg_i))
                end do
              end if
!
              Call Cspline(fft_press,ptg_press%lin_val,Rad,term,k,j)
!
              do ptg_i = 1, j
                q = k_star_all(sv_elmnt,nf,ptg_i)
                k_star_all(sv_elmnt,nf,ptg_i) = q - term(ptg_i)
              end do
!
            end do
!
          end do
!
! *** DEBUG
!
!           if (ich == 61.and.geophysic(geophys_i)%name == 'TEMP') then
!             sv_i = 9
!             sv_elmnt = m + sv_i - 1
!             r = -1.666667
!             write(*,910) sv_i,r
!             do nf = 1, no_phi_g(geophys_i)
!               write(*,919) nf,j
!               write(*,909) (k_star_all(sv_elmnt,nf,ptg_i),ptg_i=1,j)
!             end do
!           end if
!
!909  Format(6(1x,1pe12.5))
!919  Format('di_dtemp_phi_',i1,'\\',i2,'n')
!
!910  Format(/,21x,'CONVOLVED k_star_geophys dump',/,4x,
!    &'Derivatrives of Radiance with respect to TEMP',i2.2,
!    &' (Zeta=',f7.4,')',/)
!
! *** END DEBUG
!
        else
!
! Use standard processing for everything else
!
          sv_elmnt = m + geophysic(geophys_i)%no_lin_values - 1
          if (sv_elmnt > msvd) then
            Ier = 1
            print 900,' variable: sv_elmnt',msvd,sv_elmnt,msvd
            Return
          end if
!
          do nf = 1, no_phi_g(geophys_i)
!
            do sv_i = 1, geophysic(geophys_i)%no_lin_values
!
! run through representation basis coefficients
!
              sv_elmnt = m + sv_i - 1
!
              Rad(1:no_conv_hts) = &
                k_star_geophys(1:no_conv_hts,sv_i,nf,geophys_i)
!
! Now convolve:
!
              fft_angles(1:no_conv_hts) = ptg_angles
              Call fov_convolve(fft_angles,Rad,center_angle,           &
   &               1,no_conv_hts,band,fft_pts,In_Dir,Aaap,Ier)
              if (Ier /= 0) Return
!
              if (fft_index(1) > 0) then
                do ptg_i = 1, k
                  Rad(ptg_i) = Rad(fft_index(ptg_i))
                end do
              end if
!
              Call Cspline(fft_press,ptg_press%lin_val,Rad, &
                k_star_all(sv_elmnt,nf,1:j),k,j)
!
            end do
!
          end do
!
        end if
!
      end if
!
    end do
!
!********************** geometric derivatives ***************
! Need to consider fov effects due to observer to tangent point
! stretching impact. This affects earth_radius and observer position.
!
    do geom_i = 1, no_geom
!
! check to determine if derivative is desired for this parameter
!
      if (geometric(geom_i)%der_calc(band)) then
!
! Derivatives needed continue to process
! Find index that matches name in header
!
        comp_ndx = 1
        sv_i = max_no_sv_components
        do while (geometric(geom_i)%name  /=                           &
   &      header1%sv_components(comp_ndx).and.comp_ndx < sv_i)
          comp_ndx = comp_ndx + 1
        end do
!
        sv_elmnt = header1%sv_component_first_elmnt_index(comp_ndx)
        if (sv_elmnt > msvd) then
          Ier = 1
          print 900,' variable: sv_elmnt',msvd,sv_elmnt,msvd
          Return
        end if
!
! Search for specific element
!
        if (geometric(geom_i)%name == 'GEOCERAD') then
!
! compute the derivative per pointing
!
          Rad(1:no_conv_hts) = i_raw(1:no_conv_hts) / &
            (Cos(ptg_angles(1:no_conv_hts))**2 * conv_hts(1:no_conv_hts))
!
! Now convolve:
!
          fft_angles(1:no_conv_hts) = ptg_angles
          Call fov_convolve(fft_angles,Rad,center_angle,               &
   &           1,no_conv_hts,band,fft_pts,In_Dir,Aaap,Ier)
          if (Ier /= 0) Return
!
          if (fft_index(1) > 0) then
            do ptg_i = 1, k
              Rad(ptg_i) = Rad(fft_index(ptg_i))
            end do
          end if
!
          Call Cspline(fft_press,ptg_press%lin_val,Rad,term,k,j)
!
          do nf = 1, mnp
            k_star_all(sv_elmnt,nf,1:j) = term(1:j)
          end do
!
! Now do differential antenna pattern contribution first term
!
          Rad(1:no_conv_hts) = i_raw(1:no_conv_hts)
!
! Now convolve:
!
          fft_angles(1:no_conv_hts) = ptg_angles
          Call fov_convolve(fft_angles,Rad,center_angle,               &
   &           2,no_conv_hts,band,fft_pts,In_Dir,Aaap,Ier)
          if (Ier /= 0) Return
!
! Interpolate to output grid
!
          Call Lintrp(conv_press,ptg_press%lin_val,ptg_angles,         &
   &                  fft_angles,no_conv_hts,j)
!
          if (fft_index(1) > 0) then
            do ptg_i = 1, k
              Rad(ptg_i) = Rad(fft_index(ptg_i))
            end do
          end if
!
          Call Cspline(fft_press,ptg_press%lin_val,Rad,term,k,j)
!
          do nf = 1, mnp
            do ptg_i = 1, j
              c_e = tan(fft_angles(ptg_i)) / ptg_hts(ptg_i)
              q = k_star_all(sv_elmnt,nf,ptg_i)
              k_star_all(sv_elmnt,nf,ptg_i) = q + c_e * term(ptg_i)
            end do
          end do
!
! Now do differential antenna pattern contribution 2nd term
!
          Rad(1:no_conv_hts) = i_raw(1:no_conv_hts) * &
            (tan(ptg_angles(1:no_conv_hts)) / conv_hts(1:no_conv_hts))
!
! Now convolve:
!
          fft_angles(1:no_conv_hts) = ptg_angles
          Call fov_convolve(fft_angles,Rad,center_angle,               &
   &           2,no_conv_hts,band,fft_pts,In_Dir,Aaap,Ier)
          if (Ier /= 0) Return
!
! Interpolate to output grid
!
          if (fft_index(1) > 0) then
            do ptg_i = 1, k
              Rad(ptg_i) = Rad(fft_index(ptg_i))
            end do
          end if
!
          Call Cspline(fft_press,ptg_press%lin_val,Rad,term,k,j)
!
          do nf = 1, mnp
            do ptg_i = 1, j
              q = k_star_all(sv_elmnt,nf,ptg_i)
              k_star_all(sv_elmnt,nf,ptg_i) = q - term(ptg_i)
            end do
          end do
!
        else if (geometric(geom_i)%name == 'GEOCSRAD') then
!
! compute the derivative per pointing
!
          Rad(1:no_conv_hts) = - i_raw(1:no_conv_hts) / &
            (Cos(ptg_angles(1:no_conv_hts))**2 * geocsrad)
!
! Now convolve:
!
          fft_angles(1:no_conv_hts) = ptg_angles
          Call fov_convolve(fft_angles,Rad,center_angle,               &
   &           1,no_conv_hts,band,fft_pts,In_Dir,Aaap,Ier)
          if (Ier /= 0) Return
!
          if (fft_index(1) > 0) then
            do ptg_i = 1, k
              Rad(ptg_i) = Rad(fft_index(ptg_i))
            end do
          end if
!
          Call Cspline(fft_press,ptg_press%lin_val,Rad,term,k,j)
!
          do nf = 1, mnp
            k_star_all(sv_elmnt,nf,1:j) = term(1:j)
          end do
!
! Now do differential antenna pattern contribution first term
!
! average the sideband values
!
          Rad(1:no_conv_hts) = i_raw(1:no_conv_hts)
!
! Now convolve:
!
          fft_angles(1:no_conv_hts) = ptg_angles
          Call fov_convolve(fft_angles,Rad,center_angle,               &
   &           2,no_conv_hts,band,fft_pts,In_Dir,Aaap,Ier)
          if (Ier /= 0) Return
!
! Interpolate to output grid
!
          Call Lintrp(conv_press,ptg_press%lin_val,ptg_angles,         &
   &                  fft_angles,no_conv_hts,j)
!
          if (fft_index(1) > 0) then
            do ptg_i = 1, k
              Rad(ptg_i) = Rad(fft_index(ptg_i))
            end do
          end if
!
          Call Cspline(fft_press,ptg_press%lin_val,Rad,term,k,j)
!
          do nf = 1, mnp
            do ptg_i = 1, j
              c_e = tan(fft_angles(ptg_i)) / geocsrad
              q = k_star_all(sv_elmnt,nf,ptg_i)
              k_star_all(sv_elmnt,nf,ptg_i) = q -  c_e * term(ptg_i)
            end do
          end do
!
! Now do differential antenna pattern contribution 2nd term
!
          Rad(1:no_conv_hts) = -i_raw(1:no_conv_hts) * &
            (tan(ptg_angles(1:no_conv_hts)) / geocsrad)
!
! Now convolve:
!
          fft_angles(1:no_conv_hts) = ptg_angles
          Call fov_convolve(fft_angles,Rad,center_angle,               &
   &           2,no_conv_hts,band,fft_pts,In_Dir,Aaap,Ier)
          if (Ier /= 0) Return
!
! Interpolate to output grid
!
          if (fft_index(1) > 0) then
            do ptg_i = 1, k
              Rad(ptg_i) = Rad(fft_index(ptg_i))
            end do
          end if
!
          Call Cspline(fft_press,ptg_press%lin_val,Rad,term,k,j)
!
          do nf = 1, mnp
            do ptg_i = 1, j
              q = k_star_all(sv_elmnt,nf,ptg_i)
              k_star_all(sv_elmnt,nf,ptg_i) = q - term(ptg_i)
            end do
          end do
!
        end if
!
      end if
!
    end do
!
! Establish dI/dX, X, X^63 data
! Get X^63
!
    Call Lintrp(conv_press,ptg_press%lin_val,ptg_angles,fft_angles,    &
   &            no_conv_hts,j)
!
! Get dI/dX at 63
!
    ptg_i = 1
!   s_e = 0.0
!   c_e = fft_angles(1)
!   Call csinterp(fft_angles,c_e,i_star_all,s_e,j,ptg_i,j,ptg_i,       &
!  &              sc1,sc2,sc3,sc4)
    tx = 0.0
    ty = fft_angles(1)
    ier = csinterp(fft_angles,tx,i_star_all,ty,j,ptg_i,j,ptg_i,       &
   &              sc1,sc2,sc3,sc4)
    if ( ier /= 0 ) then
      print *, '-E- In I_AND_K_STAR, CSINTERP failed'
      return
    end if
!
! Where sc2 is dI/dX at 63
! Now we get derivative premultipliers -1/sqrt(1.0 - e_z^2)
! place absolute in sc3 and 63 in sc4
!
    do ptg_i = 1, j
!
! Compute equivalent epsilon 63 for the band
!
      e_z = Cos(fft_angles(ptg_i))
      root = sqrt(a63b2 - e_z*e_z)
      c_e = (e_z * a_63 + b * root) / a63b2
      s_e = (e_z * b - a_63 * root) / a63b2
!
! Compute derivative premultipliers
!
      c_d_elev = cosd(elev_offst)
      s_d_elev = sind(elev_offst)
      e_z = c_e * (a * c_d_elev + b * s_d_elev)                        &
   &      - s_e * (a * s_d_elev - b * c_d_elev)
      sc3(ptg_i) = -1.0 / sqrt(1.0 - e_z * e_z)
      sc4(ptg_i) = -1.0 / sin(fft_angles(ptg_i))
      sc1(ptg_i) = Acosd(c_e)          ! elevation for 63 GHz
!
    end do
!
    do geom_i = 1, no_geom
!
! check to determine if derivative is desired for this parameter
!
      if (geometric(geom_i)%der_calc(band)) then
!
! Derivatives needed continue to process
!
! Find index that matches name in header
!
        comp_ndx = 1
        do while (geometric(geom_i)%name  /=                           &
   &      header1%sv_components(comp_ndx))
          comp_ndx = comp_ndx + 1
        end do
!
        sv_elmnt = header1%sv_component_first_elmnt_index(comp_ndx)
        if (sv_elmnt > msvd) then
          Ier = 1
          print 900,' variable: sv_elmnt',msvd,sv_elmnt,msvd
          Return
        end if
!
! Search for specific element
!
        if (geometric(geom_i)%name == 'GEOCERAD') then
!
! compute the derivative per pointing
!
          do nf = 1, mnp
            do ptg_i = 1, j
              q = sc3(ptg_i) * sc3(ptg_i) - 1.0
              if (q > 0.0) then
                c_d = k_star_all(sv_elmnt,nf,ptg_i)
                e_z = sc4(ptg_i) * sc4(ptg_i) - 1.0
                k_star_all(sv_elmnt,nf,ptg_i) = c_d + sc2(ptg_i) *     &
   &              (sc3(ptg_i)/sqrt(q) - sc4(ptg_i)/sqrt(e_z))/geocsrad
              end if
            end do
          end do
!
        else if (geometric(geom_i)%name == 'GEOCSRAD') then
!
          do nf = 1, mnp
            do ptg_i = 1, j
              q = sc3(ptg_i) * sc3(ptg_i) - 1.0
              if (q > 0.0) then
                c_d = k_star_all(sv_elmnt,nf,ptg_i)
                e_z = sc4(ptg_i) * sc4(ptg_i) - 1.0
                k_star_all(sv_elmnt,nf,ptg_i) = c_d - sc2(ptg_i) *     &
   &               (1.0/sqrt(q) - 1.0/sqrt(e_z)) / geocsrad
              end if
            end do
          end do
!
        else if ( geometric(geom_i)%name == 'ELEV_183'                 &
   &        .or. geometric(geom_i)%name == 'ELEV_205' ) then
!
          do nf = 1, mnp
            do ptg_i = 1, j
              s_e = sc1(ptg_i) + elev_offst
              de_dx = b * cosd(s_e) - a * sind(s_e)
              k_star_all(sv_elmnt,nf,ptg_i) =                          &
   &                    deg2rad * sc2(ptg_i) * sc3(ptg_i) * de_dx
            end do
          end do
!
        else if (geometric(geom_i)%name == 'AZIM_REF' ) then
!
          e_z = s_pitch * c_yaw
          c_d = s_pitch * s_yaw * c_roll + c_pitch * s_roll
!
          do nf = 1, mnp
            do ptg_i = 1, j
              s_e = sc1(ptg_i) + elev_offst
              de_dx = cosd(azim_angle) * cosd(s_e) * c_d               &
   &                + sind(azim_angle) * cosd(s_e) * e_z
              de_dx_63 = cosd(azim_ref) * cosd(sc1(ptg_i)) * c_d       &
   &                   + sind(azim_ref) * cosd(sc1(ptg_i)) * e_z
              k_star_all(sv_elmnt,nf,ptg_i) = sc2(ptg_i) * deg2rad     &
   &                * (sc3(ptg_i) * de_dx - sc4(ptg_i) * de_dx_63)
            end do
          end do
!
        else if ( geometric(geom_i)%name == 'AZIM_183'                 &
   &        .or. geometric(geom_i)%name == 'AZIM_205' ) then
!
          e_z = s_pitch * c_yaw
          c_d = s_pitch * s_yaw * c_roll + c_pitch * s_roll
!
          do nf = 1, mnp
            do ptg_i = 1, j
              s_e = sc1(ptg_i) + elev_offst
              de_dx = cosd(azim_angle) * cosd(s_e) * c_d               &
   &                + sind(azim_angle) * cosd(s_e) * e_z
              k_star_all(sv_elmnt,nf,ptg_i) =                          &
   &                    deg2rad * sc2(ptg_i) * sc3(ptg_i) * de_dx
            end do
          end do
!
        else if (geometric(geom_i)%name == 'ROLL') then
!
          c_d = c_pitch * c_roll - s_pitch * s_yaw * s_roll
          s_d = c_pitch * s_roll + s_pitch * s_yaw * c_roll
!
          do nf = 1, mnp
            do ptg_i = 1, j
              s_e = sc1(ptg_i) + elev_offst
              de_dx = sind(azim_angle) * cosd(s_e) * c_d -             &
   &                  sind(s_e) * s_d
              de_dx_63 = sind(azim_ref) * cosd(sc1(ptg_i)) * c_d -     &
   &                     sind(sc1(ptg_i)) * s_d
              k_star_all(sv_elmnt,nf,ptg_i) = sc2(ptg_i) * deg2rad     &
   &                * (sc3(ptg_i) * de_dx - sc4(ptg_i) * de_dx_63)
            end do
          end do
!
        else if (geometric(geom_i)%name == 'PITCH') then
!
          e_z = c_pitch * c_yaw
          c_d = c_pitch * s_yaw * c_roll - s_pitch * s_roll
!
          do nf = 1, mnp
            do ptg_i = 1, j
              s_d = cosd(sc1(ptg_i))
              s_e = sc1(ptg_i) + elev_offst
              de_dx = c_d * (sind(azim_angle) * cosd(s_e) + sind(s_e)) &
   &                  - cosd(azim_angle) * cosd(s_e) * e_z
              de_dx_63 = c_d * (sind(azim_ref)*s_d + sind(sc1(ptg_i))) &
   &                     - cosd(azim_ref) * s_d * e_z
              k_star_all(sv_elmnt,nf,ptg_i) = sc2(ptg_i) * deg2rad     &
   &                * (sc3(ptg_i) * de_dx - sc4(ptg_i) * de_dx_63)
            end do
          end do
!
        else if (geometric(geom_i)%name == 'YAW') then
!
          c_d = s_pitch * s_yaw
          q = s_pitch * c_yaw
          e_z = q * s_roll
          s_d = q * c_roll
!
          do nf = 1, mnp
            do ptg_i = 1, j
              s_e = sc1(ptg_i) + elev_offst
              de_dx = sind(azim_angle) * cosd(s_e)*s_d - sind(s_e)*e_z &
   &                + cosd(azim_angle) * cosd(s_e) * c_d
              de_dx_63 = sind(azim_ref) * cosd(sc1(ptg_i)) * s_d       &
   &                   + cosd(azim_ref) * cosd(sc1(ptg_i)) * c_d       &
   &                   - sind(sc1(ptg_i)) * e_z
              k_star_all(sv_elmnt,nf,ptg_i) = sc2(ptg_i) * deg2rad     &
   &                * (sc3(ptg_i) * de_dx - sc4(ptg_i) * de_dx_63)
            end do
          end do
!
        end if
!
      end if
!
    end do
!
900 format('** Error in subroutine: i_and_k_star ..',a,/,              &
   &       '   Exceeds maximum number of derivatives:',i5,/,           &
   &       '   New value should be at least:',i4,' large',/,           &
   &       '   Need to change it in: l2pc_file_paramters.inc',/,       &
   &       '   parameter (max_no_sv_derivatives =',i4,')')
!
    Return
  End Subroutine I_AND_K_STAR
end module I_AND_K_STAR_M

! $Log$
