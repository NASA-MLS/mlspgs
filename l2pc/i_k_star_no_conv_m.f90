module I_K_STAR_NO_CONV_M
  use HYDROSTATIC_INTRP, only: GET_PRESSURES
  use L2PC_FILE_PARAMETERS, only: DEG2RAD, MAX_NO_SV_COMPONENTS, &
                                  MSVD => max_no_sv_derivatives, &
                                  MXCO => max_no_elmnts_per_sv_component, &
                                  RAD2DEG
  use L2PC_FILE_STRUCTURES, only: L2PC_HEADER_ONE
  use L2PC_PFA_STRUCTURES, only: ATMOS_COMP, GEOM_PARAM, GEOPHYS_PARAM, &
                                 LIMB_PRESS
  use L2PCDIM, only: MNP => max_no_phi, NLVL
  use MLSCommon, only: I4, R4, R8
  use D_CSINTERP_M, only: CSINTERP
  use D_CSPLINE_M, only: CSPLINE
  use D_LINTRP_M, only: LINTRP
  use DCSPLINE_DER_M, only: CSPLINE_DER
  implicit NONE
  private
  public :: I_K_STAR_NO_CONV
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName = &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
!----------------------------------------------------------------------
! This module is the same as i_and_k_star subroutine but WITHOUT
! CONVOLVING THE DATA !!
!
  Subroutine I_K_STAR_NO_CONV ( header1, ptg_press, geometric, geophysic, &
 &           atmospheric, no_geom, no_geophys, no_sps_tbl, sps_tbl,       &
 &           conv_press, conv_hts, ptg_angles, band, i_raw,               &
 &           k_star_geophys, k_star_atmos, a_grid, t_grid, z_grid,        &
 &           N_lvls, no_conv_hts, c_yaw, s_yaw, c_roll, s_roll, c_pitch,  &
 &           s_pitch, elev_183, elev_205, azim_183, azim_205, azim_ref,   &
 &           geocsrad, i_star_all, k_star_all, no_phi_g, no_phi_f, Ier )
!
    type(l2pc_header_one), intent(in) :: HEADER1
    type(limb_press), intent(in) :: PTG_PRESS
    type(geom_param), intent(in) :: GEOMETRIC(no_geom)
    type(geophys_param), intent(in) :: GEOPHYSIC(no_geophys)
    type(atmos_comp), intent(in) :: ATMOSPHERIC(*)
    integer(i4), intent(in) :: NO_GEOM
    integer(i4), intent(in) :: NO_GEOPHYS
    integer(i4), intent(in) :: NO_SPS_TBL
    integer(i4), intent(in) :: SPS_TBL(no_sps_tbl)
    real(r8), intent(in) :: CONV_PRESS(no_conv_hts)
    real(r8), intent(in) :: CONV_HTS(no_conv_hts)
    real(r8), intent(in) :: PTG_ANGLES(no_conv_hts)
    integer(i4), intent(in) :: BAND
    real(r8), intent(in) :: I_RAW(no_conv_hts)
    real(r4), intent(in) :: K_STAR_GEOPHYS(Nlvl,mxco,mnp,no_geophys)
    real(r4), intent(in) :: K_STAR_ATMOS(Nlvl,mxco,mnp,*)
    real(r8), intent(in) :: A_GRID(N_lvls)
    real(r8), intent(in) :: T_GRID(N_lvls)
    real(r8), intent(in) :: Z_GRID(N_lvls)
    integer(i4), intent(in) :: N_LVLS
    integer(i4), intent(in) :: NO_CONV_HTS
    real(r8), intent(in) :: C_YAW, S_YAW
    real(r8), intent(in) :: C_ROLL, S_ROLL
    real(r8), intent(in) :: C_PITCH, S_PITCH
    real(r8), intent(in) :: ELEV_183, ELEV_205
    real(r8), intent(in) :: AZIM_183, AZIM_205, AZIM_REF
    real(r8), intent(in) :: GEOCSRAD
    real(r8), intent(out) :: I_STAR_ALL(ptg_press%no_lin_values)
    real(r4), intent(out) :: K_STAR_ALL(msvd,mnp,ptg_press%no_lin_values)
    integer(i4), intent(in) :: NO_PHI_G(no_geophys)
    integer(i4), intent(in) :: NO_PHI_F(no_sps_tbl)
!   integer(i4), intent(inout) :: IER    ! For Debugging
    integer(i4), intent(out) :: IER
!
! -----     Local Variables     ----------------------------------------
!
    real(r8) :: A, A2B2, A63B2, A_63
    real(r8) :: ANG(Nlvl)
    real(r8) :: AZIM_ANGLE, B, B2, C_D, C_D_ELEV, C_E
    integer(i4) :: COMP_NDX
    real(r8) :: DE_DX, DE_DX_63, ELEV_OFFST, E_Z
    integer(i4) :: GEOM_I, GEOPHYS_I
!   integer(i4) :: ICH                   ! FOR DEBUGGING
    integer(i4) :: IS, J, K, KCONV, M, NF, PTG_I
    real(r8) :: PRESS(Nlvl),PtP(Nlvl)
    real(r8) :: Q
    real(r8) :: R, RAD(Nlvl)
    real(r8) :: ROOT
    real(r8) :: SC1(Nlvl), SC2(Nlvl), SC3(Nlvl), SC4(Nlvl)
    real(r8) :: S_D, S_D_ELEV, S_E
    integer(i4) :: SPS_I, SURNDX, SV_ELMNT, SV_I
    real(r8) :: TERM(Nlvl), TX(1), TY(1)
!
! =====     Private procedures     =====================================
!
! These functions compensate for the use of ACOSD, COSD and SIND,
! which are intrinsic in some compilers, but are not part of any standard.
!
  real(r8) :: DACOSD,DCOSD,DSIND,X
    dacosd(x) = rad2deg * dacos(x)
    dcosd(x)  = dcos(x*deg2rad)
    dsind(x)  = dsin(x*deg2rad)
!
  real(r8) :: ACOSD,COSD,SIND,SX
    acosd(sx) = rad2deg * acos(sx)
    cosd(sx)  = cos(sx*deg2rad)
    sind(sx)  = sin(sx*deg2rad)
!
  real(r8) :: SACOSD,SCOSD,SSIND
    sacosd(sx) = rad2deg * acos(sx)
    scosd(sx)  = cos(sx*deg2rad)
    ssind(sx)  = sin(sx*deg2rad)
!
! -----     Begin the code     -----------------------------------------
!
! This subroutine is called by channel
!
!   Ich = Ier                   ! DEBUG
    Ier = 0
!
    j = ptg_press%no_lin_values
    PtP(1:j) = dble(ptg_press%lin_val(1:j))
!
    kconv = 0
    surndx = -1
    do ptg_i = 1, no_conv_hts
      r = conv_hts(ptg_i)
      if (r > -0.01) then
        kconv = kconv + 1
        Rad(kconv) = i_raw(ptg_i)
        if (surndx < 1) surndx = ptg_i
      end if
    end do
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
! At this stage we have I vs X,now we want I vs X^63 our common radiometer
! grid.
!
    e_z = s_pitch * c_yaw
    b = c_pitch * c_roll - s_pitch * s_yaw * s_roll
    q = s_pitch * s_yaw * c_roll + c_pitch * s_roll
    a = -dcosd(azim_angle) * e_z + dsind(azim_angle) * q
!
    b2 = b * b
    a2b2 = a * a + b2
!
! These are the coresponding reference values
!
    a_63 = -cosd(azim_ref) * e_z + sind(azim_ref) * q
    a63b2 = a_63 * a_63 + b2
!
    do ptg_i = 1, kconv
!
! Compute equivalent epsilon for the band
!
      e_z = cos(ptg_angles(ptg_i+surndx-1))
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
! Adjust the 'Ang' array accordingly
!
      c_d_elev = dcosd(elev_offst)
      s_d_elev = dsind(elev_offst)
      q = c_e * (a_63 * c_d_elev - b * s_d_elev) +                     &
   &      s_e * (a_63 * s_d_elev + b * c_d_elev)
      if (abs(q) > 1.0) then
        Ier = 1
        Print *,'** Error in subroutine: i_k_star_no_conv'
        Print *,'   Taking Acos(Arg) for abs(Arg) > 1'
        Return
      end if
!
      Ang(ptg_i) = DAcos(q)
!
    end do
!
!  Get 'kconv' pressures associated with the 'Ang' array:
!
    Call get_pressures('a',a_grid,t_grid,z_grid,N_lvls,Ang,            &
   &                   Press,kconv,Ier)
    if (Ier /= 0) Return
!
! interpolate the output values and store the radiances in i_star_all
!
    Call Cspline_der(Press,PtP,Rad,i_star_all,sc3,kconv,j)
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
      if (sv_elmnt > msvd) then
        Ier = 1
        print 900,' variable: sv_elmnt',msvd,sv_elmnt,msvd
        Return
      end if
!
      do nf = 1, mnp
        k_star_all(sv_elmnt,nf,1:j) = sc3(1:j)
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
      if (atmospheric((is))%der_calc(band)) then
!
! Derivatives needed continue to process
! Find index that matches name in header
!
        comp_ndx = 1
        sv_i = max_no_sv_components
        do while (atmospheric((is))%name  /=                             &
   &      header1%sv_components(comp_ndx).and.comp_ndx < sv_i)
          comp_ndx = comp_ndx + 1
        end do
!
        m = header1%sv_component_first_elmnt_index(comp_ndx)
        sv_i = m + atmospheric((is))%no_lin_values - 1
        if (sv_i > msvd) then
          Ier = 1
          print 900,' variable: sv_elmnt',msvd,sv_i,msvd
          Return
        end if
!
        do nf = 1, no_phi_f(sps_i)
!
          do sv_i = 1, atmospheric((is))%no_lin_values
!
! run through representation basis coefficients
!
            sv_elmnt = m + sv_i - 1
!
            ptg_i = surndx+kconv-1
            Rad(1:kconv) = k_star_atmos(surndx:ptg_i,sv_i,nf,is)
!
! Interpolate onto the output grid, and store in k_star_all ..
!
            Call Lintrp(conv_press(surndx:),PtP,Rad,Sc1,kconv,j)
            k_star_all(sv_elmnt,nf,1:j) = Sc1(1:j)
!
          end do
!
        end do
!
! *** DEBUG
!
!         if (ich == 61.and.atmospheric((is))%name == 'H2O') then
!           sv_i = 9
!           r = -1.666667
!           sv_elmnt = m + sv_i - 1
!           write(*,911) sv_i,r
!           do nf = 1, no_phi_f(sps_i)
!             write(*,918) nf,char(92),j
!             write(*,908) (k_star_all(sv_elmnt,nf,l),l=1,j)
!           end do
!         end if
!
!908  Format(6(1x,1pe12.5))
!918  Format('di_dh2o_phi_',i1,a1,i2,'n')
!
!911  Format(/,21x,'UN-CONVOLVED k_star_atmos_dump',/,4x,
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
      if (geophysic((geophys_i))%der_calc(band)) then
!
! Derivatives needed continue to process:
!
! Find index that matches name in header
!
        comp_ndx = 1
        sv_i = max_no_sv_components
        do while (geophysic((geophys_i))%name  /=                        &
   &      header1%sv_components(comp_ndx).and.comp_ndx < sv_i)
          comp_ndx = comp_ndx + 1
        end do
!
        m = header1%sv_component_first_elmnt_index(comp_ndx)
        sv_i = m + geophysic((geophys_i))%no_lin_values - 1
        if (sv_i > msvd) then
          Ier = 1
          print 900,' variable: sv_elmnt',msvd,sv_i,msvd
          Return
        end if
!
        do nf = 1, no_phi_g(geophys_i)
!
          do sv_i = 1, geophysic((geophys_i))%no_lin_values
!
! run through representation basis coefficients
!
            sv_elmnt = m + sv_i - 1
!
            ptg_i = kconv-surndx+1
            Rad(1:kconv)=k_star_geophys(surndx:ptg_i,sv_i,nf,geophys_i)
!
            Call Cspline(conv_press(surndx:),PtP,Rad,Sc1,kconv,j)
            k_star_all(sv_elmnt,nf,1:j) = Sc1(1:j)
!
          end do
!
        end do
!
! *** DEBUG
!
!         if (ich == 61.and.geophysic((geophys_i))%name == 'TEMP') then
!           sv_i = 9
!           sv_elmnt = m + sv_i - 1
!           r = -1.666667
!           write(*,910) sv_i,r
!           do nf = 1, no_phi_g(geophys_i)
!             write(*,919) nf,j
!             write(*,909) (k_star_all(sv_elmnt,nf,l),l=1,j)
!           end do
!         end if
!
!909  Format(6(1x,1pe12.5))
!919  Format('di_dtemp_phi_',i1,'\\',i2,'n')
!
!910  Format(/,21x,'UN-CONVOLVED k_star_geophys dump',/,4x,
!    &'Derivatrives of Radiance with respect to TEMP',i2.2,
!    &' (Zeta=',f7.4,')',/)
!
! *** END DEBUG
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
      if (geometric((geom_i))%der_calc(band)) then
!
! Derivatives needed continue to process
! Find index that matches name in header
!
        comp_ndx = 1
        sv_i = max_no_sv_components
        do while (geometric((geom_i))%name  /=                           &
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
        if (geometric((geom_i))%name  ==  'GEOCERAD') then
!
! We set the derivative for NO CONVOLUTION CASE to zero in this application:
!
          k_star_all(sv_elmnt,1:mnp,1:j) = 0.0
!
        else if (geometric((geom_i))%name  ==  'GEOCSRAD') then
!
! We set the derivative for NO CONVOLUTION CASE to zero in this application:
!
          k_star_all(sv_elmnt,1:mnp,1:j) = 0.0
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
    Call Lintrp(conv_press(surndx:),PtP,ptg_angles(surndx:),Ang,kconv,j)
!
! Get dI/dX at 63
!
    m = 1
    ty = 0.0
    tx = Ang(1)
    ier = csinterp(Ang,tx,i_star_all,ty,j,m,j,m,sc1,sc2,sc3,sc4)
    if ( ier /= 0 ) then
      print *, '-E- In I_K_STAR_NO_CONV, CSINTERP failed'
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
      e_z = cos(Ang(ptg_i))
      root = sqrt(a63b2 - e_z*e_z)
      c_e = (e_z * a_63 + b * root) / a63b2
      s_e = (e_z * b - a_63 * root) / a63b2
!
! Compute derivative premultipliers
!
      c_d_elev = dcosd(elev_offst)
      s_d_elev = dsind(elev_offst)
      e_z = c_e * (a * c_d_elev + b * s_d_elev)                        &
   &      - s_e * (a * s_d_elev - b * c_d_elev)
      q = -1.0 / sqrt(1.0 - e_z * e_z)
      sc3(ptg_i) = q
      sc4(ptg_i) = -1.0 / sin(Ang(ptg_i))
      sc1(ptg_i) = DAcosd(c_e)          ! elevation for 63 GHz
!
    end do
!
    do geom_i = 1, no_geom
!
! check to determine if derivative is desired for this parameter
!
      if (geometric((geom_i))%der_calc(band)) then
!
! Derivatives needed continue to process
!
! Find index that matches name in header
!
        comp_ndx = 1
        do while (geometric((geom_i))%name  /=                           &
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
        if ( geometric((geom_i))%name  ==  'ELEV_183'                    &
   &        .or. geometric((geom_i))%name  ==  'ELEV_205' ) then
!
          do nf = 1, mnp
            do ptg_i = 1, j
              s_e = sc1(ptg_i) + elev_offst
              de_dx = b * dcosd(s_e) - a * dsind(s_e)
              k_star_all(sv_elmnt,nf,ptg_i) =                          &
   &                    deg2rad * sc2(ptg_i) * sc3(ptg_i) * de_dx
            end do
          end do
!
        else if (geometric((geom_i))%name  ==  'AZIM_REF' ) then
!
          e_z = s_pitch * c_yaw
          c_d = s_pitch * s_yaw * c_roll + c_pitch * s_roll
          do nf = 1, mnp
            do ptg_i = 1, j
              s_e = sc1(ptg_i) + elev_offst
              de_dx = dcosd(azim_angle) * dcosd(s_e) * c_d               &
   &                + dsind(azim_angle) * dcosd(s_e) * e_z
              de_dx_63 = cosd(azim_ref) * cosd(sc1(ptg_i)) * c_d       &
   &                   + sind(azim_ref) * cosd(sc1(ptg_i)) * e_z
              k_star_all(sv_elmnt,nf,ptg_i) = sc2(ptg_i) * deg2rad     &
   &                * (sc3(ptg_i) * de_dx - sc4(ptg_i) * de_dx_63)
            end do
          end do
!
        else if ( geometric((geom_i))%name  ==  'AZIM_183'               &
   &        .or. geometric((geom_i))%name  ==  'AZIM_205' ) then
!
          e_z = s_pitch * c_yaw
          c_d = s_pitch * s_yaw * c_roll + c_pitch * s_roll
          do nf = 1, mnp
           do ptg_i = 1, j
             s_e = sc1(ptg_i) + elev_offst
             de_dx = dcosd(azim_angle) * dcosd(s_e) * c_d                &
   &               + dsind(azim_angle) * dcosd(s_e) * e_z
             k_star_all(sv_elmnt,nf,ptg_i) =                           &
   &                   deg2rad * sc2(ptg_i) * sc3(ptg_i) * de_dx
            end do
          end do
!
        else if (geometric((geom_i))%name  ==  'ROLL') then
!
          c_d = c_pitch * c_roll - s_pitch * s_yaw * s_roll
          s_d = c_pitch * s_roll + s_pitch * s_yaw * c_roll
          do nf = 1, mnp
           do ptg_i = 1, j
             s_e = sc1(ptg_i) + elev_offst
             de_dx = dsind(azim_angle) * dcosd(s_e) * c_d -            &
   &                 dsind(s_e) * s_d
             de_dx_63 = sind(azim_ref) * cosd(sc1(ptg_i)) * c_d -      &
   &                    sind(sc1(ptg_i)) * s_d
              k_star_all(sv_elmnt,nf,ptg_i) = sc2(ptg_i) * deg2rad     &
   &                * (sc3(ptg_i) * de_dx - sc4(ptg_i) * de_dx_63)
            end do
          end do
!
        else if (geometric((geom_i))%name  ==  'PITCH') then
!
          e_z = c_pitch * c_yaw
          c_d = c_pitch * s_yaw * c_roll - s_pitch * s_roll
          do nf = 1, mnp
           do ptg_i = 1, j
             s_d = cosd(sc1(ptg_i))
             s_e = sc1(ptg_i) + elev_offst
             de_dx = c_d * (dsind(azim_angle) * dcosd(s_e) + dsind(s_e))  &
   &               - dcosd(azim_angle) * dcosd(s_e) * e_z
             de_dx_63 = c_d * (sind(azim_ref)*s_d + sind(sc1(ptg_i)))  &
   &                  - cosd(azim_ref) * s_d * e_z
             k_star_all(sv_elmnt,nf,ptg_i) = sc2(ptg_i) * deg2rad      &
   &                  * (sc3(ptg_i) * de_dx - sc4(ptg_i) * de_dx_63)
            end do
          end do
!
        else if (geometric((geom_i))%name  ==  'YAW') then
!
          c_d = s_pitch * s_yaw
          q = s_pitch * c_yaw
          e_z = q * s_roll
          s_d = q * c_roll
          do nf = 1, mnp
            do ptg_i = 1, j
              s_e = sc1(ptg_i) + elev_offst
              de_dx = dsind(azim_angle) * dcosd(s_e)*s_d - dsind(s_e)*e_z &
   &                + dcosd(azim_angle) * dcosd(s_e) * c_d
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
900 format('** Error in subroutine: i_k_star_no_conv ..',a,/,          &
   &       '   Exceeds maximum number of derivatives:',i5,/,           &
   &       '   New value should be at least:',i4,' large',/,           &
   &       '   Need to change it in: l2pc_file_paramters.inc',/,       &
   &       '   parameter (max_no_sv_derivatives =',i4,')')
!
    Return
  End subroutine I_K_STAR_NO_CONV
end module I_K_STAR_NO_CONV_M
! $Log$
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
!
