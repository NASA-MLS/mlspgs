module GET_CHI_ANGLES_M
  use MLSCOMMON, only: I4, R8
  use L2PCDIM, only: Nlvl
  use L2PC_FILE_PARAMETERS, only: mxco => MAX_NO_ELMNTS_PER_SV_COMPONENT, &
                                  DEG2RAD
  use PATH_ENTITIES_M, only: PATH_INDEX, PATH_VECTOR
  use GET_ETA_M, only: GET_ETA
  use D_LINTRP_M, only: LINTRP
  Implicit NONE
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains

!----------------------------------------------------------------------
! Set up array of pointing angles

SUBROUTINE get_chi_angles(ndx_path,n_path,tan_press,tan_hts,tan_temp, &
           phi_tan,RoC,h_obs,elev_offset,dh_dt_grid,N_lvls,no_tan_hts,&
           no_t,t_z_basis,z_grid,si,cen_ang,ptg_angle,dx_dt,d2x_dxdt,ier)

!  ===============================================================
!  Declaration of variables for sub-program: get_chi_angles
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: no_tan_hts,no_t,N_lvls,si

Real(r8), INTENT(IN) :: t_z_basis(:), z_grid(:), dh_dt_grid(:,:)
Real(r8), INTENT(IN) :: tan_hts(:), tan_temp(:), tan_press(:)
Real(r8), INTENT(IN) :: RoC, phi_tan, h_obs, elev_offset

Type(path_index), INTENT(IN)  :: ndx_path(:)
Type(path_vector), INTENT(IN) :: n_path(:)

Integer(i4), INTENT(OUT) :: ier
Real(r8), INTENT(OUT) :: ptg_angle(:), cen_ang
Real(r8), INTENT(OUT) :: dx_dt(:,:), d2x_dxdt(:,:)

!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: i, brkpt, ptg_i, sv_i, k, is, h_i

! Real(r8), PARAMETER :: ampl = 38.91
! Real(r8), PARAMETER :: phas = 51.6152 * deg2rad

Real(r8), PARAMETER :: ampl = 38.9014
Real(r8), PARAMETER :: phas = 51.6814 * deg2rad

Real(r8) :: r,t,dh,tanx,cse,ht,Rs_eq,ngrid,schi
Real(r8) :: Eta(Nlvl,mxco), temp_dh_dt(Nlvl,mxco)

  Rs_eq = h_obs + ampl * Sin(2.0*(phi_tan-phas))    ! ** Experimental

! Get 'No_tan_hts' convolution angles

  DO i = 1, no_tan_hts
    ht = tan_hts(i)
    brkpt = ndx_path(i)%break_point_index
    ngrid = n_path(i)%values(brkpt) - n_path(i)%values(1)
    schi = (1.0 + ngrid) * (ht + RoC) / Rs_eq
    IF(ABS(schi) > 1.0) THEN
      ier = 1
      PRINT *,'*** Error in get_chi_angles subroutine !'
      PRINT *,'    arg > 1.0 in ArcSin(arg) ..'
      RETURN
    END IF
    ptg_angle(i) = DAsin(schi) + elev_offset
  END DO

  cen_ang = ptg_angle(si)
!
! Set up: dx_dt, d2x_dxdt arrays for temperature derivative computations
! (NOTE: These entities has NO PHI dimension, so take the center Phi in dh_dt)
!
!  First: Get table of temperature basis functions
!
  Call get_eta(tan_press,t_z_basis,no_tan_hts,no_t,Nlvl,Eta)
!
  temp_dh_dt(1:si+1,1:no_t) = 0.0
!
  ptg_i = no_tan_hts - si + 1
  do sv_i = 1, no_t
    Call Lintrp(z_grid,tan_press(si:no_tan_hts),dh_dt_grid(1:,sv_i), &
   &            temp_dh_dt(si:no_tan_hts,sv_i),N_lvls,ptg_i)
  end do
!
  k = 2 * N_lvls
  do sv_i = 1, no_t
!
    is = 1
    r = dh_dt_grid(1,sv_i)
!
    do while (tan_hts(is)  <  0.0)
      t = tan_temp(is)
      dh = tan_hts(is) + RoC
      eta(is,sv_i) = 0.0
      tanx = Tan(ptg_angle(is))
      cse = tanx * tanx
      dx_dt(is,sv_i) = r * tanx / dh
      d2x_dxdt(is,sv_i) = (2.0+cse)*r/dh + Eta(is,sv_i)/t
      is = is + 1
    end do
!
    do h_i = is, no_tan_hts
      t = tan_temp(h_i)
      dh = tan_hts(h_i) + RoC
      tanx = tan(ptg_angle(h_i))
      cse = tanx * tanx
      r = temp_dh_dt(h_i,sv_i)
      dx_dt(h_i,sv_i) = r * tanx / dh
      d2x_dxdt(h_i,sv_i) = (2.0+cse)*r/dh + Eta(h_i,sv_i)/t
    end do
!
  end do

  RETURN
END SUBROUTINE get_chi_angles
end module GET_CHI_ANGLES_M
! $Log$
! Revision 1.1 2000/06/09 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
