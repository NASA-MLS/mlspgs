module GET_BETA_PATH_M
  use MLSCommon, only: I4, R8
  use L2PC_PFA_STRUCTURES, only: PFA_SLAB, MAXLINES
  use PATH_ENTITIES_M, only: PATH_VECTOR, PATH_BETA
  use SLABS_SW_M, only: SLABS_PREP_WDER, SLABS_PREP
  use CREATE_BETA_M, only: CREATE_BETA
  Implicit NONE
  private
  public :: get_beta_path
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
  "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
  "$RCSfile$"
!---------------------------------------------------------------------------
contains
!----------------------------------------------------------------------

 SUBROUTINE get_beta_path(ptg_i,pfa_spectrum,no_ele,no_ptg_frq, &
                     &    ptg_frq_grid,z_path,t_path,beta_path,ier)

!  ===============================================================
!  Declaration of variables for sub-program: get_beta_path
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: ptg_i, no_ele
Integer(i4), INTENT(IN) :: no_ptg_frq(*)

Integer(i4), INTENT(OUT) :: ier

Type(path_vector), INTENT(IN) :: ptg_frq_grid(*)
Type(path_vector), INTENT(IN) :: z_path, t_path

Type (pfa_slab), INTENT(IN) :: PFA_SPECTRUM(*)

Type(path_beta), INTENT(OUT) :: beta_path(:,:)  ! (sps_i,frq_i)

!  ----------------
!  Local variables:
!  ----------------

Integer(i4) :: nl, sps_i, no_sps, spectag, h_i, frq_i

Real(r8) :: Qlog(3), mass, z, p, t, Frq
Real(r8) :: values,t_power,dbeta_dw,dbeta_dn,dbeta_dnu

Real(r8) :: v0s(MAXLINES), x1(MAXLINES), y(MAXLINES), yi(MAXLINES), &
            slabs1(MAXLINES), dx1_dv0(MAXLINES), dy_dv0(MAXLINES),  &
            dslabs1_dv0(MAXLINES)

Real(r8) :: v0sp(MAXLINES), x1p(MAXLINES), yp(MAXLINES), yip(MAXLINES), &
            slabs1p(MAXLINES)

Real(r8) :: v0sm(MAXLINES), x1m(MAXLINES), ym(MAXLINES), yim(MAXLINES), &
            slabs1m(MAXLINES)

! Begin code:

  ier = 0
  no_sps = pfa_spectrum(1)%no_sps
!
! Allocate all the needed space for beta..
!
  do sps_i = 1, no_sps
    do frq_i = 1, no_ptg_frq(ptg_i)
      DEALLOCATE(beta_path(sps_i,frq_i)%values,      &
  &              beta_path(sps_i,frq_i)%t_power,     &
  &              beta_path(sps_i,frq_i)%dbeta_dw,    &
  &              beta_path(sps_i,frq_i)%dbeta_dn,    &
  &              beta_path(sps_i,frq_i)%dbeta_dnu, STAT = h_i)
      ALLOCATE(beta_path(sps_i,frq_i)%values(no_ele),    &
  &            beta_path(sps_i,frq_i)%t_power(no_ele),   &
  &            beta_path(sps_i,frq_i)%dbeta_dw(no_ele),  &
  &            beta_path(sps_i,frq_i)%dbeta_dn(no_ele),  &
  &            beta_path(sps_i,frq_i)%dbeta_dnu(no_ele), &
  &            STAT = ier)
      if(ier /= 0) then
        PRINT *,'** Allocation error in routine: get_beta_path ..'
        PRINT *,'   IER =',ier
        goto 99
      endif
    end do
  end do
!
  DO sps_i = 1, no_sps
!
    Spectag = pfa_spectrum(sps_i)%sps_spectag
    mass = Real(Spectag) / 1000.0
!
    nl = pfa_spectrum(sps_i)%NO_LINES
    Qlog(1:3) = pfa_spectrum(sps_i)%SPS_QLOG(1:3)

    do h_i = 1, no_ele
!
      z = z_path%values(h_i)
      if(z <= -4.5) CYCLE
!
      p = 10.0d0**(-z)
      t = t_path%values(h_i)
!
      Call Slabs_Prep_Arrays
!
      do frq_i = 1, no_ptg_frq(ptg_i)
!
        Frq = ptg_frq_grid(ptg_i)%values(frq_i)
!
        Call Create_beta (Spectag,p,t,Frq,nl,pfa_spectrum(sps_i), &
       &     v0s,x1,y,yi,slabs1,dx1_dv0,dy_dv0,dslabs1_dv0,v0sp,  &
       &     x1p,yp,yip,slabs1p,v0sm,x1m,ym,yim,slabs1m,values,   &
       &     t_power,dbeta_dw,dbeta_dn,dbeta_dnu,Ier)
        if(Ier /= 0) goto 99
!
        beta_path(sps_i,frq_i)%values(h_i) = values
        beta_path(sps_i,frq_i)%t_power(h_i) = t_power
        beta_path(sps_i,frq_i)%dbeta_dw(h_i) = dbeta_dw
        beta_path(sps_i,frq_i)%dbeta_dn(h_i) = dbeta_dn
        beta_path(sps_i,frq_i)%dbeta_dnu(h_i) = dbeta_dnu
!
      end do          ! On frq_i
!
    end do            ! On h_i
!
  END DO              ! On sps_i

  Return
!
! Cleanup cycle ...
!
 99  do sps_i = 1, no_sps
       do frq_i = 1, no_ptg_frq(ptg_i)
         DEALLOCATE(beta_path(sps_i,frq_i)%values,   &
        &           beta_path(sps_i,frq_i)%t_power,  &
        &           beta_path(sps_i,frq_i)%dbeta_dw, &
        &           beta_path(sps_i,frq_i)%dbeta_dn, &
        &           beta_path(sps_i,frq_i)%dbeta_dnu, STAT = h_i)
       end do
     end do

  Return

! *****     Internal procedures     **********************************
  contains
!
! --------------------------------  slabs_prep_arrays   -----
  Subroutine Slabs_Prep_Arrays
!
  Integer(i4) :: j
  Real(r8) :: v0,el,log_i,w,ps,n,n1,n2,gamma,delta,dslabs1,tp,tm

  if(Spectag==18999 .or. Spectag==28964 .or. Spectag==28965) Return
!
! Check for anything but liquid water and dry air:
!
  do j = 1, nl
!
    n = pfa_spectrum(sps_i)%SPS_N(j)
    w = pfa_spectrum(sps_i)%SPS_W(j)
    v0 = pfa_spectrum(sps_i)%SPS_V0(j)
    el = pfa_spectrum(sps_i)%SPS_EL(j)
    ps = pfa_spectrum(sps_i)%SPS_PS(j)
    n1 = pfa_spectrum(sps_i)%SPS_N1(j)
    n2 = pfa_spectrum(sps_i)%SPS_N2(j)
    log_i = pfa_spectrum(sps_i)%SPS_STR(j)
    gamma = pfa_spectrum(sps_i)%SPS_GAMMA(j)
    delta = pfa_spectrum(sps_i)%SPS_DELTA(j)
!
! Prepare the temperature weighted coefficients:
!
    Call Slabs_prep_wder(t,mass,v0,el,w,ps,p,n,log_i,Qlog,delta, &
   &           gamma,n1,n2,v0s(j),x1(j),y(j),yi(j),slabs1(j),&
   &           dx1_dv0(j),dy_dv0(j),dslabs1_dv0(j))
!
    tp = t + 10.0
    Call slabs_prep(tp,mass,v0,el,w,ps,p,n,log_i,Qlog,delta,gamma, &
   &           n1,n2,v0sp(j),x1p(j),yp(j),yip(j),slabs1p(j),dslabs1)
!
    tm = t - 10.0
    Call slabs_prep(tm,mass,v0,el,w,ps,p,n,log_i,Qlog,delta,gamma, &
   &           n1,n2,v0sm(j),x1m(j),ym(j),yim(j),slabs1m(j),dslabs1)
!
  end do
!
  End Subroutine Slabs_Prep_Arrays
!
 END SUBROUTINE get_beta_path
end module GET_BETA_PATH_M
! $Log$
! Revision 1.1 2001/02/01 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
