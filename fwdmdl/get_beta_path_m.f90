module GET_BETA_PATH_M
  use MLSCommon, only: I4, R8
  use L2PC_PFA_STRUCTURES, only: PFA_SLAB, MAXLINES
  use PATH_ENTITIES_M, only: PATH_VECTOR, PATH_BETA
  use SLABS_SW_M, only: SLABS_PREP_WDER, SLABS_PREP
  use CREATE_BETA_M, only: CREATE_BETA
  use output_m,only:output
  use Dump_0,only:dump
  Implicit NONE
  private
  public :: get_beta_path
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
  "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
  "$RCSfile$"
!---------------------------------------------------------------------------
contains
!----------------------------------------------------------------------

 SUBROUTINE get_beta_path(frequencies,pfs,no_ele, &
                       &  z_path,t_path,beta_path,vel_z,ier)

!  ===============================================================
!  Declaration of variables for sub-program: get_beta_path
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
real(r8), dimension(:), intent(in) :: frequencies
Integer(i4), INTENT(IN) :: no_ele
Real(r8),    INTENT(IN) :: vel_z

Integer(i4), INTENT(OUT) :: ier

Type(path_vector), INTENT(IN) :: z_path, t_path

Type (pfa_slab), INTENT(IN) :: pfs(:)

Type(path_beta), POINTER :: beta_path(:,:)  ! (sps_i,frq_i)

!  ----------------
!  Local variables:
!  ----------------

Real(r8), PARAMETER :: c = 299792.4583d0     ! Speed of Light Km./Sec.

Integer(i4) :: nl, i, no_sps, mnf, spectag, h_i, frq_i

Real(r8) :: Qlog(3), mass, z, p, t, Frq, Vel_z_correction
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
  no_sps = pfs(1)%no_sps
  mnf =  size(frequencies)

! call output('In get_beta_path_m',advance='yes')
! call dump(frequencies)
!
  Vel_z_correction = 1.0_r8 + vel_z / c
!
! Allocate all the needed space for beta..
!
  if ( associated(beta_path) ) then
    do i = 1, size(beta_path,1)
      do frq_i = 1, size(beta_path,2)
        deallocate ( beta_path(i,frq_i)%values, beta_path(i,frq_i)%t_power, &
          & beta_path(i,frq_i)%dbeta_dw, beta_path(i,frq_i)%dbeta_dn, &
          & beta_path(i,frq_i)%dbeta_dnu, STAT=h_i )
      end do
    end do
  end if

  DEALLOCATE(beta_path,STAT=h_i)
  ALLOCATE(beta_path(no_sps,mnf),STAT=h_i)
!
  do i = 1, no_sps
    do frq_i = 1, mnf
      ALLOCATE(beta_path(i,frq_i)%values(no_ele),    &
  &            beta_path(i,frq_i)%t_power(no_ele),   &
  &            beta_path(i,frq_i)%dbeta_dw(no_ele),  &
  &            beta_path(i,frq_i)%dbeta_dn(no_ele),  &
  &            beta_path(i,frq_i)%dbeta_dnu(no_ele), STAT = ier)
      if(ier /= 0) then
        PRINT *,'** Allocation error in routine: get_beta_path ..'
        PRINT *,'   no_ele,i,frq_i:',no_ele,i,frq_i
        PRINT *,'   STAT =',ier
        goto 99
      endif
    end do
  end do
!
  DO i = 1, no_sps
!
    Spectag = pfs(i)%sps_spectag
    mass = Real(Spectag) / 1000.0
!
    nl = pfs(i)%NO_LINES
    Qlog(1:3) = pfs(i)%SPS_QLOG(1:3)

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
! Apply velocity corrections:
!
      v0s(1:nl)  = v0s(1:nl)  * Vel_z_correction
      v0sp(1:nl) = v0sp(1:nl) * Vel_z_correction
      v0sm(1:nl) = v0sm(1:nl) * Vel_z_correction
!
      do frq_i = 1, mnf
!
        Frq = frequencies(frq_i)
!
        Call Create_beta (Spectag,p,t,Frq,nl,pfs(i),v0s,x1,y,yi,&
       &     slabs1,dx1_dv0,dy_dv0,dslabs1_dv0,v0sp,x1p,yp,yip, &
       &     slabs1p,v0sm,x1m,ym,yim,slabs1m,values,t_power,    &
       &     dbeta_dw,dbeta_dn,dbeta_dnu,Ier)
        if(Ier /= 0) goto 99
!
        beta_path(i,frq_i)%values(h_i) = values
        beta_path(i,frq_i)%t_power(h_i) = t_power
        beta_path(i,frq_i)%dbeta_dw(h_i) = dbeta_dw
        beta_path(i,frq_i)%dbeta_dn(h_i) = dbeta_dn
        beta_path(i,frq_i)%dbeta_dnu(h_i) = dbeta_dnu
!
      end do          ! On frq_i
!
    end do            ! On h_i
!
  END DO              ! On i

  Return
!
! Cleanup cycle ...
!
 99  do i = 1, no_sps
       do frq_i = 1, mnf
         DEALLOCATE(beta_path(i,frq_i)%values,beta_path(i,frq_i)%t_power,&
        &           beta_path(i,frq_i)%dbeta_dw,beta_path(i,frq_i)%dbeta_dn,&
        &           beta_path(i,frq_i)%dbeta_dnu, STAT=h_i)
       end do
     end do
     DEALLOCATE(beta_path, STAT=h_i)
     stop
  Return

! *****     Internal procedures     **********************************
  contains
!
! --------------------------------  slabs_prep_arrays   -----
  Subroutine Slabs_Prep_Arrays
!
  Integer(i4) :: j
  Real(r8) :: dslabs1,tp,tm

  if(Spectag==18999 .or. Spectag==28964 .or. Spectag==28965) Return
!
! Check for anything but liquid water and dry air:
!
  do j = 1, nl
!
! Prepare the temperature weighted coefficients:
!
    Call Slabs_prep_wder(t,mass,pfs(i)%SPS_V0(j),pfs(i)%SPS_EL(j),&
   &     pfs(i)%SPS_W(j),pfs(i)%SPS_PS(j),p,pfs(i)%SPS_N(j),      &
   &     pfs(i)%SPS_STR(j),Qlog,pfs(i)%SPS_DELTA(j),pfs(i)%SPS_GAMMA(j),& 
   &     pfs(i)%SPS_N1(j),pfs(i)%SPS_N2(j),v0s(j),x1(j),y(j),yi(j),&
   &     slabs1(j),dx1_dv0(j),dy_dv0(j),dslabs1_dv0(j))
!
    tp = t + 10.0
    Call slabs_prep(tp,mass,pfs(i)%SPS_V0(j),pfs(i)%SPS_EL(j), &
   &     pfs(i)%SPS_W(j),pfs(i)%SPS_PS(j),p,pfs(i)%SPS_N(j),   &
   &     pfs(i)%SPS_STR(j),Qlog,pfs(i)%SPS_DELTA(j),pfs(i)%SPS_GAMMA(j),&
   &     pfs(i)%SPS_N1(j),pfs(i)%SPS_N2(j),v0sp(j),x1p(j),yp(j),&
   &     yip(j),slabs1p(j),dslabs1)
!
    tm = t - 10.0
    Call slabs_prep(tm,mass,pfs(i)%SPS_V0(j),pfs(i)%SPS_EL(j), &
   &     pfs(i)%SPS_W(j),pfs(i)%SPS_PS(j),p,pfs(i)%SPS_N(j),   &
   &     pfs(i)%SPS_STR(j),Qlog,pfs(i)%SPS_DELTA(j),pfs(i)%SPS_GAMMA(j),&
   &     pfs(i)%SPS_N1(j),pfs(i)%SPS_N2(j),v0sm(j),x1m(j),ym(j),&
   &     yim(j),slabs1m(j),dslabs1)
!
  end do
!
  End Subroutine Slabs_Prep_Arrays
!
 END SUBROUTINE get_beta_path
end module GET_BETA_PATH_M
! $Log$
! Revision 1.12  2001/03/20 23:22:40  zvi
! Change to new geoc_geod routine..
!
! Revision 1.11  2001/03/20 11:03:16  zvi
! Fixing code for "real" data run, increase dim. etc.
!
! Revision 1.10  2001/03/20 02:29:26  livesey
! Interim version, gets same numbers as zvi
!
! Revision 1.9  2001/03/15 12:18:03  zvi
! Adding the Velocity effect on Line center frequency
!
! Revision 1.8  2001/03/09 02:26:11  vsnyder
! More work on deallocation
!
! Revision 1.7  2001/03/09 02:11:28  vsnyder
! Repair deallocating
!
! Revision 1.6  2001/03/05 21:37:20  zvi
! New filter format
!
! Revision 1.1 2001/02/01 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
