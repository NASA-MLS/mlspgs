! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module GET_BETA_PATH_M
  use MLSCommon, only: I4, R8
  use GLNP, only: NG
  use SpectroscopyCatalog_m, only: Catalog_T, Lines
  use L2PC_PFA_STRUCTURES, only: PFA_SLAB, SLABS_STRUCT, MAXLINES
  use PATH_ENTITIES_M, only: PATH_VECTOR, PATH_BETA, PATH_INDEX
  use CREATE_BETA_M, only: CREATE_BETA
  Implicit NONE
  private
  public :: GET_COARSE_BETA_PATH, GET_GLBETA_PATH, GET_BETA_PATH

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
  "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
  "$RCSfile$"
!---------------------------------------------------------------------------
contains
!----------------------------------------------------------------------

 SUBROUTINE get_coarse_beta_path(ptg_i,Catalog,ndx_path,gl_slabs, &
         &  no_frq,mnf,ht,frq_grid,z_path,h_path,t_path,Frq_Gap,  &
         &  temp_der, spect_der, beta_path, Ier)

!  ===============================================================
!  Declaration of variables for sub-program: get_coarse_beta_path
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: ptg_i, no_frq, mnf

Real(r8), INTENT(IN) :: frq_grid(:)
Real(r8), INTENT(IN) :: ht, Frq_Gap
Logical,  INTENT(IN) :: temp_der, spect_der

Integer(i4), INTENT(OUT) :: ier

Type(path_vector), INTENT(IN) :: z_path, h_path, t_path

Type(path_index) :: ndx_path
Type (slabs_struct), DIMENSION(:,:), POINTER :: gl_slabs

Type(Catalog_T), INTENT(IN) :: Catalog(:)

Type(path_beta), POINTER :: beta_path(:,:)  ! (sps_i,frq_i)

!  ----------------
!  Local variables:
!  ----------------

Integer(i4) :: nl,i,j,k,mp,no_sps,spectag,frq_i,Ngp1,no_ele,brkpt

Real(r8) :: z, p, h, t, Frq
Real(r8) :: values,t_power,dbeta_dw,dbeta_dn,dbeta_dnu

! Begin code:

  ier = 0
  Ngp1 = Ng + 1
  no_sps = Size(Catalog)
!
  brkpt = ndx_path%break_point_index
  no_ele = ndx_path%total_number_of_elements
!
! Allocate all the needed space for beta.
! NOTE: This is done onl once per MAF, at the lowest pointing,
! where the lenght of the path is maximal, also, allocate beta
! with the maximum number of frequencies (mnf) over all pointings..
!
  if(ptg_i == 1) then
!
    if ( associated(beta_path) ) then
      do i = 1, size(beta_path,1)
        do frq_i = 1, size(beta_path,2)
          DEALLOCATE(beta_path(i,frq_i)%values,   &
            &        beta_path(i,frq_i)%t_power,  &
            &        beta_path(i,frq_i)%dbeta_dw, &
            &        beta_path(i,frq_i)%dbeta_dn, &
            &        beta_path(i,frq_i)%dbeta_dnu, STAT=j )
        end do
      end do
      DEALLOCATE(beta_path,STAT=j)
    end if
!
    ALLOCATE(beta_path(no_sps,mnf),STAT=j)
!
    do i = 1, size(beta_path,1)
      do frq_i = 1, size(beta_path,2)
        ALLOCATE(beta_path(i,frq_i)%values(no_ele),    &
           &     beta_path(i,frq_i)%t_power(no_ele),   &
           &     beta_path(i,frq_i)%dbeta_dw(no_ele),  &
           &     beta_path(i,frq_i)%dbeta_dn(no_ele),  &
           &     beta_path(i,frq_i)%dbeta_dnu(no_ele), STAT = ier)
        if(ier /= 0) then
          PRINT *,'** Allocation error in routine: get_beta_path ..'
          PRINT *,'   no_ele,ptg_i,i,frq_i:',no_ele,ptg_i,i,frq_i
          PRINT *,'   STAT =',ier
          Return
        endif
        beta_path(i,frq_i)%values = 0.0
        beta_path(i,frq_i)%t_power = 0.0
        beta_path(i,frq_i)%dbeta_dw = 0.0
        beta_path(i,frq_i)%dbeta_dn = 0.0
        beta_path(i,frq_i)%dbeta_dnu = 0.0
      end do
    end do
!
  endif
!
  DO i = 1, no_sps
!
    Spectag = Catalog(i)%spec_tag
!
    nl = gl_slabs(1,i)%no_lines
!
    mp = 1 - Ngp1
    do
!
      mp = mp + Ngp1
      if (mp > brkpt) EXIT
!
      z = z_path%values(mp)
      if(z <= -4.5) CYCLE
!
      j = mp
      t = t_path%values(mp)
      p = 10.0d0**(-z)
!
      do frq_i = 1, no_frq
!
        Frq = frq_grid(frq_i)
!
        Call Create_beta (Spectag,p,t,Frq,nl,Catalog(i),gl_slabs(j,i)%v0s, &
       &     gl_slabs(j,i)%x1,gl_slabs(j,i)%y,gl_slabs(j,i)%yi,            &
       &     gl_slabs(j,i)%slabs1,gl_slabs(j,i)%dx1_dv0,                   &
       &     gl_slabs(j,i)%dy_dv0,gl_slabs(j,i)%dslabs1_dv0,               &
       &     gl_slabs(j,i)%v0sp,gl_slabs(j,i)%x1p,gl_slabs(j,i)%yp,        &
       &     gl_slabs(j,i)%yip,gl_slabs(j,i)%slabs1p,gl_slabs(j,i)%v0sm,   &
       &     gl_slabs(j,i)%x1m,gl_slabs(j,i)%ym,gl_slabs(j,i)%yim,         &
       &     gl_slabs(j,i)%slabs1m,values,t_power,dbeta_dw,dbeta_dn,       &
       &     dbeta_dnu, Frq_Gap, temp_der, spect_der, Ier)
        if(Ier /= 0) Return
!
        beta_path(i,frq_i)%values(mp) = values
        beta_path(i,frq_i)%t_power(mp) = t_power
        beta_path(i,frq_i)%dbeta_dw(mp) = dbeta_dw
        beta_path(i,frq_i)%dbeta_dn(mp) = dbeta_dn
        beta_path(i,frq_i)%dbeta_dnu(mp) = dbeta_dnu
!
      end do          ! On frq_i
!
    end do
!
    mp = brkpt + 1
    h = h_path%values(mp)
    do while (h < ht - 0.0001)
      mp = mp + Ngp1
      h = h_path%values(mp)
    end do
!
    mp = mp - Ngp1
!
    do
!
      mp = mp + Ngp1
      if(mp > no_ele) Return
!
      z = z_path%values(mp)
      if(z <= -4.5) CYCLE
!
      j = mp
      t = t_path%values(mp)
      p = 10.0d0**(-z)
!
      do frq_i = 1, no_frq
!
        Frq = frq_grid(frq_i)
!
        Call Create_beta (Spectag,p,t,Frq,nl,Catalog(i),gl_slabs(j,i)%v0s, &
       &     gl_slabs(j,i)%x1,gl_slabs(j,i)%y,gl_slabs(j,i)%yi,            &
       &     gl_slabs(j,i)%slabs1,gl_slabs(j,i)%dx1_dv0,                   &
       &     gl_slabs(j,i)%dy_dv0,gl_slabs(j,i)%dslabs1_dv0,               &
       &     gl_slabs(j,i)%v0sp,gl_slabs(j,i)%x1p,gl_slabs(j,i)%yp,        &
       &     gl_slabs(j,i)%yip,gl_slabs(j,i)%slabs1p,gl_slabs(j,i)%v0sm,   &
       &     gl_slabs(j,i)%x1m,gl_slabs(j,i)%ym,gl_slabs(j,i)%yim,         &
       &     gl_slabs(j,i)%slabs1m,values,t_power,dbeta_dw,dbeta_dn,       &
       &     dbeta_dnu, Frq_Gap, temp_der, spect_der, Ier)
!
        beta_path(i,frq_i)%values(mp) = values
        beta_path(i,frq_i)%t_power(mp) = t_power
        beta_path(i,frq_i)%dbeta_dw(mp) = dbeta_dw
        beta_path(i,frq_i)%dbeta_dn(mp) = dbeta_dn
        beta_path(i,frq_i)%dbeta_dnu(mp) = dbeta_dnu
!
      end do          ! On frq_i
!
    end do
!
  END DO              ! On i

  Return
!
 END SUBROUTINE get_coarse_beta_path

!----------------------------------------------------------------------

 SUBROUTINE get_glbeta_path(ptg_i,frq_i,Catalog,gl_ndx,no_gl_ndx,gl_slabs, &
         &  Frq,z_path,t_path,Frq_Gap,temp_der,spect_der,beta_path,Ier)

!  ===============================================================
!  Declaration of variables for sub-program: get_glbeta_path
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: ptg_i,frq_i, no_gl_ndx
Integer(i4), INTENT(IN) :: gl_ndx(:,:)
Logical,     INTENT(IN) :: temp_der, spect_der

Real(r8), INTENT(IN) :: Frq, Frq_Gap

Integer(i4), INTENT(OUT) :: ier

Type (slabs_struct), DIMENSION(:,:), POINTER :: gl_slabs

Type(path_vector), INTENT(IN) :: z_path, t_path

Type(Catalog_T), INTENT(IN) :: Catalog(:)

Type(path_beta), POINTER :: beta_path(:,:)  ! (sps_i,frq_i)

!  ----------------
!  Local variables:
!  ----------------

Integer(i4) :: nl, i, j, no_sps, spectag, Ngp1, mp, jj

Real(r8) :: z, p, t,values,t_power,dbeta_dw,dbeta_dn,dbeta_dnu

! Begin code:

  ier = 0
  if(no_gl_ndx < 1) Return
!
  Ngp1 = Ng + 1
  no_sps = Size(Catalog)
!
  DO i = 1, no_sps
!
    Spectag = Catalog(i)%spec_tag
!
    nl = gl_slabs(1,i)%no_lines
!
    do jj = 1, no_gl_ndx
!
      mp = gl_ndx(jj,2)
      j = mp

      do                     !  GL Inner loop
!
        j = j + 1
        if(j >= gl_ndx(jj+1,2)) EXIT

        z = z_path%values(j)
        if(z <= -4.5) CYCLE
!
        p = 10.0d0**(-z)
        t = t_path%values(j)
!
        Call Create_beta (Spectag,p,t,Frq,nl,Catalog(i),gl_slabs(j,i)%v0s, &
       &     gl_slabs(j,i)%x1,gl_slabs(j,i)%y,gl_slabs(j,i)%yi,            &
       &     gl_slabs(j,i)%slabs1,gl_slabs(j,i)%dx1_dv0,                   &
       &     gl_slabs(j,i)%dy_dv0,gl_slabs(j,i)%dslabs1_dv0,               &
       &     gl_slabs(j,i)%v0sp,gl_slabs(j,i)%x1p,gl_slabs(j,i)%yp,        &
       &     gl_slabs(j,i)%yip,gl_slabs(j,i)%slabs1p,gl_slabs(j,i)%v0sm,   &
       &     gl_slabs(j,i)%x1m,gl_slabs(j,i)%ym,gl_slabs(j,i)%yim,         &
       &     gl_slabs(j,i)%slabs1m,values,t_power,dbeta_dw,dbeta_dn,       &
       &     dbeta_dnu, Frq_Gap, temp_der, spect_der, Ier)
        if(Ier /= 0) Return
!
        beta_path(i,frq_i)%values(j) = values
        beta_path(i,frq_i)%t_power(j) = t_power
        beta_path(i,frq_i)%dbeta_dw(j) = dbeta_dw
        beta_path(i,frq_i)%dbeta_dn(j) = dbeta_dn
        beta_path(i,frq_i)%dbeta_dnu(j) = dbeta_dnu
!
      end do          ! On GL inner loop
!
    end do            ! On jj
!
  END DO              ! On i
!
  Return
!
! Cleanup cycle ...
!
 99  Print *,'** Error in get_glbeta routine ..'

  Return

 END SUBROUTINE get_glbeta_path

!----------------------------------------------------------------------
! Computes beta over the entire GL path without referring to the gl_ndx at all.

 subroutine Get_beta_path ( frequencies, Catalog, no_ele, z_path, t_path, &
      &     beta_path, vel_z, Frq_Gap, temp_der, spect_der, Ier)

 use SLABS_SW_M, only: SLABS_PREP_ARRAYS

  !  ===============================================================
  !  Declaration of variables for sub-program: get_beta_path
  !  ===============================================================
  !  ---------------------------
  !  Calling sequence variables:
  !  ---------------------------

  Real(r8), dimension(:), intent(in) :: frequencies
  Integer(i4), intent(in) :: no_ele
  Logical, intent(in) :: temp_der, spect_der

  Real(r8),    intent(in) :: vel_z, Frq_Gap

  Integer(i4), intent(out) :: ier

  Type(path_vector), intent(in) :: z_path, t_path

  Type(Catalog_T), dimension(:), intent(in) :: Catalog

  Type(path_beta), pointer :: beta_path(:,:)  ! (sps_i,frq_i)

  !  ----------------
  !  Local variables:
  !  ----------------

  Real(r8), parameter :: C = 299792.4583_r8    ! Speed of Light Km./Sec.

  Integer(i4) :: nl, i, no_sps, mnf, spectag, h_i, frq_i

  Real(r8) :: Qlog(3), mass, z, p, t, Frq, Vel_z_correction
  Real(r8) :: values,t_power,dbeta_dw,dbeta_dn,dbeta_dnu

  Real(r8) :: y(50), ym(50), yp(50)
  Real(r8) :: yi(50), yim(50), yip(50)
  Real(r8) :: x1(50), x1m(50), x1p(50)
  Real(r8) :: v0s(50), v0sm(50), v0sp(50)
  Real(r8) :: slabs1(50), slabs1m(50), slabs1p(50)
  Real(r8) :: dy_dv0(50), dx1_dv0(50), dslabs1_dv0(50)

  ! Begin code:

  ier = 0
  no_sps = Size(Catalog)
  mnf =  size(frequencies)

  ! call output('In get_beta_path_m',advance='yes')
  ! call dump(frequencies)

  Vel_z_correction = 1.0_r8 - vel_z / c

  ! Allocate all the needed space for beta..

  if ( associated(beta_path) ) then
    do i = 1, size(beta_path,1)
      do frq_i = 1, size(beta_path,2)
        deallocate ( beta_path(i,frq_i)%values, beta_path(i,frq_i)%t_power, &
          & beta_path(i,frq_i)%dbeta_dw, beta_path(i,frq_i)%dbeta_dn, &
          & beta_path(i,frq_i)%dbeta_dnu, STAT=h_i )
      end do
    end do
  end if

  deallocate ( beta_path, stat=h_i )
  allocate ( beta_path(no_sps,mnf), stat=h_i )

  do i = 1, no_sps
    do frq_i = 1, mnf
      allocate(beta_path(i,frq_i)%values(no_ele),    &
  &            beta_path(i,frq_i)%t_power(no_ele),   &
  &            beta_path(i,frq_i)%dbeta_dw(no_ele),  &
  &            beta_path(i,frq_i)%dbeta_dn(no_ele),  &
  &            beta_path(i,frq_i)%dbeta_dnu(no_ele), STAT = ier)
      if ( ier /= 0 ) then
        print *,'** Allocation error in routine: get_beta_path ..'
        print *,'   no_ele,i,frq_i:',no_ele,i,frq_i
        print *,'   STAT =',ier
        goto 99
      end if
    end do
  end do

  do i = 1, no_sps

    Spectag = Catalog(i)%spec_tag
    mass = Real(Spectag) / 1000.0

    nl = Size(Catalog(i)%lines)
    Qlog(1:3) = Catalog(i)%QLOG(1:3)

    do h_i = 1, no_ele

      z = z_path%values(h_i)
      if(z <= -4.5) CYCLE

      p = 10.0**(-z)
      t = t_path%values(h_i)
!
      Call  Slabs_Prep_Arrays(Spectag,nl,t,p,mass,Catalog(i),Qlog,v0s, &
        &   x1,y,yi,slabs1,dx1_dv0,dy_dv0,dslabs1_dv0,v0sp,x1p,yp,yip, &
        &   slabs1p,v0sm,x1m,ym,yim,slabs1m)

! Apply velocity corrections:

      v0s(1:nl)  = v0s(1:nl)  * Vel_z_correction
      v0sp(1:nl) = v0sp(1:nl) * Vel_z_correction
      v0sm(1:nl) = v0sm(1:nl) * Vel_z_correction

      do frq_i = 1, mnf

        Frq = frequencies(frq_i)

        Call Create_beta (Spectag,p,t,Frq,nl,Catalog(i),v0s,x1, &
       &     y,yi,slabs1,dx1_dv0,dy_dv0,dslabs1_dv0,v0sp,x1p,yp,&
       &     yip,slabs1p,v0sm,x1m,ym,yim,slabs1m,values,t_power,&
       &     dbeta_dw,dbeta_dn,dbeta_dnu,Frq_Gap,temp_der,spect_der,Ier)
        if(Ier /= 0) goto 99

        beta_path(i,frq_i)%values(h_i) = values
        beta_path(i,frq_i)%t_power(h_i) = t_power
        beta_path(i,frq_i)%dbeta_dw(h_i) = dbeta_dw
        beta_path(i,frq_i)%dbeta_dn(h_i) = dbeta_dn
        beta_path(i,frq_i)%dbeta_dnu(h_i) = dbeta_dnu

      end do          ! On frq_i

    end do            ! On h_i

  end do              ! On i

  Return

! Cleanup cycle ...

 99  do i = 1, no_sps
       do frq_i = 1, mnf
         deallocate(beta_path(i,frq_i)%values,beta_path(i,frq_i)%t_power, &
           &        beta_path(i,frq_i)%dbeta_dw,beta_path(i,frq_i)%dbeta_dn, &
           &        beta_path(i,frq_i)%dbeta_dnu,STAT=h_i )
       end do
     end do

     deallocate ( beta_path, STAT=h_i )

    Return

  end subroutine get_beta_path
!
End module GET_BETA_PATH_M
! $Log$
! Revision 1.21  2001/05/16 00:42:10  livesey
! Fixed velocity correction.  Make sure this is OK with people later.
!
! Revision 1.20  2001/05/15 03:47:26  zvi
! Adding derivative flag to beta calculations
!
! Revision 1.19  2001/05/14 23:14:54  zvi
! Added Freq. Gap test..
!
! Revision 1.18  2001/05/03 22:17:32  vsnyder
! Changed some d0 constants to _r8, cosmetic changes
!
! Revision 1.17  2001/05/02 20:49:23  zvi
! Cleaning up code
!
! Revision 1.16  2001/04/05 21:58:47  zvi
! Implementing l2cf inputs for FilterShape & Spectroscopy instead of FMI
!
! Revision 1.15  2001/04/03 07:32:45  zvi
! Modify the spectral structure - eliminating sps_ from the names
!
! Revision 1.14  2001/03/31 23:40:55  zvi
! Eliminate l2pcdim (dimension parameters) move to allocatable ..
!
! Revision 1.13  2001/03/29 08:51:01  zvi
! Changing the (*) toi (:) everywhere
!
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
