MODULE convolve_all_m

  USE MLSCommon, ONLY: I4, r4, R8, rp
  use intrinsic, only: L_VMR
  use Molecules, only: L_EXTINCTION
  USE Allocate_Deallocate, only: allocate_test, deallocate_test
  USE ForwardModelConfig, only: ForwardModelConfig_T
  use FOV_CONVOLVE_M, only: FOV_CONVOLVE
  USE VectorsModule, only: Vector_T, VectorValue_T, GetVectorQuantityByType
  USE AntennaPatterns_m, only: AntennaPattern_T
  USE MLSMessageModule, only: MLSMessage, MLSMSG_Error
  USE MatrixModule_0, only: M_ABSENT, M_BANDED, M_FULL
  USE MatrixModule_1, only: CREATEBLOCK, FINDBLOCK, MATRIX_T
  USE Load_sps_data_m, only: Grids_T

  IMPLICIT none
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
     "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName = &
     "$RCSfile$"
!---------------------------------------------------------------------------
 CONTAINS
! ============================================  convolve_all =====
! This subprogram adds the effects of antenna smearing to the radiance.
!
  SUBROUTINE convolve_all(FwdMdlConfig,FwdMdlIn,FwdMdlExtra,maf,channel,&
           & winStart,winFinish,mol_cat_index,temp,ptan,radiance,chi_in,&
           & rad_in,chi_out,dhdz_out,dx_dh_out,sbRatio,AntennaPattern,  &
           & t_deriv_flag,Grids_f,Jacobian,rowFlags,req,rsc,earth_frac, &
           & surf_angle,di_dt,dx_dt,d2x_dxdt,dxdt_tan,dxdt_surface,di_df)
!
! inputs
!
  Type (ForwardModelConfig_T), intent(in) :: FwdMdlCONFIG
  Type (Vector_t), intent(in) :: FwdMdlIN, FwdMdlEXTRA
!
  Integer, intent(in) :: MAF
  Integer, intent(in) :: CHANNEL
  Integer, intent(in) :: WINSTART
  Integer, intent(in) :: WINFINISH
  Integer, intent(IN) :: mol_cat_index(:)
!
  Type (VectorValue_T), intent(in) :: TEMP
  Type (VectorValue_T), intent(in) :: PTAN
!
  REAL(rp), INTENT(in) :: chi_in(:)    ! inputted pointing angles radians
  REAL(rp), INTENT(in) :: rad_in(:)    ! inputted radiances
  REAL(rp), INTENT(in) :: chi_out(:)   ! outputted pointing angles radians
  REAL(rp), INTENT(in) :: dhdz_out(:)  ! dhdz on the outputted pointing grid
  REAL(rp), INTENT(in) :: dx_dh_out(:) ! dx/dh on the outputted pointing grid
!
  Real(r8), intent(in) :: SbRatio

  Type(antennaPattern_T), intent(in) :: AntennaPattern
!
  Logical, dimension(:), pointer :: t_deriv_flag
!
  Type (Grids_T), intent(in) :: Grids_f
  Type (Matrix_T), OPTIONAL, intent(inout) :: Jacobian
!
  Logical, intent(inout) :: rowFlags(:) ! Flag to calling code
!
  REAL(rp), OPTIONAL, INTENT(in) :: req ! equivalent earth radius
  REAL(rp), OPTIONAL, INTENT(in) :: rsc ! spacecraft radius
  REAL(rp), OPTIONAL, INTENT(in) :: earth_frac ! fraction of earth in total
!                                   filled out pattern
! note req, rsc and earth_frac are non critical parameters and don't
! really need to be supplied externally. They are used to partition the
! full fft field between earth and space components.
!
! stuff for temperature derivatives
!
  REAL(rp), OPTIONAL, INTENT(in) :: surf_angle ! An angle that defines the
!                                   Earth surface.
  REAL(rp), OPTIONAL, INTENT(in) :: di_dt(:,:) ! derivative of radiance wrt
!                                   temperature on chi_in
  REAL(rp), OPTIONAL, INTENT(in) :: dx_dt(:,:,:) ! derivative of angle wrt
!                                   temperature on chi_in
  REAL(rp), OPTIONAL, INTENT(in) :: d2x_dxdt(:,:,:) ! 2nd derivative wrt angle
!                                   temperature on chi_in
  REAL(rp), OPTIONAL, INTENT(in) :: dxdt_tan(:,:,:) ! derivative of angle wrt
!                                   temperature on chi_in
  REAL(rp), OPTIONAL, INTENT(in) :: dxdt_surface(:,:,:) ! derivative of angle
!                                   wrt temperature at the surface
!
! stuff for atmospheric derivatives
!
  REAL(rp), OPTIONAL, INTENT(in) :: di_df(:,:) ! mixing ratio derivatives or
!                                   any parameter where a simple convolution
!                                   will suffice
! outputs
!
  Type (VectorValue_T), intent(inout) :: RADIANCE   ! Outputed radiances
!
! Internal stuff
!
  type (VectorValue_T), pointer :: F   ! An arbitrary species

  INTEGER(i4) :: j, k, Row, Col, ind, jz, sps_i, ptg_i, &
               & no_tan_hts, noPtan, noChans, jf

  Integer :: n_t_zeta, no_sv_p_t, sv_t_len, sv_f, f_len, no_mol
!
  REAL(r8) :: r, q
!
  REAL(r8), POINTER :: p(:)
  REAL(r8), POINTER :: dp(:)
  REAL(r8), POINTER :: angles(:)
  REAL(r8), POINTER :: rad_fft(:)
  REAL(r8), POINTER :: rad_fft1(:)
  REAL(r8), POINTER :: drad_dt_temp(:)

  REAL(r8), POINTER :: test1(:,:)
  REAL(r8), POINTER :: test2(:,:)
  REAL(r8), POINTER :: test3(:,:)
  REAL(r8), POINTER :: drad_dt_out(:,:)
  REAL(r8), POINTER :: drad_df_out(:,:)

  Real(r8) :: SRad(ptan%template%noSurfs), di_dx(ptan%template%noSurfs)
!
  n_t_zeta = temp%template%noSurfs
  no_sv_p_t = winFinish - winStart + 1
  sv_t_len = n_t_zeta * no_sv_p_t

  no_tan_hts = SIZE(chi_in)
  no_mol = SIZE(mol_cat_index)
  noChans = Radiance%template%noChans
  noPtan = ptan%template%nosurfs

  f_len = 0
  if ( PRESENT(di_df) ) f_len = SIZE(di_df,dim=2)
!
! Load the Radiance values into the Radiance structure:
!
  SRad = 0.0_r8
  CALL fov_convolve(antennaPattern,chi_in,rad_in,chi_out,SRad, &
                 &  drad_dx_out=di_dx)

  do ptg_i = 1, noPtan
    ind = channel + noChans * (ptg_i - 1)
    Radiance%values(ind,maf) = Radiance%values(ind,maf) + &
                                     &  sbRatio * SRad(ptg_i)
  end do
!
  j = 0
  IF ( PRESENT(Jacobian) ) j = 2
  if(j < 1) Return
!
! Compute dI/dPtan using the chain rule:
!
  SRad = sbRatio * di_dx * dx_dh_out * dhdz_out
!
! Now, load the dI/dPtan values into the Jacobian:
! (First, find index location Jacobian and write the derivative)
!
  row = FindBlock( Jacobian%row, radiance%index, maf )
  rowFlags(row) = .TRUE.
  col = FindBlock ( Jacobian%col, ptan%index, maf )

! Of course, we might not care about ptan

  if (col > 0) then

    select case (jacobian%block(Row,col)%kind)
      case (m_absent)
        call CreateBlock (Jacobian, row, col, m_banded, noPtan*noChans, &
                        & bandHeight=noChans)
        jacobian%block(row,col)%values = 0.0_r8
      case (m_banded)
      case default
        call MLSMessage (MLSMSG_Error, ModuleName,&
                      & 'Wrong matrix type for ptan derivative')
    end select

    do ptg_i = 1, noPtan
      ind = channel + noChans * (ptg_i-1)
      jacobian%block(row,col)%values(ind, 1) = &
              & jacobian%block(row,col)%values(ind, 1) + SRad(ptg_i)
      jacobian%block(row,col)%r1(ptg_i) = 1 + noChans * (ptg_i - 1)
      jacobian%block(row,col)%r2(ptg_i) = noChans * ptg_i
    end do

  endif

  IF (.not. ANY( (/FwdMdlConfig%temp_der, FwdMdlConfig%atmos_der, &
                &  FwdMdlConfig%spect_der/)) ) RETURN
!
  nullify (p,dp,angles,rad_fft,rad_fft1,drad_dt_temp,drad_dt_out, &
         & drad_df_out,test1,test2,test3) 
!
  IF (FwdMdlConfig%atmos_der .AND. .not. FwdMdlConfig%temp_der) THEN
!
    SRad = 0.0_r8
    CALL ALLOCATE_TEST(drad_df_out,noPtan,f_len,'drad_df_out',ModuleName)
    CALL fov_convolve(antennaPattern,chi_in,Rad_in,chi_out,SRad, &
                        & DI_DF=di_df, DRAD_DF_OUT=drad_df_out)
!
  ELSE IF (FwdMdlConfig%temp_der) THEN
!
    CALL ALLOCATE_TEST(test1,no_tan_hts,sv_t_len,'test1',Modulename)
    CALL ALLOCATE_TEST(test2,no_tan_hts,sv_t_len,'test2',Modulename)
    CALL ALLOCATE_TEST(test3,noPtan,    sv_t_len,'test3',Modulename)
    CALL ALLOCATE_TEST(drad_dt_out,noPtan,sv_t_len,'drad_dt_out',ModuleName)
!
    DO j = 1, no_sv_p_t
      DO k = 1, n_t_zeta
        r = dxdt_surface(1,j,k)
        test1(:,k+n_t_zeta*(j-1)) = dx_dt(:,j,k)
        test2(:,k+n_t_zeta*(j-1)) = d2x_dxdt(:,j,k)
! This arithmatic makes a surface value adjustment:
        test3(:,k+n_t_zeta*(j-1)) = dxdt_tan(:,j,k) - r
      end do
    end do
!
    IF (FwdMdlConfig%atmos_der) then
!
      CALL ALLOCATE_TEST(drad_df_out,noPtan,f_len,'drad_df_out',ModuleName)
      CALL fov_convolve(antennaPattern,chi_in,Rad_in,chi_out,SRad, &
         & SURF_ANGLE=surf_angle,DI_DT=di_dt,DX_DT=test1,DDX_DXDT=test2,&
         & DX_DT_OUT=test3,DRAD_DT_OUT=drad_dt_out,DI_DF=di_df, &
         & DRAD_DF_OUT=drad_df_out)
!
    ELSE
!
      CALL fov_convolve(antennaPattern,chi_in,Rad_in,chi_out,SRad, &
         & SURF_ANGLE=surf_angle,DI_DT=di_dt,DX_DT=test1, &
         & DDX_DXDT=test2,DX_DT_OUT=test3,DRAD_DT_OUT=drad_dt_out)
!
    ENDIF
!
! Load the Temp. derivative values into the Jacobian
!
    row = FindBlock( Jacobian%row, Radiance%index, maf )
    rowFlags(row) = .TRUE.

    sv_t_len = 0
    do jf = winStart, winFinish

      col = FindBlock ( Jacobian%col, temp%index, jf )
      select case ( Jacobian%block(row,col)%kind )
        case ( m_absent )
                call CreateBlock ( Jacobian, row, col, m_full )
                jacobian%block(row,col)%values = 0.0_r8
        case ( m_full )
        case default
                call MLSMessage ( MLSMSG_Error, ModuleName, &
                  & 'Wrong type for temperature derivative matrix' )
      end select

      do k = 1, n_t_zeta

! Check if derivatives are needed for this (zeta & phi) :

        sv_t_len = sv_t_len + 1
        if(.NOT. t_deriv_flag(sv_t_len)) CYCLE

! run through representation basis coefficients

        do ptg_i = 1, noPtan
          r = drad_dt_out(ptg_i,sv_t_len)
          ind = channel + noChans * (ptg_i-1)
          jacobian%block(row,col)%values(ind,k) =  &
               & jacobian%block(row,col)%values(ind,k) + sbRatio * r
        end do

      end do

    end do
!
    CALL DEALLOCATE_TEST(drad_dt_out,'drad_dt_out',ModuleName)
    CALL DEALLOCATE_TEST(test1,'test1',Modulename)
    CALL DEALLOCATE_TEST(test2,'test2',Modulename)
    CALL DEALLOCATE_TEST(test3,'test3',Modulename)

  ENDIF

  IF (.not. FwdMdlConfig%atmos_der) Return
!
! load Atmospheric derivatives into jacobian
!
  row = FindBlock( Jacobian%row, Radiance%index, maf )
  rowFlags(row) = .TRUE.
!
  sv_f = 0
  do sps_i = 1, no_mol
!
    jz = mol_cat_index(sps_i)
    if (FwdMdlConfig%molecules(jz) == l_extinction ) then
      f => GetVectorQuantityByType ( FwdMdlIn, FwdMdlExtra,quantityType= &
          l_extinction, radiometer=fwdMdlConfig%signals(1)%radiometer)
    else
      f => GetVectorQuantityByType ( FwdMdlIn, FwdMdlExtra, &
         & quantityType=l_vmr, molecule=FwdMdlConfig%molecules(jz))
    endif

    if(.not. associated(f) ) then
      jf = Grids_f%windowfinish(sps_i)-Grids_f%windowStart(sps_i)+1
      k = Grids_f%no_f(sps_i) * Grids_f%no_z(sps_i)
      sv_f = sv_f + jf * k
      CYCLE
    endif
!
    DO jf = Grids_f%windowStart(sps_i), Grids_f%windowfinish(sps_i)
!
      col = FindBlock ( Jacobian%col, f%index, jf)
      select case ( Jacobian%block(row,col)%kind )
        case ( m_absent )
          call CreateBlock ( Jacobian, row, col, m_full )
          jacobian%block(row,col)%values = 0.0_r8
        case ( m_full )
        case default
          call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Wrong type for temperature derivative matrix' )
      end select

      DO k = 1, Grids_f%no_f(sps_i) * Grids_f%no_z(sps_i)

! Check if derivatives are needed for this (zeta & phi) :

        sv_f = sv_f + 1
        if(.NOT. Grids_f%deriv_flags(sv_f)) CYCLE

! run through representation basis coefficients

        do ptg_i = 1, noPtan
          r = drad_df_out(ptg_i,sv_f)
          ind = channel + noChans * (ptg_i-1)
          jacobian%block(row,col)%values(ind,k) = &
              & jacobian%block(row,col)%values(ind,k) + sbRatio * r
        end do
!
      end do
!
    end do
!
  end do
!
  Call Deallocate_test ( drad_df_out, 'drad_df_out', ModuleName )
!
 END SUBROUTINE convolve_all

END MODULE convolve_all_m
! $Log$
! Revision 2.12  2002/06/19 11:00:31  zvi
! Removing debug statements
!
! Revision 2.10  2002/06/17 23:22:36  bill
! Add zvis modifications, some name changing
!
! Revision 2.9  2002/06/07 23:22:36  bill
! debugging study--wgr
!
! Revision 2.8  2002/06/07 04:50:25  bill
! fixes and improvements--wgr
!
! Revision 2.7  2002/06/04 10:28:02  zvi
! Encorporate deriv. flag into convolution, fixing a bug with 
! species ruuning index
!
! Revision 2.6  2002/05/22 19:42:44  zvi
! Fix a bug in the mol. index loop
!
! Revision 2.5  2002/02/15 22:52:16  livesey
! Bug fix for case where no ptan derivative
!
! Revision 2.4  2002/02/06 08:31:55  zvi
! Adding Temp. Deriv. correction
!
! Revision 2.3  2002/02/02 11:20:17  zvi
! Code to overwrite the l2cf integration & tanget grids
!
! Revision 2.2  2002/01/27 08:37:47  zvi
! Adding Users selected coefficients for derivatives
!
! Revision 2.1  2001/11/08 00:10:36  livesey
! Updated for extinction stuff
!
! Revision 2.0  2001/09/17 20:26:26  livesey
! New forward model
!
! Revision 1.29.2.1  2001/09/13 11:18:20  zvi
! Fix temp. derv. bug
!
! Revision 1.29  2001/08/24 03:42:26  jonathan
! change Ntr to Ntr-1 in do while loop
!
! Revision 1.28  2001/05/18 00:01:19  livesey
! Zero out some arrays to start with (mainly to make them safe to dump).
!
! Revision 1.27  2001/05/09 19:46:49  vsnyder
! Use new bandHeight argument of createBlock
!
! Revision 1.26  2001/05/03 02:03:16  vsnyder
! Insert copyright notice
!
! Revision 1.25  2001/05/02 20:49:23  zvi
! Cleaning up code
!
! Revision 1.24  2001/05/01 17:48:33  vsnyder
! Cosmetic changes -- put dummy arg declarations in same order as in header
!
! Revision 1.23  2001/05/01 00:42:54  zvi
! Fixing phi window bug
!
! Revision 1.22  2001/04/28 17:48:08  livesey
! Now accepts and sets rowFlags
!
! Revision 1.21  2001/04/27 22:37:54  vsnyder
! Don't compute derivatives if Jacobian isn't present
!
! Revision 1.20  2001/04/27 00:13:10  zvi
! Fixing some more phiwindow bug
!
! Revision 1.19  2001/04/26 22:54:41  zvi
! Fixing some phiwindow bug
!
! Revision 1.18  2001/04/20 23:09:13  livesey
! Cleaned up multi-channel case, also does folding in place
!
! Revision 1.17  2001/04/20 02:57:09  livesey
! Writes derivatives in matrix_t
!
! Revision 1.16  2001/04/19 23:56:52  livesey
! New parameters
!
! Revision 1.15  2001/04/10 10:14:16  zvi
! Fixing bug in convolve routines
!
! Revision 1.14  2001/04/10 02:25:14  livesey
! Tidied up some code
!
! Revision 1.13  2001/04/10 01:16:34  livesey
! Tidied up convolution
!
! Revision 1.12  2001/04/05 22:54:39  vsnyder
! Use AntennaPatterns_M
!
! Revision 1.11  2001/03/31 23:40:55  zvi
! Eliminate l2pcdim (dimension parameters) move to allocatable ..
!
! Revision 1.10  2001/03/29 12:08:17  zvi
! Fixing bugs
!
! Revision 1.9  2001/03/28 01:32:12  livesey
! Working version
!
! Revision 1.8  2001/03/28 00:40:01  zvi
! Fixing up convolution code, some minor changes in geoc_geod
!
! Revision 1.7  2001/03/26 17:56:14  zvi
! New codes to deal with dh_dt_path issue.. now being computed on the fly
!
! Revision 1.6  2001/03/21 01:10:31  livesey
! Now gets Ptan from vector
!
! Revision 1.5  2001/03/07 23:45:14  zvi
! Adding logical flags fro Temp, Atmos & Spect. derivatives
!
! Revision 1.1  2000/06/21 21:56:14  zvi
! First version D.P.
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
