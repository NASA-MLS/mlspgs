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
  USE MLSNumerics, ONLY: interpolatevalues, hunt
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
! ============================================  new_convolve_all =====
! This subprogram adds the effects of antenna smearing to the radiance.
!
  SUBROUTINE convolve_all(FwdMdlConfig,FwdMdlIn,FwdMdlExtra,maf,channel,&
           & winStart,winFinish,mol_cat_index,temp,ptan,radiance,chi_in,   &
           & rad_in,chi_out,sbRatio,AntennaPattern,t_deriv_flag,Grids_f,    &
           & Jacobian,rowFlags,req,rsc,earth_frac,surf_angle,di_dt, &
           & dx_dt,d2x_dxdt,dxdt_tan,dxdt_surface,di_df)
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
  REAL(rp), INTENT(in) :: chi_in(:)! inputted pointing angles radians
  REAL(rp), INTENT(in) :: rad_in(:)! inputted radiances
  REAL(rp), INTENT(in) :: chi_out(:)! outputted pointing angles radians
!
  Real(r8), intent(in) :: SBRATIO

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

  Integer :: n_t_zeta, no_sv_p_t, sv_t_len, sv_f, sv_i, f_len, no_mol
!
  Logical :: Want_Deriv
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

  Real(r8) :: SRad(ptan%template%noSurfs), Term(ptan%template%noSurfs)
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
  Term = 0.0_r8
  CALL fov_convolve(antennaPattern,chi_in,rad_in,chi_out,Term)
  SRad = sbRatio * Term
  do ptg_i = 1, noPtan
    ind = channel + noChans * (ptg_i - 1)
    Radiance%values(ind,maf) = Radiance%values(ind,maf) + SRad(ptg_i)
  end do
!
  want_deriv = PRESENT(jacobian) .AND. ANY( (/FwdMdlConfig%temp_der, &
               & FwdMdlConfig%atmos_der,FwdMdlConfig%spect_der/) )

! Find out if user wants pointing derivatives
    if ( .not. want_deriv ) Return
!
! Compute Ptan derivatives:
!
  Term(:) = 0.0_r8
  do ptg_i = 1, noPtan
    j = 1
    if(ptg_i == noPtan) j = -1
    q = Ptan%values(ptg_i+j,maf) -  Ptan%values(ptg_i,maf)
    Term(ptg_i) = (SRad(ptg_i+j) - SRad(ptg_i) ) / q
  end do
!
! Find index location Jacobian and write the derivative
!
  row = FindBlock( Jacobian%row, radiance%index, maf )
  rowFlags(row) = .true.
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
              & jacobian%block(row,col)%values(ind, 1) + Term(ptg_i)
      jacobian%block(row,col)%r1(ptg_i) = 1 + noChans * (ptg_i - 1)
      jacobian%block(row,col)%r2(ptg_i) = noChans * ptg_i
    end do
  end if
!
  IF (.not. FwdMdlConfig%atmos_der .AND. .not. FwdMdlConfig%temp_der) RETURN

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
    rowFlags(row) = .true.

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
  rowFlags(row) = .true.
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

      sv_i = 0
      DO k = 1, Grids_f%no_f(sps_i) * Grids_f%no_z(sps_i)

! Check if derivatives are needed for this (zeta & phi) :

        sv_f = sv_f + 1
        if(.NOT. Grids_f%deriv_flags(sv_f)) CYCLE

! *** ZEBUG
!       if(channel == 13) then
!         if(ANY(drad_df_out(:,sv_f) /= 0.0)) then
!           ind=Grids_f%no_f(sps_i)
!           Print *,' ch,sps_i,kf,sv_f:',channel,sps_i,ind,sv_f
!           Print *,' Window width:', &
!             Grids_f%windowfinish(sps_i)-Grids_f%windowStart(sps_i)+1
!           Print *,' drad_df_out(ptg_i,sv_f),ptg_i=1,noPtan'
!           Print 905,(drad_df_out(ptg_i,sv_f),ptg_i=1,noPtan)
!         endif
!       endif
!905    format ( 5(2x, 1pg13.6) )
! *** END ZEBUG

! run through representation basis coefficients

        sv_i = sv_i + 1
        do ptg_i = 1, noPtan
          ind = channel + noChans * (ptg_i-1)
          jacobian%block(row,col)%values(ind,sv_i) = &
                      jacobian%block(row,col)%values(ind,sv_i) + &
                          & sbRatio * drad_df_out(ptg_i,sv_f)
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
