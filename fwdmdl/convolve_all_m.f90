! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Convolve_All_m

  implicit NONE
  private
  public :: CONVOLVE_ALL

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
 contains
! ============================================  convolve_all =====
! This subprogram adds the effects of antenna smearing to the radiance.

  subroutine Convolve_All ( FwdMdlConfig, FwdMdlIn, FwdMdlExtra, maf, channel, &
           & winStart, winFinish, Qtys, temp, ptan, radiance, update, chi_in,  &
           & rad_in, chi_out, dhdz_out, dx_dh_out, sbRatio, AntennaPattern,    &
           & t_deriv_flag, Grids_f, Jacobian, rowFlags, req, rsc, earth_frac,  &
           & surf_angle, di_dt, dx_dt, d2x_dxdt, dxdt_tan, dxdt_surface, di_df,&
           & ptan_Der )

    use Allocate_Deallocate, only: allocate_test, deallocate_test
    use AntennaPatterns_m, only: AntennaPattern_T
    use ForwardModelConfig, only: ForwardModelConfig_T
    use ForwardModelVectorTools, only: QtyStuff_T
    use Fov_Convolve_m, only: Fov_Convolve
    use Load_sps_data_m, only: Grids_T
    use MatrixModule_0, only: M_ABSENT, M_BANDED, M_FULL, CHECKFORSIMPLEBANDEDLAYOUT
    use MatrixModule_1, only: CREATEBLOCK, FINDBLOCK, MATRIX_T
    use MLSCommon, only: I4, R8, RP, RM
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use VectorsModule, only: Vector_T, VectorValue_T

! inputs

    type (ForwardModelConfig_T), intent(in) :: FwdMdlCONFIG
    type (Vector_t), intent(in) :: FwdMdlIN, FwdMdlEXTRA

    integer, intent(in) :: MAF
    integer, intent(in) :: CHANNEL
    integer, intent(in) :: WINSTART
    integer, intent(in) :: WINFINISH
    type(QtyStuff_T), intent(in) :: Qtys(:)

    type (vectorvalue_t), intent(in) :: TEMP
    type (vectorvalue_t), intent(in) :: PTAN

    real(rp), intent(in) :: chi_in(:)    ! inputted pointing angles radians
    real(rp), intent(in) :: rad_in(:)    ! inputted radiances
    real(rp), intent(in) :: chi_out(:)   ! outputted pointing angles radians
    real(rp), intent(in) :: dhdz_out(:)  ! dhdz on the outputted pointing grid
    real(rp), intent(in) :: dx_dh_out(:) ! dx/dh on the outputted pointing grid

    real(r8), intent(in) :: SbRatio

    type(antennaPattern_T), intent(in) :: AntennaPattern

    logical, dimension(:), pointer :: t_deriv_flag

    type (Grids_T), intent(in) :: Grids_f
    type (Matrix_t), optional, intent(inout) :: Jacobian

    logical, intent(inout) :: rowFlags(:) ! Flag to calling code

    real(rp), optional, intent(in) :: Req ! equivalent earth radius
    real(rp), optional, intent(in) :: Rsc ! spacecraft radius
    real(rp), optional, intent(in) :: Earth_frac ! fraction of earth in total
!                                   filled out pattern
    logical, intent(in), optional :: Ptan_der ! Flag
    logical, intent(in) :: Update       ! If set, just add to radiance don't overwrite

! note req, rsc and earth_frac are non critical parameters and don't
! really need to be supplied externally. They are used to partition the
! full fft field between earth and space components.

! stuff for temperature derivatives

    real(rp), optional, intent(in) :: surf_angle ! An angle that defines the
!                                   Earth surface.
    real(rp), optional, intent(in) :: di_dt(:,:) ! derivative of radiance wrt
!                                   temperature on chi_in
    real(rp), optional, intent(in) :: dx_dt(:,:) ! derivative of angle wrt
!                                   temperature on chi_in
    real(rp), optional, intent(in) :: d2x_dxdt(:,:) ! 2nd derivative wrt angle
!                                   temperature on chi_in
    real(rp), optional, intent(in) :: dxdt_tan(:,:) ! derivative of angle wrt
!                                   temperature on chi_in
    real(rp), optional, intent(in) :: dxdt_surface(:,:) ! derivative of angle
!                                   wrt temperature at the surface

! stuff for atmospheric derivatives

    real(rp), optional, intent(in) :: di_df(:,:) ! mixing ratio derivatives or
!                                   any parameter where a simple convolution
!                                   will suffice
! outputs

    type (VectorValue_T), intent(inout) :: RADIANCE   ! Outputed radiances

! Internal stuff

    integer(i4) :: k, Row, Col, ind, sps_i, ptg_i, noPtan, noChans, jf, nfz

    integer :: n_t_zeta, no_sv_p_t, sv_t_len, sv_f, f_len

    real(r8) :: R

    real(r8), pointer :: drad_dt_out(:,:)
    real(r8), pointer :: drad_df_out(:,:)

    real(r8), pointer :: temp_dxdt_tan(:,:)

    real(r8) :: SRad(ptan%template%noSurfs), di_dx(ptan%template%noSurfs)

    logical :: my_ptan_der

    my_ptan_der = .false.
    if ( present ( ptan_der ) ) my_ptan_der = ptan_der
    n_t_zeta = temp%template%noSurfs
    no_sv_p_t = winFinish - winStart + 1
    sv_t_len = n_t_zeta * no_sv_p_t

    noChans = Radiance%template%noChans
    noPtan = ptan%template%nosurfs

    f_len = 0
    if ( present(di_df) ) f_len = SIZE(di_df,dim=2)

! Load the Radiance values into the Radiance structure:

    SRad = 0.0_r8
    call fov_convolve ( antennaPattern, chi_in, rad_in, chi_out, SRad,  &
                     &  drad_dx_out=di_dx )

    if ( update ) then
      do ptg_i = 1, noPtan
        ind = channel + noChans * (ptg_i - 1)
        Radiance%values(ind,maf) = Radiance%values(ind,maf) + &
          &  sbRatio * SRad(ptg_i)
      end do
    else
      do ptg_i = 1, noPtan
        ind = channel + noChans * (ptg_i - 1)
        Radiance%values(ind,maf) = sbRatio * SRad(ptg_i)
      end do
    end if

    if ( .not. present(Jacobian) ) return

! Compute dI/dPtan using the chain rule:

    SRad = sbRatio * di_dx * dx_dh_out * dhdz_out

! Now, load the dI/dPtan values into the Jacobian:
! (First, find index location Jacobian and write the derivative)

    row = FindBlock( Jacobian%row, radiance%index, maf )
    rowFlags(row) = .TRUE.

! Of course, we might not care about ptan

    if ( my_ptan_der ) then
      col = FindBlock ( Jacobian%col, ptan%index, maf )

      select case (jacobian%block(Row,col)%kind)
        case (m_absent)
          call CreateBlock ( Jacobian, row, col, m_banded, noPtan*noChans, &
                           & bandHeight=noChans )
          jacobian%block(row,col)%values = 0.0_rm
        case (m_banded)
          call CheckForSimpleBandedLayout ( jacobian%block(row,col), noChans, &
            & 'd[radiance]/d[ptan] in convolution' )
        case default
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Wrong matrix type for ptan derivative' )
      end select

      if ( update ) then 
        do ptg_i = 1, noPtan
          ind = channel + noChans * (ptg_i-1)
          r = jacobian%block(row,col)%values(ind,1)
          Jacobian%block(row,col)%values(ind, 1) = r + SRad(ptg_i)
        end do
      else
        do ptg_i = 1, noPtan
          ind = channel + noChans * (ptg_i-1)
          Jacobian%block(row,col)%values(ind, 1) = SRad(ptg_i)
        end do
      end if
    end if

    if ( .not. ANY( (/FwdMdlConfig%temp_der, FwdMdlConfig%atmos_der, &
                  &  FwdMdlConfig%spect_der/)) ) return

    nullify ( drad_dt_out, drad_df_out, temp_dxdt_tan )

    if ( FwdMdlConfig%atmos_der .AND. .not. FwdMdlConfig%temp_der ) then

      SRad = 0.0_r8
      call allocate_test ( drad_df_out, noPtan, f_len, 'drad_df_out', ModuleName )
      call fov_convolve ( antennaPattern, chi_in, Rad_in, chi_out, SRad,  &
                        & DI_DF=di_df,  DI_DF_FLAG=grids_f%deriv_flags,  &
                        & DRAD_DF_OUT=drad_df_out )

    else if ( FwdMdlConfig%temp_der ) then

      call allocate_test ( drad_dt_out, noPtan, sv_t_len, 'drad_dt_out', &
                         & ModuleName )

      call allocate_test ( temp_dxdt_tan, noPtan, sv_t_len, 'temp_dxdt_tan', &
                        &  ModuleName )

      temp_dxdt_tan = dxdt_tan - SPREAD(dxdt_surface(1,:),1,noPtan)

      if ( FwdMdlConfig%atmos_der ) then

        call allocate_test ( drad_df_out, noPtan, f_len, 'drad_df_out', ModuleName)
        call fov_convolve ( antennaPattern, chi_in, Rad_in, chi_out, SRad,  &
           & SURF_ANGLE=surf_angle, DI_DT=di_dt, DX_DT=dx_dt,               &
           & DDX_DXDT=d2x_dxdt, DX_DT_OUT=temp_dxdt_tan,                    &
           & DRAD_DT_OUT=drad_dt_out, DI_DF=di_df, DI_DT_FLAG=t_deriv_flag, &
           & DI_DF_FLAG=grids_f%deriv_flags, DRAD_DF_OUT=drad_df_out )
             

      else

        call fov_convolve ( antennaPattern, chi_in, Rad_in, chi_out, SRad,  &
           & SURF_ANGLE=surf_angle, DI_DT=di_dt, DX_DT=dx_dt,               &
           & DDX_DXDT=d2x_dxdt, DX_DT_OUT=temp_dxdt_tan,                    &
           & DI_DT_FLAG=t_deriv_flag, DRAD_DT_OUT=drad_dt_out )

      end if

! Load the Temp. derivative values into the Jacobian

      row = FindBlock( Jacobian%row, Radiance%index, maf )
      rowFlags(row) = .TRUE.

      sv_t_len = 0
      do jf = winStart, winFinish

        col = FindBlock ( Jacobian%col, temp%index, jf )
        select case ( Jacobian%block(row,col)%kind )
          case ( m_absent )
            call CreateBlock ( Jacobian, row, col, m_full )
            jacobian%block(row,col)%values = 0.0_rm
          case ( m_full )
          case default
            call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Wrong type block for temperature derivative matrix' )
        end select

        do k = 1, n_t_zeta

! Check if derivatives are needed for this (zeta & phi) :

          sv_t_len = sv_t_len + 1
          if ( .NOT. t_deriv_flag(sv_t_len)) cycle

! run through representation basis coefficients
        
          if ( update ) then
            do ptg_i = 1, noPtan
              r = drad_dt_out(ptg_i,sv_t_len)
              ind = channel + noChans * (ptg_i-1)
              jacobian%block(row,col)%values(ind,k) =  &
                & jacobian%block(row,col)%values(ind,k) + sbRatio * r
            end do
          else
            do ptg_i = 1, noPtan
              r = drad_dt_out(ptg_i,sv_t_len)
              ind = channel + noChans * (ptg_i-1)
              jacobian%block(row,col)%values(ind,k) = sbRatio * r
            end do
          end if
        end do

      end do

      call deallocate_test ( drad_dt_out, 'drad_dt_out', ModuleName )
      call deallocate_test ( temp_dxdt_tan, 'temp_dxdt_tan', ModuleName )

    end if

    if ( .not. FwdMdlConfig%atmos_der) Return

! load Atmospheric derivatives into jacobian

    row = FindBlock( Jacobian%row, Radiance%index, maf )
    rowFlags(row) = .TRUE.

    do sps_i = 1, size(qtys)

      if ( .not. qtys(sps_i)%foundInFirst ) cycle

      sv_f = grids_f%l_v(sps_i-1)
      nfz = (Grids_f%l_f(sps_i) - Grids_f%l_f(sps_i-1)) * &
          & (Grids_f%l_z(sps_i) - Grids_f%l_z(sps_i-1))

      do jf = Grids_f%windowStart(sps_i), Grids_f%windowfinish(sps_i)

        col = FindBlock ( Jacobian%col, qtys(sps_i)%qty%index, jf)
        select case ( Jacobian%block(row,col)%kind )
          case ( m_absent )
            call CreateBlock ( Jacobian, row, col, m_full )
            jacobian%block(row,col)%values = 0.0_rm
          case ( m_full )
          case default
            call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Wrong type block for atmospheric derivative matrix' )
        end select

        do k = 1, nfz
! Check if derivatives are needed for this (zeta & phi) :

          sv_f = sv_f + 1
          if ( .NOT. Grids_f%deriv_flags(sv_f) ) cycle

! run through representation basis coefficients

          if ( update ) then 
            do ptg_i = 1, noPtan
              r = drad_df_out(ptg_i,sv_f)
              ind = channel + noChans * (ptg_i-1)
              jacobian%block(row,col)%values(ind,k) = &
                & jacobian%block(row,col)%values(ind,k) + sbRatio * r
            end do
          else
            do ptg_i = 1, noPtan
              r = drad_df_out(ptg_i,sv_f)
              ind = channel + noChans * (ptg_i-1)
              jacobian%block(row,col)%values(ind,k) = sbRatio * r
            end do
          end if
        end do

      end do

    end do

    call deallocate_test ( drad_df_out, 'drad_df_out', ModuleName )

  end subroutine Convolve_All

  logical function NOT_USED_HERE()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function NOT_USED_HERE

end module Convolve_All_m

! $Log$
! Revision 2.30  2003/05/19 19:58:07  vsnyder
! Remove USEs for unreferenced symbols, remove unused local variables
!
! Revision 2.29  2003/05/05 23:00:25  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.28.2.2  2003/03/21 02:47:03  vsnyder
! Use an array of pointers to quantities instead of GetQuantityForForwardModel
!
! Revision 2.28.2.1  2003/03/20 01:42:26  vsnyder
! Revise Grids_T structure
!
! Revision 2.28  2002/11/13 17:07:17  livesey
! Bug fix, now takes FwdMdlExtra
!
! Revision 2.27  2002/10/10 19:36:22  vsnyder
! Add di_dt_flag
!
! Revision 2.26  2002/10/10 01:00:10  vsnyder
! Delete unreferenced entities from USE list
!
! Revision 2.25  2002/10/10 00:54:16  vsnyder
! Give Id a value, sort USEs, cosmetic changes
!
! Revision 2.24  2002/10/08 17:08:01  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.23  2002/09/26 18:01:52  livesey
! Now uses GetQuantityForForwardModel.
!
! Revision 2.22  2002/09/11 17:43:39  pwagner
! Began changes needed to conform with matrix%values type move to rm from r8
!
! Revision 2.21  2002/09/10 17:05:45  livesey
! Added update argument
!
! Revision 2.20  2002/09/07 00:52:24  vsnyder
! Cosmetic changes, copyright notice
!
! Revision 2.19  2002/09/05 20:48:13  livesey
! Fixed handling of vmr derivative flags.
!
! Revision 2.18  2002/08/20 22:36:39  livesey
! Moved uses inside routine
!
! Revision 2.17  2002/07/23 22:26:38  livesey
! Added ptan_der handling
!
! Revision 2.16  2002/07/16 08:47:11  mjf
! Nullified temp_dxdt_tan along with drad_dt_out, etc.
!
! Revision 2.15  2002/07/08 17:45:38  zvi
! Remove unnecessary variables
!
! Revision 2.14  2002/07/05 07:52:46  zvi
! Coor. switch (phi,z) -> (z,phi)
!
! Revision 2.13  2002/06/28 11:06:47  zvi
! compute dI/dPtan using chain rule
!
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
