! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module NO_CONV_AT_ALL_M

  implicit NONE
  private
  public :: No_Conv_At_All

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
  !-------------------------------------------------------------------------
  ! This subroutine transfers the derivatives over from the internal
  ! convolution grid to the users specified points. This module uses
  ! cubic spline interpolation to do the job.

  Subroutine No_Conv_At_All ( FwmConf, ForwardModelIn, ForwardModelExtra, maf, &
           & Channel, WindowStart, WindowFinish, Temp, Ptan, Radiance, update, &
           & t_deriv_flag, ptg_angles, chi_out, dhdz_out, dx_dh_out, Grids_f,  &
           & I_raw, sbRatio, qtys, rowFlags, Jacobian, di_dt, di_df, ptan_Der )

    use ForwardModelConfig, only: ForwardModelConfig_T
    use ForwardModelVectorTools, only: QtyStuff_T
    use Load_sps_data_m, only: Grids_T
    use MatrixModule_0, only: M_ABSENT, M_BANDED, M_FULL, CHECKFORSIMPLEBANDEDLAYOUT
    use MatrixModule_1, only: CREATEBLOCK, FINDBLOCK, MATRIX_T
    use MLSCommon, only: R8, RM
    use MLSMessageModule, only: MLSMSG_Error, MLSMessage
    use MLSNumerics, ONLY: INTERPOLATEVALUES
    use VectorsModule, only: Vector_T, VectorValue_T

    type (ForwardModelConfig_T) :: FWMCONF
    type (Vector_T), intent(in) :: FORWARDMODELIN
    type (Vector_T), intent(in) :: FORWARDMODELEXTRA

    integer, intent(in) :: MAF
    integer, intent(in) :: CHANNEL
    integer, intent(in) :: WINDOWSTART
    integer, intent(in) :: WINDOWFINISH
    type(qtyStuff_T), intent(in) :: qtys(:)

    type (VectorValue_T), intent(in) :: TEMP
    type (VectorValue_T), intent(in) :: PTAN
    type (VectorValue_T), intent(inout) :: RADIANCE

    type (Grids_T), Intent(in) :: Grids_f

    real(r8), intent(in) :: sbRatio
    real(r8), intent(in) :: i_raw(:), ptg_angles(:), chi_out(:), dhdz_out(:), &
                         &  dx_dh_out(:)
    logical, intent(in), optional :: ptan_Der     ! Flag
!
! derivative of radiance w.r.t. temperature on chi_in
    real(r8), optional, intent(in) :: di_dt(:,:)

! mixing ratio derivatives or any parameter which behaves like VMR
    real(r8), optional, intent(in) :: di_df(:,:)

    type (matrix_t), intent(inout), optional :: Jacobian
!
    logical, dimension(:), pointer :: t_deriv_flag
    logical, intent(in) :: Update       ! If set just add to radiacnes/derivatives

    logical, dimension(:), intent(inout) :: rowFlags ! Flag to calling code

    ! -----     Local Variables     ------------------------------------

    integer:: No_t, No_tan_hts

    integer :: jf, jz, nfz
    integer :: is, k, nf, sv_f, sv_t_len
    integer :: Row, col                     ! Matrix row & column indices
    integer :: ptg_i,noPtan,noChans,Ind     ! Indices

    real(r8) :: Rad( size(ptg_angles)), q
    real(r8) :: SRad(ptan%template%noSurfs)
    real(r8) :: di_dx(ptan%template%noSurfs)
    real(r8) :: I_star_all(ptan%template%noSurfs)
    logical :: my_ptan_der

    ! -----  Begin the code  -------------------------------------------

    my_ptan_der = .false.
    if ( present ( ptan_der ) ) my_ptan_der = ptan_der
    no_t = temp%template%noSurfs
    no_tan_hts = size(ptg_angles)

    noPtan = ptan%template%noSurfs
    noChans = radiance%template%noChans

! Ptan derivative

    col = 0
    if ( present (Jacobian) ) then
      row = FindBlock ( Jacobian%row, radiance%index, maf )
      rowFlags(row) = .TRUE.
    end if

! Of course, we might not care about ptan

    if ( my_ptan_der ) then
      col = FindBlock ( Jacobian%col, ptan%index, maf )

      call InterpolateValues ( ptg_angles, i_raw, chi_out, i_star_all, &
                             & METHOD='S', dyByDx=di_dx )

! Use the chain rule to compute dI/dPtan on the output grid:

      SRad(1:noPtan) = di_dx(1:noPtan) * dx_dh_out(1:noPtan)
      SRad(1:noPtan) = SRad(1:noPtan) * dhdz_out(1:noPtan)

      select case ( Jacobian%block(Row,col)%kind )
        case ( m_absent )
          call CreateBlock ( Jacobian, row, col, m_banded, noPtan*noChans, &
                           & bandHeight=noChans, init=0.0_rm )
        case ( m_banded )
          call CheckForsimpleBandedLayout ( jacobian%block(row,col), noChans, &
            & 'd[Radiance]/d[ptan] in no convolution case' )
        case default
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Wrong matrix type for ptan derivative' )
      end select

      if ( update ) then
        do ptg_i = 1, noPtan
          ind = channel + noChans*(ptg_i-1)
          q = Jacobian%block(row,col)%values(ind,1)
          Jacobian%block(row,col)%values(ind,1) = q + sbRatio * SRad(ptg_i)
        end do
      else
        do ptg_i = 1, noPtan
          ind = channel + noChans*(ptg_i-1)
          Jacobian%block(row,col)%values(ind,1) = sbRatio * SRad(ptg_i)
        end do
      end if
    else

      call InterpolateValues ( ptg_angles,i_raw,chi_out,i_star_all,METHOD='S', &
                             & extrapolate = 'C' )

    end if

    if ( update ) then
      do ptg_i = 1, noPtan
        ind = channel + noChans*(ptg_i-1)
        radiance%values( ind, maf ) = &
          & radiance%values(ind,maf) + sbRatio * i_star_all(ptg_i)
      end do
    else
      do ptg_i = 1, noPtan
        ind = channel + noChans*(ptg_i-1)
        radiance%values( ind, maf ) = sbRatio * i_star_all(ptg_i)
      end do
    end if

    if ( .not. present(Jacobian) ) return

    if ( .not. any((/fwmConf%temp_der, &
                   & fwmConf%atmos_der, &
                   & fwmConf%spect_der/)) ) return

    ! Now transfer the other fwd_mdl derivatives to the output pointing
    ! values

    ! ********************* Temperature derivatives ******************

    ! check to determine if derivative is desired for this parameter

    if ( fwmConf%temp_der ) then

    ! Derivatives needed continue to process

      sv_t_len = 0
      Rad(1:) = 0.0
      SRad(1:) = 0.0
      k = no_tan_hts

      do nf = windowStart, windowFinish

        col = FindBlock ( Jacobian%col, temp%index, nf )
        select case ( Jacobian%block(row,col)%kind )
        case ( m_absent )
          call CreateBlock ( Jacobian, row, col, m_full, init=0.0_rm )
        case ( m_full )
        case default
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Wrong type for temperature derivative matrix' )
        end select

        do jz = 1, no_t

! Check if derivatives are needed for this (zeta & phi) :

          sv_t_len = sv_t_len + 1
          if ( .NOT. t_deriv_flag(sv_t_len) ) cycle

          Rad(1:k) = di_dt(1:k,sv_t_len)
          call InterpolateValues ( ptg_angles, Rad, chi_out, Srad, 'S', &
                                 & extrapolate = 'C' )
          if ( update ) then
            do ptg_i = 1, noPtan
              ind = channel + noChans*(ptg_i-1)
              q = Jacobian%block(row,col)%values(ind,jz)
              Jacobian%block(row,col)%values(ind,jz) = q + sbRatio*Srad(ptg_i)
            end do
          else
            do ptg_i = 1, noPtan
              ind = channel + noChans*(ptg_i-1)
              Jacobian%block(row,col)%values(ind,jz) = sbRatio*Srad(ptg_i)
            end do
          end if
        end do

      end do

    end if

    if ( fwmConf%atmos_der ) then

      ! ****************** atmospheric derivatives ******************

      do is = 1, size(qtys)
        
        if ( .not. qtys(is)%foundInFirst ) cycle

        sv_f = grids_f%l_v(is-1)
        nfz = (Grids_f%l_f(is) - Grids_f%l_f(is-1)) * &
            & (Grids_f%l_z(is) - Grids_f%l_z(is-1))

        do jf = Grids_f%windowStart(is), Grids_f%windowfinish(is)

          col = FindBlock ( Jacobian%col, qtys(is)%qty%index, jf)
          select case ( Jacobian%block(row,col)%kind )
            case ( m_absent )
              call CreateBlock ( Jacobian, row, col, m_full, init=0.0_rm )
            case ( m_full )
            case default
              call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Wrong type for atmospheric derivative matrix' )
          end select

          do k = 1, nfz

! Check if derivatives are needed for this (zeta & phi) :

            sv_f = sv_f + 1
            if ( .NOT. Grids_f%deriv_flags(sv_f) ) cycle

            Rad(1:no_tan_hts) = di_df(1:no_tan_hts,sv_f)
            call InterpolateValues  ( ptg_angles, Rad, chi_out, Srad, 'L', &
                                    & extrapolate = 'C' )
            if ( update ) then
              do ptg_i = 1, noPtan
                ind = channel + noChans*(ptg_i-1)
                q = Jacobian%block(row,col)%values(ind,k)
                Jacobian%block(row,col)%values(ind,k) = q + sbRatio*Srad(ptg_i)
              end do
            else
              do ptg_i = 1, noPtan
                ind = channel + noChans*(ptg_i-1)
                Jacobian%block(row,col)%values(ind,k) = sbRatio*Srad(ptg_i)
              end do
            end if

          end do

        end do

      end do

    end if

  end subroutine No_Conv_At_All

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module NO_CONV_AT_ALL_M
! $Log$
! Revision 2.23  2003/10/09 22:17:30  livesey
! Added call to CheckForSimpleBandedLayout
!
! Revision 2.22  2003/05/17 01:17:03  vsnyder
! Remove unused names, futzing
!
! Revision 2.21  2003/05/16 23:52:36  livesey
! Removed obsolete spectag stuff
!
! Revision 2.20  2003/05/05 23:00:26  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.19.2.2  2003/03/21 02:47:03  vsnyder
! Use an array of pointers to quantities instead of GetQuantityForForwardModel
!
! Revision 2.19.2.1  2003/03/20 01:42:26  vsnyder
! Revise Grids_T structure
!
! Revision 2.19  2002/11/13 17:07:28  livesey
! Bug fix, now takes forwardModelExtra
!
! Revision 2.18  2002/10/08 17:08:05  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.17  2002/10/04 23:46:21  vsnyder
! Cosmetic changes
!
! Revision 2.16  2002/09/26 18:02:03  livesey
! Now uses GetQuantityForForwardModel.
!
! Revision 2.15  2002/09/11 17:43:39  pwagner
! Began changes needed to conform with matrix%values type move to rm from r8
!
! Revision 2.14  2002/09/10 17:05:52  livesey
! Added update argument
!
! Revision 2.13  2002/08/21 23:38:56  bill
!  added no extrapolate to interpolation calls
!
! Revision 2.12  2002/08/20 22:36:47  livesey
! Moved uses inside routine
!
! Revision 2.11  2002/07/29 23:16:32  bill
! got rid of debugging write
!
! Revision 2.10  2002/07/29 21:42:02  bill
! no changes, just debugging
!
! Revision 2.9  2002/07/23 22:26:52  livesey
! Added ptan_der handling
!
! Revision 2.8  2002/07/05 07:52:51  zvi
! Fixing bug in filling the Jacobian for atmos
!
! Revision 2.7  2002/06/28 11:06:49  zvi
! compute dI/dPtan using chain rule
!
! Revision 2.6  2002/06/19 11:00:36  zvi
! changing from Cspline to InterpolateValues routine
!
! Revision 2.5  2002/05/22 19:43:03  zvi
! Fix a bug in the mol. index loop
!
! Revision 2.4  2002/02/16 06:50:01  zvi
! Some cosmetic code changes
!
! Revision 2.3  2002/02/15 22:51:58  livesey
! Bug fix for case where no ptan derivative wanted
!
! Revision 2.2  2002/01/27 08:37:50  zvi
! Adding Users selected coefficients for derivatives
!
! Revision 2.1  2001/11/08 00:10:49  livesey
! Updated to include extinction
!
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model
!
! Revision 1.21  2001/05/09 19:46:49  vsnyder
! Use new bandHeight argument of createBlock
!
! Revision 1.20  2001/05/03 22:26:29  vsnyder
! Insert copyright notice, some cosmetic changes
!
! Revision 1.19  2001/05/02 20:49:23  zvi
! Cleaning up code
!
! Revision 1.18  2001/05/01 00:42:54  zvi
! Fixing phi window bug
!
! Revision 1.17  2001/04/28 17:47:57  livesey
! Now accepts and sets rowFlags
!
! Revision 1.16  2001/04/27 22:02:13  vsnyder
! Don't compute derivatives if Jacobian isn't present
!
! Revision 1.15  2001/04/27 00:13:29  zvi
! Fixing some phiwindow bug
!
! Revision 1.14  2001/04/26 22:54:41  zvi
! Fixing some phiwindow bug
!
! Revision 1.13  2001/04/24 21:32:45  zvi
! fixing a dimension bug..
!
! Revision 1.12  2001/04/20 23:09:29  livesey
! Now folds in place
!
! Revision 1.11  2001/04/20 02:57:00  livesey
! Writes derivatives in matrix_t
!
! Revision 1.10  2001/04/19 23:56:52  livesey
! New parameters
!
! Revision 1.9  2001/04/10 10:14:16  zvi
! Fixing bug in convolve routines
!
! Revision 1.8  2001/04/10 02:25:14  livesey
! Tidied up some code
!
! Revision 1.7  2001/03/31 23:40:55  zvi
! Eliminate l2pcdim (dimension parameters) move to allocatable ..
!
! Revision 1.6  2001/03/29 02:54:49  livesey
! Removed print statements
!
! Revision 1.5  2001/03/29 02:54:29  livesey
! Changed assumed size to assumed shape
!
! Revision 1.4  2001/03/26 17:56:14  zvi
! New codes to deal with dh_dt_path issue.. now being computed on the fly
!
! Revision 1.3  2001/03/21 01:10:38  livesey
! Now gets Ptan from vector
!
! Revision 1.2  2001/03/07 23:45:15  zvi
! Adding logical flags fro Temp, Atmos & Spect. derivatives
!
! Revision 1.1  2000/06/21 21:56:14  zvi
! First version D.P.
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
