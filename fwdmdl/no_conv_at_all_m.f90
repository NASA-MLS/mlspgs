! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module NO_CONV_AT_ALL_M
  use MLSCommon, only: I4, R4, R8
  use L2PC_PFA_STRUCTURES, only: K_MATRIX_INFO
  USE MLSNumerics, ONLY: INTERPOLATEVALUES
  use dump_0,only:dump
  use VectorsModule, only: Vector_T, VectorValue_T, GETVECTORQUANTITYBYTYPE
  use ForwardModelConfig, only: ForwardModelConfig_T
  use Intrinsic, only: L_VMR
  use String_Table, only: GET_STRING
  use MatrixModule_0, only: M_ABSENT, M_BANDED, M_FULL, DUMP
  use MatrixModule_1, only: CREATEBLOCK, FINDBLOCK, MATRIX_T, DUMP
  use Molecules, only: L_EXTINCTION
  use MLSMessageModule, only: MLSMSG_Error, MLSMessage
  USE Load_sps_data_m, only: Grids_T

  implicit NONE
  private

  public :: NO_CONV_AT_ALL

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
CONTAINS
  !-------------------------------------------------------------------------
  ! This subroutine transfers the derivatives over from the internal
  ! convolution grid to the users specified points. This module uses
  ! cubic spline interpolation to do the job.

  Subroutine no_conv_at_all ( ForwardModelConfig, ForwardModelIn, maf, &
           & Channel, WindowStart, WindowFinish, Temp, Ptan, Radiance, &
           & t_deriv_flag,ptg_angles,chi_out,dhdz_out,dx_dh_out,Grids_f,&
           & I_raw,sbRatio,mol_cat_indx, rowFlags, Jacobian, di_dt, di_df )

    Type (ForwardModelConfig_T) :: FORWARDMODELCONFIG
    Type (Vector_T), intent(in) :: FORWARDMODELIN

    Integer, INTENT(IN) :: maf
    Integer, INTENT(IN) :: CHANNEL
    Integer, INTENT(IN) :: WINDOWSTART
    Integer, INTENT(IN) :: WINDOWFINISH
    Integer, INTENT(IN) :: mol_cat_indx(:)

    Type (VectorValue_T), INTENT(IN) :: TEMP
    Type (VectorValue_T), INTENT(IN) :: PTAN
    Type (VectorValue_T), INTENT(INOUT) :: RADIANCE

    Type (Grids_T), INTENT(IN) :: Grids_f

    Real(r8), INTENT(IN) :: sbRatio
    Real(r8), INTENT(IN) :: i_raw(:),ptg_angles(:),chi_out(:),dhdz_out(:), &
                         &  dx_dh_out(:)
!
! derivative of radiance w.r.t. temperature on chi_in
    Real(r8), OPTIONAL, INTENT(IN) :: di_dt(:,:)

! mixing ratio derivatives or any parameter which behaves like VMR
    Real(r8), OPTIONAL, INTENT(IN) :: di_df(:,:)

    Type (Matrix_T), INTENT(INOUT), OPTIONAL :: Jacobian
!
    Logical, DIMENSION(:), pointer :: t_deriv_flag

    Logical, DIMENSION(:), INTENT(INOUT) :: rowFlags ! Flag to calling code

    ! -----     Local Variables     ------------------------------------

    Integer:: No_t, No_tan_hts

    Type (VectorValue_T), pointer :: F  ! VMR quantity

    Integer :: jf, jz, no_mol, l
    Integer :: is, j, k, nf, sv_f, sv_t_len
    Integer :: Row, col                     ! Matrix row & column indices
    Integer :: ptg_i,noPtan,noChans,Ind     ! Indices

    Real(r8) :: Rad( size(ptg_angles)), q
    Real(r8) :: SRad(ptan%template%noSurfs)
    Real(r8) :: di_dx(ptan%template%noSurfs)
    Real(r8) :: I_star_all(ptan%template%noSurfs)

    ! -----  Begin the code  -------------------------------------------

    no_t = temp%template%noSurfs
    no_tan_hts = size(ptg_angles)

    noPtan = ptan%template%noSurfs
    noChans = radiance%template%noChans

! Ptan derivative

    col = 0
    if ( PRESENT (Jacobian) ) &
      & col = FindBlock ( Jacobian%col, ptan%index, maf )

    if ( col > 0 ) then

      Call InterpolateValues(ptg_angles, i_raw, chi_out, i_star_all, &
                           & METHOD='S',dyByDx=di_dx)

! Use the chain rule to compute dI/dPtan on the output grid:

      SRad(1:noPtan) = di_dx(1:noPtan) * dx_dh_out(1:noPtan)
      SRad(1:noPtan) = SRad(1:noPtan) * dhdz_out(1:noPtan)

      row = FindBlock ( Jacobian%row, radiance%index, maf )
      rowFlags(row) = .TRUE.

      select case ( jacobian%block(Row,col)%kind )
        case ( m_absent )
          call CreateBlock ( Jacobian, row, col, m_banded, &
            & radiance%template%noSurfs*noChans,bandHeight=noChans )
          jacobian%block(row,col)%values = 0.0_r8
        case ( m_banded )
        case default
          call MLSMessage ( MLSMSG_Error, ModuleName,&
            & 'Wrong matrix type for ptan derivative' )
      end select

      do ptg_i = 1, noPtan
        ind = channel + noChans*(ptg_i-1)
        jacobian%block(row,col)%values(ind,1) = &
           & jacobian%block(row,col)%values(ind,1) + sbRatio * SRad(ptg_i)
        jacobian%block(row,col)%r1(ptg_i) = 1 + noChans * (ptg_i - 1)
        jacobian%block(row,col)%r2(ptg_i) = noChans * ptg_i
      end do

    else

      Call InterpolateValues(ptg_angles,i_raw,chi_out,i_star_all,METHOD='S')

    endif

    do ptg_i = 1, noPtan
      ind = channel + noChans*(ptg_i-1)
      radiance%values( ind, maf ) = &
             & radiance%values ( ind, maf ) + sbRatio * i_star_all(ptg_i)
    end do

    if ( .not. PRESENT(Jacobian) ) Return

    if ( .not. ANY((/forwardModelConfig%temp_der, &
                   & forwardModelConfig%atmos_der, &
                   & forwardModelConfig%spect_der/)) ) Return

    ! Now transfer the other fwd_mdl derivatives to the output pointing
    ! values

    ! ********************* Temperature derivatives ******************

    ! check to determine if derivative is desired for this parameter

    if ( forwardModelConfig%temp_der ) then

    ! Derivatives needed continue to process

      sv_t_len = 0
      Rad(1:) = 0.0
      SRad(1:) = 0.0
      k = no_tan_hts

      do nf = windowStart, windowFinish

        col = FindBlock ( Jacobian%col, temp%index, nf )
        select case ( Jacobian%block(row,col)%kind )
        case ( m_absent )
          call CreateBlock ( Jacobian, row, col, m_full )
          jacobian%block(row,col)%values = 0.0_r8
        case ( m_full )
        case default
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Wrong type for temperature derivative matrix' )
        end select

        do jz = 1, no_t

! Check if derivatives are needed for this (zeta & phi) :

          sv_t_len = sv_t_len + 1
          if(.NOT. t_deriv_flag(sv_t_len)) CYCLE

          Rad(1:k) = di_dt(1:k,sv_t_len)
          Call InterpolateValues ( ptg_angles, Rad, chi_out, Srad, 'S')
          do ptg_i = 1, noPtan
            ind = channel + noChans*(ptg_i-1)
            q = jacobian%block(row,col)%values(ind,jz)
            jacobian%block(row,col)%values(ind,jz) = q + sbRatio*Srad(ptg_i)
          end do
        end do

      end do

    end if

    if ( forwardModelConfig%atmos_der ) then

      ! ****************** atmospheric derivatives ******************

      sv_f = 0
      no_mol = size(mol_cat_indx)

      do is = 1, no_mol

        jz = mol_cat_indx(is)
        l = forwardModelConfig%molecules(jz)
        if ( l == l_extinction ) then
          f => GetVectorQuantityByType(forwardModelIn, &
                & quantityType=l_extinction,radiometer = &
                & radiance%template%radiometer, noError=.true. )
        else
          f => GetVectorQuantityByType ( forwardModelIn, quantityType=l_vmr,&
                & molecule=l, noError=.true. )
        endif

        if(.not. associated(f) ) then
          jf = Grids_f%windowfinish(is)-Grids_f%windowStart(is)+1
          k = Grids_f%no_f(is) * Grids_f%no_z(is)
          sv_f = sv_f + jf * k
          CYCLE
        endif
!
        DO jf = Grids_f%windowStart(is), Grids_f%windowfinish(is)
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

          DO k = 1, Grids_f%no_f(is) * Grids_f%no_z(is)

! Check if derivatives are needed for this (zeta & phi) :

            sv_f = sv_f + 1
            if(.NOT. Grids_f%deriv_flags(sv_f)) CYCLE

            Rad(1:no_tan_hts) = di_df(1:no_tan_hts,sv_f)
            Call InterpolateValues (ptg_angles, Rad, chi_out, Srad, 'L')
            do ptg_i = 1, noPtan
              ind = channel + noChans*(ptg_i-1)
              q = jacobian%block(row,col)%values(ind,sv_f)
              jacobian%block(row,col)%values(ind,sv_f) = &
                                               &  q + sbRatio*Srad(ptg_i)
            end do

          end do
!
        end do
!
      end do

    endif

    Return

  End Subroutine NO_CONV_AT_ALL

END module NO_CONV_AT_ALL_M
! $Log$
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
