! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module CONVOLVE_ALL_M
  use AntennaPatterns_m, only: AntennaPattern_T
  use Load_sps_data_m, only: Grids_T
  use DCSPLINE_DER_M, only: CSPLINE_DER
  use DUMP_0, only: DUMP
  use D_LINTRP_M, only: LINTRP
  use D_CSPLINE_M, only: CSPLINE
  use FOV_CONVOLVE_M, only: FOV_CONVOLVE
  use HYDROSTATIC_INTRP, only: GET_PRESSURES
  use L2PC_PFA_STRUCTURES, only: K_MATRIX_INFO
  use VectorsModule, only: Vector_T, VectorValue_T, GetVectorQuantityByType
  use ForwardModelConfig, only: ForwardModelConfig_T
  use MLSCommon, only: I4, R4, R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use Molecules, only: L_EXTINCTION
  use Intrinsic, only: L_VMR
  use String_table, only: GET_STRING
  use MatrixModule_0, only: M_ABSENT, M_BANDED, M_FULL, DUMP
  use MatrixModule_1, only: CREATEBLOCK, FINDBLOCK, MATRIX_T, DUMP

  implicit NONE
  private
  public :: CONVOLVE_ALL
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
  !---------------------------------------------------------------------------
  ! This subroutine transfers the derivatives over from the internal
  ! convolution grid to the users specified points. This module uses
  ! cubic spline interpolation to do the job.

  Subroutine convolve_all (ForwardModelConfig, ForwardModelIn, maf, channel, &
  & windowStart, windowFinish, mafTInstance, temp, ptan, radiance, tan_press,&
  & ptg_angles,tan_temp,dx_dt,d2x_dxdt, si,center_angle,i_raw, k_temp,       &
  & k_atmos, sbRatio, t_deriv_flag, Grids_f, Jacobian, rowFlags,             &
  & AntennaPattern, mol_cat_indx, Ier, lu_debug)

    ! Dummy arguments
    type (ForwardModelConfig_T), intent(in) :: FORWARDMODELCONFIG
    type (Vector_t), intent(in) :: FORWARDMODELIN
    integer, intent(in) :: MAF
    integer, intent(in) :: CHANNEL
    integer, intent(in) :: WINDOWSTART
    integer, intent(in) :: WINDOWFINISH
    integer, intent(in) :: MAFTINSTANCE
    integer, intent(IN) :: mol_cat_indx(:)
    type (VectorValue_T), intent(in) :: TEMP
    type (VectorValue_T), intent(in) :: PTAN
    type (VectorValue_T), intent(inout) :: RADIANCE
    real(r8), intent(IN) :: TAN_PRESS(:), PTG_ANGLES(:), TAN_TEMP(:)
    real(r8), intent(IN) :: DX_DT(:,:), D2X_DXDT(:,:)
    integer, intent(IN) :: si
    real(r8), intent(IN) :: CENTER_ANGLE
    real(r8), intent(IN) :: I_RAW(:)
    Real(r4), intent(in) :: k_temp(:,:,WindowStart:)
    Real(r4), intent(in) :: k_atmos(:,:,:,WindowStart:,:)
    real(r8), intent(in) :: SBRATIO
    logical, dimension(:), pointer :: t_deriv_flag
    type (Grids_T), intent(in) :: Grids_f
    type (Matrix_T), intent(inout), optional :: Jacobian
    logical, dimension(:), intent(inout) :: rowFlags ! Flag to calling code
    type(antennaPattern_T), intent(in) :: AntennaPattern
    integer, intent(out) :: ier         ! Flag

    ! -----     Local Variables     ----------------------------------------

    type (VectorValue_t), pointer :: f

    integer :: FFT_INDEX(size(antennaPattern%aaap))
    integer :: Ind, ht_ind(1),kp,kz,kf  ! Indecies and various sizes
    integer :: Row, Col                 ! Matrix entries
    integer :: FFT_pts, Is, J, Ktr, Nf, Ntr, Ptg_i, Sv_i
    integer :: No_tan_hts, Lk, Uk, no_mol, jf, jz, l, sv_t, sv_f

    logical :: Want_Deriv

    real(r8) :: K_star_tmp(ptan%template%noSurfs)
    Real(r8) :: Q, R
    Real(r8) :: SRad(ptan%template%noSurfs), Term(ptan%template%noSurfs)
    Real(r8), dimension(size(fft_index)) :: FFT_PRESS, FFT_ANGLES, RAD
! bills debug
    INTEGER(i4), intent(in) :: lu_debug

 ! -----  Begin the code  -----------------------------------------

    k_star_tmp = 0.0
    no_tan_hts = size(tan_press)

    ! Compute the ratio of the strengths

    ! This subroutine is called by channel

    Ier = 0
    Ntr = Size(antennaPattern%aaap)

    Rad = 0.0
    Rad(1:no_tan_hts) = i_raw(1:no_tan_hts)

    j = ptan%template%noSurfs

    ! Compute the convolution of the mixed radiances

    fft_pts = nint(Log(Real(Ntr))/log(2.0))

    fft_angles = 0.0
    fft_angles(1:no_tan_hts) = ptg_angles(1:no_tan_hts)

    Call fov_convolve ( fft_angles, Rad, center_angle, 1, no_tan_hts, &
      &                 fft_pts, AntennaPattern, Ier )
    if ( Ier /= 0) Return

    !  Get 'Ntr' pressures associated with the fft_angles:

    Call get_pressures ( 'a', ptg_angles, tan_temp, tan_press, no_tan_hts, &
      &                   fft_angles, fft_press, Ntr, Ier )
! bills debug
    WRITE(lu_debug,'(a)') 'fft_angles, fft_press, rad'
    DO  is = 1 , ntr 
      WRITE(lu_debug,'(f11.8,1x,f10.5,1x,f9.4)') fft_angles(is), &
           fft_press(is), rad(is)
    ENDDO
    if ( Ier /= 0) Return

    ! Make sure the fft_press array is MONOTONICALY increasing:

    is = 1
    fft_index = 0
    do while (is < Ntr-1  .and.  fft_press(is) >= fft_press(is+1))
      is = is + 1
    end do

    Ktr = 1
    Rad(Ktr) = Rad(is)
    fft_index(Ktr) = is
    fft_press(Ktr) = fft_press(is)

    do ptg_i = is+1, Ntr
      q = fft_press(ptg_i)
      if ( q > fft_press(Ktr)) then
        Ktr = Ktr + 1
        fft_press(Ktr) = q
        Rad(Ktr) = Rad(ptg_i)
        fft_index(Ktr) = ptg_i
      end if
    end do

    if(Ktr == Ntr) fft_index(1) = -2

! Interpolate the output values
! (Store the radiances derivative with respect to pointing pressures in: Term)

    want_deriv = present(jacobian) .and. any( (/ &
      & forwardModelConfig%temp_der, &
      & forwardModelConfig%atmos_der,&
      & forwardModelConfig%spect_der/) )

!??? Can we use cspline instead of cspline_der if .not. want_deriv ???

    Call Cspline_der(fft_press,Ptan%values(:,maf),Rad,SRad,Term,Ktr,j) 
    do ptg_i = 1, j
      ind = channel + radiance%template%noChans*(ptg_i-1)
      radiance%values( ind, maf ) = radiance%values ( ind, maf ) + &
        & sbRatio*SRad(ptg_i)
    end do

    ! Find out if user wants pointing derivatives
    if ( .not. want_deriv ) return

    ! Derivatives wanted,find index location Jacobian and write the derivative
    row = FindBlock( Jacobian%row, radiance%index, maf )
    rowFlags(row) = .true.
    col = FindBlock ( Jacobian%col, ptan%index, maf )

    ! Of course, we might not care about ptan
    if ( col > 0 ) then
      select case (jacobian%block(Row,col)%kind)
      case (m_absent)
        call CreateBlock ( Jacobian, row, col, m_banded, &
          & radiance%template%noSurfs*radiance%template%noChans, &
          & bandHeight=radiance%template%noChans )
        jacobian%block(row,col)%values = 0.0_r8
      case (m_banded)
      case default
        call MLSMessage ( MLSMSG_Error, ModuleName,&
          & 'Wrong matrix type for ptan derivative')
      end select
      do ptg_i = 1, j
        ind = channel + radiance%template%noChans*(ptg_i-1)
        jacobian%block(row,col)%values( ind, 1 ) = &
          & jacobian%block(row,col)%values( ind, 1 ) + sbRatio*term(ptg_i)
        jacobian%block(row,col)%r1(ptg_i) = &
          & 1 + radiance%template%noChans*(ptg_i-1)
        jacobian%block(row,col)%r2(ptg_i) = &
          & radiance%template%noChans*ptg_i
      end do
    end if

    ! Now transfer the other fwd_mdl derivatives to the output pointing
    ! values

    ! ********************* Temperature derivatives ******************

    ! check to determine if derivative is desired for this parameter

    if ( forwardModelConfig%temp_der) then

    ! Derivatives needed continue to process

      sv_t = 0
      ht_ind = 0

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

        do jz = 1, temp%template%noSurfs

! Check if derivatives are needed for this (zeta & phi) :

          sv_t = sv_t + 1
          if(.NOT. t_deriv_flag(sv_t)) CYCLE

          ! run through representation basis coefficients
          ! Integrand over temperature derivative plus pointing differential

!  *** Bill's correction
!
          IF ( ANY(d2x_dxdt(:,jz) > 0.0)) THEN
!
            ht_ind = MAXLOC(d2x_dxdt(:,jz)) 

            do ptg_i = 1, no_tan_hts
              q = d2x_dxdt(ptg_i,jz)
              Rad(ptg_i) = (i_raw(ptg_i) - i_raw(ht_ind(1))) * q + &
                         &  k_temp(ptg_i,jz,nf)
            end do

          ! Now, Convolve:

            fft_angles(1:no_tan_hts) = ptg_angles(1:no_tan_hts)
            Call fov_convolve ( fft_angles, Rad, center_angle, 1, &
                             &  no_tan_hts, fft_pts, AntennaPattern, Ier )
            if ( Ier /= 0) Return

            if ( fft_index(1).gt.0) then
              do ptg_i = 1, Ktr
                Rad(ptg_i) = Rad(fft_index(ptg_i))
              end do
            end if

            Call Cspline ( fft_press, Ptan%values(:,maf), Rad, SRad, Ktr, j )
            k_star_tmp(1:j) = SRad(1:j)

          ! Now the convolution of radiance with the derivative antenna field

            Rad(1:no_tan_hts) = i_raw(1:no_tan_hts)

          ! Now, Convolve:

            fft_angles(1:no_tan_hts) = ptg_angles(1:no_tan_hts)
            Call fov_convolve ( fft_angles, Rad, center_angle, 2, &
                             &  no_tan_hts, fft_pts, AntennaPattern, Ier )
            if ( Ier /= 0) Return

            if ( fft_index(1).gt.0) then
              do ptg_i = 1, Ktr
                Rad(ptg_i) = Rad(fft_index(ptg_i))
              end do
            end if

            Call Cspline ( fft_press, Ptan%values(:,maf), Rad, Term, Ktr, j )

          ! Transfer dx_dt from convolution grid onto the output grid

            Call Lintrp (tan_press, Ptan%values(:,maf), dx_dt(1:,jz), SRad, &
                       & no_tan_hts, j )

            k_star_tmp = k_star_tmp + srad*term

          ! the convolution of the radiance weighted hydrostatic derivative
          ! with the antenna derivative

            Rad(1:no_tan_hts) = dx_dt(1:no_tan_hts,jz) * &
                            & (i_raw(1:no_tan_hts) - i_raw(ht_ind(1)))

          ! Now, convolve:

            fft_angles(1:no_tan_hts) = ptg_angles(1:no_tan_hts)
            Call fov_convolve ( fft_angles, Rad, center_angle, 2, &
                             &  no_tan_hts, fft_pts, AntennaPattern, Ier )
            if ( Ier /= 0) Return

            if ( fft_index(1).gt.0) then
              do ptg_i = 1, Ktr
                Rad(ptg_i) = Rad(fft_index(ptg_i))
              end do
            end if

            Call Cspline(fft_press,Ptan%values(:,maf),Rad,Term,Ktr,j)

            do ptg_i = 1, j
              q = k_star_tmp(ptg_i)
              ind = channel + radiance%template%noChans*(ptg_i-1)
              r = jacobian%block(row,col)%values( ind, jz)
              jacobian%block(row,col)%values( ind, jz) = r + &
                                          &    sbRatio*(q - Term(ptg_i))
            end do
!
          ELSE          ! On IF(ANY(d2x_dxdt(:,jz) > 0.0)) THEN
!
! If All the d2x_dxdt(:,jz) = 0.0 then it is a simple cycle, like VMR:
!
! run through the Temp. representation basis coefficients
!
            Rad(1:no_tan_hts) = k_temp(1:no_tan_hts,jz,nf)

! Now Convolve the derivative

            fft_angles(1:no_tan_hts) = ptg_angles(1:no_tan_hts)
            Call fov_convolve ( fft_angles, Rad, center_angle, 1, &
                             &  no_tan_hts, fft_pts, AntennaPattern, Ier )
            if ( Ier /= 0) Return

            if ( fft_index(1).gt.0) then
              do ptg_i = 1, Ktr 
                Rad(ptg_i) = Rad(fft_index(ptg_i))
              end do
            end if

            Call Lintrp(fft_press,Ptan%values(:,maf),Rad,SRad,Ktr,j)

            do ptg_i = 1, j
              ind = channel + radiance%template%noChans*(ptg_i-1)
              q = jacobian%block(row,col)%values( ind, jz )
              jacobian%block(row,col)%values( ind, jz ) = &
                                         &    q + sbRatio*Srad(ptg_i)
            end do

          ENDIF    ! On IF(ANY(d2x_dxdt(:,jz) > 0.0)) THEN

        end do

      end do

    end if

    ! ****************** atmospheric derivatives ******************

    if ( forwardModelConfig%atmos_der) then

      sv_f = 0
      no_mol = size(mol_cat_indx)

      lk = lbound(k_atmos,4)   ! The lower Phi dimension
      uk = ubound(k_atmos,4)   ! The upper Phi dimension

      do is = 1, no_mol

        kp = Grids_f%no_p(is)
        kz = Grids_f%no_z(is)
        kf = Grids_f%no_f(is)

        jz = mol_cat_indx(is)
        l = forwardModelConfig%molecules(jz)
        if ( l == l_extinction ) then
          f => GetVectorQuantityByType(forwardModelIn, &
            &                          quantityType=l_extinction, &
            &  radiometer=radiance%template%radiometer, noError=.true. )
        else
          f => GetVectorQuantityByType ( forwardModelIn, quantityType=l_vmr, &
            &  molecule=l, noError=.true. )
        endif

        if ( associated(f)) then

          do nf = 1, kp

          ! run through phi representation basis coefficients

            if ( nf+lk-1 > uk) then
              sv_f = sv_f + kf * kz * (kp-nf+1)
              EXIT
            endif

            col = FindBlock ( Jacobian%col, f%index, nf+windowStart-1 )
            select case ( Jacobian%block(row,col)%kind )
            case ( m_absent )
              call CreateBlock ( Jacobian, row, col, m_full )
              jacobian%block(row,col)%values = 0.0_r8
            case ( m_full )
            case default
              call MLSMessage ( MLSMSG_Error, ModuleName, &
                & 'Wrong type for vmr derivative matrix' )
            end select

            sv_i = 0
            do jz = 1, kz

            ! run through zeta representation basis coefficients

              do jf = 1, kf

! Check if derivatives are needed for this molecule (zeta, phi & channel) :

                sv_f = sv_f + 1
                IF(.NOT. Grids_f%deriv_flags(sv_f)) CYCLE

                ! run through Frequencies basis coefficients

                sv_i = sv_i + 1
                Rad(1:no_tan_hts) = k_atmos(1:no_tan_hts,jf,jz,nf+lk-1,is)

                ! Now Convolve the derivative

                fft_angles(1:no_tan_hts) = ptg_angles(1:no_tan_hts)
                Call fov_convolve(fft_angles,Rad,center_angle,1,no_tan_hts, &
                  &               fft_pts,AntennaPattern,Ier)
                if ( Ier /= 0) Return

                if ( fft_index(1).gt.0) then
                  do ptg_i = 1, Ktr
                    Rad(ptg_i) = Rad(fft_index(ptg_i))
                  end do
                end if

                ! Interpolate onto the output grid, and store in Jacobian ..

                Call Lintrp(fft_press,Ptan%values(:,maf),Rad,SRad,Ktr,j)

                do ptg_i = 1, j
                  ind = channel + radiance%template%noChans*(ptg_i-1)
                  q = jacobian%block(row,col)%values( ind, sv_i )
                  jacobian%block(row,col)%values( ind, sv_i ) = &
                                            &    q + sbRatio*Srad(ptg_i)
                end do                  ! Pointing
!
              end do                    ! F channels (Frequencies)
!
            end do                      ! F surfs
!
          end do                        ! F instances
!
        else                            ! On: if ( associated(f)) ...
!
          sv_f = sv_f + kp * kz * kf
!
        end if                          ! Want derivatives for this species
!
      end do                            ! Loop over species
!
    end if                              ! Any derivatives
!
10  CONTINUE

    Return
!
  End Subroutine CONVOLVE_ALL
!
end module CONVOLVE_ALL_M
! $Log$
! Revision 2.7  2002/06/04 10:28:02  zvi
! Encorporate deriv. flag into convolution, fixing a bug with species ruuning index
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
