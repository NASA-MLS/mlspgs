module CONVOLVE_ALL_M
  use AntennaPatterns_m, only: AntennaPattern_T
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
    windowStart, windowFinish, mafTInstance, temp, ptan, radiance, &
    tan_press,ptg_angles,tan_temp,dx_dt,d2x_dxdt, si,center_angle,i_raw, &
    k_temp, k_atmos, sbRatio, Jacobian, rowFlags, AntennaPattern,Ier)

    ! Dummy arguments
    type (ForwardModelConfig_T), intent(in) :: FORWARDMODELCONFIG
    type (Vector_t), intent(in) :: FORWARDMODELIN
    integer, intent(in) :: MAF
    integer, intent(in) :: CHANNEL
    integer, intent(in) :: WINDOWSTART
    integer, intent(in) :: WINDOWFINISH
    integer, intent(in) :: MAFTINSTANCE
    type (VectorValue_T), intent(in) :: TEMP
    type (VectorValue_T), intent(in) :: PTAN
    type (VectorValue_T), intent(inout) :: RADIANCE
    type (Matrix_T), intent(inout), optional :: Jacobian

    integer, intent(IN) :: si
    real(r8), intent(IN) :: CENTER_ANGLE
    real(r8), intent(IN) :: I_RAW(:)
    real(r8), intent(IN) :: TAN_PRESS(:), PTG_ANGLES(:), TAN_TEMP(:)
    real(r8), intent(IN) :: DX_DT(:,:), D2X_DXDT(:,:)
    real(r8), intent(in) :: SBRATIO
    type(antennaPattern_T), intent(in) :: AntennaPattern
!
    Real(r4), intent(in) :: k_temp(:,:,WindowStart:)   ! (Nptg,mxco,mnp)
    Real(r4), intent(in) :: k_atmos(:,:,WindowStart:,:) ! (Nptg,mxco,mnp,Nsps)

    real(r8) :: k_star_tmp(ptan%template%noSurfs)
    logical, dimension(:), intent(inout) :: rowFlags ! Flag to calling code

    integer, intent(out) :: ier         ! Flag

    ! -----     Local Variables     ----------------------------------------

    type (VectorValue_t), pointer :: f

    integer :: FFT_INDEX(size(antennaPattern%aaap))
    integer :: Ind                      ! Index
    integer :: N, I, J, Is, Ktr, Nf, Ntr, Ptg_i, Sv_i, Spectag, Ki, Kc, FFT_pts
    integer :: Row, Col                 ! Matrix entries
    integer :: No_t, No_tan_hts, No_phi_t, Lk, Uk
    logical :: Want_Deriv

    Real(r8) :: Q, R
    Real(r8) :: SRad(ptan%template%noSurfs), Term(ptan%template%noSurfs)
    Real(r8), dimension(size(fft_index)) :: FFT_PRESS, FFT_ANGLES, RAD

    Character :: CA

    ! -----  Begin the code  -----------------------------------------
    no_t = temp%template%noSurfs
    no_phi_t = temp%template%noInstances
    no_tan_hts = size(tan_press)
    k_star_tmp = 0.0

    ! Compute the ratio of the strengths

    ! This subroutine is called by channel

    Ier = 0
    ntr = size(antennaPattern%aaap)

    Rad(1:no_tan_hts) = i_raw(1:no_tan_hts)

    kc = 0
    ki = 0
    j = ptan%template%noSurfs

    ! Compute the convolution of the mixed radiances

    fft_pts = nint(log(real(size(AntennaPattern%aaap)))/log(2.0))
    fft_angles(1:size(tan_press)) = ptg_angles(1:size(tan_press))
    Call fov_convolve ( fft_angles, Rad,center_angle, 1, no_tan_hts, &
      &                 fft_pts, AntennaPattern, Ier )
    if ( Ier /= 0) Return

    !  Get 'Ntr' pressures associated with the fft_angles:

    Call get_pressures ( 'a', ptg_angles, tan_temp, tan_press, no_tan_hts, &
      &                   fft_angles, fft_press, Ntr, Ier )
    if ( Ier /= 0) Return

    ! Make sure the fft_press array is MONOTONICALY increasing:

    is = 1
    do while (is < Ntr  .and.  fft_press(is) >= fft_press(is+1))
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

    ! Interpolate the output values
    ! (Store the radiances derivative with respect to pointing pressures in: Term)

    want_deriv = present(jacobian) .and. any( (/ &
      & forwardModelConfig%temp_der, &
      & forwardModelConfig%atmos_der,&
      & forwardModelConfig%spect_der/) )
!??? Can we use cspline instead of cspline_der if .not. want_deriv ???
    Call Cspline_der ( fft_press, Ptan%values(:,maf), Rad, SRad, Term, Ktr, j )
    do ptg_i = 1, j
      ind = channel + radiance%template%noChans*(ptg_i-1)
      radiance%values( ind, maf ) = radiance%values ( ind, maf ) + &
        & sbRatio*SRad(ptg_i)
    end do

    ! Find out if user wants pointing derivatives
    if ( .not. want_deriv ) return
    ! Derivatives wanted,find index location k_star_all and write the derivative
    row = FindBlock( Jacobian%row, radiance%index, maf )
    rowFlags(row) = .true.
    col = FindBlock ( Jacobian%col, ptan%index, maf )
    select case (jacobian%block(Row,col)%kind)
    case (m_absent)
      call CreateBlock ( Jacobian, row, col, m_banded, &
        &    radiance%template%noSurfs*radiance%template%noChans )
      jacobian%block(row,col)%values = 0.0_r8
      do ptg_i = 1, j
        jacobian%block(row,col)%r1(ptg_i) = &
          & 1 + radiance%template%noChans*(ptg_i-1)
        jacobian%block(row,col)%r2(ptg_i) = &
          & radiance%template%noChans*ptg_i
      end do
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

    ! Now transfer the other fwd_mdl derivatives to the output pointing
    ! values

    ! ********************* Temperature derivatives ******************

    ! check to determine if derivative is desired for this parameter

    if ( forwardModelConfig%temp_der) then

      ! Derivatives needed continue to process

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
        do sv_i = 1, temp%template%noSurfs

          ! run through representation basis coefficients
          ! Integrand over temperature derivative plus pointing differential

          do ptg_i = 1, no_tan_hts
            q = 0.0
            if ( nf == mafTInstance) q = d2x_dxdt(ptg_i,sv_i)
            Rad(ptg_i) = i_raw(ptg_i) * q + k_temp(ptg_i,sv_i,nf)
          end do

          ! Now, Convolve:

          fft_angles(1:no_tan_hts) = ptg_angles(1:no_tan_hts)
          Call fov_convolve ( fft_angles, Rad, center_angle, 1, no_tan_hts, &
            &                 fft_pts, AntennaPattern, Ier )
          if ( Ier /= 0) Return

          if ( fft_index(1).gt.0) then
            do ptg_i = 1, Ktr
              Rad(ptg_i) = Rad(fft_index(ptg_i))
            end do
          end if

          Call Cspline ( fft_press, Ptan%values(:,maf), Rad, SRad, Ktr, j )
          k_star_tmp(1:j) = SRad(1:j)

          !  For any index off center Phi, skip the rest of the phi loop ...

          if ( nf /= mafTInstance) then
            do ptg_i = 1, j
              ind = channel + radiance%template%noChans*(ptg_i-1)
              jacobian%block(row,col)%values( ind, sv_i ) = &
                & jacobian%block(row,col)%values( ind, sv_i ) + &
                &   sbRatio*k_star_tmp(ptg_i)
            end do
            cycle
          end if

          ! Now the convolution of radiance with the derivative antenna field

          Rad(1:no_tan_hts) = i_raw(1:no_tan_hts)

          ! Now, Convolve:

          fft_angles(1:no_tan_hts) = ptg_angles(1:no_tan_hts)
          Call fov_convolve ( fft_angles, Rad, center_angle, 2, no_tan_hts, &
            &                 fft_pts, AntennaPattern, Ier )
          if ( Ier /= 0) Return

          if ( fft_index(1).gt.0) then
            do ptg_i = 1, Ktr
              Rad(ptg_i) = Rad(fft_index(ptg_i))
            end do
          end if

          Call Cspline ( fft_press, Ptan%values(:,maf), Rad, Term, Ktr, j )

          ! Transfer dx_dt from convolution grid onto the output grid

          Call Lintrp (tan_press, Ptan%values(:,maf), dx_dt(1:,sv_i), SRad, &
            & no_tan_hts, j )

          k_star_tmp = k_star_tmp + srad*term

          ! the convolution of the radiance weighted hydrostatic derivative
          ! with the antenna derivative

          Rad(1:no_tan_hts) = &
            dx_dt(1:no_tan_hts,sv_i) * i_raw(1:no_tan_hts)

          ! Now, convolve:

          fft_angles(1:no_tan_hts) = ptg_angles(1:no_tan_hts)
          Call fov_convolve ( fft_angles, Rad, center_angle, 2, no_tan_hts, &
            &                 fft_pts, AntennaPattern, Ier )
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
            jacobian%block(row,col)%values( ind, sv_i ) = &
              jacobian%block(row,col)%values( ind, sv_i ) + sbRatio*(q - Term(ptg_i))
          end do

        end do

      end do

    end if

    ! ****************** atmospheric derivatives ******************

    if ( forwardModelConfig%atmos_der) then

      lk = lbound(k_atmos,3)
      uk = ubound(k_atmos,3)
      do is = 1, size(ForwardModelConfig%molecules) ! What about derivatives!???NJL

        f => GetVectorQuantityByType ( forwardModelIn, quantityType=l_vmr, &
          & molecule=forwardModelConfig%molecules(is), noError=.true. )

        if ( associated(f)) then

          ! Derivatives needed continue to process

          do nf = 1, f%template%noInstances
            if ( nf+lk-1 > uk) EXIT
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

            do sv_i = 1, f%template%noSurfs

              ! run through representation basis coefficients

              Rad(1:no_tan_hts) = k_atmos(1:no_tan_hts,sv_i,nf+lk-1,is)

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

              ! Interpolate onto the output grid, and store in k_star_all ..

              Call Lintrp(fft_press,Ptan%values(:,maf),Rad,SRad,Ktr,j)
              do ptg_i = 1,j
                ind = channel+ radiance%template%noChans*(ptg_i-1)
                jacobian%block(row,col)%values( ind, sv_i ) = &
                  & jacobian%block(row,col)%values( ind, sv_i ) + sbRatio*Srad(ptg_i)
              end do                    ! Pointing
            end do                      ! F surfs
          end do                        ! F instances
        end if                          ! Want derivatives for this species
      end do                            ! Loop over species
    end if                              ! Any derivatives
!
    !     ! ****************** Spectroscopic derivatives ******************
    ! !
    !     if ( spect_der) then
    ! !
    !       do is = 1, n_sps
    ! !
    !         i = spect_atmos(is)
    !         if ( i < 1) CYCLE
    !         if ( .not.spectroscopic(i)%DER_CALC(band)) CYCLE
    ! !
    !         ! Derivatives needed continue to process
    ! !
    !         Spectag = atmospheric(is)%spectag
    ! !
    !         DO
    ! !
    !           if ( spectroscopic(i)%Spectag /= Spectag) EXIT
    !           n = spectroscopic(i)%no_phi_values
    !           nz = spectroscopic(i)%no_zeta_values
    !           CA = spectroscopic(i)%type
    !           ki = ki + 1
    !           kc = kc + 1
    !           k_star_info(kc)%name = spectroscopic(i)%NAME
    !           k_star_info(kc)%first_dim_index = ki
    !           k_star_info(kc)%no_phi_basis = n
    !           k_star_info(kc)%no_zeta_basis = nz
    !           k_star_info(kc)%zeta_basis(1:nz) = &
    !             &  spectroscopic(i)%zeta_basis(1:nz)
    ! !
    !           do nf = 1, n
    ! !
    !             do sv_i = 1, nz
    ! !
    !               select case ( CA )
    !               case ( 'W' )
    !                 Rad(1:no_tan_hts) = k_spect_dw(1:no_tan_hts,sv_i,nf,i)
    !               case ( 'N' )
    !                 Rad(1:no_tan_hts) = k_spect_dn(1:no_tan_hts,sv_i,nf,i)
    !               case ( 'V' )
    !                 Rad(1:no_tan_hts) = k_spect_dnu(1:no_tan_hts,sv_i,nf,i)
    !               case default
    !                 Ier = -99
    !                 Print *,'** Unknown Spectroscopic element !'
    !                 Return
    !               end select
    ! !
    !               ! Now Convolve the derivative
    ! !
    !               fft_angles(1:no_tan_hts) = ptg_angles(1:no_tan_hts)
    !               Call fov_convolve(fft_angles,Rad,center_angle,1,no_tan_hts, &
    !                 &               fft_pts,AntennaPattern,Ier)
    !               if ( Ier /= 0) Return
    ! !
    !               if ( fft_index(1).gt.0) then
    !                 do ptg_i = 1, Ktr
    !                   Rad(ptg_i) = Rad(fft_index(ptg_i))
    !                 end do
    !               end if
    ! !
    !               ! Interpolate onto the output grid, and store in k_star_all ..
    ! !
    !               Call Lintrp(fft_press,Ptan,Rad,SRad,Ktr,j)
    !               k_star_all(ki,sv_i,nf,1:j) = SRad(1:j)
    ! !
    !             end do        ! sv_i loop
    ! !
    !           end do          ! nf loop
    ! !
    !           i = i + 1
    !           if ( i > 3 * n_sps) EXIT
    ! !
    !         END DO
    ! !
    !       end do
    ! !
    !     end if
!
10  CONTINUE ! K_INFO_COUNT = kc
    Return
!
  End Subroutine CONVOLVE_ALL
!
end module CONVOLVE_ALL_M
! $Log$
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
