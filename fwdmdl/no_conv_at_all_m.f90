! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module NO_CONV_AT_ALL_M
  use MLSCommon, only: I4, R4, R8
  use L2PC_PFA_STRUCTURES, only: K_MATRIX_INFO
  use D_LINTRP_M, only: LINTRP
  use D_CSPLINE_M, only: CSPLINE
  use DCSPLINE_DER_M, only: CSPLINE_DER
  use dump_0,only:dump
  use VectorsModule, only: Vector_T, VectorValue_T, GETVECTORQUANTITYBYTYPE
  use ForwardModelConfig, only: ForwardModelConfig_T
  use Intrinsic, only: L_VMR
  use String_Table, only: GET_STRING
  use MatrixModule_0, only: M_ABSENT, M_BANDED, M_FULL, DUMP
  use MatrixModule_1, only: CREATEBLOCK, FINDBLOCK, MATRIX_T, DUMP
  use Molecules, only: L_EXTINCTION
  use MLSMessageModule, only: MLSMSG_Error, MLSMessage

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
contains
  !-------------------------------------------------------------------------
  ! This subroutine transfers the derivatives over from the internal
  ! convolution grid to the users specified points. This module uses
  ! cubic spline interpolation to do the job.

  Subroutine no_conv_at_all ( ForwardModelConfig, ForwardModelIn, MAF, &
    & Channel, WindowStart, WindowFinish, Temp, Ptan, Radiance, &
    & Tan_press, I_raw, K_temp, K_atmos, SbRatio, Jacobian, rowFlags )

    type (ForwardModelConfig_T) :: FORWARDMODELCONFIG
    type (Vector_T), intent(in) :: FORWARDMODELIN
    integer, intent(in) :: MAF
    integer, intent(in) :: CHANNEL
    integer, intent(in) :: WINDOWSTART
    integer, intent(in) :: WINDOWFINISH
    type (VectorValue_T), intent(in) :: TEMP
    type (VectorValue_T), intent(in) :: PTAN
    type (VectorValue_T), intent(inout) :: RADIANCE

    real(r8), intent(IN) :: I_RAW(:),TAN_PRESS(:)
    real(r8), intent(IN) :: SBRATIO
!
    Real(r4), intent(in) :: k_temp(:,:,WindowStart:)
    Real(r4), intent(in) :: k_atmos(:,:,:,WindowStart:,:)

    type (Matrix_T), intent(inout), optional :: Jacobian
    logical, dimension(:), intent(inout) :: rowFlags ! Flag to calling code

    ! -----     Local Variables     ------------------------------------

    integer:: No_t, No_tan_hts, No_phi_t
    type (VectorValue_T), pointer :: F  ! vmr quantity
    integer :: Lk, Uk, jf, jz, n_sps
    integer :: IS, J, K, NF, SV_I, Spectag
    integer :: Row, col                 ! Indices
    integer :: Ptg                      ! Index
    integer :: Ind                      ! Index

    real(r8) :: RAD( size(tan_press)), Q
    real(r8) :: SRad(ptan%template%noSurfs)
    real(r8) :: Der_all(ptan%template%noSurfs)
    real(r8) :: I_star_all(ptan%template%noSurfs)

    ! -----  Begin the code  -------------------------------------------

    no_t = temp%template%noSurfs
    no_phi_t = temp%template%noInstances
    no_tan_hts = size(tan_press)

    ! Compute the ratio of the strengths

    ! This subroutine is called by channel

    k = no_tan_hts
    j = ptan%template%noSurfs
    if ( present(Jacobian) ) then
      Call Cspline_der ( tan_press, Ptan%values(:,maf), i_raw, i_star_all, &
        & der_all, k, j )

      row = FindBlock ( Jacobian%row, radiance%index, maf )
      col = FindBlock ( Jacobian%col, ptan%index, maf )
      rowFlags(row) = .true.
      select case ( jacobian%block(Row,col)%kind )
      case ( m_absent )
        call CreateBlock ( Jacobian, row, col, m_banded, &
          & radiance%template%noSurfs*radiance%template%noChans, &
          & bandHeight=radiance%template%noChans )
        jacobian%block(row,col)%values = 0.0_r8
      case ( m_banded )
      case default
        call MLSMessage ( MLSMSG_Error, ModuleName,&
          & 'Wrong matrix type for ptan derivative' )
      end select
      do ptg = 1, j
        ind = channel + radiance%template%noChans*(ptg-1)
        jacobian%block(row,col)%values(ind ,1 ) = &
          & jacobian%block(row,col)%values( ind ,1 ) + sbRatio*der_all(ptg)
      end do
    else
      Call Cspline ( tan_press, Ptan%values(:,maf), i_raw, i_star_all, k, j )
    end if
    do ptg = 1, j
      ind = channel + radiance%template%noChans*(ptg-1)
      radiance%values( ind, maf ) = &
        & radiance%values ( ind, maf ) + sbRatio*i_star_all(ptg)
    end do


    if ( .not. ANY((/forwardModelConfig%temp_der, forwardModelConfig%atmos_der, &
      & forwardModelConfig%spect_der/)) ) Return
    if ( .not. present(jacobian) ) return


    ! Now transfer the other fwd_mdl derivatives to the output pointing
    ! values

    ! ********************* Temperature derivatives ******************

    ! check to determine if derivative is desired for this parameter

    if ( forwardModelConfig%temp_der ) then

    ! Derivatives needed continue to process

      Rad(1:) = 0.0
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
          Rad(1:k) = k_temp(1:k,jz,nf)
          Call Cspline ( tan_press, Ptan%values(:,maf), Rad, SRad, k, j )
          do ptg = 1,j
            ind = channel + radiance%template%noChans*(ptg-1)
            jacobian%block(row,col)%values( ind, jz) = &
              & jacobian%block(row,col)%values( ind, jz) + sbRatio*Srad(ptg)
          end do
        end do

      end do


    end if

    if ( forwardModelConfig%atmos_der ) then

      ! ****************** atmospheric derivatives ******************

      lk = lbound(k_atmos,4)   ! The lower Phi dimension
      uk = ubound(k_atmos,4)   ! The upper Phi dimension
      n_sps = size(ForwardModelConfig%molecules)

      do is = 1, n_sps

        if ( forwardModelConfig%molecules(is) == l_extinction ) then
          f => GetVectorQuantityByType(forwardModelIn, quantityType=l_extinction, &
            &  radiometer=radiance%template%radiometer, noError=.true. )
        else
          f => GetVectorQuantityByType ( forwardModelIn, quantityType=l_vmr, &
            &  molecule=forwardModelConfig%molecules(is), noError=.true. )
        endif

        if ( associated(f) ) then

          Rad(1:) = 0.0

          ! Derivatives needed continue to process

          do nf = 1, f%template%noInstances

          ! run through phi representation basis coefficients

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

            sv_i = 0
            do jz = 1, f%template%noSurfs

            ! run through zeta representation basis coefficients

              do jf = 1, f%template%noChans

                ! run through Frequencies basis coefficients

                sv_i = sv_i + 1
                Rad(1:k) = k_atmos(1:k,jf,jz,nf+lk-1,is)
                Call Lintrp(tan_press,Ptan%values(:,maf),Rad,SRad,k,j)
                do ptg = 1, j
                  ind = channel + radiance%template%noChans*(ptg-1)
                  q = jacobian%block(row,col)%values( ind, sv_i)
                  jacobian%block(row,col)%values( ind, sv_i) = &
                                                       & q + sbRatio*Srad(ptg)
                end do

              end do

            end do

          end do

        endif

      end do

    endif

    Return

  End Subroutine NO_CONV_AT_ALL

end module NO_CONV_AT_ALL_M
! $Log$
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

