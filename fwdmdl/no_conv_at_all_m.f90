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

  Subroutine no_conv_at_all ( ForwardModelConfig, ForwardModelIn, maf, &
           & Channel, WindowStart, WindowFinish, Temp, Ptan, Radiance, &
           & ptg_angles, chi_out, I_raw, K_temp, K_atmos, SbRatio,     &
           & Jacobian, rowFlags, mol_cat_indx )

    Type (ForwardModelConfig_T) :: FORWARDMODELCONFIG
    Type (Vector_T), intent(in) :: FORWARDMODELIN

    Integer, intent(in) :: maf
    Integer, intent(in) :: CHANNEL
    Integer, intent(in) :: WINDOWSTART
    Integer, intent(in) :: WINDOWFINISH
    Integer, intent(IN) :: mol_cat_indx(:)

    Type (VectorValue_T), intent(in) :: TEMP
    Type (VectorValue_T), intent(in) :: PTAN
    Type (VectorValue_T), intent(inout) :: RADIANCE

    Real(r8), intent(IN) :: SBRATIO
    Real(r8), intent(IN) :: i_raw(:),ptg_angles(:),chi_out(:)
!
    Real(r4), intent(in) :: k_temp(:,:,WindowStart:)
    Real(r4), intent(in) :: k_atmos(:,:,:,WindowStart:,:)

    Type (Matrix_T), intent(inout), optional :: Jacobian

    Logical, dimension(:), intent(inout) :: rowFlags ! Flag to calling code

    ! -----     Local Variables     ------------------------------------

    Integer:: No_t, No_tan_hts

    Type (VectorValue_T), pointer :: F  ! vmr quantity

    Integer :: Lk, Uk, jf, jz, no_mol, l
    Integer :: is, j, k, nf, sv_i
    Integer :: Row, col                     ! Matrix row & column indices
    Integer :: ptg_i,noPtan,noChans,Ind     ! Indices

    Real(r8) :: RAD( size(ptg_angles)), q
    Real(r8) :: SRad(ptan%template%noSurfs)
    Real(r8) :: Der_all(ptan%template%noSurfs)
    Real(r8) :: I_star_all(ptan%template%noSurfs)

    ! -----  Begin the code  -------------------------------------------

    no_t = temp%template%noSurfs
    no_tan_hts = size(ptg_angles)

    ! Compute the ratio of the strengths

    ! This subroutine is called by channel

    col = 0
    noPtan = ptan%template%noSurfs
    noChans = radiance%template%noChans

    ! Ptan derivative
    if ( present(Jacobian) ) &
      & col = FindBlock ( Jacobian%col, ptan%index, maf )

    if ( col > 0 ) then

      Call InterpolateValues(ptg_angles, i_raw, chi_out, i_star_all, 'S')
      do k = 1, noPtan
        j = 1
        if(k == noPtan) j = -1
        q = Ptan%values(k+j,maf) - Ptan%values(k,maf)
        der_all(k) = (i_star_all(k+j) - i_star_all(k) ) / q
      end do

      row = FindBlock ( Jacobian%row, radiance%index, maf )
      rowFlags(row) = .true.

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
        jacobian%block(row,col)%values(ind ,1 ) = &
          & jacobian%block(row,col)%values( ind ,1 ) + sbRatio*der_all(ptg_i)
      end do
    else
      Call InterpolateValues ( ptg_angles, i_raw, chi_out, i_star_all, &
                            &  METHOD='S')
    end if

    do ptg_i = 1, noPtan
      ind = channel + noChans*(ptg_i-1)
      radiance%values( ind, maf ) = &
        & radiance%values ( ind, maf ) + sbRatio*i_star_all(ptg_i)
    end do

    if ( .not. ANY((/forwardModelConfig%temp_der, &
                   & forwardModelConfig%atmos_der, &
                   & forwardModelConfig%spect_der/)) ) Return
    if ( .not. PRESENT(jacobian) ) Return

    ! Now transfer the other fwd_mdl derivatives to the output pointing
    ! values

    ! ********************* Temperature derivatives ******************

    ! check to determine if derivative is desired for this parameter

    if ( forwardModelConfig%temp_der ) then

    ! Derivatives needed continue to process

      Rad(1:) = 0.0
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
          Rad(1:k) = k_temp(1:k,jz,nf)
          Call InterpolateValues ( ptg_angles, Rad, chi_out, Srad, 'S')
          do ptg_i = 1,noPtan
            ind = channel + noChans*(ptg_i-1)
            jacobian%block(row,col)%values( ind, jz) = &
              & jacobian%block(row,col)%values( ind, jz) + sbRatio*Srad(ptg_i)
          end do
        end do

      end do

    end if

    if ( forwardModelConfig%atmos_der ) then

      ! ****************** atmospheric derivatives ******************

      lk = lbound(k_atmos,4)   ! The lower Phi dimension
      uk = ubound(k_atmos,4)   ! The upper Phi dimension
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

        if ( associated(f) ) then

          Rad(1:) = 0.0
          k = no_tan_hts

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
                Call InterpolateValues (ptg_angles, Rad, chi_out, Srad, 'L')
                do ptg_i = 1, noPtan
                  ind = channel + noChans*(ptg_i-1)
                  q = jacobian%block(row,col)%values( ind, sv_i)
                  jacobian%block(row,col)%values( ind, sv_i) = &
                                                  & q + sbRatio*Srad(ptg_i)
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
