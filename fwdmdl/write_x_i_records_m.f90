module WRITE_X_I_RECORDS_M
  use L2PC_FILE_PARAMETERS, only: MAX_NO_BANDS, MAX_NO_SV_ELMNTS, &
                                  MXCO => max_no_elmnts_per_sv_component
  use L2PC_FILE_STRUCTURES, only: I_STAR, K_STAR, L2PC_HEADER_ONE, &
                                  L2PC_HEADER_TWO, X_STAR
  use L2PC_PFA_STRUCTURES, only: ATMOS_COMP, GEOM_PARAM, GEOPHYS_PARAM, &
                                 LIMB_PRESS
  use L2PCDim, only: MNP => max_no_phi
  use MLSCommon, only: I4, R4
  use WRITE_ONE_RECORD_M, only: WRITE_ONE_RECORD
  implicit NONE
  private
  public :: WRITE_X_I_RECORDS

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains
!
  Subroutine WRITE_X_I_RECORDS ( jch, time_i, l2pc_lu, recn, mxbin_nor, &
 &           lrun, rec_nos, band, no_geom, no_atmos, no_geophys,        &
 &           mr_g, mr_f, Keys, runf, ptg_press, geometric, geophysic,   &
 &           atmospheric, header1, header2, jkey, i_star_all,           &
 &           atmos_index, geom_index, geophys_index,                    &
 &           no_phi_g, no_phi_f, Ier )
!
    integer(i4), intent(in) :: JCH
    integer(i4), intent(in) :: TIME_I
    integer(i4), intent(in) :: L2PC_LU
    integer(i4), intent(inout) :: RECN
    integer(i4), intent(in) :: MXBIN_NOR
    integer(i4), intent(in) :: LRUN
    integer(i4), intent(inout) :: REC_NOS(*)
    integer(i4), intent(in) :: BAND
    integer(i4), intent(in) :: NO_GEOM
    integer(i4), intent(in) :: NO_ATMOS
    integer(i4), intent(in) :: NO_GEOPHYS
    real(r4), intent(in) :: MR_G(mxco,mnp,*)
    real(r4), intent(in) :: MR_F(mxco,mnp,*)
    character(len=*), intent(inout) :: KEYS(*)
    character(len=*), intent(in) :: RUNF
    type(limb_press), intent(in) :: PTG_PRESS
    type(geom_param), intent(in) :: GEOMETRIC(no_geom)
    type(geophys_param), intent(in) :: GEOPHYSIC(no_geophys)
    type(atmos_comp), intent(in) :: ATMOSPHERIC(no_atmos)
    type(l2pc_header_one), intent(in) :: HEADER1
    type(l2pc_header_two), intent(in) :: HEADER2
    integer(i4), intent(in) :: JKEY
    real(r4), intent(in) :: I_STAR_ALL(*)
    integer(i4), intent(in) :: ATMOS_INDEX(*)
    integer(i4), intent(in) :: GEOM_INDEX(*)
    integer(i4), intent(in) :: GEOPHYS_INDEX(*)
    integer(i4), intent(in) :: NO_PHI_G(*)
    integer(i4), intent(in) :: NO_PHI_F(*)
    integer(i4), intent(out) :: IER
!
    character(len=2) :: ABAND
    Character(len=40) :: AKEY
    integer(i4) :: I, J, JP, K
    Character(len=40) :: KEY
    integer(i4) :: M, NO_CHANNELS_PER_BAND, NO_POINTINGS, SV_I
    type(x_star) :: STATE_VECTOR
    type(k_star) :: THE_DERIVATIVES
    type(i_star), save :: THE_RADIANCES
    Character(len=40), save :: XKEY = '@@@@@@@@@@'
!
! Write the records by pointing and description
! Load in the state vector
!
    Ier = 0
    aband = '_b'
    write(aband(2:2),'(i1)') band
!
! Set a general key format:
!
    Key(1:) = ' '
    Key(01:08) = header1%avail_keys(time_i)
    Key(09:23) = 'NA_THE_C_X_STAR'
!
    no_pointings = header1%no_pointings
    no_channels_per_band = header1%no_channels_per_band
!
    if (Key /= XKey) then
!
      j = header1%no_pointings
      state_vector%no_pointings = j
      state_vector%pointings(1:j) = ptg_press%lin_val(1:j)
!
      state_vector%sv_elmnts(1:max_no_sv_elmnts) = 0.0
!
      do j = 1, no_geom
        if ( any(geometric(j)%der_calc(1:max_no_bands)) ) then
          i = geom_index(j)
          sv_i = header1%sv_component_first_elmnt_index(i)
          state_vector%sv_elmnts(sv_i) = geometric(j)%lin_val
        end if
      end do
!
      do j = 1, no_geophys
        if ( any(geophysic(j)%der_calc(1:max_no_bands)) ) then
          i = geophys_index(j)
          sv_i = header1%sv_component_first_elmnt_index(i) - 1
          jp = (no_phi_g(j) + 1) / 2           ! Temporary code ..
          do k = 1, geophysic(j)%no_lin_values
            sv_i = sv_i + 1
            state_vector%sv_elmnts(sv_i) = mr_g(k,jp,j)
          end do
        end if
      end do
!
      do j = 1, no_atmos
        if ( any(atmospheric(j)%der_calc(1:max_no_bands)) ) then
          i = atmos_index(j)
          sv_i = header1%sv_component_first_elmnt_index(i) - 1
          jp = (no_phi_f(j) + 1) / 2           ! Temporary code ..
          do k = 1, atmospheric(j)%no_lin_values
            sv_i = sv_i + 1
            state_vector%sv_elmnts(sv_i) = mr_f(k,jp,j)
          end do
        end if
      end do
!
      state_vector%no_sv_elmnts = header2%no_sv_elmnts
!
! Set the key for: x_star
!
      XKey = Key
      AKey = Key
      Call write_one_record(AKey,jkey,state_vector,the_radiances,      &
   &       the_derivatives,Keys,rec_nos,recn,mxbin_nor,time_i,l2pc_lu, &
   &       header1%no_pointings,header1%no_channels_per_band,          &
   &       'state_vector',Ier)
      if (Ier /= 0) Return
!
      write (*,'(1x,a)') AKey
      write (86,'(1x,a)') AKey
!
! Evaluate the initial channel index
!
      j = (band - 1) * no_channels_per_band
      the_radiances%first_channel_number = j + 1
!
    end if
!
! Second, i_star:
! Writing the /i_star/ record for channel: jch
! Set up the radiance vector, write records per channels
!
    the_radiances%radiances(1:no_pointings,jch) = i_star_all(1:no_pointings)
!
    if (jch == no_channels_per_band) then
!
! Write the /i_star/ records using no_sv_lin_val,sv_lin_val,
! and i_star. Set the key:
!
      Key(18:23) = 'I_STAR'
      Key(24:25) = aband
!
      AKey = Key
      Call write_one_record(AKey,jkey,state_vector,the_radiances,      &
   &       the_derivatives,Keys,rec_nos,recn,mxbin_nor,time_i,l2pc_lu, &
   &       header1%no_pointings,header1%no_channels_per_band,          &
   &       'the_radiances',Ier)
      if (Ier /= 0) Return
!
      write (*,'(1x,a)') AKey
      write (86,'(1x,a)') AKey
!
    end if
!
!  Close run file & re-open for append (this is done to ensure information
!  up to this point is safe).
!
    Close (86,iostat=m)
    Open (86,file=runf(1:lrun),position='APPEND',iostat=m)
!
    Return
  End Subroutine WRITE_X_I_RECORDS
end module WRITE_X_I_RECORDS_M

! $Log$
