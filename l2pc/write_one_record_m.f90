module WRITE_ONE_RECORD_M
  use FORMAT_KEY_M, only: FORMAT_KEY
  use L2PC_FILE_PARAMETERS, only: MAX_NO_KEY_ADDR
  use L2PC_FILE_STRUCTURES, only: I_STAR, K_STAR, X_STAR
  use MACHINE, only: IO_ERROR
  use MLSCommon, only: I4, R4, R8
  implicit NONE
  private
  public :: WRITE_ONE_RECORD
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
! This writes the linearized state vector and the radiometer spectral power.
!
  Subroutine WRITE_ONE_RECORD ( Key, jkey, state_vector, radiances, &
 &           derivatives, Keys, rec_nos, recn, mbr, bin_no, l2u,    &
 &           no_pointings, no_channels_per_band, Which, Ier)
!
    character(len=*), intent(inout) :: KEY
    integer(i4), intent(in) :: JKEY
    type(x_star), intent(in) :: STATE_VECTOR
    type(i_star), intent(in) :: RADIANCES
    type(k_star), intent(in) :: DERIVATIVES
    character(len=*), intent(inout) :: KEYS(*)
    integer(i4), intent(inout) :: REC_NOS(*)
    integer(i4), intent(inout) :: RECN
    integer(i4), intent(in) :: MBR
    integer(i4), intent(in) :: BIN_NO
    integer(i4), intent(in) :: L2U
    integer(i4), intent(in) :: NO_POINTINGS
    integer(i4), intent(in) :: NO_CHANNELS_PER_BAND
    character(len=*), intent(in) :: WHICH
    integer(i4), intent(out) :: IER
!     Include 'inc/l2pcdim.inc'
!     Include 'inc/l2pc_file_parameters.inc'
!     Include 'inc/l2pc_file_structures.inc'
!
    integer(i4) :: i, io, irec, j, k, l
    character(len=*), parameter :: MSG = &
      'I had a problem writing a record,sorry !'
    integer(i4) :: n
!
    io = 0
    if (recn == max_no_key_addr) then
      Print *,'** Maximum number of records reached .. recn =',recn
      goto 99
    end if
!
    Call format_key(Key,Ier)
    if (Ier /= 0) goto 99
!
    recn = recn + 1
    irec = 3 + (bin_no - 1) * mbr + recn
!
! Transfer the input record into the write out record
!
    n = no_pointings
    if (Which(1:12) == 'state_vector') then               ! x_star
      k = state_vector%no_sv_elmnts
      write(l2u,rec=irec,iostat=io)  Key,n,k,                          &
   &      (state_vector%pointings(j),j=1,n),                           &
   &      (state_vector%sv_elmnts(l),l=1,k)
    else if (Which(1:13) == 'the_radiances') then         ! i_star
      k = no_channels_per_band
      i = radiances%first_channel_number
      write(l2u,rec=irec,iostat=io)  Key,i,                            &
   &      ((radiances%radiances(j,l),l=1,k),j=1,n)
    else if (Which(1:14) == 'the_derivative') then        ! k_star
      k = no_channels_per_band
      i = derivatives%first_channel_number
      write(l2u,rec=irec,iostat=io)  Key,i,derivatives%sv_name1,       &
   &   derivatives%sv_name2,derivatives%log_psv1,derivatives%log_psv2, &
   &   ((derivatives%derivatives(j,l),l=1,k),j=1,n)
    else
      Print *, ' ** Unknown type of record: ', Which
      goto 99
    end if
!
    if (io /= 0) goto 99
!
! Store the key (only on first pass)
!
    if (jkey < 2) then
      rec_nos(recn) = recn
      Keys(recn)(1:32) = Key(9:40)
    end if
!
    Return
!
! Error Message
!
99  Ier = 1
!   Call ErrMsg('From subroutine: write_one_record',io)
    call io_error ( 'While writing in WRITE_ONE_RECORD', io )
    Print *, Msg
!
    Return
  End Subroutine WRITE_ONE_RECORD
end module WRITE_ONE_RECORD_M
! $Log$
! Revision 1.1  2000/05/04 18:12:06  vsnyder
! Initial conversion to Fortran 90
!
