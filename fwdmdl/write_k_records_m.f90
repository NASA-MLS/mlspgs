module WRITE_K_RECORDS_M
  use L2PC_FILE_STRUCTURES, only: DUMP_K_STAR, I_STAR, K_STAR, &
                                  L2PC_HEADER_ONE, X_STAR
  use MACHINE, only: IO_ERROR
  use MLSCommon, only: I4, R4, R8
  use WRITE_ONE_RECORD_M, only: WRITE_ONE_RECORD
  implicit NONE
  private
  public :: WRITE_K_RECORDS

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

  Subroutine WRITE_K_RECORDS ( time_i, l2pc_lu, recn, mxbin_nor,         &
 &           lrun, rec_nos, band, Keys, runf, header1, jkey, cycle, Ier )
!
    integer(i4), intent(in) :: TIME_I
    integer(i4), intent(in) :: L2PC_LU
    integer(i4), intent(inout) :: RECN
    integer(i4), intent(in) :: MXBIN_NOR
    integer(i4), intent(in) :: LRUN
    integer(i4), intent(inout) :: REC_NOS(*)
    integer(i4), intent(in) :: BAND               ! Not used
    character(len=*), intent(inout) :: KEYS(*)
    character(len=*), intent(in) :: RUNF
    type(l2pc_header_one), intent(in) :: HEADER1
    integer(i4), intent(in) :: JKEY
    integer(i4), intent(in) :: CYCLE
    integer(i4), intent(out) :: IER

    type(dump_k_star) :: DUMP_K
    character(len=*), parameter :: FDUM = 'k_star_dummy.dat'
    integer(i4) :: IO, J, JCH, K
    character(len=len(keys)) :: KEY
    character(len=76) :: Line
    integer(i4) :: M
    real(r4) :: Q
    type(x_star) :: STATE_VECTOR
    type(k_star) :: THE_DERIVATIVES
    type(i_star) :: THE_RADIANCES ! Not assigned a value?

!     Character*40 Key
!
    inquire ( iolength=j ) dump_k
    close (97,iostat=io)
    open (97,file=Fdum,status='OLD',access='DIRECT',                   &
   &      form='UNFORMATTED',ACTION='READ',recl=j,iostat=Ier)
    if (Ier /= 0) then
      call io_error ( 'Opening file in WRITE_K_RECORDS', ier, fdum )
      Return
    end if
!
    do m = 1, cycle
      k = m - cycle
      do jch = 1, header1%no_channels_per_band
        k = k + cycle
        Read(97,rec=k,iostat=io) dump_k
        if (io /= 0) then
          Ier = 1
          Print *,'** Error in: write_k_recors subroutine ..'
          Line = 'Read(97,rec=k,iostat=io) dump_k'
          Print *,'   Last statement executed:'
!         Call ErrMsg(Line,io)
          call io_error ( line, io, Fdum )
          Return
        end if
        if (jch == 1) then
          Key = dump_k%Key
          j = dump_k%first_channel_number
          the_derivatives%first_channel_number = j
          the_derivatives%sv_name1 = dump_k%sv_name1
          the_derivatives%sv_name2 = dump_k%sv_name2
          the_derivatives%log_psv1 = dump_k%log_psv1
          the_derivatives%log_psv2 = dump_k%log_psv2
        end if
        do j = 1, header1%no_pointings
          q = dump_k%derivatives(j)
          the_derivatives%derivatives(j,jch) = q
        end do
      end do
      Call write_one_record(Key,jkey,state_vector,the_radiances,       &
   &       the_derivatives,Keys,rec_nos,recn,mxbin_nor,time_i,l2pc_lu, &
   &       header1%no_pointings,header1%no_channels_per_band,          &
   &       'the_derivatives',Ier)
      if (Ier /= 0) Return
      write (*,'(1x,a)') Key
      write (86,'(1x,a)') Key
    end do
!
    close (97, iostat=m, status='delete')
!   Call Get_Rid(Fdum)
!
!  Close run file & re-open for append (this is done to ensure information
!  up to this point is safe).
!
    close (86,iostat=m)
    open (86,file=runf(1:lrun),position='APPEND',iostat=m)
!
    Return
  End Subroutine WRITE_K_RECORDS
end module WRITE_K_RECORDS_M

! $Log$
