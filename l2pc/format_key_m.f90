module FORMAT_KEY_M
  use MACHINE, only: IO_ERROR
  use MLSCommon, only: I4
  use STRINGS, only: LEFTJ, SQZSTR, STRUPR
  implicit NONE
  private
  public :: FORMAT_KEY
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
contains
!  Parse & Format l2pc Key = 'D_10N_A_NA_THE_C_SVNAMEXXEL_PHI+NF_B'
!  or: (With Spectroscopy variable):           SVNAME_XEL_PHI+NF_B'
!
  Subroutine FORMAT_KEY ( l2pc_key, ier )
    character(len=*), intent(inout) :: L2PC_KEY
    Integer(i4), intent(out) :: IER
!
! -----     Local variables     ----------------------------------------
!
    character(len=40) :: AKEY
    character(len=8) :: AN1
    character(len=20) :: AX
    character :: CH
    Integer(i4) :: BAND, COEFF, I, IO, L, PHI_INDX
!
! -----     Executable statements     ----------------------------------
!
    ier = 0
    Akey = l2pc_key
!
    i = index(Akey,'__')
    do while (i > 0)
      Akey(i:i) = ' '
      i = index(Akey,'__')
    end do
!
    do l = len_trim(Akey), 1, -1
      if (Akey(l:l) /= '_') exit
      Akey(l:l) = ' '
    end do
!
    Call SqzStr(Akey)
    Call StrUpr(Akey)
!
    i = Index(Akey,'X_STAR')
    if (i < 1) i = Index(Akey,'I_STAR')
!
    if (i > 0) then
      l2pc_key=' '
      l2pc_key = Akey
      Return
    end if
!
!  Key = 'D_10N_A_NA_THE_C_SVNAMEXXEL_PHI+NF_B'   or:
!        'D_10N_A_NA_THE_C_SVNAME_XEL_PHI+NF_B'   (Spectroscopy type)
!
    l = len_trim(Akey)
    read (Akey(l:l),*,iostat=io) band
    if (io /= 0) goto 99
    Akey(l-1:l)=' '
    l = l - 2
!
    do i = l, 1, -1
      Ch = Akey(i:i)
      if (Ch == '+' .or. Ch == '-') then
        An1(1:)=' '
        An1 = Akey(i:l)
        read (An1,*,iostat=io) phi_indx
        if (io /= 0) goto 99
        Akey(i-4:l) = ' '
        l = len_trim(Akey)
        exit
      end if
    end do
!
    io = 5
    An1(1:)=' '
    i = l + 1
    do i = l, 1, -1
      if (io <= 3) exit
      Ch = Akey(i:i)
      Akey(i:i) = ' '
      if (Index('@0123456789',Ch) <= 0) exit
      io = io - 1
      An1(io:io) = Ch
    end do
!
    read (An1,*,iostat=io) coeff
    if (io /= 0) goto 99
!
    Ax = Akey(1:17)
    Akey(1:17)=' '
    Call Leftj(Akey)
    l = len_trim(Akey)
    if (Akey(l:l) == '_') l = l - 1
    An1 = Akey(1:l)
!
    Akey(1:17) = Ax(1:17)
    Akey(18:25) = An1(1:8)
    write (Akey(26:27),'(i2.2)') coeff
    Akey(28:32)='_PHI+'
    i = iabs(phi_indx)
    if (phi_indx < 0) Akey(32:32) = '-'
    write (Akey(33:34),'(i2.2)') i
    Akey(35:35)='_'
    write (Akey(36:36),'(i1)') band
    Akey(37:)=' '
    Call SqzStr(Akey)
!
    l2pc_key(1:)=' '
    l2pc_key = Akey
    Return
!
 99 ier = io
    Print *,' Input l2pc_key: ',l2pc_key
!   Call ErrMsg('Error in subroutine: format_key',ier)
    call io_error ( 'In FORMAT_KEY', ier, ' internal file ' )
!
    Return
  End Subroutine FORMAT_KEY
end module FORMAT_KEY_M
! $Log$
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
!
