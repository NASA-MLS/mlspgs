module Test_Parse_Signals_m

! After ingesting the l2cf, it will print "Enter radiometer strings:"
! When you enter a radiometer string, it will parse it, and dump all
! of the signals in the database that match it.

! When you're bored with this, enter end-of-file to continue with the
! real work.

  use Allocate_Deallocate, only: Deallocate_Test
  use Init_Tables_Module, only: Spec_Indices
  use Parse_Signal_m, only: Parse_Signal
  use MLSSignals_m, only: GetSignalName

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains
  subroutine Test_Parse_Signals

    logical, pointer, dimension(:) :: Channels => NULL()
    integer :: I
    character(len=127) :: Line
    integer :: Sideband
    integer, pointer, dimension(:) :: Signal_Indices => NULL()

    print *, 'Enter radiometer strings: '
    do
      read ( *, '(a)', end=99 ) line
      call parse_signal ( line, signal_indices, spec_indices, &
        & sideband=sideband, channels=channels )
      if ( associated(signal_indices) ) then
        print '(10i5)', signal_indices
        do i = 1, size(signal_indices)
          call getSignalName ( signal_indices(i), line )
          print *, trim(line)
        end do
        print *, 'Sideband =', sideband
        print '(" Channels = ", 50l1:(/12x,50l1))', channels
      end if
    end do
99  return
  end subroutine Test_Parse_Signals
end module Test_Parse_Signals_m

! $Log$
! Revision 2.2  2001/04/06 20:12:14  vsnyder
! Print Sideband and Channels
!
! Revision 2.1  2001/03/15 21:05:23  vsnyder
! Set up to test Parse_Signal
!
