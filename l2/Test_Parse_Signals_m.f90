! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Test_Parse_Signals_m

! After ingesting the l2cf, it will print "Enter radiometer strings:"
! When you enter a radiometer string, it will parse it, and dump all
! of the signals in the database that match it.

! When you're bored with this, enter end-of-file to continue with the
! real work.

  use Allocate_Deallocate, only: Deallocate_Test
  use Parse_Signal_m, only: Parse_Signal
  use MLSSignals_m, only: GetSignalName

  implicit NONE

  private

  public :: Test_Parse_Signals

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains
  subroutine Test_Parse_Signals

    logical, pointer, dimension(:) :: Channels
    integer :: I
    character(len=127) :: Line
    integer :: Sideband
    integer, pointer, dimension(:) :: Signal_Indices

    nullify ( channels, signal_Indices )

    print *, 'Enter radiometer strings: '
    do
      read ( *, '(a)', end=99 ) line
      call parse_signal ( line, signal_indices, &
        & sideband=sideband, channels=channels )
      if ( associated(signal_indices) ) then
        print '(10i5)', signal_indices
        do i = 1, size(signal_indices)
          call getSignalName ( signal_indices(i), line )
          print *, trim(line)
        end do
        print *, 'Sideband =', sideband
        if ( associated(channels) ) &
          & print '(" Channels = ", 50l1:(/12x,50l1))', channels
      end if
    end do
99  return
  end subroutine Test_Parse_Signals
end module Test_Parse_Signals_m

! $Log$
! Revision 2.7  2002/03/12 23:43:26  vsnyder
! Added IMPLICIT NONE and PRIVATE statements
!
! Revision 2.6  2001/06/07 21:58:28  pwagner
! Added Copyright statement
!
! Revision 2.5  2001/04/26 02:44:17  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.4  2001/04/10 22:27:47  vsnyder
! Nullify explicitly instead of with <initialization> so as not to give
! pointers the SAVE attribute.  <initialization> is NOT executed on each
! entry to a procedure.
!
! Revision 2.3  2001/04/10 17:58:57  vsnyder
! Handle possibility that 'channels' is not associated
!
! Revision 2.2  2001/04/06 20:12:14  vsnyder
! Print Sideband and Channels
!
! Revision 2.1  2001/03/15 21:05:23  vsnyder
! Set up to test Parse_Signal
!
