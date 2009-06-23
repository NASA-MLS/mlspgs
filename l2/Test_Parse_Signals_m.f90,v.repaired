! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

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

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
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
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Test_Parse_Signals_m

! $Log$
! Revision 2.10  2009/06/23 18:46:19  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.9  2005/06/22 18:57:02  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.8  2002/10/08 17:36:23  pwagner
! Added idents to survive zealous Lahey optimizer
!
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
