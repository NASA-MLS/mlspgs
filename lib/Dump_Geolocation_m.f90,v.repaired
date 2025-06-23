! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module Dump_Geolocation_m
!=============================================================================

  implicit none
  private

  public :: Dump     ! Generic for Dump_H_t and Dump_ZOT
  public :: Dump_H_t ! Dump (lon,lat) array
  public :: Dump_ZOT ! Dump ZOT coordinates array

  interface Dump
    module procedure Dump_H_t, Dump_ZOT
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Dump_H_t ( G, Name, Format )

    ! Dump an array of lon,lat pairs enclosed in parentheses.

    use Geolocation_0, only: H_t
    use Output_m, only: NewLine, Output

    type(h_t), intent(in) :: G(:)
    character(*), intent(in), optional :: Name
    character(*), intent(in), optional :: Format ! Default '(f8.3)'

    integer :: I, J, N
    character(31) :: MyFmt

    myFmt = '(f8.3)'
    if ( present(format) ) myFmt = format
    if ( present(name) ) call output ( name, advance='yes' )
    n = size(g)
    do i = 1, n, 5
      call output ( i, 4, after=':' )
      do j = i, min(i+4,n)
        call output ( g(j)%lon%d, before=' (', format=myFmt )
        call output ( g(j)%lat, before=',', format=myFmt, after=')' )
      end do
      call NewLine
    end do

  end subroutine Dump_H_t

  subroutine Dump_ZOT ( G, Name, Format )

    ! Dump an array of x,y pairs enclosed in parentheses.

    use QTM_m, only: ZOT_t
    use Output_m, only: NewLine, Output

    type(ZOT_t), intent(in) :: G(:)
    character(*), intent(in), optional :: Name
    character(*), intent(in), optional :: Format ! Default '(f8.5)'

    integer :: I, J, N
    character(31) :: MyFmt

    myFmt = '(f8.5)'
    if ( present(format) ) myFmt = format
    if ( present(name) ) call output ( name, advance='yes' )
    n = size(g)
    do i = 1, n, 5
      call output ( i, 4, after=':' )
      do j = i, min(i+4,n)
        call output ( g(j)%x, before=' (', format=myFmt )
        call output ( g(j)%y, before=',', format=myFmt, after=')' )
      end do
      call NewLine
    end do

  end subroutine Dump_ZOT

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Dump_Geolocation_m

! $Log$
! Revision 2.3  2016/08/23 00:39:23  vsnyder
! Add generic Dump for Dump_H_t and Dump_ZOT
!
! Revision 2.2  2016/03/25 00:25:06  vsnyder
! Lon component now needs to acces its %d component
!
! Revision 2.1  2016/02/23 00:46:49  pwagner
! First commit
!
