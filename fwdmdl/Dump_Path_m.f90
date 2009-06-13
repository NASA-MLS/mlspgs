! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Dump_Path_m

! Dump path quantities in a readable format


  implicit NONE
  private
  public :: Dump_Path

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains

  ! --------------------------------------------------  Dump_Path  -----
  subroutine Dump_Path ( Config, I_start, I_end, Phi, Zeta, VMR, Beta, &
    &                    Alpha, Incoptdepth, Frq_i, Tau, IncRadPath )

    use Constants, only: Rad2Deg
    use ForwardModelConfig, only: ForwardModelConfig_t
    use Intrinsic, only: Lit_indices
    use MLSKinds, only: RP
    use Output_m, only: Blanks, NewLine, Output
    use String_Table, only: Display_String, String_Length
    use Tau_m, only: Tau_t

    type(ForwardModelConfig_t), intent(in) :: Config
    integer, intent(in) :: I_start, I_end
    real(rp), intent(in) :: Phi(:)
    real(rp), intent(in) :: Zeta(:)
    real(rp), intent(in) :: Beta(:,:)
    real(rp), intent(in) :: VMR(:,:)
    real(rp), intent(in) :: Alpha(:)
    real(rp), intent(in) :: Incoptdepth(:)
    integer, intent(in) :: Frq_i
    type(Tau_t), intent(in) :: Tau
    real(rp), intent(in) :: IncRadPath(:)

    integer :: I, J
    character(80) :: Line ! of output

    ! Dump VMR and Beta

!   call blanks ( len("         Phi     Zeta") )
    call blanks ( 3 )
    do j = 1, size(vmr,2)
      call blanks ( 4 )
      call display_string ( lit_indices(config%beta_group(j)%molecule) )
      call blanks ( 1 )
      call output ( repeat('-',23-string_length(lit_indices(config%beta_group(j)%molecule))) )
    end do
    call NewLine
!   call output (     "         Phi     Zeta" )
    call blanks ( 3 )
    do j = 1, size(vmr,2)
      call output ( "    VMR           Beta" )
      if ( j /= size(vmr,2) ) call blanks ( 6 )
    end do
    call NewLine
    do i = i_start, i_end
      call output ( i, format='(i3)' )
!     call output ( rad2deg*phi(i), format='(f10.4)' )
!     call output ( zeta(i), format='(f8.4)' )
      do j = 1, size(vmr,2)
        call output ( vmr(i,j), format='(1pg14.6)' )
        call output ( beta(i,j), format='(1pg14.6)' )
      end do
      call NewLine
    end do

    ! Dump Alpha, Incoptdepth, Tau, IncRadPath

    call output ( "         Phi     Zeta    Alpha        IncOptDepth     Tau       IncRadPath", advance="yes" )
    write ( line, 1 ) i_start, rad2deg*phi(i_start), zeta(i_start), &
      & alpha(i_start)
1   format ( i3, f10.4, f8.4, 1p, 4g14.6 )
    call output ( trim(line), advance='yes' )
    do i = i_start+1, tau%i_stop(frq_i)
      write ( line, 1 ) i, rad2deg*phi(i), zeta(i), alpha(i), incoptdepth(i), &
        & tau%tau(i,frq_i), incradpath(i)
      call output ( trim(line), advance='yes' )
    end do

  end subroutine Dump_Path

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, not_used_here ! .mod files sometimes change if PRINT is added
  end function not_used_here

end module Dump_Path_m

! $Log$
! Revision 2.2  2009/06/13 00:53:16  vsnyder
! Different meaning for i_start, simplify I/O
!
! Revision 2.1  2008/12/18 02:58:29  vsnyder
! Initial commit
!
