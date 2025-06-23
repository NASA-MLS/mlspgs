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
  public :: Dump_Path, Sps_List

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains

  ! --------------------------------------------------  Dump_Path  -----
  subroutine Dump_Path ( Details, Config, I_start, I_stop, I_end, Phi, Zeta, &
    &                    VMR, Beta, Alpha, Incoptdepth, IncRadPath, T, Frq_i,&
    &                    Tau, Frq )

    use Constants, only: Rad2Deg
    use Dump_0, only: Dump
    use ForwardModelConfig, only: ForwardModelConfig_t
    use Intrinsic, only: Lit_indices
    use MLSKinds, only: RP, R8
    use Output_m, only: Blanks, NewLine, Output
    use String_Table, only: Display_String, String_Length
    use Tau_m, only: Tau_t

    integer, intent(in) :: Details ! 0 => zeta and phi only
                                   ! 1 => phi, zeta, alpha, T, IncRadPath, Tau,
                                   !      IncOptDepth
                                   ! >1 => 1 + betas
    type(ForwardModelConfig_t), intent(in) :: Config
    integer, intent(in) :: I_start, I_Stop, I_end
    real(rp), intent(in) :: Phi(:)
    real(rp), intent(in) :: Zeta(:)
    real(rp), intent(in) :: Beta(:,:)
    real(rp), intent(in) :: VMR(:,:)
    real(rp), intent(in) :: Alpha(:)
    real(rp), intent(in) :: Incoptdepth(:)
    real(rp), intent(in) :: IncRadPath(:)
    real(rp), intent(in) :: T(:)
    integer, intent(in) :: Frq_i
    type(Tau_t), intent(in) :: Tau
    real(r8), intent(in) :: Frq

    integer :: I, J
    character(87) :: Line ! of output

    call output ( frq, before='Frequency (MHz) = ', advance='yes' )

    if ( details == 0 ) then
      call dump ( rad2deg*phi(i_start:i_end), format='(f10.4)', name='Phi', lbound=i_start )
      call dump ( zeta(i_start:i_end), format='(f10.4)', name='Zeta', lbound=i_start )
      if ( i_stop < i_end ) then
        call output ( i_stop+1, before='Blacked out from ' )
        call output ( i_end, before=' to ', advance='yes' )
      end if
      return
    end if

    if ( details > 1 ) then
      ! Dump VMR and Beta
      call blanks ( 3 )
      do j = 1, size(vmr,2)
        call blanks ( 4 )
        call display_string ( lit_indices(config%beta_group(j)%molecule) )
        call blanks ( 1 )
        call output ( repeat('-',23-string_length(lit_indices(config%beta_group(j)%molecule))) )
      end do
      call NewLine
      call blanks ( 3 )
      do j = 1, size(vmr,2)
        call output ( "    VMR           Beta" )
        if ( j /= size(vmr,2) ) call blanks ( 6 )
      end do
      call NewLine
      do i = i_start, i_stop
        call output ( i, format='(i3)' )
        do j = 1, size(vmr,2)
          call output ( vmr(i,j), format='(1pg14.6)' )
          call output ( beta(i,j), format='(1pg14.6)' )
        end do
        call NewLine
      end do
    end if

    ! Dump Alpha, Incoptdepth, Tau, IncRadPath, T

    call output ( "         Phi     Zeta    Alpha          T    IncRadPath     Tau       IncOptDepth", advance="yes" )
    write ( line, 1 ) i_start, rad2deg*phi(i_start), zeta(i_start), &
      & alpha(i_start), T(i_start), incradpath(i_start), tau%tau(i_start,frq_i)
1   format ( i3, f10.4, f8.4, 1p, g14.6, 0p, f7.2, 1p, 3g14.6 )
    call output ( trim(line), advance='yes' )
    do i = i_start+1, tau%i_stop(frq_i)
      write ( line, 1 ) i, rad2deg*phi(i), zeta(i), alpha(i), T(i), &
        & incradpath(i), tau%tau(i,frq_i), incoptdepth(i)
      call output ( trim(line), advance='yes' )
    end do

    if ( i_stop < i_end ) then
      call output ( i_stop+1, before='Blacked out from ' )
      call output ( i_end, before=' to ', advance='yes' )
    end if

  end subroutine Dump_Path

  subroutine Sps_List ( Config )
  ! Print a list of active species
    use ForwardModelConfig, only: ForwardModelConfig_t
    use Intrinsic, only: Lit_indices
    use Output_m, only: Blanks, NewLine, Output
    use String_Table, only: Display_String, String_Length

    type(ForwardModelConfig_t), intent(in) :: Config

    integer :: I, J

    j = 0
    call output ( 'Species:', advance='yes' )
    call blanks ( 7 )
    do i = 1, size(config%beta_group)
      call output ( i, after=': ', format='(i3)' )
      call display_string ( lit_indices(config%beta_group(i)%molecule) )
      j = j + 1
      if ( j < 5 ) then
        call blanks ( max(0,9 - string_length(lit_indices(config%beta_group(i)%molecule))) )
      else
        call newLine
        j = 0
      end if
    end do
    if ( j /= 0 ) call newLine

  end subroutine Sps_List

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Dump_Path_m

! $Log$
! Revision 2.4  2011/11/09 00:30:05  vsnyder
! Add details argument
!
! Revision 2.3  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.2  2009/06/13 00:53:16  vsnyder
! Different meaning for i_start, simplify I/O
!
! Revision 2.1  2008/12/18 02:58:29  vsnyder
! Initial commit
!
