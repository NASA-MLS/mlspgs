! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Read_Mie_m

  ! Read the Mie tables produced by the Mie_Tables program.

  use MLSKinds, only: R8

  implicit NONE
  private
  public :: Destroy_Mie, Read_Mie

  ! Coordinates
  real(r8), public, allocatable :: F_s(:), IWC_s(:), T_s(:), Theta_s(:)
  ! Phase
  real(r8), public, allocatable :: P(:,:,:,:)       ! T X IWC X Theta X F
  ! Derivatives
  real(r8), public, allocatable :: dP_dIWC(:,:,:,:) ! T X IWC X Theta X F
  real(r8), public, allocatable :: dP_dT(:,:,:,:)   ! T X IWC X Theta X F

!---------------------------- RCS Ident Info -------------------------------
  character(len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Destroy_Mie
  ! Deallocate the Mie tables

    use Allocate_Deallocate, only: E_dp, MEMORY_UNITS, Test_DeAllocate

    integer :: Status
    real :: S

    if ( allocated(f_s) ) then
      s = (e_dp * size(f_s) ) / MEMORY_UNITS
      deallocate ( f_s, stat=status )
      call test_deallocate ( status, moduleName, 'F_s', s )
      s = (e_dp * size(iwc_s) ) / MEMORY_UNITS
      deallocate ( iwc_s, stat=status )
      call test_deallocate ( status, moduleName, 'iwc_s', s )
      s = (e_dp * size(t_s) ) / MEMORY_UNITS
      deallocate ( t_s, stat=status )
      call test_deallocate ( status, moduleName, 't_s', s )
      s = (e_dp * size(theta_s) ) / MEMORY_UNITS
      deallocate ( theta_s, stat=status )
      call test_deallocate ( status, moduleName, 'theta_s', s )
      s = (e_dp * size(P) ) / MEMORY_UNITS
      deallocate ( P, stat=status )
      call test_deallocate ( status, moduleName, 'P', s )
      if ( allocated(dP_dIWC) ) then
        s = (e_dp * size(dP_dIWC) ) / MEMORY_UNITS
        deallocate ( dP_dIWC, stat=status )
        call test_deallocate ( status, moduleName, 'dP_dIWC', s )
        s = (e_dp * size(dP_dT) ) / MEMORY_UNITS
        deallocate ( dP_dT, stat=status )
        call test_deallocate ( status, moduleName, 'dP_dT', s )
      end if
    end if

  end subroutine Destroy_Mie

  subroutine Read_Mie ( FileName )

    use Allocate_Deallocate, only: Test_Allocate
    use IO_stuff, only: Get_Lun
    use Machine, only: IO_Error
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error

    character(len=*), intent(in) :: FileName

    integer :: Lun, Status

    integer :: N_f, N_IWC, N_T, N_theta, N_cut
    real(r8) :: R_min, R_max
    logical :: Derivs

    ! Destroy the old database
    call destroy_Mie

    ! Open the file
    call get_lun ( lun, msg=.false. )
    if ( lun < 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "No logical unit numbers available" )
    open ( unit=lun, file=filename, status='old', form='unformatted', &
      & access='sequential', iostat=status )
    if ( status /= 0 ) then
      call io_error ( "Unable to open Mie tables file ", status, filename )
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Unable to open Mie tables file " // Filename )
    end if

    ! Read some sizes and scalars
    read ( lun ) n_f, n_IWC, n_T, n_theta, r_min, r_max, n_cut, derivs

    ! Allocate and read the coordinate arrays
    allocate ( f_s(n_f), stat=status )
    call test_allocate ( status, moduleName, 'F_s', (/ 1 /), (/ n_f /) )
    allocate ( iwc_s(n_iwc), stat=status )
    call test_allocate ( status, moduleName, 'IWC_s', (/ 1 /), (/ n_iwc /) )
    allocate ( T_s(n_T), stat=status )
    call test_allocate ( status, moduleName, 'T_s', (/ 1 /), (/ n_T /) )
    allocate ( theta_s(n_theta), stat=status )
    call test_allocate ( status, moduleName, 'theta_s', (/ 1 /), (/ n_theta /) )
    read ( lun ) iwc_s, t_s, theta_s, f_s

    ! Allocate the phase array and its derivatives
    allocate ( p(n_t, n_iwc, n_theta, n_f), stat=status )
    call test_allocate ( status, moduleName, 'P', (/ 1,1,1,1 /), &
      & (/ n_t,n_iwc,n_theta,n_f /) )
    if ( derivs ) then
      allocate ( dP_dIWC(n_t, n_iwc, n_theta, n_f), stat=status )
      call test_allocate ( status, moduleName, 'dP_dIWC', (/ 1,1,1,1 /), &
        & (/ n_t,n_iwc,n_theta,n_f /) )
      allocate ( dP_dT(n_t, n_iwc, n_theta, n_f), stat=status )
      call test_allocate ( status, moduleName, 'dP_dT', (/ 1,1,1,1 /), &
        & (/ n_t,n_iwc,n_theta,n_f /) )
    end if

    ! Read everything
    call Read_Mie_Auto
    if ( derivs ) call Read_Mie_Derivs_Auto

    close ( lun )

  contains

    ! These are internal subroutines to avoid a bunch of explicit allocations
    subroutine Read_Mie_Auto

      integer, parameter :: I_c_e = 1, I_c_s = 2

!     real(r8) :: Beta(n_t,n_iwc,n_f,2) ! Final dimension: 1 = Beta(c_e), 2 = Beta(c_s)
!     real(r8) :: Eest(n_t,n_iwc,n_f,2) ! Error estimate
!     integer :: NFunc(2,n_t,n_iwc,n_f,2)
!     integer :: MaxOrd(n_t,n_iwc,n_f,2)
      real(r8) :: E_P(n_t,n_iwc,n_theta,n_f) ! Error estimate
      integer :: NFuncP(2, n_t, n_iwc, n_theta, n_f, 3)
      integer :: MaxOrdP(n_t, n_iwc, n_theta, n_f, 3)

!     read ( lun ) beta(:,:,:,i_c_e), eest(:,:,:,i_c_e), nFunc(:,:,:,:,i_c_e), maxOrd(:,:,:,i_c_e)
!     read ( lun ) beta(:,:,:,i_c_s), eest(:,:,:,i_c_s), nFunc(:,:,:,:,i_c_s), maxOrd(:,:,:,i_c_s)
      read ( lun ) ! Skip beta(:,:,:,i_c_e) etc
      read ( lun ) ! Skip beta(:,:,:,i_c_s) etc
      read ( lun ) p, e_p, nFuncP(:,:,:,:,:,1), maxOrdP(:,:,:,:,1)
    end subroutine Read_Mie_Auto

    subroutine Read_Mie_Derivs_Auto

      integer, parameter :: I_c_e = 1, I_c_s = 2

!     real(r8) :: dBeta_dIWC(n_t, n_iwc, n_f, 2)
!     real(r8) :: E_dBeta_dIWC(n_t, n_iwc, n_f, 2) ! Error estimate
!     real(r8) :: dBeta_dT(n_t, n_iwc, n_f, 2)
!     real(r8) :: E_dBeta_dT(n_t, n_iwc, n_f, 2) ! Error estimate
!     integer :: NFunc(2,n_t,n_iwc,n_f,3:6)
!     integer :: MaxOrd(n_t,n_iwc,n_f,3:6)
      real(r8) :: e_dP_dIWC(n_t,n_iwc,n_theta,n_f) ! Error estimate
      integer :: NFuncP(2,n_t,n_iwc,n_theta,n_f,3)
      integer :: MaxOrdP(n_t,n_iwc,n_theta,n_f,3)
      real(r8) :: e_dP_dT(n_t,n_iwc,n_theta,n_f) ! Error estimate

!     read ( lun ) dBeta_dIWC(:,:,:,i_c_e),  e_dBeta_dIWC(:,:,:,i_c_e), nFunc(:,:,:,:,3), maxOrd(:,:,:,3)
!     read ( lun ) dBeta_dIWC(:,:,:,i_c_s),  e_dBeta_dIWC(:,:,:,i_c_s), nFunc(:,:,:,:,4), maxOrd(:,:,:,4)
!     read ( lun ) dBeta_dT(:,:,:,i_c_e),  e_dBeta_dT(:,:,:,i_c_e), nFunc(:,:,:,:,5), maxOrd(:,:,:,5)
!     read ( lun ) dBeta_dT(:,:,:,i_c_s),  e_dBeta_dT(:,:,:,i_c_s), nFunc(:,:,:,:,6), maxOrd(:,:,:,6)
      read ( lun ) ! skip dBeta_dIWC(:,:,:,i_c_e) etc
      read ( lun ) ! skip dBeta_dIWC(:,:,:,i_c_s) etc
      read ( lun ) ! skip dBeta_dT(:,:,:,i_c_e) etc
      read ( lun ) ! skip dBeta_dT(:,:,:,i_c_s) etc
      read ( lun ) dP_dIWC, e_dP_dIWC, nFuncP(:,:,:,:,:,2), maxOrdP(:,:,:,:,2)
      read ( lun ) dP_dT, e_dP_dT, nFuncP(:,:,:,:,:,3), maxOrdP(:,:,:,:,3)
    end subroutine Read_Mie_Derivs_Auto

  end subroutine Read_Mie

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Read_Mie_m

! $Log$
! Revision 2.1  2008/05/20 00:26:40  vsnyder
! Initial commit
!
