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

  use MLSKinds, only: R8, RP

  implicit NONE
  private
  public :: Destroy_Log_Mie, Destroy_Mie, Dump_Mie, Log_Mie, Read_Mie

  ! Coordinates
  real(rp), public, target, allocatable :: F_s(:)     ! MHz
  real(rp), public, target, allocatable :: IWC_s(:)   ! log10(IWC)
  real(rp), public, target, allocatable :: T_s(:)     ! K
  real(rp), public, target, allocatable :: Theta_s(:) ! Radians above horizon
  ! Betas
  real(r8), public, target, allocatable :: Beta_c_a(:,:,:)      ! T X IWC X F
  real(r8), public, target, allocatable :: Beta_c_e(:,:,:)      ! T X IWC X F
  real(r8), public, target, allocatable :: Beta_c_s(:,:,:)      ! T X IWC X F
  real(r8), public, target, allocatable :: Log_Beta_c_a(:,:,:)  ! T X IWC X F
  real(r8), public, target, allocatable :: Log_Beta_c_e(:,:,:)  ! T X IWC X F
  real(r8), public, target, allocatable :: Log_Beta_c_s(:,:,:)  ! T X IWC X F
  ! Beta derivatives
  real(r8), public, target, allocatable :: dBeta_dIWC_c_a(:,:,:)     ! T X IWC X F
  real(r8), public, target, allocatable :: dBeta_dT_c_a(:,:,:)       ! T X IWC X F
  real(r8), public, target, allocatable :: dBeta_dIWC_c_e(:,:,:)     ! T X IWC X F
  real(r8), public, target, allocatable :: dBeta_dT_c_e(:,:,:)       ! T X IWC X F
  real(r8), public, target, allocatable :: dBeta_dIWC_c_s(:,:,:)     ! T X IWC X F
  real(r8), public, target, allocatable :: dBeta_dT_c_s(:,:,:)       ! T X IWC X F
  ! Phase
  real(r8), public, target, allocatable :: P(:,:,:,:)       ! T X IWC X Theta X F
  ! Phase Derivatives
  real(r8), public, target, allocatable :: dP_dIWC(:,:,:,:) ! T X IWC X Theta X F
  real(r8), public, target, allocatable :: dP_dT(:,:,:,:)   ! T X IWC X Theta X F

!---------------------------- RCS Ident Info -------------------------------
  character(len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! --------------------------------------------- Destroy_Log_Mie  -----
  subroutine Destroy_Log_Mie
  ! Deallocate the Log Mie tables

    use Allocate_Deallocate, only: Byte_Size, Test_DeAllocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc

    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: Status
    integer :: S

    if ( allocated(log_beta_c_a) ) then
      s = byte_size(log_beta_c_a)
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(log_beta_c_a(1,1,1)), addr)
      deallocate ( log_beta_c_a, stat=status )
      call test_deallocate ( status, moduleName, 'log_beta_c_a', s, address=addr )
    end if
    if ( allocated(log_beta_c_e) ) then
      s = byte_size(log_beta_c_e)
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(log_beta_c_e(1,1,1)), addr)
      deallocate ( log_beta_c_e, stat=status )
      call test_deallocate ( status, moduleName, 'log_beta_c_e', s, address=addr )
    end if
    if ( allocated(log_beta_c_s) ) then
      s = byte_size(log_beta_c_s)
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(log_beta_c_s(1,1,1)), addr)
      deallocate ( log_beta_c_s, stat=status )
      call test_deallocate ( status, moduleName, 'log_beta_c_s', s, address=addr )
    end if
  end subroutine Destroy_Log_Mie
  
  ! ------------------------------------------------  Destroy_Mie  -----
  subroutine Destroy_Mie
  ! Deallocate the Mie tables

    use Allocate_Deallocate, only: Byte_Size, Test_DeAllocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc

    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: Status
    integer :: S

    if ( allocated(f_s) ) then
      s = byte_size(f_s)
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(f_s(1)), addr)
      deallocate ( f_s, stat=status )
      call test_deallocate ( status, moduleName, 'F_s', s, address=addr )
      s = byte_size(iwc_s)
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(iwc_s(1)), addr)
      deallocate ( iwc_s, stat=status )
      call test_deallocate ( status, moduleName, 'iwc_s', s, address=addr )
      s = byte_size(t_s)
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(t_s(1)), addr)
      deallocate ( t_s, stat=status )
      call test_deallocate ( status, moduleName, 't_s', s, address=addr )
      s = byte_size(theta_s)
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(theta_s(1)), addr)
      deallocate ( theta_s, stat=status )
      call test_deallocate ( status, moduleName, 'theta_s', s, address=addr )
      s = byte_size(beta_c_a)
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(beta_c_a(1,1,1)), addr)
      deallocate ( beta_c_a, stat=status )
      call test_deallocate ( status, moduleName, 'beta_c_a', s, address=addr )
      s = byte_size(beta_c_e)
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(beta_c_e(1,1,1)), addr)
      deallocate ( beta_c_e, stat=status )
      call test_deallocate ( status, moduleName, 'beta_c_e', s, address=addr )
      s = byte_size(beta_c_s)
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(beta_c_s(1,1,1)), addr)
      deallocate ( beta_c_s, stat=status )
      call test_deallocate ( status, moduleName, 'beta_c_s', s, address=addr )
      s = byte_size(P)
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(P(1,1,1,1)), addr)
      deallocate ( P, stat=status )
      call test_deallocate ( status, moduleName, 'P', s, address=addr )
      if ( allocated(dP_dIWC) ) then
        s = byte_size(dBeta_dIWC_c_a)
        addr = 0
        if ( s > 0 ) addr = transfer(c_loc(dBeta_dIWC_c_a(1,1,1)), addr)
        deallocate ( dBeta_dIWC_c_a, stat=status )
        call test_deallocate ( status, moduleName, 'dBeta_dIWC_c_a', s, address=addr )
        s = byte_size(dBeta_dT_c_a)
        addr = 0
        if ( s > 0 ) addr = transfer(c_loc(dBeta_dT_c_a(1,1,1)), addr)
        deallocate ( dBeta_dT_c_a, stat=status )
        call test_deallocate ( status, moduleName, 'dBeta_dT_c_a', s, address=addr )
        s = byte_size(dBeta_dIWC_c_e)
        addr = 0
        if ( s > 0 ) addr = transfer(c_loc(dBeta_dIWC_c_e(1,1,1)), addr)
        deallocate ( dBeta_dIWC_c_e, stat=status )
        call test_deallocate ( status, moduleName, 'dBeta_dIWC_c_e', s, address=addr )
        s = byte_size(dBeta_dT_c_e)
        addr = 0
        if ( s > 0 ) addr = transfer(c_loc(dBeta_dT_c_e(1,1,1)), addr)
        deallocate ( dBeta_dT_c_e, stat=status )
        call test_deallocate ( status, moduleName, 'dBeta_dT_c_e', s, address=addr )
        s = byte_size(dBeta_dIWC_c_s)
        addr = 0
        if ( s > 0 ) addr = transfer(c_loc(dBeta_dIWC_c_s(1,1,1)), addr)
        deallocate ( dBeta_dIWC_c_s, stat=status )
        call test_deallocate ( status, moduleName, 'dBeta_dIWC_c_s', s, address=addr )
        s = byte_size(dBeta_dT_c_s)
        addr = 0
        if ( s > 0 ) addr = transfer(c_loc(dBeta_dT_c_s(1,1,1)), addr)
        deallocate ( dBeta_dT_c_s, stat=status )
        call test_deallocate ( status, moduleName, 'dBeta_dT_c_s', s, address=addr )
        s = byte_size(dP_dIWC)
        addr = 0
        if ( s > 0 ) addr = transfer(c_loc(dP_dIWC(1,1,1,1)), addr)
        deallocate ( dP_dIWC, stat=status )
        call test_deallocate ( status, moduleName, 'dP_dIWC', s, address=addr )
        s = byte_size(dP_dT)
        addr = 0
        if ( s > 0 ) addr = transfer(c_loc(dP_dT(1,1,1,1)), addr)
        deallocate ( dP_dT, stat=status )
        call test_deallocate ( status, moduleName, 'dP_dT', s, address=addr )
      end if
    end if

    call destroy_log_mie

  end subroutine Destroy_Mie

  ! ---------------------------------------------------  Dump_Mie  -----
  subroutine Dump_Mie ( Details )
  ! Dump the Mie tables
  ! Details <= 0 => F_s, T_s, IWC_s, Theta_s
  !          = 1 => 0 + P
  !         >= 2 => 1 + dP_dIWC + dP_dT
    use Dump_0, only: Dump
    use Output_m, only: Output
    use Constants, only: Rad2Deg
    integer, intent(in) :: Details
    integer :: I_f, I_theta
    if ( .not. allocated(f_s) ) return
    call dump ( f_s, name='F_s (MHz)' )
    call dump ( T_s, name='T_s (K)' )
    call dump ( IWC_s, name='IWC_s (actually log10 IWC)' )
    call dump ( Rad2Deg*theta_s, name='Theta_s (Degrees)' )
    if ( details <= 0 ) return
    do i_f = 1, size(f_s)
      call output ( f_s(i_f), before='Beta_c_a(T X IWC), F = ', advance='yes' )
      call dump ( beta_c_a(:,:,i_f) )
      call output ( f_s(i_f), before='Beta_c_e(T X IWC), F = ', advance='yes' )
      call dump ( beta_c_e(:,:,i_f) )
      if ( allocated(log_beta_c_a) ) then
        call output ( f_s(i_f), before='log_beta_c_a(T X IWC), F = ', advance='yes' )
        call dump ( log_beta_c_a(:,:,i_f) )
      end if
      if ( allocated(log_beta_c_e) ) then
        call output ( f_s(i_f), before='Log_Beta_c_e(T X IWC), F = ', advance='yes' )
        call dump ( log_beta_c_e(:,:,i_f) )
      end if
      call output ( f_s(i_f), before='Beta_c_s(T X IWC), F = ', advance='yes' )
      call dump ( beta_c_s(:,:,i_f) )
      if ( allocated(log_beta_c_s) ) then
        call output ( f_s(i_f), before='Log_Beta_c_s(T X IWC), F = ', advance='yes' )
        call dump ( log_beta_c_s(:,:,i_f) )
      end if
      if ( allocated(dBeta_dIWC_c_a) .and. details >= 2 ) then
        call output ( f_s(i_f), before='dBeta_dIWC_c_a(T X IWC), F = ', advance='yes' )
        call dump ( dBeta_dIWC_c_a(:,:,i_f) )
        call output ( f_s(i_f), before='dBeta_dT_c_a(T X IWC), F = ', advance='yes' )
        call dump ( dBeta_dT_c_a(:,:,i_f) )
        call output ( f_s(i_f), before='dBeta_dIWC_c_e(T X IWC), F = ', advance='yes' )
        call dump ( dBeta_dIWC_c_e(:,:,i_f) )
        call output ( f_s(i_f), before='dBeta_dT_c_e(T X IWC), F = ', advance='yes' )
        call dump ( dBeta_dT_c_e(:,:,i_f) )
        call output ( f_s(i_f), before='dBeta_dIWC_c_s(T X IWC), F = ', advance='yes' )
        call dump ( dBeta_dIWC_c_s(:,:,i_f) )
        call output ( f_s(i_f), before='dBeta_dT_c_s(T X IWC), F = ', advance='yes' )
        call dump ( dBeta_dT_c_s(:,:,i_f) )
      end if
      do i_theta = 1, size(theta_s)
        call output ( Rad2Deg*theta_s(i_theta), before='P(T X IWC,' )
        call output ( f_s(i_f), before=') F = ', advance='yes' )
        call dump ( p(:,:,i_theta,i_f) )
        if ( allocated(dP_dIWC) .and. details >= 2 ) then
          call output ( Rad2Deg*theta_s(i_theta), before='dP_dIWC(T X IWC,' )
          call output ( f_s(i_f), before=') F = ', advance='yes' )
          call dump ( dP_dIWC(:,:,i_theta,i_f) )
          call output ( Rad2Deg*theta_s(i_theta), before='dP_dT(T X IWC,' )
          call output ( f_s(i_f), before=') F = ', advance='yes' )
          call dump ( dP_dT(:,:,i_theta,i_f) )
        end if
      end do
    end do
  end subroutine Dump_Mie

  ! ----------------------------------------------------  Log_Mie  -----
  subroutine Log_Mie
  ! If any of Log_Beta_c_a, Log_Beta_c_e, Log_Beta_c_s are not associated,
  ! allocate them and compute them.  Read_Mie has to be called first.

    use Allocate_Deallocate, only: Bytes, Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: Status

    if ( .not. allocated(beta_c_e) ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & 'Cannot compute log(beta...) before allocating beta.')

    if ( .not. allocated(log_beta_c_a) ) then
      allocate ( &
        & log_beta_c_a(ubound(beta_c_e,1),ubound(beta_c_e,2),ubound(beta_c_e,3)), &
        & stat=status )
      addr = 0
      if ( status == 0 ) then
        if ( size(log_beta_c_a) > 0 ) addr = transfer(c_loc(log_beta_c_a(1,1,1)), addr)
      end if
      call test_allocate ( status, moduleName, 'log_beta_c_a', (/ 1,1,1 /), &
          & ubound(beta_c_e), bytes(beta_c_e), address=addr )
      log_beta_c_a = log(beta_c_e)
    end if
    if ( .not. allocated(log_beta_c_e) ) then
      allocate ( &
        & log_beta_c_e(ubound(beta_c_e,1),ubound(beta_c_e,2),ubound(beta_c_e,3)), &
        & stat=status )
      addr = 0
      if ( status == 0 ) then
        if ( size(log_beta_c_e) > 0 ) addr = transfer(c_loc(log_beta_c_e(1,1,1)), addr)
      end if
      call test_allocate ( status, moduleName, 'Log_Beta_c_e', (/ 1,1,1 /), &
          & ubound(beta_c_e), bytes(log_beta_c_e), address=addr )
      log_beta_c_e = log(beta_c_e)
    end if
    if ( .not. allocated(log_beta_c_s) ) then
      allocate ( &
        & log_beta_c_s(ubound(beta_c_s,1),ubound(beta_c_s,2),ubound(beta_c_s,3)), &
        & stat=status )
      addr = 0
      if ( status == 0 ) then
        if ( size(log_beta_c_s) > 0 ) addr = transfer(c_loc(log_beta_c_s(1,1,1)), addr)
      end if
      call test_allocate ( status, moduleName, 'Log_Beta_c_s', (/ 1,1,1 /), &
          & ubound(beta_c_s), bytes(beta_c_s), address=addr )
      log_beta_c_s = log(beta_c_s)
    end if
  end subroutine Log_Mie

  ! ---------------------------------------------------  Read_Mie  -----
  subroutine Read_Mie ( FileName )

    use Allocate_Deallocate, only: Bytes, Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    use Machine, only: IO_Error
    use MLSHDF5, only: IsHDF5DSPresent, loadAllocFromHDF5DS
    use HDF5, only: H5FClose_f, H5FOpen_f, H5F_ACC_RDONLY_F
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSStrings, only: Capitalize

    character(len=*), intent(in) :: FileName

    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: FileID, Lun, Status

    integer :: N_f, N_IWC, N_T, N_theta, N_cut
    real(r8) :: R_min, R_max
    logical :: Derivs

    ! Destroy the old database
    call destroy_Mie

    if ( index(capitalize(fileName),'.HDF') /= 0 ) then
      ! Read tables from HDF5 file
      call H5FOpen_f ( trim(fileName), H5F_ACC_RDONLY_F, fileID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Unable to open Mie tables file " // trim(fileName) )
      ! Get the index arrays
      call loadAllocFromHDF5DS ( fileID, 'IWC_s', iwc_s )
      call loadAllocFromHDF5DS ( fileID, 'T_s', T_s )
      call loadAllocFromHDF5DS ( fileID, 'Theta_s', theta_s )
      call loadAllocFromHDF5DS ( fileID, 'F_s', f_s )
      f_s = 1000.0 * f_s ! F_s are GHz, standard unit is MHz
      ! Get the betas
      call loadAllocFromHDF5DS ( fileID, 'Beta_c_a', Beta_c_a )
      call loadAllocFromHDF5DS ( fileID, 'Beta_c_e', Beta_c_e )
      call loadAllocFromHDF5DS ( fileID, 'Beta_c_s', Beta_c_s )
      ! Get the beta derivatives
      call loadAllocFromHDF5DS ( fileID, 'dBeta_dIWC_c_a', dBeta_dIWC_c_a )
      call loadAllocFromHDF5DS ( fileID, 'dBeta_dT_c_a', dBeta_dT_c_a )
      call loadAllocFromHDF5DS ( fileID, 'dBeta_dIWC_c_e', dBeta_dIWC_c_e )
      call loadAllocFromHDF5DS ( fileID, 'dBeta_dT_c_e', dBeta_dT_c_e )
      call loadAllocFromHDF5DS ( fileID, 'dBeta_dIWC_c_s', dBeta_dIWC_c_s )
      call loadAllocFromHDF5DS ( fileID, 'dBeta_dT_c_s', dBeta_dT_c_s )
      ! Get the phase function
      call loadAllocFromHDF5DS ( fileID, 'P', P )
      ! Get the phase function derivatives
      if ( IsHDF5DSPresent ( fileID, 'dP_dIWC' ) ) then
        call loadAllocFromHDF5DS ( fileID, 'dP_dIWC', dP_dIWC )
        call loadAllocFromHDF5DS ( fileID, 'dP_dT', dP_dT )
      end if
      call H5FClose_f ( fileID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Unable to close Mie tables file " // trim(fileName) )
    else
      ! Read tables from Fortran unformatted file
      ! Open the file
      open ( newunit=lun, file=filename, status='old', form='unformatted', &
        & access='sequential', iostat=status )
      if ( status /= 0 ) then
        call io_error ( "Unable to open Mie tables file ", status, filename )
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & "Unable to open Mie tables file " // trim(fileName) )
      end if

      ! Read some sizes and scalars
      read ( lun ) n_f, n_IWC, n_T, n_theta, r_min, r_max, n_cut, derivs

      ! Allocate and read the coordinate arrays
      allocate ( f_s(n_f), stat=status )
      addr = 0
      if ( status == 0 .and. n_f > 0 ) addr = transfer(c_loc(f_s(1)), addr)
      call test_allocate ( status, moduleName, 'F_s', (/ 1 /), (/ n_f /), &
        & bytes(f_s), address=addr )
      allocate ( iwc_s(n_iwc), stat=status )
      addr = 0
      if ( status == 0 .and. n_iwc > 0 ) addr = transfer(c_loc(iwc_s(1)), addr)
      call test_allocate ( status, moduleName, 'IWC_s', (/ 1 /), (/ n_iwc /), &
        & bytes(iwc_s), address=addr )
      allocate ( T_s(n_T), stat=status )
      addr = 0
      if ( status == 0 .and. n_T > 0 ) addr = transfer(c_loc(T_s(1)), addr)
      call test_allocate ( status, moduleName, 'T_s', (/ 1 /), (/ n_T /), &
        & bytes(T_s), address=addr )
      allocate ( theta_s(n_theta), stat=status )
      addr = 0
      if ( status == 0 .and. n_theta > 0 ) addr = transfer(c_loc(theta_s(1)), addr)
      call test_allocate ( status, moduleName, 'theta_s', (/ 1 /), (/ n_theta /), &
        & bytes(theta_s), address=addr )
      read ( lun ) iwc_s, t_s, theta_s, f_s

      ! Allocate the phase array and its derivatives
      allocate ( p(n_t, n_iwc, n_theta, n_f), stat=status )
      addr = 0
      if ( status == 0 ) then
        if ( size(p) > 0 ) addr = transfer(c_loc(p(1,1,1,1)), addr)
      end if
      call test_allocate ( status, moduleName, 'P', (/ 1,1,1,1 /), &
        & (/ n_t,n_iwc,n_theta,n_f /), bytes(p), address=addr )
      if ( derivs ) then
        allocate ( dP_dIWC(n_t, n_iwc, n_theta, n_f), stat=status )
        addr = 0
        if ( status == 0 ) then
          if ( size(dP_dIWC) > 0 ) addr = transfer(c_loc(dP_dIWC(1,1,1,1)), addr)
        end if
        call test_allocate ( status, moduleName, 'dP_dIWC', (/ 1,1,1,1 /), &
          & (/ n_t,n_iwc,n_theta,n_f /), bytes(dP_dIWC), address=addr )
        allocate ( dP_dT(n_t, n_iwc, n_theta, n_f), stat=status )
        addr = 0
        if ( status == 0 ) then
          if ( size(dP_dT) > 0 ) addr = transfer(c_loc(dP_dT(1,1,1,1)), addr)
        end if
        call test_allocate ( status, moduleName, 'dP_dT', (/ 1,1,1,1 /), &
          & (/ n_t,n_iwc,n_theta,n_f /), bytes(dP_dT), address=addr )
      end if

      ! Read everything
      call Read_Mie_Auto
      if ( derivs ) call Read_Mie_Derivs_Auto

      close ( lun )
    end if

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

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Read_Mie_m

! $Log$
! Revision 2.14  2018/05/15 03:26:25  vsnyder
! Change Mie tables from pointer to allocatable
!
! Revision 2.13  2015/03/28 02:03:48  vsnyder
! Added stuff to trace allocate/deallocate addresses.  Use NewUnit=
! specifier in OPEN statement instead of Get_Lun.
!
! Revision 2.12  2014/09/05 20:52:32  vsnyder
! More complete and accurate allocate/deallocate size tracking
!
! Revision 2.11  2013/06/12 02:21:57  vsnyder
! Use BYTE_SIZE from Allocate_Deallocate
!
! Revision 2.10  2011/01/26 02:50:43  vsnyder
! Add more support for beta_c_a separately from beta_c_e
!
! Revision 2.9  2010/06/07 23:13:02  vsnyder
! Added ability to compute and store Log_Beta_c_e and Log_Beta_c_s
!
! Revision 2.8  2010/01/14 02:24:06  vsnyder
! Add some comments
!
! Revision 2.7  2009/08/20 19:47:26  vsnyder
! Replace THETA_s by Theta_s to conform to Mie_Tables
!
! Revision 2.6  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.5  2009/05/13 20:03:01  vsnyder
! Get constants from Constants, kinds from MLSKinds
!
! Revision 2.4  2008/10/20 23:23:20  vsnyder
! Convert F_s from GHz to MHz, since MHz is the standard unit
!
! Revision 2.3  2008/07/31 18:00:13  vsnyder
! Change coordinates from R8 to RP
!
! Revision 2.2  2008/06/05 02:15:26  vsnyder
! Added HDF, spiffed up dump
!
! Revision 2.1  2008/05/20 00:26:40  vsnyder
! Initial commit
!
