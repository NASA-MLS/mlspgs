module D_FFT
!>> 2000-04-24 D_FFT WV Snyder Converted to Fortran 90
!--D replaces "?": ?FFT, ?RFT1, '?'
  use ERMSG_M, only: ERMSG, ERM1
  implicit NONE
  private
  public :: CFT, FFT, RFT, RFT1, TCST      ! Generic
  public :: DCFT, DFFT, DRFT, DRFT1, DTCST ! Specific
  interface CFT;  module procedure DCFT;  end interface
  interface FFT;  module procedure DFFT;  end interface
  interface RFT;  module procedure DRFT;  end interface
  interface RFT1; module procedure DRFT1; end interface
  interface TCST; module procedure DTCST; end interface

! -----     Private declarations     -----------------------------------
  integer, private, parameter :: RK = kind(1.0d0)
  character, private, parameter :: PREC = 'D' ! For error messages
  integer, private, parameter :: KEDIM=30
  logical, save :: NEEDST       ! .TRUE. if the sine table must be computed.
  integer, save :: ILAST        ! KS * 2**MM
  integer, save :: KE(KEDIM)    ! (/ KS * 2**(MM-L), L=1, MM /)
  integer, save :: KEE(KEDIM+1) ! Equivalenced to KE -- see below
  integer, save :: KS           ! distance in memory between successive
                                ! points.  The i-th coefficient, a(i), is
                                ! given by AR((I+1)*KS)+AI((I+1)*KS)*
                                ! sqrt(-1), i=0, 1, ..., (2**MM)-1.
  integer, save :: MM           ! base 2 log(number of complex fourier
                                ! coefficients)
  integer, save :: MT           ! base 2 log(NT)
  integer, save :: N1           ! Equivalenced to ILAST
  integer, save :: NT           ! number of entries in the sine table
  equivalence (ILAST, KEE(1), N1)
  equivalence (KE(1), KEE(2))

  real(rk), parameter :: HALF = 0.5_rk
!     PI4 = PI / 4
  real(rk), parameter :: PI4 = 0.7853981633974483096156608458198757_rk
!     SPI4 = SIN(PI/4) = .5 * SQRT(2)
  real(rk), parameter :: SPI4 = 0.7071067811865475244008443621048490_rk
  real(rk), parameter :: TWO = 2.0_rk

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains ! -------------------------------------------------------------


  subroutine DCFT (A, MODE, M, ND, MS, S)
    include 'cft.f9h'
  end subroutine DCFT
  subroutine DFFT ( AR, AI, S )
    include 'fft.f9h'
  end subroutine DFFT
  subroutine DRFT (A, MODE, M, ND, MS, S)
    include 'rft.f9h'
  end subroutine DRFT
  subroutine DRFT1 (A, MODE, M, MS, S)
    include 'rft1.f9h'
  end subroutine DRFT1
  subroutine DTCST (A, TCS, MODE, M, ND, MS, S)
    include 'tcst.f9h'
  end subroutine DTCST

end module D_FFT

! $Log$
