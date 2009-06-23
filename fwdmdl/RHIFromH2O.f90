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
module RHIFromH2O                     ! H2O <-> RHI Conversions
  !=============================================================================

  ! This module gathers into one place subroutines converting
  ! between H2O  concentration and RHI.
  ! Also their respective precisions

  use MLSCommon, only: R8
  implicit none
  private
  public :: RHIFromH2O_Factor, RHIPrecFromH2O
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
  real, parameter ::    UNDEFINED_VALUE = -999.99 ! Same as %template%badvalue

contains ! =====     Public Procedures     =============================

  !--------------------------------------------  RHIFromH2O_Factor  -----

  function RHIFromH2O_Factor ( T, zeta, vmr_unit_cnv, invert)
  ! Factor relating rhi to H2O vmr
  ! (See Eq. 9 from "UARS Microwave Limb Sounder upper tropospheric
  !  humidity measurement: Method and validation" Read et. al. 
  !  J. Geoph. Res. Dec. 2001 (106) D23)
  ! Arguments
  real(r8), intent(in) :: T               ! Temperature
  real(r8), intent(in) :: zeta            ! Surface in log pressure units
  integer, intent(in)  :: vmr_unit_cnv   ! E.g., 6 for ppmv
  logical, optional, intent(in) :: invert ! If TRUE, H2O = factor RHI
                                          ! else, RHI = factor H2O (default)
  real(r8)             :: RHIFromH2O_Factor
  ! Local variables
  integer :: invs
  
  ! Executable
  invs = -1
  if ( present(invert) ) then
    if ( invert ) invs = 1
  endif
  RHIFromH2O_Factor = &
     & exp(invs*( &
     & (C(T)+zeta+vmr_unit_cnv) * log(10.) &
     & + &
     & 3.56654*log(T/273.16) &
     & ))
  
  end function RHIFromH2O_Factor

  !--------------------------------------------  RHIPrecFromH2O  -----

  subroutine RHIPrecFromH2O ( H2Ovmr, T, zeta, vmr_unit_cnv, &
    & H2OPrecision, TPrecision, RHIPrecision, NEGATIVETOO )
  ! Calculate RHI Precision based on that of H2O and Temperature
  ! Eq. 9 from "UARS Microwave Limb Sounder upper tropospheric
  !  humidity measurement: Method and validation" Read et. al. 
  !  J. Geoph. Res. Dec. 2001 (106) D23
  ! recoded by Mark Filipiak
  ! Arguments
  real(r8), intent(in)  :: H2Ovmr          ! vmr of H2O
  real(r8), intent(in)  :: T               ! Temperature
  real(r8), intent(in)  :: zeta            ! Surface in log pressure units
  integer, intent(in)   :: vmr_unit_cnv    ! E.g., 6 for ppmv
  real(r8), intent(in)  :: H2OPrecision    ! Precision of H2O
  real(r8), intent(in)  :: TPrecision      ! Precision of Temperature
  real(r8), intent(out) :: RHIPrecision    ! Precision of RHI
  logical, intent(in), optional   :: NEGATIVETOO  ! Set RHI Precision negative 
  ! Local variables                          if either T or H2O Precsisions are
  integer :: invs
  real(r8) :: df_db       ! RHi deriv wrt H2O
  real(r8) :: df_dT       ! RHi deriv wrt T
  logical :: isNegative
  
  ! Executable
  isNegative = .false.
  if ( present(negativeToo) ) &
    & isNegative = negativeToo .and. &
    & ( H2OPrecision < 0.d0 .or. TPrecision < 0.d0 )
  invs = -1
  df_db = exp(invs*( &
   & (C(T)+zeta+vmr_unit_cnv) * log(10.) &
   & + &
   & 3.56654*log(T/273.16) &
   & ))
  df_dT = H2Ovmr * exp(invs*( &
   & (C(T)+zeta+vmr_unit_cnv) * log(10.) &
   & + &
   & 3.56654*log(T/273.16) &
   & )) &
   & * invs * ( dC_dT(T) * log(10.) + 3.56654 / T )
  RHIPrecision = sqrt (&
   & ( H2OPrecision * df_db )**2 &
   & + ( TPrecision * df_dT )**2 &
   & )
  if ( isNegative ) RHIPrecision = -RHIPrecision
  contains
    function dC_dT ( T )
      ! As found in ref.
      real(r8), intent(in)   :: T
      real(r8)               :: dC_dT
      ! Local
      real(r8), parameter    :: a0 = -1.2141649d0
      real(r8), parameter    :: a1 = 9.09718d0
      real(r8), parameter    :: a2 = 0.876793d0
      real, parameter        :: ILLEGALTEMP = UNDEFINED_VALUE
      !
      if ( T > 0.d0 ) then
        dC_dT = a1*(273.16/T**2) - a2/273.16
      else
        dC_dT = ILLEGALTEMP
      end if
    end function dC_dT

  end subroutine RHIPrecFromH2O

! =====     Private Procedures     =============================
    function C ( T )
      ! As found in ref.
      real(r8), intent(in)   :: T
      real(r8)               :: C
      ! Local
      real(r8), parameter    :: a0 = -1.2141649d0
      real(r8), parameter    :: a1 = 9.09718d0
      real(r8), parameter    :: a2 = 0.876793d0
      real, parameter        :: ILLEGALTEMP = UNDEFINED_VALUE
      !
      if ( T > 0.d0 ) then
        C = a0 - a1*(273.16/T -1.0d0) + a2*(1.0d0 - T/273.16)
      else
        C = ILLEGALTEMP
      end if
    end function C

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module RHIFromH2O
!=============================================================================

!
! $Log$
! Revision 2.4  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.3  2005/10/17 17:25:48  pwagner
! May set RHIPrecision negative if either T or H2O is
!
! Revision 2.2  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.1  2003/01/28 21:53:07  pwagner
! RHI H2O conversions moved to fwdmdl from l2/Fill
!
