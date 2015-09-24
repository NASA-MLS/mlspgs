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

  use MLSCommon, only: defaultUndefinedValue
  use MLSKinds, only: r8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Warning
  implicit none
  private
  public :: RHIFROMH2O_FACTOR, RHIPRECFROMH2O, H2OPRECFROMRHI
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  real(r8), parameter :: TooBig = 300.
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
  real(r8)             :: itsLog
  ! Local variables
  integer :: invs
  
  ! Executable
  invs = -1
  if ( present(invert) ) then
    if ( invert ) invs = 1
  endif
  RHIFromH2O_Factor = defaultUndefinedValue
  if ( C(T) /= defaultUndefinedValue ) then
    itsLog = (C(T)+zeta+vmr_unit_cnv) * log(10.) &
      & + &
      & 3.56654*log(T/273.16) 
    if ( abs(itsLog) < TooBig ) then
      RHIFromH2O_Factor = &
      & exp(invs*itsLog )
    else
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & "Cannot recover RHi from H2O" )
    endif
  endif
  
  end function RHIFromH2O_Factor

  !--------------------------------------------  RHIPrecFromH2O  -----
  ! Idea:
  ! Assuming the errors in {H2O] and Temperature are uncorrelated, the errors
  ! in a function jointly of x and y can be expressed as
  ! (d f)^2 = (d x (ds f / ds x)^2 + (d y (ds f / ds y)^2
  
  ! d[RHi]^2 = ( d [H2O] (ds [RHi]/ ds [H2O]) )^2 + ( [H2O] dT (ds F / ds T) )^2
  ! where ds y / ds x means the (partial) derivative of y w.r.t. x
  ! d x means precision or uncertainty in x
  ! and F(T) is the factor computed in RHIFromH2O_Factor
  ! namely, [RHi] = F(T) [H2O]; so finally
  ! d[H2O]^2 = ( d [H2O] F )^2 + ( [H2O] dT (ds F / ds T) )^2
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
  real(r8)             :: itsLog
  
  ! Executable
  isNegative = .false.
  if ( present(negativeToo) ) &
    & isNegative = negativeToo .and. &
    & ( H2OPrecision < 0.d0 .or. TPrecision < 0.d0 )
  invs = -1
  RHIPrecision = defaultUndefinedValue
  itsLog = (C(T)+zeta+vmr_unit_cnv) * log(10.) &
   & + &
   & 3.56654*log(T/273.16)
  if ( abs(itsLog) > TooBig ) then
    call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & "Cannot recover RHi Precision from H2O" )
    return
  endif
  df_db = exp(invs*( &
   & itsLog  ))
  df_dT = H2Ovmr * exp(invs*( &
   & itsLog  )) &
   & * invs * ( dC_dT(T) * log(10.) + 3.56654 / T )
  RHIPrecision = sqrt (&
   & ( H2OPrecision * df_db )**2 &
   & + ( TPrecision * df_dT )**2 &
   & )
  if ( isNegative ) RHIPrecision = -RHIPrecision
  end subroutine RHIPrecFromH2O

  !--------------------------------------------  H2OPrecFromRhI  -----
  ! Idea:
  ! Assuming the errors in {RHi] and Temperature are uncorrelated, the errors
  ! in a function jointly of x and y can be expressed as
  ! (d f)^2 = (d x (ds f / ds x)^2 + (d y (ds f / ds y)^2
  
  ! d[H2O]^2 = ( d [RHi] (ds [H2O]/ ds [RHi]) )^2 + ( [RHi] dT (ds F^-1 / ds T) )^2
  ! where ds y / ds x means the (partial) derivative of y w.r.t. x
  ! d x means precision or uncertainty in x
  ! and F(T) is the factor computed in RHIFromH2O_Factor
  ! namely, [RHi] = F(T) [H2O]; so finally
  ! d[H2O]^2 = ( d [RHi] / F )^2 + ( [H2O] dT (ds F / ds T) / F )^2
  subroutine H2OPrecFromRhI ( H2Ovmr, T, zeta, vmr_unit_cnv, &
    & RHIPrecision, TPrecision, H2OPrecision, NEGATIVETOO )
  ! Calculate H2O Precision based on that of RHi and Temperature
  ! Eq. 9 from "UARS Microwave Limb Sounder upper tropospheric
  !  humidity measurement: Method and validation" Read et. al. 
  !  J. Geoph. Res. Dec. 2001 (106) D23
  ! recoded by PAW
  ! Arguments
  real(r8), intent(in)  :: H2Ovmr          ! vmr of H2O
  real(r8), intent(in)  :: T               ! Temperature
  real(r8), intent(in)  :: zeta            ! Surface in log pressure units
  integer, intent(in)   :: vmr_unit_cnv    ! E.g., 6 for ppmv
  real(r8), intent(in)  :: TPrecision      ! Precision of Temperature
  real(r8), intent(in) :: RHIPrecision    ! Precision of RHI
  real(r8), intent(out)  :: H2OPrecision    ! Precision of H2O
  logical, intent(in), optional   :: NEGATIVETOO  ! Set RHI Precision negative 
  ! Local variables                          if either T or H2O Precsisions are
  integer :: invs
  real(r8) :: Eff         ! F in the above
  real(r8) :: df_dT       ! ds F / ds T in the above
  logical :: isNegative
  real(r8)             :: itsLog
  
  ! Executable
  isNegative = .false.
  if ( present(negativeToo) ) &
    & isNegative = negativeToo .and. &
    & ( RHIPrecision < 0.d0 .or. TPrecision < 0.d0 )
  invs = -1
  H2OPrecision = defaultUndefinedValue
  itsLog = (C(T)+zeta+vmr_unit_cnv) * log(10.) &
   & + &
   & 3.56654*log(T/273.16)
  if ( abs(itsLog) > TooBig ) then
    call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & "Cannot recover H2O Precision from RHi" )
    return
  endif
  Eff = exp(invs*( &
   & itsLog  ))
  df_dT = exp(invs*( &
   & itsLog  )) &
   & * invs * ( dC_dT(T) * log(10.) + 3.56654 / T )
  H2OPrecision = sqrt (&
   & ( RHIPrecision / Eff )**2 &
   & + ( TPrecision * H2Ovmr * df_dT / Eff )**2 &
   & )
  if ( isNegative ) H2OPrecision = -H2OPrecision
  end subroutine H2OPrecFromRhI

! =====     Private Procedures     =============================
    function C ( T )
      ! As found in ref.
      real(r8), intent(in)   :: T
      real(r8)               :: C
      ! Local
      real(r8), parameter    :: a0 = -1.2141649d0
      real(r8), parameter    :: a1 = 9.09718d0
      real(r8), parameter    :: a2 = 0.876793d0
      real, parameter        :: ILLEGALTEMP = DEFAULTUNDEFINEDVALUE
      !
      if ( T > 0.d0 ) then
        C = a0 - a1*(273.16/T -1.0d0) + a2*(1.0d0 - T/273.16)
      else
        C = ILLEGALTEMP
      end if
    end function C

    function dC_dT ( T )
      ! As found in ref.
      real(r8), intent(in)   :: T
      real(r8)               :: dC_dT
      ! Local
      real(r8), parameter    :: a1 = 9.09718d0
      real(r8), parameter    :: a2 = 0.876793d0
      real, parameter        :: ILLEGALTEMP = DEFAULTUNDEFINEDVALUE
      !
      if ( T > 0.d0 ) then
        dC_dT = a1*(273.16/T**2) - a2/273.16
      else
        dC_dT = ILLEGALTEMP
      end if
    end function dC_dT

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
! Revision 2.8  2015/09/24 23:11:34  pwagner
! Avoid crashing when conversion would overflow otherwise
!
! Revision 2.7  2013/06/12 02:21:15  vsnyder
! Cruft removal
!
! Revision 2.6  2013/03/01 01:06:57  pwagner
! Get R8 from MLSKinds; DEFAULTUNDEFINEDVALUE from MLSCommon
!
! Revision 2.5  2009/08/24 20:08:55  pwagner
! May Fill H2O precision from RHI precision
!
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
