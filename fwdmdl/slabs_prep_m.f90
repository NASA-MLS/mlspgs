module SLABS_PREP_M
  use D_Q_LOG_M, only: Q_LOG
  use MLSCommon, only: R8
  implicit NONE
  private
  public :: SLABS_PREP

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

!p|
!p|
!p|  ############################### PROLOGUE ###############################
!p|  #                                                                      #
!p|  #                          S L A B S _ P R E P                         #
!p|  #                                                                      #
!p|  ########################################################################
!p|
!p|
!p|  PURPOSE: This function will compute all the parameters needed before
!p|           calling the sigle-line-absorption computation routine (slabs)
!p|
!p|
!p|  ABSTRACT:
!p|
!p|
!p|  CORRESPONDING REQUIREMENTS:
!p|
!p|
!p|  INVOCATION METHOD: call slabs_prep(t,m,v0,el,w,ps,p,n,i,q,delta,gamma,
!p|                                     n1,n2,v0s,x1,slabs1,y,yi)
!p|
!p|
!p|  ARGUMENTS:
!p|
!p|   Name                   Type    I/O    Purpose
!p|  --------------------   ------   ---   ----------------------------------
!p|   t                       R*8     I     Temperature K
!p|   m                       R*8     I     Molecular mass amu
!p|   v0                      R*8     I     Line center frequency MHz
!p|                                           (MHz)
!p|   el                      R*8     I     Lower state energy cm-1
!p|   w                       R*8     I     Collision broadening parameter
!p|                                           (MHz/mbar at 300 K)
!p|   ps                      R*8     I     Pressure shift parameter
!p|                                           (MHz/mbar)
!p|   p                       R*8     I     Pressure mbar
!p|   n                       R*8     I     Temperature power dependence of w
!p|   i                       R*8     I     Integrated spectral intensity
!p|                                           log(nm**2 MHz) at 300 K
!p|   q(3)                    R*8     I     Logarithm of the partition function
!p|                                           at 300 , 225 , and 150 K
!p|   delta                   R*8     I     Line mixing coefficient
!p|   gamma                   R*8     I     Line mixing coefficient
!p|   n1                      R*8     I     Temperature dependency of delta
!p|   n2                      R*8     I     Temperature dependency of gamma
!p|   v0s                     R*8     O     Pressure shifted line position
!p|   x1                      R*8     O     Sqrt(Ln(2))/Doppler half width
!p|                                           (MHz)
!p|   slabs1                  R*8     O     Frequency independent piece of
!p|                                           slabs
!p|   y                       R*8     O     [collision width/doppler width] *
!p|                                           sqrt(Ln(2))
!p|   yi                      R*8     O     Interference contribution
!p|
!p|
!p|  LOGICAL/FILE REFERENCES: None
!p|
!p|
!p|  INCLUDED PARAMETERS & COMMON BLOCKS: None
!p|
!p|
!p|  PARENT: z_slabs
!p|
!p|
!p|  EXTERNAL REFERENCES: None
!p|
!p|
!p|  NOTES:  The above outputs along with frequency offset are used with
!p|          routine SLABSWINT to compute a Single Line ABSorption in
!p|          1/Km units, with unit mixing ratio.
!p|
!p|
!p|  ERROR HANDLING: None
!p|
!p|
!p|  CONFIGURATION HISTORY:
!p|
!p|   Author        Date      Comments
!p|  ----------   --------   ------------------------------------------------
!p| Z. Shippony   09/30/97   Initial design
!p| T. Lungu      09/30/97   Initial VMS release
!p|
!p|
!p|  ________________________________________________________________________
!p|  alterations:
!p|  Z. Shippony  09/30/97   Added Prologue
!p|
!p|
!d|
!d|    _____________________________________________________________________
!d|   /                                                                     \
!d|  <                         -    D E S I G N    -                         >
!d|   \_____________________________________________________________________/
!d|
!d|
!d|        \\\\\\\\  SUBROUTINE  LINES_O_CODE  \\\\\\\\
!d|
!d|
!c|
!c| ==========================================================================
!c| CODE CODE CODE CODE CODE CODE CODE CODE CODE CODE CODE CODE CODE CODE CODE
!c| ==========================================================================
!
!---------------------------------------------------------------------
!
  Subroutine SLABS_PREP ( t, m, v0, el, w, ps, p, n, i, q, delta, gamma, &
 &                        n1, n2, v0s, x1, slabs1, y, yi )
    real(r8), intent(in) :: T          ! Temperature K
    real(r8), intent(in) :: M          ! Molecular mass amu
    real(r8), intent(in) :: V0         ! Line center frequency MHz
                                       !   (MHz)
    real(r8), intent(in) :: EL         ! Lower state energy cm-1
    real(r8), intent(in) :: W          ! Collision broadening parameter
                                       !   (MHz/mbar at 300 K)
    real(r8), intent(in) :: PS         ! Pressure shift parameter
                                       !   (MHz/mbar)
    real(r8), intent(in) :: P          ! Pressure mbar
    real(r8), intent(in) :: N          ! Temperature power dependence of w
    real(r8), intent(in) :: I          ! Integrated spectral intensity
                                       !   log(nm**2 MHz) at 300 K
    real(r8), intent(in) :: Q(3)       ! Logarithm of the partition function
                                       !   at 300 , 225 , and 150 K
    real(r8), intent(in) :: DELTA      ! Line mixing coefficient
    real(r8), intent(in) :: GAMMA      ! Line mixing coefficient
    real(r8), intent(in) :: N1         ! Temperature dependency of delta
    real(r8), intent(in) :: N2         ! Temperature dependency of gamma
    real(r8), intent(out) :: V0S       ! Pressure shifted line position
    real(r8), intent(out) :: X1        ! Sqrt(Ln(2))/Doppler half width
                                       !   (MHz)
    real(r8), intent(out) :: SLABS1    ! Frequency independent piece of
                                       !   slabs
    real(r8), intent(out) :: Y         ! [collision width/doppler width] *
                                       !   sqrt(Ln(2))
    real(r8), intent(out) :: YI        ! Interference contribution
!
! Internal constants:
!
    ! i2abs           ! converts intensity into absorption
    ! dc              ! sqrt(amu/K) used to calculate doppler
                      ! width
    ! boltzcm         ! boltzmann constant cm-1/K
    ! boltzmhz        ! boltzmann constant MHz/K
    ! sqrtln2         ! sqrt(ln(2))
    ! oned300         ! 1.0 / 300.0
    ! loge            ! log10 e
!
    real(r8), parameter :: i2abs = 3.402136078e9_r8
    real(r8), parameter :: dc = 3.58117369e-7_r8
    real(r8), parameter :: boltzcm = 0.6950387e0_r8
    real(r8), parameter :: boltzmhz = 20836.74e0_r8
    real(r8), parameter :: sqrtln2 = 8.32554611157698e-1_r8
    real(r8), parameter :: loge = 4.34294481903251828e-1_r8
    real(r8), parameter :: oned300 = 1.0_r8/300.0_r8
!
! Internal data:
!
    real(r8) :: DFD, IP, BETAE, BETAV, T3T, ONEDT, NS
!
! The action begins here
!
    onedt = 1.0_r8 / t
    t3t = 300.0_r8 * onedt
    ns = 0.25_r8 + 1.5_r8 * n
    v0s = v0 + ps * p * (t3t**ns)
    betae = el / boltzcm
    betav = v0s / boltzmhz
    dfd = v0s * sqrt(t/m) * dc
    x1 = sqrtln2 / dfd
    y = x1 * w * p * (t3t**n)
    ip = i - Q_Log(q,t) + loge *  betae * (oned300 - onedt)
    slabs1 = i2abs * p * (10.0_r8**ip) * (1.0_r8 - exp(-betav*onedt)) &
   &         / (dfd * t *(1.0_r8 - exp(-betav*oned300)))
    yi = p * (delta*(t3t**n1) + gamma*(t3t**n2))
!
    Return
  End Subroutine SLABS_PREP
end module SLABS_PREP_M

! $Log$

