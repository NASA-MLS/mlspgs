module D_Q_LOG_M
  implicit NONE
  private
  public :: Q_LOG   ! generic
  public :: D_Q_LOG ! specific
  interface Q_LOG; module procedure D_Q_LOG; end interface
  integer, private, parameter :: RK = kind(0.0d0)
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
!p|
!p|
!p|
!p|  ############################### PROLOGUE ###############################
!p|  #                                                                      #
!p|  #                               Q _ L O G                              #
!p|  #                                                                      #
!p|  ########################################################################
!p|
!p|
!p|  PURPOSE: LOG OF PARTITION FUNCTION NORMALIZED TO Q(T=300K)
!p|
!p|      Returns the base 10 logarithm of the ratio of partition function
!p|              at temperature T to that at temperature 300K
!p|
!p|      Linear interpolation in log(T) and log(Q) is used.  JPL catalog
!p|              accuracy is maintained over temperature range 150-300K.
!p|
!p|
!p|  ABSTRACT:
!p|
!p|
!p|  CORRESPONDING REQUIREMENTS:
!p|
!p|
!p|  INVOCATION METHOD:  r = q_log(q0,t)
!p|
!p|
!p|  ARGUMENTS:
!p|
!p|   Name            Type    I/O    Purpose
!p|  -------------   ------   ---    ---------------------------------
!p|   q0               R*8     I     3-elements array of base 10 logarithms of
!p|                                  partition function at T=300, 225, 150K
!p|                                  from JPL catalog. [Poynter and Pickett,
!p|                                  Appl. Opt. 24, 2235, July 1985]
!p|   t                R*8     I     temperature (K)
!p|   q_log            R*8     O     output
!p|
!p|
!p|  LOGICAL/FILE REFERENCES: None
!p|
!p|
!p|  INCLUDED PARAMETERS & COMMON BLOCKS: None
!p|
!p|
!p|  PARENT: bin_pqm_intrp
!p|
!p|
!p|  EXTERNAL REFERENCES: None
!p|
!p|
!p|  NOTES:
!p|    Written by J.Waters.   JPL  2 June 1986
!p|
!p|
!p|  ERROR HANDLING: None
!p|
!p|
!p|  CONFIGURATION HISTORY:
!p|
!p|   Author        Date      Comments
!p|  ----------   --------   ------------------------------------------------
!p| J.Waters      02/06/86   Initial design
!p| Z. Shippony   09/30/97   Adaptation to SID/L2PC
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
  Function D_Q_Log(Q0,T) result ( RESULT )
    real(rk), intent(in) :: Q0(3), T
    real(rk) :: result
!
    real(rk) :: SLOPE, TLOG
!
!   tlog0 = Log10(T) for T = 300.0, 225.0, 150.0
!
    real(rk), parameter :: tlog0(3) = &
   &  (/ 2.47712125471966_rk, 2.35218251811136_rk, 2.17609125905568_rk /)
!
    tlog=log10(T)
!
    if(tlog < tlog0(2)) then
!
      slope=(q0(2)-q0(3))/(tlog0(2)-tlog0(3))
      result=q0(2)-q0(1)+slope*(tlog-tlog0(2))
!
    else
!
      slope=(q0(1)-q0(2))/(tlog0(1)-tlog0(2))
      result=slope*(tlog-tlog0(1))
!
    endif
!
    Return
  End Function D_Q_Log
end module D_Q_LOG_M
! $Log$
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
!
