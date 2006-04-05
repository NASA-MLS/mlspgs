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
MODULE BrightObjects_m
!=============================================================================

  USE Init_MLSSignals_m
  USE INTRINSIC
  USE MLSL1Common, ONLY: MaxMIFs

  IMPLICIT NONE

  PUBLIC

  PRIVATE :: add_ident, make_tree

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

! Enumeration types:

  INTEGER, PARAMETER :: t_module = last_signal_type + 1
  INTEGER, PARAMETER :: t_name = t_module + 1
  INTEGER, PARAMETER :: last_BrightObject_type = t_name

! Field indices:

  INTEGER, PARAMETER :: f_angle = last_Signal_Field + 1
  INTEGER, PARAMETER :: f_name = f_angle + 1
  INTEGER, PARAMETER :: f_negate = f_name + 1
  INTEGER, PARAMETER :: last_BrightObject_Field = f_negate

! Enumeration literals:

  INTEGER, PARAMETER :: l_mercury =  last_signal_lit + 1
  INTEGER, PARAMETER :: l_venus =  l_mercury + 1
  INTEGER, PARAMETER :: l_earth =  l_venus + 1
  INTEGER, PARAMETER :: l_mars =  l_earth + 1
  INTEGER, PARAMETER :: l_jupiter =  l_mars + 1
  INTEGER, PARAMETER :: l_saturn =  l_jupiter + 1
  INTEGER, PARAMETER :: l_uranus =  l_saturn + 1
  INTEGER, PARAMETER :: l_neptune =  l_uranus + 1
  INTEGER, PARAMETER :: l_pluto =  l_neptune + 1
  INTEGER, PARAMETER :: l_moon =  l_pluto + 1
  INTEGER, PARAMETER :: l_sun =  l_moon + 1
  INTEGER, PARAMETER :: l_sbaryctr =  l_sun + 1
  INTEGER, PARAMETER :: l_ebaryctr =  l_sbaryctr + 1
  INTEGER, PARAMETER :: l_gcenter =  l_ebaryctr + 1
  INTEGER, PARAMETER :: last_BrightObject_lit =  l_gcenter

! Specification indices:

  INTEGER, PARAMETER :: s_BrightObject = last_Signal_Spec + 1
  INTEGER, PARAMETER :: s_Name = s_BrightObject + 1
  INTEGER, PARAMETER :: last_BrightObject_Spec = s_Name

! Using the same definitions as the Toolkit values (except for the Galactic
!  Center value)

  INTEGER, PARAMETER :: BO_MaxNum = 14  ! planets plus moon plus GC

  CHARACTER (len=14), PARAMETER :: BO_name(BO_MaxNum)= (/&
       "MERCURY       ", "VENUS         ", "EARTH         ", &
       "MARS          ", "JUPITER       ", "SATURN        ", &
       "URANUS        ", "NEPTUNE       ", "PLUTO         ", &
       "MOON          ", "SUN           ", "SOLAR_BARYCNTR", &
       "EARTH_BARYCNTR", "GALACTICCENTER" /)
  INTEGER :: BO_NumGHz = 0, BO_Index_GHz(BO_MaxNum) = 0  ! GHz BO's to check
  INTEGER :: BO_NumTHz = 0, BO_Index_THz(BO_MaxNum) = 0  ! THz BO's to check
  LOGICAL :: BO_Negate_GHz(BO_MaxNum) = .FALSE.  ! GHz BO's to negate
  LOGICAL :: BO_Negate_THz(BO_MaxNum) = .FALSE.  ! THz BO's to negate

  REAL :: BO_Angle_GHz(BO_MaxNum) = -999.9, &               ! BO angles for FOV
       BO_Angle_THz(BO_MaxNum) = -999.9
  INTEGER, PARAMETER :: GC_Def = l_gcenter - l_mercury + 1  ! GC def for tests

! Matched BO(s) structure:

  TYPE BO_Match_T
     INTEGER :: num
     CHARACTER (len=LEN(BO_name)) :: Name(BO_MaxNum)
     REAL :: PrecScale(BO_MaxNum)
     LOGICAL :: InFOV(0:(MaxMIFs-1),BO_MaxNum)
  END TYPE BO_Match_T
  TYPE (BO_Match_T) :: BO_Match

! BO_stat for THz processing (read from L1BOA file):

  INTEGER, DIMENSION(:,:), ALLOCATABLE :: THz_BO_stat

CONTAINS

!=============================================================================
  SUBROUTINE Init_BrightObjects
!=============================================================================

    USE TREE
    USE TREE_TYPES, ONLY: N_DT_DEF, N_FIELD_TYPE, N_SPEC_DEF
 
  ! Put nonintrinsic predefined identifiers into the symbol table.

    ! Put enumeration literals into the symbol table:

    lit_indices (l_mercury) =                add_ident (TRIM(BO_name(1)))
    lit_indices (l_venus) =                  add_ident (TRIM(BO_name(2)))
    lit_indices (l_earth) =                  add_ident (TRIM(BO_name(3)))
    lit_indices (l_mars) =                   add_ident (TRIM(BO_name(4)))
    lit_indices (l_jupiter) =                add_ident (TRIM(BO_name(5)))
    lit_indices (l_saturn) =                 add_ident (TRIM(BO_name(6)))
    lit_indices (l_uranus) =                 add_ident (TRIM(BO_name(7)))
    lit_indices (l_neptune) =                add_ident (TRIM(BO_name(8)))
    lit_indices (l_pluto) =                  add_ident (TRIM(BO_name(9)))
    lit_indices (l_moon) =                   add_ident (TRIM(BO_name(10)))
    lit_indices (l_sun) =                    add_ident (TRIM(BO_name(11)))
    lit_indices (l_sbaryctr) =               add_ident (TRIM(BO_name(12)))
    lit_indices (l_ebaryctr) =               add_ident (TRIM(BO_name(13)))
    lit_indices (l_gcenter) =                add_ident (TRIM(BO_name(14)))

    ! Put enumeration type names into the symbol table

    data_type_indices (t_module) =           add_ident ( 'module' )
    data_type_indices (t_name) =             add_ident ( 'name' )

    ! Put field names into the symbol table

    field_indices (f_angle) =                add_ident ( 'angle' )
    field_indices (f_name) =                 add_ident ( 'name' )
    field_indices (f_negate) =               add_ident ( 'negate' )
  
    ! Put spec names into the symbol table

    spec_indices (s_BrightObject) =          add_ident ( 'BrightObject' )
    spec_indices (s_name) =                  add_ident ( 'name' )

    CALL make_tree ( (/ &
      begin, t+t_module, l+l_ghz, l+l_thz, n+n_dt_def /))

    CALL make_tree ( (/ &
      begin, s+s_name, &
             begin, f+f_name, t+t_boolean, n+n_field_type, &
             np+n_spec_def /))

    CALL make_tree ( (/ &
      begin, t+t_name, l+l_mercury, l+l_venus, l+l_earth, l+l_mars, &
      l+l_jupiter, l+l_saturn, l+l_uranus, l+l_neptune, l+l_pluto, &
      l+l_moon, l+l_sun, l+l_sbaryctr, l+l_ebaryctr, l+l_gcenter, n+n_dt_def /))

    CALL make_tree ( (/ &
      begin, s+s_BrightObject, &
             begin, f+f_angle, t+t_numeric, nr+n_field_type, &
             begin, f+f_negate, t+t_boolean, n+n_field_type, &
             begin, f+f_module, s+s_module, nr+n_field_type, &
             begin, f+f_name, s+s_name, nr+n_field_type, &
             ndp+n_spec_def /) )

  END SUBROUTINE Init_BrightObjects
    
  ! --------------------------------------------------  MAKE_TREE  -----
  INCLUDE "make_tree.f9h"

!=============================================================================
  SUBROUTINE Test_BO_stat (BO_stat)
!=============================================================================

    INTEGER, INTENT(in) :: BO_stat(:)

    LOGICAL :: BOinFOV(SIZE(BO_stat),BO_MaxNum)
    LOGICAL :: BO_Flag(0:15)         ! Flag for each B.O.
    REAL :: BO_scale(0:15)           ! Scale for precisions (+/- 1.0)

    INTEGER :: i, bitno, matchno, last_MIF

    BO_Match%Num = 0                 ! No matches yet
    IF (ALL (BO_stat <= 1)) RETURN   ! Nothing in limb FOV

    PRINT *, 'BO in FOV...'

    last_MIF = SIZE (BO_stat) - 1    ! BO_stat starts at 1
    BO_Flag = .FALSE.        ! Nothing yet
    BOinFOV = .FALSE.        ! Nothing yet
    BO_scale = 1.0
    DO i = 1, SIZE (BO_stat)
       DO bitno = 1, 15      ! test for BO in FOV (skip bit 0 since it's space)
          IF (BTEST (BO_stat(i), bitno)) THEN
             BOinFOV(i,bitno) = .TRUE.
             BO_Flag(bitno) = .TRUE.
             IF (BTEST (BO_stat(i), (bitno+16))) THEN   ! test for negating
                BO_scale(bitno) = -1.0                  ! negative scale factor
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    BO_Match%Num = COUNT (BO_Flag)
    matchno = 1
    DO i = 1, (SIZE (BO_Flag) - 1)
       IF (BO_Flag(i)) THEN
          BO_Match%Name(matchno) = BO_Name(i)
          BO_Match%PrecScale(matchno) = BO_scale(i)
          BO_Match%InFov(0:last_MIF,matchno)= BOinFOV(:,i)
          matchno = matchno + 1
       ENDIF
    ENDDO

  END SUBROUTINE Test_BO_stat

  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here

END MODULE BrightObjects_m

! $Log$
! Revision 2.3  2006/04/05 18:09:23  perun
! Remove unused variables
!
! Revision 2.2  2006/03/24 15:06:53  perun
! Corrected bit number index for BOinFOV
!
! Revision 2.1  2005/12/06 19:32:36  perun
! Initial release for bright object CF fields and for testing
!
