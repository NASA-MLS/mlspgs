! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE MLSL1Config  ! Level 1 Configuration
!=============================================================================

  USE MLSCommon, ONLY: TAI93_Range_T
  USE MLSL1Common, ONLY: MaxMIFs, BandSwitch, NumBands
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Info
  USE Init_tables_module, ONLY: First_Parm, Last_Parm
  USE Intrinsic, ONLY: parm_indices
  USE Output_m, ONLY: prunit, Output
  USE STRING_TABLE, ONLY: Get_string

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: L1Config_T, L1Config, Globals_T, Calib_T, Output_T, GHz_seq, &
       GHz_seq_use, THz_seq, THz_seq_use, MIFsGHz, MIFsTHz
  PUBLIC :: GetL1Config

  INTEGER, PARAMETER :: MIFsGHz = 125   ! Length (MIFs) of GHz module data
  INTEGER, PARAMETER :: MIFsTHz = 125   ! Length (MIFs) of THz module data

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  TYPE Globals_T
     CHARACTER(LEN=80) :: OutputVersionString
     CHARACTER(LEN=80) :: VersionComment
     LOGICAL :: ProduceL1BOA = .TRUE.
     LOGICAL :: SimOA = .TRUE.
  END TYPE Globals_T

  TYPE Calib_T
     INTEGER :: CalWindow
     INTEGER :: MIFsPerMAF
     REAL :: GHzSpaceTemp, GHzTargetTemp
     REAL :: THzSpaceTemp, THzTargetTemp
     REAL :: THzSpaceAngle, THzMaxBias
     REAL :: MIF_duration, MIF_DeadTime
     REAL :: MoonToSpaceAngle, MoonToLimbAngle
     LOGICAL :: UseDefaultGains = .FALSE.
     LOGICAL :: CalibDACS = .TRUE.
     CHARACTER(LEN=1) :: GHz_seq(0:MaxMIFs-1), THz_seq(0:MaxMIFs-1)
     CHARACTER(LEN=1) :: GHz_seq_use, THz_seq_use
  END TYPE Calib_T

  TYPE Output_T
     LOGICAL :: RemoveBaseline = .TRUE.             ! For GHz Baseline removal
     LOGICAL :: EnableChi2Err(NumBands) = .FALSE.   ! For RadErr calculation
  END TYPE Output_T

  TYPE L1Config_T
     TYPE (TAI93_Range_T) :: Input_TAI, Expanded_TAI
     TYPE (Globals_T) :: Globals
     TYPE (Calib_T) :: Calib
     TYPE (Output_T) :: Output
  END TYPE L1Config_T

  TYPE (L1Config_T), SAVE, TARGET :: L1Config

  CHARACTER(LEN=1), POINTER, DIMENSION(:) :: GHz_seq, THz_seq
  CHARACTER(LEN=1), POINTER :: GHz_seq_use, THz_seq_use

  LOGICAL :: GotParm(First_Parm:Last_Parm) = .FALSE.
  CHARACTER(LEN=32) :: ParmName

  CONTAINS

!=============================================================================
    SUBROUTINE GetL1Config
!=============================================================================

      USE Declaration_Table, ONLY: Allocate_Decl
      USE Init_tables_module, ONLY: Init_tables, z_globalsettings, &
           z_calibration, z_output
      USE Lexer_Core, ONLY: Init_Lexer
      USE MLSPCF1, ONLY: mlspcf_l1cf_start
      USE Parser, ONLY: Configuration
      USE SDPToolkit, ONLY: PGS_PC_GetReference, PGS_S_SUCCESS, &
           PGSd_IO_Gen_RSeqFrm, PGS_IO_Gen_openF, PGS_IO_Gen_closeF
      USE STRING_TABLE, ONLY: l1cf_unit => inunit
      USE Tree, ONLY: Allocate_Tree, Decoration, Nsons, Subtree
      USE Tree_checker, ONLY: Check_tree
      USE Units, ONLY: Init_units

      CHARACTER (LEN=132) :: physicalFilename
      INTEGER :: i, returnStatus, version
      INTEGER :: error           ! error return from tree checker
      INTEGER :: first_section   ! index of son of root of first n_cf_node
      INTEGER :: root            ! of the abstract syntax tree
      INTEGER :: son             ! of root

!! Open config file:

      version = 1
      returnStatus = PGS_PC_getReference (mlspcf_l1cf_start, version, &
           & physicalFilename)

      version = 1

      returnStatus = PGS_IO_Gen_openF (mlspcf_l1cf_start, PGSd_IO_Gen_RSeqFrm, &
           0, l1cf_unit, version)
      IF (returnstatus /= PGS_S_SUCCESS) THEN
         CALL MLSMessage (MLSMSG_Error, ModuleName, &
              & "Could not open L1 Config file: " // physicalFilename)
      END IF

      CALL MLSMessage (MLSMSG_Info, ModuleName, &
           & "Opened L1 Config file: " // physicalFilename)

!! Initialize the lexer, symbol table and tree checker's tables:

      CALL Init_Lexer (n_chars=10000, n_symbols=1000, hash_table_size=1003)
      CALL Allocate_Decl (ndecls=1000)
      CALL Allocate_Tree (n_tree=10000)
      CALL Init_tables
      CALL Init_units

!! Produce the abstract syntax tree

      CALL Configuration (root)

      IF (root <= 0) THEN

         CALL MLSMessage (MLSMSG_Error, ModuleName, &
              "Syntax error -- no abstract syntax tree")

      ENDIF

!! Close config file:

      returnStatus = PGS_IO_Gen_CloseF (l1cf_unit)

      IF (returnstatus /= PGS_S_SUCCESS) THEN
         CALL MLSMessage (MLSMSG_Error, ModuleName, &
              & "Could not close L1 Config file")
      END IF

      CALL MLSMessage (MLSMSG_Info, ModuleName, &
           & "Closed L1 Config file")

!! Check tree syntax

      prunit = -9                  ! to output to MLSMessage (and terminal)

      CALL Check_tree (root, error, first_section)

      IF (error /= 0 .OR. first_section <= 0) THEN

         CALL MLSMessage (MLSMSG_Error, ModuleName, "L1 CF Syntax check error")

      ENDIF

      DO i = first_section, Nsons (root)

         son = Subtree (i, root)

         SELECT CASE (Decoration (Subtree (1, son)))

         CASE (z_globalsettings)

            CALL Set_globalsettings (son)

         CASE (z_calibration)

            CALL Set_calibration (son)

         CASE (z_output)

            CALL Set_output (son)

         CASE DEFAULT

            CALL MLSMessage (MLSMSG_Error, ModuleName, 'Unknown section')

         END SELECT

      ENDDO


! Check for any missing required parameters:

      IF (ANY (.NOT. (GotParm))) THEN  ! All parameters are required!

         DO i = First_Parm, Last_Parm
            IF (.NOT. GotParm(i)) THEN
               CALL Get_String (Parm_indices(i), ParmName)
               CALL Output ('Missing L1CF parameter: '//TRIM(ParmName), &
                    advance='yes')
            ENDIF
         ENDDO

         CALL MLSMessage (MLSMSG_Error, ModuleName, 'Missing L1CF parameter(s)')

      ENDIF

    END SUBROUTINE GetL1Config

!=============================================================================
    SUBROUTINE Set_globalsettings (root)
!=============================================================================

      USE INIT_TABLES_MODULE, ONLY: p_output_version_string, &
           p_version_comment, p_produce_l1boa, p_simoa
      USE TREE, ONLY: Decoration, Nsons, Subtree, Sub_rosa
      USE MoreTree, ONLY: Get_Boolean

      INTEGER :: root

      INTEGER :: i, son

      DO i = 2, Nsons (root) - 1

         son = Subtree (i, root)

         GotParm(Decoration (Subtree (1,son))) = .TRUE.

         SELECT CASE (Decoration (Subtree (1,son)))

         CASE (p_output_version_string)

            CALL Get_string (Sub_rosa (Subtree(2,son)), &
                 L1Config%Globals%OutputVersionString, strip=.TRUE.)

         CASE (p_version_comment)

            CALL Get_string (Sub_rosa (Subtree(2,son)), &
                 L1Config%Globals%VersionComment, strip=.TRUE.)

         CASE (p_produce_l1boa)

            L1Config%Globals%ProduceL1BOA = Get_Boolean (son)

         CASE (p_simoa)

            L1Config%Globals%SimOA = Get_Boolean (son)

         END SELECT

      ENDDO

    END SUBROUTINE Set_globalsettings

!=============================================================================
   SUBROUTINE Set_output (root)
!=============================================================================

      USE EXPR_M, ONLY: Expr
      USE INIT_TABLES_MODULE, ONLY: p_removebaseline, s_chi2err, f_bandno
      USE TREE, ONLY: Decoration, Nsons, Subtree, Sub_rosa, Node_id
      USE TREE_TYPES
      USE MLSStrings, ONLY: lowercase
      USE MoreTree, ONLY: Get_Boolean

      INTEGER :: root, i, j, k, son, key, spec
      INTEGER :: expr_units(2)
      DOUBLE PRECISION :: expr_value(2)
      LOGICAL :: chi2entry

      chi2entry = .FALSE.   ! initialize to no entry yet

      DO i = 2, Nsons (root) - 1

         son = Subtree (i, root)

         SELECT CASE (node_id (son))

         CASE (n_equal)

            GotParm(Decoration (Subtree (1,son))) = .TRUE.

            SELECT CASE (Decoration (Subtree (1,son)))

            CASE (p_removebaseline)

               L1Config%Output%RemoveBaseline = Get_Boolean (son)

            END SELECT

         CASE (n_spec_args)

            key = son
            spec = decoration (subtree (1, decoration(subtree (1, key))))

            SELECT CASE (spec)

            CASE (s_chi2err)

               chi2entry = .TRUE.

               DO j = 2, nsons (key)

                  son = subtree (j, key)

                  SELECT CASE (decoration (subtree(1,son)))   ! field

                  CASE (f_bandno)

                     DO k = 2, nsons (son)
                        CALL Expr (subtree (k, son), expr_units, expr_value)
                        IF (MINVAL (INT (expr_value)) < 0 .OR. &
                             MAXVAL (INT (expr_value)) > 34) THEN
                           CALL MLSMessage (MLSMSG_Error, ModuleName, &
                            'Input band number out of range!')
                        ENDIF
                        L1Config%Output%EnableChi2Err(INT(expr_value(1)): &
                             INT(expr_value(2))) = .TRUE.
                     ENDDO

                  END SELECT

               ENDDO

            END SELECT

         END SELECT

      ENDDO

      IF (chi2entry) THEN
         CALL MLSMessage (MLSMSG_Info, ModuleName, &
              'Chi2 adjustments to Radiance precisions are ENABLED.')
      ELSE
         CALL MLSMessage (MLSMSG_Info, ModuleName, &
              'Chi2 adjustments to Radiance precisions are NOT ENABLED.')
      ENDIF

    END SUBROUTINE Set_output

!=============================================================================
    SUBROUTINE Set_calibration (root)
!=============================================================================

      USE EXPR_M, ONLY: Expr
      USE INIT_TABLES_MODULE, ONLY: p_calwindow, s_spaceMIFs, s_targetMIFs, &
           s_limbMIFs, s_discardMIFs, f_mifs, f_use, l_match, l_override, &
           f_module, f_secondary, p_usedefaultgains, p_GHzSpaceTemp, &
           p_GHzTargetTemp, p_THzSpaceTemp, p_THzTargetTemp, p_mif_duration, &
           p_mif_dead_time, p_mifspermaf, p_calibDACS, p_THzMaxBias, s_switch, &
           p_thzspaceangle, f_s, f_bandno, p_MoonToSpaceAngle, p_MoonToLimbAngle
      USE INTRINSIC, ONLY: l_ghz, l_thz, phyq_mafs, phyq_temperature, &
           phyq_mifs, phyq_time, phyq_angle
      USE TREE, ONLY: Decoration, Nsons, Subtree, Sub_rosa, Node_id
      USE TREE_TYPES
      USE MoreTree, ONLY: Get_Boolean

      INTEGER :: root

      CHARACTER(LEN=1), POINTER, DIMENSION(:) :: scan_seq
      CHARACTER(LEN=1), POINTER :: scan_use
      CHARACTER(LEN=80) :: identifier
      INTEGER :: i, j, k, son, key, spec, swno, bandno
      INTEGER :: expr_units(2)
      DOUBLE PRECISION :: expr_value(2)
      LOGICAL :: GHz_mod, sec_tgt

      TYPE Scan_T
         CHARACTER(LEN=1) :: Use    ! M or O
         CHARACTER(LEN=1) :: Type   ! L, S, T, or D
         INTEGER :: MIF(0:MaxMIFs-1)
      END TYPE Scan_T

      TYPE (Scan_T) :: scan

      CHARACTER(LEN=1), PARAMETER :: ignore = "I", unknown = "U", overlap = "O"

! Initialize desired scan sequences

      L1Config%Calib%GHz_seq = unknown
      L1Config%Calib%GHz_seq_use = ignore

      L1Config%Calib%THz_seq = unknown
      L1Config%Calib%THz_seq_use = ignore

! Initialize Target Temp to not input (and not required!)

      L1Config%Calib%GHzTargetTemp = -1.0
      L1Config%Calib%THzTargetTemp = -1.0

      GotParm(p_GHzTargetTemp) = .TRUE.
      GotParm(p_THzTargetTemp) = .TRUE.

      DO i = 2, Nsons (root) - 1

         son = Subtree (i, root)

         SELECT CASE (node_id (son))

         CASE (n_equal)

            GotParm(Decoration (Subtree (1,son))) = .TRUE.

            SELECT CASE (decoration (subtree (1,son)))

            CASE (p_calwindow)

               CALL Expr (subtree (2, son), expr_units, expr_value)
               L1Config%Calib%CalWindow = expr_value(1)
               IF (expr_units(1) /= phyq_mafs) THEN
                  CALL Get_string (Sub_rosa (Subtree(1,son)), identifier)
                  CALL MLSMessage (MLSMSG_Error, ModuleName, &
                       TRIM (identifier)//' is not input as MAFs')
               ENDIF

            CASE (p_mifspermaf)

               CALL Expr (subtree (2, son), expr_units, expr_value)
               L1Config%Calib%MIFsPerMAF = expr_value(1)
               IF (expr_units(1) /= phyq_mifs) THEN
                  CALL Get_string (Sub_rosa (Subtree(1,son)), identifier)
                  CALL MLSMessage (MLSMSG_Error, ModuleName, &
                       TRIM (identifier)//' is not input as MIFs')
               ENDIF

            CASE (p_GHzSpaceTemp)

               CALL Expr (subtree (2, son), expr_units, expr_value)
               L1Config%Calib%GHzSpaceTemp = expr_value(1)
               IF (expr_units(1) /= phyq_temperature) THEN
                  CALL Get_string (Sub_rosa (Subtree(1,son)), identifier)
                  CALL MLSMessage (MLSMSG_Error, ModuleName, &
                       TRIM (identifier)//' is not input as K')
               ENDIF

            CASE (p_GHzTargetTemp)

               CALL Expr (subtree (2, son), expr_units, expr_value)
               L1Config%Calib%GHzTargetTemp = expr_value(1)
               IF (expr_units(1) /= phyq_temperature) THEN
                  CALL Get_string (Sub_rosa (Subtree(1,son)), identifier)
                  CALL MLSMessage (MLSMSG_Error, ModuleName, &
                       TRIM (identifier)//' is not input as K')
               ENDIF

            CASE (p_THzSpaceTemp)

               CALL Expr (subtree (2, son), expr_units, expr_value)
               L1Config%Calib%THzSpaceTemp = expr_value(1)
               IF (expr_units(1) /= phyq_temperature) THEN
                  CALL Get_string (Sub_rosa (Subtree(1,son)), identifier)
                  CALL MLSMessage (MLSMSG_Error, ModuleName, &
                       TRIM (identifier)//' is not input as K')
               ENDIF

            CASE (p_THzTargetTemp)

               CALL Expr (subtree (2, son), expr_units, expr_value)
               L1Config%Calib%THzTargetTemp = expr_value(1)
               IF (expr_units(1) /= phyq_temperature) THEN
                  CALL Get_string (Sub_rosa (Subtree(1,son)), identifier)
                  CALL MLSMessage (MLSMSG_Error, ModuleName, &
                       TRIM (identifier)//' is not input as K')
               ENDIF

            CASE (p_THzSpaceAngle)

               CALL Expr (subtree (2, son), expr_units, expr_value)
               L1Config%Calib%THzSpaceAngle = expr_value(1)
               IF (expr_units(1) /= phyq_angle) THEN
                  CALL Get_string (Sub_rosa (Subtree(1,son)), identifier)
                  CALL MLSMessage (MLSMSG_Error, ModuleName, &
                       TRIM (identifier)//' is not input as deg[rees]')
               ENDIF

            CASE (p_THzMaxBias)

               CALL Expr (subtree (2, son), expr_units, expr_value)
               L1Config%Calib%THzMaxBias = expr_value(1)
!!$               IF (expr_units(1) /= phyq_temperature) THEN ! Add Volts later?
!!$                  CALL Get_string (Sub_rosa (Subtree(1,son)), identifier)
!!$                  CALL MLSMessage (MLSMSG_Error, ModuleName, &
!!$                       TRIM (identifier)//' is not input as K')
!!$               ENDIF

            CASE (p_MoonToSpaceAngle)

               CALL Expr (subtree (2, son), expr_units, expr_value)
               L1Config%Calib%MoonToSpaceAngle = expr_value(1)
               IF (expr_units(1) /= phyq_angle) THEN
                  CALL Get_string (Sub_rosa (Subtree(1,son)), identifier)
                  CALL MLSMessage (MLSMSG_Error, ModuleName, &
                       TRIM (identifier)//' is not input as deg[rees]')
               ENDIF

            CASE (p_MoonTolimbAngle)

               CALL Expr (subtree (2, son), expr_units, expr_value)
               L1Config%Calib%MoonToLimbAngle = expr_value(1)
               IF (expr_units(1) /= phyq_angle) THEN
                  CALL Get_string (Sub_rosa (Subtree(1,son)), identifier)
                  CALL MLSMessage (MLSMSG_Error, ModuleName, &
                       TRIM (identifier)//' is not input as deg[rees]')
               ENDIF

            CASE (p_mif_duration)

               CALL Expr (subtree (2, son), expr_units, expr_value)
               L1Config%Calib%MIF_duration = expr_value(1)
               IF (expr_units(1) /= phyq_time) THEN
                  CALL Get_string (Sub_rosa (Subtree(1,son)), identifier)
                  CALL MLSMessage (MLSMSG_Error, ModuleName, &
                       TRIM (identifier)//' is not input as seconds')
               ENDIF

            CASE (p_mif_dead_time)

               CALL Expr (subtree (2, son), expr_units, expr_value)
               L1Config%Calib%MIF_DeadTime = expr_value(1)
               IF (expr_units(1) /= phyq_time) THEN
                  CALL Get_string (Sub_rosa (Subtree(1,son)), identifier)
                  CALL MLSMessage (MLSMSG_Error, ModuleName, &
                       TRIM (identifier)//' is not input as seconds')
               ENDIF

            CASE (p_usedefaultgains)

               L1Config%Calib%UseDefaultGains = Get_Boolean (son)

            CASE (p_calibDACS)

               L1Config%Calib%calibDACS = Get_Boolean (son)

            END SELECT

         CASE (n_spec_args)

            key = son
            spec = decoration (subtree (1, decoration(subtree (1, key))))

            SELECT CASE (spec)

            CASE (s_spaceMIFs, s_targetMIFs, s_limbMIFs, s_discardMIFs)

               SELECT CASE (spec)

               CASE (s_spaceMIFs)
                  scan%Type = "S"
               CASE (s_targetMIFs)
                  scan%Type = "T"
                  sec_tgt = .FALSE.
               CASE (s_limbMIFs)
                  scan%Type = "L"
               CASE (s_discardMIFs)
                  scan%Type = "D"
               END SELECT

               scan%MIF = 0

               DO j = 2, nsons (key)

                  son = subtree (j, key)

                  SELECT CASE (decoration (subtree(1,son)))   ! field

                  CASE (f_mifs)
                     DO k = 2, nsons (son)
                        CALL Expr (subtree (k, son), expr_units, expr_value)
                        IF (MINVAL (INT (expr_value)) < 0 .OR. &
                             MAXVAL (INT (expr_value)) > 149) THEN
                           CALL MLSMessage (MLSMSG_Error, ModuleName, &
                            'Input MIF number out of range for scan type "'// &
                            scan%Type//'"!')
                        ENDIF
                      scan%MIF(INT(expr_value(1)):INT(expr_value(2))) = 1
                     ENDDO
                  CASE (f_use)
                     SELECT CASE (decoration (subtree (2, son)))

                     CASE (l_match)
                        scan%Use = "M"
                     CASE (l_override)
                        scan%Use = "O"
                     END SELECT

                  CASE (f_secondary)

                     sec_tgt = Get_Boolean (son)

                  CASE (f_module)

                     SELECT CASE (decoration (subtree (2, son)))

                     CASE (l_ghz)
                        scan_use => L1Config%Calib%GHz_seq_use
                        scan_seq => L1Config%Calib%GHz_seq
                        GHz_mod = .TRUE.
                     CASE (l_thz)
                        scan_use => L1Config%Calib%THz_seq_use
                        scan_seq => L1Config%Calib%THz_seq
                        GHz_mod = .FALSE.
                     END SELECT

                  END SELECT

               ENDDO

               IF (scan_use == ignore) THEN
                  scan_use = scan%Use
               ELSE IF (scan_use /= scan%Use) THEN
                  CALL MLSMessage (MLSMSG_Error, ModuleName, &
                       'USE fields do not match in L1CF Calib section!')
               ENDIF

!! GHz module Target type:

               IF (scan%type == "T" .AND. GHz_mod) THEN

!! Reverse target type if secondary target

                  IF (sec_tgt) scan%type = "t"

               ENDIF

               WHERE (scan_seq == unknown .AND. scan%MIF == 1)
                  scan_seq = scan%type
               ELSEWHERE
                  WHERE (scan_seq /= unknown .AND. scan%MIF == 1)
                     scan_seq = overlap
                  END WHERE
               END WHERE
               IF (ANY (scan_seq == overlap)) THEN
                  CALL MLSMessage (MLSMSG_Error, ModuleName, &
                       'Scan sequence overlap in L1CF Calib section!')
               ENDIF

            CASE (s_switch)

               DO j = 2, nsons (key)

                  son = subtree (j, key)

                  SELECT CASE (decoration (subtree(1,son)))   ! field

                  CASE (f_s)

                     CALL Expr (subtree (2, son), expr_units, expr_value)
                     swno = expr_value(1)
                     IF (swno < 1 .OR. swno > 5) THEN
                        CALL MLSMessage (MLSMSG_Error, ModuleName, &
                             'Switch number out of range (1-5)!')
                     ENDIF
                     
                  CASE (f_bandno)
                     CALL Expr (subtree (2, son), expr_units, expr_value)
                     bandno = expr_value(1)

                  END SELECT

               ENDDO

               BandSwitch(swno) = bandno

            CASE DEFAULT

               PRINT *, 'unknown spec!'

         END SELECT

         CASE DEFAULT

            PRINT *, 'cal default son', son

         END SELECT

      ENDDO

! Save pointers for sort/qualify:

      GHz_seq_use => L1Config%Calib%GHz_seq_use
      GHz_seq => L1Config%Calib%GHz_seq
      THz_seq_use => L1Config%Calib%THz_seq_use
      THz_seq => L1Config%Calib%THz_seq

      IF (ANY (L1Config%Calib%GHz_seq == unknown)) THEN
         CALL MLSMessage (MLSMSG_Error, ModuleName, &
              'Undefined MIF(s) for GHz module in L1CF Calib section!')
      ENDIF

      IF (ANY (L1Config%Calib%THz_seq == unknown)) THEN
         CALL MLSMessage (MLSMSG_Error, ModuleName, &
              'Undefined MIF(s) for THz module in L1CF Calib section!')
      ENDIF

!! Check switches

      DO swno = 1, 5
         IF (BandSwitch(swno) == 0) THEN
            WRITE (identifier, "(I1)") swno
            CALL MLSMessage (MLSMSG_Error, ModuleName, &
                 'Undefined switch no '//TRIM(identifier)//&
                 ' in L1CF Calib section!')
         ENDIF
      ENDDO

    END SUBROUTINE Set_calibration

END MODULE MLSL1Config

! $Log$
! Revision 2.14  2004/08/12 13:51:50  perun
! Version 1.44 commit
!
! Revision 2.13  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
! Revision 2.12  2004/01/09 17:46:22  perun
! Version 1.4 commit
!
! Revision 2.11  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.10  2002/11/21 16:57:46  perun
! Moved default HDFversion number to declaration
!
! Revision 2.9  2002/11/20 18:30:50  perun
! Restated HDFVersionString error message
!
! Revision 2.8  2002/11/19 20:35:16  perun
! Convert HDFVersionstring to HDFversion number
!
! Revision 2.7  2002/11/14 16:51:13  perun
! Split space & target temps between GHz & THz
!
! Revision 2.6  2002/11/07 21:35:03  jdone
! Added HDF4/HDF5 switch.
!
! Revision 2.5  2002/04/04 19:30:49  perun
! Check for correct input MIFs
!
! Revision 2.4  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.3  2001/04/27 14:00:32  perun
! For the latest parser version
!
! Revision 2.2  2001/03/22 16:45:06  perun
! Changed call to Get_string to strip "'s from globals
!
! Revision 2.1  2001/02/23 20:50:54  perun
! Version 0.5 commit
!
