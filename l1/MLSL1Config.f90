! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE MLSL1Config  ! Level 1 Configuration
!=============================================================================

  USE MLSCommon, ONLY: TAI93_Range_T
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Info, &
       MLSMSG_Warning

  IMPLICIT NONE

  INTEGER :: MIFsGHz   ! Length (MIFs) of GHz module data
  INTEGER :: MIFsTHz   ! Length (MIFs) of THz module data

  PRIVATE :: Id, ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  TYPE Globals_T
     CHARACTER(LEN=80) :: OutputVersionString
     CHARACTER(LEN=80) :: VersionComment
  END TYPE Globals_T

  TYPE L1Config_T
     TYPE (TAI93_Range_T) :: Input_TAI 
     TYPE (TAI93_Range_T) :: Expanded_TAI
     TYPE (Globals_T) :: Globals
  END TYPE L1Config_T

  TYPE (L1Config_T) :: L1Config

  CONTAINS

    SUBROUTINE GetL1Config

      USE Declaration_Table, ONLY: Allocate_Decl
      USE Init_tables_module, ONLY: Init_tables, lit_indices, &
           z_globalsettings, z_calibration
      USE Lexer_Core, ONLY: Init_Lexer
      USE MLSPCF1, ONLY: mlspcf_l1cf_start
      USE Output_m, ONLY: prunit
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
              & "Could not open L1 Config file" // physicalFilename)
      END IF

      CALL MLSMessage (MLSMSG_Info, ModuleName, &
           & "Opened L1 Config file" // physicalFilename)

!! Initialize the lexer, symbol table and tree checker's tables:

      CALL Init_Lexer (n_chars=10000, n_symbols=1000, hash_table_size=1003)
      CALL Allocate_Decl (ndecls=1000)
      CALL Allocate_Tree (n_tree=10000)
      CALL Init_tables
      CALL Init_units (lit_indices)

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

            !CALL Set_calibration (son)

         CASE DEFAULT

            CALL MLSMessage (MLSMSG_Error, ModuleName, 'Unknown section')

         END SELECT

      ENDDO

!! Will get the User Inputs from the L1CF file here

      MIFsGHz = 120
      MIFsTHz = 114

    END SUBROUTINE GetL1Config

    SUBROUTINE Set_globalsettings (root)

      use INIT_TABLES_MODULE, ONLY: p_output_version_string, p_version_comment
      use STRING_TABLE, ONLY: Get_string
      use TREE, ONLY: Decoration, Nsons, Subtree, Sub_rosa

      INTEGER :: root

      INTEGER :: i, son

      CHARACTER(LEN=80) :: line

      DO i = 2, Nsons (root) - 1

         son = Subtree (i, root)

         SELECT CASE (Decoration (Subtree (1,son)))

         CASE (p_output_version_string)

            CALL Get_string (Sub_rosa (Subtree(2,son)), &
                 L1Config%Globals%OutputVersionString)

         CASE (p_version_comment)

            CALL Get_string (Sub_rosa (Subtree(2,son)), &
                 L1Config%Globals%VersionComment)

         END SELECT


      ENDDO

    END SUBROUTINE Set_globalsettings

END MODULE MLSL1Config

! $Log$
! Revision 2.1  2001/02/23 20:50:54  perun
! Version 0.5 commit
!
