! Copyright 2006, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
MODULE MLSL1Config  ! Level 1 Configuration
!=============================================================================

  USE MLSCommon, ONLY: TAI93_Range_T
  USE MLSL1Common, ONLY: MaxMIFs, NumBands, BandChanBad
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Info
  USE Init_tables_module, ONLY: First_Parm, Last_Parm
  USE INTRINSIC, ONLY: parm_indices
  USE Output_m, ONLY: OutputOptions, Output, BOTHPRUNIT, STDOUTPRUNIT
  USE STRING_TABLE, ONLY: Get_string

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: L1Config_T, L1Config, Globals_T, Calib_T, Output_T, GHz_seq, &
       GHz_seq_use, THz_seq, THz_seq_use, MIFsGHz, MIFsTHz
  PUBLIC :: GetL1Config

  INTEGER, PARAMETER :: MIFsGHz = 125   ! Length (MIFs) of GHz module data
  INTEGER, PARAMETER :: MIFsTHz = 125   ! Length (MIFs) of THz module data

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

  TYPE Globals_T
     CHARACTER(LEN=80) :: OutputVersionString
     LOGICAL :: ProduceL1BOA = .TRUE.
     LOGICAL :: SimOA = .FALSE.
  END TYPE Globals_T

  TYPE Calib_T
     INTEGER :: CalWindow
     INTEGER :: MIFsPerMAF
     INTEGER :: MAFexpandNum             ! number of MAFs to expand at both ends
     INTEGER :: MaxDataGaps              ! maximum allowed data gaps
     INTEGER :: MaxErroneousCounterMAFs  ! maximum error logs count
     INTEGER :: DiffBeginEndEng          ! Beging and End for Engineering times difference
     INTEGER :: MinSpaceLimbs            ! minimum "Space" views per MAF
     INTEGER :: DACSwindow
     REAL :: GHzSpaceTemp, GHzTargetTemp
     REAL :: THzSpaceTemp, THzTargetTemp
     REAL :: THzSpaceAngle, THzMaxBias
     REAL :: MIF_duration, MIF_DeadTime
     REAL :: MoonToSpaceAngle
     REAL :: AntOffsetsScale = 1.0       ! scale factor for antenna offsets
     LOGICAL :: UseDefaultGains = .FALSE.
     LOGICAL :: CalibDACS = .TRUE.
     LOGICAL :: TPdigital = .TRUE.
     LOGICAL :: THzColdCal = .TRUE.
     LOGICAL :: Do_Slimb = .FALSE.
     CHARACTER(LEN=1) :: GHz_seq(0:MaxMIFs-1), THz_seq(0:MaxMIFs-1)
     CHARACTER(LEN=1) :: GHz_seq_use, THz_seq_use
  END TYPE Calib_T

  TYPE Output_T
     LOGICAL :: RemoveBaseline = .TRUE.             ! For GHz Baseline removal
     LOGICAL :: DeconvolveDACS = .FALSE.            ! For DACS deconvolution
     LOGICAL :: DisableRadOut(NumBands) = .FALSE.   ! To output band radiances
     LOGICAL :: SubtractBinnedBaseline(NumBands) = .FALSE. ! To adjust baseline
     LOGICAL :: EnableChi2Err(NumBands) = .FALSE.   ! For RadErr calculation
     
     ! Limits on certain logged messages
     !        'Found Attenuation WALL at MAF time'
     ! (not yet configurable by the user)
     integer :: MaxAttenuationWalls             = 10000 

     integer :: NumAttenuationWalls             = 0
     
     ! Keep a close eye on UpdateCalWindow in case of .. trouble
     logical :: DebugUpdateCalWindow            = .false.
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

      !! Opens and parses the L1CF file.

      USE Declaration_Table, ONLY: Allocate_Decl
      USE Init_tables_module, ONLY: Init_tables, z_globalsettings, &
           z_calibration, z_output
      USE Lexer_Core, ONLY: Init_Lexer
      USE MLSPCF1, ONLY: mlspcf_l1cf_start
      USE Parser, ONLY: Clean_Up_Parser, Configuration
      USE Parser_Table_m, ONLY:  Destroy_Parser_Table, Parser_Table_t
      USE Parser_Tables_L2CF, ONLY: Init_Parser_Table
      USE SDPToolkit, ONLY: PGS_PC_GetReference, PGS_S_SUCCESS, &
           PGSd_IO_Gen_RSeqFrm, PGS_IO_Gen_openF, PGS_IO_Gen_closeF
      USE STRING_TABLE, ONLY: AddInUnit
      USE Tree, ONLY: Allocate_Tree, Decoration, Nsons, Subtree
      USE Tree_checker, ONLY: Check_tree

      CHARACTER (LEN=132) :: physicalFilename
      INTEGER :: i, returnStatus, version
      INTEGER :: error           ! error return from tree checker
      INTEGER :: first_section   ! index of son of root of first n_cf_node
      INTEGER :: root            ! of the abstract syntax tree
      INTEGER :: son             ! of root
      integer :: l1cf_unit
      type(Parser_Table_t) :: Parser_Table

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
!! Initialize the lexer, symbol table and tree checker's tables:

      CALL Init_Lexer (n_chars=10000, n_symbols=1000, hash_table_size=1003)
      CALL Allocate_Decl (ndecls=1000)
      CALL Allocate_Tree (n_tree=10000)
      CALL Init_tables

      call AddInUnit(l1cf_unit)
      CALL MLSMessage (MLSMSG_Info, ModuleName, &
           & "Opened L1 Config file: " // physicalFilename)

!! Produce the abstract syntax tree

      CALL init_parser_table ( parser_table )
      CALL Configuration (root, parser_table)
      CALL destroy_parser_table ( parser_table )
      CALL clean_up_parser

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

      ! This may be causing problems, so I'm resetting it to STDOUTPRUNIT
      !OutputOptions%prunit = BOTHPRUNIT  ! to output to MLSMessage (and terminal)
      OutputOptions%prunit = STDOUTPRUNIT

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

      USE INIT_TABLES_MODULE, ONLY: p_output_version_string, p_produce_l1boa, &
           p_simoa
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
      USE INIT_TABLES_MODULE, ONLY: p_removebaseline, p_deconvolveDACS, &
           s_chi2err, f_bandno, s_subtractbinnedbaseline, s_disableradout
      USE TREE, ONLY: Decoration, Nsons, Subtree, Node_id
      USE TREE_TYPES
      USE MoreTree, ONLY: Get_Boolean

      INTEGER :: root, i, j, k, son, key, spec
      INTEGER :: expr_units(2)
      DOUBLE PRECISION :: expr_value(2)

      DO i = 2, Nsons (root) - 1

         son = Subtree (i, root)

         SELECT CASE (node_id (son))

         CASE (n_equal)

            GotParm(Decoration (Subtree (1,son))) = .TRUE.

            SELECT CASE (Decoration (Subtree (1,son)))

            CASE (p_removebaseline)

               L1Config%Output%RemoveBaseline = Get_Boolean (son)

            CASE (p_deconvolveDACS)

               L1Config%Output%DeconvolveDACS = Get_Boolean (son)

            END SELECT

         CASE (n_spec_args)

            key = son
            spec = decoration (subtree (1, decoration(subtree (1, key))))

            SELECT CASE (spec)

            CASE (s_chi2err, s_subtractbinnedbaseline, s_disableradout)

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
                        IF (spec == s_chi2err) THEN
                           L1Config%Output%EnableChi2Err(INT(expr_value(1)): &
                            INT(expr_value(2))) = .TRUE.
                        ELSE IF (spec == s_disableradout) THEN
                           L1Config%Output%DisableRadOut(INT(expr_value(1)): &
                            INT(expr_value(2))) = .TRUE.
                        ELSE IF (spec == s_subtractbinnedbaseline) THEN
                           L1Config%Output%SubtractBinnedBaseline &
                            (INT(expr_value(1)): INT(expr_value(2))) = .TRUE.
                        ENDIF
                     ENDDO

                  END SELECT

               ENDDO

            END SELECT

         END SELECT

      ENDDO

      IF (ANY(L1Config%Output%EnableChi2Err)) THEN
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

      !parse the 'Calibration' section of the .cf file.


      USE EXPR_M, ONLY: Expr
      USE INIT_TABLES_MODULE, ONLY: p_calwindow, s_spaceMIFs, s_targetMIFs, &
           s_limbMIFs, s_discardMIFs, f_mifs, f_use, l_match, l_override, &
           f_module, f_secondary, p_usedefaultgains, p_GHzSpaceTemp, &
           p_GHzTargetTemp, p_THzSpaceTemp, p_THzTargetTemp, p_mif_duration, &
           p_mif_dead_time, p_mifspermaf, p_calibDACS, p_THzMaxBias, &
           p_thzspaceangle, f_bandno, f_chan, s_markchanbad, p_thzcoldcal, &
           p_MoonToSpaceAngle, p_DACSwindow, p_UseAntOffsets, p_MinSpaceLimbs, &
           p_MAFexpandNum, p_MaxDataGaps, p_MaxErroneousCounterMAFs, p_DiffBeginEndEng, p_TPdigital, p_Do_Slimb, f_yrdoy
      USE BrightObjects_m, ONLY: s_BrightObject, f_angle, f_name, f_negate, &
           l_mercury, BO_Angle_GHz, BO_Angle_THz, BO_NumGHz, BO_NumTHz, &
           BO_Index_GHz, BO_Index_THz, BO_Negate_GHz, BO_Negate_THz
      USE INTRINSIC, ONLY: l_ghz, l_thz, phyq_mafs, phyq_temperature, &
           phyq_mifs, phyq_time, phyq_angle
      USE TREE, ONLY: Decoration, Nsons, Subtree, Sub_rosa, Node_id
      USE TREE_TYPES
      USE MoreTree, ONLY: Get_Boolean
      USE InitPCFs, ONLY: L1PCF
      USE SDPToolkit, ONLY: PGS_S_SUCCESS

      INTEGER :: root

      CHARACTER(LEN=1), POINTER, DIMENSION(:) :: scan_seq
      CHARACTER(LEN=1), POINTER :: scan_use
      CHARACTER(LEN=80) :: identifier
      CHARACTER(LEN=17) :: chanbad_UTC = '20YY-DOYT00:00:00'
      INTEGER :: i, j, k, son, key, spec, stat, bandno, channo, yrdoy
      INTEGER :: expr_units(2), BO_index
      REAL :: BO_angle
      DOUBLE PRECISION :: expr_value(2), start_TAI, chanbad_TAI
      LOGICAL :: GHz_mod, sec_tgt, Negate

      ! This type is used to store calibration directives from the .cf file. The
      ! TYPE= 'L', 'S', 'T' 'D', or 't'. 't' means 'secondary target' -- not
      ! listed in the user type documentatin here, but you see it in the code --
      ! however, this would only come into play if 1) it's the GHz module and
      ! 2)the line defining the MIFs that are of a particular type had the
      ! declaration `secondary' in them. No extant .cf file has this
      ! characteristic, and Paul Wagner assures me that no .cf file in the past,
      ! nor any in the future will.

      ! if USE='M'atch, the telemetry must agree with what's in the .cf ifle
      TYPE Scan_T
         CHARACTER(LEN=1) :: USE    ! 'M'atch or 'O'verride, 
         CHARACTER(LEN=1) :: TYPE   ! 'L'imb, 'S'pace, 'T'arget, or 'D'iscard
         INTEGER :: MIF(0:MaxMIFs-1) ! The MIFs numbers that are of type `TYPE'
      END TYPE Scan_T

      TYPE (Scan_T) :: scan

      CHARACTER(LEN=1), PARAMETER :: ignore = "I", unknown = "U", overlap = "O"

      INTEGER, EXTERNAL :: PGS_TD_UTCtoTAI

! Current input time for comparisons:

      stat = PGS_TD_UTCtoTAI (L1PCF%startUTC, start_TAI)

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

               ! remove the 'p_' and look for the remaining string in the .cf
               ! file
            CASE (p_calwindow)

               CALL Expr (subtree (2, son), expr_units, expr_value)
               L1Config%Calib%CalWindow = expr_value(1)
               IF (expr_units(1) /= phyq_mafs) THEN
                  CALL Get_string (Sub_rosa (Subtree(1,son)), identifier)
                  CALL MLSMessage (MLSMSG_Error, ModuleName, &
                       TRIM (identifier)//' is not input as MAFs')
               ENDIF

            CASE (p_MAFexpandNum)

               CALL Expr (subtree (2, son), expr_units, expr_value)
               L1Config%Calib%MAFexpandNum = expr_value(1)
               IF (expr_units(1) /= phyq_mafs) THEN
                  CALL Get_string (Sub_rosa (Subtree(1,son)), identifier)
                  CALL MLSMessage (MLSMSG_Error, ModuleName, &
                       TRIM (identifier)//' is not input as MAFs')
               ENDIF
	       
            CASE (p_MaxDataGaps)

               CALL Expr (subtree (2, son), expr_units, expr_value)
               L1Config%Calib%MaxDataGaps = expr_value(1)
          
	    CASE (p_MaxErroneousCounterMAFs)

               CALL Expr (subtree (2, son), expr_units, expr_value)
               L1Config%Calib%MaxErroneousCounterMAFs = expr_value(1)
	       
	    CASE (p_DiffBeginEndEng)

               CALL Expr (subtree (2, son), expr_units, expr_value)
               L1Config%Calib%DiffBeginEndEng = expr_value(1)           

            CASE (p_MinSpaceLimbs)

               CALL Expr (subtree (2, son), expr_units, expr_value)
               L1Config%Calib%MinSpaceLimbs = expr_value(1)
               IF (expr_units(1) /= phyq_mifs) THEN
                  CALL Get_string (Sub_rosa (Subtree(1,son)), identifier)
                  CALL MLSMessage (MLSMSG_Error, ModuleName, &
                       TRIM (identifier)//' is not input as MIFs')
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

            CASE (p_THzColdCal)

               L1Config%Calib%THzColdCal = Get_Boolean (son)

            CASE (p_MoonToSpaceAngle)

               CALL Expr (subtree (2, son), expr_units, expr_value)
               L1Config%Calib%MoonToSpaceAngle = expr_value(1)
               IF (expr_units(1) /= phyq_angle) THEN
                  CALL Get_string (Sub_rosa (Subtree(1,son)), identifier)
                  CALL MLSMessage (MLSMSG_Error, ModuleName, &
                       TRIM (identifier)//' is not input as deg[rees]')
               ENDIF

            CASE (p_dacswindow)

               CALL Expr (subtree (2, son), expr_units, expr_value)
               L1Config%Calib%DACSwindow = expr_value(1)
               IF (expr_units(1) /= phyq_mafs) THEN
                  CALL Get_string (Sub_rosa (Subtree(1,son)), identifier)
                  CALL MLSMessage (MLSMSG_Error, ModuleName, &
                       TRIM (identifier)//' is not input as MAFs')
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

            CASE (p_UseAntOffsets)

               IF (Get_Boolean (son)) THEN
                  L1Config%Calib%AntOffsetsScale = 1.0
               ELSE
                  L1Config%Calib%AntOffsetsScale = 0.0
               ENDIF

            CASE (p_calibDACS)

               L1Config%Calib%calibDACS = Get_Boolean (son)

            CASE (p_TPdigital)

               L1Config%Calib%TPdigital = Get_Boolean (son)

            CASE (p_Do_Slimb)

               L1Config%Calib%Do_Slimb = Get_Boolean (son)

            END SELECT

         CASE (n_spec_args)
            ! to handle constructs like 
            !  `limbMIFs, MIFs=[2:122], use=match, module=GHz'

            key = son
            spec = decoration (subtree (1, decoration(subtree (1, key))))

            SELECT CASE (spec)

            CASE (s_spaceMIFs, s_targetMIFs, s_limbMIFs, s_discardMIFs)

               SELECT CASE (spec)

               CASE (s_spaceMIFs)
                  scan%Type = "S"
               CASE (s_targetMIFs)
                  scan%Type = "T" ! though, can be changed to 't' if sec_tgt
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
                  ! if scan_use (i.e. L1Config%Calib%{G,T}Hz) == ignore, it
                  ! means it hasn't been set yet, so we set it to the value
                  ! we've just read from the .cf file.
                  scan_use = scan%Use
               ELSE IF (scan_use /= scan%Use) THEN
                  ! how can this possibly happen?!? 
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

            CASE (s_markchanbad)

               chanbad_TAI = 0.0d00   ! No current bad time

               DO j = 2, nsons (key)

                  son = subtree (j, key)

                  SELECT CASE (decoration (subtree(1,son)))   ! field

                  CASE (f_chan)

                     CALL Expr (subtree (2, son), expr_units, expr_value)
                     channo = expr_value(1)
                     
                  CASE (f_bandno)
                     CALL Expr (subtree (2, son), expr_units, expr_value)
                     bandno = expr_value(1)

                  CASE (f_yrdoy)

                     CALL Expr (subtree (2, son), expr_units, expr_value)
                     yrdoy = expr_value(1)
                     write (chanbad_UTC(3:4), '(I2.2)') (yrdoy / 1000)
                     write (chanbad_UTC(6:8), '(I3.3)') MOD (yrdoy, 1000)
                     stat = PGS_TD_UTCtoTAI (chanbad_UTC, chanbad_TAI)
                     IF (stat /= PGS_S_SUCCESS) THEN
                        CALL MLSMessage (MLSMSG_Error, ModuleName, &
                             'MarkChanBad YRDOY input field is incorrect!')
                     ENDIF

                   END SELECT

               ENDDO

               IF (bandno < 1 .OR. bandno > NumBands) THEN
                  CALL MLSMessage (MLSMSG_Error, ModuleName, &
                       'Bandno number out of range!')
               ENDIF
               IF (channo < 1 .OR. channo > BandChanBad%MaxChan(bandno)) THEN
                  CALL MLSMessage (MLSMSG_Error, ModuleName, &
                       'ChanNo number out of range!')
               ENDIF

               IF (start_TAI >= chanbad_TAI) &
                    BandChanBad%Sign(bandno,channo) = -1.0  ! Mark as "Bad"

            CASE (s_BrightObject)   ! Bright Objects

               Negate = .FALSE.  ! Initialize to do not negate
               DO j = 2, nsons (key)

                  son = subtree (j, key)

                  SELECT CASE (decoration (subtree(1,son)))   ! field

                  CASE (f_module)

                     SELECT CASE (decoration (subtree (2, son)))

                     CASE (l_ghz)
                        GHz_mod = .TRUE.
                     CASE (l_thz)
                        GHz_mod = .FALSE.
                     END SELECT

                  CASE(f_name)
                     CALL Expr (subtree (2, son), expr_units, expr_value)
                     BO_index = INT(expr_value(1)) - l_mercury + 1
                  CASE (f_angle)
                     CALL Expr (subtree (2, son), expr_units, expr_value)
                     BO_angle = expr_value(1)
                  CASE (f_negate)
                     Negate = Get_Boolean (son)
                  ENDSELECT

               ENDDO

               IF (GHz_mod) THEN   ! Save appropriate BO info
                  BO_NumGHz = BO_NumGHz + 1
                  BO_Index_GHz(BO_NumGHz) = BO_index
                  BO_Angle_GHz(BO_index) = BO_angle
                  BO_Negate_GHz(BO_NumGHz) = Negate
               ELSE
                  BO_NumTHz = BO_NumTHz + 1
                  BO_Index_THz(BO_NumTHz) = BO_index
                  BO_Angle_THz(BO_index) = BO_angle
                  BO_Negate_THz(BO_NumTHz) = Negate
                ENDIF

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

    END SUBROUTINE Set_calibration

  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here
END MODULE MLSL1Config

! $Log$
! Revision 2.37  2023/06/06 22:33:03  pwagner
! Reduce routine printing
!
! Revision 2.36  2018/04/09 22:15:02  whdaffer
! STDOUTPRUNIT
!
! Revision 2.35  2016/05/10 20:30:57  mmadatya
! To get the error-checking parameters from the l1 configuration file instead of them being hard-coded into the source code
!
! Revision 2.34  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.33.4.1  2015/10/09 10:21:38  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.33  2014/05/20 23:57:14  vsnyder
! New parser gets its tables from an argument instead of an include
!
! Revision 2.32  2010/08/06 17:53:43  pwagner
! Moved call to AddInUnit after init_lexer
!
! Revision 2.31  2010/05/23 04:13:48  honghanh
! Use AddInunit instead of inunit due to change in string_table
!
! Revision 2.30  2008/03/04 20:01:54  perun
! Use optional YRDOY field in MarkChanBad entry to determine when to mark channel data bad.
!
! Revision 2.29  2008/01/15 19:54:35  perun
! Add DisableRadOut to disable outputting unwanted bands.
!
! Revision 2.28  2007/02/09 15:05:24  perun
! Add Do_Slimb flag
!
! Revision 2.27  2006/09/28 16:15:37  perun
! Remove WriteDiagOffsets
!
! Revision 2.26  2006/08/02 19:23:50  pwagner
! prunit now a component of OutputOptions
!
! Revision 2.25  2006/08/02 18:55:22  perun
! Added SubtractBinnedBaseline field
!
! Revision 2.24  2006/06/14 13:47:00  perun
! Handle TPdigital input
!
! Revision 2.23  2006/04/05 18:10:57  perun
! Remove unused variables
!
! Revision 2.22  2006/03/24 15:12:19  perun
! Add MAFexpandNum, MinSpaceLimbs, THzColdCal, WriteDiagOffsets and remove Switch
!
! Revision 2.21  2005/12/06 19:27:12  perun
! Removed MoonToLimbAngles parsing and added Bright Object parsing
!
! Revision 2.20  2005/10/10 19:06:27  perun
! Add DeconvolveDACS field
!
! Revision 2.19  2005/06/23 18:41:35  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.18  2005/05/02 16:04:15  perun
! Add UseAntOffsets field and AntOffsetsScale factor
!
! Revision 2.17  2005/01/28 17:00:29  perun
! Split MoonToLimbAngle into GHz and THz
!
! Revision 2.16  2004/12/01 17:10:07  perun
! Remove VersionComment and add DACSwindow
!
! Revision 2.15  2004/11/10 15:39:17  perun
! Add case to set BandChanBad value based on user input
!
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
