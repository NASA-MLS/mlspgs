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
MODULE CalibWeightsFlags
!=============================================================================

  USE MLSL1Common, ONLY: MaxMIFs, R8, DACSnum, FileNameLen
  USE Calibration, ONLY: WeightsFlags_T


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: DetermineWeightsFlags

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

  TYPE (WeightsFlags_T), DIMENSION(:), POINTER :: WeightsFlags

  INTEGER :: nom_MIFs, TotalMAFs

CONTAINS

!=============================================================================
  SUBROUTINE ProcessMAFdata
!=============================================================================

    ! <whd>Read the data which was assembled by ExamineData ->
    ! Examine{Sci,Eng}Data and process it. These data are read from two
    ! temporary files named {sci,eng}MAF_tmp.dat (PCFids=92{0,1}).
    ! 
    !</whd>


    USE MLSL1Config, ONLY: GHz_seq
    USE MLSL1Common, ONLY: L1BFileInfo
    USE EngTbls, ONLY: EngMAF, EngPkt
    USE L0_sci_tbls, ONLY: SciMAF
    USE SciUtils, ONLY: SwMirPos, GetScAngles
    USE MLSL1Utils, ONLY: QNan
    USE DACsUtils, ONLY: TPz
    USE L1LogUtils, ONLY: MAF_dur
    USE SDPToolkit, ONLY: PGS_TD_TAItoUTC



    INTEGER :: i, ios, sci_MAFno, n, sum_S, sum_T
    INTEGER :: EngMAF_unit, SciMAF_unit, MAF_data_unit
    INTEGER :: mindx(0:(MaxMIFs-1)), ngood(DACSnum)
    INTEGER, PARAMETER :: GHz_sw_indx(0:(MaxMIFs-1)) = (/ (i, i=1,MaxMIFs) /)
    CHARACTER(len=1) :: GHz_sw_pos(0:(MaxMIFs-1))
    CHARACTER(len=70) :: PLL_DN(0:(MaxMIFs-1)), PLL_dflt_DN(0:(MaxMIFs-1))
    REAL :: LLO_EU(16,0:(MaxMIFs-1)), LLO_dflt_EU(16,0:(MaxMIFs-1))
    INTEGER, PARAMETER :: MaxMIFno = (MaxMIFs - 1)
    REAL :: swFac(0:MaxMIFno)
    REAL(r8) :: TP_ana(0:MaxMIFno,DACSnum), TP_dig(0:MaxMIFno,DACSnum), &
         engTAI, sciTAI
    REAL(r8) :: Sum_ana(DACSnum) = 0.0, Sum_dig(DACSnum) = 0.0, &
         Sum_dig_dig(DACSnum) = 0.0, Sum_dig_ana(DACSnum) = 0.0
    CHARACTER(len=1), PARAMETER :: discard = "D"
    REAL :: TP_dig_max = 4.0     ! Maximum allowable TP digital

    CHARACTER(len=27) :: asciiUTC
    CHARACTER (LEN=FileNameLen) :: msg

    real(8) :: TAI
    integer :: status

    TYPE Sci_pos_T
       REAL :: APE(2,0:(MaxMIFs-1))
       REAL :: ASE(2,0:(MaxMIFs-1))
       REAL :: GME(2,0:(MaxMIFs-1))
       REAL :: TSE(2,0:(MaxMIFs-1))
    END TYPE Sci_pos_T

    TYPE (Sci_pos_T) :: Sci_pos, Sci_dflt_pos

    PRINT *, ' Reading MAFdata...'

! Init default positions:

    Sci_dflt_pos%APE = QNan()
    Sci_dflt_pos%ASE = QNan()
    Sci_dflt_pos%GME = QNan()
    Sci_dflt_pos%TSE = QNan()

    DO n = 0, (MaxMIFs - 1)
       PLL_dflt_DN(n) = CHAR(0)
    ENDDO

    LLO_dflt_EU = QNan()

    MAF_data_unit = L1BFileInfo%MAF_data_unit
    EngMAF_unit = L1BFileInfo%EngMAF_unit
    SciMAF_unit = L1BFileInfo%SciMAF_unit

    ! rewind to the beginning of the file each time through
    REWIND (engMAF_unit)
    REWIND (sciMAF_unit)


! Determine override sequence indexes to use for comparisons

    !in IDL, this is mindx=where(GHz_seq eq 'S', sum_S) GHz_seq is a character
    !string, 147 characters long, with a 'D' for any MIF that should be
    !discarded, a 'L' with a limb scan in it, a 'S' or 'T' for space and
    !calibration target, respectively. 

    !GHz_seq => L1Config%Calib%GHz_seq. This is read from the .cf file and tells
    !what MIFs are what type: "L" = limb, 'S' = space, 'T' = target and 'D' =
    !discard. 

    !GHz_seq_use => L1Config%Calib%GHz_seq_use: also read from .cf file, tells
    !you what to do with them. If 'M'atch then the telemetry must match what's
    !in the .cf. If 'O'verride, the data will override the .cf file.

    ! These two quantities are *nominal* values! *NOT* what's in the
    ! telemetry. That's checked later.

    ! Find and count the number of 'S'pace and 'T'arget MIFs
    mindx = 0
    WHERE (GHz_seq == "S")
       mindx = GHz_sw_indx
    ENDWHERE
    sum_S = SUM (mindx)

    mindx = 0
    WHERE (GHz_seq == "T")
       mindx = GHz_sw_indx
    ENDWHERE
    sum_T = SUM (mindx)

    ngood = 0
    i = 0

    ! <whd:comment>
    ! Loop over the data in the two temporary files, engMAF_tmp.dat and
    ! sciMAF_tmp.dat. These two loops try to match eng and sci MAFs and write
    ! data to a permanent file described before, which I'll call the L1BENG
    ! file. The outer reads an engineering MAF, the inner tries to find its
    ! matching science MAF. At each step it writes data to the L1BENG
    ! file. There are three cases: 1) The MAF numbers agree and the eng MAF is
    ! within 0.5 and 1.5 MAF durations away from the science MAF time; 2) the
    ! MAF numbers don't agree (most likely) or the times are out of range, and
    ! the engineering time is <= science time; 3) the MAF numbers don't agree
    ! (or the time is out of spec) and engineering time is > science time. In 1)
    ! it declares a match, writes some science data to the permanent file, exits
    ! the inner loop and reads the next eng MAF and exits the inner loop. In 2)
    ! is writes some default science values to L1BENG file and reads the next
    ! Eng MAF and writes _that_ data to the L1BENG file and continues to the
    ! next iteration of the inner loop (which will read the next science
    ! MAF). In 3) it reads the next science MAF, calculates position info and
    ! continues to next iteration of the inner loop.

    ! When a match has been found, the following is calculated.  
    !
    ! 1) The switching mirror 'type' (i.e. which MIFs are 'L', 'S', 'T' or 'D'
    ! MIFs)
    ! 2) The DACS TPs
    ! 3) The scan angles

    ! 4) Calibration weight flags.(See the WeightsFlags user type)

    !    4.1) These flags tell whether the weights for entire MAF has to be
    !    recalculated because the number of MIFs has changed. (recomp_MAF)
    ! 
    !    4.2) Whether the MIFs declared to be 'S'pace MIFs need to be
    !    recalculated because there's a different number of Space MIFs than the
    !    previous MAF (recomp_S)
    !
    !    4.3) Whether the MIFs declared to be 'T'arget MIFs need to be
    !    recalculated because there's a different number of Target MIFs than the
    !    previous MAF (recomp_T)
    !
    !  5) Calculates the TPz values (for MLSL1G)
    ! </whd:comment>



    outer: DO

       ! <whd:comment> Read the engineering data from the temporary file.
       !
       ! engMAF_unit is the temporary file with PCF id 921, nominally named
       ! engMAF_tmp.dat. This file is written by L1LogUtils::ExamineEngData
       ! EngMAF is a user defined type defined in EngTbls, ~line 155. EngPkt is
       ! a 6*char*256 char*1 array.
       ! </whd:comment>

       READ (engMAF_unit, iostat=ios) EngMAF, EngPkt
       engTAI = EngMAF%secTAI
       IF (ios /= 0) EXIT outer

       !<whd:comment>
       ! L1BFileInfo%EngId points to a fortran 90 binary file with PCF id
       ! 30003, nominally named MLS-Aura_L1BENG_<version>_<auraday>.dat. As far
       ! as I can see, it's only for diagnostic purposes and isn't used by any
       ! other part of the Level 1 code. This file written as a fortran
       ! unformated sequential file. EngMAF%Eng%value is the EU values for the
       ! counts contained in EngPkt. EngPkt is just the 6 256 byte arrays, as
       ! read from telemetry. EngMAF contains the converted data. It's a user
       ! type (Eng_MAF_T) defined in EngTbls: ~line 150. EngMAF%Eng%Value stores
       ! the converted counts. Unfortunately, the mnemonic isn't stored in this
       ! file, so you'll have to reconstruct that yourself. See
       ! EngUtils::ConvertCounts to help you.
       ! </whd:comment>


       WRITE (L1BFileInfo%EngId, iostat=ios) EngPkt, EngMAF%Eng%value

       ! <whd:comment> Read the Science data.  sciMAF_unit is the temporary file
       ! with PCF id 922, nominally named sciMAF_tmp.dat. Thie file is written
       ! by MLSL1Log->L1LogUtils::ExamineSciData SciMAF is an array
       ! (0:MaxMIFs-1) of the user defined type (Sci_pkt_t) defined in
       ! L0_sci_tbls.  
       ! </whd:comment>

       READ (sciMAF_unit, iostat=ios) SciMAF
       sciTAI = sciMAF(1)%secTAI
       IF (ios /= 0) EXIT outer
       sci_MAFno = SciMAF(0)%MAFno

       ! store the positions of various items as read from the
       ! telementry. (APE=Atenna Position Encoder), ASA is ???, GSM = GHz
       ! switching mirror, TSSM = THz switching, scanning mirror
       ! each of these sub-fields are (2,MaxMIFs)
       ! LLO_EU is (16,MaxMIFs)
       DO n = 0, (MaxMIFs - 1)
          Sci_pos%APE(:,n) = SciMAF(n)%APE_pos
          Sci_pos%ASE(:,n) = SciMAF(n)%ASA_pos
          Sci_pos%GME(:,n) = SciMAF(n)%GSM_pos
          Sci_pos%TSE(:,n) = SciMAF(n)%TSSM_pos
          PLL_DN(n) = SciMAF(n)%PLL_DN
          LLO_EU(:,n) = SciMAF(n)%LLO_EU
       ENDDO

       DO ! inner, find the Science MAF that goes with the Engineering MAF

          IF (EngMAF%MAFno == sci_MAFno .AND. &
               NINT(ABS(sciTAI - engTAI) / MAF_dur) == 1) THEN 

             ! vp original comment: Same MAF data.

             ! <whd>:hmmm... This code hits when the MAFno agree and engMAF time
             ! is between 0.5 and 1.5 MAFs away from the sciMAF, whichever
             ! direction. Dominick thinks this makes sense, it could be that the
             ! engineering data for a particular MAF is reported somewhere
             ! around 1 MAF before or after the science data. When I go through
             ! this with the debugger, it's almost always the case that the Eng
             ! data is first and the science data with the same MAF number comes
             ! about 1 MAF later.
             ! </whd>

             ! Write the positions, the PLL and LLO data. 
             WRITE (L1BFileInfo%EngId, iostat=ios) Sci_pos
             WRITE (L1BFileInfo%EngId, iostat=ios) PLL_DN
             WRITE (L1BFileInfo%EngId, iostat=ios) LLO_EU
             EXIT   !Don't need to go any further with this inner loop
          ENDIF

          IF (engTAI <= sciTAI) THEN ! vp original comment: Catch the Sci
             
             ! <whd>:observation: It's my observation that engTAI < sciTAI is
             ! almost always true when the MAF nums agree. This branch will only
             ! hit when the block above is passed. That can happen when
             ! engMAF%MAFno != sci_MAFno, or the times are less than 0.5 or more
             ! than 1.5 MAF durations away, or both. So, I guess what's
             ! happening here is that the code thinks the file pointer in the
             ! science temp file is behind by at least 1 MAF, so it writes some
             ! default values (to match the eng MAF it's already written?) and
             ! reads the next Eng MAF and tries again.
             !
             ! sci is ahead, write some default values, then read the next
             !Eng MAF
             ! </whd>

             WRITE (L1BFileInfo%EngId, iostat=ios) Sci_dflt_pos
             WRITE (L1BFileInfo%EngId, iostat=ios) PLL_dflt_DN
             WRITE (L1BFileInfo%EngId, iostat=ios) LLO_dflt_EU
             READ (engMAF_unit, iostat=ios) EngMAF, EngPkt 
             IF (ios /= 0) EXIT outer
             engTAI = EngMAF%secTAI
             WRITE (L1BFileInfo%EngId, iostat=ios) EngPkt, EngMAF%Eng%value

          ELSE ! vp original comment: Catch the Eng

             ! <whd> science is before eng. read the next Science MAF, store the
             ! newly calculated position info and try again.</whd>


             READ (sciMAF_unit, iostat=ios) SciMAF
             sciTAI = sciMAF(2)%secTAI
             IF (ios /= 0) EXIT outer
             sci_MAFno = SciMAF(0)%MAFno
             DO n = 0, (MaxMIFs - 1)
                Sci_pos%APE(:,n) = SciMAF(n)%APE_pos
                Sci_pos%ASE(:,n) = SciMAF(n)%ASA_pos
                Sci_pos%GME(:,n) = SciMAF(n)%GSM_pos
                Sci_pos%TSE(:,n) = SciMAF(n)%TSSM_pos
             ENDDO

          ENDIF
       ENDDO ! inner loop, matching sci/eng MAFs

! Reset GHz SW mirror pos:

       ! Here we find what the data says about the GHz swiching mirror
       DO n = 0, (MaxMIFs - 1)
          SciMAF(n)%GHz_sw_pos = SwMirPos ("G", SciMAF(n)%GSM_theta)
          GHz_sw_pos(n) = SciMAF(n)%GHz_sw_pos
       ENDDO

! Accumulate DACS TPs to calculate daily TP zeros:

       swFac = 1.0        ! init switch factor to non-discard MIFs
       WHERE (SciMAF%GHz_sw_pos == discard)   ! Don't use discards
          swFac = 0.0
       ENDWHERE
!       ngood = ngood + COUNT (swFac == 1.0)
       DO n = 1, DACSNUM
          TP_ana(:,n) = SciMAF%TP(n) * swFac
          TP_dig(:,n) = (SciMAF%TPdigP(n) + SciMAF%TPdigN(n)) * swFac
          WHERE (TP_dig(:,n) > TP_dig_max)
             TP_dig(:,n) = 0.0
             TP_ana(:,n) = 0.0
          ENDWHERE
!##########################################################

 

! DEBUG: Fix accumulation of sum_ana to be only over values where TP_dig /= 0.0, consistent

! with the other sums, so the mean of TP_ana can be calculated as Sum_ana(i)/ngood(i)         

          WHERE (TP_dig(:,n) == 0.0)

             TP_ana(:,n) = 0.0

          ENDWHERE        

 

!##########################################################
          ngood(n) = ngood(n) + COUNT (TP_dig(:,n) /= 0.0)
       ENDDO
       DO n = 1, DACSNUM
          Sum_dig(n) = Sum_dig(n) + SUM (TP_dig(:,n))
          Sum_ana(n) = Sum_ana(n) + SUM (TP_ana(:,n))
          Sum_dig_dig(n) = Sum_dig_dig(n) + SUM (TP_dig(:,n)*TP_dig(:,n))
          Sum_dig_ana(n) = Sum_dig_ana(n) + SUM (TP_dig(:,n)*TP_ana(:,n))
       ENDDO

! Reset scAngles to be used for L1BOA

       CALL GetScAngles

       i = i + 1

       WeightsFlags(i)%MAFno = EngMAF%MAFno

! Check MAF duration. If the length of the MAF (in MIFs) changed, we have to
! recompute the whole MAF. Set the flag accordingly. This might come into play
! in Calibration::SetComVecs

       WeightsFlags%recomp_MAF  = (EngMAF%MIFsPerMAF /= nom_MIFs)

! Check switching mirror positions:

       !GHz_seq==`predicted value' from .cf
       !GHz_sw_pos=`telemetry value' 

       ! Mark MIFs that should be S or T, but which the telemetry says are out
       ! of range, i.e. 'D'iscards. Set a flag for when the S or T MIFs must be
       ! recomputed because they don't agree with the count from the last MAF

       mindx = 0
       WHERE (GHz_seq == "S" .AND. GHz_sw_pos /= "D")
          mindx = GHz_sw_indx
       ENDWHERE
       WeightsFlags%recomp_S  = (SUM (mindx) /= sum_S)

       mindx = 0
       WHERE (GHz_seq == "T" .AND. GHz_sw_pos /= "D")
          mindx = GHz_sw_indx
       ENDWHERE
       WeightsFlags%recomp_T  = (SUM (mindx) /= sum_T)

       PRINT *, 'i, MAFno, GHz_sw_pos: ', i, EngMAF%MAFno, SciMAF%GHz_sw_pos
       PRINT *, 'Calib Weights flags: ', WeightsFlags(i)
       if (WeightsFlags(i)%recomp_S    .OR. &
            & WeightsFlags(i)%recomp_T .OR.&
            & WeightsFlags(i)%recomp_MAF) THEN 
          status=PGS_TD_TAItoUTC(SciMAF(0)%secTAI,asciiUTC)
          msg='Weights recomp needed at ' // asciiUTC
          
          IF (WeightsFlags(i)%recomp_S)   msg = trim(msg) // ': Space'
          IF (WeightsFlags(i)%recomp_T)   msg = trim(msg) // ': Target' 
          IF (WeightsFlags(i)%recomp_MAF) msg = trim(msg) // ': MAF'
          PRINT *,trim(msg)
          print *,"Expected values:",GHz_seq
       ENDIF

       ! write data to the file pointed to by PCF id=922, nominally named
       ! MAF_data_tmp.dat

       WRITE (MAF_data_unit) WeightsFlags(i)
       WRITE (MAF_data_unit) EngMAF
       WRITE (MAF_data_unit) SciMAF
    ENDDO outer

    TotalMAFs = i
    PRINT *, 'MAFs: ', TotalMAFs
    ENDFILE MAF_data_unit

! Calculate TPz values to pass on the MLSL1G program:

    DO i = 1, DACSNUM
       TPz(i) = ((Sum_dig_dig(i) / ngood(i) * Sum_ana(i) / ngood(i)) - &
            (Sum_dig_ana(i) / ngood(i) * Sum_dig(i) / ngood(i))) / &
            (Sum_dig_dig(i) / ngood(i) - (Sum_dig(i) / ngood(i))**2)
    ENDDO

  END SUBROUTINE ProcessMAFdata

!=============================================================================
  SUBROUTINE DetermineWeightsFlags
!=============================================================================

    ! Pretty much a wrapper to ProcessMAFdata

    USE MLSL1Config, ONLY: L1Config
    USE L1LogUtils, ONLY: EngMAFs, SciMAFs

    INTEGER :: maxMAFs

    PRINT *, 'Determining weights flags...'

    nom_MIFs = L1Config%Calib%MIFsPerMAF

    maxMAFs = MAX (EngMAFs, SciMAFs)
    ALLOCATE (WeightsFlags(maxMAFs))

    CALL ProcessMAFdata

  END SUBROUTINE DetermineWeightsFlags

!=============================================================================
  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here
END MODULE CalibWeightsFlags
!=============================================================================
! $Log$
! Revision 2.13  2023/03/01 18:30:42  pwagner
! Fixed bug in calibrating DACS
!
! Revision 2.12  2018/04/09 22:12:29  whdaffer
! Documentation
!
! Revision 2.11  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.10.4.1  2015/10/09 10:21:37  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.10  2011/01/27 15:33:53  perun
! Change TAI time from MIF 2 to MIF 1
!
! Revision 2.9  2007/06/21 20:57:53  perun
! Exclude TP_dig greater than 4.0 in calculating TPz
!
! Revision 2.8  2006/09/26 16:01:33  perun
! Correct testing to insure alignment of Eng and Sci data
!
! Revision 2.7  2006/08/02 18:52:17  perun
! Accumulate DACS TPs to calculate daily TP zeros
!
! Revision 2.6  2006/04/05 18:10:08  perun
! Remove unused variables
!
! Revision 2.5  2005/06/23 18:41:35  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.4  2004/08/12 13:51:49  perun
! Version 1.44 commit
!
! Revision 2.3  2004/01/09 17:46:22  perun
! Version 1.4 commit
!
! Revision 2.2  2003/09/15 17:15:53  perun
! Version 1.3 commit
!
! Revision 2.1  2003/08/15 14:54:55  perun
! Version 1.2 commit
!
! Revision 2.1  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
!
