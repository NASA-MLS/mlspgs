! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!============================================================================
PROGRAM MLSL2                   ! MLS Level 2 software
!============================================================================

  ! External modules

  USE MLSCommon
  USE MLSMessageModule
  USE L2GPData
  USE L2AUXData
  USE QuantityTemplates
  USE VectorsModule
  USE MLSCF
  USE OpenInit
  USE ScanDivide
  USE Construct
  USE Fill
  USE Join
  USE OutputAndClose

  IMPLICIT NONE

  !--------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=130) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !--------------------------------------------------------------------------

  ! This is the Level 2 processing software for the EOS Microwave Limb Sounder.
  ! For more information see:

  ! EOS MLS Retrieval Processes Algorithm Theoretical Basis Document
  !   JPL D-16159
  ! Functional Requirements for the EOS MLS Level 2 software
  !   JPL D-18027
  ! EOS MLS Level 2 file description document
  !   JPL D-18028
  ! EOS MLS Level 2 software users' guide
  !   JPL D-18029
  ! Design of the EOS MLS Level 2 software for version 0.1
  !   JPL D-18030.

  !--------------------------------------------------------------------------

  ! Local parameters

  !--------------------------------------------------------------------------

  ! Local variables roughly ordered increasing triviality

  TYPE (TAI93_Range_T) :: processingRange ! Data processing range
  TYPE (L1BInfo_T) :: l1bInfo   ! File handles etc. for L1B dataset
  TYPE (MLSCF_T) :: l2cf        ! The information from the l2cf file
  TYPE (GriddedData_T), DIMENSION(:), POINTER :: aprioriData 
  ! Input a priori database
  TYPE (MLSChunk_T), DIMENSION(:), POINTER :: chunks ! Chunks for data

  TYPE (QuantityTemplate_T), DIMENSION(:), POINTER :: qtyTemplates
  TYPE (QuantityTemplate_T), DIMENSION(:), POINTER :: mifGeolocation
  TYPE (VectorTemplate_T), DIMENSION(:), POINTER :: vectorTemplates
  TYPE (Vector_T), DIMENSION(:), POINTER :: vectors

  TYPE (L2GPData_T), DIMENSION(:), POINTER :: l2gpDatabase ! L2GP products
  TYPE (L2AUXData_T), DIMENSION(:), POINTER :: l2auxDatabase ! L2AUX products

  INTEGER :: chunkNo,instrumentModule, n            ! Loop counter
  CHARACTER (LEN=132) :: message ! Line of message output

  !--------------------------------------------------------------------------

  ! Executable code


  CALL MLSMessage(MLSMSG_Info,ModuleName, &
       & "EOS MLS Level 2 data processing commenced")

  ! Before we do anything we nullify the database pointers.  This is needed
  ! because apparently f90 does not set pointers to .not. associated when
  ! they are first created.

  NULLIFY (aprioriData,chunks,qtyTemplates,mifGeolocation,vectorTemplates,&
       & vectors,l2gpDatabase,l2auxDatabase)

  ! The first thing to do is to open the input Level 1 files, l2cf files and do
  ! other initialization stuff.


  CALL OpenAndInitialize(processingRange, l1bInfo, l2cf, &
       & aprioriData)

  ! The next stage is to divide the input L1B dataset into chunks.

  CALL ScanAndDivide(processingRange, l1bInfo, l2cf, chunks)
  ! Now we loop over the chunks and do the data processing

  n= SIZE(chunks)
  DO chunkNo=1,n
     WRITE (UNIT=message,FMT=*) "Processing chunk #",chunkNo
     CALL MLSMessage(MLSMSG_Info,ModuleName,TRIM(message))
     ! This will need to change a lot in 0.5 and above
     ! First construct the vector quantity templates and vector templates

     CALL MLSL2Construct(l2cf,l1bInfo,chunks(chunkNo),qtyTemplates, &
          & vectorTemplates,mifGeolocation)

     ! Now fill those templates with the a priori and/or l1b data
     CALL MLSL2Fill(l2cf,l1bInfo,aprioriData, vectorTemplates, vectors)

     ! In later versions there will be retrieve etc. and also repeats, but for
     ! the moment we'll leave it at this.

     CALL MLSL2Join(l2cf,vectors,l2gpDatabase,l2auxDatabase,chunks,chunkNo)
     ! At the end of each chunk, destroy the vectorTemplate and vector
     ! information.

     CALL DestroyVectorDatabase(vectors)
     CALL DestroyVectorTemplateDatabase(vectorTemplates)
     CALL DestroyQuantityTemplateDatabase(qtyTemplates)
     CALL DestroyQuantityTemplateDatabase(mifGeolocation)
     CALL DestroyGridTemplateDatabase(aprioriData)
  END DO
  ! Now we write out the data
 CALL Ouput_Close(l2cf,l2gpDatabase,l2auxDatabase)
  ! Now we tidy up any remaining `pointer' data, do this in order of variable
  ! declaration above.

  ! processingRange needs no deallocation
  DEALLOCATE(l1bInfo%L1BRADIDs)
  ! CALL DestroyMLSCFInfo(l2cf) MLSCF is currently fixed arrays. 
  CALL DestroyGridTemplateDatabase(aprioriData)
  DEALLOCATE(chunks)
  ! vectors, vectorTemplates and qtyTemplates destroyed at the 
  ! end of each chunk
  CALL DestroyL2GPDatabase(l2gpDatabase)
  CALL DestroyL2AUXDatabase(l2auxDatabase)
  CALL MLSMessage(MLSMSG_Info,ModuleName, &
       & "EOS MLS Level 2 data processing ended")
   STOP
!=============================================================================
END PROGRAM MLSL2
!=============================================================================

!
! $Log$
! Revision 2.0  2000/09/05 18:57:03  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.16  2000/06/29 23:55:35  lungu
! Uncommented:
!   CALL DestroyL2GPDatabase(l2gpDatabase)
!   CALL DestroyL2AUXDatabase(l2auxDatabase).
!

