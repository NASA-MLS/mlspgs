

! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.


!=============================================================================
MODULE ReadParseL2cf         ! Read and Parse the l2 Configuration File
!=============================================================================

USE SDPToolkit
USE Hdf
USE MLSMessageModule
USE MLSStrings


IMPLICIT NONE
PUBLIC


PRIVATE :: Id, ModuleName


!------------------- RCS Ident Info -----------------------
CHARACTER(LEN=130) :: Id = &                                                    
"$Id$"
CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!----------------------------------------------------------
 
INTEGER, PARAMETER :: MaxKeyLen = 25
INTEGER, PARAMETER :: MaxTypeLen = 6
INTEGER, PARAMETER :: MaxNOL2cfKeys = 20
INTEGER, PARAMETER :: MaxCharValueLen =25
INTEGER, PARAMETER :: UnitsLen = 10
INTEGER, PARAMETER :: L2cfEntryLen = 25
INTEGER, PARAMETER :: MaxNoKeysPerEntry = 10
INTEGER, PARAMETER :: MaxNoEntriesPerSection = 15

TYPE l2cfKey
  CHARACTER (len=MaxKeyLen) :: Keyword
  CHARACTER (len=MaxTypeLen):: TYPE
  INTEGER :: Keylen
END TYPE l2cfKey



TYPE L2cfCell

   ! This datatype defines a L2cf "cell", i.e. a contruction of type:
   ! "keyword=value" or "keyword=LowerBound..UpperBound[units]"

   ! First the "keyword" holder:

   CHARACTER (len=MaxKeyLen) :: Keyword
   ! Value holder for "real" type

   DOUBLE PRECISION :: RealValue

   ! Value holder for integer type

   INTEGER :: IntValue

   ! Value holder for "character" type

   CHARACTER (len= MaxCharValueLen) :: CharValue

   ! Value holder for "units" field

   CHARACTER (len=UnitsLen) :: Units

   ! Value holder for upper and lower bounds fields

   DOUBLE PRECISION :: RangeLowerBound
   DOUBLE PRECISION :: RangeUpperBound


END TYPE L2cfCell


TYPE L2cfEntry

   ! This datatype defines am L2cf entry, i.e. a name followed by
   ! a list of "keyword=value" cells
 
   ! Entry name, e.g. aprioriTemp

   CHARACTER (len=L2cfEntryLen) L2cfEntryName

   ! Number of "keyword=value" cells

   INTEGER L2cfEntryNoKeys

   !Cells holder 

   TYPE (L2cfCell), DIMENSION(MaxNoKeysPerEntry) :: Cells

END TYPE L2cfEntry


TYPE L2cfSection


   ! For each section:
   ! Name, Number of entries within section
   CHARACTER (len=L2cfEntryLen) L2cfSectionName
   INTEGER NoSectionEntries

   ! The actual entries:

   TYPE (L2cfEntry), DIMENSION(MaxNoEntriesPerSection) :: Entries ! (NoSectionEntries)


END TYPE L2cfSection




TYPE L2cf

! This datatype defines the l2cf as an array of L2cfSections

  INTEGER NoSections

  TYPE (L2cfSection), POINTER, DIMENSION (:) :: Sections ! (NoSections) 

END TYPE l2cf



! --------------------------------------------------------------------------

CONTAINS


! This function extracts a named section from the l2cf

!========================================
FUNCTION GetL2CFSection(l2cfInfo,sectionName)
!========================================

  ! Dummy arguments
  TYPE (l2cf), INTENT(IN) :: l2cfInfo
  CHARACTER (LEN=*), INTENT(IN) :: sectionName

  ! Function result
  TYPE (l2cfSection) :: GetL2CFSection

  ! Local variables
  INTEGER :: index

  ! Executable code

  index=LinearSearchStringArray(l2cfInfo%sections%l2cfSectionName,sectionName,&
       & caseInsensitive=.TRUE.)
  IF (index==0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
       & "No such l2cf section:"//sectionName)

  GetL2CFSection=l2cfInfo%sections(index)
END FUNCTION GetL2CFSection
  
!========================================
SUBROUTINE Strip_Blanks (line, nc, ncnb)
!========================================
 

CHARACTER (LEN=*) :: line
INTEGER :: nc, ncnb
INTEGER :: i, ncc

LOGICAL :: done = .FALSE.

ncc = nc

! check for blank line
line = TRIM(line)
IF (LEN(line) .EQ.  0)THEN
  ncnb = 0
  RETURN
END IF


! ignore trailing comments

i= INDEX(line, ';')
IF(  i /= 0)THEN
  ncc = i-1
END IF

! replace double blanks by single ones

DO WHILE (.NOT. done)
  i = INDEX(line(1:ncc), '  ')
  IF (i .EQ. 0)THEN
    done =.TRUE.
    ncnb = ncc
  ELSE 
    line(i:ncc-1)= line (i+1:ncc)
    ncc = ncc-1
  END IF
END DO

RETURN
END SUBROUTINE Strip_Blanks


!===================================
SUBROUTINE read_parse_l2cf (L2CFUnit, L2cf_data) 
!===================================

! This subroutine reads and parses L2CF and does some preliminary checking on 
! tokens.

! Arguments

INTEGER, INTENT (IN):: L2CFUnit 
TYPE (l2cf), INTENT (OUT) :: L2cf_data






! Parameters


CHARACTER (LEN=*), PARAMETER :: lineFmt = "(q,<nc>a1)"



! Functions


! Variables

INTEGER ::  returnStatus

TYPE (l2cfKey), DIMENSION(MaxNoL2cfKeys) :: L2cfTable

CHARACTER (LEN=256) :: line
CHARACTER (LEN=256) :: msg
CHARACTER (LEN=2048) :: buff
CHARACTER (LEN=L2cfEntryLen) :: mnemonic, section, AnEntry
CHARACTER (LEN=maxKeyLen) :: key
CHARACTER (LEN=MaxCharValueLen)Value
INTEGER :: nc, ncnb, ncb, version, isection, i, ib, ie, igs, ios, j, k, processL2CF
INTEGER :: CellIndex, ValueLen, keyLen
REAL :: tmp
LOGICAL :: eof = .FALSE.
LOGICAL :: found = .FALSE.
LOGICAL :: flag_section = .FALSE.

L2cfTable(1)%Keyword = 'InputVersionString       '
L2cfTable(1)%Type = 'string'
L2cfTable(1)%Keylen =  18
L2cfTable(2)%Keyword = 'OutputVersionString      '
L2cfTable(2)%Type = 'string'
L2cfTable(2)%Keylen =  19
L2cfTable(3)%Keyword = 'AllowClimatologyOverloads'
L2cfTable(3)%Type = 'string'
L2cfTable(3)%Keylen =  25
L2cfTable(4)%Keyword = 'source                   '
L2cfTable(4)%Type = 'string'
L2cfTable(4)%Keylen =  6
L2cfTable(5)%Keyword = 'length                   '
L2cfTable(5)%Type = 'real  '
L2cfTable(5)%Keylen =  6
L2cfTable(6)%Keyword = 'versionRange             '
L2cfTable(6)%Type = 'string'
L2cfTable(6)%Keylen =  12
L2cfTable(7)%Keyword = 'species                  '
L2cfTable(7)%Type = 'string'
L2cfTable(7)%Keylen =  7
L2cfTable(8)%Keyword = 'range                    '
L2cfTable(8)%Type = 'range '
L2cfTable(8)%Keylen =  5
L2cfTable(9)%Keyword = 'method                   '
L2cfTable(9)%Type = 'string'
L2cfTable(9)%Keylen =  6
L2cfTable(10)%Keyword = 'scale                    '
L2cfTable(10)%Type = 'real  '
L2cfTable(10)%Keylen =  5
L2cfTable(11)%Keyword ='height                   '
L2cfTable(11)%Type ='real  '
L2cfTable(11)%Keylen = 6
L2cfTable(12)%Keyword ='Overlap                  '
L2cfTable(12)%Type ='int   '
L2cfTable(12)%Keylen = 7
L2cfTable(13)%Keyword ='IdealLength              '
L2cfTable(13)%Type ='real  '
L2cfTable(13)%Keylen = 11
L2cfTable(14)%Keyword ='HomeLatitude             '
L2cfTable(14)%Type ='real  '
L2cfTable(14)%Keylen = 12
L2cfTable(15)%Keyword ='CriticalScanningModules  '
L2cfTable(15)%Type ='string'
L2cfTable(15)%Keylen = 23
L2cfTable(16)%Keyword ='CriticalBands            '
L2cfTable(16)%Type = 'string'
L2cfTable(16)%Keylen = 13
L2cfTable(17)%Keyword ='MaxGap                   '
L2cfTable(17)%Type ='real  '
L2cfTable(17)%Keylen = 6
L2cfTable(18)%Keyword = 'type                    '
L2cfTable(18)%Type ='string' 
L2cfTable(18)%Keylen = 4
L2cfTable(19)%Keyword ='fraction                 '
L2cfTable(19)%Type ='real  '
L2cfTable(19)%Keylen = 8
L2cfTable(20)%Keyword ='interpolationfactor      '
L2cfTable(20)%Type ='real  '
L2cfTable(20)%Keylen = 19

version = 1
isection = 0
igs=0
! Open the L2CF as a generic file for reading

returnStatus = Pgs_io_gen_openF (L2CFUnit, PGSd_IO_Gen_RSeqFrm, 0, &
                                processL2CF, version)

IF (returnStatus /= PGS_S_SUCCESS) THEN

  CALL Pgs_smf_getMsg(returnStatus, mnemonic, msg)
  CALL MLSMessage (MLSMSG_Error, &
                  ModuleName, "Error opening L2CF:  "//mnemonic//" "//msg)

ENDIF

! Read  input
DO WHILE (.NOT. eof)

  READ(UNIT=processL2CF, IOSTAT=ios, FMT=lineFmt) nc, line
  IF (ios .eq. 0) THEN

! strip leading, trailing and multiple blanks from line
 
    CALL Strip_Blanks(line, nc, ncnb)

    IF (INDEX(line, 'BEGIN' ) /= 0 ) THEN
! Is this a BEGIN Section Line?
       isection = isection + 1
       section(1:) = line (7:ncnb)
       flag_section =.false.
    ELSE IF (isection .eq. 0 .and. ncnb .gt. 0)Then
! If not and there was no BEGIN before and this is not a comment line, log error and exit
        CALL MLSMessage (MLSMSG_Error, ModuleName, 'Missing BEGIN')
     
    END IF

    ib = 1
    READ(UNIT=processL2CF, IOSTAT=ios, FMT=lineFmt) nc, line
    
    IF (ios .eq.  0) THEN

      CALL Strip_Blanks(line, nc, ncnb)

! ignore comment lines

      IF(INDEX(line, ';' ) /= 1) THEN

        IF(INDEX(line, 'END') /= 1) THEN

          buff(ib:ib+ncnb-1) = line(1:ncnb)
          ib = ib + ncnb

! copy continuation lines into buff

          DO WHILE (line(ncnb:ncnb) .eq. '$')
            ib = ib- 1

            READ(UNIT=processL2CF, IOSTAT=ios, FMT=lineFmt) nc, line
            IF (ios /= 0) THEN

              CALL Strip_Blanks(line, nc, ncnb)
              buff(ib:ib+ncnb-1) = line(1:ncnb)
              ib = ib+ncnb
            ELSE
              eof = .TRUE.
            END IF

          END DO
          ncb = ib

! Proces l2cf entry

          AnEntry(1:L2cfEntryLen) =' '
          ib = INDEX(buff, ',')

          IF (ib .eq.  0) THEN
            CALL MLSMessage (MLSMSG_Error, ModuleName, "Comma expected, : "//buff(1:ncb))
          END IF          

          AnEntry(1:ib-1)=buff(1:ib-1)
          ib = ib + 1
          igs = igs + 1
          L2cf_data%Sections(isection)%Entries(igs)%L2cfEntryName = AnEntry
          CellIndex = 1

! Process keyword = value pairs

          DO WHILE (ib < ncb)
            ie = INDEX(buff(ib:ncb), '=')

            IF (ie .eq. 0) THEN
              CALL MLSMessage (MLSMSG_Error, ModuleName, "= expected, : "//buff(1:ncb))
            END IF
 
            key(1:maxKeyLen) = ' '
            Key = buff(ib: ie-1) 
            keyLen = ie-ib
            ib = ie + 1
            ie = INDEX(buff(ib:ncb), ',')

            IF (ie .eq. 0) THEN
              CALL MLSMessage (MLSMSG_Error, ModuleName, "Comma expected, : "//buff(ib:ncb))
            END IF 

            Value (1:MaxCharValueLen) = ' '
            Value = buff(ib: ie-1)
            ValueLen = ie-ib
            ib = ie + 1

            found = .FALSE.
            i = 1

            DO WHILE ((.NOT. found) .AND. (i <= MaxNoL2cfKeys))

              IF(key(1:KeyLen) .EQ. l2cfTable(i)%Keyword(1:l2cfTable(i)%Keylen))THEN
                found = .TRUE.
              ELSE
                i = i + 1
              END IF
            END DO

            IF(found) THEN

                
               L2cf_data%Sections(isection)%Entries(igs)%L2cfEntryNoKeys = &
               L2cf_data%Sections(isection)%Entries(igs)%L2cfEntryNoKeys + 1
               L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%Keyword= Key(1:Keylen)

               IF(l2cfTable(i)%Type .eq. 'real')THEN
                 j = 1     
                 DO WHILE((j <= ValueLen) .and. (Value(j:j) /= ' '))
                   j = j + 1
                 END DO
                 READ (unit=Value(1:j-1), fmt=*)L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RealValue
                 L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%units=Value(j:ValueLen)
               ELSE IF(l2cfTable(i)%Type .eq. 'int')THEN
                 j = 1     
                 DO WHILE((j <= ValueLen) .AND. (Value(j:j) /= ' '))
                   j = j + 1
                 END DO
                 READ (unit=Value(1:j-1), fmt=*)L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%IntValue
                 L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%units=Value(j:ValueLen)
               ELSE IF(l2cfTable(i)%Type .eq. 'range')THEN
                 j = INDEX(Value, '..')
                 IF (j > 0)THEN
                   READ (unit=Value(1:j-1), fmt=*) &
                   L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeLowerBound
                   k=j+2

                   DO WHILE((k <= ValueLen) .AND. (Value(k:k) /= ' '))
                    k = k + 1
                   END DO  
            
                   READ (unit=Value(1:j-1), fmt=*) &
                   L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeUpperBound
                   L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%units=Value(k:ValueLen)

                   IF(L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeUpperBound < &
                      L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeUpperBound)THEN

                     tmp = L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeUpperBound
                     L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeUpperBound = &
                     L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeLowerBound
                     L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeLowerBound = tmp

                   END IF ! (L2cf_data%Sections(isection)%Entries
                 ELSE
                   CALL MLSMessage (MLSMSG_Error, ModuleName, "Range expected, : "//Value)
                 END IF ! (j > 0)

               ELSE IF(l2cfTable(i)%Type .EQ. 'string')THEN

                 j = INDEX(Value, '"')
                 IF(j > 0)THEN
                   k = INDEX(Value(j+1:ValueLen), '"')
                   IF(k > 0)THEN
                      L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%CharValue=Value(j:k)
                   ELSE
                     CALL MLSMessage (MLSMSG_Error, ModuleName, "Quote expected, : "//Value(1:ValueLen))
                   END IF
                 ELSE
                   CALL MLSMessage (MLSMSG_Error, ModuleName, "Quote expected, : "//Value(1:ValueLen))
                 END IF ! j > 0
    
               END IF ! l2cfTable(i)%Type

            ELSE            
              CALL MLSMessage (MLSMSG_Error, ModuleName, "Illegal Keyword : " //key(1:KeyLen))
            END IF ! (found)                      
            CellIndex = CellIndex + 1  
          END DO ! ib < ncb
          igs = igs + 1
        ELSE
          if (  L2cf_data%Sections(isection)%Entries(igs)%L2cfEntryName .ne. line(5:ncnb))THEN
             CALL MLSMessage (MLSMSG_Warning, ModuleName, "BEGIN and END Section Names don't match " // &
                              L2cf_data%Sections(isection)%Entries(igs)%L2cfEntryName//' , '//line(5:ncnb))
          end if
          isection = isection + 1
          flag_section = .true.
! Check whether END Section matches BEGIN
        END IF !INDEX(line, 'END') /= 1  
      ELSE
        eof = .TRUE.
      END IF ! (INDEX(line, ';' ) /= 1)
    ELSE

      eof = .true.
      if(.not. flag_section) then
        CALL MLSMessage (MLSMSG_Warning, ModuleName, 'Premature EOF encountered')
      end if
    END IF ! (ios .eq.  0)

  ELSE

    eof = .true.
    if(.not. flag_section) then
       CALL MLSMessage (MLSMSG_Warning, ModuleName, 'Premature EOF encountered')
    end if
  END IF ! (ios .eq.  0)
 

END DO ! (.not. eof)

! Close L2CF

returnStatus = Pgs_io_gen_closeF (processL2CF)

IF (returnStatus /= PGS_S_SUCCESS) THEN
  CALL Pgs_smf_getMsg(returnStatus, mnemonic, msg)
  CALL MLSMessage (MLSMSG_Error, ModuleName, 'Error closing L2CF:  '//mnemonic//' '//msg)
 
ENDIF




!==============================
END SUBROUTINE  read_parse_l2cf
!==============================

!=======================
END MODULE ReadParseL2cf
!=======================

! $Log$
! Revision 1.7  2000/01/07 19:30:04  lungu
! Revoved returnStatus from argument list.
!
! Revision 1.6  2000/01/06 23:39:00  lungu
! Finalized error handling.
!
! Revision 1.5  2000/01/06 02:37:22  lungu
! Made a module out of read_parse_l2cf, strip_blanks and data structures.
!
! Revision 1.4  2000/01/04 19:46:23  lungu
! Changed error handling using MLSMessage
!
