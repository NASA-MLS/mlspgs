! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.


!===================================
SUBROUTINE read_parse_l2cf (L2cf_data, returnStatus) 
!===================================

USE l2cf
USE MLSMessageModule

IMPLICIT NONE


!------------------- RCS Ident Info -----------------------
CHARACTER(LEN=130) :: Id = &                                                    
"$Id$"
CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This subroutine reads and parses L2CF and does some preliminary checking on 
! tokens.

! Arguments

TYPE (l2cf), INTENT (OUT) :: L2cf_data
! Parameters


CHARACTER (LEN=*), PARAMETER :: lineFmt = "(q,<nc>a1)"

INTEGER, PARAMETER :: L2CF = 3000

! Functions

INTEGER, EXTERNAL :: Pgs_io_gen_openF, Pgs_io_gen_closeF

! Variables


CHARACTER (LEN=256) :: line
CHARACTER (LEN=256) :: msg
CHARACTER (LEN=2048) :: buff



INTEGER :: nc, returnStatus, version, isection, igs

LOGICAL :: eof = .FALSE.

version = 1
isection = 0
igs=0
! Open the L2CF as a generic file for reading

returnStatus = Pgs_io_gen_openF (L2CF, PGSd_IO_Gen_RSeqFrm, 0, &
                                processL2CF, version)

IF (returnStatus /= PGS_S_SUCCESS) THEN

  CALL Pgs_smf_getMsg(returnStatus, mnemonic, msg)
  CALL MSMessage (MLSMSG_Error,
         & ModuleName, "Error opening L2CF:  "//mnemonic//" "//msg)

ENDIF

! Read  input
DO WHILE (.NOT. eof)

  READ(UNIT=processL2CF, IOSTAT=ios, FMT=lineFmt) nc, line
  IF (ios /= 0) THEN

! strip leading, trailing and multiple blanks from line
 
    CALL Strip_Blanks(line, nc, ncnb)

    IF (INDEX(line, 'BEGIN') /= 0) THEN
       isection = isection + 1
       section(1:) = line (7:ncnb)
    END IF

    ib = 1
    READ(UNIT=processL2CF, IOSTAT=ios, FMT=lineFmt) nc, line
    IF (ios /= 0) THEN

      CALL Strip_Blanks(line, nc, ncnb)

! ignore comment lines

      IF(INDEX(line, ';' ) /= 1) THEN

        IF(INDEX(line, 'END') /= 1) THEN

          buff(ib:ib+ncnb-1) = line(1:ncnb)
          ib = ib + ncnb

! copy continuation lines into buff

          DO WHILE (line(ncnb:ncnb) = '$')
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

          IF (ib = 0) THEN
            CALL MSMessage (MLSMSG_Error, ModuleName, "Comma expected, : "//buff(1:ncb))
          END IF          

          AnEntry(1:ib-1)=buff(1:ib-1)
          ib = ib + 1
          igs = igs + 1
          L2cf_data%Sections(isection)%Section(igs)%L2cfEntryName = AnEntry
          CellIndex = 1

! Process keyword = value pairs

          DO WHILE (ib < ncb)
            ie = INDEX(buff(ib:ncb), '=')

            IF (ie = 0) THEN
              CALL MSMessage (MLSMSG_Error, ModuleName, "= expected, : "//buff(1:ncb))
            END IF
 
            key(1:maxKeyLen) = ' '
            Key = buff(ib, ie-1) 
            keyLen = ie-ib
            ib = ie + 1
            ie = INDEX(buff(ib:ncb), ',')

            IF (ie = 0) THEN
              CALL MSMessage (MLSMSG_Error, ModuleName, "Comma expected, : "//buff(ib:ncb))
            END IF 

            Value (1:MaxCharValueLen) = ' '
            Value = buff(ib, ie-1)
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
               L2cf_data%Section(isection)%Entries(igs)%L2cfEntryNoKeys + 1
               L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%Keyword= Key(1:Keylen)

               IF(l2cfTable(i)%Type = "real")THEN
                 j = 1     
                 DO WHILE((j <= ValueLen) = (Value(j:j) /= ' '))
                   j = j + 1
                 END DO
                 READ (unit=Value(1:j-1), fmt='*') L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RealValue
                 L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%units=Value(j:ValueLen)
               ELSE IF(l2cfTable(i)%Type = "int")THEN
                 j = 1     
                 DO WHILE((j <= ValueLen) .AND. (Value(j:j) /= ' '))
                   j = j + 1
                 END DO
                 READ (unit=Value(1:j-1), fmt='*') L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%IntValue
                 L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%units=Value(j:ValueLen)
               ELSE IF(l2cfTable(i)%Type = "range")THEN
                 j = INDEX(Value, '..')
                 IF (j > 0)THEN
                   READ (unit=Value(1:j-1), fmt='*') &
                   L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeLowerBound
                   k=j+2
                   DO WHILE((k <= ValueLen) .AND. (Value(k:k) /= ' '))
                    k = k + 1
                   END DO              
                   READ (unit=Value(1:j-1), fmt='*') &
                   L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeUpperBound
                   L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%units=Value(k:ValueLen)

                   IF(L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeUpperBound < &
                      L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeUpperBound)THEN

                     tmp = L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeUpperBound
                     L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeUpperBound = 
                     L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeLowerBound
                     L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeLowerBound = tmp

                   END IF
                 ELSE
                   CALL MSMessage (MLSMSG_Error, ModuleName, "Range expected, : "//Value
                 END IF

               ELSE IF(l2cfTable(i)%Type .EQ. 'string')THEN

                 j = INDEX(Value, '"')
                 IF(j > 0)THEN
                   k = INDEX(Value(j+1:ValueLen)
                   IF(k > 0)THEN
                      L2cf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%CharValue=Value(j:k)
                   ELSE
                     CALL MSMessage (MLSMSG_Error, ModuleName, "Quote expected, : "//Value(1:ValueLen))
                   END IF
                 ELSE
                   CALL MSMessage (MLSMSG_Error, ModuleName, "Quote expected, : "//Value(1:ValueLen))
                 END IF
    
               END IF
               
                                  
               
         END DO
       ELSE
! Check whether END Section matches BEGIN
       END IF  
     ELSE
       eof = .TRUE.
     END IF

   END IF

 END IF

END DO

! Close L2CF

returnStatus = Pgs_io_gen_closeF (processL2CF)

IF (returnStatus /= PGS_S_SUCCESS) THEN
  CALL Pgs_smf_getMsg(returnStatus, mnemonic, msg)
  CALL MSMessage (MLSMSG_Error, ModuleName, 'Error closing L2CF:  ', mnemonic//' '//msg)
 
ENDIF




!===================
END subroutine  read_parse_l2cf
!===================

! $Log$
